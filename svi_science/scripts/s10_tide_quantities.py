#!/usr/bin/env python
"""s10_tide_quantities -- S10 tide/sub-daily follow-up, step A1 (rules: results/s10_dvv/tide_subday/PLAN.md).
Per pilot row (data and c1): M2 phase of UDSS and of sigma_n (s5_lib, harmonic fit per era table), UDSS / sigma_n
values and hourly rates, fake-cycle phases. Also writes the merged hourly UDSS / sigma_n series of the three families
(for the random-time-shift control) and the 2014-01 continuity check of the two S5 tables.
"""
import os, sys, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s5_lib as L

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/s10_dvv/tide_subday"
FAMS = (10, 87, 17)
COMPS = ["tau_updip_Pa", "sigma_n_Pa"]
FAKE_H = (12.0, 13.5, 10.5)
T_FAKE0 = np.datetime64("2003-01-01T00:00:00", "ns")


def main():
    os.makedirs(f"{OUT}/tables", exist_ok=True)
    X = pd.read_parquet(f"{ROOT}/svi_science/results/s10_dvv/pilot/tables/rows.parquet")
    X["tt"] = pd.to_datetime(X.t, unit="s").astype("datetime64[ns]")
    Bo = pd.read_parquet(f"{ROOT}/svi_science/results/s5/tidal_stress_bostock_hourly.parquet", columns=["fid", "time"] + COMPS)
    Bo = Bo[Bo.fid.isin(FAMS)]
    Se = pd.read_parquet(f"{ROOT}/svi_science/results/s5_setup/tidal_stress_svi_hourly.parquet", columns=["fid", "time"] + COMPS)
    Se = Se[Se.fid.isin(FAMS)]
    # continuity check 2014-01
    m = Bo[(Bo.time >= "2014-01-01") & (Bo.time < "2014-02-01")].merge(Se, on=["fid", "time"], suffixes=("_b", "_s"))
    chk = {c: float(np.abs(m[c + "_b"] - m[c + "_s"]).max()) for c in COMPS}
    chk["n"] = int(len(m))
    json.dump(chk, open(f"{OUT}/tables/continuity_2014_01.json", "w"), indent=1)
    print("continuity", chk)
    assert chk["tau_updip_Pa"] < 1.0 and chk["sigma_n_Pa"] < 1.0, "S5 tables disagree in 2014-01"
    Sb = L.Stress(Bo[Bo.time < "2014-01-01"], COMPS)
    Ss = L.Stress(Se, COMPS)
    parts = []
    for era, St, msk in (("b", Sb, X.tt < np.datetime64("2014-01-01")), ("s", Ss, X.tt >= np.datetime64("2014-01-01"))):
        ev = X[msk].rename(columns={"fam": "fid", "tt": "t_"}).copy()
        ev["t"] = ev["t_"]
        es = L.EventSet(ev[["fid", "t", "row_id"]], St, COMPS, phase_comp="tau_updip_Pa")
        q = es.quantities(np.zeros(len(es.ev), np.int64))
        o = es.ev[["row_id"]].copy()
        o["th_udss"] = np.degrees(np.arctan2(q["sin_m2"], q["cos_m2"]))
        o["th_sn"] = np.degrees(np.arctan2(q["sin_sn_m2"], q["cos_sn_m2"]))
        o["udss"] = q["tau_updip_Pa"]; o["sn"] = q["sigma_n_Pa"]
        # M2 UDSS amplitude per family (for per-kPa conversion)
        hc = St.hc["tau_updip_Pa"]["m2"]
        o["m2_udss_amp_Pa"] = np.abs(hc)[es.fi]
        o["era_tab"] = era
        parts.append(o)
    Q = pd.concat(parts, ignore_index=True)
    # merged hourly series for value lookups (Bostock table < 2014, setup table >= 2014)
    H = pd.concat([Bo[Bo.time < "2014-01-01"], Se], ignore_index=True).sort_values(["fid", "time"])
    t0 = H.time.min(); nh = int((H.time.max() - t0) / pd.Timedelta("1h")) + 1
    arr = {}
    for c in COMPS:
        A = np.full((len(FAMS), nh), np.nan, np.float32)
        for k, f in enumerate(FAMS):
            g = H[H.fid == f]
            A[k, ((g.time - t0) / pd.Timedelta("1h")).astype(int).to_numpy()] = g[c].to_numpy()
        arr[c] = A
    np.savez(f"{OUT}/tables/hourly_3fam.npz", fams=np.array(FAMS), t0_ns=np.int64(pd.Timestamp(t0).value),
             udss=arr["tau_updip_Pa"], sn=arr["sigma_n_Pa"])
    X = X.merge(Q, on="row_id", how="left")
    # rates (Pa/h) from the merged hourly series, centred difference, linear interpolation
    tns = X.tt.to_numpy().astype(np.int64)
    fi = X.fam.map({f: k for k, f in enumerate(FAMS)}).to_numpy()
    h = (tns - pd.Timestamp(t0).value) / 3.6e12
    for nm, A in (("udss", arr["tau_updip_Pa"]), ("sn", arr["sigma_n_Pa"])):
        D = np.full_like(A, np.nan); D[:, 1:-1] = (A[:, 2:] - A[:, :-2]) / 2.0
        h0 = np.floor(h).astype(int); w = h - h0
        X[f"{nm}_rate"] = D[fi, h0] * (1 - w) + D[fi, h0 + 1] * w
        X[f"{nm}_chk"] = A[fi, h0] * (1 - w) + A[fi, h0 + 1] * w
    X["chk_diff"] = np.abs(X.udss - X.udss_chk)
    print("value lookup vs s5_lib max |diff| Pa:", float(X.chk_diff.max()))
    for P in FAKE_H:
        X[f"th_fake{P:g}"] = 360.0 * (((tns - T_FAKE0.astype(np.int64)) / (P * 3.6e12)) % 1.0)
    X["t_pt"] = pd.to_datetime(X.t, unit="s", utc=True).dt.tz_convert("America/Los_Angeles").dt.strftime("%Y-%m-%d %H:%M")
    X.drop(columns=["udss_chk", "sn_chk"]).to_parquet(f"{OUT}/tables/rows_tide.parquet", index=False)
    print(X.groupby(["fam", "series"]).agg(n=("row_id", "size"), udss_sd=("udss", "std"), m2amp=("m2_udss_amp_Pa", "mean"),
                                         nan=("th_udss", lambda s: int(s.isna().sum()))))
    open(f"{OUT}/QUANT_DONE", "w").write("ok\n")


if __name__ == "__main__":
    main()
