#!/usr/bin/env python
"""g_trigger_eq (2026-10-11) -- Study G step 1: earthquake list, arrival windows, estimated peak dynamic stress.
Rules: svi_science/results/g_trigger/PLAN.md §0-§2.  CPU only, small.

    python svi_science/scripts/g_trigger_eq.py   -> results/g_trigger/tables/eq_list_raw.csv
"""
import numpy as np
import pandas as pd
from obspy.taup import TauPyModel

ROOT = "/home/jovyan/tremorferometry"
G = f"{ROOT}/svi_science/results/g_trigger"


def refs():
    F = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv")
    W = pd.read_csv(f"{ROOT}/svi_science/results/a_wa_between/families_wa.csv")
    return {"SVI": (F.lat.mean(), F.lon.mean()), "WA": (W.lat.mean(), W.lon.mean())}


def gc_deg(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    c = np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(lo2 - lo1)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def est_sigma_kpa(M, ddeg):
    d = np.maximum(ddeg, 2.0)
    AT = 10 ** (M - 1.66 * np.log10(d) - 3.3)          # µm/s  (A/T)
    pgv = 2 * np.pi * AT * 1e-6                           # m/s
    eps = pgv / 3500.0
    return pgv, eps, 30e9 * eps / 1e3                     # kPa


def main():
    A = pd.read_csv(f"{G}/dl/global_m7.csv"); A["src_query"] = "global_M7"
    B = pd.read_csv(f"{G}/dl/regional_m55.csv"); B["src_query"] = "regional_M5.5"
    E = pd.concat([A, B]).drop_duplicates("id").copy()
    E = E[E.type == "earthquake"]
    E["t0"] = pd.to_datetime(E.time, utc=True).dt.tz_convert(None).astype("datetime64[ns]")
    E = E.sort_values("t0").reset_index(drop=True)
    R = refs()
    m = TauPyModel("iasp91")
    out = []
    for _, e in E.iterrows():
        row = dict(id=e.id, t0_utc=e.t0, lat=e.latitude, lon=e.longitude, depth_km=e.depth, mag=e.mag,
                   magType=e.magType, place=e.place, src_query=e.src_query)
        for reg, (la, lo) in R.items():
            dd = gc_deg(e.latitude, e.longitude, la, lo)
            dkm = dd * 111.195
            arr = m.get_travel_times(source_depth_in_km=max(0.0, min(float(e.depth), 700.0)), distance_in_degree=dd,
                                     phase_list=["ttbasic"])
            tp = min(a.time for a in arr) if arr else dkm / 6.0
            pgv, eps, sig = est_sigma_kpa(e.mag, dd)
            row.update({f"{reg}_ddeg": dd, f"{reg}_dkm": dkm, f"{reg}_tP_s": tp, f"{reg}_tL_s": dkm / 4.5,
                        f"{reg}_tR_s": dkm / 2.5, f"{reg}_pgv_est": pgv, f"{reg}_sig_est_kpa": sig,
                        f"{reg}_near": dd < 2.0})
        out.append(row)
    O = pd.DataFrame(out)
    O["t0_PT"] = O.t0_utc.dt.tz_localize("UTC").dt.tz_convert("America/Los_Angeles").dt.strftime("%Y-%m-%d %H:%M %Z")
    O.to_csv(f"{G}/tables/eq_list_raw.csv", index=False)
    print(len(O), "earthquakes;", "refs", R)
    print(O[["t0_PT", "mag", "place", "SVI_ddeg", "SVI_sig_est_kpa"]].sort_values("SVI_sig_est_kpa").tail(25).to_string())


if __name__ == "__main__":
    main()
