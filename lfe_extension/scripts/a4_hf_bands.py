#!/usr/bin/env python
"""Track A-2 follow-up (1b) — EMPIRICAL band scan + noise-floor diagnostic for the PGC HF-ratio veto.

Recomputes the veto statistic R = E(hi)/E(2-8 Hz) at CN.PGC, from the NATIVE-rate traces, for a set of
candidate `hi` bands, on detections of a BH epoch and an HH epoch, using exactly the a1_hf_veto.py
placement (t + ttP + sp - 0.5 + moveout, window = L s, 3 components summed, cumulative-sum energies).

It also measures, per day, the BAND ENERGY OF THE WHOLE TRACE divided by the band-averaged squared
instrument response, i.e. the ground-motion-referred power in 2-8 and in the hi band.  If the epoch
shift of R were an anti-alias/band effect it would disappear once the response is divided out; if it is
a sensitivity/noise-floor effect (the 2017 sensor swap cut PGC's gain by 10.6x, so the digitiser noise
floor sits 10.6x higher in ground-motion units) it will NOT.

Output: results/a4_hf_band_scan.csv, results/a4_hf_noisefloor.csv
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import d3_vel1d as v1d                                        # noqa: E402
from a1_hf_veto import load_native, bp, epi_km, PRE           # noqa: E402

EXT = os.path.dirname(HERE)
REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
IN = f"{REP}/results"
BANDS = [(10., 19.), (10., 16.), (10., 14.), (8., 16.), (12., 16.)]


def native_any(netsta, year, doy):
    """a1_hf_veto.load_native, but searching BOTH trace roots (lfe_extension first)."""
    from obspy import read
    import e4_chanorder as chanorder
    for root in (f"{EXT}/traces", f"{REP}/traces"):
        p = f"{root}/{netsta}/{year}/{doy:03d}.mseed"
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            continue
        st = read(p)
        for tr in st:
            r = float(tr.stats.sampling_rate)
            if r > 0 and abs(r - round(r)) / r < 1e-4 and r != round(r):
                tr.stats.sampling_rate = float(round(r))
        st.merge(method=1, fill_value=0)
        order = chanorder.resolve_order({tr.stats.channel for tr in st}, station=netsta, strict=True)
        keep = [c for c in order if c is not None]
        npts = min(min(tr.stats.npts for tr in st.select(channel=c)) for c in keep)
        X = np.vstack([st.select(channel=c)[0].data[:npts].astype(np.float64) for c in keep])
        return X, float(st[0].stats.sampling_rate), float(st.select(channel=keep[0])[0].stats.starttime.timestamp)
    return None


def placements(fams, sta, L):
    F = pd.read_csv(f"{IN}/e5r_families_2012.csv").set_index("bkey")
    S = pd.read_csv(f"{IN}/e5r_stations.csv").set_index("sta")
    z = np.load(f"{IN}/e5r_templates_2012.npz", allow_pickle=True); sp = z["sp"].item()
    MV = pd.read_csv(f"{IN}/e4_moveout_e5r_mv_build.csv")
    mv = {(float(r.fam), r.sta, float(r.L)): float(r.lag_median) for r in MV.itertuples() if np.isfinite(r.lag_median)}
    off = {}
    for fk in fams:
        fk = float(fk)
        if (fk, sta) not in sp:
            continue
        m = F.loc[int(fk)]
        dx = epi_km(float(m.lat), float(m.lon), float(S.loc[sta, "lat"]), float(S.loc[sta, "lon"]))
        off[fk] = float(v1d.tt(dx, float(m.dep), "P")) + float(sp[(fk, sta)]) - PRE + mv.get((fk, sta, L), 0.0)
    return off


def kband(resp_f, resp_a, f1, f2):
    m = (resp_f >= f1) & (resp_f <= f2)
    return float(np.mean(resp_a[m] ** 2)) if m.any() else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--legs", default="2015-12:CN.PGC:2015,2018-06:CN.PGC:2018",
                    help="comma list <trackC month>:<netsta>:<year>")
    ap.add_argument("--sta", default="PGC")
    ap.add_argument("--L", type=float, default=12.0)
    ap.add_argument("--max-days", type=int, default=8)
    a = ap.parse_args()

    from obspy.clients.fdsn import Client
    inv = Client("EARTHSCOPE", timeout=180).get_stations(network="CN", station=a.sta, level="response",
                                                         starttime="2010-01-01", endtime="2026-09-16")
    resp = {}
    for n in inv:
        for s in n:
            for ch in s:
                if ch.code[1] == "H" and ch.code[-1] == "Z":
                    f = np.linspace(0.5, min(0.49 * ch.sample_rate, 49.0), 1500)
                    resp[(ch.code[:2], str(ch.start_date)[:10])] = (
                        f, np.abs(ch.response.get_evalresp_response_for_frequencies(f, output="VEL")))

    rows, nf_rows = [], []
    for leg in a.legs.split(","):
        month, netsta, year = leg.split(":")
        year = int(year)
        det = pd.read_parquet(f"{EXT}/results/trackC/detections_{month}.parquet")
        off = placements(det.fam.unique(), a.sta, a.L)
        # representative days: the busiest ones, since the catalog median is dominated by them
        doys = list(det.groupby("doy").size().sort_values(ascending=False).index[:a.max_days])
        doys = sorted(int(d) for d in doys)
        R = {b: [] for b in BANDS}
        for doy in doys:
            got = native_any(netsta, year, int(doy))
            if got is None:
                print(f"  [{month} doy {doy}] no native trace", flush=True)
                continue
            X, fs, t0 = got
            g = det[det.doy == doy]
            key = [k for k in resp if k[0] == ("BH" if fs == 40 else "HH")]
            key = sorted(key, key=lambda k: k[1])
            key = [k for k in key if k[1] <= f"{year}-12-31"][-1]
            rf, ra = resp[key]
            elo_full = (bp(X, fs, 2.0, 8.0) ** 2).sum(0)
            n = int(round(a.L * fs))
            clo = np.concatenate([[0.0], np.cumsum(elo_full)])
            ok = np.array([f in off for f in g.fam.to_numpy()])
            i0 = np.full(len(g), -1, np.int64)
            i0[ok] = np.rint((g.t.to_numpy()[ok] + np.array([off[f] for f in g.fam.to_numpy()[ok]]) - t0) * fs).astype(np.int64)
            good = ok & (i0 >= 0) & (i0 + n <= len(elo_full))
            klo = kband(rf, ra, 2.0, 8.0)
            nf = dict(month=month, doy=int(doy), fs=fs, chan=key[0], resp_start=key[1],
                      med_e_lo_counts=float(np.median(elo_full)), k_lo=klo,
                      med_e_lo_ground=float(np.median(elo_full)) / klo)
            for (f1, f2) in BANDS:
                if f2 > 0.49 * fs:
                    continue
                ehi_full = (bp(X, fs, f1, f2) ** 2).sum(0)
                chi = np.concatenate([[0.0], np.cumsum(ehi_full)])
                num = chi[i0[good] + n] - chi[i0[good]]
                den = clo[i0[good] + n] - clo[i0[good]]
                r = np.where(den > 0, num / np.maximum(den, 1e-30), np.nan)
                R[(f1, f2)].append(pd.DataFrame(dict(
                    month=month, doy=int(doy), fs=fs, band=f"{f1:g}-{f2:g}",
                    fam=g.fam.to_numpy()[good], t=g.t.to_numpy()[good],
                    net_mad=g.net_mad.to_numpy()[good], n_sites=g.n_sites.to_numpy()[good],
                    tremor=g.tremor_flag.to_numpy()[good] if "tremor_flag" in g else True, hfr=r)))
                khi = kband(rf, ra, f1, f2)
                nf[f"med_e_hi_counts_{f1:g}-{f2:g}"] = float(np.median(ehi_full))
                nf[f"k_hi_{f1:g}-{f2:g}"] = khi
                nf[f"med_e_hi_ground_{f1:g}-{f2:g}"] = float(np.median(ehi_full)) / khi
            nf_rows.append(nf)
            print(f"  [{month} doy {doy}] fs {fs:g}, {int(good.sum())}/{len(g)} windows", flush=True)
        for b, lst in R.items():
            if lst:
                rows.append(pd.concat(lst, ignore_index=True))
    D = pd.concat(rows, ignore_index=True)
    D.to_csv(f"{EXT}/results/a4_hf_band_scan.csv", index=False)
    N = pd.DataFrame(nf_rows)
    N.to_csv(f"{EXT}/results/a4_hf_noisefloor.csv", index=False)

    pd.set_option("display.width", 220)
    print("\n[band scan] median HF ratio of the detections, by band and epoch")
    piv = D.groupby(["band", "month"]).hfr.median().unstack()
    piv["shift"] = piv.iloc[:, 1] / piv.iloc[:, 0]
    print(piv.round(3).to_string())
    print("\n[noise floor] whole-trace band power referred to GROUND MOTION, log10(counts^2 / <|Resp|^2>)")
    cols = [c for c in N.columns if c.startswith("med_e_lo_ground") or c.startswith("med_e_hi_ground")]
    G = N.groupby(["month", "chan", "fs"])[cols].median()
    print(np.log10(G).round(3).to_string())
    print("\n[noise floor] HF/LF ground-motion power ratio per epoch (whole trace, response removed)")
    print(pd.DataFrame({c.replace("med_e_hi_ground_", ""): G[c] / G["med_e_lo_ground"]
                        for c in cols if c != "med_e_lo_ground"}).round(4).to_string())
    print("\n[write] results/a4_hf_band_scan.csv, results/a4_hf_noisefloor.csv")


if __name__ == "__main__":
    main()
