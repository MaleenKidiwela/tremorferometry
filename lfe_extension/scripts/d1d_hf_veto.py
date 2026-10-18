#!/usr/bin/env python
"""Track D copy of Track A step 1(b) (`a1_hf_veto.py`) with ONE change: the native traces are looked up in
BOTH `lfe_catalog_repair/traces` and `lfe_extension/traces`, so the frozen veto can be applied to the
production months of the D1(d) withheld-station runs.  Nothing else is modified.

Original docstring -- the HF-RATIO VETO test (Merlin plan-review F14, SPECULATIVE): are the quiet-day detections
site transients / uncatalogued micro-earthquakes, which carry far more 10-20 Hz energy than an LFE?

For every mergeB detection (existing run, no GPU pass) compute at the PGC SITE the energy ratio
    R = E(HF band) / E(2-8 Hz)
over the same 8-s window the matched filter scored, i.e. starting at
    t_detection + ttP(family, station) + sp(family, station) - PRE + moveout(family, station, L)
(identical placement formula to e4_template_test.py, PRE = 0.5 s), summed over the three components of the
NATIVE-RATE day trace (the 20 Hz cache cannot carry 10-20 Hz).  CN.PGC in 2012 is BH 40 Hz -> Nyquist 20 Hz, so the
PGC HF band is 10-19 Hz (0.95 Nyquist; a band that touches Nyquist makes the Butterworth unstable); PB.B011 is EH
100 Hz and is computed alongside with the literal 10-20 Hz band as a cross-check of the band choice.
Cut = the 99th percentile of R over BUILD-half detections matched to an own-family Bostock event (the veto is
calibrated on events we believe, on the half the templates were built from; the test half is never used to set it).
Reported: quiet-day (doy 228-242) FA reduction and TEST-half own-recall cost at each threshold.
Adopt only if quiet FA falls >= 30 % at <= 2 % recall cost (Merlin step 1)."""
from __future__ import annotations
import argparse, os, sys
import numpy as np, pandas as pd
from obspy import read
from scipy.signal import butter, sosfiltfilt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d3_vel1d as v1d
import e4_chanorder as chanorder

REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
IN = f"{REP}/results"
OUT = "/home/jovyan/tremorferometry/lfe_extension/results"
PRE = 0.5
QUIET = set(range(228, 243))
R_EARTH_KM_DEG = 111.195


def epi_km(lat1, lon1, lat2, lon2):
    cs = np.cos(np.radians(48.6))
    return np.hypot((np.asarray(lon1) - lon2) * R_EARTH_KM_DEG * cs, (np.asarray(lat1) - lat2) * R_EARTH_KM_DEG)


def bp(x, fs, f1, f2):
    return sosfiltfilt(butter(4, [f1 / (fs / 2), f2 / (fs / 2)], btype="band", output="sos"), x, axis=-1)


TRACE_ROOTS = [f"{REP}/traces", "/home/jovyan/tremorferometry/lfe_extension/traces"]


def load_native(netsta, year, doy):
    # Track D: production months live in lfe_extension/traces, the referee years in lfe_catalog_repair.
    p = next((q for q in (f"{r}/{netsta}/{year}/{doy:03d}.mseed" for r in TRACE_ROOTS)
              if os.path.exists(q) and os.path.getsize(q) > 0), f"{REP}/traces/{netsta}/{year}/{doy:03d}.mseed")
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return None
    st = read(p); st.merge(method=1, fill_value=0)
    chans = sorted({tr.stats.channel for tr in st})
    order = chanorder.resolve_order(chans, station=f"{netsta}/{year}.{doy}", strict=True)
    keep = [c for c in order if c is not None]
    npts = min(min(tr.stats.npts for tr in st.select(channel=c)) for c in keep)
    X = np.vstack([st.select(channel=c)[0].data[:npts].astype(np.float64) for c in keep])
    fs = float(st[0].stats.sampling_rate)
    t0 = float(st.select(channel=keep[0])[0].stats.starttime.timestamp)
    return X, fs, t0


def one2one(a, b_sorted, tol=2.0):
    a = np.asarray(a, float); out = np.full(len(a), np.inf)
    if len(a) == 0 or len(b_sorted) == 0:
        return out
    j = 0; nb = len(b_sorted)
    for i in np.argsort(a, kind="stable"):
        ti = a[i]
        while j < nb and b_sorted[j] < ti - tol:
            j += 1
        if j < nb and abs(b_sorted[j] - ti) <= tol:
            out[i] = abs(b_sorted[j] - ti); j += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--det", default=f"{IN}/e5r_detections_e6_mergeB.parquet")
    ap.add_argument("--tag", default="e6_mergeB")
    ap.add_argument("--L", type=float, default=8.0)
    ap.add_argument("--year", type=int, default=2012)
    ap.add_argument("--sta", default="PGC"); ap.add_argument("--netsta", default="CN.PGC")
    ap.add_argument("--hf", default="10,19")
    ap.add_argument("--xsta", default="B011"); ap.add_argument("--xnetsta", default="PB.B011")
    ap.add_argument("--xhf", default="10,20")
    ap.add_argument("--block-s", type=float, default=7200.0); ap.add_argument("--guard-s", type=float, default=60.0)
    ap.add_argument("--out", default=f"{OUT}/a1_hf_ratio")
    a = ap.parse_args()
    det = pd.read_parquet(a.det); det = det[det.L == a.L].reset_index(drop=True)
    F = pd.read_csv(f"{IN}/e5r_families_2012.csv").set_index("bkey")
    S = pd.read_csv(f"{IN}/e5r_stations.csv").set_index("sta")
    z = np.load(f"{IN}/e5r_templates_2012.npz", allow_pickle=True); sp = z["sp"].item()
    MV = pd.read_csv(f"{IN}/e4_moveout_e5r_mv_build.csv")
    mv = {(float(r.fam), r.sta, float(r.L)): float(r.lag_median) for r in MV.itertuples() if np.isfinite(r.lag_median)}
    res = {}
    for sta, ns, hfs in ((a.sta, a.netsta, a.hf), (a.xsta, a.xnetsta, a.xhf)):
        if not sta:
            continue
        f1, f2 = [float(x) for x in hfs.split(",")]
        off = {}
        for fk in det.fam.unique():
            fk = float(fk)
            if (fk, sta) not in sp:
                continue
            m = F.loc[int(fk)]
            dx = epi_km(float(m.lat), float(m.lon), float(S.loc[sta, "lat"]), float(S.loc[sta, "lon"]))
            off[fk] = float(v1d.tt(dx, float(m.dep), "P")) + float(sp[(fk, sta)]) - PRE + mv.get((fk, sta, a.L), 0.0)
        R = np.full(len(det), np.nan)
        for doy, g in det.groupby("doy"):
            L_ = load_native(ns, a.year, int(doy))
            if L_ is None:
                continue
            X, fs, t0 = L_
            xlo = bp(X, fs, 2.0, 8.0); xhi = bp(X, fs, f1, min(f2, 0.98 * fs / 2))
            elo = (xlo ** 2).sum(0); ehi = (xhi ** 2).sum(0)
            n = int(round(a.L * fs))
            ok = np.array([f in off for f in g.fam.to_numpy()])
            i0 = np.full(len(g), -1, np.int64)
            i0[ok] = np.rint((g.t.to_numpy()[ok] + np.array([off[f] for f in g.fam.to_numpy()[ok]]) - t0) * fs).astype(np.int64)
            good = ok & (i0 >= 0) & (i0 + n <= elo.shape[0])
            clo = np.concatenate([[0.0], np.cumsum(elo)]); chi = np.concatenate([[0.0], np.cumsum(ehi)])
            num = chi[i0[good] + n] - chi[i0[good]]; den = clo[i0[good] + n] - clo[i0[good]]
            R[g.index.to_numpy()[good]] = np.where(den > 0, num / np.maximum(den, 1e-30), np.nan)
            print(f"  [{sta} doy {doy}] {int(good.sum())}/{len(g)} windows, fs {fs:g}", flush=True)
        res[sta] = R
        det[f"hfr_{sta}"] = R
    det.to_parquet(f"{a.out}_{a.tag}.parquet", index=False)

    # ---- half / matched labels
    bb = np.floor(det.t.to_numpy() / a.block_s).astype(np.int64); pp = det.t.to_numpy() - bb * a.block_s
    is_test = (bb % 2 == 1) & (pp >= a.guard_s) & (pp <= a.block_s - a.guard_s)
    is_build = (bb % 2 == 0) & (pp >= a.guard_s) & (pp <= a.block_s - a.guard_s)
    gt = pd.read_parquet(f"{IN}/e5r_gt_bostock2012.parquet").reset_index(drop=True)
    test_ev = gt[(gt.half == "test") & (~gt.guard)].copy().reset_index(drop=True)
    rows = []
    for sta in res:
        Rv = det[f"hfr_{sta}"].to_numpy()
        cal = Rv[is_build & det.match_own.to_numpy() & np.isfinite(Rv)]
        cut = float(np.percentile(cal, 99)) if len(cal) else np.nan
        keep = ~(np.isfinite(Rv) & (Rv > cut))          # undefined ratio -> never vetoed
        for T in (7.5, 8.0, 8.5, 9.0):
            sel = det.net_mad.to_numpy() >= T
            q = sel & det.doy.isin(QUIET).to_numpy()
            fa0 = int(q.sum()); fa1 = int((q & keep).sum())
            rec = {}
            for lbl, kp in (("before", np.ones(len(det), bool)), ("after", keep)):
                D = det[sel & kp]
                hit = np.zeros(len(test_ev), bool)
                for fk, g in test_ev.groupby("bkey"):
                    td = np.sort(D[D.fam == float(fk)].t.to_numpy())
                    hit[test_ev.index.get_indexer(g.index)] = one2one(g.t_origin.to_numpy(), td) <= 2.0
                rec[lbl] = float(hit.mean())
            rows.append(dict(sta=sta, cut_p99=cut, n_calib=len(cal), T=T, quiet_det_before=fa0, quiet_det_after=fa1,
                             quiet_fa_reduction=1 - fa1 / max(fa0, 1), recall_before=rec["before"], recall_after=rec["after"],
                             recall_cost=rec["before"] - rec["after"],
                             recall_cost_rel=(rec["before"] - rec["after"]) / max(rec["before"], 1e-9),
                             adopt=bool((1 - fa1 / max(fa0, 1)) >= 0.30 and (rec["before"] - rec["after"]) / max(rec["before"], 1e-9) <= 0.02)))
    D = pd.DataFrame(rows); D.to_csv(f"{a.out}_{a.tag}_decision.csv", index=False)
    pd.set_option("display.width", 250); print(D.round(4).to_string(index=False))
    for sta in res:
        Rv = det[f"hfr_{sta}"]
        q = det.doy.isin(QUIET) & (det.net_mad >= 8.0)
        m = det.match_own & (det.net_mad >= 8.0)
        print(f"[{sta}] HF ratio medians: own-matched {Rv[m].median():.4f} (n {int(m.sum())}), "
              f"quiet-day {Rv[q].median():.4f} (n {int(q.sum())}), all {Rv.median():.4f}; defined {np.isfinite(Rv).mean():.1%}")
    print(f"[write] {a.out}_{a.tag}.parquet, {a.out}_{a.tag}_decision.csv")


if __name__ == "__main__":
    main()
