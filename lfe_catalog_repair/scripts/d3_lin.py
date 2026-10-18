#!/usr/bin/env python
"""D3 stage 3c — GROUND TRUTH: PGC re-detections vs the Lin/Bostock LFE catalogue.

Ground truth on disk: `data/detections_lin_ets_2010_vi.parquet` = 106,644 Lin LFE times in 217
families, 2010-05-17..2010-10-14, southern Vancouver Island (from Zenodo 10016020, clustered into
0.1-deg families by src/tremorferometry/lin_catalog.py). Family locations:
`catalogs/lin_families_ets_2010_vi.csv`.

Lin times are LFE ORIGIN times; a re-detection time is the START of the 2 s template window at the
STATION, so the two differ by (S travel time + 1 s). The lag is therefore ESTIMATED, and to keep the
recall/precision honest it is estimated on a HALF of the Lin events (even index) and evaluated on
the OTHER half (odd index).

Outputs results/d3_lin_groundtruth.csv (per PGC family x rule) and a printed summary.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

LIN_DET = "data/detections_lin_ets_2010_vi.parquet"
LIN_FAM = "catalogs/lin_families_ets_2010_vi.csv"
KS = (8, 9, 10)
TOL = 2.0            # +-2 s match window (the brief)
# The project already calibrated this offset: `archive/lfe_features_thread/lfe_features/
# refine_lin_picks.py` uses ANCHOR = 11.0 s, i.e. the S arrival at B011 (0.3 km from PGC) sits at
# Lin's origin time + ~11 s (+-3 s jitter, refined there by cross-correlation). A re-detection time
# is the START of the [-1,+1) s template window, so the expected lag here is ~10 s. The scan is
# therefore restricted to a physically plausible band instead of a blind +-40 s.
LAG_LO, LAG_HI = -5.0, 35.0
ANCHOR_PROJECT = 11.0 - 1.0


def km(la1, lo1, la2, lo2):
    return np.hypot((la1 - la2) * 111.2, (lo1 - lo2) * 111.2 * np.cos(np.radians((la1 + la2) / 2)))


def est_lag(t_det, t_lin, lo=LAG_LO, hi=LAG_HI, binsec=0.25):
    """Modal offset of nearest re-detection after each Lin time, from a histogram of differences."""
    if len(t_det) < 5 or len(t_lin) < 5:
        return np.nan, 0
    d = []
    s = np.sort(t_det)
    for tl in t_lin:
        i0 = np.searchsorted(s, tl + lo)
        i1 = np.searchsorted(s, tl + hi)
        if i1 > i0:
            d.append(s[i0:i1] - tl)
    if not d:
        return np.nan, 0
    d = np.concatenate(d)
    bins = np.arange(lo, hi + binsec, binsec)
    h, _ = np.histogram(d, bins=bins)
    if h.max() < 3:
        return np.nan, int(h.max())
    return float(bins[h.argmax()] + binsec / 2), int(h.max())


def match(t_det, t_lin, lag, tol=TOL):
    """Returns (recall, precision, n_lin, n_det)."""
    if len(t_det) == 0 or len(t_lin) == 0 or not np.isfinite(lag):
        return np.nan, np.nan, len(t_lin), len(t_det)
    s = np.sort(t_det) - lag
    x = np.sort(t_lin)
    i = np.searchsorted(s, x)
    d = np.full(len(x), np.inf)
    ok = i < len(s)
    d[ok] = np.abs(s[i[ok]] - x[ok])
    ok2 = i > 0
    d[ok2] = np.minimum(d[ok2], np.abs(s[i[ok2] - 1] - x[ok2]))
    recall = float((d <= tol).mean())
    j = np.searchsorted(x, s)
    e = np.full(len(s), np.inf)
    ok = j < len(x)
    e[ok] = np.abs(x[j[ok]] - s[ok])
    ok2 = j > 0
    e[ok2] = np.minimum(e[ok2], np.abs(x[j[ok2] - 1] - s[ok2]))
    prec = float((e <= tol).mean())
    return recall, prec, len(x), len(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sta", default="PGC")
    ap.add_argument("--det", required=True, help="results/d3_detections_PGC2010.parquet")
    ap.add_argument("--max-km", type=float, default=5.0)
    ap.add_argument("--min-det", type=int, default=5)
    ap.add_argument("--out", default="lfe_catalog_repair/results/d3_lin_groundtruth.csv")
    ap.add_argument("--rules", default="", help="comma-separated rule columns (default old + new8/9/10)")
    a = ap.parse_args()

    lin = pd.read_parquet(LIN_DET)
    lin["time"] = pd.to_datetime(lin["time"])
    # NB: the parquet stores datetime64[us]; .astype("int64") therefore yields MICROseconds, so the
    # unit must be pinned to ns before dividing (this silently produced epoch ~1.27e6 and an empty
    # intersection with the detections).
    lin["ep"] = lin["time"].astype("datetime64[ns]").astype("int64") / 1e9
    lf = pd.read_csv(LIN_FAM).set_index("family_id")
    det = pd.read_parquet(a.det)
    det["fam"] = det["template"].astype(str)
    fams = sorted(det.fam.unique())
    print(f"[lin] {len(lin):,} Lin events / {lin.family_id.nunique()} Lin families "
          f"{lin.time.min().date()}..{lin.time.max().date()}; {len(fams)} PGC templates", flush=True)
    # family -> nearest Lin family
    cell = np.array([[float(f.split('__')[0].split('_')[0]), float(f.split('__')[0].split('_')[1])]
                     for f in fams])
    D = km(cell[:, 0:1], cell[:, 1:2], lf.lat.values[None, :], lf.lon.values[None, :])
    near = lf.index.values[D.argmin(1)]
    dkm = D.min(1)
    keep = dkm <= a.max_km
    print(f"  {keep.sum()}/{len(fams)} PGC templates within {a.max_km} km of a Lin family "
          f"({len(set(near[keep]))} distinct Lin families)", flush=True)

    rules = ([r.strip() for r in a.rules.split(",") if r.strip()] if a.rules
             else ["old"] + [f"new{kk}" for kk in KS])
    rules = [r for r in rules if r in det.columns]
    print(f"  rules present: {rules}", flush=True)
    t0, t1 = lin.ep.min(), lin.ep.max()
    lin_by_fam = {k: np.sort(v.ep.values) for k, v in lin.groupby("family_id")}
    det = det[(det.t >= t0 - 60) & (det.t <= t1 + 60)]
    det_by_fam = {k: v for k, v in det.groupby("fam", sort=False)}
    rows = []
    for f, lfam, d0, k in zip(fams, near, dkm, keep):
        if not k:
            continue
        tl_all = lin_by_fam.get(lfam)
        if tl_all is None or len(tl_all) < 10:
            continue
        tl_fit, tl_ev = tl_all[0::2], tl_all[1::2]
        sub = det_by_fam.get(f)
        if sub is None:
            continue
        for rule in rules:
            s = sub[sub[rule]] if rule in sub else sub
            td = s["t"].values
            if len(td) < a.min_det:
                rows.append(dict(template=f, lin_family=lfam, dist_km=d0, rule=rule,
                                 n_det=len(td), n_lin=len(tl_all), lag=np.nan,
                                 recall=np.nan, precision=np.nan, peak=0))
                continue
            lag, peak = est_lag(td, tl_fit)
            rec, prec, nl, nd = match(td, tl_ev, lag)
            rows.append(dict(template=f, lin_family=lfam, dist_km=d0, rule=rule,
                             n_det=nd, n_lin=len(tl_all), lag=lag, recall=rec,
                             precision=prec, peak=peak))
    out = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, index=False)

    # ---------------- POOLED test (the per-template rates are far too low) ----------------
    # Per template the new rules yield O(1-30) detections over the whole window, so a per-template
    # recall/precision is dominated by counting noise. The pooled test uses ONE lag, estimated from
    # the even-indexed Lin events of every Lin-proximate template together, and scores the
    # odd-indexed ones.
    print("\n POOLED over all Lin-proximate templates (one global lag, fit on the even-index Lin "
          "events, scored on the odd-index ones):")
    print(" rule     n_det   n_lin    lag_s  peak/flat   recall   precision  recall_shift600  prec_shift600")
    pooled = []
    for rule in rules:
        dall, lfit, lev = [], [], []
        for f, lfam, d0, k in zip(fams, near, dkm, keep):
            if not k:
                continue
            tl = lin_by_fam.get(lfam)
            sub = det_by_fam.get(f)
            if tl is None or len(tl) < 10 or sub is None:
                continue
            s_ = sub[sub[rule]] if rule in sub else sub
            if not len(s_):
                continue
            dall.append(s_["t"].values)
            lfit.append(tl[0::2])
            lev.append(tl[1::2])
        if not dall:
            print(f" {rule:7s}  (no detections)")
            continue
        td = np.concatenate(dall)
        tf = np.concatenate(lfit)
        te = np.concatenate(lev)
        lag, peak = est_lag(td, tf)
        # flatness reference: mean histogram height away from the peak
        d = []
        ss = np.sort(td)
        for tl_ in tf:
            i0 = np.searchsorted(ss, tl_ + LAG_LO)
            i1 = np.searchsorted(ss, tl_ + LAG_HI)
            if i1 > i0:
                d.append(ss[i0:i1] - tl_)
        h, _ = np.histogram(np.concatenate(d) if d else np.zeros(1),
                            bins=np.arange(LAG_LO, LAG_HI + 0.25, 0.25))
        flat = (h.sum() - h.max()) / max(len(h) - 1, 1)
        rec, prec, nl, nd = match(td, te, lag)
        recs, precs, _, _ = match(td + 600.0, te, lag)
        pooled.append(dict(rule=rule, n_det=len(td), n_lin=len(te), lag=lag,
                           peak_over_flat=(peak / flat if flat > 0 else np.nan),
                           recall=rec, precision=prec,
                           recall_shift=recs, precision_shift=precs))
        print(f" {rule:7s} {len(td):8,} {len(te):7,}  {lag:+7.2f}  {peak/flat if flat>0 else np.nan:9.2f}"
              f"  {rec:8.4f}   {prec:8.4f}       {recs:10.4f}     {precs:10.4f}")
    if pooled:
        pd.DataFrame(pooled).to_csv(str(a.out).replace(".csv", "_pooled.csv"), index=False)
    print("\n rule    n_tmpl  lag_med   recall_med  recall_p90  prec_med  prec_p90  "
          "frac_rec>0.1  frac_prec>0.5")
    for rule, g in out.groupby("rule"):
        gg = g.dropna(subset=["recall"])
        if not len(gg):
            print(f" {rule:7s}  {len(g):5d}   (no template reached the detection minimum)")
            continue
        print(f" {rule:7s}  {len(gg):5d}  {gg.lag.median():+7.2f}  {gg.recall.median():10.4f}  "
              f"{gg.recall.quantile(0.9):10.4f}  {gg.precision.median():8.4f}  "
              f"{gg.precision.quantile(0.9):8.4f}  {(gg.recall > 0.1).mean():12.3f}  "
              f"{(gg.precision > 0.5).mean():13.3f}")
    gg = out.dropna(subset=["lag"])
    if len(gg):
        print(f"\n lag vs the project's calibrated offset ({ANCHOR_PROJECT:.1f} s = ANCHOR 11 s "
              f"minus the 1 s window offset): median {gg.lag.median():+.2f} s, "
              f"IQR {gg.lag.quantile(.25):+.2f}..{gg.lag.quantile(.75):+.2f}, "
              f"frac within +-3 s of it {(np.abs(gg.lag - ANCHOR_PROJECT) <= 3).mean():.3f}")
    print(f"  -> {a.out}")


if __name__ == "__main__":
    main()
