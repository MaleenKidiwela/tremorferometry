#!/usr/bin/env python
"""E3 step 5 — PER-FAMILY Bostock time base (Merlin F7).

E2 measured ONE offset for the whole Slow-Earthquake-DB `Bostock2015-LFE` catalogue: the peak of
the (PGC S pick - Bostock time - ttS_PGC) residual histogram at +12.38 s.  The audit's F7 shows the
offset is station/distance dependent (+12.1..+13.1 s at 8 stations, PFB +9.4, LRIV +7.9; and
+12.1/+12.1/+7.1 s in PGC distance bins), which smears every Bostock-referenced number by ~5 s --
wider than the +-2 s scoring window.  If the offset is a property of the FAMILY (a reference
arrival at a per-family reference station is the obvious candidate) then a per-family table fixes
it exactly.

Method, per family with >= `--min-events` events in the pgc2010 window (2010-05-01..2010-11-01):
  residuals d = (S pick at the reference station) - (Bostock time + ttS_station), all pairs with
  |d| <= 40 s, histogrammed in 0.25 s bins; the offset is the argmax bin, the background is the
  median bin count further than 5 s from the peak, and the reported height is peak/background.
  A per-family offset is ACCEPTED only if peak/background >= `--min-peak` (default 3.0) and the
  Poisson significance of the peak over the background is >= 3 sigma.
The same is done per (family, station) over the 10-station vi2010 window so the offset can be
tested for being a source property (constant across stations) or a station property.

Outputs
  results/e3_bostock_family_offsets.csv   per family at the reference station (default PGC)
  results/e3_bostock_famsta_offsets.csv   per (family, station) over vi2010
  results/e3_bostock_offset_pattern.csv   the correlations against distance / depth / magnitude
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e2_resid_hist as RH  # noqa: E402

MAXLAG = 40.0
BIN = 0.25
BG_EXCL = 5.0
VI2010 = ["PGC", "SNB", "VGZ", "LZB", "PFB", "YOUB", "NLLB", "B011", "B926", "LRIV"]


def epoch(s):
    return pd.to_datetime(s).astype("datetime64[ns]").astype("int64").to_numpy() / 1e9


def peak_offset(d, maxlag=MAXLAG, half=1.0, step=BIN, bg_excl=BG_EXCL, nsigma=4.0):
    """Sliding-window peak of the residual distribution, with a Poisson significance.

    BUGFIX 2026-09-11 (E3 step 5, first run): the first version histogrammed in 0.25 s bins and
    accepted a family when peak/background >= 3.  Per family there are only ~100-350 Bostock events
    and ~0.8 PGC picks within +-40 s of each, i.e. ~100-300 pairs spread over 320 bins -> the
    background bin count is 0 for most families, `peak / max(bg, 1e-9)` returned 3e9, and a
    3-count noise bin was "accepted": it reported 96 of 100 families with offsets scattered over
    59 s.  Replaced by a 2 s-wide sliding window stepped 0.25 s (still 0.25 s offset resolution,
    but ~8x the counts per test), an expected count from the flat background OUTSIDE the peak, and
    a Poisson z = (n - mu)/sqrt(mu) with a 4 sigma bar -- 4 sigma because the window is scanned
    over ~320 positions per family, so a 3 sigma bar would admit ~0.4 false peaks per family.
    """
    out = dict(offset_s=np.nan, peak=0, bg_per_s=np.nan, mu=np.nan, peak_over_bg=np.nan,
               sigma=np.nan, n_pairs=int(len(d)), excess=np.nan, accept=False)
    if len(d) < 20:
        return out
    d = np.sort(d)
    cen = np.arange(-maxlag + half, maxlag - half + step / 2, step)
    lo = np.searchsorted(d, cen - half)
    hi = np.searchsorted(d, cen + half)
    n = hi - lo
    k = int(np.argmax(n))
    far = np.abs(cen - cen[k]) > bg_excl
    # background DENSITY from the residuals outside the peak neighbourhood (pairs per second)
    n_far = int(((d < cen[k] - bg_excl) | (d > cen[k] + bg_excl)).sum())
    span = 2 * maxlag - 2 * bg_excl
    bg_per_s = n_far / span if span > 0 else np.nan
    mu = bg_per_s * 2 * half
    z = (n[k] - mu) / np.sqrt(mu) if mu > 0 else np.nan
    out.update(offset_s=float(cen[k]), peak=int(n[k]), bg_per_s=float(bg_per_s), mu=float(mu),
               peak_over_bg=float(n[k] / mu) if mu > 0 else np.nan,
               sigma=float(z) if np.isfinite(z) else np.nan,
               excess=float(n[k] - mu),
               accept=bool(np.isfinite(z) and z >= nsigma and mu >= 1.0))
    return out


def load_picks(setname, model, thr, sta):
    p = f"{RES}/e2_picks_{setname}_{model}.parquet"
    df = pd.read_parquet(p)
    df = df[(df.phase == "S") & (df.sta == sta) & (df.thr_first <= thr) & (df.prob >= thr)]
    return np.sort(epoch(df.time))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="cascadia")
    ap.add_argument("--threshold", type=float, default=0.1)
    ap.add_argument("--ref-station", default="PGC")
    ap.add_argument("--min-events", type=int, default=50)
    ap.add_argument("--min-peak", type=float, default=3.0,
                    help="kept for the record; acceptance is now the 4 sigma Poisson test")
    args = ap.parse_args()

    bo = pd.read_parquet(f"{RES}/e2_gt_bostock.parquet")
    bo["t"] = epoch(bo.time)
    n_by_fam = bo.bkey.value_counts()
    fams = sorted(n_by_fam[n_by_fam >= args.min_events].index)
    print(f"[gt] {len(bo)} Bostock events in pgc2010, {bo.bkey.nunique()} families, "
          f"{len(fams)} with >= {args.min_events} events")

    # ---- per family at the reference station, over the whole pgc2010 window
    pk = load_picks("pgc2010", args.model, args.threshold, args.ref_station)
    print(f"[picks] {args.ref_station} pgc2010 {args.model} thr {args.threshold}: {len(pk)} S picks")
    ttcol = f"ttS_{args.ref_station}"
    rows = []
    for f in fams:
        g = bo[bo.bkey == f]
        pred = g.t.to_numpy() + g[ttcol].to_numpy()
        d = RH.pairs_within(pred, pk, MAXLAG)
        r = peak_offset(d)
        r.update(bkey=int(f), fam_id=float(g.fam_id.iloc[0]), n_events=len(g),
                 lat=float(g.lat.iloc[0]), lon=float(g.lon.iloc[0]), dep=float(g.dep.iloc[0]),
                 epi_ref=float(g[f"epi_{args.ref_station}"].iloc[0]),
                 mag_median=float(g.mag.median()), station=args.ref_station)
        rows.append(r)
    fam = pd.DataFrame(rows)
    fam.to_csv(f"{RES}/e3_bostock_family_offsets.csv", index=False)
    acc = fam[fam.accept]
    print(f"\n[per family @ {args.ref_station}] {len(acc)} of {len(fam)} accepted "
          f"(2 s window, Poisson z >= 4 over a flat background, mu >= 1)")
    if len(acc):
        print(f"  offset: median {acc.offset_s.median():+.2f} s, mean {acc.offset_s.mean():+.2f}, "
              f"sd {acc.offset_s.std():.2f}, min {acc.offset_s.min():+.2f}, "
              f"max {acc.offset_s.max():+.2f}, "
              f"IQR {acc.offset_s.quantile(.25):+.2f}..{acc.offset_s.quantile(.75):+.2f} "
              f"(range {acc.offset_s.max()-acc.offset_s.min():.1f} s)")
        print(f"  peak/bg: median {acc.peak_over_bg.median():.2f}, "
              f"sigma median {acc.sigma.median():.1f}")

    # ---- pattern: offset vs distance / depth / magnitude
    from scipy.stats import spearmanr, pearsonr
    prows = []
    for xc in ("epi_ref", "dep", "mag_median", "n_events", "lat", "lon"):
        if len(acc) >= 8:
            rs, ps = spearmanr(acc[xc], acc.offset_s)
            rp, pp = pearsonr(acc[xc], acc.offset_s)
            prows.append(dict(x=xc, n=len(acc), spearman=rs, p_spearman=ps,
                              pearson=rp, p_pearson=pp))
    pat = pd.DataFrame(prows)
    if len(pat):
        print("\n[pattern of the per-family offset]")
        print(pat.round(4).to_string(index=False))

    # ---- per (family, station) over vi2010 -> is the offset a SOURCE or a STATION property?
    rows2 = []
    for sta in VI2010:
        try:
            pk2 = load_picks("vi2010", args.model, args.threshold, sta)
        except Exception as e:
            print(f"[warn] {sta}: {e}")
            continue
        if len(pk2) == 0:
            continue
        t0, t1 = pk2.min() - 60, pk2.max() + 60
        sub = bo[(bo.t >= t0) & (bo.t <= t1)]
        # all events at this station (the station-level offset, F7's table)
        pred = sub.t.to_numpy() + sub[f"ttS_{sta}"].to_numpy()
        r = peak_offset(RH.pairs_within(pred, pk2, MAXLAG))
        r.update(station=sta, bkey=-1, n_events=len(sub))
        rows2.append(r)
        for f in fams:
            g = sub[sub.bkey == f]
            if len(g) < 20:                     # vi2010 is only 34 d: relax the per-family floor
                continue
            pred = g.t.to_numpy() + g[f"ttS_{sta}"].to_numpy()
            r = peak_offset(RH.pairs_within(pred, pk2, MAXLAG))
            r.update(station=sta, bkey=int(f), n_events=len(g),
                     epi=float(g[f"epi_{sta}"].iloc[0]), dep=float(g.dep.iloc[0]))
            rows2.append(r)
    fs = pd.DataFrame(rows2)
    fs.to_csv(f"{RES}/e3_bostock_famsta_offsets.csv", index=False)
    st = fs[fs.bkey == -1]
    print("\n[per station, all families, vi2010 window] (F7's table, recomputed)")
    print(st[["station", "n_events", "n_pairs", "offset_s", "peak_over_bg", "sigma"]]
          .round(2).to_string(index=False))
    fsa = fs[(fs.bkey >= 0) & fs.accept.fillna(False)]
    if len(fsa) >= 10:
        # variance decomposition: is the offset explained by the family or by the station?
        gm = fsa.offset_s.mean()
        ss_tot = float(((fsa.offset_s - gm) ** 2).sum())
        ss_fam = float(sum(len(g) * (g.offset_s.mean() - gm) ** 2
                           for _, g in fsa.groupby("bkey")))
        ss_sta = float(sum(len(g) * (g.offset_s.mean() - gm) ** 2
                           for _, g in fsa.groupby("station")))
        print(f"\n[family vs station] {len(fsa)} accepted (family, station) offsets, "
              f"{fsa.bkey.nunique()} families x {fsa.station.nunique()} stations; "
              f"between-FAMILY SS/total = {ss_fam/ss_tot:.3f}, "
              f"between-STATION SS/total = {ss_sta/ss_tot:.3f}, "
              f"sd {fsa.offset_s.std():.2f} s")
        prows.append(dict(x="ss_family_frac", n=len(fsa), spearman=ss_fam / ss_tot,
                          p_spearman=np.nan, pearson=np.nan, p_pearson=np.nan))
        prows.append(dict(x="ss_station_frac", n=len(fsa), spearman=ss_sta / ss_tot,
                          p_spearman=np.nan, pearson=np.nan, p_pearson=np.nan))
    pd.DataFrame(prows).to_csv(f"{RES}/e3_bostock_offset_pattern.csv", index=False)
    print(f"\n[write] {RES}/e3_bostock_family_offsets.csv, e3_bostock_famsta_offsets.csv, "
          f"e3_bostock_offset_pattern.csv")


if __name__ == "__main__":
    main()
