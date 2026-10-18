#!/usr/bin/env python
"""E2 step 2b (added 2026-09-10, resumed run) — EXCESS-based validation at PGC 2010.

Why: the +-3 s window recall/precision of e2_validate_pgc2010.py is dominated by chance
coincidence (900-1200 picks/day at thr 0.1 = one pick per 70-95 s), so the shifted-time null is
of the same order as the measurement.  The decisive quantity is the SHAPE of the
(nearest S pick - predicted arrival) residual distribution: a real detection population is a
narrow peak at 0 s sitting on the flat chance background, and its area is a chance-free estimate
of how many ground-truth events the picker actually found.

Method: for every ground-truth event take ALL S picks within +-30 s of the predicted arrival,
histogram the residuals (0.25 s bins), fit the background as the mean density in 10 < |d| < 30 s,
and integrate the excess over |d| <= w (w = 1.0 s default, and 0.5 / 2.0 s reported).
  recall_excess  = excess area / n_gt        (fraction of GT events the picker really found)
  precision_excess = excess area / n_picks   (fraction of picks that are real, upper bound: a pick
                     matched to several GT events is counted once per event, so this is an upper
                     bound on precision when the GT is denser than the picks)
Outputs results/e2_resid_hist.csv (the histograms) and results/e2_resid_excess.csv (the numbers).
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

MAXLAG = 40.0
BIN = 0.25
BG_EXCL = 5.0     # background = median bin count further than this from the peak
SWEEP = (0.1, 0.2, 0.3, 0.4, 0.5)
WIDTHS = (0.5, 1.0, 2.0)


def epoch(s):
    return pd.to_datetime(s).astype("datetime64[ns]").astype("int64").to_numpy() / 1e9


def pairs_within(gt, pk, maxlag):
    """all (gt, pick) residuals with |pick - gt| <= maxlag, via a sorted two-pointer sweep."""
    gt = np.sort(gt); pk = np.sort(pk)
    lo = np.searchsorted(pk, gt - maxlag)
    hi = np.searchsorted(pk, gt + maxlag)
    out = []
    for g, a, b in zip(gt, lo, hi):
        if b > a:
            out.append(pk[a:b] - g)
    return np.concatenate(out) if out else np.zeros(0)


def main():
    import e2_validate_pgc2010 as V
    pk = V.load_picks()
    bo, lin, C = V.gt_tables()
    linq = lin[(lin.residual < 0.5) & (lin.N >= 4)]
    gts = {
        "lin_hq": linq.pred_1d.dropna().to_numpy(),
        "lin_1d": lin.pred_1d.dropna().to_numpy(),
        "lin_near40": lin[lin.epi_PGC < 40].pred_1d.dropna().to_numpy(),
        "bostock_stack": bo.pred_stack.dropna().to_numpy(),
        "bostock_near40": bo[bo.epi_PGC < 40].pred_stack.dropna().to_numpy(),
    }
    edges = np.arange(-MAXLAG, MAXLAG + BIN / 2, BIN)
    ctr = 0.5 * (edges[1:] + edges[:-1])
    hrows, erows = [], []
    for model in sorted(pk.model.unique()):
        for thr in SWEEP:
            sub = pk[(pk.model == model) & (pk.phase == "S")
                     & (pk.thr_first <= thr) & (pk.prob >= thr)]
            t = sub.t.to_numpy()
            for gk, g in gts.items():
                d = pairs_within(g, t, MAXLAG)
                h, _ = np.histogram(d, bins=edges)
                # BUGFIX 2026-09-10 (resumed run): the excess used to be integrated around lag 0
                # with the background taken from 10 < |lag| < 30 s.  The Bostock predictions are
                # offset by +12.4 s (their catalogue time convention, see e2_pgc2010_offsets.csv),
                # so the real peak sat INSIDE that background band: the background was inflated and
                # the reported excess at 0 s came out NEGATIVE.  The excess is now integrated around
                # the measured PEAK, with the background taken >5 s away from it (median = robust to
                # the peak itself).  For Lin the peak is at -0.25 s, so its numbers are unchanged.
                k = int(np.argmax(h))
                peak_lag = float(ctr[k])
                bgm = np.abs(ctr - peak_lag) > BG_EXCL
                bgd = float(np.median(h[bgm]))
                sbg = float(h[bgm].std(ddof=1))
                for c, hh in zip(ctr, h):
                    hrows.append(dict(model=model, threshold=thr, gt=gk, lag_s=c, count=int(hh),
                                      bg_per_bin=bgd, peak_lag_s=peak_lag))
                for w in WIDTHS:
                    m = np.abs(ctr - peak_lag) <= w
                    exc = float(h[m].sum() - bgd * m.sum())
                    # sigma of the excess: Poisson on the peak + background estimate uncertainty
                    sig = float(np.sqrt(max(h[m].sum(), 1) + m.sum() ** 2 * sbg ** 2 / bgm.sum()))
                    erows.append(dict(model=model, threshold=thr, gt=gk, half_width_s=w,
                                      peak_lag_s=peak_lag,
                                      n_gt=len(g), n_picks=len(t),
                                      counts_in_peak=int(h[m].sum()),
                                      bg_expected=round(bgd * m.sum(), 1),
                                      excess=round(exc, 1), excess_sigma=round(sig, 1),
                                      snr=round(exc / max(sig, 1e-9), 2),
                                      recall_excess=exc / len(g),
                                      precision_excess=exc / max(len(t), 1)))
            print(f"[{model} thr {thr}] {len(t)} S picks", flush=True)
    pd.DataFrame(hrows).to_csv(f"{RES}/e2_resid_hist.csv", index=False)
    ex = pd.DataFrame(erows)
    ex.to_csv(f"{RES}/e2_resid_excess.csv", index=False)
    print(f"[write] {RES}/e2_resid_hist.csv, {RES}/e2_resid_excess.csv")
    show = ex[ex.half_width_s == 1.0]
    print(show.pivot_table(index=["gt", "threshold"], columns="model",
                           values=["peak_lag_s", "excess", "snr", "recall_excess",
                                   "precision_excess"]).round(4).to_string())


if __name__ == "__main__":
    main()
