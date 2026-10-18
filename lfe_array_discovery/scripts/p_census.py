#!/usr/bin/env python
"""PER-STATION NOISE-REPEATER CENSUS on agent N's 2008 pair stage.  Prefix `p_`.

WHY.  The 51-station CAFE+PB Puget/S-Olympic 2008 bed returned a pair-stage tremor/quiet ratio
of **0.89** (`results/WA_2008_RESULTS.md` §3): the network autocorrelation fires as often on a
literature-quiet day as inside the audited ETS, so the pair population is dominated by repeating
cultural/site noise rather than by LFE pairs.  On the north (Vancouver Island POLARIS) the same
unchanged recipe gives 11x.  Royer & Bostock built "composite arrays" out of the stations that
worked; station QUALITY, not station count, is the lever this script measures.

WHAT IS MEASURED.  The pair stage stores only the NETWORK function's detections (i, j, cc, nmad),
not the per-station contributions, so the per-station cc is recomputed here at exactly the stored
(i, j) pairs, with k_auto's own window definition and its own degenerate-window QC:

  * windows are 6 s (n = 120 at 20 Hz); row i starts at sample i*hop (hop = 10), column j starts
    at sample j; both are the indices k_auto wrote into the parquet;
  * per component the cc is the normalised dot product of those two windows (k_auto's a @ H.T / E);
  * a station's contribution c_s = the MEAN of its three component cc's, which is exactly the
    term station s adds to the network mean;
  * the station's OWN bar is med_s + 8 * MAD_s of c_s over REF_PAIRS random (i, j) pairs of the
    same block with the same |lag| >= 6 s rule -- i.e. the 8 x MAD rule of the recipe applied to
    the single station instead of to the network.

Outputs per station, separately over the AUDITED TREMOR hours and the literature-QUIET hours:
  n_contrib          detections of the network where this station is above its own 8 x MAD bar
  contrib_per_hour   that count divided by the blocks the station took part in
  ratio              contrib_per_hour(tremor) / contrib_per_hour(quiet)      <-- the census number
  top3_share         fraction of the block's detections in which the station is among the TOP 3
                     contributors by its own nmad_s                          <-- the other one

A "noise repeater station" (the brief's definition) has ratio < 2 and a high top3_share on quiet
hours: it manufactures the same pair over and over whether or not there is tremor, and it is one
of the three stations that put the network function over threshold when it does.

  python scripts/p_census.py --hours tremor --auto results/n_puget2008/2008/auto \
      --out results/p_2008/census_tremor.parquet
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import n_common as L                                                        # noqa: E402

REF_PAIRS = 4000        # random (i, j) pairs per block for the per-station med/MAD reference
TOPK = 3                # "top-3 contributors"
SEED = 20080505


def block_cc(Db, det_i, det_j, fs, hop_s=0.5, win_s=6.0, min_lag_s=6.0, rng=None):
    """Per-station 3-C mean cc at the given (row i, sample j) pairs, plus the per-station
    reference median/MAD from random pairs of the same block.  Returns
    (C (P, S) float32, Z (P, S) float32 = (c - med)/mad, stas, med, mad)."""
    n = int(round(win_s * fs))
    hop = int(round(hop_s * fs))
    minlag = int(round(min_lag_s * fs))
    stas = sorted(Db)
    NH = min(x.shape[1] for x in Db.values())
    NV = NH - n + 1
    NW = (NH - n) // hop + 1
    det_i = np.asarray(det_i, np.int64)
    det_j = np.asarray(det_j, np.int64)
    ok = (det_i >= 0) & (det_i < NW) & (det_j >= 0) & (det_j < NV)
    det_i, det_j = det_i[ok], det_j[ok]
    # reference pairs: same |lag| >= one window rule, uniform over the block
    rng = np.random.default_rng(SEED) if rng is None else rng
    ri = rng.integers(0, NW, size=REF_PAIRS * 2)
    rj = rng.integers(0, NV, size=REF_PAIRS * 2)
    keep = np.abs(rj - ri * hop) >= minlag
    ri, rj = ri[keep][:REF_PAIRS], rj[keep][:REF_PAIRS]
    I = np.concatenate([det_i, ri])
    J = np.concatenate([det_j, rj])
    nd, nr = len(det_i), len(ri)
    off = np.arange(n, dtype=np.int64)
    Ii = (I * hop)[:, None] + off[None, :]                                  # (P, n)
    Jj = J[:, None] + off[None, :]
    C = np.full((len(I), len(stas)), np.nan, np.float32)
    for k, s in enumerate(stas):
        X = Db[s][:, :NH]
        a = X[:, Ii]                                                        # (3, P, n)
        b = X[:, Jj]
        ea = np.sqrt((a.astype(np.float64) ** 2).sum(-1))                   # (3, P)
        eb = np.sqrt((b.astype(np.float64) ** 2).sum(-1))
        # k_auto's degenerate-window QC, per component: a window whose energy is < 1e-4 of the
        # component's median window energy is a gap / dead channel and contributes nothing.
        emed_a = np.median(ea, axis=1, keepdims=True)
        emed_b = np.median(eb, axis=1, keepdims=True)
        good = (ea > 1e-4 * emed_a) & (eb > 1e-4 * emed_b)
        num = (a.astype(np.float64) * b.astype(np.float64)).sum(-1)
        cc = np.where(good, num / np.maximum(ea * eb, 1e-30), np.nan)       # (3, P)
        with np.errstate(invalid="ignore"):
            C[:, k] = np.nanmean(cc, axis=0).astype(np.float32)
    ref = C[nd:nd + nr]
    med = np.nanmedian(ref, axis=0)
    mad = np.nanmedian(np.abs(ref - med[None, :]), axis=0)
    mad = np.where(np.isfinite(mad) & (mad > 0), mad, np.nan)
    Z = (C[:nd] - med[None, :]) / mad[None, :]
    return C[:nd], Z, stas, med, mad, ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hours", default="tremor", choices=["tremor", "quiet"])
    ap.add_argument("--auto", required=True, help="agent N's per-hour pair-stage parquet dir")
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--stas", type=int, default=0, help="smoke: keep only the first N stations")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    K = L.activate()
    import k_auto as KA
    S = pd.read_csv(f"{a.auto}/_stats.csv")
    S = S[S.nmad_in_sd > 3]                                                 # k_group's block QC
    blocks = [(str(r.date), int(r.hour), str(r.stas)) for r in S.itertuples()]
    if a.limit:
        blocks = blocks[:a.limit]
    veto = L.health_veto()
    all_stas = L.stations()
    dc = K.DayCache(all_stas)
    rng = np.random.default_rng(SEED)
    rows = []
    L.log(f"[census {a.hours}] {len(blocks)} blocks, {len(all_stas)} admitted stations")
    for (date, hh, stas_str) in blocks:
        p = f"{a.auto}/{date}_{hh:02d}.parquet"
        if not os.path.exists(p):
            continue
        R = pd.read_parquet(p)
        want = [s for s in stas_str.split(",") if s]
        if a.stas:
            want = want[:a.stas]
        bad = veto.get((date, hh), set())
        D0 = dc.day(date)
        t_h = pd.Timestamp(f"{date} {hh:02d}:00:00").timestamp()
        Db = {}
        for s in want:
            if s in bad or s not in D0:
                continue
            t0, X = D0[s]
            i0 = int(round((t_h - t0) * K.FS))
            seg = X[:, i0:i0 + int(round(3600 * K.FS))]
            if seg.shape[1] < int(3000 * K.FS) or not np.isfinite(seg).all():
                continue
            Db[s] = seg
        if len(Db) < 4:
            continue
        C, Z, stas, med, mad, ok = block_cc(Db, R.i.to_numpy(), R.j.to_numpy(), K.FS,
                                            hop_s=KA.HOP_S, win_s=KA.WIN_S,
                                            min_lag_s=KA.MIN_LAG_S, rng=rng)
        nd = Z.shape[0]
        if nd == 0:
            continue
        above = Z >= 8.0                                                    # its own 8 x MAD
        above6 = Z >= 6.0
        above4 = Z >= 4.0                                                   # denser, same shape
        # top-3 contributors by the station's own nmad (NaN stations rank last)
        Zf = np.where(np.isfinite(Z), Z, -1e9)
        k = min(TOPK, Zf.shape[1])
        thr_k = np.partition(Zf, -k, axis=1)[:, -k][:, None]
        intop = (Zf >= thr_k) & np.isfinite(Z)
        for m, s in enumerate(stas):
            rows.append(dict(sta=s, date=date, hour=hh, n_det=int(nd),
                             n_above=int(np.nansum(above[:, m])),
                             n_above6=int(np.nansum(above6[:, m])),
                             n_above4=int(np.nansum(above4[:, m])),
                             n_top3=int(intop[:, m].sum()),
                             z_mean=float(np.nanmean(Z[:, m]))
                             if np.isfinite(Z[:, m]).any() else np.nan,
                             med_cc=float(med[m]) if np.isfinite(med[m]) else np.nan,
                             mad_cc=float(mad[m]) if np.isfinite(mad[m]) else np.nan,
                             z_p90=float(np.nanpercentile(Z[:, m], 90))
                             if np.isfinite(Z[:, m]).any() else np.nan))
        L.log(f"[{a.hours}] {date} {hh:02d}h  {len(stas)} sta  {nd} det  "
              f"above8mad/det={np.nansum(above) / max(nd, 1):.2f}")
    B = pd.DataFrame(rows)
    B.to_parquet(a.out, index=False)
    L.log(f"{len(B)} station-block rows -> {a.out}")


if __name__ == "__main__":
    main()
