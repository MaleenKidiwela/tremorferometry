#!/usr/bin/env python
"""Score autocorrelation detections against the KEPT families (development truth).

A detection (t_i, t_j) is a TRUE PAIR of kept family k if there are two catalogued events of k,
at their Step-1 corrected times, with
    |(t_j - t_i) - (t_e2 - t_e1)| <= tol            (the inter-event lag invariant)
and |t_i - (t_e1 + tt_med(k))| <= anchor            (the window really sits on the first event).
Chance is measured by shifting the catalogue by +300 s and rescoring identically.
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import k_common as K                                                        # noqa: E402


def load_det(d):
    fs = sorted(glob.glob(f"{d}/*.parquet"))
    A = [pd.read_parquet(f) for f in fs]
    A = [x for x in A if len(x)]
    return pd.concat(A, ignore_index=True) if A else pd.DataFrame()


def score(R, C, TTM, tol=0.6, anchor=8.0, shift=0.0):
    """Return (matched detections DataFrame, per-family counts)."""
    out = []
    for b, g in C.groupby("bkey"):
        t = np.sort(g.t_corr.to_numpy() + shift)
        if len(t) < 2:
            continue
        ttm = TTM[b]
        # detections whose t_i sits on some event of this family
        ii = np.searchsorted(t, R.t_i.to_numpy() - ttm)
        cand = []
        for k in (-1, 0):
            idx = np.clip(ii + k, 0, len(t) - 1)
            cand.append(idx)
        best = None
        for idx in cand:
            d1 = np.abs(R.t_i.to_numpy() - ttm - t[idx])
            m = d1 <= anchor
            if not m.any():
                continue
            tj = R.t_j.to_numpy()[m]
            t1 = t[idx][m]
            # the partner event must reproduce the lag
            want = tj - R.t_i.to_numpy()[m] + t1
            jj = np.searchsorted(t, want)
            ok = np.zeros(len(want), bool)
            dmin = np.full(len(want), np.inf)
            for kk in (-1, 0):
                j2 = np.clip(jj + kk, 0, len(t) - 1)
                dd = np.abs(t[j2] - want)
                dmin = np.minimum(dmin, dd)
            ok = dmin <= tol
            if ok.any():
                sub = R[m].iloc[np.nonzero(ok)[0]].copy()
                sub["bkey"] = b
                out.append(sub)
    if not out:
        return pd.DataFrame(columns=list(R.columns) + ["bkey"])
    M = pd.concat(out, ignore_index=True)
    return M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--det", default=f"{K.OUT}/auto11")
    ap.add_argument("--times", default=f"{K.OUT}/s1_times.csv")
    ap.add_argument("--tol", type=float, default=0.6)
    ap.add_argument("--anchor", type=float, default=8.0)
    ap.add_argument("--top", type=int, default=0, help="score only the N strongest detections")
    ap.add_argument("--stas", default="")
    a = ap.parse_args()
    R = load_det(a.det)
    if not len(R):
        K.log("no detections")
        return
    if a.top:
        R = R.nlargest(a.top, "cc")
    C = pd.read_csv(a.times)
    C["t_corr"] = C.t_origin + C.dt
    stas = a.stas.split(",") if a.stas else K.CLEAN11
    TT = K.fam_tt(stas)
    TTM = {int(b): float(np.median([TT[(int(b), s)] for s in stas])) for b in C.bkey.unique()}
    K.log(f"{len(R)} detections, {C.bkey.nunique()} kept families with corrected times")
    M = score(R, C, TTM, a.tol, a.anchor)
    N = score(R, C, TTM, a.tol, a.anchor, shift=300.0)
    K.log(f"MATCHED {len(M)} detections ({len(M)/len(R):.3%}) over "
          f"{M.bkey.nunique() if len(M) else 0} families | "
          f"CHANCE(+300 s) {len(N)} ({len(N)/len(R):.3%}) over "
          f"{N.bkey.nunique() if len(N) else 0}")
    if len(M):
        g = M.groupby("bkey").agg(n=("cc", "size"), cc_max=("cc", "max"), cc_med=("cc", "median"))
        K.log("\nper family\n" + str(g.sort_values("n", ascending=False).head(25)))
        M.to_csv(f"{a.det}/_matched.csv", index=False)
        # enrichment in the strongest detections
        for q in (100, 500, 2000, 10000):
            if len(R) >= q:
                Rt = R.nlargest(q, "cc")
                Mt = score(Rt, C, TTM, a.tol, a.anchor)
                Nt = score(Rt, C, TTM, a.tol, a.anchor, shift=300.0)
                K.log(f"top {q:6d}: matched {len(Mt):5d} ({len(Mt)/q:.2%}) fams "
                      f"{Mt.bkey.nunique() if len(Mt) else 0:3d} | chance {len(Nt):5d} "
                      f"({len(Nt)/q:.2%})")


if __name__ == "__main__":
    main()
