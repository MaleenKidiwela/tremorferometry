#!/usr/bin/env python
"""P4-0b de-duplication (Merlin P4-2 step 5, second half).

Two templates are ONE family if
  (a) their declustered independent hit sets overlap >= 50 % (+-2 s), compared on a COMMON time
      axis (each template's origin times converted to its own predicted S at the reference
      station, so templates cut at different cells are comparable), or
  (b) their stacks correlate above the per-station CROSS-TEMPLATE p95 at >= half the sites where
      both register.  The p95 table is recomputed from the grown bank itself, exactly as
      g_crossfam.py recomputed it from the 2012 stacks.
Keep the larger (most independent hits).
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import v_common as V                                                        # noqa: E402

LAG_S = 1.0


def cc_pair(a, b, lag):
    """max over +-lag samples of the 3-C mean normalised cc of two unit-L2 templates."""
    n = a.shape[1]
    best = -2.0
    for k in range(-lag, lag + 1):
        if k >= 0:
            x, y = a[:, k:], b[:, :n - k] if k else b
        else:
            x, y = a[:, :n + k], b[:, -k:]
        if x.shape[1] < 10:
            continue
        num = (x * y).sum(1)
        den = np.sqrt((x ** 2).sum(1) * (y ** 2).sum(1))
        v = float(np.mean(np.where(den > 0, num / np.maximum(den, 1e-30), 0.0)))
        best = max(best, v)
    return best


def crossfam_p95(W, stas):
    """(p95 per station of |cc| over all template pairs, n pairs)."""
    lag = int(round(LAG_S * V.FS))
    out = {}
    for si, s in enumerate(stas):
        rows = [b for b in range(W.shape[0]) if np.isfinite(W[b, si]).all()]
        if len(rows) < 5:
            continue
        cc = []
        for i in range(len(rows)):
            for j in range(i + 1, len(rows)):
                cc.append(abs(cc_pair(W[rows[i], si].astype(np.float64),
                                      W[rows[j], si].astype(np.float64), lag)))
        if cc:
            out[s] = float(np.quantile(cc, 0.95))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bank", required=True)
    ap.add_argument("--indep", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    tids, stas, W, off, meta = V.load_bank(a.bank)
    H = pd.read_parquet(a.indep)
    pre = float(meta.get("pre_s", V.TPL1_PRE_S))
    ref = "PGC" if "PGC" in stas else stas[0]
    ri = stas.index(ref)
    T = {}
    for bi, tid in enumerate(tids):
        g = H[H.tid == tid]
        T[tid] = np.sort(g.t.to_numpy(float) + off[bi, ri] + pre) if len(g) else np.zeros(0)
    p95 = crossfam_p95(W, stas)
    pd.Series(p95).to_csv(a.out.replace(".csv", "_crossfam_p95.csv"))
    lag = int(round(LAG_S * V.FS))
    n = len(tids)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            ti, tj = T[tids[i]], T[tids[j]]
            ov = 0.0
            if len(ti) and len(tj):
                m, _ = V.one2one(ti, tj, V.DEDUPE_TOL_S)
                ov = m.sum() / min(len(ti), len(tj))
            nboth = nhi = 0
            for si, s in enumerate(stas):
                if s not in p95:
                    continue
                if not (np.isfinite(W[i, si]).all() and np.isfinite(W[j, si]).all()):
                    continue
                nboth += 1
                c = abs(cc_pair(W[i, si].astype(np.float64), W[j, si].astype(np.float64), lag))
                nhi += int(c > p95[s])
            same = bool(ov >= V.DEDUPE_OVERLAP or (nboth and nhi >= np.ceil(nboth / 2)))
            if same:
                parent[find(i)] = find(j)
                pairs.append(dict(a=tids[i], b=tids[j], overlap=float(ov),
                                  n_both=nboth, n_hi=nhi))
    grp = {}
    for i in range(n):
        grp.setdefault(find(i), []).append(i)
    rows = []
    for g in grp.values():
        sizes = [len(T[tids[i]]) for i in g]
        k = g[int(np.argmax(sizes))]
        for i in g:
            rows.append(dict(tid=tids[i], group=tids[k], kept=bool(i == k),
                             n_indep=int(len(T[tids[i]]))))
    R = pd.DataFrame(rows)
    R.to_csv(a.out, index=False)
    pd.DataFrame(pairs).to_csv(a.out.replace(".csv", "_pairs.csv"), index=False)
    V.log(f"[dedupe] {n} templates -> {R.kept.sum()} distinct ({len(pairs)} merged pairs)")


if __name__ == "__main__":
    main()
