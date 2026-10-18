#!/usr/bin/env python
"""Stage 5 -- DEDUPE.  Two templates are ONE family when

  (a) their independent HIT SETS have Jaccard overlap >= 0.5 (+- 2 s), or
  (b) their stacks correlate above the per-station CROSS-TEMPLATE p95 at >= 3 sites where both
      register.

Keep the larger (most independent hits).  Merlin's rule (PLAN sec.10); the p95 table is recomputed
from the bank itself, as g_crossfam did from the 2012 stacks.  `v_dedupe.cc_pair` /
`v_dedupe.crossfam_p95` are imported unchanged; only the overlap definition differs from P4-0b
(Jaccard, not containment) and the cc clause counts SITES, not stations.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w_common as W                                                        # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bank", required=True)
    ap.add_argument("--indep", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    V, E = W.patch_v()
    import v_dedupe as DD                                                   # noqa: E402
    mem, of, order, det = W.sites()
    tids, stas, Wt, off, meta = V.load_bank(a.bank)
    H = (pd.read_parquet(a.indep) if os.path.exists(a.indep)
         else pd.DataFrame(columns=W.HIT_COLS + ["tid"]))
    n = len(tids)
    if n == 0:
        pd.DataFrame(columns=W.DEDUPE_COLS).to_csv(a.out, index=False)
        pd.DataFrame(columns=["a", "b", "jaccard", "n_both", "n_hi_sites"]).to_csv(
            a.out.replace(".csv", "_pairs.csv"), index=False)
        W.log("[dedupe] empty bank")
        return
    pre = float(meta.get("pre_s", W.TPL1_PRE_S))
    ref = "PGC" if "PGC" in stas else stas[0]
    ri = stas.index(ref)
    T = {}
    for bi, tid in enumerate(tids):
        g = H[H.tid == tid]
        o = off[bi, ri]
        o = o if np.isfinite(o) else float(np.nanmean(off[bi]))
        T[tid] = np.sort(g.t.to_numpy(float) + o + pre) if len(g) else np.zeros(0)
    p95 = DD.crossfam_p95(Wt, stas)
    pd.Series(p95).to_csv(a.out.replace(".csv", "_crossfam_p95.csv"))
    lag = int(round(W.DEDUPE_LAG_S * W.FS))
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
            jac = 0.0
            if len(ti) and len(tj):
                m, _ = V.one2one(ti, tj, W.DEDUPE_TOL_S)
                k = int(m.sum())
                jac = k / max(len(ti) + len(tj) - k, 1)
            hi_sites = set()
            nboth = 0
            for si, s in enumerate(stas):
                if s not in p95 or s not in of:
                    continue
                if not (np.isfinite(Wt[i, si]).all() and np.isfinite(Wt[j, si]).all()):
                    continue
                nboth += 1
                c = abs(DD.cc_pair(Wt[i, si].astype(np.float64), Wt[j, si].astype(np.float64), lag))
                if c > p95[s]:
                    hi_sites.add(of[s])
            same = bool(jac >= W.DEDUPE_JACCARD or len(hi_sites) >= W.DEDUPE_CC_MIN_SITES)
            if same:
                parent[find(i)] = find(j)
                pairs.append(dict(a=tids[i], b=tids[j], jaccard=float(jac), n_both=nboth,
                                  n_hi_sites=len(hi_sites)))
    grp = {}
    for i in range(n):
        grp.setdefault(find(i), []).append(i)
    rows = []
    for g in grp.values():
        k = g[int(np.argmax([len(T[tids[i]]) for i in g]))]
        for i in g:
            rows.append(dict(tid=tids[i], group=tids[k], kept=bool(i == k),
                             n_indep=int(len(T[tids[i]]))))
    R = pd.DataFrame(rows)
    R.to_csv(a.out, index=False)
    pd.DataFrame(pairs).to_csv(a.out.replace(".csv", "_pairs.csv"), index=False)
    W.log(f"[dedupe] {n} templates -> {int(R.kept.sum())} distinct ({len(pairs)} merged pairs)")


if __name__ == "__main__":
    main()
