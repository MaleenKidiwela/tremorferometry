#!/usr/bin/env python
"""Merge duplicate grown families the way Plourde et al. (2015) sec.2 describe:

  "the number of shared detection times was computed ... the templates were considered to
   originate from the same family if >10% of detections were shared"

Shared = within +- 2 s (our sample grid is 0.05 s; Plourde allowed a lag because their template was
60 s long, ours is 12 s and the detection axis is the template's own reference).  Single-linkage
merge; the survivor of a merged set is the template with the most detections.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import k_common as K                                                        # noqa: E402

TOL = 0.10   # Plourde et al. 2015: shared detection times must agree "to the nearest
             # sample (dt = 0.025 s)".  Ours is 0.05 s at 20 Hz; 0.10 s = +-1 sample.
SHARE = 0.10
LAGMAX = 40.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--det", required=True)
    ap.add_argument("--tol", type=float, default=TOL)
    ap.add_argument("--share", type=float, default=SHARE)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    D = pd.read_parquet(a.det)
    T = {int(g): np.sort(gg.t.to_numpy()) for g, gg in D.groupby("grp")}
    gs = sorted(T, key=lambda g: -len(T[g]))
    parent = {g: g for g in gs}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for ii, g in enumerate(gs):
        for h in gs[ii + 1:]:
            if find(g) == find(h):
                continue
            ta, tb = T[g], T[h]
            # best constant lag between the two detection axes, then shared fraction
            d = []
            i = np.searchsorted(tb, ta)
            for k in (-1, 0):
                j = np.clip(i + k, 0, len(tb) - 1)
                d.append(ta - tb[j])
            d = np.concatenate(d)
            d = d[np.abs(d) <= LAGMAX]
            if d.size < 2:
                continue
            hh, e = np.histogram(d, bins=np.arange(-LAGMAX, LAGMAX + 0.05, 0.05))
            off = float((e[hh.argmax()] + e[hh.argmax() + 1]) / 2)
            ts = ta - off
            i = np.searchsorted(tb, ts)
            sh = 0
            for k in (-1, 0):
                j = np.clip(i + k, 0, len(tb) - 1)
                sh = max(sh, int((np.abs(tb[j] - ts) <= a.tol).sum()))
            if sh / min(len(ta), len(tb)) > a.share:
                parent[find(h)] = find(g)
    root = {g: find(g) for g in gs}
    surv = {}
    for g in gs:
        r = root[g]
        if r not in surv or len(T[g]) > len(T[surv[r]]):
            surv[r] = g
    keep = set(surv.values())
    K.log(f"{len(gs)} grown families -> {len(keep)} after the >10% shared-detection merge")
    D[D.grp.isin(keep)].to_parquet(a.out, index=False)
    pd.DataFrame([dict(grp=g, root=root[g], survivor=surv[root[g]]) for g in gs]).to_csv(
        a.out.replace(".parquet", "_map.csv"), index=False)


if __name__ == "__main__":
    main()
