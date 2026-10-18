#!/usr/bin/env python
"""Stage 2 of the literal recipe: group the autocorrelation detections into candidate families.

Bostock et al. (2012) grouped 5,775 autocorrelation detections into 140 families by "waveform
cluster analysis".  The most literal automatic reading: each detection is an edge between two
6-s windows that are the SAME waveform, so the windows are nodes, the detections are edges, and a
candidate family is a connected component of that graph.  Windows within one declustering interval
(4 s, Bostock 2015 §2) are the same node.

No minimum group size is stated in the papers; we report the whole size spectrum and carry every
component with >= --min-members (default 3, i.e. more than a bare pair).
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
import k_common as K                                                        # noqa: E402

MERGE_S = 4.0


def load_detections(d, stats=None):
    fs = sorted(glob.glob(f"{d}/*.parquet"))
    A = [pd.read_parquet(f) for f in fs]
    A = [x for x in A if len(x)]
    R = pd.concat(A, ignore_index=True)
    if "nmad" not in R.columns or R.nmad.isna().any():
        S = pd.read_csv(stats or f"{d}/_stats.csv")
        m = {f"{r.date}_{r.hour}": (r.med, r.mad) for r in S.itertuples()}
        k = (R.date + "_" + R.hour.astype(str)).to_numpy()
        med = np.array([m[x][0] for x in k])
        mad = np.array([m[x][1] for x in k])
        R["nmad"] = (R.cc.to_numpy() - med) / mad
    # a block whose 8xMAD bar is not ~5 sd is a data-quality failure, not a detection set
    S = pd.read_csv(stats or f"{d}/_stats.csv")
    bad = {(r.date, int(r.hour)) for r in S.itertuples() if not (r.nmad_in_sd > 3)}
    if bad:
        R = R[~pd.Series(list(zip(R.date, R.hour)), index=R.index).isin(bad)]
    return R.reset_index(drop=True)


class UF:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--det", default=f"{K.OUT}/auto_full")
    ap.add_argument("--detfile", default="")
    ap.add_argument("--budget", type=int, default=20000, help="strongest detections carried")
    ap.add_argument("--min-members", type=int, default=3)
    ap.add_argument("--stas", default="")
    ap.add_argument("--times", default=f"{K.OUT}/s1all_times.csv")
    ap.add_argument("--out", default=f"{K.OUT}/groups")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    R = (pd.read_parquet(a.detfile) if a.detfile else load_detections(a.det))
    K.log(f"{len(R)} clean detections; nmad {R.nmad.min():.2f}..{R.nmad.max():.2f}")
    R = R.nlargest(a.budget, "nmad").reset_index(drop=True)
    t = np.concatenate([R.t_i.to_numpy(), R.t_j.to_numpy()])
    o = np.argsort(t)
    ts = t[o]
    nid = np.zeros(len(t), np.int64)
    cur = 0
    nid[o[0]] = 0
    last = ts[0]
    for k in range(1, len(ts)):
        if ts[k] - last > MERGE_S:
            cur += 1
        nid[o[k]] = cur
        last = ts[k]
    nn = cur + 1
    ni = nid[:len(R)]
    nj = nid[len(R):]
    uf = UF(nn)
    for x, y in zip(ni, nj):
        uf.union(int(x), int(y))
    root = np.array([uf.find(i) for i in range(nn)])
    # node representative time = median of its member window times
    tn = pd.Series(t).groupby(nid).median().to_numpy()
    comp = pd.DataFrame(dict(node=np.arange(nn), root=root, t=tn))
    sz = comp.groupby("root").size()
    K.log(f"{nn} nodes, {sz.size} components; sizes: "
          f"max {sz.max()} p99 {np.quantile(sz, .99):.0f} median {sz.median():.0f}; "
          f">= {a.min_members}: {(sz >= a.min_members).sum()}")
    keep = set(sz[sz >= a.min_members].index)
    comp = comp[comp.root.isin(keep)].copy()
    ren = {r: i for i, r in enumerate(sorted(keep))}
    comp["grp"] = comp.root.map(ren)
    comp["date"] = pd.to_datetime(comp.t, unit="s").dt.strftime("%Y-%m-%d")
    # per-component detection strength
    edge = pd.DataFrame(dict(root=root[ni], nmad=R.nmad.to_numpy()))
    st = edge[edge.root.isin(keep)].groupby("root").nmad.agg(["size", "max", "median"])
    st.index = [ren[i] for i in st.index]
    G = comp.groupby("grp").agg(n_nodes=("node", "size"), t0=("t", "min"), t1=("t", "max"),
                                n_days=("date", "nunique"))
    G = G.join(st.rename(columns=dict(size="n_edges", max="nmad_max", median="nmad_med")))
    G = G.sort_values(["n_nodes", "nmad_max"], ascending=False)
    comp.to_csv(f"{a.out}_nodes.csv", index=False)
    G.to_csv(f"{a.out}_summary.csv")
    K.log("\n" + str(G.head(15)))

    # ---- purity against the kept families (development truth)
    if os.path.exists(a.times):
        C = pd.read_csv(a.times)
        C["t_corr"] = C.t_origin + C.dt
        DET = a.stas.split(",") if a.stas else K.CLEAN11
        TT = K.fam_tt(DET)
        rows = []
        for g, gg in comp.groupby("grp"):
            tn_ = np.sort(gg.t.to_numpy())
            best = (None, 0)
            for b, cc in C.groupby("bkey"):
                ttm = float(np.median([TT[(int(b), s)] for s in DET]))
                te = np.sort(cc.t_corr.to_numpy() + ttm)
                idx = np.searchsorted(te, tn_)
                d = np.full(len(tn_), np.inf)
                for k in (-1, 0):
                    j = np.clip(idx + k, 0, len(te) - 1)
                    d = np.minimum(d, np.abs(te[j] - tn_))
                nm = int((d <= 8.0).sum())
                if nm > best[1]:
                    best = (int(b), nm)
            rows.append(dict(grp=int(g), n=len(tn_), bkey=best[0], n_match=best[1],
                             purity=best[1] / len(tn_)))
        P = pd.DataFrame(rows)
        P.to_csv(f"{a.out}_purity.csv", index=False)
        pure = P[P.purity >= 0.5]
        K.log(f"\ngroups {len(P)}; purity>=0.5 {len(pure)} covering "
              f"{pure.bkey.nunique()} kept families; purity>=0.8 {(P.purity>=0.8).sum()}")
        K.log("\n" + str(P.sort_values("n", ascending=False).head(20)))
    json.dump(dict(budget=a.budget, min_members=a.min_members, merge_s=MERGE_S,
                   n_nodes=int(nn), n_components=int(sz.size),
                   n_kept=int(len(keep))), open(f"{a.out}_meta.json", "w"), indent=1)


if __name__ == "__main__":
    main()
