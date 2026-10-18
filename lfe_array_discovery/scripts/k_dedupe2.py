#!/usr/bin/env python
"""Merge duplicate grown families using BOTH of Plourde et al. (2015) sec.2's uniqueness tests.

  test 1  "the number of shared detection times was computed, allowing for a lag ... As this
           method requires times to be shared TO THE NEAREST SAMPLE, the templates were considered
           to originate from the same family if >10% of detections were shared."
  test 2  "cross-correlation coefficients were computed for each template pair, using only
           stations for which a phase was picked ... A correlation coefficient >0.3 was considered
           to suggest templates originated from the same family."

Bostock's families sit 0.1-1.6 km apart and co-detect each other's events, so test 1 alone chains
neighbours into one blob.  Here test 1 only proposes an edge; test 2 (the waveform cc of the two
STACKS at a common lag over the stations they share) decides it.  The cc threshold is reported
against the development truth so the value used is a measured one, not an assumed 0.3.
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
import k_common as K                                                        # noqa: E402

SHARE = 0.10
TOL = 0.10
LAGMAX = 40.0
CCLAG_S = 3.0


def load_tpl(p):
    z = np.load(p, allow_pickle=True)
    meta = json.loads(str(z["meta"]))
    T = {}
    for k in z.files:
        if k == "meta":
            continue
        g, s = k.split("_", 1)
        T.setdefault(int(g), {})[s] = (z[k], meta[k][0])
    return T


def tpl_cc(Ta, Tb, fs=K.FS, lagmax=CCLAG_S):
    """Network-mean 3-C cc of two template sets at the single best COMMON lag."""
    sh = [s for s in Ta if s in Tb]
    if len(sh) < 3:
        return np.nan, 0
    nl = int(round(lagmax * fs))
    n = min(min(Ta[s][0].shape[1] for s in sh), min(Tb[s][0].shape[1] for s in sh))
    acc = np.zeros(2 * nl + 1)
    cnt = 0
    for s in sh:
        a, oa = Ta[s][0][:, :n], Ta[s][1]
        b, ob = Tb[s][0][:, :n], Tb[s][1]
        pad = np.zeros((3, n + 2 * nl))
        pad[:, nl:nl + n] = b
        V = np.lib.stride_tricks.sliding_window_view(pad, n, axis=1)[:, :2 * nl + 1, :]
        na = np.sqrt((a ** 2).sum())
        num = np.einsum("clt,ct->l", V, a)
        den = np.sqrt(np.einsum("clt,clt->l", V, V))
        if na <= 0:
            continue
        cc = np.where(den > 0, num / (na * np.maximum(den, 1e-30)), 0.0)
        # shift so the lag axis is the same PHYSICAL lag for every station
        k = int(round((ob - oa) * fs))
        cc = np.roll(cc, -k)
        acc += cc
        cnt += 1
    if cnt < 3:
        return np.nan, cnt
    return float((acc / cnt).max()), cnt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--det", required=True)
    ap.add_argument("--tpl", required=True)
    ap.add_argument("--score", default="", help="pre-merge score csv, for calibration only")
    ap.add_argument("--cc-thr", type=float, default=-1.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    D = pd.read_parquet(a.det)
    T = load_tpl(a.tpl)
    G = {int(g): np.sort(gg.t.to_numpy()) for g, gg in D.groupby("grp") if int(g) in T}
    gs = sorted(G, key=lambda g: -len(G[g]))
    K.log(f"{len(gs)} grown templates")

    # ---- test 1: candidate edges
    edges = []
    for ii, g in enumerate(gs):
        ta = G[g]
        for h in gs[ii + 1:]:
            tb = G[h]
            i = np.searchsorted(tb, ta)
            d = []
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
                sh = max(sh, int((np.abs(tb[j] - ts) <= TOL).sum()))
            f = sh / min(len(ta), len(tb))
            if f > SHARE:
                edges.append((g, h, f))
    K.log(f"test 1 proposes {len(edges)} edges")

    # ---- test 2: waveform cc of the stacks
    rows = []
    for g, h, f in edges:
        cc, nsh = tpl_cc(T[g], T[h])
        rows.append(dict(g=g, h=h, share=f, cc=cc, n_sh=nsh))
    E = pd.DataFrame(rows)
    E.to_csv(a.out.replace(".parquet", "_edges.csv"), index=False)

    thr = a.cc_thr
    if a.score and os.path.exists(a.score):
        S = pd.read_csv(a.score).set_index("grp")
        lab = {g: (int(S.bkey.get(g, -1)) if S.frac.get(g, 0) >= 0.30 else -1) for g in gs}
        E["same"] = [lab.get(r.g, -1) == lab.get(r.h, -2) and lab.get(r.g, -1) >= 0
                     for r in E.itertuples()]
        E["diff"] = [lab.get(r.g, -1) >= 0 and lab.get(r.h, -2) >= 0 and
                     lab.get(r.g) != lab.get(r.h) for r in E.itertuples()]
        A = E[E.same & np.isfinite(E.cc)].cc
        B = E[E["diff"] & np.isfinite(E.cc)].cc
        K.log(f"labelled edges: SAME family n={len(A)} cc med {A.median():.3f} "
              f"p10 {A.quantile(.1):.3f} | DIFFERENT n={len(B)} cc med {B.median():.3f} "
              f"p90 {B.quantile(.9):.3f} p95 {B.quantile(.95):.3f}")
        if thr < 0 and len(B) > 5:
            thr = float(B.quantile(0.95))
            K.log(f"cc threshold = cross-family p95 = {thr:.3f} "
                  f"(keeps {(A >= thr).mean():.0%} of same-family edges)")
    if thr < 0:
        thr = 0.30
    parent = {g: g for g in gs}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    nm = 0
    for r in E.itertuples():
        if np.isfinite(r.cc) and r.cc >= thr and find(r.g) != find(r.h):
            parent[find(r.h)] = find(r.g)
            nm += 1
    root = {g: find(g) for g in gs}
    surv = {}
    for g in gs:
        rt = root[g]
        if rt not in surv or len(G[g]) > len(G[surv[rt]]):
            surv[rt] = g
    keep = set(surv.values())
    K.log(f"{len(gs)} -> {len(keep)} families after BOTH tests (cc thr {thr:.3f}, {nm} merges)")
    D[D.grp.isin(keep)].to_parquet(a.out, index=False)


if __name__ == "__main__":
    main()
