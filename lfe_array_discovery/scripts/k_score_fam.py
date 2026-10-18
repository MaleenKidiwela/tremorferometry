#!/usr/bin/env python
"""Score grown families against the 100 KEPT families with BOSTOCK'S OWN criterion.

A grown family G is "kept family k" if >= --min-frac (0.30) of G's detections match catalogued
events of k ONE-TO-ONE at +- 2 s, after removing the constant offset between G's detection axis and
k's origin-time axis (estimated from the difference histogram; both axes are arbitrary origins).
Recall of k is reported in the three denominators: all catalogued events in the bed, the
60-s-independent subset, and the loud subset (top-decile network S-energy excess).
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

TOL = 2.0
OFF_MAX = 40.0
BIN = 0.25


def best_offset(td, te):
    """Mode of (detection - event) differences within +- OFF_MAX, on a BIN grid."""
    i = np.searchsorted(te, td)
    d = []
    for k in (-2, -1, 0, 1):
        j = np.clip(i + k, 0, len(te) - 1)
        d.append(td - te[j])
    d = np.concatenate(d)
    d = d[np.abs(d) <= OFF_MAX]
    if d.size < 3:
        return 0.0, 0
    h, e = np.histogram(d, bins=np.arange(-OFF_MAX, OFF_MAX + BIN, BIN))
    k = int(np.argmax(h))
    return float((e[k] + e[k + 1]) / 2), int(h[k])


def one_to_one(td, te, tol=TOL):
    """Greedy nearest one-to-one matching.  Returns (n_match, matched detection mask, ev idx)."""
    pairs = []
    i = np.searchsorted(te, td)
    for a in range(len(td)):
        for k in (-1, 0):
            j = int(np.clip(i[a] + k, 0, len(te) - 1))
            dd = abs(te[j] - td[a])
            if dd <= tol:
                pairs.append((dd, a, j))
    pairs.sort()
    ud, ue = set(), set()
    md, me = [], []
    for dd, a, j in pairs:
        if a in ud or j in ue:
            continue
        ud.add(a)
        ue.add(j)
        md.append(a)
        me.append(j)
    return len(md), np.array(md, int), np.array(me, int)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--det", required=True)
    ap.add_argument("--times", default=f"{K.OUT}/s1all_times.csv")
    ap.add_argument("--loud", default=f"{K.OUT}/s1all_loud.csv")
    ap.add_argument("--min-frac", type=float, default=0.30)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    D = pd.read_parquet(a.det) if a.det.endswith(".parquet") else pd.read_csv(a.det)
    C = pd.read_csv(a.times)
    C["t_corr"] = C.t_origin + C.dt
    L = pd.read_csv(a.loud)
    loud = set(zip(L[L.loud].bkey, L[L.loud].t_origin))
    fams = sorted(C.bkey.unique())
    EV, EV60, EVL = {}, {}, {}
    for b, g in C.groupby("bkey"):
        t = np.sort(g.t_corr.to_numpy())
        EV[b] = t
        keep = np.concatenate([[True], np.diff(t) > 60.0])
        EV60[b] = t[keep]
        EVL[b] = np.sort(g[[(int(r.bkey), float(r.t_origin)) in loud
                            for r in g.itertuples()]].t_corr.to_numpy())
    rows = []
    for g, gg in D.groupby("grp"):
        td = np.sort(gg.t.to_numpy())
        best = None
        for b in fams:
            te = EV[b]
            off, hh = best_offset(td, te)
            if hh < 2:
                continue
            n, md, me = one_to_one(td - off, te)
            if best is None or n > best[1]:
                best = (b, n, off)
        if best is None:
            rows.append(dict(grp=int(g), n_det=len(td), bkey=-1, n_match=0, frac=0.0))
            continue
        b, n, off = best
        n60, _, _ = one_to_one(td - off, EV60[b])
        nl, _, _ = one_to_one(td - off, EVL[b]) if len(EVL[b]) else (0, None, None)
        rows.append(dict(grp=int(g), n_det=len(td), bkey=int(b), n_match=int(n),
                         frac=n / len(td), offset_s=off,
                         recall_all=n / len(EV[b]),
                         recall_60=n60 / max(len(EV60[b]), 1),
                         recall_loud=nl / max(len(EVL[b]), 1),
                         n_ev=len(EV[b]), n_ev60=len(EV60[b]), n_evloud=len(EVL[b])))
    R = pd.DataFrame(rows)
    ok = R[R.frac >= a.min_frac]
    K.log(f"{len(R)} grown groups; {len(ok)} meet Bostock's >= {a.min_frac:.0%} one-to-one bar; "
          f"**{ok.bkey.nunique()} DISTINCT kept families rediscovered** of {len(fams)}")
    if len(ok):
        K.log("\n" + str(ok.sort_values("n_match", ascending=False)
                         [["grp", "n_det", "bkey", "n_match", "frac", "recall_all",
                           "recall_60", "recall_loud"]].head(30).round(3)))
        per = ok.sort_values("n_match", ascending=False).drop_duplicates("bkey")
        K.log(f"\nper-family medians: recall_all {per.recall_all.median():.3f} "
              f"recall_60 {per.recall_60.median():.3f} recall_loud {per.recall_loud.median():.3f} "
              f"purity {per.frac.median():.3f}")
    if a.out:
        R.to_csv(a.out, index=False)


if __name__ == "__main__":
    main()
