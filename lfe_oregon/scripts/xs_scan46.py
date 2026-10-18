#!/usr/bin/env python
"""Templates that pass the >=100-detection rule -> 46-day network cross-correlation -> 10 %/5 s
duplicate elimination (Thomas & Bostock 2015 sec.2):

  "LFE templates that make the cut are cross-correlated through 46 days of continuous data from
   August 3 through September 17, 2009.  We eliminate duplicate families by ensuring that no family
   shares more than 10% of detection times, defined as a difference between detection times of
   <5s, with any other family."

"Eliminate" is read literally: families are ranked by 46-day detection count; walking down the
ranking, a family is DROPPED if it shares > 10 % of its detections (within 5 s) with an already kept
family (the smaller family's fraction, as ab_grow.merge_duplicates defines it).  --mode merge uses
ab_grow.merge_duplicates (union of detections into the representative) as an ablation column.
Detection settings are the growth run's (args.json).

  python scripts/xs_scan46.py --grow results/grow --out results/scan46 [--dates ...]
-> <out>/det46.parquet (all passing templates), families.csv, det46_families.parquet, dup_pairs.csv
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re

import numpy as np
import pandas as pd

import xs_common as C

A = C.patch_ab_common()
import ab_grow as G                                                          # noqa: E402

POOL = [s for s in C.STATIONS.netsta if os.path.isdir(f"{C.CACHE}/{s}")]


def last_iter(gdir):
    its = sorted(int(re.findall(r"state_it(\d+)", p)[0]) for p in glob.glob(f"{gdir}/state_it*.csv"))
    return its[-1]


def load_tpl(npz):
    z = np.load(npz)
    T = {}
    for k in z.files:
        g, s = k.split("|")
        T.setdefault(int(g), {})[s] = (z[k].astype(np.float64), 0.0, np.nan)
    return T


def shared_fraction(det, win=5.0):
    """{(a, b): fraction of the SMALLER family's detections with one of the other within win}."""
    det = det.sort_values("t")
    t = det.t.to_numpy()
    g = det.grp.to_numpy()
    n = det.groupby("grp").size()
    hits = {}
    j0 = 0
    for i in range(len(t)):
        while t[i] - t[j0] >= win:
            j0 += 1
        for j in range(j0, i):
            if g[j] != g[i]:
                hits.setdefault((g[i], g[j]), set()).add(i)
                hits.setdefault((g[j], g[i]), set()).add(j)
    rows = []
    for (a, b), s in hits.items():
        rows.append(dict(a=a, b=b, n_a=int(n[a]), n_b=int(n[b]), shared_a=len(s), frac_a=len(s) / n[a]))
    return pd.DataFrame(rows)


def eliminate(det, frac=0.10, win=5.0):
    P = shared_fraction(det, win)
    n = det.groupby("grp").size().sort_values(ascending=False)
    kept, dropped = [], {}
    for g in n.index:
        dup = None
        if len(P):
            pp = P[(P.a == g) | (P.b == g)]
            for k in kept:
                q = pp[((pp.a == g) & (pp.b == k)) | ((pp.a == k) & (pp.b == g))]
                if len(q):
                    small = g if n[g] <= n[k] else k
                    f = q[q.a == small].frac_a.max() if (q.a == small).any() else 0.0
                    if f > frac:
                        dup = k
                        break
        if dup is None:
            kept.append(g)
        else:
            dropped[g] = dup
    fam = pd.DataFrame(dict(grp=n.index, n_det46=n.values))
    fam["kept"] = fam.grp.isin(kept)
    fam["dup_of"] = fam.grp.map(dropped)
    return fam, P


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grow", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--dates", default=",".join(C.ALL_DAYS))
    ap.add_argument("--mode", default="eliminate", choices=["eliminate", "merge"])
    ap.add_argument("--min-det", type=int, default=0, help="0 = the growth run's --min-det")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    ga = json.load(open(f"{a.grow}/args.json"))
    it = last_iter(a.grow)
    st = pd.read_csv(f"{a.grow}/state_it{it}.csv")
    md = a.min_det or ga["min_det"]
    keep = st[st.n_prev >= md].grp.tolist()
    T = load_tpl(f"{a.grow}/tpl_it{it}.npz")
    T = {g: T[g] for g in keep if g in T}
    A.log(f"{len(T)} templates with >= {md} detections in the 3-day window (it{it})")
    dc = A.DayCache(POOL)
    det = G.matched_filter(dc, T, a.dates.split(","), nmad=ga["nmad"], net=ga["net"],
                           min_sta_valid=ga["min_sta_valid"], decl_s=ga["tpl_len"],
                           tpl_len=ga["tpl_len"])
    det.to_parquet(f"{a.out}/det46.parquet", index=False)
    if a.mode == "eliminate":
        fam, P = eliminate(det)
        P.to_csv(f"{a.out}/dup_pairs.csv", index=False)
        out = det[det.grp.isin(fam[fam.kept].grp)]
    else:
        out, fam, stats = G.merge_duplicates(det)
        A.log(stats)
    fam.to_csv(f"{a.out}/families.csv", index=False)
    out.to_parquet(f"{a.out}/det46_families.parquet", index=False)
    A.log(f"families kept {int(fam.kept.sum()) if 'kept' in fam else fam.rep.nunique()}, "
          f"detections {len(out)}")
    open(f"{a.out}/SCAN46_DONE", "w").write("ok\n")


if __name__ == "__main__":
    main()
