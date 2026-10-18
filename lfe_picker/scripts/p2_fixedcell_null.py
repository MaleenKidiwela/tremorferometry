#!/usr/bin/env python
"""Supplementary null: does the TRUE moveout matter at a FIXED cell?

The pre-registered moveout-shuffle null compares detection COUNTS, and a detection is the maximum of
the back-projection over 1,444 cells.  With that many cells a wrong travel-time table can still find
some cell whose alignment happens to be good, so the max-over-cells count is a weak test of the
geometry.  This test removes the maximisation: at each Bostock event, read the site-summed CF at the
moveout of that event's OWN family cell, under (a) the true table, (b) the shuffled table, (c) the
true table at random times on the same day.  If the physics carries information, (a) > (b) ~ (c).

  python scripts/p2_fixedcell_null.py --cf results/p2/cf_p2nnpu_2011 --year 2011
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P, p2_backproj as BP


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cf", required=True)
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--nctl", type=int, default=20)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    bp = BP.BackProjector(seed=a.seed)
    bs = BP.BackProjector(shuffle=True, seed=a.seed)
    st = P.DayStore()
    FG = pd.read_csv(f"{P.RES}/p2/family_on_grid.csv").set_index("bkey")
    B = P.bostock(); B = B[B.year == a.year]
    cnt = B.groupby("doy").size(); ED = set(cnt[cnt >= 100].index)
    B = B[B.doy.isin(ED)]
    rng = np.random.default_rng(3)
    rows = []
    for doy, g in B.groupby("doy"):
        T0 = (pd.Timestamp(f"{a.year}-01-01") + pd.Timedelta(days=int(doy) - 1)).timestamp()
        L, sok, used, val = BP.load_site_series(a.cf, P.S15, a.year, doy, st, T0)
        if len(sok) < 3:
            continue
        Lf = np.zeros((len(bp.sites), P.NDAY), np.float32)
        for j, s in enumerate(sok):
            Lf[bp.sites.index(s)] = L[j]
        n = P.NDAY
        for tab, name in ((bp.shift, "true"), (bs.shift, "shuffled")):
            c = np.array([int(FG.cell[b]) for b in g.bkey])
            j0 = np.rint((g.t_origin.values - T0) * P.FS).astype(np.int64)
            idx = (j0[:, None] + tab[c]).clip(0, n - 1)
            v = np.take_along_axis(Lf.T[:, None, :].repeat(1, 1), np.zeros((1, 1, 1), int), 0) \
                if False else Lf[np.arange(Lf.shape[0])[None, :], idx].sum(1)
            rows.append(pd.DataFrame(dict(doy=doy, kind=name, z=v)))
            if name == "true":                       # random-time control with the true table
                for _ in range(1):
                    jc = rng.integers(0, n - int(tab.max()) - 1, len(g) * a.nctl)
                    cc = c[rng.integers(0, len(c), len(jc))]
                    ii = (jc[:, None] + tab[cc]).clip(0, n - 1)
                    vv = Lf[np.arange(Lf.shape[0])[None, :], ii].sum(1)
                    rows.append(pd.DataFrame(dict(doy=doy, kind="random_time", z=vv)))
    R = pd.concat(rows, ignore_index=True)
    out = a.out or f"{P.RES}/p2/fixedcell_null_{a.year}.csv"
    s = R.groupby("kind").z.describe()[["count", "mean", "50%", "75%", "max"]]
    s.to_csv(out)
    print(s.to_string())
    t = R[R.kind == "true"].z.values
    sh = R[R.kind == "shuffled"].z.values
    rt = R[R.kind == "random_time"].z.values
    for nm, x in (("shuffled", sh), ("random_time", rt)):
        q = np.quantile(x, 0.99)
        print(f"  fraction of TRUE above the {nm} p99 ({q:.2f}): {(t > q).mean():.4f} "
              f"(chance 0.01)  -> enrichment {(t > q).mean()/0.01:.1f}x")
    print(f"  median lift true - shuffled: {np.median(t)-np.median(sh):+.3f} site-MAD units; "
          f"true - random_time: {np.median(t)-np.median(rt):+.3f}")


if __name__ == "__main__":
    main()
