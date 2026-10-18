#!/usr/bin/env python
"""The CEILING: templates built directly from the corrected-time catalogued events of each kept
family, run through exactly the same matched filter at 8 x MAD over the same 24 ETS days.

This is what the chain could achieve if the seeding were perfect, and it is the reference the
recall of the grown families must be judged against (11-15 stations, not Bostock's full network).
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import k_common as K                                                        # noqa: E402
import k_grow as GR                                                         # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--times", default=f"{K.OUT}/s1all_times.csv")
    ap.add_argument("--stas", default=",".join(K.CLEAN11))
    ap.add_argument("--nev", type=int, default=200)
    ap.add_argument("--iters", type=int, default=1)
    ap.add_argument("--out", default=f"{K.OUT}/ceiling")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    stas = a.stas.split(",")
    C = pd.read_csv(a.times)
    C["t_corr"] = C.t_origin + C.dt
    groups = {}
    for b, g in C.groupby("bkey"):
        t = np.sort(g.t_corr.to_numpy())
        if len(t) > a.nev:                      # spread evenly over the window, not the busiest hour
            t = t[np.linspace(0, len(t) - 1, a.nev).astype(int)]
        if len(t) >= 10:
            groups[int(b)] = t
    K.log(f"{len(groups)} kept families, {sum(len(v) for v in groups.values())} stacked events")
    dc = K.DayCache(stas)
    dates = K.bed_days()
    for it in range(a.iters + 1):
        TPL, info = GR.build_templates(dc, stas, groups)
        import json as _j
        np.savez_compressed(f"{a.out}/tpl_it{it}.npz",
                            **{f"{g}_{s}": TPL[g][s][0].astype(np.float32)
                               for g in TPL for s in TPL[g]},
                            meta=np.array(_j.dumps({f"{g}_{s}": [TPL[g][s][1], TPL[g][s][2]]
                                                    for g in TPL for s in TPL[g]})))
        D = GR.matched_filter(dc, stas, TPL, dates)
        D.to_parquet(f"{a.out}/det_it{it}.parquet", index=False)
        n = D.groupby("grp").size()
        K.log(f"[ceiling iter {it}] {len(TPL)} templates -> {len(D)} det; median {n.median():.0f}")
        if it == a.iters:
            break
        groups = {int(g): np.sort(gg.nlargest(GR.MAX_STACK, "cc").t.to_numpy())
                  for g, gg in D.groupby("grp") if len(gg) >= 3}


if __name__ == "__main__":
    main()
