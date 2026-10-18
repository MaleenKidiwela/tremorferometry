#!/usr/bin/env python
"""Run one validation proposal set through the gates and write results/gates/g_val_<set>.{csv,json}."""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import g_common as G                                             # noqa: E402
import g_gates as GG                                             # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("set")
    ap.add_argument("--gates", default="0,1,2")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--eval-days", default="")   # "YYYY-MM-DD,YYYY-MM-DD" for gate 3
    a = ap.parse_args()
    ctx = GG.GateContext(GG.GateConfig())
    props = GG.load_proposals(f"{G.GATES}/proposals_{a.set}.parquet")
    if a.limit:
        props = props[:a.limit]
    gates = tuple(int(x) for x in a.gates.split(","))
    ed = None
    if a.eval_days:
        lo, hi = a.eval_days.split(",")
        ed = np.arange(np.datetime64(lo), np.datetime64(hi) + np.timedelta64(1, "D"))
    rows, full = [], []
    t0 = time.time()
    for k, p in enumerate(props):
        try:
            r = GG.evaluate(ctx, p, gates=gates, evaluated_days=ed)
        except Exception as e:  # never lose a whole set to one bad proposal
            r = dict(pid=p.pid, lat=p.lat, lon=p.lon, dep=p.dep, n_members=len(p.times),
                     meta=dict(p.meta), verdict="ERROR", error=repr(e), runtime_s=0.0)
        rows.append(GG.flatten(r))
        full.append(r)
        if (k + 1) % 5 == 0 or k == 0:
            print(f"[{a.set} {k+1}/{len(props)}] {r['verdict']} "
                  f"({(time.time()-t0)/(k+1):.1f} s/proposal)", flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(f"{G.GATES}/g_val_{a.set}.csv", index=False)
    json.dump(full, open(f"{G.GATES}/g_val_{a.set}.json", "w"), indent=1, default=str)
    print(f"[{a.set}] {len(R)} proposals in {time.time()-t0:.0f} s")
    for c in ("verdict", "g0", "g1", "g2", "g3"):
        if c in R:
            print(c, dict(R[c].value_counts()))


if __name__ == "__main__":
    main()
