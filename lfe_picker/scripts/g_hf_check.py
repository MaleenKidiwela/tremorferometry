#!/usr/bin/env python
"""Why the Gate-0 HF class statistic is the MEDIAN OVER MEMBERS and not a ratio measured on the
coherent stack: at 10-19 Hz the per-member alignment error (>= 0.1-0.5 s) is many half-periods, so the
coherent stack destroys the HF band and its ratio is biased LOW by a large factor, which would make the
veto unable to fire.  Measured here on the positive proposals.  -> results/gates/g_hf_stack_vs_member.csv"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import g_common as G                                             # noqa: E402
import g_gates as GG                                             # noqa: E402


def main():
    ctx = GG.GateContext(GG.GateConfig())
    rows = []
    for setname in ("pos_bost_N200", "neg_hf_N50"):
        props = GG.load_proposals(f"{G.GATES}/proposals_{setname}.parquet")
        for p in props[:12]:
            sta = "PGC"
            _d, ttS, _ = ctx.pred_S(sta, p.lat, p.lon, p.dep)
            t = p.times[:100] + ttS
            ts = pd.to_datetime(t, unit="s")
            by = {}
            for k, tt in zip(zip(ts.year.to_numpy().tolist(), ts.dayofyear.to_numpy().tolist()), t):
                by.setdefault(k, []).append(float(tt))
            r = G.hf_ratio_members(ctx.nstore, ctx.netsta(sta), by)
            # coherent native-rate stack at the predicted S
            segs = []
            for (year, doy), times in by.items():
                v = ctx.nstore.get(ctx.netsta(sta), year, doy)
                if v is None:
                    continue
                X, fs, t0 = v
                n = int(round(12.0 * fs))
                for x in times:
                    i0 = int(round((x - 1.0 - t0) * fs))
                    if i0 < 0 or i0 + n > X.shape[1]:
                        continue
                    sg = X[:, i0:i0 + n].astype(float)
                    rms = np.sqrt((sg ** 2).mean())
                    if rms > 0:
                        segs.append(sg / rms)
            if len(segs) < 10 or r.size < 10:
                continue
            st = np.stack(segs).mean(0)
            eh = float((G._bp(st, fs, 10.0, 19.0) ** 2).sum())
            el = float((G._bp(st, fs, 2.0, 8.0) ** 2).sum())
            rows.append(dict(set=setname, pid=p.pid, n=len(segs),
                             member_median=float(np.median(r)), stack_ratio=eh / el,
                             bias=float(np.median(r)) / (eh / el)))
    R = pd.DataFrame(rows)
    R.to_csv(f"{G.GATES}/g_hf_stack_vs_member.csv", index=False)
    print(R.to_string(index=False))
    print("\nmedian bias member/stack by set:")
    print(R.groupby("set").bias.median().to_string())


if __name__ == "__main__":
    main()
