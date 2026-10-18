#!/usr/bin/env python
"""Internal coherence of the western stacks vs random times (coordinator point 5, 2026-10-16): each member's channel-mean CC
with the stack of the other members (best of +-5 s, catalogue alignment range) compared with the same statistic for 300
random 2018 ETS-window times against the same leave-one-out stacks; RC4 / RC0 / RC1 relay members as reference.
-> tables/nv_stackw_coh.csv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import nv_stackw as SW
P = SW.P
B = np.load(f'{SW.OUT}/bank.npz')
REL = pd.read_csv(f'{SW.CAT}/tables/prod2018_relay_events.csv')
groups = {'W6': list(B['starts__W6']), 'W17': list(B['starts__W17'])}
for c in (0, 1, 4):
    g = REL[REL.cluster == c]; groups[f'RC{c}'] = list((g.start0 + g.shift_s).astype(int))
rng = np.random.default_rng(7)
a, b = P.tidx('2018-06-11'), P.tidx('2018-07-20')
rand = [int(x) for x in rng.integers(a, b, 300)]
out = []
for gname, st in groups.items():
    for k, s0 in enumerate(st):
        ref = SW.stack_frame([s for j, s in enumerate(st) if j != k])
        def best(s):
            v = [SW.sumcc(P.cut_frame(s + d), ref) for d in range(-5, 6)]
            return max(x / n for x, n in v if n)
        m = best(s0)
        r = np.array([best(t) for t in rand])
        out.append(dict(group=gname, member=k, cc_member=m, rand_median=float(np.median(r)), rand_p95=float(np.percentile(r, 95)),
                        pct=float((r < m).mean())))
    print(gname, pd.DataFrame([o for o in out if o['group'] == gname]).round(3).to_string())
pd.DataFrame(out).to_csv(f'{SW.NV}/tables/nv_stackw_coh.csv', index=False)
