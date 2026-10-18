#!/usr/bin/env python
"""Selection-matched reference for the western-stack coherence: scan candidates chosen by the same synthetic grid at the same
western nodes (48.45-48.95 N, west of 124.2 W, 2018, gate + single-station test passed), between ETS (real) and in the ETS
polarity controls' timeline is not usable (different template), so: real between-ETS rows; their channel-mean CC with the
W17 and W6 stacks (best of +-5 s) vs the members' leave-one-out CC. -> tables/nv_stackw_ref.csv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import nv_stackw as SW
P = SW.P
B = np.load(f'{SW.OUT}/bank.npz')
D = pd.read_parquet(f'{SW.NV}/tables/nv_rows_scored.parquet')
x = D[(D.run == 'real') & (D.year == 2018) & (D.period == 'between') & D.gate & D.pass_f_light.fillna(False)
      & (D.lat >= 48.45) & (D.lat <= 48.95) & (D.lon <= -124.2)]
print('reference rows', len(x), 'of which B0', int(x.B0.sum()), 'B1', int(x.B1.sum()))
out = []
for nm in ('W6', 'W17'):
    st = {ch: (B[nm][k] if np.isfinite(B[nm][k]).all() else None) for k, ch in enumerate(P.CH)}
    for r in x.itertuples():
        s = int(round(r.epoch - P.EP0 - P.PRE))
        v = [SW.sumcc(P.cut_frame(s + d), st) for d in range(-5, 6)]
        v = [a / n for a, n in v if n]
        out.append(dict(stack=nm, epoch=r.epoch, z=r.z_best, B0=r.B0, B1=r.B1, cc=max(v) if v else np.nan))
R = pd.DataFrame(out); R.to_csv(f'{SW.NV}/tables/nv_stackw_ref.csv', index=False)
print(R.groupby('stack').cc.describe(percentiles=[0.5, 0.9, 0.95]).round(3))
print(R[R.B1].groupby('stack').cc.describe().round(3))
