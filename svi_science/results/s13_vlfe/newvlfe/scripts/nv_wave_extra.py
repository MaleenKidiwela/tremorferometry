#!/usr/bin/env python
"""Raw waveform figures for selected new-ETS real candidates and polarity-control rows of similar z (coordinator request
2026-10-16 ~10:50 am PT). Uses nv_survivors.wave_fig (stations > 10x the median peak clipped). -> figures/nv_wave_*.png,
figures/nv_ctrlwave_*.png"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import nv_survivors as V
import scan_refs as SR
cc = SR.all_coords()
R = pd.read_csv(f'{V.NV}/tables/nv_new_candidates.csv'); C = pd.read_csv(f'{V.NV}/tables/nv_new_controls.csv')
R = R[R.period == 'ETS'].sort_values('z_best', ascending=False)
sel = pd.concat([R.head(6), R[R.patch_p < 0.05]]).drop_duplicates('epoch')
for e in sel.itertuples():
    V.wave_fig(e, f"{V.NV}/figures/nv_wave_{e.time_pt[:19].replace(' ', '_').replace(':', '')}.png", cc)
    print('real', e.time_pt)
Cp = C[(C.period == 'ETS') & (C.run == 'polB')]
for z in R.head(6).z_best:
    e = Cp.iloc[int(np.argmin(np.abs(Cp.z_best.values - z)))]
    V.wave_fig(e, f"{V.NV}/figures/nv_ctrlwave_{e.time_pt[:19].replace(' ', '_').replace(':', '')}_v{int(e.variant)}.png", cc)
    print('ctrl', e.time_pt, e.z_best)
