#!/usr/bin/env python
"""Raw waveform figures (0.02-0.05 Hz vertical displacement: scan near network, H&G network, far references;
t0 - 300 .. t0 + 400 s) for the E2 between-ETS battery' events, using hg_r3_summary.wave_fig unchanged.
-> figures/e2_wave_<time>.png"""
import os, sys
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import cat_common as K
import hg_r3_summary as RS

E = pd.read_csv(f'{K.TB}/e2_real_battery_events.csv')
E = E[E.period == 'between'].drop_duplicates('i')
D = pd.read_parquet(f'{K.DATA}/e2/e2_rows_scored.parquet')
for k, e in enumerate(E.itertuples()):
    d = D[(D.source == 'real') & (D.i == e.i)].iloc[0].copy()
    d['cand_id'] = k
    d['time_utc'] = str(pd.Timestamp(d.epoch, unit='s'))[:19]
    d['time_pt'] = e.time_pt
    d['class'] = f"battery{'' if d.battery else ' (fails e)'} e' pass, agree {e.agree_z:.1f} MAD, {e.lfe_status}"
    tag = e.time_pt[:19].replace(' ', '_').replace(':', '')
    fn = f'{K.FIG}/e2_wave_{tag}.png'
    RS.wave_fig(d, fn)
    print(fn)
