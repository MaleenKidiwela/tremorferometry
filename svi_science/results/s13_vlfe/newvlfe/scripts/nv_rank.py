#!/usr/bin/env python
"""Rank curves (coordinator request 2026-10-16 ~10:50 am PT): cumulative B1 count above z, real vs expected from each
control type (count x real live / control live), new and known, ETS, per year; LRIV-broken flag (|Z| > 1 um in the fit
window). -> figures/nv_rank_ets.png, tables/nv_rank_ets.csv, tables/nv_ets_b1_rows.csv"""
import json, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
S13 = '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe'; NV = f'{S13}/newvlfe'
D = pd.read_parquet(f'{NV}/tables/nv_rows_scored.parquet')
X = D[D.B1 & (D.period == 'ETS')].copy()
A = {y: np.load(f'{S13}/data/scan_{y}/bp_{y}v2/UW.LRIV.npy', mmap_mode='r') for y in (2013, 2018)}
X['lriv_pk_nm'] = [np.nanmax(np.abs(np.asarray(A[int(r.year)][2, int(r.tsec) - 60:int(r.tsec) + 121]) * 1e9)) for r in X.itertuples()]
X['lriv_broken'] = X.lriv_pk_nm > 1000
X['ctype'] = np.where(X.run == 'real', 'real', np.where(X.run == 'polB', 'polB', 'shift'))
live = {}
for f in glob.glob(f'{NV}/data/live_*.json'):
    L = json.load(open(f)); ct = 'real' if L['run'] == 'real' else ('polB' if L['run'] == 'polB' else 'shift')
    live[(L['year'], ct)] = live.get((L['year'], ct), 0) + L['live_yr_per_catalogue']['ETS'] * L['n_catalogues']
X[['time_pt', 'run', 'variant', 'year', 'loc', 'lat', 'lon', 'z_best', 'z_loso', 'ratio_far', 'vr_thrust', 'lriv_pk_nm', 'lriv_broken',
   'dist_tpl_km', 'nearest_tpl']].to_csv(f'{NV}/tables/nv_ets_b1_rows.csv', index=False)
zs = np.r_[5:30.5:0.5]; rows = []
fig, ax = plt.subplots(2, 3, figsize=(14, 8))
for i, loc in enumerate(('new', 'known')):
    for j, yr in enumerate((2013, 2018, 'both')):
        a = ax[i, j]; ys = (2013, 2018) if yr == 'both' else (yr,)
        lr = sum(live[(y, 'real')] for y in ys)
        x = X[(X['loc'] == loc) & X.year.isin(ys)]
        cr = [int(((x.ctype == 'real') & (x.z_best >= z)).sum()) for z in zs]
        a.plot(zs, cr, 'o-', color='tab:red', ms=3, label='real')
        for ct, col in (('polB', 'tab:blue'), ('shift', 'tab:gray')):
            lc = sum(live[(y, ct)] for y in ys)
            ex = [((x.ctype == ct) & (x.z_best >= z)).sum() / lc * lr for z in zs]
            a.plot(zs, ex, 's--', color=col, ms=3, mfc='none', label=f'expected from {ct}')
            for z, c_, e_ in zip(zs, cr, ex):
                rows.append(dict(loc=loc, year=yr, z=z, ctype=ct, real=c_, expected=e_))
        crl = [int(((x.ctype == 'real') & (x.z_best >= z) & ~x.lriv_broken).sum()) for z in zs]
        a.plot(zs, crl, 'o:', color='darkred', ms=2, label='real, LRIV not broken')
        a.set_yscale('symlog', linthresh=1); a.set_title(f'B1, ETS, {loc}, {yr}'); a.legend(fontsize=7)
        a.set_xlabel('z threshold'); a.set_ylabel('events ≥ z')
fig.suptitle('Rank curves: real vs expected from scrambled controls (ETS, tier B1)')
fig.tight_layout(); fig.savefig(f'{NV}/figures/nv_rank_ets.png', dpi=120)
pd.DataFrame(rows).to_csv(f'{NV}/tables/nv_rank_ets.csv', index=False)
print({k: round(v, 3) for k, v in live.items()})
