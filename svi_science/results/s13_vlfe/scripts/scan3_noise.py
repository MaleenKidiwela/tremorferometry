#!/usr/bin/env python
"""VLFE scan pilot 3 (PLAN_SCAN.md S8.2): noise table for the scan year (SCAN_YEAR, e.g. 2012), same rule as scan_calib.noise
(PLAN_SCAN §3): per station-component-month median of the daily median 400-s band RMS of the v3 day files of that year;
a station-month is excluded if its Z median > 3 x the station's all-days median (tables/scan_noise_ref.csv, computed in 2011
from every v2+v3 day on disk) or > 5 x the network median of Z medians that month.
Writes tables/scan3_noise_month.csv (S.tab) and data/scan_<year>/noise_days.csv. usage: SCAN_YEAR=2012 scan3_noise.py"""
import os, sys, glob
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import scan_common as S
from scan_calib import day_noise, STAS

if __name__ == '__main__':
    assert S.YEAR != 2011, 'use scan_calib.py noise for 2011'
    from multiprocessing import Pool
    files = []
    for s in STAS:
        files += sorted(glob.glob(f'{S.V3}/{s}/{S.YEAR}/*.npz'))
    rows = []
    with Pool(int(os.environ.get('NWORK', '5'))) as p:
        for f, o in p.imap_unordered(day_noise, [(f,) for f in files], chunksize=8):
            parts = f.split('/')
            rows.append(dict(sta=parts[-3], year=int(parts[-2]), doy=int(parts[-1][:3]), N=o[0], E=o[1], Z=o[2]))
    D = pd.DataFrame(rows)
    D['date'] = pd.to_datetime(D.year.astype(str) + '-' + D.doy.astype(str), format='%Y-%j')
    D.to_csv(f'{S.SD}/noise_days.csv', index=False)
    ref = pd.read_csv(f'{C.TB}/scan_noise_ref.csv').set_index('sta').Z_ref_all_days
    D['month'] = D.date.dt.month
    M = D.groupby(['sta', 'month'])[['N', 'E', 'Z']].median().reset_index()
    netmed = M[M.sta.isin(S.PILOT + S.NOISE_ONLY)].groupby('month').Z.median().rename('Z_net_median')
    M = M.join(ref, on='sta').join(netmed, on='month')
    M['n_days'] = D.groupby(['sta', 'month']).size().values
    M['exclude'] = (M.Z > 3 * M.Z_ref_all_days) | (M.Z > 5 * M.Z_net_median)
    for c in ('N', 'E', 'Z', 'Z_ref_all_days', 'Z_net_median'):
        M[c] = M[c] * 1e9
    M.to_csv(S.tab('scan_noise_month.csv'), index=False)
    print(M.round(1).to_string())
