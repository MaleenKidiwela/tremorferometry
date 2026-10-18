#!/usr/bin/env python
"""Amendment W1 diagnostics: per station-component-month median 10-min RMS and glitch rate of the MF velocity data.
usage: cat_diag.py YEAR [YEAR..] -> tables/diag_noise_<YEAR>.csv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

for y in map(int, sys.argv[1:]):
    T = pd.read_csv(f'{K.TB}/prod_stations.csv'); st = T[T.year == y].sta.tolist()
    rows = []
    for s in st:
        cd = f'{K.DATA}/prodf{y}/cont' if os.path.isdir(f'{K.DATA}/prodf{y}/cont') else f'{K.DATA}/prod{y}/cont'
        X = np.load(f'{cd}/contv_{s}.npy', mmap_mode='r')
        n = X.shape[1] // 600 * 600
        for c, cn in enumerate('NEZ'):
            x = np.asarray(X[c, :n], np.float64).reshape(-1, 600)
            ok = np.isfinite(x).all(1)
            r = np.where(ok, np.sqrt(np.nanmean(x * x, 1)), np.nan)
            ymed = np.nanmedian(r)
            day = np.arange(len(r)) // 144
            mon = (pd.Timestamp(f'{y}-01-01') + pd.to_timedelta(day, unit='D')).month
            for m in range(1, 13):
                rm = r[(mon == m) & np.isfinite(r)]
                rows.append(dict(year=y, sta=s, comp=cn, month=m, n10=len(rm), med_rms=float(np.median(rm)) if len(rm) else np.nan,
                                 glitch=float(np.mean(rm > 10 * ymed)) if len(rm) else np.nan, year_med=ymed))
    D = pd.DataFrame(rows); D.to_csv(f'{K.TB}/diag_noise_{y}.csv', index=False)
    z = D[D.comp == 'Z'].pivot_table(index='sta', columns='month', values='med_rms')
    g = D[D.comp == 'Z'].pivot_table(index='sta', columns='month', values='glitch')
    pd.set_option('display.width', 250)
    print(y, 'Z median 10-min RMS (nm/s)'); print((z * 1e9).round(1).to_string()); print('glitch rate'); print((g * 100).round(1).to_string())
