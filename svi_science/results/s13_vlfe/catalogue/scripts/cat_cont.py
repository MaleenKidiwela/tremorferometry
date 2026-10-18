#!/usr/bin/env python
"""E3/E4 data (PLAN_CAT §8): continuous 0.02-0.05 Hz velocity in the hg_repro frame (sample 0 = 2011-06-01 UTC,
1 sps, 914 days) for stations not in the H&G set (CN.PGC, CN.SNB, CN.PFB, CN.YOUB), same recipe as hg_cont.py +
hg_r2.build_vel (v3 day files -> per-run detrend/taper/band-pass -> central-difference velocity), and the 90-s window
energy (den). Outputs data/cont/contv_STA.npy (3, N) float32 and data/cont/denv_STA.npy (3, M) float32.
The nine H&G stations are read from ../hg_repro/data/r2/contv_*.npy (read-only)."""
import os, sys
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K
import hg_common as H
import hg_cont

OUT = f'{K.DATA}/cont'
os.makedirs(OUT, exist_ok=True)


def build(sta):
    f = f'{OUT}/contv_{sta}.npy'
    if os.path.isfile(f) and os.path.isfile(f'{OUT}/denv_{sta}.npy'):
        return sta, 'exists'
    X = np.full((3, H.N), np.nan, np.float32)
    nd = 0
    for k in range(H.NDAY):
        d = H.T0 + pd.Timedelta(days=k)
        p = H.day_file(sta, d)
        if p is None:
            continue
        X[:, k * 86400:(k + 1) * 86400] = np.load(p)['d']
        nd += 1
    for c in range(3):
        x = hg_cont.filt(X[c]).astype(np.float64)
        v = np.full(x.shape, np.nan)
        v[1:-1] = 0.5 * (x[2:] - x[:-2])
        X[c] = v.astype(np.float32)
    np.save(f + '.part.npy', X); os.replace(f + '.part.npy', f)
    D = hg_cont.dens(X)
    np.save(f'{OUT}/denv_{sta}.part.npy', D); os.replace(f'{OUT}/denv_{sta}.part.npy', f'{OUT}/denv_{sta}.npy')
    return sta, f'days {nd}/{H.NDAY} live N/E/Z {(D > 0).mean(axis=1).round(3)}'


if __name__ == '__main__':
    from multiprocessing import Pool
    stas = sys.argv[1].split(',')
    with Pool(min(4, len(stas))) as p:
        for s, r in p.imap_unordered(build, stas):
            K.log(f'cont {s}: {r}', 'e3')
