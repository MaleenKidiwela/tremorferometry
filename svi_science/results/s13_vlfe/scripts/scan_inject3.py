#!/usr/bin/env python
"""VLFE scan pilot 2 (PLAN_SCAN.md S7.8): extra injection set inj3. 40 at Mw 3.8 and 40 at Mw 3.6, standard GF (free surface
2, Vp 6.9, Vs 3.95 km/s), footprint-family locations (family depth), half 1-s / half 20-s source (tau = D/2.33, S6).
Time rules as scan_inject.py (outside the headline mask v1 'head', >= 3 pilot stations live, >= 2 h apart, >= 2 h from
H&G template times) and >= 2 h from every inj2 time. Seed 20261016. Writes data/scan/inject3_2011.parquet."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import scan_common as S
from scan_detect import gain_table, day_cfg

rng = np.random.default_rng(20261016)
TP = np.load(f'{S.SD}/templates_grid.npz')
stas = list(TP['stations'])
Gt, excl = gain_table(); g = day_cfg(Gt, excl, stas, None)
BP = {s: np.load(f'{S.SD}/bp_{S.YEAR}g/{s}.npy', mmap_mode='r') for s in stas}     # SCAN_YEAR (S8)
mask = S.mask_vector('head')
F = pd.read_csv(f'{C.TB}/scan_footprint_families.csv'); F = F[F.in_footprint_set]
hg = [pd.Timestamp(t).value / 1e9 for t in pd.read_csv(f'{C.TB}/hg2018_table1.csv', comment='#').time]
J2 = pd.read_parquet(f'{S.SD}/inject_{S.YEAR}.parquet')
y0 = S.Y0.value / 1e9
specs = []
for mw in (3.6, 3.8):
    for k in range(40):
        specs.append(dict(mw=mw, tau=1.0 if k < 20 else 20.0 / 2.33, dur_s=1.0 if k < 20 else 20.0, fs=2.0, vs=3.95, vp=6.9,
                          gf='standard'))
rng.shuffle(specs)
times = []
avoid = list(J2.tsec.values)
while len(times) < len(specs):
    t = rng.uniform(1800, S.NDAY * 86400 - 1800)
    ti = int(t)
    if mask[ti - 600:ti + 600].any():
        continue
    if any(abs(t - x) < 7200 for x in times + avoid) or any(abs(y0 + t - h) < 7200 for h in hg):
        continue
    d = ti // 86400
    live = 0
    for i, s in enumerate(stas):
        if not (g[i, :, d] > 0).any():
            continue
        x = np.asarray(BP[s][:, ti - 660:ti + 721])
        if any((g[i, c, d] > 0) and np.isfinite(x[c]).all() for c in range(3)):
            live += 1
    if live >= 3:
        times.append(t)
rows = []
for k, (sp, t) in enumerate(zip(specs, times)):
    f = F.iloc[rng.integers(len(F))]
    rows.append(dict(inj=1000 + k, t=y0 + t, tsec=t, lat=f.lat, lon=f.lon, dep=f.depth_km, loc_kind='family', fam=int(f.fam),
                     m0=10 ** (1.5 * sp['mw'] + 9.1), **sp))
J = pd.DataFrame(rows)
J['time_utc'] = pd.to_datetime(J.t, unit='s')
J['season'] = S.season(J.time_utc)
J.to_parquet(f'{S.SD}/inject3_{S.YEAR}.parquet')
print(J.groupby(['mw', 'dur_s']).size(), J.season.value_counts())
