#!/usr/bin/env python
"""VLFE scan 6d (PLAN_SCAN.md §6d): injection list. Writes data/scan/inject_2011.parquet.
Mw 2.8..4.0 (step 0.2) x 40 (20 with 1-s, 20 with 20-s source; Amendment S6: tau = D/2.33), standard GF (free surface 2, Vp 6.9, Vs 3.95 km/s);
plus 40 at Mw 3.6 with a perturbed GF (free surface 1.0, Vs 3.4, Vp 6.0 km/s, tau=30 s). Times uniform over 2011, outside
the headline mask, >= 3 pilot stations live (gain kept, band-passed data finite over the fit window +-600 s), >= 2 h apart,
>= 2 h from any H&G template time. Half at random footprint-family locations (family depth), half at grid nodes 30-50 km
from the nearest family (Slab2 depth). Seed 20261015."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import scan_common as S
from scan_detect import gain_table, day_cfg

rng = np.random.default_rng(20261015)
TP = np.load(f'{S.SD}/templates_grid.npz')
stas = list(TP['stations'])
Gt, excl = gain_table(); g = day_cfg(Gt, excl, stas, None)
BP = {s: np.load(f'{S.SD}/bp_{S.YEAR}g/{s}.npy', mmap_mode='r') for s in stas}     # SCAN_YEAR (S8)
mask = S.mask_vector('head')
F = pd.read_csv(f'{C.TB}/scan_footprint_families.csv')
F = F[F.in_footprint_set]
N = pd.read_csv(f'{C.TB}/scan_grid_nodes.csv').iloc[:-1]
dmin = np.array([S.dist_km(r.lat, r.lon, F.lat.values, F.lon.values).min() for r in N.itertuples()])
far = N[(dmin >= 30) & (dmin <= 50)].reset_index(drop=True)
hg = [pd.Timestamp(t).value / 1e9 for t in pd.read_csv(f'{C.TB}/hg2018_table1.csv', comment='#').time]
y0 = S.Y0.value / 1e9
specs = []
for mw in (2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0):
    for k in range(40):
        specs.append(dict(mw=mw, tau=1.0 if k < 20 else 20.0 / 2.33, dur_s=1.0 if k < 20 else 20.0, fs=2.0, vs=3.95, vp=6.9, gf='standard'))
for k in range(40):
    specs.append(dict(mw=3.6, tau=30.0 / 2.33, dur_s=30.0, fs=1.0, vs=3.4, vp=6.0, gf='perturbed'))
rng.shuffle(specs)
times = []
while len(times) < len(specs):
    t = rng.uniform(1800, S.NDAY * 86400 - 1800)
    ti = int(t)
    if mask[ti - 600:ti + 600].any():
        continue
    if any(abs(t - x) < 7200 for x in times) or any(abs(y0 + t - h) < 7200 for h in hg):
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
    if k % 2 == 0:
        f = F.iloc[rng.integers(len(F))]; loc = (f.lat, f.lon, f.depth_km, 'family', int(f.fam))
    else:
        n = far.iloc[rng.integers(len(far))]; loc = (n.lat, n.lon, n.dep, 'off_family', -1)
    rows.append(dict(inj=k, t=y0 + t, tsec=t, lat=loc[0], lon=loc[1], dep=loc[2], loc_kind=loc[3], fam=loc[4],
                     m0=10 ** (1.5 * sp['mw'] + 9.1), **sp))
J = pd.DataFrame(rows)
J['time_utc'] = pd.to_datetime(J.t, unit='s')
J['season'] = S.season(J.time_utc)
J.to_parquet(f'{S.SD}/inject_{S.YEAR}.parquet')
print(J.groupby(['gf', 'mw', 'tau']).size(), J.season.value_counts(), J.loc_kind.value_counts())
