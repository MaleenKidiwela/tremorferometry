#!/usr/bin/env python
"""Diagnostics after Tests 1-4 (NOTES ~10:25 am PT, post-hoc, disclosed): (i) spread of the LFEs themselves inside the
+-10 min window (how sharp is the 'truth'); (ii) offsets vs live western stations (PFB, YOUB, NLLB, LZB) per event.
-> tables/t5_diag.csv, tables/t5_diag_summary.csv"""
import os, sys
import numpy as np, pandas as pd
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import t1_lfe_centroids as T1
E = pd.read_csv(f'{OUT}/tables/t2_zsurface_all.csv')
E = E[(E.period == 'ETS') & np.isfinite(E.c10_lat)]
L = T1.lfe_table(); L = L[np.isfinite(L.lat)]
rows = []
WEST = ('CN.PFB', 'CN.YOUB', 'CN.NLLB')
Wd = {y: np.load(f'{OUT}/data/win_{y}.npz') for y in sorted(E.year.unique())}
for r in E.itertuples():
    a, b = np.searchsorted(L.t.values, r.epoch - 600), np.searchsorted(L.t.values, r.epoch + 600, 'right')
    x = L.iloc[a:b]
    lx, ly = T1.xy(x.lat.values, x.lon.values); cx, cy = T1.xy(r.c10_lat, r.c10_lon)
    d = np.hypot(lx - cx, ly - cy)
    W = Wd[r.year]; j = list(W['eid']).index(r.eid); st = list(W['stas'])
    nwest = int(sum(W['won'][j][st.index(s)].sum() for s in WEST if s in st))
    nall = int(W['won'][j].sum())
    rows.append(dict(eid=r.eid, year=r.year, n_lfe=len(x), lfe_rms_km=float(np.sqrt((d ** 2).mean())), lfe_max_km=float(d.max()),
                     n_west_ch=nwest, n_ch=nall, c10_dist=r.c10_dist, c10_dE=r.c10_dE, c10_dN=r.c10_dN, z=r.z_stored))
D = pd.DataFrame(rows); D.to_csv(f'{OUT}/tables/t5_diag.csv', index=False)
S = []
for nm, g in (('all', D), ('n_lfe>=5', D[D.n_lfe >= 5]), ('west on (>=3 ch)', D[D.n_west_ch >= 3]), ('west off (<3 ch)', D[D.n_west_ch < 3])):
    S.append(dict(subset=nm, n=len(g), lfe_rms_median=g.lfe_rms_km.median(), offset_median=g.c10_dist.median(),
                  mean_dE=g.c10_dE.mean(), mean_dN=g.c10_dN.mean(), years=','.join(map(str, sorted(g.year.unique())))))
S = pd.DataFrame(S); S.to_csv(f'{OUT}/tables/t5_diag_summary.csv', index=False)
pd.set_option('display.width', 250); print(S.round(2).to_string())
print(D.groupby('year')[['n_west_ch', 'n_ch', 'lfe_rms_km', 'c10_dist']].median())
