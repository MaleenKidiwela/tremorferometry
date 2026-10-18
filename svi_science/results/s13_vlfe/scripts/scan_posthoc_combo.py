#!/usr/bin/env python
"""POST-HOC exploration (after the NO-GO verdict; not a rule): headline-filtered positive clusters with z >= 7 and
z_loso >= T, (i) under the literal Merlin m45 mask and (ii) with the best node inside the LFE footprint polygon.
Rates per live year of the corresponding mask. Writes tables/scan_posthoc_combo.csv."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scan_common as S
from scan_validate import between_mask, live_years, R_
from matplotlib.path import Path
P = pd.read_csv(f'{S.C.TB}/scan_footprint_polygon.csv'); poly = Path(np.c_[P.lon, P.lat])
btw = between_mask(); nav = np.load(f'{R_}/real/navail.npy')
M = {t: S.mask_vector(t) for t in ('head', 'm45')}
real = pd.read_parquet(f'{R_}/real/loso.parquet')
A = {s: pd.read_parquet(f'{R_}/ctrlA_s{s}/loso.parquet') for s in (3, 5, 7, 11, 13)}
nA = {s: np.load(f'{R_}/ctrlA_s{s}/navail.npy') for s in A}
inj = pd.read_parquet(f'{R_}/inj2/loso.parquet'); J = pd.read_parquet(f'{S.SD}/inject_2011.parquet')
rows = []
for tier in ('head', 'm45'):
    for infp in (False, True):
        for T in (5, 7, 10):
            def sel(D):
                d = D[(D.z_loso >= T) & ~M[tier][D.tsec.values]]
                if infp:
                    d = d[poly.contains_points(np.c_[d.lon, d.lat])]
                return d
            r = sel(real); r0 = r[r.variant == 0]; rb = r[r.variant > 0]
            lyb = live_years(nav, M[tier], btw)
            a = np.mean([len(sel(A[s])) / live_years(nA[s], M[tier]) for s in A])
            ii = sel(inj); t = ii.tsec.values
            JJ = J[(J.gf == 'standard') & ~M[tier][J.tsec.astype(int).values]]
            if infp:
                JJ = JJ[JJ.loc_kind == 'family']
            rec = {mw: float(np.mean([np.any(np.abs(t - x) <= 60) for x in JJ[JJ.mw == mw].tsec])) for mw in (3.4, 3.6, 3.8, 4.0)}
            rows.append(dict(mask=tier, node_in_footprint=infp, z_loso_thr=T, live_yr_between=lyb,
                             real_between=int(btw[r0.tsec.values].sum()), real_between_per_yr=btw[r0.tsec.values].sum() / lyb,
                             ctrlB_between_per_yr=btw[rb.tsec.values].sum() / 50 / lyb, ctrlA_per_yr=a,
                             **{f'inj_rec_{k}': v for k, v in rec.items()}))
T = pd.DataFrame(rows); T.to_csv(f'{S.C.TB}/scan_posthoc_combo.csv', index=False); print(T.round(2).to_string())
