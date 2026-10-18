#!/usr/bin/env python
"""Diagnostic (not a pre-registered rule): separation curve extended beyond z = 7 (headline mask + headline filter),
real positive / negative, control (a) mean of 5, control (b) between-ETS mean of 50 draws; plus ETS-only counts.
Writes tables/scan_separation_tail.csv."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scan_common as S
from scan_validate import headline_ok, between_mask, ets_vec, live_years, R_
use = ['n_sta', 'n_Z_agree', 'share_excess', 'vr_ratio', 'vr_gt_ss', 'vr_gt_steep', 'isolation', 'not_hf_loud']
mh = S.mask_vector('head'); btw = between_mask(); ets = ets_vec()
real = pd.read_parquet(f'{R_}/real/cand_cols.parquet'); nav = np.load(f'{R_}/real/navail.npy')
def prep(D):
    D = D[~mh[D.tsec.values]]
    return D[headline_ok(D, use)]
r = prep(real); r0 = r[r.variant == 0]; rb = r[r.variant > 0]; nd = int(real.variant.max())
A = {s: prep(pd.read_parquet(f'{R_}/ctrlA_s{s}/cand_cols.parquet')) for s in (3, 5, 7, 11, 13)}
nA = {s: np.load(f'{R_}/ctrlA_s{s}/navail.npy') for s in A}
ly = live_years(nav, mh); lyb = live_years(nav, mh, btw); lye = live_years(nav, mh, ets)
rows = []
for thr in (7, 8, 10, 12, 15, 20, 25, 30):
    p = r0[(r0.sign > 0) & (r0.z >= thr)]; n = r0[(r0.sign < 0) & (-r0.z >= thr)]
    a = np.mean([((A[s].sign > 0) & (A[s].z >= thr)).sum() / live_years(nA[s], mh) for s in A])
    bb = rb[(rb.sign > 0) & (rb.z >= thr) & btw[rb.tsec.values]]
    be = rb[(rb.sign > 0) & (rb.z >= thr) & ets[rb.tsec.values]]
    rows.append(dict(thr=thr, real_pos_per_yr=len(p) / ly, real_neg_per_yr=len(n) / ly, ctrlA_per_yr=a,
                     ctrlB_between_per_yr=len(bb) / nd / lyb, real_pos_between_per_yr=btw[p.tsec.values].sum() / lyb,
                     real_pos_ets=int(ets[p.tsec.values].sum()), real_neg_ets=int(ets[n.tsec.values].sum()),
                     ctrlB_ets_per_draw=len(be) / nd, real_pos_between=int(btw[p.tsec.values].sum()),
                     ets_live_days=lye * 365.25))
T = pd.DataFrame(rows); T.to_csv('tables/scan_separation_tail.csv', index=False); print(T.round(2).to_string())
