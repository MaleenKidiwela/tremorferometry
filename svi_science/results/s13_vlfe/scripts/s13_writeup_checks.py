#!/usr/bin/env python
"""S13 write-up checks (before gap fill), 2026-10-14 evening PT. Reads the existing scale_* tables and stacks only.
1. Flat-channel units: control SD < 1e10 N m (a channel with sigma_ctrl = 0, e.g. UW.LRIV N/Z 2019-04/05, UW.OPC 2009-04)
   -> left out of the cleaned family / zone tables (s13_scale_summary.py does not guard against this).
2. Re-fit of every pooled / between-ETS stack from data/stacks/*.pkl with the stored 1/sigma^2 weights, (a) all channels
   (must reproduce thrust_m0) and (b) without UW.SP2 E (near-dead 2012-2019: sigma 0.04-0.35 nm vs 7-20 nm on SP2 N)
   and without sigma = 0 channels. M0 only: control draws are not stored, so no new z.
Outputs: tables/writeup_bgf_refit.csv, tables/writeup_bgf_{svi,wa}_families_clean.csv, tables/writeup_bgf_{svi,wa}_zones_clean.csv"""
import os, sys, glob, pickle
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
from s13_stack import fit

out = []
for f in sorted(glob.glob(f'{C.DATA}/stacks/scale_*.pkl')):
    tag = os.path.basename(f)[:-4]
    st = pd.read_csv(f'{C.TB}/{tag}_stacomp.csv'); un = pd.read_csv(f'{C.TB}/{tag}_units.csv').set_index('unit')
    D = pickle.load(open(f, 'rb'))['units']
    for u in [x for x in D if x == 'pooled' or 'BE' in x]:
        d = D[u]; s = st[st.unit == u]
        sig = {k: np.array([s[(s.sta == k) & (s.comp == c)].sigma_ctrl_nm.iloc[0] for c in 'NEZ']) * 1e-9 for k in d['stk']}
        def run(excl, floor):
            w = {k: np.array([0.0 if (sig[k][i] <= 0 and not floor) or (k, 'NEZ'[i]) in excl
                              else 1 / max(sig[k][i], 1e-15) ** 2 for i in range(3)]) for k in sig}
            w = {k: v for k, v in w.items() if v.sum() > 0}
            return fit({k: d['stk'][k] for k in w}, {k: d['tpl'][k] for k in w}, w)
        a = run(set(), True); b = run({('UW.SP2', 'E')}, False)
        out.append(dict(tag=tag, unit=u, n=un.loc[u, 'n'], thrust_z=un.loc[u, 'thrust_z'], m0_stored=un.loc[u, 'thrust_m0'],
                        m0_refit_all=a[0], m0_no_sp2e_no_flat=b[0], vr_no_sp2e_no_flat=b[2]))
pd.DataFrame(out).to_csv(f'{C.TB}/writeup_bgf_refit.csv', index=False)

for reg in ('svi', 'wa'):
    e = pd.read_csv(f'{C.TB}/scale_{reg}_episodes.csv')
    U = pd.concat([pd.read_csv(f'{C.TB}/scale_{reg}_{ep}_units.csv').assign(episode=ep) for ep in e[e.built].episode])
    U = U[U.thrust_ctrl_sd >= 1e10]
    for kind in ('fam', 'zone'):
        X = U[U.unit.str.startswith(kind)].copy(); X['w'] = 1 / X.thrust_ctrl_sd ** 2
        T = X.groupby('unit').apply(lambda g: pd.Series(dict(
            n_ep=len(g), n_lfe=g.n.sum(), m0_wmean=(g.thrust_m0 * g.w).sum() / g.w.sum(), m0_sd=1 / np.sqrt(g.w.sum()),
            chi2=((g.thrust_m0 - (g.thrust_m0 * g.w).sum() / g.w.sum()) ** 2 * g.w).sum(), dof=len(g) - 1,
            n_z3=(g.thrust_z >= 3).sum(), n_zm3=(g.thrust_z <= -3).sum(),
            lat=g.src_lat.mean(), lon=g.src_lon.mean(), dep=g.src_dep.mean())))
        T['z_wmean'] = T.m0_wmean / T.m0_sd
        T.to_csv(f'{C.TB}/writeup_bgf_{reg}_{"families" if kind == "fam" else "zones"}_clean.csv')
print('ok')
