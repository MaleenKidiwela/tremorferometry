#!/usr/bin/env python
"""S13 check column (2026-10-14 PT, not a verdict): Amendment 4 v2 drops healthy verticals whose horizontals are tilt-noisy
(clause sigma < 0.05 x the station's other two). Re-fit M0 of every SVI pooled / between unit from the stored v2 stacks
with those Z channels added back (1/sigma^2), all other v2 drops kept. M0 only (control draws are not stored, so no z).
Output tables/guard2_zvariant.csv."""
import os, sys, pickle
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
from s13_stack import fit
import s13_scale_summary as S

out = []
E = pd.read_csv(f'{C.TB}/svi_episodes.csv')
for b in [f'{e}' for e in E[E.built].episode] + ['between']:
    f, pre = S.pick('svi', b)
    if not f:
        continue
    U = pd.read_csv(f).set_index('unit')
    if 'guard_rule' not in U.columns:
        continue                                   # first-pass table with no v2 drop: variant = as run
    st = pd.read_csv(f.replace('_units.csv', '_stacomp.csv'))
    D = pickle.load(open(f'{C.DATA}/stacks/{os.path.basename(f)[:-10]}.pkl', 'rb'))['units']
    for u in [x for x in D if x == 'pooled' or 'BE' in x]:
        d = D[u]; s = st[st.unit == u]; w = {k: v.copy() for k, v in d['wts'].items()}
        back = []
        for r in s.itertuples():
            if r.comp == 'Z' and r.drop_reason == '<0.05xsta':
                w[r.sta][2] = 1 / (r.sigma_ctrl_nm * 1e-9) ** 2; back.append(r.sta)
        m = fit(d['stk'], d['tpl'], w)[0] if back else U.loc[u, 'thrust_m0']
        out.append(dict(table=pre, rowfile=b, unit=u, n=U.loc[u, 'n'], m0_v2=U.loc[u, 'thrust_m0'], z_v2=U.loc[u, 'thrust_z'],
                        m0_zback=m, z_added_back=' '.join(back)))
pd.DataFrame(out).to_csv(f'{C.TB}/guard2_zvariant.csv', index=False)
print(pd.DataFrame(out).assign(r=lambda x: x.m0_zback / x.m0_v2).r.describe())
