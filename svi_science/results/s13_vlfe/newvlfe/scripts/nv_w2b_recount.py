#!/usr/bin/env python
"""W2b 2018 re-count (coordinator 2026-10-16 ~10:55 am PT): new-ETS B1/B2 counts by z bin, real vs polarity-control
expectation, before (real/polB) and after W2b (realW/polBW); western group (48.5-48.9 N, lon <= -124.2) LFE coincidence.
-> tables/nv_w2b_zbins.csv, tables/nv_w2b_west.csv"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import nv_summary as SM, nv_survivors as SV
import cat_common as K
NV = SM.NV
V = pd.read_csv(f'{NV}/tables/nv_validation.csv'); valid = set(V[V.valid == True].test)   # noqa: E712
b1 = [c for c in ('t_nfar3', 't_share', 't_vr012', 't_zloso5') if c in valid]
L = SM.template_locs()


def score(run):
    D = pd.read_parquet(f'{NV}/data/bat_2018_{run}.parquet')
    for c in ('travelling', 'grid_edge', 'hf_loud'):
        if c in D:
            D[c] = D[c].astype(object).where(D[c].notna(), np.nan)
    D = K.apply_rules(D); D = pd.concat([D, SM.added_tests(D)], axis=1)
    D['B1'] = D.battery_e2.fillna(False).astype(bool) & D[b1].all(axis=1)
    D['B2'] = D.B1 & D.t_ring
    dd = K.dkm(D.lat.values[:, None], D.lon.values[:, None], L.lat.values[None, :], L.lon.values[None, :])
    D['dist_tpl_km'] = dd.min(1); D['nearest_tpl'] = L.src.values[dd.argmin(1)]
    D['loc'] = np.where(D.dist_tpl_km <= 15, 'known', 'new')
    D['time_pt'] = pd.to_datetime(D.epoch, unit='s').dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.strftime('%Y-%m-%d %H:%M:%S %Z')
    lv = json.load(open(f'{NV}/data/live_2018_{run}.json'))
    return D, lv['live_yr_per_catalogue']['ETS'] * lv['n_catalogues']


rows = []
S = {}
for tag, (r, c) in {'before W2b': ('real', 'polB'), 'after W2b': ('realW', 'polBW')}.items():
    R, lr = score(r); C, lc = score(c); S[tag] = (R, C)
    for tier in ('B1', 'B2'):
        for loc in ('new', 'known'):
            for lo, hi in ((5, 7), (7, 10), (10, 1e9), (5, 1e9)):
                f = lambda X: int((X[tier] & (X.period == 'ETS') & (X['loc'] == loc) & (X.z_best >= lo) & (X.z_best < hi)).sum())
                nr, nc = f(R), f(C); ex = nc / lc * lr
                rows.append(dict(run=tag, tier=tier, loc=loc, z=f'{lo}-{hi if hi < 1e8 else "inf"}', real=nr, polB_n=nc,
                                 expected=round(ex, 1), excess=round(nr - ex, 1), real_live_yr=round(lr, 3), ctrl_live_yr=round(lc, 3)))
Z = pd.DataFrame(rows); Z.to_csv(f'{NV}/tables/nv_w2b_zbins.csv', index=False)
pd.set_option('display.width', 220); print(Z.to_string())
# western group LFE coincidence
X, LF, Rs, WA, F, P = SV.lfe_tables()
rng = np.random.default_rng(20261016)
out = []
for tag, (R, C) in S.items():
    for kind, T in (('real', R), ('polB', C)):
        w = T[T.B2 & (T.period == 'ETS') & (T['loc'] == 'new') & (T.lat >= 48.45) & (T.lat <= 48.95) & (T.lon <= -124.2)]
        for e in w.itertuples():
            ep = float(e.epoch)
            Fn = F[K.dkm(F.lat.values, F.lon.values, e.lat, e.lon) <= 20]
            sv = set(pd.to_numeric(Fn[Fn.region == 'S_Vancouver_Is'].family_id, errors='coerce').dropna().astype(int))
            wa = set(Fn[Fn.region == 'Royer_NW_WA'].family_id.astype(str))
            tt = np.sort(np.r_[Rs[Rs.fam.isin(sv)].t.values, WA[WA.fam.isin(wa)].t.values])
            obs = int(X.cnt(tt, [ep], 600)[0]); null = X.cnt(tt, np.floor(ep / 86400) * 86400 + rng.uniform(0, 86400, 1000), 600)
            n_any = int(X.cnt(LF['svi075'], [ep], 600)[0] + X.cnt(LF['wa_nw'], [ep], 600)[0])
            out.append(dict(run=tag, kind=kind, variant=e.variant, time_pt=e.time_pt, lat=e.lat, lon=e.lon, z_best=e.z_best, z_loso=e.z_loso,
                            ratio_far=e.ratio_far, vr_thrust=e.vr_thrust, n_lfe_any_10m=n_any, patch_lfe_10m=obs,
                            patch_null_mean=float(null.mean()), patch_p=float(np.mean(null >= obs)) if obs else 1.0, epoch=ep))
Wd = pd.DataFrame(out); Wd.to_csv(f'{NV}/tables/nv_w2b_west.csv', index=False)
for (tag, kind), g in Wd.groupby(['run', 'kind'], sort=False):
    print(tag, kind, 'n', len(g), 'with LFEs', int((g.n_lfe_any_10m > 0).sum()), 'patch p<=0.05', int((g.patch_p <= 0.05).sum()))
print(Wd[(Wd.kind == 'real')].round(3).drop(columns='epoch').to_string())
