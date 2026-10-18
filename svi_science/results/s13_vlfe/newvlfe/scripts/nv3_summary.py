#!/usr/bin/env python
"""NV3 (PLAN Amendment NV3): ladder (8 checks), earthquake-wave flag, known/new, separation pooled over 2011, 2013, 2016,
2018, 2019 (all W2b), raw figures. -> tables/nv3_rows.parquet, nv3_separation.csv, nv3_ladder_counts.csv, nv3_summary.json,
figures/nv3_separation.png, nv3_rank.png, nv3_map_level8.png"""
import os, sys, json
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import nv_summary as SM
import cat_common as K
from scipy.stats import poisson
NV = SM.NV; S13 = SM.S13
YR = {y: ['realW', 'polBW', 'shA3W', 'shA7W', 'shA11W'] for y in (2011, 2013, 2016, 2019)}
YR[2018] = ['realW', 'polBW', 'shA3', 'shA7', 'shA11']
CHECKS = ['pass_a', 'pass_b', 'pass_c', 'pass_e2', 'pass_f', 't_nfar3', 't_share', 't_ring']
ZS = [5, 7, 10, 15, 20]


def eq_table(y):
    name = {2011: 'scan2_eq_events.csv', 2012: 'scan3_eq_events.csv'}.get(y, f'scan{y}_eq_events.csv')
    A = pd.read_csv(f'{S13}/tables/{name}')
    return np.sort(A[A.mag >= 4.5].t.values.astype(float))


def load():
    parts, live = [], []
    for y, runs in YR.items():
        for r in runs:
            f = f'{NV}/data/bat_{y}_{r}.parquet'
            if not os.path.isfile(f):
                print('missing', f); continue
            D = pd.read_parquet(f); D['run'] = r; D['year'] = y; parts.append(D)
            L = json.load(open(f'{NV}/data/live_{y}_{r}.json'))
            for p, v in L['live_yr_per_catalogue'].items():
                live.append(dict(year=y, run=r, period=p, live_yr=v * L['n_catalogues'], n_cat=L['n_catalogues']))
    return pd.concat(parts, ignore_index=True), pd.DataFrame(live)


def main():
    D, LV = load()
    for c in ('travelling', 'grid_edge', 'hf_loud'):
        if c in D:
            D[c] = D[c].astype(object).where(D[c].notna(), np.nan)
    D = K.apply_rules(D)
    D = pd.concat([D, SM.added_tests(D)], axis=1)
    full = D.n_sta.notna()
    fl = D.pass_f_light.fillna(False).astype(bool) if 'pass_f_light' in D else full
    D['level'] = D[CHECKS].fillna(False).astype(bool).sum(axis=1)
    D['n_missing'] = np.where(full, 0, np.where(D.gate & ~fl, 5, 6))
    D['full_eval'] = full
    # earthquake-wave flag
    eqw = (D.emerg_far_med >= D.emerg_med_near).fillna(False).values
    m45 = np.zeros(len(D), bool)
    for y in YR:
        t = eq_table(y); idx = np.where(D.year.values == y)[0]; ep = D.epoch.values[idx]
        lo = np.searchsorted(t, ep - 5400); hi = np.searchsorted(t, ep + 60)
        m45[idx] = hi > lo
    D['eq_far_ge_near'] = eqw; D['eq_m45_90min'] = m45; D['eq_wave'] = eqw | m45
    L = SM.template_locs()
    dd = K.dkm(D.lat.values[:, None], D.lon.values[:, None], L.lat.values[None, :], L.lon.values[None, :])
    D['dist_tpl_km'] = dd.min(1); D['nearest_tpl'] = L.src.values[dd.argmin(1)]
    D['loc'] = np.where(D.dist_tpl_km <= 15, 'known', 'new')
    D['kind'] = np.where(D.run == 'realW', 'real', 'control')
    D['ctype'] = np.where(D.run == 'realW', 'real', np.where(D.run == 'polBW', 'polB', 'shift'))
    D['time_pt'] = pd.to_datetime(D.epoch, unit='s').dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.strftime('%Y-%m-%d %H:%M:%S %Z')
    D.to_parquet(f'{NV}/tables/nv3_rows.parquet')
    LV['ctype'] = np.where(LV.run == 'realW', 'real', np.where(LV.run == 'polBW', 'polB', 'shift'))
    rows = []
    for y in list(YR) + ['pooled']:
        Dy = D if y == 'pooled' else D[D.year == y]; Ly = LV if y == 'pooled' else LV[LV.year == y]
        for per in ('ETS', 'between', 'other'):
            lt = Ly[Ly.period == per].groupby('ctype').live_yr.sum()
            for lev in (6, 7, 8):
                for z in ZS:
                    for loc in ('known', 'new'):
                        for eqx in (False, True):
                            m = (Dy.level >= lev) & (Dy.z_best >= z) & (Dy.period == per) & (Dy['loc'] == loc)
                            if eqx:
                                m &= ~Dy.eq_wave
                            o = dict(year=y, period=per, level=lev, z=z, loc=loc, eq_wave_excluded=eqx)
                            nr = int((m & (Dy.ctype == 'real')).sum()); lr = lt.get('real', 0)
                            o.update(real_n=nr, real_live_yr=lr, real_rate=nr / lr if lr else np.nan)
                            for ct in ('polB', 'shift', 'control'):
                                n = int((m & ((Dy.ctype != 'real') if ct == 'control' else (Dy.ctype == ct))).sum())
                                l = lt.drop('real', errors='ignore').sum() if ct == 'control' else lt.get(ct, 0)
                                o.update({f'{ct}_n': n, f'{ct}_live_yr': l, f'{ct}_rate': n / l if l else np.nan})
                            rc = max(o['control_n'], 0.5) / o['control_live_yr'] if o['control_live_yr'] else np.nan
                            rp = max(o['polB_n'], 0.5) / o['polB_live_yr'] if o['polB_live_yr'] else np.nan
                            o.update(ratio_pooled=o['real_rate'] / rc if rc else np.nan, ratio_polB=o['real_rate'] / rp if rp else np.nan,
                                     excess_vs_pooled=nr - rc * lr, excess_vs_polB=nr - rp * lr,
                                     p_pooled=float(poisson.sf(nr - 1, rc * lr)) if nr else 1.0,
                                     p_polB=float(poisson.sf(nr - 1, rp * lr)) if nr else 1.0,
                                     separates=bool(nr >= 3 and o['real_rate'] >= 3 * rc and o['control_rate'] <= 1))
                            rows.append(o)
    SP = pd.DataFrame(rows); SP.to_csv(f'{NV}/tables/nv3_separation.csv', index=False)
    LC = D.groupby(['year', 'ctype', 'period', 'loc', 'level']).size().rename('n').reset_index()
    LC.to_csv(f'{NV}/tables/nv3_ladder_counts.csv', index=False)
    po = SP[(SP.year == 'pooled') & (SP.period.isin(['ETS', 'between']))]
    summ = dict(n_rows=len(D), live=LV.groupby(['ctype', 'period']).live_yr.sum().round(3).reset_index().to_dict('records'),
                separating=po[po.separates][['period', 'level', 'z', 'loc', 'eq_wave_excluded', 'real_n', 'control_n', 'ratio_pooled']].to_dict('records'),
                level8_new=po[(po.level == 8) & (po['loc'] == 'new') & (po.z == 5)][['period', 'eq_wave_excluded', 'real_n', 'polB_n', 'shift_n', 'real_rate',
                                                                                     'polB_rate', 'shift_rate', 'excess_vs_polB', 'p_polB']].round(3).to_dict('records'))
    json.dump(summ, open(f'{NV}/tables/nv3_summary.json', 'w'), indent=1, default=str)
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 300)
    print(po[(po.eq_wave_excluded == False) & (po.level.isin([7, 8]))][['period', 'level', 'z', 'loc', 'real_n', 'real_rate', 'polB_n', 'polB_rate', 'shift_n',   # noqa: E712
                                                                         'shift_rate', 'ratio_polB', 'excess_vs_polB', 'p_polB', 'separates']].round(3).to_string())
    figs(D, SP, LV)


def figs(D, SP, LV):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    po = SP[(SP.year == 'pooled')]
    fig, ax = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    for i, loc in enumerate(('known', 'new')):
        for j, lev in enumerate((6, 7, 8)):
            a = ax[i, j]
            for per, col in (('ETS', 'tab:red'), ('between', 'tab:blue')):
                q = po[(po['loc'] == loc) & (po.level == lev) & (po.period == per) & (po.eq_wave_excluded == False)]   # noqa: E712
                a.plot(q.z, q.real_rate, 'o-', color=col, label=f'real {per}')
                a.plot(q.z, q.polB_rate, 's--', color=col, mfc='none', label=f'polarity controls {per}')
                a.plot(q.z, q.shift_rate, '^:', color=col, mfc='none', label=f'shift controls {per}')
            a.axhline(1, color='k', lw=0.5, ls=':'); a.set_yscale('symlog', linthresh=0.1); a.set_title(f'level ≥ {lev}, {loc}')
            if i == 1:
                a.set_xlabel('z threshold')
            if j == 0:
                a.set_ylabel('events per live year')
            a.legend(fontsize=6)
    fig.suptitle('NV3: real vs scrambled controls, 2011/2013/2016/2018/2019 pooled (W2b), ladder levels')
    fig.tight_layout(); fig.savefig(f'{NV}/figures/nv3_separation.png', dpi=120); plt.close(fig)
    fig, ax = plt.subplots(2, 5, figsize=(20, 7))
    for j, y in enumerate(sorted(YR)):
        for i, per in enumerate(('ETS', 'between')):
            a = ax[i, j]; x = D[(D.year == y) & (D.level == 8) & (D['loc'] == 'new') & (D.period == per)]
            lt = LV[(LV.year == y) & (LV.period == per)].groupby('ctype').live_yr.sum()
            zs = np.r_[5:30.5:0.5]
            a.plot(zs, [((x.ctype == 'real') & (x.z_best >= z)).sum() for z in zs], 'o-', ms=3, color='tab:red', label='real')
            for ct, col in (('polB', 'tab:blue'), ('shift', 'tab:gray')):
                if lt.get(ct, 0):
                    a.plot(zs, [((x.ctype == ct) & (x.z_best >= z)).sum() / lt[ct] * lt.get('real', 0) for z in zs], 's--', ms=3, mfc='none',
                           color=col, label=f'expected from {ct}')
            a.set_yscale('symlog', linthresh=1); a.set_title(f'{y} {per}: level 8, new'); a.legend(fontsize=6); a.set_xlabel('z ≥')
    fig.tight_layout(); fig.savefig(f'{NV}/figures/nv3_rank.png', dpi=110); plt.close(fig)
    fig, ax = plt.subplots(1, 2, figsize=(13, 6.5), sharex=True, sharey=True)
    F = K.families()
    for a, kind in zip(ax, ('real', 'control')):
        a.scatter(F.lon, F.lat, s=3, c='0.8')
        x = D[(D.level == 8) & (D.kind == kind)]
        jit = np.random.default_rng(1).normal(0, 0.012, (len(x), 2))
        sc = a.scatter(x.lon + jit[:, 0], x.lat + jit[:, 1], s=14, c=x.year, cmap='viridis', vmin=2011, vmax=2019)
        a.set_title(f'level-8 {kind} rows (n={len(x)}), coloured by year')
    fig.colorbar(sc, ax=ax); fig.savefig(f'{NV}/figures/nv3_map_level8.png', dpi=110); plt.close(fig)


if __name__ == '__main__':
    main()
