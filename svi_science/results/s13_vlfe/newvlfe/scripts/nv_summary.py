#!/usr/bin/env python
"""newvlfe stage 3 (PLAN.md §2-3): tiers, validation of the added tests on the catalogue's production tier events,
known/new location class, separation tables (real vs scrambled controls, ETS / between), raw figures.
-> tables/nv_rows_scored.parquet, nv_validation.csv, nv_separation.csv, nv_separation_pooled.csv, nv_verdict.json,
   figures/nv_z_vs_ratio.png, nv_map.png, nv_separation.png"""
import os, sys, json, glob
os.environ['OMP_NUM_THREADS'] = '1'
S13 = '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe'
NV = f'{S13}/newvlfe'
sys.path.insert(0, f'{S13}/catalogue/scripts'); sys.path.insert(0, f'{S13}/scripts'); sys.path.insert(0, f'{S13}/hg_repro/scripts')
import numpy as np, pandas as pd
import cat_common as K

ZS = [5, 7, 10, 15, 20]
TIERS = ['B0', 'B1', 'B2']
RUNS = ['real', 'polB', 'shA3', 'shA7', 'shA11']
YEARS = [2013, 2018]
CAT = f'{S13}/catalogue'


def template_locs():
    import hg_common as H
    T = H.templates()
    E = pd.read_csv(f'{CAT}/tables/e3_events.csv')
    R = pd.read_csv(f'{CAT}/tables/prod2018_relay_events.csv')
    L = pd.concat([pd.DataFrame(dict(src='HG_T' + T.template.astype(str), lat=T.lat, lon=T.lon)),
                   pd.DataFrame(dict(src='C' + E.cluster.astype(str), lat=E.lat_best, lon=E.lon_best)),
                   pd.DataFrame(dict(src='RC' + R.cluster.astype(str), lat=R.lat_best, lon=R.lon_best))], ignore_index=True)
    return L


def added_tests(D):
    o = pd.DataFrame(index=D.index)
    o['t_nfar3'] = D.n_far >= 3
    o['t_share'] = D.share_excess <= 0.25
    o['t_vr012'] = D.vr_thrust >= 0.12
    o['t_zloso5'] = D.z_loso >= 5
    o['t_ring'] = ~((D.pw2_n >= 8) & (D.pw2_vr >= 0.5)).fillna(False)
    return o.fillna(False).astype(bool)


def validation():
    """PLAN §2: each added test must keep >= 80 % of the production tier events of 2013 + 2018 (B0 pass by construction)."""
    out = []
    for y, d in ((2013, 'prodq2013'), (2018, 'prodr2018')):
        f = f'{CAT}/data/{d}/rows_scored_v2.parquet'
        R = pd.read_parquet(f)
        x = R[(R.source == 'real') & R.tier].sort_values('epoch')
        x = x.assign(g=(x.i.diff().fillna(1e9) > 30).cumsum()).sort_values('zmad').groupby('g').tail(1)
        x['year'] = y; x['src_file'] = d
        out.append(x)
    X = pd.concat(out, ignore_index=True)
    T = added_tests(X)
    V = pd.DataFrame([dict(test=c, n=len(T), kept=int(T[c].sum()), frac=float(T[c].mean())) for c in T.columns])
    V['valid'] = V.frac >= 0.8
    ok = [c for c in ('t_nfar3', 't_share', 't_vr012', 't_zloso5') if bool(V.set_index('test').loc[c, 'valid'])]
    allB1 = T[ok].all(axis=1) if ok else pd.Series(True, index=T.index)
    ring = T.t_ring if bool(V.set_index('test').loc['t_ring', 'valid']) else True
    V = pd.concat([V, pd.DataFrame([dict(test='B1 (valid tests)', n=len(T), kept=int(allB1.sum()), frac=float(allB1.mean()), valid=np.nan),
                                    dict(test='B2 (valid tests)', n=len(T), kept=int((allB1 & ring).sum()),
                                         frac=float((allB1 & ring).mean()), valid=np.nan)])], ignore_index=True)
    return V, X


def load():
    parts, live = [], []
    for y in YEARS:
        for r in RUNS:
            f = f'{NV}/data/bat_{y}_{r}.parquet'
            if not os.path.isfile(f):
                continue
            parts.append(pd.read_parquet(f))
            L = json.load(open(f'{NV}/data/live_{y}_{r}.json'))
            for p, v in L['live_yr_per_catalogue'].items():
                live.append(dict(year=y, run=r, period=p, live_yr=v * L['n_catalogues'], n_cat=L['n_catalogues']))
    return pd.concat(parts, ignore_index=True), pd.DataFrame(live)


def main():
    V, PT = validation()
    V.to_csv(f'{NV}/tables/nv_validation.csv', index=False)
    valid = set(V[(V.valid == True)].test)              # noqa: E712
    D, LV = load()
    for c in ('travelling', 'grid_edge', 'hf_loud'):
        if c in D:
            D[c] = D[c].astype(object).where(D[c].notna(), np.nan)
    D = K.apply_rules(D)
    T = added_tests(D)
    D = pd.concat([D, T], axis=1)
    b1 = [c for c in ('t_nfar3', 't_share', 't_vr012', 't_zloso5') if c in valid]
    D['B0'] = D.battery_e2.fillna(False).astype(bool)
    D['B1'] = D.B0 & D[b1].all(axis=1)
    D['B2'] = D.B1 & (D.t_ring if 't_ring' in valid else True)
    L = template_locs()
    dd = K.dkm(D.lat.values[:, None], D.lon.values[:, None], L.lat.values[None, :], L.lon.values[None, :])
    D['dist_tpl_km'] = dd.min(axis=1); D['nearest_tpl'] = L.src.values[dd.argmin(axis=1)]
    D['loc'] = np.where(D.dist_tpl_km <= 15, 'known', 'new')
    D['loc30'] = np.where(D.dist_tpl_km <= 30, 'known30', 'new30')
    D['kind'] = np.where(D.run == 'real', 'real', 'control')
    D['time_pt'] = pd.to_datetime(D.epoch, unit='s').dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.strftime('%Y-%m-%d %H:%M:%S %Z')
    D.to_parquet(f'{NV}/tables/nv_rows_scored.parquet')
    # ---- separation
    rows = []
    for y in YEARS + ['pooled']:
        Dy = D if y == 'pooled' else D[D.year == y]
        Ly = LV if y == 'pooled' else LV[LV.year == y]
        for per in ('ETS', 'between', 'other'):
            lr = Ly[(Ly.period == per) & (Ly.run == 'real')].live_yr.sum()
            lc = Ly[(Ly.period == per) & (Ly.run != 'real')].live_yr.sum()
            ncat = Ly[(Ly.period == per) & (Ly.run != 'real')].n_cat.sum()
            for tier in TIERS:
                for z in ZS:
                    for loc in ('known', 'new', 'all'):
                        m = Dy[tier] & (Dy.z_best >= z) & (Dy.period == per)
                        if loc != 'all':
                            m &= Dy['loc'] == loc
                        nr = int((m & (Dy.kind == 'real')).sum()); nc = int((m & (Dy.kind == 'control')).sum())
                        rr = nr / lr if lr > 0 else np.nan
                        rc = nc / lc if lc > 0 else np.nan
                        rcf = max(nc, 0.5) / lc if lc > 0 else np.nan
                        from scipy.stats import poisson
                        p = float(poisson.sf(nr - 1, rcf * lr)) if nr > 0 else 1.0
                        nhf = int((m & (Dy.kind == 'real') & Dy.hf_loud.fillna(False).astype(bool)).sum())
                        rows.append(dict(year=y, period=per, tier=tier, z=z, loc=loc, real_n=nr, real_live_yr=lr, real_rate=rr,
                                         ctrl_n=nc, ctrl_catalogues=ncat, ctrl_live_yr=lc, ctrl_rate=rc, ctrl_rate_floor=rcf,
                                         ratio=rr / rcf if rcf and rcf > 0 else np.nan, poisson_p=p, real_hf_loud=nhf,
                                         passes=bool(nr >= 3 and rr >= 3 * rcf and rc <= 1)))
    SP = pd.DataFrame(rows)
    # next-z robustness column
    SP['passes_next_z'] = False
    for k, r in SP.iterrows():
        if r.z == ZS[-1]:
            continue
        zn = ZS[ZS.index(r.z) + 1]
        q = SP[(SP.year == r.year) & (SP.period == r.period) & (SP.tier == r.tier) & (SP['loc'] == r['loc']) & (SP.z == zn)]
        SP.loc[k, 'passes_next_z'] = bool(q.passes.iloc[0]) if len(q) else False
    SP.to_csv(f'{NV}/tables/nv_separation.csv', index=False)
    PO = SP[SP.year == 'pooled']
    PO.to_csv(f'{NV}/tables/nv_separation_pooled.csv', index=False)
    verdict = {}
    for loc in ('known', 'new'):
        for per in ('ETS', 'between'):
            q = PO[(PO['loc'] == loc) & (PO.period == per)]
            pas = q[q.passes]
            verdict[f'{loc}_{per}'] = dict(passes=bool(len(pas)), passing_tiers=[f'{r.tier} z>={r.z}' for r in pas.itertuples()],
                                           best=q.sort_values(['ratio', 'real_n'], ascending=False).head(3)[
                                               ['tier', 'z', 'real_n', 'real_rate', 'ctrl_n', 'ctrl_rate', 'ratio', 'poisson_p']].round(3).to_dict('records'))
    verdict['added_tests_valid'] = sorted(valid); verdict['B1_tests_used'] = b1
    verdict['validation'] = V.round(3).to_dict('records')
    verdict['live'] = LV.groupby(['run', 'period']).live_yr.sum().round(3).reset_index().to_dict('records')
    json.dump(verdict, open(f'{NV}/tables/nv_verdict.json', 'w'), indent=1, default=str)
    pd.set_option('display.width', 250); pd.set_option('display.max_rows', 400)
    print(V.round(3).to_string())
    print(PO[PO['loc'] != 'all'][['period', 'tier', 'z', 'loc', 'real_n', 'real_rate', 'ctrl_n', 'ctrl_rate', 'ratio', 'poisson_p', 'real_hf_loud', 'passes']].round(3).to_string())
    print(json.dumps({k: v for k, v in verdict.items() if k.endswith(('ETS', 'between'))}, indent=1, default=str))
    figs(D, PO, L)


def figs(D, PO, L):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    G = D[D.gate & D.pass_f_light.fillna(False) & D.ratio_far.notna()]
    # 1. raw: z vs ratio_far, real vs controls, per period x loc
    fig, ax = plt.subplots(2, 2, figsize=(12, 9), sharex=True, sharey=True)
    for i, per in enumerate(('ETS', 'between')):
        for j, loc in enumerate(('known', 'new')):
            a = ax[i, j]
            x = G[(G.period == per) & (G['loc'] == loc)]
            c = x[x.kind == 'control']; r = x[x.kind == 'real']
            a.scatter(c.z_best, c.ratio_far, s=5, c='0.6', alpha=0.4, label=f'scrambled controls (14 catalogues) n={len(c)}')
            sc = a.scatter(r.z_best, r.ratio_far, s=9, c=r.vr_thrust.clip(-0.2, 0.8), cmap='viridis', label=f'real n={len(r)}')
            b2 = x[x.B2]
            a.scatter(b2[b2.kind == 'real'].z_best, b2[b2.kind == 'real'].ratio_far, s=60, facecolors='none', edgecolors='r', label='real, B2')
            a.scatter(b2[b2.kind == 'control'].z_best, b2[b2.kind == 'control'].ratio_far, s=60, marker='s', facecolors='none', edgecolors='k', label='control, B2')
            a.axhline(0.551, color='k', lw=0.6, ls='--')
            a.set_xscale('log'); a.set_yscale('log'); a.set_title(f'{per}, {loc} (≤/> 15 km from a template location)')
            a.legend(fontsize=7, loc='upper right')
            if i == 1:
                a.set_xlabel('scan z (best node)')
            if j == 0:
                a.set_ylabel('ratio_far (median far LHZ / near amplitude)')
    fig.colorbar(sc, ax=ax, label='vr_thrust (real)')
    fig.suptitle('Scan-v2 candidates 2013 + 2018 passing location, earthquake-time and single-station gates')
    fig.savefig(f'{NV}/figures/nv_z_vs_ratio.png', dpi=130); plt.close(fig)
    # 2. map
    fig, ax = plt.subplots(1, 2, figsize=(13, 6.5), sharex=True, sharey=True)
    F = K.families()
    for j, tier in enumerate(('B0', 'B2')):
        a = ax[j]
        a.scatter(F.lon, F.lat, s=3, c='0.8', label='LFE families')
        th = np.linspace(0, 2 * np.pi, 60)
        for r in L.itertuples():
            a.plot(r.lon + 15 / (111.19 * np.cos(np.radians(48.3))) * np.cos(th), r.lat + 15 / 111.19 * np.sin(th), color='tab:blue', lw=0.3)
        x = D[D[tier]]
        c = x[x.kind == 'control']; r = x[x.kind == 'real']
        jit = lambda n: np.random.default_rng(0).normal(0, 0.01, n)
        a.scatter(c.lon + jit(len(c)), c.lat + jit(len(c)), s=12, marker='s', c='k', alpha=0.5, label=f'controls n={len(c)}')
        sc = a.scatter(r.lon + jit(len(r)), r.lat + jit(len(r)), s=18, c=np.log10(r.z_best), cmap='plasma', label=f'real n={len(r)}')
        a.set_title(f'{tier} candidates (z ≥ 5), 2013 + 2018; blue = 15 km around template locations')
        a.legend(fontsize=7); a.set_xlabel('lon')
    ax[0].set_ylabel('lat'); fig.colorbar(sc, ax=ax, label='log10 z (real)')
    fig.savefig(f'{NV}/figures/nv_map.png', dpi=130); plt.close(fig)
    # 3. separation curves
    fig, ax = plt.subplots(2, 3, figsize=(14, 8), sharex=True)
    for i, loc in enumerate(('known', 'new')):
        for j, tier in enumerate(TIERS):
            a = ax[i, j]
            for per, col in (('ETS', 'tab:red'), ('between', 'tab:blue')):
                q = PO[(PO['loc'] == loc) & (PO.tier == tier) & (PO.period == per)]
                a.plot(q.z, q.real_rate, 'o-', color=col, label=f'real {per}')
                a.plot(q.z, q.ctrl_rate, 's--', color=col, mfc='none', label=f'controls {per}')
            a.axhline(1, color='k', lw=0.5, ls=':')
            a.set_yscale('symlog', linthresh=0.1); a.set_title(f'{tier}, {loc}')
            if i == 1:
                a.set_xlabel('z threshold')
            if j == 0:
                a.set_ylabel('events per live year')
            a.legend(fontsize=7)
    fig.suptitle('Real vs scrambled-control rates, scan-only candidates, 2013 + 2018 pooled')
    fig.savefig(f'{NV}/figures/nv_separation.png', dpi=130); plt.close(fig)


if __name__ == '__main__':
    main()
