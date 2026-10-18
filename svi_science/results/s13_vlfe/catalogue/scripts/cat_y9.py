#!/usr/bin/env python
"""Amendment Y9 (Merlin review 5): Q2 completeness, transparent battery' ladder, emergence / amplitude-distance columns,
figures of the level-5 between-ETS rows.
  step1          -> tables/y9_q2_completeness.csv, tables/y9_q2_known.json
  emerge Y       -> data/ranked/y9_emerge_<Y>.parquet (every agreement row, all classes)
  ladder         -> data/ranked/y9_all_<Y>.parquet, tables/catalogue_ladder_<Y>.csv, catalogue_ladder_2010_2026.csv,
                    tables/y9_ladder_table.csv, tables/y9_ladder_completeness.csv, tables/y9_checks.json
  figs           -> figures/y9_l5_<k>.png, tables/y9_level5_between.csv (labels filled by rule; eye labels added by hand)"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

RD = f'{K.DATA}/ranked'
YEARS = list(range(2010, 2019)) + [int(x) for x in os.environ.get('LADDER_EXTRA', '').split(',') if x]   # Y11: 2019+
RUNS = {2010: 'prod', 2011: 'prod', 2012: 'prod', 2013: 'prodq', 2014: 'prod', 2015: 'prod', 2016: 'prodq', 2017: 'prodq',
        2018: 'prodq', **{y: 'prodf' for y in range(2019, 2027)}}
RUNS.update({int(a): b for a, b in (x.split('=') for x in os.environ.get('LADDER_RUNS', '').split(',') if x)})   # e.g. 2021=prodfy (Y12)
NEAR = ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.PFB', 'CN.YOUB', 'UW.LRIV', 'CN.LZB']
REG = ['UW.GNW', 'UW.DOSE', 'UW.SP2', 'TA.D03D', 'TA.B05D', 'UW.WISH', 'US.NLWA']
FAR = ['IU.COR', 'US.NEW', 'US.HAWA', 'US.MSO', 'US.WVOR']
MWS = (3.6, 3.8, 4.0)


def known_lists():
    V = pd.read_csv(f'{K.TB}/e3_events.csv')
    return np.sort(np.r_[V[V.role == 'vetted'].epoch.values, pd.read_csv(f'{K.TB}/pilot2013_ets_events_v2.csv').epoch.values])


def near(ep, ref, w=30):
    k = np.searchsorted(ref, ep)
    d = np.minimum(np.abs(ep - ref[np.clip(k - 1, 0, len(ref) - 1)]), np.abs(ep - ref[np.clip(k, 0, len(ref) - 1)]))
    return d <= w


def inj(y):
    J = pd.read_parquet(f'{K.DATA}/{RUNS[y]}{y}/inj_list.parquet').merge(pd.read_parquet(f'{K.DATA}/{RUNS[y]}{y}/inj_v2.parquet'), on='inj_id', how='left')
    J = K.rescore_inj(J)                                     # Amendment Y14: f'
    J['mf2'] = (J.mf_best_z >= 5.5) & (J.mf_second_z >= 5.5)
    J['bat'] = J.battery_e2.fillna(False).astype(bool)
    J['level'] = sum(J[c].fillna(False).astype(bool).astype(int) for c in ('pass_a', 'pass_b', 'pass_c', 'pass_e2', 'pass_f'))
    J['year'] = y
    return J


# ------------------------------------------------------------------------------------------------ step 1
def step1():
    rows = []
    for y in YEARS:
        J = inj(y); R = J[J.mf2 & J.bat]
        for per in ('ETS', 'between', 'all'):
            for mw in MWS:
                x = R[(R.Mw == mw) & ((R.period == per) if per != 'all' else True)]
                n0 = J[(J.Mw == mw) & ((J.period == per) if per != 'all' else True)]
                rows.append(dict(year=y, run=RUNS[y], period=per, Mw=mw, n_inj=len(n0), n_battery_rec=len(x),
                                 frac_vr_ge_0p12=float((x.vr_thrust >= 0.12).mean()) if len(x) else np.nan,
                                 completeness_battery=len(x) / max(len(n0), 1),
                                 completeness_battery_and_vr=float(((x.vr_thrust >= 0.12)).sum() / max(len(n0), 1))))
    T = pd.DataFrame(rows); T.to_csv(f'{K.TB}/y9_q2_completeness.csv', index=False)
    known = known_lists(); kq = []
    for y in (2011, 2012, 2013):
        A = pd.read_parquet(f'{RD}/r3_all_{y}.parquet'); A = A[A.ctrl == 'real']
        kq.append(A[near(A.epoch.values, known)][['epoch', 'q2', 'vr_thrust', 'far_over_near2', 'ratio_GNW', 'lat_best', 'zmad', 'battery_e2']])
    kq = pd.concat(kq)
    S = dict(known_n=len(kq), q2_pass=int(kq.q2.sum()), vr_ge_0p12=int((kq.vr_thrust >= 0.12).sum()),
             fon2_lt_0p5=int((kq.far_over_near2 < 0.5).sum()),
             gnw_ok=int(np.where(kq.lat_best > 48.1, kq.ratio_GNW <= 0.551, True).sum()), battery_prime=int(kq.battery_e2.fillna(False).sum()))
    json.dump(S, open(f'{K.TB}/y9_q2_known.json', 'w'), indent=1, default=float)
    pd.set_option('display.width', 250)
    print(json.dumps(S)); print(T[T.period == 'all'].pivot_table(index='year', columns='Mw', values='frac_vr_ge_0p12').round(2).to_string())
    print(T[T.period != 'all'].groupby(['period', 'Mw'])[['n_battery_rec']].sum().join(
        T[T.period != 'all'].assign(k=lambda d: d.frac_vr_ge_0p12 * d.n_battery_rec).groupby(['period', 'Mw']).k.sum().rename('n_vr')).assign(
        frac=lambda d: d.n_vr / d.n_battery_rec).round(3).to_string())


# ------------------------------------------------------------------------------------------------ step 5 (emergence)
def emerge(y):
    from scipy.signal import hilbert
    from scipy.stats import spearmanr
    import scan_refs as SR
    cc = SR.all_coords()
    A = pd.read_parquet(f'{RD}/r3_all_{y}.parquet')
    A = A[A.agree_bin == 1]
    SD = f'{K.S13}/data/scan' if y == 2011 else f'{K.S13}/data/scan_{y}'
    refbp = f'{SD}/bp2_ref_v2' if y in (2011, 2012) else f'{SD}/bp2_ref'
    arr = {}
    for s in NEAR + REG + FAR:
        f = f'{SD}/bp_{y}g/{s}.npy' if s in NEAR else f'{refbp}/{s}.npy'
        arr[s] = np.load(f, mmap_mode='r') if os.path.isfile(f) else None
    y0 = pd.Timestamp(f'{y}-01-01').value / 1e9
    out = []
    for r in A.itertuples():
        ts = int(round(r.epoch - y0)); em, amp, dist = {}, {}, {}
        for s, X in arr.items():
            if X is None or ts - 1800 < 0 or ts + 1800 > X.shape[1]:
                continue
            x = np.asarray(X[2, ts - 1800:ts + 1800], float) * 1e9
            if np.isfinite(x).sum() < 0.9 * len(x) or not np.isfinite(x[1800:1980]).all():
                continue
            x = np.nan_to_num(x - np.nanmean(x))
            env = np.abs(hilbert(x)); med = float(np.median(env))
            pk = float(np.max(np.abs(x[1800:1980])))
            amp[s] = pk; em[s] = pk / med if med > 0 else np.nan
            dist[s] = float(K.dkm(cc.loc[s, 'lat'], cc.loc[s, 'lon'], r.lat_best, r.lon_best)) if (s in cc.index and np.isfinite(r.lat_best)) else np.nan
        nf = [s for s in em if s not in FAR]
        netw = [s for s in nf if np.isfinite(dist.get(s, np.nan)) and dist[s] <= 150] or [s for s in nf if s in NEAR]
        ev = np.array([em[s] for s in netw]); allnf = sorted([em[s] for s in nf if np.isfinite(em[s])], reverse=True)
        lead = max(nf, key=lambda s: em[s]) if nf else ''
        # section 12 T3 quantities
        Rk = {s: dist[s] for s in amp if np.isfinite(dist.get(s, np.nan))}
        nm = [s for s in Rk if s not in FAR]
        nr = [s for s in nm if Rk[s] <= 40]
        if len(nr) < 2:
            nr = sorted(nm, key=Rk.get)[:2]
        mid = [s for s in nm if 40 < Rk[s] <= 150 and s not in nr]
        far = [s for s in FAR if s in amp]
        r1 = np.median([amp[s] for s in nr]) / np.median([amp[s] for s in mid]) if nr and mid else np.nan
        r2 = np.median([amp[s] for s in far]) / max(amp[s] for s in nr) if far and nr else np.nan
        allS = [s for s in nr + mid + far if s in Rk]
        rho = spearmanr([Rk[s] for s in allS], [amp[s] for s in allS]).correlation if len(allS) >= 4 else np.nan
        out.append(dict(i=int(r.i), ctrl=r.ctrl, source=r.source, draw=int(r.draw), epoch=r.epoch,
                        emerg_med_near=float(np.nanmedian(ev)) if len(ev) else np.nan, n_near_emerg_ge2=int((ev >= 2).sum()), n_near=len(ev),
                        emerg_max=allnf[0] if allnf else np.nan, emerg_max_over_second=allnf[0] / allnf[1] if len(allnf) > 1 and allnf[1] > 0 else np.nan,
                        emerg_lead_sta=lead, emerg_far_med=float(np.nanmedian([em[s] for s in far])) if far else np.nan,
                        t3_near_mid=r1, t3_far_over_nearmax=r2, t3_spearman=rho, t3_pass=bool(r1 >= 1.5 and r2 <= 0.5 and rho <= -0.4),
                        emerg_json=json.dumps({s: round(em[s], 2) for s in em})))
    pd.DataFrame(out).to_parquet(f'{RD}/y9_emerge_{y}.parquet')
    K.log(f'Y9 emerge {y}: {len(out)} rows', 'y9')


# ------------------------------------------------------------------------------------------------ ladder
def ladder_year(y):
    A = pd.read_parquet(f'{RD}/r3_all_{y}.parquet')
    b = pd.read_parquet(f'{RD}/r2bat_{y}.parquet'); b = K.apply_rules(b)
    cols = ['pass_a', 'pass_b', 'pass_c', 'pass_e', 'pass_e2', 'pass_f', 'a_testable', 'in_head2', 'in_bw5', 'in_reg', 'in_reg3', 'in_m45',
            'battery_e2', 'grid_edge_node']
    b = b[['i', 'travelling', 'z_loso', 'z_best'] + cols].rename(columns={c: c + '_b' for c in ['travelling', 'z_loso', 'z_best'] + cols})
    A = A.merge(b, on='i', how='left')
    ag = A.agree_bin == 1
    meas = dict(a=A.a_testable_b.fillna(False).astype(bool), b=A.travelling_b.notna(), c=A.z_best_b.notna() & np.isfinite(A.z_best_b.astype(float)),
                e=pd.Series(True, index=A.index), f=A.z_loso_b.notna())
    pas = dict(a=A.pass_a_b, b=A.pass_b_b, c=A.pass_c_b, e=A.pass_e2_b, f=A.pass_f_b)
    lvl = np.zeros(len(A), int); miss = np.zeros(len(A), int)
    for k in 'abcef':
        p = (meas[k] & pas[k].fillna(False).astype(bool) & ag).values
        A[f'chk_{k}'] = p; lvl += p.astype(int); miss += (~meas[k].values & ag.values).astype(int) if k != 'e' else 0
    A['level'] = np.where(ag, lvl, -1); A['n_checks_missing'] = np.where(ag, miss, np.nan)
    A['eq_in_window'] = np.where(ag, ~A.pass_e_b.fillna(True).astype(bool), A.get('pass_e').fillna(True).astype(bool) == False)
    A['in_head2_any'] = np.where(ag, A.in_head2_b.fillna(False).astype(bool), A.in_head2.fillna(False).astype(bool))
    st = pd.read_csv(f'{K.TB}/{RUNS[y]}{y}_R_stormdays_v2.csv')
    st2 = pd.concat([pd.read_csv(f'{K.TB}/{RUNS[z]}{z}_R_stormdays_v2.csv') for z in (y - 1, y, y + 1) if z in RUNS and os.path.isfile(f'{K.TB}/{RUNS[z]}{z}_R_stormdays_v2.csv')]).drop_duplicates('day')
    pm = st2.set_index(pd.to_datetime(st2.day)).pct.sort_index()
    mx = pm.rolling('3D', center=True).max()
    A['storm_pct_max_pm1d'] = pd.to_datetime(A.day).map(mx).values
    # Y8c in-sample score with loso fill -1
    import cat_ranked2 as R2
    import cat_ranked3 as R3  # noqa
    A = R2.feats(A)
    M = json.load(open(f'{RD}/score3_model.json'))
    A['score_logit_insample'] = (A[R2.FEATS].values.astype(float) - np.array(M['mean'])) / np.array(M['scale']) @ np.array(M['coef']) + M['intercept']
    e = f'{RD}/y9_emerge_{y}.parquet'
    if os.path.isfile(e):
        E = pd.read_parquet(e).drop(columns=['epoch'])
        A = A.merge(E, on=['i', 'ctrl', 'source', 'draw'], how='left')
    A.to_parquet(f'{RD}/y9_all_{y}.parquet')
    return A


def garwood(n, a=0.05):
    from scipy.stats import chi2
    lo = 0.0 if n == 0 else chi2.ppf(a / 2, 2 * n) / 2
    return lo, chi2.ppf(1 - a / 2, 2 * n + 2) / 2


def ladder():
    import cat_ychecks as YC
    known = known_lists()
    TB = {y: YC.between_years(y, RUNS[y]) for y in YEARS}
    allA, rows, comp = [], [], []
    keep = ['time_pt', 'period', 'episode', 'level', 'n_checks_missing', 'chk_a', 'chk_b', 'chk_c', 'chk_e', 'chk_f', 'zmad', 'template', 'agree',
            'n_templates', 'z_best', 'ratio_far', 'n_far', 'ratio_GNW', 'vr_thrust', 'far_over_near2', 'z_loso', 'loso_frac', 'travelling',
            'eq_in_window', 'in_head2_any', 'head2_maxmag', 'grid_edge', 'storm_pct_max_pm1d', 'lat_best', 'lon_best', 'dep_best', 'nearest_family',
            'dist_family_km', 'lfe_status', 'n_lfe_10m', 'emerg_med_near', 'n_near_emerg_ge2', 'n_near', 'emerg_max', 'emerg_max_over_second',
            'emerg_lead_sta', 'emerg_far_med', 't3_near_mid', 't3_far_over_nearmax', 't3_spearman', 't3_pass', 'strong_tier', 'q2',
            'score_logit_insample', 'emerg_json', 'epoch', 'year', 'run']
    for y in YEARS:
        A = ladder_year(y); A['year'] = y; A['run'] = RUNS[y]; allA.append(A)
        R = A[A.ctrl == 'real'].sort_values(['level', 'zmad'], ascending=False)
        R[[c for c in keep if c in R]].to_csv(f'{K.TB}/catalogue_ladder_{y}.csv', index=False)
        nn = A[A.ctrl == 'noise_template'].template.nunique()
        J = inj(y)
        for per in ('ETS', 'between'):
            for L in range(6):
                r = R[(R.period == per) & (R.level >= L)]
                n_ = A[(A.ctrl == 'noise_template') & (A.period == per) & (A.level >= L)]
                s_ = A[(A.ctrl == 'scrambled') & (A.period == per) & (A.level >= L)]
                rows.append(dict(year=y, period=per, level_ge=L, real=len(r), real_days=r.day.nunique(),
                                 noise_x11=11 * len(n_) / nn, noise_raw=len(n_), n_noise_tpl=nn,
                                 A1_per_draw=float((s_.source == 'A1').sum() / 2), A2_per_draw=float((s_.source == 'A2').sum() / 2),
                                 live_yr=TB[y] if per == 'between' else np.nan))
                for mw in MWS:
                    jj = J[(J.period == per) & (J.Mw == mw)]
                    comp.append(dict(year=y, period=per, level_ge=L, Mw=mw, n=len(jj), rec=int((jj.mf2 & (jj.level >= L)).sum())))
    T = pd.DataFrame(rows); C = pd.DataFrame(comp)
    P = T.groupby(['period', 'level_ge']).sum(numeric_only=True).reset_index().assign(year='pooled')
    P['real_days'] = [pd.concat(allA)[lambda d: (d.ctrl == 'real') & (d.period == p) & (d.level >= L)].assign(dd=lambda d: d.year.astype(str) + d.day).dd.nunique()
                      for p, L in zip(P.period, P.level_ge)]
    T = pd.concat([T, P], ignore_index=True)
    Cp = C.groupby(['period', 'level_ge', 'Mw'])[['n', 'rec']].sum().reset_index().assign(year='pooled')
    C = pd.concat([C, Cp], ignore_index=True); C['completeness'] = C.rec / C.n
    for mw in MWS:
        T = T.merge(C[C.Mw == mw][['year', 'period', 'level_ge', 'completeness']].rename(columns={'completeness': f'compl_{mw}'}),
                    on=['year', 'period', 'level_ge'], how='left')
    gw = [garwood(int(n)) for n in T.real]
    T['real_lo95'] = [g[0] for g in gw]; T['real_hi95'] = [g[1] for g in gw]
    T['excess_vs_noise'] = T.real - T.noise_x11
    T['excess_lo95'] = T.real_lo95 - T.noise_x11; T['excess_hi95'] = T.real_hi95 - T.noise_x11
    T['excess_vs_maxctrl'] = T.real - T[['noise_x11', 'A1_per_draw', 'A2_per_draw']].max(axis=1)
    T['excess_rate_per_yr'] = T.excess_vs_noise / T.live_yr
    T['excess_rate_over_compl38'] = T.excess_rate_per_yr / T['compl_3.8']
    T.to_csv(f'{K.TB}/y9_ladder_table.csv', index=False); C.to_csv(f'{K.TB}/y9_ladder_completeness.csv', index=False)
    ALL = pd.concat(allA, ignore_index=True)
    R = ALL[ALL.ctrl == 'real'].sort_values(['level', 'zmad'], ascending=False)
    R[[c for c in keep if c in R]].to_csv(f'{K.TB}/catalogue_ladder_2010_2026.csv', index=False)
    rb = R[R.period == 'between'].head(50)
    kn = R[near(R.epoch.values, known)]
    S = dict(top50_between_in_head2=int(rb.in_head2_any.sum()), top50_between_grid_edge=int(rb.grid_edge.fillna(False).astype(bool).sum()),
             top50_between_eq_in_window=int(rb.eq_in_window.sum()),
             known_n=len(kn), known_at_max_allowed=int((kn.level == 5 - kn.n_checks_missing).sum()),
             known_levels=kn.level.value_counts().to_dict(), known_missing=kn.n_checks_missing.value_counts().to_dict(),
             known_not_at_max=kn[kn.level < 5 - kn.n_checks_missing][['time_pt', 'level', 'n_checks_missing', 'chk_a', 'chk_b', 'chk_c', 'chk_e', 'chk_f', 'zmad']].astype(str).to_dict('records'))
    if 'emerg_med_near' in kn:
        thr = float(np.nanquantile(kn.emerg_med_near, 0.15))
        d23 = R[R.time_pt.str.startswith('2014-01-23')]
        S.update(emerg_threshold_known_p15=thr, known_emerg=kn.emerg_med_near.round(2).tolist(),
                 known_frac_ge_thr=float((kn.emerg_med_near >= thr).mean()),
                 day_2014_01_23_rows=len(d23), day_2014_01_23_frac_ge_thr=float((d23.emerg_med_near >= thr).mean()),
                 day_2014_01_23_level5_frac_ge_thr=float((d23[d23.level == 5].emerg_med_near >= thr).mean()) if (d23.level == 5).any() else None,
                 known_t3_pass=int(kn.t3_pass.sum()), known_spearman_med=float(kn.t3_spearman.median()),
                 known_max_over_second_med=float(kn.emerg_max_over_second.median()))
    json.dump(S, open(f'{K.TB}/y9_checks.json', 'w'), indent=1, default=str)
    pd.set_option('display.width', 250)
    print(json.dumps({k: v for k, v in S.items() if k != 'known_not_at_max'}, indent=1, default=str))
    for p in ('between', 'ETS'):
        print(T[(T.year == 'pooled') & (T.period == p)][['level_ge', 'real', 'real_days', 'noise_x11', 'noise_raw', 'A1_per_draw', 'A2_per_draw',
                                                         'compl_3.6', 'compl_3.8', 'compl_4.0', 'real_lo95', 'real_hi95', 'excess_vs_noise',
                                                         'excess_rate_per_yr', 'excess_rate_over_compl38']].round(2).to_string(index=False))
    print(T[(T.year != 'pooled') & (T.period == 'between') & (T.level_ge == 5)][['year', 'real', 'real_days', 'noise_x11', 'noise_raw', 'A1_per_draw',
                                                                                 'A2_per_draw', 'compl_3.8']].round(2).to_string(index=False))


# ------------------------------------------------------------------------------------------------ figures of level-5 between rows
def figs(years=None, pre='y9_l5', outname='y9_level5_between.csv'):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import scan_refs as SR
    cc = SR.all_coords()
    R = pd.read_csv(f'{K.TB}/catalogue_ladder_2010_2026.csv')
    t = R[(R.period == 'between') & (R.level == 5)].sort_values('epoch').reset_index(drop=True)
    if years:
        t = t[t.year.isin(years)].reset_index(drop=True)
    thr = json.load(open(f'{K.TB}/y9_checks.json')).get('emerg_threshold_known_p15', 2.0)
    lab = []
    for k, e in t.iterrows():
        sp = (e.emerg_max_over_second >= 3)
        wf = (e.emerg_med_near < thr)
        lab.append('spike' if sp else 'wavefield' if wf else 'candidate')
    t['label_rule'] = lab
    t.to_csv(f'{K.TB}/{outname}', index=False)
    for k, e in t.iterrows():
        y = int(e.year); SD = f'{K.S13}/data/scan' if y == 2011 else f'{K.S13}/data/scan_{y}'
        refbp = f'{SD}/bp2_ref_v2' if y in (2011, 2012) else f'{SD}/bp2_ref'
        ts = int(round(e.epoch - pd.Timestamp(f'{y}-01-01').value / 1e9)); a, n = ts - 1800, 3600
        tr = []
        for s in NEAR + REG + FAR:
            f = f'{SD}/bp_{y}g/{s}.npy' if s in NEAR else f'{refbp}/{s}.npy'
            if not os.path.isfile(f):
                continue
            X = np.load(f, mmap_mode='r')
            if a < 0 or a + n > X.shape[1]:
                continue
            x = np.asarray(X[2, a:a + n], float) * 1e9
            if np.isfinite(x).sum() > 0.8 * n:
                dk = K.dkm(cc.loc[s, 'lat'], cc.loc[s, 'lon'], e.lat_best, e.lon_best) if s in cc.index else np.nan
                tr.append(('far' if s in FAR else 'near', s, x, dk))
        tr.sort(key=lambda z: (z[0] == 'far', z[3]))
        em = json.loads(e.emerg_json) if isinstance(e.get('emerg_json', None), str) else {}
        fig, axs = plt.subplots(1, 2, figsize=(13, 0.36 * len(tr) + 1.6), gridspec_kw=dict(width_ratios=[1, 1.6]))
        for j, (g, s, x, dk) in enumerate(tr):
            w = x[1500:2400]; pk = np.nanmax(np.abs(w[300:480]))
            axs[0].plot(np.arange(900) - 300, w / (2.2 * np.nanmax(np.abs(w))) - j, color='k' if g == 'near' else 'tab:red', lw=0.7)
            axs[0].text(605, -j, f'{s} {dk:.0f} km  {pk:.0f} nm  em {em.get(s, np.nan):.1f}', fontsize=6.5)
            axs[1].plot(np.arange(n) - 1800, x / (2.2 * np.nanmax(np.abs(x))) - j, color='k' if g == 'near' else 'tab:red', lw=0.5)
        for ax in axs:
            ax.axvspan(0, 180, color='0.9', zorder=0); ax.set_yticks([])
        axs[0].set_xlim(-300, 900); axs[0].set_xlabel('s after t0 (own-peak scaled)'); axs[1].set_xlabel('s after t0 (±30 min)')
        fig.suptitle(f'#{k} {e.time_pt}  MF {e.zmad:.1f}  n_tpl {e.n_templates}  z {e.z_best:.1f}  ratio_far {e.ratio_far:.2f}  loso {e.loso_frac:.2f}  '
                     f'emerg med {e.emerg_med_near:.1f} (n≥2: {e.n_near_emerg_ge2}/{e.n_near})  max/2nd {e.emerg_max_over_second:.1f} {e.emerg_lead_sta}  '
                     f'Spearman {e.t3_spearman:.2f}  storm {e.storm_pct_max_pm1d:.2f}  {e.lfe_status}  rule: {e.label_rule}', fontsize=7.5)
        fig.tight_layout(); fig.savefig(f'{K.FIG}/{pre}_{k:02d}.png', dpi=80); plt.close(fig)
    pd.set_option('display.width', 250)
    print(t[['time_pt', 'zmad', 'z_best', 'ratio_far', 'loso_frac', 'emerg_med_near', 'n_near_emerg_ge2', 'n_near', 'emerg_max_over_second',
             'emerg_lead_sta', 't3_spearman', 't3_pass', 'storm_pct_max_pm1d', 'lfe_status', 'label_rule']].round(2).to_string())


if __name__ == '__main__':
    m = sys.argv[1]
    if m == 'emerge':
        emerge(int(sys.argv[2]))
    else:
        {'step1': step1, 'ladder': ladder, 'figs': figs}[m]()
