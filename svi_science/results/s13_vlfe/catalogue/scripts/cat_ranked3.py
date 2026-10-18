#!/usr/bin/env python
"""Amendment Y8c: score fix of Y8b (saturated, separable logistic) + RESULTS 20 diagnostics.
  score  -> data/ranked/r3_all_<Y>.parquet (every real + control row with features and scores), data/ranked/score3_model.json,
            tables/catalogue_ranked3_<Y>.csv, catalogue_ranked3_2010_2026.csv, ranked3_curves_2010_2026.csv,
            ranked3_calibration.csv, ranked3_sanity.json, figures/ranked3_curves_pooled.png
  diag   -> tables/r20_tierA_between.csv (24 Y8b tierA between rows, all columns, conditions, EQ beyond mask, storm pct),
            tables/r20_condition_matrix.csv, tables/r20_topdays.csv, figures/r20_wave_<k>.png (top 5)"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K
import cat_ranked2 as R2

RD = R2.RD
YEARS = R2.YEARS
RUNS = R2.RUNS
FEATS = R2.FEATS
L2GRID = [1e-3, 1e-2, 1e-1, 1, 10, 100, 1000]          # x n (per-row penalty scale in R2.fit)


def known_lists():
    V = pd.read_csv(f'{K.TB}/e3_events.csv')
    known = np.sort(np.r_[V[V.role == 'vetted'].epoch.values, pd.read_csv(f'{K.TB}/pilot2013_ets_events_v2.csv').epoch.values])
    return known, np.sort(V[V.role == 'near-miss'].epoch.values)


def near(ep, ref, w=30):
    if not len(ref):
        return np.zeros(len(ep), bool)
    k = np.searchsorted(ref, ep)
    d = np.minimum(np.abs(ep - ref[np.clip(k - 1, 0, len(ref) - 1)]), np.abs(ep - ref[np.clip(k, 0, len(ref) - 1)]))
    return d <= w


def balanced_logloss(yv, t):
    w = np.where(yv == 1, 0.5 / yv.mean(), 0.5 / (1 - yv.mean()))
    return float((w * (np.logaddexp(0, t) - yv * t)).sum() / len(yv))


def cv_fit(X, yv, tag):
    from sklearn.model_selection import StratifiedKFold
    from sklearn.metrics import roc_auc_score
    rows = []
    oof = {}
    for l2 in L2GRID:
        lam = l2 * len(yv); t_oof = np.zeros(len(yv))
        for tr, te in StratifiedKFold(5, shuffle=True, random_state=1).split(X, yv):
            m = R2.fit(X[tr], yv[tr], l2=lam)
            t_oof[te] = m['intercept'] + ((X[te] - np.array(m['mean'])) / np.array(m['scale'])) @ np.array(m['coef'])
        rows.append(dict(neg=tag, l2_per_row=l2, cv_balanced_logloss=balanced_logloss(yv, t_oof), cv_auc=float(roc_auc_score(yv, t_oof))))
        oof[l2] = t_oof
    C = pd.DataFrame(rows)
    best = float(C.sort_values('cv_balanced_logloss').l2_per_row.iloc[0])
    M = R2.fit(X, yv, l2=best * len(yv)); M.update(l2_per_row=best, features=FEATS)
    return M, C, oof[best]


def logit(A, M):
    Z = (A[FEATS].values.astype(float) - np.array(M['mean'])) / np.array(M['scale'])
    return M['intercept'] + Z @ np.array(M['coef'])


def run_score():
    all_ = {}
    for y in YEARS:
        f = f'{RD}/r3_all_{y}.parquet'
        A = R2.feats(R2.assemble(y)) if not os.path.isfile(f) else pd.read_parquet(f)
        all_[y] = A; print('assembled', y, len(A), flush=True)
    known, nearmiss = known_lists()
    T = pd.concat([all_[y] for y in (2011, 2012, 2013)], ignore_index=True)
    real = T[T.ctrl == 'real']; P = real[near(real.epoch.values, known)]
    N = T[T.ctrl == 'noise_template']; S = T[T.ctrl == 'scrambled']
    X = np.r_[P[FEATS].values, N[FEATS].values].astype(float); yv = np.r_[np.ones(len(P)), np.zeros(len(N))]
    M, CV, oof = cv_fit(X, yv, 'noise_template')
    # variant: negatives = noise-template + scrambled rows, equal total weight (agreement on the real-bank footing)
    Ns = pd.concat([N, S.sample(n=min(len(S), len(N)), random_state=1)])
    Xv = np.r_[P[FEATS].values, Ns[FEATS].values].astype(float); yvv = np.r_[np.ones(len(P)), np.zeros(len(Ns))]
    Mv, CVv, _ = cv_fit(Xv, yvv, 'noise_template+scrambled')
    M.update(n_pos=int(len(P)), n_neg=int(len(N)), variant=dict(intercept=Mv['intercept'], coef=Mv['coef'], l2_per_row=Mv['l2_per_row'],
                                                                   mean=Mv['mean'], scale=Mv['scale']))
    json.dump(M, open(f'{RD}/score3_model.json', 'w'), indent=1)
    pd.concat([CV, CVv]).to_csv(f'{K.TB}/ranked3_cv_l2.csv', index=False)
    # calibration 1: out-of-fold probability (class-balanced) by decile
    p_oof = 1 / (1 + np.exp(-oof))
    bins = np.unique(np.quantile(p_oof, np.linspace(0, 1, 11)))
    cal = pd.DataFrame(dict(p=p_oof, y=yv)).assign(bin=lambda d: pd.cut(d.p, bins, include_lowest=True))
    cal1 = cal.groupby('bin', observed=True).agg(n=('y', 'size'), n_pos=('y', 'sum'), p_mean=('p', 'mean')).reset_index()
    cal1['bin'] = cal1.bin.astype(str); cal1['kind'] = 'oof_prob_decile'
    for y, A in all_.items():
        A['score_logit'] = logit(A, M); A['score_p'] = 1 / (1 + np.exp(-A.score_logit))
        A['score_v_logit'] = logit(A, Mv)
        A.to_parquet(f'{RD}/r3_all_{y}.parquet')
    ALL = pd.concat(all_.values(), ignore_index=True)
    # between-ETS live time per year (cat_ychecks.between_years) and catalogue counts
    import cat_ychecks as YC
    TB = {y: YC.between_years(y, RUNS[y]) for y in YEARS}
    nn = {y: int(all_[y][all_[y].ctrl == 'noise_template'].template.nunique()) for y in YEARS}
    ntot = sum(nn[y] * TB[y] for y in YEARS); bank_yr = sum(TB.values())
    # calibration 2: noise-template FAR at logit thresholds, even templates predict odd templates
    nb = ALL[(ALL.ctrl == 'noise_template') & (ALL.period == 'between')].copy()
    nb['tn'] = nb.template.str[1:].astype(int)
    grid = np.quantile(ALL[ALL.ctrl == 'real'].score_logit, [0.5, 0.9, 0.99, 0.995, 0.999, 0.9995])
    cal2 = []
    for g in grid:
        for par, lab in ((0, 'even'), (1, 'odd')):
            x = nb[nb.tn % 2 == par]
            ty = sum(TB[y] * all_[y][(all_[y].ctrl == 'noise_template')].template.str[1:].astype(int).pipe(lambda s: s[s % 2 == par]).nunique() for y in YEARS)
            cal2.append(dict(kind='noise_FAR_split', logit_ge=g, half=lab, n=int((x.score_logit >= g).sum()), template_years=ty,
                             rate_per_template_yr=float((x.score_logit >= g).sum() / ty)))
    cal2 = pd.DataFrame(cal2)
    pd.concat([cal1, cal2]).to_csv(f'{K.TB}/ranked3_calibration.csv', index=False)
    # curves on several footings (between ETS)
    real = ALL[ALL.ctrl == 'real']
    q = np.quantile(real.score_logit, np.r_[np.linspace(0, 0.99, 34), [0.993, 0.995, 0.997, 0.998, 0.999, 0.9995]])
    cur = []
    for per in ('ETS', 'between'):
        r = real[real.period == per]; n_ = ALL[(ALL.ctrl == 'noise_template') & (ALL.period == per)]
        s_ = ALL[(ALL.ctrl == 'scrambled') & (ALL.period == per)]
        for g in q:
            cur.append(dict(period=per, logit_ge=g, real=int((r.score_logit >= g).sum()), real_tierA=int(((r.score_logit >= g) & r.tierA).sum()),
                            noise_tpl_per_catalogue=float((n_.score_logit >= g).sum() / np.mean(list(nn.values()))),
                            noise_tpl_x11=float(11 * (n_.score_logit >= g).sum() / np.mean(list(nn.values()))),
                            scrambled_A1_per_draw=float(((s_.source == 'A1') & (s_.score_logit >= g)).sum() / 2),
                            scrambled_A2_per_draw=float(((s_.source == 'A2') & (s_.score_logit >= g)).sum() / 2)))
    C = pd.DataFrame(cur); C.to_csv(f'{K.TB}/ranked3_curves_2010_2026.csv', index=False)
    # per-row FAR columns (between-ETS live time)
    nsort = np.sort(nb.score_logit.values)
    sb = ALL[(ALL.ctrl == 'scrambled') & (ALL.period == 'between')]; ssort = np.sort(sb.score_logit.values)
    keep = [c for c in pd.read_csv(f'{K.TB}/catalogue_ranked2_2010.csv', nrows=1).columns if c != 'score'] + ['score_logit', 'score_p', 'score_v_logit']
    combo = []
    for y in YEARS:
        R = all_[y][all_[y].ctrl == 'real'].assign(run=RUNS[y])
        R['far_noise_tpl_per_tpl_yr'] = (len(nsort) - np.searchsorted(nsort, R.score_logit.values)) / ntot
        R['far_noise_tpl_bank11_yr'] = 11 * R.far_noise_tpl_per_tpl_yr
        R['far_scrambled_per_bank_yr'] = (len(ssort) - np.searchsorted(ssort, R.score_logit.values)) / (4 * bank_yr)
        R = R[[c for c in keep + ['far_noise_tpl_per_tpl_yr', 'far_noise_tpl_bank11_yr', 'far_scrambled_per_bank_yr'] if c in R]]
        R.sort_values('epoch').to_csv(f'{K.TB}/catalogue_ranked3_{y}.csv', index=False); combo.append(R)
    Rall = pd.concat(combo).sort_values('epoch'); Rall.to_csv(f'{K.TB}/catalogue_ranked3_2010_2026.csv', index=False)
    Rall['rank_pct'] = Rall.score_logit.rank(pct=True)
    kn = Rall[near(Rall.epoch.values, known)]; nm = Rall[near(Rall.epoch.values, nearmiss)]
    top = Rall[Rall.period == 'between'].sort_values('score_logit', ascending=False).head(10)
    top.to_csv(f'{K.TB}/ranked3_top10_between.csv', index=False)
    p_sat = float((Rall.score_p >= 0.999).mean())
    S = dict(l2_per_row=M['l2_per_row'], intercept=M['intercept'], coef=dict(zip(FEATS, np.round(M['coef'], 3).tolist())),
             variant_coef=dict(zip(FEATS, np.round(Mv['coef'], 3).tolist())), variant_l2=Mv['l2_per_row'],
             cv=CV.round(4).to_dict('records'), cv_variant=CVv.round(4).to_dict('records'),
             share_real_rows_p_ge_0p999=p_sat, n_distinct_logit_top100=int(Rall.score_logit.nlargest(100).round(3).nunique()),
             known_n=int(len(kn)), known_rank_pct_median=float(kn.rank_pct.median()), known_rank_pct_min=float(kn.rank_pct.min()),
             nearmiss_rank_pct_median=float(nm.rank_pct.median()) if len(nm) else None,
             between_live_yr=float(bank_yr), noise_template_years=float(ntot))
    json.dump(S, open(f'{K.TB}/ranked3_sanity.json', 'w'), indent=1, default=float)
    figs(C)
    pd.set_option('display.width', 250)
    print(json.dumps({k: v for k, v in S.items() if k not in ('cv', 'cv_variant')}, indent=1, default=float))
    print(pd.concat([CV, CVv]).round(4).to_string(index=False))
    print(C[C.period == 'between'].round(2).tail(12).to_string(index=False))
    print(cal1.round(3).to_string(index=False)); print(cal2.round(4).to_string(index=False))


def figs(C):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    for ax, per in zip(axs, ('ETS', 'between')):
        c = C[C.period == per]
        ax.semilogy(c.logit_ge, c.real.clip(lower=0.1), 'r-', label='real (bank of 11)')
        ax.semilogy(c.logit_ge, c.real_tierA.clip(lower=0.1), 'r:', label='real, Q2 tierA')
        ax.semilogy(c.logit_ge, c.noise_tpl_per_catalogue.clip(lower=0.01), 'k-', label='noise template, per catalogue')
        ax.semilogy(c.logit_ge, c.noise_tpl_x11.clip(lower=0.01), 'k--', label='noise template x 11')
        ax.semilogy(c.logit_ge, c.scrambled_A1_per_draw.clip(lower=0.01), color='tab:purple', ls='--', label='scrambled A1 (day shift), per draw')
        ax.semilogy(c.logit_ge, c.scrambled_A2_per_draw.clip(lower=0.01), color='tab:blue', ls='--', label='scrambled A2 (polarity), per draw')
        ax.set_xlabel('score (logit) ≥'); ax.set_ylabel('count, 2010-2018'); ax.set_title(per, fontsize=9); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f'{K.FIG}/ranked3_curves_pooled.png', dpi=120); plt.close(fig)


# ------------------------------------------------------------------------------------------------ RESULTS 20 diagnostics
def run_diag():
    import cat_ychecks as YC
    R = pd.read_csv(f'{K.TB}/catalogue_ranked2_2010_2026.csv')
    ALL = pd.concat([pd.read_parquet(f'{RD}/r3_all_{y}.parquet') for y in YEARS], ignore_index=True)
    # battery' items recomputed on the Y8b columns (cat_common.apply_rules) for the real rows
    bat = []
    for y in YEARS:
        b = pd.read_parquet(f'{RD}/r2bat_{y}.parquet'); k = pd.read_parquet(f'{RD}/r2keys_{y}.parquet')
        b = k[k.ctrl == 'real'][['i', 'epoch']].merge(b.drop(columns=['epoch']), on='i')
        b = K.apply_rules(b); b['year'] = y; bat.append(b[['epoch', 'year', 'pass_a', 'pass_b', 'pass_c', 'pass_e', 'pass_e2', 'pass_f',
                                                              'battery_e2', 'n_far_excited', 'in_head2', 'in_bw5', 'in_reg', 'in_reg3', 'in_m45']])
    B = pd.concat(bat).drop_duplicates('epoch')
    t = R[(R.period == 'between') & R.tierA]
    t = t.merge(B.drop(columns=['year'] + [c for c in B.columns if c in t.columns and c != 'epoch']), on='epoch', how='left')
    sc = ALL[ALL.ctrl == 'real'][['epoch', 'score_logit']].drop_duplicates('epoch')
    t = t.merge(sc, on='epoch', how='left')
    # storm percentile, EQs beyond the mask
    st = pd.concat([pd.read_csv(f'{K.TB}/{RUNS[y]}{y}_R_stormdays_v2.csv') for y in YEARS]).drop_duplicates('day')
    t['day'] = pd.to_datetime(t.epoch, unit='s').dt.date.astype(str)
    t = t.merge(st[['day', 'pct', 'pgc_5_10s_rms_nm']].rename(columns={'pct': 'storm_pct'}), on='day', how='left')
    EQ = pd.concat([K.eq_table(y) for y in YEARS + [2019]]).drop_duplicates('id').sort_values('t')
    eqc = []
    for r in t.itertuples():
        e = EQ[(EQ.t <= r.epoch) & (EQ.t >= r.epoch - 6 * 3600)]
        big = e.sort_values('a20_pred2_nm', ascending=False).head(1)
        e6 = EQ[(EQ.t <= r.epoch) & (EQ.t >= r.epoch - 24 * 3600) & (EQ.mag >= 6.5)]
        eqc.append(dict(eq6h_n=len(e), eq6h_max_pred_nm=float(big.a20_pred2_nm.iloc[0]) if len(big) else np.nan,
                        eq6h_desc=(f"M{big.mag.iloc[0]:.1f} {big.dist_km.iloc[0]:.0f} km {(r.epoch - big.t.iloc[0]) / 60:.0f} min before") if len(big) else '',
                        eq24h_M65=';'.join(f"M{x.mag:.1f} {x.dist_km:.0f}km -{(r.epoch - x.t) / 3600:.1f}h" for x in e6.itertuples())))
    t = pd.concat([t.reset_index(drop=True), pd.DataFrame(eqc)], axis=1).sort_values('score_logit', ascending=False)
    t['y8b_score_ge_0p9'] = t.score >= 0.9
    t.to_csv(f'{K.TB}/r20_tierA_between.csv', index=False)
    t_all = t
    def cnd(t): return dict(n=len(t), agreement=int((t.n_templates >= 2).sum()), mf_ge_5p5=int((t.zmad >= 5.5).sum()),
                pass_a_far=int(t.pass_a.sum()), pass_b_not_travelling=int(t.pass_b.sum()), pass_c_located=int(t.pass_c.sum()),
                pass_e_strict=int(t.pass_e.sum()), pass_e2_prime=int(t.pass_e2.sum()), pass_f_loso=int(t.pass_f.sum()),
                battery_prime_all=int(t.battery_e2.sum()), strong_tier=int(t.strong_tier.sum()), storm_day=int(t.storm_day.sum()),
                cluster_day=int(t.cluster_day.sum()), grid_edge=int(t.grid_edge.sum()), lfe_quiet=int((t.lfe_status == 'quiet').sum()),
                storm_pct_ge_0p7=int((t.storm_pct >= 0.7).sum()), eq24h_M65_any=int((t.eq24h_M65 != '').sum()))
    cond = pd.DataFrame(dict(all_tierA_between=cnd(t_all), y8b_score_ge_0p9=cnd(t_all[t_all.y8b_score_ge_0p9])))
    cond.to_csv(f'{K.TB}/r20_condition_matrix.csv')
    t = t_all[t_all.y8b_score_ge_0p9]
    # top days
    hi = ALL[(ALL.period == 'between') & (ALL.ctrl == 'real')]
    thr = float(hi.score_logit.quantile(0.98))
    days = []
    for d in ['2014-01-23', '2013-01-18', '2014-09-14', '2011-10-04']:
        y = int(d[:4]); A = ALL[(ALL.year == y)]
        A = A[pd.to_datetime(A.epoch, unit='s').dt.date.astype(str) == d]
        p0, p1 = pd.Timestamp(d).value / 1e9, pd.Timestamp(d).value / 1e9 + 86400
        e = EQ[(EQ.t >= p0 - 6 * 3600) & (EQ.t < p1)]
        r8 = R[(R.score >= 0.9) & (pd.to_datetime(R.epoch, unit='s').dt.date.astype(str) == d)]
        days.append(dict(day_utc=d, real_rows=int((A.ctrl == 'real').sum()), real_score2_ge_0p9=len(r8),
                         real_logit_ge_p98=int(((A.ctrl == 'real') & (A.score_logit >= thr)).sum()),
                         noise_tpl_rows=int((A.ctrl == 'noise_template').sum()), noise_tpl_agree=int(((A.ctrl == 'noise_template') & (A.agree_bin == 1)).sum()),
                         noise_tpl_logit_ge_p98=int(((A.ctrl == 'noise_template') & (A.score_logit >= thr)).sum()),
                         scr_rows_per_draw=float((A.ctrl == 'scrambled').sum() / 4), scr_logit_ge_p98_per_draw=float(((A.ctrl == 'scrambled') & (A.score_logit >= thr)).sum() / 4),
                         storm_pct=float(st[st.day == d].pct.iloc[0]) if (st.day == d).any() else np.nan,
                         pgc_5_10s_rms_nm=float(st[st.day == d].pgc_5_10s_rms_nm.iloc[0]) if (st.day == d).any() else np.nan,
                         n_eq_ge_m6=int((e.mag >= 6).sum()), eqs_ge_m6=';'.join(f"{x.time_pt} M{x.mag:.1f} {x.dist_km:.0f}km" for x in e[e.mag >= 6].itertuples()),
                         real_hours_pt=','.join(sorted(set(pd.to_datetime(r8.epoch, unit='s').dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.strftime('%H')))),
                         real_lat_range=f"{r8.lat_best.min():.2f}-{r8.lat_best.max():.2f}" if len(r8) else ''))
    D = pd.DataFrame(days); D.to_csv(f'{K.TB}/r20_topdays.csv', index=False)
    pd.set_option('display.width', 300); pd.set_option('display.max_columns', 80)
    print(cond.to_string()); print(D.to_string(index=False))
    cols = ['time_pt', 'score', 'score_logit', 'zmad', 'n_templates', 'z_best', 'ratio_far', 'n_far', 'vr_thrust', 'loso_frac', 'far_over_near2', 'pass_a',
            'pass_c', 'pass_e2', 'pass_f', 'battery_e2', 'storm_pct', 'storm_day', 'cluster_day', 'lat_best', 'lon_best', 'lfe_status', 'eq6h_desc', 'eq24h_M65']
    print(t[cols].round(2).to_string(index=False))
    waves(t.head(5))


def waves(t, pre='r20_wave'):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import scan_refs as SR
    cc = SR.all_coords()
    NEAR, REG, FAR = R2.NEAR, ['UW.GNW', 'UW.DOSE', 'UW.SP2', 'TA.D03D', 'TA.B05D', 'UW.WISH', 'US.NLWA'], R2.FAR
    out = []
    for k, e in enumerate(t.itertuples()):
        y = int(e.year); SD = f'{K.S13}/data/scan' if y == 2011 else f'{K.S13}/data/scan_{y}'
        refbp = f'{SD}/bp2_ref_v2' if y in (2011, 2012) else f'{SD}/bp2_ref'
        ts = int(round(e.epoch - pd.Timestamp(f'{y}-01-01').value / 1e9)); a, n = ts - 1800, 3600
        tr = []
        for s in NEAR + REG + FAR:
            f = f'{SD}/bp_{y}g/{s}.npy' if s in NEAR else f'{refbp}/{s}.npy'
            if os.path.isfile(f):
                X = np.load(f, mmap_mode='r')
                if a < 0 or a + n > X.shape[1]:
                    continue
                x = np.asarray(X[2, a:a + n], float) * 1e9
                if np.isfinite(x).sum() > 0.8 * n:
                    dk = K.dkm(cc.loc[s, 'lat'], cc.loc[s, 'lon'], e.lat_best, e.lon_best) if s in cc.index else np.nan
                    tr.append(('far' if s in FAR else 'near', s, x, dk))
        tr.sort(key=lambda z: (z[0] == 'far', z[3]))
        fig, axs = plt.subplots(1, 2, figsize=(14, 0.42 * len(tr) + 1.8), gridspec_kw=dict(width_ratios=[1, 2]))
        for j, (g, s, x, dk) in enumerate(tr):
            w = x[1800 - 300:1800 + 600]; pk = np.nanmax(np.abs(w[300:480])); pk0 = np.nanmax(np.abs(w[:300]))
            axs[0].plot(np.arange(900) - 300, w / (2.2 * np.nanmax(np.abs(w))) - j, color='k' if g == 'near' else 'tab:red', lw=0.7)
            axs[0].text(605, -j, f'{s} {dk:.0f} km  win {pk:.0f} / pre {pk0:.0f} nm', fontsize=6.5)
            axs[1].plot(np.arange(n) - 1800, x / (2.2 * np.nanmax(np.abs(x))) - j, color='k' if g == 'near' else 'tab:red', lw=0.5)
            out.append(dict(k=k, time_pt=e.time_pt, sta=s, dist_km=dk, win_peak_nm=pk, pre300_peak_nm=pk0, hour_rms_nm=float(np.nanstd(x))))
        axs[0].axvspan(0, 180, color='0.9', zorder=0); axs[0].set_xlim(-300, 900); axs[0].set_yticks([]); axs[0].set_xlabel('s after t0 (own-peak scaled)')
        axs[1].axvspan(0, 180, color='0.9', zorder=0); axs[1].set_yticks([]); axs[1].set_xlabel('s after t0 (±30 min)')
        fig.suptitle(f'{e.time_pt}  MF {e.zmad:.1f}  n_tpl {e.n_templates}  scan z {e.z_best:.1f}  vr {e.vr_thrust:.2f}  far/near2 {e.far_over_near2:.2f}  '
                     f'storm pct {e.storm_pct:.2f}  {e.lfe_status}  {e.eq6h_desc}\nZ, 0.02–0.05 Hz displacement; near (black) by distance, far refs (red)', fontsize=8)
        fig.tight_layout(); fig.savefig(f'{K.FIG}/{pre}_{k}.png', dpi=105); plt.close(fig)
    pd.DataFrame(out).to_csv(f'{K.TB}/{pre}_amps.csv', index=False)


if __name__ == '__main__':
    {'score': run_score, 'diag': run_diag}[sys.argv[1]]()
