#!/usr/bin/env python
"""Amendment Y8b: ranked catalogue with full check columns and a monotone score.
  battery Y   (scan-v2 env of year Y) -> data/ranked/r2bat_<Y>.parquet: expensive battery columns for every agreement row
                (real; noise templates; scrambled draws A1_0, A1_1, A2_0, A2_1)
  score       -> data/ranked/score2_model.json, tables/catalogue_ranked2_<Y>.csv, catalogue_ranked2_2010_2026.csv,
                 ranked2_curves_<Y>.csv, ranked2_curves_2010_2026.csv, figures/ranked2_curves_<Y>.png, tables/ranked2_sanity.json"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

RD = f'{K.DATA}/ranked'
YEARS = list(range(2010, 2019)) + [int(x) for x in os.environ.get('LADDER_EXTRA', '').split(',') if x]   # Y11: 2019+
RUNS = {2010: 'prod', 2011: 'prod', 2012: 'prod', 2013: 'prodq', 2014: 'prod', 2015: 'prod', 2016: 'prodq', 2017: 'prodq',
        2018: 'prodq', **{y: 'prodf' for y in range(2019, 2027)}}
RUNS.update({int(a): b for a, b in (x.split('=') for x in os.environ.get('LADDER_RUNS', '').split(',') if x)})   # e.g. 2021=prodfy (Y12)
SCR_DRAWS = [('A1', 0), ('A1', 1), ('A2', 0), ('A2', 1)]
NEAR = ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.PFB', 'CN.YOUB', 'UW.LRIV', 'CN.LZB']
FAR = ['IU.COR', 'US.NEW', 'US.HAWA', 'US.MSO', 'US.WVOR']
EXP = ['ratio_far', 'n_far', 'ratio_GNW', 'vr_thrust', 'z_loso', 'travelling', 'n_Z_agree', 'share_excess', 'm0_cols']


def load(y):
    R = pd.read_csv(f'{K.TB}/catalogue_ranked_{y}.csv'); R['ctrl'] = 'real'; R['source'] = 'real'; R['draw'] = -1
    C = pd.read_parquet(f'{K.DATA}/ranked/controls_{y}.parquet')
    C = C[(C.ctrl == 'noise_template') | C.apply(lambda r: (r.source, int(r.draw)) in SCR_DRAWS, axis=1)]
    A = pd.concat([R, C], ignore_index=True)
    A['year'] = y
    A['agree_bin'] = (A.n_templates >= 2).astype(int)
    return A


def run_battery(y):
    import cat_e2
    A = load(y)
    B = A[A.agree_bin == 1].copy()
    B['i'] = np.arange(len(B)); B['template'] = B.template.astype(str)
    B[['i', 'epoch', 'year', 'template']].to_parquet(f'{RD}/r2rows_{y}.parquet')
    B[['i', 'epoch', 'ctrl', 'source', 'draw', 'template']].to_parquet(f'{RD}/r2keys_{y}.parquet')
    cat_e2.run_battery(y, rows_file=f'{RD}/r2rows_{y}.parquet', out_file=f'{RD}/r2bat_{y}.parquet', tag='rank2', gate_all=True)


def amp_fn(y):
    SD = f'{K.S13}/data/scan' if y == 2011 else f'{K.S13}/data/scan_{y}'
    refbp = (f'{SD}/bp2_ref_v2' if y in (2011, 2012) else f'{SD}/bp2_ref')
    cache = {}

    def amp(s, ep):
        if s not in cache:
            f = f'{SD}/bp_{y}g/{s}.npy' if s in NEAR else f'{refbp}/{s}.npy'
            cache[s] = np.load(f, mmap_mode='r') if os.path.isfile(f) else None
        X = cache[s]
        if X is None:
            return np.nan
        ts = int(round(ep - pd.Timestamp(f'{y}-01-01').value / 1e9))
        if ts < 0 or ts + 180 > X.shape[1]:
            return np.nan
        x = np.asarray(X[2, ts:ts + 180], float) * 1e9
        return float(np.nanmax(np.abs(x))) if np.isfinite(x).sum() > 150 else np.nan
    return amp


def assemble(y):
    import scan_refs as SR
    cc = SR.all_coords()
    A = load(y)
    A['i'] = -1; ia = np.where(A.agree_bin == 1)[0]; A.loc[A.index[ia], 'i'] = np.arange(len(ia))
    keys = pd.read_parquet(f'{RD}/r2keys_{y}.parquet'); bat = pd.read_parquet(f'{RD}/r2bat_{y}.parquet')
    chk = A[A.i >= 0].merge(keys[['i', 'epoch']], on='i', suffixes=('', '_k'))
    assert len(chk) == len(ia) and np.allclose(chk.epoch, chk.epoch_k), f'row keys changed for {y}'
    cols = [c for c in EXP + ['z_best', 'lat_best', 'lon_best', 'grid_edge_node', 'in_footprint_poly', 'dist_family_km', 'nearest_family'] if c in bat]
    A = A.merge(bat[['i'] + cols].drop_duplicates('i').rename(columns={c: c + '_r2' for c in cols}), on='i', how='left')
    for c in cols:
        if c == 'grid_edge_node':
            continue
        A[c] = np.where(A.agree_bin == 1, A[c + '_r2'], A[c] if c in A else np.nan)
    ge_old = pd.to_numeric(A['grid_edge'].astype(object).map({True: 1, False: 0, 'True': 1, 'False': 0}), errors='coerce') if 'grid_edge' in A else np.nan
    ge_new = pd.to_numeric(A['grid_edge_node_r2'].astype(object).map({True: 1, False: 0}), errors='coerce')
    A['grid_edge'] = np.where(A.agree_bin == 1, ge_new, ge_old) > 0
    A = A[[c for c in A.columns if not c.endswith('_r2')]]
    # pass flags recomputed from the (now complete) columns
    A['pass_far'] = (A.n_far >= 2) & (A.ratio_far <= 0.551)
    A['not_travelling'] = A.travelling.notna() & ~A.travelling.fillna(False).astype(bool)
    A['loso_frac'] = A.z_loso / A.z_best
    A['pass_c'] = (A.z_best >= 5) & ~A.grid_edge & (A.in_footprint_poly.fillna(False).astype(bool) | (A.dist_family_km <= 20))
    A['n_checks_missing'] = A[['ratio_far', 'vr_thrust', 'z_loso', 'travelling']].isna().sum(1)
    # Q2 / Tier A' for agreement rows
    amp = amp_fn(y)
    run = RUNS[y]
    st = pd.read_csv(f'{K.TB}/{run}{y}_R_stormdays_v2.csv'); sd = dict(zip(st.day, st.storm_day))
    A['day'] = pd.to_datetime(A.epoch, unit='s').dt.date.astype(str)
    A['storm_day'] = A.day.map(sd).fillna(False).astype(bool)
    fo = np.full(len(A), np.nan)
    for k in np.where(A.agree_bin == 1)[0]:
        r = A.iloc[k]
        if not np.isfinite(r.lat_best):
            continue
        live = [s for s in NEAR if s in cc.index and np.isfinite(amp(s, r.epoch))]
        n2 = sorted(live, key=lambda s: K.dkm(cc.loc[s, 'lat'], cc.loc[s, 'lon'], r.lat_best, r.lon_best))[:2]
        nm = max([amp(s, r.epoch) for s in n2] or [np.nan]); fm = np.nanmax([amp(s, r.epoch) for s in FAR] + [np.nan])
        fo[k] = fm / nm if np.isfinite(nm) and nm > 0 else np.nan
    A['far_over_near2'] = fo
    A['q2'] = (A.vr_thrust >= 0.12) & (A.far_over_near2 < 0.5) & np.where(A.lat_best > 48.1, A.ratio_GNW <= 0.551, True)
    # cluster day: >= 2 Q2 events (30-s declustered) on the same day in the same catalogue (real; one noise template; one draw)
    A['cat_key'] = np.where(A.ctrl == 'noise_template', A.template.astype(str), A.ctrl + A.source.astype(str) + A.draw.astype(str))
    q = A[A.q2 & (A.agree_bin == 1)].sort_values('epoch')
    q = q[q.groupby('cat_key').epoch.diff().fillna(1e9) > 30]
    nq = q.groupby(['cat_key', 'day']).size()
    nqd = nq.to_dict()
    A['cluster_day'] = np.array([nqd.get(t, 0) for t in zip(A.cat_key, A.day)]) >= 2
    A['tierA'] = A.q2 & (A.agree_bin == 1) & ~(A.storm_day & A.cluster_day)
    return A


FEATS = ['agree_bin', 'zm', 'lzb', 'pass_far_f', 'neg_ratio_far', 'vr_f', 'loso_f', 'pass_e2_f', 'pass_c_f', 'nottrav_f', 'notedge_f']


def feats(A):
    A = A.copy()
    A['zm'] = (A.zmad - 5.5).clip(0, 20)
    A['lzb'] = np.log10(A.z_best.clip(lower=0.5).fillna(0.5))
    A['pass_far_f'] = A.pass_far.astype(float)
    A['neg_ratio_far'] = -A.ratio_far.fillna(3.0).clip(0, 3)
    A['vr_f'] = A.vr_thrust.fillna(0.0).clip(-1, 1)
    A['loso_f'] = A.loso_frac.fillna(-1.0).clip(-1, 1.5)   # Y9: missing -> clip floor
    A['pass_e2_f'] = A.pass_e2.fillna(False).astype(float)
    A['pass_c_f'] = A.pass_c.astype(float)
    A['nottrav_f'] = A.not_travelling.astype(float)
    A['notedge_f'] = (~A.grid_edge).astype(float)
    return A


def fit(X, yv, l2=1.0):
    from scipy.optimize import minimize
    mu, sd = X.mean(0), X.std(0) + 1e-9
    Z = (X - mu) / sd
    w = np.where(yv == 1, 0.5 / yv.mean(), 0.5 / (1 - yv.mean()))

    def f(p):
        b, c = p[0], p[1:]
        t = b + Z @ c
        ll = w * (np.logaddexp(0, t) - yv * t)
        g = w * (1 / (1 + np.exp(-t)) - yv)
        return ll.sum() / len(yv) + 0.5 * l2 * (c @ c) / len(yv), np.r_[g.sum(), Z.T @ g] / len(yv) + np.r_[0, l2 * c / len(yv)]
    p0 = np.zeros(Z.shape[1] + 1)
    res = minimize(f, p0, jac=True, method='L-BFGS-B', bounds=[(None, None)] + [(0, None)] * Z.shape[1])
    return dict(mean=mu.tolist(), scale=sd.tolist(), intercept=float(res.x[0]), coef=res.x[1:].tolist(), ok=bool(res.success))


def apply(A, M):
    Z = (A[FEATS].values.astype(float) - np.array(M['mean'])) / np.array(M['scale'])
    return 1 / (1 + np.exp(-(M['intercept'] + Z @ np.array(M['coef']))))


def run_score():
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold
    all_ = {}
    for y in YEARS:
        all_[y] = feats(assemble(y)); print('assembled', y, len(all_[y]), flush=True)
    V = pd.read_csv(f'{K.TB}/e3_events.csv')
    nearmiss = np.sort(V[V.role == 'near-miss'].epoch.values)
    known = np.sort(np.r_[V[V.role == 'vetted'].epoch.values, pd.read_csv(f'{K.TB}/pilot2013_ets_events_v2.csv').epoch.values])
    T = pd.concat([all_[y] for y in (2011, 2012, 2013)], ignore_index=True)
    real = T[T.ctrl == 'real']
    k = np.searchsorted(known, real.epoch.values)
    d = np.minimum(np.abs(real.epoch.values - known[np.clip(k - 1, 0, len(known) - 1)]), np.abs(real.epoch.values - known[np.clip(k, 0, len(known) - 1)]))
    P = real[d <= 30]; N = T[T.ctrl == 'noise_template']
    X = np.r_[P[FEATS].values, N[FEATS].values].astype(float); yv = np.r_[np.ones(len(P)), np.zeros(len(N))]
    aucs = []
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=1).split(X, yv):
        m = fit(X[tr], yv[tr]); z = (X[te] - np.array(m['mean'])) / np.array(m['scale'])
        aucs.append(roc_auc_score(yv[te], m['intercept'] + z @ np.array(m['coef'])))
    M = fit(X, yv); M.update(features=FEATS, n_pos=int(len(P)), n_neg=int(len(N)), cv_auc=aucs,
                             note='monotone (coef >= 0) logistic; missing set to worst value; positives = known 2011-13 VLFEs')
    json.dump(M, open(f'{RD}/score2_model.json', 'w'), indent=1)
    grid = np.round(np.arange(0, 1.0001, 0.02), 2)
    combo, curves = [], []
    for y, A in all_.items():
        A['score'] = apply(A, M)
        R = A[A.ctrl == 'real']; nn = int(pd.read_parquet(f'{K.DATA}/{RUNS[y]}{y}/rows_scored_v2.parquet', columns=['kind', 'template']).query('kind == "noise"').template.nunique())
        keep = ['time_pt', 'period', 'episode', 'template', 'agree', 'n_templates', 'zmad', 'z_best', 'ratio_far', 'n_far', 'ratio_GNW',
                'vr_thrust', 'z_loso', 'loso_frac', 'travelling', 'pass_far', 'pass_e', 'pass_e2', 'in_head2', 'head2_maxmag', 'pass_c',
                'grid_edge', 'in_footprint_poly', 'lat_best', 'lon_best', 'nearest_family', 'dist_family_km', 'lfe_status', 'n_lfe_10m',
                'far_over_near2', 'storm_day', 'cluster_day', 'strong_tier', 'tierA', 'n_checks_missing', 'score', 'epoch', 'year', 'run']
        R = R.assign(run=RUNS[y])
        R[[c for c in keep if c in R]].sort_values('epoch').to_csv(f'{K.TB}/catalogue_ranked2_{y}.csv', index=False)
        combo.append(R[[c for c in keep if c in R]])
        for per in ('ETS', 'between'):
            r = R[R.period == per]; n_ = A[(A.ctrl == 'noise_template') & (A.period == per)]; s_ = A[(A.ctrl == 'scrambled') & (A.period == per)]
            for g in grid:
                curves.append(dict(year=y, period=per, score_ge=g, real=int((r.score >= g).sum()), real_tierA=int(((r.score >= g) & r.tierA).sum()),
                                   noise_template_per_catalogue=float((n_.score >= g).sum() / max(nn, 1)),
                                   noise_template_tierA_per_catalogue=float(((n_.score >= g) & n_.tierA).sum() / max(nn, 1)),
                                   scrambled_per_draw=float((s_.score >= g).sum() / len(SCR_DRAWS))))
    pd.concat(combo).sort_values('epoch').to_csv(f'{K.TB}/catalogue_ranked2_2010_2026.csv', index=False)
    C = pd.DataFrame(curves); C.to_csv(f'{K.TB}/ranked2_curves_2010_2026.csv', index=False)
    for y in YEARS:
        C[C.year == y].to_csv(f'{K.TB}/ranked2_curves_{y}.csv', index=False)
    figs(C)
    # sanity
    Rall = pd.concat(combo); Rall['rank_pct'] = Rall.score.rank(pct=True)
    k = np.searchsorted(known, Rall.epoch.values)
    d = np.minimum(np.abs(Rall.epoch.values - known[np.clip(k - 1, 0, len(known) - 1)]), np.abs(Rall.epoch.values - known[np.clip(k, 0, len(known) - 1)]))
    kn = Rall[d <= 30]
    k2 = np.searchsorted(nearmiss, Rall.epoch.values)
    d2 = np.minimum(np.abs(Rall.epoch.values - nearmiss[np.clip(k2 - 1, 0, len(nearmiss) - 1)]), np.abs(Rall.epoch.values - nearmiss[np.clip(k2, 0, len(nearmiss) - 1)]))
    nmr = Rall[d2 <= 30]
    pooled = C.groupby(['period', 'score_ge'])[['real', 'real_tierA', 'noise_template_per_catalogue', 'noise_template_tierA_per_catalogue',
                                                  'scrambled_per_draw']].sum().reset_index()
    S = dict(known_found=int(len(kn)), known_rank_pct_median=float(kn.rank_pct.median()), known_rank_pct_min=float(kn.rank_pct.min()),
             known_best_rank_per_event_pct=float(kn.assign(ev=np.round(kn.epoch / 60)).groupby('ev').rank_pct.max().median()),
             nearmiss_rank_pct_median=float(nmr.rank_pct.median()) if len(nmr) else None, n_real=int(len(Rall)),
             pooled_at_0p9={p: pooled[(pooled.period == p) & (pooled.score_ge == 0.9)].iloc[0].drop(['period', 'score_ge']).to_dict() for p in ('ETS', 'between')},
             pooled_at_0p5={p: pooled[(pooled.period == p) & (pooled.score_ge == 0.5)].iloc[0].drop(['period', 'score_ge']).to_dict() for p in ('ETS', 'between')})
    json.dump(S, open(f'{K.TB}/ranked2_sanity.json', 'w'), indent=1, default=float)
    top = Rall[Rall.period == 'between'].sort_values('score', ascending=False).head(10)
    top.to_csv(f'{K.TB}/ranked2_top10_between.csv', index=False)
    pd.set_option('display.width', 250)
    print(json.dumps({k: (np.round(v, 3).tolist() if k in ('coef', 'cv_auc') else v) for k, v in M.items() if k not in ('mean', 'scale')}, indent=1))
    print(json.dumps(S, indent=1, default=float))
    print(top[['time_pt', 'score', 'zmad', 'n_templates', 'z_best', 'ratio_far', 'vr_thrust', 'loso_frac', 'pass_c', 'grid_edge', 'tierA',
               'lfe_status', 'nearest_family']].round(3).to_string(index=False))


def figs(C):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for y in YEARS + ['pooled']:
        c0 = C if y == 'pooled' else C[C.year == y]
        if y == 'pooled':
            c0 = c0.groupby(['period', 'score_ge']).sum(numeric_only=True).reset_index()
        fig, axs = plt.subplots(1, 2, figsize=(11, 4))
        for ax, per in zip(axs, ('ETS', 'between')):
            c = c0[c0.period == per]
            ax.semilogy(c.score_ge, c.real.clip(lower=0.1), 'r-', label='real')
            ax.semilogy(c.score_ge, c.real_tierA.clip(lower=0.1), 'r:', label='real, Tier A′')
            ax.semilogy(c.score_ge, c.noise_template_per_catalogue.clip(lower=0.01), 'k-', label='noise template (per catalogue)')
            ax.semilogy(c.score_ge, c.noise_template_tierA_per_catalogue.clip(lower=0.01), 'k:', label='noise template, Tier A′')
            ax.semilogy(c.score_ge, c.scrambled_per_draw.clip(lower=0.01), color='tab:purple', ls='--', label='scrambled (per draw)')
            ax.set_xlabel('score ≥'); ax.set_ylabel('count'); ax.set_title(f'{y} {per}', fontsize=9); ax.legend(fontsize=7)
        fig.tight_layout(); fig.savefig(f'{K.FIG}/ranked2_curves_{y}.png', dpi=120); plt.close(fig)


if __name__ == '__main__':
    if sys.argv[1] == 'battery':
        run_battery(int(sys.argv[2]))
    elif sys.argv[1] == 'smoke':                       # tiny end-to-end test in a scratch dir
        RD = sys.argv[3]; y = int(sys.argv[2]); os.makedirs(RD, exist_ok=True)
        import cat_e2
        _load = load
        def load(y, _l=_load):
            A = _l(y); g = A[A.agree_bin == 1].groupby('ctrl').head(4).index
            return A[(A.agree_bin == 0) | A.index.isin(g)].reset_index(drop=True)
        run_battery(y); A = feats(assemble(y))
        print(A[A.agree_bin == 1][['ctrl', 'epoch', 'z_best', 'ratio_far', 'n_far', 'vr_thrust', 'loso_frac', 'grid_edge', 'pass_c', 'far_over_near2', 'q2', 'tierA', 'n_checks_missing']].to_string())
        print(A[FEATS].describe().T[['mean', 'min', 'max']])
    else:
        run_score()
