#!/usr/bin/env python
"""Amendment Y8: ranked catalogue product (rank, do not discard).
  fit                 -> data/ranked/score_model.json (logistic regression on 2011-2013; positives = known VLFEs)
  year Y [RUN]        -> tables/catalogue_ranked_<Y>.csv, data/ranked/controls_<Y>.parquet, tables/ranked_curves_<Y>.csv,
                         figures/ranked_curves_<Y>.png; appends to tables/catalogue_ranked_2010_2026.csv
Default runs per year (final 2011-bank / adopted bank): RUNS below; override with the RUN argument."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

RD = f'{K.DATA}/ranked'
os.makedirs(RD, exist_ok=True)
RUNS = {2010: 'prod', 2011: 'prod', 2012: 'prod', 2013: 'prodq', 2014: 'prod', 2015: 'prod', 2016: 'prodq', 2017: 'prodq',
        2018: 'prodq'}
VALID = [f'T{k}' for k in range(1, 6)] + [f'C{k}' for k in range(6)]
LEVEL = 5.5
FEATS = ['zmad', 'n_templates', 'log_zbest', 'ratio_far_f', 'vr_thrust_f', 'loso_frac_f', 'pass_e2_f', 'pass_c_f',
         'far_missing', 'travelling_f']
COLS = ['time_pt', 'period', 'episode', 'template', 'agree', 'n_templates', 'zmad', 'z_best', 'ratio_far', 'n_far', 'ratio_GNW',
        'vr_thrust', 'z_loso', 'loso_frac', 'in_head2', 'head2_maxmag', 'in_bw5', 'in_reg', 'in_reg3', 'pass_e', 'pass_e2',
        'pass_c', 'in_footprint_poly', 'grid_edge', 'travelling', 'lat_best', 'lon_best', 'dep_best', 'nearest_family',
        'dist_family_km', 'lfe_status', 'n_lfe_10m', 'battery_e2', 'strong_tier', 'score', 'epoch']


def episodes_of(y, ep):
    E = pd.read_csv(f'{K.ROOT}/svi_science/results/b_ets_record/episodes_b.csv')
    E = E[E.region == 'SVI']
    out = np.full(len(ep), '', dtype=object)
    t = pd.to_datetime(ep, unit='s')
    for r in E.itertuples():
        a = pd.Timestamp(r.start) - pd.Timedelta(days=2); b = pd.Timestamp(r.end) + pd.Timedelta(days=3)
        out[(t >= a) & (t < b)] = r.episode
    return out


def dedupe(x, keys):
    """one row per event (+-30 s) within each group of keys; representative = max zmad; n_templates, agree union."""
    out = []
    for _, g in x.sort_values('epoch').groupby(keys):
        g = g.assign(grp=(g.epoch.diff().fillna(1e9) > 30).cumsum())
        for _, h in g.groupby('grp'):
            r = h.loc[h.zmad.idxmax()].copy()
            tpls = set(h.template)
            for a in h.agree.astype(str):
                tpls |= {t for t in a.split(',') if t}
            r['n_templates'] = len(tpls); r['agree'] = ','.join(sorted(tpls - {r.template}))
            r['strong_tier'] = bool(h.tier.any()) if 'tier' in h else False
            out.append(r)
    return pd.DataFrame(out)


def feats(D):
    D = D.copy()
    D['log_zbest'] = np.log10(D.z_best.clip(lower=0.1).fillna(0.1))
    D['far_missing'] = D.ratio_far.isna().astype(float)
    D['ratio_far_f'] = D.ratio_far.fillna(2.0).clip(0, 3)
    D['vr_thrust_f'] = D.vr_thrust.fillna(0.0)
    D['loso_frac'] = D.z_loso / D.z_best
    D['loso_frac_f'] = D.loso_frac.fillna(0.0).clip(-1, 1.5)
    D['pass_e2_f'] = D.pass_e2.astype(float)
    D['pass_c_f'] = D.pass_c.astype(float)
    D['travelling_f'] = D.travelling.fillna(False).astype(float)
    return D


def load_year(y, run):
    D = pd.read_parquet(f'{K.DATA}/{run}{y}/rows_scored_v2.parquet')
    if 'n_lfe_10m' not in D:
        D['n_lfe_10m'] = np.nan
    real = D[(D.source == 'real') & D.template.isin(VALID + [t for t in D.template.unique() if str(t).startswith('R')]) & (D.zmad >= LEVEL)]
    noise = D[(D.kind == 'noise') & (D.zmad >= LEVEL)]
    scr = D[D.source.isin(['A1', 'A2']) & (D.zmad >= LEVEL)]
    R = dedupe(real, ['source'])
    N = dedupe(noise.assign(src=noise.template), ['template'])
    S = dedupe(scr, ['source', 'draw'])
    for X in (R, N, S):
        if len(X):
            X['episode'] = episodes_of(y, X.epoch.values)
    n_noise = int(D[D.kind == 'noise'].template.nunique())
    return feats(R), feats(N), feats(S), n_noise


def run_fit():
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_score
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    known = []
    V = pd.read_csv(f'{K.TB}/e3_events.csv'); known += list(V.epoch)
    P13 = pd.read_csv(f'{K.TB}/pilot2013_ets_events_v2.csv'); known += list(P13.epoch)
    known = np.sort(np.array(known, float))
    pos, neg = [], []
    for y in (2011, 2012, 2013):
        R, N, S, nn = load_year(y, RUNS[y])
        k = np.searchsorted(known, R.epoch.values)
        d = np.minimum(np.abs(R.epoch.values - known[np.clip(k - 1, 0, len(known) - 1)]), np.abs(R.epoch.values - known[np.clip(k, 0, len(known) - 1)]))
        pos.append(R[d <= 30]); neg.append(N)
    Pp, Nn = pd.concat(pos), pd.concat(neg)
    X = pd.concat([Pp[FEATS], Nn[FEATS]]).values.astype(float); yv = np.r_[np.ones(len(Pp)), np.zeros(len(Nn))]
    model = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, class_weight='balanced', max_iter=2000))
    auc = cross_val_score(model, X, yv, cv=5, scoring='roc_auc')
    model.fit(X, yv)
    sc, lr = model.named_steps['standardscaler'], model.named_steps['logisticregression']
    M = dict(features=FEATS, mean=sc.mean_.tolist(), scale=sc.scale_.tolist(), coef=lr.coef_[0].tolist(), intercept=float(lr.intercept_[0]),
             n_pos=int(len(Pp)), n_neg=int(len(Nn)), cv_auc_mean=float(auc.mean()), cv_auc_folds=auc.tolist(),
             note='positives = 2011-12 vetted/near-miss + 2013 pilot ETS tier events found at >= 5.5 MAD; negatives = noise-template detections 2011-2013')
    json.dump(M, open(f'{RD}/score_model.json', 'w'), indent=1)
    print(json.dumps({k: (np.round(v, 3).tolist() if isinstance(v, list) and k in ('coef', 'cv_auc_folds') else v) for k, v in M.items()
                      if k not in ('mean', 'scale')}, indent=1))


def score(D, M):
    X = D[M['features']].values.astype(float)
    z = (X - np.array(M['mean'])) / np.array(M['scale'])
    return 1 / (1 + np.exp(-(z @ np.array(M['coef']) + M['intercept'])))


def run_year(y, run=None):
    run = run or RUNS[y]
    M = json.load(open(f'{RD}/score_model.json'))
    R, N, S, nn = load_year(y, run)
    for X in (R, N, S):
        if len(X):
            X['score'] = score(X, M)
            X['time_pt'] = pd.to_datetime(X.epoch, unit='s').dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.strftime('%Y-%m-%d %H:%M:%S %Z')
    R['year'] = y; R['run'] = run
    cols = [c for c in COLS if c in R] + ['year', 'run']
    R.sort_values('epoch')[cols].to_csv(f'{K.TB}/catalogue_ranked_{y}.csv', index=False)
    pd.concat([N.assign(ctrl='noise_template'), S.assign(ctrl='scrambled')], ignore_index=True)[
        [c for c in cols if c in N] + ['ctrl', 'source', 'draw']].to_parquet(f'{RD}/controls_{y}.parquet')
    # combined file (replace this year's rows)
    fc = f'{K.TB}/catalogue_ranked_2010_2026.csv'
    old = pd.read_csv(fc) if os.path.isfile(fc) else pd.DataFrame()
    if len(old):
        old = old[old.year != y]
    pd.concat([old, R[cols]], ignore_index=True).sort_values('epoch').to_csv(fc, index=False)
    # curves
    grid = np.round(np.arange(0.0, 1.0001, 0.02), 2)
    rows = []
    ndraw = 16
    for per in ('ETS', 'between'):
        r = R[R.period == per]; n_ = N[N.period == per] if len(N) else N; s_ = S[S.period == per] if len(S) else S
        for g in grid:
            rows.append(dict(year=y, period=per, score_ge=g, real=int((r.score >= g).sum()),
                             noise_template_per_catalogue=float((n_.score >= g).sum() / max(nn, 1)) if len(n_) else 0.0,
                             scrambled_per_catalogue=float((s_.score >= g).sum() / ndraw) if len(s_) else 0.0))
    C = pd.DataFrame(rows); C.to_csv(f'{K.TB}/ranked_curves_{y}.csv', index=False)
    fcur = f'{K.TB}/ranked_curves_2010_2026.csv'
    oc = pd.read_csv(fcur) if os.path.isfile(fcur) else pd.DataFrame()
    if len(oc):
        oc = oc[oc.year != y]
    pd.concat([oc, C], ignore_index=True).to_csv(fcur, index=False)
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(11, 4), sharey=False)
    for ax, per in zip(axs, ('ETS', 'between')):
        c = C[C.period == per]
        ax.semilogy(c.score_ge, c.real.clip(lower=0.1), 'r-', label='real events')
        ax.semilogy(c.score_ge, c.noise_template_per_catalogue.clip(lower=0.01), 'k-', label='noise-template control (per catalogue)')
        ax.semilogy(c.score_ge, c.scrambled_per_catalogue.clip(lower=0.01), color='tab:purple', ls='--', label='scrambled control (per draw)')
        ax.set_xlabel('combined score ≥'); ax.set_ylabel('count'); ax.set_title(f'{y} {per}', fontsize=9); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f'{K.FIG}/ranked_curves_{y}.png', dpi=120); plt.close(fig)
    print(y, run, 'real events', len(R), 'noise-template events', len(N), 'per', nn, 'templates; scrambled', len(S))


if __name__ == '__main__':
    if sys.argv[1] == 'fit':
        run_fit()
    else:
        run_year(int(sys.argv[2]), sys.argv[3] if len(sys.argv) > 3 else None)
