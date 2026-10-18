#!/usr/bin/env python
"""Y13 verification (RESULTS 27): C0/C1 as single templates on the same footing, overlap, top detections, timing.
  fair     -> tables/y13b_fair.csv (level-5 single-template detections ETS/between: g1, g2, C0, C1, noise), data/y13/levels_all.parquet
  overlap  -> tables/y13b_overlap.json, tables/y13b_top_dedup.csv (top-10 per g1 template, deduplicated, with HF and catalogue)
  figs     -> figures/y13b_top_<k>.png
  timing   -> figures/y13b_timing.png, tables/y13b_timing.json"""
import os, sys, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

Y13 = f'{K.DATA}/y13'
RUNMAP = {2010: 'prod', 2011: 'prod', 2012: 'prod', 2013: 'prodq', 2014: 'prod', 2015: 'prod', 2016: 'prodq', 2017: 'prodq',
          2018: 'prodq', 2019: 'prodf', 2020: 'prodf', **{y: 'prodfy' for y in range(2021, 2027)}}


def levels(R):
    meas = dict(a=R.a_testable.fillna(False).astype(bool), b=R.travelling.notna(), c=R.z_best.notna(), e=pd.Series(True, index=R.index), f=R.z_loso.notna())
    pas = dict(a=R.pass_a, b=R.pass_b, c=R.pass_c, e=R.pass_e2, f=R.pass_f)
    return sum((meas[k] & pas[k].fillna(False).astype(bool)).astype(int) for k in 'abcef')


def load_y13():
    A = []
    for y in sorted(RUNMAP):
        f = f'{Y13}/bat_{y}.parquet'
        if not os.path.isfile(f):
            continue
        R = pd.read_parquet(f'{Y13}/rows_{y}.parquet'); b = K.apply_rules(pd.read_parquet(f))
        R = R.merge(b.drop(columns=['epoch', 'year', 'template'], errors='ignore').rename(columns={'i': 'i_row'}), on='i_row', how='left')
        R['level'] = levels(R); A.append(R)
    return pd.concat(A, ignore_index=True)


def load_c0c1():
    A = []
    for y, run in RUNMAP.items():
        d = f'{K.DATA}/{run}{y}'
        if not (os.path.isfile(f'{d}/rows.parquet') and os.path.isfile(f'{d}/bat_v2.parquet')):
            continue
        R = pd.read_parquet(f'{d}/rows.parquet'); R = R[(R.source == 'real') & R.template.isin(['C0', 'C1']) & (R.zmad >= 5.5)]
        B = pd.read_parquet(f'{d}/bat_v2.parquet').drop(columns=['epoch'])
        B['grid_edge'] = B.grid_edge.where(B.grid_edge.notna(), B.grid_edge_node) if 'grid_edge' in B else B.grid_edge_node
        R = R.merge(B, on='i', how='left')
        for c in ('n_far', 'n6'):
            R[c] = R[c].fillna(0)
        R = K.apply_rules(R); R['level'] = levels(R); R['kind'] = 'C0C1'; R['year'] = y
        A.append(R)
    return pd.concat(A, ignore_index=True)


def fair():
    A = load_y13(); C = load_c0c1()
    A.to_parquet(f'{Y13}/levels_y13.parquet'); C.to_parquet(f'{Y13}/levels_c0c1.parquet')
    rows = []
    for t, g in A.groupby('template'):
        rows.append(dict(template=t, kind=g.kind.iloc[0], years=int(g.year.nunique()), L5_ETS=int(((g.level == 5) & (g.period == 'ETS')).sum()),
                         L5_between=int(((g.level == 5) & (g.period == 'between')).sum()), det_between=int((g.period == 'between').sum())))
    for t, g in C.groupby('template'):
        rows.append(dict(template=t, kind='C0C1 (production rows)', years=int(g.year.nunique()), L5_ETS=int(((g.level == 5) & (g.period == 'ETS')).sum()),
                         L5_between=int(((g.level == 5) & (g.period == 'between')).sum()), det_between=int((g.period == 'between').sum())))
    T = pd.DataFrame(rows)
    nz = T[T.kind == 'noise']; nyr = A[A.kind == 'noise'].groupby('template').year.nunique().sum()
    T.loc[len(T)] = dict(template='noise per template (17 yr)', kind='noise', years=17, L5_ETS=nz.L5_ETS.sum() / nyr * 17,
                         L5_between=nz.L5_between.sum() / nyr * 17, det_between=nz.det_between.sum() / nyr * 17)
    T.to_csv(f'{K.TB}/y13b_fair.csv', index=False)
    py = pd.concat([A[(A.level == 5) & (A.period == 'between')].groupby(['template', 'year']).size().unstack(0),
                    C[(C.level == 5) & (C.period == 'between')].groupby(['template', 'year']).size().unstack(0)], axis=1).fillna(0).astype(int)
    py.to_csv(f'{K.TB}/y13b_fair_by_year.csv')
    pd.set_option('display.width', 250); print(T.round(1).to_string(index=False)); print(py.to_string())


def overlap():
    A = pd.read_parquet(f'{Y13}/levels_y13.parquet'); C = pd.read_parquet(f'{Y13}/levels_c0c1.parquet')
    g1 = A[(A.template.str.startswith('g1')) & (A.level == 5) & (A.period == 'between')].sort_values('epoch')
    g1 = g1.assign(grp=(g1.epoch.diff().fillna(1e9) > 30).cumsum())
    ev = g1.groupby('grp').agg(epoch=('epoch', 'first'), n_tpl=('template', 'nunique'), tpls=('template', lambda s: ','.join(sorted(set(s)))),
                               zmad=('zmad', 'max'), z_best=('z_best', 'max')).reset_index(drop=True)
    c5 = np.sort(C[(C.level == 5)].epoch.values); call = np.sort(C.epoch.values)
    def near(e, ref):
        k = np.searchsorted(ref, e); d = np.minimum(np.abs(e - ref[np.clip(k - 1, 0, len(ref) - 1)]), np.abs(e - ref[np.clip(k, 0, len(ref) - 1)])); return d <= 30
    ev['with_C0C1_L5'] = near(ev.epoch.values, c5); ev['with_C0C1_any'] = near(ev.epoch.values, call)
    S = dict(g1_rows=len(g1), g1_distinct_events=len(ev), by_n_templates=ev.n_tpl.value_counts().sort_index().to_dict(),
             with_C0C1_level5=int(ev.with_C0C1_L5.sum()), with_C0C1_any_MF=int(ev.with_C0C1_any.sum()))
    json.dump(S, open(f'{K.TB}/y13b_overlap.json', 'w'), indent=1, default=int)
    ev.to_csv(f'{K.TB}/y13b_g1_between_events.csv', index=False)
    print(json.dumps(S, indent=1, default=int))
    # top-10 per template, deduplicated
    top = []
    for t, g in g1.groupby('template'):
        top.append(g.sort_values('zmad', ascending=False).head(10))
    top = pd.concat(top).sort_values('epoch'); top = top.assign(grp=(top.epoch.diff().fillna(1e9) > 30).cumsum())
    top = top.sort_values('zmad', ascending=False).groupby('grp').head(1).sort_values('epoch').reset_index(drop=True)
    top['n_tpl_all'] = [int(ev.loc[(ev.epoch - e).abs().idxmin(), 'n_tpl']) for e in top.epoch]
    top['with_C0C1_any'] = near(top.epoch.values, call); top['with_C0C1_L5'] = near(top.epoch.values, c5)
    top['time_pt'] = pd.to_datetime(top.epoch, unit='s').dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.strftime('%Y-%m-%d %H:%M:%S %Z')
    # catalogue and high frequency
    import cat_y9b as B
    EQ = pd.concat([K.eq_table(y) for y in sorted(top.year.unique())]).drop_duplicates('id')
    eqs, hfx, hfn = [], [], []
    for r in top.itertuples():
        q = EQ[(EQ.t >= r.epoch - 3600) & (EQ.t <= r.epoch + 300)].copy()
        q['tR0'] = q.t + q.dist_km / 4.5; q['tR1'] = q.t + q.dist_km / 2.8 + 600
        hit = q[((q.tR0 <= r.epoch + 180) & (q.tR1 >= r.epoch)) | ((q.t >= r.epoch - 300) & (q.t <= r.epoch + 300))].sort_values('a20_pred2_nm', ascending=False)
        eqs.append(f"M{hit.mag.iloc[0]:.1f} {hit.dist_km.iloc[0]:.0f}km {(r.epoch - hit.t.iloc[0]) / 60:.0f}min" if len(hit) else '')
        X = [hfext(s, r.epoch) for s in B.HFSTA]; X = np.array([x for x in X if np.isfinite(x)])
        hfx.append(float(X.max()) if len(X) else np.nan); hfn.append(int((X >= 3).sum()) if len(X) else 0)
    top['eq_in_window_catalogue'] = eqs; top['hf28_ext_max'] = hfx; top['hf28_ext_n_ge3'] = hfn
    top.to_csv(f'{K.TB}/y13b_top_dedup.csv', index=False)
    print(top[['time_pt', 'template', 'n_tpl_all', 'zmad', 'z_best', 'ratio_far', 'in_head2', 'with_C0C1_any', 'eq_in_window_catalogue', 'hf28_ext_max', 'hf28_ext_n_ge3']].round(2).to_string())


def hfext(sta, ep):
    import s13_common as C
    from obspy import read, UTCDateTime
    t0 = UTCDateTime(ep); d = pd.Timestamp(ep, unit='s'); fs = C.day_files(sta, d.strftime('%Y-%m-%d'))
    if not fs:
        return np.nan
    try:
        st = read(fs[0], starttime=t0 - 1900, endtime=t0 + 1900).select(component='Z'); st.merge(fill_value=0); tr = st[0]
    except Exception:
        return np.nan
    sr = tr.stats.sampling_rate
    if sr < 20:
        return np.nan
    tr.detrend('demean'); tr.taper(0.02); tr.filter('bandpass', freqmin=2, freqmax=8, corners=4, zerophase=True)
    y = tr.data.astype(float); n0 = int((t0 - 1800 - tr.stats.starttime) * sr); w = int(30 * sr)
    seg = np.array([np.sqrt(np.mean(y[n0 + k * w:n0 + (k + 1) * w] ** 2)) for k in range(120) if n0 + (k + 1) * w <= len(y)])
    a = int((t0 - 300 - tr.stats.starttime) * sr)
    win = [np.sqrt(np.mean(y[a + k * w:a + (k + 1) * w] ** 2)) for k in range(16) if a + (k + 1) * w <= len(y)]
    return float(np.max(win) / np.median(seg)) if len(seg) and win and np.median(seg) > 0 else np.nan


def figs():
    import cat_ranked3 as C3
    t = pd.read_csv(f'{K.TB}/y13b_top_dedup.csv')
    t['storm_pct'] = np.nan; t['lfe_status'] = ''; t['n_templates'] = t.n_tpl_all; t['far_over_near2'] = np.nan
    t['eq6h_desc'] = t.eq_in_window_catalogue.fillna('') + '  HF ' + t.hf28_ext_max.round(1).astype(str) + ' (' + t.hf28_ext_n_ge3.astype(str) + ' st>=3)  ' + t.template
    C3.waves(t, 'y13b_top')


def timing():
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ev = pd.read_csv(f'{K.TB}/y13b_g1_between_events.csv')
    E = pd.read_csv(f'{K.ROOT}/svi_science/results/b_ets_record/episodes_b.csv'); E = E[E.region == 'SVI']
    ws = [((pd.Timestamp(r.start) - pd.Timedelta(days=2)).value / 1e9, (pd.Timestamp(r.end) + pd.Timedelta(days=3)).value / 1e9) for r in E.itertuples()]
    def edge_d(e):
        return min(min(abs(e - a), abs(e - b)) for a, b in ws) / 86400
    ev['days_to_ets_window'] = [edge_d(e) for e in ev.epoch]
    P = pd.read_csv(f'{K.ROOT}/catalogs/pnsn_tremor_cascadia_full.csv', usecols=['time', 'lat', 'lon'])
    P['t'] = pd.to_datetime(P.time).astype('datetime64[ns]').astype('int64') / 1e9
    pc = pd.read_csv(f'{K.TB}/inj_sources.csv').set_index('source').loc[['C0', 'C1'], ['lat', 'lon']].mean()
    P = P[np.hypot((P.lat - pc.lat) * 111.19, (P.lon - pc.lon) * 111.19 * np.cos(np.radians(pc.lat))) <= 50]
    tr = np.sort(P.t.values)
    ev['tremor50km_pm1d'] = np.searchsorted(tr, ev.epoch + 86400) - np.searchsorted(tr, ev.epoch - 86400)
    ev['tremor50km_pm7d'] = np.searchsorted(tr, ev.epoch + 7 * 86400) - np.searchsorted(tr, ev.epoch - 7 * 86400)
    # null: random between-ETS times in the same years
    rng = np.random.default_rng(3); yrs = pd.to_datetime(ev.epoch, unit='s').dt.year
    lo, hi = pd.Timestamp(f'{yrs.min()}-01-01').value / 1e9, pd.Timestamp(f'{yrs.max() + 1}-01-01').value / 1e9
    rnd = rng.uniform(lo, hi, 20000); rnd = rnd[[not any(a <= x < b for a, b in ws) for x in rnd]][:3000]
    S = dict(n_events=len(ev), median_days_to_ets_window=float(ev.days_to_ets_window.median()),
             median_days_random=float(np.median([edge_d(x) for x in rnd[:1000]])),
             frac_within_7d=float((ev.days_to_ets_window <= 7).mean()), frac_within_7d_random=float(np.mean([edge_d(x) <= 7 for x in rnd[:1000]])),
             frac_tremor_pm1d=float((ev.tremor50km_pm1d > 0).mean()),
             frac_tremor_pm1d_random=float(np.mean((np.searchsorted(tr, rnd + 86400) - np.searchsorted(tr, rnd - 86400)) > 0)),
             per_year=pd.to_datetime(ev.epoch, unit='s').dt.year.value_counts().sort_index().to_dict(),
             top_days=pd.to_datetime(ev.epoch, unit='s').dt.date.astype(str).value_counts().head(8).to_dict())
    json.dump(S, open(f'{K.TB}/y13b_timing.json', 'w'), indent=1, default=float)
    ev.to_csv(f'{K.TB}/y13b_g1_between_events.csv', index=False)
    fig, ax = plt.subplots(2, 1, figsize=(14, 6), sharex=True)
    t = pd.to_datetime(ev.epoch, unit='s')
    ax[0].hist(t, bins=pd.date_range('2010-01-01', '2026-11-01', freq='MS'), color='tab:red')
    for a, b in ws:
        if a > pd.Timestamp('2010-01-01').value / 1e9:
            for x in ax:
                x.axvspan(pd.Timestamp(a, unit='s'), pd.Timestamp(b, unit='s'), color='0.85', zorder=0)
    ax[0].set_ylabel('g1 level-5 between-ETS events / month')
    ax[1].hist(pd.to_datetime(tr, unit='s'), bins=pd.date_range('2010-01-01', '2026-11-01', freq='MS'), color='tab:orange'); ax[1].set_yscale('log')
    ax[1].set_ylabel('PNSN tremor <= 50 km / month'); ax[0].set_title('grey = SVI ETS windows (start - 2 d .. end + 3 d)', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{K.FIG}/y13b_timing.png', dpi=110); plt.close(fig)
    print(json.dumps(S, indent=1, default=float))


if __name__ == '__main__':
    {'fair': fair, 'overlap': overlap, 'figs': figs, 'timing': timing}[sys.argv[1]]()
