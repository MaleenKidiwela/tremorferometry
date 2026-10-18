#!/usr/bin/env python
"""Amendment Y13: do the between-ETS candidates repeat?
  cut Y        (year-Y run env)  -> data/y13/ev_<label>.npz (event windows on the year's MF channels)
  cluster                         -> data/y13/templates.npz, tables/y13_templates.csv (clusters, stacks, C0/C1 similarity)
  mf Y         (year-Y run env)  -> data/y13/rows_<Y>.parquet (template + noise-template single-template detections), LOO
  battery Y    (scan env Y)      -> data/y13/bat_<Y>.parquet
  summary                         -> tables/y13_summary.csv, y13_verdict.json, y13_top_between.csv, figures/y13_top_*.png"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

Y13 = f'{K.DATA}/y13'
os.makedirs(Y13, exist_ok=True)
EVENTS = [('g1', '2016-06-13 07:41'), ('g1', '2020-06-11 13:42'), ('g1', '2020-09-11 05:00'), ('g1', '2020-10-25 19:25'),
          ('g2', '2014-05-05 08:30'), ('g2', '2016-06-13 19:58'), ('g2', '2017-01-31 04:22')]
RUNMAP = {2010: 'prod', 2011: 'prod', 2012: 'prod', 2013: 'prodq', 2014: 'prod', 2015: 'prod', 2016: 'prodq', 2017: 'prodq',
          2018: 'prodq', 2019: 'prodf', 2020: 'prodf', **{y: 'prodfy' for y in range(2021, 2027)}}
NOISE = ['N00', 'N01', 'N02']
THR = 0.7


def events():
    L = pd.read_csv(f'{K.TB}/catalogue_ladder_2010_2026.csv', low_memory=False)
    out = []
    for g, ts in EVENTS:
        r = L[L.time_pt.str.startswith(ts)].iloc[0]
        lab = f"{g}_{ts[:10].replace('-', '')}_{ts[11:].replace(':', '')}"
        out.append(dict(group=g, label=lab, time_pt=r.time_pt, epoch=float(r.epoch), year=int(r.year)))
    return pd.DataFrame(out)


def cut(y):
    import cat_pilot as P
    E = events(); E = E[E.year == y]
    for r in E.itertuples():
        start = int(round(r.epoch - P.EP0)) - P.PRE
        w = P.cut_frame(start); wt = P.rms_ratio_frame(start)
        W = np.array([w[ch] if w[ch] is not None else np.full(P.L, np.nan) for ch in P.CH])
        np.savez(f'{Y13}/ev_{r.label}.npz', chans=np.array(P.CH), W=W, w=np.array([wt[ch] for ch in P.CH]), start=start, year=y)
        print('cut', r.label, int(np.isfinite(W).all(1).sum()), 'channels')


def sim(A, ca, B, cb, lag=30):
    """Y9 step-4 metric on common channels: mean over channels of the max normalised CC within +-lag s."""
    v = []
    for k, ch in enumerate(ca):
        if ch not in cb:
            continue
        a, b = A[k], B[list(cb).index(ch)]
        if not (np.isfinite(a).all() and np.isfinite(b).all()) or np.linalg.norm(a) == 0 or np.linalg.norm(b) == 0:
            continue
        c = np.correlate(a, b, mode='full') / (np.linalg.norm(a) * np.linalg.norm(b)); m = len(a) - 1
        v.append(c[m - lag:m + lag + 1].max())
    return (float(np.mean(v)) if v else np.nan), len(v)


def cluster():
    E = events()
    ev = {r.label: np.load(f'{Y13}/ev_{r.label}.npz') for r in E.itertuples()}
    rows, T = [], {}
    for g, Eg in E.groupby('group'):
        labs = list(Eg.label); n = len(labs)
        S = np.full((n, n), np.nan)
        for i in range(n):
            for j in range(i + 1, n):
                S[i, j] = S[j, i] = sim(ev[labs[i]]['W'], list(ev[labs[i]]['chans']), ev[labs[j]]['W'], list(ev[labs[j]]['chans']))[0]
        pd.DataFrame(S, index=labs, columns=labs).to_csv(f'{K.TB}/y13_similarity_{g}.csv')
        comp = list(range(n))                                  # single linkage at THR
        for i in range(n):
            for j in range(i + 1, n):
                if np.isfinite(S[i, j]) and S[i, j] >= THR:
                    a, b = comp[i], comp[j]; comp = [a if c == b else c for c in comp]
        for c in sorted(set(comp)):
            mem = [labs[i] for i in range(n) if comp[i] == c]
            name = mem[0] if len(mem) == 1 else f'{g}_cl{c}'
            base = ev[mem[0]]; chans = list(base['chans'])
            if len(mem) == 1:
                W, w = base['W'], base['w']
            else:                                              # stack on the first member's channels (members of other years: common channels)
                W = np.full_like(base['W'], np.nan); w = np.full(len(chans), np.nan)
                for k, ch in enumerate(chans):
                    xs, ws = [], []
                    for m in mem:
                        cm = list(ev[m]['chans'])
                        if ch in cm and np.isfinite(ev[m]['W'][cm.index(ch)]).all():
                            xs.append(ev[m]['W'][cm.index(ch)]); ws.append(ev[m]['w'][cm.index(ch)])
                    if xs:
                        s = np.mean(xs, 0); W[k] = s / np.linalg.norm(s); w[k] = np.nanmean(ws)
            T[name] = dict(W=W, w=w, chans=chans, members=mem)
            # similarity with the 2011 C0/C1 templates (bank of the first member's run/year)
            y = int(base['year']); B = np.load(f"{K.DATA}/{RUNMAP[y]}{y}/bank.npz"); bc = list(B['chans'])
            s0 = sim(W, chans, B['C0'], bc); s1 = sim(W, chans, B['C1'], bc)
            rows.append(dict(template=name, group=g, members=','.join(mem), n_members=len(mem), n_ch=int(np.isfinite(W).all(1).sum()),
                             sim_C0=s0[0], sim_C1=s1[0], n_common_C0=s0[1], year=y))
    np.savez(f'{Y13}/templates.npz', **{f'{k}__W': v['W'] for k, v in T.items()}, **{f'{k}__w': v['w'] for k, v in T.items()},
             **{f'{k}__chans': np.array(v['chans']) for k, v in T.items()})
    M = pd.DataFrame(rows); M.to_csv(f'{K.TB}/y13_templates.csv', index=False)
    pd.set_option('display.width', 250); print(M.round(3).to_string(index=False))
    for g in ('g1', 'g2'):
        print(pd.read_csv(f'{K.TB}/y13_similarity_{g}.csv', index_col=0).round(2).to_string())


def mf_one(args):
    import cat_pilot as P
    import cat_e3 as E3
    name, W, w, chans = args
    use = [(k, ch) for k, ch in enumerate(chans) if ch in P.CH and np.isfinite(W[k]).all() and np.isfinite(w[k])]
    if len(use) < 6:
        return name, None
    wk = np.array([w[k] for k, _ in use]); wk = wk * len(wk) / wk.sum()
    S = np.zeros(P.M, np.float32); NL = np.zeros(P.M, np.int16)
    for j, (k, ch) in enumerate(use):
        cc, live = P.ccser(ch, W[k]); S += np.float32(wk[j]) * cc; NL += live
    pi, pv, pn, med, mad = E3.peaks(S, NL)
    return name, dict(pi=pi, z=(pv - med) / mad, n_ch=len(use))


def mf(y):
    import cat_pilot as P
    from multiprocessing import Pool
    Z = np.load(f'{Y13}/templates.npz'); names = sorted({k.split('__')[0] for k in Z.files})
    M = pd.read_csv(f'{K.TB}/y13_templates.csv').set_index('template')
    E = events().set_index('label')
    jobs = [(n, Z[f'{n}__W'], Z[f'{n}__w'], list(Z[f'{n}__chans'])) for n in names]
    # LOO templates for clusters (member years only): stack of the other members
    for n in names:
        mem = M.loc[n, 'members'].split(',')
        if len(mem) < 2:
            continue
        for m in mem:
            if int(E.loc[m, 'year']) != y:
                continue
            others = [np.load(f'{Y13}/ev_{o}.npz') for o in mem if o != m]
            chans = list(others[0]['chans']); W = np.full_like(others[0]['W'], np.nan); w = np.full(len(chans), np.nan)
            for k, ch in enumerate(chans):
                xs = [o['W'][list(o['chans']).index(ch)] for o in others if ch in list(o['chans']) and np.isfinite(o['W'][list(o['chans']).index(ch)]).all()]
                ws = [o['w'][list(o['chans']).index(ch)] for o in others if ch in list(o['chans'])]
                if xs:
                    s = np.mean(xs, 0); W[k] = s / np.linalg.norm(s); w[k] = np.nanmean(ws)
            jobs.append((f'LOO|{n}|{m}', W, w, chans))
    rows, loo = [], []
    with Pool(P.NW) as p:
        res = dict(p.map(mf_one, jobs))
    own = {r.label: int(round(r.epoch - P.EP0)) - P.PRE for r in events().itertuples() if r.year == y}
    for name, o in res.items():
        if o is None:
            continue
        if name.startswith('LOO|'):
            _, n, m = name.split('|'); i0 = own[m]
            hit = o['z'][(np.abs(o['pi'] - i0) <= 30)]
            loo.append(dict(template=n, held_out=m, year=y, z_at_member=float(hit.max()) if len(hit) else np.nan,
                            detected=bool(len(hit) and hit.max() >= 5.5)))
            continue
        excl = np.array([own[m] for m in M.loc[name, 'members'].split(',') if m in own] or [-10 ** 9])
        ok = (o['z'] >= 5.5) & (np.min(np.abs(o['pi'][:, None] - excl[None, :]), axis=1) > 90)
        rows.append(pd.DataFrame(dict(template=name, kind='y13', i=o['pi'][ok] + P.PRE, zmad=o['z'][ok], n_ch=o['n_ch'])))
    Mt = pd.read_csv(f'{K.TB}/{P.NAME}{y}_bank.csv').set_index('name')
    for nz in NOISE:
        f = f'{P.MFD}/{nz}.npz'
        if not os.path.isfile(f):
            continue
        z = np.load(f); zz = (z['real_pv'] - float(z['real_med'])) / float(z['real_mad']); pi = z['real_pi']
        ok = (zz >= 5.5) & (np.abs(pi - int(Mt.loc[nz, 'i0'])) > 90)
        rows.append(pd.DataFrame(dict(template=nz, kind='noise', i=pi[ok], zmad=zz[ok], n_ch=0)))
    R = pd.concat(rows, ignore_index=True)
    R['epoch'] = P.EP0 + R.i; R['year'] = y; R['period'] = P.periods(R.epoch.values)
    R = R.reset_index(drop=True); R['i_row'] = np.arange(len(R))
    R.to_parquet(f'{Y13}/rows_{y}.parquet')
    pd.DataFrame(loo).to_csv(f'{Y13}/loo_{y}.csv', index=False)
    K.log(f'Y13 mf {y}: {len(R)} rows; LOO {len(loo)}', 'y13')


def battery(y):
    import cat_e2
    R = pd.read_parquet(f'{Y13}/rows_{y}.parquet')
    B = R.rename(columns={'i': 'i_frame', 'i_row': 'i'})[['i', 'epoch', 'year', 'template']]
    B.to_parquet(f'{Y13}/brows_{y}.parquet')
    cat_e2.run_battery(y, rows_file=f'{Y13}/brows_{y}.parquet', out_file=f'{Y13}/bat_{y}.parquet', tag='y13')


def summary():
    from scipy.stats import chi2, poisson
    import cat_ychecks as YC
    def gw(n):
        return (0.0 if n == 0 else float(chi2.ppf(.025, 2 * n) / 2), float(chi2.ppf(.975, 2 * n + 2) / 2))
    A = []
    live = {}
    for y in sorted(RUNMAP):
        f = f'{Y13}/bat_{y}.parquet'
        if not os.path.isfile(f):
            continue
        R = pd.read_parquet(f'{Y13}/rows_{y}.parquet'); b = K.apply_rules(pd.read_parquet(f))
        R = R.merge(b[['i', 'pass_a', 'pass_b', 'pass_c', 'pass_e2', 'pass_f', 'a_testable', 'travelling', 'z_best', 'z_loso', 'lat_best', 'lon_best',
                       'ratio_far', 'vr_thrust', 'in_head2']].rename(columns={'i': 'i_row'}), on='i_row', how='left')
        meas = dict(a=R.a_testable.fillna(False).astype(bool), b=R.travelling.notna(), c=R.z_best.notna(), e=pd.Series(True, index=R.index), f=R.z_loso.notna())
        pas = dict(a=R.pass_a, b=R.pass_b, c=R.pass_c, e=R.pass_e2, f=R.pass_f)
        R['level'] = sum((meas[k] & pas[k].fillna(False).astype(bool)).astype(int) for k in 'abcef')
        A.append(R)
        nm = RUNMAP[y]
        live[y] = dict(between=YC.between_years(y, nm + ('q' if nm in ('prodf', 'prodfy') else '')), year_len=(366 if y % 4 == 0 else 365) / 365.25)
    A = pd.concat(A, ignore_index=True)
    yrs = sorted(live)
    T_b = sum(live[y]['between'] for y in yrs); T_all = sum(live[y]['year_len'] for y in yrs)
    nz_ty = A[A.kind == 'noise'].groupby('year').template.nunique().sum()
    ctrl = {p: int(((A.kind == 'noise') & (A.period == p) & (A.level == 5)).sum()) / nz_ty for p in ('ETS', 'between')}
    M = pd.read_csv(f'{K.TB}/y13_templates.csv')
    loo = pd.concat([pd.read_csv(f) for f in sorted(__import__('glob').glob(f'{Y13}/loo_*.csv')) if os.path.getsize(f) > 2] or [pd.DataFrame()])
    rows = []
    for r in M.itertuples():
        x = A[(A.template == r.template) & (A.level == 5)]
        ny = int(A[A.template == r.template].year.nunique())
        o = dict(template=r.template, group=r.group, members=r.members, sim_C0=r.sim_C0, sim_C1=r.sim_C1, years_run=ny)
        for p in ('ETS', 'between'):
            n = int((x.period == p).sum()); exp = ctrl[p] * ny; lo, hi = gw(n)
            o[f'level5_{p}'] = n; o[f'level5_{p}_lo'] = lo; o[f'level5_{p}_hi'] = hi; o[f'ctrl_expect_{p}'] = exp
            o[f'p_{p}'] = float(poisson.sf(n - 1, exp)) if n > 0 and exp > 0 else (0.0 if n > 0 else 1.0)
        if len(loo) and r.template in set(loo.template):
            l = loo[loo.template == r.template]; o['loo'] = f"{int(l.detected.sum())}/{len(l)}"
        o['outcome'] = 'A' if o['p_between'] < 0.01 else ('B' if o['p_ETS'] < 0.01 else 'C')
        rows.append(o)
    S = pd.DataFrame(rows); S.to_csv(f'{K.TB}/y13_summary.csv', index=False)
    top = A[(A.kind == 'y13') & (A.level == 5) & (A.period == 'between')].sort_values('zmad', ascending=False).groupby('template').head(5)
    top = top.assign(time_pt=pd.to_datetime(top.epoch, unit='s').dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.strftime('%Y-%m-%d %H:%M:%S %Z'))
    top.to_csv(f'{K.TB}/y13_top_between.csv', index=False)
    V = dict(noise_template_years=int(nz_ty), ctrl_level5_per_template_year=ctrl, between_live_yr=T_b, years=yrs,
             outcomes=S.set_index('template').outcome.to_dict())
    json.dump(V, open(f'{K.TB}/y13_verdict.json', 'w'), indent=1, default=float)
    pd.set_option('display.width', 250); print(S.round(3).to_string(index=False)); print(json.dumps(V, indent=1, default=float))
    print(top[['template', 'time_pt', 'zmad', 'z_best', 'ratio_far', 'vr_thrust', 'in_head2']].round(2).to_string(index=False))
    # figures of the top between-ETS detections (eye labels added afterwards)
    import cat_ranked3 as C3
    t = top.copy(); t['storm_pct'] = np.nan; t['eq6h_desc'] = ''; t['lfe_status'] = ''; t['n_templates'] = 1; t['far_over_near2'] = np.nan
    if len(t):
        C3.waves(t.reset_index(drop=True), 'y13_top')


if __name__ == '__main__':
    m = sys.argv[1]
    if m in ('cut', 'mf', 'battery'):
        {'cut': cut, 'mf': mf, 'battery': battery}[m](int(sys.argv[2]))
    else:
        {'cluster': cluster, 'summary': summary}[m]()
