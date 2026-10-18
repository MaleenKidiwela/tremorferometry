#!/usr/bin/env python
"""Catalogue E3/E4 (PLAN_CAT §8, Merlin step 4): empirical-template matched filter, 2011-06-01 - 2013-01-01.
Modes
  bank        cluster the 8 vetted + 6 near-miss ETS events by best-node location (complete linkage 15 km), cut 90-s
              velocity windows (origin - 9 s) on all 13 stations, align (+-5 s) and stack per cluster; leave-one-out
              stacks; H&G T2/T4 re-cut; 25 noise templates -> data/e3/bank.npz, tables/e3_bank.csv, e3_events.csv
  mf          matched filter for every (template, station set) (8 workers): primary templates with 8 day-shift + 8
              polarity scrambled controls, LOO and noise templates real only -> data/e3/mf/*.npz
  detect      thresholds per set and year from the noise templates (0.25 / 0.5 / 1 per template per day), detections,
              leave-one-out re-find -> data/e3/rows.parquet, tables/e3_thresholds.csv, e3_refind.csv
  battery Y   battery on every detection time (SCAN_YEAR=Y) -> data/e3/bat_Y.parquet
  summary     -> tables/e3_summary.csv, figures/e3_*.png"""
import os, sys, json
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K
import hg_common as H

E3D = f'{K.DATA}/e3'
MFD = f'{E3D}/mf'
os.makedirs(MFD, exist_ok=True)
HG9 = ['CN.NLLB', 'CN.VGZ', 'UW.LRIV', 'UW.DOSE', 'UW.GNW', 'UW.SP2', 'TA.A04D', 'TA.B05D', 'TA.D03D']
EXTRA = ['CN.PGC', 'CN.SNB', 'CN.PFB', 'CN.YOUB']
ALL = HG9 + EXTRA
SETS = {'all13': ALL, 'hg9': HG9, 'cn4': ['CN.PGC', 'CN.VGZ', 'CN.NLLB', 'CN.SNB'],
        'cn4gnw': ['CN.PGC', 'CN.VGZ', 'CN.NLLB', 'CN.SNB', 'UW.GNW']}
CH = [f'{s}.{c}' for s in ALL for c in 'NEZ']
NF = H.tidx('2013-01-01')                    # frame [2011-06-01, 2013-01-01)
NDF = NF // 86400
MFN = NF - H.L + 1
PRE = 9                                       # window starts at origin - 9 s
NDRAW = 8
NNOISE = 25
RATES = [0.25, 0.5, 1.0]
SEED = 20261015 + 500
_C = {}


def cont(sta):
    if sta not in _C:
        d = f'{K.HGR}/data/r2' if sta in HG9 else f'{K.DATA}/cont'
        _C[sta] = (np.load(f'{d}/contv_{sta}.npy', mmap_mode='r'), np.load(f'{d}/denv_{sta}.npy', mmap_mode='r'))
    return _C[sta]


def cut(start):
    """dict ch -> unit 90-s window at window start, or None if dead."""
    out = {}
    for ch in CH:
        s, c = ch.rsplit('.', 1); ci = 'NEZ'.index(c)
        X, D = cont(s)
        w = np.array(X[ci, start:start + H.L], np.float64)
        if not np.all(np.isfinite(w)) or D[ci, start] <= 0:
            out[ch] = None; continue
        w = w - w.mean(); n = np.linalg.norm(w)
        out[ch] = w / n if n > 0 else None
    return out


def draws():
    rng = np.random.default_rng(SEED)
    offs = []
    for _ in range(NDRAW):
        while True:
            o = rng.choice(np.arange(30, NDF - 30), size=len(ALL), replace=False)
            if np.min(np.diff(np.sort(o))) >= 2:
                break
        offs.append(o.astype(int))
    signs = []
    for _ in range(NDRAW):
        while True:
            sg = rng.choice([-1.0, 1.0], size=len(CH))
            if 0 < (sg > 0).sum() < len(sg):
                break
        signs.append(sg)
    return np.array(offs), np.array(signs)


# ------------------------------------------------------------------------------------------------ bank
def complete_linkage(lat, lon, dmax):
    n = len(lat)
    D = K.dkm(np.asarray(lat)[:, None], np.asarray(lon)[:, None], np.asarray(lat)[None, :], np.asarray(lon)[None, :])
    cl = [[i] for i in range(n)]
    while True:
        best = None
        for a in range(len(cl)):
            for b in range(a + 1, len(cl)):
                d = max(D[i, j] for i in cl[a] for j in cl[b])
                if d <= dmax and (best is None or d < best[0]):
                    best = (d, a, b)
        if best is None:
            break
        _, a, b = best
        cl[a] = cl[a] + cl[b]; del cl[b]
    lab = np.zeros(n, int)
    for k, c in enumerate(sorted(cl, key=lambda c: -len(c))):
        lab[c] = k
    return lab


def stack(windows):
    """windows: list of dict ch -> unit window/None. Channel kept if live in >= half the members."""
    out = {}
    for ch in CH:
        ws = [w[ch] for w in windows if w[ch] is not None]
        if len(ws) == 0 or len(ws) < 0.5 * len(windows):
            out[ch] = None; continue
        s = np.mean(ws, 0); s = s - s.mean(); n = np.linalg.norm(s)
        out[ch] = s / n if n > 0 else None
    return out


def sumcc(a, b):
    return sum(float(a[ch] @ b[ch]) for ch in CH if a[ch] is not None and b[ch] is not None)


def aligned_stack(starts, others_excluded=None):
    """starts: window starts of members; align each to the stack of the others (+-5 s), return (stack, shifts)."""
    W = [cut(s) for s in starts]
    st = stack(W)
    shifts = np.zeros(len(starts), int)
    if len(starts) >= 2:
        for it in range(2):
            for k, s in enumerate(starts):
                ref = stack([W[j] for j in range(len(starts)) if j != k])
                best = max(range(-5, 6), key=lambda d: sumcc(cut(s + d), ref))
                shifts[k] = best
            W = [cut(s + d) for s, d in zip(starts, shifts)]
        st = stack(W)
    return st, shifts, W


def run_bank():
    E = pd.read_csv(f'{K.HGR}/tables/r3_ets_events.csv')
    pa = (E.n_far >= 2) & (E.ratio_far <= 0.551)
    vet = E[pa & E.pass_b & E.pass_c & E.pass_e & E.pass_f].copy(); vet['role'] = 'vetted'
    nm = E[pa & E.pass_b & E.pass_c & E.pass_f & ~E.pass_e].copy(); nm['role'] = 'near-miss'
    V = pd.concat([vet, nm], ignore_index=True)[['cand_id', 'role', 'time_pt', 'i', 'epoch', 'lat_best', 'lon_best',
                                                  'dep_best', 'z_best', 'ratio_far', 'Mw_tpl', 'Mw_scan', 'lfe_status']]
    V['cluster'] = complete_linkage(V.lat_best.values, V.lon_best.values, 15.0)
    V['start0'] = V.i - PRE
    TT = H.templates().set_index('template')
    bank = {}
    meta = []
    # cluster stacks
    for c, g in V.groupby('cluster'):
        st, sh, W = aligned_stack(list(g.start0.values))
        V.loc[g.index, 'shift_s'] = sh
        nm_ = f'C{c}'
        bank[nm_] = st
        meta.append(dict(name=nm_, kind='cluster', members=','.join(map(str, g.cand_id)), n_members=len(g),
                         lat=g.lat_best.mean(), lon=g.lon_best.mean(), n_ch=sum(v is not None for v in st.values()),
                         primary=True))
        # leave-one-out stacks for vetted members
        for k, r in enumerate(g.itertuples()):
            if r.role != 'vetted' or len(g) < 2:
                continue
            idx = [j for j in range(len(g)) if j != k]
            st2 = stack([W[j] for j in idx])
            nm2 = f'C{c}_loo{r.cand_id}'
            bank[nm2] = st2
            meta.append(dict(name=nm2, kind='loo', members=','.join(str(g.cand_id.values[j]) for j in idx),
                             n_members=len(idx), lat=g.lat_best.mean(), lon=g.lon_best.mean(),
                             n_ch=sum(v is not None for v in st2.values()), primary=False, left_out=r.cand_id))
    # self-consistency: member windows vs stack (all channels)
    for r in V.itertuples():
        w = cut(int(r.start0 + r.shift_s))
        V.loc[r.Index, 'cc_with_own_stack'] = sumcc(w, bank[f'C{r.cluster}']) / max(1, sum(
            w[ch] is not None and bank[f'C{r.cluster}'][ch] is not None for ch in CH))
    # H&G T2, T4 re-cut on all stations
    for k in (2, 4):
        io = H.tidx(TT.loc[k].time)
        st = cut(io - PRE)
        bank[f'HG_T{k}'] = st
        meta.append(dict(name=f'HG_T{k}', kind='hg_recut', members='', n_members=1, lat=TT.loc[k].lat, lon=TT.loc[k].lon,
                         n_ch=sum(v is not None for v in st.values()), primary=True))
    # noise templates (first 25 R1 noise times; window start = the noise time, as in R1)
    NT = pd.read_csv(f'{H.MF}/noise_times.csv')
    for j in range(NNOISE):
        st = cut(int(NT.i[j]))
        bank[f'N{j:02d}'] = st
        meta.append(dict(name=f'N{j:02d}', kind='noise', members=str(int(NT.i[j])), n_members=1, lat=np.nan, lon=np.nan,
                         n_ch=sum(v is not None for v in st.values()), primary=False, i0=int(NT.i[j])))
    M = pd.DataFrame(meta)
    # save
    arr = {}
    for nm_, st in bank.items():
        arr[nm_] = np.array([st[ch] if st[ch] is not None else np.full(H.L, np.nan) for ch in CH])
    np.savez(f'{E3D}/bank.npz', chans=np.array(CH), **arr)
    M.to_csv(f'{K.TB}/e3_bank.csv', index=False)
    V.to_csv(f'{K.TB}/e3_events.csv', index=False)
    live = pd.DataFrame({nm_: [st[ch] is not None for ch in CH] for nm_, st in bank.items()}, index=CH)
    live.to_csv(f'{K.TB}/e3_bank_live_channels.csv')
    K.log(f'bank: {len(V)} events in {V.cluster.nunique()} clusters; templates {M.kind.value_counts().to_dict()}', 'e3')
    print(V[['cand_id', 'role', 'time_pt', 'cluster', 'lat_best', 'lon_best', 'shift_s', 'cc_with_own_stack']].to_string())
    print(M.to_string())


# ------------------------------------------------------------------------------------------------ matched filter
def ccser(ch, w):
    from scipy.signal import oaconvolve
    s, c = ch.rsplit('.', 1); ci = 'NEZ'.index(c)
    X, D = cont(s)
    x = np.nan_to_num(np.asarray(X[ci, :NF], np.float32))
    num = oaconvolve(x, w[::-1].astype(np.float32), mode='valid')
    del x
    d = np.asarray(D[ci, :MFN])
    live = d > 0
    cc = np.zeros(MFN, np.float32)
    np.divide(num, d, out=cc, where=live)
    np.clip(cc, -1, 1, out=cc)
    return cc, live


def peaks(S, NL, hw=90):
    from scipy.ndimage import maximum_filter1d
    v = S[NL > 0]
    med = float(np.median(v)); mad = float(np.median(np.abs(v - med)))
    mf = maximum_filter1d(S, 2 * hw + 1, mode='constant', cval=-np.inf)
    i = np.where((S == mf) & (S > med + 3 * mad) & (NL > 0))[0]
    if i.size > 1:
        i = i[np.concatenate([[True], np.diff(i) > hw])]
    return i.astype(np.int64), S[i].astype(np.float32), NL[i].astype(np.int16), med, mad


def mf_job(job):
    name, setname, controls = job
    f = f'{MFD}/{name}__{setname}.npz'
    if os.path.isfile(f):
        return f'{name} {setname} exists'
    B = np.load(f'{E3D}/bank.npz')
    W = B[name]
    chs = [ch for k, ch in enumerate(CH) if ch.rsplit('.', 1)[0] in SETS[setname] and np.isfinite(W[k]).all()]
    names = ['real'] + ([f'A1_{r}' for r in range(NDRAW)] + [f'A2_{r}' for r in range(NDRAW)] if controls else [])
    S = {n: np.zeros(MFN, np.float32) for n in names}
    NL = {n: np.zeros(MFN, np.int16) for n in names}
    offs, signs = draws()
    for ch in chs:
        k = CH.index(ch); si = ALL.index(ch.rsplit('.', 1)[0])
        cc, live = ccser(ch, W[k])
        S['real'] += cc; NL['real'] += live
        if controls:
            for r in range(NDRAW):
                sh = int(offs[r, si]) * 86400
                S[f'A1_{r}'] += np.roll(cc, sh); NL[f'A1_{r}'] += np.roll(live, sh)
                S[f'A2_{r}'] += cc * np.float32(signs[r, k]); NL[f'A2_{r}'] += live
        del cc, live
    out = {'name': name, 'set': setname, 'chans': np.array(chs)}
    for n in names:
        pi, pv, pn, med, mad = peaks(S[n], NL[n])
        out[f'{n}_pi'] = pi; out[f'{n}_pv'] = pv; out[f'{n}_pnl'] = pn; out[f'{n}_med'] = med; out[f'{n}_mad'] = mad
        S[n] = None; NL[n] = None
    np.savez_compressed(f + '.part.npz', **out); os.replace(f + '.part.npz', f)
    return f'{name} {setname} done ({len(chs)} ch, controls {controls})'


def run_mf():
    M = pd.read_csv(f'{K.TB}/e3_bank.csv')
    jobs = []
    for r in M.itertuples():
        for s in SETS:
            jobs.append((r.name, s, bool(r.primary)))
    jobs.sort(key=lambda j: (not j[2], j[1] != 'all13'))      # controls first (heaviest)
    from multiprocessing import Pool
    with Pool(8, maxtasksperchild=4) as p:
        for k, r in enumerate(p.imap_unordered(mf_job, jobs)):
            K.log(f'mf {k + 1}/{len(jobs)}: {r}', 'e3')


# ------------------------------------------------------------------------------------------------ detect
def load_mf(name, s):
    return np.load(f'{MFD}/{name}__{s}.npz')


def year_of(pi):
    t = H.T0 + pd.to_timedelta(np.asarray(pi), unit='s')
    return np.asarray(t.year)


def run_detect():
    M = pd.read_csv(f'{K.TB}/e3_bank.csv')
    V = pd.read_csv(f'{K.TB}/e3_events.csv')
    ndays = {2011: (pd.Timestamp('2012-01-01') - H.T0).days, 2012: 366}
    TH, rows, refind = [], [], []
    noise = M[M.kind == 'noise']
    for s in SETS:
        # thresholds from the 25 noise templates (own median/MAD), per year
        zs = {2011: [], 2012: []}
        for r in noise.itertuples():
            z = load_mf(r.name, s)
            med, mad = float(z['real_med']), float(z['real_mad'])
            pi, pv = z['real_pi'], z['real_pv']
            m = np.abs(pi - int(r.i0)) > 90
            yy = year_of(pi)
            for y in (2011, 2012):
                zs[y].append(((pv[m & (yy == y)] - med) / mad))
        zstar = {}
        for y in (2011, 2012):
            allz = np.sort(np.concatenate(zs[y]))[::-1]
            for rate in RATES:
                n = int(round(rate * len(noise) * ndays[y]))
                zstar[(y, rate)] = float(allz[n - 1]) if 0 < n <= len(allz) else float(allz[-1])
                TH.append(dict(set=s, year=y, rate_per_tpl_day=rate, z_star_MAD=zstar[(y, rate)], n_noise=len(noise)))
        zlow = {y: zstar[(y, 1.0)] for y in (2011, 2012)}
        # detections: every template (real + controls), at the lowest threshold (rate 1/day); rate columns
        for r in M.itertuples():
            z = load_mf(r.name, s)
            med, mad = float(z['real_med']), float(z['real_mad'])
            srcs = ['real'] + ([f'A1_{q}' for q in range(NDRAW)] + [f'A2_{q}' for q in range(NDRAW)] if r.primary else [])
            for src in srcs:
                pi, pv, pn = z[f'{src}_pi'], z[f'{src}_pv'], z[f'{src}_pnl']
                zz = (pv - med) / mad
                yy = year_of(pi)
                thr = np.where(yy == 2011, zlow[2011], zlow[2012])
                m = zz >= thr
                if r.kind == 'noise':
                    m &= np.abs(pi - int(r.i0)) > 90
                if not m.any():
                    continue
                d = pd.DataFrame(dict(set=s, template=r.name, kind=r.kind, source=src if src == 'real' else src[:2],
                                      draw=-1 if src == 'real' else int(src[3:]), p=pi[m], zmad=zz[m], n_live=pn[m],
                                      nch=len(z['chans']), year=yy[m]))
                for rate in RATES:
                    d[f'det_{rate}'] = d.zmad >= np.where(d.year == 2011, zstar[(2011, rate)], zstar[(2012, rate)])
                rows.append(d)
        # leave-one-out re-find (vetted events only)
        for e in V[V.role == 'vetted'].itertuples():
            ye = int(year_of([e.i])[0])
            cand_tpl = [t for t in M.itertuples() if t.kind in ('cluster', 'loo', 'hg_recut')
                        and str(e.cand_id) not in str(t.members).split(',')
                        and not (t.kind == 'loo' and str(getattr(t, 'left_out', '')) not in ('', 'nan') and int(float(t.left_out)) != e.cand_id)]
            best = (-np.inf, '')
            per = {}
            for t in cand_tpl:
                z = load_mf(t.name, s)
                med, mad = float(z['real_med']), float(z['real_mad'])
                pi, pv = z['real_pi'], z['real_pv']
                m = np.abs(pi + PRE - e.i) <= 30
                zz = float(((pv[m] - med) / mad).max()) if m.any() else np.nan
                per[t.name] = zz
                if np.isfinite(zz) and zz > best[0]:
                    best = (zz, t.name)
            o = dict(set=s, cand_id=e.cand_id, time_pt=e.time_pt, cluster=e.cluster, best_template=best[1],
                     best_z=best[0], templates_tried=len(cand_tpl))
            for rate in RATES:
                o[f'refound_{rate}'] = bool(best[0] >= zstar[(ye, rate)])
            o['per_template_z'] = json.dumps({k: (round(v, 2) if np.isfinite(v) else None) for k, v in per.items()})
            refind.append(o)
        K.log(f'detect {s}: z* {json.dumps({f"{y}_{r}": round(v, 2) for (y, r), v in zstar.items()})}', 'e3')
    R = pd.concat(rows, ignore_index=True)
    EP0 = (H.T0 - pd.Timestamp('1970-01-01')).total_seconds()
    R['i'] = R.p + PRE
    R['epoch'] = EP0 + R.i
    R['year'] = pd.to_datetime(R.epoch, unit='s').dt.year
    R['period'] = K.period_of(R.epoch.values)
    R.to_parquet(f'{E3D}/rows.parquet')
    pd.DataFrame(TH).to_csv(f'{K.TB}/e3_thresholds.csv', index=False)
    F = pd.DataFrame(refind)
    F.to_csv(f'{K.TB}/e3_refind.csv', index=False)
    K.log(f'detect: {len(R)} rows ({R.groupby(["set", "source"]).size().to_dict()}); unique times {R.i.nunique()}', 'e3')
    print(F.drop(columns='per_template_z').to_string())


def run_battery(year):
    import cat_e2
    cat_e2.run_battery(year, rows_file=f'{E3D}/rows.parquet', out_file=f'{E3D}/bat_{year}.parquet', tag='e3')


# ------------------------------------------------------------------------------------------------ summary
def dedupe_count(i):
    t = np.sort(np.asarray(i))
    return int(1 + np.sum(np.diff(t) > 30)) if len(t) else 0


def run_summary():
    import cat_e2
    R = pd.read_parquet(f'{E3D}/rows.parquet')
    B = pd.concat([pd.read_parquet(f'{E3D}/bat_{y}.parquet') for y in (2011, 2012)], ignore_index=True).drop(columns=['epoch'])
    B['grid_edge'] = B.grid_edge.where(B.grid_edge.notna(), B.grid_edge_node) if 'grid_edge' in B else B.grid_edge_node
    D = R.merge(B, on='i', how='left')
    for c in ('n_far', 'n6'):
        D[c] = D[c].fillna(0)
    D = K.apply_rules(D)
    D['lfe_status'], D['n_lfe_10m'] = cat_e2.lfe_status(D.epoch.values)
    D.to_parquet(f'{E3D}/e3_rows_scored.parquet')
    F = pd.read_csv(f'{K.TB}/e3_refind.csv')
    yr = K.BETWEEN_YR
    out = []
    for s in SETS:
        for rate in RATES:
            for bcol in ('battery_e2', 'battery'):
                x = D[(D.set == s) & D[f'det_{rate}'] & D[bcol]]
                xp = x[x.kind.isin(['cluster', 'hg_recut'])]
                o = dict(set=s, noise_rate_per_tpl_day=rate, battery="e'" if bcol == 'battery_e2' else 'e')
                for per in ('ETS', 'between'):
                    y = xp[xp.period == per]
                    o[f'real_{per}'] = dedupe_count(y[y.source == 'real'].i)
                    cc = [dedupe_count(y[(y.source == a) & (y.draw == q)].i) for a in ('A1', 'A2') for q in range(NDRAW)]
                    o[f'ctrl_mean_{per}'] = float(np.mean(cc)); o[f'ctrl_max_{per}'] = int(np.max(cc))
                    yn = x[(x.kind == 'noise') & (x.period == per)]
                    o[f'noise_per_tpl_{per}'] = float(np.mean([dedupe_count(yn[yn.template == f'N{j:02d}'].i) for j in range(NNOISE)]))
                o['ctrl_rate_between_per_yr'] = o['ctrl_mean_between'] / yr
                o['noise_between_per_tpl_yr'] = o['noise_per_tpl_between'] / yr
                f = F[F.set == s]
                o['loo_refound'] = int(f[f'refound_{rate}'].sum()); o['loo_n'] = len(f)
                o['detections_real_between_noBattery'] = dedupe_count(D[(D.set == s) & D[f'det_{rate}'] & (D.source == 'real')
                                                                         & D.kind.isin(['cluster', 'hg_recut']) & (D.period == 'between')].i)
                out.append(o)
    T = pd.DataFrame(out)
    T['accept_all13'] = np.where(T.set == 'all13', T.loo_refound >= 7, np.nan)
    T['accept_4sta_rule'] = (T.loo_refound >= 5) & (T.ctrl_rate_between_per_yr <= 1)
    T.to_csv(f'{K.TB}/e3_summary.csv', index=False)
    pd.set_option('display.width', 250)
    print(T.to_string(index=False))
    ev = D[(D.source == 'real') & D.kind.isin(['cluster', 'hg_recut']) & D['det_0.5'] & D.battery_e2].sort_values('epoch')
    ev = ev.assign(time_pt=pd.to_datetime(ev.epoch, unit='s').dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.strftime('%Y-%m-%d %H:%M:%S %Z'))
    ev[['set', 'template', 'time_pt', 'period', 'zmad', 'z_best', 'ratio_far', 'n_far', 'ratio_GNW', 'lat_best', 'lon_best',
        'nearest_family', 'lfe_status', 'n_lfe_10m', 'head2_maxmag', 'pass_e', 'i']].to_csv(f'{K.TB}/e3_real_battery_events.csv', index=False)
    figs(D, F)


def figs(D, F):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 4, figsize=(18, 8), sharey=True)
    for b, s in enumerate(SETS):
        for a, per in enumerate(('ETS', 'between')):
            ax = axs[a, b]
            x = D[(D.set == s) & (D.period == per) & D['det_1.0']]
            for src, c, sz in (('noise', '0.7', 4), ('ctrl', 'tab:purple', 4), ('real', 'tab:red', 10)):
                y = x[x.kind == 'noise'] if src == 'noise' else (x[x.source.isin(['A1', 'A2'])] if src == 'ctrl'
                                                                   else x[(x.source == 'real') & (x.kind != 'noise') & (x.kind != 'loo')])
                ax.scatter(y.zmad, y.z_best.clip(0.5, 100), s=sz, c=c, lw=0, alpha=0.6, label=f'{src} {len(y)}')
                yb = y[y.battery_e2]
                ax.scatter(yb.zmad, yb.z_best, s=40, facecolors='none', edgecolors=c if src != 'noise' else 'k', lw=1)
            ax.set_yscale('log'); ax.axhline(5, color='k', lw=0.5)
            ax.set_title(f'{s} — {per} (rings = battery′)', fontsize=9); ax.legend(fontsize=7)
            if a == 1:
                ax.set_xlabel('MF z (MAD)')
        axs[0, 0].set_ylabel('scan z_best'); axs[1, 0].set_ylabel('scan z_best')
    fig.tight_layout(); fig.savefig(f'{K.FIG}/e3_mf_vs_scanz.png', dpi=120)
    fig, ax = plt.subplots(figsize=(8, 4))
    for k, s in enumerate(SETS):
        f = F[F.set == s]
        ax.scatter(np.arange(len(f)) + 0.15 * k, f.best_z, label=s, s=30)
    ax.set_xticks(range(len(F[F.set == 'all13']))); ax.set_xticklabels(F[F.set == 'all13'].time_pt.str[:16], rotation=60, fontsize=7)
    ax.set_ylabel('best leave-one-out MF z (MAD)'); ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig(f'{K.FIG}/e3_loo_refind.png', dpi=130)


if __name__ == '__main__':
    mode = sys.argv[1]
    {'bank': run_bank, 'mf': run_mf, 'detect': run_detect, 'summary': run_summary}.get(mode, lambda: None)()
    if mode == 'battery':
        run_battery(int(sys.argv[2]))
