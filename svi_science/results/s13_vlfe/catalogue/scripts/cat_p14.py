#!/usr/bin/env python
"""Pilot 2014 (PLAN_CAT Amendment P2; Merlin step 6). Frame 2014-01-01 - 2015-01-01, H&G 9 stations, V4 recipe.
  cont        -> data/p14/cont/contv_STA.npy, denv_STA.npy (velocity 0.02-0.05 Hz, from disp1sps_v3 2014 day files)
  bank        -> data/p14/bank.npz, tables/p14_bank.csv (H&G T1-T5 2011 windows, E3 C0-C5, H&G T6-T13, 25 noise templates)
  mf          -> data/p14/mf/NAME.npz (bank: real + 8 day-shift + 8 polarity; noise: real only)
  detect      -> data/p14/rows.parquet, tables/p14_noise_threshold.csv
  battery     (SCAN_YEAR=2014) -> data/p14/bat_2014.parquet
  truth       -> tables/p14_truth.csv, p14_counts.csv, p14_verdict.json (needs inj)
  inj_list / inj (SCAN_YEAR=2014) -> data/p14/inj_list.parquet, inj.parquet
  summary     -> tables/p14_*, figures/p14_*"""
import os, sys, json
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

P14 = f'{K.DATA}/p14'
CONT = f'{P14}/cont'
MFD = f'{P14}/mf'
for d in (P14, CONT, MFD):
    os.makedirs(d, exist_ok=True)
HG9 = ['CN.NLLB', 'CN.VGZ', 'UW.LRIV', 'UW.DOSE', 'UW.GNW', 'UW.SP2', 'TA.A04D', 'TA.B05D', 'TA.D03D']
CH = [f'{s}.{c}' for s in HG9 for c in 'NEZ']
T0 = pd.Timestamp('2014-01-01')
NDAY = 365
N = NDAY * 86400
L = 90
M = N - L + 1
EP0 = (T0 - pd.Timestamp('1970-01-01')).total_seconds()
ETS = ('2014-11-12', '2014-12-17')            # SVI 2014-11 episode (11-14 .. 12-14) +- 2 d
NDRAW, NNOISE, NW = 8, 25, 8
SEED = 20261015 + 900
LEVEL = 5.5
VALID = [f'T{k}' for k in range(1, 6)] + [f'C{k}' for k in range(6)]     # Q0: validated bank; T6-T13 columns only


def tidx(t):
    return int(round((pd.Timestamp(t) - T0).total_seconds()))


def periods(epoch):
    t = pd.to_datetime(np.asarray(epoch, float), unit='s')
    out = np.full(len(t), 'between', dtype=object)
    E = pd.read_csv(f'{K.ROOT}/svi_science/results/b_ets_record/episodes_b.csv')
    E = E[E.region.isin(['SVI', 'NW'])]
    for r in E.itertuples():
        a = pd.Timestamp(r.start) - pd.Timedelta(days=2); b = pd.Timestamp(r.end) + pd.Timedelta(days=3)
        out[(t >= a) & (t < b)] = 'other_episode'
    out[(t >= pd.Timestamp(ETS[0])) & (t < pd.Timestamp(ETS[1]))] = 'ETS'
    return out


def between_years():
    days = pd.date_range(T0, periods=NDAY, freq='D') + pd.Timedelta(hours=12)
    p = periods((days - pd.Timestamp('1970-01-01')).total_seconds().values)
    return float((p == 'between').sum() / 365.25)


# ------------------------------------------------------------------------------------------------ cont
def dens(X):
    from hg_cont import den_chunked
    import hg_common as H
    Mh = H.M
    H.M = M                                                  # den_chunked reads hg_common.M; frame-local override
    try:
        D = np.zeros((3, M), np.float32)
        for c in range(3):
            den = den_chunked(X[c])
            med = np.median(den[den > 0]) if np.any(den > 0) else 0
            den[den < 1e-4 * med] = 0.0
            D[c] = den.astype(np.float32)
    finally:
        H.M = Mh
    return D


def build_cont(sta):
    import hg_cont
    f = f'{CONT}/contv_{sta}.npy'
    if os.path.isfile(f) and os.path.isfile(f'{CONT}/denv_{sta}.npy'):
        return sta, 'exists'
    X = np.full((3, N), np.nan, np.float32)
    nd = 0
    for k in range(NDAY):
        d = T0 + pd.Timedelta(days=k)
        p = f'{K.S13}/data/disp1sps_v3/{sta}/{d.year}/{d.dayofyear:03d}.npz'
        if os.path.isfile(p):
            X[:, k * 86400:(k + 1) * 86400] = np.load(p)['d']; nd += 1
    for c in range(3):
        x = hg_cont.filt(X[c]).astype(np.float64)
        v = np.full(x.shape, np.nan); v[1:-1] = 0.5 * (x[2:] - x[:-2])
        X[c] = v.astype(np.float32)
    np.save(f + '.part.npy', X); os.replace(f + '.part.npy', f)
    D = dens(X)
    np.save(f'{CONT}/denv_{sta}.part.npy', D); os.replace(f'{CONT}/denv_{sta}.part.npy', f'{CONT}/denv_{sta}.npy')
    return sta, f'days {nd}/{NDAY} live N/E/Z {(D > 0).mean(axis=1).round(3)}'


def run_cont():
    from multiprocessing import Pool
    with Pool(5) as p:
        for s, r in p.imap_unordered(build_cont, HG9):
            K.log(f'cont {s}: {r}', 'p14')


_C = {}


def cont(sta):
    if sta not in _C:
        _C[sta] = (np.load(f'{CONT}/contv_{sta}.npy', mmap_mode='r'), np.load(f'{CONT}/denv_{sta}.npy', mmap_mode='r'))
    return _C[sta]


def cut14(start):
    out = {}
    for ch in CH:
        s, c = ch.rsplit('.', 1); ci = 'NEZ'.index(c)
        X, D = cont(s)
        w = np.array(X[ci, start:start + L], np.float64)
        if not np.all(np.isfinite(w)) or D[ci, start] <= 0:
            out[ch] = None; continue
        w = w - w.mean(); n = np.linalg.norm(w)
        out[ch] = w / n if n > 0 else None
    return out


# ------------------------------------------------------------------------------------------------ bank
def run_bank():
    import hg_common as H
    import hg_r2 as R2
    import cat_e3 as E3
    TT = H.templates() if False else pd.read_csv(f'{K.S13}/tables/hg2018_table1.csv', comment='#')
    TT['time'] = pd.to_datetime(TT.time)
    bank, meta = {}, []
    # H&G T1-T5: V4 windows [origin, origin + 90] in the 2011 hg_repro frame (velocity), off 0
    for k in range(1, 6):
        io = H.tidx(TT[TT.template == k].time.iloc[0])
        T = R2.cut(io, 'vel')
        bank[f'T{k}'] = {ch: (T[ch][0] if T.get(ch) is not None else None) for ch in CH}
        meta.append(dict(name=f'T{k}', kind='hg2011', off=0, members=f'T{k}'))
    # E3 cluster stacks (window origin - 9), off 9
    B = np.load(f'{E3.E3D}/bank.npz')
    EB = pd.read_csv(f'{K.TB}/e3_bank.csv')
    for r in EB[EB.kind == 'cluster'].itertuples():
        bank[r.name] = {ch: (B[r.name][E3.CH.index(ch)] if np.isfinite(B[r.name][E3.CH.index(ch)]).all() else None) for ch in CH}
        meta.append(dict(name=r.name, kind='cluster', off=E3.PRE, members=str(r.members)))
    # H&G T6-T13 (2014, window origin - 9), off 9
    for k in range(6, 14):
        o = tidx(TT[TT.template == k].time.iloc[0])
        bank[f'T{k}'] = cut14(o - 9)
        meta.append(dict(name=f'T{k}', kind='hg2014', off=9, members=f'T{k}', origin_i=o))
    # noise templates: random times in the 2014 ETS window; channels live in >= half of the T6-T13 windows must be live
    need = [ch for ch in CH if sum(bank[f'T{k}'][ch] is not None for k in range(6, 14)) >= 4]
    rng = np.random.default_rng(SEED)
    oi = [tidx(TT[TT.template == k].time.iloc[0]) for k in range(6, 14)]
    got = []
    a, b = tidx('2014-11-14'), tidx('2014-12-15')
    while len(got) < NNOISE:
        i = int(rng.integers(a, b - L))
        if min(abs(i - x) for x in oi) < 180 or any(abs(i - x) < 180 for x in got):
            continue
        w = cut14(i)
        if all(w[ch] is not None for ch in need):
            got.append(i)
            bank[f'N{len(got) - 1:02d}'] = w
            meta.append(dict(name=f'N{len(got) - 1:02d}', kind='noise', off=0, members='', i0=i))
    arr = {nm: np.array([st[ch] if st[ch] is not None else np.full(L, np.nan) for ch in CH]) for nm, st in bank.items()}
    np.savez(f'{P14}/bank.npz', chans=np.array(CH), **arr)
    Mt = pd.DataFrame(meta)
    Mt['n_ch'] = [int(np.isfinite(arr[n]).all(1).sum()) for n in Mt.name]
    Mt.to_csv(f'{K.TB}/p14_bank.csv', index=False)
    print(Mt.to_string()); print('noise channels required:', len(need))


# ------------------------------------------------------------------------------------------------ mf
def draws():
    rng = np.random.default_rng(SEED + 1)
    offs = []
    for _ in range(NDRAW):
        while True:
            o = rng.choice(np.arange(20, NDAY - 20), size=len(HG9), replace=False)
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


def ccser(ch, w):
    from scipy.signal import oaconvolve
    s, c = ch.rsplit('.', 1); ci = 'NEZ'.index(c)
    X, D = cont(s)
    x = np.nan_to_num(np.asarray(X[ci], np.float32))
    num = oaconvolve(x, w[::-1].astype(np.float32), mode='valid')
    d = np.asarray(D[ci]); live = d > 0
    cc = np.zeros(M, np.float32)
    np.divide(num, d, out=cc, where=live); np.clip(cc, -1, 1, out=cc)
    return cc, live


def mf_job(job):
    import cat_e3 as E3
    name, controls = job
    f = f'{MFD}/{name}.npz'
    if os.path.isfile(f):
        return f'{name} exists'
    W = np.load(f'{P14}/bank.npz')[name]
    chs = [ch for k, ch in enumerate(CH) if np.isfinite(W[k]).all()]
    names = ['real'] + ([f'A1_{r}' for r in range(NDRAW)] + [f'A2_{r}' for r in range(NDRAW)] if controls else [])
    S = {n: np.zeros(M, np.float32) for n in names}; NL = {n: np.zeros(M, np.int16) for n in names}
    offs, signs = draws()
    for ch in chs:
        k = CH.index(ch); si = HG9.index(ch.rsplit('.', 1)[0])
        cc, live = ccser(ch, W[k])
        S['real'] += cc; NL['real'] += live
        if controls:
            for r in range(NDRAW):
                sh = int(offs[r, si]) * 86400
                S[f'A1_{r}'] += np.roll(cc, sh); NL[f'A1_{r}'] += np.roll(live, sh)
                S[f'A2_{r}'] += cc * np.float32(signs[r, k]); NL[f'A2_{r}'] += live
    out = {'name': name, 'chans': np.array(chs)}
    for n in names:
        pi, pv, pn, med, mad = E3.peaks(S[n], NL[n])
        out.update({f'{n}_pi': pi, f'{n}_pv': pv, f'{n}_pnl': pn, f'{n}_med': med, f'{n}_mad': mad})
        S[n] = None; NL[n] = None
    np.savez_compressed(f + '.part.npz', **out); os.replace(f + '.part.npz', f)
    return f'{name} done ({len(chs)} ch)'


def run_mf():
    Mt = pd.read_csv(f'{K.TB}/p14_bank.csv')
    jobs = [(r.name, r.kind != 'noise') for r in Mt.itertuples()]
    jobs.sort(key=lambda j: not j[1])
    from multiprocessing import Pool
    with Pool(NW, maxtasksperchild=3) as p:
        for k, r in enumerate(p.imap_unordered(mf_job, jobs)):
            K.log(f'mf {k + 1}/{len(jobs)}: {r}', 'p14')


# ------------------------------------------------------------------------------------------------ detect
def run_detect():
    Mt = pd.read_csv(f'{K.TB}/p14_bank.csv')
    noise = Mt[Mt.kind == 'noise']
    zs = []
    for r in noise.itertuples():
        z = np.load(f'{MFD}/{r.name}.npz')
        pi, pv = z['real_pi'], z['real_pv']
        m = np.abs(pi - int(r.i0)) > 90
        zs.append((pv[m] - float(z['real_med'])) / float(z['real_mad']))
    allz = np.sort(np.concatenate(zs))[::-1]
    TH = {}
    for rate in (0.25, 0.5, 1.0):
        n = int(round(rate * len(noise) * NDAY)); TH[rate] = float(allz[n - 1])
    pd.DataFrame([dict(rate_per_tpl_day=k, z_star_MAD=v) for k, v in TH.items()]).to_csv(f'{K.TB}/p14_noise_threshold.csv', index=False)
    lo = min(LEVEL, TH[1.0])
    rows = []
    for r in Mt.itertuples():
        z = np.load(f'{MFD}/{r.name}.npz')
        med, mad = float(z['real_med']), float(z['real_mad'])
        srcs = ['real'] + ([f'A1_{q}' for q in range(NDRAW)] + [f'A2_{q}' for q in range(NDRAW)] if r.kind != 'noise' else [])
        for src in srcs:
            pi, pv = z[f'{src}_pi'], z[f'{src}_pv']
            zz = (pv - med) / mad
            m = zz >= lo
            if r.kind == 'noise':
                m &= np.abs(pi - int(r.i0)) > 90
            if m.any():
                rows.append(pd.DataFrame(dict(template=r.name, kind=r.kind, source=src if src == 'real' else src[:2],
                                              draw=-1 if src == 'real' else int(src[3:]), i=pi[m] + int(r.off), zmad=zz[m])))
    R = pd.concat(rows, ignore_index=True)
    R['epoch'] = EP0 + R.i
    R['year'] = 2014
    R['period'] = periods(R.epoch.values)
    # agreement: other bank templates (same source/draw) >= LEVEL within +-30 s
    R['agree'] = ''
    bankR = R[R.template.isin(VALID) & (R.zmad >= LEVEL)]
    key = {}
    for (src, dr, t), g in bankR.groupby(['source', 'draw', 'template']):
        key[(src, dr, t)] = np.sort(g.i.values)
    for (src, dr, t), g in R[R.kind != 'noise'].groupby(['source', 'draw', 'template']):
        others = [k for k in key if k[0] == src and k[1] == dr and k[2] != t]
        ag = []
        for i in g.i.values:
            names = [k[2] for k in others if (lambda a: np.any(np.abs(a[np.searchsorted(a, i - 30):np.searchsorted(a, i + 31)] - i) <= 30))(key[k])]
            ag.append(','.join(names))
        R.loc[g.index, 'agree'] = ag
    R.to_parquet(f'{P14}/rows.parquet')
    K.log(f'detect: thresholds {TH}; rows {len(R)} {R.groupby(["kind", "source"]).size().to_dict()}', 'p14')


def run_battery():
    import cat_e2
    cat_e2.run_battery(2014, rows_file=f'{P14}/rows.parquet', out_file=f'{P14}/bat_2014.parquet', tag='p14')


# ------------------------------------------------------------------------------------------------ injection
def bank14(exclude):
    """templates (dict) for the local MF, bank minus templates whose members intersect `exclude`."""
    Mt = pd.read_csv(f'{K.TB}/p14_bank.csv')
    B = np.load(f'{P14}/bank.npz')
    out = []
    for r in Mt[Mt.name.isin(VALID)].itertuples():
        if set(str(r.members).split(',')) & set(exclude):
            continue
        z = np.load(f'{MFD}/{r.name}.npz')
        chs = list(z['chans'])
        out.append(dict(name=r.name, W={ch: B[r.name][CH.index(ch)] for ch in chs}, off=int(r.off),
                        med=float(z['real_med']), mad=float(z['real_mad'])))
    return out


def getter14():
    def get(ch, a, b):
        s, c = ch.rsplit('.', 1)
        X, D = cont(s)
        return np.array(X['NEZ'.index(c), a:b], np.float64)
    return get


def run_inj_list():
    import cat_inj as I
    R = pd.read_parquet(f'{P14}/rows.parquet')
    det = np.sort(R[(R.source == 'real') & (R.kind != 'noise') & (R.zmad >= LEVEL)].i.unique())
    Z = np.load(f'{K.S13}/data/scan_2014/masks2_2014.npz')
    mk = Z['bw5'] | Z['reg'] | Z['reg3']
    rng = np.random.default_rng(SEED + 2)
    rows, used = [], []
    for src in ('C0', 'C1'):
        for mw in (3.4, 3.6, 3.8, 4.0, 4.2):
            got, tries = 0, 0
            while got < 100 and tries < 50000:
                tries += 1
                i = int(rng.integers(86400, N - 86400))
                ep = EP0 + i
                if periods([ep])[0] != 'between':
                    continue
                nl = 0
                for ch in CH:
                    s, c = ch.rsplit('.', 1)
                    X, D = cont(s)
                    nl += np.isfinite(X['NEZ'.index(c), i - 400:i + 400]).all()
                if nl < 0.5 * len(CH) - 1:
                    continue
                if mk[max(i - 60, 0):i + 121].any() or K.head2_maxmag([ep], [2014])[0] >= 5.5:
                    continue
                k = np.searchsorted(det, i)
                if (k < len(det) and det[k] - i < 600) or (k > 0 and i - det[k - 1] < 600):
                    continue
                if any(abs(i - j) < 600 for j in used):
                    continue
                used.append(i); got += 1
                rows.append(dict(source=src, Mw=mw, i=i, epoch=ep, n_live_ch=int(nl)))
    Lst = pd.DataFrame(rows); Lst['inj_id'] = np.arange(len(Lst))
    Lst.to_parquet(f'{P14}/inj_list.parquet')
    print(Lst.groupby(['source', 'Mw']).size())


_P = {}


def inj_init():
    import cat_inj as I
    I.bat_init(early=False)
    W = I._B['SC'].W
    I._B['base_ref_all'] = {r['sta']: r['X'] for r in W['refs']}
    S, SM = I.load_sources()
    _P['I'], _P['S'], _P['SM'] = I, S, SM
    _P['bank'] = {src: bank14(str(SM.loc[src, 'members']).split(',') + [src]) for src in ('C0', 'C1')}
    _P['get'] = getter14()


def inj_row(r):
    I, S, SM = _P['I'], _P['S'], _P['SM']
    src = r['source']; sm = SM.loc[src]
    sc = 10 ** (1.5 * (r['Mw'] - sm.Mw_ref))
    res = I.mf_local(_P['bank'][src], _P['get'], int(r['i']), inj=S[src]['mf'], scale=sc)
    zs = sorted(((v[0], k) for k, v in res.items() if np.isfinite(v[0])), reverse=True)
    o = dict(inj_id=int(r['inj_id']), mf_best_z=zs[0][0] if zs else np.nan, mf_best_tpl=zs[0][1] if zs else '',
             mf_second_z=zs[1][0] if len(zs) > 1 else np.nan)
    tsec = int(round(r['epoch'] - I._B['X'].EV['y0']))
    far = I.far_synth(sm, r['Mw'])
    I.set_injection(tsec, sc_sig=S[src]['sc'], ref_sig=far, scale=sc)
    for rr in I._B['SC'].W['refs']:
        if isinstance(rr['X'], I.Proxy):
            rr['X'].inj = [(tsec + I.W0, far[rr['sta']])]
    b = I.battery_local(r['epoch']); b['epoch'] = r['epoch']
    I.set_injection(tsec)
    p = I.score(b, 2014)
    o.update(p)
    for k in ('z_best', 'ratio_far', 'n_far', 'z_loso', 'lat_best', 'lon_best'):
        o[k] = b.get(k, np.nan)
    return o


def inj_chunk(rows):
    out = []
    for r in rows:
        try:
            out.append(inj_row(r))
        except Exception as e:  # noqa
            out.append(dict(inj_id=int(r['inj_id']), err=f'{type(e).__name__}: {str(e)[:120]}'))
    return out


def run_inj():
    assert os.environ.get('SCAN_YEAR') == '2014'
    Lst = pd.read_parquet(f'{P14}/inj_list.parquet')
    rows = Lst.to_dict('records')
    from multiprocessing import Pool
    nch = NW * 6
    chunks = [rows[k::nch] for k in range(nch)]
    out = []
    with Pool(NW, initializer=inj_init) as p:
        for k, r in enumerate(p.imap_unordered(inj_chunk, chunks)):
            out += r
            if k % 12 == 0:
                K.log(f'inj chunk {k + 1}/{nch}', 'p14')
    pd.DataFrame(out).to_parquet(f'{P14}/inj.parquet')
    K.log(f'inj done {len(out)}', 'p14')


# ------------------------------------------------------------------------------------------------ truth / summary
def dedupe(i):
    t = np.sort(np.asarray(i))
    return int(1 + np.sum(np.diff(t) > 30)) if len(t) else 0


def run_summary():
    import cat_e2
    R = pd.read_parquet(f'{P14}/rows.parquet')
    B = pd.read_parquet(f'{P14}/bat_2014.parquet').drop(columns=['epoch'])
    B['grid_edge'] = B.grid_edge.where(B.grid_edge.notna(), B.grid_edge_node) if 'grid_edge' in B else B.grid_edge_node
    D = R.merge(B, on='i', how='left')
    for c in ('n_far', 'n6'):
        D[c] = D[c].fillna(0)
    D = K.apply_rules(D)
    D['lfe_status'], D['n_lfe_10m'] = cat_e2.lfe_status(D.epoch.values)
    TH = pd.read_csv(f'{K.TB}/p14_noise_threshold.csv').set_index('rate_per_tpl_day').z_star_MAD
    D['tier'] = D.template.isin(VALID) & (D.zmad >= LEVEL) & (D.agree != '') & D.battery_e2
    D['tier_noise05'] = D.template.isin(VALID) & (D.zmad >= TH[0.5]) & D.battery_e2
    D['vr_thrust_ge012'] = D.vr_thrust >= 0.12
    D.to_parquet(f'{P14}/p14_rows_scored.parquet')
    TT = pd.read_csv(f'{K.S13}/tables/hg2018_table1.csv', comment='#')
    # truth
    tr = []
    for k in range(6, 14):
        o = tidx(TT[TT.template == k].time.iloc[0])
        x = D[(D.source == 'real') & D.template.isin(VALID) & ((D.i - o).abs() <= 30)]
        ag = x.agree.apply(lambda s: [a for a in s.split(',') if a and a != f'T{k}'])
        tier = x[(x.zmad >= LEVEL) & (ag.str.len() > 0) & x.battery_e2]
        mfonly = x[(x.zmad >= LEVEL) & (ag.str.len() > 0)]
        tr.append(dict(template=f'T{k}', time_pt=TT[TT.template == k].time_PT.iloc[0], Mw=TT[TT.template == k].Mw.iloc[0],
                       best_other_tpl=x.sort_values('zmad').template.iloc[-1] if len(x) else '',
                       best_other_z=float(x.zmad.max()) if len(x) else np.nan,
                       refound_tier=len(tier) > 0, refound_mf_2tpl=len(mfonly) > 0,
                       refound_noise05=bool(((x.zmad >= TH[0.5]) & x.battery_e2).any()),
                       battery_e2_at_time=bool(x.battery_e2.any()) if len(x) else False,
                       z_best=float(x.z_best.max()) if len(x) else np.nan,
                       ratio_far=float(x.ratio_far.min()) if len(x) else np.nan))
    TR = pd.DataFrame(tr); TR.to_csv(f'{K.TB}/p14_truth.csv', index=False)
    # counts
    yb = between_years()
    rows = []
    for tiercol in ('tier', 'tier_noise05'):
        x = D[D[tiercol]]
        r = dict(tier=tiercol, between_years=yb)
        for per in ('ETS', 'between', 'other_episode'):
            y = x[x.period == per]
            r[f'real_{per}'] = dedupe(y[y.source == 'real'].i)
            cc = [dedupe(y[(y.source == a) & (y.draw == q)].i) for a in ('A1', 'A2') for q in range(NDRAW)]
            r[f'ctrl_mean_{per}'] = float(np.mean(cc)); r[f'ctrl_max_{per}'] = int(np.max(cc))
        r['ctrl_rate_between_per_yr'] = r['ctrl_mean_between'] / yb
        rows.append(r)
    C = pd.DataFrame(rows); C.to_csv(f'{K.TB}/p14_counts.csv', index=False)
    # completeness
    Lst = pd.read_parquet(f'{P14}/inj_list.parquet')
    J = Lst.merge(pd.read_parquet(f'{P14}/inj.parquet'), on='inj_id', how='left')
    J['rec_tier'] = (J.mf_best_z >= LEVEL) & (J.mf_second_z >= LEVEL) & J.battery_e2.fillna(False).astype(bool)
    J['rec_mf2'] = (J.mf_best_z >= LEVEL) & (J.mf_second_z >= LEVEL)
    J['rec_noise05'] = (J.mf_best_z >= TH[0.5]) & J.battery_e2.fillna(False).astype(bool)
    J.to_parquet(f'{P14}/inj_scored.parquet')
    G = J.groupby('Mw').agg(n=('inj_id', 'count'), rec_tier=('rec_tier', 'mean'), rec_mf2=('rec_mf2', 'mean'),
                            battery_e2=('battery_e2', 'mean'), rec_noise05=('rec_noise05', 'mean'),
                            mf_best_z_median=('mf_best_z', 'median'), mf_second_z_median=('mf_second_z', 'median'),
                            z_best_median=('z_best', 'median')).reset_index()
    G.to_csv(f'{K.TB}/p14_completeness.csv', index=False)
    c = C[C.tier == 'tier'].iloc[0]
    comp38 = float(G[G.Mw == 3.8].rec_tier.iloc[0])
    V = dict(truth_refound=int(TR.refound_tier.sum()), truth_n=len(TR), between_real=int(c.real_between),
             between_ctrl_mean=float(c.ctrl_mean_between), ctrl_rate_per_yr=float(c.ctrl_rate_between_per_yr),
             completeness_Mw38=comp38)
    V['go_items'] = dict(truth=V['truth_refound'] >= 5, ctrl=V['ctrl_rate_per_yr'] <= 1, completeness=comp38 >= 0.8)
    V['between_is_reported_not_gated'] = True                                   # Amendment Q0
    V['GO'] = bool(all(V['go_items'].values()))
    json.dump(V, open(f'{K.TB}/p14_verdict.json', 'w'), indent=1)
    pd.set_option('display.width', 250)
    print(TR.to_string()); print(C.round(3).to_string()); print(G.round(3).to_string()); print(V)
    figs(D, G)


def figs(D, G):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8), sharey=True)
    for ax, per in zip(axs, ('ETS', 'between', 'other_episode')):
        x = D[(D.period == per) & (D.kind != 'noise')]
        for src, c, s in (('A', 'tab:purple', 4), ('real', 'tab:red', 10)):
            y = x[x.source.str.startswith(src)]
            ax.scatter(y.zmad, y.z_best.clip(0.5, 300), s=s, c=c, lw=0, alpha=0.6, label=f'{src} {len(y)}')
            yb = y[y.tier]
            ax.scatter(yb.zmad, yb.z_best, s=50, facecolors='none', edgecolors=c, lw=1.2)
        ax.set_yscale('log'); ax.axhline(5, color='k', lw=0.5); ax.axvline(LEVEL, color='k', lw=0.5)
        ax.set_title(f'2014 {per}: rings = tier (5.5 MAD + 2 templates + battery′)', fontsize=9); ax.legend(fontsize=7)
        ax.set_xlabel('MF z (MAD)')
    axs[0].set_ylabel('scan z_best')
    fig.tight_layout(); fig.savefig(f'{K.FIG}/p14_mf_vs_scanz.png', dpi=120)
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(G.Mw, G.rec_tier, 'o-', color='tab:red', label='tier (MF 2 templates + battery′)')
    ax.plot(G.Mw, G.rec_mf2, 's--', color='0.4', label='MF 2 templates only')
    ax.plot(G.Mw, G.battery_e2, '^:', color='tab:blue', label='battery′ only')
    ax.plot(G.Mw, G.rec_noise05, 'd-.', color='tab:green', label='noise-rate 0.5/day + battery′')
    ax.axhline(0.8, color='k', lw=0.5); ax.set_ylim(0, 1.02); ax.set_xlabel('injected Mw'); ax.set_ylabel('fraction recovered')
    ax.set_title('Pilot 2014, H&G 9 stations, cross-template injection (C0, C1)', fontsize=9); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f'{K.FIG}/p14_completeness.png', dpi=130)


if __name__ == '__main__':
    m = sys.argv[1]
    {'cont': run_cont, 'bank': run_bank, 'mf': run_mf, 'detect': run_detect, 'battery': run_battery,
     'inj_list': run_inj_list, 'inj': run_inj, 'summary': run_summary}[m]()
