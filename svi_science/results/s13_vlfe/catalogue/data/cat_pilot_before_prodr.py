#!/usr/bin/env python
"""Pilots 2013-09 / 2018-06 (PLAN_CAT Amendment Q3/Q4). Environment: PILOT=2013|2018 (sets frame, stations, ETS).
  cont | bank | mf | detect | battery (SCAN_YEAR=PILOT) | inj_list | inj (SCAN_YEAR=PILOT) | summary
Bank = validated templates only: H&G T1-T5 re-cut (window origin - 9 s) + E3 clusters C0-C5 (E3 alignment), cut from
2011-12 data on the pilot's station set; statistic plain or weighted (Q1 verdict, tables/q1_verdict.json).
Tier = MF >= 5.5 MAD + two-template agreement + battery'. Noise templates: 25 random ETS-window times."""
import os, sys, json
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

PROD = os.environ.get('PROD_YEAR')                     # Amendment S production mode (pilot behaviour unchanged without it)
PILOT = int(PROD or os.environ.get('PILOT', '2013'))
NAME = 'prod' if PROD else 'pilot'
CFG = {2013: dict(stations=['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.PFB', 'CN.LZB', 'UW.LRIV', 'UW.DOSE', 'UW.GNW',
                               'UW.SP2', 'TA.A04D', 'TA.D03D', 'UW.WISH', 'US.NLWA'],
                  setname='dense2013', ets=('2013-09-10', '2013-10-16')),   # Q1': dense2011 stations with >= 180 bb days in 2013
       2008: dict(stations=['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'UW.GNW'], setname='era2007',
                  ets=('2008-05-08', '2008-06-05')),                  # SVI 2008-05 (Bostock list 05-10..06-02) +- 2 d
       2018: dict(stations=['CN.PGC', 'CN.VGZ', 'UW.LRIV', 'CN.SNB', 'CN.PFB', 'CN.NLLB', 'UW.DOSE', 'UW.GNW', 'UW.SP2', 'UW.WISH'],
                  setname='era2016', ets=('2018-06-11', '2018-07-20'))}.get(PILOT)   # production years: see PROD below
if PROD:
    def _prod_cfg(year):
        E = pd.read_csv(f'{K.ROOT}/svi_science/results/b_ets_record/episodes_b.csv')
        E = E[(E.region == 'SVI') & (pd.to_datetime(E.end) >= pd.Timestamp(f'{year}-01-01'))
              & (pd.to_datetime(E.start) < pd.Timestamp(f'{year + 1}-01-01'))]
        wins = [(str((pd.Timestamp(r.start) - pd.Timedelta(days=2)).date()), str((pd.Timestamp(r.end) + pd.Timedelta(days=3)).date()), r.episode)
                for r in E.itertuples()]
        T = pd.read_csv(f'{K.TB}/prod_stations.csv')
        st = T[T.year == year].sta.tolist()
        majors = [r.episode for r in E.itertuples() if bool(r.major)]
        return dict(stations=st, setname=f'prod{year}', ets=wins[0][:2] if wins else (f'{year}-01-01', f'{year}-03-01'), wins=wins,
                    majors=majors)
    CFG = _prod_cfg(PILOT)
STAS = CFG['stations']
CH = [f'{s}.{c}' for s in STAS for c in 'NEZ']
PD = f'{K.DATA}/{NAME}{PILOT}'
CONT, MFD = f'{PD}/cont', f'{PD}/mf'
for d in (PD, CONT, MFD):
    os.makedirs(d, exist_ok=True)
T0 = pd.Timestamp(f'{PILOT}-01-01'); NDAY = 366 if PILOT % 4 == 0 else 365; N = NDAY * 86400; L = 90; M = N - L + 1
EP0 = (T0 - pd.Timestamp('1970-01-01')).total_seconds()
NDRAW, NNOISE, NW, LEVEL, PRE = 8, 25, 8, 5.5, 9
SFX = os.environ.get('PILOT_SFX', '')                # e.g. '_v2' for the scan-v2 battery (Amendment R)
NINJ = int(os.environ.get('NINJ', '50'))           # injections per source x Mw x period (Amendment R: 400)
SEED = 20261015 + PILOT
VALID = [f'T{k}' for k in range(1, 6)] + [f'C{k}' for k in range(6)]


def stat():
    f = f'{K.TB}/q1_verdict.json'
    return json.load(open(f)).get('stat_for_later', 'plain') if os.path.isfile(f) else 'plain'


def tidx(t):
    return int(round((pd.Timestamp(t) - T0).total_seconds()))


def periods(epoch):
    t = pd.to_datetime(np.asarray(epoch, float), unit='s')
    out = np.full(len(t), 'between', dtype=object)
    E = pd.read_csv(f'{K.ROOT}/svi_science/results/b_ets_record/episodes_b.csv')
    for r in E[E.region.isin(['SVI', 'NW'])].itertuples():
        a = pd.Timestamp(r.start) - pd.Timedelta(days=2); b = pd.Timestamp(r.end) + pd.Timedelta(days=3)
        out[(t >= a) & (t < b)] = 'other_episode'
    for a_, b_ in ([w[:2] for w in CFG['wins']] if CFG.get('wins') else [CFG['ets']]):
        out[(t >= pd.Timestamp(a_)) & (t < pd.Timestamp(b_))] = 'ETS'
    return out


def period_years():
    days = pd.date_range(T0, periods=NDAY, freq='D') + pd.Timedelta(hours=12)
    p = periods((days - pd.Timestamp('1970-01-01')).total_seconds().values)
    return {k: float((p == k).sum() / 365.25) for k in ('ETS', 'between', 'other_episode')}


# ------------------------------------------------------------------------------------------------ cont (frame)
def build_cont(sta):
    import hg_cont
    import cat_p14 as P
    f = f'{CONT}/contv_{sta}.npy'
    if os.path.isfile(f) and os.path.isfile(f'{CONT}/denv_{sta}.npy'):
        return sta, 'exists'
    X = np.full((3, N), np.nan, np.float32); nd = 0
    for k in range(NDAY):
        d = T0 + pd.Timedelta(days=k)
        p = f'{K.S13}/data/disp1sps_v3/{sta}/{d.year}/{d.dayofyear:03d}.npz'
        if os.path.isfile(p):
            X[:, k * 86400:(k + 1) * 86400] = np.load(p)['d']; nd += 1
    for c in range(3):
        x = hg_cont.filt(X[c]).astype(np.float64)
        v = np.full(x.shape, np.nan); v[1:-1] = 0.5 * (x[2:] - x[:-2]); X[c] = v.astype(np.float32)
    np.save(f + '.part.npy', X); os.replace(f + '.part.npy', f)
    D = dens_local(X)
    np.save(f'{CONT}/denv_{sta}.part.npy', D); os.replace(f'{CONT}/denv_{sta}.part.npy', f'{CONT}/denv_{sta}.npy')
    return sta, f'days {nd}/{NDAY} live {(D > 0).mean(axis=1).round(3)}'


def dens_local(X):
    from hg_cont import den_chunked
    import hg_common as H
    Mh = H.M; H.M = M
    try:
        D = np.zeros((3, M), np.float32)
        for c in range(3):
            den = den_chunked(X[c]); med = np.median(den[den > 0]) if np.any(den > 0) else 0
            den[den < 1e-4 * med] = 0.0; D[c] = den.astype(np.float32)
    finally:
        H.M = Mh
    return D


def run_cont():
    from multiprocessing import Pool
    with Pool(5) as p:
        for s, r in p.imap_unordered(build_cont, STAS):
            K.log(f'{NAME}{PILOT} cont {s}: {r}', 'pilot')


_C = {}


def cont(sta):
    if sta not in _C:
        _C[sta] = (np.load(f'{CONT}/contv_{sta}.npy', mmap_mode='r'), np.load(f'{CONT}/denv_{sta}.npy', mmap_mode='r'))
    return _C[sta]


def cut_frame(start):
    out = {}
    for ch in CH:
        s, c = ch.rsplit('.', 1); ci = 'NEZ'.index(c)
        X, D = cont(s)
        w = np.array(X[ci, start:start + L], np.float64)
        if not np.all(np.isfinite(w)) or D[ci, start] <= 0:
            out[ch] = None; continue
        w = w - w.mean(); n = np.linalg.norm(w); out[ch] = w / n if n > 0 else None
    return out


def rms_ratio_frame(start):
    out = {}
    for ch in CH:
        s, c = ch.rsplit('.', 1); ci = 'NEZ'.index(c)
        X, D = cont(s)
        w = np.array(X[ci, start:start + L], np.float64); pre = np.array(X[ci, start - 600:start], np.float64)
        ok = np.isfinite(w).all() and np.isfinite(pre).sum() > 300 and np.nanmean(pre ** 2) > 0
        out[ch] = float(np.sqrt(np.mean(w ** 2)) / np.sqrt(np.nanmean(pre ** 2))) if ok else np.nan
    return out


# ------------------------------------------------------------------------------------------------ bank
def run_bank():
    import cat_q1 as Q                                       # E3 machinery on the extended 2011-12 station list
    import cat_e3 as E3
    import hg_common as H
    V = pd.read_csv(f'{K.TB}/e3_events.csv')
    TT = H.templates().set_index('template')
    bank, wts, meta = {}, {}, []

    def add(name, starts, kind, members):
        Wn = [E3.cut(int(s)) for s in starts]
        st = E3.stack(Wn) if len(Wn) > 1 else Wn[0]
        rr = [Q.rms_ratio(int(s)) for s in starts]
        bank[name] = {ch: st.get(ch) for ch in CH}
        wts[name] = {ch: (float(np.nanmean([r[ch] for r in rr])) if np.isfinite([r[ch] for r in rr]).any() else np.nan) for ch in CH}
        meta.append(dict(name=name, kind=kind, members=members, off=PRE))
    for k in range(1, 6):
        add(f'T{k}', [H.tidx(TT.loc[k].time) - PRE], 'hg2011', f'T{k}')
    for c, g in V.groupby('cluster'):
        add(f'C{c}', (g.start0 + g.shift_s).astype(int).values, 'cluster', ','.join(map(str, g.cand_id)))
    rng = np.random.default_rng(SEED)
    a, b = tidx(CFG['ets'][0]), tidx(CFG['ets'][1])
    if PROD:
        wins = [(max(tidx(a_), 700), min(tidx(b_), N - 700)) for a_, b_ in [w[:2] for w in CFG['wins']]]
        wins = [w for w in wins if w[1] - w[0] > 3600] or [(700, 60 * 86400)]
    need = [ch for ch in CH if sum(bank[t][ch] is not None for t in VALID) >= 6]
    got = []
    while len(got) < NNOISE:
        if PROD:
            a, b = wins[int(rng.integers(0, len(wins)))]
        i = int(rng.integers(a + 600, b - L))
        if any(abs(i - x) < 180 for x in got):
            continue
        w = cut_frame(i)
        if all(w[ch] is not None for ch in need):
            got.append(i); nm = f'N{len(got) - 1:02d}'
            bank[nm] = w; wts[nm] = rms_ratio_frame(i)
            meta.append(dict(name=nm, kind='noise', members='', off=0, i0=i))
    arr = {nm: np.array([st[ch] if st[ch] is not None else np.full(L, np.nan) for ch in CH]) for nm, st in bank.items()}
    warr = {f'w__{nm}': np.array([wts[nm][ch] for ch in CH]) for nm in bank}
    np.savez(f'{PD}/bank.npz', chans=np.array(CH), **arr, **warr)
    Mt = pd.DataFrame(meta); Mt['n_ch'] = [int(np.isfinite(arr[n]).all(1).sum()) for n in Mt.name]
    Mt.to_csv(f'{K.TB}/{NAME}{PILOT}_bank.csv', index=False)
    print(Mt.to_string()); print('noise channels required', len(need), 'statistic', stat())


# ------------------------------------------------------------------------------------------------ mf
def draws():
    rng = np.random.default_rng(SEED + 1)
    offs, signs = [], []
    for _ in range(NDRAW):
        while True:
            o = rng.choice(np.arange(20, NDAY - 20), size=len(STAS), replace=False)
            if np.min(np.diff(np.sort(o))) >= 2:
                break
        offs.append(o.astype(int))
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
    cc = np.zeros(M, np.float32); np.divide(num, d, out=cc, where=live); np.clip(cc, -1, 1, out=cc)
    return cc, live


def chan_weights(name, chs):
    B = np.load(f'{PD}/bank.npz')
    wv = B[f'w__{name}']
    wk = np.array([wv[CH.index(ch)] for ch in chs])
    if stat() == 'plain' or not np.isfinite(wk).all() or wk.sum() <= 0:
        return np.ones(len(chs))
    return wk * len(wk) / wk.sum()


def usable_chans(name):
    B = np.load(f'{PD}/bank.npz')
    W = B[name]; wv = B[f'w__{name}']
    return [ch for k, ch in enumerate(CH) if np.isfinite(W[k]).all() and (stat() == 'plain' or np.isfinite(wv[k]))]


def mf_job(job):
    import cat_e3 as E3
    name, controls = job
    f = f'{MFD}/{name}.npz'
    if os.path.isfile(f):
        return f'{name} exists'
    W = np.load(f'{PD}/bank.npz')[name]
    chs = usable_chans(name); wk = chan_weights(name, chs)
    names = ['real'] + ([f'A1_{r}' for r in range(NDRAW)] + [f'A2_{r}' for r in range(NDRAW)] if controls else [])
    S = {n: np.zeros(M, np.float32) for n in names}; NL = {n: np.zeros(M, np.int16) for n in names}
    offs, signs = draws()
    for j, ch in enumerate(chs):
        k = CH.index(ch); si = STAS.index(ch.rsplit('.', 1)[0]); ww = np.float32(wk[j])
        cc, live = ccser(ch, W[k])
        S['real'] += ww * cc; NL['real'] += live
        if controls:
            for r in range(NDRAW):
                sh = int(offs[r, si]) * 86400
                S[f'A1_{r}'] += ww * np.roll(cc, sh); NL[f'A1_{r}'] += np.roll(live, sh)
                S[f'A2_{r}'] += (ww * np.float32(signs[r, k])) * cc; NL[f'A2_{r}'] += live
    out = {'name': name, 'chans': np.array(chs), 'weights': wk}
    for n in names:
        pi, pv, pn, med, mad = E3.peaks(S[n], NL[n])
        out.update({f'{n}_pi': pi, f'{n}_pv': pv, f'{n}_med': med, f'{n}_mad': mad}); S[n] = None; NL[n] = None
    np.savez_compressed(f + '.part.npz', **out); os.replace(f + '.part.npz', f)
    return f'{name} done ({len(chs)} ch, {stat()})'


def run_mf():
    Mt = pd.read_csv(f'{K.TB}/{NAME}{PILOT}_bank.csv')
    jobs = sorted([(r.name, r.kind != 'noise') for r in Mt.itertuples()], key=lambda j: not j[1])
    from multiprocessing import Pool
    with Pool(NW, maxtasksperchild=3) as p:
        for k, r in enumerate(p.imap_unordered(mf_job, jobs)):
            K.log(f'{NAME}{PILOT} mf {k + 1}/{len(jobs)}: {r}', 'pilot')


# ------------------------------------------------------------------------------------------------ detect
def agree_col(R, names, key_fn):
    R['agree'] = ''
    for grp, g in R.groupby(['source', 'draw', 'template']):
        others = key_fn(grp)
        ag = []
        for i in g.i.values:
            ag.append(','.join(t for t, a in others if np.any(np.abs(a[np.searchsorted(a, i - 30):np.searchsorted(a, i + 31)] - i) <= 30)))
        R.loc[g.index, 'agree'] = ag
    return R


def run_detect():
    Mt = pd.read_csv(f'{K.TB}/{NAME}{PILOT}_bank.csv')
    noise = Mt[Mt.kind == 'noise']
    zs = []
    for r in noise.itertuples():
        z = np.load(f'{MFD}/{r.name}.npz')
        m = np.abs(z['real_pi'] - int(r.i0)) > 90
        zs.append((z['real_pv'][m] - float(z['real_med'])) / float(z['real_mad']))
    allz = np.sort(np.concatenate(zs))[::-1]
    TH = {rate: float(allz[int(round(rate * len(noise) * NDAY)) - 1]) for rate in (0.25, 0.5, 1.0)}
    pd.DataFrame([dict(rate=k, z_star=v) for k, v in TH.items()]).to_csv(f'{K.TB}/{NAME}{PILOT}_noise_threshold.csv', index=False)
    rows = []
    for r in Mt.itertuples():
        z = np.load(f'{MFD}/{r.name}.npz'); med, mad = float(z['real_med']), float(z['real_mad'])
        srcs = ['real'] + ([f'A1_{q}' for q in range(NDRAW)] + [f'A2_{q}' for q in range(NDRAW)] if r.kind != 'noise' else [])
        for src in srcs:
            zz = (z[f'{src}_pv'] - med) / mad; pi = z[f'{src}_pi']
            m = zz >= min(LEVEL, TH[1.0])
            if r.kind == 'noise':
                m &= np.abs(pi - int(r.i0)) > 90
            if m.any():
                rows.append(pd.DataFrame(dict(template=r.name, kind=r.kind, source=src if src == 'real' else src[:2],
                                              draw=-1 if src == 'real' else int(src[3:]), i=pi[m] + int(r.off), zmad=zz[m])))
    R = pd.concat(rows, ignore_index=True)
    R['epoch'] = EP0 + R.i; R['year'] = PILOT; R['period'] = periods(R.epoch.values)
    key = {}
    for (src, dr, t), g in R[(R.zmad >= LEVEL)].groupby(['source', 'draw', 'template']):
        key[(src, dr, t)] = np.sort(g.i.values)
    nn = sorted(noise.name)

    def others(grp):
        src, dr, t = grp
        if t in VALID:
            return [(k[2], a) for k, a in key.items() if k[0] == src and k[1] == dr and k[2] in VALID and k[2] != t]
        if t in nn:                                           # noise pair: noise j with noise j+1
            j = nn.index(t); t2 = nn[(j + 1) % len(nn)]
            return [(t2, key[('real', -1, t2)])] if ('real', -1, t2) in key else []
        return []
    R = agree_col(R, VALID, others)
    R.to_parquet(f'{PD}/rows.parquet')
    K.log(f'{NAME}{PILOT} detect: thresholds {TH}; rows {len(R)}', 'pilot')


def run_battery():
    import cat_e2
    cat_e2.run_battery(PILOT, rows_file=f'{PD}/rows.parquet', out_file=f'{PD}/bat{SFX}.parquet', tag=f'{NAME}{PILOT}')


# ------------------------------------------------------------------------------------------------ injections
def bank_local(exclude):
    Mt = pd.read_csv(f'{K.TB}/{NAME}{PILOT}_bank.csv'); B = np.load(f'{PD}/bank.npz')
    out = []
    for r in Mt[Mt.name.isin(VALID)].itertuples():
        if set(str(r.members).split(',')) & set(exclude):
            continue
        z = np.load(f'{MFD}/{r.name}.npz'); chs = list(z['chans']); wk = z['weights']
        out.append(dict(name=r.name, W={ch: wk[j] * B[r.name][CH.index(ch)] for j, ch in enumerate(chs)}, off=int(r.off),
                        med=float(z['real_med']), mad=float(z['real_mad'])))
    return out


def getter():
    def get(ch, a, b):
        s, c = ch.rsplit('.', 1); X, D = cont(s)
        return np.array(X['NEZ'.index(c), a:b], np.float64)
    return get


def run_inj_list():
    R = pd.read_parquet(f'{PD}/rows.parquet')
    det = np.sort(R[(R.source == 'real') & R.template.isin(VALID) & (R.zmad >= LEVEL)].i.unique())
    SDY = f'{K.S13}/data/scan' if PILOT == 2011 else f'{K.S13}/data/scan_{PILOT}'
    Z = np.load(f'{SDY}/masks2_{PILOT}.npz'); mk = Z['bw5'] | Z['reg'] | Z['reg3']
    rng = np.random.default_rng(SEED + 2)
    rows, used = [], []
    for per in ('ETS', 'between'):
        for src in ('C0', 'C1'):
            for mw in (3.4, 3.6, 3.8, 4.0):
                got, tries = 0, 0
                if NINJ > 50:
                    used = []                                 # Amendment R: >= 10 min spacing within each (source, Mw) set
                while got < NINJ and tries < 400000:
                    tries += 1
                    i = int(rng.integers(86400, N - 86400)); ep = EP0 + i
                    if periods([ep])[0] != per:
                        continue
                    nl = sum(np.isfinite(cont(ch.rsplit('.', 1)[0])[0]['NEZ'.index(ch[-1]), i - 400:i + 400]).all() for ch in CH)
                    if nl < 0.5 * len(CH) - 1 or mk[max(i - 60, 0):i + 121].any() or K.head2_maxmag([ep], [PILOT])[0] >= 5.5:
                        continue
                    k = np.searchsorted(det, i)
                    if (k < len(det) and det[k] - i < 600) or (k > 0 and i - det[k - 1] < 600) or any(abs(i - j) < 600 for j in used):
                        continue
                    used.append(i); got += 1
                    rows.append(dict(period=per, source=src, Mw=mw, i=i, epoch=ep, n_live_ch=int(nl)))
    Lst = pd.DataFrame(rows); Lst['inj_id'] = np.arange(len(Lst)); Lst.to_parquet(f'{PD}/inj_list.parquet')
    print(Lst.groupby(['period', 'source', 'Mw']).size())


_P = {}


def inj_init():
    import cat_inj as I
    I.bat_init(early=False)
    S, SM = I.load_sources()
    _P.update(I=I, S=S, SM=SM, get=getter(),
              bank={src: bank_local(str(SM.loc[src, 'members']).split(',') + [src]) for src in ('C0', 'C1')})


def inj_row(r):
    I, S, SM = _P['I'], _P['S'], _P['SM']
    src = r['source']; sm = SM.loc[src]; sc = 10 ** (1.5 * (r['Mw'] - sm.Mw_ref))
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
    o.update(I.score(b, PILOT))
    for k in ('z_best', 'ratio_far', 'n_far', 'z_loso', 'vr_thrust'):
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
    assert os.environ.get('SCAN_YEAR') == str(PILOT)
    rows = pd.read_parquet(f'{PD}/inj_list.parquet').to_dict('records')
    from multiprocessing import Pool
    nch = NW * 6; out = []
    with Pool(NW, initializer=inj_init) as p:
        for r in p.imap_unordered(inj_chunk, [rows[k::nch] for k in range(nch)]):
            out += r
    pd.DataFrame(out).to_parquet(f'{PD}/inj{SFX}.parquet')
    K.log(f'{NAME}{PILOT} inj done {len(out)}', 'pilot')


# ------------------------------------------------------------------------------------------------ summary
def dedupe_events(x):
    t = np.sort(np.asarray(x.i))
    return int(1 + np.sum(np.diff(t) > 30)) if len(t) else 0


def run_summary():
    import cat_e2
    from scipy.stats import poisson
    R = pd.read_parquet(f'{PD}/rows.parquet')
    B = pd.read_parquet(f'{PD}/bat{SFX}.parquet').drop(columns=['epoch'])
    B['grid_edge'] = B.grid_edge.where(B.grid_edge.notna(), B.grid_edge_node) if 'grid_edge' in B else B.grid_edge_node
    D = R.merge(B, on='i', how='left')
    for c in ('n_far', 'n6'):
        D[c] = D[c].fillna(0)
    D = K.apply_rules(D)
    D['lfe_status'], D['n_lfe_10m'] = cat_e2.lfe_status(D.epoch.values)
    TH = pd.read_csv(f'{K.TB}/{NAME}{PILOT}_noise_threshold.csv').set_index('rate').z_star
    D['tier'] = D.template.isin(VALID) & (D.zmad >= LEVEL) & (D.agree != '') & D.battery_e2
    D['noise_pair_tier'] = (D.kind == 'noise') & (D.zmad >= LEVEL) & (D.agree != '') & D.battery_e2
    D['tier_noise05'] = D.template.isin(VALID) & (D.zmad >= TH[0.5]) & D.battery_e2
    D['vr_thrust_ge012'] = D.vr_thrust >= 0.12
    D.to_parquet(f'{PD}/rows_scored{SFX}.parquet')
    PY = period_years()
    C = {}
    for per in ('ETS', 'between', 'other_episode'):
        x = D[D.tier & (D.period == per)]
        C[f'real_{per}'] = dedupe_events(x[x.source == 'real'])
        cc = [dedupe_events(x[(x.source == a) & (x.draw == q)]) for a in ('A1', 'A2') for q in range(NDRAW)]
        C[f'ctrl_mean_{per}'] = float(np.mean(cc)); C[f'ctrl_max_{per}'] = int(np.max(cc))
        xn = D[D.noise_pair_tier & (D.period == per)]
        C[f'noise_pair_{per}'] = float(np.mean([dedupe_events(xn[xn.template == t]) for t in sorted(D[D.kind == 'noise'].template.unique())] or [0]))
    C.update({f'years_{k}': v for k, v in PY.items()})
    C['ctrl_rate_between_per_yr'] = C['ctrl_mean_between'] / PY['between']
    C['noise_pair_rate_between_per_yr'] = C['noise_pair_between'] / PY['between']
    # T1 Poisson vs scrambled in the same ETS live time
    lam = max(C['ctrl_mean_ETS'], 1e-3)
    C['T1_poisson_p'] = float(poisson.sf(C['real_ETS'] - 1, lam)) if C['real_ETS'] > 0 else 1.0
    # T2 LFE share and LFE rate vs random same-day times
    ev = D[D.tier & (D.period == 'ETS') & (D.source == 'real')].sort_values('epoch')
    ev = ev.assign(g=(ev.i.diff().fillna(1e9) > 30).cumsum()).sort_values('zmad').groupby('g').tail(1).sort_values('epoch')
    C['T2_share_with_LFEs'] = float((ev.lfe_status == 'with LFEs').mean()) if len(ev) else np.nan
    if len(ev):
        import hg_r3 as X
        LF = X.lfe_lists()
        rng = np.random.default_rng(SEED + 3)
        obs = (X.cnt(LF['svi075'], ev.epoch.values, 600) + X.cnt(LF['wa_nw'], ev.epoch.values, 600)).mean()
        rt = (np.floor(np.repeat(ev.epoch.values, 200) / 86400) * 86400 + rng.integers(0, 86400, 200 * len(ev)))
        rnd = (X.cnt(LF['svi075'], rt, 600) + X.cnt(LF['wa_nw'], rt, 600)).mean()
        C['T2_lfe_rate_ratio'] = float(obs / rnd) if rnd > 0 else np.inf
    else:
        C['T2_lfe_rate_ratio'] = np.nan
    ev.assign(time_pt=pd.to_datetime(ev.epoch, unit='s').dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles')
              .dt.strftime('%Y-%m-%d %H:%M:%S %Z'))[['time_pt', 'template', 'agree', 'zmad', 'z_best', 'ratio_far', 'n_far',
                                                     'lat_best', 'lon_best', 'nearest_family', 'vr_thrust', 'lfe_status',
                                                     'n_lfe_10m', 'i', 'epoch']].to_csv(f'{K.TB}/{NAME}{PILOT}_ets_events{SFX}.csv', index=False)
    bt = D[D.tier & (D.period == 'between') & (D.source == 'real')]
    bt.assign(time_pt=pd.to_datetime(bt.epoch, unit='s').dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles')
              .dt.strftime('%Y-%m-%d %H:%M:%S %Z'))[['time_pt', 'template', 'agree', 'zmad', 'z_best', 'ratio_far', 'vr_thrust',
                                                     'lfe_status', 'i', 'epoch']].to_csv(f'{K.TB}/{NAME}{PILOT}_between_events{SFX}.csv', index=False)
    # T3 near-station pulse on the top ETS event
    C['T3_pulse'] = None
    if len(ev):
        top = ev.sort_values('zmad').iloc[-1]
        C['T3_top_event'] = str(pd.Timestamp(top.epoch, unit='s').tz_localize('UTC').tz_convert('America/Los_Angeles'))
        C['T3_pulse'], C['T3_detail'] = t3_check(top)
    # completeness
    J = pd.read_parquet(f'{PD}/inj_list.parquet').merge(pd.read_parquet(f'{PD}/inj{SFX}.parquet'), on='inj_id', how='left')
    J['rec_tier'] = (J.mf_best_z >= LEVEL) & (J.mf_second_z >= LEVEL) & J.battery_e2.fillna(False).astype(bool)
    J['rec_mf2'] = (J.mf_best_z >= LEVEL) & (J.mf_second_z >= LEVEL)
    G = J.groupby('Mw').agg(n=('inj_id', 'count'), rec_tier=('rec_tier', 'mean'), rec_mf2=('rec_mf2', 'mean'),
                            battery_e2=('battery_e2', 'mean'), mf_best_z_median=('mf_best_z', 'median'),
                            z_best_median=('z_best', 'median')).reset_index()
    G2 = J.groupby(['period', 'Mw']).rec_tier.mean().unstack(0).reset_index()
    G.to_csv(f'{K.TB}/{NAME}{PILOT}_completeness{SFX}.csv', index=False); G2.to_csv(f'{K.TB}/{NAME}{PILOT}_completeness_by_period{SFX}.csv', index=False)
    comp = {mw: float(G[G.Mw == mw].rec_tier.iloc[0]) for mw in (3.4, 3.6, 3.8)}
    items = dict(T1=bool(C['real_ETS'] >= 5 and C['T1_poisson_p'] < 0.01),
                 T2=bool(np.isfinite(C['T2_share_with_LFEs']) and C['T2_share_with_LFEs'] >= 0.7 and C['T2_lfe_rate_ratio'] >= 3),
                 T3=bool(C['T3_pulse']), controls=bool(C['ctrl_rate_between_per_yr'] <= 1 and C['noise_pair_rate_between_per_yr'] <= 1),
                 completeness=bool(comp[3.4] >= 0.4 and comp[3.6] >= 0.7 and comp[3.8] >= 0.8))
    V = dict(pilot=PILOT, stations=STAS, statistic=stat(), counts=C, completeness=comp, items=items,
             GO=bool(all(items.values())), between_is_reported_not_gated=True)
    json.dump(V, open(f'{K.TB}/{NAME}{PILOT}_verdict{SFX}.json', 'w'), indent=1, default=str)
    pd.set_option('display.width', 250)
    print(json.dumps(V, indent=1, default=str)); print(G.round(3).to_string()); print(G2.round(3).to_string())
    figs(D, G, ev)


def t3_check(top):
    """largest near |Z| in [t0, t0 + 180] at one of the two scan stations nearest the best node, and far max < 0.5 near max."""
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import scan_refs as SR
    cc = SR.all_coords()
    SD = f'{K.S13}/data/scan' if PILOT == 2011 else f'{K.S13}/data/scan_{PILOT}'
    y0 = pd.Timestamp(f'{PILOT}-01-01').value / 1e9
    ts = int(round(top.epoch - y0)); a, n = ts - 300, 800
    near, far = [], []
    for s in ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.PFB', 'CN.YOUB', 'UW.LRIV']:
        f = f'{SD}/bp_{PILOT}g/{s}.npy'
        if os.path.isfile(f):
            x = np.asarray(np.load(f, mmap_mode='r')[2, a:a + n], float) * 1e9
            if np.isfinite(x[300:480]).any():
                near.append((s, x))
    for s in ['UW.GNW', 'IU.COR', 'US.NEW', 'US.HAWA', 'US.MSO', 'US.WVOR']:
        f = f"{os.environ.get('SCAN_REFBP', SD + '/bp2_ref')}/{s}.npy"
        if os.path.isfile(f):
            x = np.asarray(np.load(f, mmap_mode='r')[2, a:a + n], float) * 1e9
            if np.isfinite(x[300:480]).any():
                far.append((s, x))
    if not near or not np.isfinite(top.lat_best):
        return False, 'no data/location'
    d = {s: K.dkm(cc.loc[s, 'lat'], cc.loc[s, 'lon'], top.lat_best, top.lon_best) for s, _ in near}
    amp = {s: np.nanmax(np.abs(x[300:480])) for s, x in near}
    famp = {s: np.nanmax(np.abs(x[300:480])) for s, x in far if s != 'UW.GNW'}
    nearest2 = sorted(d, key=d.get)[:2]
    ok = (max(amp, key=amp.get) in nearest2) and (max(famp.values()) < 0.5 * max(amp.values()) if famp else False)
    fig, ax = plt.subplots(figsize=(10, 0.5 * (len(near) + len(far)) + 1.5)); t = np.arange(n) - 300
    for k, (s, x) in enumerate(sorted(near, key=lambda q: d[q[0]]) + far):
        pk = np.nanmax(np.abs(x)); ax.plot(t, x / (2.2 * pk) - k, color='k' if (s, x) in near else 'tab:red', lw=0.7)
        ax.text(402, -k, f'{s} peak {pk:.0f} nm' + (f' {d[s]:.0f} km' if s in d else ''), fontsize=7)
    ax.axvspan(0, 180, color='0.9', zorder=0); ax.set_yticks([]); ax.set_xlim(-300, 600)
    ax.set_title(f'pilot {PILOT} top ETS tier event {pd.Timestamp(top.epoch, unit="s").tz_localize("UTC").tz_convert("America/Los_Angeles")}'
                 f'  T3 near-pulse = {ok}', fontsize=8)
    fig.tight_layout(); fig.savefig(f'{K.FIG}/{NAME}{PILOT}_top_event{SFX}.png', dpi=120); plt.close(fig)
    return bool(ok), f'max at {max(amp, key=amp.get)} (nearest {nearest2}); far max {max(famp.values()) if famp else np.nan:.0f} vs near {max(amp.values()):.0f} nm'


def figs(D, G, ev):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8), sharey=True)
    for ax, per in zip(axs, ('ETS', 'between', 'other_episode')):
        x = D[(D.period == per) & D.template.isin(VALID)]
        for src, c, s in (('A', 'tab:purple', 4), ('real', 'tab:red', 10)):
            y = x[x.source.str.startswith(src)]
            ax.scatter(y.zmad, y.z_best.clip(0.5, 300), s=s, c=c, lw=0, alpha=0.6, label=f'{src} {len(y)}')
            yb = y[y.tier]; ax.scatter(yb.zmad, yb.z_best, s=50, facecolors='none', edgecolors=c, lw=1.2)
        ax.set_yscale('log'); ax.axhline(5, color='k', lw=0.5); ax.axvline(LEVEL, color='k', lw=0.5)
        ax.set_title(f'{PILOT} {per}: rings = tier', fontsize=9); ax.legend(fontsize=7); ax.set_xlabel('MF z (MAD)')
    axs[0].set_ylabel('scan z_best'); fig.tight_layout(); fig.savefig(f'{K.FIG}/{NAME}{PILOT}_mf_vs_scanz{SFX}.png', dpi=120)
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(G.Mw, G.rec_tier, 'o-', color='tab:red', label='tier'); ax.plot(G.Mw, G.rec_mf2, 's--', color='0.4', label='MF 2 templates')
    ax.plot(G.Mw, G.battery_e2, '^:', color='tab:blue', label='battery′')
    ax.axhline(0.8, color='k', lw=0.5); ax.set_ylim(0, 1.02); ax.set_xlabel('injected Mw'); ax.set_ylabel('fraction recovered')
    ax.set_title(f'Pilot {PILOT} cross-template injection (C0, C1; ETS + between)', fontsize=9); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f'{K.FIG}/{NAME}{PILOT}_completeness{SFX}.png', dpi=130)


if __name__ == '__main__':
    {'cont': run_cont, 'bank': run_bank, 'mf': run_mf, 'detect': run_detect, 'battery': run_battery,
     'inj_list': run_inj_list, 'inj': run_inj, 'summary': run_summary}[sys.argv[1]]()
