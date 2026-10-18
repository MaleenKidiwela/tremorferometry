#!/usr/bin/env python
"""Amendment R3 (PLAN.md): individual between-ETS Table-S1 candidates, per-event tests (a)-(g) and the selection-matched
aggregate test.
Modes
  rows               candidate set (R3.1), scrambled-control times, selection-matched noise-template detections, random
                     times for the plane-wave base rate -> tables/r3_candidates.csv, data/r3/rows.parquet
  events YEAR        scan-machinery columns (a, b, c, e, f) for every row in calendar YEAR (run with SCAN_YEAR=YEAR)
                     -> data/r3/events_YEAR.parquet
  aux                H&G-network side: (d) catalogue counts + 2-8 Hz, H&G one-station share, template-scaling Mw,
                     (g) held-out-station CC with random same-day percentiles -> data/r3/aux.parquet
  agg                R3.3 aggregate test -> tables/r3_aggregate.csv, data/r3/agg_noise.parquet
Conventions: H frame index i = seconds after 2011-06-01 (hg_common); t0 = origin-aligned index of our matched V4
detection; epoch seconds UTC are carried for the scan (calendar-year) frame."""
import os, sys, json
os.environ['OMP_NUM_THREADS'] = '1'
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/scripts')
import numpy as np
import pandas as pd

import hg_common as H

R3 = f'{H.HG}/data/r3'
os.makedirs(R3, exist_ok=True)
A_WIN, B_WIN = '2011-09-08', '2012-08-28'
# R3_TAG='ets': POST-HOC positive control (not in Amendment R3): the same pipeline on Table-S1 detections inside H&G's
# 2011 and 2012 ETS windows (template self-detections within +-300 s of a Table-1 origin removed).
TAG = os.environ.get('R3_TAG', '')
SFX = f'_{TAG}' if TAG else ''
WINS = {'': [(A_WIN, B_WIN)], 'ets': [('2011-07-23', '2011-09-06'), ('2012-08-30', '2012-10-11')]}[TAG]


def inwin(i):
    i = np.asarray(i); m = np.zeros(i.shape, bool)
    for a_, b_ in WINS:
        m |= (i >= H.tidx(a_)) & (i < H.tidx(b_))
    return m
OFF = {1: 9, 2: 12, 3: 9, 4: 9, 5: 9}
NSCR, NNOISE_CLS, NNOISE_AGG, NRND_HF, NPW_BASE = 5, 10, 100, 500, 1000
EPOCH0 = (H.T0 - pd.Timestamp('1970-01-01')).total_seconds()
HELD = ['CN.PGC', 'CN.SNB', 'CN.PFB', 'CN.YOUB']
REF6 = ['GNW', 'COR', 'NEW', 'HAWA', 'OMAK', 'I02D']
FAR5 = ['COR', 'NEW', 'HAWA', 'OMAK', 'I02D']
S13D = f'{H.S13}/data'


def log(msg):
    line = f'{pd.Timestamp.now(tz="America/Los_Angeles"):%Y-%m-%d %H:%M:%S} PT r3 {msg}'
    print(line, flush=True)
    open(f'{H.LOGD}/r3.log', 'a').write(line + '\n')


def thr_table():
    T = pd.read_csv(f'{H.TB}/r2_thresholds.csv')
    return T[(T.variant == 'V4') & (T.target_set == 'S1')].set_index('template')


# ------------------------------------------------------------------------------------------------ rows
def run_rows():
    a, b = H.tidx(WINS[0][0]), H.tidx(WINS[-1][1])
    S1 = pd.read_csv(f'{H.TB}/hg2019_tableS1.csv')
    S1 = S1[inwin(S1.i.values)].copy()
    if TAG == 'ets':
        io_all = np.array([H.tidx(t) for t in H.templates().time])
        S1 = S1[np.min(np.abs(S1.i.values[:, None] - io_all[None, :]), 1) > 300]
    ours = pd.read_csv(f'{H.TB}/r2_catalogue_chosen.csv')
    TH = thr_table()
    ours['zmad'] = (ours.S - ours.template.map(TH['median'])) / ours.template.map(TH['MAD'])
    rows = []
    for K in range(1, 6):
        o = ours[ours.template == K].sort_values('i')
        oi = o.i.values
        for r in S1[S1.template == K].itertuples():
            lo, hi = np.searchsorted(oi, r.i - 90), np.searchsorted(oi, r.i + 90)
            if hi <= lo:
                continue
            cand = o.iloc[lo:hi]
            d = cand.i.values - r.i
            k = int(np.argmin(np.abs(d - OFF[K])))
            rows.append(dict(template=K, i_S1=r.i, t_S1=r.t, i=int(cand.i.values[k]), S=float(cand.S.values[k]),
                             zmad=float(cand.zmad.values[k]), d_ours_minus_theirs=int(d[k]),
                             matched_off30=bool(abs(d[k] - OFF[K]) <= 30), matched_raw30=bool(abs(d[k]) <= 30)))
    M = pd.DataFrame(rows)
    n_S1 = len(S1)
    M = M[M.matched_off30].copy()
    n_match_rows = len(M)
    # one row per our detection (their spacing can be < 30 s): keep the closest Table-S1 row
    M['absd'] = (M.d_ours_minus_theirs - M.template.map(OFF)).abs()
    M = M.sort_values('absd').drop_duplicates(['template', 'i']).drop(columns='absd')
    # rank within our V4 between-ETS catalogue of the same template (same window)
    for K in range(1, 6):
        o = ours[(ours.template == K) & inwin(ours.i.values)].sort_values('S', ascending=False)
        rk = {int(i): r + 1 for r, i in enumerate(o.i.values)}
        m = M.template == K
        M.loc[m, 'rank_in_ours'] = M.loc[m, 'i'].map(rk)
        M.loc[m, 'n_ours'] = len(o)
    M = M.sort_values('zmad', ascending=False).reset_index(drop=True)
    M['cand_id'] = np.arange(len(M))
    M['top30'] = M.cand_id < 30
    M['epoch'] = EPOCH0 + M.i
    t = pd.to_datetime(M.epoch, unit='s')
    M['time_utc'] = t.dt.strftime('%Y-%m-%d %H:%M:%S')
    M['time_pt'] = t.dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.strftime('%Y-%m-%d %H:%M:%S %Z')
    # SVI / NW episodes +-2 d (any kind)
    E = pd.read_csv(f'{H.ROOT}/svi_science/results/b_ets_record/episodes_b.csv')
    for reg in ('SVI', 'NW'):
        e = E[E.region == reg]
        inn = np.zeros(len(M), bool)
        for x in e.itertuples():
            s0 = pd.Timestamp(x.start) - pd.Timedelta(days=2); s1 = pd.Timestamp(x.end) + pd.Timedelta(days=3)
            inn |= ((t >= s0) & (t < s1)).values
        M[f'in_{reg}_episode_pm2d'] = inn
    M.to_csv(f'{H.TB}/r3{SFX}_candidates.csv', index=False)
    log(f'rows: Table S1 in window {n_S1}; matched rows {n_match_rows}; unique candidates {len(M)} '
        f'(per template {M.template.value_counts().sort_index().to_dict()})')
    rng = np.random.default_rng(H.SEED + 300)
    R = [dict(kind='cand', cand_id=c.cand_id, template=c.template, i=c.i, draw=-1) for c in M.itertuples()]
    for c in M.itertuples():                                       # scrambled control: same UTC day, random time
        day = (c.i // 86400) * 86400
        for q in range(NSCR):
            R.append(dict(kind='scr', cand_id=c.cand_id, template=c.template, i=int(day + rng.integers(0, 86400 - 200)), draw=q))
    # selection-matched noise-template detections (rank-fraction matched)
    sel = noise_selection(M, a, b, TH)
    sel.to_parquet(f'{R3}/noise_sel{SFX}.parquet')
    for r in sel[sel.j < NNOISE_CLS].itertuples():
        R.append(dict(kind='noise', cand_id=r.cand_id, template=r.template, i=int(r.i_noise), draw=int(r.j)))
    for q, i in enumerate(rng.integers(a, b - 200, NPW_BASE)):
        R.append(dict(kind='pwbase', cand_id=-1, template=int(rng.integers(1, 6)), i=int(i), draw=q))
    R = pd.DataFrame(R)
    R['row_id'] = np.arange(len(R))
    R['epoch'] = EPOCH0 + R.i
    R['year'] = pd.to_datetime(R.epoch, unit='s').dt.year
    R.to_parquet(f'{R3}/rows{SFX}.parquet')
    log(f'rows: {R.kind.value_counts().to_dict()}; noise_sel {len(sel)} rows over {sel.j.nunique()} noise templates')


def noise_selection(M, a, b, TH):
    import hg_r2 as R2
    out = []
    for j in range(NNOISE_AGG):
        f = f'{R2.R2}/noise_{R2.tag("vel,0,sum,all,90")}_{j:03d}.npz'
        if not os.path.isfile(f):
            continue
        z = np.load(f, allow_pickle=True)
        i0 = int(z['i0'])
        for K in range(1, 6):
            pi = pv = None
            for q in range(int(z['ncfg'])):
                if K in list(z[f'c{q}_templates']):
                    pi, pv = z[f'c{q}_pi'], z[f'c{q}_pv']
            if pi is None:
                continue
            s = inwin(pi) & (np.abs(pi - i0) > H.L) & (pv >= TH.loc[K, 'thr'])
            pi_s, pv_s = pi[s], pv[s]
            o = np.argsort(pv_s)[::-1]
            pi_s, pv_s = pi_s[o], pv_s[o]
            nn = pi_s.size
            for c in M[M.template == K].itertuples():
                if nn == 0:
                    continue
                rr = int(np.clip(round(c.rank_in_ours * nn / c.n_ours), 1, nn))
                out.append(dict(j=j, i0_noise=i0, template=K, cand_id=c.cand_id, top30=c.top30, rank_noise=rr,
                                n_noise=nn, i_noise=int(pi_s[rr - 1]), S_noise=float(pv_s[rr - 1])))
    return pd.DataFrame(out)


# ------------------------------------------------------------------------------------------------ events (scan side)
def run_events(year):
    assert os.environ.get('SCAN_YEAR') == str(year), 'set SCAN_YEAR'
    R = pd.read_parquet(f'{R3}/rows{SFX}.parquet')
    R = R[R.year == year].copy()
    from multiprocessing import Pool
    import scan_common as S
    import types
    rows = [types.SimpleNamespace(**d) for d in R.to_dict('records')]
    nch = max(1, min(len(rows), 4 * 10))
    chunks = [rows[k::nch] for k in range(nch)]
    out = []
    with Pool(4, initializer=ev_init) as p:
        for k, r in enumerate(p.imap_unordered(ev_chunk, chunks)):
            out += r
            if k % 10 == 0:
                log(f'events {year}: chunk {k + 1}/{nch}')
    D = pd.DataFrame(out)
    D.to_parquet(f'{R3}/events{SFX}_{year}.parquet')
    log(f'events {year}: {len(D)} rows done')


EV = {}


def ev_init():
    import scan_common as S
    import scan_columns2 as SC
    import scan_masks2 as M2
    SC.init('real', f'{S.SD}/bp_{S.YEAR}g', f'{S.SD}/bp2_ref')
    EV['S'], EV['SC'], EV['M2'] = S, SC, M2
    EV['zmax'] = np.load(f'{S.SD}/runs/real/zmax.npy', mmap_mode='r')
    EV['argmax'] = np.load(f'{S.SD}/runs/real/argmax.npy', mmap_mode='r')
    EV['y0'] = (S.Y0 - pd.Timestamp('1970-01-01')).total_seconds()
    N = SC.W['nodes']
    TT = H.templates().set_index('template')
    EV['fixnode'] = {}
    for K in range(1, 6):
        d = S.dist_km(TT.loc[K].lat, TT.loc[K].lon, N.lat.values[:-1], N.lon.values[:-1])
        EV['fixnode'][K] = int(np.argmin(d))
    A = pd.read_csv(S.tab('scan2_eq_events.csv'))
    EV['eq'] = A
    EV['masks'] = {t: M2.mask(t) for t in ('head2', 'bw5', 'reg', 'reg3', 'm45', 'v2')}


def zlite(tsec, node):
    """z at (tsec, node), variant 0 (the z_chk part of scan_columns2.cols)."""
    S, W = EV['S'], EV['SC'].W
    d, u = divmod(int(tsec), 86400)
    stas = W['stas']; sv = W['signs'][0]; g = W['g'][d]
    srcday = (d + W['offs']) % S.NDAY
    num = den = 0.0
    for i, s in enumerate(stas):
        ts = int(srcday[i]) * 86400 + u
        x = EV['SC'].cleaned_window(i, s, ts, tsec, S.FIT0, S.FIT1 + 1)
        sig = W['DS'][s]['sig'][srcday[i]] * 1e9
        for c in range(3):
            if W['won'][d, i, c] and np.isfinite(x[c]).all() and np.isfinite(sig[c]) and sig[c] > 0:
                T = W['T']['thrust'][node][i, c] * g[i, c]
                w = 1.0 / sig[c] ** 2
                num += w * (T * x[c] * sv[i, c]).sum(); den += w * (T * T).sum()
    mad = W['mad'][0, node, d]
    return num / np.sqrt(den) / mad if den > 0 else np.nan


def plane_wave2(tsec, node):
    """S7.3 plane-wave fit with station distances (rule R3.2 b)."""
    S, SC = EV['S'], EV['SC']; W = SC.W
    d, u = divmod(int(tsec), 86400)
    srcday = (d + W['offs']) % S.NDAY
    nla, nlo = W['nodes'].lat[node], W['nodes'].lon[node]
    tr = []
    for i, s in enumerate(W['stas']):
        if not W['won'][d, i, 2]:
            continue
        ts = int(srcday[i]) * 86400 + u
        x = SC.cleaned_window(i, s, ts, tsec, S.FIT0 - 150, S.FIT1 + 151)[2]
        tr.append((s, x, True))
    for ref in W['refs']:
        dr = (d + W['shift'] * ref['k']) % S.NDAY
        t0 = int(dr) * 86400 + u + S.FIT0 - 150
        X = ref['X']
        if t0 < 0 or t0 + S.NT + 300 > X.shape[1]:
            continue
        tr.append((ref['sta'], np.asarray(X[2, t0:t0 + S.NT + 300], dtype=np.float64) * 1e9, False))
    near = [(s, x) for s, x, nr in tr if nr and np.isfinite(x[150:150 + S.NT]).all()]
    if not near:
        return dict(pw2_n=0)
    rms = [np.sqrt(np.mean(x[150:150 + S.NT] ** 2)) for _, x in near]
    piv = near[int(np.argmax(rms))][1][150:150 + S.NT]
    pn = piv - piv.mean(); pnn = np.sqrt((pn ** 2).sum())
    lag, pos, dist, names = [], [], [], []
    for s, x, _ in tr:
        if not np.isfinite(x).all() or pnn == 0:
            continue
        win = np.lib.stride_tricks.sliding_window_view(x, S.NT)          # (301, NT)
        y = win - win.mean(1, keepdims=True)
        ny = np.sqrt((y ** 2).sum(1))
        c = np.where(ny > 0, (y @ pn) / (pnn * np.where(ny > 0, ny, 1)), -2)
        k = int(np.argmax(c))
        if c[k] >= 0.6:
            la, lo = W['cc'].loc[s, 'lat'], W['cc'].loc[s, 'lon']
            pos.append(S.kmxy(la, lo)); lag.append(k - 150); names.append(s.split('.')[1])
            dist.append(float(S.dist_km(nla, nlo, la, lo)))
    o = dict(pw2_n=len(lag), pw2_nfar250=int(np.sum(np.array(dist) >= 250)) if dist else 0, pw2_list=','.join(names))
    if len(lag) >= 4:
        P = np.array(pos, float); t = np.array(lag, float)
        A = np.c_[np.ones(len(t)), P[:, 0], P[:, 1]]
        sol, *_ = np.linalg.lstsq(A, t, rcond=None)
        res = t - A @ sol
        den = ((t - t.mean()) ** 2).sum()
        sl = float(np.hypot(sol[1], sol[2]))
        o.update(pw2_vr=1 - (res ** 2).sum() / den if den > 0 else np.nan, pw2_slow=sl,
                 pw2_vapp=1 / sl if sl > 0 else np.inf,
                 pw2_baz=float((np.degrees(np.arctan2(-sol[1], -sol[2])) + 360) % 360))
    o['travelling'] = bool(o['pw2_n'] >= 8 and o.get('pw2_vr', 0) >= 0.5 and o['pw2_nfar250'] >= 2)
    return o


def eq_names(ep):
    A = EV['eq']; out = []
    for tier in ('head2', 'bw5', 'reg', 'reg3'):
        a0, a1 = A[f'{tier}_t0'].values.astype(float), A[f'{tier}_t1'].values.astype(float)
        hit = np.isfinite(a0) & (a1 > ep - 60) & (a0 < ep + 120)
        for r in A[hit].itertuples():
            out.append(f'{tier}:M{r.mag:.1f} {str(r.time)[:19]} {r.dist_km:.0f}km')
    return '; '.join(out[:4])


def ev_row(r):
    S, SC = EV['S'], EV['SC']
    ep = float(r.epoch)
    tsec = int(round(ep - EV['y0']))
    o = dict(row_id=int(r.row_id), tsec_t0=tsec)
    n = len(EV['zmax'])
    if r.kind == 'pwbase':
        k0 = EV['argmax'][tsec] if 0 <= tsec < n else -1
        node = int(k0) if k0 >= 0 else EV['fixnode'][int(r.template)]
        o.update(plane_wave2(tsec, node))
        return o
    lo, hi = max(tsec - 30, 0), min(tsec + 31, n)
    zw = np.asarray(EV['zmax'][lo:hi], np.float64)
    if np.isfinite(zw).any():
        k = int(np.nanargmax(zw)); ts = lo + k; node = int(EV['argmax'][ts]); zb = float(zw[k])
    else:
        ts, node, zb = tsec, EV['fixnode'][int(r.template)], np.nan
    o.update(tsec_best=ts, node_best=node, z_best=zb, scan_dt=ts - tsec)
    c1 = SC.cols(dict(variant=0, tsec=ts, node=node, m0=1e15))
    m0 = c1.get('m0_cols', np.nan)
    c2 = SC.cols(dict(variant=0, tsec=ts, node=node, m0=m0 if np.isfinite(m0) else 1e15)) if np.isfinite(m0) else c1
    keep = ['n_sta', 'n_sc', 'fp_interior', 'grid_edge', 'm0_cols', 'vr_thrust', 'cc', 'max_share', 'max_share_ch',
            'share_excess', 'vr_free', 'z_chk', 'z_loso', 'loso_dropped', 'isolation', 'hf1', 'hf2', 'hf_loud', 'A_near_nm',
            'n_ref', 'ref_excess', 'ref_ratio', 'n_ref_nf', 'ref_ratio_nf', 'pw_n', 'pw_vr', 'pw_slow', 'pw_baz']
    for kk in keep:
        o[kk] = c2.get(kk, np.nan)
    for tg in REF6:
        o[f'ratio_{tg}'] = c2.get(f'ratio_{tg}', np.nan); o[f'ex_{tg}'] = c2.get(f'ex_{tg}', np.nan)
    r6 = [o[f'ratio_{tg}'] for tg in REF6 if np.isfinite(o[f'ratio_{tg}'])]
    rf = [o[f'ratio_{tg}'] for tg in FAR5 if np.isfinite(o[f'ratio_{tg}'])]
    o['n6'] = len(r6); o['ratio6'] = float(np.median(r6)) if r6 else np.nan
    o['n_far'] = len(rf); o['ratio_far'] = float(np.median(rf)) if rf else np.nan
    o.update(plane_wave2(ts, node))
    fn = EV['fixnode'][int(r.template)]
    zf = [zlite(t, fn) for t in range(tsec - 10, tsec + 11)]
    o['z_fixnode'] = float(np.nanmax(zf)) if np.isfinite(zf).any() else np.nan
    o['fixnode'] = fn
    a0, a1 = max(tsec - 60, 0), min(tsec + 121, n)
    for tname, m in EV['masks'].items():
        o[f'in_{tname}'] = bool(m[a0:a1].any())
    o['eq_names'] = eq_names(ep)
    N = SC.W['nodes']
    o['lat_best'], o['lon_best'], o['dep_best'] = float(N.lat[node]), float(N.lon[node]), float(N.dep[node])
    return o


def ev_chunk(rows):
    out = []
    for r in rows:
        try:
            out.append(ev_row(r))
        except Exception as e:  # noqa
            out.append(dict(row_id=int(r.row_id), ev_error=f'{type(e).__name__}: {str(e)[:80]}'))
    return out


# ------------------------------------------------------------------------------------------------ aux (H&G side)
def lfe_lists():
    import s13_common as C
    Rw = C.svi_rows()
    Rw = Rw[(Rw.series == 'data') & Rw.counted_once]
    t0, t1 = EPOCH0 - 86400, EPOCH0 + H.tidx('2013-01-01')
    Rw = Rw[(Rw.t >= t0) & (Rw.t < t1)]
    out = {'svi075': np.sort(Rw[Rw.P >= 0.75].t.values), 'svi_c010_075': np.sort(Rw[(Rw.P >= 0.1) & (Rw.P < 0.75)].t.values)}
    c = pd.read_parquet(f'{H.ROOT}/lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet',
                        columns=['time', 'region', 'selected_v2'])
    c = c[c.selected_v2 & (c.region == 'NW')]
    tw = c.time.astype('datetime64[ns]').astype('int64').values / 1e9
    out['wa_nw'] = np.sort(tw[(tw >= t0) & (tw < t1)])
    P = pd.read_csv(f'{H.ROOT}/catalogs/pnsn_tremor_cascadia_full.csv', usecols=['time', 'lat', 'lon'])
    P['t'] = pd.to_datetime(P.time).astype('datetime64[ns]').astype('int64') / 1e9
    P = P[(P.t >= t0) & (P.t < t1)]
    TT = H.templates().set_index('template')
    for K in range(1, 6):
        d = np.hypot((P.lat - TT.loc[K].lat) * 111.19, (P.lon - TT.loc[K].lon) * 111.19 * np.cos(np.radians(TT.loc[K].lat)))
        out[f'pnsn100_T{K}'] = np.sort(P.t.values[d <= 100])
    return out


def cnt(L, t, w):
    t = np.asarray(t, np.float64)
    return np.searchsorted(L, t + w, side='right') - np.searchsorted(L, t - w, side='left')


class BP:
    """scan bp arrays (displacement 0.02-0.05 Hz, m) per calendar year, for held-out stations and hfenv."""
    def __init__(self):
        self.X = {}
        self.y0 = {2011: (pd.Timestamp('2011-01-01') - pd.Timestamp('1970-01-01')).total_seconds(),
                   2012: (pd.Timestamp('2012-01-01') - pd.Timestamp('1970-01-01')).total_seconds()}

    def arr(self, y, s):
        k = (y, s)
        if k not in self.X:
            f = f'{S13D}/scan/bp_2011g/{s}.npy' if y == 2011 else f'{S13D}/scan_2012/bp_2012g/{s}.npy'
            self.X[k] = np.load(f, mmap_mode='r') if os.path.isfile(f) else None
        return self.X[k]

    def win(self, s, ep, n=H.L):
        """(3, n) window starting at epoch ep (NaN if unavailable)."""
        y = pd.Timestamp(ep, unit='s').year
        X = self.arr(y, s)
        out = np.full((3, n), np.nan)
        if X is None:
            return out
        a = int(round(ep - self.y0[y]))
        if a < 0 or a + n > X.shape[1]:
            return out
        out[:] = X[:, a:a + n]
        return out

    def wins(self, s, eps, n=H.L):
        return np.stack([self.win(s, e, n) for e in eps])


def unit(x):
    x = x - np.nanmean(x, axis=-1, keepdims=True)
    nx = np.sqrt(np.nansum(x * x, axis=-1, keepdims=True))
    ok = np.isfinite(x).all(-1, keepdims=True) & (nx > 0)
    return np.where(ok, x / np.where(nx > 0, nx, 1), np.nan)


def heldout_ref(bp, ep_ref):
    """template (or noise) window at the held-out channels: dict ch -> unit vector or None."""
    out = {}
    for s in HELD:
        w = unit(bp.win(s, ep_ref))
        for c in range(3):
            out[(s, c)] = w[c] if np.isfinite(w[c]).all() else None
    return out


def heldout_cc(bp, ref, eps):
    """per-row mean CC over live held-out channels; also per-channel unit windows for stacking."""
    num = np.zeros(len(eps)); nl = np.zeros(len(eps))
    U = {}
    for s in HELD:
        Wn = unit(bp.wins(s, eps))                                    # (n, 3, L)
        for c in range(3):
            r = ref[(s, c)]
            if r is None:
                continue
            x = Wn[:, c]
            ok = np.isfinite(x).all(1)
            cc = np.where(ok, np.nan_to_num(x) @ r, 0.0)
            num += cc; nl += ok
            U[(s, c)] = (x, ok)
    return np.where(nl > 0, num / np.maximum(nl, 1), np.nan), nl, U


def heldout_stack(ref, U):
    vals = []
    for k, (x, ok) in U.items():
        if ok.sum() >= 2 and ref[k] is not None:
            sk = x[ok].sum(0); n = np.linalg.norm(sk)
            if n > 0:
                vals.append(float(sk @ ref[k] / n))
    return float(np.mean(vals)) if vals else np.nan


def hfenv_ratio(sta, ep, w):
    y = pd.Timestamp(ep, unit='s')
    root = f'{S13D}/scan/hfenv' if y.year == 2011 else f'{S13D}/scan_2012/hfenv'
    f = f'{root}/{sta}/{y.year}/{y.dayofyear:03d}.npy'
    if not os.path.isfile(f):
        return None
    e = np.load(f)
    med = np.nanmedian(e)
    if not np.isfinite(med) or med <= 0:
        return None
    return e, med


def run_aux():
    import hg_r2 as R2
    import hg_r2b as B
    import scan_refs as SR
    R = pd.read_parquet(f'{R3}/rows{SFX}.parquet')
    R = R[R.kind != 'pwbase'].copy()
    TT = H.templates().set_index('template')
    TH = thr_table()
    rng = np.random.default_rng(H.SEED + 310)
    LF = lfe_lists()
    cc_all = SR.all_coords()
    hf_sta = [s for s in ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.PFB', 'CN.YOUB', 'UW.LRIV']]
    bp = BP()
    out = []
    for K in range(1, 6):
        RK = R[R.template == K]
        if RK.empty:
            continue
        io = H.tidx(TT.loc[K].time); M0K = float(TT.loc[K].M0_Nm)
        # template-scaling Mw (displacement, R2c recipe) and H&G one-station share (velocity V4 template)
        Td = R2.cut(io, 'disp'); chd = R2.chans_for(K, Td, 'all')
        dd = B.Data(chd, 'disp')
        Wd = np.array([dd.X[c][io:io + H.L] for c in chd], np.float64); Wd -= Wd.mean(1, keepdims=True)
        sc = scale(dd, chd, Wd, RK.i.values)
        del dd
        iov, chv, Wv, _ = B.tmpl(K, 'vel')
        dv = B.Data(chv, 'vel')
        st = np.array([c.rsplit('.', 1)[0] for c in chv])
        ccm = np.zeros((len(RK), len(chv)))
        for k, ch in enumerate(chv):
            X = dv.cut(ch, RK.i.values, H.L).astype(np.float64)
            X = X - X.mean(1, keepdims=True); nx = np.linalg.norm(X, axis=1)
            ccm[:, k] = np.where(np.isfinite(nx) & (nx > 0), (np.nan_to_num(X) @ Wv[k]) / np.where(nx > 0, nx, 1), 0.0)
        del dv
        stas_u = sorted(set(st))
        per_sta = np.stack([ccm[:, st == s].sum(1) for s in stas_u], 1)
        Stot = per_sta.sum(1); jmax = per_sta.argmax(1)
        # held-out stations
        ref = heldout_ref(bp, EPOCH0 + io)
        g, gn, _ = heldout_cc(bp, ref, RK.epoch.values)
        # nearest two hfenv stations to the template epicentre
        dk = {s: np.hypot((cc_all.loc[s, 'lat'] - TT.loc[K].lat) * 111.19,
                          (cc_all.loc[s, 'lon'] - TT.loc[K].lon) * 111.19 * np.cos(np.radians(TT.loc[K].lat))) for s in hf_sta}
        hs = sorted(dk, key=dk.get)
        for j, r in enumerate(RK.itertuples()):
            o = dict(row_id=r.row_id, scale_tpl=sc[j], Mw_tpl=(np.log10(sc[j] * M0K) - 9.1) / 1.5 if sc[j] > 0 else np.nan,
                     hg_S_tot=Stot[j], hg_top_sta=stas_u[jmax[j]], hg_top_share=per_sta[j, jmax[j]] / Stot[j] if Stot[j] != 0 else np.nan,
                     hg_S_without_top=Stot[j] - per_sta[j, jmax[j]],
                     hg_zmad_without_top=(Stot[j] - per_sta[j, jmax[j]] - TH.loc[K, 'median']) / TH.loc[K, 'MAD'],
                     g_cc=g[j], g_nch=gn[j])
            ep = float(r.epoch)
            day0 = np.floor(ep / 86400) * 86400
            if r.kind == 'cand':
                rt = day0 + rng.integers(0, 86400 - 200, NRND_HF)
                gr, _, _ = heldout_cc(bp, ref, rt)
                o['g_pct_random'] = float(np.mean(gr[np.isfinite(gr)] < g[j])) if np.isfinite(g[j]) and np.isfinite(gr).any() else np.nan
                o['g_random_median'] = float(np.nanmedian(gr))
            else:
                rt = day0 + rng.integers(0, 86400 - 200, 50)
            # (d) catalogue counts
            for nm, L in (('svi075', LF['svi075']), ('wa_nw', LF['wa_nw']), ('c010_075', LF['svi_c010_075']),
                          ('pnsn100', LF[f'pnsn100_T{K}'])):
                for w, wl in ((60, '60s'), (600, '10m')):
                    ob = int(cnt(L, [ep], w)[0]); o[f'n_{nm}_{wl}'] = ob
                    if r.kind == 'cand':
                        rc = cnt(L, rt, w)
                        o[f'rnd_{nm}_{wl}'] = float(rc.mean()); o[f'p_{nm}_{wl}'] = float(np.mean(rc >= ob))
            # 2-8 Hz
            got = 0
            for s in hs:
                hv = hfenv_ratio(s, ep, 60)
                if hv is None:
                    continue
                e, med = hv
                u = int(ep - day0)
                got += 1
                for w, wl in ((60, '60s'), (600, '10m')):
                    seg = e[max(u - w, 0):u + w + 1]
                    v = float(np.nanmax(seg) / med) if np.isfinite(seg).any() else np.nan
                    o[f'hf{got}_{wl}'] = v
                    if r.kind == 'cand' and np.isfinite(v):
                        us = (rt - day0).astype(int)
                        rv = np.array([np.nanmax(e[max(x - w, 0):x + w + 1]) / med if np.isfinite(e[max(x - w, 0):x + w + 1]).any()
                                       else np.nan for x in us])
                        o[f'hf{got}_{wl}_pct'] = float(np.mean(rv[np.isfinite(rv)] < v)) if np.isfinite(rv).any() else np.nan
                o[f'hf{got}_sta'] = s
                if got == 2:
                    break
            out.append(o)
        log(f'aux T{K}: {len(RK)} rows')
    D = pd.DataFrame(out)
    D.to_parquet(f'{R3}/aux{SFX}.parquet')
    log(f'aux: {len(D)} rows done')


def scale(data, chs, W, starts):
    st = np.asarray(starts, np.int64)
    num = np.zeros(st.size); den = np.zeros(st.size)
    for k, ch in enumerate(chs):
        X = data.cut(ch, np.mod(st, H.N - H.L), H.L).astype(np.float64)
        ok = np.all(np.isfinite(X), 1)
        X = np.where(np.isfinite(X), X, 0.0); X = X - X.mean(1, keepdims=True)
        num += np.where(ok, X @ W[k], 0.0); den += np.where(ok, W[k] @ W[k], 0.0)
    return np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)


# ------------------------------------------------------------------------------------------------ agg (R3.3)
def run_agg():
    import hg_r2 as R2
    import hg_r2b as B
    C = pd.read_csv(f'{H.TB}/r3{SFX}_candidates.csv')
    SEL = pd.read_parquet(f'{R3}/noise_sel{SFX}.parquet')
    TT = H.templates().set_index('template')
    bp = BP()
    rng = np.random.default_rng(H.SEED + 320)
    res = []
    sets = {'TOP30': C[C.top30], 'ALL': C}
    # per template: data in RAM for H (velocity)
    for setname, CS in sets.items():
        acc = {'cand': [], 'random': {q: [] for q in range(20)}, 'noise': {}}
        for K in sorted(CS.template.unique()):
            ck = CS[CS.template == K]
            io = H.tidx(TT.loc[K].time)
            iov, chv, Wv, _ = B.tmpl(K, 'vel')
            data = B.Data(chv, 'vel')
            ref = heldout_ref(bp, EPOCH0 + io)
            g, gn, U = heldout_cc(bp, ref, EPOCH0 + ck.i.values)
            Hm = B.metrics(data, chv, Wv, ck.i.values) if len(ck) >= 3 else {'stack_cc': np.nan}
            acc['cand'].append((K, len(ck), np.nanmean(g), heldout_stack(ref, U), Hm['stack_cc'], np.nansum(g), np.sum(np.isfinite(g))))
            for q in range(20):
                tr = (ck.i.values // 86400) * 86400 + rng.integers(0, 86400 - 200, len(ck))
                gr, _, Ur = heldout_cc(bp, ref, EPOCH0 + tr)
                Hr = B.metrics(data, chv, Wv, tr) if len(ck) >= 3 else {'stack_cc': np.nan}
                acc['random'][q].append((K, len(ck), np.nanmean(gr), heldout_stack(ref, Ur), Hr['stack_cc'], np.nansum(gr), np.sum(np.isfinite(gr))))
            sk = SEL[SEL.cand_id.isin(ck.cand_id)]
            for j, sj in sk.groupby('j'):
                i0 = int(sj.i0_noise.iloc[0])
                refn = heldout_ref(bp, EPOCH0 + i0)
                gnn, _, Un = heldout_cc(bp, refn, EPOCH0 + sj.i_noise.values)
                T = R2.cut(i0, 'vel')
                Wn = np.array([T[c][0] if T[c] is not None else np.zeros(H.L) for c in chv])
                Hn = B.metrics(data, chv, Wn, sj.i_noise.values) if len(sj) >= 3 else {'stack_cc': np.nan}
                acc['noise'].setdefault(j, []).append((K, len(sj), np.nanmean(gnn), heldout_stack(refn, Un), Hn['stack_cc'],
                                                       np.nansum(gnn), np.sum(np.isfinite(gnn))))
            del data
            log(f'agg {setname} T{K} done ({len(ck)} candidates)')

        def pool(lst):
            n = np.array([x[1] for x in lst], float)
            gsum = np.nansum([x[5] for x in lst]); gcnt = np.nansum([x[6] for x in lst])
            gst = np.array([x[3] for x in lst]); hh = np.array([x[4] for x in lst])
            okg = np.isfinite(gst); okh = np.isfinite(hh)
            return dict(G=gsum / gcnt if gcnt > 0 else np.nan,
                        G_stack=float(np.sum(gst[okg] * n[okg]) / n[okg].sum()) if okg.any() else np.nan,
                        H=float(np.sum(hh[okh] * n[okh]) / n[okh].sum()) if okh.any() else np.nan, n=int(n.sum()))
        res.append(dict(set=setname, source='candidates', draw=-1, **pool(acc['cand'])))
        for q, lst in acc['random'].items():
            res.append(dict(set=setname, source='random same-day', draw=q, **pool(lst)))
        for j, lst in acc['noise'].items():
            res.append(dict(set=setname, source='selection-matched noise template', draw=j, **pool(lst)))
    D = pd.DataFrame(res)
    D.to_parquet(f'{R3}/agg_all{SFX}.parquet')
    rows = []
    for setname in sets:
        d = D[D.set == setname]
        c = d[d.source == 'candidates'].iloc[0]
        for src in ['random same-day', 'selection-matched noise template']:
            x = d[d.source == src]
            r = dict(set=setname, control=src, n_candidates=int(c.n), n_draws=len(x))
            for st in ['G', 'G_stack', 'H']:
                v = x[st].dropna().values
                r[f'{st}_cand'] = c[st]
                r[f'{st}_ctrl_median'] = float(np.median(v)) if v.size else np.nan
                r[f'{st}_ctrl_p16'], r[f'{st}_ctrl_p84'] = (np.percentile(v, [16, 84]) if v.size else (np.nan, np.nan))
                r[f'{st}_ctrl_p975'] = float(np.percentile(v, 97.5)) if v.size else np.nan
                r[f'{st}_cand_pct'] = float(np.mean(v < c[st])) if v.size else np.nan
            rows.append(r)
    A = pd.DataFrame(rows)
    A.to_csv(f'{H.TB}/r3{SFX}_aggregate.csv', index=False)
    pd.set_option('display.width', 250)
    print(A.round(3).to_string(index=False))
    log('agg done')


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'rows':
        run_rows()
    elif mode == 'events':
        run_events(int(sys.argv[2]))
    elif mode == 'aux':
        run_aux()
    elif mode == 'agg':
        run_agg()
