#!/usr/bin/env python
"""Cross-template injection machinery (PLAN_CAT Amendment P1) + early-era completeness run (P3).
  sources                 -> data/inj/sources.npz, tables/inj_sources.csv  (C0, C1 amplitude-calibrated stacks)
  check YEAR              -> local scan statistic vs stored zmax, local MF vs stored peaks (code check, SCAN_YEAR=YEAR)
  early_list              -> data/inj/early_list.parquet (injection times, 2011-12)
  early YEAR              -> data/inj/early_YEAR.parquet (SCAN_YEAR=YEAR; 8 workers)
  early_summary           -> tables/inj_early_*.csv, figures/inj_early_completeness.png
Conventions: MF data = velocity in the hg_repro frame (index i = s after 2011-06-01, cat_e3.cont); scan data = scan bp
(displacement, m) per calendar year; source windows span origin-60 .. origin+120 (181 s, the scan FIT window)."""
import os, sys, json
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

INJ = f'{K.DATA}/inj'
os.makedirs(INJ, exist_ok=True)
PIL7 = ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.PFB', 'CN.YOUB', 'UW.LRIV']
W0, W1 = -60, 121            # source window relative to origin
NW = int(os.environ.get('CAT_NW', '8'))


def epoch0_H():
    import hg_common as H
    return (H.T0 - pd.Timestamp('1970-01-01')).total_seconds()


# ------------------------------------------------------------------------------------------------ sources
def bp_year(sta, year):
    f = f'{K.S13}/data/scan/bp_2011g/{sta}.npy' if year == 2011 else f'{K.S13}/data/scan_{year}/bp_{year}g/{sta}.npy'
    return np.load(f, mmap_mode='r') if os.path.isfile(f) else None


def run_sources():
    import cat_e3 as E3
    V = pd.read_csv(f'{K.TB}/e3_events.csv')
    R3E = pd.read_csv(f'{K.HGR}/tables/r3_ets_events.csv').set_index('cand_id')
    EP0 = epoch0_H()
    out, meta = {}, []
    N = K.node_table()
    for c in (0, 1):
        g = V[V.cluster == c]
        orig = (g.i + g.shift_s).astype(int).values
        mw = R3E.loc[g.cand_id, 'Mw_tpl'].values
        # MF channels (velocity)
        mf = {}
        for ch in E3.CH:
            s, cc = ch.rsplit('.', 1); ci = 'NEZ'.index(cc)
            X, D = E3.cont(s)
            ws = []
            for o in orig:
                w = np.array(X[ci, o + W0:o + W1], np.float64)
                if np.isfinite(w).all():
                    ws.append(w - w.mean())
            if len(ws) < 2:
                continue
            U = [w / np.linalg.norm(w) for w in ws]
            st = np.mean(U, 0); st /= np.linalg.norm(st)
            a = float(np.mean([w @ st for w in ws]))
            mf[ch] = a * st
        # scan stations (displacement, scan bp of the event's calendar year)
        sc = {}
        for s in PIL7:
            ws = []
            for o in orig:
                ep = EP0 + o
                y = pd.Timestamp(ep, unit='s').year
                Xb = bp_year(s, y)
                if Xb is None:
                    continue
                a0 = int(round(ep - pd.Timestamp(f'{y}-01-01').value / 1e9))
                w = np.array(Xb[:, a0 + W0:a0 + W1], np.float64)
                if np.isfinite(w).all():
                    ws.append(w - w.mean(1, keepdims=True))
            if len(ws) < 2:
                continue
            sig = np.zeros((3, W1 - W0))
            for k in range(3):
                U = [w[k] / np.linalg.norm(w[k]) for w in ws if np.linalg.norm(w[k]) > 0]
                st = np.mean(U, 0); st /= np.linalg.norm(st)
                sig[k] = np.mean([w[k] @ st for w in ws]) * st
            sc[s] = sig
        lat, lon = g.lat_best.mean(), g.lon_best.mean()
        d = K.dkm(N.lat.values[:-1], N.lon.values[:-1], lat, lon)
        node = int(np.argmin(d))
        m0 = 10 ** (1.5 * mw[np.isfinite(mw)] + 9.1)                   # stack amplitude = mean projection ~ mean M0
        mwref = float((np.log10(m0.mean()) - 9.1) / 1.5)
        for ch, v in mf.items():
            out[f'C{c}|mf|{ch}'] = v
        for s, v in sc.items():
            out[f'C{c}|sc|{s}'] = v
        meta.append(dict(source=f'C{c}', members=','.join(map(str, g.cand_id)), n=len(g), Mw_ref=mwref,
                         Mw_members=','.join(f'{x:.2f}' for x in mw), node=node, lat=lat, lon=lon,
                         n_mf_ch=len(mf), scan_stations=','.join(sc)))
    np.savez(f'{INJ}/sources.npz', **out)
    M = pd.DataFrame(meta); M.to_csv(f'{K.TB}/inj_sources.csv', index=False)
    print(M.to_string())


def load_sources():
    z = np.load(os.environ.get('INJ_SOURCES', f'{INJ}/sources.npz'))         # Y12: alternative source file (env)
    S = {}
    for k in z.files:
        src, kind, ch = k.split('|')
        S.setdefault(src, {'mf': {}, 'sc': {}})[kind][ch] = z[k]
    M = pd.read_csv(f'{K.TB}/inj_sources.csv').set_index('source')
    return S, M


# ------------------------------------------------------------------------------------------------ local MF
def sliding_cc(x, w):
    """x (n,), w unit zero-mean (L,) -> cc for every window start (n-L+1,), 0 where the window has NaN or no energy."""
    L = len(w)
    bad = ~np.isfinite(x)
    xz = np.where(bad, 0.0, x)
    num = np.correlate(xz, w, mode='valid')
    c1 = np.concatenate([[0.0], np.cumsum(xz)]); c2 = np.concatenate([[0.0], np.cumsum(xz * xz)])
    cb = np.concatenate([[0], np.cumsum(bad)])
    s1 = c1[L:] - c1[:-L]; s2 = c2[L:] - c2[:-L]; nb = cb[L:] - cb[:-L]
    den = np.sqrt(np.maximum(s2 - s1 * s1 / L, 0.0))
    ok = (nb == 0) & (den > 0)
    return np.where(ok, num / np.where(ok, den, 1), 0.0)


def mf_local(tpls, get_data, t0, inj=None, scale=1.0, pad=300):
    """tpls: list of dict(name, W: {ch: unit (90,)}, off, med, mad). get_data(ch, a, b) -> float64 array (NaN dead).
    inj: {ch: signal (181,)} added at t0+W0. Returns {name: (zmax within +-30 s of t0, detection index)}."""
    a, b = t0 - pad, t0 + pad + 90
    chans = sorted({ch for t in tpls for ch in t['W']})
    X = {}
    for ch in chans:
        x = get_data(ch, a, b)
        if inj is not None and ch in inj:
            x = x.copy(); x[pad + W0:pad + W1] += scale * inj[ch]
        X[ch] = x
    out = {}
    for t in tpls:
        S = None
        for ch, w in t['W'].items():
            cc = sliding_cc(X[ch], w)
            S = cc if S is None else S + cc
        if S is None:
            out[t['name']] = (np.nan, -1); continue
        z = (S - t['med']) / t['mad']
        p = np.arange(a, a + len(z))
        m = np.abs(p + t['off'] - t0) <= 30
        k = int(np.argmax(np.where(m, z, -np.inf)))
        out[t['name']] = (float(z[k]), int(p[k] + t['off']))
    return out


# ------------------------------------------------------------------------------------------------ local scan battery
class Proxy:
    """array proxy adding injected samples: base (3, n); inj = list of (start, arr (3, m))."""
    def __init__(self, base, inj):
        self.base, self.inj, self.shape = base, inj, base.shape

    def __getitem__(self, key):
        c, t = key
        out = np.array(self.base[c, t], dtype=np.float64)
        t0 = t.start or 0
        t1 = t.stop if t.stop is not None else self.shape[1]
        for s, arr in self.inj:
            lo, hi = max(t0, s), min(t1, s + arr.shape[1])
            if hi <= lo:
                continue
            a = arr[c, lo - s:hi - s]
            if out.ndim == 1:
                out[lo - t0:hi - t0] += a
            else:
                out[:, lo - t0:hi - t0] += a
        return out


_B = {}


def bat_init(early=False):
    import hg_r3 as X
    X.ev_init()
    SC = X.EV['SC']; W = SC.W
    _B['X'], _B['SC'], _B['S'] = X, SC, X.EV['S']
    _B['base_BP'] = dict(W['BP'])
    _B['base_ref'] = [r['X'] for r in W['refs']]
    _B['early'] = early
    if early:
        W['won'] = W['won'].copy()
        if 'UW.LRIV' in W['stas']:
            W['won'][:, W['stas'].index('UW.LRIV'), :] = False
        W['refs'] = [r for r in W['refs'] if r['sta'] in ('IU.COR', 'US.NEW', 'US.HAWA', 'UW.GNW')]
        _B['base_ref'] = [r['X'] for r in W['refs']]
    TP = W['T']['thrust']                                   # (nodes, sta, 3, 181)
    _B['T'] = TP


def set_injection(tsec, sc_sig=None, ref_sig=None, scale=1.0):
    """install proxies: scan station s gets scale*sc_sig[s] (m) at tsec+W0 (day offsets of the real run are 0);
    references get ref_sig[sta] (m)."""
    W = _B['SC'].W
    for s in W['stas']:
        inj = [(tsec + W0, scale * sc_sig[s])] if sc_sig is not None and s in sc_sig else []
        W['BP'][s] = Proxy(_B['base_BP'][s], inj) if inj else _B['base_BP'][s]
    for r, base in zip(W['refs'], _B['base_ref']):
        inj = [(tsec + W0, ref_sig[r['sta']])] if ref_sig is not None and r['sta'] in ref_sig else []
        r['X'] = Proxy(base, inj) if inj else base


def local_z(tsec, half=30):
    """scan statistic z(t, node) for t in tsec +- half (variant 0), all nodes but the last; mirrors hg_r3.zlite."""
    SC, S = _B['SC'], _B['S']; W = SC.W
    d, u = divmod(int(tsec), 86400)
    sv = W['signs'][0]; g = W['g'][d]
    srcday = (d + W['offs']) % S.NDAY
    T = _B['T']
    nn = T.shape[0] - 1
    nt = 2 * half + 1
    num = np.zeros((nt, nn)); den = np.zeros((nt, nn))
    for i, s in enumerate(W['stas']):
        ts = int(srcday[i]) * 86400 + u
        x = SC.cleaned_window(i, s, ts, tsec, S.FIT0 - half, S.FIT1 + 1 + half)        # (3, 181 + 2 half)
        sig = W['DS'][s]['sig'][srcday[i]] * 1e9
        for c in range(3):
            if not (W['won'][d, i, c] and np.isfinite(sig[c]) and sig[c] > 0):
                continue
            xc = x[c]
            win = np.lib.stride_tricks.sliding_window_view(xc, S.NT)                 # (nt, 181)
            okt = np.isfinite(win).all(1)
            Tg = T[:nn, i, c, :] * g[i, c]                                            # (nn, 181)
            w = 1.0 / sig[c] ** 2
            num += np.where(okt[:, None], w * sv[i, c] * (np.nan_to_num(win) @ Tg.T), 0.0)
            den += np.where(okt[:, None], w * (Tg * Tg).sum(1)[None, :], 0.0)   # zlite: channel used only if finite
    mad = W['mad'][0, :nn, d]
    z = np.where(den > 0, num / np.sqrt(np.where(den > 0, den, 1)) / mad[None, :], np.nan)
    return z


def battery_local(epoch):
    """battery columns at an (injected or real) time: best point from the local scan statistic, then cols / plane wave."""
    X, SC, S = _B['X'], _B['SC'], _B['S']
    y0 = X.EV['y0']
    tsec = int(round(epoch - y0))
    z = local_z(tsec)
    if not np.isfinite(z).any():
        return dict(z_best=np.nan)
    k = np.unravel_index(np.nanargmax(z), z.shape)
    ts = tsec - 30 + int(k[0]); node = int(k[1]); zb = float(z[k])
    o = dict(z_best=zb, tsec_best=ts, node_best=node)
    c1 = SC.cols(dict(variant=0, tsec=ts, node=node, m0=1e15))
    m0 = c1.get('m0_cols', np.nan)
    c2 = SC.cols(dict(variant=0, tsec=ts, node=node, m0=m0)) if np.isfinite(m0) else c1
    for kk in ['n_sta', 'fp_interior', 'grid_edge', 'm0_cols', 'vr_thrust', 'share_excess', 'vr_free', 'vr_strikeslip',
               'vr_steepdip', 'n_Z_agree', 'z_loso', 'isolation', 'A_near_nm', 'ref_ratio',
               'loso_dropped', 'max_share', 'max_share_ch']:                 # Y11 step 1: which station carries the scan
        o[kk] = c2.get(kk, np.nan)
    rf, r6 = [], []
    for tg in K.REF_ALL:
        o[f'ratio_{tg}'] = c2.get(f'ratio_{tg}', np.nan); o[f'ex_{tg}'] = c2.get(f'ex_{tg}', np.nan)
        if np.isfinite(o[f'ratio_{tg}']):
            r6.append(o[f'ratio_{tg}'])
            if tg in K.FAR5:
                rf.append(o[f'ratio_{tg}'])
    o['n6'] = len(r6); o['ratio6'] = float(np.median(r6)) if r6 else np.nan
    o['n_far'] = len(rf); o['ratio_far'] = float(np.median(rf)) if rf else np.nan
    o.update(X.plane_wave2(ts, node))
    a0, a1 = max(tsec - 60, 0), tsec + 121
    for t in ('head2', 'bw5', 'reg', 'reg3'):
        o[f'in_{t}'] = bool(X.EV['masks'][t][a0:a1].any())
    N = K.node_table()
    o['lat_best'], o['lon_best'] = float(N.lat[node]), float(N.lon[node])
    o['dist_family_km'] = float(N.dist_family_km[node]); o['in_footprint_poly'] = bool(N.fp_interior[node])
    o['grid_edge'] = bool(N.grid_edge[node]) if not np.isfinite(o.get('grid_edge', np.nan)) else o['grid_edge']
    return o


def score(o, year):
    """battery' pass etc. for one battery_local dict."""
    D = pd.DataFrame([o])
    for c in ('n_far', 'n6', 'z_loso', 'ratio_far', 'ratio6', 'travelling', 'grid_edge', 'dist_family_km', 'n_sta',
              'n_Z_agree', 'vr_thrust', 'vr_free', 'vr_strikeslip', 'vr_steepdip', 'isolation', 'share_excess'):
        if c not in D:
            D[c] = np.nan
    D['n_far'] = D.n_far.fillna(0); D['n6'] = D.n6.fillna(0)
    D['head2_maxmag'] = K.head2_maxmag([o['epoch']], [year]) if 'epoch' in o else np.nan
    D = K.apply_rules(D)
    return D.iloc[0][['pass_a', 'pass_b', 'pass_c', 'pass_e', 'pass_e2', 'pass_f', 'battery', 'battery_e2', 'pass_S5']].to_dict()


def far_synth(src_meta, mw, scale_m0=True):
    """thrust synthetic at the far references (m) for the source's node and Mw: m0/1e12 * T (nm) -> m."""
    W = _B['SC'].W
    m0 = 10 ** (1.5 * mw + 9.1)
    out = {}
    for r in W['refs']:
        if r['far'] or r['sta'] == 'UW.GNW':
            out[r['sta']] = (m0 / 1e12) * r['T'][int(src_meta['node'])] * 1e-9          # (3, 181) m
    return out


# ------------------------------------------------------------------------------------------------ code checks
def run_check(year):
    import hg_r3 as X
    bat_init()
    EV = X.EV
    rng = np.random.default_rng(5)
    n = len(EV['zmax'])
    rows = []
    for ts in rng.integers(86400, n - 86400, 40):
        z = local_z(int(ts), half=0)
        if not np.isfinite(z).any():
            rows.append(dict(tsec=int(ts), z_local=np.nan, z_stored=float(EV['zmax'][ts]))); continue
        rows.append(dict(tsec=int(ts), z_local=float(np.nanmax(z)), z_stored=float(EV['zmax'][ts]),
                         node_local=int(np.nanargmax(z[0])), node_stored=int(EV['argmax'][ts])))
    D = pd.DataFrame(rows)
    D.to_csv(f'{K.TB}/inj_check_scan_{year}.csv', index=False)
    print(D.round(3).to_string())
    print('median |dz|', float(np.nanmedian(np.abs(D.z_local - D.z_stored))), 'node match', float((D.node_local == D.node_stored).mean()))


# ------------------------------------------------------------------------------------------------ early-era run (P3)
def early_bank(src):
    """cn4 templates of the E3 bank not containing src's members, with full-frame median/MAD and thresholds."""
    import cat_e3 as E3
    M = pd.read_csv(f'{K.TB}/e3_bank.csv')
    S, SM = load_sources()
    mem = set(str(SM.loc[src, 'members']).split(','))
    B = np.load(f'{E3.E3D}/bank.npz')
    out = []
    for r in M[M.kind.isin(['cluster', 'loo', 'hg_recut'])].itertuples():
        if mem & set(str(r.members).split(',')):
            continue
        z = np.load(f'{E3.MFD}/{r.name}__cn4.npz')
        chs = list(z['chans'])
        Wt = {ch: B[r.name][E3.CH.index(ch)] for ch in chs}
        out.append(dict(name=r.name, W=Wt, off=E3.PRE, med=float(z['real_med']), mad=float(z['real_mad'])))
    return out


def mf_getter():
    import cat_e3 as E3

    def get(ch, a, b):
        s, c = ch.rsplit('.', 1)
        X, D = E3.cont(s)
        return np.array(X['NEZ'.index(c), a:b], np.float64)
    return get


def run_early_list():
    """injection times 2011-12 (H frame): summer Jun-Sep, winter Nov-Feb; live = cn4 channels finite over the window,
    e' time masks off; >= 10 min apart and >= 10 min from real cn4 detections at the 1/day threshold."""
    import cat_e3 as E3
    import hg_common as H
    EP0 = epoch0_H()
    rng = np.random.default_rng(20261015 + 700)
    R = pd.read_parquet(f'{E3.E3D}/rows.parquet')
    det = np.sort(R[(R.set == 'cn4') & (R.source == 'real')].i.unique())
    masks = {}
    for y in (2011, 2012):
        f = f'{K.S13}/data/scan/masks2_2011.npz' if y == 2011 else f'{K.S13}/data/scan_2012/masks2_2012.npz'
        Z = np.load(f)
        masks[y] = Z['bw5'] | Z['reg'] | Z['reg3']
    A = {y: K.eq_table(y) for y in (2011, 2012)}
    rows = []
    S, SM = load_sources()
    for src in ('C0', 'C1'):
        for season, months in (('summer', (6, 7, 8, 9)), ('winter', (11, 12, 1, 2))):
            for mw in (3.4, 3.6, 3.8, 4.0, 4.2):
                got, tries = [], 0
                while len(got) < 50 and tries < 20000:
                    tries += 1
                    i = int(rng.integers(H.tidx('2011-06-02'), H.tidx('2012-12-30')))
                    t = H.T0 + pd.Timedelta(seconds=i)
                    if t.month not in months:
                        continue
                    ok = True
                    for ch in [f'{s}.{c}' for s in ['CN.PGC', 'CN.VGZ', 'CN.NLLB', 'CN.SNB'] for c in 'NEZ']:
                        s, c = ch.rsplit('.', 1)
                        X, D = E3.cont(s)
                        if not np.isfinite(X['NEZ'.index(c), i - 400:i + 400]).all():
                            ok = False; break
                    if not ok:
                        continue
                    ep = EP0 + i; y = t.year
                    ts = int(ep - pd.Timestamp(f'{y}-01-01').value / 1e9)
                    if masks[y][max(ts - 60, 0):ts + 121].any():
                        continue
                    if K.head2_maxmag([ep], [y])[0] >= 5.5:
                        continue
                    k = np.searchsorted(det, i)
                    if (k < len(det) and det[k] - i < 600) or (k > 0 and i - det[k - 1] < 600):
                        continue
                    if any(abs(i - j) < 600 for j in got) or any(abs(i - r['i']) < 600 for r in rows):
                        continue
                    got.append(i)
                for i in got:
                    rows.append(dict(source=src, season=season, Mw=mw, i=i, epoch=EP0 + i,
                                     year=(H.T0 + pd.Timedelta(seconds=i)).year))
    L = pd.DataFrame(rows)
    L['inj_id'] = np.arange(len(L))
    L.to_parquet(f'{INJ}/early_list.parquet')
    print(L.groupby(['source', 'season', 'Mw']).size())


_E = {}


def early_init():
    bat_init(early=False)
    S, SM = load_sources()
    _E['S'], _E['SM'] = S, SM
    _E['bank'] = {src: early_bank(src) for src in ('C0', 'C1')}
    _E['get'] = mf_getter()
    import cat_e3 as E3
    T = pd.read_csv(f'{K.TB}/e3_thresholds.csv')
    _E['zstar'] = {(int(r.year), float(r.rate_per_tpl_day)): r.z_star_MAD for r in T[T.set == 'cn4'].itertuples()}
    W = _B['SC'].W
    _E['full_won'] = W['won']; _E['full_refs'] = list(W['refs'])


def set_variant(early):
    W = _B['SC'].W
    if early:
        won = _E['full_won'].copy()
        if 'UW.LRIV' in W['stas']:
            won[:, W['stas'].index('UW.LRIV'), :] = False
        W['won'] = won
        W['refs'] = [r for r in _E['full_refs'] if r['sta'] in ('IU.COR', 'US.NEW', 'US.HAWA', 'UW.GNW')]
    else:
        W['won'] = _E['full_won']; W['refs'] = list(_E['full_refs'])
    _B['base_ref'] = [_B['base_ref_all'][r['sta']] for r in W['refs']]


def early_row(r):
    S, SM = _E['S'], _E['SM']
    src = r['source']; sm = SM.loc[src]
    sc = 10 ** (1.5 * (r['Mw'] - sm.Mw_ref))
    o = dict(inj_id=int(r['inj_id']))
    res = mf_local(_E['bank'][src], _E['get'], int(r['i']), inj=S[src]['mf'], scale=sc)
    zs = sorted(((v[0], k) for k, v in res.items() if np.isfinite(v[0])), reverse=True)
    o['mf_best_z'] = zs[0][0] if zs else np.nan; o['mf_best_tpl'] = zs[0][1] if zs else ''
    o['mf_second_z'] = zs[1][0] if len(zs) > 1 else np.nan
    y = int(r['year'])
    tsec = int(round(r['epoch'] - _B['X'].EV['y0']))
    for variant, early in (('early', True), ('full', False)):
        set_variant(early)
        far = far_synth(sm, r['Mw'])
        set_injection(tsec, sc_sig={s: v for s, v in S[src]['sc'].items()}, ref_sig=far, scale=sc)
        # far_synth already carries the Mw; undo the scale for refs
        for rr in _B['SC'].W['refs']:
            if isinstance(rr['X'], Proxy):
                rr['X'].inj = [(tsec + W0, far[rr['sta']])]
        b = battery_local(r['epoch']); b['epoch'] = r['epoch']
        p = score(b, y)
        for k, v in p.items():
            o[f'{variant}_{k}'] = v
        o[f'{variant}_z_best'] = b.get('z_best', np.nan); o[f'{variant}_ratio_far'] = b.get('ratio_far', np.nan)
        o[f'{variant}_n_far'] = b.get('n_far', np.nan); o[f'{variant}_z_loso'] = b.get('z_loso', np.nan)
        o[f'{variant}_dt'] = (b.get('tsec_best', np.nan) - tsec) if 'tsec_best' in b else np.nan
        o[f'{variant}_lat'] = b.get('lat_best', np.nan); o[f'{variant}_lon'] = b.get('lon_best', np.nan)
        set_injection(tsec)                                   # remove
    return o


def early_chunk(rows):
    out = []
    for r in rows:
        try:
            out.append(early_row(r))
        except Exception as e:  # noqa
            import traceback
            out.append(dict(inj_id=int(r['inj_id']), err=f'{type(e).__name__}: {str(e)[:120]}'))
    return out


def early_pool_init():
    early_init()
    W = _B['SC'].W
    _B['base_ref_all'] = {r['sta']: r['X'] for r in W['refs']}
    _B['base_ref'] = [r['X'] for r in W['refs']]


def run_early(year):
    assert os.environ.get('SCAN_YEAR') == str(year)
    L = pd.read_parquet(f'{INJ}/early_list.parquet')
    L = L[L.year == year]
    rows = L.to_dict('records')
    from multiprocessing import Pool
    nch = max(1, min(len(rows), NW * 6))
    chunks = [rows[k::nch] for k in range(nch)]
    out = []
    with Pool(NW, initializer=early_pool_init) as p:
        for k, r in enumerate(p.imap_unordered(early_chunk, chunks)):
            out += r
            if k % 12 == 0:
                K.log(f'early {year}: chunk {k + 1}/{nch}', 'inj')
    D = pd.DataFrame(out)
    D.to_parquet(f'{INJ}/early_{year}.parquet')
    K.log(f'early {year}: {len(D)} rows, errors {int(D.get("err", pd.Series(dtype=str)).notna().sum())}', 'inj')


def run_early_summary():
    L = pd.read_parquet(f'{INJ}/early_list.parquet')
    D = pd.concat([pd.read_parquet(f'{INJ}/early_{y}.parquet') for y in (2011, 2012)], ignore_index=True)
    D = L.merge(D, on='inj_id', how='left')
    T = pd.read_csv(f'{K.TB}/e3_thresholds.csv')
    zs = {(int(r.year), float(r.rate_per_tpl_day)): r.z_star_MAD for r in T[T.set == 'cn4'].itertuples()}
    for rate in (0.25, 0.5, 1.0):
        thr = np.array([zs[(int(y), rate)] for y in D.year])
        D[f'mf_{rate}'] = D.mf_best_z >= thr
        D[f'mf2_{rate}'] = D.mf_second_z >= thr
        D[f'rec_early_{rate}'] = D[f'mf_{rate}'] & D.early_battery_e2.fillna(False).astype(bool)
        D[f'rec_full_{rate}'] = D[f'mf_{rate}'] & D.full_battery_e2.fillna(False).astype(bool)
        D[f'rec_early2_{rate}'] = D[f'rec_early_{rate}'] & D[f'mf2_{rate}']
    D.to_parquet(f'{INJ}/early_scored.parquet')
    agg = {c: 'mean' for c in D.columns if c.startswith(('mf_', 'mf2_', 'rec_')) and D[c].dtype == bool}
    agg.update({'mf_best_z': 'median', 'mf_second_z': 'median', 'early_z_best': 'median'})
    agg.update({'early_pass_a': 'mean', 'early_pass_c': 'mean', 'early_pass_f': 'mean', 'early_pass_e2': 'mean',
                'early_pass_b': 'mean', 'inj_id': 'count'})
    G = D.groupby(['Mw', 'season']).agg(agg).rename(columns={'inj_id': 'n'}).reset_index()
    G.to_csv(f'{K.TB}/inj_early_completeness.csv', index=False)
    G2 = D.groupby(['Mw', 'season', 'source']).agg(agg).rename(columns={'inj_id': 'n'}).reset_index()
    G2.to_csv(f'{K.TB}/inj_early_completeness_by_source.csv', index=False)
    pd.set_option('display.width', 250)
    print(G[['Mw', 'season', 'n', 'mf_0.5', 'rec_early_0.5', 'rec_full_0.5', 'rec_early2_0.5', 'rec_early_0.25', 'rec_early_1.0',
             'early_pass_a', 'early_pass_c', 'early_pass_f', 'early_pass_e2']].round(2).to_string(index=False))
    # rule
    mwc = None
    for mw in (3.6, 3.8, 4.0, 4.2):
        g = G[G.Mw == mw]
        if len(g) == 2 and (g['rec_early_0.5'] >= 0.8).all():
            mwc = mw; break
    json.dump(dict(Mw_c_early=mwc, rule='lowest tested Mw with >= 80 % recovery (0.5/day tier, early-era battery) in both seasons'),
              open(f'{K.TB}/inj_early_verdict.json', 'w'))
    print('Mw_c(early) =', mwc)
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for season, mk in (('summer', 'o'), ('winter', 's')):
        g = G[G.season == season]
        ax.plot(g.Mw, g['mf_0.5'], ls=':', marker=mk, color='0.5', label=f'{season}: MF only')
        ax.plot(g.Mw, g['rec_full_0.5'], ls='--', marker=mk, color='tab:blue', label=f'{season}: MF + battery′ (2011-12 scan)')
        ax.plot(g.Mw, g['rec_early_0.5'], ls='-', marker=mk, color='tab:red', label=f'{season}: MF + battery′ (early-era scan)')
    ax.axhline(0.8, color='k', lw=0.5); ax.set_ylim(0, 1.02)
    ax.set_xlabel('injected Mw'); ax.set_ylabel('fraction recovered')
    ax.set_title('4-station {PGC, VGZ, NLLB, SNB} completeness, 2011-12 data, 0.5/day threshold', fontsize=9)
    ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(f'{K.FIG}/inj_early_completeness.png', dpi=130)


if __name__ == '__main__':
    m = sys.argv[1]
    if m == 'sources':
        run_sources()
    elif m == 'check':
        run_check(int(sys.argv[2]))
    elif m == 'early_list':
        run_early_list()
    elif m == 'early':
        run_early(int(sys.argv[2]))
    elif m == 'early_summary':
        run_early_summary()


# ------------------------------------------------------------------------------------------------ relay sources (Amendment U)
def run_relay_sources(year=2018, name='prod'):
    """amplitude-calibrated stacks of the accepted relay clusters of `year` -> appended to sources.npz / inj_sources.csv."""
    import shutil
    if not os.path.isfile(f'{INJ}/sources_pre_relay.npz'):
        shutil.copy(f'{INJ}/sources.npz', f'{INJ}/sources_pre_relay.npz')
        shutil.copy(f'{K.TB}/inj_sources.csv', f'{K.TB}/inj_sources_pre_relay.csv')
    if year == 2018:
        z = dict(np.load(f'{INJ}/sources_pre_relay.npz')); M = pd.read_csv(f'{K.TB}/inj_sources_pre_relay.csv')
    else:                                                     # Amendment X: append to the current set, keep a copy first
        shutil.copy(f'{INJ}/sources.npz', f'{INJ}/sources_before_{year}.npz')
        shutil.copy(f'{K.TB}/inj_sources.csv', f'{K.TB}/inj_sources_before_{year}.csv')
        z = dict(np.load(f'{INJ}/sources.npz')); M = pd.read_csv(f'{K.TB}/inj_sources.csv')
    pre, mp = ('', '') if year == 2018 else (f'R{year}', f'{year}:')
    V = pd.read_csv(f'{K.TB}/{name}{year}_relay_events.csv')
    Ev = pd.read_csv(f'{K.TB}/{name}{year}_events_v2.csv')
    mw_of = dict(zip(Ev.i, Ev.Mw_scan))
    PDY = f'{K.DATA}/{name}{year}'
    CONTY = f'{K.DATA}/prodf{year}/cont' if name.startswith('prodf') else f'{K.DATA}/prod{year}/cont'   # QC/relay runs reuse prod cont
    st_list = pd.read_csv(f'{K.TB}/prod_stations.csv').query('year == @year').sta.tolist()
    y0 = pd.Timestamp(f'{year}-01-01').value / 1e9
    N = K.node_table()
    # Amendment Y10b: W2b applies to injection sources too -- a member window in a W2b-excluded station-component-month
    # is not used for that channel (bug fix: the 2018 relay sources carried the broken UW.LRIV at 5-14 um)
    QX = pd.read_csv(f'{K.TB}/qc_exclusions.csv'); QX = QX[QX.year == year]
    excl = set(zip(QX.sta, QX.comp, QX.month))
    rows = []
    for c, g in V.groupby('cluster'):
        orig = (g.i + g.shift_s).astype(int).values
        mon = [int(pd.Timestamp(y0 + o, unit='s').month) for o in orig]
        mw = np.array([mw_of.get(int(i), np.nan) for i in g.i])
        for s in st_list:
            X = np.load(f'{CONTY}/contv_{s}.npy', mmap_mode='r')
            for ci, cc in enumerate('NEZ'):
                ws = [np.array(X[ci, o + W0:o + W1], np.float64) for o, mo in zip(orig, mon) if (s, cc, mo) not in excl]
                ws = [w - w.mean() for w in ws if np.isfinite(w).all()]
                if len(ws) < max(1, 0.5 * len(orig)):
                    continue
                U = [w / np.linalg.norm(w) for w in ws if np.linalg.norm(w) > 0]
                stk = np.mean(U, 0); stk /= np.linalg.norm(stk)
                z[f'{pre}RC{c}|mf|{s}.{cc}'] = float(np.mean([w @ stk for w in ws])) * stk
        for f in sorted(os.listdir(f'{K.S13}/data/scan_{year}/bp_{year}v2')):
            if not f.endswith('.npy'):
                continue
            s = f[:-4]
            Xb = np.load(f'{K.S13}/data/scan_{year}/bp_{year}v2/{f}', mmap_mode='r')
            ws = []
            for o, mo in zip(orig, mon):
                a0 = int(round(o))                                    # prod frame starts at 01-01 like the scan year
                w = np.array(Xb[:, a0 + W0:a0 + W1], np.float64)
                if np.isfinite(w).all():
                    ok = np.array([(s, cc, mo) not in excl for cc in 'NEZ'])     # Y10b: W2b-excluded components
                    ws.append((w - w.mean(1, keepdims=True), ok))
            if len(ws) < max(1, 0.5 * len(orig)):
                continue
            sig = np.zeros((3, W1 - W0)); nk = 0
            for k in range(3):
                wk = [w[k] for w, ok in ws if ok[k]]
                if len(wk) < max(1, 0.5 * len(orig)):
                    continue                                          # component masked: injects nothing
                U = [x / np.linalg.norm(x) for x in wk if np.linalg.norm(x) > 0]
                stk = np.mean(U, 0); stk /= np.linalg.norm(stk)
                sig[k] = np.mean([x @ stk for x in wk]) * stk; nk += 1
            if nk == 0:
                continue
            z[f'{pre}RC{c}|sc|{s}'] = sig
        m0 = 10 ** (1.5 * mw[np.isfinite(mw)] + 9.1)
        lat, lon = g.lat_best.mean(), g.lon_best.mean()
        node = int(np.argmin(K.dkm(N.lat.values[:-1], N.lon.values[:-1], lat, lon)))
        rows.append(dict(source=f'{pre}RC{c}', members=','.join(mp + x for x in g.rid), n=len(g), Mw_ref=float((np.log10(m0.mean()) - 9.1) / 1.5),
                         Mw_members=','.join(f'{x:.2f}' for x in mw), node=node, lat=lat, lon=lon, n_mf_ch=0, scan_stations=''))
    np.savez(f'{INJ}/sources.npz', **z)
    M = M[~M.source.isin([r['source'] for r in rows])]
    pd.concat([M, pd.DataFrame(rows)], ignore_index=True).to_csv(f'{K.TB}/inj_sources.csv', index=False)
    print(pd.DataFrame(rows).to_string())


if __name__ == '__main__' and sys.argv[1] == 'relay_sources':
    run_relay_sources(int(sys.argv[2]) if len(sys.argv) > 2 else 2018, sys.argv[3] if len(sys.argv) > 3 else 'prod')
