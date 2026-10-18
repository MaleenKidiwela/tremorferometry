#!/usr/bin/env python
"""VLFE scan pilot 2 (PLAN_SCAN.md S7.3, S7.5, S7.7): columns for POSITIVE clusters of a scan_detect2.py run, using the
run's cleaned data (clean_iv.parquet -> NaN), plus the reference-station veto columns, z_loso, plane-wave fit,
footprint/edge flags and the pilot-2 masks.

  scan_columns2.py templates                    -> data/scan/templates_ref.npz (thrust synthetic, g = 1, all refs with
                                                   coordinates, all grid nodes)
  scan_columns2.py RUN --bp NEARBP --refbp REFBP [--minz 5] [--workers 8] [--inj inject.parquet]
                                                -> data/scan/runs/RUN/cand2.parquet
Rows: positive candidates (sign > 0) with z >= minz, not in BOTH mask v1 and mask v2; for injection runs only rows within
+-60 s of an injection. Control (a) runs: references shifted by shift*k_ctrl days (S7.3)."""
import os, sys, json, argparse
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import scan_common as S
import scan_refs as R
import scan_masks2 as M2
import s13_green as G

W = {}
TREF = f'{S.SD}/templates_ref.npz'


def _tpl_node(args):
    src, coords = args
    M = G.mt_ned(*S.MECHS['thrust'], 1.0)
    return np.stack([S.template(M, src, c) for c in coords])


def templates():
    T = R.load_table(); T = T[T.lat.notna()]
    N = pd.read_csv(f'{C.TB}/scan_grid_nodes.csv')
    coords = [(r.lat, r.lon) for r in T.itertuples()]
    from multiprocessing import Pool
    with Pool(8) as p:
        out = p.map(_tpl_node, [((r.lat, r.lon, r.dep), coords) for r in N.itertuples()])
    np.savez_compressed(TREF, stations=np.array([str(x) for x in T.netsta], dtype='U16'), thrust=np.stack(out))       # (nodes, refs, 3, 181)
    print('templates_ref', np.stack(out).shape)


def init(run, bp, refbp):
    W['run'] = run
    W['meta'] = json.load(open(f'{S.SD}/runs/{run}/meta.json'))
    W['stas'] = W['meta']['stations']; W['offs'] = np.array(W['meta']['offsets']); W['shift'] = int(W['meta']['shift'])
    W['signs'] = np.load(f'{S.SD}/runs/{run}/signs.npy')
    W['mad'] = np.load(f'{S.SD}/runs/{run}/madscale.npy')
    TP = np.load(f'{S.SD}/templates_grid.npz')
    W['T'] = {k: TP[k] * 1e21 for k in ('thrust', 'strikeslip', 'steepdip')}
    W['elem'] = TP['elem'] * 1e21
    W['nodes'] = pd.DataFrame(dict(lat=TP['lat'], lon=TP['lon'], dep=TP['dep']))
    W['BP'] = {s: np.load(f'{bp}/{s}.npy', mmap_mode='r') for s in W['stas']}
    W['DS'] = {s: dict(np.load(f'{bp}/{s}_daystats.npz')) for s in W['stas']}
    cfg = pd.read_parquet(f'{S.SD}/runs/{run}/daycfg.parquet')
    cfg['ci'] = cfg.comp.map({'N': 0, 'E': 1, 'Z': 2})
    si = {s: i for i, s in enumerate(W['stas'])}
    W['g'] = np.zeros((S.NDAY, len(W['stas']), 3)); W['won'] = np.zeros((S.NDAY, len(W['stas']), 3), bool)
    W['g'][cfg.day.values, cfg.sta.map(si).values, cfg.ci.values] = cfg.g.values
    W['won'][cfg.day.values, cfg.sta.map(si).values, cfg.ci.values] = cfg.w_on.values
    fiv = f'{S.SD}/runs/{run}/clean_iv.parquet'                     # pilot-1 runs (no cleaner): empty list (S7.12)
    iv = pd.read_parquet(fiv) if os.path.isfile(fiv) else pd.DataFrame(columns=['sta', 't0', 't1'])
    W['clean'] = {s: iv[iv.sta == s][['t0', 't1']].values for s in W['stas']}
    cc = R.all_coords()
    W['cc'] = cc
    W['sll'] = np.array([(cc.loc[s, 'lat'], cc.loc[s, 'lon']) for s in W['stas']])
    # references
    RT = R.load_table()
    from scan_refqc import usable_month
    um = usable_month()
    TR = np.load(TREF)
    tri = {s: i for i, s in enumerate(TR['stations'])}
    W['refs'] = []
    for r in RT.itertuples():
        f = f'{refbp}/{r.netsta}.npy'
        if not os.path.isfile(f) or r.netsta not in tri or r.netsta not in um:
            continue
        W['refs'].append(dict(sta=r.netsta, role=r.role, k=int(r.k_ctrl), far=str(r.set) == 'far_lhz', lat=r.lat, lon=r.lon, X=np.load(f, mmap_mode='r'),
                              sig=np.load(f'{refbp}/{r.netsta}_daystats.npz')['sig'], month_ok=um[r.netsta],
                              T=TR['thrust'][:, tri[r.netsta]] * 1e21))
    dv = json.load(open(f'{C.TB}/scan2_decay_verdict.json')) if os.path.isfile(f'{C.TB}/scan2_decay_verdict.json') else {}
    W['cfac'] = float(dv.get('veto_pred_factor_c', 1.0))
    W['refs_veto_ok'] = not dv.get('gnw_cannot_separate', False)
    from matplotlib.path import Path
    P = pd.read_csv(f'{C.TB}/scan_footprint_polygon.csv')
    W['fp'] = Path(np.c_[P.lon, P.lat]).contains_points(np.c_[W['nodes'].lon, W['nodes'].lat])
    la, lo = W['nodes'].lat.values[:-1], W['nodes'].lon.values[:-1]
    edge = np.zeros(len(W['nodes']), bool)
    # outer row/column of the grid: no grid neighbour at +-10 km in some direction among the Slab2-valid nodes
    x, y = S.kmxy(la, lo)
    for i in range(len(la)):
        nb = 0
        for dx, dy in ((10, 0), (-10, 0), (0, 10), (0, -10)):
            nb += np.any((np.abs(x - x[i] - dx) < 1) & (np.abs(y - y[i] - dy) < 1))
        edge[i] = nb < 4
    W['edge'] = edge
    W['hf'] = {}


def cleaned_window(i, s, ts_src, t_line, n0, n1):
    """band-passed station window (3, n1-n0) in nm at source time ts_src+n0.., NaN in the run's cleaned intervals
    (timeline t_line+n0..)."""
    X = W['BP'][s]
    a, b = ts_src + n0, ts_src + n1
    out = np.full((3, n1 - n0), np.nan)
    aa, bb = max(a, 0), min(b, X.shape[1])
    if bb > aa:
        out[:, aa - a:bb - a] = X[:, aa:bb] * 1e9
    iv = W['clean'][s]
    if len(iv):
        l0, l1 = t_line + n0, t_line + n1
        hit = iv[(iv[:, 1] > l0) & (iv[:, 0] < l1)]
        for x0, x1 in hit:
            out[:, max(x0 - l0, 0):min(x1 - l0, n1 - n0)] = np.nan
    return out


def hfenv(sta, day):
    k = (sta, day)
    if k not in W['hf']:
        d = S.Y0 + pd.Timedelta(days=int(day))
        f = f'{S.SD}/hfenv/{sta}/{d.year}/{d.dayofyear:03d}.npy'
        W['hf'][k] = np.load(f) if os.path.isfile(f) else None
    return W['hf'][k]


def rolling_rms_max(x, L, n=S.NT):
    """x (comps, 2L+n): max over start offsets 0..2L of the RMS over comps and n samples (comps all finite)."""
    c = np.concatenate([np.zeros((x.shape[0], 1)), np.cumsum(x * x, axis=1)], axis=1)
    s = (c[:, n:] - c[:, :-n]).sum(axis=0) / (n * x.shape[0])
    return float(np.sqrt(s.max())), int(np.argmax(s))


def cols(r):
    stas, offs = W['stas'], W['offs']
    v, tsec, node, m0 = int(r['variant']), int(r['tsec']), int(r['node']), float(r['m0'])
    d, u = divmod(tsec, 86400)
    sv = W['signs'][v]
    D = np.full((len(stas), 3, S.NT), np.nan); act = np.zeros((len(stas), 3), bool)
    wts = np.zeros((len(stas), 3)); g = W['g'][d]
    srcday = (d + offs) % S.NDAY
    for i, s in enumerate(stas):
        ts = int(srcday[i]) * 86400 + u
        x = cleaned_window(i, s, ts, tsec, S.FIT0, S.FIT1 + 1)
        sig = W['DS'][s]['sig'][srcday[i]] * 1e9
        for c in range(3):
            if W['won'][d, i, c] and np.isfinite(x[c]).all() and np.isfinite(sig[c]) and sig[c] > 0:
                act[i, c] = True; D[i, c] = x[c] * sv[i, c]; wts[i, c] = 1.0 / sig[c] ** 2
    out = dict(n_sta=int(act.any(axis=1).sum()), n_sc=int(act.sum()), fp_interior=bool(W['fp'][node]),
               grid_edge=bool(W['edge'][node]) if node < len(W['edge']) else False)
    if out['n_sc'] == 0:
        return out
    Tt = W['T']['thrust'][node] * g[:, :, None]
    w3 = wts[:, :, None] * act[:, :, None]
    Dz = np.where(act[:, :, None], D, 0.0)

    def fitm(T):
        ss = (w3 * T * T).sum(); sd = (w3 * T * Dz).sum()
        m = sd / ss if ss > 0 else np.nan
        dd = (w3 * Dz * Dz).sum()
        vr = 1 - (w3 * (Dz - m * T) ** 2).sum() / dd if dd > 0 else np.nan
        return m, vr, sd, ss, dd
    m, vr, sd, ss, dd = fitm(Tt)
    out.update(m0_cols=m * 1e12, vr_thrust=vr, cc=sd / np.sqrt(ss * dd) if ss > 0 and dd > 0 else np.nan)
    contrib = (w3 * Tt * Dz).sum(axis=2)
    shw = (w3 * Tt * Tt).sum(axis=2); shw = shw / shw.sum() if shw.sum() > 0 else shw
    if sd != 0:
        sh = contrib / sd
        k = np.unravel_index(np.argmax(np.where(act, sh, -np.inf)), sh.shape)
        ex = np.where(act, sh - shw, -np.inf)
        out.update(max_share=float(sh[k]), max_share_ch=f'{stas[k[0]]}.{"NEZ"[k[1]]}', share_excess=float(ex.max()))
    zc = (w3[:, 2] * Tt[:, 2] * Dz[:, 2]).sum(axis=1)
    out['n_Z_agree'] = int(((zc > 0) & act[:, 2]).sum())
    out['vr_strikeslip'] = fitm(W['T']['strikeslip'][node] * g[:, :, None])[1]
    out['vr_steepdip'] = fitm(W['T']['steepdip'][node] * g[:, :, None])[1]
    sw = np.sqrt(w3)
    A = np.stack([(sw * W['elem'][node, e] * g[:, :, None])[act].ravel() for e in range(5)], axis=1)
    b = (sw * Dz)[act].ravel()
    mm, *_ = np.linalg.lstsq(A, b, rcond=None)
    out['vr_free'] = 1 - ((b - A @ mm) ** 2).sum() / (b ** 2).sum() if (b ** 2).sum() > 0 else np.nan
    # z_loso (scan_loso.py rule): drop the station with the largest sum_c w T' d at the best node
    nums = (w3 * Tt * Dz).sum(axis=(1, 2)); dens = (w3 * Tt * Tt).sum(axis=(1, 2))
    mad = W['mad'][v, node, d]
    kd = int(np.argmax(nums)); keep = np.arange(len(stas)) != kd
    out['z_chk'] = nums.sum() / np.sqrt(dens.sum()) / mad if dens.sum() > 0 else np.nan
    out['z_loso'] = nums[keep].sum() / np.sqrt(dens[keep].sum()) / mad if dens[keep].sum() > 0 else np.nan
    out['loso_dropped'] = stas[kd]
    # isolation (S5 definition), on cleaned data
    isos = []
    for side in ((-900, -300), (300, 900)):
        num = 0.0; wsum = 0.0
        for i, s in enumerate(stas):
            if not act[i].any():
                continue
            ts = int(srcday[i]) * 86400 + u
            x = cleaned_window(i, s, ts, tsec, side[0], side[1])
            med = W['DS'][s]['med400'][srcday[i]] * 1e9
            for c in range(3):
                if act[i, c] and np.isfinite(med[c]) and med[c] > 0 and np.isfinite(x[c]).any():
                    num += shw[i, c] * np.sqrt(np.nanmean(x[c] ** 2)) / med[c]; wsum += shw[i, c]
        isos.append(num / wsum if wsum > 0 else np.nan)
    out['isolation'] = np.nanmax(isos) if np.isfinite(isos).any() else np.nan
    # HF ratio (class only)
    nd = S.dist_km(W['nodes'].lat[node], W['nodes'].lon[node], W['sll'][:, 0], W['sll'][:, 1])
    hfr = []
    for i in np.argsort(nd):
        s = stas[i]
        if s not in S.PILOT:
            continue
        e = hfenv(s, srcday[i])
        if e is None or not np.isfinite(e).any():
            continue
        seg = e[max(u - 60, 0):u + 61]
        if not np.isfinite(seg).any():
            continue
        hfr.append((s, np.nanmax(seg) / np.nanmedian(e)))
        if len(hfr) == 2:
            break
    out['hf1'] = hfr[0][1] if hfr else np.nan; out['hf2'] = hfr[1][1] if len(hfr) > 1 else np.nan
    out['hf_loud'] = bool(len(hfr) == 2 and min(h[1] for h in hfr) >= 10)
    # near amplitude
    out['A_near_nm'] = float(np.sqrt((shw * np.where(act, np.nanmean(np.where(act[:, :, None], D, 0.0) ** 2, axis=2), 0)).sum()))
    # reference veto columns
    nla, nlo = W['nodes'].lat[node], W['nodes'].lon[node]
    exs, rats, rr, exn, ran = [], [], {}, [], []
    for ref in W['refs']:
        dr = (d + W['shift'] * ref['k']) % S.NDAY
        if not ref['month_ok'][(S.Y0 + pd.Timedelta(days=int(dr))).month - 1]:
            continue
        dist = float(S.dist_km(nla, nlo, ref['lat'], ref['lon']))
        L = int(max(60, dist / 2.5))
        t0 = int(dr) * 86400 + u + S.FIT0 - L
        n = S.NT + 2 * L
        X = ref['X']
        if t0 < 0 or t0 + n > X.shape[1]:
            continue
        x = np.asarray(X[:, t0:t0 + n], dtype=np.float64) * 1e9
        cm = np.isfinite(x).all(axis=1)
        if not cm[2]:
            continue
        sg = ref['sig'][dr] * 1e9
        cm &= np.isfinite(sg) & (sg > 0)
        if not cm[2]:
            continue
        Ar, off = rolling_rms_max(x[cm], L)
        Pr = float(np.sqrt(np.mean((m0 / 1e12 * ref['T'][node][cm]) ** 2))) * W['cfac']
        Nr = float(np.sqrt(np.mean(sg[cm] ** 2)))
        ex = (Ar - np.sqrt(Pr ** 2 + Nr ** 2)) / Nr
        ra = Ar / out['A_near_nm'] if out['A_near_nm'] > 0 else np.nan
        tag = ref['sta'].split('.')[1]
        rr.update({f'ex_{tag}': ex, f'ratio_{tag}': ra, f'pred_{tag}': Pr, f'A_{tag}': Ar, f'lag_{tag}': off - L})
        if ref['role'] == 'veto' and (ref['sta'] != 'UW.GNW' or W['refs_veto_ok']):
            exs.append(ex); rats.append(ra)
            if not ref['far']:
                exn.append(ex); ran.append(ra)
        if ref['sta'] == 'UW.DOSE':
            out['dose_excess'] = ex; out['dose_ratio'] = ra
    out.update(rr)
    out['n_ref'] = len(exs)
    out['ref_excess'] = float(np.median(exs)) if exs else np.nan
    out['ref_ratio'] = float(np.median(rats)) if rats else np.nan
    out['ref_list'] = ','.join(k[3:] for k in rr if k.startswith('ex_'))
    out['n_ref_nf'] = len(exn)                                          # without the far set (S7.13)
    out['ref_excess_nf'] = float(np.median(exn)) if exn else np.nan
    out['ref_ratio_nf'] = float(np.median(ran)) if ran else np.nan
    # plane-wave fit (variant 0 only; column)
    if v == 0:
        out.update(plane_wave(d, u, tsec, srcday, act, D))
    return out


def plane_wave(d, u, tsec, srcday, act, D):
    stas = W['stas']
    zr = np.where(act[:, 2], np.nanmean(np.where(act[:, 2:3], D[:, 2], 0) ** 2, axis=1), -1)
    if (zr > 0).sum() < 1:
        return {}
    ip = int(np.argmax(zr))
    piv = D[ip, 2]
    tr, xy = [], []
    for i, s in enumerate(stas):
        ts = int(srcday[i]) * 86400 + u
        x = cleaned_window(i, s, ts, tsec, S.FIT0 - 150, S.FIT1 + 151)[2]
        tr.append((s, x))
    for ref in W['refs']:
        dr = (d + W['shift'] * ref['k']) % S.NDAY
        t0 = int(dr) * 86400 + u + S.FIT0 - 150
        X = ref['X']
        if t0 < 0 or t0 + S.NT + 300 > X.shape[1]:
            continue
        tr.append((ref['sta'], np.asarray(X[2, t0:t0 + S.NT + 300], dtype=np.float64) * 1e9))
    lag, cc_, pos = [], [], []
    pn = piv - piv.mean(); pnn = np.sqrt((pn ** 2).sum())
    for s, x in tr:
        if not np.isfinite(x).all() or pnn == 0:
            continue
        best = (-2, 0)
        for L in range(-150, 151):
            y = x[150 + L:150 + L + S.NT]
            y = y - y.mean(); ny = np.sqrt((y ** 2).sum())
            if ny > 0:
                c = (pn * y).sum() / (pnn * ny)
                if c > best[0]:
                    best = (c, L)
        if best[0] >= 0.6:
            la, lo = W['cc'].loc[s, 'lat'], W['cc'].loc[s, 'lon']
            pos.append(S.kmxy(la, lo)); lag.append(best[1]); cc_.append(best[0])
    o = dict(pw_n=len(lag))
    if len(lag) >= 4:
        P = np.array(pos, float); t = np.array(lag, float)
        A = np.c_[np.ones(len(t)), P[:, 0], P[:, 1]]
        sol, *_ = np.linalg.lstsq(A, t, rcond=None)
        res = t - A @ sol
        den = ((t - t.mean()) ** 2).sum()
        o.update(pw_vr=1 - (res ** 2).sum() / den if den > 0 else np.nan, pw_slow=float(np.hypot(sol[1], sol[2])),
                 pw_baz=float((np.degrees(np.arctan2(-sol[1], -sol[2])) + 360) % 360))
    return o


def run_all(rows):
    out = []
    for i, r in rows:
        try:
            out.append(dict(idx=i, **cols(r)))
        except Exception as e:  # noqa
            out.append(dict(idx=i, col_error=f'{type(e).__name__}: {str(e)[:60]}'))
    return out


def eval_rows(run, bp, refbp, Rw, workers=8):
    """Columns for an arbitrary DataFrame of rows (variant, tsec, node, m0) on run's data."""
    init(run, bp, refbp)
    from multiprocessing import Pool
    rows = list(Rw.iterrows())
    nch = max(1, min(len(rows), workers * 8))
    chunks = [rows[k::nch] for k in range(nch)]
    res = []
    with Pool(workers, initializer=init, initargs=(run, bp, refbp)) as p:
        for r in p.imap_unordered(run_all, chunks):
            res += r
    return Rw.join(pd.DataFrame(res).set_index('idx'))


def add_meta(Cd):
    t = pd.to_datetime(Cd.tsec, unit='s', origin=S.Y0)
    Cd['time_utc'] = t
    Cd['time_pt'] = t.dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.strftime('%Y-%m-%d %H:%M:%S %Z')
    N = pd.read_csv(f'{C.TB}/scan_grid_nodes.csv')
    Cd['lat'] = N.lat.values[Cd.node]; Cd['lon'] = N.lon.values[Cd.node]; Cd['dep'] = N.dep.values[Cd.node]
    Cd['mw'] = (np.log10(Cd.m0.abs()) - 9.1) / 1.5
    for tname in ('head', 'm45', 't58', 'reg', 'head2', 'bw5', 'reg3', 'v1', 'v2'):
        Cd[f'in_{tname}'] = M2.mask(tname)[Cd.tsec.values.astype(int)]
    Cd['vr_ratio'] = Cd.vr_thrust / Cd.vr_free
    return Cd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run'); ap.add_argument('--bp', default=f'{S.SD}/bp_2011g'); ap.add_argument('--refbp', default=f'{S.SD}/bp2_ref')
    ap.add_argument('--minz', type=float, default=5.0); ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--inj', default='')
    a = ap.parse_args()
    Cd = pd.read_parquet(f'{S.SD}/runs/{a.run}/cand.parquet')
    Cd = Cd[(Cd.sign > 0) & (Cd.z >= a.minz)]
    both = M2.mask('v1') & M2.mask('v2')
    Cd = Cd[~both[Cd.tsec.values]]
    if a.inj:
        J = pd.read_parquet(a.inj); tj = np.sort(J.tsec.values)
        k = np.searchsorted(tj, Cd.tsec.values)
        dmin = np.minimum(np.abs(Cd.tsec.values - tj[np.clip(k - 1, 0, len(tj) - 1)]), np.abs(Cd.tsec.values - tj[np.clip(k, 0, len(tj) - 1)]))
        Cd = Cd[dmin <= 60]
    C.log(f'scan_columns2 {a.run}: {len(Cd)} rows', f'{C.LOGD}/scan2_detect.log')
    Cd = eval_rows(a.run, a.bp, a.refbp, Cd, a.workers)
    Cd = add_meta(Cd)
    Cd.to_parquet(f'{S.SD}/runs/{a.run}/cand2.parquet')
    C.log(f'scan_columns2 {a.run}: done {len(Cd)}', f'{C.LOGD}/scan2_detect.log')


if __name__ == '__main__':
    if sys.argv[1] == 'templates':
        templates()
    else:
        main()
