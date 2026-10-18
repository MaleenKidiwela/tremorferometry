#!/usr/bin/env python
"""VLFE scan (PLAN_SCAN.md §4-5): candidate columns for one detector run (real, control a/b, injection) and the
headline flags. Writes data/scan/runs/RUN/cand_cols.parquet.  usage: scan_columns.py RUN [--bp DIR] [--workers 8]

Control (b) columns use the run's sign draw applied to the DATA per station-component (equivalent to scrambling the
template signs) with unscrambled templates for every fit. Negative candidates are evaluated with the sign-flipped
data (thrust sense), as PLAN_SCAN §5 requires for the ETS sign test."""
import os, sys, json, argparse
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import scan_common as S

THRS = [4, 4.5, 5, 5.5, 6, 7]
W = {}


def init(run, bp):
    W['run'] = run
    W['meta'] = json.load(open(f'{S.SD}/runs/{run}/meta.json'))
    W['stas'] = W['meta']['stations']; W['offs'] = np.array(W['meta']['offsets'])
    W['signs'] = np.load(f'{S.SD}/runs/{run}/signs.npy')
    W['mad'] = np.load(f'{S.SD}/runs/{run}/madscale.npy')
    TP = np.load(f'{S.SD}/templates_grid.npz')
    W['T'] = {k: TP[k] * 1e21 for k in ('thrust', 'strikeslip', 'steepdip', 'thrust20')}
    W['elem'] = TP['elem'] * 1e21                                           # (nodes, 5, sta, 3, 181)
    W['nodes'] = pd.DataFrame(dict(lat=TP['lat'], lon=TP['lon'], dep=TP['dep']))
    W['BP'] = {s: np.load(f'{bp}/{s}.npy', mmap_mode='r') for s in W['stas']}
    W['DS'] = {s: dict(np.load(f'{bp}/{s}_daystats.npz')) for s in W['stas']}
    cfg = pd.read_parquet(f'{S.SD}/runs/{run}/daycfg.parquet')
    cfg['ci'] = cfg.comp.map({'N': 0, 'E': 1, 'Z': 2})
    W['g'] = np.zeros((S.NDAY, len(W['stas']), 3)); W['won'] = np.zeros((S.NDAY, len(W['stas']), 3), bool)
    si = {s: i for i, s in enumerate(W['stas'])}
    W['g'][cfg.day.values, cfg.sta.map(si).values, cfg.ci.values] = cfg.g.values
    W['won'][cfg.day.values, cfg.sta.map(si).values, cfg.ci.values] = cfg.w_on.values
    W['mask'] = {t: S.mask_vector(t) for t in ('head', 'm45', 't58')}
    cc = C.station_coords()
    W['sll'] = np.array([(cc.loc[s, 'lat'], cc.loc[s, 'lon']) for s in W['stas']])
    ep = S.episodes_svi()
    W['major'] = [(pd.Timestamp(r.start) - pd.Timedelta(days=2), pd.Timestamp(r.end) + pd.Timedelta(days=3))
                  for r in ep[ep.major.astype(bool)].itertuples()]
    W['minor'] = [(pd.Timestamp(r.start), pd.Timestamp(r.end) + pd.Timedelta(days=1))
                  for r in ep[~ep.major.astype(bool)].itertuples()]
    W['hf'] = {}
    W['full'] = not run.startswith('ctrlA')


def hfenv(sta, day):
    k = (sta, day)
    if k not in W['hf']:
        d = S.Y0 + pd.Timedelta(days=int(day))
        f = f'{S.SD}/hfenv/{sta}/{d.year}/{d.dayofyear:03d}.npy'
        W['hf'][k] = np.load(f) if os.path.isfile(f) else None
    return W['hf'][k]


def window(arr, t0, t1):
    n = arr.shape[1]
    out = np.full((arr.shape[0], t1 - t0), np.nan)
    a, b = max(t0, 0), min(t1, n)
    if b > a:
        out[:, a - t0:b - t0] = arr[:, a:b]
    return out


def cols(r):
    stas, offs = W['stas'], W['offs']
    v, sgn, tsec, node = int(r.variant), int(r.sign), int(r.tsec), int(r.node)
    d, u = divmod(tsec, 86400)
    sv = W['signs'][v] * sgn                                       # (sta, 3): data multiplier
    D = np.full((len(stas), 3, S.NT), np.nan); act = np.zeros((len(stas), 3), bool)
    wts = np.zeros((len(stas), 3)); g = W['g'][d]; iso_num = []
    near_gap = 1e9
    srcday = (d + offs) % S.NDAY
    for i, s in enumerate(stas):
        ts = int(srcday[i]) * 86400 + u
        x = window(W['BP'][s], ts + S.FIT0, ts + S.FIT1 + 1) * 1e9
        sig = W['DS'][s]['sig'][srcday[i]] * 1e9
        for c in range(3):
            if W['won'][d, i, c] and np.isfinite(x[c]).all() and np.isfinite(sig[c]) and sig[c] > 0:
                act[i, c] = True; D[i, c] = x[c] * sv[i, c]; wts[i, c] = 1.0 / sig[c] ** 2
        if act[i].any():
            big = window(W['BP'][s], ts - 3600, ts + 3601)
            nanpos = np.where(~np.isfinite(big[act[i]]).any(axis=0))[0]
            if nanpos.size:
                near_gap = min(near_gap, np.abs(nanpos - 3600).min() + S.DIL)
    out = dict(n_sta=int(act.any(axis=1).sum()), n_sc=int(act.sum()))
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
    shw = (w3 * Tt * Tt).sum(axis=2); shw = shw / shw.sum()
    if sd != 0:
        sh = contrib / sd
        k = np.unravel_index(np.argmax(np.where(act, sh, -np.inf)), sh.shape)
        ex = np.where(act, sh - shw, -np.inf)
        out.update(max_share=float(sh[k]), max_share_ch=f'{stas[k[0]]}.{"NEZ"[k[1]]}', share_excess=float(ex.max()))
    zc = (w3[:, 2] * Tt[:, 2] * Dz[:, 2]).sum(axis=1)
    out['n_Z'] = int(act[:, 2].sum()); out['n_Z_agree'] = int(((zc > 0) & act[:, 2]).sum())
    out['vr_strikeslip'] = fitm(W['T']['strikeslip'][node] * g[:, :, None])[1]
    out['vr_steepdip'] = fitm(W['T']['steepdip'][node] * g[:, :, None])[1]
    # free deviatoric fit (same lag)
    sw = np.sqrt(w3)
    A = np.stack([(sw * W['elem'][node, e] * g[:, :, None])[act].ravel() for e in range(5)], axis=1)
    b = (sw * Dz)[act].ravel()
    mm, *_ = np.linalg.lstsq(A, b, rcond=None)
    out['vr_free'] = 1 - ((b - A @ mm) ** 2).sum() / (b ** 2).sum()
    # 20-s template, lags -10..10 (real and injection runs only)
    T20 = W['T']['thrust20'][node] * g[:, :, None]
    best = (-1, np.nan, 0)
    for L in (range(-10, 11) if W['full'] and v == 0 else []):
        Dl = np.zeros_like(Dz); okl = act.copy()
        for i, s in enumerate(stas):
            if not act[i].any():
                continue
            ts = int(srcday[i]) * 86400 + u + L
            x = window(W['BP'][s], ts + S.FIT0, ts + S.FIT1 + 1) * 1e9
            for c in range(3):
                if act[i, c]:
                    if np.isfinite(x[c]).all():
                        Dl[i, c] = x[c] * sv[i, c]
                    else:
                        okl[i, c] = False
        w3l = wts[:, :, None] * okl[:, :, None]
        s2 = (w3l * T20 * T20).sum(); sdl = (w3l * T20 * Dl).sum()
        if s2 > 0 and sdl * sdl / s2 > best[0]:
            best = (sdl * sdl / s2, sdl / s2, L)
    out.update(m0_20s=best[1] * 1e12, lag_20s=best[2])
    # isolation
    isos = []
    for side in ((-900, -300), (300, 900)):
        num = 0.0; wsum = 0.0
        for i, s in enumerate(stas):
            if not act[i].any():
                continue
            ts = int(srcday[i]) * 86400 + u
            x = window(W['BP'][s], ts + side[0], ts + side[1]) * 1e9
            med = W['DS'][s]['med400'][srcday[i]] * 1e9
            for c in range(3):
                if act[i, c] and np.isfinite(med[c]) and med[c] > 0:
                    rr = np.sqrt(np.nanmean(x[c] ** 2)) / med[c] if np.isfinite(x[c]).any() else np.nan
                    if np.isfinite(rr):
                        num += shw[i, c] * rr; wsum += shw[i, c]
        isos.append(num / wsum if wsum > 0 else np.nan)
    out['isolation'] = np.nanmax(isos) if np.isfinite(isos).any() else np.nan
    out['near_gap_s'] = near_gap
    # HF ratio at the two nearest pilot stations with envelope data
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
    out['hf_sta'] = ','.join(h[0] for h in hfr)
    out['hf1'] = hfr[0][1] if hfr else np.nan; out['hf2'] = hfr[1][1] if len(hfr) > 1 else np.nan
    out['class_eq'] = bool(len(hfr) == 2 and min(h[1] for h in hfr) >= 10)
    # z over all nodes at this time -> region of nodes within dz <= 1 of the best (real and injection runs only)
    if not (W['full'] and v == 0):
        return out
    Ta = W['T']['thrust'] * g[None, :, :, None]
    num = (Ta * (w3 * Dz)[None]).sum(axis=(1, 2, 3)); den = (Ta * Ta * w3[None]).sum(axis=(1, 2, 3))
    zn = num / np.sqrt(np.where(den > 0, den, np.nan)) / W['mad'][v, :, d]
    zn = zn[:-1]
    if np.isfinite(zn).any():
        kb = int(np.nanargmax(zn))
        reg = np.where(zn >= zn[kb] - 1)[0]
        out.update(z_cols_best=float(zn[kb]), node_cols_best=kb,
                   region_km=float(S.dist_km(W['nodes'].lat[kb], W['nodes'].lon[kb], W['nodes'].lat[reg], W['nodes'].lon[reg]).max()),
                   region_n=int(len(reg)))
    return out


def run_all(rows):
    return [dict(idx=i, **cols(r)) for i, r in rows]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run'); ap.add_argument('--bp', default=f'{S.SD}/bp_2011g'); ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--minz', type=float, default=4.0)
    a = ap.parse_args()
    Cd = pd.read_parquet(f'{S.SD}/runs/{a.run}/cand.parquet')
    Cd = Cd[Cd.z.abs() >= a.minz]
    # Amendment S5.1: joint-sign clusters (+-150 s, keep the largest |z|), per variant
    keep = []
    for v, g_ in Cd.groupby('variant'):
        g_ = g_.sort_values('tsec')
        t = g_.tsec.values; zz = g_.z.abs().values; idx = g_.index.values
        order = np.argsort(-zz); taken = np.zeros(len(g_), bool); kk = []
        for j in order:
            if taken[j]:
                continue
            kk.append(idx[j])
            lo, hi = np.searchsorted(t, t[j] - 150, 'left'), np.searchsorted(t, t[j] + 150, 'right')
            taken[lo:hi] = True
        keep += kk
    Cd = Cd.loc[sorted(keep)]
    m58 = S.mask_vector('t58')
    Cd = Cd[~m58[Cd.tsec.values]]                      # columns only where at least the narrowest tier is unmasked
    init(a.run, a.bp)
    rows = list(Cd.iterrows())
    from multiprocessing import Pool
    chunks = [rows[k::a.workers * 4] for k in range(a.workers * 4)]
    res = []
    with Pool(a.workers, initializer=init, initargs=(a.run, a.bp)) as p:
        for r in p.imap_unordered(run_all, chunks):
            res += r
    R = pd.DataFrame(res).set_index('idx')
    Cd = Cd.join(R)
    t = pd.to_datetime(Cd.tsec, unit='s', origin=S.Y0)
    Cd['time_utc'] = t
    Cd['time_pt'] = t.dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.strftime('%Y-%m-%d %H:%M:%S %Z')
    N = W['nodes']
    Cd['lat'] = N.lat.values[Cd.node]; Cd['lon'] = N.lon.values[Cd.node]; Cd['dep'] = N.dep.values[Cd.node]
    Cd['m0_abs'] = Cd.m0.abs()
    Cd['mw'] = (np.log10(Cd.m0_abs) - 9.1) / 1.5
    for tname in ('head', 'm45', 't58'):
        Cd[f'in_mask_{tname}'] = W['mask'][tname][Cd.tsec.values]
    Cd['minutes_to_midnight'] = np.minimum(Cd.tsec % 86400, 86400 - Cd.tsec % 86400) / 60
    Cd['season'] = S.season(t)
    Cd['in_ets'] = (t >= S.ETS[0]) & (t < S.ETS[1])
    Cd['in_svi_major_pm2d'] = np.any([(t >= a0) & (t < a1) for a0, a1 in W['major']], axis=0)
    Cd['in_minor_svi_episode'] = np.any([(t >= a0) & (t < a1) for a0, a1 in W['minor']], axis=0)
    Cd['vr_ratio'] = Cd.vr_thrust / Cd.vr_free
    Cd['filt_ok'] = ((Cd.n_sta >= 3) & (Cd.n_Z_agree >= 2) & (Cd.max_share <= 0.5) & (Cd.vr_ratio >= 0.8)
                     & (Cd.vr_thrust > Cd.vr_strikeslip) & (Cd.vr_thrust > Cd.vr_steepdip) & (Cd.isolation < 2)
                     & ~Cd.class_eq.fillna(False).astype(bool))
    Cd.to_parquet(f'{S.SD}/runs/{a.run}/cand_cols.parquet')
    C.log(f'scan_columns {a.run}: {len(Cd)} candidates, filt_ok {int(Cd.filt_ok.sum())}', f'{C.LOGD}/scan_detect.log')


if __name__ == '__main__':
    main()
