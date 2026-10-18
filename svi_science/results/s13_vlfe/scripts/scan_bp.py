#!/usr/bin/env python
"""VLFE scan (PLAN_SCAN.md §4): band-passed 2011 year per station (v3, per-run detrend + 120-s taper + 0.02-0.05 Hz
filtfilt, 300-s mask dilation, vetoes), plus per-day per-component sigma (1.4826 MAD, finite, outside the headline mask)
and daily median 400-s RMS. Optional injection list (scan_inject.py) added to the raw displacement BEFORE band-passing.
Writes <out>/NET.STA.npy (3, 365*86400 float32), <out>/NET.STA_daystats.npz (sig (365,3), med400 (365,3)).
usage: scan_bp.py OUTDIR [--inject inj.parquet] [--stations ...] [--workers 8]"""
import os, sys, argparse
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import s13_green as G
import scan_common as S

STAS = S.PILOT + S.NOISE_ONLY
GLITCH_K = 200.0
GLITCH_DIL = 600


def make_injector(inj_file):
    if not inj_file:
        return None
    J = pd.read_parquet(inj_file)
    cc = C.station_coords()
    y0 = S.Y0.value / 1e9

    def inject(netsta, X):
        delta = np.zeros_like(X)
        sll = (cc.loc[netsta, 'lat'], cc.loc[netsta, 'lon'])
        gains = pd.read_csv(f'{C.TB}/sc_gain.csv')
        gains = gains[gains.sta == netsta]
        for r in J.itertuples():
            t0 = r.t - y0
            i0 = int(np.floor(t0)) - 400
            tt = np.arange(i0, i0 + 1600) - t0                     # seconds after origin
            M = G.mt_ned(*S.MECHS['thrust'], r.m0)
            kw = dict(free_surface=r.fs, vs=r.vs * 1e3, vp=r.vp * 1e3)
            u = S.raw_disp(M, (r.lat, r.lon, r.dep), sll, tt, tau=r.tau, **kw)
            u = u - u[:, :1]                                         # zero before the event (step after)
            d = pd.Timestamp(r.t, unit='s')
            g = np.ones(3)
            for c, cn in enumerate('NEZ'):
                gg = gains[(gains.comp == cn) & (gains.epoch_start <= str(d)) & (gains.epoch_end >= str(d))]
                g[c] = gg.g.iloc[0] if len(gg) and gg.keep.iloc[0] else 1.0   # true Earth response = g x synthetic
            a, b = max(i0, 0), min(i0 + 1600, X.shape[1])
            if b <= a:
                continue
            seg = u[:, a - i0:b - i0] * g[:, None]
            # the static step persists after the window (to the end of the year, via a cumulative offset)
            X[:, a:b] += seg
            if b < X.shape[1]:
                delta[:, b] += u[:, -1] * g
        X += np.cumsum(delta, axis=1)
    return inject


def one(args):
    netsta, out, inj = args
    f = f'{out}/{netsta}.npy'
    if os.path.isfile(f):
        return netsta, 'exists'
    X = S.load_year(netsta, inject=make_injector(inj))
    B = S.bandpass_year(X)
    del X
    mask = S.mask_vector('head')
    sig = np.full((S.NDAY, 3), np.nan); med = np.full((S.NDAY, 3), np.nan)
    for d in range(S.NDAY):
        sl = slice(d * 86400, (d + 1) * 86400)
        for c in range(3):
            x = B[c, sl]
            ok = np.isfinite(x) & ~mask[sl]
            if ok.mean() < 0.1:
                ok = np.isfinite(x)
            if ok.sum() > 3600:
                v = x[ok]
                sig[d, c] = 1.4826 * np.median(np.abs(v - np.median(v)))
            w = x.reshape(-1, 400) if x.size % 400 == 0 else x[:86400 // 400 * 400].reshape(-1, 400)
            okw = np.isfinite(w).mean(axis=1) >= 0.95
            if okw.any():
                med[d, c] = np.median(np.sqrt(np.nanmean(w[okw] ** 2, axis=1)))
    # Amendment S4: glitch mask, |x| > 200 sigma_day -> NaN +-600 s (per component)
    from scipy.ndimage import binary_dilation
    gl = np.zeros((S.NDAY, 3))
    for c in range(3):
        sg = np.repeat(sig[:, c], 86400)
        with np.errstate(invalid='ignore'):
            hit = np.abs(B[c]) > GLITCH_K * sg
        if hit.any():
            m = binary_dilation(hit, structure=np.ones(2 * GLITCH_DIL + 1, bool)) & np.isfinite(B[c])
            gl[:, c] = m.reshape(S.NDAY, 86400).sum(axis=1)
            B[c, m] = np.nan
    np.save(f + '.part.npy', B.astype(np.float32)); os.replace(f + '.part.npy', f)
    np.savez(f'{out}/{netsta}_daystats.npz', sig=sig, med400=med, glitch_s=gl)
    return netsta, f'ok finite {np.isfinite(B[2]).mean():.3f}'


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('--inject', default=''); ap.add_argument('--stations', default=','.join(STAS))
    ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    from multiprocessing import Pool
    with Pool(a.workers) as p:
        for s, r in p.imap_unordered(one, [(s, a.out, a.inject) for s in a.stations.split(',')]):
            C.log(f'scan_bp {a.out}: {s} {r}', f'{C.LOGD}/scan_bp.log')
