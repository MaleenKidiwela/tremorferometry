#!/usr/bin/env python
"""Continuous 0.02-0.05 Hz displacement per station for 1 Jun 2011 - 1 Dec 2013 (PLAN.md rows 3, 5).
Per station: day files (v3 recipe, N,E,Z m) -> (3, N) float32; per component, each finite run >= 1800 s is demeaned,
linearly detrended, cosine-tapered 120 s per end, and band-passed (Butterworth order-2 band = 4 poles, zero phase,
sosfiltfilt); the 120-s taper zones are set NaN. Outputs data/mf/cont_STA.npy (3, N) float32 with NaN = dead and
data/mf/den_STA.npy (3, M) float32 = sqrt(sum of squared demeaned samples) over each 90-s window [t, t+90), 0 if the
window touches a dead sample (then that channel's CC is 0 and it is not live).
usage: hg_cont.py STA[,STA..] [--workers 3]"""
import os, sys, argparse
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt
import hg_common as H

SOS = butter(2, H.BAND, btype='band', fs=1.0, output='sos')
TAP = 120
MINSEG = 1800


def runs(ok):
    m = np.concatenate([[False], ok, [False]]).astype(np.int8)
    d = np.diff(m)
    return np.where(d == 1)[0], np.where(d == -1)[0]


def filt(x):
    y = np.full(x.shape, np.nan, np.float32)
    a, b = runs(np.isfinite(x))
    w = np.hanning(2 * TAP)
    for i, j in zip(a, b):
        if j - i < MINSEG:
            continue
        s = x[i:j].astype(np.float64)
        t = np.arange(s.size)
        s = s - np.polyval(np.polyfit(t, s, 1), t)
        s[:TAP] *= w[:TAP]; s[-TAP:] *= w[TAP:]
        f = sosfiltfilt(SOS, s)
        f[:TAP] = np.nan; f[-TAP:] = np.nan
        y[i:j] = f
    return y


CH = 3600


def den_chunked(x):
    """sqrt of the demeaned energy of every 90-s window [t, t+90), computed with cumulative sums restarted every hour
    (2026-10-15 fix: a single global cumsum lost all precision after the 0.8-0.9 m glitches at TA.A04D/D03D on
    2012-06-07 and zeroed most later windows). 0 where the window touches a NaN."""
    out = np.zeros(H.M, np.float64)
    for a in range(0, H.M, CH):
        b = min(a + CH, H.M)
        seg = x[a:b + H.L - 1].astype(np.float64)
        bad = ~np.isfinite(seg)
        seg = np.where(bad, 0.0, seg)
        mu = seg.mean()
        seg = seg - mu                                   # local demean for precision
        cs = np.concatenate([[0.0], np.cumsum(seg)]); c2 = np.concatenate([[0.0], np.cumsum(seg * seg)])
        cb = np.concatenate([[0], np.cumsum(bad)])
        n = b - a
        s1 = cs[H.L:H.L + n] - cs[:n]; s2 = c2[H.L:H.L + n] - c2[:n]; nb = cb[H.L:H.L + n] - cb[:n]
        e = np.sqrt(np.maximum(s2 - s1 * s1 / H.L, 0.0))
        e[nb > 0] = 0.0
        out[a:b] = e
    return out


def dens(X):
    D = np.zeros((3, H.M), np.float32)
    for c in range(3):
        den = den_chunked(X[c])
        med = np.median(den[den > 0]) if np.any(den > 0) else 0
        den[den < 1e-4 * med] = 0.0                      # flat / zeroed data -> dead
        D[c] = den.astype(np.float32)
    return D


def build_den(sta):
    f = f'{H.MF}/den_{sta}.npy'
    if os.path.isfile(f):
        return sta, 'exists'
    X = np.load(f'{H.MF}/cont_{sta}.npy', mmap_mode='r')
    D = dens(X)
    np.save(f'{H.MF}/den_{sta}.part.npy', D); os.replace(f'{H.MF}/den_{sta}.part.npy', f)
    live = (D > 0).mean(axis=1)
    return sta, f'den rebuilt live N/E/Z {live[0]:.3f}/{live[1]:.3f}/{live[2]:.3f}'


def build(sta):
    out = f'{H.MF}/cont_{sta}.npy'
    if os.path.isfile(out) and os.path.isfile(f'{H.MF}/den_{sta}.npy'):
        return sta, 'exists'
    X = np.full((3, H.N), np.nan, np.float32)
    nd = 0
    for k in range(H.NDAY):
        d = H.T0 + pd.Timedelta(days=k)
        p = H.day_file(sta, d)
        if p is None:
            continue
        z = np.load(p)
        X[:, k * 86400:(k + 1) * 86400] = z['d']
        nd += 1
    for c in range(3):
        X[c] = filt(X[c])
    np.save(out + '.part.npy', X); os.replace(out + '.part.npy', out)
    D = dens(X)
    np.save(f'{H.MF}/den_{sta}.part.npy', D); os.replace(f'{H.MF}/den_{sta}.part.npy', f'{H.MF}/den_{sta}.npy')
    live = (D > 0).mean(axis=1)
    return sta, f'days {nd}/{H.NDAY} live N/E/Z {live[0]:.3f}/{live[1]:.3f}/{live[2]:.3f}'


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('stations'); ap.add_argument('--workers', type=int, default=3)
    ap.add_argument('--den-only', action='store_true')
    a = ap.parse_args()
    os.makedirs(H.MF, exist_ok=True)
    from multiprocessing import Pool
    with Pool(min(a.workers, 4)) as p:
        for s, r in p.imap_unordered(build_den if a.den_only else build, a.stations.split(',')):
            line = f'{pd.Timestamp.now(tz="America/Los_Angeles"):%Y-%m-%d %H:%M:%S} PT cont {s} {r}'
            print(line, flush=True)
            open(f'{H.LOGD}/cont.log', 'a').write(line + '\n')
