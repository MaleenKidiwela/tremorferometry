#!/usr/bin/env python
"""VLFE scan prep v3 (PLAN_SCAN.md §2; Merlin MERLIN_VLFE_SCAN.md step 2, flags F2/F3).

Same as s13_prep.py (imported read-only: read_span, join_short_gaps, pick_band, inventory, staged decimation,
sample-rate snapping, atomic write) with two changes:
  (a) F3: the 5 % taper (tr.taper(0.05) = 4,500 s of a 90,000-s padded segment) is replaced by a FIXED-length cosine
      taper of 600 s at each segment end, and remove_response is called with taper=False. The 600 s edge NaN
      (EDGE) covers the whole taper, so no taper reaches inside the retained data.
  (b) F2: per-day gap intervals are saved: 'gap_iv' = (k, 2) int array of [start, end) seconds-of-day where the RAW data
      are missing (after the <= S13_MAXGAP s linear join), 'nan_iv' = (k, 2) intervals where the stored Z is NaN
      (raw gaps + the 600 s segment-edge trim). The scanner tapers 120 s into every edge and dilates masks by 300 s.
  (c) S13_MAXGAP stays 15 s.
Writes $S13_DISP (default data/disp1sps_v3)/NET.STA/YYYY/DOY.npz  (d: (3, 86400) float32 N,E,Z [m]; chans; band;
gap_iv; nan_iv; prep='v3').  usage: scan_prep_v3.py NET.STA[,..] START END [--workers 8]
"""
import os, sys, argparse, warnings
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('S13_DISP', '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/data/disp1sps_v3')
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import s13_prep as P

TAPER_S = 600.0


def runs(mask):
    """[start, end) intervals where mask is True."""
    m = np.concatenate([[False], mask, [False]]).astype(np.int8)
    d = np.diff(m)
    return np.stack([np.where(d == 1)[0], np.where(d == -1)[0]], axis=1).astype(np.int32)


def process(netsta, day):
    from obspy import Stream
    day = pd.Timestamp(day)
    out = f'{C.DISP}/{netsta}/{day.year}/{day.dayofyear:03d}.npz'
    if os.path.isfile(out):
        return 'exists'
    st, t0 = P.read_span(netsta, day)
    if len(st) == 0:
        return 'nodata'
    st, band = P.pick_band(st)
    if st is None:
        return 'noband'
    inv = P.inventory(netsta.split('.')[0])
    if inv is None:
        return 'noresp'
    locs = pd.Series([tr.stats.location for tr in st]).value_counts()
    st = st.select(location=locs.index[0])
    for tr in st:
        sr = tr.stats.sampling_rate
        if sr > 1.5 and abs(sr - round(sr)) < 1e-3 * sr:
            tr.stats.sampling_rate = float(round(sr))
    st.merge(method=1, fill_value=None)
    st = st.split()
    if P.MAXGAP > 0:
        st = P.join_short_gaps(st)
    # raw coverage on the day grid (any channel of the chosen band, per channel then AND over channels)
    cov = {}
    for tr in st:
        g = cov.setdefault(tr.stats.channel, np.zeros(86400, bool))
        i0 = int(np.ceil(tr.stats.starttime - t0)); i1 = int(np.floor(tr.stats.endtime - t0)) + 1
        a, b = max(i0, 0), min(i1, 86400)
        if b > a:
            g[a:b] = True
    D = np.full((3, 86400), np.nan, np.float32)
    seg_out = Stream()
    for tr in st:
        if tr.stats.npts * tr.stats.delta < 3600:
            continue
        tr = tr.copy()
        tr.data = tr.data.astype(np.float64)
        tr.detrend('demean'); tr.detrend('linear')
        tr.taper(max_percentage=0.5, type='cosine', max_length=TAPER_S, side='both')     # (a) fixed 600 s
        sr = tr.stats.sampling_rate
        fac = int(round(sr))
        for f in ([8, 5] if fac == 40 else [10, 10] if fac == 100 else [5, 4] if fac == 20 else P._stages(fac, sr)):
            tr.decimate(f, no_filter=False, strict_length=False)
        try:
            tr.remove_response(inventory=inv, output='DISP', pre_filt=P.PRE_FILT, zero_mean=True, taper=False)
        except Exception:
            continue
        tr.data[:P.EDGE] = np.nan; tr.data[-P.EDGE:] = np.nan
        seg_out += tr
    if len(seg_out) == 0:
        return 'noseg'
    chans = sorted({tr.stats.channel for tr in seg_out})
    grid = {}
    for ch in chans:
        g = np.full(86400, np.nan)
        for tr in seg_out.select(channel=ch):
            i0 = int(round(tr.stats.starttime - t0))
            a, b = max(i0, 0), min(i0 + tr.stats.npts, 86400)
            if b > a:
                g[a:b] = tr.data[a - i0:b - i0]
        grid[ch] = g
    vecs = []
    for ch in chans:
        try:
            o = inv.get_orientation(f'{netsta}.{seg_out.select(channel=ch)[0].stats.location}.{ch}', t0 + 43200)
        except Exception:
            return f'noorient {ch}'
        az, dip = np.radians(o['azimuth']), np.radians(o['dip'])
        vecs.append([np.cos(dip) * np.cos(az), np.cos(dip) * np.sin(az), -np.sin(dip)])
    A = np.array(vecs)
    X = np.vstack([grid[c] for c in chans])
    ok = np.all(np.isfinite(X), axis=0)
    if len(chans) == 3:
        Y = np.full((3, 86400), np.nan)
        Y[:, ok] = np.linalg.inv(A) @ X[:, ok]
        D[:] = Y
    else:
        zi = [i for i, c in enumerate(chans) if c.endswith('Z')]
        if not zi:
            return 'noZ'
        D[2] = X[zi[0]] * np.sign(A[zi[0], 2])
    rawcov = np.ones(86400, bool)
    for ch in chans:
        rawcov &= cov.get(ch, np.zeros(86400, bool))
    gap_iv = runs(~rawcov)
    nan_iv = runs(~np.isfinite(D[2]))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    tmp = out + '.part'
    with open(tmp, 'wb') as fh:
        np.savez_compressed(fh, d=D, chans=','.join(chans), band=band, gap_iv=gap_iv, nan_iv=nan_iv, prep='v3')
    os.replace(tmp, out)
    return f'ok {band} {len(chans)}ch {np.isfinite(D[2]).mean():.2f} gaps {len(gap_iv)}'


def _job(a):
    try:
        return a, process(*a)
    except Exception as e:  # noqa
        return a, f'ERR {type(e).__name__}: {str(e)[:80]}'


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('stations'); ap.add_argument('start'); ap.add_argument('end')
    ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()
    days = pd.date_range(a.start, a.end, freq='D')
    jobs = [(s, d) for d in days for s in a.stations.split(',')]
    from multiprocessing import Pool
    lf = f'{C.LOGD}/scan_prep_v3.log'
    with Pool(a.workers) as p:
        for (s, d), r in p.imap_unordered(_job, jobs):
            C.log(f'{s} {d:%Y-%m-%d} {r}', lf)
