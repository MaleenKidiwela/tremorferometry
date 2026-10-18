#!/usr/bin/env python
"""VLFE scan (PLAN_SCAN.md §4): 1-s RMS envelope of the 2-8 Hz vertical (raw counts, no response), per station-day.
Used for the HF ratio column (earthquake class) and later for test 1c. Writes data/scan/hfenv/NET.STA/YYYY/DOY.npy
(86400 float32, NaN where missing). usage: scan_hfenv.py NET.STA[,..] START END [--workers 8]"""
import os, sys, argparse, warnings
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, pandas as pd
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C

OUTD = f'{C.DATA}/scan/hfenv'


def process(netsta, day):
    from obspy import read, UTCDateTime
    day = pd.Timestamp(day)
    out = f'{OUTD}/{netsta}/{day.year}/{day.dayofyear:03d}.npy'
    if os.path.isfile(out):
        return 'exists'
    fs = C.day_files(netsta, day.strftime('%Y-%m-%d'))
    if not fs:
        return 'nodata'
    t0 = UTCDateTime(str(day.date()))
    st = read(fs[0], starttime=t0 - 60, endtime=t0 + 86460)
    st = st.select(channel='BHZ') or st.select(channel='HHZ')
    if len(st) == 0:
        return 'noZ'
    for tr in st:
        sr = tr.stats.sampling_rate
        if sr > 1.5 and abs(sr - round(sr)) < 1e-3 * sr:
            tr.stats.sampling_rate = float(round(sr))
    st.merge(method=1, fill_value=None)
    st = st.split()
    E = np.full(86400, np.nan, np.float32)
    for tr in st:
        if tr.stats.npts * tr.stats.delta < 120:
            continue
        tr.data = tr.data.astype(np.float64)
        tr.detrend('demean'); tr.taper(max_percentage=0.5, max_length=5.0)
        tr.filter('bandpass', freqmin=2.0, freqmax=8.0, corners=4, zerophase=True)
        sr = int(round(tr.stats.sampling_rate))
        i0 = int(np.ceil(tr.stats.starttime - t0))          # first whole second
        k0 = int(round((t0 + i0 - tr.stats.starttime) * sr))
        x = tr.data[k0:]
        n = len(x) // sr
        if n <= 0:
            continue
        r = np.sqrt(np.mean(x[:n * sr].reshape(n, sr) ** 2, axis=1))
        a, b = max(i0, 0), min(i0 + n, 86400)
        if b > a:
            E[a:b] = r[a - i0:b - i0]
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.save(out + '.part.npy', E); os.replace(out + '.part.npy', out)
    return f'ok {np.isfinite(E).mean():.2f}'


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
    jobs = [(s, d) for d in pd.date_range(a.start, a.end, freq='D') for s in a.stations.split(',')]
    from multiprocessing import Pool
    with Pool(a.workers) as p:
        for (s, d), r in p.imap_unordered(_job, jobs):
            C.log(f'{s} {d:%Y-%m-%d} {r}', f'{C.LOGD}/scan_hfenv.log')
