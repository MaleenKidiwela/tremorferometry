#!/usr/bin/env python
"""S13 prep: raw day file -> 1-sps displacement N/E/Z (float32, NaN where missing), per station-day.

usage: s13_prep.py NET.STA[,NET.STA...] START END [--workers 8] [--days-file F]
Reads day d plus 1800 s padding from d-1 / d+1 (same roots), splits on gaps, per contiguous segment >= 1 h:
demean, detrend, 5 % taper, decimate to 1 sps (FIR/cheby via obspy), remove response to DISP with pre-filt
(0.004, 0.008, 0.25, 0.45) Hz, rotate ->ZNE with the StationXML; 600 s at every segment edge -> NaN; trim to the day.
Writes data/disp1sps/NET.STA/YYYY/DOY.npz (d: (3, 86400) float32 N,E,Z; meta: channels, sr). Skips existing outputs.
Channel preference: BH? then HH? (EH*/HN*/borehole geophones never used)."""
import os, sys, argparse, warnings
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C

PRE_FILT = (0.004, 0.008, 0.25, 0.45)
PAD = 1800.0
EDGE = 600
MAXGAP = float(os.environ.get('S13_MAXGAP', '15'))
_INV = {}



def _stages(fac, sr):
    """Split a decimation factor into stages <= 10 (obspy refuses > 16 in one step). 2026-10-14: GNW 2007-11 failed."""
    if abs(sr - fac) > 1e-6:
        raise ValueError(f"non-integer sampling rate {sr}")
    out, n = [], fac
    for p in (10, 8, 5, 4, 2, 3, 7):
        while n % p == 0 and n > 1:
            out.append(p); n //= p
    if n != 1:
        raise ValueError(f"cannot stage decimation factor {fac}")
    return out

def inventory(net):
    from obspy import read_inventory
    if net not in _INV:
        f = f'{C.RESP}/{net}_bb.xml'
        _INV[net] = read_inventory(f) if os.path.isfile(f) else None
    return _INV[net]


def read_span(netsta, day):
    from obspy import read, Stream, UTCDateTime
    t0 = UTCDateTime(str(day)); t1 = t0 + 86400
    st = Stream()
    for dd, a, b in ((day - pd.Timedelta(days=1), t0 - PAD, t0), (day, t0, t1), (day + pd.Timedelta(days=1), t1, t1 + PAD)):
        fs = C.day_files(netsta, dd.strftime('%Y-%m-%d'))
        if not fs:
            continue
        try:
            s = read(fs[0], starttime=a, endtime=b)
        except Exception:
            continue
        st += s
    return st, t0


def join_short_gaps(st, maxgap=MAXGAP):
    """Per channel: join consecutive traces whose gap is <= maxgap s by linear interpolation (Amendment 2)."""
    from obspy import Stream, Trace
    out = Stream()
    for ch in sorted({tr.id for tr in st}):
        trs = sorted(st.select(id=ch), key=lambda t: t.stats.starttime)
        cur = trs[0].copy(); cur.data = cur.data.astype(np.float64)
        for tr in trs[1:]:
            dt = cur.stats.delta
            gap = tr.stats.starttime - (cur.stats.endtime + dt)
            if -dt < gap <= maxgap:
                n = int(round(gap / dt))
                fill = np.linspace(cur.data[-1], tr.data[0], n + 2)[1:-1] if n > 0 else np.array([])
                cur.data = np.concatenate([cur.data, fill, tr.data.astype(np.float64)])
            elif gap <= -dt:                     # overlap: drop overlapping samples of the later trace
                k = int(round(-gap / dt))
                if k < tr.stats.npts:
                    cur.data = np.concatenate([cur.data, tr.data[k:].astype(np.float64)])
            else:
                out += cur
                cur = tr.copy(); cur.data = cur.data.astype(np.float64)
        out += cur
    return out


def pick_band(st):
    for band in ('BH', 'HH'):
        s = st.select(channel=f'{band}?')
        comps = {tr.stats.channel[-1] for tr in s}
        if len(comps) >= 3 or ('Z' in comps and len(comps) >= 1):
            return s, band
    return None, None


def process(netsta, day):
    from obspy import Stream
    day = pd.Timestamp(day)
    out = f'{C.DISP}/{netsta}/{day.year}/{day.dayofyear:03d}.npz'
    if os.path.isfile(out):
        return 'exists'
    st, t0 = read_span(netsta, day)
    if len(st) == 0:
        return 'nodata'
    st, band = pick_band(st)
    if st is None:
        return 'noband'
    inv = inventory(netsta.split('.')[0])
    if inv is None:
        return 'noresp'
    # keep one location code (the most common)
    locs = pd.Series([tr.stats.location for tr in st]).value_counts()
    st = st.select(location=locs.index[0])
    for tr in st:                                # 2026-10-14: CN.PGC files carry 99.99999/100.00001 sps -> merge refused
        sr = tr.stats.sampling_rate
        if sr > 1.5 and abs(sr - round(sr)) < 1e-3 * sr:
            tr.stats.sampling_rate = float(round(sr))
    st.merge(method=1, fill_value=None)
    st = st.split()
    if MAXGAP > 0:
        st = join_short_gaps(st)
    D = np.full((3, 86400), np.nan, np.float32)
    seg_out = Stream()
    for tr in st:
        if tr.stats.npts * tr.stats.delta < 3600:
            continue
        tr = tr.copy()
        tr.data = tr.data.astype(np.float64)
        tr.detrend('demean'); tr.detrend('linear'); tr.taper(0.05)
        sr = tr.stats.sampling_rate
        fac = int(round(sr))
        for f in ([8, 5] if fac == 40 else [10, 10] if fac == 100 else [5, 4] if fac == 20 else _stages(fac, sr)):
            tr.decimate(f, no_filter=False, strict_length=False)
        try:
            tr.remove_response(inventory=inv, output='DISP', pre_filt=PRE_FILT, zero_mean=True, taper=True)
        except Exception:
            continue
        tr.data[:EDGE] = np.nan; tr.data[-EDGE:] = np.nan
        seg_out += tr
    if len(seg_out) == 0:
        return 'noseg'
    # place each component on the day grid (by original channel), then rotate per sample block
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
    # rotation to N, E, Z from StationXML azimuth/dip
    vecs = []
    for ch in chans:
        try:
            o = inv.get_orientation(f'{netsta}.{seg_out.select(channel=ch)[0].stats.location}.{ch}', t0 + 43200)
        except Exception:
            return f'noorient {ch}'
        az, dip = np.radians(o['azimuth']), np.radians(o['dip'])
        vecs.append([np.cos(dip) * np.cos(az), np.cos(dip) * np.sin(az), -np.sin(dip)])  # N, E, Up
    A = np.array(vecs)
    X = np.vstack([grid[c] for c in chans])
    ok = np.all(np.isfinite(X), axis=0)
    if len(chans) == 3:
        Ainv = np.linalg.inv(A)
        Y = np.full((3, 86400), np.nan)
        Y[:, ok] = Ainv @ X[:, ok]
        D[:] = Y
    else:                                        # Z only
        zi = [i for i, c in enumerate(chans) if c.endswith('Z')]
        if not zi:
            return 'noZ'
        D[2] = X[zi[0]] * np.sign(A[zi[0], 2])
    os.makedirs(os.path.dirname(out), exist_ok=True)
    tmp = out + '.part'                          # atomic write (2026-10-14): a killed worker leaves no half file
    with open(tmp, 'wb') as fh:
        np.savez_compressed(fh, d=D, chans=','.join(chans), band=band)
    os.replace(tmp, out)
    return f'ok {band} {len(chans)}ch {np.isfinite(D[2]).mean():.2f}'


def _job(a):
    try:
        return a, process(*a)
    except Exception as e:  # noqa
        return a, f'ERR {type(e).__name__}: {str(e)[:80]}'


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('stations'); ap.add_argument('start'); ap.add_argument('end')
    ap.add_argument('--workers', type=int, default=8)
    ap.add_argument('--days-file', default='', help='CSV with a date column; overrides START/END')
    a = ap.parse_args()
    days = (pd.to_datetime(pd.read_csv(a.days_file).date) if a.days_file else pd.date_range(a.start, a.end, freq='D'))
    jobs = [(s, d) for s in a.stations.split(',') for d in days]
    from multiprocessing import Pool
    lf = f'{C.LOGD}/prep.log'
    with Pool(a.workers) as p:
        for (s, d), r in p.imap_unordered(_job, jobs):
            C.log(f'{s} {d:%Y-%m-%d} {r}', lf)
