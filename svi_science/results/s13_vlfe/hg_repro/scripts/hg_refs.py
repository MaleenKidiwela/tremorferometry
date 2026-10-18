#!/usr/bin/env python
"""Far-reference vertical displacement, 0.02-0.05 Hz, on the hg_repro time grid (PLAN.md row A3), read-only from
../data/raw_ref2011 and raw_ref2012 (LHZ, 1 sps; IU.COR, US.NEW, US.HAWA, UW.OMAK, TA.I02D; 300-700 km) with
responses ../data/resp_scan/NET.STA[.2012].xml. Per day (+-1800 s padding from neighbours): merge, response -> DISP
(pre_filt 0.004/0.008/0.25/0.45 Hz as the prep), placed on the grid; then hg_cont.filt (same band-pass, 120-s tapers).
Grid covers 1 Jun 2011 - 31 Dec 2012 (no 2013 references). Output data/mf/ref_STA.npy (N,) float32, NaN = dead.
usage: hg_refs.py"""
import os, sys
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import hg_common as H
from hg_cont import filt

REFS = ['IU.COR', 'US.NEW', 'US.HAWA', 'UW.OMAK', 'TA.I02D']
PRE = (0.004, 0.008, 0.25, 0.45)


def dayfile(sta, d):
    p = f'{H.S13}/data/raw_ref{d.year}/{sta}/{d.year}/{d.dayofyear:03d}.mseed'
    return p if os.path.isfile(p) else None


def build(sta):
    from obspy import read, read_inventory, UTCDateTime, Stream
    out = f'{H.MF}/ref_{sta}.npy'
    if os.path.isfile(out):
        return sta, 'exists'
    inv = read_inventory(f'{H.S13}/data/resp_scan/{sta}.xml')
    f2 = f'{H.S13}/data/resp_scan/{sta}.2012.xml'
    if os.path.isfile(f2):
        inv += read_inventory(f2)
    X = np.full(H.N, np.nan, np.float32)
    nd = 0
    for k in range((pd.Timestamp('2013-01-01') - H.T0).days):
        d = H.T0 + pd.Timedelta(days=k)
        if dayfile(sta, d) is None:
            continue
        t0 = UTCDateTime(str(d.date()))
        st = Stream()
        for dd, a, b in ((d - pd.Timedelta(days=1), t0 - 1800, t0), (d, t0, t0 + 86400),
                         (d + pd.Timedelta(days=1), t0 + 86400, t0 + 86400 + 1800)):
            p = dayfile(sta, dd)
            if p:
                try:
                    st += read(p, starttime=a, endtime=b)
                except Exception:
                    pass
        st = st.select(channel='LHZ')
        st.merge(method=1, fill_value=None)
        st = st.split()
        for tr in st:
            if tr.stats.npts < 3600:
                continue
            tr.data = tr.data.astype(np.float64)
            tr.detrend('demean'); tr.detrend('linear')
            tr.taper(max_percentage=0.5, type='cosine', max_length=600, side='both')
            try:
                tr.remove_response(inventory=inv, output='DISP', pre_filt=PRE, zero_mean=True, taper=False)
            except Exception:
                continue
            tr.data[:600] = np.nan; tr.data[-600:] = np.nan
            i0 = int(round(tr.stats.starttime - t0))
            a, b = max(i0, 0), min(i0 + tr.stats.npts, 86400)
            if b > a:
                X[k * 86400 + a:k * 86400 + b] = tr.data[a - i0:b - i0]
        nd += 1
    X = filt(X)
    np.save(out + '.part.npy', X); os.replace(out + '.part.npy', out)
    return sta, f'days {nd} live {np.isfinite(X).mean():.3f}'


if __name__ == '__main__':
    from multiprocessing import Pool
    with Pool(3) as p:
        for s, r in p.imap_unordered(build, REFS):
            line = f'{pd.Timestamp.now(tz="America/Los_Angeles"):%Y-%m-%d %H:%M:%S} PT ref {s} {r}'
            print(line, flush=True)
            open(f'{H.LOGD}/refs.log', 'a').write(line + '\n')
