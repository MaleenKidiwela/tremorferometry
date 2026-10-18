#!/usr/bin/env python
"""VLFE scan pilot 2 (PLAN_SCAN.md S7.1): prep v3 for reference stations. Identical recipe to scan_prep_v3.py (imported,
unchanged); this wrapper only (1) looks up the response per station in the S7.1 order (data/resp/{NET}_bb.xml ->
lfe_array_discovery/results/stationxml/{NET}.xml -> data/resp_scan/{NET}.xml) and (2) accepts LH? (1-sps) channels when
no BH/HH exist (for the pending far references; no decimation needed).
usage: scan_prep_ref.py NET.STA[,..] START END [--workers 8]   (writes data/disp1sps_v3/NET.STA/YYYY/DOY.npz)"""
import os, sys, argparse
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('S13_DISP', '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/data/disp1sps_v3')
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import scan_prep_v3 as V3
import scan_refs as R

P = V3.P
C.TRACE_ROOTS.append('svi_science/results/s13_vlfe/data/raw_ref2011')     # far references (S7.13)
C.TRACE_ROOTS.append('svi_science/results/s13_vlfe/data/raw_ref2012')     # far references 2012 (S8.2)
C.TRACE_ROOTS.append('svi_science/results/s13_vlfe/data/raw_ref_all')    # far references 2003-2026 (catalogue Amendment P)


def pick_band_lh(st):
    for band in ('BH', 'HH', 'LH'):
        s = st.select(channel=f'{band}?')
        comps = {tr.stats.channel[-1] for tr in s}
        if len(comps) >= 3 or ('Z' in comps and len(comps) >= 1):
            return s, band
    return None, None


def process_ref(netsta, day):
    from obspy import UTCDateTime
    inv, f = R.inventory_for(netsta, UTCDateTime(str(pd.Timestamp(day).date())) + 43200)
    P.inventory = lambda net: inv              # noqa: E731  (per-call station inventory)
    P.pick_band = pick_band_lh
    return V3.process(netsta, day)


def _job(a):
    try:
        return a, process_ref(*a)
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
    lf = f'{C.LOGD}/scan2_prep_ref.log'
    with Pool(a.workers) as p:
        for (s, d), r in p.imap_unordered(_job, jobs):
            C.log(f'{s} {d:%Y-%m-%d} {r}', lf)
