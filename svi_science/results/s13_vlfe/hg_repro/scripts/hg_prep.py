#!/usr/bin/env python
"""1-sps displacement day files for the H&G (2019) stations (PLAN.md step 1). Same recipe as ../scripts/scan_prep_v3.py
(imported read-only; its process() is called unchanged) with three in-process patches (nothing in s13_vlfe/ is modified):
  - raw roots: s13_common.TRACE_ROOTS + hg_repro/data/raw (fetched days);
  - band choice: BH, HH, then LH (TA 1-sps days); a 1-sps trace is not decimated;
  - responses: ../data/resp/{NET}_bb.xml merged with hg_repro/data/resp/{NET}.*.xml (TA.A04D, LH channels).
Skips a day if hg_repro/data/disp1sps or ../data/disp1sps_v3 already holds it (v3 files are read-only shared inputs).
Writes hg_repro/data/disp1sps/NET.STA/YYYY/DOY.npz (d: (3, 86400) float32 N,E,Z [m]; gap_iv; nan_iv; prep='v3').
usage: hg_prep.py NET.STA[,..] START END [--workers 4]   (END inclusive)"""
import os, sys, glob, argparse
os.environ['OMP_NUM_THREADS'] = '1'
HG = '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/hg_repro'
os.environ['S13_DISP'] = f'{HG}/data/disp1sps'
sys.path.insert(0, '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/scripts')
import pandas as pd
import s13_common as C
import s13_prep as P
import scan_prep_v3 as V

C.TRACE_ROOTS.append('svi_science/results/s13_vlfe/hg_repro/data/raw')
V3 = f'{C.DATA}/disp1sps_v3'
_INV = {}


def pick_band(st):
    for band in ('BH', 'HH', 'LH'):
        s = st.select(channel=f'{band}?')
        comps = {tr.stats.channel[-1] for tr in s}
        if len(comps) >= 3 or ('Z' in comps and len(comps) >= 1):
            return s, band
    return None, None


def inventory(net):
    from obspy import read_inventory, Inventory
    if net not in _INV:
        inv = Inventory(networks=[], source='hg_repro')
        f = f'{C.RESP}/{net}_bb.xml'
        if os.path.isfile(f):
            inv += read_inventory(f)
        for g in sorted(glob.glob(f'{HG}/data/resp/{net}.*.xml')):
            inv += read_inventory(g)
        _INV[net] = inv if len(inv.networks) else None
    return _INV[net]


P.pick_band = pick_band
P.inventory = inventory


def job(a):
    s, d = a
    rel = f'{s}/{d.year}/{d.dayofyear:03d}.npz'
    if os.path.isfile(f'{V3}/{rel}') or os.path.isfile(f'{HG}/data/disp1sps/{rel}'):
        return a, 'exists'
    try:
        return a, V.process(s, d)
    except Exception as e:  # noqa
        return a, f'ERR {type(e).__name__}: {str(e)[:80]}'


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('stations'); ap.add_argument('start'); ap.add_argument('end')
    ap.add_argument('--workers', type=int, default=4)
    a = ap.parse_args()
    days = pd.date_range(a.start, a.end, freq='D')
    jobs = [(s, d) for d in days for s in a.stations.split(',')]
    from multiprocessing import Pool
    lf = f'{HG}/logs/prep.log'
    with Pool(min(a.workers, 4)) as p:
        for (s, d), r in p.imap_unordered(job, jobs, chunksize=4):
            C.log(f'{s} {d:%Y-%m-%d} {r}', lf)
