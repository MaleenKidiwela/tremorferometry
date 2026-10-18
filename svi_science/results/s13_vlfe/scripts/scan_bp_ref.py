#!/usr/bin/env python
"""VLFE scan pilot 2 (S7.1): band-passed year, daily sigma and S4 glitch mask for reference stations — scan_bp.py
unchanged, called with station coordinates that include the reference stations (scan_refs.all_coords), so that the
injection callback can put the synthetic (g = 1) into reference data too.
usage: scan_bp_ref.py OUTDIR --stations A,B [--inject inj.parquet] [--workers 8]"""
import os, sys, argparse
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import scan_refs as R
_cc = R.all_coords()
C.station_coords = lambda: _cc                    # noqa: E731
import scan_bp as B

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('out'); ap.add_argument('--inject', default=''); ap.add_argument('--stations', required=True)
    ap.add_argument('--workers', type=int, default=8)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    from multiprocessing import Pool
    with Pool(a.workers) as p:
        for s, r in p.imap_unordered(B.one, [(s, a.out, a.inject) for s in a.stations.split(',')]):
            C.log(f'scan_bp_ref {a.out}: {s} {r}', f'{C.LOGD}/scan2_bp.log')
