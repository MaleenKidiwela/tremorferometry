#!/usr/bin/env python
"""DUC_FIX prep -- per-day burst mute / line notches / dead-channel floors / channel reference energies
(duc_fix_clean.prep_day, the sh_v2 fixes) for every station his families name, on cache_duc.
  python scripts/duc_fix_prep.py --t0 2008-04-01 --t1 2008-05-12 [--workers 3] [--tag repro]
Writes results/duc_fix/prep/<date>.npz and results/duc_fix/prep_summary/<tag>.csv.
"""
import os, sys, time, argparse, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from duc_fix_common import *
import duc_fix_clean as CL

NREF = int(round(20.0 * FS))       # 20 s reference windows (guard scales by L / NREF)
_C = None; _NETOF = None


def _init():
    global _C, _NETOF
    _C = CL.coords()
    fam = families(); used = set(",".join(fam.stations).split(","))
    H = his_rows(); _NETOF = {s: n for s, n in zip(H.sta, H.net) if s in used}


def one(d):
    warnings.filterwarnings("ignore")
    stas = [f"{n}.{s}" for s, n in sorted(_NETOF.items())]
    try:
        return CL.prep_day(d, stas, _C, NREF)
    except Exception as e:
        return [dict(day=str(pd.Timestamp(d).date()), error=str(e)[:200])]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--t0", required=True); ap.add_argument("--t1", required=True)
    ap.add_argument("--workers", type=int, default=3); ap.add_argument("--tag", default="")
    a = ap.parse_args()
    os.makedirs(CL.PREP, exist_ok=True); os.makedirs(f"{FIX}/prep_summary", exist_ok=True)
    days = [d for d in pd.date_range(a.t0, a.t1, freq="D", inclusive="left")
            if not os.path.exists(f"{CL.PREP}/{d.date()}.npz")]
    t0 = time.time(); rows = []
    import multiprocessing as mp
    with mp.Pool(a.workers, initializer=_init, maxtasksperchild=20) as pool:
        for r in pool.imap_unordered(one, days):
            rows += r
    D = pd.DataFrame(rows)
    tag = a.tag or f"{a.t0}_{a.t1}"
    if len(D): D.to_csv(f"{FIX}/prep_summary/{tag}.csv", index=False)
    nerr = int(D["error"].notna().sum()) if "error" in D else 0
    log(f"[prep {tag}] {len(days)} days, {len(D)} station-days, errors {nerr}, "
        f"muted mean {D.muted_frac.mean() if 'muted_frac' in D else float('nan'):.4f}, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
