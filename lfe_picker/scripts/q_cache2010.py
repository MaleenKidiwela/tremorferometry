#!/usr/bin/env python
"""Build the 20 Hz cache for the 2010 PB boreholes that the picker's stores do not have.

The pooled 2010 + 2012 repeat of Step 1a needs the SAME 13-station discovery network in both years,
but the picker's memory-mapped day store (results/mm, built by p_build_mm.py) carries no B001-B007
for 2010 — those station-days were never cached, although the raw day traces DO exist read-only in
../lfe_extension/traces (B001, B003, B004, B005, B006, B007, doy 222-256).

This calls scripts/o_cache_month.cache_one AS A LIBRARY, with explicit jobs, so the frozen caching
chain (merge -> component order -> Butterworth-4 1-8 Hz at the native rate -> integer decimation to
20 Hz) is byte-for-byte the one every other year went through.  Output goes to the picker's own
writable cache, lfe_picker/cache20hz/.

    python scripts/q_cache2010.py --workers 12
"""
from __future__ import annotations
import argparse, glob, os, sys, time
import multiprocessing as mp
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import o_common as C
import o_cache_month as CM

WANT = ["PB.B001", "PB.B003", "PB.B004", "PB.B005", "PB.B006", "PB.B007"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", type=int, default=2010)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--stations", default=",".join(WANT))
    a = ap.parse_args()
    jobs = []
    for netsta in a.stations.split(","):
        for p in sorted(glob.glob(f"{C.TRACES}/{netsta}/{a.year}/*.mseed")):
            if os.path.getsize(p) == 0:
                continue
            doy = int(os.path.basename(p)[:3])
            dst = f"{C.CACHE}/{netsta}/{a.year}/{doy:03d}.npz"
            if os.path.exists(dst):
                continue
            jobs.append((netsta, a.year, doy, p, C.CACHE, False, np.nan))
    print(f"[q_cache2010] {len(jobs)} station-days to build -> {C.CACHE}", flush=True)
    t0 = time.time()
    rows = []
    if jobs:
        with mp.Pool(a.workers, maxtasksperchild=20) as pool:
            for i, r in enumerate(pool.imap_unordered(CM.cache_one, jobs, chunksize=1), 1):
                rows.append(r)
                if i % 25 == 0 or i == len(jobs):
                    print(f"  {i}/{len(jobs)}  {time.time()-t0:.0f}s", flush=True)
    if rows:
        D = pd.DataFrame(rows)
        D.to_parquet(f"{C.PRES}/p3/q_cache2010_log.parquet", index=False)
        print(D.groupby(["netsta", "status"]).size().to_string())
        for r in D[D.status == "error"].itertuples():
            print(f"  ERROR {r.netsta} {r.year}.{r.doy:03d}: {r.err}")
    print(f"[q_cache2010] done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
