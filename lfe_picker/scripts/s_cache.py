#!/usr/bin/env python
"""P3D: build the frozen 20 Hz / 1-8 Hz cache for every delivered P3D station-day.

Reads the downloader's manifest READ-ONLY and caches any delivered station-day that has no npz yet,
by calling `r_cache_p0.dispatch` (which is `o_cache_month.cache_one` for 3-C stations and the
Z-only twin for STW/UWFH) -- the frozen chain, unmodified.  Writes only cache20hz/ and its own log.

  python scripts/s_cache.py [--workers 14] [--min-cov 0.5] [--stations CN.SYMB,...]
"""
from __future__ import annotations
import os, sys, csv, time, argparse, warnings
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r_cache_p0 as RC                                                    # noqa: E402
import s_common as S                                                       # noqa: E402

TRACES = f"{S.PICK}/traces"
CACHE = f"{S.PICK}/cache20hz"
MANIFEST = f"{S.P3D}/download_manifest.csv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--min-cov", type=float, default=0.5)
    ap.add_argument("--stations", default="")
    ap.add_argument("--out", default=f"{S.P3D}/s_cache_log.parquet")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    want = set(a.stations.split(",")) if a.stations else None
    jobs = []
    with open(MANIFEST) as f:
        for r in csv.DictReader(f):
            if r.get("status") != "data":
                continue
            netsta = f"{r['net']}.{r['sta']}"
            if want and netsta not in want:
                continue
            try:
                cov = float(r["coverage_frac"] or 0)
            except ValueError:
                continue
            if cov < a.min_cov:
                continue
            y, doy = int(r["date"][:4]), int(r["doy"])
            if y in (2011, 2013):
                continue                                    # sealed referee years
            dst = f"{CACHE}/{netsta}/{y}/{doy:03d}.npz"
            if os.path.exists(dst):
                continue
            src = f"{TRACES}/{netsta}/{y}/{doy:03d}.mseed"
            if not os.path.exists(src):
                continue
            jobs.append((netsta, y, doy, src, CACHE, False, cov))
    jobs = sorted(set(jobs))
    print(f"[s_cache] {len(jobs)} station-days to cache with {a.workers} workers", flush=True)
    if not jobs:
        return
    import multiprocessing as mp
    rows, t0 = [], time.time()
    with mp.Pool(a.workers, maxtasksperchild=40) as pool:
        for i, r in enumerate(pool.imap_unordered(RC.dispatch, jobs, chunksize=1), 1):
            rows.append(r)
            if i % 250 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)}  {time.time()-t0:.0f}s", flush=True)
    D = pd.DataFrame(rows)
    if os.path.exists(a.out):
        D = pd.concat([pd.read_parquet(a.out), D], ignore_index=True)
    D.to_parquet(a.out, index=False)
    print(D.groupby(["netsta", "status"]).size().to_string())
    bad = D[D.status == "error"]
    for r in bad.head(10).itertuples():
        print(f"  ERROR {r.netsta} {r.year}.{r.doy:03d}: {r.err}")


if __name__ == "__main__":
    main()
