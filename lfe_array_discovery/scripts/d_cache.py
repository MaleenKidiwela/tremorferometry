#!/usr/bin/env python
"""20 Hz / 1-8 Hz cache for the station-days of one tranche (FROZEN chain + rational twin).

Integer-decimation stations go through `lfe_picker/scripts/o_cache_month.cache_one` VERBATIM;
the non-integer ones (XG/XU at 50 Hz, GNW BH at 50 Hz, CC at 99.92 Hz) go through
`lfe_picker/scripts/u_cache.cache_one_rational`, the documented twin whose ONLY difference is
resample_poly(up, down) for the exact rational fs/20.  Every row carries `chain`.

  python scripts/d_cache.py T2A [--workers 24] [--min-cov 0.2]
-> cache20hz/<NET.STA>/<YYYY>/<DOY>.npz and results/cache_<TRANCHE>.parquet with the per
   station-day coverage_frac (from the download manifest) and frac_nonzero (measured).
The SEALED tranche is never cached (traces_sealed is not read).
"""
from __future__ import annotations
import os
import sys
import time
import argparse
import warnings
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "lfe_picker", "scripts")))
import d_common as D                                                   # noqa: E402
import o_cache_month as OC                                             # noqa: E402
import u_cache as UC                                                   # noqa: E402


def dispatch(job):
    netsta, year, doy, src, out, force, cov, chain = job
    r = OC.cache_one(job[:7]) if chain == "frozen" else UC.cache_one_rational(job[:7])
    r["chain"] = chain
    return r


def build_jobs(tranche, force, min_cov):
    M = pd.read_csv(D.MANIFEST, keep_default_na=False)
    M = M[(M.tranche == tranche) & (M.status == "data")]
    M = M[M.outroot != D.TRACES_SEALED]
    M = M.drop_duplicates(["net", "sta", "day"], keep="last")
    jobs, skipped = [], []
    for r in M.itertuples():
        netsta = f"{r.net}.{r.sta}"
        try:
            cov = float(r.coverage_frac or 0)
        except (TypeError, ValueError):
            cov = 0.0
        if cov < min_cov:
            skipped.append((netsta, r.day, f"coverage {cov:.3f} < {min_cov}"))
            continue
        d = pd.Timestamp(r.day)
        src = os.path.join(r.outroot, netsta, str(d.year), f"{int(d.dayofyear):03d}.mseed")
        if not os.path.exists(src):
            skipped.append((netsta, r.day, "no trace file"))
            continue
        try:
            fs = float(r.fs or 0)
        except (TypeError, ValueError):
            fs = 0.0
        chain = D.chain_for(fs) if fs > 0 else (r.chain or "frozen")
        jobs.append((netsta, d.year, int(d.dayofyear), src, D.CACHE, force, cov, chain))
    return sorted(set(jobs)), skipped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tranche")
    ap.add_argument("--workers", type=int, default=24)       # pod cgroup cap is 32 CPUs
    ap.add_argument("--min-cov", type=float, default=0.2)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    jobs, skipped = build_jobs(a.tranche, a.force, a.min_cov)
    D.log(f"[cache {a.tranche}] {len(jobs)} station-days ({len(skipped)} skipped), "
          f"{a.workers} workers")
    t0, rows = time.time(), []
    if jobs:
        import multiprocessing as mp
        with mp.Pool(a.workers, maxtasksperchild=20) as pool:
            for i, r in enumerate(pool.imap_unordered(dispatch, jobs, chunksize=1), 1):
                rows.append(r)
                if i % 250 == 0 or i == len(jobs):
                    D.log(f"  {i}/{len(jobs)} {time.time()-t0:.0f}s")
    if rows:
        C = pd.DataFrame(rows)
        C["date"] = (pd.to_datetime(C.year.astype(str), format="%Y")
                     + pd.to_timedelta(C.doy - 1, "D")).dt.strftime("%Y-%m-%d")
        C["coverage_frac"] = C.cov_manifest
        C = C.sort_values(["netsta", "date"])
        C.to_parquet(os.path.join(D.RES, f"cache_{a.tranche}.parquet"), index=False)
        D.log(C.groupby(["chain", "status"]).size().to_string())
        for r in C[C.status == "error"].head(15).itertuples():
            D.log(f"  ERROR {r.netsta} {r.date}: {r.err}")
        if "fs_native" in C:
            D.log(f"  native sample rates seen: "
                  f"{sorted(set(np.round(C.fs_native.dropna().values, 3)))}")
        D.log(f"[cache {a.tranche}] {len(C)} days in {time.time()-t0:.0f}s")
    for s in skipped[:10]:
        D.log(f"  SKIP {s[0]} {s[1]}: {s[2]}")
    if len(skipped) > 10:
        D.log(f"  ... {len(skipped)-10} more skipped")


if __name__ == "__main__":
    main()
