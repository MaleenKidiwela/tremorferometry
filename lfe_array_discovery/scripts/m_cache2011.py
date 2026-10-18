#!/usr/bin/env python
"""Cache the SEALED 2011 station-days with the FROZEN chain.  Prefix `m_`.

`scripts/d_cache.py` verbatim in everything that touches a waveform -- integer-decimation
stations go through `lfe_picker/scripts/o_cache_month.cache_one`, the non-integer ones (XG/XU at
50 Hz) through `lfe_picker/scripts/u_cache.cache_one_rational`, the documented twin whose only
difference is resample_poly(up, down) -- with ONE change: d_cache refuses the sealed tranche
(`M.outroot != D.TRACES_SEALED`), which is exactly the line the seal exists to enforce.  The seal
has been opened and recorded (results/SEALED_2011_OPENED_*), so that line is lifted here and
nowhere else.  traces_sealed/ is only ever READ.

Output: cache20hz/<NET.STA>/2011/<DOY>.npz  +  results/cache_T3.parquet (per station-day
coverage_frac from the download manifest and the measured frac_nonzero).

  python scripts/m_cache2011.py [--workers 24] [--min-cov 0.2] [--days 2011-08-15,...]
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import d_common as D                                                        # noqa: E402
import d_cache as DC                                                        # noqa: E402

TRANCHE = "T3"


def build_jobs(force, min_cov, days=None):
    """d_cache.build_jobs for the sealed tranche (the `outroot != TRACES_SEALED` line lifted)."""
    M = pd.read_csv(D.MANIFEST, keep_default_na=False)
    M = M[(M.tranche == TRANCHE) & (M.status == "data")]
    M = M.drop_duplicates(["net", "sta", "day"], keep="last")
    if days:
        M = M[M.day.isin(set(days))]
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
    ap.add_argument("--workers", type=int, default=24)       # pod cgroup cap is 32 CPUs
    ap.add_argument("--min-cov", type=float, default=0.2)
    ap.add_argument("--days", default="")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    import glob
    if not glob.glob(f"{D.RES}/SEALED_2011_OPENED_*"):
        raise SystemExit("the 2011 seal has not been opened/recorded -- run m_seal_open.py first")
    days = [x for x in a.days.split(",") if x] or None
    jobs, skipped = build_jobs(a.force, a.min_cov, days)
    D.log(f"[cache {TRANCHE}] {len(jobs)} station-days ({len(skipped)} skipped), "
          f"{a.workers} workers")
    t0, rows = time.time(), []
    if jobs:
        import multiprocessing as mp
        with mp.Pool(a.workers, maxtasksperchild=20) as pool:
            for i, r in enumerate(pool.imap_unordered(DC.dispatch, jobs, chunksize=1), 1):
                rows.append(r)
                if i % 250 == 0 or i == len(jobs):
                    D.log(f"  {i}/{len(jobs)} {time.time()-t0:.0f}s")
    if rows:
        C = pd.DataFrame(rows)
        C["date"] = (pd.to_datetime(C.year.astype(str), format="%Y")
                     + pd.to_timedelta(C.doy - 1, "D")).dt.strftime("%Y-%m-%d")
        C["coverage_frac"] = C.cov_manifest
        C = C.sort_values(["netsta", "date"])
        out = a.out or os.path.join(D.RES, f"cache_{TRANCHE}.parquet")
        if os.path.exists(out) and days:
            P = pd.read_parquet(out)
            C = (pd.concat([P, C], ignore_index=True)
                   .drop_duplicates(["netsta", "date"], keep="last")
                   .sort_values(["netsta", "date"]))
        C.to_parquet(out, index=False)
        D.log(C.groupby(["chain", "status"]).size().to_string())
        for r in C[C.status == "error"].head(15).itertuples():
            D.log(f"  ERROR {r.netsta} {r.date}: {r.err}")
        if "fs_native" in C:
            D.log("  native sample rates seen: "
                  f"{sorted(set(np.round(C.fs_native.dropna().values, 3)))}")
        D.log(f"[cache {TRANCHE}] {len(C)} rows in {time.time()-t0:.0f}s -> {out}")
    for s in skipped[:10]:
        D.log(f"  SKIP {s[0]} {s[1]}: {s[2]}")
    if len(skipped) > 10:
        D.log(f"  ... {len(skipped)-10} more skipped")


if __name__ == "__main__":
    main()
