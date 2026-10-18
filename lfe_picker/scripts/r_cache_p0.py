#!/usr/bin/env python
"""P3D step 5 - cache the P0 (SSE-window) station-days to the FROZEN 20 Hz / 1-8 Hz chain.

The chain itself is NOT reimplemented: 3-C station-days go through `o_cache_month.cache_one`
verbatim (merge(method=1, fill_value=0) -> component order via o_chanorder -> Butterworth-4
zero-phase 1-8 Hz at the native rate -> integer decimation to 20 Hz with resample_poly).

The only addition is the two Z-ONLY helpers (UW.STW, UW.UWFH, EHZ): `cache_one` resolves the
component order with strict=True and would raise on a single-channel day, so those days are run
through `cache_one_z` below, which calls the SAME o_cache_month.bp / med_win_energy / decimation
with strict=False and stores a (1, N) array with order ["EHZ"].  Nothing about the 3-C path
changes.

Output: cache20hz/<NET.STA>/<YYYY>/<DOY>.npz   x (n,N) float32, t0 float64, order (n,), fs 20.0
Log:    results/p3d/cache_p0.parquet (one row per station-day)

  python scripts/r_cache_p0.py [--priority P0] [--workers 12] [--force] [--stations CN.SYMB,...]
"""
from __future__ import annotations
import os
import sys
import csv
import time
import argparse
import warnings
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))
import o_cache_month as OC                                            # noqa: E402
import o_chanorder as chanorder                                       # noqa: E402

RESULTS = os.path.join(HERE, "results", "p3d")
TRACES = os.path.join(HERE, "traces")
CACHE = os.path.join(HERE, "cache20hz")
MANIFEST = os.path.join(RESULTS, "download_manifest.csv")
QUEUE = os.path.join(RESULTS, "download_queue.csv")
FS = OC.FS
WIN_S = OC.WIN_S


def cache_one_z(job):
    """Single-component (EHZ) twin of o_cache_month.cache_one - same filter, same decimation."""
    netsta, year, doy, src, out, force, cov_man = job
    t0w = time.time()
    dst = f"{out}/{netsta}/{year}/{doy:03d}.npz"
    try:
        if os.path.exists(dst) and not force:
            z = np.load(dst)
            X = z["x"].astype(np.float64)
            me = OC.med_win_energy(X, int(round(WIN_S * FS)))
            return dict(netsta=netsta, year=year, doy=doy, status="cached", n=X.shape[1],
                        t0=float(z["t0"]), order=",".join(list(z["order"])),
                        e8_Z=me[0], e8_N=np.nan, e8_E=np.nan,
                        frac_nonzero=float((np.abs(X).sum(0) > 0).mean()),
                        cov_manifest=cov_man, secs=time.time() - t0w, err="")
        from obspy import read
        st = read(src)
        n_snap = 0
        for tr in st:
            fs_ = float(tr.stats.sampling_rate)
            fr = round(fs_)
            if fr > 0 and fs_ != fr and abs(fs_ - fr) / fr < 1e-4:
                tr.stats.sampling_rate = float(fr)
                n_snap += 1
        st.merge(method=1, fill_value=0)
        chans = sorted({tr.stats.channel for tr in st})
        order = chanorder.resolve_order(chans, station=f"{netsta}/{year}.{doy}", strict=False)
        keep = [c for c in order if c is not None]
        if not keep or not keep[0].endswith("Z"):
            raise RuntimeError(f"{netsta}: no Z channel in {chans}")
        keep = keep[:1]                                     # Z-only helper: vertical only
        fs0 = float(st.select(channel=keep[0])[0].stats.sampling_rate)
        npts = min(tr.stats.npts for tr in st.select(channel=keep[0]))
        X = np.vstack([st.select(channel=keep[0])[0].data[:npts].astype(np.float64)])
        X = OC.bp(X, fs0)
        q = fs0 / FS
        if abs(q - round(q)) > 1e-9:
            raise RuntimeError(f"{netsta} fs={fs0} not an integer multiple of {FS}")
        q = int(round(q))
        if q > 1:
            X = OC.resample_poly(X, 1, q, axis=-1)
        t0 = float(st.select(channel=keep[0])[0].stats.starttime.timestamp)
        os.makedirs(f"{out}/{netsta}/{year}", exist_ok=True)
        tmp = dst + f".tmp{os.getpid()}.npz"
        np.savez(tmp, x=X.astype(np.float32), t0=np.float64(t0), order=np.array(keep),
                 fs=np.float64(FS))
        os.replace(tmp, dst)
        me = OC.med_win_energy(X, int(round(WIN_S * FS)))
        return dict(netsta=netsta, year=year, doy=doy, status="built", n=X.shape[1], t0=t0,
                    order=",".join(keep), fs_native=fs0, n_fs_snapped=n_snap,
                    e8_Z=me[0], e8_N=np.nan, e8_E=np.nan,
                    frac_nonzero=float((np.abs(X).sum(0) > 0).mean()), cov_manifest=cov_man,
                    secs=time.time() - t0w, err="")
    except Exception as e:
        return dict(netsta=netsta, year=year, doy=doy, status="error", n=0, t0=np.nan,
                    order="", e8_Z=np.nan, e8_N=np.nan, e8_E=np.nan, frac_nonzero=np.nan,
                    cov_manifest=cov_man, secs=time.time() - t0w, err=f"{type(e).__name__}: {e}")


ZONLY = {"UW.STW", "UW.UWFH"}


def dispatch(job):
    return cache_one_z(job) if job[0] in ZONLY else OC.cache_one(job)


def build_jobs(prio, force, stations, min_cov):
    prios = set(prio.split(","))
    want = set()
    with open(QUEUE) as f:
        for r in csv.DictReader(f):
            if r["priority"] in prios and r["status"] == "pending":
                want.add((r["sta"], r["date"]))
    jobs, skipped = [], []
    with open(MANIFEST) as f:
        for r in csv.DictReader(f):
            if r.get("status") != "data" or (r["sta"], r["date"]) not in want:
                continue
            netsta = f"{r['net']}.{r['sta']}"
            if stations and netsta not in stations:
                continue
            cov = float(r["coverage_frac"] or 0)
            if cov < min_cov:
                skipped.append((netsta, r["date"], f"coverage {cov:.3f} < {min_cov}"))
                continue
            y, doy = int(r["date"][:4]), int(r["doy"])
            src = os.path.join(TRACES, netsta, str(y), f"{doy:03d}.mseed")
            if not os.path.exists(src):
                skipped.append((netsta, r["date"], "no trace file"))
                continue
            jobs.append((netsta, y, doy, src, CACHE, force, cov))
    jobs = sorted(set(jobs))
    return jobs, skipped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--priority", default="P0")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--min-cov", type=float, default=0.5)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--stations", default="")
    ap.add_argument("--out", default=os.path.join(RESULTS, "cache_p0.parquet"))
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    stations = set(a.stations.split(",")) if a.stations else None
    jobs, skipped = build_jobs(a.priority, a.force, stations, a.min_cov)
    print(f"[cache {a.priority}] {len(jobs)} station-days to cache "
          f"({len(skipped)} skipped), {a.workers} workers", flush=True)
    t0 = time.time()
    rows = []
    if jobs:
        import multiprocessing as mp
        with mp.Pool(a.workers, maxtasksperchild=40) as pool:
            for i, r in enumerate(pool.imap_unordered(dispatch, jobs, chunksize=1), 1):
                rows.append(r)
                if i % 200 == 0 or i == len(jobs):
                    print(f"  {i}/{len(jobs)} {time.time()-t0:.0f}s", flush=True)
    D = pd.DataFrame(rows)
    if len(D):
        D["date"] = (pd.to_datetime(D.year.astype(str), format="%Y")
                     + pd.to_timedelta(D.doy - 1, "D")).dt.strftime("%Y-%m-%d")
        D = D.sort_values(["netsta", "date"])
        D.to_parquet(a.out, index=False)
        print(D.groupby(["netsta", "status"]).size().to_string())
        bad = D[D.status == "error"]
        for r in bad.head(20).itertuples():
            print(f"  ERROR {r.netsta} {r.year}.{r.doy:03d}: {r.err}", flush=True)
        print(f"[cache {a.priority}] {len(D)} days in {time.time()-t0:.0f}s -> {a.out}")
    for s in skipped[:15]:
        print(f"  SKIP {s[0]} {s[1]}: {s[2]}")
    if len(skipped) > 15:
        print(f"  ... {len(skipped)-15} more skipped")


if __name__ == "__main__":
    main()
