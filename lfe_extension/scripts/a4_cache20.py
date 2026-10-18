#!/usr/bin/env python
"""a4 (Track A-2 copy of e4_cache20.py) — cache 3-C, 1-8 Hz, 20 Hz data as .npy so the GPU correlator never touches obspy.

Component order is resolved EXPLICITLY (Z, N|1, E|2) by e4_chanorder -- see that module for the
bug this replaces.  One file per station-day: cache20hz/<NET.STA>/<year>/<doy>.npz with
  x  (3, N) float32   t0 (scalar float64, epoch seconds)   order (3,) str
"""
from __future__ import annotations
import os, sys, glob, argparse
import numpy as np
from obspy import read
from scipy.signal import butter, sosfiltfilt, resample_poly
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e4_chanorder as chanorder

ROOT = "/home/jovyan/tremorferometry/lfe_catalog_repair"
FS = 20.0


def bp(x, fs, f1=1.0, f2=8.0):
    sos = butter(4, [f1 / (fs / 2), f2 / (fs / 2)], btype="band", output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def cache_day(netsta, year, doy, outdir, one_c=False, root=None):
    p = f"{root or ROOT}/traces/{netsta}/{year}/{doy:03d}.mseed"
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return None
    st = read(p)
    # Track A-2 fix (2026-09-16): CN day-files in 2021 carry per-FRAGMENT sampling rates that differ
    # in the 7th digit (99.99999237, 100.00000763, ... instead of 100.0).  obspy then refuses to
    # merge them, and the integer-decimation check below rejects the day.  The relative time-base
    # difference is ~1e-7 (0.009 s over a whole day, 1/5 of a 20 Hz sample), so snapping every trace
    # to the nearest integer rate is harmless and is the only way these days enter the chain.
    for tr in st:
        r = float(tr.stats.sampling_rate)
        if r > 0 and abs(r - round(r)) / r < 1e-4 and r != round(r):
            tr.stats.sampling_rate = float(round(r))
    st.merge(method=1, fill_value=0)
    chans = sorted({tr.stats.channel for tr in st})
    order = chanorder.resolve_order(chans, station=f"{netsta}/{year}.{doy}", strict=not one_c)
    keep = [c for c in order if c is not None]
    fs0 = float(st[0].stats.sampling_rate)
    npts = min(tr.stats.npts for tr in st.select(channel=keep[0]))
    npts = min(min(tr.stats.npts for tr in st.select(channel=c)) for c in keep)
    X = np.vstack([st.select(channel=c)[0].data[:npts].astype(np.float64) for c in keep])
    X = bp(X, fs0)
    q = fs0 / FS
    if abs(q - round(q)) > 1e-9:
        raise RuntimeError(f"{netsta} fs={fs0} not an integer multiple of {FS}")
    q = int(round(q))
    if q > 1:
        X = resample_poly(X, 1, q, axis=-1)
    t0 = float(st.select(channel=keep[0])[0].stats.starttime.timestamp)
    od = f"{outdir}/{netsta}/{year}"
    os.makedirs(od, exist_ok=True)
    np.savez(f"{od}/{doy:03d}.npz", x=X.astype(np.float32), t0=np.float64(t0),
             order=np.array(keep), fs=np.float64(FS))
    return X.shape, keep


def _job(t):
    ns, year, d, outdir, one_c, root = t
    try:
        r = cache_day(ns, year, d, outdir, one_c=one_c, root=root)
        return (ns, year, d, "ok" if r else "skip", None if not r else list(r[1]))
    except Exception as e:                      # never let one bad day kill the pass
        return (ns, year, d, f"ERR {type(e).__name__}: {e}", None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", required=True)          # NET.STA,NET.STA,...
    ap.add_argument("--year", type=int, default=2010)
    ap.add_argument("--outdir", default=f"{ROOT}/cache20hz")
    ap.add_argument("--one-c", action="store_true")
    ap.add_argument("--traces-root", default=ROOT, help="lfe_extension: read <root>/traces/<NET.STA>/<year>/*.mseed")
    ap.add_argument("--doy-from", type=int, default=0)
    ap.add_argument("--doy-to", type=int, default=999)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--skip-existing", action="store_true")
    a = ap.parse_args()
    from multiprocessing import Pool
    jobs = []
    for ns in a.stations.split(","):
        days = sorted(int(os.path.basename(f)[:3])
                      for f in glob.glob(f"{a.traces_root}/traces/{ns}/{a.year}/*.mseed"))
        days = [d for d in days if a.doy_from <= d <= a.doy_to]
        for d in days:
            if a.skip_existing and os.path.exists(f"{a.outdir}/{ns}/{a.year}/{d:03d}.npz"):
                continue
            jobs.append((ns, a.year, d, a.outdir, a.one_c, a.traces_root))
    print(f"[cache] {len(jobs)} station-days, {a.workers} workers -> {a.outdir}", flush=True)
    res = []
    with Pool(a.workers) as pool:
        for r in pool.imap_unordered(_job, jobs, chunksize=1):
            res.append(r)
            if r[3].startswith("ERR"):
                print(f"  [err] {r[0]} {r[1]}.{r[2]:03d}: {r[3]}", flush=True)
    import collections
    for ns in sorted({r[0] for r in res}):
        rr = [r for r in res if r[0] == ns]
        ok = sum(1 for r in rr if r[3] == "ok")
        orders = {tuple(r[4]) for r in rr if r[4]}
        print(f"[cache] {ns} {a.year}: {ok}/{len(rr)} days, orders {orders}", flush=True)


if __name__ == "__main__":
    main()
