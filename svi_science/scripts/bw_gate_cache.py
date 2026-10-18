#!/usr/bin/env python
"""bw_gate_cache -- 20 Hz / 1-8 Hz caches of the backward raw (lfe_extension/traces_backward) with the FORWARD chain:
lfe_extension/scripts/a4_cache20.py cache_day() line for line (integer-rate snap, merge fill 0, explicit Z/N|1/E|2 order
by e4_chanorder, 4-pole 1-8 Hz zero-phase band-pass at the native rate, resample_poly to 20 Hz), only the raw path is
a parameter (a4_cache20 hard-wires <root>/traces/).  Output: lfe_extension/cache20hz_backward/<NET.STA>/<YYYY>/<DOY>.npz.
Identity is checked by --check: re-cache a day that already has a forward cache and compare arrays.

  python bw_gate_cache.py --stations CN.PGC,PB.B009 --start 2008-05-01 --end 2008-06-04 [--workers 10]
"""
from __future__ import annotations
import os, sys, glob, argparse
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bw_gate_common as G
sys.path.insert(0, G.LFE_SCRIPTS)
import a4_cache20 as A
from obspy import read
from scipy.signal import resample_poly


def cache_day(p, outdir, netsta, year, doy):
    st = read(p)
    for tr in st:
        r = float(tr.stats.sampling_rate)
        if r > 0 and abs(r - round(r)) / r < 1e-4 and r != round(r):
            tr.stats.sampling_rate = float(round(r))
    st.merge(method=1, fill_value=0)
    chans = sorted({tr.stats.channel for tr in st})
    order = A.chanorder.resolve_order(chans, station=f"{netsta}/{year}.{doy}", strict=True)
    keep = [c for c in order if c is not None]
    fs0 = float(st[0].stats.sampling_rate)
    npts = min(min(tr.stats.npts for tr in st.select(channel=c)) for c in keep)
    X = np.vstack([st.select(channel=c)[0].data[:npts].astype(np.float64) for c in keep])
    X = A.bp(X, fs0)
    q = fs0 / A.FS
    if abs(q - round(q)) > 1e-9:
        raise RuntimeError(f"{netsta} fs={fs0} not an integer multiple of {A.FS}")
    q = int(round(q))
    if q > 1:
        X = resample_poly(X, 1, q, axis=-1)
    t0 = float(st.select(channel=keep[0])[0].stats.starttime.timestamp)
    od = f"{outdir}/{netsta}/{year}"
    os.makedirs(od, exist_ok=True)
    np.savez(f"{od}/{doy:03d}.npz", x=X.astype(np.float32), t0=np.float64(t0), order=np.array(keep), fs=np.float64(A.FS))
    return X.shape, keep


def _job(t):
    p, outdir, ns, y, d = t
    try:
        r = cache_day(p, outdir, ns, y, d)
        return (ns, y, d, "ok", r[0][1])
    except Exception as e:
        return (ns, y, d, f"ERR {type(e).__name__}: {str(e)[:100]}", 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", required=True); ap.add_argument("--start", required=True); ap.add_argument("--end", required=True)
    ap.add_argument("--root", default=G.TRACES_BW); ap.add_argument("--outdir", default=G.CACHE_BW)
    ap.add_argument("--status-dir", default=f"{G.OUT}/tables")
    ap.add_argument("--workers", type=int, default=10); ap.add_argument("--check", default="",
                    help="forward cache dir holding the same station-days: compare instead of writing")
    a = ap.parse_args()
    jobs = []
    for ns in a.stations.split(","):
        for x in pd.date_range(a.start, a.end, freq="D"):
            y, d = x.year, x.dayofyear
            p = f"{a.root}/{ns}/{y}/{d:03d}.mseed"
            if not os.path.exists(p) or os.path.getsize(p) == 0:
                continue
            if not a.check and os.path.exists(f"{a.outdir}/{ns}/{y}/{d:03d}.npz"):
                continue
            jobs.append((p, a.check and f"{G.OUT}/tmp_cachecheck" or a.outdir, ns, y, d))
    G.log(f"[cache] {len(jobs)} station-days -> {a.outdir}")
    from multiprocessing import Pool
    with Pool(a.workers) as pool:
        res = list(pool.imap_unordered(_job, jobs, chunksize=1))
    bad = [r for r in res if r[3] != "ok"]
    for r in bad:
        G.log(f"  [err] {r}")
    G.log(f"[cache] ok {len(res) - len(bad)} / {len(res)}")
    if a.check:
        for p, od, ns, y, d in jobs:
            A1 = np.load(f"{od}/{ns}/{y}/{d:03d}.npz"); B1 = np.load(f"{a.check}/{ns}/{y}/{d:03d}.npz")
            n = min(A1["x"].shape[1], B1["x"].shape[1])
            G.log(f"[check] {ns} {y}.{d:03d}: shapes {A1['x'].shape} {B1['x'].shape} t0 diff {float(A1['t0'] - B1['t0']):.4f} s "
                  f"max|diff| {float(np.abs(A1['x'][:, :n] - B1['x'][:, :n]).max()):.3g} (max|x| {float(np.abs(B1['x']).max()):.3g})")
    sp_ = f"{a.status_dir}/cache_status_{a.start}_{a.end}.csv"
    pd.DataFrame(res, columns=["netsta", "year", "doy", "status", "npts20"]).to_csv(sp_, mode="a", index=False,
                                                                                    header=not os.path.exists(sp_))


if __name__ == "__main__":
    main()
