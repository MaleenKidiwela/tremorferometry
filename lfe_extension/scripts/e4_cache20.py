#!/usr/bin/env python
"""E4 — cache 3-C, 1-8 Hz, 20 Hz data as .npy so the GPU correlator never touches obspy.

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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", required=True)          # NET.STA,NET.STA,...
    ap.add_argument("--year", type=int, default=2010)
    ap.add_argument("--outdir", default=f"{ROOT}/cache20hz")
    ap.add_argument("--one-c", action="store_true")
    ap.add_argument("--traces-root", default=ROOT, help="lfe_extension: read <root>/traces/<NET.STA>/<year>/*.mseed")
    a = ap.parse_args()
    for ns in a.stations.split(","):
        days = sorted(int(os.path.basename(f)[:3])
                      for f in glob.glob(f"{a.traces_root}/traces/{ns}/{a.year}/*.mseed"))
        ok = 0
        first = None
        for d in days:
            r = cache_day(ns, a.year, d, a.outdir, one_c=a.one_c, root=a.traces_root)
            if r:
                ok += 1
                first = first or r
        print(f"[cache] {ns} {a.year}: {ok}/{len(days)} days, shape {first[0] if first else None}, "
              f"order {first[1] if first else None}", flush=True)


if __name__ == "__main__":
    main()
