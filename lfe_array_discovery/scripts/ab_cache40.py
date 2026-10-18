#!/usr/bin/env python
"""STAGE 2: the 40 sps / 1-8 Hz cache Royer & Bostock's method needs, WITH A VALIDITY MASK.

sec.3: "All data were divided into 24 hour-long segments, band-pass filtered between 1 and 8 Hz,
resampled to 40 sps".  The existing project cache is 20 Hz, so this is a new cache, not a re-use.

THE MASK IS THE POINT.  The 20 Hz cache is written after `merge(method=1, fill_value=0)` with no
mask (lfe_picker/scripts/o_cache_month.py:70).  Two things then go wrong downstream:
  * the band-pass turns a zero-filled gap into a very small but NOT exactly zero signal, so
    `frac_nonzero` in the cache manifest reads 1.000000 on station-days that are 58 % gap
    (UW.GNW 2010-08-12, coverage_frac 0.418) -- the cache's own gap indicator is blind;
  * the correlation engines therefore defend themselves with a RELATIVE energy guard
    (`ew > 1e-4 * median(ew)`, k_auto.py:70, k_grow.py:167), which silently degenerates to
    `ew > 0` on any channel more than half dead, and masked windows still enter the MAD pool.
Here the gap is known exactly, from the trace spans BEFORE the fill, and is written as `m`.

  python scripts/ab_cache40.py --sub NW [--workers 24] [--force]
-> cache40hz/<NET.STA>/<YYYY>/<DOY>.npz   x (3,N) float32 ZNE 1-8 Hz @40 sps, t0, m (N,) bool,
                                          order, fs_native, rotated, frac_valid
-> results/ab_royer/cache40_<SUB>.csv     one row per station-day
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from scipy.signal import butter, resample_poly, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import ab_common as A                                                       # noqa: E402

# A 4th-order Butterworth applied twice (sosfiltfilt) at a 1 Hz corner rings for a few cycles; the
# resample_poly anti-alias kernel adds its own support.  Both are handled by widening every gap by
# TRANSIENT_S at each edge in the MASK only -- the waveform itself is left exactly as the frozen
# chain would produce it, so the mask is a pure addition and can be ignored to reproduce the old
# behaviour for an ablation column.
TRANSIENT_S = 3.0


def _order_and_mask(st, sta, t_day_fixed=None):
    """(order, X (C,N) float64 zero-filled, valid (N,) bool, fs, t0).

    `valid` is built from the trace spans BEFORE the merge fill, AND-ed over components, so a sample
    is valid only where every component of the station actually recorded."""
    import o_chanorder as chanorder
    for tr in st:
        fs_ = float(tr.stats.sampling_rate)
        fr = round(fs_)
        if fr > 0 and fs_ != fr and abs(fs_ - fr) / fr < 1e-4:
            tr.stats.sampling_rate = float(fr)
    fs = float(st[0].stats.sampling_rate)
    if any(abs(float(tr.stats.sampling_rate) - fs) > 1e-9 for tr in st):
        raise RuntimeError(f"{sta}: mixed sampling rates {sorted({tr.stats.sampling_rate for tr in st})}")
    chans = sorted({tr.stats.channel for tr in st})
    order = chanorder.resolve_order(chans, station=sta, strict=False)
    if order[0] and not order[1] and not order[2]:
        # vertical-only station-day (UW.SQM after 2011-04-22 has EHZ + accelerometers only).  The
        # horizontals are stored as exact zeros: every downstream stage drops a zero-energy
        # channel per component, so the station contributes its one real channel (recorded as
        # `z_only`) instead of being dropped whole.
        order = [order[0], "", ""]
    elif order[0] and not all(order) and all((c == "" or c is None or c[-1] in "NE")
                                            for c in order[1:]):
        # one N/E horizontal missing (e.g. UW.OPC 2009-05: BHZ + BHE only): keep what exists,
        # the missing channel is stored as zeros and dropped per component downstream
        order = [c or "" for c in order]
    elif not all(order):
        raise RuntimeError(f"{sta}: not three components ({chans})")
    # one absolute day grid so every station shares a time base
    t_day = min(float(tr.stats.starttime.timestamp) for tr in st)
    t_day = float(pd.Timestamp(t_day, unit="s").floor("D").timestamp())
    if t_day_fixed is not None:
        # MENDO (2026-09-29): lfe_south raw day files are padded +-1 min (23:59:00 of the previous day
        # .. 00:01:00 of the next), so the grid is the REQUESTED day's midnight, not the earliest sample
        t_day = float(t_day_fixed)
    n = int(round(86400 * fs))
    X = np.zeros((3, n))
    V = np.zeros((3, n), bool)
    for ci, c in enumerate(order):
        if not c:
            continue
        for tr in st.select(channel=c):
            i0 = int(round((float(tr.stats.starttime.timestamp) - t_day) * fs))
            a, b = max(i0, 0), min(i0 + int(tr.stats.npts), n)
            if b <= a:
                continue
            X[ci, a:b] = np.asarray(tr.data[a - i0:b - i0], float)
            V[ci, a:b] = True
    valid = V[[i for i, c in enumerate(order) if c]].all(0)
    return order, X, valid, fs, t_day


def _widen(valid, k):
    """Mark k extra samples invalid on each side of every invalid run (filter transient guard)."""
    if k <= 0 or valid.all():
        return valid
    bad = (~valid).astype(np.int16)
    c = np.concatenate([[0], np.cumsum(bad)])
    n = valid.size
    lo = np.maximum(np.arange(n) - k, 0)
    hi = np.minimum(np.arange(n) + k + 1, n)
    return (c[hi] - c[lo]) == 0


def _mask_to_fs(valid, fs0, fs1, pad_s):
    """Carry a native-rate mask to `fs1` conservatively: an output sample is valid only if the whole
    native interval it draws from, widened by pad_s, is valid."""
    n1 = int(round(valid.size * fs1 / fs0))
    idx = (np.arange(n1) * (fs0 / fs1)).astype(np.int64)
    k = int(np.ceil(pad_s * fs0)) + int(np.ceil(fs0 / fs1))
    w = _widen(valid, k)
    idx = np.clip(idx, 0, valid.size - 1)
    return w[idx]


def cache_one(job):
    sta, date, force = job
    t0w = time.time()
    dst = A.cache_path(sta, date)
    src = A.raw_path(sta, date)
    base = dict(sta=sta, netsta=A.netsta(sta), date=date)
    try:
        if os.path.exists(dst) and not force:
            z = np.load(dst)
            return dict(base, status="cached", n=int(z["x"].shape[1]),
                        frac_valid=float(z["m"].mean()) if "m" in z.files else np.nan,
                        fs_native=float(z["fs_native"]) if "fs_native" in z.files else np.nan,
                        order=",".join(list(z["order"])), rotated=bool(z["rotated"])
                        if "rotated" in z.files else False,
                        secs=round(time.time() - t0w, 1), err="")
        if not os.path.exists(src):
            return dict(base, status="noraw", n=0, frac_valid=0.0, fs_native=np.nan, order="",
                        rotated=False, secs=0.0, err="raw miniSEED absent")
        from obspy import read
        st = read(src)
        st, band_tag = A.pick_band(st)
        order, X, valid, fs0, t_day = _order_and_mask(
            st, sta, pd.Timestamp(date).timestamp() if sta in A.MENDO_SET else None)
        # ---- rotation to ZNE, the same convention as the frozen 20 Hz cache (k_common.load_day)
        rotated = False
        if order[1] and order[1][-1] in "12":
            import s_common as S
            az = A.AZ12_EXTRA.get(sta) or S.AZ12.get(sta)
            if not az:
                try:
                    az = A.az12(sta, date, [order[1], order[2]])
                except Exception:
                    az = None
            if az:
                a1, a2 = az
                r1, r2 = np.radians(a1), np.radians(a2)
                M = np.array([[np.cos(r1), np.cos(r2)], [np.sin(r1), np.sin(r2)]])
                X[1:] = M @ X[1:]
                rotated = True
            else:
                return dict(base, status="noazimuth", n=0, frac_valid=0.0, fs_native=fs0,
                            order=",".join(order), rotated=False, secs=round(time.time() - t0w, 1),
                            err="1/2 channels with no StationXML azimuth: cannot make N/E")
        # ---- band-pass at the NATIVE rate, then rational resample to exactly 40 sps
        X = X - X.mean(1, keepdims=True)
        sos = butter(4, [A.BAND[0] / (fs0 / 2), A.BAND[1] / (fs0 / 2)], btype="band", output="sos")
        X = sosfiltfilt(sos, X, axis=-1)
        if abs(fs0 - A.FS) > 1e-9:
            from fractions import Fraction
            fr = Fraction(fs0 / A.FS).limit_denominator(1000)
            X = resample_poly(X, fr.denominator, fr.numerator, axis=-1)
        m = _mask_to_fs(valid, fs0, A.FS, TRANSIENT_S)
        n = min(X.shape[1], m.size)
        X, m = X[:, :n], m[:n]
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        tmp = dst + f".tmp{os.getpid()}.npz"
        np.savez(tmp, x=X.astype(np.float32), m=m, t0=np.float64(t_day),
                 order=np.array(order), fs=np.float64(A.FS), fs_native=np.float64(fs0),
                 rotated=np.bool_(rotated))
        os.replace(tmp, dst)
        return dict(base, status=("built" if all(order) else "built_z_only"
                                  if not (order[1] or order[2]) else "built_partial"), n=int(n), frac_valid=float(m.mean()), fs_native=fs0,
                    order=",".join(order), rotated=rotated, secs=round(time.time() - t0w, 1),
                    err="")
    except Exception as e:
        return dict(base, status="error", n=0, frac_valid=0.0, fs_native=np.nan, order="",
                    rotated=False, secs=round(time.time() - t0w, 1),
                    err=f"{type(e).__name__}: {e}"[:300])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", default="NW,SW,NVI")
    ap.add_argument("--workers", type=int, default=24)         # pod cgroup cap is 32 CPUs
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--ext", action="store_true", help="the EXTENSION station-days")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    jobs, tag = [], a.sub.replace(",", "")
    for sub in a.sub.split(","):
        if a.ext:
            SD = pd.read_csv(f"{A.EXT_DIR}/station_days_{sub}.csv")
            jobs += [(s, d, a.force) for s, d in zip(SD.sta, SD.date)]
        else:
            jobs += [(s, d, a.force) for s, d in A.station_days(sub)]
    jobs = sorted(set(jobs))
    A.log(f"[cache40] {len(jobs)} station-days, {a.workers} workers")
    rows, t0 = [], time.time()
    import multiprocessing as mp
    with mp.Pool(min(a.workers, 30), maxtasksperchild=20) as pool:
        for i, r in enumerate(pool.imap_unordered(cache_one, jobs, chunksize=1), 1):
            rows.append(r)
            if i % 50 == 0 or i == len(jobs):
                A.log(f"  {i}/{len(jobs)} {time.time()-t0:.0f}s")
    C = pd.DataFrame(rows).sort_values(["netsta", "date"])
    C.to_csv(f"{A.EXT_DIR if a.ext else A.RES}/cache40_{tag}.csv", index=False)
    A.log(C.groupby("status").size().to_string())
    ok = C[C.status.isin(["built", "built_z_only", "built_partial", "cached"])]
    if len(ok):
        A.log(f"  frac_valid: min {ok.frac_valid.min():.4f} p05 {ok.frac_valid.quantile(.05):.4f} "
              f"median {ok.frac_valid.median():.4f}")
        A.log(f"  station-days with any gap (frac_valid < 1): {(ok.frac_valid < 1).sum()}")
        A.log(f"  station-days more than half gap: {(ok.frac_valid < 0.5).sum()}")
        A.log("  native rates: " + str(sorted(set(ok.fs_native.dropna()))))
    for r in C[C.status.isin(["error", "noazimuth", "noraw"])].head(20).itertuples():
        A.log(f"  {r.status.upper()} {r.netsta} {r.date}: {r.err}")


if __name__ == "__main__":
    main()
