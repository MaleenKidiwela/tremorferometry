#!/usr/bin/env python
"""1-10 Hz / 40 sps day cache for the Thomas & Bostock (2015) central-Oregon reproduction.

Thomas & Bostock 2015 sec.2: "Data were demeaned, detrended, band-pass filtered between 1 and 10 Hz,
resampled to 40 sps, and stored in 24 h data files."  No instrument correction is mentioned (and the
Savard & Bostock detector normalises each hour per station, matlab_original/main_scc.m:285-287), so
none is applied here: the cache is in counts.

Adapted (copied, not imported) from lfe_array_discovery/scripts/ab_cache40.py: exact validity mask
built from the trace spans BEFORE the zero fill, widened by TRANSIENT_S at every gap edge.

Deviations the data force (logged in NOTES.md):
  * IU.COR is BH at 20 sps: its Nyquist is 10 Hz, so its band is 1-9 Hz (0.9 x Nyquist) and it is
    then up-sampled to 40 sps.  Every other station has native fs >= 40.
  * PB boreholes (EH1/EH2) are rotated to N/E with the StationXML azimuths (meta/pb_azimuths.json);
    Z with dip +90 (B031) is sign-flipped to up-positive.  Orientation accuracy is NOT checked here.
  * EHZ-only UW short-periods (BRO, HSO, MPO, RNO, SSO) are stored with N/E = exact zeros.

  python scripts/or_cache40.py [--workers 12] [--force] [--sta XA.ALP0,...] [--dates 2009-09-03,...]
-> cache40/<NET.STA>/<YYYY-MM-DD>.npz : x (3,N) float32 Z,N,E ; m (N,) bool valid ; t0 ; fs ; order
-> results/cache40_manifest.csv
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import time
import warnings
from fractions import Fraction

import numpy as np
import pandas as pd
from scipy.signal import butter, resample_poly, sosfiltfilt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = f"{ROOT}/raw"
CACHE = f"{ROOT}/cache40"
FS = 40.0
BAND = (1.0, 10.0)
TRANSIENT_S = 3.0
BAND_PREF = ["HH", "BH", "EH", "SH"]
AZ = {}
for r in json.load(open(f"{ROOT}/meta/pb_azimuths.json")):
    AZ[(r["sta"], r["ch"])] = (float(r["az"]), float(r["dip"]))


def _widen(valid, k):
    if k <= 0 or valid.all():
        return valid
    bad = (~valid).astype(np.int32)
    c = np.concatenate([[0], np.cumsum(bad)])
    n = valid.size
    lo = np.maximum(np.arange(n) - k, 0)
    hi = np.minimum(np.arange(n) + k + 1, n)
    return (c[hi] - c[lo]) == 0


def _mask_to_fs(valid, fs0, fs1, pad_s):
    n1 = int(round(valid.size * fs1 / fs0))
    idx = np.clip((np.arange(n1) * (fs0 / fs1)).astype(np.int64), 0, valid.size - 1)
    k = int(np.ceil(pad_s * fs0)) + int(np.ceil(fs0 / fs1))
    return _widen(valid, k)[idx]


def _pick_band(st):
    groups = {}
    for tr in st:
        groups.setdefault((tr.stats.location, tr.stats.channel[:2]), set()).add(tr.stats.channel[-1])
    ok = [k for k, v in groups.items() if "Z" in v]
    ok.sort(key=lambda k: (-len(groups[k]), BAND_PREF.index(k[1]) if k[1] in BAND_PREF else 99, k[0] != ""))
    loc, band = ok[0]
    return st.select(location=loc, channel=band + "?")


def cache_one(job):
    netsta, date, force = job
    t0w = time.time()
    src = f"{RAW}/{netsta}/{date}.mseed"
    dst = f"{CACHE}/{netsta}/{date}.npz"
    base = dict(netsta=netsta, date=date)
    try:
        if os.path.exists(dst) and not force:
            z = np.load(dst)
            return dict(base, status="cached", frac_valid=float(z["m"].mean()), order=",".join(z["order"]),
                        fs_native=float(z["fs_native"]), band_hi=float(z["band"][1]), secs=0.0, err="")
        if not os.path.exists(src):
            return dict(base, status="noraw", frac_valid=0.0, order="", fs_native=np.nan, band_hi=np.nan,
                        secs=0.0, err="")
        from obspy import read
        st = _pick_band(read(src))
        for tr in st:
            f_, fr = float(tr.stats.sampling_rate), round(float(tr.stats.sampling_rate))
            if fr > 0 and f_ != fr and abs(f_ - fr) / fr < 1e-4:
                tr.stats.sampling_rate = float(fr)
        # XA.BWO0 2009-09-09/10 interleave 40 Hz packets with the 100 Hz HH stream under one
        # channel code: keep the rate holding the most samples (the minority packets are dropped
        # and show up as gaps in the mask)
        from collections import Counter
        cnt = Counter()
        for tr in st:
            cnt[float(tr.stats.sampling_rate)] += int(tr.stats.npts)
        fs0 = max(cnt, key=cnt.get)
        if len(cnt) > 1:
            for tr in [tr for tr in st if float(tr.stats.sampling_rate) != fs0]:
                st.remove(tr)
        if any(abs(tr.stats.sampling_rate - fs0) > 1e-9 for tr in st):
            raise RuntimeError("mixed sampling rates")
        chans = sorted({tr.stats.channel for tr in st})
        z = [c for c in chans if c[-1] == "Z"]
        n_ = [c for c in chans if c[-1] in "N1"]
        e_ = [c for c in chans if c[-1] in "E2"]
        order = [z[0] if z else "", n_[0] if n_ else "", e_[0] if e_ else ""]
        t_day = pd.Timestamp(date).timestamp()
        n = int(round(86400 * fs0))
        X = np.zeros((3, n))
        V = np.zeros((3, n), bool)
        for ci, c in enumerate(order):
            if not c:
                continue
            for tr in st.select(channel=c):
                i0 = int(round((float(tr.stats.starttime.timestamp) - t_day) * fs0))
                a, b = max(i0, 0), min(i0 + int(tr.stats.npts), n)
                if b > a:
                    X[ci, a:b] = np.asarray(tr.data[a - i0:b - i0], float)
                    V[ci, a:b] = True
        live = [i for i, c in enumerate(order) if c]
        valid = V[live].all(0)
        sta = netsta.split(".")[1]
        rotated = False
        if order[1].endswith("1"):
            (a1, _), (a2, _) = AZ[(sta, order[1])], AZ[(sta, order[2])]
            r1, r2 = np.radians(a1), np.radians(a2)
            M = np.array([[np.cos(r1), np.cos(r2)], [np.sin(r1), np.sin(r2)]])
            X[1:] = M @ X[1:]
            rotated = True
        if (sta, order[0]) in AZ and AZ[(sta, order[0])][1] > 0:     # dip +90 = down-positive
            X[0] *= -1
        # demean + detrend over valid samples only, per channel (zeros stay zeros outside data)
        tt = np.arange(n, dtype=float)
        for ci in live:
            vv = V[ci]
            if vv.sum() > 10:
                p = np.polyfit(tt[vv], X[ci, vv], 1)
                X[ci, vv] -= np.polyval(p, tt[vv])
            X[ci, ~vv] = 0.0
        hi = min(BAND[1], 0.45 * fs0)
        sos = butter(4, [BAND[0] / (fs0 / 2), hi / (fs0 / 2)], btype="band", output="sos")
        X = sosfiltfilt(sos, X, axis=-1)
        if abs(fs0 - FS) > 1e-9:
            fr = Fraction(fs0 / FS).limit_denominator(1000)
            X = resample_poly(X, fr.denominator, fr.numerator, axis=-1)
        m = _mask_to_fs(valid, fs0, FS, TRANSIENT_S)
        nn = int(86400 * FS)
        Xo = np.zeros((3, nn), np.float32)
        k = min(nn, X.shape[1])
        Xo[:, :k] = X[:, :k]
        mo = np.zeros(nn, bool)
        mo[:min(nn, m.size)] = m[:nn]
        Xo[:, ~mo] = 0.0
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        tmp = dst + f".tmp{os.getpid()}.npz"
        np.savez(tmp, x=Xo, m=mo, t0=np.float64(t_day), order=np.array(order), fs=np.float64(FS),
                 fs_native=np.float64(fs0), band=np.array([BAND[0], hi]), rotated=np.bool_(rotated))
        os.replace(tmp, dst)
        return dict(base, status="built", frac_valid=float(mo.mean()), order=",".join(order), fs_native=fs0,
                    band_hi=hi, secs=round(time.time() - t0w, 1), err="")
    except Exception as e:
        return dict(base, status="error", frac_valid=0.0, order="", fs_native=np.nan, band_hi=np.nan,
                    secs=round(time.time() - t0w, 1), err=f"{type(e).__name__}: {e}"[:300])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--sta", default="")
    ap.add_argument("--dates", default="")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    stas = a.sta.split(",") if a.sta else sorted(os.path.basename(d) for d in glob.glob(f"{RAW}/*"))
    dates = a.dates.split(",") if a.dates else [str(d.date()) for d in pd.date_range("2009-08-03", "2009-09-17")]
    jobs = [(s, d, a.force) for s in stas for d in dates]
    print(f"[cache40] {len(jobs)} station-days, {a.workers} workers", flush=True)
    import multiprocessing as mp
    rows = []
    with mp.Pool(min(a.workers, 12), maxtasksperchild=10) as pool:
        for i, r in enumerate(pool.imap_unordered(cache_one, jobs), 1):
            rows.append(r)
            if i % 50 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)}", flush=True)
    C = pd.DataFrame(rows).sort_values(["netsta", "date"])
    out = f"{ROOT}/results/cache40_manifest.csv"
    if os.path.exists(out) and (a.sta or a.dates):
        old = pd.read_csv(out)
        C = pd.concat([old[~old.set_index(["netsta", "date"]).index.isin(C.set_index(["netsta", "date"]).index)], C])
        C = C.sort_values(["netsta", "date"])
    C.to_csv(out, index=False)
    print(C.groupby("status").size().to_string(), flush=True)
    for r in C[C.status == "error"].head(20).itertuples():
        print("  ERROR", r.netsta, r.date, r.err, flush=True)


if __name__ == "__main__":
    main()
