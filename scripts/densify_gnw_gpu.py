"""GPU matched-filter densification (NVIDIA L40S) — fast + crash-resilient.

Validated: GPU cc is bit-identical to the CPU fftconvolve path (max|Δcc|=0.0); GPU correlation
is ~180x faster (44 ms/day vs 7.9 s/day). See notes/NISQUALLY_GNW_RESUME.md.

Architecture (built for the 32-CPU + 1x L40S pod, scales to 100s of stations):
  * CPU worker pool: read + resample_poly(->40Hz) + bandpass  (the cheap part, parallel)
  * GPU (main process only): one rfft(day) + batched template multiply + batched irfft,
    GPU sliding-norm, cc>1 fix, threshold + sparse candidate extraction (only detections
    cross PCIe back, never the full cc array)
  * forward-chronological, YEARLY chunks, RESUME by skipping years whose csv exists
    (atomic .tmp+rename so a crash never leaves a partial year)
  * per-day QC (anomalous/glitch days dropped) BEFORE the top-N cap
  * top-N cap per (family, day): keep the N highest-cc detections

Run (GPU lives only in main; workers never import cupy):
  python scripts/densify_gnw_gpu.py            # all years, resumable
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import sys
import time
import warnings
from math import gcd
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, "src")
from scipy.signal import resample_poly  # noqa: E402

# ---------------- CPU preprocessing (pool workers; NO cupy here) ----------------
_G = {}


def _init_worker(root, station, fs, bandpass):
    _G.update(root=root, station=station, fs=fs, bandpass=bandpass)


def _preprocess(day):
    """read + resample_poly->fs + detrend + bandpass. Returns (day, data_f32, start_epoch)."""
    from obspy import read
    root, station, fs, bp = _G["root"], _G["station"], _G["fs"], _G["bandpass"]
    cands = list(root.glob(f"*.{station}/{day.year}/{day.timetuple().tm_yday:03d}.mseed"))
    if not cands:
        return day, None, 0.0
    try:
        st = read(str(cands[0])).select(component="Z")
    except Exception:
        return day, None, 0.0
    if len(st) == 0:
        return day, None, 0.0
    for tr in st:
        fs0 = tr.stats.sampling_rate
        if abs(fs0 - fs) > 1e-6:
            r = int(round(fs0))
            if r <= 0:
                return day, None, 0.0
            g = gcd(int(fs), r)
            try:
                tr.data = resample_poly(tr.data.astype(np.float64), int(fs)//g, r//g).astype(np.float32)
                tr.stats.sampling_rate = float(fs)
            except Exception:
                try:
                    tr.resample(fs)
                except Exception:
                    return day, None, 0.0
    try:
        st.merge(fill_value=0)
    except Exception:
        return day, None, 0.0
    if len(st) == 0:
        return day, None, 0.0
    st.detrend("demean")
    st.filter("bandpass", freqmin=bp[0], freqmax=bp[1], corners=4, zerophase=True)
    tr = st[0]
    return day, tr.data.astype(np.float32), float(tr.stats.starttime.timestamp)


# ---------------- CPU peak-pick over sparse GPU candidates ----------------
def _dedup_min_gap(ti, si, vals, n_gap):
    """ti row-major sorted; within each template si ascending. Greedy min-gap peak-pick.
    Returns {template_index: (sample_idx_array, cc_array)}."""
    out = {}
    if ti.size == 0:
        return out
    # boundaries between templates (ti is non-decreasing)
    btmpl = np.unique(ti)
    starts = np.searchsorted(ti, btmpl, side="left")
    ends = np.searchsorted(ti, btmpl, side="right")
    for t, a, b in zip(btmpl.tolist(), starts.tolist(), ends.tolist()):
        s = si[a:b]; c = vals[a:b]   # s ascending
        picks = []
        i = 0; L = s.size
        while i < L:
            j = s[i]; w_end = j + n_gap
            best = i; k = i + 1
            while k < L and s[k] < w_end:
                if c[k] > c[best]:
                    best = k
                k += 1
            picks.append(best)
            nxt = s[best] + n_gap
            while i < L and s[i] < nxt:
                i += 1
        pi = np.asarray(picks)
        out[int(t)] = (s[pi], c[pi])
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--templates-npz", default="data/gnw_pnsn_families.npz")
    p.add_argument("--summary-csv", default="data/gnw_pnsn_families.summary.csv")
    p.add_argument("--min-snr", type=float, default=15.0)
    p.add_argument("--wfdir", default="data/waveforms")
    p.add_argument("--network", default="UW")
    p.add_argument("--station", default="GNW")
    p.add_argument("--fs", type=float, default=40.0)
    p.add_argument("--fmin", type=float, default=2.0)
    p.add_argument("--fmax", type=float, default=8.0)
    p.add_argument("--threshold", type=float, default=0.8)
    p.add_argument("--min-gap-s", type=float, default=6.0)
    p.add_argument("--workers", type=int, default=12)   # CPU preprocessors (feed the GPU)
    p.add_argument("--top-n", type=int, default=100)
    p.add_argument("--qc-median-count", type=int, default=2000)
    p.add_argument("--qc-median-cc", type=float, default=0.96)
    p.add_argument("--start-year", type=int, default=0)
    p.add_argument("--end-year", type=int, default=9999)
    p.add_argument("--out-dir", default="data")
    p.add_argument("--out-prefix", default="mf_gnw_")
    args = p.parse_args()

    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor, as_completed
    import cupy as cp  # GPU only in main

    print(f"[gpu] {cp.cuda.runtime.getDeviceProperties(0)['name'].decode()}", flush=True)
    d = np.load(args.templates_npz, allow_pickle=True)
    s = pd.read_csv(args.summary_csv).set_index("family_id")
    keys = [k for k in s[s["snr"] >= args.min_snr].index if k in d.files]
    m = int(np.asarray(d[keys[0]]).size)
    Tn = np.empty((len(keys), m), np.float32)        # normalized, reversed
    for i, k in enumerate(keys):
        t0 = np.asarray(d[k], np.float64); t0 = t0 - t0.mean()
        nr = np.linalg.norm(t0)
        Tn[i] = (t0 / nr if nr else t0)[::-1]
    Tn_g = cp.asarray(Tn)
    print(f"[load] {len(keys)} templates, m={m}", flush=True)

    # FIXED day length -> constant N -> ONE template-FFT matrix + all GPU buffers reused by
    # the cupy pool -> bounded GPU memory (~4 GB). (Variable N from gappy days cached a 0.73 GB
    # Tk per distinct length and OOM'd the 46 GB GPU.)
    import scipy.fft as _sf
    N_FIX = 3_456_016                 # >= any full day (40Hz day=3,456,000; 100->40=3,456,001)
    valid_fix = N_FIX - m + 1
    N = _sf.next_fast_len(N_FIX + m - 1)
    Tk = cp.fft.rfft(Tn_g, n=N, axis=1).astype(cp.complex64)
    print(f"[gpu] N_FIX={N_FIX:,} N={N:,}, Tk resident {Tk.nbytes/1e9:.2f} GB", flush=True)

    thr = np.float32(args.threshold)
    n_gap = max(1, int(round(args.min_gap_s * args.fs)))
    qc_count = args.qc_median_count or None
    qc_cc = args.qc_median_cc or None
    _empty = (np.empty(0, np.int64), np.empty(0, np.int64), np.empty(0, np.float32))

    def gpu_candidates(data):
        n_real = int(data.size)
        if n_real < m:
            return _empty
        # pad (or truncate) to N_FIX so every day uses the same N
        if n_real >= N_FIX:
            dg = cp.asarray(data[:N_FIX]); n_real = N_FIX
        else:
            dg = cp.zeros(N_FIX, dtype=cp.float32); dg[:n_real] = cp.asarray(data)
        D = cp.fft.rfft(dg, n=N).astype(cp.complex64)
        conv = cp.fft.irfft(D[None, :] * Tk, n=N, axis=1)[:, m - 1:m - 1 + valid_fix]
        # GPU sliding norm (float64 cumsum for accuracy, cast result to f32 to halve mem)
        x = dg.astype(cp.float64)
        c1 = cp.concatenate((cp.zeros(1), cp.cumsum(x))); wm = (c1[m:] - c1[:-m]) / m
        c2 = cp.concatenate((cp.zeros(1), cp.cumsum(x * x)))
        ws = cp.sqrt(cp.maximum((c2[m:] - c2[:-m]) - m * wm * wm, 1e-12)).astype(cp.float32)
        cc = conv / ws[None, :]
        cc = cp.where(cp.isfinite(cc) & (cp.abs(cc) <= 1.0), cc, cp.float32(0.0))
        # only windows fully inside REAL samples (padding region gives cc~0 anyway)
        real_valid = max(0, n_real - m + 1)
        mask = cc[:, :real_valid] >= thr
        ti, si = cp.nonzero(mask)
        vals = cc[ti, si]
        out = (cp.asnumpy(ti), cp.asnumpy(si), cp.asnumpy(vals).astype(np.float32))
        del dg, D, conv, x, c1, c2, ws, cc, mask, ti, si, vals
        return out

    # enumerate by year
    sdir = Path(args.wfdir) / f"{args.network}.{args.station}"
    by_year = {}
    for ypath in sorted(sdir.glob("[0-9]*")):
        try:
            year = int(ypath.name)
        except ValueError:
            continue
        if not (args.start_year <= year <= args.end_year):
            continue
        days = []
        for mfile in sorted(ypath.glob("*.mseed")):
            try:
                jday = int(mfile.stem)
            except ValueError:
                continue
            days.append(datetime(year, 1, 1) + timedelta(days=jday - 1))
        if days:
            by_year[year] = sorted(days)
    years = sorted(by_year)
    print(f"[plan] {len(years)} years {years[0]}-{years[-1]}, "
          f"{sum(len(v) for v in by_year.values())} day-files; "
          f"cpu_workers={args.workers}, top_n={args.top_n}, cc>={args.threshold}", flush=True)

    out_dir = Path(args.out_dir); out_dir.mkdir(exist_ok=True)
    ctx = mp.get_context("spawn")
    station = args.station; fs = args.fs

    pool = ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx,
                               initializer=_init_worker,
                               initargs=(Path(args.wfdir), station, fs, (args.fmin, args.fmax)))
    try:
        for year in years:
            out_csv = out_dir / f"{args.out_prefix}{year}.csv"
            if out_csv.exists():
                print(f"[skip] {year}: exists ({out_csv.stat().st_size/1e6:.1f} MB)", flush=True)
                continue
            days = by_year[year]
            t0 = time.time()
            rows = []; dropped = 0
            futs = [pool.submit(_preprocess, dd) for dd in days]
            for f in as_completed(futs):
                day, data, start_ep = f.result()
                if data is None:
                    continue
                ti, si, vals = gpu_candidates(data)
                res = _dedup_min_gap(ti, si, vals, n_gap)   # {tmpl_idx: (samp, cc)}
                if not res:
                    continue
                # QC on full-day post-dedup detections (BEFORE cap)
                counts = [v[0].size for v in res.values()]
                med_count = int(np.median(counts)) if counts else 0
                all_cc = np.concatenate([v[1] for v in res.values()]) if res else np.array([0.0])
                med_cc = float(np.median(all_cc))
                if ((qc_count is not None and med_count > qc_count) or
                        (qc_cc is not None and med_cc > qc_cc)):
                    dropped += 1
                    continue
                # top-N cap per (family, day), keep highest cc
                for tix, (samp, ccv) in res.items():
                    if args.top_n and ccv.size > args.top_n:
                        sel = np.argpartition(ccv, -args.top_n)[-args.top_n:]
                    else:
                        sel = np.arange(ccv.size)
                    tkey = keys[tix]
                    tt = start_ep + samp[sel] / fs
                    for e, c in zip(tt.tolist(), ccv[sel].tolist()):
                        rows.append((tkey, e, c))
            # atomic write
            df = pd.DataFrame(rows, columns=["template", "time", "cc"])
            df["time"] = pd.to_datetime(df["time"], unit="s")
            df["station"] = station
            tmp = out_csv.with_suffix(".csv.tmp")
            df.to_csv(tmp, index=False)
            os.replace(tmp, out_csv)
            dt = time.time() - t0
            print(f"[done] {year}: {len(days)} days in {dt/60:.2f} min "
                  f"({len(days)/max(dt,1e-9):.1f}/s), {dropped} QC-dropped, "
                  f"{len(rows):,} rows -> {out_csv.name}", flush=True)
    finally:
        pool.shutdown(wait=True)
    print("[ALL YEARS DONE]", flush=True)


if __name__ == "__main__":
    main()
