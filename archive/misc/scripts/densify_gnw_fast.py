"""FAST, crash-resilient matched-filter densification of the GNW seed templates.

Designed for the 32-CPU pod (see notes/NISQUALLY_GNW_RESUME.md):
  * --workers 30  (pod cgroup cap is 32 CPUs; 48 thrashed -> ~13h. 30 -> ~1h)
  * thread env pinned to 1 (set at top, before numpy) so 30 workers don't oversubscribe
  * resample_poly (25x faster than obspy FFT resample); keeps fs=40 (NLLB-comparable)
  * day FFT computed ONCE per day, reused across all templates (was recomputed 53x/day)
  * forward-chronological, YEARLY chunks: writes data/mf_gnw_<YEAR>.csv per year and
    RESUMES by skipping years whose file already exists (atomic write via .tmp + rename)
  * per-day QC (anomalous/glitch days dropped) BEFORE the cap
  * cc>1 normalization fix
  * top-N cap per (family, day): keep the N highest-cc detections (stacks only need ~20)

Run:
  export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
  python scripts/densify_gnw_fast.py            # all years, resumable
"""
from __future__ import annotations

# ---- pin threads BEFORE numpy import (spawned workers re-import this module) ----
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
import scipy.fft as sf  # noqa: E402

# ----------------------------- per-day scan (fast) -----------------------------

def _load_day_fast(root: Path, station: str, day: datetime, bandpass, fs):
    from obspy import read
    cands = list(root.glob(f"*.{station}/{day.year}/{day.timetuple().tm_yday:03d}.mseed"))
    if not cands:
        return None
    try:
        st = read(str(cands[0]))
    except Exception:
        return None
    st = st.select(component="Z")
    if len(st) == 0:
        return None
    for tr in st:
        fs0 = tr.stats.sampling_rate
        if abs(fs0 - fs) > 1e-6:
            r = int(round(fs0))
            if r <= 0:
                return None
            g = gcd(int(fs), r)
            up, down = int(fs) // g, r // g
            try:
                tr.data = resample_poly(tr.data.astype(np.float64), up, down).astype(np.float32)
                tr.stats.sampling_rate = float(fs)
            except Exception:
                # fall back to obspy's resampler for odd/jittery rates
                try:
                    tr.resample(fs)
                except Exception:
                    return None
    try:
        st.merge(fill_value=0)
    except Exception:
        return None
    if len(st) == 0:
        return None
    st.detrend("demean")
    st.filter("bandpass", freqmin=bandpass[0], freqmax=bandpass[1], corners=4, zerophase=True)
    tr = st[0]
    return tr.data.astype(np.float32), tr.stats.starttime


def _sliding_norm(x, m):
    c1 = np.concatenate(([0.0], np.cumsum(x.astype(np.float64))))
    wm = (c1[m:] - c1[:-m]) / m
    c2 = np.concatenate(([0.0], np.cumsum(x.astype(np.float64) ** 2)))
    return np.sqrt(np.maximum((c2[m:] - c2[:-m]) - m * wm ** 2, 1e-12))


_TK_CACHE = {}  # per-process: (N, m, ntmpl) -> complex template-FFT matrix


def _scan_day(root, station, day, templates, fs, bandpass, threshold, min_gap_s):
    """Return {template_key: (times, ccs)} for one day, or None if no data."""
    loaded = _load_day_fast(root, station, day, bandpass, fs)
    if loaded is None:
        return None
    data, start = loaded
    n = data.size
    keys = list(templates.keys())
    m = templates[keys[0]].size
    if n == 0 or m > n:
        return None
    win_std = _sliding_norm(data, m)
    N = sf.next_fast_len(n + m - 1)
    ck = (N, m, len(keys))
    Tk = _TK_CACHE.get(ck)
    if Tk is None:
        # complex64 keeps the per-worker cache small (~0.7 GB vs 1.5 GB at c128);
        # cc precision to ~1e-3 is plenty for a >=0.8 threshold.
        Tk = np.empty((len(keys), N // 2 + 1), dtype=np.complex64)
        for i, k in enumerate(keys):
            t0 = templates[k].astype(np.float64)
            t0 = t0 - t0.mean()
            nr = np.linalg.norm(t0)
            if nr:
                t0 = t0 / nr
            Tk[i] = sf.rfft(t0[::-1], N).astype(np.complex64)
        _TK_CACHE.clear()  # only ever need the current-length matrix
        _TK_CACHE[ck] = Tk
    D = sf.rfft(data, N).astype(np.complex64)  # day FFT ONCE
    valid = n - m + 1
    n_gap = max(1, int(round(min_gap_s * fs)))
    out = {}
    for i, k in enumerate(keys):
        num = sf.irfft(D * Tk[i], N)[m - 1: m - 1 + valid]
        cc = (num / win_std).astype(np.float32)
        cc[~np.isfinite(cc)] = 0.0
        cc[np.abs(cc) > 1.0] = 0.0  # cc>1 fix (zero-fill-gap artifacts)
        picks = []
        j = 0
        sz = cc.size
        while j < sz:
            if cc[j] >= threshold:
                je = min(sz, j + n_gap)
                p = int(np.argmax(cc[j:je])) + j
                picks.append(p)
                j = p + n_gap
            else:
                j += 1
        if picks:
            out[k] = ([(start + (p / fs)).datetime for p in picks],
                      [float(cc[p]) for p in picks])
    return out


# module-level globals so spawned workers share the (pickled-once) templates
_G = {}


def _init_worker(root, station, templates, fs, bandpass, threshold, min_gap_s):
    _G.update(root=root, station=station, templates=templates, fs=fs,
              bandpass=bandpass, threshold=threshold, min_gap_s=min_gap_s)


def _worker(day):
    res = _scan_day(_G["root"], _G["station"], day, _G["templates"], _G["fs"],
                    _G["bandpass"], _G["threshold"], _G["min_gap_s"])
    return day, res


# ----------------------------- QC + cap -----------------------------

def _is_bad_day(res, qc_count, qc_cc):
    """Drop glitch days: all families saturate near-equally at near-perfect cc."""
    if qc_count is None and qc_cc is None:
        return False
    counts = [len(t) for t, _ in res.values()]
    if not counts:
        return False
    med_count = int(np.median(counts))
    all_cc = [c for _, ccs in res.values() for c in ccs]
    med_cc = float(np.median(all_cc)) if all_cc else 0.0
    return ((qc_count is not None and med_count > qc_count) or
            (qc_cc is not None and med_cc > qc_cc))


def _cap_rows(res, station, top_n):
    """Per (family, day): keep top_n highest-cc detections. Returns list of row dicts."""
    rows = []
    for k, (times, ccs) in res.items():
        if top_n and len(ccs) > top_n:
            idx = np.argpartition(np.asarray(ccs), -top_n)[-top_n:]
            sel = sorted(idx.tolist())
        else:
            sel = range(len(ccs))
        for j in sel:
            rows.append({"template": k, "time": times[j], "cc": ccs[j], "station": station})
    return rows


# ----------------------------- driver -----------------------------

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
    p.add_argument("--workers", type=int, default=30)      # pod cap is 32 CPUs
    p.add_argument("--top-n", type=int, default=100)        # cap per family-day
    p.add_argument("--qc-median-count", type=int, default=2000)
    p.add_argument("--qc-median-cc", type=float, default=0.96)
    p.add_argument("--start-year", type=int, default=0)
    p.add_argument("--end-year", type=int, default=9999)
    p.add_argument("--out-dir", default="data")
    p.add_argument("--out-prefix", default="mf_gnw_")
    args = p.parse_args()

    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor, as_completed

    print(f"[load] templates from {args.templates_npz} (SNR>={args.min_snr})", flush=True)
    d = np.load(args.templates_npz, allow_pickle=True)
    s = pd.read_csv(args.summary_csv).set_index("family_id")
    hi = s[s["snr"] >= args.min_snr].index
    T = {k: np.asarray(d[k], dtype=np.float32) for k in hi if k in d.files}
    print(f"[load] {len(T)} seed templates, len={next(iter(T.values())).size}", flush=True)

    wfdir = Path(args.wfdir)
    sdir = wfdir / f"{args.network}.{args.station}"
    # enumerate day-files grouped by year
    by_year: dict[int, list[datetime]] = {}
    for ypath in sorted(sdir.glob("[0-9]*")):
        try:
            year = int(ypath.name)
        except ValueError:
            continue
        if year < args.start_year or year > args.end_year:
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
    total_days = sum(len(v) for v in by_year.values())
    print(f"[plan] {len(years)} years {years[0]}-{years[-1]}, {total_days} day-files; "
          f"workers={args.workers}, top_n={args.top_n}, cc>={args.threshold}", flush=True)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)
    ctx = mp.get_context("spawn")
    qc_count = args.qc_median_count or None
    qc_cc = args.qc_median_cc or None

    for year in years:
        out_csv = out_dir / f"{args.out_prefix}{year}.csv"
        if out_csv.exists():
            print(f"[skip] {year}: {out_csv.name} exists ({out_csv.stat().st_size/1e6:.1f} MB)", flush=True)
            continue
        days = by_year[year]
        t0 = time.time()
        rows = []
        dropped = 0
        done = 0
        with ProcessPoolExecutor(
            max_workers=args.workers, mp_context=ctx,
            initializer=_init_worker,
            initargs=(wfdir, args.station, T, args.fs, (args.fmin, args.fmax),
                      args.threshold, args.min_gap_s),
        ) as ex:
            futs = [ex.submit(_worker, dd) for dd in days]
            for f in as_completed(futs):
                day, res = f.result()
                done += 1
                if res:
                    if _is_bad_day(res, qc_count, qc_cc):
                        dropped += 1
                    else:
                        rows.extend(_cap_rows(res, args.station, args.top_n))
        # atomic write: .tmp then rename, so a crash never leaves a partial year file
        tmp = out_csv.with_suffix(".csv.tmp")
        df = pd.DataFrame(rows, columns=["template", "time", "cc", "station"])
        df.to_csv(tmp, index=False)
        os.replace(tmp, out_csv)
        dt = time.time() - t0
        print(f"[done] {year}: {len(days)} days in {dt/60:.1f} min "
              f"({len(days)/max(dt,1e-9)*60:.0f}/min), {dropped} QC-dropped, "
              f"{len(rows):,} rows -> {out_csv.name}", flush=True)

    print("[ALL YEARS DONE]", flush=True)


if __name__ == "__main__":
    main()
