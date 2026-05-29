"""Faster batched matched-filter: load each day once, scan many templates.

The original `matched_filter.scan_day` reloads + filters the day file for every
(template, day) pair. With T templates, that's T-x redundant I/O per day.

This module loads the day, computes its FFT once, then runs each template via
cheap freq-domain multiply + irfft + sliding normalization. For T=35 templates
on N=2700 days, expect ~10-20x speedup over the per-task approach.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import fftconvolve


def _load_day_filt(
    root: Path, station: str, day: datetime, bandpass: tuple[float, float], fs: float,
):
    from obspy import read

    candidates = list(root.glob(f"*.{station}/{day.year}/{day.timetuple().tm_yday:03d}.mseed"))
    if not candidates:
        return None
    try:
        st = read(str(candidates[0]))
    except Exception:
        return None
    st = st.select(component="Z")
    if len(st) == 0:
        return None
    # Normalize sampling rates BEFORE merging (some files have slight FS jitter
    # like 99.999... vs 100.0 across segments which breaks merge).
    for tr in st:
        if abs(tr.stats.sampling_rate - fs) > 1e-6:
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


def _sliding_norm(continuous: np.ndarray, m: int) -> np.ndarray:
    """Sliding-window std of length-m windows on `continuous`.

    Returns array of length n-m+1.
    """
    n = continuous.size
    if n < m:
        return np.zeros(0, dtype=np.float64)
    c1 = np.concatenate(([0.0], np.cumsum(continuous.astype(np.float64))))
    win_sum = c1[m:] - c1[:-m]
    win_mean = win_sum / m
    c2 = np.concatenate(([0.0], np.cumsum(continuous.astype(np.float64) ** 2)))
    win_sumsq = c2[m:] - c2[:-m]
    win_var = win_sumsq - m * win_mean ** 2
    return np.sqrt(np.maximum(win_var, 1e-12))


def scan_day_multi(
    waveform_root: Path,
    station: str,
    day: datetime,
    templates: dict[str, np.ndarray],
    fs: float = 40.0,
    bandpass: tuple[float, float] = (2.0, 8.0),
    threshold: float = 0.7,
    min_gap_s: float = 6.0,
) -> dict[str, tuple[list[datetime], list[float]]] | None:
    """Apply many templates to one station-day in a single I/O + FFT pass.

    Returns {template_key: (times, ccs)} for each template that produced
    detections, or None if no data for that day.
    """
    loaded = _load_day_filt(waveform_root, station, day, bandpass, fs)
    if loaded is None:
        return None
    data, start = loaded
    n = data.size
    if n == 0:
        return None

    # Sliding norm is shared across templates of the same length.
    # We assume all templates have the same length here. If not, group by length.
    template_lens = {tkey: t.size for tkey, t in templates.items()}
    norms_by_m = {}
    out: dict[str, tuple[list[datetime], list[float]]] = {}
    n_gap = max(1, int(round(min_gap_s * fs)))

    for tkey, template in templates.items():
        m = template.size
        if m > n:
            continue
        if m not in norms_by_m:
            norms_by_m[m] = _sliding_norm(data, m)
        win_std = norms_by_m[m]
        # template -> zero-mean, unit norm
        t0 = template - template.mean()
        nrm_t = float(np.linalg.norm(t0))
        if nrm_t == 0:
            continue
        t0 = t0 / nrm_t
        # numerator
        num = fftconvolve(data, t0[::-1], mode="valid")
        cc = (num / win_std).astype(np.float32)
        # Reject degenerate windows. A normalized (Pearson) CC satisfies
        # |cc| <= 1 by Cauchy-Schwarz; any |cc| > 1 is a numerical artifact
        # from a near-flat window (e.g. a zero-filled data gap), where the
        # true win_std ~ 0 was clamped to the 1e-12 floor while the numerator
        # was not exactly zero. These flooded gappy stations (GNW 1995-2026:
        # 99.9% of cc>=0.8 "detections" were cc>1 artifacts). Set them to 0 so
        # they never pass the threshold.
        cc[~np.isfinite(cc)] = 0.0
        cc[np.abs(cc) > 1.0] = 0.0
        # peak-pick
        picks = []
        i = 0
        while i < cc.size:
            if cc[i] >= threshold:
                j_end = min(cc.size, i + n_gap)
                j = int(np.argmax(cc[i:j_end])) + i
                picks.append(j)
                i = j + n_gap
            else:
                i += 1
        times = [(start + (k / fs)).datetime for k in picks]
        ccs = [float(cc[k]) for k in picks]
        if times:
            out[tkey] = (times, ccs)

    return out


def _worker_task(args):
    """Module-level worker for ProcessPoolExecutor (must be pickleable)."""
    (waveform_root, station, day, templates, fs, bandpass, threshold, min_gap_s) = args
    return day, scan_day_multi(
        waveform_root, station, day, templates,
        fs=fs, bandpass=bandpass, threshold=threshold, min_gap_s=min_gap_s,
    )


def scan_many_days_multi(
    waveform_root: Path,
    station: str,
    days: list[datetime],
    templates: dict[str, np.ndarray],
    fs: float = 40.0,
    bandpass: tuple[float, float] = (2.0, 8.0),
    threshold: float = 0.7,
    min_gap_s: float = 6.0,
    n_workers: int = 24,
    progress_every: int = 50,
    out_path: "Path | str | None" = None,
    flush_every: int = 50,
    qc_median_count_max: "int | None" = None,
    qc_median_cc_max: "float | None" = None,
):
    """Parallel: load each day once, run all templates.

    Two modes:
      * out_path is None (default): accumulate all detections in memory and
        return a DataFrame. Fine for small scans; can OOM on large ones.
      * out_path given: STREAM detections to that CSV as each day's future
        completes, flushing every `flush_every` days and dropping them from
        RAM. Peak memory stays flat (~one batch of days) regardless of total
        detection count. Returns a per-template detection-count dict (small),
        not the full DataFrame.

    Submission is BOUNDED: at most `max_inflight` (= 2*n_workers) tasks are
    submitted at once, and a new one is submitted only as each completes. This
    caps how many completed-but-undrained worker results can buffer in RAM —
    without it, submitting all days up front lets fast workers outrun the
    single-threaded consumer and pile up results until OOM (seen on GNW
    1995-2026, where noisy early years fire ~200K detections/day).
    """
    import multiprocessing as mp
    from collections import Counter
    from concurrent.futures import ProcessPoolExecutor, as_completed

    # Use "spawn" (fresh interpreter per worker) rather than the Linux default
    # "fork". numpy/scipy/obspy start internal thread pools on import; forking
    # while one of those threads holds a lock gives the child a locked mutex
    # with no thread to release it -> the worker deadlocks on its first FFT
    # (seen as workers stuck in futex_wait, 0% CPU). spawn avoids this.
    mp_ctx = mp.get_context("spawn")

    rows = []
    done = 0
    dropped = 0
    total = len(days)
    tmpl_counts: "Counter[str]" = Counter()
    header_written = False

    def _is_bad_day(res):
        """Per-day QC: artifact days saturate ALL templates near-equally with
        near-perfect cc (e.g. GNW 1995 telemetry glitch: ~9,100 det/template,
        median cc ~0.98). Real tremor has high MAX but moderate MEDIAN counts and
        cc spread down to the threshold. Flag a day if its median per-template
        count or its median cc is implausibly high."""
        if qc_median_count_max is None and qc_median_cc_max is None:
            return False, 0, 0.0
        counts = [len(times) for times, _ in res.values()]
        all_cc = [c for _, ccs in res.values() for c in ccs]
        if not counts:
            return False, 0, 0.0
        med_count = int(np.median(counts))
        med_cc = float(np.median(all_cc)) if all_cc else 0.0
        bad = ((qc_median_count_max is not None and med_count > qc_median_count_max)
               or (qc_median_cc_max is not None and med_cc > qc_median_cc_max))
        return bad, med_count, med_cc

    def _flush():
        nonlocal rows, header_written
        if not rows:
            return
        df = pd.DataFrame(rows)
        df.to_csv(out_path, mode="a", header=not header_written, index=False)
        header_written = True
        rows = []

    # Submit ALL tasks up front and drain with as_completed. This keeps the
    # worker pool continuously saturated (workers never wait on the main thread
    # to hand them the next day), which is what makes this fast. An earlier
    # bounded-submission throttle was added to cap RAM, but that was only needed
    # because the (now-fixed) cc>1 bug flooded memory with bogus detections;
    # with the cc fix + per-day QC the detection volume is normal, the streaming
    # flush below keeps `rows` bounded, so submit-all is safe again.
    args_list = [
        (waveform_root, station, d, templates, fs, bandpass, threshold, min_gap_s)
        for d in days
    ]
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=mp_ctx) as ex:
        futs = [ex.submit(_worker_task, a) for a in args_list]
        for f in as_completed(futs):
            day, res = f.result()
            if res is not None:
                bad, med_count, med_cc = _is_bad_day(res)
                if bad:
                    dropped += 1
                    print(f"[{station}] QC DROP {day.date()}: median "
                          f"{med_count}/template, median cc {med_cc:.3f} "
                          f"(artifact day)", flush=True)
                    res = None
            if res is not None:
                for tkey, (times, ccs) in res.items():
                    for t, cc in zip(times, ccs):
                        rows.append({"template": tkey, "time": t, "cc": cc, "station": station})
                    if out_path is not None:
                        tmpl_counts[tkey] += len(times)
            done += 1
            if out_path is not None and done % flush_every == 0:
                _flush()
            if done % progress_every == 0:
                seen = (sum(tmpl_counts.values()) if out_path is not None else len(rows))
                print(f"[{station}] day {done}/{total}: cumulative {seen} "
                      f"detections ({dropped} days QC-dropped)", flush=True)

    if qc_median_count_max is not None or qc_median_cc_max is not None:
        print(f"[{station}] QC: dropped {dropped}/{done} days as artifact days",
              flush=True)
    if out_path is not None:
        _flush()
        return dict(tmpl_counts)
    return pd.DataFrame(rows)
