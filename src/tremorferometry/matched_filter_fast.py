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
    st = read(str(candidates[0]))
    st = st.select(component="Z")
    if len(st) == 0:
        return None
    st.merge(fill_value=0)
    st.detrend("demean")
    st.filter("bandpass", freqmin=bandpass[0], freqmax=bandpass[1], corners=4, zerophase=True)
    if abs(st[0].stats.sampling_rate - fs) > 1e-6:
        st.resample(fs)
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
) -> pd.DataFrame:
    """Parallel: load each day once, run all templates."""
    from concurrent.futures import ProcessPoolExecutor, as_completed

    def task(day):
        return day, scan_day_multi(
            waveform_root, station, day, templates,
            fs=fs, bandpass=bandpass, threshold=threshold, min_gap_s=min_gap_s,
        )

    rows = []
    done = 0
    total = len(days)
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = {ex.submit(task, d): d for d in days}
        for f in as_completed(futs):
            day, res = f.result()
            if res is not None:
                for tkey, (times, ccs) in res.items():
                    for t, cc in zip(times, ccs):
                        rows.append({"template": tkey, "time": t, "cc": cc, "station": station})
            done += 1
            if done % progress_every == 0:
                import time
                print(f"[{station}] day {done}/{total}: cumulative {len(rows)} detections", flush=True)
    return pd.DataFrame(rows)
