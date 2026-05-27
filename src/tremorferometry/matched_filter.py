"""Sliding-window matched filter for finding template hits in continuous data.

Given a template and a continuous trace at the same sample rate, compute the
normalized cross-correlation at every possible lag and return peaks above a
threshold. The standard tool for finding repeating-event family members
without relying on a pre-existing trigger list.

Network matched filter: do this independently at each station, then keep
detections where the network-summed CC exceeds the threshold (i.e. the
template hits at the same time at multiple stations).
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.signal import fftconvolve


def normalized_cc(continuous: np.ndarray, template: np.ndarray) -> np.ndarray:
    """Return normalized CC trace of `template` sliding through `continuous`.

    Output length = len(continuous) - len(template) + 1.
    Each value is the Pearson CC of the template against the window starting
    at that index.
    """
    n = len(continuous)
    m = len(template)
    if m > n:
        return np.zeros(0, dtype=np.float32)

    template = template - template.mean()
    nrm_t = float(np.linalg.norm(template))
    if nrm_t == 0:
        return np.zeros(n - m + 1, dtype=np.float32)
    template = template / nrm_t

    # Numerator: fft-based cross-correlation
    num = fftconvolve(continuous, template[::-1], mode="valid")  # length n-m+1

    # Sliding window mean and sum-of-squares for normalization
    c1 = np.concatenate(([0], np.cumsum(continuous.astype(np.float64))))
    win_sum = c1[m:] - c1[:-m]
    win_mean = win_sum / m

    c2 = np.concatenate(([0], np.cumsum((continuous.astype(np.float64)) ** 2)))
    win_sumsq = c2[m:] - c2[:-m]
    win_var = win_sumsq - m * win_mean ** 2
    # numerical floor
    win_var = np.maximum(win_var, 1e-12)
    win_std = np.sqrt(win_var)

    # Numerator was correlation against zero-mean unit-norm template;
    # The window's effective amplitude is win_std. So CC = num / win_std.
    cc = num / win_std
    return cc.astype(np.float32)


def pick_peaks(cc: np.ndarray, fs: float, threshold: float,
               min_gap_s: float = 6.0) -> np.ndarray:
    """Local-maximum peak picker on a CC trace with minimum separation.

    Returns sample indices of picks (each is the local max within a
    `min_gap_s` window where CC > threshold).
    """
    n_gap = max(1, int(round(min_gap_s * fs)))
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
    return np.array(picks, dtype=np.int64)


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


def scan_day(
    waveform_root: Path,
    station: str,
    day: datetime,
    template: np.ndarray,
    fs: float = 40.0,
    bandpass: tuple[float, float] = (2.0, 8.0),
    threshold: float = 0.5,
    min_gap_s: float = 6.0,
) -> tuple[list[datetime], list[float]] | tuple[None, None]:
    """Apply one template to one station-day; return (peak_times, peak_ccs)."""
    loaded = _load_day_filt(waveform_root, station, day, bandpass, fs)
    if loaded is None:
        return None, None
    data, start = loaded
    cc = normalized_cc(data, template)
    if cc.size == 0:
        return [], []
    picks = pick_peaks(cc, fs, threshold, min_gap_s)
    times = [(start + (i / fs)).datetime for i in picks]
    ccs = [float(cc[i]) for i in picks]
    return times, ccs


def network_scan_days(
    waveform_root: Path,
    stations: list[str],
    days: Iterable[datetime],
    templates: dict[str, dict[str, np.ndarray]],
    fs: float = 40.0,
    bandpass: tuple[float, float] = (2.0, 8.0),
    threshold_per_station: float = 0.4,
    network_threshold: float = 0.5,
    min_gap_s: float = 6.0,
    time_match_tol_s: float = 1.5,
) -> pd.DataFrame:
    """For each family in `templates` (dict mapping family_id -> {station: template}),
    scan each station-day and aggregate detections.

    A detection is recorded at time T if the SAME family template has a peak
    above `threshold_per_station` at every station within `time_match_tol_s`
    of T, AND the mean CC across stations exceeds `network_threshold`.

    Returns DataFrame with columns family_id, time, station_count, network_cc.
    """
    days = list(days)
    rows = []
    for family_id, tmpl_per_sta in templates.items():
        for day in days:
            # Get peaks per station
            station_peaks = {}
            for sta in stations:
                if sta not in tmpl_per_sta:
                    continue
                tmpl = tmpl_per_sta[sta]
                times, ccs = scan_day(
                    waveform_root, sta, day, tmpl, fs=fs,
                    bandpass=bandpass, threshold=threshold_per_station,
                    min_gap_s=min_gap_s,
                )
                if times is None:
                    continue
                station_peaks[sta] = list(zip(times, ccs))
            if len(station_peaks) < 2:
                continue
            # Aggregate: find times where peaks appear at >=2 stations within tolerance
            # Use the first station's picks as anchors
            ref_sta = sorted(station_peaks.keys())[0]
            ref_picks = station_peaks[ref_sta]
            other_stas = [s for s in station_peaks if s != ref_sta]
            for ref_time, ref_cc in ref_picks:
                cc_list = [ref_cc]
                hits = {ref_sta: ref_cc}
                for sta in other_stas:
                    best = None
                    for t, cc in station_peaks[sta]:
                        if abs((t - ref_time).total_seconds()) <= time_match_tol_s:
                            if best is None or cc > best:
                                best = cc
                    if best is not None:
                        cc_list.append(best)
                        hits[sta] = best
                if len(hits) >= 2:
                    ncc = float(np.mean(cc_list))
                    if ncc >= network_threshold:
                        rows.append({
                            "family_id": family_id,
                            "time": ref_time,
                            "n_stations": len(hits),
                            "network_cc": ncc,
                            **{f"cc_{s}": float(v) for s, v in hits.items()},
                        })
    return pd.DataFrame(rows)
