"""Detect candidate LFE-like events from continuous data inside tremor windows.

Used when there is no pre-existing LFE catalog for a region. The PNSN tremor
catalog gives 5-minute windows where tremor is active anywhere in Cascadia;
inside each such window we look for local envelope peaks at a station,
treating each peak as a candidate detection.

The candidate list then feeds `repeater.py` exactly as Lin's catalog did for
southern V.I., so the rest of the family-discovery pipeline is unchanged.
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import hilbert


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


def envelope_peaks_in_windows(
    waveform_root: Path,
    station: str,
    tremor_windows: list[tuple[datetime, datetime]],
    fs: float = 40.0,
    bandpass: tuple[float, float] = (2.0, 8.0),
    snr_threshold: float = 3.0,
    min_peak_separation_s: float = 6.0,
) -> pd.DataFrame:
    """Inside each tremor window, find envelope peaks above SNR threshold.

    Returns DataFrame with columns ['time', 'envelope', 'background'] where
    `time` is the peak time (datetime), `envelope` is the peak amplitude, and
    `background` is the median envelope amplitude in the window.

    The SNR is defined as envelope_peak / median_envelope_in_window. We pick
    local maxima above `snr_threshold` (default 3.0) separated by at least
    `min_peak_separation_s` (default 6 s, slightly wider than an LFE).
    """
    from obspy import UTCDateTime

    rows = []
    cur_day = None
    cur_loaded = None
    n_sep = int(min_peak_separation_s * fs)
    for t0, t1 in tremor_windows:
        day = datetime(t0.year, t0.month, t0.day)
        if day != cur_day:
            cur_loaded = _load_day_filt(waveform_root, station, day, bandpass, fs)
            cur_day = day
        if cur_loaded is None:
            continue
        data, start = cur_loaded
        t0_utc = UTCDateTime(t0)
        t1_utc = UTCDateTime(t1)
        end_utc = start + (data.size - 1) / fs
        if t0_utc < start or t1_utc > end_utc:
            continue
        i0 = int(round((t0_utc - start) * fs))
        i1 = int(round((t1_utc - start) * fs))
        seg = data[i0:i1]
        if seg.size < int(min_peak_separation_s * fs * 2):
            continue
        env = np.abs(hilbert(seg))
        bg = float(np.median(env))
        if bg == 0:
            continue
        threshold = bg * snr_threshold
        # find all peaks above threshold; pick local maxima with separation
        above = env > threshold
        # naive peak picking: where derivative changes sign and value above threshold
        if not above.any():
            continue
        # scan in a sliding-max fashion
        i = 0
        picks = []
        while i < env.size:
            if env[i] >= threshold:
                # find max in next 2*n_sep window
                j_end = min(env.size, i + n_sep)
                j = int(np.argmax(env[i:j_end])) + i
                picks.append(j)
                i = j + n_sep
            else:
                i += 1
        for j in picks:
            t = (start + j / fs).datetime
            rows.append({"time": t, "envelope": float(env[j]), "background": bg,
                         "snr": float(env[j] / bg)})
    return pd.DataFrame(rows)
