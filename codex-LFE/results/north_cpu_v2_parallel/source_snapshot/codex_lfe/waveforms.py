"""Gap-aware native-coordinate processing with an explicit common UTC sample grid.

Matching native components to templates from the SAME channels needs no guessed
Z12->ZNE rotation. Geographic rotations / response removal are separate prerequisites
for cross-instrument transfer or a ZNE neural picker, not silently assumed here.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
from scipy.signal import butter, detrend, resample_poly, sosfiltfilt

from .provenance import digest, sha256, write_json


@dataclass
class Waveform:
    x: np.ndarray                 # channels x samples, zeros outside mask
    valid: np.ndarray             # channels x samples
    start_s: float
    fs: float
    channels: tuple[str, ...]     # FULL seed IDs, explicit location + channel
    site_id: str


def runs(mask):
    changes = np.diff(np.r_[False, np.asarray(mask, bool), False].astype(np.int8))
    return zip(np.flatnonzero(changes == 1), np.flatnonzero(changes == -1))


def invalidate_flat_runs(data, valid, min_samples):
    valid = valid.copy()
    for lo, hi in runs(np.r_[False, np.diff(data) == 0] & valid):
        if hi - lo + 1 >= min_samples:
            valid[max(0, lo - 1):hi] = False
    return valid


def preprocess_stream(stream, channels, site_id, start_s, duration_s, fs=20.0,
                      band=(1.0, 8.0), edge_s=10.0):
    """Conservatively mask gaps/conflicting overlaps and filter each valid run.

    Resample each run with anti-alias filtering, then interpolate its fractional
    timestamp to a fixed UTC grid. Fractional sample alignment is preserved, not
    rounded away. Linear interpolation is a documented first implementation.
    Native rate changes or sub-sample phase changes within one channel are rejected.
    """
    from obspy import UTCDateTime
    if len(channels) != 3 or len(set(channels)) != 3:
        raise ValueError("exactly three distinct full channel IDs required")
    if not (0 < band[0] < band[1] < fs / 2) or duration_s <= 0:
        raise ValueError("invalid passband or duration")
    n = int(round(duration_s * fs))
    out = np.zeros((3, n), np.float32)
    mask = np.zeros((3, n), bool)
    for c, channel in enumerate(channels):
        selected = stream.select(id=channel).copy()
        selected.trim(UTCDateTime(start_s - edge_s - 2), UTCDateTime(start_s + duration_s + edge_s + 2))
        if not selected:
            continue
        rates = {float(tr.stats.sampling_rate) for tr in selected}
        if len(rates) != 1:
            raise ValueError(f"rate epoch changes in {channel}")
        rate = rates.pop()
        origin = float(selected[0].stats.starttime)
        for tr in selected:
            offset = (float(tr.stats.starttime) - origin) * rate
            if abs(offset - round(offset)) > 1e-3:
                raise ValueError(f"sample phase changes in {channel}; split epochs first")
        # method=0 masks conflicting overlaps, fill_value=None preserves missing data.
        selected.merge(method=0, fill_value=None)
        if len(selected) != 1:
            raise ValueError(f"unresolved channel epochs in {channel}")
        tr = selected[0]
        a = np.ma.asarray(tr.data, dtype=np.float64)
        data = np.asarray(a.filled(0), dtype=np.float64)
        good = ~np.ma.getmaskarray(a) & np.isfinite(data)
        good = invalidate_flat_runs(data, good, max(2, int(rate)))
        if band[1] >= rate / 2:
            raise ValueError("native sampling rate too low")
        ratio = Fraction(fs / rate).limit_denominator(10000)
        if abs(float(ratio) - fs / rate) > 1e-10:
            raise ValueError("unsupported sample-rate ratio")
        sos = butter(4, band, fs=rate, btype="bandpass", output="sos")
        for lo, hi in runs(good):
            if hi - lo < (2 * edge_s + 2) * rate:
                continue
            y = sosfiltfilt(sos, detrend(data[lo:hi], type="linear"))
            y = resample_poly(y, ratio.numerator, ratio.denominator)
            seg_start = float(tr.stats.starttime) + lo / rate
            lower = seg_start + edge_s
            upper = seg_start + (hi - lo - 1) / rate - edge_s
            i0 = max(0, int(np.ceil((lower - start_s) * fs)))
            i1 = min(n, int(np.floor((upper - start_s) * fs)) + 1)
            if i1 <= i0:
                continue
            source_index = ((start_s - seg_start) * fs) + np.arange(i0, i1)
            out[c, i0:i1] = np.interp(source_index, np.arange(len(y)), y)
            mask[c, i0:i1] = True
    return Waveform(out, mask, float(start_s), float(fs), tuple(channels), site_id)


def load_day(path, cache_dir, channels, site_id, start_s, fs=20.0, duration_s=86400, edge_s=10.0):
    """Cache key includes raw bytes, code bytes and preprocessing parameters."""
    from obspy import read
    path = Path(path)
    config = dict(source=str(path.resolve()), source_sha256=sha256(path), channels=list(channels),
                  site_id=site_id, start_s=float(start_s), duration_s=duration_s, fs=fs,
                  edge_s=edge_s, band=[1.0, 8.0], code_sha256=sha256(__file__),
                  coordinates="native", units="instrument_counts_filtered")
    key = digest(config)
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = cache_dir / f"{key}.npz"
    if dest.exists():
        with np.load(dest, allow_pickle=False) as z:
            return Waveform(z["x"], z["valid"], start_s, fs, tuple(channels), site_id), config
    wave = preprocess_stream(read(str(path)), channels, site_id, start_s, duration_s, fs, edge_s=edge_s)
    tmp = dest.with_suffix(".tmp.npz")
    np.savez(tmp, x=wave.x, valid=wave.valid)
    tmp.replace(dest)
    write_json(dest.with_suffix(".json"), config | {"valid_fraction": wave.valid.mean(axis=1).tolist()})
    return wave, config
