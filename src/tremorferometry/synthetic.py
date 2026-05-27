"""Synthetic LFE-coda data with a known dv/v pattern.

Used by `scripts/00_smoke_synthetic.py` and the end-to-end pytest. Generates
per-family HDF5 stack files in exactly the same schema as `stacking.py` would
produce on real data, so everything downstream of stacking can be exercised.

We deliberately skip the FDSN/template-matching upstream stages here: the
upstream code is independently testable, and the value of this synthetic is to
shake out integration bugs in the dv/v measurement + plotting end of the
pipeline against a known signal.
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import h5py
import numpy as np


def master_template(
    fs: float = 100.0,
    duration: float = 35.0,
    t_arrival: float = 2.0,
    freq_band: tuple[float, float] = (2.0, 8.0),
    coda_q: float = 12.0,
    seed: int = 0,
) -> np.ndarray:
    """A clean LFE-like trace: short direct arrival at `t_arrival` + decaying coda."""
    rng = np.random.default_rng(seed)
    n = int(round(duration * fs))
    t = np.arange(n) / fs

    # bandpassed white noise, then enveloped: zero before arrival, exp decay after.
    x = rng.standard_normal(n)
    X = np.fft.rfft(x)
    f = np.fft.rfftfreq(n, d=1.0 / fs)
    band = (f >= freq_band[0]) & (f <= freq_band[1])
    X[~band] = 0.0
    x = np.fft.irfft(X, n=n)
    env = np.where(t < t_arrival, 0.0, np.exp(-(t - t_arrival) / coda_q))

    # mild boost at the direct arrival to mimic the LFE pulse
    pulse_width = 0.4
    pulse = np.exp(-((t - t_arrival) ** 2) / (2 * pulse_width**2))
    pulse *= np.cos(2 * np.pi * 4.0 * (t - t_arrival))
    pulse[t < t_arrival - 3 * pulse_width] = 0.0

    return x * env + 0.6 * pulse * np.where(t >= t_arrival - 3 * pulse_width, 1.0, 0.0)


def stretch_waveform(x: np.ndarray, fs: float, dvv: float) -> np.ndarray:
    """Resample x onto a stretched time axis matching a homogeneous dv/v change."""
    from scipy.interpolate import CubicSpline

    n = x.size
    t = np.arange(n) / fs
    spline = CubicSpline(t, x, extrapolate=False)
    out = spline(t / (1.0 + dvv))
    return np.nan_to_num(out, copy=False)


def dvv_ets_pattern(
    t_center: datetime,
    ets_start: datetime,
    ets_end: datetime,
    dvv_min: float = -0.005,
    drop_days: float = 10.0,
    recover_days: float = 15.0,
) -> float:
    """Synthetic dv/v vs time: flat 0 -> ramp down before ETS -> hold -> ramp up after.

    Mirrors the canonical Mexico-style SSE dv/v: a few-tenths-of-a-% drop with
    onset slightly before the ETS, recovery afterward.
    """
    days_before_ets = (ets_start - t_center).total_seconds() / 86400.0
    days_after_ets = (t_center - ets_end).total_seconds() / 86400.0

    if t_center < ets_start:
        if days_before_ets >= drop_days:
            return 0.0
        # ramp from 0 to dvv_min over drop_days
        frac = 1.0 - days_before_ets / drop_days
        return dvv_min * frac
    if t_center <= ets_end:
        return dvv_min
    if days_after_ets >= recover_days:
        return 0.0
    frac = 1.0 - days_after_ets / recover_days
    return dvv_min * frac


def make_bin_edges(t_start: datetime, t_end: datetime, bin_days: int) -> list[datetime]:
    edges = []
    cur = t_start
    while cur <= t_end:
        edges.append(cur)
        cur += timedelta(days=bin_days)
    return edges


def write_synthetic_family(
    out_path: Path,
    family_id: str,
    stations: list[str],
    bin_edges: list[datetime],
    template: np.ndarray,
    fs: float,
    ets_start: datetime,
    ets_end: datetime,
    dvv_min: float = -0.005,
    noise_level: float = 0.02,
    n_det_per_bin: int = 50,
    seed: int = 0,
) -> Path:
    """Write one family's stack HDF5 with the imposed dv/v pattern baked in.

    Per (station, bin): stretch the template by the dv/v at the bin center,
    add small white noise to simulate residual after stacking n_det_per_bin
    detections, write to /STATION/bin_{idx:04d} with t_center attr.
    """
    rng = np.random.default_rng(seed)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    bin_centers = [(bin_edges[i] + (bin_edges[i + 1] - bin_edges[i]) / 2) for i in range(len(bin_edges) - 1)]

    with h5py.File(out_path, "w") as h5:
        h5.attrs["family_id"] = family_id
        h5.attrs["fs"] = fs
        h5.attrs["bin_edges"] = np.array([e.timestamp() for e in bin_edges])
        for sta in stations:
            g = h5.create_group(sta)
            for i, t_center in enumerate(bin_centers):
                dvv = dvv_ets_pattern(t_center, ets_start, ets_end, dvv_min=dvv_min)
                stretched = stretch_waveform(template, fs, dvv)
                noise = noise_level * rng.standard_normal(stretched.size) / np.sqrt(n_det_per_bin)
                stack = stretched + noise
                d = g.create_dataset(f"bin_{i:04d}", data=stack.astype(np.float32), compression="gzip")
                d.attrs["t_center"] = t_center.timestamp()
                d.attrs["n_det"] = n_det_per_bin
                d.attrs["dvv_true"] = float(dvv)
    return out_path
