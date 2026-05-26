"""Stretching method for dv/v on coda waveforms.

Forward model: a coda waveform recorded after a homogeneous fractional velocity
change dv/v is a time-stretched copy of the reference,

    cur(t) = ref( t / (1 + dv/v) )    equivalently    ref(t) = cur( t * (1 + dv/v) ).

To recover dv/v we scan a grid of trial stretches eps, resample the current trace
onto the reference's time axis assuming that stretch, and pick the eps that
maximizes correlation with the reference on the coda window.

Convention: a negative dv/v means the medium slowed down, so coda arrivals
appear later in `cur` than in `ref`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline


@dataclass(frozen=True)
class StretchResult:
    dvv: float
    cc_max: float
    cc_curve: np.ndarray
    eps_grid: np.ndarray
    dvv_err: float


def stretch_dvv(
    ref: np.ndarray,
    cur: np.ndarray,
    fs: float,
    t_min: float,
    t_max: float,
    eps_max: float = 0.02,
    n_eps: int = 201,
) -> StretchResult:
    """Estimate dv/v between reference and current waveforms by stretching.

    Parameters
    ----------
    ref, cur : 1D arrays, same length, same sample rate `fs`.
    t_min, t_max : coda window in seconds, relative to sample 0.
    eps_max : half-width of the stretch search (e.g. 0.02 = +/- 2 %).
    n_eps : number of trial stretches.

    Returns
    -------
    StretchResult with dvv (peak), cc_max, the full cc curve over eps_grid,
    and a parabolic-fit uncertainty estimate.
    """
    if ref.shape != cur.shape or ref.ndim != 1:
        raise ValueError("ref and cur must be matching 1D arrays")
    if fs <= 0 or t_max <= t_min:
        raise ValueError("invalid sampling rate or window")

    n = ref.size
    t = np.arange(n) / fs
    mask = (t >= t_min) & (t < t_max)
    if mask.sum() < 8:
        raise ValueError("coda window contains fewer than 8 samples")

    t_coda = t[mask]
    ref_coda = ref[mask]
    ref_demean = ref_coda - ref_coda.mean()
    ref_norm = np.sqrt(np.sum(ref_demean * ref_demean))
    if ref_norm == 0:
        raise ValueError("reference coda is identically zero")

    spline = CubicSpline(t, cur, extrapolate=False)
    eps_grid = np.linspace(-eps_max, eps_max, n_eps)
    cc_curve = np.empty(n_eps)

    for i, eps in enumerate(eps_grid):
        cur_stretched = spline(t_coda * (1.0 + eps))
        cur_stretched = np.nan_to_num(cur_stretched, copy=False)
        cur_demean = cur_stretched - cur_stretched.mean()
        denom = ref_norm * np.sqrt(np.sum(cur_demean * cur_demean))
        cc_curve[i] = float(np.dot(ref_demean, cur_demean) / denom) if denom > 0 else 0.0

    i_max = int(np.argmax(cc_curve))
    dvv_peak, dvv_err = _parabolic_refine(eps_grid, cc_curve, i_max)

    return StretchResult(
        dvv=dvv_peak,
        cc_max=float(cc_curve[i_max]),
        cc_curve=cc_curve,
        eps_grid=eps_grid,
        dvv_err=dvv_err,
    )


def _parabolic_refine(eps_grid: np.ndarray, cc: np.ndarray, i: int) -> tuple[float, float]:
    """Sub-grid peak by parabolic interpolation; width gives a crude uncertainty."""
    if i == 0 or i == len(cc) - 1:
        return float(eps_grid[i]), float(eps_grid[1] - eps_grid[0])
    y0, y1, y2 = cc[i - 1], cc[i], cc[i + 1]
    denom = y0 - 2.0 * y1 + y2
    if denom == 0:
        return float(eps_grid[i]), float(eps_grid[1] - eps_grid[0])
    delta = 0.5 * (y0 - y2) / denom
    step = float(eps_grid[1] - eps_grid[0])
    peak = float(eps_grid[i]) + delta * step
    curvature = abs(denom)
    err = step * float(np.sqrt(max(1.0 - y1, 0.0) / max(curvature, 1e-12)))
    return peak, err


def apply_stretch(x: np.ndarray, fs: float, dvv: float) -> np.ndarray:
    """Synthesize a stretched copy of x with the given dv/v (for tests / forward modelling)."""
    n = x.size
    t = np.arange(n) / fs
    t_query = t / (1.0 + dvv)
    spline = CubicSpline(t, x, extrapolate=False)
    out = spline(t_query)
    return np.nan_to_num(out, copy=False)
