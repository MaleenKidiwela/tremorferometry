"""Synthetic recovery test for the stretching dv/v estimator."""

from __future__ import annotations

import numpy as np
import pytest

from tremorferometry.dvv import apply_stretch, stretch_dvv


def _synthetic_coda(fs: float = 100.0, duration: float = 40.0, seed: int = 0) -> np.ndarray:
    """Band-limited, exponentially decaying noise resembling a real coda."""
    rng = np.random.default_rng(seed)
    n = int(duration * fs)
    t = np.arange(n) / fs
    # white noise
    x = rng.standard_normal(n)
    # bandpass via FFT to 2-8 Hz
    X = np.fft.rfft(x)
    freqs = np.fft.rfftfreq(n, d=1.0 / fs)
    band = (freqs >= 2.0) & (freqs <= 8.0)
    X[~band] = 0.0
    x = np.fft.irfft(X, n=n)
    # coda decay envelope: zero before "S arrival" at t=2s, decay with Q~50 thereafter
    env = np.where(t < 2.0, 0.0, np.exp(-(t - 2.0) / 12.0))
    return x * env


@pytest.mark.parametrize("true_dvv", [-0.010, -0.005, 0.000, 0.003, 0.008])
def test_recovers_imposed_dvv(true_dvv: float) -> None:
    fs = 100.0
    ref = _synthetic_coda(fs=fs)
    cur = apply_stretch(ref, fs=fs, dvv=true_dvv)
    # window past the direct arrival
    result = stretch_dvv(ref, cur, fs=fs, t_min=4.0, t_max=20.0, eps_max=0.02, n_eps=401)
    assert result.cc_max > 0.9, f"cc too low: {result.cc_max}"
    assert abs(result.dvv - true_dvv) < 5e-4, (
        f"recovery error {result.dvv - true_dvv:.5f} > 5e-4 (true={true_dvv})"
    )


def test_flat_signal_returns_low_cc() -> None:
    fs = 100.0
    ref = _synthetic_coda(fs=fs, seed=1)
    cur = np.zeros_like(ref)
    with pytest.raises(ValueError):
        stretch_dvv(cur, cur, fs=fs, t_min=4.0, t_max=20.0)


def test_window_bounds_error() -> None:
    fs = 100.0
    ref = _synthetic_coda(fs=fs)
    with pytest.raises(ValueError):
        stretch_dvv(ref, ref, fs=fs, t_min=5.0, t_max=5.0)
