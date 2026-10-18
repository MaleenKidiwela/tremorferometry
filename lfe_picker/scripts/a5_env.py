#!/usr/bin/env python
"""Envelope-onset + AIC picker — verbatim copy of the three functions of
lfe_extension/scripts/a5_eq_corrections.py (smooth, aic_pick, env_onset), which are themselves the
d3_stage4_locate recipe.  Copied so lfe_picker never imports from the read-only trees."""
from __future__ import annotations

import numpy as np


def smooth(x, n):
    return np.convolve(x, np.ones(n) / n, mode="same") if n >= 2 else x


def aic_pick(x, i0, i1):
    seg = np.asarray(x[i0:i1], float)
    n = seg.size
    if n < 12:
        return (i0 + i1) // 2
    aic = np.full(n, np.inf)
    for k in range(5, n - 5):
        va, vb = seg[:k].var(), seg[k:].var()
        if va > 0 and vb > 0:
            aic[k] = k * np.log(va) + (n - k) * np.log(vb)
    return i0 + (int(np.argmin(aic)) if np.isfinite(aic).any() else n // 2)


def env_onset(x, lo, hi, nnoise, fs):
    """Envelope-onset pick inside [lo, hi) -> (index, peak/noise-rms).  d3_stage4_locate recipe."""
    from scipy.signal import hilbert
    x = np.asarray(x, float)
    env = smooth(np.abs(hilbert(x - x.mean())), max(2, int(round(0.25 * fs))))
    noise = float(np.median(env[:nnoise])) if nnoise > 2 else 0.0
    nrms = float(np.sqrt(np.mean(x[:nnoise] ** 2))) if nnoise > 2 else 0.0
    lo, hi = max(lo, nnoise), min(hi, len(x) - 2)
    if hi <= lo or nrms <= 0:
        return None, 0.0
    j = lo + int(np.argmax(env[lo:hi]))
    snr = float(env[j] / nrms)
    lvl = noise + 0.25 * (env[j] - noise)
    k = j
    while k > max(lo - int(2.0 * fs), nnoise) and env[k] > lvl:
        k -= 1
    w = int(round(0.75 * fs))
    return aic_pick(x, max(k - w, nnoise), min(k + w, len(x) - 1)), snr
