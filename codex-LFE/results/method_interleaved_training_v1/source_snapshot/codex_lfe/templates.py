"""Reference-aligned templates from training days only, in native channel coordinates."""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter1d


def cut(waveform, start_s, n, alignment='linear'):
    if alignment not in ('linear', 'sinc'):
        raise ValueError('unknown template alignment')
    if alignment == 'sinc':
        from .resampling import fractional_samples
        position = (start_s - waveform.start_s) * waveform.fs
        first = int(np.floor(position))
        width = 32
        lo, hi = first - width, first + n + width + 1
        if lo < 0 or hi > waveform.x.shape[1]:
            return None
        if not waveform.valid[:, lo:hi].all():
            return None
        # Only convolve the padded window, not the entire day for every event.
        return np.asarray([fractional_samples(x[lo:hi], position - lo, n, width)
                           for x in waveform.x])
    idx = (start_s - waveform.start_s) * waveform.fs + np.arange(n)
    lo = np.floor(idx).astype(int)
    hi = np.ceil(idx).astype(int)
    if lo[0] < 0 or hi[-1] >= waveform.x.shape[1]:
        return None
    if not (waveform.valid[:, lo] & waveform.valid[:, hi]).all():
        return None
    fraction = idx - lo
    return waveform.x[:, lo] * (1 - fraction) + waveform.x[:, hi] * fraction


def stack_windows(windows):
    """One common RMS per event preserves native inter-component amplitude ratios."""
    w = np.asarray(windows, float)
    w = w - w.mean(axis=2, keepdims=True)
    rms = np.sqrt(np.mean(w * w, axis=(1, 2)))
    if (rms <= 0).any() or not np.isfinite(w).all():
        raise ValueError("invalid training window")
    return (w / rms[:, None, None]).mean(axis=0)


def vector_cc(a, b):
    a = a - a.mean(axis=1, keepdims=True)
    b = b - b.mean(axis=1, keepdims=True)
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    return float(np.sum(a * b) / denom) if denom else 0.0


def build_station_template(windows, fs, search_start_s=-10, template_s=8,
                           pre_peak_s=2, min_events=20, min_split_cc=0.25):
    """Choose the peak using half A, verify half B, then combine training halves.

    No per-event maximum-lag alignment: avoiding a hidden noise-selection step is
    more valuable than optimizing first-pilot sensitivity. Validation performance
    and absolute timestamp accuracy remain separate empirical questions.
    """
    if len(windows) < min_events:
        return None, {"usable": False, "reason": "insufficient_training_windows", "n_windows": len(windows)}
    a, b = stack_windows(windows[::2]), stack_windows(windows[1::2])
    energy = gaussian_filter1d(np.sum(a * a, axis=0), sigma=0.25 * fs)
    begin = int(round((0 - search_start_s) * fs))
    end = min(len(energy) - int(template_s * fs), int(round((50 - search_start_s) * fs)))
    peak = begin + int(np.argmax(energy[begin:end]))
    first = peak - int(round(pre_peak_s * fs))
    last = first + int(round(template_s * fs))
    split_cc = vector_cc(a[:, first:last], b[:, first:last])
    qc = {"usable": split_cc >= min_split_cc, "n_windows": len(windows),
          "split_stack_cc": split_cc, "peak_offset_from_published_s": search_start_s + peak / fs,
          "template_start_from_published_s": search_start_s + first / fs,
          "reason": "passed_training_split" if split_cc >= min_split_cc else "incoherent_training_halves"}
    template = stack_windows(windows)[:, first:last]
    return (template if qc["usable"] else None), qc
