"""Signed, demeaned multichannel correlation with explicit valid-window masks."""
from __future__ import annotations

import numpy as np
from scipy.signal import fftconvolve, find_peaks


def window_sum(x, n):
    a = np.cumsum(np.r_[0.0, np.asarray(x, dtype=np.float64)])
    return a[n:] - a[:-n]


def normalized_cc(x, template, valid=None, backend="numpy", device="cpu"):
    """One coefficient per possible TEMPLATE START. Equal weight per channel.

    A coefficient is invalid if any channel has any missing sample in that window.
    Signed correlation is retained; a polarity reversal is not an LFE match.
    Float64 rolling variance avoids single-precision cancellation on day-long data.
    CPU and torch share exactly the same normalization and validity implementation.
    """
    x, template = np.asarray(x, float), np.asarray(template, float)
    if x.ndim != 2 or template.ndim != 2 or x.shape[0] != template.shape[0]:
        raise ValueError("channel-by-time arrays with matching channels required")
    n = template.shape[1]
    if n < 2 or x.shape[1] < n:
        raise ValueError("template needs at least 2 samples and cannot exceed data")
    if not np.isfinite(template).all():
        raise ValueError("nonfinite template")
    valid = np.isfinite(x) if valid is None else np.asarray(valid, bool) & np.isfinite(x)
    if valid.shape != x.shape:
        raise ValueError("valid mask shape mismatch")
    clean = np.where(valid, x, 0.0)
    t = template - template.mean(axis=1, keepdims=True)
    et = (t * t).sum(axis=1)
    if (et <= 0).any():
        raise ValueError("constant template channel")
    coefficients, masks = [], []
    for channel in range(x.shape[0]):
        data = clean[channel]
        if backend == "numpy":
            dot = fftconvolve(data, t[channel, ::-1], mode="valid")
        elif backend == "torch":
            import torch
            a = torch.as_tensor(data, dtype=torch.float64, device=device)
            b = torch.as_tensor(t[channel, ::-1].copy(), dtype=torch.float64, device=device)
            length = 1 << (len(data) + n - 2).bit_length()
            dot = torch.fft.irfft(torch.fft.rfft(a, n=length) * torch.fft.rfft(b, n=length), n=length)
            dot = dot[n - 1:len(data)].cpu().numpy()
        else:
            raise ValueError(f"unknown backend {backend}")
        sums = window_sum(data, n)
        squares = window_sum(data * data, n)
        energy = squares - sums * sums / n
        # Explicitly reject near-constant windows instead of amplifying numerical error.
        good = (window_sum(valid[channel], n) == n) & (energy > np.maximum(squares * 1e-10, 1e-30))
        cc = np.zeros(len(dot), float)
        cc[good] = dot[good] / np.sqrt(energy[good] * et[channel])
        coefficients.append(np.clip(cc, -1, 1))
        masks.append(good)
    good = np.all(masks, axis=0)
    cc = np.mean(coefficients, axis=0)
    cc[~good] = np.nan
    return cc, good


def align_to_anchor(cc, valid, offset_s, fs, n_output):
    """Station template starts at anchor_time + offset_s. No circular wrapping."""
    offset = offset_s * fs
    nearest = int(round(offset))
    if abs(offset - nearest) > 1e-6:
        raise ValueError("template moveout must be on the output sample grid")
    out = np.full(n_output, np.nan)
    good = np.zeros(n_output, bool)
    lo, hi = max(0, -nearest), min(n_output, len(cc) - nearest)
    if hi > lo:
        out[lo:hi] = cc[lo + nearest:hi + nearest]
        good[lo:hi] = valid[lo + nearest:hi + nearest]
    return out, good


def combine_sites(station_series, station_sites):
    """Equal mean within a site, then equal mean across sites.

    All declared station series must be valid: no untracked denominator changes.
    Stations at the same site cannot increase the site's total network weight.
    """
    groups = {}
    for station, series in station_series.items():
        if station not in station_sites:
            raise ValueError(f"missing site identity for {station}")
        groups.setdefault(station_sites[station], []).append(series)
    return {site: np.mean(values, axis=0) for site, values in groups.items()}


def shifted_sites(sites, shifts):
    """Independent site shifts, shared by ALL channels/families at that site.

    NaNs shift with the data. No wraparound: exposed edges stay invalid.
    Positive shift means later output; sample t reads original sample t-shift.
    """
    output = {}
    for site, series in sites.items():
        shift = int(shifts[site])
        x = np.asarray(series)
        y = np.full(x.shape, np.nan)
        if abs(shift) < len(x):
            if shift > 0:
                y[shift:] = x[:-shift]
            elif shift < 0:
                y[:shift] = x[-shift:]
            else:
                y[:] = x
        output[site] = y
    return output


def robust_scale(x):
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    if len(x) < 100:
        return None
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    return (med, mad) if mad > 1e-12 else None


def extract_candidates(sites, fs, start_s, threshold, support_mad=2.0, min_sites=3,
                       refractory_s=2.0, common_valid=None):
    """Same operation for observed and controls, including support and declustering.

    This function expects one normalization block, not a full multi-day record.
    Every accepted event requires positive independent-site support at the same
    moveout-aligned sample. Availability records include valid subthreshold samples.
    """
    if len(sites) < min_sites:
        return [], np.zeros(len(next(iter(sites.values()))) if sites else 0, bool)
    names = sorted(sites)
    x = np.asarray([sites[s] for s in names])
    good = np.isfinite(x).all(axis=0)
    if common_valid is not None:
        good &= common_valid
    stats = [robust_scale(row[good]) for row in x]
    scale = robust_scale(x[:, good].mean(axis=0))
    if scale is None or any(s is None for s in stats):
        return [], np.zeros(x.shape[1], bool)
    network = x.mean(axis=0)
    z = (network - scale[0]) / scale[1]
    z[~good] = -np.inf
    support = np.asarray([(row - s[0]) / s[1] for row, s in zip(x, stats)])
    # Gate first, so a strong single-site artifact cannot suppress a nearby event.
    gate = (support >= support_mad).sum(axis=0) >= min_sites
    gated = np.where(gate, z, -np.inf)
    peaks, _ = find_peaks(gated, height=threshold, distance=max(1, int(round(refractory_s * fs))))
    rows = []
    for i in peaks:
        rows.append({"time_s": float(start_s + i / fs), "score_mad": float(z[i]),
                     "network_cc": float(network[i]), "n_sites": len(names),
                     "n_supporting_sites": int(gate[i] and (support[:, i] >= support_mad).sum()),
                     "supporting_sites": ",".join(s for s, v in zip(names, support[:, i]) if v >= support_mad)})
    return rows, good


def resolve_candidates(candidates, tolerance_s=1.0, ambiguity_margin=0.5):
    """Greedy temporal event grouping on a common physical clock, with alternatives.

    Close competing families produce an UNASSIGNED event unless the score margin
    exceeds the configured bar. This is a conservative provisional resolution rule,
    not a validated family classifier. Competing family scores are retained for audit.
    """
    import pandas as pd
    if candidates.empty:
        return candidates.copy()
    if candidates.time_basis.nunique() != 1:
        raise ValueError("cannot group incompatible clocks")
    c = candidates.sort_values(["score_mad", "time_s", "family_id"], ascending=[False, True, True], kind="stable")
    remaining = np.ones(len(c), bool)
    times = c.time_s.to_numpy()
    families = c.family_id.to_numpy()
    scores = c.score_mad.to_numpy()
    records = c.to_dict('records')
    out = []
    for i in range(len(c)):
        if not remaining[i]:
            continue
        near = remaining & (np.abs(times - times[i]) <= tolerance_s)
        indices = np.flatnonzero(near)
        row = records[i].copy()
        others = indices[families[indices] != families[i]]
        margin = float(scores[i] - scores[others].max()) if len(others) else None
        row["family_alternatives"] = ",".join(sorted(set(families[indices])))
        row["family_score_margin"] = margin
        row["assignment_status"] = "ambiguous" if margin is not None and margin < ambiguity_margin else "provisional"
        if row["assignment_status"] == "ambiguous":
            row["family_id"] = None
        out.append(row)
        remaining[near] = False
    return pd.DataFrame(out).sort_values("time_s").reset_index(drop=True)
