"""REDPy-style incremental network-CC family discovery for LFEs.

REDPy (Hotovec-Ellis & Jeffries) builds repeating-event families by, for each
new candidate event, cross-correlating against the template of every existing
family and assigning it to the best-matching family if `CC > threshold`, or
seeding a new family otherwise. The CC is computed at the *network* level —
summed across stations — which averages out station-specific noise that
sinks single-station CC.

This module adapts that approach for LFEs:
  * candidates come from Lin (2023)'s detection times (we don't re-detect),
  * for each candidate we cut a short window at every station,
  * pair-wise CC against a family template uses `numpy.correlate` with a
    small ±lag allowance to absorb sub-cell timing uncertainty,
  * the family template is the running mean of accepted member waveforms
    (per station).

The output is a label per detection: -1 for unclassifiable, else the family
ID. Families with too few members can be discarded.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.signal import fftconvolve


@dataclass
class _Family:
    template: dict[str, np.ndarray]  # station -> L2-normalized waveform
    member_idx: list[int] = field(default_factory=list)


def _load_day_z(
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
    return tr.data, tr.stats.starttime


def _cut_window(
    data: np.ndarray, start_utc, t_utc, win_s: tuple[float, float], fs: float
) -> np.ndarray | None:
    cut_start = t_utc + win_s[0]
    cut_end = t_utc + win_s[1]
    n_samples = int(round((win_s[1] - win_s[0]) * fs))
    if cut_start < start_utc or cut_end > start_utc + (data.size - 1) / fs:
        return None
    i_start = int(round((cut_start - start_utc) * fs))
    seg = data[i_start : i_start + n_samples]
    if seg.size != n_samples:
        return None
    seg = seg - seg.mean()
    nrm = np.linalg.norm(seg)
    if nrm == 0 or not np.isfinite(nrm):
        return None
    return (seg / nrm).astype(np.float32)


def _shifted_cc_max(a: np.ndarray, b: np.ndarray, max_shift: int) -> float:
    """Max Pearson CC of two L2-normalized 1-D signals over lags in ±max_shift."""
    cc_full = fftconvolve(a, b[::-1], mode="full")
    center = cc_full.size // 2
    lo = center - max_shift
    hi = center + max_shift + 1
    return float(cc_full[lo:hi].max())


def cut_all_stations(
    detections: pd.DataFrame,
    family_id: str,
    waveform_root: Path,
    stations: list[str],
    win_s: tuple[float, float] = (5.0, 15.0),
    bandpass: tuple[float, float] = (2.0, 8.0),
    fs: float = 40.0,
) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    """Cut windows for every detection in a 0.1° family, at every station.

    Returns (waveforms, det_df). waveforms is {station: array shape (n, L)}.
    det_df has the same `time` column as `detections` filtered to this family
    plus a `good` mask of detections successfully cut at >= 1 station.
    """
    from obspy import UTCDateTime

    fam = detections[detections["family_id"] == family_id].sort_values("time").reset_index(drop=True)
    n_samples = int(round((win_s[1] - win_s[0]) * fs))
    n_det = len(fam)
    waveforms = {sta: np.zeros((n_det, n_samples), dtype=np.float32) for sta in stations}
    good_per = {sta: np.zeros(n_det, dtype=bool) for sta in stations}

    for sta in stations:
        cur_day = None
        cur_data = None
        cur_start = None
        for i, row in fam.iterrows():
            t = row["time"]
            day = datetime(t.year, t.month, t.day)
            if day != cur_day:
                loaded = _load_day_z(waveform_root, sta, day, bandpass, fs)
                if loaded is None:
                    cur_day = day
                    cur_data = None
                    continue
                cur_data, cur_start = loaded
                cur_day = day
            if cur_data is None:
                continue
            t_utc = UTCDateTime(t)
            w = _cut_window(cur_data, cur_start, t_utc, win_s, fs)
            if w is None:
                continue
            waveforms[sta][i] = w
            good_per[sta][i] = True

    fam = fam.copy()
    fam["good"] = np.logical_or.reduce(list(good_per.values()))
    fam["n_stations_good"] = np.sum(np.stack(list(good_per.values())), axis=0)
    return waveforms, fam


def network_cluster_redpy(
    waveforms: dict[str, np.ndarray],
    cc_threshold: float = 0.5,
    min_stations: int = 3,
    max_shift_samples: int = 60,
) -> np.ndarray:
    """REDPy-style streaming network-CC clustering.

    For each detection in row order, compute network CC = mean over stations
    where both the candidate and the family template have a valid window, of
    the max shifted CC. Assign to the best-matching family if it passes
    `cc_threshold` and there are at least `min_stations` overlapping
    stations; else seed a new family.

    Returns labels of shape (n_det,) with values in {0, 1, ...} per family,
    or -1 for unclassifiable detections.
    """
    stations = list(waveforms.keys())
    n_det = next(iter(waveforms.values())).shape[0]
    labels = np.full(n_det, -1, dtype=np.int64)
    families: list[_Family] = []

    # precompute which detections have a valid window per station
    valid = {sta: np.any(waveforms[sta] != 0, axis=1) for sta in stations}

    for i in range(n_det):
        good_stations = [sta for sta in stations if valid[sta][i]]
        if not good_stations:
            continue
        candidate = {sta: waveforms[sta][i] for sta in good_stations}

        best_fam = -1
        best_cc = -np.inf
        for fid, fam in enumerate(families):
            overlap = [sta for sta in good_stations if sta in fam.template]
            if len(overlap) < min_stations:
                continue
            ccs = [
                _shifted_cc_max(candidate[sta], fam.template[sta], max_shift_samples)
                for sta in overlap
            ]
            mean_cc = float(np.mean(ccs))
            if mean_cc > best_cc:
                best_cc = mean_cc
                best_fam = fid

        if best_fam >= 0 and best_cc >= cc_threshold:
            labels[i] = best_fam
            fam = families[best_fam]
            fam.member_idx.append(i)
            # update template = running mean (re-L2-normalized) of members
            for sta in good_stations:
                if sta in fam.template:
                    new_t = (fam.template[sta] * (len(fam.member_idx) - 1) + candidate[sta]) / len(fam.member_idx)
                else:
                    new_t = candidate[sta]
                new_t = new_t - new_t.mean()
                nrm = np.linalg.norm(new_t)
                if nrm > 0:
                    fam.template[sta] = (new_t / nrm).astype(np.float32)
        else:
            # seed a new family
            labels[i] = len(families)
            families.append(_Family(template={sta: candidate[sta].copy() for sta in good_stations},
                                    member_idx=[i]))

    return labels


def family_sizes(labels: np.ndarray) -> pd.Series:
    """Return a Series of family sizes sorted descending, excluding -1."""
    s = pd.Series(labels)
    return s[s >= 0].value_counts().sort_values(ascending=False)
