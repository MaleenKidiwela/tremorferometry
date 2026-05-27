"""Waveform-similarity reclustering of LFE detections.

The 0.1-deg grid clustering used by `lin_catalog.cluster_into_families` is too
coarse: each cell mixes several LFE sub-sources whose waveforms don't align,
so the resulting "families" don't satisfy the repeating-source assumption
underlying coda-wave interferometry.

This module implements the standard Bostock/Shelly recipe: cut a short window
around the direct phase at each detection, cross-correlate every pair, and
hierarchically cluster on (1 - CC) distance with a `cc_threshold` cutoff.
Each cluster is a true "family" — same source patch, repeating waveform.

For Cascadia southern V.I.:
- direct S arrives ~13 s after OT at PGC for a ~30 km LFE -> cut window in
  [10, 18] s post-cut-start (with pre_s=5, post_s=35 frame)
- bandpass 2-8 Hz
- cc_threshold = 0.5 is the typical Bostock criterion
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


def _load_day_band(
    root: Path, station: str, day: datetime, bandpass: tuple[float, float], fs: float,
    component: str = "Z",
):
    """Open a station's day MSEED, filter, return the trace data + start UTC."""
    from obspy import read

    candidates = list(root.glob(f"*.{station}/{day.year}/{day.timetuple().tm_yday:03d}.mseed"))
    if not candidates:
        return None
    st = read(str(candidates[0]))
    st = st.select(component=component)
    if len(st) == 0:
        return None
    st.merge(fill_value=0)
    st.detrend("demean")
    st.filter("bandpass", freqmin=bandpass[0], freqmax=bandpass[1], corners=4, zerophase=True)
    if abs(st[0].stats.sampling_rate - fs) > 1e-6:
        st.resample(fs)
    tr = st[0]
    return tr.data, tr.stats.starttime


def extract_detection_windows(
    detections: pd.DataFrame,
    family_id: str,
    waveform_root: Path,
    station: str,
    win_s: tuple[float, float] = (10.0, 18.0),
    bandpass: tuple[float, float] = (2.0, 8.0),
    fs: float = 40.0,
    component: str = "Z",
) -> tuple[np.ndarray, pd.DataFrame]:
    """Return (X, det_df) where X has shape (n_dets, n_samples) of bandpassed
    waveforms cut at (OT-5 + win_s[0], OT-5 + win_s[1]).

    Skips detections with no available day file or short-cut situations.
    Each row is L2-normalized (zero mean, unit norm) ready for dot-product CC.
    """
    from obspy import UTCDateTime

    fam = detections[detections["family_id"] == family_id].sort_values("time").reset_index(drop=True)
    pre_s = 5.0
    n_samples = int(round((win_s[1] - win_s[0]) * fs))
    X = np.zeros((len(fam), n_samples), dtype=np.float32)
    good = np.zeros(len(fam), dtype=bool)
    cur_day = None
    cur_data = None
    cur_start = None
    for i, row in fam.iterrows():
        t = row["time"]
        day = datetime(t.year, t.month, t.day)
        if day != cur_day:
            loaded = _load_day_band(waveform_root, station, day, bandpass, fs, component)
            if loaded is None:
                cur_day = day
                cur_data = None
                continue
            cur_data, cur_start = loaded
            cur_day = day
        if cur_data is None:
            continue
        t_utc = UTCDateTime(t)
        # cut start = t - pre_s + win_s[0] = t + (win_s[0] - pre_s)
        cut_start = t_utc + (win_s[0] - pre_s)
        cut_end = t_utc + (win_s[1] - pre_s)
        if cut_start < cur_start or cut_end > cur_start + (cur_data.size - 1) / fs:
            continue
        i_start = int(round((cut_start - cur_start) * fs))
        seg = cur_data[i_start : i_start + n_samples]
        if seg.size != n_samples:
            continue
        # demean + L2-normalize for dot-product CC
        seg = seg - seg.mean()
        nrm = np.linalg.norm(seg)
        if nrm == 0 or not np.isfinite(nrm):
            continue
        X[i] = (seg / nrm).astype(np.float32)
        good[i] = True
    return X[good], fam[good].reset_index(drop=True)


def pairwise_cc(X: np.ndarray) -> np.ndarray:
    """Compute the full pairwise CC matrix using dot product on L2-normalized rows."""
    return X @ X.T


def hier_cluster_by_cc(cc: np.ndarray, cc_threshold: float = 0.5) -> np.ndarray:
    """Hierarchical (average linkage) cluster, cut at distance = 1 - cc_threshold."""
    from scipy.cluster.hierarchy import fcluster, linkage
    from scipy.spatial.distance import squareform

    dist = np.clip(1.0 - cc, 0.0, 2.0)
    np.fill_diagonal(dist, 0.0)
    cond = squareform(dist, checks=False)
    Z = linkage(cond, method="average")
    return fcluster(Z, t=1.0 - cc_threshold, criterion="distance")
