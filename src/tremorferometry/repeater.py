"""Properly discover repeating LFE waveforms at a station across time.

This is the Shelly-Beroza / Brown-Beroza-Shelly recipe adapted to our setup:

  1. For each candidate event (e.g. every Lin detection in a known LFE hotspot),
     cut a wide window from filtered continuous data, identify the envelope
     peak inside the expected direct-phase arrival region, and cut a tight
     out-window centered on that peak. This removes Lin's OT-timing slop.
  2. L2-normalize each aligned waveform.
  3. Compute all-pairs cross-correlation with a small shift allowance via
     batched FFT (per-row IFFT). Stack across stations for network CC.
  4. Threshold high to suppress random matches (>= 0.7 for confirmed repeats).
  5. Cluster matched pairs into families via transitive closure on the match
     graph; each family is candidate evidence for a repeating LFE patch.

Methodological notes (lessons logged in METHODS.md):
- Pairwise CC of broad band-limited noise gives ~0.5 baseline; thresholds
  below ~0.65 admit many noise coincidences. Use >= 0.7.
- Envelope-peak alignment is essential -- raw OT alignment lets ~1-2 s of
  jitter scatter the direct phase across the window and kill CC.
- We use single station (or pair) here; "same signal at same station over
  time" is the criterion. Verification across more stations is downstream.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import hilbert


def _load_day_filt(
    root: Path, station: str, day: datetime, bandpass: tuple[float, float], fs: float,
):
    """Load + filter one day at one station to component Z. Returns (data, start_UTC)."""
    from obspy import read

    candidates = list(root.glob(f"*.{station}/{day.year}/{day.timetuple().tm_yday:03d}.mseed"))
    if not candidates:
        return None
    try:
        st = read(str(candidates[0]))
    except Exception:
        return None
    st = st.select(component="Z")
    if len(st) == 0:
        return None
    # Normalize sampling rate per trace before merging -- some files have
    # slight FS jitter (e.g. 99.999... vs 100.0) across segments which
    # otherwise breaks st.merge.
    for tr in st:
        if abs(tr.stats.sampling_rate - fs) > 1e-6:
            try:
                tr.resample(fs)
            except Exception:
                return None
    try:
        st.merge(fill_value=0)
    except Exception:
        return None
    if len(st) == 0:
        return None
    st.detrend("demean")
    st.filter("bandpass", freqmin=bandpass[0], freqmax=bandpass[1], corners=4, zerophase=True)
    tr = st[0]
    return tr.data.astype(np.float32), tr.stats.starttime


def cut_aligned_window(
    data: np.ndarray, start_utc, t_ot, fs: float,
    search_window: tuple[float, float] = (5.0, 17.0),
    out_window: tuple[float, float] = (-1.0, 1.0),
) -> tuple[np.ndarray, float] | tuple[None, None]:
    """Find envelope-peak inside [t_ot+search_window], cut +/- out_window around it.

    Returns (L2-normalized waveform, peak_time_relative_to_OT) or (None, None).
    """
    n_search = int(round((search_window[1] - search_window[0]) * fs))
    n_out = int(round((out_window[1] - out_window[0]) * fs))
    t_start = t_ot + search_window[0]
    t_end = t_ot + search_window[1]
    if t_start < start_utc:
        return None, None
    fs_int = float(fs)
    end_utc = start_utc + (data.size - 1) / fs_int
    if t_end > end_utc:
        return None, None
    i_search0 = int(round((t_start - start_utc) * fs_int))
    search_seg = data[i_search0 : i_search0 + n_search]
    if search_seg.size != n_search:
        return None, None
    env = np.abs(hilbert(search_seg))
    i_peak = int(np.argmax(env))
    peak_time = search_window[0] + i_peak / fs_int  # relative to OT

    # Cut out_window around peak; bounds-check
    i_out0 = i_search0 + i_peak + int(round(out_window[0] * fs_int))
    if i_out0 < 0 or i_out0 + n_out > data.size:
        return None, None
    seg = data[i_out0 : i_out0 + n_out].copy()
    seg = seg - seg.mean()
    nrm = float(np.linalg.norm(seg))
    if nrm == 0 or not np.isfinite(nrm):
        return None, None
    return (seg / nrm).astype(np.float32), peak_time


def cut_all_detections(
    detections: pd.DataFrame,
    waveform_root: Path,
    station: str,
    fs: float = 40.0,
    bandpass: tuple[float, float] = (2.0, 8.0),
    search_window: tuple[float, float] = (5.0, 17.0),
    out_window: tuple[float, float] = (-1.0, 1.0),
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """For every detection, cut envelope-aligned window at the given station.

    Returns:
      X: array (n_ok, n_samples) of normalized windows (rows where successful)
      peak_times: array (n_ok,) of envelope-peak times relative to OT (seconds)
      df_ok: DataFrame of the surviving detections (same row order as X)
    """
    from obspy import UTCDateTime

    df = detections.sort_values("OT").reset_index(drop=True).copy()
    n_samples = int(round((out_window[1] - out_window[0]) * fs))
    n = len(df)
    X = np.zeros((n, n_samples), dtype=np.float32)
    peaks = np.zeros(n, dtype=np.float32)
    ok = np.zeros(n, dtype=bool)

    cur_day = None
    cur_loaded = None
    for i, row in df.iterrows():
        t = pd.Timestamp(row["OT"]).to_pydatetime()
        day = datetime(t.year, t.month, t.day)
        if day != cur_day:
            cur_loaded = _load_day_filt(waveform_root, station, day, bandpass, fs)
            cur_day = day
        if cur_loaded is None:
            continue
        data, start = cur_loaded
        t_ot = UTCDateTime(t)
        w, pk = cut_aligned_window(
            data, start, t_ot, fs,
            search_window=search_window, out_window=out_window,
        )
        if w is None:
            continue
        X[i] = w
        peaks[i] = pk
        ok[i] = True

    return X[ok], peaks[ok], df[ok].reset_index(drop=True)


def all_pairs_cc_max_shifted(
    X: np.ndarray, max_shift_samples: int = 20,
) -> np.ndarray:
    """All-pairs max-shifted CC on L2-normalized rows of X.

    Returns symmetric (n, n) matrix of max CC over lags in [-max_shift, +max_shift].
    Uses batched FFT-IFFT per row to avoid materializing the full O(n^2 * pad) tensor.
    """
    n, m = X.shape
    pad = 1 << int(np.ceil(np.log2(2 * m)))
    F = np.fft.rfft(X, n=pad, axis=1)
    # We'll accumulate per-row max over lags in [-s..s]
    # Use np.fft.irfft of (F[i].conj() * F) for each i
    out = np.zeros((n, n), dtype=np.float32)
    s = max_shift_samples
    for i in range(n):
        prod = np.conj(F[i:i + 1]) * F          # shape (1, n_freq) broadcast on F (n, n_freq)
        # Note: numpy broadcasting; (1,nf) * (n,nf) -> (n, nf)
        cc_full = np.fft.irfft(prod, n=pad, axis=1)
        # cc_full[j, k] is correlation X[i] vs X[j] at lag k (0..pad-1)
        # lag 0..s on the right, lag pad-s..pad-1 on the wraparound (negative lags)
        cc_pos = cc_full[:, : s + 1]
        cc_neg = cc_full[:, -s:]
        cc_max = np.maximum(cc_pos.max(axis=1), cc_neg.max(axis=1))
        out[i, :] = cc_max
    np.fill_diagonal(out, 0.0)
    # symmetrize (small floating-point asymmetries can happen)
    out = np.maximum(out, out.T)
    return out


def network_cc_all_pairs(
    Xs: list[np.ndarray], valid: list[np.ndarray], max_shift_samples: int = 20,
) -> np.ndarray:
    """Network CC across multiple stations, averaging where both detections have valid data.

    Xs: list of per-station (n, m) arrays (rows for invalid have all zeros)
    valid: list of per-station (n,) bool arrays indicating which rows are valid
    """
    n = Xs[0].shape[0]
    accum = np.zeros((n, n), dtype=np.float32)
    denom = np.zeros((n, n), dtype=np.float32)
    for X, v in zip(Xs, valid):
        cc = all_pairs_cc_max_shifted(X, max_shift_samples=max_shift_samples)
        pair_valid = np.outer(v, v).astype(np.float32)
        accum += cc * pair_valid
        denom += pair_valid
    out = np.where(denom > 0, accum / np.maximum(denom, 1), 0.0)
    np.fill_diagonal(out, 0.0)
    return out


def cluster_matches(network_cc: np.ndarray, threshold: float) -> np.ndarray:
    """Cluster events by transitive closure of (CC > threshold) match graph.

    Returns labels array of shape (n,). -1 for singletons (no match above threshold).
    Cluster IDs start at 0.
    """
    n = network_cc.shape[0]
    adj = (network_cc >= threshold).astype(np.int8)
    np.fill_diagonal(adj, 0)
    labels = np.full(n, -1, dtype=np.int64)
    cluster_id = 0
    # connected components via BFS
    visited = np.zeros(n, dtype=bool)
    for start in range(n):
        if visited[start]:
            continue
        if not adj[start].any():
            visited[start] = True
            continue
        queue = [start]
        visited[start] = True
        members = []
        while queue:
            v = queue.pop()
            members.append(v)
            neighbors = np.where(adj[v] & ~visited)[0]
            for nb in neighbors:
                visited[nb] = True
                queue.append(int(nb))
        if len(members) >= 2:
            for m in members:
                labels[m] = cluster_id
            cluster_id += 1
    return labels
