"""Inter-station cross-correlation driven by PNSN tremor windows.

The pipeline:

  1. From the PNSN tremor catalog, take every 5-minute event time as the
     start of a "tremor-active" window.
  2. For each (station-pair, tremor window): cut the same window from each
     station's continuous data, preprocess (bandpass + 1-bit + whitening),
     FFT-cross-correlate to get a CC trace from -max_lag to +max_lag.
  3. Stack CC traces over each downstream time bin (e.g. 2-day) into a single
     (pair, bin) CC stack. Write to HDF5.
  4. (Downstream) build a pre-ETS reference CC stack per pair and run
     `dvv.stretch_dvv` on the coda of the bin CC traces vs that reference.

This deliberately avoids the LFE-as-repeating-source assumption: tremor
windows are used as a quasi-noise wavefield emitted from a known patch of
the plate interface, so the inter-station CC approximates the Green's
function for paths illuminated by that source region.
"""

from __future__ import annotations

import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta
from itertools import combinations
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

DEFAULT_WINDOW_S = 300.0
DEFAULT_MAX_LAG_S = 60.0
DEFAULT_BANDPASS = (2.0, 8.0)
DEFAULT_FS = 40.0


def select_tremor_windows(
    pnsn: pd.DataFrame,
    t_start: datetime,
    t_end: datetime,
    bbox: tuple[float, float, float, float] | None = None,
    window_s: float = DEFAULT_WINDOW_S,
) -> list[tuple[datetime, datetime]]:
    """Each PNSN event time becomes the start of a tremor-active window of `window_s` seconds."""
    df = pnsn.copy()
    df["time"] = pd.to_datetime(df["time"])
    df = df[(df["time"] >= t_start) & (df["time"] < t_end)]
    if bbox is not None:
        lat_min, lat_max, lon_min, lon_max = bbox
        df = df[
            (df["lat"] >= lat_min) & (df["lat"] <= lat_max)
            & (df["lon"] >= lon_min) & (df["lon"] <= lon_max)
        ]
    windows = []
    for t in df["time"]:
        t_py = t.to_pydatetime() if hasattr(t, "to_pydatetime") else t
        windows.append((t_py, t_py + timedelta(seconds=window_s)))
    return windows


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
    return tr.data.astype(np.float32), tr.stats.starttime


def _preprocess_window(x: np.ndarray, fs: float, do_whiten: bool = False,
                       do_onebit: bool = False) -> np.ndarray | None:
    """Demean + per-window L2 normalize. Optional 1-bit / whitening (default off).

    For tremor inter-station CC, one-bit + whitening tend to destroy coherent
    direct-phase amplitude and replace it with noise; the empirical compare
    (`figures/smoke_tremor_cc_preprocess_compare.png`) shows raw-bandpass +
    L2-normalize gives by far the strongest CC peak.
    """
    if x.size == 0 or not np.isfinite(x).all():
        return None
    x = x - x.mean()
    nrm = np.linalg.norm(x)
    if nrm == 0:
        return None
    x = (x / nrm).astype(np.float32)
    if do_onebit:
        x = np.sign(x).astype(np.float32)
    if do_whiten:
        X = np.fft.rfft(x)
        mag = np.abs(X)
        mag[mag < 1e-6] = 1e-6
        X = X / mag
        x = np.fft.irfft(X, n=x.size).astype(np.float32)
    return x


def cross_correlate_window(a: np.ndarray, b: np.ndarray, fs: float,
                           max_lag_s: float = DEFAULT_MAX_LAG_S) -> np.ndarray | None:
    """Return CC trace for lags in [-max_lag_s, +max_lag_s], normalized."""
    if a.size != b.size or a.size == 0:
        return None
    a = a - a.mean()
    b = b - b.mean()
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    if norm_a == 0 or norm_b == 0:
        return None
    # FFT-based: cc[k] = sum_n a[n] * b[n+k]
    n = a.size
    pad = 1 << int(np.ceil(np.log2(2 * n)))
    A = np.fft.rfft(a, pad)
    B = np.fft.rfft(b, pad)
    cc_full = np.fft.irfft(np.conj(A) * B, pad)
    # Reorder to negative lags first, then positive
    cc_full = np.concatenate([cc_full[-(n - 1):], cc_full[:n]])
    cc_full /= (norm_a * norm_b)
    max_lag = int(round(max_lag_s * fs))
    center = n - 1
    return cc_full[center - max_lag : center + max_lag + 1].astype(np.float32)


def _bin_index(t: datetime, edges: list[datetime]) -> int | None:
    if t < edges[0] or t >= edges[-1]:
        return None
    for i in range(len(edges) - 1):
        if edges[i] <= t < edges[i + 1]:
            return i
    return None


def stack_pair(
    station_a: str,
    station_b: str,
    windows: list[tuple[datetime, datetime]],
    waveform_root: Path,
    bin_edges: list[datetime],
    fs: float = DEFAULT_FS,
    bandpass: tuple[float, float] = DEFAULT_BANDPASS,
    max_lag_s: float = DEFAULT_MAX_LAG_S,
) -> tuple[dict[int, np.ndarray], dict[int, int]]:
    """For one station pair, compute and stack CCs per bin.

    Returns (stacks, counts): stacks[bin_idx] = mean CC, counts[bin_idx] = n contributions.
    """
    from obspy import UTCDateTime

    n_lag_samples = int(round(max_lag_s * fs))
    cc_len = 2 * n_lag_samples + 1

    stacks_acc: dict[int, np.ndarray] = {}
    counts: dict[int, int] = {}

    # cache day streams per station to avoid repeated reads
    cache_a: tuple[datetime | None, np.ndarray | None, object] = (None, None, None)
    cache_b: tuple[datetime | None, np.ndarray | None, object] = (None, None, None)

    for t0, t1 in windows:
        bin_idx = _bin_index(t0, bin_edges)
        if bin_idx is None:
            continue
        day = datetime(t0.year, t0.month, t0.day)
        if cache_a[0] != day:
            loaded = _load_day_z(waveform_root, station_a, day, bandpass, fs)
            cache_a = (day, loaded[0] if loaded else None, loaded[1] if loaded else None)
        if cache_b[0] != day:
            loaded = _load_day_z(waveform_root, station_b, day, bandpass, fs)
            cache_b = (day, loaded[0] if loaded else None, loaded[1] if loaded else None)
        if cache_a[1] is None or cache_b[1] is None:
            continue

        # window may span two days; we keep simple: skip if either station starts later or ends earlier
        t0_utc = UTCDateTime(t0)
        t1_utc = UTCDateTime(t1)
        # data A
        start_a = cache_a[2]
        start_b = cache_b[2]
        end_a = start_a + (cache_a[1].size - 1) / fs
        end_b = start_b + (cache_b[1].size - 1) / fs
        if t0_utc < start_a or t1_utc > end_a or t0_utc < start_b or t1_utc > end_b:
            continue
        ia = int(round((t0_utc - start_a) * fs))
        ib = int(round((t0_utc - start_b) * fs))
        nsamp = int(round((t1_utc - t0_utc) * fs))
        sega = cache_a[1][ia : ia + nsamp]
        segb = cache_b[1][ib : ib + nsamp]
        if sega.size != nsamp or segb.size != nsamp:
            continue
        sega = _preprocess_window(sega, fs)
        segb = _preprocess_window(segb, fs)
        if sega is None or segb is None:
            continue
        cc = cross_correlate_window(sega, segb, fs, max_lag_s=max_lag_s)
        if cc is None or cc.size != cc_len:
            continue
        if bin_idx not in stacks_acc:
            stacks_acc[bin_idx] = np.zeros(cc_len, dtype=np.float32)
            counts[bin_idx] = 0
        stacks_acc[bin_idx] += cc
        counts[bin_idx] += 1

    # Convert accumulators to means
    stacks = {bi: arr / counts[bi] for bi, arr in stacks_acc.items() if counts[bi] > 0}
    return stacks, counts


def write_pair_hdf5(
    out_path: Path,
    station_pair: tuple[str, str],
    stacks: dict[int, np.ndarray],
    counts: dict[int, int],
    bin_edges: list[datetime],
    fs: float,
    max_lag_s: float,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out_path, "w") as h5:
        h5.attrs["station_a"] = station_pair[0]
        h5.attrs["station_b"] = station_pair[1]
        h5.attrs["fs"] = fs
        h5.attrs["max_lag_s"] = max_lag_s
        h5.attrs["bin_edges"] = np.array([e.timestamp() for e in bin_edges])
        for bi, arr in stacks.items():
            d = h5.create_dataset(f"bin_{bi:04d}", data=arr, compression="gzip")
            d.attrs["t_center"] = (bin_edges[bi] + (bin_edges[bi + 1] - bin_edges[bi]) / 2).timestamp()
            d.attrs["n_windows"] = counts.get(bi, 0)


def stack_all_pairs(
    stations: list[str],
    windows: list[tuple[datetime, datetime]],
    waveform_root: Path,
    bin_edges: list[datetime],
    out_dir: Path,
    fs: float = DEFAULT_FS,
    bandpass: tuple[float, float] = DEFAULT_BANDPASS,
    max_lag_s: float = DEFAULT_MAX_LAG_S,
    n_workers: int = 8,
) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    pairs = list(combinations(sorted(stations), 2))
    paths: list[Path] = []
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futures = {
            ex.submit(
                stack_pair, a, b, windows, waveform_root, bin_edges,
                fs, bandpass, max_lag_s,
            ): (a, b)
            for a, b in pairs
        }
        for f in as_completed(futures):
            a, b = futures[f]
            stacks, counts = f.result()
            if stacks:
                out_path = out_dir / f"{a}_{b}.h5"
                write_pair_hdf5(out_path, (a, b), stacks, counts, bin_edges, fs, max_lag_s)
                paths.append(out_path)
                log.info("wrote %s (%d bins)", out_path.name, len(stacks))
    return paths
