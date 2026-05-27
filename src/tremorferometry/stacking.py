"""Time-binned LFE-aligned waveform stacks.

For each (family, station) pair and each time bin, cut a fixed-length window
around every LFE detection (aligned on detect time, which approximates the S
arrival in Bostock's LFE templates) and average them. Output is an HDF5 file
per family with one group per station and one dataset per bin.

Parallelism: ProcessPoolExecutor across (family, station) pairs.
"""

from __future__ import annotations

import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def make_bin_edges(t_start: datetime, t_end: datetime, bin_days: int) -> list[datetime]:
    edges = []
    cur = t_start
    while cur <= t_end:
        edges.append(cur)
        cur += timedelta(days=bin_days)
    return edges


def stack_family_station(
    detections: pd.DataFrame,
    family_id: str,
    station: str,
    waveform_root: Path,
    bin_edges: list[datetime],
    pre_s: float = 2.0,
    post_s: float = 30.0,
    bandpass: tuple[float, float] = (2.0, 8.0),
    fs: float = 100.0,
) -> dict[int, np.ndarray]:
    """Return {bin_index: stacked waveform} for one (family, station) pair."""
    from obspy import read

    if "station" in detections.columns:
        fam_det = detections[
            (detections["family_id"] == family_id) & (detections["station"] == station)
        ]
    else:
        # Catalog-direct mode: detection times apply to every station, no per-
        # station filter. Used when stacking from Lin (2023) or any other LFE
        # catalog where origin times are the detection record.
        fam_det = detections[detections["family_id"] == family_id]
    if fam_det.empty:
        return {}

    stacks: dict[int, list[np.ndarray]] = {i: [] for i in range(len(bin_edges) - 1)}
    cur_day = None
    cur_stream = None
    for _, row in fam_det.sort_values("time").iterrows():
        t = row["time"]
        bin_idx = _bin_index(t, bin_edges)
        if bin_idx is None:
            continue
        day = datetime(t.year, t.month, t.day)
        if day != cur_day:
            cur_day = day
            cur_stream = _load_day_stream(waveform_root, station, day, bandpass, fs)
        if cur_stream is None:
            continue
        cut = _cut_window(cur_stream, t, pre_s, post_s, fs)
        if cut is not None:
            stacks[bin_idx].append(cut)

    return {i: np.mean(np.stack(arrs), axis=0) for i, arrs in stacks.items() if len(arrs) > 0}


def _bin_index(t: datetime | pd.Timestamp, edges: list[datetime]) -> int | None:
    t = t if isinstance(t, datetime) else t.to_pydatetime()
    if t < edges[0] or t >= edges[-1]:
        return None
    # binary search would be nicer; linear is fine for ~hundreds of bins
    for i in range(len(edges) - 1):
        if edges[i] <= t < edges[i + 1]:
            return i
    return None


def _load_day_stream(
    root: Path, station: str, day: datetime, bandpass: tuple[float, float], fs: float
):
    from obspy import read

    candidates = list(root.glob(f"*.{station}/{day.year}/{day.timetuple().tm_yday:03d}.mseed"))
    if not candidates:
        return None
    st = read(str(candidates[0]))
    st.merge(fill_value=0)
    st.detrend("demean")
    st.filter("bandpass", freqmin=bandpass[0], freqmax=bandpass[1], corners=4, zerophase=True)
    if abs(st[0].stats.sampling_rate - fs) > 1e-6:
        st.resample(fs)
    return st


def _cut_window(st, t: datetime | pd.Timestamp, pre: float, post: float, fs: float):
    from obspy import UTCDateTime

    t_utc = UTCDateTime(t)
    try:
        cut = st.slice(t_utc - pre, t_utc + post).copy()
    except Exception:  # noqa: BLE001
        return None
    if len(cut) == 0:
        return None
    tr = cut[0]
    n_expected = int(round((pre + post) * fs))
    if tr.stats.npts < n_expected - 2:
        return None
    return tr.data[:n_expected]


def write_stacks_hdf5(
    stacks_by_station: dict[str, dict[int, np.ndarray]],
    bin_edges: list[datetime],
    out_path: Path,
    family_id: str,
    fs: float,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(out_path, "w") as h5:
        h5.attrs["family_id"] = family_id
        h5.attrs["fs"] = fs
        h5.attrs["bin_edges"] = np.array([e.timestamp() for e in bin_edges])
        for station, bins in stacks_by_station.items():
            g = h5.create_group(station)
            for bin_idx, stack in bins.items():
                d = g.create_dataset(f"bin_{bin_idx:04d}", data=stack, compression="gzip")
                d.attrs["t_center"] = bin_edges[bin_idx].timestamp()


def stack_all_parallel(
    detections: pd.DataFrame,
    families: Iterable[str],
    stations: Iterable[str],
    waveform_root: Path,
    bin_edges: list[datetime],
    out_dir: Path,
    fs: float = 100.0,
    n_workers: int = 32,
) -> list[Path]:
    """Run `stack_family_station` over all pairs in parallel; write one HDF5 per family."""
    families = list(families)
    stations = list(stations)
    paths: list[Path] = []
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        for family_id in families:
            futures = {
                ex.submit(
                    stack_family_station,
                    detections,
                    family_id,
                    sta,
                    waveform_root,
                    bin_edges,
                ): sta
                for sta in stations
            }
            stacks_by_station: dict[str, dict[int, np.ndarray]] = {}
            for f in as_completed(futures):
                sta = futures[f]
                result = f.result()
                if result:
                    stacks_by_station[sta] = result
            if stacks_by_station:
                out_path = out_dir / f"{family_id}.h5"
                write_stacks_hdf5(stacks_by_station, bin_edges, out_path, family_id, fs)
                paths.append(out_path)
    return paths
