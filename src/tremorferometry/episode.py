"""PNSN tremor catalog ingestion and episode definition.

PNSN serves the public Wech-style tremor catalog (envelope-correlation
locations, 5 min windows) at https://pnsn.org/tremor — typically as CSV
exports filterable by time range and bbox.

This module:
  * fetches or loads a tremor CSV for a date range
  * defines an ETS "episode" by clustering tremor in time (and optionally space)
  * returns a refined bbox + (t_start, t_end) you can write back into the config
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class TremorEpisode:
    t_start: datetime
    t_end: datetime
    bbox: tuple[float, float, float, float]  # lat_min, lat_max, lon_min, lon_max
    n_detections: int


def load_pnsn_tremor(path: str | Path) -> pd.DataFrame:
    """Load a PNSN tremor CSV. Expected columns: time, lat, lon, depth (optional)."""
    df = pd.read_csv(path, parse_dates=["time"])
    required = {"time", "lat", "lon"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"PNSN tremor CSV missing columns {missing}")
    return df.sort_values("time").reset_index(drop=True)


def _runs(active: pd.Series) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
    """Yield (start, end) of contiguous True runs in a bool series indexed by time."""
    edges = active.astype(int).diff().fillna(active.iloc[0].astype(int))
    starts = list(active.index[edges == 1])
    ends = list(active.index[edges == -1])
    if active.iloc[0]:
        starts = [active.index[0]] + starts
    if active.iloc[-1]:
        ends = ends + [active.index[-1]]
    # dedupe / pair
    return list(zip(starts, ends))


def list_episodes(
    tremor: pd.DataFrame,
    rate_window_hours: float = 24.0,
    rate_threshold: int = 20,
    bbox_padding_deg: float = 0.2,
    gap_tolerance_days: float = 2.0,
) -> list[TremorEpisode]:
    """Return every contiguous active span as a TremorEpisode, sorted by n_detections desc.

    `gap_tolerance_days`: short below-threshold gaps (e.g. weekends, instrument
    dropouts) are bridged so a single ETS isn't split into many sub-episodes.
    """
    if tremor.empty:
        return []
    s = tremor.set_index("time").sort_index()
    counts = s.resample(f"{rate_window_hours}h")["lat"].count()
    active = counts >= rate_threshold
    if not active.any():
        return []

    # bridge short gaps
    n_bridge = max(1, int(round(gap_tolerance_days * 24.0 / rate_window_hours)))
    bridged = active.copy()
    rolled = active.rolling(window=2 * n_bridge + 1, center=True, min_periods=1).max()
    bridged = (rolled > 0) & (counts.rolling(2 * n_bridge + 1, center=True, min_periods=1).sum() >= rate_threshold)

    runs = _runs(bridged)
    episodes: list[TremorEpisode] = []
    for t0, t1 in runs:
        t_start = t0.to_pydatetime()
        t_end = (t1 + pd.Timedelta(f"{rate_window_hours}h")).to_pydatetime()
        inside = tremor[(tremor["time"] >= t_start) & (tremor["time"] < t_end)]
        if inside.empty:
            continue
        bbox = (
            float(inside["lat"].min()) - bbox_padding_deg,
            float(inside["lat"].max()) + bbox_padding_deg,
            float(inside["lon"].min()) - bbox_padding_deg,
            float(inside["lon"].max()) + bbox_padding_deg,
        )
        episodes.append(
            TremorEpisode(t_start=t_start, t_end=t_end, bbox=bbox, n_detections=len(inside))
        )
    episodes.sort(key=lambda e: e.n_detections, reverse=True)
    return episodes


def detect_episode(
    tremor: pd.DataFrame,
    rate_window_hours: float = 24.0,
    rate_threshold: int = 20,
    bbox_padding_deg: float = 0.2,
    select: str = "largest",
    gap_tolerance_days: float = 2.0,
) -> TremorEpisode:
    """Pick one ETS episode out of the catalog.

    `select`:
      - "largest"  : episode with the most detections (default)
      - "first"    : earliest active span
      - "latest"   : most recent active span
    """
    episodes = list_episodes(
        tremor,
        rate_window_hours=rate_window_hours,
        rate_threshold=rate_threshold,
        bbox_padding_deg=bbox_padding_deg,
        gap_tolerance_days=gap_tolerance_days,
    )
    if not episodes:
        raise RuntimeError("no time window exceeded the tremor rate threshold")
    if select == "largest":
        return episodes[0]
    if select == "first":
        return min(episodes, key=lambda e: e.t_start)
    if select == "latest":
        return max(episodes, key=lambda e: e.t_start)
    raise ValueError(f"unknown select={select!r}")
