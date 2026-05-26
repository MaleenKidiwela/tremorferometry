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


def detect_episode(
    tremor: pd.DataFrame,
    rate_window_hours: float = 24.0,
    rate_threshold: int = 20,
    bbox_padding_deg: float = 0.2,
) -> TremorEpisode:
    """Find one ETS episode: contiguous span where tremor rate exceeds threshold.

    rate is detections per `rate_window_hours`. bbox covers all tremor detections
    inside the active span, padded by `bbox_padding_deg`.
    """
    if tremor.empty:
        raise ValueError("tremor catalog is empty")
    s = tremor.set_index("time").sort_index()
    counts = s.resample(f"{rate_window_hours}h")["lat"].count()
    active = counts >= rate_threshold
    if not active.any():
        raise RuntimeError("no time window exceeded the tremor rate threshold")

    # take the first contiguous run of `active`
    first = active.idxmax()
    after_first = active.loc[first:]
    end = after_first[~after_first].index.min() if (~after_first).any() else after_first.index.max()
    t_start = first.to_pydatetime()
    t_end = end.to_pydatetime()

    inside = tremor[(tremor["time"] >= t_start) & (tremor["time"] < t_end)]
    if inside.empty:
        raise RuntimeError("episode interval is empty")
    bbox = (
        float(inside["lat"].min()) - bbox_padding_deg,
        float(inside["lat"].max()) + bbox_padding_deg,
        float(inside["lon"].min()) - bbox_padding_deg,
        float(inside["lon"].max()) + bbox_padding_deg,
    )
    return TremorEpisode(t_start=t_start, t_end=t_end, bbox=bbox, n_detections=len(inside))
