"""Bostock LFE catalog ingestion and spatial filtering.

Bostock et al. (2012, 2015) publish per-family information as supplementary
material — typically a `families.txt` table (id, lat, lon, depth, n_detections)
and per-family template waveforms (multi-station, multi-channel) plus a list of
detection times.

This module is intentionally I/O-tolerant: the exact format of the supplement
varies between papers and personal-communication packages. `load_families`
expects a normalized CSV after `scripts/02_ingest_bostock.py` has done the
parsing.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

FAMILY_COLUMNS = ("family_id", "lat", "lon", "depth_km", "n_detections")


def load_families(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = set(FAMILY_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(f"family CSV missing columns {missing}")
    return df


def filter_by_bbox(
    families: pd.DataFrame, bbox: tuple[float, float, float, float]
) -> pd.DataFrame:
    lat_min, lat_max, lon_min, lon_max = bbox
    mask = (
        (families["lat"] >= lat_min)
        & (families["lat"] <= lat_max)
        & (families["lon"] >= lon_min)
        & (families["lon"] <= lon_max)
    )
    return families.loc[mask].reset_index(drop=True)
