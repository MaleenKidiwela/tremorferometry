"""Ingest the Lin (2023) LFE catalog into our family + detection schema.

Source: https://doi.org/10.5281/ZENODO.10016020
File:   EQloc_001_0.1_3_S.csv (~106 MB, 1.05 M detections, 2005-01..2017-02,
        southern Vancouver Island).

Schema of the raw CSV:
    starttime, OT, lon, lat, depth, residual, dt, N

Lin's catalog is per-event (no native family labels), so we cluster events
spatially into proto-families by grid binning at the configured cell size
(default 0.1 deg, ~10 km). Each cell with >= min_n events becomes one family.

Outputs of `ingest_lin_catalog`:
  - families  : DataFrame matching `catalog.FAMILY_COLUMNS`
                (family_id, lat, lon, depth_km, n_detections)
  - detections: DataFrame with columns (family_id, time)
                — `time` is the LFE origin time. station column is intentionally
                empty; downstream stacking can drop the same template times into
                each FDSN station's continuous data.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def load_raw(path: str | Path) -> pd.DataFrame:
    """Load the Lin Zenodo CSV. Returns columns: time, lat, lon, depth_km."""
    df = pd.read_csv(path, parse_dates=["OT"], usecols=["OT", "lat", "lon", "depth", "N", "residual"])
    df = df.rename(columns={"OT": "time", "depth": "depth_km"})
    # Lin reports depth as a negative number (below sea level); flip sign so it
    # matches the rest of our catalog convention (depth_km positive downward).
    df["depth_km"] = -df["depth_km"]
    return df


def filter_window(
    df: pd.DataFrame,
    t_start,
    t_end,
    bbox: tuple[float, float, float, float],
) -> pd.DataFrame:
    """Restrict to time + bbox window. bbox = lat_min, lat_max, lon_min, lon_max."""
    t_start = pd.Timestamp(t_start)
    t_end = pd.Timestamp(t_end)
    lat_min, lat_max, lon_min, lon_max = bbox
    mask = (
        (df["time"] >= t_start)
        & (df["time"] < t_end)
        & (df["lat"] >= lat_min)
        & (df["lat"] <= lat_max)
        & (df["lon"] >= lon_min)
        & (df["lon"] <= lon_max)
    )
    return df.loc[mask].reset_index(drop=True)


def cluster_into_families(
    df: pd.DataFrame,
    cell_deg: float = 0.1,
    min_n: int = 50,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Bin detections onto a lat/lon grid; each populated cell is one family.

    Returns (families, detections), where detections has a `family_id` column
    referring to families.family_id. Family IDs are zero-padded "L{n:04d}".
    """
    if df.empty:
        cols = ("family_id", "lat", "lon", "depth_km", "n_detections")
        return pd.DataFrame(columns=cols), pd.DataFrame(columns=("family_id", "time"))

    lat_bin = (df["lat"] / cell_deg).round().astype(int)
    lon_bin = (df["lon"] / cell_deg).round().astype(int)
    df = df.assign(lat_b=lat_bin, lon_b=lon_bin)

    sizes = df.groupby(["lat_b", "lon_b"]).size()
    kept = sizes[sizes >= min_n].index
    if len(kept) == 0:
        return cluster_into_families(df.drop(columns=["lat_b", "lon_b"]), cell_deg, max(1, min_n // 2))

    fams_centroid = (
        df.set_index(["lat_b", "lon_b"])
          .loc[kept]
          .groupby(level=[0, 1])
          .agg(lat=("lat", "mean"), lon=("lon", "mean"), depth_km=("depth_km", "mean"),
               n_detections=("time", "count"))
          .reset_index()
    )
    fams_centroid = fams_centroid.sort_values("n_detections", ascending=False).reset_index(drop=True)
    fams_centroid["family_id"] = [f"L{i:04d}" for i in range(len(fams_centroid))]
    families = fams_centroid[["family_id", "lat", "lon", "depth_km", "n_detections"]].copy()

    key_to_id = {(int(r.lat_b), int(r.lon_b)): r.family_id for r in fams_centroid.itertuples()}
    df_keep = df[df.set_index(["lat_b", "lon_b"]).index.isin(kept)].copy()
    df_keep["family_id"] = [
        key_to_id[(int(a), int(b))] for a, b in zip(df_keep["lat_b"], df_keep["lon_b"])
    ]
    detections = df_keep[["family_id", "time"]].reset_index(drop=True)
    return families, detections


def ingest_lin_catalog(
    csv_path: str | Path,
    t_start,
    t_end,
    bbox: tuple[float, float, float, float],
    cell_deg: float = 0.1,
    min_n: int = 50,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = load_raw(csv_path)
    win = filter_window(raw, t_start, t_end, bbox)
    return cluster_into_families(win, cell_deg=cell_deg, min_n=min_n)
