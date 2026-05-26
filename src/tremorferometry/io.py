"""Schema helpers for the on-disk parquet / HDF5 artifacts."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

DETECTIONS_COLUMNS = ("family_id", "station", "channel", "time", "cc", "shift")
DVV_COLUMNS = ("family_id", "station", "t_center", "dvv", "dvv_err", "cc_max", "n_det")


def write_detections(df: pd.DataFrame, path: str | Path) -> None:
    missing = set(DETECTIONS_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(f"detections df missing columns {missing}")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)


def write_dvv(df: pd.DataFrame, path: str | Path) -> None:
    missing = set(DVV_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(f"dvv df missing columns {missing}")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)


def read_parquet(path: str | Path) -> pd.DataFrame:
    return pd.read_parquet(path)
