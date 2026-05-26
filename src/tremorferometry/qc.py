"""Sanity checks on the dv/v result.

These run on the final parquet from `07_measure_dvv.py`. Each function returns
a (passed: bool, info: dict) tuple; the script aggregates into a small report.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def spatial_coherence(dvv: pd.DataFrame, min_cc: float = 0.5) -> tuple[bool, dict]:
    """For each (t_center) snapshot, do dv/v values across stations agree?

    Crude version: compute per-time-bin std-of-dv/v across stations within each
    family. If most bins have a tight spread (< 2x the median per-bin error),
    consider it coherent.
    """
    if dvv.empty:
        return False, {"reason": "empty"}
    grp = dvv.groupby(["family_id", "t_center"])
    stats = grp.agg(std=("dvv", "std"), err=("dvv_err", "median"), n=("dvv", "count"))
    stats = stats.dropna()
    coherent = (stats["std"] < 2.0 * stats["err"]) & (stats["n"] >= 2)
    frac = float(coherent.mean()) if len(coherent) else 0.0
    return frac > 0.5, {"coherent_fraction": frac, "n_bins": int(len(coherent))}


def detection_count_independence(dvv: pd.DataFrame) -> tuple[bool, dict]:
    """dv/v should not be a simple function of n_det. Report Pearson r."""
    if dvv.empty:
        return False, {"reason": "empty"}
    x = np.log(dvv["n_det"].clip(lower=1))
    y = dvv["dvv"]
    if x.std() == 0 or y.std() == 0:
        return True, {"r": 0.0}
    r = float(np.corrcoef(x, y)[0, 1])
    return abs(r) < 0.5, {"r": r}


def quality_filter(
    dvv: pd.DataFrame, min_cc: float = 0.6, min_n_det: int = 20
) -> pd.DataFrame:
    return dvv[(dvv["cc_max"] >= min_cc) & (dvv["n_det"] >= min_n_det)].copy()
