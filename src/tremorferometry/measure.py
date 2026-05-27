"""dv/v measurement from per-family HDF5 stacks.

This is the driver invoked by `scripts/07_measure_dvv.py`; it lives in the
package (rather than the script) so the smoke / integration tests can call
it directly.

For each family file the routine builds a reference stack from the bins whose
centers fall inside the reference window, then runs `stretch_dvv` on every bin
against that reference. Embarrassingly parallel across families with joblib.
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from .dvv import stretch_dvv


def _bins_with_centers(h5: h5py.File, station: str) -> list[tuple[int, float, np.ndarray]]:
    g = h5[station]
    out = []
    for key in sorted(g.keys()):
        ds = g[key]
        idx = int(key.split("_")[1])
        t_center = float(ds.attrs["t_center"])
        out.append((idx, t_center, ds[...]))
    return out


def _build_reference(
    bins: list[tuple[int, float, np.ndarray]],
    ref_start: datetime,
    ref_end: datetime,
) -> np.ndarray | None:
    t0, t1 = ref_start.timestamp(), ref_end.timestamp()
    members = [arr for _, tc, arr in bins if t0 <= tc < t1]
    if not members:
        return None
    L = min(a.size for a in members)
    return np.mean(np.stack([a[:L] for a in members]), axis=0)


def measure_family(
    h5_path: Path,
    coda_window: tuple[float, float],
    fs: float,
    eps_max: float,
    n_eps: int,
    ref_start: datetime,
    ref_end: datetime,
    min_cc: float,
) -> list[dict]:
    rows: list[dict] = []
    with h5py.File(h5_path, "r") as h5:
        family_id = str(h5.attrs.get("family_id", h5_path.stem))
        stations = list(h5.keys())
        for station in stations:
            bins = _bins_with_centers(h5, station)
            if not bins:
                continue
            ref = _build_reference(bins, ref_start, ref_end)
            if ref is None:
                continue
            for _, t_center, arr in bins:
                L = min(ref.size, arr.size)
                try:
                    res = stretch_dvv(
                        ref[:L], arr[:L], fs=fs,
                        t_min=coda_window[0], t_max=coda_window[1],
                        eps_max=eps_max, n_eps=n_eps,
                    )
                except ValueError:
                    continue
                if res.cc_max < min_cc:
                    continue
                rows.append({
                    "family_id": family_id,
                    "station": station,
                    "t_center": datetime.fromtimestamp(t_center),
                    "dvv": res.dvv,
                    "dvv_err": res.dvv_err,
                    "cc_max": res.cc_max,
                    "n_det": int(arr.size),
                })
    return rows


def measure_many(
    h5_paths: list[Path],
    coda_window: tuple[float, float],
    fs: float,
    eps_max: float,
    n_eps: int,
    ref_start: datetime,
    ref_end: datetime,
    min_cc: float,
    n_workers: int = 1,
) -> pd.DataFrame:
    if not h5_paths:
        return pd.DataFrame(columns=["family_id", "station", "t_center", "dvv", "dvv_err", "cc_max", "n_det"])

    if n_workers <= 1:
        results = [
            measure_family(p, coda_window, fs, eps_max, n_eps, ref_start, ref_end, min_cc)
            for p in h5_paths
        ]
    else:
        results = Parallel(n_jobs=n_workers, prefer="processes")(
            delayed(measure_family)(
                p, coda_window, fs, eps_max, n_eps, ref_start, ref_end, min_cc
            )
            for p in h5_paths
        )
    rows = [r for batch in results for r in batch]
    return pd.DataFrame(rows)


def ref_window_bounds(t_start: datetime, ref_window_days: tuple[int, int]) -> tuple[datetime, datetime]:
    a, b = ref_window_days
    return t_start + timedelta(days=a), t_start + timedelta(days=b)
