"""Stretching dv/v on per-bin stacks against a pre-ETS reference.

For each (family, station): build the reference as the mean of bin-stacks whose
center falls in the reference window; then for every bin measure dv/v vs that
reference. Embarrassingly parallel across families x stations x bins.

Usage:
    python scripts/07_measure_dvv.py \\
        --config configs/ets_2010_vi.yaml \\
        --stacks data/stacks \\
        --out data/dvv/<event_id>.parquet \\
        --workers 64
"""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from tremorferometry.config import load_config
from tremorferometry.dvv import stretch_dvv
from tremorferometry.io import write_dvv


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


def _process_family(
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
            for idx, t_center, arr in bins:
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
                    "n_det": int(arr.size),  # placeholder; refine to true count later
                })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--stacks", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=64)
    ap.add_argument("--fs", type=float, default=100.0)
    args = ap.parse_args()

    cfg = load_config(args.config)
    ref_pre_d, ref_post_d = cfg.dvv.reference_window
    ref_start = cfg.episode.t_start + timedelta(days=ref_pre_d)
    ref_end = cfg.episode.t_start + timedelta(days=ref_post_d)

    h5_paths = sorted(args.stacks.glob("*.h5"))
    if not h5_paths:
        raise SystemExit(f"no stack files under {args.stacks}")

    results = Parallel(n_jobs=args.workers, prefer="processes")(
        delayed(_process_family)(
            p,
            cfg.dvv.coda_window,
            args.fs,
            cfg.dvv.stretch_range,
            cfg.dvv.stretch_steps,
            ref_start,
            ref_end,
            cfg.dvv.min_cc,
        )
        for p in h5_paths
    )
    rows = [r for batch in results for r in batch]
    df = pd.DataFrame(rows)
    write_dvv(df, args.out)
    print(f"wrote {len(df)} dv/v rows -> {args.out}")


if __name__ == "__main__":
    main()
