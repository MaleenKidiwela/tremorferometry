"""Stretching dv/v on per-bin stacks against a pre-ETS reference.

Driver lives in `tremorferometry.measure`; this script is a thin CLI.

Usage:
    python scripts/07_measure_dvv.py \\
        --config configs/ets_2010_vi.yaml \\
        --stacks data/stacks \\
        --out data/dvv/ets_2010_vi.parquet \\
        --workers 64
"""

from __future__ import annotations

import argparse
from pathlib import Path

from tremorferometry.config import load_config
from tremorferometry.io import write_dvv
from tremorferometry.measure import measure_many, ref_window_bounds


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--stacks", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=64)
    ap.add_argument("--fs", type=float, default=100.0)
    args = ap.parse_args()

    cfg = load_config(args.config)
    ref_start, ref_end = ref_window_bounds(cfg.episode.t_start, cfg.dvv.reference_window)

    h5_paths = sorted(args.stacks.glob("*.h5"))
    if not h5_paths:
        raise SystemExit(f"no stack files under {args.stacks}")

    df = measure_many(
        h5_paths=h5_paths,
        coda_window=cfg.dvv.coda_window,
        fs=args.fs,
        eps_max=cfg.dvv.stretch_range,
        n_eps=cfg.dvv.stretch_steps,
        ref_start=ref_start,
        ref_end=ref_end,
        min_cc=cfg.dvv.min_cc,
        n_workers=args.workers,
    )
    write_dvv(df, args.out)
    print(f"wrote {len(df)} dv/v rows -> {args.out}")


if __name__ == "__main__":
    main()
