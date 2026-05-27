"""Stack LFE detections in time bins per (family, station).

Writes data/stacks/{family_id}.h5 with one group per station and one dataset
per bin (compressed gzip).

Usage:
    python scripts/06_stack_bins.py \\
        --config configs/ets_2010_vi.yaml \\
        --detections data/detections \\
        --waveforms data/waveforms \\
        --selected catalogs/selected_families.csv \\
        --out data/stacks --workers 32
"""

from __future__ import annotations

import argparse
from datetime import timedelta
from pathlib import Path

import pandas as pd

from tremorferometry.config import load_config
from tremorferometry.stacking import make_bin_edges, stack_all_parallel


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--detections", type=Path, required=True)
    ap.add_argument("--waveforms", type=Path, required=True)
    ap.add_argument("--selected", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--fs", type=float, default=40.0,
                    help="Target sampling rate (matches BHZ at 40 Hz for 2010 V.I.)")
    ap.add_argument("--pre-s", type=float, default=5.0,
                    help="Seconds before each detection time (origin time) to start the cut")
    ap.add_argument("--post-s", type=float, default=35.0,
                    help="Seconds after each detection time to end the cut")
    args = ap.parse_args()

    cfg = load_config(args.config)
    selected = pd.read_csv(args.selected)

    if args.detections.is_dir():
        # Template-matching layout: one parquet per family
        all_det = []
        for fid in selected["family_id"].astype(str):
            p = args.detections / f"{fid}.parquet"
            if p.exists():
                all_det.append(pd.read_parquet(p))
        if not all_det:
            raise SystemExit(f"no detection parquets under {args.detections}")
        detections = pd.concat(all_det, ignore_index=True)
    else:
        # Catalog-direct layout: one parquet with all families
        if not args.detections.exists():
            raise SystemExit(f"detections file not found: {args.detections}")
        detections = pd.read_parquet(args.detections)
        # restrict to selected families
        keep = set(selected["family_id"].astype(str))
        detections = detections[detections["family_id"].astype(str).isin(keep)].copy()
    detections["time"] = pd.to_datetime(detections["time"])

    ref_pre, _ = cfg.dvv.reference_window
    t0 = cfg.episode.t_start + timedelta(days=min(ref_pre, -90))
    t1 = cfg.episode.t_end + timedelta(days=30)
    edges = make_bin_edges(t0, t1, cfg.stack.bin_days)

    paths = stack_all_parallel(
        detections=detections,
        families=selected["family_id"].astype(str).tolist(),
        stations=list(cfg.stations.list),
        waveform_root=args.waveforms,
        bin_edges=edges,
        out_dir=args.out,
        n_workers=args.workers,
        fs=args.fs,
        pre_s=args.pre_s,
        post_s=args.post_s,
        bandpass=cfg.dvv.freq_band,
    )
    print(f"wrote {len(paths)} family HDF5 files to {args.out}")


if __name__ == "__main__":
    main()
