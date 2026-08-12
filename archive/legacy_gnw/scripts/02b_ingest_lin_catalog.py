"""Ingest the Lin (2023) Cascadia LFE catalog into family + detection tables.

Lin's catalog is per-event with no family labels; this script clusters detections
by 0.1-deg grid binning (configurable) within the configured episode window +
reference window, and writes:

    catalogs/lin_families_<event_id>.csv         (FAMILY_COLUMNS schema)
    data/detections_lin_<event_id>.parquet       (family_id, time)

Download the raw file first:

    mkdir -p data/raw_lfe
    curl -fsSL -o data/raw_lfe/lin2023_lfe.csv \\
        'https://zenodo.org/records/10016020/files/EQloc_001_0.1_3_S.csv?download=1'

Usage:
    python scripts/02b_ingest_lin_catalog.py --config configs/ets_2010_vi.yaml \\
        --raw data/raw_lfe/lin2023_lfe.csv
"""

from __future__ import annotations

import argparse
from datetime import timedelta
from pathlib import Path

import pandas as pd

from tremorferometry.config import load_config
from tremorferometry.lin_catalog import ingest_lin_catalog


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--raw", type=Path, required=True)
    ap.add_argument("--cell-deg", type=float, default=0.1)
    ap.add_argument("--min-n", type=int, default=50,
                    help="minimum detections per cell to keep it as a family")
    ap.add_argument("--out-families", type=Path, default=None)
    ap.add_argument("--out-detections", type=Path, default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    # ingest covers ref window + ETS + post; we want detections for the full
    # span the dv/v measurement will use.
    ref_pre_d, _ = cfg.dvv.reference_window
    t_start = cfg.episode.t_start + timedelta(days=min(ref_pre_d, -90))
    t_end = cfg.episode.t_end + timedelta(days=30)

    families, detections = ingest_lin_catalog(
        csv_path=args.raw,
        t_start=t_start,
        t_end=t_end,
        bbox=cfg.episode.bbox,
        cell_deg=args.cell_deg,
        min_n=args.min_n,
    )

    out_fam = args.out_families or Path(f"catalogs/lin_families_{cfg.event_id}.csv")
    out_det = args.out_detections or Path(f"data/detections_lin_{cfg.event_id}.parquet")
    out_fam.parent.mkdir(parents=True, exist_ok=True)
    out_det.parent.mkdir(parents=True, exist_ok=True)
    families.to_csv(out_fam, index=False)
    detections.to_parquet(out_det, index=False)
    print(f"families: {len(families)} -> {out_fam}")
    print(f"detections: {len(detections):,} -> {out_det}")
    print(f"detections per family (median): {families['n_detections'].median():.0f}")
    print(f"top 5 families by detections:")
    print(families.head(5).to_string(index=False))


if __name__ == "__main__":
    main()
