"""Subset the Bostock LFE catalog to families inside the ETS slipping patch.

Usage:
    python scripts/03_select_families.py \\
        --config configs/ets_2010_vi.yaml \\
        --catalog catalogs/bostock_lfe_families.csv \\
        --out catalogs/selected_families.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

from tremorferometry.catalog import filter_by_bbox, load_families
from tremorferometry.config import load_config


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--catalog", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    fams = load_families(args.catalog)
    selected = filter_by_bbox(fams, cfg.episode.bbox)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.out, index=False)
    print(f"selected {len(selected)}/{len(fams)} families -> {args.out}")


if __name__ == "__main__":
    main()
