"""Parse the Bostock LFE supplementary material into our normalized layout.

Expected raw inputs (place under data/raw_bostock/ before running):
  * families.txt           — family_id, lat, lon, depth_km, n_detections (whitespace or comma)
  * templates/{fam}.mseed  — multi-station template waveform per family

This script writes:
  * catalogs/bostock_lfe_families.csv
  * data/templates/{family_id}.mseed  (just copies from data/raw_bostock/templates/)

The exact format of Bostock's supplement varies; if your copy looks different,
edit `parse_families` below.

Usage:
    python scripts/02_ingest_bostock.py \\
        --raw-dir data/raw_bostock \\
        --out-catalog catalogs/bostock_lfe_families.csv \\
        --out-templates data/templates
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd

from tremorferometry.catalog import FAMILY_COLUMNS


def parse_families(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(path)
    except Exception:
        df = pd.read_csv(path, sep=r"\s+", header=None)
        df.columns = ["family_id", "lat", "lon", "depth_km", "n_detections"][: df.shape[1]]
    missing = set(FAMILY_COLUMNS) - set(df.columns)
    if missing:
        raise SystemExit(
            f"families file is missing columns {missing}; expected {FAMILY_COLUMNS}. "
            "Edit parse_families() to match your local Bostock supplement."
        )
    return df[list(FAMILY_COLUMNS)].copy()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw-dir", type=Path, required=True)
    ap.add_argument("--out-catalog", type=Path, required=True)
    ap.add_argument("--out-templates", type=Path, required=True)
    args = ap.parse_args()

    fams = parse_families(args.raw_dir / "families.txt")
    args.out_catalog.parent.mkdir(parents=True, exist_ok=True)
    fams.to_csv(args.out_catalog, index=False)
    print(f"wrote {len(fams)} families -> {args.out_catalog}")

    args.out_templates.mkdir(parents=True, exist_ok=True)
    raw_templates = args.raw_dir / "templates"
    for src in raw_templates.glob("*.mseed"):
        shutil.copy2(src, args.out_templates / src.name)
    print(f"copied {len(list(raw_templates.glob('*.mseed')))} templates -> {args.out_templates}")


if __name__ == "__main__":
    main()
