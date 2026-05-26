"""Run EQcorrscan template matching per family, write detections to parquet.

Uses fast-matched-filter (GPU) when the optional [gpu] extra is installed.

Usage:
    python scripts/05_match_lfe.py \\
        --config configs/ets_2010_vi.yaml \\
        --selected catalogs/selected_families.csv \\
        --templates data/templates \\
        --waveforms data/waveforms \\
        --out data/detections
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from tremorferometry.config import load_config
from tremorferometry.io import write_detections
from tremorferometry.matching import detect_family


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--selected", type=Path, required=True)
    ap.add_argument("--templates", type=Path, required=True)
    ap.add_argument("--waveforms", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--no-gpu", action="store_true")
    args = ap.parse_args()

    cfg = load_config(args.config)
    selected = pd.read_csv(args.selected)
    args.out.mkdir(parents=True, exist_ok=True)

    waveform_paths = sorted(args.waveforms.rglob("*.mseed"))
    if not waveform_paths:
        raise SystemExit(f"no continuous waveforms found under {args.waveforms}")

    for _, fam in selected.iterrows():
        family_id = str(fam["family_id"])
        template = args.templates / f"{family_id}.mseed"
        if not template.exists():
            print(f"skip {family_id}: no template")
            continue
        det = detect_family(
            template_path=template,
            waveform_paths=waveform_paths,
            cc_threshold=cfg.match.cc_threshold,
            mad_multiplier=cfg.match.mad_multiplier,
            use_gpu=not args.no_gpu,
        )
        out_path = args.out / f"{family_id}.parquet"
        write_detections(det, out_path)
        print(f"{family_id}: {len(det)} detections -> {out_path}")


if __name__ == "__main__":
    main()
