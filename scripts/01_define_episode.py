"""Define an ETS episode from a PNSN tremor CSV.

Inputs: a PNSN tremor catalog CSV (columns: time, lat, lon, depth).
Outputs: prints a refined (t_start, t_end, bbox) you can paste into the config.

Usage:
    python scripts/01_define_episode.py \\
        --tremor catalogs/pnsn_tremor_2010.csv \\
        --config configs/ets_2010_vi.yaml
"""

from __future__ import annotations

import argparse
from pathlib import Path

from tremorferometry.config import load_config
from tremorferometry.episode import detect_episode, load_pnsn_tremor


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tremor", type=Path, required=True, help="PNSN tremor CSV")
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--rate-threshold", type=int, default=20)
    ap.add_argument("--rate-window-hours", type=float, default=24.0)
    args = ap.parse_args()

    cfg = load_config(args.config)
    tremor = load_pnsn_tremor(args.tremor)
    ep = detect_episode(
        tremor,
        rate_window_hours=args.rate_window_hours,
        rate_threshold=args.rate_threshold,
    )
    print(f"event_id: {cfg.event_id}")
    print(f"episode.t_start: {ep.t_start.isoformat()}")
    print(f"episode.t_end:   {ep.t_end.isoformat()}")
    print(f"episode.bbox:    {list(ep.bbox)}")
    print(f"n_detections:    {ep.n_detections}")


if __name__ == "__main__":
    main()
