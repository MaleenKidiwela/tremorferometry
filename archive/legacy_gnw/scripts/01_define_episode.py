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
    ap.add_argument("--select", choices=("largest", "first", "latest"), default="largest")
    ap.add_argument("--list-all", action="store_true",
                    help="Print every episode sorted by n_detections instead of picking one.")
    args = ap.parse_args()

    cfg = load_config(args.config)
    tremor = load_pnsn_tremor(args.tremor)

    if args.list_all:
        from tremorferometry.episode import list_episodes
        eps = list_episodes(
            tremor,
            rate_window_hours=args.rate_window_hours,
            rate_threshold=args.rate_threshold,
        )
        print(f"found {len(eps)} episodes:")
        for i, ep in enumerate(eps):
            print(f"  [{i:2d}] {ep.t_start.date()} -> {ep.t_end.date()}  "
                  f"n={ep.n_detections:5d}  "
                  f"bbox=lat[{ep.bbox[0]:.2f},{ep.bbox[1]:.2f}] lon[{ep.bbox[2]:.2f},{ep.bbox[3]:.2f}]")
        return

    ep = detect_episode(
        tremor,
        rate_window_hours=args.rate_window_hours,
        rate_threshold=args.rate_threshold,
        select=args.select,
    )
    print(f"event_id: {cfg.event_id}")
    print(f"episode.t_start: {ep.t_start.isoformat()}")
    print(f"episode.t_end:   {ep.t_end.isoformat()}")
    print(f"episode.bbox:    {list(ep.bbox)}")
    print(f"n_detections:    {ep.n_detections}")


if __name__ == "__main__":
    main()
