"""Fetch the PNSN tremor catalog over a date range into a normalized CSV.

Output CSV columns: time, lat, lon, depth, magnitude, num_stas, duration, energy, id.
The first four columns satisfy `episode.load_pnsn_tremor` so this drops in to
`scripts/01_define_episode.py` directly.

Usage:
    python scripts/00b_fetch_pnsn_tremor.py \\
        --start 2024-01-01 --end 2024-07-01 \\
        --bbox 47.5 50.0 -125.5 -122.0 \\
        --out catalogs/pnsn_tremor_2024_h1.csv
"""

from __future__ import annotations

import argparse
import logging
from datetime import datetime
from pathlib import Path

from tremorferometry.pnsn import fetch_events, filter_bbox, write_csv

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--start", type=str, required=True, help="YYYY-MM-DD")
    ap.add_argument("--end", type=str, required=True, help="YYYY-MM-DD (exclusive)")
    ap.add_argument(
        "--bbox", type=float, nargs=4, default=None,
        metavar=("LAT_MIN", "LAT_MAX", "LON_MIN", "LON_MAX"),
        help="Optional spatial filter, lat_min lat_max lon_min lon_max",
    )
    ap.add_argument("--chunk-days", type=int, default=14)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    t0 = datetime.fromisoformat(args.start)
    t1 = datetime.fromisoformat(args.end)
    df = fetch_events(t0, t1, chunk_days=args.chunk_days)
    print(f"fetched {len(df)} events from PNSN tremor API")
    if args.bbox is not None:
        df = filter_bbox(df, tuple(args.bbox))
        print(f"after bbox filter: {len(df)}")
    write_csv(df, args.out)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
