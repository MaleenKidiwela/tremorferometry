"""Backfill continuous broadband waveforms for one station.

Generic version of backfill_pgc_2005-2013.py: pass --station and --network.
Skips day-files already on disk. ThreadPool, I/O-bound.

Typical use:
    python scripts/backfill_waveforms.py --network CN --station NLLB \
        --start 2005-01-01 --end 2026-12-31 --workers 24 --provider IRIS
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

from tremorferometry.waveforms import fetch_many


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--network", default="CN")
    ap.add_argument("--station", required=True)
    ap.add_argument("--start", default="2005-01-01")
    ap.add_argument("--end", default="2026-12-31")
    ap.add_argument("--root", default="data/waveforms")
    ap.add_argument("--provider", default="IRIS",
                    help="IRIS, NRCAN URL, NCEDC, etc.")
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--channel", default="BHZ,HHZ,EHZ")
    args = ap.parse_args()

    t0 = datetime.fromisoformat(args.start)
    t1 = datetime.fromisoformat(args.end)

    paths = fetch_many(
        root=Path(args.root),
        network=args.network,
        stations=[args.station],
        t_start=t0,
        t_end=t1,
        provider=args.provider,
        max_workers=args.workers,
        channel=args.channel,
    )
    print(f"wrote/verified {len(paths)} day-files to {args.root}")


if __name__ == "__main__":
    main()
