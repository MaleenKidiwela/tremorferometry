"""Backfill PGC continuous waveforms for 2005-01-01 .. 2014-01-01.

Pre-2014 we only fetched ETS-active months. Now fetch every day so the
matched-filter + daily-dv/v products are continuous back to 2005.

Skips day-files we already have. About 9 years x 365 = ~3000 new files,
~17 GB. ThreadPool, I/O bound.
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

from tremorferometry.waveforms import fetch_many


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", default="2005-01-01")
    ap.add_argument("--end", default="2014-01-01")
    ap.add_argument("--root", default="data/waveforms")
    ap.add_argument("--provider", default="IRIS")
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--channel", default="BHZ,HHZ,EHZ")
    args = ap.parse_args()

    t0 = datetime.fromisoformat(args.start)
    t1 = datetime.fromisoformat(args.end)

    paths = fetch_many(
        root=Path(args.root),
        network="CN",
        stations=["PGC"],
        t_start=t0,
        t_end=t1,
        provider=args.provider,
        max_workers=args.workers,
        channel=args.channel,
    )
    print(f"wrote/verified {len(paths)} day-files to {args.root}")


if __name__ == "__main__":
    main()
