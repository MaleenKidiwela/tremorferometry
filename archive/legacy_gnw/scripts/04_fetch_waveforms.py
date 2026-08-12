"""Fetch continuous FDSN waveforms for the ETS window +/- reference padding.

Pulls one MSEED per station per UTC day into data/waveforms/{net}.{sta}/{year}/{jday}.mseed.
Parallel by default (ThreadPool, I/O-bound).

Usage:
    python scripts/04_fetch_waveforms.py --config configs/ets_2010_vi.yaml --workers 16
"""

from __future__ import annotations

import argparse
import logging
from datetime import timedelta
from pathlib import Path

from tremorferometry.config import load_config
from tremorferometry.waveforms import fetch_many

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--root", type=Path, default=Path("data/waveforms"))
    ap.add_argument("--provider", type=str, default="IRIS")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--channel", type=str, default="HH?,BH?,EH?")
    args = ap.parse_args()

    cfg = load_config(args.config)
    ref_pre, ref_post = cfg.dvv.reference_window  # days, both negative-or-zero typically
    t0 = cfg.episode.t_start + timedelta(days=min(ref_pre, ref_post, -90))
    t1 = cfg.episode.t_end + timedelta(days=30)

    paths = fetch_many(
        root=args.root,
        network=cfg.stations.network,
        stations=cfg.stations.list,
        t_start=t0,
        t_end=t1,
        provider=args.provider,
        max_workers=args.workers,
        channel=args.channel,
    )
    print(f"wrote {len(paths)} day-files to {args.root}")


if __name__ == "__main__":
    main()
