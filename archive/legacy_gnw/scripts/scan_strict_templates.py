"""Run matched-filter for the 16 STRICT new LFE templates at PGC.

The original 35 templates' MF detections were already saved at
`data/mf_pgc_2005-2026_cc08.csv`. The 16 STRICT templates were used to
build the canonical dv/v product but their per-detection times weren't
saved. We need them to build the long-window daily stacks for the
coda-window dv/v redo, so scan them here.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, "src")
from tremorferometry.matched_filter_fast import scan_many_days_multi  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--templates-npz", default="data/strict_new_templates.npz")
    p.add_argument("--wfdir", default="data/waveforms")
    p.add_argument("--station", default="PGC")
    p.add_argument("--threshold", type=float, default=0.70)
    p.add_argument("--min-gap-s", type=float, default=6.0)
    p.add_argument("--fmin", type=float, default=2.0)
    p.add_argument("--fmax", type=float, default=8.0)
    p.add_argument("--fs", type=float, default=40.0)
    p.add_argument("--workers", type=int, default=32)
    p.add_argument("--out", default="data/mf_pgc_strict.csv")
    return p.parse_args()


def main():
    args = parse_args()

    tdict = np.load(args.templates_npz, allow_pickle=True)
    templates = {k: np.asarray(tdict[k], dtype=np.float32) for k in tdict.keys()
                 if str(k).startswith("STRICT")}
    print(f"Loaded {len(templates)} STRICT templates: {list(templates)[:5]}...")

    # Discover days with PGC data
    wfdir = Path(args.wfdir) / f"CN.{args.station}"
    days = []
    for ypath in sorted(wfdir.glob("[0-9]*")):
        try:
            year = int(ypath.name)
        except ValueError:
            continue
        for mfile in sorted(ypath.glob("*.mseed")):
            try:
                jday = int(mfile.stem)
            except ValueError:
                continue
            d = datetime(year, 1, 1) + timedelta(days=jday - 1)
            days.append(d)
    print(f"PGC has {len(days)} day-files spanning {days[0].date()} - {days[-1].date()}")

    df = scan_many_days_multi(
        waveform_root=Path(args.wfdir),
        station=args.station,
        days=days,
        templates=templates,
        fs=args.fs,
        bandpass=(args.fmin, args.fmax),
        threshold=args.threshold,
        min_gap_s=args.min_gap_s,
        n_workers=args.workers,
        progress_every=200,
    )
    print(f"Got {len(df):,} detections at cc>={args.threshold}")
    df.to_csv(args.out, index=False)
    print(f"saved {args.out}")


if __name__ == "__main__":
    main()
