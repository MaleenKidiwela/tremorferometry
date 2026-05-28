"""Matched-filter scan of newly-backfilled PGC day-files only.

After the 2005-2013 backfill we have ~1,800 new day-files that have not
been scanned. Diff the on-disk day-list against the existing detections
in mf_pgc_all51_cc08.csv, run the multi-template fast matched filter
across all 51 templates on just the new days, filter to cc in [0.8, 1.1],
and append.
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
    p.add_argument("--wfdir", default="data/waveforms")
    p.add_argument("--station", default="PGC")
    p.add_argument("--orig-npz", default="data/patch_templates.npz")
    p.add_argument("--strict-npz", default="data/strict_new_templates.npz")
    p.add_argument("--existing-csv", default="data/mf_pgc_all51_cc08.csv")
    p.add_argument("--out-csv", default="data/mf_pgc_new_days.csv")
    p.add_argument("--combined-csv", default="data/mf_pgc_all51_cc08.csv")
    p.add_argument("--threshold", type=float, default=0.7)
    p.add_argument("--cc-keep-lo", type=float, default=0.8)
    p.add_argument("--cc-keep-hi", type=float, default=1.1)
    p.add_argument("--workers", type=int, default=36)
    return p.parse_args()


def main():
    args = parse_args()

    print("[1/5] Loading templates...")
    orig = np.load(args.orig_npz, allow_pickle=True)
    strict = np.load(args.strict_npz, allow_pickle=True)
    templates: dict[str, np.ndarray] = {}
    # patch_templates.npz contains BOTH PGC_* and LZB_* templates (35 each, for
    # two-station discovery). For a PGC scan, only the PGC_* templates are
    # physically meaningful -- LZB templates against PGC waveforms produce
    # spurious detections we'd have to filter out downstream.
    for k in orig.keys():
        if not str(k).startswith(f"{args.station}_"):
            continue
        templates[k] = np.asarray(orig[k], dtype=np.float32)
    for k in strict.keys():
        if str(k).startswith("STRICT"):
            templates[k] = np.asarray(strict[k], dtype=np.float32)
    print(f"  {len(templates)} templates total ({sum(1 for k in templates if not k.startswith('STRICT'))} {args.station} + "
          f"{sum(1 for k in templates if k.startswith('STRICT'))} strict)")

    print("[2/5] Enumerating on-disk day-files...")
    wfdir = Path(args.wfdir) / f"CN.{args.station}"
    disk_days: set[datetime] = set()
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
            disk_days.add(datetime(year, 1, 1) + timedelta(days=jday - 1))
    print(f"  {len(disk_days)} day-files on disk")

    print("[3/5] Finding days already scanned...")
    exist = pd.read_csv(args.existing_csv, usecols=["time"])
    exist["time"] = pd.to_datetime(exist["time"], format="mixed")
    scanned_days = set(exist["time"].dt.normalize().dt.to_pydatetime())
    print(f"  {len(scanned_days)} scanned days in existing CSV")

    new_days = sorted(disk_days - scanned_days)
    print(f"  {len(new_days)} NEW days to scan")
    if not new_days:
        print("  nothing to do.")
        return

    print(f"[4/5] Scanning {len(new_days)} days x {len(templates)} templates "
          f"(workers={args.workers})...")
    df = scan_many_days_multi(
        waveform_root=Path(args.wfdir),
        station=args.station,
        days=new_days,
        templates=templates,
        fs=40.0,
        bandpass=(2.0, 8.0),
        threshold=args.threshold,
        min_gap_s=6.0,
        n_workers=args.workers,
        progress_every=200,
    )
    print(f"  {len(df):,} raw detections at cc>={args.threshold}")
    df.to_csv(args.out_csv, index=False)
    print(f"  saved raw {args.out_csv}")

    print("[5/5] Filter to cc in [0.8, 1.1] and append to combined CSV...")
    keep = df[(df["cc"] >= args.cc_keep_lo) & (df["cc"] <= args.cc_keep_hi)]
    print(f"  keep: {len(keep):,} at cc in [{args.cc_keep_lo}, {args.cc_keep_hi}]")

    existing = pd.read_csv(args.combined_csv)
    combined = pd.concat([existing, keep], ignore_index=True)
    combined.to_csv(args.combined_csv, index=False)
    print(f"  combined: {len(combined):,} rows in {args.combined_csv}")
    print("done")


if __name__ == "__main__":
    main()
