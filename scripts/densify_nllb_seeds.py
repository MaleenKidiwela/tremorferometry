"""Matched-filter densification of NLLB high-quality seed templates.

Inputs:
  - data/nllb_pnsn_families.npz : seed templates (NLLB-side, 2 s @ 40 Hz)
  - data/nllb_pnsn_families.summary.csv : with `snr` column
  - data/waveforms/CN.NLLB/* : continuous NLLB day-files

For seeds with template SNR >= --min-snr, run the batched-fast matched filter
across every NLLB day-file. Save raw detections at the lower CC threshold
(default 0.7); downstream filtering to cc>=0.8 happens at the dv/v step.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
from tremorferometry.matched_filter_fast import scan_many_days_multi  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--templates-npz", default="data/nllb_pnsn_families.npz")
    p.add_argument("--summary-csv", default="data/nllb_pnsn_families.summary.csv")
    p.add_argument("--min-snr", type=float, default=10.0)
    p.add_argument("--wfdir", default="data/waveforms")
    p.add_argument("--network", default="CN")
    p.add_argument("--station", default="NLLB")
    p.add_argument("--fs", type=float, default=40.0)
    p.add_argument("--fmin", type=float, default=2.0)
    p.add_argument("--fmax", type=float, default=8.0)
    p.add_argument("--threshold", type=float, default=0.7)
    p.add_argument("--min-gap-s", type=float, default=6.0)
    p.add_argument("--workers", type=int, default=32)
    p.add_argument("--out", default="data/mf_nllb_seeds_top.csv")
    return p.parse_args()


def main():
    args = parse_args()

    print(f"[1/3] Loading NLLB seed templates and filtering to SNR>={args.min_snr}...")
    d = np.load(args.templates_npz, allow_pickle=True)
    s = pd.read_csv(args.summary_csv).set_index("family_id")
    s_hi = s[s["snr"] >= args.min_snr]
    print(f"  {len(s_hi)} seed templates pass SNR threshold "
          f"(out of {len(s)} total)")
    if len(s_hi) == 0:
        raise SystemExit("no templates")

    templates = {}
    for fid in s_hi.index:
        if fid in d.files:
            templates[fid] = np.asarray(d[fid], dtype=np.float32)
    print(f"  loaded {len(templates)} templates")
    for k in list(templates)[:5]:
        print(f"    sample: {k}, shape {templates[k].shape}, "
              f"snr {s_hi.loc[k, 'snr']:.1f}")

    print(f"[2/3] Enumerating {args.station} day-files...")
    wfdir = Path(args.wfdir) / f"{args.network}.{args.station}"
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
            days.append(datetime(year, 1, 1) + timedelta(days=jday - 1))
    print(f"  {len(days)} {args.station} day-files spanning "
          f"{days[0].date()} - {days[-1].date()}")

    print(f"[3/3] Matched filter across all days "
          f"(threshold={args.threshold}, workers={args.workers})...")
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
    print(f"  raw detections at cc>={args.threshold}: {len(df):,}")
    df.to_csv(args.out, index=False)
    print(f"  saved {args.out}")

    # Per-template count summary
    counts = df["template"].value_counts()
    print()
    print(f"detections per template at cc>={args.threshold}:")
    print(f"  min={int(counts.min())}, median={int(counts.median())}, "
          f"max={int(counts.max())}, mean={int(counts.mean())}")
    for thr in [100, 1000, 10000, 100000]:
        n = (counts >= thr).sum()
        print(f"  >= {thr:,} detections: {n} families")


if __name__ == "__main__":
    main()
