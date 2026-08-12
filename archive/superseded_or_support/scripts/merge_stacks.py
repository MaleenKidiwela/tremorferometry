"""Merge several long_window_daily_*.npz stacks into one (concatenate the
per-(patch,date) rows, dedupe overlapping family-days, keep the better-sampled).
Avoids re-stacking families that are already stacked.

Usage:
  python scripts/merge_stacks.py --out data/long_window_daily_GNW_merged.npz \
      data/long_window_daily_GNW.npz data/long_window_daily_GNWcircle.npz \
      data/long_window_daily_GNWgap.npz
"""
from __future__ import annotations

import argparse

import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True)
    p.add_argument("inputs", nargs="+")
    args = p.parse_args()

    stacks, patches, dates, n_det = [], [], [], []
    t_ref = fs_ref = None
    for f in args.inputs:
        d = np.load(f, allow_pickle=True)
        if t_ref is None:
            t_ref, fs_ref = d["t"], float(d["fs"])
        else:
            assert d["stacks"].shape[1] == len(t_ref), f"{f}: sample-count mismatch"
            assert abs(float(d["fs"]) - fs_ref) < 1e-6, f"{f}: fs mismatch"
        stacks.append(d["stacks"]); patches.append(d["patches"])
        dates.append(d["dates"]); n_det.append(d["n_det"])
        print(f"  {f}: {d['stacks'].shape[0]:,} rows, {len(set(d['patches']))} patches")

    stacks = np.concatenate(stacks); patches = np.concatenate(patches)
    dates = np.concatenate(dates); n_det = np.concatenate(n_det)
    print(f"concatenated: {len(stacks):,} rows, {len(set(patches))} distinct patches")

    # dedupe (patch,date): keep the row with the most detections
    order = np.argsort(-n_det, kind="stable")            # high n_det first
    key = np.char.add(np.char.add(patches[order], "@"), dates[order])
    _, keep_in_order = np.unique(key, return_index=True)
    keep = np.sort(order[keep_in_order])
    print(f"after dedupe (patch,date): {len(keep):,} rows "
          f"({len(stacks)-len(keep):,} duplicate family-days dropped)")

    np.savez(args.out, stacks=stacks[keep], patches=patches[keep],
             dates=dates[keep], n_det=n_det[keep], t=t_ref, fs=fs_ref)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
