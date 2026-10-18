#!/usr/bin/env python
"""MENDO (2026-09-29) day lists for the blind Royer-chain discovery at Mendocino.

scan days  = the 572 days of lfe_south/results/aa_shelly/stage2_days_{build,test}.json (union over the
             Shelly2025 keys, 2018-01-06 .. 2024-11-30).  Used as-is: no further selection.
auto days  = the autocorrelation subset (Royer sec.3.1 autocorrelates a few episodes, sec.3.3 scans all
             dates).  Autocorrelating all 572 days would give ~60k detections -> a 60k x 60k similarity
             matrix and ~20k initial templates (see notes/2026-09-25_Notes.md MENDO entry), so a
             CONTIGUOUS block of 80 of the 572 listed days (NW used 75) is taken: the block with the most
             three-component S1 stations with raw data (BRIC / RVIT not counted), i.e. chosen by station
             coverage only, never by tremor or LFE activity.
-> results/ab_royer/tremor/scan_days_MENDO.csv, auto_days_MENDO.csv (column `date`)
"""
import json
import os

import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
SOUTH = f"{ROOT}/lfe_south"
OUT = f"{ROOT}/lfe_array_discovery/results/ab_royer/tremor"
NBLOCK = 80
S3 = ["KMR", "B049", "SPRL", "B933", "JEWT", "PRDS", "GHOP", "BLAS", "B934", "HALS", "B932", "KMPB",
      "HUNT"]
NET = {s: ("PB" if s[0] == "B" and s[1:].isdigit() else "NC" if s[0] == "K" else "BK") for s in S3}


def main():
    u = set()
    for h in ("build", "test"):
        for k, v in json.load(open(f"{SOUTH}/results/aa_shelly/stage2_days_{h}.json")).items():
            if k.split("|")[0] == "Shelly2025":
                u |= set(v)
    days = sorted(u)
    os.makedirs(OUT, exist_ok=True)
    pd.DataFrame(dict(date=days)).to_csv(f"{OUT}/scan_days_MENDO.csv", index=False)
    M = np.array([[os.path.exists(f"{SOUTH}/traces/{NET[s]}.{s}/{d[:4]}/"
                                  f"{pd.Timestamp(d).dayofyear:03d}.mseed") for s in S3] for d in days])
    n = M.sum(1).astype(float)
    cs = np.convolve(n, np.ones(NBLOCK), "valid") / NBLOCK
    i = int(np.argmax(cs))
    blk = days[i:i + NBLOCK]
    pd.DataFrame(dict(date=blk)).to_csv(f"{OUT}/auto_days_MENDO.csv", index=False)
    print(f"scan days {len(days)} ({days[0]}..{days[-1]}); auto block {blk[0]}..{blk[-1]} "
          f"({len(blk)} days, mean {cs[i]:.2f} / min {int(n[i:i + NBLOCK].min())} 3-C S1 stations)")


if __name__ == "__main__":
    main()
