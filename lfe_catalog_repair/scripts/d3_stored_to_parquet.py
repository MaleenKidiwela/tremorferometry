#!/usr/bin/env python
"""Convert the STORED densify detections (data/mf_*.csv, Z-only cc>=0.80, 6 s gap; read-only) into
the same parquet schema d3_lin.py / d3_network.py consume, so the ground-truth test can be run on
the stored catalogue as the "before" leg."""
import argparse
import glob
from pathlib import Path

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--glob", required=True)
ap.add_argument("--start", required=True)
ap.add_argument("--end", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
d0, d1 = pd.Timestamp(a.start), pd.Timestamp(a.end)
parts = []
for f in sorted(glob.glob(a.glob)):
    for ch in pd.read_csv(f, usecols=["template", "time", "cc"], chunksize=4_000_000,
                          parse_dates=["time"]):
        ch = ch[(ch.time >= d0) & (ch.time < d1)]
        if len(ch):
            parts.append(ch)
d = pd.concat(parts, ignore_index=True)
d["t"] = d["time"].astype("int64") / 1e9
d["cc_stack"] = d["cc"].astype("float32")
d["old"] = True                       # every stored row already satisfies the fixed cc >= 0.80 rule
d["snr"] = float("nan")
out = d[["template", "t", "cc_stack", "old", "snr"]]
Path(a.out).parent.mkdir(parents=True, exist_ok=True)
out.to_parquet(a.out, index=False)
print(f"{len(out):,} stored detections {a.start}..{a.end}, {out.template.nunique()} templates -> {a.out}")
