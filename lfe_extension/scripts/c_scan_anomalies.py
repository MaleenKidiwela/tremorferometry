#!/usr/bin/env python
"""Track C — scan the whole processed record for impossible detections.

A normalised network cc cannot exceed 1 (Cauchy-Schwarz), so |net_cc| > 1 is proof that the window
energy in the denominator was wrong — the failure Track E found on 2026-09-17.  net_mad > 50 is the
same thing seen through the MAD normalisation.  Both are reported here, per month and per date.

  python scripts/c_scan_anomalies.py                 # table + results/trackC/anomaly_scan.csv
  python scripts/c_scan_anomalies.py --cc 1.0 --mad 50
"""
from __future__ import annotations
import os, sys, glob, argparse
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cc", type=float, default=1.0)
    ap.add_argument("--mad", type=float, default=50.0)
    ap.add_argument("--out-dir", default="")
    a = ap.parse_args()
    OUT = a.out_dir or C.TRC
    rows, tot, n_months = [], 0, 0
    for f in sorted(glob.glob(f"{OUT}/detections_*.parquet")):
        ym = os.path.basename(f).rsplit("_", 1)[1][:-8]
        d = pd.read_parquet(f, columns=["fam", "date", "t", "net_mad", "net_cc",
                                        "cat_core", "cat_ext"])
        n_months += 1
        tot += len(d)
        b = d[(d.net_cc.abs() > a.cc) | (d.net_mad > a.mad)]
        if len(b):
            rows.append(b.assign(ym=ym))
    print(f"[scan] {n_months} months, {tot:,} detections, thresholds |net_cc| > {a.cc} or "
          f"net_mad > {a.mad}")
    if not rows:
        print("[scan] CLEAN — no impossible detection anywhere in the record")
        pd.DataFrame(columns=["ym", "date", "n", "max_cc", "max_mad", "cat_core", "cat_ext"]).to_csv(
            f"{OUT}/anomaly_scan.csv", index=False)
        return
    B = pd.concat(rows, ignore_index=True)
    g = B.groupby(["ym", "date"]).agg(n=("t", "size"), max_cc=("net_cc", lambda x: float(x.abs().max())),
                                      max_mad=("net_mad", "max"), cat_core=("cat_core", "sum"),
                                      cat_ext=("cat_ext", "sum")).reset_index()
    pd.set_option("display.width", 250)
    print(g.to_string(index=False))
    print(f"\n[scan] {len(B)} anomalous rows on {g.date.nunique()} dates in {g.ym.nunique()} months; "
          f"{int(B.cat_core.sum())} are in cat_core, {int(B.cat_ext.sum())} in cat_ext")
    g.to_csv(f"{OUT}/anomaly_scan.csv", index=False)
    B.to_csv(f"{OUT}/anomaly_scan_rows.csv", index=False)
    print(f"[scan] -> {OUT}/anomaly_scan.csv, anomaly_scan_rows.csv")


if __name__ == "__main__":
    main()
