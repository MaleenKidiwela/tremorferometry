#!/usr/bin/env python
"""Merlin test T1 setup: a merged WA_SETUP for the all-templates matched-filter test.

families.csv  = the trial's 329 representatives (results/setup), val_windows + A2011 for families whose own
                array run scanned A2011 days (NW runs); A2011 truth = results/step4b/setup/truth.parquet
stations.csv  = the trial's 42 permanent stations (same sites)
truth.parquet = trial truth (V2008, V2010) + step-4b truth rows of A2011 (window column added)
scan_days.csv = union of both scan tables
  python scripts/merlin_t1_setup.py  -> results/merlin_tests/t1_setup/
"""
import os, sys, shutil
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

R = W.RES
OD = f"{R}/merlin_tests/t1_setup"
os.makedirs(OD, exist_ok=True)
F = pd.read_csv(f"{R}/setup/families.csv")
T0 = pd.read_parquet(f"{R}/setup/truth.parquet")
T4 = pd.read_parquet(f"{R}/step4b/setup/truth.parquet")
s0 = pd.read_csv(f"{R}/setup/scan_days.csv"); s4 = pd.read_csv(f"{R}/step4b/setup/scan_days.csv")
S = pd.concat([s0, s4]).drop_duplicates()
a11 = {d for d, _ in W.window_days("A2011")}
src11 = set(S[S.date.isin(a11)].det_source)
F["val_windows"] = [v + (",A2011" if s in src11 else "") for v, s in zip(F.val_windows, F.det_source)]
wmap = {}
for w in ("V2008", "V2010", "A2011"):
    for d, _ in W.window_days(w):
        wmap[d] = w
T0 = T0.copy(); T0["window"] = T0.date.map(wmap)
TA = T4[T4.window == "A2011"][["fid", "t", "date", "det_source", "half", "window"]]
T = pd.concat([T0, TA], ignore_index=True).drop_duplicates(["fid", "t"])
F.to_csv(f"{OD}/families.csv", index=False)
shutil.copy(f"{R}/setup/stations.csv", f"{OD}/stations.csv")
T.to_parquet(f"{OD}/truth.parquet", index=False)
S.to_csv(f"{OD}/scan_days.csv", index=False)
print("families", len(F), F.val_windows.value_counts().to_dict())
print("truth by window", T.groupby("window", dropna=False).size().to_dict())
