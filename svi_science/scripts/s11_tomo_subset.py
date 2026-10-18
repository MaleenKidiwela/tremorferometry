#!/usr/bin/env python
"""S11 tomo (Amendment 8): WA borehole-only subset (PB boreholes within 35 km epicentral) of the pick tables.
  python s11_tomo_subset.py   -> tables/picks_{BV3,CV3,B,C}bh.csv"""
import os, sys
import numpy as np, pandas as pd
ROOT = "/home/jovyan/tremorferometry"
T = f"{ROOT}/svi_science/results/s11_depth/tomo/tables"
sys.path.insert(0, f"{ROOT}/lfe_array_discovery/scripts")
cwd = os.getcwd(); os.chdir(f"{ROOT}/lfe_array_discovery")
import a4f_common as AC  # noqa: E402
os.chdir(cwd)
for st in ("BV3", "CV3", "B", "C"):
    P = pd.read_csv(f"{T}/picks_{st}.csv", low_memory=False)
    P["net"] = P.sta.map(AC.NET)
    P["slat"] = P.sta.map(lambda s: AC.COORD[s][0]); P["slon"] = P.sta.map(lambda s: AC.COORD[s][1])
    P["epi_km"] = np.hypot((P.slon - P.lon) * 111.19 * np.cos(np.radians(P.lat)), (P.slat - P.lat) * 111.19)
    Q = P[(P.net == "PB") & (P.epi_km <= 35)]
    ok = Q[(Q.rP >= 1.2) & (Q.rS >= 1.2)]
    Q.to_csv(f"{T}/picks_{st}bh.csv", index=False)
    print(st, "boreholes:", sorted(Q.sta.unique()), "rows ok", len(ok), "families with >=3 ok rows",
          int((ok.groupby("fam").size() >= 3).sum()), ">=4:", int((ok.groupby("fam").size() >= 4).sum()))
