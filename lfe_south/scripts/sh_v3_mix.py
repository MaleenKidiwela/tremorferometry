#!/usr/bin/env python
"""SH_V3: Text S1 "for six families with poor signal-to-noise ... an additional iteration".
The 6 families with the lowest build-day stack quality (median over stations of cc50 = mean
leave-one-out cc of the kept 50 events, tpl_it0) take their it1 detections; the other 21 keep it0.
Writes mf_mix_test.parquet + mix_families.csv.  Uses build-day information only."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sh_v3_common import *
T = pd.read_csv(f"{V3}/tpl_it0.csv"); T["fam"] = T.fam.astype(str)
q = T.groupby("fam").cc50.median().sort_values()
low = list(q.index[:6])
pd.DataFrame(dict(fam=q.index, cc50_median=q.values, it=[1 if f in low else 0 for f in q.index])).to_csv(f"{V3}/mix_families.csv", index=False)
A = pd.read_parquet(f"{V3}/mf_it0_test.parquet"); A["fam"] = A.fam.astype(str)
B = pd.read_parquet(f"{V3}/mf_it1_test.parquet"); B["fam"] = B.fam.astype(str)
M = pd.concat([A[~A.fam.isin(low)], B[B.fam.isin(low)]], ignore_index=True).sort_values("t")
M.to_parquet(f"{V3}/mf_mix_test.parquet")
log(f"[mix] low-SNR families (it1): {low}; {len(M)} detections")
