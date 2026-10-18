#!/usr/bin/env python
"""Quiet-day waveform test, step 0: pull NW core-network candidates (all series) + all core data candidates (timing test)
from the production score tables; add the post-stage day class (wa_prod_post.day_class: quiet = no PNSN tremor <= 60 km
that UTC day, before 2009-08-06 outside the WECC ETS windows). Read-only on production outputs."""
import os, sys, glob
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_extension_wa/scripts")
import wa_prod_common as P
import wa_prod_post as PP
OD = f"{P.OUT}/quietwave"

def one(p):
    D = pd.read_parquet(p, columns=["date", "net", "series", "fid", "k", "i", "t", "net_cc", "net_mad", "n_sta", "n_sites",
                                    "hfr_rel", "reg", "role", "P_own_WAVE", "P_nbr3_WAVE", "tpl_insample"])
    D = D[(D.net == "core") & ~D.role.astype(str).str.startswith("sealed")]
    A = D[D.series == "data"][["date", "fid", "reg", "t", "i", "net_mad", "P_own_WAVE"]].copy()
    N = D[(D.reg == "NW")].copy()
    N = PP.day_class(N)
    return N, A

if __name__ == "__main__":
    fs = sorted(glob.glob(f"{P.OUT}/conf/P/prod/P_*.parquet"))
    with Pool(16) as pool:
        R = pool.map(one, fs)
    N = pd.concat([r[0] for r in R], ignore_index=True)
    A = pd.concat([r[1] for r in R], ignore_index=True)
    N.to_parquet(f"{OD}/nw_core_candidates.parquet", index=False)
    A.to_parquet(f"{OD}/core_data_all_regions.parquet", index=False)
    print(len(N), len(A))
    for cls in ("quiet", "tremor30"):
        h = N[N[cls].astype(bool)]
        for c in (0.75, 0.85, 0.9):
            print(cls, c, h[h.P_own_WAVE >= c].series.value_counts().to_dict())
