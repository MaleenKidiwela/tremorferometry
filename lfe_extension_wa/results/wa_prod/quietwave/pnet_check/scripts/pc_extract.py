#!/usr/bin/env python
"""PNET check, step 0: core-network candidates of NW / SW / Arm B with both P_own_WAVE and P_own_WAVE_PNET, plus the v1
day class (wa_prod_post.day_class; same as the quiet-day waveform test).  Read-only on production outputs."""
import os, sys, glob
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_extension_wa/scripts")
import wa_prod_common as P
import wa_prod_post as PP
OD = f"{P.OUT}/quietwave/pnet_check"

def one(p):
    D = pd.read_parquet(p, columns=["date", "net", "series", "fid", "i", "t", "net_cc", "net_mad", "hfr_rel", "reg", "role",
                                    "P_own_WAVE", "P_own_WAVE_PNET"])
    D = D[(D.net == "core") & ~D.role.astype(str).str.startswith("sealed") & D.reg.isin(["NW", "SW", "ArmB"])]
    D = D[np.fmax(D.P_own_WAVE.fillna(0), D.P_own_WAVE_PNET.fillna(0)) >= 0.5].copy()
    return PP.day_class(D)

if __name__ == "__main__":
    fs = sorted(glob.glob(f"{P.OUT}/conf/P/prod/P_*.parquet"))
    with Pool(12) as pool:
        N = pd.concat(pool.map(one, fs), ignore_index=True)
    N = N.drop(columns=["net", "role"])
    N.to_parquet(f"{OD}/cand.parquet", index=False)
    cov = (N.date >= "2009-08-06") & (N.date <= "2026-05-27")
    N["cls"] = np.where(N.tremor30, "tremor30", np.where(N.quiet & cov, "quiet_cov", np.where(N.quiet & ~cov, "unknown", "other")))
    w = N.P_own_WAVE >= 0.75; p = N.P_own_WAVE_PNET >= 0.75
    N["sel"] = np.select([w & p, p & ~w, w & ~p], ["both", "pnet_only", "wave_only"], "none")
    print(N[N.sel != "none"].groupby(["reg", "cls", "series", "sel"]).size().unstack("sel").fillna(0).astype(int).to_string())
