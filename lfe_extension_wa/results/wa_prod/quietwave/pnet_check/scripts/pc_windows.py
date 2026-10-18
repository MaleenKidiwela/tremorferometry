#!/usr/bin/env python
"""PNET check, step 1: rows to score + cut their station windows with the quiet-day test's exact replay
(quietwave/scripts/qw_windows.do_day, unchanged; only the output folder differs).

Rows (core network; sel = both / pnet_only / wave_only at 0.75 on P_own_WAVE and P_own_WAVE_PNET):
  NW quiet_cov (no PNSN tremor <= 60 km, 2009-08-06..2026-05-27), NW unknown (no catalogue), NW tremor30:
      every data / scrambled / reversed row selected by either model at 0.75
  SW, Arm B tremor30: every data row with either P >= 0.5 (the 0.5 rows feed the reference stacks), controls at 0.75
  random times: one per data row (same family, same day), NW tremor30: one per 4 data rows."""
import os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_extension_wa/results/wa_prod/quietwave/scripts")
import qw_windows as QW
OD = f"{QW.P.OUT}/quietwave/pnet_check"
QW.WK = f"{OD}/work"
CAT0, CAT1 = "2009-08-06", "2026-05-27"


def build_rows():
    N = pd.read_parquet(f"{OD}/cand.parquet")
    cov = (N.date >= CAT0) & (N.date <= CAT1)
    N["cls"] = np.where(N.tremor30, "tremor30", np.where(N.quiet & cov, "quiet_cov", np.where(N.quiet & ~cov, "unknown", "other")))
    w = N.P_own_WAVE.fillna(0); p = N.P_own_WAVE_PNET.fillna(0)
    N["sel"] = np.select([(w >= .75) & (p >= .75), (p >= .75) & (w < .75), (w >= .75) & (p < .75)], ["both", "pnet_only", "wave_only"], "none")
    m_nw = (N.reg == "NW") & N.cls.isin(["quiet_cov", "unknown", "tremor30"]) & (N.sel != "none")
    m_s = N.reg.isin(["SW", "ArmB"]) & (N.cls == "tremor30") & ((N.series == "data") | (N.sel != "none"))
    S = N[m_nw | m_s].copy()
    rng = np.random.default_rng(20261006)
    dd = S[S.series == "data"]
    dd = dd[~((dd.reg == "NW") & (dd.cls == "tremor30")) | (rng.random(len(dd)) < 0.25)]
    Rr = dd[["date", "fid", "reg", "cls", "sel", "P_own_WAVE", "P_own_WAVE_PNET", "tremor30", "quiet"]].copy()
    Rr["i"] = rng.integers(0, QW.NV, len(Rr)); Rr["series"] = "random"
    Rr["t"] = pd.to_datetime(Rr.date).values.astype("datetime64[ns]").astype(np.int64) / 1e9 + Rr.i / QW.FS
    S = pd.concat([S, Rr], ignore_index=True)
    S["row_id"] = np.arange(len(S))
    S.to_parquet(f"{OD}/rows.parquet", index=False)
    return S


if __name__ == "__main__":
    os.makedirs(QW.WK, exist_ok=True)
    S = build_rows()
    print(S.groupby(["reg", "cls", "series", "sel"]).size().to_string(), flush=True)
    jobs = [(ds, g[["row_id", "fid", "series", "i"]]) for ds, g in S.groupby("date")]
    jobs.sort(key=lambda x: -len(x[1]))
    print(len(jobs), "days", flush=True)
    t0 = time.time()
    with Pool(int(os.environ.get("QW_NPROC", "20")), initializer=QW.init) as pool:
        for k, (ds, res) in enumerate(pool.imap_unordered(QW.do_day, jobs)):
            if k % 100 == 0:
                print(f"{k}/{len(jobs)} {ds} {res} {time.time() - t0:.0f}s", flush=True)
    print("done", time.time() - t0, flush=True)
