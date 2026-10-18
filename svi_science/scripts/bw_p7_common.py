#!/usr/bin/env python
"""bw_p7_common -- shared loaders for P7 (rules: results/backward/p7/P7.md)."""
import os, sys, glob, json
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

ROOT = "/home/jovyan/tremorferometry"
BW = f"{ROOT}/svi_science/results/backward"
P6 = f"{BW}/p6"; P7 = f"{BW}/p7"
CAT = f"{ROOT}/data/catalogue/svi_catalogue_v1_4_backward.parquet"
PROD = [f"p5f_{y}" for y in range(2003, 2014)] + ["p5e_2003_kelb", "p5e_2004", "p5e_2005", "p5e_2006"]
ERA_DET = ("2003-02-01", "2006-06-30")


def run_dir(run):
    if run.startswith("lad_"): return f"{BW}/p4b/mf/feat/{run}"
    if run.startswith("p4e0_"): return f"{BW}/p4/mf/feat/{run}"
    return f"{BW}/mf/feat/{run}"


def _fd(args):
    run, p = args
    if os.path.getsize(p) == 0: return None
    f = pd.read_parquet(p, columns=["series", "fam", "date", "observable"])
    f["observable"] = f.observable.astype(bool); f["run"] = run
    return f


def load_fd(runs):
    jobs = [(r, p) for r in runs for p in sorted(glob.glob(f"{run_dir(r)}/*_fd.parquet"))]
    with ProcessPoolExecutor(8) as ex:
        out = [x for x in ex.map(_fd, jobs, chunksize=16) if x is not None]
    return pd.concat(out, ignore_index=True)


def day_flags():
    import bw_p4_common as P4
    Q = pd.read_csv(f"{BW}/p4/envflag/tables/quiet_flag_env.csv", dtype={"date": str})
    Q["ets"] = P4.bostock_near(Q.date, pad=0); Q["between"] = ~P4.bostock_near(Q.date, pad=3)
    Q["env_quiet"] = Q.quiet_env3_CN4.fillna(False).astype(bool)
    Q["pnsn_quiet_b"] = (Q.pnsn_quiet == True) & (Q.date >= "2010-01-10")       # noqa: E712
    return Q.set_index("date")[["ets", "between", "env_quiet", "pnsn_quiet_b", "n_tremor_day"]]


def k_chain(rows="all"):
    T = pd.read_csv(f"{P7}/tables/k_chain.csv")
    T = T[T.rows == rows].set_index("net_set")
    return T


BOUND = [("2005-07-01", "B1"), ("2005-11-02", "B2"), ("2007-11-13", "B3"), ("2009-08-27", "S15")]


def set_of_date(d):
    s = "CN4"
    for b, n in BOUND:
        if d >= b: s = n
    return s


def daily(A, FD):
    """per (run, date) over EVERY day of the run: ΣP data/control (all rows, multi-site n_z2 >= 2 rows), observable family-days,
    burst flag, detector, net_set.  A needs run, date, series, P, multi_site_z2, burst_day."""
    A = A[["run", "date", "series", "P", "multi_site_z2", "burst_day"]].copy()
    A["Pd"] = np.where(A.series == "data", A.P, 0.0); A["Pc"] = np.where(A.series == "c1", A.P, 0.0)
    ms = A.multi_site_z2.astype(bool).to_numpy()
    A["Pd_ms"] = A.Pd * ms; A["Pc_ms"] = A.Pc * ms
    g = A.groupby(["run", "date"]).agg(spd=("Pd", "sum"), spc=("Pc", "sum"), spd_ms=("Pd_ms", "sum"), spc_ms=("Pc_ms", "sum"),
                                       burst=("burst_day", "max")).reset_index()
    m = FD[FD.observable].groupby(["run", "date", "series"]).size().unstack(fill_value=0).reset_index()
    m = m.rename(columns={"data": "m_data", "c1": "m_c1"})
    for c in ("m_data", "m_c1"):
        if c not in m: m[c] = 0
    D = m.merge(g, on=["run", "date"], how="left")
    for c in ("spd", "spc", "spd_ms", "spc_ms"): D[c] = D[c].fillna(0.0)
    D["burst"] = D.burst.fillna(False).astype(bool)
    D["detector"] = np.where(D.run.str.startswith("p5e_"), "era", "frozen")
    D["net_set"] = np.where(D.detector == "era", "ERA", D.date.map(set_of_date))
    D["primary"] = np.where(D.detector == "era", (D.date >= ERA_DET[0]) & (D.date <= ERA_DET[1]),
                            ~((D.date >= ERA_DET[0]) & (D.date <= ERA_DET[1])))
    return D
