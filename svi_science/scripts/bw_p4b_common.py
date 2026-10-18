#!/usr/bin/env python
"""bw_p4b_common -- burst-robust control-rate matching for P4-2b (rules: results/backward/p4b/P4B.md).
Control candidates = c1 rows with net_mad >= T (no HF veto); control family-days = observable c1 rows of _fd.
MED  = median over days of the per-day control rate;  TRIM = pooled rate after dropping given days;
POOLED = pooled rate over all days (the P4-2 statistic).  Matched T = smallest T on the 0.01 grid with statistic <= target."""
import os, glob, json
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

BW = "/home/jovyan/tremorferometry/svi_science/results/backward"
TG = np.round(np.arange(8.0, 14.0 + 1e-9, 0.01), 2)
ZC = None


def _day(args):
    d, ds, zcols = args
    p, q = f"{d}/{ds}.parquet", f"{d}/{ds}_fd.parquet"
    c = np.zeros(0); Z = None; data = np.zeros(0)
    if os.path.exists(p) and os.path.getsize(p) > 0:
        x = pd.read_parquet(p)
        if len(x):
            m = (x.series == "c1").to_numpy()
            c = x.loc[m, "net_mad"].to_numpy(float)
            data = x.loc[~m & (x.series == "data").to_numpy(), "net_mad"].to_numpy(float)
            if zcols:
                Z = x.loc[m, [k for k in zcols + ["dom_top"] if k in x.columns]].reset_index(drop=True)
    mc = md = 0; sta = None
    if os.path.exists(q) and os.path.getsize(q) > 0:
        f = pd.read_parquet(q, columns=["series", "observable"])
        ob = f.observable.astype(bool)
        mc = int(((f.series == "c1") & ob).sum()); md = int(((f.series == "data") & ob).sum())
    js = f"{d}/{ds}.json"
    if os.path.exists(js):
        sta = json.load(open(js)).get("stations")
    return ds, np.sort(c), mc, Z, sta, np.sort(data), md


def load_ctl(featdir, dates=None, zcols=None, workers=6):
    """-> DataFrame per day: date, m_c1, m_data, stations; dict date -> sorted c1 net_mad; (opt) per-day c1 z tables."""
    ds = sorted(os.path.basename(p)[:10] for p in glob.glob(f"{featdir}/*.json"))
    if dates is not None:
        s = set(dates); ds = [x for x in ds if x in s]
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(_day, [(featdir, x, zcols) for x in ds], chunksize=8))
    D = pd.DataFrame(dict(date=[r[0] for r in res], m_c1=[r[2] for r in res], m_data=[r[6] for r in res],
                          stations=[",".join(r[4]) if r[4] is not None else "" for r in res]))
    C = {r[0]: r[1] for r in res}; Zs = {r[0]: r[3] for r in res}; DA = {r[0]: r[5] for r in res}
    return D, C, Zs, DA


def count_matrix(D, C, grid=TG):
    """days x grid: control candidates >= T."""
    M = np.zeros((len(D), len(grid)))
    for i, ds in enumerate(D.date):
        c = C[ds]
        M[i] = len(c) - np.searchsorted(c, grid, side="left")
    return M


def stats(D, C, keep_days=None, grid=TG):
    """MED(T), POOLED(T) over D (m_c1 > 0), restricted to keep_days if given (TRIM = POOLED on kept days)."""
    E = D[D.m_c1 > 0]
    if keep_days is not None:
        E = E[E.date.isin(set(keep_days))]
    M = count_matrix(E, C, grid)
    r = M / E.m_c1.to_numpy()[:, None]
    return dict(MED=np.median(r, axis=0), POOLED=M.sum(0) / E.m_c1.sum(), n_days=len(E), n_fd=int(E.m_c1.sum()))


def at(curve, T, grid=TG):
    return float(curve[int(np.argmin(np.abs(grid - T)))])


def matched(curve, target, grid=TG):
    ok = np.where(curve <= target + 1e-12)[0]
    return float(grid[ok[0]]) if len(ok) else np.nan


def up25(T):
    return float(np.ceil(np.round(T / 0.25, 6)) * 0.25) if np.isfinite(T) else np.nan


def day_rate(D, C, T):
    E = D[D.m_c1 > 0]
    return pd.Series([(len(C[d]) - np.searchsorted(C[d], T, side="left")) / m for d, m in zip(E.date, E.m_c1)],
                     index=E.date.to_numpy())


def top_days(rate, frac=0.05):
    k = int(np.ceil(frac * len(rate)))
    return list(rate.sort_values(ascending=False, kind="mergesort").index[:k])
