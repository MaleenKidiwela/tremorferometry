#!/usr/bin/env python
"""P4-2b step 3: sidecar flag tables for the P5 (and ladder) feature files -- the feature files are NOT rewritten.
  burst_day      = the day's control rate (all stored c1 candidates >= 8 MAD / observable c1 family-days) >= 5x the
                   median daily control rate of that run's calendar year (mean if the median is 0; ref_stat says which)
  multi_site_z2  = n_z2 >= 2 (single-site column; n_z2 copied)
Key: run, date, row (0-based row of <date>.parquet); series, fam, t, net_mad copied for checking.
Only FINISHED runs (year marker present); rerun is idempotent (skips a run whose day count is unchanged).
Rules: results/backward/p4b/P4B.md.

  python bw_p4b_flags.py [--force] [run ...]   -> p4b/sidecar/<run>_flags.parquet, <run>_days.csv, <run>_meta.json
"""
import os, sys, glob, json, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd

BW = "/home/jovyan/tremorferometry/svi_science/results/backward"
OUT = f"{BW}/p4b/sidecar"
COLS = ["series", "fam", "t", "net_mad", "n_z2"]


def runs_finished():
    out = {}
    for d in sorted(glob.glob(f"{BW}/mf/feat/*")):
        r = os.path.basename(d)
        if r == "p5e_2003_kelb":
            m = f"{BW}/P5E_2003_KELB_DONE"
        elif r.startswith("p5f_"):
            m = f"{BW}/P5F_{r[4:]}_DONE"
        elif r.startswith("p5e_"):
            m = f"{BW}/P5E_{r[4:]}_DONE"
        else:
            continue
        if os.path.exists(m):
            out[r] = d
    for d in sorted(glob.glob(f"{BW}/p4b/mf/feat/lad_*")):
        r = os.path.basename(d)
        if os.path.exists(f"{BW}/p4b/LAD_{r.split('_')[1]}_DONE"):
            out[r] = d
    return out


def one_day(args):
    d, ds = args
    p, q = f"{d}/{ds}.parquet", f"{d}/{ds}_fd.parquet"
    A = pd.DataFrame(columns=COLS)
    if os.path.exists(p) and os.path.getsize(p) > 0:
        x = pd.read_parquet(p)
        if len(x):
            A = x[[c for c in COLS if c in x.columns]].copy()
    m_c1 = m_data = 0
    if os.path.exists(q) and os.path.getsize(q) > 0:
        f = pd.read_parquet(q, columns=["series", "observable"])
        ob = f.observable.astype(bool)
        m_c1 = int(((f.series == "c1") & ob).sum()); m_data = int(((f.series == "data") & ob).sum())
    A = A.reset_index(drop=True); A.insert(0, "row", np.arange(len(A), dtype=np.int32)); A.insert(0, "date", ds)
    n_c1 = int((A.series == "c1").sum()) if len(A) else 0
    n_data = int((A.series == "data").sum()) if len(A) else 0
    return A, dict(date=ds, n_data=n_data, n_c1=n_c1, m_data=m_data, m_c1=m_c1)


def do_run(run, d, workers=6):
    ds = sorted(os.path.basename(p)[:10] for p in glob.glob(f"{d}/*.json"))
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(one_day, [(d, x) for x in ds], chunksize=8))
    A = pd.concat([r[0] for r in res if len(r[0])], ignore_index=True) if any(len(r[0]) for r in res) else pd.DataFrame()
    D = pd.DataFrame([r[1] for r in res])
    D["ctl_rate"] = np.where(D.m_c1 > 0, D.n_c1 / D.m_c1.clip(lower=1), np.nan)
    D["year"] = D.date.str[:4]
    ref, stat = {}, {}
    for y, g in D[D.m_c1 > 0].groupby("year"):
        med = float(g.ctl_rate.median())
        if med > 0:
            ref[y], stat[y] = med, "median"
        else:
            ref[y], stat[y] = float(g.ctl_rate.mean()), "mean"
    D["ctl_rate_ref"] = D.year.map(ref); D["ref_stat"] = D.year.map(stat)
    D["burst_ratio"] = D.ctl_rate / D.ctl_rate_ref
    D["burst_day"] = (D.ctl_rate_ref > 0) & (D.ctl_rate >= 5 * D.ctl_rate_ref)
    D.insert(0, "run", run)
    if len(A):
        A.insert(0, "run", run)
        A = A.merge(D[["date", "burst_day", "burst_ratio"]], on="date", how="left")
        A["n_z2"] = A.n_z2.astype("Int16")
        A["multi_site_z2"] = (A.n_z2 >= 2).fillna(False).astype(bool)
        A["burst_day"] = A.burst_day.fillna(False).astype(bool)
    A.to_parquet(f"{OUT}/{run}_flags.parquet", index=False)
    D.drop(columns="year").to_csv(f"{OUT}/{run}_days.csv", index=False)
    meta = dict(run=run, featdir=d, n_days=len(ds), n_rows=len(A), n_burst_days=int(D.burst_day.sum()),
                burst_days=D.loc[D.burst_day, "date"].tolist(), ref_per_year=ref, ref_stat=stat,
                frac_data_multi_site_z2=float(A.loc[A.series == "data", "multi_site_z2"].mean()) if len(A) else None,
                frac_c1_multi_site_z2=float(A.loc[A.series == "c1", "multi_site_z2"].mean()) if len(A) else None,
                written=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
    json.dump(meta, open(f"{OUT}/{run}_meta.json", "w"), indent=1)
    return meta


def main():
    os.makedirs(OUT, exist_ok=True)
    import fcntl
    lk = open(f"{OUT}/.lock", "w"); fcntl.flock(lk, fcntl.LOCK_EX)          # one flags pass at a time (queue + manual)
    force = "--force" in sys.argv
    want = [a for a in sys.argv[1:] if not a.startswith("--")]
    fin = runs_finished()
    for run, d in fin.items():
        if want and run not in want:
            continue
        mp = f"{OUT}/{run}_meta.json"
        nd = len(glob.glob(f"{d}/*.json"))
        if not force and os.path.exists(mp) and json.load(open(mp)).get("n_days") == nd:
            continue
        t0 = time.time(); m = do_run(run, d)
        print(f"[flags] {run}: {m['n_days']} days, {m['n_rows']} rows, burst days {m['n_burst_days']}, "
              f"multi_site_z2 data {m['frac_data_multi_site_z2']}, c1 {m['frac_c1_multi_site_z2']} ({time.time()-t0:.0f} s)",
              flush=True)


if __name__ == "__main__":
    main()
