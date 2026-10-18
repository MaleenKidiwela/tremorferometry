#!/usr/bin/env python
"""bw_firstlook_counts -- first look (sanity check, not science) at the SVI backward-extension matched-filter output.

Per UTC day and run: candidates with net_mad >= T for series 'data' and 'c1' (scrambled control: each station site
shifted to another time of the same day), plus observable family-days per series.  No HF veto (same as the P4
control-rate rule), no confidence model P (not yet applied), no catalogue matching (Bostock events are NOT read; only
the ETS window dates from backward_scope/truth_episodes.csv are used, in the plotting/table script).
Read-only on mf/feat; skips any day whose files are missing or unreadable (runs still being written).
Light: <= 4 worker processes.
"""
import os, glob, json
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor

BW = "/home/jovyan/tremorferometry/svi_science/results/backward"
FEAT = f"{BW}/mf/feat"
OUT = f"{BW}/first_look"
RUNS = {  # run -> detector label
    **{f"p5f_{y}": "frozen" for y in range(2003, 2014)},
    "p5e_2003_kelb": "era", "p5e_2004": "era", "p5e_2005": "era", "p5e_2006": "era",
}
# frozen detector: per-station-set T_core (p4/tables/p4_2_decision.json); boundaries = station on-dates (bw_gate_common)
SETS = [("CN4", "2003-01-01", 10.75), ("B1", "2005-07-01", 10.5), ("B2", "2005-11-02", 9.0),
        ("B3", "2007-11-13", 9.0), ("S15fix", "2009-08-27", 8.75)]
ERA_T = (8.75, 9.0)  # provisional: P4-3 era threshold not derived (p4_3_era.csv empty)
TGRID = sorted({t for _, _, t in SETS} | set(ERA_T))


def one(args):
    run, d = args
    p, q = f"{FEAT}/{run}/{d}.parquet", f"{FEAT}/{run}/{d}_fd.parquet"
    try:
        A = pd.read_parquet(p, columns=["series", "net_mad"]) if os.path.getsize(p) > 0 else pd.DataFrame(columns=["series", "net_mad"])
        F = pd.read_parquet(q, columns=["series", "observable"])
    except Exception:
        return None
    r = dict(run=run, date=d)
    for s in ("data", "c1"):
        m = A.net_mad[A.series == s].to_numpy(float)
        for T in TGRID:
            r[f"n_{s}_{T:g}"] = int((m >= T).sum())
        r[f"fd_{s}"] = int((F.series == s).astype(bool).mul(F.observable.astype(bool)).sum())
    return r


def main():
    jobs = []
    for run in RUNS:
        for j in sorted(glob.glob(f"{FEAT}/{run}/*.json")):
            jobs.append((run, os.path.basename(j)[:10]))
    with ProcessPoolExecutor(4) as ex:
        R = [r for r in ex.map(one, jobs, chunksize=32) if r is not None]
    D = pd.DataFrame(R)
    D["detector"] = D.run.map(RUNS)
    dd = pd.to_datetime(D.date)
    D["sset"] = "era(S15ERA)"
    D["T_core"] = np.nan
    fz = D.detector == "frozen"
    for name, start, T in SETS:
        k = fz & (dd >= pd.Timestamp(start))
        D.loc[k, "sset"] = name; D.loc[k, "T_core"] = T
    for s in ("data", "c1"):
        D[f"n_{s}_Tcore"] = np.nan
        for name, start, T in SETS:
            k = fz & (D.sset == name)
            D.loc[k, f"n_{s}_Tcore"] = D.loc[k, f"n_{s}_{T:g}"]
    D.to_csv(f"{OUT}/daily_counts.csv", index=False)
    print(D.groupby("run").size().to_string(), "\nrows", len(D))


if __name__ == "__main__":
    main()
