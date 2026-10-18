#!/usr/bin/env python
"""q4_core9 (2026-10-07) -- Q4 column: the fixed 9-borehole run (wa_fixed10/post_v2_core9), WAVE-only P_own >= 0.75,
raw counting (no dedup flags exist for it), headline deep families, big-ETS clock, D1, all majors dropped.
Uses the functions of q4_recharge_wa.py unchanged.  Rules: q4_recharge_wa/PLAN.md ("Fixed 9-borehole run").
    python svi_science/scripts/q4_core9.py
"""
import sys, json
import numpy as np
import pandas as pd
from multiprocessing import Pool
sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import q4_recharge_wa as Q

C9 = f"{Q.ROOT}/lfe_extension_wa/results/wa_fixed10/post_v2_core9"


def main():
    deep = json.load(open(f"{Q.OUT}/q4_summary.json"))["deep_sets"]["deep"]
    R = pd.read_parquet(f"{C9}/rates_daily_v2.parquet", columns=["unit_kind", "label", "date", "fid", "region",
                                                                "n075_data", "n075_rev", "n075_c1"])
    R = R[(R.unit_kind == "family") & (R.region == "NW") & (R.label == "own") & R.fid.isin(deep)
          & (R.date >= Q.P0) & (R.date <= Q.P1)]
    days = np.arange(np.datetime64(Q.P0), np.datetime64(Q.P1) + 1)
    ds = pd.Series(days.astype(str))
    g = R.groupby("date")
    n = g.size().reindex(ds, fill_value=0).to_numpy(float)
    Epi = pd.read_csv(f"{Q.EPR}/episodes_nw_final.csv", parse_dates=["start", "end"])
    maj = Epi[Epi.major].sort_values("start").reset_index(drop=True)
    Bp = pd.read_csv(Q.B_EP, parse_dates=["start", "end"])
    big = Bp[(Bp.region == "NW") & (Bp.big == True)].sort_values("start").reset_index(drop=True)
    jobs = []
    for col, lab in (("n075_data", "core9 deep (raw)"), ("n075_rev", "core9 time-reversed control (raw)"),
                     ("n075_c1", "core9 scrambled control (raw)")):
        y = g[col].sum().reindex(ds, fill_value=0).to_numpy(float)
        Q.log(f"{lab}: events per era all days { {e: int(y[(ds >= a) & (ds <= b)].sum()) for e, (a, b) in Q.ERAS.items()} }")
        C = Q.cycle_table(days, y, n, big, maj)
        nb = Q.NB if col == "n075_data" else 0
        for e in Q.ERAS:
            jobs.append((f"{lab} | clock big | D1 | {e}", Q.era_sub(C, e), True, nb))
        jobs.append((f"{lab} | clock big | D1 | pooled-normalised", Q.pooled_norm(C), True, nb))
    with Pool(12) as p:
        res = p.map(Q.analyse, jobs, chunksize=1)
    O = pd.DataFrame([{k: (json.dumps(v) if isinstance(v, list) else v) for k, v in r[0].items()} for r in res])
    O.to_csv(f"{Q.OUT}/recharge_fits_core9.csv", index=False)
    pd.concat([r[1] for r in res]).to_csv(f"{Q.OUT}/recharge_curves_core9.csv", index=False)
    for r in res:
        Q.log(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r[0].items()}, default=str))


if __name__ == "__main__":
    main()
