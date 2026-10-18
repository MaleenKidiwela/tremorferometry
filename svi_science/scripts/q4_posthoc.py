#!/usr/bin/env python
"""q4_posthoc (2026-10-07) -- NOT pre-registered, descriptive: (1) the headline Q4 tables with each day counted as
active (>= 1 deep event) or not, so a burst counts once per day; (2) the busiest E3 days.
    python svi_science/scripts/q4_posthoc.py
"""
import sys, json
import numpy as np
import pandas as pd
from multiprocessing import Pool
sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import q4_recharge_wa as Q


def main():
    deep = json.load(open(f"{Q.OUT}/q4_summary.json"))["deep_sets"]["deep"]
    E = pd.read_parquet(f"{Q.OUT}/events_p075.parquet")
    E = E[(E.series == "data") & E.fid.isin(deep) & E.once_D]
    days = np.arange(np.datetime64(Q.P0), np.datetime64(Q.P1) + 1)
    ds = pd.Series(days.astype(str))
    y = E.groupby("date").size().reindex(ds, fill_value=0).to_numpy(float)
    Epi = pd.read_csv(f"{Q.EPR}/episodes_nw_final.csv", parse_dates=["start", "end"])
    maj = Epi[Epi.major].sort_values("start").reset_index(drop=True)
    Bp = pd.read_csv(Q.B_EP, parse_dates=["start", "end"])
    big = Bp[(Bp.region == "NW") & (Bp.big == True)].sort_values("start").reset_index(drop=True)
    occ = (y > 0).astype(float)
    n1 = np.ones(days.size)
    C = Q.cycle_table(days, occ, n1, big, maj)
    Cy = Q.cycle_table(days, y, n1, big, maj)
    jobs = [(f"POSTHOC active days | clock big | D1 | {e}", Q.era_sub(C, e), True, 300) for e in Q.ERAS]
    jobs.append(("POSTHOC active days | clock big | D1 | pooled-normalised", Q.pooled_norm(C), True, 300))
    # SVI, same day counting (P3 data; big clock; all majors dropped = Q4 definition)
    Fs = pd.read_csv(f"{Q.S0}/families.csv")
    sdeep = Fs.fid[Fs.group.isin([1, 4, 10, 12])].to_numpy()
    Es = pd.read_csv(f"{Q.EPR}/episodes_svi_final.csv", parse_dates=["start", "end"])
    maj39 = Es[Es.major].sort_values("start").reset_index(drop=True)
    Bs = Bp[Bp.region == "SVI"]
    sbig = Bs[(Bs.big == True) & (Bs.start >= "2013-01-01")].sort_values("start").reset_index(drop=True)
    T = pd.read_parquet(f"{Q.S0}/analysis_table.parquet", columns=["date", "series", "fam", "P_same", "keep_xfam8"])
    T = T[(T.date >= "2014-01-01") & (T.date <= Q.P1) & (T.P_same >= 0.75) & T.keep_xfam8 & (T.series == "data")]
    sdays = np.arange(np.datetime64("2014-01-01"), np.datetime64(Q.P1) + 1)
    ys = T[T.fam.isin(sdeep)].groupby("date").size().reindex(pd.Series(sdays.astype(str)), fill_value=0).to_numpy(float)
    both = pd.concat([maj39[["start", "end"]], sbig[["start", "end"]]])
    jobs.append(("POSTHOC SVI active days | clock big | D1 | drop all majors",
                 Q.cycle_table(sdays, (ys > 0).astype(float), np.ones(sdays.size), sbig, both), True, 300))
    jobs.append(("POSTHOC SVI active days | clock big | D1 | P3 as coded",
                 Q.cycle_table(sdays, (ys > 0).astype(float), np.ones(sdays.size), sbig, sbig), True, 300))
    with Pool(4) as p:
        res = p.map(Q.analyse, jobs, chunksize=1)
    O = pd.DataFrame([{k: (json.dumps(v) if isinstance(v, list) else v) for k, v in r[0].items()} for r in res])
    O.to_csv(f"{Q.OUT}/recharge_fits_posthoc_daycount.csv", index=False)
    pd.concat([r[1] for r in res]).to_csv(f"{Q.OUT}/recharge_curves_posthoc_daycount.csv", index=False)
    for r in res:
        Q.log(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r[0].items()}, default=str))
    # busiest usable days (E3), with t
    top = Cy[Cy.era == "E3"].sort_values("y", ascending=False).head(12)
    top["date"] = top.date.astype(str)
    tot = Cy[Cy.era == "E3"].y.sum()
    top["share_of_E3"] = top.y / tot
    top.to_csv(f"{Q.OUT}/busiest_days_E3.csv", index=False)
    Q.log("POSTHOC busiest E3 usable days:\n" + top.to_string(index=False))


if __name__ == "__main__":
    main()
