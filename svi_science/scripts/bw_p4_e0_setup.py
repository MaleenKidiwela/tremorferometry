#!/usr/bin/env python
"""P4-2 setup: E0 (2014-01..2017-02) PNSN-quiet days + one exclusion file per reduced station set
(= forward R2 exclusions U every station outside the set, every E0 month), so the forward detector run with
--sset eff differs from the forward production only by the station set.  Writes p4/tables/."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import bw_p4_common as P
import bw_gate_common as G

Q = P.pnsn_quiet("2014-01-01", "2017-02-28")
days = sorted(Q[Q.quiet == True].date)                                    # noqa: E712
os.makedirs(P.TAB, exist_ok=True)
json.dump(days, open(f"{P.TAB}/p4_e0_quiet_days.json", "w"))
# cross-check against Track C's stored quiet_flag (must be the same rule)
import glob
fs = [f for f in sorted(glob.glob(f"{P.LFE}/results/trackC/family_days_*.parquet")) if "2014-01" <= f[-15:-8] <= "2017-02"]
tc = pd.concat([pd.read_parquet(f, columns=["date", "quiet_flag"]) for f in fs]).groupby("date").quiet_flag.max()
tcq = set(tc[tc].index)
print(f"E0 quiet days: {len(days)}; Track C quiet_flag days {len(tcq)}; symmetric difference {len(set(days) ^ tcq)}")
R2 = pd.read_csv(f"{P.REFINE}/r2/r2_exclusions.csv", dtype=str).fillna("")
months = sorted({d[:7] for d in days})
S15 = ["B001", "B003", "B004", "B005", "B006", "B007", "B009", "B010", "B011", "B926", "LRIV", "NLLB", "PGC", "SNB", "VGZ"]
for name, mem in dict(CN4=G.CN4, B1=G.B1, B2=G.B2, B3=G.B3).items():
    add = pd.DataFrame([dict(sta=s, month=m, reason=f"p4_set_{name}") for s in S15 if s not in mem for m in months])
    E = pd.concat([R2, add], ignore_index=True).drop_duplicates(["sta", "month"])
    E.to_csv(f"{P.TAB}/p4_excl_{name}.csv", index=False)
    print(name, "stations", mem, "exclusion rows", len(E))
