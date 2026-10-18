#!/usr/bin/env python
"""P4-2b step 5 setup: completeness-ladder day lists (60 d from each station-set boundary, split by year so each year
waits for its own P2 marker) and exclusion files (every S15 station outside the OLD set, every month touched), for the
frozen detector run with --sset eff (the P4-2 exclusion-file mechanism).  Rules: results/backward/p4b/P4B.md."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
import bw_gate_common as G

O = "/home/jovyan/tremorferometry/svi_science/results/backward/p4b/tables"
S15 = ["PGC", "VGZ", "NLLB", "SNB", "LRIV", "B001", "B003", "B004", "B005", "B006", "B007", "B009", "B010", "B011", "B926"]
LAD = [("CN4", G.CN4, "2005-07-01"), ("B1", G.B1, "2005-11-02"), ("B2", G.B2, "2007-11-13"), ("B3", G.B3, "2009-08-27")]
for name, mem, b in LAD:
    days = [str(x.date()) for x in pd.date_range(b, periods=60)]
    run = f"lad_{name}_{b.replace('-', '')}"
    for y in sorted({d[:4] for d in days}):
        json.dump([d for d in days if d[:4] == y], open(f"{O}/{run}_days_{y}.json", "w"))
    months = sorted({d[:7] for d in days})
    E = pd.DataFrame([dict(sta=s, month=m, reason=f"p4b_ladder_old_{name}") for s in S15 if s not in mem for m in months])
    E.to_csv(f"{O}/{run}_excl.csv", index=False)
    print(run, days[0], days[-1], "years", sorted({d[:4] for d in days}), "kept", [s for s in S15 if s in mem], "excl rows", len(E))
