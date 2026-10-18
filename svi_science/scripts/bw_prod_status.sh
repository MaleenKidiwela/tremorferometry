#!/usr/bin/env bash
# one-shot progress of the SVI backward production queues (no waiting)
O=/home/jovyan/tremorferometry/svi_science/results/backward
TZ=America/Los_Angeles date
echo "markers: $(cd $O && ls P*_DONE P*_FAILED 2>/dev/null | tr '\n' ' ')"
for q in P1A P1B P2 P3 P5; do p=$(cat $O/${q}_PID 2>/dev/null); echo -n "$q pid $p $(ps -p $p >/dev/null 2>&1 && echo alive || echo stopped); "; done; echo
python - <<PY
import pandas as pd, os
for q in "AB":
    f = "$O/tables/ledger_%s.csv" % q
    if os.path.exists(f):
        L = pd.read_csv(f, dtype=str); L = L[L.chan_log != "resumed"]
        print(q, len(L), "station-days;", L.groupby("status").size().to_dict(), "; last", L.date.iloc[-1], L.sta.iloc[-1])
PY
du -sh /home/jovyan/tremorferometry/lfe_extension/traces_backward /home/jovyan/tremorferometry/lfe_extension/cache20hz_backward 2>/dev/null
tail -n 1 $O/logs/p1a.log $O/logs/p1b.log 2>/dev/null | cut -c1-160
