#!/bin/bash
# Mechanism test (2026-09-04): re-measure 2-4 s dv/v WITHOUT the SVD-Wiener filter at 6 strong-ETS-response stations.
# NEW files data/daily_dvv_<TAG>_Z_2to4_nosvd.csv only; production files untouched.
cd /home/jovyan/tremorferometry
for T in ANGUS B004p90f40 LOPEZ B024p90f40 B928p90f40 B011p90f40; do
  s=$(echo $T | tr A-Z a-z | sed 's/p90f40//'); C=data/${s}_causality_cert.csv; [ -f "$C" ] || C=/tmp/claude-1000/-home-jovyan-tremorferometry/03499166-a677-42d8-9397-e2cf5e28621f/scratchpad/${s}_causality_cert.csv
  bash scripts/run_capped.sh 100 python scripts/dvv_roll30cal.py --station $T --npz data/long_window_daily_${T}_Z.npz --window 2 4 --origin-anchor --no-svd \
     --cert-csv $C --workers 14 --out data/daily_dvv_${T}_Z_2to4_nosvd.csv > logs/dvv_nosvd_${s}.log 2>&1 && echo "$(date +%H:%M) $T ok" || echo "$T FAILED"
done
echo NOSVD_DONE
