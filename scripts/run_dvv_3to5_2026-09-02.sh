#!/bin/bash
# 3-5 s dv/v for every station of the last joint inversion (user request 2026-09-02). Same estimator + S-anchor as the
# 2-4 s production files; certified families only; eps-max 0.03 / n-eps 301 (cycle-skip-safe at 3-5 s, same 0.02% step).
# NEW files data/daily_dvv_<TAG>_Z_3to5.csv only — never touches 2-4 s results or the stacks.
cd /home/jovyan/tremorferometry
SC=/tmp/claude-1000/-home-jovyan-tremorferometry/03499166-a677-42d8-9397-e2cf5e28621f/scratchpad
one() {
  T=$1; s=$(echo $T | tr A-Z a-z | sed 's/p90f40//'); OUT=data/daily_dvv_${T}_Z_3to5.csv
  [ -s "$OUT" ] && { echo "$T exists, skip"; return; }
  C=data/${s}_causality_cert.csv; [ -f "$C" ] || C=$SC/${s}_causality_cert.csv
  bash scripts/run_capped.sh 100 python scripts/dvv_roll30cal.py --station $T --npz data/long_window_daily_${T}_Z.npz \
     --window 3 5 --origin-anchor --eps-max 0.03 --n-eps 301 --cert-csv $C --workers 14 --out $OUT > logs/dvv35_${s}.log 2>&1 \
     && echo "$(date +%H:%M) $T ok $(wc -l < $OUT) rows" || echo "$(date +%H:%M) $T FAILED"
}
export -f one; export SC
cat data/window35_station_list.txt | xargs -P 2 -I{} bash -c 'one {}'
echo "$(date) ALL DONE"
