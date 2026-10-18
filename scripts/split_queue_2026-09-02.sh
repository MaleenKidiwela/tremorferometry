#!/bin/bash
# Merlin split-ruling steps 4-5, queued behind the three gate chains (2026-09-02). Gate 0.3 population.
cd /home/jovyan/tremorferometry
RUN="bash scripts/run_capped.sh 200 python fault_tomography/inversion/run_subset_chain.py --lam 0.464 --drop-fe data/family_drift_gate_D0.3.csv"
until grep -q -E "DONE ->|stage failed" logs/split_gateD0.2.log 2>/dev/null && grep -q -E "DONE ->|stage failed" logs/split_gateD0.5.log 2>/dev/null; do sleep 120; done
echo "$(date) gate chains finished; starting design (E) family split"
$RUN --name E_odd  --parity odd  --workers 14 > logs/split_E_odd.log 2>&1 &
$RUN --name E_even --parity even --workers 14 > logs/split_E_even.log 2>&1 &
wait
echo "$(date) (E) done; starting leave-one-out (top-5 gated q16 support contributors)"
for pair in "B943p90f40 B018p90f40" "NLWA STOR" "B005p90f40"; do
  for s in $pair; do
    echo "$s" > data/loo_exclude_$s.txt
    $RUN --name loo_$s --exclude data/loo_exclude_$s.txt --workers 14 > logs/split_loo_$s.log 2>&1 &
  done
  wait
done
echo "$(date) leave-one-out done"
