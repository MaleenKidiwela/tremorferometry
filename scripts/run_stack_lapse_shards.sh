#!/bin/bash
cd /home/jovyan/tremorferometry
N=6
for k in $(seq 0 $((N-1))); do SHARD=$k NSHARD=$N nohup bash scripts/run_capped.sh 60 python fault_tomography/inversion/sse_stack_lapse_fleet.py > logs/stack_lapse_s$k.log 2>&1 & done
wait
echo "$(date) all shards finished" >> logs/stack_lapse_fleet.log; echo FLEET_STACK_DONE >> logs/stack_lapse_fleet.log
