#!/bin/bash
# Pipelined driver: GPU densify of station N+1 overlaps the CPU tail (stacks/dv/v/finalize) of station N.
# Only the GPU-densify stages are serialized (wait_gpu_free inside gpu_phase). CPU tails run in background.
# B004/B927/B928 (B001 handled by its own orphaned pipeline; B003 dropped). Skip-and-continue on abort.
cd /home/jovyan/tremorferometry
STATUS=logs/rollout_status.log
echo "$(date +%H:%M) PIPELINED driver start: B004 -> B927 -> B928 (densify overlaps prev CPU tail)" >> "$STATUS"
for STA in B004 B927 B928; do
  echo "$(date +%H:%M) === $STA GPU PHASE begin ===" >> "$STATUS"
  bash scripts/borehole_gpu_phase.sh "$STA"
  rc=$?
  if [ $rc -eq 0 ]; then
    echo "$(date +%H:%M) $STA -> launching CPU tail in background" >> "$STATUS"
    nohup bash scripts/borehole_cpu_tail.sh "$STA" > logs/tail_${STA}.log 2>&1 &
  else
    echo "$(date +%H:%M) $STA GPU phase EXITED rc=$rc — skipping to next" >> "$STATUS"
  fi
done
echo "$(date +%H:%M) all GPU phases done; waiting for background CPU tails to finalize" >> "$STATUS"
wait
echo "$(date +%H:%M) PIPELINED BATCH complete" >> "$STATUS"
