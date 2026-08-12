#!/bin/bash
# B004 re-run with the P>=0.55 (22k) candidates, forward-only, AFTER B928's GPU phase finishes
# (no 3rd concurrent densify). Then CPU tail (fwd stacks -> fwd Z dv/v -> causality finalize + low-count flag).
cd /home/jovyan/tremorferometry
STATUS=logs/rollout_status.log
echo "$(date +%H:%M) B004 RE-RUN armed: waiting for B928 GPU PHASE done, then forward-only pipeline" >> "$STATUS"
until grep -q 'B928 GPU PHASE done' "$STATUS" 2>/dev/null; do sleep 120; done
echo "$(date +%H:%M) === B004 RE-RUN begin (forward-only, 22k cand) ===" >> "$STATUS"
bash scripts/borehole_gpu_phase_fwd.sh B004
rc=$?
if [ $rc -eq 0 ]; then
  bash scripts/borehole_cpu_tail.sh B004
else
  echo "$(date +%H:%M) B004 RE-RUN GPU phase EXITED rc=$rc" >> "$STATUS"
fi
echo "$(date +%H:%M) B004 RE-RUN complete" >> "$STATUS"
