#!/bin/bash
# Takeover driver: let B001 finish serially (orphaned pipeline), then run B004/B927/B928 with the
# CONCURRENT densify variant. B003 already dropped (no data). Skip-and-continue on abort.
cd /home/jovyan/tremorferometry
STATUS=logs/rollout_status.log
echo "$(date +%H:%M) CONC driver armed: waiting for B001 to finish, then B004/B927/B928 concurrent" >> "$STATUS"
# wait for B001 pipeline to fully finish (finalize or abort)
until grep -qE 'B001 (FINALIZE done|ABORT)' "$STATUS" 2>/dev/null; do sleep 60; done
echo "$(date +%H:%M) B001 finished -> starting concurrent batch B004/B927/B928" >> "$STATUS"
for STA in B004 B927 B928; do
  echo "$(date +%H:%M) === $STA pipeline BEGIN (concurrent densify) ===" >> "$STATUS"
  bash scripts/borehole_pipeline_conc.sh "$STA"
  rc=$?
  [ $rc -ne 0 ] && echo "$(date +%H:%M) $STA pipeline EXITED rc=$rc — skipping to next" >> "$STATUS"
done
echo "$(date +%H:%M) CONC BATCH complete" >> "$STATUS"
