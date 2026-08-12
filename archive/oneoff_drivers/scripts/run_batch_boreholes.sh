#!/bin/bash
# Batch driver: run the borehole 3-comp dv/v pipeline sequentially for this batch.
# B003 skipped (no waveforms / no cand_filtered on disk — logged blocker).
cd /home/jovyan/tremorferometry
STATUS=logs/rollout_status.log
echo "$(date +%H:%M) BATCH start: B001 B004 B927 B928 (B003 blocked)" >> "$STATUS"
for STA in B001 B004 B927 B928; do
  echo "$(date +%H:%M) === $STA pipeline BEGIN ===" >> "$STATUS"
  bash scripts/borehole_pipeline.sh "$STA"
  rc=$?
  if [ $rc -ne 0 ]; then
    echo "$(date +%H:%M) $STA pipeline EXITED rc=$rc — skipping to next" >> "$STATUS"
  fi
done
echo "$(date +%H:%M) BATCH complete" >> "$STATUS"
