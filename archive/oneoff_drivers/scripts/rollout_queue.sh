#!/bin/bash
# Trust-battery rollout queue: for each PB station, full re-download -> battery (auto-cleans waveforms) -> next.
# Sequential (one station's waveforms at a time; one GPU job at a time). Runs unattended.
cd /home/jovyan/tremorferometry
export PYTHONPATH=src
PY=/home/jovyan/envs/tremorferometry/bin/python
L=logs/rollout_queue.log; : > $L
STAS="B935 B017 B201 B022 B001 B928 B014 B040"   # on-band / natural-dominated PB boreholes (B935 = positive control)
for STA in $STAS; do
  s=$(echo $STA|tr A-Z a-z)
  # skip if already certified this session
  if [ -f data/family_trust_tier1_${STA}.csv ]; then echo "$(date +%H:%M) $STA already done, skip" | tee -a $L; continue; fi
  echo "=== $(date +%H:%M) $STA: download ===" | tee -a $L
  $PY scripts/download_station.py --network PB --station $STA --start 2005-06-21 --end 2026-06-11 \
      --workers 8 --client EARTHSCOPE >> logs/${s}_redl.log 2>&1
  echo "=== $(date +%H:%M) $STA: battery ===" | tee -a $L
  bash scripts/run_trust_battery.sh PB $STA "2015 2016" >> $L 2>&1
  grep -E "REAL fams|FAKES|verdicts" logs/trust_${s}.log 2>/dev/null | tail -3 | tee -a $L
done
echo "=== $(date +%H:%M) ROLLOUT QUEUE DONE ===" | tee -a $L