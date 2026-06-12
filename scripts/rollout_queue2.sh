#!/bin/bash
# Trust-battery rollout — REMAINING 25 densified stations (anchors first), each: download -> battery
# (auto-cleans waveforms) -> next. Sequential (one station's waveforms at a time; one GPU job at a time).
cd /home/jovyan/tremorferometry
export PYTHONPATH=src
PY=/home/jovyan/envs/tremorferometry/bin/python
L=logs/rollout_queue2.log; : > $L
# per-station network / client / start  (anchors on CN/UW/IU via IRIS; PB boreholes via EARTHSCOPE)
declare -A NET CLI START
for S in PGC NLLB; do NET[$S]=CN; CLI[$S]=IRIS; START[$S]=2004-01-01; done
for S in GNW HDW COLT;  do NET[$S]=UW; CLI[$S]=IRIS; START[$S]=2004-01-01; done
NET[COR]=IU; CLI[COR]=IRIS; START[COR]=2004-01-01
ORDER="B004 B005 B011 B013 B020 B023 B026 B028 B030 B032 B033 B036 B039 B045 B049 B204 B927 B932 B941 PGC GNW NLLB HDW COLT COR"
for STA in $ORDER; do
  s=$(echo $STA|tr A-Z a-z)
  net=${NET[$STA]:-PB}; cli=${CLI[$STA]:-EARTHSCOPE}; st=${START[$STA]:-2005-06-21}
  if [ -f data/family_trust_tier1_${STA}.csv ]; then echo "$(date +%H:%M) $STA already done, skip" | tee -a $L; continue; fi
  echo "=== $(date +%H:%M) $STA ($net, $cli, from $st): download ===" | tee -a $L
  $PY scripts/download_station.py --network $net --station $STA --start $st --end 2026-06-11 \
      --workers 8 --client $cli >> logs/${s}_redl.log 2>&1
  echo "=== $(date +%H:%M) $STA: battery ===" | tee -a $L
  bash scripts/run_trust_battery.sh $net $STA "2015 2016" >> $L 2>&1
  grep -E "REAL fams|FAKES|verdicts" logs/trust_${s}.log 2>/dev/null | tail -3 | tee -a $L
done
echo "=== $(date +%H:%M) ROLLOUT QUEUE 2 DONE (25 stations) ===" | tee -a $L
