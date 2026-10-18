#!/usr/bin/env bash
# Family ACTIVITY by stacking — queue (design: PLAN.md §Stack activity; results/stack_activity/STACK_ACTIVITY_RESULTS.md)
#   nohup setsid -f bash scripts/stack_activity_queue.sh > logs/stack_activity_queue.out 2>&1 < /dev/null
#   results/STACK_ACT_PID  queue pid      results/STACK_ACT_HEARTBEAT  "<utc> <stage>" every 60 s
#   results/STACK_ACT_DONE / results/STACK_ACT_FAILED (stage + log tail)
#   results/stack_activity/stage_<name>.done  per-stage guard (delete one to re-run it); resumable
#   stop: kill $(cat results/STACK_ACT_PID)   (exact pid; the heartbeat child dies with the trap)
# Phase A (now): valA = V2008 + V2010 validation on the 42 trial stations (cache already on disk).
# Phase B waits for step 4b (results/STEP4B_DONE or STEP4B_FAILED or its pid gone), then:
#   valB = A2011 (independent episode, build templates), valC = V2010 + step-4b stations,
#   tr2012 / tr2013 / tr2014 trial episodes, selection of the continuous-series families + stations,
#   fetch + cache 2012-06..2014-12 for those stations, valS / valSB (validation of that exact configuration),
#   cont (continuous series), report.
set -u
S=/home/jovyan/tremorferometry/lfe_extension_wa
cd "$S" || exit 1
R=$S/results; L=$S/logs; G=$R/stack_activity
mkdir -p "$G" "$L"
exec 9>"$R/.stack_act.lock"
if ! flock -n 9; then echo "another stack_activity_queue holds the lock"; exit 1; fi
echo $$ > "$R/STACK_ACT_PID"
rm -f "$R/STACK_ACT_DONE" "$R/STACK_ACT_FAILED"
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
( while true; do echo "$(date -u +%FT%TZ) $(cat $G/.stage 2>/dev/null)" > "$R/STACK_ACT_HEARTBEAT"; sleep 60; done ) &
HB=$!
trap 'kill $HB 2>/dev/null' EXIT

run() {  # run <name> <cmd...>
  local name=$1; shift
  if [ -f "$G/stage_$name.done" ]; then echo "[skip] $name"; return 0; fi
  echo "$name" > "$G/.stage"
  echo "$(date -u +%FT%TZ) [stage] $name: $*" | tee -a "$L/stack_activity_queue.log"
  if "$@" >> "$L/stack_activity_$name.log" 2>&1; then
    date -u +%FT%TZ > "$G/stage_$name.done"
    echo "$(date -u +%FT%TZ) [ok] $name" | tee -a "$L/stack_activity_queue.log"
  else
    { echo "stage $name failed $(date -u +%FT%TZ)"; tail -40 "$L/stack_activity_$name.log"; } > "$R/STACK_ACT_FAILED"
    echo "$(date -u +%FT%TZ) [FAILED] $name" | tee -a "$L/stack_activity_queue.log"
    exit 1
  fi
}

df -h /home/jovyan | tail -1 | tee -a "$L/stack_activity_queue.log"
SA="python scripts/stack_activity.py"
SC="python scripts/stack_activity_score.py"
ST0=$R/setup; ST4=$R/step4b/setup
# ---------------- phase A
run sa_valA  $SA --tag valA --windows V2008,V2010 --tpl build --stations $ST0/stations.csv --scored test
run sc_valA  $SC --run valA --mode validation --truth $ST0/truth.parquet --scan $ST0/scan_days.csv
run report_A python scripts/stack_activity_report.py

# ---------------- wait for step 4b (its data: A2011 / T2013 / T2014 caches, step4b_a templates)
echo "wait_step4b" > "$G/.stage"
while [ ! -f "$R/STEP4B_DONE" ] && [ ! -f "$R/STEP4B_FAILED" ] && kill -0 3740832 2>/dev/null; do sleep 120; done
echo "$(date -u +%FT%TZ) step 4b finished: $(ls $R/STEP4B_DONE $R/STEP4B_FAILED 2>/dev/null)" | tee -a "$L/stack_activity_queue.log"
if [ -f "$R/templates/tpl_step4b_a.npz" ]; then T4=build,step4b_a; else T4=build; fi
echo "templates for step-4b runs: $T4" | tee -a "$L/stack_activity_queue.log"

run sa_valB  $SA --tag valB --windows A2011 --tpl build --stations $ST4/stations.csv --scored all
run sc_valB  $SC --run valB --mode validation --truth $ST4/truth.parquet --scan $ST4/scan_days.csv
run sa_valC  $SA --tag valC --windows V2010 --tpl $T4 --stations $ST4/stations.csv --scored test
run sc_valC  $SC --run valC --mode validation --truth $ST4/truth.parquet --scan $ST4/scan_days.csv
for E in 2012 2013 2014; do
  run sa_tr$E $SA --tag tr$E --windows T$E --tpl $T4 --stations $ST4/stations.csv --scored all
  run sc_tr$E $SC --run tr$E --mode trial
done
run report_B python scripts/stack_activity_report.py

# ---------------- continuous monthly series 2012-06 -> 2014-12 for the best-validated families
run select   python scripts/stack_activity_select.py
df -h /home/jovyan | tail -1 | tee -a "$L/stack_activity_queue.log"
run fetch_cont  env WA_SETUP=$G/setup_cont python scripts/wa_fetch.py --workers 8
run fetch_cont2 env WA_SETUP=$G/setup_cont python scripts/wa_fetch.py --workers 4
run cache_cont  env WA_SETUP=$G/setup_cont WA_EXTRA_NODES=$ST4/extra_nodes.csv python scripts/wa_cache.py --workers 12
FS="--fids $G/cont_fids.csv --only-stations $G/cont_stations.csv --stations $ST4/stations.csv"
run sa_valS  $SA --tag valS --windows V2008,V2010 --tpl $T4 $FS --scored test
run sc_valS  $SC --run valS --mode validation --truth $ST0/truth.parquet --scan $ST0/scan_days.csv
run sa_valSB $SA --tag valSB --windows A2011 --tpl build $FS --scored all
run sc_valSB $SC --run valSB --mode validation --truth $ST4/truth.parquet --scan $ST4/scan_days.csv
run sa_cont  $SA --tag cont --date-range 2012-06-01:2014-12-31 --tpl $T4 $FS --scored all
run sc_cont  $SC --run cont --mode trial
run report   python scripts/stack_activity_report.py
echo "done" > "$G/.stage"
date -u +%FT%TZ > "$R/STACK_ACT_DONE"
echo "$(date -u +%FT%TZ) [DONE]" | tee -a "$L/stack_activity_queue.log"
