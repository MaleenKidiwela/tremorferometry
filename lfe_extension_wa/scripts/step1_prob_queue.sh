#!/usr/bin/env bash
# STEP 1 (Merlin, 2026-09-30): per-repeat probability by local FDR, scrambled control = empirical noise model.
# Guarded, resumable stages:  mf_c34 (GPU, 2 extra site-rolled scrambled controls) -> feat (CPU) -> fit (CPU) -> report.
#   nohup setsid -f bash scripts/step1_prob_queue.sh > logs/step1_prob_queue.out 2>&1 < /dev/null
#   results/STEP1_PID  results/STEP1_HEARTBEAT ("<utc> <stage>" every 60 s)  results/STEP1_DONE / results/STEP1_FAILED
#   results/step1_prob/stage_<name>.done  (delete one to re-run that stage)
#   stop: kill $(cat results/STEP1_PID)   (exact pid; the heartbeat is killed by the trap)
set -u
S=/home/jovyan/tremorferometry/lfe_extension_wa
cd "$S" || exit 1
R=$S/results; L=$S/logs; G=$R/step1_prob
mkdir -p "$G" "$L"
exec 9>"$R/.step1.lock"
if ! flock -n 9; then echo "another step1_prob_queue holds the lock"; exit 1; fi
echo $$ > "$R/STEP1_PID"
rm -f "$R/STEP1_DONE" "$R/STEP1_FAILED"
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
export WA_SETUP=$R/merlin_tests/t1_setup
echo start > "$G/.stage"
( while true; do echo "$(date -u +%FT%TZ) $(cat $G/.stage 2>/dev/null)" > "$R/STEP1_HEARTBEAT"; sleep 60; done ) &
HB=$!
trap 'kill $HB 2>/dev/null' EXIT
run() {
  local name=$1; shift
  if [ -f "$G/stage_$name.done" ]; then echo "[skip] $name"; return 0; fi
  echo "$name" > "$G/.stage"
  echo "$(date -u +%FT%TZ) [stage] $name: $*" | tee -a "$L/step1.log"
  if "$@" >> "$L/step1_$name.log" 2>&1; then
    date -u +%FT%TZ > "$G/stage_$name.done"; echo "$(date -u +%FT%TZ) [ok] $name" | tee -a "$L/step1.log"
  else
    { echo "stage $name failed $(date -u +%FT%TZ)"; tail -40 "$L/step1_$name.log"; } > "$R/STEP1_FAILED"
    echo "$(date -u +%FT%TZ) [FAILED] $name" | tee -a "$L/step1.log"; exit 1
  fi
}
W3="V2008,V2010,A2011"
run mf_c34 python scripts/wa_mf.py --tag step1c34 --windows $W3 --tpl build --n-control 2 --drop-dead 0.05 --all-templates --seed 20261001
run feat   python scripts/step1_features.py --workers 12 --windows $W3
run fit    env OMP_NUM_THREADS=16 python scripts/step1_fit.py
run report python scripts/step1_report.py
rm -f "$G/stage_report.done"
date -u +%FT%TZ > "$R/STEP1_DONE"
echo "$(date -u +%FT%TZ) [DONE]" | tee -a "$L/step1.log"
