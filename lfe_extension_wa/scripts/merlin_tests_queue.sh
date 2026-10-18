#!/usr/bin/env bash
# Merlin's two deciding tests (2026-09-30): T1 all-templates matched filter, T2 withheld-station template recovery.
# Two INDEPENDENT chains run side by side (T2 failing never stops T1, and vice versa); the report runs at the end
# with whatever exists.  Guarded, resumable stages.
#   nohup setsid -f bash scripts/merlin_tests_queue.sh > logs/merlin_tests_queue.out 2>&1 < /dev/null
#   results/MERLIN_T_PID        queue pid            results/MERLIN_T_HEARTBEAT  "<utc> T1:<stage> T2:<stage>" every 60 s
#   results/MERLIN_T_DONE       both chains + report results/MERLIN_T_FAILED     which chain / stage failed + log tail
#   results/merlin_tests/stage_<name>.done           per-stage guard (delete one to re-run that stage)
#   stop: kill $(cat results/MERLIN_T_PID)   (exact pid; chain subshells + heartbeat are killed by the trap)
set -u
S=/home/jovyan/tremorferometry/lfe_extension_wa
cd "$S" || exit 1
R=$S/results; L=$S/logs; G=$R/merlin_tests
mkdir -p "$G" "$L"
exec 9>"$R/.merlin_t.lock"
if ! flock -n 9; then echo "another merlin_tests_queue holds the lock"; exit 1; fi
echo $$ > "$R/MERLIN_T_PID"
rm -f "$R/MERLIN_T_DONE" "$R/MERLIN_T_FAILED" "$G/T1_FAILED" "$G/T2_FAILED"
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
echo start > "$G/.stage_T1"; echo start > "$G/.stage_T2"
( while true; do echo "$(date -u +%FT%TZ) T1:$(cat $G/.stage_T1 2>/dev/null) T2:$(cat $G/.stage_T2 2>/dev/null)" > "$R/MERLIN_T_HEARTBEAT"; sleep 60; done ) &
HB=$!
C1=; C2=
trap 'kill $HB 2>/dev/null; [ -n "$C1" ] && kill $C1 2>/dev/null; [ -n "$C2" ] && kill $C2 2>/dev/null' EXIT

run() {  # run <chain> <name> <cmd...>   (returns 1 on failure; the caller stops ITS chain only)
  local ch=$1 name=$2; shift 2
  if [ -f "$G/stage_$name.done" ]; then echo "[skip] $name"; return 0; fi
  echo "$name" > "$G/.stage_$ch"
  echo "$(date -u +%FT%TZ) [stage] $ch $name: $*" | tee -a "$L/merlin_tests.log"
  if "$@" >> "$L/merlin_$name.log" 2>&1; then
    date -u +%FT%TZ > "$G/stage_$name.done"
    echo "$(date -u +%FT%TZ) [ok] $name" | tee -a "$L/merlin_tests.log"
  else
    { echo "chain $ch stage $name failed $(date -u +%FT%TZ)"; tail -40 "$L/merlin_$name.log"; } > "$G/${ch}_FAILED"
    echo "$(date -u +%FT%TZ) [FAILED] $ch $name" | tee -a "$L/merlin_tests.log"
    return 1
  fi
}

chain_t1() {
  export WA_SETUP=$G/t1_setup
  W="V2008,V2010,A2011"
  run T1 t1_setup     python scripts/merlin_t1_setup.py || return 1
  run T1 t1_mf_before python scripts/wa_mf.py --tag mt1_before --windows $W --tpl build --n-control 2 --drop-dead 0.05 || return 1
  run T1 t1_mf_all    python scripts/wa_mf.py --tag mt1_all --windows $W --tpl build --n-control 2 --drop-dead 0.05 --all-templates || return 1
  run T1 t1_hf_before python scripts/wa_hf.py --mf mt1_before --tpl build --workers 8 || return 1
  run T1 t1_hf_all    python scripts/wa_hf.py --mf mt1_all --tpl build --workers 8 --all-templates || return 1
  run T1 t1_sc_val    python scripts/merlin_t1_score.py --mf val --tplmode clear || return 1
  run T1 t1_sc_before python scripts/merlin_t1_score.py --mf mt1_before --tplmode clear || return 1
  run T1 t1_sc_all    python scripts/merlin_t1_score.py --mf mt1_all --tplmode all || return 1
  echo done > "$G/.stage_T1"
}

chain_t2() {
  unset WA_SETUP                     # T2 uses the trial setup (42 stations, V2008/V2010 truth)
  run T2 t2_setup python scripts/merlin_t2.py --stage setup --cap 4 --workers 8 || return 1
  run T2 t2_r0    python scripts/merlin_t2.py --stage round --round 0 || return 1
  run T2 t2_r1    python scripts/merlin_t2.py --stage round --round 1 || return 1
  run T2 t2_r2    python scripts/merlin_t2.py --stage round --round 2 || return 1
  run T2 t2_score python scripts/merlin_t2.py --stage score || return 1
  echo done > "$G/.stage_T2"
}

chain_t1 & C1=$!
chain_t2 & C2=$!
wait $C1; wait $C2
C1=; C2=
export WA_SETUP=$G/t1_setup
run R report python scripts/merlin_tests_report.py
rm -f "$G/stage_report.done"          # the report is always re-written on a resume
if [ -f "$G/T1_FAILED" ] || [ -f "$G/T2_FAILED" ] || [ -f "$G/R_FAILED" ]; then
  cat "$G"/T1_FAILED "$G"/T2_FAILED "$G"/R_FAILED 2>/dev/null > "$R/MERLIN_T_FAILED"
  echo "$(date -u +%FT%TZ) [FAILED] see $R/MERLIN_T_FAILED" | tee -a "$L/merlin_tests.log"
  exit 1
fi
date -u +%FT%TZ > "$R/MERLIN_T_DONE"
echo "$(date -u +%FT%TZ) [DONE]" | tee -a "$L/merlin_tests.log"
