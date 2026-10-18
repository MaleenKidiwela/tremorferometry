#!/usr/bin/env bash
# STEP 2 (Merlin, 2026-09-30): detect <= 3-km GROUPS of families, then attribute which member fired (probability),
# combined with a group-level P(real) (step-1 lfdr refitted on the group candidates).
# Guarded, resumable stages:
#   groups -> tpl -> { A: mf (GPU) -> hf -> attrib_cc (GPU) | B: mf_c34 (GPU, 2 more scrambled controls) }
#          -> { score | feat -> fit } -> attrib_eval -> report
#   nohup setsid -f bash scripts/step2_groups_queue.sh > logs/step2_queue.out 2>&1 < /dev/null
#   results/STEP2_PID  results/STEP2_HEARTBEAT ("<utc> A:<stage> B:<stage> M:<stage>" every 60 s)
#   results/STEP2_DONE / results/STEP2_FAILED (which stage + log tail)
#   results/step2_groups/stage_<name>.done   (delete one to re-run that stage)
#   stop: kill $(cat results/STEP2_PID)   (exact pid; chain subshells + heartbeat are killed by the trap)
set -u
S=/home/jovyan/tremorferometry/lfe_extension_wa
cd "$S" || exit 1
R=$S/results; L=$S/logs; G=$R/step2_groups
mkdir -p "$G" "$L"
exec 9>"$R/.step2.lock"
if ! flock -n 9; then echo "another step2_groups_queue holds the lock"; exit 1; fi
echo $$ > "$R/STEP2_PID"
rm -f "$R/STEP2_DONE" "$R/STEP2_FAILED" "$G"/*_FAILED
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
export WA_SETUP=$G/setup
echo start > "$G/.stage_M"; echo - > "$G/.stage_A"; echo - > "$G/.stage_B"
( while true; do echo "$(date -u +%FT%TZ) A:$(cat $G/.stage_A 2>/dev/null) B:$(cat $G/.stage_B 2>/dev/null) M:$(cat $G/.stage_M 2>/dev/null)" > "$R/STEP2_HEARTBEAT"; sleep 60; done ) &
HB=$!
C1=; C2=
trap 'kill $HB 2>/dev/null; [ -n "$C1" ] && kill $C1 2>/dev/null; [ -n "$C2" ] && kill $C2 2>/dev/null' EXIT

run() {  # run <chain> <name> <cmd...>
  local ch=$1 name=$2; shift 2
  if [ -f "$G/stage_$name.done" ]; then echo "[skip] $name"; return 0; fi
  echo "$name" > "$G/.stage_$ch"
  echo "$(date -u +%FT%TZ) [stage] $ch $name: $*" | tee -a "$L/step2_queue.log"
  if "$@" >> "$L/step2_$name.log" 2>&1; then
    date -u +%FT%TZ > "$G/stage_$name.done"
    echo "$(date -u +%FT%TZ) [ok] $name" | tee -a "$L/step2_queue.log"
  else
    { echo "chain $ch stage $name failed $(date -u +%FT%TZ)"; tail -40 "$L/step2_$name.log"; } > "$G/${ch}_FAILED"
    echo "$(date -u +%FT%TZ) [FAILED] $ch $name" | tee -a "$L/step2_queue.log"
    return 1
  fi
}
fail() { cat "$G"/*_FAILED 2>/dev/null > "$R/STEP2_FAILED"; echo "$(date -u +%FT%TZ) [FAILED] see $R/STEP2_FAILED" | tee -a "$L/step2_queue.log"; exit 1; }
W3="V2008,V2010,A2011"
run M groups python scripts/step2_groups.py --stage groups || fail
run M tpl    python scripts/step2_groups.py --stage tpl --workers 12 || fail
echo "wait A+B" > "$G/.stage_M"
chain_a() {
  run A mf        python scripts/wa_mf.py --tag step2grp --windows $W3 --tpl step2grp --n-control 2 --drop-dead 0.05 --all-templates || return 1
  run A hf        python scripts/wa_hf.py --mf step2grp --tpl step2grp --workers 8 --all-templates || return 1
  run A attrib_cc python scripts/step2_attrib.py --stage cc || return 1
  echo done > "$G/.stage_A"
}
chain_b() {
  run B mf_c34 python scripts/wa_mf.py --tag step2grp_c34 --windows $W3 --tpl step2grp --n-control 2 --drop-dead 0.05 --all-templates --seed 20261001 || return 1
  echo done > "$G/.stage_B"
}
chain_a & C1=$!
chain_b & C2=$!
wait $C1; wait $C2
C1=; C2=
ls "$G"/*_FAILED >/dev/null 2>&1 && fail
# score (1 CPU) beside feat (12 workers)
( run A score python scripts/step2_score.py --mf step2grp ) & C1=$!
export STEP1_MF1=step2grp STEP1_MF2=step2grp_c34 STEP1_TPL=step2grp STEP1_FEATDIR=$G/step1/feat STEP1_OUT=$G/step1
mkdir -p "$G/step1"
run M feat python scripts/step1_features.py --workers 12 --windows $W3 || { wait $C1; fail; }
run M fit  env OMP_NUM_THREADS=16 python scripts/step1_fit.py || { wait $C1; fail; }
wait $C1; C1=
ls "$G"/*_FAILED >/dev/null 2>&1 && fail
run M attrib_eval python scripts/step2_attrib.py --stage eval || fail
run M report python scripts/step2_report.py || fail
rm -f "$G/stage_report.done"          # the report is always re-written on a resume
echo done > "$G/.stage_M"
date -u +%FT%TZ > "$R/STEP2_DONE"
echo "$(date -u +%FT%TZ) [DONE]" | tee -a "$L/step2_queue.log"
