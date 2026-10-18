#!/usr/bin/env bash
# PNET as independent per-repeat evidence for the WA extension (Merlin step 6, 2026-09-30).
# Stages (guarded, resumable; results/pnet/stage_<name>.done; a failed stage does NOT stop independent later stages):
#   A  cf_w05 (GPU) -> featA -> aucA -> refitA (step-1 lfdr + off-the-shelf PNET features)
#   B  labels -> finetune (GPU) -> cf_ft (GPU) -> evalB (both CFs)
#   C  featB -> aucB -> refitB (step-1 lfdr + fine-tuned PNET features)
#   D  t2_w05 (reference, always) ; t2_ft (only if stageB_verdict.json says works) -> score each
#   report -> results/PNET_RESULTS.md
#   nohup setsid -f bash scripts/pnet_queue.sh > logs/pnet_queue.out 2>&1 < /dev/null
#   results/PNET_PID  results/PNET_HEARTBEAT ("<utc> <stage>" every 60 s)  results/PNET_DONE / results/PNET_FAILED
#   stop: kill $(cat results/PNET_PID)   (exact pid)
set -u
S=/home/jovyan/tremorferometry/lfe_extension_wa
cd "$S" || exit 1
R=$S/results; L=$S/logs; G=$R/pnet
mkdir -p "$G" "$L"
exec 9>"$R/.pnet.lock"
if ! flock -n 9; then echo "another pnet_queue holds the lock"; exit 1; fi
echo $$ > "$R/PNET_PID"
rm -f "$R/PNET_DONE" "$R/PNET_FAILED" "$G/failed_stages.txt"
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
export WA_SETUP=$R/merlin_tests/t1_setup
echo start > "$G/.stage"
( while true; do echo "$(date -u +%FT%TZ) $(cat $G/.stage 2>/dev/null)" > "$R/PNET_HEARTBEAT"; sleep 60; done ) &
HB=$!
trap 'kill $HB 2>/dev/null' EXIT
done_() { [ -f "$G/stage_$1.done" ]; }
run() {   # run <name> <deps (comma list or -)> cmd...
  local name=$1 deps=$2; shift 2
  if done_ "$name"; then echo "[skip] $name"; return 0; fi
  if [ "$deps" != "-" ]; then
    for d in ${deps//,/ }; do
      if ! done_ "$d"; then echo "$(date -u +%FT%TZ) [blocked] $name (needs $d)" | tee -a "$L/pnet.log"; echo "$name blocked by $d" >> "$G/failed_stages.txt"; return 1; fi
    done
  fi
  echo "$name" > "$G/.stage"
  echo "$(date -u +%FT%TZ) [stage] $name: $*" | tee -a "$L/pnet.log"
  if "$@" >> "$L/pnet_$name.log" 2>&1; then
    date -u +%FT%TZ > "$G/stage_$name.done"; echo "$(date -u +%FT%TZ) [ok] $name" | tee -a "$L/pnet.log"
  else
    { echo "stage $name failed $(date -u +%FT%TZ)"; tail -30 "$L/pnet_$name.log"; } >> "$G/failed_stages.txt"
    echo "$(date -u +%FT%TZ) [FAILED] $name" | tee -a "$L/pnet.log"; return 1
  fi
}
W05=/home/jovyan/tremorferometry/lfe_picker/results/p1/models/pnet_w05.pt
FT=$G/models/ft1.pt
EXTRA="pnet_max,pnet_2nd,n_sites_pnet03,n_sites_pnet05"
# ---- A: off-the-shelf
run cf_w05  -          python scripts/pnet_cf.py --model $W05 --tag w05 --windows V2008,V2010,A2011
run featA   cf_w05     python scripts/pnet_feat.py --cf w05 --workers 12
run aucA    featA      python scripts/pnet_auc.py --tag w05
# ---- B: fine-tune (independent of A's features)
run labels  -          python scripts/pnet_labels.py --workers 16
run finetune labels    python scripts/pnet_finetune.py --tag ft1 --steps 20000 --val-every 1000
run cf_ft   finetune   python scripts/pnet_cf.py --model $FT --tag ft1 --windows V2008,V2010,A2011
run evalB   cf_w05,cf_ft python scripts/pnet_eval.py --tags w05,ft1 --workers 16
# ---- D (reference with off-the-shelf CF; GPU-free)
run t2_w05  cf_w05     python scripts/pnet_t2.py --cf w05 --od $G/t2_w05 --workers 8
run t2s_w05 t2_w05     env MT2_OD=$G/t2_w05 python scripts/merlin_t2.py --stage score
# ---- C: step-1 refits (CPU heavy; after the GPU stages so they overlap less)
run refitA  featA      env OMP_NUM_THREADS=16 STEP1_FEATDIR=$G/feat_w05 STEP1_OUT=$G/step1_w05 STEP1_EXTRA_FEATS=$EXTRA python scripts/step1_fit.py
run repA    refitA     env STEP1_OUT=$G/step1_w05 python scripts/step1_report.py
run featB   cf_ft      python scripts/pnet_feat.py --cf ft1 --workers 12
run aucB    featB      python scripts/pnet_auc.py --tag ft1
run refitB  featB      env OMP_NUM_THREADS=16 STEP1_FEATDIR=$G/feat_ft1 STEP1_OUT=$G/step1_ft1 STEP1_EXTRA_FEATS=$EXTRA python scripts/step1_fit.py
run repB    refitB     env STEP1_OUT=$G/step1_ft1 python scripts/step1_report.py
# ---- D with the fine-tuned CF, only if B works
if done_ evalB && python -c "import json,sys; sys.exit(0 if json.load(open('$G/stageB_verdict.json')).get('works') else 1)"; then
  run t2_ft1  cf_ft    python scripts/pnet_t2.py --cf ft1 --od $G/t2_ft1 --workers 8
  run t2s_ft1 t2_ft1   env MT2_OD=$G/t2_ft1 python scripts/merlin_t2.py --stage score
else
  echo "$(date -u +%FT%TZ) [skip] t2_ft1: stage B verdict 'works' is false or evalB missing" | tee -a "$L/pnet.log"
fi
rm -f "$G/stage_report.done"
run report  -          python scripts/pnet_report.py
if [ -s "$G/failed_stages.txt" ]; then
  cp "$G/failed_stages.txt" "$R/PNET_FAILED"; echo "$(date -u +%FT%TZ) [FINISHED WITH FAILURES]" | tee -a "$L/pnet.log"
else
  date -u +%FT%TZ > "$R/PNET_DONE"; echo "$(date -u +%FT%TZ) [DONE]" | tee -a "$L/pnet.log"
fi
