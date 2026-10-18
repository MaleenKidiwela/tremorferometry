#!/usr/bin/env bash
# WA FIXED-10 re-run queue (2026-10-07; design + every deviation: results/wa_fixed10/PLAN.md).  Launch-and-exit:
#   nohup setsid -f bash scripts/wa_f10_queue.sh > results/wa_fixed10/logs/queue.log 2>&1
# Production code, run with WA_PROD_DIR=results/wa_fixed10 (own setup / mf / hf / conf / post); the production caches,
# templates and raw traces are only READ; nothing is written under results/wa_prod/; nothing deleted; no git.
# Markers OUT/WA_F10_{PID,HEARTBEAT,DONE,FAILED}; lock locks/wa_f10.lock; stage markers OUT/<STAGE>_DONE (resume = relaunch).
# Smoke: F10_OUT=<scratch> F10_TRUTH_DAYS=<json> F10_PROD_DAYS=<json> F10_NSHARD=1 WA_PROD_FIT_MAX=.. F10_LOCK=locks/wa_f10_smoke.lock
set -u -o pipefail
WS=/home/jovyan/tremorferometry/lfe_extension_wa
cd "$WS" || exit 1
export WA_PROD_DIR=${F10_OUT:-$WS/results/wa_fixed10}
O=$WA_PROD_DIR
export WA_PROD_LOG=$O/logs/wa_f10.log
export WA_PROD_CUT_DIR=$O/no_cut                      # never used (no cut stage); keeps the production cut path out of reach
NSH=${F10_NSHARD:-2}
LOCKF=${F10_LOCK:-locks/wa_f10.lock}
mkdir -p "$O" "$O/logs" locks
LG=$O/logs
exec 9>$LOCKF
flock -n 9 || { echo "[f10] another instance holds $LOCKF; exiting"; exit 0; }
echo $$ > $O/WA_F10_PID
STAGE="init"; FINISHED=0; KIDS=""; HB=""
beat(){ echo "stage=$STAGE pid=$$ $(date -u +%FT%TZ) gpu=$(nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader 2>/dev/null | head -1)" > $O/WA_F10_HEARTBEAT; }
hb_start(){ ( exec 9>&-; while sleep 120; do
                g=$(nvidia-smi --query-gpu=memory.used,memory.total,utilization.gpu --format=csv,noheader 2>/dev/null | head -1)
                echo "$(date -u +%FT%TZ) $g $(ls $O/mf/prod/run_*.json 2>/dev/null | wc -l) prod-days $(ls $O/hf/prod/hf_*.parquet 2>/dev/null | wc -l) hf-days" >> $O/gpu_mem.log
                [ -f $O/WA_F10_HEARTBEAT ] && touch $O/WA_F10_HEARTBEAT
              done ) & HB=$!; }
killtree(){ local c; for c in $(pgrep -P "$1" 2>/dev/null); do killtree "$c"; done; kill "$1" 2>/dev/null; }
cleanup(){ [ -n "$HB" ] && kill "$HB" 2>/dev/null; for p in $KIDS; do killtree "$p"; done
           for p in $(pgrep -P $$ 2>/dev/null); do killtree "$p"; done; }
fail(){ echo "[f10] FAILED stage=$STAGE rc=$1 $(date -u +%FT%TZ)"; echo "stage=$STAGE rc=$1 $(date -u +%FT%TZ)" > $O/WA_F10_FAILED; cleanup; exit "$1"; }
onexit(){ rc=$?; cleanup; [ "$FINISHED" = 1 ] && return 0; [ -f $O/WA_F10_FAILED ] && return 0
          echo "stage=$STAGE rc=$rc unexpected-exit $(date -u +%FT%TZ)" > $O/WA_F10_FAILED; }
trap 'fail 143' TERM INT HUP
trap onexit EXIT
run(){ STAGE="$1"; shift; beat; echo "[f10] === $STAGE $(date -u +%FT%TZ)"; "$@" 9>&- || fail $?; }
done_or(){ [ -f "$O/$1" ]; }
mark(){ date -u +%FT%TZ > "$O/$1"; }
bg(){ local n=$1; shift
  if done_or BG_${n}_DONE; then echo "[f10] bg $n already done"; return 0; fi
  rm -f $O/BG_${n}_FAILED
  ( exec 9>&-; set -o pipefail; "$@" > $LG/bg_$n.log 2>&1 && date -u +%FT%TZ > $O/BG_${n}_DONE || { echo "rc=$? $(date -u +%FT%TZ)" > $O/BG_${n}_FAILED; } ) &
  KIDS="$KIDS $!"; echo "[f10] bg $n launched pid $! ($LG/bg_$n.log)"; }
bgwait(){ STAGE="wait_$1"; beat
  while ! done_or BG_$1_DONE; do
    [ -f $O/BG_$1_FAILED ] && { tail -15 $LG/bg_$1.log; fail 7; }
    sleep 60
  done; }
rm -f $O/WA_F10_FAILED $O/WA_F10_DONE
hb_start
export PYTHONPATH=$WS/scripts:${PYTHONPATH:-}
export PYTHONUNBUFFERED=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONWARNINGS=ignore
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
PY=python
echo "[f10] start $(date -u +%FT%TZ) OUT=$O nshard=$NSH"
gpu_wait(){ STAGE="gpu_wait"; beat
  while true; do
    held=0
    for f in /home/jovyan/tremorferometry/*/locks/*.lock /home/jovyan/tremorferometry/lfe_extension_wa/results/wa_pre/gpu.lock; do
      [ -e "$f" ] || continue
      case "$f" in *wa_f10*.lock) continue;; esac
      flock -n "$f" true 2>/dev/null || { held=1; echo "$(date -u +%FT%TZ) waiting for GPU lock $f"; }
    done
    busy=$(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader,nounits 2>/dev/null | awk -F', ' '$2>4096' | wc -l)
    [ "$held" = 0 ] && [ "$busy" -le 0 ] && break
    sleep 120
  done; }
TDAYS=${F10_TRUTH_DAYS:-$O/setup/days_truth.json}
PDAYS=${F10_PROD_DAYS:-$O/setup/days_prod.json}

# ================================================================== F0 setup (copies of the production setup; 3 restricted tables)
done_or SETUP_DONE || { run setup $PY scripts/wa_f10_setup.py; mark SETUP_DONE; }
# ================================================================== F1 matched filter: production pass (bg, shards) + truth run
mkdir -p $O/mf/prod
if ! done_or PRODMF_DONE || ! done_or TRUTH_MF_DONE; then gpu_wait; fi
prodmf(){ local bad=0; rm -f $O/gpujob_prodmf_*.rc; : > $O/gpu_prod.txt
  for i in $(seq 0 $((NSH-1))); do
    echo "$PY scripts/wa_prod_mf.py --run prod --days-file $PDAYS --nets core,all --series data,c1,rev --shard $i --nshard $NSH" >> $O/gpu_prod.txt
  done
  nl -ba -w1 -s' ' $O/gpu_prod.txt | xargs -P $NSH -d '\n' -n 1 bash -c 'i=${0%% *}; cmd=${0#* }; ( eval "$cmd" ) > '"$LG"'/prodmf_$i.log 2>&1; echo $? > '"$O"'/gpujob_prodmf_$i.rc'
  for r in $O/gpujob_prodmf_*.rc; do [ "$(cat $r)" = "0" ] || bad=1; done
  [ $bad = 0 ] && $PY -c "import wa_prod_common as P; P.guard('mf_prod', 'mf_prod_s*')" && date -u +%FT%TZ > $O/mf/prod/MF_DONE; }
if done_or PRODMF_DONE; then echo "[f10] prod mf already done"; else bg prodmf prodmf; fi
if ! done_or TRUTH_MF_DONE; then
  run truthmf $PY scripts/wa_prod_mf.py --run truth --days-file $TDAYS --nets core,all --series data,c1,c2,rev
  run guard_mf_truth $PY -c "import wa_prod_common as P; P.guard('mf_truth', 'mf_truth_s*')"
  mark TRUTH_MF_DONE
fi
# HF ratio for every production candidate, as the day files appear (CPU)
bg hfprod $PY scripts/wa_prod_hf.py --run prod --workers ${F10_HF_WORKERS:-16} --wait-mf $O/mf/prod/MF_DONE
# ================================================================== F2 confidence: HF on truth, refit (WAVE recipe), sealed one look
done_or HF_TRUTH_DONE || { run hf_truth $PY scripts/wa_prod_hf.py --run truth --workers 8
                           run guard_hf_truth $PY -c "import wa_prod_common as P; P.guard('hf', 'hf_truth*')"; mark HF_TRUTH_DONE; }
export OMP_NUM_THREADS=12
done_or FIT_DONE || { run fit $PY scripts/wa_prod_fit.py fit; mark FIT_DONE; }
done_or SEALED_DONE || { run open_sealed $PY scripts/wa_prod_fit.py open_sealed; mark SEALED_DONE; }
export OMP_NUM_THREADS=2
# ================================================================== F3 apply to every production candidate
if ! done_or PRODMF_DONE; then bgwait prodmf; mark PRODMF_DONE; fi
bgwait hfprod
run guard_hf_prod $PY -c "import wa_prod_common as P; P.guard('hf', 'hf_prod')"
apply_par(){ local n=${F10_APPLY_PROCS:-4} bad=0 i pids=""
  for i in $(seq 0 $((n-1))); do
    OMP_NUM_THREADS=6 $PY scripts/wa_prod_fit.py apply --run prod --shard $i --nshard $n > $LG/apply_prod_$i.log 2>&1 9>&- &
    pids="$pids $!"
  done
  for i in $pids; do wait $i || bad=1; done
  return $bad; }
done_or APPLY_PROD_DONE || { run apply_prod apply_par; mark APPLY_PROD_DONE; }
# ================================================================== F4 catalogue + rates (post_v2 logic incl. the k and zero-day fixes), both networks
export WA_POST_HEAD_LAB=own WA_POST_SCANNED=famday
for NETX in core all; do
  case $NETX in core) PD=$O/post_v2_core9;; all) PD=$O/post_v2_all10;; esac
  if ! done_or POST_${NETX}_DONE; then
    STAGE=post_$NETX; beat; echo "[f10] === post $NETX -> $PD $(date -u +%FT%TZ)"
    WA_POST_NET=$NETX WA_POST_KNET=$NETX WA_POST_V2_DIR=$PD $PY scripts/wa_prod_post_v2.py all > $LG/post_$NETX.log 2>&1 9>&- || { tail -20 $LG/post_$NETX.log; fail 8; }
    WA_POST_NET=$NETX WA_POST_KNET=$NETX WA_POST_V2_DIR=$PD $PY scripts/wa_f10_strict.py > $LG/strict_$NETX.log 2>&1 9>&- || { tail -20 $LG/strict_$NETX.log; fail 9; }
    mark POST_${NETX}_DONE
  fi
done
date -u +%FT%TZ > $O/WA_F10_DONE
STAGE="done"; FINISHED=1; beat
echo "[f10] DONE $(date -u +%FT%TZ)"
