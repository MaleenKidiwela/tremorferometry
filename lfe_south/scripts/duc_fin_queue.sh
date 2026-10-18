#!/usr/bin/env bash
# DUC_FIN (2026-10-01): finish the Ducellier & Creager (2022) extension on the data already on disk (2008-01 .. 2012-07).
# No downloads, nothing deleted (step 5 writes a delete-CANDIDATE list only).  Launch-and-exit; guarded, resumable stages:
#   0 day table + arrival offsets            1 2011 diagnosis part A (results/duc_fin/D2011.md; decision: his data gap)
#   2 re-template the new epochs on disk     3 feature MF + scrambled control (GPU, NMF shards) -> station-day guard
#   4 mfcheck vs duc_fix mf1 + final catalogue + re-joins   + D2011 part B (scrambled control, his station set)
#   5 in parallel: [cut windows -> read guard -> verify 20 windows -> stack check] || [HF ratio from raw, guard]
#                  || [PNET labels -> fine-tune -> PNET features (GPU, 3 shards), guard]
#   6 confidence fit / evaluation / application   7 delete-candidate list   8 report results/DUC_FIN_RESULTS.md
# Every data-reading stage FAILS (marker + exit != 0) if > 3 % of its station-days errored or are silently missing.
# Markers results/DUC_FIN_{HEARTBEAT,DONE,FAILED,PID}; lock locks/duc_fin.lock (+ locks/duc_fix.lock: no DUC_FIX run in parallel);
# log logs/duc_fin_queue.log; GPU memory every 2 min -> results/duc_fin/gpu_mem.log.
# Smoke: DUC_FIN_TAG=SMOKE DUC_FIN_DIR=<scratch> DUC_FIN_CUT_DIR=<scratch>/cut DUC_FIN_REPORT=<scratch>/R.md DUC_FIN_DAYS=d1,d2,..
#        DUC_FIN_RETPL_ARGS="--epochs LDB#2" DUC_FIN_CUT_ARGS="--months 2009-06 --max-days 12" DUC_FIN_FT_STEPS=200 bash scripts/duc_fin_queue.sh
set -u -o pipefail
WS=/home/jovyan/tremorferometry/lfe_south
cd "$WS" || exit 1
TAGM=${DUC_FIN_TAG:-DUC_FIN}
export DUC_FIN_DIR=${DUC_FIN_DIR:-$WS/results/duc_fin}
export DUC_FIN_CUT_DIR=${DUC_FIN_CUT_DIR:-$WS/results/cut/duc}
export DUC_FIN_REPORT=${DUC_FIN_REPORT:-$WS/results/DUC_FIN_RESULTS.md}
O=$DUC_FIN_DIR
MK=results; [ "$TAGM" != "DUC_FIN" ] && MK=$O
mkdir -p "$O" "$DUC_FIN_CUT_DIR" logs locks
exec 9>locks/${TAGM,,}.lock
flock -n 9 || { echo "[fin] another instance holds the lock; exiting"; exit 0; }
exec 8>locks/duc_fix.lock
flock -n 8 || { echo "[fin] a DUC_FIX / DUC_DIAG queue holds locks/duc_fix.lock; exiting"; exit 0; }
echo $$ > $MK/${TAGM}_PID
STAGE="init"; FINISHED=0; HB=""; KIDS=""
beat(){ echo "stage=$STAGE pid=$$ $(date -u +%FT%TZ) gpu=$(nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader 2>/dev/null | head -1)" > $MK/${TAGM}_HEARTBEAT; }
hb_start(){ ( exec 9>&- 8>&-; while sleep 120; do
                g=$(nvidia-smi --query-gpu=memory.used,memory.total,utilization.gpu --format=csv,noheader 2>/dev/null | head -1)
                echo "$(date -u +%FT%TZ) $g" >> $O/gpu_mem.log
                [ -f $MK/${TAGM}_HEARTBEAT ] && sed -i "s#gpu=.*#gpu=$g#" $MK/${TAGM}_HEARTBEAT && touch $MK/${TAGM}_HEARTBEAT
              done ) & HB=$!; }
cleanup(){ [ -n "$HB" ] && kill "$HB" 2>/dev/null; for p in $KIDS; do pkill -P "$p" 2>/dev/null; kill "$p" 2>/dev/null; done; }
fail(){ echo "[fin] FAILED stage=$STAGE rc=$1 $(date -u +%FT%TZ)"; echo "stage=$STAGE rc=$1 $(date -u +%FT%TZ)" > $MK/${TAGM}_FAILED; cleanup; exit "$1"; }
onexit(){ rc=$?; cleanup; [ "$FINISHED" = 1 ] && return 0; [ -f $MK/${TAGM}_FAILED ] && return 0
          echo "[fin] EXITED unexpectedly stage=$STAGE rc=$rc $(date -u +%FT%TZ)"
          echo "stage=$STAGE rc=$rc unexpected-exit $(date -u +%FT%TZ)" > $MK/${TAGM}_FAILED; }
trap 'fail 143' TERM INT HUP
trap onexit EXIT
run(){ STAGE="$1"; shift; beat; echo "[fin] === $STAGE $(date -u +%FT%TZ)"; "$@" || fail $?; }
rm -f $MK/${TAGM}_FAILED $MK/${TAGM}_DONE $O/CUT_FAILED $O/HF_FAILED
hb_start
export PYTHONUNBUFFERED=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
export DUC_MUTE=$(python -c "import json; print(json.load(open('results/duc_fix/detector_choice.json'))['mute'])") || fail 2
DAYS=${DUC_FIN_DAYS:+--days $DUC_FIN_DAYS}
NMF=${DUC_FIN_NMF:-5}; NPN=${DUC_FIN_NPN:-3}; FTSTEPS=${DUC_FIN_FT_STEPS:-6000}
echo "[fin] start $(date -u +%FT%TZ) OUT=$O CUT=$DUC_FIN_CUT_DIR days=${DUC_FIN_DAYS:-all} DUC_MUTE=$DUC_MUTE"; df -h /home/jovyan | tail -1
PY=python

# ---- 0. day table, arrival offsets
[ -f $O/day_table.csv ] || run days $PY scripts/duc_fin_common.py days
[ -f $O/arrival_offsets.csv ] || run offsets $PY scripts/duc_fin_common.py offsets
# ---- 1. 2011 diagnosis, part A (catalogue files)
[ -f $O/D2011_A_DONE ] || { run d2011_events $PY scripts/duc_fin_d2011.py events; date -u +%FT%TZ > $O/D2011_A_DONE; }
# ---- 2. re-templating of the new epochs that start inside the data on disk
[ -f $O/RETPL_DONE ] || { run retpl $PY scripts/duc_fin_retpl.py --workers 6 ${DUC_FIN_RETPL_ARGS:-}; date -u +%FT%TZ > $O/RETPL_DONE; }
# ---- 3. feature MF + scrambled control (GPU shards)
if [ ! -f $O/MF_DONE ]; then
  STAGE=mf; beat; echo "[fin] === mf ($NMF shards) $(date -u +%FT%TZ)"
  pids=()
  for i in $(seq 0 $((NMF-1))); do
    $PY scripts/duc_fin_mf.py --shard $i --nshard $NMF $DAYS > logs/duc_fin_mf_$i.log 2>&1 &
    pids+=($!); KIDS="$KIDS $!"
  done
  rcs=0; for p in "${pids[@]}"; do wait $p || rcs=$?; done
  [ $rcs = 0 ] || { tail -5 logs/duc_fin_mf_*.log; fail $rcs; }
  run mf_guard $PY scripts/duc_fin_common.py guard mf
  date -u +%FT%TZ > $O/MF_DONE
fi
# ---- 4. mfcheck + final catalogue + joins; D2011 part B
[ -f $O/JOIN_DONE ] || { run join $PY scripts/duc_fin_join.py $DAYS; date -u +%FT%TZ > $O/JOIN_DONE; }
[ -f $O/D2011_B_DONE ] || { run d2011_control $PY scripts/duc_fin_d2011.py control; date -u +%FT%TZ > $O/D2011_B_DONE; }
# ---- 5. parallel: cut || HF || PNET
CPID=""; HPID=""
if [ ! -f $O/CUT_DONE ]; then
  ( set -o pipefail
    { $PY scripts/duc_fin_cut.py cut --workers 6 ${DUC_FIN_CUT_ARGS:-} && $PY scripts/duc_fin_cut.py guard \
      && $PY scripts/duc_fin_cut.py verify --n 20 && $PY scripts/duc_fin_cut.py stackcheck; } > logs/duc_fin_cut.log 2>&1 \
      && date -u +%FT%TZ > $O/CUT_DONE || { echo "rc=$? $(date -u +%FT%TZ)" > $O/CUT_FAILED; } ) &
  CPID=$!; KIDS="$KIDS $CPID"; echo "[fin] cut launched pid $CPID (logs/duc_fin_cut.log)"
fi
if [ ! -f $O/HF_DONE ]; then
  ( $PY scripts/duc_fin_hf.py --workers 10 $DAYS > logs/duc_fin_hf.log 2>&1 && date -u +%FT%TZ > $O/HF_DONE \
      || { echo "rc=$? $(date -u +%FT%TZ)" > $O/HF_FAILED; } ) &
  HPID=$!; KIDS="$KIDS $HPID"; echo "[fin] hf launched pid $HPID (logs/duc_fin_hf.log)"
fi
if [ ! -f $O/PNET_FT_DONE ]; then
  STAGE=pnet_finetune; beat; echo "[fin] === pnet labels + fine-tune $(date -u +%FT%TZ)"
  LDAYS=${DUC_FIN_LDAYS:+--days $DUC_FIN_LDAYS}
  if $PY scripts/duc_fin_pnet.py labels --workers 8 $LDAYS > logs/duc_fin_pnet_labels.log 2>&1; then :; else
    rc=$?; tail -5 logs/duc_fin_pnet_labels.log; fail $rc; fi
  if $PY scripts/duc_fin_pnet.py finetune --steps $FTSTEPS > logs/duc_fin_pnet_ft.log 2>&1; then
    date -u +%FT%TZ > $O/PNET_FT_DONE
  else
    echo "[fin] WARN PNET fine-tune failed -- continuing with the off-the-shelf model only"; echo failed > $O/PNET_FT_DONE
    [ -f $O/pnet/models/ft1.pt ] && mv $O/pnet/models/ft1.pt $O/pnet/models/ft1_failed.pt
  fi
fi
if [ ! -f $O/PNET_DONE ]; then
  STAGE=pnet_feat; beat; echo "[fin] === pnet feat ($NPN shards) $(date -u +%FT%TZ)"
  pids=()
  for i in $(seq 0 $((NPN-1))); do
    $PY scripts/duc_fin_pnet.py feat --shard $i --nshard $NPN $DAYS > logs/duc_fin_pnet_$i.log 2>&1 &
    pids+=($!); KIDS="$KIDS $!"
  done
  rcs=0; for p in "${pids[@]}"; do wait $p || rcs=$?; done
  [ $rcs = 0 ] || fail $rcs
  run pnet_guard $PY scripts/duc_fin_common.py guard pnet
  date -u +%FT%TZ > $O/PNET_DONE
fi
STAGE=wait_cut_hf; beat; echo "[fin] === waiting for cut / hf $(date -u +%FT%TZ)"
[ -n "$HPID" ] && wait $HPID
[ -f $O/HF_FAILED ] && { tail -5 logs/duc_fin_hf.log; STAGE=hf; fail 6; }
[ -n "$CPID" ] && wait $CPID
[ -f $O/CUT_FAILED ] && { tail -8 logs/duc_fin_cut.log; STAGE=cut; fail 6; }
[ -f $O/HF_DONE ] && [ -f $O/CUT_DONE ] || { STAGE=cut_hf_markers; fail 6; }
# ---- 6. confidence
export OMP_NUM_THREADS=8
[ -f $O/FIT_DONE ] || { run fit $PY scripts/duc_fin_fit.py; date -u +%FT%TZ > $O/FIT_DONE; }
# ---- 7. delete-candidate list (NOTHING deleted)
[ -f $O/DELETE_LIST_DONE ] || { run delete_list $PY scripts/duc_fin_delete.py; date -u +%FT%TZ > $O/DELETE_LIST_DONE; }
# ---- 8. report
run report $PY scripts/duc_fin_report.py
date -u +%FT%TZ > $MK/${TAGM}_DONE
STAGE="done"; FINISHED=1; beat
echo "[fin] DONE $(date -u +%FT%TZ)"
