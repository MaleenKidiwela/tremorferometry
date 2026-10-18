#!/usr/bin/env bash
# SH_CONF (2026-10-01): per-repeat confidence for Shelly (2025)'s 27 Mendocino families, fitted against his catalogue, applied to
# every candidate of the forward run 2023-12 .. 2026-09.  Everything that needs continuous data is computed here and saved.
# Stages (guarded, resumable): days -> mf (GPU, 5 shards: feature re-run + site-rolled scrambled control) -> [hf (CPU) || pnet labels ->
# fine-tune -> pnet feat (GPU, 3 shards)] -> fit -> report.  Each data-reading stage fails if > 3 % of its station-days errored.
# Markers results/SH_CONF_{DONE,FAILED,HEARTBEAT}; pid results/SH_CONF_PID; lock locks/sh_conf.lock; log logs/sh_conf_queue.log
# Smoke: SH_CONF_DIR=<scratch> SH_CONF_DAYS=d1,d2 SH_CONF_BDAYS=b1,b2 SH_CONF_FT_STEPS=200 SH_CONF_TAG=SMOKE bash scripts/sh_conf_queue.sh
set -u -o pipefail
WS=/home/jovyan/tremorferometry/lfe_south
cd "$WS" || exit 1
TAGM=${SH_CONF_TAG:-SH_CONF}
export SH_CONF_DIR=${SH_CONF_DIR:-$WS/results/sh_conf}
O=$SH_CONF_DIR
MK=results; [ "$TAGM" != "SH_CONF" ] && MK=$O
mkdir -p "$O" logs locks
exec 9>locks/${TAGM,,}.lock
flock -n 9 || { echo "[conf] another instance holds the lock; exiting"; exit 0; }
echo $$ > $MK/${TAGM}_PID
STAGE="init"; FINISHED=0; HB=""
beat(){ echo "stage=$STAGE pid=$$ $(date -u +%FT%TZ)" > $MK/${TAGM}_HEARTBEAT; }
hb_start(){ ( exec 9>&-; while sleep 300; do [ -f $MK/${TAGM}_HEARTBEAT ] && touch $MK/${TAGM}_HEARTBEAT; done ) & HB=$!; }
cleanup(){ [ -n "$HB" ] && kill "$HB" 2>/dev/null; }
fail(){ echo "[conf] FAILED stage=$STAGE rc=$1 $(date -u +%FT%TZ)"; echo "stage=$STAGE rc=$1 $(date -u +%FT%TZ)" > $MK/${TAGM}_FAILED; cleanup; exit "$1"; }
onexit(){ rc=$?; cleanup; [ "$FINISHED" = 1 ] && return 0; [ -f $MK/${TAGM}_FAILED ] && return 0
          echo "[conf] EXITED unexpectedly stage=$STAGE rc=$rc $(date -u +%FT%TZ)"
          echo "stage=$STAGE rc=$rc unexpected-exit $(date -u +%FT%TZ)" > $MK/${TAGM}_FAILED; }
trap 'fail 143' TERM INT HUP
trap onexit EXIT
run(){ STAGE="$1"; shift; beat; echo "[conf] === $STAGE $(date -u +%FT%TZ)"; "$@" || fail $?; }
rm -f $MK/${TAGM}_FAILED $MK/${TAGM}_DONE
hb_start
export PYTHONUNBUFFERED=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
CUPYLIB="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')"
DAYS=${SH_CONF_DAYS:+--days $SH_CONF_DAYS}
BDAYS=${SH_CONF_BDAYS:+--days $SH_CONF_BDAYS}
FTSTEPS=${SH_CONF_FT_STEPS:-6000}
NMF=${SH_CONF_NMF:-5}; NPN=${SH_CONF_NPN:-3}
echo "[conf] start $(date -u +%FT%TZ) OUT=$O days=${SH_CONF_DAYS:-all}"; df -h /home/jovyan | tail -1

# ---- 1. day table
[ -f $O/day_table.csv ] || run days python scripts/sh_conf_common.py days

# ---- 2. MF features + scrambled control (GPU, NMF shards in parallel)
if [ ! -f $O/MF_DONE ]; then
  STAGE=mf; beat; echo "[conf] === mf ($NMF shards) $(date -u +%FT%TZ)"
  pids=()
  for i in $(seq 0 $((NMF-1))); do
    LD_LIBRARY_PATH="$CUPYLIB${LD_LIBRARY_PATH:-}" python scripts/sh_conf_mf.py --shard $i --nshard $NMF $DAYS > logs/sh_conf_mf_$i.log 2>&1 &
    pids+=($!)
  done
  rcs=0; for p in "${pids[@]}"; do wait $p || rcs=$?; done
  [ $rcs = 0 ] || fail $rcs
  run mf_guard python scripts/sh_conf_common.py guard mf
  run mf_check python scripts/sh_conf_common.py mfcheck ${SH_CONF_DAYS:-}
  date -u +%FT%TZ > $O/MF_DONE
fi

# ---- 3. HF (CPU, background) || PNET labels -> fine-tune -> features (GPU)
HFPID=""
if [ ! -f $O/HF_DONE ]; then
  ( python scripts/sh_conf_hf.py --workers 6 $DAYS > logs/sh_conf_hf.log 2>&1 && date -u +%FT%TZ > $O/HF_DONE ) &
  HFPID=$!
  echo "[conf] hf launched in background pid $HFPID (log logs/sh_conf_hf.log)"
fi
if [ ! -f $O/PNET_FT_DONE ]; then
  STAGE=pnet_finetune; beat; echo "[conf] === pnet labels + fine-tune $(date -u +%FT%TZ)"
  if python scripts/sh_conf_pnet.py labels --workers 6 $BDAYS && python scripts/sh_conf_pnet.py finetune --steps $FTSTEPS; then
    date -u +%FT%TZ > $O/PNET_FT_DONE
  else
    echo "[conf] WARN PNET fine-tune failed -- continuing with the off-the-shelf model only"; echo failed > $O/PNET_FT_DONE
    [ -f $O/pnet/models/ft1.pt ] && mv $O/pnet/models/ft1.pt $O/pnet/models/ft1_failed.pt
  fi
fi
if [ ! -f $O/PNET_DONE ]; then
  STAGE=pnet_feat; beat; echo "[conf] === pnet feat ($NPN shards) $(date -u +%FT%TZ)"
  pids=()
  for i in $(seq 0 $((NPN-1))); do
    python scripts/sh_conf_pnet.py feat --shard $i --nshard $NPN $DAYS > logs/sh_conf_pnet_$i.log 2>&1 &
    pids+=($!)
  done
  rcs=0; for p in "${pids[@]}"; do wait $p || rcs=$?; done
  [ $rcs = 0 ] || fail $rcs
  run pnet_guard python scripts/sh_conf_common.py guard pnet
  date -u +%FT%TZ > $O/PNET_DONE
fi
if [ -n "$HFPID" ]; then
  STAGE=hf_wait; beat; echo "[conf] === waiting for hf $(date -u +%FT%TZ)"
  wait $HFPID; rc=$?
  [ $rc = 0 ] || { tail -20 logs/sh_conf_hf.log; fail $rc; }
fi
[ -f $O/HF_DONE ] || fail 6

# ---- 4. fit + evaluation + application, report
export OMP_NUM_THREADS=8
[ -f $O/FIT_DONE ] || { run fit python scripts/sh_conf_fit.py; date -u +%FT%TZ > $O/FIT_DONE; }
run report python scripts/sh_conf_report.py
date -u +%FT%TZ > $MK/${TAGM}_DONE
STAGE="done"; FINISHED=1; beat
echo "[conf] DONE $(date -u +%FT%TZ)"
