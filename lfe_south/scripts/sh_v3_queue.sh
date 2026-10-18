#!/usr/bin/env bash
# SH_V3 (2026-09-29): Shelly (2025) Text S1 literal templates -- see results/SH_V3_LOG.md.
# tpl_it0 (his times, own build days, top-100/500, <=0.1 s shifts, best 50, S pick, S-1..S+5) ->
# MF test it0 + Text S1 rescore -> MF build it0 -> tpl_it1 (our build detections >=8xMAD & sum>=2.75)
# -> MF test it1 + rescore -> mix (6 low-SNR families on it1) + rescore -> report.
# Outputs only under results/sh_v3/.  Markers results/SH_V3_{DONE,FAILED,HEARTBEAT}; pid results/SH_V3_PID.
set -u -o pipefail
WS=/home/jovyan/tremorferometry/lfe_south
cd "$WS" || exit 1
mkdir -p results/sh_v3 logs locks
exec 9>locks/sh_v3.lock
flock -n 9 || { echo "[v3] another instance holds the lock; exiting"; exit 0; }
echo $$ > results/SH_V3_PID
STAGE="init"; FINISHED=0; HB=""
beat(){ echo "stage=$STAGE pid=$$ $(date -u +%FT%TZ)" > results/SH_V3_HEARTBEAT; }
hb_start(){ ( exec 9>&-; while sleep 300; do [ -f results/SH_V3_HEARTBEAT ] && touch results/SH_V3_HEARTBEAT; done ) & HB=$!; }
cleanup(){ [ -n "$HB" ] && kill "$HB" 2>/dev/null; }
fail(){ echo "[v3] FAILED stage=$STAGE rc=$1 $(date -u +%FT%TZ)"; echo "stage=$STAGE rc=$1 $(date -u +%FT%TZ)" > results/SH_V3_FAILED; cleanup; exit "$1"; }
onexit(){ rc=$?; cleanup; [ "$FINISHED" = 1 ] && return 0; [ -f results/SH_V3_FAILED ] && return 0
          echo "[v3] EXITED unexpectedly stage=$STAGE rc=$rc $(date -u +%FT%TZ)"
          echo "stage=$STAGE rc=$rc unexpected-exit $(date -u +%FT%TZ)" > results/SH_V3_FAILED; }
trap 'fail 143' TERM INT HUP
trap onexit EXIT
run(){ STAGE="$1"; shift; beat; echo "[v3] === $STAGE $(date -u +%FT%TZ)"; "$@" || fail $?; }
soft(){ STAGE="$1"; shift; beat; echo "[v3] === $STAGE (non-fatal) $(date -u +%FT%TZ)"; "$@" || echo "[v3] WARN $STAGE rc=$? -- continuing"; }
rm -f results/SH_V3_FAILED results/SH_V3_DONE
hb_start
export PYTHONUNBUFFERED=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
V=results/sh_v3
TEST=results/aa_shelly/stage2_days_test.json
BUILD=results/aa_shelly/stage2_days_build.json
T0=2018-01-01; T1=2024-12-01

[ -f $V/tpl_it0.csv ] || run tpl0 python scripts/sh_v3_restack.py --mode his --days-file $BUILD --workers 4 --out $V/tpl_it0.npz
for k in 0 1; do
  [ -f $V/mf_it${k}_test.parquet ] || run mf${k}_test python scripts/t11_mf_v2.py --templates $V/tpl_it$k \
      --days-file $TEST --truth shelly --t0 $T0 --t1 $T1 --out $V/mf_it${k}_test.parquet
  [ -f $V/rescore_it$k.csv ] || run rescore$k python scripts/sh_v3_rescore.py --ours $V/mf_it${k}_test.parquet --tag it$k --out $V/rescore_it$k.csv
  soft report$k python scripts/sh_v3_report.py
  if [ $k -eq 0 ]; then
    [ -f $V/mf_it0_build.parquet ] || run mf0_build python scripts/t11_mf_v2.py --templates $V/tpl_it0 \
        --days-file $BUILD --t0 $T0 --t1 $T1 --out $V/mf_it0_build.parquet
    [ -f $V/tpl_it1.csv ] || run tpl1 python scripts/sh_v3_restack.py --mode det --det $V/mf_it0_build.parquet \
        --picks $V/tpl_it0.csv --days-file $BUILD --thr 8 --workers 4 --out $V/tpl_it1.npz
  fi
done
run mix python scripts/sh_v3_mix.py
[ -f $V/rescore_mix.csv ] || run rescore_mix python scripts/sh_v3_rescore.py --ours $V/mf_mix_test.parquet --tag mix --out $V/rescore_mix.csv
run report python scripts/sh_v3_report.py
date -u +%FT%TZ > results/SH_V3_DONE
STAGE="done"; FINISHED=1; beat
echo "[v3] DONE $(date -u +%FT%TZ)"
