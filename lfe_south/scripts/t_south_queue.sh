#!/usr/bin/env bash
# Agent T — southern LFE REPRODUCTION + existence + extension queue.
# ONE detached queue.  Every stage is guarded by its own output, so a relaunch resumes.
# Order: the reproductions first (they carry a known answer), then existence, then extension.
set -u -o pipefail
WS=/home/jovyan/tremorferometry/lfe_south
cd "$WS" || exit 1
mkdir -p results logs locks traces cache20hz cache_duc resp results/ext
exec 9>locks/t_south.lock
flock -n 9 || { echo "[queue] another instance holds the lock; exiting"; exit 0; }
STAGE="init"; FINISHED=0
beat(){ echo "stage=$STAGE pid=$$ $(date -u +%FT%TZ)" > results/T_SOUTH_HEARTBEAT; }
fail(){ echo "[queue] FAILED at stage=$STAGE rc=$1 $(date -u +%FT%TZ)"; echo "stage=$STAGE rc=$1 $(date -u +%FT%TZ)" > results/T_SOUTH_FAILED; exit "$1"; }
onexit(){ rc=$?; [ "$FINISHED" = 1 ] && return 0; [ -f results/T_SOUTH_FAILED ] && return 0
          echo "[queue] EXITED unexpectedly at stage=$STAGE rc=$rc $(date -u +%FT%TZ)"
          echo "stage=$STAGE rc=$rc unexpected-exit $(date -u +%FT%TZ)" > results/T_SOUTH_FAILED; }
trap 'fail 143' TERM INT HUP
trap onexit EXIT
run(){ STAGE="$1"; shift; beat; echo "[queue] === $STAGE $(date -u +%FT%TZ)"; "$@" || fail $?; }
soft(){ STAGE="$1"; shift; beat; echo "[queue] === $STAGE (non-fatal) $(date -u +%FT%TZ)"; "$@" || echo "[queue] WARN $STAGE rc=$? — continuing"; }
rm -f results/T_SOUTH_FAILED
PY=python
export PYTHONUNBUFFERED=1
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"

# ---------- 0. inventory / geometry / plans (fast, idempotent) ----------
[ -f results/stage1_download_plan.csv ] || run geometry  $PY scripts/t0_geometry.py
[ -f results/stage1_rosters.csv ]       || run rosters   $PY scripts/t1_rosters.py
[ -f results/duc_download_plan.csv ]    || run duc_plan  $PY scripts/t5_duc_plan.py
run inventory $PY scripts/t0b_inventory.py

# ---------- 1. REPRODUCTION of Ducellier & Creager 2022 — his templates, his stations ----------
if [ ! -f results/duc_repro_bedA_score.csv ]; then
  run duc_fetch  $PY scripts/t2_fetch.py results/duc_download_plan.csv --workers 8 \
      --manifest results/duc_download_manifest.csv
  run duc_cache  $PY scripts/t6_duc_cache.py results/duc_download_plan.csv --workers 16
  run duc_mfA    $PY scripts/t7_duc_mf.py --bed A --t0 2008-04-01 --t1 2008-05-10
  run duc_scoreA $PY scripts/t8_score.py --ours results/duc_repro_bedA_detections.parquet \
      --truth ducellier --t0 2008-04-01 --t1 2008-05-11 --out results/duc_repro_bedA_score.csv
fi
date -u +%FT%TZ > results/T_SOUTH_REPRO_A_DONE
soft report_A $PY scripts/t13_report.py
# bed B — his set PLUS every station installed since, each given a template from its own stack
if [ ! -f results/duc_repro_bedB_score.csv ]; then
  soft bedB_tpl  $PY scripts/t7b_bedB_templates.py --t0 2008-04-01 --t1 2008-05-10
  soft duc_mfB   $PY scripts/t7_duc_mf.py --bed B --t0 2008-04-01 --t1 2008-05-10 \
      --extra-templates results/duc_bedB_templates.npz
  soft duc_scoreB $PY scripts/t8_score.py --ours results/duc_repro_bedB_detections.parquet \
      --truth ducellier --t0 2008-04-01 --t1 2008-05-11 --out results/duc_repro_bedB_score.csv
fi
soft report_B $PY scripts/t13_report.py

# ---------- 2. STAGE 1 — existence test on the modern permanent network ----------
if [ ! -f results/stage1_existence_family.csv ]; then
  run s1_fetch  $PY scripts/t2_fetch.py results/stage1_download_plan.csv --workers 8
  run s1_cache  $PY scripts/t3_cache.py results/stage1_download_plan.csv --workers 24
  run s1_test   $PY scripts/t4_stack_test.py --nnull 200
fi
date -u +%FT%TZ > results/T_SOUTH_STAGE1_DONE
soft report_1 $PY scripts/t13_report.py

# ---------- 3. templates from our own stacks (build half, and all) ----------
[ -f results/templates_south_build.csv ] || run tpl_build $PY scripts/t9_templates.py \
    --half build --out results/templates_south_build.npz
[ -f results/templates_south.csv ]       || run tpl_all   $PY scripts/t9_templates.py \
    --half all   --out results/templates_south.npz

# ---------- 4. REPRODUCTION of Shelly 2025 — our templates, his times, HELD-OUT days ----------
if [ ! -f results/shelly_repro_score.csv ]; then
  run shelly_mf $PY scripts/t11_mf.py --source Shelly2025 \
      --templates results/templates_south_build --days-file results/stage2_days_test.json \
      --t0 2018-01-01 --t1 2024-12-01 --out results/shelly_repro_detections.parquet
  run shelly_score $PY scripts/t8_score.py --ours results/shelly_repro_detections.parquet \
      --truth shelly --t0 2018-01-01 --t1 2024-12-01 --out results/shelly_repro_score.csv
fi
date -u +%FT%TZ > results/T_SOUTH_REPRO_DONE
soft report_2 $PY scripts/t13_report.py

# ---------- 5. blind discovery on the southern permanent bed (reproduced Bostock recipe) ----------
soft port_k    $PY scripts/t14_port_k.py
soft discovery bash scripts/t15_discovery.sh

# ---------- 6. EXTENSION to 2026 (long tail; per-month guards) ----------
soft extend bash scripts/t16_extend.sh

run report $PY scripts/t13_report.py
date -u +%FT%TZ > results/T_SOUTH_DONE
STAGE="done"; FINISHED=1; beat
echo "[queue] DONE $(date -u +%FT%TZ)"
