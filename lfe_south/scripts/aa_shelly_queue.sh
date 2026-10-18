#!/usr/bin/env bash
# Agent AA -- re-run of the SHELLY (2025) stage ONLY, after the two detector fixes.
# Never touches the Ducellier reproduction, the 15,400 downloaded station-days, or any path the
# running t_south_queue.sh (pid 1509294) writes: every output goes to results/aa_shelly/.
# t11_mf.py / t9_templates.py / t8_score.py are COPIED (t11_mf_aa.py, t9_templates_aa.py,
# t8_score_aa.py), never edited, because the running queue's extension stage uses t11_mf.py.
set -u -o pipefail
WS=/home/jovyan/tremorferometry/lfe_south
cd "$WS" || exit 1
mkdir -p results/aa_shelly logs locks
exec 9>locks/aa_shelly.lock
flock -n 9 || { echo "[aa] another instance holds the lock; exiting"; exit 0; }
STAGE="init"; FINISHED=0
beat(){ echo "stage=$STAGE pid=$$ $(date -u +%FT%TZ)" > results/AA_SHELLY_HEARTBEAT; }
fail(){ echo "[aa] FAILED stage=$STAGE rc=$1 $(date -u +%FT%TZ)"; echo "stage=$STAGE rc=$1 $(date -u +%FT%TZ)" > results/AA_SHELLY_FAILED; exit "$1"; }
onexit(){ rc=$?; [ "$FINISHED" = 1 ] && return 0; [ -f results/AA_SHELLY_FAILED ] && return 0
          echo "[aa] EXITED unexpectedly stage=$STAGE rc=$rc $(date -u +%FT%TZ)"
          echo "stage=$STAGE rc=$rc unexpected-exit $(date -u +%FT%TZ)" > results/AA_SHELLY_FAILED; }
trap 'fail 143' TERM INT HUP
trap onexit EXIT
run(){ STAGE="$1"; shift; beat; echo "[aa] === $STAGE $(date -u +%FT%TZ)"; "$@" || fail $?; }
soft(){ STAGE="$1"; shift; beat; echo "[aa] === $STAGE (non-fatal) $(date -u +%FT%TZ)"; "$@" || echo "[aa] WARN $STAGE rc=$? -- continuing"; }
rm -f results/AA_SHELLY_FAILED results/AA_SHELLY_DONE
export PYTHONUNBUFFERED=1
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
AA=results/aa_shelly
DAYS=results/stage2_days_test.json
T0=2018-01-01; T1=2024-12-01
SW=6,7,8,9,10,12,15

# ---- ablation first (cheap, 2018 only): normalisation fixed, moveout NOT applied ----------------
[ -f $AA/mfN.parquet ] || soft mfN env PYTHONHASHSEED=0 python scripts/t11_mf_aa.py --source Shelly2025 \
    --templates results/templates_south_build --days-file $DAYS --truth shelly --tmin 6 --gap 4 \
    --align none --t0 2018-01-01 --t1 2019-01-01 --out $AA/mfN.parquet
[ -f $AA/mfN_score.csv ] || soft scoreN env PYTHONHASHSEED=0 python scripts/t8_score_aa.py --ours $AA/mfN.parquet \
    --truth shelly --t0 2018-01-01 --t1 2019-01-01 --days-file $DAYS --days-mode source --sweep $SW \
    --out $AA/mfN_score.csv

# ---- row A: the fix, existing 12 s templates, cache band -----------------------------------------
[ -f $AA/mfA.parquet ] || run mfA env PYTHONHASHSEED=0 python scripts/t11_mf_aa.py --source Shelly2025 \
    --templates results/templates_south_build --days-file $DAYS --truth shelly --tmin 6 --gap 4 \
    --t0 $T0 --t1 $T1 --out $AA/mfA.parquet
[ -f $AA/mfA_score.csv ] || run scoreA env PYTHONHASHSEED=0 python scripts/t8_score_aa.py --ours $AA/mfA.parquet \
    --truth shelly --t0 $T0 --t1 $T1 --days-file $DAYS --days-mode source --sweep $SW --out $AA/mfA_score.csv
[ -f $AA/mfA_score_fam.csv ] || soft scoreA_fam env PYTHONHASHSEED=0 python scripts/t8_score_aa.py --ours $AA/mfA.parquet \
    --truth shelly --t0 $T0 --t1 $T1 --days-file $DAYS --days-mode fam --sweep $SW --out $AA/mfA_score_fam.csv
[ -f $AA/mfA_score_all.csv ] || soft scoreA_all env PYTHONHASHSEED=0 python scripts/t8_score_aa.py --ours $AA/mfA.parquet \
    --truth shelly --t0 $T0 --t1 $T1 --days-mode all --sweep $SW --out $AA/mfA_score_all.csv
soft reportA python scripts/aa_report.py

# ---- row B: Shelly's 2-9 Hz band ----------------------------------------------------------------
[ -f $AA/mfB.parquet ] || run mfB env PYTHONHASHSEED=0 python scripts/t11_mf_aa.py --source Shelly2025 \
    --templates results/templates_south_build --days-file $DAYS --truth shelly --tmin 6 --gap 4 --hp2 \
    --t0 $T0 --t1 $T1 --out $AA/mfB.parquet
[ -f $AA/mfB_score.csv ] || run scoreB env PYTHONHASHSEED=0 python scripts/t8_score_aa.py --ours $AA/mfB.parquet \
    --truth shelly --t0 $T0 --t1 $T1 --days-file $DAYS --days-mode source --sweep $SW --out $AA/mfB_score.csv
soft reportB python scripts/aa_report.py

# ---- row C: Shelly's 6 s template length, 2-9 Hz -------------------------------------------------
[ -f $AA/templates6_build.csv ] || run tpl6 env PYTHONHASHSEED=0 python scripts/t9_templates_aa.py \
    --half build --pre 2.0 --len 6.0 --out $AA/templates6_build.npz
[ -f $AA/mfC.parquet ] || run mfC env PYTHONHASHSEED=0 python scripts/t11_mf_aa.py --source Shelly2025 \
    --templates $AA/templates6_build --days-file $DAYS --truth shelly --tmin 6 --gap 4 --hp2 \
    --t0 $T0 --t1 $T1 --out $AA/mfC.parquet
[ -f $AA/mfC_score.csv ] || run scoreC env PYTHONHASHSEED=0 python scripts/t8_score_aa.py --ours $AA/mfC.parquet \
    --truth shelly --t0 $T0 --t1 $T1 --days-file $DAYS --days-mode source --sweep $SW --out $AA/mfC_score.csv
soft reportC python scripts/aa_report.py

# ---- row D: Shelly's template ITERATION, built on BUILD days only -------------------------------
[ -f $AA/mfA_build.parquet ] || soft mfA_build env PYTHONHASHSEED=0 python scripts/t11_mf_aa.py --source Shelly2025 \
    --templates results/templates_south_build --days-file results/stage2_days_build.json --tmin 8 --gap 4 \
    --t0 $T0 --t1 $T1 --out $AA/mfA_build.parquet
[ -f $AA/templatesD_build.csv ] || soft restackD env PYTHONHASHSEED=0 python scripts/aa_restack.py \
    --det $AA/mfA_build.parquet --templates results/templates_south_build --nstack 200 \
    --out $AA/templatesD_build.npz
[ -f $AA/mfD.parquet ] || soft mfD env PYTHONHASHSEED=0 python scripts/t11_mf_aa.py --source Shelly2025 \
    --templates $AA/templatesD_build --days-file $DAYS --truth shelly --tmin 6 --gap 4 \
    --t0 $T0 --t1 $T1 --out $AA/mfD.parquet
[ -f $AA/mfD_score.csv ] || soft scoreD env PYTHONHASHSEED=0 python scripts/t8_score_aa.py --ours $AA/mfD.parquet \
    --truth shelly --t0 $T0 --t1 $T1 --days-file $DAYS --days-mode source --sweep $SW --out $AA/mfD_score.csv

run report python scripts/aa_report.py
date -u +%FT%TZ > results/AA_SHELLY_DONE
STAGE="done"; FINISHED=1; beat
echo "[aa] DONE $(date -u +%FT%TZ)"
