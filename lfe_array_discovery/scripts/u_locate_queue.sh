#!/usr/bin/env bash
# =====================================================================================
# Agent U -- rebuild the LFE family LOCATION step and calibrate it against Bostock's
# PUBLISHED family coordinates, then apply it to Washington.
#
#   cd /home/jovyan/tremorferometry/lfe_array_discovery
#   nohup setsid -f bash scripts/u_locate_queue.sh >> logs/u_locate_queue.log 2>&1 < /dev/null
#
# Markers: results/U_LOCATE_BASELINE_DONE (interim, after the baseline measurement),
#          results/U_LOCATE_DONE | results/U_LOCATE_FAILED.
# Every stage is guarded by its own output, so a relaunch resumes instead of redoing hours.
# Detection and growth are NEVER touched; this queue only re-reads their outputs.
# =====================================================================================
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery || exit 1
mkdir -p logs locks results/u_locate

exec 9>locks/u_locate.lock
if ! flock -n 9; then echo "[U] another u_locate_queue is running; exiting"; exit 0; fi

STAGE="init"
fail() { echo "[U] FAILED at stage $STAGE ($1)"; date -u +"%Y-%m-%dT%H:%M:%SZ" > results/U_LOCATE_FAILED; echo "stage=$STAGE" >> results/U_LOCATE_FAILED; }
on_exit() { rc=$?; if [ ! -s results/U_LOCATE_DONE ] && [ ! -s results/U_LOCATE_FAILED ]; then fail "exit rc=$rc"; fi; }
trap 'fail SIGTERM; exit 143' TERM
trap 'fail SIGINT;  exit 130' INT
trap 'fail SIGHUP;  exit 129' HUP
trap on_exit EXIT

rm -f results/U_LOCATE_FAILED results/U_LOCATE_DONE
U=results/u_locate
REPRO=results/repro
O=results/o_oly2010
say() { echo "[U $(date -u +%H:%M:%S)] $*"; }

# ---------------------------------------------------------------- 1. NORTH long stacks --------
STAGE=north_stacks
if [ ! -s "$U/north_stacks.npz" ]; then
  say "rebuilding the 44-s long stacks of every grown 2005 family from its own detections"
  env PYTHONUNBUFFERED=1 python scripts/uu_stacks.py --bed north \
      --det "$REPRO/grow15/det_it0.parquet" --out "$U/north_stacks.npz" \
      > logs/u_north_stacks.log 2>&1 || { fail "uu_stacks north"; exit 1; }
fi
tail -2 logs/u_north_stacks.log

# ---------------------------------------------------------------- 2. BASELINE (interim marker)
STAGE=baseline
if [ ! -s "$U/north_eval_base/ladder.csv" ]; then
  say "measuring the CURRENT locator against Bostock's published coordinates (never done before)"
  env PYTHONUNBUFFERED=1 python scripts/uu_eval.py --locdir "$U/__none__" \
      --out "$U/north_eval_base" > logs/u_baseline.log 2>&1 || { fail "baseline eval"; exit 1; }
fi
cat logs/u_baseline.log
date -u +"%Y-%m-%dT%H:%M:%SZ" > results/U_LOCATE_BASELINE_DONE
cat "$U/north_eval_base/ladder.csv" >> results/U_LOCATE_BASELINE_DONE

# ---------------------------------------------------------------- 3. NORTH picks ---------------
STAGE=north_picks
if [ ! -s "$U/north/picks.parquet" ]; then
  say "picking S (envelope + AIC onset), P and S-P on every long stack"
  env PYTHONUNBUFFERED=1 python scripts/uu_picks.py --stacks "$U/north_stacks.npz" \
      --out "$U/north/picks.parquet" > logs/u_north_picks.log 2>&1 || { fail "uu_picks"; exit 1; }
fi
tail -2 logs/u_north_picks.log

STAGE=north_resid
if [ ! -s "$U/north/pick_residuals.csv" ]; then
  say "per-station pick residuals against the published locations"
  env PYTHONUNBUFFERED=1 python scripts/uu_resid.py --picks "$U/north/picks.parquet" \
      --out "$U/north/pick_residuals.csv" > logs/u_north_resid.log 2>&1 \
      || { fail "uu_resid"; exit 1; }
fi
tail -3 logs/u_north_resid.log

# ---------------------------------------------------------------- 4. the LADDER ----------------
STAGE=north_ladder
# uu_locate.py already skips any variant whose loc_*.csv exists, so this is idempotent; the
# sentinel additionally spares the relaunch the 267 MB stack load once the ladder is complete.
if [ ! -s "$U/north/.ladder_done" ]; then
  say "running the locator ladder on all 1,437 grown northern templates"
  env PYTHONUNBUFFERED=1 python scripts/uu_locate.py --bed north \
      --picks "$U/north/picks.parquet" --stacks "$U/north_stacks.npz" --outdir "$U/north" \
      >> logs/u_north_ladder.log 2>&1 || { fail "uu_locate ladder"; exit 1; }
  date -u +"%Y-%m-%dT%H:%M:%SZ" > "$U/north/.ladder_done"
else
  say "ladder already complete ($(cat "$U/north/.ladder_done")); skipping"
fi
tail -4 logs/u_north_ladder.log

STAGE=north_eval
env PYTHONUNBUFFERED=1 python scripts/uu_eval.py --locdir "$U/north" --out "$U/north_eval" \
    > logs/u_north_eval.log 2>&1 || { fail "uu_eval"; exit 1; }
cat logs/u_north_eval.log

# ---------------------------------------------------------------- 5. adopt ---------------------
STAGE=adopt
env PYTHONUNBUFFERED=1 python scripts/uu_adopt.py --ladder "$U/north_eval/ladder.csv" \
    --out "$U/north_eval/adopted.json" > logs/u_adopt.log 2>&1 || { fail "uu_adopt"; exit 1; }
cat logs/u_adopt.log
VAR=$(tail -1 logs/u_adopt.log)
case "$VAR" in (L*) ;; (*) fail "adopt returned '$VAR'"; exit 1;; esac
VAR2=L12_comb     # the maximum-COVERAGE companion: identical estimator plus the absolute-pick
                  # term, so families with too few correlating station pairs still locate.  Both
                  # are carried through Washington; neither filters anything.
say "adopted variant = $VAR ; coverage companion = $VAR2"

# ---------------------------------------------------------------- 6. uncertainty + control ----
STAGE=north_boot
if [ ! -s "$U/north_boot/loc_${VAR2}.csv" ]; then
  say "station bootstrap + shuffled-station control for $VAR and $VAR2"
  env PYTHONUNBUFFERED=1 python scripts/uu_locate.py --bed north \
      --picks "$U/north/picks.parquet" --stacks "$U/north_stacks.npz" \
      --outdir "$U/north_boot" --only "$VAR,$VAR2" --boot 60 --nshuf 20 \
      > logs/u_north_boot.log 2>&1 || { fail "uu_locate boot"; exit 1; }
fi
tail -2 logs/u_north_boot.log
cp -f "$U/north_boot/loc_${VAR}.csv" "$U/north/loc_${VAR}.csv"
cp -f "$U/north_boot/loc_${VAR2}.csv" "$U/north/loc_${VAR2}.csv"
date -u +"%Y-%m-%dT%H:%M:%SZ" > "$U/north/.ladder_done"   # the ladder outputs now carry boot cols

# ---------------------------------------------------------------- 7. double difference --------
STAGE=north_dd
if [ ! -s "$U/north/loc_ZZdd.csv" ]; then
  say "double-difference pass between neighbouring families"
  env PYTHONUNBUFFERED=1 python scripts/uu_dd.py --bed north --loc "$U/north/loc_${VAR}.csv" \
      --picks "$U/north/picks.parquet" --stacks "$U/north_stacks.npz" \
      --out "$U/north/loc_ZZdd.csv" > logs/u_north_dd.log 2>&1 || say "dd failed (continuing)"
fi
tail -4 logs/u_north_dd.log 2>/dev/null || true

STAGE=north_eval2
env PYTHONUNBUFFERED=1 python scripts/uu_eval.py --locdir "$U/north" --out "$U/north_eval" \
    > logs/u_north_eval2.log 2>&1 || { fail "uu_eval 2"; exit 1; }
cat logs/u_north_eval2.log
env PYTHONUNBUFFERED=1 python scripts/uu_adopt.py --ladder "$U/north_eval/ladder.csv" \
    --out "$U/north_eval/adopted.json" >> logs/u_adopt.log 2>&1 || true
VAR=$(tail -1 logs/u_adopt.log)
case "$VAR" in (L*) ;; (*) VAR=$(python -c "import json;print(json.load(open('$U/north_eval/adopted.json'))['variant'])");; esac
say "adopted (final) = $VAR"

STAGE=report_partial
env PYTHONUNBUFFERED=1 python scripts/uu_report.py > logs/u_report.log 2>&1 || true

# ---------------------------------------------------------------- 8. WASHINGTON ----------------
wa_pop() {   # $1 = tag, $2 = template npz, $3 = detections, $4 = extra env, $5 = bed
  local TAG="$1" TPL="$2" DET="$3" E="$4" BED="$5"
  STAGE="wa_${TAG}_stacks"
  if [ ! -s "$U/wa_${TAG}_stacks.npz" ]; then
    say "WA $TAG: long stacks"
    env $E PYTHONUNBUFFERED=1 python scripts/uu_stacks.py --bed "$BED" --det "$DET" \
        --tpl "$TPL" --out "$U/wa_${TAG}_stacks.npz" \
        > "logs/u_wa_${TAG}_stacks.log" 2>&1 || { fail "wa $TAG stacks"; return 1; }
  fi
  tail -2 "logs/u_wa_${TAG}_stacks.log"
  STAGE="wa_${TAG}_picks"
  if [ ! -s "$U/wa_${TAG}/picks.parquet" ]; then
    env $E PYTHONUNBUFFERED=1 python scripts/uu_picks.py --stacks "$U/wa_${TAG}_stacks.npz" \
        --out "$U/wa_${TAG}/picks.parquet" \
        > "logs/u_wa_${TAG}_picks.log" 2>&1 || { fail "wa $TAG picks"; return 1; }
  fi
  tail -2 "logs/u_wa_${TAG}_picks.log"
  STAGE="wa_${TAG}_locate"
  if [ ! -s "$U/wa_${TAG}/loc_${VAR2}.csv" ]; then
    say "WA $TAG: relocating with $VAR and $VAR2 (+ bootstrap + shuffled-station control)"
    env $E PYTHONUNBUFFERED=1 python scripts/uu_locate.py --bed "$BED" \
        --picks "$U/wa_${TAG}/picks.parquet" --stacks "$U/wa_${TAG}_stacks.npz" \
        --outdir "$U/wa_${TAG}" --only "$VAR,$VAR2" --boot 40 --nshuf 20 \
        > "logs/u_wa_${TAG}_locate.log" 2>&1 || { fail "wa $TAG locate"; return 1; }
  fi
  tail -2 "logs/u_wa_${TAG}_locate.log"
  return 0
}

wa_pop real "$O/2010/grow/tpl_it0.npz" "$O/2010/grow/det_it0.parquet" \
       "O_BED=2010 O_OUT=$O/2010" wa || exit 1
env PYTHONUNBUFFERED=1 python scripts/uu_report.py >> logs/u_report.log 2>&1 || true
wa_pop null "$O/2010_null/grow/tpl_it0.npz" "$O/2010_null/grow/det_it0.parquet" \
       "O_BED=2010 O_OUT=$O/2010_null O_NULL_SHIFT=1 O_NULL_SEED=11" wa || exit 1
wa_pop armb "$O/chestler_arm/tpl_final.npz" "$O/chestler_arm/det_final.parquet" \
       "O_BED=2010 O_OUT=$O/chestler_arm" wa_mem || exit 1

STAGE=wa_chestler
say "Chestler-coordinate comparison: real vs null vs Arm B, ranked, nothing filtered"
for V in "$VAR" "$VAR2"; do
  OUTD="$U/wa"; [ "$V" = "$VAR2" ] && OUTD="$U/wa_cov"
  env PYTHONUNBUFFERED=1 python scripts/uu_wa.py \
      --real "$U/wa_real/loc_${V}.csv" --null "$U/wa_null/loc_${V}.csv" \
      --armb "$U/wa_armb/loc_${V}.csv" \
      --det-real "$O/2010/grow/det_it0.parquet" \
      --det-null "$O/2010_null/grow/det_it0.parquet" \
      --det-armb "$O/chestler_arm/det_final.parquet" \
      --sf-real results/sf_final/table_oly2010_real.csv \
      --sf-null results/sf_final/table_oly2010_null.csv \
      --sf-col dedupe_n_partners --sf-min 20 \
      --out "$OUTD" > "logs/u_wa_chest_${V}.log" 2>&1 || { fail "uu_wa $V"; exit 1; }
  cat "logs/u_wa_chest_${V}.log"
done

# ---------------------------------------------------------------- 9. report --------------------
STAGE=report
env PYTHONUNBUFFERED=1 python scripts/uu_report.py > logs/u_report.log 2>&1 \
    || { fail "uu_report"; exit 1; }
tail -2 logs/u_report.log

date -u +"%Y-%m-%dT%H:%M:%SZ" > results/U_LOCATE_DONE
echo "adopted=$VAR" >> results/U_LOCATE_DONE
say "DONE"
