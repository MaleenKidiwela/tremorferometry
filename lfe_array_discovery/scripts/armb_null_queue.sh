#!/usr/bin/env bash
# ARMB_NULL -- scrambled control (site-day shift) for Arm B, Chestler's own chain on the 2010
# N-Olympic arrays.  Identical code + parameters to scripts/o_oly2010_queue.sh sec.9
# (o_chestler.py seeds -> grow --max-seeds 1200 -> dedupe) and to the Arm-B long re-stack in
# scripts/u_locate_queue.sh (uu_stacks.py --bed wa_mem --tpl tpl_final.npz), except that every
# python stage runs through scripts/armb_null_shift.py with O_NULL_SHIFT=1 O_NULL_SEED=11: each
# SITE's day-series is circularly shifted by arm A's whole-day offsets
# (results/o_oly2010/2010_null/null_offsets.json), members of a site shift together.
#
#   nohup setsid -f bash scripts/armb_null_queue.sh >> logs/armb_null_queue.log 2>&1 < /dev/null
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
mkdir -p locks logs results/o_oly2010 results/u_locate results/ab_royer/OLY2010_cvm

exec 9>locks/armb_null.lock
flock -n 9 || { echo "[AN] already running"; exit 0; }
rm -f results/ARMB_NULL_DONE results/ARMB_NULL_FAILED
echo $$ > results/ARMB_NULL_PID
STAGE=start
fail() { [ -f results/ARMB_NULL_DONE ] || \
  printf 'stage=%s signal=%s at=%s pid=%s\n' "$STAGE" "${1:-exit}" "$(date -u +%FT%TZ)" "$$" \
  > results/ARMB_NULL_FAILED; }
trap 'rc=$?; [ -n "${HB:-}" ] && kill "$HB" 2>/dev/null; echo "[AN] TRAP EXIT rc=$rc stage=$STAGE"; [ $rc -ne 0 ] && fail "exit$rc"' EXIT
trap 'echo "[AN] TRAP TERM stage=$STAGE"; fail TERM; exit 143' TERM
trap 'echo "[AN] TRAP INT  stage=$STAGE"; fail INT;  exit 130' INT
trap 'echo "[AN] TRAP HUP  stage=$STAGE"; fail HUP;  exit 129' HUP

# heartbeat (own child, killed on exit)
( while true; do echo "$(date -u +%FT%TZ) pid=$$" > results/ARMB_NULL_HEARTBEAT; sleep 60; done ) &
HB=$!

CU=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')
export LD_LIBRARY_PATH="${CU}${LD_LIBRARY_PATH:-}"

O=results/o_oly2010
ARMBN=$O/chestler_arm_null
U=results/u_locate
mkdir -p "$ARMBN"
E="O_BED=2010 O_OUT=$ARMBN O_NULL_SHIFT=1 O_NULL_SEED=11 PYTHONUNBUFFERED=1"
W="python scripts/armb_null_shift.py"
echo "[AN] start $(date -u +%FT%TZ) pid=$$"

STAGE=seeds
if [ ! -s "$ARMBN/seeds.csv" ]; then
  env $E $W scripts/o_chestler.py --stage seeds --out "$ARMBN" \
    > logs/armb_null_seeds.log 2>&1 || exit 1
fi
tail -2 logs/armb_null_seeds.log

STAGE=grow
if [ ! -s "$ARMBN/det_it3.parquet" ]; then
  env $E $W scripts/o_chestler.py --stage grow --out "$ARMBN" --max-seeds 1200 \
    > logs/armb_null_grow.log 2>&1 || exit 1
fi
tail -3 logs/armb_null_grow.log

STAGE=dedupe
if [ ! -s "$ARMBN/keep.csv" ]; then
  env $E $W scripts/o_chestler.py --stage dedupe --out "$ARMBN" \
    > logs/armb_null_dedupe.log 2>&1 || exit 1
fi
tail -2 logs/armb_null_dedupe.log
cat "$ARMBN/dedupe.json" 2>/dev/null

STAGE=stacks
if [ ! -s "$U/wa_armb_null_stacks.npz" ]; then
  env $E $W scripts/uu_stacks.py --bed wa_mem --det "$ARMBN/det_final.parquet" \
    --tpl "$ARMBN/tpl_final.npz" --out "$U/wa_armb_null_stacks.npz" \
    > logs/armb_null_stacks.log 2>&1 || exit 1
fi
tail -2 logs/armb_null_stacks.log

STAGE=score
python scripts/armb_site_score.py --stacks "$U/wa_armb_null_stacks.npz" \
  --out results/ab_royer/OLY2010_cvm/site_score_armb_null_siteavg.csv \
  > logs/armb_null_score.log 2>&1 || exit 1
# the REAL Arm-B stacks through the SAME scorer (like-for-like; new file, nothing overwritten)
python scripts/armb_site_score.py --stacks "$U/wa_armb_stacks.npz" \
  --out results/ab_royer/OLY2010_cvm/site_score_armb_real_rescored.csv \
  >> logs/armb_null_score.log 2>&1 || exit 1
cat logs/armb_null_score.log

STAGE=done
echo "OK $(date -u +%FT%TZ)" > results/ARMB_NULL_DONE
echo "[AN] DONE $(date -u +%FT%TZ)"
