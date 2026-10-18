#!/usr/bin/env bash
# S_FINAL -- Washington, re-scored under RANK-DON'T-DISCARD (user rule, 2026-09-23).
#
#   nohup setsid -f bash scripts/s_final_queue.sh >> logs/s_final_queue.log 2>&1 < /dev/null
#
# Nothing is filtered.  Every grown template of every chain keeps its row; every criterion is a
# column; the whole-chain site-shift null gets the byte-identical table; the evidence is the
# real-vs-null separation as a function of each criterion.  The only list-reducing steps are the
# two the published chains actually perform -- the dedupe and the location -- and both are
# recorded as flags with the rows retained.
#
# Stage 1  (done earlier in this queue's first incarnation) the confirmation rule calibrated on
#          Bostock's own families, north 2005 half-split -- reported, never used as a gate.
# Stage 2  the published dedupe (a flag) for every chain.
# Stage 3  Arm B: the site-shift null it never had, and the locations it never got.
# Stage 4  sealed-window re-fire of EVERY template of EVERY chain (2011 for the 2010 beds,
#          2007 for the 2008 bed), so the re-fire outcome is a column on every row.
# Stage 5  per-template tables + real-vs-null separation curves + the report.
#
# Every stage pins its own environment with `env VAR=... python`; every stage has an output guard,
# so the queue resumes where it stopped.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
mkdir -p locks logs results/sf_final results/sf_final/north

exec 9>locks/s_final.lock
flock -n 9 || { echo "[S] already running"; exit 0; }
rm -f results/S_FINAL_DONE results/S_FINAL_FAILED
STAGE=start
fail() { [ -f results/S_FINAL_DONE ] || \
  printf 'stage=%s signal=%s at=%s pid=%s\n' "$STAGE" "${1:-exit}" "$(date -u +%FT%TZ)" "$$" \
  > results/S_FINAL_FAILED; }
trap 'rc=$?; echo "[S] TRAP EXIT rc=$rc stage=$STAGE"; fail "exit$rc"' EXIT
trap 'echo "[S] TRAP TERM stage=$STAGE"; fail TERM; exit 143' TERM
trap 'echo "[S] TRAP INT  stage=$STAGE"; fail INT;  exit 130' INT
trap 'echo "[S] TRAP HUP  stage=$STAGE"; fail HUP;  exit 129' HUP

CU=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')
export LD_LIBRARY_PATH="${CU}${LD_LIBRARY_PATH:-}"

F=results/sf_final
N=$F/north
O10=results/o_oly2010
P08=results/n_puget2008
ARMB=$O10/chestler_arm
ARMBN=$F/armb_null
SEAL=$F/sealed
mkdir -p "$ARMBN" "$SEAL"
echo "[S] start $(date -u +%FT%TZ) pid=$$ ppid=$PPID"

# =============================================================== STAGE 1 =======================
# the north half-split calibration (reported as information; never a gate)
STAGE=n_select
[ -s "$N/set.csv" ] || python scripts/sf_north_split.py --stage select \
  > logs/s_n_select.log 2>&1 || exit 1
STAGE=n_h1
[ -s "$N/tpl_h1.npz" ] || python scripts/sf_north_split.py --stage h1 \
  > logs/s_n_h1.log 2>&1 || exit 1
STAGE=n_refire_h2
[ -s "$N/det_h2.parquet" ] || python scripts/sf_north_split.py --stage refire --which h2 \
  > logs/s_n_refire_h2.log 2>&1 || exit 1
STAGE=n_refire_quiet
[ -s "$N/det_quiet.parquet" ] || python scripts/sf_north_split.py --stage refire --which quiet \
  > logs/s_n_refire_quiet.log 2>&1 || exit 1
STAGE=n_h2stack
[ -s "$N/tpl_h2.npz" ] || python scripts/sf_north_split.py --stage h2stack \
  > logs/s_n_h2stack.log 2>&1 || exit 1
for H in 1 2; do
  STAGE="n_locate_h$H"
  [ -s "$N/loc_h$H.csv" ] || python scripts/sf_north_split.py --stage locate --half "$H" \
    > "logs/s_n_loc_h$H.log" 2>&1 || exit 1
done
STAGE=n_calib
[ -s "$F/calib.json" ] || python scripts/sf_calib.py > logs/s_n_calib.log 2>&1 || exit 1
cat "$F/calib.json"

# =============================================================== STAGE 2 =======================
# the published dedupe, as a FLAG, for every chain (rows are never deleted)
dedupe() {                       # $1 key  $2 chain  $3 det  $4 tpl  $5 loc
  STAGE="dedupe_$1_$2"
  [ -s "$F/$1_$2_dedupe.csv" ] && { echo "[S] dedupe $1 $2 already done"; return 0; }
  python scripts/sf_accept.py --det "$3" --tpl "$4" --loc "$5" --out "$F/$1_$2" \
    --tag "$1_$2" > "logs/s_dedupe_$1_$2.log" 2>&1 || exit 1
  grep -E "dedupe:" "logs/s_dedupe_$1_$2.log" || tail -5 "logs/s_dedupe_$1_$2.log"
}
dedupe oly2010 real "$O10/2010/grow/det_it0.parquet" "$O10/2010/grow/tpl_it0.npz" \
  "$O10/2010/loc_it0.csv"
dedupe oly2010 null "$O10/2010_null/grow/det_it0.parquet" "$O10/2010_null/grow/tpl_it0.npz" \
  "$O10/2010_null/loc_it0.csv"
dedupe puget2008 real "$P08/2008/grow/det_it0.parquet" "$P08/2008/grow/tpl_it0.npz" \
  "$P08/2008/loc_it0.csv"
dedupe puget2008 null "$P08/2008_null/grow/det_it0.parquet" "$P08/2008_null/grow/tpl_it0.npz" \
  "$P08/2008_null/loc_it0.csv"

# =============================================================== STAGE 3 =======================
echo "[S] ===================== ARM B: the null it never had ================="
STAGE=armb_null_seeds
[ -s "$ARMBN/seeds.csv" ] || env O_BED=2010 O_NULL_SHIFT=1 O_NULL_SEED=11 \
  python scripts/sf_armb_null.py --stage seeds --out "$ARMBN" \
  > logs/s_armb_null_seeds.log 2>&1 || exit 1
tail -3 logs/s_armb_null_seeds.log

STAGE=armb_null_grow
[ -s "$ARMBN/det_it3.parquet" ] || env O_BED=2010 O_NULL_SHIFT=1 O_NULL_SEED=11 \
  python scripts/sf_armb_null.py --stage grow --out "$ARMBN" --max-seeds 1200 \
  > logs/s_armb_null_grow.log 2>&1 || exit 1
tail -6 logs/s_armb_null_grow.log

STAGE=armb_null_dedupe_rcom
[ -s "$ARMBN/keep.csv" ] || env O_BED=2010 O_NULL_SHIFT=1 O_NULL_SEED=11 \
  python scripts/sf_armb_null.py --stage dedupe --out "$ARMBN" \
  > logs/s_armb_null_rcom.log 2>&1 || exit 1
tail -3 logs/s_armb_null_rcom.log

STAGE=armb_null_fa
[ -s "$ARMBN/det_fa.parquet" ] || env O_BED=2010 O_NULL_SHIFT=1 O_NULL_SEED=11 \
  python scripts/sf_armb_null.py --stage fa --out "$ARMBN" \
  > logs/s_armb_null_fa.log 2>&1 || exit 1
tail -3 logs/s_armb_null_fa.log

echo "[S] ===================== ARM B: locate both populations =============="
STAGE=armb_locate_real
[ -s "$F/armb_real_loc.csv" ] || env O_BED=2010 python scripts/sf_locate_member.py \
  --tpl "$ARMB/tpl_final.npz" --det "$ARMB/det_final.parquet" --out "$F/armb_real_loc.csv" \
  > logs/s_armb_loc_real.log 2>&1 || exit 1
tail -2 logs/s_armb_loc_real.log
STAGE=armb_locate_null
[ -s "$F/armb_null_loc.csv" ] || env O_BED=2010 python scripts/sf_locate_member.py \
  --tpl "$ARMBN/tpl_final.npz" --det "$ARMBN/det_final.parquet" --out "$F/armb_null_loc.csv" \
  > logs/s_armb_loc_null.log 2>&1 || exit 1
tail -2 logs/s_armb_loc_null.log

dedupe armb real "$ARMB/det_final.parquet" "$ARMB/tpl_final.npz" "$F/armb_real_loc.csv"
dedupe armb null "$ARMBN/det_final.parquet" "$ARMBN/tpl_final.npz" "$F/armb_null_loc.csv"

STAGE=depthcheck
[ -s "$F/depth_check.json" ] || python scripts/sf_depthcheck.py \
  > logs/s_depthcheck.log 2>&1 || exit 1
grep -E "verdict|rejected_genuinely" logs/s_depthcheck.log | head -6

# =============================================================== STAGE 4 =======================
# sealed-window re-fire of EVERY template of EVERY chain -- no --grps, nothing pre-selected
TD11=$(cat "$O10/tremor_dates_2011.txt"); QD11=$(cat "$O10/quiet_dates_2011.txt")
TD07=$(cat "$P08/2007/refire_dates.txt"); QD07=$(cat "$P08/2007/quiet_dates.txt")

refire_all() {   # $1 key  $2 chain  $3 tpl  $4 TD  $5 QD  $6 env  $7 script
  local OD="$SEAL/$1_$2"
  mkdir -p "$OD"
  for W in tremor fa; do
    local DD="$4"; [ "$W" = "fa" ] && DD="$5"
    STAGE="seal_$1_$2_$W"
    [ -s "$OD/det_all_${W}.parquet" ] && { echo "[S] seal $1 $2 $W done"; continue; }
    env $6 python scripts/$7 --tpl "$3" --dates "$DD" --out "$OD" --tag "all_${W}" \
      > "logs/s_seal_$1_$2_$W.log" 2>&1 || exit 1
    tail -5 "logs/s_seal_$1_$2_$W.log"
  done
}
echo "[S] ===================== SEALED-WINDOW RE-FIRE (all templates) ======="
refire_all oly2010 real "$O10/2010/grow/tpl_it0.npz" "$TD11" "$QD11" "O_BED=2011" \
  "o_run.py n_refire"
refire_all oly2010 null "$O10/2010_null/grow/tpl_it0.npz" "$TD11" "$QD11" "O_BED=2011" \
  "o_run.py n_refire"
refire_all puget2008 real "$P08/2008/grow/tpl_it0.npz" "$TD07" "$QD07" "N_BED=2007" \
  "n_refire.py"
refire_all puget2008 null "$P08/2008_null/grow/tpl_it0.npz" "$TD07" "$QD07" "N_BED=2007" \
  "n_refire.py"
refire_all armb real "$ARMB/tpl_final.npz" "$TD11" "$QD11" "O_BED=2011" "sf_armb_refire.py"
refire_all armb null "$ARMBN/tpl_final.npz" "$TD11" "$QD11" "O_BED=2011" "sf_armb_refire.py"

# =============================================================== STAGE 5 =======================
echo "[S] ===================== TABLES + SEPARATION ========================="
HOURS10=results/l_wa2010/candidate_hours.csv
HOURS08=$P08/tremor_hours_2008.csv
TD10=$(cat "$O10/tremor_dates_2010.txt"); QD10=$(cat "$O10/quiet_dates_2010.txt")
TD08=$(cat "$P08/tremor_dates_2008.txt"); QD08=$(cat "$P08/quiet_dates_2008.txt")

table() {  # $1 key $2 chain $3 det $4 detfa $5 loc $6 screen $7 hours $8 TD $9 QD $10 sealTD $11 sealQD
  STAGE="table_$1_$2"
  [ -s "$F/table_$1_$2.csv" ] && { echo "[S] table $1 $2 done"; return 0; }
  local OD="$SEAL/$1_$2"
  python scripts/sf_table.py --key "$1" --chain "$2" --det "$3" --det-fa "$4" --loc "$5" \
    --dedupe "$F/$1_$2_dedupe.csv" --screen "$6" --hours "$7" \
    --tremor-dates "$8" --quiet-dates "$9" \
    --sealed-det "$OD/det_all_tremor.parquet" --sealed-fa "$OD/det_all_fa.parquet" \
    --sealed-tremor-dates "${10}" --sealed-quiet-dates "${11}" --chestler \
    --out "$F/table_$1_$2.csv" > "logs/s_table_$1_$2.log" 2>&1 || exit 1
  tail -22 "logs/s_table_$1_$2.log"
}
table oly2010 real "$O10/2010/grow/det_it0.parquet" "$O10/2010/det_fa.parquet" \
  "$O10/2010/loc_it0.csv" "$O10/2010/screen_final.csv" "$HOURS10" "$TD10" "$QD10" "$TD11" "$QD11"
table oly2010 null "$O10/2010_null/grow/det_it0.parquet" "$O10/2010_null/det_fa.parquet" \
  "$O10/2010_null/loc_it0.csv" "$O10/2010_null/screen_final.csv" "$HOURS10" "$TD10" "$QD10" \
  "$TD11" "$QD11"
table puget2008 real "$P08/2008/grow/det_it0.parquet" "$P08/2008/det_fa.parquet" \
  "$P08/2008/loc_it0.csv" "$P08/2008/screen_final.csv" "$HOURS08" "$TD08" "$QD08" "$TD07" "$QD07"
table puget2008 null "$P08/2008_null/grow/det_it0.parquet" "$P08/2008_null/det_fa.parquet" \
  "$P08/2008_null/loc_it0.csv" "$P08/2008_null/screen_final.csv" "$HOURS08" "$TD08" "$QD08" \
  "$TD07" "$QD07"
table armb real "$ARMB/det_final.parquet" "$ARMB/det_fa.parquet" "$F/armb_real_loc.csv" "" \
  "$HOURS10" "$TD10" "$QD10" "$TD11" "$QD11"
table armb null "$ARMBN/det_final.parquet" "$ARMBN/det_fa.parquet" "$F/armb_null_loc.csv" "" \
  "$HOURS10" "$TD10" "$QD10" "$TD11" "$QD11"

for B in oly2010 puget2008 armb; do
  STAGE="sep_$B"
  [ -s "$F/sep_$B.json" ] || python scripts/sf_separation.py --real "$F/table_${B}_real.csv" \
    --null "$F/table_${B}_null.csv" --key "$B" --out "$F/sep_$B" \
    > "logs/s_sep_$B.log" 2>&1 || exit 1
  tail -30 "logs/s_sep_$B.log"
done

STAGE=report
python scripts/sf_report.py --out results/WA_FINAL_ASSESSMENT.md \
  > logs/s_report.log 2>&1 || exit 1
tail -30 logs/s_report.log

echo "OK $(date -u +%FT%TZ)" > results/S_FINAL_DONE
echo "[S] DONE $(date -u +%FT%TZ)"
