#!/usr/bin/env bash
# M_REFIRE2011 -- confirmation of the WA 2010 candidate families by re-firing the UNCHANGED
# templates on the SEALED 2011 window, plus a wider-grid re-location of both years, as ONE
# detached queue.
#
#   nohup setsid -f bash scripts/m_refire2011_queue.sh >> logs/m_refire2011_queue.log 2>&1 < /dev/null
#
# Every stage pins its OWN environment with `env VAR=... python` -- nothing is exported across
# stages (the K_YEAR leak that voided a blind run must not be repeatable), and every stage has an
# output guard so the queue resumes from where it stopped.  Never `pkill -f` this.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
mkdir -p locks logs results/m_2011

exec 9>locks/m_refire2011.lock
flock -n 9 || { echo "[M] already running"; exit 0; }
rm -f results/M_REFIRE2011_DONE results/M_REFIRE2011_FAILED
STAGE=start
fail() { [ -f results/M_REFIRE2011_DONE ] || \
  printf 'stage=%s signal=%s at=%s pid=%s\n' "$STAGE" "${1:-exit}" "$(date -u +%FT%TZ)" "$$" \
  > results/M_REFIRE2011_FAILED; }
trap 'rc=$?; echo "[M] TRAP EXIT rc=$rc stage=$STAGE"; fail "exit$rc"' EXIT
trap 'echo "[M] TRAP TERM stage=$STAGE"; fail TERM; exit 143' TERM
trap 'echo "[M] TRAP INT  stage=$STAGE"; fail INT;  exit 130' INT
trap 'echo "[M] TRAP HUP  stage=$STAGE"; fail HUP;  exit 129' HUP

# CUDA-12 runtime libs (memory: torch dies with libcufft.so.11 missing after a pod restart)
CU=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')
export LD_LIBRARY_PATH="${CU}${LD_LIBRARY_PATH:-}"

O=results/m_2011
R=results/l_wa2010
N=results/l_wa2010_null
echo "[M] start $(date -u +%FT%TZ) pid=$$ ppid=$PPID"

# ---------------------------------------------------------------- 1. open the seal --------------
STAGE=seal
if ! ls results/SEALED_2011_OPENED_* >/dev/null 2>&1; then
  env python scripts/m_seal_open.py --workers 16 > logs/m_seal.log 2>&1 || exit 1
fi
tail -8 logs/m_seal.log 2>/dev/null || true
ls results/SEALED_2011_OPENED_* >/dev/null 2>&1 || { echo "[M] seal not opened"; exit 1; }

# ---------------------------------------------------------------- 2. cache 2011 -----------------
STAGE=cache
if [ ! -s "$O/CACHE_DONE" ]; then
  env python scripts/m_cache2011.py --workers 24 --min-cov 0.2 \
    > logs/m_cache2011.log 2>&1 || exit 1
  echo "OK $(date -u +%FT%TZ)" > "$O/CACHE_DONE"
fi
tail -12 logs/m_cache2011.log

# ---------------------------------------------------------------- 3. health veto 2011 -----------
STAGE=health
if [ ! -s results/wa_station_hour_veto_2011.csv ]; then
  env python scripts/b_wa_health.py --year 2011 --workers 24 \
    > logs/m_health2011.log 2>&1 || exit 1
fi
tail -20 logs/m_health2011.log

# ---------------------------------------------------------------- 4. the 2011 bed ---------------
STAGE=bed
if [ ! -s "$O/bed.json" ]; then
  env M_OUT="$O" python scripts/m_bed2011.py > logs/m_bed2011.log 2>&1 || exit 1
fi
tail -40 logs/m_bed2011.log
[ -s "$O/tremor_dates.txt" ] || { echo "[M] no tremor_dates.txt"; exit 1; }
TREMOR=$(cat "$O/tremor_dates.txt")
FA=$(cat "$O/nontremor_dates.txt")
echo "[M] 2011 tremor days: $(echo "$TREMOR" | tr ',' '\n' | wc -l); "\
"non-tremor (FA) days: $(echo "$FA" | tr ',' '\n' | wc -l)"

# ---------------------------------------------------------------- 5. RE-FIRE --------------------
# $1 tag (real|null)  $2 template npz  $3 candidate csv  $4 dates  $5 suffix
refire() {
  STAGE="refire_$1_$5"
  if [ ! -s "$O/det_$1_$5.parquet" ]; then
    env M_OUT="$O" python scripts/m_refire.py --tpl "$2" --grps "$3" --dates "$4" \
      --out "$O" --tag "$1_$5" > "logs/m_refire_$1_$5.log" 2>&1 || exit 1
  fi
  tail -14 "logs/m_refire_$1_$5.log"
}
refire real "$R/grow/tpl_it0.npz" "$R/verify_it0.csv" "$TREMOR" tremor
refire null "$N/grow/tpl_it0.npz" "$N/verify_it0.csv" "$TREMOR" tremor
refire real "$R/grow/tpl_it0.npz" "$R/verify_it0.csv" "$FA"     fa
refire null "$N/grow/tpl_it0.npz" "$N/verify_it0.csv" "$FA"     fa

# ---------------------------------------------------------------- 6. WIDE-GRID RE-LOCATION ------
# 2010, from the 2010 detections (the null chain must be re-stacked on its own shifted network)
STAGE=loc2010_real
if [ ! -s "$O/loc2010_real.csv" ]; then
  env M_OUT="$O" python scripts/m_relocate.py --bed 2010 --det "$R/grow/det_it0.parquet" \
    --grps "$R/verify_it0.csv" --out "$O/loc2010_real" --tag loc2010_real \
    > logs/m_loc2010_real.log 2>&1 || exit 1
fi
tail -6 logs/m_loc2010_real.log

STAGE=loc2010_null
if [ ! -s "$O/loc2010_null.csv" ]; then
  env M_OUT="$O" L_NULL_SHIFT=1 L_NULL_SEED=11 python scripts/m_relocate.py --bed 2010 \
    --det "$N/grow/det_it0.parquet" --grps "$N/verify_it0.csv" --out "$O/loc2010_null" \
    --tag loc2010_null > logs/m_loc2010_null.log 2>&1 || exit 1
fi
tail -6 logs/m_loc2010_null.log

# 2011, from the 2011 re-fire detections (both template sets fire on the SAME real 2011 data)
STAGE=loc2011_real
if [ ! -s "$O/loc2011_real.csv" ]; then
  env M_OUT="$O" python scripts/m_relocate.py --bed 2011 --det "$O/det_real_tremor.parquet" \
    --grps "$R/verify_it0.csv" --out "$O/loc2011_real" --tag loc2011_real \
    > logs/m_loc2011_real.log 2>&1 || exit 1
fi
tail -6 logs/m_loc2011_real.log

STAGE=loc2011_null
if [ ! -s "$O/loc2011_null.csv" ]; then
  env M_OUT="$O" python scripts/m_relocate.py --bed 2011 --det "$O/det_null_tremor.parquet" \
    --grps "$N/verify_it0.csv" --out "$O/loc2011_null" --tag loc2011_null \
    > logs/m_loc2011_null.log 2>&1 || exit 1
fi
tail -6 logs/m_loc2011_null.log

# ---------------------------------------------------------------- 7. report ---------------------
STAGE=report
env M_OUT="$O" python scripts/m_report.py --out results/WA_2011_REFIRE_RESULTS.md \
  > logs/m_report.log 2>&1 || exit 1
cat logs/m_report.log

echo "OK $(date -u +%FT%TZ)" > results/M_REFIRE2011_DONE
echo "[M] DONE $(date -u +%FT%TZ)"
