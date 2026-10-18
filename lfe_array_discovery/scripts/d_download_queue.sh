#!/usr/bin/env bash
# lfe_array_discovery - agent D download + metadata queue (ONE background queue, no idle waiting).
# Tranche order is Merlin's (results/merlin_review_plan_2026-09-20.md "(d) Ordering so nobody
# idles"): metadata -> N-Olympic 2010 thin -> 2010 quiet -> XG rest -> XU -> SEALED 2011 ->
# Puget 2008 -> CAFE 2007 -> intraslab EQ windows -> extension tail 2007-2026.
# Markers: results/markers/DL_<TRANCHE>_DONE|FAILED, final results/DL_ALL_DONE.
set -u
cd "$(dirname "$0")/.." || exit 1
ROOT="$PWD"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY=python
MK="$ROOT/results/markers"
mkdir -p "$MK" "$ROOT/logs"

step () {   # step <TRANCHE> <command...>
  local tag="$1"; shift
  if [ -f "$MK/DL_${tag}_DONE" ]; then
    echo "[queue] $tag already done, skipping"; return 0
  fi
  echo "[queue] ===== $tag start $(date -u +%FT%TZ) ====="
  if "$@"; then
    date -u +%FT%TZ > "$MK/DL_${tag}_DONE"
    echo "[queue] ===== $tag DONE $(date -u +%FT%TZ) ====="
  else
    date -u +%FT%TZ > "$MK/DL_${tag}_FAILED"
    echo "[queue] ===== $tag FAILED (rc=$?) $(date -u +%FT%TZ) ====="
  fi
}

# ---- (1) metadata + queues (idempotent; StationXML is cached on disk)
[ -f "$MK/DL_META_DONE" ] || $PY scripts/d_meta.py
[ -f "$ROOT/results/queues/q_T7_extension.csv" ] || $PY scripts/d_queues.py

# ---- (2) N-Olympic 2010: thin set -> quiet days -> remaining XG -> XU  (cache after each)
step T2A  $PY scripts/d_fetch.py results/queues/q_T2A_thin.csv
step T2A_CACHE $PY scripts/d_cache.py T2A --workers 24
step T2B  $PY scripts/d_fetch.py results/queues/q_T2B_quiet.csv
step T2B_CACHE $PY scripts/d_cache.py T2B --workers 24
step T2C  $PY scripts/d_fetch.py results/queues/q_T2C_xg_rest.csv
step T2C_CACHE $PY scripts/d_cache.py T2C --workers 24
step T2D  $PY scripts/d_fetch.py results/queues/q_T2D_xu.csv
step T2D_CACHE $PY scripts/d_cache.py T2D --workers 24

# ---- (3) SEALED 2011 (traces_sealed/, hashed, never cached, never read)
step T3   $PY scripts/d_fetch.py results/queues/q_T3_sealed.csv
step T3_SEAL $PY scripts/d_seal.py

# ---- (4) Puget / S-Olympic 2008 + literature-quiet days
step T4   $PY scripts/d_fetch.py results/queues/q_T4_puget2008.csv
step T4_CACHE $PY scripts/d_cache.py T4 --workers 24

# ---- (5) CAFE 2007 + literature-quiet days
step T5   $PY scripts/d_fetch.py results/queues/q_T5_cafe2007.csv
step T5_CACHE $PY scripts/d_cache.py T5 --workers 24

# ---- (6) intraslab earthquake windows (+-3 min), for the station corrections
step T6   $PY scripts/d_eq.py

# ---- (7) extension roster background tail, ETS months first then oldest-first
step T7   $PY scripts/d_fetch.py results/queues/q_T7_extension.csv

date -u +%FT%TZ > "$ROOT/results/DL_ALL_DONE"
echo "[queue] ALL TRANCHES COMPLETE $(date -u +%FT%TZ)"
