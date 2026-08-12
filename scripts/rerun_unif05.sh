#!/usr/bin/env bash
# Re-run dv/v for the stations whose ORIGINAL (eps_max=0.02) output is most contaminated by the
# +/-2% search rail (fleet diagnostic 2026-08-12: these carry 4.5-37% railed rows among cc>0.7 input,
# vs a fleet median of 0.28%). UNIFORM eps_max=0.05 with n_eps=501 -> grid step 0.0002 (0.02%),
# identical to the original 0.02/201, so only the search RANGE changes, not its resolution.
# Writes *_unif05.csv -- does NOT touch the originals (still the inversion's input) or the
# Jul-28 *_eps05.csv files (which used inconsistent per-station ranges 0.023-0.049).
set -u
PY=/home/jovyan/envs/tremorferometry/bin/python
LOG=logs/rerun_unif05.log
STATIONS="NEMA B036p90f40 B035p90f40 BABR KBO KRP B031p90f40 B017p90f40 B012p90f40 LCM"

mkdir -p logs
echo "=== uniform eps_max=0.05 / n_eps=501 re-run started $(date -u +%FT%TZ) ===" >> "$LOG"
for T in $STATIONS; do
  NPZ="data/long_window_daily_${T}_Z.npz"
  OUT="data/daily_dvv_${T}_Z_2to4_unif05.csv"
  if [ ! -f "$NPZ" ]; then echo "[$T] SKIP - no npz" >> "$LOG"; continue; fi
  if [ -f "$OUT" ]; then echo "[$T] SKIP - $OUT already exists" >> "$LOG"; continue; fi
  # only the causality-certified families are used downstream (map + inversion both filter to them);
  # per-family independence makes this exactly equivalent for the kept families. ~4.3x less work.
  STEM=$(echo "$T" | tr 'A-Z' 'a-z' | sed 's/p90f40//')
  CERT="data/${STEM}_causality_cert.csv"
  CERTARG=""; [ -f "$CERT" ] && CERTARG="--cert-csv $CERT"
  echo "--- [$T] start $(date -u +%FT%TZ) ${CERTARG:-（all families）} ---" >> "$LOG"
  S=$(date +%s)
  "$PY" scripts/dvv_roll30cal.py --station "$T" --npz "$NPZ" --out "$OUT" \
        --eps-max 0.05 --n-eps 501 --origin-anchor --workers 28 $CERTARG >> "$LOG" 2>&1
  RC=$?
  echo "[$T] rc=$RC elapsed=$(( $(date +%s) - S ))s" >> "$LOG"
done
echo "=== ALL DONE $(date -u +%FT%TZ) ===" >> "$LOG"
