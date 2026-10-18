#!/usr/bin/env bash
# STEP-1 POST-MORTEM DIAGNOSTIC QUEUE (agent G, 2026-09-20).
#
# Runs, in one nohup queue: the FIXED w_diag stage for Step 1b and Step 1c (D1/D2/D3), then the
# D4-D7 diagnostics of scripts/y_step1diag.py, then composes results/STEP1_RESULTS.md (v2) and
# appends the NOTES.md entry.
#
# DIAGNOSTIC ONLY.  No frozen verdict is changed, no engine constant is retuned, nothing is
# downloaded (no FDSN request anywhere in this queue), results/WA_GO is NEVER created, and no
# process is ever killed.  Single instance via flock; every stage writes its own output file and
# is skipped on a relaunch if that file already exists; TERM/INT/HUP/EXIT always leave either
# results/STEP1DIAG_DONE or results/STEP1DIAG_FAILED behind.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
S=scripts
RES=results
LOG=logs/step1diag_queue.log
MDONE=$RES/STEP1DIAG_DONE
MFAIL=$RES/STEP1DIAG_FAILED
mkdir -p logs locks "$RES/markers"

export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export PYTHONUNBUFFERED=1

# ---------------------------------------------------------------- single instance --------------
exec 9>locks/step1diag.lock
if ! flock -n 9; then
  echo "=== [$(date -u +%FT%TZ)] another y_step1diag_queue.sh holds locks/step1diag.lock -- exiting" >&2
  exit 0
fi

STAGE="startup"
FINISHED=0

fail_marker() {
  { echo "step: $STAGE";
    echo "when: $(date -u +%FT%TZ)";
    echo "reason: $1";
    echo "--- last 120 log lines ---";
    tail -120 "$LOG" 2>/dev/null; } > "$MFAIL"
}

on_signal() {
  local sig="$1"
  if [ "$FINISHED" != "1" ]; then
    fail_marker "received SIG$sig during stage $STAGE"
    echo "!!! SIG$sig during $STAGE" >&2
  fi
  exit 143
}

on_exit() {
  local c=$?
  if [ "$FINISHED" != "1" ] && [ ! -f "$MFAIL" ]; then
    fail_marker "queue exited with code $c during stage $STAGE"
  fi
}

trap 'on_signal TERM' TERM
trap 'on_signal INT'  INT
trap 'on_signal HUP'  HUP
trap on_exit EXIT

rm -f "$MDONE" "$MFAIL"
T0=$(date +%s)
echo "=== [$(date -u +%FT%TZ)] y_step1diag_queue START (pid $$)" >&2

# ---------------------------------------------------------------- stage runner -----------------
# step <name> <guard-file> <command...>
step() {
  local name="$1"; shift
  local guard="$1"; shift
  STAGE="$name"
  if [ -n "$guard" ] && [ -f "$guard" ]; then
    echo "=== [$(date -u +%FT%TZ)] SKIP $name (guard $guard exists)" >&2
    return 0
  fi
  echo "=== [$(date -u +%FT%TZ)] STEP $name" >&2
  "$@"
  local c=$?
  if [ $c -ne 0 ]; then
    fail_marker "exit_code $c from: $*"
    echo "!!! FAILED at $name ($c)" >&2
    exit $c
  fi
  if [ -n "$guard" ] && [ ! -f "$guard" ]; then
    fail_marker "stage $name returned 0 but did not write $guard"
    echo "!!! $name produced no output" >&2
    exit 3
  fi
  return 0
}

# ---------------------------------------------------------------- 0. self-test -----------------
STAGE="selftest"
if [ ! -f "$RES/markers/STEP1DIAG_SELFTEST_OK" ]; then
  echo "=== [$(date -u +%FT%TZ)] STEP selftest" >&2
  W_BED=c python - <<'PY' || { STAGE=selftest; fail_marker "selftest failed"; exit 4; }
import sys, json
sys.path.insert(0, "scripts")
import pandas as pd
import y_step1diag as Y, w_cluster as C
Q = pd.read_parquet("results/step1c/ac_quiet.parquet")
for T in (5.0, 6.0, 6.5, 7.0):
    S = Q[Q.net_mad >= T]
    a = sorted([len(t) for _, t in C.components_of(S)], reverse=True)[:8]
    b = sorted([int(x) for x in Y.comp_sizes(S)], reverse=True)[:8]
    assert a == b, (T, a, b)
Tf, binds = Y.fast_threshold(Q, 40)
frozen = json.load(open("results/step1c/T.json"))["T"]
assert abs(Tf - frozen) < 1e-9, (Tf, frozen)
print(f"[selftest] comp_sizes == w_cluster.components_of; fast_threshold reproduces T = {Tf}")
PY
  touch "$RES/markers/STEP1DIAG_SELFTEST_OK"
fi

# ---------------------------------------------------------------- 1. D0: fixed w_diag ----------
STAGE="diag_1c"
if [ ! -f "$RES/step1c_diag.json" ]; then
  echo "=== [$(date -u +%FT%TZ)] STEP diag_1c" >&2
  W_BED=c python $S/w_diag.py --out $RES/step1c_diag --pairs $RES/step1c/ac.parquet \
    || { fail_marker "w_diag (bed c) failed"; exit 5; }
  [ -f "$RES/step1c_diag.json" ] || { fail_marker "diag_1c wrote no json"; exit 5; }
fi

STAGE="diag_1b"
if [ ! -f "$RES/step1b_diag.json" ]; then
  echo "=== [$(date -u +%FT%TZ)] STEP diag_1b" >&2
  W_BED=b python $S/w_diag.py --out $RES/step1b_diag --pairs $RES/step1b/ac.parquet \
    || { fail_marker "w_diag (bed b) failed"; exit 6; }
  [ -f "$RES/step1b_diag.json" ] || { fail_marker "diag_1b wrote no json"; exit 6; }
fi

# ---------------------------------------------------------------- 2. D4/D7 (cheap first) -------
STAGE="d7"
if [ ! -f "$RES/step1diag_D7.json" ]; then
  echo "=== [$(date -u +%FT%TZ)] STEP d7" >&2
  W_BED=c python $S/y_step1diag.py d7 --out $RES/step1diag --run $RES/step1c \
    || { fail_marker "d7 failed"; exit 7; }
fi

STAGE="d4"
if [ ! -f "$RES/step1diag_D4.json" ]; then
  echo "=== [$(date -u +%FT%TZ)] STEP d4" >&2
  W_BED=c python $S/y_step1diag.py d4 --out $RES/step1diag --run $RES/step1c \
    || { fail_marker "d4 failed"; exit 8; }
fi

STAGE="d4b"
if [ ! -f "$RES/step1diag_D4b.json" ]; then
  echo "=== [$(date -u +%FT%TZ)] STEP d4b" >&2
  W_BED=c python $S/y_step1diag.py d4b --out $RES/step1diag --run $RES/step1c \
    || { fail_marker "d4b failed"; exit 12; }
fi

STAGE="d4b_1b"
if [ ! -f "$RES/step1bdiag_D4b.json" ]; then
  echo "=== [$(date -u +%FT%TZ)] STEP d4b_1b" >&2
  if ! W_BED=b python $S/y_step1diag.py d4b --out $RES/step1bdiag --run $RES/step1b; then
    echo "!!! d4b on the 1b bed failed (non-fatal); continuing" >&2
  fi
fi

# ---------------------------------------------------------------- 3. D5 ------------------------
STAGE="d5"
if [ ! -f "$RES/step1diag_D5.json" ]; then
  echo "=== [$(date -u +%FT%TZ)] STEP d5" >&2
  W_BED=c python $S/y_step1diag.py d5 --out $RES/step1diag --min-events 50 --budget 1200 \
    || { fail_marker "d5 failed"; exit 9; }
fi

# ---------------------------------------------------------------- 4. D6 ------------------------
STAGE="d6"
if [ ! -f "$RES/step1diag_D6.json" ]; then
  echo "=== [$(date -u +%FT%TZ)] STEP d6" >&2
  W_BED=c python $S/y_step1diag.py d6 --out $RES/step1diag \
    || { fail_marker "d6 failed"; exit 10; }
fi

# D6 native-rate: raw mseed already on disk; NEVER downloads.  A failure here is NOT fatal --
# it is one sub-clause of one diagnostic, and the write-up reports its absence.
STAGE="d6native"
if [ ! -f "$RES/step1diag_D6native.json" ]; then
  echo "=== [$(date -u +%FT%TZ)] STEP d6native" >&2
  if ! W_BED=c python $S/y_step1diag.py d6native --out $RES/step1diag; then
    echo "!!! d6native failed (non-fatal); recording it and continuing" >&2
    python - <<'PY'
import json
json.dump(dict(status="d6native stage failed; see logs/step1diag_queue.log"),
          open("results/step1diag_D6native.json", "w"), indent=1)
PY
  fi
fi

# ---------------------------------------------------------------- 5. write-up ------------------
STAGE="write"
if [ ! -f "$RES/markers/STEP1DIAG_WRITE_DONE" ]; then
  echo "=== [$(date -u +%FT%TZ)] STEP write" >&2
  python $S/y_write_step1_results.py || { fail_marker "write-up failed"; exit 11; }
  touch "$RES/markers/STEP1DIAG_WRITE_DONE"
fi

# ---------------------------------------------------------------- done -------------------------
STAGE="done"
FINISHED=1
{
  echo "STEP1 DIAGNOSTICS D0-D7 -- DONE $(date -u +%FT%TZ)"
  echo "elapsed_s: $(( $(date +%s) - T0 ))"
  echo "stages: selftest diag_1c diag_1b d7 d4 d4b d4b_1b d5 d6 d6native write"
  echo "write-up: results/STEP1_RESULTS.md (v2); appendix results/STEP1_RESULTS_v1_superseded.md"
  echo "NOTE: no frozen verdict changed, no constant retuned, results/WA_GO NOT created."
  for f in step1c_diag.json step1b_diag.json step1diag_D4.json step1diag_D4b.json \
           step1bdiag_D4b.json step1diag_D5.json \
           step1diag_D6.json step1diag_D6native.json step1diag_D7.json; do
    echo "--- $f"; head -c 1200 "$RES/$f" 2>/dev/null; echo
  done
} > "$MDONE"
echo "=== [$(date -u +%FT%TZ)] y_step1diag_queue DONE ($(( $(date +%s) - T0 ))s)" >&2
exit 0
