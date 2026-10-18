#!/usr/bin/env bash
# STEP 1c' -- THE ONE PERMITTED REPAIR of the north positive control (Merlin's prescription 1-10,
# results/merlin_review_step1_2026-09-20.md; PLAN.md sec.11 is binding).
#
#   1  D8 diagnostic (time base / physics / chain) over the 50 busiest kept family-hours
#   2  station-hour health veto over every station-hour of the bed
#   3  hour-matched, earthquake-screened quiet set (catalog query only -- NO waveform download)
#   4  FREEZE results/freeze_step1c_prime.json, then the engine fixes + window-transient veto and
#      the quiet/tremor autocorrelation
#   5  T_pair by the frozen FA rule with the mandatory leave-one-hour-out column
#   6  ORACLE CEILING and the 2005 scoring reference -- written BEFORE any discovery output is read
#   7  pair seeds (every >= 2-node component, ranked by max net_mad, budget 2000) + handoff
#   8  growth / survival / 2 iterations / dedupe / locate / Gates 0-2, exactly as frozen
#   9  whole-chain site-shift null x3, counting gated families AND survivors
#  10  score against the oracle reference and write the verdict line
#
# RULES OBSERVED: no FDSN waveform download anywhere; results/WA_GO is NEVER created; no process is
# ever killed (the T7 download queue, the parked Step-5 queue and any Arm-B GPU job keep running
# and share the GPU); single instance via flock; every stage is guarded by its own output file so a
# relaunch resumes; TERM/INT/HUP/EXIT always leave results/STEP1CP_DONE or results/STEP1CP_FAILED.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
S=scripts
RES=results
RUN=$RES/step1c_prime
LOG=logs/step1cp_queue.log
MDONE=$RES/STEP1CP_DONE
MFAIL=$RES/STEP1CP_FAILED
mkdir -p logs locks "$RUN" "$RES/markers"

export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export PYTHONUNBUFFERED=1
export W_BED=c
unset W_MAX_DAYS

# ---- the Step 1c' repair flags (their values are recorded in freeze_step1c_prime.json) ---------
export W_HEALTH_VETO="$PWD/$RUN/station_hour_veto.csv"
export W_QUIET_HOURS_CSV="$PWD/$RUN/quiet_hours.csv"
export W_WIN_TRANSIENT=10
export W_DROP_ZERO_MEDIAN=1
export W_SEED_RULE=pairs2
export W_SEED_BUDGET=2000
export W_HANDOFF_MIN_SEG=2
export W_NEAR_SEED_MODE=either

# ---------------------------------------------------------------- single instance ---------------
exec 9>locks/step1cp.lock
if ! flock -n 9; then
  echo "=== [$(date -u +%FT%TZ)] another z_step1cp_queue.sh holds locks/step1cp.lock -- exiting" >&2
  exit 0
fi

STAGE="startup"
FINISHED=0

fail_marker() {
  { echo "step: $STAGE";
    echo "when: $(date -u +%FT%TZ)";
    echo "reason: $1";
    if [ -f "$RUN/ABANDON_REASON.txt" ]; then echo "--- abandon reason ---"; cat "$RUN/ABANDON_REASON.txt"; fi
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

rm -f "$MDONE" "$MFAIL" "$RUN/ABANDON_REASON.txt"
T0=$(date +%s)
echo "=== [$(date -u +%FT%TZ)] z_step1cp_queue START (pid $$)" >&2
echo "=== flags: veto=$W_HEALTH_VETO quiet=$W_QUIET_HOURS_CSV transient=$W_WIN_TRANSIENT seed=$W_SEED_RULE/$W_SEED_BUDGET handoff_min=$W_HANDOFF_MIN_SEG near=$W_NEAR_SEED_MODE" >&2

# step <name> <guard-file ("" = none)> <command...>
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
  if [ $c -eq 20 ]; then
    fail_marker "ABANDON condition fired in $name"
    echo "!!! ABANDON at $name" >&2
    exit 20
  fi
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

# ================================ 0. self-test =================================================
STAGE="selftest"
if [ ! -f "$RES/markers/STEP1CP_SELFTEST_OK" ]; then
  echo "=== [$(date -u +%FT%TZ)] STEP selftest" >&2
  python - <<'PY' || { STAGE=selftest; fail_marker "selftest failed"; exit 4; }
import os, sys
sys.path.insert(0, "scripts")
import torch
import w_common as W, w_cluster as C, w_autocorr as A, w_grow as G   # noqa: F401
import z_common as Z, z_d8, z_health, z_quiet, z_tpair, z_oracle, z_score, z_null, z_checks  # noqa
assert torch.cuda.is_available(), "no CUDA device"
f = W.repair_flags()
assert f["W_SEED_RULE"] == "pairs2" and f["W_SEED_BUDGET"] == 2000, f
assert f["W_DROP_ZERO_MEDIAN"] and f["W_WIN_TRANSIENT"] == 10.0, f
assert f["W_HANDOFF_MIN_SEG"] == 2 and f["W_NEAR_SEED_MODE"] == "either", f
assert W.BED_MODE == "c" and W.STEP1_YEAR == 2005, (W.BED_MODE, W.STEP1_YEAR)
print("[selftest] flags", f)
print("[selftest] sites", W.sites()[2])
PY
  touch "$RES/markers/STEP1CP_SELFTEST_OK"
fi

# ================================ 1. D8 ========================================================
step d8        "$RUN/d8.json"                   python $S/z_d8.py --out $RUN/d8 --n-hours 50 --chain-hours 5
step d8_check  ""                               python $S/z_checks.py d8

# ================================ 2. station-hour health veto ==================================
step health       "$RUN/station_hour_veto.csv"  python $S/z_health.py --out $RUN/station_hour_veto.csv --workers 12
step health_check ""                            python $S/z_checks.py health

# ================================ 3. quiet set (catalog query, no download) ====================
step quiet        "$RUN/quiet_hours.csv"        python $S/z_quiet.py --out $RUN/quiet_hours.csv --veto $RUN/station_hour_veto.csv
step quiet_check  ""                            python $S/z_checks.py quiet

# ================================ 4. FREEZE, then the repaired autocorrelation =================
step freeze    "$RES/freeze_step1c_prime.json"  python $S/z_freeze.py --out $RES/freeze_step1c_prime.json --run $RUN
step ac_quiet  "$RUN/ac_quiet.parquet"          python $S/w_autocorr.py --hours quiet  --out $RUN/ac_quiet --t-record 4.0
step ac_tremor "$RUN/ac.parquet"                python $S/w_autocorr.py --hours tremor --out $RUN/ac      --t-record 5.0
step ac_check  ""                               python $S/z_checks.py ac

# ================================ 5. T_pair + leave-one-hour-out ===============================
step tpair       "$RUN/T.json"                  python $S/z_tpair.py --quiet-pairs $RUN/ac_quiet.parquet --out $RUN/T
step tpair_check ""                             python $S/z_checks.py tpair

# ================================ 6. oracle ceiling (before any discovery output) ==============
step oracle       "$RUN/oracle_reference_2005.csv" python $S/z_oracle.py --dir $RUN/oracle --out $RUN/oracle_reference_2005.csv
step oracle_check ""                            python $S/z_checks.py oracle

# ================================ 7. pair seeds and handoff ====================================
step cluster     "$RUN/cl_seeds.csv"            python $S/w_cluster.py cluster --pairs $RUN/ac.parquet --thr $RUN/T.json --out $RUN/cl
step handoff     "$RUN/bank0.npz"               python $S/w_grow.py handoff --clusters $RUN/cl_clusters.parquet --seeds $RUN/cl_seeds.csv --out $RUN/bank0.npz
step seeds_check ""                             python $S/z_checks.py seeds

# ================================ 8. growth / survival / iterations / locate ===================
step quiet0   "$RUN/p0_thresholds.csv" python $S/w_grow.py quiet    --bank $RUN/bank0.npz --out $RUN/p0
step grow0    "$RUN/p0_hits.parquet"   python $S/w_grow.py grow     --bank $RUN/bank0.npz --out $RUN/p0 --thr $RUN/p0_thresholds.csv
step scramble "$RUN/bank0scr.npz"      python $S/w_grow.py scramble --clusters $RUN/cl_clusters.parquet --seeds $RUN/cl_seeds.csv --out $RUN/bank0scr.npz
step grow0scr "$RUN/p0scr_hits.parquet" python $S/w_grow.py grow    --bank $RUN/bank0scr.npz --out $RUN/p0scr --thr $RUN/p0_thresholds.csv
step survive0 "$RUN/p0_survival.csv"   python $S/w_grow.py survive  --bank $RUN/bank0.npz --out $RUN/p0 --seeds $RUN/cl_seeds.csv \
                                              --scramble-hits $RUN/p0scr_hits.parquet --scramble-bank $RUN/bank0scr.npz

for it in 1 2; do
  prev=$((it-1))
  step iterate$it "$RUN/bank$it.npz"          python $S/w_iterate.py --bank $RUN/bank$prev.npz --indep $RUN/p${prev}_indep.parquet \
                                                     --survivors $RUN/p0_survival.csv --out-bank $RUN/bank$it.npz --out-diag $RUN/iter${it}_diag.csv
  step quiet$it   "$RUN/p${it}_thresholds.csv" python $S/w_grow.py quiet --bank $RUN/bank$it.npz --out $RUN/p$it
  step grow$it    "$RUN/p${it}_hits.parquet"   python $S/w_grow.py grow  --bank $RUN/bank$it.npz --out $RUN/p$it --thr $RUN/p${it}_thresholds.csv
  step surv$it    "$RUN/p${it}_survival.csv"   python $S/w_grow.py survive --bank $RUN/bank$it.npz --out $RUN/p$it --seeds $RUN/cl_seeds.csv
done

step dedupe "$RUN/dedupe.csv"          python $S/w_dedupe.py --bank $RUN/bank2.npz --indep $RUN/p2_indep.parquet --out $RUN/dedupe.csv
step locate "$RUN/lg.json"             python $S/w_locate_gates.py --bank $RUN/bank2.npz --indep $RUN/p2_indep.parquet \
                                              --dedupe $RUN/dedupe.csv --survival $RUN/p2_survival.csv --out $RUN/lg

# ================================ 9. whole-chain site-shift null x3 ============================
step nulls "$RUN/null_summary.csv"     python $S/z_null.py --out-root $RUN --n 3

# ================================ 10. score ====================================================
step score "$RES/step1c_prime_score.json" python $S/z_score.py --run $RUN --oracle $RUN/oracle_reference_2005.csv \
                                                 --out $RES/step1c_prime_score

# ================================ done =========================================================
STAGE="done"
FINISHED=1
python - <<'PY' > "$MDONE"
import json, os, sys
import pandas as pd
sys.path.insert(0, "scripts")
os.environ["W_BED"] = "c"
import w_common as W
RUN = "results/step1c_prime"


def j(p, k=None, d=None):
    try:
        v = json.load(open(p))
        return v if k is None else v.get(k, d)
    except Exception:
        return d


def n(p, col=None, val=None):
    try:
        x = pd.read_csv(p) if p.endswith("csv") else pd.read_parquet(p)
        return int(len(x) if col is None else (x[col] == val).sum())
    except Exception:
        return -1


S = j("results/step1c_prime_score.json") or {}
O = j(f"{RUN}/oracle_reference_2005.json") or {}
D = j(f"{RUN}/d8.json") or {}
H = j(f"{RUN}/station_hour_veto.json") or {}
Q = j(f"{RUN}/quiet_hours.json") or {}
T = j(f"{RUN}/T.json") or {}
print("STEP 1c' -- the one permitted repair -- DONE", pd.Timestamp.now("UTC").isoformat())
print(f"flags {W.repair_flags()}")
print(f"D8 {D.get('verdict')}: peaks within 3 s at {D.get('n_sites_within_3s')} sites, control "
      f"{D.get('control_excess_median')}, chain ratio {D.get('chain_ratio_median')}")
print(f"health veto {H.get('n_vetoed')}/{H.get('n_station_hours')} station-hours; SHVB flagged "
      f"{H.get('shvb_veto_frac_on_bad_doys')} on doy 241-251/255-261")
print(f"quiet set {Q.get('n_quiet_hours')} hours, {Q.get('clock_hours_covered')}/24 clock hours, "
      f"{Q.get('n_hours_eq_vetoed')} earthquake-vetoed")
print(f"T_pair {T.get('T')} ({T.get('binds')}); LOO range {T.get('loo_range')}")
print(f"oracle C_iso {O.get('C_iso')}/10, C_30 {O.get('C_30')}/30, fallback "
      f"{O.get('fallback_active')}, bar {(O.get('pass_bar') or {}).get('value')}")
print(f"pairs {n(f'{RUN}/ac.parquet')}; seeds {n(f'{RUN}/cl_seeds.csv')}; survivors p0 "
      f"{n(f'{RUN}/p0_survival.csv','survives',True)}; final "
      f"{n(f'{RUN}/p2_survival.csv','survives',True)}; located {j(f'{RUN}/lg.json','n_located')}; "
      f"gated {j(f'{RUN}/lg.json','n_gate_pass')}")
print(f"nulls gated {S.get('null_gated_per_realisation')}; null survivors "
      f"{S.get('null_survivors_p0')} vs real {S.get('real_survivors_p0')} "
      f"(fraction {S.get('null_survivor_fraction')})")
print(f"isolated rediscovered {S.get('n_isolated_rediscovered')}/10; all-30 "
      f"{S.get('n_hidden_rediscovered_all30')}/30; split {S.get('spurious_split_rate')}")
print(f"VERDICT: {S.get('verdict_line')}")
print(f"CONSEQUENCE: {S.get('consequence')}")
print("NOTE: results/WA_GO was NOT created by this queue; freeze results/freeze_step1c_prime.json")
PY
echo "elapsed_s: $(( $(date +%s) - T0 ))" >> "$MDONE"
echo "=== [$(date -u +%FT%TZ)] z_step1cp_queue COMPLETE in $(( $(date +%s) - T0 ))s" >&2
exit 0
