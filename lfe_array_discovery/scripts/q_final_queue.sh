#!/usr/bin/env bash
# =================================================================================================
# THE FINAL SPECIFICATION of the north positive control -- ONE chained queue.
#
# Merlin's ruling, verbatim in results/merlin_ruling_step1cprime_2026-09-20.md and PLAN.md sec.12.
# Declared "the last change of any kind": if anything is found wanting after S5, STOP.
#
#   S0  self-test (flags, CUDA, the reused 1c' stage-1..5 products)
#   S7' the stage-7 seed set (every >= 2-node component of pairs >= T_pair = 8.0, top 2,000 by max
#       pair net_mad) -- 1c' stopped before stage 7, so it is built here, once, and hashed into v3
#   S1  quantisation A/B and re-alignment (CPU)   -> results/final/s1_quantisation.csv
#                                                    results/final/corrected_times_2005.parquet
#   S2  the corrected oracle ceiling (GPU)        -> results/final/oracle_reference_2005_v3.csv
#   S3  truth-seeded bootstrap-convergence test   -> results/final/s3/s3_verdict.json
#                                                    results/final/SELECTED_ORDER.txt
#   S4  results/freeze_step1c_prime_v3.json (NEW file, written BEFORE S5)
#   S5  THE ONE BLIND RUN + whole-chain site-shift null x3 + the frozen scorer
#
# MARKERS
#   results/FINAL_DONE    the verdict line (PASS / FAIL / FA_BROKEN) and the counts
#   results/FINAL_STOP    which pre-registered abandon fired, with its evidence lines
#   results/FINAL_FAILED  a crash or a signal (signal, stage, last 120 log lines)
#
# RULES OBSERVED
#   * NO FDSN waveform download anywhere (Step 6 is not this queue's);
#   * results/WA_GO is NEVER created, whatever the verdict;
#   * NO process is ever killed -- d_download_queue (T7 tail) and the parked x_step5_queue keep
#     running and share the GPU; single instance via flock; every stage guarded by its own output
#     file so a relaunch resumes instead of repeating GPU hours;
#   * TERM/INT/HUP/EXIT always leave one of the three markers.
# =================================================================================================
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
S=scripts
RES=results
RUN=$RES/step1c_prime
FIN=$RES/final
LOG=logs/final_queue.log
MDONE=$RES/FINAL_DONE
MSTOP=$RES/FINAL_STOP
MFAIL=$RES/FINAL_FAILED
mkdir -p logs locks "$FIN" "$RES/markers"

export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export PYTHONUNBUFFERED=1
export W_BED=c
unset W_MAX_DAYS

# ---- the 1c' repair flags, UNCHANGED (values recorded in freeze_step1c_prime.json) -------------
export W_HEALTH_VETO="$PWD/$RUN/station_hour_veto.csv"
export W_QUIET_HOURS_CSV="$PWD/$RUN/quiet_hours.csv"
export W_WIN_TRANSIENT=10
export W_DROP_ZERO_MEDIAN=1
export W_SEED_RULE=pairs2
export W_SEED_BUDGET=2000
export W_HANDOFF_MIN_SEG=2
export W_NEAR_SEED_MODE=either
# ---- the ONE new engine flag of the final spec (Merlin F5), recorded in v3 ----------------------
export W_SURV_ETS_ONLY=1

# ---------------------------------------------------------------- single instance ---------------
exec 9>locks/final.lock
if ! flock -n 9; then
  echo "=== [$(date -u +%FT%TZ)] another q_final_queue.sh holds locks/final.lock -- exiting" >&2
  exit 0
fi

STAGE="startup"
FINISHED=0

_tail() { echo "--- last 120 log lines ---"; tail -120 "$LOG" 2>/dev/null; }

_evidence() {
  for f in "$FIN"/S1_ABANDON.txt "$FIN"/S1_BELOW_ACCEPT.txt "$FIN"/S2_ABANDON.txt \
           "$FIN"/S3_ABANDON.txt; do
    [ -f "$f" ] && { echo "--- $(basename "$f") ---"; cat "$f"; }
  done
  for f in "$FIN"/s1_quantisation.json "$FIN"/oracle_reference_2005_v3.json \
           "$FIN"/s3/labelling.json "$FIN"/s3/s3_verdict.json; do
    [ -f "$f" ] && { echo "--- $f ---"; head -60 "$f"; }
  done
}

fail_marker() {
  { echo "step: $STAGE";
    echo "when: $(date -u +%FT%TZ)";
    echo "reason: $1";
    _tail; } > "$MFAIL"
}

stop_marker() {
  { echo "FINAL SPECIFICATION -- STOP (a pre-registered abandon fired)";
    echo "step: $STAGE";
    echo "when: $(date -u +%FT%TZ)";
    echo "reason: $1";
    _evidence;
    _tail; } > "$MSTOP"
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
  if [ "$FINISHED" != "1" ] && [ ! -f "$MFAIL" ] && [ ! -f "$MSTOP" ]; then
    fail_marker "queue exited with code $c during stage $STAGE"
  fi
}

trap 'on_signal TERM' TERM
trap 'on_signal INT'  INT
trap 'on_signal HUP'  HUP
trap on_exit EXIT

rm -f "$MDONE" "$MSTOP" "$MFAIL"
rm -f "$FIN"/S1_ABANDON.txt "$FIN"/S1_BELOW_ACCEPT.txt "$FIN"/S2_ABANDON.txt "$FIN"/S3_ABANDON.txt
T0=$(date +%s)
echo "=== [$(date -u +%FT%TZ)] q_final_queue START (pid $$)" >&2
echo "=== flags: veto=$W_HEALTH_VETO quiet=$W_QUIET_HOURS_CSV transient=$W_WIN_TRANSIENT" >&2
echo "=== seed=$W_SEED_RULE/$W_SEED_BUDGET handoff_min=$W_HANDOFF_MIN_SEG near=$W_NEAR_SEED_MODE ets_only=$W_SURV_ETS_ONLY" >&2

# step <name> <guard-file ("" = none)> <command...>
#   exit 20 = a pre-registered ABANDON fired      -> FINAL_STOP
#   exit 21 = BELOW_ACCEPT, report before Step 2  -> FINAL_STOP
#   any other non-zero                            -> FINAL_FAILED
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
  if [ $c -eq 20 ] || [ $c -eq 21 ]; then
    stop_marker "$([ $c -eq 20 ] && echo 'ABANDON' || echo 'BELOW_ACCEPT') condition fired in $name"
    echo "!!! STOP at $name (rc $c)" >&2
    exit $c
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

# ================================ S0. self-test =================================================
STAGE="selftest"
if [ ! -f "$RES/markers/FINAL_SELFTEST_OK" ]; then
  echo "=== [$(date -u +%FT%TZ)] STEP selftest" >&2
  python - <<'PY' || { STAGE=selftest; fail_marker "selftest failed"; exit 4; }
import os, sys
sys.path.insert(0, "scripts")
import torch
import w_common as W, w_grow as G, w_cluster as C                            # noqa: F401
import z_common as Z, z_score, z_oracle                                      # noqa: F401
import q_common as Q, q_chain, q_s1_quant, q_s2_oracle, q_s3_bootstrap       # noqa: F401
import q_s4_freeze, q_s5_blind                                               # noqa: F401
assert torch.cuda.is_available(), "no CUDA device"
f = W.repair_flags()
assert f["W_SEED_RULE"] == "pairs2" and f["W_SEED_BUDGET"] == 2000, f
assert f["W_HANDOFF_MIN_SEG"] == 2 and f["W_NEAR_SEED_MODE"] == "either", f
assert f["W_SURV_ETS_ONLY"] is True, f
assert W.BED_MODE == "c" and W.STEP1_YEAR == 2005, (W.BED_MODE, W.STEP1_YEAR)
# F5: the ETS window must be a STRICT subset of the cached days (24 of 44 on this bed)
bd, ed = len(W.bed_days()), len(W.ets_days())
assert 0 < ed < bd, (ed, bd)
# the 1c' stages 1-5 products that are REUSED, never recomputed
for p in ("results/step1c_prime/station_hour_veto.csv", "results/step1c_prime/quiet_hours.csv",
          "results/step1c_prime/ac.parquet", "results/step1c_prime/T.json",
          "results/freeze_step1c_prime.json", "results/freeze_step1c_prime_v2.json"):
    assert os.path.exists(p), p
import json
assert float(json.load(open("results/step1c_prime/T.json"))["T"]) == 8.0
assert not os.path.exists("results/WA_GO"), "results/WA_GO must not exist"
print("[selftest] flags", f)
print(f"[selftest] bed days {bd}, ETS days {ed}, sites {W.sites()[2]}")
PY
  touch "$RES/markers/FINAL_SELFTEST_OK"
fi

# ================================ S7'. the stage-7 seed set =====================================
step seeds "$RUN/cl_seeds.csv" python $S/w_cluster.py cluster --pairs $RUN/ac.parquet \
                                      --thr $RUN/T.json --out $RUN/cl

# ================================ S1. quantisation A/B (CPU) ====================================
step s1 "$FIN/S1_OK" python $S/q_s1_quant.py --out $FIN/s1_quantisation.csv \
                            --times-out $FIN/corrected_times_2005.parquet

# ================================ S2. corrected oracle ceiling (GPU) ============================
step s2 "$FIN/S2_OK" python $S/q_s2_oracle.py --dir $FIN/s2_oracle \
                            --times $FIN/corrected_times_2005.parquet \
                            --out $FIN/oracle_reference_2005_v3.csv

# ================================ S3. truth-seeded bootstrap test (GPU) =========================
step s3 "$FIN/SELECTED_ORDER.txt" python $S/q_s3_bootstrap.py --run $RUN --dir $FIN/s3 \
                            --oracle $FIN/oracle_reference_2005_v3.csv

# ================================ S4. the v3 freeze (before S5) =================================
step freeze "$RES/freeze_step1c_prime_v3.json" python $S/q_s4_freeze.py \
                            --out $RES/freeze_step1c_prime_v3.json --run $RUN --final $FIN

# ================================ S5. THE ONE BLIND RUN (GPU) ===================================
step s5 "$RES/final_score.json" python $S/q_s5_blind.py --run $RUN --dir $FIN/s5 \
                            --oracle $FIN/oracle_reference_2005_v3.csv --out $RES/final_score

# ================================ done ==========================================================
STAGE="done"
FINISHED=1
python - <<'PY' > "$MDONE"
import json, os, sys
import pandas as pd
sys.path.insert(0, "scripts")
os.environ["W_BED"] = "c"
import w_common as W


def j(p, d=None):
    try:
        return json.load(open(p))
    except Exception:
        return d or {}


S1 = j("results/final/s1_quantisation.json")
S2 = j("results/final/oracle_reference_2005_v3.json")
S3 = j("results/final/s3/s3_verdict.json")
FR = j("results/freeze_step1c_prime_v3.json")
SC = j("results/final_score.json")
print("FINAL SPECIFICATION of the north positive control -- DONE",
      pd.Timestamp.now("UTC").isoformat())
print(f"flags {W.repair_flags()}")
print(f"S1 quantisation {S1.get('verdict')}: 2012 coherent "
      f"{S1.get('n_families_2012_coherent')}/5 families (best {S1.get('best_2012_coherence')}); "
      f"2005 quantised median {S1.get('median_2005_quantised_coherence')}; lag spread "
      f"{S1.get('median_lag_spread_p90_p10_s')} s; POLARIS sites corrected "
      f"{S1.get('polaris_sites_corrected_ge2')}")
print(f"S2 corrected oracle: C_30 {S2.get('C_30')}/30, C_iso {S2.get('C_iso')}/10, fallback "
      f"{S2.get('fallback_active')}, bar {(S2.get('pass_bar') or {}).get('value')}; median recall "
      f"all/indep/loud {S2.get('median_recall_all')}/{S2.get('median_recall_indep')}/"
      f"{S2.get('median_recall_loud')}; median purity {S2.get('median_purity')}")
print(f"S3 bootstrap {S3.get('verdict')}: K {S3.get('K')} families from {S3.get('n_matched')} "
      f"true seeds (chance {S3.get('n_chance')}, {S3.get('chance_fraction')}); bar "
      f"{S3.get('bar')}; frozen {(S3.get('frozen') or {}).get('n_rediscovered')}, staged "
      f"{(S3.get('staged') or {}).get('n_rediscovered')} "
      f"({(S3.get('staged') or {}).get('n_isolated_rediscovered')} isolated); twins "
      f"{S3.get('twin_survivors')}/{S3.get('real_survivors_staged')} "
      f"({S3.get('twin_fraction')}) -> ORDER {S3.get('selected_order')}")
print(f"S4 freeze results/freeze_step1c_prime_v3.json (order {FR.get('selected_order')}, "
      f"ETS-day-only survival, seeds {(FR.get('seeds') or {}).get('n')})")
print(f"S5 blind: clusters {SC.get('n_clusters')}, survivors p0 {SC.get('real_survivors_p0')} / "
      f"final {SC.get('real_survivors_final')}; isolated rediscovered "
      f"{SC.get('n_isolated_rediscovered')}/10; all-30 {SC.get('n_hidden_rediscovered_all30')}/30; "
      f"split {SC.get('spurious_split_rate')}; null gated {SC.get('null_gated_per_realisation')}; "
      f"null survivors {SC.get('null_survivors_p0')} vs real {SC.get('real_survivors_p0')} "
      f"(fraction {SC.get('null_survivor_fraction')})")
print(f"VERDICT: {SC.get('verdict_line')}")
print(f"CONSEQUENCE: {SC.get('consequence')}")
print("NOTE: results/WA_GO was NOT created by this queue; no FDSN download was made; "
      "Step 6 (the 2004-07 re-fire) is NOT this queue's and runs only after a PASS.")
PY
echo "elapsed_s: $(( $(date +%s) - T0 ))" >> "$MDONE"
echo "=== [$(date -u +%FT%TZ)] q_final_queue COMPLETE in $(( $(date +%s) - T0 ))s" >&2
exit 0
