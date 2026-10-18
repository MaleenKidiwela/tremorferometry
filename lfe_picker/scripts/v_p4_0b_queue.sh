#!/usr/bin/env bash
# P4-0b -- the TRUE Arm B kill test on the PQ-era north.  ONE unattended queue.
# Launch:  nohup bash scripts/v_p4_0b_queue.sh > logs/p4_0b_queue.log 2>&1 &
# Markers: results/p4/P4_0B_DONE  |  results/p4/P4_0B_FAILED
# Every step is resumable and idempotent; re-running the queue continues where it stopped.
set -u
cd /home/jovyan/tremorferometry/lfe_picker
export PYTHONUNBUFFERED=1
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
D=results/p4/p4_0b
LOG=logs/p4_0b_queue.log
mkdir -p "$D" results/p4 logs

ts(){ date -u +'%Y-%m-%dT%H:%M:%SZ'; }
STEP=""
fail(){
  local code=$1
  { echo "step: $STEP"; echo "exit_code: $code"; echo "--- last 120 log lines ---";
    tail -120 "$LOG"; } > results/p4/P4_0B_FAILED
  echo "[$(ts)] FAILED at $STEP (exit $code)"
  exit "$code"
}
run(){  # run <step-name> <guard-file-or-> <command...>
  STEP="$1"; local guard="$2"; shift 2
  if [ "$guard" != "-" ] && [ -e "$guard" ]; then
    echo "[$(ts)] SKIP  $STEP (have $guard)"; return 0
  fi
  echo "[$(ts)] START $STEP :: $*"
  "$@" || fail $?
  echo "[$(ts)] DONE  $STEP"
}

T0=$(date +%s)
rm -f results/p4/P4_0B_FAILED
echo "[$(ts)] P4-0b queue start, pid $$"

# ---------------------------------------------------------------- seeds ------------------------
run seeds_given  "$D/seeds_given.parquet" python scripts/v_seeds.py --arm given
run hp_seed_band "$D/hp_seed_z.json"      python scripts/v_seeds.py --arm hpband

# ---------------------------------------------------------------- ARM (i): GIVEN seeds ---------
run bank_given0    "$D/bank_given0.npz"  python scripts/v_bank.py pass0 --seeds "$D/seeds_given.parquet" --out "$D/bank_given0.npz"
run given0_quiet   "$D/run_given0_thresholds.csv" python scripts/v_run.py quiet --bank "$D/bank_given0.npz" --out "$D/run_given0"
run given0_grow    "$D/run_given0_hits.parquet"   python scripts/v_run.py grow  --bank "$D/bank_given0.npz" --out "$D/run_given0"

run bank_given_nm  "$D/bank_given_nm.npz" python scripts/v_bank.py nmoveout --seeds "$D/seeds_given.parquet" --out "$D/bank_given_nm.npz"
run given_nm_quiet "$D/run_given_nm_thresholds.csv" python scripts/v_run.py quiet --bank "$D/bank_given_nm.npz" --out "$D/run_given_nm"
run given_nm_grow  "$D/run_given_nm_hits.parquet"   python scripts/v_run.py grow  --bank "$D/bank_given_nm.npz" --out "$D/run_given_nm"

run given0_survive "$D/run_given0_survival.csv" python scripts/v_run.py survive --bank "$D/bank_given0.npz" --out "$D/run_given0" --seeds "$D/seeds_given.parquet" --nmoveout-out "$D/run_given_nm" --nmoveout-bank "$D/bank_given_nm.npz"

run given_it1_bank "$D/bank_given_it1.npz" python scripts/v_iterate.py iterate --bank "$D/bank_given0.npz" --indep "$D/run_given0_indep.parquet" --survivors "$D/run_given0_survival.csv" --out-bank "$D/bank_given_it1.npz" --out-diag "$D/it1_given_diag.csv"
run given_it1_quiet "$D/run_given_it1_thresholds.csv" python scripts/v_run.py quiet --bank "$D/bank_given_it1.npz" --out "$D/run_given_it1"
run given_it1_grow  "$D/run_given_it1_hits.parquet"   python scripts/v_run.py grow  --bank "$D/bank_given_it1.npz" --out "$D/run_given_it1"
run given_it1_surv  "$D/run_given_it1_survival.csv"   python scripts/v_run.py survive --bank "$D/bank_given_it1.npz" --out "$D/run_given_it1" --seeds "$D/seeds_given.parquet"

run given_it2_bank "$D/bank_given_it2.npz" python scripts/v_iterate.py iterate --bank "$D/bank_given_it1.npz" --indep "$D/run_given_it1_indep.parquet" --survivors "$D/run_given_it1_survival.csv" --use-all --out-bank "$D/bank_given_it2.npz" --out-diag "$D/it2_given_diag.csv"
run given_it2_quiet "$D/run_given_it2_thresholds.csv" python scripts/v_run.py quiet --bank "$D/bank_given_it2.npz" --out "$D/run_given_it2"
run given_it2_grow  "$D/run_given_it2_hits.parquet"   python scripts/v_run.py grow  --bank "$D/bank_given_it2.npz" --out "$D/run_given_it2"
run given_it2_surv  "$D/run_given_it2_survival.csv"   python scripts/v_run.py survive --bank "$D/bank_given_it2.npz" --out "$D/run_given_it2" --seeds "$D/seeds_given.parquet"

run given_dedupe "$D/dedupe_given.csv" python scripts/v_dedupe.py --bank "$D/bank_given_it2.npz" --indep "$D/run_given_it2_indep.parquet" --out "$D/dedupe_given.csv"
run given_screen "$D/screen_given.csv" python scripts/v_iterate.py screen --bank "$D/bank_given_it2.npz" --indep "$D/run_given_it2_indep.parquet" --out "$D/screen_given.csv"

# ---------------------------------------------------------------- arm (i) REFERENCE ------------
run bank_ref  "$D/bank_ref.npz"            python scripts/v_bank.py ref --out "$D/bank_ref.npz"
run ref_quiet "$D/run_ref_thresholds.csv"  python scripts/v_run.py quiet --bank "$D/bank_ref.npz" --out "$D/run_ref"
run ref_grow  "$D/run_ref_hits.parquet"    python scripts/v_run.py grow  --bank "$D/bank_ref.npz" --out "$D/run_ref"
run ref_surv  "$D/run_ref_survival.csv"    python scripts/v_run.py survive --bank "$D/bank_ref.npz" --out "$D/run_ref" --seeds "$D/seeds_given.parquet"

# arm (i) is the DECISIVE arm: score it as soon as it is complete, before the free arm runs.
STEP=score_given_partial; python scripts/v_score.py > "$D/score_partial.txt" 2>&1 || true

# ---------------------------------------------------------------- ARM (ii): TEMPLATE-FREE ------
run seeds_free  "$D/seeds_free.parquet" python scripts/v_seeds.py --arm free
run bank_free0  "$D/bank_free0.npz"     python scripts/v_bank.py pass0 --seeds "$D/seeds_free.parquet" --out "$D/bank_free0.npz"
run free0_quiet "$D/run_free0_thresholds.csv" python scripts/v_run.py quiet --bank "$D/bank_free0.npz" --out "$D/run_free0"
run free0_grow  "$D/run_free0_hits.parquet"   python scripts/v_run.py grow  --bank "$D/bank_free0.npz" --out "$D/run_free0"

run bank_free_nm  "$D/bank_free_nm.npz" python scripts/v_bank.py nmoveout --seeds "$D/seeds_free.parquet" --out "$D/bank_free_nm.npz"
run free_nm_quiet "$D/run_free_nm_thresholds.csv" python scripts/v_run.py quiet --bank "$D/bank_free_nm.npz" --out "$D/run_free_nm"
run free_nm_grow  "$D/run_free_nm_hits.parquet"   python scripts/v_run.py grow  --bank "$D/bank_free_nm.npz" --out "$D/run_free_nm"

run free0_survive "$D/run_free0_survival.csv" python scripts/v_run.py survive --bank "$D/bank_free0.npz" --out "$D/run_free0" --seeds "$D/seeds_free.parquet" --nmoveout-out "$D/run_free_nm" --nmoveout-bank "$D/bank_free_nm.npz"

run free_it1_bank "$D/bank_free_it1.npz" python scripts/v_iterate.py iterate --bank "$D/bank_free0.npz" --indep "$D/run_free0_indep.parquet" --survivors "$D/run_free0_survival.csv" --out-bank "$D/bank_free_it1.npz" --out-diag "$D/it1_free_diag.csv"
run free_it1_quiet "$D/run_free_it1_thresholds.csv" python scripts/v_run.py quiet --bank "$D/bank_free_it1.npz" --out "$D/run_free_it1"
run free_it1_grow  "$D/run_free_it1_hits.parquet"   python scripts/v_run.py grow  --bank "$D/bank_free_it1.npz" --out "$D/run_free_it1"
run free_it1_surv  "$D/run_free_it1_survival.csv"   python scripts/v_run.py survive --bank "$D/bank_free_it1.npz" --out "$D/run_free_it1" --seeds "$D/seeds_free.parquet"

run free_it2_bank "$D/bank_free_it2.npz" python scripts/v_iterate.py iterate --bank "$D/bank_free_it1.npz" --indep "$D/run_free_it1_indep.parquet" --survivors "$D/run_free_it1_survival.csv" --use-all --out-bank "$D/bank_free_it2.npz" --out-diag "$D/it2_free_diag.csv"
run free_it2_quiet "$D/run_free_it2_thresholds.csv" python scripts/v_run.py quiet --bank "$D/bank_free_it2.npz" --out "$D/run_free_it2"
run free_it2_grow  "$D/run_free_it2_hits.parquet"   python scripts/v_run.py grow  --bank "$D/bank_free_it2.npz" --out "$D/run_free_it2"
run free_it2_surv  "$D/run_free_it2_survival.csv"   python scripts/v_run.py survive --bank "$D/bank_free_it2.npz" --out "$D/run_free_it2" --seeds "$D/seeds_free.parquet"

run free_dedupe "$D/dedupe_free.csv" python scripts/v_dedupe.py --bank "$D/bank_free_it2.npz" --indep "$D/run_free_it2_indep.parquet" --out "$D/dedupe_free.csv"
run free_screen "$D/screen_free.csv" python scripts/v_iterate.py screen --bank "$D/bank_free_it2.npz" --indep "$D/run_free_it2_indep.parquet" --out "$D/screen_free.csv"

# ---------------------------------------------------------------- scoring ----------------------
run score - python scripts/v_score.py

# ---------------------------------------------------------------- done -------------------------
T1=$(date +%s)
python - "$D" $((T1-T0)) > results/p4/P4_0B_DONE <<'PY'
import json, os, sys
import pandas as pd
D, secs = sys.argv[1], int(sys.argv[2])
def n(p, col=None):
    try:
        if p.endswith(".parquet"):
            d = pd.read_parquet(p)
        else:
            d = pd.read_csv(p)
        return int(d[col].sum()) if col else len(d)
    except Exception:
        return "n/a"
S = {}
try:
    S = json.load(open(f"{D}/P4_0B_SCORE.json"))
except Exception:
    pass
L = []
L.append(f"P4-0b COMPLETE  runtime {secs/3600:.2f} h")
L.append(f"eligible kept families           : {S.get('n_eligible','n/a')}")
L.append(f"seeds  given / free              : {n(D+'/seeds_given.parquet')} / {n(D+'/seeds_free.parquet')}")
L.append(f"pass-0 survivors given / free    : {n(D+'/run_given0_survival.csv','survives')} / {n(D+'/run_free0_survival.csv','survives')}")
L.append(f"iterated templates given / free  : {n(D+'/run_given_it2_survival.csv')} / {n(D+'/run_free_it2_survival.csv')}")
L.append(f"distinct after dedupe g / f      : {n(D+'/dedupe_given.csv','kept')} / {n(D+'/dedupe_free.csv','kept')}")
L.append(f"screen passes given / free       : {n(D+'/screen_given.csv','passes')} / {n(D+'/screen_free.csv','passes')}")
a1 = S.get("arm_i_given", {}); a2 = S.get("arm_ii_free", {})
L.append(f"ARM (i)  assembled               : {a1.get('n_assembled','n/a')} of {a1.get('n_testable','n/a')} testable"
         f"  (no-reference {a1.get('n_no_reference','n/a')})   PASS={a1.get('pass_i','n/a')}")
L.append(f"ARM (ii) rediscovered families   : {a2.get('n_rediscovered_families','n/a')}   PASS={a2.get('pass_ii','n/a')}")
L.append(f"burst convergence (>50% near seed): {S.get('burst_convergence',{}).get('frac_families_gt50pct_near_seed','n/a')}")
hz = S.get("hp_seed_z_band", {})
L.append(f"TRUE hp-event seed-z band        : median {hz.get('z_med_percentiles',{}).get('50','n/a')}"
         f"  in-band frac {hz.get('frac_in_band','n/a')}  below-4 {hz.get('frac_zmed_below_4','n/a')}")
fs = S.get("free_seed_stats", {})
L.append(f"free seed histogram              : max/hour {fs.get('max_seeds_per_hour','n/a')},"
         f" top-10-hours {fs.get('top10_hours_fraction','n/a')}")
L.append(f"VERDICT                          : {S.get('verdict','n/a')}")
L.append(f"score json                       : {D}/P4_0B_SCORE.json")
L.append(f"freeze                           : results/p4/p4_0b_freeze.json")
print("\n".join(L[:20]))
PY
echo "[$(ts)] P4-0b queue COMPLETE in $(( (T1-T0)/60 )) min"
cat results/p4/P4_0B_DONE
