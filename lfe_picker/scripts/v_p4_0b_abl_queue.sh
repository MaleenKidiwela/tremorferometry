#!/usr/bin/env bash
# P4-0b SENSITIVITY ABLATION -- NOT THE FROZEN RESULT.
# The frozen given-seed verdict (0/18 assembled) stands; this queue only brackets how much of it
# is the >= 20-independent-hit survival floor, set for a multi-year record and applied to 203 days.
# Launch:  nohup bash scripts/v_p4_0b_abl_queue.sh > logs/p4_0b_abl_queue.log 2>&1 &
# Markers: results/p4/P4_0B_ABL_DONE | results/p4/P4_0B_ABL_FAILED
# Writes ONLY under results/p4/p4_0b_abl/; results/p4/p4_0b/ is opened READ-ONLY.
# Shares the GPU with the running frozen free-arm queue (v_engine.wait_gpu).
set -u
cd /home/jovyan/tremorferometry/lfe_picker
export PYTHONUNBUFFERED=1
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=6
export MKL_NUM_THREADS=6
A=results/p4/p4_0b_abl
LOG=logs/p4_0b_abl_queue.log
mkdir -p "$A" logs

ts(){ date -u +'%Y-%m-%dT%H:%M:%SZ'; }
STEP=""
fail(){ local c=$1
  { echo "SENSITIVITY ABLATION -- not the frozen result"; echo "step: $STEP"; echo "exit_code: $c";
    echo "--- last 120 log lines ---"; tail -120 "$LOG"; } > results/p4/P4_0B_ABL_FAILED
  echo "[$(ts)] FAILED at $STEP (exit $c)"; exit "$c"; }
run(){ STEP="$1"; local guard="$2"; shift 2
  if [ "$guard" != "-" ] && [ -e "$guard" ]; then echo "[$(ts)] SKIP  $STEP"; return 0; fi
  echo "[$(ts)] START $STEP :: $*"; "$@" || fail $?; echo "[$(ts)] DONE  $STEP"; }

T0=$(date +%s)
rm -f results/p4/P4_0B_ABL_FAILED
echo "[$(ts)] P4-0b ABLATION queue start, pid $$  -- SENSITIVITY ABLATION, NOT THE FROZEN RESULT"

run prereg "$A/../p4_0b_abl_prereg.json" python scripts/v_abl.py prereg

# ---- A3 / A4: no GPU, land the diagnostics first -----------------------------------------------
run A3_diag     "$A/A3_diagnostics.json"  python scripts/v_abl.py diag
run A4_refcheck "$A/A4_reference.json"    python scripts/v_abl.py refcheck

# ---- A1 (floor >= 8, >= 3 days) and A2 (floor >= 5, >= 2 days) ----------------------------------
for L in A1 A2; do
  run ${L}_surv0    "$A/${L}_run0_survival.csv"        python scripts/v_abl.py survive --label $L --level 0
  run ${L}_it1bank  "$A/${L}_bank_it1.npz"             python scripts/v_abl.py iterate --label $L --level 1
  run ${L}_it1quiet "$A/${L}_run_it1_thresholds.csv"   python scripts/v_run.py quiet --bank "$A/${L}_bank_it1.npz" --out "$A/${L}_run_it1"
  run ${L}_it1grow  "$A/${L}_run_it1_hits.parquet"     python scripts/v_run.py grow  --bank "$A/${L}_bank_it1.npz" --out "$A/${L}_run_it1"
  run ${L}_it1surv  "$A/${L}_run_it1_survival.csv"     python scripts/v_abl.py survive --label $L --level 1
  run ${L}_it2bank  "$A/${L}_bank_it2.npz"             python scripts/v_abl.py iterate --label $L --level 2
  run ${L}_it2quiet "$A/${L}_run_it2_thresholds.csv"   python scripts/v_run.py quiet --bank "$A/${L}_bank_it2.npz" --out "$A/${L}_run_it2"
  run ${L}_it2grow  "$A/${L}_run_it2_hits.parquet"     python scripts/v_run.py grow  --bank "$A/${L}_bank_it2.npz" --out "$A/${L}_run_it2"
  run ${L}_it2surv  "$A/${L}_run_it2_survival.csv"     python scripts/v_abl.py survive --label $L --level 2
  run ${L}_dedupe   "$A/${L}_dedupe.csv"               python scripts/v_dedupe.py --bank "$A/${L}_bank_it2.npz" --indep "$A/${L}_run_it2_indep.parquet" --out "$A/${L}_dedupe.csv"
  run ${L}_score    "$A/${L}_score.json"               python scripts/v_abl.py score --label $L
done

run report - python scripts/v_abl.py report

T1=$(date +%s)
python - "$A" $((T1-T0)) > results/p4/P4_0B_ABL_DONE <<'PY'
import json, sys, os
import pandas as pd
A, secs = sys.argv[1], int(sys.argv[2])
def J(p):
    try: return json.load(open(f"{A}/{p}"))
    except Exception: return {}
L = ["SENSITIVITY ABLATION -- NOT THE FROZEN RESULT",
     "FROZEN VERDICT (unchanged): given-seed arm 0/18 families assembled",
     f"ablation runtime {secs/3600:.2f} h"]
for lab in ("A1", "A2"):
    s = J(f"{lab}_score.json")
    L.append(f"{lab} floor>={s.get('clauses',{}).get('floor','?')} days>="
             f"{s.get('clauses',{}).get('min_days','?')}: assembled by level "
             f"{s.get('assembled_by_level','n/a')} of {s.get('n_eligible','?')}; "
             f"near-misses {len(s.get('near_miss',[]))}")
d3 = J("A3_diagnostics.json")
L.append(f"A3 pass-0 templates {d3.get('n_pass0_templates','?')}: median indep hits "
         f"{d3.get('median_n_indep','?')} (p75 {d3.get('p75_n_indep','?')}, max "
         f"{d3.get('max_n_indep','?')}), median own-family purity "
         f"{d3.get('median_purity_own_all','?')}, Spearman(purity, n_indep) "
         f"{d3.get('spearman_purity_vs_hits','?')}")
for b in d3.get("purity_by_hitcount", []):
    L.append(f"   n_indep {b['n_indep_bin']:>6s}: {b['n_templates']:3d} tpl, median purity_own "
             f"{b['median_purity_own']}")
d4 = J("A4_reference.json")
for k, v in d4.get("by_threshold", {}).items():
    L.append(f"A4 reference @ {k}: median recall {v['median_recall']:.3f}, purity "
             f"{v['median_purity']:.3f}, n_indep {v['median_n_indep']:.0f}")
L.append(f"summary json : {A}/P4_0B_ABL_SUMMARY.json")
L.append(f"pre-registration: results/p4/p4_0b_abl_prereg.json")
print("\n".join(L[:20]))
PY
echo "[$(ts)] ABLATION COMPLETE in $(( (T1-T0)/60 )) min"
cat results/p4/P4_0B_ABL_DONE
