#!/usr/bin/env bash
# STEP 1c -- the 2005-09 POLARIS bed: the positive control with the geometry the WA arrays offer.
#
# PHASE A (no GPU, starts immediately): download + frozen-chain cache -> array-energy hour flag ->
#          full freeze (availability / quiet-hour / candidate-hour lists).
# PHASE B: wait for Step 1b to release the GPU (results/STEP1B_DONE or STEP1B_FAILED, poll 300 s).
# PHASE C (GPU): quiet T -> autocorrelation -> clusters -> handoff -> growth -> survival +
#          template-scramble null -> 2 iterations -> dedupe -> locate/gates -> 3 site-shift nulls
#          -> scoring vs the 10 isolated hidden families -> engine diagnostics.
#
# NO engine constant is retuned: window, hop, thresholds, survival, dedupe and gate rules are
# byte-identical to freeze_step1.json.  Only the BED changes (year, window, hours, quiet days,
# sites, withheld), and the SPEC was frozen before launch in results/freeze_step1c.json.
# Markers: results/STEP1C_DONE / results/STEP1C_FAILED.  Log: logs/step1c_queue.log.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
S=scripts
O=results/step1c
LOG=logs/step1c_queue.log
mkdir -p "$O" logs
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export W_BED=c
unset W_MAX_DAYS
rm -f results/STEP1C_DONE results/STEP1C_FAILED

step() {
  local name="$1"; shift
  echo "=== [$(date -u +%FT%TZ)] STEP $name" >&2
  "$@"
  local c=$?
  if [ $c -ne 0 ]; then
    { echo "step: $name"; echo "exit_code: $c"; echo "cmd: $*"; echo "--- last 120 log lines ---";
      tail -120 "$LOG"; } > results/STEP1C_FAILED
    echo "!!! FAILED at $name ($c)" >&2
    exit $c
  fi
}
T0=$(date +%s)

# ================================ PHASE A -- no GPU, runs now ================================
step spec_check python - <<'PY'
import json, os, sys
sys.path.insert(0, "scripts")
os.environ["W_BED"] = "c"
import w_common as W
F = json.load(open("results/freeze_step1c.json"))
bad = [os.path.basename(p) for p in W.w_scripts()
       if F["sha256"]["w_scripts"].get(os.path.basename(p)) != W.sha256_file(p)]
print("[spec] scripts changed since the pre-launch freeze:", bad or "none")
assert not bad, f"re-freeze before running: {bad}"
assert os.path.exists(W.STEP1C_BED), "step1c_bed.json missing"
print("[spec] bed:", F["bed_config"]["sites"], "withheld", F["bed_config"]["withheld"])
PY

step fetch_bed  python $S/w_fetch_bed2005.py --workers 4 --cache-workers 10
step energy     python $S/w_energy_hours.py
step freeze_full python $S/w_freeze.py --out results/freeze_step1c_full.json --tag step1c \
  --deviation "STEP 1c, coordinator-authorised deviation for Merlin's audit. The 2012 window has no dense array (POLARIS BC ran 2002-2007); the 2005-09 ETS is the only north data with array-class geometry. Window from the literature (WTD catalogue not retrievable: the Slow Earthquake Database serves catalogues only through php/sqlIntrTrem.php, whose parameter shape returns [] even for a catalogue inside its own date range). HOUR selection from the ARRAY-ENERGY FLAG, which Merlin's Q3 sanctions for hour selection only; quiet days are literature-quiet. NO Bostock event time is used for window, hour or day selection - only inside the scorer. NO engine constant retuned."

# ================================ PHASE B -- wait for the GPU ================================
echo "=== [$(date -u +%FT%TZ)] PHASE B: waiting for Step 1b to finish (poll 300 s)" >&2
while [ ! -f results/STEP1B_DONE ] && [ ! -f results/STEP1B_FAILED ]; do sleep 300; done
echo "=== [$(date -u +%FT%TZ)] Step 1b released the GPU ($(ls results/STEP1B_DONE results/STEP1B_FAILED 2>/dev/null | tr '\n' ' '))" >&2

# ================================ PHASE C -- GPU ================================
step ac_quiet   python $S/w_autocorr.py --hours quiet  --out $O/ac_quiet --t-record 4.0
step threshold  python $S/w_cluster.py  threshold --quiet-pairs $O/ac_quiet.parquet --out $O/T
step ac_tremor  python $S/w_autocorr.py --hours tremor --out $O/ac --t-record 5.0
step cluster    python $S/w_cluster.py  cluster --pairs $O/ac.parquet --thr $O/T.json --out $O/cl

step handoff    python $S/w_grow.py handoff  --clusters $O/cl_clusters.parquet --seeds $O/cl_seeds.csv --out $O/bank0.npz
step quiet0     python $S/w_grow.py quiet    --bank $O/bank0.npz --out $O/p0
step grow0      python $S/w_grow.py grow     --bank $O/bank0.npz --out $O/p0 --thr $O/p0_thresholds.csv
step scramble   python $S/w_grow.py scramble --clusters $O/cl_clusters.parquet --seeds $O/cl_seeds.csv --out $O/bank0scr.npz
step grow0scr   python $S/w_grow.py grow     --bank $O/bank0scr.npz --out $O/p0scr --thr $O/p0_thresholds.csv
step survive0   python $S/w_grow.py survive  --bank $O/bank0.npz --out $O/p0 --seeds $O/cl_seeds.csv \
                       --scramble-hits $O/p0scr_hits.parquet --scramble-bank $O/bank0scr.npz

for it in 1 2; do
  prev=$((it-1))
  step iterate$it python $S/w_iterate.py --bank $O/bank$prev.npz --indep $O/p${prev}_indep.parquet \
                         --survivors $O/p0_survival.csv --out-bank $O/bank$it.npz --out-diag $O/iter${it}_diag.csv
  step quiet$it   python $S/w_grow.py quiet --bank $O/bank$it.npz --out $O/p$it
  step grow$it    python $S/w_grow.py grow  --bank $O/bank$it.npz --out $O/p$it --thr $O/p${it}_thresholds.csv
  step surv$it    python $S/w_grow.py survive --bank $O/bank$it.npz --out $O/p$it --seeds $O/cl_seeds.csv
done

step dedupe     python $S/w_dedupe.py --bank $O/bank2.npz --indep $O/p2_indep.parquet --out $O/dedupe.csv
step locate     python $S/w_locate_gates.py --bank $O/bank2.npz --indep $O/p2_indep.parquet \
                       --dedupe $O/dedupe.csv --survival $O/p2_survival.csv --out $O/lg
step nulls      python $S/w_null_siteshift.py --out-root $O --n 3
step score      python $S/w_score_step1.py --clusters $O/lg_clusters.parquet --gates $O/lg_gates.csv \
                       --nulls $O/null_summary.csv --out results/step1c_score --year 2005
step diag       python $S/w_diag.py --out results/step1c_diag --pairs $O/ac.parquet --year 2005

python - <<'PY' > results/STEP1C_DONE
import json, os, sys
import pandas as pd
sys.path.insert(0, "scripts")
os.environ["W_BED"] = "c"
import w_common as W
O = "results/step1c"


def j(p, k, d=None):
    try:
        return json.load(open(p)).get(k, d)
    except Exception:
        return d


def n(p, col=None, val=None):
    try:
        d = pd.read_csv(p) if p.endswith("csv") else pd.read_parquet(p)
        return int(len(d) if col is None else (d[col] == val).sum())
    except Exception:
        return -1


S = json.load(open("results/step1c_score.json"))
cfg = W.bed_cfg()
print("STEP1c 2005-09 POLARIS bed -- DONE", pd.Timestamp.now("UTC").isoformat())
print(f"window {cfg['day0']}..{cfg['day1']} + {len(cfg['quiet_days'])} literature-quiet days; "
      f"hours from the ARRAY-ENERGY flag")
print(f"sites {cfg['sites']} (spread {cfg['moveout_spread_s_final']:.1f} s); withheld {cfg['withheld']}")
print(f"T = {j(f'{O}/T.json','T')} ({j(f'{O}/T.json','binds')}), "
      f"{j(f'{O}/T.json','n_quiet_hours')} quiet hours")
print(f"tremor hours {j(f'{O}/ac_meta.json','n_hours')}, pairs {j(f'{O}/ac_meta.json','n_pairs')}, "
      f"{j(f'{O}/ac_meta.json','median_s_per_hour')} s/hour")
print(f"seed clusters {n(f'{O}/cl_seeds.csv')}; survivors "
      f"{n(f'{O}/p0_survival.csv','survives',True)}/{n(f'{O}/p0_survival.csv')}; "
      f"deduped kept {n(f'{O}/dedupe.csv','kept',True)}")
print(f"located {j(f'{O}/lg.json','n_located')}; Gates 0-2 PASS {j(f'{O}/lg.json','n_gate_pass')}")
print(f"isolated hidden rediscovered {S['n_isolated_rediscovered']}/10 (ceiling 5/10); "
      f"all-30 {S['n_hidden_rediscovered_all30']}/30")
print(f"centroid {S['median_centroid_km']} km; recall {S['median_recall']}; "
      f"purity {S['median_purity']} (3-km truth {S['median_purity_cluster3km']}); "
      f"split {S['spurious_split_rate']}")
print(f"gate-PASS matching no Bostock family {S['gate_pass_no_bostock']}; "
      f"site-shift null per realisation {S['null_gated_per_realisation']}")
print(f"VERDICT: {S['verdict']} (accept >=3/10, split <=0.1, null <=1)")
D = {}
try:
    D = json.load(open("results/step1c_diag.json"))
except Exception:
    pass
print(f"diag D1 z_self {D.get('D1_median_z_self')}; D2 true-pair z {D.get('D2_median_z_true')} vs "
      f"null {D.get('D2_median_z_null')} (frac>nullmax {D.get('D2_frac_true_gt_null_max')})")
print("   COMPARE  Step 1a (46-137 km): 14.4 / 2.13 / 1.98 / 0.16")
print(f"diag D3 pair-invariant {D.get('D3_true_matches')}/{D.get('D3_n_pairs_tested')} vs null "
      f"{D.get('D3_null_matches')}")
print("files: results/freeze_step1c.json (spec), freeze_step1c_full.json, step1c_score.json, "
      f"step1c_diag.json, {O}/")
PY
T1=$(date +%s)
echo "elapsed_s: $((T1-T0))" >> results/STEP1C_DONE
echo "=== [$(date -u +%FT%TZ)] STEP1C QUEUE COMPLETE in $((T1-T0)) s" >&2
