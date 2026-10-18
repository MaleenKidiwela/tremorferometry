#!/usr/bin/env bash
# STEP 1 -- north positive control of the array-discovery engine.  ONE queue, no idle waiting.
#
#   quiet-hour T -> autocorrelation on the 2012 tremor hours -> clusters -> handoff -> growth ->
#   survival + template-scramble null -> 2 iterations -> dedupe -> locate/gates ->
#   3 whole-chain site-shift null realisations -> scoring + engine diagnostics
#
# Markers: results/STEP1_DONE (20-line summary) or results/STEP1_FAILED (step, code, last 120 log
# lines).  Nothing here touches traces/, cache20hz/ or results/markers/ in this study folder and
# nothing makes an FDSN request.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
S=scripts
O=results/step1
LOG=logs/step1_queue.log
mkdir -p "$O" logs
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
unset W_MAX_DAYS
rm -f results/STEP1_DONE results/STEP1_FAILED

step() {
  local name="$1"; shift
  echo "=== [$(date -u +%FT%TZ)] STEP $name" >&2
  "$@"
  local c=$?
  if [ $c -ne 0 ]; then
    { echo "step: $name"; echo "exit_code: $c"; echo "cmd: $*"; echo "--- last 120 log lines ---";
      tail -120 "$LOG"; } > results/STEP1_FAILED
    echo "!!! FAILED at $name ($c)" >&2
    exit $c
  fi
}

T0=$(date +%s)

# ---- 0. the freeze must already exist and its hashes must match the scripts on disk
step freeze_check python - <<'PY'
import json, os, sys
sys.path.insert(0, "scripts")
import w_common as W
F = json.load(open("results/freeze_step1.json"))
bad = [os.path.basename(p) for p in W.w_scripts()
       if F["sha256"]["w_scripts"].get(os.path.basename(p)) != W.sha256_file(p)]
print("[freeze] scripts changed since the freeze:", bad or "none")
assert not bad, f"re-freeze before running: {bad}"
PY

# ---- 1. quiet hours -> the FA-only autocorrelation threshold T (no 8.0 floor)
step ac_quiet   python $S/w_autocorr.py --hours quiet  --out $O/ac_quiet --t-record 4.0
step threshold  python $S/w_cluster.py  threshold --quiet-pairs $O/ac_quiet.parquet --out $O/T

# ---- 2. autocorrelation on the 2012 tremor hours -> clusters
step ac_tremor  python $S/w_autocorr.py --hours tremor --out $O/ac --t-record 5.0
step cluster    python $S/w_cluster.py  cluster --pairs $O/ac.parquet --thr $O/T.json --out $O/cl

# ---- 3. handoff -> growth -> scramble null -> survival
step handoff    python $S/w_grow.py handoff  --clusters $O/cl_clusters.parquet --seeds $O/cl_seeds.csv --out $O/bank0.npz
step quiet0     python $S/w_grow.py quiet    --bank $O/bank0.npz --out $O/p0
step grow0      python $S/w_grow.py grow     --bank $O/bank0.npz --out $O/p0 --thr $O/p0_thresholds.csv
step scramble   python $S/w_grow.py scramble --clusters $O/cl_clusters.parquet --seeds $O/cl_seeds.csv --out $O/bank0scr.npz
step grow0scr   python $S/w_grow.py grow     --bank $O/bank0scr.npz --out $O/p0scr --thr $O/p0_thresholds.csv
step survive0   python $S/w_grow.py survive  --bank $O/bank0.npz --out $O/p0 --seeds $O/cl_seeds.csv \
                       --scramble-hits $O/p0scr_hits.parquet --scramble-bank $O/bank0scr.npz

# ---- 4. two iterations (12-s templates from S - 2 s)
for it in 1 2; do
  prev=$((it-1))
  step iterate$it python $S/w_iterate.py --bank $O/bank$prev.npz --indep $O/p${prev}_indep.parquet \
                         --survivors $O/p0_survival.csv --out-bank $O/bank$it.npz --out-diag $O/iter${it}_diag.csv
  step quiet$it   python $S/w_grow.py quiet --bank $O/bank$it.npz --out $O/p$it
  step grow$it    python $S/w_grow.py grow  --bank $O/bank$it.npz --out $O/p$it --thr $O/p${it}_thresholds.csv
  step surv$it    python $S/w_grow.py survive --bank $O/bank$it.npz --out $O/p$it --seeds $O/cl_seeds.csv
done

# ---- 5. dedupe -> locate -> Gates 0-2
step dedupe     python $S/w_dedupe.py --bank $O/bank2.npz --indep $O/p2_indep.parquet --out $O/dedupe.csv
step locate     python $S/w_locate_gates.py --bank $O/bank2.npz --indep $O/p2_indep.parquet \
                       --dedupe $O/dedupe.csv --survival $O/p2_survival.csv --out $O/lg

# ---- 6. whole-chain site-shift null, 3 realisations
step nulls      python $S/w_null_siteshift.py --out-root $O --n 3

# ---- 7. scoring against the 10 isolated hidden families + engine diagnostics
step score      python $S/w_score_step1.py --clusters $O/lg_clusters.parquet --gates $O/lg_gates.csv \
                       --nulls $O/null_summary.csv --out results/step1_score
step diag       python $S/w_diag.py --out results/step1_diag --pairs $O/ac.parquet

# ---- 8. marker
python - <<'PY' > results/STEP1_DONE
import json, os, sys
import pandas as pd
sys.path.insert(0, "scripts")
O = "results/step1"


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


S = json.load(open("results/step1_score.json"))
print("STEP1 north positive control -- DONE", pd.Timestamp.now("UTC").isoformat())
print(f"bed: 12 stations / 6 sites / 52 days; withheld VGZ SNB B926")
print(f"autocorrelation T = {j(f'{O}/T.json','T')} ({j(f'{O}/T.json','binds')}), "
      f"{j(f'{O}/T.json','n_quiet_hours')} quiet hours")
print(f"tremor hours: {j(f'{O}/ac_meta.json','n_hours')}, pairs >= 5.0: "
      f"{j(f'{O}/ac_meta.json','n_pairs')}, median {j(f'{O}/ac_meta.json','median_s_per_hour')} s/hour")
print(f"seed clusters: {n(f'{O}/cl_seeds.csv')}")
print(f"pass-0 survivors: {n(f'{O}/p0_survival.csv','survives',True)} of {n(f'{O}/p0_survival.csv')}")
print(f"iter-2 templates: {n(f'{O}/p2_survival.csv')}; deduped kept: {n(f'{O}/dedupe.csv','kept',True)}")
print(f"located: {j(f'{O}/lg.json','n_located')}; Gates 0-2 PASS: {j(f'{O}/lg.json','n_gate_pass')}")
print(f"isolated hidden rediscovered: {S['n_isolated_rediscovered']}/10  "
      f"(ceiling: frozen template 5/10)")
print(f"all-30 hidden rediscovered: {S['n_hidden_rediscovered_all30']}/30 (ceiling 9/30)")
print(f"median centroid distance: {S['median_centroid_km']} km; median recall {S['median_recall']}; "
      f"purity {S['median_purity']} (3-km cluster truth {S['median_purity_cluster3km']})")
print(f"spurious-split rate: {S['spurious_split_rate']} (bar <= 0.10)")
print(f"gate-PASS clusters matching no Bostock family: {S['gate_pass_no_bostock']}")
print(f"site-shift null gated families per realisation: {S['null_gated_per_realisation']}")
print(f"VERDICT: {S['verdict']}  (accept: >= 3/10 isolated, split <= 0.1, null <= 1)")
D = {}
try:
    D = json.load(open("results/step1_diag.json"))
except Exception:
    pass
print(f"diag D1 median z_self (perfect-repeat ceiling): {D.get('D1_median_z_self')}")
print(f"diag D2 true-pair z {D.get('D2_median_z_true')} vs jitter null {D.get('D2_median_z_null')} "
      f"(frac > null max {D.get('D2_frac_true_gt_null_max')})")
print(f"diag D3 pair-invariant matches {D.get('D3_true_matches')} of {D.get('D3_n_pairs_tested')} "
      f"vs null {D.get('D3_null_matches')}")
print(f"files: results/freeze_step1.json, results/step1_score.json, results/step1_diag.json, {O}/")
PY
T1=$(date +%s)
echo "elapsed_s: $((T1-T0))" >> results/STEP1_DONE
echo "=== [$(date -u +%FT%TZ)] STEP1 QUEUE COMPLETE in $((T1-T0)) s" >&2
