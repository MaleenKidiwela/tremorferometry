#!/usr/bin/env bash
# STEP 1b -- the same 2012 window, the CLOSE stations that exist then.
#
#   download + frozen-chain cache -> quiet-hour T -> autocorrelation -> clusters -> handoff ->
#   growth -> survival + template-scramble null -> 2 iterations -> dedupe -> locate/gates ->
#   3 whole-chain site-shift null realisations -> scoring vs the 10 isolated hidden families ->
#   engine diagnostics
#
# NO engine constant is retuned: window, hop, thresholds, survival, dedupe and gate rules are
# byte-identical to freeze_step1.json.  Only the BED changes (sites, withheld, availability).
#
# ROBUSTNESS (added 2026-09-20 after two silent deaths).  The deaths were NOT OOM and NOT a code
# fault: an interactive `pkill -f w_step1b_queue.sh`, run while fixing an unrelated bug, matched
# the detached instances too, and the queue had no signal trap, so SIGTERM left no marker and no
# diagnosis.  Three guards now:
#   * SINGLE INSTANCE -- flock on locks/step1b.lock; a second launch exits 0 with a note instead
#     of racing (w_step1a_refresh.sh execs this queue at its end, so races are expected);
#   * SIGNAL TRAP -- TERM/INT/HUP and any non-zero exit write results/STEP1B_FAILED with the
#     signal name, the stage that was running and the last 120 log lines;
#   * RESUMABLE -- every stage is guarded by its own output file, so a relaunch continues instead
#     of re-fetching and re-computing.  Delete the output of a stage to force it to re-run.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
S=scripts
O=results/step1b
LOG=logs/step1b_queue.log
mkdir -p "$O" logs locks
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export W_BED=b
unset W_MAX_DAYS

exec 9>locks/step1b.lock
if ! flock -n 9; then
  echo "=== [$(date -u +%FT%TZ)] another Step-1b queue holds locks/step1b.lock -- exiting 0" >&2
  exit 0
fi

CURRENT_STAGE="startup"
finish() {
  local rc=$? sig="${1:-exit}"
  if [ ! -f results/STEP1B_DONE ] && [ ! -f results/STEP1B_FAILED ]; then
    { echo "step: $CURRENT_STAGE"; echo "signal: $sig"; echo "exit_code: $rc";
      echo "pid: $$"; echo "time: $(date -u +%FT%TZ)";
      echo "--- last 120 log lines ---"; tail -120 "$LOG" 2>/dev/null; } \
      > results/STEP1B_FAILED
    echo "!!! [$(date -u +%FT%TZ)] terminated during '$CURRENT_STAGE' ($sig, rc=$rc)" >&2
  fi
}
trap 'finish SIGTERM; exit 143' TERM
trap 'finish SIGINT;  exit 130' INT
trap 'finish SIGHUP;  exit 129' HUP
trap 'finish exit' EXIT

rm -f results/STEP1B_DONE results/STEP1B_FAILED

step() {
  local name="$1"; shift
  CURRENT_STAGE="$name"
  echo "=== [$(date -u +%FT%TZ)] STEP $name" >&2
  "$@"
  local c=$?
  if [ $c -ne 0 ]; then
    { echo "step: $name"; echo "exit_code: $c"; echo "cmd: $*"; echo "--- last 120 log lines ---";
      tail -120 "$LOG"; } > results/STEP1B_FAILED
    echo "!!! FAILED at $name ($c)" >&2
    exit $c
  fi
}

# step_if <output-file> <name> <cmd...> -- skip the stage when its output already exists
step_if() {
  local out="$1"; shift
  local name="$1"
  if [ -s "$out" ]; then
    echo "=== [$(date -u +%FT%TZ)] SKIP $name (have $out)" >&2
    return 0
  fi
  step "$@"
}
T0=$(date +%s)

step_if results/step1b_cache_manifest.csv fetch_bed \
        python $S/w_fetch_bed2012b.py --workers 4 --cache-workers 10

step_if results/freeze_step1b.json freeze \
        python $S/w_freeze.py --out results/freeze_step1b.json --tag step1b \
  --deviation "BED CHANGED from permanent-only to permanent+close-2012 because the permanent bed has no single-pair SNR margin (Step 1a: true-pair z 2.13 vs null 1.98, moveout spread 13.8 s). MERLIN TO AUDIT. C8 SOKB and CN PFB are metadata-open but deliver NO waveforms for this window at either IRIS or NRCan (results/step1b_not_delivered.csv), so the only NEW close site is CN LZB at 28.0 km epicentral (54/130 families within 25 km); PGCq gains five C8 PA0x members. Sites are matched by MEMBER station, never by the data-driven site id. Site rule: maximal deliverable subset with median S-moveout spread <= AC_WIN_S. NO engine constant retuned."

step_if $O/ac_quiet.parquet ac_quiet   python $S/w_autocorr.py --hours quiet  --out $O/ac_quiet --t-record 4.0
step_if $O/T.json           threshold  python $S/w_cluster.py  threshold --quiet-pairs $O/ac_quiet.parquet --out $O/T
step_if $O/ac.parquet       ac_tremor  python $S/w_autocorr.py --hours tremor --out $O/ac --t-record 5.0
step_if $O/cl_seeds.csv     cluster    python $S/w_cluster.py  cluster --pairs $O/ac.parquet --thr $O/T.json --out $O/cl

step_if $O/bank0.npz        handoff    python $S/w_grow.py handoff  --clusters $O/cl_clusters.parquet --seeds $O/cl_seeds.csv --out $O/bank0.npz
step_if $O/p0_thresholds.csv quiet0    python $S/w_grow.py quiet    --bank $O/bank0.npz --out $O/p0
step_if $O/p0_hits.parquet  grow0      python $S/w_grow.py grow     --bank $O/bank0.npz --out $O/p0 --thr $O/p0_thresholds.csv
step_if $O/bank0scr.npz     scramble   python $S/w_grow.py scramble --clusters $O/cl_clusters.parquet --seeds $O/cl_seeds.csv --out $O/bank0scr.npz
step_if $O/p0scr_hits.parquet grow0scr python $S/w_grow.py grow     --bank $O/bank0scr.npz --out $O/p0scr --thr $O/p0_thresholds.csv
step_if $O/p0_survival.csv  survive0   python $S/w_grow.py survive  --bank $O/bank0.npz --out $O/p0 --seeds $O/cl_seeds.csv \
                       --scramble-hits $O/p0scr_hits.parquet --scramble-bank $O/bank0scr.npz

for it in 1 2; do
  prev=$((it-1))
  step_if $O/bank$it.npz          iterate$it python $S/w_iterate.py --bank $O/bank$prev.npz --indep $O/p${prev}_indep.parquet \
                         --survivors $O/p0_survival.csv --out-bank $O/bank$it.npz --out-diag $O/iter${it}_diag.csv
  step_if $O/p${it}_thresholds.csv quiet$it  python $S/w_grow.py quiet --bank $O/bank$it.npz --out $O/p$it
  step_if $O/p${it}_hits.parquet   grow$it   python $S/w_grow.py grow  --bank $O/bank$it.npz --out $O/p$it --thr $O/p${it}_thresholds.csv
  step_if $O/p${it}_survival.csv   surv$it   python $S/w_grow.py survive --bank $O/bank$it.npz --out $O/p$it --seeds $O/cl_seeds.csv
done

step_if $O/dedupe.csv dedupe python $S/w_dedupe.py --bank $O/bank2.npz --indep $O/p2_indep.parquet --out $O/dedupe.csv
step_if $O/lg.json    locate python $S/w_locate_gates.py --bank $O/bank2.npz --indep $O/p2_indep.parquet \
                       --dedupe $O/dedupe.csv --survival $O/p2_survival.csv --out $O/lg
step_if $O/null_summary.csv nulls python $S/w_null_siteshift.py --out-root $O --n 3
step_if results/step1b_score.json score python $S/w_score_step1.py --clusters $O/lg_clusters.parquet --gates $O/lg_gates.csv \
                       --nulls $O/null_summary.csv --out results/step1b_score
step_if results/step1b_diag.json diag python $S/w_diag.py --out results/step1b_diag --pairs $O/ac.parquet

CURRENT_STAGE="marker"
python - <<'PY' > results/STEP1B_DONE
import json, os, sys
import pandas as pd
sys.path.insert(0, "scripts")
os.environ["W_BED"] = "b"
import w_common as W
O = "results/step1b"


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


S = json.load(open("results/step1b_score.json"))
mem, of, order, det = W.sites()
print("STEP1b close-2012 bed -- DONE", pd.Timestamp.now("UTC").isoformat())
print(f"sites: {order}  ({len(det)} stations); withheld {list(W.WITHHELD)}")
print(f"autocorrelation T = {j(f'{O}/T.json','T')} ({j(f'{O}/T.json','binds')}), "
      f"{j(f'{O}/T.json','n_quiet_hours')} quiet hours")
print(f"tremor hours {j(f'{O}/ac_meta.json','n_hours')}, pairs {j(f'{O}/ac_meta.json','n_pairs')}, "
      f"{j(f'{O}/ac_meta.json','median_s_per_hour')} s/hour")
print(f"seed clusters {n(f'{O}/cl_seeds.csv')}; pass-0 survivors "
      f"{n(f'{O}/p0_survival.csv','survives',True)} of {n(f'{O}/p0_survival.csv')}")
print(f"deduped kept {n(f'{O}/dedupe.csv','kept',True)}; located {j(f'{O}/lg.json','n_located')}; "
      f"Gates 0-2 PASS {j(f'{O}/lg.json','n_gate_pass')}")
print(f"isolated hidden rediscovered {S['n_isolated_rediscovered']}/10 (ceiling 5/10); "
      f"all-30 {S['n_hidden_rediscovered_all30']}/30 (ceiling 9/30)")
print(f"centroid {S['median_centroid_km']} km; recall {S['median_recall']}; purity {S['median_purity']} "
      f"(3-km truth {S['median_purity_cluster3km']}); split {S['spurious_split_rate']}")
print(f"gate-PASS matching no Bostock family: {S['gate_pass_no_bostock']}")
print(f"site-shift null gated per realisation: {S['null_gated_per_realisation']}")
print(f"VERDICT: {S['verdict']} (accept >=3/10, split <=0.1, null <=1)")
D = {}
try:
    D = json.load(open("results/step1b_diag.json"))
except Exception:
    pass
print(f"diag D1 z_self ceiling {D.get('D1_median_z_self')}")
print(f"diag D2 true-pair z {D.get('D2_median_z_true')} vs null {D.get('D2_median_z_null')} "
      f"(frac > null max {D.get('D2_frac_true_gt_null_max')})  <- COMPARE with Step 1a 2.13/1.98/0.16")
print(f"diag D3 pair-invariant {D.get('D3_true_matches')}/{D.get('D3_n_pairs_tested')} vs null "
      f"{D.get('D3_null_matches')}")
print("files: results/freeze_step1b.json, results/step1b_score.json, results/step1b_diag.json, "
      f"{O}/")
PY
T1=$(date +%s)
echo "elapsed_s: $((T1-T0))" >> results/STEP1B_DONE
echo "=== [$(date -u +%FT%TZ)] STEP1B QUEUE COMPLETE in $((T1-T0)) s" >&2
