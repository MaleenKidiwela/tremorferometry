#!/usr/bin/env bash
# ARM B -- the north control with LFEDetect(`jms`) on the 2005 POLARIS bed.
#
# Merlin's Step-1 audit (Q3): `cascadia`/`cms`/`jcms`/`jcs` are trained on Bostock 2015 stacks, so
# only `jms` (Guerrero + San Andreas + Nankai, no Cascadia) makes this a real control.  This is the
# only path that survives a STOP of Arm A, and WA must later run the SAME weights.
#
#   stage 1  model card -> results/armB/jms_model_card.json + the freeze
#   stage 2  LFEDetect(jms) picks over the 2005 bed (14 detection stations, ETS + quiet days)
#   stage 3  PyOcto association, Crosson 1-D, Slab2 28-45 km depth prior, at 0.30 and 0.15
#   stage 4  DBSCAN 5 km / whole window, clusters >= 20 events
#   stage 5  frozen handoff -> quiet T_seed -> ONE growth pass -> scramble -> survive -> dedupe
#            -> locate + Gates 0-2
#   stage 6  scoring vs the 10 isolated hidden families (absolute bars; relative if the 2005
#            oracle reference exists)
#   stage 7  ONE site-shift null realisation of stages 2-5 (runs BEFORE the score so the null
#            column is filled)
#
# The engine is a SNAPSHOT of scripts/ taken at queue start (scripts_armB/), because scripts/w_*.py
# are being edited concurrently for Step 1c'.  The GPU is SHARED with the Step-1c' queue.
# Markers: results/ARMB_DONE / results/ARMB_FAILED.  Log: logs/armB_queue.log.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
ROOT=/home/jovyan/tremorferometry/lfe_array_discovery
S=$ROOT/scripts_armB
O=results/armB
LOG=logs/armB_queue.log
LOCK=locks/armB_queue.lock
mkdir -p "$O" logs locks results

# ------------------------------------------------------------------ single instance (flock) ----
exec 9>"$LOCK"
if ! flock -n 9; then
  echo "=== [$(date -u +%FT%TZ)] another Arm-B queue holds $LOCK -- exiting" >&2
  exit 0
fi

export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
export W_BED=c
# The Step-1c' repair flags must NEVER be active here: Arm B hands off to the FROZEN 1c chain.
unset W_HEALTH_VETO W_WIN_TRANSIENT W_DROP_ZERO_MEDIAN W_SEED_RULE W_SEED_BUDGET \
      W_HANDOFF_MIN_SEG W_NEAR_SEED_MODE W_QUIET_HOURS_CSV W_MAX_DAYS 2>/dev/null || true

STAGE="init"
T0=$(date +%s)

fail() {
  local c="$1"
  { echo "stage: $STAGE"; echo "exit_code: $c"; echo "utc: $(date -u +%FT%TZ)";
    echo "elapsed_s: $(( $(date +%s) - T0 ))"; echo "--- last 150 log lines ---";
    tail -150 "$LOG" 2>/dev/null; } > results/ARMB_FAILED
  echo "!!! ARM B FAILED at $STAGE ($c)" >&2
}

on_exit() {
  local c=$?
  if [ ! -f results/ARMB_DONE ]; then fail "$c"; fi
}
trap on_exit EXIT
trap 'echo "=== signal received at stage $STAGE" >&2; exit 143' TERM INT HUP

step() {                       # step <name> <cmd...>
  STAGE="$1"; shift
  echo "=== [$(date -u +%FT%TZ)] STAGE $STAGE" >&2
  "$@" || { local c=$?; fail "$c"; exit "$c"; }
}

guard() {                      # guard <marker-file> -- true when the stage may be skipped
  [ -e "$1" ] && { echo "--- [$(date -u +%FT%TZ)] skip (exists: $1)" >&2; return 0; } || return 1
}

rm -f results/ARMB_DONE results/ARMB_FAILED

# ------------------------------------------------------------------ stage 0: engine snapshot ---
if [ ! -d "$S" ]; then
  mkdir -p "$S"
  cp $ROOT/scripts/*.py "$S"/
  echo "=== [$(date -u +%FT%TZ)] engine snapshot -> $S ($(ls $S/*.py | wc -l) files)" >&2
fi
PY=$(command -v python)

# ------------------------------------------------------------------ stage 1: card + freeze -----
step freeze $PY $S/b_armB_freeze.py

# ------------------------------------------------------------------ stage 2: picks -------------
step picks $PY $S/b_armB_pick.py --out $O/picks --days all

# ------------------------------------------------------------------ stage 3: association -------
guard $O/events_t30.json || step assoc30 $PY $S/b_armB_assoc.py --picks "$O/picks_*.parquet" \
  --thr 0.30 --out $O/events_t30
guard $O/events_t15.json || step assoc15 $PY $S/b_armB_assoc.py --picks "$O/picks_*.parquet" \
  --thr 0.15 --out $O/events_t15
# DIAGNOSTIC only (never a gate): the 4-station variant of the default level
guard $O/events_t30_n4.json || step assoc30n4 $PY $S/b_armB_assoc.py --picks "$O/picks_*.parquet" \
  --thr 0.30 --n-picks 4 --min-stations 4 --out $O/events_t30_n4

# ------------------------------------------------------------------ stage 4: clustering --------
guard $O/cl30_meta.json || step cluster30 $PY $S/b_armB_cluster.py --events $O/events_t30.csv \
  --out $O/cl30
guard $O/cl15_meta.json || step cluster15 $PY $S/b_armB_cluster.py --events $O/events_t15.csv \
  --out $O/cl15

# pre-registered carry rule: the DEFAULT level, unless it yields zero clusters
CARRY=cl30
if [ "$($PY -c "import json;print(json.load(open('$O/cl30_meta.json'))['n_clusters'])")" = "0" ]; then
  CARRY=cl15
fi
if [ "$CARRY" = "cl30" ]; then NTHR=0.30; else NTHR=0.15; fi
echo "carry_level: $CARRY (threshold $NTHR)" | tee $O/carry_level.txt >&2

# ------------------------------------------------------------------ stage 5: frozen growth -----
guard $O/bank0.npz  || step handoff  $PY $S/w_grow.py handoff --clusters $O/${CARRY}_clusters.parquet \
  --seeds $O/${CARRY}_seeds.csv --out $O/bank0.npz
guard $O/p0_thresholds.csv || step quiet0 $PY $S/w_grow.py quiet --bank $O/bank0.npz --out $O/p0
guard $O/p0_hits.parquet   || step grow0  $PY $S/w_grow.py grow  --bank $O/bank0.npz --out $O/p0 \
  --thr $O/p0_thresholds.csv
guard $O/bank0scr.npz || step scramble $PY $S/w_grow.py scramble \
  --clusters $O/${CARRY}_clusters.parquet --seeds $O/${CARRY}_seeds.csv --out $O/bank0scr.npz
guard $O/p0scr_hits.parquet || step grow0scr $PY $S/w_grow.py grow --bank $O/bank0scr.npz \
  --out $O/p0scr --thr $O/p0_thresholds.csv
guard $O/p0_survival.csv || step survive0 $PY $S/w_grow.py survive --bank $O/bank0.npz --out $O/p0 \
  --seeds $O/${CARRY}_seeds.csv --scramble-hits $O/p0scr_hits.parquet \
  --scramble-bank $O/bank0scr.npz
guard $O/dedupe.csv || step dedupe $PY $S/w_dedupe.py --bank $O/bank0.npz \
  --indep $O/p0_indep.parquet --out $O/dedupe.csv
guard $O/lg.json || step locate $PY $S/w_locate_gates.py --bank $O/bank0.npz \
  --indep $O/p0_indep.parquet --dedupe $O/dedupe.csv --survival $O/p0_survival.csv --out $O/lg

# ------------------------------------------------------------------ stage 7: site-shift null ---
mkdir -p $O/null1
if ! guard $O/null1/lg.json; then
  STAGE="null1"
  SHIFT=$($PY -c "
import json,os,sys
sys.path.insert(0,'$S'); os.environ['W_BED']='c'
import w_common as W
print(json.dumps(W.site_shift(1)))")
  echo "null1 shift: $SHIFT" >&2
  step null1_assoc   $PY $S/b_armB_assoc.py --picks "$O/picks_*.parquet" --thr $NTHR \
    --shift "$SHIFT" --out $O/null1/events
  step null1_cluster $PY $S/b_armB_cluster.py --events $O/null1/events.csv --out $O/null1/cl
  step null1_handoff $PY $S/w_grow.py handoff --clusters $O/null1/cl_clusters.parquet \
    --seeds $O/null1/cl_seeds.csv --out $O/null1/bank0.npz --shift "$SHIFT"
  step null1_quiet   $PY $S/w_grow.py quiet --bank $O/null1/bank0.npz --out $O/null1/p0 \
    --shift "$SHIFT"
  step null1_grow    $PY $S/w_grow.py grow --bank $O/null1/bank0.npz --out $O/null1/p0 \
    --thr $O/null1/p0_thresholds.csv --shift "$SHIFT"
  step null1_surv    $PY $S/w_grow.py survive --bank $O/null1/bank0.npz --out $O/null1/p0 \
    --seeds $O/null1/cl_seeds.csv
  step null1_dedupe  $PY $S/w_dedupe.py --bank $O/null1/bank0.npz \
    --indep $O/null1/p0_indep.parquet --out $O/null1/dedupe.csv
  step null1_locate  $PY $S/w_locate_gates.py --bank $O/null1/bank0.npz \
    --indep $O/null1/p0_indep.parquet --dedupe $O/null1/dedupe.csv \
    --survival $O/null1/p0_survival.csv --out $O/null1/lg --shift "$SHIFT"
fi
step null_summary $PY - <<PY
import json, os
import pandas as pd
rows = []
for k in (1,):
    p = "$O/null%d/lg.json" % k
    if os.path.exists(p):
        j = json.load(open(p)); j["realisation"] = k
        rows.append(j)
pd.DataFrame(rows).to_csv("$O/null_summary.csv", index=False)
print("[null] ", rows)
PY

# ------------------------------------------------------------------ stage 6: scoring -----------
step score $PY $S/b_armB_score.py --clusters $O/lg_clusters.parquet --gates $O/lg_gates.csv \
  --nulls $O/null_summary.csv --out results/armB_score

# ------------------------------------------------------------------ done ----------------------
STAGE="report"
$PY - <<PY > results/ARMB_DONE
import json, os, glob
import pandas as pd

O = "$O"


def j(p, k=None, d=None):
    try:
        v = json.load(open(p))
        return v if k is None else v.get(k, d)
    except Exception:
        return d


def n(p, col=None, val=None):
    try:
        d = pd.read_csv(p) if p.endswith("csv") else pd.read_parquet(p)
        return int(len(d) if col is None else (d[col] == val).sum())
    except Exception:
        return -1


print("ARM B -- LFEDetect(jms) north control on the 2005 POLARIS bed -- DONE",
      pd.Timestamp.now("UTC").isoformat())
P = sorted(glob.glob(f"{O}/picks_*.parquet"))
tot = 0
per = []
for p in P:
    try:
        d = pd.read_parquet(p)
        tot += len(d)
        if len(d):
            per.append((d.sta.iloc[0], len(d), d.doy.nunique(), len(d) / max(d.doy.nunique(), 1)))
    except Exception:
        pass
print(f"picks: {tot} over {sum(x[2] for x in per)} station-days, {len(per)} stations")
for s, k, nd, r in sorted(per, key=lambda x: -x[3]):
    print(f"   {s:6s} {k:8d} picks / {nd:3d} station-days = {r:7.0f} per station-day")
for lvl in ("t30", "t15", "t30_n4"):
    print(f"assoc {lvl}: {j(f'{O}/events_{lvl}.json', 'n_events')} events from "
          f"{j(f'{O}/events_{lvl}.json', 'n_picks')} picks; median stations "
          f"{j(f'{O}/events_{lvl}.json', 'median_n_stations')}; depth median "
          f"{j(f'{O}/events_{lvl}.json', 'depth_median')}")
for lvl in ("cl30", "cl15"):
    print(f"clusters {lvl}: {j(f'{O}/{lvl}_meta.json', 'n_clusters')} "
          f"(median size {j(f'{O}/{lvl}_meta.json', 'median_cluster_size')})")
print("carried level:", open(f"{O}/carry_level.txt").read().strip())
print(f"templates handed off: {n(f'{O}/bank0_diag.csv')} station-rows; survivors "
      f"{n(f'{O}/p0_survival.csv', 'survives', True)}/{n(f'{O}/p0_survival.csv')}; "
      f"deduped kept {n(f'{O}/dedupe.csv', 'kept', True)}")
print(f"located {j(f'{O}/lg.json', 'n_located')}; Gates 0-2 PASS {j(f'{O}/lg.json', 'n_gate_pass')}")
print(f"null1: clusters {j(f'{O}/null1/cl_meta.json', 'n_clusters')}, survivors "
      f"{n(f'{O}/null1/p0_survival.csv', 'survives', True)}, located "
      f"{j(f'{O}/null1/lg.json', 'n_located')}, gated {j(f'{O}/null1/lg.json', 'n_gate_pass')}")
A = j("results/armB_score_armB.json") or {}
print(f"ISOLATED rediscovered (absolute bars): {A.get('n_isolated_rediscovered_abs')}/10 "
      f"(relative {A.get('n_isolated_rediscovered_relative')}/10)")
print("relative reference:", A.get("relative_reference"))
print(f"median centroid {A.get('median_centroid_km')} km; median recall {A.get('median_recall')}; "
      f"median purity(3 km) {A.get('median_purity_cluster3km')}")
print(f"VERDICT: {A.get('verdict')}  (bars {A.get('bars')})")
print("STOP rule:", A.get("stop_rule"))
print("files: results/freeze_armB.json, results/armB_score{,_armB,_isolated_detail}.json/csv, "
      f"{O}/")
PY
echo "elapsed_s: $(( $(date +%s) - T0 ))" >> results/ARMB_DONE
echo "=== [$(date -u +%FT%TZ)] ARM B QUEUE COMPLETE in $(( $(date +%s) - T0 )) s" >&2
