#!/usr/bin/env bash
# STEP 5 -- N-OLYMPIC 2010 DISCOVERY (Washington).  Merlin prescription Step 5 + Step 5b.
#
#   hash check -> GPU GATE -> quiet-hour T -> autocorrelation on the candidate hours -> clusters
#   -> handoff -> growth -> template-scramble null -> survival -> 2 (+1) iterations -> dedupe
#   -> locate + Gates 0-2 with the withheld sites -> whole-chain site-shift null x3
#   -> Step-5b injections (build, same chain on the injected cache, score) -> candidate table.
#
# NO engine constant is retuned: every `w_` stage runs byte-identical to Step 1c through
# `x_stage.py`, which binds the Washington bed in process.  The spec was frozen BEFORE launch in
# results/freeze_wa_v1.json and stage 0 refuses to run if any script has changed since.
#
# GPU GATE (stage 1).  The queue does not touch the GPU until BOTH
#     (a) results/WA_GO exists -- the orchestrator creates it after the north control passes, AND
#     (b) no `w_step1` process is running (Step 1a/1b/1c own the GPU until they are done).
# It polls; it NEVER kills anything and it never runs `pkill` on a shared name.
#
# SMOKE MODE.  X_SMOKE=1 skips the GPU gate, but ONLY when results/SMOKE_OK exists AND
# results/WA_GO does not.  It runs one tremor hour, one cluster, two growth days, a fake-proposal
# location/gate pass and a one-station-day injection into results/step5_smoke, and writes
# results/STEP5_SMOKE_DONE|FAILED instead of the real markers.
#
# Markers: results/STEP5_DONE / results/STEP5_FAILED.  Log: logs/step5_queue.log.
# Single instance via flock on locks/step5.lock; TERM/INT/HUP/EXIT always write a FAILED marker
# naming the signal, the stage and the last 120 log lines; every stage is skipped when its own
# output already exists, so a relaunch resumes where it stopped.
set -u
cd /home/jovyan/tremorferometry/lfe_array_discovery
S=scripts
LOG=logs/step5_queue.log
mkdir -p logs locks results

SMOKE=${X_SMOKE:-0}
if [ "$SMOKE" = "1" ]; then
  O=results/step5_smoke
  MDONE=results/STEP5_SMOKE_DONE
  MFAIL=results/STEP5_SMOKE_FAILED
  LOCK=locks/step5_smoke.lock
else
  O=results/step5
  MDONE=results/STEP5_DONE
  MFAIL=results/STEP5_FAILED
  LOCK=locks/step5.lock
fi
mkdir -p "$O"

# ---------------------------------------------------------------- single instance
exec 9>"$LOCK"
if ! flock -n 9; then
  echo "=== [$(date -u +%FT%TZ)] another Step-5 queue holds $LOCK -- exiting 0" >&2
  exit 0
fi

export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8
unset W_MAX_DAYS
unset X_CACHE

STAGE="startup"
FINISHED=0

fail_marker() {
  { echo "step: $STAGE";
    echo "reason: $1";
    echo "time: $(date -u +%FT%TZ)";
    echo "pid: $$";
    echo "--- last 120 log lines ---";
    tail -120 "$LOG" 2>/dev/null; } > "$MFAIL"
}

on_signal() {
  STAGE="$STAGE (killed by SIG$1)"
  fail_marker "queue received SIG$1 while running stage '$STAGE'"
  echo "!!! [$(date -u +%FT%TZ)] SIG$1 during $STAGE -- wrote $MFAIL" >&2
  exit 143
}
trap 'on_signal TERM' TERM
trap 'on_signal INT'  INT
trap 'on_signal HUP'  HUP
trap '[ "$FINISHED" = "1" ] || [ -f "$MFAIL" ] || fail_marker "queue exited without finishing stage \"$STAGE\""' EXIT

step() {
  local name="$1"; shift
  STAGE="$name"
  echo "=== [$(date -u +%FT%TZ)] STAGE $name" >&2
  "$@"
  local c=$?
  if [ $c -ne 0 ]; then
    fail_marker "stage exited $c: $*"
    echo "!!! FAILED at $name ($c)" >&2
    exit $c
  fi
}

# resume guard: run the stage only if its output file is missing
guard() {
  local out="$1"; local name="$2"; shift 2
  if [ -e "$out" ]; then
    echo "=== [$(date -u +%FT%TZ)] STAGE $name -- $out exists, skip" >&2
    return 0
  fi
  step "$name" "$@"
}

rm -f "$MDONE" "$MFAIL"
T0=$(date +%s)
XS="python $S/x_stage.py"

# ================================================================ 0. hash check vs the freeze ===
step hash_check python - <<'PY'
import json, os, sys
sys.path.insert(0, "scripts")
import x_common as X
F = json.load(open(X.FREEZE))
bad = []
for group, paths in (("x_scripts", X.x_scripts()), ("w_scripts", X.w_scripts())):
    for p in paths:
        b = os.path.basename(p)
        want = F["sha256"][group].get(b)
        if want is None:
            bad.append(f"{b} (not in the freeze)")
        elif want != X.sha256_file(p):
            bad.append(f"{b} (changed)")
for p, want in F["sha256"]["inputs"].items():
    pass
print("[hash] scripts differing from results/freeze_wa_v1.json:", bad or "none")
assert not bad, f"re-freeze before running: {bad}"
print("[hash] sub-region", F["sub_region"]["name"], F["window"]["day0"], "..",
      F["window"]["day1"], "|", F["sites"]["n_sites"], "sites /",
      F["sites"]["n_detection_stations"], "stations | withheld",
      F["withheld"]["stations"], "| spread",
      round(F["moveout"]["spread_s_at_centre"], 2), "s | HF veto",
      F["gates"]["gate0"]["hf_station"], round(F["gates"]["gate0"]["hf_veto"], 3))
PY

# ================================================================ 1. GPU GATE ===================
STAGE="gpu_gate"
if [ "$SMOKE" = "1" ]; then
  if [ -f results/WA_GO ] || [ ! -f results/SMOKE_OK ]; then
    fail_marker "smoke mode refused: it needs results/SMOKE_OK present AND results/WA_GO absent"
    echo "!!! smoke mode refused (SMOKE_OK missing or WA_GO present)" >&2
    exit 2
  fi
  echo "=== [$(date -u +%FT%TZ)] SMOKE MODE: GPU gate skipped (SMOKE_OK present, WA_GO absent)" >&2
else
  echo "=== [$(date -u +%FT%TZ)] STAGE gpu_gate: waiting for results/WA_GO and for no w_step1 process" >&2
  n=0
  while true; do
    if [ -f results/WA_GO ] && ! pgrep -f "[w]_step1" >/dev/null 2>&1; then
      echo "=== [$(date -u +%FT%TZ)] gate open: WA_GO present, no w_step1 process" >&2
      break
    fi
    n=$((n+1))
    if [ $((n % 30)) -eq 1 ]; then
      echo "    [$(date -u +%FT%TZ)] waiting: WA_GO=$([ -f results/WA_GO ] && echo yes || echo no) w_step1=$(pgrep -cf '[w]_step1' 2>/dev/null || echo 0)" >&2
    fi
    sleep 60
  done
fi

# ================================================================ 2. quiet-hour T ===============
if [ "$SMOKE" = "1" ]; then QLIM="--limit 4"; else QLIM=""; fi
guard $O/ac_quiet.parquet ac_quiet $XS w_autocorr --hours quiet --out $O/ac_quiet --t-record 4.0 $QLIM
guard $O/T.json threshold $XS w_cluster threshold --quiet-pairs $O/ac_quiet.parquet --out $O/T
step threshold_check python - "$O" <<'PY'
import json, sys
import pandas as pd
O = sys.argv[1]
T = json.load(open(f"{O}/T.json"))
G = pd.read_csv(f"{O}/T_grid.csv")
ok = G[(G["T"] <= 20.0) & (G.n_components_ge5 == 0)]
print(f"[T] T = {T['T']} ({T['binds']}) on {T['n_quiet_hours']} quiet hours; "
      f"grid rows <= 20 with zero >=5-member components: {len(ok)}")
if not len(ok):
    raise SystemExit("SITE-COHERENT NOISE: the quiet hours contain >= 5-member components at "
                     "every T <= 20 -- Merlin Step-3 failure signature")
PY

# ================================================================ 3. autocorrelation -> clusters
if [ "$SMOKE" = "1" ]; then ALIM="--limit 1"; else ALIM=""; fi
guard $O/ac.parquet   ac_tremor $XS w_autocorr --hours tremor --out $O/ac --t-record 5.0 $ALIM
guard $O/cl_seeds.csv cluster   $XS w_cluster cluster --pairs $O/ac.parquet --thr $O/T.json --out $O/cl

# ================================================================ 4. handoff -> growth -> survival
if [ "$SMOKE" = "1" ]; then export W_MAX_DAYS=2; fi
guard $O/bank0.npz        handoff   $XS w_grow handoff --clusters $O/cl_clusters.parquet --seeds $O/cl_seeds.csv --out $O/bank0.npz
guard $O/p0_thresholds.csv quiet0   $XS w_grow quiet   --bank $O/bank0.npz --out $O/p0
guard $O/p0_hits.parquet  grow0     $XS w_grow grow    --bank $O/bank0.npz --out $O/p0 --thr $O/p0_thresholds.csv
guard $O/bank0scr.npz     scramble  $XS w_grow scramble --clusters $O/cl_clusters.parquet --seeds $O/cl_seeds.csv --out $O/bank0scr.npz
guard $O/p0scr_hits.parquet grow0scr $XS w_grow grow   --bank $O/bank0scr.npz --out $O/p0scr --thr $O/p0_thresholds.csv
guard $O/p0_survival.csv  survive0  $XS w_grow survive --bank $O/bank0.npz --out $O/p0 --seeds $O/cl_seeds.csv \
      --scramble-hits $O/p0scr_hits.parquet --scramble-bank $O/bank0scr.npz

# ================================================================ 5. iterations -> dedupe =======
for it in 1 2; do
  prev=$((it-1))
  guard $O/bank$it.npz          iterate$it $XS w_iterate --bank $O/bank$prev.npz --indep $O/p${prev}_indep.parquet \
        --survivors $O/p0_survival.csv --out-bank $O/bank$it.npz --out-diag $O/iter${it}_diag.csv
  guard $O/p${it}_thresholds.csv quiet$it  $XS w_grow quiet   --bank $O/bank$it.npz --out $O/p$it
  guard $O/p${it}_hits.parquet   grow$it   $XS w_grow grow    --bank $O/bank$it.npz --out $O/p$it --thr $O/p${it}_thresholds.csv
  guard $O/p${it}_survival.csv   surv$it   $XS w_grow survive --bank $O/bank$it.npz --out $O/p$it --seeds $O/cl_seeds.csv
done

# third iteration ONLY if the survivors' S-SNR is still rising (the frozen rule: median
# `s_snr_db` of the `__global__` rows rose by >= 0.5 dB from iteration 1 to iteration 2)
step iter_snr_rule python - "$O" <<'PY'
import json, os, sys
import numpy as np
import pandas as pd
O = sys.argv[1]


def med(p):
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return float("nan")
    try:
        d = pd.read_csv(p)
    except Exception:
        return float("nan")
    if not len(d) or "s_snr_db" not in d:
        return float("nan")
    g = d[d.sta == "__global__"] if "sta" in d else d
    v = (g if len(g) else d).s_snr_db.astype(float)
    return float(np.nanmedian(v))


m1, m2 = med(f"{O}/iter1_diag.csv"), med(f"{O}/iter2_diag.csv")
rise = (m2 - m1) if (np.isfinite(m1) and np.isfinite(m2)) else float("nan")
run3 = bool(np.isfinite(rise) and rise >= 0.5)
json.dump(dict(iter1_median_s_snr_db=m1, iter2_median_s_snr_db=m2, rise_db=rise,
               rule="a third iteration runs only if the rise is >= 0.5 dB", run_third=run3),
          open(f"{O}/iter_snr.json", "w"), indent=1)
print(f"[iter] median S-SNR {m1} -> {m2} (rise {rise}); third iteration: {run3}")
PY
if python -c "import json,sys; sys.exit(0 if json.load(open('$O/iter_snr.json'))['run_third'] else 1)"; then
  guard $O/bank3.npz           iterate3 $XS w_iterate --bank $O/bank2.npz --indep $O/p2_indep.parquet \
        --survivors $O/p0_survival.csv --out-bank $O/bank3.npz --out-diag $O/iter3_diag.csv
  guard $O/p3_thresholds.csv   quiet3   $XS w_grow quiet   --bank $O/bank3.npz --out $O/p3
  guard $O/p3_hits.parquet     grow3    $XS w_grow grow    --bank $O/bank3.npz --out $O/p3 --thr $O/p3_thresholds.csv
  guard $O/p3_survival.csv     surv3    $XS w_grow survive --bank $O/bank3.npz --out $O/p3 --seeds $O/cl_seeds.csv
  LAST=3
else
  LAST=2
fi
echo "$LAST" > $O/last_iteration.txt
guard $O/dedupe.csv dedupe $XS w_dedupe --bank $O/bank$LAST.npz --indep $O/p${LAST}_indep.parquet --out $O/dedupe.csv

# ================================================================ 6. locate + Gates 0-2 ========
guard $O/lg_gates.csv locate_gates python $S/x_locate_gates.py --bank $O/bank$LAST.npz \
      --indep $O/p${LAST}_indep.parquet --dedupe $O/dedupe.csv --survival $O/p${LAST}_survival.csv \
      --out $O/lg

# ================================================================ 7. whole-chain site-shift null
if [ "$SMOKE" = "1" ]; then
  # plumbing check of the location + gate wiring, which the smoke bank is too small to reach
  guard $O/fake.csv fake_proposal python $S/x_smoke.py fake-proposal --n 40 --out $O/fake
  echo "=== [$(date -u +%FT%TZ)] SMOKE: site-shift null skipped" >&2
else
  guard $O/null_summary.csv nulls python $S/x_null_siteshift.py --out-root $O --n 3
fi

# ================================================================ 8. Step-5b injections ========
if [ "$SMOKE" = "1" ]; then
  guard $O/inject/planted.csv inject_build python $S/x_inject.py build --out $O/inject \
        --n-families 1 --n-events 8 --n-days 1 --amps 1.0 --only-stations B003
  echo "=== [$(date -u +%FT%TZ)] SMOKE: injected chain not run (builder only)" >&2
else
  guard $O/inject/planted.csv inject_build python $S/x_inject.py build --out $O/inject
  export X_CACHE="$PWD/$O/inject/cache20hz"
  I=$O/inject
  guard $I/ac.parquet   inj_ac      $XS w_autocorr --hours file --hours-file $I/hours.csv --out $I/ac --t-record 5.0
  guard $I/cl_seeds.csv inj_cluster $XS w_cluster cluster --pairs $I/ac.parquet --thr $O/T.json --out $I/cl
  guard $I/bank0.npz    inj_handoff $XS w_grow handoff --clusters $I/cl_clusters.parquet --seeds $I/cl_seeds.csv --out $I/bank0.npz
  guard $I/p0_thresholds.csv inj_quiet0 $XS w_grow quiet --bank $I/bank0.npz --out $I/p0
  guard $I/p0_hits.parquet   inj_grow0  $XS w_grow grow  --bank $I/bank0.npz --out $I/p0 --thr $I/p0_thresholds.csv
  guard $I/bank0scr.npz      inj_scr    $XS w_grow scramble --clusters $I/cl_clusters.parquet --seeds $I/cl_seeds.csv --out $I/bank0scr.npz
  guard $I/p0scr_hits.parquet inj_grow0scr $XS w_grow grow --bank $I/bank0scr.npz --out $I/p0scr --thr $I/p0_thresholds.csv
  guard $I/p0_survival.csv   inj_surv0  $XS w_grow survive --bank $I/bank0.npz --out $I/p0 --seeds $I/cl_seeds.csv \
        --scramble-hits $I/p0scr_hits.parquet --scramble-bank $I/bank0scr.npz
  for it in 1 2; do
    prev=$((it-1))
    guard $I/bank$it.npz          inj_iterate$it $XS w_iterate --bank $I/bank$prev.npz --indep $I/p${prev}_indep.parquet \
          --survivors $I/p0_survival.csv --out-bank $I/bank$it.npz --out-diag $I/iter${it}_diag.csv
    guard $I/p${it}_thresholds.csv inj_quiet$it  $XS w_grow quiet   --bank $I/bank$it.npz --out $I/p$it
    guard $I/p${it}_hits.parquet   inj_grow$it   $XS w_grow grow    --bank $I/bank$it.npz --out $I/p$it --thr $I/p${it}_thresholds.csv
    guard $I/p${it}_survival.csv   inj_surv$it   $XS w_grow survive --bank $I/bank$it.npz --out $I/p$it --seeds $I/cl_seeds.csv
  done
  guard $I/dedupe.csv   inj_dedupe python $S/x_stage.py w_dedupe --bank $I/bank2.npz --indep $I/p2_indep.parquet --out $I/dedupe.csv
  guard $I/lg_gates.csv inj_gates  python $S/x_locate_gates.py --bank $I/bank2.npz --indep $I/p2_indep.parquet \
        --dedupe $I/dedupe.csv --survival $I/p2_survival.csv --out $I/lg
  guard $I/recovery.csv inj_score  python $S/x_inject.py score --root $I --lg $I/lg
  unset X_CACHE
fi

# ================================================================ 9. candidate table ============
step candidates python $S/x_candidates.py --root $O

# ================================================================ marker =========================
STAGE="marker"
python - "$O" "$SMOKE" <<'PY' > "$MDONE"
import json, os, sys
import pandas as pd
O, SMOKE = sys.argv[1], sys.argv[2]


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


S = {}
try:
    S = json.load(open(f"{O}/step5_summary.json"))
except Exception:
    pass
print(("STEP5 SMOKE" if SMOKE == "1" else "STEP5 N-Olympic 2010 discovery"), "-- DONE",
      pd.Timestamp.now("UTC").isoformat())
F = json.load(open("results/freeze_wa_v1.json"))
print(f"bed: {F['sites']['n_sites']} sites / {F['sites']['n_detection_stations']} stations; "
      f"withheld {F['withheld']['stations']}; spread "
      f"{F['moveout']['spread_s_at_centre']:.2f} s; {F['sub_region']['n_active_cells']} cells")
print(f"T = {j(f'{O}/T.json','T')} ({j(f'{O}/T.json','binds')}), "
      f"{j(f'{O}/T.json','n_quiet_hours')} quiet hours")
print(f"tremor hours {j(f'{O}/ac_meta.json','n_hours')}, pairs {j(f'{O}/ac_meta.json','n_pairs')}, "
      f"{j(f'{O}/ac_meta.json','median_s_per_hour')} s/hour")
print(f"seed clusters {n(f'{O}/cl_seeds.csv')}; survivors "
      f"{n(f'{O}/p0_survival.csv','survives',True)}/{n(f'{O}/p0_survival.csv')}; "
      f"iterations to {open(f'{O}/last_iteration.txt').read().strip() if os.path.exists(f'{O}/last_iteration.txt') else '?'}; "
      f"deduped kept {n(f'{O}/dedupe.csv','kept',True)}")
print(f"located {j(f'{O}/lg.json','n_located')}; Gates 0-2 PASS {j(f'{O}/lg.json','n_gate_pass')}")
print(f"candidates {S.get('n_candidates')}; GATED FAMILIES {S.get('n_gated_families')}; "
      f"site-shift null {S.get('null_counts')} (mean {S.get('null_mean')}), "
      f"condition met {S.get('null_condition_met')}")
print(f"injection recovered {S.get('n_injected_recovered')}/{S.get('n_injected')}; "
      f"median location error {S.get('injection_median_loc_err_km')} km")
ch = (S.get("chestler") or {})
if ch.get("available"):
    print(f"Chestler nearest {ch['median_nearest_km']:.2f} km vs random-cell null "
          f"{ch['null_median_km']:.2f} [{ch['null_p05_km']:.2f}, {ch['null_p95_km']:.2f}] "
          f"km, p={ch['p_closer']:.3f} (descriptive)")
print(f"ACCEPT BAR: {S.get('accept_bar')}")
print(f"sealed 2011: {S.get('sealed_2011')}")
print(f"files: results/freeze_wa_v1.json, {O}/CANDIDATES.csv, {O}/step5_summary.json, {O}/")
PY
T1=$(date +%s)
echo "elapsed_s: $((T1-T0))" >> "$MDONE"
FINISHED=1
echo "=== [$(date -u +%FT%TZ)] STEP5 QUEUE COMPLETE in $((T1-T0)) s" >&2
