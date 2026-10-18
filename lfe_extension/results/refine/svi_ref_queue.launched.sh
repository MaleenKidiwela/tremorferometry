#!/usr/bin/env bash
# SVI_REF (2026-10-02): southern Vancouver Island refinement queue -- results/SVI_REFINE_PLAN.md as AMENDED by Merlin.
# Merlin's order: R0 -> R1 -> R8 pilot -> R2 -> R3 (-> R4 if gated) -> R8 full (parallel) -> PNET gate -> R5 -> inter-ETS
# validation -> R7 -> R9 -> report.  R3's GPU runs start while R2's CPU pass runs (no dependency between them); R2 must finish
# before the production pass (it supplies the effective station sets).
# No downloads; NOTHING deleted (R9 writes a delete-CANDIDATE list only); no git; v1.1 never overwritten (outputs under
# results/refine/ and results/cut/svi_ext/, rebuilt caches under cache20hz_r1/).  Every data-reading stage FAILS (marker +
# exit != 0) when > 3 % of its station-days error or are silently missing (svi_ref_common.guard).
# Markers results/SVI_REF_{HEARTBEAT,DONE,FAILED,PID}; lock locks/svi_ref.lock; log logs/svi_ref_queue.log (caller redirects);
# GPU memory every 2 min -> results/refine/gpu_mem.log; report results/SVI_REF_RESULTS.md.
# Smoke: SVI_REF_TAG=SMOKE SVI_REF_DIR=<scratch> SVI_REF_CUT_DIR=<scratch>/cut SVI_REF_REPORT=<scratch>/R.md SVI_REF_CACHE_R1=<scratch>/c1
#        SVI_REF_SEALED_YEAR=2011 SVI_REF_PROD_DAYS=.. SVI_REF_TRUTH_DAYS_<Y>=.. SVI_REF_R3M_DAYS=.. SVI_REF_S13_DAYS=..
#        SVI_REF_R2_ARGS="--months 2021-10 --stations LRIV,B926" SVI_REF_PILOT_ARGS="--max-days 1" SVI_REF_CUT_ARGS="--months 2021-10 --max-days 1"
set -u -o pipefail
WS=/home/jovyan/tremorferometry/lfe_extension
cd "$WS" || exit 1
TAGM=${SVI_REF_TAG:-SVI_REF}
export SVI_REF_DIR=${SVI_REF_DIR:-$WS/results/refine}
export SVI_REF_CUT_DIR=${SVI_REF_CUT_DIR:-$WS/results/cut/svi_ext}
export SVI_REF_REPORT=${SVI_REF_REPORT:-$WS/results/SVI_REF_RESULTS.md}
export SVI_REF_SEALED_YEAR=${SVI_REF_SEALED_YEAR:-2013}
O=$SVI_REF_DIR
MK=results; [ "$TAGM" != "SVI_REF" ] && MK=$O
mkdir -p "$O" "$SVI_REF_CUT_DIR" logs locks "$O/logs"
LG=$O/logs
exec 9>locks/${TAGM,,}.lock
flock -n 9 || { echo "[ref] another instance holds the lock; exiting"; exit 0; }
echo $$ > $MK/${TAGM}_PID
STAGE="init"; FINISHED=0; HB=""; KIDS=""
beat(){ echo "stage=$STAGE pid=$$ $(date -u +%FT%TZ) gpu=$(nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader 2>/dev/null | head -1)" > $MK/${TAGM}_HEARTBEAT; }
hb_start(){ ( exec 9>&-; while sleep 120; do
                g=$(nvidia-smi --query-gpu=memory.used,memory.total,utilization.gpu --format=csv,noheader 2>/dev/null | head -1)
                echo "$(date -u +%FT%TZ) $g" >> $O/gpu_mem.log
                [ -f $MK/${TAGM}_HEARTBEAT ] && sed -i "s#gpu=.*#gpu=$g#" $MK/${TAGM}_HEARTBEAT && touch $MK/${TAGM}_HEARTBEAT
              done ) & HB=$!; }
killtree(){ local c; for c in $(pgrep -P "$1" 2>/dev/null); do killtree "$c"; done; kill "$1" 2>/dev/null; }
cleanup(){ [ -n "$HB" ] && kill "$HB" 2>/dev/null; for p in $KIDS; do killtree "$p"; done
           for p in $(pgrep -P $$ 2>/dev/null); do killtree "$p"; done
           sleep 2; pkill -s 0 -f "python scripts/svi_ref_" 2>/dev/null; }       # orphaned pool workers of this session
fail(){ echo "[ref] FAILED stage=$STAGE rc=$1 $(date -u +%FT%TZ)"; echo "stage=$STAGE rc=$1 $(date -u +%FT%TZ)" > $MK/${TAGM}_FAILED; cleanup; exit "$1"; }
onexit(){ rc=$?; cleanup; [ "$FINISHED" = 1 ] && return 0; [ -f $MK/${TAGM}_FAILED ] && return 0
          echo "[ref] EXITED unexpectedly stage=$STAGE rc=$rc $(date -u +%FT%TZ)"
          echo "stage=$STAGE rc=$rc unexpected-exit $(date -u +%FT%TZ)" > $MK/${TAGM}_FAILED; }
trap 'fail 143' TERM INT HUP
trap onexit EXIT
run(){ STAGE="$1"; shift; beat; echo "[ref] === $STAGE $(date -u +%FT%TZ)"; "$@" 9>&- || fail $?; }
done_or(){ [ -f "$O/$1" ]; }
mark(){ date -u +%FT%TZ > "$O/$1"; }
rm -f $MK/${TAGM}_FAILED $MK/${TAGM}_DONE $O/R2_FAILED $O/PREP_FAILED $O/CUT_FAILED
hb_start
export PYTHONUNBUFFERED=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
PY=python
echo "[ref] start $(date -u +%FT%TZ) OUT=$O CUT=$SVI_REF_CUT_DIR sealed=$SVI_REF_SEALED_YEAR"; df -h /home/jovyan | tail -1

# GPU: never share the card with another GPU queue (DUC_FIN must be done; no foreign compute process holding > 4 GB)
gpu_wait(){ STAGE="gpu_wait"; beat
  while true; do
    [ -f /home/jovyan/tremorferometry/lfe_south/results/DUC_FIN_DONE ] || { sleep 120; continue; }
    busy=$(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader,nounits 2>/dev/null | awk -v me="$KIDS" '$2>4096' | wc -l)
    [ "$busy" -le 0 ] && break
    [ "${1:-}" = "own" ] && break
    sleep 120
  done; }
# up to 3 matched-filter passes at once on the H100 (Track A-2 OOM): every line of $1 is one command; all must succeed
gpu_par(){ local f=$1; local tag=$2
  STAGE="$tag"; beat; echo "[ref] === $tag ($(wc -l < $f) GPU jobs, <= 3 at once) $(date -u +%FT%TZ)"
  rm -f $O/gpujob_${tag}_*.rc
  nl -ba -w1 -s' ' $f | xargs -P 3 -d '\n' -n 1 bash -c 'i=${0%% *}; cmd=${0#* }; ( eval "$cmd" ) > '"$LG"'/'"$tag"'_$i.log 2>&1; echo $? > '"$O"'/gpujob_'"$tag"'_$i.rc' 9>&-
  local bad=0
  [ "$(ls $O/gpujob_${tag}_*.rc 2>/dev/null | wc -l)" = "$(wc -l < $f)" ] || { echo "[ref] $tag: missing job exit codes"; bad=1; }
  for r in $O/gpujob_${tag}_*.rc; do [ "$(cat $r)" = "0" ] || { bad=1; echo "[ref] $tag job $r rc=$(cat $r)"; }; done
  [ $bad = 0 ] || { tail -5 $LG/${tag}_*.log; fail 4; }
}
mfcmd(){ # run kind year sset clean daysfile shard nshard [extra]
  echo "$PY scripts/svi_ref_mf.py --run $1 --kind $2 --year $3 --sset $4 --clean $5 --days-file $6 --shard $7 --nshard $8 ${9:-}"; }
guard_runs(){ for r in "$@"; do run guard_$r $PY scripts/svi_ref_common.py guard mf_$r "mf_${r}_s*"; done; }

# ---- 0. frozen day table (before any fit), R0
run days $PY scripts/svi_ref_common.py days
run r0_watcher $PY scripts/svi_ref_r0.py check --wait-s ${SVI_REF_R0_WAIT:-1200}
run r0_snapshot $PY scripts/svi_ref_r0.py snapshot
# ---- 1. R1 completeness (+ rebuild of failing cache days into cache20hz_r1)
done_or R1_DONE || { run r1 $PY scripts/svi_ref_r1.py --workers 16; mark R1_DONE; }
run r3_months $PY scripts/svi_ref_r3.py months
run daylists $PY scripts/svi_ref_common.py daylists
# ---- 2. R8 pilot (one busy month, 20/20 verify, stack check, size extrapolation)
PM=$($PY scripts/svi_ref_common.py busiest)
done_or PILOT_DONE || { run r8_pilot $PY scripts/svi_ref_cut.py pilot --month $PM --workers 6 ${SVI_REF_PILOT_ARGS:-}; mark PILOT_DONE; }
CUT_OK=$($PY -c "import json; print(int(json.load(open('$O/cut_pilot/PILOT.json'))['extrapolation']['pass_lt_0p5TB']))")
echo "[ref] pilot month $PM; extrapolated cut < 0.5 TB: $CUT_OK"
# ---- 3. CPU in the background: R2 (raw health + stacks, guard, A-2 reproduction, summary) || cleaning prep
R2PID=""; PPID_=""
if ! done_or R2_DONE; then
  ( exec 9>&-; set -o pipefail
    { $PY scripts/svi_ref_r2.py run --workers ${SVI_REF_R2_WORKERS:-18} ${SVI_REF_R2_ARGS:-} && \
      $PY scripts/svi_ref_common.py guard r2raw r2raw_all && $PY scripts/svi_ref_common.py guard r2cache r2cache_all && \
      $PY scripts/svi_ref_r2.py pc_a2 && $PY scripts/svi_ref_r2.py summary --workers 12; } > $LG/r2.log 2>&1 \
      && mark R2_DONE || { echo "rc=$? $(date -u +%FT%TZ)" > $O/R2_FAILED; } ) &
  R2PID=$!; KIDS="$KIDS $R2PID"; echo "[ref] R2 launched pid $R2PID ($LG/r2.log)"
fi
done_or PREP1_DONE || { run prep_first $PY scripts/svi_ref_prep.py --days-file $O/days/prep_first.json --workers 10 --tag first; \
                        run prep_first_guard $PY scripts/svi_ref_common.py guard prep_first prep_first; mark PREP1_DONE; }
if ! done_or PREP_DONE; then
  ( exec 9>&-; { $PY scripts/svi_ref_prep.py --days-file $O/days/prep_rest.json --workers ${SVI_REF_PREP_WORKERS:-8} --tag rest && \
      $PY scripts/svi_ref_common.py guard prep_rest prep_rest; } > $LG/prep_rest.log 2>&1 && mark PREP_DONE \
      || { echo "rc=$? $(date -u +%FT%TZ)" > $O/PREP_FAILED; } ) &
  PPID_=$!; KIDS="$KIDS $PPID_"; echo "[ref] prep (rest) launched pid $PPID_"
fi
# ---- 4. R3 A/B on the GPU (truth years 2010-2012 + flagged months), both variants
gpu_wait
if ! done_or R3_DONE; then
  : > $O/gpu_r3.txt
  for v in nc cl; do c=$([ $v = cl ] && echo 1 || echo 0)
    mfcmd r3t2010_$v truth 2010 S9 $c $O/days/truth_2010.json 0 1 >> $O/gpu_r3.txt
    mfcmd r3t2011_$v truth 2011 S15 $c $O/days/truth_2011.json 0 1 >> $O/gpu_r3.txt
    mfcmd r3t2012_$v truth 2012 S15 $c $O/days/truth_2012.json 0 1 >> $O/gpu_r3.txt
    mfcmd r3m_$v prod 0 eff $c $O/days/r3m.json 0 1 >> $O/gpu_r3.txt
  done
  gpu_par $O/gpu_r3.txt r3mf
  guard_runs r3t2010_nc r3t2010_cl r3t2011_nc r3t2011_cl r3t2012_nc r3t2012_cl r3m_nc r3m_cl
  run r3_eval $PY scripts/svi_ref_r3.py evaluate
  mark R3_DONE
fi
VAR=$($PY -c "import json; print('cl' if json.load(open('$O/r3/r3_decision.json'))['adopt_cleaning'] else 'nc')")
CL=$([ $VAR = cl ] && echo 1 || echo 0)
echo "[ref] R3 decision: variant $VAR (cleaning applied = $CL); R4 = the production pass with this variant"
# ---- 5. truth runs with the chosen variant: S14 / S13 on 2011-2012, S15 / S14 / S13 on the sealed year (no labels read here)
if ! done_or TRUTH_DONE; then
  : > $O/gpu_truth.txt
  for y in 2011 2012; do for s in S14 S13; do mfcmd t${y}_${s}_$VAR truth $y $s $CL $O/days/truth_$y.json 0 1 >> $O/gpu_truth.txt; done; done
  for s in S15 S14 S13; do mfcmd t2013_${s}_$VAR truth 2013 $s $CL $O/days/truth_2013.json 0 1 >> $O/gpu_truth.txt; done
  gpu_par $O/gpu_truth.txt truthmf
  guard_runs t2011_S14_$VAR t2011_S13_$VAR t2012_S14_$VAR t2012_S13_$VAR t2013_S15_$VAR t2013_S14_$VAR t2013_S13_$VAR
  mark TRUTH_DONE
fi
# ---- 6. wait for R2 (effective station sets) and the prep of every production day (cleaning or flags-as-columns)
STAGE=wait_r2_prep; beat; echo "[ref] === waiting for R2 / prep $(date -u +%FT%TZ)"
[ -n "$R2PID" ] && wait $R2PID
[ -f $O/R2_FAILED ] && { tail -8 $LG/r2.log; STAGE=r2; fail 6; }
[ -n "$PPID_" ] && wait $PPID_
[ -f $O/PREP_FAILED ] && { tail -8 $LG/prep_rest.log; STAGE=prep_rest; fail 6; }
done_or R2_DONE && done_or PREP_DONE || { STAGE=r2_prep_markers; fail 6; }
# ---- 7. production pass (= R4 when the cleaning is adopted) with the R2 effective station sets, + the fixed-S13 re-run
gpu_wait
if ! done_or PROD_DONE; then
  : > $O/gpu_prod.txt
  for i in 0 1 2; do mfcmd prod prod 0 eff $CL $O/days/prod.json $i 3 "--excl-file $O/r2/r2_exclusions.csv" >> $O/gpu_prod.txt; done
  gpu_par $O/gpu_prod.txt prodmf
  guard_runs prod
  mark PROD_DONE
fi
# R8 full cut in the background (needs the production candidates for the v1.1 U v2 union)
CPID=""
if [ "$CUT_OK" = 1 ] && ! done_or CUT_DONE; then
  ( exec 9>&-; set -o pipefail
    { $PY scripts/svi_ref_cut.py cut --workers ${SVI_REF_CUT_WORKERS:-8} ${SVI_REF_CUT_ARGS:-} && $PY scripts/svi_ref_cut.py guard && \
      $PY scripts/svi_ref_cut.py verify --n 20 && $PY scripts/svi_ref_cut.py stackcheck; } > $LG/cut.log 2>&1 \
      && mark CUT_DONE || { echo "rc=$? $(date -u +%FT%TZ)" > $O/CUT_FAILED; } ) &
  CPID=$!; KIDS="$KIDS $CPID"; echo "[ref] R8 full cut launched pid $CPID ($LG/cut.log)"
elif [ "$CUT_OK" != 1 ]; then echo "[ref] R8 full cut NOT run: the pilot extrapolates to >= 0.5 TB (decision for the user)"; echo skipped > $O/CUT_SKIPPED; fi
if ! done_or S13_DONE; then
  : > $O/gpu_s13.txt
  for i in 0 1 2; do mfcmd prod_S13 prod 0 S13 $CL $O/days/prodS13.json $i 3 >> $O/gpu_s13.txt; done
  gpu_par $O/gpu_s13.txt s13mf
  guard_runs prod_S13
  mark S13_DONE
fi
# ---- 8. PNET gate features (2011 + sealed year only)
cat > $O/conf_runs.tmp.py <<EOF
import json, os
v = "$VAR"; o = "$O"
cfg = dict(variant=v, prod="prod", prod_S13="prod_S13",
           truth={"2010": {"S9": f"r3t2010_{v}"},
                  "2011": {"S15": f"r3t2011_{v}", "S14": f"t2011_S14_{v}", "S13": f"t2011_S13_{v}"},
                  "2012": {"S15": f"r3t2012_{v}", "S14": f"t2012_S14_{v}", "S13": f"t2012_S13_{v}"},
                  "2013": {"S15": f"t2013_S15_{v}", "S14": f"t2013_S14_{v}", "S13": f"t2013_S13_{v}"}})
os.makedirs(f"{o}/conf", exist_ok=True); json.dump(cfg, open(f"{o}/conf/runs.json", "w"), indent=1)
EOF
run runs_cfg $PY $O/conf_runs.tmp.py
if ! done_or PNET_GATE_FEAT_DONE; then
  PR="r3t2011_$VAR,t2011_S14_$VAR,t2011_S13_$VAR,t2013_S15_$VAR,t2013_S14_$VAR,t2013_S13_$VAR"
  [ "$SVI_REF_SEALED_YEAR" = 2011 ] && PR="r3t2012_$VAR,t2012_S14_$VAR,t2012_S13_$VAR,r3t2011_$VAR,t2011_S14_$VAR,t2011_S13_$VAR"
  printf "%s\n%s\n" "$PY scripts/svi_ref_pnet.py --runs $PR --shard 0 --nshard 2 --tag gate" \
                    "$PY scripts/svi_ref_pnet.py --runs $PR --shard 1 --nshard 2 --tag gate" > $O/gpu_pnetg.txt
  gpu_par $O/gpu_pnetg.txt pnetgate
  run pnet_gate_guard $PY scripts/svi_ref_common.py guard pnet_gate "pnet_gate_s*"
  mark PNET_GATE_FEAT_DONE
fi
# ---- 9. R5: fit + freeze, the ONE sealed look (also the PNET gate), production PNET only if the gate passes, apply
export OMP_NUM_THREADS=8
done_or FIT_DONE || { run fit $PY scripts/svi_ref_fit.py fit; mark FIT_DONE; }
done_or SEALED_DONE || { run open_sealed $PY scripts/svi_ref_fit.py open_sealed; mark SEALED_DONE; }
export OMP_NUM_THREADS=2
ADOPT=$($PY -c "import json; print(int(bool(json.load(open('$O/conf/pnet_gate.json')).get('adopt'))))")
if [ "$ADOPT" = 1 ] && ! done_or PNET_PROD_DONE; then
  gpu_wait
  printf "%s\n" "$PY scripts/svi_ref_pnet.py --runs prod --shard 0 --nshard 3 --tag prod" "$PY scripts/svi_ref_pnet.py --runs prod --shard 1 --nshard 3 --tag prod" \
               "$PY scripts/svi_ref_pnet.py --runs prod --shard 2 --nshard 3 --tag prod" > $O/gpu_pnetp.txt
  gpu_par $O/gpu_pnetp.txt pnetprod
  run pnet_prod_guard $PY scripts/svi_ref_common.py guard pnet_prod "pnet_prod_s*"
  mark PNET_PROD_DONE
fi
export OMP_NUM_THREADS=8
done_or APPLY_DONE || { run apply $PY scripts/svi_ref_fit.py apply; mark APPLY_DONE; }
# ---- 10. inter-ETS validation, R7 rates
done_or VAL_DONE || { run validation $PY scripts/svi_ref_val.py; mark VAL_DONE; }
done_or RATES_DONE || { run rates $PY scripts/svi_ref_rates.py; mark RATES_DONE; }
export OMP_NUM_THREADS=2
# ---- 11. wait for the cut; 2013 Bostock windows (only now that the sealed check is done); R9 delete-CANDIDATE list
STAGE=wait_cut; beat
[ -n "$CPID" ] && wait $CPID
[ -f $O/CUT_FAILED ] && { tail -8 $LG/cut.log; STAGE=cut; fail 6; }
if [ "$CUT_OK" = 1 ] && [ "$SVI_REF_SEALED_YEAR" = 2013 ] && ! done_or CUT2013_DONE; then
  run cut_bostock2013 $PY scripts/svi_ref_cut.py cut --bostock2013 --workers 6
  run cut_bostock2013_verify $PY scripts/svi_ref_cut.py verify --bostock2013 --n 20
  mark CUT2013_DONE
fi
done_or R9_DONE || { run r9_delete_candidates $PY scripts/svi_ref_delete.py; mark R9_DONE; }
# ---- 12. v1.1 integrity, report
run r0_verify $PY scripts/svi_ref_r0.py verify
run report $PY scripts/svi_ref_report.py
date -u +%FT%TZ > $MK/${TAGM}_DONE
STAGE="done"; FINISHED=1; beat
echo "[ref] DONE $(date -u +%FT%TZ)"
