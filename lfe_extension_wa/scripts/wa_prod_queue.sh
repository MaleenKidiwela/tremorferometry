#!/usr/bin/env bash
# WA PRODUCTION queue (2026-10-03): results/WA_PROD_PLAN.md W0-W13 as amended by Merlin's review and the user's decisions
# of 2026-10-03 (see results/wa_prod/DEVIATIONS.md for every departure, with its reason).  Launch-and-exit:
#   nohup setsid -f bash scripts/wa_prod_queue.sh > logs/wa_prod_queue.log 2>&1
# Markers results/wa_prod/WA_PROD_{PID,HEARTBEAT,DONE,FAILED}; instance lock locks/wa_prod.lock; report results/WA_PROD_RESULTS.md.
# No downloads, NOTHING deleted (W13 writes a delete-CANDIDATE list only), no git; outputs only under results/wa_prod/,
# results/cut/wa_ext/ and the new cache20hz_prod/; wa_common.TRACE_ROOTS gets traces_dl2 appended once WA_DL2_DONE exists
# (user instruction; backup scripts/wa_common.py.bak_pre_dl2).  Every data-reading stage FAILS when > 3 % of its station-days
# error or are silently missing (wa_prod_common.guard).  GPU: waits for other projects' GPU locks / foreign GPU processes;
# at most 3 matched-filter passes at once; CUDA-12 libraries on LD_LIBRARY_PATH; cuFFT plan cache bounded + cleared daily.
# Smoke: WA_PROD_SMOKE=1 WA_PROD_DIR=<scratch> WA_PROD_CUT_DIR=<scratch>/cut WA_PROD_CACHE=<scratch>/cache
#        WA_PROD_REPORT=<scratch>/R.md WA_PROD_LOCK=locks/wa_prod_smoke.lock WA_PROD_TRUTH_DAYS=.. WA_PROD_PROD_DAYS=..
#        WA_PROD_SKIP_DL2=1 WA_PROD_FIT_MAX=.. WA_PROD_CUT_ARGS="--months YYYY-MM --max-days 1" WA_PROD_NSHARD=1
set -u -o pipefail
WS=/home/jovyan/tremorferometry/lfe_extension_wa
cd "$WS" || exit 1
export WA_PROD_DIR=${WA_PROD_DIR:-$WS/results/wa_prod}
export WA_PROD_CUT_DIR=${WA_PROD_CUT_DIR:-$WS/results/cut/wa_ext}
export WA_PROD_CACHE=${WA_PROD_CACHE:-$WS/cache20hz_prod}
export WA_PROD_REPORT=${WA_PROD_REPORT:-$WS/results/WA_PROD_RESULTS.md}
export WA_EXTRA_NODES=$WS/results/step4b/setup/extra_nodes.csv
O=$WA_PROD_DIR
NSH=${WA_PROD_NSHARD:-3}
SMOKE=${WA_PROD_SMOKE:-0}
LOCKF=${WA_PROD_LOCK:-locks/wa_prod.lock}
mkdir -p "$O" "$O/logs" "$WA_PROD_CUT_DIR" logs locks
LG=$O/logs
exec 9>$LOCKF
flock -n 9 || { echo "[wa_prod] another instance holds $LOCKF; exiting"; exit 0; }
echo $$ > $O/WA_PROD_PID
STAGE="init"; FINISHED=0; KIDS=""; HB=""
beat(){ echo "stage=$STAGE pid=$$ $(date -u +%FT%TZ) gpu=$(nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader 2>/dev/null | head -1)" > $O/WA_PROD_HEARTBEAT; }
hb_start(){ ( exec 9>&-; while sleep 120; do
                g=$(nvidia-smi --query-gpu=memory.used,memory.total,utilization.gpu --format=csv,noheader 2>/dev/null | head -1)
                echo "$(date -u +%FT%TZ) $g" >> $O/gpu_mem.log
                [ -f $O/WA_PROD_HEARTBEAT ] && touch $O/WA_PROD_HEARTBEAT
              done ) & HB=$!; }
killtree(){ local c; for c in $(pgrep -P "$1" 2>/dev/null); do killtree "$c"; done; kill "$1" 2>/dev/null; }
cleanup(){ [ -n "$HB" ] && kill "$HB" 2>/dev/null; for p in $KIDS; do killtree "$p"; done
           for p in $(pgrep -P $$ 2>/dev/null); do killtree "$p"; done; }
fail(){ echo "[wa_prod] FAILED stage=$STAGE rc=$1 $(date -u +%FT%TZ)"; echo "stage=$STAGE rc=$1 $(date -u +%FT%TZ)" > $O/WA_PROD_FAILED; cleanup; exit "$1"; }
onexit(){ rc=$?; cleanup; [ "$FINISHED" = 1 ] && return 0; [ -f $O/WA_PROD_FAILED ] && return 0
          echo "stage=$STAGE rc=$rc unexpected-exit $(date -u +%FT%TZ)" > $O/WA_PROD_FAILED; }
trap 'fail 143' TERM INT HUP
trap onexit EXIT
run(){ STAGE="$1"; shift; beat; echo "[wa_prod] === $STAGE $(date -u +%FT%TZ)"; "$@" 9>&- || fail $?; }
done_or(){ [ -f "$O/$1" ]; }
mark(){ date -u +%FT%TZ > "$O/$1"; }
bg(){ # bg <name> <cmd...>: background stage, markers BG_<name>_DONE / BG_<name>_FAILED
  local n=$1; shift
  if done_or BG_${n}_DONE; then echo "[wa_prod] bg $n already done"; return 0; fi
  rm -f $O/BG_${n}_FAILED
  ( exec 9>&-; set -o pipefail; "$@" > $LG/bg_$n.log 2>&1 && date -u +%FT%TZ > $O/BG_${n}_DONE || { echo "rc=$? $(date -u +%FT%TZ)" > $O/BG_${n}_FAILED; } ) &
  KIDS="$KIDS $!"; echo "[wa_prod] bg $n launched pid $! ($LG/bg_$n.log)"; }
bgwait(){ # bgwait <name>
  STAGE="wait_$1"; beat
  while ! done_or BG_$1_DONE; do
    [ -f $O/BG_$1_FAILED ] && { tail -15 $LG/bg_$1.log; fail 7; }
    sleep 60
  done; }
rm -f $O/WA_PROD_FAILED $O/WA_PROD_DONE
hb_start
export PYTHONPATH=$WS/scripts:${PYTHONPATH:-}
export PYTHONUNBUFFERED=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 PYTHONWARNINGS=ignore
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
PY=python
echo "[wa_prod] start $(date -u +%FT%TZ) OUT=$O CUT=$WA_PROD_CUT_DIR CACHE=$WA_PROD_CACHE smoke=$SMOKE"; df -h /home/jovyan | tail -1

gpu_wait(){ STAGE="gpu_wait"; beat
  while true; do
    held=0
    for f in /home/jovyan/tremorferometry/*/locks/*.lock /home/jovyan/tremorferometry/lfe_extension_wa/results/wa_pre/gpu.lock; do
      [ -e "$f" ] || continue
      case "$f" in *wa_prod*.lock) continue;; esac
      flock -n "$f" true 2>/dev/null || { held=1; echo "$(date -u +%FT%TZ) waiting for GPU lock $f"; }
    done
    mine=$(pgrep -d, -f "scripts/wa_prod_" || true)
    busy=$(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader,nounits 2>/dev/null | \
           awk -F', ' -v me=",$mine," 'index(me, ","$1",")==0 && $2>4096' | wc -l)
    [ "$held" = 0 ] && [ "$busy" -le 0 ] && break
    sleep 120
  done; }
gpu_par(){ # gpu_par <file> <tag>: <= 3 GPU jobs at once; all must succeed
  local f=$1 tag=$2; STAGE="$tag"; beat
  echo "[wa_prod] === $tag ($(wc -l < $f) GPU jobs) $(date -u +%FT%TZ)"
  rm -f $O/gpujob_${tag}_*.rc
  nl -ba -w1 -s' ' $f | xargs -P 3 -d '\n' -n 1 bash -c 'i=${0%% *}; cmd=${0#* }; ( eval "$cmd" ) > '"$LG"'/'"$tag"'_$i.log 2>&1; echo $? > '"$O"'/gpujob_'"$tag"'_$i.rc' 9>&-
  local bad=0
  [ "$(ls $O/gpujob_${tag}_*.rc 2>/dev/null | wc -l)" = "$(wc -l < $f)" ] || bad=1
  for r in $O/gpujob_${tag}_*.rc; do [ "$(cat $r)" = "0" ] || { bad=1; echo "[wa_prod] $tag $r rc=$(cat $r)"; }; done
  [ $bad = 0 ] || { tail -8 $LG/${tag}_*.log; fail 4; }; }

# ================================================================== W0 / W1 / W4a setup (frozen before any run)
if ! done_or SETUP_DONE; then
  run freeze    $PY scripts/wa_prod_setup.py freeze
  run roster    $PY scripts/wa_prod_setup.py roster
  run families  $PY scripts/wa_prod_setup.py families
  run truth     $PY scripts/wa_prod_setup.py truth
  run days      $PY scripts/wa_prod_setup.py days
  run daylists  $PY scripts/wa_prod_setup.py daylists
  run excl      $PY scripts/wa_prod_setup.py excl
  run rawinv    $PY scripts/wa_prod_setup.py rawinv
  run templates $PY scripts/wa_prod_setup.py templates --workers 12
  run core      $PY scripts/wa_prod_setup.py core
  mark SETUP_DONE
fi
# ================================================================== W3 caches (phase 1 now; phase 2 = traces_dl2 after WA_DL2_DONE)
done_or DET_DONE || { run determinism $PY scripts/wa_prod_cache.py determinism; mark DET_DONE; }
done_or PLAN1_DONE || { run cache_plan1 $PY scripts/wa_prod_cache.py plan --phase 1; mark PLAN1_DONE; }
bg cache1 $PY scripts/wa_prod_cache.py run --phase 1 --workers ${WA_PROD_CACHE_WORKERS:-12}
dl2_phase(){
  if [ "${WA_PROD_SKIP_DL2:-0}" = 1 ]; then
    $PY - <<EOF
import wa_prod_common as P, os
P.ensure(f"{P.OUT}/cache")
for m in sorted({d[:7] for d in P.days()}):
    if m >= "2012-01": P.marker(f"cache/done_p2_{m}", "smoke: dl2 skipped")
EOF
    return 0
  fi
  t0=$(date +%s)
  until [ -f $WS/results/wa_dl2/WA_DL2_DONE ] || [ -f $WS/results/wa_dl2/STOPPED_STORAGE ]; do
    [ $(( $(date +%s) - t0 )) -gt 172800 ] && { echo "dl2 not done after 48 h"; return 3; }
    sleep 300
  done
  # user instruction: append traces_dl2 to wa_common.TRACE_ROOTS once WA_DL2_DONE exists
  if ! grep -q 'TRACES_DL2' scripts/wa_common.py; then
    cp -n scripts/wa_common.py scripts/wa_common.py.bak_pre_dl2
    $PY - <<'EOF'
p = "scripts/wa_common.py"; s = open(p).read()
old = 'TRACE_ROOTS = [f"{ROOT}/lfe_array_discovery/traces", f"{ROOT}/lfe_extension/traces",\n               f"{ROOT}/lfe_catalog_repair/traces", TRACES, TRACES_DL]'
new = ('TRACES_DL2 = f"{STUDY}/traces_dl2"   # tiers C, D, D2 + temporary networks (scripts/wa_dl2.py); appended once WA_DL2_DONE\n'
       '                                     # existed (WA production queue, 2026-10-03); only days missing from the roots above\n'
       'TRACE_ROOTS = [f"{ROOT}/lfe_array_discovery/traces", f"{ROOT}/lfe_extension/traces",\n'
       '               f"{ROOT}/lfe_catalog_repair/traces", TRACES, TRACES_DL, TRACES_DL2]')
assert old in s, "TRACE_ROOTS line not found"
open(p, "w").write(s.replace(old, new))
print("traces_dl2 appended to wa_common.TRACE_ROOTS")
EOF
  fi
  $PY scripts/wa_prod_setup.py rawinv && $PY scripts/wa_prod_cache.py plan --phase 2 && \
  $PY scripts/wa_prod_cache.py run --phase 2 --workers ${WA_PROD_CACHE2_WORKERS:-8}
}
bg cache2 dl2_phase
# ================================================================== W5 truth runs (all networks, 2 scrambled + time-reversed controls)
STAGE=wait_truth_caches; beat
until $PY -c "
import wa_prod_common as P, os, sys
sys.exit(0 if all(os.path.exists(f'{P.OUT}/cache/done_p1_{m}') for m in P.TRUTH_MONTHS) else 1)"; do
  [ -f $O/BG_cache1_FAILED ] && { tail -15 $LG/bg_cache1.log; fail 7; }
  sleep 60
done
if ! done_or TRUTH_MF_DONE; then
  gpu_wait
  : > $O/gpu_truth.txt
  echo "$PY scripts/wa_prod_mf.py --run truth --days-file $O/setup/days_truth.json --nets E3,all --series data,c1,c2,rev" >> $O/gpu_truth.txt
  echo "$PY scripts/wa_prod_mf.py --run truth_e12 --days-file $O/setup/days_truth.json --nets E1,E2 --series data,c1,rev" >> $O/gpu_truth.txt
  gpu_par $O/gpu_truth.txt truthmf
  run guard_mf_truth $PY -c "import wa_prod_common as P; P.guard('mf_truth', 'mf_truth_s*'); P.guard('mf_truth_e12', 'mf_truth_e12_s*')"
  mark TRUTH_MF_DONE
fi
# ================================================================== W6 production pass (W4a templates; core + all networks) in the background
gpu_wait
: > $O/gpu_prod.txt
for i in $(seq 0 $((NSH-1))); do
  echo "$PY scripts/wa_prod_mf.py --run prod --days-file $O/setup/days_prod.json --nets core,all --series data,c1,rev --shard $i --nshard $NSH --wait-cache" >> $O/gpu_prod.txt
done
prodmf(){ local bad=0; rm -f $O/gpujob_prodmf_*.rc
  nl -ba -w1 -s' ' $O/gpu_prod.txt | xargs -P 3 -d '\n' -n 1 bash -c 'i=${0%% *}; cmd=${0#* }; ( eval "$cmd" ) > '"$LG"'/prodmf_$i.log 2>&1; echo $? > '"$O"'/gpujob_prodmf_$i.rc'
  for r in $O/gpujob_prodmf_*.rc; do [ "$(cat $r)" = "0" ] || bad=1; done
  [ $bad = 0 ] && $PY -c "import wa_prod_common as P; P.guard('mf_prod', 'mf_prod_s*')" && date -u +%FT%TZ > $O/mf/prod/MF_DONE; }
mkdir -p $O/mf/prod
bg prodmf prodmf
# ================================================================== W7 / W8 confidence on the truth runs (CPU; PNET CF on the GPU, small)
done_or HF_TRUTH_DONE || { run hf_truth $PY scripts/wa_prod_hf.py --run truth --workers 8
                           run hf_truth_e12 $PY scripts/wa_prod_hf.py --run truth_e12 --workers 8
                           run guard_hf_truth $PY -c "import wa_prod_common as P; P.guard('hf', 'hf_truth*')"; mark HF_TRUTH_DONE; }
bg hfprod $PY scripts/wa_prod_hf.py --run prod --workers ${WA_PROD_HF_WORKERS:-8} --wait-mf $O/mf/prod/MF_DONE
done_or PNET_TRUTH_DONE || { run pnet_truth $PY scripts/wa_prod_pnet.py --run truth
                             run guard_pnet_truth $PY -c "import wa_prod_common as P; P.guard('pnet_truth', 'pnet_truth_s*')"; mark PNET_TRUTH_DONE; }
export OMP_NUM_THREADS=12
done_or FIT_DONE || { run fit $PY scripts/wa_prod_fit.py fit; mark FIT_DONE; }
done_or SEALED_DONE || { run open_sealed $PY scripts/wa_prod_fit.py open_sealed; mark SEALED_DONE; }
export OMP_NUM_THREADS=2
# ================================================================== production P
bgwait prodmf
bgwait hfprod
run guard_hf_prod $PY -c "import wa_prod_common as P; P.guard('hf', 'hf_prod')"
ADOPT=$($PY -c "import json; print(int(bool(json.load(open('$O/conf/pnet_gate.json')).get('adopt'))))" 2>/dev/null || echo 0)
if [ "$ADOPT" = 1 ] && ! done_or PNET_PROD_DONE; then
  gpu_wait; : > $O/gpu_pnetp.txt
  for i in $(seq 0 $((NSH-1))); do echo "$PY scripts/wa_prod_pnet.py --run prod --shard $i --nshard $NSH" >> $O/gpu_pnetp.txt; done
  gpu_par $O/gpu_pnetp.txt pnetprod
  run guard_pnet_prod $PY -c "import wa_prod_common as P; P.guard('pnet_prod', 'pnet_prod_s*')"
  mark PNET_PROD_DONE
fi
export OMP_NUM_THREADS=16
done_or APPLY_PROD_DONE || { run apply_prod $PY scripts/wa_prod_fit.py apply --run prod; mark APPLY_PROD_DONE; }
done_or REFIT_RISK_DONE || { run refit_risk $PY scripts/wa_prod_fit.py refit_risk; mark REFIT_RISK_DONE; }
export OMP_NUM_THREADS=2
# ================================================================== W4b re-templates + T2B new nodes, then the rich pass
bgwait cache2
run guard_cache $PY -c "import wa_prod_common as P; P.guard('cache_p1', 'cache_p1_*'); P.guard('cache_p2', 'cache_p2_*', fail=False)"
done_or W4B_DONE || { run w4b_targets $PY scripts/wa_prod_retpl.py targets
                      run w4b_build $PY scripts/wa_prod_retpl.py build --workers 12
                      run w4b_richdays $PY scripts/wa_prod_retpl.py richdays; mark W4B_DONE; }
NRICH=$($PY -c "import json; print(len(json.load(open('$O/templates/rich_days.json'))))" 2>/dev/null || echo 0)
if [ "$NRICH" -gt 0 ] && ! done_or RICH_MF_DONE; then
  gpu_wait; : > $O/gpu_rich.txt
  for i in $(seq 0 $((NSH-1))); do
    echo "$PY scripts/wa_prod_mf.py --run rich --days-file $O/setup/days_prod.json --rich --rich-days $O/templates/rich_days.json --nets rich --series data,c1,rev --shard $i --nshard $NSH --wait-cache" >> $O/gpu_rich.txt
  done
  gpu_par $O/gpu_rich.txt richmf
  run guard_mf_rich $PY -c "import wa_prod_common as P; P.guard('mf_rich', 'mf_rich_s*')"
  mark RICH_MF_DONE
fi
if [ "$NRICH" -gt 0 ]; then
  done_or HF_RICH_DONE || { run hf_rich $PY scripts/wa_prod_hf.py --run rich --workers 14
                            run guard_hf_rich $PY -c "import wa_prod_common as P; P.guard('hf', 'hf_rich')"; mark HF_RICH_DONE; }
  export OMP_NUM_THREADS=16
  done_or APPLY_RICH_DONE || { run apply_rich $PY scripts/wa_prod_fit.py apply --run rich; mark APPLY_RICH_DONE; }
  export OMP_NUM_THREADS=2
fi
# ================================================================== W9 / W10 / W11
done_or POST_DONE || { run tremor $PY scripts/wa_prod_post.py tremor
                       run member $PY scripts/wa_prod_post.py member
                       run cutoff $PY scripts/wa_prod_post.py cutoff
                       run rates  $PY scripts/wa_prod_post.py rates; mark POST_DONE; }
# ================================================================== W12 windows (pilot month first; size limit 3 TB, pre-registered)
if ! done_or CUT_DONE; then
  STAGE=cut_pilot; beat
  $PY scripts/wa_prod_cut.py pilot --workers 8 ${WA_PROD_PILOT_ARGS:-} > $LG/cut_pilot.log 2>&1 9>&-; prc=$?
  if [ $prc = 0 ]; then
    run cut        $PY scripts/wa_prod_cut.py cut --workers ${WA_PROD_CUT_WORKERS:-10} ${WA_PROD_CUT_ARGS:-}
    run cut_guard  $PY scripts/wa_prod_cut.py guard
    run cut_verify $PY scripts/wa_prod_cut.py verify --n 20
    STAGE=cut_stackcheck; beat
    $PY scripts/wa_prod_cut.py stackcheck > $LG/cut_stackcheck.log 2>&1 9>&- || echo "[wa_prod] stack check did not pass (reported, not fatal)"
    mark CUT_DONE
  elif [ $prc = 6 ]; then
    echo "[wa_prod] cut pilot: size or verify bar not met -> full cut NOT run (decision for the user)"; mark CUT_DECISION_NEEDED
  else
    tail -15 $LG/cut_pilot.log; STAGE=cut_pilot; fail $prc
  fi
fi
# ================================================================== W13 delete-candidate list (nothing deleted), verify, report
done_or DELETE_LIST_DONE || { run delete_list $PY scripts/wa_prod_delete.py; mark DELETE_LIST_DONE; }
run verify_inputs $PY scripts/wa_prod_setup.py verify
run report $PY scripts/wa_prod_report.py
date -u +%FT%TZ > $O/WA_PROD_DONE
STAGE="done"; FINISHED=1; beat
echo "[wa_prod] DONE $(date -u +%FT%TZ)"
