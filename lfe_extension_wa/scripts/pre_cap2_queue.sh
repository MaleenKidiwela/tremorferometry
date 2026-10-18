#!/bin/bash
# NW dominance-cap repeat (2026-10-03): V2008 (test half) + A2011, 21-core network, uncapped / cap 0.3 / dominant-site-out,
# 2 scrambled controls + time-reversed templates.  Rules: results/wa_pre/cap2/PREREGISTERED.md (written before this ran).
# Outputs only under results/wa_pre/cap2/ (+ the lock file locks/wa_cap2.lock); never reads traces_dl*/traces_sealed;
# no downloads, deletions or git.  GPU block: waits while any project GPU lock (*/locks/*.lock, results/wa_pre/gpu.lock)
# is held, then holds locks/wa_cap2.lock (visible to the other queues' */locks/*.lock check).
set -u
cd /home/jovyan/tremorferometry/lfe_extension_wa
R=results/wa_pre/cap2
LOG=logs/wa_pre_cap2.log
echo $$ > $R/CAP2_PID
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')${LD_LIBRARY_PATH:-}
export PYTHONWARNINGS=ignore
hb() { echo "$(date -u +%FT%TZ) $1" > $R/CAP2_HEARTBEAT; echo "$(date -u +%FT%TZ) [cap2] $1" >> $LOG; }
fail() { echo "$(date -u +%FT%TZ) stage $1 failed (rc $2)" > $R/CAP2_FAILED; hb "FAILED $1 rc $2"; exit 1; }
gpu_wait() {
  while true; do
    held=0
    for f in /home/jovyan/tremorferometry/*/locks/*.lock /home/jovyan/tremorferometry/lfe_extension_wa/results/wa_pre/gpu.lock; do
      [ -e "$f" ] || continue
      case "$f" in *wa_cap2.lock) continue;; esac
      flock -n "$f" true 2>/dev/null || { held=1; echo "$(date -u +%FT%TZ) waiting for $f" >> $LOG; }
    done
    [ $held -eq 0 ] && break
    sleep 120
  done
}
run() {  # run <stage> <cmd...>
  st=$1; shift
  [ -e $R/stage_$st.done ] && { hb "skip $st (done)"; return 0; }
  hb "start $st"
  "$@" >> $LOG 2>&1; rc=$?
  [ $rc -eq 0 ] || fail $st $rc
  touch $R/stage_$st.done; hb "done $st"
}
PR="python scripts/pre_run.py --res $R --setup results/wa_gonogo/setup mfcap --"
MFARGS="--tpl build --n-control 2 --drop-dead 0.05 --all-templates"
REVARGS="--tpl buildrev --no-control --drop-dead 0.05 --all-templates"
C2="python scripts/pre_cap2.py"

# ---- setup: templates (symlinks) and the V2010 days of the A2011 score runs (symlinks to the existing V2010 runs,
#      needed only so A2011 is scored exactly as the go/no-go: HF cuts from the V2010 build half)
mkdir -p $R/templates
for f in tpl_build.npz tpl_qc_build.csv; do [ -e $R/templates/$f ] || ln -s ../../../templates/$f $R/templates/$f; done
for f in tpl_buildrev.npz tpl_qc_buildrev.csv; do [ -e $R/templates/$f ] || ln -s ../../ctl/res/templates/$f $R/templates/$f; done
link_v2010() {  # link_v2010 <srcdir> <dstdir>
  mkdir -p $2
  for p in $1/det_2010-*.parquet $1/famday_2010-*.parquet $1/run_2010-*.json; do
    [ -e $2/$(basename $p) ] || ln -s $(realpath $p) $2/$(basename $p)
  done
}
link_v2010 results/wa_gonogo/mf_core21 $R/mf_frozen_va
link_v2010 results/wa_pre/cap/mf_cap30 $R/mf_cap30_va
link_v2010 results/wa_pre/cap/mf_lodo $R/mf_lodo_va

run guard $C2 guard

# ---- GPU block (one lock for all matched-filter runs)
mkdir -p locks
if [ ! -e $R/stage_gpu_all.done ]; then
  hb "gpu_wait"; gpu_wait; hb "gpu lock taken"
  (
    flock 9
    run mf_frozen_v08 $PR --tag frozen_v08 --windows V2008 $MFARGS --weight-mode frozen
    run mf_cap30_v08  $PR --tag cap30_v08  --windows V2008 $MFARGS --weight-mode cap --cap 0.3
    run mf_lodo_v08   $PR --tag lodo_v08   --windows V2008 $MFARGS --weight-mode lodo
    run mf_frozen_va  $PR --tag frozen_va  --windows A2011 $MFARGS --weight-mode frozen
    run mf_cap30_va   $PR --tag cap30_va   --windows A2011 $MFARGS --weight-mode cap --cap 0.3
    run mf_lodo_va    $PR --tag lodo_va    --windows A2011 $MFARGS --weight-mode lodo
    run mf_revfrozen_v08 $PR --tag revfrozen_v08 --windows V2008 $REVARGS --weight-mode frozen
    run mf_revcap30_v08  $PR --tag revcap30_v08  --windows V2008 $REVARGS --weight-mode cap --cap 0.3
    run mf_revlodo_v08   $PR --tag revlodo_v08   --windows V2008 $REVARGS --weight-mode lodo
    run mf_revfrozen_va  $PR --tag revfrozen_va  --windows A2011 $REVARGS --weight-mode frozen
    run mf_revcap30_va   $PR --tag revcap30_va   --windows A2011 $REVARGS --weight-mode cap --cap 0.3
    run mf_revlodo_va    $PR --tag revlodo_va    --windows A2011 $REVARGS --weight-mode lodo
  ) 9>locks/wa_cap2.lock || exit 1
  touch $R/stage_gpu_all.done
fi

# ---- checks: completeness + the uncapped reruns equal the go/no-go matched filter
for t in frozen_v08 cap30_v08 lodo_v08 frozen_va cap30_va lodo_va revfrozen_v08 revcap30_v08 revlodo_v08 revfrozen_va revcap30_va revlodo_va; do
  run postguard_$t $C2 postguard $R/mf_$t
done
run equal_v08 $C2 equal $R/mf_frozen_v08 results/wa_gonogo/mf_core21_v2008 V2008
run equal_a11 $C2 equal $R/mf_frozen_va results/wa_gonogo/mf_core21 A2011

# ---- HF column (CPU, raw traces from the non-download roots)
for t in frozen_v08 cap30_v08 lodo_v08 frozen_va cap30_va lodo_va; do
  run hf_$t python scripts/pre_run.py --res $R --setup results/wa_gonogo/setup hf -- --mf $t --tpl build --workers 6 --all-templates
  run hfcheck_$t $C2 hfcheck $t
done

# ---- scoring (gonogo_score.py unchanged)
SC="python scripts/gonogo_score.py --tplmode all --stations results/wa_gonogo/setup/stations.csv"
for v in frozen cap30 lodo; do
  run score_${v}_v08 env WA_SETUP=$PWD/results/wa_gonogo/setup $SC --wins V2008 --mf ${v}_v08 --res $PWD/$R --od $PWD/$R/score_${v}_v08
  run score_${v}_va  env WA_SETUP=$PWD/results/wa_gonogo/setup $SC --wins V2010,A2011 --mf ${v}_va --res $PWD/$R --od $PWD/$R/score_${v}_va
done
run report $C2 report
hb "ALL DONE"
date -u +%FT%TZ > $R/CAP2_DONE
