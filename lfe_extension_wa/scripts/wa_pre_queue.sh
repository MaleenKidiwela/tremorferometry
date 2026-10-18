#!/bin/bash
# WA pre-production queue (2026-10-02). Markers: results/wa_pre/WA_PRE_{PID,HEARTBEAT,DONE,FAILED}, stage_<name>.done.
# Read-only on existing results; outputs under results/wa_pre/ only; never reads traces_dl*/traces_sealed; no downloads,
# no deletions, no git. GPU stages hold results/wa_pre/gpu.lock and wait while any other project GPU lock is held.
set -u
cd /home/jovyan/tremorferometry/lfe_extension_wa
R=results/wa_pre
LOG=logs/wa_pre_queue.log
echo $$ > $R/WA_PRE_PID
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')${LD_LIBRARY_PATH:-}
export PYTHONWARNINGS=ignore
hb() { echo "$(date -u +%FT%TZ) $1" > $R/WA_PRE_HEARTBEAT; echo "$(date -u +%FT%TZ) [queue] $1" >> $LOG; }
fail() { echo "$(date -u +%FT%TZ) stage $1 failed (rc $2)" > $R/WA_PRE_FAILED; hb "FAILED $1 rc $2"; exit 1; }
gpu_wait() {
  while true; do
    held=0
    for f in /home/jovyan/tremorferometry/*/locks/*.lock; do
      [ -e "$f" ] || continue
      flock -n "$f" true 2>/dev/null || { held=1; echo "$(date -u +%FT%TZ) waiting for $f" >> $LOG; }
    done
    [ $held -eq 0 ] && break
    sleep 120
  done
}
run() {  # run <stage> <gpu 0|1> <cmd...>
  st=$1; g=$2; shift 2
  [ -e $R/stage_$st.done ] && { hb "skip $st (done)"; return 0; }
  hb "start $st"
  if [ "$g" = "1" ]; then gpu_wait; flock $R/gpu.lock "$@" >> $LOG 2>&1; rc=$?; else "$@" >> $LOG 2>&1; rc=$?; fi
  [ $rc -eq 0 ] || fail $st $rc
  touch $R/stage_$st.done; hb "done $st"
}
PR="python scripts/pre_run.py"
GOS=results/wa_gonogo/setup
MFARGS="--tpl build --n-control 2 --drop-dead 0.05 --all-templates"

# ---- GPU matched-filter runs (cheap CPU; may overlap the W1 audit)
run cap_mf_cap30 1 $PR --res $R/cap --setup $GOS mfcap -- --tag cap30 --windows V2010 $MFARGS --weight-mode cap --cap 0.3
run cap_mf_lodo  1 $PR --res $R/cap --setup $GOS mfcap -- --tag lodo  --windows V2010 $MFARGS --weight-mode lodo
D60=$(python -c "import pandas as pd; print(','.join(pd.read_csv('$R/ctl/days60.csv').date))")
run ctl_mf_n60   1 $PR --res $R/ctl/res --setup $GOS --inject-window D60=$R/ctl/days60_kind.csv mf -- --tag n60 --windows D60 $MFARGS
run ctl_mf_rev60 1 $PR --res $R/ctl/res --setup $GOS --inject-window D60=$R/ctl/days60_kind.csv mf -- --tag rev60 --windows D60 --tpl buildrev --no-control --drop-dead 0.05 --all-templates

# ---- wait for the W1 audit (CPU) before the CPU-heavy stages
while [ ! -e $R/W1_EXIT ]; do hb "waiting for W1 audit"; sleep 60; done
[ "$(cat $R/W1_EXIT)" = "0" ] || fail w1_complete "$(cat $R/W1_EXIT)"
touch $R/stage_w1_complete.done

# ---- step 3: HF column, scoring, report
run cap_hf_cap30 0 $PR --res $R/cap --setup $GOS hf -- --mf cap30 --tpl build --workers 8 --all-templates
run cap_hf_lodo  0 $PR --res $R/cap --setup $GOS hf -- --mf lodo  --tpl build --workers 8 --all-templates
SC="python scripts/gonogo_score.py --tplmode all --stations $GOS/stations.csv --wins V2010"
run cap_score_frozen 0 env WA_SETUP=$PWD/$GOS $SC --mf core21 --res $PWD/results/wa_gonogo --od $PWD/$R/cap/score_frozen
run cap_score_cap30  0 env WA_SETUP=$PWD/$GOS $SC --mf cap30  --res $PWD/$R/cap --od $PWD/$R/cap/score_cap30
run cap_score_lodo   0 env WA_SETUP=$PWD/$GOS $SC --mf lodo   --res $PWD/$R/cap --od $PWD/$R/cap/score_lodo
run cap_report 0 python scripts/pre_cap_report.py

# ---- step 4: per-site z, contamination, report
CC="python scripts/pre_ctlclean.py"
run ctl_sitez_n60    0 $CC sitez --mf $R/ctl/res/mf_n60 --setup $GOS --tpl build --out $R/ctl/sitez_n60 --workers 12
run ctl_sitez_core21 0 $CC sitez --mf results/wa_gonogo/mf_core21 --setup $GOS --tpl build --out $R/ctl/sitez_core21 --workers 12
run ctl_sitez_mt1    0 $CC sitez --mf results/mf_mt1_all --mf2 results/mf_step1c34 --setup results/merlin_tests/t1_setup --tpl build --out $R/ctl/sitez_mt1all --workers 12
run ctl_contam_n60    0 $CC contam --mf $R/ctl/res/mf_n60 --setup $GOS --sitez $R/ctl/sitez_n60 --out $R/ctl/contam_n60.parquet
run ctl_contam_core21 0 $CC contam --mf results/wa_gonogo/mf_core21 --setup $GOS --sitez $R/ctl/sitez_core21 --out $R/ctl/contam_core21.parquet
run ctl_contam_mt1    0 $CC contam --mf results/mf_mt1_all --setup results/merlin_tests/t1_setup --sitez $R/ctl/sitez_mt1all --out $R/ctl/contam_mt1all.parquet
run ctl_report 0 python scripts/pre_ctl_report.py

# ---- W2 station health
W2="python scripts/pre_w2.py"
run w2_cache 0 $W2 cache --workers 16
run w2_mf    1 $PR --res $R/w2/res --setup $R/w2/setup --cache-dirs $R/w2/cache20hz --inject-window W2=$R/w2/mf_days.csv mf -- --tag w2 --windows W2 --tpl build,step4b_a --all-templates --drop-dead 0.05 --no-control
run w2_lo    0 $W2 lo --workers 12
run w2_stack 0 $W2 stack --workers 16
run w2_raw   0 $W2 raw --workers 24
run w2_summary 0 $W2 summary --workers 16
run w0_verify 0 python scripts/pre_w01.py verify
run report 0 python scripts/pre_report.py
hb "ALL DONE"
date -u +%FT%TZ > $R/WA_PRE_DONE
