#!/usr/bin/env bash
# WA permanent-station TRIAL extension queue.  Guarded, resumable stages; markers in results/.
#   nohup bash scripts/wa_trial_queue.sh > logs/queue.out 2>&1 &
#   results/WA_TRIAL_PID        queue pid          results/WA_TRIAL_HEARTBEAT  "<utc> <stage>" every 60 s
#   results/WA_TRIAL_DONE       all stages done    results/WA_TRIAL_FAILED     stage + log tail on failure
#   results/stage_<name>.done   per-stage guard (delete one to re-run that stage)
#   stop: kill $(cat results/WA_TRIAL_PID)   (exact pid; the heartbeat child is killed by the trap)
set -u
S=/home/jovyan/tremorferometry/lfe_extension_wa
cd "$S" || exit 1
R=$S/results; L=$S/logs
exec 9>"$R/.wa_trial.lock"
if ! flock -n 9; then echo "another wa_trial_queue holds the lock"; exit 1; fi
echo $$ > "$R/WA_TRIAL_PID"
rm -f "$R/WA_TRIAL_DONE" "$R/WA_TRIAL_FAILED"
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
STAGE=start
( while true; do echo "$(date -u +%FT%TZ) $(cat $R/.stage 2>/dev/null)" > "$R/WA_TRIAL_HEARTBEAT"; sleep 60; done ) &
HB=$!
trap 'kill $HB 2>/dev/null' EXIT

run() {  # run <name> <cmd...>
  local name=$1; shift
  if [ -f "$R/stage_$name.done" ]; then echo "[skip] $name"; return 0; fi
  echo "$name" > "$R/.stage"
  echo "$(date -u +%FT%TZ) [stage] $name: $*" | tee -a "$L/queue.log"
  if "$@" >> "$L/stage_$name.log" 2>&1; then
    date -u +%FT%TZ > "$R/stage_$name.done"
    echo "$(date -u +%FT%TZ) [ok] $name" | tee -a "$L/queue.log"
  else
    { echo "stage $name failed $(date -u +%FT%TZ)"; tail -40 "$L/stage_$name.log"; } > "$R/WA_TRIAL_FAILED"
    echo "$(date -u +%FT%TZ) [FAILED] $name" | tee -a "$L/queue.log"
    exit 1
  fi
}

# the pre-launch fetch (started during development) must finish before the queue's own fetch pass
if [ -f "$L/fetch_pre.pid" ]; then
  P=$(cat "$L/fetch_pre.pid")
  while [ -d /proc/$P ] && grep -q wa_fetch /proc/$P/cmdline 2>/dev/null; do echo "waitfetch" > "$R/.stage"; sleep 30; done
fi

[ -f "$R/setup/families.csv" ] || run setup python scripts/wa_setup.py
run fetch      python scripts/wa_fetch.py --workers 8
run cache      python scripts/wa_cache.py --workers 16
run tpl_build  python scripts/wa_templates.py --tag build --workers 12
run mf_val     python scripts/wa_mf.py --tag val --windows V2008,V2010 --tpl build
run hf_val     python scripts/wa_hf.py --mf val --tpl build --workers 16
run score_val  python scripts/wa_score.py --mode validation --mf val --tpl build --cuts $R/hf_cuts.csv
run mf_t1      python scripts/wa_mf.py --tag t2012p1 --windows T2012 --tpl build
run hf_t1      python scripts/wa_hf.py --mf t2012p1 --tpl build --workers 16 --cuts $R/hf_cuts.csv
run score_t1   python scripts/wa_score.py --mode trial --mf t2012p1 --tpl build --cuts $R/hf_cuts.csv
run newnodes   python scripts/wa_newnodes.py --mf t2012p1 --tpl build
NNEW=$(python -c "import pandas as pd; print(len(pd.read_csv('$R/trial/new_nodes.csv')))")
if [ "$NNEW" -gt 0 ]; then
  run tpl_new  python scripts/wa_templates.py --tag t2012new --events $R/trial/pass1_events.parquet \
                 --only-new-nodes $R/trial/new_nodes.csv --weights-ref $R/templates/tpl_qc_build.csv --workers 12
  run mf_t2    python scripts/wa_mf.py --tag t2012p2 --windows T2012 --tpl build,t2012new
  run hf_t2    python scripts/wa_hf.py --mf t2012p2 --tpl build,t2012new --workers 16 --cuts $R/hf_cuts.csv
  run score_t2 python scripts/wa_score.py --mode trial --mf t2012p2 --tpl build,t2012new --cuts $R/hf_cuts.csv
fi
run report     python scripts/wa_report.py
echo "done" > "$R/.stage"
date -u +%FT%TZ > "$R/WA_TRIAL_DONE"
echo "$(date -u +%FT%TZ) [DONE]" | tee -a "$L/queue.log"
