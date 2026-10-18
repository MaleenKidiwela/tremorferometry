#!/usr/bin/env bash
# Step 4b queue: ADD STATIONS for the NW + Olympic families, bridge templating, re-run the trial.
#   nohup bash scripts/step4b_queue.sh > logs/step4b_queue.out 2>&1 &
#   results/STEP4B_PID         queue pid        results/STEP4B_HEARTBEAT  "<utc> <stage>" every 60 s
#   results/STEP4B_DONE        all stages done  results/STEP4B_FAILED     stage + log tail on failure
#   results/step4b/stage_<name>.done   per-stage guard (delete one to re-run that stage); resumable
#   stop: kill $(cat results/STEP4B_PID)   (exact pid; the heartbeat child is killed by the trap)
# Runs:  V2010 before / after (validation, known answer, test half)
#        T2012 / T2013 / T2014: before (42-station templates) / p1 (+ array-era bridge templates and earlier
#        episodes) / p2 (+ templates of this episode's pass-1 detections of the separating families)
set -u
S=/home/jovyan/tremorferometry/lfe_extension_wa
cd "$S" || exit 1
R=$S/results; L=$S/logs; G=$R/step4b
mkdir -p "$G" "$L"
exec 9>"$R/.step4b.lock"
if ! flock -n 9; then echo "another step4b_queue holds the lock"; exit 1; fi
echo $$ > "$R/STEP4B_PID"
rm -f "$R/STEP4B_DONE" "$R/STEP4B_FAILED"
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}
export PYTHONUNBUFFERED=1 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
export WA_SETUP=$G/setup WA_EXTRA_NODES=$G/setup/extra_nodes.csv STEP4B_PREFIX=step4b_
( while true; do echo "$(date -u +%FT%TZ) $(cat $G/.stage 2>/dev/null)" > "$R/STEP4B_HEARTBEAT"; sleep 60; done ) &
HB=$!
trap 'kill $HB 2>/dev/null' EXIT

run() {  # run <name> <cmd...>
  local name=$1; shift
  if [ -f "$G/stage_$name.done" ]; then echo "[skip] $name"; return 0; fi
  echo "$name" > "$G/.stage"
  echo "$(date -u +%FT%TZ) [stage] $name: $*" | tee -a "$L/step4b_queue.log"
  if "$@" >> "$L/step4b_$name.log" 2>&1; then
    date -u +%FT%TZ > "$G/stage_$name.done"
    echo "$(date -u +%FT%TZ) [ok] $name" | tee -a "$L/step4b_queue.log"
  else
    { echo "stage $name failed $(date -u +%FT%TZ)"; tail -40 "$L/step4b_$name.log"; } > "$R/STEP4B_FAILED"
    echo "$(date -u +%FT%TZ) [FAILED] $name" | tee -a "$L/step4b_queue.log"
    exit 1
  fi
}

df -h /home/jovyan | tail -1 | tee -a "$L/step4b_queue.log"
[ -f "$G/setup/families.csv" ] || run setup python scripts/step4b_setup.py
run fetch     python scripts/wa_fetch.py --workers 8
run fetch2    python scripts/wa_fetch.py --workers 4          # retry pass for transient errors
run cache     python scripts/wa_cache.py --workers 12
run events_a  python scripts/step4b_events.py --mode array
run tpl_a     python scripts/wa_templates.py --tag step4b_a --events $G/events_array.parquet \
                --only-new-nodes $G/nodes_array.csv --weights-ref $R/templates/tpl_qc_build.csv --workers 12

B=build,t2012new
A=$B,step4b_a
# ---------------- V2010 validation: before (42-station templates) and after (+ added stations / epochs)
run mf_v2010b    python scripts/wa_mf.py --tag step4b_v2010b --windows V2010 --tpl $B --n-control 2
run hf_v2010b    python scripts/wa_hf.py --mf step4b_v2010b --tpl $B --workers 16 --cuts $R/hf_cuts.csv
run score_v2010b python scripts/wa_score.py --mode validation --mf step4b_v2010b --tpl $B --cuts $G/hf_cuts_v2010b.csv
run mf_v2010     python scripts/wa_mf.py --tag step4b_v2010 --windows V2010 --tpl $A --n-control 2
run hf_v2010     python scripts/wa_hf.py --mf step4b_v2010 --tpl $A --workers 16 --cuts $R/hf_cuts.csv
run score_v2010  python scripts/wa_score.py --mode validation --mf step4b_v2010 --tpl $A --cuts $G/hf_cuts_v2010.csv
# trial cuts = the frozen trial cuts (V2008 + V2010) + cuts of HF nodes new in the step-4b V2010 validation
run cuts python -c "
import pandas as pd
o = pd.read_csv('$R/hf_cuts.csv'); n = pd.read_csv('$G/hf_cuts_v2010.csv')
pd.concat([o.assign(src='trial V2008+V2010'), n[~n.hf_node.isin(o.hf_node)].assign(src='step4b V2010')]).to_csv('$G/hf_cuts_merged.csv', index=False)"
CUT=$G/hf_cuts_merged.csv

# ---------------- trial episodes, bridge rule iterated episode by episode
PREV=""
for E in T2012 T2013 T2014; do
  e=$(echo $E | tr 'A-Z' 'a-z')
  run mf_${e}b     python scripts/wa_mf.py --tag step4b_${e}b --windows $E --tpl $B --n-control 2
  run hf_${e}b     python scripts/wa_hf.py --mf step4b_${e}b --tpl $B --workers 16 --cuts $CUT
  run score_${e}b  python scripts/wa_score.py --mode trial --mf step4b_${e}b --tpl $B --cuts $CUT
  P1=$A$PREV
  run mf_${e}p1    python scripts/wa_mf.py --tag step4b_${e}p1 --windows $E --tpl $P1 --n-control 2
  run hf_${e}p1    python scripts/wa_hf.py --mf step4b_${e}p1 --tpl $P1 --workers 16 --cuts $CUT
  run score_${e}p1 python scripts/wa_score.py --mode trial --mf step4b_${e}p1 --tpl $P1 --cuts $CUT
  run events_${e}  python scripts/step4b_events.py --mode episode --episode $E
  run tpl_${e}     python scripts/wa_templates.py --tag step4b_${e} --events $G/events_$E.parquet \
                     --only-new-nodes $G/nodes_$E.csv --weights-ref $R/templates/tpl_qc_build.csv --workers 12
  P2=$P1,step4b_${e}
  run mf_${e}p2    python scripts/wa_mf.py --tag step4b_${e}p2 --windows $E --tpl $P2 --n-control 2
  run hf_${e}p2    python scripts/wa_hf.py --mf step4b_${e}p2 --tpl $P2 --workers 16 --cuts $CUT
  run score_${e}p2 python scripts/wa_score.py --mode trial --mf step4b_${e}p2 --tpl $P2 --cuts $CUT
  PREV=$PREV,step4b_${e}
done
run report python scripts/step4b_report.py
echo "done" > "$G/.stage"
date -u +%FT%TZ > "$R/STEP4B_DONE"
echo "$(date -u +%FT%TZ) [DONE]" | tee -a "$L/step4b_queue.log"
