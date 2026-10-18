#!/usr/bin/env bash
# P4-2 GPU queue: forward detector (svi_ref_mf.py via bw_gate_mf.py, unchanged) with the reduced station sets CN4/B1/B2/B3
# on every E0 PNSN-quiet day (2014-01..2017-02, 509 days), forward R2 exclusions applied (--sset eff + per-set exclusion
# file, see P4.md), then the extra S15 fixed set without exclusions.  ONE GPU process at a time (sequential).
# Markers p4/P4_2_<SET>_DONE|FAILED, then p4/P4_2_DONE after scoring (bw_p4_sets.py).  Nothing deleted.
set -u -o pipefail
source /home/jovyan/tremorferometry/svi_science/scripts/bw_prod_env.sh
P=$O/p4
echo $$ > $P/Q2_PID
export SVI_REF_DIR=$P/mf
D=$P/tables/p4_e0_quiet_days.json
PREP=/home/jovyan/tremorferometry/lfe_extension/results/refine/prep
guard(){ python -c "import sys; sys.path.insert(0,'$LFE/scripts'); import svi_ref_common as R; R.guard('mf_$1','mf_$1_s*')" >> $P/logs/q2_guard.log 2>&1; }
for SET in CN4 B1 B2 B3 S15fix; do
  [ -f $P/P4_2_${SET}_DONE ] && continue
  if [ $SET = S15fix ]; then ARGS="--sset S15"; else ARGS="--sset eff --excl-file $P/tables/p4_excl_$SET.csv"; fi
  echo "[q2] $SET start $(date -u +%FT%TZ)" >> $P/logs/q2.log
  python $S/bw_gate_mf.py -- --run p4e0_$SET --kind prod $ARGS --clean 0 --days-file $D --prep-dir $PREP \
     --shard 0 --nshard 1 > $P/logs/p4e0_$SET.log 2>&1 && guard p4e0_$SET \
     && date -u +%FT%TZ > $P/P4_2_${SET}_DONE || date -u +%FT%TZ > $P/P4_2_${SET}_FAILED
  echo "[q2] $SET end $(date -u +%FT%TZ)" >> $P/logs/q2.log
done
python $S/bw_p4_sets.py > $P/logs/p4_sets.log 2>&1 && date -u +%FT%TZ > $P/P4_2_DONE || date -u +%FT%TZ > $P/P4_2_FAILED
