#!/usr/bin/env bash
# G3: era templates (build blocks, re-timed) -> detector (a) CN4 frozen vs (b) CN4 + POLARIS/LZB era templates on
# 2004-06-04..06-17 (quiet) + 2004-07-01..07-28 (ETS), same scrambled-control site list -> score on test blocks.
set -u -o pipefail
S=/home/jovyan/tremorferometry/svi_science/scripts
O=/home/jovyan/tremorferometry/svi_science/results/backward_gates
export SVI_REF_DIR=$O/mf PYTHONUNBUFFERED=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=2 BW_SITELIST_EXT=1
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
rm -f $O/G3_DONE $O/G3_FAILED
fail(){ echo "FAILED $1 $(date)"; date > $O/G3_FAILED; exit 1; }
while [ ! -f $O/G3_DATA_DONE ] && [ ! -f $O/G3_DATA_FAILED ]; do sleep 30; done
[ -f $O/G3_DATA_DONE ] || fail data
OMP_NUM_THREADS=4 python $S/bw_gate_g3_templates.py > $O/logs/g3_templates.log 2>&1 || fail templates
python - <<PY
import json, pandas as pd
d = [str(x.date()) for x in pd.date_range("2004-06-04", "2004-06-17")] + [str(x.date()) for x in pd.date_range("2004-07-01", "2004-07-28")]
json.dump(d, open("$O/tables/g3_days.json", "w"))
PY
CB=/home/jovyan/tremorferometry/lfe_extension/cache20hz_backward
J=$O/logs/g3_jobs.txt; : > $J
for i in 0 1; do
  echo "g3a_CN4 CN4 $i none" >> $J
  echo "g3b_CN4ERA CN4ERA $i $O/tables/g3_era_templates.npz" >> $J
done
cat $J | xargs -P 3 -L 1 bash -c '[ "$3" = none ] && E= || E=$3; BW_ERA_TPL=$E python '"$S"'/bw_gate_mf.py --backward -- --run $0 --kind truth --year 2012 \
   --sset $1 --clean 0 --days-file '"$O"'/tables/g3_days.json --cache-dirs '"$CB"' --shard $2 --nshard 2 > '"$O"'/logs/$0_s$2.log 2>&1; echo $? > '"$O"'/logs/$0_s$2.rc'
for r in $O/logs/g3a_CN4_s*.rc $O/logs/g3b_CN4ERA_s*.rc; do [ "$(cat $r)" = 0 ] || fail "mf $r"; done
python $S/bw_gate_g3_score.py > $O/logs/g3_score.log 2>&1 || fail score
date > $O/G3_DONE
