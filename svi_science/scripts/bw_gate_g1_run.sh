#!/usr/bin/env bash
# G1 (no download): the forward detector with the CN4 / B1 / B2 station sets on the forward truth days on disk.
# <= 3 GPU passes at once (A-2 OOM lesson).  Markers OUT/G1_RUN_{DONE,FAILED}.  Nothing deleted, nothing downloaded.
set -u -o pipefail
S=/home/jovyan/tremorferometry/svi_science/scripts
O=/home/jovyan/tremorferometry/svi_science/results/backward_gates
D=/home/jovyan/tremorferometry/lfe_extension/results/refine/days
PREP=/home/jovyan/tremorferometry/lfe_extension/results/refine/prep
export SVI_REF_DIR=$O/mf PYTHONUNBUFFERED=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
mkdir -p $O/logs; rm -f $O/G1_RUN_DONE $O/G1_RUN_FAILED
J=$O/logs/g1_jobs.txt; : > $J
for y in 2010 2011 2012; do echo "g1t${y}_CN4 $y CN4" >> $J; done
for y in 2011 2012; do for s in B2 B1; do echo "g1t${y}_${s} $y $s" >> $J; done; done
cat $J | xargs -P 3 -L 1 bash -c 'python '"$S"'/bw_gate_mf.py -- --run $0 --kind truth --year $1 --sset $2 --clean 0 \
   --days-file '"$D"'/truth_$1.json --prep-dir '"$PREP"' > '"$O"'/logs/$0.log 2>&1; echo $? > '"$O"'/logs/$0.rc'
bad=0
while read r y s; do [ "$(cat $O/logs/$r.rc 2>/dev/null)" = 0 ] || { bad=1; echo "FAILED $r"; }
  python -c "import sys; sys.path.insert(0,'/home/jovyan/tremorferometry/lfe_extension/scripts'); import svi_ref_common as R; R.guard('mf_$r','mf_${r}_s*')" || bad=1
done < $J
[ $bad = 0 ] && date > $O/G1_RUN_DONE || date > $O/G1_RUN_FAILED
