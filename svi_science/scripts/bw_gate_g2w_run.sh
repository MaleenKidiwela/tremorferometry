#!/usr/bin/env bash
# G2, A-2 original protocol: WITHHELD-STATION detections (forward detector, the target station left out of the network)
# as the reference events.  Needed because Bostock's 2003-2005 event times are ~0.3-0.6 s later and broader than in
# 2008-2012 at the never-changed PGC (bw_gate_g2_lagdiag.py), so Bostock-referenced stacks smear at every station.
# New legs: 2004-07 (CN4 minus target), 2005-09 (B1 minus target), 2008-05 (B3 minus target, cross-check);
# control legs: the SAME withheld network on the 2012 truth fit days (stacked on test blocks only).  <= 3 GPU passes.
set -u -o pipefail
S=/home/jovyan/tremorferometry/svi_science/scripts
O=/home/jovyan/tremorferometry/svi_science/results/backward_gates
export SVI_REF_DIR=$O/mf PYTHONUNBUFFERED=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=2
export LD_LIBRARY_PATH="$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib 2>/dev/null | tr '\n' ':')${LD_LIBRARY_PATH:-}"
rm -f $O/G2W_DONE
J=$O/logs/g2w_jobs.txt; : > $J
CB=/home/jovyan/tremorferometry/lfe_extension/cache20hz_backward
for x in PGC VGZ NLLB SNB; do
  echo "g2w2004_CN4no$x CN4no$x $O/tables/g2w_days_2004.json $CB" >> $J
  echo "g2w2005_B1no$x B1no$x $O/tables/g2w_days_2005.json $CB,/home/jovyan/tremorferometry/lfe_array_discovery/cache20hz" >> $J
  echo "g2w2008_B3no$x B3no$x $O/tables/g2w_days_2008.json $CB,/home/jovyan/tremorferometry/lfe_picker/cache20hz" >> $J
  for b in CN4 B1 B3; do echo "g2w2012_${b}no$x ${b}no$x $O/tables/g2w_days_2012.json /home/jovyan/tremorferometry/lfe_catalog_repair/cache20hz" >> $J; done
done
cat $J | xargs -P 3 -L 1 bash -c 'python '"$S"'/bw_gate_mf.py --backward -- --run $0 --kind truth --year 2012 --sset $1 --clean 0 \
   --days-file $2 --cache-dirs $3 --no-hfr > '"$O"'/logs/$0.log 2>&1; echo $? > '"$O"'/logs/$0.rc'
date > $O/G2W_DONE
