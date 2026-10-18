#!/bin/bash
# E2 step 1 driver — 4 parallel GPU queues, 2 torch threads each (8 CPU threads total).
cd /home/jovyan/tremorferometry
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages
L=lfe_catalog_repair/logs
R=lfe_catalog_repair/results
P="python lfe_catalog_repair/scripts/e2_pick.py --torch-threads 2"

q1() {
  for M in cascadia jcms; do
    $P --set pgc2010 --model $M --out $R/e2_picks_pgc2010_$M.parquet
    $P --set vi2010  --model $M --out $R/e2_picks_vi2010_$M.parquet
  done
  $P --set offgt --model cascadia --stations UW.WISH --out $R/e2_picks_offgt_WISH_cascadia.parquet
  $P --set offgt --model jcms     --stations UW.WISH --out $R/e2_picks_offgt_WISH_jcms.parquet
}
q2() {
  for M in cascadia jcms; do $P --set vi2017 --model $M --out $R/e2_picks_vi2017_$M.parquet; done
  $P --set offgt --model cascadia --stations UW.YACT --out $R/e2_picks_offgt_YACT_cascadia.parquet
  $P --set offgt --model jcms     --stations UW.YACT --out $R/e2_picks_offgt_YACT_jcms.parquet
}
q3() { $P --set offgt --model cascadia --stations PB.B943 --out $R/e2_picks_offgt_B943_cascadia.parquet; }
q4() { $P --set offgt --model jcms     --stations PB.B943 --out $R/e2_picks_offgt_B943_jcms.parquet; }

q1 > $L/e2_pick_q1.log 2>&1 &
q2 > $L/e2_pick_q2.log 2>&1 &
q3 > $L/e2_pick_q3.log 2>&1 &
q4 > $L/e2_pick_q4.log 2>&1 &
wait
echo "[e2 pick ALL DONE]"
