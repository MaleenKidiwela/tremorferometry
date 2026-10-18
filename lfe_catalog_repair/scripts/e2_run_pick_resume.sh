#!/bin/bash
# E2 step 1 RESUME driver (after the 2026-09-10 21:03 pod crash) — only the unfinished runs.
# Same env, same 4-queue / 2-torch-thread layout and same output paths as e2_run_pick.sh.
# Already on disk (NOT redone): pgc2010 {cascadia,jcms}, vi2010 cascadia, vi2017 cascadia.
# e2_pick.py writes one parquet per SET at the end, so partially-finished sets are re-run whole.
cd /home/jovyan/tremorferometry
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages
L=lfe_catalog_repair/logs
R=lfe_catalog_repair/results
P="python lfe_catalog_repair/scripts/e2_pick.py --torch-threads 2"

q1() {
  $P --set vi2010 --model jcms --out $R/e2_picks_vi2010_jcms.parquet
  $P --set offgt --model cascadia --stations UW.WISH --out $R/e2_picks_offgt_WISH_cascadia.parquet
  $P --set offgt --model jcms     --stations UW.WISH --out $R/e2_picks_offgt_WISH_jcms.parquet
}
q2() {
  $P --set vi2017 --model jcms --out $R/e2_picks_vi2017_jcms.parquet
  $P --set offgt --model cascadia --stations UW.YACT --out $R/e2_picks_offgt_YACT_cascadia.parquet
  $P --set offgt --model jcms     --stations UW.YACT --out $R/e2_picks_offgt_YACT_jcms.parquet
}
q3() { $P --set offgt --model cascadia --stations PB.B943 --out $R/e2_picks_offgt_B943_cascadia.parquet; }
q4() { $P --set offgt --model jcms     --stations PB.B943 --out $R/e2_picks_offgt_B943_jcms.parquet; }

q1 > $L/e2_pick_resume_q1.log 2>&1 &
q2 > $L/e2_pick_resume_q2.log 2>&1 &
q3 > $L/e2_pick_resume_q3.log 2>&1 &
q4 > $L/e2_pick_resume_q4.log 2>&1 &
wait
echo "[e2 pick resume ALL DONE]"
