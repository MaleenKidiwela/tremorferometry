#!/bin/bash
# Track A step 3 -- SELECTION passes for template iteration (Merlin plan-review step 3).
#   a3_sel_gated     : the frozen config (mergeB, 12 s, 15 stations) with the two-site coincidence gate
#                      (>= 2 sites each >= 2.5 MAD of their OWN site series) and the catalog written from T_core = 8.75.
#   a3_sel_loso_PGC  : the same, with the PGC SITE (PGC, B009, B010, B011) EXCLUDED from the network sum, so the events
#                      used to re-stack the PGC-site templates were never selected by those stations (leave-one-site-out;
#                      the PGC site carries > 40 % of most families' weight).
cd /home/jovyan/tremorferometry/lfe_extension
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages CUDA_PATH=/opt/conda/targets/x86_64-linux OMP_NUM_THREADS=4
I=/home/jovyan/tremorferometry/lfe_catalog_repair/results
O=/home/jovyan/tremorferometry/lfe_extension/results
ST=PGC,VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
ST_noPGC=VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B926
NS="PGC=CN.PGC,VGZ=CN.VGZ,NLLB=CN.NLLB,SNB=CN.SNB,LRIV=UW.LRIV,B001=PB.B001,B003=PB.B003,B004=PB.B004,B005=PB.B005,B006=PB.B006,B007=PB.B007,B009=PB.B009,B010=PB.B010,B011=PB.B011,B926=PB.B926"
SITES="PGC:PGC,B009,B010,B011;LRIV:LRIV,B005,B006,B007"
BASE="--year 2012 --lengths 12 --n-random 200 --doy-from 228 --doy-to 279 --netsta $NS --no-seisbench \
  --families-csv $I/e5r_families_2012.csv --stations-csv $I/e5r_stations.csv --gt-parquet $I/e5r_gt_bostock2012.parquet \
  --offsets-csv $I/e5r_family_offsets_2012.csv --extra-tpl $I/e5r_templates_2012.npz --skip-z --half all \
  --moveout-in $I/e4_moveout_e5r_mv_build.csv --weights $I/e5r_weights_snr_2012.csv --sites $SITES --site-merge B \
  --budget --budget-thr 8.75 --catalog-out --amp-ref-csv $I/e5r_amp_ref_2012.csv --out-dir $O --coinc-gate 2,2.5"
echo "[$(date -u +%FT%TZ)] selection pass 1 (gated, all 15 stations)"
python scripts/e4_template_test.py $BASE --stations $ST --tag a3_sel_gated > logs/a3_sel_gated.log 2>&1
echo "s1=$? $(date -u +%FT%TZ)" >> logs/a3_sel_gated.log
echo "[$(date -u +%FT%TZ)] selection pass 2 (gated, PGC site left out of the sum)"
python scripts/e4_template_test.py $BASE --stations $ST_noPGC --tag a3_sel_loso_PGC > logs/a3_sel_loso_PGC.log 2>&1
echo "s2=$? $(date -u +%FT%TZ)" >> logs/a3_sel_loso_PGC.log
echo "[$(date -u +%FT%TZ)] SELECTION DONE"
