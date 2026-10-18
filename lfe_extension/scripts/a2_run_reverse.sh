#!/bin/bash
# Track A step 2(c) -- the REVERSE frozen transfer: templates built from the 2010 events only (scripts/e5r_templates.py
# --year 2010, 9 persistent stations, results/e5r_templates_2010rev.npz, 716 usable pairs / 90 families, q-control cc to
# the SeisBench stacks 0.79-0.87) applied to 2012 at the SAME frozen thresholds.  Pass 1 = the 2010 build-half moveout
# (method parity with the forward direction, which used the 2012 build-half moveout); pass 2 = the 2012 run.
cd /home/jovyan/tremorferometry/lfe_extension
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages CUDA_PATH=/opt/conda/targets/x86_64-linux OMP_NUM_THREADS=4
I=/home/jovyan/tremorferometry/lfe_catalog_repair/results
O=/home/jovyan/tremorferometry/lfe_extension/results
ST=PGC,VGZ,SNB,NLLB,B926,LRIV,B009,B010,B011
NS="PGC=CN.PGC,VGZ=CN.VGZ,SNB=CN.SNB,NLLB=CN.NLLB,B926=PB.B926,LRIV=UW.LRIV,B009=PB.B009,B010=PB.B010,B011=PB.B011"
SITES="PGC:PGC,B009,B010,B011"
BASE="--lengths 12 --n-random 600 --stations $ST --netsta $NS --no-seisbench --skip-z \
  --families-csv $I/e5r_families_2012.csv --stations-csv $I/e5r_stations.csv \
  --extra-tpl $O/e5r_templates_2010rev.npz --sites $SITES --site-merge B --out-dir $O \
  --cache-dirs /home/jovyan/tremorferometry/lfe_extension/cache20hz"
echo "[$(date -u +%FT%TZ)] (c1) 2010 build-half moveout with the 2010 templates"
python scripts/e4_template_test.py $BASE --year 2010 --doy-from 222 --doy-to 256 --half build \
  --gt-parquet $O/a2_gt_bostock2010.parquet --offsets-csv $I/e3_bostock_family_offsets.csv \
  --tag a2_rev_mv > logs/a2_rev_mv.log 2>&1
echo "c1=$? $(date -u +%FT%TZ)" >> logs/a2_rev_mv.log
echo "[$(date -u +%FT%TZ)] (c2) 2010-built templates applied to 2012"
python scripts/e4_template_test.py $BASE --year 2012 --doy-from 228 --doy-to 279 --half all \
  --gt-parquet $I/e5r_gt_bostock2012.parquet --offsets-csv $I/e5r_family_offsets_2012.csv \
  --weights $O/e5r_weights_snr_2010rev.csv --moveout-in $O/e4_moveout_a2_rev_mv.csv \
  --amp-ref-csv $I/e5r_amp_ref_2012.csv --budget --budget-thr 7,7.5,8,8.25,8.5,8.75,9,10 --catalog-out \
  --tag a2_rev2012 > logs/a2_rev2012.log 2>&1
echo "c2=$? $(date -u +%FT%TZ)" >> logs/a2_rev2012.log
echo "[$(date -u +%FT%TZ)] STEP 2c DONE"
