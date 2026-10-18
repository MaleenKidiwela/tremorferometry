#!/bin/bash
# Track A step 2 -- the FROZEN cross-year transfer test (Merlin plan-review step 2, make-or-break).
#   (a) a2_abl2012_7p : 2012, the FROZEN config (mergeB, 12-s templates, 2012 weights/moveout) but only the persistent
#                       stations PGC + its co-located boreholes B009/B010/B011 (one site), VGZ, SNB, NLLB, B926, LRIV
#                       -> the ablation PREDICTION for what 2010 should give.
#   (b) a2_xfer2010   : the SAME 2012-built templates, weights, moveout and MAD thresholds applied to 2010 doy 222-256.
#                       Ground truth = Bostock's 2010 events with the 90 ACCEPTED 2010 per-family offsets
#                       (e3_bostock_family_offsets.csv, accept == True), fam_id set to bkey (results/a2_gt_bostock2010.parquet).
#                       B009/B010 2010 were delivered by Track B and cached into lfe_extension/cache20hz; the other seven
#                       station-years are already in lfe_catalog_repair/cache20hz.
# Budget thresholds include the frozen T_ext = 8.25 and T_core = 8.75 exactly, so the site-rolled control is measured there.
cd /home/jovyan/tremorferometry/lfe_extension
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages CUDA_PATH=/opt/conda/targets/x86_64-linux OMP_NUM_THREADS=4
I=/home/jovyan/tremorferometry/lfe_catalog_repair/results
O=/home/jovyan/tremorferometry/lfe_extension/results
ST=PGC,VGZ,SNB,NLLB,B926,LRIV,B009,B010,B011
NS="PGC=CN.PGC,VGZ=CN.VGZ,SNB=CN.SNB,NLLB=CN.NLLB,B926=PB.B926,LRIV=UW.LRIV,B009=PB.B009,B010=PB.B010,B011=PB.B011"
SITES="PGC:PGC,B009,B010,B011"
THR=7,7.5,8,8.25,8.5,8.75,9,10
COMMON="--lengths 12 --n-random 600 --stations $ST --netsta $NS --no-seisbench --skip-z --half all \
  --families-csv $I/e5r_families_2012.csv --stations-csv $I/e5r_stations.csv \
  --extra-tpl $I/e5r_templates_2012.npz --moveout-in $I/e4_moveout_e5r_mv_build.csv \
  --weights $I/e5r_weights_snr_2012.csv --sites $SITES --site-merge B \
  --budget --budget-thr $THR --catalog-out --out-dir $O \
  --cache-dirs /home/jovyan/tremorferometry/lfe_extension/cache20hz"
echo "[$(date -u +%FT%TZ)] (a) 2012 ablation, 7 persistent sites"
python scripts/e4_template_test.py $COMMON --year 2012 --doy-from 228 --doy-to 279 \
  --gt-parquet $I/e5r_gt_bostock2012.parquet --offsets-csv $I/e5r_family_offsets_2012.csv \
  --amp-ref-csv $I/e5r_amp_ref_2012.csv --tag a2_abl2012_7p > logs/a2_abl2012_7p.log 2>&1
echo "a=$? $(date -u +%FT%TZ)" >> logs/a2_abl2012_7p.log
echo "[$(date -u +%FT%TZ)] (b) frozen transfer to 2010"
python scripts/e4_template_test.py $COMMON --year 2010 --doy-from 222 --doy-to 256 \
  --gt-parquet $O/a2_gt_bostock2010.parquet --offsets-csv $I/e3_bostock_family_offsets.csv \
  --tag a2_xfer2010 > logs/a2_xfer2010.log 2>&1
echo "b=$? $(date -u +%FT%TZ)" >> logs/a2_xfer2010.log
echo "[$(date -u +%FT%TZ)] STEP 2 a+b DONE"
