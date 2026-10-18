#!/bin/bash
# Track A step 1(a) -- mergeB (adopted variant: within-site weighted-mean cc, site weight sqrt(sum w^2); no PGC cap, no
# coincidence gate) with 12-s templates on 2012, otherwise IDENTICAL to the E6 mergeB pass (scripts/e6_run_gpu1.sh):
# same templates/weights/moveout/thresholds/days/stations, 20 Hz 1-8 Hz chain frozen (Merlin F6).  Outputs -> lfe_extension/results.
cd /home/jovyan/tremorferometry/lfe_extension
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages CUDA_PATH=/opt/conda/targets/x86_64-linux OMP_NUM_THREADS=4
I=/home/jovyan/tremorferometry/lfe_catalog_repair/results
O=/home/jovyan/tremorferometry/lfe_extension/results
ST=PGC,VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
NS="PGC=CN.PGC,VGZ=CN.VGZ,NLLB=CN.NLLB,SNB=CN.SNB,LRIV=UW.LRIV,B001=PB.B001,B003=PB.B003,B004=PB.B004,B005=PB.B005,B006=PB.B006,B007=PB.B007,B009=PB.B009,B010=PB.B010,B011=PB.B011,B926=PB.B926"
SITES="PGC:PGC,B009,B010,B011;LRIV:LRIV,B005,B006,B007"
python scripts/e4_template_test.py --year 2012 --lengths 12 --n-random 600 --doy-from 228 --doy-to 279 \
  --stations $ST --netsta "$NS" --no-seisbench --families-csv $I/e5r_families_2012.csv --stations-csv $I/e5r_stations.csv \
  --gt-parquet $I/e5r_gt_bostock2012.parquet --offsets-csv $I/e5r_family_offsets_2012.csv \
  --extra-tpl $I/e5r_templates_2012.npz --skip-z --half all --moveout-in $I/e4_moveout_e5r_mv_build.csv \
  --budget --budget-thr 7,7.5,8,8.5,9,10 --catalog-out --amp-ref-csv $I/e5r_amp_ref_2012.csv \
  --weights $I/e5r_weights_snr_2012.csv --sites "$SITES" --site-merge B --out-dir $O --tag a1_mergeB_L12 \
  > logs/a1_mergeB_L12.log 2>&1
echo "exit=$? $(date -u +%FT%TZ)" >> logs/a1_mergeB_L12.log
