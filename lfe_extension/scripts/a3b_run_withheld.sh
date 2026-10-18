#!/bin/bash
# Track A step 3b -- the WITHHELD-STATION network: the frozen configuration (mergeB, 12-s templates, 2012 weights and
# moveout) with CN.VGZ, CN.SNB and PB.B926 EXCLUDED from the network sum, so those three stations can be used to confirm
# (or refute) the quiet-day detections without having contributed to them.  113 of the 116 templated families have a
# usable template at one of the three.
cd /home/jovyan/tremorferometry/lfe_extension
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages CUDA_PATH=/opt/conda/targets/x86_64-linux OMP_NUM_THREADS=3
I=/home/jovyan/tremorferometry/lfe_catalog_repair/results
O=/home/jovyan/tremorferometry/lfe_extension/results
ST=PGC,NLLB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011
NS="PGC=CN.PGC,VGZ=CN.VGZ,NLLB=CN.NLLB,SNB=CN.SNB,LRIV=UW.LRIV,B001=PB.B001,B003=PB.B003,B004=PB.B004,B005=PB.B005,B006=PB.B006,B007=PB.B007,B009=PB.B009,B010=PB.B010,B011=PB.B011,B926=PB.B926"
SITES="PGC:PGC,B009,B010,B011;LRIV:LRIV,B005,B006,B007"
python scripts/e4_template_test.py --year 2012 --lengths 12 --n-random 200 --doy-from 228 --doy-to 279 \
  --stations $ST --netsta "$NS" --no-seisbench --families-csv $I/e5r_families_2012.csv --stations-csv $I/e5r_stations.csv \
  --gt-parquet $I/e5r_gt_bostock2012.parquet --offsets-csv $I/e5r_family_offsets_2012.csv \
  --extra-tpl $I/e5r_templates_2012.npz --skip-z --half all --moveout-in $I/e4_moveout_e5r_mv_build.csv \
  --budget --budget-thr 8,8.25,8.5,8.75,9 --catalog-out --amp-ref-csv $I/e5r_amp_ref_2012.csv \
  --weights $I/e5r_weights_snr_2012.csv --sites "$SITES" --site-merge B --out-dir $O --tag a3b_noVSB \
  > logs/a3b_noVSB.log 2>&1
echo "x=$? $(date -u +%FT%TZ)" >> logs/a3b_noVSB.log
echo "[$(date -u +%FT%TZ)] 3b WITHHELD PASS DONE"
