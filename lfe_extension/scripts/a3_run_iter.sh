#!/bin/bash
# Track A step 3 -- apply the ITERATED templates: pass 1 the moveout on Bostock's 2012 BUILD half (method parity with
# the frozen point, which used the same table), pass 2 the full 2012 run, pass 3 the HELD-OUT 2010 year on the nine
# persistent stations (the only honest recall test, since the iterated templates were built from detections on both
# 2012 halves).  $1 = template tag (default a3_it1).
cd /home/jovyan/tremorferometry/lfe_extension
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages CUDA_PATH=/opt/conda/targets/x86_64-linux OMP_NUM_THREADS=3
V=${1:-a3_it1}
I=/home/jovyan/tremorferometry/lfe_catalog_repair/results
O=/home/jovyan/tremorferometry/lfe_extension/results
ST=PGC,VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
ST9=PGC,VGZ,SNB,NLLB,B926,LRIV,B009,B010,B011
NS="PGC=CN.PGC,VGZ=CN.VGZ,NLLB=CN.NLLB,SNB=CN.SNB,LRIV=UW.LRIV,B001=PB.B001,B003=PB.B003,B004=PB.B004,B005=PB.B005,B006=PB.B006,B007=PB.B007,B009=PB.B009,B010=PB.B010,B011=PB.B011,B926=PB.B926"
S15="PGC:PGC,B009,B010,B011;LRIV:LRIV,B005,B006,B007"
S9="PGC:PGC,B009,B010,B011"
B="--lengths 12 --n-random 200 --netsta $NS --no-seisbench --skip-z --families-csv $I/e5r_families_2012.csv \
  --stations-csv $I/e5r_stations.csv --extra-tpl $O/e5r_templates_${V}.npz --site-merge B --out-dir $O \
  --cache-dirs /home/jovyan/tremorferometry/lfe_extension/cache20hz"
echo "[$(date -u +%FT%TZ)] (1) moveout, 2012 build half"
python scripts/e4_template_test.py $B --year 2012 --doy-from 228 --doy-to 279 --stations $ST --sites "$S15" --half build \
  --gt-parquet $I/e5r_gt_bostock2012.parquet --offsets-csv $I/e5r_family_offsets_2012.csv --tag ${V}_mv > logs/${V}_mv.log 2>&1
echo "m=$? $(date -u +%FT%TZ)" >> logs/${V}_mv.log
echo "[$(date -u +%FT%TZ)] (2) full 2012"
python scripts/e4_template_test.py $B --year 2012 --doy-from 228 --doy-to 279 --stations $ST --sites "$S15" --half all \
  --gt-parquet $I/e5r_gt_bostock2012.parquet --offsets-csv $I/e5r_family_offsets_2012.csv \
  --weights $O/e5r_weights_snr_${V}.csv --moveout-in $O/e4_moveout_${V}_mv.csv --amp-ref-csv $I/e5r_amp_ref_2012.csv \
  --budget --budget-thr 7,7.5,8,8.25,8.5,8.75,9,10 --catalog-out --tag ${V}_2012 > logs/${V}_2012.log 2>&1
echo "f=$? $(date -u +%FT%TZ)" >> logs/${V}_2012.log
echo "[$(date -u +%FT%TZ)] (3) held-out 2010, nine persistent stations"
python scripts/e4_template_test.py $B --year 2010 --doy-from 222 --doy-to 256 --stations $ST9 --sites "$S9" --half all \
  --gt-parquet $O/a2_gt_bostock2010.parquet --offsets-csv $I/e3_bostock_family_offsets.csv \
  --weights $O/e5r_weights_snr_${V}.csv --moveout-in $O/e4_moveout_${V}_mv.csv \
  --budget --budget-thr 7,7.5,8,8.25,8.5,8.75,9,10 --catalog-out --tag ${V}_2010 > logs/${V}_2010.log 2>&1
echo "h=$? $(date -u +%FT%TZ)" >> logs/${V}_2010.log
echo "[$(date -u +%FT%TZ)] ITERATION RUNS DONE"
