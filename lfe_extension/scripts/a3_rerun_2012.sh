#!/bin/bash
# Track A step 3 -- re-run of the iterated-template 2012 pass after a CUDA OOM (the H100 is shared; three other
# processes held 88 GB).  Identical to a3_run_iter.sh pass 2 except --fam-chunk 8 (families are chunked, not
# throttled) and expandable_segments, so the peak allocation is a quarter of the default.
cd /home/jovyan/tremorferometry/lfe_extension
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages CUDA_PATH=/opt/conda/targets/x86_64-linux OMP_NUM_THREADS=3
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
V=${1:-a3_it1}
I=/home/jovyan/tremorferometry/lfe_catalog_repair/results
O=/home/jovyan/tremorferometry/lfe_extension/results
ST=PGC,VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
NS="PGC=CN.PGC,VGZ=CN.VGZ,NLLB=CN.NLLB,SNB=CN.SNB,LRIV=UW.LRIV,B001=PB.B001,B003=PB.B003,B004=PB.B004,B005=PB.B005,B006=PB.B006,B007=PB.B007,B009=PB.B009,B010=PB.B010,B011=PB.B011,B926=PB.B926"
python scripts/e4_template_test.py --year 2012 --lengths 12 --n-random 200 --doy-from 228 --doy-to 279 \
  --stations $ST --netsta "$NS" --no-seisbench --families-csv $I/e5r_families_2012.csv --stations-csv $I/e5r_stations.csv \
  --gt-parquet $I/e5r_gt_bostock2012.parquet --offsets-csv $I/e5r_family_offsets_2012.csv \
  --extra-tpl $O/e5r_templates_${V}.npz --skip-z --half all --moveout-in $O/e4_moveout_${V}_mv.csv \
  --weights $O/e5r_weights_snr_${V}.csv --amp-ref-csv $I/e5r_amp_ref_2012.csv \
  --sites "PGC:PGC,B009,B010,B011;LRIV:LRIV,B005,B006,B007" --site-merge B --out-dir $O \
  --budget --budget-thr 7,7.5,8,8.25,8.5,8.75,9,10 --catalog-out --fam-chunk 8 --tag ${V}_2012 \
  > logs/${V}_2012_retry.log 2>&1
echo "f=$? $(date -u +%FT%TZ)" >> logs/${V}_2012_retry.log
echo "[$(date -u +%FT%TZ)] 2012 RETRY DONE"
