#!/bin/bash
# E5 step 5 DIAGNOSTIC (not an operating point): the quiet-day excess over the rolled control -- site-coherent noise amplified by the
# co-located clusters (PGC/B009/B010/B011 within 0.3 km; LRIV/B005/B006/B007 within 0.3 km) or genuine moveout-consistent (inter-ETS LFE)
# signals?  Same weights, templates, moveout and code; ONE station per site; quiet window + first two SSE days (doy 228-244), 8 s.
cd /home/jovyan/tremorferometry/lfe_catalog_repair
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages CUDA_PATH=/opt/conda/targets/x86_64-linux OMP_NUM_THREADS=4
R=results
ST=PGC,VGZ,NLLB,SNB,LRIV,B001,B003,B004,B926
NS="PGC=CN.PGC,VGZ=CN.VGZ,NLLB=CN.NLLB,SNB=CN.SNB,LRIV=UW.LRIV,B001=PB.B001,B003=PB.B003,B004=PB.B004,B926=PB.B926"
python scripts/e4_template_test.py --year 2012 --lengths 8 --n-random 600 --doy-from 228 --doy-to 244 --stations $ST --netsta $NS --no-seisbench \
  --families-csv $R/e5r_families_2012.csv --stations-csv $R/e5r_stations.csv --gt-parquet $R/e5r_gt_bostock2012.parquet \
  --offsets-csv $R/e5r_family_offsets_2012.csv --extra-tpl $R/e5r_templates_2012.npz --skip-z --half all \
  --moveout-in $R/e4_moveout_e5r_mv_build.csv --budget --budget-thr 7,8,9,10 --weights $R/e5r_weights_snr_2012.csv --tag e5r_snr_onepersite > logs/e5r_mf_onepersite.log 2>&1
echo "[$(date -u +%FT%TZ)] onepersite DONE"
