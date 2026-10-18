#!/bin/bash
# E5 step 5 STEP D -- NETWORK MATCHED FILTER on the 2012 data, from-scratch templates only (e4_template_test.py behind flags).
# pass 1: per-(family, station, L) residual moveout = median best lag at the BUILD-half events (+-3 s scan) -- build half only.
# pass 2: all 52 days, --half all (events of both halves scored and labelled), SNR weights from the BUILD-half stack S-SNR,
#         --budget (detections / rolled control / +-2 s matches per family-day-threshold at 7,8,9,10 MAD) and --catalog-out
#         (every declustered detection >= 7 MAD with network MAD, cc, relative amplitude).  8 and 12 s.  Equal weights as reference.
cd /home/jovyan/tremorferometry/lfe_catalog_repair
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages CUDA_PATH=/opt/conda/targets/x86_64-linux OMP_NUM_THREADS=4
set -e
R=results
ST=PGC,VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
NS="PGC=CN.PGC,VGZ=CN.VGZ,NLLB=CN.NLLB,SNB=CN.SNB,LRIV=UW.LRIV,B001=PB.B001,B003=PB.B003,B004=PB.B004,B005=PB.B005,B006=PB.B006,B007=PB.B007,B009=PB.B009,B010=PB.B010,B011=PB.B011,B926=PB.B926"
COMMON="--year 2012 --lengths 8,12 --n-random 600 --doy-from 228 --doy-to 279 --stations $ST --netsta $NS --no-seisbench \
  --families-csv $R/e5r_families_2012.csv --stations-csv $R/e5r_stations.csv --gt-parquet $R/e5r_gt_bostock2012.parquet \
  --offsets-csv $R/e5r_family_offsets_2012.csv --extra-tpl $R/e5r_templates_2012.npz --skip-z"
T="python scripts/e4_template_test.py $COMMON"
echo "[$(date -u +%FT%TZ)] pass 1 (build-half moveout)"
$T --half build --tag e5r_mv_build > logs/e5r_mf_mv_build.log 2>&1
MV=$R/e4_moveout_e5r_mv_build.csv
echo "[$(date -u +%FT%TZ)] pass 2 SNR-weighted"
$T --half all --moveout-in $MV --budget --budget-thr 7,8,9,10 --catalog-out --amp-ref-csv $R/e5r_amp_ref_2012.csv \
   --weights $R/e5r_weights_snr_2012.csv --tag e5r_snr > logs/e5r_mf_snr.log 2>&1
echo "[$(date -u +%FT%TZ)] pass 2 equal-weight reference"
$T --half all --moveout-in $MV --budget --budget-thr 7,8,9,10 --catalog-out --amp-ref-csv $R/e5r_amp_ref_2012.csv \
   --tag e5r_eq > logs/e5r_mf_eq.log 2>&1
echo "[$(date -u +%FT%TZ)] STEP D DONE"
