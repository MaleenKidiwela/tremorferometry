#!/bin/bash
# Track D step D1(a) -- EPOCH / STATION-SET ABLATION (Merlin Step 6, component 1).
# Re-runs the FROZEN configuration (mergeB, 12-s templates, 2012 templates/weights/moveout, one-to-one
# scoring) on the referee years with each EXTENSION EPOCH's station set, so that per-family recall can be
# predicted for every epoch of the 2014-2026 production run.
#
# Extension epochs actually realised by Track C (results/trackC/templates_by_epoch.csv + station_epochs.csv):
#   E1 2014-01-01..2019-05-29  15 stations  (= the frozen 2012 network; already run as a1_mergeB_L12 / a2)
#   E2 2019-05-30..2021-03-27  14 stations  (LRIV DROPPED by the epoch policy)
#   E3 2021-03-28..2026-09-15  13 stations  (LRIV dropped + SNB dead 2021-03-28)
# CN.PFB (ends 2018-12-20) is NOT in the network at all -- it has no 2012 template -- so the "post-PFB"
# boundary of the brief has no effect on our detector and needs no ablation.
#
# 2012 = the referee year (test half scored); 2010 = the held-out year (9-station base network, the only
# stations cached for 2010), ablated the same way so that the EPOCH EFFECT (a ratio) can be validated
# out of sample.
cd /home/jovyan/tremorferometry/lfe_extension
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages CUDA_PATH=/opt/conda/targets/x86_64-linux OMP_NUM_THREADS=4
I=/home/jovyan/tremorferometry/lfe_catalog_repair/results
O=/home/jovyan/tremorferometry/lfe_extension/results
NS="PGC=CN.PGC,VGZ=CN.VGZ,NLLB=CN.NLLB,SNB=CN.SNB,LRIV=UW.LRIV,B001=PB.B001,B003=PB.B003,B004=PB.B004,B005=PB.B005,B006=PB.B006,B007=PB.B007,B009=PB.B009,B010=PB.B010,B011=PB.B011,B926=PB.B926"
SITES="PGC:PGC,B009,B010,B011;LRIV:LRIV,B005,B006,B007"
THR=7,7.5,8,8.25,8.5,8.75,9,10
COMMON="--lengths 12 --n-random 600 --netsta $NS --no-seisbench --skip-z --half all \
  --families-csv $I/e5r_families_2012.csv --stations-csv $I/e5r_stations.csv \
  --extra-tpl $I/e5r_templates_2012.npz --moveout-in $I/e4_moveout_e5r_mv_build.csv \
  --weights $I/e5r_weights_snr_2012.csv --sites $SITES --site-merge B \
  --budget --budget-thr $THR --catalog-out --out-dir $O \
  --cache-dirs /home/jovyan/tremorferometry/lfe_extension/cache20hz"

run12 () {  # $1 tag  $2 stations
  python scripts/e4_template_test.py $COMMON --year 2012 --doy-from 228 --doy-to 279 --stations "$2" \
    --gt-parquet $I/e5r_gt_bostock2012.parquet --offsets-csv $I/e5r_family_offsets_2012.csv \
    --amp-ref-csv $I/e5r_amp_ref_2012.csv --tag "$1" > logs/$1.log 2>&1
  echo "$1 exit=$? $(date -u +%FT%TZ)" >> logs/d1_ablation.log
}
run10 () {  # $1 tag  $2 stations
  python scripts/e4_template_test.py $COMMON --year 2010 --doy-from 222 --doy-to 256 --stations "$2" \
    --gt-parquet $O/a2_gt_bostock2010.parquet --offsets-csv $I/e3_bostock_family_offsets.csv \
    --tag "$1" > logs/$1.log 2>&1
  echo "$1 exit=$? $(date -u +%FT%TZ)" >> logs/d1_ablation.log
}

E2_15=PGC,VGZ,NLLB,SNB,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
E3_15=PGC,VGZ,NLLB,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
E2_09=PGC,VGZ,SNB,NLLB,B926,B009,B010,B011
E3_09=PGC,VGZ,NLLB,B926,B009,B010,B011

echo "[$(date -u +%FT%TZ)] D1(a) ablation START" >> logs/d1_ablation.log
run12 d1_abl2012_E2 "$E2_15" &
run12 d1_abl2012_E3 "$E3_15" &
wait
run10 d1_abl2010_E2 "$E2_09" &
run10 d1_abl2010_E3 "$E3_09" &
wait
echo "[$(date -u +%FT%TZ)] D1(a) ablation DONE" >> logs/d1_ablation.log
