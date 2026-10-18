#!/bin/bash
# E6 GPU chain 1 (2026-09-11, Fable agent; GPU goes offline in a few hours -> all GPU passes first, 8-s templates only, checkpointed).
#   pass 1  step 1  : UNCHANGED 15-station SNR-weighted network + SITE-rolled control (--sites), per-site contributions   -> tag e6_sites
#   pass 2  step 2A : site-merged network, site weight = max member w                                                       -> tag e6_mergeA
#   pass 3  step 2B : site-merged network, site weight = sqrt(sum w^2)                                                      -> tag e6_mergeB
#   pass 4  step 1b : the 5 worst quiet-FA families, VGZ and SNB EXCLUDED from the sum, all 52 days (quiet = test, ETS = +control) -> e6_1b_noVGZSNB
# Same templates / weights / moveout / thresholds as E5r (e5r_run_mf2.sh) plus --budget-thr 7,7.5,8,8.5,9,10.  Each pass writes its own
# budget csv + detections parquet; results/e6_checkpoint.txt lists the completed passes.  A CUDA failure is logged and the chain moves on.
cd /home/jovyan/tremorferometry/lfe_catalog_repair
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages CUDA_PATH=/opt/conda/targets/x86_64-linux OMP_NUM_THREADS=4
R=results; CK=$R/e6_checkpoint.txt
ST=PGC,VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
ST13=PGC,NLLB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
NS="PGC=CN.PGC,VGZ=CN.VGZ,NLLB=CN.NLLB,SNB=CN.SNB,LRIV=UW.LRIV,B001=PB.B001,B003=PB.B003,B004=PB.B004,B005=PB.B005,B006=PB.B006,B007=PB.B007,B009=PB.B009,B010=PB.B010,B011=PB.B011,B926=PB.B926"
SITES="PGC:PGC,B009,B010,B011;LRIV:LRIV,B005,B006,B007"
COMMON="--year 2012 --lengths 8 --n-random 600 --doy-from 228 --doy-to 279 --netsta $NS --no-seisbench \
  --families-csv $R/e5r_families_2012.csv --stations-csv $R/e5r_stations.csv --gt-parquet $R/e5r_gt_bostock2012.parquet \
  --offsets-csv $R/e5r_family_offsets_2012.csv --extra-tpl $R/e5r_templates_2012.npz --skip-z --half all \
  --moveout-in $R/e4_moveout_e5r_mv_build.csv --budget --budget-thr 7,7.5,8,8.5,9,10 --catalog-out --amp-ref-csv $R/e5r_amp_ref_2012.csv \
  --weights $R/e5r_weights_snr_2012.csv --sites $SITES"
run () {  # tag, then extra args
  tag=$1; shift
  echo "[$(date -u +%FT%TZ)] START $tag" | tee -a logs/e6_gpu_chain.log
  if python scripts/e4_template_test.py $COMMON --tag $tag "$@" > logs/e6_mf_$tag.log 2>&1; then
    echo "$tag DONE $(date -u +%FT%TZ)" >> $CK; echo "[$(date -u +%FT%TZ)] DONE $tag" | tee -a logs/e6_gpu_chain.log
  else
    echo "$tag FAILED $(date -u +%FT%TZ) (see logs/e6_mf_$tag.log)" >> $CK; echo "[$(date -u +%FT%TZ)] FAILED $tag" | tee -a logs/e6_gpu_chain.log; tail -20 logs/e6_mf_$tag.log >> logs/e6_gpu_chain.log
  fi
}
run e6_sites  --stations $ST
run e6_mergeA --stations $ST --site-merge A
run e6_mergeB --stations $ST --site-merge B
run e6_1b_noVGZSNB --stations $ST13 --only-fams 103,32,75,128,74
echo "[$(date -u +%FT%TZ)] CHAIN 1 DONE" | tee -a logs/e6_gpu_chain.log
