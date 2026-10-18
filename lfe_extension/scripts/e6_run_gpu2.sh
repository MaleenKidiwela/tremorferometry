#!/bin/bash
# E6 GPU chain 2 (2026-09-11, Fable agent): step-3 recall LEVERS at the step-2 operating variant ($1 = A or B), then the step-2 FALLBACKS
# as labelled extras.  8-s templates only.  Each lever = the site-merged network of the adopted variant with ONE change:
#   (a) e6_lev_a : templates rebuilt at LRIV/B001-B007 with the +-8 s residual-centred S pick (e6_templates_2012_a.npz + weights);
#                  its own pass-1 moveout (e6_lev_a_mv, build half) for method parity with E5r
#   (b) e6_lev_b : network restricted to the 8 stations with median template S-SNR >= 5 dB (PGC B010 B009 B011 SNB VGZ B926 NLLB)
#   (c) e6_lev_c : 2010+2012 combined templates at the 7 persistent stations (e6_templates_2012_c.npz + weights); own pass-1 moveout
#   extras       : e6_merge${V}_cap25 (PGC site <= 25 % of the total weight); e6_merge${V}_coinc (>= 2 sites each >= 3 MAD)
# results/e6_checkpoint.txt lists completed passes; a CUDA failure is logged and the chain moves on.
V=${1:-A}
cd /home/jovyan/tremorferometry/lfe_catalog_repair
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages CUDA_PATH=/opt/conda/targets/x86_64-linux OMP_NUM_THREADS=4
R=results; CK=$R/e6_checkpoint.txt
ST=PGC,VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926
ST8=PGC,VGZ,NLLB,SNB,B009,B010,B011,B926
NS="PGC=CN.PGC,VGZ=CN.VGZ,NLLB=CN.NLLB,SNB=CN.SNB,LRIV=UW.LRIV,B001=PB.B001,B003=PB.B003,B004=PB.B004,B005=PB.B005,B006=PB.B006,B007=PB.B007,B009=PB.B009,B010=PB.B010,B011=PB.B011,B926=PB.B926"
SITES="PGC:PGC,B009,B010,B011;LRIV:LRIV,B005,B006,B007"
BASE="--year 2012 --lengths 8 --n-random 600 --doy-from 228 --doy-to 279 --netsta $NS --no-seisbench \
  --families-csv $R/e5r_families_2012.csv --stations-csv $R/e5r_stations.csv --gt-parquet $R/e5r_gt_bostock2012.parquet \
  --offsets-csv $R/e5r_family_offsets_2012.csv --skip-z"
FULL="--half all --budget --budget-thr 7,7.5,8,8.5,9,10 --catalog-out --amp-ref-csv $R/e5r_amp_ref_2012.csv --sites $SITES --site-merge $V"
run () {  # tag, then extra args
  tag=$1; shift
  echo "[$(date -u +%FT%TZ)] START $tag" | tee -a logs/e6_gpu_chain.log
  if python scripts/e4_template_test.py $BASE --tag $tag "$@" > logs/e6_mf_$tag.log 2>&1; then
    echo "$tag DONE $(date -u +%FT%TZ)" >> $CK; echo "[$(date -u +%FT%TZ)] DONE $tag" | tee -a logs/e6_gpu_chain.log
  else
    echo "$tag FAILED $(date -u +%FT%TZ) (see logs/e6_mf_$tag.log)" >> $CK; echo "[$(date -u +%FT%TZ)] FAILED $tag" | tee -a logs/e6_gpu_chain.log; tail -20 logs/e6_mf_$tag.log >> logs/e6_gpu_chain.log
  fi
}
echo "chain 2 variant $V" >> $CK
# (a) pass 1 (moveout, build half, no budget) then pass 2
run e6_lev_a_mv --stations $ST --extra-tpl $R/e6_templates_2012_a.npz --half build
run e6_lev_a    --stations $ST --extra-tpl $R/e6_templates_2012_a.npz --weights $R/e6_weights_snr_2012_a.csv --moveout-in $R/e4_moveout_e6_lev_a_mv.csv $FULL
# (b) >= 5 dB stations, E5r templates / weights / moveout
run e6_lev_b    --stations $ST8 --extra-tpl $R/e5r_templates_2012.npz --weights $R/e5r_weights_snr_2012.csv --moveout-in $R/e4_moveout_e5r_mv_build.csv $FULL
# (c) combined templates
run e6_lev_c_mv --stations $ST --extra-tpl $R/e6_templates_2012_c.npz --half build
run e6_lev_c    --stations $ST --extra-tpl $R/e6_templates_2012_c.npz --weights $R/e6_weights_snr_2012_c.csv --moveout-in $R/e4_moveout_e6_lev_c_mv.csv $FULL
# step-2 fallbacks as labelled extras
run e6_merge${V}_cap25 --stations $ST --extra-tpl $R/e5r_templates_2012.npz --weights $R/e5r_weights_snr_2012.csv --moveout-in $R/e4_moveout_e5r_mv_build.csv $FULL --site-cap PGC=0.25
run e6_merge${V}_coinc --stations $ST --extra-tpl $R/e5r_templates_2012.npz --weights $R/e5r_weights_snr_2012.csv --moveout-in $R/e4_moveout_e5r_mv_build.csv $FULL --coinc-gate 2,3
echo "[$(date -u +%FT%TZ)] CHAIN 2 DONE" | tee -a logs/e6_gpu_chain.log
