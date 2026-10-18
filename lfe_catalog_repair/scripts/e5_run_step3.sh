#!/bin/bash
# E5 step 3 -- SNR-WEIGHTED network on 2010, HELD OUT (alternate 2-h blocks; moveout + templates + weights from the
# BUILD half only).  Weights = results/e5_weights_snr.csv (build-half stack S-SNR, 10^(dB/20), capped 3x median,
# median-normalised); the rolled control is recomputed WITH the weights inside e4_template_test.py.
# Operating subset (PGC/VGZ/SNB + B011/B926/LRIV/NLLB) primary; full ring + transfers reference.
cd /home/jovyan/tremorferometry/lfe_catalog_repair
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages
set -e
R=results; X="$R/e4_transfer_stacks_e5.npz"; MV=$R/e4_moveout_e5_mv_build.csv; W=$R/e5_weights_snr.csv
OP=PGC,VGZ,SNB,B011,B926,LRIV,NLLB; ALL=PGC,LZB,VGZ,SNB,YOUB,PFB,B011,B926,LRIV,NLLB; XS=B011,B926,LRIV,NLLB
NS="PGC=CN.PGC,LZB=CN.LZB,VGZ=CN.VGZ,SNB=CN.SNB,YOUB=CN.YOUB,PFB=CN.PFB,B011=PB.B011,B926=PB.B926,LRIV=UW.LRIV,NLLB=CN.NLLB"
T="python scripts/e4_template_test.py --lengths 8,12 --n-random 600 --doy-from 222 --doy-to 255 --netsta $NS --half test --moveout-in $MV --budget --extra-tpl $X --extra-tpl-stations $XS"
$T --tag e5_w_op_eq    --stations $OP
$T --tag e5_w_op_snr   --stations $OP  --weights $W
$T --tag e5_w_ring_eq  --stations $ALL
$T --tag e5_w_ring_snr --stations $ALL --weights $W
echo "[e5 step3 DONE]"
