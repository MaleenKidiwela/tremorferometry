#!/bin/bash
# E5 step 1 -- DETECTION-BUDGET table on the 2010 data (Merlin E4-audit prescription 1).
#   pass 1  empirical moveout measured on the BUILD half only (alternate 2-h blocks, 60-s guard), all
#           10 stations: ring published stacks + build-half transfers B011/B926/LRIV/NLLB (e5 npz, 4/8/12 s)
#   pass 2  --budget on the FULL window (quiet days 222-229 and 249-255 have <= 50 Bostock events) for
#           (a) the OPERATING SUBSET PGC/VGZ/SNB + B011/B926/LRIV/NLLB (primary),
#           (b) the full ring + transfers (reference), (c) the E4 6-station ring alone (ties to E4's F1 numbers)
cd /home/jovyan/tremorferometry/lfe_catalog_repair
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages
set -e
R=results
X="$R/e4_transfer_stacks_e5.npz"
RING=PGC,LZB,VGZ,SNB,YOUB,PFB
OP=PGC,VGZ,SNB,B011,B926,LRIV,NLLB
ALL=$RING,B011,B926,LRIV,NLLB
XS=B011,B926,LRIV,NLLB
NS="PGC=CN.PGC,LZB=CN.LZB,VGZ=CN.VGZ,SNB=CN.SNB,YOUB=CN.YOUB,PFB=CN.PFB,B011=PB.B011,B926=PB.B926,LRIV=UW.LRIV,NLLB=CN.NLLB"
T="python scripts/e4_template_test.py --lengths 4,8,12 --n-random 600 --doy-from 222 --doy-to 255 --netsta $NS"
$T --tag e5_mv_build --half build --stations $ALL --extra-tpl $X --extra-tpl-stations $XS
MV=$R/e4_moveout_e5_mv_build.csv
$T --tag e5_budget_op   --stations $OP   --extra-tpl $X --extra-tpl-stations $XS --moveout-in $MV --budget
$T --tag e5_budget_ring --stations $ALL  --extra-tpl $X --extra-tpl-stations $XS --moveout-in $MV --budget
$T --tag e5_budget_ring6 --stations $RING --moveout-in $MV --budget
echo "[e5 step1 DONE]"
