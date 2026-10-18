#!/bin/bash
# E4 step 2 evaluation — held-out (alternate 2-h blocks) registration.
#   A  ring, PUBLISHED SeisBench templates                      -> the reference
#   B  ring, OUR 2010 re-stacks built on the BUILD half only    -> the measured q
#   C  ring published + TRANSFERRED B011/B926/LRIV/STW          -> the network gain
# Each leg measures its empirical moveout on the BUILD half and is scored on the TEST half, so
# nothing that sets the alignment has seen a test-half event.
cd /home/jovyan/tremorferometry/lfe_catalog_repair
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages
set -e
R=results
X="$R/e4_transfer_stacks_t2.npz"
RING=PGC,LZB,VGZ,SNB,YOUB,PFB
ALL=$RING,B011,B926,LRIV,STW
NS="PGC=CN.PGC,LZB=CN.LZB,VGZ=CN.VGZ,SNB=CN.SNB,YOUB=CN.YOUB,PFB=CN.PFB,B011=PB.B011,B926=PB.B926,LRIV=UW.LRIV,STW=UW.STW"
T="python scripts/e4_template_test.py --lengths 4,6,8 --n-random 600 --doy-from 222 --doy-to 255"

# A
$T --tag A_build --half build --stations $RING
$T --tag A_test  --half test  --stations $RING --moveout-in $R/e4_moveout_A_build.csv
# B
$T --tag B_build --half build --stations $RING --drop-stations $RING --extra-tpl $X --extra-tpl-stations $RING
$T --tag B_test  --half test  --stations $RING --drop-stations $RING --extra-tpl $X --extra-tpl-stations $RING --moveout-in $R/e4_moveout_B_build.csv
# C
$T --tag C_build --half build --stations $ALL --netsta $NS --sta-weight STW=1 --extra-tpl $X --extra-tpl-stations B011,B926,LRIV,STW
$T --tag C_test  --half test  --stations $ALL --netsta $NS --sta-weight STW=1 --extra-tpl $X --extra-tpl-stations B011,B926,LRIV,STW --moveout-in $R/e4_moveout_C_build.csv
echo "[e4 step2 eval DONE]"
