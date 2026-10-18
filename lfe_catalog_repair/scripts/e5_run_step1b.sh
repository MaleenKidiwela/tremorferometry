#!/bin/bash
# E5 step 1, pass 2 only (pass 1 = e5_mv_build finished; the first pass-2 attempt died on KeyError (6.0,'B011'), fixed).
cd /home/jovyan/tremorferometry/lfe_catalog_repair
export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages
set -e
R=results; X="$R/e4_transfer_stacks_e5.npz"; MV=$R/e4_moveout_e5_mv_build.csv
RING=PGC,LZB,VGZ,SNB,YOUB,PFB; OP=PGC,VGZ,SNB,B011,B926,LRIV,NLLB; ALL=$RING,B011,B926,LRIV,NLLB; XS=B011,B926,LRIV,NLLB
NS="PGC=CN.PGC,LZB=CN.LZB,VGZ=CN.VGZ,SNB=CN.SNB,YOUB=CN.YOUB,PFB=CN.PFB,B011=PB.B011,B926=PB.B926,LRIV=UW.LRIV,NLLB=CN.NLLB"
T="python scripts/e4_template_test.py --lengths 4,8,12 --n-random 600 --doy-from 222 --doy-to 255 --netsta $NS --moveout-in $MV --budget"
$T --tag e5_budget_op    --stations $OP   --extra-tpl $X --extra-tpl-stations $XS
$T --tag e5_budget_ring  --stations $ALL  --extra-tpl $X --extra-tpl-stations $XS
$T --tag e5_budget_ring6 --stations $RING
echo "[e5 step1 DONE]"
