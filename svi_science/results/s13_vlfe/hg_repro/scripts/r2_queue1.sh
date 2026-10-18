#!/bin/bash
# Amendment R2 stage 1: velocity arrays, then the pre-listed variants V0-V6 for T1-T5 (5 workers). Marker R2_VAR_DONE.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe/hg_repro
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python scripts/hg_r2.py vel CN.NLLB,CN.VGZ,UW.LRIV,UW.DOSE,UW.GNW,UW.SP2,TA.A04D,TA.B05D,TA.D03D || exit 1
for K in 1 2 3 4 5; do python scripts/hg_r2.py var $K & done
wait
n=$(ls data/r2/var_T?.npz 2>/dev/null | wc -l)
[ "$n" = "5" ] && touch R2_VAR_DONE
echo "$(TZ=America/Los_Angeles date) queue1 end, var files $n"
