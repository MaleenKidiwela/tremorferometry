#!/bin/bash
# VLFE scan pilot 2 (S7.1): prep v3 of the on-disk reference stations. Marker SCAN2_PREP_REF_DONE.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe
export OMP_NUM_THREADS=1 S13_DISP=$PWD/data/disp1sps_v3
python scripts/scan_prep_ref.py UW.GNW,UW.DOSE 2011-01-01 2011-12-31 --workers 8 > logs/scan2_prep_ref.out 2>&1
python scripts/scan_prep_ref.py CC.OBSR,UW.LEBA,UW.RATT,UW.ELW,TA.D04D,UW.WISH,UW.SP2,TA.B05D,US.NLWA,TA.D03D 2011-07-30 2011-08-27 --workers 8 >> logs/scan2_prep_ref.out 2>&1
touch SCAN2_PREP_REF_DONE
