#!/bin/bash
# S7.13: prep v3 + band-pass (real, inj2, inj3) for the far references. Marker SCAN2_FAR_PREP_DONE.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe
export OMP_NUM_THREADS=1 S13_DISP=$PWD/data/disp1sps_v3
F=IU.COR,US.NEW,US.HAWA,UW.OMAK,TA.I02D
python scripts/scan_prep_ref.py $F 2011-01-01 2011-12-31 --workers 8 > logs/scan2_far_prep.out 2>&1
python scripts/scan_bp_ref.py data/scan/bp2_ref --stations $F --workers 5 >> logs/scan2_far_prep.out 2>&1
python scripts/scan_bp_ref.py data/scan/bp2_ref_inj2 --stations $F --inject data/scan/inject_2011.parquet --workers 5 >> logs/scan2_far_prep.out 2>&1
python scripts/scan_bp_ref.py data/scan/bp2_ref_inj3 --stations $F --inject data/scan/inject3_2011.parquet --workers 5 >> logs/scan2_far_prep.out 2>&1
touch SCAN2_FAR_PREP_DONE
