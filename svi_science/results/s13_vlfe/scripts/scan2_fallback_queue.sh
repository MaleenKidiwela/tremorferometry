#!/bin/bash
# S7.12/S7.13: columns (with the far references) on the uncleaned fallback runs, then validate2. Marker SCAN2_FALLBACK_DONE.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe
export OMP_NUM_THREADS=1
L=logs/scan2_fallback.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') $*" >> $L; }
log start
python scripts/scan_columns2.py real --bp data/scan/bp_2011g --refbp data/scan/bp2_ref --workers 8 >> $L 2>&1; log "real rc $?"
for s in 3 5 7 11 13; do python scripts/scan_columns2.py ctrlA_s$s --bp data/scan/bp_2011g --refbp data/scan/bp2_ref --workers 8 >> $L 2>&1; log "ctrlA_s$s rc $?"; done
python scripts/scan_columns2.py inj2 --bp data/scan/bp_2011inj2 --refbp data/scan/bp2_ref_inj2 --inj data/scan/inject_2011.parquet --workers 8 >> $L 2>&1; log "inj2 rc $?"
python scripts/scan_columns2.py inj3_u --bp data/scan/bp_2011inj3 --refbp data/scan/bp2_ref_inj3 --inj data/scan/inject3_2011.parquet --workers 8 >> $L 2>&1; log "inj3_u rc $?"
python scripts/scan_validate2.py > logs/scan2_validate.out 2>&1; log "validate2 rc $?"
touch SCAN2_FALLBACK_DONE
