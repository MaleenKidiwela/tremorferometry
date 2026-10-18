#!/bin/bash
# Q0: prep + bp + reference templates + refqc for US.MSO/US.WVOR in a scan year. usage: q0_refs.sh YEAR
Y=$1; D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D
export OMP_NUM_THREADS=1 SCAN_YEAR=$Y S13_DISP=$D/data/disp1sps_v3
SD=$D/data/scan_$Y; L=$D/catalogue/logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') q0_refs $Y $*" >> $L; }
fail(){ log "FAIL $*"; exit 1; }
python scripts/scan_prep_ref.py US.MSO,US.WVOR $Y-01-01 $Y-12-31 --workers 8 > /dev/null 2>&1 || fail prep
python scripts/scan_bp_ref.py $SD/bp2_ref --stations US.MSO,US.WVOR --workers 2 >> $L 2>&1 || fail bp_ref
if [ -L $SD/templates_ref.npz ]; then mv $SD/templates_ref.npz $SD/templates_ref.npz.symlink_2011; fi
[ -f $SD/templates_ref.npz ] || python scripts/scan_columns2.py templates >> $L 2>&1 || fail templates
python scripts/scan_refqc.py >> $L 2>&1 || fail refqc
log done
