#!/bin/bash
# Scan year Y for the catalogue pilots (PLAN_CAT Amendments P/Q). usage: scan_year_queue.sh YEAR
# Idempotent (STEP_*_DONE in data/scan_Y). Marker data/scan_Y/SCANY_DONE. FAIL lines -> catalogue/logs/q_main.log.
Y=$1; D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D
export OMP_NUM_THREADS=1 SCAN_YEAR=$Y S13_DISP=$D/data/disp1sps_v3
unset SCAN_GRID SCAN_GAIN_FALLBACK SCAN_RUN SCAN_BP SCAN_REFBP   # v1 scan year: no v2 overrides (a leaked SCAN_GRID broke 2020)
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
SD=$D/data/scan_$Y; mkdir -p $SD/runs
L=$D/logs/scan${Y}_queue.log; M=$D/catalogue/logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') $*" >> $L; }
fail(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') scan$Y FAIL $*" >> $M; exit 1; }
PIL=CN.PGC,CN.VGZ,CN.SNB,CN.NLLB,CN.PFB,CN.YOUB,UW.LRIV,CN.LZB
REG=UW.GNW,UW.DOSE,UW.SP2,TA.D03D,TA.B05D,TA.A04D,UW.WISH,US.NLWA${SCAN_REG_EXTRA}
FAR=IU.COR,US.NEW,US.HAWA,US.MSO,US.WVOR
REFS=UW.GNW,UW.DOSE,UW.SP2,TA.D03D,TA.B05D,UW.WISH,US.NLWA${SCAN_REG_EXTRA},$FAR
W=${CAT_NW:-8}
[ -e $SD/templates_grid.npz ] || ln -s $D/data/scan/templates_grid.npz $SD/templates_grid.npz
[ -e $SD/hfenv ] || ln -s $D/data/scan/hfenv $SD/hfenv
step(){ local n=$1; shift; [ -f $SD/STEP_${n}_DONE ] && return 0
  log "start $n"; "$@" >> $L 2>&1; local rc=$?; log "end $n rc $rc"; [ $rc -eq 0 ] && touch $SD/STEP_${n}_DONE; return $rc; }
step prep_pilot python scripts/scan_prep_v3.py $PIL $Y-01-01 $Y-12-31 --workers $W || fail prep_pilot
step prep_reg python scripts/scan_prep_ref.py $REG $Y-01-01 $Y-12-31 --workers $W || fail prep_reg
step prep_far python scripts/scan_prep_ref.py $FAR $Y-01-01 $Y-12-31 --workers $W || fail prep_far
step noise python scripts/scan3_noise.py || fail noise
step masks python scripts/scan_masks2.py || fail masks
step bp python scripts/scan_bp.py $SD/bp_${Y}g --workers $W || fail bp
step bp_ref python scripts/scan_bp_ref.py $SD/bp2_ref --stations $REFS --workers $W || fail bp_ref
step templates_ref python scripts/scan_columns2.py templates || fail templates_ref
step refqc python scripts/scan_refqc.py || fail refqc
[ -f $SD/runs/real/DONE ] || { python scripts/scan_detect.py real --bp $SD/bp_${Y}g --ndraw 0 --shift 0 >> $L 2>&1; log "detect real rc $?"; }
[ -f $SD/runs/real/DONE ] || fail detect
touch $SD/SCAN${Y}_DONE; log "queue end"
