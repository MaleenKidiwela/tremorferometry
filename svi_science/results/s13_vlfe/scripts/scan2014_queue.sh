#!/bin/bash
# Scan year 2014 for the catalogue pilot (catalogue/PLAN_CAT.md Amendment P). Idempotent (STEP_*_DONE). <= 8 workers.
#   nohup bash scripts/scan2014_queue.sh > logs/scan2014_queue.out 2>&1 &     marker data/scan_2014/SCAN2014_DONE
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe
cd $D
export OMP_NUM_THREADS=1 SCAN_YEAR=2014 S13_DISP=$D/data/disp1sps_v3
for p in /opt/conda/lib/python3.13/site-packages/nvidia/*/lib; do export LD_LIBRARY_PATH=$p:$LD_LIBRARY_PATH; done
SD=$D/data/scan_2014; mkdir -p $SD/runs
L=$D/logs/scan2014_queue.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') $*" >> $L; }
PIL=CN.PGC,CN.VGZ,CN.SNB,CN.NLLB,CN.PFB,CN.YOUB,UW.LRIV,CN.LZB
REG=UW.GNW,UW.DOSE,UW.SP2,TA.D03D,TA.B05D,TA.A04D
FAR=IU.COR,US.NEW,US.HAWA
REFS=UW.GNW,UW.DOSE,UW.SP2,TA.D03D,TA.B05D,$FAR
W=8
log "queue start"
[ -e $SD/templates_grid.npz ] || ln -s $D/data/scan/templates_grid.npz $SD/templates_grid.npz
[ -e $SD/templates_ref.npz ] || ln -s $D/data/scan/templates_ref.npz $SD/templates_ref.npz
[ -e $SD/hfenv ] || ln -s $D/data/scan/hfenv $SD/hfenv
step(){ local n=$1; shift; [ -f $SD/STEP_${n}_DONE ] && return 0
  log "start $n"; "$@" >> $L 2>&1; local rc=$?; log "end $n rc $rc"; [ $rc -eq 0 ] && touch $SD/STEP_${n}_DONE; return $rc; }
step prep_pilot python scripts/scan_prep_v3.py $PIL 2014-01-01 2014-12-31 --workers $W || exit 1
step prep_reg python scripts/scan_prep_ref.py $REG 2014-01-01 2014-12-31 --workers $W || exit 1
step prep_far python scripts/scan_prep_ref.py $FAR 2014-01-01 2014-12-31 --workers $W || exit 1
step noise python scripts/scan3_noise.py || exit 1
step masks python scripts/scan_masks2.py || exit 1
step bp python scripts/scan_bp.py $SD/bp_2014g --workers $W || exit 1
step bp_ref python scripts/scan_bp_ref.py $SD/bp2_ref --stations $REFS --workers $W || exit 1
step refqc python scripts/scan_refqc.py || exit 1
[ -f $SD/runs/real/DONE ] || { python scripts/scan_detect.py real --bp $SD/bp_2014g --ndraw 0 --shift 0 >> $L 2>&1; log "detect real rc $?"; }
[ -f $SD/runs/real/DONE ] && touch $SD/SCAN2014_DONE
log "queue end"
