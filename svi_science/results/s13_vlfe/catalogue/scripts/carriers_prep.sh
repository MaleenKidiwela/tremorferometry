#!/bin/bash
# Amendment S (b): make the 2011 carriers usable. Marker CARRIERS_DONE. FAIL lines -> logs/q_main.log.
D=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe; cd $D
export OMP_NUM_THREADS=1 S13_DISP=$D/data/disp1sps_v3
L=$D/catalogue/logs/q_main.log
log(){ echo "$(TZ=America/Los_Angeles date '+%F %H:%M:%S PT') carriers $*" >> $L; }
fail(){ log "FAIL $*"; exit 1; }
python - >> $D/catalogue/logs/carriers.out 2>&1 <<'PY' || fail link
import os, glob
D = '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe'
n = 0
for s in ('UW.SP2', 'TA.A04D', 'TA.D03D'):
    for f in glob.glob(f'{D}/hg_repro/data/disp1sps/{s}/*/*.npz'):
        y, doy = f.split('/')[-2], f.split('/')[-1]
        t = f'{D}/data/disp1sps_v3/{s}/{y}/{doy}'
        if not os.path.exists(t):
            os.makedirs(os.path.dirname(t), exist_ok=True); os.symlink(f, t); n += 1
print('linked', n)
PY
log "linked hg_repro v3-recipe day files (see carriers.out)"
python scripts/scan_prep_ref.py UW.SP2,TA.A04D,TA.D03D,UW.WISH,US.NLWA 2011-01-01 2012-12-31 --workers 8 > /dev/null 2>&1 || fail prep
log "prep done"
for Y in 2011 2012; do
  SD=$D/data/scan_$Y; [ $Y = 2011 ] && SD=$D/data/scan
  mkdir -p $SD/bp2_ref_v2
  for f in $SD/bp2_ref/*; do b=$(basename $f); case $b in UW.SP2*|TA.D03D*|UW.WISH*|US.NLWA*) ;; *) [ -e $SD/bp2_ref_v2/$b ] || ln -s $f $SD/bp2_ref_v2/$b;; esac; done
  SCAN_YEAR=$Y python scripts/scan_bp_ref.py $SD/bp2_ref_v2 --stations UW.SP2,TA.D03D,UW.WISH,US.NLWA --workers 4 >> $D/catalogue/logs/carriers.out 2>&1 || fail bp_ref_$Y
  log "bp2_ref_v2 $Y done"
done
TZ=America/Los_Angeles date > $D/catalogue/CARRIERS_DONE; log done
