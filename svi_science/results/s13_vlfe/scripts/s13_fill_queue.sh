#!/bin/bash
# S13 gap fill (PLAN.md Amendment 3; user OK 2026-10-13). Idempotent. Run: nohup bash scripts/s13_fill_queue.sh > logs/fill_queue.out 2>&1 &
# 1 fetch raw day files (<=3 parallel, resumable) -> data/raw/   2 wait for the scale-up queue (SCALE_QUEUE_DONE)
# 3 prep fetched station-days -> 1-sps displacement data/disp1sps_v2/   4 re-stack episodes/between sets touching them (tag scalef_*)
# 5 summary (prefers scalef_ tables). Markers: data/raw/FETCH_fill_DONE, PREP_FILL_DONE, STACK_scalef_*_DONE, FILL_QUEUE_DONE.
cd /home/jovyan/tremorferometry/svi_science/results/s13_vlfe
export OMP_NUM_THREADS=1
export S13_DISP=/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/data/disp1sps_v2
LOG=logs/fill_queue.log
ts(){ TZ=America/Los_Angeles date '+%Y-%m-%d %H:%M PT'; }
echo "$(ts) fill queue start" >> $LOG
until [ -f data/raw/FETCH_fill_DONE ]; do
  python scripts/s13_fetch.py tables/fetch_fill.csv --parallel 3 --tag fill >> logs/fetch_fill.out 2>&1
  [ -f data/raw/FETCH_fill_DONE ] || { echo "$(ts) fetch incomplete, retry in 10 min" >> $LOG; sleep 600; }
done
echo "$(ts) fetch done" >> $LOG
until [ -f SCALE_QUEUE_DONE ]; do sleep 120; done
SVI_STA=CN.PGC,CN.VGZ,CN.SNB,CN.NLLB,CN.LZB,CN.PFB,CN.YOUB,CN.SYMB,CN.CLRS,CN.PTRF,C8.SOKB,C8.SHVB,C8.KHVB,C8.MGCB,C8.TWBB,C8.SHDB,C8.BPCB,UW.LRIV,UW.OPC,TA.B04A,UW.SQM,UW.ODUC,UW.OSQM,UW.SLDQ,UW.ODEER,UW.OTR,XU.BS11,CN.GOBB
WA_STA=UW.DOSE,UW.GNW,UW.LRIV,UW.SQM,UW.OSQM,UW.ODEER,UW.NATEM,UW.CHIMA,UW.EQUIL,UW.HDW,UW.HANS,UW.COYL,UW.KNGST,UW.KTSAP,UW.SKOKO,UW.MORSE,UW.WISH,UW.SP2,UW.OPC,UW.PRIN,UW.WINDI,UW.EAGLE,US.NLWA,TA.D03D,TA.D04E,TA.D04D,TA.B04A,YW.FACB,XU.W020,XU.W030,XU.W040,XU.C04A,XU.BS11,XU.PL11,XU.W060,XU.N040,XU.N050,XU.N060,CN.PGC,CN.VGZ
if [ ! -f PREP_FILL_DONE ]; then
  python - <<'PY'
import pandas as pd
L = pd.read_csv('tables/fetch_fill.csv'); L[['date']].drop_duplicates().to_csv('tables/fill_days.csv', index=False)
open('tables/fill_stations.txt', 'w').write(','.join(sorted(L.netsta.unique())))
PY
  python scripts/s13_prep.py $(cat tables/fill_stations.txt) x x --days-file tables/fill_days.csv --workers 8 > logs/prep_fill.out 2>&1 && ! grep -q "ModuleNotFoundError" logs/prep_fill.out && touch PREP_FILL_DONE
  echo "$(ts) prep fill done" >> $LOG
fi
# affected row files: any needed day (+-1) among the fetched days with a non-'nodata' result
python - <<'PY' > tables/fill_affected.txt
import pandas as pd, glob, numpy as np, os
ok = []
for l in open('logs/fetch_fill.log'):
    p = l.split()
    if len(p) >= 5 and p[4] == 'ok':
        ok.append(p[3])
ok = set(ok)
for f in sorted(glob.glob('data/rows/*.parquet')):
    R = pd.read_parquet(f, columns=['date'])
    d = set(pd.to_datetime(R.date).dt.strftime('%Y-%m-%d'))
    d |= {(pd.Timestamp(x) + pd.Timedelta(days=k)).strftime('%Y-%m-%d') for x in list(d) for k in (-1, 1)}
    if d & ok:
        print(f)
PY
echo "$(ts) affected row files: $(wc -l < tables/fill_affected.txt)" >> $LOG
n=0
for f in $(cat tables/fill_affected.txt); do
  b=$(basename $f .parquet); tag=scalef_$b
  [ -f STACK_${tag}_DONE ] && continue
  STA=$([[ $b == svi_* ]] && echo $SVI_STA || echo $WA_STA)
  ( python scripts/s13_stack2.py $f $tag $STA > logs/${tag}.out 2>&1 && echo "$(ts) $tag done" >> $LOG || echo "$(ts) $tag FAILED" >> $LOG ) &
  n=$((n+1)); if [ $n -ge 4 ]; then wait -n; n=$((n-1)); fi
done
wait
miss=0; for f in $(cat tables/fill_affected.txt); do [ -f STACK_scalef_$(basename $f .parquet)_DONE ] || miss=$((miss+1)); done   # guard added 2026-10-14
python scripts/s13_scale_summary.py > logs/scale_summary_fill.out 2>&1 && [ $miss -eq 0 ] && touch FILL_QUEUE_DONE
[ $miss -eq 0 ] || echo "$(ts) NOT done: $miss scalef stacks missing" >> $LOG
echo "$(ts) fill queue end" >> $LOG
