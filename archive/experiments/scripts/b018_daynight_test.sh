#!/bin/bash
# DECISIVE B018 test: day-only vs night-only detection stacks -> two independent dv/v series.
# Real medium signal => the two series match. Noise-source-carried => they diverge.
cd /home/jovyan/tremorferometry
export PYTHONPATH=src OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY=/home/jovyan/envs/tremorferometry/bin/python
L=logs/b018_daynight.log; : > $L
say(){ echo "=== $(date +%H:%M:%S) $* ===" | tee -a $L; }

say "split detections (local day 15-23 UTC / night 5-13 UTC)"
$PY - <<'EOF' >> $L 2>&1
import pandas as pd
m=pd.read_csv('data/mf_b018_all.csv')
hr=pd.to_datetime(m.time).dt.hour
m[hr.between(15,23)].to_csv('data/mf_b018_dayonly.csv',index=False)
m[hr.between(5,13)].to_csv('data/mf_b018_nightonly.csv',index=False)
print('day',hr.between(15,23).sum(),'night',hr.between(5,13).sum())
EOF

for ARM in dayonly nightonly; do
  say "stack $ARM"
  [ -f data/long_window_daily_B018_${ARM}.npz ] || $PY scripts/build_long_window_resp.py \
    --mf-csv data/mf_b018_${ARM}.csv --network PB --station B018 --no-deconv --min-det 8 \
    --despike-mad 8 --workers 16 --out data/long_window_daily_B018_${ARM}.npz >> $L 2>&1
  say "dv/v $ARM (2-4s, S-anchored)"
  [ -f data/daily_dvv_B018_2to4_${ARM}.csv ] || $PY scripts/dvv_roll30cal.py --station B018 \
    --npz data/long_window_daily_B018_${ARM}.npz --window 2.0 4.0 --origin-anchor \
    --out data/daily_dvv_B018_2to4_${ARM}.csv --workers 16 >> $L 2>&1
done

say "compare arms"
$PY - <<'EOF' 2>&1 | tee -a $L
import pandas as pd, numpy as np
D=pd.read_csv('data/daily_dvv_B018_2to4_dayonly.csv',parse_dates=['date']).groupby('date').dvv.median()*100
N=pd.read_csv('data/daily_dvv_B018_2to4_nightonly.csv',parse_dates=['date']).groupby('date').dvv.median()*100
a,b=D.align(N,join='inner')
r=np.corrcoef(a,b)[0,1]; slope=np.sum(a*b)/np.sum(a*a)
print(f'day-vs-night dv/v: overlap {len(a)} days | corr {r:+.3f} | slope {slope:.2f}')
print(f'  day std {a.std():.3f}%  night std {b.std():.3f}%')
# the 2016-17 drop in both?
for t0,t1 in [('2016-09-01','2017-03-31')]:
    da=a[(a.index>=t0)&(a.index<=t1)].min(); na=b[(b.index>=t0)&(b.index<=t1)].min()
    print(f'  2016-17 drop: day-only min {da:+.3f}%  night-only min {na:+.3f}%')
print('VERDICT: corr>~0.7 + drop in BOTH = medium-carried (REAL). corr low / drop one-sided = source-carried (artifact).')
EOF
say "DAYNIGHT TEST DONE"
