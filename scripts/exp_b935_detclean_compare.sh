#!/bin/bash
# Wait for the cleaned-stack build, run dv/v (2-4s, origin-anchor true-scale) on RAW and CLEANED,
# then compare GOLD-sharpens vs FAIL-stays-noise. Autonomous tail of the B935 detection-cleaning test.
set -u
cd /home/jovyan/tremorferometry
export PYTHONPATH=src OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY=/home/jovyan/envs/tremorferometry/bin/python
L=logs/b935_detclean_compare.log; : > $L
echo "$(date +%H:%M) waiting for BUILD_DONE..." | tee -a $L
for i in $(seq 1 240); do grep -q BUILD_DONE logs/b935_detclean_build.log 2>/dev/null && break; sleep 30; done
grep -q BUILD_DONE logs/b935_detclean_build.log || { echo "BUILD never finished — abort" | tee -a $L; exit 1; }
echo "$(date +%H:%M) build done; running dv/v on RAW and CLEANED" | tee -a $L
for KIND in RAW CLEAN; do
  $PY scripts/dvv_roll30cal.py --station B935 --npz data/long_window_daily_B935_${KIND}exp.npz \
     --window 2.0 4.0 --origin-anchor --workers 24 --out data/daily_dvv_B935_${KIND}exp_2to4_calT.csv >> $L 2>&1
done
echo "$(date +%H:%M) comparing" | tee -a $L
$PY - <<'EOF' 2>&1 | tee -a $L
import pandas as pd, numpy as np
GOLD=['40.800_-123.450__c18','40.050_-122.600__c130']; FAIL=['40.550_-122.800__c406','40.050_-122.600__c193']
raw=pd.read_csv('data/daily_dvv_B935_RAWexp_2to4_calT.csv'); cln=pd.read_csv('data/daily_dvv_B935_CLEANexp_2to4_calT.csv')
def seasonal(g):
    # fraction of dv/v variance explained by an annual sinusoid (clarity of seasonal signal)
    d=pd.to_datetime(g.date); doy=d.dt.dayofyear.values; y=g.dvv.values
    if len(y)<40: return np.nan, np.nan
    X=np.c_[np.ones_like(doy,float),np.cos(2*np.pi*doy/365.25),np.sin(2*np.pi*doy/365.25)]
    b,_,_,_=np.linalg.lstsq(X,y,rcond=None); pred=X@b; ss=1-np.sum((y-pred)**2)/np.sum((y-y.mean())**2)
    amp=np.hypot(b[1],b[2]); return ss, amp
print(f"\n{'family':24s} {'set':5s} {'npts':>5s} {'cc_max':>7s} {'dvv_std':>8s} {'seasR2':>7s} {'seasAmp':>8s}")
for fam in GOLD+FAIL:
    tag='GOLD' if fam in GOLD else 'FAIL'
    for nm,df in [('RAW',raw),('CLN',cln)]:
        g=df[df.patch==fam]
        if not len(g): print(f"{fam:24s} {nm:5s}  (none)"); continue
        r2,amp=seasonal(g)
        print(f"{fam:24s} {nm:5s} {len(g):5d} {g.cc_max.mean():7.3f} {g.dvv.std()*100:7.3f}% {r2:7.2f} {amp*100:7.3f}%   [{tag}]")
print("\nVERDICT GUIDE: GOLD should gain cc_max & seasonal-R2 and/or lose scatter after cleaning;")
print("FAIL should NOT gain a seasonal (else cleaning manufactures signal = circular). COMPARE_DONE")
EOF
echo "$(date +%H:%M) EXPERIMENT COMPLETE" | tee -a $L