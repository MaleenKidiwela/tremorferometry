#!/bin/bash
# After the baseline 35-station inversion finishes, run the CONTINUOUS-ONLY variant and A/B compare.
cd /home/jovyan/tremorferometry
export PYTHONPATH=src OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY=/home/jovyan/envs/tremorferometry/bin/python
L=logs/ab_continuous.log; : > $L

# wait for the consolidate35 chain (baseline inversion) to finish
while pgrep -f "[c]onsolidate35.sh|[i]nvert_multiwin" >/dev/null 2>&1; do sleep 60; done
echo "=== baseline done; running CONTINUOUS-ONLY inversion ===" | tee -a $L

CONTINUOUS_ONLY=1 SFX=_calT_des DESEASON=1 \
  OUT=fault_tomography/inversion/fault_4d_multiwin_calT35_cont.npz \
  $PY fault_tomography/inversion/invert_multiwin.py >> $L 2>&1

echo "=== A/B COMPARISON (all-families vs continuous-only, true-scale 35-sta) ===" | tee -a $L
$PY - <<'EOF' 2>&1 | tee -a $L
import numpy as np
def stats(f,tag):
    try: d=np.load(f,allow_pickle=True)
    except Exception as e: print(tag,'MISSING',e); return
    ok=d['ok'].astype(bool); well=d['well'].astype(bool); MF=d['MF']
    perc=np.nanstd(MF[:,ok],axis=1)
    fidx=np.nanmean(MF[np.ix_(well,ok)],axis=0)
    vr=float(np.nanmean(d['VR'][ok])) if 'VR' in d.files else float('nan')
    print(f'{tag:16s}: resolved cells {int(well.sum())} | median per-cell RMS {100*np.nanmedian(perc[well]):.4f}% '
          f'| network-index std {100*np.nanstd(fidx):.4f}% | mean VR {vr:.0f}%')
stats('fault_tomography/inversion/fault_4d_multiwin_calT35.npz','ALL families')
stats('fault_tomography/inversion/fault_4d_multiwin_calT35_cont.npz','CONTINUOUS only')
EOF
echo "=== AB DONE ===" | tee -a $L
