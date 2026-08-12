#!/bin/bash
# B003 EXPANSION PILOT (the gate before the fleet campaign).
# Densify the 346 genuine Stage-1 candidates, merge with the existing 81, restack, re-dv/v (true-scale),
# rerun the 35-station inversion with B003 EXPANDED (others fixed), and measure the resolved-cell gain.
cd /home/jovyan/tremorferometry
export PYTHONPATH=src OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export CUDA_PATH=/opt/conda/targets/x86_64-linux
PY=/home/jovyan/envs/tremorferometry/bin/python
L=logs/expand_b003.log; : > $L
say(){ echo "=== $(date +%H:%M:%S) $* ===" | tee -a $L; }

say "1. densify 346 genuine candidates (GPU)"
if [ $(ls data/mf_b003_exp_[12]*.csv 2>/dev/null|wc -l) -eq 0 ]; then
  $PY scripts/densify_launcher.py --templates-npz data/b003_pnsn_families_100km.npz \
    --summary-csv data/b003_expansion_candidates.summary.csv --min-snr 0 --network PB --station B003 \
    --out-prefix mf_b003_exp_ --workers 20 --top-n 100 --max-raw-det 3000000 --despike-mad 8 >> $L 2>&1
fi

say "2. merge existing 81 + new 346 -> expanded detection set"
$PY - <<'EOF' 2>&1 | tee -a $L
import pandas as pd, glob
base=pd.read_csv('data/mf_b003_all.csv')
exp=pd.concat([pd.read_csv(f) for f in sorted(glob.glob('data/mf_b003_exp_[12]*.csv'))],ignore_index=True)
alld=pd.concat([base,exp],ignore_index=True).drop_duplicates(['template','time'])
alld.to_csv('data/mf_b003_expanded_all.csv',index=False)
print(f'base {base.template.nunique()} fam, exp {exp.template.nunique()} fam, merged {alld.template.nunique()} fam, {len(alld)} det')
EOF

say "3. restack expanded (min-det 5, sparse station)"
$PY scripts/build_long_window_resp.py --mf-csv data/mf_b003_expanded_all.csv --network PB --station B003 \
  --no-deconv --min-det 5 --despike-mad 8 --workers 16 --out data/long_window_daily_B003_expanded.npz >> $L 2>&1

say "4. re-dv/v true-scale (3 windows) + deseason"
for W in "1.0 3.0 1to3" "2.0 4.0 2to4" "3.0 5.0 3to5"; do set -- $W
  $PY scripts/dvv_roll30cal.py --station B003 --npz data/long_window_daily_B003_expanded.npz --window $1 $2 \
    --origin-anchor --out data/daily_dvv_B003_${3}_calT_exp.csv --workers 16 >> $L 2>&1
  $PY scripts/deseason_cal.py --glob "data/daily_dvv_B003_${3}_calT_exp.csv" >> $L 2>&1
done

say "5. rerun 35-sta inversion with B003 EXPANDED (others fixed)"
EXP_STA=B003 EXP_SFX=_calT_exp_des SFX=_calT_des DESEASON=1 \
  OUT=fault_tomography/inversion/fault_4d_multiwin_calT35_b003exp.npz \
  $PY fault_tomography/inversion/invert_multiwin.py >> $L 2>&1

say "6. A/B: B003 baseline (81 fam) vs expanded"
$PY - <<'EOF' 2>&1 | tee -a $L
import numpy as np, pandas as pd
# B003 patch->cell coverage before/after
def cells(csv):
    d=pd.read_csv(csv,usecols=['patch']);
    return set((round(round(float(p.split('__')[0].split('_')[0])/0.1)*0.1,3),
                round(round(float(p.split('__')[0].split('_')[1])/0.1)*0.1,3)) for p in d.patch.unique())
cb=cells('data/daily_dvv_B003_2to4_calT_des.csv'); ce=cells('data/daily_dvv_B003_2to4_calT_exp_des.csv')
print(f'B003 distinct 0.1-deg cells imaged: baseline {len(cb)} -> expanded {len(ce)} (+{len(ce)-len(cb)})')
def res(f):
    d=np.load(f,allow_pickle=True); return int(d['well'].astype(bool).sum()), \
        100*np.nanmedian(np.nanstd(d['MF'][:,d['ok'].astype(bool)],axis=1)[d['well'].astype(bool)])
rb,mb=res('fault_tomography/inversion/fault_4d_multiwin_calT35.npz')
re_,me=res('fault_tomography/inversion/fault_4d_multiwin_calT35_b003exp.npz')
print(f'WHOLE-NETWORK resolved cells: baseline {rb} -> B003-expanded {re_} (+{re_-rb})')
print(f'median per-cell RMS: baseline {mb:.4f} -> {me:.4f} (note already %)')
EOF
say "EXPAND_B003 PILOT DONE"
