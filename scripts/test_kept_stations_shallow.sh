#!/bin/bash
# Kept-station test program (B941, B927): compute the template-free ambient-ACF shallow monitor on each,
# report its coherence (is the shallow signal real here?), then run the T2a site-filter (family dv/v vs ACF)
# to extend the B018 depth/shallow split to 2 more live stations. Non-destructive (waveforms retained).
set -u
cd /home/jovyan/tremorferometry
export PYTHONPATH=src OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
PY=/home/jovyan/envs/tremorferometry/bin/python
L=logs/test_kept_shallow.log; : > $L
for STA in B941 B927; do
  s=$(echo $STA|tr A-Z a-z)
  echo "=== $(date +%H:%M) $STA: ambient-ACF (template-free shallow monitor) ===" | tee -a $L
  if [ ! -f data/autocorr_dvv_${STA}.csv ]; then
    $PY scripts/autocorr_dvv.py --net PB --sta $STA --y0 2005 --y1 2026 --workers 24 \
        --out data/autocorr_dvv_${STA}.csv >> $L 2>&1
  fi
  echo "=== $(date +%H:%M) $STA: ACF quality + T2a site-filter ===" | tee -a $L
  $PY - "$STA" <<'EOF' >> $L 2>&1
import pandas as pd, numpy as np, sys
sta=sys.argv[1]
a=pd.read_csv(f'data/autocorr_dvv_{sta}.csv')
ok = a.cc_max.median()>=0.6
print(f"[{sta}] ACF monitor: n={len(a)} cc_max median {a.cc_max.median():.2f} mean {a.cc_max.mean():.2f} "
      f"| cc>0.6 {(a.cc_max>0.6).mean()*100:.0f}% | dvv std {a.dvv.std()*100:.2f}% -> "
      f"{'USABLE shallow monitor' if ok else 'WEAK/dead monitor (T2a will be uninformative, like B935)'}")
EOF
  $PY scripts/exp_t2a_fixed.py $STA >> $L 2>&1
  grep -E "ACF monitor|class x tier|GOLD summary" $L | tail -4
done
echo "=== $(date +%H:%M) KEPT-STATION TEST DONE ===" | tee -a $L