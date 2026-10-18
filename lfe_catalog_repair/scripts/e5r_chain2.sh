#!/bin/bash
# E5 step 5: after the pass-2 re-run -> score both runs (SNR-weighted writes the catalog) -> Thomas harness on the final catalogs -> figures
cd /home/jovyan/tremorferometry/lfe_catalog_repair
export PYTHONPATH=/opt/conda/lib/python3.13/site-packages OMP_NUM_THREADS=8
until grep -q "STEP D DONE" logs/e5r_run_mf2.log || grep -lq "Traceback" logs/e5r_mf_snr.log logs/e5r_mf_eq.log 2>/dev/null; do sleep 15; done
if ! grep -q "STEP D DONE" logs/e5r_run_mf2.log; then echo "MF re-run FAILED"; exit 1; fi
echo "[$(date -u +%FT%TZ)] scoring"
python scripts/e5r_score.py --tag e5r_snr --write-catalog > logs/e5r_score_snr.log 2>&1 &
python scripts/e5r_score.py --tag e5r_eq > logs/e5r_score_eq.log 2>&1 &
wait
echo "[$(date -u +%FT%TZ)] thomas"
bash scripts/e5r_run_thomas.sh > logs/e5r_run_thomas.log 2>&1
echo "[$(date -u +%FT%TZ)] figures"
python scripts/e5r_figures.py > logs/e5r_figures.log 2>&1
echo "[$(date -u +%FT%TZ)] chain2 DONE"
