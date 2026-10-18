#!/bin/bash
cd /home/jovyan/tremorferometry; export RESDIR=fault_tomography/inversion/splits/Gdiff24 WORKERS=20
bash scripts/run_capped.sh 200 python fault_tomography/inversion/invert_dvv_4d_daily.py --lam 0 --workers 20 --max-sweeps 30 > logs/gdiff24_solver.log 2>&1 && \
bash scripts/run_capped.sh 200 python fault_tomography/inversion/transmission_T.py > logs/gdiff24_T.log 2>&1 && \
bash scripts/run_capped.sh 200 python fault_tomography/inversion/transmission_deconv.py > logs/gdiff24_deconv.log 2>&1
echo "GDIFF_DONE rc=$?" >> logs/gdiff24_solver.log
