#!/bin/bash
cd /home/jovyan/tremorferometry; export RESDIR=fault_tomography/inversion/coarse05x05_gated WORKERS=20
bash scripts/run_capped.sh 200 python fault_tomography/inversion/invert_dvv_4d_daily.py --lam 0.100 --workers 20 --max-sweeps 40 > logs/coarse05x05_gated_solver.log 2>&1; echo "COARSE_DONE rc=$?" >> logs/coarse05x05_gated_solver.log
