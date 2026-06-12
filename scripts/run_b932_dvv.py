#!/home/jovyan/envs/tremorferometry/bin/python
"""Run B932 all-time dv/v then calendar rolling windows."""
import os, sys, subprocess

os.environ['PYTHONPATH'] = 'src'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

PYBIN = '/home/jovyan/envs/tremorferometry/bin/python'

steps = [
    # All-time 1-4s dv/v
    [PYBIN, 'scripts/dvv_coda_parallel.py',
     '--npz', 'data/long_window_daily_B932.npz',
     '--window', '1.0', '4.0',
     '--station', 'B932',
     '--out-csv', 'data/daily_dvv_B932_coda_1to4.csv',
     '--out-fig', 'figures/smoke_dvv_B932_raw.png',
     '--workers', '20'],
    # Calendar rolling 1-3
    [PYBIN, 'scripts/dvv_roll30cal.py',
     '--station', 'B932',
     '--npz', 'data/long_window_daily_B932.npz',
     '--window', '1.0', '3.0',
     '--out', 'data/daily_dvv_B932_1to3_cal.csv',
     '--workers', '16'],
    # Calendar rolling 2-4
    [PYBIN, 'scripts/dvv_roll30cal.py',
     '--station', 'B932',
     '--npz', 'data/long_window_daily_B932.npz',
     '--window', '2.0', '4.0',
     '--out', 'data/daily_dvv_B932_2to4_cal.csv',
     '--workers', '16'],
    # Calendar rolling 3-5
    [PYBIN, 'scripts/dvv_roll30cal.py',
     '--station', 'B932',
     '--npz', 'data/long_window_daily_B932.npz',
     '--window', '3.0', '5.0',
     '--out', 'data/daily_dvv_B932_3to5_cal.csv',
     '--workers', '16'],
]

for cmd in steps:
    print(f"[run] {' '.join(cmd)}", flush=True)
    r = subprocess.run(cmd, env={**os.environ})
    if r.returncode != 0:
        print(f"[FAIL] exit {r.returncode}", flush=True)
        sys.exit(r.returncode)
    print(f"[ok]", flush=True)

print("[B932 dvv+rolling DONE]")
