#!/home/jovyan/envs/tremorferometry/bin/python
"""Run B933 stage-1 PNSN candidate discovery."""
import os, sys, subprocess

os.environ['PYTHONPATH'] = 'src'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

PYBIN = sys.executable  # Use this script's own interpreter

cmd = [PYBIN, 'scripts/discover_nllb_pnsn_driven.py',
       '--station', 'B933',
       '--pnsn', 'catalogs/pnsn_tremor_cascadia_full.csv',
       '--bbox', '39.1561', '40.9581', '-125.1490', '-122.7966',
       '--candidates-out', 'data/b933_pnsn_candidates_100km.parquet',
       '--candidates-only',
       '--workers', '16']

print(f"[run] {' '.join(cmd)}", flush=True)
r = subprocess.run(cmd, env={**os.environ})
print(f"[exit] {r.returncode}", flush=True)
sys.exit(r.returncode)
