#!/usr/bin/env python3
"""Launcher for densify_gnw_gpu.py that sets CUDA_PATH before any cupy import."""
import os
import sys

# Must set CUDA_PATH BEFORE cupy is imported anywhere
os.environ['CUDA_PATH'] = '/opt/conda/targets/x86_64-linux'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['PYTHONPATH'] = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'src')

# Add src to path
sys.path.insert(0, os.environ['PYTHONPATH'])

# Now run densify_gnw_gpu main
import runpy
sys.argv[0] = os.path.join(os.path.dirname(__file__), 'densify_gnw_gpu.py')
runpy.run_path(sys.argv[0], run_name='__main__')
