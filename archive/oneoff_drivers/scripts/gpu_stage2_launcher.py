#!/usr/bin/env python3
"""Launcher for discover_gpu.py that sets CUDA_PATH before any cupy import."""
import os
import sys

# Must set CUDA_PATH BEFORE cupy is imported anywhere
os.environ['CUDA_PATH'] = '/opt/conda/targets/x86_64-linux'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'src'))

# Now run discover_gpu main
import runpy
sys.argv[0] = os.path.join(os.path.dirname(__file__), 'discover_gpu.py')
runpy.run_path(sys.argv[0], run_name='__main__')
