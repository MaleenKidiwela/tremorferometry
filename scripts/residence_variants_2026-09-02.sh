#!/bin/bash
# Merlin Step-0 mandatory checks (2026-09-02 reply): seed reproducibility, HG forward-scattering bracket, slab
# thickness bracket, voxel convergence. Subsample = PAIRS=200 (RandomState(0) choice, same pairs every run).
cd /home/jovyan/tremorferometry
M="bash scripts/run_capped.sh 100 python fault_tomography/kernels/residence_kernel_mc.py"
PAIRS=200 WORKERS=4 OUTTAG=_sub200_seed1 $M
PAIRS=200 WORKERS=4 OUTTAG=_sub200_seed2 SEED_OFFSET=7777 $M
PAIRS=200 WORKERS=4 OUTTAG=_sub200_g05 HG_G=0.5 $M
PAIRS=200 WORKERS=4 OUTTAG=_sub200_g08 HG_G=0.8 $M
PAIRS=200 WORKERS=4 OUTTAG=_sub200_H3 HSLAB=3.0 $M
PAIRS=200 WORKERS=4 OUTTAG=_sub200_H5 HSLAB=5.0 $M
PAIRS=100 WORKERS=4 OUTTAG=_sub100_fine DX=2.5 DZ=0.5 DS=0.5 $M
PAIRS=100 WORKERS=4 OUTTAG=_sub100_coarse $M
echo VARIANTS_DONE
