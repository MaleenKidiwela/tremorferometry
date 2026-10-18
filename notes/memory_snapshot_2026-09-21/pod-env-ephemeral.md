---
name: pod-env-ephemeral
description: "Pod env after rebuild — use envs/tremorferometry python (not base), and set CUDA_PATH so cupy JIT finds headers"
metadata: 
  node_type: memory
  type: reference
  originSessionId: 400c26f9-17ff-46fc-914c-9c6a7671906d
---

This Jupyter pod's `/opt/conda` lives on an **ephemeral overlay** filesystem — a pod crash/reschedule resets it to the base image. Only `/home/jovyan` (NFS PVC) persists.

**As of 2026-06-05 the pod was REBUILT and the env layout changed:**
- **USE the project env: `/home/jovyan/envs/tremorferometry/bin/python` (Python 3.11.15).** It lives on persistent `/home/jovyan/envs`, so it survives crashes. It has the full stack: obspy 1.5.0, cupy-cuda12x 13.6.0, pyarrow, cartopy 0.25, h5py, scipy, pandas, **numpy 2.4.6** (numba NOT installed here — that's why numpy 2 is fine; the old "numpy<2 / cupy<14" constraint does NOT apply to this env).
- **Do NOT use base conda `/opt/conda/bin/python`** — it's now Python 3.13 with numpy 2 and **no obspy/cupy/pyarrow/cartopy**; every pipeline import fails.
- There's also `/home/jovyan/envs/pnw` (unrelated).

**cupy CUDA-header gotcha (cost an hour 2026-06-05):** after a rebuild, cupy's NVRTC JIT fails on the first GPU kernel with `NVRTCError: NVRTC_ERROR_COMPILATION ... cannot open source file "cuda_fp16.h"`. FIX: **`export CUDA_PATH=/opt/conda/targets/x86_64-linux`** (a CUDA 12.4 cudart-dev include dir with cuda_fp16.h, matching cupy-cuda12x). Set it for every GPU job (densify, discover_gpu). Verified: `cp.fft.rfft` then compiles.

**Standard invocation for any pipeline script:**
```
export PYTHONPATH=src OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export CUDA_PATH=/opt/conda/targets/x86_64-linux   # GPU jobs only
/home/jovyan/envs/tremorferometry/bin/python scripts/<...>.py ...
```
GPU now reports 4 devices (was 1); 47 GB free. Memory cap ~187 GB. See [[pod-cpu-cap]], [[crash-recovery-context]], [[nisqually-precursor-project]].

**UPDATE 2026-09-10 (crash at ~21:03):** the current base image `/opt/conda/bin/python` (3.13) DOES ship obspy 1.5.0, cupy-cuda12x, pyarrow, the nvidia cu12 libs — the per-station fleet used it. What the crash wiped was everything pip-installed into `/opt/conda` afterwards: **the CUDA torch build** (e2 picker ran on `cuda` pre-crash; post-crash only `torch 2.12.0+cpu` in `~/.local` remained → `torch.cuda.is_available()` False, LFEDetect falls to CPU). `pip install --user torch --index-url .../cu128` is a NO-OP ("already satisfied" by +cpu) — must `pip uninstall -y torch` then `pip install --user torch==2.12.0+cu130 --index-url https://download.pytorch.org/whl/cu130 --extra-index-url https://pypi.org/simple` (driver 580 / CUDA 13.0; cu128 index has no 2.12). seisbench 0.12.5, pyocto 0.2.0, NonLinLoc (`~/.local/bin`) are in `~/.local` and survive. Rule: **anything GPU/DL must be `pip install --user`**; after every crash run the import check in `lfe_catalog_repair/scripts/e2_run_pick.sh`'s env before launching GPU queues.
