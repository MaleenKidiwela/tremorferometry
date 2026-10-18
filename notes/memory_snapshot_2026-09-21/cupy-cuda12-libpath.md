---
name: cupy-cuda12-libpath
description: GPU steps (discover_gpu/densify) fail with libcufft.so.11 missing — fix by prepending the base-conda CUDA-12 nvidia libs to LD_LIBRARY_PATH
metadata: 
  node_type: memory
  type: reference
  originSessionId: bfadfcae-1eeb-458c-a71d-45a42210b672
---

cupy in the env (`/home/jovyan/envs/tremorferometry/bin/python`, cupy 13.6.0 = the **cuda12x** build)
needs the **CUDA-12** runtime libs (`libcufft.so.11`, `libcublas.so.12`, `libcusolver.so.11`, …), but the
env itself only ships **CUDA-13** libs (`nvidia/cu13/lib/libcufft.so.12`) — so any GPU step dies with
`ImportError: libcufft.so.11: cannot open shared object file`. This recurs after pod restarts (see
[[pod-env-ephemeral]]).

**Fix (export before every GPU job — discover_gpu.py, densify):**
```
NVLIBS=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')
export LD_LIBRARY_PATH="${NVLIBS}${LD_LIBRARY_PATH}"
export CUDA_PATH=/opt/conda/targets/x86_64-linux
```
The base-conda py3.13 site-packages carries a COMPLETE CUDA-12 nvidia set (cublas, cufft, cusolver,
cusparse, curand, cuda_nvrtc, …) even though the env's own python is 3.11 — that's the donor. Verify with a
cupy FFT smoke test (`cp.fft.rfft`) before launching. `CUDA_PATH` alone is NOT enough (that was the old note);
the loader path to the CUDA-12 `.so`s is the actual missing piece. See also [[pod-cpu-cap]] (workers ≤30).
