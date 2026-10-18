---
name: pod-cpu-cap
description: This pod is CPU-capped at 32 cores via cgroup; os.cpu_count()/nproc lie (report 176)
metadata: 
  node_type: memory
  type: reference
  originSessionId: be83f309-c83d-4336-b8a2-99a689098c27
  modified: 2026-08-22T07:51:09.033Z
---

The Jupyter pod's **cgroup CPU quota is 32 CPUs**, far below the node's 176 cores.
Check: `cat /sys/fs/cgroup/cpu.max` → `3200000 100000` (quota/period = 32).

**`os.cpu_count()`, `nproc`, `len(os.sched_getaffinity(0))` all report 176** — they see the
*node*, not the pod's CFS quota. So code that auto-sizes a worker pool from `cpu_count()`
will massively oversubscribe.

**Why it matters:** a ProcessPool with `--workers 48` (or anything >32) thrashes the 32-CPU
quota → CFS throttling + context-switch overhead → ~2.4× effective speedup instead of ~30×.
The GNW densify projected ~13 h for this reason; **`--workers 30` (≤32) is the single biggest
speedup.** Also pin per-process threads so 30 workers don't each spawn FFT/BLAS threads:
`export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`
before importing numpy.

Memory cap `memory.max` ≈ 187 GB. See [[nisqually-precursor-project]], [[pod-env-ephemeral]].

2026-08-22 additions: unpinned BLAS doesn't just slow the job — it makes the whole pod so
laggy the user restarts it (looks like a crash; memory/disk were fine). Prefer setting the
thread env vars IN-SCRIPT via `os.environ.setdefault` before `import numpy` so no launch can
forget them. For CPU python jobs, launch under `scripts/run_capped.sh [GB] <cmd>` (RLIMIT_AS,
default 140 GB) — the only OOM protection fast enough for a runaway allocation (a 212 GB
einsum OOM-killed the pod; watchdogs poll too slowly). NOT for GPU/cupy jobs (AS cap breaks
CUDA); densify keeps using oom_guard.sh. Leave headroom: `--workers 24` keeps typing usable.
