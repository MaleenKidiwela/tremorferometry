#!/bin/bash
# Hard per-process memory cap — THE robust OOM protection (2026-08-22, after invert_dvv_4d_daily
# materialized a 212 GB einsum operand and OOM-killed the whole pod; watchdogs poll at 30-60 s and
# could not react to a 90-s blowup, and oom_guard's pause action only knows densify jobs).
#
# RLIMIT_AS makes a runaway allocation fail INSTANTLY with MemoryError inside the job (traceback in
# its own log) instead of taking the pod down. Cap well below the ~187 GB pod cgroup limit.
#
# Usage:  nohup bash scripts/run_capped.sh [GB] <command...> > log 2>&1 &
#         default cap 140 GB; child processes inherit the limit (each worker individually).
#
# CAVEAT: do NOT use for GPU/cupy jobs (CUDA maps huge virtual ranges; AS cap breaks it).
# For densify, keep using oom_guard.sh (SIGSTOP/CONT). This wrapper is for CPU python jobs
# (inversion, assembly, stacking, catalogs).
if [[ "$1" =~ ^[0-9]+$ ]]; then CAP_GB=$1; shift; else CAP_GB=140; fi
ulimit -v $((CAP_GB * 1024 * 1024))
exec "$@"
