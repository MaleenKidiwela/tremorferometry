---
name: gpu-family-discovery
description: discover_gpu.py is the fast/memory-safe family-discovery stage-2; obspy resample is the recurring OOM culprit
metadata: 
  node_type: memory
  type: project
  originSessionId: 5cecf374-37e4-4bcc-aea9-79b2b6cdc569
---

`scripts/discover_gpu.py` is the production stage-2 (clustering) for NLLB-style family
discovery — built this session, validated bit-equivalent to the CPU path. Three wins:
load-each-day-once (group cuts by day), `resample_poly` not obspy FFT-resample, and GPU
all-pairs CC (`src/tremorferometry/repeater_gpu.py`, 197× faster, ZERO cluster-flips at
cc≥0.8). HDW: 29,721 families in 89 s at anon 2 G (CPU was ~100 min / 54 G).

**Why:** the CPU discovery OOM'd (anon→117 G) and was slow because obspy's FFT `tr.resample`
on 100 Hz day-files × many workers blows memory, and the N²-CC-per-worker piled on. The fix
is single-GPU-process + CPU workers that only cut waveforms (KB each).

**How to apply:** (1) stage-1 via `discover_nllb_pnsn_driven.py --candidates-only` (now uses
`_load_day_poly`/resample_poly — without it 24-worker stage-1 hit anon 99 G). (2) GPU stage-2
via discover_gpu.py. The obspy-resample → `resample_poly` swap is the recurring fix for BOTH
memory and speed across this codebase. discover_gpu.py writes the **snr column** the committed
discover script silently dropped (env-peak / pre-pulse-RMS) — needed by select_coverage_families /
plot_family_map. Always run a memory watchdog (kill if anon>120 G). See [[pod-cpu-cap]],
[[pnsn-catalog-cut]], [[nisqually-precursor-project]].
