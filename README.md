# tremorferometry

LFE coda-wave interferometry for dv/v across Cascadia ETS events.

## Idea

Use Bostock-style Low-Frequency Earthquake (LFE) families as **repeating
sources** for coda-wave interferometry, and track temporal velocity changes
(dv/v) at the plate interface across ETS cycles. LFEs sit on the megathrust /
transition zone, so the dv/v sensitivity kernel is biased toward the slipping
patch — distinct from ambient-noise dv/v, which is dominated by shallow crust.

Catalogs used:
- **PNSN tremor catalog** (Wech envelope-correlation; pnsn.org/tremor) — defines
  ETS episode timing and the slipping bbox.
- **Bostock LFE catalog (2012/2015)** — repeating LFE families with template
  waveforms. The catalog template itself is the "source" for CWI.

## Pipeline

```
01 PNSN tremor    -> episode (t_start, t_end, bbox)
02 Bostock raw    -> normalized family table + template MSEED
03 catalog + bbox -> selected families
04 FDSN           -> continuous waveforms (parallel ThreadPool)
05 templates + data -> EQcorrscan + fast-matched-filter (GPU) -> detections
06 detections     -> time-binned LFE-aligned stacks (HDF5)
07 stacks + ref   -> dv/v(t) per family/station via stretching (parallel joblib)
08 dv/v parquet   -> figures
```

Scripts under `scripts/` are numbered and idempotent; each takes `--config`.

## Repo layout

```
configs/        # one YAML per ETS episode
catalogs/       # ingested catalog tables
scripts/        # numbered thin CLIs
src/tremorferometry/  # the package
tests/          # pytest
data/           # gitignored; waveforms, templates, detections, stacks, dvv
figures/        # output plots
```

## Install

A dedicated conda env at `/home/jovyan/envs/tremorferometry` already has the
full stack installed (obspy, eqcorrscan, fast-matched-filter built against
CUDA 12.4, plus this package in editable mode). Activate with

```bash
mamba activate /home/jovyan/envs/tremorferometry
```

To recreate from scratch on another host:

```bash
mamba create -p /path/to/env -c conda-forge python=3.11 \
  numpy scipy pandas matplotlib h5py pyarrow pyyaml \
  obspy eqcorrscan joblib tqdm pytest
mamba activate /path/to/env

# fast-matched-filter (GPU template matching, needs nvcc on PATH)
pip install --no-build-isolation \
  "git+https://github.com/beridel/fast_matched_filter.git"
# FMF ships sources only; build the CPU + GPU .so files:
FMF=$(python -c "import fast_matched_filter, os; print(os.path.dirname(fast_matched_filter.__file__))")
mkdir -p $FMF/lib && cd $FMF/src
gcc -O3 -fopenmp -fPIC -march=native -shared matched_filter.c \
    -o $FMF/lib/matched_filter_CPU.so
nvcc -O3 -Xcompiler "-fPIC -fopenmp" -shared matched_filter.cu \
     -o $FMF/lib/matched_filter_GPU.so

pip install --no-deps --no-build-isolation -e .
```

## v1 target

One well-recorded southern Vancouver Island ETS (e.g. 2010-08), ~5 PNSN
broadbands, a handful of Bostock LFE families inside the slipping patch,
+/- 90 d window, 2-day bins, 2-8 Hz coda 5-25 s post-S. Output: a dv/v(t)
figure and a parquet of values, with QC checks passing.

## Status

Scaffold + core dv/v measurement + tests in place. External-data ingestion
(PNSN CSV, Bostock supplement) requires manual file placement; see the
docstrings in scripts 01 and 02.
