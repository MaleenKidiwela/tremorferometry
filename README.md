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

```bash
pip install -e .[dev]
# Optional GPU template matching (requires CUDA toolchain):
pip install -e .[gpu]
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
