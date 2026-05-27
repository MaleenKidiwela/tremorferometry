# tremorferometry

LFE coda-wave interferometry for dv/v across Cascadia ETS events.

## Idea

Use Low-Frequency Earthquake (LFE) families as **repeating sources** for
coda-wave interferometry, and track temporal velocity changes (dv/v) at the
plate interface across ETS cycles. LFEs sit on the megathrust / transition
zone, so the dv/v sensitivity kernel is biased toward the slipping patch —
distinct from ambient-noise dv/v, which is dominated by shallow crust.

## v1 result (2010-08 V.I. ETS)

![Phase D dv/v figure](figures/smoke_phaseD_dvv.png)

- **217 LFE families × 6 V.I./Olympic broadbands × 76 two-day bins → 5,948 dv/v measurements**
- ETS-window mean dv/v = **+0.04 ± 0.05 %** (1σ standard error on the mean)
- Gaussian-smoothed envelope stays within ±0.2 % through the ETS
- Spatial coherence QC: 89 % of bins consistent across stations
- Coda-window stability QC: three different coda windows all return ETS means consistent with zero
- **Headline: no detectable sustained ETS-wide dv/v signature above ~0.1 %**
  in the 2010-08 V.I. ETS at our station coverage.

See `notes/METHODS.md` for the full methods, results, and limitations.

## Catalogs

- **PNSN tremor catalog** (Wech envelope-correlation method;
  `https://tremorapi.pnsn.org/api/v3.0/events`) — used to identify and time
  ETS episodes. Implementation: `src/tremorferometry/pnsn.py` and
  `scripts/00b_fetch_pnsn_tremor.py`.
- **Lin (2023) LFE catalog** (Zenodo `10.5281/ZENODO.10016020`) — 1.05 M LFE
  detection times + locations for southern Vancouver Island 2005–2017.
  Ingest: `src/tremorferometry/lin_catalog.py` and
  `scripts/02b_ingest_lin_catalog.py`. We cluster Lin's per-event records
  into 217 "proto-families" by 0.1° grid binning (limitation: see
  per-family analysis section in METHODS).

## Pipeline

```
00b PNSN API         -> tremor CSV
01  tremor CSV       -> episode (t_start, t_end, bbox)
02b Lin Zenodo CSV   -> families + detection-times parquet
04  FDSN             -> per-day MSEED waveform cache (parallel ThreadPool)
06  detections + MSEED -> per-family HDF5 stacks
                          (flat (family, station) ProcessPool)
07  stacks + ref     -> dv/v(t) parquet via stretching
09  dv/v parquet     -> aggregate figure + QC

Optional / for time periods outside Lin's coverage:
05  templates + data -> EQcorrscan + fast-matched-filter (GPU) detections
```

Scripts under `scripts/` are numbered and idempotent; each takes `--config`.

## Repo layout

```
configs/        # one YAML per ETS episode
catalogs/       # ingested catalog tables (gitignored, reproducible)
scripts/        # numbered thin CLIs
src/tremorferometry/  # the package
tests/          # pytest, 12 tests passing
notes/          # METHODS.md and other write-ups
data/           # gitignored; waveforms, detections, stacks, dvv
figures/        # tracked smoke + result figures
```

## Install

A dedicated conda env at `/home/jovyan/envs/tremorferometry` has the full
stack installed (obspy, eqcorrscan, fast-matched-filter built against
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

# fast-matched-filter (GPU template matching, needs nvcc on PATH; optional
# for v1 since we use Lin's detection times directly).
pip install --no-build-isolation \
  "git+https://github.com/beridel/fast_matched_filter.git"
FMF=$(python -c "import fast_matched_filter, os; print(os.path.dirname(fast_matched_filter.__file__))")
mkdir -p $FMF/lib && cd $FMF/src
gcc -O3 -fopenmp -fPIC -march=native -shared matched_filter.c \
    -o $FMF/lib/matched_filter_CPU.so
nvcc -O3 -Xcompiler "-fPIC -fopenmp" -shared matched_filter.cu \
     -o $FMF/lib/matched_filter_GPU.so

pip install --no-deps --no-build-isolation -e .
```

## Reproduce v1 end-to-end

```bash
mkdir -p data/raw_lfe
curl -fsSL -o data/raw_lfe/lin2023_lfe.csv \
    'https://zenodo.org/records/10016020/files/EQloc_001_0.1_3_S.csv?download=1'

python scripts/02b_ingest_lin_catalog.py --config configs/ets_2010_vi.yaml \
    --raw data/raw_lfe/lin2023_lfe.csv

python scripts/04_fetch_waveforms.py --config configs/ets_2010_vi.yaml --workers 24

python scripts/06_stack_bins.py --config configs/ets_2010_vi.yaml \
    --detections data/detections_lin_ets_2010_vi.parquet \
    --selected catalogs/lin_families_ets_2010_vi.csv \
    --waveforms data/waveforms --out data/stacks --workers 48 \
    --fs 40 --pre-s 5 --post-s 35

python scripts/07_measure_dvv.py --config configs/ets_2010_vi.yaml \
    --stacks data/stacks --out data/dvv/v1.parquet --workers 16 --fs 40

python scripts/09_aggregate_plot.py --config configs/ets_2010_vi.yaml \
    --dvv data/dvv/v1.parquet \
    --detections data/detections_lin_ets_2010_vi.parquet \
    --out figures/v1.png
```

## Status

v1 complete: end-to-end pipeline on real 2010-08 V.I. ETS data, multi-station,
multi-family, with QC checks. The remaining open work (waveform-similarity
family reclustering, frequency-band split, margin-wide extension to other LFE
catalogs) is captured in `notes/METHODS.md` § 6.
