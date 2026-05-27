# tremorferometry

LFE coda-wave interferometry for dv/v across Cascadia ETS events.

## Idea

Use Low-Frequency Earthquake (LFE) families as **repeating sources** for
coda-wave interferometry, and track temporal velocity changes (dv/v) at the
plate interface across ETS cycles. LFEs sit on the megathrust / transition
zone, so the dv/v sensitivity kernel is biased toward the slipping patch —
distinct from ambient-noise dv/v, which is dominated by shallow crust.

## v1 results (S. V.I., 2010-2013 ETS cycles)

**Two independent methods, consistent null.**

**Method A — tremor-windowed inter-station CC** (Phase E.6):
    dv/v = +0.0003 +/- 0.0075 %  (within the 2010-08 ETS)
    2-sigma upper bound: |dv/v| < 0.015 %

**Method B — repeating LFE coda-wave interferometry** (Phase E.7e-g, the
breakthrough): proper Shelly-Beroza-style network autocorrelation discovered
**16 cross-year LFE families** at southern V.I. (CC ≥ 0.7 at PGC + LZB),
with members spanning **2005–2013** (up to 8-year baselines). Multi-seed
matched-filter then grew the catalog to **~17,000 LFE detections across 4 ETS
cycles**. Stretching per-family year-stacks against 2010 reference:

    2011 vs 2010 ETS:  dv/v = -0.054 +/- 0.291 %
    2012 vs 2010 ETS:  dv/v = +0.033 +/- 0.406 %
    2013 vs 2010 ETS:  dv/v = -0.246 +/- 0.297 %

All consistent with zero within 1σ. Per-family CC matches are 0.93+ — the
templates are real, the measurement is between actual same-source events.

**Conclusion:** the southern V.I. plate interface does not show a
detectable dv/v signal at ~0.3 % across 2010-2013 ETS cycles. Either the
real change is below that floor, or it's spatially localized below our
per-family resolution. Two methods agree.

### The path (LFE-CWI eventually worked — with the right methodology)

The first LFE-CWI attempt (using Lin's catalog binned at 0.1°) **failed**
because at that clustering grain, detections within a "family" are NOT
waveform-repeating sources. The early Phase A–D numbers were stretching
on source-mixture variation, not medium change.

The pivot to **tremor-windowed inter-station CC** (Method A above) gave
the first clean null bound — does not require repeating sources, biases
sensitivity to the megathrust depth via the tremor source field.

Then the proper LFE-CWI was rescued (Method B) with the Shelly-Beroza
recipe in `src/tremorferometry/repeater.py`:
1. Envelope-peak alignment (don't trust Lin's OT — find the direct phase).
2. Tight 2-sec window (not 6 s of band-limited noise).
3. Network CC at PGC + LZB with shift allowance.
4. Strict threshold CC ≥ 0.7.
5. Hotspot focus.

Result: 16 cross-year families found in seed analysis; matched-filter then
grew them to ~17,000 detections. See `notes/METHODS.md` § Phase E for the
full arc.

### Margin-wide and 2026 extension

The pipeline is region-agnostic. To extend everywhere in Cascadia / to 2026:

- **Southern V.I.** (have): use Lin (2023) catalog.
- **Olympic / SW Washington**: Sweet et al. (2019) catalog (9 LFE families).
- **N California**: Plourde (2015) and Ducellier (2022).
- **Central Oregon gap** and **2014-onward** beyond Lin's coverage: use
  `src/tremorferometry/detect.py` (envelope peaks inside PNSN tremor
  windows) to seed `repeater.py` on continuous data. PNSN tremor catalog
  covers all of Cascadia continuously through 2026.

### Figures

- `figures/smoke_tremor_cc_dvv_v2.png` — main result (aggregate dv/v vs time with LFE-rate panel)
- `figures/smoke_tremor_cc_per_pair_distance.png` — per-pair robustness
- `figures/smoke_tremor_cc_preprocess_compare.png` — preprocessing tradeoff
- `figures/smoke_family_similarity.png`, `figures/smoke_reclust_L0000_PGC.png` —
  diagnostics that walked back the LFE-CWI approach

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
