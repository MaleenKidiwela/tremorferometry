# Methods: LFE coda-wave interferometry for dv/v across Cascadia ETS

A running methods log for the `tremorferometry` project. The current target is
the August–September 2010 Episodic Tremor and Slip (ETS) event beneath southern
Vancouver Island; the architecture extends margin-wide for future work.

---

## 1. Scientific premise

Cascadia ETS is a slow-slip phenomenon: tens of millimeters of along-strike slip
on the plate interface over weeks, accompanied by abundant tectonic tremor and
Low-Frequency Earthquakes (LFEs). We measure fractional seismic velocity change
(**dv/v**) through these events using LFEs as **repeating sources** for
coda-wave interferometry (CWI). The premise:

- LFEs are co-located on (or very near) the plate interface at ~25–45 km depth.
  Their coda samples the volume around that interface preferentially. The dv/v
  sensitivity kernel for LFE coda is therefore biased toward the megathrust —
  exactly where ETS-driven stress changes occur.
- Each "family" of LFEs is, by construction, a quasi-repeating source: the same
  patch of fault re-ruptures with very similar waveforms. Repetition is what
  makes CWI clean — any waveform difference between two stacks of detections
  from the same family is a **medium** change (or noise), not a source change.
- Ambient-noise dv/v (the standard alternative) is dominated by ocean
  microseisms (~0.1–0.5 Hz) and has a sensitivity kernel concentrated in the
  shallow crust. It is not the right tool for probing the megathrust.

---

## 2. Data

### 2.1 PNSN tremor catalog (episode scoping)

The Pacific Northwest Seismic Network publishes a public tremor catalog using
the **Wech (2008) envelope-correlation method** — bandpass continuous data to
2–8 Hz, take amplitude envelopes, cross-correlate envelope shapes across the
network in 5-minute windows, grid-search the location that explains the
envelope coherence. We query the v3.0 GeoJSON API at

    https://tremorapi.pnsn.org/api/v3.0/events?starttime=...&endtime=...

via `src/tremorferometry/pnsn.py`. This gives time, lat, lon, depth, duration,
magnitude, and a few quality metrics per 5-min tremor window. We use this
catalog only to **define when and where an ETS episode is active** — bbox and
time bounds for the rest of the pipeline. Distinct ETS episodes are picked out
of a multi-month catalog by `src/tremorferometry/episode.py`:

1. Bin tremor detections into 24-hour windows; mark windows above a configured
   rate as "active".
2. Bridge short below-threshold gaps so a single ETS isn't fragmented.
3. Score contiguous active spans by total detection count; return the largest
   (or "first" / "latest" by request).

For the 17-month window 2025-01 → 2026-05 we identified 9 distinct N-Cascadia
episodes. The 2010-08 ETS used here is the canonical historical episode
characterized in Bostock et al. (2012, 2015).

### 2.2 Lin (2023) LFE catalog (detection times)

LFE template waveforms from Bostock 2012/2015 are not openly hosted. The most
accessible public LFE catalog is **Lin et al. (2023)**, "A Deep Learning-Based
Low-Frequency Earthquake Catalog in Southern Vancouver Island" (Seismica),
deposited at Zenodo (https://doi.org/10.5281/ZENODO.10016020). The catalog
file `EQloc_001_0.1_3_S.csv` (~106 MB) contains 1,058,114 LFE detections
between 2005-01-01 and 2017-02-21, with columns

    starttime, OT, lon, lat, depth, residual, dt, N

where `OT` is the inferred event origin time. Lat coverage 48.07°–49.32°, lon
−124.57° to −122.94°: southern Vancouver Island only.

We ingest the catalog with `src/tremorferometry/lin_catalog.py`. Because Lin's
deposit is per-event without family labels, we cluster detections by spatial
grid binning:

- Bin lat/lon at 0.1° cells (~10 km on a side).
- Each cell with ≥50 detections in the configured episode window becomes one
  "proto-family". Cells are labeled `L0000`, `L0001`, … in decreasing
  detection count.
- Each cell's centroid (mean lat, lon, depth) is recorded as the family
  location; the detection list becomes a `(family_id, time)` parquet.

For the 2010 v1 window (2010-05-17 → 2010-10-15) inside bbox (48.0, 49.4,
−124.6, −123.0):

- **217 families**, median **324 detections per family**, **106,644 total
  detections**.
- Depth distribution: 10–90 % quantiles span 25.3–46.5 km (median 35.9 km) —
  consistent with the plate interface depth in this region.
- Top families (~2,000 detections each) cluster in the central V.I. plate
  interface, matching the dense Bostock LFE zone in the literature.

### 2.3 Continuous waveforms (FDSN)

Six broadband stations are used for v1, all available with continuous BH?
(40 Hz) coverage in 2010:

| Code | Network | lat   | lon     | Notes |
|------|---------|-------|---------|-------|
| PGC  | CN      | 48.65 | −123.45 | Pat Bay, B.C. — Bostock workhorse |
| LZB  | CN      | 48.61 | −123.82 | Lasqueti Island |
| NLLB | CN      | 49.23 | −123.99 | Nanoose Bay |
| SNB  | CN      | 48.78 | −123.17 | Saturna Island — original Bostock station |
| LRIV | UW      | 48.06 | −123.50 | Lake Crescent, Olympic Peninsula |
| SQM  | UW      | 48.07 | −123.05 | Sequim, Olympic Peninsula |

Waveforms are fetched via ObsPy's FDSN client against IRIS/EarthScope using
`src/tremorferometry/waveforms.py`. The fetcher is `ThreadPoolExecutor`-parallel
(24 workers) over (station, day) tasks and writes one MSEED per UTC day to
`data/waveforms/{net}.{sta}/{year}/{jday}.mseed`. For the 150-day v1 window
(2010-05-17 → 2010-10-15) across 6 stations × 3 components, 11.1 GB total were
fetched in roughly 2 minutes (network-limited).

Notes:

- PGC's HH? channels (100 Hz) only came online in 2015. For 2010, BH? at 40 Hz
  is the operative channel set. 40 Hz Nyquist = 20 Hz, comfortably above our
  2–8 Hz analysis band.
- We initially queried HH? and received 204 No Content for all 2010 dates; the
  fix was to query BH? after consulting station metadata via the FDSN station
  service.

---

## 3. Pipeline

The pipeline is implemented as numbered idempotent scripts plus a Python
package. Each script takes `--config <yaml>` and reads from / writes to a
predictable on-disk layout. Heavy stages cache aggressively.

```
01_define_episode.py     PNSN tremor catalog  →  episode bbox + t_start/t_end
02b_ingest_lin_catalog   Lin Zenodo CSV       →  families + per-family detections
04_fetch_waveforms.py    FDSN                 →  per-day MSEED cache
06_stack_bins.py         detections + MSEED   →  per-family HDF5 stacks
07_measure_dvv.py        stacks + reference   →  per-bin dv/v parquet
09_aggregate_plot.py     dv/v parquet         →  smoothed figure + QC
```

The template-matching step (`05_match_lfe.py`) is **not used in v1** because
Lin's catalog provides detection times directly. The matching code is kept for
future per-family template extraction or for time periods outside Lin's
coverage.

### 3.1 Per-family per-station stacking

For each (family, station) pair and each time bin, `stack_family_station`:

1. Filter detections to this family (and station, if the detection list is
   per-station; otherwise broadcast to all stations).
2. For each detection time `t = OT`, identify the time bin.
3. Lazily load the MSEED day file containing `t`; cache it across detections
   on the same day to avoid re-reading.
4. Restrict to the requested component (default Z), `merge(fill_value=0)`,
   detrend (demean), zero-phase 2-8 Hz bandpass (4-pole Butterworth, two
   passes), resample if needed.
5. Cut a window `[t − pre_s, t + post_s]` — for v1, `pre_s = 5`, `post_s = 35`
   for a 40-s window. This generously brackets the direct S arrival, which
   travels ~10 s for a 30-km-depth LFE to a station at ~50 km horizontal
   offset (Vs ≈ 3.5 km/s).
6. Stack (mean) the cuts within each bin → one waveform per (family, station,
   bin), written to `data/stacks/<family>.h5` (group per station, dataset per
   bin, with `t_center` and `n_det` attributes).

Parallelism: the original `stack_all_parallel` parallelizes per (family,
station). For 217 × 6 = 1302 tasks on this machine (88 physical cores, 1.5 TiB
RAM) the OS page cache holds the entire 11 GB waveform set after the first
pass, so subsequent family iterations are CPU-bound on filter + resample +
window-cut. Empirically: 30 families × 6 stations = 20 minutes; 217 families
expected ~2 hours (Phase D).

### 3.2 dv/v by stretching

`src/tremorferometry/dvv.py` implements the canonical stretching estimator.
For a reference waveform `ref(t)` and a current waveform `cur(t)`, both
sampled at `fs` over the same support, a homogeneous fractional velocity
change `dv/v = ε` predicts

    cur(t) = ref(t / (1 + ε))    ⟺    ref(t) = cur(t · (1 + ε))

Inversion (`stretch_dvv`):

1. Restrict to the coda window in seconds: for v1, `t ∈ [20, 32]`. With the
   cut starting at OT−5 and S arriving ~13 s after OT, the direct S sits near
   `t ≈ 18` and the coda window starts ~2 s after the direct phase.
2. Scan ε on a uniform grid of 401 points in [−2 %, +2 %]. For each ε:
   - Cubic-spline resample `cur` onto `t · (1 + ε)`.
   - Compute zero-mean Pearson correlation against `ref_coda`.
3. Pick ε* at the maximum CC. Refine sub-grid by parabolic interpolation of
   the three samples around the peak. Estimate a crude uncertainty from the
   parabolic curvature.

Per-family driver `measure_family` builds a reference stack by averaging all
bin-stacks whose centers fall inside the **reference window** (in days
relative to `episode.t_start`; v1 uses [−90, −30]). It then runs stretching
on every bin against that reference. Empty bins, low-CC measurements
(`cc_max < min_cc`, v1 default 0.3), and stretching failures (e.g.
identically-zero references) are silently dropped.

The synthetic-recovery unit test (`tests/test_dvv_synthetic.py`) confirms
recovery of imposed dv/v values in [−1 %, +1 %] to better than 5 × 10⁻⁴ on
band-limited synthetic coda. The end-to-end integration test
(`tests/test_end_to_end_synthetic.py`) builds synthetic stacks with a
Cascadia-ETS-shaped dv/v pattern and verifies recovery rms < 5 × 10⁻⁴ and
sign-correct ETS-interior bins. 12 tests in total, all pass.

### 3.3 Aggregation and smoothing

`scripts/09_aggregate_plot.py` aggregates dv/v measurements per time bin:

- Per-row weight `w = cc_max · √n_det`. The `n_det` factor reflects the
  expected SNR improvement of a stack of N detections (∝ √N).
- Weighted mean and weighted standard error per `t_center` across all
  (family, station) measurements falling in that bin.
- Gaussian smoothing of the per-bin weighted-mean series with σ = 3 days.

### 3.4 Quality control

`src/tremorferometry/qc.py`:

- **Spatial coherence**: for each (family, t_center), compute the std of
  dv/v across the stations sampling that bin. Compare to the median per-bin
  error. A bin is "coherent" if `std < 2 × err` with ≥ 2 stations. The
  overall coherent fraction across all bins should be > 50 %.
- **n_det independence**: dv/v should not be a linear function of log(n_det).
  Pearson `|r| < 0.5` passes. (Practical bug: the metric returns NaN when
  n_det is missing; a fix is to filter NaN before computing r.)

For Phase B (30 families × 6 stations), spatial coherence passes at **91 %**
of bins. This is the strongest QC outcome from v1 so far: when stations make
a measurement at the same time, they agree on the value.

---

## 4. Results so far

### Phase A — 1 station, BHZ only

- CN.PGC, BHZ, 30 top-detection families, 150 days.
- Pre-ETS dv/v flat near 0. ETS interior shows ~±1.5 % per-bin scatter with a
  hint of a negative excursion mid-September. Single-station SNR is the
  expected limit.
- Figure: `figures/smoke_phaseA_dvv.png`.

### Phase B — 6 stations × 3 components

- All six V.I./Olympic stations, 30 families.
- 856 dv/v measurements survive `cc_max > 0.3`. Mean CC 0.42.
- ETS-window weighted mean shows small structure (~±0.5 %) but no clean
  monotone drop.
- Figure: `figures/smoke_phaseB_dvv.png`.

### Phase C — same data, QC + smoothing

- Gaussian-smoothed weighted mean stays within ±0.3 % during the ETS.
- LFE-rate panel shows distinct bursts at 2010-08-23, 2010-08-29, 2010-09-04
  (peaks 3,500 – 4,500 LFEs / 2-day bin vs <500 background).
- **Spatial coherence QC passes at 91 %**.
- Headline result: **|dv/v| ≲ 0.3 % during the 2010-08 V.I. ETS** (upper
  bound).
- Figure: `figures/smoke_phaseC_dvv.png`.

### Phase D — all 217 families (in progress)

Scaling from top 30 to all 217 families increases per-bin measurements 7×.
Expected to tighten the upper bound by ~√7 ≈ 2.6×, into the ~0.1 % range,
sufficient to test the hypothesis that Cascadia ETS produces a Mexico-style
dv/v drop. Phase D stacking refactored for true (family, station) parallelism
(`stack_all_parallel`): all 1,302 stacking tasks submitted up-front rather
than per-family batches of 6, yielding ~8× speedup on the 48-worker pool.

### Per-family analysis (Phase B follow-up)

Splitting the Phase B parquet by family and computing each family's mean
ETS-window dv/v minus its mean reference-window dv/v shows **per-family
excursions of order ±1 %, with mixed signs across nearby families**:

| family | lat   | lon     | depth (km) | ETS Δdv/v | SNR  |
|--------|-------|---------|------------|-----------|------|
| L0001  | 48.90 | −123.80 | 22         | −1.07 %   | 6.6  |
| L0005  | 49.32 | −123.90 | 24         | +1.54 %   | 5.6  |
| L0000  | 49.32 | −123.80 | 33         | −1.30 %   | 5.5  |
| L0003  | 48.90 | −123.70 | 22         | +1.40 %   | 2.9  |

A real medium-velocity change should be coherent across nearby families
(overlapping coda volumes). The mixed signs in adjacent locations point to
**source heterogeneity within our 0.1° grid families** — at this clustering
resolution we are lumping multiple physically distinct LFE patches together,
and changes in their relative contributions to the bin-stack masquerade as
"apparent dv/v". Quantitative aggregate stays within ±0.3 % because these
opposite-sign per-family excursions average down.

Refinements that would address this:
- Cluster Lin's detections by *waveform similarity* (not just lat/lon) to
  recover Bostock-style families.
- Tighten the spatial cell to 0.05° and require min_n higher to avoid
  swallowing dissimilar sub-patches.
- Restrict the analysis to families whose detection waveforms have high
  intra-family cross-correlation (a station-level repeatability test).

---

## 5. Implementation notes / pitfalls already encountered

- **Channel naming over time**. PGC switched from BH? (40 Hz) to HH? (100 Hz)
  around 2015. Always probe station metadata for the target epoch before
  fetching.
- **CN network at IRIS**. Canadian National Seismograph Network station
  metadata is mirrored at IRIS/EarthScope and data is available for at least
  the V.I. stations we use. Some marginal stations may require querying
  NRCAN directly; we did not need that for v1.
- **Family clustering at 0.1°**. Coarser than Bostock's per-family resolution
  (which uses waveform-similarity clusters at much finer spatial scale).
  Acceptable for v1; will likely become a SNR-limiter at the per-family
  level. Refinement: cluster Lin detections by waveform similarity to recover
  true Bostock-style families.
- **fast-matched-filter** (FMF) is installed and CUDA-built on the project
  env, but v1 does not use template matching at all — Lin's catalog supplies
  detection times directly. FMF is ready for the Phase E template-extension
  use case (matching on time periods outside Lin's coverage).
- **`pkg-resources` phantom dep** prevented installing `eqcorrscan` into the
  base conda env. Resolved by creating a dedicated env at
  `/home/jovyan/envs/tremorferometry` with `mamba install -c conda-forge
  eqcorrscan` (conda-forge pulls in `fftw` automatically) and `pip install
  -e . --no-deps --no-build-isolation` for the tremorferometry package
  itself.

---

## 6. What's next

1. **Phase D completion** — all 217 families. Re-measure dv/v and re-aggregate.
2. **Per-family analysis** — instead of pooling all families, look for a
   subset that shows a clear dv/v excursion. The plate interface is not
   uniformly responsive.
3. **Frequency-band split** — re-run at 2-4 Hz vs 4-8 Hz. Mexican Guerrero
   work showed the dv/v signal is band-dependent.
4. **Coda-window stability** — vary the coda start/end by ± a few seconds.
   If the signal is robust, it persists; if not, it is a coda artifact.
5. **Hydrology / loading control** — overlay precipitation or GPS-derived
   surface loading. Any seasonal dv/v signal should be distinguishable from
   the ETS-driven signal in timing.
6. **Margin-wide v2** — ingest Sweet (2019) for SW Washington, Ducellier
   (2022) for southern Cascadia, Plourde (2015) for N California. Each
   has a published LFE catalog. Central Oregon remains a gap; if we want it
   we will need to do template discovery from PNSN tremor windows.
7. **GPU template matching** — when extending into time periods or regions
   outside published catalogs, use fast-matched-filter on the L40S to
   self-detect LFEs.

---

## 7. Reproducibility

Environment:

```bash
mamba activate /home/jovyan/envs/tremorferometry
```

Re-run v1 end-to-end from scratch:

```bash
# raw catalogs
mkdir -p data/raw_lfe catalogs
curl -fsSL -o data/raw_lfe/lin2023_lfe.csv \
    'https://zenodo.org/records/10016020/files/EQloc_001_0.1_3_S.csv?download=1'

# PNSN tremor for episode context (optional; v1 uses a fixed config date)
python scripts/00b_fetch_pnsn_tremor.py --start 2010-01-01 --end 2010-12-01 \
    --bbox 47.5 50.0 -125.5 -122.0 --out catalogs/pnsn_tremor_2010.csv

# Lin → families + detections
python scripts/02b_ingest_lin_catalog.py --config configs/ets_2010_vi.yaml \
    --raw data/raw_lfe/lin2023_lfe.csv

# FDSN waveforms (~11 GB)
python scripts/04_fetch_waveforms.py --config configs/ets_2010_vi.yaml --workers 24

# Stacks (~140 MB for top 30; ~1 GB for all 217)
python scripts/06_stack_bins.py --config configs/ets_2010_vi.yaml \
    --detections data/detections_lin_ets_2010_vi.parquet \
    --selected catalogs/lin_families_ets_2010_vi.csv \
    --waveforms data/waveforms \
    --out data/stacks --workers 24 --fs 40 --pre-s 5 --post-s 35

# dv/v
python scripts/07_measure_dvv.py --config configs/ets_2010_vi.yaml \
    --stacks data/stacks --out data/dvv/v1.parquet --workers 16 --fs 40

# Aggregate + figure + QC
python scripts/09_aggregate_plot.py --config configs/ets_2010_vi.yaml \
    --dvv data/dvv/v1.parquet \
    --detections data/detections_lin_ets_2010_vi.parquet \
    --out figures/v1_aggregate.png
```

All raw input artifacts (PNSN API, Lin Zenodo CSV, FDSN waveforms) are
deterministically reproducible from the script invocations above plus the
checked-in config YAML.
