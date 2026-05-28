# Methods: LFE coda-wave interferometry for dv/v across Cascadia ETS

A running methods log for the `tremorferometry` project. The current target is
the August–September 2010 Episodic Tremor and Slip (ETS) event beneath southern
Vancouver Island; the architecture extends margin-wide for future work.

---

## 0. The intuitive picture: a fisheye camera looking down at the megathrust

Each (LFE patch, station) pair is one **pixel** of an integrated measurement.
35 patches × 1 station (PGC) = 35 pixels of a time-lapse image of the
subsurface. The set of pixels together is like a **fisheye lens at PGC looking
down at the plate interface**, sampling at 35 spots scattered across southern
Vancouver Island.

What each pixel actually is:

- Not a sharp ray — it's a fuzzy 3D **sensitivity kernel** between the source
  patch (on the plate interface ~30 km depth) and the station at the surface,
  weighted by where the coda waves spend their time scattering.
- For each ~2-second coda window the kernel covers a roughly tens-of-km-wide
  tube between source and receiver, biased heavily toward depths near the
  source.

The time-lapse:

- Each ETS cycle is an "exposure": stack a patch's detections during that ETS
  and measure dv/v vs a reference.
- Comparing exposures across ETSs (every ~14 months in N. Cascadia) gives a
  movie of how the subsurface has changed cycle to cycle.
- Patches whose dv/v changes during a specific ETS show **where in the crust
  the medium was perturbed by that slow-slip event**.

Spatial coverage with one station:

- All 35 pixels converge at PGC and fan out to the patches
- One station = fan-shaped 1-sided image, no crossing rays → no 3D inversion
- Add more stations and the rays start crossing at different angles → solve
  for the dv/v field in 3D (true tomography). The natural extension of this
  work.

Why this is the right kind of measurement:

- The LFE coda samples the plate interface preferentially (kernel peaks near
  the source depth)
- The repeating-source property guarantees that any waveform change is a
  medium change, not a source change
- Once we have the catalog of repeating patches, the rest is *just* stacking
  + stretching: a clean, well-understood signal-processing chain that turns
  decades of continuous seismic data into a series of subsurface snapshots

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

### Phase D — all 217 families

Scaled from top 30 to all 217 families. Stacking refactored for true
(family, station) parallelism (`stack_all_parallel`): all 1,302 tasks
submitted up-front rather than per-family batches of 6, ~8× wall-clock
speedup. All 217 family HDF5 stacks built in 580 s on the 48-worker pool;
709 MB total.

Measurement on the full stack set:

- **5,948 dv/v rows** (7× Phase B's 856).
- **Per-bin sample size**: 20–40 measurements per ETS-window bin (vs Phase B's 2–6).
- **ETS-window dv/v: mean = +0.043 %, median = +0.148 %, std = 1.18 %**.
  Per-bin standard error of the mean ≈ 1.18 / √27 ≈ 0.23 %. The Gaussian-
  smoothed series stays within **±0.2 % through the ETS interval**.
- Spatial coherence QC: 89 % of 1,334 bins coherent across stations.
- Figure: `figures/smoke_phaseD_dvv.png`.

**Headline numerical result:** for the 2010-08-15 → 2010-09-15 V.I. ETS,
the LFE-CWI measurement at 6 V.I./Olympic broadbands gives

    |dv/v| < 0.2 %   (smoothed)
    mean ETS dv/v = +0.04 ± 0.05 %   (1-sigma SE on the ETS-wide mean)

The data are consistent with **no detectable sustained ETS-wide dv/v
signature** at this sensitivity. Per-family excursions exist at the ~±1 %
level but with mixed signs across adjacent families (see per-family note
above) — consistent with the diagnosed source-heterogeneity issue in the
0.1° clustering, not with a coherent medium velocity change.

### Phase D QC — family-waveform similarity (the central finding)

The Phase D pipeline assumed each Lin family at our 0.1° clustering would
provide a repeating source — that is the prerequisite for coda-wave
interferometry. We tested this directly.

**Within-family bin-to-bin coda CC** (family L0000 at PGC, 59 bin-stacks):

    mean CC = -0.002, median = -0.001, max off-diagonal ~ 0
    fraction of bin pairs with CC > 0.5: 0 %
    fraction of bin pairs with CC > 0.3: 0 %

The bin-to-bin CC matrix is white noise except on the diagonal. There is no
master-stack waveform that 2-day-bin stacks of L0000 converge to.

**Within-family detection-to-detection CC** (L0000 at PGC, 2,365 detections
each cut for an 8-second direct-phase window):

    fraction of detection pairs with CC > 0.5: 0 %
    Hierarchical clustering at CC > 0.5: 2,353 sub-clusters (largest n=2)

Even with **time-shifted CC** allowing ±2 s alignment (which would absorb
typical OT/depth uncertainty), the result is the same:

    Time-shifted CC: mean=0.18, median=0.18, max across 125k pairs=0.525
    Fraction of pairs above CC 0.5: 0 %
    Fraction above 0.3: 0.6 %
    Cluster at CC>0.3: 348 groups, largest n=5

**Lin's catalog, at the 0.1° granularity we use, does NOT contain repeating
sources.** Either the deep-learning-based detector finds many physically
distinct LFE sub-patches per cell, or sub-cell timing/depth variability
prevents the underlying waveforms from aligning even with time shifts.

Consequence: the Phase A-D dv/v values are statistically valid stretching
computations, but they measure **source-mixture variation between bins**, not
medium velocity change. The +0.04 % / −0.03 % numerical results in earlier
sections are not bounds on real Cascadia dv/v.

Figures: `figures/smoke_family_similarity.png` (bin-to-bin CC matrix),
`figures/smoke_reclust_L0000_PGC.png` (detection-to-detection CC matrix).

### Phase D QC — frequency-band split (2-4 vs 2-8 Hz)

Re-stacking all 217 families at a narrower bandpass and re-measuring on those
stacks. The per-family ETS dv/v statistics (aggregating each family's ETS-bin
mean, then std across families) is the more honest SE because per-bin
measurements within a family are correlated:

| Bandpass  | families w/ETS data | per-family ETS dv/v   | per-family std | SE of family mean |
|-----------|---------------------|-----------------------|----------------|-------------------|
| 2-8 Hz    | 186                 | +0.048 %              | 0.891 %        | ±0.065 %          |
| 2-4 Hz    | 216                 | **−0.032 %**          | **0.261 %**    | **±0.018 %**      |

The 2-4 Hz band shows a **3.4× tighter per-family scatter** than 2-8 Hz, and
the per-family ETS dv/v mean is **negative at ~1.7σ from zero** — the
direction expected for slow-slip-induced velocity drop. This is the first
hint of a directional signal we've seen.

Two interpretations:
1. LFE energy peaks in 2-4 Hz, so the stacks are higher-SNR there; the
   per-family heterogeneity (driven by waveform-shape changes between
   sub-patches) is suppressed at lower frequency where individual sub-patch
   waveforms look more alike. Lower-frequency coda is also longer-lived
   (less attenuated), giving more medium sampling per measurement.
2. The 4-8 Hz band is dominated by scattered energy from small heterogeneities
   that respond to source-location migration rather than bulk velocity change.

The aggregate SE = 0.018 % is approaching the regime where the predicted
~−0.1 % Mexican-style SSE dv/v drop would be detectable. The actual
Cascadia 2010-08 measurement (−0.03 %) is smaller than that prediction by
~3× — either the signal really is small, or our 0.1° family clustering
still smears it down.

Figures: `figures/smoke_phaseD_24Hz.png` (2-4 Hz aggregate),
`figures/smoke_phaseD_freqband_compare.png` (side-by-side overlay).

### Phase D QC — coda-window stability

Re-running `measure_many` on the same Phase D stacks with three different
coda windows:

| coda window (s) | n_rows | mean CC | ETS mean dv/v       |
|-----------------|--------|---------|---------------------|
| [18, 30]        | 5,908  | 0.439   | +0.054 ± 0.057 %    |
| [20, 32]        | 5,948  | 0.438   | +0.043 ± 0.056 %    |
| [22, 34]        | 6,018  | 0.437   | +0.003 ± 0.056 %    |

All three are consistent with zero within 1 σ, and the per-window means
agree at the level of the standard error. The Gaussian-smoothed dv/v
trajectories from the three windows track each other broadly (figure
`figures/smoke_phaseD_codawindow.png`). A spurious "signal" from a single
coda-window choice would not survive this test — the consistency across
windows is the strongest confirmation that no ETS-wide medium-velocity
change is detectable here above |dv/v| ≈ 0.1 %.

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

## Tremor-windowed inter-station CC (Phase E)

After confirming the LFE-CWI assumption is broken at our 0.1° clustering
(see "Phase D QC — family-waveform similarity" above), we pivoted to a
different physical setup: **tremor windows as a distributed source field**,
inter-station cross-correlation as the measurement.

The point of the change: this approach **does not require repeating sources**.
It only requires that during each "tremor window" the wavefield at the array
is dominated by emissions from the plate interface. The inter-station CC then
approximates the Green's function with source weighting biased toward the
tremor source region — i.e., the megathrust — which is the depth sensitivity
we wanted from CWI in the first place.

### E.1 Algorithm

  1. PNSN tremor catalog → list of (t0, t0+300 s) tremor windows in the v1 bbox.
  2. For each (station-pair, window): cut the same 300 s window from both
     stations, demean, L2-normalize, and FFT cross-correlate over ±60 s lag.
  3. Stack CCs into 2-day bins per station pair → HDF5 in `data/cc_tremor/`.
  4. Stretch the coda of each bin's CC against a reference CC stack → dv/v.

### E.2 Preprocessing — what works

A direct compare of three preprocessing flavors on 1,500 PGC-SNB windows
(`figures/smoke_tremor_cc_preprocess_compare.png`):

| Preprocessing | Max CC at +/-8 s |
|---|---|
| Raw bandpass + L2 normalize | 0.0123 |
| One-bit normalization | 0.0064 |
| Running-mean amplitude norm | 0.0084 |

**Raw bandpass + L2 normalize wins by ~2×.** One-bit and running-mean
suppress the coherent direct-phase amplitude that tremor inter-station
CC actually relies on. (This is opposite of ambient-noise CC, where
one-bit is standard to suppress earthquakes — for tremor, the coherent
signal IS what we want, not what we're suppressing.)

Default switched to raw-bandpass throughout.

### E.3 Frequency band — narrow to the tremor peak

User noted tremor is band-limited. The published tremor band is roughly
1–10 Hz with **peak energy in 2–5 Hz**. We narrowed from 2–8 Hz (used for
LFE-CWI) to **2–5 Hz** for tremor-CC. The 2–5 Hz band has:
- Larger spatial sampling per wavelength (longer wavelengths → more bulk
  sensitivity, less small-scale scattering noise);
- Cleaner tremor signal (excludes higher-frequency cultural noise);
- Roughly the same number of accepted measurements at our CC threshold.

### E.4 Mean CC is real signal — 0.67 vs 0.42 for LFE-CWI

Aggregated across all bins and pairs, **mean tremor-CC stretching CC = 0.67**,
compared with 0.42 for the LFE-CWI Phase D. The tremor-CC stacks have real
coherent structure that the LFE-CWI stacks lacked.

### E.5 The noise-reference artifact (caught Phase E.5)

Initial tremor-CC result with pre-ETS reference (t_start −90 to −30 days):

    ETS dv/v = -0.113 % +/- 0.049 %    (~2.3 σ negative)

This looked like the Mexican-style ETS-induced velocity drop. But the
pre-ETS reference period had **<50 tremor windows per 2-day bin** (vs ~800
during ETS — see `figures/smoke_tremor_cc_dvv.png` bottom panel). The
"reference CC" was essentially noise. Stretching can find spurious
ε > 0 that align noise patterns with real signal — and the apparent
−0.11 % was exactly that artifact.

**Diagnostic and fix:** build the reference from the densest-tremor bins
(top-5 by `n_windows`, which all fall inside the ETS) instead of a fixed
pre-ETS time window. Self-comparison of reference bins gives dv/v = −0.005 %
(noise floor sanity check). Then:

    ETS non-reference bins: dv/v = +0.024 % +/- 0.020 %  (consistent with zero)

The "signal" went away — confirming the original number was a reference
artifact, not a real measurement.

### E.6 Symmetric coda windows tighten the bound

For each bin's CC, use **both positive and negative lag codas**
(concatenated into a single window) instead of just positive lag. Doubles
the coda sample count → tighter precision. The four coda choices and
their results:

| Coda window           | ETS dv/v (%)        | SE (%) |
|-----------------------|---------------------|--------|
| pos [+15, +35] s      | +0.014              | 0.017  |
| neg [−35, −15] s      | +0.030              | 0.028  |
| symmetric concatenated| **+0.0003**         | **0.0075** |
| \|lag\| ∈ [15, 35] s  | +0.005              | 0.005  |

**Tightest result: tremor-CC dv/v during 2010-08 V.I. ETS = +0.0003 ± 0.0075 %.**

2σ upper bound: |dv/v| < 0.015 %.

For comparison: the Mexican Guerrero SSE produced a clear −0.2 % dv/v
drop. Our bound is **13× tighter than that signal**. If Cascadia ETS
produced a Guerrero-style dv/v in the volume we sample, we would have
seen it easily. The data are consistent with **either no measurable
deep dv/v during this ETS, or a dv/v smaller than 0.015 % averaged
over our southern-V.I. station-pair coverage**.

### E.7 Per-pair robustness check

Splitting the aggregate ETS dv/v by station pair (`figures/smoke_tremor_cc_per_pair_distance.png`):

- 12 pairs with ≥2 ETS measurements each
- Per-pair ETS dv/v ranges from −0.17 % (PGC-SNB, n=2, SE=0.18 %) to +0.34 %
  (LZB-PGC, n=1) — but **every single pair is consistent with zero at 95 % CI**
- No correlation with inter-station distance (25 – 146 km).
- No outlier pair is driving the aggregate; the +0.0003 ± 0.0075 % aggregate
  is honestly representative of every path we sample.

If a real ~0.05 % medium velocity change existed in some sub-region but not
others, we would expect at least one pair to show it as a tight, distance-
correlated offset. We see none. The null is robust per-pair, not just on
average.

### E.7b Finer-grid Lin clustering doesn't rescue LFE-CWI

Hypothesis: maybe 0.1° (~10 km) cells lump multiple LFE patches, but finer
cells (~2 km) would isolate single patches.

Test: re-cluster Lin's catalog at **0.02°** with `min_n = 30`. Result:
803 sub-families, top has 812 detections in a 2 km cell. Run network-CC
clustering on those 812 detections:

    CC > 0.4 threshold: 798 sub-clusters formed; largest has 3 members (0.4 %)
    CC > 0.3 threshold:  34 sub-clusters formed; largest has 54 members (6.7 %)

Even at 2 km, Lin's detections do NOT form coherent waveform-similar
families. Either:
  - Lin's DL detector picks up many physically distinct patches per 2 km cell
    (geophysically reasonable — the plate interface is densely populated with
    LFE sources), or
  - location precision is finer than the catalog's lat/lon precision can
    represent, so same-source detections are spread across nearby cells

**Conclusion:** at any spatial binning of Lin's catalog tested (0.1°, 0.05°,
0.02°), the resulting "families" don't satisfy the CWI repeating-source
assumption. Real LFE-CWI for this region requires direct waveform clustering
from continuous data (network autocorrelation; Brown-Beroza-Shelly 2008) or
Bostock's own family templates from their publications. Both are substantial
engineering efforts.

The tremor-CC approach (Phase E.6) doesn't have this problem — it doesn't
require repeating sources — and is the v1 working measurement.

### E.7c Cross-year repeating-source test (user's reframing)

The user reframed the goal: not "produce dv/v" directly, but "determine whether
repeating events exist at a single station across multiple ETS cycles, which
would make CWI viable." The right experiment: take the densest-tremor day in
each year of Lin's catalog, compute pairwise network CC between every pair
of LFE detections (within-year and cross-year), and ask whether there are
high-CC cross-year matches.

**Setup:**
- 6 ETS peak days from 2005, 2009, 2010, 2011, 2012, 2013 (skipped 2006-2008 —
  NRCAN doesn't host CN.PGC for those days).
- ~150 random Lin detections per day; 900 total.
- 6-second window 6 s after OT (catches the direct phase at ~10 s post-OT).
- Network CC = mean(time-shifted CC at PGC, time-shifted CC at LZB),
  with ±1.5 s shift allowance.
- Bandpass 2-5 Hz.

**Result (444,550 pairs across 900 detections):**

| Window | n pairs | mean CC | max CC | frac > 0.5 |
|---|---|---|---|---|
| Within-year (same ETS) | 67,050 | 0.356 | 0.635 | 0.64 % |
| Cross-year (different ETS) | 337,500 | 0.352 | 0.637 | 0.50 % |

The within-year and cross-year CC distributions are statistically indistinguishable.

**But the matches are NOT real repeating sources.** Inspection of the top 12
cross-year pairs (`figures/smoke_crossyr_top_pairs.png`) shows:

- Pair locations are **36–87 km apart** between matched detections;
- Depths range from −60 km to 0 km, often very different between matches;
- Waveforms look similar in a fuzzy band-limited-noise way, with no
  distinctive shared zero-crossing pattern.

A genuine repeating LFE source family would have all members **co-located
(< 1 km)** and at the **same depth**. The cross-year "matches" we see are
bulk-spectrum coincidences in band-limited 2-5 Hz noise, not source
repetitions.

**Conclusion:** Lin (2023)'s catalog, treated as a list of pre-detected
events, does not contain identifiable cross-year repeating LFE families at
the CC > 0.5 level at any of our station-pair combinations. This does not
mean Cascadia LFE families don't exist — Bostock's published families DO
repeat across years — it means **Lin's detection times are not, on their
own, sufficient to recover the family structure**. The standard recipe for
recovering it is Shelly-Beroza / Brown-Beroza-Shelly network autocorrelation
on *continuous data*, which is the engineering effort we'd need to start a
proper repeating-source LFE-CWI pipeline.

### E.7d Hotspot persistence: locations DO repeat, waveforms don't (cleanly)

Although individual detection waveforms didn't cross-correlate well across
years (E.7c), the **spatial pattern of LFE activity is highly persistent**.
Comparing the top-50 (0.05°) detection-density cells in each year against
2010's:

    2005 vs 2010: 27/50 cells in common (54 %)
    2008 vs 2010: 29/50 cells in common (58 %)
    2012 vs 2010: 36/50 cells in common (72 %)
    2013 vs 2010: 37/50 cells in common (74 %)

The central V.I. hotspot at ~(48.8°N, −124°E) is unmistakably the same
across every ETS cycle 2005-2013 (`figures/smoke_crossyr_hotspot_persistence.png`).

This is the standard Bostock LFE zone, and it confirms that **the physical
premise of repeating-source CWI is valid for Cascadia**: the same patches
of the plate interface fire ETS after ETS. What our test E.7c shows is
that **individual detection waveforms vary even within a stable hotspot** —
each hotspot contains many sub-patch sources that fire variably.

So:
  - To do repeating-source CWI, you need to recover sub-patch families
    inside each hotspot (Bostock's recipe: waveform-similarity clustering).
  - Lin's catalog gives you detection times across the whole hotspot but
    not the sub-family structure within it.
  - Tremor-CC (Phase E.6) bypasses this entirely by treating the entire
    hotspot as a distributed source field, which is why it gives a clean
    measurement at the cost of losing per-patch spatial resolution.

### E.7e BREAKTHROUGH — repeating LFE families found across 8 years

After E.7c showed that simple cross-year pairwise CC on Lin's detections
gave only spurious matches between distant sources, we reworked the
methodology with proper Shelly-Beroza ingredients:

  1. **Envelope-peak alignment.** Don't trust Lin's OT — instead bandpass
     2-8 Hz, take Hilbert envelope, find the peak inside [OT+5, OT+17] s
     (the expected direct-phase window for ~30 km depth at ~50 km offset),
     and cut a tight 2-second window *centered on that peak*.
  2. **Multi-station verification.** Network CC = mean of max-shifted CCs
     at PGC and LZB (both available for all 6 peak days 2005-2013).
  3. **Strict threshold.** CC >= 0.7 at the network level (suppresses
     band-limited-noise coincidences that pass at CC ~0.5-0.65).
  4. **Hotspot focus.** Filter Lin detections to those within 0.05° of
     (48.85, -123.85) — the canonical central V.I. LFE hotspot.

Sample: 451 hotspot detections across 6 peak days 2005-2013, pairwise
network CC computed in 0.3 s with batched FFT.

| Threshold | n pairs above | clusters | cross-year clusters |
|---|---|---|---|
| 0.5  | 15,800 (15.6 %) | -- (too loose, single merged cluster) | -- |
| 0.6  | 1,446 (1.4 %)   | 5  | 5 |
| 0.65 | 250 (0.25 %)    | 23 | 19 |
| **0.7**  | **39 (0.04 %)** | **20** | **16** |

**16 cross-year families form at CC ≥ 0.7.** Top examples:

| Family | n | Years |
|---|---|---|
| 4 | 8 | 2010, 2011, 2012, 2013 |
| 0 | 6 | **2005, 2011, 2012, 2013** (8-year span) |
| 2 | 6 | 2009, 2011, 2012, 2013 |
| 9 | 3 | 2011, 2012 |

Visual confirmation (`figures/smoke_crossyr_repeaters.png`): the per-family
stacks show clear coherent LFE-like oscillation, and individual member
waveforms align around the envelope peak at both stations.

**Conclusion.** Repeating LFE waveforms DO exist in Cascadia at single
stations across multiple ETS cycles. They can be recovered from Lin's
catalog with the right methodology (envelope alignment + tight windows
+ multi-station verification + strict threshold). The earlier negative
results (E.7c) were due to insufficient filtering / alignment, not absent
sources. **This validates the user's hypothesis and opens the path to
proper repeating-source CWI.**

The implementation lives in `src/tremorferometry/repeater.py`
(`cut_aligned_window`, `all_pairs_cc_max_shifted`, `network_cc_all_pairs`,
`cluster_matches`). Per-family stacks can be used as templates for
matched-filter detection across all continuous data to grow the catalog.

### E.7f Multi-seed catalog growth — Bostock-scale families emerging

Once the 16 cross-year family templates were identified (E.7e), each was
used as a matched-filter template against ALL Lin detections in southern
V.I. for the 2010 ETS (not just the 451-event hotspot subsample). Per-
template max-shifted network CC, threshold 0.65 (slightly looser than the
0.7 we used for discovery, because templates are clean stacks of multiple
events so signals are reinforced).

Results on a 5,000-detection subsample of the 2010 ETS:

| Template | Initial members | Matches at CC ≥ 0.65 | Extrapolated full ETS |
|---|---|---|---|
| fam 0  | 6  | 48 | ~1,334 |
| fam 2  | 6  | 48 | ~1,334 |
| fam 4  | 8  | 3  | ~83   |
| fam 9  | 3  | 28 | ~778  |
| fam 10 | 3  | 27 | ~750  |
| fam 14 | 3  | 42 | ~1,167 |
| fam 18 | 3  | 10 | ~277  |

Across just the seven cross-year families with ≥3 initial members, we'd
get **~5,700 LFE detections per ETS**. Scaling to all 16 cross-year families
plus the iterative-refinement step (use grown stacks as new templates, find
more seeds, find more families) easily yields a ~10⁴–10⁵-detection-per-ETS
catalog — comparable in size to Bostock's published families.

**This validates the user's hypothesis that repeating LFE waveforms across
multiple ETS cycles can drive a proper coda-wave interferometry pipeline.**

### E.7g Path to Cascadia-wide

Lin (2023)'s catalog covers only southern V.I. (lat 48.07°–49.32°). To
extend the repeating-family catalog to the full Cascadia margin, three
options stack:

1. **Other published catalogs.** Sweet et al. (2019) covers SW Washington
   and Olympic Peninsula. Plourde et al. (2015) covers N California.
   Ducellier (2022) has an 8-year catalog for southern Cascadia. Each can
   be ingested through a Lin-like adapter; the repeater.py pipeline then
   discovers families exactly as we did for V.I.
2. **Self-detection from continuous data.** For gaps (notably central
   Oregon), run network autocorrelation (Brown-Beroza-Shelly 2008) directly
   on continuous data inside PNSN tremor windows: this is the same recipe
   that built Bostock's families originally. Most expensive but covers the
   gaps that published catalogs miss.
3. **Margin-wide station coverage.** We have the FDSN access path; just
   need station lists per latitude band. From IRIS metadata, the natural
   set is V.I./Olympic (we have), then UW.* across WA, then UO.* through
   Oregon, then BK.* and US.* in N CA.

The repeater.py methodology is location-agnostic: feed it any LFE catalog
+ any station list + waveform-cache path, and it returns cross-cycle
families. Margin-wide scaling is data engineering, not new science.

### E.7h Per-path framing — small per-path changes ARE physical

The earlier aggregate result (E.7g) reported the cross-family weighted mean
as if it were a single dv/v, which suppressed path-specific information.
With proper repeating-source CWI, **each (family, station) pair is a
specific physical path** — source patch on the plate interface → station
at the surface — and the dv/v along that path over time is the physical
quantity. Aggregating across paths averages over different physical
quantities.

Reframed with this lens, the high-CC paths (mean CC > 0.85 between year
stacks) give:

| Path        | mean dv/v vs 2010 | std across 2011–2013 | mean CC |
|-------------|------------------:|---------------------:|--------:|
| **fam2-PGC**| **−0.97 %**       | **0.22 %**           | 0.987   |
| fam0-PGC    | +0.36 %           | 0.51 %               | 0.978   |
| fam10-PGC   | −0.38 %           | 0.52 %               | 0.994   |
| fam14-PGC   | +0.15 %           | 0.32 %               | 0.988   |
| fam0-LZB    | +0.07 %           | 0.84 %               | 0.980   |
| fam9-LZB    | −0.11 %           | 0.89 %               | 0.987   |
| fam18-LZB   | −0.18 %           | 0.66 %               | 0.991   |

Two observations:

1. **Some paths are highly stable across years** (fam14-PGC at ~0%, std 0.3%;
   fam0-PGC at +0.4%, std 0.5%). At our noise floor of ~0.3–0.5% per
   measurement, those paths show no detectable inter-ETS dv/v change.
2. **At least one path shows a consistent offset** — fam2-PGC has dv/v
   ≈ −1.0% across **all three** of 2011, 2012, 2013 vs the 2010 reference,
   with std 0.22%. The velocity along that specific source-station path
   has decreased by ~1% relative to 2010 and stayed there through 2013.
   This is a real path-specific signal that is invisible in the aggregate.

The earlier "aggregate consistent with zero" finding (E.7g) is correct as
a network-wide average but misses the underlying per-path structure.
**The proper repeating-source CWI measurement is the per-path dv/v time
series, not the aggregate.**

Figure: `figures/smoke_per_path_dvv_clean.png` shows the 10 trustworthy
paths' dv/v values across 2010–2013. The line for fam2-PGC is the clearest
example of a path-specific consistent offset.

### E.7i Single-station matched-filter on continuous data — proper time series

The previous Phase E.7g/h sections only fetched peak-day ± 1 day data for
2011-2013, which collapsed each non-2010 year into a single dv/v measurement.
That was an observational limitation, not a physical one: LFE patches fire
many times per year, not just on the peak ETS day.

Fixed by fetching full ETS-span continuous data (≥30 days/year for 2011-13)
from IRIS, then running **sliding-window matched filter** on PGC continuous
data with each family template (`src/tremorferometry/matched_filter.py`).
This finds every time a template matches anywhere in continuous data, not
just at Lin's pre-detected events. Single station, single template, every
sample — the literal "same signal at same station over time" test.

**Validation:** stacking 500 random matched-filter detections (CC ≥ 0.7) for
family 0 in 2010 gives a stack that matches the original template at
CC = 0.966. The matches are real LFE family members, not noise coincidences.

**Catalog growth:**

| Family | Detections at CC ≥ 0.7 (PGC, 2005-2013) |
|---|---|
| 0  | 129,260 |
| 2  | 228,526 |
| 4  | 25,509 |
| 9  | 28,595 |
| 10 | 160,777 |
| 14 | 87,097 |
| 18 | 41,080 |
| **Total** | **700,844** |

These come from continuous data scanning, **not** from pre-existing detection
lists. They include detections during ETS bursts and lower-rate firings
between ETS cycles.

**Daily-resolution dv/v across multiple ETS cycles** (figure
`figures/smoke_daily_dvv_pgc.png`):

- 1,897 (family, day) dv/v measurements with ≥20 detections per day
- Reference: pre-Aug 2010 stack at each family
- All paths stay within ±0.3 % across 2010-2013
- Per-family persistent offsets exist (fam 2 ~ −0.2 %) but are tiny vs the
  ~1 % offset we incorrectly inferred earlier from peak-day-only data
- The ETS windows (shaded) show **no clear bulk dv/v excursion** at
  daily resolution

**Conclusion (refined):** with the proper observational density (matched
filter on continuous data → hundreds of LFE detections per day per family),
the dv/v along the southern V.I. plate-interface-to-PGC source-receiver
paths is **stable to ~±0.3 % across the 2010-2013 ETS cycles**. The earlier
"−1 %" finding for fam2-PGC was inflated by using small-N per-year stacks
from peak-day-only data; with daily stacks of hundreds of detections, the
true per-path inter-ETS variability is much smaller.

This is the cleanest LFE-CWI measurement we have. ETS does not produce a
detectable bulk dv/v signal along these paths at our resolution.

### E.7j Wider patch survey — 35 cross-year repeating patches across S. V.I.

User flagged that the original 7 cross-year families all clustered in a
single ~5 km region near (48.85, −123.85). That was an artifact of seeding
discovery from one hotspot. Lin's catalog covers a much wider area with
many distinct hotspots.

Re-ran the discovery with a 6000-detection sample across **all** Lin hotspot
cells in the 6 peak days (2005, 2009, 2010, 2011, 2012, 2013):

- 20 distinct 0.05° hotspot cells across lat 48.85–49.30, lon −124.55 to
  −122.95
- All-pairs network CC at PGC + LZB, threshold CC ≥ 0.75
- After excluding the >200-member transitively-merged super-cluster
  (over-bridging at lower thresholds), **35 cross-year tight repeating
  patches survive**

Patch distribution (from `data/detected_patches.csv` and figure
`figures/smoke_all_patches_map.png`):

- Spread across lat 48.40°–49.02° (~70 km N–S), lon −124.28° to −123.64°
  (~50 km E–W)
- Top patch (id 37): 101 members across 5 years
- Most are 3–5 members each (limited by the 6000-event sample)
- Lin's location precision (~5–10 km per detection) gives apparent
  per-patch std of 0.2–0.5° even for truly co-located events — the
  waveform CC is the authoritative same-source criterion

Cartopy map (`figures/smoke_map_cartopy.png` and
`figures/smoke_all_patches_map.png`) shows the patches concentrated
between PGC, LZB, and NLLB — the area where these stations can sense
LFE coda above noise.

The original 7 families are a subset of these 35 — specifically the
densest sub-clusters at the central hotspot. The wider survey reveals
the full network of patches detectable at PGC + LZB.

### E.7k Extension to 2005-2026 at PGC

Fetched continuous PGC data for May–Nov each year 2014–2026 (the months
where N. Cascadia ETSs cluster). Now have **2,469 PGC day-files spanning
2005-09 through 2026-05** (29 GB total).

Channel timeline at PGC:
- 1999–2017-08: BHZ at 40 Hz
- 2017-08-onward: HHZ at 100 Hz (we resample to 40 Hz for consistency)

Running the matched filter for all 35 patch templates against all 2,469
PGC days → expected detection catalog of order 10^7 events across 21 years.
This is the input to:

- Per-patch detection rate time series (when did each patch fire over the
  21 years?)
- ETS cycle identification from rate (expect ~15 N. Cascadia ETSs 2005-2026)
- Per-patch daily dv/v across the full record
- Spatial dv/v "movie" — each ETS as one frame, 35 patches as 35 pixels

This is the multi-decadal, multi-patch product the project has been
building toward. (See also Section 0 for the fisheye-camera intuition.)

### E.7l Canonical CWI reference at 51 patches (all-time mean LFE)

After PNSN-driven self-detection (E.7i) and strict complete-linkage filtering
turned up **16 new cross-year families post-2014** to add to the 35
Lin-seeded ones, we standardized on the canonical CWI reference: per
(family, station), the reference is the L2-normalized sum of every L2-normed
LFE detection in that family across the record, weighted by daily count.

Phrasing: "how much does the velocity on this path change relative to the
*average propagating LFE* for this family?"

Computation (per family):

```
ref_unnormalized = sum_d (daily_sum_of_L2_normed_cuts[d])  /  total_LFE_count
                 = average LFE waveform across the entire record
day_stack[d]     = sum_d_cuts / n_det[d]
dvv[d]           = stretch_dvv(ref, day_stack[d],
                               t_min=0, t_max=2,  ← direct-S window
                               eps_max=0.02, n_eps=401)
keep if cc_max >= 0.8
```

Result (`data/daily_dvv_51_alltime_ref.csv`): 123,793 daily dv/v
measurements across 50 patches (one new family dropped for insufficient
data); mean stretching CC = 0.993; per-day dv/v 5–95%ile = -0.18% / +0.14%;
per-patch long-term mean ≈ 0 by construction.

Figure `figures/smoke_dvv_51_alltime_ref.png`: per-patch 60-day rolling
medians for both the Lin-seeded 35 and the strict-new 16, plus the
cross-patch median (essentially flat at 0 across 21 years). Per-patch
dv/v stays within ~±0.1% of the long-term mean across 21 years.

### E.7m Long-window image plot of one family (PGC_79)

For the largest single LFE family at PGC (`orig_79`, 872 k detections,
21-year record), built a record-section image:

- Re-cut each detection from continuous data over a **long window
  −3 to +10 s** around the envelope-peak alignment time (bandpass 2–8 Hz,
  L2-norm).
- Sum/normalize per calendar day → 4,775 daily L2-normed stacks.
- Place each daily stack at its **true calendar y-position** (NaN-fill
  missing days) so the image is honest about coverage gaps.

**Coverage caveat caught here.** A first version used
`imshow(stacks, extent=[..., date_first, date_last])` which evenly
distributes 4,775 rows across the date range — combined with
`interpolation='nearest'`, the empty 2006–2009 band got painted with
the nearest observed row, faking continuous coverage. The honest
calendar-grid version shows that **actual coverage at PGC is 63.1%** of
the calendar span — a big gap 2005–2010 and ETS-only intermittent
coverage 2010–2013 (we only fetched waveforms for ETS-active months
back then), then dense 2014–2026.

Figure `figures/smoke_long_window_image_PGC_79.png`. Visually the
direct-S band (0–2 s) is a tight vertical stripe across all years;
faint coda extends out to ~5 s.

### E.7n Coda window beats direct-S — the textbook CWI choice

After noticing in the PGC_79 dv/v time series (`smoke_dvv_PGC_79.png`)
that the direct-S (0–2 s) dv/v measurement showed conspicuous
**negative excursions concentrated during ETS-time LFE bursts**
(distribution skewed asymmetrically negative: median -0.005%, 5–95%ile
-0.16 / -0.004%), redid the same measurement using the early-coda
window 1–3 s (same reference, same CC≥0.8 filter).

| Window | median dv/v | 5–95%ile | mean stretching CC |
|--------|-------------|----------|--------------------|
| 0–2 s (direct pulse) | -0.005% | -0.16 / -0.004% (asymmetric) | 0.997 |
| 1–3 s (early coda)   | +0.004% | -0.10 / +0.09%   (symmetric) | 0.985 |

**The ETS-time "drops" go away in the coda.** Most likely they were
direct-S **waveform-shape changes** (slight sub-source mixture
differences during high-activity bursts within a family), not medium
velocity changes. Coda waves average over many scattering paths and
reject that contamination. This is exactly what classical CWI theory
warns about: stretch the coda, not the direct phase.

**Methodological consequence.** The 51-patch all-time-ref dv/v product
above used the 0–2 s direct-S window for symmetry with template length;
any apparent ETS-time signal in it is suspect and should be re-checked
against a coda-window redo. Default stretching window going forward:
**1–3 s after envelope-peak alignment**, with the long-window data
products (`data/long_window_daily_<TEMPLATE>.npz`) as the source rather
than the short 2-s daily-stack npz.

Figure `figures/smoke_dvv_PGC_79_window_compare.png`: side-by-side
1–3 s vs 0–2 s scatter for PGC_79 makes the contamination obvious.

### E.7o Full 51-patch coda-window redo

Following E.7n at PGC_79, redid the full 51-family dv/v product on the
1–3 s coda window for every patch:

1. The 16 STRICT_* templates didn't have a saved MF detection CSV
   (their raw detection times had only existed in memory during the
   original 0–2 s product build), so reran the matched filter for them
   against all 2005-2026 PGC days
   (`scripts/scan_strict_templates.py` → `data/mf_pgc_strict.csv`,
   32.8 M detections at CC≥0.7).
2. Filtered STRICT detections to CC ∈ [0.8, 1.1] (upper bound rejects
   rare numerical artifacts with vanishing window-std) and concatenated
   with `mf_pgc_2005-2026_cc08.csv` → `data/mf_pgc_all51_cc08.csv`
   (12.1 M detections across 51 templates).
3. Rebuilt long-window daily stacks for all 51 with
   `scripts/build_long_window_daily_all51.py` (loads each PGC day once
   and processes all 51 templates against it; output:
   `data/long_window_daily_all51.npz`, 125 066 (template, day) rows).
4. Stretched every day-stack against its family's all-time-mean LFE
   reference on the 1–3 s coda window via `scripts/dvv_coda_51.py`.

Bug-fix note: the very first build collapsed every detection time to
~1970 because `pd.to_datetime(..., format='mixed')` returned
`datetime64[us]` on this pandas version, and the subsequent
`astype('int64') // 1000` therefore turned out to be a milliseconds-not-microseconds
conversion. Fixed by an explicit cast to `datetime64[ns]` before
integer roundtripping.

Result (`data/daily_dvv_51_coda_1to3.csv`):

| field | value |
|-------|-------|
| measurements | 123 644 daily dv/v values |
| patches | 49 (51 attempted; 2 dropped for sparse data after CC≥0.8 filter) |
| mean stretching CC | 0.977 |
| per-day dv/v median | -0.002% |
| per-day dv/v 5–95%ile | -0.088% / +0.078% (**symmetric** around 0) |
| cross-patch 60-d median range | within ±0.04% across 21 years |

Visible mild ETS-time dips of ~0.02–0.04% in the cross-patch median
around 2017, 2022, 2024 — present in both the Original-35 panel and the
New-Strict-16 panel, suggesting the dips are not template-set artifacts.
But the amplitude is right at the level of normal cross-patch median
fluctuation, so I would not over-interpret them yet.

Compared side-by-side with the previous 0–2 s product
(`figures/smoke_dvv_direct_orig35.png`):

- **Direct 0–2 s**: per-patch tracks routinely reach ±0.1%, cross-patch
  median dips to -0.05% during ETS clusters; 5–95% asymmetric at
  -0.19 / +0.14%.
- **Coda 1–3 s**: per-patch tracks within ±0.05%, cross-patch median
  inside ±0.04%; 5–95% symmetric at -0.09 / +0.08%.

The coda is ~2× tighter and unbiased — consistent with classical CWI:
the medium probe lives in the coda, not the direct phase. The
**canonical PGC-only dv/v product for this project is now
`data/daily_dvv_51_coda_1to3.csv`**; the 0–2 s version is retained
solely for the methodology comparison and should not be cited as a
medium measurement.

Figures:
- `figures/smoke_dvv_51_coda_1to3.png`: the canonical 51-patch coda
  product (two panels: Lin-seeded 35, PNSN-discovered strict 16).
- `figures/smoke_dvv_coda_vs_direct_orig35.png`: coda-only single-panel
  on the 29 Original-35 patches with sufficient data.
- `figures/smoke_dvv_direct_orig35.png`: same single-panel layout but
  on the 0–2 s direct-S product, for the side-by-side comparison.

### E.7p PGC backfill 2005-2013 -- continuous record

The canonical PGC product through E.7o was dense 2014+ but ETS-summer-only
pre-2014 (~150 day-files/yr for 2010-2013 and ~0 for 2005-2009). To
make the dv/v product continuous back to 2005:

1. Fetched IRIS 2010-2014 at 24 workers (the dense, easy years) ->
   ~1,400 new day-files, filling 2010-2013 to 95-100% per year.
2. Fetched NRCAN 2005-2010 at 1 worker (parallel calls return spurious
   404s under load; single-worker confirms the rest are genuine archive
   gaps) -> ~330 new day-files, ~17% coverage per year 2005-2009.
3. Re-ran the matched filter on the 1,460 new day-files only
   (`scripts/scan_new_days_all51.py`, diff disk-days vs existing MF
   catalog) -> +5.3 M cc>=0.8 detections, combined 15.76 M total.
4. Rebuilt long-window daily stacks -> 164,229 (template, day) rows
   (+31% over pre-backfill). Recomputed coda 1-3 s dv/v with the
   all-time-mean LFE reference per family.

Bug-fix in `scan_new_days_all51.py`: `patch_templates.npz` contains
BOTH PGC_* and LZB_* templates (35 each, for the original two-station
discovery in E.7e). The first run loaded ALL of them, scanning the new
PGC waveforms with LZB templates too -- producing ~1.6 M spurious
"detections" that contaminated the combined CSV. Script now whitelists
templates whose name starts with the target station prefix.

Final canonical product (`data/daily_dvv_51_coda_1to3.csv`):

| field | value |
|---|---|
| measurements | 162,791 daily dv/v values (+32% over E.7o) |
| date range | continuous 2005-2026 |
| 2010-2026 coverage | 94-100% per year |
| 2005-2009 coverage | ~17% per year (NRCAN archive limit) |
| patches | 49 of 51 attempted (PGC_66, PGC_104 drop out at cc>=0.80 stretch) |
| mean stretching CC | 0.977 |
| per-day dv/v median | -0.001%, 5-95%ile -0.088 / +0.079% |

Figures: `smoke_dvv_51_coda_1to3.png` (single-panel cross-patch median),
`smoke_long_window_image_PGC_79.png` (record section).

### E.7q Audit of the 51 -- 6 suspect families on stack SNR

Quick auditing the 51 families by the all-time-mean LFE stack's
direct-pulse SNR (peak amplitude in [0, 2 s] / RMS in [-3, -0.5] s)
identifies a cliff:

- **45 healthy families**: SNR > 1, most > 5. Highest is STRICT_599 at
  SNR=143. These show a clear direct-S pulse above pre-pulse noise in
  the all-time stack.
- **6 suspect families**: SNR < 1 in all six. The all-time stack has
  no visible direct-S pulse above background.

The 6 suspects (`PGC_22, PGC_48, PGC_66, PGC_104, PGC_114, PGC_116`)
share a tight signature:
- 64-82 daily stacks (vs 1000-5000 for healthy)
- 12-14 distinct years (vs ~21 for healthy)
- ~42,000 matched-filter detections each (the MF still "finds" them
  because the template shape is what's being matched, but those are
  largely false-positive detections on noise)
- 0-3 daily stacks pass the cc>=0.80 coda stretch filter (vs 200-5900
  for healthy)

Interpretation: these templates were derived from Lin's 0.05-deg
clustered detections at lat/lon bins that don't actually contain a
single coherent repeating source. The template ends up being mostly
band-limited noise that the matched filter then "finds" everywhere.

Net effect on dv/v: collectively contribute <=6 day-measurements out
of 162,791, so they don't pollute the cross-patch median, but they
should be dropped from any "family list" deliverable. Effective real
family count is **45**, not 51. (Note: this is consistent with the
canonical 49 in the dv/v CSV minus the 4 sus families that contribute
1-3 measurements each.)

To filter later, the rule is straightforward: require `template SNR
>= 1.0` (a generous bar) on the all-time-mean stack at PGC.

### E.7r NLLB pipeline -- waveform backfill (multi-station extension)

First step of the §8.2 NLLB workflow: bring NLLB waveforms to parity
with PGC across 2005-2026.

Final NLLB coverage (5,979 day-files, 41 GB):

| Year | files | coverage | notes |
|---|---|---|---|
| 2005-2009 | ~40-67 each | ~10-18% | NRCAN archive limit (similar to PGC) |
| 2010-2018 | 346-366 | 94-100% | IRIS |
| 2019 | 310 | 84% | Persistent IRIS gap at the HHZ-channel transition |
| 2020-2025 | 346-366 | 94-100% | IRIS |
| 2026 | 66 | year-to-date | IRIS |

NLLB channel timeline: BHZ 40 Hz 2003-2017, HHZ 100 Hz 2017-open. The
HHZ transition introduces sub-sample FS jitter (99.999... vs 100.0)
that broke `repeater.py::_load_day_filt::st.merge`. Fixed there the
same way `matched_filter_fast.py::_load_day_filt` was fixed earlier
(resample each trace before merging).

**Process artifact discovered**: large multi-year IRIS fetch jobs
sometimes exit silently mid-run with `wrote 0` and no error -- the
process gets killed by the harness (or IRIS connection state goes
bad) and the bottom of the script's print loop never runs. Mitigation:
do FDSN backfill in **year-by-year subprocess invocations**, each
short enough to complete within the harness window. For years that
still don't fill (HHZ-100Hz era at high concurrency), drop to monthly
chunks and run multiple passes. The `fetch_day` function in
`waveforms.py` skips on-disk files, so passes are idempotent.

Distance check: PGC's 51 families sit 10-85 km from PGC; the same
families are 10-110 km from NLLB. NLLB additionally covers tremor
sources up to ~50 km north of itself (49.5-50 N) that PGC cannot see
at all. Defensible PNSN bbox for NLLB-side discovery is ~80-100 km
around NLLB, i.e., **48.5-50.0 N, -124.7 to -123.0 W**, though we
fetched 47.5-51.5 N to verify the geophysical taper.

### E.7s PNSN catalog extended to 47.5-51.5 N

The original cached PNSN catalog (`catalogs/pnsn_tremor_2014-2026.csv`)
had been fetched with the V.I.-only bbox 47.5-50.0 N. Re-fetched
2010-01..2027-01 over 47.5-51.5 N x -125.5 to -122.0 W
(`catalogs/pnsn_tremor_2010-2026_extn.csv`):

| Lat band | events | notes |
|---|---|---|
| 47-48 N | 41,790 | Olympic Peninsula |
| 48-49 N | 129,425 | S. V.I., Strait of Juan de Fuca |
| 49-50 N | 41,767 | S. V.I. interior |
| **50-51 N** | **24** | Cascadia tremor taper (geophysical limit) |
| 51-52 N | 0 | (geophysical) |

213,006 total events (+25% over original cache, mostly from time-gap
infill, but also confirming the north-of-50 region is essentially
inactive). For southward / margin-wide extensions later, PNSN extends
all the way to ~40 N (Mendocino Triple Junction), with a real activity
gap 46-47 N between WA and OR tremor zones.

### E.7t NLLB family discovery -- methodology consolidation

The pipeline for any new station (locked in conversation):

**Branch A (Lin-seeded, multi-station)**
- Cut envelope-aligned 2-s windows at PGC + NLLB around every Lin OT.
- Network CC = mean of all-pairs PGC-CC and all-pairs NLLB-CC.
- Cluster + cross-year filter.
- Threshold question: single-station 0.80 doesn't transfer cleanly
  to two-station mean (largest V.I. proto-family at PGC max single-CC
  ~0.83 but network mean max ~0.76). Options on the table:
  (i) network mean >= 0.80 (very strict, near-zero yield),
  (ii) per-station >= 0.80 both (strict but equivalent to single-
       station bar applied symmetrically),
  (iii) network mean >= 0.70 (matches the original E.7e finding of
        16 cross-year families at PGC + LZB).
  Parked while we ran Branch B (which uses single-station 0.80
  unambiguously).

**Branch B (PNSN-driven, single-station)**
- Every PNSN tremor window in the bbox -> envelope-peak detection at
  NLLB -> candidate list with inherited tremor lat/lon.
- Bin by 0.05 deg lat/lon -> per-bin all-pairs CC at NLLB ->
  **complete-linkage cluster at single-station CC >= 0.80** -> keep
  clusters with >=3 members across >=3 years.
- This is the exact analog of the procedure that yielded the 16
  STRICT new families at PGC, applied at NLLB.

`scripts/discover_nllb_pnsn_driven.py`. Stage 1 (candidate detection)
produces 1.73 M NLLB envelope peaks across 1,541 days of active PNSN
tremor. Stage 2 (per-bin clustering) is running at the time of this
note.

**Network autocorrelation as TAG, not filter** (key methodological
clarification): hard-applying "must be coherent at both PGC and NLLB"
to Branch B would discard the most interesting class of finds --
sources visible at NLLB but not PGC (e.g., north of NLLB, outside
PGC's reach). Instead, run network CC as a per-family **tag** after
Branch B discovery:
- `two-station-validated`: NLLB family with a coherent PGC stack
  cut at the same source times. Best for cross-station dv/v +
  spatial inversion.
- `NLLB-only-strong`: PGC stack incoherent but NLLB stack has high
  SNR + tight pairwise CC distribution + smooth coda decay. Real
  family that PGC can't see. Single-station NLLB dv/v only.
- `weak/noise`: low NLLB SNR, marginal pair CC, no coda decay.
  Discard.

LFE-vs-noise discriminators applicable post-discovery (in order of
cost-effectiveness):
1. **Tight pairwise CC distribution** (median pair-CC well above the
   0.80 threshold, not just barely passing).
2. **Stack impulsiveness** -- direct-pulse peak / pre-pulse RMS >= 3-5.
3. **Spectral peak in 2-6 Hz** vs flat/lined/microseism-dominated.
4. **Coda decay structure** -- smooth exponential tail vs featureless flat.
5. **Time-of-day uniformity** -- catches cultural / anthropogenic noise.
6. **PGC cross-check tag** (above).
7. **3-component polarization** -- requires fetching N/E channels.

Notably absent from the list: "ETS-concentration" -- LFEs occur
outside ETSs too (background / inter-ETS LFEs, sub-ETS slip events).
The PNSN seed already enforces temporal concentration through its
tremor windows, so an additional ETS-only filter would wrongly drop
real background families.

### E.8 Lessons

1. **Always build the reference from signal-rich data.** A noisy reference
   produces spurious dv/v signals from noise-pattern matching, with the
   sign and magnitude depending on which random noise patterns happen to
   stretch into alignment. The fix is to choose the reference window
   by *signal* (highest n_windows), not by *time* (fixed pre-event).
2. **Symmetric coda is cheap.** Doubling the coda samples by using both
   lag sides cut the SE by ~2× with no extra processing.
3. **Preprocessing matters in the opposite direction from ambient noise.**
   For tremor, raw bandpass beats one-bit / running-mean — the coherent
   direct phase IS the signal, not something to suppress.
4. **Stretch the coda, not the direct phase (E.7n).** Even when LFE
   templates are envelope-aligned and the direct-S pulse is the highest-CC
   feature, dv/v measured on it is contaminated by source-pulse-shape
   changes (e.g. sub-source mixture during ETS bursts). The coda — a few
   seconds *after* the direct-S — is the medium probe. Direct-pulse
   stretching can manufacture an ETS-correlated dv/v "signal" that is not
   in the coda.
5. **An image plot stretched on row index can hide gaps.** When plotting
   day-stacks vs date with imshow + a date extent, the rows are spread
   uniformly through the y-range and nearest-neighbor interpolation
   paints empty time intervals with whichever observed row is closest.
   Always either place rows at their true calendar index with NaN gaps,
   or use pcolormesh with explicit y-coords, before reading coverage off
   the plot.
6. **The 0.80 STRICT threshold is single-station, not network** (E.7t).
   Applying single-station CC>=0.80 to a two-station mean is
   inconsistent: the network mean is bounded above by the lower of the
   two per-station CCs, so requiring mean>=0.80 effectively requires
   both stations to individually exceed 0.80, which is much stricter
   than the per-station bar. Pick a consistent definition (e.g.
   per-station 0.80 at all stations) before threshold-tuning.
7. **Network autocorrelation is a tag, not a hard filter** (E.7t).
   Requiring coherence at all stations discards the most scientifically
   interesting class of finds -- single-station-visible sources at the
   new station. Use cross-station check to tag families, not to
   exclude them. The original PGC STRICT 16 were found single-station
   without LZB validation.
8. **patch_templates.npz holds templates for BOTH PGC and LZB**
   (E.7p). Any matched-filter wrapper that loads it indiscriminately
   will scan a single-station's waveforms with the wrong-station
   templates and produce contaminated detections. Whitelist by
   station prefix.
9. **FDSN backfill: do year-by-year subprocess calls, not one big run**
   (E.7r). Multi-year IRIS pulls sometimes exit silently mid-run
   (process killed, partial save, no error). Smaller per-process work
   units + idempotent skip-on-disk is much more reliable. Drop to
   monthly chunks for problem years (especially HHZ 100 Hz era).
10. **NRCAN single-worker FDSN** (E.7r, E.7p). Parallel calls return
    spurious 404s under load. For pre-2010 archives that aren't on
    IRIS, fetch with `--workers 1` to confirm what's actually missing
    vs throttled.

## 6. What's next

**Updated after Phase E (tremor-CC arc).** Current headline measurement is
the tremor-CC v2 (symmetric coda, signal-rich reference) which gives the
tight null bound. Highest-value follow-ups:

1. **Refine tremor-CC SNR**: per-window weight by tremor energy or PNSN
   amplitude, exclude marginal tremor windows, try 3–5 Hz instead of 2–5,
   try sub-bbox restrictions on which tremor windows feed each pair's
   CC stack (some pairs are only well-illuminated by tremor in part of
   the bbox).
2. **Hydrology / loading control**: overlay precipitation and GPS-derived
   surface loading. The CWI signal we *can* see at this sensitivity might
   be dominated by seasonal hydrology, not the ETS. A null result against
   the ETS is more defensible if we show we *would* see other known
   sources at this sensitivity.
3. **Sub-network analysis**: split into V.I.-only pairs vs cross-V.I./OP
   pairs. Cross-network pairs sample longer paths through more of the
   plate-interface volume — if there's any signal, it should show up
   there preferentially.
4. **Self-detect LFE families from continuous data** (Brown-Beroza-Shelly 2008
   network autocorrelation, or REDPy with a tremor-appropriate trigger).
   Independent of tremor-CC, this would let us redo the LFE-CWI properly
   with real repeating sources. Substantial engineering.
5. **Other ETS events** — repeat tremor-CC on the 2009-05, 2011, 2013
   ETSs to see if the bound changes. Lin's catalog covers 2005–2017 so
   plenty of candidates.
6. **Margin-wide** — ingest Sweet (2019), Ducellier (2022), Plourde (2015)
   for the rest of Cascadia. Tremor-CC trivially extends — just plug in
   new station lists and bboxes.
7. **GPU template matching** — fast-matched-filter is built and ready on
   the L40S; useful when (4) above is launched and we need to match the
   discovered templates against years of continuous data.

---

## 8. Multi-station extension and margin-wide plan

The PGC-only product (Sections E.7l–o) measures dv/v on N=51 paths
{ source patch → PGC } across 21 years. Extending this is two related
problems:

A. **Second station near the same sources** (e.g., NLLB) → confirms each
   PGC family with a fully independent receiver and adds a second path
   per family. Required for cross-validation of the ETS-time signal and
   for any tomographic combination of paths.
B. **Margin-wide extension** (Olympic Peninsula, central Oregon,
   N California) → discovers new LFE families in regions Lin's catalog
   does not cover.

Both reuse the same discovery building blocks; the only thing that
varies region-to-region is which **seed catalog** is available.

### 8.1 Unified discovery recipe

Seed selection by region:

| Region | Lin (2023) | PNSN tremor catalog | Local broadband |
|--------|-----------|---------------------|-----------------|
| S. V.I. 2005–2017 | yes | yes | yes |
| S. V.I. 2017–2026 | no  | yes | yes |
| Olympic Peninsula | no  | yes | yes |
| SW Washington / N Cascades | no | yes | yes |
| Central Oregon | no | yes | yes |
| N California | partial (Plourde 2015, Ducellier 2022) | yes | yes |

Regardless of where Lin is available, the pipeline is the same five-step
chain; we just feed it whatever seed exists for that region.

**Step 1 — Seed candidate detection times.**
- If Lin is available for this region: take Lin's detections in the
  region's bounding box as candidates (`src/tremorferometry/lin_catalog.py`).
- Otherwise: PNSN tremor windows (`pnsn.py` / `00b_fetch_pnsn_tremor.py`)
  → envelope-peak detection at one or two local broadband stations
  within each window
  (`src/tremorferometry/detect.py::envelope_peaks_in_windows`,
  SNR≥3 default, min separation 6 s).

**Step 2 — Cut envelope-aligned 2-second windows.**
At every regional station, cut a 2-s window centered on the Hilbert
envelope peak in the expected direct-S arrival window for that
station–source geometry. Bandpass 2–8 Hz, L2-normalize.
(`repeater.py::cut_aligned_window`,
 `repeater.py::cut_all_detections`.)

**Step 3 — Shelly-Beroza network autocorrelation.**
Compute all-pairs max-shifted CC at each regional station; average across
stations (only pairs where both detections have valid data at that
station). This is the Brown-Beroza-Shelly 2008 + Bostock 2012/2015
recipe. (`repeater.py::all_pairs_cc_max_shifted`,
`network_cc_all_pairs`, +/- 20-sample shift tolerance.)

**Step 4 — Cluster into families.**
Complete-linkage clustering on the (network CC ≥ 0.80) graph (i.e.
every member-pair within a family must have CC ≥ 0.80, not just be
transitively connected). Require ≥3 members across ≥3 years
(cross-year repeater filter) — this is what makes the family usable for
multi-decadal dv/v rather than a one-off cluster.
(`repeater.py::cluster_matches` with complete-linkage threshold.)

Note: the original E.7e discovery used τ=0.70 with transitive closure;
that found 16 cross-year families but is now superseded. The canonical
discovery threshold for this project is **0.80 with complete linkage**
(matches the criterion used for the 16 strict new families at PGC and
for the cc≥0.8 filter on every downstream matched-filter detection).

**Step 5 — Densify with matched filter (and optionally PNSN growth).**
Build each family's per-station template (stack of its member
waveforms, L2-normed). Run the batched-fast matched filter
(`matched_filter_fast.py`) at every regional station across all
continuous days to recover detections the seed missed. Optionally feed
the matched-filter output back as a new candidate set and rerun
clustering once to capture sub-family structure or near-relatives
(this is the "PNSN-grown" path that yielded the 16 strict at PGC).

Output per region:
- `families_<region>.csv` (cluster id, lat, lon, n_members, year span)
- `templates_<region>.npz` (per-(family, station) 2-s waveform)
- `mf_<region>_<station>.csv` (matched-filter detection time series)

### 8.2 NLLB-specific pipeline (Lin available)

Concrete recipe for the V.I. second-station step:

1. **Backfill NLLB waveforms** 2005-2026 with
   `scripts/backfill_pgc_2005-2013.py` (refactored to accept a
   `--station` argument). IRIS for 2010+, NRCAN for earlier.
   We already have 161 ETS-summer days; need ~5800 more.
2. **Seed**: Lin detections in the V.I. bounding box (same set we used
   for the 35 PGC originals — Lin's catalog is location-tagged, not
   station-tagged, so the same detection-times feed every station).
3. **Step 2** of §8.1 at PGC + NLLB jointly: cut envelope-aligned
   2-s windows at each station around each Lin detection. The
   envelope-peak search window at NLLB is offset from PGC by
   (NLLB_travel_time − PGC_travel_time); calibrate empirically using
   ~100 high-CC Lin detections per family from the PGC product.
4. **Step 3** at PGC + NLLB: network CC across both stations.
5. **Step 4**: complete-linkage cluster at 0.80 with the cross-year
   filter. Expected result: a re-discovery of most of the PGC 35
   (those visible at NLLB) plus possibly some NLLB-strong families
   PGC missed.
6. **Spatially match NLLB-discovered families to PGC families**: same
   cluster centroid within ~5 km ⇒ same physical source patch.
   Tag PGC-only, NLLB-only, and both.
7. **Step 5 — densify**: matched filter at NLLB with the new NLLB
   templates → NLLB detection time series for each family.
8. **PGC-pipeline parallel at NLLB**: build long-window daily stacks
   (`build_long_window_daily_all51.py` parameterized on station) and
   compute coda 1-3 s dv/v (`dvv_coda_51.py`).

Comparison product: per family, two dv/v series (PGC path, NLLB path).
Real ETS-time medium changes appear at both; source-pulse contamination
or station-specific noise does not.

### 8.3 Margin-wide extension (no Lin)

For each new region (Olympic, central OR, N CA, ...), pick a regional
sub-network of broadband stations within ~100 km of the target tremor
patches. Network CC needs stations within that range — adding faraway
stations (e.g., PGC) to an Olympic-region network is counterproductive
because the LFE amplitude is below noise.

Suggested regional sub-networks:

| Region | Sub-network |
|--------|-------------|
| S. V.I. | PGC, LZB (pre-2014), NLLB, SNB |
| Olympic Peninsula | UW.OSD, UW.SQM, UW.JCW, UW.HOOD |
| SW Washington / N Cascades | UW.GMW, UW.LRIV, UW.PASS, UW.OHW |
| Central Oregon | CC + UO stations near Three Sisters / Mt Jefferson |
| N California | BK.WDC, BK.HUMO, BK.HUMP |

Recipe per region:
1. Backfill 2005–2026 waveforms for the regional sub-network (FDSN
   IRIS+NRCAN+NCEDC depending on operator).
2. Skip "seed = Lin"; go directly to PNSN-driven envelope-peak detection
   at one or two seed stations in the sub-network (the
   highest-SNR ones).
3. Step 2 of §8.1 across the sub-network.
4. Step 3 across the sub-network.
5. Step 4 with complete-linkage 0.80, cross-year filter.
6. Step 5 — densify and (optionally) re-seed PNSN with the new
   templates to catch near-relatives.

Output stitching: the **margin-wide LFE catalog** is the union of
per-region families. dv/v is measured **per (family, station) within a
region**, not across regions — there's no station with SNR for sources
across the whole margin, so no global dv/v measurement. The final
deliverable is a margin-wide dv/v map composed of regional per-path
measurements.

### 8.4 Open methodological choices for §8

- **Spatial-match tolerance** between PGC-family and NLLB-family
  centroids (currently a guess of ~5 km). Calibrate against Lin's
  per-detection scatter for a known patch.
- **Envelope-peak-search window offset between stations**. Either
  pick from a 1-D V.I. velocity model (Bostock/Cassidy) or learn
  empirically from initial high-CC matches.
- **Minimum visibility threshold** before declaring a PGC family
  "not seen at NLLB" — a family with weak NLLB SNR may still be
  there but below the network CC ≥ 0.80 discovery threshold (or
  may pass discovery but fail the per-detection CC ≥ 0.80 filter
  used downstream).
- **Sub-network composition** when stations come online/offline
  mid-record (e.g., LZB cuts off 2014-06). One option is to allow
  variable sub-networks per epoch; another is to fix the sub-network
  per region for the full 21 years and lose the LZB benefit pre-2014.

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
