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
