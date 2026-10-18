# PLAN — tremor vs SSE slip: relative position through 2010-2026 (drafted 2026-09-10; user: no dv/v yet; execution started 2026-09-10)

## 0. Question and what "new" would look like
Tremor (PNSN) and the Gualandi (2025) vbICA slip model locate the same process independently. Their RELATIVE geometry
(across-strike offset, depth edges, along-strike lead/lag, slip-per-tremor) and its EVOLUTION over 16 years has not been
measured with a continuous slip model. Candidate findings: (a) a secular drift of the tremor–slip across-strike offset or of
either depth edge; (b) a change in tremor efficiency (slip per detection) after particular events or after ~2015; (c) a
systematic lead/lag of the slip front relative to the tremor front that changes across cycles; (d) the ETS slip centroid
wandering across cycles while the tremor band stays fixed (or vice versa).

## 1. Data and resolution limits (set before looking)
- Tremor: catalogs/pnsn_tremor_cascadia_full.csv, 748k detections 2010-01..2026-05, epicentres only (depth column is a
  placeholder). Location error ~5-10 km; detection completeness grows with the network (esp. south, and step at 2019+).
  -> normalise counts per band per year when comparing epochs; never compare raw rates across years.
- Slip: 742-triangle 0.25° mesh (~28 km), fixed patterns x 13 amplitudes. A slip CENTROID can only move by re-weighting
  13 patterns -> within-event "migration" of slip is a mixture effect; cross-event drift of the centroid is meaningful only
  above ~15 km (half a triangle) and must be tested against the mirror catalog.
- Frame: slip_thrust_inc30 is slip within (t-30, t]; tremor is binned the same way for band/epoch work; for FRONTS use the
  raw daily slip rate (d/dt of A_raw x Ld, 7-d smoothed) and daily tremor centroids.
- Depth: tremor depth = Slab2 interface depth at the epicentre (finite for 99 %); slip depth = slip-weighted mean of the
  triangle depths. Both are "on-interface" by assumption.
- Coordinates: local along-strike / across-strike from the Slab2 depth gradient at each point (NOT E-W / N-S; the strike
  rotates from ~N in Oregon to ~NW off Vancouver Island).
- Scope: north (>=46N) and central (44-46N) confirmed patches only; the south is at the model's noise floor (tremor-confirmed
  23 % vs 4 % chance). Mw>=6.0 confirmed patches: ~30 north + 6 central. Small n -> trends need the per-event null.

## 2. Quantities (exact definitions)
For each confirmed patch p (DUR window [ts, te], footprint F):
  Q1 across-strike offset  = (tremor centroid during DUR) - (slip-weighted centroid of raw slip accumulated ts-30..te), in km,
     projected on the local down-dip direction (+ = tremor down-dip of slip).
  Q2 along-strike offset   = same, projected along strike (+ = north).
  Q3 depth offset          = median Slab2 depth of DUR tremor - slip-weighted interface depth (km).
  Q4 edges                 = tremor 10th/90th percentile depth; slip up-dip edge = shallowest triangle with accumulated slip
                             >= 3 mm; down-dip edge likewise (km).
  Q5 efficiency            = accumulated slip (potency, km^3) / tremor detections during DUR on F; also per triangle per epoch.
  Q6 fronts (north/central Mw>=6.2 only): daily tremor centroid along-strike position vs daily slip-rate centroid; lead/lag
     by cross-correlation of the two tracks; propagation speeds (Theil-Sen km/d) of each.
Per year y and segment: tremor depth distribution (p10/p50/p90), slip-weighted depth, band-normalised tremor rate.
Per triangle and epoch (2010-14 / 2015-19 / 2020-26): tremor detections, accumulated thrust slip, ratio.

## 3. Hypotheses (pre-registered)
  H1 tremor sits down-dip of the slip maximum (literature: tremor ~ at/below the slip peak); Q1 > 0, stable.
  H2 Q1 / Q3 / edges drift over 2010-2026 (secular deepening/shallowing) — Theil-Sen slope vs per-event null.
  H3 efficiency Q5 changes with time or after the largest events (2016, 2017, 2024) — epoch maps and per-event series.
  H4 slip front lags tremor front (model needs slip to accumulate) by a lag that is constant across cycles; a CHANGING lag
     would be new.
  H5 the ETS slip centroid drifts across cycles at a site while the tremor centroid does not.

## 4. Controls and nulls
- Mirror (back-slip) catalog: same Q1-Q5 -> distribution of offsets for non-events.
- Time-shifted tremor (rigid shift >= 1 yr): offsets/efficiency under no association.
- Network-growth control: recompute everything with tremor counts normalised per band-year; report the 2019 boundary.
- Fixed-pattern control: synthesise slip with a single pattern and constant amplitude history -> what centroid "motion" the
  13-pattern mixture alone can produce (upper bound on artefactual drift).
- Multiple testing: every trend reported with FWER over the quantities tested in that block.

## 5. Steps and decision points
  S1 Co-registration table (Q1-Q5 per confirmed patch, both catalogs) + per-event figure vs time by segment.   [~30 min]
     DECISION: if mirror offsets are as tight as thrust offsets, the offsets are mesh geometry, stop the geometry line.
  S2 Depth occupancy per year (Q4 + tremor depth percentiles) with normalised rates; trends + null.                 [~20 min]
  S3 Efficiency maps per epoch + change maps; per-event efficiency series.                                           [~20 min]
     DECISION: an efficiency change confined to the network-growth boundary is a detection artefact -> drop.
  S4 Fronts (Q6) for the ~20 north/central Mw>=6.2 events; lead/lag and speeds vs cycle.                             [~45 min]
  S5 Site recurrence: tremor-episode vs slip-patch recurrence intervals per site, phase drift.                       [~15 min]
  S6 (model-only checks) inter-ETS slip: modeled slip between ETS events per site and whether tremor confirms it; the southern
     +-5 mm oscillation: spectrum, seasonal phase, coherence across bands (station-driven vs random).                 [~20 min]
  (dv/v bridge DEFERRED by user 2026-09-10: 'lets not add anything with dv/v yet')
Deliverables: res_catalog_sse/tremor_slip_geometry_{events,years,epochs,fronts,sites}.csv; figures/sse_tremor_slip_geometry_*.png;
ledger in notes/2026-09-10_Notes.md; a results note if anything survives the nulls.

## 6. Known pitfalls to avoid
- Union footprints of chained events bias centroids north (already seen: N-S offset +21 km in the north) -> use slip-weighted
  centroids of accumulated slip, not footprint centroids, and per-day tracks for fronts.
- Tremor onset thresholds trigger on precursory bursts -> report lags with two thresholds.
- The 30-d trailing frame shifts every slip time by 0-30 d; never compare tremor daily timing with inc30 timing.
- Do not read absolute stress/slip amplitudes; only relative geometry and ratios.

## Round 2 (2026-09-10, user: answer everything; push for mechanism)
  S7 occupancy + terminations (Q5, Q17): per triangle x epoch classify slip-with-tremor / slip-without / tremor-without-slip
     (completeness-normalised tremor-days), persistence across epochs; per event along-strike extents of tremor vs slip and the
     termination latitudes -> persistent barriers 2010-2026?
  S8 big events + predictability (Q7, Q18): before/after 2016-01, 2017-03, 2024-09 at their sites (offset, depth, efficiency,
     recurrence, event size); time- vs slip-predictable recurrence per site using catalog Mw AND tremor counts as size proxies.
  S9 within-event mechanics (Q14, Q13, Q2): tremor depth migration within events with a null; tremor rate vs slip rate scaling
     (exponent) per event/segment/epoch from the daily tracks; up-dip/down-dip edges referenced to concurrent margin-wide tremor.

## Round 3 (2026-09-10, user: "find a publishable result") — settle the three candidates
  S10 literature check (WebSearch agent): novelty of (a) sub-linear tremor-rate vs slip-rate scaling constant over 16 yr,
      (b) persistent central-Oregon tremor without resolvable slip / weak co-location, (c) within-event tremor deepening,
      (d) 16-yr stationarity bounds of the tremor-slip system; plus what Gualandi (2025) itself showed vs tremor.
  S11 scaling robustness: tremor ENERGY and DURATION (populated 2017+; 64.5 % of rows), unique-location counts, unsmoothed
      slip rate, saturation tests -> is the exponent < 1 a detection artefact?
  S12 central-Oregon representability: can the 13 fixed vbICA patterns represent a localized slip patch at 44.3-45.4N, 27-38 km?
      Projection of unit patches onto the pattern span (map), GNSS station density vs latitude from the solution files, misfit
      by latitude -> model limitation vs physical decoupling.

## Round 4 (2026-09-10, after the S10 literature check) — make the top candidates claimable
  S11b model-free scaling: exponent of tremor count/energy vs RAW GNSS displacement rate (no slip patterns) and vs the ICA
       reconstruction; agreement with the vbICA-based 0.40 / 0.64 removes the last blocker on the scaling result.
  S9b  reconcile within-event deepening with Wech & Creager 2011 (up-dip migration): onset ramp vs peak vs tail; their
       down-dip-distance statistic; spatial control on the intensity-depth relation.
  S12b temporal loophole for central Oregon (component time functions vs central tremor; data-space check at nearest stations).
  S13  test Meng et al. 2026's quantities (nucleation latitude, propagation direction vs time) from a model-free tremor-event
       catalog with completeness controls (bright-only, >= 8 stations, Feb-2018 algorithm step); power analysis of the S1 bounds.
  S10b literature: predictability candidate; obtain Bartlow 2020 data; central-Oregon statements in Brudzinski & Allen 2007 and
       Schmidt & Gao 2010.
  S14  LFE-based exponent: the project's own LFE detections (data/daily_dvv_<STA>_Z_2to4.csv n_stack = 30-d family counts; data/mf_*
       per-detection catalogs for daily counts in event windows) as a third independent tremor measure vs slip rate.

## Round 5 (2026-09-10) — paper readiness
  S15 (A) amplitude snapshot check: rebuild daily amplitudes from the 2026-08-18 solution, confirm the 27 confirmed events and the
       count exponent are unchanged; (B) third tremor measure from the project's picker-ranked candidates (data/<sta>_cand_filtered.parquet,
       envelope peaks inside PNSN tremor windows, ~195 stations): exponent vs slip rate, paired against PNSN counts.
  Next: Frank & Brodsky aggregation match; figure set; stationarity power table.
  S16 (2026-09-10): Frank & Brodsky aggregation match (median per-detection amplitude vs slip rate, 2017+ events) + the paper figure set
       (Figs 1-4, Tables 1-2, caption drafts).
  S16 DONE 2026-09-10 (resumed after the 21:03 pod crash; PART A from the crashed run's CSVs, PART B new): F&B match ->
       results/s16_fb_summary.md (verdict: reproduces in sign/sub-linearity; value matches 1/3.17 only at matched granularity
       and estimator — daily forward 0.121, 7-d 0.227, cross-event 0.296, RMA 0.271, F&B-direction beta 1.01); figure set ->
       figures/s16_paper_fig1..4.png, results/s16_table1.csv, s16_table2.csv, s16_captions.md; daily tracks (never saved by
       S9/S11/S16) regenerated + validated to 8e-17 -> scripts/s16_daily_tracks.py, results/s16_daily_tracks.csv,
       s16_track_validation.csv. Ledger: NOTES.md "2026-09-10 S16 (resumed after crash)"; results: RESULTS.md R13/R14.
