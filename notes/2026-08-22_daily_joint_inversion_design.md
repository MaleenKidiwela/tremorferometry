# 2026-08-22 — Daily (30-day rolling) 4-D joint inversion: fixed-effects design (PRE-BUILD, for Merlin audit)

Goal: replace the annual-window inversion with one that consumes the DAILY 30-day-trailing-stack dv/v
series (exactly what the HTML map shows, gaps preserved), inverting for m(x,t) = interface δβ/β movie +
per-station-day site terms, with per-pair reference offsets treated as explicit unknowns.

Prereq applied 2026-08-22: the `G = -d["G"]` sign flip in `invert_dvv_4d.py` was removed (data column is
eps = -(dv/v), forward-tested; G.npz stores -K; eps = G_stored·m as stored). SITES/idx_data/closure_mean
negated at save so all outputs are true velocity-change sign. Results on disk not yet re-run.

## The user's formulation (verbatim math)

Write the field as m(x,t); pair p has kernel K_p(x), Σ_x K_p(x) = 1. A single measurement:

    d_p(t) = K_pᵀ m(·,t) + c_p + ε_p(t)

c_p exists because dv/v is measured against a reference stack over window R_p:
c_p = −K_pᵀ m̄_{R_p}. It becomes a free parameter because pairs have different usable spans (the
production reference is each pair's ALL-TIME mean, `dvv_roll30cal.py` ref=Rf.mean(0)) and because of
per-pair processing offsets (response changes etc.).

Augmented system: d = G m + S c + ε, S[(p,t),q] = δ_pq. Null space: any static field f with
m(·,t)=f, c_p=−K_pᵀf gives Gm+Sc=0 → M-dimensional. Gauge: (1/|R|) Σ_{t∈R} m(x,t) = 0 ∀x supplies
exactly M constraints (the reference window IS the gauge).

Eliminating c is a rank-P update (SᵀWS diagonal):
    ĉ = (SᵀWS)⁻¹SᵀW(d − Gm̂);   (GᵀW̃G + Λ)m̂ = GᵀW̃d,  W̃ = W − WS(SᵀWS)⁻¹SᵀW
= the panel fixed-effects "within" estimator. Keeping c explicit gives its posterior covariance and
handles short pairs gracefully (loose offsets instead of noisy sample-mean demeans).

## Key insight (why this beats the current demean)

The CURRENT code demeans d only (per pair, 2019–24 common window else own mean) and solves each window
independently with raw K_p. The within-equivalence requires demeaning BOTH d and G per pair over its
observed epochs — which couples epochs. Demean-d-only ≡ asserting K_pᵀm̄_{R_p}=0 PER PAIR (thousands of
gauge constraints, not one M-dim gauge) → any real slow drift of the field near pair p is removed by its
own mean. This is the common-mode-pre-subtraction slow-field killer, in time. The explicit-c
formulation preserves the slow differential field up to ONE static field (the best a reference-based
measurement can do).

## Proposed build

1. **Assembler** (`assemble_res_catalog.py --daily` or sibling): daily pair-table from the per-station
   CSVs (preference `_unif05` > `_eps05` > original, same as map), certified families only, one row per
   (pair, date), NO aggregation/gap-filling. Keep cc_max, n_stack, source-eps-range; flag railed rows
   (|eps| at its file's search edge). ~25M rows → pair_days.parquet.
2. **Solver** (`invert_dvv_4d_daily.py`, sibling; annual script untouched):
   - unknowns: m(cells × days), s_j(t) (station × day site terms), c_p (per pair×era offsets;
     era boundaries from era_table.csv, ±35 d straddle drop kept).
   - solve by block coordinate descent: (i) per-day solves for m,s given c (existing machinery,
     data d−Sc); (ii) closed-form c update; iterate to convergence. (Full sparse CG as fallback /
     if temporal smoothing on m is added later.)
   - gauges: mean_{t∈R} m(x,t)=0 per cell; mean_{t∈R} s_j(t)=0 per station (removes the c_p-vs-
     constant-site null space, one dim per station). R = 2019–2024 (densely illuminated); cells with no
     illumination in R flagged as floating-level.
   - NO deseason: per-DAY site terms absorb station-common seasonality frame-by-frame; cell-coherent
     (mostly shallow-leak) seasonal wiggle stays in the movie, separated at interpretation.
   - weights: per-datum σ = pair residual scatter (no √n, no N_eff within a day); λ re-oracled at
     single-day noise (NOT transplanted from annual 0.433).
   - gates: MIN_ROWS/MIN_SITES per day → NaN frame (gaps stay gaps).
   - output: MODEL (days × cells), SITES (days × stations), C (+ its σ), per-frame σ_m.
3. **Deferred**: significance machinery (year-scramble meaningless here; per-pair circular-shift nulls
   for episodic features, split-network coherence for slow field), temporal smoothing on m, mirror v2.

## Open questions for Merlin

A. Formulation/gauges correct & complete? Any additional null spaces (limited illumination × free c)?
B. c_p vs per-day s_j(t): is the per-station gauge on s sufficient, or does the s×c interplay need more
   (e.g. can time-varying site terms + offsets conspire to absorb real slow field)?
C. Dropping deseason given per-day site terms — any leak path that argues for keeping it?
D. Block coordinate descent convergence on this structure (weak coupling claim) — safe, or CG needed?
E. Weights/λ at daily noise; the 30× oversampling (consecutive frames share 29/30 of data) — traps?
F. Railed-row policy (drop vs keep) and the mixed eps-range input (unif05/eps05/clipped across stations).
G. Anything that invalidates comparison against the annual gated product.

---

## MERLIN AUDIT RESULT (round 1, 2026-08-22)

Formulation sound; BCD = exact alternating minimization, safe. Core "demean kills slow field" claim
OVERSTATED: 91% of pairs use the common 2019-24 window (single-gauge-equivalent up to ~0.11% per-pair
offset-sampling noise); only 603/6726 own-mean-fallback pairs are genuinely slow-field-killed. Expected
outcome: new slow field SIMILAR to annual = success, not failure (recorded to prevent later "fixing").

BLOCKERS: B1 offset unit = FAMILY not (cell,station) (census below); B2 gauge site terms per
(station,ERA) not per station (era x per-day-site null space); B3 drop railed rows per-file edge +
cycle-skip gate + robust MAD gate; B4 per-epoch transmission T(y,patch) normalization (roster growth
13->188 makes fixed-lambda response time-varying -> fake "accelerating strain").
IMPORTANT: I1 R-gauge as exact post-hoc reparameterization every sweep; I2 reduced-problem BCD-vs-direct
validation; I3 gated-out days excluded from c updates too; I4 seasonal x sampling rectification -> fit
harmonics JOINTLY with slow spline at curve stage + deseason-on variant; I5 lambda oracle with
station-common noise + re-oracle LAM_RATIO; I6 30-d trailing = 15-d centroid lag, gap jitter fakes
propagation for episodic features (movie resolution = 30 d hard); I7 annual-comparison protocol (sign
negation, same gauge, same aggregation, vs ERA_SPLIT variant; bar r~0.8 not 0.99).

## FAMILY-OFFSET CENSUS (Merlin step 1, run 2026-08-22)

2741 (cell,station) pairs / 181 stations. Within-cell spread of per-family mean eps (%):
std p25 0.009 | MEDIAN 0.024 | p75 0.052 | p90 0.105; span median 0.056, p90 0.248.
By generation: orig 0.016, eps05 0.035, unif05 0.177 (selection-confounded: unif05 = the 10 stations
chosen FOR rail contamination). Median passes the <0.05 bar; tail fails. Sent back to Merlin for the
pair-vs-family-vs-hybrid ruling. csv: res_catalog/family_offset_census.csv

## MERLIN ROUND-2 RULING (B1, 2026-08-22)

**(b) FAMILY-level c fleet-wide.** Median-passes-bar is the wrong summary (tail ~25% of pairs at/above
the 0.043% class, concentrated at high-scrutiny stations); hybrid REJECTED (data-dependent estimator
switch correlated with geography = uninterpretable seam; census statistic itself noisy ~0.03-0.1%);
family-level over-parameterizes harmlessly where offsets are equal. Guard: >=30 rows per family-era,
else DROP those rows (loose offset would soak field signal). unif05 0.177 NOT yet the failure signature
(selection confound + orig-clipping shrinks means toward zero = un-censoring explains part); decisive
PAIRED check during build: for the 10 unif05 stations, same within-cell spread from earlier-generation
file on COMMON (family,date) rows excluding rows railed in either -> ~2x match = proceed; 5-10x = STOP,
generations disagree on identical data. Assembler spec settled: family-day rows, c per family x era,
exact asymmetric era straddle [b, b+30d), B3 gates as drops.

## BUILD LOG (2026-08-22)

- FAMILY-DAY assembler (`assemble_res_catalog_daily.py`, rewritten to family granularity): 189 stations,
  32.4M family-days, 14,630 family-eras, 11,947 families, 872 cells, span 2009-01..2026-07. Drops:
  rail 87k, cycle-skip 61k, outlier(8MAD) 22k, era-straddle[b,b+30d) 39k, fam-era<30 1.5k. Family
  sigma (robust MAD) median 0.296%. -> res_catalog/family_days.parquet + families_daily.csv
- PAIRED GENERATION CHECK (Merlin's decisive B1 follow-up): on common non-railed (family,date) rows the
  unif05 and orig files are BIT-IDENTICAL (row diff 0.0000%, identical fraction 1.00, spread ratio 1.00
  all 10 stations) — same stretch grid, only search range differs. unif05 census spread 0.177% = station
  quality + un-censoring. Generations CONSISTENT -> PROCEED. (scratchpad/paired_generation_check.py)
- Solver `invert_dvv_4d_daily.py` written (BCD + per-sweep exact gauge (I1/B2) + closed-form family-era
  offsets + iid+station-common lambda oracle (I5b) + I3 day-gating + true-sign outputs (I7a)).
- Validator `validate_daily_bcd.py` (step 3, BCD vs direct sparse joint solve, 3 mo x 30 stations)
  PASSED 07:04 (max|dM| 9.7e-4% << 0.1 sigma_m, monotone objective; threads capped 16 after a
  176-thread BLAS thrash).
- CRASH LOG: production solver launch #1 OOM-killed the pod (212 GB einsum in predict(); fixed to
  per-day loop 07:19, `run_capped.sh` 140 GB ulimit wrapper written 07:23). Launch #2 (07:31, capped)
  made the pod laggy -> user manually restarted 07:43: NOT memory (flat 11 GB) or disk (79 TB free) —
  unpinned BLAS threads (pod advertises 176 CPUs, capped at 32). Fix: OMP/MKL/OPENBLAS/NUMEXPR=1 set
  in-script before numpy import. Launch #3 07:50: run_capped 140 GB, --workers 24 (8 CPUs headroom),
  nlwp=4, load ~1. Log `logs/invert_daily.log`.
- **PRODUCTION RUN v1 COMPLETE 2026-08-22 ~09:00** (~70 min): lam_f oracle 0.215 (100-km recovery
  0.58 at single-day noise — NOT the annual 0.433). 30 sweeps, objective monotone throughout
  (5.139e7 -> 5.121e7, final rel 2.0e-7). max|dc|/sig_c plateaued EXACTLY at 2.320 for sweeps 13-29:
  diagnosed as gauge shuttle, not non-convergence (ridge-optimal s has nonzero era-mean; gauge()
  moves it into c every sweep, so cn-vs-post-gauge-c contains a constant b; obj contraction 0.78/sweep
  is inconsistent with a real mode that slow). FIXED in script: dC now compares successive un-gauged
  optima (prev_cn). Outputs `res_catalog/inversion_4d_daily.npz` + `daily_offsets.csv` (13,826
  family-eras; 5,013 with |c|>2sig — offsets doing real work; WISH/LCCR era steps clearly captured,
  e.g. LCCR era0 +0.38% vs era1 ~0 = the instrument-leak blocker CLOSED). Field sanity: per-frame
  posterior sigma_m median 0.69%; ~no cell-day below 0.3% — single frames are noise-dominated as
  expected (0.043% class), all signal claims must come from temporal aggregation (slow-curve stage /
  composites). |m|>5% extremes = 1.3% of cell-days, median SIGM 1.24 (poorly-illuminated, maskable).
  NEXT: I7 annual-comparison protocol (bar r~0.8), then slow-curve fitting + split-network coherence.
- USER DECISIONS 2026-08-22 (afternoon): SKIP I7 annual comparison (annual product inherits its own
  flaws — no era offsets, deseason, old-sign npz); go straight to B4. Merlin B4 audit done (see
  transmission_T.py docstring; key: c/BCD channel = 2nd shrinkage pathway, amplitude claims need
  end-to-end BCD injection; crosstalk transfer matrix required). transmission_T.py steps 1-3 RUNNING.
  B4 outputs = separate files in res_catalog (nothing overwritten).
- B4 STEPS 1-3 RESULTS (transmission_T.py, run 2026-08-22): checks PASS (null=0 exact, 2x-lin=0
  exact, mux-vs-single max 9% at weak-dia patches). Deep-band mean T: 0.49 (2009) -> plateau
  0.63-0.67 (2011-2026) — WELL above 0.1 kill floor, deep program viable; roster artifact = 2009-10
  ramp + ~6% drift. BUT crosstalk: 9/10 deep rows fail the 30% bar at ALL epochs — dominant leak =
  directly-overlying shallow patch (down-dip smearing, med 0.45-0.64x dia, max 2.7x at lat 40.4/
  44.8/46.0); secondary adjacent-deep 0.3-0.9x. Leak pattern epoch-STABLE; cond(Tqp)=40. Clean:
  lat 41.6/42.5. Sent to Merlin round 2: deconvolve vs mask vs coarsen vs per-segment 2x2 unmix.
- MERLIN ROUND-2 (B4 crosstalk): 30% bar MISCALIBRATED for the non-orthogonal Gaussian basis —
  adjacent-deep "crosstalk" ~= Gram overlap (identity-estimator baseline 0.188 vs observed 0.193);
  only shallow->deep down-dip leak is real estimator smearing (bad: lat 40.4/44.8/46.0). RULINGS:
  UNDAMPED 20x20 deconvolution (mask on propagated noise, not leak ratios); factorization C (epoch-
  mean shape) x daily diag T APPROVED w/ 3-epoch systematic; masking/coarsening/2x2 REJECTED; rebuild
  daily T single-patch (multiplex 9% fails bar; archived _multiplexed.npz); add 2009 epoch; aliasing
  probe; end-to-end BCD runs THROUGH deconv pipeline.
- B4 STEPS 1-5 v2 RESULTS (transmission_T.py v2 + transmission_deconv.py, 2026-08-22 evening):
  single-patch daily T: deep 0.487(2009)->plateau 0.62-0.66; checks exact 0 (4 epochs). 2009 shape
  deviates -> 2009-10 deep MASKED. Noise propagation (full cov, Neff=span/75=85): C^-1 amplification
  benign (med ~1, max 1.95) — masking driven by raw patch posterior. SURVIVORS at 0.05% agg bar:
  deep 47.2 (0.040%) + deep 48.2 (0.036%) ONLY; near-miss <2x: deep 41.6/42.5, shallow 41.4/47.2/
  48.4; leak-bad rows priced themselves out; ALL shallow masked (LFE kernels are deep — network is
  a deep instrument). Aliasing probe PASS 0.83 (bar 0.8-1.2). 18/20 masked. Sent Merlin round 3:
  2-survivor partial acceptance, per-claim-type Neff (slow-curve promotes near-misses?), mirror-v2
  before freezing mask, aggregate construction, nan-column handling. Products:
  res_catalog/transmission_deconv.npz/.png.
- MERLIN ROUND-3 RULINGS (B4): (A) hybrid accepted + merged-tile eval from full coeff covariance;
  (B) Neff promotion REJECTED (k-dof curve sigma = sig/sqrt(same Neff); promotions only via
  empirical circular-shift nulls or mirror v2); (C) mirror v2 BEFORE freezing mask, but needs USER
  sign-off (supervisor paused mirror cleaning 2026-07-22); mask = PROVISIONAL-v1; v2/v1 sig ratio
  <1.3x => near-miss tier dead, hybrid final; (D) aggregate = precision-weighted mean of ALL 10 deep
  deconvolved coeffs (frozen w, published); (E) 2009 stat repaired (nan-safe, enumerated, restricted
  headline). Probe: ~20% ONE-SIDED LOW bias systematic for unaligned structure.
- ROUND-3 AMENDMENTS RUN (transmission_deconv.py v2): merged tiles CONFIRM anticorrelation
  (-0.2..-0.6): **q11+q12 (lat 41.6/42.5, ~250 km) PROMOTED at 0.047%**; WA merge 0.022%. CLAIMABLE
  DEEP ROSTER (PROVISIONAL-v1): deep 47.2 + 48.2 individually + merged 41.6/42.5. Deep-margin
  precision-weighted aggregate sigma_agg 0.0168% (w: 70% WA, published). 2009 restricted dev 0.757 >
  spread 0.268 -> masked for SHAPE INSTABILITY (confirmed, not just illumination). All in
  transmission_deconv.npz (SDEC covariance, merged table, w_deep, agg series, systematics).
  AWAITING: user sign-off on mirror v2 (blocks mask freeze -> end-to-end BCD -> empirical nulls).
- USER DECISION 2026-08-23: the mirror-corrected ETS value (+0.043%) is RETRACTED as a signal-class
  anchor ("forget the mirror corrected value, its not good"). Consequences: the 0.05% masking bar in
  transmission_deconv is a PLACEHOLDER convention, NOT a validated threshold; tile sigma RANKING
  stands but NO claimable/unclaimable verdicts may be stated until the empirical per-claim nulls
  (Merlin r3-B) set data-derived bars. Mirror-as-control (certification) unaffected. Mirror-v2
  enthusiasm should be considered LOW pending user word. Memory: mirror-corrected-ets-retracted.md.
- USER DECISION 2026-08-22 (late): EXCLUDE 2009-2010 from the deep transmission-corrected product
  (ratifies the shape-instability mask; already enforced in transmission_deconv.npz A/agg pre-2011
  NaN). Mirror-v2 sign-off still PENDING.
- 2026-08-23: OLD-SIGN/SUPERSEDED PRODUCT PURGE (user-approved): deleted annual old-sign set
  (inversion_4d.npz + annual/ets/mirror/monthly figs, anom_corrected.parquet, common_mode movie/fig,
  step_scan*), raw daily movie, res_catalog_g20/g40/g50, June calT35 npz+null caches,
  inversion_4d_meta.csv. KEPT: G.npz, era_table, pairs/cells, parquets, daily products, July
  checkerboard/atlas figs (until redo lands). PENDING user OK: June fault_4d*.npz family at
  inversion/ level. FULL SIGN AUDIT of the daily chain: all hops consistent (see memory
  dvv-sign-convention).
- CHECKERBOARD REDO ACTIVATED (user 2026-08-23) — Merlin audit r1 rulings: noise-scaling surrogate
  for temporal averaging APPROVED (multi-frame rejected — asserts tau=1d, worse); mean response for
  W=365/1825 from 3/5 sampled days (roster drift); tau BRACKET {30,75,150} d (75 unvalidated; verdict
  flip across bracket => split-half measurement gates the axis); injection: per-band = HEADLINE,
  joint-opposite = adversarial bound, joint-same = leak-inflation demo ONLY (map-view checker gives
  bands same sign -> leak flatters deep; my joint proposal was backwards); bar 0.7 + July ladder
  {44,70,100,140,180,250}; R=8 seeded; caveats: corr amplitude-blind, 1% SNR legend, long-W excludes
  c-channel (B4 end-to-end); BRIDGE ASSERTION deep/W=365/tau75/2023 within one step of July 70-140
  else HALT; ATLAS = STALE-FOR-DAILY (TICKET: PSF rows under daily operator, deferred); July 3-D
  VOLUME verdict OUT OF SCOPE (stands). Script checkerboard_daily.py RUNNING.
- CHECKERBOARD REDO RESULTS (2026-08-23, checkerboard_daily.{npz,csv,png}): **BRIDGE PASS** (deep,
  W=1y, tau=75, 2023 -> crossing 100 km, inside July 70-140; July figs may be archived). Headline
  (per-band, tau=75, bar 0.7): single frames + 90-d resolve NOTHING at any scale/either band
  (movie = qualitative, as stated); W=1y deep resolves only 2019+ (180/100/250 km), shallow 140 km
  all epochs; W=5y deep 100 km everywhere (70 at 2023), shallow 70 km. Leak check: joint-same
  FLATTERS deep (70 km, corr .80) vs honest per-band (100 km, .74); joint-opposite never crosses
  (.50@100) — reality between per-band and opposite. Anchor: 1% deep 100-km checker = 0.0294%
  data-space. **CRITICAL: tau BRACKET FLIPS the W=1y verdict** (tau=30: 70-140 everywhere; tau=75:
  2019+ only; tau=150: never) -> per Merlin r1(b) the AVERAGING-AXIS ROWS ARE PROVISIONAL and the
  split-half empirical noise-decorrelation measurement (odd/even family difference series, variance
  vs averaging window) is now GATING for resolution-with-averaging claims. NEXT: split-half tau
  measurement. — the July
  resolution study certified a different estimator (lam 0.433 vs 0.215, no per-day site terms, no
  era offsets/new gates, epoch-stacked noise). Design: scale sweep 50-250 km through day_solve,
  iid+station-common noise, several roster epochs, evaluated per-frame AND vs temporal averaging
  (resolution as f(averaging window, year) = what the slow-curve stage needs). Merlin-check first.
