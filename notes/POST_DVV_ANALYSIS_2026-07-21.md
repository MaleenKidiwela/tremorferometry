# Post-dv/v Analysis Log — Resolution + 4-D Inversion (2026-07-21, LIVING DOC)

Everything AFTER dv/v completion: the deep-resolution characterization, the 4-D δβ/β inversion, and the full
adversarial-audit trail (every test, every Merlin round, every bug caught, every fix). Companion: the pipeline/
coverage methodology is in `METHODOLOGY_END_OF_RUN_2026-07-21.md` (§1–6). Scripts live in
`fault_tomography/inversion/`. **STATUS: resolution DONE (Merlin-signed); inversion IN the adversarial loop
(round 3b) — not yet signed off. Numbers below marked [prelim] until Merlin agrees.**

Working principle for this phase (user directive + `dont-over-conclude` memory): **loop with Merlin as an
adversarial referee — fix → re-run → re-audit — until he agrees the result is correct.** No dramatic conclusion
is stated before the setup is audited. Every surprising/negative result is treated as a probable bug first.

---

## 1. Deep-resolution characterization (SIGNED OFF)

**Goal:** how well is the DEEP plate interface resolved, and how does it improve as coverage grows 2009→2026?
NO velocity inversion here — only resolution/information characterization.

**Operator:** δβ/β on 0.10° interface cells at Slab2 depth; single-scattering early-coda ellipsoid kernel
(2–4 s window); per-station SITE terms; graph-Laplacian smoothing. Combined network = ALL 12,703 certified
families / 187 stations → 6,726 (cell,station) data → 872 cells (456 deep), 173 distinct sites.

**Captured-fraction correction (retiring the old 0.77).** The coda kernel is VOLUMETRIC; a thin interface
captures only a **median ~5%** (±1.5 km slab, grid-converged) of the 3-D shell mass. The historic checkerboard
used a unit-sum kernel that renormalized every ray to full interface sensitivity — inflating interface
amplitudes ~12×. Fix (Merlin, adopted as the PRODUCTION operator): multiply each row by its captured fraction f
(calibrated Gc = unit-sum × f). This is calibration + de-inflation in one.

**Gates A & B (are the family-specific residuals removable noise or a real floor?):**
- Gate A (mirror confound): on the 7 mirror-equipped stations, raw vs mirror-corrected post-site residual moves
  only ~5% (0.331%→0.316%) → the family-specific floor is REAL, not the removable mirror artifact.
- Gate B (statics): post-site residual is 99% TIME-VARYING, 1% static; lag-1 month autocorr 0.60.

**λ AUDIT — the bug the user caught.** My first checkerboards reused the atlas PSF-floor λ=2.069 (tuned on a
different, monthly-whitened system) → over-smoothing → artificially low recovery (I reported 0.06/0.69 and
prematurely concluded "index, not map"). Merlin's code audit found the λ transplant. CORRECTED with a per-scale
ORACLE λ:

| interface recovery (deep, oracle λ) | noise-free (geometry) | σ=0.13% (measurement) | σ=0.25% (residual) |
|---|---|---|---|
| 44 km  | ~1.00 | 0.66 | 0.43 |
| 250 km | ~1.00 | 0.84 | 0.73 |
| crossing scale (corr>0.7) 2011→2025 | 44 km always | 180→70 km | none→~140 km |

**Result:** the geometry resolves down to the 44 km cell scale noise-free (NO geometric obstruction); the limit
is SNR. The fleet sharpened the finest resolved scale ~180 km (2011) → ~70 km (2021+) at measurement noise. A
±1% interface checker produces ~0.035–0.05% data — matching the observed +0.043% ETS transient (calibration is
right, not pessimistic).

**INTERFACE vs VOLUME (which old number survives):** ablation on the SAME network — under OLD assumptions
(unit-sum, no site, inverse-crime noise) the new bigger network recovers 0.92/0.96, HIGHER than the old 0.77/
0.86 (the network improved). The honest drop is almost entirely CAPTURE-WEIGHTING (thin interface sees ~5% of
volumetric sensitivity). The **3-D VOLUME checkerboard reproduces ~0.86 honestly** (0.85 @ 250 km incl. depth) —
so **only the 2-D INTERFACE 0.77 is retired; the VOLUME finding stands.** Broad volumetric deep δβ/β is
resolvable at large scale; interface-confined change is the harder, more specific hypothesis.

**Atlas products (λ-independent, stand):** resolved deep modes ~3 (2010–2020) → ~24 (2021+); deep-index
precision 0.64%→0.12%.

Scripts: `build_G_captured.py`, `sensitivity_atlas.py`, `checkerboard_corrected.py`, `checkerboard_grid_corrected.py`,
`volume3d_checkerboard_new.py`, `checkerboard_ablation.py`. Figures in `res_catalog/`.

---

## 2. The 4-D δβ/β inversion — adversarial Merlin loop

Estimator (mirrors production `invert_epoch.py`): per window m = argmin ‖W(d−[Gc S]m)‖² + λ_f²‖Lm‖² +
λ_s²‖m_site‖²; +K sign convention; anomaly (deseason + demean); annual cadence; 0.10° grid; masking + null.

**Round 1 audit — NOT presentable, 12 issues** incl. FATAL: (1) SIGN bug — used −K, production is +K → every
map inverted; (2) NULL rigged — compared mean-of-|cells| (noise magnitude) to |mean-of-cells| (signal), so
"within null" was true by construction; (3) demean over each pair's own record mixes baselines → can fake the
roster-growth drift. Plus λ tuned on the wrong (pooled) system, closure one-sided, σ ignored autocorrelation,
ETS margin-wide/rate-nonstationary, VR from n=1 permutation, etc.

**Round 2 (all 12 fixed) — structurally sound; 6 residual defects; the "bound" not actually computed.** Results:
injection recovers a planted +1% patch at +0.89% (correct sign) → machinery validated; deep index within the
MATCHED null and in the long-record control → velocity-stable holds; but with autocorrelation-corrected σ (ρ=0.54)
and honest λ, ~0 cells pass the annual detectability bar. Residual defects: composite σ skipped the N_eff
correction; recent-year index biased to zero by short-pair fallback; VR/closure/null under-sampled; 0/456 is the
wrong summary (report per-element MDA); the injection was a TRANSMISSION test, not a detectability BOUND.

**Round 3 (6 fixes + the amplitude-sweep bound)** — running/complete; the bound = smallest coherent-patch
amplitude whose recovery exceeds its own year-scramble null envelope (per deep patch N/C/S).

**Round 3b (user-prompted seasonality fix + 2026 drop)** — see §3.5. RUNNING.

**[prelim] current inversion result:** deep interface velocity-STABLE within resolution — the spatial-mean index
stays within the honest matched null (incl. the long-record roster-growth control); a planted ~1% coherent deep
patch IS recovered, so this is a genuine BOUND not a coverage failure; no individually-resolved annual cell maps
under honest noise; ETS composite consistent with ZERO (expected for the RAW tensor — the ETS signal needs the
mirror-corrected v2). Everything labeled **v1/raw-tensor**; the fleet-wide mirror pass (v2) is the biggest lever
(it lowers the noise floor and, since geometry already resolves 44 km, buys resolution directly).

Scripts: `invert_dvv_4d.py`, `invert_4d_figures.py`.

---

## 3. Tests & investigations (the audit trail)

3.1 **Captured fraction** — median ~5% (±1.5 km), grid-converged 0.99; offset- and station-type-independent →
    intrinsic surface-vs-volume geometry, not a bug.
3.2 **Gate A / B** — family-specific floor is real (mirror removes ~5%), 99% time-varying (§1).
3.3 **Sign convention** — production `invert_epoch.py:16` flips forward −K to **+K** (dvv=∫K·δβ/β); G.npz stores
    −K; inversion must flip. (Round-1 FATAL bug.)
3.4 **Matched-statistic null** — scramble YEAR labels within each pair (≥50), compare the SAME statistic
    (|deep spatial mean|, and per-patch mean) real vs null; report the real value's percentile, not a binary.
3.5 **Seasonality (user-prompted, 2026-07-21) — key findings:**
    - Raw dv/v has a real ~0.12% seasonal cycle (winter +, late-summer/fall −).
    - 2-harmonic per-pair deseason removes ~70% (residual 0.038%); MORE harmonics don't collapse it (0.025% at
      4 harmonics) → the residual is **non-stationary** seasonal, not higher-harmonic shape.
    - Sampling season is STABLE across years (mean sample-month 6.4–6.9) → the residual does NOT drift → does
      NOT alias into a fake secular trend (the velocity-stable conclusion is not contaminated by it).
    - **2026 is a partial year** (sample-month 3.9, winter/spring only) → its index carries a different seasonal
      bias → **DROPPED from the index/null**.
    - **27% of pairs have NO significant seasonality** (F-test p≥0.05); unconditional deseason over-fit them
      (~5.7% variance spuriously removed). Seasonal amplitude is location-dependent (S 0.24% > C/N 0.18%; no
      depth dependence). → FIX: **conditional, location-dependent deseason** — per-pair harmonics subtracted
      ONLY where F-test-significant; short pairs use their LAT-BAND regional seasonal (not one global signal),
      also conditional. Implemented in `invert_dvv_4d.py`; re-running (round 3b).
3.6 **Autocorrelation** — monthly anomalies have ρ≈0.54 (30-day rolling stacks) → N_eff × (1−ρ)/(1+ρ); inflates
    σ_m ~1.7× and correctly deflates the resolved-cell count.
3.7 **Synthetic injection** — plant a known +1% deep patch into REAL data, invert, confirm sign+amplitude
    (recovered +0.89, correct sign) → validates operator + catches the sign bug (round-1).
3.8 **Instrument/metadata-change contamination (user-prompted, 2026-07-21) — OPEN, with Merlin:**
    - Instrument/metadata changes (sensor swap, response, gain) mid-record create SPURIOUS dv/v STEPS (not
      velocity changes). Across such a change the two eras are not the same station.
    - The dv/v uses an ALL-TIME reference (`dvv_roll30cal.py:53` ref=Rf.mean(0)) — NOT per-era — so the step is
      NOT corrected upstream; it enters the tensor.
    - Prevalent: fleet stations have within-band response-epoch changes (UW.WISH BHZ @2011-06, UW.LCCR @2013-06,
      most stations multiple Z epochs). Anchors PGC/SHB sensor swaps WERE handled as seams; fleet stations (one
      band via pick_band) were NOT checked for within-band changes.
    - Risk: a per-station-per-YEAR site term absorbs the common-mode annual-mean step, but the FAMILY-SPECIFIC
      part of a frequency-dependent response change leaks into the FAULT (interface) term as a spurious anomaly;
      a per-pair demean straddling a step uses one wrong baseline for a two-level series.
    - **Merlin ruling (round-4 blocker for the INDEX/stability claim only; resolution products unaffected):**
      - **GAIN changes are PROVABLY HARMLESS** — stacks are peak-normalized (`dvv_roll30cal.py` pk-division) and
        the stretch is correlation-based → pure gain/sensitivity changes have ZERO effect. Only PHASE-response
        (waveform-shape) changes bite (rarer; most FDSN epochs are gain/metadata).
      - Severity = "second-order of a potentially huge number": a sensor swap shifting group delay 5–10 ms in
        2–8 Hz → 0.17–0.33% apparent step (4–25× the signal budget). The per-station-per-YEAR site term absorbs
        the station-COMMON annual-mean step fully; only the FAMILY-SPECIFIC leak (coda-spectral-shape dependent)
        reaches the fault term — must be MEASURED, not assumed.
      - **The gate is a one-day MEASUREMENT (step-scan), not a rebuild:** `step_scan.py` — fit a step at each FDSN
        response boundary per pair; COMMON=median(b over families), LEAK=MAD(b); PLACEBO at random dates; +cc_max
        step + boundary-date histogram vs the index (clustered-upgrade discriminator). THRESHOLDS (data-space %):
        common median |b| > 0.05 → per-era ref mandatory (v2); LEAK > 0.02 → era-split rerun (v1) before present.
      - **Direction: both bad.** Incoherent steps widen the null → "within null" anti-conservatively easier.
        Clustered upgrade waves (2010–2015 UW/PB) → coherent spurious index drift that could FAKE the elevated-
        2010s pattern (or mask a real one). Since the index sits ~0.89 of the null, either can flip the verdict.
      - **Fix if triggered:** site tag = station×era (within year), per-(pair,era) demean, exclude ~35 days
        straddling a boundary; v2 = per-era reference in `dvv_roll30cal.py` (CHEAP — re-runs only dv/v from
        existing npz, no re-densify) folded with the fleet mirror pass.
      - Root cause: `scripts/plot_dvv_metadata.py` (metadata-overlay QC) was **never run fleet-wide** — precedent
        was CPW ROCK1 2011–2019 (a bad-sensor era of exactly this failure class).
    - **STEP-SCAN RESULT (`step_scan.py`, 156 boundaries / 76 stations):** COMMON step median 0.142% (placebo
      0.073%), family LEAK median 0.192% (placebo 0.167%) → placebo-corrected instrument leak ≈ 0.095% quadrature
      = **~5× the 0.02% bar**. My scan under-counts (all-band boundaries + gain-only revisions dilute toward
      placebo → the medians are LOWER bounds). Silver lining: USED-band physical boundaries are LOW in 2010–2019
      (2–9/yr) → the elevated-2010s index is NOT boundary-driven (an argument, not a clearance).
    - **MERLIN VERDICT: TRIGGER FIRED → ERA-SPLIT NOW. Inversion NOT presentable as "velocity-stable" until an
      AFTER-fix scan is clean.** Round-4 plan (in progress):
      1. `build_era_table.py` — per-station USED-band physical boundaries (poles/zeros change only; gain-only &
         unused-band excluded; merge <30 d). **RUNNING.**
      2. Refined before-scan: ≥5 placebos/station → per-station placebo distribution; attribute a boundary if
         common-|b| > its station's 95th placebo pctile; leak = quadrature excess of attributed boundaries only.
      3. Era-split inversion: site×era columns (incl. the ETS composite, `invert_dvv_4d.py:157`), per-(pair,era)
         demean (common-window WITHIN era, ≥6 mo else era-mean+flag), drop ±35-day boundary straddle, drop eras
         <6 mo. Split ONLY at used-band physical boundaries.
      4. After-scan (mechanical GATE): corrected-anomaly common |b| within ~1.3× placebo AND leak excess <0.02%.
         If fail → the era table is incomplete; iterate the TABLE, not the thresholds.
      5. Re-run index + matched null + percentile; before/after index curves side by side.
      - Presentable only after the after-scan passes, with a QUANTIFIED caveat ("instrument-era systematics
        bounded by X% common / Y% leak post-correction"). v2 upstream fix = per-era reference in
        `dvv_roll30cal.py` (cheap) folded with the fleet mirror pass.
      - If the used-band table shows <10 physical swaps AND attributed leak already <0.02% → era-split is a
        near-no-op and the blocker dissolves (do it anyway, ~free).

---

## 4. Honest conclusions (current, [prelim] until Merlin signs off)

- **Resolution:** the deep interface is a resolvable MAP at large scale (~70–140 km), sharpening with coverage;
  geometry resolves 44 km noise-free; the limit is SNR, not geometry. Volumetric deep structure resolves better
  (~0.85 @ 250 km) than interface-confined; fine scale (<~70 km) is noise-limited.
- **Inversion:** deep megathrust velocity-STABLE within resolution — a null-gated bound. No resolvable secular
  (annual) or ETS deep change at current (raw-tensor) sensitivity. The network WOULD detect a ~1%-class coherent
  deep anomaly (injection-validated); none is present. This reinforces (now spatially) the project's standing
  ETS-null / deep-velocity-stable result.
- **The lever:** fleet-wide MIRROR correction (v2) lowers the noise floor → directly sharpens both the resolution
  and the inversion detectability. Durable stacks survive on disk; it's a CPU pass.

---

## 5. Open items
- [ ] **Instrument/metadata-change handling (§3.8)** — Merlin ruling pending; likely per-response-epoch site
      terms + per-era demean, or a step-test to flag/split. BLOCKER for sign-off.
- [x] Conditional location-dependent deseason (§3.5) — implemented, round 3b running.
- [ ] Merlin sign-off on the inversion (round 3b → round 4 with the instrument fix if needed).
- [ ] Update this doc + figures with the FINAL signed-off numbers (headline detectability bound; index null
      percentile; MDA).
- [ ] Fleet-wide mirror pass (v2) → re-run inversion (the noise-floor lever).
- [ ] Two-depth shallow-leakage probe (quantify off-interface contamination).
- [ ] Drain re-run queue (6); NRCan retry for CN-sparse (3).

## 6. Audit-loop ledger (rounds)
| round | what | verdict |
|---|---|---|
| Resolution R1 | captured-fraction fork | proceed, capture-weight (production fix) |
| Resolution R2 | interface under-recovers? | λ-transplant BUG found (user caught); oracle-λ fix |
| Inversion R1 | first inversion | NOT presentable — 12 bugs (sign, rigged null, demean drift) |
| Inversion R2 | 12 fixes | structurally sound; 6 residual defects; bound not computed |
| Inversion R3 | 6 fixes + amplitude-sweep bound | running/complete |
| Inversion R3b | conditional deseason (§3.5) + drop 2026 | done (superseded by R4) |
| Inversion R4 | instrument/metadata era-split (§3.8) | DONE (pending Merlin sign-off): era-table 34 boundaries/32 stns; era-split AFTER-scan PASSES (attributed 6→0, common 0.14→0.02%, leak<0.02%); index exceedance PERSISTS @94th pctile after era-split + long-record control — BUT the apples-to-apples MIRROR test on the 7 mirror stations COLLAPSES it (coda−mirror decline −0.002%). +0.27% model-space = ~0.013% data-space = noise level. **Verdict: v1 exceedance = uncorrected NOISE-FIELD (mirror) drift, not deep signal; velocity-stable within resolution stands; v2 fleet-wide mirror pass = the definitive test.** |
| MIRROR check | 7-station coda vs own mirror | apples-to-apples: coda decline fully explained by mirror → collapses. Mismatched-network version looked encouraging (a trap — avoided by running apples-to-apples first). |
| Spatial-detector R6 | is the network-mean hiding a structured deep signal? | User instinct (correct): the network-mean deep index is orthogonal to zero-mean spatial patterns + dilutes asynchronous N/S ETS → ≥10× SNR lost (Merlin); "index within null" only excludes a margin-uniform change. Ran CONFIRMATORY per-cell LOCAL-ETS-locked composite (`ets_composite_confirm.py`, PRE-REGISTERED: primary=deep-cell mean, predicted sign NEGATIVE, detection⇔one-sided p<0.01 on ≥500 per-pair circular-shift nulls + gates G1 year-shift/G2 both-grids/G3 LOSO/G4 no-mixed/G5 n_days). RESULT both grids mirror-free: sign CORRECT (g50 −0.044%, g20 −0.027%), circular-shift **p=0.058 / 0.054 → FAILS p<0.01, misses p<0.05**; PASSES LOSO (sign holds at ALL ~180 stns, ≤1% atten), G4, G5. Merlin exploratory p≈0.03 (30 nulls, multi-stat) = FORKING TAX, did NOT survive 500-null confirmation. G1 year-shift gate BROKEN (all-NaN, shift-lookup bug — moot, primary fails; not counted as evidence). Grids share data → consistent NOT independent. Monthly 0.37% index DISPOSED: 82nd pctile under corrected max-over-months null (`plot_deep_index_monthly.py` null-construction bug fixed + hardcoded title removed). **VERDICT: sub-threshold (~1.6σ), sign-correct, spatially-robust LEAN — the strongest positive-leaning result in the project, but NOT a mirror-free detection (the ~60% outcome Merlin predicted). Mirror v2 (−65% noise) is the lever, blocked pending supervisor's mirror-method decision.** |
| Coarse-grid R5 | cell-depth / kernel source position | **FAMILY-CENTROID FIX** (user caught: 0.2° south had deep cells, 0.5° didn't). Root cause = cell depth AND kernel source were the *snapped geometric grid center*, not where the LFEs sit. Two bugs: (1) center-based deep(>30km) label is grid-dependent near the slab edge; (2) kernel horizontal source = grid center is up to **14.7 km median / 30 km max** off the family centroid on coarse grids → mispositions the 2–4 s coda shell by up to ~8.6 s (Merlin flag #1, the load-bearing one). FIX: represent each cell by its LFE **family centroid** (median lon/lat + family-median Slab2 depth), source each pair from its own family centroid; both grids. `max`-footprint (my first stab) REJECTED — balloons 32→70 deep, 3× area, kernel-inconsistent. VALIDATION: Slab2-at-family-position vs Lin 2023 located depths (n=1868, 48–49.3°N) → unbiased (+2 km) but ±13 km scatter; the >30 km label agrees only ~50% (chance) → deep boundary is **inherently fuzzy**, so add `depth_mixed` flag (frac_deep 0.3–0.7) and report conclusions with/without. Result: g50 deep 32→39 (+12 mixed), g20 deep 147→158 (+45 mixed); all deep cells illuminated (g20 all-zero rows 9.5% vs g50 44.6% → confirms the coarse-grid all-zero is thin-shell/coarse-cell discretization, NOT a bug — finer grid, fewer misses). Headline deep-index story UNCHANGED (still borderline: max|idx| at 100th pctile ~0.2% vs null ~0.11–0.13%, driven by 2010 endpoint) — this is a geometry/classification refinement, not a signal change. Files: assemble_res_catalog.py (src_lat/src_lon/src_dep, family-median depth_km, depth_center diagnostic), build_G_captured.py (CXY/CZ + per-pair source from family centroid). |
