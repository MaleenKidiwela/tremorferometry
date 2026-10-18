# S11 pilot (SVI) — separating deep from shallow dv/v with S-pinned residence kernels: rules fixed first

Written Mon 2026-10-12, ~1:12–1:17 pm PT, before any kernel, measurement or inversion number was computed.
Binding design: `MERLIN_DESIGN.md` (D1–D8). Decision log: `NOTES.md`. ε = −dv/v (positive ε = slower); every plot shows
dv/v = −ε. "Scrambled control" = the detector's c1 rows (each station site's data shifted to another time of the same
day) read at frame time (`c2` windows of S10 multi): noise through the identical path. Scope of this pilot: SVI only,
76 families × 7 stations (PB.B009/B010/B011/B926, CN.PGC/VGZ/SNB) on the existing windows
`svi_science/results/s10_dvv/multi/win/<NET.STA>.npz`. Step 7 (POLARIS/LZB/NLLB + WA) is NOT run until review.

## 0. Status (updated as stages finish)
- 2026-10-12 1:17 pm PT — PLAN written. Nothing computed yet.
- 1:20 pm PT — step 2 gate (a) part 1 done (tables/gate_a.csv; NOTES). Sensor depths: tables/sensor_depths.csv.
- 1:33 pm PT — step 3 kernels PASS the unit test (max |Σa−1| 2e-14; a_d < 0 for 0 % of pairs). KERNEL_DONE.
- 1:32 pm PT — queue running steps 4–6: `svi_science/scripts/s11_queue.sh` (PID in logs/queue.pid). Check: `cat logs/queue.log`;
  stage markers MEASURE_DONE, SVD_DONE, INVERT_DONE, INJECT_DONE, GATES_DONE, FIGURES_DONE, QUEUE_DONE in this folder.
  Results: tables/headline.json, inv_summary.csv, inject_summary.csv, gates.csv, gate_a.csv; figures/.
- 1:37 pm PT — QUEUE_DONE (steps 4–6 complete; re-run once after an SVD-stage KeyError fix). Headline path = raw.
  Post-hoc column cc0 (no cc cut) added 1:38 pm PT. Summary: RESULTS.md. Step 7 NOT started (awaits review).
- ~2:10 pm PT — Revision 1 (§3) run done (s11_v2.py, V2_DONE); RESULTS.md top section. No scale-up.
- 2:20 pm PT — Revision 2 (§4) done: queue s11_queue_q3.sh (ext inventory/windows/measure/kernels + s11_q3.py), QUEUE_Q3_DONE. Results in RESULTS.md top section.
- 2:50 pm PT — Step 7 WA-NW (§5) running: `S11_REGION=WA_NW` queue `svi_science/scripts/s11_wa_queue.sh` (log `wa_nw/logs/queue.log`,
  markers in `wa_nw/`, QUEUE_DONE at the end); WA-SW column chained after it (`s11_wa_sw_chain.sh`, `wa_sw/logs/queue.log`).
- 3:52 pm PT — Step 7 WA-NW done (QUEUE_DONE, post-hoc KERR_DONE); RESULTS.md 'Step 7' section. WA-SW closed: 0 measurable series.
- 10:41 pm PT — ROUND 5 (§9 + Amendment 1) queue running: `svi_science/scripts/s11_r5_queue.sh` (PID round5/queue.pid, log round5/queue.log, markers round5/*_DONE, QUEUE_DONE at the end).
- 2026-10-13 9:42 am PT — ROUND 7 (§11 + Amendment 1) QUEUE_DONE; round7/RESULTS.md.

## 1. Decisions D1–D8 as adopted (with executor choices made before any number)
**D1 Kernel** — new script `svi_science/scripts/s11_kernel_mc.py` (adapted; `fault_tomography/kernels/residence_kernel_mc.py`
is not edited). Monte-Carlo radiative transfer of S energy, isotropic scattering, isotropic source, perfectly reflecting
free surface, next-event estimation to the receiver at its **sensor depth** plus its **free-surface image** (equal
weights), straight legs.
- Velocity: regional USGS CVM v1.7 S profile `data/velocity_models/usgs_cascadia_v1.7/profiles/south_VI.csv`. Every
  straight segment is cut **exactly** at the union of the CVM layer edges and the reporting layer edges; time residence
  = length / β(layer); total path times (lapse) are the sums of these exact times.
- Reporting layers (km): 0–0.1, 0.1–0.3, 0.3–1, 1–3, 3–5, 5–10, 10–25, 25–35, 35–45, > 45. Zones (exclusive, priority
  in this order): near-source sphere r ≤ 10 km around the hypocentre; near-station vertical cylinder r ≤ 10 km; rest.
  Zone intervals along each segment are also exact (line–sphere and line–cylinder quadratics). 30 categories.
- Lapse bins 0.5 s past the direct S, edges 0 … 6.5 s (13 bins). Direct S = straight ray, time and residence exact;
  for a borehole receiver the anchor is the mean of the direct ray and its surface ghost (equal weights).
- Differential kernel per bin b: R_b = contribution-weighted mean time-residence vector (s) of coda paths in bin b,
  r_S = residence vector (s) of the direct S, τ_b = mean lapse past S in bin b. A stretch fitted about S over a window
  weights lapse τ by energy × τ², so the window coefficient vector is
  **a = Σ_b W_b τ_b (R_b − r_S) / Σ_b W_b τ_b²** (Σ over categories of a = 1 exactly for a homogeneous change;
  dv/v_apparent = Σ_k a_k δ_k with δ = δβ/β). (This replaces the fleet script's t̄-only form; logged in NOTES.)
- Lapse weights W_b from **each series' own grand stack** (`multi/tables/grandstacks.npz` G1, envelope² per 0.5-s bin
  past its tS, Hilbert envelope; for 3C each component is normalised by its own window energy, as the joint z-scored
  3C stretch does). Z-only and H-only columns use their own component weights.
- Windows past S: **2–5 s (headline)**; 3–6, 2–3.5, 3.5–6 s as columns. The fleet "w24" (1–3 s past S) is not used.
- Variants: ell = 20 / 40 (primary) / 80 km; two-layer ell (5 km for z < 3 km, 40 km below); second seed at ell 40.
  8,000 particles per pair and variant.
- Groupings used by the inversion: a_s = shallow 0–5 km (all zones), a_m = mid 5–25 km, a_d = deep 25–45 km (all
  zones), a_b = below 45 km. Columns: near-station-cylinder share of a_s, near-source share of a_d, every category.
- **Unit test (pass rule):** per pair and variant, |Σ_k a_k − 1| ≤ 0.02. **Stop rule:** any failure, or a_d < 0 for
  > 10 % of pairs (primary variant, 2–5 s) → stop, report, re-define the deep layer before continuing.
- Reported: ell-spread of a_s (max/min over the ell variants, per pair; expected ~×2) and family-depth spread of a_d
  at each site (max/min over families; expected ≥ ×2).

**D2 Deep unknown** = volumetric 25–45 km layer, SVI region, one value D(e) per episode. Not the ±1.5 km slab.

**D3 Data unit (raw path, Merlin)** — per (family, station, response epoch, episode): P-weighted (P ≥ 0.5) stack of all
repeats, per-event 3C realignment (`s10_measure.realign`, ±0.5 s, on tS−1…tS+2 s against the epoch grand stack),
leave-episode-out reference within the epoch, **joint 3C S-anchored stretch** (each component z-scored on the window,
cc = mean of the component cc's = cc of the concatenated z-scored windows, as `s10_fam_3c_26.stretch`), ±1 % grid
0.01 %, parabolic peak, rail |ε| ≥ 0.99 % = missing. Columns: 2–5 (headline), 3–6, 2–3.5, 3.5–6 s; Z-only and H-only
(N+E joint) at 2–5 s; 2-parameter fit (shift ±50 ms × stretch) at 2–5 s; cc; 200-draw event bootstrap σ for every
3C window and for Z/H 2–5 s; pre-arrival (pre-P) RMS of the episode stack relative to the series median; c2 scrambled-
control episode stacks (equal weights, ≤ n data windows of the episode, same days) through the identical path.
- **Headline cells: n ≥ 50 repeats and cc3c(2–5) ≥ 0.8**; n 20–49 as a column set.

**Data path (ii), the user's recipe (SVD path)** — same rows, realignment and grand stack; P-weighted 3C **day**
stacks (days with ΣP ≥ 3), SVD-Wiener filter per component per epoch (`s10_svd_measure.svdw`, noise columns τ < tS − 8 s),
joint 3C stretch 2–5 s anchored at the picked S against the realigned epoch grand stack (as `s10_fam_3c_26`), days with
CC > 0.8 and |ε| < 0.99 %, **CC-weighted mean per episode**. Same window columns. Episode σ = CC-weighted day scatter /
√n_days (n_days ≥ 2), else the series' median σ. Cells: same episode, n ≥ 50 repeats, ≥ 1 day with CC > 0.8.
- The two paths share the cell definition (family, station, epoch, episode), the kernels, the inversion code, the
  injection code and the gate code. Only the measured value and σ differ.

**Headline rule (fixed before numbers):** statistic = median over series of the waveform-level step recovery of the
SVD path, 3C 2–5 s, at injected ε steps of 0.1 % and 0.3 % (each amplitude separately; step at an episode boundary
near the middle of each series; recovery = injected minus natural contrast / amplitude). **SVD path is the headline if
both medians lie in 0.9–1.1; otherwise the raw path is the headline** (< 0.8 is Merlin's rule; 0.8–0.9 and > 1.1 are
undefined there and are assigned to the raw path because it is the calibrated one) and the SVD path is reported beside
it with its recovery. The raw path's own recovery is measured with the same code and reported.

**D4 Inversion** — per path, window 2–5 s headline. Data d = dv/v (%) per headline cell. Model
d = o_series + a_s(series) S_site(e) + a_d(series) D(e). Unknowns: offset per (family, station, epoch) series (free,
absorbs the epoch zero and any instrument step; series with < 3 cells dropped); S_site(e) per site per episode
(sites: PGC = {PGC, B009, B010, B011}, B926, VGZ, SNB); D(e) per episode. Mid-crust and > 45 km not free; reported:
the factor by which a uniform 5–25 km change projects onto D (invert synthetic data a_m × 1 through the same operator).
Variant: S_site = two annual harmonics + per-year step instead of per-episode values. Weighted least squares, weights
1/σ² with σ inflated per series by √max(1, χ²_red of its scrambled-control cells about zero) (fewer than 3 control
cells → the station median inflation); Gaussian prior σ = 1 % on S and D (conditioning only); offsets eliminated
exactly (Schur complement). Reported: formal covariance, resolution-matrix diagonal, eigenvectors of the data normal
matrix (S, D block after offsets) whose data-only σ along the direction exceeds 0.1 % ("below the noise level"), and
their overlap with D ≈ −mean(S). **Prior test:** prior × 2 and × 0.5 must change D by < 10 % (RMS of the change /
RMS of D); else D is reported as unresolved. No λ scans, no oracle, nothing transplanted.

**D5 Time bins** — SVI episode units of S10 multi (`multi/tables/active_days.csv`, majors + between-episode bursts);
D per episode; S per site per episode.

**D6 Inclusion** — pilot: the 76 S10-multi families (Bostock hypocentres via `data/catalogue/families_v1_3.csv`,
identical to `lfe_catalog_repair/results/e5r_families_2012.csv`; depths 29–46 km, all inside 25–45 km except none
need replacement — checked at kernel time and flagged if not). Stations: the 7 with windows; per series causality
(coda/pre-P RMS > 1.5) and realign cc median ≥ 0.5 are applied as inclusion; failures kept as columns.

**D7 Gates (pass rules)**
- (a) Co-located, kernel-free (first): annual cycle fitted (A cos + B sin + const per series-demeaned episode values,
  WLS) per instrument PGC, B009, B010, B011 and per component (Z, H, 3C) on the raw path 2–5 s; first with the
  existing Z 2–4 s episode stacks of S10 multi, then the new columns. Report amplitude ± σ and phase per instrument and
  component, the χ² test of equal amplitudes across the 4 instruments, and the sensor depths (FDSN, table
  `tables/sensor_depths.csv`). No pass/fail; any outcome is informative.
- (b) Injection suite, differential, waveform level (day stacks stretched by the series' kernel coefficient × the
  injected amplitude about tS before any stacking/SVD; whole chain re-run; result = injected minus natural):
  - A: deep step +0.3 % from 2015-01-01 → D recovery 0.7–1.3, leak into every S_site < 30 %.
  - B: common-mode shallow annual 0.1 % (peak day 233) at all sites → leak into D < 25 % (D's annual amplitude /
    0.1 %). If > 50 % the headline quantity becomes "deep + regional shallow".
  - C: one-site shallow step +0.3 % from 2015-01-01, run once per site (PGC site, B926, VGZ, SNB) → D step / 0.3 %
    < 20 % for each.
  - D: kernel error: B's signal generated with a_s × 1.5 or ÷ 1.5 (random per pair, seed fixed), inverted with
    nominal kernels → leak into D < 30 %.
- (c) Scrambled control: D from the control cells (same headline cells, control value where not railed; no cc cut,
  noise cannot pass it) within 2σ of zero in ≥ 90 % of measurable episodes (≥ 5 control cells).
- (d) Split network: {PGC site} vs {B926, VGZ, SNB}, and families odd vs even: Theil–Sen slope of D_B vs D_A over
  common episodes in 0.7–1.3 and Pearson r > 0.5 (each split).
- (e) Station-specific systematic: leave-one-site-out D vs that site's own D: same sign in ≥ 70 % of common episodes
  (mean over sites); and across headline cells (series-demeaned) |Spearman ρ| < 0.2 between ε and the pre-arrival RMS
  change, and between ε and (1 − cc).
- (f) Lapse: D from 3.5–6 s vs D from 2–3.5 s (each inverted with its own kernels; the kernel-predicted ratio is
  therefore 1.0): Theil–Sen slope in 0.7–1.3. The raw ratio a_d(3.5–6)/a_d(2–3.5) is reported.
- (g) Seasonal: annual amplitude of D < 2σ while the S_site amplitudes carry 0.07–0.08 %; otherwise adjudicate with (a).
- Step-4 check: raw 3C 2–5 s per-cell bootstrap σ (median, cells with n ≥ 100) at the boreholes > 0.25 % → the study
  is a bound only; say so and continue.
- Failure of (d) or (e) → the D series is a station systematic; report the bound from the scrambled control and
  injection A only.

**D8 Staging** — kernels (8 workers) → measurement (12 workers, one station's windows in memory at a time) → SVD path
+ SVD recovery → inversion (both paths) → injection suite (both paths) → gates. One nohup queue
`svi_science/scripts/s11_queue.sh`, DONE markers in this folder. ≤ 16 CPU workers, no GPU, RAM < 100 GB.

## 2. Files
- Scripts: `svi_science/scripts/s11_{common,kernel_mc,measure,svdpath,invert,inject,gates,figures}.py`.
- Outputs: `tables/`, `figures/`, `logs/`; DONE markers `KERNEL_DONE`, `MEASURE_DONE`, `SVD_DONE`, `INVERT_DONE`,
  `INJECT_DONE`, `GATES_DONE`.
- Plot rules (user): dv/v points coloured by CC, no error bars on dv/v plots (σ in tables/text), CC-weighted averages,
  line graphs welcome.

## 3. REVISION 1 — pre-registered before any v2 number (written 2026-10-12, ~1:58 pm PT)
Source: `MERLIN_REVIEW_1.md`, NOTES decisions 9–12. No new measurement; same `tables/episode_raw.parquet`, kernels and
store. Raw path only (headline per §1 rule). Script `svi_science/scripts/s11_v2.py`; outputs `tables/v2_*`,
`figures/fig_v2_*`; marker `V2_DONE`.

**Cell rule v2.** n ≥ 50 AND cc3c(2–5 s) ≥ 0.35 (99th percentile of scrambled-control cc) AND series causality > 1.5,
finite ε and bootstrap σ; weights 1/σ² with σ inflated per series by the scrambled-control χ² (as §1 D4); series with
≥ 3 cells. The cc ≥ 0.8 rule is kept as a column.

**Parameterisation** (offsets per (family, station, epoch) series always free; weak prior σ = 1 (% or %/yr) on every
coefficient for conditioning only; t in decimal years; slow basis = linear (t − 2015) + hinges max(0, t − k),
k = 2007, 2011, 2015, 2019, 2023):
- `episode`: S per site per episode + D per episode (as §1) → per-episode D reported as a bound (2σ) only.
- `smooth` (HEADLINE): S per site = 2 annual harmonics (cos/sin 1 yr, cos/sin 0.5 yr) + slow basis (10 dof per site);
  D = slow basis + 1 annual harmonic (8 dof).
- `slow+annual`: as smooth but S per site with 1 annual harmonic only (8 dof per site) — the minimal symmetric model
  (our reading of Merlin's label; logged).
- `joint-lapse`: smooth model fitted jointly to the 2–3.5 s and 3.5–6 s cells, each with its own kernels and its own
  series offsets (column).
- Kernel variants as columns (ell20, ell80, two5_40); cc ≥ 0.8 cells as a column.

**Headline quantities.** Q1 = deep annual amplitude √(c² + s²) of D with σ (delta method on the posterior covariance,
scaled by max(1, χ²_red)) and its peak day. Q2 = decade contrast of the slow D curve: mean over 2015-01-01…2025-12-31
minus mean over 2004-01-01…2014-12-31 (monthly grid), σ from gᵀCg.

**Injections (waveform level, raw path, v2 cells, smooth model, injected − natural):**
- A deep step +0.3 % at 2015-01-01 (true Q2 = +0.3 %): Q2 recovery 0.7–1.3; |ΔQ1| < 0.03 %.
- B shallow annual 0.1 % (peak day 233) at all sites: Q1 leak |ΔQ1| / 0.1 % < 0.25.
- Coherent kernel error: A and B generated with the two5_40 and with the ell80 coefficients (a_s, a_d of that variant
  for every pair), inverted with ell40: same pass rules (A Q2 recovery 0.7–1.3, B leak < 0.25).

**Q1 gates (pass rules):**
- G1 leave-one-site-out: for every site, |Q1₋ₛ − Q1| < 2σ(Q1₋ₛ) and, where Q1₋ₛ/σ ≥ 1, peak within ±60 d of the full peak.
- G2 family split (odd vs even): |Q1_odd − Q1_even| < 2√(σ²_odd + σ²_even); peaks within ±60 d where both ≥ 1σ.
- G3 scrambled control: smooth model on the control values of the v2 cells (non-railed, own bootstrap σ, same
  inflation): Q1_ctrl < 2σ. Column (not "scrambled control"): episode-date permutation of the data cells, 500 draws,
  p of Q1.
- G4 station-amplitude pattern: per station, the kernel-free in-phase annual amplitude (cells series-demeaned, cos at
  the pooled peak, WLS); χ² against A ∝ mean a_d of the station's cells and against A ∝ mean a_s (one free scale each).
  "Deep pattern supported" iff p(∝a_d) > 0.05 AND χ²(∝a_s) − χ²(∝a_d) > 4.
- G5 lapse: kernel-free pooled in-phase annual amplitude ratio 3.5–6 / 2–3.5 s; predicted 0.83 (deep, ratio of a_d)
  vs 1.14 (shallow, ratio of a_s) from the kernels at these cells; reported as z to each; "deep supported" iff the
  observation is ≥ 1σ closer to the deep prediction than to the shallow one.
- SNB tension investigated (per-epoch amplitude, phase, Q1 without SNB) — diagnostic, no pass rule.

**Failure signatures / reporting.** Q1 is called a deep annual change only if Q1/σ ≥ 3 AND G1–G3 pass AND B and
coherent-B leak < 0.25; G4/G5 are supporting evidence. Otherwise Q1 is reported as a measured value with its 2σ upper
bound and "not attributable". Q2 is a detection only if |Q2|/σ ≥ 3 and A (incl. coherent) recovery is in 0.7–1.3;
otherwise a bound. Prior ×2/×0.5 must move Q1 and Q2 by < 10 % of their σ-scaled values (|Δ| < 0.1σ), else flagged.

## 4. REVISION 2 — Q3 = ΔD(e), pre-registered before any Q3 number (written 2026-10-12, ~2:15 pm PT)
Source: `MERLIN_FORECAST.md` ("Prescription" 1–3), NOTES decision 14. Raw path. Script `svi_science/scripts/s11_q3.py`
(+ `s11_ext_inventory.py`, `s11_ext_windows.py` for the 39 added families); outputs `tables/q3_*`, `figures/fig_q3_*`.
Correction carried from the forecast: per-episode deep/shallow ATTRIBUTION (episode model) is unidentified (series
offsets absorb a_d); the earlier per-episode σ_D was dominated by the unidentifiable common level of D.

**Model `temporal`.** d = o_series + a_s S_site(t) + a_d D(e); S per site = 2 annual harmonics + linear + hinges at
2007/2011/2015/2019/2023 (as v2 smooth); D free per episode (majors and bursts); offsets free; prior σ = 1 % on all
coefficients (conditioning); σ scaled by max(1, χ²_red); weights 1/σ² with control inflation (B010 overflow fixed:
control cells with σ < 0.01 % are dropped from the χ² and the station fallback is the median of its series χ²).
**Q3(e) = ΔD(e) = D(e) − mean D(e′) over the other episodes with |mid(e′) − mid(e)| ≤ 3 yr.** Composite = equal-weight
mean of ΔD over the 15 fixed majors: 2003-02, 2004-07, 2005-09, 2008-05, 2010-08b1, 2011-08, 2012-08, 2013-09, 2015-12,
2018-06, 2019-08, 2019-11, 2020-11, 2023-04, 2024-08 (big majors with ≥ 20 v2 cells; fixed now, not re-chosen).

**Configurations.** (a) 76 families, v2 cells (n ≥ 50, cc ≥ 0.35) — reproduction of the forecast; (b) 76 families,
n ≥ 20, B010 fixed; (c) HEADLINE: 115 families (76 + 39 added with the identical inventory/window/measure/kernel
chain; their dates map to the existing episode units, days outside every unit form new bursts by the same 7-day rule),
n ≥ 20, cc ≥ 0.35, B010 fixed.

**Acceptance:** σ_ΔD ≤ 0.10 % for ≥ 5 of the 15 majors in (c).

**Gates on (c) (pass rules):**
- prior ×2 and ×0.5: every major's ΔD and the composite move by < 0.1σ.
- scrambled control: same model on the control values of the same cells (non-railed, own σ, same inflation):
  ΔD_ctrl within 2σ of 0 for ≥ 90 % of the majors with ≥ 5 control cells; composite within 2σ.
- one-episode deep pulse: +0.3 % × a_d injected into the day stacks of ONE major (waveform level, full raw chain),
  once per major: recovery Δ(ΔD)/0.3 median over the 15 in 0.8–1.1.
- common shallow pulse: 0.1 % × a_s at every site in ONE major, once per major: leak |Δ(ΔD)|/0.1 median ≤ 0.35; the
  all-15-majors-pulsed composite leak is reported (forecast 0.41).
- leave-one-site-out: for each site, ΔD of the majors within 2σ (σ of the leave-out estimate) for ≥ 90 % of the
  majors measurable without that site; leave-PGC-out reported separately.
- Failure of the control or the leave-PGC-out gate → ΔD is a station systematic (reported as such).
- Wording: a major's ΔD is called a deep-weighted episode change only if |ΔD| ≥ 3σ and all gates pass; otherwise
  the values are reported with their 2σ bounds. Declared shallow leak: the measured shallow-pulse leak.

## 5. STEP 7 — WA-NW, the validated pipeline on an independent region (pre-registered before any WA number)
Written 2026-10-12, ~2:30–2:50 pm PT, before any WA window, kernel, dv/v or inversion number. Source: MERLIN_DESIGN D6/D8,
NOTES decision 15. The SVI pipeline runs unchanged; a region switch (`S11_REGION=WA_NW`, `s11_common.py`) sends every
WA input/output to `svi_science/results/s11_depth/wa_nw/{tables,logs,figures}` (markers there too). No SVI file is
written. Final figures are copied to `figures/fig_wa_*.png`. Raw path only (SVD path retired at SVI, recovery 0.26).

**Families.** `data/catalogue/families_v1_3.csv`, region Royer_NW_WA, location_quality A or B (113), that exist in the WA
production family table `lfe_extension_wa/results/wa_prod/setup/families.csv` (expected 108; the rest have no repeats).
Location = v1_3 lat/lon/depth. Depth outside 25–45 km → Slab2 depth at the epicentre
(`data/slab_models/slab2/cas_slab2_dep_02.24.18.xyz`, bilinear), flagged `depth_replaced`. C-quality families are not
measured in this run (no kernel; could not enter any quantity); their count is reported.

**Repeats (headline).** `selected_v2 == True` rows of `lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet`,
region NW (rule PNET ≥ 0.75 on every day class), 2007-01-01 → 2026-09-15. This already contains the E1 2007–2011 rows
(roles fit / prod / sealed / excluded-month); all are kept. Weight P = `P_own_pnet` (the WA confidence). Active family-
day = ΣP ≥ 3 (identical to SVI); only rows on active family-days are used (rows lost are reported).
**Scrambled control.** c1 rows of the production detector (`wa_prod/mf/prod/det_<date>.parquet`, net == core, series == c1),
same family, active family-days; at most n_data rows per family-day kept (seeded random; the measurement caps the control
at n data windows per day anyway); read at frame time t + S offset (the SVI "c2" convention).

**Stations, S offsets, epochs.** S offset of (family, station) = `s_off` of the WA template with the most events (template
library `wa_prod/setup/tpl_library.csv` + `templates/tpl_qc_w4b.csv` use == True; n_events ≥ 30, finite s_off). Spread of
s_off across a station's template nodes is a column (> 2 s flagged). S−P for the pre-arrival window from the royer_NW_WA
CVM (straight ray). Station set: every station with a template S for the family and hypocentral distance ≤ 80 km.
Response epochs: boundaries = union of the node start dates (`data/extension_stations/{NW_WA,SW_WA,OLY2010}_bridge_nodes_strict.csv`,
`step4b/setup/extra_nodes.csv`) and of `wa_pre/w2/epoch_changes.csv` (sensor, channel, gain and node changes); series =
(family, station, epoch), never referenced across epochs. Station-days in `wa_pre/w2/r2_bad_for_dvv.csv` with reason 'dead' are
not read. [Amended 2:36 pm PT, after the inventory count and before any window: 'flat' (≥ 1 % of samples repeating the
previous one) marks 6,760 GNW and 7,035 HDW days, i.e. their old quantised digitisers, not dead data; the window
reader already rejects zero runs. 'flat' is therefore not used.] Sites = stations < 3 km apart (single linkage; wa_common SITE_KM). Series gates as SVI: causality (Z 2–5 s coda /
pre-P RMS) > 1.5 is the inclusion rule; realign-cc median ≥ 0.5 is a column (SVI decision 1:25 pm PT: that cc is a single-
event SNR and removed 86 % of SVI series).

**Windows.** `s10_windows.task` unchanged (2–8 Hz Butterworth-4 zero-phase, 40 Hz, despike 8 × MAD, RMS-normalised, 3C,
−14..+10 s about S, read span S−15..S+30 s) with the Raw reader subclassed to read the WA raw roots
(`wa_prod_common.prod_roots()` = `wa_common.TRACE_ROOTS` + traces_dl2 + the production extra roots; first non-empty file
wins, as in production). Data and c2 passes stay separate tasks (same despike reference as SVI) sharing a per-worker day
cache. Window sizes recorded.

**Kernels.** `s11_kernel_mc.py` unchanged physics, CVM v1.7 `royer_NW_WA.csv`, sensor depths from EarthScope FDSN
(channel level; `wa_nw/tables/sensor_depths.csv`), lapse weights from the S11 store grand stacks, 100,000 particles, the 5
variants, only pairs with a measured series. Unit test |Σa − 1| ≤ 0.02; stop if a_d < 0 for > 10 % of pairs.

**Episodes.** Majors = `lfe_extension_wa/results/episode_rates/episodes_nw_final.csv` with major == True (69), unit window
start − 2 d … end + 2 d; other active days form bursts by the 7-day rule (as SVI).

**Cells.** v2 rule: cc3c(2–5 s) ≥ 0.35, causality > 1.5, finite ε and σ, series ≥ 3 cells, 1/σ² weights with the control
inflation incl. the 2:13 pm fix. **n ≥ 20 is the WA headline for every quantity** (as instructed for step 7; at SVI it was
the Q3 headline); n ≥ 50 (SVI's v2 rule for Q1/Q2) and cc ≥ 0.8 and the realign rule are columns.

**Models and quantities** (as §3/§4): episode (per-episode bound only), smooth (Q1 deep annual, Q2), slow+annual, joint-
lapse, temporal (Q3 = ΔD(e) vs episodes within ±3 yr). Same basis and knots. **Q2 for WA = mean 2015-01-01…2025-12-31 minus
mean 2007-01-01…2014-12-31** (WA repeats start in 2007; SVI's 2004 start would be extrapolation).
**Q3 major list (rule fixed now, list written by counting cells before any ΔD):** every major == True episode with ≥ 20
headline cells; composite = equal-weight mean of their ΔD. Column: "co-timed" composite = WA majors whose unit window lies
within ±30 d of one of the SVI 15 majors.

**Gates (pass rules identical to §3 and §4):** prior ×2/×0.5 (< 0.1σ); scrambled control (Q1 < 2σ; Q3 ≥ 90 % of measurable
majors within 2σ and composite within 2σ); deep pulse (median recovery 0.8–1.1) and common shallow pulse (median leak
≤ 0.35) per major; injections A (Q2 recovery 0.7–1.3, Q1 leak < 0.03 %) and B (Q1 leak < 0.25), each also generated with
the two5_40 and ell80 kernels (coherent kernel error); kernel-variant columns; leave-one-site-out (Q1 G1; Q3 ≥ 90 % within
2σ); family split (Q1 G2; plus Q3: odd/even composites within 2√(σ₁² + σ₂²) and ≥ 90 % of common majors within that
bound); station-amplitude pattern G4; lapse G5. Wording rules as §3/§4 (Q1 "deep annual change" only if Q1/σ ≥ 3 and
G1–G3 and B pass; Q3 episode "deep-weighted change" only if |ΔD| ≥ 3σ and all Q3 gates pass).
**Single-site dominance (W060 lesson):** per site, share of the total cell weight Σ1/σ²; a site above 50 % is flagged and
its leave-out result is the guard; per-station kernel-free annual amplitude more than 3 robust σ from the station median
(MAD-scaled) is flagged as a possible bursty / instrument station (column, no removal).

**Comparison with SVI (no pass rule).** Q2 and Q3 composites side by side (SVI headline values from RESULTS; SVI Q1/Q2 also
recomputed read-only at n ≥ 20 for a like-for-like column); calendar-year inverse-variance means of ΔD per region with σ;
annual cycle "matches" if the Q1 peak days are within ±60 d and the amplitudes within 2σ combined (also for the kernel-free
pooled in-phase amplitude).

**WA-SW (column, only if cheap).** `S11_REGION=WA_SW`, Royer_SW_WA A/B families, catalogue rows region SW with
P_own_pnet ≥ 0.75 and day_class == tremor, NW episode units, the same chain. Reported as counts; Q1/Q3 only if ≥ 100 cells;
no gates.

**Staging.** `svi_science/scripts/s11_wa_queue.sh` (nohup): inventory → windows → measure → kernels → v2 → major list → q3
→ WA extras/comparison/figures. Markers in `wa_nw/`. ≤ 14 CPU workers, no GPU, RAM < 100 GB.

## 6. ROUND 3 — kernel-robust deep term, REDPy estimators in Q3, coarse along-strike zones (pre-registered before any round-3 number)
Written Mon 2026-10-12, ~7:25 pm PT, before any round-3 kernel, dv/v, injection or forecast number. Both regions (SVI, WA-NW).
Source: user via main session ("keep trying until something works"; coarse resolution fine if new), NOTES decision 20,
`redpy/RESULTS.md`. Outputs `round3/{svi,wa_nw}/{tables,logs}`, `round3/figures`; scripts `svi_science/scripts/s11_r3_*.py`;
queue `s11_r3_queue.sh` (markers `round3/<STAGE>_DONE`). No earlier result file is changed. Design by the executor (no
Merlin call: Fable-quota rule); every choice below is fixed now.

**Common.** Raw path. Cells as the headline of each region: WA-NW = §5 (n ≥ 20, cc ≥ 0.35, causality > 1.5, control
inflation, series ≥ 3; 2,466 cells); SVI = §4 configuration (c), 115 families, n ≥ 20 (3,531 cells) — like-for-like with WA.
SVI v2 (76 fam, n ≥ 50) is a column for Q1/Q2. Models unchanged: smooth (Q1, Q2) and temporal (Q3, ±3-yr neighbours);
majors: SVI the 15 of §4, WA the 32 of `wa_nw/tables/maj_list.csv`. Q2 decades as before (SVI 2004–14, WA 2007–14 vs 2015–25).

### 6.1 Task 1 — a kernel-robust deep term
- **Kernel variants.** Existing ell20, ell40, ell80, two5_40 (ell40_s2 = MC-noise column, never averaged). Two new variants
  outside their range, same physics, same random numbers (`s11_r3_kernel.py` imports `s11_kernel_mc`): **ell10** (10 km at all
  depths) and **two3_20** (3 km for z < 3 km, 20 km below). V6 = {ell10, ell20, ell40, ell80, two5_40, two3_20}. Stop rules:
  new variants |Σa − 1| ≤ 0.02 and a_d < 0 for ≤ 10 % of pairs; the ell40 re-run (made for the zone tally, §6.3) must reproduce
  the stored ell40 a_s, a_d within 1e-6.
- **Setups** (per-series a_s, a_d, each window its own): NOM = ell40; AVG4 = mean of ell20/40/80/two5_40 (= M2, post hoc in
  round 2); **LOKO_k** = mean of V6 \ {k}, k ∈ V6; AVG6 = mean of V6; and three lapse-chosen setups:
  - Joint-lapse smooth fit per variant v ∈ V6 on the 2–3.5 s and 3.5–6 s cells (cell rule as the headline with each window's own
    cc; own series offsets per window; each window with v's coefficients for that window; shared S_site and D); χ²_v = Σ w r² on
    the identical cell set. **L-SEL** = argmin χ²_v. **L-W** = Σ_v w_v a_v with w_v ∝ exp(−(χ²_v − χ²_min)/2). **L-STA** = per site
    the variant with the smallest χ² contribution of that site's cells (sites with < 30 lapse cells take L-SEL).
  - Diagnostic only: per station the kernel-free in-phase annual amplitude ratio A(3.5–6)/A(2–3.5) with σ, beside each variant's
    deep (ratio of a_d) and shallow (ratio of a_s) predictions.
- **Injections** (waveform level, existing P-weighted stores, full raw chain), generated with each g ∈ V6 plus AVG6 (in-sample
  check of the averaged setup): A = deep step 0.3 % × a_d^g from 2015-01-01; B = shallow annual 0.1 % × a_s^g, peak day 233;
  Pd(e) = 0.3 % × a_d^g in one major; Ps(e) = 0.1 % × a_s^g in one major (each major in turn); Ps(all) = all majors.
- **The six injection rules** (for a setup and a generating kernel): R1 A: Q2 recovery 0.7–1.3; R2 A: |ΔQ1| < 0.03 %;
  R3 B: Q1 leak < 0.25; R4 Pd: median Q3 recovery over majors 0.8–1.1; R5 Ps: median |leak| ≤ 0.35; R6 B: |Δ(Q3 composite)|/0.1
  ≤ 0.35 (new: a mis-scaled shallow annual cycle must not reach the episode composite).
- **Pass.** A setup passes when R1–R6 hold for every generating kernel outside it (LOKO_k: kernel k; L-setups: all of V6; AVG4:
  ell10 and two3_20). Verdicts are given separately for Q1/Q2 (R1–R3) and Q3 (R4–R6). If R5 (or R4) fails also in-sample (NOM
  generated with ell40), that failure is not a kernel error: the Q3 kernel verdict then requires |R5_g − R5_in-sample| ≤ 0.10
  (and |R4_g − R4_in| ≤ 0.10), and the absolute failure is carried to task 2.
- **Winner (order fixed now):** AVG6 if every LOKO fold passes; else L-W, L-STA, L-SEL, AVG4 (first that passes); else none
  (NOM stays, with the measured obstacle). Chosen separately for Q1/Q2 and for Q3.
- **Re-report under the winner**, both regions: Q1, Q2 with G1 leave-one-site-out, G2 family split, G3 scrambled control
  (+ 500-draw date permutation), G4 station pattern (∝ the winner's a_d vs a_s), G5 lapse (predictions from the winner's window
  coefficients), prior ×2/×0.5; Q3 with control, leave-one-site-out, family split, prior. Pass rules and wording as §3–§5.

### 6.2 Task 2 — REDPy findings folded into Q3
- **Measurement** `s11_r3_measure.py`: same rows, windows, 3C realignment and S re-pick as `s11_measure` (series keys identical;
  only series with causality > 1.5 in the existing series table are run). Event weight **w_cc = max(cc_i, 0)**, cc_i = zero-lag 3C
  cc (per component demeaned, unit norm, mean over components) of the realigned event with the leave-self-out P-weighted grand
  stack on S−1…S+5 s (REDPy PLAN §2); scrambled-control windows get the same rule against the grand stack.
- **Estimators.** B0 = the existing raw path (P weights, leave-episode-out reference), recomputed. **B1 (ccw-ref)** = B0 with
  w_cc; 200-draw event bootstrap σ; control episode stacks with w_cc. **B2 (ccw-graph)** = nodes are episode stacks with ≥ 50 data
  events (w_cc); every node pair of the series measured (joint 3C S-anchored stretch 2–5 s, node j against node i); solved with the
  REDPy solve (Weaver a-priori weights scaled by misclosure, 4 Huber passes, zero mean); **σ = 50-draw event bootstrap** (events
  resampled inside every node, all pairs re-measured, re-solved); misclosure σ reported only. Control = scrambled-control episode
  nodes with ≥ 50 control events (≤ n data per episode), same graph and bootstrap. Columns: **B2p** = graph with P weights;
  **B0n50** = B0 at n ≥ 50.
- **Cells.** B1: n ≥ 20, its own cc ≥ 0.35. B2/B2p: graph nodes whose reference-path cc (B1 / B0) ≥ 0.35. All: causal series,
  finite σ, series ≥ 3 cells, σ inflated by the estimator's own scrambled-control χ² (series with ≥ 3 control values, else the
  station median), NOM kernel (the task-1 Q3 winner as a column).
- **Q3 gates per estimator** (rules §4): acceptance σ_ΔD ≤ 0.10 % for ≥ 5 majors; scrambled control; leave-one-site-out; prior
  ×2/×0.5; deep pulse (median recovery 0.8–1.1) and common shallow pulse (median |leak| ≤ 0.35) injected into the estimator's own
  day stacks (ell40 coefficients) and run through that estimator. Reported: σ change vs B0, the deep-pulse over-read (SVI 1.14,
  WA 1.06 in B0), the shallow leak (SVI 0.45, WA 0.60).
- **Adoption rule:** an estimator replaces B0 as the Q3 headline if it passes the deep-pulse rule, the control and leave-one-site-
  out gates, and its 15/32-major composite σ ≤ 1.1 × B0's. Order if several: B2, B1, B2p.

### 6.3 Task 3 — coarse along-strike deep zones (design first)
- **Along-strike kernels.** The ell40 MC is re-run with identical random numbers and additionally tallies the 25–45 km time
  residence of every path (all zones) in 5-km along-strike bins (exact segment–plane intervals; bins −120…+120 km about the
  family centroid, outer bins open). Check: Σ_bins a_d,bin = a_d (≤ 1e-6).
- **Axis:** mean Slab2 strike at the region's family epicentres (SVI N38°W, WA-NW N15°W; computed from the Slab2 depth gradient,
  design only), origin = family centroid.
- **Zones:** 3, boundaries at the 1/3 and 2/3 quantiles of the family along-strike coordinate weighted by each family's Σ1/σ²
  over the headline cells (σ only, no dv/v); outer zones open-ended. a_d,z = Σ of the bins in z (for averaged setups:
  a_d(setup) × the ell40 share of z).
- **Zone model:** the temporal model with D(e) replaced by D_z(e) (3 per episode; prior 1 %). ΔD_z(e) as Q3. Zone composite =
  weighted mean over the region's majors with weights 1/σ²_forecast (fixed by the forecast, before data).
- **Forecast (design only: cell σ and kernels, no dv/v):** σ of each zone composite, median per-major σ, per-year composite σ
  (majors within one calendar year), σ of the zone contrasts.
- **Run rule:** run on data if all three zone composites have forecast σ ≤ 0.10 %. Else forecast 2 zones (halves by weight) and run
  if both ≤ 0.10 %. Else report the forecast and stop.
- **When run** (B0 cells, NOM; task-1 winner and adopted task-2 estimator as columns): zone composites, contrasts, per-year values
  (column). Test "does one part behave differently": χ² of the zone composites about their common weighted mean (full
  covariance, dof = zones − 1).
- **Gates:** scrambled control (zone composites within 2σ AND its χ² test p > 0.05); leave-one-site-out (each zone composite
  within 2σ of the leave-out σ, for every site); waveform-level injections (P store, ell40): deep pulse 0.3 % × a_d,z in all majors
  in zone z only → recovery of z's composite 0.7–1.3 and |leak| into every other zone composite ≤ 0.30 (as a fraction of 0.3 %);
  common shallow pulse 0.1 % × a_s in all majors → |leak| into every zone contrast ≤ 0.35.
- **Wording:** "one part of the fault behaves differently" only if χ² p < 0.01, all zone gates pass and the contrast keeps its sign
  under leave-one-site-out; otherwise zone contrasts are reported with 2σ bounds.

**Resources.** ≤ 12 workers, 1 BLAS thread each, RAM < 80 GB, no GPU. Queue survives the reboot via `scripts/resume_after_reboot.sh`.

**§6 Amendment 1 (Mon 2026-10-12 ~8:00 pm PT, after a 25-series code test on one SVI station, before the full run).** The
in-sample clause of §6.1 applies to every rule R1–R6, not only R4/R5: a rule that already fails when the data are generated and
inverted with the same kernel (NOM, ell40) is not a kernel error, and is then judged on |R_g − R_in-sample| ≤ 0.10 (R2: ≤ 0.01 %).
Why: R2 (A-step leak into Q1 < 0.03 %) already failed at SVI in Revision 1 with the nominal kernel (0.06–0.07 %, ETS-timing
correlation), so the absolute rule could never separate kernel robustness from that known leak. The strict (all-absolute)
verdict is reported beside it in `kerr_verdict.json`; the winner order is unchanged.

## 7. LEVER 4 — distant stations (80–150 km)
Rules in `distant/PLAN.md` (§1–§9 + Amendment 1, written 2026-10-12 8:20–8:45 pm PT); results `distant/RESULTS.md`.

## 8. ROUND 4 — combined model, the 2019-11 candidate, lapse-selection recovery (pre-registered before any round-4 number)
Written Mon 2026-10-12, ~9:30 pm PT, before any round-4 measurement, kernel, fit or injection. Both regions (SVI, WA-NW).
Source: main session ("combine what worked, fix what failed, chase 2019-11"; coarse is fine if new), `round3/RESULTS.md`,
`distant/RESULTS.md` ("what failed / next"). Outputs `round4/{svi,wa_nw}/{tables,logs}`, `round4/figures`; scripts
`svi_science/scripts/s11_r4_*.py`; queue `s11_r4_queue.sh` (markers `round4/<STAGE>_DONE`). No earlier file is changed. No
Merlin call (Fable-quota rule). ≤ 12 workers, 1 BLAS thread each, RAM < 80 GB.

### 8.0 Stages
- **S1 measure** (`s11_r4_measure.py`). (a) Near: the round-3 measurement (`s11_r3_measure.run_series`, same rows, windows, seeds)
  re-run with ONE change: a scrambled-control graph node needs ≥ 10 control events (was 50; data nodes stay ≥ 50 events). Why:
  with 50, the B2 control existed for 1 SVI / 7 WA majors only (median control events per episode is 5). Code check: g_eps,
  g_sig, b1_* reproduce round 3 to ≤ 1e-8 (the data bootstrap draws precede the control draws). (b) Far: the same measurement on
  the far windows (`distant/<region>/tables/win`, far grid −28…+10 s, pre-P window grid-relative as `s11_dist_measure`), only
  far series with causality > 1.5; writes the cc-weighted far store (needed for waveform-level injections on far cells).
- **S2 far kernels** (`s11_r4_farkernel.py`): far pairs, variants ell10 and two3_20 plus an ell40 re-run (same seeds as
  `s11_dist_kernel`), lapse weights from the far P store. Stop rules: |Σa − 1| ≤ 0.02; a_d < 0 for ≤ 10 % of pairs per
  variant; ell40 re-run reproduces the stored far ell40 a_s, a_d within 1e-6.
- **S3 task 3** (`s11_r4_lapse.py`), **S4 task 1** (`s11_r4_combined.py`), **S5 task 2** (`s11_r4_dive.py`, SVI),
  **S6 figures** (`s11_r4_figures.py`).

### 8.1 Task 1 — the pre-registered combined model (both regions)
**Cells.** Near = round-3 B2 (cc-weighted all-pairs graph on episode nodes ≥ 50 events; cell rule b1_cc ≥ 0.35, finite value
and bootstrap σ, series ≥ 3 cells). Far = the same B2 rule on the far series (S1b); far sites as `s11_dist_model.site_map`.
σ = bootstrap σ × √max(1, mean z² of the series' scrambled-control values) (≥ 3 control values; else the station median;
else 1), control values = the S1 B2 control nodes. **Headline H = near B2 + far B2.** Episode terms only for episodes with
cells at ≥ 2 sites; cells of other episodes are dropped (count reported).
**Model.** d = o_series + a_d D(e) + a_s Sc(e) + a_s S_site(t) + R_sta(t) + a_m M(t) + a_b B(t):
- D(e) deep-weighted per episode; Sc(e) one regional shallow term per episode (distant Amendment 1 "T+Sc");
- S_site = 2 annual harmonics + linear + hinges 2007/11/15/19/23 (10 per site; as before);
- **R_sta = receiver term per station, coefficient 1 (not scaled by a_s): annual cos/sin + linear (3 per station), for
  stations with ≥ 20 cells** (a cycle that follows the station is absorbed here);
- M, B = regional mid-crust (5–25 km) and > 45 km terms, slow basis + annual (8 each), free;
- **no prior on D and Sc; instead Σ_e D(e) = 0 and Σ_e Sc(e) = 0 exactly** (null-space solve). Smooth terms (S, R, M, B)
  keep the conditioning prior σ = 1. Numerically singular restricted normal matrix → Moore–Penrose (logged).
- weights 1/σ², σ scaled by max(1, χ²_red); offsets per (family, station, epoch) series eliminated exactly.
**Quantities.** ΔD(e) = D(e) − mean D over the model's other episodes within ±3 yr; per-major σ (median; number ≤ 0.10 %);
composite = equal-weight mean of ΔD over the region's majors present (SVI the 15 of §4; WA the 32 of `maj_list.csv`); χ² of
majors about zero. Kernel setups: **NOM (ell40) = headline**; L-SEL (round-3 lapse choice, two3_20 in both regions) = column
with the identical suite. Columns (natural values only): near only; near B2 + far B0 (n ≥ 20); H without R; H without Sc;
H with the old 1 % prior on D and Sc instead of the constraint.
**Injection suite** (waveform level: day stacks of the near and far cc-weighted stores stretched about tS before stacking, the
graph re-measured; injected − natural). Generating kernels g ∈ V6 ∪ {AVG6} (V6 = ell10, ell20, ell40, ell80, two5_40,
two3_20):
- Pd(e): +0.3 % × a_d^g in ONE major (each in turn) → recovery Δ(ΔD)/0.3;
- Ps(e): +0.1 % × a_s^g at every station in ONE major → leak |Δ(ΔD)|/0.1; Psall: all majors → composite leak;
- one-site: the Ps(e) runs restricted to the stations of ONE site (the injected difference is zero elsewhere), for every site
  with ≥ 5 % of the cell weight;
- B: shallow annual 0.1 % × a_s^g (peak day 233) → |Δ composite|/0.1 and median major |Δ(ΔD)|/0.1;
- A: deep step 0.3 % × a_d^g from 2015-01-01 → composite response reported (a real change; no rule);
- Rann: receiver annual 0.1 % with coefficient 1 (no kernel), all stations, and one site at a time (masked).
**Rules** (NOM; in-sample = g = ell40):
- T1 deep pulse: median recovery over majors 0.8–1.1.
- T2 common shallow pulse: median leak ≤ 0.35 (≤ 0.10 is reported as "leak fixed").
- T3 one-site pulse: median leak ≤ 0.35 for every site with ≥ 5 % weight.
- T4 receiver annual: |Δ composite|/0.1 ≤ 0.35 and median major leak ≤ 0.35, all-station and every one-site run.
- T5 B shallow annual → |Δ composite|/0.1 ≤ 0.35.
- T6 kernel error: T1, T2, T5 for every g ∈ V6 ∪ {AVG6}; a rule that fails in-sample is judged on |R_g − R_in| ≤ 0.10 (§6
  Amendment 1). T6 failing does not block wording; it sets the stated size uncertainty (the range of T1 over g).
**Gates.** G-ctrl scrambled control (same model, constraints, cells and σ, control values of the S1 B2 control nodes): ≥ 90 % of
majors with ≥ 5 control cells within 2σ, and composite within 2σ (column: B1 reference-path control values). G-loso: for every
site, ≥ 90 % of majors within 2σ of the leave-out σ. G-fam: odd/even family composites within 2√(σ²₁ + σ²₂) and ≥ 90 % of
common majors within their combined bound. G-prior: smooth-term prior ×2 and ×0.5 moves every major's ΔD and the composite by
< 0.1σ. G-site: no site > 50 % of the weight (flag). Acceptance: σ_ΔD ≤ 0.10 % for ≥ 5 majors.
**Wording.** A major's ΔD is a "deep-weighted episode change" only if |ΔD| ≥ 3σ and T1–T5, G-ctrl, G-loso, G-fam and G-prior
pass. Otherwise values with 2σ bounds. Composite: 2σ bound unless |composite| ≥ 3σ with the same gates.

### 8.2 Task 2 — SVI 2019-11 deep dive
On H unless stated; round-3 B2 (near, temporal model) and distant F1_B2 values quoted beside.
- (a) **Carriers:** per-cell contributions c_i = (g_2019-11ᵀ Op)_i d_i summed by family, station, site, near/far.
- (b) **Subsets** (same model; terms without data drop out): near only, far only, borehole only (all PB stations, near and
  far), PGC site only, everything except the PGC site.
- (c) **Own scrambled control:** ΔD_ctrl(2019-11) from the B2 control (S1), the B1 control and the B0 control values (same
  cells where they exist). Look-elsewhere: z = ΔD/σ over all model episodes; number with |z| ≥ z(2019-11) vs the Gaussian
  expectation.
- (d) **Leave-one-family-out** (every family with a 2019-11 cell) and leave-one-station-out: range of ΔD(2019-11).
- (e) **Lapse:** graph re-measured at 2–3.5 s and 3.5–6 s (natural, same nodes, cc-weighted store, SVI near). Apparent model
  (all coefficients 1: per-episode term + smooth per site, offsets free, sum-to-zero) → r = ΔE(3.5–6)/ΔE(2–3.5) for 2019-11
  vs the deep prediction (ratio of the contribution-weighted mean a_d, 3.5–6 over 2–3.5 s) and the shallow prediction (same
  with a_s). "deep-like" iff r is ≥ 1σ closer to the deep prediction; "shallow-like" mirror; else undecided. B0 path
  (eps3c_2035/3560) as a column.
- (f) **Arrival times:** ballistic round-1/2/3 values in 2019-11 (clock terms, station S lags, pair dS, region S−P common
  term) and the expected S−P change if ΔD were real: fraction of the direct-S time spent in 25–45 km (straight ray, CVM) ×
  ΔD.
- (g) **Instrument / timing:** (1) clock terms in 2019-11 (round 2/3); (2) per-event realignment lags (`s10_measure.realign`,
  the same windows): median lag of 2019-11 events minus the series median, per station → is the PGC −6 ms removed before
  stretching?; (3) the shift-free estimator: ΔD(2019-11) from eps2p_25 (stretch with a free ±50 ms shift) vs eps3c_25 in the
  B0 temporal model, and the 2019-11 shift2p_25 per station; (4) stretch bias per ms of residual shift in 2–5 s (kernel-free
  arithmetic, lapse-weighted); (5) s10 response epochs and FDSN channel epochs within ±120 d of 2019-11-15 for every station
  carrying ≥ 5 % of the contribution; (6) data days per station 2019-10…12.
  **Timing-artefact rule:** (i) leaving out every station with a 2019-11 clock term ≥ 3 ms gives |ΔD| < 2σ AND within 1σ
  of zero, OR (ii) the shift-free estimator gives |ΔD_2p| < 1σ while ΔD_3c ≥ 3σ (B0 path) AND the stations carrying ≥ 50 %
  of the contribution show a 2019-11 residual shift ≥ 3 ms.
- (h) **Seismicity and tremor:** PNSN tremor within 30 km of the carrier families, per day 2019-10-15…12-15, vs the other 14
  majors (count, days); earthquakes (`g2_local/tables/eq_local.parquet`) M ≥ 2 within 100 km and M ≥ 3 within 150 km of the
  carrier centroid, 2019-10-01…12-31; expected coseismic sign = velocity drop.
- (i) **The fast burst:** q8_rtr streaks of 2019-11 that start 2019-11-08…11-11 (UTC in tables, PT quoted); carrier vs burst
  families; sub-episode split from the cc-weighted day stacks: UTC days ≤ 2019-11-10 ("burst phase") vs ≥ 2019-11-11, each
  measured against the same leave-2019-11-out reference (2–5 s, joint 3C), minus the mean of the series' neighbouring episodes
  (±3 yr, same reference); inverse-variance pooled (σ from a 100-draw day bootstrap), per region near cells.
- **Verdict (fixed now).** "Real deep-weighted change": |ΔD_H| ≥ 3σ AND the H gates pass AND own control within 2σ AND
  leave-one-family-out minimum ≥ 2σ AND near-only and far-only have the same sign, each ≥ 1.5σ AND lapse not shallow-like AND
  not a timing artefact AND no metadata event within ±30 d at stations carrying ≥ 50 % of the contribution. "Instrument or
  timing artefact": the timing rule fires, or a metadata event coincides at stations carrying ≥ 50 % of the contribution.
  Otherwise "undecided", with the blocking tests named.

### 8.3 Task 3 — does the lapse fit pick the right kernel? (both regions)
- Data: the round-3 lapse cells (B0, n ≥ 20, 2–3.5 and 3.5–6 s, own offsets per window, cells usable under all six variants),
  smooth model, χ²_v per variant as round 3; L-SEL = argmin.
- Synthetic: for each generating g ∈ V6, signal s_g = the smooth model fitted to the real lapse cells with g's coefficients
  (prediction, offsets free). Noise (N1) Gaussian σ_i √max(1, χ²_red of that fit); (N2) bootstrap of that fit's residuals
  within station × window. 300 draws per g and noise model. Control: noise only (signal × 0).
- Output: confusion matrix P(sel = v | true = g); uniform-prior confidence P(true = v* | sel = v*), v* = the data's L-SEL;
  class K = {v: Δχ²_data(v) ≤ 4}: P(sel ∈ K | true ∈ K) and max_{g ∉ K} P(sel ∈ K | true = g).
- **Adoption (Q1 only):** adopt v* if P(sel = v* | true = v*) ≥ 0.8 AND max_{g ≠ v*} P(sel = v* | true = g) ≤ 0.2, in both
  noise models. Else adopt the class K (mean coefficients of K) if P(sel ∈ K | true ∈ K) ≥ 0.8 AND max_{g ∉ K} P(sel ∈ K | g)
  ≤ 0.2 (both noise models). Else nothing is adopted. Confidence = the uniform-prior posterior above.
- Under the adopted setup: Q1, Q2 on the round-3 near cells with the round-3 gates (`s11_r3_kerr.rereport`: G1–G5, prior,
  date permutation, Q3) and the round-3 kernel-error rows of `kerr_rules.csv` for that setup; expected recovery and leak =
  posterior-weighted mean over the generating kernels. Column: F1 (near + far, M/B free) Q1 under the adopted kernel.
  Wording for Q1 as §3 (deep annual change only if Q1/σ ≥ 3, G1–G3 pass, in-sample B leak < 0.25; G4/G5 supporting).

**§8 Amendment 1 (Mon 2026-10-12 ~9:45 pm PT; after S1/S2 and a natural-value smoke test of H on SVI, before any injection,
lapse-selection or dive number).** The smoke test shows that with the exact constraint (no prior) the regional shallow term
Sc(e) trades off with D(e): SVI σ per major 0.199 %, composite σ 0.177 % (H) vs 0.111 / 0.041 % without Sc, and 0.149 /
0.067 % with the old 1 % prior on D and Sc. So the precision of the round-3/distant "T+Sc" model came largely from its 1 %
prior on Sc, which is also why its prior test failed. Additions (H stays the pre-registered headline; nothing removed):
(1) the full task-1 suite (T1–T6, every gate) is also run for **HnoSc = H without Sc** (receiver term, M/B, constraint kept),
under NOM and L-SEL; (2) the 2019-11 dive runs sections (a)–(d) for H and for HnoSc and states the verdict for both (rule
unchanged); the carrier families/stations used in (f)–(i) come from the form whose ΔD(2019-11) has the candidate's
positive sign (H if positive, else HnoSc), and the timing rule (i) is evaluated in that form.

**§8 Amendment 2 (Mon 2026-10-12 ~10:00 pm PT; after a code test of the dive on the real SVI data — its numbers were seen,
so these changes are POST-TEST and are labelled so in the results).** (1) §8.2(e) apparent lapse model: the coefficient-1
smooth-per-site term is exactly degenerate with the coefficient-1 episode term (code test: σ ≈ 1 %); the apparent model is
now the episode term + series offsets + sum-to-zero only. (2) §8.2(i) sub-episode split: a day bootstrap of a 2–3-day part
is degenerate; σ of a part = the 2019-11 cell's event-bootstrap σ × √(events in 2019-11 / events in the part) (day
bootstrap kept as a column). (3) New timing test (g7), because the code test showed residual 2019-11 shifts after the
per-event realignment (2-parameter fit: PGC −5 ms, B009 −2.5 ms; stretch bias −0.029 % dv/v per ms of delay in 2–5 s):
every episode-node stack is realigned on the direct S (s10 realign: tS−1…tS+2 s, 3C, sub-sample) to the series grand stack
before the graph stretch, near and far; ΔD(2019-11) with and without (same cells, H and HnoSc). Timing-artefact rule
extended: (iii) node realignment brings |ΔD| below 1σ while the same cells without it give ≥ 2σ (candidate's model form).

## 9. ROUND 5 — close the failing round-4 gates, settle the WA yearly cycle and SVI 2019-11 (pre-registered before any round-5 number)
Written Mon 2026-10-12, ~10:25–10:30 pm PT, before any round-5 fit, injection response or measurement. Both regions (SVI,
WA-NW). Source: main session (round-4 "what failed / next"; user: persevere, coarse is fine if new). Outputs
`round5/{svi,wa_nw}/{tables,logs}`, `round5/figures`; scripts `svi_science/scripts/s11_r5_*.py`; queue `s11_r5_queue.sh`
(markers `round5/<STAGE>_DONE`, `round5/QUEUE_DONE`). No earlier file is changed. Reused as is: the round-4 cells (near B2 +
far B2 = H, `load_H`), the round-4 waveform-level injections (`round4/<r>/tables/r4_injected.npz`, every run of §8.1, re-aligned
to H by (series, episode) and checked), all kernel tables. New waveform measurement only in task 4 (graph re-measured on the
round-3 cc-weighted stores at extra lapse windows, PGC-site series only) and new kernel windows only for those series (same
MC seeds as `s11_r3_kernel`; the 2–5 / 2–3.5 / 3.5–6 s coefficients must reproduce the stored ones within 1e-6). ≤ 12
workers, 1 BLAS thread each, RAM < 80 GB. No Merlin call (Fable-quota rule).

### 9.0 Common definitions
- **Units (one shallow term per site-station).** A station with ≥ 20 cells in the model's cell set is its own unit; the
  other stations of a site share one "site-rest" unit. Unit basis by size (fixed now): ≥ 40 cells → the full smooth basis
  STERMS (c1, s1, c2, s2, lin, hinges 2007/11/15/19/23); 10–39 cells → (c1, s1, lin); < 10 cells → no unit term (series
  offsets only). Within each unit block, exactly-null directions of the weighted, series-demeaned columns are removed (SVD,
  singular values < 1e-6 × the block maximum); this does not change any fit. **No S_site AND R_sta pair any more.**
- **Unit coefficient.** "Ua" = a_s × basis (kernel shallow term); "U1" = 1 × basis (receiver-type term).
- **D(e).** Per episode, coefficient a_d, exact Σ_e D(e) = 0 ("D0"). "Dp" adds the exact constraints
  Σ_e D(e) b_k(t_e) = 0 for b_k ∈ {lin, h2007, h2011, h2015, h2019, h2023, c1, s1} evaluated at the episode mid-dates of
  the model (dependent/empty rows dropped by the null-space solve): D(e) then carries no slow or annual part (round-4 "next
  3"; the smooth terms carry them).
- **M, B** regional mid-crust / > 45 km terms, a_m / a_b × (slow + annual), as round 4.
- **Exact basis = no ridge on any term.** Solve on the null space of the constraints; condition number of the restricted
  normal matrix logged; > 1e12 → pseudo-inverse with relative tolerance 1e-10 (logged, then G-prior includes tolerance
  ×100 and ÷100).
- **G-prior (r5)** for an exact form: add a ridge σ_p = 1 % and σ_p = 2 % on every non-episode term; every major's ΔD
  and the composite must move < 0.1σ (σ of the exact form). For a penalised form: σ_p ×2 and ×0.5 as before.
- **Calibrated-penalty column ("Pcal")**: the selected exact form with a ridge per term group (units, M, B, Sc-spline)
  σ_p(group) = max(RMS of that group's exact-fit coefficients on the data, RMS on the B2 scrambled-control cells). Judged by
  ×2 / ×0.5.
- **Suite.** T1–T5 exactly as §8.1 (from the round-4 injections; in-sample generator = the fit setup); T3 for every site
  with ≥ 5 % weight; T4 all-station + one-site receiver annual; G-ctrl (B2 control, B1 column), G-loso, G-fam, G-site as §8.1;
  acceptance σ_ΔD ≤ 0.10 % for ≥ 5 majors. **T6 is replaced by the band of task 5** (no pass/fail on one variant).
- **Setups.** NOM (ell40) and L-SEL (two3_20) both reported for every form. DECISION 27: the round-5 headline setup is
  L-SEL for Q3 as well as Q1 (the lapse fit's choice was validated in round 4 task 3 for both regions); NOM stays a full column.

### 9.1 Task 1 — prior test with one shallow term per unit and an exact basis
- **Forms:** E1 = Ua + D0; E2 = Ua + Dp; E3 = U1 + D0; E4 = U1 + Dp (each with M, B; no Sc; exact basis). Columns:
  r4-HnoSc reproduced (must equal `round4/<r>/tables/r4_gates.csv` HnoSc σ/major and composite to 1e-6), Pcal of the selected
  form.
- **Each form:** natural ΔD per episode, σ per major (median, number ≤ 0.10 %), composite; full suite in both regions under
  NOM and L-SEL.
- **Selection rule (fixed now):** score = number of passes of {T1, T2, T3, T5, G-prior, G-ctrl, G-loso, G-fam} summed over
  the 4 (region × setup) combinations; highest score is "selected"; tie → smaller mean (over regions) L-SEL σ per major;
  still tied → E2. If no exact form passes G-prior in a combination, Pcal is judged there by ×2 / ×0.5.
- **Pass (task 1):** the selected form passes G-prior < 0.1σ (majors and composite) in both regions under L-SEL, with
  T1–T3, T5, G-ctrl, G-loso, G-fam also passing. T4 reported (Dp is expected to help; it is not required).

### 9.2 Task 2 — regional shallow term smooth in time (no per-episode trade-off)
- **Form:** the task-1 selected form + Sc_spl = a_s × regional cubic B-spline in time, interior knots every k years over
  the episode span (k = 2, 3, 4), its columns first made orthogonal (unweighted, at the model's episode mid-dates) to the
  unit basis functions STERMS, so it carries only regional 2–4-yr shallow content the unit terms cannot; exact basis.
  ("Tied to the shared annual + slow parts" reduces to this: the shared annual and 4-yr slow parts are already in every
  unit term, so a regional copy would be an exact duplicate.)
- **Goal (both regions, L-SEL and NOM):** median T2 leak ≤ 0.35 AND σ per major ≤ round-4 HnoSc of the same setup
  (SVI NOM 0.111, L-SEL 0.082; WA NOM 0.128, L-SEL 0.124 — from `round4/*/tables/r4_gates.csv`). Full suite reported.
- **Adoption:** among k meeting the goal in all 4 combinations, the one with the smallest mean L-SEL σ per major; none →
  report each k's failing item. The round-4 per-episode Sc (H) is reproduced as a column.

### 9.3 Task 3 — the WA yearly cycle (Q1): deep, shallow or mixed
- **Cells:** (a) round-4 Q1 cells (near head cells, B0 path, `R.head_cells`) — the subject of the round-4 "deep by rule";
  (b) near + far (`s11_dist_model.all_cells`, the F1 set: contains the PGC-site stations as far stations). Kernels two3_20
  (far pairs: `kernel_coef_far_r4`). SVI run as a column (same code).
- **Models** (smooth; D = a_d × (slow + c1, s1)):
  Q1-r4 = round-4 form (D + S per site, ridge 1; must reproduce WA 0.151 ± 0.043 on (a));
  **Q1-R (headline) = D + one U1 term per unit (receiver-type, coefficient 1, unit basis by the 9.0 size rule) + M, B,
  exact basis**; Q1-Ua = the same with Ua units; Q1-RL = Q1-R on the joint-lapse cells (2–3.5 and 3.5–6 s, own offsets per
  window, U1 shared by both windows of a station, D with each window's a_d).
- **Groups:** P = PGC site (PB.B009, PB.B010, PB.B011, CN.PGC); B = other PB boreholes; O = all other stations.
  Per group: Q1-R on the group's cells; contribution share of the group to the pooled in-phase Q1 (operator rows);
  **lapse test per group** = kernel-free in-phase annual amplitude (at the pooled Q1-r4 peak day) at 3.5–6 s over 2–3.5 s
  (lapse cells, near + far), vs the group's predictions: deep (mean a_d ratio, two3_20), shallow (mean a_s ratio) and
  receiver-type (1.0). Call per group: "deep-like" iff ≥ 1σ closer to deep than to both others; "shallow-like" iff ≥ 1σ
  closer to shallow or receiver than to deep; else undecided. Pooled call the same way. Per-station amplitudes as a
  column (figure).
- **Gates on Q1-R:** G1 (leave-one-site-out), G2 (family split), G3 (scrambled control < 2σ), date permutation (500), G-prior
  (r5 definition).
- **Verdict (fixed now, per cell set; headline = (b) for WA because only (b) contains the PGC site):**
  DEEP iff Q1-R/σ ≥ 3 AND G1–G3 pass AND the pooled lapse call is not shallow-like AND no group with ≥ 25 % share is
  shallow-like AND the cycle is ≥ 2σ in ≥ 2 groups with peaks within 60 d of the pooled peak.
  MIXED iff not DEEP and (one group with ≥ 25 % share is deep-like while another is shallow-like, OR Q1-R/σ ≥ 3 with a
  shallow-like pooled lapse call).
  SHALLOW iff not DEEP/MIXED and (Q1-R/σ < 2 — the cycle is absorbed by per-station receiver terms — OR the pooled lapse
  call is shallow-like) and no group with ≥ 25 % share is deep-like.
  Otherwise UNDECIDED, blocking items named.

### 9.4 Task 4 — SVI 2019-11 per instrument at the PGC site (co-located: CN.PGC vault, PB.B009/B010/B011 boreholes)
- **Per-instrument values** (SVI near B2 cells, 2–5 s): (i) apparent ΔE_i(2019-11) = episode model with coefficient 1, series
  offsets, sum-to-zero, on instrument i's cells only (σ inflated as H); (ii) model ΔD_i = task-1 selected form (L-SEL and NOM)
  on all non-PGC-site cells + instrument i only; (iii) sub-episode parts per instrument (day stacks, `split` machinery,
  same reference): whole episode, UTC ≤ 11-10, 11-11…11-21, 11-22…11-25 (carrier families fired), ≥ 11-26; σ by the
  Amendment-2 event-count scaling. Scrambled control per instrument (B2 control values, same model).
- **Call (fixed now)** on (i): common-value fit over the 4 instruments → m ± σ_m, χ² (3 dof).
  SITE-WIDE ("the site or real") iff χ² p ≥ 0.05 AND m ≥ 2.5σ_m AND ≥ 3 of 4 instruments > 0.
  ONE INSTRUMENT j iff ΔE_j ≥ 2.5σ_j AND the common value of the other three is within 1.5σ of 0 AND
  |ΔE_j − m_{−j}| ≥ 2.5σ_diff. Otherwise MIXED/UNDECIDED, naming the instruments ≥ 2σ.
  INSTRUMENT ARTEFACT iff ONE INSTRUMENT and that instrument has a metadata event within ±30 d.
  Column: vault − borehole-mean contrast (PGC − mean(B009, B010, B011)) ± σ (a near-surface change would load the vault).
- **Metadata per instrument:** FDSN channel epochs (EarthScope for PB, NRCan for CN; start/end/sensor/sample rate) and
  s10 response epochs within ±120 d of 2019-11-15; data days per instrument 2019-10…2019-12.
- **Lapse-window profile at PGC:** graph re-measured (same nodes, same stores, 50-draw event bootstrap) at 2–3, 3–4, 4–5,
  5–6 s for the PGC-site series; apparent ΔE(2019-11) per window, pooled over the site and per instrument. Kernel
  coefficients for these windows from re-simulated kernels (ell40, ell10, two3_20; s11_r3_kernel seeds; reproduction check
  above). Four one-parameter shapes fitted by weighted least squares: deep k·ā_d(w), shallow k·ā_s(w), uniform k
  (receiver/whole-medium), timing k/t̄(w) (t̄ = mean lapse after S of the window; a residual delay gives apparent dv/v ∝ 1/t).
  Best shape = lowest χ²; "decisive" iff Δχ² ≥ 4 to every other shape. Headline kernels two3_20; ell40 column.

### 9.5 Task 5 — kernel error as a systematic band (replaces the T6 pass/fail)
- **Plausible set P per region** (fixed now from the round-4 data Δχ², `round4/*/tables/lapse_verdict.json`):
  P = {v ∈ V6 : Δχ²_data(v) ≤ 9} → SVI {two3_20, two5_40, ell40}; WA {two3_20, two5_40}. Full band over V6 ∪ {AVG6} as a
  column.
- **Q3:** for the task-1 selected form (and the task-2 adopted form if any): composite ± σ and per-major ΔD under every fit
  setup v; the deep-pulse recovery T1(headline setup, generator g) for g ∈ P. Band = [min, max] over P of the composite;
  size band = 2σ bound / min_{g∈P} T1 … / max T1.
- **Q1:** Q1-r4 and Q1-R under every v ∈ P (cell sets a and b); band [min, max] of Q1, Q1/σ and peak day.
- **"Conclusion holds across the band"** iff for every v ∈ P the same statement is true: Q3 — composite within 2σ of zero
  (or the same sign at ≥ 3σ); Q1 — the same side of 3σ, and the task-3 call unchanged when the task-3 verdict is recomputed
  with v's coefficients.

### 9.6 Stages (queue `s11_r5_queue.sh`, markers in round5/)
S1 `s11_r5_model.py` per region (tasks 1, 2, 5-Q3) → MODEL_<R>_DONE; S2 `s11_r5_q1.py` per region (task 3, 5-Q1) →
Q1_<R>_DONE; S3 `s11_r5_pgc.py` (task 4, SVI) → PGC_DONE; S4 `s11_r5_figures.py` → FIGURES_DONE; QUEUE_DONE.
Figures: dv/v points coloured by weight or cc, no error bars, line graphs.

**§9 Amendment 1 (Mon 2026-10-12 ~10:35 pm PT; POST-CODE-TEST, labelled so).** A code test of stages forms/sc on the first
40 SVI series (268 cells, PGC site only; outputs in `round5/svi/tables/test/`, not results) showed that every exact form
still moves by 0.17–0.56σ when a 1 % ridge is added, and a one-group-at-a-time diagnostic on that subset (E2, two3_20) put
the cause in the regional > 45 km term B (ridge on B alone 0.17σ; on M alone 0.085σ; on the unit terms 0.003σ): with
a_b ≈ 0.02 the data cannot pin B (σ_θ ≈ 0.6 %), so any prior on it moves ΔD. Additions (nothing removed, numbers of the
full run not seen):
(1) a per-group G-prior diagnostic (ridge 1 % / 2 % on U, M, B, Scs one at a time) is reported for every exact form;
(2) new forms **E2b, E4b** = E2, E4 with the > 45 km layer folded into the deep term (coefficient a_d + a_b on D(e), no B
term; i.e. "deep" = everything below 25 km) — eligible for the §9.1 selection with the same rule (tie order E2, E2b, …);
**E2m** = E2b also without M (column only, not eligible: no mid-crust injection exists to test the leak this would open);
(3) task 3: Q1-R is the round-4 Q1 form with ONLY the per-station change (D smooth + U1 per unit; no M/B, as Q1-r4 has
none); Q1-R+MB (M, B free; or the "b" treatment if the selected form is a "b" form) is a column. G-prior for Q1-R is the
same r5 definition.

## 10. ROUND 6 — the WA deep-pulse recovery, and does the per-episode deep change scale with episode size? (pre-registered before any round-6 number)
Written Mon 2026-10-12, ~11:04 pm PT, before any round-6 fit, regression or injection response. Both regions (SVI,
WA-NW). Source: main session (round-5 "next 1"; user: persevere, coarse is fine if new; the σ ≈ 0.06 % per major makes the
size question answerable). Outputs `round6/{svi,wa_nw}/{tables,logs}`, `round6/figures`; scripts
`svi_science/scripts/s11_r6_*.py`; queue `s11_r6_queue.sh` (markers `round6/<STAGE>_DONE`, `round6/QUEUE_DONE`). No earlier
file is changed. **No new waveform measurement:** the round-5 cell set H (`load_H`), the round-4 waveform-level injections
(`round4/<r>/tables/r4_injected.npz`, re-aligned as round 5) and the round-5 model class (`s11_r5_common.Model5`) are reused
as is. ≤ 12 workers, 1 BLAS thread each, RAM < 80 GB. No Merlin call (Fable-quota rule).

### 10.0 Common definitions
- Everything in §9.0 holds (units, U1, exact basis, G-prior r5, suite T1–T5, G-ctrl/G-loso/G-fam/G-site, acceptance).
- Setups: L-SEL (two3_20) = headline, NOM (ell40) = full column (decision 27 unchanged).
- New D constraint sets (each with the exact sum-to-zero Σ_e D(e) = 0, at the model's episode mid-dates):
  - **Dy** = D ⊥ {c1, s1} only (the yearly pair; round-5 "next 1"). Form **E4y** = U1 + Dy + M + B.
  - **Dl** = D ⊥ {lin, c1, s1} (drops the five hinges). Form **E4l**.
  - **Dh** = D ⊥ {lin, h2011, h2015, h2019, c1, s1} (drops the two end hinges h2007, h2023). Form **E4h**.
  - **E4y+Scs2** = E4y + the round-5 adopted smooth regional shallow spline (k = 2 yr).
  - Round-5 E4 (= Dp, all 8) and E3 (= D0) are reproduced as columns (must equal round-5 σ/major, composite, T1 to 1e-6).

### 10.1 Task 1 — why WA L-SEL T1 is 0.75, and the fix
- **Cause analysis (diagnostic, every form, both regions, both setups), per major e:**
  (i) h_ee = leverage of e in the constraint basis [1, PROJ-set] evaluated at the model's episode mids; the projection-only
  prediction of the ΔD recovery = (δ_e − Pδ_e)(e) − mean over ±3-yr neighbours of (δ_e − Pδ_e);
  (ii) R_struct = recovery of a model-consistent pulse (cell value a_d(setup) × 0.3 % on e's cells only) through the form's
  operator;
  (iii) R_wave = recovery of the round-4 waveform-level pulse Pd|gen=setup|e (the T1 input), per major.
  The cause is named by which part carries the E4 − E3 shortfall: projection (ii ≈ iii and ii(E3) ≈ 1), or not.
  Spearman ρ(h_ee, R_wave) over majors is reported.
- **Pass rule (task 1, round 6), per form, in all 4 region × setup combinations:** T1 ∈ [0.8, 1.1] AND G-prior < 0.1σ AND
  T2, T3, T4 (all-station and one-site, composite and major median ≤ 0.35), T5, G-ctrl (B2), G-loso, G-fam pass AND
  acceptance (≥ 5 majors with σ ≤ 0.10 %). G-site and the B1 column are reported, not required (as round 5).
- **Adoption (fixed now):** the first form in the order **E4y, E4l, E4h, E4y+Scs2** that meets the pass rule is adopted as
  **E4′**. None → the form with the most passes (tie: order above) is reported as "best, not passing", its failing items named;
  round-5 E4 then stays the headline for task 2, with every size quoted divided by its measured T1.
- Reported for every form: natural ΔD per episode, σ per major (median, count ≤ 0.10 %), composite ± σ, all gates and
  rules, per-major T1, and for E4′ the Q3 band over the plausible kernels (§9.5 definition: composite per fit kernel in P,
  T1 over generators in P, size band).

### 10.2 Task 2 — does ΔD(e) scale with episode size?
- **ΔD source (fixed now):** headline = E4′ under L-SEL; columns: E4′ NOM, E4′+Scs2 L-SEL (the smooth-regional-shallow form;
  if E4′ already contains Scs2, the column is E4′ without it), round-5 E4 L-SEL. ΔD(e), σ(e) for every model episode with a
  ±3-yr neighbour (as the round-5 tables).
- **Size measures** (only episodes of the model with a finite value enter that measure's test):
  - SVI: **G = S1 Tier A GNSS offset D (mm)** (`s1/episode_geodetic.csv` `D`; 2014+) — **SVI headline**;
    M0 = S1 Tier B regional moment (`tierB_M0_Nm`, per 10^18 N·m; regional gate PASS); T = PNSN tremor epicentres within
    10 km of a family (`b_ets_record/episodes_b.csv` `n_trem10`, per 1,000; 2010+); L = S3 LFE count (`s3/episodes_all.csv`
    `n_headline`, per 1,000; 2014+); Du = tremor duration (`episodes_b.csv` `duration_d`, per 10 d); Pc = S3 patch coverage
    (share of the 12 S0 groups with ≥ 3 headline LFEs in the episode, from `s3/events_svi.parquet`); full-record columns
    from the dv/v family set itself (`active_days.csv`, main + ext, restricted to the families in H): Fs = Σ sumP (per 1,000)
    and Fc = share of H's families with sumP ≥ 1 in the episode.
  - WA: **T = PNSN tremor `n_trem10` (per 1,000; 2009-08+) — WA headline** (no per-WA-episode GNSS series exists);
    G_WA = S1 Tier A WA-station offset `D_WA` (mm) of the S1 (SVI) episode whose window overlaps the WA episode window
    (column; partial); L = NW LFE count counted once P ≥ 0.75 (`episodes_b.csv` `lfe_nw_once075`, = S3 for E3; all eras);
    Du; Pc = share of NW 10-km along-strike bins (bins holding ≥ 1 % of all S3 NW headline events) with ≥ 3 headline events
    (`s3/events_nw.parquet`, E3 only); Fs, Fc from `wa_nw/tables/inv/active_days.csv`.
- **Regression (per region × ΔD source × measure):** ΔD(e) = a + b·size(e), weighted least squares, w = 1/σ(e)²;
  σ_b from the WLS covariance × max(χ²_red, 1). Columns: GLS with the full ΔD covariance (g_eᵀ C g_e'), majors only,
  standardised slope b·SD(size) (% per SD), leave-one-episode-out (max |Δb| / |b|, sign kept?).
- **Scrambled control:** sizes shuffled among the episodes of that test, 5,000 draws, seed 0; p = share with |b_perm| ≥ |b|
  (two-sided); the 95th and 99th percentiles of |b_perm| are reported.
- **Recovery (waveform-level):** a deep change proportional to size, 0.3 % × size/size_max on every major with a size, is
  built by superposing the round-4 Pd|gen=setup|e injections (input − natural); recovery = fitted slope / the slope of the
  true ΔD pattern (same dvec functional, same episodes, same weights). Slopes and bounds are quoted as measured and divided
  by this recovery when it lies outside 0.9–1.1.
- **What would have been seen:** b_det = 3σ_b (and the scrambled-control 99th percentile of |b_perm|), translated into a
  ΔD difference between a large and a small episode = b_det × (size q90 − size q10).
- **Top third vs rest (the change at the largest episodes):** among the episodes with a size, the top third by size
  (ceil(n/3)) against the rest:
  (a) **during-episode contrast** = weighted mean ΔD(top) − weighted mean ΔD(rest) (w = 1/σ²), σ from the full covariance;
  (b) **before → after step** at each episode = mean D of model episodes in (e, e + 3 yr] − mean D in [e − 3 yr, e)
  (both sides non-empty), weighted mean over top minus rest, σ from the full covariance.
  Scrambled control for both: random top sets of the same size, 5,000 draws, seed 0, two-sided p.
- **Verdict per region (headline measure, headline ΔD source; fixed now):**
  **SCALES** iff |b| ≥ 3σ_b AND scrambled-control p < 0.01 AND the E4′+Scs2 column has the same sign at ≥ 2σ AND
  leave-one-episode-out keeps the sign; **HINT** iff not SCALES and |b| ≥ 2σ_b and p < 0.05;
  else **NO SCALING SEEN**, quoted as the 2σ bound |b| + 2σ_b (÷ recovery) and the b_det translation. Column measures get the
  same numbers; they are 7–8 tests per region, so a column at p < 0.05 alone is a hint at most (stated as such).
- **Figure (fixed now):** ΔD against size per region × measure, points coloured by weight 1/σ², a thin line joining the
  points in size order, the WLS line; no error bars.

### 10.3 Stages (queue `s11_r6_queue.sh`, markers in round6/)
S1 `s11_r6_model.py` per region (task 1: forms, cause analysis, suite) → MODEL_<R>_DONE; S2 `s11_r6_model.py --stage adopt`
→ ADOPT_DONE (`round6/adopted.json`); S3 `s11_r6_size.py` per region (task 2) → SIZE_<R>_DONE; S4 `s11_r6_figures.py` →
FIGURES_DONE; QUEUE_DONE.

**§10 Amendment 1 (Mon 2026-10-12 ~11:20 pm PT; POST-RESULT for task 1, labelled so; written before any task-2 number was
looked at).** Task-1 result: no candidate passes in all 4 combinations (E4y/E4l lift WA L-SEL T1 to 0.82/0.81 but fail
G-prior 0.15/0.14σ there, the shift sitting in the unit terms U; E4h keeps the prior but T1 0.77). Cause measured: the D
constraints remove the part of a one-episode pulse that the slow basis can fit (E4 per-major WA recovery tracks the
leverage h_ee, ρ −0.81; model-consistent recovery 0.876 = the projection-only prediction 0.861), on top of a waveform-level
loss already present without any constraint (E3 0.885). Since the prior shift comes from D trading against the unit slow
terms, the next adaptation removes the shared slow directions on BOTH sides instead of only from D:
- **E4lu** = units with {c1, s1, c2, s2, lin} (≥ 40 cells) / {c1, s1, lin} (10–39), D ⊥ {lin, c1, s1};
- **E4yu** = units with {c1, s1, c2, s2} (≥ 40) / {c1, s1} (10–39), D ⊥ {c1, s1}.
(The ±3-yr neighbour difference in ΔD already removes most of a slow drift; T2/T3 measure whatever leaks.) Same pass rule
(§10.1). Adoption order E4lu, E4yu → "E4′ (amendment)". The pre-registered task-2 headline stays round-5 E4 (the §10.1
"none" branch, fixed before this amendment); the amendment form, if adopted, is added to task 2 as a labelled column
(L-SEL), with the same tests. Stages: `s11_r6_model.py --stage amend` per region → AMEND_<R>_DONE; `--stage adopt_amend`
→ ADOPT_AMEND_DONE (`round6/adopted_amend.json`); `s11_r6_size.py --amend` per region → SIZE_AMEND_<R>_DONE.

## 11. ROUND 7 — front timing, per-episode zones with E4, joint arrival-time + coda, the WA large-episode hint (pre-registered before any round-7 number)
Written Tue 2026-10-13, ~9:26 am PT, before any round-7 measurement, fit, injection or contrast. Both regions (SVI, WA-NW).
Source: main session (user: persevere, coarse is fine if new; finish the open ideas: Merlin test 3, per-episode zones,
idea-4 middle ground, independent test of the round-6 WA hint). Outputs `round7/{svi,wa_nw}/{tables,logs}`,
`round7/figures`; scripts `svi_science/scripts/s11_r7_*.py`; queue `s11_r7_queue.sh` (markers `round7/<STAGE>_DONE`,
`round7/QUEUE_DONE`; PID `round7/queue.pid`). No earlier file is changed. **No new download; no new event-level waveform
measurement:** the round-3 / round-4 cc-weighted day-stack stores (`round3/<r>/tables/store_ccw[_ext]`,
`round4/<r>/tables/store_ccw_far`), the round-5 cell set H (`s11_r5_common.load_H`), the round-4 injections, the round-3
along-strike kernel bins and the ballistic round-3 lag tables are reused. ≤ 12 workers, 1 BLAS thread each, RAM < 80 GB.
No Merlin call (Fable-quota rule). dv/v in %; "scrambled control" as defined per task; L-SEL (two3_20) = headline setup,
NOM (ell40) = column; E4 = round-5 form (U1 units, D ⊥ PROJ, M, B, exact basis); majors = the region lists (`R4.majors()`).

### 11.1 Task 1 — front timing (Merlin test 3)
- **Series:** every series of H (near and far; SVI main + ext). **Day values:** for each series, reference G = the
  weight-summed stack of all its stored days; each stored day stack is measured against G with the 2–5 s joint-3C stretch
  (`Ref3`); dv/v_day = −100 ε; |ε| ≥ RAIL → missing. Weight w_day = mean_c Ws[day, c] × max(cc_day, 0)². The stored
  time-of-day scrambled-control day stacks (Sc) are measured the same way (**time-of-day control**).
- **Front time per (family, major):** (A) S3: SVI unit = family, time = w0 + lfe_t_first + onset (`s3/family_onsets.csv`,
  `episodes_kinematics.csv`, UTC days); NW: the S3 source onset given to each member family (S3 `unit_fid` → S11 fam via
  `wa_nw/tables/inv/families.csv` family_id). (B) own: the date at which the family's cumulative sumP in that major
  (`active_days.csv`) reaches 10 % (the S3 10 %-quantile definition on the S11 detection days). **Headline = A where it
  exists, else B.** Columns: A only, B only. S3 episode names are matched to the S11 major names; where they differ the
  S11 major whose [start, end] contains the S3 window mid is used.
- **Lag** = (day − UTC date of the front time) in days. Bins: ≤ −1, 0, 1, 2, 3, 4–6, ≥ 7. Only days of majors enter the
  curves; the between baseline b_s = weighted mean of the series' non-major days.
- **Curves:** (i) absolute = weighted mean over (series, major, day) of dv/v_day − b_s per bin; (ii) **within-episode
  shape** = dv/v_day − weighted mean of the same (series, major) days, per bin (immune to any per-episode level, slow or
  yearly term). **Headline contrast C** = shape(lag ≥ 3) − shape(lag 0–1); **slope b** = WLS slope of the day shape values
  on lag clipped to [−1, 10] (% per day).
- **σ:** bootstrap over majors (1,000 draws, seed 20261012), recomputing everything.
- **Front-scrambled control (headline control):** 500 draws; each (family, major) front time replaced by
  start(e) + (front(f, e′) − start(e′)), e′ drawn from the other majors with a front for f (start = episode unit start);
  p = share of |C*| ≥ |C| (and for b). Time-of-day control: C and b on Sc days with the true fronts.
- **Bias on ΔD:** L̄_e = weighted mean lag of major e's days (clipped [−1, 10]); bias_e = b_up × (L̄_e − L̄_all majors), with
  b_up = |b| + 2σ_b; compared with round-5 E4 L-SEL σ_e and σ_composite. **"Front timing biases ΔD"** iff any |bias_e| ≥
  0.5 σ_e or the composite bias ≥ 0.5 σ_comp; else "no bias at the 0.5σ level".
- **Verdict per region:** **FRONT-LOCKED CHANGE** iff |C| ≥ 3σ AND front-scrambled p < 0.01 AND the time-of-day control
  |C_c| < 2σ_c; else **NONE SEEN**, bound |C| + 2σ (and |b| + 2σ_b). Both regions FRONT-LOCKED with the same sign →
  "seen in both". Columns (reported, not gates): front source A / B; near vs far; series with L-SEL a_d above / below the
  region median; a bin-stack check (days of each (series, bin) stacked over majors and measured against G) vs the day
  pooling.

### 11.2 Task 2 — along-strike zones with E4, one major at a time; slipped vs not slipped
- **Zones:** round-3 bounds (`round3/<r>/tables/zones_def.json`: SVI 2, WA 3). Zone shares from the round-3 ell40
  along-strike bins (`s11_r3_zones.zone_coefs`, setup a_d × ell40 share); cells without bins (far) get the mean share of the
  same family's near cells, else are dropped (count reported).
- **Model E4Z:** E4 with D replaced by D_z(e) (coefficient a_d × share_z); each zone's D sum-to-zero and ⊥ PROJ at the
  episode mids; U1 units, M, B as E4. Setup L-SEL headline, NOM column.
- **Per major:** ΔD_z(e) ± σ (±3-yr neighbour functional per zone), contrasts z_i − z_j, equal-zones χ² and p.
  **Pooled:** IVW composite per zone, contrasts, equal-zones χ² p; and Σ per-major χ² (dof n_major × (n_z − 1)).
- **Slipped labels:** (S3) patches = SVI S0 groups / NW 10-km bins (the round-6 Pc patches), each assigned to the zone of
  the median along-strike position of its S3 headline events projected on the round-3 axis; zone z slipped in e iff
  ≥ 50 % of its patches have ≥ 3 S3 headline events in e. Column (F) for every major: ≥ 50 % of the zone's H families
  have sumP ≥ 1 in e. **Slipped − not** per major with both classes = mean ΔD over slipped zones − mean over the rest
  (σ from the covariance); pooled IVW; label-permutation p (labels permuted among zones within each major, 5,000 draws).
- **Controls and injections:** time-of-day scrambled control cells (`R4.control_cells`) through E4Z: ≥ 90 % of per-major
  contrasts within 2σ and pooled contrasts within 2σ. Model-consistent zone pulses (0.3 % in zone z only, one major at a
  time, cell value a_d share_z × 0.3) → recovery ΔD_z(e)/0.3 (median) and cross-zone leak; round-4 waveform-level Pd
  (all-zone deep) and Ps (shallow) per major → |zone contrast| /0.3 and /0.1. Rule: recovery 0.8–1.1, cross leak ≤ 0.2,
  Ps leak into contrasts ≤ 0.35. Leave-one-site-out of the pooled contrasts as a column.
- **Verdict:** ZONES DIFFER IN AN EPISODE iff a major has p < 0.01 / n_major AND the control passes; ZONES DIFFER POOLED iff
  pooled p < 0.01, a contrast ≥ 3σ, control and injection rules pass; SLIP-LOCKED iff pooled slipped − not ≥ 3σ AND
  permutation p < 0.01 AND its control within 2σ. Else 2σ bounds.

### 11.3 Task 3 — joint arrival-time + coda inversion (idea-4 middle ground)
- **Ray tracing:** 1-D layered region profile (`C.VEL`: ztop_km, vp, vs), source at the family (`s11_kernel_mc.geometry`),
  receiver at the surface; P and S traced separately by ray-parameter bisection; travel time per depth group s 0–5,
  m 5–25, d 25–45, b > 45 km. Check (not a gate): traced tS − tP vs observed SP (`series_geom.csv`), median ratio.
- **Ballistic rows (headline):** `lags_station` bin EP, kind data, not railed (S or P3), finite sigSmP3 > 0;
  y = lagS − lagP3 (ms); per-series (fam, netsta, epoch) offsets demeaned; series need ≥ 3 episodes in the model set;
  kernel (constant Vp/Vs) K_g = −10 (tS_g − tP_g) ms per % dv/v; σ = sigSmP3 × √max(1, χ²_red of the ballistic-only
  episode fit) (one factor per region). Columns: Vs only (K_g = −10 tS_g); + far-pair dS rows (`lags_pair` EP, pairs
  between a co-located reference site and a far station; K_g = −10 (tS_g^B − tS_g^A)); + a source-locus nuisance (one
  constant ms per (family, episode) on S−P rows with ≥ 2 stations).
- **Per-episode model** (episode set = H's episodes): unknowns S(e), M(e), D(e), each sum-to-zero and ⊥ PROJ; coda rows
  d = o + a_s S + a_m M + a_d D + U1 + a_b B(t) (L-SEL coefficients); ballistic rows y = o + Σ_g K_g X_g(e) (K_b × B(t)).
  Fits: coda-only, ballistic-only, joint. Per major ΔX_g(e) ± σ; σ per major (median); composite (mean of majors) ± σ.
- **Decade model:** coda d = o + Σ_g a_g X_g·[t ≥ 2015-01-01] + per-unit seasonal (c1 s1 c2 s2 ≥ 40 cells / c1 s1 ≥ 10);
  ballistic: DH rows (E/L) per series, demeaned, kernel K_g·step; PB.B013 excluded (documented instrument step). Fits
  coda-only, ballistic-only, joint; values ± σ and 2σ bounds.
- **Model-consistent injections** (noise-free, each fit's operator): +0.3 % D pulse at each major → ΔD recovery and leak
  into ΔS, ΔM; +0.1 % S pulse → leak into ΔD; decade +0.1 % deep step → recovery and leaks.
- **Reading (fixed now):** a group is DETECTED iff joint |value| ≥ 3σ AND coda-only and ballistic-only have the same sign;
  else 2σ bounds. A data type "adds" to group g iff σ_joint ≤ 0.8 × σ of the other type alone; reported as σ ratios and the
  ballistic information share 1 − σ²_joint/σ²_coda.

### 11.4 Task 4 — independent test of the WA large-episode hint
- **ΔD source:** round-5 E4 L-SEL (`s11_r6_size.fit_source(H, inj, "E4", "two3_20")`), as round 6. Sizes: round-6
  `r6_sizes.csv`. Measures: duration Du (headline), tremor T (column). **Top set:** the round-6 definition (top
  ceil(n/3) by size among the full-fit episodes with a finite size), fixed on the full WA fit and reused for the halves.
- **Splits (WA):** family halves (fam odd / even); station halves: sites ordered by their total cell weight in the full fit
  and assigned alternately to A / B (co-located instruments share a site, so they stay together); each half refit after the
  ≥ 2-site episode rule. Contrast = round-6 "during" (IVW top − rest, σ from the covariance); scrambled control 5,000
  random top sets (seed 20261012).
- **SVI replication:** the same definition on SVI (Du headline, T column). Caveat stated now: round 6 already showed SVI
  T top third (+0.030 ± 0.029) and the SVI Du slope (−0.009 ± 0.023); the SVI Du top-third contrast was not computed.
- **Stack-length column:** per episode n̄ = weight-averaged events per cell; corr(Du, log n̄) and the WLS of ΔD on
  [Du, log n̄].
- **Verdict (Du):** WA-CONSISTENT iff in both splits each half's contrast > 0 and the halves differ by < 2σ_diff;
  SVI-REPLICATES iff the SVI contrast > 0 with z ≥ 2 and scrambled p < 0.05. **REPLICATED** = both; **WA-ONLY** =
  WA-CONSISTENT only; **NOT ROBUST** = WA not consistent. **STACK-LENGTH** flag iff the Du partial slope drops below 1σ
  with log n̄ included. The halves test internal consistency only; SVI is the only independent replication.

### 11.5 Figures (fixed now; line graphs, colour by weight or cc, no error bars on dv/v)
`fig_r7_front.png` (shape and absolute curves vs lag, both regions; front-scrambled draws as thin grey lines; points
coloured by weight); `fig_r7_zones.png` (per-major ΔD_z lines per zone; slipped − not per major coloured by weight);
`fig_r7_joint.png` (σ per major and decade σ for S/M/D: coda, ballistic, joint); `fig_r7_hint.png` (top-third contrasts:
WA full, halves, SVI; points coloured by weight, scrambled 95 % lines).

### 11.6 Stages (queue `s11_r7_queue.sh`)
Per region: `s11_r7_front.py` → FRONT_<R>_DONE; `s11_r7_zones.py` → ZONES_<R>_DONE; `s11_r7_joint.py` → JOINT_<R>_DONE;
`s11_r7_hint.py` → HINT_<R>_DONE; then `s11_r7_hint.py --stage verdict` → HINT_VERDICT_DONE; `s11_r7_figures.py` →
FIGURES_DONE; QUEUE_DONE. Code tests run first on subsets (tag `_test`, outputs in `round7/<r>/tables/test/`, not results).

**§11 Amendment 1 (Tue 2026-10-13 ~9:38 AM PT; POST-CODE-TEST, labelled so; before the queue run).** Code tests (subsets, tag
_test, outputs in `round7/<r>/tables/test/`, not results) ran every stage end to end.
- (a) Task 3 ray check: in SVI the traced tS − tP is 1.18 × the observed S−P (q10–q90 1.02–1.43; WA 0.98). The
  ballistic bound scales with the kernel, so every ballistic kernel is rescaled by f = observed SP / traced (tS − tP) per
  (family, station) (median over epochs; layer shares kept from the ray; each dS leg by its own f). The unscaled kernel is
  kept as the columns ballistic_raw / joint_raw. Nothing else in §11.3 changes.
- (b) Task 1: the front-scrambled p of the slope b is uninformative by construction (a per-(family, major) lag shift leaves
  a within-episode slope unchanged except through the clipping); it is reported, not read. The verdict uses C (as written).
- (c) Seen before the queue: the WA zone code test (`--test 250`) exceeded the WA series count, so the full WA zone fit
  (both setups) ran in the test: L-SEL pooled equal-zones p 0.47, Σ per-major χ² 109 on 48 dof, Ps leak into contrasts
  0.75 (injection rule fails), control passes. No §11.2 rule was changed after that. SVI front / joint / zone tests ran on
  40–250 series subsets (numbers seen, not results). The hint test ran with the sizes shuffled first.
