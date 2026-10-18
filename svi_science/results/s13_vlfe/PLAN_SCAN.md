# S13 scan — between-ETS VLFE scan ("VLFEs without LFEs"), PILOT: S Vancouver Island, calendar 2011

Written 2026-10-14 (Wed) ~10:30–10:40 pm PT, before any detector, control, injection or calibration number exists.
Design: `MERLIN_VLFE_SCAN.md` (Merlin, 2026-10-14 ~3:30 pm PT), steps 1–7; flags F1–F9. This file fixes every rule.
Any later change is a dated PT amendment written before the numbers it affects. All new files are prefixed `scan_`.

What existed when this was written (no detector numbers):
- prep v3 started 10:23 pm PT (`scripts/scan_prep_v3_run.sh`); its acceptance test on PGC 2011-08-20 was run first (§2).
- truth tables and earthquake-mask event lists (`scripts/scan_truth.py`, §1), including the masked-time fraction per
  year (a property of the earthquake catalogues only). That fraction forced amendment S1 below.

## 0. Question and pilot
Are there discrete 20–50 s slow events (VLFEs) on the S Vancouver Island plate interface BETWEEN ETS, and do they come
with LFEs? Pilot = SVI, calendar 2011 (2011-08 ETS; H&G 2011 templates 2–5; winter + summer). Pilot stations: CN.PGC,
CN.VGZ, CN.SNB, CN.NLLB, CN.PFB, CN.YOUB, UW.LRIV; CN.LZB enters the noise table and is excluded only by the noise /
gain rules, never by hand. Existing waveform data only.

## 1. Truth tables and masks (`scripts/scan_truth.py`)
- `tables/hg2018_table1.csv`: Hutchison & Ghosh (2019, JGR 124, 554–565, accepted 2018) Table 1, 13 templates, with
  a source line and the paper's between-ETS statements (no numeric between-ETS rate is given in the paper; counts are in
  Figs 6–7 only). 2011 SVI-area templates 2–5: Mw 3.7, 3.3, 3.5, 3.4 at depths 54, 52, 30, 55 km. Template 1
  (47.825 N) is WA-side (outside our grid's station coverage but inside the grid; reported, not a pass item).
- `tables/scan_footprint_polygon.csv`: 10-km-buffered convex hull of the SVI families (`families_v1_3.csv`) with ≥ 30
  catalogued P ≥ 0.75 counted-once LFEs 2003–2026 (`scan_footprint_families.csv`: 111 of 115 families).
- `tables/station_epochs.csv`: BH/HH channel epochs from `data/resp/*.xml` (2011: PGC one epoch 2000–2015; SNB BH until
  2011-03-17 18:15 UTC then HH (same CMG3ESP sensor, new digitizer); VGZ HH 2009–2017; NLLB 2007–2017; PFB 2007–2018;
  YOUB 2009–2014; LRIV BH 2009–2019; LZB 2001–2014).
- `tables/scan_eq_events.csv` + `tables/scan_mask_fraction.csv`. Sources: ComCat global M ≥ 4.5 2003–2026 (event
  metadata only, `data/eq/`), `data/rf_teleseis/events_M58_1993_2026.csv`, ComCat M ≥ 2.5 within 300 km and M ≥ 4.0
  within 1,500 km of (48.37 N, −123.61), `lfe_extension/results/trackC/eq_catalog.csv`. Distance Δ from (48.37, −123.61).
  - Regional (all tiers): M ≥ 2.5 within 300 km: origin ± 5 min; M ≥ 4 within 1,500 km: origin ± 20 min.
  - Teleseism (Δ > 1,500 km), three tiers:
    - `m45` (Merlin literal): every M ≥ 4.5: origin → origin + Δ/2.5 km/s + 20 min (M < 6), + 60 min (6–7), + 3 h
      (7–8), + 24 h (≥ 8). Column `in_mask_m45`.
    - `t58`: Merlin literal on M ≥ 5.8 only. Column `in_mask_t58`.
    - `head` (HEADLINE, amendment S1): M ≥ 6, Merlin's full window; M < 6 only if the predicted 20-s amplitude
      A20 ≥ 20 nm, window origin + Δ/4.5 km/s − 300 s → origin + Δ/2.8 km/s + 900 s (Rayleigh/Love group window at
      20–50 s plus coda). A20 from the Prague Ms formula at T = 20 s, Ms = M for Mw/Ms-type magnitudes and
      Ms = (mb − 2.5)/0.63 (Gutenberg–Richter 1956) for mb; Δ clipped to 20–160°.
  - **Amendment S1 (written ~10:35 pm PT, after the mask fractions, before any detector number).** Merlin's literal
    M ≥ 4.5 mask covers 69 % of 2011 and 53–69 % of every year 2003–2026 (`scan_mask_fraction.csv`), not the 3–8 % he
    expected (his range matches the M ≥ 5.8 list: 5.8 % in 2011). A scan on 31 % of the time would make completeness in
    time the dominant loss. The headline tier keeps every large event and every small teleseism whose surface waves are
    predicted to reach the band noise (20 nm ≈ the quietest pilot verticals' 400-s RMS, half a Mw 3.3 VLFE peak) and
    restricts M < 6 windows to the surface-wave group window: 34 % of 2011 masked (2003–2026: 21–34 %; 2011 is the
    Tohoku year). The literal `m45` and `t58` masks stay as columns; separation curves are reported under all three, so
    the effect of the choice is visible. Rate denominators always use the unmasked live time of the tier in use.
- Season: summer = April–September, winter = October–March (UTC dates).
- ETS window (6a): SVI major 2011-08 = 2011-08-12 → 2011-09-03 (`episodes_b.csv`). Between-ETS (step 8) = outside every
  SVI major ±2 d; minor SVI episodes in 2011 (01, 04, 06, 06b1, 09) are a column `in_minor_svi_episode`, not removed.

## 2. Prep v3 (`scripts/scan_prep_v3.py`, root `data/disp1sps_v3/`; v1 and v2 untouched)
- As `s13_prep.py` (staged decimation, sample-rate snapping, atomic writes, S13_MAXGAP 15 s, pre-filt
  0.004/0.008/0.25/0.45 Hz, 600 s NaN at segment edges) with (F3) a fixed 600-s cosine taper per segment end instead of
  `taper(0.05)`, and `remove_response(taper=False)`; (F2) per-day raw gap intervals `gap_iv` and NaN intervals `nan_iv`.
- Acceptance (run 10:23 pm PT): PGC 2011-08-20 v3 vs v2. Band-passed (0.02–0.05 Hz) RMS ratio over 2,700 s → 83,700 s:
  1.0000 on N, E, Z (relative difference 0.0000); raw mid-day relative difference ≤ 0.3 %; |v3|/|v2| over the first 300 s
  of the day 2.4 / 3.0 / 2.6 (N/E/Z) → the F3 ramp, removed. PASS.

## 3. Templates and calibration (`scripts/scan_calib.py`)
- Grid: nodes every 10 km (flat-earth at 48.37 N) over 47.8–49.3 N, −124.8 to −122.6 E; depth = Slab2
  (`cas_slab2_dep_02.24.18.xyz`, nearest grid value) at the node; nodes without Slab2 are dropped. One extra node at the
  exact S13 2011-08 pooled centroid (48.3744, −123.6109, 34.67 km) for 6b.
- Mechanisms: thrust 315/12/90 (headline); strike-slip 315/90/0, steep dip 315/80/90 (competitors); 5-element
  deviatoric basis (free fit). Source time function: 1-s smooth step (headline scan), 20-s (column at candidates).
  Whole-space GF with free-surface factor 2 (`s13_green.synth`), band-passed 0.02–0.05 Hz (4-pole zero-phase), cut to
  the fit window τ ∈ [−60, +120] s about origin.
- Calibration: grand stack per station-component-epoch = n-weighted mean of the per-episode pooled stacks (SVI majors,
  file preference `scalef_` > `refit_` > `scale_` as in `s13_scale_summary.py`), using only episodes inside the
  station's 2011 channel epoch and only station-components not dropped by guard2 in that episode. σ of the grand stack
  per component = sqrt(Σ n_e² σ_e²)/Σ n_e (σ_e = control SD from `*_stacomp.csv`). The thrust synthetic at each
  episode's pooled centroid with that episode's fitted lag is n-weighted the same way. M0_grand = joint weighted fit
  (1/σ²) of the grand stacks onto the grand synthetic (all kept components). Gain g_sc = ⟨T_sc, S_sc⟩ / (M0_grand ⟨T_sc,
  T_sc⟩); g_snr = ⟨T,S⟩/(σ_sc |T_sc|). Rule (Merlin): keep a station-component-epoch if g > 0 and |ln g| ≤ ln 2, with
  template scaled by g; otherwise exclude it from the scan. A common factor (e.g. the ×2 free surface) is absorbed in
  M0_grand and cannot be tested by this calibration; relative gains are the product.
- Acceptance: g in [0.5, 2] and positive for ≥ 12 of 21 pilot station-components (7 stations × 3); Z gains at PGC, VGZ,
  LRIV in [0.7, 1.4]. Failure: Z gains all ≈ 0.5 or ≈ 2 → carry the factor into the moment scale, continue.
- Noise table `tables/scan_noise_month.csv`: per station-component-month, median band RMS over non-overlapping 400-s
  windows (< 5 % NaN). Exclude a station-month from the scan if its Z median > 3 × that station's median over all its
  v2 + v3 days on disk (the available 2003–2026 days) or > 5 × the network median of Z medians that month.
- Empirical-template option: if any 2011 ETS-time detection or H&G event has Mw ≥ 3.8 with ≥ 500 nm at PGC/VGZ/LRIV
  and SNR ≥ 30 per component, cut it and run it as a second detector. H&G 2011 templates are Mw ≤ 3.7, so this fires
  only on a scan detection.
- Output `tables/sc_gain.csv`, `data/scan/templates_grid.npz`.

## 4. Detector (`scripts/scan_detect.py`, cupy FFT)
- Data: v3, per station, per component; per contiguous finite run: linear detrend, 120-s cosine taper into each edge,
  0.02–0.05 Hz 4-pole zero-phase band-pass, then NaN within 300 s of every NaN sample (mask dilation 300 s). Runs
  < 1,800 s dropped. Instrument vetoes (`svi_science/instrument_vetoes.csv`) blank their component-days first.
- Per UTC day, per station-component: σ_sc = 1.4826 × MAD of the band-passed component over the day (finite, outside
  the headline earthquake mask). Effective template T'_sc = g_sc T_sc; weight w_sc = 1/σ_sc²; excluded components w = 0.
- M̂0(t, node) = Σ_s a_s(t) Σ_c w_sc (T'_sc ⋆ d_sc)(t) / Σ_s a_s(t) Σ_c w_sc |T'_sc|², where a_s(t) = 1 if station s
  is finite over the whole fit window [t − 60, t + 120] s, else 0 (prep NaNs are common to the 3 components).
- z(t, node) = M̂0 / (1.4826 × MAD of M̂0(·, node) over the UTC day, finite and outside the headline mask; MAD on every
  4th sample; whole day if < 10 % is unmasked). Network CC(t) = Σ w T' d / sqrt(Σ w T'² · Σ w d²) over the fit window.
- Z-only variant zZ (verticals only, same rules) computed in the same pass.
- Candidates: per second, zmax = max over nodes (argmax = best node); local maxima of zmax ≥ 4 with M̂0 > 0, declustered
  ±150 s keeping the largest; negative exceedances (min over nodes of z ≤ −4, same declustering) kept as the sign control.
- Columns per candidate (never cuts): time (UTC and PT), best node lat/lon/depth, z, M̂0, Mw = (log10 M̂0 − 9.1)/1.5,
  CC, zZ at the best node, n_sta (stations with a_s = 1 and weight > 0), n_Z_agree (verticals with own ⟨T'_Z, d_Z⟩ > 0),
  max_share (largest single station-component share of Σ w T' d) and its channel, VR_thrust, VR_strikeslip,
  VR_steepdip, VR_free (5-element deviatoric, same lag), M̂0 with the 20-s template (best of lags −10…+10 s),
  isolation, near_gap_s (distance to the nearest raw gap at any contributing station), in_mask_head / m45 / t58 /
  regional, HF ratio at the two nearest pilot stations, class_eq (HF ratio ≥ 10 at both), minutes_to_midnight, season,
  in_ets (2011-08-12…09-03), in_svi_major_pm2d, in_minor_svi_episode, Δz to the second-best node region (Δz ≤ 1 region
  size in km).
  - isolation = max over the windows [−900, −300] s and [+300, +900] s of the template-weight-share-weighted mean over
    station-components of (band RMS in the window / that station-component's daily median 400-s RMS).
  - HF ratio = peak of the 1-s 2–8 Hz RMS envelope (Z, raw counts, `scripts/scan_hfenv.py`) within ±60 s of the
    candidate / that station-day's median envelope.
- Headline filter (a detection): z ≥ threshold, M̂0 > 0, not in the headline earthquake mask (head ∪ regional),
  n_sta ≥ 3, n_Z_agree ≥ 2, max_share ≤ 0.5, VR_thrust / VR_free ≥ 0.8, VR_thrust > VR_strikeslip and
  > VR_steepdip, isolation < 2, not class_eq. Z-only check: headline detections with z ≥ 1.2 × threshold must also have
  zZ ≥ threshold (agreement fraction reported; rule: ≥ 80 %).

## 5. Controls (identical detector, identical columns and headline filter, thresholds 4, 4.5, 5, 5.5, 6, 7)
- (a) Station-day-shifted scrambled control: station k (pilot order PGC, VGZ, SNB, NLLB, PFB, YOUB, LRIV, k = 0…6)
  uses its own data from day (d + s·k) mod 365 of 2011 (same epoch by construction); five realisations
  s ∈ {3, 5, 7, 11, 13}; output timeline and masks = the unshifted day d. Rate = mean over the five.
- (b) Polarity-scrambled template control: unshifted data, random sign per station-component (50 draws, seed 20261014),
  gains unchanged.
- (c) Sign control: negative-M̂0 exceedances in the real data.
- Time reversal is NOT used (Merlin).
- Rates per year = count / (unmasked live time in years); live time = seconds with ≥ 3 stations available outside the
  mask tier in use.
- Headline threshold = smallest listed threshold at which (a) ≤ 1/yr AND (b) ≤ 2/yr after masks and headline filter.
  Separation curve = real positive, real negative, (a), (b) per threshold; reported for raw exceedances and for headline
  detections, under each mask tier.
- Acceptance: at the headline threshold, ETS-window positives ≥ 10 × ETS-window negatives (headline filter applied to both
  signs, negative with the sign-flipped template for the VR rules). Failure: positive ≈ negative in ETS → check
  polarity vs the grand stack; stop.

## 6. Validation inside 2011 (`scripts/scan_validate.py`)
- (a) ETS re-find: headline detections per unmasked live day in 2011-08-12…09-03 vs the rest of 2011: ratio ≥ 10.
  Also Σ M̂0 of ETS-time detections / 7.7e15 N m (S13 time-locked total, 5.3e12 × 1,445) = fraction of the LFE-locked
  VLF moment in discrete detectable VLFEs (a new number).
- (b) Mechanical consistency: mean M̂0(t) at the 2011-08 pooled-centroid node, sampled at the 1,445 P ≥ 0.75
  counted-once LFE origins of 2011-08 (S13 origins, + S13's fitted lag 1 s), dropping an LFE if its 400-s window fails
  S13 window QC (RMS > 3 × the episode median of that station-component's windows, or > 5 % masked) at any
  station-component carrying ≥ 10 % of the template weight. Must equal 5.3e12 within 20 % (rule on the headline gains;
  g = 1 version as a column). Failure > 2× → implementation bug; stop.
- (c) H&G templates 2–5 (2011-08-21 04:09:06, 08:43:29; 08-24 10:47:38; 08-26 13:44:34 UTC) must be candidates with
  z ≥ the headline threshold within ±60 s, best node within 40 km of the H&G epicentre (masks and other filters
  reported, not applied for this test). Expect ≥ 3/4; compare our Mw with Table 1. Failure ≤ 1/4 while (d) recovers
  ≥ 80 % at Mw ≥ 3.5 → truth table or polarity/timing error; debug before any between-ETS number.
- (d) Injection (`scripts/scan_inject.py`): synthetic thrusts (g_sc × whole-space GF, no noise added) added to the v3
  displacement before band-passing: Mw 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0 × 40 each (20 with 1-s, 20 with 20-s source
  duration), plus 40 at Mw 3.6 with a perturbed GF (free surface 1.0, Vs 3.4 km/s, Vp 6.0 km/s, 30-s duration). Times:
  uniform random over 2011, outside the headline mask, with ≥ 3 stations live, ≥ 2 h apart and ≥ 2 h from any H&G
  template time. Locations: half at random footprint-family locations (family depth), half at grid nodes whose distance
  to the nearest family is 30–50 km (Slab2 depth). One injected pass of the same detector; a recovery = headline
  detection within ±60 s. Outputs: recovery vs Mw by season and duration (completeness Mw at 50 % and 80 %, linear
  interpolation in Mw), M̂0 bias (median log10 ratio), location error (best node vs truth, km): median and 95 % radius
  vs Mw.
- (e) ETS-stack location test: grid search of the VR-maximising node (thrust, ±8 s lag scan, S13 stored weights) for the
  2012-08 pooled stack (`refit_svi_2012-08`) and 2018-06 (`scalef_svi_2018-06`). Target = the stack's own count-weighted
  centroid (2012-08: 48.440, −123.752; 2018-06: 48.447, −123.797); distance to Merlin's (48.37, −123.61) as a column.
  Acceptance: within 20 km. Failure > 30 km or a ridge along station geometry → Okada half-space location templates.

## 7. Pilot go/no-go (Merlin step 7, verbatim)
- GO if all: 6a ≥ 10×; 6b within 20 %; 6c ≥ 3/4; 6d ≥ 80 % at Mw 3.6 and ≥ 50 % at 3.3 (summer), location median
  ≤ 30 km at 3.6; controls ≤ 1/yr (a) and ≤ 2/yr (b); 6e within 20 km.
- The between-ETS count is NOT a criterion: zero is a result (rate bound at measured completeness) if the above hold.
- NO-GO signatures → adaptations: (i) H&G not re-found though injection says they should be → truth-table or polarity
  error, debug; (ii) control (b) > 1/month at the Mw-3.6 threshold → coherent long-period noise: require HF-quiet +
  isolation in the headline, raise the low corner to 25 s, add era stations (C8 POLARIS 2003–06, XU CAFE 2006–10, UW
  compact 2021+ in `inventory_bb_near.csv`); if still failing, restrict to Mw ≥ 3.8 and fall back to the
  stack-at-candidate-times path (S13 A5 extended to P 0.5–0.75 candidates and to PNSN tremor times between ETS) —
  answers the complementary question at coarse resolution; (iii) completeness worse than Mw 3.8 in every season → same
  fallback.
- Interpolation for "≥ 50 % at 3.3": recovery at Mw 3.3 = linear interpolation between the 3.2 and 3.4 summer bins.
  "Location median ≤ 30 km at 3.6" uses the recovered Mw-3.6 standard-GF injections (both durations).

## 8. Step 8 (only if GO, on 2011 between-ETS headline detections) — as `MERLIN_VLFE_SCAN.md` step 8, unchanged.

## Outputs
`tables/scan_*`, `data/scan/`, `figures/scan_*` (raw first: z(t) for 2011 with masks and ETS; 2011-08 panel with H&G
ticks; separation curves; points coloured by z/CC, no error bars), RESULTS.md "Scan pilot (2011)", NOTES.md entries.

## Amendment S2 — calibration weights (written 2026-10-14 ~10:43 pm PT, after the first gain run, before any detector number)
- First gain run (rule as in §3, n-weighted, guard2 drops respected) is kept as `tables/scan_sc_gain_run1_nweighted.csv`:
  7 of 21 kept; Z gains PGC 0.91, VGZ 1.35, LRIV 19.6, SNB 4.9 → acceptance FAILED. Two mechanical causes, both visible in
  the inputs, not in any detector output:
  1. n-weighting lets broken-sensor episodes dominate: UW.LRIV BHZ control σ is 89–411 nm in 2015-12 … 2018-06 (vs
     0.3–2 nm in other episodes) inside the same 2009–2019 channel epoch, so the n-weighted LRIV grand stack is noise.
  2. guard2's `<0.05xsta` clause drops healthy verticals when the horizontals are tilt-noisy (already noted in NOTES.md
     2026-10-14 ~7:15 pm PT and `tables/guard2_zvariant.csv`): SNB Z is dropped in 19 episodes, leaving 1 for the HH epoch.
- Rule change: per-episode weights 1/σ_e² (inverse variance; same weights for the stacks and the synthetic),
  σ_grand = (Σ 1/σ_e²)^−1/2; for calibration only, guard2 drops whose sole reason is `<0.05xsta` are ignored
  (drops for n_used < 5, σ ≤ 0, `<0.02xunit` still respected). Everything else in §3 unchanged (same keep rule,
  same acceptance). The scan itself (§4) never uses guard2.

## Amendment S3 — z normalisation (written 2026-10-14 ~10:46 pm PT, after a 4-day code test, before the year run)
- Code test `data/scan/runs/test_d232` (days 2011-08-20, 21, 24, 26; seen: candidate lists incl. the H&G-time values;
  not used for any rule or result; all numbers are recomputed by the year run). It exposed a pathology of §4's
  z = M̂0 / MAD(M̂0): when station availability changes within a day, the denominator Σ a w g²|T|² changes, and so does
  the noise of M̂0 (var M̂0 ∝ 1/den for 1/σ² weights). At the far edge node 170 (48.70 N, −124.80) z reached 1e6 when the
  near stations dropped out. Fix: z(t, node) = s(t)/(1.4826 MAD of s over the UTC day, unmasked) with
  s = M̂0·sqrt(den) = num/sqrt(den) (the normalised matched-filter statistic). Identical to §4 when availability is
  constant through the day. M̂0 = num/den is unchanged and still reported. The Z-only zZ gets its own normalisation the
  same way.

## Amendment S4 — glitch mask (written 2026-10-14 ~10:49 pm PT, after the 4-day code test, before the year run)
- The code test showed instrument glitches that S13 never saw (S13 averages random windows and its window QC drops
  outliers): CN.PGC E has a ~10 mm excursion within seconds of every UTC midnight (raw v3 values, present in both
  adjacent independently processed day files); CN.YOUB Z has bursts of 2e7 nm (band-passed); UW.LRIV N/E and CN.SNB
  have many. Counts of band-passed samples > 200 × σ_day outside the earthquake mask in 2011: PGC E 71,494 (328 days),
  YOUB Z 286,803 (47 days), LRIV N 29,781 (52 days), SNB E 31,021 (57 days); 1–2 orders fewer elsewhere.
- Rule: per station-component-day, every band-passed sample with |x| > 200 × σ_day (σ_day = 1.4826 × MAD, §4) is a
  glitch; the component is set to NaN within ±600 s of every glitch sample (covers the band-pass ringing). Applied in
  `scan_bp.py` (new cache `data/scan/bp_2011g/`; `bp_2011/` kept) before the detector, the controls and the injections,
  identically. Why 200: a Mw 4.3 VLFE (largest H&G size) gives ≤ 900 nm on the nearest vertical (≈ 0.5 nm per 1e12 N m)
  = 45–90 σ at σ 10–20 nm, so no real VLFE in the H&G range reaches 200 σ; the glitches sit at 10^3–10^6 σ. Masked time
  per station-component is reported (`tables/scan_glitch_mask.csv`).

## Amendment S5 — sign side lobes, polarity draws, filter validation (written 2026-10-14 ~10:55 pm PT, after the 4-day code test `test4_d232` (glitch-masked), before the year run)
What the code test showed (code check only; numbers are recomputed by the year run):
- The band-passed step template is antisymmetric in time (a zero-phase band-pass of a step is odd about the step), so
  every event has an opposite-sign match 13–17 s away with 60–110 % of the main |z|: at the H&G 2011-08-21 04:09 time
  +18.0 at 04:09:09 and −15.6 at 04:08:55; at 2011-08-24 10:47 +18.8 at 10:47:38 and −20.8 at 10:47:53. Counting
  negative peaks separately (§4–5) would count the side lobe of every real event as "sign control".
- The polarity-scrambled template control (b) keeps part of every real signal: with ~12 live station-components a
  random sign draw typically shares a large fraction of the true template's weighted energy, so real VLFEs (and the
  ETS) leak into (b).
- Headline filter values at the four H&G times: max_share 0.58–0.64 (all fail ≤ 0.5: PGC Z carries most of the fit
  weight on this network, as in S13 where PGC Z held 46–78 %), VR ratio 0.96–0.97, isolation 0.7–1.4; at 04:09 the
  2–8 Hz envelope peaks 9–24 s after the VLFE origin at all five stations (PGC 15×, LRIV 33× the day median; no
  catalogued earthquake within ±10 min) → class_eq.
Rule changes (all mechanical, fixed before the year run):
1. Joint-sign clusters: candidates of both signs are merged and declustered together (±150 s, keep the largest |z|);
   the cluster's sign is the sign of that peak. Rates, controls (a)/(b), separation curves, 6a and the ETS sign
   acceptance (positives ≥ 10 × negatives) use clusters. The per-sign lists stay (6c uses the positive list, as written).
2. Control (b) draws: generated sequentially from seed 20261014; a draw is used only if the weighted correlation of the
   scrambled with the true thrust template, ρ = Σ s_sc w̄_sc g_sc² |T_sc|² / Σ w̄_sc g_sc² |T_sc|² (S13 centroid node,
   w̄ = median 2011 daily weight), satisfies |ρ| ≤ 0.3; the first 50 that qualify. (b) rates are reported for all of 2011
   and for between-ETS time (outside SVI majors ±2 d); the headline-threshold rule uses the between-ETS rate (the scan's
   target), the all-year rate is a column.
3. Filter validation (project rule: a cut must keep ≥ 85 % of known events): each headline criterion is checked on the
   recovered standard-GF injections with Mw ≥ 3.4 (recovered = cluster with z ≥ 4 within ±60 s). A criterion keeping
   < 85 % of them is replaced by its pre-declared fallback if that keeps ≥ 85 %, otherwise dropped. Fallbacks:
   max_share ≤ 0.5 → share_excess = max over station-components of (share of Σ w T' d − expected share
   w T'² / Σ w T'²) ≤ 0.25. No fallback for n_sta ≥ 3, n_Z_agree ≥ 2, VR ratio ≥ 0.8, VR thrust > competitors,
   isolation < 2. class_eq cannot be checked on injections (the HF envelope is the real data's); it stays as Merlin
   wrote it, and HF-loud clusters (HF ratio ≥ 10 at both nearest stations) are reported as their own class.
   The final headline filter is then applied unchanged to real data and both controls.

## Amendment S6 — source-duration convention bug (written 2026-10-14 ~11:40 pm PT, after the first validation pass)
- The first validation pass (`tables/scan_verdict.json` as of 11:40 pm PT, kept as `tables/scan_verdict_pass1.json`)
  recovered 0 of 140 standard-GF injections with the "20-s" source and 1 of 40 perturbed-GF ("30-s") injections, while
  the 1-s injections at Mw ≥ 3.8 gave z 20–105. Cause: the code passed τ = D to the erf moment function
  M(t) = ½(1 + erf(t/τ)), whose moment rate has spectrum exp(−(π f τ)²): at 0.03 Hz that is 0.03 for τ = 20 s and
  4e-4 for τ = 30 s. A τ of 20 s is a ~50-s rise (10–90 %), not a 20-s source.
- Fix (implementation only, the stated design is unchanged): duration D = 10–90 % rise time, τ = D / 2.33. 20-s
  source: τ = 8.6 s; 30-s (perturbed set): τ = 12.9 s. The headline template keeps τ = 1 s (S13's "1-s" step, 2.3-s
  rise; at 20–50 s indistinguishable from shorter). Same injection times, locations and Mw (same seed); the old list
  and run are kept as `data/scan/inject_2011_tauD.parquet` and `data/scan/runs/inj/`. The 20-s template column
  (`m0_20s`) in the real-run table still uses τ = 20 s; it enters no rule and is flagged as such.

## Amendment S7 — second pilot, SVI 2011 (written 2026-10-15 (Thu) ~0:15–0:50 am PT, before any second-pilot number)
Source: Merlin's NO-GO review `MERLIN_SCAN_NOGO_REVIEW.md` (steps 1–7). When this was written no reference station had been
prepped for the scan, no decay stack, veto column, cleaner run, new mask fraction, control rate or classification
existed. Reused unchanged: v3 prep of the pilot stations, gains (`sc_gain.csv`), template grid, noise table, S3 z, S4
glitch mask, S5 filter rules, S6 injection list, H&G truth, footprint polygon. Pilot-1 runs and tables are kept as they
are; every new file is prefixed `scan2_` (scripts: `scan_refs.py`, `scan_prep_ref.py`, `scan_refqc.py`, `scan_decay.py`,
`scan_masks2.py`, `scan_bwlist.py`, `scan_detect2.py`, `scan_columns2.py`, `scan_validate2.py`, `scan_pilot2_queue.sh`).

### S7.1 Reference stations (Merlin step 1)
- On disk and used now: UW.GNW (all 2011), UW.DOSE (2011, column only), and the ETS-days set, 2011-07-30 … 08-27 (DOY
  211–239): CC.OBSR, UW.LEBA, UW.RATT, UW.ELW, TA.D04D, UW.WISH, UW.SP2 (E vetoed from 2011-09-05 by
  `instrument_vetoes.csv`, outside these days), TA.B05D, US.NLWA, TA.D03D.
- Asked of the user, NOT approved, NOT fetched: 1-sps LHZ 2011 for IU.COR, US.NEW, US.HAWA, CN.PNT, CN.BBB (~300–540 km).
  Code is written so they drop in: `scan_prep_ref.py` accepts LH? channels (no decimation; responses from
  `data/resp_scan/{NET}.xml`), and every later script uses whatever reference band-passed files exist. Re-running
  prep → bp → columns → validate adds them. Until then the year-round veto rests on GNW alone, which is weak above
  Mw 3.8 (a local event that large puts a non-negligible amplitude at 100–150 km, so a local event and a wave from
  outside look alike at GNW); every result states this.
- Coordinates: `lfe_extension_wa/results/wa_temp_networks/stations.csv` where listed (it holds the TA stations only),
  otherwise `lfe_array_discovery/results/stationxml/{NET}.xml`. Responses: `data/resp/{NET}_bb.xml` if it has the
  station, else `lfe_array_discovery/results/stationxml/{NET}.xml`, else `data/resp_scan/{NET}.xml`. Table
  `tables/scan2_ref_stations.csv`.
- Role: veto reference if ≥ 100 km from (48.37 N, −123.61): GNW 107, D03D 101, NLWA 110, B05D 112, SP2 135, WISH 140,
  D04D 146, ELW 161, RATT 170, LEBA 203, OBSR 210 km. DOSE (87 km) = columns only, never a rule.
- Prep: the `scan_prep_v3.py` recipe (wrapper only changes the response lookup and allows LH). Band-pass, daily σ and
  S4 glitch mask: `scan_bp.py` unchanged.
- Reference QC and the GNW quality check (`scan_refqc.py` → `tables/scan2_refqc.csv`), fixed now:
  - per station-month, Z median 400-s band RMS; a month is excluded if > 3 × that station's 2011 median (same rule as
    the pilot noise table); glitch-mask seconds and finite fraction per month reported;
  - amplitude check against CN.PGC on 2011 events M ≥ 6.8 (Δ ≥ 1,500 km) with both stations live: ratio of the peak
    |Z| in the band (0.02–0.05 Hz) inside the event's head-mask surface-wave window. A reference is used only if the
    median ratio is within [0.5, 2]; no event in its days → used, flagged "unchecked". Failing → excluded.

### S7.2 Decay measurement (Merlin step 2)
- S13 stacker `s13_stack2_v2.py` unchanged, pooled unit only (row file filtered to unit == 'pooled'), `--kmax 1`.
  - `scan2dec_svi_2011-08`: S13_DISP = v3; stations PGC, VGZ, SNB, NLLB, PFB, YOUB, LRIV + GNW, DOSE, WISH, SP2, D04D,
    B05D, D03D, NLWA. OBSR, LEBA, RATT and ELW are not in `data/resp`, whose coordinate table the stacker needs, so
    they are not stacked.
  - `scan2dec_svi_2012-08`: S13_DISP = v2; the `refit_svi_2012-08` station list + GNW, DOSE, WISH, SP2, D03D, D04D, NLWA.
- Per station, Z: A = sqrt(max(RMS² − σ², 0)), RMS of the pooled Z stack over the stacker's fit window at its fitted
  lag, σ = that component's control SD (`*_stacomp.csv`). Decay ratio R_obs = A / A_PGC. If A ≤ 2σ the station is "not
  detected" and R_up = RMS / A_PGC (upper bound) is reported. Predicted R_pred = RMS(template Z) / RMS(PGC template Z).
- Rule: PASS if R_obs (or R_up) ≤ 0.5 at GNW in both episodes. If not, the veto's predicted reference amplitude is
  multiplied by c = max over the two episodes of R_obs/R_pred at GNW; if R_obs > 0.8 at GNW in either episode, GNW
  cannot tell local from travelling: the veto then uses only the ETS-days set (no year-round veto), stated.

### S7.3 Veto columns and their validation (Merlin step 3)
Computed (`scan_columns2.py`) for every positive cluster with z ≥ 5 in every run, and at the validation times:
- A_near = sqrt(Σ_sc shw_sc · mean(d_sc²)) over the near-station components active in the fit window [t − 60, t + 120] s;
  shw = template-weight share (S5).
- Per live veto reference r (QC-passed month, finite data): A_r = largest RMS (over r's finite components, 3-comp) of a
  181-s window starting within ± L_r of t − 60, L_r = max(60 s, Δ_r / 2.5 km/s), Δ_r = distance best node → r (allows
  any arrival-time difference of a travelling wave); P_r = RMS of M̂0 · T_r (whole-space thrust synthetic at r for the
  best node, g = 1, band-passed, same components, × c of S7.2 if it applies); N_r = RMS over components of σ_day.
  - ref_excess_r = (A_r − sqrt(P_r² + N_r²)) / N_r; ref_ratio_r = A_r / A_near.
  - Summary: n_ref, ref_excess = median over live veto references, ref_ratio = median. DOSE: dose_excess, dose_ratio
    (columns only).
  - Controls: in (a) references are shifted like the near stations (station index k = 7, 8, … in the reference-table
    order, day (d + s·k) mod 365); in (b) references are unscrambled (they enter only as amplitudes); in the injection
    runs the synthetic is injected into the reference data too (g = 1).
- Plane-wave fit (column only, never a rule): Z cross-correlation (lags ± 150 s) of every live station (near +
  references) with the near station of largest Z RMS, window [t − 60, t + 120] s; stations with CC ≥ 0.6 → least-squares
  plane wave t_i = t0 + s·x_i → pw_n, pw_vr, pw_slow (s/km), pw_baz; NaN if < 4 stations.
- fp_interior = best node inside the footprint polygon; grid_edge = best node on the outer row or column of the grid.
- Veto variants, used in this order, the first that passes validation is the rule: V1 ref_excess ≥ 3 AND ref_ratio ≥ 0.5;
  V2 ref_excess ≥ 3; V3 ref_excess ≥ 2 AND ref_ratio ≥ 0.3 (all need n_ref ≥ 1; n_ref = 0 → not flagged).
- Validation (all three): (i) H&G 2 and 4 not flagged (positive peaks within ± 60 s of the H&G times); (ii) ≥ 80 % of the
  known body-wave list (S7.4) flagged — 9 of 11 if the list has 11 members; (iii) ≥ 85 % of the recovered standard-GF
  family-location injections with Mw ≥ 3.6 not flagged (project 85 % rule). None passes → the veto stays a column, the
  verdict says "veto unvalidated", and the pilot is NO-GO.

### S7.4 Known body-wave list (reconstruction of Merlin's 11; `scan_bwlist.py` → `tables/scan2_bw_known_list.csv`)
- Population: the 47 pilot-1 survivors (real, variant 0, headline mask, S5 filter, z ≥ 7, z_loso ≥ 7;
  `data/scan/runs/real/loso.parquet`).
- A survivor is a known body-wave cluster if its time lies in [t_P − 60 s, t_S + 600 s] of a ComCat event with M ≥ 5.0 at
  Δ ≥ 1,500 km; t_P, t_S = first P-type and first S-type arrival (iasp91, `obspy.taup`, phases P/Pdiff/PKP*, S/Sdiff/SKS*)
  at (48.37 N, −123.61) for the event depth. The list from this rule is used whatever its size (the 80 % rule applies to
  it) and its difference from Merlin's 11 is reported.

### S7.5 Single-station step cleaner (Merlin step 4; `scan_detect2.py`)
- Pass 1 per run and UTC day: z_s(t, n) = Σ_c w g (T ⋆ d) / sqrt(Σ_c w g² |T|²) for one station s at node n (the S3
  statistic restricted to one station), normalised by its own daily MAD (unmasked, every 4th sample). For each station and
  second, n* = argmax over nodes of |z_s|; flag if |z_s(t, n*)| ≥ 8 while every other available station has
  |z_s'(t, n*)| < 2.5, with ≥ 2 other stations available. Flags are dilated ± 600 s and the station's three components set
  to NaN there, on the run's timeline (so control (a) and the injections are cleaned the same way). All variants of a run
  use the cleaned data. Intervals saved per run (`clean_iv.parquet`).
- Pass 2 = the pilot-1 detector (§4, S3, S4) on the cleaned data; also z maximised over footprint-interior nodes only
  (column). Runs: `real_c` (+ 50 polarity draws, same seed and |ρ| rule), `ctrlA_s{3,5,7,11,13}_c`, `inj2_c`, `inj3_c`.
- Clusters (pilot 2): positive peaks only (z ≥ 4), declustered ± 150 s keeping the largest z. S5.1 joint-sign clusters are
  dropped: the sign test is no longer a criterion (Merlin) and the antisymmetric side lobe made the joint sign a coin flip
  (H&G 4 came out negative). Negative peaks are still saved.
- Acceptance: control (a) ≤ 20/yr at z ≥ 7 (pilot-1 headline mask and S5 filter, no z_loso) AND recovery of standard-GF
  Mw ≥ 3.6 injections (z ≥ 7, ± 60 s, S5 filter) ≥ 85 % of the pilot-1 recovery. Injections failing → the cleaner is not
  used (pilot 2 falls back to the pilot-1 runs); (a) failing alone → the cleaner is used and the miss stated.

### S7.6 Masks as columns (Merlin step 5; `scan_masks2.py` → `tables/scan2_eq_events.csv`, `tables/scan2_mask_fraction.csv`)
- head2 = S1 head with Ms = mb − 0.3 for mb-type magnitudes (replaces (mb − 2.5)/0.63).
- bw5 = [t_P − 60 s, t_S + 600 s] of every M ≥ 5.0 event at Δ ≥ 1,500 km (the S7.4 window).
- reg3 = ComCat M ≥ 3.0 within 2,500 km of (48.37, −123.61), 2011 (event metadata only, fetched for this amendment):
  [origin, origin + Δ / 2.5 km/s + 600 s].
- Headline mask v2 = head2 ∪ bw5 ∪ regional (S1) ∪ reg3, provided reg3 alone covers ≤ 15 % of 2011; otherwise reg3 is a
  column only (stated). head, m45, t58 stay columns. Rates use live time outside the mask in use; separation curves are
  reported under v1 and v2.

### S7.7 Headline detection (pilot 2)
Positive cluster; not in mask v2; S5 filter as validated (re-checked on the pilot-2 injections, Mw ≥ 3.4, 85 % rule,
same fallback) WITHOUT class_eq (HF-loud is a class, never a cut — H&G 2 is HF-loud); fp_interior; the validated veto
does not fire; z ≥ z_thr and z_loso ≥ zl_thr. z_loso as `scan_loso.py` (drop the station with the largest Σ_c w T' d
at the best node; same daily MAD scale).

### S7.8 Injections and the threshold pair (Merlin step 6)
- Extra set `inj3` (`scan_inject.py` logic, seed 20261016): 40 at Mw 3.8 and 40 at Mw 3.6, standard GF, footprint-family
  locations, half 1-s / half 20-s; same time rules (outside head, ≥ 3 live, ≥ 2 h from each other, from H&G times and from
  `inj2` times). Completeness at Mw 3.8 = recovered fraction of all family-location standard-GF Mw 3.8 injections
  (inj2 20 + inj3 40) outside mask v2. Off-family injections are reported (they fail fp_interior by design).
- Grid z_thr ∈ {5, 6, 7, 8, 10, 12, 15, 20} × zl_thr ∈ {3, 4, 5, 6, 7, 8, 10}. Chosen pair = maximum Mw-3.8 recovery with
  (a) ≤ 0.3/yr and (b) ≤ 0.7/yr between ETS; ties → smaller (a) + (b), then larger z_thr. (a) = mean of 5 shifts per
  live year (all year, mask v2); (b) = mean of 50 draws per between-ETS live year. No qualifying pair → NO-GO on
  controls; the pair with the smallest (a) + (b) among those with ≥ 50 % Mw-3.8 recovery is reported.

### S7.9 GO / NO-GO (Merlin step 7; replaces §7 for pilot 2)
GO if all: (a) + (b) ≤ 1/yr between ETS at the chosen pair (strong GO ≤ 0.5); ≥ 80 % recovery at Mw 3.8 with median
location error ≤ 15 km; H&G 2 AND 4 are headline detections at the chosen pair; ETS-window detections ≥ 3 with Poisson
p < 0.01 against the between-ETS real rate (p = P(N ≥ n_ETS | λ = r_between × ETS live time), r_between from the real
headline detections between ETS, floored at 0.5 / between live time if zero); the veto validated (S7.3). The between-ETS
count is not a criterion. NO-GO → Merlin's fallback (stacks at candidate times) is the coordinator's next decision; not
run here.

### S7.10 Individual classification: 2011-08-18 21:19:46 UTC (14:19:46 PDT; pilot 1: z 42.9, Mw 3.86, best node
49.24 N −123.18 on the north grid edge, LRIV-dominant, z_loso 30.2) and the 2011-08-22 triple (08:21:49, 08:25:25,
08:28:22 UTC = 01:21–01:28 PDT). Each, first match wins:
1. instrument — removed by the step cleaner, or z_loso < 5 in the cleaned run;
2. catalogued earthquake — inside head2, bw5, regional or reg3 (event named); m45 membership reported;
3. travelling wave — the validated veto fires with every live reference of that day (ETS-days set included); the plane-wave
   fit and DOSE are reported as support, not as the rule;
4. local candidate — none of the above, n_sta ≥ 3, ref_excess < 3 → "candidate VLFE". For 08-18 the statistic is then
   recomputed on an extended grid to 49.6 N (10-km Slab2 nodes, same templates, gains and weights; `scan_ext_grid.py`)
   at that time: best node, z, and whether the node is interior; still on the edge → "source north of the grid,
   unresolved";
5. otherwise "unresolved".

### S7.11 Outputs and cost
`tables/scan2_*`, `figures/scan2_*` (raw first: z(t) with masks, cleaner intervals and references; separation curves;
veto validation; points coloured by z or CC, no error bars), `tables/scan2_verdict.json`, RESULTS.md "Scan pilot 2
(2011)". Compute ~1–1.5 h (GPU detector ~25 min, columns ~30 min), ≤ 8 CPU workers, < 100 GB.

### S7.12 Clarification (written 2026-10-15 ~1:10 am PT, after the cleaner acceptance numbers, before any fallback number)
- The step cleaner FAILED its injection acceptance (S7.5): standard-GF Mw ≥ 3.6 recovery 0.23 vs 0.82 in pilot 1 (keep
  29 %, rule ≥ 85 %); control (a) at z ≥ 7 was 1,061/yr (rule ≤ 20). As S7.5 prescribes, the cleaner is NOT used: pilot 2
  falls back to the uncleaned pilot-1 detector runs (`real`, `ctrlA_s*`, `inj2`), and `inj3` gets an uncleaned run
  (`inj3_u`, `scan_detect2.py --noclean`, otherwise identical). All S7.3–S7.10 rules then apply to these runs unchanged
  (columns from `scan_columns2.py`, with an empty cleaner interval list). Single-station events are handled by z_loso
  (S7.7/S7.8), which was Merlin's other route.
- fp_interior is a cut and so falls under the project 85 % rule (S5.3) like every other cut: it is part of the headline
  only if it keeps ≥ 85 % of the recovered standard-GF family-location injections with Mw ≥ 3.6 in the runs used;
  otherwise it stays a column. (Implied by S5.3 but not written in S7.7; the cleaned runs gave 0.66, which is not the
  fallback value.)

### S7.13 Far references added after user approval (written 2026-10-15 ~1:20 am PT)
- The user approved the 1-sps download at ~0:30 am PT (relayed by the coordinator ~1:15 am PT). Fetched from EarthScope
  FDSN (`scripts/scan_fetch_ref2011.py`, log `logs/ref_download.log`, raw `data/raw_ref2011/`, StationXML
  `data/resp_scan/NET.STA.xml`): LHZ 2011 for IU.COR, US.NEW, US.HAWA. LHN/LHE not fetched (3× the volume).
- CN.PNT and CN.BBB have no LH channel at NRCan or EarthScope (BH 40 Hz only; a year of BHZ ≈ 2 GB per station, far over
  the approved ~150 MB). Substitutes with LHZ at EarthScope, as the approval allows: UW.OMAK for PNT (same azimuth,
  ~315 km) and TA.I02D for BBB (no LH station exists north of the network in 2011; I02D is coastal Oregon, ~470 km S).
  The northern azimuth therefore has no far reference.
- Status of veto numbers at this moment: veto COLUMNS had been computed (not inspected) for the cleaned runs and the
  pilot-1 `real` run with the on-disk references; no veto-validation, rate or classification number had been computed or
  seen (the first validation pass crashed before the veto step). The far set is added to the veto references (same rules,
  role 'veto', Z only), all columns are recomputed, and every veto result is reported with and without the far set.

## Amendment S8 — pilot 3: the ratio veto fixed on 2011, tested on SVI calendar 2012 (written 2026-10-15 (Thu) ~1:55–2:15 am PT, before anything is run on 2012)
Source: coordinator instruction (track A, 2026-10-15 ~1:50 am PT) after the pilot-2 NO-GO; Merlin's fallback paragraph in
`MERLIN_SCAN_NOGO_REVIEW.md` runs in parallel as PLAN.md Amendment 7. When this was written no 2012 scan file existed:
no 2012 prep v3, band-pass, mask, detector run, column, injection, control rate or reference QC. The only 2012 facts used
are catalogue/inventory facts (raw day counts on disk, `episodes_b.csv` dates, the H&G paper text). Every rule below
that depends on a number takes that number from 2011 files that already existed (`tables/scan2_validation_events.csv`,
`tables/scan2_injection.csv`). New files: prefix `scan3_` (tables, figures), `data/scan_2012/` (runs, band-passed data,
masks, injections). The 2011 files are not touched; the shared scripts read `SCAN_YEAR` (default 2011, so every 2011
call is unchanged).

### S8.1 The ratio veto R1 (fixed from 2011 only)
- Statistic: `ref_ratio` exactly as S7.3 computes it — median over the live veto references (role 'veto', ≥ 100 km,
  QC-passed month, far set included) of A_r / A_near — with `n_ref` the number of live veto references.
- 2011 values (pilot-2 tables, existing): H&G 2 0.345, H&G 4 0.455 (n_ref 13, 14); the 8 known body-wave clusters
  0.668–1.115 (n_ref 6 each). Threshold = geometric midpoint between the larger H&G value and the smallest body-wave
  value: sqrt(0.4553 × 0.6678) = **0.551**.
- **R1: flag (travelling wave from outside) if n_ref ≥ 2 AND ref_ratio > 0.551.**
- 2011 in-sample behaviour (post hoc, stated as such, not a validation): flags 8/8 body-wave clusters, 0/2 H&G 2/4,
  5/130 (3.8 %) of the recovered standard-GF family injections with Mw ≥ 3.6 (median ratio 0.17).
- Headline also requires fp_interior (best node inside the footprint polygon), subject to the S7.12 85 % rule
  re-checked on the 2012 injections.
- Columns (never rules): R1 on the six fixed references only (GNW + the five far LHZ); R1 without the far set
  (`ref_ratio_nf`, `n_ref_nf` ≥ 2); the old V1–V3 excess variants.

### S8.2 Year 2012 set-up
- Pilot stations as 2011 (CN.PGC, CN.VGZ, CN.SNB, CN.NLLB, CN.PFB, CN.YOUB, UW.LRIV; CN.LZB noise table only); prep v3
  (`scan_prep_v3.py`, unchanged) for calendar 2012.
- Gains and templates: unchanged from 2011 (`tables/sc_gain.csv`, `templates_grid.npz`). Every pilot channel epoch that
  covers 2011 also covers 2012 (`tables/station_epochs.csv`). Note: the gains came from grand stacks that include the
  2012-08 ETS stack; they are a calibration, not a detection.
- Noise table 2012 (`tables/scan3_noise_month.csv`): the §3 rule on 2012 v3 days (Z median > 3 × the station's
  all-days median in `tables/scan_noise_ref.csv`, or > 5 × the network median that month → month excluded).
- References (same table `scan2_ref_stations.csv`, same roles and control index k): GNW, DOSE (column), CC.OBSR, UW.LEBA,
  UW.RATT, UW.ELW, TA.D04D, UW.WISH, UW.SP2 (E vetoed 2011-09-05 → 2019-05-22), TA.B05D, US.NLWA, TA.D03D — prep v3 for
  every 2012 day with raw data on disk (in 2012 GNW, DOSE, WISH, SP2, NLWA, D03D, D04D are on disk most of the year;
  the others ~49 days) — and the five far LHZ stations IU.COR, US.NEW, US.HAWA, UW.OMAK, TA.I02D, downloaded for 2012
  (user-approved reference class, ~300 MB; `scripts/scan_fetch_ref2012.py` → `data/raw_ref2012/`, StationXML
  `data/resp_scan/NET.STA.2012.xml`). Reference QC = S7.1 rules on 2012 (monthly Z noise > 3 × the station's 2012 median →
  month excluded; amplitude vs PGC on 2012 M ≥ 6.8 teleseisms in [0.5, 2], else excluded; no event → 'unchecked', used).
  Tables `scan3_refqc*.csv`.
- Decay factor c = 1 (S7.2 verdict, which already includes the 2012-08 episode).
- Masks: S7.6 tiers on 2012 (`scan_masks2.py` with SCAN_YEAR=2012; reg3 from a 2012 ComCat M ≥ 3.0 within 2,500 km
  metadata pull, `data/eq/comcat_reg2500_M30_2012.csv`); headline mask v2 = head2 ∪ bw5 ∪ reg ∪ reg3 (reg3 only if
  ≤ 15 % of 2012).
- ETS window (6a, Poisson) = the SVI major **2012-08: 2012-08-31 → 2012-09-28** inclusive (`episodes_b.csv`; the strongest
  SVI episode). Between-ETS = outside every SVI major ±2 d (2012-04, 2012-07 and 2012-08 are all majors).
- H&G 2012 events: the paper lists none. Table 1 has the 2011 and 2014 templates only; the 2011-07 → 2012-12 matched-filter
  catalogue appears only as counts in Figs 6–7; "Table S1" is referenced in the Supporting Information but is not on
  disk. The H&G item is therefore not applicable in 2012.
- HF envelopes (`scan_hfenv.py`) for the pilot stations, 2012 (class column only).

### S8.3 Runs (rules S3, S4, S5, S6, S7.3, S7.7 unchanged; pilot-1 detector without the cleaner, as the S7.12 fallback)
- `scan_bp.py` (pilot) and `scan_bp_ref.py` (references) on 2012 → `data/scan_2012/bp_2012g`, `bp2_ref`.
- `scan_detect.py`: `real` (+ 50 polarity draws, |ρ| ≤ 0.3, seed 20261014), `ctrlA_s{3,5,7,11,13}` (station k uses day
  (d + s·k) mod 366), `inj2` (`scan_inject.py` rules on 2012, seed 20261015), `inj3` (`scan_inject3.py` rules, seed
  20261016); injections are added to near and reference data.
- `scan_columns2.py` on every run (z ≥ 5 positive clusters; references shifted in control (a) as S7.3).

### S8.4 Known body-wave list 2012 (validation population, S7.4 rule)
Real run, variant 0, positive clusters with z ≥ 7 and z_loso ≥ 7, S5 filter as validated in 2012 (no class_eq), outside
mask v1; a member is "known body wave" if its time lies in a bw5 window (M ≥ 5.0, Δ ≥ 1,500 km, [t_P − 60 s, t_S + 600 s]).
The list is used whatever its size; < 5 members → result stated with the small-number caveat.

### S8.5 Validation of R1 in 2012 (out of sample)
- (ii) ≥ 80 % of the 2012 known body-wave list flagged (ceil);
- (iii) ≥ 85 % of the recovered standard-GF family-location injections with Mw ≥ 3.6 not flagged;
- (i) H&G: not applicable (no 2012 events listed).
Both (ii) and (iii) → "R1 validated". Otherwise "R1 unvalidated" and the veto item fails.

### S8.6–S8.8 Filter, threshold pair, GO / NO-GO
- S5 filter re-validated on the 2012 injections (S7.7: Mw ≥ 3.4, 85 % rule, share_excess fallback); fp_interior by S7.12.
- Threshold pair: S7.8 grid and selection rule (max Mw-3.8 family recovery with (a) ≤ 0.3/yr, (b) ≤ 0.7/yr between ETS;
  ties → smaller (a) + (b), then larger z_thr), **with R1 in the headline**. The same grid is also computed without R1
  (column) → separation curves (real between-ETS, real in ETS, (a), (b) vs z_thr at each zl) with and without R1.
- **GO if all:** (a) + (b) ≤ 1/yr between ETS at the chosen pair (strong ≤ 0.5); ≥ 80 % recovery of the family Mw 3.8
  injections with median location error ≤ 15 km; ETS-window headline detections ≥ 3 with Poisson p < 0.01 against
  the between-ETS real rate (S7.9 formula, floor 0.5 / between live time); R1 validated (S8.5). H&G: N/A.
  The between-ETS count is not a criterion. NO-GO → the fallback (PLAN.md Amendment 7, running) answers the question.
- Columns reported: the 2012 rate of R1 flags among control (a) clusters (how often noise looks "from outside"), real
  headline detections with PT times, the 2011-08-28 01:08 UTC 2011 candidate's R1 status (ratio 0.36 in 2011: not flagged).
- Outputs: `tables/scan3_*` (verdict `scan3_verdict.json`), `figures/scan3_*` (raw first: z(t) 2012 with masks, ETS and
  detections coloured by z; separation curves with/without R1; veto plane coloured by z; injection recovery), RESULTS
  "Scan pilot 3 (2012, ratio veto)", NOTES.
- Compute: ≤ 5 CPU workers for this track (the other 3 go to PLAN.md Amendment 7), OMP_NUM_THREADS=1, < 100 GB.
