# MERLIN — design for a between-ETS VLFE scan ("VLFEs without LFEs")
Written by Merlin (read-only advisor, Fable 5.1) 2026-10-14 ~3:30 pm PT; saved verbatim by the main session.
Status: design only. Pilot (SVI 2011) to start after the S13 per-episode scale-up finishes (step 3 needs its grand stacks).

## Verdict

The scan is viable and worth doing, but **not with the S13 stacks as matched-filter templates**: per-component stack SNR is only 2–3 (`tables/pilot_2012-08_stacomp.csv`: VGZ Z peak 1.24 nm vs control-p99 0.93; LRIV Z 1.80 vs 0.92; SNB Z 1.98 vs 0.96; `figures/s13_pilot_2012-08_pooled.png`), and the thrust synthetic already explains everything a free moment tensor explains (VR 0.39 vs 0.40, `RESULTS.md` line 13). The one decision that matters: **detect with the synthetic thrust template on a source grid; use the ETS stacks only to calibrate per-station-epoch gain and polarity; set thresholds from two scrambled controls plus the sign-symmetric control.** Single-event completeness will be about Mw 3.2–3.4 in summer and 3.6–3.8 in winter (numbers in step 6d), i.e. the upper half of the Mw 3.1–4.3 range PLAN.md §0 attributes to H&G. The data already hint at the answer: at the four 2011 SVI-area H&G template times our catalogue has 3, 1, 1, 0 LFEs within ±60 s (expected 0.44, 0.44, 0.17, 0.55; `tables/hg2018_crosscheck.csv`), whereas the ETS partition (5e12 N m per LFE) predicts ~20–40 LFEs per Mw 3.3–3.5 VLFE. Large discrete VLFEs and the LFE-locked VLF moment are probably two populations even inside ETS; the scan quantifies this between ETS.

## Flaws / flags (severity-ordered)

**F1 — MAJOR (fatal for "stacks as templates"; harmless for S13 itself).** Stack SNR 2–3 per component (above). A template whose noise rivals its signal loses most matched-filter gain and fires on its own wiggles. Synthetic explains the stacks fully (thrust VR = free-MT VR). Fix: synthetic templates; stacks only for calibration (step 3).

**F2 — MAJOR for a continuous scan.** Gap-edge transients. `s13_stack.py` line 23 `DIL = 150` s; `load_station` lines 33–40 zero-fill NaNs before the zero-phase band-pass. Measured on `G.bandpass` at 1 sps: impulse response > 1e-3 of peak spans ±240 s; a unit step at a gap edge leaves 0.44 % at 150 s, 0.15 % at 200 s, 0.017 % at 300 s. Prep (`s13_prep.py` line 134) NaNs 600 s at segment edges but does not taper the retained data into the NaN, so the edge value is the detrended long-period drift (100s of nm to µm on horizontals) → 1–10 nm transients just outside the 150-s mask = a Mw 3.0–3.3 VLFE. Averages out in 1,400-window stacks; a false-alarm generator for single events (PGC 2012: gaps every few minutes, Amendment 2). Fix in step 2.

**F3 — MAJOR for a scan, MINOR for S13.** Day-edge attenuation. `s13_prep.py` line 125 `tr.taper(0.05)` = 5 % at EACH end (obspy semantics verified) = 4,500 s of the 90,000-s padded segment; line 131 `remove_response(taper=True, taper_fraction=0.05)` = obspy `cosine_taper` p=0.05 → 2.5 % each end = 2,250 s. Only 1,800 s padding + 600 s NaN cover them, so the first and last ~2,700 s of every UTC day are multiplied by a ramp from ≈0.31 at midnight to 1. For S13 stacks ≈2 % mean amplitude bias (ignore); for single events, ~6 % of each day at 1/3 sensitivity and biased amplitudes near midnight. Fix in step 2.

**F4 — MAJOR for S13's interpretation.** Window QC (`stack_windows` line 59: drop if window RMS > 3 × median) rejects any window containing a VLFE above ≈ Mw 3.6 at VGZ / 3.8 at PGC (a one-cycle pulse of amplitude A adds ≈0.22 A to a 400-s window RMS; thresholds 39 nm VGZ, 69 nm PGC from `thr_rms_nm`). S13's "VLF moment per LFE" is therefore the SMALL, CONTINUOUS accompaniment with large discrete VLFEs excised — a lower bound on time-locked VLF moment; bears on A3 (vs GNSS). State in RESULTS.md. S13 numbers cannot predict discrete-VLFE counts.

**F5 — MAJOR for location, MINOR for detection.** Whole-space GF with flat ×2 (`s13_green.py` lines 24–53). At r 30–60 km, T 20–50 s, near/intermediate and far-field terms are comparable (λ/2πr ≈ 0.3–1), so free surface and layering change RELATIVE station amplitudes by order-1 factors — exactly what location uses. Polarity/shape are validated by the stacks; the amplitude pattern is not. Calibrate location empirically (6e) and if biased switch location templates to half-space Okada (`pyrocko.modelling.okada` imports; `cutde` installed; no pyrocko GF store: `gf_store_superdirs = []`).

**F6 — MINOR.** Template source is a 1-s step (`synth`, `tau=1.0`); Mw 3–4 VLFEs last 10–50 s → M̂0 biased low by tens of % for long sources. Use 1-s (headline) and 20-s (column) durations.

**F7 — MINOR, mandatory in pipeline.** Epoch noise jumps: SNB Z median window RMS 135 nm in 2011-08 (`pilot_2011-08_stacomp.csv`) vs 12 nm in 2012-08; LZB ~1,000 nm. Per-station-month noise table and per-epoch gain calibration required.

**F8 — SPECULATIVE (paper not on disk).** All H&G numbers come from PLAN.md §0 and the 12 template times in `s13_hg_crosscheck.py` lines 8–13; their Mw, locations, method and between-ETS rate are not on disk. Transcribe Table 1 to `tables/hg2018_table1.csv` with a source line before using as truth (ask the user for the PDF if needed). Template 1 (47.825 N, −123.075) is WA-side; templates 6–12 (2014-12) at 48.7–49.1 N are north of the SVI family cloud (families 48.12–48.93 N, −124.34 to −123.10, 28.8–45.9 km; `families_v1_3.csv`).

**F9 — MINOR/known.** `scale_svi_episodes.csv` and `scale_verdict.json` on disk are from the broken 12:24 pm run (NaN from 2004-07 on) although unit files exist: 2004-07 1.95e12 N m z 4.9; 2005-09 3.7e12 z 6.6; 2012-08 4.7e12 z 7.8; 2013-09 5.0e12 z 5.8; 2018-06 6.2e12 z 9.5; 2024-08 5.0e12 z 9.1. `scale_svi_between_units.csv`: SVI_BE_broad (688 LFEs 2014–26) 4.3e12 per LFE, z 3.2, VR 0.07, ±1-day −1.8e12/−0.7e12 → between-ETS LFEs DO carry the VLF accompaniment (provisional until `scalef_*`). SVI_BE (76) z −1.9 (min detectable ≈1.4e13 per LFE).

## Prescription

Design principle: at 20–50 s a VLFE is a band-passed step whose waveform is nearly identical at every station; the discriminating information is (i) the sign/amplitude pattern across station-components (uplift on all nearby verticals, trenchward horizontals; `NOTES.md` line 47), (ii) isolation in time, (iii) absence of 2–8 Hz energy. Build the detector around those, not waveform shape.

### 1. Pre-register and build truth tables (half a day, no compute)
- `PLAN_SCAN.md` in `svi_science/results/s13_vlfe/` with every rule below fixed before any number exists; amendments dated PT.
- `tables/hg2018_table1.csv` (F8). LFE-catalogue footprint polygon: 10-km-buffered hull of SVI families with ≥ 30 catalogued LFEs at P ≥ 0.75 over 2003–2026 → fixed "inside/outside" test for step 8.
- `tables/station_epochs.csv` from `data/resp/*.xml` BH/HH channel epochs (LFE-channel version exists at `lfe_extension/results/cut/svi_ext/channel_epochs.csv`).
- Earthquake masks: teleseismic from `data/rf_teleseis/events_M58_1993_2026.csv` (8,086 events, M ≥ 5.8) now, plus a ComCat global M ≥ 4.5 pull 2003–2026 (tens of MB, inside the download rule). Why M 4.5: by the Ms relation an Ms 4.5 at 60–90° gives 20-s Rayleigh amplitudes ~200–350 nm, 5–10× a Mw 3.3 VLFE (40–55 nm at PGC/VGZ/LRIV). Regional: `lfe_extension/results/trackC/eq_catalog.csv` (26,464 events 2009–2026, 46.5–50 N, M ≤ 4.9) plus ComCat regional M ≥ 2.5 for 2003–2008.
- Mask lengths: origin → origin + Δ(km)/2.5 km/s + 20 min (M < 6); + 60 min (6–7); + 3 h (7–8); + 24 h (≥ 8, multi-orbit). Regional M ≥ 2.5 within 300 km ±5 min; M ≥ 4 within 1,500 km ±20 min. Report masked fraction (expect 3–8 %).
- Acceptance: four files exist; PLAN_SCAN.md lists thresholds, windows, controls, pass rules verbatim.

### 2. Prep v3 for the pilot year (≈25 min)
- Pilot: SVI, calendar 2011 (2011-08 ETS; 11 months after 2010-08b1; H&G templates 2–5; winter+summer). Stations CN.PGC, CN.VGZ, CN.SNB, CN.NLLB, CN.PFB, CN.YOUB, UW.LRIV (all 2011 on disk per `tables/inventory_bb_near.csv`); CN.LZB enters the noise table and is excluded by the step-3 noise rule, not by hand.
- New root `data/disp1sps_v3/` (env `S13_DISP`); v1/v2 untouched. Prep changes: (a) replace `tr.taper(0.05)` by a fixed-length taper ≤ 600 s and call `remove_response` with `taper=False`, so no taper reaches inside the day (F3); (b) save per-day gap intervals so the scanner cosine-tapers 120 s into every gap edge before band-passing and dilates masks by 300 s not 150 (F2); (c) keep `S13_MAXGAP` 15 s.
- Throughput from `logs/prep.log`: median 112 station-days/min at 8 workers → 2,555 station-days ≈ 25 min, ~1 MB each.
- Acceptance: PGC 2011-08-20 (gap-free) v3 vs v2 agree < 1 % except the first/last 2,700 s where v3/v2 falls from ~3 to 1 (that IS the F3 ramp). Failure: ratio ≠ 1 mid-day → stop.

### 3. Templates and station calibration
- Grid: nodes every 10 km, 47.8–49.3 N, −124.8 to −122.6 E (~200), depth = Slab2 at node (`data/slab_models/slab2/cas_slab2_dep_02.24.18.xyz`). Thrust 315/12/90 headline; strike-slip and steep-dip competitors; 5-element deviatoric basis for the post-detection free fit; durations 1 s and 20 s. Use existing `templates()`/`elem_templates()`; store once.
- Calibration (the legitimate use of S13's output): grand stack per station-component-epoch = n-weighted sum of per-episode pooled stacks in `data/stacks/scale_svi_*.pkl` (~25,000 LFEs; per-component SNR ≈ 8–12 at PGC/VGZ/LRIV/SNB). Project onto the synthetic at the pooled centroid → gain g and sign per station-component-epoch. Rules: |log g| ≤ ln 2 and positive → keep, weight scaled by g; otherwise exclude that station-component-epoch (LZB, SNB-2011 expected). Write `tables/sc_gain.csv`.
- Noise table per station-component-month (median band RMS over 400-s windows); exclude a month if Z median > 3 × that station's 2003–2026 median or > 5 × the network median that month.
- Empirical-template option (strongest form of the user's idea): if any 2011 ETS-time detection or H&G event has Mw ≥ 3.8 (≥ 500 nm at PGC/VGZ/LRIV, SNR ≥ 30 per component), cut it as an empirical all-station 3-component template (true Earth response) and run it as a second detector; report both lists. Below Mw 3.8 the synthetic is better.
- Acceptance: g in [0.5, 2] and positive for ≥ 12 of 21 pilot station-components; Z gains at PGC/VGZ/LRIV in [0.7, 1.4]. Failure: Z gains all ≈0.5 or ≈2 → the ×2 free-surface factor is off; carry the factor into the moment scale, continue.

### 4. Detector (sliding `fit()`)
- Data: v3, gap-tapered, 0.02–0.05 Hz zero-phase on contiguous runs, masks dilated 300 s.
- M̂0(t) per node = Σ_sc w_sc (template_sc ⋆ data_sc)(t) / Σ_sc w_sc |template_sc|² over the −60…+120 s fit window (`FIT0/FIT1`), w_sc = g_sc²/σ_sc², σ_sc = 1.4826 × MAD of the band-passed component over the UTC day. Implement as FFT cross-correlation per node per component (~200 × 21 × 3.2e7 samples/yr: minutes–1 h CPU per year; trivial on cupy).
- z(t) = M̂0(t) / (1.4826 × MAD of M̂0 over the UTC day, unmasked). Store network CC(t) = Σ w t d / sqrt(Σ w t² Σ w d²) as a column.
- Candidate = local max of z over nodes and time with M̂0 > 0 (thrust sense; negative exceedances kept as sign control). Decluster ±150 s keep max (template side-lobes: opposite sign at ±30–40 s, ~20 % positive at ±60–70 s → only Mw ≥ 4 could self-trigger; covered).
- Columns per candidate (never cuts; headline rules in step 7): n_sta finite (headline ≥ 3); n_Z_agree = verticals with own positive CC (headline ≥ 2); max single station-component share of Σ w·t·d (headline ≤ 0.5); VR_thrust/VR_free (headline ≥ 0.8); VR_thrust − VR_strikeslip and − VR_steepdip (> 0); isolation = band RMS in [−900,−300] and [+300,+900] s over daily median (headline < 2, kills surface-wave trains); near_gap_s (flag < 300); in_tele_mask, in_regional_mask; HF ratio = peak 2–8 Hz envelope within ±60 s at the two nearest stations over daily median (≥ 10 at ≥ 2 stations → class "earthquake", never VLFE); minutes to midnight; season.
- Horizontals stay in with 1/σ² weights; Z-only variant as a column; lists must agree above 1.2× threshold.

### 5. Controls and threshold (all three at thresholds 4, 4.5, 5, 5.5, 6, 7)
- (a) Station-day-shifted scrambled control: each station a different integer-day offset (e.g. +3k days for station k, same epoch), identical detector. Kills network coherence; keeps single-station artefacts.
- (b) Polarity-scrambled template control: unshifted data, random sign per station-component (50 draws). Keeps network-coherent noise (teleseisms, storms); tests dependence on the thrust sign pattern.
- (c) Sign control: count of negative-M̂0 exceedances in real data per threshold.
- Do NOT headline time reversal: under a zero-phase band-pass the near-field (step) part of the template is odd in time and just flips sign, while the far-field part is even; the reversed statistic is neither independent nor equivalent.
- Headline threshold = smallest listed at which (a) ≤ 1/yr and (b) ≤ 2/yr after masks. Report real positive, negative, (a), (b) at every threshold = the separation curve (the evidence).
- Acceptance: at headline threshold, ETS-window positives ≥ 10 × negatives. Failure: positive ≈ negative in ETS → check polarity vs grand stack; stop.

### 6. Validation against known answers (inside the 2011 pilot)
- (a) ETS re-find: detection density in 2011-08-12…09-03 vs rest of year ≥ 10×. Also total moment of ETS-time detections vs S13 time-locked total (5.3e12 × 1,445 = 7.7e15 N m) → fraction of LFE-locked VLF moment in discrete detectable VLFEs (a new number).
- (b) Mechanical consistency: mean of M̂0(t) at the pooled-centroid node, sampled at the 1,445 P ≥ 0.75 LFE origins of 2011-08 with S13 window QC, must equal 5.3e12 within 20 %. Failure > 2× → implementation bug; stop.
- (c) H&G templates 2–5 (2011-08-21 04:09:06, 08:43:29; 08-24 10:47:38; 08-26 13:44:34 UTC) must be detections within ±60 s at a node within 40 km of their location; compare Mw once Table 1 exists. Expect ≥ 3/4. Failure ≤ 1/4 while (d) says Mw ≥ 3.5 recovered ≥ 80 % → truth table or our polarity/timing wrong; debug before any between-ETS number (reproduce-the-known-answer rule).
- (d) Injection: add synthetic thrusts (whole-space GF, 1-s and 20-s) into v3 displacement before band-passing: Mw 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0 × 40 each, random unmasked times over the year, half at random family locations, half at nodes 30–50 km from every family; plus 40 at Mw 3.6 with a perturbed GF (free-surface 1.0, vs 3.4 km/s, 30-s duration). Outputs: recovery vs Mw by season (completeness), M̂0 bias, location-error median and 95 % radius vs Mw (= resolution for test 2; the only honest way with 5–7 stations). Expected from pilot noise (Z window RMS PGC 23, VGZ 13, LRIV 16, SNB 12 good-year, PFB 13 nm; template nm per 1e12 N m: PGC 0.54, VGZ 0.41, LRIV 0.55, SNB 0.39, PFB 0.23): Mw 3.3 → per-Z SNR 2.3/3.2/3.4/3.3/1.8, network ≈ 6 → ~50 % at Mw 3.2–3.3, ≥ 80 % at 3.5–3.6 in summer; winter 2–4× noisier → +0.3–0.5 Mw. Location: expect 20–30 km (68 %) at Mw 3.3–3.6, 10–15 km at ≥ 4; depth unresolvable (fix to Slab2).
- (e) ETS-stack location test (settles F5): grid-search the VR-maximising node for the 2012-08 pooled stack (z 7.8) and 2018-06 (z 9.5). Acceptance: within 20 km of the count-weighted centroid (48.37 N, −123.61). Failure > 30 km or a ridge along station geometry → switch LOCATION templates to half-space Okada static × step response; re-run.

### 7. Pilot go/no-go (SVI 2011)
- GO if all: 6a ≥ 10×; 6b within 20 %; 6c ≥ 3/4; 6d ≥ 80 % at Mw 3.6 and ≥ 50 % at 3.3 (summer), location median ≤ 30 km at 3.6; controls ≤ 1/yr (a) and ≤ 2/yr (b); 6e within 20 km.
- The between-ETS count is NOT a criterion: zero is a result (rate bound at measured completeness) if the above hold.
- NO-GO signatures → adaptations: (i) H&G not re-found though injection says they should be → truth-table or polarity error, debug; (ii) control (b) > 1/month at the Mw-3.6 threshold → coherent long-period noise: require HF-quiet + isolation in the headline, raise the low corner to 25 s, add era stations (C8 POLARIS 2003–06, XU CAFE 2006–10, UW compact 2021+ in `inventory_bb_near.csv`); if still failing, restrict to Mw ≥ 3.8 and fall back to the stack-at-candidate-times path (S13 A5 extended to P 0.5–0.75 candidates and to PNSN tremor times between ETS) — answers the complementary question at coarse resolution; (iii) completeness worse than Mw 3.8 in every season → same fallback.

### 8. The three explanations — pre-set tests (every between-ETS detection; between-ETS = outside every `episodes_b.csv` SVI major ±2 d)
- Branch 0 first (false detection): headline drops anything in a mask, HF ratio ≥ 10, isolation ≥ 2, or single share > 0.5; all kept in the table with flags.
- Test 1 — LFEs present below the catalogue bar (each vs 500 same-day random-time draws):
  - 1a catalogue counts: SVI candidates (`data/catalogue/svi_catalogue_v1_4_backward.parquet` 2003–13: 309k rows P < 0.5, 37k at 0.5–0.75, 52k ≥ 0.75; `svi_science/results/s0/analysis_table.parquet` prod/keep_xfam8 2014+: 480k P < 0.5, 10k ≥ 0.75) within ±60 s and ±300 s, three P classes, summed over detections. Pass "present" if excess ≥ 3 σ_control in any class/window. ETS partition would predict ~20–40 P ≥ 0.75 LFEs per Mw 3.3–3.5 VLFE; H&G cross-check gives ~1 → expect small/absent.
  - 1b sub-threshold matched-filter stack: run the SVI family templates (`lfe_extension/results/a4_templates_merged_*.npz`; matched-filter machinery `lfe_extension/scripts/a4_run_matched.sh`, `c_common.py`) over ±600 s of 2–8 Hz data at the 5 nearest families; stack CC² aligned on the VLFE origin over all detections vs 200 scrambled-time stacks. Pass if ±60 s mean > control p99. Gain ~1/√N below the single-event bar.
  - 1c HF envelope stack: 2–8 Hz log-envelope (daily-median normalised), 3 nearest stations, ±600 s; same rule. Catches tremor without family match.
  - Report LFE-to-VLFE moment ratio = (excess P ≥ 0.75 count × 4.5e11) / Σ M̂0_VLFE vs the ETS value 1/12 (`RESULTS.md` line 28).
- Test 2 — off the LFE families: per detection, best node and Δz ≤ 1 region; 68 %/95 % radii from injection at that Mw and season. Classes: OUTSIDE footprint (95 % region entirely outside polygon); INSIDE (best node inside and 68 % region overlaps a ≥ 30-LFE family); UNDETERMINED. State plainly: inside the footprint no single family can be assigned at 20–30 km; the resolvable statement is footprint in/out and along-strike tercile (S13 zones). OUTSIDE (e.g. north of 48.7 N where H&G 2014 templates sit) says nothing about LFE occurrence — no templates there.
- Test 3 — genuinely smooth slip at a family location: pooled over INSIDE, branch-0-clean detections: pass if 1a, 1b, 1c all < 2 σ AND the pooled bound (3 σ excess × 4.5e11 / Σ M̂0) < 1/48 (quarter of the ETS ratio). Statement: "between-ETS VLFEs of Mw ≥ completeness, inside the LFE footprint, carry ≤ X of the ETS LFE accompaniment". Else UNDETERMINED with the bound. Expected: either passes, or 1a shows ~1 LFE per VLFE (≈1/200 of the ETS ratio).
- Present both directions: S13 A5 (between-ETS LFEs carry 4.3e12 per LFE, z 3.2, provisional) and this scan (discrete between-ETS VLFEs carry [no/few] LFEs).

### 9. Scale-up and costs
- Full SVI: PGC/VGZ/NLLB 2003–26, SNB 2003–21, LRIV 2009–26, PFB 2009–18, YOUB 2010–14, C8 2003–06, UW compact 2021+: ≈34,000 new station-days ≈ 5 h prep, ~35 GB; scan + 2 controls + injection: minutes–1 h per year CPU.
- WA-NW: GNW/DOSE/WISH/SP2 + XU 2006–10 + YW: ≈23,000 station-days ≈ 3.5 h; H&G template 1 (2011-08-17, 47.825 N) is the WA truth; families `lfe_extension_wa/results/wa_prod/setup/families.csv` (lat_use/lon_use/depth_km); candidates `lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet` (60 M rows; selected_v2 42,696; P_own_pnet ≥ 0.75: 39,122).
- Final resolution: per region, per year / inter-ETS interval: counts and moment of between-ETS VLFEs by class (INSIDE/OUTSIDE/UNDETERMINED × LFE-accompanied/not) with seasonal completeness. Coarse, and new.

### 10. Reporting
- Raw first: z(t) for the pilot year with masks shaded and ETS marked; 2011-08 panel with H&G times as ticks; detection map coloured by class with the footprint polygon; separation curves (real +, real −, control a, control b vs threshold). Points coloured by CC, no error bars.

## Experiments Merlin could not run
1. Grand-stack gain table (step 3): |log g| ≤ ln 2 everywhere → proceed; systematic 0.5/2 → moment-scale factor only.
2. ETS-stack grid search (6e): within 20 km → whole-space GF adequate for 20–30 km location; > 30 km → Okada half-space for location.
3. Injection completeness by season (6d): ≥ 50 % at Mw 3.3 summer → H&G upper half accessible; ≥ 50 % only at ≥ 3.8 → stack fallback.
4. H&G re-find (6c): ≥ 3/4 → validated against external truth; ≤ 1/4 → stop.
5. Control (b) rate: ≤ 2/yr at the Mw-3.6 threshold → single-event detection viable; ≫ → coherent noise dominates, step-7 adaptations.

## What would change Merlin's mind
- Grand stack reaching per-component SNR ≥ 10 at PGC/VGZ/LRIV/SNB with VR_free − VR_thrust > 0.1 → an empirical grand-stack template becomes the better detector; F1 lapses.
- Control (b) ≤ 2/yr at headline threshold with no earthquake masking → relax the mask to M ≥ 5.5.
- Injection location medians ≤ 10 km at Mw 3.3 → test 2 can be per family instead of per footprint/tercile.

Key paths: `svi_science/results/s13_vlfe/{PLAN.md,RESULTS.md,NOTES.md}`, `.../scripts/{s13_stack.py,s13_stack2.py,s13_prep.py,s13_green.py,s13_common.py,s13_hg_crosscheck.py}`, `.../tables/{pilot_2011-08_stacomp.csv,pilot_2012-08_stacomp.csv,hg2018_crosscheck.csv,inventory_bb_near.csv,scale_svi_between_units.csv}`, `.../data/stacks/scale_svi_*.pkl`, `data/rf_teleseis/events_M58_1993_2026.csv`, `lfe_extension/results/trackC/eq_catalog.csv`, `data/slab_models/slab2/cas_slab2_dep_02.24.18.xyz`, `svi_science/results/b_ets_record/episodes_b.csv`, `svi_science/results/s7/between_events.csv`.
