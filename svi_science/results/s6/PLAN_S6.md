# S6 — LFE size per detection: PLAN (frozen before any S6 number)

Written 2026-10-07 ~1:10 pm PT, before any amplitude was measured. Only counts were looked at while writing it
(data P ≥ 0.75: 13,466 raw / 10,378 counted once; P ≥ 0.5: 30,411 / 22,378; all match a cut-store detection exactly).
Scripts `svi_science/scripts/s6_*.py`; outputs `svi_science/results/s6/`, figures `svi_science/figures/s6_*.png`.
Times shown to readers in Pacific time; joins in UTC.

## Question
How big is each southern Vancouver Island LFE, and does LFE size change with (1) position in the local episode
(at the front vs behind it), (2) episode size, (3) between-episode vs ETS, (4) the 12 years, (5) depth?

## Unit
- Headline: `s0/analysis_table.parquet`, run prod, series data, P_same ≥ 0.75, keep_xfam8 (counted once).
- Columns: P ≥ 0.5 and P ≥ 0.85; raw (not counted once). Nothing is dropped; every row keeps its flags.
- Non-event reference (no scrambled-control windows were cut): 6,000 random data candidates with P_same < 0.05
  (seed 0). They are measured with the same code and give the amplitude a noise "detection" reaches.

## Measurement (per detection × station)
- Windows from `lfe_extension/results/cut/svi_ext/` (−15…+30 s around the expected S; S at +15 s from the template
  offset; `read_window` of `lfe_south/scripts/cut_windows.py`). All 3 components of the station's band.
- Demean, 2–8 Hz Butterworth-4 zero-phase at native rate, divide by the channel-epoch sensitivity
  (`channel_epochs.csv` gain, counts per m/s) → ground velocity. Gapped samples (zeros) are masked.
- Stored: 3-C mean-square envelope in 0.1 s bins over the whole 45 s (so windows can be audited later).
- A_S = 3-C RMS in the S window [−0.5, +3.5] s about the expected S. A_N = 3-C RMS in [−14.5, −5] s (the first 0.5 s
  is dropped for the filter edge). S/N = A_S / A_N.
- S-window check (per family × station): the median-envelope of the family's headline detections must peak inside
  [−1.5, +4.5] s; family-stations that fail are not eligible as best stations (flag, counted in the results).

## Best stations and relative size
- Per family, eligible stations ranked by the median S/N over its headline detections; best = the top 5 with
  median S/N ≥ 2 and present for ≥ 50 % of the family's headline detections; if fewer than 3 qualify, the top 3 by
  median S/N (family flagged `weak_family`).
- Reference per family × station: ref = median log10 A_S over the family's headline detections (all years, after the
  gain correction). r(d,s) = log10 A_S(d,s) − ref(f,s).
- **ΔM(d) = median of r(d,s) over the family's best stations present** (n_used column; ΔM is set only when n_used ≥ 2).
  ΔM is a relative amplitude magnitude (log10 amplitude ratio to the family's median repeat).
- Noise level ΔN(d) = the same median of log10 A_N − refN(f,s). Flag `low_sn` when the median S/N over the used
  stations is < 2 (kept, flagged).
- Absolute (approximate): M ≈ mag_med(family, Bostock, `e5r_families_2012.csv`) + ΔM. Stated as approximate:
  Bostock's median refers to his 2003–2013 catalogue, ours to our 2014–2026 repeats.

## Station gain epochs
- Sensitivity correction from `channel_epochs.csv` (NLLB BH→HH 2017-02-23 ×1.77; PGC 2015-01-12 and BH→HH 2017-08-10;
  VGZ 2017-03-03; SNB 2017-06-08; LRIV 2019-05-30).
- Check 1 (as asked): median r(d,s) of all headline detections at that station in the 30 d before vs after each change.
- Check 2 (cleaner): station residual e(d,s) = r(d,s) − ΔM(d) computed from the OTHER best stations; per station ×
  epoch median. A residual step > 0.05 log units (12 %) at a change = sensitivity correction incomplete.
- Correction rule: the per-station-epoch residual median (Check 2) is ALWAYS subtracted and ΔM recomputed (one pass).
  It is relative to the other stations, so it cannot remove a size change seen by the whole network. Both the
  sensitivity-only and the epoch-corrected ΔM are kept; the epoch-corrected one is headline.

## Winner's curse and completeness
- Test: Spearman of ΔM vs net_mad within family (ranks within family), headline rows; the same for the non-event
  reference.
- Completeness: per family, z = ΔM − ΔN (size relative to the local noise). Mc_z(f) = the maximum-curvature bin
  (mode of the 0.1-wide histogram of z) + 0.1, families with < 30 headline detections use their group's pooled z.
  Per family-month completeness amplitude Mc(f,m) = Mc_z(f) + median ΔN of that family-month's detections.
  `complete` = ΔM ≥ Mc(f,m). All rows kept; every test is run on complete rows (headline) and on all rows (column).

## Tests (headline = complete, P ≥ 0.75, counted once; families weighted equally)
Effect sizes in ΔM units (log10 amplitude). Bootstrap = 2,000 resamples of families. Permutation = 2,000.
1. **Front vs behind.** The 10 SVI S3 episodes (`s3/episodes_kinematics.csv`). Front arrival at a family = the
   10th-percentile time of PNSN tremor in the family's 5 km along-strike bin (≥ 20 epicentres; the s3b rule).
   lag = event time − arrival. Classes: front = −1 ≤ lag < 2 d; behind = lag ≥ 2 d. Effect = mean over family × episode
   cells with ≥ 3 events in each class of [median ΔM behind − median ΔM front]. Permutation: class labels shuffled
   within each family × episode cell. Column: Spearman ΔM vs lag (within-cell ranks).
2. **Big vs small episodes.** The 39 MAJOR episodes; size = PNSN family_exposure_days (`episodes_svi_final.csv`,
   independent of our detector); big = above the median. Effect = mean over families with events in both of
   [median ΔM big − median ΔM small]. Permutation: big/small labels shuffled across episodes. Column: Spearman of
   each episode's median within-family ΔM vs log size; and vs S1 Tier A geodetic offset D on the clean episodes.
3. **Between-episode vs ETS.** Between = the 76 S7 events (`s7/between_events.csv`); ETS = events in MAJOR episodes.
   Effect = mean over families with ≥ 1 between event of [median ΔM between − median ΔM ETS]. Permutation: each family's
   between events replaced by the same number of its own ETS events drawn at random. Column: all 734 events outside
   episodes.
4. **Twelve years.** Per family with ≥ 20 events spanning ≥ 5 years: Theil–Sen slope of ΔM vs time (per decade).
   Effect = mean slope. Permutation: event times shuffled among the family's episodes (whole episodes moved, so the
   within-episode structure stays). Noise check (no control needed): the same slope for ΔN must not trend; the step
   checks above must pass. Figure: yearly median ΔM per group.
5. **Depth.** Across families: Spearman of (a) Bostock mag_med and (b) our distance-corrected family size M_f vs
   depth. M_f from log10 A_S(f,s) (family-station medians) = M_f + S_s − n·log10 R_fs (station terms, one spreading
   exponent n, least squares). Permutation: depths shuffled among families within 10 km along-strike bins.
   Column: within-family spread (IQR of complete ΔM) vs depth.

## Pass rule (per test, fixed here)
"Seen" = permutation p < 0.01 AND bootstrap 95 % CI excludes 0 AND the effect keeps its sign on all rows AND the same
test on ΔN (noise level) is smaller than the ΔM effect. Otherwise "not seen", reported with its CI as a bound.
Surprising results are provisional until Washington-NW repeats them.

## Not done here
Tidal phase (needs the S5 tide setup). Moment / b-value fitting beyond the completeness cut.

## Deviation 1 (2026-10-07 ~1:45 pm PT, before any test statistic was computed)
Seen while checking the measurement (on 52 of 153 months): the 4 s RMS S window is buried in the tremor/noise —
median S/N per family × station = 1.17 (IQR 1.08–1.32; no family has 3 stations with S/N ≥ 2, so all families would be
`weak_family`). An RMS amplitude at S/N ≈ 1.2 measures mostly the noise level, not the LFE. The median envelope does
peak at the expected S (0–0.5 s after it, power ×2–4.5 at PGC/B009/VGZ/SNB), so the alignment is right; the S pulse is
~1 s long. Change (estimator only; units, tests, pass rules unchanged):
- **Headline estimator = matched-filter amplitude ratio** to the family reference waveform: per family × station the
  reference T = linear stack of the family's headline detections (3-C, 2–8 Hz velocity, resampled to 20 Hz, window
  [−0.5, +3.5] s about the expected S). Per detection a(d,s) = Σ x·T / Σ T·T over the 3 components (best lag within
  ±0.25 s). ΔM_mf = median over the best stations of log10 a − ref(f,s) (a ≤ 0 → that station's ratio is set to the
  floor log10(σ_noise) and flagged). Noise: the same projection on two pre-arrival windows [−14.5, −10.5] and
  [−9.0, −5.0] s → σ_noise(d,s) = RMS of those projections; S/N_mf = a / σ_noise; `low_sn` when the median S/N_mf < 2.
- Best stations for ΔM_mf: ranked by the median S/N_mf of the family's headline detections (same ≥ 2 / ≥ 50 % / top-5
  / fallback-3 rule), using only family-stations that pass the S-window check.
- Station-epoch correction, completeness (z = ΔM − ΔN with ΔN from σ_noise), and all five tests: unchanged, run on ΔM_mf.
- The pre-registered RMS ΔM is kept as a column (`dM_rms`) and every test is also reported on it.
- The stack includes the detection itself (leave-one-out not needed: ≥ 20 events per stack; families with < 20
  headline detections use all their P ≥ 0.5 detections for the stack, flagged).
