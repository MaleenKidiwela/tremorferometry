# S10 pilot — daily active-day LFE stacks for coda dv/v (rules fixed before any number)

Written Sun 2026-10-11, 8:40 am PT, before any stack or stretch was computed. Base design:
`../MERLIN_DESIGN_2026-10-11.md`; its rules hold except where the user's revised design (2026-10-11) changes them:
no N ≥ 100 / N_eff ≤ 150 acceptance gates; the unit is one stack per (family, station, active day), weighted, never
gated; gates are columns. POLARIS addition (user, same day) in §9. "Scrambled control" = the same detector run with each
station site's data shifted to another time of the same day. eps = −dv/v throughout (positive eps = slower).

## 1. Rows (catalogue, read only)
- 2014-01 → 2026-09: `svi_science/results/s0/analysis_table.parquet`, run `prod`, `keep_xfam8` (counted once),
  P = `P_same`.
- 2003-01 → 2013-12: `data/catalogue/svi_catalogue_v1_4_backward.parquet`, `primary` & `keep_xfam8`, P = `P` (era
  detector rows 2003-02 → 2006-06 carry the era P; flagged, used as weights as the user asked).
- Families 10, 87, 17. Day = UTC date (as in the catalogue); times reported in PT.
- **Active day**: Σ P over counted-once data rows with P ≥ 0.5 ≥ 3 for that family-date.
- Data windows on an active day: every counted-once data row with P ≥ 0.5, weight P.
- P ≥ 0.75 stack (era 0.90): same day, those rows only, equal weights (gate c).
- Scrambled-control stack: every counted-once c1 row of that family-date (any P, because c1 rows with P ≥ 0.5 are
  ~2 % of the catalogue), equal weights, randomly subsampled (seed fixed) to at most the day's number of data rows.
  Days with < 3 c1 windows at a station get no control value (counted in the coverage column).

## 2. Windows
- Expected S time at station s: data row `t + offset(fam, s)`; offset = frozen template placement + 0.5 s
  (`lfe_extension/results/cut/svi_ext/offsets_used.csv`; for the POLARIS/LZB stations the same placement computed from
  the era templates `backward/templates/kelb/g3_era_templates.npz`, which reproduces offsets_used at PGC/B-stations).
- Scrambled-control row: the detector frame index i of the c1 row (`feat/<run>/<date>.parquet`), the day's tref/NVo
  (`<date>_fd.parquet`) and the site shift (`<date>.json`): S time at s = t + (((i + shift_site) mod NVo) − i)/20 +
  offset. B009/B010/B011 share the PGC site shift.
- Source: data rows 2010-08 → 2026 from the cut store `lfe_extension/results/cut/svi_ext` when the row is in it
  (merge on family and time ± 0.1 s); every other window (pre-2010-08 rows, rows not cut, all scrambled-control rows)
  from the raw miniSEED day files (`lfe_extension/traces_backward` ≤ 2013, `lfe_extension/traces` ≥ 2014) read at the
  same offsets. The cut store is sample-identical to the raw files (its VERIFY.txt); checked again here on ≥ 50 windows.
  No new cut stage, nothing written outside this folder.
- Processing (identical for data and control): window S−15 … S+10 s, demean, 5 % taper, 2–8 Hz 4-pole zero-phase
  Butterworth at native rate, resample to 40 Hz; windows with gaps/zeros in S−15 … S+10 dropped; despike = clip at
  ± 8 × MAD, MAD = median over that station-day-component's data windows (the window analogue of the day-level
  `--despike-mad 8`); per-window RMS normalisation over S−15 … S+10 s.
- Grand stack G0 per (family, station, epoch, component): P-weighted mean of all data windows of the epoch.
- S re-pick per (family, station, epoch): peak of the 0.25 s-smoothed 3-component envelope of G0 within ±2 s of the
  expected S (Z only if horizontals are missing). τ = time after the picked S.
- cc-realign: each window (data and control alike) is shifted by the lag (±0.5 s, 40 Hz, parabolic refinement) that
  maximises the 3-C summed cc of its τ ∈ [−1, +2] s segment with G0; then G1 is rebuilt, S re-picked once more, done.
- Stretch uses Z (headline). N/E are stacked and measured as columns if the run time allows.

## 3. Response epochs (never reference across a boundary)
From `lfe_catalog_repair/results/e4_station_epochs_full.csv` (and bw_gate_common BW_EPOCHS/INSTR for POLARIS):
B009, B010, B011 2005-11-02 →; B926 2007-11-13 →; PGC 2000-06-28 → 2015-01-12 → 2017-08-10 →;
VGZ 2003-02-27 → 2007-02-28 → 2009-07-31 → 2017-03-03 →; SNB 2000-05-25 → 2006-04-20 → 2008-12-11 → 2011-03-17 →
2017-06-08 → 2021-03-28. Each (station, epoch) is its own series with its own zero.

## 4. Estimator
- S-anchored stretch: the reference r(τ) is stretched, r(τ/(1+ε)), and correlated with the day stack in the window;
  ε on a grid ±1 % (0.01 % step) with parabolic refinement; ε = −dv/v. Windows τ ∈ [2, 4] s (headline) and [3, 5] s.
- Reference = the (family, station, epoch) P-weighted grand stack **leaving the measured day out**.
- Per-day uncertainty σ = SD of ε over 200 bootstrap resamples of the day's windows (weights kept).
- Pairwise doublet graph between all pairs of day stacks of one series: ε_ij measured with day i as reference;
  v_j − v_i = ε_ij, weight 1/(σ_i² + σ_j²), 4 Huber IRLS passes, zero mean per series; σ from the weighted normal
  equations.
- Per-episode stacks: all data windows of the family in one episode (`r23/tables/episodes_2003_2026.csv`, ±2 d),
  measured the same way against the leave-episode-out reference.
- Mirror window (noise only): τ ∈ [−15, −(S−P) − 1] s, S−P from the template tables. For the mirror stretch a segment of
  the same length as each coda window, ending at −(S−P) − 1 s, is time-reversed and stretched as if it lay at the
  coda lapse (same estimator, same length, pure noise).

## 5. Measured, not gated — precision vs number of repeats
In the busiest days (top 3 per family-station series by data-window count; consecutive days of the same episode pooled
when one day has < 2n windows), draw two disjoint random halves of n equal-weight P ≥ 0.5 windows, n = 10, 20, 50, 100,
200 (where available), 50 draws; σ(n) = SD(ε_A − ε_B)/√2 in 2–4 and 3–5 s. Fit σ = c·n^−β; report the n giving 0.1 %
and 0.05 %.

## 6. Gates (columns; pass rules fixed now)
- (a) causality: coda RMS (τ 2–4 s) / mirror RMS > 1.5 on each epoch grand stack (also 3–5 s as a column).
- (b) scrambled control: ≥ 90 % of control day values have |ε_c/σ_c| < 2 and the mean of ε_c/σ_c is within ±0.5;
  and |ε_mirror| < |ε_coda| on ≥ 70 % of days.
- (c) P ≥ 0.75-only vs P-weighted: z = Δ/√|σ75² − σw²| (nested samples); pass if ≥ 90 % of days have |z| < 2.
- (d) lapse: Theil–Sen slope of ε(3–5) on ε(2–4) over the days of a series ≥ 0.6 and sign agreement ≥ 70 %.
- (e) differential step injection: boundary = middle active day of each series; ε stretch a = ±0.1 % and ±0.3 %
  injected into every later day stack (the leave-day-out reference is rebuilt from the injected stacks). Recovery =
  [(post − pre)_inj − (post − pre)_natural]/a; leakage = [(last pre-boundary episode − earlier pre days)_inj − (same)
  _natural]/a. Pass: recovery 0.8–1.2 and |leakage| < 0.2 for all four amplitudes, for the daily and the graph series.
- (f) metadata overlay: every time-series panel shows the response-epoch boundaries and the known instrument / network
  changes (svi_ref_common CHANGES, cut channel_epochs.csv, BW INSTR); column = overlay drawn + any change date that
  falls inside a series.
- A series "passes" when (a)–(e) pass. Nothing is dropped for failing; every value keeps its gate columns.

## 7. Change tests (reported, provisional)
- Within-episode: per series and episode, day ε minus the episode's weighted mean, against days since the front reached
  the family (front = first date in the episode, ±2 d, with a counted-once data row P ≥ 0.75, era 0.90). Pooled weighted
  slope (%/day) ± error, same for the scrambled control and for the mirror.
- Long-term: weighted LS slope of daily ε against time per series (σ scaled by √χ²_red), and the per-episode series.
- "Above the errors" = |value| > 3σ AND outside the scrambled-control / mirror spread AND reproduced at a second station
  or window. One series alone is never a finding.

## 8. Outputs
`tables/`: rows, windows inventory, daily series (ε, σ, control, P ≥ 0.75, mirror, graph), episode series, gates,
precision-vs-N, injection, POLARIS tables. `figures/`: daily ε vs time with ETS shading and metadata lines; within-
episode ε vs days since front; σ vs n. `RESULTS.md` in plain words.

## 9. POLARIS addition (user, 2026-10-11, received before any number)
- Stations 2003–2006: C8 KHVB, MGCB, SHDB, SHVB, TWBB; PO KELB/KLNB, SILB, SSIB, TSJB, TWGB, TWKB; CN LZB; raw data
  from `lfe_extension/traces_backward`. VGZ and SNB are added per epoch for the whole record (cheap: same reader).
- Rows: era-detector rows (primary 2003-02 → 2006-06), P-weighted; P ≥ 0.90 for the equal-weight check (c).
- Epochs: KELB 2002-09-12 → 2003-07-31 (gain 6.17e7) → 2004-06-04 (gain 3.07e8) and KLNB 2004-06-04 → 2006-06-22
  are separate series first. Gain-normalised check: RMS-normalised grand stacks of the KELB(b) and KLNB series (Z);
  merge into one series only if their τ 2–4 s coda cc ≥ (the median half-split cc within each series − 0.05), their
  stretch differs from 0 by < 2σ (bootstrap), and the coda/S RMS ratios agree within 20 %. Otherwise separate.
- Reported: active days each POLARIS station adds per family in 2003–2006 (days with a stack; and days where none of
  PGC/VGZ/SNB has one); agreement with PGC, VGZ and SNB on the same days: χ² of (ε_a − ε_b − weighted mean offset)
  /√(σ_a² + σ_b²), agree if p(χ²) > 0.01, plus the weighted correlation.

## Revision 1 (2026-10-11, ~8:45 am PT — after a 3-series test run on B011, before the full run or any gate verdict)
What the test showed that the plan did not foresee:
1. **The site-shifted scrambled-control windows are real LFEs at the PGC site.** The scrambled detections are carried
   by one site: on 15 forward fam-87 days the median per-station detection z of the c1 rows is 6.0 at B011 and 5.6 at
   PGC but 0.3–0.7 at B926, VGZ, SNB, NLLB (data rows: 5.7 / 5.8 / 3.3 / 2.4 / 3.0 / 0.8). At B011 the c1 windows
   realign as well as the data (median cc 0.55 vs 0.48) and their stacks correlate with the coda reference at 0.62.
   So at the carrying station the §1 control is a second, single-site-selected sample of real LFEs, not noise; at the
   other stations it is noise. **Kept unchanged** (column `c_*`) and **added** a second scrambled control `c2_*`: the same
   c1 rows read at their detector frame time t + offset (no site shift), i.e. a same-day time unrelated to what the
   station contributed. Same reading, despike (MAD from those windows), realign and stacking path. Gate (b) is
   reported for both; the c2 version is the noise control.
2. **Injection test design check.** Against a fixed reference a single noisy day recovers the injected stretch
   (median 1.00, mean 0.96; 0.67–0.77 for the 4 days with cc < 0.5). The pre-registered leave-day-out reference mixes
   the injected and natural halves; when the stacked waveform differs between the halves (e.g. across detector eras),
   recovery drops (B011 test: fam 87 0.97, fam 10 0.50, fam 17 0.17). This is what the gate is meant to catch, so the
   rule is unchanged. Added diagnostic columns `prepost_cc24/35` (cc between the pre- and post-boundary grand stacks)
   to say why. Injection now uses cubic interpolation (was linear; 0.97 vs 1.00 on clean stacks).
3. Raw day files are missing for the 2010-08-09 → 09-12 truth window at the B stations (data windows there come from
   the cut store; scrambled-control windows cannot be read on those days).

## Revision 2 (2026-10-11, ~8:58 am PT — after the first full run, before any gate verdict was read)
- 20 % of the day values (37 % of days with ≤ 5 repeats, 8 % with > 50) sit on the ±1 % grid edge: the cc has no peak
  inside the search range. These are not measurements. Rule: |ε| ≥ 0.99 % → value kept as `*_raw`, flag `rail*`, and
  the value and its σ are set to missing for every statistic, the graph (railed days and railed pairs dropped) and the
  injection contrasts (railed natural days get zero weight). Same for the controls, the P ≥ 0.75 stacks, the
  horizontals and the episode stacks. The mirror keeps its raw value for gate (b) (a railed mirror counts as
  |mirror| > |coda|) and has its own rail flag. Rail fractions are reported per series.
- The bootstrap σ includes railed draws (conservative); it is capped by the grid (a value uniform in ±1 % has SD
  0.58 %), so σ ≳ 0.4 % means "barely constrained".
- Precision vs n: the half-split σ also saturates near 0.6–0.7 %; the fit uses only n levels whose draws are < 5 %
  railed (both the plain and this fit are reported).
- Checked and rejected: re-anchoring each day stack's S on the reference (±0.15 s) before stretching changed nothing
  (day S lags 3–5 ms SD; rail fraction and day-to-day SD unchanged at B011/PGC/SNB).
