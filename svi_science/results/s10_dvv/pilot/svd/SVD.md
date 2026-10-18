# S10 pilot — SVD-Wiener denoising of the daily stacks (rules fixed before any number)

Written Sun 2026-10-11, ~9:30 am PT, before any SVD stack or stretch was computed. User request (2026-10-11): add
the original pipeline's SVD-Wiener filter to the S10 daily-stack pilot and measure what it buys. Base rules: `../PLAN.md`
(+ Revisions 1–2). ε = −dv/v (positive ε = slower). "Scrambled control" = windows that are not LFEs of the family
(here the Revision-1 frame-time control `c2`: the scrambled detector rows read at their frame time, no site shift).

## 0. Scope
Families 10, 87, 17 × PB.B011, PB.B009, PB.B010, PB.B926, CN.PGC (E1, E2, E3 separately) = 21 series. Same windows,
same S re-pick and cc-realignment, same P-weighted day stacks (active day ΣP ≥ 3), Z component, as `s10_measure.py`
(functions imported, not re-written). Stretch grid ±1 %, 0.01 % step, parabolic refinement; rail rule of Revision 2
(|ε| ≥ 0.99 % = no peak = not a measurement) for both raw and SVD.

## 1. SVD-Wiener filter (reused unchanged from `scripts/dvv_roll30cal.py::_family`)
- Matrix R: rows = the active-day stacks of one (family, station, epoch) in date order; columns = every sample of the
  stored window (τ = −14 … +10 s about the expected S, 40 Hz, 961 columns).
- Each row divided by its peak |value| (as in the original); σ_n = SD of all entries of the row-normalised matrix in
  the "pre-noise" columns; nmf = σ_n (√n_rows + √n_cols); U, s, Vᵀ = svd(R); weights w_k = s_k² / (s_k² + nmf²);
  R_f = U diag(s·w) Vᵀ. No hard rank: the original uses this soft Wiener weight on every singular value.
- Pre-noise columns: the original stack started 3 s before the detection time with the S at +1 s, and used t < −1 s,
  i.e. **τ − t_S ∈ [−4, −2) s**. Same position relative to the (re-picked) S here. (It contains P coda at these
  stations, S−P 4.4–5.3 s at PGC/B009–B011; the original had the same content.) No other noise column is tried
  for the headline.
- "Rank chosen" is reported as (i) the number of singular values with s_k > nmf (w_k > 0.5) and (ii) Σ w_k.
- Rows are scaled back by their peak after filtering (stretch is scale-free; this only keeps the P-weighted reference
  definition identical to the pilot).
- Reference for day i = Σ_{j≠i} sumP_j R_f,j / Σ_{j≠i} sumP_j (the pilot's leave-day-out epoch reference, built from the
  filtered rows). The filtered row i is stretched against it in τ − t_S ∈ [2, 4] and [3, 5] s.

## 2. Per-day errors
- **Primary (σ)**: 200 bootstrap draws per day: resample day i's repeats with replacement (P weights kept), replace
  row i, re-run the full SVD-Wiener (σ_n and nmf re-estimated) on the whole matrix, rebuild the leave-day-out reference
  from the new filtered matrix, measure row i. σ = SD over draws (railed draws included, as in the pilot). The raw σ is
  recomputed in the same code with the filter switched off (should match the pilot σ up to the random draws).
- **Secondary (σ_joint)**: 200 draws in which every day is resampled at once, one SVD per draw.
- **Independent check (half-split)**: 20 random splits of every day's repeats into two disjoint halves → two matrices,
  each filtered and measured on its own. σ_half,emp = 1.4826·MAD(ε_A − ε_B)/√2 over all days of the series; predicted
  = √2 × median σ (a half stack has half the repeats). Ratio emp/pred reported for raw and SVD. A ratio > 1.5 for SVD
  means the bootstrap understates the SVD error, and the σ comparison in §3 is then not trusted for that series.

## 3. Raw vs SVD comparison (columns, per series and pooled)
- median σ (2–4 and 3–5 s), ratio SVD/raw;
- median cc of the day stack with its reference;
- lapse test (pilot gate d): Theil–Sen slope of ε(3–5) on ε(2–4) ≥ 0.6 AND sign agreement ≥ 70 %;
- fraction of days at the grid edge (rail).
- "SVD shrinks the error" is stated only for series where σ_SVD/σ_raw < 0.8 AND the §2 half-split ratio ≤ 1.5 AND the
  §4 injection test passes.

## 4. Leakage / recovery test (mandatory, pre-registered)
- Step dates: the three cluster (episode) boundaries of the series closest to day index n/3, n/2, 2n/3 (distinct;
  "between two episodes"); series with < 6 days or < 3 boundaries use what exists.
- Amplitudes a = +0.1 % and +0.3 %: every raw day stack on or after the step date is stretched by a about S (cubic
  interpolation, as the pilot). The SVD filter is re-run on the injected matrix; references rebuilt from it.
- **Recovery** (pilot definition) = [(wmean_post − wmean_pre)_inj − (same)_natural] / a, weights 1/σ_raw² of
  natural non-railed days. **Pass 0.8–1.2.**
- **Leakage, fixed reference** (the SVD-specific test): every filtered injected pre-step row is measured against the
  natural filtered leave-day-out reference of that day; leak = wmean_pre(ε_inj − ε_nat)/a. For the raw path this is
  0 by construction. **Pass |leak| < 0.2.** Post-step recovery against the same fixed reference reported as a column.
- **Leakage, pilot contrast** = [(last pre-step episode − earlier pre-step days)_inj − (same)_natural]/a. **Pass |·| < 0.2.**
- A series passes §4 when all three pass for both amplitudes, all step dates and both windows. Run for raw and SVD.
- Columns only (no pass rule, read with the same 0.8–1.2 yardstick): a **transient** of +0.3 % injected (i) into every
  day of the one episode containing day n/2 only, (ii) into day n/2 only. Recovery = [(that episode/day − the rest)_inj
  − (same)_natural]/a. A filter that keeps only patterns shared by many rows may pass the long step and still damp
  short changes; this column says whether it does.

## 5. Scrambled control (noise) through the same path
- Control day stack = the pilot's c2 stack (same-day c2 windows, equal weights, ≤ the day's number of data windows,
  same seed rule), on days with ≥ 3 c2 windows.
- (own) the control day stacks form their own matrix (≥ 3 control days), filtered by the identical SVD-Wiener, each
  measured against the filtered data leave-day-out reference of that day; σ by the §2 per-day bootstrap with SVD re-run.
- (ins) the control stack replaces that day's data row in the data matrix, SVD-Wiener re-run, the filtered row measured
  against the leave-day-out reference of that filtered matrix; σ the same way (100 draws).
- Pass (pilot gate b, control part): ≥ 90 % of non-railed control values have |ε/σ| < 2 AND |mean ε/σ| ≤ 0.5
  (series with ≥ 3 control values). Also reported: median cc of control vs reference, raw vs SVD. If SVD raises the
  control cc toward the data cc, the filter is shaping noise like the reference; that is written down.

## 6. Figures
- `figures/s10_svd_stacks_fam{F}.png`: packed day-stack image (as `s10_daily_stacks_fam10.png`) raw | SVD side by side
  for B011 and PGC. Raw here = the realigned day stacks that go into the SVD (the original figure had no realignment).
- `figures/s10_svd_daily_fam{F}.png`: daily dv/v (= −ε) vs calendar time with ± σ, raw and SVD, ETS shading,
  epoch boundaries.

## 7. Reading
Provisional pilot. A smaller σ alone is not a gain: it counts only if the injected steps are still recovered, nothing
leaks, the half-split scatter agrees with σ, and the scrambled control stays at zero.

## Run record (added after the run, ~9:45 am PT)
- Run 9:32–9:40 am PT, 8 workers, 21 series, no errors. Raw path reproduces the pilot's daily ε and cc exactly
  (max difference 1e-15 on the test series), so any raw-vs-SVD difference comes from the filter alone.
- Outputs (`tables/`): `daily_svd.parquet` (per day: raw_/svd_ ε, cc, σ, σ_joint, both windows), `series_svd.csv`
  (rank, Σw, half-split, lapse), `compare_raw_vs_svd.csv`, `injection_svd.csv` + `injection_summary.csv`,
  `control_svd.csv` + `control_summary.csv`, `summary_svd.json`, `stacks_svd.npz` (raw and filtered day matrices),
  `ablation_noise_cols.csv`. Figures `figures/s10_svd_stacks_fam{10,17,87}.png`, `figures/s10_svd_daily_fam{10,17,87}.png`.
- Post-hoc additions (not in the rules above, labelled as such): (i) the column σ/recovery ("error in units of a true
  change"); (ii) `s10_svd_ablation.py`: the pre-noise columns replaced by the pre-P noise window, to see whether the
  shrinkage hangs on that one choice. Fam 17 at B926 has no pre-P window inside the stored span (S−P 13.1 s), so 20 series.
- Lapse sign agreement here is counted over days where both windows give a value. The pilot's gates.csv counted a
  railed day as a disagreement, so the raw lapse numbers here are less strict than the pilot's (7/21 pass here vs
  5/48 there). Raw and SVD use the same definition, so the comparison holds.
- Verdicts: see `../RESULTS.md`, section "SVD".
