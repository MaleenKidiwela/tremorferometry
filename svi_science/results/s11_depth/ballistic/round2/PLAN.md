# S11 ballistic timing, round 2: four adaptations (rules fixed first)

Written Mon 2026-10-12, 3:35 pm PT, before any round-2 number was computed. Round 1 (`../PLAN.md`, `../RESULTS.md`) found no
path signal, a co-located S floor of 10 ms SD (3.5–4.1 ms robust), station-timing glitches at PGC (2019-11, 2021-01) and B010
(2009-03, 2022-03, 2024-08), far-pair scatter 18–21 ms, and P cc 0.24–0.42. User (10-12): "science requires perseverance,
keep trying until something works". Four adaptations, run in order, each judged by its own rule below. "Scrambled control"
= the S10 c2 windows (each station's data shifted to another time of the same day), as in round 1.

Scripts `svi_science/scripts/s11_ball2_*.py`; outputs here (`tables/`, `figures/`, `logs/`). CPU only, ≤ 8 workers,
RAM < 40 GB; nothing deleted; the WA-NW S11 queue and REDPy pilot files are not touched.

## 0. One new measurement pass (`s11_ball2_measure.py`)
Same data, windows, families (76), stations (7 core), S pick, S window, ±100 ms search, parabolic refinement, 200 bootstrap
draws and leave-bin-out references as round 1 (`s11_ball_measure.py`, functions reused unchanged). New:
- **Bins:** `EP` = episode unit (n ≥ 20 rows, as round 1); `YR` = calendar year, all data rows of the series in that year
  (also rows on non-active days, which round 1 used only in grand stacks), n ≥ 20; `Y2` = two-year bins (2003–04, 2005–06, …),
  n ≥ 20. Reference = the station-epoch grand stack minus the bin's rows. An epoch is never crossed.
- **P estimators** (window = full P wavetrain [tP − 0.5, min(tP + 2.5, tS − 1.0)] s, tP as round 1 Amendment 1):
  - `P3` (**headline**): projection on the 3C principal direction (largest eigenvector of the 3×3 Z/1/2 covariance of the
    station-epoch grand stack in the P window);
  - `P3j`: joint 3C cc (Z, 1, 2 summed, no projection, as S is measured on 1+2);
  - `Pz`: Z only; `Pzr`: round-1 estimator (Z + polarised horizontal).
  Diagnostics per series: split-half cc of the grand stack (odd vs even rows) for S and each P estimator; angle between the
  P3 direction and the vertical (incidence proxy).
- **Pairs on common rows** (exact row match, common bootstrap draws), for EP, YR and Y2: all 6 quartet pairs (PGC, B009, B010,
  B011) and each far station (B926, SNB, VGZ) against each quartet instrument (12 pairs). S and P3 differences.
- **Scrambled control** (c2) for station bins and pairs, n matched to the data, equal weights, as round 1.

## 1. Clock and timing-glitch correction (`s11_ball2_clock.py`)
Model, per (family f, bin b): station-level lag = O(f,b) [detector frame / origin, common to all stations] + K(s,b)
[station term] + noise. On common-row pairs O cancels exactly.
- **Instrument term** K_inst(s,b) for the quartet: per (f, b) solve the 6 quartet pair differences for 4 instrument lags (LS),
  gauge = **median of the instruments present** ("majority clock"; mean if only 2). K_inst(s,b) = σ-weighted mean over
  families; **for correcting family f, K_inst is estimated without family f** (leave-family-out), so the corrected floor is
  not fitted to itself. Requires ≥ 3 other families in the bin, else no correction (value kept uncorrected, flagged).
- **Site term** (common to all instruments at a site, i.e. path near the site, or a clock error shared by all four): not
  identifiable from the quartet; it is what remains in far-pair differences after the quartet instrument terms are removed.
- **Far-station term** K(X,b) − K(PGC site,b): clock and near-station path cannot be separated at single-instrument sites;
  reported as "clock + near-site path" and labelled clock only if the metadata flag it.
- **Detector-frame drift** O(b) = σ-weighted mean over families and stations of station-level all-row lags (common to all
  stations); reported as a series, no rule.
- **Glitch call:** |K_inst(s,b)| ≥ 3σ (σ from the family scatter / √n_fam) with sign agreement in ≥ 70 % of families →
  "instrument timing glitch".
- **Metadata check:** EarthScope MUSTANG daily metrics (`timing_quality`, `clock_locked`, `suspect_time_tag`,
  `timing_correction`, `num_gaps`) for the 7 stations where available, compared between glitch bins and the rest; the
  channel band used per row (cut-store index `cha`) for PGC/VGZ/SNB (an HH/BH/EH switch inside an epoch would shift timing).
- **Corrected differences:** quartet pairs dS − [K_inst(B) − K_inst(A)]; far pairs X−A: dS + K_inst(A); and a **site-mean
  reference**: mean over available A of (dS(X−A) + K_inst(A)).
- **Pass rule 1:** (a) robust co-located floor (1.4826·MAD with the round-1 dof factor, demeaned per series, round-1 pair
  sets: PGC vs B009/B010/B011 and B011 vs B009/B010) **≤ 3 ms after correction, for both reference sets**, EP bins; and
  (b) far-pair scatter (dof-corrected SD and robust) **lower after correction** than before on the same values for ≥ 2 of
  the 3 far stations, with PGC and with B011 as reference. The site-mean reference is reported as a column. Both (a) and (b)
  needed for PASS; (b) alone = partial.

## 2. Better P (`s11_ball2_measure.py` + `s11_ball2_ppass.py`)
- **Pass rule 2:** for the headline P3, the median P cc of a station's non-railed stacks ≥ 0.6 on EP, YR or Y2 stacks
  (n ≥ 20) for **≥ 3 stations**. Other estimators are columns; if only a column passes, it is reported as "passes on a
  column (4 estimators tried)". The scrambled control's median P cc on the same bins must stay below 0.3, else the cc is not
  P (withdrawn). P timing precision = median bootstrap σ_P per bin type.
- Multi-station P alignment is left for round 3 unless rule 2 fails on all estimators.

## 3. Pooled long-term constraints (`s11_ball2_longterm.py`)
Items (YR bins headline, Y2 column): S−P per station (P3, or the best estimator by split-half cc if P3 fails rule 2);
inter-station dS per far station (headline: site-mean reference after clock correction; columns PGC-ref, B011-ref,
uncorrected); dP3 per far station; co-located dS (control of the method at year scale).
- Demean per series ((fam, station, epoch) or (fam, pair, epoch pair)); series with ≥ 3 years.
- **Family-common component per year** = one-way decomposition by year as round 1 (MSB, MSW, variance split), **permutation
  control**: year labels shuffled within series, 10,000 draws, seed 20261012.
- **Decade contrast** per series: mean of years ≥ 2015 minus mean of years ≤ 2014 (within one epoch; ≥ 2 years each side);
  common contrast = inverse-variance mean (σ from bootstrap, inflated by √(χ²/dof) if > 1); permutation p from the same
  within-series year shuffle; scrambled-control contrast reported. Linear trend (ms/yr) the same way.
- **Detection** = |common contrast| ≥ 3σ AND permutation p < 0.01 AND scrambled control p ≥ 0.01 AND (for dS) the same
  sign and p < 0.05 with both PGC and B011 references. Otherwise a **2σ bound**.
- **Conversion** with observed S−P (SP) per (family, station): T_S = 1.75/0.75·SP, T_P = SP/0.75. dS on the far leg only:
  δVs/Vs = −Δ/T_S,X; uniform on both legs: −Δ/(T_S,X − T_S,ref). dP: δVp/Vp = −Δ/T_P,X (leg only). S−P: constant
  Vp/Vs: −Δ/SP; Vs only: −Δ/T_S; Vp only: +Δ/T_P. If rule 2 passes: Vp/Vs change on the far leg =
  δVp/Vp − δVs/Vs from (dP, dS).
- **Success rule 3:** any item yielding a 2σ bound on δVs/Vs or δVp/Vp tighter than 0.1 % (decade contrast), or a detection.

## 4. Source-locus nuisance (`s11_ball2_locus.py`)
Data: pilot windows (`s10_dvv/pilot/win`), families 10, 87, 17, 2003-02 to 2006-12; stations KHVB, MGCB, SHDB, SHVB, TWBB,
LZB, TWKB, TSJB, SILB, SSIB, TWGB + the 7 core stations (boreholes from 2005-11). EP stacks, n ≥ 20 (n ≥ 10 column), S lag
on 1+2 as §0, leave-episode-out reference per station-epoch, demeaned per series.
- **Model** per (family, episode) with ≥ 7 stations from ≥ 5 sites: lag_s = O + g_s·Δx − T_S,s·c + ε_s, unknowns O, Δx
  (E, N, Z) and c (common δVs/Vs), weighted by 1/σ²_s. g_s = ∂T_S/∂(source) = −(unit vector source→station)/Vs(source)
  (straight ray, Vs at the source depth from USGS CVM v1.7 south_VI); T_S,s = 1.75/0.75 × observed SP. Locations:
  `families_v1_3.csv`; stations from FDSN. Quartet K_inst corrections applied where available.
- **Wander detected** if (i) pooled over all family-episodes, the χ² of the horizontal Δx about zero (formal covariance) has
  p < 0.01, AND (ii) a geometry-scramble control (g vectors permuted among stations, 1,000 draws) gives that χ² in < 1 % of
  draws. Centroid motion = RMS horizontal Δx with the formal noise removed in quadrature (m); otherwise a 2σ bound (m).
- Reported: c with and without the locus terms (does the nuisance absorb a velocity term?), depth (expected to trade with O).

## 5. Figures (user rules: points coloured by cc, no error bars, line graphs)
`fig_clock_terms.png` (K_inst per instrument and far-station terms vs time), `fig_floor_before_after.png`, `fig_P_cc.png`
(median cc per station vs bin type, data and control), `fig_longterm_years.png` (year series, grey per family, black =
common), `fig_locus.png` (Δx per episode).

## 6. Order and conduct
Measure (nohup queue, DONE marker) → 1 → 2 → 3 → 4. A failed rule leads to the next adaptation, not to a stop; each
failure is written up with the next thing to try. Amendments are added below with time stamps, before the affected number
is seen.

## Amendment 1 (Mon 2026-10-12, 3:41 pm PT — after a family-10-only test of the measurement; written just after the 76-family run was launched at 3:40 pm, none of its output seen; no clock, floor, long-term or locus number seen)
The test (`tables/*_test.*`) showed only per-station P cc (family 10: P3 median 0.72–0.86 at PGC, B010, B011, SNB on
EP/YR stacks; 0.15–0.45 at B926 and VGZ; control 0.24) and the +5 ms self-check (S 4.91, P3 4.94 ms). No rule changed.
Noted: B009's P3 direction is 86° from vertical (its Z channel carries little P; P3 then rides on the horizontals) — kept,
reported as a station-quality column. Bins that hold all but < 10 of a series' rows are skipped (their reference would be
too thin); the round-1 rule had no such case because episodes are small.

## Amendment 2 (Mon 2026-10-12, 3:44 pm PT — while the 76-family measurement runs; no output of it seen)
§3 headline reference for the long-term dS changed, because the PGC epoch boundaries (2015-01-12, 2017-08-10) demean the
PGC pair series inside each epoch, so the PGC pairs carry no decade contrast and would dilute a site mean that includes
them:
- **Long-term dS headline = borehole-site mean**: mean over A ∈ {B009, B010, B011} (one epoch, 2005-11 onward) of the
  clock-corrected X−A values. Columns: B011 alone, PGC alone (year decomposition only), site mean incl. PGC.
- The detection clause "same sign and p < 0.05 with both references" now means: with **B011 alone** and with the
  **B009+B010 mean** (B011 excluded), i.e. two independent instruments.
- dP3 pairs get the same clock correction as dS (an instrument clock error shifts P and S alike; K from S).
- Column "split contrast": per series, the later half of its years minus the earlier half (split at the series' median
  year), so that series inside one epoch (PGC E1, VGZ, SNB) also contribute; the decade contrast stays the headline.

## Amendment 3 (Mon 2026-10-12, 3:47 pm PT — before any locus number; locus measurement not yet run)
§4: T_S,s for the velocity term is the straight-ray S time through the USGS CVM v1.7 south_VI Vs profile (family
location → station), because P picks at the POLARIS stations are untested; the observed S−P pick is kept as a column.
The 2σ bound on centroid motion = √(motion² + 2·formal²/√N_fits) (formal = RMS formal horizontal error, N_fits the number
of family-episodes; the √N term is the 2σ sampling spread of a χ²-based variance estimate).

## Amendment 4 (Mon 2026-10-12, 3:57 pm PT — after rule 1 was judged on the headline S estimator: (a) FAILED, robust floor 4.8 / 3.8 ms after correction; written before any variant below is computed)
Next adaptation for rule 1(a) (the floor is set by per-value S noise, not by clocks — see RESULTS): two pre-declared S
estimator variants, measured on the quartet pairs, EP bins only, otherwise identical (`s11_ball2_measure.py` with
`S11B2_SVAR`): `s3c25` = joint 3C (Z, 1, 2) cc in [tS − 0.5, tS + 2.5] s (more S energy and coda); `sh07` = 1+2 in
[tS − 0.3, tS + 0.7] s (direct S only). Same clock correction and floor definitions. Rule 1(a) on a variant counts as
"passes on a column (3 S estimators tried)". The headline long-term results keep the round-1 S estimator.

## Amendment 5 (Mon 2026-10-12, 3:59 pm PT — after the Amendment-4 variants: s3c25 robust floor after correction 3.59 / 3.11 ms, sh07 5.50 / 4.45 ms, both fail ≤ 3 ms; written before the variant below is computed)
The longer 3C window helped (4.8/3.8 → 3.6/3.1 ms). One more pre-declared variant in the same direction, for the
CLOCK/FLOOR purpose only (it reaches into the early coda, where a velocity change adds lag ∝ lapse time, so it is not
used for path inferences): `s3c40` = joint 3C in [tS − 0.5, tS + 4.0] s. Same rule 1(a), "passes on a column (4 S
estimators tried)" if it passes.
