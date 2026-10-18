# S15 — VLFE science, Q1 (timing vs the local LFE sequence) and Q2 (tides) — pre-registration

Written 2026-10-17 ~2:45 pm PDT, **before any Q1/Q2 number was computed** (only file layouts and column names were read).
Source of the tests: `MERLIN_SCIENCE_PLAN.md` (Merlin, 2026-10-17 ~6:30 am PT), questions 1 and 2, copied below verbatim
where they define a test or a rule. Operational choices that Merlin left open are listed in §3 and fixed here. Any later
change is appended in §6 with time and reason; nothing above §6 is edited after numbers exist.

Conventions: times in PT (PDT = UTC−7, PST = UTC−8); every criterion is a column, no hard gates; "random-draw control" /
"scrambled control", never the other word; dv/v-style rules for plots do not apply here, but raw plots (Δt histograms vs
control, phase histograms) come before any verdict. Resources: ≤ 4 workers, < 40 GB, no downloads.

## 0. Coordinator / user instruction received 2026-10-17 ~2:40 pm PT (before any number)
"Use the closed VLFE catalogue as it is. Do not re-audit, re-vet or re-flag its events. State its limits once, in the
scope section, as given: per-episode completeness and ±16 km locations. Merlin's detection-recovery check in Q1 stays,
but only as a quick column, not an abandon gate. Focus the write-up on what the timing and tides show."
→ Merlin's abandon rule "injection recovery differs by > 30 % between days 0–2 and later days" is computed and reported
as a column; it does not stop the study. The other two report-back conditions (variant sign disagreement, one front
direction only) are reported as columns/caveats with the verdict.

## 1. Data (fixed)
- **VLFEs: the 540 ETS events of the closed catalogue** = rows of `s13_vlfe/catalogue/tables/<run>_ets_events_v2.csv`
  for runs prod2010, prod2011, prod2012, prodq2013, prod2014, prod2015, prodq2016, prodq2017, prodrq2018, prodf2019,
  prodf2020, prodfy2021, prodfy2022, prodfy2023, prodfy2024, prodfy2025, prodfy2026 (11+41+60+20+1+4+4+11+85+38+28+5+2+
  63+138+29+0 = 540, the findings-page table). Episode label from the same run's `<run>_events_v2.csv` (join on `epoch`);
  episode window `start`–`end` from `<run>_episodes_v2.csv` (an episode split over two year frames = union of its windows).
  Columns used: `time_pt`, `epoch` (UTC s), `lat_best`, `lon_best`, `n_lfe_10m`, `lfe_status`.
- **LFEs** = exactly the pool `cat_e2.lfe_status()` counts (`hg_r3.lfe_lists()`):
  SVI = `s13_common.svi_rows()` (forward `s0/analysis_table.parquet` run prod & keep_xfam8; backward
  `data/catalogue/svi_catalogue_v1_4_backward.parquet` primary & keep_xfam8), series `data`, counted_once, **P ≥ 0.75**;
  WA = `lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet`, `selected_v2` & region NW, with family ids.
  Positions: SVI `s0/families.csv` (lat, lon, `along_strike_km`, + = NW); WA NW families from `data/catalogue/families_v1_3.csv`
  via `family_uid` (lat/lon, else pos_/approx_ position), projected onto the S0 along-strike axis by a linear fit of
  `along_strike_km` on local x/y of the SVI families.
- Check before use: recomputed ±600 s counts must reproduce `n_lfe_10m` (expect 431 events with ≥ 1, 109 without).

## 2. Tests copied from Merlin (verbatim where quoted)
### Q1 — "Do VLFEs fire at the local LFE front, behind it, or in re-activations?"
- **Day scale (single best first test).** "for each VLFE with ≥ 1 P ≥ 0.75 LFE within ±10 min (431 events), Δt = VLFE time −
  first LFE of those families in that episode (count-weighted mean; min/median as columns). Control: 1,000 draws of a random
  LFE of the same families and episode, same Δt."
  - **Rule (verbatim):** "PASS-late if the observed median Δt (count-weighted, concurrent families) exceeds the 99th
    percentile of the 1,000 null medians, the sign holds in ≥ 6 of the 9 episodes with ≥ 10 VLFEs (2011-08, 2012-08,
    2013-09, 2018-06, 2019-08, 2020-11, 2023-04, 2024-08, 2025-09b1) and separately in SE-ward and NW-ward groups;
    PASS-early symmetric at the 1st percentile; NULL otherwise." ("null medians" = the random-draw control medians.)
  - Report-back conditions (columns here, §0): "injection recovery differs by > 30 % between days 0–2 and later days of the
    local sequence; concurrent-LFE and VLFE-location variants disagree in sign; result holds in one front direction only."
  - Outcomes: "late = VLFEs are a post-peak phenomenon of the slipped zone (same regime as rising tidal sensitivity and
    RTRs) — new; early = VLFEs mark front arrival — new; indistinguishable = sum-of-LFEs holds at the day scale,
    minutes/hours decide."
- **Minutes scale.** "LFE counts in [−10,−2], [−2,0], [0,+2], [+2,+10] min, null symmetric."
- **Hours scale.** "fraction of 2010–13 VLFEs inside the 30 Bostock-density RTR windows of 2010-08/2011-08/2012-08/2013-09
  (`q8_rtr/tables/rtr_list_frozen_d2.csv`: start, end, x_start, x_end; within ±20 km along strike) vs the LFE-drawn null";
  "plus S3 day-scale returns (`s3/bursts.csv`, `return_streaks.csv`; 2023-04 Apr 18–19, 2018-06 Jul 8–10)".
- **The 109 events without LFEs.** "require the same sign in SE-ward (2015-12, 2018-06, 2019-11, 2021-01, 2022-03) and
  NW-ward (2019-08, 2020-11, 2023-04, 2024-08) episodes (s3 `episodes_kinematics.csv` direction; 2010–13 from
  q6_front_arrest/episodes_pnsn.csv)."
- **Detection control.** "injection recovery vs injection day in the episode from data/<run>/inj_list.parquet (epoch) ×
  inj_v2.parquet (battery); MF thresholds are per template per run, so no daily self-suppression is expected — verify."
  Column only (§0).
### Q2 — "Are VLFEs tidally modulated, and more or less than the LFEs around them?"
- "M2/K1 phase of 540 VLFEs, Rayleigh + modulation amplitude; control = LFE-time draws from the same families/episode
  matched on days-since-local-onset". Stresses: `s5_setup/tidal_stress_svi_hourly.parquet` (2014–26) and
  `s5/tidal_stress_bostock_hourly.parquet` (2003–2014-01); `tau_updip_Pa`, `sigma_n_Pa`, `coulomb_0p4`.
- "Power: at ±35 % modulation N = 540 gives Rayleigh z ≈ 16; a 0.2 difference in modulation is ~3σ."
- "Limits: 2016–17 (LRIV fault) flagged." Outcomes: "same = VLFEs are LFE-like slip; stronger = large slow events more
  tide-sensitive (friction/fluid statement); none = distinct, loosely coupled process."

## 3. Operational choices fixed now (Merlin left these open)
### Q1
1. **Concurrent families** of a VLFE = families (SVI or WA-NW) with ≥ 1 pool LFE within ±600 s of `epoch`; weight w_f = that
   count. **First LFE** of f = earliest pool LFE of f inside the episode window. Δt = t_VLFE − Σ w_f·first_f / Σ w_f (days).
   Columns: Δt_min = t − latest first_f (smallest per-family Δt), Δt_median = median over families of (t − first_f),
   Δt_earliest = t − earliest first_f.
2. **Control draw**: for each VLFE, one LFE drawn uniformly from all pool LFEs of its concurrent families inside the episode
   window; same Δt formula with the VLFE's own weights. 1,000 draws (seed 15). Statistic = median over the 431 events.
   Percentile of the observed median among the 1,000 control medians.
3. **Per-episode sign** = sign of (observed median Δt − median over draws of the control episode median), for the 9 listed
   episodes. "Sign holds" = same sign as the pooled (obs − control) difference.
4. **Direction groups**: SE-ward / NW-ward from s3 `episodes_kinematics.csv` (SVI rows). For episodes not in that table
   (2010-08b1, 2011-08, 2012-08, 2013-09, 2025-09b1, …): sign of the LFE front speed from our own pool (Spearman of family
   first-LFE time vs along-strike position, + = NW-ward), cross-checked against s3 where both exist; q6 PNSN start position
   as a column. Group test = pooled median (obs − control) sign and the group's percentile against its own 1,000 draws.
5. **Extra columns (not part of the rule)**: per-event rank u = fraction of the event's control pool earlier than t_VLFE
   (uniform under "VLFE = random LFE"); Δt with first_f replaced by the 5th-percentile LFE time of f (robust onset);
   flag when first_f falls on the first day of the episode window (window truncation).
6. **VLFE-location variant** (all 540, and the only one for the 109): families within 20 km of (lat_best, lon_best) with
   ≥ 1 pool LFE in the episode (30 km if none; else no value); weights = their episode LFE counts; same Δt and control.
7. **Minutes scale**: pool LFEs of any family within ±10 min of `epoch`; counts per bin and per minute in [−10,−2), [−2,0),
   [0,+2), [+2,+10]. Symmetry tests: pooled after/before ratio for the inner pair and outer pair, binomial p (expected 0.5).
   Random-draw control: the same bins around randomly drawn LFEs of the concurrent families/episode (the drawn LFE itself
   excluded), 1,000 draws → expected after/before ratio given LFE clustering. Raw plot: ±10 min histogram in 20-s bins
   (shows the VLFE-time offset relative to LFEs; `epoch` is the matched-filter time, not an origin time).
8. **Hours scale**: 2010–13 VLFEs (the four 2010–2013 runs). Inside = UTC time in [start, end] of one of the 30 RTRs and
   along-strike x within [min(x_start, x_end) − 20, max + 20] km. VLFE x = count-weighted along-strike of its concurrent
   families (VLFE-location x for events without LFEs; VLFE-location x for all as a variant). Control: the drawn LFE's time
   and its family's x, 1,000 draws. S3 returns: VLFE PT date = a return day of an SVI `bursts.csv` row and |x − bin_center|
   ≤ 10 km; same control. Report fractions and control percentiles.
9. **Detection column**: C0/C1 injections (Mw 3.6 and 3.8, ETS period) of the same 17 runs; recovered = mf_best_z ≥ 5.5 &
   mf_second_z ≥ 5.5 & battery_e2 after `cat_common.rescore_inj` (f′, the catalogue's completeness definition). Days since
   local onset = injection epoch − first pool LFE of families within 20 km of the source in that episode. Recovery in
   days 0–2 vs ≥ 2 (and < 0, i.e. before local onset); relative difference column ("> 30 %" flag).
### Q2
10. Phases with the S5 machinery (`svi_science/scripts/s5_lib.py`: harmonic fit per family per era, θ = V + u + arg(hc),
    θ = 0 at the constituent's peak thrust-promoting stress, −90° at peak stress rate). **Headline component
    `tau_updip_Pa`** (Royer's up-dip shear); `tau_Pa` and `coulomb_0p4` as columns; constituents **M2** and **K1** (K1 is
    locked to the solar day: flagged as confounded with any time-of-day detection effect, as in S5). 2010–2013 events use
    the Bostock-era file, 2014+ the SVI file (separate fits, as S5).
11. Stress family of a VLFE = SVI family nearest to its count-weighted concurrent-family position (VLFE location if no
    LFEs); of a control LFE = its own family if SVI, else the SVI family nearest to its family position.
12. Statistics: R = |mean e^{iθ}|, Rayleigh z = N R², p (S5 formula), modulation amplitude a = 2R (rate ∝ 1 + a cos(θ − θ0)),
    θ0 = arg mean e^{iθ}; C = mean cos θ. Stress-level columns: mean z-scored τ_updip and Coulomb at event times.
13. **Matched LFE control**: for each VLFE, one LFE drawn from the same families (concurrent; location-based for the 109)
    within ±24 h of the VLFE time (= matched on days since local onset to ±1 d, since the onset is shared; ±24 h spans two
    K1 and four M2 cycles); widen to ±48 h, then the whole episode, if empty (flag column). 1,000 draws → distributions of
    a, θ0, C. **Random-time control**: uniform random times within ±24 h of each VLFE, 1,000 draws (absolute significance
    of R given the sampling).
14. **Verdict rules (Q2)**: VLFE modulation *detected* if M2 Rayleigh p < 0.01 and R above the 99th percentile of the
    random-time control. Relative to LFEs: *stronger* if a_VLFE > 99th percentile of the matched-LFE a distribution,
    *weaker* if < 1st percentile, *same* otherwise (when detected). If not detected: *none* if a_VLFE < 1st percentile of
    the matched-LFE distribution, otherwise *not resolved at this N*. Columns: with/without 2016–17 (runs prodq2016,
    prodq2017 flagged), declustered (first VLFE of each 1-h cluster), with-LFE vs no-LFE subsets, era 2010–13 vs 2014–26.

## 4. Outputs
`scripts/` (s15_common.py, s15_q1.py, s15_q2.py, s15_detect.py, s15_figs.py), `tables/` (per-event and summary CSV/JSON),
`figures/` (raw first: Δt and u histograms vs control, minute histogram, phase histograms vs control), `RESULTS.md`.

## 5. Scope (stated once)
The catalogue is used as closed. Its limits: per-episode completeness (C0/C1 injections; 0.14–0.89 at Mw 3.8 by episode)
and ±16 km locations (11.5 km NW bias). Nothing in this study re-vets or re-labels catalogue events.

## 6. Amendments (append only, with time and reason)

### Amendment 1 — 2026-10-17 ~2:50 pm PDT (after the first Q1 day-scale run; before any re-anchored number)
Reason: auditing the first run, the fixed-weight control of §3.2 is **biased toward "late" under the very hypothesis it
tests**. A VLFE's concurrent families are, by construction, all active at the VLFE time, so the VLFE time lies after the
*latest* of their onsets, while the control draws from the *union* of their LFEs over the whole episode, including times
before some of them had started. So even if VLFEs were random LFE times, Δt_obs would exceed the control whenever ≥ 2
families are concurrent. The same conditioning (an LFE within ±10 min is required) also affects the RTR control.
Change: the pre-registered fixed-weight control stays as reported (run A), and a **re-anchored control** is added and is
the one the verdict is read from: each control draw is a pseudo-VLFE placed at a random LFE time t* from the same pool
(concurrent families, same episode) that has ≥ 1 *other* pool LFE within ±600 s; its own concurrent set, weights,
first-LFE times and Δt are recomputed around t* exactly as for a real VLFE (anchor LFE excluded from its own counts;
including it is a column). Merlin's rule (99th/1st percentile, ≥ 6 of 9 episodes, both directions) is applied unchanged
to the re-anchored control. Column checks: events with exactly one concurrent family (no intersection effect) under the
fixed-weight control. The RTR fraction is also re-computed against re-anchored draws; the location variant (§3.6) is
re-computed with re-anchored draws (families within 20 km of the pseudo-VLFE's own family).

### Amendment 2 — 2026-10-17 ~2:55 pm PDT (after the first Q2 run; checks only, the §3.14 verdict is unchanged)
The first Q2 run gives a very strong VLFE M2 modulation (a ≈ 1.0 against ≈ 0.5 for the matched LFEs). Because the
VLFEs come in bursts and the detector works in a band where coastal long-period noise could follow the ocean tide, these
robustness columns are added before interpreting it: (i) cluster bootstrap of a_VLFE (resampling 1-h and 6-h VLFE
clusters, 2,000 resamples) and declustering at 6 h; (ii) **detection vs tidal phase**: recovery of the catalogue's C0/C1
injections (17 runs, ETS and between, Mw 3.6/3.8) in M2 τ_updip phase bins, with the fitted detection modulation;
(iii) per-episode a and θ0 for the 9 big episodes, VLFE vs matched LFEs; (iv) stress-level curve (VLFE rate per τ_updip
decile relative to random times, LFEs the same).

### Amendment 3 — 2026-10-17 ~3:15 pm PDT (Merlin audit `MERLIN_AUDIT_Q1Q2.md`, its "Amendment 2"; written before any number below)
Merlin's audit keeps Q1 and Q2 and asks for one gate test, one cross-check and two Q1 columns. Fixed now:
1. **Poisson regression gate (Q2).** Episodes: all with ≥ 5 catalogue VLFEs. Per episode: VLFE centroid = median of the
   VLFEs' concurrent-LFE anchor positions (VLFE location where no LFEs); zone families = pool families within 20 km of it;
   stress point = SVI family nearest the centroid. Unit = UTC hour inside the episode window. Columns per hour: VLFE count
   (catalogue events), zone-LFE count n_L, M2 phase θ of tidal up-dip shear at the hour midpoint, days since local onset
   (first zone LFE of the episode) in factor bins {before onset, 0–1, 1–3, 3–6, ≥ 6 d}, episode factor (offsets).
   - VLFE model: log E[N_V] = α_ep + β·log(1 + n_L) + γ_bin + c·cos θ + s·sin θ.
   - LFE model: log E[n_L] = α_ep + γ_bin + c·cos θ + s·sin θ.
   - M2 amplitude A = √(c² + s²) (log-rate amplitude, rate ∝ exp(A cos(θ − θ0))), phase θ0 = atan2(s, c).
   - Cluster bootstrap: resample episode-days (24-h blocks within episodes) with replacement, refit both models, 500 resamples.
   - **Acceptance:** A_V ≥ 2·A_L and the bootstrap 95 % interval of A_V / A_L excludes 1 → the "more tide-sensitive than
     the LFEs" claim stands with the LFE-rate term in. **Failure signature:** A_V consistent with β·A_L (bootstrap interval of
     A_V − β·A_L includes 0) and A_V / A_L interval includes 1 → VLFEs are a sample of tidally modulated slip; drop
     "more tide-sensitive". Anything else = intermediate, reported as such.
   - Columns: without the β term; zone LFEs replaced by all pool LFEs; hours with n_L ≥ 1 only; without 2016–17 runs.
2. **O1 cross-check (all events):** O1 a for VLFEs vs matched LFE draws vs random times (same machinery as §3.12–13);
   Merlin's exponential-scaling expectation VLFE ≈ 0.33, LFE ≈ 0.13, random p99 ≈ 0.21.
3. **Q2 reporting columns:** declustered pair as the headline; von Mises κ and peak/trough ratio exp(2κ); LFE baseline with
   the VLFE's own ±1 h excluded from the matched draws; VLFE a by days since local onset (< 1, 1–3, ≥ 3 d), by ±10-min LFE
   count tercile, and above/below the median Mw_scan.
4. **Q1 columns:** the NW-ward single-family subset stated explicitly; the re-anchored day-scale test repeated with
   concurrency, onsets and pools defined from P ≥ 0.5 SVI LFEs (WA-NW unchanged), reported for all events, by direction,
   and for the events that were single-family at P ≥ 0.75.
5. Interpretation wording becomes "VLFEs are a late-phase, strongly tide-gated part of local slow slip; they are not a
   proportional (or convex) sample of LFE activity" unless the gate fails.


### Amendment 4 — 2026-10-17 ~3:52 pm PDT (Q3 partition, Q5 migration rate; Merlin `MERLIN_PARTITION_MIGRATION_DESIGN.md`; pasted verbatim before any number)
Q3/Q5 use the closed catalogue (540 ETS VLFEs, `q1_events.csv`), the S15 LFE pool, S13 `scale_svi_episodes.csv`, S1 `episode_geodetic.csv`, the vbICA regional slip increment for 2010–13, and S3 onset kinematics; nothing is re-vetted. Before any partition number we re-stack S13 for 2012-08, 2018-06 and 2024-08 with LFEs within ±60 s of a catalogued VLFE removed; additivity holds if the per-LFE moment stays within ±20 % of `thrust_m0`, else the partition is LFE : all-slow : GNSS. The partition table lists, per episode, LFE self-moment (family `mag_med`), continuous (`m0_total`, zone variant), discrete (raw, ÷compl_3.4, Mw ≥ 3.5 bin-corrected, 2018 relay self-detections removed) and GNSS (Tier A `D`, Tier B M0, regional slip), raw and corrected side by side; the size test is Spearman of VLFE-per-zone-LFE vs `D` on the compl_3.4 0.46–0.79 subset (7 episodes), reported with its Spearman against compl_3.4 and `n_sta_day`, and claims only a ≥ 3× effect with |ρ| ≥ 0.79 and completeness |ρ| < 0.6. For migration rate, each anchored VLFE and each pool LFE gets a local passage speed (Theil–Sen of family onsets within ±10 km, S3 onset rule) and a daily front speed; the verdict is read from VLFE-per-LFE by speed tercile against 1,000 re-anchored draws matched to ±24 h, requiring the outer-tercile contrast ≥ 1.5× above the 99th percentile, the sign in ≥ 6 of 9 big episodes and both directions; the RTR fraction is recomputed against the ±24 h-matched draws; tidal M2 `a` by tercile and a speed × tide interaction in the S15 Poisson model are reported with a claimable difference ≥ 0.3; injection recovery by speed tercile is a column with a 30 % flag. Times PT, scrambled/random-draw controls, raw plots before verdicts, no error bars on dv/v-style plots.

Operational notes fixed now (where the design leaves a choice): episodes for the partition = the 11 with ≥ 300 LFEs and
≥ 5 VLFEs plus 2015-12, 2021-01, 2021-10, 2022-03 as near-zero-VLFE rows, as in the design. Speed terciles: edges from the
pooled pool-LFE distribution of local passage speed (folded so + = along the episode's front direction). Matched draws for
Q5 = re-anchored pseudo-VLFEs (s15_q1b eligibility) within ±24 h of each VLFE, 1,000 draws. Detection-vs-speed: column with
a 30 % flag, not a gate.

### Amendment 5 — 2026-10-17 ~5:23 pm PDT (Merlin plan items 4, 5 and the b-value substitute; written before any number)
1. **Q4 detection map (is the Strait segment the only source we can see?).** Green's-function (GF) thrust sources (scan
   mechanism 315/12/90, the scan's own synthetics `scan_common.raw_disp`/`templates_grid_v2.npz`) at 7 scan-grid nodes:
   5 along the SVI LFE cloud 48.3–48.7 N at along-strike quantiles 0.1/0.3/0.5/0.7/0.9 of the SVI families there (nodes
   112, 109, 126, 143, 176) and 2 under the Olympic families (quantiles 0.25/0.75: nodes 14, 80); plus GF references at
   the C0/C1 nodes (77, 94) so the comparison is like for like. MF channels = the run's velocity frame (zero-phase
   2nd-order 0.02–0.05 Hz band-pass of displacement, central difference), scan channels = the grid thrust templates.
   Injected through the production runs' own detector, bank and battery (`cat_pilot.inj_row`, unchanged): 2012 = prod2012
   (2011 bank), 2024 = prodfy2024 (2011 ∪ 2018 relay bank, Y12 stations), Mw 3.6 and 3.8, at the times of the run's first
   150 ETS and 150 between-ETS C0 Mw 3.6 injections (paired times across nodes). Recovered = MF ≥ 5.5 on two templates and
   battery e′ with f′ (the catalogue's completeness definition). Bank = full bank (GF sources are not templates).
   - Check first: GF reference at C0/C1 vs the empirical C0/C1 sources (waveform correlation per channel; recovery).
   - Outcomes (Merlin): recovery ≥ 0.5 at nodes away from C0/C1 with no real events there → physical confinement (then
     Slab2 geometry 48.0–48.5 N); < 0.2 → the catalogue cannot speak about other patches. Ratios to the GF C0/C1 reference
     are reported as columns. Run under marker `S15_Q4_DONE`.
   - Then per episode: along-strike centroid of VLFE locations vs of their concurrent LFEs (median, S0 frame), and the
     difference by front direction (the ~11.5 km NW location bias is the expectation).
2. **Original Q5 (column): VLFEs in the regional slip-rate history.** Regional mean `slip_rate_mm_per_day` (115 families,
   daily) from `slip_svi_family_daily.parquet`; per episode, percentile of the rate on VLFE days within the distribution
   on patch-LFE days (zone families, PLAN Amendment 3), and the median day offset VLFE − slip-rate peak. Resolution ±5 d.
3. **b-value substitute:** per episode the largest Mw_scan and the count of events ≥ Mw 3.6, divided by compl_3.6
   (span-weighted episode completeness), with Poisson 95 % ranges.
