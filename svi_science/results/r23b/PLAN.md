# R23b — findings 1 and 6 re-done on our own 23-year SVI LFE record (2003–2026): PLAN (pre-registered)

Written Sun 2026-10-11 ~5:40 am PT, **before any R23b number was computed**. Only inventories were looked at while
writing it: catalogue / rate-table columns, station sets per year in the rate tables (forward S15 2014–18, S14 2019–20,
S13 2021–26; backward CN4/ERA/B2/B3/S15), the r23 episode list (145 episodes, 52 majors, 23 big ETS), the CN channel
epochs (`backward_scope/inventory_channels.csv`), and one timing check (backward S15f detection times equal the
forward cut-store detection times for 2011-08, 1,443/1,443 to < 1 µs, so the cut-store template offsets give the
expected S time for backward rows). Changes after this point go to "Revisions" with a PT time and the reason.
Science only: no catalogue or rate file is changed. Scripts `svi_science/scripts/r23b_*.py`; outputs here (`tables/`,
`figures/`, `logs/`). CPU ≤ 8 cores, RAM < 60 GB. "Scrambled control" = the same detector on the same day with each
station site's data shifted to another time of that day. Times shown to readers are PT; joins are UTC days.
No Merlin review (Fable quota).

## 0. Usage rules (from `../backward/RESULTS_BACKWARD.md`, MUST)
- Primary detector only; single repeats: era rows P ≥ 0.90, frozen rows P ≥ 0.75; `counted_once`. Forward 2014–26:
  `s0/analysis_table.parquet`, run prod, P_same ≥ 0.75, `keep_xfam8`. These are exactly the r23 event rows
  (`r23/tables/events_2003_2026.parquet`, 115 families), reused unchanged.
- Episode sizes across networks only through the per-family Bostock-equivalent `rate` (backward primary rows +
  forward `refine/rates/famday.parquet`); `rate_kc` (completeness-chained) is a column.
- No long-term between-episode rate trend across a network change. Nothing in R23b uses between-episode rates.
- Network segments (r23): ERA 2003-02 → 2006-06 (era detector), B2f/B3f 2006-07 → 2009-08 (frozen), S15f 2009-08 →
  2013 (frozen forward model on the forward station set), FWD 2014 →. Forward station sets: S15 2014–2018, S14
  2019–20, S13 2021–26. **Like-for-like pair = S15f (2009-08-27 → 2013-12-31) vs FWD-S15 (2014–2018 days whose rate-table
  `sset` is S15)**: same model, same station set.

## 1. Shared foundation (`r23b_base.py` → `tables/`)
- Episodes = the r23 list (145; 52 majors; 23 big ETS = Σ rate split at 1,640, 10 in 2003–13, 13 in 2014–26). Event
  windows = S3/Q1 padded windows (`s3_build.padded_windows`: ±3 d, capped at half the gap) on the full 145 list.
- Patches = the 12 S0 location groups (position = group centroid, S0 Slab2 frame), as Q1.
- Per patch × episode and per family × episode:
  - `n_once` = counted-once data events in the padded window (r23 rows); `n_c1` = scrambled-control rows, same rule;
  - `R` = Σ `rate` over the window's days (a day counts if day + 12 h is inside the window), clipped at 0 per
    family-episode; `R_kc` the same with `rate_kc`;
  - **`C_eq` (headline unit) = N_once(e) × R(p,e) / Σ_p R(p,e)**: Bostock-equivalent (network-corrected) shares put
    on the counted-once scale of that episode, so residuals, χ²/df and multinomial controls keep real counting noise
    while the family detection differences between detectors (k varies ×10 between families within one set) cancel.
    Columns: raw counted-once counts (`once`), `R_kc` shares (`kc`).
- Tremor per episode (Q1 rule): epicentres ≤ 10 km of a family inside the padded window, along-strike p05 / p95 and the
  Theil–Sen speed (direction). Source by date as r23: env CN4 windows (`r23/tables/tremor_points.parquet`, src env)
  before 2010-01-06, PNSN from 2010-01-06. Column: env windows also for 2010–2013.
- Detectability column (D1) per patch × episode = log(mean day station count × fraction of observable family-days),
  minus the patch mean (rate-table `n_sta_day`, r23 `famdays_observable`).

**Acceptance (foundation):** the rebuilt counted-once counts per group × S3-big episode (2015–2024) must equal Q1's
`counts_long.parquet` (`once075`, group) exactly; tremor p05/p95 for those episodes within 1 km of Q1's
`tremor_ep.csv`. Fail → stop and report.

## 2. Finding 1 — "Each ETS slips the fault in a different pattern, set mostly by how far the front gets" (`r23b_f1.py`)
Code reused unchanged by import: `q1_tests` (Design, run_budget, kin_tests, h4_and_l2, episode_level, fe_slope,
fe_multi_r2, calibrate), `s3_patches` (pair_stats, excursion, persistence/catch-up controls), `episode_rates`
(`share_tests`, `mantel`, `set_jaccard_svi`). Shuffle counts as the originals (2,000; per-patch 1,000; calibration 100
simulations × 200 draws instead of 200 × 200, to keep run time down).

**Designs.** Targets = big ETS. History for "time since last active" = every listed episode of the same period.
Residual pool (Q1 `era`): headline `C_eq` = one pool over all targets of the design; column `once` = pools by network
segment group (ERA, B2f+B3f, S15f, FWD), order shuffles inside each group.
- R0 reproduction: 2014–26, S3's 10 big, counted-once — must give Q1's numbers (K1 in-reach β +1.37, R² 0.53; H2 −1.11,
  p 0.001; unevenness vs tremor extent ρ −1.00) and S3 rule 7 (+0.103; consecutive −0.45). Mismatch → reported.
- A (headline): `C_eq`, 23 big ETS; reported as 2003–13 (10 targets, 2003–13 history), 2014–26 (13), 2003–26 (23).
- Columns: `once` (segment pools); `kc`; targets with ≥ 300 counted-once events only; env tremor for 2010–13 too.

**Tests (Q1 / S3 rules, unchanged).** K1 in tremor extent (+) and km outside (−) vs shuffled tremor extents; K2 / K3
direction; D1; H1 days since last active (+), H2 previous big-ETS residual ≥ 90 d (−, "turn-taking"), H3 cumulative
(−) vs order shuffle; H4 budget beyond front terms; L2; partial R² table; E1 episode unevenness (χ²/df of `C_eq`) vs
tremor extent width, LFE count, front speed; S3 rule 7a persistence (mean residual correlation of pairs ≥ 90 d apart
vs patch-label shuffle and multinomial control) and 7b catch-up (consecutive pairs vs order shuffle + cumulative
excursion).

**Neighbours (figure analogue of `f1_neighbours.png`).** Family × MAJOR episode `C_eq` (52 majors), `episode_rates`
`share_tests` (40-km neighbourhoods; exposure D = days with ≥ 1 tremor epicentre/window ≤ 30 km, r23 n30; X = Σ n30)
and `mantel` (log distance + carrier Jaccard, 1,000 permutations), residual correlation in bins < 3, 3–6, 6–10,
10–20, 20–40, 40–80, > 80 km. Blocks 2003–13 (13 majors), 2014–26 (39), 2003–26 (52). Columns: counted-once counts
inside one network (S15f 2010–13 majors; FWD 2014–26 majors).

**Tremor validation (check, frozen):** on the 4 big ETS of 2010–13, env extents vs PNSN extents: median |Δp05|,
|Δp95| ≤ 15 km → env extents validated. If not, every 2003–09 K1 / E1 / K2 / K3 number is flagged `unvalidated` and
the 2003–13 verdict uses the 2010–13 targets only.

**Verdict (Finding 1), on the headline design:**
- **holds** if all of: (i) patterns differ beyond counting noise (patch χ²/df > 2 in the joint design, and family
  over-dispersion q < 0.05 for ≥ 1/4 of tested families); (ii) K1 in-reach β > 0 with p < 0.01 both in 2003–13 alone
  and in 2003–26, and in-reach is the largest single within-episode R² (≥ 0.3) in 2003–26; (iii) no lasting
  strong/weak patches: rule 7a persistence not passed (p ≥ 0.01) in 2003–13 and 2003–26.
- **changes** if any of: K1 β ≤ 0 or R² < 0.1 with p ≥ 0.05 in 2003–13 or 2003–26; persistence passes (p < 0.01) in
  2003–26.
- **undecided** otherwise.
- Turn-taking (H2) is reported by period, not part of the verdict (it was SVI-2015–26-only before).
- Neighbours: "neighbours move together" holds in a block if the distance coefficient is negative with permutation
  p < 0.01 AND the mean correlation is > 0 in both 3–6 and 6–10 km bins (S0b criterion).

## 3. Finding 6 — "LFE sizes and rates did not change over the years"
### 3a. ETS-episode sizes per family (`r23b_f6_rates.py`)
- S_fe = Σ `rate` over the episode's days (start…end, r23 `S_rate` definition) per family f and big ETS e. Columns:
  `rate_kc`; counted-once counts; Bostock's own per-family counts (his 10 windows, 2003–13, fixed method).
- **E-L4L (like-for-like):** big ETS of S15f (2010-08b1, 2011-08, 2012-08, 2013-09) vs big ETS of FWD-S15 (those of
  2014–2018 with ≥ 90 % of family-days on S15). Statistic: median over families (≥ 20 Bostock-equivalent events in each
  side) of log10(mean S_fe late / mean S_fe early); also the episode-total ratio. Exact permutation of the episode
  labels between the two sides; bootstrap over families and episodes (2,000) for the 95 % CI. Column: counted-once.
- **E-DEC (whole record):** same statistic 2014–26 vs 2003–13 (13 vs 10 big ETS), with `rate` and with `rate_kc`;
  permutation of episode labels (10,000); bootstrap CI. Episode-level: Mann–Whitney of S_e, and release per year
  (Σ S over all majors ÷ years) per decade (descriptive).
- **E-TREM (normalised):** S_e per tremor family-exposure-day — PNSN 2010–26 (S15f vs FWD-S15, S15f vs all FWD) and
  env (fixed CN4 envelope network) within 2003–13; never across the env/PNSN switch.
- **E-B (Bostock check):** within his catalogue 2003–13, per-family log-count Theil–Sen slope over his 10 windows
  (episode-order permutation); and Spearman of our per-family-episode S_fe vs his counts (validation of `rate`).
- Family shares: Spearman of family mean shares 2003–13 vs 2014–26 (descriptive).
- **Verdict (3a):** **holds (no change)** if the E-L4L CI includes 1 (no change on the same model and stations) AND
  the E-DEC CI includes 1 under `rate` or `rate_kc`; **changes** if E-L4L excludes 1 with p < 0.01 AND E-DEC excludes 1
  in the same direction under both `rate` and `rate_kc`; **undecided** otherwise. The between-episode rate is not
  trended (usage rule).

### 3b. LFE size per detection 2003–2026 within ETS (S6 ΔM; `r23b_f6_wave.py`, `r23b_f6_size.py`)
- **Station set (fixed for the whole record):** CN.PGC, CN.SNB, CN.VGZ, CN.NLLB (the only S6 stations present
  2003–2026; NLLB has no 2003 data). PO/C8 era stations and KELB are NOT used, so the KELB gain problem and the era-only
  amplitude feature (`rel_amp`) do not enter. Comparisons are on this same station set throughout ("same station set"),
  and within or across detector segments only above a shared per-family completeness floor (below).
- **Events:** data, counted once (single-repeat rule), inside the padded windows of the 52 majors (headline); column:
  big ETS only.
- **Waveforms:** 2014–26 = the S6 store (`s6/env/wave_*.npz`); 2003–13 = raw miniSEED (`lfe_extension/traces_backward`),
  cut −15…+30 s about the expected S (t + cut-store template offset of the family × station), demeaned, 2–8 Hz
  Butterworth-4 zero-phase at the native rate, divided by the channel-epoch sensitivity (NRCan inventory rows of
  `backward_scope/inventory_channels.csv`), spline to 20 Hz on −14.5…+4.5 s; a component with a gap in that span is
  dropped (the S6 rules).
  **Reader check (frozen):** on 2010-08 → 2013 events both the raw reader and the S6 cut-store reader are run; median
  |Δ log10 a| ≤ 0.02 and 90 % of |Δ| ≤ 0.05 required, else the raw reader is fixed before any test.
- **Amplitude:** S6 matched-filter ratio a = Σ x·T / Σ T·T against the S6 family × station reference stack T (2014–26
  headline; P ≥ 0.5 fallback), best lag ±0.25 s; noise = the same projection on the two pre-arrival windows. Forward
  values are S6's `mf_station.parquet` rows (identical code and stacks).
- **ΔM (fixed CN set):** best stations = the S6 rule restricted to the 4 CN stations (eligible = S6 window check passed;
  ranked by forward median S/N; ≥ 2 and present ≥ 50 % else top 3 with presence ≥ 50 %, flagged weak); ref(f,s) =
  median log10 a over the 2014–26 headline; ΔM = median over best stations present (n_used ≥ 2). Station-epoch
  correction (S6 rule, residual vs the other best stations, ≥ 20 rows per epoch) with epochs at every CN change:
  PGC 2015-01-12, 2017-08-10; SNB 2006-04-20, 2008-12-11, 2011-03-17, 2017-06-08; VGZ 2007-02-28, 2009-07-31, 2017-03-03;
  NLLB 2007-01-27, 2007-02-27, 2017-02-23. ΔN the same on the noise projections.
- **S-window check per segment:** the median normalised envelope of each segment's events (per family × station) must
  peak inside −1.5…+4.5 s (S6 rule); the share of family × stations passing is reported per segment; a segment with
  < 70 % passing is reported but not used for tests.
- **Completeness (shared floor):** z = ΔM − ΔN; Mc_z(f) from the 2014–26 headline (S6 max-curvature rule, group
  fallback); per family and segment Mc(f, seg) = Mc_z(f) + median ΔN of that family's events in the segment; the
  shared floor of a comparison = the maximum of the compared segments' Mc(f, seg). Only events with ΔM ≥ floor count;
  ≥ 3 such events per family per compared group.
- **Tests (S6 statistics; families weighted equally; bootstrap 2,000 over families; permutation 2,000):**
  - T6a whole-record trend 2003–2026: per family (≥ 20 above-floor events, span ≥ 10 yr) Theil–Sen slope of ΔM vs time
    (per decade), floor = max over all its segments; mean over families; permutation = event times shuffled among the
    family's episodes (whole episodes move).
  - T6b like-for-like: S15f majors vs FWD-S15 majors, mean over families of [median ΔM late − median ΔM early] above
    the pair's shared floor; permutation = episode labels swapped between sides.
  - T6c PGC alone, 2003-01-01 → 2015-01-11 (one instrument, one gain, no correction at all): r = log10 a_PGC − ref,
    the T6a statistic within that span.
  - Columns: segment medians (ERA, B2f+B3f, S15f, FWD) relative to FWD above each pair's shared floor; all rows (no
    floor); big-ETS-only; ΔN version of every test.
- **"Seen"** (S6 rule) = p < 0.01 AND bootstrap 95 % CI excludes 0 AND the same sign on all rows AND the ΔN effect is
  smaller than the ΔM effect. **Verdict (3b): holds (no change)** if T6a is not seen with its 95 % CI inside ±0.05 log
  units per decade AND T6b is not seen; **changes** if T6a is seen AND (T6b or T6c) is seen with the same sign;
  **undecided** otherwise.

## 4. Outputs
`RESULTS.md` (per finding: plain verdict first, then numbers, then what is new vs the earlier version), `tables/`,
`figures/r23b_f1_neighbours.png`, `r23b_f1b_why.png`, `r23b_f6_size.png`, `r23b_f6b_ets_sizes.png`; notes lines in
`svi_science/NOTES.md` and `notes/2026-09-25_Notes.md`.

## Revisions
### Revision 1 — Sun 2026-10-11 ~5:50 am PT, after a 50-draw smoke run of `r23b_f1.py fast` (POST-HOC; flagged)
- Seen in the smoke run: in 2003–13 every big ETS had PNSN/env tremor over the whole SVI box (in-reach = 1 for 119 of
  120 patch × episode cells; the only exception is 1 end patch in 2003-02). K1 then has no variation, so its slope and
  R² are ~0 by construction (smoke: β +0.08, R² 0.02). The §2 "changes" branch ("K1 β ≤ 0 or R² < 0.1 with p ≥ 0.05")
  did not foresee a period without out-of-reach patches.
- Added (post-hoc, stated as such in RESULTS): a K1 test is **not testable** in a period with < 10 out-of-reach
  patch × episode cells; the pre-registered letter verdict is still printed next to the revised reading.
- Added (post-hoc column): do 2003–13 patterns cluster by detector segment? Mean residual correlation (leave-one-out
  log residuals of `C_eq`) for target pairs in the same network segment group vs different groups, with 10,000
  segment-label permutations. A clear excess (p < 0.01) would mean network effects survive the `rate` correction.
- The R0 reproduction matched Q1/S3 exactly in the smoke run (K1 +1.372, R² 0.529; H2 −1.112; 7a +0.103; 7b −0.451).
### Revision 2 — Sun 2026-10-11 ~6:30 am PT, after the full runs of r23b_f1 / r23b_f6_rates / r23b_f6_size (POST-HOC columns; no rule changed)
- **F1 persistence (rule 7a) passes in the joint 2003–26 design (p 0.002) but in neither decade alone.** Post-hoc
  `r23b_f1_posthoc.py`: rule 7a rerun (same code) on (a) full-box big ETS only (tremor in-reach ≥ 0.9), (b) without the
  three era-detector ETS (their large counted-once N dominates the pooled share), (c) both, (d) all 23 with every ETS
  weighted equally in the pool. Columns only; the pre-registered letter verdict is printed next to the reading.
- **F6a like-for-like:** the exact permutation of 4 vs 4 big ETS has a smallest possible two-sided p of 2/70 = 0.029,
  so the pre-registered "changes" branch (p < 0.01) could never be met; this design flaw is stated, not repaired.
  Post-hoc columns (`r23b_f6_rates.py posthoc`): S15f vs FWD-S15 without the two ETS with few LFEs (2014-11:
  139 counted-once LFEs under 1,620 tremor family-days; 2017-03: B010 out all 18 days + VGZ swap, the S0 flag);
  all MAJOR episodes (not only big) of S15f vs FWD-S15 per PNSN tremor family-exposure-day.
- **F6b T6c** (PGC alone) is "seen" (+0.032 log units/decade, p 0.002) while T6a is not: post-hoc column T6c on the
  frozen-detector span only (2006-07-01 → 2015-01-11, no era detector), to test whether the era detector's extra small
  events leak above the shared floor.
