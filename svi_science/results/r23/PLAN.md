# R23 — science reruns on the 23-year SVI LFE record (2003–2026): PLAN (pre-registered)

Written Sun 2026-10-11 ~12:10 am PT, **before any R23 number was computed**. Only inventories were looked at
while writing it: catalogue columns, the date ranges of each detector/station set, the per-year day coverage of the
v1.4-backward rate table (every day 2003–2013 has observable families), the env-flag coverage per year, and the
2010–2013 PNSN episode rows of `b_ets_record/episodes_b.csv`. Changes after this point go to "Revisions" at the end with
a PT time and the reason. Science only: no catalogue file is changed. Scripts `svi_science/scripts/r23_*.py`; outputs in
this folder (`tables/`, `figures/`). CPU ≤ 8 cores, RAM < 60 GB. "Scrambled control" = the same detector on the same day
with each station site shifted to another time of that day. Times shown to the reader are PT; all joins are UTC days.

## 0. Usage rules carried from `../backward/RESULTS_BACKWARD.md` (MUST)
- Primary detector only (`primary == True`): era 2003-02-01 → 2006-06-30, frozen otherwise.
- Single repeats: era rows P ≥ 0.90 (P\*), frozen rows P ≥ 0.75; `counted_once` (data) and the same flags for the
  scrambled control (`series == c1 & keep_xfam8 & primary`, same P rule).
- Between-episode activity is a control-subtracted bracket, never a rate. No long-term between-episode trend across a
  network change (2005-07, 2005-11, 2006-07, 2007-11, 2009-08, 2009-11 LRIV, 2014-01 join).
- `burst_day` and `n_z2` are reported as columns (with / without burst days; n_z2 ≥ 2 subset); never used to drop.
- Episode sizes across networks only through `rate` (Bostock-equivalent).

## 1. Shared foundation (`r23_base.py` → `tables/`)

**F1 Families.** The 115 production families (`s0/families.csv`: group, Bostock depth, along-strike / down-dip km).
The 11 families present only in the backward catalogue are left out of every headline (their counted-once count is
reported once). Deep = P3 deep groups 1, 4, 10, 12 (centroid down-dip ≥ +10 km, 32 families); column: depth ≥ 40 km.

**F2 Events.** 2003-01-01 … 2013-12-31: `svi_catalogue_v1_4_backward.parquet`, rule §0. 2014-01-01 … 2026-09-15:
`s0/analysis_table.parquet` run prod, `keep_xfam8`, P_same ≥ 0.75 (the forward rule, as S7/P3). Data and scrambled
control rows, same steps. Network segment column `seg`: CN4f 2003-01-01…01-31 · ERA 2003-02-01…2006-06-30 · B2f
2006-07-01…2007-11-12 · B3f 2007-11-13…2009-08-26 · S15f 2009-08-27…2013-12-31 · FWD 2014-01-01 →. (S15f uses the
frozen forward model on the forward station set, so S15f vs FWD is the closest like-for-like pair.)

**F3 Observable family-days.** Backward rate table, primary rows with `m_data > 0`; forward
`s0/famday_detectability.parquet` (run prod, observable).

**F4 Tremor.** Before 2010-01-06: the env-CC tremor windows of `backward/p4/envflag/win/env_<Y>.parquet` that make the
frozen flag's tremor window (sset CN4, series data, coh ≥ 0.42, s3 ≥ 0.42, located ≤ 30 km from a family); from
2010-01-06 to 2026-05-27: PNSN epicentres (`catalogs/pnsn_tremor_cascadia_full.csv`). Per family-day: n30 (≤ 30 km
epicentral), n30_100. Env for 2010–2013 is kept as a column (overlap check). The env CN4 detector is one fixed
4-station network for all of 2003–2013 (3 stations in 2003: NLLB missing).

**F5 Episodes 2003–2026.**
- 2014-01-01 →: `episodes_svi_final.csv` (79; 39 majors). 2010-01-06 … 2013-12-31: the PNSN extension in
  `b_ets_record/episodes_b.csv` (same rule; the 2013-12-25 episode crossing 2014-01-01 is kept from that file).
- 2003-01-01 … 2010-01-05: the same rule (`episode_rates.find_episodes`: ≥ 5 % of the 115 families exposed, gaps ≤ 5 d
  bridged, span ≥ 5 d, MAJOR if family-exposure-days ≥ 241.5) with env exposure (family-day exposed if n30 ≥ 1 env
  tremor window). Variant (b): a family-day is exposed only if the day is also env-active (N_CN4 ≥ 12, the flag's k).
- **Validation (frozen):** each variant is also run on 2010-01-06 … 2013-12-31 and compared with the PNSN list.
  Recall = PNSN majors (7) that overlap an env major; precision = env majors that overlap any PNSN episode; and every
  Bostock window 2003-02, 2004-07, 2005-09, 2007-01, 2008-05, 2009-05 must overlap an env major of 2003–2009.
  PASS = recall ≥ 6/7 AND precision ≥ 0.8 AND 6/6 Bostock windows. Use (a) if it passes, else (b) if it passes, else
  the fallback: majors = the Bostock windows (start/end from his catalogue) and minors = variant (b) episodes that do
  not overlap them, flagged `unvalidated`. Start/end differences vs PNSN on 2010–13 reported.
- **Episode size** `S_rate` = Σ over the 115 families and the UTC days start…end of the primary `rate` (backward table;
  forward `refine/rates/famday.parquet`). Columns: Bostock catalogue events (his 10 windows), forward S3 counted-once
  count, S1 GNSS D (2014+), `rate_kc` sum (backward).
- **Big ETS (2003–2026)** = upper class of the optimal two-class split of log10(S_rate) over all majors 2003–2026
  (B's method, applied to the rate size). Columns: all majors; B big (2010+); S3 big (2014+).
- Onset (for intervals and clocks) = list start day; column: 5th percentile of counted-once event times in the padded
  window.

**Acceptance (foundation):** forward event counts must reproduce S7 (76 BE / 3 control, 2014-01-01…2026-05-27) and
the backward counted-once totals must reproduce RESULTS_BACKWARD (51,539 / 305 over all families; per year table).
Fail → stop and report.

## 2. Study 1 — front arrest (finding 8, `q6_front_arrest`)
- Code reused: `q6_lfe.episode_tables` and `run_check` (walk = continue/stop hazard; enter = active or not; β of
  z(log10 T_any) and z(log10 T_big) with distance, log size and bin-propensity covariates; episode-order shuffle,
  2,000 draws, seed 0). Bins 5 km of family along-strike; bin active = ≥ 3 events in the episode; history = every
  listed episode 2003–2026 (padded S3 windows); start = S3 `lfe_start10` where it exists, else the median position of
  the first 10 % of events.
- **Targets (primary): rate-big ETS 2003–2026**; columns: all majors; S3-big 10 (2014+ only; reproduction check: must
  give the original's 6 stops and β −0.35 to within the shuffle draw noise — if not, the difference is reported).
- **Shuffle groups = network segment** (start dates permuted only among targets of the same segment family: ERA+CN4f,
  B2f+B3f, S15f, FWD), so detector sensitivity moves with the slot.
- Tremor-history column: T_trem = days from the episode's first event back to the last 5-km bin-day with ≥ 3 F4
  tremor windows/epicentres (≤ 30 km of a family) before the padded window; same model and same order shuffle (the
  window moves with the permuted date; tremor inside the moved window is excluded).
- Segments: 2003–13 targets alone, 2014–26 alone, joint 2003–26 (headline).
- **Verdict (joint, walk, T_any):** at least 20 stop trials needed, else **undecided (too few stops)**. With ≥ 20:
  **supported** β > 0 and one-sided p < 0.01; **suggestive** β > 0, p < 0.05; **contradicted** β < 0, p_neg < 0.01;
  otherwise **not supported**. Descriptive: stop share for T_any < 30 d vs 30–90 d vs > 90 d.

## 3. Study 2 — between-episode deep bursts (finding 3; s7, q5_bursts, c_decay)
- BE (S7, frozen): event day outside every listed episode ± 2 d AND n30 = 0 for its family that day (F4 source of the
  date). BE_broad = outside episodes ± 2 d. Eligible family-days = F3 ∩ the same condition. 2014–26 period =
  2014-01-01…2026-05-27 (S7); 2003–13 = 2003-01-01…2013-12-31.
- Per segment (2003–13, its sub-segments ERA / B2f+B3f / S15f, 2014–26):
  1. **Exist:** data vs scrambled control counts. Exist if data ≥ 20 AND ≥ 5× control.
  2. **Bursty:** distinct days vs family-preserving day shuffle *within the event's network segment* (1000×); p < 0.01.
  3. **Deep:** mean depth vs 10,000 draws from the families' ETS shares (ETS = events inside the segment's majors);
     p < 0.01; per burst as a column.
  4. **Same patches:** per-group BE share 2003–13 vs 2014–26, Spearman over the 12 groups with a group-label
     permutation (10,000); same if ρ > 0 and p < 0.01. Top families and G10 share listed. Column: BE_broad.
  5. **Burst sizes** (events chained by gaps < 24 h; S7/Q5): events per burst, Bostock-equivalent size (Σ k of the
     family-set per event, from the rate table), duration (h), families, max family separation (km), median depth.
     Compared S15f vs FWD (same model) and ERA / B2f+B3f vs FWD: Mann–Whitney p and medians; "no detectable difference"
     if p ≥ 0.05, "different" if p < 0.01. Event counts per burst are not compared across the era/frozen switch as a
     claim (sensitivity ×8–10 on ETS days), only the Bostock-equivalent size and the shape columns.
  6. **Early-2006 bursts:** each date of `backward/first_look/burst2006/NOTE.md` (2006-01-15/16, 03-27…30, 04-24, 05-01,
     05-03, 05-11, 05-25/26, 06-18): counted-once P\* ≥ 0.90 events (data / control), families, groups, depth, env N,
     BE / BE_broad status, `burst_day`.
  7. Counts per eligible family-day within each segment are listed, never trended across segments.
- Columns everywhere: without `burst_day` rows, n_z2 ≥ 2 rows, env tremor source for 2010–13.
- **Verdict:** "exist in 2003–13 … same patches … same sizes" each answered supported / not supported / undecided by
  the rules above.

## 4. Study 3 — deep quiet after big ETS (finding 4, provisional; s7, p3, q4)
2003–2013 is new, independent data: this is the **confirmatory** test of the 1–2-month quiet.
- (A) **S7 quiet month:** BE events, τ = days since the previous major's end (D1; P3's D2 = end + 7 d and D3 = the
  family's own last tremor day as columns), F30 = share with τ < 30 d vs the family-preserving, segment-preserving day
  shuffle (1000×); p_low. Clock = all majors (S7); column = big ETS. BE_broad as a column.
- (B) **P3 recharge:** deep-group events outside the reference episodes, daily counts per observable deep family-day;
  D30 = 1 − rate(0–30 d) / rate(30 d … next start); time-shift surrogates within each cycle (1000×; a cycle that
  crosses a network-segment boundary is cut there and each part shifted separately); fit r∞(1 − a e^(−t/τ)), bootstrap
  over cycles for τ and a. Clock = **big ETS (primary)**; column = all majors. Day-count version (days with ≥ 1 deep
  event) as a column (P3 post-hoc).
- (C) Detection-free column: env CN4 deep tremor (windows ≤ 15 km of a deep family), 2003–2013, same D30 and surrogate.
- **Verdict (2003–13):** **confirmed** if (A) p_low < 0.01 and (B) D30 > 99th percentile of the surrogates
  (p < 0.01) with the bootstrap CI of a excluding 0; **partly** if one passes; **rejected** if neither passes AND both
  point estimates are on the wrong side (F30 ≥ shuffle median, D30 ≤ 0); **undecided** otherwise. Number of big ETS
  and cycles in 2003–13 stated. 2014–26 recomputed with the same code (must match S7/P3 numbers); joint 2003–26 given.

## 5. Study 4 — cycle regularity and "bigger episode → longer wait" (b_ets_record; GNSS ρ 0.80)
- Big ETS 2003–2026 (F5); intervals between consecutive onsets (the 2005-09 → 2007-01 gap is in by construction).
  Columns: all majors; 2010+ B big list for comparison with B.
- **Regularity (B Q1):** CV of intervals vs the uniform-onset scrambled control (same number of onsets over the same
  span, 10,000); quasi-periodic if CV < 1st percentile. Drift: Spearman(interval, mid-year), permutation p.
- **Time-predictable (B Q2 iii):** Spearman(S_rate_i, interval after i), permutation p (10,000); detrended column (rank
  residual after a linear trend in onset year). Slip-predictable (size vs interval before) alongside. Columns: Bostock
  counts (2003–13), GNSS D (2014+), S_rate restricted to the 115 families' quiet-independent ETS days.
- **Verdict:** **supported** p < 0.01 (both raw and detrended same sign); **suggestive** p < 0.05; else **not
  supported**. 2003–13, 2014–26 and 2003–26 reported.

## 6. Study 5 — pre-ETS ramp-up
- Events: deep-group counted-once events outside every listed episode ± 2 d (BE_broad, deep groups). Columns: BE
  (tremor-free) all families as S7; scrambled-control rows through the same steps.
- τ_next = days to the next major start. **Primary statistic F30** = share with τ_next < 30 d (F14 column).
- **Control 1 (S7):** family-preserving, segment-preserving day shuffle (1000×). **Control 2 (scrambled dates):** each
  major's 30-day pre-window is moved to a random position inside its own inter-episode gap (the same gap, so the same
  network in nearly every case), 1000×; statistic = events inside the windows.
- Also S7's "events ≥ 30 d after the previous end" condition and the burst-level version (each burst counted once).
- **Ramp supported** if F30 > 99th percentile of control 1 AND the ≥ 30 d version passes AND control 2 p < 0.01;
  **not supported** if F30 ≤ control-1 median; **undecided** otherwise. Superposed-epoch table per segment.
  Detection-free column: env/PNSN deep tremor counts in the same pre-windows vs control 2.

## 7. Network-change caveat (every study)
Sensitivity steps: era detector 8–10× the frozen network on ETS days (similar on quiet days); B1/B2/B3/S15 steps;
LRIV 2009-11-20; tremor source env (CN4, 4 fixed stations) before 2010-01-06, PNSN after. Every shuffle and surrogate is
kept inside a network segment; nothing is trended across segments; across-era statements use `rate` only.

## 8. Outputs
`RESULTS.md` (plain verdict first per study, then numbers, caveat, new vs 2014–26-only), `tables/*.csv|json`,
`figures/r23_*.png`; notes lines in `svi_science/NOTES.md` and `notes/2026-09-25_Notes.md`.

## Revisions
### Revision 1 — Sun 2026-10-11 ~12:13 am PT, after the foundation build (`r23_base.py`), before any study statistic
- **Acceptance number corrected.** RESULTS_BACKWARD's 51,539 / 305 is the total at P ≥ 0.75 for every detector. Under
  the mandated era P\* ≥ 0.90 rule the reference is its "data at era P\* 0.90" row (2003 3,730 + 3 CN4-frozen January
  rows, 2004 10,891, 2005 6,620) and the data row for 2007–2013. All match exactly (34,093 data / 85 control over all
  families). 2006 (859) has no printed reference (834 era + 25 frozen).
- **Episode validation rule strengthened (it passed a broken list).** Variant (a) met the written bars (recall 7/7,
  precision 1.0, 6/6 Bostock windows) only because it merged 2003-03 → 2010-01 into ONE episode: env windows located
  within 30 km of some family occur almost every day. The rule lacked a start/end criterion, which was reported only.
  Added: ≥ 6 of the 7 matched PNSN majors must start and end within ± 7 d. (a) fails it (Δstart −70…−312 d);
  (b) meets it (6/7; Δ ≤ 5 d except 2012-07 end +10 d) but fails precision (0.74 < 0.8; 23 env majors vs 7 PNSN majors).
  So the **pre-registered fallback** is used: 2003-01-01…2010-01-05 majors = the 6 Bostock windows (dates from his
  catalogue), minors = the 39 variant-(b) episodes that do not overlap them (`validated = False`).
  Consequence (caveat, not a rule change): before 2010 every major is a Bostock ETS, so mid-size majors such as
  2010-03 / 2012-04 have no pre-2010 counterpart in the "all majors" clock (they are minors there); the big-ETS clock is
  the like-for-like one. Episode rate is similar in every era (5–8 listed episodes per year).
- Result of the build: 145 episodes 2003–2026, 52 majors, **23 big** (split S_rate 1,640 Bostock-equivalent events):
  10 in 2003–2013 (6 Bostock + 2010-08, 2011-08, 2012-08, 2013-09), 13 in 2014–2026. PNSN n30 reproduces
  `exposure_svi` except 8 of 60,076 family-days (distance formula).
### Revision 2 — Sun 2026-10-11 ~12:30 am PT, after the first study results (POST-HOC columns; no verdict rule changed)
- **Study 3 diagnostic `--union`** (`tables/study3_union/`): the fallback's Bostock-window majors end where his scan ended,
  which is before the tremor stopped in 2003-02 (env tremor to 03-18 vs window 03-10) and 2008-05 (06-09 vs 06-02).
  Events in those tails count as "between episodes, just after the end". The diagnostic widens each 2003–2009 Bostock
  major to the union with its overlapping env-(b) tremor episode (start/end) and reruns study 3 unchanged. Reported beside
  the pre-registered numbers; the pre-registered verdict is not replaced by it.
- Study 3 (A): burst-level F30 column (each tremor-free burst once; S7's deviation column).
- Study 2: burst sizes of BE_broad events (between episodes, local tremor allowed) per network segment, because the
  early-2006 bursts and most deep between-ETS activity are of this kind (Q5: 53 % of SVI broad bursts have tremor ≤ 30 km).
- Code fix (no rule change): the forward catalogue has no `burst_day` flag (NA); "without burst-day rows" columns now
  treat NA as not-a-burst-day, so for 2014–26 they equal all rows.
### Revision 3 — Sun 2026-10-11 ~12:35 am PT (POST-HOC column, no verdict rule changed)
- `r23_posthoc.py`: the S7 quiet-month F30 restricted to the deep groups (tremor-free and broad), fallback and union lists,
  because the 2003–13 failure of the all-family version comes from shallow-family ETS-tail events. Column only.
- Bookkeeping: importing `r23_quiet` with `--union` for this column re-opened (truncated) the in-progress log
  `tables/study3/study3_union.log`; the full union output is in `tables/study3_union_run.out`.
