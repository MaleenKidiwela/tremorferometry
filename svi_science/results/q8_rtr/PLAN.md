# Q8 — Rapid tremor reversals (RTRs) in the LFE catalogues

Written 2026-10-08 ~11:20 am PT, BEFORE any detector was run on any catalogue. Analysis only (no catalogue
development, no downloads). CPU ≤ 8 cores, RAM < 30 GB. Scripts `svi_science/scripts/q8_rtr_*.py`.

## 0. Known answer (the positive-control target)

Royer, Thomas & Bostock (2015 JGR) = Royer PhD thesis 2014 ch. 4 (on disk: `notes/royer_sources/royer2014_phd_thesis.txt`,
§4.3.4, Tables 4.1–4.2). They hand-picked ("subjectively") **64 RTRs** in 7 SVI episodes of the Bostock catalogue
(2003, 2004, 2005, 2008, 2010, 2011, 2012) as "linear streaks of LFE detections moving at high apparent velocities in
the opposite direction to the main slip front". Table values: duration 1:33–8:13 h (mean ~4 h), along-strike distance
4.7–43 km (5 of 64 < 10 km), speed 1.3–9.1 km/h (mean ~5 km/h); 6.8 % of detections per episode are in RTRs; 72 % of
RTRs occur "exclusively during periods of positive UDSS" (up-dip shear stress, + = thrust sense).
- Their list is transcribed into `tables/royer2015_rtr_list.csv` (start/end times assumed **UTC**, like the catalogue).
  Obvious typos (end before start: #11, #20; #44 month "02/08" → 02/09; #23 "22.12") are fixed by start + listed
  duration and flagged in a column. No other change.
- Their per-RTR family membership is not public; only times.

## 1. Data (headline; all columns below are carried, nothing discarded)

| set | events | positions | episodes |
|---|---|---|---|
| Bostock 2003–13 (SVI) | `p3_bostock/bostock_events_mapped.parquet` (115 of our families; times UTC; mag) | S0 frame along-strike (`along_km`), down-dip | `p3_bostock/episodes_bostock_final.csv` (10), window = start − 3 d … end + 3 d |
| ours SVI 2014–26 | `s3/events_svi.parquet`, series data, P_same ≥ 0.75, keep_xfam8 (counted once). Columns: P ≥ 0.5 once; scrambled-control series c1 | same S0 frame (same 115 families) | S3 big episodes (10, ≥ 300 headline) |
| ours WA-NW 2015–26 | `s3/events_nw.parquet`, selected_v2 (P_own ≥ 0.75) & counted_once (dedup D). Column: P ≥ 0.5 once | S3 NW source positions (100 array-confirmed sources) | S3 big NW episodes (14) |
| ours WA-SW | catalogue_v2 SW selected_v2: only 2008-05 has > 100 events (774); run as one extra episode if positions map; otherwise reported "too few" | families_v1_3 pos_lat/lon, own Slab2 frame | 2008-05 only |

## 2. Detector (fixed now)

Per episode, events = (t, x) with x = along-strike km (+ = NW, S0 frame for SVI).
- **Candidate lines**: x(t) = x0 + v (t − t0), t ∈ [t0, t0 + D]. Grid: |v| ∈ {1, 1.5, 2, 3, 4, 5, 6, 8, 10, 15, 20} km/h,
  both signs; D ∈ {1.5, 2, 3, 4, 5, 6, 8} h; span |v|·D ∈ [8, 80] km; t0 every 10 min; x0 every 1 km; the whole
  line must lie inside the 2nd–98th percentile range of the episode's event positions. Corridor half-width w = 3 km.
  Time bins 10 min, position bins 1 km.
- **Score**: N = events in the corridor (parallelogram). Expected μ = the count predicted if, inside the bounding
  rectangle [t0, t0 + D] × [x range of the line ± w], event times and positions were independent
  (μ = Σ_cells n_t(t)·n_x(x)/N_rect over the parallelogram cells). This removes stationary bursts and rate changes that
  hit the whole region at once. Score S = −log10 Poisson P(≥ N | μ).
- **Streak shape** (required): N ≥ 5; events in ≥ 3 of the 4 time-quarters of the line; Spearman ρ between time and
  position-in-direction-of-motion of corridor events ≥ 0.5.
- **Overlap removal**: candidates sorted by S; a candidate is dropped if > 30 % of its events (of the smaller set)
  belong to an already-kept candidate.
- **Direction relative to the main front (local, handles two-way episodes)**: onset map = 10th-percentile event time
  in 4-km along-strike bins (bins with ≥ 5 episode events). A kept streak is a **reversal (RTR)** if
  (i) behind the front: ≥ 60 % of the path's bins with an onset have onset < t0, and
  (ii) back-propagating: median onset of the last quarter of the path is ≥ 12 h EARLIER than that of the first quarter.
  Streaks with (i) but the opposite onset order = "forward streak"; others = "other". All classes are kept in the table.
- **Threshold S\***: chosen on Bostock 2003–13 only, as the smallest S at which the RTR count in the scrambled controls
  (mean of 5 copies, worst of the two control types) is ≤ 10 % of the RTR count in the data. Frozen and applied
  unchanged to all other catalogues; the control-based false fraction is reported for every catalogue, and the full
  count-vs-S curve (data vs control) is tabulated (rank, do not discard). Because S is a significance and sparse
  catalogues cannot reach the same S, a second column applies the same 10 % rule to each catalogue's OWN controls
  (S\*_own); both are reported, S\* (frozen) is the headline. Streak search keeps everything with S ≥ 2 (after
  local-maximum reduction in (t0, x0) per (v, D) over ± 1 h × ± 3 km).

**Scrambled controls** (same detector, same threshold):
- C1 "time-shuffled": each event's time redrawn uniformly within its own family-day (UTC day). Keeps each family's
  daily count (the daily front), destroys hour-scale structure.
- C2 "position-permuted": within each UTC day, family-day blocks swap positions at random (times kept). Keeps all
  bursts and tidal/diurnal timing, destroys within-day space–time coupling.

## 3. Step 1 — reproduce Royer (pass rule, fixed now)

On Bostock, Royer's 7 episodes. **PASS** if (a) RTR count over the 7 episodes is within 43–96 (64 ± 50 %), (b) the
median RTR start is ≥ 2 days after episode onset (onset = first day of the frozen episode), and (c) ≥ 70 % of Royer's 64
RTRs have a one-to-one match among our RTRs (start within ± 6 h; greedy by smallest |Δt|). Also reported: our RTRs
matched (precision), chance match rate (our list shifted by ±1…3 days), and a UTC vs UTC−7/−8 diagnostic (not used
for the pass). FAIL → debug the detector against Royer's list (logged as deviations) before going further.

## 4. Step 3 questions (each with a test and a control)

a. **Counts / timing / kinematics** per episode; days after onset; speed, span, duration. Regions and eras compared
   with Mann–Whitney. **Density matching**: our catalogues are ~20–50× sparser than Bostock's. Bostock episodes are
   thinned at random to the median headline count of our SVI big episodes (10 copies) and re-run → expected RTR count
   per episode at our density and recall of the full-density RTRs. Era comparison uses the thinned numbers.
b. **Same places?** Statistic = mean distance from each RTR start (and end) to the nearest RTR start (end) of a
   DIFFERENT episode. Control: each RTR's start replaced by the position of a random event of the same episode within
   ± 12 h of its start (2,000 draws). One-sided p.
c. **Start at the main front?** Front position at t0 = median position of the bins whose onset falls in
   [t0 − 24 h, t0]. Statistic = median |start − front|. Control: random back-propagating grid lines of the same v and D
   in the same episode within ± 24 h of t0 meeting criterion (i) (non-circular: same geometry, random place).
d. **Tides** (SVI only; S5 hourly τ_updip at each family; Bostock 130-family file 2003–13, setup file 2014–26):
   RTR "exclusively positive" = event-weighted mean UDSS of its families > 0 at every hour of [start, end].
   Expected = for each RTR, the probability that the same-length window at a random start inside the episode's active
   days (days with ≥ 1 % of the episode's events) is all-positive; Poisson-binomial p. Applied to (1) Royer's own list
   (families = all families with events in the window), (2) our Bostock RTRs, (3) our 2014–26 RTRs.
e. **Sizes and depth**: Bostock mag of RTR events minus the median mag of the same family in the same episode (non-RTR
   events) → mean; S6 per-detection size dM for 2014–26 the same way. Control: the same corridor moved to random times
   in the same episode (≥ 5 events), 1,000 draws. "Stronger updip": Spearman of each family's RTR share vs its
   down-dip position (Bostock; ours if n allows); control = the same with control-type C2 detections at the same S\*
   and with random-time corridors.

Columns everywhere: raw (no de-dup) and P ≥ 0.5 for ours; WA tides not available (no WA tide setup) → not tested.

## 5. Deviations log (filled during the run)

**V0 (as pre-registered) — FAIL** (run 11:19–12:01 pm PT; outputs moved to `tables/v0_prereg/`). S\* came out 48
(the position-permuted control C2 made MORE streaks than the data below S ≈ 15: 6,727 vs 1,635 at S ≥ 2). At S\*:
25 RTRs in Royer's 7 episodes (rule 43–96 ✘), median 13.9 d after onset (✔), 13/64 Royer RTRs matched = 20 %
(chance with the list shifted ±1–3 d: 13 %) ✘. Eye check (`/tmp` debug plot of 6 windows, UTC): Royer's windows sit
on dense diagonal bursts (time zone UTC confirmed: recall UTC 20 % vs UTC−7/−8 9–11 %), but V0 drew lines that hop
between scattered events and between stationary single-family bursts, and several crossing lines of both signs through
the same burst. Royer windows are mild rate surges: median 2.2× the ±24 h median hourly rate, median 42 families active.

**D1 (debug iteration 1, ~12:03 pm PT, before any re-run):**
- Expected count μ = the same positions' event rate smoothed over ± 24 h (and ± 1 km), not the space/time-independence
  product → a streak must be a rate SURGE moving along a line.
- Streak must involve ≥ 5 distinct families, no family > 50 % of its events, family median times ordered along the
  motion (Spearman over families ≥ 0.6), and ≥ 50 % of the 2-km bins along the path hit.
- One streak per space–time box: a candidate is dropped if it overlaps a kept streak by > 50 % of the shorter duration
  AND > 50 % of the smaller along-strike range (in addition to the 30 % shared-event rule).
- Everything else unchanged (grid, W, classes, S\* rule, pass rule). Bostock re-run first; our catalogues only after
  step 1 is settled.

**D1 result — FAIL** (`tables/*_bostock_only.*`, 12:28 pm PT): the position-permuted control still gives ~1/3–1/2 of the
data's streak count at any score (S ≥ 50: 151 vs 70 per 10 episodes) → S\* = ∞ under the 10 % rule. With ~40 families
firing per hour in Bostock's catalogue, some ≥ 5-family subset always lines up. Royer-matched streaks are stronger
(median S 64 vs 23, 16 vs 9 families) but no (score, family-count) threshold gives ≥ 70 % recall at ~64 detections:
e.g. ≥ 15 families & S ≥ 50 → 76 streaks, recall 63 % vs 39 % with our list shifted ±1–2 days.

**D2 (debug iteration 2, ~12:33 pm PT) — back-zone surge detector** (`scripts/q8_rtr_backzone.py`): Royer picked by eye
on whole-episode along-strike-vs-time plots, where an RTR is a near-vertical dense column BEHIND the front. Rule:
onset map as before; slipped zone at time τ = 4-km bins with onset < τ − 12 h; count slipped-zone events per 30 min;
surge bin = count ≥ 3 and Poisson p < 0.01 against the same slipped set's ± 24 h mean; surge bins merged (gap ≤ 1 bin);
run ≥ 1.5 h; speed = Theil–Sen of the run's slipped-zone events; span = 10–90 % range; **RTR = span ≥ 8 km, |v| ≥ 1 km/h,
direction opposite to the episode's main front direction** (Theil–Sen of daily median position over days with ≥ 1 % of
the episode's events) — Royer's literal definition. Score S = Poisson significance of the run.
Prototype numbers on Royer's 7 episodes (seen before adopting it — so this is calibration on the known answer, not a
blind test): 72 RTRs, 58 % of Royer's matched within ± 6 h vs 32–37 % with our list shifted by ±12–48 h.
Variants tried in the same session (for the record, all below the 70 % bar at counts ≤ 96): local "away from front"
direction 48 RTRs / 41 %; surge + D1 best line for direction 74–106 / 59–67 % vs 41–50 % chance.
**Decision:** step 1 stays FAILED on the 70 % match rule (count and timing pass). The remaining questions are run with
D2 unchanged and labelled "detector reproduces Royer only partly". The S\* 10 % control rule is applied to D2's S;
everything else in §4 unchanged (front distance for (c) uses D2's own front estimate: reference = positions of
slipped-zone events in the same 30-min bins).

**D3 (~12:37 pm PT, after the D2 run, before the questions were computed):** the 10 % control rule cannot be met by D2
either: the position-permuted control keeps every rate surge (it keeps times) and assigns a random direction, so it makes
RTR-class and forward-class surges about equally (Bostock 46 vs 46 per 10 episodes) — it measures DIRECTION coherence,
not detection. The time-shuffled control almost never makes a surge (1.6 vs 88). Decision: no score threshold (every D2
RTR, S ≥ 2) for all catalogues; both controls reported as rates; per-catalogue own-control thresholds kept as a column.
Post-hoc additions (labelled as such in RESULTS): thinned-Bostock precision; start-histogram correlation for (b);
tides of non-migrating surges.

---

# Q8b — the same RTR finder on PNSN tremor 2010–2026 (pre-registered 2026-10-08 ~12:50 pm PT, before any PNSN number)

**Detector:** final D2 (`q8_rtr_backzone.py`), every parameter unchanged (4-km onset bins, ≥ 5 events per onset bin,
slipped zone = onset < τ − 12 h, 30-min bins, surge = count ≥ 3 and Poisson p < 0.01 vs ± 24 h, run ≥ 1.5 h,
RTR = span ≥ 8 km, |v| ≥ 1 km/h, opposite to the episode's main direction). No score cut (D3).

**Input mapping:** `catalogs/pnsn_tremor_cascadia_full.csv` (master; starts 2010-01-06, so "2009" is not available).
One "event" = one PNSN tremor window (epicentre, ~2.5-min windows). Position x = epicentre projected along strike;
down-dip coordinate likewise. Location error: PNSN tremor epicentres scatter by several km (envelope cross-correlation,
typical epicentral uncertainty ~3–5 km; depths not used). That is comparable to the 8-km minimum span and the 4-km onset
bins, so short RTRs will be blurred and some spurious migration from scatter is possible. The controls measure the
latter.

**Boxes and frames (same as S3):**
- SVI: epicentres ≤ 10 km from any of the 115 SVI families, S0 frame.
- WA-NW: ≤ 10 km from the S3 NW quality A/B source positions, S3 NW frame.
- WA-SW: ≤ 10 km from the 81 `families_v1_3` Royer_SW_WA positions (pos_lat/lon), frame by the S3 rule (circular mean
  Slab2 down-dip azimuth − 90°, origin = mean position).
- A window can fall in two boxes (as in S3); each box is analysed separately.

**Episodes:** the PNSN-era big-ETS lists already in use (`b_ets_record/episodes_b.csv`, `big == True`, windows w0–w1,
episode start = list start date).
- SVI: 18 episodes 2010-08 … 2025-09 (SVI list extended to 2010 by B with the same rule).
- NW: the 28 NW big episodes from 2010 on.
- SW: no SW episode list exists; the NW big-episode windows are used for the SW box, a stated limitation.
- Column: the S3 big subset.

**Controls (exactly as before):** C1 time-shuffled (time redrawn within its unit-day) and C2 position-permuted
(unit-days swap positions within each UTC day), 5 copies each. The PNSN "unit" (no families) = the 2-km along-strike bin.

**Cross-check (2010–2013 overlap, SVI):**
- PNSN-SVI RTRs vs the Bostock-based D2 RTRs inside the Bostock episode windows 2010-08, 2011-08, 2012-08, 2013-09.
- One-to-one match of start times within ± 6 h (greedy), both ways (recall of Bostock RTRs; precision of PNSN RTRs).
- Chance = our PNSN list shifted by ± 1, 2, 3 days.
- Secondary: the same vs Royer's 2010–2012 list. No pass bar is set (the catalogues differ). Stated as agreement vs chance.

**Questions (same tests and controls as Q8 §4):**
- a. Count per episode, days after onset, speed, length, duration; backward (RTR) vs forward counts against C2's 50/50.
- b. Repeat sites: nearest start in another episode vs random event within ± 12 h (2,000 draws), plus the 10-km
  start-histogram correlation (pre-registered here as secondary).
- c. Start-to-front distance (as Q8).
- d. Tides, SVI only: each epicentre gets the S5 hourly up-dip shear of its nearest SVI family (Bostock-era file before
  2014-01-15, setup file after). Same "whole RTR positive" and "mean positive" tests. WA has no tide setup → not tested.
- e. Up-dip: "families" replaced by 4 × 4 km (along × down-dip) cells with ≥ 50 windows. Spearman of RTR share vs cell
  down-dip position, with the random-time window control as before. No sizes (PNSN tremor has no usable amplitude).
- SVI vs WA: Mann–Whitney on speed, length, duration, timing; RTRs per episode by Poisson rate comparison.
- Eras: PNSN-SVI RTRs per episode 2010–2013 vs 2014–2026 (conditional binomial test on counts with episode numbers as
  exposure); same for NW. This is the clean era comparison, because it is one catalogue.

---

# Q8c — the same RTR finder on ALL our candidates (no P cut) (pre-registered 2026-10-08 ~5:35 pm PT, before any Q8c number)

**Detector:** final D2 unchanged (as Q8b). No score cut.

**Inputs (every ≥ 8 MAD candidate):**
- **SVI 2014–26:** `s0/analysis_table.parquet`. series `data` = real candidates, series `c1` = scrambled-control
  candidates (stations shifted to different times). Positions = S0 family along-strike/down-dip.
- **WA-NW and WA-SW 2015–26:** `wa_prod/post_v2/dedup/flags/flags_*.parquet` (every candidate, series data/c1, family
  id, time, dup_D). P_own_pnet is joined from `catalogue_v2` for data rows. Scrambled-control rows have no P.
  - NW position = the S3 source position of the family (families with no source dropped, as S3).
  - SW position = `wa_prod/setup/families.csv` lat_use/lon_use (reg == SW), frame by the S3 rule.
- **Calibration 2010–13 (SVI):** the backward production's frozen-network forward-detector candidates
  (`backward/mf/feat/p5f_<year>/*.parquet`, series data/c1, keep_xfam8), inside the Bostock windows 2010-08, 2011-08,
  2012-08, 2013-09. No P exists yet for these (P4 scoring not finished), so only variants (a) and (b).

**Episodes:** the same as Q8 — S3 big episodes (SVI 10, NW 14); SW uses the NW windows (as Q8b).

**Variants (all reported):**
- (a) raw: every candidate.
- (b) counted once: SVI/2010–13 keep_xfam8 = True; WA dup_D = False.
- (c) P-weighted: D2 counts events, so weights are applied by random thinning — each candidate is kept with
  probability P (5 draws, mean reported). This keeps the Poisson logic of D2 valid. SVI scrambled-control rows have P
  (P_same), so (c) is also run on them. WA scrambled-control rows have no P, so (c) has no scrambled-control run for WA.
- Q8's P ≥ 0.75 counted-once results are the reference column.

**Controls:**
- C1 time-shuffled and C2 position-permuted, 5 copies each (as before), for (a) and (b).
- **Scrambled-control candidates:** the same finder on series c1 with the same variant rule. Per-episode RTRs there =
  the false-RTR rate from noise.
- **Pass rule (per variant and region):** data RTRs per episode exceed scrambled-control RTRs per episode by ≥ 2×,
  AND the paired episode-level comparison (data − scrambled per episode) has one-sided Wilcoxon p < 0.05. A variant that
  fails is reported, but its science answers are labelled "not separable from noise".

**Calibration 2010–13:**
- RTRs from (a) and (b) on the backward candidates vs the Bostock-based D2 RTRs (Q8) and vs Royer's 2010–12 list.
- One-to-one ± 6 h; chance = list shifted ± 1–3 d. Also the scrambled-control series there.

**Thinned-Bostock expectation at the new density:** Bostock is thinned to the median per-episode count of each
variant (10 copies each). Reported: expected RTRs per episode, precision and recall vs full-density RTRs.

**Questions (same tests as Q8/Q8b):** count per episode vs the thinned expectation; speed, length, duration; backward
vs forward (binomial, with C2 as the 50/50 check); tides (SVI; RTR, forward and non-migrating bursts separately);
up-dip (family RTR share vs down-dip, random-window control); repeat sites (nearest-start and histogram correlation);
start-to-front; SVI vs WA (Mann–Whitney, rate test); 2014–26 vs 2003–13 (counts per episode vs thinned-Bostock
expectation at the same density, plus the 2010–13 backward-candidate calibration run).

---

# Q8d — P for the WA scrambled-control candidates, then the WA P-weighted RTR test (pre-registered 2026-10-08 ~5:55 pm PT)

**Source of P:** no refit, no recomputation. The frozen WA model (`wa_prod/conf/MODELS_FROZEN.json`, sha256
e184734b…, HGB + isotonic, label P_own_WAVE_PNET) was already applied in production to every series (data, c1
scrambled control, rev), with the same WAVE + PNET features: `wa_prod/conf/P/prod/P_YYYY-MM.parquet`.
- Check done before use: the join key (date, series, fid, k) is unique in both the P files and the dedup flag files.
  On 2018-05, data rows reproduce catalogue_v2 P_own_pnet to 3 × 10⁻⁸, and every c1 row has a P (0 missing).
- The scrambled-control rows of NW/SW families 2015–26 with their P are written to a NEW file,
  `q8_rtr/tables/wa_scrambled_P_2015_2026.parquet`. No WA catalogue file is edited.

**Test:** exactly Q8c variant (c).
- Each candidate is kept with probability P, 5 draws, for data and scrambled-control rows.
- Pass rule: data ≥ 2× scrambled AND per-episode Wilcoxon p < 0.05.
- Same episodes (S3 NW big; SW uses the NW windows).

**Questions:**
- backward vs forward (binomial; plus the position-permuted control on the P-weighted draws);
- up-dip (family RTR share vs down-dip, random-window control);
- repeat sites and start-to-front, as Q8c;
- per-episode counts vs SVI P-weighted (2.54 per episode), by Poisson rate comparison, scaled also per 1,000 P-sum.
- Tides: not testable for WA (no WA tide setup); stated, not improvised.

---

# Q8e — why WA-NW P-weighted bursts show no backward preference (pre-registered 2026-10-08 ~10:05 pm PT, before any Q8e number)

**Base:** Q8d WA-NW P-weighted test.
- D2 unchanged, P-thinning, 5 draws.
- Data and scrambled-control candidates, plus the position-permuted control on the draws (5).
- S3 big NW episodes (14) unless stated.

**Variants (all run, all ranked; none chosen in advance):**
- **V0 baseline:** S3 source positions, all families with a source, S3 NW frame (strike 344.73°), LFE-candidate onset map,
  LFE episode main direction.
- **P-a one family per source group:** only each source's `best_fid` (the family whose location is the source
  position). 100 sources.
- **P-b well-located only:** sources with quality A in `s3/sources_nw.csv` (erh-based, 55 of 100), all their families.
  The criterion is fixed now, before results.
- **P-c both:** best_fid of quality-A sources.
- **F-1 tremor front:**
  - Main direction per episode = sign of the S3 PNSN tremor Theil–Sen speed (`trem_speed`).
  - Slipped zone and front = onset map (same rule: 4-km bins, 10th percentile, ≥ 5 windows) built from PNSN tremor
    windows within 10 km of the A/B sources (S3 box), same frame, same episode window.
  - Positions and candidates as V0.
- **F-2 single-clear-front episodes:** V0 restricted to episodes whose S3 `front_class` starts with "one front"
  (11 of 14; drops 2017-03, 2018-05, 2020-10).
- **F-3:** F-1 restricted to the same 11 episodes.
- **A-1 slab-contour axis:**
  - Strike = PCA direction of the Slab2 contour at the sources' median slab depth (41.8 km ± 1 km) inside the NW source
    box = 345.24° (computed before results; 0.5° from the current 344.73°).
  - Positions re-projected on this axis, same origin. A curvilinear contour coordinate is not run, because the contour
    is straight across the box (Slab2 longitude at 41.8 km moves smoothly from −122.88 to −123.33 between 47.3° and
    48.4° N).

**Tides for WA-NW:** the S5 machinery unchanged (`s5_bostock_stress.py` functions):
- geometry = Slab2 plane within 10 km of each source, 236° slip direction;
- body tide = pyTMD IERS SET → plane stress;
- ocean = Boussinesq on the saved S5 EOT20 + gauge cell fields. The S5 inner polygon includes Puget Sound and Hood
  Canal gauges, and all NW sources are within the 330-km inner zone.
- Run at the 100 NW source positions. Depth = **Slab2 depth at the source** (Hyp2000 depths are unreliable for some
  sources).
- Hourly 2015-01 to 2026-09 UTC → `q8_rtr/tables/tidal_stress_nw_sources_hourly.parquet` (new file).
- Tide test as before (whole burst in positive up-dip shear vs random same-length windows on active days; Poisson-
  binomial). Run for data RTR, forward and non-migrating bursts and for scrambled-control RTRs, for every variant.
- Pooled draws: p is also given with n divided by 5.

**"Fixed" rule (per variant; all three must hold):**
1. Backward share (RTR / (RTR + forward)) exceeds the position-permuted control's backward share.
2. One-sided binomial p < 0.01 of the backward count against the control's share. Counts are pooled over the 5 draws
   and divided by 5 (effective n).
3. Data RTRs exceed scrambled-control RTRs by the Q8c/Q8d rule (≥ 2× and per-episode Wilcoxon p < 0.05).

**Also reported per variant:** RTR and forward counts, burst length and speed vs the position-permuted control,
family location table (location method, erh/erz, quality, copy/source-group flags from `data/catalogue/families_v1_3.csv`).

**Q8e post-hoc (added after seeing the 8 pre-registered variants, ~10:20 pm PT; labelled post-hoc everywhere):**
- Pd = P-weighted counted-once candidates (dup_D == False, for data and scrambled rows alike).
- Pe = Pd restricted to quality-A sources.
- Reason: S0b found 31 % of NW repeats are cross-family duplicates with the winner a median 19.5 km away, which would
  smear positions. Same tests and "fixed" rule.
