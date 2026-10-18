# S3 + S4: what each big ETS episode did (SVI and Washington-NW)

## 0. Pre-registration (written 2026-10-07 ~1:15 pm PT, before any kinematic or patch number was computed)

**Data (headline).**
- SVI: `svi_science/results/s0/analysis_table.parquet`, run prod, series data, P_same ≥ 0.75, keep_xfam8 (counted once).
  Columns: raw (no de-dup) and P ≥ 0.5. Scrambled-control rows (series c1) go through the same steps.
- NW: `catalogue_v2.parquet`, region NW, selected_v2 (own_pnet ≥ 0.75), 2015-01-01 → 2026-09-15 (era E3); counted once =
  dedup variant D == False ("same event", any distance + same arrival). Columns: raw, P ≥ 0.5.
  Location = the array-confirmed source group (`source_groups_E3_array_confirmed.csv`, 100 sources); a source's position
  = its best-located member (smallest erh). A column flags sources whose best member is quality C.
- Event-to-episode assignment by time only: frozen episode window padded by 3 days each side (capped at half the gap to
  the neighbouring episode). This deliberately does NOT use the S0 "family exposure" rule, because that rule only admits
  an LFE when PNSN tremor was near the family, which would make LFE onsets copy the tremor onsets.
- Big episodes: counted-once headline ≥ 300 events in the padded window.
- Along-strike / down-dip frame: Slab2 (as S0); along-strike positive toward the NW strike azimuth.

**Per-episode quantities.** start = median along-strike position of the first 10 % of events (and of the first 2 days);
end = median of the last 10 %; speed = Theil–Sen slope of along-strike position vs time (km/day), 90 % CI by
day-block bootstrap; same for down-dip; secondary bursts = a 5-km along-strike bin (SVI and NW alike) that had ≥ 3 events
on a day, then ≥ 3 consecutive days with 0 events, then again ≥ 3 events on a day. PNSN tremor: epicentres within 10 km of
a family / source position (NW: quality A/B sources only), same frame, same quantities.

**Pass / fail rules (fixed now).**
1. *Migration is real* in an episode if the Theil–Sen speed lies outside the central 99 % of 500 time-shuffled copies
   (event times permuted within the episode) — i.e. p < 0.01, two-sided. Region-level: "episodes migrate" if ≥ half of the big
   episodes pass.
2. *Agreement with tremor*: LFE and PNSN speeds have the same sign and the LFE start lies within 15 km of the tremor start
   in ≥ 2/3 of the big episodes where both pass rule 1.
3. *Detectability check*: thinning every family to the episode's median per-family count changes the start by < 10 km and
   leaves the speed sign unchanged in ≥ 80 % of big episodes; family onset delay vs family count: Spearman ρ reported; if
   ρ < −0.3 with p < 0.01 the onset-delay map is called detection-driven.
4. *Start positions cluster* if their standard deviation across big episodes is below the 5th percentile of 10,000 sets of
   uniform draws over the region's along-strike span (5th–95th percentile of event positions). Secondary control (stricter):
   each episode's start replaced by the median of a random 10 % of its own events — clustering beyond that control is
   required to call it a "repeatable start site".
5. *Same direction*: binomial test of NW-ward vs SE-ward signs among rule-1-passing episodes, p < 0.05.
6. *Start relates to previous end*: Spearman ρ(start_i, end_{i−1}) over consecutive big episodes, permutation p < 0.05.
7. *Patch strength (S4)*: share of events per group (SVI: the 12 S0 groups; NW: 10-km along-strike bins) per big episode.
   Residual = Pearson residual against the group's share in all OTHER big episodes (leave-pair-out). Persistence if the
   mean residual correlation over episode pairs ≥ 90 days apart is > 0 with permutation p < 0.01 (control: group labels
   shuffled within each episode, 2,000 times). Catch-up if the mean over consecutive pairs (≥ 90 d apart) is < 0 with
   episode-order-shuffle p < 0.01, and the variance ratio of each group's cumulative residual is below the order-shuffle
   5th percentile. Neither = "patches vary episode to episode without memory".

---

## 1. The answer in plain words (written after computing, 2026-10-07 ~2:30 pm PT)

- **Every big episode is a front that sweeps along the coast at about 5 km per day** (SVI 2.7–6.7 km/day, NW 2–12.5 km/day).
  It moves either toward the NW or toward the SE; it is about 50/50 which (SVI 4 vs 5, NW 5 vs 6). PNSN tremor, which
  is detected completely independently, shows the same sweep in 9 of 10 SVI episodes (7 of 12 in Washington-NW).
- **Episodes do not keep starting in one place inside either box.** They start at one END of the box or the other —
  where the front came in from outside — never in the middle. In SVI: fronts moving SE start near Lake Cowichan
  (+20 to +40 km), fronts moving NW start near Victoria / the Strait (−24 to −47 km). In Washington-NW: SE-moving fronts
  start near Sequim / the Strait (+20 to +28 km) and NW-moving ones near southern Hood Canal (−23 to −33 km).
- **Seen across both regions together, there are three kinds of episode:** one sweep from Washington up Vancouver
  Island (2024-08, and 2020-11 and 2021-10 start in Washington 4–8 days earlier), one sweep from Vancouver Island down
  into Washington (2015-12), or a start near the Strait of Juan de Fuca that splits both ways (2019-08, 2023-04). Three
  Washington episodes (2019-03, 2024-01, 2025-04b1) began near the Strait and ran SE fast (9–13 km/day) while SVI stayed quiet.
- **Where an episode stops is the main thing that decides which patches look "strong" or "weak".** Several SVI
  episodes only reached half the island (2019-08 only the SE half; 2019-11, 2021-01, 2022-03 only the NW half), and
  the patches they never reached are empty. Beyond that, there are **no patches that are always strong or always weak**,
  and no reliable "catch-up" (pre-registered tests fail; one hint, below).
- **Activity sometimes comes back** after ≥ 3 quiet days, more often than chance in 3 SVI episodes (2015-12, 2018-06,
  2023-04). The two clear cases are second fronts running back fast: 2018-06 activity re-appeared at the NW end on
  Jul 8–10 PT moving NW at ~15 km/day; 2023-04 activity re-appeared at the SE end on Apr 18–19 PT moving SE at ~24 km/day
  (opposite to the 4 km/day main front).
- **What is new vs known:** the sweep speeds, the both-ways directions and the back-running fronts confirm what was
  described from tremor (Wech & Creager; Houston et al. 2011; Rubin & Armbruster 2013; Peng & Rubin). New here: the same
  picture from LFE families over 10–12 years (2015–2026, 24 big episodes, with controls), and the negative result that
  start places and strong patches do not repeat; episodes differ mostly by which end they come from and how far they get.

Figures: `svi_science/figures/s3_svi_episodes.png`, `s3_nw_episodes.png` (one panel per episode; single-episode files
`s3_ep_<region>_<episode>.png`), `s3_svi_summary.png`, `s3_nw_summary.png` (map of starts + arrows; start→end vs date),
`s3_cross_region.png` (SVI and NW on one along-margin axis), `s3_patch_heatmap.png`.

## 2. Each big episode (headline: counted once, P ≥ 0.75; positions = along-strike km, + = NW; times PT)

Place names: SVI −40 km ≈ Strait south of Victoria, −25 ≈ Sooke, 0 ≈ Sooke Hills, +20 ≈ Shawnigan/Cowichan, +38 ≈ Lake
Cowichan. NW −33 ≈ southern Hood Canal, 0 ≈ central Hood Canal (Quilcene), +28 ≈ Sequim.

| region | episode | first LFE (PT) | LFEs | start → end (km) | speed km/day (90 % CI) | tremor start / speed | what happened |
|---|---|---|---|---|---|---|---|
| SVI | 2015-12 | 2015-12-18 | 1138 | +20 → −34 | −6.7 (−7.9…−4.7) | +30 / −6.4 | full sweep to SE; small returns at the start area; continued into WA (NW starts 5.8 d later) |
| SVI | 2018-06 | 2018-06-25 | 1498 | +35 → −36 | −5.6 (−6.1…−5.1) | +46 / −4.6 | full sweep to SE; second front back toward NW on Jul 8–10 (~15 km/day) |
| SVI | 2019-08 | 2019-08-12 | 662 | −47 → −8 | +4.5 (3.6…5.2) | −33 / +4.0 | started at the Strait, went NW, stopped mid-island; WA side went SE (split) |
| SVI | 2019-11 | 2019-11-08 | 613 | +25 → 0 | −2.7 (−5.1…−1.6) | +40 / −1.8 | NW half only, slow; returns not above chance |
| SVI | 2020-11 | 2020-11-06 | 617 | −38 → 0 | +6.2 (5.5…6.8) | −41 / +5.4 | came from WA (WA 4.3 d earlier), stopped mid-island |
| SVI | 2021-01 | 2021-02-03 | 394 | +40 → 0 | −6.3 (−6.9…−5.6) | +49 / −4.9 | NW half only |
| SVI | 2021-10 | 2021-10-04 | 344 | +24 → −10 | +0.8 (−2.6…5.4) | +18 / +0.6 | **no single front**: an early burst at the NW end, then a late NW-moving sweep from Oct 22; WA started 8 d earlier |
| SVI | 2022-03 | 2022-03-16 | 381 | +35 → 0 | −4.8 (−5.5…−3.5) | +50 / −4.7 | 11 LFEs at the SE end on Mar 17, main front from Mar 29 at the NW end, NW half only |
| SVI | 2023-04 | 2023-04-05 | 881 | −24 → +35 | +4.4 (3.1…5.4) | −17 / +4.1 | started near the Strait (split with WA), full sweep NW; second front back toward SE on Apr 18–19 (~24 km/day) |
| SVI | 2024-08 | 2024-08-15 | 1225 | −38 → +35 | +6.0 (5.6…6.7) | −36 / +5.5 | came from WA (4.0 d earlier), full sweep NW |
| NW | 2015-12 | 2015-12-29 | 864 | +28 → −23 | −4.4 (−6.1…−3.2) | +38 / −5.6 | continuation of the SVI front toward SE |
| NW | 2017-02 | 2017-02-18 | 706 | −23 → −1 | +4.9 (3.5…6.6) | −30 / +1.3 | went NW |
| NW | 2017-03 | 2017-03-17 | 372 | +18 → +28 | −1.0 | +21 / +3.6 | no single front (NW end only) |
| NW | 2018-05 | 2018-05-01 | 946 | −32 → −33 | −0.6 (p 0.014) | −25 / −3.4 | no clear migration (early SE-end burst, then whole box active) |
| NW | 2019-03 | 2019-03-11 | 628 | +20 → −23 | −8.7 (−13…−5.1) | +34 / −7.6 | started near the Strait, fast SE; SVI quiet |
| NW | 2019-08 | 2019-08-14 | 529 | +1 → −23 | −8.6 (−13.6…−1.5) | +40 / −8.7 | SE half of the split with SVI |
| NW | 2020-09 | 2020-10-05 | 821 | −33 → +25 | +6.2 (5.5…7.0) | −24 / +2.3 | went NW; then SVI 2020-11 |
| NW | 2020-10 | 2020-11-02 | 335 | +18 → +28 | −0.2 (p 0.68) | +22 / +3.4 | no clear migration (NW end, as the front moved on into SVI) |
| NW | 2021-09 | 2021-09-26 | 907 | −23 → +28 | +2.0 (1.2…2.6) | −30 / +2.3 | slow, patchy NW-ward (32 days); returns not above chance |
| NW | 2023-03 | 2023-04-05 | 766 | +6 → −32 | −8.0 (−10.6…−4.5) | +31 / −6.7 | SE half of the 2023 split |
| NW | 2024-01 | 2024-01-02 | 353 | +26 → −23 | −11.3 (−12.8…−7.0) | +36 / −7.8 | started near the Strait, fast SE; SVI quiet |
| NW | 2024-08 | 2024-08-10 | 539 | −4 → +28 | +3.1 (0.5…4.1) | −18 / −0.6 | went NW into SVI |
| NW | 2025-04b1 | 2025-05-15 | 341 | +26 → −23 | −12.5 (−16…−9) | −44 / −6.1 | started near the Strait, fast SE (tremor starts with a separate SE cluster 6 days earlier) |
| NW | 2025-09 | 2025-09-18 | 1188 | −32 → +28 | +4.2 (3.0…5.6) | −37 / +5.6 | went NW (SVI 2025-09 not big) |

Full columns (raw and P ≥ 0.5 variants, down-dip speeds, onset-based front speeds, thinning, bursts, tremor):
`episodes_kinematics.csv`; activity returns `bursts.csv`, `return_streaks.csv`; family onsets `family_onsets.csv`;
SVI/NW overlaps `cross_region_overlaps.csv`.

**Up/down-dip:** not resolvable as a consistent motion. SVI down-dip speeds are −1.4…+1.6 km/day with both signs; NW
−4.7…+3.4 km/day, but NW down-dip positions of single sources are uncertain (several best-located depths < 10 km).

## 3. Pre-registered checks (rules of §0)

| rule | SVI | NW |
|---|---|---|
| 1 migration real (time-shuffle p < 0.01) | **PASS** 10/10 (p = 0.002, the floor) | **PASS** 12/14 (2018-05 p 0.014, 2020-10 p 0.68) |
| 2 agrees with PNSN tremor (same sign, start ≤ 15 km) in ≥ 2/3 | **PASS** 9/10 | **FAIL** 7/12 (2017-03 and 2024-08 opposite sign; 2019-08, 2023-03, 2025-04b1 starts 24–70 km apart: NW LFE positions are lumpy, ~29 % of events come from one source at −23 km) |
| 3 detectability: thinning keeps start (< 10 km) and sign in ≥ 80 % | **PASS** 9/10 (2022-03 fails: thinning picks the 11-LFE early burst at the SE end) | **PASS** 12/14 (2019-08, 2021-09 starts move 34 / 12 km; signs kept) |
| 3b onset delay vs family count (ρ < −0.3, p < 0.01 = detection-driven) | 0/10 episodes | 0/14 episodes |
| 4 start positions cluster | **FAIL**: spread 35 km vs uniform 5th pct 16.7 km (and random-10 % control 14.8 km). The starts are MORE spread than uniform (p ≈ 0.0001): they sit at the two ends | **FAIL**: 24 km vs 14.7 km; more spread than uniform (p ≈ 0.02) |
| 5 same direction (binomial) | **FAIL** 5 NW-ward / 5 SE-ward, p = 1.0 | **FAIL** 5 / 7, p = 0.77 |
| 6 start relates to previous end | **FAIL** ρ = −0.12, p = 0.77 | **FAIL** ρ = −0.12, p = 0.69 |
| 7a persistent strong patches | **FAIL**: mean residual correlation (pairs ≥ 90 d apart) +0.103 vs patch-label shuffle 0.092 (p 0.094) and independent-episode simulation 0.101 (p 0.37) | **FAIL** (at the threshold): +0.087 vs 0.068 (p 0.0105); simulation 0.070 (p 0.068) |
| 7b catch-up | **FAIL (half)**: consecutive pairs −0.45 vs order shuffle +0.09 (p 0.001) ✔, but cumulative-residual test p 0.058 ✘ | **FAIL (half)**: consecutive −0.22 (p 0.030) ✘, cumulative-residual test p 0.0065 ✔ |

Raw (no de-dup) and P ≥ 0.5 give the same answers: speed sign identical in 24/24 episodes, speeds differ by a median
3 % (raw) / 7 % (P ≥ 0.5), at most 1.1 / 1.7 km/day; starts within 5 km in 21/24 and 20/24 episodes. Patch tests on raw
counts (`patch_tests.json`, `raw075`): same verdicts.

**Scrambled control** (stations shifted to different days, SVI only): 0–9 events per big episode at P ≥ 0.75 against
344–1,498 real → far too few to make a front; the fronts are not a detector artefact. NW has no scrambled rows in
catalogue_v2; its time-reversed control has 15 rows in all of 2015–26.

**Activity returns** (bin active, ≥ 3 quiet days, active again) vs the same count after shuffling event times: SVI
2015-12 3 vs 0.6 (p 0.01), 2018-06 5 vs 0.1 (p 0.005), 2023-04 4 vs 0.2 (p 0.005); 2019-11 3 vs 3.1 and 2021-10 2 vs
6.2 are chance. NW: 2023-03 1 vs 0.0 (p 0.01), 2024-08 3 vs 0.5 (p 0.02), 2025-09 3 vs 1.4 (p 0.16), 2015-12 5 vs 4.2,
2018-05 3 vs 7.1, 2021-09 10 vs 11.4 (chance).

## 4. Post-hoc (not pre-registered; treat as leads)

- **SVI patch pattern follows direction.** Episodes moving the same way have similar patch patterns (mean residual
  correlation +0.41), opposite-way episodes differ (−0.14); direction-label shuffle p = 0.018 (`post_direction_modes.json`).
  So the consecutive-episode anti-correlation (rule 7b) is mostly "the next episode came from the other end": direction
  changed in 7 of 9 consecutive SVI pairs (shuffle p 0.16, not significant). NW shows no such link (+0.10 vs +0.06, p 0.29).
- **In Washington-NW, SE-moving fronts were faster** (median 8.7 km/day, 4.4–12.5; n = 6) than NW-moving ones (4.2,
  2.0–6.2; n = 5). Three of the fast ones started near the Strait while SVI was quiet. Not tested; small n.
- **Start time between regions (6 overlaps):** Washington first by 4.0–8.0 days in 2020-10/11, 2021-09/10, 2024-08 (front
  moving NW); SVI first by 5.8 d in 2015-12 (front moving SE) and by 3.2–3.8 d in the two split episodes 2019-08 and 2023-04.

## 5. Changes from the pre-registration (stated)

- Theil–Sen uses only pairs of events at different positions. In NW one source holds ~29 % of the events and same-position
  pairs (slope 0) forced the median to exactly 0 (first run: NW 2018-05 speed 0.00, p = 1.0). Same change applied to the
  shuffle control and to SVI (SVI speeds changed by ≤ 0.6 km/day).
- Activity return: the return day may be the first or second day after the quiet gap (a day with 1–2 events starting the
  return no longer hides it).
- "no single front" class (start→end opposite to the fitted sign, or |speed| < 1.5 km/day) added for figures and the
  table; it reclassifies SVI 2021-10 and NW 2017-03, which pass rule 1.
- Theil–Sen and shuffles use a random 1,500-event subsample when an episode has more (PNSN tremor windows).

## 6. Prior work and what is new (honest)

- Known (confirmed here): ETS tremor in northern Cascadia migrates along strike at ~5–10 km/day, in either direction,
  sometimes bilaterally (Wech & Creager 2008 and later); secondary fronts run back along the already-slipped zone much
  faster than the main front (rapid tremor reversals: Houston et al. 2011 in Washington; Rubin & Armbruster 2013 and Peng
  & Rubin for southern Vancouver Island, 2003–2013). Our day-resolution returns (15–24 km/day over 20–25 km) are slower
  and longer than the hour-scale reversals those studies imaged; we cannot see the hour-scale ones at this binning.
- New: (a) the same front picture from de-duplicated LFE families, 2015–2026, 24 big episodes, with a time-shuffle
  control, a detectability check and an independent tremor check; (b) **no repeatable start site inside SVI or
  Washington-NW, no preferred direction, no link between where one episode stops and the next starts**; (c) patch
  strength is set mainly by how far the front reaches, with no lasting strong or weak patches (and, in SVI, a direction
  "mode" in the patch pattern, post-hoc). These are mostly negative results; the positive picture confirms known behaviour.
- Not done: hour-scale reversals; S1 GNSS timing cross-check; front speed vs episode size.
