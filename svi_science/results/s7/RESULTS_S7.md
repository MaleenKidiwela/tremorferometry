# S7 results — LFEs between ETS episodes (southern Vancouver Island)

Run 2026-10-07 ~4:00 pm PT. Rules frozen beforehand in `PLAN_S7.md`. Script `svi_science/scripts/s7_between.py`
(CPU, ~15 s). Log `s7.log`; all numbers in `s7_summary.json`. No Merlin review (run directly at the caller's request);
provisional until Washington-NW repeats it.

## Plain answer
- **The LFEs that fire between ETS episodes come from the deep families.** Their mean depth is 39.2 km. If they
  followed each family's share of ETS activity, it would be 35.4 km (99 % of draws ≤ 36.6 km, p < 0.0001). Half of
  them (51 %) come from families deeper than 40 km, which carry only 22 % of the ETS events. The deep northern group
  G10 alone gives 30 % of them (6 % of ETS events). The busiest ETS group, G2 (shallow, south-east), gives none.
  The result holds when each burst counts once, on the fixed-13-station run, and for the scrambled control.
- **They come in bursts, not as a steady trickle.** The 76 events fall on 44 days. Random days for the same
  families would fill 75 (p = 0.001). Over half of the events sit in 7 bursts of 3 or more. The two biggest:
  - 2015-04-12/13 PT: 11 events in 3 deep G10 families.
  - 2018-10-21/22 PT: 11 events in 8 families of groups 1, 5, 7, 9 and 11, with PNSN tremor 30–100 km away that day.
- **They do not announce the next ETS.** Activity does not pick up before the next major episode. Counting events,
  there are fewer than expected in the last 30 days (13 % vs 27 %). Counting bursts, the timing is random (21 % vs
  26 %, p = 0.29). There are too few events before a mapped ETS (10 events, 3 episodes) to test whether they sit where
  the next ETS starts. The tremor-start version of that test shows nothing (p = 0.29).
- **No link to slip shows up.**
  - GNSS on those days: −0.12 mm vs −0.06 ± 0.07 mm on random days. That bounds slip at < 0.2 mm per ~6 days,
    about 4 % of a big ETS.
  - The same stack does not see the 40 minor episodes either, so this is a bound, not evidence of no slip.
  - Distant tremor (30–100 km) on the same day: 17 % vs 9 % (p = 0.02), which misses the p < 0.01 bar. Suggestive only.
- **A provisional extra:** LFEs are scarce in the first month after a major ETS ends (14 % vs 27 % expected,
  p = 0.008). The broad set shows the same (7 % vs 27 %, p = 0.001).

## Definition (frozen)
- A BE event:
  - P ≥ 0.75;
  - counted once across neighbouring families (keep_xfam8);
  - its UTC day lies outside all 79 SVI episode windows ± 2 days;
  - no PNSN epicentre lies within 30 km of its family that day.
- Period 2014-01-01 … 2026-05-27, which is where the PNSN catalogue ends. Another 33 events after that date are not used.
- All rates are per eligible family-day (family observable, same condition).

| definition (prod run) | data P ≥ 0.75 | scrambled control | data P ≥ 0.5 | control P ≥ 0.5 |
|---|---|---|---|---|
| **BE (headline)** | **76** | **3** | 763 | 403 |
| outside episodes ± 2 d (BE_broad, any local tremor) | 688 | 11 | 2,268 | 639 |
| no PNSN epicentre in the whole box that day | 161 | 5 | 1,009 | 446 |
| quiet flag (box quiet, ≥ 3 d from tremor) | 51 | 2 | 475 | 229 |
| outside the 79 episodes, no pad | 713 | 12 | 2,368 | 678 |
| BE, not de-duplicated | 118 | 9 | 1,228 | 583 |
| BE, fixed-13-station run (to 2021-03) | 51 | 7 | 513 | 300 |

P ≥ 0.5 is unusable here: the scrambled control is 53 % of the data. At 0.75 it is 4 %.

## 1. Who and when
- **Families.** 37 of 115 families have at least one BE event, and 5 families hold 34 % of them:
  - 119 (G10, 45.9 km): 9 events
  - 114 (G10, 42.2 km): 6
  - 33 (G1, 43.3 km): 4
  - 82 (G7, 36.4 km): 4
  - 76 (G7, 34.4 km): 3
- **Groups** (`groups_between.csv`). BE share vs ETS share by group:
  - G10 (43.6 km): 30 % vs 6 %
  - G1 (41.3 km): 14 % vs 9 %
  - G7: 14 % vs 12 %
  - G2 (33.5 km): 0 % vs 17 %
  - G5: 3 % vs 10 %

  The broad set (any local tremor) is even more one-sided: G10 42 %, G4 36 %, G1 13 %, which is 89 % from families
  deeper than 40 km. These are the small deep tremor bursts that did not make the episode list.
- **Time** (`bursts_between.csv`, `figures/s7_timeline.png`). Events per year:

  | 2014 | 2015 | 2016 | 2017 | 2018 | 2019 | 2020 | 2021 | 2022 | 2023 | 2024 | 2025 | 2026 |
  |---|---|---|---|---|---|---|---|---|---|---|---|---|
  | 6 | 17 | 9 | 2 | 12 | 3 | 4 | 6 | 1 | 0 | 10 | 1 | 5 |

  - Bursts are runs of events chained by gaps under 24 h. There are 40 of them; 7 have ≥ 3 events and together hold
    54 % of all events.
  - Biggest bursts: 2018-10-21 22:44 → 10-22 07:04 PT (11 events), 2015-04-12 11:13 → 04-13 08:34 PT (11),
    2016-11-24 PT (5, G10/G12), 2014-04-26–28 PT (4, G1), 2026-04-27/28 PT (4, G4, deep; 16 d after the 2026-04 episode).
- **Detectability.** The median n_sta_day is 14 for BE events and 13 for all eligible days. The extra station is
  small and cannot make the deep families stand out.

## 2. Relation to the next ETS (`superposed_epoch.csv`, `figures/s7_superposed_epoch.png`)
- Event counts by days to the next major start (shuffle mean in brackets):

  | days to next major | 0–7 | 7–14 | 14–30 | 30–60 | 60–120 | 120–240 | > 240 |
  |---|---|---|---|---|---|---|---|
  | events (expected) | 1 (3.4) | 4 (5.2) | 4 (10.1) | 11 (15.2) | 13 (15.9) | 19 (13.0) | 19 (7.3) |

  - F30 (share of events in the last 30 days) is 0.13 against 0.27 in the shuffle. The ramp test **FAILS**: it is in
    the wrong direction, p_low = 0.005.
  - The same holds for events ≥ 30 d after the previous ETS: 0.12 vs 0.27.
- **Counting each burst once** (38 bursts), F30 is 0.21 vs 0.26 (p_low = 0.29). The event-level deficit comes from a
  few big bursts that fell far from an ETS.
  - Share of bursts more than 240 d before the next major: 21 % vs 11 % (p = 0.04).
  - Median time to the next major: 106 d vs 62 d for random timing.
  - In plain words: the timing is consistent with random, maybe slightly favouring long quiet gaps.
- **After the previous ETS ends** (the ± 2 d pad already removes the first 2 days):
  - BE: 0 events in 0–7 d (expected 3.4), F30 = 0.14 vs 0.27 (p_low = 0.008).
  - BE_broad: 0 in 0–7 d (expected 30), F30 = 0.07 vs 0.27 (p_low = 0.001).
  - The deep families seem to go quiet for a few weeks after a big ETS. This is provisional: the definition depends
    on the episode window, and the end of an episode is set by tremor dying out.
- **Where the next ETS starts.**
  - Against the S3 LFE start points: only 10 BE events (3 episodes) fall ≤ 60 d before one of the 10 mapped episodes
    → **too few to test** (rule: ≥ 20).
  - Column, PNSN-tremor start for all majors: 20 events, median distance 37 km vs 43 km when re-paired (p = 0.29) → no.
  - The broad set "passes" against S3 starts (p = 0.001), but 45 of its 49 events precede a single episode (2024-08).
    It is not evidence.
  - The recomputed PNSN start disagrees with S3 for 2021-10 (−63 vs +18 km) because the box takes in Washington
    tremor. It is a column only.

## 3. Small slip, tremor and earthquakes
- **Tremor 30–100 km away on the same day:** 17.1 % of BE events vs 9.2 % for random eligible days of the same
  families. The 99th percentile is 17.1 % and p = 0.022, so the test is **not passed** at the frozen bar (suggestive).
  - Within ± 6 h: 18.4 %.
  - Scrambled control: 0 of 3.
  - Most of the excess is the 2018-10-21/22 burst, which had 6 epicentres 30–100 km away.
- **Catalogued earthquakes within 60 s:** 2 of 76 events (chance 0.6 %).
  - 2024-07-30 16:48 PT, family 17: an M 0.5 earthquake 35 km away, 7 s apart.
  - 2021-08-17 11:12 PT, family 108: an M 1.5 earthquake 250 km away.
- **Same-minute multi-family hits:** 13 events have another family firing within 120 s. Example: 2024-07-05 09:37 PT,
  three families 45 km apart. These may be one event seen by several families, or not an LFE. Flagged in
  `between_events.csv` (via time) and not removed.
- **GNSS stack** (Merlin S2 step 7; 14 SVI stations; own ETS azimuths; offset = 5-day median after minus before):

  | day set | days used | mean offset (mm) | random-day control mean ± σ | z | north z |
  |---|---|---|---|---|---|
  | BE burst days | 30 of 38 | −0.12 | −0.06 ± 0.066 | −0.9 | +0.3 |
  | BE_broad burst days | 64 of 74 | −0.04 | −0.06 ± 0.046 | +0.5 | +1.7 |
  | minor-episode mid-days (reference) | 40 | −0.02 | −0.06 ± 0.057 | +0.7 | +0.2 |

  - No slip-like signal: the 3σ bound is about 0.2 mm per ~6 days (BE) and 0.14 mm (broad), against ~5 mm for a
    big ETS.
  - The reference row shows the stack cannot see minor-episode slip either. A null here therefore says nothing
    beyond the bound.
  - The control mean of −0.06 mm is the normal loading drift between ETS.

## 4. Down-dip pattern
| set | n | mean depth (km) | expected from ETS shares | 99th pct | p (deeper) | mean down-dip km (expected) |
|---|---|---|---|---|---|---|
| **BE** | 76 | **39.2** | 35.4 | 36.6 | < 0.0001 | +11.6 (+0.3) |
| BE, one per burst | 40 | 38.5 | 35.5 | 37.0 | < 0.0001 | — |
| BE, fixed-13 run (to 2021-03) | 51 | 38.5 | 35.3 | 36.6 | < 0.0001 | — |
| BE_broad | 688 | 41.9 | 35.4 | 35.8 | < 0.0001 | +24.0 (+0.3) |
| scrambled control (BE) | 3 | 33.3 | 35.5 | — | 0.78 | — |

- Per family, the ratio (BE + 0.5)/(ETS + 0.5) rises with depth: Spearman 0.32 (p = 0.003, 84 families with ≥ 10
  ETS events). Against down-dip distance it is 0.22 (p = 0.05).
- Mean depth by cut-off:

  | cut-off | data | scrambled control |
  |---|---|---|
  | P ≥ 0.5 | 36.2 km | 35.2 km |
  | P ≥ 0.75 | 39.2 km (76 events) | — |
  | P ≥ 0.85 | 40.6 km (15 events) | 33.1 km (1 event) |

  The deep preference grows with confidence and the scrambled control is not deep. Noise is not causing it.
- Map and depth histograms: `figures/s7_map_depth.png`.
- About half (49 %) of BE events still come from families shallower than 40 km, mostly in bursts (2018-10, 2021-12,
  2024-07).

## Failed or not-passed checks (as pre-registered)
- Ramp before the next ETS: **FAIL**. The event level points the other way; the burst level is random.
- Location vs the next ETS start: **too few** events (10 < 20). The PNSN-start column is not significant.
- Tied to tremor 30–100 km away: **not passed** (p = 0.022 > 0.01).
- GNSS slip-like signal: **none**; bound only. The reference stack does not see the minor episodes either.
- Passed: **deeper than ETS** (all variants) and **clustered in time** (p = 0.001).

## What is new vs prior work
- Already known:
  - Tremor between ETS and smaller, more frequent episodes concentrate down-dip (Wech & Creager 2011, Cascadia).
  - Deep SVI LFE families are active more often and more continuously (Bostock et al. 2015, 2003–2013 catalogue).
  - S7 confirms this with LFE families over 2014–2026, with de-duplication and a scrambled control.
- New here:
  - The deep preference is measured against each family's own ETS share. G10 is over-represented ~5× and the main
    ETS group G2 is never active between ETS.
  - It holds per burst, on fixed stations, and grows with confidence.
  - Between-ETS LFEs come in bursts of hours, some of them spanning ~40 km and several groups (2018-10).
  - They give no warning of the next ETS: there is no ramp, no tie to the start location, and no GNSS signal above
    ~0.2 mm.
  - A possible post-ETS quiet period of a few weeks, provisional.

## Files
- `between_events.csv` (76 events; PT times; τ to next/previous major; tremor and earthquake columns),
  `between_events_control.csv`, `between_broad_events.parquet`
- `families_between.csv`, `groups_between.csv`, `bursts_between.csv`, `inventory_definitions.csv`
- `superposed_epoch.csv`, `majors_tremor_start.csv`, `gnss_stack.npz`, `s7_summary.json`, `s7.log`
- Figures: `svi_science/figures/s7_timeline.png`, `s7_map_depth.png`, `s7_superposed_epoch.png`

Deviations from the plan, added after the first run as checks (columns only, no rule changed):
- the burst-level versions of the next-ETS and depth tests;
- the earthquake ± 60 s and same-minute multi-family columns;
- the lower-tail p for F30.

The plan's pool estimate (109) included post-PNSN days; the frozen period gives 76.
