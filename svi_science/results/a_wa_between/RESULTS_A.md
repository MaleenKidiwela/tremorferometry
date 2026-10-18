# A results — Washington-NW repeat of S7 (LFEs between ETS episodes)

Run 2026-10-07 ~1:25–1:35 pm PT. Rules frozen beforehand in `PLAN_A.md` (S7 rules + a listed set of adaptations).
Scripts `svi_science/scripts/a_wa_build.py` (tables) and `a_wa_between.py` (tests, figures; CPU, ~20 s). Logs
`a_build.log`, `a.log`; all numbers in `a_summary.json`. No Merlin review (run directly at the caller's request); provisional.
"Time-reversed control" = every template reversed in time (main); "scrambled control" = station data shifted to other
times of the same day.

## Plain answer: does each SVI result repeat in Washington-NW?

| SVI result (S7) | Washington-NW, 2015–2026 (primary) | verdict |
|---|---|---|
| 1. Between-episode LFEs come from deeper families than their ETS share predicts | 40.2 km vs 32.3 km expected (p < 0.0001); by position too: Slab2 depth 47.6 vs 45.5 km, down-dip +27 vs +17 km (both p < 0.0001). Holds per burst, at P ≥ 0.85, without de-duplication, by source group. 2012–14 and 2009–11 agree by position. | **Yes** |
| 2. They come in bursts | 220 events on 63 days; random days for the same families fill 211 (p = 0.001). 73 % sit in 14 bursts of ≥ 3. Same in 2012–14 and 2009–11. | **Yes** (even more bursty than SVI) |
| 3. No ramp before the next ETS | Last 30 days before a major hold 27 % vs 35 % expected (fewer, as in SVI); per burst 31 % vs 34 % (p = 0.30). | **Yes** (no ramp) |
| 4. (provisional, pre-registered here) Fewer in the month after a major ETS ends | 32 % vs 35 % expected, p_low = 0.17. The frozen rule says **does not repeat**. | **No** (in the primary era) |

- The confirmatory test fails in 2015–26. The older eras point the SVI way: 2012–14 has 18 % vs 39 % (p = 0.001) and
  2009–11 has 13 % vs 23 % (p = 0.012). Both are secondary, so they cannot rescue the test.
- Looking only at the first two weeks (NOT pre-registered, descriptive) there are fewer events in both regions:
  WA 15 vs 37 expected, SVI 1 vs 9. In WA the shortfall is cancelled by an excess on days 14–30 (56 vs 41). At most this
  says "quiet for about two weeks, not a month". It needs its own pre-registered test before it can count.
- **How many events:** WA has 220 between-episode LFEs in 2015-01 … 2026-05-27, against 9 in the time-reversed control
  (4 %) and 4 in the scrambled control (2 %). SVI had 76 against 3. **The controls are low at 0.75.** At P ≥ 0.5 they are
  not: 998 vs 135 / 196 (14 % / 20 %). The cut stays at 0.75.
- **Differences from SVI:**
  - **The WA bursts often go with PNSN tremor 30–100 km away.** 55 % of events and 18 % of bursts have it, against
    7 % / 6 % on random days (p ≤ 0.002). SVI had 17 % vs 9 %, which did not pass. So many WA "between-episode"
    bursts are small tremor bursts just outside the 30-km radius, not isolated activity.
  - **They concentrate even more:** 5 families give 75 % of the events. The leader is source G019 (fid 75/126, 46 km,
    along −27 km, southern Hood Canal): 93 events on 4 days.
  - **The start-location test passes the frozen event rule but rests on two bursts.** Before the next big ETS, events
    sit 10.5 km from its start vs 27 km when re-paired (p = 0.0004, 33 events). But 28 of the 33 come from two bursts.
    Counting each burst once (8 bursts, 5 episodes) gives p = 0.07. **Can't tell.**

## Definition (frozen) and what changed from S7
- **Event:** core-network NW candidate with P_own_WAVE_PNET ≥ 0.75, counted once by dedup variant D.
- **Between-episode:** its UTC day lies outside all 118 NW episodes ± 2 d, and post_v2 calls the day 'quiet' for the
  family (no PNSN epicentre within 30 km).
- **Exposure:** scanned family-days that meet the same condition (rates_daily_v2, zero rows included).
- **Positions:** a copy takes its source's primary position (`families_v1_3` `nw_copy_of`).
- **Depth:** Hyp2000, quality A/B only.
- **Majors:** the 69 `major` episodes of `episodes_nw_final.csv`.
- Full adaptation list: `PLAN_A.md`. The P-file events reproduce `rates_daily_v2_dedup` exactly (data, every era/class;
  `a_build.log`).

| definition (P ≥ 0.75, counted once D) | E3 2015–26 data / time-rev / scrambled | E2 2012–14 | E1 2009-08–2011 |
|---|---|---|---|
| **BE (headline)** | **220 / 9 / 4** | 72 / 4 / 2 | 100 / 1 / 2 |
| outside episodes ± 2 d, any local tremor (BE_broad) | 934 / 9 / 8 | 402 / 4 / 2 | 254 / 2 / 2 |
| outside episodes, no pad | 965 / 10 / 8 | 421 / 4 / 2 | 267 / 4 / 2 |
| BE and no tremor within 60 km | 105 / 9 / 4 | 70 / 4 / 2 | 75 / 1 / 2 |
| BE, source-group once | 234 / 10 / 4 | 79 / 7 / 2 | 102 / 1 / 3 |
| BE, raw (not de-duplicated) | 279 / 10 / 4 | 109 / 7 / 2 | 124 / 1 / 3 |
| BE, P ≥ 0.85 | 110 / 0 / 1 | 39 / 0 / 0 | 40 / 0 / 0 |
| BE, P ≥ 0.5 | 998 / 135 / 196 | 386 / 59 / 134 | 531 / 57 / 121 |

- Rate: 0.64 BE events per 1,000 eligible family-days in E3. Another 229 counted-once events after 2026-05-27 are
  not used.
- Full table (incl. P ≥ 0.5 / 0.85 per counting): `inventory_definitions.csv`.

## 1. Who and when
- **Families.** 33 of 132 families (23 of 93 source groups) fire between episodes in E3. Top source groups, BE share vs
  ETS share (`groups_between_E3.csv`):
  - G019 (46.4 km): 43 % vs 11 %
  - G006 (42.5 km): 16 % vs 3 %
  - G005: 10 % vs 15 %
  - G032 (45.0 km): 8 % vs 2 %
- The same deep sources (G006, G032, G019) lead in 2009–11 (76 %) and, with G005, in 2012–14.
- **Depth bins** (A/B, E3): 74 % of BE events come from families deeper than 40 km, which hold 40 % of ETS events.
  - 45–60 km: 47 % of BE vs 17 % of ETS.
  - ≤ 25 km: 11 % vs 29 %.
- **Time** (`bursts_between.csv`, `figures/a_timeline.png`). Events per year:

  | 2015 | 2016 | 2017 | 2018 | 2019 | 2020 | 2021 | 2022 | 2023 | 2024 | 2025 | 2026 |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | 17 | 3 | 67 | 69 | 6 | 26 | 7 | 7 | 1 | 4 | 8 | 5 |

- **Bursts.** 62 bursts (gaps < 24 h); 14 have ≥ 3 events and hold 73 % of the events. Distinct days 63 vs 211 for the
  shuffle (1st percentile 205; p = 0.001). The largest:
  - 2018-11-18 22:18 → 11-19 14:05 PT: 53 events, G019 + G005. Two of them fall within 60 s of an M 4.1 earthquake
    74 km away (03:08 PT).
  - 2017-11-10 23:29 → 11-11 07:41 PT: 42 events, G019.
  - 2017-01-24 01:53–05:37 PT: 17 events, G006 + G032.
  - 2020-09-13 PT: 9 events, G005.
  - E1/E2 (`bursts_between_E1E2.csv`): 2013-05-25/26 PT (43 events), 2010-02-19/20 PT (25), 2011-03-12/13 PT (22).
- **Detectability:** live core stations, median 18 on BE events and 18 on eligible days (E2 34 / 35, E1 11 / 11). The
  BE events are not just the better-detected days.
- **Core-network changes:** 69 % of BE events lie within ± 7 d of a core-set change, vs 75 % of eligible family-days.
  Changes are so frequent (1,444) that this check cannot separate anything.

## 2. Relation to the ETS cycle (`superposed_epoch.csv`, `figures/a_superposed_epoch.png`)
- **Before the next major (E3), events (shuffle mean):**

  | days to next major | 0–7 | 7–14 | 14–30 | 30–60 | 60–120 | 120–240 |
  |---|---|---|---|---|---|---|
  | events (expected) | 2 (15.5) | 27 (22.7) | 30 (36.6) | 58 (51.2) | 79 (66.5) | 23 (24.3) |

  - F30 = 0.27 vs 0.35 (p_low = 0.01). The ramp test **FAILS**, in the same wrong direction as SVI.
  - For events ≥ 30 d after the previous end: 0.32 vs 0.36.
  - Per burst: 0.31 vs 0.34 (p = 0.30).
  - The near-empty 0–7 d bin may partly be the definition: tremor that starts early near a family makes the day 'tremor'.
- **After the previous major ends (confirmatory, E3):**

  | days after end | 0–7 | 7–14 | 14–30 | 30–60 | 60–120 | 120–240 |
  |---|---|---|---|---|---|---|
  | events (expected) | 7 (15.4) | 8 (21.6) | 56 (41.1) | 48 (55.4) | 85 (62.1) | 16 (24.3) |

  - F30_prev = 0.32 vs 0.35, p_low = 0.17 → **DOES NOT REPEAT** (pre-registered rule).
  - Per burst: 0.26 vs 0.35, p_low = 0.07.
  - BE_broad goes the other way: 0.44 vs 0.36 (more right after an ETS; SVI broad was 0.07).
  - E2: 0.18 vs 0.39 (p_low = 0.001, "repeats"). E1: 0.13 vs 0.23 (p_low = 0.012, same direction, not confirmed).
- **Secondary eras, ramp:** E2 fails (0.19 vs 0.38, wrong direction). E1 "passes" at event level (0.36 vs 0.17;
  ≥ 30 d after the previous end 0.41 vs 0.22). This rests on one burst (2010-02-19, 25 events, 24 d before the next major):
  per burst 0.21 vs 0.17, p = 0.38. Not a ramp.
- **Where the next ETS starts** (14 S3 big NW episodes, A/B families):
  - Event rule: 33 events (5 episodes), 10.5 km vs 26.6 km re-paired, p = 0.0004 → formally **pass**.
  - But 17 events are the 2017-01-24 burst (12.7 km before the 2017-02 start) and 11 the 2020-09-13 burst (before 2020-09).
  - Per burst (column added after the run): 8 bursts, p = 0.07 → **can't tell**.
  - PNSN-start column: 116 events, p = 0.0002 at event level, but per burst p = 0.91 (32 bursts, 21 majors). The
    recomputed PNSN starts disagree with S3 for several episodes (2015-12, 2017-02, 2023-03) because the box takes in
    distant tremor. Column only, as in SVI.

## 3. Tremor elsewhere, earthquakes, GNSS
- **PNSN tremor 30–100 km on the same day:**
  - Events: 54.5 % of BE events vs 6.8 % on shuffled days (99th percentile 10.9 %, p = 0.001) → **tied**. Within ± 6 h: 53.6 %.
  - Bursts: 17.7 % vs 6.5 % (p = 0.002).
  - Controls: 0 of 9 time-reversed and 0 of 4 scrambled.
  - SVI did not pass (17 % vs 9 %, p = 0.022).
- **Earthquakes within ± 60 s:** 10 of 220 (4.5 %; chance 0.6 %). These include:
  - 2 next to the M 4.1 of 2018-11-19 03:08 PT (74 km);
  - 3 in 2025-03-19..28 PT (M 1.2–2.5, 38–77 km);
  - 2017-05-25 16:27 PT, an M 1.6 7.5 km from fid 115.

  Controls: 0. Flagged in `between_events.csv`, not removed.
- **Same-minute multi-family hits:** 129 of 220 events have another family firing within 120 s. That is expected, since
  bursts span neighbouring sources. They are not removed.
- **GNSS stack** (S7 estimator; 12 WA stations; own ETS azimuths from S1):

  | day set | days used | mean offset (mm) | season-matched random days, mean ± σ | z | north z (CHCM, P437) |
  |---|---|---|---|---|---|
  | BE burst days | 41 of 58 | −0.05 | −0.09 ± 0.062 | +0.7 | +0.2 |
  | BE_broad burst days | 58 of 77 | −0.13 | −0.08 ± 0.052 | −0.8 | −0.3 |
  | NW minor-episode mid-days (reference) | 40 | 0.00 | −0.07 ± 0.059 | +1.2 | +1.2 |

  - No slip-like signal. The 3σ bound is about 0.19 mm per ~6 days.
  - As in SVI, the reference stack does not see the minor episodes either, so this is a bound only.

## 4. Down-dip (A/B families; ETS shares from the same families and counting)

| set (E3) | n | metric | observed | expected from ETS shares | 99th pct | p (deeper) |
|---|---|---|---|---|---|---|
| **BE** | 215 | Hyp2000 depth (km) | **40.2** | 32.3 | 34.5 | < 0.0001 |
| BE | 215 | Slab2 depth at epicentre (km) | 47.6 | 45.5 | 46.3 | < 0.0001 |
| BE | 215 | down-dip (km) | +27.5 | +17.4 | +21.5 | < 0.0001 |
| BE, all families (Slab2) | 220 | Slab2 depth | 47.4 | 45.2 | 46.0 | < 0.0001 |
| BE, one per burst | 60 | Hyp2000 depth | 39.1 | 32.3 | — | < 0.0001 |
| BE, P ≥ 0.85 | 110 | Hyp2000 depth | 42.0 | 32.3 | — | < 0.0001 |
| BE, source-group once / raw | 229 / 273 | Hyp2000 depth | 40.0 / 39.3 | 32.3 | — | < 0.0001 |
| BE_broad | 929 | Hyp2000 depth | 31.0 | 32.3 | — | 0.997 (not deeper) |
| time-reversed control | 9 | Hyp2000 depth | 30.0 | 32.3 | — | 0.70 |
| scrambled control | 1 | Hyp2000 depth | 42.5 | — | — | — |

- **Mean depth rises with the cut-off:** P ≥ 0.5 37.9, 0.75 40.2, 0.85 42.0, 0.9 43.2 km.
  The time-reversed control at P ≥ 0.5 is 35.5 km.
- **Older eras:**
  - E1: Hyp2000 40.5 vs 33.4; Slab2 47.8 vs 43.4; down-dip +29 vs +7 km (all p < 0.0001; per burst p = 0.07).
  - E2: by position deeper (Slab2 48.6 vs 45.0, down-dip +31 vs +15, p < 0.0001), but not by Hyp2000 depth (32.2 vs 32.6).
- **Caveat on Hyp2000 depth.** 23 A/B families have Hyp2000 depths of 0–25 km, while the slab lies at 29–54 km below
  their epicentres. Copies of one source disagree too: the G005 primary is at 12 km, its copies at 43 km; fid 75 is at
  16.6 km, its G019 primary at 46.4 km.
  - This pulls the ETS-share expectation down, which inflates the Hyp2000 gap.
  - In E2 it also hides the effect: G005, the most down-dip source by position (Slab2 54 km), carries a 12 km depth.
  - The position-based columns (Slab2 depth, down-dip km) are the more robust reading, and they agree in all three eras.
- **Per-family ratio** (BE + 0.5)/(ETS + 0.5) vs depth: Spearman 0.06 (p = 0.62, 85 families). Unlike SVI (0.32), WA has
  no smooth per-family trend. The deep excess comes from a few deep sources (G019, G006, G032), not from depth in general.
- **BE_broad is not deeper in E3.** BE_broad means "outside episodes, but with tremor ≤ 30 km". That is unlike SVI, where
  it was strongly deeper. Its shallow events are G005 / Hood Canal activity during small tremor bursts.

## Failed / not passed (as pre-registered)
- **Failed:** the post-ETS quiet month (confirmatory, E3), p_low = 0.17. The ramp before the next ETS also fails.
- **Passes the frozen rule, but driven by bursts:** start location (event level). Per burst p = 0.07, so treat it as
  "can't tell".
- **Passed:**
  - deeper than ETS share;
  - clustered in time;
  - tied to tremor 30–100 km away, at both event and burst level. This is a WA-only result; SVI missed the bar.

## Files
`PLAN_A.md`, `RESULTS_A.md`, `a_summary.json`, `a.log`, `a_build.log`, `inventory_definitions.csv`, `between_events.csv`
(220 E3 events, PT times, τ, tremor / earthquake / station-change columns), `between_events_control.csv`,
`between_events_E1E2.csv`, `between_broad_events.parquet`, `bursts_between.csv`, `bursts_between_E1E2.csv`,
`families_between_E{1,2,3}.csv`, `groups_between_E{1,2,3}.csv`, `depth_bins_E{1,2,3}.csv`, `superposed_epoch.csv`,
`majors_tremor_start.csv`, `gnss_stack.npz`, build tables `families_wa.csv`, `events_all.parquet`,
`exposure_famday.parquet`, `live_core_by_day.csv`.

Figures:
- `svi_science/figures/a_timeline.png`: along-strike vs time 2009–2026, depth colour, both controls, monthly counts.
- `svi_science/figures/a_superposed_epoch.png`: events vs days after the previous ETS end and before the next start, with
  the shuffle band; observed/expected for WA E1/E2/E3 next to SVI.

Deviations added after the first run (columns only, no rule changed): the burst-level location and tremor tests, the E1/E2
burst table, and the descriptive first-14-days comparison.
