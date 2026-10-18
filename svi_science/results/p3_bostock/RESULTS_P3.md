# P3 results — Bostock 2003–2013 baseline, the quiet month after ETS, and the recharge curve

Run 2026-10-07 ~1:45–2:00 pm PT. Rules frozen beforehand in `PLAN_P3.md` (two bookkeeping revisions at its end, no rule
changed). Scripts `svi_science/scripts/p3_bostock_map.py`, `p3_bostock_episodes.py`, `p3_quiet_month.py`,
`p3_recharge.py`, `p3_figures.py` (CPU, < 5 GB). No Merlin review (run directly at the caller's request); provisional.

## Plain answer
- **Bostock's catalogue maps cleanly onto our families.**
  - All 115 of our families are found, and every per-year count is reproduced exactly.
  - His ten ETS episodes come out of our own episode rule unchanged.
- **But Bostock only searched during the ten big ETS episodes.**
  - His 269,423 LFEs fall on 230 days, all inside 10 windows of 11–31 days. Between them, nothing was searched.
  - His paper says the same: the catalogue was "assembled from 10 major ETS episodes".
  - So the quiet month after ETS, the between-episode LFEs and the recharge curve **cannot be measured in 2003–2013**.
    The blank months are "not searched", not "quiet".
  - The 23-year test of these patterns that Merlin proposed is therefore not possible with this catalogue. The test
    can neither kill the result nor confirm it.
- **The quiet month after a big ETS does not replicate cleanly.**

  | dataset | regional end | end + 7 days | family's own end | verdict (frozen rule) |
  |---|---|---|---|---|
  | S. Vancouver Island 2014–26 (where it was found) | 14 % vs 29 % expected, p = 0.004 | 17 % vs 30 %, p = 0.007 | 14 % vs 27 %, p = 0.012 | holds only at the regional end (third version misses the 0.01 bar narrowly) |
  | Bostock 2003–13 | — | — | — | not testable |
  | Washington 2015–26 (confirmatory) | 32 % vs 36 %, p = 0.17 | 42 % vs 38 %, p = 0.91 | 32 % vs 30 %, p = 0.68 | **does not hold** |
  | Washington 2012–14 (older, secondary) | 18 % vs 39 %, p = 0.001 | 15 % vs 41 %, p = 0.001 | 18 % vs 36 %, p = 0.002 | holds under all three |
  | Washington 2009–11 (older, secondary) | 13 % vs 23 %, p = 0.009 | 12 % vs 27 %, p = 0.002 | 12 % vs 22 %, p = 0.009 | holds under all three |

  - The percentages are the share of between-episode LFEs that fall in the first 30 days after a major ETS ends.
    The "expected" value comes from moving each LFE to a random eligible day of its own family.
  - In Vancouver Island the quiet month survives moving the clock a week later. So it is not just the few calm days
    that the episode rule builds into every episode end.
  - It fails in the main Washington era. Washington's two older eras show it, but they are secondary.
  - Status: **found in one dataset, failed its pre-registered confirmation in the second, untestable in the third.**
    Merlin's minimum ("replicated in three datasets") is not met.
- **Deep activity builds up slowly through the whole gap between ETS, rather than "recharging" in a month.**
  - In 2014–2026, the deep families (down-dip groups 1, 4, 10, 12) fire about 3.6 times less often in the first month
    after a major ETS than later: 3.4 vs 11.9 events per 1,000 family-days.
  - Their rate then keeps rising for the whole ~400 days we can follow; it never levels off. The fitted recovery
    time hits the 400-day limit of the fit (95 % range 11–400 d).
  - **The frozen test does not pass with the primary clock (all 39 majors): p = 0.08 against time-shifted controls.**
    Part of the low early rate is expected anyway: long gaps between ETS are busier than short ones.
  - With the clock reset only by the ~14 big ETS (a column), the build-up passes (p = 0.002). The rate is near zero for
    a month, then rises over ~300 days (τ 322 d, 95 % range 85–400 d).
  - The up-dip families are almost silent between episodes (25 events), so they show no curve.
  - Scrambled control: 13 deep events in 12 years.
- **Is the cycle the same in both eras? We can't tell.**
  - The LFE curve cannot be built for 2003–2013.
  - The two substitutes have too little data to decide:
    - our own first vs second half (2014–19 vs 2020–26);
    - PNSN deep tremor, 2010–13 vs 2014–26.
  - Both give "no detectable difference". The uncertainty is larger than the effect itself (a change in the early
    deficit of ~±1.5–1.8 could not be detected), so this is no evidence that the eras are the same.
  - Deep tremor in 2014–26 shows the same early deficit (rate 0–30 d = 0.32× later, τ ≈ 35 d, 95 % range 11–115 d).
    It just misses the bar (p = 0.014).
- **What Bostock's catalogue can say about stationarity (shape only, inside ETS windows):**
  - The deep groups' share of each ETS is steady at about 19 % (between 14 % and 26 %) over 2003–2013.
  - There is no trend with year (ρ −0.26, p = 0.47).
  - There is no clear tie to how long the previous gap was (ρ +0.53, p = 0.15; 9 episodes).
  - Our 2014–26 big-ETS shares scatter far more (CV 1.0 vs 0.2) because the counts are small. They also show no trend.

## Step 1 — mapping (PASS)
- Key = exact (lat, lon, depth) → `e5r_families_2012.csv` bkey (= our fid).
  - 0 of 269,423 rows unmapped. All 130 Bostock hypocentres are matched.
  - **115/115 production families** receive events: 256,920 events.
- Per-year totals are reproduced exactly (`step1_checks.json`):
  - 2003: 32,978
  - 2004: 59,659
  - 2005: 57,504
  - 2006: 0
  - 2007: 3,147
  - 2008: 14,224
  - 2009: 22,112
  - 2010: 18,114
  - 2011: 18,663
  - 2012: 28,149
  - 2013: 14,873
- Cross-checks:
  - Per family, 2010–2012 agrees with `svi_ref_common.bostock()` (0 of 115 families differ).
  - 2012 per family equals `n_all` of `e5r_families_2012.csv` (0 differ).
- Gaps:
  - **2006 has no events.** By the ~14-month cycle no big SVI ETS is expected in 2006 (2005-09 → 2007-01 is 16 months),
    so this is most likely no episode rather than missing data.
  - **2007 is sparse**: 2,847 LFEs in 11 days, vs 13,600–55,300 for the other episodes. This is an incomplete
    episode or network (Bostock 2015 notes the VGZ change in 2007), not a gap between episodes.
  - Per-year × group counts: `bostock_year_by_group.csv`.
- Coverage (`bostock_scanned_windows.csv`):
  - Events fall on 230 of the 3,883 days between 2003-02-26 and 2013-10-13 (5.9 %).
  - Every day inside the 10 windows has events. No day outside them does.
  - Bostock et al. 2015, §2: "assembled from 10 major ETS episodes between 2003 and 2013". 2003–2005 were the densest
    temporary networks; the templates were then applied to 2007–2013 with subsets of common stations.
  - **Detection ran on ETS windows only, so no between-episode analysis is possible on Bostock.**
  - No scrambled control exists for his catalogue. Time-shift surrogates were not computed, because there are no
    between-episode days to shift into.

## Step 2 — Bostock-era episodes (PASS)
Exposure = the family has ≥ 1 Bostock LFE that UTC day. The final ladder rule (5 %, gaps ≤ 5 d, ≥ 5 d, MAJOR 0.15 → 241.5)
finds 10 episodes, all MAJOR, one per target, with no splits (`episodes_bostock_final.csv`):

| episode | start | end | days | LFEs (115 fams) | max families same day |
|---|---|---|---|---|---|
| 2003-02 | 2003-02-26 | 03-10 | 13 | 31,095 | 93 |
| 2004-07 | 2004-07-02 | 07-27 | 26 | 55,327 | 109 |
| 2005-09 | 2005-09-03 | 09-25 | 23 | 55,316 | 100 |
| 2007-01 | 2007-01-24 | 02-03 | 11 | 2,847 | 80 |
| 2008-05 | 2008-05-10 | 06-02 | 24 | 13,649 | 69 |
| 2009-05 | 2009-05-10 | 06-07 | 29 | 21,507 | 81 |
| 2010-08 | 2010-08-16 | 09-08 | 24 | 17,691 | 81 |
| 2011-08 | 2011-08-19 | 09-09 | 22 | 17,989 | 83 |
| 2012-08 | 2012-08-30 | 09-25 | 27 | 27,445 | 95 |
| 2013-09 | 2013-09-13 | 10-13 | 31 | 14,054 | 101 |

- The frozen first rung (10 %, gap 3, 0.25) gives the same 10 episodes. 2004-07 then starts a day later.
- Recovery is near-automatic, because the catalogue only holds these windows.
- Against the PNSN-tremor list for 2010–2013:
  - Starts agree within 1 day for 2010, 2012 and 2013, and ends within 3 days.
  - For 2011, Bostock's LFEs start 7 days after the tremor (08-19 vs 08-12) and end 6 days later.

## Step 3 — the quiet month (`quiet_month_three_ends.csv`, `figures/p3_quiet_month.png`)
- Statistic and control as in S7/A:
  - F30 = share of between-episode LFEs (P ≥ 0.75, counted once, no tremor within 30 km) within 30 d after the
    previous major ends.
  - Control: the same LFEs on random eligible days of their own family (1000×).
- **Three end definitions:**
  - D1 = the regional end.
  - D2 = end + 7 d, which skips the calm days the episode rule builds in.
  - D3 = the last day that family itself had tremor within 30 km in that major. This is the most local "slip stopped
    here" date that does not use the LFE detector. It is a median of 1 day before the regional end in SVI and 0 days
    in WA.
- **Reproduction:** D1 reproduces S7 for SVI (14 % vs 29 %; S7 printed 27 %, most likely because its shuffle counted
  days before the first major as "not < 30"; p 0.004 vs 0.008) and A for WA (0.323 vs 0.355, p 0.17). Table in the plain answer.
- **Columns, any local tremor allowed (BE_broad):**
  - SVI holds under all three (7 % vs 28 %, p 0.001).
  - WA 2015–26 goes the **opposite** way (44 % vs 36 %): the month after a WA ETS has *more* activity when local
    tremor days are included.
  - WA 2012–14 is the same (76 % vs 39 %). WA 2009–11 holds (9 % vs 23 %).
- Washington item A is final (RESULTS_A.md): it fails its confirmatory test in 2015–26 (p_low 0.17). Descriptively,
  both regions are short of events in the first two weeks (WA 15 vs 37, SVI 1 vs 9), but in WA days 14–30 make up for
  it. **Three-dataset status of the quiet month:** SVI yes (discovery) / Bostock untestable / WA no (main era);
  WA 2009–14 yes (secondary).

## Step 4 — recharge curve (`recharge_fits.csv`, `recharge_curves.csv`, `recharge_era_compare.csv`, `figures/p3_recharge_curves.png`)
- Activity counted: ALL deep LFEs outside majors (P ≥ 0.75, counted once), whether between episodes, in minor
  episodes, or on local-tremor days.
- Rate is per observable family-day. t = days since the previous reference major ended.
- Fit: rate = r∞ (1 − a e^(−t/τ)). Bootstrap over cycles; control = time-shift surrogates (each cycle's daily series
  rotated).

| set | events | cycles | rate 0–30 d / 30 d+ (per 1,000 fam-days) | D30 deficit | surrogate 99th pct | p | τ (95 % range) | recovery present? |
|---|---|---|---|---|---|---|---|---|
| **deep groups 1/4/10/12, clock = 39 majors (primary)** | 1,099 | 39 | 3.4 / 11.9 | 0.72 | 0.82 | 0.077 | 400* (11–400) | **no** (by the rule) |
| families ≥ 40 km, 39 majors | 1,053 | 39 | 4.6 / 16.3 | 0.72 | 0.84 | 0.080 | 400* | no |
| deep groups, end + 7 d | 1,097 | 38 | 5.9 / 12.0 | 0.51 | 0.80 | 0.20 | 400* | no |
| deep groups, family's own end | 1,072 | 39 | 2.1 / 11.8 | 0.82 | — | — | 303 (11–400) | (no control) |
| deep groups, clock = big ETS only | 1,912 | 15 | 5.1 / 16.3 | 0.69 | 0.65 | **0.002** | 322 (85–400) | **yes** (column) |
| group 10 alone, big ETS | 623 | 15 | 1.2 / 28.6 | 0.96 | — | — | 233 (66–400) | (no control) |
| up-dip groups 3/6/8/11 | 25 | 39 | 0.19 / 0.16 | — | — | 0.86 | — | no (too few) |
| scrambled control, deep | 13 | 39 | 0 / 0.16 | — | — | — | — | — |
| PNSN deep tremor 2014–26 (per day) | 11,798 | 40 | 1.2 / 3.7 | 0.68 | 0.69 | 0.014 | 35 (11–115) | no (just misses) |
| PNSN deep tremor 2010–13 | 3,893 | 7 | 2.7 / 3.3 | 0.18 | 0.77 | 0.55 | 400* | no |

\* at the fit's upper limit, i.e. no levelling-off within 400 days.

- **Shape.**
  - With the 39-major clock, the deep rate rises roughly steadily: 0.03× its mean in days 0–7, 0.6× at 14–30 d,
    ~1× at 90–150 d and 2.1× at 240–400 d. A single exponential "recharge" fits poorly; it is a long build-up.
  - The surrogates show a deficit of D30 = 0.33 even with timing scrambled. Longer gaps are busier per day, so part of
    the observed 0.72 is built into the cycle-length mix.
- **Era comparisons (pre-registered rule): both "no detectable difference".**
  - Ours 2014–19 vs 2020–26: τ 400* vs 14 d, overlapping ranges; D30 0.56 vs 0.86, p 0.47.
  - PNSN deep tremor 2010–13 vs 2014–26: D30 0.18 vs 0.68, p 0.24 (end + 7 d: p 0.57).
  - The detectable D30 difference is ±1.5–1.8, which is uninformative. 2010–2013 has only 7 cycles, and 2020–26 has
    almost no long gaps.
- **Bostock-era substitute (deep share inside each ETS, `deep_share_by_episode.csv`):**

  | era | episodes | mean | CV | vs year | vs previous gap |
  |---|---|---|---|---|---|
  | Bostock 2003–13 | 10 | 0.19 | 0.21 | ρ −0.26 (p 0.47) | ρ +0.53 (p 0.15, n 9) |
  | ours 2014–26, big ETS | 14 | 0.24 | 1.04 | ρ −0.14 (p 0.63) | ρ 0.00 (p 1.0) |
  | ours 2014–26, 39 majors (≥ 30 LFEs) | 27 | 0.54 | 0.72 | ρ −0.10 (p 0.61) | ρ +0.24 (p 0.23) |

  - Absolute shares are not compared across eras, because the detector and recall differ.
  - Bostock's network also changed between episodes, so his within-era stability is a bound, not a pure measurement.

## What this means for the paper plan (Merlin Part 3)
- The "23 years" headline cannot rest on Bostock's catalogue for anything between ETS (quiet month, recharge,
  between-episode bursts). What it can carry:
  - within-episode shape observables: tides (step 5, the Royer 2015 positive control is inside his windows), front
    behaviour, and deep share per ETS;
  - our 12 years for the between-ETS behaviour.
- **The quiet month is not a replicated result.** It is one-region (SVI) and fails the confirmatory Washington test.
- **The deep build-up through the inter-ETS gap is suggestive.** The deep rate rises over ~300 days, and it passes the
  surrogate test only with the big-ETS clock. It is a candidate for its own pre-registered test, which could use deep
  tremor in both regions (PNSN 2010–2026). It is not a result yet.
- The idea "the cycle is the same in both eras" is untested: the data cannot resolve it.

## Files
- `PLAN_P3.md`
- Step 1: `step1_checks.json`, `bostock_events_mapped.parquet`, `bostock_scanned_windows.csv`,
  `bostock_year_by_group.csv`
- Step 2: `episodes_bostock_{final,frozen}.csv`, `step2_*`
- Step 3: `quiet_month_three_ends.csv`, `quiet_month_verdicts.json`, `p3_quiet_month.log`
- Step 4: `recharge_fits.csv`, `recharge_curves.csv`, `recharge_era_compare.csv`, `deep_share_by_episode.csv`,
  `recharge_summary.json`, `p3_recharge.log`
- Figures: `svi_science/figures/p3_bostock_timeline.png`, `p3_quiet_month.png`, `p3_recharge_curves.png`
