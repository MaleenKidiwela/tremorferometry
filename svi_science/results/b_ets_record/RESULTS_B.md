# B. ETS over the long record — southern Vancouver Island (SVI) and Washington-NW (2007/2010 → 2026)

Rules: `PLAN_B.md` (frozen 2026-10-07 ~1:25 pm PT; Revision 1 ~1:45 pm PT, before any test). Script
`svi_science/scripts/b_ets_record.py` (build / tests / figures). Times PT. Verdict words as pre-registered: **robust** = p < 0.01
with the primary big-ETS list AND ≥ 2 of the alternatives a1/a2/a4; **suggestive** = p < 0.05 with the primary only;
otherwise **not supported**. Every permutation: 10,000 draws, seed 0.

## 1. Plain-words answer

- **Both regions run on a regular clock, and it is the same clock.** Big ETS are far more regular than random
  (interval CV 0.20–0.61 vs ~0.93 for random timing; every list p ≤ 0.003, **robust**). The main ETS come every ~12–15 months
  in both regions, depending on how you count the smaller pieces. Taking only the largest (post hoc, ≥ 2,000/2,500 epicentres),
  the median is **13.5 months (NW) and 14.3 months (SVI)**, the same as the known ~14 months.
  Over 2010–2026, 13 SVI and NW big ETS pair up (4.7 expected by chance, p = 0.008, **robust**). The delay between them
  does not drift (ρ = 0.01).
- **What differs is the extra pieces.** Washington-NW has extra ETS-sized episodes of its own (2012-04, 2014-08,
  2019-03/04/05, 2022-05, 2024-01, 2025-04b1), and SVI has its own (2019-11, 2021-01, 2022-03). This is why the primary lists
  give a shorter mean wait in NW (6.7 mo) than in SVI (10.7 mo) (p = 0.002, **suggestive**, not robust).
- **No clear trend in size or timing 2009–2026.** ETS sizes (tremor, duration, LFEs, GNSS) do not trend (all p > 0.1).
  The NW intervals shorten in three alternative lists (ρ −0.50 to −0.65, best p = 0.002) but not in the primary list
  (p = 0.07). They are **not supported**, and they could come from the PNSN catalogue detecting more tremor after 2017.
- **A longer wait does not give a bigger ETS (not slip-predictable).** There is a hint that a bigger ETS is followed by a
  longer wait (time-predictable). For SVI, the GNSS offset vs the next interval gives ρ = 0.80, p = 0.002, and NW "cycles" give
  ρ = 0.69–0.89, p ≤ 0.002. These do not hold across the lists, so the hint is **suggestive only**. It comes mostly from
  SVI's 2019–2022 phase: six small ETS (GNSS 1.4–3.8 mm), mostly followed by short waits (3–12 months). In the phases before and
  after, ETS were ~5 mm and came every 13–16 months.
- **No seasonality.** ETS onset months are not clustered (Rayleigh p 0.11–0.94). Tremor between episodes shows no annual
  cycle beyond the scrambled control (SVI p 0.12, NW p 0.69). SVI between-episode LFEs show none either (p 0.86).
- **Neither side of the border reliably goes first.** Washington starts first in 8 of 13 pairs, SVI in 5 (binomial p 0.58).
  The median delay is 4 days, Washington first (90 % CI −3 to +8 d), and |delay| is 1–14 days. The order PNSN tremor gives
  matches the S3 LFE order in sign in all 6 pairs S3 measured. It also matches S3's front types: WA→SVI sweeps are
  Washington-first; SVI→WA sweeps and Strait-start splits are SVI-first or about equal.

## 2. Data, footprint, big-ETS lists

- Episodes: `episodes_nw_final.csv` (118; three WECC windows 2007-01, 2008-04, 2009-04 + PNSN 2009-08 → 2026-05) and
  `episodes_svi_final.csv` (79, 2014+). The SVI list was extended to 2010 with the identical code and rule (21 more episodes,
  7 major). The rerun reproduces every 2014+ episode exactly (`build_info.json`).
- Tremor footprint: PNSN epicentres ≤ 10 km from a family (SVI 115 Bostock families, NW 132 families). 9,032 epicentres
  within 10 km of both regions were given to the region of the nearest family. In total, SVI has 85,289 and NW 73,323 epicentres.
- LFE columns: SVI = S3 counted-once P ≥ 0.75. NW = catalogue_v2 `selected_v2`, counted once (dedup D), recomputed for all eras.
  For the 47 E3 majors it matches S3 exactly (max relative difference 0). Size tests use NW LFE counts from E3 only.
  GNSS = S1 Tier A `D` (SVI, 2014+).
- **Big ETS.** Primary = optimal two-class split of log10(tremor count) among MAJOR episodes. This gives SVI 18 (threshold
  1,171 epicentres) and NW 28 + 3 WECC (threshold 655). Alternatives: a1 = family-exposure ≥ half the families × 14 d
  (SVI 17, NW 21 + 3); a2 = the 16 largest, one per year (16 + 3 WECC); a3 = S3 LFE "big" (SVI 10, NW 14, 2014/15+ only,
  shown only); a4 = primary with onsets < 120 d apart merged into one cycle (SVI 15, NW 20 + 3).
  Revision 1: the plan's first primary rule split over ALL episodes. It only separated tiny bursts, giving 51/55 "big",
  about 3 per year. It was replaced before any test.

| Big ETS (primary) | onset PT | tremor | LFE once | GNSS D mm |
|---|---|---|---|---|
| SVI | 2010-08-20, 2011-08-19, 2012-09-06, 2013-09-15, 2014-11-18, 2015-12-23, 2017-03-21, 2018-06-25, 2019-08-12, 2019-11-08, 2020-11-06, 2021-02-02, 2021-10-04, 2022-03-27, 2023-04-05, 2024-08-16, 2025-08-17, 2025-10-01 | 1,171–5,911 | 42–1,498 (2017-03 and 2025-09b1 low = known station problems) | 4.3, 5.9, 5.6, 5.5, 1.6, 2.2, 1.9, 3.8, 1.9, 1.4, 5.0, 5.1, 0.0, 5.1 (2014-11 on) |
| NW | WECC 2007-01, 2008-04, 2009-04; then episodes (onset PT for the largest): 2010-08 (08-08), 2011-07 (2011-08-06), 2012-04, 2012-08, 2012-09 (09-19), 2013-09 (09-07), 2014-08, 2014-11 (11-22), 2015-07, 2015-12 (2016-01-02), 2017-02 (02-18), 2017-03, 2018-05 (05-02), 2019-03, 2019-04, 2019-05, 2019-08, 2020-05, 2020-09 (2020-10-04), 2020-10, 2021-09 (09-27), 2022-05, 2022-10b1, 2023-03 (2023-04-09), 2024-01, 2024-08 (08-09), 2025-04b1, 2025-09 (09-18) | 655–4,128 | E3: 211–1,188 | — |

Full table: `episodes_b.csv` (every listed episode, all columns, all big flags).

## 3. Q1 — intervals (`intervals.csv`, `tests.json` Q1)

| list | n | mean (mo) | median (d) | CV | random-timing CV (median / 1st pct) | p (CV ≤ random) | drift ρ (p) | before / after 2018 (d), MW p |
|---|---|---|---|---|---|---|---|---|
| SVI primary | 18 | 10.7 | 374 | 0.44 | 0.93 / 0.58 | < 1e-4 | −0.31 (0.23) | 409 / 266, 0.055 |
| SVI a1 | 17 | 11.3 | 379 | 0.37 | 0.93 / 0.57 | < 1e-4 | −0.06 (0.84) | 409 / 295, 0.17 |
| SVI a2 | 16 | 12.1 | 375 | 0.49 | 0.92 / 0.56 | 0.003 | −0.10 (0.71) | 478 / 295, 0.15 |
| SVI a4 cycles | 15 | 12.8 | 392 | 0.20 | 0.92 / 0.55 | < 1e-4 | 0.00 (1.0) | 409 / 373, 0.54 |
| NW primary + WECC | 31 | 7.5 | 244 | 0.61 | 0.95 / 0.67 | 0.002 | −0.35 (0.065) | 275 / 180, 0.068 |
| NW primary, PNSN only | 28 | 6.7 | 220 | 0.61 | 0.95 / 0.65 | 0.003 | −0.12 (0.57) | 235 / 180, 0.29 |
| NW a1 + WECC | 24 | 9.7 | 313 | 0.45 | 0.94 / 0.63 | < 1e-4 | −0.58 (0.005) | 344 / 245, 0.053 |
| NW a2 + WECC | 19 | 12.5 | 384 | 0.30 | 0.93 / 0.60 | < 1e-4 | −0.50 (0.035) | 413 / 337, 0.10 |
| NW a4 + WECC | 23 | 10.2 | 328 | 0.39 | 0.94 / 0.63 | < 1e-4 | −0.65 (0.002) | 375 / 245, 0.013 |

- Regularity: **robust** in both regions.
- Drift: SVI **not supported** under any list. NW **not supported**: the primary list has p 0.065. Two alternatives pass
  (a1 0.005, a4 0.002); a2 does not (0.035). Caveat: footprint tremor per year rose from ~1,500–3,800 (2015–16) to ~5,500–6,700
  (2024–25). Size-threshold lists therefore admit more mid-size episodes in later years, which shortens intervals without any
  change at the fault.
- Post hoc, not tested: the 13 largest per region (NW ≥ 2,000, SVI ≥ 2,500 epicentres) recur at median 13.5 / 14.3 months
  (mean 15.1 both, CV 0.32 / 0.27). Each region has one long gap: NW 2018-05 → 2020-10, 29 months, which held the split 2019
  episodes. SVI 2013-09 → 2015-12, 27 months, because the 2014-11 SVI episode had only 1,489 epicentres even though its
  GNSS offset was 4.3 mm.

## 4. Q2 — size vs time and vs the wait (`tests.json` Q2; figure `b_interval_size.png`)

Spearman ρ (permutation p); "before" = slip-predictable, "after" = time-predictable; detrended = sizes with the year trend removed.

| list / size | trend with year | vs interval before | vs interval after | after, detrended |
|---|---|---|---|---|
| SVI primary, tremor | 0.01 (0.97) | 0.11 (0.67) | 0.47 (0.054) | 0.47 (0.060) |
| SVI primary, GNSS D | −0.35 (0.22) | 0.39 (0.17) | **0.80 (0.002)** | 0.68 (0.013) |
| SVI a1, GNSS D | −0.24 (0.44) | 0.56 (0.050) | **0.81 (0.002)** | 0.79 (0.004) |
| SVI a2, GNSS D | −0.30 (0.34) | 0.67 (0.020) | 0.75 (0.012) | 0.79 (0.005) |
| SVI a4, GNSS D | −0.29 (0.38) | 0.53 (0.10) | 0.30 (0.40) | 0.25 (0.49) |
| SVI primary, LFE | −0.09 (0.77) | 0.33 (0.25) | 0.21 (0.49) | 0.21 (0.48) |
| NW primary, tremor | 0.05 (0.82) | 0.29 (0.13) | 0.34 (0.084) | 0.34 (0.082) |
| NW a4, tremor | −0.21 (0.37) | −0.01 (0.98) | **0.69 (0.002)** | 0.60 (0.008) |
| NW a4, LFE (E3) | −0.04 (0.90) | −0.23 (0.42) | **0.89 (< 1e-4)** | 0.89 (< 1e-4) |
| NW a2, tremor | −0.15 (0.57) | −0.06 (0.83) | 0.18 (0.51) | 0.13 (0.64) |

- Sizes do not trend with time in either region (all |ρ| ≤ 0.48, p > 0.03; none pass) → **not supported**.
- Slip-predictable (longer wait → bigger ETS): **not supported** (best SVI a2 GNSS 0.67, p 0.02; primary 0.39, p 0.17).
- Time-predictable (bigger ETS → longer wait): SVI GNSS passes with the primary list and a1, fails a2 by a hair (0.012), and
  fails a4 → **suggestive**. NW passes only with a4 cycles → **not supported** by the rule. Note that merging a cycle's pieces
  makes its size larger and its onset earlier, which by itself favours a positive ρ with the next interval.
- What drives the SVI hint: 2019-08 → 2022-03 had six ETS of 1.4–3.8 mm, 3–12 months apart. Before (2014-11 → 2018-06) and
  after (2023-04 → 2025-09) they were 4.3–5.9 mm and 13–16 months apart. Post hoc, the mean GNSS offset released per year is
  similar in the three phases: ~4.4, ~3.8 and ~4.1 mm/yr (sum of D / phase length). The total release stays about the same;
  it was split into smaller, more frequent pieces in 2019–2022.

## 5. Q3 — seasonality (`tests.json` Q3, `seasonal_monthly.csv`; figure `b_seasonal.png`)

- Onsets (Rayleigh / interval-shuffle scrambled control): SVI primary R 0.35, p 0.11 / 0.20 (mean onset ~late Sep; many
  Aug–Sep onsets in 2010–2013 then drifting through the year). NW primary R 0.21, p 0.30 / 0.44. Every alternative p ≥ 0.11
  → **not supported**.
- Between-episode tremor (rate per between-day, annual harmonic, each year's months rotated as the scrambled control):
  SVI amp/mean 0.56, peak ~Aug, p 0.115 (4,077 between-days, 6,497 epicentres). NW 0.18, p 0.69 (3,767 days, 6,131).
  Share of between-days with any tremor: SVI p 0.31, NW p 0.57. SVI between-episode LFEs (S7, 74 in these days):
  p 0.86 → **not supported**.
- Detection caveat (diagnostic): epicentres per day on big-ETS days show no significant annual cycle (SVI p 0.88, NW p 0.18;
  NW peaks in Dec/Feb, the opposite of what winter station noise would do). The PNSN record therefore gives no sign that
  seasonal noise distorts these numbers, but with 270–370 ETS days per region this check is weak.

## 6. Q4 / Q5 — cross-border order and clocks (`pairs.csv`, `tests.json` Q4/Q5; figure `b_crossborder.png`)

Primary pairs (delay = SVI onset − NW onset; > 0 = Washington first; 5th-percentile tremor onset, PT):

| SVI | NW | delay d | 10th-pct delay | S3 LFE NW lead d | S3 fronts (SVI / NW) |
|---|---|---|---|---|---|
| 2010-08b1 | 2010-08 | +11.9 | +11.7 | — | — |
| 2011-08 | 2011-07 | +13.0 | +12.5 | — | — |
| 2012-08 | 2012-09 | −13.8 | −13.7 | — | — |
| 2013-09 | 2013-09 | +7.8 | +10.0 | — | — |
| 2014-11 | 2014-11 | −4.6 | −5.1 | — | — |
| 2015-12 | 2015-12 | −10.0 | −9.6 | −5.8 | SE / SE (SVI→WA sweep) |
| 2017-03 | 2017-03 | +3.9 | +4.4 | — | — / no single front |
| 2019-08 | 2019-08 | −1.4 | −1.0 | −3.8 | NW / SE (Strait split) |
| 2020-11 | 2020-10 | +4.0 | +4.3 | +4.3 | NW / none |
| 2021-10 | 2021-09 | +7.7 | +6.6 | +8.0 | — / NW |
| 2023-04 | 2023-03 | −3.3 | −3.6 | −3.2 | NW / SE (Strait split) |
| 2024-08 | 2024-08 | +7.4 | +7.0 | +4.0 | NW / NW (WA→SVI sweep) |
| 2025-09b1 | 2025-09 | +12.5 | +12.8 | — | — / NW |

- Pairing beyond chance: 13 pairs vs 4.7 expected for randomly shifted SVI onsets (p 0.008). a1 10 vs 3.5 (p 0.007),
  a2 9 vs 2.6 (p 0.010), a4 7 vs 3.0 (p 0.041) → **robust**.
- Order: Washington first 8 / 13 (primary), 6 / 10 (a1), 6 / 9 (a2), 10 / 13 (a4, binomial p 0.09) → **no consistent order**.
  The median delay is +4.0 d (90 % CI −3.3 … +7.8). With a4 cycles it is +11.9 d (CI 7.4 … 32.8), because a cycle's onset
  is its first Washington piece (e.g. NW 2020-09 starts 33 d before SVI 2020-11).
- Unpaired big ETS: SVI-only 2018-06 (NW 2018-05 ended ~5 weeks before), 2019-11, 2021-01, 2022-03, 2025-08. NW-only
  2012-04, 2012-08, 2014-08, 2015-07, 2017-02 (the first piece of the 2017 pair), 2018-05, 2019-03/04/05, 2020-05, 2020-09,
  2022-05, 2022-10b1, 2024-01, 2025-04b1.
- Q5 relative clocks: the delay does not trend with year (primary ρ 0.01, p 0.99; a1 −0.09, a2 −0.13, a4 −0.47, p ≥ 0.11)
  → **no relative drift**. Mean intervals differ (primary SVI 325 d vs NW 204 d, p 0.002; a4 p 0.005; a1 p 0.047; a2
  identical by construction) → **suggestive**: same main beat, more extra pieces in Washington.

## 7. What is known vs new

- Known: ~14-month ETS recurrence under southern Vancouver Island / northern Washington (Rogers & Dragert 2003). Recurrence
  differs between along-strike segments (Brudzinski & Allen 2007). The 2010–2026 tremor record reproduces the ~14 months
  for the largest episodes in both regions. Lowry (2006) proposed a seasonal modulation of slow slip by hydrological loading;
  the claim has been debated since (we did not re-check the later literature here). This record shows no seasonal phase
  in onsets or in between-episode tremor. PNSN detection depends on seasonal station noise; our check found no sign of it
  in ETS-day rates.
- New here, from one detection-independent footprint for both regions: (i) the two regions share the main beat with no
  relative drift 2010–2026, but each has its own extra mid-size episodes; (ii) the cross-border order changes from episode to
  episode (|delay| 1–14 d), and the PNSN order matches S3's LFE front types; (iii) SVI 2019–2022: six smaller, more frequent
  ETS with about the same long-term release rate, and a suggestive time-predictable link between GNSS size and the next wait.
  None of (iii) is robust by the pre-registered rule, and n is 13–18.

## 8. Files
`PLAN_B.md`, `build_info.json`, `episodes_b.csv`, `tremor_footprints.parquet`, `intervals.csv`, `pairs.csv`,
`seasonal_monthly.csv`, `tests.json`, `b.log`; figures `svi_science/figures/b_timeline.png`, `b_interval_size.png`,
`b_seasonal.png`, `b_crossborder.png`.
