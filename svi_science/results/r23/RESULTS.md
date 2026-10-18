# R23 — five findings re-tested on the 23-year SVI LFE record (2003–2026)

Written Sun 2026-10-11, ~12:15–12:40 am PT. Rules: `PLAN.md` (pre-registered before any number; Revisions 1–3 at its end).
Scripts `svi_science/scripts/r23_{base,common,front,bursts,quiet,cycle,ramp,posthoc}.py`; tables `tables/`, figures
`figures/`. Science only: no catalogue file was changed. "Scrambled control" = the same detector on the same day with
each station site shifted to another time of that day. Times are PT.

## Verdicts at a glance

| study | 2003–13 (new) | 2014–26 (before) | 2003–26 together | verdict |
|---|---|---|---|---|
| 1 Fronts stop where the fault was recently active (SVI LFE check) | 1 stop in 10 big ETS | 8 stops | 9 stops, β +0.24 (p 0.41) | **undecided** — still too few stops |
| 2 Deep between-episode bursts | exist, deep, bursty, same patches | yes | yes | **supported** (sizes: no detectable difference) |
| 3 Deep quiet 1–2 months after big ETS | same shape; p 0.059 (0.009 with tremor-based ends, diagnostic) | passes | passes (p 0.001) | **undecided** (pre-registered); leans towards yes |
| 4a Regular ETS cycle | CV 0.13 | CV 0.44 | CV 0.32 | **supported** |
| 4b Bigger episode → longer wait (LFE sizes) | ρ −0.07 | ρ +0.38 | ρ +0.31 (p 0.17) | **not supported** |
| 5 Deep LFEs ramp up before ETS | no | no | no | **not supported** |

## Foundation (what every study uses)

- **Events.** Single repeats on the 115 production families:
  - 2003–13 from v1.4-backward: primary detector, counted once, era rows at P ≥ 0.90, frozen rows at P ≥ 0.75. This
    gives 34,093 data events and 85 scrambled control over all families; per year it matches the P\* table of
    RESULTS_BACKWARD exactly.
  - 2014–26 from the forward catalogue (P ≥ 0.75, counted once).
  - Network segments carried on every row: era 2003-02 → 2006-06, frozen B2 / B3 → 2009-08, frozen S15 → 2013, forward.
- **Tremor.** PNSN from 2010-01-06. Before that, the frozen envelope tremor detector on the 4 CN stations (one fixed
  network for 2003–13). PNSN distances reproduce the S7 exposure on all but 8 of 60,076 family-days.
- **Episodes 2003–2026: 145 listed, 52 majors, 23 big ETS.**
  - 2010–26 use the PNSN lists unchanged.
  - Before 2010 the envelope-based list failed its validation (Revision 1). One variant merged 2003-03 → 2010-01 into a
    single episode. The other found 23 "majors" in 2010–13 where PNSN has 7.
  - So the pre-registered fallback is used: the 6 Bostock ETS windows are the majors (2003-02, 2004-07, 2005-09,
    2007-01, 2008-05, 2009-05), and 39 envelope tremor episodes are the minors. That is 5–8 listed episodes per year,
    as in every other era.
  - Big ETS = the upper class of episode size in Bostock-equivalent LFEs (Σ `rate`, split at 1,640). That gives 10 big
    ETS in 2003–13 (the 6 Bostock windows plus 2010-08, 2011-08, 2012-08, 2013-09) and 13 in 2014–26.
- **Caveat that matters below.** The Bostock-window majors end where his scan ended. In 2003-02 and 2008-05 that is
  before the tremor stopped (envelope tremor to 03-18 vs 03-10; 06-09 vs 06-02). The start of 2008-05 is also late
  (tremor from 04-30). Events in those tails count as "between episodes". Study 3 therefore has a diagnostic that
  widens those majors to the tremor episode (`--union`).

---

## 1. Front arrest (finding 8): do fronts stop where the fault was active in the previous month?

**Undecided: adding 2003–13 does not help.** Big-ETS fronts in the LFE record almost never stop inside the
~100-km SVI family box. They sweep through it. The 10 new big ETS add **one** stop (2007-01, at −55 km). Joint
2003–2026 has 9 stops, below the 20 the plan required for any verdict. The whole-margin PNSN tremor result (finding 8,
p 0.0005) is neither confirmed nor contradicted by SVI LFEs.

Numbers (`tables/study1/`; walk = continue/stop at the next 5-km bin; Q6 code and order shuffle unchanged; shuffle
within network segment):

| targets, period | stops | β(log T_any) | p (one-sided, > 0) | β(log T_big) | p |
|---|---|---|---|---|---|
| **big ETS, 2003–26 (headline)** | **9** | **+0.24** | **0.41** | +0.62 | 0.022 |
| big ETS, 2003–13 | 1 | — (1 stop) | 0.62 | +0.64 | 0.33 |
| big ETS, 2014–26 | 8 | +0.06 | 0.42 | +0.68 | 0.11 |
| all 47 majors, 2003–26 (column) | 42 | −0.13 | 0.54 | −0.04 | 0.46 |
| S3's 10 big, 2014–26 (reproduction of Q6) | 6 | −0.35 | 0.61 | +0.17 | 0.48 |

- The reproduction matches Q6 exactly: 6 stops, β −0.35, enter β +0.30.
- Stop share by time since the bin's last LFE activity (big ETS, 2003–26): < 30 d 0 of 2, 30–90 d 3 of 44, > 90 d 6
  of 264. Fronts almost never reach a bin that was active in the previous month (2 of 310 steps), so the question
  cannot be asked of SVI LFEs.
- Columns, not verdicts:
  - **Enter model, joint, big ETS:** bins that had a big ETS longer ago are entered more often (β +0.83, p 0.0015;
    T_any +0.45, p 0.016). This is the same direction as the tremor result. In 2014–26 alone it gives +0.83, p 0.005.
  - **Tremor-history version** (PNSN / envelope bin-days ≥ 3): β −1.23, p (negative) 0.048. It points the other way and
    does not pass.
- **Network caveat:** the era detector (2003–06) marks more bins active. The shuffle is within segment, and the
  propensity covariate is the Q6 one.
- **New vs 2014–26 only:** the 2003–13 fronts behave the same way, sweeping the box. The SVI LFE check stays without
  power: 9 stops instead of 6.

## 2. Between-episode deep bursts (finding 3)

**Supported: they exist in 2003–13, in the same deep patches, with burst sizes not detectably different.**
Between ETS, with no tremor within 30 km, the LFEs come from deep families and cluster on a few days, as in 2014–26.

| period / segment | tremor-free between-episode LFEs (scrambled control) | broad, local tremor allowed (control) | mean depth vs ETS-share expectation | distinct days vs shuffle |
|---|---|---|---|---|
| **2003–13** | **109 (9)** | 1,280 (14) | **39.0 vs 35.0 km, p < 0.0001** | **68 vs 106, p 0.001** |
| era 2003-02 → 2006-06 | 24 (2) | 857 (2) | 38.9 vs 35.0 km, p < 0.0001 | 15 vs 24, p 0.001 |
| frozen B2/B3 2006-07 → 2009-08 | 8 (2) | 147 (4) | 33.8 vs 34.6 (n.s.) | 8 vs 8 |
| frozen S15 2009-08 → 2013 | 77 (5) | 276 (8) | 39.5 vs 35.0 km, p < 0.0001 | 45 vs 75, p 0.001 |
| 2014–26 (S7, reproduced exactly: 76 / 3 / 688) | 76 (3) | 688 (11) | 39.2 vs 35.4 km, p < 0.0001 | 44 vs 75, p 0.001 |

- **Same patches.** The group shares of tremor-free between-episode LFEs agree between the decades: Spearman ρ 0.74
  over the 12 groups (group-label permutation p 0.003). The excess over the ETS share agrees too (ρ 0.71, p 0.007). The
  broad set agrees almost perfectly (ρ 0.93, p 0.0001).
  - The deep groups dominate in both decades. 54 % of the 2003–13 events come from families ≥ 40 km deep, against 16 %
    of ETS events.
  - The mix shifts between decades. In 2003–13, G4 (fam 87, 88; 41–45 km, central down-dip) holds 29 % and G10 (NW,
    44 km) 15 %. In 2014–26 these are 8 % and 30 %. The busiest ETS group, G2, has 2 % and 0 %.
  - Top families: 2003–13 fam 87 (18 events), 88 (14), 119 (7); 2014–26 fam 119 (9), 114 (6), 33 (4).
- **Burst sizes** (bursts = events chained by gaps < 24 h; `tables/study2/burst_size_compare.csv`):
  - The median burst is 1 event in every segment. The share of events in bursts of ≥ 3 is 54 % (era), 43 % (S15) and
    54 % (forward). Durations, families per burst, separations and depths show no detectable difference
    (S15 vs forward: all p ≥ 0.2).
  - Only the Bostock-equivalent size of era-detector bursts is smaller (median 2.3 vs 5.5, p 0.0001). The era detector
    picks up more small ones per Bostock-equivalent unit, so this is a detector difference.
  - The largest tremor-free burst of 2003–13 is 2013-04-21 17:02–23:26 PT: 21 LFEs from fam 65/87/88, ~45 km. In
    2014–26 the largest is 11.
  - Post-hoc, broad bursts (local tremor allowed) are the same size in both decades. The largest Bostock-equivalent
    bursts are 2003-06-12 (G10, 224), 2013-06-26 (G10/12, 188), 2006-05-25 (187), 2007-05-19 (G10, 168) and
    2008-10-22 (G1, 163). In 2014–26 they are 2020-09-12 (404), 2015-10-07 (291) and 2014-03-22 (162).
- **Early-2006 bursts** (`tables/study2/bursts_2006.csv`; counted once, era P\* ≥ 0.90; none has a scrambled-control
  event):
  - 2006-01-14 17:58 → 01-16 05:05 PT: 91 LFEs, G10 (fam 116/112/114/110/119), 42 km.
  - 2006-05-25/26: 121 LFEs, G9/10/12, 42 km.
  - 2006-06-17/18: 33 LFEs, G12 (fam 128/125/126), 41 km.

  All three are deep between-episode bursts *with* local envelope tremor (broad set, not tremor-free): the same patch,
  G10, as the main 2014–26 bursts. The 2006-04-24 (G4, 60 LFEs), 05-01 (G1/G2, 203), 05-03 (G2, 122) and 05-11 (G10, 5)
  bursts fall inside a listed envelope tremor episode (2006-04-23 → 05-13), so they count as a minor episode, not as
  between-episode. **2006-03-27…30: only 2 LFEs at P\* ≥ 0.90.** The new catalogue does not support it as a burst.
- **Rates are not compared across segments** (usage rule). Within-segment counts per 1,000 eligible family-days are in
  `inventory.csv` for reference only.
- **Columns:**
  - Without burst-day rows: 74 of 109 tremor-free events in 2003–13 remain. The forward catalogue has no burst-day flag.
  - Rows with n_z2 ≥ 2: 108 of 109.
  - Envelope tremor as the source for 2010–13 too: 71.
- **Network caveat:** B2/B3 (2006-07 → 2009-08) sees only 8 tremor-free events. That network is the least sensitive
  between episodes, so this is no evidence of an absence.
- **New vs 2014–26 only:** the deep between-ETS LFE bursts are a 23-year feature, not a property of the forward
  network. The same deep patches (G10, G4, G1, G12) fire in both decades. The share between the two deepest clusters
  changes between decades.

## 3. Deep quiet after big ETS (finding 4, provisional): confirm or reject the 1–2-month quiet?

**Undecided by the pre-registered rule, but leaning towards confirmed.** In 2003–13 the deep families are almost
silent for about two months after a big ETS, as in 2014–26. Neither pre-registered test clears its bar in 2003–13
alone:
- the burst-driven control is wide;
- the tremor-free all-family test is spoiled by shallow ETS-tail events after the Bostock window ends.

With both decades together, the recharge test passes clearly.

**(B) Deep-group activity vs days since the previous big ETS ended** (P3 recharge; rate per deep family-day,
normalised by the fitted later level; `tables/study3/recharge_curves.csv`):

| days after the end | 0–7 | 7–14 | 14–30 | 30–60 | 60–90 | 90–150 | 150–240 | 240–400 |
|---|---|---|---|---|---|---|---|---|
| 2003–13 (10 big ETS) | 0.13 | 0.12 | 0.03 | 0.01 | 0.43 | 0.34 | 0.69 | 0.69 |
| 2014–26 (13 big ETS) | 0.01 | 0.01 | 0.00 | 0.06 | 0.14 | 0.23 | 0.43 | 0.55 |
| envelope deep tremor 2003–13 (fixed CN network) | 0.98 | 0.61 | 0.29 | 0.26 | 0.73 | 0.80 | 0.78 | 0.96 |

| test | 2003–13 | 2014–26 | 2003–26 |
|---|---|---|---|
| D30 (deficit of the first 30 d) | 0.86 (95 % 0.79–0.98) | 0.98 (0.96–1.0) | 0.90 (0.83–0.99) |
| time-shift control p (bar < 0.01) | **0.059** | 0.002 | **0.001** |
| recovery time τ | 241 d (134–400) | 400 d (pinned; 112–400) | 321 d (166–400) |
| D30, days counted instead of events | 0.49 (p 0.11) | 0.84 (p 0.005) | 0.67 (p 0.006) |
| D2 (end + 7 d) | 0.92 (p 0.039) | 0.98 (p 0.003) | 0.94 (p 0.001) |

**(A) S7 quiet month** (tremor-free between-episode LFEs of all families, share within 30 d after the previous major;
`tables/study3/quiet_month_F30.csv`):
- 2003–13: 0.20 vs 0.13 for shuffled days. That is the wrong side (p_low 0.98); the big-ETS clock gives 0.17 vs 0.09.
- 2014–26: 0.13 vs 0.28 (p 0.002), which reproduces S7.
- 2003–26: 0.17 vs 0.19 (p 0.29).
- Why 2003–13 fails: the early events are mostly *shallow* families right after a Bostock window end. Examples:
  2003-03-18 (6 LFEs of G11/G12, 8 d after his window but inside the envelope tremor episode), 2004-08-09/11/23,
  2008-06-18/23, 2012-10-03.
- Deep families only (post-hoc, Revision 3; `tables/posthoc/deep_quiet_month_F30.csv`):
  - Deep between-episode LFEs (local tremor allowed) in 2003–13: F30 0.069 vs 0.107 (p 0.001), F60 0.076 vs 0.220
    (p 0.001). With the big-ETS clock: 0.050 vs 0.094 (p 0.001).
  - Tremor-free deep only: F30 0.13 vs 0.13 (p 0.49), F60 0.15 vs 0.28 (p 0.024).
- Widening the Bostock majors to the tremor episode (`--union`, diagnostic, `tables/study3_union/`):
  - Tremor-free, all families: 0.16 vs 0.13 (still not passed).
  - Broad: 0.043 vs 0.109 (p 0.001). The deep broad set: 0.033 vs 0.108 (p 0.001).
  - **Recharge (B) passes:** D30 0.96, time-shift p 0.009, a 95 % CI 0.996–1.0, τ 240 d (130–400). D2 gives 0.97
    (p 0.003) and day-count 0.70 (p 0.021).
  - So with tremor-based episode ends, 2003–13 would be "partly": (B) yes, (A) no.

**Pre-registered verdict for 2003–13:** (A) fails on the wrong side and (B) misses p < 0.01 (0.059) with D30 > 0, so
the verdict is **undecided**. In plain words, the deep zone in 2003–13 shows the same post-ETS gap in its point
estimates: 1–13 % of the later rate for ~60 days, and an envelope-tremor dip at 14–60 d (0.26–0.29, not significant,
p 0.21). But 10 cycles of very bursty deep activity cannot reach the bar alone. Joint 2003–26 passes (p 0.001).
When the Bostock-window ends are moved to where the tremor stopped (diagnostic), 2003–13 passes alone (p 0.009).
- **Network caveat:** the test compares days inside one cycle, so slow network drifts cancel. Cycles crossing a segment
  boundary were cut and shifted separately. The 2005-09 → 2007-01 cycle crosses the era → frozen switch at day ~280,
  after the window that matters.
- **New vs 2014–26 only:**
  - The quiet first two months after big ETS are seen in a second, independent decade, at the same point-estimate
    strength.
  - The S7 version (tremor-free LFEs of all families) does not replicate, because shallow ETS tails dominate it.
  - The finding should be stated for the **deep families**, not for all between-episode LFEs.
  - The recovery time is long (τ ~ 240–320 d, fits pinned near the 400-d limit). This is not the ~48-d recovery of
    PNSN deep tremor.

## 4. Cycle regularity and "bigger episode → longer wait"

**(a) Regular cycle: supported, and 2003–13 is even more regular.** Big ETS under SVI recur every 12–17 months in
2003–13: intervals 362–508 d, mean 428 d = 14.1 months, CV 0.13 against 0.90 for random timing (p 0.0001). The
2005-09 → 2007-01 gap is 508 d (16.7 months), the longest of the record. 2014–26 is less regular (CV 0.44, p 0.002),
because of the short 2019–22 intervals of 3, 3, 8 and 5 months. Together: CV 0.32 (p 0.0001). The drift towards shorter
intervals over 2003–26 (ρ −0.43, p 0.048) misses the bar. It comes from 2019–22.

**(b) Bigger → longer wait: not supported with our own LFE sizes.**

| size measure | period | ρ(size, wait after) | p | detrended |
|---|---|---|---|---|
| Σ rate (Bostock-equivalent) | 2003–13 | −0.07 | 0.84 | −0.39 (0.27) |
| Σ rate | 2014–26 | +0.38 | 0.22 | +0.43 (0.16) |
| **Σ rate** | **2003–26** | **+0.31** | **0.17** | +0.15 (0.51) |
| Bostock's own counts | 2003–13 | +0.36 | 0.30 | — |
| GNSS D (S1), our big list | 2014–26 | +0.81 | 0.002 | +0.76 (0.006) |
| GNSS D, B's list (reproduction of ρ 0.80) | 2014–26 | +0.80 | 0.002 | +0.66 (0.014) |
| all 52 majors, Σ rate (column) | 2003–26 | +0.70 | 0.0001 | +0.43 (0.002) |

- The GNSS-size result is reproduced (ρ 0.80), but the LFE sizes of the same 2014–26 episodes do not show it (ρ 0.38).
  It does not appear in 2003–13 at all. For example, the largest episode, 2004-07 (49,800), was followed by 14 months,
  and 2007-01, one of the smallest (5,400), by 15.5.
- The all-majors column mixes small majors that come just before a big one. It is a property of the list, not the
  time-predictable test.
- **Network caveat:** Σ rate is the only size allowed across networks. Pre-2006 sizes come from the era detector, and
  the ETS-day k is uncertain by ×2–3 (P7 brackets). Rank tests within each decade do not depend on it.
- **New vs 2014–26 only:** the ~14-month cycle holds over 23 years and was very regular in 2003–13. The "bigger →
  longer wait" hint is a GNSS-only, 2014–26-only result that our own LFE sizes do not reproduce.

## 5. Pre-ETS ramp-up

**Not supported in any period.** Deep families do not become more active in the weeks before an ETS starts.

| period | deep between-episode LFEs in the last 30 d before a major (F30) | shuffled days | p | ≥ 30 d after previous end: F30 (p) | scrambled pre-window position: rate per 1,000 family-days (p) |
|---|---|---|---|---|---|
| 2003–13 | 0.085 (n 1,128) | 0.096 | 0.89 | 0.071 (0.95) | 9.8 vs 9.6 (0.49) |
| 2014–26 | 0.25 (n 612) | 0.27 | 0.81 | 0.24 (0.97) | 5.5 vs 6.2 (0.64) |
| 2003–26 | 0.14 (n 1,740) | 0.16 | 0.95 | 0.13 (0.99) | 6.6 vs 7.1 (0.57) |

- The same holds with the big-ETS clock, for tremor-free events, for deep tremor-free events, burst by burst
  (`tables/study5/ramp_F30.csv`), and in deep tremor (scrambled pre-window p 0.74 / 0.19 / 0.49).
- One column hint: in 2003–13, F14 (last 14 days) in the ≥ 30-d version reaches 4.9 % vs 3.3 % (p 0.006). It is a
  single burst, 2008-05-06/07 (49 LFEs). That burst is inside the 2008-05 ETS by tremor (envelope episode from 04-30)
  and only before Bostock's window start (05-10), so it is a window-start artefact, not a ramp.
- In 2014–26 the last two weeks before a major are, if anything, emptier than random (2 % vs 12 %).
- **Network caveat:** the shuffles stay within network segment, and control 2 stays inside the same inter-ETS gap.
- **New vs 2014–26 only:** the S7 "no ramp" result now holds on 23 years and in a fixed-network tremor check.

## Files
- `tables/events_2003_2026.parquet`, `famdays_observable.parquet`, `famday_tremor.parquet`, `tremor_points.parquet`,
  `episodes_2003_2026.csv` (sizes, big flags, kinds), `base_summary.json`, `env_episodes_2010_2013_variant_{a,b}.csv`.
- `tables/study1/` (front_arrest_table.csv, walk rows, summary), `study2/` (inventory, bursts, burst_size_compare,
  group_shares, bursts_2006, between-episode event lists 2003–13), `study3/` and `study3_union/` (quiet_month_F30,
  recharge_fits, recharge_curves), `study4/` (cycle_tests, big_ets_2003_2026), `study5/` (ramp_F30,
  ramp_scrambled_dates, superposed_epoch), `posthoc/`.
- Figures:
  - `figures/r23_between_timeline.png`: 2003–26 between-episode LFEs and group shares;
  - `r23_between_burst_sizes.png`;
  - `r23_quiet_after_ets.png` (and `_union`);
  - `r23_cycle.png`;
  - `r23_preets_ramp.png`.
