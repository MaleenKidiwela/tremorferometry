# Q8 — Rapid tremor reversals (RTRs) in our LFE catalogues: results

Written 2026-10-08 ~12:40 pm PT. Rules: `PLAN.md` (§0–4 pre-registered 11:20 am PT; §5 = deviations V0 → D1 → D2 → D3).
Scripts `svi_science/scripts/q8_rtr_{lib,data,backzone,run,analyze,questions,extra,figures}.py`; CPU only, ≤ 8 cores,
< 2 GB. Times in tables are UTC; times quoted below are PT (PDT = UTC − 7). Final detector = **D2** (bursts of activity
behind the front that move against the main front). Every number below is from D2 unless marked otherwise.

## 1. In plain words

- **The step-1 reproduction of Royer 2015 PASSES PARTLY.** Royer, Thomas & Bostock picked their RTRs by eye. With an
  automatic detector we find about the same number in the same 7 episodes: 72 vs their 64. They come at the same time
  in the episode: a median of 12.5 vs 11.7 days after onset. Speeds, lengths and durations are alike: median 3.4 km/h vs
  4.2, 15 km vs 17, 2.5 h vs 4 h.
  - **But the one-to-one match fails the pre-registered 70 % bar:** 37 of their 64 (58 %) have one of ours within ±6 h.
    Chance gives 34 % (our list shifted by ±1–3 days). The matched pairs agree closely in time (median 0.9 h).
  - Three detector designs all stall at 58–67 % against 35–50 % chance. So hand picks and automatic picks agree on about
    half of the events, and we could not do better without tuning to their list.
- **Hour-scale reversals are real in Bostock's 2003–13 catalogue and they are not a detector artefact.**
  - Shuffling event times within each family-day leaves almost no such bursts: 1.6 vs 88 per 10 episodes.
  - The bursts that move have a strong preference for moving BACK against the main front: 88 backward vs 31 forward.
    When positions are scrambled among the families active that day, the split is 46 vs 46.
  - Bostock has 0–18 RTRs per episode (8.8 on average), mostly 4–20 days after onset.
- **Our own catalogues (2014–26) are far too sparse to see most RTRs**, but the few we do see behave the same way.
  - **We see few:** 8 RTRs in the 10 SVI episodes and 6 in the 14 WA-NW episodes. Our catalogues have ~30× fewer events
    per episode than Bostock's (median 640 vs 19,700).
  - **Bostock thinned to our density** detects 0.45 per episode. 82 % of those are real full-density RTRs, but only ~4 %
    of all RTRs are recovered. So most RTRs are invisible at our density. **The absence of an RTR in any 2014–26
    episode means nothing.**
  - **The count is what thinning predicts.** SVI 8 vs 4.5 expected (p 0.09); WA-NW 6 vs 6.3. Also their speed, length
    and duration (4.6 km/h, 11 km, 1.6 h) are exactly what thinned Bostock gives (4.8, 11.3, 1.4). So RTRs are no
    rarer in 2014–26 than in 2003–13, as far as we can tell.
  - The day-scale "back-running fronts" of S3 (SVI 2018-06, 2023-04) show no hour-scale RTR at our density, which is
    expected given the 4 % recall.
- **Tides: confirmed, and the strongest result.**
  - RTRs happen while the tide pushes the fault in the slip direction (positive up-dip shear, Royer's "UDSS").
  - With our own S5 tide model, **76 % of our 2003–13 RTRs lie entirely within positive up-dip shear**, against 31 %
    expected for random windows of the same length (p < 0.0001). Royer's own list gives 56 % vs 22 % with our tide
    model; they reported 72 % with theirs.
  - **2014–26 says the same: 7 of 8 SVI RTRs vs 40 % expected (p 0.009). This is new.**
  - Bursts that do not migrate are also tide-favoured, but less so (53 % vs 36 %), as are forward-moving bursts
    (65 % vs 36 %). So the tide favours all hour-scale bursts behind the front, and the backward ones most.
- **RTRs start near the main front and run back into the slipped zone.** In 2003–13 the median start is 5.7 km from the
  front and the median end is 14.4 km from it. The start is closer than typical activity in the slipped zone (8.5 km,
  p < 10⁻⁶), and 44 % start within 5 km of the front. This matches Houston et al. 2011. It is partly built in, because a
  streak that runs back from the front must start nearer to it than where it ends.
- **RTRs are stronger up-dip (confirms Bostock 2015 and Royer's "region A").**
  - Up-dip families put more of their LFEs into RTRs (Spearman −0.56 with down-dip position vs −0.26 for random time
    windows, p 0.01). RTR LFEs sit 3.3 km up-dip of the episode average.
  - **2014–26 gives the same, despite only 8 RTRs:** −0.26 vs −0.04 (p 0.01) and 3.3 km up-dip.
- **"RTRs hold the largest LFEs" is not specific to RTRs (new, a correction to the usual reading).**
  - In 2003–13, RTR LFEs are +0.09 magnitude units above their family's median, against +0.02 at random times
    (p 0.001).
  - But every other LFE firing in the same hours is just as big (+0.08; RTR minus same-hours difference +0.009, p 0.22).
    The big LFEs belong to the busy hours, not to the reversal.
  - Ours (S6 size, 8 RTRs): +0.08 vs +0.02, p 0.14, not significant. S6 already showed that our size measure rises with
    the tremor background.
- **Do RTRs repeat in the same places? Only a hint.**
  - 2003–13 RTR start positions per episode are more alike than chance (histogram correlation 0.54 vs 0.41, p 0.03).
    Their end positions are not (p 0.23), and the nearest-start test is not significant (p 0.24).
  - Starts spread over the whole SVI box (10–90 %: −28 to +35 km along strike). Our 2014–26 RTRs are too few to test.

**What is new vs Royer/Bostock/Houston:**
1. The first automatic, control-tested RTR detection in this catalogue.
2. The tidal preference of RTRs holds in a second, independent decade (2014–26).
3. The up-dip preference holds in 2014–26.
4. The "large-LFE" signature belongs to busy hours, not to the reversals themselves.
5. The tidal preference extends, more weakly, to every hour-scale burst behind the front.

## 2. Step 1 numbers (Royer 2015 reproduction; Bostock 2003–13, Royer's 7 episodes)

| detector | RTRs (7 ep.) | rule 43–96 | median days after onset | Royer matched ±6 h | chance (shift ±1–3 d / ±12–60 h) | verdict |
|---|---|---|---|---|---|---|
| V0 pre-registered (line scan, space/time independence) | 25 | ✘ | 13.9 ✔ | 13/64 = 20 % ✘ | 13 % | FAIL (S\* = 48 from a broken control) |
| D1 (rate surge along a line, ≥ 5 ordered families) | 0 at S\* = ∞; at ≥ 15 fam & S ≥ 50: 76 | – | – | 63 % | 39 % | FAIL |
| **D2 back-zone surge (final)** | **72** | ✔ | **12.5 ✔** (Royer 11.7) | **37/64 = 58 % ✘** | 34 % / 32 % | **count + timing pass, match fails** |

- Per episode (ours / Royer): 2003 5/11, 2004 10/7, 2005 18/10, 2008 11/6, 2010 8/10, 2011 7/8, 2012 13/12. Royer
  RTRs matched per episode: 36 %, 29 %, 80 %, 67 %, 40 %, 75 %, 75 %.
- Time zone check: Royer's times are UTC (recall 58 % at UTC vs 23–34 % if read as PT).
- Matched pairs: |Δt| median 0.9 h (quartiles 0.4–1.2 h). Their per-event speed/length/duration do not correlate with
  ours (Spearman 0.03 / 0.24 / 0.00): both sides measure noisily, and 12 of Royer's 64 lengths are exactly 16.7 km.
- Of Royer's RTRs, 95 % have SOME back-zone surge of ours (any direction) within ± 6 h (chance 70 %). The loss is in
  the direction and speed test.
- Both lists prefer the same UTC hours (Royer: 38/64 start 21:00–05:00 UTC). This inflates the day-shifted chance level.
- Tables: `tables/royer2015_rtr_list.csv` (transcribed; typos #11, #20, #23, #44 fixed and flagged),
  `royer_match_d2.csv`, `royer_recall_vs_S_d2.csv`, `results_step1_d2.json`. V0/D1 outputs are kept in
  `tables/v0_prereg/` and `*_bostock_only.*`.

## 3. Controls (all catalogues, D2, every RTR = S ≥ 2; per run of all episodes)

| catalogue | RTRs data | time-shuffled control | position-permuted control | forward surges (data) |
|---|---|---|---|---|
| Bostock 2003–13 (10 ep.) | 88 | 1.6 | 46 (forward 46) | 31 |
| Bostock thinned to 639 ev./ep. (10 copies) | 4.5 | 0.1 | 3.1 | 1.5 |
| ours SVI head (P ≥ 0.75, once; 10 ep.) | 8 | 0.8 | 2.8 | 2 |
| ours SVI P ≥ 0.5 once | 13 | 0 | 1.4 | 7 |
| ours SVI raw P ≥ 0.75 | 13 | 0 | 2.8 | 5 |
| ours WA-NW head (14 ep.) | 6 | 0.2 | 5.4 | 7 |
| ours WA-NW P ≥ 0.5 once | 15 | 0 | 4.4 | 11 |

- Time-shuffled = each event's time redrawn within its family-day. Position-permuted = within each UTC day, the
  families swap positions, keeping their times.
- The position-permuted control keeps the surges and scrambles the direction. So data minus control ≈ "direction
  beyond chance": Bostock ~42 of 88, SVI ~5 of 8.
- **WA-NW is not distinguishable from its position-permuted control** (6 vs 5.4). Its source positions are lumpy (S3:
  29 % of events from one source), and one NW "RTR" is 3 sources over 86 km in 1.3 h (2025-10-04). **Treat the WA-NW
  RTRs as unconfirmed.**
- SVI's own scrambled-control catalogue (stations shifted to different days) has 0–9 events per episode, too few to run.
- Pre-registered S\* (10 % false fraction) could not be met by any detector (deviation D3), so there is no score cut.
  Count-vs-score curves are in `tables/threshold_curves_d2.csv` and `figures/q8_summary.png`.

## 4. Questions (3a–e)

**a. Counts, timing, kinematics** (`tables/a_summary_d2.csv`, `a_per_episode_d2.csv`, `rtr_list_frozen_d2.csv`)

| | RTRs | per episode | days after onset (median) | speed km/h | span km | duration h | events in RTRs |
|---|---|---|---|---|---|---|---|
| Bostock 2003–13 | 88 in 9/10 ep. | 8.8 | 12.0 (3–30) | 3.4 | 15.4 | 2.5 | 17 % (Royer 6.8 %; ours counts all slipped-zone events in the burst) |
| Bostock thinned to our density | 4.5 / 10 ep. | 0.45 | – | 4.8 | 11.3 | 1.4 | – |
| ours SVI 2014–26 | 8 in 5/10 | 0.8 | 13.6 | 4.6 | 11.3 | 1.6 | 2.6 % |
| ours WA-NW 2015–26 | 6 in 4/14 | 0.43 | 14.9 | 9.8 | 23.7 | 1.4 | 1.1 % |

- Ours vs Bostock differ in span (p 0.015) and duration (p 0.02), but thinned Bostock shows the same shift. That is
  density, not physics.
- Timing does not differ (p 0.40–0.95). Speed: no difference (p 0.12–0.41).
- Regions: WA-NW spans are longer than SVI's (p 0.004). Given the lumpy NW positions, this is not interpreted.
- Our SVI RTRs (start, PT):
  - 2016-01-04 03:36 PST
  - 2020-11-10 00:32 PST and 11-14 16:51 PST
  - 2021-02-10 04:25 PST
  - 2022-04-05 21:13 PDT
  - 2024-08-19 17:33, 08-22 20:09 and 08-24 22:10 PDT

**b. Same places?** (pre-registered: nearest-start in another episode vs a random event of the same episode within
± 12 h of the start)
- Bostock: starts 0.43 vs 0.54 km, p 0.24; ends p 0.61.
- SVI 2014–26: p 0.57 / 0.40. WA-NW: p 0.88 / 0.10.
- Post-hoc (Bostock, start histograms in 10-km bins): correlation between episodes 0.54 vs control 0.41, p 0.028; ends
  0.46 vs 0.41, p 0.23.
- Verdict: **no clear repeat sites**; at most a hint that starts favour the same stretches.

**c. Start at the main front?** (D2 front = median of 4-km bins that switched on in the previous 24 h)

| set | n | start → front (median) | end → front | slipped-zone activity → front | percentile of start | p | starts ≤ 5 km |
|---|---|---|---|---|---|---|---|
| Bostock 2003–13 | 63 | 5.7 km | 14.4 km | 8.5 km | 0.33 | 9 × 10⁻⁷ | 44 % |
| ours SVI | 5 | 3.9 | 10.1 | 10.8 | 0.20 | 0.03 | 60 % |
| ours WA-NW | 4 | 12.4 | 13.4 | 15.9 | 0.35 | 0.31 | 0 % |

The plan's random-line reference (D1 geometry) does not exist for D2, so this test is partly built in (see §1).

**d. Tides (SVI; S5 hourly up-dip shear at the RTR's families, event-weighted)**

| set | n | whole RTR in positive up-dip shear: observed / expected | p | mean shear > 0: observed / expected | p |
|---|---|---|---|---|---|
| Royer's list (their times) | 64 | 56 % / 22 % | < 10⁻⁴ | 92 % / 55 % | < 10⁻⁴ |
| our RTRs 2003–13 | 88 | **76 % / 31 %** | < 10⁻⁴ | 92 % / 54 % | < 10⁻⁴ |
| our RTRs 2014–26 | 8 | **88 % / 40 %** | 0.009 | 100 % / 54 % | 0.009 |
| forward surges 2003–13 | 31 | 65 % / 36 % | 0.001 | 77 % / 55 % | 0.009 |
| non-migrating surges 2003–13 (post-hoc) | 135 | 53 % / 36 % | – | 72 % / 55 % | – |

- "Expected" is the chance that a window of the same length, started at a random time on the episode's active days, is
  all-positive.
- Royer's 72 % is a fraction that is "exclusively" positive. With our tide model, Royer's own list gives 56 %, and our
  list gives 76 %. The S5 model was validated on rates, not phases (S5 setup README).
- Caveat: a surge is a rate excess, and tides modulate the LFE rate (S5: ~35 % M2 modulation). Part of the preference is
  therefore expected for any surge. The RTR excess (76 vs 31) is larger than for non-migrating surges (53 vs 36).

**e. Sizes and depth**
- **Bostock magnitudes** (RTR LFE minus its family's median in the episode):
  - RTR LFEs: +0.086.
  - Random-time windows in the same stretch: +0.020 (97.5 % 0.041) → p 0.001.
  - Other LFEs firing in the same hours outside the RTR: +0.079. Paired RTR − same hours: +0.009, p 0.22.
- **Ours (S6 dM):**
  - P ≥ 0.75: +0.084 vs +0.021 shifted, p 0.14 (n 8).
  - P ≥ 0.5: +0.060 vs +0.030, p 0.27. Same hours: +0.028, paired −0.017, p 0.56.
- **Depth:**
  - Bostock: share of each family's LFEs that fall in RTRs vs down-dip position, Spearman −0.56; random-window control
    −0.26 (95 % −0.48…−0.02), p 0.010. RTR LFEs are 3.3 km up-dip of the average.
  - Ours: −0.26 vs −0.04 (95 % −0.24…+0.13), p 0.010; 3.3 km up-dip.
  - Per-family table: `tables/e_family_rtr_share_{bostock,svi}_d2.csv`.

## 5. Figures (`svi_science/results/q8_rtr/figures/`)
- `q8_examples.png`: space–time plots (along strike vs time, PT) for Bostock 2003, 2012 and 2005 with our RTRs (blue),
  forward bursts (yellow) and Royer's RTRs (pink bands); ours SVI 2024-08 and 2020-11; WA-NW 2015-12.
- `q8_bostock_royer_episodes.png`: all 7 Royer episodes, ours vs his.
- `q8_summary.png`:
  - count vs score, data vs the two controls (three catalogues);
  - Royer recall vs score threshold with chance;
  - speed and timing distributions.
- `q8_tides.png`: positive up-dip shear fractions, observed vs expected.

## 6. Deviations (full text in PLAN.md §5)
1. **V0 (as pre-registered) failed step 1.** The position-permuted control made more "streaks" than the data, and the
   detector drew lines between unrelated single-family bursts.
2. **D1 (rate surge + ≥ 5 ordered families + one streak per space–time box) also failed.** In Bostock's catalogue ~40
   families fire per hour, so some ordered subset always exists.
3. **D2 (bursts of activity behind the front, Theil–Sen direction, opposite to the episode's main direction) is final.**
   It was chosen after seeing its recall on Royer's list. This is calibration on the known answer, not a blind test.
4. **D3: no score threshold.** The position-permuted control keeps the bursts and tests only their direction.
5. **(c) uses D2's own front estimate and a slipped-zone reference.** The planned random-line reference needs line
   geometry, which D2 does not have.
6. **Post-hoc additions** (not pre-registered): thinned precision; start-histogram correlation; tides of non-migrating
   surges.
7. **Not done:**
   - WA tides: no WA tide setup exists.
   - WA-SW: only 2008-05 has > 100 events, so it was not run.
   - Bostock 2007-01 (2,847 events) and 2013-09 contribute 0 and 2 RTRs.

## 7. What would make this stronger
- A denser 2014–26 SVI catalogue (P ≥ 0.5 doubles the events, 15,312 vs 7,753, and gives 13 RTRs; its tides and depth
  were not run) or a tremor-based check of the 8 SVI RTRs.
- A blind test of D2 on an independent hand-picked RTR list (e.g. Houston 2011 WA tremor RTRs on PNSN tremor).

---

# Q8b — the same RTR finder on PNSN tremor, 2010–2026 (SVI, WA-NW, WA-SW)

Run 2026-10-08, ~1:00–3:05 pm PT. Rules: PLAN.md section Q8b, pre-registered ~12:50 pm PT before any PNSN number.
- Detector: the final Q8 finder (D2), unchanged. One event = one PNSN tremor window. Position = epicentre projected along
  strike.
- Script `scripts/q8b_rtr_pnsn.py` and figures `q8b_figures.py`. Outputs: `tables/q8b_*.csv`, `results_q8b.json`,
  `streaks_all_abs_pnsn_d2.parquet`.
- The PNSN catalogue starts 2010-01-06, so there is no 2009.

## Q8b.1 In plain words

- **Tremor shows hour-scale bursts in every ETS, and they are real.** Shuffling times within each day leaves almost none
  (1–3 per run vs 36–61). The finder labels 42 RTRs in 18 SVI episodes, 61 in 28 WA-NW and 36 in 28 WA-SW: about 1–2 per
  episode, in 15/18, 25/28 and 16/28 episodes.
- **In tremor, the backward preference is weak.** The LFE catalogue had 88 backward vs 31 forward (Q8).
  - SVI: 42 vs 30 (p 0.19). WA-NW: 61 vs 37 (p 0.02). WA-SW: 36 vs 39.
  - The position-permuted control makes about half as many migrating bursts as the data in SVI and NW (22 and 29 RTRs).
    In SW it makes as many (32 vs 36), so **SW "RTRs" are not separable from chance**.
  - The likely reason is tremor epicentre scatter (several km), which blurs short streaks and adds random slopes. So in
    tremor, "RTR" means "a migrating hour-scale burst", with only a mild backward bias outside NW.
- **Tremor and LFE RTRs agree only weakly where both exist (SVI, 2010–2013).**
  - 8 of the 30 LFE-based RTRs have a tremor RTR within ±6 h (27 %; chance 14 %). 8 of the 15 tremor RTRs in those
    windows match (53 %).
  - Against Royer's 2010–12 list: 27 % vs 14 % chance.
  - About twice chance, but most RTRs seen in one catalogue are missed by the other.
- **Tides: the strongest and most reproducible result holds in tremor too.**
  - SVI tremor RTRs lie wholly in positive up-dip shear 67 % of the time vs 38 % expected (p 0.0001).
  - By era: 2010–13 87 % vs 36 %. 2014–26 56 % vs 39% (p 0.055); the mean shear is positive in 89 % vs 55 %
    (p 0.0002).
  - Forward-moving tremor bursts are just as tidal (77 % vs 39 %); bursts that do not migrate less so (59 % vs 38 %).
  - **So, across two independent catalogues and 23 years, hour-scale migrating bursts in an ETS happen at positive tidal
    up-dip shear.** In tremor this is not specific to the backward ones.
- **Up-dip preference is not seen in tremor** (all boxes ≈ control). The LFE effect is a ~3 km up-dip shift, smaller than
  tremor location scatter, so this is a resolution limit, not a contradiction.
- **They start nearer the front than the typical slipped-zone activity**:
  - WA-NW p 0.0005, WA-SW p 0.001, SVI p 0.08;
  - medians 12–15 km from the front, blurred compared with 6 km in the LFEs.
- **No repeat sites** in any box (nearest-start p 0.10–0.34; histogram p 0.10–0.69).
- **SVI vs WA: the same bursts.**
  - Speed 4.1 / 4.9 / 4.1 km/h, length ~20–22 km and duration ~1.5 h, all equal (p ≥ 0.16).
  - SVI and NW have the same rate per episode (2.3 vs 2.2). SW has fewer (1.3, p 0.01); its box borrows the NW episode
    windows.
  - SVI bursts come later in the episode than NW ones (median 13.4 vs 10.4 days, p 0.02).
- **Are 2014–26 counts consistent with 2003–13? Yes, within what can be resolved; no consistent change.**
  - LFE catalogue: 8 SVI RTRs vs 4.5 expected at our density (Q8, p 0.09).
  - PNSN SVI, one catalogue across both eras: 3.8 per episode in 2010–13 (4 episodes) vs 1.9 in 2014–26 (14 episodes),
    p 0.04. Per 1,000 tremor windows: 0.96 vs 0.55.
  - But WA-NW goes the other way (1.3 → 2.4, p 0.12) and SW is flat (0.8 → 1.4, p 0.32). PNSN tremor detection also
    changed over 16 years.
  - So the SVI drop is a single marginal result, not a change in RTR behaviour.

## Q8b.2 Numbers

| box | episodes | tremor windows | RTRs (eps with ≥ 1) | per ep | forward | time-shuffled RTR | position-permuted RTR / fwd | median days after onset | speed km/h | length km | duration h |
|---|---|---|---|---|---|---|---|---|---|---|---|
| SVI | 18 (2010-08…2025-09) | 64,873 | 42 (15) | 2.3 | 30 | 2.8 | 22.2 / 17.0 | 13.4 | 4.1 | 21.4 | 1.5 |
| WA-NW | 28 (2010-08…2025-09) | 53,343 | 61 (25) | 2.2 | 37 | 1.2 | 28.8 / 28.0 | 10.4 | 4.9 | 20.3 | 1.6 |
| WA-SW (NW windows) | 28 | 32,000 | 36 (16) | 1.3 | 39 | 1.6 | 32.2 / 25.8 | 13.0 | 4.1 | 21.5 | 1.5 |
| (Q8 LFE, Bostock 2003–13) | 10 | 256,920 | 88 (9) | 8.8 | 31 | 1.6 | 46 / 46 | 12.0 | 3.4 | 15.4 | 2.5 |

- Per-episode counts and PT start times: `tables/q8b_per_episode.csv`. Full list: `tables/q8b_rtr_list.csv`.

**Cross-check, SVI 2010-08 / 2011-08 / 2012-08 / 2013-09 (Bostock windows), one-to-one ±6 h:**
- 30 LFE-based RTRs and 15 tremor RTRs; matched 8.
- Recall of LFE RTRs 27 % (chance 14 %, tremor list shifted ±1–3 d); precision of tremor RTRs 53 %.
- Any tremor burst (any class) within ±6 h of an LFE RTR: 40 % vs 27 % chance.
- Royer 2010–12 (30 RTRs): 8 matched, 27 % vs 14 % chance.

**Tides (SVI, S5 up-dip shear at the nearest family; whole RTR positive / expected):**

| set | n | whole RTR positive / expected | p | mean shear > 0 / expected | p |
|---|---|---|---|---|---|
| tremor RTRs, all | 42 | 67 % / 38 % | 0.0001 | 93 % / 55 % | < 0.0001 |
| tremor RTRs 2010–13 | 15 | 87 % / 36 % | 0.0001 | 100 % / 55 % | 0.0002 |
| tremor RTRs 2014–26 | 27 | 56 % / 39 % | 0.055 | 89 % / 55 % | 0.0002 |
| tremor forward bursts | 30 | 77 % / 39 % | < 0.0001 | 87 % / 55 % | 0.0002 |
| tremor non-migrating bursts | 22 | 59 % / 38 % | 0.04 | 77 % / 56 % | 0.03 |

- WA tides: not tested (no WA tide setup).

**Up-dip (4 × 4 km cells, ≥ 50 windows; RTR share vs down-dip position vs the random-window control):**
- SVI: ρ +0.03 vs −0.01 (p 0.66); RTR windows +0.3 km down-dip of the average.
- NW: +0.03 vs −0.02 (p 0.73); −0.2 km.
- SW: −0.11 vs −0.05 (p 0.22); −1.0 km.

**Start vs front:**

| box | n | start → front | end → front | slipped-zone reference | p (start closer) | starts ≤ 5 km |
|---|---|---|---|---|---|---|
| SVI | 21 | 12.3 km | 18.1 | 12.8 | 0.08 | 14 % |
| NW | 37 | 13.2 | 22.8 | 17.2 | 0.0005 | 19 % |
| SW | 23 | 15.5 | 21.3 | 18.0 | 0.001 | 9 % |

**Repeat sites:**
- Nearest start in another episode vs random-activity control:
  - SVI 1.42 vs 1.75 km (p 0.14);
  - NW 0.84 vs 1.01 (p 0.10);
  - SW 2.17 vs 2.55 (p 0.34).
- Start-histogram correlation: SVI p 0.32, NW p 0.10, SW p 0.69. None significant.

**Eras (PNSN, RTRs per episode):**
- SVI: 3.75 (2010–13, 4 eps) vs 1.93 (2014–26, 14 eps), p 0.04.
- NW: 1.33 (6) vs 2.41 (22), p 0.12.
- SW: 0.83 vs 1.41, p 0.32.

**Regions:**
- Speed, length and duration are the same in all pairs (p ≥ 0.16).
- Timing: SVI is later than NW (p 0.02).
- Rate: SW is lower than SVI and NW (p 0.01).

## Q8b.3 Figures
- `figures/q8b_pnsn_examples.png`: tremor along strike vs time (PT) with tremor RTRs (blue) and forward bursts (yellow).
  - SVI 2012-09 and 2011-08 also show the LFE-based RTRs (violet bands) and Royer's starts (pink dashes).
  - SVI 2023-04 and 2020-11, WA-NW 2018-05, WA-SW 2021-09.
- `figures/q8b_tides.png`: fraction wholly in positive up-dip shear, observed vs expected, for LFE RTRs (both eras),
  tremor RTRs (both eras), and tremor forward / non-migrating bursts.

## Q8b.4 Deviations and limits
- **No 2009:** the PNSN catalogue starts in 2010.
- **Episodes:** the PNSN-era big-ETS lists of B (`episodes_b.csv`, which extends the SVI list to 2010 with the same
  rule) were used. They include episodes not in the S3 list; the S3 subset is a column (`s3_big` in
  `q8b_per_episode.csv`). The SW box has no episode list of its own and borrows the NW windows.
- **Control "unit":** a 2-km along-strike bin, because tremor has no families. In the position-permuted control, a
  bin's windows take the position of another bin's first window.
- **No sizes:** PNSN tremor has no usable amplitude.
- **No WA tides:** no WA tide setup exists.

---

# Q8c — the same RTR finder on ALL our candidates (no P cut)

Run 2026-10-08, ~5:35–5:45 pm PT. Rules: PLAN.md section Q8c, pre-registered before computing.
- Scripts: `q8c_rtr_run.py`, `q8c_rtr_analyze.py`, `q8c_pw.py`, `q8c_figures.py`.
- Outputs: `tables/q8c_*.csv`, `results_q8c.json`, `results_q8c_pw.json`, `streaks_all_abs_q8c_d2.parquet`.
- Finder: D2, unchanged.

## Q8c.1 In plain words

- **Raw candidates do not help.** Taking every ≥ 8 MAD candidate gives SVI ~9,300 per big episode, about half of
  Bostock's density. But the scrambled-control candidates (stations shifted to different times) are ~60 % as numerous,
  and the RTR finder finds nearly as many "RTRs" in them as in the data:
  - SVI: 34 vs 29 (raw) and 30 vs 20 (counted once);
  - WA-NW: 5 vs 5 and 16 vs 20;
  - WA-SW: 3 vs 3.
  - **Raw and counted-once fail the pre-registered pass rule in every region.**
- **The scrambled-control candidates are not pure noise.** In 2010–13 their "RTRs" also line up with Bostock's RTRs
  above chance (27 % vs 14 %). The scrambled detector still responds when strong tremor hits one or two stations, so its
  output rises and falls with the real bursts.
- **P-weighting fixes it for SVI: this is the useful new catalogue result.** Each candidate is kept with probability P
  (5 draws).
  - SVI gives 25.4 RTRs per run (2.5 per episode) vs 5.2 in the scrambled-control candidates (ratio 4.9, per-episode
    Wilcoxon p 0.001) → **PASS**.
  - These RTRs strongly prefer running backward: 25 vs 11 forward (p 3 × 10⁻⁸).
  - Their rate, 2.5 per episode, is what thinned Bostock predicts at the same effective density (2.3 per episode at
    2,116 per episode). **So RTRs are as common in 2014–26 as in 2003–13.**
  - WA has no P for its scrambled-control rows, so the WA P-weighted RTRs (NW 23 per run, 28 forward; SW 0.8) cannot be
    tested and stay unconfirmed.
- **Tides separate real RTRs from noise RTRs.**
  - SVI P-weighted RTRs: 78 % lie wholly in positive up-dip shear vs 36 % expected. Forward bursts are equally tidal
    (89 % vs 42 %), non-migrating bursts less so (61 % vs 39 %), as in tremor (Q8b).
  - The scrambled-control RTRs are not tidal: 35 % vs 42 % (P-weighted), and in 2010–13 12 % and 0 % vs 33–35 %.
  - Even in the failing raw and counted-once variants, the data RTRs are tidal (56–75 % vs 33–35 %, p ≤ 0.008), while
    the counted-once and 2010–13 scrambled-control ones are not.
- **Calibration 2010–13 (our forward-detector candidates in Bostock's ETS windows):**
  - 23–24 RTRs, matching 43–47 % of the Bostock-based RTRs (chance 21–23 %; precision 54–61 %).
  - Royer: 30–33 % matched (chance 19–23 %).
  - Backward 23–24 vs forward 9–10 (p 0.02).
  - The scrambled-control series gives 10–17 RTRs, so the pass rule fails here too. With only 4 episodes, the Wilcoxon
    part of the rule cannot go below p 0.0625; the counted-once ratio is 2.4.
- **Other answers (counted-once SVI, labelled "not separable from noise" by the pass rule):**
  - RTR LFEs are more up-dip (Spearman −0.48 vs −0.16 in the control, p 0.005; 3.9 km up-dip). Same in 2010–13
    (p 0.025).
  - Starts are nearer the front than ends: 16 vs 32 km (p 0.003).
  - No repeat sites (p 0.18–0.27).
  - Speed 4–5 km/h, length 17–28 km, duration ~2 h, median 15 days after onset.
- **SVI vs WA cannot be compared on candidates.**
  - WA-NW candidates are ~93 % matched by scrambled-control candidates (52k vs 48k per episode).
  - The NW "RTRs" do not exceed the scrambled control, and NW shows 5–16 backward vs 121–161 forward bursts, a noise
    pattern.
  - WA-SW has only ~375 candidates per episode.

## Q8c.2 Numbers

**Pass rule (per region and variant): data RTRs ≥ 2× scrambled-control RTRs AND per-episode Wilcoxon p < 0.05**

| set | variant | candidates / ep (scrambled) | RTRs (per ep) | scrambled RTRs | ratio | Wilcoxon p | forward | position-permuted RTR/fwd | PASS |
|---|---|---|---|---|---|---|---|---|---|
| SVI 2014–26 | raw | 9,316 (5,479) | 34 (3.4) | 29 | 1.2 | 0.32 | 38 | 30 / 35 | ✘ |
| SVI 2014–26 | counted once | 6,151 (3,802) | 30 (3.0) | 20 | 1.5 | 0.12 | 33 | 23 / 29 | ✘ |
| SVI 2014–26 | **P-weighted (5 draws)** | 2,116 (P-sum) | **25.4 (2.5)** | **5.2** | **4.9** | **0.001** | 10.6 | – | **✔** |
| WA-NW | raw | 51,917 (48,045) | 5 (0.36) | 5 | 1.0 | 1.0 | 161 | 47 / 67 | ✘ |
| WA-NW | counted once | 36,347 (33,734) | 16 (1.1) | 20 | 0.8 | 0.92 | 121 | 46 / 68 | ✘ |
| WA-NW | P-weighted | – | 23 (1.6) | no P for scrambled rows | – | – | 27.8 | – | untestable |
| WA-SW | raw / once | 375 (337) | 3 / 2 | 3 / 2 | 1.0 | 0.5 | 4 / 4 | – | ✘ |
| 2010–13 fwd-det. | raw | 13,671 (6,264) | 23 (5.8) | 17 | 1.4 | 0.13 | 9 | 12 / 12 | ✘ |
| 2010–13 fwd-det. | counted once | 9,950 (4,919) | 24 (6.0) | 10 | 2.4 | 0.063 (floor with 4 eps) | 10 | 13 / 16 | ✘ |

- Time-shuffled control: 0–3 RTRs per run everywhere.

**Thinned Bostock at the new densities (10 copies; RTRs per episode, precision / recall vs full-density RTRs):**
- 9,316 per episode: 6.5 (89 % / 66 %).
- 6,150: 5.8 (86 % / 56 %).
- 2,116: 2.3 (82 % / 21 %).
- 375: 0.13 (77 % / 1 %).
- Full density: 8.8.
- Raw and counted-once SVI (3.4 and 3.0 per episode) fall below the 6.5 and 5.8 expected. That comparison is not fair,
  because ~60 % of the raw candidate count is matched by scrambled-control candidates; only the P-weighted comparison
  (2.5 vs 2.3) is like for like.

**Calibration 2010–13 (Bostock windows 2010-08, 2011-08, 2012-08, 2013-09; ±6 h one-to-one; chance = list shifted
±1–3 d):**

| candidates | RTRs | Bostock-RTRs matched (of 30) | chance | precision | Royer matched (of 30) | chance |
|---|---|---|---|---|---|---|
| raw data | 23 | 14 (47 %) | 21 % | 61 % | 9 (30 %) | 19 % |
| counted-once data | 24 | 13 (43 %) | 23 % | 54 % | 10 (33 %) | 23 % |
| raw scrambled | 17 | 8 (27 %) | 14 % | 47 % | 10 (33 %) | 18 % |
| counted-once scrambled | 10 | 6 (20 %) | 8 % | 60 % | 6 (20 %) | 9 % |

**Kinematics (medians):**

| set | speed km/h | length km | duration h | days after onset |
|---|---|---|---|---|
| SVI P-weighted | 3.9 | 16.8 | 1.9 | 15.2 |
| SVI raw | 4.3 | 25.5 | 1.9 | 14.4 |
| SVI counted once | 4.6 | 28.0 | 2.0 | 16.0 |
| 2010–13 raw | 3.1 | 30.2 | 2.5 | 16.1 |
| 2010–13 counted once | 3.5 | 22.4 | 2.5 | 15.8 |
| (Bostock LFE, Q8) | 3.4 | 15.4 | 2.5 | 12.0 |

**Tides (SVI families; whole RTR in positive up-dip shear, observed vs expected):**

| set | data RTRs | forward | non-migrating | scrambled-control RTRs |
|---|---|---|---|---|
| SVI P-weighted (5 draws pooled; p with n/5) | 78 % vs 36 % (p 8 × 10⁻⁶) | 89 % vs 42 % (p 0.008) | 61 % vs 39 % (p 0.003) | 35 % vs 42 % (p 0.69) |
| SVI raw | 56 % vs 35 % (p 0.008) | 50 % vs 36 % (p 0.04) | 40 % vs 32 % (p 0.03) | 55 % vs 39 % (p 0.055) |
| SVI counted once | 60 % vs 35 % (p 0.004) | 42 % vs 35 % (p 0.25) | 42 % vs 31 % (p 0.008) | 35 % vs 41 % (p 0.76) |
| 2010–13 raw | 61 % vs 33 % (p 0.004) | 67 % vs 34 % (p 0.05) | 41 % vs 31 % (p 0.04) | 12 % vs 35 % (p 0.99) |
| 2010–13 counted once | 75 % vs 33 % (p 5 × 10⁻⁵) | 50 % vs 32 % (p 0.18) | 42 % vs 32 % (p 0.03) | 0 % vs 33 % (p 1.0) |

- WA tides: not tested (no WA tide setup).

**Up-dip (family RTR share vs down-dip; random-window control):**
- SVI counted once: −0.48 vs −0.16 (p 0.005); RTR LFEs 3.9 km up-dip.
- SVI raw: −0.17 vs −0.10 (p 0.25).
- 2010–13 counted once: −0.36 vs −0.12 (p 0.025). 2010–13 raw: −0.34 vs −0.11 (p 0.06).
- NW raw: p 0.005, but NW does not pass the noise rule.
- P-weighted up-dip not computed.

**Start vs front (start → front / end → front):**
- SVI raw: 18 / 26 km (p 0.01).
- SVI counted once: 16 / 32 km (p 0.003).
- 2010–13: 11–12 / 19–22 km (p 0.03 / 0.001).

**Repeat sites:**
- SVI raw: nearest-start p 0.043, histogram p 0.24.
- SVI counted once: 0.27 / 0.18.
- 2010–13: 0.22–0.71 / 0.48–0.57.
- NW counted once: 0.95.
- Verdict: no repeat sites.

**SVI vs WA:** WA does not pass the noise rule, so it is not compared. For the record, NW lengths are longer (46–49 vs
26–28 km, p 0.02).

## Q8c.3 Figures
- `figures/q8c_examples.png`:
  - SVI 2018-06, 2023-04, 2024-08 and 2015-12: all candidates, with P-weighted draw-0 RTRs (blue), forward bursts
    (yellow) and scrambled-control RTRs (red dashed).
  - 2012 and 2010 forward-detector candidates (counted once), with Bostock-LFE RTRs as violet bands.
- `figures/q8c_tides.png`: tidal fraction, data vs scrambled-control RTRs, observed vs expected, for each variant.

## Q8c.4 Deviations and limits
- **Wilcoxon floor.** With 4 calibration episodes, the Wilcoxon part of the pass rule cannot reach p < 0.05 (floor
  0.0625).
- **P-weighted pass-rule test.** It was computed in a separate add-on script (`q8c_pw.py`) after the main summary,
  which had only computed the counts. Same rule, no change.
- **No 2010–13 P-weighted variant.** The backward production has no P yet, so variant (c) is not available for
  2010–13.
- **No WA scrambled-control check of P-weighted RTRs.** The WA scrambled-control rows have no P.
- **Tide p-values for pooled P-weighted draws.** The 5 draws share events, so p is also given with the effective n
  divided by 5.

---

# Q8d — WA scrambled-control candidates get their P; the WA P-weighted RTR test

Run 2026-10-08 ~5:55–6:00 pm PT. Rules: PLAN.md section Q8d.
- Scripts: `q8d_wa_scrambled_P.py`, `q8d_rtr_run.py`, `q8d_analyze.py`.
- Outputs: `tables/wa_scrambled_P_2015_2026.parquet` (new file, 6,661,852 scrambled-control rows of NW/SW families),
  `wa_scrambled_P_join_check.csv`, `q8d_per_episode_{NW,SW}.csv`, `results_q8d.json`.

## Q8d.1 In plain words

- **No new scoring was needed, and nothing was improvised.** The frozen WA confidence model had already been applied in
  production to the scrambled-control candidates, with the same WAVE + PhaseNet features as the data rows
  (`wa_prod/conf/P/prod`). I only joined those values to the scrambled rows.
  - All 6.66 M rows have a P (0 missing, times identical).
  - The same join reproduces the catalogue's data-row P to 3 × 10⁻⁸.
  - No WA file was edited.
- **WA-NW passes the noise rule narrowly.** P-weighted data give 23.8 RTRs per run (1.7 per episode) vs 11.0 in the
  scrambled-control candidates (ratio 2.2, per-episode Wilcoxon p 0.018). Time-shuffling leaves 0.6.
- **But WA-NW bursts have no backward preference:** 23.8 backward vs 26.2 forward. Scrambling positions gives the same
  split (19.6 / 23.8). So in Washington-NW we see hour-scale migrating bursts in ETS, but not a reversal signature.
  SVI had 25 vs 11, and PNSN NW tremor had a weak 61 vs 37.
- **WA-NW has fewer such bursts than SVI**: 1.7 vs 2.5 per episode (p 0.002). Per 1,000 units of P-sum (the expected
  number of real LFEs) it is 0.44 vs 0.94 (p 4 × 10⁻⁹), so the deficit is not from detection density.
- **WA-NW bursts start nearer the front than they end** (15 vs 28 km, p 10⁻⁶).
- **Up-dip in WA-NW is mixed:** the family-share ranking leans up-dip relative to the control (p 0.02), but the burst LFEs
  sit 3.7 km down-dip of the average. Not interpreted.
- **No repeat sites** (p 0.69). Speed ~6.7 km/h, length ~37 km, duration ~1.7 h, median 14 days after onset.
- **WA-SW cannot be tested.** Its confident candidates are too few (P-sum 66 per episode): 1.0 RTR per run vs 0 in the
  scrambled control, p 0.09 → fails.
- **Tides not tested for WA:** no WA tide setup exists.

## Q8d.2 Numbers (P-weighted, 5 draws; per run = all 14 episodes)

| | WA-NW | WA-SW | SVI (Q8c, reference) |
|---|---|---|---|
| P-sum per episode (scrambled) | 3,844 (1,848) | 66 (43) | 2,116 (–) |
| RTRs per run (per episode) | 23.8 (1.7) | 1.0 (0.07) | 25.4 (2.5) |
| scrambled-control RTRs per run | 11.0 | 0.0 | 5.2 |
| ratio / Wilcoxon p / PASS | 2.2 / 0.018 / ✔ | – / 0.09 / ✘ | 4.9 / 0.001 / ✔ |
| forward per run (scrambled) | 26.2 (35.0) | 1.4 (0) | 10.6 |
| backward vs forward binomial p | 0.49 | 0.77 | 3 × 10⁻⁸ |
| time-shuffled RTRs | 0.6 | 0 | – |
| position-permuted RTR / forward | 19.6 / 23.8 | 1.4 / 0.8 | – |
| RTRs per 1,000 P-sum | 0.44 | 1.08 | 0.94 |
| vs SVI: per episode p / per P-sum p | 0.002 / 4 × 10⁻⁹ | ≪ 0.001 / 0.63 | – |
| speed / length / duration | 6.7 km/h / 37 km / 1.7 h | 11.8 / 16 / 1.1 | 3.9 / 17 / 1.9 |
| start → front / end → front | 15 / 28 km (p 1 × 10⁻⁶) | – | 16–18 / 26–32 km |
| up-dip (Spearman vs control; RTR LFE − all) | −0.03 vs +0.16, p 0.02; +3.7 km down-dip | −0.38 vs −0.18, p 0.045 (1 RTR per run) | −0.48 vs −0.16 (counted once) |
| repeat sites (draw 0) | p 0.69 | – | none |

## Q8d.3 Limits
- **Pooled draws.** The 5 P-weighted draws share events, so pooled up-dip/front p-values are optimistic; repeat sites
  use draw 0 only.
- **WA-SW** uses the NW episode windows (no SW list).

---

# Q8e — why WA-NW bursts show no backward preference (diagnosis; not fixed)

Run 2026-10-08, ~10:05–10:25 pm PT. Rules: PLAN.md section Q8e, pre-registered. The two counted-once variants (Pd, Pe)
were added afterwards and are labelled post-hoc.
- Scripts: `q8e_nw_tides.py`, `q8e_rtr_run.py`, `q8e_analyze.py`.
- Outputs:
  - `tables/q8e_variants_ranked.csv`, `results_q8e.json`;
  - `q8e_nw_family_locations.csv`;
  - new tide files `tidal_stress_nw_sources_hourly.parquet`, `nw_fault_geometry.csv`, `nw_ocean_constituents.csv`.
- Figure: `figures/q8e_wa_nw_variants.png`.
- No WA file was edited.

## Q8e.1 In plain words

- **Not fixed.** None of the 10 versions gives WA-NW a backward preference that passes the rule (backward share above the
  position-permuted control with p < 0.01, AND data clearly above the scrambled control).
- **The best lead is location quality.** Keeping only families at well-located (quality A) sources raises the backward
  share to 0.62 (18 backward, 11 forward per run) against 0.47 in the position-permuted control. But p is only 0.07,
  and these data RTRs are just 1.6× the scrambled-control count, below the 2× bar.
- **Ruled out as the explanation, each by its own test:**
  - **Duplicate families** (132 families = 100 sources). One family per source: share 0.38 vs control 0.43.
  - **Front definition.** Taking the direction and the slipped zone from PNSN tremor: 0.55 vs 0.52. Single-front
    episodes only: 0.51 vs 0.51; with both: 0.52 vs 0.49.
  - **Axis.** Our strike axis is already the slab-contour strike: the 41.8-km Slab2 contour across the box runs at
    345.2°, 0.5° from the 344.7° in use. Rerun: 0.49 vs 0.45.
  - **Cross-family co-detections** (post-hoc). Counting each event once: 0.38 vs 0.44; with quality-A sources:
    0.50 vs 0.47.
- **What WA-NW bursts do show:**
  - They are more compact along strike than chance (data ~30–35 km vs ~44–49 km with positions permuted).
  - Time-shuffling destroys them (0–1 per run).
  - They are about as frequent as the scrambled-control count allows.
  - **But their direction is random.**
- **The pass in Q8d was fragile.** The same baseline with new random draws gives data / scrambled = 1.8 (p 0.035),
  below the 2× bar; Q8d's draws had given 2.2 (p 0.018).
- **Tides for WA-NW (new set-up):**
  - The S5 machinery at the 100 NW sources gives up-dip shear with std 0.24–2.6 kPa per source (median 0.70).
  - WA-NW bursts lie wholly in positive up-dip shear more often than expected: baseline RTRs 66 % vs 44 % (p 0.03);
    non-migrating 58 % vs 46 %.
  - **The scrambled-control bursts are just as tidal** (62–94 % vs ~48 %), unlike SVI. In WA, tides cannot tell
    real bursts from scrambled-control ones.
  - That fits the main suspicion: the WA candidates, including the scrambled-control ones, are driven by the real
    tremor (scrambled control ≈ 93 % of the data candidate count).
- **Most likely reason** (a judgement, not a test result): the WA-NW detections do not place an event at the right
  source position often enough to trace an hour-scale streak.
  - Location quality is the only factor that moves the share upward.
  - 38 of 132 families have erh > 10 km, and 20 of 100 sources sit > 30 km off the down-dip centre of the box.
  - The tremor shows a weak backward bias (Q8b, 61 vs 37), so WA-NW RTRs probably exist.
  - A real fix needs better WA-NW locations (or detections tied to single sources), not a different RTR rule.

## Q8e.2 Ranked variants (P-weighted, 5 draws; per run = all episodes of the variant)

| rank | variant | families / sources / eps | RTR / forward per run | backward share | position-permuted share | p (share > control) | scrambled RTR per run | data / scrambled, Wilcoxon p | span data / permuted km | fixed? |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Pb quality-A sources | 79 / 55 / 14 | 18.0 / 11.2 | 0.62 | 0.47 | 0.073 | 11.2 | 1.6, 0.029 | 32 / 43 | ✘ |
| 2 | A1 slab-contour axis 345.2° | 132 / 100 / 14 | 23.6 / 24.2 | 0.49 | 0.45 | 0.28 | 10.8 | 2.2, 0.005 | 35 / 49 | ✘ |
| 3 | F1 tremor front | 132 / 100 / 14 | 24.6 / 20.0 | 0.55 | 0.52 | 0.38 | 18.6 | 1.3, 0.004 | 35 / 46 | ✘ |
| 4 | F3 tremor front, single-front eps | 132 / 100 / 11 | 20.2 / 18.8 | 0.52 | 0.49 | 0.43 | 16.0 | 1.3, 0.25 | 33 / 46 | ✘ |
| 5 | Pe (post-hoc) counted once + quality A | – / 55 / 14 | 13.8 / 13.8 | 0.50 | 0.47 | 0.47 | 4.8 | 2.9, 0.002 | 29 / 45 | ✘ |
| 6 | Pc one family per quality-A source | 55 / 55 / 14 | 12.8 / 15.6 | 0.45 | 0.46 | 0.54 | 6.2 | 2.1, 0.013 | 30 / 44 | ✘ |
| 7 | F2 single-front episodes | 132 / 100 / 11 | 22.4 / 21.6 | 0.51 | 0.51 | 0.62 | 7.6 | 2.9, 0.006 | 35 / 50 | ✘ |
| 8 | Pa one family per source | 100 / 100 / 14 | 16.6 / 27.4 | 0.38 | 0.43 | 0.76 | 6.2 | 2.7, 0.045 | 35 / 49 | ✘ |
| 9 | Pd (post-hoc) counted once | – / 100 / 14 | 16.6 / 27.6 | 0.38 | 0.44 | 0.79 | 5.8 | 2.9, 0.043 | 35 / 48 | ✘ |
| 10 | V0 baseline (fresh draws) | 132 / 100 / 14 | 23.4 / 27.0 | 0.46 | 0.52 | 0.85 | 12.8 | 1.8, 0.035 | 35 / 49 | ✘ |

- Burst speeds are ~5.5–6.8 km/h in the data and 5.4–6.6 km/h in the position-permuted control: no difference.
- Time-shuffled control: 0–1 RTR per run in every variant.

**Tides (WA-NW; whole burst in positive up-dip shear, observed vs expected; p with effective n = pooled / 5):**

| variant | RTR | forward | non-migrating | scrambled-control RTR |
|---|---|---|---|---|
| V0 | 66 % vs 44 % (p 0.03) | 53 % vs 47 % (0.38) | 58 % vs 46 % (0.03) | 62 % vs 48 % (0.25) |
| Pb | 69 % vs 46 % (0.06) | 64 % vs 47 % (0.22) | 54 % vs 47 % (0.10) | 71 % vs 49 % (0.11) |
| F1 | 66 % vs 43 % (0.03) | 63 % vs 49 % (0.11) | 59 % vs 46 % (0.01) | 76 % vs 48 % (0.02) |
| A1 | 68 % vs 44 % (0.02) | 58 % vs 47 % (0.18) | 56 % vs 46 % (0.02) | 80 % vs 48 % (0.02) |
| Pd (post-hoc) | 70 % vs 45 % (0.03) | 60 % vs 47 % (0.10) | 52 % vs 46 % (0.24) | 55 % vs 48 % (0.61) |

All variants: `results_q8e.json`.

**Family locations (diagnostic 1; `tables/q8e_nw_family_locations.csv`):**
- 132 setup NW families map onto 100 array-confirmed sources (1–4 families each).
- Methods: 131 Hyp2000 (GSC 1-D, Vp/Vs 1.73, free depth); 1 approximate (clear-arrival station centroid).
- Family location quality A/B/C = 64/44/24; erh median 5.1 km, 38 families > 10 km.
- Source quality A/B/C = 55/30/15. Sources span −147 … +68 km down-dip; 20 sit > 30 km off the box centre.
- Flags in families_v1_3: 32 families have `nw_copy_of` set; `nw_link_type` same_source_array 54, inseparable 4,
  none 74.

## Q8e.3 Notes and limits
- **WA-NW tides use the S5 SVI elastic constants and gauge field unchanged.** The inner gauge polygon covers Puget Sound
  and Hood Canal. Depth = Slab2 depth at the source. Phases were not checked against an independent WA tidal-stress
  result.
- **Pooled draws:** p-values use n divided by 5. The baseline's ratio moved from 2.2 to 1.8 between two sets of 5 draws,
  so P-weighted ratios near 2 are not stable.
