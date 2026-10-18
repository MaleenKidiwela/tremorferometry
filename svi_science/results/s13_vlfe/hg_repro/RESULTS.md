# hg_repro RESULTS — Hutchison & Ghosh (2019) reproduced, then our controls

Written 2026-10-15 ~4:40 am PT. Rules: `PLAN.md` (written ~2:15 am PT before any number, Amendments 1–2 before the
matched-filter numbers). Log: `NOTES.md`. Periods are the paper's UTC dates: 2011 ETS 23 Jul–6 Sep 2011, between ETS
6 Sep 2011–30 Aug 2012, 2012 ETS 30 Aug–11 Oct 2012. "Detections" = H&G-style matched-filter detections; all five 2011
templates; all nine Fig.-1 stations.

## Plain answer
- **The paper's literal threshold finds almost nothing.** It is the mean third-highest noise-template peak over the full
  2.5-year noise run. With it we get 7 detections in 18 months, against ~4,400 in the paper.
- **A threshold set to the paper's own false rate reproduces the shape of their catalogue, at about half the count.**
  - The threshold is set where noise templates fire 0.79 times per template per day, the paper's Fig.-6 background
    (median + 6×MAD gives nearly the same).
  - We get 4.4 detections/day, against 8.5 in the paper.
  - The rate rises clearly in both ETS (7.3 and 7.8/day) and drops between ETS (3.8/day, paper 8.2).
- **Their stack-and-invert check does not tell real repeats from noise.**
  - Windows picked because they correlate with any 90-s window stack into that window.
  - Stacks of noise-template detections match their own noise window at CC 0.60–0.79 (4 of 5 noise templates). Real
    stacks match their templates at 0.22–0.82.
  - With our whole-space Green's functions, neither real nor noise-template stacks consistently pass a VR or location test.
- **In the ETS, many detections look real; between ETS, the answer depends on which false-alarm reference is used.**
  - In the ETS, 70–74 % of detections survive our scrambled controls, and detections carry our LFEs at 2.5–3.8× the
    random rate.
  - Between ETS, 45 % survive our scrambled controls, but about 0 % survive the paper's own noise-template background.
  - Between-ETS detections carry no more LFEs than random times.

## 1. Known-answer check (K1–K4; `tables/known_answer.csv`, `period_counts.csv`)
Raw figures first: `figures/fig6_repro_thr_rate079.png` (ours on top, the paper's digitised Fig. 6 below, same tremor
histogram) and `figures/fig7_repro_thr_rate079.png` (cumulative per template, ours red, paper dashed).
Literal versions: `fig6_repro_thrA_3rd_incl_self.png`, `fig7_repro_thrA_3rd_incl_self.png`.

| quantity (all 5 templates) | paper (digitised) | literal 3rd-highest | median + 6×MAD | noise rate = 0.79/day |
|---|---|---|---|---|
| threshold, mean CC per channel | — | 0.43 | 0.30–0.32 | 0.31–0.33 |
| detections/day, 1 Jul 2011–1 Dec 2012 (K1) | 8.5 | 0.013 | 3.7 | **4.4 (pass)** |
| 2011 ETS /day | 10.7 | 0.02 | 5.9 | 7.3 |
| between ETS /day | 8.2 | 0.003 | 3.1 | 3.8 |
| 2012 ETS /day | 10.3 | 0.12 | 7.6 | 7.8 |
| ETS / between-ETS rate (K2, 1.0–1.8) | 1.27 | 25 | 2.2 (fail) | 2.0 (fail) |
| noise templates, per template per day (K3) | 0.6 text / 0.79 Fig. 6 | 0.02 (fail) | **0.87 (pass)** | **0.80 (pass)** |
| count 1 Jul–1 Oct 2011 | 1,394 (text run) / ~776 (Fig. 6) | 2 | 412 | 501 |
| largest 2011-ETS rise (K4) | T1 | T3 | **T1 (pass)** | T2 (T1 second, 1.88) |

Per template, rate relative to its own 18-month mean (2011 ETS / between / 2012 ETS), noise rate = 0.79/day threshold:

| | T1 | T2 | T3 | T4 | T5 |
|---|---|---|---|---|---|
| paper | 2.14 / 0.88 / 1.36 | 0.99 / 1.01 / 1.09 | 0.63 / 1.08 / 0.57 | 1.06 / 0.97 / 0.82 | 1.46 / 0.93 / 0.91 |
| ours | 1.88 / 0.83 / 1.30 | 2.31 / 0.72 / 2.51 | 1.10 / 0.99 / 1.35 | 1.85 / 0.75 / 3.07 | 1.53 / 0.82 / 2.05 |

- **T1 reproduces closely** (2011 ETS 1.9× vs 2.1×; 2012 ETS 1.3× vs 1.4×), and its cumulative curve lies on the paper's.
- **T2, T4 and T5 are more ETS-locked in ours than in the paper.** We find 2–3× rises in the 2012 ETS (T2, T4) where the
  paper's Fig. 7 shows ~1×.
- **Our between-ETS rate is half of theirs**, at the same noise-template false rate. So their extra between-ETS
  detections come from something in their setup we cannot see:
  - the station sets (not published);
  - the declustering (not stated);
  - the separate runs (text 1,394 vs Fig. 6 ~776 for the same window).
- **Threshold reading.** The paper cannot have used the literal third-highest reading: with it, a noise template fires
  ~0.003 times per day by construction. The rate-calibrated and 6×MAD readings both match their stated background, so
  below we use **noise rate = 0.79/day** (the K1 winner, `tables/best_threshold.json`) and report 6×MAD alongside.

## 2. False-alarm rates: theirs vs ours (`tables/noise_rate.csv`, `controls_counts_*.csv`)
Raw figures: `figures/controls_daily_thr_rate079.png` (daily real vs both scrambled controls vs noise templates, with
tremor) and `figures/controls_separation_thr_rate079.png` (detections/day vs threshold, ETS and between ETS, all 16
control draws drawn individually).

| per template per day, noise rate = 0.79/day threshold | 2011 ETS | between ETS | 2012 ETS |
|---|---|---|---|
| real detections | 1.46 | 0.75 | 1.56 |
| H&G-style noise templates (their reference) | 1.01 | 0.75 | 0.95 |
| station-day-shift scrambled control (8 draws) | 0.36 | 0.36 | 0.34 |
| polarity-scrambled control (8 draws) | 0.29 | 0.31 | 0.32 |

- **Our scrambled controls fire ~2.5× less than the paper's noise templates.** The noise templates are real 90-s network
  windows, so they keep the true station-to-station coherence of the 20–50 s wavefield. That makes them match other
  coherent-wave times more often. Our controls break that coherence. As a false-alarm reference here, the noise templates
  are the stricter and more realistic one.
- **The noise templates also fire more in the ETS** (1.0 vs 0.75 per day). They were drawn from 2011 ETS data, so they
  may themselves hold some ETS signal or ETS-season noise. The paper's statement of "no clear trend" in its noise
  detections does not hold for ours.
- **Between ETS, the real rate (0.75) equals the noise-template rate (0.75).** In the ETS, the real rate is 1.5× the
  noise templates.

## 3. Stack-and-invert: is their validation discriminating? (`tables/stackinv_*.csv`, `stack_selfsimilarity_*.csv`)
Raw figures:
- `figures/fig2_repro_thr_rate079.png`: the paper's Fig.-2 panels redone. They look like theirs: a clear stack over
  noisy grayscale rows.
- `figures/fig2_selection_bias.png`: the same display for T2 detections, for the detections of a noise template under the
  T2 configuration, and for random times.
- `figures/stackinv_map_thr_rate079.png`: map of all solutions.

**Waveform self-similarity** (CC over 0–90 s, all channels), noise rate = 0.79/day threshold:

| stack | CC with the window it was selected by |
|---|---|
| real detections, T1–T5 (all) | 0.34, 0.54, 0.25, 0.62, 0.22 |
| real, ETS only | 0.82, 0.53, 0.44, 0.66, 0.48 |
| real, between ETS only | 0.75, 0.42, 0.25, 0.53, 0.70 |
| noise-template detections (5 noise templates, T2 config) | 0.70, 0.66, 0.76, 0.60, 0.15 |
| random times (vs the template) | −0.22 to +0.14 |

- A stack of noise-template detections resembles its noise window as strongly as a real stack resembles its template.
  So "the stack looks like the template" (their Fig. 2) holds for noise too; it cannot confirm repeating events.
- Why T1 "all" (0.34) is below its ETS (0.82) and between-ETS (0.75) parts: the linear stack includes the 2011-07-21
  pre-ETS days. That is when UW.LRIV has a 3.6-cm glitch, which dominates a linear mean.

**Grid-search CMT** (whole-space GF, τ = 10 s chosen on the templates; coarse 0.1° × 5 km, fine 0.025° × 1 km):
- **Template events themselves.** VR is 0.13–0.47. Only T2 (13 km) and T4 (11 km) relocate near Table 1. T1, T3 and T5
  land 49–104 km away, mostly on grid edges.
- **Comparison with the paper.** It reports VR > 40 % for every stacked inversion; we reach that for none of the
  template events except T4 (0.47). So our inversion tool cannot reproduce their stack-inversion match, and the stack
  VR/location numbers below are not a fair test of their claim.
- **Stacks under the same tool.**
  - Real "all" stacks: VR 0.02–0.39.
  - Noise-template stacks: VR 0.09–0.26.
  - Random-time stacks: VR 0.01–0.72. The highest single VR (0.72) came from a random stack.
  - Only T1 stacks (ETS and between ETS) invert close to the T1 template solution (14 km, Kagan 9–11°). One noise-template
    stack also lands within Kagan 4° of T1.
- **Verdict.** The waveform part of their validation is not discriminating (selection bias). The inversion part we
  could not reproduce with whole-space GFs; a 1-D layered GF would be needed to judge it.

## 4. Survival of their-style detections under our controls (`tables/controls_survival_*.csv`, `controls_counts_*.csv`)
Rule (PLAN §3): surviving fraction = (real detections not seen far away − the larger of the two scrambled-control means)
/ real detections; per template, then summed.

| noise rate = 0.79/day threshold | 2011 ETS | between ETS | 2012 ETS |
|---|---|---|---|
| real detections | 329 | 1,353 | 328 |
| station-day-shift control (mean of 8) | 80 | 643 | 71 |
| polarity-scrambled control (mean of 8) | 66 | 548 | 67 |
| **surviving fraction, our scrambled controls** | **0.70** | **0.45** | **0.74** |
| surviving fraction vs H&G-style noise templates | 0.31 | 0.01 | 0.39 |
| same with median + 6×MAD: ours / vs noise templates | 0.77 / 0.08 | 0.58 / 0.00 | 0.80 / 0.31 |

Per template, between ETS (noise rate = 0.79/day threshold), surviving fraction vs our scrambled controls:
T1 0.02, T2 0.70, T3 0.70, T4 0.50, T5 0.67.
- T1 between-ETS detections (445, a third of all) are indistinguishable from its station-day-shift control (442/draw).
- **Outside the teleseism masks (v2), survival barely changes** (between ETS: T2 0.64, T3 0.69, T4 0.48, T5 0.72). Masked
  earthquake waves do not explain the excess over our controls.

**Far/near amplitude** (raw figure `figures/controls_farnear_thr_rate079.png`; `tables/controls_calibration.csv`).
Far references are IU.COR, US.NEW, US.HAWA, UW.OMAK and TA.I02D, 300–700 km away.
- The five template events are not seen far away: median excess 0.49–1.14, far/near ratio 0.66–0.77.
- Random times are flagged 1.6–2.3 % of the time (far/near ratio 0.82–0.84).
- Real detections are flagged 0–17 % per template (4–6 % overall); scrambled-control detections ~2 %.
- The flag removes few detections. The far/near ratio separates the template events from random times only weakly
  (0.7 vs 0.83), so it is a weak referee for events this small.

## 5. Do between-ETS detections behave like events? (`tables/controls_clustering_*.csv`, `lfe_coincidence_*.csv`)
**Clustering in time, between ETS** (real vs mean of the scrambled-control draws):

| | T1 | T2 | T3 | T4 | T5 |
|---|---|---|---|---|---|
| Fano factor of daily counts, real / controls | 1.22 / 1.20–1.30 | 1.20 / 1.06–1.08 | 1.36 / 1.06–1.17 | 1.11 / 1.04–1.07 | 1.04 / 1.01–1.02 |
| fraction of gaps < 1 h, real / Poisson for that rate | 0.063 / 0.050 | 0.028 / 0.017 | 0.124 / 0.056 | 0.019 / 0.018 | 0.020 / 0.012 |
| Spearman ρ with tremor ≤ 100 km, real / controls | 0.05 / 0.00–0.03 | 0.13 / 0.00–0.02 | 0.06 / 0.00–0.02 | 0.07 / 0.01–0.03 | 0.03 / −0.03–0.00 |

- Between-ETS detections cluster slightly more than the controls (T2, T3) and follow tremor bursts weakly (ρ ≤ 0.13).
  That is a small signal, not an event-like catalogue.

**LFEs within ±60 s** (our SVI catalogue at P ≥ 0.5 plus our WA-NW catalogue; raw figure `figures/lfe_coincidence.png`):

| noise rate = 0.79/day threshold | 2011 ETS | between ETS | 2012 ETS |
|---|---|---|---|
| detections with ≥ 1 LFE within ±60 s | 14.3 % | 0.4 % | 28.7 % |
| polarity-scrambled detections | 7.0 % | 0.4 % | 13.4 % |
| random times | 6.5 % | 0.5 % | 11.1 % |
| mean LFEs within ±60 s, detections / random | 2.5× | 1.9× (1 T2 detection) | 3.8× |

- **In the ETS, T2 and T4 detections carry LFEs at 4.6–5.5× the random rate** (29–42 % have an LFE within 60 s). T3 does
  so in the 2012 ETS (4.6×). T1 does not (0.9–1.6×), which fits T1 being the noisiest template.
- **Between ETS, the detections carry LFEs at the random rate** (0.4 % vs 0.5 %).
- Scope: our between-ETS LFE catalogue is sparse (base rate 0.5 % per 2-min window). So "no LFE excess" means the
  between-ETS detections are not tied to the LFEs we catalogue. It does not rule out VLFEs without LFEs.

## 6. Caveats
- **Station sets.** They are not published, so we used all nine stations for every template. T1 has 24 channels
  (TA.D03D dead at its origin). The SNR-selected variant was empty for T1, T2 and T5 (no channel with template SNR ≥ 3),
  so it carries no information.
- **Inversion.** Whole-space Green's functions, not Crosson's 1-D model. The template events themselves relocate poorly
  (T1, T3, T5 by 49–104 km), so we cannot judge their stacked-inversion claim.
- **Our counts are half of theirs.** We reproduce their background rate and their T1, but not their between-ETS rate.
  The remaining difference sits in unreported choices (stations, declustering, run windows).
- **Bug fixed before any reported number.** One global running sum lost precision after ~0.9-m TA glitches on
  2012-06-07. Fixed with hourly sums; the stale outputs are kept in `data/mf_stale_den_v1/`.
- **Data.** 3.20 GB fetched from EarthScope (TA LH; SP2/DOSE BH decimated to 1 sps); 1.1 GB stored.

## Reverse-engineered reproduction (R2 / R2a / R2b / R2c, 2026-10-15 ~11:00 am–12:08 pm PT)
Rules: PLAN.md Amendments R2 (11:00), R2a (11:15, Table S1 targets), R2b (11:40, debugging), R2c (11:48, moment budget),
each written before its numbers. Table S1 (4,355 detections) was downloaded by the user (the publisher blocks scripts).
Which targets: every number below uses the Table S1 counts; the earlier equal-split R2 targets produced no number.
Raw figures first: `figures/r2_fig6_chosen.png`, `r2_fig7_chosen.png`, `r2_heldout_counts.png`, `r2_match.png`,
`r2b_pertrace.png`, `r2b_stack.png`, `r2b_lag_hist.png`, `r2_s1_season.png`, `r2_s1_controls.png`,
`r2_controls_daily.png`, `r2c_moment.png`.

### Plain answer
- **We cannot reproduce their catalogue detection by detection.** The best setting (velocity, all 9 stations, chosen on
  1 Jul–6 Sep 2011 only) matches only 12 % of their 2011-ETS detections within ±30 s, 5 % between ETS, 15 % in the 2012
  ETS. The held-out counts fail the ±25 % rule (between ETS 2,094 vs 2,923; 2012 ETS 530 vs 419).
- **It is not a timing or template bug on our side.** Our templates find themselves (top peak of their day). Their
  time stamps sit exactly 9 s (T1, T4, T5) and 12 s (T2) before the Table-1 origins, and our matches cluster there.
- **The matches are real, and the strongest detections agree.** Matches run 30–100× above chance. Their top 10 %
  (ranked by our CC at their times) are 77 % matched in the 2011 ETS, 74 % in 2012, 52 % between ETS.
- **The rest is a near-threshold tail that no pipeline reproduces.**
  - Our own processing variants, all fitted to the same counts, agree with each other only 24–67 %.
  - At their count level, our their-style noise templates fire as often as the real templates.
- **Part of it is a template/station mismatch, mainly T4 and T5.** Their T1/T2 windows stack to our templates (CC
  0.45–0.52, controls ~0). Their T4/T5 between-ETS windows do not (CC 0.03–0.06), though they are coherent with
  each other. So their T4/T5 between-ETS detections were selected by something other than our T4/T5 templates.
- **Verdict.** Both explanations hold, in different parts of the catalogue.
  - A strong core (top-ranked, ETS, T1/T2) is shared and looks like real repeats: LFE-tied, above our scrambled controls.
  - Most of the between-ETS catalogue is noise-dominated near-threshold detections. Among them, T4/T5 also show a
    station/processing mismatch we cannot fix without their station lists.
- **Moment.** Their detections carry no discrete-VLFE moment above what correlation-selected noise gives. The upper
  bound is ~6e16 N m for the 2011 ETS (< 1 % of the ETS's GNSS moment).

### 1. Calibration and held-out test (R2a; `tables/r2_calibration.csv`, `r2_heldout*.csv`, `r2_heldout_pass.json`)
- Thresholds per template were fitted to Table S1's 2011-ETS counts (34/299/82/8/62).
- Recall of their calibration-period detections within ±30 s, by variant:
  - V4 velocity 0.101 (chosen);
  - V6b 20-s declustering 0.090; V0 displacement 0.085; V5b window origin−9 s 0.084;
  - V3b 5 nearest stations 0.079; V3a SNR channels 0.051; V8 earthquake masks 0.035; V2 network normalisation 0.009.
  - No combination run was needed (only one group beat V0).
- Held out, chosen V4 (ours vs theirs):
  - between ETS 2,094 vs 2,923 (0.72, **fail**);
  - 2012 ETS 530 vs 419 (1.26, **fail**);
  - post-2012 289 vs 380.
  - Per template, only 1 of 5 is within ±25 % in each held-out period.
- Ratios of the other variants: between ETS 0.58–1.01, 2012 ETS 0.71–1.75. V8 (masks) gets between ETS right (1.01) but
  2012 ETS wrong (1.75). No variant passes both.

### 2. Debugging the overlap (R2b)
**Self-detection** (`tables/r2b_self.csv`):
- Every template, displacement and velocity, is the top peak of its own day: summed CC = channel count, next peak
  8–11 of 24–27.
- Their self-detections lie 9 s (T1, T4, T5) and 12 s (T2) before the Table-1 origin. Our CC at their exact stamp is
  negative, at +9/+12 s it is 1.0. So their stamp = our window start minus 9–12 s; ±30-s matching covers it.
- T3 has no self-detection in Table S1.
- T2 and T4 resemble each other (summed CC ~10 of 27 at each other's origin).

**Overlap vs rank, tolerance, lag** (`tables/r2b_overlap_rank.csv`, `r2b_overlap_tolerance.csv`, `figures/r2b_lag_hist.png`):

| | 2011 ETS | between ETS | 2012 ETS |
|---|---|---|---|
| their top 10 % (by our CC at their time) matched ±30 s | 0.77 | 0.52 | 0.74 |
| their top 25 % | 0.47 | 0.21 | 0.61 |
| all of theirs | 0.12 | 0.05 | 0.15 |
| our top 10 % in theirs | 0.29 | 0.18 | 0.21 |
| chance level ±30 s | 0.003 | 0.0015 | — |

- Ours − theirs lag histogram over ±600 s: one peak at 0 to +20 s, nothing elsewhere. There is no hidden time-convention
  shift.
- Unmatched Table-S1 detections fall in daytime (8 am–6 pm PT) 43 % of the time (uniform 42 %); matched ones 26 %.
- Our CC at their unmatched times: median 3.2–3.6 MAD, against our thresholds of 5.6–7.5 MAD.

**Self-consistency** (`tables/r2b_selfconsistency.csv`, 2011 ETS, ±30 s, all fitted to the same counts):

| fraction of A also in B | in V0 | in V4 | in V5b | in V6b | in V3b | in Table S1 |
|---|---|---|---|---|---|---|
| V0 | 1 | 0.49 | 0.67 | 0.96 | 0.27 | 0.10 |
| V4 | 0.49 | 1 | 0.48 | 0.49 | 0.24 | 0.11 |

- Changing displacement to velocity alone keeps only half the catalogue; changing the station set keeps a quarter.
- Near-threshold detections are not robust to processing in any implementation. Table S1 (~10 %) sits below even our
  most different variant (V3b, 24–27 %), which adds a real station/processing difference on top.

**Stacks of their detections, independent of our matched filter** (`tables/r2b_stack.csv`, `r2b_pertrace.csv`;
velocity, their time + 9 s):

| stack CC with our template | T1 | T2 | T3 | T4 | T5 |
|---|---|---|---|---|---|
| Table S1, 2011 ETS | 0.51 | 0.45 | 0.27 | 0.28 | 0.35 |
| Table S1, between ETS | 0.45 | 0.48 | 0.28 | 0.06 | 0.03 |
| Table S1, 2012 ETS | 0.47 | 0.52 | 0.23 | 0.04 | 0.08 |
| random same-day / station-day shift / polarity-scrambled | −0.06…0.08 | −0.07…0.06 | −0.05…0.07 | −0.09…0.06 | −0.05…0.06 |
| our own detections (selected by our template) | 0.78–0.88 | 0.90–0.91 | 0.83–0.85 | 0.54–0.71 | 0.76–0.83 |
| selection-matched noise template vs its own window | 0.78–0.84 | 0.84–0.86 | 0.79–0.86 | 0.64–0.83 | 0.78–0.84 |

- Their T1/T2 windows contain our template waveform well above every control.
- Stacked at their raw stamp (Δ = 0) the CC is negative (−0.07 to −0.34). With Δ = 9–11 s it is positive.
- Split by our CC at their times: the top tercile stacks at 0.59–0.68 (T1, T2), the bottom at 0.12–0.44.
- Semblance of their windows (0.07–0.20) is far above random (0.001–0.03), about as high as selection-matched noise
  (0.12–0.29). So their windows are mutually coherent, as any correlation-selected set is.
- T4/T5 between ETS are coherent among themselves but not like our T4/T5: their selection differed from ours.
- Seasons between ETS: T1, T2 and T3 stack CCs are steady (0.25–0.48), with the lowest in Dec–Feb.
- Reading: a stack resembling the selecting template does not prove real repeats. Selected noise gives 0.8 with its
  own window (R1, A6).

**Processing grid, calibration period only** (`tables/r2b_grid.csv`): all single changes from V4, each re-fitted to the Table-S1 counts.
- Code checks pass: G0 = V4 at 0.101, G10 displacement = V0 at 0.085.
- Recall within ±30 s:
  - G9 20-s declustering 0.111 (best);
  - G1 window origin−11 s 0.095; G7 zero-phase 2-pole 0.095;
  - G4 5 nearest stations 0.082; G6 causal filter 0.066; G3 horizontals 0.047;
  - G8 0.01–0.05 Hz 0.032; G2 Z only 0.022; G5 network normalisation 0.017.
- Every grid point recovers 0 % of T4 and ≤ 1 % of T5.
- So no simple processing change found here brings the catalogues together; the gap is in what is not published (station
  sets per template).

### 3. Tests on Table S1's own detections (`tables/r2_s1_*.csv`, `figures/r2_s1_controls.png`, `r2_s1_season.png`)
Thresholds are the chosen setting's (fitted to Table S1 counts).

| per template per day | 2011 ETS | between ETS | 2012 ETS |
|---|---|---|---|
| Table S1 | 2.16 | 1.63 | 2.00 |
| their-style noise templates (100, our data, same thresholds) | 2.48 | 1.89 | 2.80 |
| station-day-shift / polarity-scrambled controls | 0.28 / 0.36 | 0.26 / 0.27 | 0.25 / 0.29 |
| Table S1 excess over noise templates (mean; 16–84 % over the 100) | −0.15 (−0.92…+0.60) | −0.16 (−0.94…+0.53) | −0.40 (−1.32…+0.53) |
| Table S1 excess over our scrambled controls | 0.84 | 0.84 | 0.85 |

- **At Table S1's count level, a noise template fires about as often as a real template.**
  - Their stated ~0.6 false per template per day would need a threshold where the real templates give ~4× fewer
    detections than Table S1 holds.
  - Exception: T3 between ETS, 2.9/day vs 1.4/day from noise templates, excess 0.52 (0.13…0.86); T4 likewise 0.43.
- **Far/near:** Table S1 times are flagged "seen far away" 8–14 % of the time, against 2–4 % for random times. That is
  3–6× more earthquake (non-local) contamination than chance. The paper does not mask earthquakes.
- **LFEs within ±60 s** (SVI P ≥ 0.5 + WA-NW), against random times:
  - 1.9× in the 2011 ETS and 2.8× in 2012 (T1, T2 2.3–3.7×; T3, T5 ≤ 1.1×);
  - 1.24× between ETS (T1 5.3× on 3 of 88; T2–T5 0.6–1.2×).
- **Season:** daily counts between ETS against a 5–10 s microseism proxy (PGC, LRIV):
  - Table S1 all ρ = −0.29 (T2 −0.44), the same sign and size as the noise-template detections (−0.48);
  - T3 is the exception (+0.14): its Nov–Jan 2011 bulge sits in storm season, and our T3 has no such bulge.
  - The 20–50 s band shows no relation (|ρ| ≤ 0.2).
  - So T2's between-ETS rate behaves like noise-template detections, and T3's winter excess is storm-season timed.
- **Clustering between ETS:**
  - T2 and T3: Fano 1.9 and 2.3; gaps < 1 h 26 % (Poisson for that rate 11–15 %); scrambled controls 1.0–1.4 and 1–9 %.
  - Tremor correlation ρ is 0.21 for T2 and −0.07 for T3.
  - Clustering is stronger than in our scrambled controls, but noise bursts cluster too.

### 4. Moment budget of Table S1 detections (R2c; `tables/r2c_*.csv`, `figures/r2c_moment.png`)
- Method: each detection's amplitude is the least-squares scale of its template (displacement, fixed +9/+12 s) times
  the Table-1 M0.
- The pre-registered plain sum is dominated by a few huge windows: between ETS 7.4e18 N m from T2 with a median scale
  of 0.25. So the robust total (N × 5–95 % trimmed mean, post hoc) is the headline.

| N m, five templates | 2011 ETS | between ETS | 2012 ETS |
|---|---|---|---|
| Table S1, robust total (bootstrap 16–84 %) | 5.8e16 (4.9–7.1e16) | 2.4e17 (2.3–2.6e17) | 2.5e16 (2.1–3.0e16) |
| random same-day times, same counts (20 draws) | −0.1e16 (−0.4…+0.2e16) | 0.2e16 (−0.6…+0.8e16) | −0.1e16 (−0.4…+0.2e16) |
| selection-matched noise (10 noise templates, top-N, own window) | 10.9e16 (9.3–12.4e16) | 7.0e17 (5.7–8.3e17) | 10.7e16 (8.2–13.2e16) |
| net = Table S1 − selection-matched noise | −5.0e16 (−6.6…−3.5e16) | −4.6e17 | −8.2e16 |
| thrust-synthetic scaling (cross-check), robust | 2.1e17 | 5.3e17 | 6.8e16 |

- Table S1 windows hold more template-shaped amplitude than random windows. T2 has 72–77 % positive scales; T4/T5
  between ETS are ~50 % positive, i.e. random.
- They hold less than windows that noise templates pick for themselves. So no discrete-VLFE moment stands above the
  noise selection effect.
- The comparison is not exactly like for like: their windows were picked by their template, not ours, which dilutes the
  projection. So 5.8e16 (2011 ETS) and 2.5e16 (2012 ETS) are best read as **upper bounds**.

**Three-way split per ETS** (SVI LFE numbers from S13; GNSS not measured by us for 2011/2012, since the S1 tables start
in 2014; literature Mw 6.5–6.8 = 7e18–2e19 N m):

| N m | 2011 ETS | 2012 ETS |
|---|---|---|
| LFE moment (catalogued SVI LFEs × 4.5e11) | 6.5e14 (1,445 LFEs) | 1.3e15 (2,864) |
| LFE-locked slow moment (S13, 20–50 s) | 7.7e15 (± 1.4e15) | 1.28e16 (± 0.2e16) |
| discrete VLFE moment, Table S1 (upper bound; net above noise ≤ 0) | ≤ 6e16 | ≤ 2.5e16 |
| GNSS (literature) | 7e18–2e19 | 7e18–2e19 |

- In shares of GNSS: LFE moment ~0.003–0.02 %, LFE-locked slow moment ~0.04–0.2 %, discrete VLFEs ≤ 0.3–0.9 %, and
  probably near 0 net of noise.
- Caveats:
  - Template scaling is biased upward on any correlation-selected window (noise projects positive amplitude).
  - Noise share as above.
  - 20–50 s band only: smooth slip at longer periods is invisible.
  - The whole-space synthetic scale is uncertain by ×2: synthetic/Table-1 M0 is 0.8–2.3 for the five templates.
  - Areas differ: H&G templates span 47.8–48.4 N, our LFEs are SVI families.
  - The LFE moment uses one median magnitude.

### 5. Caveats and what would settle it
- Their per-template station/channel sets are still unknown; the T4/T5 mismatch points there.
- Asking the authors for the station list per template would settle T4/T5, but the project rule is no author contact.
- The noise templates are ours (100 random 2011-ETS windows). Their own noise background was quoted as ~0.6 per
  template per day, but at the count level of Table S1 ours fire at 1.9–2.8.
- One code test printed old-target held-out ratios before the automatic R2 choice (disclosed in NOTES); the choice was
  automatic and unchanged.

## Individual between-ETS candidates (R3, 2026-10-15 ~12:28–12:45 pm PT)
Rules: PLAN.md Amendment R3 (written ~12:28 pm PT, before any R3 number). The ETS-window run is a POST-HOC positive
control (same code and rules, `R3_TAG=ets`), not in the amendment. Merlin not consulted (flagged).
Raw figures first: `figures/r3_overview.png` (ratio vs scan z coloured by our CC; best-node map coloured by z;
held-out CC; pass fractions vs controls), `r3_aggregate.png`, `r3_wave_cand000–002.png` (top 3 between-ETS),
`r3_nearmiss_wave_cand007/049/116/128.png` (the 4 between-ETS candidates failing only one test), and for the positive
control `r3_ets_overview.png`, `r3_ets_aggregate.png`, `r3_ets_wave_cand*.png` (8 real-looking ETS events).
Tables: `tables/r3_candidates.csv`, `r3_events.csv` (every candidate, all columns), `r3_pass_counts.csv`,
`r3_real_looking.csv` (empty), `r3_aggregate.csv`, `r3_pw_baserate.csv`; ETS: `r3_ets_*`.

### Plain answer
- **No real-looking VLFE between ETS: 0 of 130 candidates pass all five tests** (top 30: 0). So there are none to
  split into with/without LFEs.
- **The tests do find VLFEs when they exist.** The same pipeline on Table-S1 detections inside the 2011 and 2012 ETS
  finds **8 of 115** (top 30: 4). All 8 are T2, Mw 3.1–3.5, within 1–6 km of an LFE family, all with catalogued LFEs.
  One is the scan-pilot-3 event (2012-09-19 19:18 PDT), found independently.
- Between-ETS 0/130 vs ETS 8/115: Fisher p = 0.002. At the ETS rate we would expect ~9 between ETS.
- **What the between-ETS candidates are instead:** mostly earthquake waves and waves from outside the network.
  - 101 of 130 sit in a catalogued-earthquake window (78 %, vs 52 % for random same-day times).
  - Of the 29 outside earthquake windows, half (15) have far stations as loud as the near ones; the rest are
    unlocated, single-station or both (no candidate passes a, b, c and f together).
- **Selection-matched stack (pre-registered R3.3):** the held-out stations (PGC, SNB, PFB, YOUB) show nothing.
  - Mean held-out CC is −0.018 for the top 30 and −0.003 for all 130.
  - That is the 27th–28th percentile of correlation-selected noise templates (97.5th percentile: 0.25 and 0.11).
  - So the set is consistent with correlation-selected noise.
  - Caveat: this test is weak. The template events themselves are barely above noise at these stations.

### Candidate set (R3.1)
- Table S1 has 2,886 detections between 8 Sep 2011 and 28 Aug 2012 (±2 d from their ETS windows removed).
- 148 match one of our V4 detections within ±30 s after the 9/12-s offset: **130 distinct events**.
  - T2 100, T3 28, T1 1, T5 1.
  - Ranked by our CC (5.6–8.5 MAD).

### Per-event tests (R3.2; `tables/r3_pass_counts.csv`)
| passes | candidates (130) | top 30 | scrambled control (650) | selection-matched noise (1,300) | ETS positive control (115) |
|---|---|---|---|---|---|
| (a) far/near ratio6 ≤ 0.551 (n ≥ 2 refs) | 30 (23 %) | 3 | 45 % | 32 % | 52 % |
| (b) not a travelling plane wave | 98 (75 %) | 23 | 90 % | 88 % | 88 % |
| (c) located (z ≥ 5), inside footprint or ≤ 20 km from a family | 45 (35 %) | 11 | 6 % | 9 % | 41 % |
| (e) no catalogued earthquake window | 29 (22 %) | 4 | 48 % | 48 % | 41 % |
| (f) not single-station (z_loso ≥ 3 and ≥ 0.5 z) | 59 (45 %) | 12 | 14 % | 19 % | 45 % |
| **all five = real-looking local VLFE** | **0** | **0** | 0 | 1 (0.08 %) | **8 (7 %)** |

- Classes, first match wins: catalogued earthquake 101; travelling wave 17; not locatable (z < 5) 10; single-station 2.
- Earthquake windows, candidates vs random same-day times:
  - teleseismic surface waves (global M ≥ 4.5) 72 % vs 47 %;
  - body waves (M ≥ 5) 22 % vs 8 %;
  - M ≥ 3 within 2,500 km 12 % vs 4 %;
  - regional M ≥ 2.5 3 % vs 1 %.
- Their strongest between-ETS detections are 1.5–3× more often inside earthquake wavetrains than chance. The scan sees
  them clearly (z ≥ 5 for 72 %, vs 22 % at random times), but as waves the far stations see too (ratio6 median 1.06
  for the located ones).
- The travelling-wave flag fires at 9 % of random between-ETS times (`r3_pw_baserate.csv`). This is the background
  rate of coherent long-period wavefields.
- The plane-wave base rate and the 0/650 scrambled-control count show the rules are not loose. The ETS 8/115 shows they
  are not too strict.

**Between ETS, all five tests passed by none. Four candidates fail only one test** (post-hoc look, raw figures
`r3_nearmiss_wave_cand*.png`):
| cand | time (PT) | template | fails | what the waveforms show |
|---|---|---|---|---|
| 7 | 2012-05-14 05:30:35 PDT | T3 | (a) ratio 1.26 | far stations as loud as near; travelling long-period wave |
| 49 | 2012-04-10 00:02:02 PDT | T2 | (a) ratio 1.34 | continuous long-period wavefield on every station incl. 300–500 km; no pulse. In the April 2012 SVI minor episode (1 LFE P ≥ 0.75, 8 candidates within ±10 min) |
| 116 | 2012-05-09 13:47:27 PDT | T2 | (f) z_loso 1.4 of 6.7 | noise level everywhere (25–55 nm); one pulse at PFB only |
| 128 | 2011-12-10 20:44:03 PST | T2 | (a) ratio 1.53 | 300–700 nm long-period wavefield on all stations incl. far set; a separate transient at NLLB |

None of the four looks like the ETS events. Those show a pulse that is largest at the nearest stations (VGZ and LRIV,
20 km), decays with distance, and is absent at the far stations (`r3_ets_wave_cand000.png`).

**LFE status (d) of the 130:** quiet 99, candidates/tremor only 25, with LFEs 6 (none of them passes).

### The 8 real-looking ETS events (positive control; `tables/r3_ets_real_looking.csv`)
| time (PT) | Mw (template / scan) | best node (lat, lon, depth) | nearest family | LFEs P ≥ 0.75 ±10 min |
|---|---|---|---|---|
| 2011-08-19 21:37:55 PDT | 3.5 / 3.3 | 47.98, −123.31, 38 km | Olympic_2010:750, 1 km | 2 SVI, 9 PNSN tremor |
| 2011-08-25 18:25:51 PDT | 3.1 / 3.3 | 48.25, −123.72, 36 km | SVI 13, 1 km | 7 SVI + 7 WA-NW |
| 2011-08-26 03:44:21 PDT | 3.3 / 3.3 | 48.25, −123.58, 37 km | Olympic_2010:907, 1 km | 5 SVI |
| 2012-09-14 20:49:19 PDT | 2.6 / 3.5 | 48.34, −123.99, 34 km | SVI 22, 4 km | 8 SVI |
| 2012-09-18 05:02:57 PDT | 3.5 / 3.3 | 48.34, −123.85, 35 km | SVI 28, 1 km | 6 SVI |
| 2012-09-18 07:15:08 PDT | 3.1 / 3.2 | 48.16, −123.45, 38 km | Olympic_2010:849, 1 km | 9 SVI + 4 WA-NW |
| 2012-09-19 19:18:43 PDT | 3.4 / 3.5 | 48.25, −123.45, 38 km | Olympic_2010_armB:419, 1 km | 16 SVI + 4 WA-NW (= scan-pilot-3 event) |
| 2012-09-19 19:21:55 PDT | 3.5 / 3.1 | 48.16, −123.58, 36 km | SVI 3, 6 km | 14 SVI + 4 WA-NW |
- Scrambled-control and noise-template rates in the ETS run: 2/575 and 1/1,150.
- The 2011-ETS dates fall in H&G's own ETS window, so these overlap their ETS result. That is what a positive control
  should do.

### Aggregate test against selection bias (R3.3; `tables/r3_aggregate.csv`, `figures/r3_aggregate.png`)
| | between ETS top 30 | between ETS all 130 | ETS top 30 (post hoc) | ETS all 115 (post hoc) |
|---|---|---|---|---|
| G = mean held-out-station CC with the template event | −0.018 | −0.003 | 0.045 | 0.021 |
| random same-day times, median (97.5 %) | 0.015 (0.037) | 0.000 (0.018) | −0.001 (0.032) | 0.007 (0.022) |
| selection-matched noise templates, median (97.5 %) | 0.023 (0.25) | 0.013 (0.11) | 0.023 (0.19) | 0.012 (0.10) |
| percentile of G among selection-matched | 27 % | 28 % | 63 % | 60 % |
| H = H&G-network stack CC (selection-dependent) | 0.54 | 0.70 | 0.71 | 0.77 |
| H, selection-matched noise median | 0.72 | 0.83 | 0.69 | 0.80 |
- By the pre-set reading, the between-ETS set is consistent with correlation-selected noise.
  - G sits below the random-time median.
  - H (their own-network stack) is below what selected noise gives itself. H cannot discriminate (R1/R2b).
- Power, post hoc:
  - The template events are barely above noise at the held-out stations (window/pre-event RMS 0.4–1.7).
  - The real neighbouring events T2 and T4 correlate at only 0.36 there.
  - The ETS set reaches G 0.045 (top of random times), but that is within the selection-matched range.
  - So G could only have caught template-sized repeats. The per-event tests (a–f), with the ETS control at 8/115,
    carry the answer.

### Caveats
- The scan network (PGC, VGZ, SNB, NLLB, PFB, YOUB, LRIV) locates sources best inside the SVI footprint. H&G's T1
  (47.8 N) and T5 (−123.2) are near its edge, but those two templates contribute only 2 of the 130 candidates.
- Test (e) uses physical arrival windows (surface waves, body waves, regional), not the pilot-1 m45 mask (82 % of
  candidates, 69 % of random times; column only).
- The class label "single-station" is given only to located events (z ≥ 5). Unlocated ones are labelled "not
  locatable". The pass columns are as pre-registered.
- Candidates are only the Table-S1 detections our run also finds (5 % of theirs between ETS). Their other 95 % include
  weaker detections our run does not reproduce. Those are even less likely to pass. R2b already showed they sit at the
  selection-matched noise level.
- Data: all existing (2011 and 2012 far references on disk). PGC Oct–Dec 2012 lies outside the window.
