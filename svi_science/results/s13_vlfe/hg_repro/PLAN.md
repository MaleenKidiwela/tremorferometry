# hg_repro PLAN — literal reproduction of Hutchison & Ghosh (2019), then our controls

Written 2026-10-15 ~2:15 am PT, before any number from our data. Paper: Hutchison, A. A. & Ghosh, A. (2019), JGR Solid Earth
124, 554–565, doi:10.1029/2018JB016138 (`notes/papers/Hutchison_Ghosh_2018_JGR_VLFE.pdf`, read in full, pp. 554–565).
Table 1 transcribed in `../tables/hg2018_table1.csv`. The paper numbers below are the known answer; our code is wrong until
it reproduces them or a stated difference explains the gap.

## 1. Known-answer targets (from the paper)
Text:
- T-a. 2011 ETS catalogue (run 1 Jul–1 Oct 2011, five 2011 templates): **1,394 new detections** (§3.1) ≈ 15.2/day.
- T-b. Threshold ≈ 6×MAD of all summed CC; stated false rate **~0.6 per template per day** (§2.3).
- T-c. 18-month catalogue: activity in all five templates between ETS, above the noise-template background (§3.3, Fig. 6).
  3 of 5 templates rise again in the 2012 ETS; 2 show no rate increase in the inter-ETS or 2012 ETS (§3.3; Fig. 7 labels:
  2012 ETS "strong" T2, "mild" T1, "negligible" T3 and T5; 2011 ETS "strong" T1 and T5, "mild" T2 and T4).
- T-d. Stacks of detections invert to solutions close to the templates (moment tensor, location, VR, CLVD), VR > 40 % for
  every stacked inversion; also separately for ETS and inter-ETS stacks (§2.3, §3, §3.3). No numbers are given.
- T-e. Fig. 2 stack sizes (period not stated): T1 SP2 BHE 150 traces, T2 B05D BHN 210, T4 A04D BHZ 180 (2011 row).
Digitised (APPROXIMATE; scripts `scripts/hg_digitize_fig6.py`, `hg_digitize_fig7.py`, pixel calibration in the docstrings):
- T-f. Fig. 6 (all five templates summed, 18-month run): mean **8.5/day** (1 Jul 2011–1 Dec 2012); 2011 ETS **10.7/day**;
  inter-ETS (6 Sep 2011–30 Aug 2012) **8.2/day**; 2012 ETS **10.3/day** → ETS / inter-ETS ≈ **1.3**. Noise-template
  background (cyan "BGN", all five summed) **3.9/day** ≈ 0.8 per template per day. Table `tables/paper_fig6_period_means.csv`.
  (Fig. 6 gives ~776 detections for 1 Jul–1 Oct 2011, vs 1,394 in the text's separate 2011 ETS run: the paper's two runs
  differ by ~1.8× in that window. We compare with both.)
- T-g. Fig. 7 per-template rate in each period relative to the template's own 18-month mean (`tables/paper_fig7_period_fractions.csv`):

| period | T1 | T2 | T3 | T4 | T5 |
|---|---|---|---|---|---|
| 2011 ETS (23 Jul–6 Sep 2011) | 2.14 | 0.99 | 0.63 | 1.06 | 1.46 |
| inter-ETS (6 Sep 2011–30 Aug 2012) | 0.88 | 1.01 | 1.08 | 0.97 | 0.93 |
| 2012 ETS (30 Aug–11 Oct 2012) | 1.36 | 1.09 | 0.57 | 0.82 | 0.91 |

  Note: the paper defines the inter-ETS as 6 Sep 2011–30 Aug 2012 (§2.1) and as 1 Oct 2011–30 Aug 2012 (§3.3). We report both
  (digitised rates differ by < 0.05).

Known-answer pass rules (set now): (K1) our pooled rate 1 Jul 2011–1 Dec 2012 within ×2 of 8.5/day; (K2) pooled ETS /
inter-ETS rate ratio in 1.0–1.8; (K3) our H&G-style noise-template rate within ×2 of 0.6–0.8 per template per day; (K4) T1
has the largest 2011-ETS / inter-ETS ratio of the five. If K1 or K3 is off by > ×3, debug before steps 4–5.

## 2. Paper steps: followed / changed / added
| # | paper step | status | detail / reason |
|---|---|---|---|
| 1 | three-component broadband UW, CN, TA (Fig. 1: CN.NLLB, CN.VGZ, UW.LRIV, UW.DOSE, UW.GNW, UW.SP2, TA.A04D, TA.B05D, TA.D03D) | followed | on-disk raw first; missing 2011–2013 days fetched (TA LH 1 sps; SP2 BH decimated to 1 sps in memory, only the 1-sps file stored); budget ≤ 5 GB, logged in `logs/fetch.log` |
| 2 | station set per template = the stations of its best CMT solution | **changed** | the per-template sets are not published (only Fig. 2 shows T1 uses SP2, T2 B05D, T4 A04D). Primary: all nine Fig.-1 stations live at the template time, 3 components. Variant V-snr (column): per template, only components whose 90-s template RMS ≥ 3× the RMS of the 30 min before the origin |
| 3 | remove response; band-pass 0.02–0.05 Hz | followed | response → displacement (the paper does not say; the CMT inversion uses displacement); 1 sps (the paper does not state a rate; 1-s CC step); 4-pole zero-phase Butterworth 0.02–0.05 Hz |
| 4 | template = 90-s window (CMT window) | followed | window [origin, origin + 90 s] at every station (moveout kept, as in a CMT window); the paper does not give the window start |
| 5 | SEC-C sliding normalised CC, summed over stations/channels, 1-s step | followed | our FFT implementation of the same statistic (not the SEC-C code); sum over LIVE channels (literal "summed"); n_live column |
| 6 | detection = time windows above threshold | followed + declustering stated | the paper gives no declustering: one detection per above-threshold excursion; peaks closer than 90 s merged (keep max); positive CC only (negative-CC excursions kept as a polarity column) |
| 7 | threshold: 100 random noise windows (2011 ETS, not detected) as noise templates, run over 1 Jun 2011–1 Dec 2013 with the template's station configuration; threshold = mean over the 100 of the 3rd-highest summed CC | followed | windows drawn from 23 Jul–6 Sep 2011 with all nine stations live, ≥ 180 s from the five template origins and from every positive summed-CC peak of the template above its 6×MAD level (stands in for "not detected" before the threshold exists; one pass, no iteration); "3rd highest" = 3rd-highest declustered peak including the self-match (primary, A); variant B = 3rd highest excluding the self-match (column). 6×MAD of all summed CC reported too |
| 8 | run 1 Jul 2011–1 Dec 2012 (also stated 1 Jun 2011–31 Dec 2012) | followed | we run 1 Jun 2011–1 Dec 2013 once (the noise period) and count in the paper's windows |
| 9 | 2011 ETS catalogue (1 Jul–1 Oct 2011) | followed | count in that window from the same run (the paper's separate run may have differed; T-a vs T-f) |
| 10 | stack detections per template per station-channel; invert by grid-search CMT (0.1° × 0.1° × 5 km, Crosson 1-D GF, 90-s window, 0.1-s step), confirm on a 0.025° × 1 km grid; VR, CLVD | **changed** | whole-space moment-tensor Green's functions (`../scripts/s13_green.py`, free surface ×2; no 1-D layered code here); coarse grid 0.1° × 0.1° × 5 km, 46.8–49.6 N, 125.5–121.5 W, depth 20–70 km; origin-time shift ±10 s at 1 s; fine grid 0.025° × 2 km around the best node; deviatoric MT by linear least squares per node; VR and CLVD. The template events are inverted with the same code, so stack-vs-template is like-for-like. Stacks: all detections, ETS only, inter-ETS only |
| 11 | spatiotemporal comparison with tremor (50 km of the stacked solution; 100 km for the summed view) | followed (plots only) | PNSN `catalogs/pnsn_tremor_cascadia_full.csv` |
| A1 | — | **added**: station-day-shift control | each station's continuous record circularly shifted by a distinct whole number of days (station k: k × 41 d mod period, draw-specific offset), all three components together; real templates; same threshold; 8 draws |
| A2 | — | **added**: polarity-scrambled control | random ±1 per station-component of the real template (not all equal), real data, same threshold; 8 draws |
| A3 | — | **added**: far/near amplitude | far references on disk (LHZ, `../data/raw_ref2011`, `raw_ref2012`: IU.COR, US.NEW, US.HAWA, UW.OMAK, TA.I02D; 300–700 km). Per detection: ref_excess = RMS(far Z, [t−60, t+330] s) / RMS(far Z, surrounding ±1 h excluding that window); ref_ratio = median far RMS / median near-Z RMS in [t, t+90]. **Flag "seen far away" if ref_excess ≥ 3 at ≥ 2 live references** (a local Mw ≤ 3.7 VLFE should not reach 3 noise units at 300 km: whole-space decay ≤ 0.2 at 100 km, measured GNW/PGC ≤ 0.35). Calibrated on the 5 templates themselves and on random times (base rate). GNW ratio a column for T2–T5 (GNW is a template station) |
| A4 | — | **added**: temporal clustering | inter-ETS detections per template, real vs A1/A2 control detections: Fano factor of daily counts, fraction of inter-event times < 1 h, Spearman correlation of daily counts with PNSN tremor counts within 100 km |
| A5 | — | **added**: earthquake masks as a column only | `../data/scan/masks2_2011.npz`, `../data/scan_2012/masks2_2012.npz` (v2), read-only; the literal run uses no mask |

## 3. Survival under our controls (set now)
Per template and period (2011 ETS, inter-ETS, 2012 ETS):
- N_real = H&G-style detections; N_A1, N_A2 = mean detections per draw of each control in the same period (by output time).
- Detections flagged by A3 are removed from both real and control counts before the subtraction.
- **Surviving fraction = max(0, N_real,unflagged − max(N_A1, N_A2)_unflagged) / N_real.** Reported per template, pooled,
  ETS vs inter-ETS; the per-draw spread is shown as the range (no error bars in figures).
- A4 answers whether inter-ETS detections cluster like real events (tremor-correlated, bursty) or like the controls.

## 4. Order, compute, rules
1. Fetch (≤ 5 GB) → prep 1-sps displacement (reuse `../data/disp1sps_v3` files read-only where they exist; new files in
   `data/disp1sps/`, same v3 recipe) → `prep` marker.
2. Templates + matched filter (real, 100 noise templates × config, A1, A2), CPU FFT, ≤ 4 workers, OMP_NUM_THREADS=1.
3. Known-answer check (K1–K4). 4. Stack-and-invert. 5. Controls/survival. RESULTS.md, raw figures first, no error bars.
Stage > 20 min → nohup + DONE marker in `hg_repro/`, stop and report the resume command. Never delete files.
Merlin not consulted (coordinator's spec is explicit; Fable quota) — flagged for the coordinator.

## Amendment 1 (2026-10-15 ~2:30 am PT, before any matched-filter or inversion number)
- A6 (added, step 4): **selection-bias control for stack-and-invert.** A stack of windows chosen because they correlate with a
  template resembles that template even if every window is noise. So the same stack-and-invert is run on (i) the detections
  of 5 noise templates (j = 0–4, each under the T-config it was thresholded with, self excluded) and (ii) 5 sets of random
  times (as many as the real detections). The paper's test (VR > 40 %, solution close to the template) is discriminating
  only if the noise-template stacks fail it. Reported side by side; the noise templates are also inverted as single windows.
- GF source duration: τ (erf step) chosen from {2, 5, 10, 20} s as the value with the highest mean VR over the five
  template events, then fixed for every stack.
- Step-1 check made at ~2:25 am PT (template RMS / RMS of the 30 min before, CN stations): single-channel template SNR is
  0.5–4 (VGZ horizontals ~300 nm noise). Recorded only; the literal primary config keeps all live channels.

## Amendment 2 (2026-10-15 ~2:50 am PT, before any matched-filter number)
Reading §2.3 against Fig. 6: with the literal threshold (mean 3rd-highest noise-template peak over 1 Jun 2011–1 Dec 2013),
a noise template exceeds the threshold only ~2–3 times in 2.5 years (≈ 0.003/day) by construction. The paper instead quotes
~0.6 false detections per template per day (its 6×MAD Gaussian estimate) and plots a noise background of 3.9/day for the five
templates (≈ 0.8 per template per day). So the literal reading and the paper's own false rate cannot both hold. Two
alternative thresholds are therefore pre-registered as columns, each run through the same detections, controls and K1–K4:
- thr_6MAD = median + 6 × MAD of the template's own summed CC over the run (the paper's stated near-equivalent);
- thr_rate079 = the summed CC at which the 100 noise templates give 0.79 detections per template per day (self excluded,
  1 Jun 2011–1 Dec 2013), i.e. calibrated to the paper's Fig.-6 background.
The literal thrA stays the primary. If thrA fails K1/K3 and an alternative passes, the known-answer comparison is reported
for that alternative as "reproduced under interpretation X", and our controls are reported for all three thresholds.
Stack-and-invert is run for thrA and for the alternative that best matches K1.

## Amendment R2 — reverse-engineered reproduction (2026-10-15 ~11:00 am PT, before any new number from our data)
User request (morning 10-15): "can't we just reverse engineer to reproduce what they got". Goal: tune the unreported
choices on a calibration period only, then predict held-out periods and test against their counts.

**R2.0 Supporting information.** Not reachable without a browser login/challenge: agupubs/onlinelibrary landing page,
full-xml and every `downloadSupplement` URL return HTTP 403 (Cloudflare challenge); the open copies (ETH Research
Collection hdl 20.500.11850/328710, sonar.ch 140124) hold no supplement files; figshare has none. Log:
`notes/papers/HG2019_SI/SOURCE_LOG.md`. So Table S1 is not used and step 3 (detection-by-detection match) cannot be run.
Targets are the digitised Figs 6–7 (paper numbers only, no data of ours).

**R2.1 Targets** (`tables/r2_targets.csv`). Pooled (Fig. 6, five templates): pre-ETS 1–23 Jul 2011 148; 2011 ETS
481; inter-ETS 6 Sep 2011–30 Aug 2012 2,958 (1 Oct–30 Aug: 2,812); 2012 ETS 431; post 11 Oct–1 Dec 2012 392.
Per template the split of the Fig-6 total between templates is not recoverable: a non-negative fit of the Fig-6 daily
rate on the five Fig-7 cumulative shapes gives totals 574/725/1,436/213/1,411 with bootstrap ranges of a factor ≥ 2, and an
equal split fits nearly as well (relative residual 0.267 vs 0.240). So:
- C1 (primary): equal 18-month totals (882 per template) × the template's Fig-7 fraction in each period. 2011-ETS targets
  T1–T5 = 164, 76, 48, 80, 111 (sum 479).
- C2 (column): the non-negative-fit totals × Fig-7 fractions (2011 ETS: 107, 62, 78, 19, 178).
- C0 (column): one common threshold level for all templates, median + k × MAD of each template's own summed CC (the
  paper's "≈ 6 × MAD" reading), k tuned so the pooled 2011-ETS count = 481.

**R2.2 Calibration vs held-out data.** Calibration = 1 Jul–6 Sep 2011 only. Thresholds: per template, the summed-CC level
at which our declustered detections in the 2011 ETS (23 Jul–6 Sep 2011, self-detection excluded) equal the target.
Variant choice (thresholds always re-fitted per variant): primary score = |ln(our pooled pre-ETS count / 148)| on
1–23 Jul 2011 (inside the calibration run but not used by any threshold); if two variants are within 0.10 of each other,
the higher Pearson r between our pooled daily counts and the digitised Fig-6 daily rate on 1 Jul–6 Sep 2011 wins.
Held out (never looked at before the choice is logged in NOTES.md): inter-ETS, 2012 ETS, post-2012, and the noise-template
false rate at the calibrated thresholds (paper 0.6 text / 0.79 Fig. 6 per template per day).
**Pass rule (set now): pooled held-out counts within ±25 % of the paper for BOTH inter-ETS (2,958) and 2012 ETS (431).**
Per template (vs C1 targets) reported as a secondary column (how many of five within ±25 %), not a pass criterion.

**R2.3 Pre-listed processing variants** (where the paper is vague; each is one change from V0):
- V0: R1 code — 9 stations × 3 components, live channels, plain sum, displacement, window [origin, origin + 90 s],
  1-s step, 90-s declustering.
- V1: per-channel normalisation of the sum: S / n_live (mean CC over live channels; changes behaviour across data gaps).
- V2: network normalisation: one normalised CC over all live channels concatenated (amplitude-weighted).
- V3a: channel set = channels with template SNR ≥ 1.5 (90-s template RMS / RMS of the 30 min before).
- V3b: channel set = the 5 stations nearest the template epicentre (15 channels).
- V4: velocity instead of displacement (time derivative of the same band-passed record; templates likewise).
- V5: alignment: template window [origin − 30 s, origin + 60 s].
- V6: declustering 30 s instead of 90 s.
- V8: detections inside our earthquake masks (v2) removed (the AGU FM15 abstract of the same work removed ANSS
  earthquake windows; the paper does not say). Thresholds re-fitted after removal.
- V7: one combination run, made only if singles from different groups (normalisation V1/V2, channel set V3a/V3b,
  others V4/V5/V6/V8) each beat V0 on the score: the best of each such group together.
For each variant also logged (calibration period only): calibrated thresholds in MAD units, pre-ETS count, daily r.

**R2.4 At the chosen setting** (same code, chosen thresholds fixed): 100 noise templates (the R1 noise times, each
noise window under each template's channel set), station-day-shift and polarity-scrambled controls (8 draws each, R1
draws), far/near flag, survival (PLAN §3 rule), LFE coincidence (±60 s, SVI P ≥ 0.5 + WA-NW, vs random and scrambled),
ETS vs inter-ETS. **Inter-ETS excess above noise** = (N_real − mean N_noise-template) / N_real, with uncertainty = the
16–84 % range of that quantity over the 100 noise templates (plus √N_real as a column). Raw figures first, no error bars.
Compute: ≤ 8 workers, OMP_NUM_THREADS=1, < 100 GB; new outputs in `data/r2/`, `tables/r2_*`, `figures/r2_*`.
Merlin not consulted (coordinator's spec is explicit) — flagged.

## Amendment R2a — Table S1 targets (2026-10-15 ~11:15 am PT, before any new number from our data)
The user downloaded the supplement from the publisher: `notes/papers/HG2019_SI/HG2019_TableS1_detections.csv` (4,355
detections 1 Jul 2011–29 Nov 2012; year, month, day, seconds of day UTC, template 1–5) and the SI text (false rate
~0.6 per template per day). The R2 queue was stopped before any variant output existed (killed by PID ~11:14 am PT; no R2
number was produced) and is restarted under this amendment. Everything below uses the Table S1 targets; the R2
equal-split targets are kept only as a comparison column (`C1_equal`).

**R2a.1 Targets (Table S1 counts; `tables/r2_targets_S1.csv`).** Per template T1/T2/T3/T4/T5 (pooled):
pre-ETS 1–23 Jul 2011 5/89/26/5/23 (148); 2011 ETS 34/299/82/8/62 (485); inter-ETS 6 Sep 2011–30 Aug 2012
88/1,417/1,036/84/298 (2,923); 1 Oct 2011–30 Aug 2012 (2,786); 2012 ETS 27/281/59/8/44 (419); post 11 Oct–1 Dec 2012
13/133/172/22/40 (380). (T2 + T3 = 82 % of the catalogue; the R2 equal split was far off.)

**R2a.2 Catalogue facts used to extend the variant list** (properties of their catalogue, not of our data):
- Within a template, detections are as close as 22–24 s (T1 27 s; T4 ≥ 490 s): their declustering is much shorter
  than our 90 s → add **V6b: declustering 20 s**.
- Their detections of the template events themselves sit 9 s (T1, T4, T5) and 12 s (T2) before the Table-1 origin
  (T3: none within 5 min) → their time stamp ≈ window start ≈ origin − 9 s → add **V5b: template window [origin − 9 s,
  origin + 81 s]**.

**R2a.3 Calibration and choice (calibration period 1 Jul–6 Sep 2011 only).** Per template, the threshold matching the
Table-S1 2011-ETS count (self-detection excluded). Variant choice, now that their detections exist:
primary = fraction of Table-S1 detections in 1 Jul–6 Sep 2011 that our catalogue matches within ±30 s (pooled over
templates); within 0.03 of the best → smaller |ln(our pre-ETS count / 148)|; then higher daily r with Table-S1 daily counts.
V7 rule unchanged (groups: normalisation V1/V2; channel set V3a/V3b; others V4/V5/V5b/V6/V6b/V8).

**R2a.4 Held-out test (unchanged periods).** Pass = pooled inter-ETS (2,923) and 2012-ETS (419) counts both within ±25 %.
Secondary: templates within ±25 % per period. **Detection-by-detection match** of our calibrated catalogue vs Table S1,
per template and period, at ±30 s and at ±5 s (±5 s after removing one constant offset = the median of
ours − theirs over the ±30-s matches in the calibration period; the raw ±5 s match is a column): fraction of theirs we
recover, fraction of ours in theirs. Unmatched vs matched Table-S1 detections compared on: our summed CC at their time
(largest of our declustered peaks within ±45 s, in MAD units), Pacific hour of day, period, month.

**R2a.5 Tests on Table S1's own detections** (thresholds and runs of the chosen setting):
- Scrambled controls and noise templates: per template and period, Table-S1 count vs station-day-shift (8 draws),
  polarity-scrambled (8 draws) and their-style noise-template (100) detections at the calibrated thresholds; excess
  fraction (N_S1 − N_control)/N_S1, spread over draws / noise templates (16–84 %).
- Far/near: flag and far/near ratio at Table-S1 times vs random times (R1 refs and rule).
- LFE coincidence ±60 s (SVI P ≥ 0.5 + WA-NW) at Table-S1 times vs random times, per template and period.
- Seasonality: daily and monthly counts per template vs a microseism/storm proxy = daily median RMS of 10-min windows
  at CN.PGC and UW.LRIV (Z) in 5–10 s (from the 1-sps day files) and 20–50 s (same files, our band); Spearman ρ on
  inter-ETS days, Table S1 vs our catalogue vs noise-template detections (which should follow the proxy if noise-driven).
- Clustering (inter-ETS): Fano factor of daily counts, fraction of gaps < 1 h vs Poisson, Spearman ρ with PNSN tremor
  ≤ 100 km — Table S1 vs our scrambled controls.
Reported for whichever setting is chosen; raw figures first, no error bars.

## Amendment R2b — debug the low overlap with Table S1 (2026-10-15 ~11:40 am PT, before any R2b number)
Seen so far (R2a, logged): at the chosen setting V4 (velocity) only ~10–12 % of Table-S1 detections in the calibration
2011 ETS are within ±30 s of ours, and every R2a variant gave 1–10 %. User rule: treat this as OUR bug until shown
otherwise. The R2a stages continue (noise templates, stage 2, Table-S1 tests); R2b runs alongside, ≤ 8 workers in total.

**R2b.1 Self-detection check** (`hg_r2b.py self`). For V0 (displacement) and V4 (velocity) templates, summed CC over the
template's day ± 1 day: value at the Table-1 origin (must equal the channel count), whether it is the highest peak of
the day, the second-highest peak of the day; our CC (best lag within ±15 s) at the Table-S1 detection nearest each
template origin (found 9 s before for T1, T4, T5, 12 s for T2; none for T3); each template's CC at the other four
template origins. A failure (self not the day's highest, or the Table-S1 self-time CC low) = fix templates / windows /
channel sets before anything else.

**R2b.2 Stack Table S1's detections, independent of our matched filter** (`hg_r2b.py stack`). Per template, all live
channels of its configuration, 0.02–0.05 Hz displacement and velocity (columns), windows [t + Δ, t + Δ + 90] at the
Table-S1 times t, Δ ∈ {0, 9, 11} s (11 = median ours − theirs offset). Per detection: mean per-channel CC with our template
at Δ and at the best lag within ±15 s. Per group: stack CC (mean over channels of CC(channel stack, template channel),
traces unit-normalised before stacking, fixed Δ alignment, no per-trace re-alignment) and semblance (mean over channels).
Groups: 2011 ETS, inter-ETS, 2012 ETS; inter-ETS by season (Sep–Nov 2011, Dec–Feb, Mar–May, Jun–Aug 2012); terciles of
our CC at their times. Controls, same counts: (a) random times on the same UTC days (5 draws); (b) station-day shift
(each station's window from the same time on another day, 5 draws); (c) polarity-scrambled template, real windows
(5 draws); (d) **selection-matched**: the top-N detections (N = their count in that group's period) of 10 noise templates
(chosen-setting runs, each under the same template configuration), stacked and compared with their own noise window —
because R1 showed that stacks of windows selected by correlation resemble the selecting window even when all are noise.
Reading: their windows correlate with OUR template above (a)–(c) → their template ≈ ours (no basic mismatch); stack
CC / per-trace CC at their times far above (d) → real repeats; ≈ (d) → consistent with selected noise.

**R2b.3 Overlap diagnostics** (`hg_r2b.py overlap`, `hg_r2b.py grid`).
- Our calibrated catalogue (V4, Table-S1 thresholds) ranked by summed CC: fraction of our top 10 %, 25 %, 50 %, all
  that are in Table S1 (±30 s; ±5 s around the +11 s offset), per period; same in reverse using our CC at their times.
- Tolerance curve: overlap for windows ±2, 5, 10, 20, 30, 60 s around +11 s and around 0, with the chance level
  (rate × window); histogram of ours − theirs over ±600 s (a time-convention error would show as a peak away from 0).
- Self-consistency: pairwise overlap (±30 s, 2011 ETS) of our own calibrated catalogues V0, V4, V5b, V6b, V3b. If our
  own variants agree with each other far better than with Table S1 → template/processing mismatch; if they also disagree
  with each other → near-threshold detections are not robust (noise-dominated).
- Short processing grid on the calibration period only (1 Jul–6 Sep 2011; thresholds fitted to the Table-S1 2011-ETS
  counts; score = recall of Table-S1 within ±30 s, then ±5 s around the offset), single changes from V4:
  G0 V4 as is (code check: must reproduce V4's 0.101); G1 window [origin − 11, origin + 79]; G2 Z only; G3 horizontals
  only; G4 5 nearest stations; G5 network normalisation; G6 causal (one-pass) 4-pole Butterworth; G7 zero-phase 2-pole;
  G8 corners 0.01–0.05 Hz; G9 declustering 20 s; G10 displacement (code check vs V0's 0.085). Report the best and
  whether the strongest of our detections overlap theirs.
Verdict rule (set now): "template/processing mismatch" if any grid point or the self-consistency test shows a stable,
reproducible catalogue that Table S1 does not share AND their windows correlate poorly with our template; "noise-dominated
near-threshold detections" if our own variants also disagree, the overlap rises strongly toward our top-ranked
detections, and their stacks/per-trace CC sit at the selection-matched control.

## Amendment R2c — moment budget of Table S1 detections (2026-10-15 ~11:48 am PT, before computing)
Coordinator request. Per Table-S1 detection (template K, time t), periods 2011 ETS, inter-ETS, 2012 ETS:
- **Template scaling (primary).** a = least-squares scale of the template's own 0.02–0.05 Hz displacement window (raw,
  not unit-normalised, all live channels of the template configuration, 90 s from the Table-1 origin) onto the data at
  t + Δ_K, with Δ_K the template's own Table-S1 self-offset (T1, T4, T5: 9 s; T2: 12 s; T3: 9 s, no self-detection):
  a = Σ_c w_c·x_c / Σ_c |w_c|², fixed lag (no lag search). M0 = a × M0_K (Table 1 Mw). Column: best lag within ±15 s
  (biased upward by the search; upper reference only).
- **Thrust-synthetic cross-check.** Same LS with our whole-space thrust Green's function (strike 315°, dip 12°, rake 90°,
  Table-1 location and depth, τ = 10 s as in R1, zero-phase 0.02–0.05 Hz like the data) → M0 directly. Applied first to
  the five template events (synthetic M0 vs Table-1 M0 = the scale factor of this tool), then to the detections.
- **Noise contribution, same counts per template and period:** (i) random times on the same UTC days (20 draws) —
  unselected noise, expected ≈ 0 net; (ii) selection-matched: top-N detections of 10 noise templates (N = Table-S1
  count), each scaled onto its own noise window, converted to template-equivalent moment by |noise window| / |template|
  — what selecting noise windows by correlation yields as apparent moment. Net VLFE moment = Table-S1 sum − (ii) mean;
  spread over noise templates / draws (16–84 %) is the uncertainty; Poisson-like scatter of the sum by bootstrap over
  detections (column).
- **Comparison (three-way split per ETS):** LFE moment = catalogued SVI LFEs × 4.5e11 N m (Bostock median magnitude
  1.7; S13 note); LFE-locked slow moment = S13 values (2011-08: 1,445 × 5.3e12; 2012-08: 2,864 × 4.5e12 refit /
  5.1e12 pilot, `../tables/scale_svi_episodes.csv`); discrete VLFE moment = net Table-S1 sum above. GNSS moment: S1
  tables start in 2014 (no 2011/2012 row), so the commonly cited Mw 6.5–6.8 for northern-Cascadia ETS (repo README) is
  used as a range, flagged as not our measurement. Caveats reported: template-scaling bias (noise projects positive
  amplitude on selected windows), noise share, 20–50 s band limits, whole-space ×2 scale, different areas (H&G templates
  vs our SVI LFE families).

## Amendment R3 — individual between-ETS candidates (2026-10-15 ~12:28 pm PT, before any R3 number)
User: "don't you have at least a few real-looking VLFEs in between?" Hints from R2b: their top-10 % between-ETS
detections match our independent run 52 %; T1/T2 between-ETS windows stack at CC 0.45–0.52 with our templates (controls
~0); T3 excess +0.52/day. Caveat: stacks of correlation-selected windows reach 0.6–0.79 even for noise (R1, R2b), so the
aggregate test must be selection-matched and use data NOT used in the selection. Merlin not consulted (coordinator's spec
is explicit) — flagged. Scripts `scripts/hg_r3.py`, queue `scripts/r3_queue.sh`; outputs `data/r3/`, `tables/r3_*`,
`figures/r3_*`; marker `R3_DONE`. ≤ 6 workers, OMP_NUM_THREADS=1, < 60 GB, existing data only.

**R3.1 Candidate set.** Table S1 detections in 2011-09-08 00:00 – 2012-08-28 00:00 UTC (H&G between-ETS window 6 Sep 2011 –
30 Aug 2012 with ±2 d from their ETS windows removed). Matched = our calibrated V4 catalogue (`tables/r2_catalogue_chosen.csv`,
same template) has a detection with |(ours − theirs) − off_K| ≤ 30 s, off_K = 9 s (T1, T3, T4, T5), 12 s (T2) (R2b
self-detection offsets); the raw ±30 s flag is a column. Rank = our CC of the matched detection in MAD units,
(S − median)/MAD (V4 S1 rows of `tables/r2_thresholds.csv`). Sets: ALL matched, and TOP30 by that rank. Origin used by
every test t0 = our matched detection's (origin-aligned) index. Column: inside an SVI minor episode ±2 d (`b_ets_record`).

**R3.2 Per-event tests** (each a column; pass rules fixed here). Scan machinery = `../scripts/scan_columns2.py` `cols()`
(SCAN_YEAR 2011 / 2012, run `real` variant 0, near network PGC VGZ SNB NLLB PFB YOUB LRIV, bp arrays `bp_{year}g`, references
`bp2_ref`). Best scan point = the second of max `zmax` within t0 ± 30 s and its `argmax` node (m0 = that node's m0_cols).
- (a) far/near ratio. ratio6 = median over the live references among {UW.GNW, IU.COR, US.NEW, US.HAWA, UW.OMAK, TA.I02D}
  of A_r / A_near (S7.3 definitions, the 2011 R1 reference set). **Pass a: n6 ≥ 2 and ratio6 ≤ 0.551** (R1 threshold);
  n6 < 2 → untestable (not a pass). Columns: ref_ratio over all live veto references (the 2012 R1 statistic), far-only
  ratio (5 far LHZ), ratio_GNW (GNW is ~50 km from T1, location-confounded), ref_excess, n_ref.
- (b) apparent velocity / plane-wave fit. S7.3 fit re-implemented to also return station distances: pivot = loudest
  live near-network Z in [t0 − 60, t0 + 120]; all near + reference Z, lag ±150 s, stations with CC ≥ 0.6 kept; plane fit
  of lags → VR, slowness, back-azimuth. Travelling = pw_n ≥ 8 AND pw_vr ≥ 0.5 AND ≥ 2 kept stations ≥ 250 km from the best
  node. **Pass b: not travelling** (pw_n < 4 = no coherent wavefield = pass). Column: the same flag at 1,000 random
  between-ETS times (its base rate), apparent velocity 1/slowness.
- (c) location. **Pass c: z_best ≥ 5 AND best node not on the grid edge AND (inside the scan LFE footprint polygon OR
  ≤ 20 km from the nearest LFE family)**; z_best < 5 → "not locatable". Families: `data/catalogue/families_v1_3.csv`,
  regions S_Vancouver_Is, Royer_NW_WA, Royer_SW_WA, Royer_N_Vancouver_Is, Olympic_2010, Olympic_2010_armB (lat/lon, else
  pos_/approx_ coordinates). Columns: best node lat/lon/depth, distance to the H&G template epicentre, distance to the
  nearest family, vr_thrust, vr_free, vr_ratio, and z at the fixed node nearest the template epicentre (max over t0 ± 10 s).
- (d) 2–8 Hz and catalogues (reported, not a pass rule). HF ratio = max 2–8 Hz envelope in t0 ± 60 s and ± 10 min over
  that day's median, at the two scan stations nearest the template epicentre (`hfenv`). Counts within ±60 s and ±10 min of
  t0: SVI LFEs P ≥ 0.75 (svi_rows, data, counted once); WA-NW selected_v2 detections; SVI candidates 0.1 ≤ P < 0.75; PNSN
  tremor epicentres ≤ 100 km from the template epicentre. Each vs 500 random times on the same UTC day (random mean, and
  p = fraction of random ≥ observed). Status: "with LFEs" (≥ 1 SVI P ≥ 0.75 or WA-NW within ±10 min); else
  "candidates/tremor only" (≥ 1 P 0.1–0.75 or PNSN tremor within ±10 min); else "quiet".
- (e) earthquakes. Mask tiers of `masks2_{2011,2012}` at any second of [t0 − 60, t0 + 120]: head2 (teleseismic surface
  waves from the ComCat global M ≥ 4.5 list, Ms = mb − 0.3), bw5 (body-wave windows, M ≥ 5), reg (regional M ≥ 2.5), reg3
  (ComCat M ≥ 3 within 2,500 km). **Pass e: none of them.** m45 (pilot-1 windows, ~69 % base rate) column only; the masking
  event is named.
- (f) single-station dominance. From cols() at the best point. **Pass f: z_loso ≥ 3 AND z_loso ≥ 0.5 × z_best.** Columns:
  max_share, share_excess, loso_dropped; H&G network: the largest one-station share of our V4 summed CC at t0 and the summed
  CC without that station.
- (g) column only: held-out stations (in the scan, not in H&G's network): CN.PGC, CN.SNB, CN.PFB, CN.YOUB (bp displacement
  0.02–0.05 Hz). Mean CC over live channels between [t0, t0 + 90] and the template event's [io, io + 90] at the same channel;
  vs 500 random same-day times (percentile).
- Mw: template scaling (R2c recipe, displacement, K's live channels, fixed lag at t0) × M0_K; and the scan m0_cols at the
  best point.

**R3.3 Aggregate test against selection bias.** For TOP30 and ALL:
- G = mean of (g) over candidates and held-out channels (data NOT used in any selection); G_stack = per-template stack of
  unit-normalised held-out windows, CC with the template's held-out window, mean over channels and templates (weighted by
  count). H = H&G-network stack CC with the template (R2b.2 metric, velocity, at t0; selection-dependent, reported).
- Controls: (i) random same-day times (20 draws). (ii) selection-matched: each of the 100 chosen-setting noise templates
  (`data/r2/noise_vel_0_sum_all_90_*.npz`, configuration of template K) — its declustered peaks in the same window above
  template K's calibrated threshold, ranked; each candidate (template K, rank r of n_K in our V4 between-ETS catalogue for K)
  is replaced by the noise template's peak at rank round(r × n_noise / n_K) (same count per template, same CC-rank
  fraction); G, G_stack and H computed against the noise template's own window (held-out channels from the scan bp arrays
  at the noise-window time).
- Reading (fixed): candidate G above the 97.5th percentile of the 100 selection-matched values → the held-out stations
  confirm template-like events in the set; inside the range → the set is consistent with correlation-selected noise.
  H cannot discriminate (R1/R2b) and is reported only.

**R3.4 Classification.** "Real-looking local VLFE" = a ∧ b ∧ c ∧ e ∧ f. For each: time (UTC, PT), Mw, location, (d)
status, (g) percentile, and a figure of its waveforms at near and far stations. Other classes, first match wins:
catalogued earthquake (e) → travelling wave (a or b) → single-station (f) → not locatable (z_best < 5) → located outside
the footprint/family zone or on the grid edge. Same classification applied, as context for the count, to (i) the candidate
times moved to random times on the same UTC day (5 draws = scrambled control) and (ii) the selection-matched detections of
the first 10 noise templates; the rate of "real-looking" among them is reported next to the candidates' rate.
Raw figures first (points coloured by z or CC, no error bars).
