# S11 round 6: the WA deep-pulse recovery, and does the per-episode deep change scale with episode size?

Written Mon 2026-10-12, ~11:30 pm PT. The rules are in `../PLAN.md` §10, written ~11:04 pm PT before any round-6 number.
- The main queue ran 11:09–11:13 pm PT and every stage had rc = 0.
- §10 Amendment 1 (~11:20 pm PT) is **post-result for task 1**. It was written after the task-1 numbers were seen and before
  any task-2 number was looked at. Its queue ran 11:16–11:17 pm PT.
- One **post-hoc descriptive** column (~11:23 pm PT, decision 31) runs task 2 with the amendment form E4lu. It cannot change
  the pre-registered verdict.

The decision log is `../NOTES.md` (decisions 29–31 and the result lines). No new waveform measurement was made: the round-5
cells, the round-4 waveform injections and the round-5 model class were reused.

**Conventions**
- dv/v is in %.
- ΔD = the deep-weighted change of one slow-slip episode, minus its neighbours within ±3 yr.
- L-SEL = two3_20 (headline); NOM = ell40.
- "Scrambled control": for task 1, each station's data shifted to another time of day; for task 2, the sizes shuffled among
  the episodes.
- Everything here is provisional.

**Checks of the machinery**
- Round-5 E3 and E4 are reproduced exactly (max |Δ| 1e-16 in σ/major, composite, T1 and prior).
- The S3 headline LFE re-count reproduces `n_headline` exactly in both regions.

## In plain words

**1. The cause of WA's low deep-pulse recovery is now known, and one form fixes it under the headline kernel.**
- **There are two losses, and they multiply.**
  - (a) The constraints that keep the slow part out of the episode term also remove part of a one-episode pulse. This loss
    is largest for episodes at the ends of the record or isolated in time:
    - 2008-04 recovers 0.41 (leverage 0.75);
    - 2025-09 recovers 0.38 (leverage 0.73);
    - 2023-03 recovers 0.28.
  - Across the WA majors, recovery falls with leverage (Spearman ρ −0.81). The model-consistent recovery (0.876) matches the
    projection-only prediction (0.861).
  - (b) A waveform-level loss is present even with no constraint: E3 recovers 0.885 in WA and 0.961 in SVI.
  - Together: 0.876 × 0.885 ≈ 0.775. The measured value is 0.748.
- **What did not work.** Constraining D only against the yearly pair (E4y) or against lin + yearly (E4l) lifts WA T1 to
  0.82 / 0.81, but the prior test then fails: 0.15σ / 0.14σ, in the unit terms. The episode term starts trading slow
  content with the per-station slow terms.
- **The fix (Amendment 1, E4lu).** Remove the shared slow directions on both sides:
  - per-station terms: seasonal + linear only;
  - episode term: orthogonal to lin + yearly only.
- **What E4lu achieves.**
  - WA L-SEL: T1 0.817, prior 0.021σ, every gate passes.
  - SVI: everything passes under both kernels.
  - Precision is unchanged or better: σ per major is 0.054 % (SVI, 15/15 ≤ 0.1) and 0.067 % (WA, 23/24).
- **The one remaining failure.** E4lu fails WA NOM on the family split. The two family halves agree on 14 of 16 majors
  (0.875; the rule is 0.9), and their composites agree.
- So E4lu meets the round-5 task-1 pass rule (L-SEL, both regions) but not the stricter round-6 rule (all four
  combinations). It was not adopted. Following the rule fixed in advance, the task-2 headline stays round-5 E4, with sizes
  divided by the measured recovery.

**2. No scaling of the deep-weighted change with episode size is seen. The bounds are now physical numbers.**

**SVI** (GNSS Tier A offset, 15 episodes, 2014+):
- The slope is +0.011 ± 0.009 % per mm (1.2σ; scrambled control p 0.18).
- 2σ bound, after the 0.89 recovery: **≤ 0.033 % per mm of GNSS offset**. For a ~5 mm ETS, that is ≤ ~0.17 %.
- To be seen at 3σ, a large episode (q90) would have needed to differ from a small one (q10) by **0.17 %**.

**WA** (PNSN tremor count, 23 episodes):
- The slope is +0.022 ± 0.013 % per 1,000 epicentres (1.7σ; p 0.095).
- 2σ bound, after the 0.69 recovery: **≤ 0.069 % per 1,000 epicentres**.
- The 3σ-visible large-minus-small difference is **0.16 %**.

**All 15–16 other measure tests.** None reaches the SCALES rule. These cover moment, tremor, LFE count, duration, patch
coverage and family-set activity, in each region and for every ΔD source.

**3. One WA thread is worth a note: longer or larger WA episodes read slightly higher during the episode. This is a hint,
not a finding.**
- Top third by tremor vs the rest: +0.065 ± 0.025 % (scrambled p 0.035).
- Top third by duration: +0.072 ± 0.027 % (p 0.014).
- Duration slope: +0.033 ± 0.015 % per 10 d (p 0.050). With E4lu it is 0.049 ± 0.015 (3.2σ, p 0.008).
- Removing 2023-03 or the two WECC windows does not remove it (2.1–2.2σ).
- **Why it stays a hint:**
  - about 24 tests are run per region and ΔD source;
  - the before → after step for the same top third goes the other way (−0.025, p 0.022);
  - SVI shows nothing comparable (top third by GNSS +0.043 ± 0.031, p 0.25).
- If it is real, it is a positive (stiffening) excess of ~0.07 % during the longest WA episodes. It does not persist past the
  episode.

**Overall reading.** The per-episode deep-weighted change is flat against slip size at ~0.01–0.02 % per mm (1σ). A big ETS
changes the deep-weighted velocity by at most ~0.17 % more than a small one (SVI) and ~0.16 % (WA, per tremor size).
- This is a new and coarse, but clean, bound.
- It turns "no average change" (round 5) into "no change that grows with slip".

## Task 1: the WA deep-pulse recovery (PLAN §10.1 + Amendment 1)

**Cause analysis** (`{svi,wa_nw}/tables/r6_cause.csv`, `r6_t1_permajor.csv`, `fig_r6_t1.png`).

Medians over majors, L-SEL:

| form (D ⊥ …) | WA recovery (waveform) | WA model-consistent | projection-only | WA h_ee median / max | ρ(h_ee, recovery) WA | SVI recovery |
|---|---|---|---|---|---|---|
| E3 (sum-to-zero only) | 0.885 | 1.000 | 1.000 | 0.037 / 0.037 | −0.11 | 0.961 |
| E4 (lin + 5 hinges + yearly) | **0.748** | 0.876 | 0.861 | 0.256 / 0.750 | −0.81 | 0.901 |
| E4h (no end hinges) | 0.770 | 0.884 | 0.877 | 0.226 / 0.750 | −0.56 | 0.915 |
| E4l (lin + yearly) | 0.813 | 0.931 | 0.924 | 0.141 / 0.249 | −0.33 | 0.930 |
| E4y (yearly only) | 0.820 | 0.947 | 0.932 | 0.115 / 0.160 | −0.25 | 0.936 |
| E4lu (lin + yearly; units seasonal + lin) | 0.817 | 0.933 | 0.924 | 0.141 / 0.249 | −0.35 | 0.916 |

**Full suite.** `r6_forms.csv`, `r6a_forms.csv`; rule ≤ 0.35 for leaks, < 0.1σ for the prior.

| form | combo | σ/major (≤ 0.1) | composite | T1 | T2 leak | T4 all comp / major | prior | fails |
|---|---|---|---|---|---|---|---|---|
| E4 (r5) | WA L-SEL | 0.065 (23/24) | −0.004 ± 0.009 | **0.748** | 0.042 | 0.043 / 0.119 | 0.034 | T1 |
| E4y | WA L-SEL | 0.075 (21/24) | −0.022 ± 0.013 | 0.820 | 0.038 | 0.023 / 0.166 | **0.150** | G-prior |
| E4y | WA NOM | 0.098 (13/24) | −0.031 ± 0.013 | 0.882 | 0.322 | 0.003 / 0.152 | 0.072 | G-ctrl, G-fam |
| E4l | WA L-SEL | 0.072 (21/24) | −0.023 ± 0.013 | 0.813 | 0.040 | 0.018 / 0.161 | **0.142** | G-prior |
| E4h | WA L-SEL | 0.068 (23/24) | −0.004 ± 0.009 | **0.770** | 0.048 | 0.047 / 0.149 | 0.029 | T1 |
| E4y+Scs2 | WA L-SEL | 0.076 (21/24) | −0.006 ± 0.013 | 0.820 | 0.029 | 0.027 / 0.159 | **0.187** | G-prior (also SVI NOM prior, WA NOM ctrl) |
| **E4lu** | SVI NOM | 0.066 (15/15) | −0.013 ± 0.028 | 0.943 | 0.320 | 0.018 / 0.175 | 0.028 | — |
| **E4lu** | SVI L-SEL | **0.054 (15/15)** | −0.002 ± 0.024 | 0.916 | 0.056 | 0.030 / 0.135 | 0.035 | — |
| **E4lu** | WA NOM | 0.086 (22/24) | −0.023 ± 0.010 | 0.874 | 0.321 | 0.020 / 0.174 | 0.023 | G-fam (0.875) |
| **E4lu** | WA L-SEL | **0.067 (23/24)** | −0.014 ± 0.011 | **0.817** | 0.067 | 0.041 / 0.103 | 0.021 | — |
| E4yu | SVI L-SEL / WA both | 0.053 / 0.067 | | 0.908 / 0.823 | | | ≤ 0.038 | SVI G-ctrl, WA G-fam (0.75) |

**E4lu, other gates.**
- WA L-SEL:
  - T3 ≤ 0.11 and T5 0.013;
  - B2 scrambled control 100 % of majors within 2σ, composite 0.000 ± 0.012;
  - leave-one-site-out 1.0; family split 0.94;
  - G-site: B013+B943 carry 31 %.
- SVI: G-site flag (PGC 65 %) as before; B1 column 79 % (as round 5).
- No major is at ≥ 3σ under any kernel. The largest are:
  - WA: 2017-02 −0.14 ± 0.06 (−2.2σ), 2008-04 +0.12 ± 0.07 (1.8σ);
  - SVI: 2019-11 +0.11 ± 0.08 (1.4σ).
- WA 2023-03 falls to +0.11 ± 0.08 (1.3σ). In round 5 it was 2.6σ under L-SEL, so it is not a candidate in E4lu.

**Q3 kernel band (E4lu).**
- The composite is within 2σ of 0 for every plausible kernel.
- 2σ bound: SVI 0.077 %, WA 0.052 %.
- T1 over plausible generators: SVI 0.67–0.92, WA 0.60–0.82.
- Size bound per average episode: **SVI 0.08–0.12 %, WA 0.06–0.09 %.** This is the same as round 5 (0.08–0.11 / 0.05–0.07),
  so the bound did not move.

**Pass rule (round 6, all four combinations): not met by any form.** E4lu is closest, failing on one family-split item in
the NOM column. **Adoption: none.**
- Under the rule fixed in advance, the task-2 headline is round-5 E4.
- E4lu runs as a post-hoc descriptive column (`r6p_*`). It does not change any verdict below.

## Task 2: ΔD against episode size (PLAN §10.2)

**Sizes** (`{svi,wa_nw}/tables/r6_sizes.csv`):

| | SVI (53 model episodes) | WA-NW (27) |
|---|---|---|
| headline | G = S1 Tier A GNSS offset (15) | T = PNSN tremor within 10 km of a family (23) |
| other measures | M0 Tier B moment (15), T (20), L = S3 LFEs (15), Du duration (20), Pc S3 patch coverage (15, 12 groups), Fs / Fc family-set sumP / coverage (53) | G_WA = S1 WA-station offset of the overlapping SVI window (15), L = NW LFEs once P ≥ 0.75 (25), Du (25), Pc (17, 9 bins of 10 km), Fs / Fc (27) |

**Headline regressions** (WLS 1/σ², round-5 E4, L-SEL; `r6_size_regress.csv`, `r6_size_verdict.json`):

| region, measure | n | slope b ± σ | z | scrambled p (99 % of shuffled) | recovery | 2σ bound (÷ recovery) | 3σ-visible change q90 − q10 | verdict |
|---|---|---|---|---|---|---|---|---|
| SVI, GNSS mm | 15 | +0.011 ± 0.009 % / mm | 1.2 | 0.18 (0.020) | 0.89 | 0.033 % / mm | 0.17 % | **NO SCALING SEEN** |
| WA, tremor / 1,000 | 23 | +0.022 ± 0.013 % / 1,000 | 1.7 | 0.095 (0.033) | 0.69 | 0.069 % / 1,000 | 0.16 % | **NO SCALING SEEN** |

**Columns.** Per SD of size, slopes are 0.006–0.03 % in every case.

- SVI:

  | measure | slope | z | scrambled p |
  |---|---|---|---|
  | moment | +0.009 ± 0.010 % per 10^18 N·m | 1.0 | 0.28 |
  | tremor | +0.009 ± 0.009 | 1.0 | 0.17 |
  | LFE count | +0.014 ± 0.045 | 0.3 | |
  | duration | −0.009 ± 0.023 | | |
  | patch coverage | +0.03 ± 0.07 | | |
  | family-set sumP | +0.004 ± 0.004 | | |
  | family-set coverage | +0.003 ± 0.033 | | |

  - NOM and E4+Scs2 give the same picture.
- WA:

  | measure | slope | z | scrambled p |
  |---|---|---|---|
  | GNSS (WA stations) | 0.000 ± 0.011 % per mm | | |
  | LFE count | +0.013 ± 0.027 | | |
  | **duration** | **+0.033 ± 0.015 % per 10 d** | **2.3** | **0.050** |
  | patch coverage | +0.05 ± 0.08 | | |
  | family-set sumP | −0.02 ± 0.08 | | |
  | family-set coverage | +0.05 ± 0.05 | | |

  - The duration slope is +0.032 ± 0.015 (p 0.039) in E4+Scs2. It is +0.049 ± 0.015 (z 3.2, p 0.008) in the E4lu
    post-hoc column, but only after its 0.61 recovery is considered.
- GLS columns agree in sign or are smaller. For WA G_WA and Pc the GLS σ is implausibly small (0.002–0.003), because the
  ΔD covariance is near-singular there. Those GLS numbers are not used.

**Top third vs the rest** (`r6_size_topthird.csv`; the headline measure and the WA hint):

| region, measure | during-episode contrast | scrambled p | before → after step contrast | scrambled p |
|---|---|---|---|---|
| SVI, GNSS | +0.043 ± 0.031 | 0.25 | +0.016 ± 0.016 | 0.57 |
| SVI, tremor | +0.030 ± 0.029 | 0.23 | +0.021 ± 0.014 | 0.35 |
| WA, tremor | **+0.065 ± 0.025** | **0.035** | −0.001 ± 0.013 | 0.94 |
| WA, duration | **+0.072 ± 0.027** | **0.014** | −0.025 ± 0.025 | 0.022 |
| WA, LFE count | +0.029 ± 0.024 | 0.34 | −0.001 ± 0.010 | 0.94 |

- The largest SVI and WA episodes show no lasting step: every step contrast is ≤ 0.02 % with |z| ≤ 1.5.
- In E4lu, SVI top-third-by-GNSS step is +0.07 ± 0.05 (p 0.034), with 1.6σ formal. It is a column only and no other
  measure agrees.

**Verdict.**
- **SVI: NO SCALING SEEN** (GNSS headline).
- **WA: NO SCALING SEEN** (tremor headline), with a hint, at the multiple-test level, that the longest/largest WA episodes
  read ~0.07 % higher during the episode.

## What failed, the measured obstacle, and the next adaptation

1. **The round-6 pass rule (all four combinations) is not met.** E4lu fails one WA NOM family split (14/16 majors).
   - Obstacle: WA NOM ΔD noise (σ 0.086) against the halved family set.
   - Next options:
     - adopt E4lu by the round-5 rule (L-SEL both regions), with NOM as a column;
     - or test E4lu on a 3-way family split.
2. **The waveform-level loss in WA (E3 0.885) is untouched.** It sits mostly in 2023-03 (0.56), 2019-03 (0.64) and 2015-07
   (0.70), even without constraints.
   - Next: per-episode check of those cells (cell count, stretch-grid limits, reference path). Divide by per-episode
     recovery instead of the median.
3. **The size bound is limited by the number of geodetically sized episodes** (15 in SVI).
   - Next: extend the S1 Tier A GNSS sizes to 2010–2013 (4 more SVI majors with cells). Use a per-WA-episode WA GNSS series
     instead of the SVI-window proxy.
4. **The WA duration/tremor hint.**
   - Next, cheapest first:
     - split the WA top third by site (are they PGC-site or borehole cells?);
     - check whether duration tracks the number of stacked days per cell. Longer episodes give more days per stack, and a
       stack-length bias would look like this.

## Files
- **Scripts:**
  - `svi_science/scripts/s11_r6_{model,size,figures}.py`;
  - queues `s11_r6_queue.sh`, `s11_r6a_queue.sh`;
  - resume blocks 11 and 11b in `scripts/resume_after_reboot.sh`.
- **Code tests (not results):** `{svi,wa_nw}/tables/test/` (forms: 30 WA series; size: sizes shuffled first).
- **Task 1:**
  - `round6/adopted.json`, `adopted_amend.json`;
  - `{svi,wa_nw}/tables/r6_forms.csv`, `r6a_forms.csv` (all gates and rules), `r6*_forms_rules.csv`, `r6*_forms_onesite.csv`,
    `r6*_forms_dD.csv`, `r6*_t1_permajor.csv`, `r6*_cause.csv`, `r6*_band_q3.csv`, `forms_reproduction_check.csv`,
    `forms_<form>_<setup>_{dD,loso,ctrl_*}.csv`.
- **Task 2:** `r6_sizes.csv`, `r6_size_{regress,topthird,dD}.csv`, `r6_size_verdict.json`; post-hoc E4lu column
  `r6p_size_*.csv`.
- **Figures:**
  - `figures/fig_r6_t1.png`: per-major recovery vs leverage;
  - `figures/fig_r6_size.png`: ΔD vs size, headline E4, coloured by weight, no error bars;
  - `figures/fig_r6_size_E4lu.png`: the post-hoc column.
