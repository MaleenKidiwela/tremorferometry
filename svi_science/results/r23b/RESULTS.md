# R23b — findings 1 and 6 re-done on our own 23-year SVI LFE record (2003–2026)

Written Sun 2026-10-11, ~6:35–7:00 am PT. Rules in `PLAN.md`, written before any number was computed. Revisions 1–2 at
its end are post-hoc and labelled as such below. Scripts: `svi_science/scripts/r23b_{base,f1,f1_posthoc,f6_rates,f6_wave,f6_size,figures}.py`.
Tables are in `tables/` and figures in `figures/`. Science only: no catalogue or rate file was changed.
"Scrambled control" = the same detector on the same day, with each station site's data shifted to another time of that day.

**Data used (the RESULTS_BACKWARD usage rules):**
- Single repeats: primary detector, counted once; era rows at P ≥ 0.90, frozen rows at P ≥ 0.75. Forward 2014–26 at
  P ≥ 0.75, counted once. These are the r23 event rows.
- Episode and patch sizes across networks use only the per-family Bostock-equivalent `rate`.
- No between-episode rates are used anywhere.
- Episodes are the r23 list: 145 episodes, 52 majors, 23 big ETS (10 in 2003–13, 13 in 2014–26).

**Foundation checks:**
- The rebuilt counted-once counts equal Q1's exactly (7,753 events, 0 cells differ).
- Tremor extents are within 0.06 km of Q1's.
- Before 2010 the envelope tremor extents (fixed CN network) agree with PNSN on the four 2010–13 big ETS: median
  difference 5.4 km, same direction 4/4. They count as validated.

---

## Finding 1 — "Each ETS slips the fault in a different pattern, set mostly by how far the front gets"

**Verdict: undecided.** Under the letter of the pre-registered rule the verdict is "changes". Revisions 1–2 show that both
triggers are artefacts of the design, not the fault:
- K1 has no leverage in 2003–13.
- The joint "lasting patches" signal does not survive like-for-like checks.

In plain words:
- **Each ETS slipping in its own pattern holds over all 23 years.**
  - In 2003–13 the patch shares differ between big ETS far beyond counting noise: χ²/df 38, against 28 in 2014–26.
  - Neighbouring families within ~10 km change together in both decades, as strongly as before.
  - There are still no patches that stay strong or weak.
- **"Set mostly by how far the front gets" is only testable in 2019–2022.**
  - Every one of the 10 big ETS of 2003–13 had tremor over the whole ~90-km box. So did the four 2014–18 ETS, 2021-10 and
    the three 2023–25 ETS (92–100 % of the patches inside the tremor extent).
  - Only the five 2019–22 half-island ETS (2019-08, 2019-11, 2020-11, 2021-01, 2022-03) did not reach the whole box.
  - The front-reach explanation (R² 0.53) rests on those five. The other 18 big ETS still differ strongly.
  - Nothing we tested explains those 18. In 2003–13, every predictor reaches within-episode R² ≤ 0.07: direction,
    position along the path, time since last activity, previous residual, station availability.

### Numbers

Headline unit `C_eq` = Bostock-equivalent rate shares put on the counted-once scale (PLAN §1). The 12 S0 patches are
used. p values are against the Q1/S3 controls, which are unchanged.

| test | 2003–13 (10 big ETS, new) | 2014–26 (13) | 2003–26 (23) | before (Q1/S3, 2015–26, 10) |
|---|---|---|---|---|
| inside tremor extent, β (R², p) | +0.08 (0.02, p 0.20) — **no out-of-reach cell except 1** | +0.97 (0.53, p 0.0005) | +0.73 (0.45, p 0.0005) | +1.37 (0.53, p 0.0005) |
| km outside extent | −0.08 (0.02) | −0.85 (0.40) | −0.63 (0.35) | −1.18 (0.39) |
| all front terms R² | 0.06 | 0.59 | 0.51 | 0.62 |
| unevenness vs tremor extent (Spearman) | −0.01 (p 1.0) | −0.88 (p 0.0002) | −0.38 (p 0.08) | −1.00 |
| days since last active, H1 (+) | +0.02 (p 0.53) | +0.32 (p 0.019) | +0.16 (p 0.027) | +0.46 (p 0.073) |
| turn-taking, H2 (previous big-ETS residual, −) | **+0.06 (p 0.67)** | −0.71 (p 0.0005) | −0.48 (p 0.0005) | −1.11 (p 0.001) |
| budget beyond front terms, H4 ΔR² | +0.01 (p 0.66) | +0.06 (p 0.009) | +0.04 (p 0.003) | +0.06 (p 0.07) |
| persistence, rule 7a: mean residual corr. vs patch-label shuffle | 0.135 vs 0.149 (p 0.75) | 0.114 vs 0.092 (p 0.012) | **0.110 vs 0.070 (p 0.002)** | 0.103 vs 0.092 (p 0.094) |
| consecutive-pair corr. (catch-up 7b) | +0.30 (p 0.92) | −0.32 (p 0.0005) | +0.01 (p 0.10) | −0.45 (p 0.001) |
| patch χ²/df | 38 | 28 | 30 | 45 |

- **Reproduction (R0).** 2014–26, S3's 10 big ETS, counted once. Same numbers as Q1/S3 to the third decimal: K1 +1.372 /
  R² 0.529, H2 −1.112 (p 0.001), 7a 0.103 (p 0.094), 7b −0.451, unevenness ρ −1.00. The code is unchanged.
- **Calibration** (100 no-memory multinomial records through the 2003–26 design): the order-shuffle test fires at
  p < 0.01 in 0–4 % of runs. H1 3 %, H2 0 %, H3 0 %.
- **Columns:**
  - counted-once counts with network-segment pools;
  - `rate_kc` shares;
  - ETS with ≥ 300 LFEs only;
  - envelope tremor for 2010–13 as well.

  All give the same picture. In 2003–13, K1 R² ≤ 0.03 and H2 > 0. In 2003–26, K1 R² 0.41–0.49 and H2 −0.48 to −0.58.
- **Do 2003–13 patterns cluster by detector segment?** (post-hoc, Revision 1) Same-segment pairs minus different-segment
  pairs: +0.25 (p 0.046). This is a hint, not passed. In 2003–26: +0.07 (p 0.13).
- **Where the joint "persistence" comes from** (post-hoc, Revision 2, `tables/f1/posthoc_persistence.json`).
  Rule 7a rerun with the same code:

  | big ETS used | n | mean residual corr. vs label shuffle | p |
  |---|---|---|---|
  | all 23 | 23 | 0.110 vs 0.070 | 0.002 |
  | full-box only (tremor reaches ≥ 90 % of patches) | 18 | 0.161 vs 0.095 | 0.001 |
  | without the 3 era-detector ETS | 20 | 0.070 vs 0.053 | 0.004 |
  | **full-box and not era detector** | **15** | **0.085 vs 0.076** | **0.13** |
  | all 23, each ETS weighted equally | 23 | 0.082 vs 0.046 | 0.0005 |

  It needs either the repeated half-island ETS of 2019–22 (three NW-half and two SE-half look alike) or the era-detector
  ETS. The era-detector ETS hold 17.6k of the 40.7k counted-once events, so they dominate the pooled share. Among
  comparable whole-box ETS on frozen detectors, no lasting strong or weak patches are seen.
- **Neighbours move together** (`episode_rates` share test + Mantel, 52 majors). Residual correlation for family pairs at
  < 3 / 3–6 / 6–10 / 10–20 / 20–40 km:

  | block | correlation by distance bin | distance coefficient (p) | over-dispersed families (q < 0.05) | median χ²/df |
  |---|---|---|---|---|
  | 2003–13 | 0.61 / 0.42 / 0.28 / 0.07 / −0.10 | −0.23 (p 0.001) | 94 / 114 | 4.4 |
  | 2014–26 | 0.59 / 0.46 / 0.36 / 0.17 / −0.06 | −0.49 (p 0.001) | 24 / 97 | 1.4 |
  | 2003–26 | 0.59 / 0.39 / 0.26 / 0.08 / −0.07 | −0.29 (p 0.001) | 77 / 114 | 2.6 |

  - On one network only (2010–13 counted once, 7 majors): 0.61 / 0.52 / 0.38 / 0.16 / −0.01, and 25 / 77
    over-dispersed. That is close to 2014–26 (FWD counted once: 29 / 73).
  - So the large 2003–13 over-dispersion count (94 / 114) partly comes from the bigger era-detector counts (more power)
    and from per-set k errors across detector segments.
  - The shared-station term is positive in 2003–13 (+0.21, p 0.001), as Washington was before. It is zero in 2014–26
    (p 0.40).

### What is new vs the earlier version
- The earlier statement, "patch strength is set mostly by how far the front gets (R² ~0.5–0.6)", depends on the five
  half-island ETS of 2019–22. All 14 big ETS before 2019, 2021-10 and the three 2023–25 big ETS crossed the whole box. Among those
  whole-box ETS the patterns still differ strongly, and none of our predictors explains it. **The "why" is open again
  for whole-box ETS.**
- The turn-taking of 2015–26 (H2) is absent in our own 2003–13 record (+0.06, p 0.67), as it was in Bostock's
  catalogue. It is a feature of the 2019–22 alternating half-island fronts, not a persistent law.
- The neighbour coupling (≤ 10 km) is a 23-year property.

**Figures:**
- `figures/r23b_f1_neighbours.png` (analogue of `f1_neighbours.png`): correlation vs distance and vs shared stations,
  per decade.
- `figures/r23b_f1b_why.png` (analogue of `f1b_why.png`): R² per factor per decade, and unevenness vs tremor extent.
  Diamonds are half-box ETS.

---

## Finding 6 — "LFE sizes and rates did not change over the years"

### 6b. Size per LFE (S6 ΔM), 2003–2026, inside ETS — **verdict: holds (no change)**
Within ETS, a family's LFEs have the same amplitude in 2003 as in 2026, to within ~±1 % per decade. One small,
provisional exception, PGC alone, is described below.

**Method:** S6's ΔM is measured on the same four CN stations for the whole record (PGC, SNB, VGZ, NLLB).
- Raw miniSEED is used for 2003–13 and the S6 store for 2014–26.
- The amplitude is the S6 matched filter against the same family × station stacks (stacks equal S6's exactly).
- Station-epoch correction (S6 rule) is applied at every CN instrument change.
- Only events above a per-family completeness floor are compared. The floor is shared by all compared segments.

**Station choice:** the KELB gain problem and the era-only `rel_amp` feature are not used: KELB is not in the station
set. "Same station set" is the comparison form throughout; no k-correction is needed for sizes.

**Reader check:**
- The raw reader reproduces S6's amplitudes on 10,541 forward station-detections: median |Δ log10 a| 1e-8.
- It matches the cut store on 11,215 station-detections of 2011–13 (median 4e-8).
- The S-window check passes in 93–97 % of family × stations in every segment.

| test | effect | 95 % CI | p | noise ΔN effect | verdict |
|---|---|---|---|---|---|
| T6a trend 2003–2026, 69 families, 7,925 LFEs above the shared floor | **+0.0005 log units / decade** | −0.004 … +0.005 | 0.67 | −0.038 | not seen; bound ±~1 % per decade |
| T6a, big ETS only | +0.0007 | −0.004 … +0.005 | 0.66 | | not seen |
| T6a, all rows (no floor) | +0.072 | +0.055 … +0.090 | 0.0005 | | detector selection, see below |
| T6b like-for-like: 2010–13 S15f vs 2014–18 FWD-S15 (63 families) | +0.018 | +0.004 … +0.031 | 0.26 | +0.020 | not seen |
| T6c PGC alone 2003 → 2015-01 (one instrument, no correction) | +0.032 / decade | +0.016 … +0.048 | 0.002 | −0.022 | **seen** by the S6 rule |
| T6c PGC, frozen detectors only 2006-07 → 2015-01 (post-hoc, 59 families) | +0.011 / decade | −0.017 … +0.038 | 0.97 | +0.038 | not seen |
| ERA (2003–06) − FWD above shared floor | +0.020 | +0.009 … +0.031 | 0.27 | −0.054 | — |
| B2/B3 (2006–09) − FWD | +0.003 | −0.013 … +0.018 | 0.83 | +0.012 | — |
| S15f (2009–13) − FWD | −0.004 | −0.016 … +0.009 | 0.83 | −0.035 | — |

- **Why "holds":** the pre-registered rule needs T6a not seen with its CI inside ±0.05, and T6b not seen. Both are met.
- **PGC alone.**
  - It shows a small rise of about 7 % per decade, 2003–2015, that passes the S6 bar.
  - Its all-rows version (+0.24 / decade) shows how strongly detector selection pushes in the same direction: the era
    detector of 2003–06 finds many more small LFEs.
  - So a leftover selection effect is the likely cause. Without the era detector (post-hoc, 2006-07 → 2015-01), PGC alone gives +0.011 per decade (p 0.97).
  - It stays provisional. It is inside the ±0.05 bound in any case.
- **Network caveat:**
  - Without the completeness floor, sizes "grow" by 0.07 log units per decade. That is the era detector finding smaller
    LFEs, not a change in the LFEs. Only comparisons above the shared floor are meaningful.
  - PGC read ~1.7× high relative to the other stations from 2000 to 2015-01-12 (station-epoch term +0.23, as S6 found).
    The correction removes this; it is a constant within that epoch.
  - The noise level ΔN falls over the record (−0.04 per decade). Noise does not drive the ΔM result.

**New vs the earlier version:**
- S6 covered 2014–26 (−0.008 per decade, ±0.03).
- The size constancy now holds over 23 years, on one fixed station set and across four detector eras. The bound is
  tighter: ±~1 % per decade.

### 6a. ETS-episode sizes per family (Bostock-equivalent `rate`) — **verdict: undecided**
Big ETS of 2014–26 contain ~3× fewer Bostock-equivalent LFEs per family than those of 2003–13. That is significant,
but it cannot be separated from network changes.
- On the same detector model and the same stations (2010–13 vs 2014–18), the two ordinary 2014–18 ETS are within ~25 %
  of 2010–13.
- The two very small 2014–18 ETS have known problems.
- The decade contrast is driven by the large 2003–05 ETS and by the small, split 2019–22 ETS.

| comparison (median over families of late / early mean ETS size) | ratio | 95 % CI | p (episode labels permuted) |
|---|---|---|---|
| 2014–26 vs 2003–13, `rate` (111 families; 13 vs 10 big ETS) | **0.34** | 0.20 … 0.57 | 0.001 |
| same, `rate_kc` | 0.20 | 0.09 … 0.48 | 0.001 |
| like-for-like 2014–18 FWD-S15 vs 2010–13 S15f (4 vs 4) | 0.42 | 0.09 … 0.88 | 0.17 (smallest possible p 0.029) |
| same, counted-once counts | 0.41 | 0.05 … 0.93 | 0.17 |
| like-for-like without 2014-11 and 2017-03 (post-hoc) | **0.74** | 0.48 … 1.28 | 0.73 |
| 2014–26 (all) vs 2010–13 S15f | 0.40 | 0.23 … 0.66 | 0.04 |

- **Rule outcome:** the like-for-like CI excludes 1, but the exact permutation can never reach p < 0.01 with 4 vs 4 ETS.
  This design flaw is stated in Revision 2. The decade CI excludes 1 under both `rate` and `rate_kc`. So the rule gives
  **undecided**.
- **Per tremor** (PNSN tremor family-exposure-days; never across the env/PNSN switch). Σ rate per tremor family-day:
  - 2010–13 big ETS: 8.0, 8.5, 14.4, 6.8.
  - 2015-12: 7.1, and 2018-06: 7.2.
  - 2014-11: 1.3 (139 LFEs under 1,620 tremor family-days, mostly at the south end).
  - 2017-03: 1.6 (B010 out all 18 days, VGZ swap: the S0 network flag).
  - 2019–25 on the S14/S13 station sets: 2.0–7.3, median ~3.9.
  - All 2014–26 vs 2010–13: ×0.48 (p 0.006). All majors, 2014–18 S15 vs 2010–13: ×0.29 (p 0.09).
- **Release per year** (all majors): 20,100 Bostock-equivalent LFEs per year in 2003–13 against 8,900 in 2014–26.
  On S15: 17,400 per year (2009-08 → 2013) against 7,600 (2014–18).
- **Bostock's own catalogue** (his method; his stations also changed):
  - His per-family counts fall over his 10 windows (−0.30 log units per decade). Episode-order permutation p 0.39: the
    big 2004/2005 ETS dominate.
  - Our `rate` tracks his per-family counts: Spearman 0.82 (era), 0.74 (B2), 0.62 (B3), 0.85 (S15f).
- **Which families carry the ETS is the same in both decades:** family-share Spearman 0.88 (p 1e-37).
- **Network caveat:**
  - `rate` and `rate_kc` differ ×1.7 for 2003–13 (P7 K brackets: ERA 0.66–1.59, B2 0.41–1.12).
  - The forward network changed station set in 2019 (S14) and 2021 (S13).
  - The k of a set is calibrated on 2010–13 truth days only.
  - A drop of ×2–3 is inside what these network steps could produce. The two clean like-for-like ETS show ×0.74.
  - Between-episode rates were not trended (usage rule).

**New vs the earlier version:**
- Finding 6's rate part was Washington-only before (fixed boreholes, no change 2007–26).
- For SVI, LFE output per big ETS is not demonstrably constant: the raw 23-year contrast is ×0.34.
- A fixed-network test does not show it, but that test has 2–4 ETS per side and cannot settle it.
- An independent size measure (GNSS slip before 2014, or a fixed-station LFE count 2003–26) would be needed.

**Figures:**
- `figures/r23b_f6_size.png` (analogue of `f6_size.png`): yearly median ΔM per group 2003–2026, fixed CN stations,
  above the shared floor, with noise.
- `figures/r23b_f6b_ets_sizes.png`: big-ETS sizes 2003–26 with Bostock's counts, per-tremor sizes 2010–26, and the
  ratio tests.

---

## Files
- `tables/`:
  - `windows.csv`, `fam_ep.parquet`, `patch_ep.parquet`, `tremor_ep.csv`, `det_patch_ep.csv`, `base_summary.json`;
  - `f1/`: `f1_results.json`, `f1_table.csv`, rows / kin per design, `neighbour_corr_bins.csv`, `pairs_*.csv`,
    `families_*.csv`, `posthoc_persistence.json`;
  - `f6/`: `fam_episode_sizes.parquet`, `episode_sizes.csv`, `f6_rates.json`, `f6_rates_posthoc.json`,
    `wave_targets.parquet`, `wave/*.npz` (3-C waveforms, 0.49 GB), `cn_stacks.npz`, `mf_cn_backward.parquet`,
    `reader_validation.json`, `sizes_2003_2026.parquet`, `f6_size.json`, `cn_station_epoch_residuals.csv`,
    `completeness_by_segment.csv`, `cn_window_check_segments.csv`, yearly medians.
- Logs: `logs/`.
