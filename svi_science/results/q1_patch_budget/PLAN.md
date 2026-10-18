# Q1 — Why does each ETS slip the fault in a different pattern? (slip budget vs front kinematics)

Pre-registered 2026-10-07 2:30 pm PT, before any Q1 number was computed. Changes after this point go in §9 with a time.
Scripts `svi_science/scripts/q1_*.py`; outputs here; figures `svi_science/figures/q1_*.png`. CPU only, < 30 GB RAM.

## 1. Question and hypotheses

A patch's share of an episode's LFEs varies between episodes beyond counting noise, and neighbours within ~10–20 km
vary together (episode_rates RESULTS, S0b, S3 rule 7). Why?

- **H-budget (main).** A patch fires more in an ETS if it has gone longer since it last slipped / slipped less before.
  - H1: longer time since the patch was last active → larger residual share (coefficient > 0).
  - H2: larger residual in the previous big episode (≥ 90 days earlier) → smaller residual now (coefficient < 0).
  - H3: larger cumulative residual (running sum over earlier big episodes ≥ 90 d earlier) → smaller now (coefficient < 0).
- **Alternatives.**
  - K1 reach: the patch lies inside the along-strike extent the PNSN tremor covered in that episode (independent catalogue).
  - K2 direction mode: patch's mean residual in OTHER big episodes moving the same way (≥ 90 d apart) predicts it now.
  - K3 position along the front path (0 = entry end, 1 = far end), from the tremor direction.
  - K4 front speed, episode size (LFE count; SVI GNSS Tier A offset `D`, S1), tremor extent width: episode-level, so they
    can only change how UNEVEN an episode is → tested on episode unevenness (§4 E1).
  - D1 detectability: the patch's station availability in that episode relative to its own average.
- **L2 (does the budget steer the front?).** Time since last active predicts whether the tremor front reaches the patch.

## 2. Data (frozen)

- **SVI** `svi_science/results/s0/analysis_table.parquet` (run prod, series data), 2014-01-01 → end. Headline = P_same ≥ 0.75
  AND keep_xfam8 (counted once). Patches = the 12 S0 location groups (`s0/groups.csv`, position = group centroid in the S0
  Slab2 frame). Episodes = `lfe_extension_wa/results/episode_rates/episodes_svi_final.csv` (79, MAJOR + minor), padded by
  3 days each side, capped at half the gap (S3 rule).
  - Variants: raw (no de-dup) P ≥ 0.75; counted-once P ≥ 0.5; **fixed-S13** (`analysis_table_S13.parquet`, same rules);
    patches = 10-km along-strike bins of family positions.
  - Scrambled control: series c1 through the same steps (§5).
- **WA-NW** `lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet`, region NW, selected_v2 (own_pnet ≥ 0.75);
  counted once = dedup variant D == False (`post_v2/dedup/flags`). Patches = 10-km along-strike bins of the array-confirmed
  source positions exactly as S3 (`s3/sources_nw.csv`, family → source by `nw_sharing/source_groups_E3*.csv`).
  Episodes = `episodes_nw_final.csv` (118, all eras), padded as above.
  - Target episodes: era E3 (2015-01-01 → 2026-09-15). History for time-since-last-active: E3 only (primary); E2 + E3
    (2012 →) and E1 + E2 + E3 (2007 →) as variants (activity judged against each era's own long-term shares).
  - Variants: raw P ≥ 0.75; counted-once P ≥ 0.5; bins shifted by 5 km; **fixed-station check** = the fixed-9 borehole
    run `lfe_extension_wa/results/wa_fixed10/post_v2_core9/catalogue_v2.parquet` (WAVE-only, selected_v2, NOT de-duplicated,
    just produced and unvalidated → check only).
- **Bostock 2003–2013** `svi_science/results/p3_bostock/bostock_events_mapped.parquet` (fid → S0 group) and its 10 episodes
  `episodes_bostock_final.csv`. Only ETS windows were searched, so "active" history uses those 10 episodes only.

## 3. Quantities

- **Target episodes.** Primary = S3 "big" (headline ≥ 300 in the padded window): SVI 10, NW 14 (E3). Secondary = every
  MAJOR episode with headline ≥ 100 (SVI 15, NW 25).
- **Count matrix** C[p, e] (patch × episode), N_e = column sum.
- **Residual share (response)** r_pe = log((C_pe + 0.5) / (E_pe + 0.5)), E_pe = N_e · s_p^(−e), s_p^(−e) = patch's pooled
  share over all OTHER target episodes (leave-one-out). Variant: Pearson residual (C − E)/√E. When a predictor uses a past
  residual r_pk, both r_pk and r_pe use the leave-PAIR-out share s_p^(−e,−k) (Merlin's regression-to-mean fix).
- **T1 time since last active** = days from the target's padded start back to the median time of the patch's events in
  the most recent earlier episode (any episode in the list, MAJOR or minor) where the patch was active:
  C_pa ≥ 3 and C_pa ≥ 0.25 · N_a · s_p (s_p = patch's long-term share over all MAJOR episodes of that catalogue/era,
  excluding the target). Predictor = log10(T1). Rows with no activity since the history start are censored → dropped
  (variant: kept at the lower bound). Threshold variants: (≥ 1 event), (≥ 5 and ≥ 0.5 · expected).
- **T2** = r_pk of the most recent earlier target episode k ≥ 90 days before e (leave-pair-out).
- **T3** = Σ r_pk over all earlier target episodes k ≥ 90 days before e (each leave-pair-out with e).
- **K1** in_reach = patch centroid inside [tremor p05 − 5 km, tremor p95 + 5 km]; also km outside that span (0 inside).
  Tremor = PNSN epicentres ≤ 10 km of a family / A-B source, in the padded window, same frame (S3 rule). Computed for every
  target episode.
- **K2** dir_mode = mean leave-pair-out r_pk over other target episodes with the same tremor direction (sign of the tremor
  Theil–Sen speed) and |t_k − t_e| ≥ 90 d. Variant: LFE direction (S3 `direction`).
- **K3** u = patch position along the direction of travel, scaled 0 (entry end) … 1 (far end) over the patch span.
- **D1** det = log of (mean stations available to the patch's families × fraction of family-days observable) over the
  episode's padded days, minus the patch's mean over target episodes. SVI: `s0/famday_detectability.parquet` (fam_n_sta,
  observable; S13 table for the fixed run). NW: `episode_rates/wa_core_famday.parquet` (n_sta, observable).
- **Episode unevenness** U_e = Pearson χ²/df of the episode's counts against leave-one-out shares.

## 4. Tests and pass rules (fixed now)

Statistic: OLS slope of r on the z-scored predictor with episode fixed effects (r and predictor demeaned within each
episode), one predictor at a time. Referee = **order shuffle**: the episodes' count vectors (with their tremor extents
and directions) are permuted among the episode dates — big/target episodes among target slots, other episodes among
other slots — and every predictor is recomputed; 2,000 draws; one-sided p in the predicted direction (two-sided also
reported). This keeps each episode's composition (shares sum to one) and any mechanical regression-to-mean, and destroys
memory. Second control: the per-patch shuffle the caller specified (each patch's episode sequence permuted on its own;
breaks the sum-to-one, reported as a check). Calibration: 200 no-memory multinomial simulations (pooled shares, real N_e
and dates) through the same pipeline → false-positive rate of the order-shuffle test at p < 0.01 must be ≤ 3 %.

- **H1 slip budget SUPPORTED** if β_T1 > 0 with order-shuffle p < 0.01 in BOTH SVI and NW (primary targets, headline),
  AND the per-patch shuffle p < 0.05 in both, AND the fixed-station checks keep β_T1 > 0 with p < 0.05 (SVI fixed-S13; NW
  fixed-9). PARTIAL if exactly one region passes the p < 0.01 bar; NOT SUPPORTED otherwise.
- **H2 / H3**: same rule with β < 0.
- **H4 budget beyond kinematics**: increase in R² (within-episode) from adding {T1, T2, T3} to {K1 in_reach, K1 km
  outside, K2, K3} has order-shuffle p < 0.01 in both regions.
- **K1 reach** explains if β > 0, p < 0.01 against shuffling the tremor extents among target episodes.
- **K2 direction mode** explains if β > 0, p < 0.01 against shuffling direction labels among target episodes.
- **K3** reported with the direction-label shuffle (two-sided p < 0.01 = pass).
- **D1 detectability** contributes if β > 0 with p < 0.01 against permuting det within each patch across episodes.
- **E1 episode-level**: Spearman(U_e, x) for x = log N_e, |tremor speed|, tremor extent width, episode mean station count,
  SVI GNSS Tier A D; exact permutation p < 0.01 = pass (n is 10–14, so only |ρ| ≳ 0.75 can pass).
- **L2**: in_reach (0/1) on log10 T1, episode fixed effects, order-shuffle p < 0.01 (one-sided, β > 0) in both regions.
- **Bostock replication (extra, not in the pass rule)**: H1–H3 on Bostock's 10 ETS alone (history = those ETS), and on
  the joint SVI sequence 2003–2026 (Bostock 10 + our big 10; residuals within each catalogue's own long-term shares;
  T1 from big episodes only in both eras for a like-for-like history). Reported with order-shuffle p (p < 0.05 noted).
- Partial R² table (within-episode) for each predictor alone and the full model: descriptive answer to "what controls
  which patches are strong".

## 5. Controls and pitfalls

- Scrambled control (SVI c1 at the headline cut, and at P ≥ 0.5): run through §3–4 only if every target episode has ≥ 30
  control events; otherwise report the counts (S3: 0–9 per big episode at 0.75) and say it is too thin to make patterns.
- Compositional / regression-to-mean: leave-one-out / leave-pair-out shares; ≥ 90-day pairs for T2/T3/K2; the order
  shuffle carries any leftover mechanical correlation; the multinomial calibration measures it.
- Close pairs: NW 2017-02/2017-03 and 2020-09/2020-10 are < 90 d apart. T1 keeps them (it is the point of T1); a variant
  drops target episodes with another target < 90 d earlier.
- Families detecting the same event: headline is counted once; raw is a variant.

## 6. Outputs

`q1_build.py` (tables) → `patch_episode_{region}_{variant}.csv`, `episode_covariates.csv`; `q1_tests.py` →
`tests.json`, `tests_table.csv`, `calibration.json`; `q1_figures.py` → `figures/q1_resid_vs_tsince.png`,
`figures/q1_example_episodes_{svi,nw}.png`, `figures/q1_partial_r2.png`. RESULTS.md: plain answer first.

## 9. Changes after pre-registration
- 2026-10-07 ~2:40 pm PT (code, before the full run; no rule changed): the activity threshold's long-term share s_p is taken
  over ALL episodes of the era (MAJOR + minor) instead of MAJOR only, so it stays fixed under the order shuffle; minor
  episodes add little. Order and per-patch shuffles keep each era separate (E1/E2/E3, Bostock vs ours).
- Per-patch shuffle control uses 1,000 draws (not 2,000), headline designs only.
- A 20-draw debug run wrote an earlier tests.json; it was overwritten by the full 2,000-draw run (log `q1_tests.log`).
- Scrambled control: too thin at P >= 0.75 (0–9 events per big episode); every big episode has 35–150 at P >= 0.5, so it was
  run there (`q1_scrambled.py`, `scrambled_control.json`) against the data at the same cut.
- D1 slopes are "per SD" of a predictor whose SD is dominated by between-episode changes (removed by the episode
  effects), so the slope looks large; the within-episode R² is the meaningful size and is reported next to the verdict.

