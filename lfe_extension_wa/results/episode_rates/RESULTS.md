# Per-family × per-ETS-episode LFE rate consistency — WA-NW and southern Vancouver Island

Run 2026-10-06, ~8:20–8:45 pm PT. Spec: `MERLIN_DESIGN_2026-10-06.md`. Frozen rules: `PLAN.md` (written before any count; one
revision, before counting, recorded there). Script: `lfe_extension_wa/scripts/episode_rates.py` (sub-commands exposure, episodes,
count_wa, count_svi, table, checks, check_a_boot, consistency, figures). Log: `episode_rates.log`. Every number below comes from
a file in this folder. All results are provisional.

"Scrambled control" = each station site's data shifted to a different time of the same day. "Time-reversed control" = the same
detector with every template reversed in time. N = repeats at P ≥ 0.75 minus the main control (WA: time-reversed; SVI: scrambled).
X = sum of PNSN tremor epicentres ≤ 30 km of the family over its episode window. D = the family's tremor days in the episode.

## The answer

1. **Within one station era, most WA-NW family rates are not consistent from episode to episode beyond counting noise.**
   Era E3 (2015–2026), 47 MAJOR episodes, 89 of 132 families testable (sum N ≥ 20). 58 of the 89 are over-dispersed
   (BH q < 0.05) and 31 are consistent. Median χ²/df of the share test is 4.0 (quartiles 2.2–7.8). Each family's share of
   its 40-km neighbourhood's total varies by about 70 % between episodes beyond counting noise (median 0.71, quartiles
   0.60–0.88). E2 and E1 look the same: 48 of 60 and 63 of 81 over-dispersed, median χ²/df 4.6 and 3.8.
   (`consistency_summary.json`, `consistency_families.csv`)
2. **SVI (fixed 13-station set, 2014–2026) is much closer to consistent.** 39 MAJOR episodes, 75 of 115 families testable:
   48 consistent, 27 over-dispersed. Median χ²/df is 1.8 (quartiles 1.1–3.8), with about 43 % excess scatter.
3. **Tremor exposure explains almost none of the over-dispersion.** "Exposure-explained" (Poisson dispersion against X
   below 2 while share χ²/df is above 3) applies to 3 of 89 NW E3 families, 1 of 75 SVI families, 2 of 60 in E2 and 12 of 81
   in E1. The median Poisson dispersion against X is itself 4.2 (NW E3) and 3.2 (SVI).
4. **The over-dispersion is spatially coherent. In WA-NW, shared stations also carry part of it.**
   - **Distance.** Residuals of nearby families rise and fall together, and the correlation dies out by 20–40 km.
     - NW E3: 0.31 at < 3 km, 0.16 at 6–10 km, 0.06 at 10–20 km, −0.03 at 20–40 km.
     - SVI: 0.73, 0.39, 0.19, −0.05 in the same bins.
     - Distance coefficient (permutation over families, 1,000 draws): NW E3 −0.11 (p = 0.001); SVI −0.52 (p = 0.001).
   - **Carrier overlap, controlling for distance.** SVI: no effect (coefficient −0.035, p = 0.48). NW E3: yes (coefficient
     +0.23, p = 0.001). In NW E3 the 653 pairs with the highest template-weight overlap still correlate at 0.06 when > 20 km
     apart, against about −0.02 for the other pairs.
   - **Merlin's mind-changer is not met.** It needed the distance term to vanish while the carrier term survived. The distance
     term does not vanish, so the NW rates are not purely station-driven, but they are partly station-coupled.
   - **Per-cell labels (|residual| ≥ 2).**
     - NW E3, 753 cells: 122 spatially coherent, 8 station-coincident, 3 both, 620 isolated.
     - SVI, 339 cells: 168 spatially coherent, 0 station-coincident, 171 isolated.
   - **Caveats.**
     - Correlation at < 3 km can also come from co-located families detecting the same events (a detection effect, not
       physics).
     - The share test normalises within 40 km, so a change that hits every family alike, such as a regional station swap,
       cancels in the residuals. It shows only at the regional level (point 5).
   - (`residual_pair_correlations.csv`, `consistency_rows.csv`, `fig3_residual_corr_vs_distance.png`)
5. **At the regional level, N tracks tremor, with two SVI exceptions that coincide with station problems.**
   - Spearman(N, X) across MAJOR episodes is 0.92 for NW E3 and 0.81 for SVI.
   - SVI 2017-03 has 37 repeats for X = 196,559, 20× below the SVI median N/X. B010 was missing on all 18 days, and the VGZ
     sensor swap (2017-03-03) and NLLB BH→HH (2017-02-23) came just before.
   - SVI 2025-10 has N/X 5× below the median. VGZ was loaded on 2 of 30 days, B926 on 12 and B001/B005 on 19
     (`svi_day_stations.csv`).
   - These two look like detectability (provisional), not a change in the families.
   - WA-NW N/X in the 15 large E3 episodes (D > 1,000) runs 0.0049–0.0108. It has been lower since 2024: 2024-01 0.0070,
     2024-08 0.0049, 2025-04 0.0053, 2025-09 0.0085, against 0.0078–0.0108 for 2015–2023. A possible slow detectability drift
     is not separated from activity here.

## Sanity checks (pre-registered; verbatim outcomes, `sanity_checks.json`)

- **(a) WA build/test ratio: FAIL by the pre-registered bar.**
  - Pooled over 64 out-of-sample MAJOR episodes: 1.094, inside [0.85, 1.15].
  - 4 of the 19 large episodes (≥ 400 half-counted repeats) have Poisson 95 % intervals that miss the band: 2012-04, 2019-08,
    2021-09, 2025-09 (`check_a_build_test.csv`).
  - Diagnostic, not pre-registered: a day bootstrap, which allows for LFE bursts clustering in time. With it, all 21 large
    episodes overlap the band, and the pooled ratio is 1.096 [1.02, 1.17] when episodes are resampled
    (`check_a_build_test_dayboot.csv`, `check_a_dayboot.json`).
  - The Poisson interval is too narrow for bursty counts. The pooled value does sit about 9 % above 1, an unexplained
    asymmetry.
  - In-sample episodes: 2008-04 2.05 and 2010-08 2.35, as Merlin found (2.11 / 2.39).
  - Test-half share on out-of-sample days: 0.468. The empirical factor is therefore 2.135, against 2.034 from block coverage.
    2.034 was used, as planned.
- **(b) Control ≤ 2 % of data in every MAJOR episode.**
  - NW: PASS. Maximum time-reversed share 0.36 %, pooled 0.03 %. The scrambled share reaches 33 % only in tiny episodes
    (1 of 3 repeats).
  - SVI fixed-S13: FAIL in 2 low-count episodes: 2017-10 (2 of 70 = 2.9 %) and 2019-05 (1 of 40 = 2.5 %). Pooled 0.56 %.
  - SVI eff: FAIL in 2018-05 (4 of 59 = 6.8 %). Pooled 0.49 %.
  - These episodes stay in the tables (control share is a column).
- **(c) Spearman(log N, log X) ≥ 0.5 within an era.**
  - NW E3: 0.918 (47 episodes) PASS. E2: 0.923 (13). E1: 0.886 (6).
  - SVI fixed-S13: 0.814 (39) PASS.
  - SVI eff by set: S15 0.452 (8, FAIL, secondary only), S14 0.791 (11), S13 0.772 (19).
  - No primary region failed, so step 5 ran for both.
- **Named cases** (`sanity_checks.json` → named):

  | case | N | X | D | N/D (era median) | N/X (era median) | N/X rank |
  |---|---|---|---|---|---|---|
  | SVI 2017-03 (2017-03-18..04-04), fixed-S13 | 37 | 196,559 | 1,627 | 0.023 (0.116) | 0.00019 (0.0037) | 5 of 39, lowest first |
  | SVI 2025-10 (episode "2025-09b1", 2025-09-20..10-19) | 226 | 290,270 | 2,223 | 0.102 (0.116) | 0.00078 (0.0037) | 9 of 39 |
  | WA 2024-08 (2024-08-08..08-26), E3 | 672 | 137,384 | 1,849 | 0.36 (0.30) | 0.0049 (0.0069) | 14 of 47 |

- **Merlin's NW reference (net per family-tremor-day).** My N/D is 0.62–0.93× his numbers, median 0.81× (`merlin_reference_compare.csv`).
  - 2012-09 0.89 vs 0.97. 2024-08 0.36 vs 0.40. 2014-12 → episode 2014-11: 0.32 vs 0.52. 2023-04 → 2023-03: 0.47 vs 0.71.
  - Part of the gap is definition: episode windows here, calendar months there.
  - **Side finding, provisional — it needs a check of `wa_prod_post_v2.summary()`.**
    - `rates_daily_v2.parquet` has no row for a family-day with zero candidates in all three series, although the detector
      scanned that day (`mf/prod/famday_*` observable = True).
    - In E3 NW, 18 % of observable tremor family-days and 24 % of observable quiet family-days are missing.
    - Per-family-day rates built from those rows, including `rates_summary_v2.csv` and the reference above, are therefore
      overstated by about 1.2–1.3×.
  - Episode 2014-11 is low in my table for a second reason: its first half (November, 0.32 per family-day in the monthly
    table) is weak, and December alone (0.69) is what Merlin's month picks up.

## Episodes and acceptance (before counting)

- Frozen rule (10 % of families exposed, gaps ≤ 3 d, ≥ 5 d, MAJOR 0.25 × N_fam × 14) missed one reference entry in each
  region (`acceptance_*_frozen.csv`):
  - NW literature 2009-08-05..08-25. The PNSN 2009 pull starts 2009-08-06 and is sparse.
  - SVI band-48.5 2020-09-13..09-29. A 4-day burst, shorter than the 5-day minimum.
- The pre-registered ladder ended at fraction 5 %, gaps ≤ 5 d, MAJOR 0.15 (`PLAN.md`, Revision 1).
  - NW: 69 MAJOR episodes, 18 of 18 reference entries accepted.
  - SVI: 39 MAJOR, 12 of 13 accepted; 2020-09 is still missing.
  - Every row carries `in_frozen_major`.
  - Robustness on the frozen subset: NW E3 72 of 88 families over-dispersed, median χ²/df 4.3. SVI 32 of 75, 2.3. The
    answer is unchanged.
- Exposure checks (`exposure_checks.json`):
  - WA `tremor_family_day.parquet` covers all 132 NW families, 2009-08-06..2026-05-27, and agrees with `rates_daily_v2` n_trem30
    on all 603,874 covered rows (0 mismatches).
  - PNSN reaches 49.19 N in every year 2010–2026 and covers all 115 SVI families every year (`svi_pnsn_coverage_by_year.csv`).
  - There is no exposure after 2026-05-27, so no episode there.
- Counts check: the WA counts I recomputed from `conf/P/prod` equal `rates_daily_v2` (own_pnet) at every cut and series,
  with 0 differing rows (`counts_wa_check_vs_rates_daily.json`).
- In-sample:
  - WA: V2008/V2010 days use the test half × 2.034 (episodes 2008-04 and 2010-08, label `insample_halves`). Other labels:
    2007-01 sealed, 2009-04 semi, fit months, 2011-04 excluded (not tested).
  - 2012-09: 43 NW families carry in-sample T2012 templates (`tpl_insample` in `conf/P`; label "t2012 templates"). Also 4
    families in 2013-09 and 1 in 2014-11, from step-4b templates.
  - **SVI: no extension day is in-sample.** `refine/day_table.csv` roles from 2014-01-01 are all "prod" (4,641 days); the
    truth/fit days are 2010–2013.

## Tables and figures (all in this folder)

- `episodes_{nw,svi}_final.csv`: episode lists. `acceptance_*.csv`: acceptance checks.
- `fam_episode_nw.csv`, `fam_episode_svi_S13fixed.csv`, `fam_episode_svi_eff.csv`: per (family, episode), identical columns.
  - Counts: N, N_scr, N050/N085/N090, sumP_net, n_data, ctl_main, ctl_scr, N_expdays.
  - Exposure: X, D, D_obs, X_obs.
  - Activity: active_days, onset_lag_d, offset_lag_d.
  - Stations: live_carrier_days, mean_sites, station_change (+ list), station_change_principal.
  - Labels: insample, t2012_templates, sample_label, excluded, era.
  - WA E3 only: k_E3 and kN_E3 (corrected k, cross-era bridging column).
- `regional_*.csv`: per-episode regional totals and control shares.
- `consistency_families.csv`: per family × block × variant. Columns: χ²/df, bootstrap p (2,000 draws), BH q, Poisson
  dispersion, excess scatter, labels. Variants: main, frozen_major, N_expdays, ctl_main, ctl_scr. No family reaches
  sum ≥ 20 on either control, so the control rows are listed but none is tested.
- `consistency_rows.csv`: per (family, episode) residuals and labels. The WA rows also carry `label_broadflag`.
- `residual_pair_correlations.csv`; `svi_fixedS13_ratio*.csv`: 28 of 310 SVI (family, episode) ratios are "set-affected".
- Figures:
  - `fig1_family_episode_heatmap.png`: per-family N/X relative to the family's own rate, NW E3 and SVI.
  - `fig2_regional_N_vs_X.png`: regional N against X, NW coloured by era, SVI alongside.
  - `fig3_residual_corr_vs_distance.png`: residual correlation against distance and carrier overlap.

## Deviations from the design / plan

1. MAJOR rule loosened by the pre-registered ladder (Revision 1). One SVI reference entry (2020-09) is still not covered.
2. WECC windows are in `wa_prod_post.WECC`; `wa_common` has none.
3. **WA carrier sets are degenerate.** Every NW family has templates at every core station (about 19 in E3), so set Jaccard
   is about 1 and the pre-registered station-change flag fires on 85–100 % of rows.
   - Used instead:
     - Carrier overlap: weighted Jaccard of the template weights `w` (`tpl_library`).
     - Station-coincident labels: "principal carriers", the core stations holding ≥ 50 % of the family's template weight,
       1–6 stations, median 4.
   - The broad flag and its labels are kept as columns.
4. The `N_span` column was not computed; `N_expdays` (repeats on exposure days only) was added instead as the robustness count.
5. Check (a): the day-bootstrap diagnostic was added because the Poisson interval ignores burst clustering. The
   pre-registered verdict stays FAIL.
6. The E1 in-sample episodes enter the binomial test with ×2.034 counts. Their variance is understated about 2× there.
7. SVI has no per-family scan record. Live carriers = 2012 usable-template stations ∩ the stations loaded that day.
8. The WA fixed-10 subset re-run was not done. It needs a new matched-filter pass and is out of scope.
9. The E1 episode 2009-08 is cut short by the catalogue start (N = 3) and is still tested in E1.

## What this does and does not show

- **It shows:** per-family repeat counts within an era scatter well beyond counting noise in WA-NW (≈ 4× Poisson variance),
  and much less in SVI. Most of the coherent part is shared by families within about 20 km.
- **It does not separate** physical patchiness of slip from co-detection between nearby families. For WA-NW it also leaves
  open how much is station-specific coupling, given the significant carrier-overlap term.
- **Next tests that would settle it:**
  - Repeat the share test after merging families within 3 km (removes co-detection).
  - The WA fixed-10 subset run (removes station turnover).
  - Fix the zero-candidate denominator in `rates_daily_v2`.
