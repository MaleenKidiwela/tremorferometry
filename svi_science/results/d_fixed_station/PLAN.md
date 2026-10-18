# D. Did NW LFE activity change across 2007–11 / 2012–14 / 2015–26 with the station network held fixed? — pre-registered plan

Written 2026-10-07 ~3:15 pm PDT, before any number below is computed. Already seen when writing (from the run's own
summary files, not from this study): post_v2 summary NW tremor-day net per family-day, WAVE ≥ 0.75 — production
0.471 / 0.306 / 0.246; core9 all days 0.374 / 0.237 / 0.267; core9 strict 0.271 / 0.260 / 0.300; all10 all days
0.465 / 0.280 / 0.295; all10 strict 0.341 / 0.290 / 0.330. These looks fix nothing below; the tests, cut-offs and pass
rule are chosen without reference to them.

Script `svi_science/scripts/d_fixed_station.py`. Inputs (read only):
- fixed network: `lfe_extension_wa/results/wa_fixed10/post_v2_core9/` (F9 = B001, B003–B007, B013, B017, B018, primary) and
  `post_v2_all10/` (F10 = F9 + GNW, secondary): `rates_daily_v2.parquet` (label `own` = WAVE model),
  `unit_day_strict.parquet` (n_fixed, n_fixed_live, strict), `conf/P/prod/P_*.parquet` (per-candidate times, for the
  once-per-source column).
- production network: `lfe_extension_wa/results/wa_prod/post_v2/rates_daily_v2.parquet`, label `own` (P_own_WAVE, the
  WAVE column; NOT the PNET headline), and `wa_prod/conf/P/prod/` for once-per-source;
  `svi_science/results/a_wa_between/live_core_by_day.csv` (n_core_live).
- families: `svi_science/results/a_wa_between/families_wa.csv` (fid → source_group = families_v1_3 `nw_source_group`,
  downdip_km, Hyp2000 dep).
- episodes: `lfe_extension_wa/results/episode_rates/episodes_nw_final.csv` (118 rows, 69 MAJOR).

## Definitions
- Unit = NW family (132 units). Eras E1 2007-01-01..2011-12-31, E2 2012..2014, E3 2015..2026-09-15 (as production).
- Tremor family-day = post_v2 `day_class == "tremor"` (≥ 1 PNSN epicentre ≤ 30 km; before 2009-08-06 the WECC ETS windows).
  Only scanned unit-days exist in rates_daily (fixed runs: family observable at ≥ 2 sites; production: ≥ 1 live node).
- Count = P_own_WAVE ≥ 0.75 (headline, user cut-off). Net = data − time-reversed control (main). Scrambled-control net
  (data − scrambled control) alongside. 0.5 and 0.85 as columns.
- Strict unit-day = every fixed station of the family's template set live (`unit_day_strict.strict`). "All days" = every
  scanned unit-day. n_sta_day covariate = live fraction n_fixed_live / n_fixed (fixed runs); n_core_live / era core size
  (production).
- Count once per source (column, not headline): per source group × day × series, candidates of different member families
  within ±2 s are merged (one LFE); each family-day count is scaled by merged/raw of its group-day. Source groups from
  families_v1_3 `nw_source_group` (via families_wa.csv).
- Cluster = run of NW tremor dates with gaps ≤ 5 days (one ETS or burst ≈ one cluster). Clusters are the resampling unit
  (days inside one ETS are not independent).

## Measures (each for production, core9, all10; each for all days and strict days; each raw and once-per-source)
1. **Net repeats per family-tremor-day** per era = Σ net / Σ unit-days. Ratios E1/E3, E2/E3, E1/E2.
   n_sta_day covariate: (a) table by live-fraction stratum (= 1, 0.8–<1, < 0.8); (b) day-level weighted least squares
   y_day = Σnet/Σunit-days ~ era + live fraction (weights unit-days), era effect reported as the ratio to E3 at the
   pooled mean live fraction.
2. **Net repeats per PNSN tremor epicentre ≤ 30 km** per era = Σ net / Σ n_trem30, tremor family-days from 2009-08-06
   (E1 = 2009-08-06..2011-12-31 only). Detection-independent normaliser (PNSN's own detection also changed; noted). Also
   measure 1 restricted to the same ≥ 2009-08-06 days for a like-for-like comparison.
3. **Deep / shallow ratio** per era = (measure 1 for deep families) / (measure 1 for shallow families). Deep = down-dip
   position > median of the NW families (position-based, primary); column: Hyp2000 depth > median.
4. **Per-ETS-episode series 2007–2026**: for each episode in episodes_nw_final.csv, measures 1 and 2 on the tremor
   family-days inside [start, end]. Trend over MAJOR episodes with ≥ 100 unit-days: Spearman ρ of rate vs episode start,
   permutation p (10,000 shuffles); slope of log rate vs year (weighted by unit-days) in %/yr with episode bootstrap 95 %
   CI. Primary series = core9 strict measure 1; others as columns.
5. **Production network, same measures** (WAVE ≥ 0.75), and the share of the production E1→E3 drop removed by holding the
   network fixed: 1 − ln(R_fixed E1/E3) / ln(R_prod E1/E3) (same for E2/E3).

## Statistics
- Uncertainty: cluster bootstrap (resample clusters within each era, 2,000 draws) → 95 % and 99 % CI of each era ratio.
- Era difference p: cluster permutation (pool the two eras' clusters, reassign era labels keeping cluster counts, 10,000
  draws), two-sided on |ln ratio|.
- Deep/shallow: same permutation on ln(D/S E_a) − ln(D/S E_b).

## Pass rule (frozen)
For each era pair (E1–E3, E2–E3, E1–E2), "**real change**" only if ALL hold:
- (i) core9 measure 2 (per tremor epicentre, all days) differs with permutation p < 0.01;
- (ii) the direction matches in core9 and all10 (measure 2);
- (iii) it survives strict days: core9 strict measure 2, same direction, p < 0.05.
Otherwise "**no change detectable; bound ±X %**", X = the core9 strict measure-2 99 % cluster-bootstrap interval of the
ratio (reported as the range of ln-ratio → %). Measure 1 gets the same rule as a secondary verdict (it covers 2007–2009
via WECC windows). E1 strict coverage is 11–33 % of days in 2007–2008 → any E1 statement is phrased as a bound, and
measure 2 E1 uses only 2009-08..2011.
The question "was the production 0.48/0.35/0.28 drop the network?" is answered by measure 5 (share removed) — descriptive,
no pass rule.

## Outputs
`svi_science/results/d_fixed_station/`: RESULTS.md, `era_rates.csv` (every network × day set × count variant × measure ×
era), `era_ratios.csv` (ratios, CIs, p), `strata_live.csv`, `wls_covariate.csv`, `deep_shallow.csv`, `episode_series.csv`,
`episode_trend.csv`, `checks.json` (reproduction of post_v2 summary numbers). Figures `svi_science/figures/d_era_rates.png`,
`d_episode_series.png`.

## Acceptance checks (before results are read)
- Measure 1 all days for each network reproduces post_v2 `rates_summary_v2.csv` (wave view, NW tremor, cut_0.75) exactly;
  strict reproduces `strict/rates_summary_v2.csv`.
- Once-per-source raw candidate counts (before merging) equal rates_daily n075 sums per network × series.
