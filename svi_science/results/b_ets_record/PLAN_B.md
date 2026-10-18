# B. ETS over the long record — southern Vancouver Island (SVI) and Washington-NW — PLAN

Written 2026-10-07 ~1:25 pm PT, before any interval, size, seasonal or lead/lag number was computed. Rules below are
frozen; any change goes into a "Revisions" section with time and reason. Script: `svi_science/scripts/b_ets_record.py`.
Times are Pacific (PT); day-of-year and months are taken from PT onset times.

## 0. Inputs and definitions

- **Episode lists** (frozen rule of `lfe_extension_wa/results/episode_rates/PLAN.md`, Revision 1: 5 % of families exposed,
  gaps ≤ 5 d bridged, span ≥ 5 d): `episodes_nw_final.csv` (2007–2026: three WECC windows 2007-01, 2008-04, 2009-04 that are
  major ETS by definition, then PNSN 2009-08 → 2026-05) and `episodes_svi_final.csv` (2014–2026).
  **SVI extension 2010–2013**: the same code (`episode_rates.find_episodes`, same rule, same 115 families, same
  `exposure_svi.parquet`) run from 2010-01-06. Check before use: rerun from 2010 must reproduce every 2014+ episode of the
  final file exactly (start, end, family-exposure days); if not, SVI stays 2014+ only and this is stated.
- **Tremor footprint** (detection-independent of the LFE work; same as S3): PNSN epicentres (`catalogs/pnsn_tremor_cascadia_full.csv`)
  within 10 km (epicentral) of any family of the region (SVI: the 115 Bostock families of S0; NW: the 132 NW families of
  `wa_prod/setup/families.csv`, the set the episode list was built from). An epicentre within 10 km of both regions is given
  to the region of its nearest family. Column: 30 km footprint (the exposure radius), same nearest-region rule.
- **Episode window** for counting: S3's padded window (start − 3 d … end + 3 d, padding capped at half the gap to the
  neighbouring episode of the same list; all listed episodes, major or not, are used as neighbours).
- **Per-episode quantities**: tremor count in the footprint (n_trem10; column n_trem30); onset = PT time of the 5th
  percentile of footprint epicentres in the window (columns: 10th percentile, first epicentre, list start day); end = 95th
  percentile; duration = end − onset (days; column span_days of the list); family-exposure days (list).
  LFE size columns: SVI = S3 `n_headline` (P_same ≥ 0.75, counted once, keep_xfam8) and `n_c1_075` (scrambled control),
  2014+; NW = catalogue_v2 `selected_v2`, counted once (dedup variant D == False), recomputed here for every era in the same
  window, with era label (E1 2007–11, E2 2012–14, E3 2015+; LFE detectability differs between eras, so LFE sizes are compared
  only within E3), checked against S3 `n_headline` for E3 majors (must match within 1 %); `n_sources` = number of
  `families_v1_3.csv` NW source groups (`nw_source_group`) with ≥ 3 counted-once events. SVI GNSS: S1 `D` (mm), `sigma`, `SN`
  (Tier A), `tierB_slip_mm` (column).
- **Big ETS (primary)**: per region, split all PNSN-era listed episodes into two classes by log10(n_trem10) with the optimal
  1-D two-class split (minimum within-class sum of squares, exhaustive); the upper class = big. WECC windows are big by
  definition (no PNSN tremor; they enter interval statistics only, as a column "with/without WECC").
  Alternatives (columns): (a1) family-exposure days ≥ 0.5 × N_fam × 14 (SVI 805, NW 924); (a2) the K largest n_trem10 with
  K = number of whole years covered (one per year); (a3) S3 "big" (LFE ≥ 300 counted once; SVI 2014+, NW 2015+);
  (a4) "cycles": primary big episodes whose onsets are < 120 d apart are merged into one cycle (onset = first, size = sum).

## 1. Questions, statistics and pass rules

A result is **robust** if it passes (p < 0.01, same sign) with the primary big list AND with at least 2 of the alternatives
a1, a2, a4 (a3 is too short for p < 0.01 and is shown only); **suggestive** if p < 0.05 with the primary only;
otherwise **not supported**. Every permutation uses 10,000 draws, fixed seed 0.

**Q1 intervals.** Interval = onset_i+1 − onset_i (days) between consecutive big ETS. Mean, median, SD, CV, range.
Regularity: CV compared with the scrambled control = the same number of onsets drawn uniformly over the same span (Poisson);
"quasi-periodic" if observed CV < 1st percentile of the control. Drift: Spearman ρ(interval, interval mid-year), permutation p;
and Mann–Whitney of intervals with onset before vs after 2018-01-01. Literature check: NW/SVI ~14 months (Rogers & Dragert
2003; Brudzinski & Allen 2007) — reported, not tested.

**Q2 sizes.** Sizes = n_trem10 (primary), duration, family-exposure days, LFE count (SVI 2014+; NW E3), GNSS D (SVI).
(i) Trend with time: Spearman(size, onset year). (ii) Slip-predictable: Spearman(size_i, interval before i).
(iii) Time-predictable: Spearman(size_i, interval after i). Permutation p (sizes shuffled = scrambled control).
Because PNSN tremor counts can drift with the network, (ii)/(iii) are repeated with sizes replaced by their rank residual after
removing a linear trend in onset year (column `detrended`); a pass needs both versions to agree in sign.

**Q3 seasonality.** (a) Onsets: Rayleigh test on PT day-of-year of big-ETS onsets per region; plus a scrambled control that keeps
the region's interval set (intervals shuffled in order, random start phase) and recomputes the mean resultant length R;
seasonal onsets need Rayleigh p < 0.01 AND control p < 0.01. (b) Between-episode activity: footprint tremor on days outside
every listed episode window (major or not); monthly rate = count / between-days per (year, month). Statistic = amplitude of the
annual harmonic of the 12-month pooled rate profile divided by its mean. Scrambled control = each year's 12 monthly values
rotated by an independent random number of months. Pass p < 0.01. Columns: share of between-days with ≥ 1 epicentre;
SVI between-episode LFEs from S7 (`results/s7/between_events.csv`, counted by day). (c) Detection caveat (diagnostic only, not
a test): the same harmonic for the daily tremor rate on big-ETS days (if detection were seasonal, winter ETS days should show
fewer epicentres per day), and for the whole PNSN catalogue's between-episode days outside both footprints.

**Q4 cross-border order.** Pair a big SVI and a big NW episode when their [onset, end] intervals overlap or onsets are within
30 d. Delay = onset_SVI − onset_NW (days; positive = Washington first). Typical delay = median, bootstrap 90 % CI. Order
"consistent" if ≥ 80 % of pairs have the same sign and two-sided binomial p < 0.05. Coincidence control: SVI onset sequence
shifted by a random offset (circular over the common span, 10,000), number of pairs within 30 d; pairing "beyond chance" if
p < 0.01. Columns: 10th-percentile onset; S3 LFE leads (`cross_region_overlaps.csv`) and S3 front classes for 2015+ pairs.
Unpaired big episodes in each region are listed.

**Q5 relative clocks.** (i) Spearman(delay, year) over pairs, permutation p < 0.01 → the clocks drift relative to each other.
(ii) Counts of big ETS per region over the common span 2010-01-06 → 2026-05-27 and Mann–Whitney of the two interval sets
(p < 0.01 → different clocks). (iii) Figure: cumulative big-ETS number vs time, both regions.

## 2. Outputs
`episodes_b.csv` (every listed episode, both regions, all columns and big flags), `intervals.csv`, `tests.json`,
`pairs.csv`, `seasonal_monthly.csv`, `RESULTS_B.md`; figures `svi_science/figures/b_timeline.png`, `b_interval_size.png`,
`b_seasonal.png`, `b_crossborder.png`.

## 3. What is known (stated, not tested)
~14-month ETS recurrence under southern Vancouver Island / northern Washington (Rogers & Dragert 2003); along-strike segments
with different recurrence (Brudzinski & Allen 2007); a seasonal (hydrological-load) modulation of slow slip was proposed by
Lowry (2006) and remains debated; PNSN tremor detection depends on seasonal station noise.

## Revisions

### Revision 1 — 2026-10-07 ~1:45 pm PT, after the build step, BEFORE any interval / size / seasonal / lead-lag test
- **SVI extension check.** The 2010 rerun reproduces every 2014+ episode of the final file exactly (start, end, family-exposure
  days). One extension episode (2013-12-25 → 2014-01-02, 75 family-exposure days, not major) runs across 2014-01-01; it does
  not overlap any final-list episode (the final list starts 2014-01-01 and its 2014 stub is < 5 d, so it never formed). The
  "no crossing episode" clause was written to catch double-listing, which does not occur; the extension (2010-01-06 → 2013-12)
  is therefore used.
- **Primary big-ETS rule.** The two-class split of log10(n_trem10) over ALL listed episodes separated tiny bursts
  (≤ ~100–180 epicentres) from everything else: 49 of 79 SVI (2014+) and 55 of 115 NW episodes came out "big", ~3 per year,
  i.e. it does not isolate ETS. Replaced (before any test) by the same optimal two-class split of log10(n_trem10) applied only
  to the episodes already MAJOR under the frozen episode rule (upper class = big). Alternatives a1–a4 unchanged; the old
  all-episode split is kept as a column (`big_split_all`) and not tested.
