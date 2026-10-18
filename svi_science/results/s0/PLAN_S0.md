# S0 — shared foundations table: PLAN (frozen before any per-group LFE count)

Written 2026-10-07 ~9:00 am PT. Nothing below has been computed per group yet; the family grouping uses family
locations only. Script: `svi_science/scripts/s0_build.py`. Outputs only in `svi_science/results/s0/` and
`svi_science/figures/s0_*.png`.

## Inputs (read-only)
- Candidates: `lfe_extension/results/refine/conf/P_by_candidate.parquet` (runs prod and prod_S13; series data and c1 =
  scrambled control).
- Episodes: `lfe_extension_wa/results/episode_rates/episodes_svi_final.csv` (79 episodes, 39 MAJOR, frozen ladder rule) and
  the per-family exposure `exposure_svi.parquet` (PNSN epicentres ≤ 30 km per family-day).
- Day station sets: `episode_rates/svi_day_stations.csv` (from `refine/feat/<run>/<date>.json`).
- Family-day scan record: `refine/feat/<run>/<date>_fd.parquet` (per family-day network MAD, live template stations / sites).
- Completeness model: `lfe_extension/results/trackD/mc_family_day.parquet` (Track D D1(b), see below).
- Family positions: `lfe_catalog_repair/results/e5r_families_2012.csv` (Bostock hypocentres, mag_med); carriers =
  `e5r_template_qc_2012.csv` usable & template stations.
- Slab2 (`data/slab_models/slab2/cas_slab2_dep_02.24.18.xyz`) for strike and depth contours.

## Tables
1. `analysis_table.parquet` — every candidate of run prod (data + scrambled control); `analysis_table_S13.parquet` — the
   same for run prod_S13 (fixed 13 stations, to 2021-03-31). Columns: time (UTC), family, P_same, series, run, station set
   (S13/S14/S15 and the station list id), n_sta_day, keep_xfam8, quiet, n_tremor_day, n_trem30 (family-day PNSN count);
   episode (assigned by the family's own window = union of [d − 2, d + 2] over its exposure days inside the episode span,
   exactly as `episode_rates.py`), major, in_frozen_major, episode start/end, family window start/end; time since the
   episode's regional start and since the family's first exposure day in that episode (days); days since the end of the
   previous MAJOR episode and to the start of the next (from the candidate time) plus the episode's own spacing to the
   previous/next MAJOR episode; detection-problem flags; detectability per family-day (below).
2. `famday_detectability.parquet` — per run × family × day: n_sta_day, family live template stations and sites
   (`_fd`), network MAD (data series) and its log10 ratio to the family's median over the same run and station set
   (u_mad; a noisier cc series needs a bigger event), Track D Mc / P(detect | M 1.8) where available.
3. `families.csv` — per family: Bostock lat / lon / depth, mag_med, Slab2 depth at the epicentre, along-strike and
   down-dip coordinates, carrier stations (all, and ∩ S13), group id, totals at P ≥ 0.75 raw and de-duplicated
   (data and scrambled control, both runs).

## Detection-problem flags (carried, never used to drop)
- `flag_ep_2017_03`: candidate in episode 2017-03 (B010 out on all 18 days, VGZ swap 2017-03-03).
- `flag_ep_2025_10`: candidate in episode 2025-09b1 (VGZ loaded 2 of 30 days, B926 12, B001/B005 19).
- Station steps as epochs and as "within ±45 d" flags (the episode_rates convention): PGC response sub-change 2015-01-12,
  PGC BH→HH 2017-08-10, NLLB BH→HH 2017-02-23 (gain step ×1.77), VGZ sensor swap 2017-03-03.

## Completeness model
`lfe_extension/scripts/d1_mc_model.py` (Track D D1(b)) is a probit P(detect | family, day, M) model fitted to the Track C
detector (T_core / T_ext MAD tiers), validated on the held-out 2010 year (bar ±0.05; met, max |Δ| 0.024), with output
`mc_family_day.parquet` for 2014-01-01 .. 2026-09-15 (production effective station set). It is NOT calibrated for the
refine P ≥ 0.75 catalogue used here and 48 % of production family-days lie outside its calibration range. Rule: carry its
Mc_T_core and Pdet_T_core_M1.8 as covariate columns for run prod only; do not use it as a completeness correction. It is
not re-run here.

## Coordinate frame (fixed here)
- Origin: mean epicentre of the 115 families. Local km: x = east, y = north (equirectangular at the origin latitude).
- Strike: circular mean (axial) of the Slab2 local strike at the 115 family epicentres (strike = down-dip azimuth − 90°,
  down-dip azimuth from the gradient of the Slab2 depth grid, bilinear at each epicentre).
- along_strike_km = projection on the strike azimuth, oriented so positive = north-west; downdip_km = projection on the
  perpendicular, positive = down-dip (north-east). The principal axis of the family cloud is reported for comparison only.

## Family groups (pre-registered; locations only)
- Complete-linkage agglomerative clustering of the 115 families on epicentral distance (km) between Bostock positions
  (horizontal only: depth errors are larger and depth follows down-dip).
- Cut D tried in order 10, 12.5, 15, 17.5, 20, 25, 30 km (every pair inside a group ≤ D km apart). For each D: any cluster
  with < 5 families is merged into the cluster with the nearest centroid (smallest first, repeated until all ≥ 5).
  Chosen D = the smallest giving 8–12 groups after merging; if none does, the D whose group count is nearest to 10
  (ties: smaller D).
- Group ids numbered by the centroid's along-strike coordinate (1 = south-east end).
- The assignment is written to `families.csv` and `family_groups.csv` before any count per group is made.

## Counts per group × MAJOR episode (after the grouping is frozen)
- Unit: de-duplicated repeats (keep_xfam8 == 1), P ≥ 0.75, data minus scrambled control, summed over each family's own
  episode window; primary series fixed-S13 (prod_S13 before 2021-04-01, prod after — as episode_rates), production
  effective set as a second table. Raw (not de-duplicated) and P ≥ 0.5 as columns.
- Floors (Revision 1): ≥ 30 per family × episode cell, ≥ 100 per group × episode, ≥ 300 per episode. Reported: how many
  cells meet each floor, and how many group × episode cells meet ≥ 100 inside episodes that meet ≥ 300.

## Acceptance
The table must reproduce `episode_rates` `fam_episode_svi_S13fixed.csv` and `fam_episode_svi_eff.csv` raw counts exactly
(n_data, n075_c1, N, n050_data, n050_c1, n085_data, n090_data for every family × episode). Any difference is stated.

## Figures
`svi_science/figures/s0_family_groups_map.png` (families coloured by group, Slab2 depth contours);
`svi_science/figures/s0_counts_group_episode.png` (heatmap, group × MAJOR episode, de-dup P ≥ 0.75 fixed-S13, floors marked).

## Grouping result (appended 2026-10-07 ~9:02 am PT, locations only, before any per-group count)
- Frame: Slab2 mean down-dip azimuth at the families 51.2° (spread 7.4°) → strike azimuth 321.2° (N39W); along-strike
  positive toward N39W, down-dip positive toward N51E. The family cloud's principal axis is 136.2° / 316.2° (5° from Slab2).
- Cut tries (clusters before → groups after merging < 5): 10 km 35→15, 12.5 km 31→17, 15 km 27→17, 17.5 km 21→14,
  20 km 19→14, **25 km 14→12 (chosen: smallest cut with 8–12 groups)**, 30 km 11→10.
- 12 groups of 5–13 families (`groups.csv`, `family_groups.csv`). Merging small clusters makes two groups wider than the cut:
  group 7 (max pair 26 km) and group 12 (40 km; north end, sparse). Numbered 1 (south-east) … 12 (north-west).

---------------------------------------------------------------------------------------------------------------------------
# Results (added 2026-10-07 ~9:05 am PT)
- `analysis_table.parquet` (run prod, 1,064,518 candidates: 658,099 data, 406,419 scrambled control) and
  `analysis_table_S13.parquet` (prod_S13, 743,746: 458,821 / 284,925); together = all 1,808,264 rows of P_by_candidate.
  `famday_detectability.parquet` 845,408 run × family × day rows (116 families scanned, 115 with candidates).
- **Acceptance: exact.** Raw counts per family × episode (all 9,085 rows = 115 × 79, cuts 0.5/0.75/0.85/0.9, data and control,
  N) equal `fam_episode_svi_S13fixed.csv` and `fam_episode_svi_eff.csv` with 0 differences (`acceptance_fam_episode.json`).
- Completeness model: Track D Mc columns carried for run prod (524,063 of 538,356 family-days matched); its Mc_T_core is
  essentially unrelated to this run's own family-day MAD ratio (Spearman −0.05), another reason it stays a column only.
- Counts (`counts_group_episode.csv`, `counts_summary.json`), fixed-S13, de-dup, P ≥ 0.75: 7,750 net repeats in the 39
  MAJOR episodes (raw 10,172). Episodes ≥ 300: 10 (2015-12, 2018-06, 2019-08, 2019-11, 2020-11, 2021-01, 2021-10, 2022-03,
  2023-04, 2024-08). Group × episode cells ≥ 100: 24 of 468, all inside ≥ 300 episodes; ≥ 30: 82. Family × episode cells
  ≥ 30: 49 (raw 91; eff raw 111 = Merlin's number). Usable cells per group: 1, 6, 2, 0, 3, 0, 5, 4, 0, 0, 3, 0 (groups 1–12);
  only groups 2 and 7 have ≥ 5 usable episodes. At P ≥ 0.5 de-dup: 12 episodes ≥ 300, 52 usable cells, 7 groups with ≥ 5.
- Merlin counted 11 episodes ≥ 300; here 10 with eff raw 0.75 net of the control in family windows (definition difference).
