# A — Washington-NW repeat of S7 (LFEs between ETS episodes): PLAN (frozen before any A number)

Written 2026-10-07 ~1:30 pm PT, before any A statistic was computed. Nothing about the WA between-episode pool was looked
at while writing it (only file schemas and the existing post_v2 / nw_sharing / S3 summaries).
Script `svi_science/scripts/a_wa_between.py`; outputs `svi_science/results/a_wa_between/`; figures
`svi_science/figures/a_*.png`. Rules = S7 (`results/s7/PLAN_S7.md`), unchanged unless listed under "Adaptations".

## Unit, counting and controls
- Event = one core-network candidate of a NW family in `lfe_extension_wa/results/wa_prod/conf/P/prod/P_*.parquet`
  (net core, reg NW) with `P_own_WAVE_PNET ≥ 0.75` (the post_v2 headline, = `catalogue_v2.parquet` `selected_v2` for the
  data series; checked by count). Series: data; **time-reversed control** (`rev`, main); **scrambled control** (`c1`).
  Every step is run on the controls with the same rules.
- **Counted once (headline):** dedup variant D == False (`post_v2/dedup/flags/`, README_dedup.md: another family within
  ±2 s sharing ≥ 50 % of stations, predicting the same arrival ≤ 0.5 s, with higher net_cc). Columns, never substituted:
  raw (no de-dup); **source-group once** (dropped only if its D winner is in the same `nw_source_group` of
  `data/catalogue/families_v1_3.csv`); P ≥ 0.5; P ≥ 0.85.
- **Period:** E3 2015-01-01 … 2026-05-27 (UTC days) = primary (PNSN catalogue end, `catalogue_span.json`).
  Secondary: E1 2009-08-06 … 2011-12-31 (PNSN coverage starts 2009-08-06) and E2 2012-01-01 … 2014-12-31, reported
  separately (inventory, clustering, depth, post-ETS month), never pooled with E3. Events after 2026-05-27 counted on a
  separate line, not used.
- Dates/joins are UTC days; times shown are Pacific.

## Definition of "between-episode" (fixed here)
- **BE (headline):** UTC day outside every one of the 118 NW episode windows widened by 2 days ([start − 2 d,
  end + 2 d], `episode_rates/episodes_nw_final.csv`, majors and minors), AND post_v2 day_class == 'quiet' for the
  event's family (PNSN catalogue covers the day and no epicentre within 30 km of the family = S7's `n_trem30 == 0`).
- Columns: BE_broad = outside episodes ± 2 d with any local tremor; outside episodes without the pad; BE with no
  epicentre within 60 km (post_v2 `quiet_beyond60` analogue, n_trem60 == 0).
- **Eligible family-days** (exposure for every rate and shuffle): NW family × UTC day in the period on which the family
  was scanned (the headline rows of `post_v2/dedup/rates_daily_v2_dedup.parquet`, unit_kind family, region NW; these
  include scanned days with zero candidates), satisfying the same BE condition (day_class and episode pad).

## Family position, depth and grouping (adaptation of S7's S0 families/groups)
- wa_prod fid → `family_uid` (`nw_sharing/nw_families_meta.csv`) → `families_v1_3.csv`.
- **Copies take the source's primary position:** if `nw_copy_of` is set, lat/lon/depth/quality are those of that
  primary family (array-confirmed same source). Otherwise the family's own.
- **Depth = Hyp2000 `depth_km`, used only where the (primary) location quality is A or B.** Quality-C families
  (erh up to ~80 km, some depths 4–9 km) are kept in every count but left out of the depth statistic and of the ETS
  shares it is drawn from (both sides restricted identically). Columns: Slab2 depth at the primary epicentre (A/B) and
  down-dip km; all families with Slab2 depth.
- Along-strike / down-dip frame = the S3 NW frame (`results/s3/build_info.json`, strike 344.7°, + = NW).
- "Groups" table (S7 per-group shares): the 93 `nw_source_group`s, and 5-km depth bins.

## Tests and pass/fail rules (as S7)
1. **Inventory**: data vs time-reversed vs scrambled per definition and era; families / source groups / depth bins; per
   year. **Clustered in time**: distinct UTC days of BE events vs family-preserving day shuffle (each event to a random
   eligible day of its own family, 1000×); clustered if observed < 1st percentile. Bursts = events chained by gaps
   < 24 h; number, sizes, share in bursts ≥ 3.
2. **Next ETS**: τ = days to the start of the next MAJOR (`major == True`, 69 episodes; censored after the last major
   start). Superposed-epoch bins 0–7, 7–14, 14–30, 30–60, 60–120, 120–240, > 240 d. **Ramp** if F30 (share τ < 30 d)
   > 99th pct of the shuffle AND the same holds for events ≥ 30 d after the previous major's end. Burst-level column
   (one unit per burst). **Start location**: events with τ ≤ 60 d before one of the 14 S3 big NW episodes with an LFE
   start (`s3/episodes_kinematics.csv` region NW, `lfe_start10`); statistic median |along(family) − start|; control =
   re-paired with a random other of those episodes (5000×); pass if p < 0.01 with ≥ 20 events, else "too few".
   Column: PNSN-tremor start (median along of the first 10 % of epicentres in the NW family box) for all majors.
3. **Tremor elsewhere**: share of BE events with ≥ 1 PNSN epicentre 30–100 km from the family on the same UTC day vs
   shuffle; tied if > 99th pct. Column ±6 h. Earthquakes within ±60 s (`lfe_extension/results/trackC/eq_catalog.csv`,
   covers WA) vs chance. Same-minute multi-family hits (≤ 120 s) counted.
   **GNSS stack**: S7 estimator unchanged (M5[d+1] − M5[d−5], 5-day medians, projected on own ETS azimuth, station
   mean) on the 12 WA stations of S1 (`s1/tierA_summary.json` `wa_azimuths`); north subset = those with along-strike
   > 0 in the NW frame. Cluster days (< 3 d merged; day with most events); days whose [d−6, d+6] touches a NW episode
   ± 2 d dropped; control = 1000 sets of season-matched (± 30 d of year) allowed days. Slip-like only if positive and
   > 3σ. Reference row: NW minor-episode mid-days.
4. **Down-dip**: ETS share per family = its headline events inside the majors' windows (same period, same counting).
   Statistic = mean depth of BE events (A/B families) vs 10,000 multinomial draws of the same number from the ETS
   shares; **deeper** if > 99th pct. Columns: one per burst; Slab2 depth; down-dip km; per-family Spearman of
   (BE + 0.5)/(ETS + 0.5) vs depth over families with ≥ 10 ETS events; time-reversed / scrambled control rows.

## Pre-registered confirmatory test (SVI provisional result 4)
- SVI (S7): fewer BE LFEs in the first 30 d after a major ETS ends: F30_prev 0.14 vs 0.27 expected, p_low 0.008.
- **WA statistic:** F30_prev = share of E3 headline BE events with 0 ≤ τ_prev < 30 d (τ_prev = days since the end of
  the previous major), compared with the family-preserving day shuffle (1000×; the ± 2 d pad applies to both).
- **Pass rule (fixed now):** REPEATS if F30_prev < shuffle median AND p_low = P(shuffle ≤ observed) < 0.01.
  DOES NOT REPEAT if F30_prev ≥ shuffle median, or p_low ≥ 0.05. 0.01 ≤ p_low < 0.05 in the same direction = "same
  direction, not confirmed". Fewer than 20 events with a finite τ_prev = "can't tell".
- Columns (not part of the rule): burst level; S3 big episodes (14) as "major"; BE_broad; E1/E2.

## Checks carried (not filters)
- Count of time-reversed and scrambled BE events vs data; if a control exceeds 10 % of the data at 0.75 it is said
  plainly (the cut is NOT lowered).
- Core live-station count per day (`mf/prod/run_<date>.json` members.core) for BE events vs eligible family-days.
- Core-network changes (`post_v2/core_set_changes_v2.csv`): days to the nearest change per event; share within ± 7 d
  vs eligible family-days.
- No Merlin review (caller asked to run directly); provisional.

## Adaptations from S7 (complete list)
1. Data source: WA post_v2 P files instead of `s0/analysis_table`; P_own_WAVE_PNET instead of P_same; dedup D instead of
   keep_xfam8 (+ source-group once column).
2. Main control = time-reversed (WA has it; SVI's S7 used the scrambled control); scrambled = second control.
3. Episodes = 118 NW episodes, 69 majors; period E3 2015-01 … 2026-05-27 primary, E1 (from 2009-08-06) / E2 secondary.
4. Day class from post_v2 (n_trem30 per family-day) instead of S0's n_trem30; same 30-km rule.
5. Exposure = scanned family-days from rates_daily_v2 (zero rows included) instead of `famday_detectability`.
6. Depth = Hyp2000 A/B (copies → primary) instead of Bostock depths; groups = source groups + depth bins.
7. Location test uses the 14 S3 big NW episodes; GNSS uses the 12 S1 WA stations.
8. Fixed-13 check (SVI only) replaced by the E1/E2 secondary eras and the station-change / live-station checks.
