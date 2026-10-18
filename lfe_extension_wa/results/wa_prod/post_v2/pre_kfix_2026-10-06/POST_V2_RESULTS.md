# Washington production — post-processing v2 (user decisions 2026-10-06)

Script `scripts/wa_prod_post_v2.py all` (CPU only, ~50 min; log `post_v2/post_v2.log`, `run_all.out`). It reads the v1
inputs and writes only into this folder. No v1 file, queue marker, GPU, download, deletion or git was touched.
"Time-reversed control" = the same detector with every template reversed in time (main false-alarm reference).
"Scrambled control" = each station site's data shifted to a different time of the same day. "False share" = control
count / data count at the same cut-off. All counts are the fixed core network of each era (E1 2007–11, E2 2012–14,
E3 2015–26), production run. The cut-off table leaves out the two sealed months (as in v1). The catalogue keeps them.

## Plain summary

- **Day classes (per family and day).** 'tremor' = at least one PNSN tremor epicentre within 30 km. 'quiet' = the PNSN
  catalogue covers the day and shows nothing within 30 km. 'unknown' = the catalogue does not cover the day. Before
  2009-08-06, days inside the WECC ETS windows stay 'tremor'. Catalogue span checked on disk (`catalogue_span.json`):
  2009-08-06 (2009 API pull) to 2026-05-27 (`catalogs/pnsn_tremor_cascadia_full.csv`). The full file starts on
  2010-01-06 and the 2009 pull ends on 2009-12-30. The 6 days in between are treated as covered, as the user's span says.
- **Headline rule (final user decision).** NW and Arm B: P_own_WAVE_PNET ≥ 0.75 on every day class. SW:
  P_own_WAVE_PNET ≥ 0.5 on tremor days only. Olympic A: position-group rates on P_nbr3_WAVE ≥ 0.75. The earlier
  'class rules' are kept as secondary columns for both models (tremor 0.75, SW tremor 0.5, NW quiet 0.85, nothing
  elsewhere).
- **Headline repeats kept, 2012–2026** (`summary_headline_false_share_2012_2026.csv`; data / time-reversed / scrambled):

  | region | tremor days | quiet days | unknown days | families ≥ 10 repeats (of kept) |
  |---|---|---|---|---|
  | NW (PNET ≥ 0.75) | 21,637 / 4 / 118 (0.02 % / 0.5 %) | 2,887 / 22 / 16 (0.8 % / 0.6 %) | 312 / 0 / 0 | 108 of 125 |
  | SW (PNET ≥ 0.5, tremor days only) | 770 / 31 / 167 (4.0 % / 21.7 %) | — (none) | — (none) | 21 of 45 |
  | Arm B (PNET ≥ 0.75) | 276 / 0 / 1 (0 % / 0.4 %) | 37 / 0 / 7 (0 % / 19 %) | 1 / 0 / 0 | 8 of 23 |
  | Olympic A groups (nbr3 WAVE ≥ 0.75) | 767 / 3 / 80 (0.4 % / 10 %) | **959 / 36 / 781 (3.8 % / 81 %)** | 11 / 0 / 6 (0 % / 55 %) | 8 of 10 groups (4 on tremor days alone) |

  Over all years (cut-off table, sealed months excluded), the headline false shares are: NW tremor 0.0 % / 0.6 %,
  quiet 0.7 % / 0.6 %, unknown 0.1 % / 0.5 %; SW tremor 2.1 % / 13.9 %; Arm B tremor 0 % / 0.7 %, quiet 0 % / 10.6 %,
  unknown 0 % / 18.8 %; Olympic A tremor 1.8 % / 10.0 %, quiet 10.0 % / 74.8 %, unknown 41.3 % / 77.5 %.
- **Olympic A off tremor days is mostly false by the scrambled control.** I read "Olympic A stays group rates on
  P_nbr3_WAVE ≥ 0.75" as applying on every day class, as for NW and Arm B. On quiet and unknown days the scrambled
  control nearly equals the data (81 % and 55 % in 2012–26). The day_class column lets you keep Olympic A to tremor
  days only. **This needs a user decision.**
- **What the 'unknown' relabel changed (NW, WAVE model).** v1 'quiet' meant no tremor within 60 km, and uncovered
  days counted as quiet.
  - At P ≥ 0.85, v1 quiet held 999 data detections. 695 of them are now 'unknown' (574 before 2009-08, 121 after
    2026-05-27). The other 304 are days with no tremor within 60 km.
  - The 30-km rule adds 1,570 detections on days with tremor 30–60 km away. So v2 'quiet' at 0.85 = 1,874 detections
    (time-reversed 82, scrambled 32 → 4.4 % / 1.7 %; 69 families).
  - On the 304 beyond 60 km, the controls give 26 % / 8 %, the same as the quiet-day waveform test.
  - At 0.75: v1 quiet 2,956 → v2 quiet 3,874 (of which 1,135 beyond 60 km) plus 1,821 unknown.
  - Quiet days within 30–60 km of tremor carry mostly real repeats. That pulls the 'quiet' false share down. For an
    inter-ETS measure, use the diagnostic class `quiet_beyond60` in `false_share_by_cutoff_v2.csv`.
- **NW quiet-day repeats kept at 0.85 (secondary WAVE class rule): 1,874** (all years, sealed months excluded; per year
  in `summary_nw_quiet_wave085_by_year.csv`). The headline PNET ≥ 0.75 keeps **3,781** NW quiet-day repeats over all
  years (time-reversed 26, scrambled 22 → 0.7 % / 0.6 %; 85 families). 678 of them fall on days with no tremor within
  60 km (3.2 % / 1.5 %).
- **Rates headline (E3 2015–2026, per family-day, k-scaled, minus the time-reversed control; 90 % day bootstrap;
  `rates_summary_v2.csv`, view `headline`, column `rule`):**
  - NW: tremor days 6.0 [5.4–6.6], quiet days 0.062 [0.053–0.072], unknown days 0.53 [0.19–0.93].
  - SW, tremor days: 11.8 [9.8–13.9] (scrambled-subtracted 9.0).
  - Arm B: tremor days 4.8 [3.6–6.0], quiet days 0.072 [0.038–0.110].
  - Olympic A per group: tremor days 92 [69–116], quiet days 5.5 [4.9–6.2]. Scrambled-subtracted, quiet drops to
    0.1 [−0.8–0.9].
  - Unscaled net counts per family-day (NW E3): tremor 0.35, quiet 0.005.
  - The k factors are large and uncertain (median k at the headline cut-off: NW 36, SW 224, Arm B 186, Olympic A 440;
    E3). Read the k-scaled rates as "array-equivalent events" with that uncertainty.
- **E1 / E2 headline rates are net counts only (no k).** The E1 / E2 truth run never got the PNET block: 0 % of its fit
  rows carry PhaseNet features (E3: 100 %). Its `P_own_WITH_PNET` is therefore the PNET model fed missing inputs, and a
  PNET k cannot be measured there. Those k are set to NaN (`k_status` in `rates_k_factors_v2.csv`). The summary gives
  net counts per family-day for every era (NW tremor: E1 0.61, E2 0.47, E3 0.35). The production PNET scores
  themselves look normal in E1 (e.g. 2008-05 NW core: WAVE ≥ 0.75 2,661 vs PNET ≥ 0.75 2,504). Only the scaling is
  missing. The WAVE secondary view keeps k for every era.
- **Step checks** (`rates_step_checks_v2.csv`; quiet-day rate and rate per tremor hour, ±21 / ±45 d, vs the same
  calendar windows in other years). 2,866 of 23,104 checks are flagged at |z| > 2: headline NW 310, SW 201, Arm B 364,
  Olympic A 507 of 2,888 each. v1 flagged 1,624 of 11,552, a similar 12–14 %. A step next to an ETS cannot be told apart
  from activity.

## Counts by region × day class (all candidates, headline view; `summary_catalogue_region_class.csv`)

| region | class | candidates | days | kept (headline) | families kept | P ≥ 0.5 | P ≥ 0.75 | P ≥ 0.85 | P ≥ 0.9 |
|---|---|---|---|---|---|---|---|---|---|
| NW | tremor | 1,149,599 | 1,582 | 31,779 | 132 | 74,947 | 31,779 | 17,530 | 8,443 |
| NW | quiet | 7,448,495 | 6,138 | 3,781 | 85 | 10,662 | 3,781 | 2,124 | 1,092 |
| NW | unknown | 1,008,856 | 966 | 1,501 | 79 | 5,031 | 1,501 | 715 | 288 |
| SW | tremor | 14,365 | 983 | 1,746 | 70 | 1,746 | 311 | 96 | 31 |
| SW | quiet | 46,861 | 5,866 | 0 | 0 | 79 | 3 | 1 | 0 |
| SW | unknown | 12,027 | 957 | 0 | 0 | 26 | 1 | 0 | 0 |
| Arm B | tremor | 161,802 | 1,134 | 618 | 40 | 1,748 | 618 | 309 | 122 |
| Arm B | quiet | 1,109,066 | 5,752 | 85 | 26 | 496 | 85 | 46 | 20 |
| Arm B | unknown | 259,012 | 964 | 32 | 18 | 248 | 32 | 16 | 3 |
| Olympic A | tremor | 4,753,962 | 1,434 | 1,203 | 23 | 18,478 | 1,203 | 629 | 346 |
| Olympic A | quiet | 38,524,050 | 6,125 | 1,355 | 25 | 140,641 | 1,355 | 141 | 38 |
| Olympic A | unknown | 5,482,927 | 934 | 596 | 21 | 36,581 | 596 | 51 | 11 |

P = P_own_WAVE_PNET for NW / SW / Arm B and P_nbr3_WAVE for Olympic A. Olympic A 'kept' counts family candidates; its
rates merge a group's members within ±2 s. The 2012–2026 rows and the two class-rule views are in the same file.

## Unvalidated columns

Every rate row also carries counts at 0.5 / 0.75 / 0.85 / 0.9 (`n050_*` … `n090_*`, data / rev / c1), the matching net
counts and k-scaled rates, and ΣP minus the time-reversed control (`sumP_scaled`). Where no rule applies (SW quiet and
unknown days; the class-rule views on unknown days and on SW / Arm B / Olympic A quiet days), these are columns only:
`status = column only` in `rates_summary_v2.csv`. The cut-off table has every region / label (own, own_pnet, nbr3, any)
/ class / cut, including 0.85, with diagnostic classes: quiet_v1, quiet_beyond60, unknown_pre2009 and
unknown_post_catalogue.

## Files (all in `results/wa_prod/post_v2/`)

- `false_share_by_cutoff_v2.csv` — data vs time-reversed vs scrambled per run / net / region / label / class / cut.
  Flags: `is_headline`, `is_classrules_wave`, `is_classrules_pnet`.
- `catalogue_v2.parquet` — every core-network data candidate (59,971,022 rows, 1.6 GB). Columns: det_id (YYYYMMDD·10⁷ +
  candidate k), time, date, family, family_uid, region, pos_group, P_own (WAVE), P_own_pnet, P_nbr3, P_any, P_max,
  net_mad, n_sta, role, day_class, unknown_kind, n_trem30, n_trem60, label_v2, cut_v2, rule_v2, **selected_v2
  (headline)**, cut_v2_classrules, rule_v2_classrules, selected_v2_wave_classrules, selected_v2_pnet_classrules.
- `rates_daily_v2.parquet` / `rates_monthly_v2.csv` — per unit (family / Olympic A group) × label × day.
  Headline = `*_v2` columns on label own_pnet + Olympic A groups. Class rules = `*_v2cr`.
  `net_count_*` subtracts the time-reversed control; `*_scr` subtracts the scrambled control.
- `rates_k_factors_v2.csv` — k per era set × cut × label × unit (shrinkage to the region median, pseudo-count 20;
  100-draw day bootstrap; `k_status`).
- `rates_summary_v2.csv`, `rates_step_checks_v2.csv`, `core_set_changes_v2.csv`, `summary_*.csv`,
  `unit_day_counts_v2.parquet`.

## Link from the main report

The v2 section with the link has been appended to `results/WA_PROD_RESULTS.md`. That copy of the report is the one
the earlier floor-0.2 run wrote at 16:08Z. The relaunched queue was still cutting windows at ~17:40Z, and it rewrites
`WA_PROD_RESULTS.md` at its very end, which will drop the section. **After that final report, the link has to be added
again.** A read-only watcher does this: `scripts/wa_prod_report_v2_watch.sh` (log `report_v2_watch.log`) waits for a
new `WA_PROD_DONE` and then runs `scripts/wa_prod_report_v2.py`, which appends the section once. If the watcher is
gone, run that script by hand. The original `wa_prod_report.py` is unchanged.
