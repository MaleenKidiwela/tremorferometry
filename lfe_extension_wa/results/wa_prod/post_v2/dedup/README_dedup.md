# WA post_v2 — "counted once" columns (cross-family duplicate flags)

Run 2026-10-07, ~10:30–10:45 am PT. Script: `lfe_extension_wa/scripts/wa_post_dedup.py`. It is new; `wa_prod_post_v2.py`,
`wa_prod_common.py` and `wa_common.py` were imported read-only and not edited or re-run. Every output is a new file in this
folder. **Nothing is deleted and no row is dropped from any file.** Duplicates are flag columns plus extra count columns.
**The raw counts stay the headline.**

## Rule
- Same rule as the SVI `keep_xfam8` flag. The shared function `_wa_share` comes from `svi_science/scripts/s0_s0b_dedup.py`.
- The rule works per UTC day × series. Data, time-reversed control and scrambled control are treated the same way.
- It runs over all families in the core-net rows of `conf/P/prod/P_*.parquet`, not only NW. These are the rows catalogue_v2
  and rates_daily_v2 were built from.
  - `catalogue_v2.parquet` holds only the data series and no `net_cc`, so it cannot carry the rule.
- A candidate is **flagged** as a duplicate when all of these hold:
  - a candidate of another family lies within ±2 s;
  - the two families share stations: |A ∩ B| / min(|A|, |B|) ≥ 0.5 over the `tpl_library` stations in that era's core set;
  - the other candidate has a higher raw `net_cc` (ties: the earlier one wins).
- A candidate at any P can be the winner. The flag records the strongest winner: its family, `k`, time offset and distance.
- Variants (all of them are columns):
  - **A, strict (upper bound):** any station-sharing family, at any distance. Known to also catch separate events on busy
    days.
  - **B, neighbours only (≤ 10 km):** only a station-sharing family within 10 km may flag. Distance is great-circle
    between `setup/families.csv` `lat_use` / `lon_use`. One Olympic A family has no position, so it never flags and is
    never flagged under B.
  - **C = B + same arrivals:** as B, and the two detections predict the same arrival at their best shared station:
    |(t_i − t_w) − (off_w(s*) − off_i(s*))| ≤ 0.5 s.
    - s* = the shared era-core station with the largest min(template weight).
    - off = `tpl_library` `offset_s` of the family's highest-weight node at that station.
  - **D = A + same arrivals:** as A, with the arrival condition of C, at any distance.
- Olympic A position groups already merge their members' candidates at ±2 s. Here a group candidate is flagged only by a
  family **outside** its group (members never flag each other). The kept member candidates are then merged as before.
- Per-candidate flags: `flags/flags_<YYYY-MM>.parquet` (all core rows, all series; 3.6 GB).
  - Columns: `det_id` (equal to catalogue_v2 `det_id` for data rows), `dup_{A,B,C,D}`, `win_fid_*`, `win_k_*`, `win_dt_*`,
    `win_km_*`, `arrival_mismatch_*`.

## Is a flagged pair the same LFE? Chance coincidences and the arrival check
At the headline cut (data series), two checks for each variant:
- **Chance share** = flags produced when the winner pool is shifted by ±60, 120, 300 and 600 s (mean of the 8 shifts),
  divided by the true flags.
- **Arrival check** = share of flagged pairs whose best-station arrivals agree within 0.5 s, compared with the same share
  for chance pairs (shifted pool). Estimated same-event share = (pass − chance pass) / (1 − chance pass).
- This is a waveform-identity **proxy**: one station's arrival time. It is not a waveform cross-correlation. WA cut windows
  do not exist (the cut was stopped), so a true waveform check was not done.

| region   | era   | variant   |   n_selected |   n_flagged |   flagged % |   chance % of flagged |   arrival-consistent % (0.5 s) |   chance pairs consistent % |   same-event % of flagged (est.) |
|:---------|:------|:----------|-------------:|------------:|------------:|----------------------:|-------------------------------:|----------------------------:|---------------------------------:|
| ArmB     | E1    | A         |          421 |         152 |        36.1 |                  40.5 |                           59.2 |                         3.7 |                             57.7 |
| ArmB     | E1    | B         |          421 |          66 |        15.7 |                   8   |                           71.2 |                         9.5 |                             68.2 |
| ArmB     | E2    | A         |          142 |          48 |        33.8 |                  32.6 |                           54.2 |                         2.4 |                             53   |
| ArmB     | E2    | B         |          142 |          17 |        12   |                   2.9 |                           94.1 |                         0   |                             94.1 |
| ArmB     | E3    | A         |          172 |          77 |        44.8 |                  33   |                           49.4 |                         3   |                             47.8 |
| ArmB     | E3    | B         |          172 |          20 |        11.6 |                   5   |                           30   |                        25   |                              6.7 |
| NW       | E1    | A         |        12225 |        3820 |        31.2 |                  22.2 |                           87.1 |                         4.6 |                             86.5 |
| NW       | E1    | B         |        12225 |         396 |         3.2 |                  25.8 |                           83.1 |                         4.7 |                             82.3 |
| NW       | E2    | A         |         6311 |        1926 |        30.5 |                  23.3 |                           84.2 |                         4.7 |                             83.4 |
| NW       | E2    | B         |         6311 |         364 |         5.8 |                  17.4 |                           86.5 |                         5.1 |                             85.8 |
| NW       | E3    | A         |        18525 |        5723 |        30.9 |                  24.1 |                           83.2 |                         5.1 |                             82.3 |
| NW       | E3    | B         |        18525 |        1277 |         6.9 |                  14.4 |                           87.4 |                         8.6 |                             86.2 |
| OlyA     | E1    | A         |         1414 |         595 |        42.1 |                  15   |                           23.7 |                         5.1 |                             19.6 |
| OlyA     | E1    | B         |         1414 |           7 |         0.5 |                  80.4 |                           14.3 |                         8.9 |                              5.9 |
| OlyA     | E2    | A         |          939 |         489 |        52.1 |                  13.2 |                           19.4 |                         4.4 |                             15.7 |
| OlyA     | E2    | B         |          939 |           0 |         0   |                 inf   |                                |                         0   |                                  |
| OlyA     | E3    | A         |          801 |         310 |        38.7 |                  31.7 |                           29.7 |                         3.4 |                             27.2 |
| OlyA     | E3    | B         |          801 |           5 |         0.6 |                 167.5 |                            0   |                         0   |                              0   |
| SW       | E1    | A         |          976 |         138 |        14.1 |                  62.8 |                           35.5 |                         5.1 |                             32.1 |
| SW       | E1    | B         |          976 |          20 |         2   |                   8.1 |                           75   |                         7.7 |                             72.9 |
| SW       | E2    | A         |          288 |          52 |        18.1 |                  63.2 |                            7.7 |                         4.9 |                              2.9 |
| SW       | E2    | B         |          288 |           5 |         1.7 |                  15   |                           40   |                         0   |                             40   |
| SW       | E3    | A         |          482 |          96 |        19.9 |                  60.7 |                           16.7 |                         5.2 |                             12.1 |
| SW       | E3    | B         |          482 |          11 |         2.3 |                   3.4 |                           72.7 |                         0   |                             72.7 |

Variants C and D (the arrival condition is built in, so every pair passes by construction):

| region   | era   | variant   |   n_selected |   n_flagged |   flagged % |   chance % of flagged |
|:---------|:------|:----------|-------------:|------------:|------------:|----------------------:|
| ArmB     | E1    | C         |          421 |          51 |        12.1 |                   1.2 |
| ArmB     | E1    | D         |          421 |          98 |        23.3 |                   4.1 |
| ArmB     | E2    | C         |          142 |          16 |        11.3 |                   0   |
| ArmB     | E2    | D         |          142 |          30 |        21.1 |                   2.5 |
| ArmB     | E3    | C         |          172 |           7 |         4.1 |                   3.6 |
| ArmB     | E3    | D         |          172 |          43 |        25   |                   3.2 |
| NW       | E1    | C         |        12225 |         330 |         2.7 |                   1.7 |
| NW       | E1    | D         |        12225 |        3424 |        28   |                   1.6 |
| NW       | E2    | C         |         6311 |         317 |         5   |                   1.2 |
| NW       | E2    | D         |         6311 |        1670 |        26.5 |                   1.7 |
| NW       | E3    | C         |        18525 |        1118 |         6   |                   1.5 |
| NW       | E3    | D         |        18525 |        4920 |        26.6 |                   1.9 |
| OlyA     | E1    | C         |         1414 |           2 |         0.1 |                  43.8 |
| OlyA     | E1    | D         |         1414 |         180 |        12.7 |                   3   |
| OlyA     | E2    | C         |          939 |           0 |         0   |                       |
| OlyA     | E2    | D         |          939 |         119 |        12.7 |                   3.2 |
| OlyA     | E3    | C         |          801 |           0 |         0   |                       |
| OlyA     | E3    | D         |          801 |         101 |        12.6 |                   4.3 |
| SW       | E1    | C         |          976 |          16 |         1.6 |                   0.8 |
| SW       | E1    | D         |          976 |          56 |         5.7 |                  12.9 |
| SW       | E2    | C         |          288 |           2 |         0.7 |                   0   |
| SW       | E2    | D         |          288 |           7 |         2.4 |                  35.7 |
| SW       | E3    | C         |          482 |           8 |         1.7 |                   0   |
| SW       | E3    | D         |          482 |          21 |         4.4 |                  28   |

Controls at the headline cut (all regions and eras):

| series   | variant   |   n_selected |   n_flagged |   flagged % |
|:---------|:----------|-------------:|------------:|------------:|
| c1       | A         |         2112 |        1001 |        47.4 |
| c1       | B         |         2112 |          11 |         0.5 |
| c1       | C         |         2112 |          10 |         0.5 |
| c1       | D         |         2112 |         329 |        15.6 |
| rev      | A         |          478 |         148 |        31   |
| rev      | B         |          478 |           4 |         0.8 |
| rev      | C         |          478 |           1 |         0.2 |
| rev      | D         |          478 |          42 |         8.8 |

**Reading.**
- **NW: A flags 31 % of headline repeats in every era, and most of those are the same event.**
  - About 83–87 % of A-flagged pairs predict the same arrival within 0.5 s. Only about 5 % of chance pairs do, so roughly
    82–86 % of A flags are genuine double counts.
  - The time-shift test puts chance at 22–24 % of A flags. Together the two tests say A slightly over-flags.
  - D (A + same arrivals) flags 27–28 %, with chance about 2 % of D flags.
- **The double counts in NW are mostly between families 10–40 km apart, not within 10 km.**
  - D winners for NW E3: 19 % within 3 km, 2 % at 3–10 km, 32 % at 10–20 km, 33 % at 20–40 km, 13 % beyond 40 km.
  - So **B (≤ 10 km) catches only about a quarter of the same-event duplicates in NW** (B 6.9 % vs D 26.6 % in E3).
    B is the conservative "counted once" column (high same-event share, few flags). D is the best estimate of all
    same-event double counting. A is the upper bound.
- **Olympic A and SW:** the A flags are mostly not the same event (arrival-consistent share 8–35 %). There B and D change
  little.

## Headline table (net = data − time-reversed control per unit-day; headline cuts)
- **Columns:** `net_per_unit_day_*` = mean per unit-day. `regional_events_per_day_*` = sum over the region's units per day,
  averaged over days. `rate_per_unit_day_*` = k-scaled.
- **k:** for raw, the existing k. For A–D, k refitted with the same `k_factors_v2` code (imported) on de-duplicated truth
  rows and de-duplicated truth events.
  - Truth events: the same variant rule; the winner is the family whose own candidate within ±2 s has the highest cc.
  - k exists only where the production k exists (E3; E1/E2 PNET k unavailable as before).
- **The existing k must not be applied to de-duplicated counts.** It was fitted on counts that were not de-duplicated.
- 90 % day-bootstrap ranges are in `rates_summary_v2_dedup.csv`.

| region   | era   | day_class   |   net_per_unit_day_raw |   net_per_unit_day_A |   net_per_unit_day_B |   net_per_unit_day_C |   net_per_unit_day_D |   regional_events_per_day_raw |   regional_events_per_day_B |   regional_events_per_day_D |   rate_per_unit_day_raw |   rate_per_unit_day_B |   rate_per_unit_day_D |
|:---------|:------|:------------|-----------------------:|---------------------:|---------------------:|---------------------:|---------------------:|------------------------------:|----------------------------:|----------------------------:|------------------------:|----------------------:|----------------------:|
| ArmB     | E1    | quiet       |              0.001459  |            0.001064  |            0.001277  |            0.001277  |            0.001155  |                      0.0564   |                    0.04935  |                    0.04465  |                         |                       |                       |
| ArmB     | E1    | tremor      |              0.04942   |            0.03136   |            0.04176   |            0.04379   |            0.03815   |                      1.644    |                    1.389    |                    1.269    |                         |                       |                       |
| ArmB     | E1    | unknown     |              0.0008843 |            0.000485  |            0.0006846 |            0.0007132 |            0.0005991 |                      0.03626  |                    0.02807  |                    0.02456  |                         |                       |                       |
| ArmB     | E2    | quiet       |              0.0005332 |            0.0002285 |            0.0004571 |            0.0004571 |            0.0003809 |                      0.01994  |                    0.01709  |                    0.01425  |                         |                       |                       |
| ArmB     | E2    | tremor      |              0.02179   |            0.0153    |            0.01927   |            0.01945   |            0.01746   |                      0.6173   |                    0.5459   |                    0.4949   |                         |                       |                       |
| ArmB     | E3    | quiet       |              0.0001066 |            5.997e-05 |            0.0001066 |            0.0001066 |            9.329e-05 |                      0.004014 |                    0.004014 |                    0.003512 |                0.009479 |              0.004159 |              0.007607 |
| ArmB     | E3    | tremor      |              0.007506  |            0.004116  |            0.006538  |            0.007167  |            0.005521  |                      0.1965   |                    0.1711   |                    0.1445   |                0.6543   |              0.256    |              0.4529   |
| ArmB     | E3    | unknown     |              0.0002197 |            0.0002197 |            0.0002197 |            0.0002197 |            0.0002197 |                      0.009009 |                    0.009009 |                    0.009009 |                0.01978  |              0.00857  |              0.0178   |
| NW       | E1    | quiet       |              0.008309  |            0.004966  |            0.008271  |            0.008271  |            0.005284  |                      1.014    |                    1.009    |                    0.6446   |                         |                       |                       |
| NW       | E1    | tremor      |              0.4816    |            0.3305    |            0.4652    |            0.4682    |            0.3465    |                     37.55     |                   36.28     |                   27.02     |                         |                       |                       |
| NW       | E1    | unknown     |              0.01053   |            0.008054  |            0.01011   |            0.01013   |            0.008249  |                      1.389    |                    1.335    |                    1.089    |                         |                       |                       |
| NW       | E2    | quiet       |              0.00752   |            0.004694  |            0.007366  |            0.007404  |            0.005088  |                      0.8887   |                    0.8704   |                    0.6013   |                         |                       |                       |
| NW       | E2    | tremor      |              0.3514    |            0.249     |            0.3289    |            0.3316    |            0.2622    |                     19.51     |                   18.26     |                   14.56     |                         |                       |                       |
| NW       | E3    | quiet       |              0.003845  |            0.002608  |            0.003818  |            0.003828  |            0.002808  |                      0.454    |                    0.4509   |                    0.3316   |                0.03533  |              0.03478  |              0.02944  |
| NW       | E3    | tremor      |              0.282     |            0.1953    |            0.2606    |            0.2632    |            0.2072    |                     14.88     |                   13.75     |                   10.93     |                3.704    |              3.548    |              2.984    |
| NW       | E3    | unknown     |              0.02129   |            0.01474   |            0.01966   |            0.01979   |            0.01563   |                      2.811    |                    2.595    |                    2.063    |                0.3374   |              0.3211   |              0.2631   |
| OlyA     | E1    | quiet       |              0.03575   |            0.0166    |            0.03575   |            0.03575   |            0.03002   |                      0.3337   |                    0.3337   |                    0.2802   |                2.525    |              2.383    |              2.728    |
| OlyA     | E1    | tremor      |              0.2726    |            0.1967    |            0.268     |            0.2713    |            0.2455    |                      1.626    |                    1.598    |                    1.465    |               13.12     |             12.2      |             14.93     |
| OlyA     | E1    | unknown     |              0.03953   |            0.01333   |            0.03953   |            0.03953   |            0.03193   |                      0.3953   |                    0.3953   |                    0.3193   |                2.592    |              2.446    |              2.689    |
| OlyA     | E2    | quiet       |              0.06706   |            0.02855   |            0.06706   |            0.06706   |            0.05751   |                      0.6086   |                    0.6086   |                    0.5219   |                5.983    |              5.69     |              5.49     |
| OlyA     | E2    | tremor      |              0.2505    |            0.1578    |            0.2505    |            0.2505    |            0.2318    |                      1.05     |                    1.05     |                    0.9711   |               11.1      |             10.52     |             10.63     |
| OlyA     | E3    | quiet       |              0.006785  |            0.002491  |            0.006785  |            0.006785  |            0.005539  |                      0.06146  |                    0.06146  |                    0.05018  |                1.928    |              1.869    |              1.954    |
| OlyA     | E3    | tremor      |              0.1304    |            0.09744   |            0.1292    |            0.1304    |            0.1179    |                      0.509    |                    0.504    |                    0.4601   |               29.94     |             28.62     |             35.04     |
| OlyA     | E3    | unknown     |              0.00991   |            0.006306  |            0.00991   |            0.00991   |            0.009009  |                      0.0991   |                    0.0991   |                    0.09009  |                2.533    |              2.45     |              2.945    |
| SW       | E1    | tremor      |              0.07689   |            0.06596   |            0.07531   |            0.07563   |            0.07246   |                      3.735    |                    3.658    |                    3.519    |                         |                       |                       |
| SW       | E2    | tremor      |              0.03676   |            0.03037   |            0.03611   |            0.0365    |            0.03585   |                      1.226    |                    1.204    |                    1.196    |                         |                       |                       |
| SW       | E3    | tremor      |              0.0166    |            0.01326   |            0.0162    |            0.01631   |            0.01587   |                      0.5607   |                    0.5472   |                    0.5362   |                1.744    |              1.669    |              1.684    |

**What changes, in plain words.**
- **NW (headline region):** counting each event once lowers repeats per family-day in E3 tremor days from 0.282 to:
  - 0.261 with B (−8 %);
  - 0.207 with D (−27 %);
  - 0.195 with A (−31 %).
  - Regional events per day: 14.9 raw, 13.7 B, 10.9 D.
  - The same holds in E1/E2 and on quiet / unknown days, so ratios between eras and day classes barely move.
- **k-scaled NW E3 tremor rate (k refitted):** 3.70 raw, 3.55 B, 2.98 D, 2.70 A. The truth catalogue is double counted as
  well (A flags 55,231 of 198,804 truth events; D 36,384), which offsets part of the drop.
- **Arm B:** −13 % with B, −26 % with D. **SW:** −2 % with B, −4 % with D. **Olympic A groups:** 0 % with B, −10 % with D.
- **Nothing replaces the headline.** The raw columns equal `rates_daily_v2` exactly (`acceptance_dedup.json`: 0 differing
  values at every cut and series, including the Olympic A merged counts). The raw rows of `rates_summary_v2_dedup.csv`
  equal `rates_summary_v2.csv` (max difference 9e-17).

## Files
- `rates_daily_v2_dedup.parquet`: 1,842,484 rows = every headline row of `rates_daily_v2`, zero rows included.
  - Units: family own_pnet for NW / SW / Arm B; Olympic A group nbr3.
  - Columns: `n_v2_{data,rev,c1}_{A,B,C,D}`, `net_count_v2_*`, `n0xx_*` at every cut.
- `rates_summary_v2_dedup.csv`: the `rates_summary_v2` layout plus `count_variant` (raw / A / B / C / D) and the
  regional events per day.
- `headline_table_dedup.csv`.
- `chance_identity_dedup.csv`: per series × era × region × day class × variant.
- `rates_k_factors_v2_dedup_{A,B,C,D}.csv`.
- `flags/`, `unit_day_counts_dedup.parquet`, `acceptance_dedup.json`, `summary_check_vs_rates_summary_v2.json`, `dedup.log`.
