## Step 2 (W2) — station health of the 21 core stations, 2007 → 2026-09-15 (auto-written by scripts/pre_report.py)

- Raw metrics: 87112 station-days read by this run (72613 ok, 5 errors, 14494 expected but no file); the 13 SVI-shared stations 2014-01-01 → 2026-09-15 come from SVI R2's daily table (same files, same code). Files: `w2/r2_health_daily.parquet`, `w2/r2_health_stationday.parquet`.
- Stack test: 261 strongest-tremor days (2007–2026; sealed 2007-01 / 2009-05 and 2011-04 never used), cache {'ok': 3513, 'wa_cache': 1318, 'no_raw': 409, 'cached': 230, 'raw_not_prod_visible': 10, 'empty': 1}; 1111 station-months tested (≥ 3 families with ≥ 5 confident repeats, ≥ 30 repeats). Confident repeats for a station = detections whose network WITHOUT that station's site still reaches 8.75 × MAD with ≥ 2 other sites at site z ≥ 3. Scrambled control = the same station read on other selected days at the same time of day.
- Station-month status (all 4859 metadata-alive station-months): {'inherit_ok': 2508, 'inherit_control_level': 808, 'ok': 684, 'control_level': 426, 'nodata': 418, 'dead': 13, 'untested': 1, 'inherit_dead': 1}.
- Tested months at control level: 426 of 1111. Excluded station-months (health reason, not plain no-data): 1248 (`w2/r2_exclusions.csv`).

| station | tested months | at control level | median stack cc at placement | median control p95 |
|---|---|---|---|---|
| CN.PGC | 56 | 26 | 0.088 | 0.088 |
| CN.VGZ | 54 | 21 | 0.105 | 0.101 |
| PB.B001 | 59 | 9 | 0.110 | 0.092 |
| PB.B003 | 58 | 16 | 0.098 | 0.092 |
| PB.B004 | 57 | 25 | 0.088 | 0.086 |
| PB.B005 | 60 | 17 | 0.095 | 0.086 |
| PB.B006 | 60 | 19 | 0.091 | 0.083 |
| PB.B007 | 60 | 17 | 0.099 | 0.087 |
| PB.B009 | 53 | 19 | 0.085 | 0.082 |
| PB.B010 | 53 | 21 | 0.088 | 0.082 |
| PB.B011 | 53 | 21 | 0.087 | 0.083 |
| PB.B013 | 23 | 0 | 0.535 | 0.087 |
| PB.B014 | 53 | 26 | 0.086 | 0.086 |
| PB.B017 | 59 | 39 | 0.088 | 0.090 |
| PB.B018 | 58 | 33 | 0.084 | 0.086 |
| PB.B926 | 51 | 30 | 0.078 | 0.079 |
| PB.B943 | 23 | 0 | 0.591 | 0.088 |
| UW.DOSE | 50 | 24 | 0.100 | 0.092 |
| UW.GNW | 60 | 15 | 0.102 | 0.085 |
| UW.HDW | 60 | 25 | 0.155 | 0.146 |
| UW.LRIV | 51 | 23 | 0.091 | 0.088 |

Swap / epoch-change checks (`w2/r2_swap_checks.csv`; carry-over = the old node's template used on the new node):

| station | date | change | cc before / ctl p95 | cc after / ctl p95 | best channel transform after | verdict |
|---|---|---|---|---|---|---|
| CN.PGC | 2015-01-12 | node change (WA_PROD_PLAN W2 / node tables) | 0.136 / 0.099 | 0.074 / 0.088 | (0, 1, 2)(1, -1, -1);(1, 2, 0)(1, 1, 1) | AT CONTROL LEVEL |
| CN.PGC | 2017-08-10 | channels  -> HHE,HHN,HHZ; sensor  -> Nanometrics Trillium 120 Sec PH R | 0.088 / 0.088 | 0.097 / 0.082 | (1, 0, 2)(-1, -1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| CN.VGZ | 2007-02-28 | sensor Guralp CMG3ESP/CNSN GD1 digitizer -> Guralp CMG3T_30sec/CNSN GD | nan / nan | 0.376 / 0.092 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| CN.VGZ | 2009-07-31 | channels  -> HHE,HHN,HHZ; sensor  -> Guralp CMG3ESP_60sec/CNSN GD2_v21 | 0.270 / 0.093 | 0.314 / 0.094 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| CN.VGZ | 2017-03-03 | sensor Guralp CMG3ESP_60sec/CNSN GD2_v21 digitizer -> Nanometrics Tril | 0.135 / 0.109 | 0.108 / 0.119 | (0, 2, 1)(1, -1, -1);(0, 1, 2)(1, 1, 1) | AT CONTROL LEVEL |
| PB.B013 | 2007-01-06 | channels  -> EH1,EH2,EHZ; sensor  -> HS-1-LT/Quanterra 330 Linear Phas | nan / nan | 0.710 / 0.092 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| PB.B014 | 2008-03-05 | channels  -> EH1,EH2,EHZ; sensor  -> HS-1-LT/Quanterra 330 Linear Phas | nan / nan | 0.112 / 0.075 | (0, 1, 2)(1, 1, 1);(0, 2, 1)(1, 1, 1) | above control (carry-over OK) |
| PB.B017 | 2007-07-12 | channels  -> EH1,EH2,EHZ; sensor  -> HS-1-LT/Quanterra 330 Linear Phas | nan / nan | 0.136 / 0.089 | (0, 1, 2)(-1, -1, -1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| PB.B926 | 2007-11-13 | channels  -> EH1,EH2,EHZ; sensor  ->  | nan / nan | 0.236 / 0.080 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| PB.B943 | 2008-02-27 | channels  -> EH1,EH2,EHZ; sensor  -> HS-1-LT/Quanterra 330 Linear Phas | nan / nan | 0.721 / 0.085 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| UW.DOSE | 2010-07-21 | channels  -> BHE,BHN,BHZ; sensor  -> TRILLIUM 120PA,Velocity Sensor,NA | nan / nan | 0.176 / 0.085 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| UW.DOSE | 2019-05-07 | channels BHE,BHN,BHZ -> HHE,HHN,HHZ | 0.098 / 0.092 | 0.087 / 0.089 | (0, 1, 2)(-1, 1, 1);(0, 2, 1)(-1, 1, 1) | AT CONTROL LEVEL |
| UW.GNW | 2008-04-18 | gain BHE x0.533;BHN x0.533;BHZ x0.533 | 0.102 / 0.086 | 0.195 / 0.083 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| UW.GNW | 2010-01-29 | gain BHE x1.878;BHN x1.878;BHZ x1.878 | 0.095 / 0.086 | 0.391 / 0.083 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| UW.GNW | 2010-12-01 | channels  -> BHE,BHN,BHZ; sensor  -> Velocity Transducer | 0.284 / 0.084 | 0.193 / 0.084 | (0, 1, 2)(1, 1, 1);(2, 0, 1)(-1, 1, 1) | above control (carry-over OK) |
| UW.GNW | 2019-05-07 | channels BHE,BHN,BHZ -> HHE,HHN,HHZ | 0.143 / 0.079 | 0.123 / 0.079 | (0, 1, 2)(1, 1, 1);(1, 2, 0)(1, -1, 1) | above control (carry-over OK) |
| UW.HDW | 2008-04-01 | gain EHZ x0.505 | 0.183 / 0.149 | 0.216 / 0.147 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| UW.HDW | 2008-10-01 | channels EHZ -> EH2; sensor L-4C VERTICAL,Velocity Sensor,SERCEL -> S- | 0.229 / 0.148 | 0.162 / nan | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | untested |
| UW.HDW | 2020-12-04 | channels EH2,EHZ -> EHZ; sensor L-4C VERTICAL,Velocity Sensor,SERCEL;S | 0.150 / 0.142 | 0.129 / 0.145 | (0, 1, 2)(-1, 1, 1);(0, 1, 2)(-1, 1, 1) | AT CONTROL LEVEL |
| UW.HDW | 2022-09-01 | channels EH2,EHZ -> HHE,HHN,HHZ; sensor L-4C VERTICAL,Velocity Sensor, | 0.135 / 0.146 | 0.165 / 0.132 | (0, 1, 2)(-1, 1, 1);(1, 0, 2)(1, 1, 1) | above control (carry-over OK) |
| UW.LRIV | 2009-08-27 | channels  -> BHE,BHN,BHZ; sensor  -> Velocity Transducer | nan / nan | 0.292 / 0.096 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(1, 1, 1) | above control (carry-over OK) |
| UW.LRIV | 2019-05-30 | channels BHE,BHN,BHZ -> HHE,HHN,HHZ | 0.111 / 0.109 | 0.082 / 0.079 | (0, 1, 2)(1, 1, 1);(0, 1, 2)(-1, -1, -1) | above control (carry-over OK) |
| PB.B004 | 2020 | gap year check | nan / nan | 0.090 / 0.080 | nan | above control (carry-over OK) |
| PB.B014 | 2023 | gap year check | nan / nan | 0.080 / 0.081 | nan | AT CONTROL LEVEL |
| PB.B014 | 2024 | gap year check | nan / nan | nan / nan | nan | untested |
| PB.B926 | 2025 | gap year check | nan / nan | 0.068 / 0.074 | nan | AT CONTROL LEVEL |

- Positive control (LRIV after 2019-05-30, HH node, old BH template carried): tested months 27, at control level 14, median cc 0.085 vs control p95 0.087; before 2019-05: 24 tested, 9 at control, median cc 0.100. Flagged: True. (SVI R2: LRIV 2021-10 0.17 vs control p95 0.09 (weak, above control); A-2 -0.216 on its own event set)
- Effective station set per month: `w2/r2_effective_sets.csv` (median 14 of 21, min 2 in 2022-02).
- Bad-for-dv/v windows (KMR style, `w2/r2_bad_for_dvv.csv`): 16779 windows; by reason {'line': 12075, 'noisy': 1455, 'month': 1248, 'flat': 625, 'dc_jump': 526, 'gap': 501, 'dead': 171, '2-8': 165, 'meta_change': 13}.
- Read with care: the stack test has little power where a station's templates are weak for the families active on the test days (e.g. SW carriers B017/B018 on NW-dominated days); a 'control_level' month there means 'no evidence the station sees these families', not 'broken'. Months between test days inherit the nearest tested month of the same instrument epoch.
- W0 re-hash at the end: 139 frozen input files, changed 2 (`w0/verify_20261003T0252.csv`).
- W1 final readability (full reads by W2 raw + SVI R2): of 131597 present raw station-days, 131590 read in full without error, 7 errors, 0 not fully read by either (`w1/w1_readability_final.csv`).
