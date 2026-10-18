# WA permanent-station trial extension — results

Written 2026-09-29 17:14Z by `scripts/wa_report.py`. Design: `PLAN.md`; log: `NOTES.md`.

## Data

- Downloaded: 2134 station-days, 26.3 GB; nodata 91, errors 0.
- 20 Hz cache: 4575 station-days (42 stations).

## Permanent templates (build half of the array-era detections)

- (family, station-epoch) stacks with >= 1 event: 13830; with >= 30 events: 11358; clear arrival (imp > 5, dur < 3 s at 2–8 Hz, fs = 20 Hz): **386**.
- Families with >= 1 clear permanent template: 180 of 329; with >= 2 clear SITES (the observability rule): **61**.

| region | families | ≥1 clear template | ≥2 clear sites |
|---|---|---|---|
| Olympic_2010 | 75 | 21 | 2 |
| Olympic_2010_armB | 43 | 19 | 2 |
| Royer_NW_WA | 132 | 115 | 55 |
| Royer_SW_WA | 79 | 25 | 2 |

Clear templates per station: PB.B943 104, PB.B013 102, UW.BLN 45, UW.HDW 39, PB.B001 35, PB.B018 13, US.NLWA 11, UW.GNW 7, CN.PGC 6, CN.VGZ 6, UW.LRIV 4, PB.B005 2, UW.SMW 2, UW.OSD 2, PB.B006 2, PB.B007 2, PB.B009 1, PB.B010 1, UW.DOSE 1, UW.OTR 1

## HF veto cuts (99th percentile of build-half own-matched ratios, per HF station node)

| HF node | n calib | cut | median matched |
|---|---|---|---|
| PB.B943..EH.100@2008-02-27 | 3332 | 6.319 | 0.090 |
| PB.B013..EH.100@2007-01-06 | 2373 | 25.037 | 0.056 |
| PB.B001..EH.100@2005-06-30 | 304 | 0.404 | 0.031 |
| PB.B007..EH.100@2005-07-22 | 234 | 1.476 | 0.035 |
| US.NLWA..BH.40@2006-07-05 | 29 | – | 0.028 |
| PB.B018..EH.100@2006-01-21 | 27 | – | 0.054 |
| UW.HDW..EH.100@1998-05-21 | 23 | – | 0.064 |
| UW.BLN..EH.100@1995-01-19 | 5 | – | 0.099 |

## Validation — the known answer (array-era detections of the same families)

Test half only (2-h blocks, 60-s guards), observable family-days (≥ 2 permanent sites) the family's own array run scanned. One-to-one ±2 s, own family. `cat` = tier & HF pass & cross-family kept; `raw` = tier only. Scrambled control = site-rolled series (each site shifted by its own 300–3000 s), detections per family-day at the same tier (HF veto / cross-family rule not applicable to it).

`control prec` = fraction of scrambled-control detections (test half) that land within ±2 s of an own-family array detection by chance (the precision baseline); `control recall` = array detections matched by the control.

| window | region | fam rows | fam obs | truth (test) | recall core | prec own core | prec any core | control recall / prec core | recall ext | prec own ext | data/famday core | control/famday core | quiet data / control core |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ALL | ALL | 461 | 116 | 16917 | 0.158 | 0.250 | 0.439 | 0.006 / 0.007 | 0.188 | 0.202 | 13.948 | 12.237 | 11.693 / 11.569 |
| V2008 | ALL | 211 | 57 | 8059 | 0.153 | 0.347 | 0.505 | 0.004 / 0.007 | 0.183 | 0.260 | 10.809 | 9.304 | 6.661 / 6.503 |
| V2010 | ALL | 250 | 59 | 8858 | 0.162 | 0.202 | 0.407 | 0.008 / 0.008 | 0.192 | 0.170 | 17.328 | 15.395 | 13.151 / 13.037 |
| V2008 | Royer_NW_WA | 132 | 55 | 7695 | 0.157 | 0.343 | 0.502 | 0.004 / 0.007 | 0.188 | 0.257 | 11.082 | 9.587 | 6.903 / 6.739 |
| V2008 | Royer_SW_WA | 79 | 2 | 364 | 0.069 | 0.714 | 0.743 | 0.000 / 0.000 | 0.080 | 0.659 | 2.175 | 0.350 | 0.000 / 0.000 |
| V2010 | Olympic_2010 | 75 | 2 | 384 | 0.432 | 0.053 | 0.312 | 0.005 / 0.001 | 0.443 | 0.047 | 120.446 | 119.161 | 137.400 / 137.450 |
| V2010 | Olympic_2010_armB | 43 | 2 | 268 | 0.045 | 0.571 | 0.762 | 0.000 / 0.000 | 0.056 | 0.395 | 0.875 | 0.393 | 0.150 / 0.150 |
| V2010 | Royer_NW_WA | 132 | 55 | 8206 | 0.154 | 0.315 | 0.479 | 0.009 / 0.012 | 0.184 | 0.238 | 12.915 | 10.876 | 9.105 / 8.982 |

**Per-network-epoch tiers** (SVI `c_tiers_by_epoch.py` rule, FA only: T_core = smallest T ≥ 8.0 on a 0.25 grid with quiet-day FA ≤ 0.1 per template-day, T_ext ≤ 1.0; quiet = the window's non-episode observable family-days; V2008's are post-episode days with no tremor catalogue):

| window | quiet fam-days | FA at 8.25 | FA at 8.75 | T_ext_ep | T_core_ep |
|---|---|---|---|---|---|
| V2008 | 183 | 12.104 | 8.508 | 16.0 | 20.75 |
| V2010 | 610 | 19.110 | 12.454 | 13.5 | 17.75 |

| window | region | recall core_ep | prec own core_ep | prec any core_ep | control recall / prec core_ep | recall ext_ep | prec own ext_ep | data/famday core_ep | control/famday core_ep |
|---|---|---|---|---|---|---|---|---|---|
| ALL | ALL | 0.003 | 0.389 | 0.850 | 0.000 / 0.000 | 0.017 | 0.314 | 0.103 | 0.077 |
| V2008 | ALL | 0.000 | 1.000 | 1.000 | 0.000 / 0.000 | 0.006 | 0.944 | 0.010 | 0.005 |
| V2010 | ALL | 0.005 | 0.367 | 0.844 | 0.000 / 0.000 | 0.027 | 0.275 | 0.204 | 0.155 |
| V2008 | Royer_NW_WA | 0.001 | 1.000 | 1.000 | 0.000 / 0.000 | 0.007 | 0.944 | 0.010 | 0.006 |
| V2008 | Royer_SW_WA | 0.000 | – | – | 0.000 / – | 0.000 | – | 0.000 | 0.000 |
| V2010 | Olympic_2010 | 0.029 | 0.145 | 0.816 | 0.000 / 0.000 | 0.174 | 0.103 | 2.946 | 2.893 |
| V2010 | Olympic_2010_armB | 0.000 | 0.000 | 1.000 | 0.000 / – | 0.004 | 0.500 | 0.018 | 0.000 |
| V2010 | Royer_NW_WA | 0.004 | 0.906 | 0.906 | 0.000 / 0.000 | 0.021 | 0.795 | 0.074 | 0.024 |

Per family: `results/score_val/validation_by_family.csv`. Families with ≥ 20 test-half array detections on observable days: 116; recall(core) median 0.136, ≥ 0.2: 27, ≥ 0.5: 0.

Top 15 by recall (core):

| family | window | med sites | truth test | recall core | prec own core | data/famday | control/famday |
|---|---|---|---|---|---|---|---|
| RNW:mute:2435 | V2010 | 4 | 140 | 0.457 | 0.552 | 13.000 | 7.750 |
| OLY:938 | V2010 | 2 | 340 | 0.447 | 0.049 | 237.750 | 238.036 |
| RNW:mute:2435 | V2008 | 4 | 143 | 0.434 | 0.602 | 8.783 | 5.000 |
| RNW:fair:542 | V2008 | 2 | 166 | 0.325 | 0.211 | 46.652 | 45.391 |
| RNW:mute:3288 | V2008 | 2 | 93 | 0.323 | 0.157 | 21.652 | 22.261 |
| OLY:127 | V2010 | 3 | 44 | 0.318 | 0.424 | 3.143 | 0.286 |
| RNW:fair:507 | V2008 | 2 | 114 | 0.316 | 0.419 | 9.609 | 5.478 |
| RNW:lag15:1236 | V2010 | 2 | 87 | 0.310 | 0.435 | 9.050 | 3.450 |
| RNW:run1:2800 | V2010 | 4 | 183 | 0.306 | 0.258 | 24.550 | 21.500 |
| RNW:fair:542 | V2010 | 2 | 173 | 0.289 | 0.181 | 51.900 | 47.350 |
| RNW:lag15:3076 | V2010 | 3 | 140 | 0.286 | 0.513 | 8.900 | 7.350 |
| RNW:mute:3288 | V2010 | 2 | 86 | 0.267 | 0.111 | 22.300 | 22.500 |
| RNW:lag15:2226 | V2008 | 2 | 169 | 0.266 | 0.243 | 44.435 | 42.391 |
| RNW:lag15:1352 | V2010 | 5 | 154 | 0.266 | 0.872 | 6.950 | 1.550 |
| RNW:qc:2653 | V2008 | 2 | 245 | 0.257 | 0.223 | 43.783 | 43.261 |

## Trial — 2012-09 ETS, after the arrays (XU ended 2012-05-31, XG 2011-09)

Window 2012-08-20..24 (quiet), 2012-09-05..10-13 (episode), 2012-10-18..22 (quiet). No truth: tremor days (PNSN, the region's tremor-day table, ≥ 1 window) are a LABEL, never a gate.

New station epochs re-templated from pass-1 detections: US.NLWA.00.BH.40@2011-05-02, UW.GNW..BH.40@2010-12-01, UW.RER..EH.100@2012-09-06, UW.SP2..BH.40@2011-08-25, UW.SQM..EH.100@2011-04-22, UW.WISH..BH.40@2011-06-14.
  → 348 stacks, 0 clear templates.

| pass | region | families | with ≥2 live permanent sites | cat core | on tremor days | other days | rate/famday tremor | rate/famday other | control rate tremor | control rate other |
|---|---|---|---|---|---|---|---|---|---|---|
| pass 1 (validation-epoch templates only) [core] | ALL | 329 | 59 | 16068 | 11642 | 4426 | 6.283 | 4.264 | 7.189 | 5.249 |
| pass 1 (validation-epoch templates only) [core] | Olympic_2010 | 75 | 2 | 9693 | 6906 | 2787 | 93.324 | 116.125 | 100.041 | 126.958 |
| pass 1 (validation-epoch templates only) [core] | Olympic_2010_armB | 43 | 2 | 30 | 24 | 6 | 0.324 | 0.250 | 0.176 | 0.167 |
| pass 1 (validation-epoch templates only) [core] | Royer_NW_WA | 132 | 55 | 6345 | 4712 | 1633 | 2.764 | 1.649 | 3.464 | 2.421 |
| pass 1 (validation-epoch templates only) [core] | Royer_SW_WA | 79 | 0 | 0 | 0 | 0 | – | – | – | – |
| pass 2 (+ re-templated new epochs) [core] | ALL | 329 | 59 | 16068 | 11642 | 4426 | 6.283 | 4.264 | 7.189 | 5.249 |
| pass 2 (+ re-templated new epochs) [core] | Olympic_2010 | 75 | 2 | 9693 | 6906 | 2787 | 93.324 | 116.125 | 100.041 | 126.958 |
| pass 2 (+ re-templated new epochs) [core] | Olympic_2010_armB | 43 | 2 | 30 | 24 | 6 | 0.324 | 0.250 | 0.176 | 0.167 |
| pass 2 (+ re-templated new epochs) [core] | Royer_NW_WA | 132 | 55 | 6345 | 4712 | 1633 | 2.764 | 1.649 | 3.464 | 2.421 |
| pass 2 (+ re-templated new epochs) [core] | Royer_SW_WA | 79 | 0 | 0 | 0 | 0 | – | – | – | – |
| pass 1 (validation-epoch templates only) [core_ep] | ALL | 329 | 59 | 264 | 196 | 68 | 0.106 | 0.066 | 0.100 | 0.063 |
| pass 1 (validation-epoch templates only) [core_ep] | Olympic_2010 | 75 | 2 | 251 | 183 | 68 | 2.473 | 2.833 | 2.514 | 2.708 |
| pass 1 (validation-epoch templates only) [core_ep] | Olympic_2010_armB | 43 | 2 | 0 | 0 | 0 | 0.000 | 0.000 | 0.000 | 0.000 |
| pass 1 (validation-epoch templates only) [core_ep] | Royer_NW_WA | 132 | 55 | 13 | 13 | 0 | 0.008 | 0.000 | 0.000 | 0.000 |
| pass 1 (validation-epoch templates only) [core_ep] | Royer_SW_WA | 79 | 0 | 0 | 0 | 0 | – | – | – | – |
| pass 2 (+ re-templated new epochs) [core_ep] | ALL | 329 | 59 | 264 | 196 | 68 | 0.106 | 0.066 | 0.100 | 0.063 |
| pass 2 (+ re-templated new epochs) [core_ep] | Olympic_2010 | 75 | 2 | 251 | 183 | 68 | 2.473 | 2.833 | 2.514 | 2.708 |
| pass 2 (+ re-templated new epochs) [core_ep] | Olympic_2010_armB | 43 | 2 | 0 | 0 | 0 | 0.000 | 0.000 | 0.000 | 0.000 |
| pass 2 (+ re-templated new epochs) [core_ep] | Royer_NW_WA | 132 | 55 | 13 | 13 | 0 | 0.008 | 0.000 | 0.000 | 0.000 |
| pass 2 (+ re-templated new epochs) [core_ep] | Royer_SW_WA | 79 | 0 | 0 | 0 | 0 | – | – | – | – |

t2012p1 epoch tiers (trial quiet days, 590 family-days): FA at 8.75 = 5.793, T_ext_ep 13.0, T_core_ep 18.25.

t2012p2 epoch tiers (trial quiet days, 590 family-days): FA at 8.75 = 5.793, T_ext_ep 13.0, T_core_ep 18.25.

Per family: `results/score_t2012p2/trial_by_family.csv`. Top 20 by catalogue count (core):

| family | region | obs days | med sites | cat core | tremor days / rate | other days / rate | control rate tremor / other |
|---|---|---|---|---|---|---|---|
| OLY:938 | Olympic_2010 | 49 | 2 | 9669 | 37 / 186.000 | 12 / 232.250 | 199.973 / 253.917 |
| RNW:lag15:1782 | Royer_NW_WA | 49 | 2 | 653 | 31 / 10.548 | 18 / 18.111 | 11.452 / 20.611 |
| RNW:run1:2800 | Royer_NW_WA | 49 | 3 | 594 | 31 / 9.387 | 18 / 16.833 | 10.161 / 18.333 |
| RNW:qc:1479 | Royer_NW_WA | 49 | 4 | 396 | 31 / 6.613 | 18 / 10.611 | 9.161 / 14.222 |
| RNW:fair:542 | Royer_NW_WA | 49 | 2 | 322 | 31 / 9.323 | 18 / 1.833 | 12.677 / 3.222 |
| RNW:lag15:2224 | Royer_NW_WA | 49 | 2 | 245 | 31 / 6.871 | 18 / 1.778 | 12.194 / 3.500 |
| RNW:mute:2603 | Royer_NW_WA | 49 | 3 | 237 | 31 / 7.065 | 18 / 1.000 | 10.581 / 1.889 |
| RNW:mute:2435 | Royer_NW_WA | 49 | 4 | 228 | 31 / 7.161 | 18 / 0.333 | 5.032 / 0.389 |
| RNW:run1:2671 | Royer_NW_WA | 49 | 2 | 221 | 31 / 5.903 | 18 / 2.111 | 9.871 / 3.611 |
| RNW:lag15:2226 | Royer_NW_WA | 49 | 2 | 213 | 31 / 6.323 | 18 / 0.944 | 11.387 / 2.111 |
| RNW:qc:2653 | Royer_NW_WA | 49 | 2 | 204 | 31 / 5.323 | 18 / 2.167 | 12.258 / 4.389 |
| RNW:lag15:3076 | Royer_NW_WA | 49 | 3 | 178 | 31 / 5.613 | 18 / 0.222 | 5.742 / 0.333 |
| RNW:mute:2607 | Royer_NW_WA | 49 | 2 | 146 | 31 / 4.226 | 18 / 0.833 | 9.484 / 1.889 |
| RNW:lag15:553 | Royer_NW_WA | 49 | 4 | 138 | 31 / 2.903 | 18 / 2.667 | 2.968 / 3.889 |
| RNW:mute:3288 | Royer_NW_WA | 49 | 2 | 136 | 31 / 3.097 | 18 / 2.222 | 3.452 / 2.778 |
| RNW:lag15:1236 | Royer_NW_WA | 49 | 2 | 128 | 31 / 3.677 | 18 / 0.778 | 1.581 / 0.722 |
| RNW:mute:1608 | Royer_NW_WA | 49 | 2 | 128 | 31 / 2.645 | 18 / 2.556 | 3.935 / 4.000 |
| RNW:mute:1144 | Royer_NW_WA | 49 | 3 | 120 | 31 / 2.742 | 18 / 1.944 | 2.742 / 2.556 |
| RNW:run1:2341 | Royer_NW_WA | 49 | 2 | 117 | 31 / 2.742 | 18 / 1.778 | 4.419 / 3.833 |
| RNW:lag15:641 | Royer_NW_WA | 49 | 2 | 113 | 31 / 2.871 | 18 / 1.333 | 2.968 / 1.611 |

