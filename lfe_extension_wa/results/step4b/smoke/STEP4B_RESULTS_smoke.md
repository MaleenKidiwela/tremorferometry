# Step 4b — added stations for the NW and Olympic families

Written 2026-09-29 19:28Z by `scripts/step4b_report.py`. Plan: `PLAN.md` §Step 4b; log `NOTES.md`.

## Stations added

- 52 stations ({'7D': 5, 'C8': 2, 'CC': 1, 'CN': 6, 'PB': 1, 'TA': 5, 'UW': 31, 'XD': 1}); tier 1 = permanent / served later arrays (42), tier 2 = TA / 7D taken only where they are the only option (10).
- Station-days in the step-4b windows: 15401 live; downloaded (this step, incl. the 42 trial stations in 2013/2014/A2011): 1825 station-days, 21.6 GB; nodata 95.
- Candidate table (every candidate, every column, exclusion reason): `results/step4b/candidates.csv`.

| station | tier | km to nearest family | families ≤30 km | comps | live windows | clear templates (families) |
|---|---|---|---|---|---|---|
| 7D.FN05C | 2 | 2.8 | 1 | 3 | T2013 | 0 |
| UW.ALKI | 1 | 3.9 | 16 | 1 | V2010,A2011,T2012,T2013,T2014 | 0 |
| TA.D03D | 2 | 6.5 | 48 | 3 | V2010,A2011,T2012,T2013,T2014 | 0 |
| 7D.FN07C | 2 | 7.9 | 1 | 3 | T2013 | 0 |
| UW.UWFH | 1 | 10.2 | 5 | 1 | V2010,T2012,T2013,T2014 | 0 |
| UW.OBC | 1 | 10.3 | 6 | 1 | V2010,A2011,T2012,T2013,T2014 | 0 |
| UW.MEAN | 1 | 11.2 | 13 | 1 | V2010,T2012,T2013,T2014 | 0 |
| UW.OOW | 1 | 11.8 | 3 | 1 | V2010,A2011,T2012,T2013,T2014 | 0 |
| UW.OFR | 1 | 12.0 | 2 | 3 | V2010,A2011,T2012 | 0 |
| CN.PFB | 1 | 15.8 | 1 | 3 | V2010,A2011,T2013,T2014 | 0 |
| UW.MEW | 1 | 16.6 | 4 | 3 | V2010,A2011,T2012,T2013 | 0 |
| C8.SOKB | 1 | 16.6 | 6 | 3 | T2013,T2014 | 0 |
| CN.CLRS | 1 | 17.0 | 1 | 3 | T2014 | 0 |
| PB.B926 | 1 | 17.0 | 1 | 3 | V2010,A2011,T2012,T2013,T2014 | 0 |
| CN.SNB | 1 | 17.8 | 2 | 3 | V2010,A2011,T2012,T2013,T2014 | 0 |
| CN.YOUB | 1 | 18.2 | 1 | 3 | V2010,A2011,T2012,T2013 | 0 |
| CN.LZB | 1 | 18.3 | 1 | 3 | V2010,A2011,T2012,T2013 | 0 |
| TA.D04E | 2 | 18.5 | 5 | 3 | T2012,T2013,T2014 | 0 |
| TA.D04D | 2 | 18.5 | 5 | 3 | V2010,A2011,T2012 | 0 |
| UW.MCW | 1 | 18.9 | 1 | 1 | V2010,A2011,T2012,T2013,T2014 | 0 |
| 7D.FN18C | 2 | 19.4 | 1 | 3 | T2013 | 0 |
| UW.OSR | 1 | 21.3 | 2 | 3 | V2010,A2011 | 0 |
| UW.OCP | 1 | 22.5 | 1 | 1 | V2010,A2011,T2012,T2013,T2014 | 0 |
| UW.FORK | 1 | 23.0 | 2 | 3 | V2010,A2011,T2012,T2013,T2014 | 0 |
| C8.BPCB | 1 | 23.9 | 1 | 3 | V2010,A2011,T2012,T2013 | 0 |
| 7D.FN01C | 2 | 26.7 | 1 | 3 | T2013 | 0 |
| UW.TOLT | 1 | 26.8 | 1 | 3 | V2010,A2011,T2012,T2013,T2014 | 0 |
| 7D.FN14C | 2 | 27.2 | 1 | 3 | T2013 | 0 |
| UW.ELW | 1 | 27.5 | 1 | 3 | V2010,A2011,T2012 | 0 |
| TA.A04D | 2 | 29.2 | 1 | 3 | V2010,A2011,T2012,T2013,T2014 | 0 |
| UW.GHW | 1 | 29.8 | 1 | 1 | V2010,A2011,T2012,T2013,T2014 | 0 |
| UW.RMW | 1 | 33.3 | 0 | 1 | V2010 | 2 |
| UW.RATT | 1 | 36.7 | 0 | 3 | A2011,T2012,T2013,T2014 | 0 |
| CN.GOBB | 1 | 37.7 | 0 | 3 | V2010,A2011,T2012,T2013,T2014 | 0 |
| TA.B05D | 2 | 39.0 | 0 | 3 | V2010,A2011,T2012,T2013,T2014 | 0 |
| UW.ATES | 1 | 40.5 | 0 | 1 | V2010 | 0 |
| UW.CMW | 1 | 43.9 | 0 | 1 | V2010,A2011,T2012,T2013,T2014 | 0 |
| UW.PCMD | 1 | 45.2 | 0 | 1 | V2010,T2012 | 0 |
| UW.ON2 | 1 | 47.9 | 0 | 1 | V2010,A2011,T2012,T2013,T2014 | 0 |
| UW.JCW | 1 | 48.3 | 0 | 1 | V2010,A2011,T2012,T2013,T2014 | 0 |
| UW.GSM | 1 | 49.3 | 0 | 1 | V2010,A2011,T2012,T2013,T2014 | 0 |
| UW.SBES | 1 | 51.0 | 0 | 1 | V2010,A2011,T2012,T2013 | 0 |
| UW.RVC | 1 | 51.6 | 0 | 1 | V2010,A2011,T2012,T2013,T2014 | 0 |
| CC.OBSR | 1 | 63.7 | 0 | 3 | V2010,A2011,T2012,T2013,T2014 | 0 |
| UW.STAR | 1 | 68.6 | 0 | 3 | V2010,A2011,T2012,T2013,T2014 | 0 |
| UW.FMW | 1 | 69.5 | 0 | 3 | V2010,A2011,T2012,T2013,T2014 | 0 |
| UW.LMW | 1 | 69.5 | 0 | 1 | V2010,A2011,T2012 | 0 |
| UW.RCS | 1 | 70.5 | 0 | 1 | V2010,A2011,T2012,T2013,T2014 | 1 |
| UW.LCW | 1 | 71.2 | 0 | 1 | V2010,A2011,T2012,T2013,T2014 | 0 |
| UW.LEBA | 1 | 72.2 | 0 | 3 | V2010,A2011,T2012,T2013,T2014 | 0 |
| UW.RCM | 1 | 73.1 | 0 | 3 | V2010,A2011,T2012,T2013,T2014 | 0 |
| XD.MB05 | 1 | 74.8 | 0 | 3 | T2014 | 0 |

## Clear-template sites per family (before = the 42 trial stations; after = + added stations / epochs)

| region | families | ≥3 sites before | ≥3 sites after | projected ≥3 (P ≥ 0.5) | ≥2 before | ≥2 after |
|---|---|---|---|---|---|---|
| Olympic_2010 | 75 | 1 | 1 | 1 | 2 | 2 |
| Olympic_2010_armB | 43 | 0 | 1 | 0 | 2 | 2 |
| Royer_NW_WA | 132 | 24 | 28 | 27 | 55 | 55 |
| ALL | 250 | 25 | 30 | 28 | 59 | 59 |

## Before / after per episode

Rates per family-day on the family's observable tremor days (label only). Separation = (data + 1) / (scrambled control + 1) per family; the same statistic on the control = (second scrambled control + 1) / (scrambled control + 1). `sep ≥ 2` counts families; `ctl` = how many the control produces by the same rule.

| window | pass | fam observable | fam ≥3 sites | sep≥2 data / ctl (core) | sep≥2 data / ctl, ≥3 sites | sep≥3 data / ctl | data / ctl / ctl2 rate (core) | sep≥2 data / ctl (core_ep) | val recall / prec (med, fam ≥20 truth) |
|---|---|---|---|---|---|---|---|---|---|
| V2010 | before | 59 | 25 | 15 / 4 | 9 / 3 | 7 / 0 | 10.92 / 9.08 / 9.19 | 15 / 4 | 0.12 / 0.46 |
| V2010 | after | 59 | 27 | 15 / 4 | 10 / 3 | 7 / 0 | 10.96 / 9.08 / 9.18 | 15 / 4 | 0.12 / 0.46 |
| T2013 | before | 59 | 21 | 9 / 2 | 7 / 1 | 3 / 0 | 10.24 / 9.22 / 8.61 | 1 / 0 | – |
| T2013 | after_p1 | 59 | 25 | 9 / 2 | 7 / 1 | 3 / 0 | 10.24 / 9.22 / 8.61 | 1 / 0 | – |
| T2013 | after | 59 | 25 | 9 / 2 | 7 / 1 | 3 / 0 | 10.24 / 9.22 / 8.61 | 1 / 0 | – |

Separation curve (families with sep ≥ x, core tier; data vs control):

| window | pass | ≥1.5 | ≥2 | ≥3 | ≥5 |
|---|---|---|---|---|---|
| V2010 | before | 18 / 9 | 15 / 4 | 7 / 0 | 0 / 0 |
| V2010 | after | 18 / 9 | 15 / 4 | 7 / 0 | 1 / 0 |
| T2013 | before | 10 / 5 | 9 / 2 | 3 / 0 | 2 / 0 |
| T2013 | after_p1 | 10 / 5 | 9 / 2 | 3 / 0 | 2 / 0 |
| T2013 | after | 10 / 5 | 9 / 2 | 3 / 0 | 2 / 0 |

## Families ranked by separation (after; core tier)

Top 30 per window. Every family is in `results/step4b/step4b_by_family.csv` (all runs, all columns).

### V2010

| family | region | sites before → after (med obs) | tpl sites b/a | data / ctl rate | sep (ctl sep) | sep before | recall / prec | flag |
|---|---|---|---|---|---|---|---|---|
| OLB:514 | Olympic_2010_armB | 2 → 3 | 2 / 3 | 8.00 / 0.50 | 8.50 (0.50) | 4.00 | 0.06 / 0.86 |  |
| RNW:lag15:2338 | Royer_NW_WA | 3 → 3 | 3 / 3 | 10.00 / 2.00 | 4.20 (1.00) | 4.20 | 0.11 / 0.90 |  |
| RNW:run1:1756 | Royer_NW_WA | 4 → 4 | 4 / 4 | 6.50 / 1.50 | 3.50 (1.25) | 3.50 | 0.24 / 1.00 |  |
| RNW:fair:1439 | Royer_NW_WA | 2 → 2 | 2 / 2 | 2.50 / 0.50 | 3.00 (1.00) | 3.00 | 0.05 / 0.33 |  |
| RNW:mute:1174 | Royer_NW_WA | 2 → 2 | 2 / 2 | 1.00 / 0.00 | 3.00 (2.00) | 3.00 | 0.00 / 0.00 |  |
| RNW:lag15:1352 | Royer_NW_WA | 5 → 5 | 5 / 5 | 7.00 / 2.00 | 3.00 (1.00) | 3.00 | 0.33 / 0.75 |  |
| RNW:lag15:1024 | Royer_NW_WA | 3 → 3 | 3 / 3 | 1.00 / 0.00 | 3.00 (2.00) | 3.00 | 0.00 / – |  |
| RNW:fair:507 | Royer_NW_WA | 2 → 2 | 2 / 2 | 7.50 / 2.50 | 2.67 (1.00) | 2.67 | 0.00 / – |  |
| RNW:mute:2435 | Royer_NW_WA | 4 → 4 | 4 / 4 | 31.00 / 11.50 | 2.62 (1.46) | 2.62 | 0.31 / 0.50 |  |
| RNW:run1:1470 | Royer_NW_WA | 3 → 3 | 3 / 3 | 2.00 / 0.50 | 2.50 (2.00) | 2.50 | 0.00 / 0.00 |  |
| RNW:lag15:1236 | Royer_NW_WA | 2 → 2 | 2 / 2 | 6.00 / 2.50 | 2.17 (1.17) | 2.17 | 0.00 / 0.00 |  |
| RNW:lag15:2110 | Royer_NW_WA | 3 → 3 | 4 / 4 | 0.50 / 0.00 | 2.00 (1.00) | 2.00 | 0.00 / – |  |
| RNW:lag15:1356 | Royer_NW_WA | 5 → 5 | 5 / 5 | 6.50 / 3.00 | 2.00 (0.86) | 2.00 | 0.06 / 0.50 |  |
| RNW:lag15:2013 | Royer_NW_WA | 2 → 2 | 2 / 2 | 4.50 / 2.00 | 2.00 (1.00) | 2.00 | 0.00 / – |  |
| RNW:lag15:1087 | Royer_NW_WA | 3 → 3 | 3 / 3 | 0.50 / 0.00 | 2.00 (2.00) | 2.00 | – / – |  |
| RNW:mute:808 | Royer_NW_WA | 3 → 3 | 3 / 3 | 6.00 / 3.00 | 1.86 (1.14) | 1.86 | 0.00 / 0.00 |  |
| RNW:lag15:3076 | Royer_NW_WA | 3 → 3 | 3 / 3 | 18.50 / 10.50 | 1.73 (1.14) | 1.73 | 0.14 / 0.38 |  |
| RNW:lag15:553 | Royer_NW_WA | 5 → 5 | 5 / 5 | 2.00 / 1.00 | 1.67 (1.33) | 1.67 | 0.50 / 0.50 |  |
| RNW:mute:1144 | Royer_NW_WA | 3 → 3 | 3 / 3 | 1.50 / 1.00 | 1.33 (1.67) | 1.33 | 0.00 / – |  |
| RNW:mute:2603 | Royer_NW_WA | 3 → 3 | 3 / 3 | 39.50 / 30.50 | 1.29 (1.00) | 1.29 | 0.30 / 0.41 |  |
| RNW:lag15:2442 | Royer_NW_WA | 2 → 2 | 2 / 3 | 4.00 / 3.00 | 1.29 (1.86) | 1.29 | 0.00 / 0.00 |  |
| RNW:mute:1476 | Royer_NW_WA | 2 → 2 | 2 / 2 | 4.00 / 3.00 | 1.29 (1.29) | 1.29 | 0.07 / 0.67 |  |
| RNW:lag15:2224 | Royer_NW_WA | 2 → 2 | 2 / 2 | 32.50 / 27.00 | 1.20 (0.96) | 1.20 | 0.10 / 0.44 |  |
| RNW:mute:2607 | Royer_NW_WA | 2 → 2 | 2 / 2 | 32.50 / 27.00 | 1.20 (0.98) | 1.20 | 0.00 / 0.00 |  |
| RNW:fair:946 | Royer_NW_WA | 2 → 2 | 2 / 2 | 2.50 / 2.00 | 1.20 (1.00) | 1.20 | 0.00 / – |  |
| RNW:lag15:1370 | Royer_NW_WA | 2 → 2 | 2 / 2 | 6.50 / 5.50 | 1.17 (1.17) | 1.17 | 0.00 / 0.00 |  |
| RNW:mute:1608 | Royer_NW_WA | 2 → 2 | 2 / 3 | 3.50 / 3.00 | 1.14 (1.57) | 1.14 | 0.00 / – |  |
| RNW:run1:2341 | Royer_NW_WA | 2 → 2 | 2 / 3 | 3.50 / 3.00 | 1.14 (1.43) | 1.14 | 0.00 / 0.00 |  |
| RNW:run1:2671 | Royer_NW_WA | 2 → 2 | 2 / 2 | 29.50 / 26.00 | 1.13 (0.98) | 1.13 | 0.12 / 0.75 |  |
| RNW:lag15:2226 | Royer_NW_WA | 2 → 2 | 2 / 2 | 36.50 / 32.50 | 1.12 (0.92) | 1.12 | 0.32 / 0.48 |  |

### T2013

| family | region | sites before → after (med obs) | tpl sites b/a | data / ctl rate | sep (ctl sep) | sep before | recall / prec | flag |
|---|---|---|---|---|---|---|---|---|
| RNW:run1:1756 | Royer_NW_WA | 4 → 4 | 4 / 4 | 12.00 / 1.00 | 6.50 (1.00) | 6.50 | – / – |  |
| RNW:lag15:1352 | Royer_NW_WA | 5 → 5 | 5 / 5 | 9.00 / 1.00 | 5.00 (1.50) | 5.00 | – / – |  |
| OLB:514 | Olympic_2010_armB | 2 → 2 | 2 / 3 | 5.00 / 1.00 | 3.00 (1.00) | 3.00 | – / – |  |
| RNW:lag15:2338 | Royer_NW_WA | 3 → 3 | 3 / 3 | 18.00 / 7.00 | 2.38 (0.62) | 2.38 | – / – |  |
| RNW:lag15:1356 | Royer_NW_WA | 5 → 5 | 5 / 5 | 10.00 / 4.00 | 2.20 (0.60) | 2.20 | – / – |  |
| RNW:lag15:1942 | Royer_NW_WA | 2 → 2 | 2 / 2 | 1.00 / 0.00 | 2.00 (1.00) | 2.00 | – / – |  |
| RNW:mute:1829 | Royer_NW_WA | 3 → 3 | 4 / 4 | 1.00 / 0.00 | 2.00 (2.00) | 2.00 | – / – |  |
| RNW:fair:1818 | Royer_NW_WA | 3 → 3 | 3 / 3 | 1.00 / 0.00 | 2.00 (1.00) | 2.00 | – / – |  |
| RNW:fair:1144 | Royer_NW_WA | 3 → 3 | 3 / 3 | 1.00 / 0.00 | 2.00 (1.00) | 2.00 | – / – |  |
| RNW:mute:1688 | Royer_NW_WA | 2 → 2 | 2 / 2 | 2.00 / 1.00 | 1.50 (1.00) | 1.50 | – / – |  |
| RNW:fair:1439 | Royer_NW_WA | 2 → 2 | 2 / 2 | 6.00 / 4.00 | 1.40 (1.40) | 1.40 | – / – |  |
| RNW:run1:707 | Royer_NW_WA | 2 → 2 | 3 / 3 | 3.00 / 2.00 | 1.33 (1.33) | 1.33 | – / – |  |
| RNW:mute:1174 | Royer_NW_WA | 2 → 2 | 2 / 2 | 3.00 / 2.00 | 1.33 (1.33) | 1.33 | – / – |  |
| RNW:run1:2671 | Royer_NW_WA | 2 → 2 | 2 / 2 | 28.00 / 22.00 | 1.26 (0.91) | 1.26 | – / – |  |
| RNW:lag15:1370 | Royer_NW_WA | 2 → 2 | 2 / 2 | 9.00 / 7.00 | 1.25 (1.12) | 1.25 | – / – |  |
| RNW:qc:2653 | Royer_NW_WA | 2 → 2 | 2 / 2 | 29.00 / 24.00 | 1.20 (1.00) | 1.20 | – / – |  |
| RNW:run1:1722 | Royer_NW_WA | 2 → 2 | 2 / 2 | 5.00 / 4.00 | 1.20 (0.80) | 1.20 | – / – |  |
| RNW:mute:2603 | Royer_NW_WA | 3 → 3 | 3 / 3 | 31.00 / 26.00 | 1.19 (0.96) | 1.19 | – / – |  |
| RNW:mute:2607 | Royer_NW_WA | 2 → 2 | 2 / 2 | 23.00 / 20.00 | 1.14 (0.90) | 1.14 | – / – |  |
| RNW:lag15:2224 | Royer_NW_WA | 2 → 2 | 2 / 2 | 30.00 / 27.00 | 1.11 (0.93) | 1.11 | – / – |  |
| RNW:fair:542 | Royer_NW_WA | 2 → 2 | 2 / 2 | 28.00 / 26.00 | 1.07 (0.85) | 1.07 | – / – |  |
| RNW:lag15:2226 | Royer_NW_WA | 2 → 2 | 2 / 2 | 28.00 / 26.00 | 1.07 (0.89) | 1.07 | – / – |  |
| RNW:mute:3288 | Royer_NW_WA | 2 → 2 | 2 / 2 | 39.00 / 38.00 | 1.03 (1.00) | 1.03 | – / – |  |
| OLY:938 | Olympic_2010 | 2 → 2 | 2 / 2 | 233.00 / 230.00 | 1.01 (1.00) | 1.01 | – / – | runaway ~200 det/day in data AND scrambled control (B013/B943 pair); kept, flagged |
| RNW:lag15:797 | Royer_NW_WA | 2 → 2 | 2 / 2 | 0.00 / 0.00 | 1.00 (1.00) | 1.00 | – / – |  |
| RNW:lag15:553 | Royer_NW_WA | 4 → 4 | 5 / 5 | 11.00 / 11.00 | 1.00 (0.58) | 1.00 | – / – |  |
| RNW:lag15:641 | Royer_NW_WA | 2 → 2 | 2 / 2 | 1.00 / 1.00 | 1.00 (1.50) | 1.00 | – / – |  |
| RNW:lag15:1087 | Royer_NW_WA | 3 → 3 | 3 / 3 | 4.00 / 4.00 | 1.00 (1.00) | 1.00 | – / – |  |
| OLY:127 | Olympic_2010 | 3 → 3 | 3 / 3 | 0.00 / 0.00 | 1.00 (1.00) | 1.00 | – / – |  |
| RNW:mute:1608 | Royer_NW_WA | 2 → 3 | 2 / 3 | 0.00 / 0.00 | 1.00 (1.00) | 1.00 | – / – |  |

## Flag

- OLY:938: runaway (~200 detections/day in data AND in the scrambled control; its sites are the B013/B943 borehole pair + neighbour). Kept in every table, flagged in `flag`.

