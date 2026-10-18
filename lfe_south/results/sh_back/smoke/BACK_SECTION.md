## Backward extension (before 2018)

Question: can the 27 families be followed back before 2018 with the stations that existed then?
A station counts for a period only if it ran the same sensor model then as during the 2018-2024 template build. Gain-only changes do not affect a normalised cc; other stations are dropped.

**Plain summary.** Only 6 of the 17 S1 stations were running before 2018: NC.KMPB (CMG-3ESPC since 2008-04-29), NC.KMR (STS-2 since 2007-01-30) and PB.B049/B932/B933/B934 (borehole HS-1, installed Sep-Nov 2008). We re-ran the detector with only these stations on the same held-out 2018-2024 days used for the v3b score. 17 stations (v3b): any-family recall 0.68, precision 0.78; S6 (pre-2018 stations): any-family recall 0.00, precision 0.88 (9×MAD, Text S1 rules). PASS: BK1 S6 at 9xMAD (S1 rules, source mode) any-family precision 0.877, recall 0.002; rule: BK2 runs only if precision >= 0.5

### Stations: full FDSN epoch history

| netsta   | role           | first_start   | last_end   | build_sensor                                    | same_sensor_since   | operating_before_2018   | usable_before_2018   | band_sensor_changes                                                                                                                                         |
|:---------|:---------------|:--------------|:-----------|:------------------------------------------------|:--------------------|:------------------------|:---------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------|
| BK.BLAS  | S1 detection   | 2022-06-29    | 2100-01-01 | not operating 2021-06                           | nan                 | False                   | False                | 2022-06 off; 2022-07 HH MBB-2 100sps                                                                                                                        |
| BK.BRIC  | S1 detection   | 2019-02-28    | 2100-01-01 | HH CMG-3T-PH                                    | 2019-02-28          | False                   | False                | 2019-02 HH CMG-3T-PH 100sps                                                                                                                                 |
| BK.GHOP  | S1 detection   | 2022-05-05    | 2100-01-01 | not operating 2021-06                           | nan                 | False                   | False                | 2022-05 off; 2022-05 HH MBB-2 100sps                                                                                                                        |
| BK.HALS  | S1 detection   | 2020-06-19    | 2100-01-01 | HH CMG-3T-PH                                    | 2020-06-19          | False                   | False                | 2020-06 HH CMG-3T-PH 100sps; 2023-09 HH MBB-2 100sps                                                                                                        |
| BK.JEWT  | S1 detection   | 2023-08-23    | 2100-01-01 | not operating 2021-06                           | nan                 | False                   | False                | 2023-08 off; 2023-09 HH STS-5A 100sps                                                                                                                       |
| BK.PRDS  | S1 detection   | 2022-09-09    | 2100-01-01 | not operating 2021-06                           | nan                 | False                   | False                | 2022-09 off; 2022-09 HH STS-5A 100sps                                                                                                                       |
| BK.SPRL  | S1 detection   | 2023-09-13    | 2100-01-01 | not operating 2021-06                           | nan                 | False                   | False                | 2023-09 off; 2023-09 HH MBB-2 100sps                                                                                                                        |
| BK.RVIT  | S1 detection   | 2019-10-01    | 2100-01-01 | HH CMG-3T-PH                                    | 2019-10-01          | False                   | False                | 2019-10 HH CMG-3T-PH 100sps                                                                                                                                 |
| BK.HUNT  | S1 detection   | 2022-06-28    | 2100-01-01 | not operating 2021-06                           | nan                 | False                   | False                | 2022-06 off; 2022-07 HH STS-5A 100sps                                                                                                                       |
| NC.KMPB  | S1 detection   | 2000-04-27    | 2100-01-01 | HH CMG-3ESPC                                    | 2008-04-29          | True                    | True                 | 2000-04 off; 2000-05 HH CMG-40T-1 100sps; 2000-06 BH CMG-40T-1 50sps; 2002-09 HH CMG-40T-1 100sps; 2008-05 HH CMG-3ESPC 100sps                              |
| NC.KMR   | S1 detection   | 2003-07-25    | 2100-01-01 | HH STS-2                                        | 2007-01-30          | True                    | True                 | 2003-07 off; 2003-08 HH CMG-40T-1 100sps; 2007-02 HH STS-2 100sps; 2025-09 HH TRILLIUM HORIZON 120 100sps; 2026-03 HH TRILLIUM COMPACT POSTHOLE 100sps      |
| NC.KPR   | S1 detection   | 2019-05-22    | 2100-01-01 | EH L-4C VERTICAL                                | 2019-05-22          | False                   | False                | 2019-05 EH L-4C VERTICAL 100sps                                                                                                                             |
| NC.KBU   | S1 detection   | 2019-05-22    | 2100-01-01 | EH L-4C VERTICAL                                | 2019-05-22          | False                   | False                | 2019-05 EH L-4C VERTICAL 100sps                                                                                                                             |
| PB.B049  | S1 detection   | 2008-11-12    | 2100-01-01 | EH HS-1-LT/Quanterra 330 Linear Phase Composite | 2008-11-12          | True                    | True                 | 2008-11 EH HS-1-LT/Quanterra 330 Linear Phase Composite 100sps                                                                                              |
| PB.B932  | S1 detection   | 2008-10-17    | 2100-01-01 | EH HS-1-LT/Quanterra 330 Linear Phase Composite | 2008-10-17          | True                    | True                 | 2008-10 EH HS-1-LT/Quanterra 330 Linear Phase Composite 100sps                                                                                              |
| PB.B934  | S1 detection   | 2008-10-14    | 2100-01-01 | EH HS-1-LT/Quanterra 330 Linear Phase Composite | 2008-10-14          | True                    | True                 | 2008-10 EH HS-1-LT/Quanterra 330 Linear Phase Composite 100sps                                                                                              |
| PB.B933  | S1 detection   | 2008-09-13    | 2100-01-01 | EH HS-1-LT/Quanterra 330 Linear Phase Composite | 2008-09-13          | True                    | True                 | 2008-09 EH HS-1-LT/Quanterra 330 Linear Phase Composite 100sps                                                                                              |
| NC.KBS   | location stage | 1989-09-14    | 2021-11-05 | EH L-4C VERTICAL                                | 2011-10-27          | True                    | False                | 1989-09 off; 1989-09 EH L-4C VERTICAL 100sps; 2002-10 SH L-4C VERTICAL 50sps; 2011-10 EH L-4C VERTICAL 100sps; 2021-11 off                                  |
| NC.KBSB  | location stage | 2024-05-15    | 2100-01-01 | not operating 2021-06                           | nan                 | False                   | False                | 2024-05 off; 2024-05 EH L-4C VERTICAL 100sps                                                                                                                |
| NC.KKP   | location stage | 1989-08-12    | 2020-11-05 | not operating 2021-06                           | nan                 | True                    | False                | 1989-08 off; 1989-08 EH L-4C VERTICAL 100sps; 2002-09 SH L-4C VERTICAL 50sps; 2011-11 EH L-4C VERTICAL 100sps; 2020-11 off                                  |
| NC.KIP   | location stage | 1989-11-30    | 2100-01-01 | EH L-4C VERTICAL                                | 2011-10-27          | True                    | False                | 1989-11 off; 1989-12 EH L-4C VERTICAL 100sps; 2002-10 SH L-4C VERTICAL 50sps; 2011-11 EH L-4C VERTICAL 100sps                                               |
| NC.KCPB  | location stage | 1999-08-03    | 2022-11-15 | HH CMG-3ESPC                                    | 2008-04-28          | True                    | False                | 1999-08 BH CMG-40T-1 50sps; 2002-10 HH CMG-40T-1 100sps; 2008-05 HH CMG-3ESPC 100sps; 2022-11 off                                                           |
| NC.KCPC  | location stage | 2023-06-12    | 2100-01-01 | not operating 2021-06                           | nan                 | False                   | False                | 2023-06 off; 2023-06 HH MBB-2 100sps                                                                                                                        |
| NC.KBN   | location stage | 1988-10-19    | 2100-01-01 | EH L-4C VERTICAL                                | 2021-02-26          | True                    | False                | 1988-10 off; 1988-10 EH L-4C VERTICAL 100sps; 2002-10 SH L-4C VERTICAL 50sps; 2011-10 EH L-4C VERTICAL 100sps; 2020-11 off; 2021-02 EH L-4C VERTICAL 100sps |
| NC.KFP   | location stage | 1988-10-19    | 2100-01-01 | EH L-4C VERTICAL                                | 2020-11-03          | True                    | False                | 1988-10 off; 1988-10 EH L-4C VERTICAL 100sps; 2002-10 SH L-4C VERTICAL 50sps; 2011-10 EH L-4C VERTICAL 100sps; 2020-10 off; 2020-11 EH L-4C VERTICAL 100sps |
| NC.KCR   | location stage | 1989-08-11    | 2100-01-01 | EH L-4C VERTICAL                                | 2020-11-05          | True                    | False                | 1989-08 off; 1989-08 EH L-4C VERTICAL 100sps; 2002-09 SH L-4C VERTICAL 50sps; 2011-11 EH L-4C VERTICAL 100sps; 2020-09 off; 2020-11 EH L-4C VERTICAL 100sps |
| NC.KSM   | location stage | 1989-08-08    | 2100-01-01 | not operating 2021-06                           | nan                 | True                    | False                | 1989-08 off; 1989-08 EH L-4C VERTICAL 100sps; 2002-08 SH L-4C VERTICAL 50sps; 2011-11 EH L-4C VERTICAL 100sps; 2020-11 off; 2023-01 EH L-4C VERTICAL 100sps |
| PB.B045  | location stage | 2008-09-17    | 2100-01-01 | EH HS-1-LT/Quanterra 330 Linear Phase Composite | 2008-09-17          | True                    | False                | 2008-09 EH HS-1-LT/Quanterra 330 Linear Phase Composite 100sps                                                                                              |
| PB.B046  | location stage | 2008-08-26    | 2100-01-01 | EH HS-1-LT/Quanterra 330 Linear Phase Composite | 2008-08-26          | True                    | False                | 2008-08 EH HS-1-LT/Quanterra 330 Linear Phase Composite 100sps                                                                                              |
| PB.B047  | location stage | 2008-09-04    | 2100-01-01 | EH HS-1-LT/Quanterra 330 Linear Phase Composite | 2008-09-04          | True                    | False                | 2008-09 EH HS-1-LT/Quanterra 330 Linear Phase Composite 100sps                                                                                              |
| PB.B935  | location stage | 2008-10-29    | 2100-01-01 | EH HS-1-LT/Quanterra 330 Linear Phase Composite | 2008-10-29          | True                    | False                | 2008-10 EH HS-1-LT/Quanterra 330 Linear Phase Composite 100sps                                                                                              |

### BK1: known-answer test (385 held-out 2018-2024 test days, scored exactly as v3b)

| stations               |   xMAD | mode   |   ours |   his |   recall_own |   precision_own |   recall_any |   precision_any |   same_fam_share |
|:-----------------------|-------:|:-------|-------:|------:|-------------:|----------------:|-------------:|----------------:|-----------------:|
| 17 stations (v3b)      |      8 | source |  37973 | 26895 |        0.562 |           0.398 |        0.809 |           0.573 |            0.695 |
| 17 stations (v3b)      |      8 | fam    |   5361 |  5021 |        0.520 |           0.487 |        0.558 |           0.523 |            0.932 |
| 17 stations (v3b)      |      9 | source |  23534 | 26895 |        0.485 |           0.555 |        0.685 |           0.782 |            0.709 |
| 17 stations (v3b)      |      9 | fam    |   3651 |  5021 |        0.453 |           0.623 |        0.483 |           0.665 |            0.937 |
| 17 stations (v3b)      |     10 | source |  16307 | 26895 |        0.392 |           0.647 |        0.541 |           0.892 |            0.725 |
| 17 stations (v3b)      |     10 | fam    |   2613 |  5021 |        0.369 |           0.709 |        0.392 |           0.753 |            0.942 |
| S6 (pre-2018 stations) |      8 | source |     98 | 26895 |        0.001 |           0.408 |        0.002 |           0.673 |            0.606 |
| S6 (pre-2018 stations) |      8 | fam    |     22 |  5021 |        0.002 |           0.545 |        0.002 |           0.545 |            1.000 |
| S6 (pre-2018 stations) |      9 | source |     57 | 26895 |        0.001 |           0.544 |        0.002 |           0.877 |            0.620 |
| S6 (pre-2018 stations) |      9 | fam    |     13 |  5021 |        0.002 |           0.769 |        0.002 |           0.769 |            1.000 |
| S6 (pre-2018 stations) |     10 | source |     39 | 26895 |        0.001 |           0.564 |        0.001 |           0.872 |            0.647 |
| S6 (pre-2018 stations) |     10 | fam    |     10 |  5021 |        0.002 |           0.800 |        0.002 |           0.800 |            1.000 |

PASS: BK1 S6 at 9xMAD (S1 rules, source mode) any-family precision 0.877, recall 0.002; rule: BK2 runs only if precision >= 0.5

### BK2: 2008-2017 scan

```
BK2 plan: 20650 station-days (3058 already cached, 17592 new) -> est. 457 GB miniSEED + ~366 GB 20-sps cache; df -h /home/jovyan: prod-jupyter-homes.intra.oceanobservatories.org:/ifs/juphome/pvc-631c5c51-19c4-4b49-867d-384217e76f89  100T   35T   66T  35% /home/jovyan
      size         min
sta                   
B049  3337  2008-11-12
B932  3363  2008-10-17
B933  3397  2008-09-13
B934  3366  2008-10-14
KMPB  3534  2008-04-29
KMR   3653  2008-01-01
```

57 S1 9× detections 2008-2017 on 2 processed days (live stations per day in `days_processed_bk2.csv`).
Per year: 2022: 57

![2008-2026](sh_back/monthly_rate_2008_2026.png)

|   fam |   n |   days |   ratio_bk1 |   rate |   rate_his_equiv |
|------:|----:|-------:|------------:|-------:|-----------------:|
|     1 |   0 |      0 |       0.002 |    nan |              nan |
|     2 |   0 |      0 |       0.000 |    nan |              nan |
|     3 |   0 |      0 |       0.000 |    nan |              nan |
|     4 |   0 |      0 |       0.000 |    nan |              nan |
|     5 |   0 |      0 |       0.000 |    nan |              nan |
|     6 |   0 |      0 |       0.001 |    nan |              nan |
|     7 |   0 |      0 |       0.000 |    nan |              nan |
|     8 |   0 |      0 |       0.004 |    nan |              nan |
|     9 |   0 |      0 |       0.000 |    nan |              nan |
|    10 |   0 |      0 |       0.003 |    nan |              nan |
|    11 |   0 |      0 |       0.000 |    nan |              nan |
|    12 |   0 |      0 |       0.002 |    nan |              nan |
|    13 |   0 |      0 |       0.001 |    nan |              nan |
|    14 |   0 |      0 |       0.000 |    nan |              nan |
|    15 |   0 |      0 |       0.000 |    nan |              nan |
|    16 |   0 |      0 |       0.005 |    nan |              nan |
|    17 |   0 |      0 |       0.002 |    nan |              nan |
|    18 |   0 |      0 |       0.005 |    nan |              nan |
|    19 |   0 |      0 |       0.000 |    nan |              nan |
|    20 |   0 |      0 |       0.009 |    nan |              nan |
|    21 |   0 |      0 |       0.000 |    nan |              nan |
|    22 |   0 |      0 |       0.002 |    nan |              nan |
|    23 |   0 |      0 |       0.004 |    nan |              nan |
|    24 |   0 |      0 |       0.005 |    nan |              nan |
|    25 |   0 |      0 |       0.000 |    nan |              nan |
|    26 |   0 |      0 |       0.000 |    nan |              nan |
|    27 |   0 |      0 |       0.005 |    nan |              nan |

Files (results/sh_back/): `catalogue_bk2_S1_9x.{parquet,csv}`, `candidates_bk2_raw8x.parquet`, `daily_counts_bk2.csv`, `monthly_bk2.csv`.
