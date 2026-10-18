# DUC_FIN: Ducellier & Creager (2022) families finished on the data on disk (2008-01 .. 2012-07)

_Auto-written by `scripts/duc_fin_report.py`, 2026-10-01 21:51 UTC; queue `scripts/duc_fin_queue.sh`; outputs `duc_fin/` and `cut/duc/`._

## Plain-language summary

**1. Why we find more events than Ducellier in May-Dec 2011 (decision point).** His catalogue lost stations; our scan did not (full diagnosis: [`D2011.md`](duc_fin/D2011.md)).

* **His catalogue lost stations from May 2011; ours did not.** His median channel count for his rule events fell from 18 (2008-01..2011-04, monthly median) to 15 (May-Dec 2011), and for all his detections from 14 to 8; ours stayed at 19 -> 18. With his full station list, his channel count is explained exactly for 38 of 60 families per month before May 2011 and for 0 of 60 after; removing a few stations explains 53 of 60. Stations inferred missing in >= 3 of the May-Dec 2011 months but in none of 2010-09..2011-04: GASB, KHBB, WDC, YBH (BK GASB / WDC / YBH changed location code to 00 on 2011-06-16 / 05-06 / 06-03); in Nov-Dec 2011 many NC analogue stations follow (their SHZ -> EHZ switch).
* **Fewer channels make his catalogue rule harder to pass** (the rule sums the correlation over channels). Of our 3260 extra events, 1035 are also in his RAW detections of the same family within 2 s with the same correlation (his 0.117 vs ours 0.115) but 12 instead of 16 channels: 92% of them would pass his rule with his correlation and our channel count.
* **The extras look like ordinary repeats, not an artefact.** They spread over 60 families (largest share 5%), 81% fall within 1 h of one of his catalogued events (matched events: 100%; extras before May 2011: 72%), median 6 per hour (matched: 9); only 9% sit in hours with >= 20 extras (matched 16%). Their rule margin (summed cc / his threshold) is 1.14 (matched 1.48): weaker events near his threshold.
* **Before May 2011 the match is one-to-one:** same-family recall 0.931, precision 0.932 over 2008-01..2011-04 (his rule both sides, +-2 s).
* **Decision (part A): his data gap, not our artefact.** Nothing is re-scanned. May-Dec 2011 is kept, evaluated separately in the confidence step (his catalogue there is incomplete) and never used to fit. Part B below re-sums our events on his remaining stations and adds the scrambled control.

**2. Re-templating of the new 2011-10..12 station epochs** ([`retpl/epoch_summary.csv`](duc_fin/retpl/epoch_summary.csv)): 26 new epochs start inside the data on disk; 23 re-templated, 3 dropped (no family-channel beat its scrambled control); 174 of 277 family-channels kept. 17 epochs start after 2012-07 and stay pending (no data on disk; no downloads in this run). Kept templates vs his template of the old channel: median cc 0.41, median |cc| 0.41; 17 of 174 best matched with reversed polarity ([`retpl/template_vs_his.csv`](duc_fin/retpl/template_vs_his.csv)). 
Templates per family-day in 2012: 11.8 before -> 14.5 after; events passing his rule 2012-01..07: 2443 before -> 4074 after ([`joins_by_month.csv`](duc_fin/joins_by_month.csv)).

Ours / his under his rule (final scan): 2008 0.999, 2009 1.014, 2010 0.992, 2011 1.498 ([`joins_by_year.csv`](duc_fin/joins_by_year.csv)).

**3. Windows for dv/v** ([`cut/duc/`](cut/duc)): 25,527,230 windows (-15..+30 s, raw counts) around 321,611 candidates, 50.1 GB cut from 708 GB raw. Arrival offsets per window: {'moveout_cvm17': 22326686, 'template': 3200544}. Check: PASS (PASS: 20 random windows re-read from HDF5 vs the raw miniSEED (obspy slice): 20 sample-identical, 58500 samples compared, gap counts consistent 20/20). Stack check: 080408.16.026 at BK.GASB 2011-05, 212 windows, stack vs his template cc 0.65 (single windows 0.085, shifted control 0.234) -> PASS.

**4. Per-repeat confidence P** (probability that a candidate is one of his catalogue events; fitted on 2008-2010 days, out-of-fold by day). Independent checks, AUC same family / any family: held-out months 0.997 / 0.993, Jan-Apr 2011 0.999 / 0.997; his rule flag alone 0.979 / 0.921, peak height (xmad) alone 0.689. Without PNET: 0.997 / 0.993. Calibration slope (target 0.8-1.2) held-out months 1.00, Jan-Apr 2011 0.99; matched / sum P 1.00 and 0.98. Scrambled control: P > 0.5 for 1.7% of control candidates vs 23.2% of data candidates (held-out months). May-Dec 2011 (his catalogue incomplete, not fitted): AUC 0.990, matched / sum P 0.44 (below 1 = we find more than his catalogue, as expected from step 1). Note: our detector reproduces his and his catalogue is his rule on the same statistic, so P is close to the his-rule flag; the gain is the ranking inside / near his threshold and a P that stays defined when channel counts change (2011-2012). Against his RAW detections (all 66 families) the FULL model reaches AUC 0.898 / 0.909 (same / any family). P is written for every candidate, 2012-01..07 included ([`P_by_candidate.parquet`](duc_fin/P_by_candidate.parquet), joinable [`conf_P_by_det_id.parquet`](cut/duc/conf_P_by_det_id.parquet)).

PNET fine-tune on his events: internal validation AUC 0.547 (off-the-shelf) -> 0.742 (kept step 500); learning curve [`figures/pnet_learning_curve.png`](duc_fin/figures/pnet_learning_curve.png).

**5. Delete candidates (nothing deleted)**: candidates: 138942 files, 1.315 TB -> /home/jovyan/tremorferometry/lfe_south/results/cut/DELETE_CANDIDATES_ducellier.csv; kept: 31095 files, 0.344 TB -> KEEP_ducellier.csv (summary KEEP_ducellier_summary.csv) ([`DELETE_CANDIDATES_ducellier.txt`](cut/DELETE_CANDIDATES_ducellier.txt)).

![monthly](duc_fin/figures/monthly_counts.png)

## 1. 2011 diagnosis

See [`D2011.md`](duc_fin/D2011.md) (one-to-one matching by month, his channel counts, stations his catalogue lost, the extras' margin / clustering / co-occurrence, part B: extras re-summed on his station set and the scrambled control).

## 2. Re-templating and the re-scan

| epoch   | sta   | start      | chans       |      sr | change                                                            | on_disk   |   family_channels |   kept |   dropped |   families_kept |   families |   median_N |   median_cc_his |   median_abscc_his | decision                                                                  | reasons                                                                         |
|:--------|:------|:-----------|:------------|--------:|:------------------------------------------------------------------|:----------|------------------:|-------:|----------:|----------------:|-----------:|-----------:|----------------:|-------------------:|:--------------------------------------------------------------------------|:--------------------------------------------------------------------------------|
| GASB#3  | GASB  | 2024-04-18 | BHE,BHN,BHZ |  40.000 | sensor:STS-2->STS-7                                               | False     |                 0 |      0 |         0 |               0 |          0 |    nan     |         nan     |            nan     | PENDING: starts after 2012-07, no data on disk (no downloads in this run) | nan                                                                             |
| GBB#3   | GBB   | 2024-04-17 | HHZ         | 100.000 | loc:->01;channels:EHZ->HHZ;sensor:L-4C VERTICAL->MBB-2            | False     |                 0 |      0 |         0 |               0 |          0 |    nan     |         nan     |            nan     | PENDING: starts after 2012-07, no data on disk (no downloads in this run) | nan                                                                             |
| GSN#2   | GSN   | 2019-03-28 | EHZ         | 100.000 | loc:->01                                                          | False     |                 0 |      0 |         0 |               0 |          0 |    nan     |         nan     |            nan     | PENDING: starts after 2012-07, no data on disk (no downloads in this run) | nan                                                                             |
| GTC#2   | GTC   | 2011-10-27 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                 3 |      2 |         1 |               2 |          3 |    150.000 |           0.121 |              0.121 | RE-TEMPLATED                                                              | {"kept": 2, "stack not above every scrambled control": 1}                       |
| HOPS#3  | HOPS  | 2024-10-17 | BHE,BHN,BHZ |  40.000 | sensor:STS-1 VERTICAL->STS-7;azimuth:(90.0, 0.0)->(0.0, 270.0)    | False     |                 0 |      0 |         0 |               0 |          0 |    nan     |         nan     |            nan     | PENDING: starts after 2012-07, no data on disk (no downloads in this run) | nan                                                                             |
| HOPS#4  | HOPS  | 2024-11-12 | BHE,BHN,BHZ |  40.000 | azimuth:(0.0, 270.0)->(90.0, 0.0)                                 | False     |                 0 |      0 |         0 |               0 |          0 |    nan     |         nan     |            nan     | PENDING: starts after 2012-07, no data on disk (no downloads in this run) | nan                                                                             |
| KBN#2   | KBN   | 2011-10-27 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                11 |      8 |         3 |               8 |         11 |     65.000 |           0.307 |              0.307 | RE-TEMPLATED                                                              | {"kept": 8, "N<20 windows": 3}                                                  |
| KBN#3   | KBN   | 2021-02-26 | EHZ         | 100.000 | loc:->01                                                          | False     |                 0 |      0 |         0 |               0 |          0 |    nan     |         nan     |            nan     | PENDING: starts after 2012-07, no data on disk (no downloads in this run) | nan                                                                             |
| KBS#2   | KBS   | 2011-10-27 | EHZ         | 100.000 | loc:->02;channels:SHZ->EHZ;sample_rate:50.0->100.0                | True      |                 2 |      0 |         2 |               0 |          2 |    150.000 |         nan     |            nan     | DROPPED (no family-channel passes the scrambled control)                  | {"stack not above every scrambled control": 2}                                  |
| KCS#2   | KCS   | 2011-11-01 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                26 |     24 |         2 |              24 |         26 |     46.500 |           0.448 |              0.448 | RE-TEMPLATED                                                              | {"kept": 24, "N<20 windows": 1, "stack not above every scrambled control": 1}   |
| KFP#2   | KFP   | 2011-10-27 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                 4 |      2 |         2 |               2 |          4 |    150.000 |           0.178 |              0.190 | RE-TEMPLATED                                                              | {"kept": 2, "N<20 windows": 1, "stack not above every scrambled control": 1}    |
| KFP#3   | KFP   | 2020-11-03 | EHZ         | 100.000 | loc:->01                                                          | False     |                 0 |      0 |         0 |               0 |          0 |    nan     |         nan     |            nan     | PENDING: starts after 2012-07, no data on disk (no downloads in this run) | nan                                                                             |
| KHBB#3  | KHBB  | 2022-08-09 | HHE,HHN,HHZ | 100.000 | sensor:CMG-3ESPC->TRILLIUM HORIZON 120                            | False     |                 0 |      0 |         0 |               0 |          0 |    nan     |         nan     |            nan     | PENDING: starts after 2012-07, no data on disk (no downloads in this run) | nan                                                                             |
| KIP#2   | KIP   | 2011-10-27 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                 6 |      3 |         3 |               3 |          6 |     89.500 |           0.234 |              0.234 | RE-TEMPLATED                                                              | {"kept": 3, "N<20 windows": 2, "stack not above every scrambled control": 1}    |
| KIP#3   | KIP   | 2019-08-21 | EHZ         | 100.000 | loc:->01                                                          | False     |                 0 |      0 |         0 |               0 |          0 |    nan     |         nan     |            nan     | PENDING: starts after 2012-07, no data on disk (no downloads in this run) | nan                                                                             |
| KKP#2   | KKP   | 2011-11-01 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                15 |     11 |         4 |              11 |         15 |     56.000 |           0.458 |              0.458 | RE-TEMPLATED                                                              | {"kept": 11, "N<20 windows": 2, "stack not above every scrambled control": 2}   |
| KOM#2   | KOM   | 2011-11-02 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                38 |     16 |        22 |              16 |         38 |     36.000 |           0.409 |              0.411 | RE-TEMPLATED                                                              | {"kept": 16, "N<20 windows": 12, "stack not above every scrambled control": 10} |
| KOM#3   | KOM   | 2013-10-30 | EHZ         | 100.000 | loc:->02                                                          | False     |                 0 |      0 |         0 |               0 |          0 |    nan     |         nan     |            nan     | PENDING: starts after 2012-07, no data on disk (no downloads in this run) | nan                                                                             |
| KOM#4   | KOM   | 2016-09-28 | EHZ         | 100.000 | nan                                                               | False     |                 0 |      0 |         0 |               0 |          0 |    nan     |         nan     |            nan     | PENDING: starts after 2012-07, no data on disk (no downloads in this run) | nan                                                                             |
| KPP#2   | KPP   | 2011-11-01 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                 5 |      0 |         5 |               0 |          5 |     49.000 |         nan     |            nan     | DROPPED (no family-channel passes the scrambled control)                  | {"stack not above every scrambled control": 4, "N<20 windows": 1}               |
| KRK#2   | KRK   | 2011-10-27 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                 5 |      2 |         3 |               2 |          5 |    150.000 |           0.183 |              0.183 | RE-TEMPLATED                                                              | {"N<20 windows": 2, "kept": 2, "stack not above every scrambled control": 1}    |
| KRK#3   | KRK   | 2019-06-04 | EHZ         | 100.000 | loc:->01                                                          | False     |                 0 |      0 |         0 |               0 |          0 |    nan     |         nan     |            nan     | PENDING: starts after 2012-07, no data on disk (no downloads in this run) | nan                                                                             |
| KRP#5   | KRP   | 2021-08-03 | HHE,HHN,HHZ | 100.000 | sensor:CMG-3ESPC->TRILLIUM HORIZON 120                            | False     |                 0 |      0 |         0 |               0 |          0 |    nan     |         nan     |            nan     | PENDING: starts after 2012-07, no data on disk (no downloads in this run) | nan                                                                             |
| KTR#2   | KTR   | 2011-11-02 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                 5 |      1 |         4 |               1 |          5 |      6.000 |           0.692 |              0.692 | RE-TEMPLATED                                                              | {"N<20 windows": 4, "kept": 1}                                                  |
| LBF#2   | LBF   | 2011-11-03 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                 3 |      3 |         0 |               3 |          3 |     83.000 |           0.323 |              0.323 | RE-TEMPLATED                                                              | {"kept": 3}                                                                     |
| LBK#2   | LBK   | 2011-12-08 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                24 |     13 |        11 |              13 |         24 |     43.500 |           0.310 |              0.331 | RE-TEMPLATED                                                              | {"kept": 13, "stack not above every scrambled control": 9, "N<20 windows": 2}   |
| LBP#2   | LBP   | 2011-12-08 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                28 |     21 |         7 |              21 |         28 |     40.000 |           0.475 |              0.475 | RE-TEMPLATED                                                              | {"kept": 21, "N<20 windows": 4, "stack not above every scrambled control": 3}   |
| LDB#2   | LDB   | 2011-12-08 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                 1 |      1 |         0 |               1 |          1 |     26.000 |           0.234 |              0.234 | RE-TEMPLATED                                                              | {"kept": 1}                                                                     |
| LGB#2   | LGB   | 2011-11-03 | EHZ         | 100.000 | loc:->02;channels:SHZ->EHZ;sample_rate:50.0->100.0                | True      |                 2 |      1 |         1 |               1 |          2 |     20.500 |           0.262 |              0.262 | RE-TEMPLATED                                                              | {"kept": 1, "N<20 windows": 1}                                                  |
| LGP#2   | LGP   | 2011-12-08 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                16 |     12 |         4 |              12 |         16 |     37.500 |           0.392 |              0.392 | RE-TEMPLATED                                                              | {"kept": 12, "N<20 windows": 2, "stack not above every scrambled control": 2}   |
| LPG#2   | LPG   | 2011-12-08 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                12 |      8 |         4 |               8 |         12 |     40.000 |           0.619 |              0.619 | RE-TEMPLATED                                                              | {"kept": 8, "N<20 windows": 3, "stack not above every scrambled control": 1}    |
| LRB#2   | LRB   | 2011-12-08 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                13 |      8 |         5 |               8 |         13 |     41.000 |           0.388 |              0.388 | RE-TEMPLATED                                                              | {"kept": 8, "stack not above every scrambled control": 3, "N<20 windows": 2}    |
| LRB#3   | LRB   | 2024-10-28 | EHZ         | 100.000 | loc:->01                                                          | False     |                 0 |      0 |         0 |               0 |          0 |    nan     |         nan     |            nan     | PENDING: starts after 2012-07, no data on disk (no downloads in this run) | nan                                                                             |
| LRR#2   | LRR   | 2011-12-08 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                 1 |      0 |         1 |               0 |          1 |     41.000 |         nan     |            nan     | DROPPED (no family-channel passes the scrambled control)                  | {"stack not above every scrambled control": 1}                                  |
| LSF#2   | LSF   | 2011-12-08 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                25 |     16 |         9 |              16 |         25 |     39.000 |           0.422 |              0.422 | RE-TEMPLATED                                                              | {"kept": 16, "stack not above every scrambled control": 6, "N<20 windows": 3}   |
| LSH#2   | LSH   | 2011-12-08 | EHZ         | 100.000 | loc:->03;channels:SHZ->EHZ;sample_rate:50.0->100.0                | True      |                 7 |      5 |         2 |               5 |          7 |     62.000 |           0.361 |              0.361 | RE-TEMPLATED                                                              | {"kept": 5, "stack not above every scrambled control": 1, "N<20 windows": 1}    |
| LSR#2   | LSR   | 2011-12-08 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                 7 |      5 |         2 |               5 |          7 |     57.000 |           0.325 |              0.325 | RE-TEMPLATED                                                              | {"kept": 5, "stack not above every scrambled control": 1, "N<20 windows": 1}    |
| LTC#2   | LTC   | 2011-12-08 | EHZ         | 100.000 | loc:->02;channels:SHZ->EHZ;sample_rate:50.0->100.0                | True      |                 3 |      3 |         0 |               3 |          3 |     41.000 |           0.283 |              0.316 | RE-TEMPLATED                                                              | {"kept": 3}                                                                     |
| LTC#3   | LTC   | 2019-11-25 | EHZ         | 100.000 | loc:02->03                                                        | False     |                 0 |      0 |         0 |               0 |          0 |    nan     |         nan     |            nan     | PENDING: starts after 2012-07, no data on disk (no downloads in this run) | nan                                                                             |
| LTC#4   | LTC   | 2023-01-25 | HHZ         | 100.000 | channels:EHZ->HHZ;sensor:L-4C VERTICAL->TRILLIUM COMPACT POSTHOLE | False     |                 0 |      0 |         0 |               0 |          0 |    nan     |         nan     |            nan     | PENDING: starts after 2012-07, no data on disk (no downloads in this run) | nan                                                                             |
| LVR#2   | LVR   | 2011-12-08 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                 8 |      5 |         3 |               5 |          8 |     62.500 |           0.614 |              0.614 | RE-TEMPLATED                                                              | {"kept": 5, "stack not above every scrambled control": 2, "N<20 windows": 1}    |
| LWH#2   | LWH   | 2011-12-08 | EHZ         | 100.000 | channels:SHZ->EHZ;sample_rate:50.0->100.0                         | True      |                 7 |      4 |         3 |               4 |          7 |     57.000 |           0.469 |              0.477 | RE-TEMPLATED                                                              | {"kept": 4, "stack not above every scrambled control": 3}                       |
| YBH#3   | YBH   | 2022-03-15 | BHE,BHN,BHZ |  40.000 | sensor:STS-1 VERTICAL->STS-7                                      | False     |                 0 |      0 |         0 |               0 |          0 |    nan     |         nan     |            nan     | PENDING: starts after 2012-07, no data on disk (no downloads in this run) | nan                                                                             |

Joins by year:

|     year |   ours_rule |   his_rule |   ours_raw |    his_raw |   ours_rule_mf1 |   ratio_rule |   precision_same |   recall_same |
|---------:|------------:|-----------:|-----------:|-----------:|----------------:|-------------:|-----------------:|--------------:|
| 2008.000 |   10729.000 |  10740.000 |  44027.000 |  39961.000 |       10682.000 |        0.999 |            0.944 |         0.943 |
| 2009.000 |    8649.000 |   8528.000 |  40863.000 |  37531.000 |        8649.000 |        1.014 |            0.928 |         0.942 |
| 2010.000 |    6560.000 |   6610.000 |  49372.000 |  52622.000 |        6559.000 |        0.992 |            0.927 |         0.920 |
| 2011.000 |    8565.000 |   5717.000 | 156982.000 | 253093.000 |        7941.000 |        1.498 |            0.523 |         0.783 |
| 2012.000 |    4074.000 |      0.000 |  30367.000 |      0.000 |        2443.000 |      nan     |          nan     |       nan     |

Joins by month (2011-2012):

| month   |   ours_rule |   his_rule |   ratio_rule |   recall_same |   precision_same |   ours_raw |   his_raw |   ours_rule_mf1 |   ours_raw_mf1 |   ntpl_famday |   ntpl_famday_mf1 |   new_tpl_famday |
|:--------|------------:|-----------:|-------------:|--------------:|-----------------:|-----------:|----------:|----------------:|---------------:|--------------:|------------------:|-----------------:|
| 2011-01 |         531 |    559.000 |        0.950 |         0.826 |            0.870 |       9108 | 14870.000 |             531 |           9108 |        15.439 |            15.439 |            0.000 |
| 2011-02 |         318 |    319.000 |        0.997 |         0.909 |            0.912 |       6516 | 10027.000 |             318 |           6516 |        15.439 |            15.439 |            0.000 |
| 2011-03 |         186 |    192.000 |        0.969 |         0.755 |            0.780 |      22082 | 37966.000 |             186 |          22077 |        15.439 |            15.454 |            0.000 |
| 2011-04 |        1036 |   1055.000 |        0.982 |         0.923 |            0.940 |       4072 |  4686.000 |            1036 |           4059 |        15.286 |            15.286 |            0.000 |
| 2011-05 |        1116 |    752.000 |        1.484 |         0.839 |            0.565 |      48624 | 78044.000 |            1112 |          48620 |        15.439 |            15.465 |            0.000 |
| 2011-06 |         977 |    656.000 |        1.489 |         0.780 |            0.524 |      31070 | 56011.000 |             976 |          31064 |        15.397 |            15.397 |            0.000 |
| 2011-07 |         226 |     95.000 |        2.379 |         0.474 |            0.199 |      16690 | 24265.000 |             226 |          16689 |        15.633 |            15.633 |            0.000 |
| 2011-08 |         726 |    345.000 |        2.104 |         0.777 |            0.369 |       5423 |  7795.000 |             726 |           5407 |        15.876 |            15.876 |            0.000 |
| 2011-09 |         602 |    406.000 |        1.483 |         0.635 |            0.429 |       4357 |  7029.000 |             602 |           4344 |        16.253 |            16.253 |            0.000 |
| 2011-10 |        1396 |    974.000 |        1.433 |         0.695 |            0.485 |       3400 |  4984.000 |            1396 |           3403 |        16.142 |            16.101 |            0.042 |
| 2011-11 |         554 |    175.000 |        3.166 |         0.651 |            0.206 |       2387 |  3228.000 |             492 |           2397 |        15.406 |            14.296 |            1.109 |
| 2011-12 |         897 |    189.000 |        4.746 |         0.545 |            0.115 |       3253 |  4188.000 |             340 |           3059 |        15.090 |            12.730 |            2.360 |
| 2012-01 |         336 |    nan     |      nan     |       nan     |          nan     |       2788 |   nan     |             185 |           2910 |        14.599 |            11.895 |            2.705 |
| 2012-02 |         864 |    nan     |      nan     |       nan     |          nan     |       2847 |   nan     |             428 |           2772 |        14.734 |            12.016 |            2.719 |
| 2012-03 |         736 |    nan     |      nan     |       nan     |          nan     |       3614 |   nan     |             499 |           3736 |        14.988 |            12.270 |            2.719 |
| 2012-04 |         300 |    nan     |      nan     |       nan     |          nan     |       4407 |   nan     |             225 |           4640 |        14.659 |            11.941 |            2.719 |
| 2012-05 |         689 |    nan     |      nan     |       nan     |          nan     |       5152 |   nan     |             359 |           5246 |        14.344 |            11.626 |            2.719 |
| 2012-06 |         846 |    nan     |      nan     |       nan     |          nan     |       5771 |   nan     |             527 |           5857 |        14.172 |            11.453 |            2.719 |
| 2012-07 |         303 |    nan     |      nan     |       nan     |          nan     |       5788 |   nan     |             220 |           5978 |        14.281 |            11.562 |            2.719 |

## 3. Windows

| month   |   n_det |   n_windows |   n_full |   n_station_days |   n_read_ok |   raw_GB |   cut_GB |   n_segments |   seg_hours |   minutes |
|:--------|--------:|------------:|---------:|-----------------:|------------:|---------:|---------:|-------------:|------------:|----------:|
| 2008-01 |    2965 |      197683 |   197226 |             1386 |        1386 |   10.903 |    0.379 |        92271 |    1408.790 |     0.800 |
| 2008-02 |    2042 |      141236 |   140805 |             1353 |        1353 |   10.502 |    0.264 |        70093 |     982.050 |     0.700 |
| 2008-03 |    2243 |      164750 |   164493 |             1457 |        1457 |   11.033 |    0.373 |       101091 |    1415.050 |     0.800 |
| 2008-04 |    7441 |      746417 |   745719 |             2070 |        2070 |   18.981 |    1.532 |       315340 |    5100.720 |     2.800 |
| 2008-05 |    3541 |      285700 |   285425 |             1677 |        1677 |   13.631 |    0.772 |       208771 |    2848.350 |     1.600 |
| 2008-06 |    4176 |      292475 |   292283 |             1371 |        1371 |   10.514 |    0.774 |       222292 |    2997.700 |     1.500 |
| 2008-07 |    5557 |      395955 |   395651 |             1436 |        1436 |   11.189 |    1.017 |       294427 |    3989.540 |     2.000 |
| 2008-08 |    5175 |      376291 |   375831 |             1457 |        1457 |   11.546 |    0.959 |       275021 |    3761.910 |     1.800 |
| 2008-09 |    3679 |      271135 |   270894 |             1416 |        1416 |   10.994 |    0.660 |       189906 |    2592.800 |     1.300 |
| 2008-10 |    1884 |      142601 |   142360 |             1475 |        1475 |   12.123 |    0.315 |        81940 |    1194.150 |     0.600 |
| 2008-11 |    2494 |      196216 |   195980 |             1466 |        1466 |   12.594 |    0.347 |        84992 |    1312.800 |     0.700 |
| 2008-12 |    2830 |      220819 |   220566 |             1511 |        1511 |   13.174 |    0.464 |       119242 |    1751.240 |     1.000 |
| 2009-01 |    1660 |      129059 |   128814 |             1504 |        1504 |   12.522 |    0.269 |        71933 |    1035.570 |     0.600 |
| 2009-02 |    1273 |       96415 |    96220 |             1359 |        1359 |   11.627 |    0.213 |        57742 |     791.880 |     0.500 |
| 2009-03 |    3152 |      246538 |   246362 |             1519 |        1519 |   12.961 |    0.443 |       110257 |    1662.690 |     0.900 |
| 2009-04 |    2358 |      186093 |   185956 |             1470 |        1470 |   12.143 |    0.489 |       140076 |    1890.280 |     1.100 |
| 2009-05 |    4444 |      343621 |   343407 |             1504 |        1504 |   12.571 |    0.892 |       253930 |    3454.810 |     1.700 |
| 2009-06 |    6349 |      490644 |   490336 |             1456 |        1456 |   11.694 |    1.053 |       284317 |    4149.230 |     2.100 |
| 2009-07 |    5693 |      458717 |   458156 |             1537 |        1537 |   12.462 |    1.075 |       304498 |    4263.330 |     2.100 |
| 2009-08 |    4555 |      374886 |   374547 |             1549 |        1549 |   13.038 |    1.004 |       293247 |    3957.850 |     2.100 |
| 2009-09 |    4450 |      370280 |   369835 |             1501 |        1501 |   12.989 |    0.909 |       255478 |    3542.990 |     1.900 |
| 2009-10 |    3129 |      263768 |   263320 |             1552 |        1552 |   14.031 |    0.538 |       135675 |    2016.960 |     1.100 |
| 2009-11 |    1254 |      103760 |   103534 |             1500 |        1500 |   13.459 |    0.228 |        60548 |     873.400 |     0.500 |
| 2009-12 |    2546 |      212437 |   212167 |             1550 |        1550 |   13.851 |    0.501 |       141534 |    1957.640 |     1.000 |
| 2010-01 |    1956 |      151648 |   151097 |             1521 |        1521 |   12.548 |    0.316 |        84805 |    1228.560 |     0.700 |
| 2010-02 |    1762 |      134748 |   134430 |             1355 |        1355 |   11.164 |    0.269 |        74340 |    1059.550 |     0.600 |
| 2010-03 |    3270 |      256880 |   256530 |             1516 |        1516 |   12.565 |    0.501 |       137937 |    1989.360 |     1.100 |
| 2010-04 |    2520 |      197975 |   197746 |             1470 |        1470 |   11.980 |    0.450 |       130933 |    1823.740 |     1.000 |
| 2010-05 |    3084 |      246161 |   245991 |             1539 |        1539 |   12.298 |    0.597 |       177986 |    2432.720 |     1.500 |
| 2010-06 |    4249 |      339569 |   339292 |             1500 |        1500 |   12.107 |    0.853 |       251497 |    3431.680 |     2.000 |
| 2010-07 |    7101 |      547621 |   547122 |             1532 |        1532 |   12.061 |    1.259 |       375111 |    5242.890 |     3.500 |
| 2010-08 |    5191 |      412547 |   411754 |             1567 |        1567 |   12.384 |    1.002 |       303530 |    4125.200 |     2.300 |
| 2010-09 |    4914 |      399464 |   398804 |             1528 |        1528 |   12.265 |    0.822 |       233277 |    3345.540 |     1.700 |
| 2010-10 |    2773 |      224395 |   224087 |             1578 |        1578 |   13.117 |    0.471 |       132208 |    1903.500 |     1.000 |
| 2010-11 |    1873 |      155852 |   154673 |             1530 |        1530 |   13.320 |    0.359 |        97721 |    1381.210 |     0.700 |
| 2010-12 |   10679 |      798091 |   796264 |             1473 |        1473 |   12.069 |    1.160 |        88338 |    4658.500 |     1.500 |
| 2011-01 |    9108 |      741072 |   740136 |             1550 |        1550 |   12.360 |    1.161 |        93848 |    4849.100 |     1.600 |
| 2011-02 |    6516 |      521642 |   520621 |             1401 |        1401 |   11.081 |    0.842 |        84492 |    3582.860 |     1.200 |
| 2011-03 |   22082 |     1769385 |  1767780 |             1550 |        1550 |   12.885 |    2.824 |        96746 |   11361.450 |     3.200 |
| 2011-04 |    4072 |      326197 |   325167 |             1473 |        1473 |   11.883 |    0.590 |       141809 |    2384.690 |     1.100 |
| 2011-05 |   48624 |     3903146 |  3901075 |             1550 |        1550 |   12.119 |    5.384 |       182580 |   22194.290 |     5.800 |
| 2011-06 |   31070 |     2471437 |  2470417 |             1499 |        1499 |   11.857 |    3.935 |       198646 |   15995.640 |     4.400 |
| 2011-07 |   16690 |     1289627 |  1288534 |             1544 |        1544 |   12.241 |    2.035 |       279717 |    8170.390 |     2.900 |
| 2011-08 |    5423 |      433194 |   432153 |             1553 |        1553 |   12.658 |    1.095 |       308677 |    4229.510 |     2.100 |
| 2011-09 |    4357 |      360627 |   359450 |             1530 |        1530 |   13.196 |    0.923 |       252700 |    3486.110 |     1.700 |
| 2011-10 |    3400 |      278795 |   276840 |             1549 |        1549 |   14.090 |    0.671 |       167050 |    2436.990 |     1.200 |
| 2011-11 |    2387 |      190346 |   188299 |             1482 |        1482 |   14.497 |    0.446 |        96262 |    1486.370 |     0.800 |
| 2011-12 |    3253 |      257685 |   255003 |             1547 |        1547 |   15.423 |    0.602 |       129031 |    1987.720 |     1.000 |
| 2012-01 |    2788 |      223679 |   221224 |             1524 |        1524 |   15.521 |    0.526 |       113652 |    1751.420 |     0.900 |
| 2012-02 |    2847 |      229870 |   227569 |             1441 |        1441 |   14.773 |    0.554 |       120390 |    1805.500 |     0.900 |
| 2012-03 |    3614 |      289959 |   287137 |             1549 |        1549 |   16.211 |    0.678 |       145535 |    2236.230 |     1.100 |
| 2012-04 |    4407 |      348838 |   346497 |             1489 |        1489 |   14.669 |    0.903 |       206955 |    3024.370 |     1.600 |
| 2012-05 |    5152 |      407140 |   404553 |             1526 |        1526 |   14.907 |    1.022 |       239719 |    3487.390 |     1.800 |
| 2012-06 |    5771 |      455853 |   453084 |             1470 |        1470 |   14.560 |    1.202 |       284896 |    4075.500 |     2.100 |
| 2012-07 |    5788 |      456291 |   453538 |             1523 |        1523 |   14.562 |    1.213 |       298225 |    4276.310 |     2.100 |

Station-day read guard: `duc_fin/status/GUARD_cut.txt`.

## 4. Confidence

Day roles fixed before fitting (duc_fin/day_table.csv): held-out months 2008-02, 2008-09, 2009-03, 2009-10, 2010-05, 2010-09; PNET training days = every 10th 2008-2010 day (never fitted / evaluated); PNET held-out stations PB.B935, NC.KHMB, NC.LBP (NC.LBP is 1-component, fed as Z, 0, 0).

AUC (selected columns; all in eval_auc.csv):

| eval_set   | label    |      n |   n_pos |   base_rate |   n_days |   auc_FULL |   ap_FULL |   auc_NO_PNET |   auc_PNET_W_ONLY |   auc_NO_FAM |   auc_MF_RULE |   auc_xmad |   auc_pass_his |   auc_margin |   auc_site_z_max |   auc_pnet_max_w |   auc_pnet_max_f |
|:-----------|:---------|-------:|--------:|------------:|---------:|-----------:|----------:|--------------:|------------------:|-------------:|--------------:|-----------:|---------------:|-------------:|-----------------:|-----------------:|-----------------:|
| fit_oof    | same     | 103189 |   18067 |       0.175 |      823 |      0.998 |     0.989 |         0.998 |             0.998 |        0.998 |         0.997 |      0.664 |          0.979 |        0.982 |            0.355 |            0.759 |            0.804 |
| fit_oof    | any      | 103189 |   20590 |       0.200 |      823 |      0.994 |     0.978 |         0.994 |             0.994 |        0.994 |         0.987 |      0.639 |          0.922 |        0.947 |            0.361 |            0.753 |            0.816 |
| fit_oof    | same_raw | 103189 |   88797 |       0.861 |      823 |      0.842 |     0.959 |       nan     |           nan     |      nan     |       nan     |      0.808 |          0.567 |        0.673 |            0.632 |            0.527 |            0.588 |
| fit_oof    | any_raw  | 103189 |   89816 |       0.870 |      823 |      0.851 |     0.964 |       nan     |           nan     |      nan     |       nan     |      0.803 |          0.569 |        0.658 |            0.625 |            0.533 |            0.609 |
| ho_month   | same     |  20000 |    4511 |       0.226 |      182 |      0.997 |     0.991 |         0.997 |             0.997 |        0.997 |         0.997 |      0.689 |          0.979 |        0.982 |            0.433 |            0.742 |            0.802 |
| ho_month   | any      |  20000 |    5154 |       0.258 |      182 |      0.993 |     0.981 |         0.993 |             0.993 |        0.992 |         0.985 |      0.666 |          0.921 |        0.950 |            0.436 |            0.737 |            0.814 |
| ho_month   | same_raw |  20000 |   16538 |       0.827 |      182 |      0.898 |     0.976 |       nan     |           nan     |      nan     |       nan     |      0.681 |          0.601 |        0.662 |            0.686 |            0.613 |            0.667 |
| ho_month   | any_raw  |  20000 |   16765 |       0.838 |      182 |      0.909 |     0.981 |       nan     |           nan     |      nan     |       nan     |      0.669 |          0.603 |        0.651 |            0.684 |            0.625 |            0.687 |
| ind2011a   | same     |  41778 |    1943 |       0.047 |      120 |      0.999 |     0.985 |         0.999 |             0.999 |        0.999 |         0.999 |      0.149 |          0.980 |        0.976 |            0.097 |            0.950 |            0.938 |
| ind2011a   | any      |  41778 |    2329 |       0.056 |      120 |      0.997 |     0.964 |         0.997 |             0.997 |        0.997 |         0.992 |      0.138 |          0.905 |        0.921 |            0.093 |            0.948 |            0.939 |
| ind2011a   | same_raw |  41778 |   39976 |       0.957 |      120 |      0.960 |     0.998 |       nan     |           nan     |      nan     |       nan     |      0.784 |          0.481 |        0.575 |            0.781 |            0.307 |            0.477 |
| ind2011a   | any_raw  |  41778 |   40141 |       0.961 |      120 |      0.968 |     0.999 |       nan     |           nan     |      nan     |       nan     |      0.765 |          0.486 |        0.565 |            0.765 |            0.328 |            0.511 |
| gap2011    | same     | 115204 |    2749 |       0.024 |      245 |      0.990 |     0.618 |         0.990 |             0.990 |        0.990 |         0.991 |      0.250 |          0.975 |        0.906 |            0.209 |            0.928 |            0.947 |
| gap2011    | any      | 115204 |    3515 |       0.031 |      245 |      0.984 |     0.611 |         0.983 |             0.984 |        0.985 |         0.972 |      0.238 |          0.908 |        0.852 |            0.217 |            0.913 |            0.938 |
| gap2011    | same_raw | 115204 |  103639 |       0.900 |      245 |      0.852 |     0.974 |       nan     |           nan     |      nan     |       nan     |      0.853 |          0.404 |        0.701 |            0.808 |            0.232 |            0.341 |
| gap2011    | any_raw  | 115204 |  105289 |       0.914 |      245 |      0.832 |     0.974 |       nan     |           nan     |      nan     |       nan     |      0.837 |          0.414 |        0.635 |            0.787 |            0.256 |            0.381 |

Calibration and scrambled control (FULL, NO_PNET):

| eval_set   | label    | model   |      n |   slope |   icpt |       sumP |   matched |   obs_over_exp |   brier |   data_fracP_gt05 |   ctl_fracP_gt05 |   data_meanP |   ctl_meanP |   n_ctl |   ctl_sumP |   ctl_over_data_sumP |
|:-----------|:---------|:--------|-------:|--------:|-------:|-----------:|----------:|---------------:|--------:|------------------:|-----------------:|-------------:|------------:|--------:|-----------:|---------------------:|
| fit_oof    | same     | FULL    | 103189 |   1.002 | -0.000 |  18072.026 |     18067 |          1.000 |   0.013 |             0.180 |            0.016 |        0.175 |       0.015 |  140649 |   2063.146 |                0.114 |
| fit_oof    | same     | NO_PNET | 103189 |   1.002 | -0.000 |  18068.202 |     18067 |          1.000 |   0.013 |             0.179 |            0.016 |        0.175 |       0.015 |  140649 |   2089.708 |                0.116 |
| fit_oof    | any      | FULL    | 103189 |   1.003 | -0.001 |  20588.338 |     20590 |          1.000 |   0.024 |             0.199 |            0.017 |        0.200 |       0.020 |  140649 |   2747.749 |                0.133 |
| fit_oof    | any      | NO_PNET | 103189 |   1.003 | -0.001 |  20592.960 |     20590 |          1.000 |   0.024 |             0.199 |            0.019 |        0.200 |       0.021 |  140649 |   2953.446 |                0.143 |
| fit_oof    | same_raw | FULL    | 103189 |   1.021 | -0.017 |  88662.724 |     88797 |          1.002 |   0.092 |             0.950 |            0.611 |        0.859 |       0.679 |  140649 |  95474.656 |                1.077 |
| fit_oof    | any_raw  | FULL    | 103189 |   1.023 | -0.019 |  89705.534 |     89816 |          1.001 |   0.086 |             0.952 |            0.621 |        0.869 |       0.672 |  140649 |  94563.406 |                1.054 |
| ho_month   | same     | FULL    |  20000 |   1.002 | -0.000 |   4510.975 |      4511 |          1.000 |   0.015 |             0.232 |            0.017 |        0.226 |       0.016 |   32484 |    524.237 |                0.116 |
| ho_month   | same     | NO_PNET |  20000 |   1.002 | -0.001 |   4513.314 |      4511 |          0.999 |   0.015 |             0.232 |            0.017 |        0.226 |       0.016 |   32484 |    528.654 |                0.117 |
| ho_month   | any      | FULL    |  20000 |   1.003 | -0.001 |   5165.150 |      5154 |          0.998 |   0.029 |             0.258 |            0.019 |        0.258 |       0.023 |   32484 |    751.891 |                0.146 |
| ho_month   | any      | NO_PNET |  20000 |   1.004 | -0.001 |   5157.884 |      5154 |          0.999 |   0.029 |             0.257 |            0.020 |        0.258 |       0.025 |   32484 |    821.328 |                0.159 |
| ho_month   | same_raw | FULL    |  20000 |   1.193 | -0.161 |  16565.691 |     16538 |          0.998 |   0.091 |             0.907 |            0.401 |        0.828 |       0.528 |   32484 |  17152.044 |                1.035 |
| ho_month   | any_raw  | FULL    |  20000 |   1.236 | -0.206 |  16896.361 |     16765 |          0.992 |   0.086 |             0.910 |            0.417 |        0.845 |       0.608 |   32484 |  19748.124 |                1.169 |
| ind2011a   | same     | FULL    |  41778 |   0.990 | -0.001 |   1988.943 |      1943 |          0.977 |   0.004 |             0.048 |            0.007 |        0.048 |       0.006 |   31465 |    194.115 |                0.098 |
| ind2011a   | same     | NO_PNET |  41778 |   0.989 | -0.001 |   1987.520 |      1943 |          0.978 |   0.004 |             0.049 |            0.007 |        0.048 |       0.006 |   31465 |    197.258 |                0.099 |
| ind2011a   | any      | FULL    |  41778 |   0.984 | -0.001 |   2410.492 |      2329 |          0.966 |   0.009 |             0.056 |            0.007 |        0.058 |       0.017 |   31465 |    541.811 |                0.225 |
| ind2011a   | any      | NO_PNET |  41778 |   0.985 | -0.001 |   2408.066 |      2329 |          0.967 |   0.009 |             0.055 |            0.007 |        0.058 |       0.022 |   31465 |    679.462 |                0.282 |
| ind2011a   | same_raw | FULL    |  41778 |   1.069 | -0.063 |  39837.971 |     39976 |          1.003 |   0.024 |             0.977 |            0.639 |        0.954 |       0.683 |   31465 |  21491.707 |                0.539 |
| ind2011a   | any_raw  | FULL    |  41778 |   1.114 | -0.099 |  39751.788 |     40141 |          1.010 |   0.022 |             0.977 |            0.696 |        0.952 |       0.730 |   31465 |  22967.258 |                0.578 |
| gap2011    | same     | FULL    | 115204 |   0.452 | -0.001 |   6249.575 |      2749 |          0.440 |   0.029 |             0.056 |            0.006 |        0.054 |       0.006 |  104707 |    588.819 |                0.094 |
| gap2011    | same     | NO_PNET | 115204 |   0.451 | -0.001 |   6252.595 |      2749 |          0.440 |   0.029 |             0.056 |            0.006 |        0.054 |       0.006 |  104707 |    594.981 |                0.095 |
| gap2011    | any      | FULL    | 115204 |   0.494 |  0.000 |   7033.818 |      3515 |          0.500 |   0.031 |             0.061 |            0.007 |        0.061 |       0.012 |  104707 |   1251.231 |                0.178 |
| gap2011    | any      | NO_PNET | 115204 |   0.495 |  0.000 |   7044.511 |      3515 |          0.499 |   0.031 |             0.061 |            0.007 |        0.061 |       0.015 |  104707 |   1520.732 |                0.216 |
| gap2011    | same_raw | FULL    | 115204 |   0.854 |  0.096 | 108454.609 |    103639 |          0.956 |   0.078 |             0.974 |            0.700 |        0.941 |       0.727 |  104707 |  76103.801 |                0.702 |
| gap2011    | any_raw  | FULL    | 115204 |   0.910 |  0.058 | 108452.260 |    105289 |          0.971 |   0.067 |             0.978 |            0.740 |        0.941 |       0.771 |  104707 |  80684.312 |                0.744 |

Within one family-day:

| eval_set   | label   | score   |   n_famdays |   auc_within_famday |
|:-----------|:--------|:--------|------------:|--------------------:|
| fit_oof    | same    | FULL    |        1353 |               0.971 |
| fit_oof    | same    | xmad    |        1353 |               0.928 |
| fit_oof    | any     | FULL    |        2095 |               0.957 |
| fit_oof    | any     | xmad    |        2095 |               0.789 |
| ho_month   | same    | FULL    |         294 |               0.977 |
| ho_month   | same    | xmad    |         294 |               0.918 |
| ho_month   | any     | FULL    |         494 |               0.950 |
| ho_month   | any     | xmad    |         494 |               0.784 |
| ind2011a   | same    | FULL    |         154 |               0.979 |
| ind2011a   | same    | xmad    |         154 |               0.917 |
| ind2011a   | any     | FULL    |         265 |               0.892 |
| ind2011a   | any     | xmad    |         265 |               0.591 |
| gap2011    | same    | FULL    |         313 |               0.866 |
| gap2011    | same    | xmad    |         313 |               0.878 |
| gap2011    | any     | FULL    |         553 |               0.921 |
| gap2011    | any     | xmad    |         553 |               0.524 |

Day-block bootstrap of AUC(FULL) - AUC(other):

| eval_set   | label   | vs       |   n_days |   dAUC_mean |   dAUC_lo |   dAUC_hi |
|:-----------|:--------|:---------|---------:|------------:|----------:|----------:|
| ho_month   | same    | xmad     |      182 |       0.307 |     0.273 |     0.354 |
| ho_month   | same    | pass_his |      182 |       0.019 |     0.015 |     0.022 |
| ho_month   | same    | NO_PNET  |      182 |      -0.000 |    -0.000 |     0.000 |
| ho_month   | same    | MF_RULE  |      182 |       0.000 |     0.000 |     0.001 |
| ho_month   | any     | xmad     |      182 |       0.331 |     0.296 |     0.377 |
| ho_month   | any     | pass_his |      182 |       0.073 |     0.059 |     0.088 |
| ho_month   | any     | NO_PNET  |      182 |      -0.000 |    -0.000 |     0.000 |
| ho_month   | any     | MF_RULE  |      182 |       0.008 |     0.006 |     0.011 |
| ind2011a   | same    | xmad     |      120 |       0.848 |     0.804 |     0.888 |
| ind2011a   | same    | pass_his |      120 |       0.020 |     0.013 |     0.028 |
| ind2011a   | same    | NO_PNET  |      120 |       0.000 |     0.000 |     0.000 |
| ind2011a   | same    | MF_RULE  |      120 |       0.000 |     0.000 |     0.000 |
| ind2011a   | any     | xmad     |      120 |       0.857 |     0.800 |     0.897 |
| ind2011a   | any     | pass_his |      120 |       0.093 |     0.071 |     0.129 |
| ind2011a   | any     | NO_PNET  |      120 |       0.000 |     0.000 |     0.000 |
| ind2011a   | any     | MF_RULE  |      120 |       0.005 |     0.003 |     0.008 |

Permutation importance (independent rows; blocks and top single features):

| label   | feature               | is_block   |   auc_drop |   logloss_gain |     n |
|:--------|:----------------------|:-----------|-----------:|---------------:|------:|
| any     | his-rule flags        | True       |      0.101 |          0.285 | 61778 |
| any     | density / competition | True       |      0.015 |          0.050 | 61778 |
| any     | family id             | True       |      0.005 |          0.037 | 61778 |
| any     | PNET (all)            | True       |      0.000 |          0.002 | 61778 |
| any     | PNET fine-tuned       | True       |      0.000 |          0.002 | 61778 |
| any     | site cc / z           | True       |      0.000 |          0.000 | 61778 |
| any     | PNET off-the-shelf    | True       |      0.000 |         -0.000 | 61778 |
| same    | his-rule flags        | True       |      0.230 |          1.160 | 61778 |
| same    | density / competition | True       |      0.000 |          0.003 | 61778 |
| same    | family id             | True       |      0.000 |          0.001 | 61778 |
| same    | site cc / z           | True       |      0.000 |          0.001 | 61778 |
| same    | PNET (all)            | True       |      0.000 |          0.000 | 61778 |
| same    | PNET off-the-shelf    | True       |      0.000 |          0.000 | 61778 |
| same    | PNET fine-tuned       | True       |      0.000 |          0.000 | 61778 |

| label   | feature      | is_block   |   auc_drop |   logloss_gain |     n |
|:--------|:-------------|:-----------|-----------:|---------------:|------:|
| same    | margin       | False      |      0.195 |          1.031 | 61778 |
| any     | margin       | False      |      0.088 |          0.241 | 61778 |
| any     | n_oth_4s     | False      |      0.006 |          0.022 | 61778 |
| any     | fam_code     | False      |      0.005 |          0.039 | 61778 |
| same    | xmad         | False      |      0.001 |          0.013 | 61778 |
| any     | win_cross_2s | False      |      0.001 |          0.007 | 61778 |
| any     | xmad         | False      |      0.001 |          0.008 | 61778 |
| any     | ccsum        | False      |      0.000 |          0.003 | 61778 |
| same    | win_cross_2s | False      |      0.000 |          0.001 | 61778 |
| same    | fam_code     | False      |      0.000 |          0.001 | 61778 |
| any     | pnet_2nd_f   | False      |      0.000 |          0.001 | 61778 |
| same    | dens_own_10m | False      |      0.000 |          0.001 | 61778 |
| same    | ccsum        | False      |      0.000 |          0.001 | 61778 |
| any     | site_cc_2nd  | False      |      0.000 |          0.000 | 61778 |
| any     | pass_his     | False      |      0.000 |          0.001 | 61778 |

PNET per station, before (off-the-shelf) / after (fine-tuned) -- AUC of the single-station value:

| netsta   | held_out_of_finetune   |   auc_any_fine-tuned |   auc_any_off-the-shelf |   auc_same_fine-tuned |   auc_same_off-the-shelf |
|:---------|:-----------------------|---------------------:|------------------------:|----------------------:|-------------------------:|
| BK.GASB  | False                  |                0.622 |                   0.563 |                 0.630 |                    0.586 |
| BK.HOPS  | False                  |                0.582 |                   0.515 |                 0.550 |                    0.509 |
| BK.JCC   | False                  |                0.575 |                   0.586 |                 0.542 |                    0.572 |
| BK.WDC   | False                  |                0.660 |                   0.658 |                 0.650 |                    0.670 |
| BK.YBH   | False                  |                0.881 |                   0.586 |                 0.875 |                    0.643 |
| NC.GBB   | False                  |                0.702 |                   0.454 |                 0.703 |                    0.429 |
| NC.GHM   | False                  |                0.520 |                   0.510 |                 0.509 |                    0.506 |
| NC.GRO   | False                  |                0.490 |                   0.388 |                 0.477 |                    0.368 |
| NC.GSN   | False                  |                0.529 |                   0.466 |                 0.507 |                    0.474 |
| NC.GTC   | False                  |                0.926 |                   0.728 |                 0.924 |                    0.734 |
| NC.GVV   | False                  |                0.490 |                   0.504 |               nan     |                  nan     |
| NC.KBN   | False                  |                0.536 |                   0.497 |                 0.531 |                    0.501 |
| NC.KBS   | False                  |                0.550 |                   0.536 |                 0.519 |                    0.512 |
| NC.KCPB  | False                  |                0.885 |                   0.926 |                 0.874 |                    0.923 |
| NC.KCS   | False                  |                0.560 |                   0.529 |                 0.531 |                    0.513 |
| NC.KFP   | False                  |                0.871 |                   0.779 |                 0.867 |                    0.784 |
| NC.KHBB  | False                  |                0.711 |                   0.528 |                 0.693 |                    0.533 |
| NC.KHMB  | True                   |                0.623 |                   0.560 |                 0.594 |                    0.573 |
| NC.KIP   | False                  |                0.859 |                   0.808 |                 0.848 |                    0.810 |
| NC.KKP   | False                  |                0.619 |                   0.515 |                 0.603 |                    0.501 |
| NC.KOM   | False                  |                0.638 |                   0.502 |                 0.620 |                    0.475 |
| NC.KPP   | False                  |                0.598 |                   0.554 |                 0.570 |                    0.542 |
| NC.KRK   | False                  |                0.817 |                   0.694 |                 0.811 |                    0.703 |
| NC.KRMB  | False                  |                0.694 |                   0.679 |                 0.699 |                    0.690 |
| NC.KRP   | False                  |                0.511 |                   0.462 |                 0.480 |                    0.459 |
| NC.KSXB  | False                  |                0.868 |                   0.550 |                 0.864 |                    0.574 |
| NC.KTR   | False                  |                0.656 |                   0.568 |                 0.651 |                    0.551 |
| NC.LAM   | False                  |                0.623 |                   0.588 |                 0.607 |                    0.563 |
| NC.LBF   | False                  |                0.519 |                   0.419 |                 0.533 |                    0.424 |
| NC.LBK   | False                  |                0.601 |                   0.374 |                 0.592 |                    0.358 |
| NC.LBP   | True                   |                0.605 |                   0.440 |                 0.574 |                    0.403 |
| NC.LDB   | False                  |                0.670 |                   0.500 |                 0.676 |                    0.482 |
| NC.LGB   | False                  |                0.470 |                   0.437 |                 0.423 |                    0.399 |
| NC.LGP   | False                  |                0.469 |                   0.436 |                 0.460 |                    0.429 |
| NC.LPG   | False                  |                0.632 |                   0.388 |                 0.607 |                    0.364 |
| NC.LRB   | False                  |                0.612 |                   0.405 |                 0.586 |                    0.383 |
| NC.LRR   | False                  |                0.603 |                   0.520 |                 0.586 |                    0.543 |
| NC.LSF   | False                  |                0.625 |                   0.350 |                 0.605 |                    0.339 |
| NC.LSH   | False                  |                0.589 |                   0.367 |                 0.566 |                    0.334 |
| NC.LSR   | False                  |                0.496 |                   0.367 |                 0.478 |                    0.352 |
| NC.LTC   | False                  |                0.537 |                   0.520 |                 0.528 |                    0.508 |
| NC.LVR   | False                  |                0.600 |                   0.494 |                 0.611 |                    0.495 |
| NC.LWH   | False                  |                0.504 |                   0.455 |                 0.483 |                    0.426 |
| PB.B039  | False                  |                0.651 |                   0.719 |                 0.643 |                    0.739 |
| PB.B040  | False                  |                0.709 |                   0.670 |                 0.690 |                    0.703 |
| PB.B933  | False                  |                0.535 |                   0.471 |                 0.535 |                    0.471 |
| PB.B935  | True                   |                0.563 |                   0.612 |                 0.539 |                    0.627 |

![pnet learning curve](duc_fin/figures/pnet_learning_curve.png)

## 5. Delete candidates

```
DELETE-CANDIDATE list for the Ducellier stations 2008-01-01..2012-07-31 (NOTHING DELETED) -- 2026-10-01 21:51 UTC
candidates: 138942 files, 1.315 TB -> /home/jovyan/tremorferometry/lfe_south/results/cut/DELETE_CANDIDATES_ducellier.csv
kept: 31095 files, 0.344 TB -> KEEP_ducellier.csv (summary KEEP_ducellier_summary.csv)
DUC_FIN stage markers: {'MF_DONE': True, 'HF_DONE': True, 'PNET_DONE': True, 'CUT_DONE': True, 'FIT_DONE': True}
days of 2008-2012 used by Mendocino discovery / Shelly sampled days: 0
by kind (cand = candidate):
          kind  cand  size          sum    TB
     cache20hz  True 14630 117103976220 0.117
     cache_duc False 15767 187703880718 0.188
     cache_duc  True 61098 704104583296 0.704
  prep_duc_fix False   347     14999587 0.000
  prep_duc_fix  True  1327     58897360 0.000
    raw_traces False  2678  17875078656 0.018
    raw_traces  True 13355  72426188800 0.072
raw_traces_duc False 12303 138894818816 0.139
raw_traces_duc  True 48532 421363882496 0.421
keep reasons (files):
keep_reason
DUC_FIX full-record extension: final 24-month re-templating of the 2011-10..12 epochs + pass-2 re-scan of 2011-10..2012-07    28062
known-answer window (duc_fix stage 1)                                                                                          1932
sh_back BK2 (S6 2008-2017                                                                                                      1406
queue halted, not finished)                                                                                                    1406
scripts that reference lfe_south raw / cache folders (for the user to check):
  /home/jovyan/tremorferometry/lfe_array_discovery/scripts/ab_common.py
  /home/jovyan/tremorferometry/lfe_array_discovery/scripts/mendo_discovery_queue.sh
  /home/jovyan/tremorferometry/lfe_south/scripts/cut_windows.py
  /home/jovyan/tremorferometry/lfe_south/scripts/duc_fin_common.py
  /home/jovyan/tremorferometry/lfe_south/scripts/duc_fin_cut.py
  /home/jovyan/tremorferometry/lfe_south/scripts/duc_fin_delete.py
  /home/jovyan/tremorferometry/lfe_south/scripts/duc_fin_pnet.py
  /home/jovyan/tremorferometry/lfe_south/scripts/duc_fix_clean.py
  /home/jovyan/tremorferometry/lfe_south/scripts/duc_fix_common.py
  /home/jovyan/tremorferometry/lfe_south/scripts/duc_fix_data.py
  /home/jovyan/tremorferometry/lfe_south/scripts/duc_fix_prep.py
  /home/jovyan/tremorferometry/lfe_south/scripts/ks_common.py
  /home/jovyan/tremorferometry/lfe_south/scripts/t14_port_k.py
  /home/jovyan/tremorferometry/lfe_south/scripts/t6_duc_cache.py
  /home/jovyan/tremorferometry/lfe_south/scripts/t7_duc_mf.py
  /home/jovyan/tremorferometry/lfe_south/scripts/t7b_bedB_templates.py
  /home/jovyan/tremorferometry/lfe_south/scripts/t_south_queue.sh
```

## Station-day read guards

* [guard:cut] station-days ok 82935, error/missing 0 (0.00%); all {'ok': 82935, 'nodata': 946}
* [guard:hf] station-days ok 14518, error/missing 0 (0.00%); all {'ok': 14518}
* [guard:mf] station-days ok 76853, error/missing 0 (0.00%); all {'ok': 76853, 'nodata': 949, 'short': 12}
* [guard:pnet] station-days ok 76853, error/missing 0 (0.00%); all {'ok': 76853}
* [guard:pnetlabels] station-days ok 3922, error/missing 0 (0.00%); all {'ok': 3922, 'nodata': 82}
* [guard:retpl] station-days ok 6649, error/missing 0 (0.00%); all {'ok': 6649, 'nodata': 18}
* [mfcheck] days 1674, feature files missing 0; days without new templates compared with mf1 1395: identical 1358, different 37; max |cc diff| 3.72e-15

