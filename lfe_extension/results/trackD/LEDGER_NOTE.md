# Track C ledger vs the catalog — 2026-09-17 11:21 UTC

`n_cat_ext` / `n_cat_core` in `results/trackC/months_done.csv` are STALE for the months below: those rows
were written by `c_watch.mark_done` BEFORE Track C re-derived the tiers per epoch and re-flagged
`cat_ext` / `cat_core`. Each month's parquet was verified to already carry `cat_*` columns matching the
CURRENT per-epoch thresholds (/home/jovyan/tremorferometry/lfe_extension/results/trackC/tiers_by_epoch.csv (Track C, FA-only rule; 4 epochs)), so the FILES are right and only the ledger is behind — and
`results/trackD/lfe_extension_catalog_v1.parquet`, which is built from the parquets, is correct.

| ym | epochs | ledger_cat_ext | ledger_cat_core | file_cat_ext | file_cat_core | file_matches_current_thresholds |
|---|---|---|---|---|---|---|
| 2014-11 | E0 | 995.0 | 579.0 | 976 | 570 | True |
| 2014-12 | E0 | 1618.0 | 880.0 | 1264 | 695 | True |
| 2017-02 | E0 | 434.0 | 164.0 | 426 | 162 | True |
| 2017-03 | E1 | 3035.0 | 2294.0 | 3035 | 1232 | True |
| 2017-04 | E1 | 923.0 | 512.0 | 923 | 130 | True |
| 2017-05 | E1 | 2154.0 | 1324.0 | 2154 | 352 | True |
| 2017-06 | E1 | 3274.0 | 2027.0 | 3274 | 618 | True |
| 2017-07 | E1 | 1248.0 | 660.0 | 1248 | 106 | True |
| 2017-08 | E1 | 1065.0 | 572.0 | 1065 | 114 | True |
| 2017-09 | E1 | 1608.0 | 1039.0 | 1608 | 322 | True |
| 2017-10 | E1 | 2379.0 | 1809.0 | 2379 | 918 | True |
| 2017-11 | E1 | 433.0 | 234.0 | 433 | 45 | True |
| 2017-12 | E1 | 418.0 | 232.0 | 418 | 51 | True |
| 2018-01 | E1 | 1092.0 | 756.0 | 1092 | 365 | True |
| 2018-02 | E1 | 1022.0 | 761.0 | 1022 | 359 | True |
| 2018-03 | E1 | 850.0 | 521.0 | 835 | 153 | True |
| 2018-04 | E1 | 508.0 | 289.0 | 499 | 66 | True |
| 2018-05 | E1 | 1965.0 | 1302.0 | 1965 | 442 | True |
| 2018-06 | E1 | 2305.0 | 1675.0 | 2305 | 796 | True |
| 2018-07 | E1 | 6897.0 | 5371.0 | 6897 | 2903 | True |
| 2018-08 | E1 | 1245.0 | 787.0 | 1245 | 281 | True |
| 2018-09 | E1 | 1201.0 | 716.0 | 1201 | 187 | True |
| 2018-10 | E1 | 2329.0 | 1610.0 | 2329 | 619 | True |
| 2018-11 | E1 | 1379.0 | 803.0 | 1379 | 190 | True |
| 2018-12 | E1 | 1028.0 | 782.0 | 1028 | 411 | True |
| 2019-01 | E2 | 604.0 | 175.0 | 604 | 110 | True |
| 2019-02 | E2 | 295.0 | 65.0 | 291 | 41 | True |
| 2019-03 | E2 | 693.0 | 162.0 | 688 | 110 | True |
| 2019-04 | E2 | 817.0 | 263.0 | 812 | 192 | True |
| 2019-05 | E2 | 1741.0 | 720.0 | 1728 | 551 | True |
| 2019-06 | E2 | 1019.0 | 311.0 | 1019 | 247 | True |
| 2019-07 | E2 | 1743.0 | 654.0 | 1743 | 507 | True |
| 2019-08 | E2 | 5607.0 | 2808.0 | 5607 | 2321 | True |
| 2019-09 | E2 | 1131.0 | 299.0 | 1131 | 193 | True |
| 2019-10 | E2 | 716.0 | 177.0 | 716 | 114 | True |
| 2019-11 | E2 | 4792.0 | 2359.0 | 4792 | 1921 | True |
| 2019-12 | E2 | 408.0 | 67.0 | 408 | 40 | True |
| 2020-01 | E2 | 247.0 | 40.0 | 247 | 25 | True |
| 2020-02 | E2 | 736.0 | 308.0 | 736 | 261 | True |
| 2020-03 | E2 | 1123.0 | 570.0 | 1123 | 456 | True |
| 2020-04 | E2 | 999.0 | 273.0 | 999 | 181 | True |
| 2020-05 | E2 | 1464.0 | 517.0 | 1464 | 388 | True |
| 2020-06 | E2 | 906.0 | 240.0 | 906 | 163 | True |
| 2020-07 | E2 | 1134.0 | 328.0 | 1134 | 232 | True |
| 2020-08 | E2 | 1004.0 | 277.0 | 1004 | 194 | True |
| 2020-09 | E2 | 2968.0 | 1413.0 | 2968 | 1156 | True |
| 2020-10 | E2 | 858.0 | 307.0 | 858 | 230 | True |
| 2020-11 | E2 | 3220.0 | 1687.0 | 3220 | 1407 | True |
| 2020-12 | E2 | 418.0 | 121.0 | 418 | 84 | True |
| 2021-01 | E2 | 357.0 | 102.0 | 357 | 78 | True |
| 2021-02 | E2 | 2571.0 | 1261.0 | 2571 | 1004 | True |
| 2021-03 | E2 | 549.0 | 145.0 | 549 | 101 | True |
| 2024-05 | E3 | 1292.0 | 189.0 | 1284 | 182 | True |

Produced by `lfe_extension/scripts/d_fix_ledger.py`. Run it with `--apply` to update those two fields in
place in the ledger (it backs the file up to `months_done.csv.pre_d_fix_ledger` first and refuses unless
every disagreeing month's parquet already matches the current thresholds).
