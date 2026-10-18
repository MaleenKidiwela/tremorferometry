# Washington LFE families -- the full ranked tables, real against the site-shift null

Agent S, 2026-09-23.  Workspace `results/sf_final/`; queue `scripts/s_final_queue.sh`; code
`scripts/sf_*.py`.

## 0. What changed, and why there is no candidate count in this document

The user's objection, 2026-09-23: *"you are discarding them, what if there are undiscovered
families"*.  It is correct and it is fatal to the design every previous Washington write-up used.
Calibrating an acceptance criterion against the 58 blind-rediscovered Bostock families can only
detect a criterion that is too strict **for that population**.  It offers no protection at all to
a real family that is smaller, deeper, noisier, or seen by fewer stations -- which is exactly what
an undiscovered Washington family would be.  Every "0 candidates" verdict in
`WA_2010_RESULTS.md`, `WA_2011_REFIRE_RESULTS.md`, `WA_2008_RESULTS.md`, `WA_2008Q_RESULTS.md` and
`WA_2010_ARRAY_RESULTS.md` was produced by exactly that mistake, compounded by clauses that have
no counterpart in any published Cascadia catalogue.

The rule now in force:

1. **Nothing is filtered.**  Every grown family keeps its row, permanently.
2. **Every criterion is a column**, carrying its measured value -- never a gate.
3. **The whole-chain site-shift null gets the byte-identical table.**  It is the only honest
   comparator.
4. **The evidence is the real-vs-null count as a function of each criterion** -- a separation
   curve -- plus a ranking by the strongest discriminators, not a pre-chosen threshold.
5. The >= 85 %-of-known-families test survives only as a **veto on a proposed cut**.  Passing it
   never licenses a cut.
6. The only steps allowed to reduce a list are the two the papers actually perform -- the
   **dedupe** (>10 % shared declustered detections, or stack cc > 0.3) and the **location** -- and
   both are reported with counts before and after, with the merged and unlocatable rows kept and
   flagged rather than deleted.

## 1. Inventory: what is in the tables

| bed | chain | grown templates (ALL kept) | dedupe: kept as distinct | merged (row retained) | located >= 5 stations | unlocatable (row retained) | on a grid face | on a depth rail | re-fired on the sealed year |
|---|---|---|---|---|---|---|---|---|---|
| 2010 N-Olympic, arrays as arrays (agent O) | real | 1485 | 780 | 705 | 971 | 514 | 858 | 800 | 1425 |
| 2010 N-Olympic, arrays as arrays (agent O) | null | 1486 | 1247 | 239 | 761 | 725 | 694 | 639 | 1391 |
| 2008 Puget / S-Olympic (agent N) | real | 160 | 129 | 31 | 156 | 4 | 128 | 122 | 78 |
| 2008 Puget / S-Olympic (agent N) | null | 36 | 35 | 1 | 36 | 0 | 32 | 31 | 12 |
| 2010 Arm B, Chestler's own chain (agent O) | real | 739 | 342 | 397 | 737 | 2 | 577 | 560 | 738 |
| 2010 Arm B, Chestler's own chain (agent O) | null | 50 | 50 | 0 | 50 | 0 | 46 | 44 | 50 |

Both reducing steps are reversible from the tables: `dedupe_kept` / `dedupe_merged_into` /
`dedupe_max_share` / `dedupe_max_cc` record what the dedupe did and why, and `located` / `n_sta`
record the location outcome.  No row was deleted.

## 2. Real against null, criterion by criterion

AUC is the probability that a randomly drawn real template scores better than a randomly drawn
null one.  **0.5 means the two populations are indistinguishable on that column.**  Lift is the
enrichment at the strictest threshold that still leaves at least 5 null templates.

### 2010 N-Olympic, arrays as arrays (agent O)

| criterion | direction | AUC real vs null | +- | separation \|AUC-0.5\| | favours | best lift | at threshold | n real | n null |
|---|---|---|---|---|---|---|---|---|---|
| `dedupe_kept` | flag | -- | -- | **0.314** | nan | 0.63 | -- | 780 | 1247 |
| `dedupe_n_partners` | lo | 0.242 | 0.011 | **0.258** | nan | 0.98 | 166.000 | 1456 | 1486 |
| `rms` | lo | 0.642 | 0.017 | **0.142** | real | 4.77 | 2.273 | 134 | 22 |
| `located` | flag | -- | -- | **0.142** | real | 1.28 | -- | 971 | 761 |
| `n_days` | hi | 0.638 | 0.012 | **0.138** | real | 8.49 | 20.000 | 280 | 33 |
| `rate_tremor` | hi | 0.631 | 0.013 | **0.131** | real | 8.75 | 2.552 | 271 | 31 |
| `n_det` | hi | 0.631 | 0.013 | **0.131** | real | 8.75 | 74.000 | 271 | 31 |
| `ctrl_margin` | hi | 0.621 | 0.018 | **0.121** | real | 3.76 | 1.064 | 144 | 30 |
| `edge` | flag | -- | -- | **0.111** | real | 1.24 | -- | 858 | 694 |
| `dep_rail` | flag | -- | -- | **0.109** | real | 1.25 | -- | 800 | 639 |
| `med_cc` | hi | 0.392 | 0.013 | **0.108** | nan | 1.00 | 0.075 | 1485 | 1485 |
| `n_sta` | hi | 0.607 | 0.013 | **0.107** | real | 3.72 | 9.000 | 253 | 68 |
| `sealed_rate` | hi | 0.594 | 0.013 | **0.094** | real | 2.06 | 0.769 | 202 | 98 |
| `sealed_n_det` | hi | 0.594 | 0.013 | **0.094** | real | 2.06 | 20.000 | 202 | 98 |
| `ctrl_pass` | flag | -- | -- | **0.093** | real | 3.88 | -- | 186 | 48 |
| `depth_window_15_58` | flag | -- | -- | **0.089** | real | 1.24 | -- | 689 | 557 |
| `co_timing` | hi | 0.413 | 0.013 | **0.087** | nan | 1.00 | 0.000 | 1485 | 1485 |
| `sealed_n_days` | hi | 0.583 | 0.013 | **0.083** | real | 6.60 | 25.000 | 33 | 5 |

Separation curve for `rms`:

| threshold | real beyond it | null beyond it | % of real | % of null | lift |
|---|---|---|---|---|---|
| 4.935 | 434 | 190 | 44.7 | 25.0 | 1.79 |
| 5.545 | 524 | 255 | 54.0 | 33.5 | 1.61 |
| 6.107 | 610 | 325 | 62.8 | 42.7 | 1.47 |
| 6.730 | 693 | 398 | 71.4 | 52.3 | 1.36 |
| 7.277 | 764 | 483 | 78.7 | 63.5 | 1.24 |
| 8.018 | 836 | 567 | 86.1 | 74.5 | 1.16 |
| 8.803 | 900 | 658 | 92.7 | 86.5 | 1.07 |
| 11.132 | 964 | 750 | 99.3 | 98.6 | 1.01 |

Reverse separations on this bed (the NULL population scores better than the real one): `dedupe_n_partners` AUC 0.242, `med_cc` AUC 0.392.  A criterion on which noise beats signal is evidence about the chain, not a candidate filter.

**2010 N-Olympic, arrays as arrays (agent O)** -- the 25 highest-ranking REAL rows on `rms`, the strongest discriminator for this bed.  `null %ile` is the fraction of the NULL population this row beats on that column; a row near 1.00 is unusual for the null, a row near 0.50 is typical of it.

| grp | rms | null %ile | n_det | n_days | lat | lon | dep | edge | dep_rail | ctrl_pass | near_chestler_km | sealed_n_det | sealed_n_days | sealed_fa |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 642 | 0.056 | 1.000 | 18 | 6 | 48.1 | -123.2 | 21.0 | False | False | True | 20.5 | 12 | 6 | 0.000 |
| 155 | 0.135 | 1.000 | 13 | 8 | 47.8 | -122.7 | 37.0 | False | False | True | 10.2 | 10 | 6 | 1.000 |
| 118 | 0.160 | 1.000 | 26 | 8 | 47.6 | -122.7 | 21.0 | False | False | True | 26.3 | 8 | 5 | 0.000 |
| 351 | 0.327 | 1.000 | 12 | 2 | 48.2 | -125.0 | 15.0 | True | True | True | 122.8 | 5 | 3 | 0.143 |
| 410 | 0.328 | 1.000 | 7 | 5 | 47.6 | -123.2 | 15.0 | True | True | True | 13.8 | 4 | 4 | 0.286 |
| 942 | 0.333 | 1.000 | 57 | 23 | 47.7 | -123.4 | 59.0 | True | True | False | 35.0 | 1 | 1 | 0.143 |
| 64 | 0.354 | 1.000 | 32 | 10 | 46.5 | -123.1 | 59.0 | True | True | True | 138.4 | 42 | 11 | 0.714 |
| 561 | 0.356 | 1.000 | 28 | 7 | 47.9 | -123.2 | 43.0 | False | False | False | 8.544 | 9 | 7 | 0.286 |
| 435 | 0.366 | 0.999 | 8 | 6 | 47.9 | -123.8 | 59.0 | True | True | True | 37.0 | 11 | 9 | 0.000 |
| 475 | 0.388 | 0.999 | 20 | 8 | 47.9 | -123.3 | 25.0 | False | False | True | 9.859 | 4 | 4 | 0.000 |
| 944 | 0.407 | 0.999 | 9 | 6 | 48.0 | -123.6 | 57.0 | False | False | True | 26.4 | 0 | 0 | 0.429 |
| 870 | 0.414 | 0.999 | 16 | 3 | 47.8 | -123.6 | 31.0 | False | False | False | 26.4 | 5 | 4 | 0.714 |
| 615 | 0.425 | 0.999 | 101 | 22 | 47.9 | -123.0 | 17.0 | False | False | True | 19.3 | 4 | 4 | 0.429 |
| 1450 | 0.427 | 0.999 | 41 | 4 | 47.9 | -122.4 | 59.0 | True | True | True | 24.3 | 1 | 1 | 0.143 |
| 777 | 0.433 | 0.999 | 16 | 11 | 48.0 | -123.8 | 31.0 | False | False | True | 27.7 | 10 | 8 | 1.000 |
| 350 | 0.434 | 0.999 | 19 | 4 | 47.7 | -123.3 | 15.0 | True | True | True | 5.982 | 3 | 3 | 0.286 |
| 890 | 0.515 | 0.999 | 6 | 3 | 47.8 | -123.2 | 53.0 | False | False | False | 23.1 | 5 | 5 | 0.000 |
| 85 | 0.550 | 0.999 | 117 | 7 | 47.7 | -122.8 | 57.0 | False | False | True | 20.5 | 37 | 4 | 0.000 |
| 1452 | 0.590 | 0.999 | 26 | 5 | 48.0 | -122.2 | 59.0 | True | True | True | 37.9 | 4 | 4 | 0.000 |
| 700 | 0.609 | 0.999 | 5 | 3 | 47.9 | -123.4 | 15.0 | True | True | True | 21.0 | 4 | 2 | 0.286 |
| 451 | 0.644 | 0.999 | 5 | 3 | 47.8 | -123.5 | 27.0 | False | False | False | 17.1 | 5 | 3 | 0.143 |
| 955 | 0.650 | 0.999 | 20 | 10 | 48.8 | -124.0 | 15.0 | True | True | True | 78.0 | 5 | 3 | 0.286 |
| 122 | 0.651 | 0.999 | 749 | 29 | 48.0 | -123.0 | 15.0 | True | True | False | 23.3 | 30 | 15 | 1.000 |
| 566 | 0.661 | 0.999 | 24 | 7 | 48.1 | -123.5 | 49.0 | False | False | True | 19.4 | 8 | 3 | 0.286 |
| 844 | 0.704 | 0.999 | 8 | 4 | 47.9 | -123.7 | 47.0 | False | False | False | 23.3 | 7 | 4 | 0.000 |

### 2008 Puget / S-Olympic (agent N)

| criterion | direction | AUC real vs null | +- | separation \|AUC-0.5\| | favours | best lift | at threshold | n real | n null |
|---|---|---|---|---|---|---|---|---|---|
| `med_cc` | hi | 0.273 | 0.074 | **0.227** | nan | 1.00 | 0.065 | 160 | 36 |
| `n_det` | hi | 0.722 | 0.075 | **0.222** | real | 3.87 | 8.000 | 86 | 5 |
| `rate_tremor` | hi | 0.722 | 0.075 | **0.222** | real | 3.87 | 0.400 | 86 | 5 |
| `dedupe_kept` | flag | -- | -- | **0.166** | nan | 0.83 | -- | 129 | 35 |
| `tremor_quiet_ratio` | hi | 0.666 | 0.081 | **0.166** | real | 2.21 | 2.217 | 74 | 8 |
| `n_days` | hi | 0.645 | 0.080 | **0.145** | real | 2.25 | 5.000 | 60 | 6 |
| `dedupe_n_partners` | lo | 0.363 | 0.080 | **0.137** | nan | 0.99 | 9.050 | 158 | 36 |
| `max_nmad` | hi | 0.364 | 0.080 | **0.136** | nan | 1.20 | 64.144 | 32 | 6 |
| `ctrl_margin` | hi | 0.608 | 0.081 | **0.108** | real | 1.85 | 0.848 | 48 | 6 |
| `sealed_rate` | hi | 0.605 | 0.081 | **0.105** | real | 1.46 | 0.038 | 78 | 12 |
| `sealed_n_det` | hi | 0.605 | 0.081 | **0.105** | real | 1.46 | 0.300 | 78 | 12 |
| `sealed_ratio` | hi | 0.602 | 0.088 | **0.102** | real | 1.62 | 0.938 | 53 | 8 |
| `sealed_fa` | lo | 0.599 | 0.082 | **0.099** | real | 1.49 | 0.000 | 33 | 5 |
| `dep_rail` | flag | -- | -- | **0.099** | nan | 0.89 | -- | 122 | 31 |
| `sealed_n_days` | hi | 0.598 | 0.082 | **0.098** | real | 1.46 | 0.300 | 78 | 12 |
| `edge` | flag | -- | -- | **0.089** | nan | 0.90 | -- | 128 | 32 |
| `rms` | lo | 0.576 | 0.082 | **0.076** | real | 1.77 | 6.374 | 46 | 6 |
| `co_excess` | hi | 0.432 | 0.083 | **0.068** | nan | 1.00 | -0.130 | 160 | 36 |

Separation curve for `n_det`:

| threshold | real beyond it | null beyond it | % of real | % of null | lift |
|---|---|---|---|---|---|
| 6.000 | 119 | 16 | 74.4 | 44.4 | 1.67 |
| 7.000 | 102 | 11 | 63.7 | 30.6 | 2.09 |
| 8.000 | 86 | 5 | 53.8 | 13.9 | 3.87 |
| 8.850 | 69 | 4 | 43.1 | 11.1 | 3.88 |
| 11.000 | 55 | 2 | 34.4 | 5.6 | 6.19 |
| 15.950 | 36 | 2 | 22.5 | 5.6 | 4.05 |
| 43.000 | 20 | 0 | 12.5 | 0.0 | -- |
| 363.750 | 2 | 0 | 1.2 | 0.0 | -- |

Reverse separations on this bed (the NULL population scores better than the real one): `med_cc` AUC 0.273, `dedupe_n_partners` AUC 0.363, `max_nmad` AUC 0.364.  A criterion on which noise beats signal is evidence about the chain, not a candidate filter.

**2008 Puget / S-Olympic (agent N)** -- the 25 highest-ranking REAL rows on `n_det`, the strongest discriminator for this bed.  `null %ile` is the fraction of the NULL population this row beats on that column; a row near 1.00 is unusual for the null, a row near 0.50 is typical of it.

| grp | n_det | null %ile | n_det | n_days | lat | lon | dep | edge | dep_rail | ctrl_pass | near_chestler_km | sealed_n_det | sealed_n_days | sealed_fa |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 85 | 597 | 1.000 | 597 | 14 | 47.1 | -125.2 | 59.0 | True | True | False | 161.2 | 0 | 0 | 0.133 |
| 103 | 549 | 1.000 | 549 | 14 | 46.8 | -122.4 | 59.0 | True | True | False | 113.6 | 0 | 0 | 0.000 |
| 141 | 354 | 1.000 | 354 | 20 | 47.7 | -123.8 | 59.0 | True | True | False | 48.2 | 2 | 2 | 0.067 |
| 88 | 345 | 1.000 | 345 | 14 | 47.5 | -123.1 | 59.0 | True | True | False | 43.8 | 1 | 1 | 0.067 |
| 47 | 342 | 1.000 | 342 | 20 | 47.4 | -122.8 | 31.0 | False | False | False | 43.6 | 3 | 2 | 0.200 |
| 83 | 243 | 1.000 | 243 | 14 | 47.5 | -122.8 | 59.0 | True | True | False | 39.3 | 1 | 1 | 0.067 |
| 115 | 195 | 1.000 | 195 | 13 | 47.4 | -122.4 | 59.0 | True | True | False | 56.1 | 0 | 0 | 0.000 |
| 36 | 195 | 1.000 | 195 | 6 | 46.5 | -122.4 | 29.0 | True | False | False | 145.9 | 41 | 1 | 18.7 |
| 82 | 113 | 1.000 | 113 | 15 | 47.5 | -122.4 | 59.0 | True | True | False | 53.2 | 0 | 0 | 0.000 |
| 84 | 103 | 1.000 | 103 | 13 | 47.3 | -123.0 | 59.0 | True | True | False | 53.5 | 2 | 1 | 0.267 |
| 81 | 88 | 1.000 | 88 | 12 | 47.9 | -122.1 | 59.0 | True | True | False | 40.3 | 0 | 0 | 0.267 |
| 107 | 85 | 1.000 | 85 | 17 | 47.3 | -122.6 | 23.0 | False | False | True | 65.6 | 0 | 0 | 0.067 |
| 25 | 62 | 1.000 | 62 | 4 | 47.6 | -122.8 | 49.0 | False | False | False | 24.5 | 1 | 1 | 0.333 |
| 80 | 60 | 1.000 | 60 | 2 | 47.3 | -121.8 | 17.0 | True | False | False | 95.4 | 0 | 0 | 0.267 |
| 75 | 54 | 1.000 | 54 | 17 | 47.5 | -123.0 | 41.0 | False | False | False | 28.7 | 0 | 0 | 0.067 |
| 134 | 51 | 1.000 | 51 | 13 | 47.4 | -122.8 | 59.0 | True | True | False | 45.5 | 1 | 1 | 0.133 |
| 159 | 48 | 1.000 | 48 | 16 | 47.4 | -122.8 | 59.0 | True | True | False | 45.5 | 0 | 0 | 0.067 |
| 73 | 45 | 1.000 | 45 | 14 | 47.0 | -122.8 | 59.0 | True | True | False | 94.9 | 1 | 1 | 0.067 |
| 102 | 45 | 1.000 | 45 | 1 | 47.5 | -122.4 | 59.0 | True | True | False | 46.2 | 0 | 0 | 0.067 |
| 20 | 44 | 1.000 | 44 | 7 | 47.4 | -122.6 | 59.0 | True | True | False | 51.6 | 70 | 3 | 0.133 |
| 15 | 42 | 1.000 | 42 | 6 | 47.4 | -122.9 | 59.0 | True | True | False | 50.8 | 25 | 3 | 0.133 |
| 158 | 41 | 1.000 | 41 | 12 | 47.8 | -124.2 | 59.0 | True | True | True | 68.9 | 1 | 1 | 0.000 |
| 76 | 40 | 1.000 | 40 | 17 | 47.3 | -122.8 | 59.0 | True | True | False | 60.6 | 1 | 1 | 0.000 |
| 63 | 38 | 1.000 | 38 | 8 | 47.5 | -122.5 | 59.0 | True | True | False | 49.9 | 1 | 1 | 0.133 |
| 59 | 33 | 1.000 | 33 | 7 | 47.1 | -125.2 | 59.0 | True | True | False | 164.0 | 1 | 1 | 0.333 |

### 2010 Arm B, Chestler's own chain (agent O)

| criterion | direction | AUC real vs null | +- | separation \|AUC-0.5\| | favours | best lift | at threshold | n real | n null |
|---|---|---|---|---|---|---|---|---|---|
| `dedupe_kept` | flag | -- | -- | **0.537** | nan | 0.46 | -- | 342 | 50 |
| `n_days` | hi | 0.024 | 0.022 | **0.476** | nan | 1.00 | 9.000 | 739 | 50 |
| `dedupe_n_partners` | lo | 0.143 | 0.050 | **0.357** | nan | 0.99 | 178.120 | 731 | 50 |
| `ctrl_pass` | flag | -- | -- | **0.350** | real | 9.74 | -- | 288 | 2 |
| `rms` | lo | 0.815 | 0.055 | **0.315** | real | 5.70 | 4.296 | 420 | 5 |
| `sealed_n_days` | hi | 0.193 | 0.056 | **0.307** | nan | 1.00 | 0.000 | 739 | 50 |
| `co_timing` | hi | 0.741 | 0.062 | **0.241** | real | 4.13 | 0.720 | 427 | 7 |
| `depth_window_15_58` | flag | -- | -- | **0.237** | real | 2.32 | -- | 308 | 9 |
| `med_cc` | hi | 0.735 | 0.062 | **0.235** | real | 4.83 | 0.109 | 428 | 6 |
| `ctrl_margin` | hi | 0.665 | 0.067 | **0.165** | real | 3.00 | 0.938 | 354 | 8 |
| `sealed_fa` | lo | 0.660 | 0.067 | **0.160** | real | 2.03 | 0.000 | 390 | 13 |
| `n_det` | hi | 0.343 | 0.067 | **0.157** | nan | 1.29 | 190.720 | 210 | 11 |
| `rate_tremor` | hi | 0.343 | 0.067 | **0.157** | nan | 1.29 | 6.577 | 210 | 11 |
| `sealed_rate` | hi | 0.348 | 0.067 | **0.152** | nan | 1.00 | 0.000 | 739 | 50 |
| `sealed_n_det` | hi | 0.348 | 0.067 | **0.152** | nan | 1.00 | 0.000 | 739 | 50 |
| `fa_quiet` | lo | 0.643 | 0.068 | **0.143** | real | 2.94 | 0.400 | 304 | 7 |
| `edge` | flag | -- | -- | **0.139** | nan | 0.85 | -- | 577 | 46 |
| `tierA` | flag | -- | -- | **0.137** | real | -- | -- | 101 | 0 |

Separation curve for `rms`:

| threshold | real beyond it | null beyond it | % of real | % of null | lift |
|---|---|---|---|---|---|
| 3.605 | 280 | 3 | 38.0 | 6.0 | 6.33 |
| 3.965 | 349 | 5 | 47.4 | 10.0 | 4.74 |
| 4.296 | 420 | 5 | 57.0 | 10.0 | 5.70 |
| 4.803 | 490 | 6 | 66.5 | 12.0 | 5.54 |
| 5.381 | 556 | 10 | 75.4 | 20.0 | 3.77 |
| 5.998 | 617 | 20 | 83.7 | 40.0 | 2.09 |
| 6.856 | 676 | 32 | 91.7 | 64.0 | 1.43 |
| 8.370 | 730 | 49 | 99.1 | 98.0 | 1.01 |

Reverse separations on this bed (the NULL population scores better than the real one): `n_days` AUC 0.024, `dedupe_n_partners` AUC 0.143, `sealed_n_days` AUC 0.193, `n_det` AUC 0.343, `rate_tremor` AUC 0.343, `sealed_rate` AUC 0.348, `sealed_n_det` AUC 0.348.  A criterion on which noise beats signal is evidence about the chain, not a candidate filter.

**2010 Arm B, Chestler's own chain (agent O)** -- the 25 highest-ranking REAL rows on `rms`, the strongest discriminator for this bed.  `null %ile` is the fraction of the NULL population this row beats on that column; a row near 1.00 is unusual for the null, a row near 0.50 is typical of it.

| grp | rms | null %ile | n_det | n_days | lat | lon | dep | edge | dep_rail | ctrl_pass | near_chestler_km | sealed_n_det | sealed_n_days | sealed_fa |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 500 | 0.583 | 1.000 | 3077 | 25 | 47.5 | -122.7 | 15.0 | True | True | True | 44.5 | 15 | 12 | 0.000 |
| 243 | 0.678 | 1.000 | 137 | 18 | 47.9 | -123.2 | 27.0 | False | False | True | 7.818 | 56 | 12 | 0.143 |
| 332 | 0.770 | 1.000 | 429 | 27 | 46.9 | -122.3 | 19.0 | False | False | True | 110.7 | 18 | 9 | 0.000 |
| 499 | 0.790 | 1.000 | 3308 | 23 | 47.3 | -122.4 | 15.0 | True | True | True | 67.6 | 14 | 10 | 0.000 |
| 673 | 0.807 | 1.000 | 1653 | 24 | 47.4 | -122.3 | 15.0 | True | True | True | 69.3 | 25 | 14 | 0.000 |
| 304 | 0.812 | 1.000 | 83 | 23 | 47.9 | -122.8 | 15.0 | True | True | True | 23.3 | 21 | 11 | 0.286 |
| 672 | 0.963 | 1.000 | 1394 | 22 | 47.4 | -122.6 | 15.0 | True | True | True | 52.5 | 22 | 12 | 0.143 |
| 671 | 0.991 | 1.000 | 2410 | 22 | 47.4 | -122.6 | 15.0 | True | True | True | 51.1 | 31 | 18 | 0.286 |
| 669 | 1.021 | 1.000 | 2491 | 24 | 47.4 | -122.6 | 15.0 | True | True | True | 51.1 | 27 | 13 | 0.286 |
| 489 | 1.032 | 1.000 | 404 | 19 | 48.2 | -123.4 | 15.0 | True | True | True | 23.4 | 368 | 18 | 0.143 |
| 206 | 1.038 | 1.000 | 82 | 27 | 47.6 | -123.7 | 49.0 | False | False | True | 45.3 | 10 | 7 | 0.000 |
| 33 | 1.039 | 1.000 | 82 | 24 | 47.4 | -122.6 | 15.0 | True | True | True | 53.8 | 18 | 11 | 0.143 |
| 10 | 1.055 | 1.000 | 79 | 19 | 47.9 | -123.0 | 15.0 | True | True | True | 22.2 | 18 | 10 | 0.000 |
| 430 | 1.060 | 1.000 | 369 | 20 | 48.3 | -123.3 | 15.0 | True | True | True | 24.7 | 274 | 13 | 0.143 |
| 240 | 1.114 | 1.000 | 113 | 18 | 47.9 | -123.1 | 27.0 | False | False | True | 12.1 | 66 | 11 | 0.000 |
| 374 | 1.127 | 1.000 | 69 | 21 | 48.8 | -123.1 | 49.0 | True | False | True | 54.1 | 14 | 7 | 0.000 |
| 402 | 1.127 | 1.000 | 393 | 20 | 48.2 | -123.4 | 15.0 | True | True | True | 23.4 | 288 | 14 | 0.286 |
| 503 | 1.131 | 1.000 | 3910 | 26 | 47.4 | -122.5 | 15.0 | True | True | True | 55.6 | 27 | 18 | 0.571 |
| 690 | 1.138 | 1.000 | 53 | 14 | 48.0 | -121.9 | 43.0 | False | False | True | 55.9 | 13 | 7 | 0.000 |
| 658 | 1.150 | 1.000 | 52 | 15 | 47.8 | -122.7 | 59.0 | True | True | True | 15.1 | 18 | 11 | 0.000 |
| 413 | 1.179 | 1.000 | 307 | 18 | 48.2 | -123.5 | 19.0 | False | False | True | 20.2 | 240 | 11 | 0.143 |
| 736 | 1.189 | 1.000 | 111 | 16 | 47.6 | -122.8 | 31.0 | False | False | True | 23.0 | 13 | 7 | 0.000 |
| 248 | 1.198 | 1.000 | 134 | 19 | 47.8 | -123.2 | 23.0 | False | False | True | 6.027 | 34 | 11 | 0.000 |
| 498 | 1.218 | 1.000 | 1563 | 25 | 47.4 | -122.6 | 15.0 | True | True | True | 51.1 | 20 | 13 | 0.286 |
| 207 | 1.249 | 1.000 | 63 | 24 | 47.6 | -123.9 | 59.0 | True | True | True | 62.6 | 17 | 9 | 0.000 |


## 3. Distance to the published Olympic coordinates, as a column

`derived_coordinate_rows_not_a_catalog.csv` is not a catalogue; it is 69 coordinates derived from
Chestler & Creager (2017) / Sweet et al.  It is used here only as one more column, with the null
answering the same question.

| bed | chain | located rows | median km to nearest published coordinate | within 5 km | within 10 km | within 20 km | of the 69 coordinates, hit within 10 km |
|---|---|---|---|---|---|---|---|
| 2010 N-Olympic, arrays as arrays (agent O) | real | 971 | 50.4 | 2 | 28 | 126 | 11 |
| 2010 N-Olympic, arrays as arrays (agent O) | null | 761 | 62.6 | 2 | 15 | 66 | 6 |
| 2008 Puget / S-Olympic (agent N) | real | 156 | 65.5 | 0 | 0 | 2 | 0 |
| 2008 Puget / S-Olympic (agent N) | null | 36 | 61.6 | 0 | 1 | 2 | 1 |
| 2010 Arm B, Chestler's own chain (agent O) | real | 737 | 37.5 | 3 | 30 | 123 | 14 |
| 2010 Arm B, Chestler's own chain (agent O) | null | 50 | 49.0 | 0 | 0 | 4 | 0 |


## 4. Two measurements that are reported, and are NOT used to cut anything

### 4.1 The confirmation rule calibrated on Bostock's own families

`results/CONFIRMATION_CALIBRATION_NORTH.md`.  The 2005 ETS window was cut into calendar halves;
every template of the 58 rediscovered families and of a matched-size noise sample was re-stacked
from **H1 detections only** and re-fired unchanged on H2 and on the quiet days.  The strictest
values that still pass >= 85 % of the testable true families are **>= 1
distinct tremor day** and **quiet-day false alarms <= 1.0 per
template-day**; the location-move clause **cannot be calibrated at all** -- no finite value passes 85 % of the true families while rejecting a majority of the noise, confirming SCREEN_VALIDATION_NORTH.md sec.3.1 on the confirmation side.

Two things this measurement says, both of which matter more than the thresholds:
* only **18 of the 58 rediscovered families could be
  tested at all** -- the rest do not have enough detections in half an episode to re-stack.  A
  confirmation protocol built on half a window silently discards two thirds of a real catalogue
  before it starts.
* at the calibrated thresholds the rule confirms 16 of the
  18 testable true families
  (89 %) **and
  60 of 84 noise templates
  (71 %)**.  A rule tuned to keep real families
  keeps almost as much noise: on this bed, sealed-window re-fire is a weak discriminator, not a
  verdict.  That is why the sealed-year outcome appears above as columns
  (`sealed_n_days`, `sealed_fa`, `sealed_ratio`, `sealed_move_km`) and not as a gate.

### 4.2 The depth window (15-58 km) -- measured, and rejected as a cut

The window was in an earlier draft of this analysis.  It is not used, for three measured reasons.

| chain | dedupe-kept and located | would be rejected by 15-58 km | of those, on the bottom depth rail (59 km) | genuinely interior depths outside the window | kept while sitting on the 15-km floor rail |
|---|---|---|---|---|---|
| armb_null | 50 | 41 | 41 | 0 | 3 |
| armb_real | 341 | 214 | 214 | 0 | 54 |
| oly2010_null | 642 | 166 | 166 | 0 | 375 |
| oly2010_real | 511 | 139 | 139 | 0 | 277 |
| puget2008_null | 35 | 27 | 27 | 0 | 3 |
| puget2008_real | 125 | 84 | 84 | 0 | 10 |

* **Every single rejection on every Washington chain is the bottom grid rail, and none is a
  genuinely interior solution outside the window.**  The window is a grid-rail clause wearing a
  physical label -- and an asymmetric one, since it deletes the 59-km rail while keeping hundreds
  of solutions on the 15-km floor rail.
* On the north bed the literal window is **vacuous**: that locator's grid is
  18-56 km, so 100 % of the located true
  templates pass it for a reason that has nothing to do with the families.
* Its transferable form -- "reject a solution that rails on a depth face" -- keeps only
  **66 % of the located
  true north families** (33 of
  50), far below the 85 % veto bar.
* The literature bar is physical, not numerical: Bostock and Chestler place families at roughly
  28-45 km, but they get there with a 3-D velocity model and double-difference relocation.  With a
  1-D model and no P picks our depths are coarse -- a 28-45 km window would keep
  24
  of 50 located true north families -- so a tight window is not
  transferable.

For completeness, the one location clause that does pass the veto bar is the locator's own
**>= 5 stations**: 50 of 58 true
families (86 %), reproducing agent Q's
86 %.  It is still reported as the `located` flag, with the unlocatable rows retained.

## 5. Superseded comparison: agent Q's tier A/B/C as columns

Tier A/B/C were a threshold design and are kept only so the numbers can be compared with the
earlier write-ups.  They are columns in every table (`tierA`, `tierB`, `tierC`) and appear in the
separation ranking above like any other column; they are not used to select anything.

## 6. Verdict

**Real and null DO separate, and here is where.**  The largest separation IN FAVOUR OF THE REAL CHAIN anywhere is on **2010 Arm B, Chestler's own chain (agent O)**, column `rms`: AUC 0.815 (0.5 = indistinguishable), i.e. |AUC-0.5| = 0.315, best lift 5.70x (favours real).  Every other column and bed is listed above with its own AUC, so the reader can see exactly which measurements carry the difference and which carry none.

Read the negative precisely.  It is **not** "nothing was found in the raw output": the published
recipe grows thousands of templates on both Washington beds, hundreds survive the literature
dedupe and locate at plausible interface depths, and they are all in the tables with their
measurements.  It is **not** "a screen rejected everything" -- no screen is applied here.  It is
the statement that on these beds, measured column by column, the population produced by the real
data and the population produced by data whose inter-site coherence has been destroyed are not
distinguishable by any statistic this project can compute. A real Washington family, if it is in
these tables, is not separable from noise by these measurements -- which is a statement about the
measurements and the beds, and leaves every row available for a future test rather than deleting
it.

## 7. Files

| what | where |
|---|---|
| per-template tables (every row kept, every criterion a column) | `results/sf_final/table_<bed>_<chain>.csv` |
| separation curves and AUC summaries | `results/sf_final/sep_<bed>_sweep.csv`, `sep_<bed>_summary.csv`, `sep_<bed>.json` |
| ranked tables | `results/sf_final/sep_<bed>_ranked.csv` |
| dedupe detail (what merged, on what) | `results/sf_final/<bed>_<chain>_dedupe.csv` |
| sealed-year re-fires (ALL templates) | `results/sf_final/sealed/<bed>_<chain>/` |
| Arm-B site-shift null chain (new) | `results/sf_final/armb_null/` |
| Arm-B locations (new -- agent O's locate stage had crashed) | `results/sf_final/armb_<chain>_loc.csv` |
| confirmation calibration | `results/CONFIRMATION_CALIBRATION_NORTH.md`, `results/sf_final/calib.json`, `results/sf_final/north/` |
| depth-window measurement | `results/sf_final/depth_check.json` |
| code | `scripts/sf_common.py`, `sf_north_split.py`, `sf_calib.py`, `sf_depthcheck.py`, `sf_accept.py`, `sf_table.py`, `sf_separation.py`, `sf_locate_member.py`, `sf_armb_null.py`, `sf_armb_refire.py`, `sf_report.py` |
| queue / log | `scripts/s_final_queue.sh`, `logs/s_final_queue.log` |
