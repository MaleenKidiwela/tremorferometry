# Washington N-Olympic — CONFIRMATION by re-fire on the SEALED 2011 window

Agent M, 2026-09-21.  The 2010 candidate templates are re-fired **unchanged** on data that no part of this project had read before: the sealed 2011-07-20…09-10 N-Olympic set.  Nothing was tuned, re-cut or re-thresholded — `k_grow.matched_filter` (network mean of the per-station normalised 3-C cc, threshold = median + 8 × MAD of that day's own function, 4 s decluster) reads the template waveforms straight out of the 2010 run's `grow/tpl_it0.npz`.  The 17 site-shift-null candidate templates are re-fired on the same 2011 data by the same code: **their confirmation count is the false-confirmation rate.**

## 0. The seal
```
SEAL OPENED 2026-09-21T07:48:19Z
by: agent M (m_seal_open.py, pid 1252286)
reason: confirmation re-fire of the WA 2010 candidate templates on the sealed 2011 window; nothing tuned -- the templates, the 8xMAD matched filter, the stacking and the locator are the unchanged 2010 code (agent M, PLAN.md sec.15 / results/WA_2010_RESULTS.md addendum).
integrity: 5062/5062 manifest entries verified sha256-clean in 67s (317 of them are symlinks into lfe_extension/traces and were followed)
scope: traces_sealed/ is READ ONLY from here on; it is cached to cache20hz/<NET.STA>/2011/<DOY>.npz with the frozen chain and never modified.
```

## 1. The 2011 bed

* window **2011-07-20 … 2011-09-10**, sub-region N-Olympic, tranche T3.
* **89 three-component stations** = the 2010 probe-admitted, 2010-health-kept set (96) minus 6 that do not exist in 2011 (XG.DR02, XG.LC11, XG.LC13, XG.LC14, XG.TB14, XU.W030) minus 1 that fail the same > 30 % vetoed-station-hour bar on the 2011 census (XG.LC03).
* by network: CN 1, PB 10, UW 2, XG 58, XU 18.
* **3052 station-days cached** with the frozen 20 Hz / 1–8 Hz chain (rational twin for the 50 Hz XG/XU sites); median 85 stations per tremor day.
* PNSN candidate (tremor) hours **429** on **26** days (2011-08-07 … 2011-09-01); the same code on the 2010 window reproduces that run's **454 hours / 29 days** exactly.
* strictly-quiet days (no PNSN tremor within 100 km of the region centre, ≥ 3 days from any such day): **7** (2011-07-22, 2011-07-23, 2011-07-24, 2011-09-07, 2011-09-08, 2011-09-09, 2011-09-10).

> **Documented deviation (bed-forced).**  The brief asks for 20 quiet days of 2011.  The sealed set is only the 53-day window around the 2011 ETS and it contains exactly 7 days that satisfy the frozen quiet rule; no more data may be downloaded.  The false-alarm reference is therefore run over **all 27 non-tremor days** of the window and the pre-stated quiet-day rate is evaluated on the 7 strictly-quiet ones.  Both are reported.  With 7 quiet days the ≤ 0.1 / template-day bar is met only by a template with **zero** quiet-day detections.

## 2. Station survival of the templates

* real: template stations present in 2011, median 58 of 63 (min 25, worst fraction 0.88).
* null: template stations present in 2011, median 47 of 51 (min 40, worst fraction 0.90).

## 3. The re-fire

* **real (50 candidates)**: 858 detections on the 26 tremor days; 49/50 templates fire at all, 40 on ≥ 3 distinct days; 2551 detections on the 27 non-tremor days (1.890 per template-day), 625 on the 7 strictly-quiet days (1.786 per template-day).
* **null (17 candidates)**: 31 detections on the 26 tremor days; 14/17 templates fire at all, 4 on ≥ 3 distinct days; 181 detections on the 27 non-tremor days (0.394 per template-day), 37 on the 7 strictly-quiet days (0.311 per template-day).

## 4. Re-location on the wider grid

Grid **46.8–49.0 N / −125.2…−121.8 / 15–60 km** at the 2010 spacing (0.03°, 0.04°, 2 km) — the 2010 box was 47.50–48.80 / −124.80…−122.20 / 20–58 km, and 46 of the 50 real candidates railed on one of its faces.  Both years are re-stacked from their own detections with `k_grow.build_templates` (unchanged) and located with `k_locate`'s picks, demeaned-residual rms and 50-realisation shuffled-station control.

* real: located from 2010 detections 50/50 (39 still on an edge), from 2011 detections 45/50 (35 on an edge); median 3-D move 79.8 km (epicentral 79.4 km).
* null: located from 2010 detections 17/17 (14 still on an edge), from 2011 detections 4/17 (1 on an edge); median 3-D move 90.2 km (epicentral 85.1 km).

## 5. CONFIRMATION (rule fixed before the seal was opened)

> ≥ 3 distinct 2011 tremor days with detections **and** quiet-day rate ≤ 0.1/template-day **and** a 2011 location that is interior to the wider grid and within 10 km (3-D) of the 2010 location.

### **REAL 0 / 50 confirmed.  NULL 0 / 17 confirmed (false-confirmation rate 0/17 = 0.00).**

With the epicentral (horizontal-only) reading of the 10 km clause: real 0, null 0.

Clause-by-clause (how many pass each on its own):

| clause | real | null |
|---|---|---|
| ≥ 3 tremor days | 40/50 | 4/17 |
| quiet rate ≤ 0.1 | 14/50 | 4/17 |
| interior + ≤ 10 km | 0/50 | 0/17 |

### 5a. The three interior 2010 candidates (the actual evidence of the 2010 run)

grp 1656, 1193 and the marginal 1152 are the interior, ≥5-day, co-timing-excess candidates of `WA_2010_RESULTS.md`'s addendum; 429 is the fourth interior one that failed the day/co-timing filter there.

| grp | 2010 lat* | lon* | dep* | 2010 lat | lon | dep | 2011 lat | lon | dep | move km | 2011 edge | n_det 2011 | days 2011 | quiet det | quiet FA/d | CONFIRMED |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 429 | 47.56 | -123.04 | 34.00 | 47.61 | -123.16 | 31.00 | 47.97 | -123.76 | 15.00 | 61.69 | **yes** | 4 | 3 | 2 | 0.29 | no |
| 1152 | 47.92 | -123.84 | 32.00 | 48.99 | -125.20 | 59.00 | 48.21 | -123.64 | 15.00 | 150.39 | **yes** | 149 | 11 | 49 | 7.00 | no |
| 1193 | 48.16 | -123.52 | 30.00 | 48.12 | -124.72 | 15.00 | 47.52 | -124.64 | 15.00 | 66.98 | **yes** | 56 | 10 | 47 | 6.71 | no |
| 1656 | 48.04 | -123.92 | 30.00 | 48.39 | -123.76 | 15.00 | 48.60 | -124.24 | 15.00 | 42.32 | **yes** | 45 | 7 | 0 | 0.00 | no |

\* `2010 lat/lon/dep` with a star = the ORIGINAL 2010 narrow-grid location from `loc_it0.csv`; the unstarred 2010 columns are the wide-grid re-location from the same 2010 detections, which is the one the 10 km clause uses.

### 5b. The 46 edge-railing 2010 candidates

| grp | 2010 lat* | lon* | dep* | 2010 lat | lon | dep | 2011 lat | lon | dep | move km | 2011 edge | n_det 2011 | days 2011 | quiet det | quiet FA/d | CONFIRMED |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 8 | 48.55 | -123.60 | 20.00 | 48.06 | -125.20 | 59.00 | -- | -- | -- | -- | -- | 1 | 1 | 7 | 1.00 | no |
| 32 | 47.56 | -122.24 | 24.00 | 47.55 | -122.80 | 35.00 | -- | -- | -- | -- | -- | 1 | 1 | 1 | 0.14 | no |
| 75 | 47.50 | -122.68 | 20.00 | 47.73 | -123.16 | 59.00 | 46.86 | -123.48 | 15.00 | 108.85 | **yes** | 16 | 5 | 1 | 0.14 | no |
| 99 | 47.50 | -123.40 | 44.00 | 46.89 | -122.04 | 15.00 | 47.79 | -123.64 | 15.00 | 154.46 | **yes** | 4 | 3 | 0 | 0.00 | no |
| 132 | 47.50 | -122.84 | 20.00 | 47.58 | -122.92 | 15.00 | 47.43 | -122.64 | 15.00 | 26.50 | **yes** | 3 | 2 | 4 | 0.57 | no |
| 179 | 47.50 | -122.64 | 20.00 | 47.40 | -121.80 | 19.00 | 47.73 | -123.28 | 15.00 | 114.92 | **yes** | 16 | 7 | 8 | 1.14 | no |
| 181 | 48.19 | -122.28 | 20.00 | 48.30 | -122.64 | 31.00 | 48.30 | -124.80 | 15.00 | 159.64 | **yes** | 6 | 4 | 2 | 0.29 | no |
| 189 | 47.50 | -122.28 | 28.00 | 47.70 | -122.00 | 15.00 | 47.76 | -122.72 | 29.00 | 55.17 | no | 28 | 3 | 6 | 0.86 | no |
| 223 | 47.50 | -122.80 | 20.00 | 47.43 | -122.76 | 29.00 | 47.64 | -122.60 | 17.00 | 28.77 | no | 20 | 4 | 0 | 0.00 | no |
| 271 | 47.50 | -123.04 | 20.00 | 47.31 | -122.96 | 59.00 | 47.85 | -122.88 | 15.00 | 74.67 | **yes** | 3 | 3 | 1 | 0.14 | no |
| 305 | 47.80 | -122.32 | 20.00 | 47.76 | -123.24 | 15.00 | -- | -- | -- | -- | -- | 2 | 2 | 0 | 0.00 | no |
| 319 | 47.83 | -122.24 | 20.00 | 48.12 | -122.12 | 15.00 | 48.15 | -121.80 | 15.00 | 23.77 | **yes** | 3 | 2 | 2 | 0.29 | no |
| 328 | 48.28 | -122.68 | 20.00 | 48.99 | -122.24 | 59.00 | 47.64 | -123.60 | 27.00 | 183.19 | no | 10 | 6 | 1 | 0.14 | no |
| 355 | 47.50 | -122.76 | 20.00 | 47.76 | -122.72 | 59.00 | 47.64 | -123.16 | 15.00 | 56.22 | **yes** | 15 | 6 | 7 | 1.00 | no |
| 369 | 47.59 | -122.80 | 20.00 | 47.19 | -122.80 | 15.00 | 47.88 | -123.36 | 17.00 | 87.10 | no | 4 | 4 | 4 | 0.57 | no |
| 410 | 47.50 | -122.76 | 20.00 | 47.61 | -123.28 | 15.00 | 47.43 | -123.28 | 59.00 | 48.34 | **yes** | 14 | 6 | 44 | 6.29 | no |
| 438 | 47.53 | -122.76 | 20.00 | 47.58 | -123.12 | 57.00 | 47.97 | -123.00 | 15.00 | 61.01 | **yes** | 29 | 4 | 2 | 0.29 | no |
| 451 | 48.04 | -123.12 | 20.00 | 48.00 | -124.00 | 33.00 | 48.00 | -123.28 | 59.00 | 58.98 | **yes** | 3 | 3 | 8 | 1.14 | no |
| 468 | 47.77 | -122.24 | 26.00 | 48.00 | -123.20 | 39.00 | 47.91 | -123.00 | 25.00 | 22.64 | no | 8 | 3 | 1 | 0.14 | no |
| 486 | 47.86 | -124.16 | 20.00 | 47.40 | -122.12 | 29.00 | 47.76 | -123.00 | 17.00 | 77.03 | no | 4 | 3 | 1 | 0.14 | no |
| 497 | 47.50 | -122.68 | 20.00 | 46.89 | -123.76 | 15.00 | 47.16 | -124.40 | 59.00 | 71.08 | **yes** | 5 | 4 | 0 | 0.00 | no |
| 505 | 47.50 | -124.24 | 26.00 | 47.79 | -123.76 | 59.00 | 47.79 | -123.16 | 15.00 | 62.31 | **yes** | 4 | 4 | 0 | 0.00 | no |
| 519 | 47.53 | -122.84 | 20.00 | 47.76 | -123.04 | 59.00 | 47.67 | -123.56 | 59.00 | 39.53 | **yes** | 123 | 6 | 21 | 3.00 | no |
| 521 | 47.50 | -123.20 | 30.00 | 47.28 | -122.68 | 23.00 | 47.88 | -122.84 | 15.00 | 68.22 | **yes** | 56 | 12 | 273 | 39.00 | no |
| 537 | 48.52 | -122.24 | 20.00 | 47.64 | -122.96 | 59.00 | 47.73 | -121.80 | 59.00 | 85.89 | **yes** | 25 | 10 | 43 | 6.14 | no |
| 569 | 48.22 | -122.96 | 20.00 | 47.31 | -125.20 | 25.00 | 48.27 | -123.12 | 53.00 | 188.61 | no | 29 | 9 | 2 | 0.29 | no |
| 629 | 47.50 | -124.08 | 28.00 | 47.58 | -124.64 | 15.00 | 47.67 | -121.84 | 23.00 | 206.30 | no | 3 | 2 | 0 | 0.00 | no |
| 844 | 48.28 | -124.80 | 32.00 | 48.48 | -123.28 | 15.00 | 47.13 | -124.32 | 15.00 | 168.47 | **yes** | 6 | 5 | 6 | 0.86 | no |
| 875 | 48.52 | -122.56 | 20.00 | 48.60 | -122.44 | 15.00 | 48.15 | -121.84 | 27.00 | 67.78 | no | 17 | 9 | 31 | 4.43 | no |
| 901 | 48.16 | -124.60 | 20.00 | 47.52 | -124.08 | 15.00 | 47.97 | -122.72 | 15.00 | 111.83 | **yes** | 6 | 4 | 19 | 2.71 | no |
| 984 | 48.01 | -124.16 | 20.00 | 46.92 | -125.20 | 59.00 | 47.52 | -123.84 | 59.00 | 120.22 | **yes** | 4 | 4 | 0 | 0.00 | no |
| 1029 | 48.13 | -123.12 | 20.00 | 47.85 | -123.64 | 15.00 | 47.16 | -123.92 | 23.00 | 79.84 | no | 5 | 4 | 0 | 0.00 | no |
| 1035 | 47.50 | -123.00 | 20.00 | 46.80 | -122.24 | 59.00 | 47.64 | -122.72 | 15.00 | 109.12 | **yes** | 13 | 9 | 0 | 0.00 | no |
| 1040 | 48.79 | -123.20 | 20.00 | 48.99 | -123.88 | 23.00 | 47.58 | -125.20 | 59.00 | 187.88 | **yes** | 5 | 5 | 0 | 0.00 | no |
| 1063 | 47.92 | -122.56 | 20.00 | 47.70 | -123.08 | 15.00 | 47.67 | -123.24 | 59.00 | 45.67 | **yes** | 3 | 2 | 3 | 0.43 | no |
| 1094 | 47.50 | -124.00 | 20.00 | 47.70 | -123.72 | 59.00 | -- | -- | -- | -- | -- | 0 | 0 | 0 | 0.00 | no |
| 1227 | 47.80 | -122.24 | 20.00 | 47.82 | -121.80 | 15.00 | 47.64 | -123.56 | 15.00 | 130.96 | **yes** | 4 | 4 | 3 | 0.43 | no |
| 1344 | 47.50 | -122.92 | 20.00 | 47.46 | -122.88 | 59.00 | 47.88 | -124.76 | 15.00 | 152.41 | **yes** | 13 | 4 | 3 | 0.43 | no |
| 1399 | 47.50 | -122.56 | 20.00 | 47.76 | -124.04 | 15.00 | 48.63 | -122.60 | 15.00 | 143.43 | **yes** | 7 | 7 | 2 | 0.29 | no |
| 1471 | 48.79 | -122.92 | 58.00 | 47.82 | -121.88 | 29.00 | 47.34 | -124.88 | 15.00 | 227.40 | **yes** | 5 | 5 | 0 | 0.00 | no |
| 1476 | 48.79 | -123.52 | 20.00 | 48.06 | -123.12 | 59.00 | 47.46 | -124.76 | 15.00 | 144.67 | **yes** | 8 | 6 | 5 | 0.71 | no |
| 1506 | 47.77 | -122.32 | 20.00 | 47.76 | -122.40 | 31.00 | 47.73 | -122.92 | 59.00 | 47.51 | **yes** | 7 | 6 | 8 | 1.14 | no |
| 1598 | 47.50 | -122.88 | 20.00 | 47.25 | -121.80 | 19.00 | 47.10 | -123.04 | 59.00 | 100.96 | **yes** | 14 | 6 | 5 | 0.71 | no |
| 1613 | 47.50 | -122.92 | 20.00 | 47.40 | -124.72 | 15.00 | -- | -- | -- | -- | -- | 1 | 1 | 0 | 0.00 | no |
| 1678 | 48.58 | -124.72 | 20.00 | 48.45 | -124.24 | 15.00 | 48.51 | -124.44 | 15.00 | 16.15 | **yes** | 48 | 6 | 2 | 0.29 | no |
| 1685 | 48.79 | -123.48 | 20.00 | 48.96 | -123.48 | 15.00 | 47.88 | -125.20 | 55.00 | 178.94 | **yes** | 3 | 2 | 1 | 0.14 | no |

### 5c. The 17 site-shift-null candidates (false-confirmation control)

| grp | 2010 lat* | lon* | dep* | 2010 lat | lon | dep | 2011 lat | lon | dep | move km | 2011 edge | n_det 2011 | days 2011 | quiet det | quiet FA/d | CONFIRMED |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 486 | 47.50 | -122.92 | 20.00 | 47.37 | -123.44 | 31.00 | -- | -- | -- | -- | -- | 2 | 2 | 0 | 0.00 | no |
| 491 | 48.55 | -122.44 | 22.00 | 47.55 | -123.16 | 15.00 | -- | -- | -- | -- | -- | 0 | 0 | 3 | 0.43 | no |
| 757 | 48.22 | -123.32 | 20.00 | 47.64 | -124.68 | 15.00 | -- | -- | -- | -- | -- | 2 | 2 | 6 | 0.86 | no |
| 802 | 47.74 | -122.80 | 20.00 | 48.06 | -122.88 | 59.00 | -- | -- | -- | -- | -- | 2 | 2 | 2 | 0.29 | no |
| 1077 | 47.53 | -122.84 | 20.00 | 47.55 | -123.12 | 59.00 | 47.91 | -123.32 | 15.00 | 61.28 | **yes** | 6 | 4 | 9 | 1.29 | no |
| 1084 | 47.89 | -122.80 | 32.00 | 47.58 | -122.80 | 59.00 | -- | -- | -- | -- | -- | 1 | 1 | 1 | 0.14 | no |
| 1592 | 47.89 | -122.56 | 20.00 | 47.73 | -124.20 | 59.00 | 47.55 | -123.00 | 25.00 | 96.66 | no | 3 | 3 | 2 | 0.29 | no |
| 1668 | 48.22 | -122.72 | 20.00 | 47.55 | -121.80 | 59.00 | -- | -- | -- | -- | -- | 0 | 0 | 0 | 0.00 | no |
| 1886 | 47.74 | -122.24 | 20.00 | 47.97 | -121.84 | 23.00 | -- | -- | -- | -- | -- | 2 | 2 | 1 | 0.14 | no |
| 1931 | 47.50 | -123.04 | 20.00 | 47.64 | -122.92 | 15.00 | -- | -- | -- | -- | -- | 1 | 1 | 2 | 0.29 | no |
| 1970 | 47.65 | -122.72 | 32.00 | 47.82 | -123.72 | 15.00 | -- | -- | -- | -- | -- | 2 | 2 | 2 | 0.29 | no |
| 1981 | 47.50 | -122.64 | 20.00 | 48.00 | -121.80 | 59.00 | -- | -- | -- | -- | -- | 1 | 1 | 0 | 0.00 | no |
| 2080 | 48.01 | -122.88 | 20.00 | 47.79 | -124.36 | 27.00 | -- | -- | -- | -- | -- | 1 | 1 | 6 | 0.86 | no |
| 2491 | 48.52 | -124.20 | 20.00 | 47.76 | -123.64 | 15.00 | 47.40 | -121.84 | 23.00 | 138.51 | no | 3 | 3 | 1 | 0.14 | no |
| 2553 | 48.16 | -122.24 | 26.00 | 47.58 | -123.84 | 59.00 | -- | -- | -- | -- | -- | 0 | 0 | 1 | 0.14 | no |
| 2739 | 47.98 | -123.44 | 40.00 | 47.88 | -123.08 | 59.00 | -- | -- | -- | -- | -- | 1 | 1 | 0 | 0.00 | no |
| 2977 | 47.92 | -122.60 | 20.00 | 47.79 | -123.72 | 15.00 | 48.30 | -122.96 | 41.00 | 83.76 | no | 4 | 4 | 1 | 0.14 | no |

## 6. Descriptive — distance to the Chestler-derived coordinates

`../codex-LFE/results/olympic_reference_discovery_v1/derived_coordinate_rows_not_a_catalog.csv` is not a catalogue and is used only to say where the confirmed families sit relative to previously published Olympic LFE coordinates.

* no confirmed family to compare.

## 7. Files

| | |
|---|---|
| bed / stations / hours / quiet days | `results/m_2011/bed.json, stations.csv, candidate_hours.csv, quiet_days.csv` |
| 2011 re-fire detections (tremor days) | `results/m_2011/det_real_tremor.parquet, det_null_tremor.parquet` |
| 2011 false-alarm detections (non-tremor days) | `results/m_2011/det_real_fa.parquet, det_null_fa.parquet` |
| template station survival | `results/m_2011/tpl_stations_*.csv` |
| wide-grid re-locations | `results/m_2011/loc2010_*.csv, loc2011_*.csv (+ _tpl.npz re-stacks)` |
| confirmation tables | `results/m_2011/confirm_real.csv, confirm_null.csv` |
| seal | `results/sealed_2011_manifest.sha256, SEALED_2011_OPENED_*` |

## 8. What this does and does not settle

* The re-fire is a genuine out-of-sample test: the 2011 data were hashed and sealed before the 2010 chain was written, and the templates crossing into 2011 were not re-fitted in any way.
* The null templates are the controlled comparison for every number above — they were built by the identical chain on a site-shifted network and are fired by the identical code on the identical 2011 data.
* There is still no Washington LFE catalogue, so nothing here is a rediscovery count and no family is claimed to be a catalogued LFE family.
* The quiet-day clause rests on 7 days (see the deviation in §1); the non-tremor false-alarm rate over the whole window is the more precise, weaker-definition version.

