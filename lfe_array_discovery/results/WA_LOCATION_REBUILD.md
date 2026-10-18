# Rebuilding the LFE family LOCATION step, calibrated on Bostock's published coordinates

Agent U, 2026-09-23.  Workspace `results/u_locate/`; queue `scripts/u_locate_queue.sh` (markers `results/U_LOCATE_DONE|FAILED`, interim `results/U_LOCATE_BASELINE_DONE`); code `scripts/uu_common.py`, `uu_stacks.py`, `uu_picks.py`, `uu_locate.py`, `uu_dd.py`, `uu_eval.py`, `uu_wa.py`.  **Detection and growth are untouched** -- every template, every detection and every stack member comes from the frozen chain on disk; this study only re-reads them and replaces what happens after.

No family is ever filtered out.  Every grown family keeps its row in every table and every criterion is a column (memory `rank-dont-discard`).


## 0. The calibration set, and the proof that it is the published answer

`lfe_catalog_repair/results/e5r_families_2012.csv` holds **130 family coordinates**, and **130 of them** are exactly among the 130 distinct (lat, lon, dep) triples carried by the event rows of Bostock et al. (2015)'s own catalogue export (`results/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv`).  That file is the published location of each family; `bkey` is this project's index into it.  Depth convention: km positive down, 28.8-46.6 km (median 35.3).

The 58 families the frozen recipe rediscovered blind on the 2005 POLARIS bed (`results/repro/classify15.csv`, column `kept_bkey`) are the DEVELOPMENT set.  The 18 located SEALED families (`hid_bkey`) are scored ONCE at the end and no choice is made on them.  One representative template per family is fixed before any variant runs (highest one-to-one match fraction), so every variant is scored on the same 58 objects.

The noise population is the 1,219 grown templates that match no catalogued family; it gets the identical table.


## 1. BASELINE -- what the current locator does, measured for the first time

| variant | n_families | n_located | epi_med | epi_p90 | dep_med | dep_p90 | dep_bias | d3_med | rail | within5 | within10 | within20 | null_d3_med |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BASELINE_on_disk(k_locate, tpl_it0) | 58 | 44 | 27.63 | 87.85 | 15.86 | 22.88 | -11.68 | 35.24 | 0.75 | 4 | 11 | 17 | 39.75 |
| BASELINE_on_disk[BLIND sealed 18] | 22 | 19 | 22.8 | 64.67 | 12.93 | 22.27 | -11.51 | 26.67 | 0.63 | 0 | 3 | 7 | 48.71 |


`epi_*` = epicentral distance from our solution to Bostock's published coordinate for THAT family, in km; `dep_*` = |depth error|; `dep_bias` = median signed depth error; `rail` = fraction of located solutions sitting on a face of the search grid; `null_d3_med` = the same solutions scored against a RANDOMLY PERMUTED family label (the honest chance level).


## 2. Where the error comes from -- per-station pick residuals against the published locations

Residual = observed pick minus the time predicted by Bostock's published coordinate for that family in the Crosson 1-D model, median-centred per family (the origin time is free).

| h_envH_min | pick | n | n_fam | med_abs_res | frac_lt_0p5 | frac_lt_1 | frac_lt_2 | frac_gt_4 |
|---|---|---|---|---|---|---|---|---|
| 0 | t_env | 3262 | 176 | 2.91 | 0.23 | 0.32 | 0.44 | 0.44 |
| 0 | t_envH | 3262 | 176 | 2.86 | 0.24 | 0.33 | 0.44 | 0.44 |
| 0 | t_onS | 3262 | 176 | 3.11 | 0.2 | 0.29 | 0.4 | 0.43 |
| 3 | t_env | 2872 | 176 | 2.3 | 0.26 | 0.37 | 0.48 | 0.4 |
| 3 | t_envH | 2872 | 176 | 2.29 | 0.26 | 0.37 | 0.48 | 0.4 |
| 3 | t_onS | 2872 | 176 | 2.68 | 0.21 | 0.31 | 0.43 | 0.39 |
| 4 | t_env | 1866 | 176 | 0.93 | 0.38 | 0.51 | 0.62 | 0.28 |
| 4 | t_envH | 1866 | 176 | 0.93 | 0.39 | 0.52 | 0.62 | 0.27 |
| 4 | t_onS | 1866 | 176 | 1.53 | 0.3 | 0.41 | 0.56 | 0.26 |
| 5 | t_env | 1161 | 174 | 0.39 | 0.56 | 0.69 | 0.78 | 0.15 |
| 5 | t_envH | 1161 | 174 | 0.36 | 0.57 | 0.7 | 0.8 | 0.14 |
| 5 | t_onS | 1161 | 174 | 0.73 | 0.44 | 0.56 | 0.7 | 0.13 |
| 6 | t_env | 799 | 151 | 0.26 | 0.68 | 0.79 | 0.88 | 0.07 |
| 6 | t_envH | 799 | 151 | 0.25 | 0.68 | 0.81 | 0.88 | 0.06 |
| 6 | t_onS | 799 | 151 | 0.43 | 0.53 | 0.66 | 0.79 | 0.06 |
| 8 | t_env | 467 | 110 | 0.15 | 0.79 | 0.9 | 0.96 | 0.01 |
| 8 | t_envH | 467 | 110 | 0.15 | 0.8 | 0.9 | 0.96 | 0.01 |
| 8 | t_onS | 467 | 110 | 0.21 | 0.7 | 0.82 | 0.91 | 0.01 |
| 10 | t_env | 289 | 72 | 0.14 | 0.87 | 0.96 | 0.98 | 0.0 |
| 10 | t_envH | 289 | 72 | 0.13 | 0.88 | 0.96 | 0.98 | 0.0 |
| 10 | t_onS | 289 | 72 | 0.17 | 0.8 | 0.91 | 0.98 | 0.0 |



| sta | median | std | count |
|---|---|---|---|
| SSIB | -0.24 | 1.45 | 64 |
| SILB | -0.18 | 2.11 | 53 |
| MGCB | -0.06 | 1.64 | 71 |
| TWKB | -0.05 | 3.17 | 84 |
| GLBC | -0.04 | 1.36 | 53 |
| PGC | -0.01 | 2.46 | 64 |
| VGZ | 0.0 | 2.81 | 40 |
| SHVB | 0.0 |  | 1 |
| LCBC | 0.02 | 1.96 | 46 |
| LZB | 0.03 | 2.6 | 68 |
| JRBC | 0.03 | 1.03 | 48 |
| TSJB | 0.13 | 2.1 | 60 |
| TWBB | 0.22 | 1.82 | 61 |
| SNB | 0.58 | 0.79 | 24 |
| B005 | 1.26 | 5.75 | 10 |
| B006 | 1.29 | 1.65 | 12 |
| B007 | 1.38 | 7.71 | 11 |
| B001 | 1.58 | 8.51 | 11 |
| B003 | 3.6 | 6.97 | 18 |


**Reading.** the travel times are NOT the problem: at stack envelope peaks of >= 8 MAD the observed S time agrees with the Crosson 1-D prediction from Bostock's published coordinate to a median of 0.15 s and 90 % of picks land within 1 s.  The PICK is: over all admitted stations the median is 2.86 s and 44 % of picks are more than 4 s out -- the envelope maximum of a low-SNR stack lands on the wrong bump, and the unweighted least-squares fit of the current locator lets those outliers drag the solution onto a grid face.  Station medians span -0.24 to +3.60 s, i.e. the model/station-term error is small beside the pick scatter.  The AIC ONSET, which Chestler uses, is WORSE than the envelope maximum at every SNR level here (the LFE stack onset is emergent), so it is discarded; the P onset and the S-P envelope cross-correlation are also measured and are too poor to help (below).


## 3. The ladder -- every change, measured on the same 58 families

| variant | n_located | epi_med | epi_p90 | dep_med | dep_bias | d3_med | rail | within5 | within10 | within20 | noise_located | noise_rail |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| L0_base_env_L2 | 56 | 25.93 | 72.95 | 16.48 | -10.37 | 31.1 | 0.7 | 3 | 15 | 22 | 1015.0 | 0.92 |
| L10_xc_refine | 46 | 4.7 | 19.32 | 5.31 | 1.29 | 7.41 | 0.09 | 24 | 35 | 42 | 452.0 | 0.47 |
| L11_xc_zprior | 46 | 4.94 | 16.4 | 5.74 | -1.63 | 10.73 | 0.0 | 23 | 36 | 43 | 452.0 | 0.05 |
| L12_comb | 56 | 6.11 | 62.28 | 8.59 | -1.12 | 11.83 | 0.16 | 23 | 39 | 45 | 995.0 | 0.59 |
| L13_comb_zprior | 56 | 5.1 | 58.44 | 7.61 | -4.49 | 11.74 | 0.02 | 28 | 40 | 46 | 995.0 | 0.08 |
| L1_wide | 56 | 26.26 | 107.6 | 27.81 | -13.65 | 39.99 | 0.55 | 2 | 14 | 23 | 1015.0 | 0.83 |
| L2_robust | 56 | 10.49 | 73.21 | 17.48 | -7.5 | 22.3 | 0.29 | 9 | 27 | 35 | 1015.0 | 0.72 |
| L3_snrw | 56 | 9.08 | 97.83 | 13.69 | -8.83 | 20.48 | 0.27 | 17 | 29 | 39 | 1015.0 | 0.68 |
| L4_snrw2 | 56 | 12.46 | 96.0 | 13.38 | -6.49 | 20.18 | 0.29 | 17 | 24 | 36 | 1015.0 | 0.67 |
| L5_horiz | 56 | 8.99 | 89.18 | 16.56 | -9.12 | 22.49 | 0.25 | 15 | 31 | 41 | 995.0 | 0.66 |
| L6_terms | 56 | 8.81 | 90.08 | 16.03 | -6.03 | 20.1 | 0.27 | 23 | 32 | 41 | 995.0 | 0.67 |
| L7_refine | 56 | 8.45 | 90.55 | 15.29 | -5.08 | 19.57 | 0.23 | 20 | 31 | 40 | 995.0 | 0.67 |
| L8_bp | 56 | 15.63 | 72.17 | 16.25 | -8.48 | 28.17 | 0.21 | 16 | 22 | 30 | 995.0 | 0.18 |
| L8b_irls | 56 | 7.86 | 87.85 | 10.6 | -4.06 | 17.21 | 0.21 | 22 | 29 | 41 | 995.0 | 0.62 |
| L8c_irls_adapt | 56 | 9.26 | 87.85 | 11.25 | -3.0 | 17.22 | 0.23 | 20 | 28 | 38 | 995.0 | 0.6 |
| L8d_irls_zprior | 56 | 8.25 | 67.56 | 9.5 | -7.44 | 15.46 | 0.04 | 18 | 32 | 40 | 995.0 | 0.07 |
| L9_xcdiff | 46 | 4.0 | 19.39 | 5.24 | 2.07 | 6.44 | 0.09 | 26 | 35 | 41 | 452.0 | 0.47 |


## 4. Ablations -- what was tried and DISCARDED because it did not improve the error

| variant | n_located | epi_med | epi_p90 | dep_med | dep_bias | d3_med | rail | within10 |
|---|---|---|---|---|---|---|---|---|
| A_L2_at_L7 | 56 | 13.3 | 75.13 | 26.41 | -7.84 | 33.12 | 0.46 | 24 |
| A_adaptive | 56 | 8.85 | 86.01 | 13.11 | -5.42 | 17.7 | 0.27 | 30 |
| A_comb_xc2 | 56 | 5.72 | 49.96 | 7.4 | -3.28 | 11.56 | 0.0 | 39 |
| A_minpeak3 | 57 | 8.73 | 79.71 | 12.67 | -6.97 | 17.16 | 0.18 | 30 |
| A_minpeak5 | 34 | 4.14 | 21.84 | 6.38 | -2.49 | 8.76 | 0.09 | 26 |
| A_no_terms | 56 | 8.86 | 90.91 | 16.03 | -9.97 | 21.17 | 0.21 | 30 |
| A_onset | 56 | 16.94 | 71.72 | 24.44 | -11.87 | 28.51 | 0.32 | 16 |
| A_slab2_depth | 56 | 7.7 | 67.6 | 9.6 | -7.02 | 14.55 | 0.02 | 32 |
| A_sp_onP | 56 | 19.68 | 74.38 | 27.36 | 23.29 | 38.51 | 0.21 | 13 |
| A_sp_xc | 56 | 12.52 | 81.2 | 18.4 | -6.82 | 24.25 | 0.18 | 23 |
| A_vpvs_172 | 56 | 8.55 | 94.24 | 18.52 | -7.52 | 23.28 | 0.25 | 31 |
| A_vpvs_185 | 56 | 8.4 | 86.81 | 13.45 | -4.05 | 17.09 | 0.23 | 31 |
| A_wpow3 | 56 | 8.07 | 81.3 | 12.82 | -6.32 | 17.4 | 0.25 | 30 |
| A_wpow4 | 56 | 9.26 | 91.6 | 14.59 | -6.32 | 17.88 | 0.2 | 28 |
| A_xc_cc50 | 29 | 2.87 | 8.11 | 5.25 | -3.11 | 6.88 | 0.0 | 27 |
| A_xc_noterms | 46 | 5.13 | 15.82 | 5.8 | -2.01 | 10.63 | 0.0 | 37 |
| A_zprior8 | 56 | 9.22 | 62.96 | 10.49 | -6.41 | 14.84 | 0.04 | 30 |


## 5. The adopted locator

**Adopted: `L9_xcdiff`** -- {'obs': 'xc', 'hcol': 'h_envH', 'robust': True, 'wpow': 2, 'grid': 'wide', 'terms': 'loo'}

| variant | n_located | epi_med | epi_p90 | dep_med | d3_med | rail | within5 | within10 | within20 | null_d3_med |
|---|---|---|---|---|---|---|---|---|---|---|
| BASELINE_on_disk(k_locate, tpl_it0) | 44 | 27.63 | 87.85 | 15.86 | 35.24 | 0.75 | 4 | 11 | 17 | 39.75 |
| L9_xcdiff | 46 | 4.0 | 19.39 | 5.24 | 6.44 | 0.09 | 26 | 35 | 41 | 34.72 |


Acceptance (stated before the run): adopt only if the median distance to Bostock's published coordinates improves materially AND the rail fraction drops, on the 58 known families.  Measured: epicentral median 27.6 -> 4.0 km, 3-D median 35.2 -> 6.4 km, rail 0.75 -> 0.09, families located 44.0 -> 46.0 of 58.  Verdict: **ADOPTED**.


## 6. Blind confirmation -- the 18 located SEALED families

| variant | n_families | n_located | epi_med | epi_p90 | dep_med | d3_med | within10 |
|---|---|---|---|---|---|---|---|
| BASELINE_on_disk[BLIND sealed 18] | 22 | 19 | 22.8 | 64.67 | 12.93 | 26.67 | 3 |


## 7. Per-family uncertainty

Station bootstrap (resample the picking stations with replacement, 60 realisations, 68th percentile of the epicentral scatter): median 8.0 km, p90 31.8 km.  Correlation with the TRUE error: 0.81 -- this is the number the Washington catalogue carries per family.

| bkey | grp | n_sta | lat | lon | dep | boot_km | boot_dz | err_epi | err_dep |
|---|---|---|---|---|---|---|---|---|---|
| 97 | 1294 | 15.0 | 48.62 | -124.08 | 36.0 | 3.34 | 2.0 | 1.33 | -2.14 |
| 81 | 1105 | 9.0 | 48.53 | -124.04 | 34.0 | 6.77 | 12.0 | 1.33 | 1.46 |
| 56 | 954 | 16.0 | 48.47 | -123.8 | 36.0 | 0.0 | 0.0 | 1.33 | 2.19 |
| 86 | 985 | 13.0 | 48.56 | -123.76 | 40.0 | 13.34 | 36.0 | 1.47 | 3.45 |
| 37 | 744 | 12.0 | 48.41 | -123.72 | 30.0 | 2.95 | 4.0 | 1.85 | -3.23 |
| 76 | 1001 | 13.0 | 48.53 | -123.84 | 32.0 | 2.95 | 4.0 | 2.21 | -2.39 |
| 104 | 1058 | 13.0 | 48.65 | -123.72 | 42.0 | 10.06 | 10.0 | 2.34 | 1.99 |
| 96 | 1036 | 14.0 | 48.62 | -123.72 | 42.0 | 3.34 | 6.0 | 2.34 | 2.16 |
| 47 | 1169 | 17.0 | 48.44 | -123.68 | 34.0 | 3.34 | 2.0 | 2.48 | -0.32 |
| 67 | 920 | 12.0 | 48.5 | -123.64 | 42.0 | 2.95 | 4.24 | 2.67 | 5.04 |
| 53 | 836 | 7.0 | 48.44 | -123.96 | 38.0 | 13.55 | 10.0 | 2.67 | 5.85 |
| 5 | 574 | 12.0 | 48.23 | -123.68 | 36.0 | 8.92 | 10.0 | 2.67 | 4.17 |
| 105 | 1069 | 14.0 | 48.65 | -123.68 | 36.0 | 0.0 | 2.0 | 3.13 | -4.21 |
| 54 | 810 | 16.0 | 48.44 | -123.92 | 38.0 | 3.0 | 4.0 | 3.14 | 5.58 |
| 6 | 502 | 16.0 | 48.23 | -123.6 | 32.0 | 3.34 | 6.24 | 3.14 | 1.36 |
| 92 | 1132 | 8.0 | 48.59 | -123.48 | 36.0 | 3.34 | 20.0 | 3.15 | -3.25 |
| 49 | 737 | 13.0 | 48.41 | -123.72 | 42.0 | 20.45 | 34.0 | 3.34 | 7.68 |
| 111 | 1250 | 18.0 | 48.8 | -124.04 | 36.0 | 4.44 | 2.0 | 3.42 | -5.42 |
| 73 | 1003 | 10.0 | 48.5 | -123.92 | 28.0 | 4.45 | 18.0 | 3.42 | -5.46 |
| 72 | 1028 | 12.0 | 48.5 | -123.96 | 38.0 | 2.95 | 6.0 | 3.65 | 4.61 |
| 39 | 1114 | 12.0 | 48.38 | -123.96 | 38.0 | 14.6 | 32.0 | 3.65 | 7.07 |
| 66 | 987 | 10.0 | 48.53 | -123.96 | 38.0 | 4.73 | 10.0 | 3.69 | 5.06 |
| 101 | 1017 | 17.0 | 48.59 | -123.52 | 36.0 | 0.0 | 2.0 | 4.0 | -4.66 |
| 41 | 151 | 12.0 | 48.44 | -123.16 | 32.0 | 8.91 | 10.24 | 4.0 | -10.37 |
| 68 | 804 | 14.0 | 48.56 | -123.56 | 38.0 | 13.67 | 10.0 | 4.69 | -1.42 |
| 21 | 481 | 13.0 | 48.35 | -123.48 | 26.0 | 12.82 | 16.0 | 4.69 | -8.35 |
| 3 | 404 | 17.0 | 48.23 | -123.56 | 30.0 | 3.34 | 4.0 | 5.35 | -0.08 |
| 116 | 383 | 8.0 | 48.77 | -123.64 | 32.0 | 8.88 | 14.0 | 5.75 | -11.12 |
| 16 | 191 | 15.0 | 48.23 | -123.08 | 40.0 | 11.63 | 10.24 | 5.75 | -1.24 |
| 0 | 705 | 11.0 | 48.17 | -123.68 | 36.0 | 4.46 | 11.36 | 5.75 | 6.91 |
| 33 | 173 | 11.0 | 48.41 | -123.12 | 34.0 | 2.95 | 2.24 | 6.15 | -9.31 |
| 91 | 1160 | 12.0 | 48.62 | -123.84 | 40.0 | 13.53 | 14.56 | 6.8 | 3.27 |
| 9 | 503 | 14.0 | 48.26 | -123.52 | 32.0 | 7.05 | 14.0 | 7.41 | 1.0 |
| 42 | 834 | 15.0 | 48.35 | -123.88 | 74.0 | 13.8 | 14.24 | 8.09 | 41.34 |
| 107 | 1164 | 17.0 | 48.65 | -123.68 | 32.0 | 3.34 | 4.0 | 8.42 | -7.85 |
| 55 | 833 | 7.0 | 48.41 | -123.68 | 12.0 | 0.0 | 2.0 | 11.09 | -21.53 |
| 127 | 1211 | 9.0 | 48.8 | -123.8 | 48.0 | 30.9 | 28.0 | 11.59 | 7.46 |
| 44 | 658 | 8.0 | 48.29 | -123.64 | 62.0 | 20.57 | 14.0 | 14.92 | 27.03 |
| 46 | 837 | 11.0 | 48.47 | -123.72 | 76.0 | 42.04 | 3.44 | 15.41 | 43.82 |
| 23 | 617 | 12.0 | 48.17 | -123.72 | 42.0 | 12.25 | 16.88 | 18.53 | 10.41 |
| 117 | 1212 | 8.0 | 48.74 | -123.8 | 38.0 | 72.4 | 38.0 | 18.76 | -1.85 |
| 65 | 791 | 9.0 | 48.35 | -123.4 | 72.0 | 18.88 | 24.0 | 20.02 | 33.85 |
| 64 | 847 | 7.0 | 48.29 | -123.92 | 76.0 | 23.42 | 20.88 | 23.46 | 42.8 |
| 74 | 1102 | 11.0 | 48.29 | -123.72 | 76.0 | 32.72 | 0.0 | 30.5 | 42.45 |
| 84 | 1485 | 7.0 | 48.35 | -123.56 | 72.0 | 59.53 | 47.76 | 49.53 | 40.55 |
| 61 | 955 | 11.0 | 47.9 | -124.32 | 76.0 | 71.29 | 34.0 | 75.14 | 41.03 |
| 14 | 582 | 4.0 |  |  |  |  |  |  |  |
| 20 | 796 | 5.0 |  |  |  |  |  |  |  |
| 108 | 1466 | 5.0 |  |  |  |  |  |  |  |
| 75 | 999 | 6.0 |  |  |  |  |  |  |  |
| 13 | 576 | 8.0 |  |  |  |  |  |  |  |
| 106 | 1152 | 6.0 |  |  |  |  |  |  |  |
| 114 | 327 | 5.0 |  |  |  |  |  |  |  |
| 71 | 281 |  |  |  |  |  |  |  |  |
| 30 | 175 | 6.0 |  |  |  |  |  |  |  |
| 38 | 631 | 7.0 |  |  |  |  |  |  |  |
| 52 | 1025 | 6.0 |  |  |  |  |  |  |  |
| 128 | 1205 | 9.0 |  |  |  |  |  |  |  |


## 8. Double-difference pass

`ZZdd` is the adopted solution after a Waldhauser-Ellsworth double-difference pass between neighbouring families, with the differential times measured by cross-correlating the two families' long stacks at the SAME station (`scripts/uu_dd.py`, observations in `results/u_locate/north/dd_ZZdd_obs.csv`).

| variant | n_located | epi_med | epi_p90 | dep_med | d3_med | rail | within5 | within10 | within20 |
|---|---|---|---|---|---|---|---|---|---|
| L9_xcdiff | 46 | 4.0 | 19.39 | 5.24 | 6.44 | 0.09 | 26 | 35 | 41 |
| ZZdd | 46 | 4.73 | 19.39 | 5.1 | 7.4 | 0.09 | 24 | 35 | 41 |


Acceptance, stated before the pass ran: adopt the double-difference solution only if it improves that table; otherwise it is reported and not used.


**Measured: NOT ADOPTED.**  Median 3-D distance 6.44 km (`L9_xcdiff`) -> 7.40 km (`ZZdd`), epicentral 4.00 -> 4.73 km, within 10 km 35 -> 35 of 58.  The pass does NOT improve the known answer, so the single-family solution stands and Washington is located with it.  The reason is visible in `logs/u_north_dd.log`: only 1,123 differential times over 179 linked pairs survive the cc >= 0.4 / >= 4 shared-station requirement out of 3,363 candidate pairs, so most families receive no differential constraint at all (median |move| 0.00 km) and the few that move do not move toward the published coordinates.


## 9. Washington -- relocated, with the identical table for the null and for Arm B

| pop | n_families | n_located | rail | med_nearest_3d | med_nearest_epi | p10_nearest_3d | boot_km_med | boot_km_p90 | ctrl_pass_frac | fam_within_5km | fam_within_5km_frac | fam_within_10km | fam_within_10km_frac | fam_within_20km | fam_within_20km_frac | coords_hit_5km | coords_hit_10km | coords_hit_20km | n_coords |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| real | 1479 | 788 | 0.53 | 40.98 | 19.23 | 11.43 | 35.64 | 113.71 | 0.37 | 12 | 0.02 | 65 | 0.08 | 175 | 0.22 | 16 | 61 | 69 | 69 |
|  | 1471 | 230 | 0.64 | 47.88 | 25.57 | 25.84 | 29.05 | 113.75 | 0.06 | 0 | 0.0 | 2 | 0.01 | 13 | 0.06 | 0 | 3 | 31 | 69 |
| armB | 739 | 730 | 0.69 | 47.21 | 22.03 | 19.82 | 29.68 | 95.87 | 0.56 | 10 | 0.01 | 25 | 0.03 | 74 | 0.1 | 9 | 46 | 69 | 69 |


`real` = the 2010 Olympic arrays-as-arrays discovery run (`results/o_oly2010/2010/`); `null` = its whole-chain site-shift null (`2010_null/`); `armB` = Chestler & Creager's OWN detection chain run literally on her own data (`chestler_arm/`, 739 families).  `chest_km_3d` = hypocentral distance to the nearest of the 69 published Olympic coordinates.  Ranked per-family tables with every criterion as a column: `results/u_locate/wa/wa_{real,null,armB}_ranked.csv`; real-vs-null separation curves in `wa_separation_curves.csv`.


## 9b. THE PRIORITY QUESTION -- the 280 richest real 2010-Olympic templates, relocated

`results/sf_final/table_oly2010_{real,null}.csv` report **280 real 2010-Olympic templates with >= 20 dedupe partners and ZERO null templates that rich** (the null maximum is 16); AUC 0.758 in the correct direction, surviving a detection-richness control at matched `n_det`.  That is the strongest candidate family population Washington has produced, so it is the population the rebuilt locator has to be pointed at.  Below they are relocated with the adopted locator and carry their own bootstrap radius, and they are placed beside the **null's best-ranked rows** -- the equal-size head of the null ordered by the SAME criterion, which is the matched comparator.  This filters nothing: all 1,485 real families keep their rows in `wa_real_ranked.csv` and `dedupe_n_partners` / `sf_top` are simply two more columns (memory `rank-dont-discard`).


**adopted locator** (`/home/jovyan/tremorferometry/lfe_array_discovery/results/u_locate/wa/`)

| pop | n_families | n_located | rail | boot_km_med | boot_km_p90 | med_nearest_3d | med_nearest_epi | p10_nearest_3d | fam_within_5km | fam_within_10km | fam_within_20km | coords_hit_5km | coords_hit_10km | coords_hit_20km | n_coords |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| real[dedupe_n_partners>=20] | 280 | 171 | 0.71 | 49.45 | 124.49 | 47.45 | 27.03 | 22.87 | 0.0 | 1.0 | 15.0 | 0.0 | 3.0 | 45.0 | 69.0 |
| null[top 280 by dedupe_n_partners] | 280 | 51 | 0.59 | 27.42 | 116.82 | 46.31 | 24.94 | 26.77 | 0.0 | 0.0 | 3.0 | 0.0 | 0.0 | 4.0 | 69.0 |
| null[dedupe_n_partners>=20] | 0 | 0 |  |  |  |  |  |  |  |  |  |  |  |  |  |
| real[all] | 1479 | 788 | 0.53 | 35.64 | 113.71 | 40.98 | 19.23 | 11.43 | 12.0 | 65.0 | 175.0 | 16.0 | 61.0 | 69.0 | 69.0 |
| null[all] | 1471 | 230 | 0.64 | 29.05 | 113.75 | 47.88 | 25.57 | 25.84 | 0.0 | 2.0 | 13.0 | 0.0 | 3.0 | 31.0 | 69.0 |


**coverage companion locator** (`/home/jovyan/tremorferometry/lfe_array_discovery/results/u_locate/wa_cov/`)

| pop | n_families | n_located | rail | boot_km_med | boot_km_p90 | med_nearest_3d | med_nearest_epi | p10_nearest_3d | fam_within_5km | fam_within_10km | fam_within_20km | coords_hit_5km | coords_hit_10km | coords_hit_20km | n_coords |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| real[dedupe_n_partners>=20] | 280 | 270 | 0.63 | 65.9 | 141.3 | 46.94 | 31.75 | 24.16 | 2.0 | 4.0 | 20.0 | 2.0 | 5.0 | 37.0 | 69.0 |
| null[top 280 by dedupe_n_partners] | 280 | 262 | 0.65 | 87.94 | 140.86 | 57.07 | 47.95 | 22.72 | 0.0 | 5.0 | 19.0 | 0.0 | 8.0 | 58.0 | 69.0 |
| null[dedupe_n_partners>=20] | 0 | 0 |  |  |  |  |  |  |  |  |  |  |  |  |  |
| real[all] | 1479 | 1412 | 0.55 | 58.46 | 138.13 | 44.28 | 29.27 | 15.31 | 17.0 | 73.0 | 218.0 | 16.0 | 63.0 | 69.0 | 69.0 |
| null[all] | 1471 | 1353 | 0.66 | 87.22 | 158.02 | 57.0 | 47.77 | 27.62 | 2.0 | 13.0 | 64.0 | 3.0 | 25.0 | 69.0 | 69.0 |


`pop` rows: `real[dedupe_n_partners>=20]` = the 280; `null[top 280 by dedupe_n_partners]` = the matched head of the null; `null[dedupe_n_partners>=20]` = empty by construction (no null template is that rich), kept as a row so the table is the identical table; `real[all]` / `null[all]` = the whole populations for scale.  `med_nearest_3d` = median hypocentral distance from a located family to the NEAREST of Chestler & Creager's 69 published Olympic coordinates; `coords_hit_k` = how many of those 69 have at least one family within k km, the reverse direction that a mere pile-up cannot fake.  `boot_km` = the per-family station-bootstrap epicentral radius, so a distance is only meaningful beside it.


Reading the three outcomes this distinguishes: if the 280 sit on Chestler's coordinates within their bootstrap radii and cover many of the 69, they ARE her families rediscovered; if they sit in the same depth/along-strike band but away from her coordinates by more than their radii, they are new families on the same interface; if their distance distribution is indistinguishable from the null's best-ranked rows, the dedupe-partner excess is not locating anything.  Per-family rows: `/home/jovyan/tremorferometry/lfe_array_discovery/results/u_locate/wa/wa_sf_top_real_located.csv`.


## 10. Files

| what | path |
|---|---|
| long stacks (44 s, 12 s pre-node) | `/home/jovyan/tremorferometry/lfe_array_discovery/results/u_locate/north_stacks.npz`, `/home/jovyan/tremorferometry/lfe_array_discovery/results/u_locate/wa_*_stacks.npz` |
| picks | `/home/jovyan/tremorferometry/lfe_array_discovery/results/u_locate/north/picks.parquet`, `/home/jovyan/tremorferometry/lfe_array_discovery/results/u_locate/wa_*/picks.parquet` |
| pick residuals vs published locations | `/home/jovyan/tremorferometry/lfe_array_discovery/results/u_locate/north/pick_residuals.csv` |
| every variant's locations | `/home/jovyan/tremorferometry/lfe_array_discovery/results/u_locate/north/loc_*.csv` |
| ladder + per-family errors | `/home/jovyan/tremorferometry/lfe_array_discovery/results/u_locate/north_eval/ladder.csv`, `per_family.csv` |
| adopted locator | `/home/jovyan/tremorferometry/lfe_array_discovery/results/u_locate/north_eval/adopted.json` |
| Washington ranked tables | `/home/jovyan/tremorferometry/lfe_array_discovery/results/u_locate/wa/wa_*_ranked.csv` |
| the 280 richest real templates, relocated | `/home/jovyan/tremorferometry/lfe_array_discovery/results/u_locate/wa/wa_sf_top_summary.csv`, `wa_sf_top_real_located.csv` |
| double-difference observations | `/home/jovyan/tremorferometry/lfe_array_discovery/results/u_locate/north/dd_ZZdd_obs.csv` |
