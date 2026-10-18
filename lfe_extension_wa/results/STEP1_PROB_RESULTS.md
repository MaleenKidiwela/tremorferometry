# STEP 1 — a probability for every repeat (local false-discovery rate, scrambled control as the noise model)

Written 2026-09-30 06:39Z by `scripts/step1_report.py` (numbers only; interpretation after review). Queue `scripts/step1_prob_queue.sh`; outputs `results/step1_prob/`; figure `figures/step1_reliability.png`.

## Summary (plain language)

Every candidate of the all-templates matched filter (≥ 8.0 MAD, nothing cut) gets a probability P(real) that it is a real repeat and not noise. The noise model is the scrambled control: the same matched filter run after shifting each site to a different time of the same day, so nothing lines up across sites except by chance. A classifier learns how data candidates differ from scrambled-control candidates (per-site correlations above each site's own noise level, HF ratio, how many candidates the family and its ≤ 3-km group have within ±1 h, …); the control's share, scaled by how many candidates the data has on tremor-free days, is the local fraction of noise. Fitted on the build half of V2008 + V2010; checked on the V2010 / V2008 test half and on A2011 (a separate episode).

- **V2010 test**: calibration slope own / ≤3 km / any = 0.65 / 0.58 / 0.54 (target 0.8–1.2); matched / ΣP own 1.29, ≤3 km 2.72. Scrambled control with P > 0.5: 0.004 (never-seen control realization 0.005); data with P > 0.5: 0.021. Quiet days: data P > 0.5 0.020 vs control 0.010. Expected count ΣP vs true array count per family-day: Spearman 0.49 (raw candidate count 0.46), ratio ΣP / truth 0.20.
- **A2011**: calibration slope own / ≤3 km / any = 0.54 / 0.55 / 0.50 (target 0.8–1.2); matched / ΣP own 1.05, ≤3 km 1.06. Scrambled control with P > 0.5: 0.014 (never-seen control realization 0.036); data with P > 0.5: 0.028. Quiet days: data P > 0.5 – vs control –. Expected count ΣP vs true array count per family-day: Spearman 0.48 (raw candidate count 0.40), ratio ΣP / truth 0.25.
- **V2008 test**: calibration slope own / ≤3 km / any = 0.56 / 0.57 / 0.53 (target 0.8–1.2); matched / ΣP own 1.55, ≤3 km 1.56. Scrambled control with P > 0.5: 0.007 (never-seen control realization 0.007); data with P > 0.5: 0.029. Quiet days: data P > 0.5 0.001 vs control 0.008. Expected count ΣP vs true array count per family-day: Spearman 0.45 (raw candidate count 0.43), ratio ΣP / truth 0.10.

Caveat for reading calibration: the array catalogue is the truth but is itself incomplete, and a real LFE of a neighbouring family counts as 'real' for P but matches only at the ≤ 3-km / any levels. Matched fraction below P can therefore mean real-but-uncatalogued; matched fraction ABOVE P means P is conservative.

## Method

- Candidates: `results/mf_mt1_all` (T1 all ≥ 30-event build templates, w-weighted, 42 stations, `--drop-dead 0.05`), record tier ≥ 8.0 MAD, series data + control + control2; plus `results/mf_step1c34` = the same run with seed 20261001 → two more independent site-rolled scrambled controls (its data series is identical to mt1_all's: checked per day). K = 4 scrambled-control realizations.
- Features (all computable after the array era): `net_mad`, `net_cc`, `n_sta`, `n_sites`, `wsum`, `site_cc_max`, `site_cc_2nd`, `site_z_max`, `site_z_2nd`, `site_z_min`, `n_z2`, `n_z3`, `n_z5`, `frac_z3`, `dom_top`, `hfr_rel`, `dens_own_1h`, `dens_own_10m`, `dens_grp_1h`, `site_noise_mad_med`, `region_code`. Site values are the within-site w-weighted station cc at the candidate's aligned position (scrambled control: at the site's rolled position, exactly as the MF built it); the recomputed network cc equals the MF's to ≤ 2e-6 on every day. Site z = (site cc − median)/(1.4826 MAD) over 1000 random valid origins of the same family-day. `hfr_rel` = log10 E(10–19)/E(2–8) at the family-day's HF station minus the median of 300 random windows of that station-day. Densities count candidates of the same series.
- Density ratio: HistGradientBoosting + isotonic (3-fold by day) classifier data-vs-control on the fit set; f0/f = [(1−p)/p]·N_data/N_ctl.
- π0 rule: on quiet days (V2010 'quiet', V2008 'post'; tremor-free by the window design) every data candidate is noise → r_q = N_data(quiet) / (N_ctl(quiet)/K) per region (global if < 500 quiet data candidates), π0 = r_q·(N_ctl/K)/N_data, so P = 1 − r_q·(1−p)/(p·K), clipped to [0, 1]. Real LFEs on quiet days count as noise → conservative.
- Fit set: BUILD half of V2008 + V2010 (every series). Held out: V2010 / V2008 TEST half, A2011 (all). Leave-one-control-out refits (K = 3) give the model spread of the rate product and a P for a control realization the model never saw. Matching: within ±2 s (nearest, not one-to-one) of an array detection of the own family / any family ≤ 3 km (lat_use/lon_use) / any family, only on days the family's own array run scanned.

| region | quiet data candidates | quiet control candidates (all K) | r_q | source |
|---|---|---|---|---|
| NW | 7633 | 29221 | 1.045 | region |
| SW | 63 | 257 | 1.020 | global |
| Olympic | 53232 | 209575 | 1.016 | region |
| Arm B | 870 | 3243 | 1.073 | region |
| ALL | 61798 | 242296 | 1.020 | global |

Feature ranking (fit set; permutation log-loss gain on 200k candidates; univariate AUC data vs control):

| feature | perm. log-loss gain | AUC data vs control | NaN frac |
|---|---|---|---|
| n_z2 | 0.0417 | 0.567 | 0.000 |
| dom_top | 0.0326 | 0.431 | 0.000 |
| net_mad | 0.0081 | 0.514 | 0.000 |
| dens_grp_1h | 0.0050 | 0.469 | 0.000 |
| region_code | 0.0048 | 0.464 | 0.000 |
| dens_own_10m | 0.0035 | 0.501 | 0.000 |
| net_cc | 0.0018 | 0.476 | 0.000 |
| dens_own_1h | 0.0016 | 0.479 | 0.000 |
| n_sta | 0.0015 | 0.500 | 0.000 |
| site_cc_max | 0.0012 | 0.474 | 0.000 |
| hfr_rel | 0.0011 | 0.482 | 0.000 |
| wsum | 0.0011 | 0.463 | 0.000 |
| n_sites | 0.0010 | 0.505 | 0.000 |
| site_noise_mad_med | 0.0006 | 0.514 | 0.000 |
| site_z_max | 0.0006 | 0.455 | 0.000 |
| site_z_min | 0.0005 | 0.544 | 0.000 |
| frac_z3 | 0.0004 | 0.550 | 0.000 |
| site_z_2nd | 0.0004 | 0.561 | 0.000 |
| site_cc_2nd | 0.0002 | 0.533 | 0.000 |
| n_z3 | 0.0001 | 0.548 | 0.000 |
| n_z5 | 0.0000 | 0.472 | 0.000 |

## (a) Calibration — fraction matched to the array catalogue by P bin

Slope = weighted least-squares slope of matched fraction on bin-mean P (bins with ≥ 20 candidates); matched/ΣP = calibration in the large. Control rows are the chance level.

| set | region | series | n | ΣP | slope own · 3km · any | intercept own · 3km · any | matched/ΣP own · 3km · any |
|---|---|---|---|---|---|---|---|
| V2010_test | ALL | data | 115301 | 6031.5 | 0.65 · 0.58 · 0.54 | 0.03 · 0.11 · 0.24 | 1.29 · 2.72 · 5.08 |
| V2010_test | ALL | control | 447017 | 15035.5 | 0.01 · 0.05 · -0.01 | 0.00 · 0.01 · 0.08 | 0.06 · 0.47 · 2.23 |
| V2010_test | NW | data | 12527 | 2471.1 | 0.56 · 0.56 · 0.47 | 0.18 · 0.18 · 0.38 | 1.48 · 1.49 · 2.39 |
| V2010_test | NW | control | 38294 | 1349.2 | 0.01 · 0.01 · 0.07 | 0.01 · 0.01 · 0.11 | 0.18 · 0.20 · 3.06 |
| V2010_test | Olympic | data | 101356 | 3481.8 | 0.21 · 0.19 · 0.06 | 0.03 · 0.12 · 0.24 | 1.15 · 3.63 · 6.96 |
| V2010_test | Olympic | control | 402955 | 13512.4 | 0.00 · 0.08 · -0.06 | 0.00 · 0.01 · 0.07 | 0.04 · 0.50 · 2.15 |
| V2010_test | Arm B | data | 1374 | 53.1 | 0.60 · 0.60 · 0.33 | 0.03 · 0.03 · 0.34 | 1.38 · 1.47 · 9.00 |
| V2010_test | Arm B | control | 5666 | 115.6 | 0.04 · 0.04 · 0.03 | 0.00 · 0.01 · 0.08 | 0.22 · 0.41 · 3.72 |
| A2011 | ALL | data | 61999 | 8464.7 | 0.54 · 0.55 · 0.50 | 0.07 · 0.07 · 0.19 | 1.05 · 1.06 · 1.90 |
| A2011 | ALL | control | 222561 | 15846.9 | 0.00 · 0.00 · -0.00 | 0.00 · 0.00 · 0.05 | 0.04 · 0.05 · 0.68 |
| A2011 | NW | data | 61999 | 8464.7 | 0.54 · 0.55 · 0.50 | 0.07 · 0.07 · 0.19 | 1.05 · 1.06 · 1.90 |
| A2011 | NW | control | 222561 | 15846.9 | 0.00 · 0.00 · -0.00 | 0.00 · 0.00 · 0.05 | 0.04 · 0.05 · 0.68 |
| V2008_test | ALL | data | 13943 | 2494.6 | 0.56 · 0.57 · 0.53 | 0.18 · 0.18 · 0.31 | 1.55 · 1.56 · 2.25 |
| V2008_test | ALL | control | 44129 | 1327.9 | 0.01 · 0.01 · 0.05 | 0.00 · 0.00 · 0.08 | 0.15 · 0.16 · 2.70 |
| V2008_test | NW | data | 13647 | 2333.4 | 0.55 · 0.56 · 0.52 | 0.18 · 0.18 · 0.31 | 1.59 · 1.60 · 2.34 |
| V2008_test | NW | control | 43591 | 1252.8 | 0.01 · 0.01 · 0.06 | 0.00 · 0.00 · 0.08 | 0.16 · 0.17 · 2.83 |
| V2008_test | SW | data | 232 | 120.7 | 0.68 · 0.70 · 0.63 | 0.10 · 0.11 · 0.20 | 0.87 · 0.90 · 1.02 |
| V2008_test | SW | control | 481 | 35.5 | – · – · – | – · – · – | 0.06 · 0.06 · 1.21 |
| fit_in_sample | ALL | data | 160191 | 40320.5 | 0.97 · 0.89 · 0.76 | 0.00 · 0.08 · 0.22 | 0.98 · 1.22 · 1.63 |
| fit_in_sample | ALL | control | 501947 | 14071.4 | 0.01 · 0.02 · -0.01 | 0.00 · 0.02 · 0.08 | 0.08 · 0.59 · 3.01 |
| fit_in_sample | NW | data | 43208 | 23724.5 | 0.85 · 0.85 · 0.69 | 0.11 · 0.11 · 0.29 | 1.05 · 1.05 · 1.22 |
| fit_in_sample | NW | control | 84142 | 2458.4 | -0.00 · -0.00 · 0.04 | 0.00 · 0.01 · 0.09 | 0.16 · 0.18 · 3.28 |
| fit_in_sample | SW | data | 6814 | 6699.1 | 0.87 · 0.87 · 0.77 | 0.13 · 0.13 · 0.22 | 1.00 · 1.00 · 1.00 |
| fit_in_sample | SW | control | 471 | 39.3 | – · – · – | – · – · – | 0.00 · 0.03 · 0.99 |
| fit_in_sample | Olympic | data | 105984 | 7085.1 | 0.96 · 0.88 · 0.73 | -0.02 · 0.08 · 0.20 | 0.73 · 2.05 · 3.79 |
| fit_in_sample | Olympic | control | 412112 | 11411.2 | 0.01 · 0.05 · -0.09 | 0.00 · 0.02 · 0.08 | 0.06 · 0.67 · 2.95 |
| fit_in_sample | Arm B | data | 4115 | 2772.9 | 0.96 · 0.95 · 0.61 | 0.04 · 0.05 · 0.39 | 1.02 · 1.02 · 1.19 |
| fit_in_sample | Arm B | control | 5093 | 78.6 | – · – · – | – · – · – | 1.08 · 1.70 · 7.77 |

Reliability table, V2010_test, ALL regions:

| P bin | data n | mean P | frac own | frac ≤3 km | frac any | control n | control frac any |
|---|---|---|---|---|---|---|---|
| 0.0-0.1 | 101221 | 0.019 | 0.050 | 0.128 | 0.254 | 404415 | 0.076 |
| 0.1-0.2 | 9132 | 0.144 | 0.076 | 0.151 | 0.244 | 34444 | 0.066 |
| 0.2-0.3 | 1930 | 0.226 | 0.124 | 0.184 | 0.298 | 6440 | 0.071 |
| 0.3-0.4 | 438 | 0.344 | 0.196 | 0.242 | 0.356 | 604 | 0.070 |
| 0.4-0.5 | 166 | 0.452 | 0.235 | 0.289 | 0.452 | 240 | 0.087 |
| 0.5-0.6 | 161 | 0.548 | 0.466 | 0.503 | 0.596 | 213 | 0.146 |
| 0.6-0.7 | 163 | 0.653 | 0.466 | 0.491 | 0.693 | 154 | 0.182 |
| 0.7-0.8 | 168 | 0.754 | 0.482 | 0.500 | 0.643 | 138 | 0.130 |
| 0.8-0.9 | 236 | 0.856 | 0.568 | 0.589 | 0.708 | 127 | 0.134 |
| 0.9-1.0 | 1686 | 0.979 | 0.730 | 0.735 | 0.835 | 242 | 0.099 |

Reliability table, A2011, ALL regions:

| P bin | data n | mean P | frac own | frac ≤3 km | frac any | control n | control frac any |
|---|---|---|---|---|---|---|---|
| 0.0-0.1 | 40454 | 0.023 | 0.103 | 0.104 | 0.231 | 164686 | 0.051 |
| 0.1-0.2 | 12075 | 0.151 | 0.090 | 0.091 | 0.186 | 42235 | 0.038 |
| 0.2-0.3 | 2717 | 0.233 | 0.130 | 0.131 | 0.230 | 7987 | 0.045 |
| 0.3-0.4 | 973 | 0.345 | 0.183 | 0.186 | 0.263 | 2745 | 0.046 |
| 0.4-0.5 | 537 | 0.455 | 0.248 | 0.250 | 0.363 | 1225 | 0.060 |
| 0.5-0.6 | 553 | 0.555 | 0.354 | 0.356 | 0.450 | 988 | 0.067 |
| 0.6-0.7 | 597 | 0.655 | 0.402 | 0.409 | 0.519 | 1002 | 0.058 |
| 0.7-0.8 | 610 | 0.752 | 0.474 | 0.475 | 0.572 | 673 | 0.091 |
| 0.8-0.9 | 569 | 0.854 | 0.503 | 0.508 | 0.603 | 518 | 0.060 |
| 0.9-1.0 | 2914 | 0.974 | 0.671 | 0.682 | 0.758 | 502 | 0.060 |

Reliability table, V2008_test, ALL regions:

| P bin | data n | mean P | frac own | frac ≤3 km | frac any | control n | control frac any |
|---|---|---|---|---|---|---|---|
| 0.0-0.1 | 10328 | 0.016 | 0.171 | 0.172 | 0.300 | 41466 | 0.080 |
| 0.1-0.2 | 924 | 0.153 | 0.398 | 0.400 | 0.543 | 1654 | 0.086 |
| 0.2-0.3 | 269 | 0.230 | 0.435 | 0.442 | 0.591 | 337 | 0.095 |
| 0.3-0.4 | 114 | 0.350 | 0.553 | 0.553 | 0.675 | 126 | 0.119 |
| 0.4-0.5 | 82 | 0.460 | 0.537 | 0.537 | 0.659 | 88 | 0.102 |
| 0.5-0.6 | 93 | 0.552 | 0.462 | 0.462 | 0.645 | 72 | 0.097 |
| 0.6-0.7 | 120 | 0.648 | 0.575 | 0.575 | 0.692 | 78 | 0.090 |
| 0.7-0.8 | 136 | 0.751 | 0.537 | 0.551 | 0.669 | 53 | 0.170 |
| 0.8-0.9 | 188 | 0.854 | 0.585 | 0.590 | 0.707 | 76 | 0.079 |
| 0.9-1.0 | 1689 | 0.981 | 0.710 | 0.721 | 0.804 | 179 | 0.140 |

![reliability](../figures/step1_reliability.png)

## (b) Scrambled control's own P, (c) quiet vs episode days, (d) per region

`ctl P>0.5 (never seen)` = P of each control realization from the leave-that-control-out model. ΣP per realization = expected 'real' count the method would wrongly claim in one scrambled network.

| set | region | days | data n | data P>0.5 | data P>0.9 | data ΣP | ctl n / real. | ctl P>0.1 | ctl P>0.5 | ctl P>0.5 (never seen) | ctl ΣP / real. | ctl ΣP / real. (never seen) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V2010_test | ALL | all | 172854 | 0.021 | 0.011 | 8795.9 | 168969 | 0.088 | 0.004 | 0.005 | 5952.3 | 6171.95 |
| V2010_test | ALL | quiet | 48576 | 0.020 | 0.001 | 2305.9 | 48451 | 0.076 | 0.010 | 0.010 | 1929.0 | 1903.21 |
| V2010_test | ALL | episode | 124278 | 0.021 | 0.015 | 6490.1 | 120518 | 0.093 | 0.002 | 0.003 | 4023.3 | 4268.74 |
| V2010_test | NW | all | 21800 | 0.106 | 0.076 | 2723.3 | 19002 | 0.063 | 0.011 | 0.015 | 561.5 | 628.03 |
| V2010_test | NW | quiet | 5429 | 0.006 | 0.002 | 133.9 | 5521 | 0.047 | 0.004 | 0.008 | 130.1 | 148.95 |
| V2010_test | NW | episode | 16371 | 0.140 | 0.101 | 2589.4 | 13482 | 0.070 | 0.013 | 0.018 | 431.3 | 479.07 |
| V2010_test | SW | all | 381 | 0.386 | 0.289 | 158.7 | 168 | 0.200 | 0.024 | 0.042 | 12.9 | 15.56 |
| V2010_test | SW | quiet | 33 | 0.061 | 0.061 | 3.7 | 44 | 0.169 | 0.023 | 0.034 | 3.1 | 3.71 |
| V2010_test | SW | episode | 348 | 0.417 | 0.310 | 155.0 | 124 | 0.211 | 0.024 | 0.045 | 9.8 | 11.85 |
| V2010_test | Olympic | all | 148866 | 0.007 | 0.000 | 5827.5 | 147917 | 0.092 | 0.003 | 0.004 | 5328.6 | 5464.23 |
| V2010_test | Olympic | quiet | 42747 | 0.022 | 0.000 | 2160.4 | 42469 | 0.080 | 0.011 | 0.010 | 1786.1 | 1738.41 |
| V2010_test | Olympic | episode | 106119 | 0.001 | 0.000 | 3667.1 | 105448 | 0.096 | 0.000 | 0.001 | 3542.5 | 3725.82 |
| V2010_test | Arm B | all | 1807 | 0.032 | 0.015 | 86.3 | 1882 | 0.065 | 0.010 | 0.016 | 49.3 | 64.14 |
| V2010_test | Arm B | quiet | 367 | 0.008 | 0.005 | 7.8 | 417 | 0.060 | 0.007 | 0.012 | 9.6 | 12.13 |
| V2010_test | Arm B | episode | 1440 | 0.037 | 0.018 | 78.5 | 1465 | 0.067 | 0.011 | 0.017 | 39.7 | 52.00 |
| A2011 | ALL | all | 294824 | 0.028 | 0.010 | 24821.5 | 287372 | 0.264 | 0.014 | 0.036 | 20152.4 | 29780.87 |
| A2011 | ALL | episode | 294824 | 0.028 | 0.010 | 24821.5 | 287372 | 0.264 | 0.014 | 0.036 | 20152.4 | 29780.87 |
| A2011 | NW | all | 61999 | 0.085 | 0.047 | 8464.7 | 55640 | 0.260 | 0.017 | 0.032 | 3961.7 | 5136.51 |
| A2011 | NW | episode | 61999 | 0.085 | 0.047 | 8464.7 | 55640 | 0.260 | 0.017 | 0.032 | 3961.7 | 5136.51 |
| A2011 | SW | all | 709 | 0.477 | 0.134 | 362.9 | 700 | 0.989 | 0.488 | 0.657 | 369.1 | 418.78 |
| A2011 | SW | episode | 709 | 0.477 | 0.134 | 362.9 | 700 | 0.989 | 0.488 | 0.657 | 369.1 | 418.78 |
| A2011 | Olympic | all | 215645 | 0.011 | 0.000 | 14677.8 | 214685 | 0.260 | 0.011 | 0.036 | 14610.2 | 22599.98 |
| A2011 | Olympic | episode | 215645 | 0.011 | 0.000 | 14677.8 | 214685 | 0.260 | 0.011 | 0.036 | 14610.2 | 22599.98 |
| A2011 | Arm B | all | 16471 | 0.023 | 0.003 | 1316.1 | 16348 | 0.298 | 0.016 | 0.028 | 1211.4 | 1625.61 |
| A2011 | Arm B | episode | 16471 | 0.023 | 0.003 | 1316.1 | 16348 | 0.298 | 0.016 | 0.028 | 1211.4 | 1625.61 |
| V2008_test | ALL | all | 93754 | 0.029 | 0.018 | 5110.6 | 90830 | 0.080 | 0.007 | 0.007 | 3437.9 | 3507.20 |
| V2008_test | ALL | quiet | 13913 | 0.001 | 0.000 | 359.3 | 13122 | 0.123 | 0.008 | 0.003 | 527.1 | 570.04 |
| V2008_test | ALL | episode | 79841 | 0.033 | 0.021 | 4751.3 | 77708 | 0.073 | 0.007 | 0.008 | 2910.8 | 2937.15 |
| V2008_test | NW | all | 16322 | 0.128 | 0.097 | 2411.3 | 13546 | 0.053 | 0.008 | 0.013 | 386.8 | 423.07 |
| V2008_test | NW | quiet | 1993 | 0.003 | 0.001 | 53.1 | 1966 | 0.041 | 0.004 | 0.008 | 56.3 | 62.65 |
| V2008_test | NW | episode | 14329 | 0.145 | 0.111 | 2358.3 | 11580 | 0.055 | 0.009 | 0.013 | 330.6 | 360.42 |
| V2008_test | SW | all | 305 | 0.518 | 0.361 | 162.5 | 145 | 0.360 | 0.071 | 0.100 | 20.1 | 22.79 |
| V2008_test | SW | quiet | 9 | 0.111 | 0.000 | 1.3 | 11 | 0.302 | 0.070 | 0.093 | 1.3 | 1.63 |
| V2008_test | SW | episode | 296 | 0.530 | 0.372 | 161.2 | 134 | 0.364 | 0.071 | 0.100 | 18.8 | 21.16 |
| V2008_test | Olympic | all | 73842 | 0.006 | 0.000 | 2477.6 | 73896 | 0.087 | 0.007 | 0.006 | 2978.9 | 2997.80 |
| V2008_test | Olympic | quiet | 11474 | 0.000 | 0.000 | 295.7 | 10724 | 0.142 | 0.009 | 0.002 | 462.8 | 497.82 |
| V2008_test | Olympic | episode | 62368 | 0.007 | 0.000 | 2181.9 | 63171 | 0.077 | 0.007 | 0.007 | 2516.1 | 2499.97 |
| V2008_test | Arm B | all | 3285 | 0.001 | 0.000 | 59.1 | 3242 | 0.027 | 0.001 | 0.002 | 52.0 | 63.54 |
| V2008_test | Arm B | quiet | 437 | 0.000 | 0.000 | 9.3 | 420 | 0.019 | 0.002 | 0.004 | 6.7 | 7.94 |
| V2008_test | Arm B | episode | 2848 | 0.001 | 0.000 | 49.9 | 2822 | 0.028 | 0.001 | 0.002 | 45.3 | 55.60 |
| fit_in_sample | ALL | all | 297796 | 0.122 | 0.110 | 47167.2 | 257759 | 0.053 | 0.002 | 0.003 | 7332.5 | 7961.96 |
| fit_in_sample | ALL | quiet | 61798 | 0.010 | 0.001 | 3148.6 | 60574 | 0.051 | 0.002 | 0.003 | 1744.5 | 1864.43 |
| fit_in_sample | ALL | episode | 235998 | 0.151 | 0.139 | 44018.5 | 197185 | 0.053 | 0.002 | 0.004 | 5588.0 | 6097.52 |
| fit_in_sample | NW | all | 55588 | 0.421 | 0.389 | 24174.2 | 32886 | 0.058 | 0.009 | 0.014 | 924.9 | 1083.58 |
| fit_in_sample | NW | quiet | 7633 | 0.006 | 0.002 | 262.7 | 7305 | 0.045 | 0.005 | 0.009 | 180.5 | 217.76 |
| fit_in_sample | NW | episode | 47955 | 0.487 | 0.451 | 23911.5 | 25581 | 0.062 | 0.011 | 0.016 | 744.4 | 865.83 |
| fit_in_sample | SW | all | 7337 | 0.939 | 0.924 | 6910.2 | 349 | 0.298 | 0.036 | 0.069 | 34.8 | 44.99 |
| fit_in_sample | SW | quiet | 63 | 0.000 | 0.000 | 4.0 | 64 | 0.253 | 0.031 | 0.070 | 5.4 | 7.87 |
| fit_in_sample | SW | episode | 7274 | 0.948 | 0.932 | 6906.2 | 285 | 0.309 | 0.037 | 0.069 | 29.3 | 37.12 |
| fit_in_sample | Olympic | all | 226905 | 0.014 | 0.008 | 13206.8 | 219660 | 0.052 | 0.001 | 0.002 | 6291.5 | 6718.20 |
| fit_in_sample | Olympic | quiet | 53232 | 0.010 | 0.000 | 2851.9 | 52394 | 0.051 | 0.001 | 0.002 | 1544.3 | 1619.36 |
| fit_in_sample | Olympic | episode | 173673 | 0.015 | 0.010 | 10354.9 | 167266 | 0.052 | 0.000 | 0.001 | 4747.2 | 5098.84 |
| fit_in_sample | Arm B | all | 7966 | 0.349 | 0.334 | 2875.9 | 4863 | 0.032 | 0.004 | 0.009 | 81.4 | 115.18 |
| fit_in_sample | Arm B | quiet | 870 | 0.013 | 0.003 | 30.0 | 811 | 0.043 | 0.003 | 0.007 | 14.3 | 19.45 |
| fit_in_sample | Arm B | episode | 7096 | 0.390 | 0.375 | 2845.9 | 4052 | 0.030 | 0.004 | 0.009 | 67.1 | 95.73 |
| ALL | ALL | all | 867908 | 0.059 | 0.046 | 86381.7 | 814130 | 0.138 | 0.007 | 0.016 | 37186.1 | 47754.87 |
| ALL | ALL | quiet | 126292 | 0.013 | 0.001 | 5904.0 | 124542 | 0.069 | 0.006 | 0.006 | 4277.8 | 4415.44 |
| ALL | ALL | episode | 741616 | 0.067 | 0.053 | 80477.7 | 689588 | 0.150 | 0.007 | 0.018 | 32908.3 | 43339.43 |
| ALL | NW | all | 156925 | 0.212 | 0.178 | 37931.7 | 122218 | 0.150 | 0.013 | 0.022 | 5869.4 | 7310.13 |
| ALL | NW | quiet | 15254 | 0.006 | 0.002 | 455.4 | 15052 | 0.046 | 0.005 | 0.008 | 373.8 | 437.30 |
| ALL | NW | episode | 141671 | 0.234 | 0.197 | 37476.3 | 107166 | 0.165 | 0.014 | 0.024 | 5495.5 | 6872.83 |
| ALL | SW | all | 8754 | 0.862 | 0.811 | 7605.3 | 1373 | 0.644 | 0.268 | 0.369 | 437.7 | 503.23 |
| ALL | SW | quiet | 106 | 0.028 | 0.019 | 9.1 | 122 | 0.224 | 0.031 | 0.060 | 10.0 | 13.38 |
| ALL | SW | episode | 8648 | 0.872 | 0.821 | 7596.2 | 1251 | 0.685 | 0.292 | 0.399 | 427.7 | 489.85 |
| ALL | Olympic | all | 672541 | 0.010 | 0.003 | 36503.4 | 664021 | 0.132 | 0.005 | 0.014 | 29481.6 | 38068.64 |
| ALL | Olympic | quiet | 109234 | 0.014 | 0.000 | 5391.9 | 107682 | 0.072 | 0.006 | 0.005 | 3862.3 | 3924.00 |
| ALL | Olympic | episode | 563307 | 0.010 | 0.003 | 31111.5 | 556340 | 0.144 | 0.005 | 0.015 | 25619.3 | 34144.65 |
| ALL | Arm B | all | 29688 | 0.108 | 0.092 | 4341.3 | 26518 | 0.198 | 0.011 | 0.020 | 1397.4 | 1872.87 |
| ALL | Arm B | quiet | 1698 | 0.008 | 0.003 | 47.5 | 1687 | 0.041 | 0.004 | 0.008 | 31.7 | 40.76 |
| ALL | Arm B | episode | 27990 | 0.114 | 0.097 | 4293.8 | 24831 | 0.208 | 0.012 | 0.021 | 1365.7 | 1832.10 |

By PNSN regional tremor-day label (label only):

| set | region | tremor day | family-days | data ΣP / fam-day | ctl ΣP / fam-day (per real.) | data P>0.5 | ctl P>0.5 |
|---|---|---|---|---|---|---|---|
| V2010_test | ALL | False | 884 | 2.619 | 2.193 | 0.020 | 0.010 |
| V2010_test | ALL | True | 2690 | 2.409 | 1.492 | 0.021 | 0.002 |
| V2010_test | NW | False | 551 | 0.259 | 0.252 | 0.006 | 0.004 |
| V2010_test | NW | True | 1596 | 1.617 | 0.265 | 0.143 | 0.013 |
| V2010_test | SW | False | 41 | 0.096 | 0.099 | 0.047 | 0.028 |
| V2010_test | SW | True | 204 | 0.759 | 0.044 | 0.429 | 0.022 |
| V2010_test | Olympic | False | 244 | 8.854 | 7.320 | 0.022 | 0.011 |
| V2010_test | Olympic | True | 720 | 5.093 | 4.920 | 0.001 | 0.000 |
| V2010_test | Arm B | False | 48 | 0.163 | 0.201 | 0.008 | 0.007 |
| V2010_test | Arm B | True | 170 | 0.462 | 0.233 | 0.037 | 0.011 |
| A2011 | ALL | False | 172 | 1.496 | 1.543 | 0.047 | 0.057 |
| A2011 | ALL | True | 3962 | 6.200 | 5.019 | 0.028 | 0.014 |
| A2011 | NW | False | 68 | 2.946 | 3.065 | 0.027 | 0.037 |
| A2011 | NW | True | 2662 | 3.105 | 1.410 | 0.086 | 0.016 |
| A2011 | SW | False | 104 | 0.547 | 0.548 | 0.186 | 0.204 |
| A2011 | SW | True | 213 | 1.436 | 1.465 | 0.573 | 0.578 |
| A2011 | Olympic | True | 724 | 20.273 | 20.180 | 0.011 | 0.011 |
| A2011 | Arm B | True | 363 | 3.626 | 3.337 | 0.023 | 0.016 |
| V2008_test | ALL | False | 2788 | 1.833 | 1.233 | 0.029 | 0.007 |
| V2008_test | NW | False | 1667 | 1.447 | 0.232 | 0.128 | 0.008 |
| V2008_test | SW | False | 186 | 0.874 | 0.108 | 0.518 | 0.071 |
| V2008_test | Olympic | False | 712 | 3.480 | 4.184 | 0.006 | 0.007 |
| V2008_test | Arm B | False | 223 | 0.265 | 0.233 | 0.001 | 0.001 |

## (e) Expected-count check — ΣP per family-day vs the true array count

Episode family-days the family's array run scanned; scored half only (test half; A2011 all). truth ≤3 km = array detections of the ≤ 3-km group merged within ±2 s. `ΣP ctl` = the same sum on one scrambled-control realization.

| set | region | family-days | ΣP | ΣP ctl | truth own | truth ≤3 km | ΣP/own | ΣP/≤3 km | Spearman own | Spearman ≤3 km | Spearman raw count vs own |
|---|---|---|---|---|---|---|---|---|---|---|---|
| V2010_test | ALL | 4572 | 6031.5 | 3754.1 | 30921 | 147454 | 0.195 | 0.041 | 0.49 | 0.55 | 0.46 |
| V2010_test | NW | 2640 | 2471.1 | 333.3 | 17082 | 37188 | 0.145 | 0.066 | 0.51 | 0.47 | 0.48 |
| V2010_test | Olympic | 784 | 3481.8 | 3381.7 | 7056 | 85199 | 0.493 | 0.041 | 0.37 | 0.31 | 0.28 |
| V2010_test | Arm B | 1148 | 78.5 | 39.1 | 6783 | 25067 | 0.012 | 0.003 | 0.18 | 0.09 | 0.19 |
| V2008_test | ALL | 4934 | 2494.6 | 318.6 | 23962 | 44153 | 0.104 | 0.057 | 0.45 | 0.40 | 0.43 |
| V2008_test | NW | 3036 | 2333.4 | 306.5 | 16129 | 25449 | 0.145 | 0.092 | 0.52 | 0.49 | 0.49 |
| V2008_test | SW | 1898 | 161.2 | 12.1 | 7833 | 18704 | 0.021 | 0.009 | 0.28 | 0.24 | 0.25 |
| A2011 | ALL | 3828 | 8464.7 | 3863.0 | 33306 | 50539 | 0.254 | 0.167 | 0.48 | 0.45 | 0.40 |
| A2011 | NW | 3828 | 8464.7 | 3863.0 | 33306 | 50539 | 0.254 | 0.167 | 0.48 | 0.45 | 0.40 |

## Rate product — expected real repeats per family and day

`results/step1_prob/rate_family_day.parquet`: one row per observable family-day (all 97 days, V2008 + V2010 + A2011): `sumP` (whole day; build-half P is in-sample), `sumP_scored_half`, `sd_total` = √(Σ P(1−P) + model spread²) with the model spread = SD of ΣP over the 4 leave-one-control-out refits, `sumP_ctl_per_real` (the same sum on one scrambled network), truth counts where the family's array run scanned the day. Tremor-day label is a column only.

| region | days | family-days | mean ΣP / fam-day | mean sd | mean ΣP ctl / fam-day | family-days with ΣP ≥ 5 |
|---|---|---|---|---|---|---|
| Olympic | episode | 2352 | 13.23 | 6.21 | 10.89 | 1571 |
| Olympic | post | 84 | 10.93 | 4.43 | 8.97 | 64 |
| Olympic | quiet | 280 | 15.98 | 5.00 | 11.10 | 221 |
| Arm B | episode | 3444 | 1.25 | 0.61 | 0.39 | 230 |
| Arm B | post | 123 | 0.12 | 0.52 | 0.08 | 0 |
| Arm B | quiet | 410 | 0.08 | 0.39 | 0.05 | 1 |
| NW | episode | 11088 | 3.38 | 0.65 | 0.48 | 1471 |
| NW | post | 396 | 0.33 | 0.38 | 0.27 | 5 |
| NW | quiet | 1320 | 0.25 | 0.32 | 0.19 | 18 |
| SW | episode | 6132 | 1.24 | 0.28 | 0.05 | 327 |
| SW | post | 219 | 0.02 | 0.23 | 0.00 | 0 |
| SW | quiet | 730 | 0.01 | 0.17 | 0.00 | 0 |

Files: `candidates.parquet` (every candidate, every feature, P, P_loco_c1..c4, matches), `calibration.csv`, `calibration_slopes.csv`, `control_and_quiet.csv`, `tremor_label_split.csv`, `expected_count.csv`, `feature_importance.csv`, `model_info.json`, `rate_family_day.parquet`.
