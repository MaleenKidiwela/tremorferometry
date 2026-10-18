# STEP 1 — a probability for every repeat (local false-discovery rate, scrambled control as the noise model)

Written 2026-09-30 06:18Z by `scripts/step1_report.py` (numbers only; interpretation after review). Queue `scripts/step1_prob_queue.sh`; outputs `results/step1_prob/`; figure `figures/step1_reliability.png`.

## Summary (plain language)

Every candidate of the all-templates matched filter (≥ 8.0 MAD, nothing cut) gets a probability P(real) that it is a real repeat and not noise. The noise model is the scrambled control: the same matched filter run after shifting each site to a different time of the same day, so nothing lines up across sites except by chance. A classifier learns how data candidates differ from scrambled-control candidates (per-site correlations above each site's own noise level, HF ratio, how many candidates the family and its ≤ 3-km group have within ±1 h, …); the control's share, scaled by how many candidates the data has on tremor-free days, is the local fraction of noise. Fitted on the build half of V2008 + V2010; checked on the V2010 / V2008 test half and on A2011 (a separate episode).

- **V2010 test**: calibration slope own / ≤3 km / any = 0.80 / 0.67 / 0.65 (target 0.8–1.2); matched / ΣP own 4.13, ≤3 km 9.12. Scrambled control with P > 0.5: 0.003 (never-seen control realization 0.005); data with P > 0.5: 0.017. Quiet days: data P > 0.5 0.001 vs control 0.004. Expected count ΣP vs true array count per family-day: Spearman 0.11 (raw candidate count 0.10), ratio ΣP / truth 0.01.
- **A2011**: calibration slope own / ≤3 km / any = 0.61 / 0.62 / 0.48 (target 0.8–1.2); matched / ΣP own 2.24, ≤3 km 2.26. Scrambled control with P > 0.5: 0.010 (never-seen control realization 0.018); data with P > 0.5: 0.047. Quiet days: data P > 0.5 – vs control –. Expected count ΣP vs true array count per family-day: Spearman 0.17 (raw candidate count 0.16), ratio ΣP / truth 0.01.
- **V2008 test**: calibration slope own / ≤3 km / any = – / – / – (target 0.8–1.2); matched / ΣP own –, ≤3 km –. Scrambled control with P > 0.5: – (never-seen control realization –); data with P > 0.5: –. Quiet days: data P > 0.5 – vs control –. Expected count ΣP vs true array count per family-day: Spearman – (raw candidate count –), ratio ΣP / truth 0.00.

Caveat for reading calibration: the array catalogue is the truth but is itself incomplete, and a real LFE of a neighbouring family counts as 'real' for P but matches only at the ≤ 3-km / any levels. Matched fraction below P can therefore mean real-but-uncatalogued; matched fraction ABOVE P means P is conservative.

## Method

- Candidates: `results/mf_mt1_all` (T1 all ≥ 30-event build templates, w-weighted, 42 stations, `--drop-dead 0.05`), record tier ≥ 8.0 MAD, series data + control + control2; plus `results/mf_step1c34` = the same run with seed 20261001 → two more independent site-rolled scrambled controls (its data series is identical to mt1_all's: checked per day). K = 4 scrambled-control realizations.
- Features (all computable after the array era): `net_mad`, `net_cc`, `n_sta`, `n_sites`, `wsum`, `site_cc_max`, `site_cc_2nd`, `site_z_max`, `site_z_2nd`, `site_z_min`, `n_z2`, `n_z3`, `n_z5`, `frac_z3`, `dom_top`, `hfr_rel`, `dens_own_1h`, `dens_own_10m`, `dens_grp_1h`, `site_noise_mad_med`, `region_code`. Site values are the within-site w-weighted station cc at the candidate's aligned position (scrambled control: at the site's rolled position, exactly as the MF built it); the recomputed network cc equals the MF's to ≤ 2e-6 on every day. Site z = (site cc − median)/(1.4826 MAD) over 1000 random valid origins of the same family-day. `hfr_rel` = log10 E(10–19)/E(2–8) at the family-day's HF station minus the median of 300 random windows of that station-day. Densities count candidates of the same series.
- Density ratio: HistGradientBoosting + isotonic (3-fold by day) classifier data-vs-control on the fit set; f0/f = [(1−p)/p]·N_data/N_ctl.
- π0 rule: on quiet days (V2010 'quiet', V2008 'post'; tremor-free by the window design) every data candidate is noise → r_q = N_data(quiet) / (N_ctl(quiet)/K) per region (global if < 500 quiet data candidates), π0 = r_q·(N_ctl/K)/N_data, so P = 1 − r_q·(1−p)/(p·K), clipped to [0, 1]. Real LFEs on quiet days count as noise → conservative.
- Fit set: BUILD half of V2008 + V2010 (every series). Held out: V2010 / V2008 TEST half, A2011 (all). Leave-one-control-out refits (K = 3) give the model spread of the rate product and a P for a control realization the model never saw. Matching: within ±2 s (nearest, not one-to-one) of an array detection of the own family / any family ≤ 3 km (lat_use/lon_use) / any family, only on days the family's own array run scanned.

| region | quiet data candidates | quiet control candidates (all K) | r_q | source |
|---|---|---|---|---|
| NW | 469 | 1576 | 1.136 | global |
| SW | 1 | 13 | 1.136 | global |
| Olympic | 3042 | 10762 | 1.131 | region |
| Arm B | 20 | 82 | 1.136 | global |
| ALL | 3532 | 12433 | 1.136 | global |

Feature ranking (fit set; permutation log-loss gain on 200k candidates; univariate AUC data vs control):

| feature | perm. log-loss gain | AUC data vs control | NaN frac |
|---|---|---|---|
| n_z2 | 0.0404 | 0.586 | 0.000 |
| dens_grp_1h | 0.0124 | 0.455 | 0.000 |
| dom_top | 0.0079 | 0.410 | 0.000 |
| hfr_rel | 0.0030 | 0.474 | 0.000 |
| region_code | 0.0029 | 0.492 | 0.000 |
| n_sta | 0.0027 | 0.506 | 0.000 |
| dens_own_10m | 0.0017 | 0.502 | 0.000 |
| dens_own_1h | 0.0016 | 0.463 | 0.000 |
| site_cc_max | 0.0016 | 0.458 | 0.000 |
| site_z_min | 0.0015 | 0.566 | 0.000 |
| site_z_2nd | 0.0012 | 0.577 | 0.000 |
| site_z_max | 0.0010 | 0.438 | 0.000 |
| n_sites | 0.0010 | 0.507 | 0.000 |
| frac_z3 | 0.0008 | 0.570 | 0.000 |
| net_cc | 0.0005 | 0.458 | 0.000 |
| site_noise_mad_med | 0.0005 | 0.503 | 0.000 |
| wsum | 0.0004 | 0.444 | 0.000 |
| site_cc_2nd | 0.0004 | 0.536 | 0.000 |
| net_mad | 0.0003 | 0.516 | 0.000 |
| n_z5 | 0.0001 | 0.458 | 0.000 |
| n_z3 | 0.0000 | 0.562 | 0.000 |

## (a) Calibration — fraction matched to the array catalogue by P bin

Slope = weighted least-squares slope of matched fraction on bin-mean P (bins with ≥ 20 candidates); matched/ΣP = calibration in the large. Control rows are the chance level.

| set | region | series | n | ΣP | slope own · 3km · any | intercept own · 3km · any | matched/ΣP own · 3km · any |
|---|---|---|---|---|---|---|---|
| V2010_test | ALL | data | 8387 | 165.7 | 0.80 · 0.67 · 0.65 | 0.07 · 0.17 · 0.33 | 4.13 · 9.12 · 17.24 |
| V2010_test | ALL | control | 31948 | 173.3 | 0.01 · -0.00 · 0.06 | 0.00 · 0.02 · 0.11 | 0.29 · 3.80 · 21.19 |
| V2010_test | NW | data | 775 | 114.0 | 0.59 · 0.59 · 0.49 | 0.24 · 0.24 · 0.46 | 2.23 · 2.23 · 3.65 |
| V2010_test | NW | control | 2439 | 4.5 | – · – · – | – · – · – | 1.76 · 1.76 · 60.01 |
| V2010_test | Olympic | data | 7462 | 25.8 | – · – · – | – · – · – | 15.26 · 47.25 · 91.28 |
| V2010_test | Olympic | control | 29099 | 116.3 | 0.03 · 0.07 · 0.24 | 0.00 · 0.02 · 0.11 | 0.37 · 5.58 · 28.93 |
| V2010_test | Arm B | data | 89 | 0.0 | – · – · – | – · – · – | 230.40 · 230.40 · 990.70 |
| V2010_test | Arm B | control | 288 | 0.3 | – · – · – | – · – · – | 0.00 · 7.53 · 90.38 |
| A2011 | ALL | data | 1926 | 319.2 | 0.61 · 0.62 · 0.48 | 0.27 · 0.27 · 0.45 | 2.24 · 2.26 · 3.20 |
| A2011 | ALL | control | 5309 | 108.9 | 0.04 · 0.04 · 0.08 | 0.00 · 0.00 · 0.10 | 0.21 · 0.27 · 4.72 |
| A2011 | NW | data | 1926 | 319.2 | 0.61 · 0.62 · 0.48 | 0.27 · 0.27 · 0.45 | 2.24 · 2.26 · 3.20 |
| A2011 | NW | control | 5309 | 108.9 | 0.04 · 0.04 · 0.08 | 0.00 · 0.00 · 0.10 | 0.21 · 0.27 · 4.72 |
| fit_in_sample | ALL | data | 10175 | 2180.6 | 0.94 · 0.84 · 0.65 | 0.05 · 0.15 · 0.34 | 1.16 · 1.54 · 2.25 |
| fit_in_sample | ALL | control | 32592 | 114.2 | 0.00 · -0.02 · -0.05 | 0.00 · 0.02 · 0.13 | 0.47 · 5.80 · 36.19 |
| fit_in_sample | NW | data | 1678 | 1022.4 | 0.81 · 0.81 · 0.56 | 0.17 · 0.17 · 0.43 | 1.09 · 1.09 · 1.26 |
| fit_in_sample | NW | control | 2540 | 1.7 | – · – · – | – · – · – | 6.44 · 7.03 · 192.06 |
| fit_in_sample | Olympic | data | 7442 | 255.0 | 0.72 · 0.69 · 0.53 | 0.04 · 0.15 · 0.34 | 1.86 · 5.17 · 10.32 |
| fit_in_sample | Olympic | control | 29609 | 91.4 | 0.01 · -0.03 · -0.16 | 0.00 · 0.02 · 0.13 | 0.47 · 7.09 · 40.92 |
| fit_in_sample | Arm B | data | 937 | 838.2 | – · – · – | – · – · – | 1.03 · 1.03 · 1.07 |
| fit_in_sample | Arm B | control | 377 | 0.1 | – · – · – | – · – · – | 0.00 · 8.93 · 464.44 |

Reliability table, V2010_test, ALL regions:

| P bin | data n | mean P | frac own | frac ≤3 km | frac any | control n | control frac any |
|---|---|---|---|---|---|---|---|
| 0.0-0.1 | 8115 | 0.002 | 0.069 | 0.170 | 0.329 | 31517 | 0.115 |
| 0.1-0.2 | 79 | 0.156 | 0.089 | 0.152 | 0.380 | 267 | 0.109 |
| 0.2-0.3 | 17 | 0.252 | 0.412 | 0.471 | 0.588 | 39 | 0.128 |
| 0.3-0.4 | 24 | 0.343 | 0.375 | 0.417 | 0.583 | 31 | 0.194 |
| 0.4-0.5 | 20 | 0.468 | 0.450 | 0.450 | 0.650 | 48 | 0.250 |
| 0.5-0.6 | 25 | 0.564 | 0.480 | 0.480 | 0.720 | 24 | 0.167 |
| 0.6-0.7 | 18 | 0.672 | 0.611 | 0.611 | 0.833 | 8 | 0.125 |
| 0.7-0.8 | 56 | 0.732 | 0.625 | 0.625 | 0.839 | 22 | 0.000 |
| 0.8-0.9 | 38 | 0.862 | 0.842 | 0.842 | 0.868 | 3 | 0.333 |
| 0.9-1.0 | 30 | 0.966 | 0.833 | 0.833 | 0.933 | 2 | 0.000 |

Reliability table, A2011, ALL regions:

| P bin | data n | mean P | frac own | frac ≤3 km | frac any | control n | control frac any |
|---|---|---|---|---|---|---|---|
| 0.0-0.1 | 1360 | 0.014 | 0.271 | 0.273 | 0.452 | 5072 | 0.096 |
| 0.1-0.2 | 18 | 0.178 | 0.444 | 0.500 | 0.667 | 26 | 0.115 |
| 0.2-0.3 | 149 | 0.246 | 0.490 | 0.497 | 0.617 | 118 | 0.127 |
| 0.3-0.4 | 41 | 0.346 | 0.585 | 0.585 | 0.683 | 35 | 0.114 |
| 0.4-0.5 | 81 | 0.469 | 0.642 | 0.642 | 0.728 | 35 | 0.114 |
| 0.5-0.6 | 93 | 0.562 | 0.516 | 0.527 | 0.645 | 23 | 0.130 |
| 0.6-0.7 | 36 | 0.685 | 0.667 | 0.667 | 0.750 | 5 | 0.200 |
| 0.7-0.8 | 53 | 0.726 | 0.679 | 0.679 | 0.792 | 2 | 0.000 |
| 0.8-0.9 | 113 | 0.843 | 0.779 | 0.796 | 0.858 | 0 | – |
| 0.9-1.0 | 1 | 0.924 | 0.000 | 0.000 | 1.000 | 0 | – |

![reliability](../figures/step1_reliability.png)

## (b) Scrambled control's own P, (c) quiet vs episode days, (d) per region

`ctl P>0.5 (never seen)` = P of each control realization from the leave-that-control-out model. ΣP per realization = expected 'real' count the method would wrongly claim in one scrambled network.

| set | region | days | data n | data P>0.5 | data P>0.9 | data ΣP | ctl n / real. | ctl P>0.1 | ctl P>0.5 | ctl P>0.5 (never seen) | ctl ΣP / real. | ctl ΣP / real. (never seen) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V2010_test | ALL | all | 11536 | 0.017 | 0.003 | 225.5 | 11354 | 0.018 | 0.003 | 0.005 | 85.7 | 81.26 |
| V2010_test | ALL | quiet | 3055 | 0.001 | 0.000 | 23.4 | 3358 | 0.029 | 0.004 | 0.013 | 39.6 | 54.16 |
| V2010_test | ALL | episode | 8481 | 0.022 | 0.004 | 202.1 | 7996 | 0.014 | 0.002 | 0.002 | 46.1 | 27.10 |
| V2010_test | NW | all | 1295 | 0.123 | 0.021 | 154.4 | 1185 | 0.043 | 0.010 | 0.008 | 23.5 | 15.92 |
| V2010_test | NW | quiet | 455 | 0.007 | 0.000 | 13.5 | 556 | 0.057 | 0.007 | 0.006 | 13.4 | 8.98 |
| V2010_test | NW | episode | 840 | 0.186 | 0.032 | 140.8 | 629 | 0.030 | 0.013 | 0.009 | 10.1 | 6.94 |
| V2010_test | SW | all | 61 | 0.344 | 0.033 | 20.9 | 10 | 0.100 | 0.075 | 0.075 | 0.6 | 0.56 |
| V2010_test | SW | quiet | 2 | 0.000 | 0.000 | 0.9 | 4 | 0.062 | 0.062 | 0.062 | 0.2 | 0.16 |
| V2010_test | SW | episode | 59 | 0.356 | 0.034 | 20.0 | 6 | 0.125 | 0.083 | 0.083 | 0.4 | 0.39 |
| V2010_test | Olympic | all | 10063 | 0.000 | 0.000 | 42.3 | 10061 | 0.015 | 0.001 | 0.004 | 59.0 | 62.58 |
| V2010_test | Olympic | quiet | 2581 | 0.000 | 0.000 | 8.8 | 2777 | 0.023 | 0.003 | 0.014 | 25.7 | 44.85 |
| V2010_test | Olympic | episode | 7482 | 0.000 | 0.000 | 33.5 | 7284 | 0.012 | 0.001 | 0.000 | 33.3 | 17.73 |
| V2010_test | Arm B | all | 117 | 0.068 | 0.026 | 7.9 | 98 | 0.059 | 0.018 | 0.031 | 2.7 | 2.20 |
| V2010_test | Arm B | quiet | 17 | 0.000 | 0.000 | 0.2 | 20 | 0.025 | 0.012 | 0.012 | 0.3 | 0.16 |
| V2010_test | Arm B | episode | 100 | 0.080 | 0.030 | 7.7 | 77 | 0.068 | 0.019 | 0.036 | 2.3 | 2.04 |
| A2011 | ALL | all | 8173 | 0.047 | 0.000 | 540.0 | 7506 | 0.071 | 0.010 | 0.018 | 241.1 | 195.46 |
| A2011 | ALL | episode | 8173 | 0.047 | 0.000 | 540.0 | 7506 | 0.071 | 0.010 | 0.018 | 241.1 | 195.46 |
| A2011 | NW | all | 1945 | 0.152 | 0.001 | 323.4 | 1329 | 0.046 | 0.006 | 0.005 | 28.4 | 21.54 |
| A2011 | NW | episode | 1945 | 0.152 | 0.001 | 323.4 | 1329 | 0.046 | 0.006 | 0.005 | 28.4 | 21.54 |
| A2011 | SW | all | 10 | 0.000 | 0.000 | 0.7 | 7 | 0.034 | 0.000 | 0.000 | 0.3 | 0.69 |
| A2011 | SW | episode | 10 | 0.000 | 0.000 | 0.7 | 7 | 0.034 | 0.000 | 0.000 | 0.3 | 0.69 |
| A2011 | Olympic | all | 5966 | 0.014 | 0.000 | 197.3 | 5927 | 0.071 | 0.010 | 0.020 | 191.3 | 159.05 |
| A2011 | Olympic | episode | 5966 | 0.014 | 0.000 | 197.3 | 5927 | 0.071 | 0.010 | 0.020 | 191.3 | 159.05 |
| A2011 | Arm B | all | 252 | 0.020 | 0.000 | 18.5 | 243 | 0.214 | 0.036 | 0.026 | 21.1 | 14.18 |
| A2011 | Arm B | episode | 252 | 0.020 | 0.000 | 18.5 | 243 | 0.214 | 0.036 | 0.026 | 21.1 | 14.18 |
| fit_in_sample | ALL | all | 13798 | 0.166 | 0.126 | 2415.2 | 11285 | 0.009 | 0.001 | 0.004 | 42.5 | 86.84 |
| fit_in_sample | ALL | quiet | 3532 | 0.022 | 0.001 | 194.5 | 3108 | 0.005 | 0.001 | 0.012 | 8.7 | 43.79 |
| fit_in_sample | ALL | episode | 10266 | 0.216 | 0.169 | 2220.7 | 8177 | 0.011 | 0.001 | 0.002 | 33.8 | 43.05 |
| fit_in_sample | NW | all | 2207 | 0.511 | 0.353 | 1068.4 | 1037 | 0.015 | 0.003 | 0.004 | 7.1 | 6.59 |
| fit_in_sample | NW | quiet | 469 | 0.023 | 0.002 | 22.8 | 394 | 0.018 | 0.003 | 0.004 | 3.3 | 2.91 |
| fit_in_sample | NW | episode | 1738 | 0.642 | 0.448 | 1045.5 | 643 | 0.013 | 0.003 | 0.003 | 3.8 | 3.68 |
| fit_in_sample | SW | all | 92 | 0.446 | 0.076 | 40.1 | 24 | 0.031 | 0.021 | 0.021 | 0.5 | 0.67 |
| fit_in_sample | SW | quiet | 1 | 0.000 | 0.000 | 0.0 | 3 | 0.077 | 0.077 | 0.077 | 0.2 | 0.19 |
| fit_in_sample | SW | episode | 91 | 0.451 | 0.077 | 40.1 | 21 | 0.024 | 0.012 | 0.012 | 0.3 | 0.48 |
| fit_in_sample | Olympic | all | 10503 | 0.022 | 0.009 | 439.8 | 10103 | 0.008 | 0.000 | 0.004 | 32.1 | 77.15 |
| fit_in_sample | Olympic | quiet | 3042 | 0.022 | 0.000 | 171.7 | 2690 | 0.003 | 0.000 | 0.013 | 5.0 | 40.39 |
| fit_in_sample | Olympic | episode | 7461 | 0.023 | 0.012 | 268.2 | 7412 | 0.010 | 0.001 | 0.001 | 27.1 | 36.76 |
| fit_in_sample | Arm B | all | 996 | 0.895 | 0.859 | 866.9 | 120 | 0.052 | 0.019 | 0.021 | 2.9 | 2.44 |
| fit_in_sample | Arm B | quiet | 20 | 0.000 | 0.000 | 0.0 | 20 | 0.024 | 0.012 | 0.024 | 0.3 | 0.31 |
| fit_in_sample | Arm B | episode | 976 | 0.913 | 0.877 | 866.9 | 100 | 0.058 | 0.020 | 0.020 | 2.6 | 2.14 |
| ALL | ALL | all | 33908 | 0.085 | 0.052 | 3188.4 | 30518 | 0.028 | 0.004 | 0.008 | 371.9 | 366.04 |
| ALL | ALL | quiet | 6651 | 0.012 | 0.000 | 220.8 | 6569 | 0.017 | 0.002 | 0.012 | 49.1 | 98.83 |
| ALL | ALL | episode | 27257 | 0.102 | 0.065 | 2967.6 | 23949 | 0.031 | 0.004 | 0.007 | 322.8 | 267.21 |
| ALL | NW | all | 5506 | 0.288 | 0.147 | 1552.1 | 3583 | 0.036 | 0.007 | 0.006 | 60.0 | 44.78 |
| ALL | NW | quiet | 941 | 0.015 | 0.001 | 38.6 | 964 | 0.041 | 0.005 | 0.005 | 17.1 | 12.13 |
| ALL | NW | episode | 4565 | 0.344 | 0.177 | 1513.6 | 2618 | 0.034 | 0.007 | 0.006 | 42.9 | 32.65 |
| ALL | SW | all | 163 | 0.380 | 0.055 | 61.7 | 42 | 0.048 | 0.030 | 0.030 | 1.4 | 1.91 |
| ALL | SW | quiet | 3 | 0.000 | 0.000 | 0.9 | 7 | 0.069 | 0.069 | 0.069 | 0.4 | 0.35 |
| ALL | SW | episode | 160 | 0.388 | 0.056 | 60.8 | 34 | 0.043 | 0.022 | 0.022 | 1.0 | 1.56 |
| ALL | Olympic | all | 26871 | 0.012 | 0.003 | 680.4 | 26430 | 0.025 | 0.003 | 0.008 | 283.8 | 300.38 |
| ALL | Olympic | quiet | 5670 | 0.012 | 0.000 | 181.1 | 5556 | 0.013 | 0.002 | 0.014 | 31.0 | 85.88 |
| ALL | Olympic | episode | 21201 | 0.012 | 0.004 | 499.2 | 20874 | 0.028 | 0.003 | 0.006 | 252.7 | 214.50 |
| ALL | Arm B | all | 1368 | 0.662 | 0.628 | 894.1 | 462 | 0.139 | 0.028 | 0.026 | 26.8 | 18.97 |
| ALL | Arm B | quiet | 37 | 0.000 | 0.000 | 0.2 | 41 | 0.024 | 0.012 | 0.018 | 0.6 | 0.47 |
| ALL | Arm B | episode | 1331 | 0.680 | 0.645 | 893.9 | 422 | 0.150 | 0.030 | 0.027 | 26.2 | 18.50 |

By PNSN regional tremor-day label (label only):

| set | region | tremor day | family-days | data ΣP / fam-day | ctl ΣP / fam-day (per real.) | data P>0.5 | ctl P>0.5 |
|---|---|---|---|---|---|---|---|
| V2010_test | ALL | False | 78 | 0.300 | 0.508 | 0.001 | 0.004 |
| V2010_test | ALL | True | 198 | 1.021 | 0.233 | 0.022 | 0.002 |
| V2010_test | NW | False | 46 | 0.294 | 0.291 | 0.007 | 0.007 |
| V2010_test | NW | True | 108 | 1.304 | 0.094 | 0.186 | 0.013 |
| V2010_test | SW | False | 2 | 0.453 | 0.092 | 0.000 | 0.062 |
| V2010_test | SW | True | 24 | 0.834 | 0.016 | 0.356 | 0.083 |
| V2010_test | Olympic | False | 25 | 0.352 | 1.029 | 0.000 | 0.003 |
| V2010_test | Olympic | True | 49 | 0.684 | 0.679 | 0.000 | 0.001 |
| V2010_test | Arm B | False | 5 | 0.036 | 0.067 | 0.000 | 0.012 |
| V2010_test | Arm B | True | 17 | 0.455 | 0.137 | 0.080 | 0.019 |
| A2011 | ALL | True | 168 | 3.214 | 1.435 | 0.047 | 0.010 |
| A2011 | NW | True | 115 | 2.812 | 0.247 | 0.152 | 0.006 |
| A2011 | SW | True | 10 | 0.073 | 0.034 | 0.000 | 0.000 |
| A2011 | Olympic | True | 25 | 7.893 | 7.650 | 0.014 | 0.010 |
| A2011 | Arm B | True | 18 | 1.030 | 1.172 | 0.020 | 0.036 |

## (e) Expected-count check — ΣP per family-day vs the true array count

Episode family-days the family's array run scanned; scored half only (test half; A2011 all). truth ≤3 km = array detections of the ≤ 3-km group merged within ±2 s. `ΣP ctl` = the same sum on one scrambled-control realization.

| set | region | family-days | ΣP | ΣP ctl | truth own | truth ≤3 km | ΣP/own | ΣP/≤3 km | Spearman own | Spearman ≤3 km | Spearman raw count vs own |
|---|---|---|---|---|---|---|---|---|---|---|---|
| V2010_test | ALL | 4572 | 182.1 | 45.7 | 30921 | 147454 | 0.006 | 0.001 | 0.11 | 0.12 | 0.10 |
| V2010_test | NW | 2640 | 140.8 | 10.1 | 17082 | 37188 | 0.008 | 0.004 | 0.09 | 0.09 | 0.05 |
| V2010_test | Olympic | 784 | 33.5 | 33.3 | 7056 | 85199 | 0.005 | 0.000 | 0.09 | 0.17 | 0.11 |
| V2010_test | Arm B | 1148 | 7.7 | 2.3 | 6783 | 25067 | 0.001 | 0.000 | 0.11 | 0.11 | 0.14 |
| V2008_test | ALL | 4934 | 0.0 | 0.0 | 23962 | 44153 | 0.000 | 0.000 | – | – | – |
| V2008_test | NW | 3036 | 0.0 | 0.0 | 16129 | 25449 | 0.000 | 0.000 | – | – | – |
| V2008_test | SW | 1898 | 0.0 | 0.0 | 7833 | 18704 | 0.000 | 0.000 | – | – | – |
| A2011 | ALL | 3828 | 323.4 | 28.2 | 33306 | 50539 | 0.010 | 0.006 | 0.17 | 0.16 | 0.16 |
| A2011 | NW | 3828 | 323.4 | 28.2 | 33306 | 50539 | 0.010 | 0.006 | 0.17 | 0.16 | 0.16 |

## Rate product — expected real repeats per family and day

`results/step1_prob/rate_family_day.parquet`: one row per observable family-day (all 97 days, V2008 + V2010 + A2011): `sumP` (whole day; build-half P is in-sample), `sumP_scored_half`, `sd_total` = √(Σ P(1−P) + model spread²) with the model spread = SD of ΣP over the 4 leave-one-control-out refits, `sumP_ctl_per_real` (the same sum on one scrambled network), truth counts where the family's array run scanned the day. Tremor-day label is a column only.

| region | days | family-days | mean ΣP / fam-day | mean sd | mean ΣP ctl / fam-day | family-days with ΣP ≥ 5 |
|---|---|---|---|---|---|---|
| Olympic | episode | 2352 | 0.21 | 3.17 | 0.11 | 42 |
| Olympic | post | 84 | 0.00 | nan | 0.00 | 0 |
| Olympic | quiet | 280 | 0.65 | 4.24 | 0.11 | 14 |
| Arm B | episode | 3444 | 0.26 | 0.57 | 0.01 | 44 |
| Arm B | post | 123 | 0.00 | nan | 0.00 | 0 |
| Arm B | quiet | 410 | 0.00 | 0.07 | 0.00 | 0 |
| NW | episode | 11088 | 0.14 | 0.72 | 0.00 | 71 |
| NW | post | 396 | 0.00 | nan | 0.00 | 0 |
| NW | quiet | 1320 | 0.03 | 0.31 | 0.01 | 2 |
| SW | episode | 6132 | 0.01 | 0.39 | 0.00 | 5 |
| SW | post | 219 | 0.00 | nan | 0.00 | 0 |
| SW | quiet | 730 | 0.00 | 0.36 | 0.00 | 0 |

Files: `candidates.parquet` (every candidate, every feature, P, P_loco_c1..c4, matches), `calibration.csv`, `calibration_slopes.csv`, `control_and_quiet.csv`, `tremor_label_split.csv`, `expected_count.csv`, `feature_importance.csv`, `model_info.json`, `rate_family_day.parquet`.
