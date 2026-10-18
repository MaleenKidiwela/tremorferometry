# STEP 1 — a probability for every repeat (local false-discovery rate, scrambled control as the noise model)

Written 2026-09-30 13:45Z by `scripts/step1_report.py` (numbers only; interpretation after review). Queue `scripts/step1_prob_queue.sh`; outputs `results/step1_prob/`; figure `figures/step1_reliability.png`.

## Summary (plain language)

Every candidate of the all-templates matched filter (≥ 8.0 MAD, nothing cut) gets a probability P(real) that it is a real repeat and not noise. The noise model is the scrambled control: the same matched filter run after shifting each site to a different time of the same day, so nothing lines up across sites except by chance. A classifier learns how data candidates differ from scrambled-control candidates (per-site correlations above each site's own noise level, HF ratio, how many candidates the family and its ≤ 3-km group have within ±1 h, …); the control's share, scaled by how many candidates the data has on tremor-free days, is the local fraction of noise. Fitted on the build half of V2008 + V2010; checked on the V2010 / V2008 test half and on A2011 (a separate episode).

- **V2010 test**: calibration slope own / ≤3 km / any = 0.62 / 0.55 / 0.52 (target 0.8–1.2); matched / ΣP own 1.08, ≤3 km 2.29. Scrambled control with P > 0.5: 0.005 (never-seen control realization 0.005); data with P > 0.5: 0.026. Quiet days: data P > 0.5 0.024 vs control 0.012. Expected count ΣP vs true array count per family-day: Spearman 0.50 (raw candidate count 0.46), ratio ΣP / truth 0.23.
- **A2011**: calibration slope own / ≤3 km / any = 0.55 / 0.56 / 0.52 (target 0.8–1.2); matched / ΣP own 0.89, ≤3 km 0.90. Scrambled control with P > 0.5: 0.008 (never-seen control realization 0.021); data with P > 0.5: 0.026. Quiet days: data P > 0.5 – vs control –. Expected count ΣP vs true array count per family-day: Spearman 0.48 (raw candidate count 0.40), ratio ΣP / truth 0.30.
- **V2008 test**: calibration slope own / ≤3 km / any = 0.58 / 0.59 / 0.57 (target 0.8–1.2); matched / ΣP own 1.33, ≤3 km 1.34. Scrambled control with P > 0.5: 0.007 (never-seen control realization 0.006); data with P > 0.5: 0.034. Quiet days: data P > 0.5 0.001 vs control 0.001. Expected count ΣP vs true array count per family-day: Spearman 0.47 (raw candidate count 0.43), ratio ΣP / truth 0.12.

Caveat for reading calibration: the array catalogue is the truth but is itself incomplete, and a real LFE of a neighbouring family counts as 'real' for P but matches only at the ≤ 3-km / any levels. Matched fraction below P can therefore mean real-but-uncatalogued; matched fraction ABOVE P means P is conservative.

## Method

- Candidates: `results/mf_mt1_all` (T1 all ≥ 30-event build templates, w-weighted, 42 stations, `--drop-dead 0.05`), record tier ≥ 8.0 MAD, series data + control + control2; plus `results/mf_step1c34` = the same run with seed 20261001 → two more independent site-rolled scrambled controls (its data series is identical to mt1_all's: checked per day). K = 4 scrambled-control realizations.
- Features (all computable after the array era): `net_mad`, `net_cc`, `n_sta`, `n_sites`, `wsum`, `site_cc_max`, `site_cc_2nd`, `site_z_max`, `site_z_2nd`, `site_z_min`, `n_z2`, `n_z3`, `n_z5`, `frac_z3`, `dom_top`, `hfr_rel`, `dens_own_1h`, `dens_own_10m`, `dens_grp_1h`, `site_noise_mad_med`, `region_code`, `pnet_max`, `pnet_2nd`, `n_sites_pnet03`, `n_sites_pnet05`. Site values are the within-site w-weighted station cc at the candidate's aligned position (scrambled control: at the site's rolled position, exactly as the MF built it); the recomputed network cc equals the MF's to ≤ 2e-6 on every day. Site z = (site cc − median)/(1.4826 MAD) over 1000 random valid origins of the same family-day. `hfr_rel` = log10 E(10–19)/E(2–8) at the family-day's HF station minus the median of 300 random windows of that station-day. Densities count candidates of the same series.
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
| n_z2 | 0.0337 | 0.567 | 0.000 |
| dom_top | 0.0318 | 0.431 | 0.000 |
| net_mad | 0.0054 | 0.514 | 0.000 |
| dens_grp_1h | 0.0053 | 0.469 | 0.000 |
| dens_own_10m | 0.0027 | 0.501 | 0.000 |
| region_code | 0.0025 | 0.464 | 0.000 |
| net_cc | 0.0019 | 0.476 | 0.000 |
| dens_own_1h | 0.0016 | 0.479 | 0.000 |
| n_sta | 0.0015 | 0.500 | 0.000 |
| pnet_max | 0.0013 | 0.509 | 0.000 |
| pnet_2nd | 0.0012 | 0.517 | 0.000 |
| hfr_rel | 0.0012 | 0.482 | 0.000 |
| n_sites | 0.0009 | 0.505 | 0.000 |
| site_cc_max | 0.0009 | 0.474 | 0.000 |
| n_sites_pnet03 | 0.0009 | 0.520 | 0.000 |
| wsum | 0.0008 | 0.463 | 0.000 |
| site_z_max | 0.0006 | 0.455 | 0.000 |
| site_noise_mad_med | 0.0005 | 0.514 | 0.000 |
| site_z_min | 0.0005 | 0.544 | 0.000 |
| frac_z3 | 0.0005 | 0.550 | 0.000 |
| site_z_2nd | 0.0004 | 0.561 | 0.000 |
| site_cc_2nd | 0.0002 | 0.533 | 0.000 |
| n_z3 | 0.0001 | 0.548 | 0.000 |
| n_sites_pnet05 | 0.0001 | 0.514 | 0.000 |
| n_z5 | 0.0000 | 0.472 | 0.000 |

## (a) Calibration — fraction matched to the array catalogue by P bin

Slope = weighted least-squares slope of matched fraction on bin-mean P (bins with ≥ 20 candidates); matched/ΣP = calibration in the large. Control rows are the chance level.

| set | region | series | n | ΣP | slope own · 3km · any | intercept own · 3km · any | matched/ΣP own · 3km · any |
|---|---|---|---|---|---|---|---|
| V2010_test | ALL | data | 115301 | 7179.2 | 0.62 · 0.55 · 0.52 | 0.03 · 0.11 · 0.23 | 1.08 · 2.29 · 4.27 |
| V2010_test | ALL | control | 447017 | 16843.1 | 0.01 · 0.04 · -0.01 | 0.00 · 0.01 · 0.08 | 0.05 · 0.42 · 1.99 |
| V2010_test | NW | data | 12527 | 3069.4 | 0.57 · 0.57 · 0.49 | 0.15 · 0.15 · 0.35 | 1.19 · 1.20 · 1.92 |
| V2010_test | NW | control | 38294 | 1571.1 | 0.01 · 0.02 · 0.06 | 0.01 · 0.01 · 0.11 | 0.16 · 0.17 · 2.63 |
| V2010_test | Olympic | data | 101356 | 4005.0 | 0.15 · 0.13 · 0.02 | 0.03 · 0.12 · 0.24 | 1.00 · 3.16 · 6.05 |
| V2010_test | Olympic | control | 402973 | 15082.2 | 0.00 · 0.06 · -0.05 | 0.00 · 0.01 · 0.07 | 0.04 · 0.45 · 1.92 |
| V2010_test | Arm B | data | 1366 | 73.5 | 0.65 · 0.66 · 0.44 | 0.02 · 0.02 · 0.32 | 0.97 · 1.03 · 6.46 |
| V2010_test | Arm B | control | 5665 | 129.9 | 0.04 · 0.04 · 0.04 | 0.00 · 0.01 · 0.07 | 0.19 · 0.35 · 3.29 |
| A2011 | ALL | data | 61999 | 9918.0 | 0.55 · 0.56 · 0.52 | 0.06 · 0.06 · 0.18 | 0.89 · 0.90 · 1.63 |
| A2011 | ALL | control | 222561 | 18464.7 | -0.00 · 0.00 · -0.00 | 0.00 · 0.00 · 0.05 | 0.04 · 0.04 · 0.58 |
| A2011 | NW | data | 61999 | 9918.0 | 0.55 · 0.56 · 0.52 | 0.06 · 0.06 · 0.18 | 0.89 · 0.90 · 1.63 |
| A2011 | NW | control | 222561 | 18464.7 | -0.00 · 0.00 · -0.00 | 0.00 · 0.00 · 0.05 | 0.04 · 0.04 · 0.58 |
| V2008_test | ALL | data | 13943 | 2898.4 | 0.58 · 0.59 · 0.57 | 0.16 · 0.16 · 0.29 | 1.33 · 1.34 · 1.94 |
| V2008_test | ALL | control | 44129 | 1489.6 | 0.01 · 0.01 · 0.08 | 0.00 · 0.00 · 0.08 | 0.14 · 0.15 · 2.41 |
| V2008_test | NW | data | 13647 | 2732.9 | 0.58 · 0.59 · 0.57 | 0.16 · 0.16 · 0.29 | 1.36 · 1.37 · 1.99 |
| V2008_test | NW | control | 43591 | 1412.4 | 0.01 · 0.01 · 0.09 | 0.00 · 0.00 · 0.08 | 0.14 · 0.15 · 2.51 |
| V2008_test | SW | data | 225 | 121.6 | 0.74 · 0.75 · 0.70 | 0.08 · 0.08 · 0.17 | 0.88 · 0.90 · 1.02 |
| V2008_test | SW | control | 476 | 34.9 | -0.03 · -0.03 · 0.04 | 0.01 · 0.01 · 0.09 | 0.06 · 0.06 · 1.23 |
| fit_in_sample | ALL | data | 160191 | 41653.2 | 0.96 · 0.88 · 0.76 | -0.00 · 0.08 · 0.21 | 0.95 · 1.18 · 1.58 |
| fit_in_sample | ALL | control | 501947 | 15353.7 | 0.00 · 0.01 · -0.01 | 0.00 · 0.02 · 0.08 | 0.07 · 0.54 · 2.76 |
| fit_in_sample | NW | data | 43208 | 24382.0 | 0.85 · 0.85 · 0.70 | 0.09 · 0.10 · 0.27 | 1.02 · 1.02 · 1.18 |
| fit_in_sample | NW | control | 84142 | 2803.7 | -0.00 · 0.00 · 0.06 | 0.00 · 0.01 · 0.09 | 0.14 · 0.16 · 2.87 |
| fit_in_sample | SW | data | 6790 | 6687.5 | 0.90 · 0.90 · 0.79 | 0.10 · 0.10 · 0.20 | 1.00 · 1.00 · 1.00 |
| fit_in_sample | SW | control | 467 | 39.6 | 0.00 · -0.02 · -0.11 | 0.00 · 0.00 · 0.09 | 0.00 · 0.03 · 0.99 |
| fit_in_sample | Olympic | data | 105984 | 7742.8 | 0.89 · 0.81 · 0.65 | -0.02 · 0.08 · 0.21 | 0.67 · 1.88 · 3.47 |
| fit_in_sample | Olympic | control | 412095 | 12300.5 | 0.00 · 0.03 · -0.09 | 0.00 · 0.02 · 0.08 | 0.05 · 0.62 · 2.73 |
| fit_in_sample | Arm B | data | 4130 | 2798.3 | 0.97 · 0.96 · 0.62 | 0.03 · 0.04 · 0.38 | 1.01 · 1.02 · 1.18 |
| fit_in_sample | Arm B | control | 5118 | 119.9 | 0.09 · 0.11 · 0.17 | 0.02 · 0.02 · 0.12 | 0.73 · 1.14 · 5.15 |

Reliability table, V2010_test, ALL regions:

| P bin | data n | mean P | frac own | frac ≤3 km | frac any | control n | control frac any |
|---|---|---|---|---|---|---|---|
| 0.0-0.1 | 95846 | 0.016 | 0.048 | 0.126 | 0.252 | 387390 | 0.076 |
| 0.1-0.2 | 11606 | 0.132 | 0.060 | 0.134 | 0.236 | 44100 | 0.067 |
| 0.2-0.3 | 3735 | 0.247 | 0.087 | 0.150 | 0.258 | 12392 | 0.070 |
| 0.3-0.4 | 793 | 0.346 | 0.139 | 0.202 | 0.299 | 1649 | 0.072 |
| 0.4-0.5 | 201 | 0.445 | 0.294 | 0.333 | 0.463 | 345 | 0.084 |
| 0.5-0.6 | 222 | 0.553 | 0.378 | 0.423 | 0.563 | 287 | 0.115 |
| 0.6-0.7 | 250 | 0.653 | 0.488 | 0.504 | 0.656 | 263 | 0.091 |
| 0.7-0.8 | 241 | 0.752 | 0.469 | 0.498 | 0.676 | 196 | 0.112 |
| 0.8-0.9 | 303 | 0.853 | 0.492 | 0.508 | 0.670 | 160 | 0.163 |
| 0.9-1.0 | 2104 | 0.982 | 0.728 | 0.735 | 0.839 | 235 | 0.106 |

Reliability table, A2011, ALL regions:

| P bin | data n | mean P | frac own | frac ≤3 km | frac any | control n | control frac any |
|---|---|---|---|---|---|---|---|
| 0.0-0.1 | 35064 | 0.022 | 0.096 | 0.097 | 0.221 | 148603 | 0.050 |
| 0.1-0.2 | 11865 | 0.140 | 0.088 | 0.089 | 0.195 | 42220 | 0.043 |
| 0.2-0.3 | 6811 | 0.240 | 0.105 | 0.106 | 0.210 | 21842 | 0.041 |
| 0.3-0.4 | 2139 | 0.345 | 0.154 | 0.157 | 0.251 | 6293 | 0.045 |
| 0.4-0.5 | 474 | 0.446 | 0.264 | 0.264 | 0.380 | 1008 | 0.064 |
| 0.5-0.6 | 519 | 0.554 | 0.389 | 0.391 | 0.524 | 753 | 0.058 |
| 0.6-0.7 | 571 | 0.652 | 0.424 | 0.429 | 0.581 | 666 | 0.080 |
| 0.7-0.8 | 549 | 0.751 | 0.448 | 0.450 | 0.585 | 478 | 0.086 |
| 0.8-0.9 | 608 | 0.851 | 0.510 | 0.515 | 0.627 | 396 | 0.098 |
| 0.9-1.0 | 3399 | 0.979 | 0.669 | 0.680 | 0.762 | 302 | 0.073 |

Reliability table, V2008_test, ALL regions:

| P bin | data n | mean P | frac own | frac ≤3 km | frac any | control n | control frac any |
|---|---|---|---|---|---|---|---|
| 0.0-0.1 | 9664 | 0.015 | 0.152 | 0.152 | 0.275 | 40510 | 0.080 |
| 0.1-0.2 | 950 | 0.140 | 0.306 | 0.312 | 0.476 | 2244 | 0.084 |
| 0.2-0.3 | 417 | 0.252 | 0.439 | 0.444 | 0.621 | 570 | 0.111 |
| 0.3-0.4 | 160 | 0.351 | 0.550 | 0.550 | 0.669 | 151 | 0.132 |
| 0.4-0.5 | 115 | 0.450 | 0.591 | 0.600 | 0.722 | 85 | 0.024 |
| 0.5-0.6 | 142 | 0.554 | 0.437 | 0.437 | 0.627 | 131 | 0.153 |
| 0.6-0.7 | 185 | 0.652 | 0.535 | 0.535 | 0.654 | 101 | 0.149 |
| 0.7-0.8 | 176 | 0.750 | 0.545 | 0.545 | 0.665 | 85 | 0.118 |
| 0.8-0.9 | 224 | 0.854 | 0.536 | 0.540 | 0.701 | 99 | 0.111 |
| 0.9-1.0 | 1910 | 0.985 | 0.724 | 0.735 | 0.826 | 153 | 0.163 |

![reliability](../figures/step1_reliability.png)

## (b) Scrambled control's own P, (c) quiet vs episode days, (d) per region

`ctl P>0.5 (never seen)` = P of each control realization from the leave-that-control-out model. ΣP per realization = expected 'real' count the method would wrongly claim in one scrambled network.

| set | region | days | data n | data P>0.5 | data P>0.9 | data ΣP | ctl n / real. | ctl P>0.1 | ctl P>0.5 | ctl P>0.5 (never seen) | ctl ΣP / real. | ctl ΣP / real. (never seen) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V2010_test | ALL | all | 172854 | 0.026 | 0.013 | 10366.1 | 168969 | 0.134 | 0.005 | 0.005 | 6806.9 | 7005.46 |
| V2010_test | ALL | quiet | 48576 | 0.024 | 0.000 | 2672.5 | 48451 | 0.139 | 0.012 | 0.010 | 2289.3 | 2235.38 |
| V2010_test | ALL | episode | 124278 | 0.027 | 0.018 | 7693.7 | 120518 | 0.131 | 0.003 | 0.004 | 4517.6 | 4770.09 |
| V2010_test | NW | all | 21800 | 0.135 | 0.094 | 3341.4 | 19002 | 0.078 | 0.012 | 0.014 | 634.9 | 696.00 |
| V2010_test | NW | quiet | 5429 | 0.004 | 0.001 | 146.2 | 5521 | 0.045 | 0.002 | 0.004 | 140.5 | 144.53 |
| V2010_test | NW | episode | 16371 | 0.179 | 0.125 | 3195.2 | 13482 | 0.091 | 0.015 | 0.019 | 494.4 | 551.47 |
| V2010_test | SW | all | 381 | 0.457 | 0.307 | 176.2 | 168 | 0.235 | 0.031 | 0.045 | 12.1 | 15.54 |
| V2010_test | SW | quiet | 33 | 0.091 | 0.030 | 3.3 | 44 | 0.181 | 0.023 | 0.034 | 2.5 | 3.04 |
| V2010_test | SW | episode | 348 | 0.491 | 0.333 | 172.9 | 124 | 0.255 | 0.034 | 0.049 | 9.6 | 12.50 |
| V2010_test | Olympic | all | 148866 | 0.009 | 0.000 | 6735.4 | 147917 | 0.141 | 0.004 | 0.004 | 6102.0 | 6217.18 |
| V2010_test | Olympic | quiet | 42747 | 0.026 | 0.000 | 2514.7 | 42469 | 0.152 | 0.013 | 0.010 | 2135.9 | 2074.64 |
| V2010_test | Olympic | episode | 106119 | 0.002 | 0.000 | 4220.7 | 105448 | 0.137 | 0.001 | 0.001 | 3966.1 | 4142.54 |
| V2010_test | Arm B | all | 1807 | 0.043 | 0.021 | 113.2 | 1882 | 0.084 | 0.011 | 0.017 | 57.8 | 76.73 |
| V2010_test | Arm B | quiet | 367 | 0.008 | 0.000 | 8.3 | 417 | 0.070 | 0.007 | 0.011 | 10.3 | 13.17 |
| V2010_test | Arm B | episode | 1440 | 0.051 | 0.026 | 104.9 | 1465 | 0.088 | 0.012 | 0.019 | 47.5 | 63.57 |
| A2011 | ALL | all | 294824 | 0.026 | 0.012 | 26879.4 | 287372 | 0.304 | 0.008 | 0.021 | 21251.5 | 25941.55 |
| A2011 | ALL | episode | 294824 | 0.026 | 0.012 | 26879.4 | 287372 | 0.304 | 0.008 | 0.021 | 21251.5 | 25941.55 |
| A2011 | NW | all | 61999 | 0.091 | 0.055 | 9918.0 | 55640 | 0.332 | 0.012 | 0.016 | 4616.2 | 4602.67 |
| A2011 | NW | episode | 61999 | 0.091 | 0.055 | 9918.0 | 55640 | 0.332 | 0.012 | 0.016 | 4616.2 | 4602.67 |
| A2011 | SW | all | 709 | 0.142 | 0.016 | 247.4 | 700 | 0.990 | 0.165 | 0.324 | 253.4 | 306.16 |
| A2011 | SW | episode | 709 | 0.142 | 0.016 | 247.4 | 700 | 0.990 | 0.165 | 0.324 | 253.4 | 306.16 |
| A2011 | Olympic | all | 215645 | 0.007 | 0.000 | 15432.2 | 214685 | 0.298 | 0.007 | 0.023 | 15244.0 | 19706.90 |
| A2011 | Olympic | episode | 215645 | 0.007 | 0.000 | 15432.2 | 214685 | 0.298 | 0.007 | 0.023 | 15244.0 | 19706.90 |
| A2011 | Arm B | all | 16471 | 0.016 | 0.006 | 1281.9 | 16348 | 0.253 | 0.006 | 0.011 | 1137.9 | 1325.82 |
| A2011 | Arm B | episode | 16471 | 0.016 | 0.006 | 1281.9 | 16348 | 0.253 | 0.006 | 0.011 | 1137.9 | 1325.82 |
| V2008_test | ALL | all | 93754 | 0.034 | 0.021 | 5889.4 | 90830 | 0.115 | 0.007 | 0.006 | 3768.5 | 3815.40 |
| V2008_test | ALL | quiet | 13913 | 0.001 | 0.000 | 346.9 | 13122 | 0.135 | 0.001 | 0.001 | 502.3 | 567.14 |
| V2008_test | ALL | episode | 79841 | 0.040 | 0.025 | 5542.5 | 77708 | 0.111 | 0.008 | 0.007 | 3266.2 | 3248.26 |
| V2008_test | NW | all | 16322 | 0.153 | 0.111 | 2818.5 | 13546 | 0.069 | 0.010 | 0.012 | 426.9 | 460.38 |
| V2008_test | NW | quiet | 1993 | 0.002 | 0.001 | 55.0 | 1966 | 0.035 | 0.003 | 0.004 | 57.8 | 57.76 |
| V2008_test | NW | episode | 14329 | 0.174 | 0.126 | 2763.5 | 11580 | 0.075 | 0.011 | 0.014 | 369.1 | 402.62 |
| V2008_test | SW | all | 305 | 0.518 | 0.364 | 166.7 | 145 | 0.418 | 0.086 | 0.120 | 20.5 | 24.21 |
| V2008_test | SW | quiet | 9 | 0.111 | 0.000 | 1.1 | 11 | 0.349 | 0.093 | 0.093 | 1.2 | 1.34 |
| V2008_test | SW | episode | 296 | 0.530 | 0.375 | 165.5 | 134 | 0.424 | 0.086 | 0.123 | 19.3 | 22.88 |
| V2008_test | Olympic | all | 73842 | 0.007 | 0.001 | 2822.2 | 73896 | 0.126 | 0.007 | 0.005 | 3254.1 | 3254.29 |
| V2008_test | Olympic | quiet | 11474 | 0.000 | 0.000 | 280.7 | 10724 | 0.158 | 0.000 | 0.000 | 435.9 | 500.05 |
| V2008_test | Olympic | episode | 62368 | 0.008 | 0.001 | 2541.5 | 63171 | 0.121 | 0.008 | 0.006 | 2818.2 | 2754.24 |
| V2008_test | Arm B | all | 3285 | 0.005 | 0.001 | 82.0 | 3242 | 0.038 | 0.001 | 0.003 | 67.1 | 76.52 |
| V2008_test | Arm B | quiet | 437 | 0.002 | 0.000 | 10.1 | 420 | 0.019 | 0.002 | 0.004 | 7.4 | 8.00 |
| V2008_test | Arm B | episode | 2848 | 0.005 | 0.001 | 71.9 | 2822 | 0.041 | 0.001 | 0.003 | 59.6 | 68.53 |
| fit_in_sample | ALL | all | 297796 | 0.124 | 0.111 | 49407.0 | 257759 | 0.095 | 0.002 | 0.003 | 8249.9 | 9014.22 |
| fit_in_sample | ALL | quiet | 61798 | 0.008 | 0.000 | 3550.5 | 60574 | 0.105 | 0.002 | 0.002 | 2068.0 | 2204.89 |
| fit_in_sample | ALL | episode | 235998 | 0.155 | 0.140 | 45856.5 | 197185 | 0.092 | 0.002 | 0.004 | 6181.8 | 6809.33 |
| fit_in_sample | NW | all | 55588 | 0.434 | 0.394 | 24868.8 | 32886 | 0.073 | 0.011 | 0.014 | 1023.6 | 1172.78 |
| fit_in_sample | NW | quiet | 7633 | 0.004 | 0.001 | 280.7 | 7305 | 0.046 | 0.003 | 0.005 | 192.8 | 206.63 |
| fit_in_sample | NW | episode | 47955 | 0.503 | 0.456 | 24588.1 | 25581 | 0.080 | 0.013 | 0.017 | 830.8 | 966.14 |
| fit_in_sample | SW | all | 7337 | 0.944 | 0.927 | 6934.2 | 349 | 0.366 | 0.039 | 0.065 | 34.9 | 46.06 |
| fit_in_sample | SW | quiet | 63 | 0.016 | 0.000 | 3.3 | 64 | 0.253 | 0.023 | 0.047 | 4.6 | 6.63 |
| fit_in_sample | SW | episode | 7274 | 0.952 | 0.935 | 6930.9 | 285 | 0.392 | 0.042 | 0.069 | 30.3 | 39.43 |
| fit_in_sample | Olympic | all | 226905 | 0.014 | 0.008 | 14676.7 | 219660 | 0.099 | 0.001 | 0.002 | 7090.9 | 7654.72 |
| fit_in_sample | Olympic | quiet | 53232 | 0.008 | 0.000 | 3234.7 | 52394 | 0.114 | 0.001 | 0.002 | 1855.9 | 1971.12 |
| fit_in_sample | Olympic | episode | 173673 | 0.016 | 0.010 | 11442.1 | 167266 | 0.094 | 0.000 | 0.001 | 5234.9 | 5683.60 |
| fit_in_sample | Arm B | all | 7966 | 0.352 | 0.335 | 2927.3 | 4863 | 0.046 | 0.005 | 0.010 | 100.5 | 140.66 |
| fit_in_sample | Arm B | quiet | 870 | 0.013 | 0.001 | 31.8 | 811 | 0.040 | 0.003 | 0.008 | 14.7 | 20.50 |
| fit_in_sample | Arm B | episode | 7096 | 0.394 | 0.376 | 2895.5 | 4052 | 0.047 | 0.006 | 0.010 | 85.8 | 120.16 |
| ALL | ALL | all | 867908 | 0.060 | 0.047 | 93102.8 | 814130 | 0.179 | 0.006 | 0.011 | 40423.6 | 46146.35 |
| ALL | ALL | quiet | 126292 | 0.013 | 0.000 | 6669.1 | 124542 | 0.122 | 0.005 | 0.005 | 4949.0 | 5095.41 |
| ALL | ALL | episode | 741616 | 0.069 | 0.055 | 86433.7 | 689588 | 0.189 | 0.006 | 0.012 | 35474.6 | 41050.94 |
| ALL | NW | all | 156925 | 0.226 | 0.187 | 41142.1 | 122218 | 0.191 | 0.011 | 0.015 | 6739.1 | 6972.65 |
| ALL | NW | quiet | 15254 | 0.004 | 0.001 | 488.0 | 15052 | 0.044 | 0.003 | 0.004 | 398.3 | 416.00 |
| ALL | NW | episode | 141671 | 0.249 | 0.207 | 40654.1 | 107166 | 0.212 | 0.012 | 0.016 | 6340.8 | 6556.65 |
| ALL | SW | all | 8754 | 0.842 | 0.805 | 7536.6 | 1373 | 0.673 | 0.107 | 0.201 | 321.9 | 393.21 |
| ALL | SW | quiet | 106 | 0.047 | 0.009 | 7.8 | 122 | 0.235 | 0.029 | 0.047 | 8.4 | 11.21 |
| ALL | SW | episode | 8648 | 0.852 | 0.815 | 7528.8 | 1251 | 0.716 | 0.115 | 0.216 | 313.5 | 381.99 |
| ALL | Olympic | all | 672541 | 0.010 | 0.003 | 40014.5 | 664021 | 0.176 | 0.004 | 0.009 | 31995.1 | 37155.54 |
| ALL | Olympic | quiet | 109234 | 0.015 | 0.000 | 6122.4 | 107682 | 0.134 | 0.006 | 0.005 | 4508.8 | 4625.30 |
| ALL | Olympic | episode | 563307 | 0.009 | 0.003 | 33892.1 | 556340 | 0.184 | 0.004 | 0.010 | 27486.3 | 32530.24 |
| ALL | Arm B | all | 29688 | 0.106 | 0.095 | 4409.7 | 26518 | 0.175 | 0.006 | 0.010 | 1367.5 | 1624.95 |
| ALL | Arm B | quiet | 1698 | 0.009 | 0.001 | 50.9 | 1687 | 0.043 | 0.004 | 0.008 | 33.5 | 42.90 |
| ALL | Arm B | episode | 27990 | 0.112 | 0.100 | 4358.8 | 24831 | 0.184 | 0.006 | 0.010 | 1334.1 | 1582.05 |

By PNSN regional tremor-day label (label only):

| set | region | tremor day | family-days | data ΣP / fam-day | ctl ΣP / fam-day (per real.) | data P>0.5 | ctl P>0.5 |
|---|---|---|---|---|---|---|---|
| V2010_test | ALL | False | 884 | 3.033 | 2.600 | 0.024 | 0.012 |
| V2010_test | ALL | True | 2690 | 2.857 | 1.676 | 0.027 | 0.003 |
| V2010_test | NW | False | 551 | 0.281 | 0.269 | 0.004 | 0.002 |
| V2010_test | NW | True | 1596 | 1.997 | 0.305 | 0.183 | 0.016 |
| V2010_test | SW | False | 41 | 0.087 | 0.081 | 0.070 | 0.028 |
| V2010_test | SW | True | 204 | 0.846 | 0.043 | 0.506 | 0.033 |
| V2010_test | Olympic | False | 244 | 10.306 | 8.754 | 0.026 | 0.013 |
| V2010_test | Olympic | True | 720 | 5.862 | 5.508 | 0.002 | 0.001 |
| V2010_test | Arm B | False | 48 | 0.174 | 0.215 | 0.008 | 0.007 |
| V2010_test | Arm B | True | 170 | 0.617 | 0.280 | 0.051 | 0.012 |
| A2011 | ALL | False | 172 | 1.668 | 1.675 | 0.025 | 0.028 |
| A2011 | ALL | True | 3962 | 6.712 | 5.291 | 0.026 | 0.008 |
| A2011 | NW | False | 68 | 3.402 | 3.405 | 0.008 | 0.006 |
| A2011 | NW | True | 2662 | 3.639 | 1.647 | 0.093 | 0.012 |
| A2011 | SW | False | 104 | 0.535 | 0.545 | 0.136 | 0.181 |
| A2011 | SW | True | 213 | 0.900 | 0.924 | 0.145 | 0.160 |
| A2011 | Olympic | True | 724 | 21.315 | 21.055 | 0.007 | 0.007 |
| A2011 | Arm B | True | 363 | 3.531 | 3.135 | 0.016 | 0.006 |
| V2008_test | ALL | False | 2788 | 2.112 | 1.352 | 0.034 | 0.007 |
| V2008_test | NW | False | 1667 | 1.691 | 0.256 | 0.153 | 0.010 |
| V2008_test | SW | False | 186 | 0.896 | 0.110 | 0.518 | 0.086 |
| V2008_test | Olympic | False | 712 | 3.964 | 4.570 | 0.007 | 0.007 |
| V2008_test | Arm B | False | 223 | 0.368 | 0.301 | 0.005 | 0.001 |

## (e) Expected-count check — ΣP per family-day vs the true array count

Episode family-days the family's array run scanned; scored half only (test half; A2011 all). truth ≤3 km = array detections of the ≤ 3-km group merged within ±2 s. `ΣP ctl` = the same sum on one scrambled-control realization.

| set | region | family-days | ΣP | ΣP ctl | truth own | truth ≤3 km | ΣP/own | ΣP/≤3 km | Spearman own | Spearman ≤3 km | Spearman raw count vs own |
|---|---|---|---|---|---|---|---|---|---|---|---|
| V2010_test | ALL | 4572 | 7179.2 | 4205.2 | 30921 | 147454 | 0.232 | 0.049 | 0.50 | 0.55 | 0.46 |
| V2010_test | NW | 2640 | 3069.4 | 388.0 | 17082 | 37188 | 0.180 | 0.083 | 0.53 | 0.47 | 0.48 |
| V2010_test | Olympic | 784 | 4005.0 | 3770.4 | 7056 | 85199 | 0.568 | 0.047 | 0.37 | 0.32 | 0.28 |
| V2010_test | Arm B | 1148 | 104.8 | 46.9 | 6783 | 25067 | 0.015 | 0.004 | 0.21 | 0.10 | 0.19 |
| V2008_test | ALL | 4934 | 2898.4 | 358.8 | 23962 | 44153 | 0.121 | 0.066 | 0.47 | 0.42 | 0.43 |
| V2008_test | NW | 3036 | 2732.9 | 346.4 | 16129 | 25449 | 0.169 | 0.107 | 0.54 | 0.51 | 0.49 |
| V2008_test | SW | 1898 | 165.5 | 12.4 | 7833 | 18704 | 0.021 | 0.009 | 0.28 | 0.24 | 0.25 |
| A2011 | ALL | 3828 | 9918.0 | 4521.4 | 33306 | 50539 | 0.298 | 0.196 | 0.48 | 0.46 | 0.40 |
| A2011 | NW | 3828 | 9918.0 | 4521.4 | 33306 | 50539 | 0.298 | 0.196 | 0.48 | 0.46 | 0.40 |

## Rate product — expected real repeats per family and day

`results/step1_prob/rate_family_day.parquet`: one row per observable family-day (all 97 days, V2008 + V2010 + A2011): `sumP` (whole day; build-half P is in-sample), `sumP_scored_half`, `sd_total` = √(Σ P(1−P) + model spread²) with the model spread = SD of ΣP over the 4 leave-one-control-out refits, `sumP_ctl_per_real` (the same sum on one scrambled network), truth counts where the family's array run scanned the day. Tremor-day label is a column only.

| region | days | family-days | mean ΣP / fam-day | mean sd | mean ΣP ctl / fam-day | family-days with ΣP ≥ 5 |
|---|---|---|---|---|---|---|
| Olympic | episode | 2352 | 14.41 | 6.05 | 11.68 | 1602 |
| Olympic | post | 84 | 10.24 | 4.40 | 8.61 | 65 |
| Olympic | quiet | 280 | 18.79 | 5.32 | 13.52 | 228 |
| Arm B | episode | 3444 | 1.27 | 0.66 | 0.38 | 233 |
| Arm B | post | 123 | 0.13 | 0.52 | 0.09 | 0 |
| Arm B | quiet | 410 | 0.08 | 0.42 | 0.05 | 1 |
| NW | episode | 11088 | 3.67 | 0.69 | 0.56 | 1550 |
| NW | post | 396 | 0.34 | 0.35 | 0.29 | 5 |
| NW | quiet | 1320 | 0.27 | 0.31 | 0.21 | 17 |
| SW | episode | 6132 | 1.23 | 0.27 | 0.04 | 330 |
| SW | post | 219 | 0.02 | 0.21 | 0.00 | 0 |
| SW | quiet | 730 | 0.01 | 0.12 | 0.00 | 0 |

Files: `candidates.parquet` (every candidate, every feature, P, P_loco_c1..c4, matches), `calibration.csv`, `calibration_slopes.csv`, `control_and_quiet.csv`, `tremor_label_split.csv`, `expected_count.csv`, `feature_importance.csv`, `model_info.json`, `rate_family_day.parquet`.
