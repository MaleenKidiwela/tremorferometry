# STEP 1 — a probability for every repeat (local false-discovery rate, scrambled control as the noise model)

Written 2026-09-30 13:39Z by `scripts/step1_report.py` (numbers only; interpretation after review). Queue `scripts/step1_prob_queue.sh`; outputs `results/step1_prob/`; figure `figures/step1_reliability.png`.

## Summary (plain language)

Every candidate of the all-templates matched filter (≥ 8.0 MAD, nothing cut) gets a probability P(real) that it is a real repeat and not noise. The noise model is the scrambled control: the same matched filter run after shifting each site to a different time of the same day, so nothing lines up across sites except by chance. A classifier learns how data candidates differ from scrambled-control candidates (per-site correlations above each site's own noise level, HF ratio, how many candidates the family and its ≤ 3-km group have within ±1 h, …); the control's share, scaled by how many candidates the data has on tremor-free days, is the local fraction of noise. Fitted on the build half of V2008 + V2010; checked on the V2010 / V2008 test half and on A2011 (a separate episode).

- **V2010 test**: calibration slope own / ≤3 km / any = 0.63 / 0.56 / 0.53 (target 0.8–1.2); matched / ΣP own 1.27, ≤3 km 2.70. Scrambled control with P > 0.5: 0.003 (never-seen control realization 0.005); data with P > 0.5: 0.019. Quiet days: data P > 0.5 0.011 vs control 0.005. Expected count ΣP vs true array count per family-day: Spearman 0.50 (raw candidate count 0.46), ratio ΣP / truth 0.20.
- **A2011**: calibration slope own / ≤3 km / any = 0.54 / 0.55 / 0.49 (target 0.8–1.2); matched / ΣP own 1.06, ≤3 km 1.07. Scrambled control with P > 0.5: 0.020 (never-seen control realization 0.044); data with P > 0.5: 0.033. Quiet days: data P > 0.5 – vs control –. Expected count ΣP vs true array count per family-day: Spearman 0.47 (raw candidate count 0.40), ratio ΣP / truth 0.25.
- **V2008 test**: calibration slope own / ≤3 km / any = 0.55 / 0.56 / 0.53 (target 0.8–1.2); matched / ΣP own 1.49, ≤3 km 1.51. Scrambled control with P > 0.5: 0.007 (never-seen control realization 0.005); data with P > 0.5: 0.031. Quiet days: data P > 0.5 0.001 vs control 0.006. Expected count ΣP vs true array count per family-day: Spearman 0.46 (raw candidate count 0.43), ratio ΣP / truth 0.11.

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
| n_z2 | 0.0392 | 0.567 | 0.000 |
| dom_top | 0.0297 | 0.431 | 0.000 |
| net_mad | 0.0074 | 0.514 | 0.000 |
| dens_grp_1h | 0.0051 | 0.469 | 0.000 |
| region_code | 0.0037 | 0.464 | 0.000 |
| dens_own_10m | 0.0033 | 0.501 | 0.000 |
| net_cc | 0.0015 | 0.476 | 0.000 |
| n_sta | 0.0015 | 0.500 | 0.000 |
| dens_own_1h | 0.0015 | 0.479 | 0.000 |
| site_cc_max | 0.0011 | 0.474 | 0.000 |
| hfr_rel | 0.0011 | 0.482 | 0.000 |
| n_sites | 0.0009 | 0.505 | 0.000 |
| wsum | 0.0009 | 0.463 | 0.000 |
| n_sites_pnet03 | 0.0007 | 0.523 | 0.000 |
| pnet_2nd | 0.0007 | 0.520 | 0.000 |
| site_noise_mad_med | 0.0006 | 0.514 | 0.000 |
| site_z_max | 0.0005 | 0.455 | 0.000 |
| site_z_min | 0.0005 | 0.544 | 0.000 |
| frac_z3 | 0.0005 | 0.550 | 0.000 |
| site_z_2nd | 0.0004 | 0.561 | 0.000 |
| pnet_max | 0.0003 | 0.516 | 0.000 |
| site_cc_2nd | 0.0002 | 0.533 | 0.000 |
| n_z3 | 0.0001 | 0.548 | 0.000 |
| n_sites_pnet05 | 0.0000 | 0.518 | 0.000 |
| n_z5 | 0.0000 | 0.472 | 0.000 |

## (a) Calibration — fraction matched to the array catalogue by P bin

Slope = weighted least-squares slope of matched fraction on bin-mean P (bins with ≥ 20 candidates); matched/ΣP = calibration in the large. Control rows are the chance level.

| set | region | series | n | ΣP | slope own · 3km · any | intercept own · 3km · any | matched/ΣP own · 3km · any |
|---|---|---|---|---|---|---|---|
| V2010_test | ALL | data | 115301 | 6093.4 | 0.63 · 0.56 · 0.53 | 0.03 · 0.11 · 0.24 | 1.27 · 2.70 · 5.03 |
| V2010_test | ALL | control | 447017 | 14366.9 | 0.01 · 0.05 · -0.00 | 0.00 · 0.01 · 0.08 | 0.06 · 0.49 · 2.34 |
| V2010_test | NW | data | 12527 | 2598.9 | 0.55 · 0.55 · 0.46 | 0.18 · 0.18 · 0.38 | 1.41 · 1.42 · 2.27 |
| V2010_test | NW | control | 38294 | 1404.7 | 0.01 · 0.01 · 0.07 | 0.01 · 0.01 · 0.11 | 0.18 · 0.19 · 2.94 |
| V2010_test | Olympic | data | 101337 | 3399.8 | 0.18 · 0.18 · 0.05 | 0.03 · 0.12 · 0.24 | 1.18 · 3.72 · 7.12 |
| V2010_test | Olympic | control | 402973 | 12810.8 | 0.00 · 0.07 · -0.06 | 0.00 · 0.01 · 0.07 | 0.04 · 0.53 · 2.27 |
| V2010_test | Arm B | data | 1385 | 59.0 | 0.61 · 0.62 · 0.37 | 0.03 · 0.03 · 0.33 | 1.29 · 1.39 · 8.19 |
| V2010_test | Arm B | control | 5681 | 110.2 | 0.01 · 0.01 · -0.03 | 0.00 · 0.01 · 0.08 | 0.23 · 0.42 · 3.90 |
| A2011 | ALL | data | 61999 | 8359.8 | 0.54 · 0.55 · 0.49 | 0.07 · 0.07 · 0.19 | 1.06 · 1.07 · 1.93 |
| A2011 | ALL | control | 222561 | 15782.6 | 0.00 · 0.00 · -0.01 | 0.00 · 0.00 · 0.05 | 0.04 · 0.05 · 0.68 |
| A2011 | NW | data | 61999 | 8359.8 | 0.54 · 0.55 · 0.49 | 0.07 · 0.07 · 0.19 | 1.06 · 1.07 · 1.93 |
| A2011 | NW | control | 222561 | 15782.6 | 0.00 · 0.00 · -0.01 | 0.00 · 0.00 · 0.05 | 0.04 · 0.05 · 0.68 |
| V2008_test | ALL | data | 13943 | 2584.7 | 0.55 · 0.56 · 0.53 | 0.17 · 0.17 · 0.31 | 1.49 · 1.51 · 2.17 |
| V2008_test | ALL | control | 44129 | 1329.5 | 0.01 · 0.01 · 0.06 | 0.00 · 0.00 · 0.08 | 0.15 · 0.16 · 2.70 |
| V2008_test | NW | data | 13647 | 2417.6 | 0.55 · 0.55 · 0.53 | 0.17 · 0.18 · 0.31 | 1.53 · 1.55 · 2.25 |
| V2008_test | NW | control | 43591 | 1253.1 | 0.01 · 0.01 · 0.07 | 0.00 · 0.00 · 0.08 | 0.16 · 0.17 · 2.83 |
| V2008_test | SW | data | 217 | 124.6 | 0.74 · 0.75 · 0.69 | 0.05 · 0.06 · 0.15 | 0.83 · 0.86 · 0.95 |
| V2008_test | SW | control | 479 | 33.7 | -0.03 · -0.03 · 0.17 | 0.01 · 0.01 · 0.08 | 0.06 · 0.06 · 1.31 |
| fit_in_sample | ALL | data | 160191 | 40446.2 | 0.96 · 0.88 · 0.76 | 0.00 · 0.08 · 0.22 | 0.98 · 1.21 · 1.62 |
| fit_in_sample | ALL | control | 501947 | 13527.7 | 0.01 · 0.01 · -0.00 | 0.00 · 0.02 · 0.08 | 0.08 · 0.61 · 3.13 |
| fit_in_sample | NW | data | 43208 | 23902.1 | 0.85 · 0.85 · 0.69 | 0.11 · 0.11 · 0.29 | 1.04 · 1.04 · 1.21 |
| fit_in_sample | NW | control | 84142 | 2470.7 | -0.00 · 0.00 · 0.04 | 0.00 · 0.01 · 0.09 | 0.16 · 0.18 · 3.26 |
| fit_in_sample | SW | data | 6810 | 6699.4 | 0.88 · 0.88 · 0.78 | 0.12 · 0.12 · 0.21 | 1.00 · 1.00 · 1.00 |
| fit_in_sample | SW | control | 465 | 35.3 | 0.00 · -0.02 · -0.25 | 0.00 · 0.00 · 0.10 | 0.00 · 0.03 · 1.11 |
| fit_in_sample | Olympic | data | 105984 | 7024.8 | 0.93 · 0.86 · 0.70 | -0.01 · 0.08 · 0.21 | 0.74 · 2.07 · 3.83 |
| fit_in_sample | Olympic | control | 412112 | 10857.8 | 0.00 · 0.04 · -0.08 | 0.00 · 0.02 · 0.08 | 0.06 · 0.71 · 3.10 |
| fit_in_sample | Arm B | data | 4092 | 2765.5 | 0.96 · 0.95 · 0.61 | 0.04 · 0.04 · 0.39 | 1.02 · 1.02 · 1.18 |
| fit_in_sample | Arm B | control | 5095 | 73.7 | 0.13 · 0.16 · 0.32 | 0.02 · 0.02 · 0.12 | 1.17 · 1.83 · 8.33 |

Reliability table, V2010_test, ALL regions:

| P bin | data n | mean P | frac own | frac ≤3 km | frac any | control n | control frac any |
|---|---|---|---|---|---|---|---|
| 0.0-0.1 | 101703 | 0.018 | 0.050 | 0.127 | 0.253 | 408611 | 0.076 |
| 0.1-0.2 | 7296 | 0.137 | 0.073 | 0.148 | 0.247 | 27093 | 0.066 |
| 0.2-0.3 | 3149 | 0.237 | 0.103 | 0.173 | 0.268 | 9471 | 0.067 |
| 0.3-0.4 | 296 | 0.351 | 0.172 | 0.223 | 0.318 | 473 | 0.072 |
| 0.4-0.5 | 249 | 0.439 | 0.273 | 0.317 | 0.482 | 327 | 0.098 |
| 0.5-0.6 | 191 | 0.554 | 0.408 | 0.440 | 0.565 | 239 | 0.126 |
| 0.6-0.7 | 163 | 0.651 | 0.460 | 0.491 | 0.669 | 206 | 0.126 |
| 0.7-0.8 | 191 | 0.748 | 0.539 | 0.550 | 0.691 | 166 | 0.193 |
| 0.8-0.9 | 278 | 0.851 | 0.536 | 0.554 | 0.683 | 148 | 0.142 |
| 0.9-1.0 | 1785 | 0.979 | 0.719 | 0.725 | 0.826 | 283 | 0.099 |

Reliability table, A2011, ALL regions:

| P bin | data n | mean P | frac own | frac ≤3 km | frac any | control n | control frac any |
|---|---|---|---|---|---|---|---|
| 0.0-0.1 | 41438 | 0.025 | 0.102 | 0.103 | 0.230 | 168461 | 0.050 |
| 0.1-0.2 | 9322 | 0.140 | 0.097 | 0.098 | 0.196 | 32631 | 0.040 |
| 0.2-0.3 | 4632 | 0.242 | 0.114 | 0.114 | 0.205 | 14219 | 0.041 |
| 0.3-0.4 | 765 | 0.352 | 0.173 | 0.174 | 0.256 | 2126 | 0.041 |
| 0.4-0.5 | 753 | 0.444 | 0.258 | 0.263 | 0.349 | 1808 | 0.054 |
| 0.5-0.6 | 678 | 0.551 | 0.351 | 0.354 | 0.476 | 1040 | 0.062 |
| 0.6-0.7 | 620 | 0.651 | 0.437 | 0.439 | 0.523 | 855 | 0.073 |
| 0.7-0.8 | 474 | 0.749 | 0.447 | 0.451 | 0.525 | 577 | 0.062 |
| 0.8-0.9 | 537 | 0.854 | 0.529 | 0.534 | 0.646 | 449 | 0.053 |
| 0.9-1.0 | 2780 | 0.972 | 0.671 | 0.683 | 0.757 | 395 | 0.071 |

Reliability table, V2008_test, ALL regions:

| P bin | data n | mean P | frac own | frac ≤3 km | frac any | control n | control frac any |
|---|---|---|---|---|---|---|---|
| 0.0-0.1 | 10267 | 0.014 | 0.170 | 0.171 | 0.299 | 41384 | 0.080 |
| 0.1-0.2 | 731 | 0.138 | 0.342 | 0.343 | 0.487 | 1449 | 0.092 |
| 0.2-0.3 | 382 | 0.244 | 0.469 | 0.474 | 0.628 | 521 | 0.109 |
| 0.3-0.4 | 85 | 0.349 | 0.541 | 0.553 | 0.647 | 85 | 0.129 |
| 0.4-0.5 | 109 | 0.446 | 0.560 | 0.560 | 0.642 | 121 | 0.099 |
| 0.5-0.6 | 136 | 0.559 | 0.471 | 0.471 | 0.654 | 125 | 0.128 |
| 0.6-0.7 | 134 | 0.651 | 0.545 | 0.545 | 0.709 | 100 | 0.140 |
| 0.7-0.8 | 142 | 0.755 | 0.542 | 0.542 | 0.683 | 69 | 0.072 |
| 0.8-0.9 | 188 | 0.855 | 0.580 | 0.580 | 0.676 | 83 | 0.072 |
| 0.9-1.0 | 1769 | 0.980 | 0.705 | 0.717 | 0.804 | 192 | 0.141 |

![reliability](../figures/step1_reliability.png)

## (b) Scrambled control's own P, (c) quiet vs episode days, (d) per region

`ctl P>0.5 (never seen)` = P of each control realization from the leave-that-control-out model. ΣP per realization = expected 'real' count the method would wrongly claim in one scrambled network.

| set | region | days | data n | data P>0.5 | data P>0.9 | data ΣP | ctl n / real. | ctl P>0.1 | ctl P>0.5 | ctl P>0.5 (never seen) | ctl ΣP / real. | ctl ΣP / real. (never seen) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V2010_test | ALL | all | 172854 | 0.019 | 0.011 | 8711.9 | 168969 | 0.080 | 0.003 | 0.005 | 5703.9 | 6403.27 |
| V2010_test | ALL | quiet | 48576 | 0.011 | 0.000 | 2157.1 | 48451 | 0.069 | 0.005 | 0.009 | 1847.4 | 2004.02 |
| V2010_test | ALL | episode | 124278 | 0.023 | 0.015 | 6554.8 | 120518 | 0.084 | 0.002 | 0.003 | 3856.5 | 4399.25 |
| V2010_test | NW | all | 21800 | 0.113 | 0.081 | 2838.9 | 19002 | 0.060 | 0.012 | 0.015 | 559.9 | 676.29 |
| V2010_test | NW | quiet | 5429 | 0.005 | 0.001 | 119.8 | 5521 | 0.034 | 0.004 | 0.006 | 118.8 | 148.16 |
| V2010_test | NW | episode | 16371 | 0.149 | 0.107 | 2719.1 | 13482 | 0.071 | 0.015 | 0.019 | 441.0 | 528.13 |
| V2010_test | SW | all | 381 | 0.417 | 0.289 | 164.1 | 168 | 0.191 | 0.030 | 0.042 | 11.9 | 15.27 |
| V2010_test | SW | quiet | 33 | 0.061 | 0.061 | 3.3 | 44 | 0.158 | 0.023 | 0.028 | 2.9 | 3.50 |
| V2010_test | SW | episode | 348 | 0.451 | 0.310 | 160.8 | 124 | 0.202 | 0.032 | 0.047 | 9.1 | 11.77 |
| V2010_test | Olympic | all | 148866 | 0.004 | 0.000 | 5621.3 | 147917 | 0.082 | 0.002 | 0.003 | 5085.5 | 5642.29 |
| V2010_test | Olympic | quiet | 42747 | 0.012 | 0.000 | 2027.0 | 42469 | 0.074 | 0.005 | 0.009 | 1717.0 | 1839.23 |
| V2010_test | Olympic | episode | 106119 | 0.002 | 0.000 | 3594.3 | 105448 | 0.086 | 0.001 | 0.001 | 3368.5 | 3803.06 |
| V2010_test | Arm B | all | 1807 | 0.033 | 0.017 | 87.5 | 1882 | 0.058 | 0.011 | 0.017 | 46.6 | 69.42 |
| V2010_test | Arm B | quiet | 367 | 0.008 | 0.000 | 7.0 | 417 | 0.059 | 0.007 | 0.010 | 8.8 | 13.13 |
| V2010_test | Arm B | episode | 1440 | 0.039 | 0.022 | 80.5 | 1465 | 0.058 | 0.012 | 0.019 | 37.9 | 56.30 |
| A2011 | ALL | all | 294824 | 0.033 | 0.010 | 28725.2 | 287372 | 0.295 | 0.020 | 0.044 | 23912.6 | 31138.05 |
| A2011 | ALL | episode | 294824 | 0.033 | 0.010 | 28725.2 | 287372 | 0.295 | 0.020 | 0.044 | 23912.6 | 31138.05 |
| A2011 | NW | all | 61999 | 0.082 | 0.045 | 8359.8 | 55640 | 0.243 | 0.015 | 0.024 | 3945.7 | 4746.19 |
| A2011 | NW | episode | 61999 | 0.082 | 0.045 | 8359.8 | 55640 | 0.243 | 0.015 | 0.024 | 3945.7 | 4746.19 |
| A2011 | SW | all | 709 | 0.203 | 0.034 | 253.5 | 700 | 0.988 | 0.210 | 0.569 | 261.0 | 373.82 |
| A2011 | SW | episode | 709 | 0.203 | 0.034 | 253.5 | 700 | 0.988 | 0.210 | 0.569 | 261.0 | 373.82 |
| A2011 | Olympic | all | 215645 | 0.020 | 0.000 | 18846.5 | 214685 | 0.310 | 0.021 | 0.049 | 18555.2 | 24527.22 |
| A2011 | Olympic | episode | 215645 | 0.020 | 0.000 | 18846.5 | 214685 | 0.310 | 0.021 | 0.049 | 18555.2 | 24527.22 |
| A2011 | Arm B | all | 16471 | 0.020 | 0.003 | 1265.4 | 16348 | 0.247 | 0.013 | 0.018 | 1150.7 | 1490.82 |
| A2011 | Arm B | episode | 16471 | 0.020 | 0.003 | 1265.4 | 16348 | 0.247 | 0.013 | 0.018 | 1150.7 | 1490.82 |
| V2008_test | ALL | all | 93754 | 0.031 | 0.019 | 5172.4 | 90830 | 0.072 | 0.007 | 0.005 | 3326.7 | 3653.06 |
| V2008_test | ALL | quiet | 13913 | 0.001 | 0.000 | 340.0 | 13122 | 0.111 | 0.006 | 0.003 | 495.1 | 582.84 |
| V2008_test | ALL | episode | 79841 | 0.036 | 0.022 | 4832.4 | 77708 | 0.065 | 0.008 | 0.006 | 2831.6 | 3070.22 |
| V2008_test | NW | all | 16322 | 0.136 | 0.101 | 2484.9 | 13546 | 0.053 | 0.010 | 0.013 | 373.8 | 450.58 |
| V2008_test | NW | quiet | 1993 | 0.002 | 0.001 | 45.7 | 1966 | 0.033 | 0.004 | 0.007 | 47.9 | 60.58 |
| V2008_test | NW | episode | 14329 | 0.154 | 0.115 | 2439.2 | 11580 | 0.056 | 0.011 | 0.014 | 325.9 | 390.00 |
| V2008_test | SW | all | 305 | 0.548 | 0.397 | 168.3 | 145 | 0.365 | 0.083 | 0.107 | 20.3 | 24.25 |
| V2008_test | SW | quiet | 9 | 0.111 | 0.000 | 1.2 | 11 | 0.349 | 0.070 | 0.070 | 1.2 | 1.52 |
| V2008_test | SW | episode | 296 | 0.561 | 0.409 | 167.1 | 134 | 0.366 | 0.084 | 0.110 | 19.1 | 22.73 |
| V2008_test | Olympic | all | 73842 | 0.006 | 0.000 | 2460.7 | 73896 | 0.077 | 0.007 | 0.004 | 2884.6 | 3108.37 |
| V2008_test | Olympic | quiet | 11474 | 0.000 | 0.000 | 286.0 | 10724 | 0.129 | 0.007 | 0.002 | 440.6 | 512.76 |
| V2008_test | Olympic | episode | 62368 | 0.008 | 0.000 | 2174.7 | 63171 | 0.068 | 0.007 | 0.004 | 2443.9 | 2595.61 |
| V2008_test | Arm B | all | 3285 | 0.003 | 0.001 | 58.5 | 3242 | 0.024 | 0.001 | 0.002 | 48.0 | 69.87 |
| V2008_test | Arm B | quiet | 437 | 0.002 | 0.000 | 7.2 | 420 | 0.012 | 0.002 | 0.003 | 5.4 | 7.99 |
| V2008_test | Arm B | episode | 2848 | 0.003 | 0.001 | 51.3 | 2822 | 0.026 | 0.001 | 0.002 | 42.7 | 61.87 |
| fit_in_sample | ALL | all | 297796 | 0.122 | 0.111 | 47081.1 | 257759 | 0.048 | 0.002 | 0.003 | 7088.3 | 8413.92 |
| fit_in_sample | ALL | quiet | 61798 | 0.007 | 0.001 | 3032.6 | 60574 | 0.046 | 0.001 | 0.002 | 1712.3 | 2004.42 |
| fit_in_sample | ALL | episode | 235998 | 0.152 | 0.140 | 44048.5 | 197185 | 0.048 | 0.002 | 0.004 | 5375.9 | 6409.49 |
| fit_in_sample | NW | all | 55588 | 0.427 | 0.392 | 24328.3 | 32886 | 0.055 | 0.011 | 0.015 | 898.6 | 1147.58 |
| fit_in_sample | NW | quiet | 7633 | 0.007 | 0.002 | 243.7 | 7305 | 0.035 | 0.005 | 0.007 | 161.9 | 215.38 |
| fit_in_sample | NW | episode | 47955 | 0.494 | 0.454 | 24084.6 | 25581 | 0.061 | 0.012 | 0.017 | 736.7 | 932.20 |
| fit_in_sample | SW | all | 7337 | 0.942 | 0.927 | 6920.2 | 349 | 0.287 | 0.038 | 0.067 | 33.1 | 44.29 |
| fit_in_sample | SW | quiet | 63 | 0.032 | 0.000 | 3.9 | 64 | 0.245 | 0.031 | 0.058 | 5.1 | 7.14 |
| fit_in_sample | SW | episode | 7274 | 0.950 | 0.935 | 6916.3 | 285 | 0.296 | 0.039 | 0.069 | 28.0 | 37.15 |
| fit_in_sample | Olympic | all | 226905 | 0.013 | 0.008 | 12950.7 | 219660 | 0.047 | 0.000 | 0.001 | 6079.6 | 7096.98 |
| fit_in_sample | Olympic | quiet | 53232 | 0.006 | 0.000 | 2759.2 | 52394 | 0.048 | 0.001 | 0.002 | 1533.6 | 1762.38 |
| fit_in_sample | Olympic | episode | 173673 | 0.015 | 0.010 | 10191.5 | 167266 | 0.046 | 0.000 | 0.001 | 4546.0 | 5334.61 |
| fit_in_sample | Arm B | all | 7966 | 0.350 | 0.336 | 2881.9 | 4863 | 0.028 | 0.005 | 0.009 | 76.9 | 125.06 |
| fit_in_sample | Arm B | quiet | 870 | 0.013 | 0.003 | 25.9 | 811 | 0.033 | 0.003 | 0.008 | 11.8 | 19.52 |
| fit_in_sample | Arm B | episode | 7096 | 0.392 | 0.377 | 2856.1 | 4052 | 0.027 | 0.005 | 0.010 | 65.1 | 105.54 |
| ALL | ALL | all | 867908 | 0.061 | 0.046 | 90174.5 | 814130 | 0.145 | 0.009 | 0.018 | 40329.3 | 49953.46 |
| ALL | ALL | quiet | 126292 | 0.008 | 0.000 | 5615.7 | 124542 | 0.062 | 0.003 | 0.005 | 4131.4 | 4673.49 |
| ALL | ALL | episode | 741616 | 0.070 | 0.054 | 84558.8 | 689588 | 0.159 | 0.010 | 0.021 | 36197.9 | 45279.97 |
| ALL | NW | all | 156925 | 0.215 | 0.179 | 38175.2 | 122218 | 0.141 | 0.013 | 0.019 | 5811.3 | 7062.09 |
| ALL | NW | quiet | 15254 | 0.006 | 0.001 | 414.3 | 15052 | 0.035 | 0.004 | 0.007 | 334.8 | 431.64 |
| ALL | NW | episode | 141671 | 0.237 | 0.198 | 37761.0 | 107166 | 0.156 | 0.014 | 0.020 | 5476.5 | 6630.45 |
| ALL | SW | all | 8754 | 0.845 | 0.807 | 7517.6 | 1373 | 0.640 | 0.129 | 0.324 | 327.1 | 458.77 |
| ALL | SW | quiet | 106 | 0.047 | 0.019 | 8.4 | 122 | 0.220 | 0.031 | 0.049 | 9.3 | 12.34 |
| ALL | SW | episode | 8648 | 0.854 | 0.817 | 7509.2 | 1251 | 0.681 | 0.139 | 0.351 | 317.9 | 446.43 |
| ALL | Olympic | all | 672541 | 0.012 | 0.003 | 40184.6 | 664021 | 0.143 | 0.008 | 0.018 | 32865.1 | 40672.65 |
| ALL | Olympic | quiet | 109234 | 0.008 | 0.000 | 5152.6 | 107682 | 0.066 | 0.003 | 0.005 | 3760.4 | 4187.69 |
| ALL | Olympic | episode | 563307 | 0.013 | 0.003 | 35032.0 | 556340 | 0.158 | 0.009 | 0.020 | 29104.7 | 36484.96 |
| ALL | Arm B | all | 29688 | 0.107 | 0.093 | 4297.1 | 26518 | 0.165 | 0.010 | 0.015 | 1325.7 | 1759.95 |
| ALL | Arm B | quiet | 1698 | 0.009 | 0.002 | 40.5 | 1687 | 0.035 | 0.004 | 0.007 | 27.0 | 41.81 |
| ALL | Arm B | episode | 27990 | 0.113 | 0.099 | 4256.6 | 24831 | 0.173 | 0.010 | 0.015 | 1298.8 | 1718.13 |

By PNSN regional tremor-day label (label only):

| set | region | tremor day | family-days | data ΣP / fam-day | ctl ΣP / fam-day (per real.) | data P>0.5 | ctl P>0.5 |
|---|---|---|---|---|---|---|---|
| V2010_test | ALL | False | 884 | 2.449 | 2.099 | 0.011 | 0.005 |
| V2010_test | ALL | True | 2690 | 2.434 | 1.431 | 0.023 | 0.002 |
| V2010_test | NW | False | 551 | 0.231 | 0.228 | 0.005 | 0.004 |
| V2010_test | NW | True | 1596 | 1.699 | 0.272 | 0.153 | 0.016 |
| V2010_test | SW | False | 41 | 0.089 | 0.090 | 0.047 | 0.028 |
| V2010_test | SW | True | 204 | 0.787 | 0.040 | 0.464 | 0.031 |
| V2010_test | Olympic | False | 244 | 8.307 | 7.037 | 0.012 | 0.005 |
| V2010_test | Olympic | True | 720 | 4.992 | 4.678 | 0.002 | 0.001 |
| V2010_test | Arm B | False | 48 | 0.145 | 0.182 | 0.008 | 0.007 |
| V2010_test | Arm B | True | 170 | 0.474 | 0.223 | 0.039 | 0.012 |
| A2011 | ALL | False | 172 | 1.690 | 1.703 | 0.056 | 0.057 |
| A2011 | ALL | True | 3962 | 7.177 | 5.962 | 0.033 | 0.020 |
| A2011 | NW | False | 68 | 3.438 | 3.459 | 0.035 | 0.037 |
| A2011 | NW | True | 2662 | 3.053 | 1.394 | 0.083 | 0.014 |
| A2011 | SW | False | 104 | 0.546 | 0.555 | 0.198 | 0.202 |
| A2011 | SW | True | 213 | 0.924 | 0.955 | 0.205 | 0.213 |
| A2011 | Olympic | True | 724 | 26.031 | 25.629 | 0.020 | 0.021 |
| A2011 | Arm B | True | 363 | 3.486 | 3.170 | 0.020 | 0.013 |
| V2008_test | ALL | False | 2788 | 1.855 | 1.193 | 0.031 | 0.007 |
| V2008_test | NW | False | 1667 | 1.491 | 0.224 | 0.136 | 0.010 |
| V2008_test | SW | False | 186 | 0.905 | 0.109 | 0.548 | 0.083 |
| V2008_test | Olympic | False | 712 | 3.456 | 4.051 | 0.006 | 0.007 |
| V2008_test | Arm B | False | 223 | 0.262 | 0.215 | 0.003 | 0.001 |

## (e) Expected-count check — ΣP per family-day vs the true array count

Episode family-days the family's array run scanned; scored half only (test half; A2011 all). truth ≤3 km = array detections of the ≤ 3-km group merged within ±2 s. `ΣP ctl` = the same sum on one scrambled-control realization.

| set | region | family-days | ΣP | ΣP ctl | truth own | truth ≤3 km | ΣP/own | ΣP/≤3 km | Spearman own | Spearman ≤3 km | Spearman raw count vs own |
|---|---|---|---|---|---|---|---|---|---|---|---|
| V2010_test | ALL | 4572 | 6093.4 | 3587.0 | 30921 | 147454 | 0.197 | 0.041 | 0.50 | 0.55 | 0.46 |
| V2010_test | NW | 2640 | 2598.9 | 347.0 | 17082 | 37188 | 0.152 | 0.070 | 0.52 | 0.47 | 0.48 |
| V2010_test | Olympic | 784 | 3414.0 | 3202.5 | 7056 | 85199 | 0.484 | 0.040 | 0.38 | 0.33 | 0.28 |
| V2010_test | Arm B | 1148 | 80.5 | 37.4 | 6783 | 25067 | 0.012 | 0.003 | 0.18 | 0.08 | 0.19 |
| V2008_test | ALL | 4934 | 2584.7 | 318.8 | 23962 | 44153 | 0.108 | 0.059 | 0.46 | 0.41 | 0.43 |
| V2008_test | NW | 3036 | 2417.6 | 306.5 | 16129 | 25449 | 0.150 | 0.095 | 0.53 | 0.50 | 0.49 |
| V2008_test | SW | 1898 | 167.1 | 12.3 | 7833 | 18704 | 0.021 | 0.009 | 0.28 | 0.25 | 0.25 |
| A2011 | ALL | 3828 | 8359.8 | 3862.6 | 33306 | 50539 | 0.251 | 0.165 | 0.47 | 0.45 | 0.40 |
| A2011 | NW | 3828 | 8359.8 | 3862.6 | 33306 | 50539 | 0.251 | 0.165 | 0.47 | 0.45 | 0.40 |

## Rate product — expected real repeats per family and day

`results/step1_prob/rate_family_day.parquet`: one row per observable family-day (all 97 days, V2008 + V2010 + A2011): `sumP` (whole day; build-half P is in-sample), `sumP_scored_half`, `sd_total` = √(Σ P(1−P) + model spread²) with the model spread = SD of ΣP over the 4 leave-one-control-out refits, `sumP_ctl_per_real` (the same sum on one scrambled network), truth counts where the family's array run scanned the day. Tremor-day label is a column only.

| region | days | family-days | mean ΣP / fam-day | mean sd | mean ΣP ctl / fam-day | family-days with ΣP ≥ 5 |
|---|---|---|---|---|---|---|
| Olympic | episode | 2352 | 14.89 | 6.85 | 12.37 | 1551 |
| Olympic | post | 84 | 10.24 | 4.31 | 8.42 | 57 |
| Olympic | quiet | 280 | 15.33 | 5.10 | 10.90 | 221 |
| Arm B | episode | 3444 | 1.24 | 0.64 | 0.37 | 218 |
| Arm B | post | 123 | 0.09 | 0.46 | 0.07 | 0 |
| Arm B | quiet | 410 | 0.07 | 0.38 | 0.04 | 0 |
| NW | episode | 11088 | 3.41 | 0.66 | 0.48 | 1453 |
| NW | post | 396 | 0.29 | 0.37 | 0.23 | 5 |
| NW | quiet | 1320 | 0.23 | 0.33 | 0.18 | 14 |
| SW | episode | 6132 | 1.22 | 0.27 | 0.04 | 327 |
| SW | post | 219 | 0.02 | 0.21 | 0.00 | 0 |
| SW | quiet | 730 | 0.01 | 0.16 | 0.00 | 0 |

Files: `candidates.parquet` (every candidate, every feature, P, P_loco_c1..c4, matches), `calibration.csv`, `calibration_slopes.csv`, `control_and_quiet.csv`, `tremor_label_split.csv`, `expected_count.csv`, `feature_importance.csv`, `model_info.json`, `rate_family_day.parquet`.
