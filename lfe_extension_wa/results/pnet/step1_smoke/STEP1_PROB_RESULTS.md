# STEP 1 — a probability for every repeat (local false-discovery rate, scrambled control as the noise model)

Written 2026-09-30 13:04Z by `scripts/step1_report.py` (numbers only; interpretation after review). Queue `scripts/step1_prob_queue.sh`; outputs `results/step1_prob/`; figure `figures/step1_reliability.png`.

## Summary (plain language)

Every candidate of the all-templates matched filter (≥ 8.0 MAD, nothing cut) gets a probability P(real) that it is a real repeat and not noise. The noise model is the scrambled control: the same matched filter run after shifting each site to a different time of the same day, so nothing lines up across sites except by chance. A classifier learns how data candidates differ from scrambled-control candidates (per-site correlations above each site's own noise level, HF ratio, how many candidates the family and its ≤ 3-km group have within ±1 h, …); the control's share, scaled by how many candidates the data has on tremor-free days, is the local fraction of noise. Fitted on the build half of V2008 + V2010; checked on the V2010 / V2008 test half and on A2011 (a separate episode).

- **V2010 test**: calibration slope own / ≤3 km / any = 0.61 / 0.54 / 0.45 (target 0.8–1.2); matched / ΣP own 1.34, ≤3 km 2.88. Scrambled control with P > 0.5: 0.006 (never-seen control realization 0.011); data with P > 0.5: 0.024. Quiet days: data P > 0.5 0.010 vs control 0.014. Expected count ΣP vs true array count per family-day: Spearman 0.11 (raw candidate count 0.09), ratio ΣP / truth 0.02.
- **A2011**: calibration slope own / ≤3 km / any = 0.71 / 0.81 / 0.80 (target 0.8–1.2); matched / ΣP own 2.15, ≤3 km 2.26. Scrambled control with P > 0.5: 0.015 (never-seen control realization 0.012); data with P > 0.5: 0.031. Quiet days: data P > 0.5 – vs control –. Expected count ΣP vs true array count per family-day: Spearman -0.00 (raw candidate count -0.02), ratio ΣP / truth 0.01.
- **V2008 test**: calibration slope own / ≤3 km / any = 0.53 / 0.53 / 0.42 (target 0.8–1.2); matched / ΣP own 2.13, ≤3 km 2.13. Scrambled control with P > 0.5: 0.003 (never-seen control realization 0.006); data with P > 0.5: 0.029. Quiet days: data P > 0.5 0.001 vs control 0.003. Expected count ΣP vs true array count per family-day: Spearman 0.15 (raw candidate count 0.15), ratio ΣP / truth 0.01.

Caveat for reading calibration: the array catalogue is the truth but is itself incomplete, and a real LFE of a neighbouring family counts as 'real' for P but matches only at the ≤ 3-km / any levels. Matched fraction below P can therefore mean real-but-uncatalogued; matched fraction ABOVE P means P is conservative.

## Method

- Candidates: `results/mf_mt1_all` (T1 all ≥ 30-event build templates, w-weighted, 42 stations, `--drop-dead 0.05`), record tier ≥ 8.0 MAD, series data + control + control2; plus `results/mf_step1c34` = the same run with seed 20261001 → two more independent site-rolled scrambled controls (its data series is identical to mt1_all's: checked per day). K = 4 scrambled-control realizations.
- Features (all computable after the array era): `net_mad`, `net_cc`, `n_sta`, `n_sites`, `wsum`, `site_cc_max`, `site_cc_2nd`, `site_z_max`, `site_z_2nd`, `site_z_min`, `n_z2`, `n_z3`, `n_z5`, `frac_z3`, `dom_top`, `hfr_rel`, `dens_own_1h`, `dens_own_10m`, `dens_grp_1h`, `site_noise_mad_med`, `region_code`, `pnet_max`, `pnet_2nd`, `n_sites_pnet03`, `n_sites_pnet05`. Site values are the within-site w-weighted station cc at the candidate's aligned position (scrambled control: at the site's rolled position, exactly as the MF built it); the recomputed network cc equals the MF's to ≤ 2e-6 on every day. Site z = (site cc − median)/(1.4826 MAD) over 1000 random valid origins of the same family-day. `hfr_rel` = log10 E(10–19)/E(2–8) at the family-day's HF station minus the median of 300 random windows of that station-day. Densities count candidates of the same series.
- Density ratio: HistGradientBoosting + isotonic (3-fold by day) classifier data-vs-control on the fit set; f0/f = [(1−p)/p]·N_data/N_ctl.
- π0 rule: on quiet days (V2010 'quiet', V2008 'post'; tremor-free by the window design) every data candidate is noise → r_q = N_data(quiet) / (N_ctl(quiet)/K) per region (global if < 500 quiet data candidates), π0 = r_q·(N_ctl/K)/N_data, so P = 1 − r_q·(1−p)/(p·K), clipped to [0, 1]. Real LFEs on quiet days count as noise → conservative.
- Fit set: BUILD half of V2008 + V2010 (every series). Held out: V2010 / V2008 TEST half, A2011 (all). Leave-one-control-out refits (K = 3) give the model spread of the rate product and a P for a control realization the model never saw. Matching: within ±2 s (nearest, not one-to-one) of an array detection of the own family / any family ≤ 3 km (lat_use/lon_use) / any family, only on days the family's own array run scanned.

| region | quiet data candidates | quiet control candidates (all K) | r_q | source |
|---|---|---|---|---|
| NW | 1191 | 4121 | 1.156 | region |
| SW | 15 | 34 | 1.056 | global |
| Olympic | 5940 | 22839 | 1.040 | region |
| Arm B | 135 | 587 | 1.056 | global |
| ALL | 7281 | 27581 | 1.056 | global |

Feature ranking (fit set; permutation log-loss gain on 200k candidates; univariate AUC data vs control):

| feature | perm. log-loss gain | AUC data vs control | NaN frac |
|---|---|---|---|
| n_z2 | 0.0724 | 0.590 | 0.000 |
| dom_top | 0.0195 | 0.408 | 0.000 |
| dens_grp_1h | 0.0151 | 0.469 | 0.000 |
| dens_own_10m | 0.0050 | 0.525 | 0.000 |
| hfr_rel | 0.0043 | 0.470 | 0.000 |
| n_sta | 0.0038 | 0.515 | 0.000 |
| pnet_2nd | 0.0021 | 0.546 | 0.000 |
| region_code | 0.0021 | 0.461 | 0.000 |
| site_cc_max | 0.0019 | 0.462 | 0.000 |
| n_sites_pnet03 | 0.0019 | 0.547 | 0.000 |
| site_z_min | 0.0018 | 0.559 | 0.000 |
| dens_own_1h | 0.0018 | 0.499 | 0.000 |
| n_sites | 0.0018 | 0.523 | 0.000 |
| site_z_2nd | 0.0017 | 0.579 | 0.000 |
| net_mad | 0.0013 | 0.516 | 0.000 |
| pnet_max | 0.0010 | 0.537 | 0.000 |
| wsum | 0.0010 | 0.460 | 0.000 |
| site_noise_mad_med | 0.0010 | 0.524 | 0.000 |
| frac_z3 | 0.0009 | 0.561 | 0.000 |
| site_cc_2nd | 0.0009 | 0.549 | 0.000 |
| net_cc | 0.0008 | 0.465 | 0.000 |
| site_z_max | 0.0006 | 0.437 | 0.000 |
| n_sites_pnet05 | 0.0001 | 0.539 | 0.000 |
| n_z5 | -0.0000 | 0.460 | 0.000 |
| n_z3 | -0.0000 | 0.565 | 0.000 |

## (a) Calibration — fraction matched to the array catalogue by P bin

Slope = weighted least-squares slope of matched fraction on bin-mean P (bins with ≥ 20 candidates); matched/ΣP = calibration in the large. Control rows are the chance level.

| set | region | series | n | ΣP | slope own · 3km · any | intercept own · 3km · any | matched/ΣP own · 3km · any |
|---|---|---|---|---|---|---|---|
| V2010_test | ALL | data | 7836 | 518.3 | 0.61 · 0.54 · 0.45 | 0.05 · 0.15 · 0.31 | 1.34 · 2.88 · 5.12 |
| V2010_test | ALL | control | 31602 | 1442.2 | 0.00 · 0.01 · -0.11 | 0.00 · 0.02 · 0.10 | 0.04 · 0.42 · 2.07 |
| V2010_test | NW | data | 1052 | 164.8 | 0.61 · 0.61 · 0.47 | 0.15 · 0.15 · 0.34 | 1.53 · 1.53 · 2.64 |
| V2010_test | NW | control | 3532 | 27.2 | – · – · – | – · – · – | 0.55 · 0.63 · 11.60 |
| V2010_test | Olympic | data | 6577 | 299.9 | 0.01 · 0.19 · 0.13 | 0.06 · 0.17 · 0.32 | 1.31 · 3.96 · 7.08 |
| V2010_test | Olympic | control | 27549 | 1378.2 | 0.00 · 0.01 · -0.11 | 0.00 · 0.02 · 0.10 | 0.03 · 0.43 · 1.91 |
| V2010_test | Arm B | data | 100 | 1.3 | – · – · – | – · – · – | 3.92 · 3.92 · 27.45 |
| V2010_test | Arm B | control | 458 | 13.2 | – · – · – | – · – · – | 0.00 · 0.08 · 2.66 |
| A2011 | ALL | data | 2250 | 221.2 | 0.71 · 0.81 · 0.80 | 0.14 · 0.14 · 0.25 | 2.15 · 2.26 · 3.35 |
| A2011 | ALL | control | 7027 | 174.4 | -0.00 · -0.00 · 0.16 | 0.00 · 0.00 · 0.05 | 0.09 · 0.09 · 2.08 |
| A2011 | NW | data | 2250 | 221.2 | 0.71 · 0.81 · 0.80 | 0.14 · 0.14 · 0.25 | 2.15 · 2.26 · 3.35 |
| A2011 | NW | control | 7027 | 174.4 | -0.00 · -0.00 · 0.16 | 0.00 · 0.00 · 0.05 | 0.09 · 0.09 · 2.08 |
| V2008_test | ALL | data | 596 | 124.6 | 0.53 · 0.53 · 0.42 | 0.33 · 0.34 · 0.47 | 2.13 · 2.13 · 2.67 |
| V2008_test | ALL | control | 1517 | 8.1 | – · – · – | – · – · – | 0.37 · 0.50 · 22.40 |
| V2008_test | NW | data | 593 | 124.3 | 0.53 · 0.53 · 0.41 | 0.34 · 0.34 · 0.48 | 2.13 · 2.14 · 2.68 |
| V2008_test | NW | control | 1508 | 7.5 | – · – · – | – · – · – | 0.40 · 0.53 · 24.06 |
| V2008_test | SW | data | 0 | 0.0 | – · – · – | – · – · – | – · – · – |
| V2008_test | SW | control | 0 | 0.0 | – · – · – | – · – · – | – · – · – |
| fit_in_sample | ALL | data | 11182 | 3802.9 | 0.93 · 0.82 · 0.64 | 0.01 · 0.12 · 0.32 | 0.96 · 1.18 · 1.59 |
| fit_in_sample | ALL | control | 30382 | 683.2 | -0.02 · 0.04 · 0.09 | 0.00 · 0.02 · 0.11 | 0.15 · 0.86 · 5.13 |
| fit_in_sample | NW | data | 3708 | 2364.9 | 0.82 · 0.83 · 0.61 | 0.12 · 0.12 · 0.36 | 1.01 · 1.01 · 1.17 |
| fit_in_sample | NW | control | 5481 | 34.5 | -0.07 · -0.09 · -0.26 | 0.01 · 0.01 · 0.14 | 1.80 · 2.21 · 22.32 |
| fit_in_sample | SW | data | 150 | 149.6 | – · – · – | – · – · – | 1.00 · 1.00 · 1.00 |
| fit_in_sample | SW | control | 0 | 0.0 | – · – · – | – · – · – | – · – · – |
| fit_in_sample | Olympic | data | 6612 | 671.9 | 0.66 · 0.62 · 0.47 | 0.01 · 0.14 · 0.33 | 0.72 · 1.99 · 3.68 |
| fit_in_sample | Olympic | control | 24491 | 630.6 | -0.01 · 0.06 · 0.14 | 0.00 · 0.02 · 0.11 | 0.07 · 0.81 · 4.28 |
| fit_in_sample | Arm B | data | 651 | 579.6 | – · – · – | – · – · – | 1.00 · 1.00 · 1.05 |
| fit_in_sample | Arm B | control | 348 | 3.1 | – · – · – | – · – · – | 0.00 · 0.33 · 9.77 |

Reliability table, V2010_test, ALL regions:

| P bin | data n | mean P | frac own | frac ≤3 km | frac any | control n | control frac any |
|---|---|---|---|---|---|---|---|
| 0.0-0.1 | 6397 | 0.017 | 0.071 | 0.170 | 0.323 | 26719 | 0.098 |
| 0.1-0.2 | 793 | 0.145 | 0.073 | 0.197 | 0.338 | 3346 | 0.073 |
| 0.2-0.3 | 333 | 0.239 | 0.090 | 0.264 | 0.387 | 1095 | 0.079 |
| 0.3-0.4 | 67 | 0.341 | 0.119 | 0.179 | 0.328 | 249 | 0.108 |
| 0.4-0.5 | 37 | 0.445 | 0.243 | 0.270 | 0.405 | 125 | 0.040 |
| 0.5-0.6 | 20 | 0.536 | 0.600 | 0.600 | 0.650 | 38 | 0.053 |
| 0.6-0.7 | 22 | 0.652 | 0.545 | 0.545 | 0.591 | 30 | 0.033 |
| 0.7-0.8 | 37 | 0.752 | 0.676 | 0.676 | 0.757 | 16 | 0.000 |
| 0.8-0.9 | 24 | 0.846 | 0.583 | 0.583 | 0.708 | 5 | 0.000 |
| 0.9-1.0 | 106 | 0.974 | 0.736 | 0.736 | 0.802 | 1 | 0.000 |

Reliability table, A2011, ALL regions:

| P bin | data n | mean P | frac own | frac ≤3 km | frac any | control n | control frac any |
|---|---|---|---|---|---|---|---|
| 0.0-0.1 | 1777 | 0.007 | 0.137 | 0.141 | 0.245 | 6468 | 0.048 |
| 0.1-0.2 | 147 | 0.143 | 0.327 | 0.333 | 0.476 | 290 | 0.083 |
| 0.2-0.3 | 54 | 0.251 | 0.463 | 0.481 | 0.611 | 56 | 0.071 |
| 0.3-0.4 | 44 | 0.352 | 0.295 | 0.341 | 0.477 | 126 | 0.151 |
| 0.4-0.5 | 41 | 0.454 | 0.439 | 0.488 | 0.610 | 65 | 0.108 |
| 0.5-0.6 | 36 | 0.539 | 0.694 | 0.722 | 0.778 | 15 | 0.200 |
| 0.6-0.7 | 28 | 0.644 | 0.786 | 0.786 | 0.857 | 22 | 0.000 |
| 0.7-0.8 | 64 | 0.755 | 0.547 | 0.625 | 0.766 | 1 | 0.000 |
| 0.8-0.9 | 27 | 0.848 | 0.852 | 0.926 | 0.926 | 0 | – |
| 0.9-1.0 | 32 | 0.956 | 0.688 | 0.844 | 0.938 | 0 | – |

Reliability table, V2008_test, ALL regions:

| P bin | data n | mean P | frac own | frac ≤3 km | frac any | control n | control frac any |
|---|---|---|---|---|---|---|---|
| 0.0-0.1 | 439 | 0.006 | 0.335 | 0.337 | 0.469 | 1494 | 0.120 |
| 0.1-0.2 | 28 | 0.151 | 0.429 | 0.429 | 0.607 | 23 | 0.043 |
| 0.2-0.3 | 15 | 0.247 | 0.600 | 0.600 | 0.733 | 7 | 0.143 |
| 0.3-0.4 | 7 | 0.355 | 0.714 | 0.714 | 0.714 | 5 | 0.200 |
| 0.4-0.5 | 14 | 0.456 | 0.714 | 0.714 | 0.714 | 8 | 0.125 |
| 0.5-0.6 | 12 | 0.529 | 0.750 | 0.750 | 0.750 | 2 | 0.000 |
| 0.6-0.7 | 10 | 0.636 | 0.500 | 0.500 | 0.600 | 1 | 0.000 |
| 0.7-0.8 | 24 | 0.745 | 0.792 | 0.792 | 0.833 | 2 | 0.000 |
| 0.8-0.9 | 20 | 0.852 | 0.850 | 0.850 | 0.900 | 0 | – |
| 0.9-1.0 | 85 | 0.975 | 0.824 | 0.824 | 0.847 | 1 | 1.000 |

![reliability](../figures/step1_reliability.png)

## (b) Scrambled control's own P, (c) quiet vs episode days, (d) per region

`ctl P>0.5 (never seen)` = P of each control realization from the leave-that-control-out model. ΣP per realization = expected 'real' count the method would wrongly claim in one scrambled network.

| set | region | days | data n | data P>0.5 | data P>0.9 | data ΣP | ctl n / real. | ctl P>0.1 | ctl P>0.5 | ctl P>0.5 (never seen) | ctl ΣP / real. | ctl ΣP / real. (never seen) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| V2010_test | ALL | all | 10940 | 0.024 | 0.011 | 644.2 | 11269 | 0.153 | 0.006 | 0.011 | 528.3 | 583.09 |
| V2010_test | ALL | quiet | 3055 | 0.010 | 0.000 | 106.5 | 3358 | 0.148 | 0.014 | 0.027 | 163.1 | 196.37 |
| V2010_test | ALL | episode | 7885 | 0.029 | 0.015 | 537.7 | 7911 | 0.155 | 0.003 | 0.005 | 365.2 | 386.72 |
| V2010_test | NW | all | 1567 | 0.127 | 0.066 | 216.7 | 1453 | 0.058 | 0.009 | 0.007 | 31.9 | 32.60 |
| V2010_test | NW | quiet | 455 | 0.015 | 0.000 | 20.9 | 556 | 0.094 | 0.012 | 0.009 | 18.8 | 17.35 |
| V2010_test | NW | episode | 1112 | 0.173 | 0.093 | 195.7 | 896 | 0.035 | 0.006 | 0.006 | 13.1 | 15.25 |
| V2010_test | SW | all | 51 | 0.333 | 0.176 | 19.6 | 9 | 0.189 | 0.054 | 0.108 | 0.7 | 1.14 |
| V2010_test | SW | quiet | 2 | 0.000 | 0.000 | 0.2 | 4 | 0.125 | 0.062 | 0.062 | 0.3 | 0.40 |
| V2010_test | SW | episode | 49 | 0.347 | 0.184 | 19.4 | 5 | 0.238 | 0.048 | 0.143 | 0.4 | 0.74 |
| V2010_test | Olympic | all | 9187 | 0.004 | 0.000 | 400.3 | 9666 | 0.168 | 0.006 | 0.012 | 489.3 | 542.75 |
| V2010_test | Olympic | quiet | 2581 | 0.010 | 0.000 | 85.1 | 2777 | 0.160 | 0.014 | 0.030 | 143.6 | 178.16 |
| V2010_test | Olympic | episode | 6606 | 0.002 | 0.000 | 315.2 | 6889 | 0.171 | 0.002 | 0.005 | 345.7 | 364.59 |
| V2010_test | Arm B | all | 135 | 0.022 | 0.015 | 7.6 | 141 | 0.147 | 0.012 | 0.021 | 6.4 | 6.60 |
| V2010_test | Arm B | quiet | 17 | 0.000 | 0.000 | 0.2 | 20 | 0.062 | 0.012 | 0.012 | 0.4 | 0.46 |
| V2010_test | Arm B | episode | 118 | 0.025 | 0.017 | 7.4 | 121 | 0.161 | 0.012 | 0.023 | 6.0 | 6.14 |
| A2011 | ALL | all | 11814 | 0.031 | 0.003 | 931.7 | 11187 | 0.197 | 0.015 | 0.012 | 695.6 | 728.00 |
| A2011 | ALL | episode | 11814 | 0.031 | 0.003 | 931.7 | 11187 | 0.197 | 0.015 | 0.012 | 695.6 | 728.00 |
| A2011 | NW | all | 2250 | 0.083 | 0.014 | 221.2 | 1761 | 0.082 | 0.005 | 0.007 | 45.8 | 53.57 |
| A2011 | NW | episode | 2250 | 0.083 | 0.014 | 221.2 | 1761 | 0.082 | 0.005 | 0.007 | 45.8 | 53.57 |
| A2011 | SW | all | 9 | 0.000 | 0.000 | 0.8 | 10 | 0.538 | 0.000 | 0.000 | 1.0 | 1.55 |
| A2011 | SW | episode | 9 | 0.000 | 0.000 | 0.8 | 10 | 0.538 | 0.000 | 0.000 | 1.0 | 1.55 |
| A2011 | Olympic | all | 9006 | 0.020 | 0.000 | 685.2 | 8892 | 0.223 | 0.017 | 0.013 | 628.1 | 646.66 |
| A2011 | Olympic | episode | 9006 | 0.020 | 0.000 | 685.2 | 8892 | 0.223 | 0.017 | 0.013 | 628.1 | 646.66 |
| A2011 | Arm B | all | 549 | 0.009 | 0.000 | 24.5 | 525 | 0.151 | 0.008 | 0.021 | 20.7 | 26.22 |
| A2011 | Arm B | episode | 549 | 0.009 | 0.000 | 24.5 | 525 | 0.151 | 0.008 | 0.021 | 20.7 | 26.22 |
| V2008_test | ALL | all | 5280 | 0.029 | 0.016 | 274.2 | 5077 | 0.117 | 0.003 | 0.006 | 166.7 | 182.97 |
| V2008_test | ALL | quiet | 3753 | 0.001 | 0.000 | 102.4 | 3681 | 0.123 | 0.003 | 0.005 | 124.6 | 135.13 |
| V2008_test | ALL | episode | 1527 | 0.099 | 0.056 | 171.8 | 1396 | 0.099 | 0.004 | 0.007 | 42.1 | 47.84 |
| V2008_test | NW | all | 1291 | 0.119 | 0.066 | 159.5 | 1065 | 0.063 | 0.006 | 0.007 | 22.4 | 19.58 |
| V2008_test | NW | quiet | 640 | 0.003 | 0.000 | 9.9 | 682 | 0.082 | 0.007 | 0.006 | 17.7 | 14.69 |
| V2008_test | NW | episode | 651 | 0.232 | 0.131 | 149.6 | 384 | 0.030 | 0.004 | 0.008 | 4.6 | 4.89 |
| V2008_test | SW | all | 6 | 0.000 | 0.000 | 0.6 | 5 | 0.333 | 0.048 | 0.143 | 0.5 | 0.82 |
| V2008_test | SW | quiet | 3 | 0.000 | 0.000 | 0.3 | 3 | 0.333 | 0.083 | 0.167 | 0.4 | 0.50 |
| V2008_test | SW | episode | 3 | 0.000 | 0.000 | 0.3 | 2 | 0.333 | 0.000 | 0.111 | 0.1 | 0.32 |
| V2008_test | Olympic | all | 3748 | 0.000 | 0.000 | 107.8 | 3782 | 0.135 | 0.003 | 0.006 | 139.9 | 157.81 |
| V2008_test | Olympic | quiet | 2943 | 0.000 | 0.000 | 88.0 | 2840 | 0.137 | 0.002 | 0.005 | 104.0 | 116.75 |
| V2008_test | Olympic | episode | 805 | 0.000 | 0.000 | 19.7 | 943 | 0.130 | 0.005 | 0.007 | 35.9 | 41.06 |
| V2008_test | Arm B | all | 235 | 0.000 | 0.000 | 6.4 | 224 | 0.056 | 0.000 | 0.001 | 3.9 | 4.76 |
| V2008_test | Arm B | quiet | 167 | 0.000 | 0.000 | 4.2 | 156 | 0.056 | 0.000 | 0.002 | 2.5 | 3.19 |
| V2008_test | Arm B | episode | 68 | 0.000 | 0.000 | 2.2 | 68 | 0.055 | 0.000 | 0.000 | 1.4 | 1.56 |
| fit_in_sample | ALL | all | 19797 | 0.183 | 0.147 | 4675.8 | 15636 | 0.070 | 0.001 | 0.007 | 381.8 | 592.70 |
| fit_in_sample | ALL | quiet | 7281 | 0.028 | 0.000 | 742.4 | 6895 | 0.090 | 0.002 | 0.008 | 190.0 | 300.02 |
| fit_in_sample | ALL | episode | 12516 | 0.273 | 0.233 | 3933.4 | 8740 | 0.054 | 0.001 | 0.006 | 191.8 | 292.69 |
| fit_in_sample | NW | all | 4899 | 0.493 | 0.423 | 2415.2 | 2413 | 0.041 | 0.004 | 0.008 | 32.3 | 40.69 |
| fit_in_sample | NW | quiet | 1191 | 0.017 | 0.002 | 50.4 | 1030 | 0.052 | 0.002 | 0.005 | 16.8 | 17.25 |
| fit_in_sample | NW | episode | 3708 | 0.646 | 0.559 | 2364.9 | 1382 | 0.033 | 0.005 | 0.009 | 15.5 | 23.44 |
| fit_in_sample | SW | all | 226 | 0.765 | 0.699 | 176.6 | 29 | 0.138 | 0.017 | 0.043 | 1.5 | 2.58 |
| fit_in_sample | SW | quiet | 15 | 0.000 | 0.000 | 0.9 | 8 | 0.088 | 0.000 | 0.029 | 0.2 | 0.51 |
| fit_in_sample | SW | episode | 211 | 0.820 | 0.749 | 175.7 | 20 | 0.159 | 0.024 | 0.049 | 1.3 | 2.07 |
| fit_in_sample | Olympic | all | 13809 | 0.032 | 0.008 | 1485.3 | 12904 | 0.076 | 0.001 | 0.006 | 343.4 | 542.62 |
| fit_in_sample | Olympic | quiet | 5940 | 0.031 | 0.000 | 687.5 | 5710 | 0.099 | 0.001 | 0.009 | 171.5 | 279.64 |
| fit_in_sample | Olympic | episode | 7869 | 0.032 | 0.013 | 797.8 | 7194 | 0.058 | 0.000 | 0.005 | 171.9 | 262.98 |
| fit_in_sample | Arm B | all | 863 | 0.687 | 0.671 | 598.7 | 290 | 0.039 | 0.005 | 0.011 | 4.6 | 6.82 |
| fit_in_sample | Arm B | quiet | 135 | 0.000 | 0.000 | 3.6 | 147 | 0.026 | 0.005 | 0.007 | 1.5 | 2.62 |
| fit_in_sample | Arm B | episode | 728 | 0.815 | 0.795 | 595.0 | 143 | 0.052 | 0.005 | 0.016 | 3.1 | 4.20 |
| ALL | ALL | all | 48368 | 0.091 | 0.065 | 6549.3 | 43722 | 0.130 | 0.006 | 0.009 | 1797.3 | 2115.29 |
| ALL | ALL | quiet | 14328 | 0.017 | 0.000 | 959.6 | 14194 | 0.113 | 0.005 | 0.012 | 487.8 | 644.38 |
| ALL | ALL | episode | 34040 | 0.122 | 0.092 | 5589.6 | 29528 | 0.139 | 0.007 | 0.008 | 1309.5 | 1470.91 |
| ALL | NW | all | 10088 | 0.293 | 0.228 | 3019.8 | 6776 | 0.059 | 0.006 | 0.007 | 134.1 | 148.37 |
| ALL | NW | quiet | 2317 | 0.013 | 0.001 | 83.7 | 2309 | 0.071 | 0.006 | 0.006 | 54.4 | 50.46 |
| ALL | NW | episode | 7771 | 0.377 | 0.295 | 2936.0 | 4468 | 0.053 | 0.005 | 0.008 | 79.8 | 97.91 |
| ALL | SW | all | 292 | 0.651 | 0.572 | 197.6 | 54 | 0.236 | 0.023 | 0.056 | 3.7 | 6.09 |
| ALL | SW | quiet | 20 | 0.000 | 0.000 | 1.5 | 16 | 0.143 | 0.032 | 0.063 | 0.9 | 1.41 |
| ALL | SW | episode | 272 | 0.699 | 0.614 | 196.1 | 38 | 0.275 | 0.020 | 0.052 | 2.9 | 4.68 |
| ALL | Olympic | all | 36197 | 0.018 | 0.003 | 2694.6 | 35700 | 0.145 | 0.007 | 0.010 | 1623.4 | 1916.10 |
| ALL | Olympic | quiet | 11669 | 0.018 | 0.000 | 866.4 | 11538 | 0.124 | 0.005 | 0.013 | 427.9 | 585.96 |
| ALL | Olympic | episode | 24528 | 0.018 | 0.004 | 1828.2 | 24162 | 0.155 | 0.008 | 0.008 | 1195.5 | 1330.14 |
| ALL | Arm B | all | 1791 | 0.336 | 0.324 | 637.3 | 1191 | 0.105 | 0.006 | 0.015 | 36.0 | 44.73 |
| ALL | Arm B | quiet | 322 | 0.000 | 0.000 | 8.0 | 332 | 0.042 | 0.004 | 0.005 | 4.7 | 6.55 |
| ALL | Arm B | episode | 1469 | 0.409 | 0.396 | 629.3 | 860 | 0.129 | 0.007 | 0.019 | 31.3 | 38.18 |

By PNSN regional tremor-day label (label only):

| set | region | tremor day | family-days | data ΣP / fam-day | ctl ΣP / fam-day (per real.) | data P>0.5 | ctl P>0.5 |
|---|---|---|---|---|---|---|---|
| V2010_test | ALL | False | 78 | 1.366 | 2.090 | 0.010 | 0.014 |
| V2010_test | ALL | True | 200 | 2.689 | 1.826 | 0.029 | 0.003 |
| V2010_test | NW | False | 46 | 0.455 | 0.408 | 0.015 | 0.012 |
| V2010_test | NW | True | 117 | 1.673 | 0.112 | 0.173 | 0.006 |
| V2010_test | SW | False | 2 | 0.120 | 0.134 | 0.000 | 0.062 |
| V2010_test | SW | True | 17 | 1.142 | 0.026 | 0.347 | 0.048 |
| V2010_test | Olympic | False | 25 | 3.405 | 5.744 | 0.010 | 0.014 |
| V2010_test | Olympic | True | 49 | 6.432 | 7.054 | 0.002 | 0.002 |
| V2010_test | Arm B | False | 5 | 0.045 | 0.080 | 0.000 | 0.012 |
| V2010_test | Arm B | True | 17 | 0.434 | 0.351 | 0.025 | 0.012 |
| A2011 | ALL | True | 137 | 6.801 | 5.077 | 0.031 | 0.015 |
| A2011 | NW | True | 93 | 2.379 | 0.493 | 0.083 | 0.005 |
| A2011 | SW | True | 7 | 0.111 | 0.146 | 0.000 | 0.000 |
| A2011 | Olympic | True | 27 | 25.378 | 23.261 | 0.020 | 0.017 |
| A2011 | Arm B | True | 10 | 2.450 | 2.068 | 0.009 | 0.008 |
| V2008_test | ALL | False | 204 | 1.344 | 0.817 | 0.029 | 0.003 |
| V2008_test | NW | False | 135 | 1.182 | 0.166 | 0.119 | 0.006 |
| V2008_test | SW | False | 6 | 0.101 | 0.088 | 0.000 | 0.048 |
| V2008_test | Olympic | False | 48 | 2.245 | 2.915 | 0.000 | 0.003 |
| V2008_test | Arm B | False | 15 | 0.424 | 0.263 | 0.000 | 0.000 |

## (e) Expected-count check — ΣP per family-day vs the true array count

Episode family-days the family's array run scanned; scored half only (test half; A2011 all). truth ≤3 km = array detections of the ≤ 3-km group merged within ±2 s. `ΣP ctl` = the same sum on one scrambled-control realization.

| set | region | family-days | ΣP | ΣP ctl | truth own | truth ≤3 km | ΣP/own | ΣP/≤3 km | Spearman own | Spearman ≤3 km | Spearman raw count vs own |
|---|---|---|---|---|---|---|---|---|---|---|---|
| V2010_test | ALL | 4572 | 518.3 | 364.5 | 30921 | 147454 | 0.017 | 0.004 | 0.11 | 0.13 | 0.09 |
| V2010_test | NW | 2640 | 195.7 | 12.9 | 17082 | 37188 | 0.011 | 0.005 | 0.07 | 0.08 | 0.03 |
| V2010_test | Olympic | 784 | 315.2 | 345.7 | 7056 | 85199 | 0.045 | 0.004 | 0.11 | 0.21 | 0.11 |
| V2010_test | Arm B | 1148 | 7.4 | 6.0 | 6783 | 25067 | 0.001 | 0.000 | 0.15 | 0.10 | 0.16 |
| V2008_test | ALL | 4934 | 149.9 | 4.7 | 23962 | 44153 | 0.006 | 0.003 | 0.15 | 0.14 | 0.15 |
| V2008_test | NW | 3036 | 149.6 | 4.6 | 16129 | 25449 | 0.009 | 0.006 | 0.19 | 0.18 | 0.18 |
| V2008_test | SW | 1898 | 0.3 | 0.1 | 7833 | 18704 | 0.000 | 0.000 | 0.02 | 0.01 | 0.02 |
| A2011 | ALL | 3828 | 221.2 | 44.9 | 33306 | 50539 | 0.007 | 0.004 | -0.00 | -0.00 | -0.02 |
| A2011 | NW | 3828 | 221.2 | 44.9 | 33306 | 50539 | 0.007 | 0.004 | -0.00 | -0.00 | -0.02 |

## Rate product — expected real repeats per family and day

`results/step1_prob/rate_family_day.parquet`: one row per observable family-day (all 97 days, V2008 + V2010 + A2011): `sumP` (whole day; build-half P is in-sample), `sumP_scored_half`, `sd_total` = √(Σ P(1−P) + model spread²) with the model spread = SD of ΣP over the 4 leave-one-control-out refits, `sumP_ctl_per_real` (the same sum on one scrambled network), truth counts where the family's array run scanned the day. Tremor-day label is a column only.

| region | days | family-days | mean ΣP / fam-day | mean sd | mean ΣP ctl / fam-day | family-days with ΣP ≥ 5 |
|---|---|---|---|---|---|---|
| Olympic | episode | 2352 | 0.78 | 5.45 | 0.51 | 82 |
| Olympic | post | 84 | 4.06 | 5.92 | 2.27 | 22 |
| Olympic | quiet | 280 | 1.88 | 5.15 | 0.85 | 22 |
| Arm B | episode | 3444 | 0.18 | 0.35 | 0.01 | 47 |
| Arm B | post | 123 | 0.06 | 0.71 | 0.03 | 0 |
| Arm B | quiet | 410 | 0.00 | 0.21 | 0.00 | 0 |
| NW | episode | 11088 | 0.26 | 0.59 | 0.01 | 103 |
| NW | post | 396 | 0.10 | 0.25 | 0.08 | 2 |
| NW | quiet | 1320 | 0.03 | 0.26 | 0.01 | 2 |
| SW | episode | 6132 | 0.03 | 0.31 | 0.00 | 13 |
| SW | post | 219 | 0.01 | 0.22 | 0.00 | 0 |
| SW | quiet | 730 | 0.00 | 0.24 | 0.00 | 0 |

Files: `candidates.parquet` (every candidate, every feature, P, P_loco_c1..c4, matches), `calibration.csv`, `calibration_slopes.csv`, `control_and_quiet.csv`, `tremor_label_split.csv`, `expected_count.csv`, `feature_importance.csv`, `model_info.json`, `rate_family_day.parquet`.
