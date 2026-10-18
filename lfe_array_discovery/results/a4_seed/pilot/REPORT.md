# A4 location-seeded positive control (pilot)

Written by `scripts/a4_report.py` from `results/a4_seed/pilot`.  Method, choices and caveats: `results/a4_seed/A4_LOG.md` and the docstring of `scripts/a4_seed.py`.

Conditions: **real** = Royer & Bostock (2014) Fig 4 Hyp2000 template locations, real data; **nullT** = same locations, every station reading a different (day, in-day offset) -- the no-family null; **nullR** = the template cloud rotated + shifted (a location-specificity control; its points are a median 2-5 km from some real template, so it is NOT family-free).

AUC = P(a real location scores better than a null location); 0.5 = no separation.  Missing values (no seed, no detections, not located) count as the worst value.  Nothing is cut.

## Counts per subarray and condition

| grp | cond | n | seeded | located | within_3km | within_5km | frac_S_1s_med | snr_med | det_per_day_med |
|---|---|---|---|---|---|---|---|---|---|
| NVI/attr | nullR | 10 | 10 | 0 | 0 | 0 | 0.158 | 2.41 | 51.4 |
| NVI/attr | nullT | 10 | 10 | 0 | 0 | 0 | 0.133 | 2.11 | 27.6 |
| NVI/attr | real | 10 | 10 | 1 | 1 | 1 | 0.2 | 2.6 | 50.7 |
| NVI/top | nullR | 10 | 10 | 0 | 0 | 0 | 0.2 | 2.62 | 54.6 |
| NVI/top | nullT | 10 | 10 | 0 | 0 | 0 | 0.2 | 2.27 | 42.7 |
| NVI/top | real | 10 | 10 | 3 | 2 | 2 | 0.167 | 2.57 | 43 |
| NW/attr | nullR | 10 | 10 | 2 | 0 | 0 | 0.179 | 5.5 | 152 |
| NW/attr | nullT | 10 | 10 | 0 | 0 | 0 | 0.0742 | 3.04 | 29.6 |
| NW/attr | real | 10 | 10 | 1 | 0 | 0 | 0.143 | 5.16 | 131 |
| NW/top | nullR | 10 | 10 | 4 | 0 | 0 | 0.11 | 2.94 | 11.3 |
| NW/top | nullT | 10 | 10 | 0 | 0 | 0 | 0.143 | 2.87 | 40.3 |
| NW/top | real | 10 | 10 | 4 | 0 | 0 | 0.143 | 3.52 | 13.9 |
| SW/attr | nullR | 10 | 10 | 0 | 0 | 0 | 0.133 | 1.79 | 14.7 |
| SW/attr | nullT | 10 | 10 | 0 | 0 | 0 | 0.194 | 2.01 | 10.9 |
| SW/attr | real | 10 | 10 | 3 | 0 | 1 | 0.258 | 1.96 | 11.1 |
| SW/top | nullR | 10 | 10 | 0 | 0 | 0 | 0.167 | 2.02 | 21.6 |
| SW/top | nullT | 10 | 10 | 0 | 0 | 0 | 0.167 | 1.9 | 7.7 |
| SW/top | real | 10 | 10 | 0 | 0 | 0 | 0.0958 | 2.02 | 14.6 |

## Real vs null, criterion by criterion (AUC)

| subarray | metric | better | auc_nullR | auc_nullT | null_median_nullR | null_median_nullT | real_median_nullR | real_median_nullT |
|---|---|---|---|---|---|---|---|---|
| NVI/attr | cc_med_fin | hi | 0.23 | 0.38 | 0.237 | 0.216 | 0.206 | 0.206 |
| NVI/attr | dS_med_abs | lo | 0.33 | 0.66 | 2.84 | 3.88 | 3.56 | 3.56 |
| NVI/attr | day_max_nmad_med_fin | hi | 0.42 | 0.53 | 11.1 | 10.8 | 10.9 | 10.9 |
| NVI/attr | det_per_day_fin | hi | 0.47 | 0.76 | 51.4 | 27.6 | 50.7 | 50.7 |
| NVI/attr | frac_P_1s | hi | 0.525 | 0.405 | 0.233 | 0.276 | 0.233 | 0.233 |
| NVI/attr | frac_S_1s | hi | 0.715 | 0.79 | 0.158 | 0.133 | 0.2 | 0.2 |
| NVI/attr | hyp_d_truth_km | lo | 0.55 | 0.55 | nan | nan | 0.662 | 0.662 |
| NVI/attr | hyp_dist_km | lo | 0.55 | 0.55 | nan | nan | 1.5 | 1.5 |
| NVI/attr | moveout_corr | hi | 0.52 | 0.52 | 0.129 | 0.162 | 0.144 | 0.144 |
| NVI/attr | n_picks_P | hi | 0.55 | 0.455 | 0.5 | 1 | 1 | 1 |
| NVI/attr | n_picks_S | hi | 0.55 | 0.46 | 0 | 0 | 0 | 0 |
| NVI/attr | nmad_med_fin | hi | 0.77 | 0.41 | 8.95 | 9.24 | 9.08 | 9.08 |
| NVI/attr | seed_excess | hi | 0.22 | 0.42 | 0.154 | 0.133 | 0.0991 | 0.0991 |
| NVI/attr | snr_med | hi | 0.62 | 0.9 | 2.41 | 2.11 | 2.6 | 2.6 |
| NVI/top | cc_med_fin | hi | 0.21 | 0.29 | 0.221 | 0.217 | 0.191 | 0.191 |
| NVI/top | dS_med_abs | lo | 0.72 | 0.71 | 3.97 | 4.11 | 3.69 | 3.69 |
| NVI/top | day_max_nmad_med_fin | hi | 0.41 | 0.54 | 10.8 | 10.8 | 10.7 | 10.7 |
| NVI/top | det_per_day_fin | hi | 0.38 | 0.5 | 54.6 | 42.7 | 43 | 43 |
| NVI/top | frac_P_1s | hi | 0.515 | 0.34 | 0.233 | 0.333 | 0.2 | 0.2 |
| NVI/top | frac_S_1s | hi | 0.435 | 0.425 | 0.2 | 0.2 | 0.167 | 0.167 |
| NVI/top | hyp_d_truth_km | lo | 0.65 | 0.65 | nan | nan | 0.825 | 0.825 |
| NVI/top | hyp_dist_km | lo | 0.65 | 0.65 | nan | nan | 1.68 | 1.68 |
| NVI/top | moveout_corr | hi | 0.61 | 0.74 | 0.255 | 0.0245 | 0.295 | 0.295 |
| NVI/top | n_picks_P | hi | 0.705 | 0.55 | 1 | 1 | 1 | 1 |
| NVI/top | n_picks_S | hi | 0.68 | 0.77 | 0 | 0 | 1 | 1 |
| NVI/top | nmad_med_fin | hi | 0.89 | 0.7 | 8.86 | 9.2 | 9.43 | 9.43 |
| NVI/top | seed_excess | hi | 0.44 | 0.72 | 0.169 | 0.124 | 0.169 | 0.169 |
| NVI/top | snr_med | hi | 0.56 | 0.76 | 2.62 | 2.27 | 2.57 | 2.57 |
| NW/attr | cc_med_fin | hi | 0.51 | 0.64 | 0.207 | 0.2 | 0.206 | 0.206 |
| NW/attr | dS_med_abs | lo | 0.5 | 0.5 | 3.76 | 3.56 | 3.51 | 3.51 |
| NW/attr | day_max_nmad_med_fin | hi | 0.43 | 0.86 | 14.2 | 10.2 | 13.8 | 13.8 |
| NW/attr | det_per_day_fin | hi | 0.45 | 0.76 | 152 | 29.6 | 131 | 131 |
| NW/attr | frac_P_1s | hi | 0.68 | 0.57 | 0.371 | 0.407 | 0.445 | 0.445 |
| NW/attr | frac_S_1s | hi | 0.47 | 0.51 | 0.179 | 0.0742 | 0.143 | 0.143 |
| NW/attr | hyp_d_truth_km | lo | 0.45 | 0.55 | 5.03 | nan | 2.41 | 2.41 |
| NW/attr | hyp_dist_km | lo | 0.46 | 0.55 | 21.9 | nan | 5.62 | 5.62 |
| NW/attr | moveout_corr | hi | 0.65 | 0.81 | 0.0916 | -0.0104 | 0.287 | 0.287 |
| NW/attr | n_picks_P | hi | 0.26 | 0.295 | 1 | 2 | 0 | 0 |
| NW/attr | n_picks_S | hi | 0.56 | 0.7 | 0 | 0 | 0 | 0 |
| NW/attr | nmad_med_fin | hi | 0.37 | 0.67 | 9.89 | 8.67 | 9.66 | 9.66 |
| NW/attr | seed_excess | hi | 0.53 | 0.89 | 0.22 | 0.0941 | 0.226 | 0.226 |
| NW/attr | snr_med | hi | 0.47 | 0.79 | 5.5 | 3.04 | 5.16 | 5.16 |
| NW/top | cc_med_fin | hi | 0.57 | 0.15 | 0.161 | 0.201 | 0.16 | 0.16 |
| NW/top | dS_med_abs | lo | 0.63 | 0.96 | 3.12 | 3.9 | 2.77 | 2.77 |
| NW/top | day_max_nmad_med_fin | hi | 0.57 | 0.61 | 10.4 | 10.4 | 10.6 | 10.6 |
| NW/top | det_per_day_fin | hi | 0.77 | 0.24 | 11.3 | 40.3 | 13.9 | 13.9 |
| NW/top | frac_P_1s | hi | 0.52 | 0.795 | 0.536 | 0.371 | 0.593 | 0.593 |
| NW/top | frac_S_1s | hi | 0.55 | 0.57 | 0.11 | 0.143 | 0.143 | 0.143 |
| NW/top | hyp_d_truth_km | lo | 0.53 | 0.7 | 7.44 | nan | 4.85 | 4.85 |
| NW/top | hyp_dist_km | lo | 0.54 | 0.7 | 31 | nan | 18.1 | 18.1 |
| NW/top | moveout_corr | hi | 0.48 | 0.71 | 0.439 | 0.179 | 0.331 | 0.331 |
| NW/top | n_picks_P | hi | 0.49 | 0.91 | 3 | 1 | 3 | 3 |
| NW/top | n_picks_S | hi | 0.445 | 0.4 | 0 | 0 | 0 | 0 |
| NW/top | nmad_med_fin | hi | 0.4 | 0.57 | 8.88 | 8.72 | 8.87 | 8.87 |
| NW/top | seed_excess | hi | 0.67 | 1 | 0.277 | 0.0849 | 0.308 | 0.308 |
| NW/top | snr_med | hi | 0.69 | 0.73 | 2.94 | 2.87 | 3.52 | 3.52 |
| SW/attr | cc_med_fin | hi | 0.24 | 0.4 | 0.18 | 0.172 | 0.167 | 0.167 |
| SW/attr | dS_med_abs | lo | 0.67 | 0.68 | 3.91 | 3.52 | 2.73 | 2.73 |
| SW/attr | day_max_nmad_med_fin | hi | 0.33 | 0.32 | 10.2 | 10.8 | 9.98 | 9.98 |
| SW/attr | det_per_day_fin | hi | 0.46 | 0.46 | 14.7 | 10.9 | 11.1 | 11.1 |
| SW/attr | frac_P_1s | hi | 0.66 | 0.4 | 0.267 | 0.419 | 0.3 | 0.3 |
| SW/attr | frac_S_1s | hi | 0.575 | 0.595 | 0.133 | 0.194 | 0.258 | 0.258 |
| SW/attr | hyp_d_truth_km | lo | 0.65 | 0.65 | nan | nan | 3.58 | 3.58 |
| SW/attr | hyp_dist_km | lo | 0.65 | 0.65 | nan | nan | 9.91 | 9.91 |
| SW/attr | moveout_corr | hi | 0.69 | 0.78 | 0.365 | 0.226 | 0.51 | 0.51 |
| SW/attr | n_picks_P | hi | 0.735 | 0.635 | 0 | 1 | 1.5 | 1.5 |
| SW/attr | n_picks_S | hi | 0.65 | 0.615 | 0 | 0 | 0 | 0 |
| SW/attr | nmad_med_fin | hi | 0.35 | 0.3 | 8.73 | 9.23 | 8.67 | 8.67 |
| SW/attr | seed_excess | hi | 0.7 | 0.8 | 0.0493 | 0.0504 | 0.0744 | 0.0744 |
| SW/attr | snr_med | hi | 0.75 | 0.48 | 1.79 | 2.01 | 1.96 | 1.96 |
| SW/top | cc_med_fin | hi | 0.59 | 0.28 | 0.166 | 0.173 | 0.168 | 0.168 |
| SW/top | dS_med_abs | lo | 0.52 | 0.68 | 3.19 | 4.21 | 3.17 | 3.17 |
| SW/top | day_max_nmad_med_fin | hi | 0.52 | 0.47 | 10.5 | 10.6 | 10.5 | 10.5 |
| SW/top | det_per_day_fin | hi | 0.41 | 0.84 | 21.6 | 7.7 | 14.6 | 14.6 |
| SW/top | frac_P_1s | hi | 0.52 | 0.845 | 0.4 | 0.2 | 0.4 | 0.4 |
| SW/top | frac_S_1s | hi | 0.21 | 0.18 | 0.167 | 0.167 | 0.0958 | 0.0958 |
| SW/top | hyp_d_truth_km | lo | 0.5 | 0.5 | nan | nan | nan | nan |
| SW/top | hyp_dist_km | lo | 0.5 | 0.5 | nan | nan | nan | nan |
| SW/top | moveout_corr | hi | 0.65 | 0.79 | 0.271 | 0.155 | 0.505 | 0.505 |
| SW/top | n_picks_P | hi | 0.76 | 0.71 | 0.5 | 0 | 1 | 1 |
| SW/top | n_picks_S | hi | 0.5 | 0.45 | 0 | 0 | 0 | 0 |
| SW/top | nmad_med_fin | hi | 0.42 | 0.07 | 8.69 | 8.9 | 8.64 | 8.64 |
| SW/top | seed_excess | hi | 0.83 | 1 | 0.181 | 0.0552 | 0.224 | 0.224 |
| SW/top | snr_med | hi | 0.52 | 0.81 | 2.02 | 1.9 | 2.02 | 2.02 |

## Metric definitions

- `hyp_dist_km` (lower = better): Hyp2000 epicentre of the final stack to the seed location (km)
- `hyp_d_truth_km` (lower = better): Hyp2000 epicentre to the NEAREST Royer template (km)
- `frac_S_1s` (higher = better): fraction of stack stations whose horizontal envelope peak is within 1 s of predicted S
- `frac_P_1s` (higher = better): same for the vertical peak before S vs predicted P
- `dS_med_abs` (lower = better): median |observed - predicted S| (s)
- `moveout_corr` (higher = better): corr(observed S, predicted S) across stations
- `n_picks_S` (higher = better): unambiguous (R>=1.5) S picks on the Hilbert-rotated stack
- `n_picks_P` (higher = better): unambiguous (R>=1.5) P picks
- `snr_med` (higher = better): median station stack SNR (envelope peak / median, 30 s stack)
- `seed_excess` (higher = better): seed coherence: mean top-5 similarity minus median candidate similarity
- `det_per_day_fin` (higher = better): final matched-filter detections per day (8 x MAD)
- `nmad_med_fin` (higher = better): median nmad of the final detections
- `cc_med_fin` (higher = better): median network cc of the final detections
- `day_max_nmad_med_fin` (higher = better): median over days of the day's maximum nmad

Figure: `separation.png` (cumulative fraction of locations at least as good as x). Curve table: `curves.csv`.  Per-location table: `all_locations.csv`.
