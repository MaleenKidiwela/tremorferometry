# smoke3

Written by `scripts/a4_report.py` from `results/a4_seed/smoke3`.  Method, choices and caveats: `results/a4_seed/A4_LOG.md` and the docstring of `scripts/a4_seed.py`.

Conditions: **real** = Royer & Bostock (2014) Fig 4 Hyp2000 template locations, real data; **nullT** = same locations, every station reading a different (day, in-day offset) -- the no-family null; **nullR** = the template cloud rotated + shifted (a location-specificity control; its points are a median 2-5 km from some real template, so it is NOT family-free).

AUC = P(a real location scores better than a null location); 0.5 = no separation.  Missing values (no seed, no detections, not located) count as the worst value.  Nothing is cut.

## Counts per subarray and condition

| grp | cond | n | seeded | located | within_3km | within_5km | frac_S_1s_med | snr_med | det_per_day_med |
|---|---|---|---|---|---|---|---|---|---|
| NVI/attr | nullR | 10 | 10 | 1 | 0 | 0 | 0.156 | 2.01 | 104 |
| NVI/attr | nullT | 10 | 10 | 0 | 0 | 0 | 0.133 | 1.84 | 31.5 |
| NVI/attr | real | 10 | 10 | 3 | 2 | 2 | 0.188 | 2.11 | 95 |
| NVI/top | nullR | 10 | 10 | 3 | 0 | 0 | 0.25 | 2.15 | 122 |
| NVI/top | nullT | 10 | 10 | 0 | 0 | 0 | 0.133 | 1.89 | 31.5 |
| NVI/top | real | 10 | 10 | 6 | 4 | 4 | 0.344 | 2.74 | 106 |

## Real vs null, criterion by criterion (AUC)

| subarray | metric | better | auc_nullR | auc_nullT | null_median_nullR | null_median_nullT | real_median_nullR | real_median_nullT |
|---|---|---|---|---|---|---|---|---|
| NVI/attr | cc_med_fin | hi | 0.41 | 0.08 | 0.236 | 0.277 | 0.223 | 0.223 |
| NVI/attr | dS_med_abs | lo | 0.65 | 0.81 | 3.88 | 4.86 | 3.3 | 3.3 |
| NVI/attr | day_max_nmad_med_fin | hi | 0.54 | 0.66 | 16.4 | 15.9 | 18.8 | 18.8 |
| NVI/attr | det_per_day_fin | hi | 0.505 | 0.77 | 104 | 31.5 | 95 | 95 |
| NVI/attr | frac_P_1s | hi | 0.57 | 0.76 | 0.219 | 0.133 | 0.281 | 0.281 |
| NVI/attr | frac_S_1s | hi | 0.615 | 0.62 | 0.156 | 0.133 | 0.188 | 0.188 |
| NVI/attr | hyp_d_truth_km | lo | 0.585 | 0.65 | 1.48 | nan | 1.96 | 1.96 |
| NVI/attr | hyp_dist_km | lo | 0.605 | 0.65 | 7.37 | nan | 2.93 | 2.93 |
| NVI/attr | moveout_corr | hi | 0.37 | 0.4 | 0.23 | 0.273 | 0.0973 | 0.0973 |
| NVI/attr | n_picks_P | hi | 0.86 | 0.86 | 1 | 1 | 2.5 | 2.5 |
| NVI/attr | n_picks_S | hi | 0.52 | 0.615 | 0 | 0 | 0 | 0 |
| NVI/attr | nmad_med_fin | hi | 0.62 | 0.49 | 9.63 | 10.9 | 10 | 10 |
| NVI/attr | seed_excess | hi | 0.46 | 0.5 | 0.0917 | 0.0807 | 0.0922 | 0.0922 |
| NVI/attr | snr_med | hi | 0.66 | 0.96 | 2.01 | 1.84 | 2.11 | 2.11 |
| NVI/top | cc_med_fin | hi | 0.34 | 0 | 0.231 | 0.291 | 0.219 | 0.219 |
| NVI/top | dS_med_abs | lo | 0.71 | 0.89 | 2.57 | 4.38 | 1.78 | 1.78 |
| NVI/top | day_max_nmad_med_fin | hi | 0.65 | 0.72 | 15.8 | 15.2 | 18.2 | 18.2 |
| NVI/top | det_per_day_fin | hi | 0.435 | 0.845 | 122 | 31.5 | 106 | 106 |
| NVI/top | frac_P_1s | hi | 0.72 | 0.84 | 0.312 | 0.2 | 0.469 | 0.469 |
| NVI/top | frac_S_1s | hi | 0.71 | 0.72 | 0.25 | 0.133 | 0.344 | 0.344 |
| NVI/top | hyp_d_truth_km | lo | 0.71 | 0.8 | 4.51 | nan | 1.27 | 1.27 |
| NVI/top | hyp_dist_km | lo | 0.72 | 0.8 | 22 | nan | 2.08 | 2.08 |
| NVI/top | moveout_corr | hi | 0.49 | 0.36 | 0.0953 | 0.196 | 0.0051 | 0.0051 |
| NVI/top | n_picks_P | hi | 0.53 | 0.76 | 1.5 | 0.5 | 2 | 2 |
| NVI/top | n_picks_S | hi | 0.725 | 0.82 | 0 | 0 | 2 | 2 |
| NVI/top | nmad_med_fin | hi | 0.71 | 0.37 | 9.64 | 11.2 | 10 | 10 |
| NVI/top | seed_excess | hi | 0.63 | 0.7 | 0.0987 | 0.0888 | 0.105 | 0.105 |
| NVI/top | snr_med | hi | 0.62 | 0.93 | 2.15 | 1.89 | 2.74 | 2.74 |

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
