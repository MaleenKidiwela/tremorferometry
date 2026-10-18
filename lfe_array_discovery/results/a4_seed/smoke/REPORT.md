# smoke

Written by `scripts/a4_report.py` from `results/a4_seed/smoke`.  Method, choices and caveats: `results/a4_seed/A4_LOG.md` and the docstring of `scripts/a4_seed.py`.

Conditions: **real** = Royer & Bostock (2014) Fig 4 Hyp2000 template locations, real data; **nullT** = same locations, every station reading a different (day, in-day offset) -- the no-family null; **nullR** = the template cloud rotated + shifted (a location-specificity control; its points are a median 2-5 km from some real template, so it is NOT family-free).

AUC = P(a real location scores better than a null location); 0.5 = no separation.  Missing values (no seed, no detections, not located) count as the worst value.  Nothing is cut.

## Counts per subarray and condition

| subarray | cond | n | seeded | located | within_3km | within_5km | frac_S_1s_med | snr_med | det_per_day_med |
|---|---|---|---|---|---|---|---|---|---|
| NVI | nullT | 10 | 10 | 0 | 0 | 0 | 0.133 | 1.89 | 31.5 |
| NVI | real | 10 | 10 | 6 | 4 | 4 | 0.344 | 2.74 | 106 |

## Real vs null, criterion by criterion (AUC)

| subarray | metric | better | auc_nullT | null_median_nullT | real_median_nullT |
|---|---|---|---|---|---|
| NVI | cc_med_fin | hi | 0 | 0.291 | 0.219 |
| NVI | dS_med_abs | lo | 0.89 | 4.38 | 1.78 |
| NVI | day_max_nmad_med_fin | hi | 0.72 | 15.2 | 18.2 |
| NVI | det_per_day_fin | hi | 0.845 | 31.5 | 106 |
| NVI | frac_P_1s | hi | 0.84 | 0.2 | 0.469 |
| NVI | frac_S_1s | hi | 0.72 | 0.133 | 0.344 |
| NVI | hyp_d_truth_km | lo | 0.8 | nan | 1.27 |
| NVI | hyp_dist_km | lo | 0.8 | nan | 2.08 |
| NVI | moveout_corr | hi | 0.36 | 0.196 | 0.0051 |
| NVI | n_picks_P | hi | 0.76 | 0.5 | 2 |
| NVI | n_picks_S | hi | 0.82 | 0 | 2 |
| NVI | nmad_med_fin | hi | 0.37 | 11.2 | 10 |
| NVI | seed_excess | hi | 0.7 | 0.0888 | 0.105 |
| NVI | snr_med | hi | 0.93 | 1.89 | 2.74 |

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
