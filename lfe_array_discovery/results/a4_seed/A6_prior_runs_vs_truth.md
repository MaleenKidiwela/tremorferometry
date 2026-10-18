# A6 -- earlier located Washington family tables vs the Royer & Bostock (2014) Hyp2000 templates

Script `scripts/a6_prior_vs_truth.py`; scorer `scripts/ab_truth_score.py` (unchanged: one-to-one
Hungarian match on epicentral distance, candidates = our families within 30 km of the subarray's
templates, null = our candidate cloud randomly rotated + shifted, 1000 draws).  Truth = Fig 4 Hyp2000
rows (NW 122, SW 54).  `frac_Xkm` = fraction of TRUTH templates matched one-to-one within X km.
Rows: every located row of each table ("all located"), plus the dedupe-kept subset where the table
has one.  The chain's own site-shift null tables get the identical scoring.
Caveat: a dense cloud of hundreds of our points inside the 30 km region matches many truth templates
by chance, which the rotation null measures -- read `frac` against `null_mean`/`p`, and real chain
against null chain.

| bed | chain | subset | subarray | n_ours_in_region | median_match_km | frac_3km | null_mean_3km | p_3km | frac_5km | null_mean_5km | null_p95_5km | p_5km | frac_10km | null_mean_10km | p_10km | n_depth_pairs | dz_median_km |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2010 N-Olympic arrays (agent O / sf_final) | real | all located | NW | 482 | 2.7 | 0.549 | 0.429 | 0.218 | 0.746 | 0.657 | 0.869 | 0.383 | 0.951 | 0.874 | 0.394 | 114 | -23.3 |
| 2010 N-Olympic arrays (agent O / sf_final) | real | all located | SW | 243 | 12.3 | 0.111 | 0.159 | 0.653 | 0.222 | 0.269 | 0.556 | 0.57 | 0.426 | 0.423 | 0.459 | 23 | -8.92 |
| 2010 N-Olympic arrays (agent O / sf_final) | real | dedupe-kept | NW | 252 | 7.7 | 0.246 | 0.198 | 0.381 | 0.369 | 0.355 | 0.598 | 0.488 | 0.615 | 0.61 | 0.531 | 73 | -23 |
| 2010 N-Olympic arrays (agent O / sf_final) | real | dedupe-kept | SW | 107 | 26.1 | 0.0556 | 0.057 | 0.509 | 0.13 | 0.119 | 0.278 | 0.415 | 0.259 | 0.246 | 0.394 | 14 | -12.1 |
| 2010 N-Olympic arrays (agent O / sf_final) | null | all located | NW | 305 | 5.9 | 0.254 | 0.228 | 0.444 | 0.451 | 0.401 | 0.656 | 0.441 | 0.721 | 0.662 | 0.474 | 86 | -23.9 |
| 2010 N-Olympic arrays (agent O / sf_final) | null | all located | SW | 165 | 14.3 | 0.111 | 0.135 | 0.573 | 0.241 | 0.248 | 0.537 | 0.46 | 0.352 | 0.431 | 0.621 | 19 | -18.1 |
| 2010 N-Olympic arrays (agent O / sf_final) | null | dedupe-kept | NW | 260 | 7.38 | 0.197 | 0.184 | 0.463 | 0.352 | 0.333 | 0.574 | 0.475 | 0.631 | 0.588 | 0.48 | 75 | -24.1 |
| 2010 N-Olympic arrays (agent O / sf_final) | null | dedupe-kept | SW | 143 | 16.4 | 0.111 | 0.117 | 0.5 | 0.222 | 0.221 | 0.481 | 0.431 | 0.296 | 0.402 | 0.688 | 16 | -17.8 |
| 2010 N-Olympic, u_locate L9 relocation | real | all located | NW | 602 | 2.31 | 0.656 | 0.43 | 0.146 | 0.803 | 0.609 | 0.91 | 0.248 | 0.967 | 0.79 | 0.191 | 116 | 1.77 |
| 2010 N-Olympic, u_locate L9 relocation | real | all located | SW | 230 | 34.7 | 0.0556 | 0.0668 | 0.579 | 0.0926 | 0.108 | 0.278 | 0.545 | 0.185 | 0.179 | 0.39 | 10 | 33.4 |
| 2010 N-Olympic, u_locate L9 relocation | real | dedupe-kept | NW | 303 | 4.19 | 0.344 | 0.24 | 0.3 | 0.574 | 0.394 | 0.73 | 0.261 | 0.828 | 0.616 | 0.281 | 99 | -0.528 |
| 2010 N-Olympic, u_locate L9 relocation | real | dedupe-kept | SW | 103 | 49.1 | 0 | 0.0237 | 1 | 0.0185 | 0.0449 | 0.13 | 0.7 | 0.0741 | 0.0917 | 0.514 | 4 | 23.9 |
| 2010 N-Olympic, u_locate L9 relocation | null | all located | NW | 170 | 24.8 | 0.041 | 0.0532 | 0.425 | 0.0656 | 0.1 | 0.344 | 0.498 | 0.213 | 0.211 | 0.366 | 25 | 34.3 |
| 2010 N-Olympic, u_locate L9 relocation | null | all located | SW | 72 | 43.7 | 0.037 | 0.0132 | 0.173 | 0.037 | 0.0284 | 0.0926 | 0.354 | 0.0556 | 0.0718 | 0.524 | 3 | 29.2 |
| 2010 N-Olympic, u_locate L9 relocation | null | dedupe-kept | NW | 134 | 28.8 | 0.0328 | 0.034 | 0.383 | 0.0492 | 0.0675 | 0.221 | 0.479 | 0.139 | 0.152 | 0.397 | 17 | 35.1 |
| 2010 N-Olympic, u_locate L9 relocation | null | dedupe-kept | SW | 60 | 50.5 | 0.037 | 0.00883 | 0.106 | 0.037 | 0.0184 | 0.0741 | 0.247 | 0.037 | 0.0509 | 0.552 | 2 | 3.2 |
| 2010 Arm B Chestler chain (sf_final) | real | all located | NW | 481 | 3.53 | 0.402 | 0.352 | 0.427 | 0.639 | 0.521 | 0.852 | 0.368 | 0.902 | 0.72 | 0.269 | 108 | 14.6 |
| 2010 Arm B Chestler chain (sf_final) | real | all located | SW | 211 | 32.3 | 0.0185 | 0.0822 | 0.897 | 0.0926 | 0.131 | 0.333 | 0.601 | 0.167 | 0.215 | 0.524 | 9 | -9.08 |
| 2010 Arm B Chestler chain (sf_final) | null | all located | NW | 30 | 25.1 | 0.0164 | 0.0221 | 0.556 | 0.0246 | 0.0386 | 0.107 | 0.595 | 0.0656 | 0.0656 | 0.459 | 7 | 19 |
| 2010 Arm B Chestler chain (sf_final) | null | all located | SW | 9 | 24.6 | 0 | 0.00826 | 1 | 0.0185 | 0.0178 | 0.0556 | 0.558 | 0.037 | 0.0402 | 0.547 | 2 | 23.7 |
| 2010 Arm B, u_locate L9 | real | all located | NW | 560 | 3.11 | 0.484 | 0.4 | 0.356 | 0.77 | 0.591 | 0.861 | 0.234 | 0.984 | 0.808 | 0.102 | 118 | 33.1 |
| 2010 Arm B, u_locate L9 | real | all located | SW | 283 | 30.5 | 0.0926 | 0.108 | 0.522 | 0.111 | 0.168 | 0.408 | 0.638 | 0.259 | 0.264 | 0.431 | 14 | 19.8 |
| 2008 Puget/S-Olympic (agent N) | real | all located | NW | 73 | 37.7 | 0.0082 | 0.0175 | 0.695 | 0.0164 | 0.0334 | 0.115 | 0.671 | 0.041 | 0.0743 | 0.632 | 4 | 18.7 |
| 2008 Puget/S-Olympic (agent N) | real | all located | SW | 99 | 11 | 0.0741 | 0.0949 | 0.676 | 0.204 | 0.197 | 0.352 | 0.5 | 0.444 | 0.422 | 0.499 | 24 | 21.1 |
| 2008 Puget/S-Olympic (agent N) | real | dedupe-kept | NW | 59 | 33.8 | 0.0164 | 0.0199 | 0.503 | 0.0328 | 0.0371 | 0.115 | 0.458 | 0.0574 | 0.075 | 0.517 | 6 | 19.5 |
| 2008 Puget/S-Olympic (agent N) | real | dedupe-kept | SW | 77 | 14.4 | 0.0556 | 0.0619 | 0.588 | 0.13 | 0.129 | 0.259 | 0.502 | 0.296 | 0.302 | 0.515 | 16 | 20.2 |
| 2008 Puget/S-Olympic (agent N) | null | all located | NW | 18 | 23.6 | 0.0082 | 0.0145 | 0.819 | 0.0082 | 0.0249 | 0.0574 | 0.94 | 0.0082 | 0.0429 | 0.991 | 1 | 18.2 |
| 2008 Puget/S-Olympic (agent N) | null | all located | SW | 20 | 12.8 | 0.0185 | 0.0353 | 0.852 | 0.037 | 0.0708 | 0.13 | 0.913 | 0.148 | 0.145 | 0.579 | 8 | 18.8 |
| 2008 Puget/S-Olympic (agent N) | null | dedupe-kept | NW | 18 | 23.6 | 0.0082 | 0.0145 | 0.819 | 0.0082 | 0.0249 | 0.0574 | 0.94 | 0.0082 | 0.0429 | 0.991 | 1 | 18.2 |
| 2008 Puget/S-Olympic (agent N) | null | dedupe-kept | SW | 20 | 12.8 | 0.0185 | 0.0353 | 0.852 | 0.037 | 0.0708 | 0.13 | 0.913 | 0.148 | 0.145 | 0.579 | 8 | 18.8 |
| 2010 WA (agent L, l_wa2010) | real | all located | NW | 874 | 2.56 | 0.566 | 0.526 | 0.473 | 0.852 | 0.688 | 0.926 | 0.236 | 0.984 | 0.829 | 0.133 | 120 | -14.1 |
| 2010 WA (agent L, l_wa2010) | real | all located | SW | 411 | 40.2 | 0.0926 | 0.118 | 0.556 | 0.0926 | 0.166 | 0.407 | 0.706 | 0.13 | 0.244 | 0.688 | 7 | -13.1 |
| 2010 WA (agent L, l_wa2010) | null | all located | NW | 884 | 2.35 | 0.598 | 0.518 | 0.391 | 0.795 | 0.662 | 0.934 | 0.327 | 0.951 | 0.796 | 0.249 | 116 | -17.8 |
| 2010 WA (agent L, l_wa2010) | null | all located | SW | 479 | 40.2 | 0.111 | 0.12 | 0.47 | 0.111 | 0.158 | 0.389 | 0.607 | 0.13 | 0.223 | 0.652 | 7 | -13.1 |
| 2010 refire (m_2011) | real | all located | NW | 25 | 9.86 | 0.0328 | 0.0332 | 0.579 | 0.0492 | 0.0566 | 0.0984 | 0.652 | 0.107 | 0.0913 | 0.398 | 12 | 8.19 |
| 2010 refire (m_2011) | real | all located | SW | 17 | 16.7 | 0.0185 | 0.0163 | 0.505 | 0.0926 | 0.0346 | 0.111 | 0.12 | 0.0926 | 0.0784 | 0.404 | 5 | -0.793 |
| 2010 refire (m_2011) | null | all located | NW | 10 | 9.7 | 0.0164 | 0.013 | 0.484 | 0.0246 | 0.0222 | 0.041 | 0.545 | 0.041 | 0.038 | 0.534 | 5 | 19.1 |
| 2010 refire (m_2011) | null | all located | SW | 5 | 3.71 | 0 | 0.00615 | 1 | 0.0741 | 0.0147 | 0.037 | 0.00899 | 0.0741 | 0.036 | 0.14 | 4 | 0.428 |
| 2011 refire (m_2011) | real | all located | NW | 26 | 7.04 | 0.041 | 0.0388 | 0.506 | 0.0738 | 0.0666 | 0.107 | 0.483 | 0.148 | 0.107 | 0.0989 | 18 | -22.1 |
| 2011 refire (m_2011) | real | all located | SW | 14 | 19.8 | 0 | 0.0132 | 1 | 0.0185 | 0.0295 | 0.0926 | 0.671 | 0.0556 | 0.064 | 0.533 | 3 | -15.5 |
| 2011 refire (m_2011) | null | all located | NW | 3 | 4.71 | 0.0082 | 0.00493 | 0.536 | 0.0164 | 0.00907 | 0.0164 | 0.274 | 0.0246 | 0.0154 | 0.143 | 3 | -14.9 |
| 2011 refire (m_2011) | null | all located | SW | 1 | 7.65 | 0 | 0.00174 | 1 | 0 | 0.00422 | 0.0185 | 1 | 0.0185 | 0.0106 | 0.575 | 1 | -11.3 |

## Reading (A6 verdict)

- **No earlier located table matches the Royer templates better than its own rotation null** at
  any radius: every p >= 0.1 (lowest: 2010 WA refire, 2011 real NW 10 km p=0.10; u_locate L9 real
  NW 3 km 0.656 vs null mean 0.430, p=0.15).  The matched fractions are high in NW (0.5-0.85 within
  3-5 km) only because hundreds of our points fill the 30 km region; the rotated cloud does almost as well.
- **Chain real vs chain null gives the same picture:** l_wa2010 NW 3 km 0.566 (real) vs 0.598 (site-shift
  null chain); sf_final oly2010 NW 0.549 vs 0.254 -- but the null table there is smaller (305 vs 482 in
  region), and against the rotation null the real table is still not significant (p=0.22).
- u_locate L9 real NW is the only relocation with depth residuals near zero (median +1.8 km, dedupe-kept
  -0.5 km); the grid locators of sf_final / l_wa2010 sit 14-23 km shallow or deep of the Royer depths.
- SW: our 2010-era tables barely reach SW (median match 11-49 km); the 2008 Puget run (agent N) does
  reach SW but matches at the null level (0.20 vs 0.20 at 5 km).
- So the earlier runs neither show nor exclude Royer's families: their location clouds carry no
  information about where the templates are, beyond occupying the same region.
