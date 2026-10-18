# E5 step 2 -- corrected prediction model

* bg = 2.728; raw p10/50/90 at M1.8-2.1, L=8: 3.39/5.89/9.51 -> sigma_log pooled 0.175 (PRIMARY), within-family robust 0.140 (sd 0.144); g_eq = 0.863
* per-station published-quality excess (L=8, M1.8-2.1): LZB 24 km 2.34, VGZ 36 km 1.15, PGC 38 km 2.15, B011 41 km 1.65, PFB 42 km 0.68, B926 47 km 2.60, LRIV 52 km 0.54, YOUB 52 km 1.87, SNB 56 km 0.98
* class model: surface <40 km 2.15, 40-60 km 0.83; borehole 2.12 (B011/B926 at 41-47 km, both bands); >60 km: far x (60/d)^2. e(12 s)/e(8 s) = 1.100. q-ratio ours(N=90)/published = 0.72 -> r1_eff = 0.088; q(25/50/100/200/500) = 0.40/0.53/0.66/0.78/0.89
* IN-SAMPLE (actual station sets, station-measured e_s, equal weight): predicted raw median 5.47 vs measured 5.89 (bar 5.9+-0.5); frac>=8 0.165 vs 0.206 (bar 0.21+-0.05) -> PASS
* IN-SAMPLE with the CLASS model instead of station values: raw 5.35, frac>=8 0.161; with the WITHIN-family spread (0.140) frac>=8 = 0.114 (under-predicts the tail: between-family variance). Excess calibration factor measured/predicted = 1.15 (reported as a labelled bracket, NOT applied to the primary numbers)
* MEASURED recall vs magnitude, 2010 ring L=8, >=8 MAD: POOLED is flat/non-monotonic (M1.35:0.18, M1.45:0.23, M1.55:0.22, M1.65:0.25, M1.75:0.23, M1.85:0.21, M1.95:0.20, M2.05:0.19, M2.15:0.19, M2.25:0.24); WITHIN-family (centred) rises 0.157 at dM -0.35 -> 0.272 at dM +0.35 (d log10 recall / dM = 0.36). The per-bin 2012 numbers inherit the model slope; alpha = 0 is the flat bracket.
* Bostock_21_templated: median 11 stations/family, e_net 5.16 -> Mc(8) = 1.98 (target 1.6-1.9; his 2003-05 modal M 1.6-1.7)
* Bostock_21_all_within60: median 19 stations/family, e_net 6.27 -> Mc(8) = 1.75 (target 1.6-1.9; his 2003-05 modal M 1.6-1.7)
* ring7_2010-13_templated (ring + SOKB): 3 stations/family, e_net 2.58 -> Mc(8) = 2.77 vs his 2010-13 modal 1.75: implementation-vs-Bostock deficit +1.02 mag = 2.44x in registration (alpha 0.38); i.e. his 7-station matcher needs 2.05x our predicted network excess to put M1.75 at 8 MAD
* ring7_2010-13_all_within60 (ring + SOKB): 5 stations/family, e_net 3.16 -> Mc(8) = 2.54 vs his 2010-13 modal 1.75: implementation-vs-Bostock deficit +0.79 mag = 1.99x in registration (alpha 0.38); i.e. his 7-station matcher needs 1.67x our predicted network excess to put M1.75 at 8 MAD
* 2012 set15: 130 families active 2012-08-15..10-15 (per-family n p10/50/90 [ 49.8 192.5 392. ]); 126 get >= 1 template (build-half n >= 20); stations within 60 km per family median 11 (p10 7, p90 12); q(build) median 0.66; e_net equal median 3.92, snrw 4.04; Mc(8) equal median 2.29, snrw 2.25
*   sensitivity set15 L=8 snrw T=8 alpha=0.38 g=g0.86: recall all 0.279 | M<1.6 0.184 M1.6-1.8 0.273 M1.8-2.0 0.355 M>=2.0 0.429
*   sensitivity set15 L=8 snrw T=8 alpha=0.38 g=g1/1.38: recall all 0.209 | M<1.6 0.135 M1.6-1.8 0.202 M1.8-2.0 0.269 M>=2.0 0.335
*   sensitivity set15 L=8 snrw T=8 alpha=0.38 g=g0.86xCAL1.15: recall all 0.347 | M<1.6 0.236 M1.6-1.8 0.342 M1.8-2.0 0.436 M>=2.0 0.511
*   sensitivity set15 L=8 snrw T=8 alpha=0.46 g=g0.86: recall all 0.264 | M<1.6 0.158 M1.6-1.8 0.253 M1.8-2.0 0.348 M>=2.0 0.444
*   sensitivity set15 L=8 snrw T=8 alpha=0.46 g=g1/1.38: recall all 0.198 | M<1.6 0.115 M1.6-1.8 0.186 M1.8-2.0 0.263 M>=2.0 0.349
*   sensitivity set15 L=8 snrw T=8 alpha=0.46 g=g0.86xCAL1.15: recall all 0.328 | M<1.6 0.204 M1.6-1.8 0.319 M1.8-2.0 0.428 M>=2.0 0.527
*   sensitivity set15 L=8 snrw T=8 alpha=0.58 g=g0.86: recall all 0.243 | M<1.6 0.125 M1.6-1.8 0.225 M1.8-2.0 0.339 M>=2.0 0.467
*   sensitivity set15 L=8 snrw T=8 alpha=0.58 g=g1/1.38: recall all 0.182 | M<1.6 0.091 M1.6-1.8 0.165 M1.8-2.0 0.256 M>=2.0 0.371
*   sensitivity set15 L=8 snrw T=8 alpha=0.58 g=g0.86xCAL1.15: recall all 0.303 | M<1.6 0.162 M1.6-1.8 0.285 M1.8-2.0 0.417 M>=2.0 0.550
*   sensitivity set15 L=8 snrw T=8 alpha=0.0 g=g0.86: recall all 0.378 | M<1.6 0.367 M1.6-1.8 0.385 M1.8-2.0 0.388 M>=2.0 0.360
*   sensitivity set15 L=8 snrw T=8 alpha=0.0 g=g1/1.38: recall all 0.289 | M<1.6 0.280 M1.6-1.8 0.294 M1.8-2.0 0.296 M>=2.0 0.275
*   sensitivity set15 L=8 snrw T=8 alpha=0.0 g=g0.86xCAL1.15: recall all 0.460 | M<1.6 0.448 M1.6-1.8 0.468 M1.8-2.0 0.472 M>=2.0 0.438
*   sensitivity set15 L=8 snrw T=8 alpha=0.2 g=g0.86: recall all 0.321 | M<1.6 0.259 M1.6-1.8 0.323 M1.8-2.0 0.370 M>=2.0 0.395
*   sensitivity set15 L=8 snrw T=8 alpha=0.2 g=g1/1.38: recall all 0.242 | M<1.6 0.192 M1.6-1.8 0.242 M1.8-2.0 0.282 M>=2.0 0.305
*   sensitivity set15 L=8 snrw T=8 alpha=0.2 g=g0.86xCAL1.15: recall all 0.397 | M<1.6 0.325 M1.6-1.8 0.399 M1.8-2.0 0.453 M>=2.0 0.476
*   set15                    L=   8s equal T=   8: recall all 0.267 | M<1.6 0.175 M1.6-1.8 0.259 M1.8-2.0 0.340 M>=2.0 0.414
*   set15                    L=   8s snrw  T=   8: recall all 0.279 | M<1.6 0.184 M1.6-1.8 0.273 M1.8-2.0 0.355 M>=2.0 0.429
*   set15                    L=  12s equal T=   8: recall all 0.310 | M<1.6 0.207 M1.6-1.8 0.303 M1.8-2.0 0.391 M>=2.0 0.469
*   set15                    L=  12s snrw  T=   8: recall all 0.324 | M<1.6 0.218 M1.6-1.8 0.318 M1.8-2.0 0.408 M>=2.0 0.483
*   set18_with_LZB_YOUB_PFB  L=   8s equal T=   8: recall all 0.298 | M<1.6 0.193 M1.6-1.8 0.289 M1.8-2.0 0.380 M>=2.0 0.477
*   set18_with_LZB_YOUB_PFB  L=   8s snrw  T=   8: recall all 0.316 | M<1.6 0.206 M1.6-1.8 0.307 M1.8-2.0 0.401 M>=2.0 0.496
*   set18_with_LZB_YOUB_PFB  L=  12s equal T=   8: recall all 0.345 | M<1.6 0.228 M1.6-1.8 0.336 M1.8-2.0 0.436 M>=2.0 0.536
*   set18_with_LZB_YOUB_PFB  L=  12s snrw  T=   8: recall all 0.364 | M<1.6 0.242 M1.6-1.8 0.356 M1.8-2.0 0.457 M>=2.0 0.556
* region-year (q at N=96 chain-start, alpha 0.38, sigma pooled): see e5_prediction_regionyear.csv
* DOWNLOAD (no fetch): 15 stations x 3 ch x 51 d (2012-08-15..10-05) at 40 Hz = 11.2 GB (F1-measured 1.416 B/sample at 40 Hz; E4's 1.17 B/sample would give 9.3 GB); native rates (12 stations at 100 Hz, 3 at 40 Hz) = 17.0 GB
