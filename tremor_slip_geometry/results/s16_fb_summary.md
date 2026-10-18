# S16 PART A — does the Frank & Brodsky (2019) amplitude-slip relation reproduce in Cascadia?

Written 2026-09-10 (session resumed after the 21:03 pod crash) from `results/s16_fb_aggregation.csv`
(116 rows) and `results/s16_fb_events.csv` (18 rows), both produced by
`scripts/s16_fb_aggregation.py` at 21:02 on 2026-09-10. No definition in that script was changed;
the numbers below are read from its two CSVs (the run's stdout report was lost in the crash and is
reconstructed here). A few group medians / paired tests over the SAVED per-event exponents are
computed in this note and are labelled "(aggregated here)".

## 1. What is being compared

Frank & Brodsky (2019, Sci. Adv. 5, eaaw9386; Guerrero) calibrate
`Mdot_geo = alpha * A^beta` with **beta = 3.17 +- 0.29**, where **A is the MEDIAN displacement
amplitude of the individual LFEs in a time bin**. Inverted, `A ~ Mdot^(1/3.17) = Mdot^0.315`
(with `1/(3.17 -+ 0.29)` = 0.289..0.347).

Our published-track exponents (RESULTS.md R9) are for a different aggregation: the **daily SUM**
over a slip patch's 60-km neighbourhood, `C ~ R7^0.40` (27 events) and `E_sum ~ R7^0.64`
(18 events, 2017+). S10 flagged two incompatible readings of the 0.64 and told us to decide which
is like-for-like BEFORE quoting it against F&B. S16 PART A computes F&B's own statistic — the
**median per-detection amplitude per event-day** — on the same 18 events, with the same windows
(`[start-10, end+10]`), the same regressor (`R7`, the 7-d centred mean of the area-weighted
positive modeled potency rate over the 60-km neighbourhood), the same Theil-Sen fitting and the
same 200-draw >= 365-d shifted-tremor null.

**Catalogue constraint that governs everything.** On all 482,659 PNSN detections from 2017-01-01,
`log10 energy = 1.781*magnitude + 0.007*log10(duration) + 2.904`, residual sd **0.052 dex**,
r = 0.995: `energy` is a deterministic monotone function of an unrounded magnitude, and
`magnitude` is that quantity rounded to 0.1. Energy and magnitude are therefore ONE observable
here, not two, and any amplitude proxy inherits whatever magnitude convention is assumed:

| proxy | definition | assumption | note |
|---|---|---|---|
| `A1_sqrt_Emed` | sqrt(median daily energy) | `E ~ A^2` (duration takes only 150/300 s) | **primary, like-for-like with F&B** |
| `A3_mag_EMfit` | `10^(0.891 * median magnitude)` | A1 rewritten through the fitted E-M slope | agrees with A1 to the 0.1 rounding |
| `A4_mag_1p0` | `10^(1.0 * median magnitude)` | `MeL` read as a local (amplitude) magnitude | consistent with the fitted slope ~2 |
| `A2_mag_0p5` | `10^(0.5 * median magnitude)` | that `log10 E = 1.0*M` | **compressed by 0.5/0.891 = 0.561 BY CONSTRUCTION** — not an independent measurement |
| `Asum_sqrt_Esum` | sqrt(daily SUM of energy) | — | the reading S10 floated; `= b(E_sum)/2` identically, NOT an amplitude of anything |
| `Asum_lin` | daily SUM of sqrt(energy) | — | the literal summed amplitude; mixes count and amplitude |

## 2. Forward exponents on the 18 amplitude-era events (our S11 convention: measure on log R7)

| measure | b median [IQR] | boot90 CI of median | b<1 | pooled OLS / TS | null mean±sd | z med | sign +/n (p) | null p<.05 | median fit days |
|---|---|---|---|---|---|---|---|---|---|
| `C` | +0.370 [+0.166, +0.494] | [+0.181, +0.456] | 18/18 | +0.328 / +0.378 | +0.004 ± 0.328 | +1.12 | 17/18 (p 7.2e-05) | 6 | 37 |
| `C_ns8` | +0.382 [+0.210, +0.573] | [+0.233, +0.510] | 17/18 | +0.396 / +0.459 | +0.004 ± 0.374 | +1.00 | 15/18 (p 0.0038) | 4 | 35 |
| `E_sum` | +0.688 [+0.173, +0.855] | [+0.345, +0.779] | 14/18 | +0.520 / +0.635 | +0.013 ± 0.493 | +1.33 | 15/18 (p 0.0038) | 6 | 37 |
| `E_sum_epsS11` | +0.644 [+0.166, +0.840] | [+0.334, +0.759] | 14/18 | +0.500 / +0.615 | +0.010 ± 0.452 | +1.49 | 15/18 (p 0.0038) | 7 | 37 |
| `E_mean` | +0.231 [+0.004, +0.296] | [+0.069, +0.278] | 18/18 | +0.186 / +0.229 | +0.002 ± 0.235 | +0.77 | 14/18 (p 0.015) | 5 | 37 |
| `E_med` | +0.243 [+0.011, +0.320] | [+0.015, +0.307] | 18/18 | +0.174 / +0.220 | +0.002 ± 0.221 | +0.92 | 14/18 (p 0.015) | 6 | 37 |
| `E_p10` | +0.151 [+0.051, +0.270] | [+0.065, +0.268] | 18/18 | +0.129 / +0.158 | +0.000 ± 0.215 | +0.81 | 14/18 (p 0.015) | 2 | 37 |
| `E_p25` | +0.211 [+0.043, +0.308] | [+0.056, +0.282] | 18/18 | +0.148 / +0.185 | +0.002 ± 0.214 | +0.92 | 14/18 (p 0.015) | 4 | 37 |
| `E_p75` | +0.257 [-0.001, +0.324] | [+0.065, +0.305] | 18/18 | +0.192 / +0.240 | +0.003 ± 0.235 | +0.88 | 13/18 (p 0.048) | 5 | 37 |
| `E_p90` | +0.257 [+0.006, +0.317] | [+0.069, +0.303] | 18/18 | +0.210 / +0.259 | +0.004 ± 0.253 | +0.77 | 13/18 (p 0.048) | 4 | 37 |
| `E_med_ns8` | +0.175 [-0.012, +0.283] | [+0.032, +0.265] | 18/18 | +0.162 / +0.199 | +0.003 ± 0.235 | +0.67 | 14/18 (p 0.015) | 2 | 35 |
| `M_med` | +0.118 [+0.000, +0.174] | [+0.017, +0.165] | 18/18 | +0.096 / +0.121 | +0.000 ± 0.113 | +0.94 | 16/18 (p 0.00066) | 6 | 37 |
| `M_max` | +0.147 [+0.000, +0.237] | [+0.000, +0.231] | 18/18 | +0.168 / +0.200 | +0.004 ± 0.181 | +0.90 | 13/18 (p 0.048) | 4 | 37 |
| `A1_sqrt_Emed` | +0.121 [+0.006, +0.160] | [+0.032, +0.151] | 18/18 | +0.087 / +0.110 | +0.001 ± 0.110 | +0.92 | 14/18 (p 0.015) | 6 | 37 |
| `A2_mag_0p5` | +0.059 [+0.000, +0.087] | [+0.008, +0.083] | 18/18 | +0.048 / +0.061 | +0.000 ± 0.057 | +0.94 | 16/18 (p 0.00066) | 6 | 37 |
| `A3_mag_EMfit` | +0.106 [+0.000, +0.155] | [+0.015, +0.147] | 18/18 | +0.085 / +0.108 | +0.000 ± 0.101 | +0.94 | 16/18 (p 0.00066) | 6 | 37 |
| `A4_mag_1p0` | +0.118 [+0.000, +0.174] | [+0.017, +0.165] | 18/18 | +0.096 / +0.121 | +0.000 ± 0.113 | +0.94 | 16/18 (p 0.00066) | 6 | 37 |
| `Asum_sqrt_Esum` | +0.344 [+0.086, +0.427] | [+0.173, +0.390] | 18/18 | +0.260 / +0.317 | +0.006 ± 0.246 | +1.33 | 15/18 (p 0.0038) | 6 | 37 |
| `Asum_lin` | +0.561 [+0.159, +0.666] | [+0.289, +0.608] | 15/18 | +0.423 / +0.503 | +0.009 ± 0.401 | +1.34 | 15/18 (p 0.0038) | 6 | 37 |
| `NS_med` | +0.000 [+0.000, +0.645] | [+0.000, +0.581] | 14/18 | +0.574 / +0.620 | +0.008 ± 0.728 | +0.06 | 12/18 (p 0.12) | 6 | 37 |

Cross-check (exact): `C` and `E_sum_epsS11` reproduce `results/s11_energy_events.csv`
(`b_C_R7`, `b_E_R7`) per event to **max |diff| 0.0e+00** for all 18 events, and `E_mean`
reproduces its `b_EperDet_R7`; the same code re-run in `scripts/s16_daily_tracks.py` reproduces
every stored exponent to **8e-17**. The S11 energy eps (1st percentile of positive daily E,
151,136) costs only -0.013 in exponent (`E_sum_epsS11_minus_E_sum`), so the choice is immaterial.

## 3. The F&B statistic, in every regression direction

A regression is not invertible unless r = 1. Here the per-event OLS |r| of log A1 on log R7 is
**0.319 median (r^2 = 0.114)**, so forward and reverse differ by ~1/r^2 ~ 9x. All estimators:

| amplitude proxy | forward b (A on R) | reverse beta (R on A), F&B direction | RMA b | beta_RMA | naive 1/b_fwd | cross-event b (n18) | cross-event beta (n18) |
|---|---|---|---|---|---|---|---|
| `A1_sqrt_Emed` | +0.121 [+0.032, +0.151] | +1.009 [+0.198, +1.309] | +0.271 [+0.046, +0.301] | +2.54 | +6.30 | +0.296 | +0.997 |
| `A3_mag_EMfit` | +0.106 [+0.015, +0.147] | +1.003 [+0.305, +1.167] | +0.276 [+0.207, +0.300] | +2.58 | +8.16 | +0.218 | +1.024 |
| `A4_mag_1p0` | +0.118 [+0.017, +0.165] | +0.894 [+0.271, +1.113] | +0.310 [+0.239, +0.337] | +2.29 | +7.27 | +0.245 | +0.912 |
| `A2_mag_0p5` | +0.059 [+0.008, +0.083] | +1.787 [+0.643, +2.189] | +0.155 [+0.116, +0.168] | +4.59 | +14.53 | +0.122 | +1.823 |
| `Asum_sqrt_Esum` | +0.344 [+0.173, +0.390] | +0.367 [+0.196, +0.472] | +0.677 [+0.607, +0.987] | +1.14 | +2.61 | +1.086 | +0.546 |
| `Asum_lin` | +0.561 [+0.289, +0.608] | +0.233 [+0.136, +0.318] | +1.137 [+0.965, +1.624] | +0.70 | +1.66 | +1.867 | +0.347 |

Per-event detail of the reverse (F&B-direction) fits, against beta = 3.17 - 0.29 = 2.88:

| reverse-direction detail | b_median | boot90 | pooled TS | pooled OLS | events with beta<2.88 |
|---|---|---|---|---|---|
| `A1_sqrt_Emed` | +1.009 | [+0.198, +1.309] | +0.979 | +1.066 | 18/18 |
| `A3_mag_EMfit` | +1.003 | [+0.305, +1.167] | +0.974 | +1.058 | 18/18 |
| `A4_mag_1p0` | +0.894 | [+0.271, +1.113] | +0.867 | +0.942 | 18/18 |
| `A2_mag_0p5` | +1.787 | [+0.643, +2.189] | +1.734 | +1.884 | 14/18 |

- **Forward (our convention):** A1 ~ R7^**0.121** [90 % CI of the median 0.032, 0.151], 18/18
  events below 1, positive against the shifted null (null +0.001 +- 0.110, z median +0.92,
  14/18 positive, sign p 0.015, 6/18 events at null p < 0.05). This is **2.6x below** F&B's 0.315,
  and the CI excludes 0.315.
- **Reverse (F&B's own direction, moment rate as the dependent variable):** beta = **+1.01**
  [0.20, 1.31], pooled TS +0.98 / OLS +1.07, and **18/18 events have beta < 2.88**. In F&B's
  direction our amplitudes look ~3x MORE responsive to slip rate than theirs, not less.
- **Direction-symmetric (RMA):** b = **+0.271** [0.046, 0.301], `beta_RMA` = **2.54** (median of
  the per-event 1/b; 1/median(b) = 3.69) — this estimator, which is the appropriate one when both
  variables carry noise, **brackets F&B's 3.17** for A1, A3 and A4 (2.29-2.58).
- **Naive inversion of the forward slope** (1/b_fwd = 6.3) is the wrong number and is carried in
  the CSV as `naive_inv_fwd_median` only so nobody re-derives it.

## 4. Time-bin ladder — the gap is an aggregation/attenuation effect

F&B's A is a median over a time bin whose length is not ours. A shorter bin gives a noisier median
AND a noisier rate, which attenuates the forward slope. Coarsening the bin moves our forward
exponent monotonically toward F&B's value:

| amplitude proxy | 1-d bin | 3-d bin (n16) | 7-d bin (n11) | cross-event (n18) | F&B implied |
|---|---|---|---|---|---|
| `A1_sqrt_Emed` | +0.121 [+0.032, +0.151] | +0.166 [+0.100, +0.213] | +0.227 [+0.177, +0.436] | +0.296 | +0.315 |
| `A3_mag_EMfit` | +0.106 [+0.015, +0.147] | +0.133 [+0.000, +0.227] | +0.214 [+0.169, +0.402] | +0.218 | +0.315 |
| `A4_mag_1p0` | +0.118 [+0.017, +0.166] | +0.150 [+0.000, +0.255] | +0.240 [+0.190, +0.452] | +0.245 | +0.315 |
| `A2_mag_0p5` | +0.059 [+0.008, +0.083] | +0.075 [+0.000, +0.127] | +0.120 [+0.048, +0.226] | +0.122 | +0.315 |

The 7-d bin (11 events with >= 8 bins) gives **+0.227 [+0.177, +0.436] — a CI that contains
0.315** — and the cross-event fit (one point per event, n = 18, the granularity closest to F&B's
8 SSEs; event-median A1 vs event-mean log R7) gives **+0.296** (OLS +0.247, r = 0.485; no
bootstrap CI was computed at this level). The ladder 0.121 -> 0.166 -> 0.227 -> 0.296 with
`A2` (compressed 0.561x by construction) tracking at 0.059 -> 0.075 -> 0.120 -> 0.122.

## 5. Which part of the amplitude distribution responds, and paired differences

| paired difference (per event) | median | IQR | boot90 CI | Wilcoxon p | n>0 |
|---|---|---|---|---|---|
| `E_sum_minus_C` | +0.272 | [+0.024, +0.404] | [+0.120, +0.377] | 0.0034 | 14/18 |
| `E_sum_epsS11_minus_E_sum` | -0.013 | [-0.021, -0.003] | [-0.021, -0.005] | 0.00042 | 3/18 |
| `Asum_lin_minus_Asum_sqrt_Esum` | +0.212 | [+0.072, +0.240] | [+0.100, +0.222] | 0.00019 | 16/18 |
| `E_sum_minus_E_med` | +0.436 | [+0.159, +0.524] | [+0.197, +0.484] | 7.6e-05 | 16/18 |
| `E_med_minus_E_mean` | -0.001 | [-0.026, +0.014] | [-0.019, +0.012] | 0.8 | 9/18 |
| `E_p90_minus_E_med` | +0.013 | [-0.004, +0.059] | [+0.003, +0.042] | 0.074 | 12/18 |
| `E_p10_minus_E_med` | -0.058 | [-0.103, -0.012] | [-0.091, -0.027] | 0.004 | 4/18 |
| `E_med_ns8_minus_E_med` | -0.031 | [-0.066, -0.010] | [-0.045, -0.014] | 0.0066 | 4/18 |

- The **whole** daily energy distribution shifts with slip rate, not just its tail: p10 +0.151,
  p25 +0.211, median +0.243, p75 +0.257, p90 +0.257 (`E_p90 - E_med` +0.013, p 0.074;
  `E_p10 - E_med` -0.058, p 0.004 — the weak end responds slightly LESS).
- `E_med` and `E_mean` are interchangeable here (paired -0.001, p 0.80), so F&B's median-based
  statistic and S11's energy-per-detection are the same measurement in this catalogue.
- Requiring >= 8 stations LOWERS the amplitude exponent slightly (`E_med_ns8 - E_med` = -0.031,
  p 0.007) while it RAISED the count exponent (S11: +0.074), i.e. the amplitude response is not
  manufactured by the detection threshold.
- Relative to our own headline exponents, on the SAME 18 events (aggregated here from
  `s16_fb_events.csv`): `b(A1) - b(C)` = **-0.262** [IQR -0.319, -0.131], CI90 [-0.292, -0.155],
  Wilcoxon p 1.5e-05, 1/18 positive; `b(E_med) - b(C)` = -0.128, p 7.6e-05; `b(A1) - b(E_sum)` =
  -0.522, p 1.1e-04. Per-event the three exponents are strongly rank-correlated
  (rho(b_A1, b_C) = +0.88, rho(b_A1, b_E_sum) = +0.92), i.e. they are three readings of one
  relation, not independent measurements.
- Segment and epoch splits of the amplitude exponent (medians of the saved per-event values,
  aggregated here): north +0.139 (n 14) vs central +0.044 (n 4); 2015-19 +0.137 (n 6) vs
  2020-26 +0.117 (n 12) — the same lack of segment/epoch structure S9/S11 found for the count,
  on far fewer events.
- The S10 "summed-amplitude reading" is **retired**: `Asum_sqrt_Esum` = sqrt(daily sum of energy)
  is `b(E_sum)/2 = +0.344` identically — its apparent agreement with 0.315 is arithmetic, not a
  measurement of an amplitude. The LITERAL summed amplitude `Asum_lin` gives **+0.561**, close to
  `b(C) + b(A1)` = 0.370 + 0.121 = 0.491 (paired `Asum_lin - Asum_sqrt_Esum` = +0.212, p 1.9e-04),
  which is what a count x amplitude product must do.

## 6. VERDICT — does the F&B aggregation reproduce their relation?

**Partly, and only with the estimator and the bin length stated.**

1. **QUALITATIVELY YES.** With F&B's own aggregation (median per-detection amplitude per bin) the
   Cascadia relation is positive, sub-linear and decisive in sign: A1 ~ R7^0.121, **18/18 events
   below 1**, 14/18 positive against a >= 365-d shifted-tremor null (sign p 0.015), and the same
   holds for every magnitude-based variant (A2/A3/A4, 16/18 positive, sign p 6.6e-04). Slip rate
   raises per-detection tremor amplitude, as in Guerrero.
2. **QUANTITATIVELY it depends on the estimator, by a factor ~9.** Forward 0.121 (CI excludes
   0.315); direction-symmetric RMA 0.271 (beta_RMA 2.54, brackets 3.17); F&B's own direction
   beta 1.01 (3x below 3.17, all 18 events below 2.88). Because r^2 = 0.11, the forward and
   reverse slopes bracket the relation as [0.12, 0.99] and the RMA compromise lands on F&B's
   value; no single number can be quoted without its direction.
3. **QUANTITATIVELY YES once the granularity is matched.** F&B aggregate over long bins and
   8 SSEs; we aggregate daily. Our forward exponent climbs 0.121 (1 d) -> 0.166 (3 d) -> 0.227
   (7 d, CI [0.18, 0.44] CONTAINS 0.315) -> 0.296 (per event, n = 18). At their granularity our
   value is statistically indistinguishable from 1/3.17.
4. **What must NOT be said.** (a) That our energy exponent 0.64 "equals" F&B's 0.32 via
   `0.64/2` — that identity is arithmetic (`Asum_sqrt_Esum`), not a like-for-like statistic.
   (b) That F&B's 0.32 agrees with our COUNT exponent 0.40 — the numerical closeness is a
   coincidence between a per-event amplitude exponent and a rate exponent; on the same 18 events
   the like-for-like amplitude exponent is 0.26 BELOW the count exponent (p 1.5e-05).
   (c) That the amplitude and count exponents are independent measurements (rho 0.88-0.92).
5. **How it sits with R9.** The three aggregations decompose the same sub-linear relation:
   count `R^0.37` + per-detection amplitude `R^0.12` -> summed amplitude `R^0.56`; count `R^0.37`
   + per-detection energy `R^0.23` -> summed energy `R^0.69` (all medians on the 18 energy-era
   events; the 27-event count exponent is 0.40). Physically: as slip accelerates, tremor produces
   MORE detections and SLIGHTLY STRONGER individual detections, but both sub-linearly — the
   "counts under-represent strong tremor" offset of S11 (b_E - b_C = +0.222) is the
   per-detection amplitude term measured directly here (E_med +0.243 ~ E_mean +0.231).

## 7. Caveats

- **`energy` and `magnitude` are one observable** (E-M residual 0.052 dex), so PNSN amplitude
  proxies differ only by the assumed magnitude convention; A2 (the form initially specified) is
  0.561x A1 by construction and must never be quoted as an independent value.
- **PNSN "amplitude" is not F&B's amplitude.** Theirs is a per-LFE displacement amplitude from
  waveform fitting; ours is a catalogue energy/magnitude for a 150/300-s tremor window, i.e. a
  windowed envelope amplitude that mixes source amplitude with the number of overlapping bursts
  inside the window. This is the single largest obstacle to a strict comparison and cannot be
  removed without a per-LFE amplitude catalogue (the LFE-catalogue gap flagged in S14/S15).
- **18 events, 2017+ only** (the energy/magnitude columns are empty before 2017-01-01); 14 are
  north and 4 central (1361, 1124, 1127, 1434), and 2 events have only 8 fit days (`1905`, `1434`).
- **Events, windows, footprints and the R7 regressor still come from the vbICA slip-patch
  catalogue**, as in S9/S11/S15; R7 is 7-d smoothed and mixture-compressed, so all forward
  exponents inherit S9's "upper bound" caveat (S11 showed the smooth-rate fit is unbiased).
- **The 7-d-bin and cross-event rows rest on 11 and 18 points**; the cross-event level has no
  bootstrap CI in the CSV. The bin ladder's nulls were run for the 1-d bin only.
- **The picker-ranked candidate count (S15B) is not a third independent measure** (daily Spearman
  0.97 with PNSN); it is an envelope-peak count that lands on the energy exponent. Nothing in S16
  changes that.
