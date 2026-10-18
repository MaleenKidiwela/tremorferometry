# CONFIDENCE — a per-repeat confidence fitted against the KNOWN ANSWER (array-era truth)

Written by `scripts/confidence_report.py` from `results/confidence/` (fit: `scripts/confidence_fit.py`, held-out-station PNET features: `scripts/confidence_pnet_ho.py`). Scrambled control = the same matched filter with each site's data shifted to a different time of the same day.

## Summary (plain language)

Every matched-filter repeat (step-1 candidate, ≥ 8.0 MAD, nothing cut) gets a confidence P = the probability that it matches the array catalogue, learned directly from the known answer rather than from the data-vs-scrambled-control contrast used in step 1. Inputs: the step-1 features, the step-1 P, PNET (off-the-shelf and fine-tuned on Washington) and the step-2 attribution. One model per label: own family (for firing rates), any family ≤ 3 km, any family (for "is it a real LFE").

**One design change from the spec, and why.** The spec asked for a fit on the build half. But the templates are stacks of build-half array detections, so on the build half a catalogued repeat partly correlates with itself: the raw feature n_z2 separates matched from unmatched with AUC 0.92 on the build half but only 0.67 on the test half, and the own-family match rate is 55–61 % there vs 27–29 % (NW). A model fitted there learns an over-optimistic mapping. It was run as specified (**D1**), and it comes out over-confident on independent data (A2011 calibration slope 0.64). The design to use is **D2**: fit on the TEST half with 5 folds grouped by day, so every test-half row is scored by a model that never saw its day. A2011, a different year, is scored by the model fitted on the whole test half and is the fully independent check. The PNET fine-tune used build-half windows at 36 stations only, so it has no test-half, A2011 or held-out-station leakage (verified from its training metadata).

- **Separation (D2, own family), vs step-1 P / fine-tuned PNET alone / net_mad alone**: A2011 (independent) AUC **0.87** vs 0.64 / 0.78 / 0.68, average precision 0.61 vs 0.38. V2010 test: NW 0.85 (step-1 P 0.71), Olympic 0.90 (0.52), Arm B 0.82 (0.79). V2008 test: NW 0.87 (0.71), SW 0.89 (0.82; only 296 candidates). The day-block bootstrap gain over step-1 P is +0.16 to +0.29 AUC, and every interval excludes 0. Label "any": A2011 0.80 vs 0.57; V2010 0.78 (NW 0.78, Olympic 0.76, Arm B 0.66).
- **Within one family-day** (does P rank the repeats of one family on one day?): 0.73–0.85 for the model in NW / Olympic / SW (Arm B 0.57) vs **0.50–0.63 for step-1 P** and 0.61–0.73 for PNET alone. Step-1 P mainly told busy family-days from quiet ones; the new P also ranks individual repeats.
- **Calibration (D2)**: A2011 slope own / any **0.93 / 0.98** (target 0.8–1.2), but calibration in the large (matched / ΣP) is 0.76 / 0.70. P runs about 25–30 % high on A2011 because its own-family match rate (14 %) is half the training rate, so a base-rate shift needs a per-epoch correction. On the test half (out-of-fold) the slope is 0.95–1.19, except Arm B (own 1.40, any 0.74; only 86 / 498 positives). Step-1 P for comparison: 0.21–0.68.
- **Scrambled control (c)**, own model: P > 0.5 for **0.7 % of control vs 2.2 % of data** (V2010 test, all regions); NW 3.7 % vs 13.6 %; A2011 NW 3.5 % vs 11.1 %; Olympic 0.3 % vs 0.4 %, i.e. no separation, as in step 1. **But the control's total ΣP is 50–65 % of the data's in NW and ~100 % in Olympic.** Much of the low-level P comes from context (a family's busy hours), not from evidence in the individual repeat. The "any family" model is worse: control P > 0.5 in 10–18 % vs data 12–24 %. It has largely learned "busy hour" and must not be used alone as a per-repeat "is it real" score. Adding the scrambled controls as known negatives (**D3**) lowers control P > 0.5 (NW own 3.6 → 2.3 %, any 13.6 → 11.1 %), costs 0.005 AUC and leaves ΣP unchanged.
- **Quiet vs episode days (d)**: on quiet days data and control look the same (NW P > 0.5 0.1 % vs 0.1 %; ΣP ratio 0.99), so there is no false excess on tremor-free days. All the separation comes from episode days (NW 18 % vs 5 %).
- **Expected count (e)**: ΣP per family-day vs the true own-family array count, Spearman V2010 0.52 / V2008 0.46 / A2011 0.47, vs step-1 ΣP 0.49 / 0.45 / 0.48 and the raw candidate count 0.46 / 0.43 / 0.40. Only a small gain, except Olympic (0.75 vs 0.37). ΣP / truth is 0.16–0.35, which reflects recall: most array detections never produce a permanent-station candidate. Subtracting the scrambled-control ΣP makes the ranking worse (ρ 0.34 at A2011), because the context part tracks activity too.
- **Does PNET contribute? (f)** Yes. The **fine-tuned PNET is the most important input**. Permuting it on A2011 drops AUC by 0.16 (own) and 0.20 (any), and removing all PNET features drops A2011 AUC from 0.870 to 0.781 (own) and from 0.804 to 0.699 (any). The off-the-shelf PNET adds nothing (≈ 0.001). **The gain is station-specific**: with PNET read only at the 6 stations the fine-tune never saw, A2011 own AUC is 0.778, no better than without PNET (0.781). Fine-tuned PNET alone at those stations gives AUC 0.58–0.69 vs 0.48–0.65 off-the-shelf. The step-2 attribution adds nothing (AUC drop 0.003); the step-1 P adds a little (0.02).
- **Caveats**: (1) the catalogue is incomplete, so P = P(matches the catalogue), a lower bound on P(real); (2) most of the P mass is family/hour context, so a per-repeat P > 0.5 is informative while a ΣP firing rate still carries a context floor, which the scrambled control measures; (3) the fit is array-era (2008–2011), and after 2012 the station roster, instrument responses, PNET on changed stations, base rate and density features all need re-checking (list at the end).


## Method

- Candidates: every step-1 candidate (T1 all-templates MF, record tier ≥ 8.0 MAD, nothing cut) on V2008 / V2010 / A2011, `results/step1_prob/candidates.parquet`; 4 site-rolled scrambled controls c1–c4.
- Features: the 21 step-1 features; the step-1 lfdr P (`P1`; on the build half cross-fitted by day, so no build row gets a P from a model fitted on its own day); PNET off-the-shelf (`_w`) and WA fine-tuned (`_f`): pnet_max, pnet_2nd, n_sites_pnet03, n_sites_pnet05, pnet_rel_max; step 2: group size, P(member) and n eligible members of the nearest group-MF candidate within ± 2 s (group frame t + δ_f; exists for data and c1 only; NaN = none / singleton). Tremor-day labels are NOT features.
- Labels (data candidates on days the family's array run scanned): matched within ± 2 s (nearest) to an array detection of the own family / any family of the ≤ 3-km group / any family. One model per label.
- Model: HistGradientBoosting (15 leaves, ≥ 300 per leaf, L2 = 2, early stopping) + isotonic calibration (3 folds grouped by day); logistic regression (C = 0.1, + isotonic) for comparison.
- **D1 (as specified)**: fit on the BUILD half of V2008 + V2010; evaluated on the V2010 / V2008 TEST half and A2011.
- **D2 (template-leak-free)**: the templates are stacks of BUILD-half array detections, so on the build half a truth-matched candidate is partly correlating with itself (table below). D2 fits on the TEST half with 5 outer folds grouped by day (every test-half row, data and scrambled control, is scored by a model that never saw its day) and scores A2011 with a model fitted on the whole test half. A2011 is independent in both designs.
- PNET leakage: the fine-tune used BUILD-half windows of V2008 + V2010 only (`pnet_labels.build_ok`: window inside one build 2-h block), 36 training stations; the test half, A2011 and the 6 held-out stations (HDW, OSD, B018, B017, NLWA, GNW) never entered it (`results/pnet/train/meta.parquet`: 36 stations, none held out; no A2011 date). `PNET_HO` = the same model with PNET read only at those 6 stations.

### Why D2: build half vs test half (own-family label, raw features, AUC matched vs unmatched)

| window   | half   | region   |      n |   frac_own |   frac_any |   auc_own_net_mad |   auc_own_n_z2 |   auc_own_site_cc_max |   auc_own_P1_orig |   auc_own_P1 |   auc_own_pnet_max_f |   auc_own_pnet_max_w |
|:---------|:-------|:---------|-------:|-----------:|-----------:|------------------:|---------------:|----------------------:|------------------:|-------------:|---------------------:|---------------------:|
| V2008    | build  | ALL      |  29471 |      0.648 |      0.710 |             0.767 |          0.941 |                 0.313 |             0.963 |        0.961 |                0.607 |                0.597 |
| V2008    | build  | NW       |  22618 |      0.548 |      0.629 |             0.738 |          0.916 |                 0.315 |             0.946 |        0.943 |                0.681 |                0.615 |
| V2008    | build  | SW       |   6853 |      0.975 |      0.978 |             0.964 |          0.948 |                 0.531 |             0.989 |        0.991 |                0.634 |                0.647 |
| V2008    | test   | ALL      |  13943 |      0.277 |      0.403 |             0.632 |          0.669 |                 0.469 |             0.711 |        0.711 |                0.777 |                0.667 |
| V2008    | test   | NW       |  13647 |      0.272 |      0.399 |             0.633 |          0.665 |                 0.463 |             0.705 |        0.705 |                0.781 |                0.667 |
| V2008    | test   | SW       |    296 |      0.503 |      0.571 |             0.800 |          0.558 |                 0.723 |             0.819 |        0.819 |                0.745 |                0.682 |
| V2010    | build  | ALL      | 130720 |      0.157 |      0.343 |             0.656 |          0.916 |                 0.440 |             0.926 |        0.920 |                0.655 |                0.647 |
| V2010    | build  | NW       |  20590 |      0.607 |      0.709 |             0.767 |          0.923 |                 0.330 |             0.949 |        0.946 |                0.666 |                0.624 |
| V2010    | build  | Olympic  | 105984 |      0.049 |      0.254 |             0.839 |          0.786 |                 0.799 |             0.806 |        0.791 |                0.681 |                0.517 |
| V2010    | build  | Arm B    |   4146 |      0.685 |      0.798 |             0.912 |          0.975 |                 0.100 |             0.992 |        0.990 |                0.598 |                0.626 |
| V2010    | test   | ALL      | 115301 |      0.067 |      0.266 |             0.639 |          0.639 |                 0.657 |             0.626 |        0.626 |                0.745 |                0.614 |
| V2010    | test   | NW       |  12527 |      0.292 |      0.471 |             0.654 |          0.677 |                 0.483 |             0.710 |        0.710 |                0.753 |                0.657 |
| V2010    | test   | Olympic  | 101356 |      0.040 |      0.239 |             0.796 |          0.531 |                 0.843 |             0.519 |        0.519 |                0.709 |                0.505 |
| V2010    | test   | Arm B    |   1418 |      0.061 |      0.351 |             0.635 |          0.761 |                 0.297 |             0.785 |        0.785 |                0.674 |                0.563 |

`P1_orig` = step-1 P as published (in-sample on the build half); `P1` = cross-fitted (used as the feature).

## (a) AUC and average precision on held-out rows — label: own family (firing-rate use)

| eval_set   | region   |      n |   n_pos |   base_rate |   auc_D2_hgb_FULL |   auc_D1_hgb_FULL |   auc_P1 |   auc_pnet_max_f |   auc_pnet_max_w |   auc_net_mad |   ap_D2_hgb_FULL |   ap_D1_hgb_FULL |   ap_P1 |   ap_pnet_max_f |   ap_net_mad |
|:-----------|:---------|-------:|--------:|------------:|------------------:|------------------:|---------:|-----------------:|-----------------:|--------------:|-----------------:|-----------------:|--------:|----------------:|-------------:|
| V2010_test | ALL      | 115301 |    7761 |       0.067 |             0.918 |             0.900 |    0.626 |            0.745 |            0.614 |         0.639 |            0.549 |            0.472 |   0.266 |           0.215 |        0.137 |
| V2010_test | NW       |  12527 |    3662 |       0.292 |             0.853 |             0.786 |    0.710 |            0.753 |            0.657 |         0.654 |            0.730 |            0.634 |   0.569 |           0.503 |        0.441 |
| V2010_test | Olympic  | 101356 |    4013 |       0.040 |             0.900 |             0.884 |    0.519 |            0.709 |            0.505 |         0.796 |            0.327 |            0.284 |   0.057 |           0.114 |        0.158 |
| V2010_test | Arm B    |   1418 |      86 |       0.061 |             0.824 |             0.835 |    0.785 |            0.674 |            0.563 |         0.635 |            0.419 |            0.381 |   0.336 |           0.227 |        0.183 |
| V2008_test | ALL      |  13943 |    3857 |       0.277 |             0.870 |             0.795 |    0.711 |            0.777 |            0.667 |         0.632 |            0.732 |            0.621 |   0.549 |           0.519 |        0.378 |
| V2008_test | NW       |  13647 |    3708 |       0.272 |             0.868 |             0.793 |    0.705 |            0.781 |            0.667 |         0.633 |            0.724 |            0.615 |   0.538 |           0.518 |        0.375 |
| V2008_test | SW       |    296 |     149 |       0.503 |             0.893 |             0.762 |    0.819 |            0.745 |            0.682 |         0.800 |            0.889 |            0.746 |   0.799 |           0.712 |        0.817 |
| A2011      | ALL      |  61999 |    8865 |       0.143 |             0.870 |             0.795 |    0.640 |            0.780 |            0.729 |         0.682 |            0.610 |            0.469 |   0.379 |           0.316 |        0.254 |
| A2011      | NW       |  61999 |    8865 |       0.143 |             0.870 |             0.795 |    0.640 |            0.780 |            0.729 |         0.682 |            0.610 |            0.469 |   0.379 |           0.316 |        0.254 |

## (a) AUC and average precision on held-out rows — label: any family ('is it a real LFE')

| eval_set   | region   |      n |   n_pos |   base_rate |   auc_D2_hgb_FULL |   auc_D1_hgb_FULL |   auc_P1 |   auc_pnet_max_f |   auc_pnet_max_w |   auc_net_mad |   ap_D2_hgb_FULL |   ap_D1_hgb_FULL |   ap_P1 |   ap_pnet_max_f |   ap_net_mad |
|:-----------|:---------|-------:|--------:|------------:|------------------:|------------------:|---------:|-----------------:|-----------------:|--------------:|-----------------:|-----------------:|--------:|----------------:|-------------:|
| V2010_test | ALL      | 115301 |   30626 |       0.266 |             0.776 |             0.790 |    0.511 |            0.656 |            0.537 |         0.632 |            0.593 |            0.603 |   0.330 |           0.432 |        0.389 |
| V2010_test | NW       |  12527 |    5903 |       0.471 |             0.779 |             0.773 |    0.636 |            0.726 |            0.626 |         0.626 |            0.776 |            0.756 |   0.651 |           0.687 |        0.598 |
| V2010_test | Olympic  | 101356 |   24225 |       0.239 |             0.764 |             0.781 |    0.477 |            0.638 |            0.504 |         0.669 |            0.535 |            0.556 |   0.237 |           0.372 |        0.389 |
| V2010_test | Arm B    |   1418 |     498 |       0.351 |             0.661 |             0.686 |    0.511 |            0.585 |            0.510 |         0.576 |            0.527 |            0.524 |   0.394 |           0.445 |        0.438 |
| V2008_test | ALL      |  13943 |    5620 |       0.403 |             0.830 |             0.818 |    0.659 |            0.783 |            0.675 |         0.569 |            0.783 |            0.749 |   0.618 |           0.700 |        0.457 |
| V2008_test | NW       |  13647 |    5451 |       0.399 |             0.828 |             0.817 |    0.654 |            0.785 |            0.675 |         0.567 |            0.778 |            0.746 |   0.610 |           0.700 |        0.453 |
| V2008_test | SW       |    296 |     169 |       0.571 |             0.894 |             0.801 |    0.803 |            0.754 |            0.694 |         0.807 |            0.914 |            0.828 |   0.832 |           0.776 |        0.863 |
| A2011      | ALL      |  61999 |   16123 |       0.260 |             0.804 |             0.770 |    0.565 |            0.767 |            0.723 |         0.598 |            0.645 |            0.576 |   0.411 |           0.524 |        0.334 |
| A2011      | NW       |  61999 |   16123 |       0.260 |             0.804 |             0.770 |    0.565 |            0.767 |            0.723 |         0.598 |            0.645 |            0.576 |   0.411 |           0.524 |        0.334 |

## (a) AUC and average precision on held-out rows — label: ≤ 3-km group

| eval_set   | region   |      n |   n_pos |   base_rate |   auc_D2_hgb_FULL |   auc_D1_hgb_FULL |   auc_P1 |   auc_pnet_max_f |   auc_pnet_max_w |   auc_net_mad |   ap_D2_hgb_FULL |   ap_D1_hgb_FULL |   ap_P1 |   ap_pnet_max_f |   ap_net_mad |
|:-----------|:---------|-------:|--------:|------------:|------------------:|------------------:|---------:|-----------------:|-----------------:|--------------:|-----------------:|-----------------:|--------:|----------------:|-------------:|
| V2010_test | ALL      | 115301 |   16425 |       0.142 |             0.834 |             0.825 |    0.542 |            0.666 |            0.550 |         0.649 |            0.513 |            0.487 |   0.243 |           0.283 |        0.248 |
| V2010_test | NW       |  12527 |    3686 |       0.294 |             0.851 |             0.800 |    0.711 |            0.754 |            0.657 |         0.652 |            0.726 |            0.648 |   0.573 |           0.507 |        0.443 |
| V2010_test | Olympic  | 101356 |   12647 |       0.125 |             0.821 |             0.824 |    0.484 |            0.631 |            0.504 |         0.686 |            0.430 |            0.432 |   0.131 |           0.228 |        0.250 |
| V2010_test | Arm B    |   1418 |      92 |       0.065 |             0.813 |             0.831 |    0.787 |            0.653 |            0.550 |         0.641 |            0.393 |            0.386 |   0.338 |           0.222 |        0.188 |
| V2008_test | ALL      |  13943 |    3892 |       0.279 |             0.871 |             0.808 |    0.713 |            0.778 |            0.668 |         0.630 |            0.738 |            0.640 |   0.559 |           0.523 |        0.380 |
| V2008_test | NW       |  13647 |    3739 |       0.274 |             0.869 |             0.806 |    0.707 |            0.782 |            0.667 |         0.632 |            0.730 |            0.633 |   0.547 |           0.522 |        0.377 |
| V2008_test | SW       |    296 |     153 |       0.517 |             0.895 |             0.767 |    0.829 |            0.741 |            0.689 |         0.802 |            0.897 |            0.774 |   0.825 |           0.717 |        0.829 |
| A2011      | ALL      |  61999 |    8959 |       0.145 |             0.872 |             0.801 |    0.641 |            0.780 |            0.729 |         0.680 |            0.624 |            0.496 |   0.385 |           0.319 |        0.254 |
| A2011      | NW       |  61999 |    8959 |       0.145 |             0.872 |             0.801 |    0.641 |            0.780 |            0.729 |         0.680 |            0.624 |            0.496 |   0.385 |           0.319 |        0.254 |

Day-block bootstrap (200 resamples of days) of AUC(model) − AUC(step-1 P), all regions:

| eval_set   | label   | design   |   n_days |   dAUC_mean |   dAUC_lo |   dAUC_hi |
|:-----------|:--------|:---------|---------:|------------:|----------:|----------:|
| V2010_test | own     | D1       |       29 |       0.273 |     0.245 |     0.301 |
| V2010_test | own     | D2       |       29 |       0.294 |     0.267 |     0.321 |
| V2010_test | any     | D1       |       29 |       0.278 |     0.253 |     0.305 |
| V2010_test | any     | D2       |       29 |       0.263 |     0.235 |     0.294 |
| V2008_test | own     | D1       |       26 |       0.084 |     0.063 |     0.104 |
| V2008_test | own     | D2       |       26 |       0.160 |     0.140 |     0.181 |
| V2008_test | any     | D1       |       26 |       0.159 |     0.134 |     0.181 |
| V2008_test | any     | D2       |       26 |       0.171 |     0.146 |     0.193 |
| A2011      | own     | D1       |       29 |       0.157 |     0.127 |     0.190 |
| A2011      | own     | D2       |       29 |       0.228 |     0.190 |     0.265 |
| A2011      | any     | D1       |       29 |       0.204 |     0.171 |     0.243 |
| A2011      | any     | D2       |       29 |       0.238 |     0.209 |     0.276 |

### Ablations (AUC; D = design, hgb / logit, feature set)

FULL = all features; NO_PNET = without any PNET feature; NO_STEP2 = without step-2 attribution; NO_P1 = without the step-1 P; PNET_HO = PNET only from the 6 never-fine-tuned stations; S1_ONLY = the 21 step-1 features only.

| eval_set   | region   | label   |   D1_hgb_FULL |   D2_hgb_FULL |   D1_logit_FULL |   D2_logit_FULL |   D2_hgb_NO_PNET |   D2_hgb_NO_STEP2 |   D2_hgb_NO_P1 |   D2_hgb_PNET_HO |   D2_hgb_S1_ONLY |   D1_hgb_NO_PNET |   D1_hgb_NO_STEP2 |
|:-----------|:---------|:--------|--------------:|--------------:|----------------:|----------------:|-----------------:|------------------:|---------------:|-----------------:|-----------------:|-----------------:|------------------:|
| V2010_test | ALL      | own     |         0.900 |         0.918 |           0.868 |           0.898 |            0.910 |             0.918 |          0.918 |            0.911 |            0.910 |            0.892 |             0.900 |
| V2010_test | ALL      | any     |         0.790 |         0.776 |           0.754 |           0.756 |            0.763 |             0.770 |          0.775 |            0.766 |            0.755 |            0.780 |             0.786 |
| V2010_test | NW       | own     |         0.786 |         0.853 |           0.742 |           0.817 |            0.819 |             0.852 |          0.849 |            0.829 |            0.816 |            0.749 |             0.784 |
| V2010_test | NW       | any     |         0.773 |         0.779 |           0.735 |           0.755 |            0.735 |             0.779 |          0.776 |            0.748 |            0.728 |            0.731 |             0.772 |
| V2008_test | ALL      | own     |         0.795 |         0.870 |           0.742 |           0.826 |            0.844 |             0.869 |          0.869 |            0.843 |            0.844 |            0.758 |             0.793 |
| V2008_test | ALL      | any     |         0.818 |         0.830 |           0.776 |           0.794 |            0.795 |             0.829 |          0.828 |            0.792 |            0.788 |            0.780 |             0.817 |
| V2008_test | NW       | own     |         0.793 |         0.868 |           0.739 |           0.823 |            0.841 |             0.867 |          0.867 |            0.839 |            0.841 |            0.756 |             0.792 |
| V2008_test | NW       | any     |         0.817 |         0.828 |           0.777 |           0.791 |            0.791 |             0.827 |          0.827 |            0.788 |            0.786 |            0.779 |             0.817 |
| A2011      | ALL      | own     |         0.795 |         0.870 |           0.700 |           0.760 |            0.781 |             0.869 |          0.876 |            0.778 |            0.805 |            0.703 |             0.795 |
| A2011      | ALL      | any     |         0.770 |         0.804 |           0.749 |           0.779 |            0.699 |             0.803 |          0.809 |            0.690 |            0.711 |            0.665 |             0.777 |
| A2011      | NW       | own     |         0.795 |         0.870 |           0.700 |           0.760 |            0.781 |             0.869 |          0.876 |            0.778 |            0.805 |            0.703 |             0.795 |
| A2011      | NW       | any     |         0.770 |         0.804 |           0.749 |           0.779 |            0.699 |             0.803 |          0.809 |            0.690 |            0.711 |            0.665 |             0.777 |

### Held-out PNET stations: families with HDW / OSD / B018 / B017 / NLWA / GNW among their 2 nearest stations

`pnet_max_ho_*` = PNET read only at the 6 held-out stations (never used in the fine-tune).

| eval_set   | region   | label   |     n |   n_pos |   auc_pnet_max_ho_f |   auc_pnet_max_ho_w |   auc_pnet_max_f |   auc_pnet_max_w |   auc_D2_hgb_PNET_HO |   auc_D2_hgb_NO_PNET |   auc_D2_hgb_FULL |   auc_P1 |
|:-----------|:---------|:--------|------:|--------:|--------------------:|--------------------:|-----------------:|-----------------:|---------------------:|---------------------:|------------------:|---------:|
| V2010_test | ALL      | own     | 13603 |    1803 |               0.623 |               0.594 |            0.764 |            0.601 |                0.863 |                0.860 |             0.871 |    0.635 |
| V2010_test | ALL      | any     | 13603 |    5032 |               0.567 |               0.555 |            0.707 |            0.542 |                0.772 |                0.769 |             0.786 |    0.530 |
| V2010_test | NW       | own     |  2117 |     819 |               0.658 |               0.624 |            0.674 |            0.592 |                0.804 |                0.788 |             0.829 |    0.700 |
| V2010_test | NW       | any     |  2117 |    1188 |               0.638 |               0.597 |            0.680 |            0.595 |                0.728 |                0.698 |             0.766 |    0.644 |
| V2010_test | Olympic  | own     | 10874 |     953 |               0.500 |               0.508 |            0.737 |            0.505 |                0.832 |                0.832 |             0.842 |    0.530 |
| V2010_test | Olympic  | any     | 10874 |    3709 |               0.522 |               0.525 |            0.693 |            0.500 |                0.772 |                0.772 |             0.784 |    0.483 |
| V2010_test | Arm B    | own     |   612 |      31 |               0.583 |               0.455 |            0.581 |            0.383 |                0.733 |                0.737 |             0.701 |    0.744 |
| V2010_test | Arm B    | any     |   612 |     135 |               0.558 |               0.535 |            0.594 |            0.502 |                0.678 |                0.654 |             0.670 |    0.539 |
| V2008_test | ALL      | own     |  2131 |     835 |               0.689 |               0.653 |            0.655 |            0.592 |                0.836 |                0.839 |             0.848 |    0.711 |
| V2008_test | ALL      | any     |  2131 |    1220 |               0.662 |               0.621 |            0.702 |            0.614 |                0.742 |                0.751 |             0.792 |    0.642 |
| V2008_test | NW       | own     |  1920 |     727 |               0.669 |               0.630 |            0.665 |            0.589 |                0.826 |                0.833 |             0.841 |    0.696 |
| V2008_test | NW       | any     |  1920 |    1100 |               0.642 |               0.600 |            0.704 |            0.607 |                0.719 |                0.733 |             0.781 |    0.633 |
| V2008_test | SW       | own     |   211 |     108 |               0.794 |               0.795 |            0.713 |            0.675 |                0.895 |                0.881 |             0.878 |    0.820 |
| V2008_test | SW       | any     |   211 |     120 |               0.814 |               0.803 |            0.728 |            0.694 |                0.896 |                0.878 |             0.881 |    0.804 |
| A2011      | ALL      | own     |  9152 |    1532 |               0.582 |               0.484 |            0.725 |            0.637 |                0.799 |                0.795 |             0.874 |    0.679 |
| A2011      | ALL      | any     |  9152 |    2935 |               0.577 |               0.489 |            0.733 |            0.664 |                0.670 |                0.672 |             0.796 |    0.590 |
| A2011      | NW       | own     |  9152 |    1532 |               0.582 |               0.484 |            0.725 |            0.637 |                0.799 |                0.795 |             0.874 |    0.679 |
| A2011      | NW       | any     |  9152 |    2935 |               0.577 |               0.489 |            0.733 |            0.664 |                0.670 |                0.672 |             0.796 |    0.590 |

## (b) Calibration on held-out rows — label own (slope target 0.8–1.2; obs/exp = matched / ΣP)

| eval_set   | region   |   slope D1_hgb_FULL |   slope D2_hgb_FULL |   slope D1_logit_FULL |   slope D2_logit_FULL |   slope D2_hgb_NO_PNET |   slope step1_lfdr_P1 |   obs_over_exp D1_hgb_FULL |   obs_over_exp D2_hgb_FULL |   obs_over_exp D1_logit_FULL |   obs_over_exp D2_logit_FULL |   obs_over_exp D2_hgb_NO_PNET |   obs_over_exp step1_lfdr_P1 |
|:-----------|:---------|--------------------:|--------------------:|----------------------:|----------------------:|-----------------------:|----------------------:|---------------------------:|---------------------------:|-----------------------------:|-----------------------------:|------------------------------:|-----------------------------:|
| V2010_test | ALL      |                0.89 |                1.01 |                  0.80 |                  0.99 |                   1.02 |                  0.65 |                       1.32 |                       1.00 |                         1.24 |                         1.00 |                          1.01 |                         1.29 |
| V2010_test | NW       |                0.76 |                1.00 |                  0.64 |                  0.99 |                   1.01 |                  0.56 |                       1.20 |                       1.00 |                         1.13 |                         0.99 |                          1.03 |                         1.48 |
| V2010_test | Olympic  |                0.90 |                1.01 |                  0.93 |                  1.01 |                   0.99 |                  0.21 |                       1.49 |                       0.99 |                         1.40 |                         1.01 |                          0.98 |                         1.15 |
| V2010_test | Arm B    |                0.64 |                1.40 |                  0.57 |                  1.09 |                   0.97 |                  0.60 |                       0.75 |                       0.96 |                         0.56 |                         0.88 |                          0.95 |                         1.10 |
| V2008_test | ALL      |                0.76 |                1.03 |                  0.63 |                  1.00 |                   1.05 |                  0.56 |                       1.15 |                       0.99 |                         1.11 |                         1.02 |                          0.99 |                         1.55 |
| V2008_test | NW       |                0.76 |                1.03 |                  0.63 |                  0.99 |                   1.05 |                  0.55 |                       1.15 |                       0.99 |                         1.13 |                         1.01 |                          0.99 |                         1.59 |
| V2008_test | SW       |                0.69 |                1.08 |                  0.70 |                  1.22 |                   1.09 |                  0.68 |                       1.10 |                       1.06 |                         0.86 |                         1.20 |                          0.99 |                         0.92 |
| A2011      | ALL      |                0.64 |                0.93 |                  0.54 |                  0.86 |                   0.86 |                  0.54 |                       0.82 |                       0.76 |                         1.02 |                         2.19 |                          0.79 |                         1.05 |
| A2011      | NW       |                0.64 |                0.93 |                  0.54 |                  0.86 |                   0.86 |                  0.54 |                       0.82 |                       0.76 |                         1.02 |                         2.19 |                          0.79 |                         1.05 |

## (b) Calibration on held-out rows — label any (slope target 0.8–1.2; obs/exp = matched / ΣP)

| eval_set   | region   |   slope D1_hgb_FULL |   slope D2_hgb_FULL |   slope D1_logit_FULL |   slope D2_logit_FULL |   slope D2_hgb_NO_PNET |   slope step1_lfdr_P1 |   obs_over_exp D1_hgb_FULL |   obs_over_exp D2_hgb_FULL |   obs_over_exp D1_logit_FULL |   obs_over_exp D2_logit_FULL |   obs_over_exp D2_hgb_NO_PNET |   obs_over_exp step1_lfdr_P1 |
|:-----------|:---------|--------------------:|--------------------:|----------------------:|----------------------:|-----------------------:|----------------------:|---------------------------:|---------------------------:|-----------------------------:|-----------------------------:|------------------------------:|-----------------------------:|
| V2010_test | ALL      |                1.09 |                1.04 |                  1.02 |                  1.01 |                   1.04 |                  0.54 |                       1.01 |                       1.01 |                         1.01 |                         1.00 |                          1.01 |                         5.08 |
| V2010_test | NW       |                0.90 |                0.95 |                  0.80 |                  0.96 |                   0.95 |                  0.47 |                       1.07 |                       1.04 |                         1.05 |                         1.04 |                          1.08 |                         2.39 |
| V2010_test | Olympic  |                1.15 |                1.06 |                  1.12 |                  1.01 |                   1.04 |                  0.06 |                       1.00 |                       1.00 |                         1.01 |                         0.99 |                          0.99 |                         6.96 |
| V2010_test | Arm B    |                0.87 |                0.74 |                  0.38 |                  0.44 |                   0.82 |                  0.33 |                       0.88 |                       0.96 |                         0.84 |                         0.98 |                          0.96 |                         6.34 |
| V2008_test | ALL      |                0.99 |                1.07 |                  0.87 |                  1.06 |                   1.09 |                  0.53 |                       0.99 |                       0.97 |                         0.95 |                         0.99 |                          0.99 |                         2.25 |
| V2008_test | NW       |                0.99 |                1.06 |                  0.89 |                  1.06 |                   1.09 |                  0.52 |                       0.99 |                       0.97 |                         0.96 |                         0.99 |                          0.98 |                         2.34 |
| V2008_test | SW       |                0.87 |                1.19 |                  1.60 |                  1.27 |                   1.19 |                  0.63 |                       0.96 |                       1.10 |                         0.68 |                         1.01 |                          1.05 |                         1.05 |
| A2011      | ALL      |                0.86 |                0.98 |                  1.25 |                  1.36 |                   0.84 |                  0.50 |                       0.75 |                       0.70 |                         0.58 |                         0.37 |                          0.69 |                         1.90 |
| A2011      | NW       |                0.86 |                0.98 |                  1.25 |                  1.36 |                   0.84 |                  0.50 |                       0.75 |                       0.70 |                         0.58 |                         0.37 |                          0.69 |                         1.90 |

## (b) Calibration on held-out rows — label nbr3 (slope target 0.8–1.2; obs/exp = matched / ΣP)

| eval_set   | region   |   slope D1_hgb_FULL |   slope D2_hgb_FULL |   slope D1_logit_FULL |   slope D2_logit_FULL |   slope D2_hgb_NO_PNET |   slope step1_lfdr_P1 |   obs_over_exp D1_hgb_FULL |   obs_over_exp D2_hgb_FULL |   obs_over_exp D1_logit_FULL |   obs_over_exp D2_logit_FULL |   obs_over_exp D2_hgb_NO_PNET |   obs_over_exp step1_lfdr_P1 |
|:-----------|:---------|--------------------:|--------------------:|----------------------:|----------------------:|-----------------------:|----------------------:|---------------------------:|---------------------------:|-----------------------------:|-----------------------------:|------------------------------:|-----------------------------:|
| V2010_test | ALL      |                0.98 |                1.02 |                  0.92 |                  1.00 |                   1.02 |                  0.58 |                       1.09 |                       1.00 |                         1.07 |                         1.00 |                          1.00 |                         2.72 |
| V2010_test | NW       |                0.79 |                1.01 |                  0.73 |                  0.99 |                   1.03 |                  0.56 |                       1.22 |                       1.02 |                         1.07 |                         1.01 |                          1.06 |                         1.49 |
| V2010_test | Olympic  |                1.06 |                1.02 |                  1.06 |                  1.00 |                   1.00 |                  0.19 |                       1.06 |                       0.99 |                         1.09 |                         1.00 |                          0.98 |                         3.63 |
| V2010_test | Arm B    |                0.67 |              nan    |                  0.63 |                nan    |                   0.71 |                  0.60 |                       0.72 |                       0.90 |                         0.45 |                         0.82 |                          0.90 |                         1.17 |
| V2008_test | ALL      |                0.79 |                1.04 |                  0.73 |                  1.01 |                   1.06 |                  0.57 |                       1.16 |                       1.00 |                         1.05 |                         1.04 |                          1.00 |                         1.56 |
| V2008_test | NW       |                0.79 |                1.03 |                  0.73 |                  1.00 |                   1.06 |                  0.56 |                       1.16 |                       0.99 |                         1.06 |                         1.03 |                          1.00 |                         1.60 |
| V2008_test | SW       |                0.68 |                1.13 |                  0.98 |                  1.22 |                   1.10 |                  0.70 |                       1.13 |                       1.13 |                         0.78 |                         1.16 |                          1.04 |                         0.95 |
| A2011      | ALL      |                0.66 |                0.90 |                  0.81 |                  0.93 |                   0.91 |                  0.55 |                       0.68 |                       0.64 |                         0.89 |                         1.21 |                          0.65 |                         1.06 |
| A2011      | NW       |                0.66 |                0.90 |                  0.81 |                  0.93 |                   0.91 |                  0.55 |                       0.68 |                       0.64 |                         0.89 |                         1.21 |                          0.65 |                         1.06 |

Reliability, V2010_test, ALL, label own:

| bin     |   frac D1_hgb_FULL |   frac D2_hgb_FULL |   frac step1_lfdr_P1 |   mean_P D1_hgb_FULL |   mean_P D2_hgb_FULL |   mean_P step1_lfdr_P1 |   n D1_hgb_FULL |   n D2_hgb_FULL |   n step1_lfdr_P1 |
|:--------|-------------------:|-------------------:|---------------------:|---------------------:|---------------------:|-----------------------:|----------------:|----------------:|------------------:|
| 0.0-0.1 |              0.027 |              0.017 |                0.050 |                0.011 |                0.017 |                  0.019 |      101234.000 |       95656.000 |        101221.000 |
| 0.1-0.2 |              0.216 |              0.145 |                0.076 |                0.141 |                0.143 |                  0.144 |        6136.000 |        7910.000 |          9132.000 |
| 0.2-0.3 |              0.309 |              0.239 |                0.124 |                0.248 |                0.245 |                  0.226 |        2380.000 |        4296.000 |          1930.000 |
| 0.3-0.4 |              0.371 |              0.359 |                0.196 |                0.344 |                0.345 |                  0.344 |        1774.000 |        2423.000 |           438.000 |
| 0.4-0.5 |              0.477 |              0.457 |                0.235 |                0.446 |                0.447 |                  0.452 |         819.000 |        1567.000 |           166.000 |
| 0.5-0.6 |              0.453 |              0.551 |                0.466 |                0.550 |                0.547 |                  0.548 |         585.000 |        1109.000 |           161.000 |
| 0.6-0.7 |              0.568 |              0.634 |                0.466 |                0.653 |                0.647 |                  0.653 |         507.000 |         786.000 |           163.000 |
| 0.7-0.8 |              0.631 |              0.722 |                0.482 |                0.747 |                0.751 |                  0.754 |         464.000 |         630.000 |           168.000 |
| 0.8-0.9 |              0.686 |              0.881 |                0.568 |                0.853 |                0.850 |                  0.856 |         564.000 |         546.000 |           236.000 |
| 0.9-1.0 |              0.809 |              0.939 |                0.730 |                0.959 |                0.931 |                  0.979 |         838.000 |         378.000 |          1686.000 |

Reliability, V2008_test, ALL, label own:

| bin     |   frac D1_hgb_FULL |   frac D2_hgb_FULL |   frac step1_lfdr_P1 |   mean_P D1_hgb_FULL |   mean_P D2_hgb_FULL |   mean_P step1_lfdr_P1 |   n D1_hgb_FULL |   n D2_hgb_FULL |   n step1_lfdr_P1 |
|:--------|-------------------:|-------------------:|---------------------:|---------------------:|---------------------:|-----------------------:|----------------:|----------------:|------------------:|
| 0.0-0.1 |              0.111 |              0.041 |                0.171 |                0.040 |                0.052 |                  0.016 |        6664.000 |        5171.000 |         10328.000 |
| 0.1-0.2 |              0.205 |              0.135 |                0.398 |                0.141 |                0.144 |                  0.153 |        2461.000 |        2451.000 |           924.000 |
| 0.2-0.3 |              0.326 |              0.245 |                0.435 |                0.249 |                0.246 |                  0.230 |        1107.000 |        1460.000 |           269.000 |
| 0.3-0.4 |              0.406 |              0.339 |                0.553 |                0.345 |                0.347 |                  0.350 |         811.000 |         984.000 |           114.000 |
| 0.4-0.5 |              0.473 |              0.473 |                0.537 |                0.448 |                0.448 |                  0.460 |         408.000 |         789.000 |            82.000 |
| 0.5-0.6 |              0.568 |              0.558 |                0.462 |                0.548 |                0.550 |                  0.552 |         347.000 |         767.000 |            93.000 |
| 0.6-0.7 |              0.587 |              0.679 |                0.575 |                0.653 |                0.647 |                  0.648 |         341.000 |         710.000 |           120.000 |
| 0.7-0.8 |              0.639 |              0.767 |                0.537 |                0.750 |                0.749 |                  0.751 |         438.000 |         644.000 |           136.000 |
| 0.8-0.9 |              0.750 |              0.844 |                0.585 |                0.853 |                0.851 |                  0.854 |         501.000 |         583.000 |           188.000 |
| 0.9-1.0 |              0.784 |              0.924 |                0.710 |                0.956 |                0.936 |                  0.981 |         865.000 |         384.000 |          1689.000 |

Reliability, A2011, ALL, label own:

| bin     |   frac D1_hgb_FULL |   frac D2_hgb_FULL |   frac step1_lfdr_P1 |   mean_P D1_hgb_FULL |   mean_P D2_hgb_FULL |   mean_P step1_lfdr_P1 |   n D1_hgb_FULL |   n D2_hgb_FULL |   n step1_lfdr_P1 |
|:--------|-------------------:|-------------------:|---------------------:|---------------------:|---------------------:|-----------------------:|----------------:|----------------:|------------------:|
| 0.0-0.1 |              0.054 |              0.027 |                0.103 |                0.038 |                0.048 |                  0.023 |       36357.000 |       32031.000 |         40454.000 |
| 0.1-0.2 |              0.128 |              0.078 |                0.090 |                0.140 |                0.141 |                  0.151 |       10764.000 |       11781.000 |         12075.000 |
| 0.2-0.3 |              0.195 |              0.161 |                0.130 |                0.246 |                0.245 |                  0.233 |        3969.000 |        5506.000 |          2717.000 |
| 0.3-0.4 |              0.260 |              0.264 |                0.183 |                0.345 |                0.346 |                  0.345 |        2475.000 |        3345.000 |           973.000 |
| 0.4-0.5 |              0.271 |              0.374 |                0.248 |                0.448 |                0.447 |                  0.455 |        1452.000 |        2469.000 |           537.000 |
| 0.5-0.6 |              0.356 |              0.482 |                0.354 |                0.549 |                0.548 |                  0.555 |        1147.000 |        2136.000 |           553.000 |
| 0.6-0.7 |              0.391 |              0.582 |                0.402 |                0.652 |                0.647 |                  0.655 |        1186.000 |        1764.000 |           597.000 |
| 0.7-0.8 |              0.492 |              0.664 |                0.474 |                0.750 |                0.749 |                  0.752 |        1125.000 |        1019.000 |           610.000 |
| 0.8-0.9 |              0.569 |              0.812 |                0.503 |                0.853 |                0.861 |                  0.854 |        1188.000 |        1342.000 |           569.000 |
| 0.9-1.0 |              0.690 |              0.927 |                0.671 |                0.956 |                0.930 |                  0.974 |        2336.000 |         606.000 |          2914.000 |

Reliability, V2010_test, ALL, label any:

| bin     |   frac D1_hgb_FULL |   frac D2_hgb_FULL |   frac step1_lfdr_P1 |   mean_P D1_hgb_FULL |   mean_P D2_hgb_FULL |   mean_P step1_lfdr_P1 |   n D1_hgb_FULL |   n D2_hgb_FULL |   n step1_lfdr_P1 |
|:--------|-------------------:|-------------------:|---------------------:|---------------------:|---------------------:|-----------------------:|----------------:|----------------:|------------------:|
| 0.0-0.1 |              0.057 |              0.066 |                0.254 |                0.072 |                0.074 |                  0.019 |       18418.000 |       19778.000 |        101221.000 |
| 0.1-0.2 |              0.128 |              0.140 |                0.244 |                0.144 |                0.144 |                  0.144 |       39540.000 |       38775.000 |          9132.000 |
| 0.2-0.3 |              0.242 |              0.253 |                0.298 |                0.243 |                0.246 |                  0.226 |       20382.000 |       20694.000 |          1930.000 |
| 0.3-0.4 |              0.365 |              0.349 |                0.356 |                0.347 |                0.345 |                  0.344 |       12776.000 |       12003.000 |           438.000 |
| 0.4-0.5 |              0.489 |              0.467 |                0.452 |                0.448 |                0.448 |                  0.452 |        9147.000 |        8083.000 |           166.000 |
| 0.5-0.6 |              0.599 |              0.553 |                0.596 |                0.549 |                0.547 |                  0.548 |        6449.000 |        6485.000 |           161.000 |
| 0.6-0.7 |              0.707 |              0.665 |                0.693 |                0.643 |                0.646 |                  0.653 |        4285.000 |        4329.000 |           163.000 |
| 0.7-0.8 |              0.792 |              0.762 |                0.643 |                0.745 |                0.743 |                  0.754 |        1940.000 |        2671.000 |           168.000 |
| 0.8-0.9 |              0.828 |              0.856 |                0.708 |                0.844 |                0.847 |                  0.856 |        1125.000 |        1619.000 |           236.000 |
| 0.9-1.0 |              0.885 |              0.954 |                0.835 |                0.962 |                0.942 |                  0.979 |        1239.000 |         864.000 |          1686.000 |

Reliability, A2011, ALL, label any:

| bin     |   frac D1_hgb_FULL |   frac D2_hgb_FULL |   frac step1_lfdr_P1 |   mean_P D1_hgb_FULL |   mean_P D2_hgb_FULL |   mean_P step1_lfdr_P1 |   n D1_hgb_FULL |   n D2_hgb_FULL |   n step1_lfdr_P1 |
|:--------|-------------------:|-------------------:|---------------------:|---------------------:|---------------------:|-----------------------:|----------------:|----------------:|------------------:|
| 0.0-0.1 |              0.057 |              0.072 |                0.231 |                0.089 |                0.089 |                  0.023 |        2013.000 |         844.000 |         40454.000 |
| 0.1-0.2 |              0.094 |              0.077 |                0.186 |                0.150 |                0.155 |                  0.151 |       18280.000 |       16479.000 |         12075.000 |
| 0.2-0.3 |              0.168 |              0.123 |                0.230 |                0.245 |                0.246 |                  0.233 |       13839.000 |       15008.000 |          2717.000 |
| 0.3-0.4 |              0.256 |              0.201 |                0.263 |                0.347 |                0.345 |                  0.345 |        9064.000 |        8513.000 |           973.000 |
| 0.4-0.5 |              0.351 |              0.321 |                0.363 |                0.447 |                0.448 |                  0.455 |        6193.000 |        5630.000 |           537.000 |
| 0.5-0.6 |              0.448 |              0.421 |                0.450 |                0.547 |                0.547 |                  0.555 |        3832.000 |        4425.000 |           553.000 |
| 0.6-0.7 |              0.525 |              0.527 |                0.519 |                0.646 |                0.648 |                  0.655 |        2616.000 |        3289.000 |           597.000 |
| 0.7-0.8 |              0.595 |              0.642 |                0.572 |                0.748 |                0.749 |                  0.752 |        1659.000 |        2892.000 |           610.000 |
| 0.8-0.9 |              0.681 |              0.745 |                0.603 |                0.851 |                0.855 |                  0.854 |        1800.000 |        2826.000 |           569.000 |
| 0.9-1.0 |              0.800 |              0.895 |                0.758 |                0.961 |                0.944 |                  0.974 |        2703.000 |        2093.000 |          2914.000 |

## (c)/(d) Scrambled control and quiet vs episode days — model D2_hgb_FULL_own

FULL models see step-2 attribution, which exists for the first control only → c1; NO_STEP2 uses all 4 controls. `ctl_over_data_sumP` = expected 'real' count the model would assign to one scrambled network ÷ the data's.

| eval_set   | region   | days    | controls   |   n_data |   n_ctl_per_real |   data_fracP_gt05 |   ctl_fracP_gt05 |   data_meanP |   ctl_meanP |   data_sumP |   ctl_sumP_per_real |   ctl_over_data_sumP |
|:-----------|:---------|:--------|:-----------|---------:|-----------------:|------------------:|-----------------:|-------------:|------------:|------------:|--------------------:|---------------------:|
| V2010_test | ALL      | all     | c1         |   172854 |       169413.000 |             0.022 |            0.007 |        0.060 |       0.047 |   10354.561 |            8015.337 |                0.774 |
| V2010_test | ALL      | quiet   | c1         |    48576 |        49192.000 |             0.003 |            0.001 |        0.039 |       0.036 |    1875.827 |            1764.190 |                0.940 |
| V2010_test | ALL      | episode | c1         |   124278 |       120221.000 |             0.030 |            0.009 |        0.068 |       0.052 |    8478.734 |            6251.148 |                0.737 |
| V2010_test | NW       | all     | c1         |    21800 |        18954.000 |             0.136 |            0.037 |        0.198 |       0.125 |    4320.688 |            2369.402 |                0.548 |
| V2010_test | NW       | quiet   | c1         |     5429 |         5599.000 |             0.001 |            0.001 |        0.060 |       0.058 |     326.191 |             323.200 |                0.991 |
| V2010_test | NW       | episode | c1         |    16371 |        13355.000 |             0.181 |            0.053 |        0.244 |       0.153 |    3994.497 |            2046.202 |                0.512 |
| V2010_test | SW       | all     | c1         |      381 |          169.000 |             0.462 |            0.006 |        0.432 |       0.086 |     164.408 |              14.450 |                0.088 |
| V2010_test | SW       | quiet   | c1         |       33 |           44.000 |             0.030 |            0.000 |        0.070 |       0.037 |       2.297 |               1.640 |                0.714 |
| V2010_test | SW       | episode | c1         |      348 |          125.000 |             0.503 |            0.008 |        0.466 |       0.102 |     162.111 |              12.810 |                0.079 |
| V2010_test | Olympic  | all     | c1         |   148866 |       148422.000 |             0.004 |            0.003 |        0.039 |       0.037 |    5769.231 |            5550.715 |                0.962 |
| V2010_test | Olympic  | quiet   | c1         |    42747 |        43134.000 |             0.003 |            0.002 |        0.036 |       0.033 |    1536.869 |            1427.508 |                0.929 |
| V2010_test | Olympic  | episode | c1         |   106119 |       105288.000 |             0.005 |            0.004 |        0.040 |       0.039 |    4232.362 |            4123.207 |                0.974 |
| V2010_test | Arm B    | all     | c1         |     1807 |         1868.000 |             0.011 |            0.000 |        0.055 |       0.043 |     100.233 |              80.770 |                0.806 |
| V2010_test | Arm B    | quiet   | c1         |      367 |          415.000 |             0.000 |            0.000 |        0.029 |       0.029 |      10.470 |              11.842 |                1.131 |
| V2010_test | Arm B    | episode | c1         |     1440 |         1453.000 |             0.014 |            0.000 |        0.062 |       0.047 |      89.764 |              68.929 |                0.768 |
| V2008_test | ALL      | all     | c1         |    93754 |        90700.000 |             0.035 |            0.009 |        0.078 |       0.058 |    7286.536 |            5269.066 |                0.723 |
| V2008_test | ALL      | quiet   | c1         |    13913 |        13285.000 |             0.002 |            0.002 |        0.056 |       0.058 |     785.839 |             776.047 |                0.988 |
| V2008_test | ALL      | episode | c1         |    79841 |        77415.000 |             0.040 |            0.010 |        0.081 |       0.058 |    6500.697 |            4493.019 |                0.691 |
| V2008_test | NW       | all     | c1         |    16322 |        13531.000 |             0.181 |            0.050 |        0.244 |       0.153 |    3977.435 |            2075.460 |                0.522 |
| V2008_test | NW       | quiet   | c1         |     1993 |         1980.000 |             0.005 |            0.004 |        0.076 |       0.076 |     151.135 |             151.054 |                0.999 |
| V2008_test | NW       | episode | c1         |    14329 |        11551.000 |             0.206 |            0.057 |        0.267 |       0.167 |    3826.300 |            1924.406 |                0.503 |
| V2008_test | SW       | all     | c1         |      305 |          153.000 |             0.472 |            0.026 |        0.461 |       0.133 |     140.654 |              20.324 |                0.144 |
| V2008_test | SW       | quiet   | c1         |        9 |           11.000 |             0.000 |            0.000 |        0.049 |       0.030 |       0.439 |               0.327 |                0.746 |
| V2008_test | SW       | episode | c1         |      296 |          142.000 |             0.486 |            0.028 |        0.474 |       0.141 |     140.215 |              19.996 |                0.143 |
| V2008_test | Olympic  | all     | c1         |    73842 |        73854.000 |             0.002 |            0.001 |        0.041 |       0.042 |    3056.928 |            3072.979 |                1.005 |
| V2008_test | Olympic  | quiet   | c1         |    11474 |        10883.000 |             0.001 |            0.001 |        0.054 |       0.056 |     623.570 |             614.590 |                0.986 |
| V2008_test | Olympic  | episode | c1         |    62368 |        62971.000 |             0.002 |            0.001 |        0.039 |       0.039 |    2433.358 |            2458.388 |                1.010 |
| V2008_test | Arm B    | all     | c1         |     3285 |         3162.000 |             0.001 |            0.000 |        0.034 |       0.032 |     111.520 |             100.304 |                0.899 |
| V2008_test | Arm B    | quiet   | c1         |      437 |          411.000 |             0.000 |            0.000 |        0.024 |       0.025 |      10.696 |              10.075 |                0.942 |
| V2008_test | Arm B    | episode | c1         |     2848 |         2751.000 |             0.001 |            0.000 |        0.035 |       0.033 |     100.824 |              90.229 |                0.895 |
| A2011      | ALL      | all     | c1         |   294824 |       287560.000 |             0.024 |            0.007 |        0.082 |       0.069 |   24122.506 |           19790.394 |                0.820 |
| A2011      | ALL      | episode | c1         |   294824 |       287560.000 |             0.024 |            0.007 |        0.082 |       0.069 |   24122.506 |           19790.394 |                0.820 |
| A2011      | NW       | all     | c1         |    61999 |        55860.000 |             0.111 |            0.035 |        0.187 |       0.133 |   11619.814 |            7414.632 |                0.638 |
| A2011      | NW       | episode | c1         |    61999 |        55860.000 |             0.111 |            0.035 |        0.187 |       0.133 |   11619.814 |            7414.632 |                0.638 |
| A2011      | SW       | all     | c1         |      709 |          714.000 |             0.003 |            0.001 |        0.143 |       0.146 |     101.562 |             104.257 |                1.027 |
| A2011      | SW       | episode | c1         |      709 |          714.000 |             0.003 |            0.001 |        0.143 |       0.146 |     101.562 |             104.257 |                1.027 |
| A2011      | Olympic  | all     | c1         |   215645 |       214694.000 |             0.000 |            0.000 |        0.055 |       0.054 |   11786.682 |           11695.536 |                0.992 |
| A2011      | Olympic  | episode | c1         |   215645 |       214694.000 |             0.000 |            0.000 |        0.055 |       0.054 |   11786.682 |           11695.536 |                0.992 |
| A2011      | Arm B    | all     | c1         |    16471 |        16292.000 |             0.001 |            0.000 |        0.037 |       0.035 |     614.449 |             575.968 |                0.937 |
| A2011      | Arm B    | episode | c1         |    16471 |        16292.000 |             0.001 |            0.000 |        0.037 |       0.035 |     614.449 |             575.968 |                0.937 |

## (c)/(d) Scrambled control and quiet vs episode days — model D2_hgb_FULL_any

FULL models see step-2 attribution, which exists for the first control only → c1; NO_STEP2 uses all 4 controls. `ctl_over_data_sumP` = expected 'real' count the model would assign to one scrambled network ÷ the data's.

| eval_set   | region   | days    | controls   |   n_data |   n_ctl_per_real |   data_fracP_gt05 |   ctl_fracP_gt05 |   data_meanP |   ctl_meanP |   data_sumP |   ctl_sumP_per_real |   ctl_over_data_sumP |
|:-----------|:---------|:--------|:-----------|---------:|-----------------:|------------------:|-----------------:|-------------:|------------:|------------:|--------------------:|---------------------:|
| V2010_test | ALL      | all     | c1         |   172854 |       169413.000 |             0.122 |            0.102 |        0.253 |       0.241 |   43701.209 |           40767.668 |                0.933 |
| V2010_test | ALL      | quiet   | c1         |    48576 |        49192.000 |             0.086 |            0.075 |        0.228 |       0.221 |   11098.350 |           10891.637 |                0.981 |
| V2010_test | ALL      | episode | c1         |   124278 |       120221.000 |             0.136 |            0.114 |        0.262 |       0.249 |   32602.859 |           29876.031 |                0.916 |
| V2010_test | NW       | all     | c1         |    21800 |        18954.000 |             0.235 |            0.133 |        0.346 |       0.279 |    7540.746 |            5295.395 |                0.702 |
| V2010_test | NW       | quiet   | c1         |     5429 |         5599.000 |             0.020 |            0.014 |        0.190 |       0.183 |    1029.915 |            1025.539 |                0.996 |
| V2010_test | NW       | episode | c1         |    16371 |        13355.000 |             0.306 |            0.183 |        0.398 |       0.320 |    6510.831 |            4269.857 |                0.656 |
| V2010_test | SW       | all     | c1         |      381 |          169.000 |             0.512 |            0.053 |        0.508 |       0.210 |     193.672 |              35.562 |                0.184 |
| V2010_test | SW       | quiet   | c1         |       33 |           44.000 |             0.061 |            0.000 |        0.184 |       0.142 |       6.067 |               6.231 |                1.027 |
| V2010_test | SW       | episode | c1         |      348 |          125.000 |             0.555 |            0.072 |        0.539 |       0.235 |     187.606 |              29.331 |                0.156 |
| V2010_test | Olympic  | all     | c1         |   148866 |       148422.000 |             0.103 |            0.098 |        0.237 |       0.234 |   35324.831 |           34790.547 |                0.985 |
| V2010_test | Olympic  | quiet   | c1         |    42747 |        43134.000 |             0.094 |            0.082 |        0.233 |       0.225 |    9945.503 |            9726.255 |                0.978 |
| V2010_test | Olympic  | episode | c1         |   106119 |       105288.000 |             0.107 |            0.104 |        0.239 |       0.238 |   25379.327 |           25064.292 |                0.988 |
| V2010_test | Arm B    | all     | c1         |     1807 |         1868.000 |             0.191 |            0.172 |        0.355 |       0.346 |     641.960 |             646.164 |                1.007 |
| V2010_test | Arm B    | quiet   | c1         |      367 |          415.000 |             0.139 |            0.128 |        0.318 |       0.322 |     116.866 |             133.612 |                1.143 |
| V2010_test | Arm B    | episode | c1         |     1440 |         1453.000 |             0.205 |            0.185 |        0.365 |       0.353 |     525.094 |             512.552 |                0.976 |
| V2008_test | ALL      | all     | c1         |    93754 |        90700.000 |             0.177 |            0.152 |        0.302 |       0.289 |   28294.442 |           26194.314 |                0.926 |
| V2008_test | ALL      | quiet   | c1         |    13913 |        13285.000 |             0.223 |            0.238 |        0.325 |       0.334 |    4526.172 |            4436.084 |                0.980 |
| V2008_test | ALL      | episode | c1         |    79841 |        77415.000 |             0.169 |            0.138 |        0.298 |       0.281 |   23768.269 |           21758.231 |                0.915 |
| V2008_test | NW       | all     | c1         |    16322 |        13531.000 |             0.294 |            0.159 |        0.378 |       0.294 |    6168.457 |            3981.615 |                0.645 |
| V2008_test | NW       | quiet   | c1         |     1993 |         1980.000 |             0.032 |            0.028 |        0.202 |       0.206 |     403.474 |             408.858 |                1.013 |
| V2008_test | NW       | episode | c1         |    14329 |        11551.000 |             0.331 |            0.181 |        0.402 |       0.309 |    5764.984 |            3572.756 |                0.620 |
| V2008_test | SW       | all     | c1         |      305 |          153.000 |             0.502 |            0.065 |        0.510 |       0.234 |     155.634 |              35.787 |                0.230 |
| V2008_test | SW       | quiet   | c1         |        9 |           11.000 |             0.000 |            0.000 |        0.156 |       0.127 |       1.403 |               1.401 |                0.998 |
| V2008_test | SW       | episode | c1         |      296 |          142.000 |             0.517 |            0.070 |        0.521 |       0.242 |     154.231 |              34.387 |                0.223 |
| V2008_test | Olympic  | all     | c1         |    73842 |        73854.000 |             0.152 |            0.153 |        0.283 |       0.286 |   20871.565 |           21113.209 |                1.012 |
| V2008_test | Olympic  | quiet   | c1         |    11474 |        10883.000 |             0.262 |            0.282 |        0.347 |       0.358 |    3982.374 |            3893.757 |                0.978 |
| V2008_test | Olympic  | episode | c1         |    62368 |        62971.000 |             0.132 |            0.130 |        0.271 |       0.273 |   16889.190 |           17219.452 |                1.020 |
| V2008_test | Arm B    | all     | c1         |     3285 |         3162.000 |             0.126 |            0.118 |        0.334 |       0.336 |    1098.786 |            1063.704 |                0.968 |
| V2008_test | Arm B    | quiet   | c1         |      437 |          411.000 |             0.087 |            0.095 |        0.318 |       0.321 |     138.921 |             132.068 |                0.951 |
| V2008_test | Arm B    | episode | c1         |     2848 |         2751.000 |             0.132 |            0.121 |        0.337 |       0.339 |     959.864 |             931.636 |                0.971 |
| A2011      | ALL      | all     | c1         |   294824 |       287560.000 |             0.187 |            0.169 |        0.319 |       0.309 |   94189.558 |           88866.601 |                0.943 |
| A2011      | ALL      | episode | c1         |   294824 |       287560.000 |             0.187 |            0.169 |        0.319 |       0.309 |   94189.558 |           88866.601 |                0.943 |
| A2011      | NW       | all     | c1         |    61999 |        55860.000 |             0.250 |            0.170 |        0.369 |       0.320 |   22891.151 |           17862.052 |                0.780 |
| A2011      | NW       | episode | c1         |    61999 |        55860.000 |             0.250 |            0.170 |        0.369 |       0.320 |   22891.151 |           17862.052 |                0.780 |
| A2011      | SW       | all     | c1         |      709 |          714.000 |             0.048 |            0.055 |        0.296 |       0.302 |     209.751 |             215.923 |                1.029 |
| A2011      | SW       | episode | c1         |      709 |          714.000 |             0.048 |            0.055 |        0.296 |       0.302 |     209.751 |             215.923 |                1.029 |
| A2011      | Olympic  | all     | c1         |   215645 |       214694.000 |             0.177 |            0.177 |        0.308 |       0.309 |   66517.978 |           66298.330 |                0.997 |
| A2011      | Olympic  | episode | c1         |   215645 |       214694.000 |             0.177 |            0.177 |        0.308 |       0.309 |   66517.978 |           66298.330 |                0.997 |
| A2011      | Arm B    | all     | c1         |    16471 |        16292.000 |             0.074 |            0.069 |        0.277 |       0.276 |    4570.677 |            4490.297 |                0.982 |
| A2011      | Arm B    | episode | c1         |    16471 |        16292.000 |             0.074 |            0.069 |        0.277 |       0.276 |    4570.677 |            4490.297 |                0.982 |

## (c)/(d) Scrambled control and quiet vs episode days — model D2_hgb_NO_STEP2_own

FULL models see step-2 attribution, which exists for the first control only → c1; NO_STEP2 uses all 4 controls. `ctl_over_data_sumP` = expected 'real' count the model would assign to one scrambled network ÷ the data's.

| eval_set   | region   | days    | controls    |   n_data |   n_ctl_per_real |   data_fracP_gt05 |   ctl_fracP_gt05 |   data_meanP |   ctl_meanP |   data_sumP |   ctl_sumP_per_real |   ctl_over_data_sumP |
|:-----------|:---------|:--------|:------------|---------:|-----------------:|------------------:|-----------------:|-------------:|------------:|------------:|--------------------:|---------------------:|
| V2010_test | ALL      | all     | c1,c2,c3,c4 |   172854 |       168969.250 |             0.022 |            0.007 |        0.060 |       0.047 |   10398.062 |            8023.406 |                0.772 |
| V2010_test | ALL      | quiet   | c1,c2,c3,c4 |    48576 |        48451.250 |             0.003 |            0.001 |        0.039 |       0.037 |    1905.313 |            1789.020 |                0.939 |
| V2010_test | ALL      | episode | c1,c2,c3,c4 |   124278 |       120518.000 |             0.029 |            0.009 |        0.068 |       0.052 |    8492.750 |            6234.386 |                0.734 |
| V2010_test | NW       | all     | c1,c2,c3,c4 |    21800 |        19002.500 |             0.138 |            0.036 |        0.198 |       0.125 |    4327.169 |            2378.183 |                0.550 |
| V2010_test | NW       | quiet   | c1,c2,c3,c4 |     5429 |         5521.000 |             0.001 |            0.001 |        0.060 |       0.059 |     326.029 |             324.804 |                0.996 |
| V2010_test | NW       | episode | c1,c2,c3,c4 |    16371 |        13481.500 |             0.184 |            0.051 |        0.244 |       0.152 |    4001.140 |            2053.379 |                0.513 |
| V2010_test | SW       | all     | c1,c2,c3,c4 |      381 |          167.750 |             0.465 |            0.028 |        0.433 |       0.094 |     164.989 |              15.851 |                0.096 |
| V2010_test | SW       | quiet   | c1,c2,c3,c4 |       33 |           44.250 |             0.030 |            0.000 |        0.071 |       0.032 |       2.339 |               1.417 |                0.606 |
| V2010_test | SW       | episode | c1,c2,c3,c4 |      348 |          123.500 |             0.506 |            0.038 |        0.467 |       0.117 |     162.651 |              14.433 |                0.089 |
| V2010_test | Olympic  | all     | c1,c2,c3,c4 |   148866 |       147917.250 |             0.004 |            0.003 |        0.039 |       0.038 |    5805.939 |            5550.447 |                0.956 |
| V2010_test | Olympic  | quiet   | c1,c2,c3,c4 |    42747 |        42469.000 |             0.003 |            0.001 |        0.037 |       0.034 |    1566.208 |            1450.088 |                0.926 |
| V2010_test | Olympic  | episode | c1,c2,c3,c4 |   106119 |       105448.250 |             0.004 |            0.003 |        0.040 |       0.039 |    4239.731 |            4100.358 |                0.967 |
| V2010_test | Arm B    | all     | c1,c2,c3,c4 |     1807 |         1881.750 |             0.011 |            0.000 |        0.055 |       0.042 |      99.965 |              78.925 |                0.790 |
| V2010_test | Arm B    | quiet   | c1,c2,c3,c4 |      367 |          417.000 |             0.000 |            0.001 |        0.029 |       0.030 |      10.738 |              12.710 |                1.184 |
| V2010_test | Arm B    | episode | c1,c2,c3,c4 |     1440 |         1464.750 |             0.013 |            0.000 |        0.062 |       0.045 |      89.227 |              66.216 |                0.742 |
| V2008_test | ALL      | all     | c1,c2,c3,c4 |    93754 |        90829.500 |             0.034 |            0.008 |        0.076 |       0.057 |    7165.579 |            5181.509 |                0.723 |
| V2008_test | ALL      | quiet   | c1,c2,c3,c4 |    13913 |        13121.500 |             0.001 |            0.001 |        0.057 |       0.058 |     794.678 |             767.564 |                0.966 |
| V2008_test | ALL      | episode | c1,c2,c3,c4 |    79841 |        77708.000 |             0.039 |            0.009 |        0.080 |       0.057 |    6370.901 |            4413.945 |                0.693 |
| V2008_test | NW       | all     | c1,c2,c3,c4 |    16322 |        13546.250 |             0.181 |            0.053 |        0.245 |       0.156 |    3995.131 |            2109.227 |                0.528 |
| V2008_test | NW       | quiet   | c1,c2,c3,c4 |     1993 |         1966.000 |             0.005 |            0.004 |        0.078 |       0.076 |     154.634 |             150.027 |                0.970 |
| V2008_test | NW       | episode | c1,c2,c3,c4 |    14329 |        11580.250 |             0.205 |            0.062 |        0.268 |       0.169 |    3840.497 |            1959.200 |                0.510 |
| V2008_test | SW       | all     | c1,c2,c3,c4 |      305 |          145.250 |             0.485 |            0.029 |        0.467 |       0.131 |     142.482 |              18.965 |                0.133 |
| V2008_test | SW       | quiet   | c1,c2,c3,c4 |        9 |           10.750 |             0.000 |            0.000 |        0.056 |       0.047 |       0.503 |               0.502 |                0.999 |
| V2008_test | SW       | episode | c1,c2,c3,c4 |      296 |          134.500 |             0.500 |            0.032 |        0.480 |       0.137 |     141.979 |              18.463 |                0.130 |
| V2008_test | Olympic  | all     | c1,c2,c3,c4 |    73842 |        73895.750 |             0.001 |            0.000 |        0.040 |       0.040 |    2929.304 |            2963.217 |                1.012 |
| V2008_test | Olympic  | quiet   | c1,c2,c3,c4 |    11474 |        10724.500 |             0.000 |            0.000 |        0.055 |       0.057 |     629.986 |             607.829 |                0.965 |
| V2008_test | Olympic  | episode | c1,c2,c3,c4 |    62368 |        63171.250 |             0.001 |            0.000 |        0.037 |       0.037 |    2299.318 |            2355.388 |                1.024 |
| V2008_test | Arm B    | all     | c1,c2,c3,c4 |     3285 |         3242.250 |             0.001 |            0.000 |        0.030 |       0.028 |      98.662 |              90.100 |                0.913 |
| V2008_test | Arm B    | quiet   | c1,c2,c3,c4 |      437 |          420.250 |             0.000 |            0.000 |        0.022 |       0.022 |       9.555 |               9.206 |                0.963 |
| V2008_test | Arm B    | episode | c1,c2,c3,c4 |     2848 |         2822.000 |             0.001 |            0.000 |        0.031 |       0.029 |      89.107 |              80.895 |                0.908 |
| A2011      | ALL      | all     | c1,c2,c3,c4 |   294824 |       287372.500 |             0.024 |            0.007 |        0.083 |       0.070 |   24477.291 |           20182.633 |                0.825 |
| A2011      | ALL      | episode | c1,c2,c3,c4 |   294824 |       287372.500 |             0.024 |            0.007 |        0.083 |       0.070 |   24477.291 |           20182.633 |                0.825 |
| A2011      | NW       | all     | c1,c2,c3,c4 |    61999 |        55640.250 |             0.112 |            0.038 |        0.189 |       0.136 |   11724.463 |            7543.417 |                0.643 |
| A2011      | NW       | episode | c1,c2,c3,c4 |    61999 |        55640.250 |             0.112 |            0.038 |        0.189 |       0.136 |   11724.463 |            7543.417 |                0.643 |
| A2011      | SW       | all     | c1,c2,c3,c4 |      709 |          699.750 |             0.004 |            0.003 |        0.155 |       0.157 |     109.640 |             110.081 |                1.004 |
| A2011      | SW       | episode | c1,c2,c3,c4 |      709 |          699.750 |             0.004 |            0.003 |        0.155 |       0.157 |     109.640 |             110.081 |                1.004 |
| A2011      | Olympic  | all     | c1,c2,c3,c4 |   215645 |       214685.000 |             0.000 |            0.000 |        0.056 |       0.056 |   12123.394 |           12053.432 |                0.994 |
| A2011      | Olympic  | episode | c1,c2,c3,c4 |   215645 |       214685.000 |             0.000 |            0.000 |        0.056 |       0.056 |   12123.394 |           12053.432 |                0.994 |
| A2011      | Arm B    | all     | c1,c2,c3,c4 |    16471 |        16347.500 |             0.001 |            0.000 |        0.032 |       0.029 |     519.794 |             475.704 |                0.915 |
| A2011      | Arm B    | episode | c1,c2,c3,c4 |    16471 |        16347.500 |             0.001 |            0.000 |        0.032 |       0.029 |     519.794 |             475.704 |                0.915 |

## (c)/(d) Scrambled control and quiet vs episode days — model D2_hgb_NO_STEP2_any

FULL models see step-2 attribution, which exists for the first control only → c1; NO_STEP2 uses all 4 controls. `ctl_over_data_sumP` = expected 'real' count the model would assign to one scrambled network ÷ the data's.

| eval_set   | region   | days    | controls    |   n_data |   n_ctl_per_real |   data_fracP_gt05 |   ctl_fracP_gt05 |   data_meanP |   ctl_meanP |   data_sumP |   ctl_sumP_per_real |   ctl_over_data_sumP |
|:-----------|:---------|:--------|:------------|---------:|-----------------:|------------------:|-----------------:|-------------:|------------:|------------:|--------------------:|---------------------:|
| V2010_test | ALL      | all     | c1,c2,c3,c4 |   172854 |       168969.250 |             0.116 |            0.098 |        0.255 |       0.243 |   43998.926 |           41010.432 |                0.932 |
| V2010_test | ALL      | quiet   | c1,c2,c3,c4 |    48576 |        48451.250 |             0.088 |            0.081 |        0.235 |       0.230 |   11422.783 |           11131.014 |                0.974 |
| V2010_test | ALL      | episode | c1,c2,c3,c4 |   124278 |       120518.000 |             0.126 |            0.104 |        0.262 |       0.248 |   32576.143 |           29879.417 |                0.917 |
| V2010_test | NW       | all     | c1,c2,c3,c4 |    21800 |        19002.500 |             0.235 |            0.136 |        0.347 |       0.283 |    7575.079 |            5371.523 |                0.709 |
| V2010_test | NW       | quiet   | c1,c2,c3,c4 |     5429 |         5521.000 |             0.021 |            0.018 |        0.190 |       0.186 |    1032.848 |            1025.206 |                0.993 |
| V2010_test | NW       | episode | c1,c2,c3,c4 |    16371 |        13481.500 |             0.306 |            0.184 |        0.400 |       0.322 |    6542.231 |            4346.318 |                0.664 |
| V2010_test | SW       | all     | c1,c2,c3,c4 |      381 |          167.750 |             0.501 |            0.066 |        0.508 |       0.214 |     193.439 |              35.974 |                0.186 |
| V2010_test | SW       | quiet   | c1,c2,c3,c4 |       33 |           44.250 |             0.061 |            0.000 |        0.180 |       0.131 |       5.953 |               5.812 |                0.976 |
| V2010_test | SW       | episode | c1,c2,c3,c4 |      348 |          123.500 |             0.543 |            0.089 |        0.539 |       0.244 |     187.486 |              30.162 |                0.161 |
| V2010_test | Olympic  | all     | c1,c2,c3,c4 |   148866 |       147917.250 |             0.096 |            0.092 |        0.239 |       0.236 |   35588.793 |           34965.119 |                0.982 |
| V2010_test | Olympic  | quiet   | c1,c2,c3,c4 |    42747 |        42469.000 |             0.096 |            0.088 |        0.240 |       0.235 |   10265.388 |            9967.053 |                0.971 |
| V2010_test | Olympic  | episode | c1,c2,c3,c4 |   106119 |       105448.250 |             0.096 |            0.093 |        0.239 |       0.237 |   25323.405 |           24998.067 |                0.987 |
| V2010_test | Arm B    | all     | c1,c2,c3,c4 |     1807 |         1881.750 |             0.193 |            0.160 |        0.355 |       0.339 |     641.615 |             637.814 |                0.994 |
| V2010_test | Arm B    | quiet   | c1,c2,c3,c4 |      367 |          417.000 |             0.163 |            0.138 |        0.323 |       0.319 |     118.594 |             132.944 |                1.121 |
| V2010_test | Arm B    | episode | c1,c2,c3,c4 |     1440 |         1464.750 |             0.200 |            0.166 |        0.363 |       0.345 |     523.020 |             504.871 |                0.965 |
| V2008_test | ALL      | all     | c1,c2,c3,c4 |    93754 |        90829.500 |             0.160 |            0.139 |        0.296 |       0.285 |   27773.352 |           25872.991 |                0.932 |
| V2008_test | ALL      | quiet   | c1,c2,c3,c4 |    13913 |        13121.500 |             0.191 |            0.209 |        0.312 |       0.322 |    4347.161 |            4228.738 |                0.973 |
| V2008_test | ALL      | episode | c1,c2,c3,c4 |    79841 |        77708.000 |             0.155 |            0.127 |        0.293 |       0.279 |   23426.191 |           21644.253 |                0.924 |
| V2008_test | NW       | all     | c1,c2,c3,c4 |    16322 |        13546.250 |             0.291 |            0.160 |        0.378 |       0.295 |    6168.510 |            3995.246 |                0.648 |
| V2008_test | NW       | quiet   | c1,c2,c3,c4 |     1993 |         1966.000 |             0.031 |            0.025 |        0.201 |       0.204 |     400.581 |             400.530 |                1.000 |
| V2008_test | NW       | episode | c1,c2,c3,c4 |    14329 |        11580.250 |             0.328 |            0.183 |        0.403 |       0.310 |    5767.929 |            3594.716 |                0.623 |
| V2008_test | SW       | all     | c1,c2,c3,c4 |      305 |          145.250 |             0.498 |            0.062 |        0.511 |       0.235 |     155.862 |              34.150 |                0.219 |
| V2008_test | SW       | quiet   | c1,c2,c3,c4 |        9 |           10.750 |             0.000 |            0.000 |        0.154 |       0.141 |       1.384 |               1.513 |                1.093 |
| V2008_test | SW       | episode | c1,c2,c3,c4 |      296 |          134.500 |             0.514 |            0.067 |        0.522 |       0.243 |     154.477 |              32.637 |                0.211 |
| V2008_test | Olympic  | all     | c1,c2,c3,c4 |    73842 |        73895.750 |             0.133 |            0.137 |        0.276 |       0.282 |   20391.608 |           20801.734 |                1.020 |
| V2008_test | Olympic  | quiet   | c1,c2,c3,c4 |    11474 |        10724.500 |             0.224 |            0.249 |        0.332 |       0.344 |    3808.152 |            3694.356 |                0.970 |
| V2008_test | Olympic  | episode | c1,c2,c3,c4 |    62368 |        63171.250 |             0.116 |            0.118 |        0.266 |       0.271 |   16583.456 |           17107.378 |                1.032 |
| V2008_test | Arm B    | all     | c1,c2,c3,c4 |     3285 |         3242.250 |             0.096 |            0.096 |        0.322 |       0.321 |    1057.372 |            1041.861 |                0.985 |
| V2008_test | Arm B    | quiet   | c1,c2,c3,c4 |      437 |          420.250 |             0.064 |            0.072 |        0.314 |       0.315 |     137.042 |             132.339 |                0.966 |
| V2008_test | Arm B    | episode | c1,c2,c3,c4 |     2848 |         2822.000 |             0.101 |            0.099 |        0.323 |       0.322 |     920.330 |             909.522 |                0.988 |
| A2011      | ALL      | all     | c1,c2,c3,c4 |   294824 |       287372.500 |             0.194 |            0.178 |        0.320 |       0.310 |   94422.941 |           89158.516 |                0.944 |
| A2011      | ALL      | episode | c1,c2,c3,c4 |   294824 |       287372.500 |             0.194 |            0.178 |        0.320 |       0.310 |   94422.941 |           89158.516 |                0.944 |
| A2011      | NW       | all     | c1,c2,c3,c4 |    61999 |        55640.250 |             0.239 |            0.160 |        0.362 |       0.313 |   22425.125 |           17411.415 |                0.776 |
| A2011      | NW       | episode | c1,c2,c3,c4 |    61999 |        55640.250 |             0.239 |            0.160 |        0.362 |       0.313 |   22425.125 |           17411.415 |                0.776 |
| A2011      | SW       | all     | c1,c2,c3,c4 |      709 |          699.750 |             0.039 |            0.041 |        0.287 |       0.293 |     203.498 |             204.749 |                1.006 |
| A2011      | SW       | episode | c1,c2,c3,c4 |      709 |          699.750 |             0.039 |            0.041 |        0.287 |       0.293 |     203.498 |             204.749 |                1.006 |
| A2011      | Olympic  | all     | c1,c2,c3,c4 |   215645 |       214685.000 |             0.193 |            0.193 |        0.313 |       0.314 |   67502.374 |           67316.518 |                0.997 |
| A2011      | Olympic  | episode | c1,c2,c3,c4 |   215645 |       214685.000 |             0.193 |            0.193 |        0.313 |       0.314 |   67502.374 |           67316.518 |                0.997 |
| A2011      | Arm B    | all     | c1,c2,c3,c4 |    16471 |        16347.500 |             0.052 |            0.046 |        0.261 |       0.259 |    4291.944 |            4225.834 |                0.985 |
| A2011      | Arm B    | episode | c1,c2,c3,c4 |    16471 |        16347.500 |             0.052 |            0.046 |        0.261 |       0.259 |    4291.944 |            4225.834 |                0.985 |

## (c)/(d) Scrambled control and quiet vs episode days — model D1_hgb_FULL_own

FULL models see step-2 attribution, which exists for the first control only → c1; NO_STEP2 uses all 4 controls. `ctl_over_data_sumP` = expected 'real' count the model would assign to one scrambled network ÷ the data's.

| eval_set   | region   | days    | controls   |   n_data |   n_ctl_per_real |   data_fracP_gt05 |   ctl_fracP_gt05 |   data_meanP |   ctl_meanP |   data_sumP |   ctl_sumP_per_real |   ctl_over_data_sumP |
|:-----------|:---------|:--------|:-----------|---------:|-----------------:|------------------:|-----------------:|-------------:|------------:|------------:|--------------------:|---------------------:|
| V2010_test | ALL      | all     | c1         |   172854 |       169413.000 |             0.021 |            0.005 |        0.046 |       0.032 |    8017.106 |            5417.364 |                0.676 |
| V2010_test | ALL      | quiet   | c1         |    48576 |        49192.000 |             0.009 |            0.004 |        0.032 |       0.027 |    1534.697 |            1313.437 |                0.856 |
| V2010_test | ALL      | episode | c1         |   124278 |       120221.000 |             0.026 |            0.006 |        0.052 |       0.034 |    6482.409 |            4103.927 |                0.633 |
| V2010_test | NW       | all     | c1         |    21800 |        18954.000 |             0.111 |            0.018 |        0.171 |       0.093 |    3718.813 |            1764.326 |                0.474 |
| V2010_test | NW       | quiet   | c1         |     5429 |         5599.000 |             0.006 |            0.006 |        0.062 |       0.060 |     334.985 |             336.237 |                1.004 |
| V2010_test | NW       | episode | c1         |    16371 |        13355.000 |             0.145 |            0.024 |        0.207 |       0.107 |    3383.828 |            1428.089 |                0.422 |
| V2010_test | SW       | all     | c1         |      381 |          169.000 |             0.318 |            0.047 |        0.345 |       0.142 |     131.479 |              23.997 |                0.183 |
| V2010_test | SW       | quiet   | c1         |       33 |           44.000 |             0.091 |            0.068 |        0.169 |       0.148 |       5.578 |               6.499 |                1.165 |
| V2010_test | SW       | episode | c1         |      348 |          125.000 |             0.339 |            0.040 |        0.362 |       0.140 |     125.901 |              17.498 |                0.139 |
| V2010_test | Olympic  | all     | c1         |   148866 |       148422.000 |             0.007 |            0.004 |        0.027 |       0.024 |    4038.412 |            3544.553 |                0.878 |
| V2010_test | Olympic  | quiet   | c1         |    42747 |        43134.000 |             0.009 |            0.004 |        0.028 |       0.022 |    1180.958 |             958.626 |                0.812 |
| V2010_test | Olympic  | episode | c1         |   106119 |       105288.000 |             0.006 |            0.004 |        0.027 |       0.025 |    2857.454 |            2585.927 |                0.905 |
| V2010_test | Arm B    | all     | c1         |     1807 |         1868.000 |             0.036 |            0.010 |        0.071 |       0.045 |     128.402 |              84.489 |                0.658 |
| V2010_test | Arm B    | quiet   | c1         |      367 |          415.000 |             0.008 |            0.005 |        0.036 |       0.029 |      13.175 |              12.075 |                0.916 |
| V2010_test | Arm B    | episode | c1         |     1440 |         1453.000 |             0.043 |            0.012 |        0.080 |       0.050 |     115.227 |              72.414 |                0.628 |
| V2008_test | ALL      | all     | c1         |    93754 |        90700.000 |             0.027 |            0.004 |        0.052 |       0.031 |    4849.326 |            2786.706 |                0.575 |
| V2008_test | ALL      | quiet   | c1         |    13913 |        13285.000 |             0.001 |            0.001 |        0.024 |       0.025 |     335.474 |             337.828 |                1.007 |
| V2008_test | ALL      | episode | c1         |    79841 |        77415.000 |             0.032 |            0.005 |        0.057 |       0.032 |    4513.852 |            2448.878 |                0.543 |
| V2008_test | NW       | all     | c1         |    16322 |        13531.000 |             0.147 |            0.025 |        0.208 |       0.108 |    3399.860 |            1457.658 |                0.429 |
| V2008_test | NW       | quiet   | c1         |     1993 |         1980.000 |             0.005 |            0.007 |        0.056 |       0.060 |     111.162 |             118.335 |                1.065 |
| V2008_test | NW       | episode | c1         |    14329 |        11551.000 |             0.166 |            0.028 |        0.230 |       0.116 |    3288.697 |            1339.323 |                0.407 |
| V2008_test | SW       | all     | c1         |      305 |          153.000 |             0.407 |            0.170 |        0.451 |       0.250 |     137.669 |              38.279 |                0.278 |
| V2008_test | SW       | quiet   | c1         |        9 |           11.000 |             0.111 |            0.000 |        0.213 |       0.149 |       1.916 |               1.644 |                0.858 |
| V2008_test | SW       | episode | c1         |      296 |          142.000 |             0.416 |            0.183 |        0.459 |       0.258 |     135.753 |              36.635 |                0.270 |
| V2008_test | Olympic  | all     | c1         |    73842 |        73854.000 |             0.000 |            0.000 |        0.017 |       0.017 |    1247.366 |            1236.794 |                0.992 |
| V2008_test | Olympic  | quiet   | c1         |    11474 |        10883.000 |             0.000 |            0.000 |        0.019 |       0.019 |     215.537 |             211.748 |                0.982 |
| V2008_test | Olympic  | episode | c1         |    62368 |        62971.000 |             0.001 |            0.000 |        0.017 |       0.016 |    1031.830 |            1025.046 |                0.993 |
| V2008_test | Arm B    | all     | c1         |     3285 |         3162.000 |             0.002 |            0.000 |        0.020 |       0.017 |      64.431 |              53.975 |                0.838 |
| V2008_test | Arm B    | quiet   | c1         |      437 |          411.000 |             0.002 |            0.000 |        0.016 |       0.015 |       6.859 |               6.100 |                0.889 |
| V2008_test | Arm B    | episode | c1         |     2848 |         2751.000 |             0.002 |            0.000 |        0.020 |       0.017 |      57.572 |              47.874 |                0.832 |
| A2011      | ALL      | all     | c1         |   294824 |       287560.000 |             0.027 |            0.009 |        0.056 |       0.038 |   16480.872 |           11056.271 |                0.671 |
| A2011      | ALL      | episode | c1         |   294824 |       287560.000 |             0.027 |            0.009 |        0.056 |       0.038 |   16480.872 |           11056.271 |                0.671 |
| A2011      | NW       | all     | c1         |    61999 |        55860.000 |             0.113 |            0.030 |        0.175 |       0.105 |   10857.238 |            5845.893 |                0.538 |
| A2011      | NW       | episode | c1         |    61999 |        55860.000 |             0.113 |            0.030 |        0.175 |       0.105 |   10857.238 |            5845.893 |                0.538 |
| A2011      | SW       | all     | c1         |      709 |          714.000 |             0.436 |            0.438 |        0.478 |       0.483 |     339.061 |             344.895 |                1.017 |
| A2011      | SW       | episode | c1         |      709 |          714.000 |             0.436 |            0.438 |        0.478 |       0.483 |     339.061 |             344.895 |                1.017 |
| A2011      | Olympic  | all     | c1         |   215645 |       214694.000 |             0.002 |            0.001 |        0.019 |       0.017 |    4044.436 |            3741.641 |                0.925 |
| A2011      | Olympic  | episode | c1         |   215645 |       214694.000 |             0.002 |            0.001 |        0.019 |       0.017 |    4044.436 |            3741.641 |                0.925 |
| A2011      | Arm B    | all     | c1         |    16471 |        16292.000 |             0.022 |            0.015 |        0.075 |       0.069 |    1240.137 |            1123.842 |                0.906 |
| A2011      | Arm B    | episode | c1         |    16471 |        16292.000 |             0.022 |            0.015 |        0.075 |       0.069 |    1240.137 |            1123.842 |                0.906 |

## (c)/(d) Scrambled control and quiet vs episode days — model D1_hgb_FULL_any

FULL models see step-2 attribution, which exists for the first control only → c1; NO_STEP2 uses all 4 controls. `ctl_over_data_sumP` = expected 'real' count the model would assign to one scrambled network ÷ the data's.

| eval_set   | region   | days    | controls   |   n_data |   n_ctl_per_real |   data_fracP_gt05 |   ctl_fracP_gt05 |   data_meanP |   ctl_meanP |   data_sumP |   ctl_sumP_per_real |   ctl_over_data_sumP |
|:-----------|:---------|:--------|:-----------|---------:|-----------------:|------------------:|-----------------:|-------------:|------------:|------------:|--------------------:|---------------------:|
| V2010_test | ALL      | all     | c1         |   172854 |       169413.000 |             0.114 |            0.095 |        0.252 |       0.239 |   43569.019 |           40523.638 |                0.930 |
| V2010_test | ALL      | quiet   | c1         |    48576 |        49192.000 |             0.080 |            0.070 |        0.226 |       0.219 |   10991.383 |           10774.652 |                0.980 |
| V2010_test | ALL      | episode | c1         |   124278 |       120221.000 |             0.128 |            0.105 |        0.262 |       0.247 |   32577.636 |           29748.986 |                0.913 |
| V2010_test | NW       | all     | c1         |    21800 |        18954.000 |             0.213 |            0.111 |        0.342 |       0.271 |    7448.927 |            5134.304 |                0.689 |
| V2010_test | NW       | quiet   | c1         |     5429 |         5599.000 |             0.019 |            0.011 |        0.186 |       0.180 |    1011.145 |            1006.230 |                0.995 |
| V2010_test | NW       | episode | c1         |    16371 |        13355.000 |             0.277 |            0.152 |        0.393 |       0.309 |    6437.782 |            4128.074 |                0.641 |
| V2010_test | SW       | all     | c1         |      381 |          169.000 |             0.391 |            0.047 |        0.460 |       0.235 |     175.351 |              39.722 |                0.227 |
| V2010_test | SW       | quiet   | c1         |       33 |           44.000 |             0.061 |            0.068 |        0.240 |       0.217 |       7.909 |               9.531 |                1.205 |
| V2010_test | SW       | episode | c1         |      348 |          125.000 |             0.422 |            0.040 |        0.481 |       0.242 |     167.442 |              30.191 |                0.180 |
| V2010_test | Olympic  | all     | c1         |   148866 |       148422.000 |             0.098 |            0.091 |        0.237 |       0.233 |   35237.937 |           34637.221 |                0.983 |
| V2010_test | Olympic  | quiet   | c1         |    42747 |        43134.000 |             0.086 |            0.076 |        0.230 |       0.223 |    9839.062 |            9608.041 |                0.977 |
| V2010_test | Olympic  | episode | c1         |   106119 |       105288.000 |             0.102 |            0.098 |        0.239 |       0.238 |   25398.875 |           25029.180 |                0.985 |
| V2010_test | Arm B    | all     | c1         |     1807 |         1868.000 |             0.229 |            0.209 |        0.391 |       0.381 |     706.804 |             712.391 |                1.008 |
| V2010_test | Arm B    | quiet   | c1         |      367 |          415.000 |             0.213 |            0.227 |        0.363 |       0.363 |     133.267 |             150.850 |                1.132 |
| V2010_test | Arm B    | episode | c1         |     1440 |         1453.000 |             0.233 |            0.204 |        0.398 |       0.386 |     573.537 |             561.540 |                0.979 |
| V2008_test | ALL      | all     | c1         |    93754 |        90700.000 |             0.073 |            0.045 |        0.237 |       0.220 |   22188.818 |           19985.394 |                0.901 |
| V2008_test | ALL      | quiet   | c1         |    13913 |        13285.000 |             0.050 |            0.046 |        0.222 |       0.224 |    3087.550 |            2981.529 |                0.966 |
| V2008_test | ALL      | episode | c1         |    79841 |        77415.000 |             0.077 |            0.044 |        0.239 |       0.220 |   19101.268 |           17003.866 |                0.890 |
| V2008_test | NW       | all     | c1         |    16322 |        13531.000 |             0.250 |            0.117 |        0.365 |       0.274 |    5954.644 |            3705.218 |                0.622 |
| V2008_test | NW       | quiet   | c1         |     1993 |         1980.000 |             0.009 |            0.009 |        0.150 |       0.153 |     298.363 |             302.248 |                1.013 |
| V2008_test | NW       | episode | c1         |    14329 |        11551.000 |             0.284 |            0.136 |        0.395 |       0.295 |    5656.282 |            3402.971 |                0.602 |
| V2008_test | SW       | all     | c1         |      305 |          153.000 |             0.544 |            0.209 |        0.587 |       0.363 |     178.935 |              55.489 |                0.310 |
| V2008_test | SW       | quiet   | c1         |        9 |           11.000 |             0.111 |            0.000 |        0.271 |       0.197 |       2.436 |               2.167 |                0.890 |
| V2008_test | SW       | episode | c1         |      296 |          142.000 |             0.557 |            0.225 |        0.596 |       0.376 |     176.499 |              53.322 |                0.302 |
| V2008_test | Olympic  | all     | c1         |    73842 |        73854.000 |             0.034 |            0.033 |        0.206 |       0.209 |   15242.783 |           15440.351 |                1.013 |
| V2008_test | Olympic  | quiet   | c1         |    11474 |        10883.000 |             0.058 |            0.054 |        0.234 |       0.237 |    2684.276 |            2579.656 |                0.961 |
| V2008_test | Olympic  | episode | c1         |    62368 |        62971.000 |             0.030 |            0.029 |        0.201 |       0.204 |   12558.507 |           12860.695 |                1.024 |
| V2008_test | Arm B    | all     | c1         |     3285 |         3162.000 |             0.012 |            0.008 |        0.247 |       0.248 |     812.455 |             784.337 |                0.965 |
| V2008_test | Arm B    | quiet   | c1         |      437 |          411.000 |             0.009 |            0.010 |        0.234 |       0.237 |     102.475 |              97.458 |                0.951 |
| V2008_test | Arm B    | episode | c1         |     2848 |         2751.000 |             0.013 |            0.008 |        0.249 |       0.250 |     709.980 |             686.878 |                0.967 |
| A2011      | ALL      | all     | c1         |   294824 |       287560.000 |             0.090 |            0.068 |        0.261 |       0.249 |   76957.606 |           71558.212 |                0.930 |
| A2011      | ALL      | episode | c1         |   294824 |       287560.000 |             0.090 |            0.068 |        0.261 |       0.249 |   76957.606 |           71558.212 |                0.930 |
| A2011      | NW       | all     | c1         |    61999 |        55860.000 |             0.203 |            0.114 |        0.345 |       0.289 |   21379.103 |           16162.934 |                0.756 |
| A2011      | NW       | episode | c1         |    61999 |        55860.000 |             0.203 |            0.114 |        0.345 |       0.289 |   21379.103 |           16162.934 |                0.756 |
| A2011      | SW       | all     | c1         |      709 |          714.000 |             0.570 |            0.611 |        0.551 |       0.564 |     390.541 |             402.503 |                1.031 |
| A2011      | SW       | episode | c1         |      709 |          714.000 |             0.570 |            0.611 |        0.551 |       0.564 |     390.541 |             402.503 |                1.031 |
| A2011      | Olympic  | all     | c1         |   215645 |       214694.000 |             0.056 |            0.053 |        0.234 |       0.234 |   50370.180 |           50310.567 |                0.999 |
| A2011      | Olympic  | episode | c1         |   215645 |       214694.000 |             0.056 |            0.053 |        0.234 |       0.234 |   50370.180 |           50310.567 |                0.999 |
| A2011      | Arm B    | all     | c1         |    16471 |        16292.000 |             0.090 |            0.081 |        0.293 |       0.287 |    4817.781 |            4682.208 |                0.972 |
| A2011      | Arm B    | episode | c1         |    16471 |        16292.000 |             0.090 |            0.081 |        0.293 |       0.287 |    4817.781 |            4682.208 |                0.972 |

## (c)/(d) Scrambled control and quiet vs episode days — model P1

FULL models see step-2 attribution, which exists for the first control only → c1; NO_STEP2 uses all 4 controls. `ctl_over_data_sumP` = expected 'real' count the model would assign to one scrambled network ÷ the data's.

| eval_set   | region   | days    | controls    |   n_data |   n_ctl_per_real |   data_fracP_gt05 |   ctl_fracP_gt05 |   data_meanP |   ctl_meanP |   data_sumP |   ctl_sumP_per_real |   ctl_over_data_sumP |
|:-----------|:---------|:--------|:------------|---------:|-----------------:|------------------:|-----------------:|-------------:|------------:|------------:|--------------------:|---------------------:|
| V2010_test | ALL      | all     | c1,c2,c3,c4 |   172854 |       168969.250 |             0.021 |            0.004 |        0.051 |       0.035 |    8795.923 |            5952.251 |                0.677 |
| V2010_test | ALL      | quiet   | c1,c2,c3,c4 |    48576 |        48451.250 |             0.020 |            0.010 |        0.047 |       0.040 |    2305.872 |            1928.991 |                0.837 |
| V2010_test | ALL      | episode | c1,c2,c3,c4 |   124278 |       120518.000 |             0.021 |            0.002 |        0.052 |       0.033 |    6490.051 |            4023.260 |                0.620 |
| V2010_test | NW       | all     | c1,c2,c3,c4 |    21800 |        19002.500 |             0.106 |            0.011 |        0.125 |       0.030 |    2723.345 |             561.451 |                0.206 |
| V2010_test | NW       | quiet   | c1,c2,c3,c4 |     5429 |         5521.000 |             0.006 |            0.004 |        0.025 |       0.024 |     133.914 |             130.144 |                0.972 |
| V2010_test | NW       | episode | c1,c2,c3,c4 |    16371 |        13481.500 |             0.140 |            0.013 |        0.158 |       0.032 |    2589.430 |             431.308 |                0.167 |
| V2010_test | SW       | all     | c1,c2,c3,c4 |      381 |          167.750 |             0.386 |            0.024 |        0.417 |       0.077 |     158.714 |              12.950 |                0.082 |
| V2010_test | SW       | quiet   | c1,c2,c3,c4 |       33 |           44.250 |             0.061 |            0.023 |        0.112 |       0.070 |       3.688 |               3.119 |                0.846 |
| V2010_test | SW       | episode | c1,c2,c3,c4 |      348 |          123.500 |             0.417 |            0.024 |        0.445 |       0.080 |     155.026 |               9.831 |                0.063 |
| V2010_test | Olympic  | all     | c1,c2,c3,c4 |   148866 |       147917.250 |             0.007 |            0.003 |        0.039 |       0.036 |    5827.517 |            5328.555 |                0.914 |
| V2010_test | Olympic  | quiet   | c1,c2,c3,c4 |    42747 |        42469.000 |             0.022 |            0.011 |        0.051 |       0.042 |    2160.446 |            1786.084 |                0.827 |
| V2010_test | Olympic  | episode | c1,c2,c3,c4 |   106119 |       105448.250 |             0.001 |            0.000 |        0.035 |       0.034 |    3667.071 |            3542.470 |                0.966 |
| V2010_test | Arm B    | all     | c1,c2,c3,c4 |     1807 |         1881.750 |             0.032 |            0.010 |        0.048 |       0.026 |      86.347 |              49.295 |                0.571 |
| V2010_test | Arm B    | quiet   | c1,c2,c3,c4 |      367 |          417.000 |             0.008 |            0.007 |        0.021 |       0.023 |       7.824 |               9.644 |                1.233 |
| V2010_test | Arm B    | episode | c1,c2,c3,c4 |     1440 |         1464.750 |             0.037 |            0.011 |        0.055 |       0.027 |      78.523 |              39.651 |                0.505 |
| V2008_test | ALL      | all     | c1,c2,c3,c4 |    93754 |        90829.500 |             0.029 |            0.007 |        0.055 |       0.038 |    5110.635 |            3437.870 |                0.673 |
| V2008_test | ALL      | quiet   | c1,c2,c3,c4 |    13913 |        13121.500 |             0.001 |            0.008 |        0.026 |       0.040 |     359.342 |             527.095 |                1.467 |
| V2008_test | ALL      | episode | c1,c2,c3,c4 |    79841 |        77708.000 |             0.033 |            0.007 |        0.060 |       0.037 |    4751.293 |            2910.775 |                0.613 |
| V2008_test | NW       | all     | c1,c2,c3,c4 |    16322 |        13546.250 |             0.128 |            0.008 |        0.148 |       0.029 |    2411.336 |             386.823 |                0.160 |
| V2008_test | NW       | quiet   | c1,c2,c3,c4 |     1993 |         1966.000 |             0.003 |            0.004 |        0.027 |       0.029 |      53.063 |              56.271 |                1.060 |
| V2008_test | NW       | episode | c1,c2,c3,c4 |    14329 |        11580.250 |             0.145 |            0.009 |        0.165 |       0.029 |    2358.273 |             330.552 |                0.140 |
| V2008_test | SW       | all     | c1,c2,c3,c4 |      305 |          145.250 |             0.518 |            0.071 |        0.533 |       0.138 |     162.531 |              20.076 |                0.124 |
| V2008_test | SW       | quiet   | c1,c2,c3,c4 |        9 |           10.750 |             0.111 |            0.070 |        0.144 |       0.120 |       1.293 |               1.289 |                0.997 |
| V2008_test | SW       | episode | c1,c2,c3,c4 |      296 |          134.500 |             0.530 |            0.071 |        0.545 |       0.140 |     161.238 |              18.788 |                0.117 |
| V2008_test | Olympic  | all     | c1,c2,c3,c4 |    73842 |        73895.750 |             0.006 |            0.007 |        0.034 |       0.040 |    2477.650 |            2978.939 |                1.202 |
| V2008_test | Olympic  | quiet   | c1,c2,c3,c4 |    11474 |        10724.500 |             0.000 |            0.009 |        0.026 |       0.043 |     295.724 |             462.799 |                1.565 |
| V2008_test | Olympic  | episode | c1,c2,c3,c4 |    62368 |        63171.250 |             0.007 |            0.007 |        0.035 |       0.040 |    2181.925 |            2516.141 |                1.153 |
| V2008_test | Arm B    | all     | c1,c2,c3,c4 |     3285 |         3242.250 |             0.001 |            0.001 |        0.018 |       0.016 |      59.119 |              52.031 |                0.880 |
| V2008_test | Arm B    | quiet   | c1,c2,c3,c4 |      437 |          420.250 |             0.000 |            0.002 |        0.021 |       0.016 |       9.261 |               6.736 |                0.727 |
| V2008_test | Arm B    | episode | c1,c2,c3,c4 |     2848 |         2822.000 |             0.001 |            0.001 |        0.018 |       0.016 |      49.857 |              45.295 |                0.908 |
| A2011      | ALL      | all     | c1,c2,c3,c4 |   294824 |       287372.500 |             0.028 |            0.014 |        0.084 |       0.070 |   24821.464 |           20152.425 |                0.812 |
| A2011      | ALL      | episode | c1,c2,c3,c4 |   294824 |       287372.500 |             0.028 |            0.014 |        0.084 |       0.070 |   24821.464 |           20152.425 |                0.812 |
| A2011      | NW       | all     | c1,c2,c3,c4 |    61999 |        55640.250 |             0.085 |            0.017 |        0.137 |       0.071 |    8464.658 |            3961.714 |                0.468 |
| A2011      | NW       | episode | c1,c2,c3,c4 |    61999 |        55640.250 |             0.085 |            0.017 |        0.137 |       0.071 |    8464.658 |            3961.714 |                0.468 |
| A2011      | SW       | all     | c1,c2,c3,c4 |      709 |          699.750 |             0.477 |            0.488 |        0.512 |       0.527 |     362.880 |             369.081 |                1.017 |
| A2011      | SW       | episode | c1,c2,c3,c4 |      709 |          699.750 |             0.477 |            0.488 |        0.512 |       0.527 |     362.880 |             369.081 |                1.017 |
| A2011      | Olympic  | all     | c1,c2,c3,c4 |   215645 |       214685.000 |             0.011 |            0.011 |        0.068 |       0.068 |   14677.801 |           14610.222 |                0.995 |
| A2011      | Olympic  | episode | c1,c2,c3,c4 |   215645 |       214685.000 |             0.011 |            0.011 |        0.068 |       0.068 |   14677.801 |           14610.222 |                0.995 |
| A2011      | Arm B    | all     | c1,c2,c3,c4 |    16471 |        16347.500 |             0.023 |            0.016 |        0.080 |       0.074 |    1316.125 |            1211.407 |                0.920 |
| A2011      | Arm B    | episode | c1,c2,c3,c4 |    16471 |        16347.500 |             0.023 |            0.016 |        0.080 |       0.074 |    1316.125 |            1211.407 |                0.920 |

### By PNSN regional tremor-day label (label only, never a feature)

| eval_set   | region   | tremor_day   | model           |   n_data |   n_c1 |   data_fracP_gt05 |   c1_fracP_gt05 |   data_meanP |   c1_meanP |
|:-----------|:---------|:-------------|:----------------|---------:|-------:|------------------:|----------------:|-------------:|-----------:|
| V2010_test | ALL      | False        | D2_hgb_FULL_own |    49015 |  49667 |             0.003 |           0.001 |        0.039 |      0.036 |
| V2010_test | ALL      | False        | D2_hgb_FULL_any |    49015 |  49667 |             0.086 |           0.074 |        0.228 |      0.221 |
| V2010_test | ALL      | False        | D1_hgb_FULL_own |    49015 |  49667 |             0.009 |           0.004 |        0.032 |      0.027 |
| V2010_test | ALL      | False        | D1_hgb_FULL_any |    49015 |  49667 |             0.079 |           0.070 |        0.226 |      0.219 |
| V2010_test | ALL      | False        | P1              |    49015 |  49667 |             0.020 |           0.011 |        0.047 |      0.041 |
| V2010_test | ALL      | True         | D2_hgb_FULL_own |   123839 | 119746 |             0.030 |           0.009 |        0.068 |      0.052 |
| V2010_test | ALL      | True         | D2_hgb_FULL_any |   123839 | 119746 |             0.136 |           0.114 |        0.262 |      0.249 |
| V2010_test | ALL      | True         | D1_hgb_FULL_own |   123839 | 119746 |             0.026 |           0.006 |        0.052 |      0.034 |
| V2010_test | ALL      | True         | D1_hgb_FULL_any |   123839 | 119746 |             0.128 |           0.105 |        0.262 |      0.248 |
| V2010_test | ALL      | True         | P1              |   123839 | 119746 |             0.021 |           0.002 |        0.052 |      0.034 |
| V2010_test | NW       | False        | D2_hgb_FULL_own |     5858 |   6062 |             0.002 |           0.001 |        0.062 |      0.060 |
| V2010_test | NW       | False        | D2_hgb_FULL_any |     5858 |   6062 |             0.022 |           0.016 |        0.192 |      0.186 |
| V2010_test | NW       | False        | D1_hgb_FULL_own |     5858 |   6062 |             0.007 |           0.006 |        0.063 |      0.061 |
| V2010_test | NW       | False        | D1_hgb_FULL_any |     5858 |   6062 |             0.020 |           0.014 |        0.189 |      0.183 |
| V2010_test | NW       | False        | P1              |     5858 |   6062 |             0.006 |           0.004 |        0.024 |      0.024 |
| V2010_test | NW       | True         | D2_hgb_FULL_own |    15942 |  12892 |             0.185 |           0.055 |        0.248 |      0.156 |
| V2010_test | NW       | True         | D2_hgb_FULL_any |    15942 |  12892 |             0.313 |           0.188 |        0.402 |      0.323 |
| V2010_test | NW       | True         | D1_hgb_FULL_own |    15942 |  12892 |             0.148 |           0.024 |        0.210 |      0.108 |
| V2010_test | NW       | True         | D1_hgb_FULL_any |    15942 |  12892 |             0.283 |           0.156 |        0.398 |      0.312 |
| V2010_test | NW       | True         | P1              |    15942 |  12892 |             0.143 |           0.015 |        0.162 |      0.034 |
| V2010_test | SW       | False        | D2_hgb_FULL_own |       43 |     56 |             0.023 |           0.000 |        0.063 |      0.046 |
| V2010_test | SW       | False        | D2_hgb_FULL_any |       43 |     56 |             0.047 |           0.000 |        0.174 |      0.151 |
| V2010_test | SW       | False        | D1_hgb_FULL_own |       43 |     56 |             0.070 |           0.071 |        0.159 |      0.151 |
| V2010_test | SW       | False        | D1_hgb_FULL_any |       43 |     56 |             0.047 |           0.071 |        0.231 |      0.222 |
| V2010_test | SW       | False        | P1              |       43 |     56 |             0.047 |           0.054 |        0.092 |      0.096 |
| V2010_test | SW       | True         | D2_hgb_FULL_own |      338 |    113 |             0.518 |           0.009 |        0.478 |      0.105 |
| V2010_test | SW       | True         | D2_hgb_FULL_any |      338 |    113 |             0.571 |           0.080 |        0.551 |      0.240 |
| V2010_test | SW       | True         | D1_hgb_FULL_own |      338 |    113 |             0.349 |           0.035 |        0.369 |      0.138 |
| V2010_test | SW       | True         | D1_hgb_FULL_any |      338 |    113 |             0.435 |           0.035 |        0.489 |      0.242 |
| V2010_test | SW       | True         | P1              |      338 |    113 |             0.429 |           0.000 |        0.458 |      0.056 |
| V2010_test | Olympic  | False        | D2_hgb_FULL_own |    42747 |  43134 |             0.003 |           0.002 |        0.036 |      0.033 |
| V2010_test | Olympic  | False        | D2_hgb_FULL_any |    42747 |  43134 |             0.094 |           0.082 |        0.233 |      0.225 |
| V2010_test | Olympic  | False        | D1_hgb_FULL_own |    42747 |  43134 |             0.009 |           0.004 |        0.028 |      0.022 |
| V2010_test | Olympic  | False        | D1_hgb_FULL_any |    42747 |  43134 |             0.086 |           0.076 |        0.230 |      0.223 |
| V2010_test | Olympic  | False        | P1              |    42747 |  43134 |             0.022 |           0.012 |        0.051 |      0.043 |
| V2010_test | Olympic  | True         | D2_hgb_FULL_own |   106119 | 105288 |             0.005 |           0.004 |        0.040 |      0.039 |
| V2010_test | Olympic  | True         | D2_hgb_FULL_any |   106119 | 105288 |             0.107 |           0.104 |        0.239 |      0.238 |
| V2010_test | Olympic  | True         | D1_hgb_FULL_own |   106119 | 105288 |             0.006 |           0.004 |        0.027 |      0.025 |
| V2010_test | Olympic  | True         | D1_hgb_FULL_any |   106119 | 105288 |             0.102 |           0.098 |        0.239 |      0.238 |
| V2010_test | Olympic  | True         | P1              |   106119 | 105288 |             0.001 |           0.000 |        0.035 |      0.034 |
| V2010_test | Arm B    | False        | D2_hgb_FULL_own |      367 |    415 |             0.000 |           0.000 |        0.029 |      0.029 |
| V2010_test | Arm B    | False        | D2_hgb_FULL_any |      367 |    415 |             0.139 |           0.128 |        0.318 |      0.322 |
| V2010_test | Arm B    | False        | D1_hgb_FULL_own |      367 |    415 |             0.008 |           0.005 |        0.036 |      0.029 |
| V2010_test | Arm B    | False        | D1_hgb_FULL_any |      367 |    415 |             0.213 |           0.227 |        0.363 |      0.363 |
| V2010_test | Arm B    | False        | P1              |      367 |    415 |             0.008 |           0.007 |        0.021 |      0.020 |
| V2010_test | Arm B    | True         | D2_hgb_FULL_own |     1440 |   1453 |             0.014 |           0.000 |        0.062 |      0.047 |
| V2010_test | Arm B    | True         | D2_hgb_FULL_any |     1440 |   1453 |             0.205 |           0.185 |        0.365 |      0.353 |
| V2010_test | Arm B    | True         | D1_hgb_FULL_own |     1440 |   1453 |             0.043 |           0.012 |        0.080 |      0.050 |
| V2010_test | Arm B    | True         | D1_hgb_FULL_any |     1440 |   1453 |             0.233 |           0.204 |        0.398 |      0.386 |
| V2010_test | Arm B    | True         | P1              |     1440 |   1453 |             0.037 |           0.014 |        0.055 |      0.029 |
| V2008_test | ALL      | False        | D2_hgb_FULL_own |    93754 |  90700 |             0.035 |           0.009 |        0.078 |      0.058 |
| V2008_test | ALL      | False        | D2_hgb_FULL_any |    93754 |  90700 |             0.177 |           0.152 |        0.302 |      0.289 |
| V2008_test | ALL      | False        | D1_hgb_FULL_own |    93754 |  90700 |             0.027 |           0.004 |        0.052 |      0.031 |
| V2008_test | ALL      | False        | D1_hgb_FULL_any |    93754 |  90700 |             0.073 |           0.045 |        0.237 |      0.220 |
| V2008_test | ALL      | False        | P1              |    93754 |  90700 |             0.029 |           0.009 |        0.055 |      0.040 |
| V2008_test | NW       | False        | D2_hgb_FULL_own |    16322 |  13531 |             0.181 |           0.050 |        0.244 |      0.153 |
| V2008_test | NW       | False        | D2_hgb_FULL_any |    16322 |  13531 |             0.294 |           0.159 |        0.378 |      0.294 |
| V2008_test | NW       | False        | D1_hgb_FULL_own |    16322 |  13531 |             0.147 |           0.025 |        0.208 |      0.108 |
| V2008_test | NW       | False        | D1_hgb_FULL_any |    16322 |  13531 |             0.250 |           0.117 |        0.365 |      0.274 |
| V2008_test | NW       | False        | P1              |    16322 |  13531 |             0.128 |           0.008 |        0.148 |      0.029 |
| V2008_test | SW       | False        | D2_hgb_FULL_own |      305 |    153 |             0.472 |           0.026 |        0.461 |      0.133 |
| V2008_test | SW       | False        | D2_hgb_FULL_any |      305 |    153 |             0.502 |           0.065 |        0.510 |      0.234 |
| V2008_test | SW       | False        | D1_hgb_FULL_own |      305 |    153 |             0.407 |           0.170 |        0.451 |      0.250 |
| V2008_test | SW       | False        | D1_hgb_FULL_any |      305 |    153 |             0.544 |           0.209 |        0.587 |      0.363 |
| V2008_test | SW       | False        | P1              |      305 |    153 |             0.518 |           0.046 |        0.533 |      0.133 |
| V2008_test | Olympic  | False        | D2_hgb_FULL_own |    73842 |  73854 |             0.002 |           0.001 |        0.041 |      0.042 |
| V2008_test | Olympic  | False        | D2_hgb_FULL_any |    73842 |  73854 |             0.152 |           0.153 |        0.283 |      0.286 |
| V2008_test | Olympic  | False        | D1_hgb_FULL_own |    73842 |  73854 |             0.000 |           0.000 |        0.017 |      0.017 |
| V2008_test | Olympic  | False        | D1_hgb_FULL_any |    73842 |  73854 |             0.034 |           0.033 |        0.206 |      0.209 |
| V2008_test | Olympic  | False        | P1              |    73842 |  73854 |             0.006 |           0.009 |        0.034 |      0.043 |
| V2008_test | Arm B    | False        | D2_hgb_FULL_own |     3285 |   3162 |             0.001 |           0.000 |        0.034 |      0.032 |
| V2008_test | Arm B    | False        | D2_hgb_FULL_any |     3285 |   3162 |             0.126 |           0.118 |        0.334 |      0.336 |
| V2008_test | Arm B    | False        | D1_hgb_FULL_own |     3285 |   3162 |             0.002 |           0.000 |        0.020 |      0.017 |
| V2008_test | Arm B    | False        | D1_hgb_FULL_any |     3285 |   3162 |             0.012 |           0.008 |        0.247 |      0.248 |
| V2008_test | Arm B    | False        | P1              |     3285 |   3162 |             0.001 |           0.001 |        0.018 |      0.015 |
| A2011      | ALL      | False        | D2_hgb_FULL_own |     1379 |   1381 |             0.007 |           0.007 |        0.116 |      0.112 |
| A2011      | ALL      | False        | D2_hgb_FULL_any |     1379 |   1381 |             0.066 |           0.061 |        0.300 |      0.296 |
| A2011      | ALL      | False        | D1_hgb_FULL_own |     1379 |   1381 |             0.053 |           0.063 |        0.163 |      0.164 |
| A2011      | ALL      | False        | D1_hgb_FULL_any |     1379 |   1381 |             0.099 |           0.113 |        0.310 |      0.311 |
| A2011      | ALL      | False        | P1              |     1379 |   1381 |             0.047 |           0.056 |        0.187 |      0.192 |
| A2011      | ALL      | True         | D2_hgb_FULL_own |   293445 | 286179 |             0.024 |           0.007 |        0.082 |      0.069 |
| A2011      | ALL      | True         | D2_hgb_FULL_any |   293445 | 286179 |             0.187 |           0.170 |        0.320 |      0.309 |
| A2011      | ALL      | True         | D1_hgb_FULL_own |   293445 | 286179 |             0.027 |           0.008 |        0.055 |      0.038 |
| A2011      | ALL      | True         | D1_hgb_FULL_any |   293445 | 286179 |             0.090 |           0.068 |        0.261 |      0.249 |
| A2011      | ALL      | True         | P1              |   293445 | 286179 |             0.028 |           0.014 |        0.084 |      0.070 |
| A2011      | NW       | False        | D2_hgb_FULL_own |     1202 |   1218 |             0.007 |           0.007 |        0.111 |      0.106 |
| A2011      | NW       | False        | D2_hgb_FULL_any |     1202 |   1218 |             0.062 |           0.057 |        0.304 |      0.298 |
| A2011      | NW       | False        | D1_hgb_FULL_own |     1202 |   1218 |             0.017 |           0.022 |        0.126 |      0.128 |
| A2011      | NW       | False        | D1_hgb_FULL_any |     1202 |   1218 |             0.055 |           0.062 |        0.285 |      0.285 |
| A2011      | NW       | False        | P1              |     1202 |   1218 |             0.027 |           0.034 |        0.167 |      0.170 |
| A2011      | NW       | True         | D2_hgb_FULL_own |    60797 |  54642 |             0.113 |           0.036 |        0.189 |      0.133 |
| A2011      | NW       | True         | D2_hgb_FULL_any |    60797 |  54642 |             0.254 |           0.172 |        0.371 |      0.320 |
| A2011      | NW       | True         | D1_hgb_FULL_own |    60797 |  54642 |             0.115 |           0.030 |        0.176 |      0.104 |
| A2011      | NW       | True         | D1_hgb_FULL_any |    60797 |  54642 |             0.206 |           0.115 |        0.346 |      0.289 |
| A2011      | NW       | True         | P1              |    60797 |  54642 |             0.086 |           0.017 |        0.136 |      0.069 |
| A2011      | SW       | False        | D2_hgb_FULL_own |      177 |    163 |             0.006 |           0.006 |        0.152 |      0.161 |
| A2011      | SW       | False        | D2_hgb_FULL_any |      177 |    163 |             0.090 |           0.086 |        0.274 |      0.288 |
| A2011      | SW       | False        | D1_hgb_FULL_own |      177 |    163 |             0.299 |           0.368 |        0.412 |      0.438 |
| A2011      | SW       | False        | D1_hgb_FULL_any |      177 |    163 |             0.401 |           0.491 |        0.483 |      0.510 |
| A2011      | SW       | False        | P1              |      177 |    163 |             0.186 |           0.227 |        0.322 |      0.355 |
| A2011      | SW       | True         | D2_hgb_FULL_own |      532 |    551 |             0.002 |           0.000 |        0.140 |      0.142 |
| A2011      | SW       | True         | D2_hgb_FULL_any |      532 |    551 |             0.034 |           0.045 |        0.303 |      0.307 |
| A2011      | SW       | True         | D1_hgb_FULL_own |      532 |    551 |             0.481 |           0.459 |        0.500 |      0.496 |
| A2011      | SW       | True         | D1_hgb_FULL_any |      532 |    551 |             0.626 |           0.646 |        0.573 |      0.580 |
| A2011      | SW       | True         | P1              |      532 |    551 |             0.573 |           0.588 |        0.575 |      0.589 |
| A2011      | Olympic  | True         | D2_hgb_FULL_own |   215645 | 214694 |             0.000 |           0.000 |        0.055 |      0.054 |
| A2011      | Olympic  | True         | D2_hgb_FULL_any |   215645 | 214694 |             0.177 |           0.177 |        0.308 |      0.309 |
| A2011      | Olympic  | True         | D1_hgb_FULL_own |   215645 | 214694 |             0.002 |           0.001 |        0.019 |      0.017 |
| A2011      | Olympic  | True         | D1_hgb_FULL_any |   215645 | 214694 |             0.056 |           0.053 |        0.234 |      0.234 |
| A2011      | Olympic  | True         | P1              |   215645 | 214694 |             0.011 |           0.012 |        0.068 |      0.068 |
| A2011      | Arm B    | True         | D2_hgb_FULL_own |    16471 |  16292 |             0.001 |           0.000 |        0.037 |      0.035 |
| A2011      | Arm B    | True         | D2_hgb_FULL_any |    16471 |  16292 |             0.074 |           0.069 |        0.277 |      0.276 |
| A2011      | Arm B    | True         | D1_hgb_FULL_own |    16471 |  16292 |             0.022 |           0.015 |        0.075 |      0.069 |
| A2011      | Arm B    | True         | D1_hgb_FULL_any |    16471 |  16292 |             0.090 |           0.081 |        0.293 |      0.287 |
| A2011      | Arm B    | True         | P1              |    16471 |  16292 |             0.023 |           0.016 |        0.080 |      0.074 |

## (e) Expected count — ΣP per family-day vs the true array count (episode family-days the family's array run scanned; scored half)

The truth counts every array detection, including the many the permanent-station MF never produced a candidate for, so ratio < 1 is recall, not mis-calibration; the rank agreement is the test.

| eval_set   | region   | model              | truth             |   n_famdays |      sumP |   truth_count |   ratio |   spearman |   spearman_ncand |
|:-----------|:---------|:-------------------|:------------------|------------:|----------:|--------------:|--------:|-----------:|-----------------:|
| V2010_test | ALL      | D2_hgb_FULL_own    | truth_own_scored  |        4572 |  7790.995 |         30921 |   0.252 |      0.524 |            0.458 |
| V2010_test | ALL      | D1_hgb_FULL_own    | truth_own_scored  |        4572 |  5875.657 |         30921 |   0.190 |      0.527 |            0.458 |
| V2010_test | ALL      | D2_hgb_FULL_nbr3   | truth_nbr3_scored |        4572 | 16416.985 |        147454 |   0.111 |      0.534 |            0.506 |
| V2010_test | ALL      | D1_hgb_FULL_nbr3   | truth_nbr3_scored |        4572 | 15036.335 |        147454 |   0.102 |      0.538 |            0.506 |
| V2010_test | ALL      | D2_hgb_NO_PNET_own | truth_own_scored  |        4572 |  7714.036 |         30921 |   0.249 |      0.518 |            0.458 |
| V2010_test | ALL      | P1                 | truth_own_scored  |        4572 |  6031.465 |         30921 |   0.195 |      0.491 |            0.458 |
| V2010_test | NW       | D2_hgb_FULL_own    | truth_own_scored  |        2640 |  3652.344 |         17082 |   0.214 |      0.518 |            0.480 |
| V2010_test | NW       | D1_hgb_FULL_own    | truth_own_scored  |        2640 |  3059.592 |         17082 |   0.179 |      0.505 |            0.480 |
| V2010_test | NW       | D2_hgb_FULL_nbr3   | truth_nbr3_scored |        2640 |  3598.981 |         37188 |   0.097 |      0.443 |            0.403 |
| V2010_test | NW       | D1_hgb_FULL_nbr3   | truth_nbr3_scored |        2640 |  3024.091 |         37188 |   0.081 |      0.441 |            0.403 |
| V2010_test | NW       | D2_hgb_NO_PNET_own | truth_own_scored  |        2640 |  3545.725 |         17082 |   0.208 |      0.508 |            0.480 |
| V2010_test | NW       | P1                 | truth_own_scored  |        2640 |  2471.141 |         17082 |   0.145 |      0.514 |            0.480 |
| V2010_test | Olympic  | D2_hgb_FULL_own    | truth_own_scored  |         784 |  4049.136 |          7056 |   0.574 |      0.749 |            0.281 |
| V2010_test | Olympic  | D1_hgb_FULL_own    | truth_own_scored  |         784 |  2701.007 |          7056 |   0.383 |      0.797 |            0.281 |
| V2010_test | Olympic  | D2_hgb_FULL_nbr3   | truth_nbr3_scored |         784 | 12715.755 |         85199 |   0.149 |      0.714 |            0.386 |
| V2010_test | Olympic  | D1_hgb_FULL_nbr3   | truth_nbr3_scored |         784 | 11884.395 |         85199 |   0.139 |      0.803 |            0.386 |
| V2010_test | Olympic  | D2_hgb_NO_PNET_own | truth_own_scored  |         784 |  4077.594 |          7056 |   0.578 |      0.752 |            0.281 |
| V2010_test | Olympic  | P1                 | truth_own_scored  |         784 |  3481.801 |          7056 |   0.493 |      0.373 |            0.281 |
| V2010_test | Arm B    | D2_hgb_FULL_own    | truth_own_scored  |        1148 |    89.514 |          6783 |   0.013 |      0.195 |            0.189 |
| V2010_test | Arm B    | D1_hgb_FULL_own    | truth_own_scored  |        1148 |   115.058 |          6783 |   0.017 |      0.194 |            0.189 |
| V2010_test | Arm B    | D2_hgb_FULL_nbr3   | truth_nbr3_scored |        1148 |   102.249 |         25067 |   0.004 |      0.082 |            0.076 |
| V2010_test | Arm B    | D1_hgb_FULL_nbr3   | truth_nbr3_scored |        1148 |   127.849 |         25067 |   0.005 |      0.082 |            0.076 |
| V2010_test | Arm B    | D2_hgb_NO_PNET_own | truth_own_scored  |        1148 |    90.717 |          6783 |   0.013 |      0.194 |            0.189 |
| V2010_test | Arm B    | P1                 | truth_own_scored  |        1148 |    78.523 |          6783 |   0.012 |      0.180 |            0.189 |
| V2008_test | ALL      | D2_hgb_FULL_own    | truth_own_scored  |        4934 |  3892.114 |         23962 |   0.162 |      0.455 |            0.428 |
| V2008_test | ALL      | D1_hgb_FULL_own    | truth_own_scored  |        4934 |  3347.190 |         23962 |   0.140 |      0.443 |            0.428 |
| V2008_test | ALL      | D2_hgb_FULL_nbr3   | truth_nbr3_scored |        4934 |  3896.783 |         44153 |   0.088 |      0.392 |            0.363 |
| V2008_test | ALL      | D1_hgb_FULL_nbr3   | truth_nbr3_scored |        4934 |  3369.527 |         44153 |   0.076 |      0.386 |            0.363 |
| V2008_test | ALL      | D2_hgb_NO_PNET_own | truth_own_scored  |        4934 |  3910.877 |         23962 |   0.163 |      0.451 |            0.428 |
| V2008_test | ALL      | P1                 | truth_own_scored  |        4934 |  2494.646 |         23962 |   0.104 |      0.452 |            0.428 |
| V2008_test | NW       | D2_hgb_FULL_own    | truth_own_scored  |        3036 |  3751.899 |         16129 |   0.233 |      0.537 |            0.494 |
| V2008_test | NW       | D1_hgb_FULL_own    | truth_own_scored  |        3036 |  3211.437 |         16129 |   0.199 |      0.517 |            0.494 |
| V2008_test | NW       | D2_hgb_FULL_nbr3   | truth_nbr3_scored |        3036 |  3761.760 |         25449 |   0.148 |      0.495 |            0.449 |
| V2008_test | NW       | D1_hgb_FULL_nbr3   | truth_nbr3_scored |        3036 |  3233.756 |         25449 |   0.127 |      0.487 |            0.449 |
| V2008_test | NW       | D2_hgb_NO_PNET_own | truth_own_scored  |        3036 |  3760.326 |         16129 |   0.233 |      0.529 |            0.494 |
| V2008_test | NW       | P1                 | truth_own_scored  |        3036 |  2333.408 |         16129 |   0.145 |      0.519 |            0.494 |
| V2008_test | SW       | D2_hgb_FULL_own    | truth_own_scored  |        1898 |   140.215 |          7833 |   0.018 |      0.259 |            0.251 |
| V2008_test | SW       | D1_hgb_FULL_own    | truth_own_scored  |        1898 |   135.753 |          7833 |   0.017 |      0.254 |            0.251 |
| V2008_test | SW       | D2_hgb_FULL_nbr3   | truth_nbr3_scored |        1898 |   135.024 |         18704 |   0.007 |      0.216 |            0.208 |
| V2008_test | SW       | D1_hgb_FULL_nbr3   | truth_nbr3_scored |        1898 |   135.771 |         18704 |   0.007 |      0.212 |            0.208 |
| V2008_test | SW       | D2_hgb_NO_PNET_own | truth_own_scored  |        1898 |   150.551 |          7833 |   0.019 |      0.258 |            0.251 |
| V2008_test | SW       | P1                 | truth_own_scored  |        1898 |   161.238 |          7833 |   0.021 |      0.278 |            0.251 |
| A2011      | ALL      | D2_hgb_FULL_own    | truth_own_scored  |        3828 | 11619.814 |         33306 |   0.349 |      0.473 |            0.404 |
| A2011      | ALL      | D1_hgb_FULL_own    | truth_own_scored  |        3828 | 10857.238 |         33306 |   0.326 |      0.479 |            0.404 |
| A2011      | ALL      | D2_hgb_FULL_nbr3   | truth_nbr3_scored |        3828 | 14004.321 |         50539 |   0.277 |      0.454 |            0.383 |
| A2011      | ALL      | D1_hgb_FULL_nbr3   | truth_nbr3_scored |        3828 | 13169.157 |         50539 |   0.261 |      0.469 |            0.383 |
| A2011      | ALL      | D2_hgb_NO_PNET_own | truth_own_scored  |        3828 | 11206.639 |         33306 |   0.336 |      0.466 |            0.404 |
| A2011      | ALL      | P1                 | truth_own_scored  |        3828 |  8464.658 |         33306 |   0.254 |      0.476 |            0.404 |
| A2011      | NW       | D2_hgb_FULL_own    | truth_own_scored  |        3828 | 11619.814 |         33306 |   0.349 |      0.473 |            0.404 |
| A2011      | NW       | D1_hgb_FULL_own    | truth_own_scored  |        3828 | 10857.238 |         33306 |   0.326 |      0.479 |            0.404 |
| A2011      | NW       | D2_hgb_FULL_nbr3   | truth_nbr3_scored |        3828 | 14004.321 |         50539 |   0.277 |      0.454 |            0.383 |
| A2011      | NW       | D1_hgb_FULL_nbr3   | truth_nbr3_scored |        3828 | 13169.157 |         50539 |   0.261 |      0.469 |            0.383 |
| A2011      | NW       | D2_hgb_NO_PNET_own | truth_own_scored  |        3828 | 11206.639 |         33306 |   0.336 |      0.466 |            0.404 |
| A2011      | NW       | P1                 | truth_own_scored  |        3828 |  8464.658 |         33306 |   0.254 |      0.476 |            0.404 |

## (f) Permutation importance — D2_hgb_FULL_own (rows: A2011, n = 61999)

Blocks (features permuted jointly):

| feature            |   auc_drop |   logloss_gain |
|:-------------------|-----------:|---------------:|
| PNET (all)         |     0.1689 |         0.1137 |
| PNET fine-tuned    |     0.1578 |         0.1001 |
| step-1 P           |     0.0219 |         0.0184 |
| step-2 attribution |     0.0027 |         0.0015 |
| PNET off-the-shelf |     0.0012 |        -0.0012 |

Top 15 single features:

| feature          |   auc_drop |   logloss_gain |
|:-----------------|-----------:|---------------:|
| pnet_max_f       |     0.0947 |         0.0585 |
| site_cc_max      |     0.0239 |         0.0159 |
| P1               |     0.0220 |         0.0185 |
| net_mad          |     0.0104 |         0.0099 |
| pnet_2nd_f       |     0.0103 |         0.0049 |
| dens_own_1h      |     0.0063 |         0.0061 |
| dens_own_10m     |     0.0048 |         0.0060 |
| dens_grp_1h      |     0.0041 |         0.0024 |
| site_z_2nd       |     0.0028 |         0.0020 |
| pnet_2nd_w       |     0.0022 |        -0.0007 |
| grp_size         |     0.0019 |         0.0006 |
| site_z_max       |     0.0018 |         0.0001 |
| dom_top          |     0.0017 |         0.0010 |
| n_sites_pnet03_f |     0.0017 |         0.0008 |
| hfr_rel          |     0.0013 |         0.0011 |

## (f) Permutation importance — D2_hgb_FULL_any (rows: A2011, n = 61999)

Blocks (features permuted jointly):

| feature            |   auc_drop |   logloss_gain |
|:-------------------|-----------:|---------------:|
| PNET fine-tuned    |     0.2013 |         0.1477 |
| PNET (all)         |     0.1908 |         0.1383 |
| step-1 P           |     0.0087 |         0.0095 |
| step-2 attribution |     0.0017 |         0.0016 |
| PNET off-the-shelf |    -0.0023 |        -0.0023 |

Top 15 single features:

| feature          |   auc_drop |   logloss_gain |
|:-----------------|-----------:|---------------:|
| pnet_max_f       |     0.1335 |         0.0929 |
| pnet_2nd_f       |     0.0209 |         0.0132 |
| dens_own_1h      |     0.0171 |         0.0115 |
| hfr_rel          |     0.0106 |         0.0055 |
| P1               |     0.0083 |         0.0098 |
| net_mad          |     0.0071 |         0.0067 |
| site_cc_max      |     0.0038 |        -0.0011 |
| n_sites_pnet03_f |     0.0027 |         0.0012 |
| dom_top          |     0.0016 |         0.0014 |
| wsum             |     0.0015 |         0.0026 |
| site_z_max       |     0.0014 |        -0.0009 |
| pnet_rel_max_f   |     0.0012 |         0.0004 |
| site_z_2nd       |     0.0009 |         0.0005 |
| n_elig           |     0.0009 |         0.0015 |
| grp_size         |     0.0007 |         0.0002 |

## (f) Permutation importance — D1_hgb_FULL_own (rows: V2010_test+V2008_test+A2011, n = 150000)

Blocks (features permuted jointly):

| feature            |   auc_drop |   logloss_gain |
|:-------------------|-----------:|---------------:|
| PNET (all)         |     0.0518 |         0.0445 |
| PNET fine-tuned    |     0.0512 |         0.0439 |
| step-1 P           |     0.0194 |         0.0153 |
| step-2 attribution |     0.0039 |         0.0023 |
| PNET off-the-shelf |     0.0009 |         0.0006 |

Top 15 single features:

| feature     |   auc_drop |   logloss_gain |
|:------------|-----------:|---------------:|
| region_code |     0.0489 |         0.0318 |
| pnet_max_f  |     0.0479 |         0.0405 |
| net_cc      |     0.0352 |         0.0354 |
| dom_top     |     0.0338 |         0.0231 |
| wsum        |     0.0285 |         0.0306 |
| P1          |     0.0186 |         0.0141 |
| net_mad     |     0.0080 |         0.0122 |
| dens_own_1h |     0.0036 |         0.0017 |
| n_z2        |     0.0022 |        -0.0021 |
| pnet_2nd_f  |     0.0020 |         0.0012 |
| dens_grp_1h |     0.0019 |         0.0000 |
| n_sta       |     0.0018 |         0.0024 |
| P_member    |     0.0018 |         0.0011 |
| site_cc_max |     0.0016 |         0.0000 |
| grp_size    |     0.0014 |         0.0004 |

## (f) Permutation importance — D1_hgb_FULL_any (rows: V2010_test+V2008_test+A2011, n = 150000)

Blocks (features permuted jointly):

| feature            |   auc_drop |   logloss_gain |
|:-------------------|-----------:|---------------:|
| PNET (all)         |     0.0865 |         0.0628 |
| PNET fine-tuned    |     0.0861 |         0.0631 |
| step-1 P           |     0.0192 |         0.0300 |
| step-2 attribution |     0.0098 |         0.0076 |
| PNET off-the-shelf |     0.0005 |         0.0003 |

Top 15 single features:

| feature     |   auc_drop |   logloss_gain |
|:------------|-----------:|---------------:|
| pnet_max_f  |     0.0711 |         0.0512 |
| P1          |     0.0198 |         0.0304 |
| n_sta       |     0.0163 |         0.0124 |
| dens_own_1h |     0.0146 |         0.0091 |
| hfr_rel     |     0.0144 |         0.0110 |
| net_mad     |     0.0118 |         0.0129 |
| wsum        |     0.0113 |         0.0083 |
| site_cc_max |     0.0100 |         0.0064 |
| dom_top     |     0.0078 |         0.0067 |
| net_cc      |     0.0061 |         0.0072 |
| pnet_2nd_f  |     0.0056 |         0.0040 |
| n_sites     |     0.0048 |         0.0034 |
| P_member    |     0.0041 |         0.0032 |
| dens_grp_1h |     0.0038 |         0.0038 |
| region_code |     0.0035 |         0.0021 |

## D3 — scrambled-control candidates added to the fit as known negatives (`scripts/confidence_d3.py`)

Same as D2 with the NO_STEP2 features, but c1–c3 candidates of the same family-days enter the fit with label 0 (weighted like the data rows); isotonic calibration on data rows only; c4 never enters any fit. `withinFD` = AUC computed inside each family-day and averaged (pairs-weighted): does P rank the repeats of one family on one day, beyond knowing which family-days are busy?

| eval_set   | region   | label   |      n |   n_pos |   auc_D3_hgb_NO_STEP2 |   auc_D2_hgb_NO_STEP2 |   ap_D3_hgb_NO_STEP2 |   ap_D2_hgb_NO_STEP2 |   slope_D3_hgb_NO_STEP2 |   slope_D2_hgb_NO_STEP2 |   obs_exp_D3_hgb_NO_STEP2 |   obs_exp_D2_hgb_NO_STEP2 |   auc_withinFD_D3_hgb_NO_STEP2 |   auc_withinFD_D2_hgb_NO_STEP2 |   auc_withinFD_P1 |   auc_withinFD_pnet_max_f |
|:-----------|:---------|:--------|-------:|--------:|----------------------:|----------------------:|---------------------:|---------------------:|------------------------:|------------------------:|--------------------------:|--------------------------:|-------------------------------:|-------------------------------:|------------------:|--------------------------:|
| V2010_test | ALL      | own     | 115301 |    7761 |                 0.916 |                 0.918 |                0.535 |                0.547 |                   0.999 |                   1.006 |                     0.985 |                     0.994 |                          0.825 |                          0.827 |             0.511 |                     0.722 |
| V2010_test | ALL      | any     | 115301 |   30626 |                 0.767 |                 0.770 |                0.571 |                0.587 |                   1.044 |                   1.042 |                     1.003 |                     1.009 |                          0.705 |                          0.708 |             0.494 |                     0.641 |
| V2010_test | NW       | own     |  12527 |    3662 |                 0.846 |                 0.852 |                0.715 |                0.728 |                   0.978 |                   0.994 |                     0.984 |                     0.999 |                          0.725 |                          0.736 |             0.511 |                     0.677 |
| V2010_test | NW       | any     |  12527 |    5903 |                 0.771 |                 0.779 |                0.763 |                0.776 |                   0.927 |                   0.948 |                     1.021 |                     1.043 |                          0.675 |                          0.685 |             0.487 |                     0.607 |
| V2010_test | Olympic  | own     | 101356 |    4013 |                 0.897 |                 0.899 |                0.314 |                0.326 |                   1.054 |                   1.023 |                     0.990 |                     0.990 |                          0.846 |                          0.846 |             0.511 |                     0.731 |
| V2010_test | Olympic  | any     | 101356 |   24225 |                 0.755 |                 0.757 |                0.510 |                0.526 |                   1.089 |                   1.067 |                     0.999 |                     1.002 |                          0.707 |                          0.710 |             0.495 |                     0.643 |
| V2010_test | Arm B    | own     |   1418 |      86 |                 0.808 |                 0.820 |                0.424 |                0.419 |                   0.982 |                   1.502 |                     0.807 |                     0.967 |                          0.541 |                          0.577 |             0.571 |                     0.429 |
| V2010_test | Arm B    | any     |   1418 |     498 |                 0.643 |                 0.657 |                0.502 |                0.527 |                   0.637 |                   0.744 |                     0.957 |                     0.968 |                          0.557 |                          0.564 |             0.485 |                     0.507 |
| V2008_test | ALL      | own     |  13943 |    3857 |                 0.865 |                 0.869 |                0.716 |                0.731 |                   1.025 |                   1.029 |                     0.993 |                     0.987 |                          0.721 |                          0.732 |             0.497 |                     0.694 |
| V2008_test | ALL      | any     |  13943 |    5620 |                 0.823 |                 0.829 |                0.765 |                0.782 |                   1.048 |                   1.062 |                     0.954 |                     0.974 |                          0.643 |                          0.662 |             0.486 |                     0.641 |
| V2008_test | NW       | own     |  13647 |    3708 |                 0.862 |                 0.867 |                0.707 |                0.723 |                   1.022 |                   1.026 |                     0.991 |                     0.985 |                          0.721 |                          0.732 |             0.496 |                     0.694 |
| V2008_test | NW       | any     |  13647 |    5451 |                 0.821 |                 0.827 |                0.760 |                0.778 |                   1.044 |                   1.057 |                     0.953 |                     0.971 |                          0.643 |                          0.662 |             0.486 |                     0.641 |
| V2008_test | SW       | own     |    296 |     149 |                 0.889 |                 0.894 |                0.879 |                0.889 |                   1.054 |                   1.087 |                     1.033 |                     1.049 |                          0.734 |                          0.755 |             0.626 |                     0.612 |
| V2008_test | SW       | any     |    296 |     169 |                 0.894 |                 0.892 |                0.911 |                0.912 |                   1.163 |                   1.186 |                     0.989 |                     1.094 |                          0.774 |                          0.790 |             0.613 |                     0.637 |
| A2011      | ALL      | own     |  61999 |    8865 |                 0.864 |                 0.869 |                0.582 |                0.605 |                   0.897 |                   0.929 |                     0.730 |                     0.756 |                          0.735 |                          0.742 |             0.555 |                     0.650 |
| A2011      | ALL      | any     |  61999 |   16123 |                 0.792 |                 0.803 |                0.611 |                0.645 |                   0.982 |                   0.976 |                     0.748 |                     0.719 |                          0.664 |                          0.676 |             0.521 |                     0.603 |
| A2011      | NW       | own     |  61999 |    8865 |                 0.864 |                 0.869 |                0.582 |                0.605 |                   0.897 |                   0.929 |                     0.730 |                     0.756 |                          0.735 |                          0.742 |             0.555 |                     0.650 |
| A2011      | NW       | any     |  61999 |   16123 |                 0.792 |                 0.803 |                0.611 |                0.645 |                   0.982 |                   0.976 |                     0.748 |                     0.719 |                          0.664 |                          0.676 |             0.521 |                     0.603 |

D3 vs D2 on the scrambled control (c1–c3: seen in training on OTHER days only; c4: never seen):

| eval_set   | region   | days   | model               | controls   |   data_fracP_gt05 |   ctl_fracP_gt05 |   data_sumP |   ctl_sumP_per_real |   ctl_over_data_sumP |
|:-----------|:---------|:-------|:--------------------|:-----------|------------------:|-----------------:|------------:|--------------------:|---------------------:|
| V2010_test | ALL      | all    | D3_hgb_NO_STEP2_own | c1,c2,c3   |             0.022 |            0.003 |   10549.618 |            7870.413 |                0.746 |
| V2010_test | ALL      | all    | D3_hgb_NO_STEP2_own | c4         |             0.022 |            0.003 |   10549.618 |            7833.811 |                0.743 |
| V2010_test | ALL      | all    | D3_hgb_NO_STEP2_any | c1,c2,c3   |             0.115 |            0.089 |   44545.919 |           41172.546 |                0.924 |
| V2010_test | ALL      | all    | D3_hgb_NO_STEP2_any | c4         |             0.115 |            0.090 |   44545.919 |           40953.837 |                0.919 |
| V2010_test | ALL      | all    | D2_hgb_NO_STEP2_own | c1,c2,c3   |             0.022 |            0.007 |   10398.062 |            8030.796 |                0.772 |
| V2010_test | ALL      | all    | D2_hgb_NO_STEP2_own | c4         |             0.022 |            0.007 |   10398.062 |            8001.236 |                0.769 |
| V2010_test | ALL      | all    | D2_hgb_NO_STEP2_any | c1,c2,c3   |             0.116 |            0.097 |   43998.926 |           41045.018 |                0.933 |
| V2010_test | ALL      | all    | D2_hgb_NO_STEP2_any | c4         |             0.116 |            0.098 |   43998.926 |           40906.672 |                0.930 |
| V2010_test | ALL      | all    | D2_hgb_FULL_own     | c1,c2,c3   |             0.022 |            0.006 |   10354.561 |            7620.072 |                0.736 |
| V2010_test | ALL      | all    | D2_hgb_FULL_own     | c4         |             0.022 |            0.005 |   10354.561 |            7411.506 |                0.716 |
| V2010_test | ALL      | all    | D2_hgb_FULL_any     | c1,c2,c3   |             0.122 |            0.087 |   43701.209 |           39101.452 |                0.895 |
| V2010_test | ALL      | all    | D2_hgb_FULL_any     | c4         |             0.122 |            0.080 |   43701.209 |           38174.316 |                0.874 |
| V2010_test | ALL      | quiet  | D3_hgb_NO_STEP2_own | c1,c2,c3   |             0.002 |            0.000 |    1971.523 |            1816.678 |                0.921 |
| V2010_test | ALL      | quiet  | D3_hgb_NO_STEP2_own | c4         |             0.002 |            0.000 |    1971.523 |            1847.066 |                0.937 |
| V2010_test | ALL      | quiet  | D3_hgb_NO_STEP2_any | c1,c2,c3   |             0.091 |            0.080 |   11756.776 |           11429.236 |                0.972 |
| V2010_test | ALL      | quiet  | D3_hgb_NO_STEP2_any | c4         |             0.091 |            0.082 |   11756.776 |           11416.550 |                0.971 |
| V2010_test | ALL      | quiet  | D2_hgb_NO_STEP2_own | c1,c2,c3   |             0.003 |            0.001 |    1905.313 |            1782.075 |                0.935 |
| V2010_test | ALL      | quiet  | D2_hgb_NO_STEP2_own | c4         |             0.003 |            0.002 |    1905.313 |            1809.854 |                0.950 |
| V2010_test | ALL      | quiet  | D2_hgb_NO_STEP2_any | c1,c2,c3   |             0.088 |            0.080 |   11422.783 |           11123.762 |                0.974 |
| V2010_test | ALL      | quiet  | D2_hgb_NO_STEP2_any | c4         |             0.088 |            0.083 |   11422.783 |           11152.771 |                0.976 |
| V2010_test | ALL      | quiet  | D2_hgb_FULL_own     | c1,c2,c3   |             0.003 |            0.001 |    1875.827 |            1675.167 |                0.893 |
| V2010_test | ALL      | quiet  | D2_hgb_FULL_own     | c4         |             0.003 |            0.001 |    1875.827 |            1662.421 |                0.886 |
| V2010_test | ALL      | quiet  | D2_hgb_FULL_any     | c1,c2,c3   |             0.086 |            0.067 |   11098.350 |           10449.429 |                0.942 |
| V2010_test | ALL      | quiet  | D2_hgb_FULL_any     | c4         |             0.086 |            0.065 |   11098.350 |           10300.682 |                0.928 |
| V2010_test | NW       | all    | D3_hgb_NO_STEP2_own | c1,c2,c3   |             0.151 |            0.023 |    4390.474 |            2213.775 |                0.504 |
| V2010_test | NW       | all    | D3_hgb_NO_STEP2_own | c4         |             0.151 |            0.024 |    4390.474 |            2213.092 |                0.504 |
| V2010_test | NW       | all    | D3_hgb_NO_STEP2_any | c1,c2,c3   |             0.250 |            0.111 |    7722.177 |            5210.512 |                0.675 |
| V2010_test | NW       | all    | D3_hgb_NO_STEP2_any | c4         |             0.250 |            0.109 |    7722.177 |            5228.638 |                0.677 |
| V2010_test | NW       | all    | D2_hgb_NO_STEP2_own | c1,c2,c3   |             0.138 |            0.036 |    4327.169 |            2375.031 |                0.549 |
| V2010_test | NW       | all    | D2_hgb_NO_STEP2_own | c4         |             0.138 |            0.036 |    4327.169 |            2387.639 |                0.552 |
| V2010_test | NW       | all    | D2_hgb_NO_STEP2_any | c1,c2,c3   |             0.235 |            0.136 |    7575.079 |            5360.978 |                0.708 |
| V2010_test | NW       | all    | D2_hgb_NO_STEP2_any | c4         |             0.235 |            0.136 |    7575.079 |            5403.160 |                0.713 |
| V2010_test | NW       | all    | D2_hgb_FULL_own     | c1,c2,c3   |             0.136 |            0.034 |    4320.688 |            2355.850 |                0.545 |
| V2010_test | NW       | all    | D2_hgb_FULL_own     | c4         |             0.136 |            0.034 |    4320.688 |            2366.779 |                0.548 |
| V2010_test | NW       | all    | D2_hgb_FULL_any     | c1,c2,c3   |             0.235 |            0.132 |    7540.746 |            5304.745 |                0.703 |
| V2010_test | NW       | all    | D2_hgb_FULL_any     | c4         |             0.235 |            0.133 |    7540.746 |            5334.180 |                0.707 |
| V2010_test | NW       | quiet  | D3_hgb_NO_STEP2_own | c1,c2,c3   |             0.001 |            0.000 |     324.885 |             320.434 |                0.986 |
| V2010_test | NW       | quiet  | D3_hgb_NO_STEP2_own | c4         |             0.001 |            0.000 |     324.885 |             326.615 |                1.005 |
| V2010_test | NW       | quiet  | D3_hgb_NO_STEP2_any | c1,c2,c3   |             0.015 |            0.011 |    1067.010 |            1060.741 |                0.994 |
| V2010_test | NW       | quiet  | D3_hgb_NO_STEP2_any | c4         |             0.015 |            0.014 |    1067.010 |            1074.510 |                1.007 |
| V2010_test | NW       | quiet  | D2_hgb_NO_STEP2_own | c1,c2,c3   |             0.001 |            0.001 |     326.029 |             323.011 |                0.991 |
| V2010_test | NW       | quiet  | D2_hgb_NO_STEP2_own | c4         |             0.001 |            0.001 |     326.029 |             330.185 |                1.013 |
| V2010_test | NW       | quiet  | D2_hgb_NO_STEP2_any | c1,c2,c3   |             0.021 |            0.017 |    1032.848 |            1020.683 |                0.988 |
| V2010_test | NW       | quiet  | D2_hgb_NO_STEP2_any | c4         |             0.021 |            0.020 |    1032.848 |            1038.775 |                1.006 |
| V2010_test | NW       | quiet  | D2_hgb_FULL_own     | c1,c2,c3   |             0.001 |            0.001 |     326.191 |             321.170 |                0.985 |
| V2010_test | NW       | quiet  | D2_hgb_FULL_own     | c4         |             0.001 |            0.001 |     326.191 |             327.400 |                1.004 |
| V2010_test | NW       | quiet  | D2_hgb_FULL_any     | c1,c2,c3   |             0.020 |            0.015 |    1029.915 |            1017.687 |                0.988 |
| V2010_test | NW       | quiet  | D2_hgb_FULL_any     | c4         |             0.020 |            0.018 |    1029.915 |            1032.886 |                1.003 |
| V2010_test | Olympic  | all    | D3_hgb_NO_STEP2_own | c1,c2,c3   |             0.002 |            0.001 |    5875.633 |            5554.342 |                0.945 |
| V2010_test | Olympic  | all    | D3_hgb_NO_STEP2_own | c4         |             0.002 |            0.001 |    5875.633 |            5521.641 |                0.940 |
| V2010_test | Olympic  | all    | D3_hgb_NO_STEP2_any | c1,c2,c3   |             0.093 |            0.086 |   35977.295 |           35285.963 |                0.981 |
| V2010_test | Olympic  | all    | D3_hgb_NO_STEP2_any | c4         |             0.093 |            0.087 |   35977.295 |           35063.860 |                0.975 |
| V2010_test | Olympic  | all    | D2_hgb_NO_STEP2_own | c1,c2,c3   |             0.004 |            0.003 |    5805.939 |            5560.947 |                0.958 |
| V2010_test | Olympic  | all    | D2_hgb_NO_STEP2_own | c4         |             0.004 |            0.003 |    5805.939 |            5518.948 |                0.951 |
| V2010_test | Olympic  | all    | D2_hgb_NO_STEP2_any | c1,c2,c3   |             0.096 |            0.092 |   35588.793 |           35007.703 |                0.984 |
| V2010_test | Olympic  | all    | D2_hgb_NO_STEP2_any | c4         |             0.096 |            0.092 |   35588.793 |           34837.368 |                0.979 |
| V2010_test | Olympic  | all    | D2_hgb_FULL_own     | c1,c2,c3   |             0.004 |            0.002 |    5769.231 |            5171.844 |                0.896 |
| V2010_test | Olympic  | all    | D2_hgb_FULL_own     | c4         |             0.004 |            0.001 |    5769.231 |            4955.041 |                0.859 |
| V2010_test | Olympic  | all    | D2_hgb_FULL_any     | c1,c2,c3   |             0.103 |            0.080 |   35324.831 |           33122.026 |                0.938 |
| V2010_test | Olympic  | all    | D2_hgb_FULL_any     | c4         |             0.103 |            0.072 |   35324.831 |           32176.067 |                0.911 |
| V2010_test | Olympic  | quiet  | D3_hgb_NO_STEP2_own | c1,c2,c3   |             0.002 |            0.000 |    1631.964 |            1481.480 |                0.908 |
| V2010_test | Olympic  | quiet  | D3_hgb_NO_STEP2_own | c4         |             0.002 |            0.000 |    1631.964 |            1504.073 |                0.922 |
| V2010_test | Olympic  | quiet  | D3_hgb_NO_STEP2_any | c1,c2,c3   |             0.100 |            0.089 |   10562.549 |           10226.227 |                0.968 |
| V2010_test | Olympic  | quiet  | D3_hgb_NO_STEP2_any | c4         |             0.100 |            0.090 |   10562.549 |           10199.934 |                0.966 |
| V2010_test | Olympic  | quiet  | D2_hgb_NO_STEP2_own | c1,c2,c3   |             0.003 |            0.001 |    1566.208 |            1445.254 |                0.923 |
| V2010_test | Olympic  | quiet  | D2_hgb_NO_STEP2_own | c4         |             0.003 |            0.002 |    1566.208 |            1464.593 |                0.935 |
| V2010_test | Olympic  | quiet  | D2_hgb_NO_STEP2_any | c1,c2,c3   |             0.096 |            0.088 |   10265.388 |            9964.568 |                0.971 |
| V2010_test | Olympic  | quiet  | D2_hgb_NO_STEP2_any | c4         |             0.096 |            0.091 |   10265.388 |            9974.506 |                0.972 |
| V2010_test | Olympic  | quiet  | D2_hgb_FULL_own     | c1,c2,c3   |             0.003 |            0.001 |    1536.869 |            1341.112 |                0.873 |
| V2010_test | Olympic  | quiet  | D2_hgb_FULL_own     | c4         |             0.003 |            0.001 |    1536.869 |            1321.038 |                0.860 |
| V2010_test | Olympic  | quiet  | D2_hgb_FULL_any     | c1,c2,c3   |             0.094 |            0.073 |    9945.503 |            9294.802 |                0.935 |
| V2010_test | Olympic  | quiet  | D2_hgb_FULL_any     | c4         |             0.094 |            0.070 |    9945.503 |            9130.617 |                0.918 |
| V2008_test | ALL      | all    | D3_hgb_NO_STEP2_own | c1,c2,c3   |             0.035 |            0.004 |    7415.499 |            5194.154 |                0.700 |
| V2008_test | ALL      | all    | D3_hgb_NO_STEP2_own | c4         |             0.035 |            0.003 |    7415.499 |            5260.127 |                0.709 |
| V2008_test | ALL      | all    | D3_hgb_NO_STEP2_any | c1,c2,c3   |             0.132 |            0.101 |   26462.467 |           24158.930 |                0.913 |
| V2008_test | ALL      | all    | D3_hgb_NO_STEP2_any | c4         |             0.132 |            0.101 |   26462.467 |           24515.818 |                0.926 |
| V2008_test | ALL      | all    | D2_hgb_NO_STEP2_own | c1,c2,c3   |             0.034 |            0.008 |    7165.579 |            5169.167 |                0.721 |
| V2008_test | ALL      | all    | D2_hgb_NO_STEP2_own | c4         |             0.034 |            0.008 |    7165.579 |            5218.536 |                0.728 |
| V2008_test | ALL      | all    | D2_hgb_NO_STEP2_any | c1,c2,c3   |             0.160 |            0.139 |   27773.352 |           25779.373 |                0.928 |
| V2008_test | ALL      | all    | D2_hgb_NO_STEP2_any | c4         |             0.160 |            0.138 |   27773.352 |           26153.844 |                0.942 |
| V2008_test | ALL      | all    | D2_hgb_FULL_own     | c1,c2,c3   |             0.035 |            0.008 |    7286.536 |            4995.657 |                0.686 |
| V2008_test | ALL      | all    | D2_hgb_FULL_own     | c4         |             0.035 |            0.008 |    7286.536 |            4902.604 |                0.673 |
| V2008_test | ALL      | all    | D2_hgb_FULL_any     | c1,c2,c3   |             0.177 |            0.127 |   28294.442 |           24921.957 |                0.881 |
| V2008_test | ALL      | all    | D2_hgb_FULL_any     | c4         |             0.177 |            0.113 |   28294.442 |           24629.901 |                0.870 |
| V2008_test | ALL      | quiet  | D3_hgb_NO_STEP2_own | c1,c2,c3   |             0.000 |            0.000 |     899.472 |             864.144 |                0.961 |
| V2008_test | ALL      | quiet  | D3_hgb_NO_STEP2_own | c4         |             0.000 |            0.000 |     899.472 |             866.658 |                0.964 |
| V2008_test | ALL      | quiet  | D3_hgb_NO_STEP2_any | c1,c2,c3   |             0.182 |            0.207 |    4238.241 |            4128.641 |                0.974 |
| V2008_test | ALL      | quiet  | D3_hgb_NO_STEP2_any | c4         |             0.182 |            0.197 |    4238.241 |            4153.954 |                0.980 |
| V2008_test | ALL      | quiet  | D2_hgb_NO_STEP2_own | c1,c2,c3   |             0.001 |            0.001 |     794.678 |             769.256 |                0.968 |
| V2008_test | ALL      | quiet  | D2_hgb_NO_STEP2_own | c4         |             0.001 |            0.001 |     794.678 |             762.488 |                0.959 |
| V2008_test | ALL      | quiet  | D2_hgb_NO_STEP2_any | c1,c2,c3   |             0.191 |            0.211 |    4347.161 |            4224.067 |                0.972 |
| V2008_test | ALL      | quiet  | D2_hgb_NO_STEP2_any | c4         |             0.191 |            0.204 |    4347.161 |            4242.750 |                0.976 |
| V2008_test | ALL      | quiet  | D2_hgb_FULL_own     | c1,c2,c3   |             0.002 |            0.001 |     785.839 |             714.597 |                0.909 |
| V2008_test | ALL      | quiet  | D2_hgb_FULL_own     | c4         |             0.002 |            0.000 |     785.839 |             684.978 |                0.872 |
| V2008_test | ALL      | quiet  | D2_hgb_FULL_any     | c1,c2,c3   |             0.223 |            0.208 |    4526.172 |            4169.100 |                0.921 |
| V2008_test | ALL      | quiet  | D2_hgb_FULL_any     | c4         |             0.223 |            0.183 |    4526.172 |            4077.698 |                0.901 |
| V2008_test | NW       | all    | D3_hgb_NO_STEP2_own | c1,c2,c3   |             0.187 |            0.023 |    3958.820 |            1875.026 |                0.474 |
| V2008_test | NW       | all    | D3_hgb_NO_STEP2_own | c4         |             0.187 |            0.022 |    3958.820 |            1870.650 |                0.473 |
| V2008_test | NW       | all    | D3_hgb_NO_STEP2_any | c1,c2,c3   |             0.302 |            0.119 |    6263.404 |            3815.519 |                0.609 |
| V2008_test | NW       | all    | D3_hgb_NO_STEP2_any | c4         |             0.302 |            0.117 |    6263.404 |            3824.558 |                0.611 |
| V2008_test | NW       | all    | D2_hgb_NO_STEP2_own | c1,c2,c3   |             0.181 |            0.053 |    3995.131 |            2110.123 |                0.528 |
| V2008_test | NW       | all    | D2_hgb_NO_STEP2_own | c4         |             0.181 |            0.054 |    3995.131 |            2106.540 |                0.527 |
| V2008_test | NW       | all    | D2_hgb_NO_STEP2_any | c1,c2,c3   |             0.291 |            0.162 |    6168.510 |            3996.260 |                0.648 |
| V2008_test | NW       | all    | D2_hgb_NO_STEP2_any | c4         |             0.291 |            0.155 |    6168.510 |            3992.201 |                0.647 |
| V2008_test | NW       | all    | D2_hgb_FULL_own     | c1,c2,c3   |             0.181 |            0.049 |    3977.435 |            2068.189 |                0.520 |
| V2008_test | NW       | all    | D2_hgb_FULL_own     | c4         |             0.181 |            0.050 |    3977.435 |            2059.654 |                0.518 |
| V2008_test | NW       | all    | D2_hgb_FULL_any     | c1,c2,c3   |             0.294 |            0.160 |    6168.457 |            3976.499 |                0.645 |
| V2008_test | NW       | all    | D2_hgb_FULL_any     | c4         |             0.294 |            0.153 |    6168.457 |            3970.510 |                0.644 |
| V2008_test | NW       | quiet  | D3_hgb_NO_STEP2_own | c1,c2,c3   |             0.002 |            0.001 |     148.951 |             146.000 |                0.980 |
| V2008_test | NW       | quiet  | D3_hgb_NO_STEP2_own | c4         |             0.002 |            0.000 |     148.951 |             139.950 |                0.940 |
| V2008_test | NW       | quiet  | D3_hgb_NO_STEP2_any | c1,c2,c3   |             0.020 |            0.012 |     397.480 |             395.100 |                0.994 |
| V2008_test | NW       | quiet  | D3_hgb_NO_STEP2_any | c4         |             0.020 |            0.014 |     397.480 |             391.764 |                0.986 |
| V2008_test | NW       | quiet  | D2_hgb_NO_STEP2_own | c1,c2,c3   |             0.005 |            0.004 |     154.634 |             152.164 |                0.984 |
| V2008_test | NW       | quiet  | D2_hgb_NO_STEP2_own | c4         |             0.005 |            0.003 |     154.634 |             143.614 |                0.929 |
| V2008_test | NW       | quiet  | D2_hgb_NO_STEP2_any | c1,c2,c3   |             0.031 |            0.026 |     400.581 |             402.480 |                1.005 |
| V2008_test | NW       | quiet  | D2_hgb_NO_STEP2_any | c4         |             0.031 |            0.023 |     400.581 |             394.681 |                0.985 |
| V2008_test | NW       | quiet  | D2_hgb_FULL_own     | c1,c2,c3   |             0.005 |            0.003 |     151.135 |             147.955 |                0.979 |
| V2008_test | NW       | quiet  | D2_hgb_FULL_own     | c4         |             0.005 |            0.001 |     151.135 |             139.230 |                0.921 |
| V2008_test | NW       | quiet  | D2_hgb_FULL_any     | c1,c2,c3   |             0.032 |            0.027 |     403.474 |             403.620 |                1.000 |
| V2008_test | NW       | quiet  | D2_hgb_FULL_any     | c4         |             0.032 |            0.023 |     403.474 |             394.246 |                0.977 |
| V2008_test | Olympic  | all    | D3_hgb_NO_STEP2_own | c1,c2,c3   |             0.001 |            0.000 |    3200.372 |            3202.592 |                1.001 |
| V2008_test | Olympic  | all    | D3_hgb_NO_STEP2_own | c4         |             0.001 |            0.000 |    3200.372 |            3270.978 |                1.022 |
| V2008_test | Olympic  | all    | D3_hgb_NO_STEP2_any | c1,c2,c3   |             0.097 |            0.101 |   19119.842 |           19430.623 |                1.016 |
| V2008_test | Olympic  | all    | D3_hgb_NO_STEP2_any | c4         |             0.097 |            0.101 |   19119.842 |           19760.373 |                1.034 |
| V2008_test | Olympic  | all    | D2_hgb_NO_STEP2_own | c1,c2,c3   |             0.001 |            0.000 |    2929.304 |            2950.607 |                1.007 |
| V2008_test | Olympic  | all    | D2_hgb_NO_STEP2_own | c4         |             0.001 |            0.000 |    2929.304 |            3001.046 |                1.024 |
| V2008_test | Olympic  | all    | D2_hgb_NO_STEP2_any | c1,c2,c3   |             0.133 |            0.137 |   20391.608 |           20712.316 |                1.016 |
| V2008_test | Olympic  | all    | D2_hgb_NO_STEP2_any | c4         |             0.133 |            0.137 |   20391.608 |           21069.988 |                1.033 |
| V2008_test | Olympic  | all    | D2_hgb_FULL_own     | c1,c2,c3   |             0.002 |            0.001 |    3056.928 |            2817.606 |                0.922 |
| V2008_test | Olympic  | all    | D2_hgb_FULL_own     | c4         |             0.002 |            0.000 |    3056.928 |            2735.270 |                0.895 |
| V2008_test | Olympic  | all    | D2_hgb_FULL_any     | c1,c2,c3   |             0.152 |            0.123 |   20871.565 |           19852.758 |                0.951 |
| V2008_test | Olympic  | all    | D2_hgb_FULL_any     | c4         |             0.152 |            0.106 |   20871.565 |           19553.981 |                0.937 |
| V2008_test | Olympic  | quiet  | D3_hgb_NO_STEP2_own | c1,c2,c3   |             0.000 |            0.000 |     737.168 |             705.580 |                0.957 |
| V2008_test | Olympic  | quiet  | D3_hgb_NO_STEP2_own | c4         |             0.000 |            0.000 |     737.168 |             712.952 |                0.967 |
| V2008_test | Olympic  | quiet  | D3_hgb_NO_STEP2_any | c1,c2,c3   |             0.216 |            0.249 |    3708.459 |            3608.116 |                0.973 |
| V2008_test | Olympic  | quiet  | D3_hgb_NO_STEP2_any | c4         |             0.216 |            0.236 |    3708.459 |            3626.133 |                0.978 |
| V2008_test | Olympic  | quiet  | D2_hgb_NO_STEP2_own | c1,c2,c3   |             0.000 |            0.000 |     629.986 |             607.574 |                0.964 |
| V2008_test | Olympic  | quiet  | D2_hgb_NO_STEP2_own | c4         |             0.000 |            0.000 |     629.986 |             608.596 |                0.966 |
| V2008_test | Olympic  | quiet  | D2_hgb_NO_STEP2_any | c1,c2,c3   |             0.224 |            0.251 |    3808.152 |            3690.627 |                0.969 |
| V2008_test | Olympic  | quiet  | D2_hgb_NO_STEP2_any | c4         |             0.224 |            0.242 |    3808.152 |            3705.544 |                0.973 |
| V2008_test | Olympic  | quiet  | D2_hgb_FULL_own     | c1,c2,c3   |             0.001 |            0.001 |     623.570 |             557.302 |                0.894 |
| V2008_test | Olympic  | quiet  | D2_hgb_FULL_own     | c4         |             0.001 |            0.000 |     623.570 |             535.938 |                0.859 |
| V2008_test | Olympic  | quiet  | D2_hgb_FULL_any     | c1,c2,c3   |             0.262 |            0.247 |    3982.374 |            3635.113 |                0.913 |
| V2008_test | Olympic  | quiet  | D2_hgb_FULL_any     | c4         |             0.262 |            0.218 |    3982.374 |            3544.046 |                0.890 |
| A2011      | ALL      | all    | D3_hgb_NO_STEP2_own | c1,c2,c3   |             0.030 |            0.007 |   28901.251 |           23770.184 |                0.822 |
| A2011      | ALL      | all    | D3_hgb_NO_STEP2_own | c4         |             0.030 |            0.007 |   28901.251 |           23763.714 |                0.822 |
| A2011      | ALL      | all    | D3_hgb_NO_STEP2_any | c1,c2,c3   |             0.141 |            0.113 |   83441.359 |           77292.663 |                0.926 |
| A2011      | ALL      | all    | D3_hgb_NO_STEP2_any | c4         |             0.141 |            0.113 |   83441.359 |           77269.363 |                0.926 |
| A2011      | ALL      | all    | D2_hgb_NO_STEP2_own | c1,c2,c3   |             0.024 |            0.007 |   24477.291 |           20182.068 |                0.825 |
| A2011      | ALL      | all    | D2_hgb_NO_STEP2_own | c4         |             0.024 |            0.008 |   24477.291 |           20184.328 |                0.825 |
| A2011      | ALL      | all    | D2_hgb_NO_STEP2_any | c1,c2,c3   |             0.194 |            0.178 |   94422.941 |           89167.822 |                0.944 |
| A2011      | ALL      | all    | D2_hgb_NO_STEP2_any | c4         |             0.194 |            0.179 |   94422.941 |           89130.599 |                0.944 |
| A2011      | ALL      | all    | D2_hgb_FULL_own     | c1,c2,c3   |             0.024 |            0.007 |   24122.506 |           19288.254 |                0.800 |
| A2011      | ALL      | all    | D2_hgb_FULL_own     | c4         |             0.024 |            0.006 |   24122.506 |           19039.330 |                0.789 |
| A2011      | ALL      | all    | D2_hgb_FULL_any     | c1,c2,c3   |             0.187 |            0.161 |   94189.558 |           87599.727 |                0.930 |
| A2011      | ALL      | all    | D2_hgb_FULL_any     | c4         |             0.187 |            0.158 |   94189.558 |           86946.463 |                0.923 |
| A2011      | NW       | all    | D3_hgb_NO_STEP2_own | c1,c2,c3   |             0.131 |            0.027 |   12141.397 |            7279.361 |                0.600 |
| A2011      | NW       | all    | D3_hgb_NO_STEP2_own | c4         |             0.131 |            0.027 |   12141.397 |            7279.134 |                0.600 |
| A2011      | NW       | all    | D3_hgb_NO_STEP2_any | c1,c2,c3   |             0.234 |            0.112 |   21544.327 |           15903.551 |                0.738 |
| A2011      | NW       | all    | D3_hgb_NO_STEP2_any | c4         |             0.234 |            0.113 |   21544.327 |           15895.286 |                0.738 |
| A2011      | NW       | all    | D2_hgb_NO_STEP2_own | c1,c2,c3   |             0.112 |            0.038 |   11724.463 |            7541.073 |                0.643 |
| A2011      | NW       | all    | D2_hgb_NO_STEP2_own | c4         |             0.112 |            0.040 |   11724.463 |            7550.446 |                0.644 |
| A2011      | NW       | all    | D2_hgb_NO_STEP2_any | c1,c2,c3   |             0.239 |            0.160 |   22425.125 |           17409.640 |                0.776 |
| A2011      | NW       | all    | D2_hgb_NO_STEP2_any | c4         |             0.239 |            0.162 |   22425.125 |           17416.742 |                0.777 |
| A2011      | NW       | all    | D2_hgb_FULL_own     | c1,c2,c3   |             0.111 |            0.033 |   11619.814 |            7328.103 |                0.631 |
| A2011      | NW       | all    | D2_hgb_FULL_own     | c4         |             0.111 |            0.033 |   11619.814 |            7292.139 |                0.628 |
| A2011      | NW       | all    | D2_hgb_FULL_any     | c1,c2,c3   |             0.250 |            0.168 |   22891.151 |           17767.871 |                0.776 |
| A2011      | NW       | all    | D2_hgb_FULL_any     | c4         |             0.250 |            0.170 |   22891.151 |           17750.455 |                0.775 |
| A2011      | Olympic  | all    | D3_hgb_NO_STEP2_own | c1,c2,c3   |             0.002 |            0.002 |   16074.110 |           15869.317 |                0.987 |
| A2011      | Olympic  | all    | D3_hgb_NO_STEP2_own | c4         |             0.002 |            0.002 |   16074.110 |           15869.004 |                0.987 |
| A2011      | Olympic  | all    | D3_hgb_NO_STEP2_any | c1,c2,c3   |             0.119 |            0.116 |   57501.715 |           57094.929 |                0.993 |
| A2011      | Olympic  | all    | D3_hgb_NO_STEP2_any | c4         |             0.119 |            0.117 |   57501.715 |           57098.348 |                0.993 |
| A2011      | Olympic  | all    | D2_hgb_NO_STEP2_own | c1,c2,c3   |             0.000 |            0.000 |   12123.394 |           12053.996 |                0.994 |
| A2011      | Olympic  | all    | D2_hgb_NO_STEP2_own | c4         |             0.000 |            0.000 |   12123.394 |           12051.737 |                0.994 |
| A2011      | Olympic  | all    | D2_hgb_NO_STEP2_any | c1,c2,c3   |             0.193 |            0.193 |   67502.374 |           67326.159 |                0.997 |
| A2011      | Olympic  | all    | D2_hgb_NO_STEP2_any | c4         |             0.193 |            0.193 |   67502.374 |           67287.597 |                0.997 |
| A2011      | Olympic  | all    | D2_hgb_FULL_own     | c1,c2,c3   |             0.000 |            0.000 |   11786.682 |           11369.365 |                0.965 |
| A2011      | Olympic  | all    | D2_hgb_FULL_own     | c4         |             0.000 |            0.000 |   11786.682 |           11204.448 |                0.951 |
| A2011      | Olympic  | all    | D2_hgb_FULL_any     | c1,c2,c3   |             0.177 |            0.168 |   66517.978 |           65224.915 |                0.981 |
| A2011      | Olympic  | all    | D2_hgb_FULL_any     | c4         |             0.177 |            0.163 |   66517.978 |           64643.818 |                0.972 |

D3 expected count:

| eval_set   | region   | model                | truth             |   n_famdays |      sumP |   truth_count |   ratio |   spearman |
|:-----------|:---------|:---------------------|:------------------|------------:|----------:|--------------:|--------:|-----------:|
| V2010_test | ALL      | D3_hgb_NO_STEP2_own  | truth_own_scored  |        4572 |  7882.586 |         30921 |   0.255 |      0.524 |
| V2010_test | ALL      | D3_hgb_NO_STEP2_nbr3 | truth_nbr3_scored |        4572 | 16448.515 |        147454 |   0.112 |      0.533 |
| V2010_test | NW       | D3_hgb_NO_STEP2_own  | truth_own_scored  |        2640 |  3723.376 |         17082 |   0.218 |      0.518 |
| V2010_test | NW       | D3_hgb_NO_STEP2_nbr3 | truth_nbr3_scored |        2640 |  3802.125 |         37188 |   0.102 |      0.443 |
| V2010_test | Olympic  | D3_hgb_NO_STEP2_own  | truth_own_scored  |         784 |  4052.668 |          7056 |   0.574 |      0.744 |
| V2010_test | Olympic  | D3_hgb_NO_STEP2_nbr3 | truth_nbr3_scored |         784 | 12531.863 |         85199 |   0.147 |      0.718 |
| V2010_test | Arm B    | D3_hgb_NO_STEP2_own  | truth_own_scored  |        1148 |   106.542 |          6783 |   0.016 |      0.197 |
| V2010_test | Arm B    | D3_hgb_NO_STEP2_nbr3 | truth_nbr3_scored |        1148 |   114.527 |         25067 |   0.005 |      0.082 |
| V2008_test | ALL      | D3_hgb_NO_STEP2_own  | truth_own_scored  |        4934 |  3884.046 |         23962 |   0.162 |      0.455 |
| V2008_test | ALL      | D3_hgb_NO_STEP2_nbr3 | truth_nbr3_scored |        4934 |  4030.206 |         44153 |   0.091 |      0.392 |
| V2008_test | NW       | D3_hgb_NO_STEP2_own  | truth_own_scored  |        3036 |  3739.827 |         16129 |   0.232 |      0.536 |
| V2008_test | NW       | D3_hgb_NO_STEP2_nbr3 | truth_nbr3_scored |        3036 |  3882.530 |         25449 |   0.153 |      0.494 |
| V2008_test | SW       | D3_hgb_NO_STEP2_own  | truth_own_scored  |        1898 |   144.219 |          7833 |   0.018 |      0.259 |
| V2008_test | SW       | D3_hgb_NO_STEP2_nbr3 | truth_nbr3_scored |        1898 |   147.675 |         18704 |   0.008 |      0.216 |
| A2011      | ALL      | D3_hgb_NO_STEP2_own  | truth_own_scored  |        3828 | 12141.397 |         33306 |   0.365 |      0.475 |
| A2011      | ALL      | D3_hgb_NO_STEP2_nbr3 | truth_nbr3_scored |        3828 | 15024.034 |         50539 |   0.297 |      0.450 |
| A2011      | NW       | D3_hgb_NO_STEP2_own  | truth_own_scored  |        3828 | 12141.397 |         33306 |   0.365 |      0.475 |
| A2011      | NW       | D3_hgb_NO_STEP2_nbr3 | truth_nbr3_scored |        3828 | 15024.034 |         50539 |   0.297 |      0.450 |

### Expected count with the scrambled control subtracted (ΣP data − mean ΣP of c1–c4, per family-day)

`rho_*` = Spearman vs the true own-family array count.

| eval_set   | region   | model               |    n |      sumP |      ctl |   excess |   truth |   rho_sumP |   rho_excess |
|:-----------|:---------|:--------------------|-----:|----------:|---------:|---------:|--------:|-----------:|-------------:|
| V2010_test | ALL      | D2_hgb_NO_STEP2_own | 4572 |  7805.469 | 5730.772 | 2074.697 |   30921 |      0.524 |        0.272 |
| V2010_test | ALL      | D3_hgb_NO_STEP2_own | 4572 |  7882.586 | 5530.829 | 2351.757 |   30921 |      0.524 |        0.288 |
| V2010_test | ALL      | P1                  | 4572 |  6031.465 | 3758.868 | 2272.597 |   30921 |      0.491 |        0.239 |
| V2010_test | NW       | D2_hgb_NO_STEP2_own | 2640 |  3663.918 | 1742.142 | 1921.776 |   17082 |      0.518 |        0.359 |
| V2010_test | NW       | D3_hgb_NO_STEP2_own | 2640 |  3723.376 | 1581.809 | 2141.567 |   17082 |      0.518 |        0.370 |
| V2010_test | NW       | P1                  | 2640 |  2471.141 |  337.308 | 2133.832 |   17082 |      0.514 |        0.337 |
| V2010_test | Olympic  | D2_hgb_NO_STEP2_own |  784 |  4052.580 | 3922.766 |  129.814 |    7056 |      0.746 |        0.201 |
| V2010_test | Olympic  | D3_hgb_NO_STEP2_own |  784 |  4052.668 | 3876.339 |  176.329 |    7056 |      0.744 |        0.238 |
| V2010_test | Olympic  | P1                  |  784 |  3481.801 | 3381.939 |   99.863 |    7056 |      0.373 |        0.214 |
| V2010_test | Arm B    | D2_hgb_NO_STEP2_own | 1148 |    88.971 |   65.865 |   23.106 |    6783 |      0.196 |        0.065 |
| V2010_test | Arm B    | D3_hgb_NO_STEP2_own | 1148 |   106.542 |   72.681 |   33.861 |    6783 |      0.197 |        0.069 |
| V2010_test | Arm B    | P1                  | 1148 |    78.523 |   39.621 |   38.902 |    6783 |      0.180 |        0.043 |
| V2008_test | ALL      | D2_hgb_NO_STEP2_own | 4934 |  3906.877 | 1915.535 | 1991.342 |   23962 |      0.455 |        0.235 |
| V2008_test | ALL      | D3_hgb_NO_STEP2_own | 4934 |  3884.046 | 1691.253 | 2192.793 |   23962 |      0.455 |        0.247 |
| V2008_test | ALL      | P1                  | 4934 |  2494.646 |  331.980 | 2162.666 |   23962 |      0.452 |        0.229 |
| V2008_test | NW       | D2_hgb_NO_STEP2_own | 3036 |  3764.898 | 1897.073 | 1867.825 |   16129 |      0.537 |        0.308 |
| V2008_test | NW       | D3_hgb_NO_STEP2_own | 3036 |  3739.827 | 1674.689 | 2065.138 |   16129 |      0.536 |        0.325 |
| V2008_test | NW       | P1                  | 3036 |  2333.408 |  313.192 | 2020.216 |   16129 |      0.519 |        0.292 |
| V2008_test | SW       | D2_hgb_NO_STEP2_own | 1898 |   141.979 |   18.463 |  123.517 |    7833 |      0.259 |        0.110 |
| V2008_test | SW       | D3_hgb_NO_STEP2_own | 1898 |   144.219 |   16.564 |  127.655 |    7833 |      0.259 |        0.111 |
| V2008_test | SW       | P1                  | 1898 |   161.238 |   18.788 |  142.450 |    7833 |      0.278 |        0.125 |
| A2011      | ALL      | D2_hgb_NO_STEP2_own | 3828 | 11724.463 | 7543.417 | 4181.046 |   33306 |      0.475 |        0.342 |
| A2011      | ALL      | D3_hgb_NO_STEP2_own | 3828 | 12141.397 | 7279.305 | 4862.092 |   33306 |      0.475 |        0.339 |
| A2011      | ALL      | P1                  | 3828 |  8464.658 | 3961.714 | 4502.943 |   33306 |      0.476 |        0.308 |
| A2011      | NW       | D2_hgb_NO_STEP2_own | 3828 | 11724.463 | 7543.417 | 4181.046 |   33306 |      0.475 |        0.342 |
| A2011      | NW       | D3_hgb_NO_STEP2_own | 3828 | 12141.397 | 7279.305 | 4862.092 |   33306 |      0.475 |        0.339 |
| A2011      | NW       | P1                  | 3828 |  8464.658 | 3961.714 | 4502.943 |   33306 |      0.476 |        0.308 |

## Caveats

- The array catalogue is incomplete: P is the probability that a candidate MATCHES THE CATALOGUE, a lower bound on P(real). A real LFE the array run missed counts as a negative.
- 'any family' matches include chance coincidences on busy tremor days (the scrambled control matches 'any' at ~5–8 %), so P(any) mixes real-LFE evidence with 'busy hour' evidence; P(own) is the cleaner label for firing rates.
- The model is fitted in the array era (2008–2011). Before using it after 2012, re-check: (1) station roster / epochs (the site-count and n_z features, dom_top and wsum depend on which stations are live; new, dead or replaced stations shift them); (2) instrument / response changes at PB boreholes and UW stations (PNET inputs, HF ratio); (3) the fine-tuned PNET on stations and years it never saw (use PNET_HO-type checks with a post-2012 known answer, e.g. PNSN tremor-day coincidence or a later temporary array); (4) the region mix and base rate (the isotonic map assumes the array-era matched fraction; a quieter or busier epoch shifts calibration in the large); (5) the density features (dens_own / dens_grp) whose meaning depends on the candidate rate per day; (6) templates — post-2012 data never contributed to them, like the test half here (D2 is the right analogue, D1 is not).
