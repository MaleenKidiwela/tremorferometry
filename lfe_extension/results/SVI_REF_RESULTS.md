# SVI_REF results (auto-written)

## Plain summary

The southern-Vancouver-Island refinement queue (results/SVI_REFINE_PLAN.md as amended by Merlin) ran on the data already on disk. Nothing was downloaded, nothing was deleted, and no v1.1 file was overwritten (re-hash of 5173 v1.1 files: 0 changed; `results/refine/r0/v11_verify.json`).

**Completeness (R1).** 67,667 expected station-days 2014-2026. v1.1 had silently lost 136 station-days at the cache stage; they were rebuilt into `cache20hz_r1/` (136 built, 0 errors). After the rebuild 0 station-months lose > 3 % at the cache stage, 0 at the raw stage. 10 station-months have > 3 % of windows masked by the energy floor (dead or gappy channels; handed to R2). (`results/refine/r1/r1_summary.json`)

**Station health (R2), LRIV positive control.** Track A-2's LRIV number is reproduced exactly from A-2's own inputs (cc median -0.216 against the stored -0.216). The same A-2 procedure on the R2 event set (v1.1 catalogue repeats inside tremor windows, 2021 doy 276-302) gives LRIV 0.250 (7 families) and B926 0.469. R2's own test (cc at the detector's fixed placement, against the station read on other days) puts LRIV after 2019-05-30 at a median 0.205 against a control 95th percentile of 0.092; it flags LRIV as at control level in 1 of 15 tested months. So LRIV is much weaker than B926 or PGC, but whether it is at control level depends on the event set. LRIV stays dropped by the v1.1 policy. (`results/refine/r2/r2_positive_control.json`, `results/refine/r2/pc_a2.json`)

**Cleaning A/B (R3).** Adopt the southern cleaning: **False** (recall at matched control rate False, known-family veto True, control falls more than data on the flagged months True: data x0.844, control x0.796). R4 was skipped; the production pass kept the v1.1 data path and stores the cleaning flags as columns. (`results/refine/r3/r3_decision.json`)

**Confidence (R5), the one-time 2013 check** (opened 2026-10-02T14:41:00.105127+00:00; nothing fitted on it). Headline model WAVE: S15_same: AUC 0.916 vs net_mad 0.771, slope 0.99, matched/sumP 0.83, control P>0.5 2.42 %; S15_any: AUC 0.871 vs net_mad 0.724, slope 0.93, matched/sumP 0.78, control P>0.5 40.77 %; S14_same: AUC 0.916 vs net_mad 0.769, slope 0.99, matched/sumP 0.84, control P>0.5 2.63 %; S14_any: AUC 0.871 vs net_mad 0.725, slope 0.93, matched/sumP 0.78, control P>0.5 40.83 %; S13_same: AUC 0.920 vs net_mad 0.775, slope 0.97, matched/sumP 0.80, control P>0.5 1.76 %; S13_any: AUC 0.869 vs net_mad 0.710, slope 0.89, matched/sumP 0.72, control P>0.5 33.09 %. (`results/refine/conf/sealed_bars.json`, `results/refine/conf/eval_sealed.csv`)

**PNET gate.** dAUC(same) 0.006 [0.003, 0.011], slope 0.97 -> not adopted: no production PNET (`results/refine/conf/pnet_gate.json`)

**Inter-ETS validation.** withheld S-energy True, control sum of P flat False, station-set invariance False -> inter-ETS rates reported as **control subtracted bracket**. (`results/refine/val/val_summary.json`)

**Rates (R7).** Bostock-equivalent repeats per day, network total, median month 9.3 (range -0.1-507.4); scrambled-control sum of P subtracted and shown. Step checks: 912 of 5832 flagged (|z| > 2 vs the same calendar windows of other years), 912 of them next to tremor days (not separable from activity). (`results/refine/rates/monthly_total.csv`, `results/refine/rates/step_checks.csv`)

**Repeat windows (R8).** 39,935,529 windows, 79.2 GB from 1578 GB raw (194 months). Verify: PASS: 20 random windows re-read from HDF5 vs the raw miniSEED (obspy slice): 20 sample-identical, 87300 samples compared. Stack check: 8/8 families pass. (`results/cut/svi_ext/`)

**Delete candidates (R9).** candidates 9025 files, 0.193 TB ({'cache20hz': 93.2, 'cache20hz_r1': 0.2, 'raw': 99.2} GB). Nothing was deleted; the list is input to a joint SVI + Washington decision. (`results/refine/r9/`)

Caveat: P is the probability of matching Bostock's ETS-only catalogue. Between ETS windows it is a lower bound on 'real', not a measure of it; the inter-ETS checks above say how far the rates can be read there. 'Scrambled control' = the same detector run on the same day with each station site's data shifted to a different time of that day.

## Stages and files

* R0 watcher: {'pid': 629057, 'still_running': False, 'ps': '', 'stop_file_present': False, 'watcher_log_last': ['[watch] /home/jovyan/tremorferometry/lfe_extension/results/trackC/STOP present — stopping after the finished month'], 'checked_utc': '2026-10-02T01:18:53.504137+00:00', 'note': 'STOP touched 2026-10-01T22:28:14Z; nothing was killed'}  (`results/refine/r0/watcher.json`); v1.1 hash snapshot `results/refine/r0/v11_hashes.csv`
* R1: failing station-months (cache) []; (mf, re-run by the SVI_REF pass) [['NLLB', '2018-06', 0.3], ['PGC', '2018-06', 0.3], ['SNB', '2018-06', 0.3], ['VGZ', '2018-06', 0.3], ['B001', '2018-06', 0.3], ['B003', '2018-06', 0.276], ['B004', '2018-06', 0.3], ['B005', '2018-06', 0.3], ['B006', '2018-06', 0.125], ['B007', '2018-06', 0.3], ['B009', '2018-06', 0.3], ['B010', '2018-06', 0.3], ['B011', '2018-06', 0.583], ['B926', '2018-06', 0.3], ['LRIV', '2017-10', 0.032], ['LRIV', '2018-06', 0.3]]; energy-floor > 3 % [['NLLB', '2014-12', 0.03], ['NLLB', '2017-02', 0.041], ['SNB', '2014-12', 0.308], ['B005', '2024-05', 0.031], ['LRIV', '2018-03', 1.0], ['LRIV', '2018-04', 0.939], ['LRIV', '2019-02', 0.983], ['LRIV', '2019-03', 0.996], ['LRIV', '2019-04', 0.987], ['LRIV', '2019-05', 0.959]]  (`results/refine/r1/`)
* R8 pilot: verify True, stack check True, extrapolation {'pilot_month': '2024-08', 'pilot_cut_GB': 1.47, 'pilot_raw_GB': 10.601, 'cut_over_raw': 0.13866616356947456, 'n_windows': 852237, 'n_det': 20358, 'minutes': 2.3, 'raw_total_TB': 1.579045300736, 'extrapolated_TB_upper': 0.21896015395546833, 'note': 'upper bound: the pilot month is the busiest v1.1 month, and the full cut also adds the v2 / truth candidates, which mostly coincide with v1.1 ones', 'pass_lt_0p5TB': True}  (`results/refine/cut_pilot/`)
* R2: station-month status {'inherit_ok': 1437, 'ok': 732, 'inherit_control_level': 85, 'nodata': 29, 'control_level': 25, 'dead': 9}; exclusions applied to the production pass: 115  (`results/refine/r2/r2_health_monthly.csv`, `r2_effective_sets.csv`, `r2_swap_checks.csv`, `r2_bad_for_dvv.csv`)

  CNSN swap checks (cc at placement, before -> after; control p95; best channel transform after):

  * PGC 2015-01-12: 0.693 -> 0.729 (control 0.095 / 0.094), lag 0.00 -> 0.00 s, identity transform best after: True, status after ok;ok;ok
  * NLLB 2017-02-23: 0.391 -> 0.471 (control 0.122 / 0.117), lag 0.00 -> 0.00 s, identity transform best after: True, status after ok;ok;ok
  * VGZ 2017-03-03: 0.516 -> 0.390 (control 0.116 / 0.133), lag 0.00 -> 0.00 s, identity transform best after: True, status after ok;ok;ok
  * SNB 2017-06-08: 0.502 -> 0.450 (control 0.136 / 0.123), lag 0.00 -> 0.00 s, identity transform best after: True, status after ok;ok;ok
  * PGC 2017-08-10: 0.461 -> 0.746 (control 0.101 / 0.091), lag 0.00 -> 0.00 s, identity transform best after: True, status after ok;ok;ok
  * LRIV 2019-05-30: 0.191 -> 0.205 (control 0.089 / 0.091), lag 0.00 -> 0.05 s, identity transform best after: True, status after ok;ok;ok

* R3 recall table: [{'year': 2011, 'T_v11': 8.25, 'ctl_rate_target': 1.6826, 'T_clean': 8.05, 'ctl_rate_clean': 1.6822, 'n_events': 18663, 'recall_v11': 0.1995, 'recall_clean': 0.1806, 'matched_v11': 3724, 'matched_clean': 3370, 'kept_of_v11': 0.8534, 'n_data_v11': 9905, 'n_data_clean': 9661, 'n_fam_rediscovered_v11': 87, 'frac_fam_kept': 0.9885}, {'year': 2011, 'T_v11': 8.75, 'ctl_rate_target': 0.9273, 'T_clean': 8.52, 'ctl_rate_clean': 0.9265, 'n_events': 18663, 'recall_v11': 0.1611, 'recall_clean': 0.1487, 'matched_v11': 3006, 'matched_clean': 2776, 'kept_of_v11': 0.856, 'n_data_v11': 6302, 'n_data_clean': 6238, 'n_fam_rediscovered_v11': 85, 'frac_fam_kept': 0.9647}, {'year': 2012, 'T_v11': 8.25, 'ctl_rate_target': 1.0548, 'T_clean': 8.1, 'ctl_rate_clean': 1.0412, 'n_events': 13381, 'recall_v11': 0.1842, 'recall_clean': 0.1707, 'matched_v11': 2465, 'matched_clean': 2284, 'kept_of_v11': 0.8771, 'n_data_v11': 6390, 'n_data_clean': 6115, 'n_fam_rediscovered_v11': 88, 'frac_fam_kept': 0.9773}, {'year': 2012, 'T_v11': 8.75, 'ctl_rate_target': 0.6113, 'T_clean': 8.57, 'ctl_rate_clean': 0.611, 'n_events': 13381, 'recall_v11': 0.1463, 'recall_clean': 0.1379, 'matched_v11': 1958, 'matched_clean': 1845, 'kept_of_v11': 0.8723, 'n_data_v11': 4031, 'n_data_clean': 3882, 'n_fam_rediscovered_v11': 83, 'frac_fam_kept': 0.988}, {'year': 2010, 'T_v11': 8.25, 'ctl_rate_target': 1.5402, 'T_clean': 8.16, 'ctl_rate_clean': 1.538, 'n_events': 18114, 'recall_v11': 0.2053, 'recall_clean': 0.1973, 'matched_v11': 3719, 'matched_clean': 3573, 'kept_of_v11': 0.9285, 'n_data_v11': 9712, 'n_data_clean': 9531, 'n_fam_rediscovered_v11': 97, 'frac_fam_kept': 0.9897}, {'year': 2010, 'T_v11': 8.75, 'ctl_rate_target': 0.8107, 'T_clean': 8.66, 'ctl_rate_clean': 0.8099, 'n_events': 18114, 'recall_v11': 0.1663, 'recall_clean': 0.1583, 'matched_v11': 3012, 'matched_clean': 2867, 'kept_of_v11': 0.9226, 'n_data_v11': 6223, 'n_data_clean': 6078, 'n_fam_rediscovered_v11': 91, 'frac_fam_kept': 1.0}]
* R3 months: [{'month': '2017-06', 'T_ext': 8.0, 'data_v11': 4486, 'data_clean': 3690, 'c1_v11': 2962, 'c1_clean': 2414}, {'month': '2017-09', 'T_ext': 8.0, 'data_v11': 6349, 'data_clean': 5313, 'c1_v11': 2992, 'c1_clean': 2293}, {'month': '2018-10', 'T_ext': 8.0, 'data_v11': 8330, 'data_clean': 7195, 'c1_v11': 3827, 'c1_clean': 3039}, {'month': '2022-01', 'T_ext': 8.0, 'data_v11': 6560, 'data_clean': 6016, 'c1_v11': 3424, 'c1_clean': 2969}, {'month': '2023-10', 'T_ext': 8.0, 'data_v11': 7816, 'data_clean': 5849, 'c1_v11': 4560, 'c1_clean': 3242}, {'month': '2023-11', 'T_ext': 8.0, 'data_v11': 6452, 'data_clean': 5682, 'c1_v11': 3832, 'c1_clean': 3231}]

* Models frozen 2026-10-02T14:40:50.329051+00:00 (sha256 686e4de45c423e11), fit rows 73078 on 49 days; features `results/refine/conf/MODELS_FROZEN.json`; never features: ['conf_z', 'conf_ze', 'tremor', 'quiet', 'date', 'season', 'ez_*', 'cc_* (per station)']

### Development (out-of-fold, fit days) -- `results/refine/conf/eval_dev.csv`

| eval    | sset   | label   | model      |     n |   n_pos |   auc |   auc_net_mad |   slope |   matched_over_sumP |   ctl_frac_P_gt05 |
|:--------|:-------|:--------|:-----------|------:|--------:|------:|--------------:|--------:|--------------------:|------------------:|
| dev_oof | S15    | same    | WAVE       | 21130 |    6931 | 0.872 |         0.738 |   1.025 |               1.01  |             0.033 |
| dev_oof | S15    | same    | FULL       | 21130 |    6931 | 0.894 |         0.738 |   1.033 |               1.012 |             0.05  |
| dev_oof | S15    | same    | NO_FAM     | 21130 |    6931 | 0.888 |         0.738 |   1.019 |               1.002 |             0.038 |
| dev_oof | S15    | same    | NO_DENSITY | 21130 |    6931 | 0.878 |         0.738 |   1.039 |               1.009 |             0.05  |
| dev_oof | S15    | same    | MF_ONLY    | 21130 |    6931 | 0.819 |         0.738 |   1.01  |               1.016 |             0.125 |
| dev_oof | S15    | same    | NO_PNET    | 21130 |    6931 | 0.875 |         0.738 |   1.026 |               1.01  |             0.016 |
| dev_oof | S15    | same    | WITH_PNET  | 21130 |    6931 | 0.88  |         0.738 |   1.024 |               1.003 |             0.007 |
| dev_oof | S15    | any     | WAVE       | 21130 |   11248 | 0.844 |         0.706 |   1.031 |               1.005 |             0.566 |
| dev_oof | S15    | any     | FULL       | 21130 |   11248 | 0.874 |         0.706 |   1.04  |               1.009 |             0.408 |
| dev_oof | S15    | any     | NO_FAM     | 21130 |   11248 | 0.871 |         0.706 |   1.029 |               1.008 |             0.393 |
| dev_oof | S15    | any     | NO_DENSITY | 21130 |   11248 | 0.85  |         0.706 |   1.04  |               1.006 |             0.558 |
| dev_oof | S15    | any     | MF_ONLY    | 21130 |   11248 | 0.814 |         0.706 |   1.027 |               1.016 |             0.586 |
| dev_oof | S15    | any     | NO_PNET    | 21130 |   11248 | 0.853 |         0.706 |   1.033 |               1.012 |             0.303 |
| dev_oof | S15    | any     | WITH_PNET  | 21130 |   11248 | 0.876 |         0.706 |   1.029 |               1.009 |             0.215 |
| dev_oof | S14    | same    | WAVE       | 21233 |    6973 | 0.871 |         0.74  |   1.025 |               1.011 |             0.032 |
| dev_oof | S14    | same    | FULL       | 21233 |    6973 | 0.893 |         0.74  |   1.029 |               1.008 |             0.051 |
| dev_oof | S14    | same    | NO_FAM     | 21233 |    6973 | 0.888 |         0.74  |   1.019 |               1.002 |             0.038 |
| dev_oof | S14    | same    | NO_DENSITY | 21233 |    6973 | 0.878 |         0.74  |   1.035 |               1.004 |             0.05  |
| dev_oof | S14    | same    | MF_ONLY    | 21233 |    6973 | 0.817 |         0.74  |   1.015 |               1.024 |             0.118 |
| dev_oof | S14    | same    | NO_PNET    | 21233 |    6973 | 0.873 |         0.74  |   1.026 |               1.012 |             0.015 |
| dev_oof | S14    | same    | WITH_PNET  | 21233 |    6973 | 0.879 |         0.74  |   1.023 |               1.005 |             0.007 |
| dev_oof | S14    | any     | WAVE       | 21233 |   11338 | 0.844 |         0.705 |   1.032 |               1.006 |             0.567 |
| dev_oof | S14    | any     | FULL       | 21233 |   11338 | 0.873 |         0.705 |   1.041 |               1.007 |             0.41  |
| dev_oof | S14    | any     | NO_FAM     | 21233 |   11338 | 0.871 |         0.705 |   1.031 |               1.01  |             0.394 |
| dev_oof | S14    | any     | NO_DENSITY | 21233 |   11338 | 0.849 |         0.705 |   1.04  |               1.006 |             0.558 |
| dev_oof | S14    | any     | MF_ONLY    | 21233 |   11338 | 0.814 |         0.705 |   1.027 |               1.014 |             0.594 |
| dev_oof | S14    | any     | NO_PNET    | 21233 |   11338 | 0.853 |         0.705 |   1.034 |               1.012 |             0.305 |
| dev_oof | S14    | any     | WITH_PNET  | 21233 |   11338 | 0.875 |         0.705 |   1.029 |               1.01  |             0.218 |
| dev_oof | S13    | same    | WAVE       | 30715 |    7164 | 0.871 |         0.734 |   1.027 |               0.983 |             0.022 |
| dev_oof | S13    | same    | FULL       | 30715 |    7164 | 0.897 |         0.734 |   1.051 |               0.957 |             0.039 |
| dev_oof | S13    | same    | NO_FAM     | 30715 |    7164 | 0.891 |         0.734 |   1.034 |               0.986 |             0.027 |
| dev_oof | S13    | same    | NO_DENSITY | 30715 |    7164 | 0.879 |         0.734 |   1.049 |               0.948 |             0.036 |
| dev_oof | S13    | same    | MF_ONLY    | 30715 |    7164 | 0.816 |         0.734 |   0.996 |               0.961 |             0.08  |
| dev_oof | S13    | same    | NO_PNET    | 30715 |    7164 | 0.875 |         0.734 |   1.029 |               0.983 |             0.011 |
| dev_oof | S13    | same    | WITH_PNET  | 30715 |    7164 | 0.881 |         0.734 |   1.02  |               0.979 |             0.005 |
| dev_oof | S13    | any     | WAVE       | 30715 |   14011 | 0.834 |         0.689 |   1.031 |               0.978 |             0.485 |
| dev_oof | S13    | any     | FULL       | 30715 |   14011 | 0.866 |         0.689 |   1.037 |               0.968 |             0.38  |
| dev_oof | S13    | any     | NO_FAM     | 30715 |   14011 | 0.863 |         0.689 |   1.023 |               0.974 |             0.361 |
| dev_oof | S13    | any     | NO_DENSITY | 30715 |   14011 | 0.841 |         0.689 |   1.047 |               0.967 |             0.49  |
| dev_oof | S13    | any     | MF_ONLY    | 30715 |   14011 | 0.804 |         0.689 |   1.008 |               0.972 |             0.515 |
| dev_oof | S13    | any     | NO_PNET    | 30715 |   14011 | 0.843 |         0.689 |   1.028 |               0.988 |             0.261 |
| dev_oof | S13    | any     | WITH_PNET  | 30715 |   14011 | 0.866 |         0.689 |   1.017 |               0.986 |             0.196 |
| s9_2010 | S9     | same    | WAVE       | 12454 |    4105 | 0.867 |         0.747 |   1.029 |               0.798 |             0.051 |
| s9_2010 | S9     | same    | FULL       | 12454 |    4105 | 0.897 |         0.747 |   1.013 |               0.798 |             0.054 |
| s9_2010 | S9     | same    | NO_FAM     | 12454 |    4105 | 0.885 |         0.747 |   1.005 |               0.805 |             0.04  |
| s9_2010 | S9     | same    | NO_DENSITY | 12454 |    4105 | 0.88  |         0.747 |   1.036 |               0.785 |             0.071 |
| s9_2010 | S9     | same    | MF_ONLY    | 12454 |    4105 | 0.754 |         0.747 |   1.028 |               1.005 |             0.073 |
| s9_2010 | S9     | any     | WAVE       | 12454 |    6277 | 0.788 |         0.708 |   1.028 |               0.788 |             0.715 |
| s9_2010 | S9     | any     | FULL       | 12454 |    6277 | 0.829 |         0.708 |   1.05  |               0.82  |             0.464 |
| s9_2010 | S9     | any     | NO_FAM     | 12454 |    6277 | 0.822 |         0.708 |   1.027 |               0.814 |             0.46  |
| s9_2010 | S9     | any     | NO_DENSITY | 12454 |    6277 | 0.799 |         0.708 |   1.041 |               0.815 |             0.648 |
| s9_2010 | S9     | any     | MF_ONLY    | 12454 |    6277 | 0.726 |         0.708 |   0.932 |               0.827 |             0.771 |

### Sealed 2013 (one look) -- `results/refine/conf/eval_sealed.csv`

| sset   | label   | model      |     n |   n_pos |   auc |   auc_net_mad |   slope |   matched_over_sumP |   ctl_frac_P_gt05 |
|:-------|:--------|:-----------|------:|--------:|------:|--------------:|--------:|--------------------:|------------------:|
| S15    | same    | WAVE       | 13077 |    3263 | 0.916 |         0.771 |   0.989 |               0.832 |             0.024 |
| S15    | same    | FULL       | 13077 |    3263 | 0.933 |         0.771 |   0.987 |               0.874 |             0.03  |
| S15    | same    | NO_FAM     | 13077 |    3263 | 0.927 |         0.771 |   0.983 |               0.858 |             0.02  |
| S15    | same    | NO_DENSITY | 13077 |    3263 | 0.921 |         0.771 |   0.994 |               0.863 |             0.032 |
| S15    | same    | MF_ONLY    | 13077 |    3263 | 0.859 |         0.771 |   0.982 |               0.867 |             0.091 |
| S15    | same    | NO_PNET    | 13077 |    3263 | 0.914 |         0.771 |   0.983 |               0.852 |             0.022 |
| S15    | same    | WITH_PNET  | 13077 |    3263 | 0.919 |         0.771 |   0.977 |               0.861 |             0.012 |
| S15    | any     | WAVE       | 13077 |    4727 | 0.871 |         0.724 |   0.934 |               0.779 |             0.408 |
| S15    | any     | FULL       | 13077 |    4727 | 0.894 |         0.724 |   0.944 |               0.819 |             0.25  |
| S15    | any     | NO_FAM     | 13077 |    4727 | 0.89  |         0.724 |   0.939 |               0.815 |             0.248 |
| S15    | any     | NO_DENSITY | 13077 |    4727 | 0.876 |         0.724 |   0.934 |               0.793 |             0.393 |
| S15    | any     | MF_ONLY    | 13077 |    4727 | 0.837 |         0.724 |   0.897 |               0.777 |             0.504 |
| S15    | any     | NO_PNET    | 13077 |    4727 | 0.863 |         0.724 |   0.891 |               0.805 |             0.414 |
| S15    | any     | WITH_PNET  | 13077 |    4727 | 0.881 |         0.724 |   0.899 |               0.814 |             0.294 |
| S14    | same    | WAVE       | 13082 |    3266 | 0.916 |         0.769 |   0.991 |               0.837 |             0.026 |
| S14    | same    | FULL       | 13082 |    3266 | 0.932 |         0.769 |   0.986 |               0.872 |             0.031 |
| S14    | same    | NO_FAM     | 13082 |    3266 | 0.927 |         0.769 |   0.984 |               0.858 |             0.024 |
| S14    | same    | NO_DENSITY | 13082 |    3266 | 0.92  |         0.769 |   0.991 |               0.862 |             0.034 |
| S14    | same    | MF_ONLY    | 13082 |    3266 | 0.858 |         0.769 |   0.992 |               0.872 |             0.086 |
| S14    | same    | NO_PNET    | 13082 |    3266 | 0.913 |         0.769 |   0.984 |               0.856 |             0.024 |
| S14    | same    | WITH_PNET  | 13082 |    3266 | 0.919 |         0.769 |   0.978 |               0.866 |             0.014 |
| S14    | any     | WAVE       | 13082 |    4717 | 0.871 |         0.725 |   0.935 |               0.777 |             0.408 |
| S14    | any     | FULL       | 13082 |    4717 | 0.893 |         0.725 |   0.938 |               0.815 |             0.255 |
| S14    | any     | NO_FAM     | 13082 |    4717 | 0.889 |         0.725 |   0.936 |               0.812 |             0.249 |
| S14    | any     | NO_DENSITY | 13082 |    4717 | 0.875 |         0.725 |   0.934 |               0.791 |             0.396 |
| S14    | any     | MF_ONLY    | 13082 |    4717 | 0.839 |         0.725 |   0.895 |               0.77  |             0.511 |
| S14    | any     | NO_PNET    | 13082 |    4717 | 0.863 |         0.725 |   0.891 |               0.802 |             0.42  |
| S14    | any     | WITH_PNET  | 13082 |    4717 | 0.881 |         0.725 |   0.896 |               0.811 |             0.292 |
| S13    | same    | WAVE       | 20396 |    3316 | 0.92  |         0.775 |   0.966 |               0.802 |             0.018 |
| S13    | same    | FULL       | 20396 |    3316 | 0.941 |         0.775 |   0.98  |               0.81  |             0.03  |
| S13    | same    | NO_FAM     | 20396 |    3316 | 0.934 |         0.775 |   0.983 |               0.832 |             0.019 |
| S13    | same    | NO_DENSITY | 20396 |    3316 | 0.927 |         0.775 |   0.98  |               0.807 |             0.026 |
| S13    | same    | MF_ONLY    | 20396 |    3316 | 0.867 |         0.775 |   0.943 |               0.798 |             0.053 |
| S13    | same    | NO_PNET    | 20396 |    3316 | 0.917 |         0.775 |   0.961 |               0.824 |             0.015 |
| S13    | same    | WITH_PNET  | 20396 |    3316 | 0.925 |         0.775 |   0.955 |               0.828 |             0.009 |
| S13    | any     | WAVE       | 20396 |    5589 | 0.869 |         0.71  |   0.888 |               0.715 |             0.331 |
| S13    | any     | FULL       | 20396 |    5589 | 0.893 |         0.71  |   0.905 |               0.76  |             0.223 |
| S13    | any     | NO_FAM     | 20396 |    5589 | 0.886 |         0.71  |   0.898 |               0.759 |             0.21  |
| S13    | any     | NO_DENSITY | 20396 |    5589 | 0.877 |         0.71  |   0.901 |               0.73  |             0.328 |
| S13    | any     | MF_ONLY    | 20396 |    5589 | 0.838 |         0.71  |   0.817 |               0.687 |             0.443 |
| S13    | any     | NO_PNET    | 20396 |    5589 | 0.862 |         0.71  |   0.843 |               0.744 |             0.349 |
| S13    | any     | WITH_PNET  | 20396 |    5589 | 0.88  |         0.71  |   0.85  |               0.755 |             0.259 |

### Inter-ETS validation

```
{
 "energy": {
  "n_quiet_data": 155789,
  "n_quiet_ctl": 105065,
  "n_data_P_ge_05": 783,
  "diff_hiP_minus_ctl": [
   1.8135953468282544,
   1.4530588950055614,
   2.132397164053484
  ],
  "spearman_bins": 1.0,
  "pass": true
 },
 "control_sumP": {
  "by_epoch": [
   {
    "epoch": "E0",
    "data": 10688.741,
    "ctl": 4676.0294,
    "days": 1153,
    "ctl_over_data": 0.4375
   },
   {
    "epoch": "E1",
    "data": 8835.2627,
    "ctl": 2972.3574,
    "days": 670,
    "ctl_over_data": 0.3364
   },
   {
    "epoch": "E2",
    "data": 11370.4015,
    "ctl": 3980.1833,
    "days": 821,
    "ctl_over_data": 0.35
   },
   {
    "epoch": "E3",
    "data": 26883.926,
    "ctl": 13434.3514,
    "days": 1993,
    "ctl_over_data": 0.4997
   }
  ],
  "spearman_time": 0.06101818717917788,
  "pass": false
 },
 "invariance": {
  "available": true,
  "k_S15": 1.0102530869484787,
  "k_S13": 0.9834255155876075,
  "n_months": 87,
  "pooled_ratio": 1.0927024452071814,
  "median_monthly_ratio": 1.236850238612612,
  "pass": false
 },
 "all_pass": false,
 "inter_ets_mode": "control_subtracted_bracket"
}
```

* k = his / sumP per (family, set), shrunk: median by set {'S13': 4.399, 'S14': 4.647, 'S15': 4.671}; tau2 {'S13': 0.71, 'S14': 0.775, 'S15': 0.786}  (`results/refine/rates/k_factors.csv`)

### Flagged steps (network rate)

| date       | change                                       |   window_d |   ratio_after_before |   z_vs_other_years |   tremor_days_in_window |
|:-----------|:---------------------------------------------|-----------:|---------------------:|-------------------:|------------------------:|
| 2015-12-01 | R2 B001 inherit_ok->ok                       |         45 |               38.736 |              2.495 |                      23 |
| 2016-02-01 | R2 B001 ok->inherit_ok                       |         45 |                0.024 |             -2.834 |                      23 |
| 2017-03-01 | R2 B001 inherit_ok->ok                       |         21 |                0.014 |             -2.31  |                      13 |
| 2017-04-01 | R2 B001 ok->inherit_ok                       |         45 |                0.093 |             -2.952 |                      24 |
| 2018-01-01 | R2 B001 inherit_control_level->control_level |         21 |              246.311 |              4.486 |                      12 |
| 2018-01-01 | R2 B001 inherit_control_level->control_level |         45 |               32.042 |              2.076 |                      19 |
| 2018-06-01 | R2 B001 inherit_ok->ok                       |         45 |               18.426 |              2.649 |                      45 |
| 2018-08-01 | R2 B001 ok->inherit_ok                       |         45 |                0.035 |             -3.444 |                      34 |
| 2019-09-01 | R2 B001 ok->inherit_ok                       |         21 |                0.021 |             -2.061 |                      14 |
| 2019-09-01 | R2 B001 ok->inherit_ok                       |         45 |                0.047 |             -2.278 |                      19 |
| 2019-11-01 | R2 B001 inherit_ok->ok                       |         21 |               57.113 |              2.039 |                      18 |
| 2019-11-01 | R2 B001 inherit_ok->ok                       |         45 |               34.789 |              2.202 |                      26 |
| 2019-12-01 | R2 B001 ok->inherit_ok                       |         21 |                0.006 |             -2.865 |                      20 |
| 2020-06-01 | R2 B001 ok->inherit_ok                       |         21 |                0.065 |             -2.595 |                      10 |
| 2021-02-01 | R2 B001 inherit_ok->ok                       |         21 |              240.713 |              4.885 |                      18 |
| 2021-02-01 | R2 B001 inherit_ok->ok                       |         45 |              147.907 |              3.243 |                      20 |
| 2021-03-01 | R2 B001 ok->inherit_ok                       |         45 |                0.048 |             -2.69  |                      24 |
| 2023-04-01 | R2 B001 inherit_ok->ok                       |         21 |               54.722 |              2.423 |                      21 |
| 2023-04-01 | R2 B001 inherit_ok->ok                       |         45 |               31.651 |              2.382 |                      35 |
| 2023-05-01 | R2 B001 ok->inherit_ok                       |         21 |                0.018 |             -3.58  |                      20 |
| 2023-05-01 | R2 B001 ok->inherit_ok                       |         45 |                0.024 |             -2.565 |                      31 |
| 2024-08-01 | R2 B001 inherit_ok->ok                       |         21 |               16.471 |              2.076 |                      19 |
| 2024-09-01 | R2 B001 ok->inherit_ok                       |         45 |                0.052 |             -2.173 |                      40 |
| 2025-10-01 | R2 B001 control_level->ok                    |         21 |               33.265 |              2.273 |                      32 |
| 2015-10-01 | R2 B003 inherit_control_level->control_level |         21 |               41.624 |              2.555 |                       7 |
| 2015-12-01 | R2 B003 inherit_control_level->ok            |         45 |               38.736 |              2.495 |                      23 |
| 2016-02-01 | R2 B003 ok->inherit_ok                       |         45 |                0.024 |             -2.834 |                      23 |
| 2017-03-01 | R2 B003 inherit_ok->ok                       |         21 |                0.014 |             -2.31  |                      13 |
| 2017-04-01 | R2 B003 ok->inherit_ok                       |         45 |                0.093 |             -2.952 |                      24 |
| 2018-01-01 | R2 B003 inherit_ok->ok                       |         21 |              246.311 |              4.486 |                      12 |
| 2018-01-01 | R2 B003 inherit_ok->ok                       |         45 |               32.042 |              2.076 |                      19 |
| 2018-06-01 | R2 B003 inherit_ok->ok                       |         45 |               18.426 |              2.649 |                      45 |
| 2018-08-01 | R2 B003 ok->inherit_ok                       |         45 |                0.035 |             -3.444 |                      34 |
| 2019-09-01 | R2 B003 ok->inherit_ok                       |         21 |                0.021 |             -2.061 |                      14 |
| 2019-09-01 | R2 B003 ok->inherit_ok                       |         45 |                0.047 |             -2.278 |                      19 |
| 2019-11-01 | R2 B003 inherit_ok->ok                       |         21 |               57.113 |              2.039 |                      18 |
| 2019-11-01 | R2 B003 inherit_ok->ok                       |         45 |               34.789 |              2.202 |                      26 |
| 2019-12-01 | R2 B003 ok->inherit_ok                       |         21 |                0.006 |             -2.865 |                      20 |
| 2021-02-01 | R2 B003 inherit_ok->ok                       |         21 |              240.713 |              4.885 |                      18 |
| 2021-02-01 | R2 B003 inherit_ok->ok                       |         45 |              147.907 |              3.243 |                      20 |
| 2021-03-01 | R2 B003 ok->inherit_ok                       |         45 |                0.048 |             -2.69  |                      24 |
| 2023-04-01 | R2 B003 inherit_ok->ok                       |         21 |               54.722 |              2.423 |                      21 |
| 2023-04-01 | R2 B003 inherit_ok->ok                       |         45 |               31.651 |              2.382 |                      35 |
| 2023-05-01 | R2 B003 ok->inherit_ok                       |         21 |                0.018 |             -3.58  |                      20 |
| 2023-05-01 | R2 B003 ok->inherit_ok                       |         45 |                0.024 |             -2.565 |                      31 |
| 2024-08-01 | R2 B003 inherit_ok->ok                       |         21 |               16.471 |              2.076 |                      19 |
| 2024-09-01 | R2 B003 ok->inherit_ok                       |         45 |                0.052 |             -2.173 |                      40 |
| 2025-10-01 | R2 B003 inherit_ok->ok                       |         21 |               33.265 |              2.273 |                      32 |
| 2015-12-01 | R2 B004 inherit_ok->ok                       |         45 |               38.736 |              2.495 |                      23 |
| 2016-02-01 | R2 B004 ok->inherit_ok                       |         45 |                0.024 |             -2.834 |                      23 |
| 2017-03-01 | R2 B004 inherit_ok->ok                       |         21 |                0.014 |             -2.31  |                      13 |
| 2017-04-01 | R2 B004 ok->inherit_ok                       |         45 |                0.093 |             -2.952 |                      24 |
| 2018-01-01 | R2 B004 inherit_ok->ok                       |         21 |              246.311 |              4.486 |                      12 |
| 2018-01-01 | R2 B004 inherit_ok->ok                       |         45 |               32.042 |              2.076 |                      19 |
| 2018-06-01 | R2 B004 inherit_ok->ok                       |         45 |               18.426 |              2.649 |                      45 |
| 2018-08-01 | R2 B004 ok->inherit_ok                       |         45 |                0.035 |             -3.444 |                      34 |
| 2019-09-01 | R2 B004 ok->inherit_ok                       |         21 |                0.021 |             -2.061 |                      14 |
| 2019-09-01 | R2 B004 ok->inherit_ok                       |         45 |                0.047 |             -2.278 |                      19 |
| 2019-11-01 | R2 B004 inherit_ok->ok                       |         21 |               57.113 |              2.039 |                      18 |
| 2019-11-01 | R2 B004 inherit_ok->ok                       |         45 |               34.789 |              2.202 |                      26 |
| 2019-12-01 | R2 B004 ok->inherit_ok                       |         21 |                0.006 |             -2.865 |                      20 |
| 2021-02-01 | R2 B004 inherit_ok->ok                       |         21 |              240.713 |              4.885 |                      18 |
| 2021-02-01 | R2 B004 inherit_ok->ok                       |         45 |              147.907 |              3.243 |                      20 |
| 2021-03-01 | R2 B004 ok->inherit_ok                       |         45 |                0.048 |             -2.69  |                      24 |
| 2023-04-01 | R2 B004 inherit_ok->ok                       |         21 |               54.722 |              2.423 |                      21 |
| 2023-04-01 | R2 B004 inherit_ok->ok                       |         45 |               31.651 |              2.382 |                      35 |
| 2023-05-01 | R2 B004 ok->inherit_ok                       |         21 |                0.018 |             -3.58  |                      20 |
| 2023-05-01 | R2 B004 ok->inherit_ok                       |         45 |                0.024 |             -2.565 |                      31 |
| 2024-08-01 | R2 B004 inherit_ok->ok                       |         21 |               16.471 |              2.076 |                      19 |
| 2024-09-01 | R2 B004 ok->inherit_ok                       |         45 |                0.052 |             -2.173 |                      40 |
| 2025-10-01 | R2 B004 inherit_ok->ok                       |         21 |               33.265 |              2.273 |                      32 |
| 2015-12-01 | R2 B005 inherit_ok->ok                       |         45 |               38.736 |              2.495 |                      23 |
| 2016-02-01 | R2 B005 ok->inherit_ok                       |         45 |                0.024 |             -2.834 |                      23 |
| 2017-03-01 | R2 B005 inherit_ok->ok                       |         21 |                0.014 |             -2.31  |                      13 |
| 2018-01-01 | R2 B005 inherit_control_level->control_level |         21 |              246.311 |              4.486 |                      12 |
| 2018-01-01 | R2 B005 inherit_control_level->control_level |         45 |               32.042 |              2.076 |                      19 |
| 2018-06-01 | R2 B005 control_level->ok                    |         45 |               18.426 |              2.649 |                      45 |
| 2018-08-01 | R2 B005 ok->inherit_ok                       |         45 |                0.035 |             -3.444 |                      34 |
| 2019-09-01 | R2 B005 ok->inherit_ok                       |         21 |                0.021 |             -2.061 |                      14 |
| 2019-09-01 | R2 B005 ok->inherit_ok                       |         45 |                0.047 |             -2.278 |                      19 |
| 2019-11-01 | R2 B005 inherit_ok->ok                       |         21 |               57.113 |              2.039 |                      18 |
| 2019-11-01 | R2 B005 inherit_ok->ok                       |         45 |               34.789 |              2.202 |                      26 |
| 2019-12-01 | R2 B005 ok->inherit_ok                       |         21 |                0.006 |             -2.865 |                      20 |
| 2021-02-01 | R2 B005 inherit_ok->ok                       |         21 |              240.713 |              4.885 |                      18 |
| 2021-02-01 | R2 B005 inherit_ok->ok                       |         45 |              147.907 |              3.243 |                      20 |
| 2021-03-01 | R2 B005 ok->inherit_ok                       |         45 |                0.048 |             -2.69  |                      24 |
| 2023-03-01 | R2 B005 inherit_control_level->inherit_ok    |         45 |               86.965 |              2.886 |                      28 |
| 2023-04-01 | R2 B005 inherit_ok->ok                       |         21 |               54.722 |              2.423 |                      21 |
| 2023-04-01 | R2 B005 inherit_ok->ok                       |         45 |               31.651 |              2.382 |                      35 |
| 2023-05-01 | R2 B005 ok->inherit_ok                       |         21 |                0.018 |             -3.58  |                      20 |
| 2023-05-01 | R2 B005 ok->inherit_ok                       |         45 |                0.024 |             -2.565 |                      31 |
| 2024-08-01 | R2 B005 inherit_ok->ok                       |         21 |               16.471 |              2.076 |                      19 |
| 2024-09-01 | R2 B005 ok->inherit_ok                       |         45 |                0.052 |             -2.173 |                      40 |
| 2025-10-01 | R2 B005 control_level->ok                    |         21 |               33.265 |              2.273 |                      32 |
| 2015-10-01 | R2 B006 inherit_ok->ok                       |         21 |               41.624 |              2.555 |                       7 |
| 2015-12-01 | R2 B006 inherit_ok->ok                       |         45 |               38.736 |              2.495 |                      23 |
| 2016-02-01 | R2 B006 ok->inherit_ok                       |         45 |                0.024 |             -2.834 |                      23 |
| 2017-03-01 | R2 B006 inherit_ok->ok                       |         21 |                0.014 |             -2.31  |                      13 |
| 2018-06-01 | R2 B006 inherit_ok->ok                       |         45 |               18.426 |              2.649 |                      45 |
| 2018-08-01 | R2 B006 ok->inherit_ok                       |         45 |                0.035 |             -3.444 |                      34 |
| 2019-09-01 | R2 B006 ok->inherit_ok                       |         21 |                0.021 |             -2.061 |                      14 |
| 2019-09-01 | R2 B006 ok->inherit_ok                       |         45 |                0.047 |             -2.278 |                      19 |
| 2019-11-01 | R2 B006 inherit_ok->ok                       |         21 |               57.113 |              2.039 |                      18 |
| 2019-11-01 | R2 B006 inherit_ok->ok                       |         45 |               34.789 |              2.202 |                      26 |
| 2019-12-01 | R2 B006 ok->inherit_ok                       |         21 |                0.006 |             -2.865 |                      20 |
| 2021-02-01 | R2 B006 inherit_ok->ok                       |         21 |              240.713 |              4.885 |                      18 |
| 2021-02-01 | R2 B006 inherit_ok->ok                       |         45 |              147.907 |              3.243 |                      20 |
| 2021-03-01 | R2 B006 ok->inherit_ok                       |         45 |                0.048 |             -2.69  |                      24 |
| 2023-03-01 | R2 B006 nodata->inherit_ok                   |         45 |               86.965 |              2.886 |                      28 |
| 2023-04-01 | R2 B006 inherit_ok->ok                       |         21 |               54.722 |              2.423 |                      21 |
| 2023-04-01 | R2 B006 inherit_ok->ok                       |         45 |               31.651 |              2.382 |                      35 |
| 2023-05-01 | R2 B006 ok->inherit_ok                       |         21 |                0.018 |             -3.58  |                      20 |
| 2023-05-01 | R2 B006 ok->inherit_ok                       |         45 |                0.024 |             -2.565 |                      31 |
| 2024-08-01 | R2 B006 inherit_ok->ok                       |         21 |               16.471 |              2.076 |                      19 |
| 2024-09-01 | R2 B006 ok->inherit_ok                       |         45 |                0.052 |             -2.173 |                      40 |
| 2015-10-01 | R2 B007 inherit_ok->ok                       |         21 |               41.624 |              2.555 |                       7 |
| 2015-12-01 | R2 B007 inherit_ok->ok                       |         45 |               38.736 |              2.495 |                      23 |
| 2016-02-01 | R2 B007 ok->inherit_ok                       |         45 |                0.024 |             -2.834 |                      23 |
| 2017-03-01 | R2 B007 inherit_ok->ok                       |         21 |                0.014 |             -2.31  |                      13 |
| 2017-04-01 | R2 B007 ok->inherit_ok                       |         45 |                0.093 |             -2.952 |                      24 |
| 2018-06-01 | R2 B007 inherit_ok->ok                       |         45 |               18.426 |              2.649 |                      45 |
| 2018-08-01 | R2 B007 ok->inherit_ok                       |         45 |                0.035 |             -3.444 |                      34 |
| 2019-09-01 | R2 B007 ok->inherit_ok                       |         21 |                0.021 |             -2.061 |                      14 |
| 2019-09-01 | R2 B007 ok->inherit_ok                       |         45 |                0.047 |             -2.278 |                      19 |
| 2019-11-01 | R2 B007 inherit_ok->ok                       |         21 |               57.113 |              2.039 |                      18 |
| 2019-11-01 | R2 B007 inherit_ok->ok                       |         45 |               34.789 |              2.202 |                      26 |
| 2019-12-01 | R2 B007 ok->inherit_ok                       |         21 |                0.006 |             -2.865 |                      20 |
| 2021-02-01 | R2 B007 inherit_ok->ok                       |         21 |              240.713 |              4.885 |                      18 |
| 2021-02-01 | R2 B007 inherit_ok->ok                       |         45 |              147.907 |              3.243 |                      20 |
| 2021-03-01 | R2 B007 ok->inherit_ok                       |         45 |                0.048 |             -2.69  |                      24 |
| 2023-04-01 | R2 B007 inherit_ok->ok                       |         21 |               54.722 |              2.423 |                      21 |
| 2023-04-01 | R2 B007 inherit_ok->ok                       |         45 |               31.651 |              2.382 |                      35 |
| 2023-05-01 | R2 B007 ok->inherit_ok                       |         21 |                0.018 |             -3.58  |                      20 |
| 2023-05-01 | R2 B007 ok->inherit_ok                       |         45 |                0.024 |             -2.565 |                      31 |
| 2024-08-01 | R2 B007 inherit_ok->ok                       |         21 |               16.471 |              2.076 |                      19 |
| 2024-09-01 | R2 B007 ok->inherit_ok                       |         45 |                0.052 |             -2.173 |                      40 |
| 2015-10-01 | R2 B009 inherit_ok->ok                       |         21 |               41.624 |              2.555 |                       7 |
| 2015-12-01 | R2 B009 inherit_ok->ok                       |         45 |               38.736 |              2.495 |                      23 |
| 2016-02-01 | R2 B009 ok->inherit_ok                       |         45 |                0.024 |             -2.834 |                      23 |
| 2017-07-01 | R2 B009 ok->inherit_ok                       |         21 |                0.155 |             -2.3   |                       6 |
| 2018-01-01 | R2 B009 inherit_ok->ok                       |         21 |              246.311 |              4.486 |                      12 |
| 2018-01-01 | R2 B009 inherit_ok->ok                       |         45 |               32.042 |              2.076 |                      19 |
| 2018-08-01 | R2 B009 ok->inherit_ok                       |         45 |                0.035 |             -3.444 |                      34 |
| 2019-09-01 | R2 B009 ok->inherit_ok                       |         21 |                0.021 |             -2.061 |                      14 |
| 2019-09-01 | R2 B009 ok->inherit_ok                       |         45 |                0.047 |             -2.278 |                      19 |
| 2019-11-01 | R2 B009 inherit_ok->ok                       |         21 |               57.113 |              2.039 |                      18 |
| 2019-11-01 | R2 B009 inherit_ok->ok                       |         45 |               34.789 |              2.202 |                      26 |
| 2019-12-01 | R2 B009 ok->inherit_ok                       |         21 |                0.006 |             -2.865 |                      20 |
| 2020-06-01 | R2 B009 ok->inherit_ok                       |         21 |                0.065 |             -2.595 |                      10 |
| 2020-07-01 | R2 B009 inherit_ok->ok                       |         45 |                0.105 |             -2.407 |                      20 |
| 2021-02-01 | R2 B009 inherit_ok->ok                       |         21 |              240.713 |              4.885 |                      18 |
| 2021-02-01 | R2 B009 inherit_ok->ok                       |         45 |              147.907 |              3.243 |                      20 |
| 2021-03-01 | R2 B009 ok->inherit_ok                       |         45 |                0.048 |             -2.69  |                      24 |
| 2023-03-01 | R2 B009 ok->inherit_ok                       |         45 |               86.965 |              2.886 |                      28 |
| 2023-04-01 | R2 B009 inherit_ok->ok                       |         21 |               54.722 |              2.423 |                      21 |
| 2023-04-01 | R2 B009 inherit_ok->ok                       |         45 |               31.651 |              2.382 |                      35 |
| 2023-05-01 | R2 B009 ok->inherit_ok                       |         21 |                0.018 |             -3.58  |                      20 |
| 2023-05-01 | R2 B009 ok->inherit_ok                       |         45 |                0.024 |             -2.565 |                      31 |
| 2024-10-01 | R2 B009 ok->inherit_ok                       |         45 |                0.042 |             -3.814 |                      32 |
| 2025-03-01 | R2 B009 ok->inherit_ok                       |         21 |                0.017 |             -2.128 |                       4 |
| 2015-10-01 | R2 B010 inherit_ok->ok                       |         21 |               41.624 |              2.555 |                       7 |
| 2015-12-01 | R2 B010 inherit_ok->ok                       |         45 |               38.736 |              2.495 |                      23 |
| 2016-02-01 | R2 B010 ok->inherit_ok                       |         45 |                0.024 |             -2.834 |                      23 |
| 2017-03-01 | R2 B010 ok->inherit_ok                       |         21 |                0.014 |             -2.31  |                      13 |
| 2017-04-01 | R2 B010 inherit_ok->nodata                   |         45 |                0.093 |             -2.952 |                      24 |
| 2017-07-01 | R2 B010 ok->inherit_ok                       |         21 |                0.155 |             -2.3   |                       6 |
| 2018-01-01 | R2 B010 inherit_ok->ok                       |         21 |              246.311 |              4.486 |                      12 |
| 2018-01-01 | R2 B010 inherit_ok->ok                       |         45 |               32.042 |              2.076 |                      19 |
| 2018-08-01 | R2 B010 ok->inherit_ok                       |         45 |                0.035 |             -3.444 |                      34 |
| 2019-09-01 | R2 B010 ok->inherit_ok                       |         21 |                0.021 |             -2.061 |                      14 |
| 2019-09-01 | R2 B010 ok->inherit_ok                       |         45 |                0.047 |             -2.278 |                      19 |
| 2019-11-01 | R2 B010 inherit_ok->ok                       |         21 |               57.113 |              2.039 |                      18 |
| 2019-11-01 | R2 B010 inherit_ok->ok                       |         45 |               34.789 |              2.202 |                      26 |
| 2019-12-01 | R2 B010 ok->inherit_ok                       |         21 |                0.006 |             -2.865 |                      20 |
| 2020-06-01 | R2 B010 ok->inherit_ok                       |         21 |                0.065 |             -2.595 |                      10 |
| 2020-07-01 | R2 B010 inherit_ok->ok                       |         45 |                0.105 |             -2.407 |                      20 |
| 2021-02-01 | R2 B010 inherit_ok->ok                       |         21 |              240.713 |              4.885 |                      18 |
| 2021-02-01 | R2 B010 inherit_ok->ok                       |         45 |              147.907 |              3.243 |                      20 |
| 2021-03-01 | R2 B010 ok->inherit_ok                       |         45 |                0.048 |             -2.69  |                      24 |
| 2023-03-01 | R2 B010 ok->inherit_ok                       |         45 |               86.965 |              2.886 |                      28 |
| 2023-04-01 | R2 B010 inherit_ok->ok                       |         21 |               54.722 |              2.423 |                      21 |
| 2023-04-01 | R2 B010 inherit_ok->ok                       |         45 |               31.651 |              2.382 |                      35 |
| 2023-05-01 | R2 B010 ok->inherit_ok                       |         21 |                0.018 |             -3.58  |                      20 |
| 2023-05-01 | R2 B010 ok->inherit_ok                       |         45 |                0.024 |             -2.565 |                      31 |
| 2024-10-01 | R2 B010 ok->inherit_ok                       |         45 |                0.042 |             -3.814 |                      32 |
| 2025-03-01 | R2 B010 ok->inherit_ok                       |         21 |                0.017 |             -2.128 |                       4 |
| 2015-10-01 | R2 B011 inherit_ok->ok                       |         21 |               41.624 |              2.555 |                       7 |
| 2015-12-01 | R2 B011 inherit_ok->ok                       |         45 |               38.736 |              2.495 |                      23 |
| 2016-02-01 | R2 B011 ok->inherit_ok                       |         45 |                0.024 |             -2.834 |                      23 |
| 2017-07-01 | R2 B011 ok->inherit_ok                       |         21 |                0.155 |             -2.3   |                       6 |
| 2018-01-01 | R2 B011 inherit_ok->ok                       |         21 |              246.311 |              4.486 |                      12 |
| 2018-01-01 | R2 B011 inherit_ok->ok                       |         45 |               32.042 |              2.076 |                      19 |
| 2018-08-01 | R2 B011 ok->inherit_ok                       |         45 |                0.035 |             -3.444 |                      34 |
| 2019-09-01 | R2 B011 ok->inherit_ok                       |         21 |                0.021 |             -2.061 |                      14 |
| 2019-09-01 | R2 B011 ok->inherit_ok                       |         45 |                0.047 |             -2.278 |                      19 |
| 2019-11-01 | R2 B011 inherit_ok->ok                       |         21 |               57.113 |              2.039 |                      18 |
| 2019-11-01 | R2 B011 inherit_ok->ok                       |         45 |               34.789 |              2.202 |                      26 |
| 2019-12-01 | R2 B011 ok->inherit_ok                       |         21 |                0.006 |             -2.865 |                      20 |
| 2020-06-01 | R2 B011 ok->inherit_ok                       |         21 |                0.065 |             -2.595 |                      10 |
| 2020-07-01 | R2 B011 inherit_ok->ok                       |         45 |                0.105 |             -2.407 |                      20 |
| 2021-02-01 | R2 B011 inherit_ok->ok                       |         21 |              240.713 |              4.885 |                      18 |
| 2021-02-01 | R2 B011 inherit_ok->ok                       |         45 |              147.907 |              3.243 |                      20 |
| 2021-03-01 | R2 B011 ok->inherit_ok                       |         45 |                0.048 |             -2.69  |                      24 |
| 2023-03-01 | R2 B011 ok->inherit_ok                       |         45 |               86.965 |              2.886 |                      28 |
| 2023-04-01 | R2 B011 inherit_ok->ok                       |         21 |               54.722 |              2.423 |                      21 |
| 2023-04-01 | R2 B011 inherit_ok->ok                       |         45 |               31.651 |              2.382 |                      35 |
| 2023-05-01 | R2 B011 ok->inherit_ok                       |         21 |                0.018 |             -3.58  |                      20 |
| 2023-05-01 | R2 B011 ok->inherit_ok                       |         45 |                0.024 |             -2.565 |                      31 |
| 2024-10-01 | R2 B011 ok->inherit_ok                       |         45 |                0.042 |             -3.814 |                      32 |
| 2025-03-01 | R2 B011 ok->inherit_ok                       |         21 |                0.017 |             -2.128 |                       4 |
| 2015-10-01 | R2 B926 inherit_ok->ok                       |         21 |               41.624 |              2.555 |                       7 |
| 2015-12-01 | R2 B926 inherit_ok->ok                       |         45 |               38.736 |              2.495 |                      23 |
| 2016-02-01 | R2 B926 ok->inherit_ok                       |         45 |                0.024 |             -2.834 |                      23 |
| 2018-01-01 | R2 B926 inherit_ok->ok                       |         21 |              246.311 |              4.486 |                      12 |
| 2018-01-01 | R2 B926 inherit_ok->ok                       |         45 |               32.042 |              2.076 |                      19 |
| 2018-08-01 | R2 B926 ok->inherit_ok                       |         45 |                0.035 |             -3.444 |                      34 |
| 2019-09-01 | R2 B926 ok->inherit_ok                       |         21 |                0.021 |             -2.061 |                      14 |
| 2019-09-01 | R2 B926 ok->inherit_ok                       |         45 |                0.047 |             -2.278 |                      19 |
| 2019-11-01 | R2 B926 inherit_ok->ok                       |         21 |               57.113 |              2.039 |                      18 |
| 2019-11-01 | R2 B926 inherit_ok->ok                       |         45 |               34.789 |              2.202 |                      26 |
| 2019-12-01 | R2 B926 ok->inherit_ok                       |         21 |                0.006 |             -2.865 |                      20 |
| 2020-06-01 | R2 B926 ok->inherit_ok                       |         21 |                0.065 |             -2.595 |                      10 |
| 2020-07-01 | R2 B926 inherit_ok->control_level            |         45 |                0.105 |             -2.407 |                      20 |
| 2021-02-01 | R2 B926 inherit_ok->ok                       |         21 |              240.713 |              4.885 |                      18 |
| 2021-02-01 | R2 B926 inherit_ok->ok                       |         45 |              147.907 |              3.243 |                      20 |
| 2021-03-01 | R2 B926 ok->inherit_ok                       |         45 |                0.048 |             -2.69  |                      24 |
| 2023-03-01 | R2 B926 control_level->inherit_control_level |         45 |               86.965 |              2.886 |                      28 |
| 2023-04-01 | R2 B926 inherit_control_level->ok            |         21 |               54.722 |              2.423 |                      21 |
| 2023-04-01 | R2 B926 inherit_control_level->ok            |         45 |               31.651 |              2.382 |                      35 |
| 2023-05-01 | R2 B926 ok->inherit_ok                       |         21 |                0.018 |             -3.58  |                      20 |
| 2023-05-01 | R2 B926 ok->inherit_ok                       |         45 |                0.024 |             -2.565 |                      31 |
| 2024-09-01 | R2 B926 ok->inherit_ok                       |         45 |                0.052 |             -2.173 |                      40 |
| 2015-10-01 | R2 LRIV inherit_control_level->control_level |         21 |               41.624 |              2.555 |                       7 |
| 2015-12-01 | R2 LRIV inherit_control_level->ok            |         45 |               38.736 |              2.495 |                      23 |
| 2016-02-01 | R2 LRIV ok->inherit_ok                       |         45 |                0.024 |             -2.834 |                      23 |
| 2017-03-01 | R2 LRIV inherit_ok->ok                       |         21 |                0.014 |             -2.31  |                      13 |
| 2017-04-01 | R2 LRIV ok->inherit_ok                       |         45 |                0.093 |             -2.952 |                      24 |
| 2018-01-01 | R2 LRIV inherit_ok->ok                       |         21 |              246.311 |              4.486 |                      12 |
| 2018-01-01 | R2 LRIV inherit_ok->ok                       |         45 |               32.042 |              2.076 |                      19 |
| 2018-06-01 | R2 LRIV inherit_ok->ok                       |         45 |               18.426 |              2.649 |                      45 |
| 2018-08-01 | R2 LRIV ok->inherit_ok                       |         45 |                0.035 |             -3.444 |                      34 |
| 2019-09-01 | R2 LRIV ok->inherit_ok                       |         21 |                0.021 |             -2.061 |                      14 |
| 2019-09-01 | R2 LRIV ok->inherit_ok                       |         45 |                0.047 |             -2.278 |                      19 |
| 2019-11-01 | R2 LRIV inherit_ok->ok                       |         21 |               57.113 |              2.039 |                      18 |
| 2019-11-01 | R2 LRIV inherit_ok->ok                       |         45 |               34.789 |              2.202 |                      26 |
| 2019-12-01 | R2 LRIV ok->inherit_ok                       |         21 |                0.006 |             -2.865 |                      20 |
| 2021-02-01 | R2 LRIV inherit_ok->ok                       |         21 |              240.713 |              4.885 |                      18 |
| 2021-02-01 | R2 LRIV inherit_ok->ok                       |         45 |              147.907 |              3.243 |                      20 |
| 2021-03-01 | R2 LRIV ok->inherit_ok                       |         45 |                0.048 |             -2.69  |                      24 |
| 2023-04-01 | R2 LRIV inherit_ok->ok                       |         21 |               54.722 |              2.423 |                      21 |
| 2023-04-01 | R2 LRIV inherit_ok->ok                       |         45 |               31.651 |              2.382 |                      35 |
| 2023-05-01 | R2 LRIV ok->inherit_ok                       |         21 |                0.018 |             -3.58  |                      20 |
| 2023-05-01 | R2 LRIV ok->inherit_ok                       |         45 |                0.024 |             -2.565 |                      31 |
| 2024-08-01 | R2 LRIV inherit_ok->ok                       |         21 |               16.471 |              2.076 |                      19 |
| 2024-09-01 | R2 LRIV ok->inherit_ok                       |         45 |                0.052 |             -2.173 |                      40 |
| 2015-10-01 | R2 NLLB inherit_ok->ok                       |         21 |               41.624 |              2.555 |                       7 |
| 2015-12-01 | R2 NLLB inherit_ok->ok                       |         45 |               38.736 |              2.495 |                      23 |
| 2016-02-01 | R2 NLLB ok->inherit_ok                       |         45 |                0.024 |             -2.834 |                      23 |
| 2017-03-01 | R2 NLLB inherit_ok->ok                       |         21 |                0.014 |             -2.31  |                      13 |
| 2017-07-01 | R2 NLLB ok->inherit_ok                       |         21 |                0.155 |             -2.3   |                       6 |
| 2018-01-01 | R2 NLLB inherit_ok->ok                       |         21 |              246.311 |              4.486 |                      12 |
| 2018-01-01 | R2 NLLB inherit_ok->ok                       |         45 |               32.042 |              2.076 |                      19 |
| 2018-08-01 | R2 NLLB ok->inherit_ok                       |         45 |                0.035 |             -3.444 |                      34 |
| 2019-09-01 | R2 NLLB nodata->inherit_ok                   |         21 |                0.021 |             -2.061 |                      14 |
| 2019-09-01 | R2 NLLB nodata->inherit_ok                   |         45 |                0.047 |             -2.278 |                      19 |
| 2019-11-01 | R2 NLLB inherit_ok->ok                       |         21 |               57.113 |              2.039 |                      18 |
| 2019-11-01 | R2 NLLB inherit_ok->ok                       |         45 |               34.789 |              2.202 |                      26 |
| 2019-12-01 | R2 NLLB ok->inherit_ok                       |         21 |                0.006 |             -2.865 |                      20 |
| 2021-02-01 | R2 NLLB inherit_ok->ok                       |         21 |              240.713 |              4.885 |                      18 |
| 2021-02-01 | R2 NLLB inherit_ok->ok                       |         45 |              147.907 |              3.243 |                      20 |
| 2021-03-01 | R2 NLLB ok->inherit_ok                       |         45 |                0.048 |             -2.69  |                      24 |
| 2023-03-01 | R2 NLLB ok->inherit_ok                       |         45 |               86.965 |              2.886 |                      28 |
| 2023-04-01 | R2 NLLB inherit_ok->ok                       |         21 |               54.722 |              2.423 |                      21 |
| 2023-04-01 | R2 NLLB inherit_ok->ok                       |         45 |               31.651 |              2.382 |                      35 |
| 2023-05-01 | R2 NLLB ok->inherit_ok                       |         21 |                0.018 |             -3.58  |                      20 |
| 2023-05-01 | R2 NLLB ok->inherit_ok                       |         45 |                0.024 |             -2.565 |                      31 |
| 2024-09-01 | R2 NLLB ok->inherit_ok                       |         45 |                0.052 |             -2.173 |                      40 |
| 2025-03-01 | R2 NLLB ok->inherit_ok                       |         21 |                0.017 |             -2.128 |                       4 |
| 2015-10-01 | R2 PGC inherit_ok->ok                        |         21 |               41.624 |              2.555 |                       7 |
| 2015-12-01 | R2 PGC inherit_ok->ok                        |         45 |               38.736 |              2.495 |                      23 |
| 2016-02-01 | R2 PGC ok->inherit_ok                        |         45 |                0.024 |             -2.834 |                      23 |
| 2017-07-01 | R2 PGC ok->inherit_ok                        |         21 |                0.155 |             -2.3   |                       6 |
| 2018-01-01 | R2 PGC inherit_ok->ok                        |         21 |              246.311 |              4.486 |                      12 |
| 2018-01-01 | R2 PGC inherit_ok->ok                        |         45 |               32.042 |              2.076 |                      19 |
| 2018-08-01 | R2 PGC ok->inherit_ok                        |         45 |                0.035 |             -3.444 |                      34 |
| 2019-09-01 | R2 PGC ok->inherit_ok                        |         21 |                0.021 |             -2.061 |                      14 |
| 2019-09-01 | R2 PGC ok->inherit_ok                        |         45 |                0.047 |             -2.278 |                      19 |
| 2019-11-01 | R2 PGC inherit_ok->ok                        |         21 |               57.113 |              2.039 |                      18 |
| 2019-11-01 | R2 PGC inherit_ok->ok                        |         45 |               34.789 |              2.202 |                      26 |
| 2019-12-01 | R2 PGC ok->inherit_ok                        |         21 |                0.006 |             -2.865 |                      20 |
| 2020-06-01 | R2 PGC ok->inherit_ok                        |         21 |                0.065 |             -2.595 |                      10 |
| 2020-07-01 | R2 PGC inherit_ok->ok                        |         45 |                0.105 |             -2.407 |                      20 |
| 2021-02-01 | R2 PGC inherit_ok->ok                        |         21 |              240.713 |              4.885 |                      18 |
| 2021-02-01 | R2 PGC inherit_ok->ok                        |         45 |              147.907 |              3.243 |                      20 |
| 2021-03-01 | R2 PGC ok->inherit_ok                        |         45 |                0.048 |             -2.69  |                      24 |
| 2023-03-01 | R2 PGC ok->inherit_ok                        |         45 |               86.965 |              2.886 |                      28 |
| 2023-04-01 | R2 PGC inherit_ok->ok                        |         21 |               54.722 |              2.423 |                      21 |
| 2023-04-01 | R2 PGC inherit_ok->ok                        |         45 |               31.651 |              2.382 |                      35 |
| 2023-05-01 | R2 PGC ok->inherit_ok                        |         21 |                0.018 |             -3.58  |                      20 |
| 2023-05-01 | R2 PGC ok->inherit_ok                        |         45 |                0.024 |             -2.565 |                      31 |
| 2024-10-01 | R2 PGC ok->inherit_ok                        |         45 |                0.042 |             -3.814 |                      32 |
| 2025-03-01 | R2 PGC ok->inherit_ok                        |         21 |                0.017 |             -2.128 |                       4 |
| 2015-10-01 | R2 SNB inherit_ok->ok                        |         21 |               41.624 |              2.555 |                       7 |
| 2015-12-01 | R2 SNB inherit_ok->ok                        |         45 |               38.736 |              2.495 |                      23 |
| 2016-02-01 | R2 SNB ok->inherit_ok                        |         45 |                0.024 |             -2.834 |                      23 |
| 2017-03-01 | R2 SNB inherit_ok->ok                        |         21 |                0.014 |             -2.31  |                      13 |
| 2017-07-01 | R2 SNB ok->inherit_ok                        |         21 |                0.155 |             -2.3   |                       6 |
| 2018-01-01 | R2 SNB inherit_ok->ok                        |         21 |              246.311 |              4.486 |                      12 |
| 2018-01-01 | R2 SNB inherit_ok->ok                        |         45 |               32.042 |              2.076 |                      19 |
| 2018-08-01 | R2 SNB ok->inherit_ok                        |         45 |                0.035 |             -3.444 |                      34 |
| 2019-09-01 | R2 SNB ok->inherit_ok                        |         21 |                0.021 |             -2.061 |                      14 |
| 2019-09-01 | R2 SNB ok->inherit_ok                        |         45 |                0.047 |             -2.278 |                      19 |
| 2019-11-01 | R2 SNB inherit_ok->ok                        |         21 |               57.113 |              2.039 |                      18 |
| 2019-11-01 | R2 SNB inherit_ok->ok                        |         45 |               34.789 |              2.202 |                      26 |
| 2019-12-01 | R2 SNB ok->inherit_ok                        |         21 |                0.006 |             -2.865 |                      20 |
| 2020-06-01 | R2 SNB ok->inherit_ok                        |         21 |                0.065 |             -2.595 |                      10 |
| 2020-07-01 | R2 SNB inherit_ok->ok                        |         45 |                0.105 |             -2.407 |                      20 |
| 2021-02-01 | R2 SNB inherit_ok->ok                        |         21 |              240.713 |              4.885 |                      18 |
| 2021-02-01 | R2 SNB inherit_ok->ok                        |         45 |              147.907 |              3.243 |                      20 |
| 2021-03-01 | R2 SNB ok->inherit_ok                        |         45 |                0.048 |             -2.69  |                      24 |
| 2015-10-01 | R2 VGZ inherit_ok->ok                        |         21 |               41.624 |              2.555 |                       7 |
| 2015-12-01 | R2 VGZ inherit_ok->ok                        |         45 |               38.736 |              2.495 |                      23 |
| 2016-02-01 | R2 VGZ ok->inherit_ok                        |         45 |                0.024 |             -2.834 |                      23 |
| 2017-07-01 | R2 VGZ ok->inherit_ok                        |         21 |                0.155 |             -2.3   |                       6 |
| 2018-01-01 | R2 VGZ inherit_ok->ok                        |         21 |              246.311 |              4.486 |                      12 |
| 2018-01-01 | R2 VGZ inherit_ok->ok                        |         45 |               32.042 |              2.076 |                      19 |
| 2018-08-01 | R2 VGZ ok->inherit_ok                        |         45 |                0.035 |             -3.444 |                      34 |
| 2019-09-01 | R2 VGZ ok->inherit_ok                        |         21 |                0.021 |             -2.061 |                      14 |
| 2019-09-01 | R2 VGZ ok->inherit_ok                        |         45 |                0.047 |             -2.278 |                      19 |
| 2019-11-01 | R2 VGZ inherit_ok->ok                        |         21 |               57.113 |              2.039 |                      18 |
| 2019-11-01 | R2 VGZ inherit_ok->ok                        |         45 |               34.789 |              2.202 |                      26 |
| 2019-12-01 | R2 VGZ ok->inherit_ok                        |         21 |                0.006 |             -2.865 |                      20 |
| 2020-06-01 | R2 VGZ ok->inherit_ok                        |         21 |                0.065 |             -2.595 |                      10 |
| 2020-07-01 | R2 VGZ inherit_ok->control_level             |         45 |                0.105 |             -2.407 |                      20 |
| 2021-02-01 | R2 VGZ inherit_ok->ok                        |         21 |              240.713 |              4.885 |                      18 |
| 2021-02-01 | R2 VGZ inherit_ok->ok                        |         45 |              147.907 |              3.243 |                      20 |
| 2021-03-01 | R2 VGZ ok->inherit_ok                        |         45 |                0.048 |             -2.69  |                      24 |
| 2023-03-01 | R2 VGZ ok->inherit_ok                        |         45 |               86.965 |              2.886 |                      28 |
| 2023-04-01 | R2 VGZ inherit_ok->ok                        |         21 |               54.722 |              2.423 |                      21 |
| 2023-04-01 | R2 VGZ inherit_ok->ok                        |         45 |               31.651 |              2.382 |                      35 |
| 2023-05-01 | R2 VGZ ok->inherit_ok                        |         21 |                0.018 |             -3.58  |                      20 |
| 2023-05-01 | R2 VGZ ok->inherit_ok                        |         45 |                0.024 |             -2.565 |                      31 |
| 2024-08-01 | R2 VGZ inherit_ok->ok                        |         21 |               16.471 |              2.076 |                      19 |
| 2024-10-01 | R2 VGZ ok->inherit_ok                        |         45 |                0.042 |             -3.814 |                      32 |
| 2025-03-01 | R2 VGZ ok->inherit_ok                        |         21 |                0.017 |             -2.128 |                       4 |
| 2025-10-01 | R2 VGZ ok->inherit_ok                        |         21 |               33.265 |              2.273 |                      32 |

* Runs used: {'variant': 'nc', 'prod': 'prod', 'prod_S13': 'prod_S13', 'truth': {'2010': {'S9': 'r3t2010_nc'}, '2011': {'S15': 'r3t2011_nc', 'S14': 't2011_S14_nc', 'S13': 't2011_S13_nc'}, '2012': {'S15': 'r3t2012_nc', 'S14': 't2012_S14_nc', 'S13': 't2012_S13_nc'}, '2013': {'S15': 't2013_S15_nc', 'S14': 't2013_S14_nc', 'S13': 't2013_S13_nc'}}}

**Correction 2026-10-07 (count columns only, no rate changed):** `monthly_total.csv` `n_quiet_famdays` and `monthly_family.csv` `quiet_days` counted only family-days that had a candidate. Renamed to `n_quiet_famdays_with_candidates` / `quiet_days_with_candidates`; new `n_quiet_days` and `n_quiet_famdays` (= quiet days × 115 families; quiet is a network-level day flag) and `quiet_days` (quiet days in the month). Totals 2014–2026: 231,035 quiet family-days (was shown as 91,769). Rates are sums divided by calendar days and were never affected. Script fixed (`svi_ref_rates.py`, backup `.bak_pre_quietcount`); CSV backups `*.bak_pre_quietcount`.
