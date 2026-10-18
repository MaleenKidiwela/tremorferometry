# SH_KMRFIX: the NC.KMR fix for the Shelly (2025) Mendocino forward extension (2023-12 → 2026-09-28)

2026-10-01. CPU only, saved files only (the raw continuous data for this run is deleted). No downloads, nothing deleted or overwritten, no git. Scripts: `scripts/sh_kmrfix_t1.py`, `sh_kmrfix_t1b.py`, `sh_kmrfix_cat.py`, `sh_kmrfix_conf.py`, `sh_kmrfix_rates.py`, `sh_kmrfix_windows.py`, `sh_kmrfix_out.py`. Outputs are in `results/sh_kmrfix/` unless a path says otherwise.

## Plain summary

**1. What happened to KMR before the swap. The saved data do not fully match the earlier description.**
- The earlier note said HHE and HHN were stuck from 2025-07-26/27. The cut windows show something different (`t1_kmr_daily_health.csv`):
  - From **2025-07-25** the DC level of all three STS-2 channels drifts by millions of counts, as if the masses were moving off centre.
  - **HHE** is stuck at a constant −8,980,480 counts from **2025-08-02/03**.
  - **HHN** is never constant. It settles at about −7.05 M counts but keeps normal-looking noise (100–200 counts std).
  - **Z** re-centres on 2025-08-06.
  - **KMR has no data at all from 2025-08-26 to 2025-09-17.** On those days segment A already ran as a 16-station network (nsta 16, `mf_A_diag.csv`).
- **The cleaning masked HHE as dead from 2025-08-03 and counted HHN and Z as live** (`t1_kmr_prep_flags.csv`: dead_ch = 1, ref of HHE = 0). The sh_v2 dead-channel rule did its job on HHE.
- **HHN and Z still carry a weakened LFE signal** (`t1_kmr_channel_stacks_summary.csv`, median over 27 families). Each value is the cc of a direct stack against the old template; the scrambled control (each event's slice taken at a random time) is in brackets.

  | period | HHE | HHN | HHZ |
  |---|---|---|---|
  | healthy, 2025-06-01 → 07-24 | 0.83 (0.23) | 0.85 (0.23) | 0.75 (0.21) |
  | drift, 07-25 → 08-02 | 0.62 (0.18) | 0.65 (0.20) | 0.55 (0.21) |
  | HHE stuck, 08-03 → 08-25 | n/a | 0.60 (0.19) | 0.43 (0.20) |

  Single-event cc fell about 3× (HHN median 0.23 → 0.08, Z 0.16 → 0.05). So KMR was failing, not dead, and it still contributed weak true signal.
- **The effect on the catalogue is small, so no correction was applied** (`t1b_effect_S1.csv`). Method: subtract KMR's per-channel cc from the network sum at every ≥ 8×MAD candidate, shrink the daily MAD by that channel's noise, and re-apply the S1 rules.
  - Dropping the live KMR horizontals changes the S1 count over 07-25 → 08-25 by **+1.2 %** (1306 → 1322): −2.6 % in the drift days, +2.8 % in the HHE-stuck days. Per family the median is +1.5 % (range −11 % … +14 %).
  - Dropping all of KMR gives +0.8 %.
  - These changes are within the error of the approximation.
  - The recompute without any change reproduces the stored catalogue exactly (1306 = 1306).

**2. Forward catalogue v2** (`catalogue_S1_9x_v2.parquet`).
- **Before 2025-09-17:** cont + A, unchanged. 24,902 events, identical to v1.
- **From 2025-09-17:** the 16-station no-KMR run (mf_NK) with the same S1 rules.
  - Check: cat_NK_B (7,289) and cat_NK_C (8,139) are reproduced row for row (`t2_checks.json`).
- **After the swap, v2 has 15,718 S1 events against 13,664 in v1 (1.150×).** The per-family median is 1.137× (range 0.82 for family 9, n = 65 → 53, up to 1.27). Total v2 = 40,620 events.
- **Each event has a det_id_v2 and a mapping to a cut window** (`results/cut/sh_ext/detections_v2.parquet`):
  - **All S1 events:** 97.4 % exact (same family, same peak); 1.6 % same family within 2 s; 0.2 % other family; 298 (0.7 %) unmapped. Unmapped means no cut window, flag `window_unmapped`.
  - **Post-swap S1 events:** 97.6 % map to the same family.
  - **All post-swap ≥ 8×MAD candidates (85,299):** 74.5 % exact, 6.8 % same family, 10.6 % other family, 8.0 % unmapped.

**3. Confidence v2 (KMR removed from every feature): it costs about 0.01 AUC and stays calibrated** (`t3_eval.csv`).
- Design: same as SH_CONF. HGB + isotonic calibration, the stored 5 day-folds, the same 132 fit days 2018–2023, the same labels.
- Training rows had their no-KMR features rebuilt from saved per-candidate values.
  - Check on the post-swap days, against the real 16-station run at the same peaks (65,049 pairs, `t3_network_approx_check.json`):
    - ccsum and mean cc are reproduced exactly (ratio 1.000).
    - net_mad: median ratio 1.000, 10–90 % range 0.99–1.05, r = 0.990. Without the fix it was 0.978 and r = 0.981.
    - At the 9× threshold, 2,541 pairs pass on the rebuilt value only and 344 on the real one only, out of 43,760 that pass on both. So the rebuilt net_mad is slightly high at the top end.
  - Site and PNET features are rebuilt exactly when KMR is included. The site z values of the combined KBU+B049 site are reproduced at r = 0.996–0.997 (`t3_reconstruction_check.csv`).

  | eval set | label | original FULL AUC | v2 AUC | v2 calibration slope | v2 matched/ΣP | control P > 0.5 |
  |---|---|---|---|---|---|---|
  | 2024 held-out (34 days) | same | 0.867 | **0.856** | 0.97 | 0.92 | 0.0 % |
  | 2024 held-out | any | 0.869 | **0.857** | 1.05 | 0.95 | 0.0 % |
  | forward overlap (244 days) | same | 0.884 | **0.870** | 0.99 | 0.94 | 0.0 % |
  | forward overlap | any | 0.883 | **0.869** | 1.05 | 0.96 | 0.3 % |
  | fit, out-of-fold | same | 0.857 | 0.829 | 1.01 | 1.00 | 0.1 % |

- **Post-swap events:**
  - Network features come from mf_NK itself. Site, PNET and HF features come from the cut candidate at the same peak (same family, |dt| ≤ 0.1 s).
  - The rest, 19,524 of 85,299 candidates (700 of 15,718 S1 events), including all unmapped ones, get P from a **fallback model** fitted the same way without site/PNET/HF features. Its AUC is 0.853 / 0.867 (same; 2024 held-out / forward overlap) and its matched/ΣP is 0.93 / 0.94.
  - The fallback carries 16.6 % of the post-swap ΣP(same) (`t5_outputs_summary.json`).
- **Two series are given** (`results/cut/sh_ext/conf_P_v2_by_det_id.parquet`):
  - `P_*_v2a`: the original model up to 2025-07-24, v2 after.
  - `P_*_v2`: the v2 model throughout.
  - They are identical from 2025-07-25 on. Before that, the v2 model gives 1.020× the original ΣP(same) (the model-only effect).

**4. Rates.**
- **After the swap (2025-09-17 → 2026-09-28), v2 / v1** (`t4_v1_vs_v2_post_swap.csv`; ratios are of unscaled counts):
  - S1 count: **1.150×**, per-family median 1.137.
  - ΣP(same): **1.227×**, per-family median 1.226, range 1.00–1.45. About 2 % of this comes from the model change and the rest from the 16-station run.
  - ΣP(any): 1.196×.
  - Join-year scaled ΣP(same): per-family median 1.191×.
  - The S1 rise fits the expected 11–18 %. The ΣP rise is larger. The post-swap 16-station run finds more and stronger candidates than the noisy-KMR run: on the same days, ΣP_v2 on the 16-station candidates is 1.106× ΣP_v2 on the with-KMR candidates (`t3_consistency_NK_vs_BC.csv`).
- **Network mean rate, scaled to his at the join year** (events/day, S1 v1 → v2):

  | period | S1 v1 → v2 | ΣP v1 → v2 |
  |---|---|---|
  | 2024 | 38.8 | 38.8 |
  | 2024-12 → 2025-07-24 | 44.5 | 45.2 → 45.6 |
  | 2025-07-25 → 09-16 | 45.3 | 46.8 → 48.7 |
  | 2025-09-17 → 2026-09-28 | 39.5 → **45.6** | 39.8 → **48.3** |

  So the apparent drop after the swap in v1 was the KMR artefact. In v2 the post-swap level is about the same as 2024-12 → 2025-07.
- **Steps** (`t4_step_checks.csv`; daily network totals, unscaled; ±21 and ±45 days; the same calendar windows in 2024 give the seasonal comparison).
  - **2025-07-25:**
    - 2025: S1 v2 ratio 0.96 (21 d) and 0.94 (45 d).
    - 2024, same calendar windows: 0.82 and 1.23.
    - No clear step.
  - **2025-09-17:**
    - 2025: S1 v2 ratio 1.42 (21 d) and 1.22 (45 d). Medians 1.32 / 1.00; quiet half 1.07 / 1.21.
    - 2024, same calendar windows: 0.83 / 0.80.
    - A tremor-active period follows 09-17 (figure, lower panel), so the rise cannot be cleanly separated from activity.
    - From 2025-08-26 to 2026-09-28 the detector setup is the same, 16 stations with the same templates: KMR was absent in segment A before 09-17 and dropped in mf_NK after. So v2 has no instrument step at 09-17 by construction.
  - **2025-08-26** (KMR lost):
    - 2025: S1 v2 ratio 1.10 (21 d) and 0.94 (45 d).
    - 2024, same calendar windows: 1.09 / 1.08.
  - **The step that remains is losing a working KMR.** Estimated from the rebuilt no-KMR features of the old-epoch candidates, a 16-station network would have found **0.937×** as many S1 events as the 17-station run over 2023-12 → 2025-07-24, and 0.928× over the join year. So the 16-station part of v2 (from 2025-08-26) is probably still **about 6–7 % low** against the 17-station join year it is scaled to. This was not corrected.
    - It is approximate. It uses only the saved ≥ 8×MAD candidates, and the rebuilt net_mad is slightly high at the top end, so the true deficit may be a little larger.
- Figures:
  - `figures/monthly_rate_total_v1_v2.png`: monthly total since 2018, and a daily close-up with the KMR dates.
  - `figures/monthly_rate_per_family_v1_v2.png`: per family.

**5. Cut windows.**
- The HDF5 files are not rewritten. `results/cut/sh_ext/window_flags_v2.parquet` has one row per window (11,590,888 rows: det_id, seed_id, h5_file, h5_start, t0).
- **`bad_for_dvv`** (183,038 windows):
  - KMR HHE/HHN from **2025-07-25**: 14,276 windows. That is one day earlier than the 07-26 in the brief, because the drift starts on 07-25.
  - All KMR channels from 2025-09-17: 168,762 windows.
- **`kmr_z_caution`** marks KMR HHZ from 2025-07-25 to 2025-09-16 (7,138 windows): mass drift and weakened signal. Decide whether to use these.
- The v2 det_id → window mapping is in `detections_v2.parquet`: `det_id` is the window id, −1 if unmapped, and `map_dt_s` is the time difference to shift by when reading the window.

**What could not be done without the raw data.**
- A true 16-station re-run before 2025-08-26 (needed to measure the loss-of-KMR step exactly). Only estimated, ~0.93×.
- New candidates that would rise above 8×MAD only without KMR, before the swap.
- Per-site PNET medians at random times: pnet_rel_max is approximated by shifting it by the change in pnet_max.
- The scrambled control for the 16-station run after the swap: the with-KMR control rows were used, with KMR removed approximately.
- Exact per-channel live counts at each candidate: nominal channel counts were used.
- HF features where the family's HF station was KMR (2,357 rows): set missing.

## Files

- `results/sh_kmrfix/`:
  - Task 1: `t1_*`, `t1b_*`.
  - Catalogue and candidates: `catalogue_S1_9x_v2.parquet`, `candidates8x_v2.parquet`, `t2_checks.json`.
  - Confidence: `conf_v2_rows.parquet` (every model row with v2 features' key columns and P), `conf_v2_BCcheck_rows.parquet`, `models_v2.joblib`, `t3_eval.csv`, `t3_network_approx_check.json`, `t3_reconstruction_check.csv`, `t3_consistency_NK_vs_BC.csv`, `t3_log.txt`.
  - Rates: `t4_daily_per_family.csv.gz` (date × family: his_n, S1_v1/v2, sumP_v1/v2a/v2b, raw and scaled), `t4_monthly_per_family.csv`, `t4_monthly_total.csv`, `t4_join_scale_his_over_ours.csv`, `t4_v1_vs_v2_post_swap.csv`, `t4_step_checks.csv`.
  - `t5_outputs_summary.json`, `figures/`.
- `results/cut/sh_ext/` (new files):
  - `detections_v2.parquet`.
  - `conf_P_v2_by_det_id.parquet`: keyed by det_id_v2, with det_id for the windows.
  - `window_flags_v2.parquet`.
- The originals (`conf_P_by_det_id.parquet`, `detections.parquet`, `results/sh_ext2026/*`, `results/sh_conf/*`) are untouched.

## 7. Correction for losing a working KMR (added 2026-10-02 by the main session)
- Old-epoch candidates (2023-12 → 2025-07-24) with KMR removed from the network sum: S1 count 21,268 vs 22,693 with KMR (0.937×); join year 12,173 vs 13,118 (0.928×). ΣP from the v2 model (fitted without KMR): 1.020× / 1.016× the original — so **ΣP_v2b has no loss-of-KMR step** and needs no correction.
- S1 counts only: from 2025-08-26 (KMR's last day with data was 2025-08-25) divided by 0.928 (band 0.918–0.937). Per-family factors are noisy (median 0.965, range 0.83–1.22; `sh_kmrfix/t6_kmrloss_factor_by_family.csv`), so the network-wide factor is used.
- Outputs: `sh_kmrfix/t6_daily_per_family_v3.csv.gz` (S1_v3_scaled, lo/hi), `sh_kmrfix/t6_monthly_total_v3.csv`.
- Recommended rate series: ΣP_v2b (scaled at the join year). S1_v3 is the count-based cross-check.
