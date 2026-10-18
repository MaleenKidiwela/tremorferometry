# Robustness of the per-episode rate result — pre-registration (written before any re-run)

Written 2026-10-06 ~8:55 pm PT. Parent run: `../RESULTS.md`, rules `../PLAN.md` (unchanged here). Outputs only in this folder.
Script: `lfe_extension_wa/scripts/episode_rates.py` with `EPR_OUT=<this folder>` (the first-run outputs are not touched).

## Variants (all at the PLAN.md final episode rule, k-free share test, same neighbourhood / bootstrap / permutation settings)
- V0 `rerun_075`: full chain re-run on the updated inputs (`rates_daily_v2.parquet` with zero rows, DEVIATIONS §29/§30); headline
  cut 0.75. Expected: identical to the first run (the episode check counts its own days from `conf/P`, k-free).
- V1 `cut_050`: N = data - main control at P >= 0.5 (WA own_pnet n050 - time-reversed; SVI P_same >= 0.5 - scrambled control).
- V2 `cut_085`: same at P >= 0.85.
- V3 `flag_excl_075`: at 0.75, drop SVI episodes 2017-03 and 2025-10 (id 2025-09b1) entirely; and remove from their own family's
  test set every (family, episode) cell that is (a) labelled station-coincident in V0 (either the principal-carrier label or, WA,
  the broad-flag label) or (b) has mean live carriers per window day < 50 % of the family's median over its MAJOR episodes of
  that block. Removed cells keep their counts in the neighbours' totals T_e (the share test's normalisation is unchanged).
Each variant reports, per block (NW E3 / E2 / E1, SVI fixed-S13, SVI eff S15/S14/S13): tested, over-dispersed (BH q < 0.05),
consistent, exposure-explained, median chi2/df, residual-correlation bins vs distance, distance and carrier-overlap coefficients
with permutation p (1,000).

## Pre-registered criterion
The finding "most NW families over-dispersed, SVI closer to consistent, scatter shared within ~20 km, not station-driven"
**holds** if in every variant, for the primary blocks (NW E3 and SVI fixed-S13):
(i) the distance coefficient is negative with permutation p < 0.01 in both regions, and
(ii) the NW E3 over-dispersed fraction (over-dispersed / tested) stays > 50 %.
Reported alongside but not part of the criterion: SVI over-dispersed fraction < NW's ("closer to consistent"), the carrier
coefficient (the "not station-driven" part is read as: distance term survives; a significant carrier term is reported as
partial station coupling, as in the first run), and the secondary blocks (E1/E2, SVI eff).

---------------------------------------------------------------------------------------------------------------------------
# Results (added after the runs, 2026-10-06 ~9:55 pm PT)

Combined table: `robustness_table.csv` (rows = variant × block). Source: `consistency_summary.json`,
`consistency_families.csv`, `consistency_rows.csv` and `residual_pair_correlations.csv` (all with a `variant` column), plus
`sanity_checks.json`. Log: `episode_rates.log`.

## Verdict: the pre-registered criterion HOLDS in all four variants

| variant | NW E3 over-dispersed / tested (fraction) | NW E3 median χ²/df | NW E3 distance coef (p) | NW E3 carrier coef (p) | SVI-S13 over-dispersed / tested | SVI median χ²/df | SVI distance coef (p) | SVI carrier coef (p) |
|---|---|---|---|---|---|---|---|---|
| V0 re-run 0.75 | 58/89 (0.65) | 4.03 | −0.111 (0.001) | +0.230 (0.001) | 27/75 (0.36) | 1.81 | −0.517 (0.001) | −0.035 (0.48) |
| V1 P ≥ 0.5 | 93/112 (0.83) | 6.43 | −0.140 (0.001) | +0.179 (0.001) | 54/93 (0.58) | 3.05 | −0.497 (0.001) | −0.025 (0.52) |
| V2 P ≥ 0.85 | 40/60 (0.67) | 3.26 | −0.137 (0.001) | +0.238 (0.005) | 22/56 (0.39) | 1.71 | −0.478 (0.001) | +0.061 (0.25) |
| V3 flagged removed, 0.75 | 60/89 (0.67) | 3.95 | −0.110 (0.001) | +0.230 (0.001) | 25/75 (0.33) | 1.83 | −0.547 (0.001) | −0.042 (0.44) |

p = permutation p from 1,000 draws (0.001 is the floor). Exposure-explained families: NW E3 1–3 and SVI 1–2 in every variant.

Mean residual correlation by distance bin (< 3 / 3–6 / 6–10 / 10–20 / 20–40 km):

| variant | NW E3 | SVI fixed-S13 |
|---|---|---|
| V0 | 0.31 / 0.11 / 0.16 / 0.06 / −0.03 | 0.73 / 0.49 / 0.39 / 0.18 / −0.05 |
| V1 | 0.29 / 0.19 / 0.15 / 0.08 / −0.03 | 0.61 / 0.48 / 0.39 / 0.20 / −0.06 |
| V2 | 0.42 / 0.15 / 0.15 / 0.07 / −0.04 | 0.68 / 0.49 / 0.38 / 0.20 / −0.07 |
| V3 | 0.31 / 0.11 / 0.16 / 0.06 / −0.03 | 0.73 / 0.51 / 0.41 / 0.20 / −0.05 |

## 1. Re-run on the updated inputs (V0): identical to the first run
- `rates_daily_v2.parquet` now holds 810,265 covered NW rows (603,874 before), because the zero rows were added.
- n_trem30 still agrees on every row (0 mismatches). My counts from `conf/P` still equal `rates_daily_v2` at every cut and
  series (0 differing rows).
- Per-(family, episode) tables (NW, SVI fixed-S13, SVI eff) and `consistency_families.csv` are cell-for-cell identical to
  `../`. The only exception is the order of entries inside the NW `station_changes` text: `core_set_changes_v2.csv` was
  regenerated, but all station-change flags are unchanged. k_E3 is unchanged.
- Sanity checks are unchanged: (a) pooled 1.094, FAIL on the Poisson bar; (b) NW pass, SVI fail in 3 tiny episodes;
  (c) NW E3 0.918, SVI 0.814.

## 2. Cut-off sensitivity (V1, V2)
- **NW over-dispersed fraction:** stays above 50 % at both cut-offs (0.83 at 0.5, 0.67 at 0.85).
- **Distance term:** stays negative with p = 0.001 in both regions.
- **SVI "closer to consistent":**
  - It holds only relative to NW.
  - At P ≥ 0.5 a majority of SVI families are over-dispersed (54 of 93, median χ²/df 3.05).
  - At 0.5 the SVI scrambled control is about 8–9 % of data (prod 2,551 / 30,411; S13 run 1,473 / 16,257 —
    `counts_svi_check.json`/`P_by_candidate`), against < 1 % at 0.75. So part of the extra scatter at 0.5 is
    false-detection noise.

## 3. Flagged cells removed (V3)
- **SVI:** episodes 2017-03 and 2025-10 dropped (37 MAJOR episodes left). The < 50 % live-carrier rule removed no cell: the
  2025-10 station gaps cut the mean live carriers from about 8.7 to 7.2, not below half.
- **NW E3:** 11 station-coincident cells removed (E2: 11; E1: 22, plus 16 low-carrier cells).
- **SVI eff S15:** 6 cells removed. None elsewhere in SVI.
- **Effect:** small. NW E3 60/89 over-dispersed. SVI 25/75, slightly more consistent without the two episodes.

## Which parts hold and which do not
- **Holds in every variant:**
  - Most NW E3 families are over-dispersed (0.65–0.83).
  - SVI has a lower over-dispersed fraction than NW (0.33–0.58).
  - Exposure explains almost none of it.
  - The extra scatter is shared by families within about 10–20 km and decays to about 0 by 20–40 km. The distance
    coefficient is p = 0.001 in both regions in all four variants.
- **Does not hold as an unqualified statement:**
  - **"Not station-driven" holds only in the pre-registered weak sense: the distance term survives.** In NW E3 the carrier
    term (overlap of template weights) is also significant in every variant (+0.18 to +0.24, p ≤ 0.005). Removing the
    station-coincident cells does not change it. So the NW over-dispersion is partly station-coupled in every variant.
    In SVI fixed-S13 the carrier term is never significant (p 0.25–0.52).
  - **"SVI closer to consistent" means relative to NW only.** At P ≥ 0.5 most SVI families are over-dispersed too.
- **Secondary blocks (not in the criterion):**
  - NW E1 loses the distance term at P ≥ 0.85 (p = 0.47) and in V3 (p = 0.09).
  - NW E2's distance term is p = 0.03 at 0.75.
  - These blocks have only 8 and 13 episodes. They are reported, not used.
- **Caveat carried over from the first run:** the correlation at < 3 km may partly be co-detection between co-located
  families. That test (merging families within 3 km) has not been run.
