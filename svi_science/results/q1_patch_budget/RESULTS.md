# Q1 results: why each big ETS slips the fault in a different pattern

Run 2026-10-07, ~2:30–3:00 pm PT. Rules frozen in `PLAN.md` before any Q1 number was computed; changes listed in its §9.
Scripts: `svi_science/scripts/q1_build.py`, `q1_tests.py`, `q1_scrambled.py`, `q1_figures.py` (CPU only, < 2 GB).
No Merlin review of this run. All results are provisional.

## The answer in plain words

- **What decides which patches are strong in a big ETS is mainly how far the slip front gets.**
  - Patches inside the stretch that the PNSN tremor covered in that episode fire far more than usual. Patches outside
    it are nearly empty. Tremor is an independent catalogue.
  - This one fact explains about half of the patch-to-patch differences within an episode (R² 0.53 SVI, 0.50 NW).
    Adding all the front terms raises this to 0.62 / 0.61.
  - The shorter the stretch the front covers, the more uneven the episode. In SVI the ranking is perfect over all 10 big
    ETS (Spearman −1.00).
- **In Vancouver Island the pattern also depends on which end the front came in from.** This was not the case in
  Washington.
  - SVI episodes running the same way look alike (R² 0.34). Position along the front's path explains 0.40.
  - Neither passes the strict p < 0.01 bar on its own (p 0.026, 0.016). Both largely repeat the reach effect.
- **The "slip budget" idea is not supported.**
  - A patch that has gone longer since it last slipped does NOT reliably fire more. The slope is positive in both
    regions but within chance: order-shuffle p 0.073 (SVI) and 0.11 (NW).
  - It explains only 5 % (SVI) and 2 % (NW) of the differences.
  - It flips with the activity threshold: at "≥ 1 event" the slope is about zero in both regions.
- **One budget-like signal, only in Vancouver Island 2015–2026: turn-taking between consecutive big ETS.**
  - A patch that was strong in the previous big ETS tends to be weak in the next one, and the reverse.
  - Slope −1.11 per SD, p 0.001, R² 0.36. It survives the per-patch shuffle and the fixed 13-station run (p 0.008).
  - It is absent in Washington (+0.14, p 0.69).
  - It is absent in Bostock's 2003–2013 SVI catalogue (+0.05, p 0.73).
  - It adds little beyond the front terms: +0.06 R², p 0.07.
  - The likely reading, consistent with S3: SVI fronts in 2015–26 often alternated their entry end and stopped halfway.
    "Strong here last time, weak here now" is then mostly "the front came from the other side and stopped".
  - We could not show that the front goes where slip is overdue. Time since last active vs whether the front reached
    the patch: p 0.033 SVI, 0.09 NW, both fail.
- **Detectability and episode size do not explain it.**
  - Station availability has a statistically clear but tiny effect (R² 0.008 SVI, 0.002 NW).
  - The LFE count, the front speed and the station count do not change how uneven an episode is.
  - In SVI, larger GNSS offsets go with more even episodes (ρ −0.76, p 0.017, misses the bar). This is the same story:
    big episodes are the ones that sweep the whole box.

## Pre-registered verdicts

The test is the slope of the residual share (log obs/expected, leave-one-out patch share) on the predictor (per SD), with
episode fixed effects. "Order shuffle" means the episodes' patch patterns are reassigned to the other episode dates, which
removes any memory and keeps each episode's shares summing to one. 2,000 draws.

| test (predicted sign) | SVI 2015–26 (10 big ETS × 12 groups) | NW E3 2015–26 (14 big ETS × 9 bins) | verdict |
|---|---|---|---|
| H1 days since last active (+) | +0.46, p 0.073 (per-patch shuffle 0.042), R² 0.05 | +0.32, p 0.11 (per-patch 0.026), R² 0.02 | **NOT SUPPORTED** |
| H2 previous big-ETS residual, ≥ 90 d (−) | **−1.11, p 0.001** (per-patch 0.001; fixed-S13 −0.87, p 0.008), R² 0.36 | +0.14, p 0.69; fixed-9 +0.17, p 0.70 | **PARTIAL (SVI only)** |
| H3 cumulative past residual (−) | −0.32, p 0.002, but per-patch 0.083 and fixed-S13 p 0.55; R² 0.01 | +0.16, p 0.92 | **PARTIAL by the letter (SVI headline only); fails SVI's own checks** |
| H4 budget adds beyond front terms | ΔR² +0.062 (0.62 → 0.68), p 0.070 | ΔR² +0.028 (0.61 → 0.64), p 0.68 | **FAIL** |
| K1 inside tremor extent (+) | +1.37, p 0.0005, R² 0.53 | +0.91, p 0.0005, R² 0.50 | **PASS both** |
| K1 km outside extent (−) | −1.18, p 0.0005, R² 0.39 | −0.91, p 0.0005, R² 0.49 | **PASS both** |
| K2 same-direction pattern (+) | +1.01, p 0.026, R² 0.34 | +0.40, p 0.074, R² 0.11 | FAIL (SVI hint) |
| K3 position along front path (two-sided) | −1.09 (entry end stronger), p 0.016, R² 0.40 | +0.10, p 0.58 | FAIL (SVI hint) |
| D1 station availability (+) | +0.45, p 0.0055, R² 0.008 | p 0.0005, R² 0.002 | passes the rule; explains < 1 % |
| L2 time since last active → front reaches patch (+) | p 0.033 | p 0.092 | FAIL |
| E1 episode unevenness vs tremor extent width | **ρ −1.00, p 0.0002** | ρ −0.60, p 0.025 | SVI PASS, NW fails the bar |
| E1 vs LFE count / front speed / station count / GNSS D | −0.50 / −0.49 / −0.15 / −0.76 (p 0.15 / 0.16 / 0.67 / 0.017) | −0.20 / −0.04 / +0.36 | FAIL |

Tremor direction labels were used for K2/K3. With the LFE direction labels from S3 the results are the same: SVI K2
p 0.026, K3 p 0.011; NW p 0.27 / 0.15.

**Calibration.** 200 simulated records with no memory (random patch counts drawn from the long-term shares, real episode
sizes and dates) were put through the same pipeline. The order-shuffle test fired at p < 0.01 in 0–1.5 % of them (H1
1.5 %, H2 0.5 % / 0 %, H3 1 % / 0.5 %; SVI / NW). The mean slopes were about 0 (H2 +0.014 / +0.009). The leave-pair-out
residuals therefore do not manufacture the SVI turn-taking signal (`tests.json` → calibration).

## Robustness of the budget tests (slope, order-shuffle p; `tests_table.csv`)

| variant | SVI H1 | SVI H2 | SVI H3 | NW H1 | NW H2 | NW H3 |
|---|---|---|---|---|---|---|
| headline (counted once, P ≥ 0.75) | +0.46, 0.073 | −1.11, 0.001 | −0.32, 0.002 | +0.32, 0.11 | +0.14, 0.69 | +0.16, 0.92 |
| fixed stations (SVI S13; NW fixed-9 boreholes, raw) | +0.47, 0.056 | −0.87, 0.008 | −0.18, 0.55 | +0.07, 0.14 | +0.17, 0.70 | +0.18, 0.73 |
| raw (no de-dup) | +0.47, 0.085 | −1.16, 0.0005 | −0.31, 0.002 | +0.24, 0.19 | +0.13, 0.70 | +0.15, 0.93 |
| P ≥ 0.5 | +0.29, 0.18 | −1.07, 0.0015 | −0.22, 0.006 | +0.15, 0.17 | +0.08, 0.61 | +0.14, 0.95 |
| 10-km bins (SVI) / bins shifted 5 km (NW) | +0.52, 0.080 | −1.08, 0.003 | −0.29, 0.004 | +0.63, 0.010 | +0.13, 0.60 | +0.19, 0.90 |
| more episodes (MAJOR ≥ 100 LFEs: 15 / 25) | +0.28, 0.062 | −0.72, 0.002 | −0.12, 0.25 | +0.38, 0.077 | +0.17, 0.73 | +0.10, 0.74 |
| Pearson residuals | p 0.058 | p 0.001 | p 0.004 | p 0.034 | p 0.44 | p 0.76 |
| activity ≥ 1 event | +0.17, 0.30 | — | — | −0.01, 0.49 | — | — |
| activity ≥ 5 and ≥ 0.5 × expected | +0.49, 0.078 | — | — | **+0.71, 0.003** | — | — |
| longer history (SVI + Bostock 2003–13; NW from 2012 / 2007) | +0.46, 0.087 | same | same | +0.33, 0.11 / 0.10 | same | same |
| NW without targets < 90 d after another | — | — | — | +0.20, 0.20 | −0.05, 0.34 | −0.02, 0.09 |

So H1 reaches p < 0.01 only in NW under two of ten variants (strict threshold, shifted bins). It is never there in SVI.

**Bostock 2003–2013 (SVI groups, his 10 ETS).** His episodes reach nearly every patch, so time since last active
hardly varies. H1 +0.13, p 0.08. H2 +0.05, p 0.73. H3 p 0.13. De-duplicated at ±2 s (drops 22 % of events): the same.
The turn-taking seen in SVI 2015–26 is not in his decade.
Joint SVI 2003–2026 sequence (20 big ETS, big-only history): H1 +1.19, p 0.043. H2 −0.77, p 0.001; this comes from the
2015–26 half. H3 −0.25, p 0.020.

**Scrambled control (stations shifted to different days, SVI).**
- At P ≥ 0.75 there are 0–9 events per big ETS, too few to make a pattern.
- At P ≥ 0.5 there are 35–150 per episode. The scrambled control gives H1 +0.25 (p 0.04), H2 −0.14 (p 0.03),
  H3 −0.04 (p 0.54).
- The data at the same cut give H2 −1.07, 8× stronger.
- The scrambled-control events also lean toward the tremor-covered patches (+0.41, p 0.0005, vs +1.34 in data). Either a
  few real LFEs leak through the scramble, or tremor-time noise raises false detections near active patches. The
  headline cut is not affected.

## Figures

- `svi_science/figures/q1_resid_vs_tsince.png`: residual share vs days since the patch was last active (SVI, NW, SVI
  2003–2026). Open orange circles are patches the tremor front did not reach. They make up the empty cloud at the
  bottom, at every time since last active.
- `q1_example_episodes_svi.png`, `q1_example_episodes_nw.png`: every big ETS as a map. Colour = residual share,
  number = days since last active, × = outside the tremor extent, title = front direction. Examples:
  - SVI 2019-08 (front from the Strait toward NW, stopped mid-island: the NW groups are empty although 400 days had
    passed);
  - SVI 2019-11 / 2021-01 / 2022-03 (NW half only: the SE groups are empty; they were last active 75–146 days
    before);
  - NW 2017-03 and 2020-10 (only the NW end).
- `q1_partial_r2.png`: within-episode R² of each predictor, and episode unevenness vs tremor extent.

## What this does and does not show

- **Shows.** Per-patch strength in a big ETS is set by the front's extent, plus (in SVI) which end it entered from. A
  simple per-patch slip budget (time since last slip, deficit) adds nothing beyond that in either region.
- **Does not show.** Why the front stops where it does. With 10–14 episodes per region, a budget effect on the front's
  reach smaller than about R² 0.1 cannot be excluded (L2 p 0.033 in SVI is a hint only).
- The SVI turn-taking of 2015–26 is real in our catalogue (fixed stations, calibration, per-patch control) but is
  missing in Washington and in Bostock's 2003–13. It is best read as a feature of how the recent SVI fronts behaved
  (alternating ends, half-island episodes), not as a persistent law.

## Files

`counts_long.parquet` (patch × episode counts, all variants), `patches.csv`, `episodes.csv`, `tremor_ep.csv`,
`det_patch_ep.csv`, `rows_<design>.csv` (per patch × big ETS: residual, days since last active, censored flag,
previous / cumulative residual), `kin_<design>.npz` (front predictors), `tests.json`, `tests_table.csv`,
`partial_r2.csv`, `episode_level_values.json`, `scrambled_control.json`, logs `q1_build.log`, `q1_tests.log`,
`q1_scrambled.log`.
