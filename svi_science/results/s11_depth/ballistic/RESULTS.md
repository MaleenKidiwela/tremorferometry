# S11 ballistic test: results (Mon 2026-10-12, 3:25 pm PT)

Rules: `PLAN.md` (written 2:45 pm PT, before any lag; Amendment 1 at 2:50 pm after a family-10-only test run; Amendment 2
at 2:58 pm, written while the full run was in flight). Scripts: `svi_science/scripts/s11_ball_measure.py`,
`s11_ball_analysis.py`, `s11_ball_posthoc.py` (post hoc, labelled). Tables are in `tables/`, figures in `figures/`.
"Scrambled control" means the S10 c2 windows: each station's data shifted to another time of the same day.

## What the data show

**There is no episode-to-episode path signal in the ballistic arrival times above what this data set can detect.
The co-located floor is not 1–2 ms; at 76-family scale it is about 3.5–4 ms (robust) or 10 ms (plain SD).**

1. **Floor (co-located pairs).** The pre-registered floor is the dof-corrected SD of the S difference between instruments at the
   same site (PGC vs B009/B010/B011, and B011 vs B009/B010). It is **10.2 ms with PGC as reference and 10.0 ms with B011 as
   reference**, against the rule of ≤ 2 ms. So timing precision is **not confirmed** at this scale. The robust (MAD)
   scatter is 4.1 / 3.5 ms. The heavy tail comes from poorly measured stacks:
   - for values whose bootstrap σ is ≤ 2 ms (435 + 424 values, about 30 %), the floor is **3.0 / 2.1 ms SD and 1.9 / 1.7 ms robust**;
   - the scatter grows with σ, bin by bin (`tables/posthoc_scatter_by_sigma.csv`), at about 1.3–2 × σ;
   - so Merlin's ~1.2 ms from family 10 at PGC E1 holds for the best stacks only (family 10 in this run: PGC E1–B011 SD 1.1 ms).

2. **The co-located pairs carry an episode-common term**, which cannot be path because the instruments share a site:
   SD 1.9 ms (PGC reference) and 2.2 ms (B011 reference), both with permutation p = 0.0001. These are **station-timing
   glitches** that hit every family in an episode (`tables/posthoc_floor_top_episodes.csv`):
   - **PGC**: +6 to +10 ms in 2019-11 and −4 to −10 ms in 2021-01 against every other station;
   - **B010**: −24 ms in the 2009-03 burst, −15 ms in 2024-08 and −10 ms in 2022-03 against B011.
   These glitches are exactly the kind of thing that would pass for a "path change" with a single reference station.

3. **Far pairs** (B926, SNB, VGZ; S cc 0.77–0.84 against 0.93–0.94 at the PGC site):
   - the scatter is 18–21 ms;
   - the median bootstrap σ per value is 9–14 ms, and the bootstrap noise alone accounts for most of the scatter (noise SD 15–22 ms);
   - episode-common SD is 0–4.6 ms, permutation p 0.006–0.69;
   - family-specific SD is 0–7 ms.

   | pair | scatter (ms) | ratio to floor | episode-common SD | family-specific SD | p | rule verdict |
   |---|---|---|---|---|---|---|
   | PGC–B926 | 20.6 | 2.03 | 4.0 | 0.0 | 0.36 | WANDER (borderline) |
   | PGC–SNB | 20.7 | 2.03 | 4.6 | 7.2 | 0.08 | WANDER (borderline) |
   | PGC–VGZ | 18.4 | 1.81 | 2.6 | 5.4 | 0.006 | FLAT |
   | B011–B926 | 20.0 | 2.00 | 2.0 | 0.0 | 0.69 | FLAT |
   | B011–SNB | 17.9 | 1.78 | 3.3 | 0.0 | 0.20 | FLAT |
   | B011–VGZ | 19.5 | 1.94 | 0.0 | 4.5 | 0.24 | FLAT |

   - The two WANDER calls sit at ratio 2.03 against the cut of 2, and B926's family-specific part is zero. They are
     noise-level calls, not evidence of locus wander.
   - The only p < 0.01 (PGC–VGZ) does not reproduce with B011 as reference (p 0.24). Its largest episode terms are the
     PGC glitches of 2019-11 and 2021-01. Under Amendment 1 it is not PATH.
   - Majors only: every far pair is FLAT (p 0.006–0.52).
   - Overall pre-registered call: PGC reference, 2 of 3 WANDER (borderline); B011 reference, 3 of 3 FLAT. No PATH anywhere.

4. **S−P is limited by P, as Merlin expected.**
   - P cc medians: 0.24–0.42 on Z plus the polarised horizontal, 0.27–0.55 on Z only. Adding the polarised horizontal
     *lowered* cc everywhere except B009.
   - Bootstrap σ(S−P) is 18–45 ms per episode stack. The S−P scatter of 26–38 ms equals the co-located d(S−P) floor
     (31 ms), so every station is FLAT.
   - One station has an episode-common S−P term at p < 0.01: B010, SD 7 ms, p 0.005. It is absent at co-located B009
     and B011 (p 0.65 and 0.63), so it is a B010 matter, not a path. PGC: p 0.035.
   - The ±100 ms search did stop the cycle skips: P railed in 4–17 % of stacks; no ±300–400 ms jumps as in the family-10 pilot.

5. **Self-checks.** A +5 ms Fourier shift reads back as 4.94 ms (median; 5–95 % 4.73–5.28). S railed in 1.3 % of the data
   stacks against 40 % of the scrambled-control stacks. Median S cc is 0.90 for the data against 0.19 for the control.

6. **Scrambled control.** Few control stacks survive the n ≥ 20 rule (1,577 station rows, 2,264 pair rows). Of those,
   40 % rail; the rest scatter 50–90 ms with no episode-common part at p < 0.01 (p 0.06–0.78). The decomposition does not
   manufacture an episode-common term from noise.

## What this means in δβ/β (post-hoc sensitivity, `tables/posthoc_sensitivity.csv`)
Method: inject a random episode-common term of SD s into the real demeaned values and re-run the 1,000-draw permutation.
The table gives the smallest s caught at p < 0.01 in ≥ 75 % of trials.
- **Inter-station S, far pairs:** s ≈ 6–8 ms (PGC–VGZ 2–4 ms; B011–B926 8 ms). With T_S ≈ 15–20 s, that is
  - **≈ 0.03–0.05 % if the change sits on the far station's leg only**;
  - ≈ 0.12–0.16 % if it is uniform on both legs, since only the 4–6 s travel-time difference then counts. The PGC–VGZ
    pair, at ΔT_S 0.4 s, has no uniform-change sensitivity.
- **S−P:** s ≈ 6 ms at PGC (0.10 % at constant Vp/Vs; 0.04 % if only Vs changes), 8–12 ms at B011, and
  12–16 ms at B926/SNB/VGZ (0.15–0.26 % constant Vp/Vs; 0.06–0.11 % Vs only).
- Measured episode-common SDs on the far pairs (0–4.6 ms) are below these detection levels. The data are consistent with
  zero path change. The bound is roughly "no episode-common S-time change larger than ~8 ms (~0.04 % on a 20 s leg)".

For comparison, S10/S11 coda episode σ is ≈ 0.09 % per episode. Along a far station's leg, the ballistic S difference is
2–3 × tighter per episode as a bound, but on a straight ray rather than on the deep coda kernel.

## Verdict against the pre-registered rules
- **Floor ≤ 2 ms: NOT met** (10 ms SD; 3.5–4.1 ms robust). It is met only by the best-measured ~30 % of values (σ ≤ 2 ms: 1.7–1.9 ms robust).
- **PATH: not found** with either reference, all units or majors only.
- **WANDER: not demonstrated.** The two rule-level WANDER calls are at the 2 × floor edge; the family-specific SD is 0–7 ms
  and within the misfit of the bootstrap σ (co-located χ²/noise ≈ 3).
- **Closest to FLAT:** the far-pair and S−P scatter is set by measurement noise at the far stations (S) and by P everywhere (S−P).

## What would move it
These are honest options, not a recommendation.
- **S at the far stations is the limit** (σ 9–14 ms against 3–5 ms at the PGC site). More rows per stack, or weighting by
  1/σ², would tighten it. With equal noise, the PGC-site level is ~3 ms.
- **P needs a different estimator** (Z only did better than Z plus polarised; per-row realignment of P is not possible
  without the S-anchored drift).
- **Station-timing glitches** (PGC 2019-11/2021-01; B010 2009-03, 2022-03, 2024-08) must be modelled as station-episode
  terms in any ballistic inversion. Otherwise they become fake path signals.
- **Locus wander:** the 3-parameter locus nuisance Merlin describes is not needed at this noise level. Nothing family-specific is resolved above noise.

## Files
- `tables/lags_station.parquet`: per (family, station, epoch, episode, data/ctrl): S, P (Z + polarised), P(Z) lags, cc,
  railed flags, bootstrap σ, inj5S.
- `tables/lags_pair.parquet`: per pair and episode on common rows: dS, dP, σ.
- `tables/series_geom.csv`: S/P picks and polarisation.
- `tables/decomposition.csv`: every item, subset (all/majors), data/ctrl: scatter, MSB/MSW, variance components, F, p,
  δβ/β conversions, verdict.
- `tables/episode_common.csv`, `tables/station_quality.csv`, `tables/summary.json`, and `tables/posthoc_*.csv`.
- Figures:
  - `figures/fig_dS_series.png`: per pair, dS per episode, cc-coloured, one grey line per family, black = episode-common;
  - `figures/fig_SmP_series.png`: the same for S−P;
  - `figures/fig_variance_decomp.png`: variance decomposition bars, data and scrambled control.
- `*_test.*` = the family-10 test run (Amendment 1 diagnosis).
