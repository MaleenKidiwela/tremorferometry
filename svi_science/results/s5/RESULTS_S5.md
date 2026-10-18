# S5 results — tidal sensitivity of the SVI LFE families, 2003–2026

Run 2026-10-07, about 2:05–3:12 pm PT. The rules were frozen in `PLAN_S5.md` before any modulation number was computed.
Two post-hoc additions are labelled there (appendices A and B). CPU only, < 5 GB per process.
No Merlin review (run directly at the caller's request), so treat everything as provisional.

## Plain answer
- **Our tide model reproduces Royer et al. (2015). Step 2 PASSES.**
  - Data: Bostock's 269,423 LFEs (2003–2013, 130 families).
  - More LFEs occur when the tidal shear stress pushes the fault toward thrust slip.
  - The rate swings by about ±35 % over an M2 tidal cycle. The scrambled control (0.3 % on 406k candidates) and random
    times show no swing.
  - On the first day of activity at a family the tidal response is about zero. It rises to a maximum about 3–4 days later.
  - Early in an episode, LFEs peak when the shear stress is rising fastest. Later they peak at the shear-stress maximum.
  - All three points match Royer's abstract.
  - So the tide model, including its phases, works. The README had flagged the phases as never checked; this result
    checks them indirectly.
- **Our own 2014–2026 catalogue shows the same response.**
  - The size, the phase and the within-episode rise all match (10,378 LFEs, counted once, P ≥ 0.75).
  - The rise holds whether the clock starts at each family's own onset or at the arrival of the PNSN tremor front.
  - Measured from the front, the response peaks later, 4–10 days after it. Activity at a family begins 1–3 days after
    the tremor front arrives, and the strong response starts after that.
- **Has tidal sensitivity changed over 23 years, or between eras? Not detectably.**
  - Some tests do give small p-values: 2014–19 vs 2020–26, 2003–08 vs 2009–13, and shallow families 2014–19 vs 2020–26.
    Those tests use only the tidal-sampling noise.
  - But single episodes differ from each other far more than tidal-sampling noise allows (χ² 99 on 9 df in Bostock's
    catalogue; 46 on 14 df in ours).
  - Tested against that episode-to-episode scatter, no era difference is significant (p 0.18–0.60).
  - What we could have seen: a change in modulation of about ±25–30 percentage points (out of ~35 %). That is roughly
    a factor of 1.8.
- **Within episodes:** the rise after onset is the same in every era and both halves of ours (all p ≤ 0.004).
  - Its shape correlates at 0.86 between Bostock's catalogue and ours.
- **Depth:** there is no depth effect in ours.
  - In Bostock's catalogue the deepest third responds about half as strongly (a 20 % vs 39 % raw).
  - After allowing for when in the episode the events fall, that becomes 20 vs 29 % (p = 0.015). This misses the
    0.01 bar. It is a hint, not a result.

## Royer et al. (2015): what we compared against
The full text was not accessible. The abstract was retrieved verbatim (Crossref, doi:10.1002/2014JB011430):

> "...Sensitivities of LFE families to tidally induced fault normal stress, updip shear stress (UDSS), and corresponding
> time derivatives are computed... We observe that sensitivity of LFEs to UDSS rises dramatically from near zero on the
> first day of strong activity to a maximum ∼4 days later. In addition, the peak LFE rate transitions from a correlation
> with peak tidal shear stress rate to a correlation with peak tidal shear stress through large slow slip events. We
> identify 64 Rapid-Tremor-Reversals (RTRs)... a majority (72%) occurs during periods of large positive UDSS."

A secondary source (Zhou et al. 2026, arXiv 2602.06703, Table 4 and text) adds three points: a tidal shear stress of
about 2 kPa, aσ ≈ 7.6 kPa (Beeler 2013 model), and a phase of about −90° early moving to ~0° late.
"Modulation of tens of %" is from memory only and is consistent with exp(2/7.6) − 1 ≈ 30 %.
The RTR result was not attempted.

## Step 1 — tidal stress 2003–2013 (PASS, exact)
- Stress was computed hourly from 2003-01-01 to 2014-01-31 UTC at all 130 Bostock hypocentres, using the S5-setup model unchanged.
  - Body tide: pyTMD IERS2010.
  - Ocean load: EOT20 plus the inland tide-gauge field, half-space, out to 700 km.
  - Fault: Slab2 plane, slip direction 236°.
  - Output: `tidal_stress_bostock_hourly.parquet` (130 × 97,176 h).
- Acceptance checks (`step1_acceptance.json`):
  - For the 115 shared families, the ocean constituent stresses equal the setup's to 2e-15 (relative).
  - The hourly τ, σn, body and ocean parts in 2014-01 equal the setup's parquet exactly (max |diff| 0 Pa).
  - The median τ std is 1.22 kPa (setup: 1.22 kPa).
- Event→family key: exact (lat, lon, depth) = e5r bkey, the same key P3 used. All 269,423 events are mapped; none fall outside the 10 windows.

## Method (as pre-registered)
- **M2 response:** θ = M2 phase of each family's shear-stress series, with θ = 0 at peak thrust-promoting shear and
  θ = −90° at the peak stress *rate*.
  - The response is excess C = mean cos θ of the LFEs minus the same for surrogates.
  - The modulation amplitude a is 2 × |excess vector|, so that rate ∝ 1 + a cos(θ − θ0). θ0 is the preferred phase.
- **Stress level:** Δτ = mean total τ at the LFE times minus the surrogate mean. The implied exponential sensitivity is
  aσ ≈ Var(τ)/Δτ.
- **Surrogates (1,000):** each family's active days are shuffled among themselves. Every LFE moves with its whole day
  and keeps its time of day.
  - This preserves within-day bursts and any time-of-day detection effect.
  - It randomises M2 (24°/day drift) and O1 (25°/day drift).
  - K1 and S2 are locked to the solar day, so they are reported but not tested. They are confounded with day/night
    detection.
- **Self-tests** (`selftest.json`):
  - The phase convention correlates at 0.9996 with the M2 part of τ.
  - The shifted-phase formula is exact to 0.007°.
  - An injected 30 % / 0° modulation was recovered as 34 % / +7°. An injected −90° lag was recovered as −102°.

## Step 2 — positive control on Bostock 2003–2013 (130 families): **PASS**
| quantity | value | control |
|---|---|---|
| M2 modulation a | **35.0 %** (surrogate SD 1.4 %) | surrogates: mean |R| 0.002 |
| preferred phase θ0 | **−32°** (between peak stress rate and peak stress) | — |
| excess C (P1) | 0.148, **p = 0.001** (one-sided; 0/1000 surrogates reach it) | — |
| Δτ (P1) | **+250 Pa**, p = 0.001 → aσ ≈ **5.1 kPa** (Royer via Zhou: ~7.6 kPa) | — |
| up-dip shear (Royer's UDSS) | +245 Pa, p = 0.001, aσ 5.2 kPa | — |
| ocean vs body part of Δτ | +238 Pa vs +12 Pa (the ocean load carries the signal) | — |
| normal stress σn (tension +) | −173 Pa: LFEs slightly favour *clamping*, not unclamping | Coulomb μ′ = 0.4 still positive (+181 Pa) because τ dominates |
| O1 (lunar diurnal) | a 5.3 % ± 1.5 %, p = 0.10, not detected | expected ~8 % if the same aσ applies (O1 τ is 0.23 × M2) |
| K1 | a 6.8 %, θ0 −22° | not testable (time-of-day locked) |
| random times on the same days | M2 Rayleigh R = 0.0046 (p 0.0035, i.e. a 0.9 % artefact; see PLAN app. A) | flat to < 1 % |
| 115 production families only | a 35.4 %, θ0 −34° | the 15 others: 29.7 %, −1° |

**Within-episode rise (P2): PASS.**
- Response on days 2–6 minus the first day: +0.454, p = 0.001 (day permutation; null mean 0.022, SD 0.017).
- With the first-event onset instead: +0.124, p = 0.001.
- With stress level instead of phase: p = 0.001.

| days since family onset | <1 | 1–2 | 2–3 | 3–4 | 4–6 | 6–10 | ≥10 |
|---|---|---|---|---|---|---|---|
| LFEs | 145,588 | 49,328 | 27,173 | 15,522 | 14,969 | 11,741 | 5,102 |
| a (%) | 13.5 | 67 | 92 | **94** | 86 | 28 | 36 |
| θ0 | **−90°** | −43° | −18° | −5° | +2° | +29° | −52° |
| excess C, p | 0.000, p 0.48 | 0.25, 0.001 | 0.44, 0.001 | 0.47, 0.001 | 0.43, 0.001 | 0.12, 0.001 | 0.11, 0.02 |
| aσ (kPa) | 32 | 3.3 | 1.9 | 1.7 | 1.8 | 5.7 | 8.4 |

- This is Royer's result: near zero on the first day, maximum ~3–4 days later, and the phase moving from the stress
  rate (−90°) to the stress (0°).
- Every one of the 10 ETS windows shows positive modulation (a 28–57 %, p = 0.001) except the sparse 2007-01 window
  (3,147 events, a 9 %, p 0.82).

## Step 3 — our 2014–2026 catalogue (115 families)
**Overall, with controls:**

| set | n | a (%) | θ0 | p | Δτ (Pa) |
|---|---|---|---|---|---|
| data, P ≥ 0.75, counted once (headline) | 10,378 | **36.8** (SD 2.4) | **−36°** | 0.001 | +261 (aσ 4.7 kPa) |
| data, P ≥ 0.5, once | 22,378 | 37.0 | −35° | 0.001 | +278 |
| data, P ≥ 0.85, once | 5,887 | 34.8 | −37° | 0.001 | +240 |
| data, P ≥ 0.75, raw | 13,466 | 37.5 | −36° | 0.001 | +266 |
| **scrambled control**, P ≥ 0.75, once | 65 | 4.2 | — | 0.65 | +38 (p 0.34) |
| **scrambled control**, P ≥ 0.5, once | 1,928 | 3.5 (SD 3.8) | — | two-sided 0.45 | −23 (p 0.75) |
| **scrambled control, all candidates** (detector check) | 406,419 | **0.3** (SD 0.8) | — | two-sided 0.89 | −9 (p 0.87) |

- Random times on the same days are flat: Rayleigh p = 0.57.
- O1: a 4.7 %, p 0.61, not detected.
- The two flagged episodes (2017-03, 2025-10, 212 LFEs) have a 45 %; the rest have 36.7 %.
- LFEs outside the major episodes (1,218) also respond: a 32 %, θ0 −9°, p 0.001.

**Within-episode rise (ours):**

| clock | first day → days 2–6 | p (day permutation) |
|---|---|---|
| days since PNSN tremor front (10 migrating episodes; headline per PLAN) | +0.228 | 0.001 |
| days since family's own onset, all major episodes | +0.457 | 0.001 |
| own onset, the same 10 episodes | +0.495 | 0.001 |
| front, 2014–2019 only / 2020–2026 only | +0.171 / +0.242 | 0.003 / 0.001 |

| days since own onset (majors) | <1 | 1–2 | 2–3 | 3–4 | 4–6 | 6–10 | ≥10 |
|---|---|---|---|---|---|---|---|
| LFEs | 6,121 | 1,533 | 586 | 350 | 316 | 171 | 83 |
| a (%) | 24 | 76 | 82 | **120** | 106 | 59 | 73 |
| θ0 | −69° | −42° | −24° | +10° | +22° | −24° | −15° |

| days since PNSN front | <1 | 1–2 | 2–3 | 3–4 | 4–6 | 6–10 | ≥10 |
|---|---|---|---|---|---|---|---|
| a (%) | 0.6 | 34 | 46 | 49 | **74** | 68 | 56 |
| θ0 | — | −79° | −59° | −40° | −15° | −5° | −50° |

- On both clocks the phase moves from about −80° toward 0°. This is the same transition Royer described, and the same
  as in Bostock's catalogue.
- Values above 100 % mean the response is no longer a simple cosine (rate close to zero at the opposite phase).
  They are not an error.

## Era comparisons (SHAPE across detectors; amplitude only within ours)
| comparison | surrogate-only test | standardised for days since onset (post-hoc) | **episode-level test (appendix B)** | could have detected (MDD, a) |
|---|---|---|---|---|
| ours 2020–26 vs 2014–19 | a 41.7 vs 32.3 %, p 0.037 | 38 vs 22 %, p 0.001 | **p 0.60** (15 episodes; active days only p 0.25) | 14 % (tidal noise) / **29 %** (episodes) / 53 % (active days) |
| Bostock 2009–13 vs 2003–08 | excess C 0.171 vs 0.133, p 0.004 | — | **p 0.25** (10 episodes; active days p 0.18) | 7 % / **25 %** / 34 % |
| Bostock vs ours (detectors differ) | a 35.0 vs 36.8 %; θ0 −32° vs −36° | 29 vs 38 %, p 0.008: **not interpretable as change** | — | — |
| shape of the rise, Bostock vs ours | normalised curves r = 0.86 over 7 bins; both peak at 3–4 d after onset; θ0 path −90° → ~+10…+30° in both | | | |

- Episodes really do differ from each other:
  - Bostock: χ² 99.5 on 9 df, between-episode SD of the excess 0.060, against a sampling SD of 0.019.
  - Ours: χ² 46 on 14 df, between-episode SD 0.077, against a sampling SD of 0.055.
  - Examples in ours: the two biggest 2014–19 episodes (2015-12, 2018-06) respond weakly (excess 0.02 / 0.07). 2014-11,
    2017-12 and 2022-03 respond strongly (0.30–0.35).
- So "is the cycle stationary?" can only be answered at episode level. The answer is "no change ≥ ~25–30 percentage
  points in a, i.e. about ×1.8, in either era".

## Depth
| comparison | raw | standardised for days since onset | MDD (a) |
|---|---|---|---|
| ours deep vs shallow third | 31.6 vs 49.1 % (excess C 0.157 vs 0.171), p 0.64 | 34.5 vs 28.6 %, p 0.41 | 17–19 % |
| Bostock deep vs shallow third | 20.5 vs 38.6 % (excess 0.082 vs 0.164), p 0.001 | 19.8 vs 29.0 %, **p 0.015** | 10 % |
| ours shallow third, 2020–26 vs 2014–19 | p 0.001 | 49.6 vs 14.7 %, p 0.001 | 26 % |

- The shallow-third era difference is post-hoc, uses tidal-sampling noise only, and was not re-tested at episode level.
  Given the episode scatter above, it is **not a detection**.
- Bostock's weaker deep response is consistent within his own catalogue. It is not reproduced in ours.
- 12 location groups (ours): a 23–57 %, all p ≤ 0.01 except G5 (p 0.04). Listed in `step3_table_data_075_once.csv`
  (label `group`).

## Caveats
- **Bursts.**
  - The surrogates keep within-day bursts together, but the p-values only measure tidal-sampling noise.
  - Episode-to-episode differences are larger. Episode-level tests are the honest ones for eras.
  - The Poisson bars in `s5_modulation_vs_stress.png` understate the noise.
- **The onset clock depends on detection.**
  - The 10th-percentile onset of each family-episode is set by the LFEs themselves.
  - The PNSN-front clock is detection-free, and it gives the same rise.
- **Model choices.**
  - The surrogates make the solar-day constituents (S2 exactly, K1 nearly) untestable.
  - The M2-phase statistic is unaffected by any day/night detection bias.
- **Normal stress.**
  - The σn result (LFEs slightly favour clamping) comes from a total-stress model. It has no pore-pressure (Skempton) response.
  - In this model σn and τ are not independent, so it is not separate evidence about friction.
- **Royer comparison.**
  - It rests on the abstract and a secondary table. Families, stress models and statistics differ, so only the sign,
    the order of magnitude and the time pattern are compared.

## Figures
- `svi_science/figures/s5_modulation_vs_stress.png`:
  - LFE rate (observed / expected) against tidal shear stress, in family-std bins: rising from ×0.64 to ×1.46 (Bostock).
  - LFE rate against M2 phase: ×1.38 at about −15°, ×0.65 at about +135°.
  - Both panels show Bostock 2003–13 and ours.
- `svi_science/figures/s5_sensitivity_vs_days.png`: amplitude and phase against days since onset or front (Bostock onset,
  ours front, ours onset).
- `svi_science/figures/s5_era_comparison.png`: a and θ0 for Bostock, 2014–19, 2020–26 and the scrambled control; the
  normalised rise shapes. Its error bars show tidal-sampling noise only.

## Files
- **Rules:** `PLAN_S5.md`.
- **Step 1:** `step1_acceptance.json`, `bostock_fault_geometry.csv`, `bostock_ocean_constituents.csv`,
  `tidal_stress_bostock_hourly.parquet`.
- **Step 2:** `step2_result.json`, `step2_modulation_table.csv`, `step2_extra_diffs.json`, `step2_hists.npz`,
  `step2_surrogates.pkl`.
- **Step 3:**
  - per set: `step3_table_<set>.csv`, `step3_result_<set>.json`, `step3_surrogates_<set>.pkl`;
  - `step3_front_arrivals.csv`, `step3_posthoc_standardised.json`, `step3_episode_heterogeneity.json`, `step3_per_episode.csv`.
- **Scripts** (`svi_science/scripts/`): `s5_bostock_stress.py`, `s5_lib.py`, `s5_selftest.py`, `s5_step2_bostock.py`,
  `s5_step3_ours.py`, `s5_step3_posthoc.py`, `s5_step3_episodes.py`, `s5_figures.py`.
