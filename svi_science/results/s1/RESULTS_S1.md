# S1 results — episode size from GNSS (Tier A) and slip-model gate (Tier B)

Run 2026-10-07 ~12:10 pm PT. Rules frozen beforehand in `PLAN_S1.md`. Script `svi_science/scripts/s1_geodetic.py`
(CPU only, ~1 min, deterministic seed). Log `s1.log`.

## Plain answer
- **Can GNSS measure the size of SVI episodes? Only for the big ones.** 10 of the 23 frozen-major episodes have a
  GNSS offset above 3σ of the random-window control. That is between the fail bar (8) and the accept bar (15), so
  **Tier A = PARTIAL**. The count is 10 with every projection tried (own ETS direction, fixed 236°, Slab2 231°), so
  the circularity of the own-direction choice does not matter. None of the 40 minor episodes is detected.
  - The 10 detections: 2014-11, 2015-12, 2017-03, 2018-06, 2019-08, 2021-01, 2021-10, 2023-04, 2024-08, 2025-09b1.
  - Eight of these are big ETS with a 3.8–6.0 mm mean offset (S/N 4–22).
  - The other two, **2019-08 and 2021-10, are signal at the southern (Olympic-side) stations only**. 2021-10 also
    carries a strong Washington signal (WA S/N 6.0) while the northern SVI stations see nothing (S/N −0.4), so its
    "SVI slip" is probably Washington slip. That leaves **8 clean SVI-wide sizes, plus 2 flagged**.
- **Many LFE "major" episodes have no geodetic size at all.** 2017-10, 2017-12, 2023-02, 2025-08 and 2026-04 have
  430–715 family-exposure days each, but their GNSS offsets are 0.0 ± 0.3–0.8 mm. Their 3σ upper bounds are about
  1–2 mm, against ~5 mm for a big ETS.
- **Can the slip model be used? Yes for whole-SVI size, no for north/south.**
  - **Tier B regional: PASS, at the edge.** At ALBH, PGC5 and SC04 the model reproduces the observed offsets of the
    10 largest episodes with slope 0.71 (95 % episode-bootstrap 0.65–0.74, r = 0.90; E/N components 0.69). The bar
    is 0.7–1.3, so the model under-predicts the observed offsets by ~30 %.
  - Model regional slip ranks the episodes like Tier A does: Spearman 0.90 over the 23 frozen majors, 0.78 over all 79.
  - **Tier B along-strike: FAIL.** The resolution trace over the 37 gate patches is 2.05, which passes (≥ 2), but it
    is only 1.64 over the 29 family patches. Neither spike test passes: north-half spikes are recovered at 0.36 and
    south-half spikes at 0.48 (both below 0.5). Cross-talk is 0.29 / 0.28 (below 0.3, so that part passes).
  - Consequence: the model's north/south shares are not used (`tierB_north_share_gated` is empty). This agrees with
    RESOLUTION.md (~1 independent blob).

## Tier A details (`episode_geodetic.csv`, `tierA_*.csv/json`)
- **Stations.** SVI-14 = Merlin's 13 stations + SC03. Each episode uses a median of 13 stations (11–14).
- **Own ETS azimuths** are 227–263°; the median is about 246° (WSW).
- **Random-window control.** Each episode gets 37–444 season-matched windows (median 236); no episode needed a
  widened match.
  - σ per episode is 0.26–0.93 mm (median 0.47).
  - The drift bias (inter-ETS loading, opposite to the ETS direction) has a median of −0.15 mm and reaches −0.69 mm
    for the 31-day 2014-11 window.
  - The control's 99th percentile sits at a median of 2.1σ, so the 3σ bar is conservative.
- **Station bootstrap.** For frozen majors σ_boot has a median of 0.38 mm. It measures station-to-station spread,
  not detection noise.

| detected (S/N > 3) | D (mm) | σ | S/N | S/N north | S/N south | S/N WA | model slip (mm) | model Mw SVI / N-Cascadia |
|---|---|---|---|---|---|---|---|---|
| 2014-11 | 4.33 | 0.77 | 5.7 | 2.9 | 7.5 | 4.3 | 8.1 | 6.26 / 6.69 |
| 2015-12 | 5.95 | 0.85 | 7.0 | 6.6 | 6.4 | 3.0 | 13.7 | 6.42 / 6.75 |
| 2017-03 | 5.59 | 0.72 | 7.7 | 7.1 | 7.3 | 4.2 | 12.7 | 6.40 / 6.76 |
| 2018-06 | 5.51 | 0.26 | 21.6 | 17.6 | 17.0 | 1.7 | 14.4 | 6.43 / 6.74 |
| 2019-08 (south only) | 1.59 | 0.40 | 4.0 | 0.4 | 7.0 | 2.5 | 3.4 | 6.01 / 6.11 |
| 2021-01 | 3.77 | 0.93 | 4.0 | 4.7 | 3.0 | 1.4 | 9.9 | 6.32 / 6.71 |
| 2021-10 (south only, WA leak) | 1.91 | 0.56 | 3.4 | −0.4 | 6.8 | 6.0 | 3.4 | 6.02 / 6.36 |
| 2023-04 | 5.05 | 0.59 | 8.6 | 6.7 | 8.0 | 9.6 | 11.2 | 6.36 / 6.72 |
| 2024-08 | 5.08 | 0.61 | 8.3 | 4.9 | 12.7 | 12.8 | 12.1 | 6.38 / 6.64 |
| 2025-09b1 | 5.13 | 0.51 | 10.0 | 6.7 | 11.4 | 9.8 | 11.1 | 6.36 / 6.66 |

**Near misses** (frozen majors with 2 < S/N ≤ 3):
- 2019-11 (2.81) and 2022-03 (1.98) are large LFE episodes, with 1,730 and 1,566 family-exposure days. Their
  signal sits on the **northern** stations only (north S/N 4.7 and 3.8).
- 2020-11 (2.91) shows signal at the **southern** stations only (south S/N 4.8).

These one-sided cases are a post-hoc observation, not a test. They fit along-strike-limited slip, but the GNSS
cannot confirm that: the Tier B along-strike gate failed.

## Leakage check (SVI vs Washington-NW; `leakage_svi_vs_wa.csv`)
- 61 of 79 SVI episodes overlap a NW episode in time.
- The WA set (12 Puget/Olympic stations) detects 13 of 34 NW frozen majors from 2014 on. Normalisation sizes:
  WA 3.1 mm, SVI 5.1 mm.
- `wa_leak_possible` is raised for 6 episodes: 2017-02, 2018-05, 2019-05, 2020-10, 2021-10 and 2022-05. In all of
  them WA S/N ≥ 3 and the normalised WA size is at least twice the SVI size.
  - Four of these are SVI majors: 2018-05, 2019-05, 2020-10 and 2021-10.
  - In each of the four, the SVI-14 offset is small or southern-only.
- **Caution on the SVI-14 set itself:** 5 of the 14 stations (P435, P436, PTAA, COUP, SC03) sit within 52 km of the
  WA-NW family centroid. Washington slip therefore enters the "SVI" mean through the southern subset. The north-subset
  column (PGC5, SC04, NEAH, PTRF, NANO, CLRS) is the cleaner SVI measure.

## Tier B details (`tierB_gate.json`, `tierB_model_vs_obs.csv`, `tierB_spike_maps.npz`)
- **B1. Resolution trace (thrust, IC-weighted R).**
  - 2.05 over the 37 patches within 15 km of a family.
  - 1.64 over the 29 family patches.
  - 4.30 counting strike-slip as well.
- **B2. Spike tests.**

  | input half | recovery | cross-talk | total potency recovered |
  |---|---|---|---|
  | north | 0.36 | 0.29 | 0.93 |
  | south | 0.48 | 0.28 | 0.95 |
  | up-dip | 0.43 | 0.31 | — |
  | down-dip | 0.38 | 0.33 | — |

  The splits are at the median along-strike coordinate (+2.0 km in the S0 frame) and at the median depth (39.7 km).
- **B3. Modelled vs observed offsets.** Slope 0.709 over 30 points (per station: ALBH 0.79, PGC5 0.68, SC04 0.68).
  The model's offsets are systematically ~30 % small. The cause could be the 13-IC truncation, the smoothing prior,
  or the missing seasonal/other ICs. Whatever the cause, absolute model slip/M0 is probably low by about that factor.
- **B4. Model Mw.**
  - Big episodes: 6.6–6.8 over 47–50.5 N, and 6.3–6.4 over the SVI patch set alone (μ = 40 GPa). Given B3, the true
    values are likely ~0.1 higher.
  - **No per-episode published comparison was made.** The Michel et al. (2019) table (2007–2017.6, overlapping only
    2014-06 … 2017-03 here) is not on disk, and the paper PDF could not be fetched (over the size limit).
    Gualandi (2025) is the model itself, so it is not an independent check.
  - The only check is that the range matches the commonly cited Mw ~6.5–6.8 for northern-Cascadia ETS and the repo
    README ("Mw 6.6–6.8, 1–2 cm"). This is a consistency check, not a validation.
- **B5. IC 26 share of SVI regional slip, 2024-01 … 2026-08.**
  - Covariance share: IC26 3.4 %, IC17 3.5 %. IC10 (66 %) and IC11 (18 %) dominate.
  - Share of summed |increments| over the 2024–2026 episodes: IC26 5.7 %. Its largest single increment is +1.0 mm
    (2025-07); the 2024-08 and 2025-09b1 totals are 14.7 and 13.2 mm.
  - IC26 does not drive the 2024–2026 SVI episode sizes.

## What this means for S2
- Under the frozen rule (8–14 detected = PARTIAL), the GNSS side of S2 uses only the detected episodes.
- That gives 8 clean episodes plus 2 flagged ones. All of them are big, so log N vs log D spans less than a factor of
  ~4 in D. That is a weak lever for the slope b.
- Two options for the user / Merlin:
  - (a) Use Tier B regional slip as the size for all episodes. It passed the regional gate and ranks like Tier A
    (ρ 0.90), but increments below ~3 mm are unresolved (RESOLUTION.md §4).
  - (b) Use the Tier A 3σ upper bounds on the null-offset LFE majors as censored points. This is the most informative
    new fact here: several LFE-rich episodes carry < ~1.5 mm of SVI-wide offset.
- Per-family or north/south slip shares from the model stay forbidden.

## Files
- `episode_geodetic.csv`: one row per episode (79). It holds D, σ, S/N in three projections, subsets, WA leakage,
  Tier B columns and flags.
- `tierA_station_offsets.csv`, `tierA_control.csv`, `tierA_summary.json`, `leakage_svi_vs_wa.csv`
- `tierB_gate.json`, `tierB_model_vs_obs.csv`, `tierB_spike_maps.npz`
- Figures: `svi_science/figures/s1_episode_size.png`, `s1_model_vs_obs.png`, `s1_spike_maps.png`

**Correction 2026-10-07 (~12:30 pm PT, Merlin review):** the "five LFE-rich majors with ~0 mm offset" (2017-10, 2017-12, 2023-02, 2025-08, 2026-04) is a labelling error, not a result. "LFE-rich" used `family_exposure_days` (tremor exposure, saturates at 115 families/day), not LFE counts. De-duplicated P ≥ 0.75 counts for these episodes are 78, 139, 0, 67, 31 (vs 880–1,498 in the clean detected ETS); linear scaling predicts 0.1–0.6 mm (LFE) or 0.5–1.0 mm (tremor), inside every 3σ bound. Four of the five are northern deep-edge episodes (group 10). 2023-02 (719 epicentres, 0 LFEs, n_sta_day 10) is a probable detection-problem episode. Also: the "D" column in episode_rates/regional_svi_eff.csv is exposure days, not GNSS offset. Merlin's S2 design: svi_science/results/s2/MERLIN_S2_DESIGN.md.
