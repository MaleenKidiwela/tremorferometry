# S1 — Episode size from GNSS (Tier A) and slip-model gate (Tier B): frozen plan

Written 2026-10-07 ~12:10 pm PT, BEFORE any S1 number is computed. Design: Merlin
(`lfe_extension/results/slip_compare/MERLIN_DESIGN_2026-10-07.md` steps 1–2), resolution facts from
`lfe_extension/results/slip_compare/RESOLUTION.md`. Script: `svi_science/scripts/s1_geodetic.py`.
Anything changed after the first number is listed under "Deviations" at the bottom.

## Inputs (read-only)
- Episodes: `lfe_extension_wa/results/episode_rates/episodes_svi_final.csv` — 79 episodes (39 `major`, 40 minor);
  the **23 `in_frozen_major` episodes are the acceptance set** (Merlin: "≥ 15 of 23"). All 79 get a row.
  Pairs < 30 d apart (start of next − end of previous) are flagged `close_pair` (their 5-day windows may
  contain the neighbour's slip).
- GNSS: `/home/jovyan/cascadia-sse-strain/data/2026-08-18/X.mat` (`Xd.ts`, mm, already detrended by the
  trajectory model; `Xd.data_mask` = 1 for real days; masked days are filled in `ts` and are NEVER used).
  Dates: Gualandi timeline converted with `slip_build.fixed_dates` (floor, true UTC day), not the repo's
  rounding function (half the epochs one day late). The repo is imported read-only; nothing in it is modified.
- Model slip: `slip_cache.npz` (G, L, Ls, Ld, A, corrected dates) and `slip_svi_family_daily.parquet`
  (thrust-positive, corrected dates), both from `lfe_extension/scripts/slip_build.py`. The repo's fig07 sign
  (−thrust) is NOT used.

## Station sets (geometry only, fixed now)
- **SVI-14** = ALBH, PGC5, SC02, SC04, NEAH, PTRF, PTAA, NANO, P435, P436, P439, CLRS, COUP (Merlin's 13 named)
  + **SC03** (the closest station within 75 km of the family centroid not named; Merlin's "+1").
- North / south subsets: sign of the along-strike coordinate in the S0 frame (origin 48.4837 N 123.7485 W,
  strike 321.15°, `svi_science/results/s0/frame_groups.json`). Stations with along-strike > 0 = north.
- **WA set** (leakage check) = stations within 65 km of the WA-NW family centroid (47.85 N, 123.02 W) that are
  not in SVI-14: CHCM, KTBW, P426, ELSR, P437, CUSH, P419, PFLD, SSHO, SEAT, SMAI, P399.
- Note (not a result): 5 of SVI-14 (P435, P436, PTAA, COUP, SC03) sit within 52 km of the WA-NW family centroid,
  i.e. above Washington LFE families. This is why the leakage check matters.

## Tier A — model-free episode offsets
1. Per station, per horizontal component: offset = median(ts over the 5 days end+2 … end+6) −
   median(ts over the 5 days start−6 … start−2), masked days excluded; ≥ 3 valid days required on each side,
   else the station is missing for that episode.
2. Projection: each station's **own ETS direction** = unit vector of the mean horizontal offset over the 10
   largest frozen-major episodes, "largest" ranked by `family_exposure_days` (an LFE quantity, so GNSS does
   not choose the episodes). For those 10 episodes the direction is **leave-one-out** (mean over the other 9)
   to limit circularity. Second variant (column): **fixed azimuth 236°** for every station = the reverse of the
   Juan de Fuca–North America convergence direction near southern Vancouver Island (≈ N56°E; approximate
   value from the literature, not looked up here). Third variant (column): the reverse of the S0 Slab2 mean
   down-dip azimuth, 231.15°.
3. D_e = equal-weight mean of the projected offsets over available stations (mm, positive = ETS sense).
   σ_boot = std of D_e over 1,000 station bootstrap resamples. `n_sta` = stations used.
   Also D_e for the north and south subsets, and a template amplitude `amp_e` (least-squares scale of the
   episode's projected offsets onto the mean projected offset pattern of the 10 largest episodes; robust to
   missing stations), as columns.
4. Random-window control (no tremor or LFE input): for each episode, pseudo-episodes of the
   same span are started on every day 2014-01-15 … 2026-08-18 such that the whole measurement range
   [start−6, end+6] is ≥ 10 d from every SVI episode (all 79) and every NW `major` episode; season-matched =
   start day-of-year within ±60 d of the episode's start. If < 30 windows qualify, widen to ±120 d, then drop
   the season match (flag `ctrl_widened`). Each window gets the identical Tier-A computation (same directions,
   same stations, same missing-day rule).
   bias_e = median of the control D; σ_e = std of the control D.
5. Headline size: **D_corr = D_e − bias_e**; **S/N = D_corr / σ_e**. Detected = S/N > 3. The uncorrected
   D_e / σ_e is a column.
6. **Accept Tier A** if ≥ 15 of the 23 frozen-major episodes have S/N > 3 (own direction). **Fail** if < 8:
   then only the LFE detection curve is used downstream (S2 GNSS side dropped). 8–14: partial — only the
   detected episodes are used in S2, flagged. Fixed-direction counts reported alongside; if the own-direction
   and fixed-direction counts fall on different sides of a threshold, the fixed-direction result governs
   (conservative w.r.t. circularity).

## Tier B gate — can the slip model be used regionally?
Model = Gualandi (2025) vbICA solution rebuilt read-only (same as RESOLUTION.md), IC-weighted model resolution
matrix R̄ (IC weights = share of SVI slip variance, `slip_resolution.ic_weights`), thrust (dip-slip) block.
- Patch sets: (a) the 29 patches containing ≥ 1 SVI family; (b) **gate set** = patches whose centroid is
  ≤ 15 km (horizontal) from any SVI family (≈ Merlin's "36 patches under SVI").
- B1 **Resolution trace** = Σ diag(R̄_thrust) over set (b) (set (a) as column).
- B2 **Spike tests**: unit thrust on the north half of set (b) (split at the median along-strike coordinate of
  its patch centroids, S0 frame), recovered = R̄_thrust · m. Recovery = area-weighted mean recovered slip over
  the input half; cross-talk = area-weighted mean over the other half ÷ recovery. Same for south input, and for
  up-dip / down-dip halves (split at the median centroid depth; reported, not gated).
- B3 **Modelled vs observed offsets** at ALBH, PGC5, SC04 for the 10 largest episodes (same 10 as above):
  modelled displacement = Σ_k (G L_k) A_k(t) (retained ICs only), Tier-A offset rule applied to both.
  Slope = least squares through the origin of modelled on observed, using the own-direction projected
  offsets (30 points); E/N components (60 points) as a column. Pass: slope 0.7–1.3.
- B4 **Per-episode Mw**: model M0 = μ Σ area × thrust increment (signed; same 5-day rule on patch slip),
  μ = 40 GPa (as repo fig09), over (i) set (b) "SVI regional" and (ii) all patches 47–50.5 N "northern
  Cascadia". Compare with Michel et al. (2019) / Gualandi (2025) per-episode values if a table can be
  obtained; otherwise state exactly what was compared (e.g. the RESOLUTION.md 2010–2017 sanity values).
- B5 **IC 26 share**: share of the variance of the SVI regional slip (mean of r25 along-rake slip over
  families) carried by IC26 in 2024-01-01 … 2026-08-18, and its share of the summed |increments| of the
  2024–2026 episodes; IC17 (the late-2024 −88 excursion per RESOLUTION.md) reported alongside.
- **Tier B for whole-SVI (regional) use** requires B3 slope 0.7–1.3 (B4/B5 reported, flags only).
  **Tier B for along-strike (north/south) use** additionally requires trace ≥ 2 AND both N and S spikes
  recovered ≥ 50 % with cross-talk < 30 %. Only if this passes are the north/south model shares used.
- Per-episode Tier B columns (always computed, but labelled `tierB_usable` per the gate): regional mean
  thrust slip increment (patch set (b), area-weighted), regional M0 and Mw, family-mean r25 along-rake
  increment, north/south model M0 share.

## Leakage check (SVI vs Washington-NW)
- An SVI episode "overlaps" a NW episode (`episodes_nw_final.csv`, any kind) if their date ranges intersect.
- For overlapping SVI episodes compute Tier A at the WA set (own directions from the 10 largest NW
  frozen-major episodes 2014–2026 by `family_exposure_days`, LOO for those; same control rule).
- Normalised sizes: D_corr divided by the median D_corr of that region's frozen-major episodes with S/N > 3.
- Flag `wa_leak_possible` if WA S/N ≥ 3 and normalised WA size ≥ 2 × normalised SVI size.
  Flag `south_only` if SVI S/N ≥ 3 but the north-subset S/N < 2 while the south-subset S/N ≥ 3
  (signal carried only by the Olympic-side stations). Flags only — nothing is dropped.

## Not done here
No LFE counts enter S1 except `family_exposure_days` for choosing the 10 direction episodes. No per-family or
per-depth slip use. "Scrambled control" in this study = the random-window control above.

## Deviations
- None to the rules. Execution notes: the first launch crashed on a pandas attribute bug before writing any output
  (log kept in `s1_run1.log`); figure labels were tidied and the script re-run, numbers identical (fixed seed).
  Gate set (b) came out with 37 patches. The north/south station split by along-strike sign put NEAH in the north
  subset and SC02/P439 in the south subset (geometry, as frozen). Michel et al. (2019) per-episode table could not be
  obtained (not on disk; PDF fetch over size limit) -> B4 compared only against the commonly cited Mw range.
