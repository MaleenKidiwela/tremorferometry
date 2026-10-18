# SH_CHECKS: two checks on the Shelly (2025) Mendocino forward extension (2023-12 → 2026-09)

2026-10-01. CPU only. Everything was read from the cut windows (`results/cut/sh_ext/`, raw counts −15…+30 s around every ≥8×MAD candidate) and from saved results. No downloads, nothing deleted, no git. Scripts: `scripts/sh_checks_kmr.py` (extract + analyze), `scripts/sh_checks_kmr_env.py`, `scripts/sh_checks_kmr_health.py`, `scripts/sh_checks_pnsn.py`, `scripts/sh_checks_pnsn_daily.py`. Outputs and figures are in `results/sh_checks/`.

## Plain summary

**Check 1 (NC.KMR). The cc of 0.13 does not come from a changed waveform. After the swap, KMR records almost no usable ground motion.**
- Before the swap, KMR is healthy. A plain stack of the 100 strongest detections per family matches the old STS-2 template at cc **0.90** (median of 27 families). It peaks at lag 0. Templates built by the same recipe from two halves of the old epoch match each other at **0.91** and the old template at **0.90 / 0.88**. So the recipe works when there is signal.
- After 2025-09-17, stacks at the strongest detections show no LFE at KMR:
  - With detections made without KMR, the stack matches the old template at cc **0.18** (Horizon) and **0.17** (Compact Posthole).
  - With the highest-P(same) detections, it is **0.16 / 0.18**.
  - A scrambled control gives **0.17–0.19**. It takes the same windows but cuts each event's slice at its own random time.
  - The stack's signal-to-noise is **0.8–1.0**, against 3.4 before the swap.
  - The energy in the 3 s after the expected S is **1.02× / 0.98×** the pre-arrival energy, against 1.63× before. No energy peak appears anywhere in −15…+30 s.
- The data themselves show the cause (`kmr_spectra.png`, `kmr_health_monthly.csv`). KMPB, measured the same way, is normal.
  - **Trillium Horizon epoch (2025-09-17 → 2026-03-12):** every channel sits at a fixed offset of about ±8.6 million counts. 2–8 Hz noise is 80→190 nm/s and rises month by month, against 5–12 nm/s on the STS-2. The spectrum is a comb of electronic lines; the 4 Hz line is 40× the band median by Jan 2026. The ocean microseism (0.1–0.5 Hz) is 50–110 nm/s, against 200–900 nm/s normally.
  - **Trillium Compact Posthole epoch (from 2026-03-12):** E and N show about 0.5 counts RMS (19–33 distinct values), 1.6 nm/s white noise and 1–2 nm/s microseism, so essentially digitiser noise only. Z has a line comb (7 Hz line 46×), 75 nm/s at 2–8 Hz and 40 nm/s microseism.
  - **The old STS-2 failed before the swap.** Its HHE and HHN were stuck at a constant value from **2025-07-26 / 27**; only Z worked until the swap.
- **Checks (a)–(g) all point the same way:**
  - (a) Channel names and order are HHE/HHN/HHZ in every epoch, sorted the same way in the cache and in the windows.
  - (b) Polarity: after the swaps, the best signed cc per channel within ±2 s is at most 0.44 in size, with either sign about equally often. A polarity flip would show up as about −0.9; none does.
  - (c) Azimuths are 0°/90° in all StationXML epochs. Rotating the horizontals to the best angle lifts cc only from 0.16 to 0.21, the level of noise.
  - (d) Time offset: the lag scan over ±2 s is flat (median peak 0.15 at a random lag), and the energy stack rules out offsets up to −15/+30 s.
  - (e) Response: phase difference against the STS-2 at 2–8 Hz is a 2.4–4.6 ms delay (Horizon / Compact), under 0.1 sample. All sensitivities are positive; gain cancels in a normalised cc.
  - (f) n_stack was 100 for most families. But the new templates' cc50 (0.21–0.30 in the original run) is matched by the scrambled control built the same way (0.49 Horizon, 0.32 Compact in this re-run, the same as the real windows: 0.48 / 0.32). The new templates are stacks of instrument noise, and the electronic lines make random slices correlate with each other.
  - (g) The old-epoch control gives 0.87–0.91 (see the second bullet).
- **Verdict: drop KMR from 2025-09-17 on.** It cannot be fixed with these data; re-templating will not help. The no-KMR run (`mf_NK`, 16 stations) already covers every day from 2025-09-17 to 2026-09-28. Also consider Z-only KMR, or dropping KMR, for 2025-07-26 → 2025-09-17, when only its Z worked.
- **Effect on rates.** On the same days, the catalogue with the new KMR templates has **0.89×** (Horizon, 6472 vs 7289) and **0.85×** (Compact, 6902 vs 8139) as many S1 9× detections as the run without KMR. Per family the median is 0.91 / 0.85 (range 0.76–1.29). 95% of the with-KMR detections are also in the without-KMR run. So the current forward series after 2025-09-17 is probably **11–15% low** (S1 count), because the noise station dilutes the network mean cc. The same probably applies to ΣP(same), though this was not measured. A 16-station run also lacks the working KMR of before 2025-07, so it is not a perfect match to 2024 either. That gap cannot be measured without continuous data from the old epoch.

**Check 2 (PNSN tremor). PNSN tremor cannot confirm or refute our post-2024 rates. Within 30 km of the families it is too sparse and does not follow Shelly's own monthly counts.**
- Over 2018-01 → 2023-11, monthly PNSN tremor windows within 30 km and Shelly's monthly event count are **not related** (Pearson r = −0.08, Spearman −0.17). The same holds at 15 km (−0.06) and 60 km (−0.23).
  - PNSN has **zero** windows within 30 km from late 2022 to mid-2024, while his catalogue had 9,862 events in 2023.
  - His catalogue grew about 5× from 2018 to 2022 while PNSN did not.
  - PNSN's median station count per tremor window near Mendocino also changes: 7 (2018–21), 6 (2022–24), 5 (2025), 7 (2026).
- So the "fit his count from PNSN, then test ours" relation has no skill. The residual tests are given in `pnsn_relation_tests.csv`, but they mean little. For example, ours from 2024-12 to 2026-05 is 1.95× the PNSN-predicted count, while his own join-year count is already 1.77× it.
- **What PNSN does support is timing.** Days with PNSN tremor within 30 km have about 3× more LFEs than other days:
  - his 2018–2023: 2.6×; his join year: 3.5×;
  - ours, join year: 3.8×; ours, 2024-12 → 2026-05: 3.3×.
  - At 15 km the figures are the same. At 60 km they are 1.7–2.8×.
  - So the forward catalogue still lights up on the same tremor days as before. This is a check of timing, not of rate level.
- **Rise vs 2024.** Our 2025/2024 ratio is 1.22× (ΣP) and 1.19× (S1). PNSN windows within 30 km went the other way, 69 → 49 (0.71×). Given the lack of any monthly relation, this does not count against our rise. In short, PNSN neither supports nor contradicts the rise. The KMR finding above also means post-2025-09 rates are, if anything, under-counted.

## Check 1 details (NC.KMR)

Event sets (`kmr_windows_events.parquet`). The top 100 per family are used.

| Set | Detections | Period | Ranked by | Windows with data |
|---|---|---|---|---|
| OLD_all | S1 9× | 2023-12 → 2025-09-15 (old STS-2) | summed cc | 7356 |
| OLD_h1 / OLD_h2 | halves of OLD_all | before / after 2024-11-01 | summed cc | 6936 / 6504 |
| P_OLD | non-build-day candidates | old epoch | P(same) | 7390 |
| NK_B / NK_C | S1 9× from the run without KMR | Horizon / Compact epochs | summed cc | 4055 / 5202 |
| P_B / P_C | candidates from the with-KMR run | Horizon / Compact epochs | P(same) | 5618 / 6602 |

Processing is the t3_cache + v2 chain: 1–8 Hz at 100 sps, then 20 sps, then a 2 Hz high-pass. Each window is 45 s long with S at +15 s, and the template slice runs from S − 1 s to S + 5 s.

Median over 27 families (`kmr_summary_by_set.csv`; per family in `kmr_per_family_set.csv`):

| | OLD_all | OLD_h1 | OLD_h2 | P_OLD | NK_B | P_B | NK_C | P_C |
|---|---|---|---|---|---|---|---|---|
| 2–8 Hz noise, nm/s (StationXML gain) | 10.5 | 10.1 | 10.2 | 9.6 | 132 | 131 | 1.6 | 1.6 |
| direct stack SNR (S window / −14…−5 s) | 3.45 | 2.92 | 3.01 | 2.99 | 0.80 | 0.91 | 0.96 | 1.01 |
| direct stack cc vs old template (±5 samples) | 0.904 | 0.895 | 0.880 | 0.883 | 0.176 | 0.159 | 0.170 | 0.180 |
| cc at lag 0 / best lag in ±2 s | 0.90 / 0.90 | 0.90 / 0.90 | 0.88 / 0.88 | 0.88 / 0.88 | 0.01 / 0.16 | −0.01 / 0.15 | 0.01 / 0.15 | 0.02 / 0.16 |
| channel permutation: diagonal / best off-diagonal \|cc\| | 0.90 / 0.57 | | | | 0.27 / 0.34 | 0.24 / 0.31 | 0.27 / 0.32 | 0.28 / 0.34 |
| horizontals: cc at 0° / at best rotation (angle) | 0.92 / 0.92 (0°) | | | | 0.16 / 0.22 (135°) | 0.14 / 0.20 | 0.18 / 0.20 (170°) | 0.17 / 0.22 |
| Text S1 recipe: cc50 / cc vs old template | 0.52 / 0.89 | 0.47 / 0.90 | 0.47 / 0.88 | 0.47 / 0.87 | 0.48 / 0.19 | 0.48 / 0.18 | 0.32 / 0.17 | 0.32 / 0.18 |
| scrambled control, Text S1 recipe: cc50 / cc vs old | 0.25 / 0.30 | 0.25 / 0.29 | 0.24 / 0.30 | 0.24 / 0.27 | 0.49 / 0.21 | 0.50 / 0.20 | 0.32 / 0.18 | 0.32 / 0.17 |
| scrambled control, direct stack cc vs old | 0.22 | 0.23 | 0.23 | 0.21 | 0.18 | 0.19 | 0.18 | 0.17 |
| energy in S…S+3 s / pre-arrival (`kmr_envelope_summary.csv`) | 1.63 | | | 1.44 | 1.02 | 1.05 | 0.99 | 0.98 |

- Pipeline check: Text S1 re-stacks from the same windows match the run's own new templates only at cc 0.26–0.30 (Horizon) and 0.42–0.48 (Compact). Stacks of noise do not reproduce each other well, which fits there being no signal. The run's templates `tpl_kmrB/C` match the old template at 0.130 / 0.131 and match each other at 0.25.
- The old-epoch direct stacks use detections made with the old KMR template in the network sum (KMR is 1 of 17 stations). So 0.90 may be slightly optimistic. The held-out-style P_OLD set gives the same 0.88.
- Raw data health per month: `kmr_health_monthly.csv`, `kmr_raw_levels.csv`, `kmr_psd_monthly.npz`. Control: on NC.KMPB HHZ with the same method, 2–8 Hz is 4.5–9 counts and 0.1–0.5 Hz is 190–470 counts in 2024-12, 2025-12 and 2026-06, so nothing like the KMR change.
- StationXML (`NC.KMR_response.xml`): same coordinates in all epochs. The Compact Posthole is listed at **71 m depth** (elevation 844 m vs 915 m), so it is a borehole install. Azimuths are nominal 0/90 and the dip of Z is −90. A borehole sensor's true orientation can differ from nominal, but rotation does not recover any signal (see the table).
- Rate table: `kmr_rate_effect.csv`.
- Figures:
  - `kmr_lag_and_cc.png`: lag curves; direct-stack cc against the scrambled control.
  - `kmr_envelope_stacks.png`: energy stacks over the whole window.
  - `kmr_example_stacks.png`: families 1, 8, 19, 23 (old template, stacks, new templates).
  - `kmr_spectra.png`: amplitude spectra per epoch.
  - `kmr_raw_examples.png`: raw counts.

Not done: there is no 16-station run before the swap (it needs continuous data, now deleted). So the size of the step from losing a working KMR is unknown. It could be estimated later from the cut windows by re-summing the network cc at candidate times with and without KMR, but the MAD threshold needs the continuous cc, so only approximately.

## Check 2 details (PNSN tremor)

- **Box, chosen before looking at any rate:** a circle of 30 km radius around the centroid of Shelly's 27 Table S1 family locations (`results/families_south.csv`, 40.062°N, −123.627°W). The families span about 15 km E–W. Also reported: 15 km and 60 km circles.
- **Catalogue:** `catalogs/pnsn_tremor_cascadia_full.csv`, 2.5-min windows, up to 2026-05-27. May 2026 is scaled to a full month.
- **Monthly series 2018-01 → 2026-05** are in `pnsn_monthly.csv`:
  - PNSN windows per circle and their median num_stas;
  - his monthly count (2018-01 → 2024-11);
  - ours: ΣP(same) and S1 9×, both as daily rate × days in month, from `results/sh_conf/monthly.csv`, forward run from 2023-12.
- Yearly totals are in `pnsn_yearly.csv`. PNSN windows within 30 km per year: 2018 204, 2019 83, 2020 90, 2021 134, 2022 101, 2023 **0**, 2024 69, 2025 49, 2026 (Jan–May) 20. His events per year: 4587, 7019, 5760, 7345, 13496, 9862, 13372 (2024 to Nov).

Fit of his monthly count on PNSN windows, 2018-01 → 2023-11 (`pnsn_fit.csv`):

| circle | slope (events per PNSN window) | intercept | Pearson r | Spearman |
|---|---|---|---|---|
| 30 km | −2.8 | 688 | −0.08 | −0.17 |
| 15 km | −2.5 | 681 | −0.06 | −0.14 |
| 60 km | −2.9 | 719 | −0.23 | −0.29 |

Tests against this (useless) fit, 30 km, ours scaled to his by the join-year ratio (ΣP × 1.03, S1 × 1.09). Full table in `pnsn_relation_tests.csv`.

| period | observed / predicted | r with PNSN |
|---|---|---|
| join year, his | 1.77 | 0.67 |
| join year, ours ΣP | 1.77 (by construction) | 0.71 |
| 2024-12 → 2026-05, ours ΣP | 1.95 | −0.05 |
| 2025 | 2.13 | 0.26 |
| 2026-01 → 05 | 1.72 | 0.11 |

The join-year correlation of 0.67–0.71 rests on one PNSN burst (2024-08). The "observed / predicted" values only repeat the rise of his own catalogue since 2018, which PNSN does not show.

Day-level coincidence (`pnsn_daily_coincidence.csv`). PNSN tremor days are days with ≥1 window in the circle. The ratio is the mean LFE count per day on PNSN tremor days divided by the mean on other days.

| series | 30 km: tremor days / ratio | 15 km ratio | 60 km ratio |
|---|---|---|---|
| his 2018-01 → 2023-11 | 86 of 2160 / 2.64 | 2.91 | 1.68 |
| his join year | 7 of 366 / 3.47 | 3.47 | 2.63 |
| ours ΣP join year | 7 of 366 / 3.78 | 3.78 | 2.82 |
| ours ΣP 2024-12 → 2026-05-27 | 10 of 543 / 3.31 | 3.31 | 1.74 |
| ours S1 2024-12 → 2026-05-27 | 10 of 543 / 3.11 | 3.11 | 1.65 |

These rest on few days (7–10 tremor days within 30 km per period), so they are indicative only.

Figures:
- `pnsn_vs_rates.png`: PNSN windows, his and our monthly counts, and observed/predicted.
- `pnsn_scatter.png`: monthly scatter per circle.

PNSN tremor detection itself depends on its station set: the median number of stations per window near Mendocino was 7 → 6 → 5 → 7 over 2018–2026. That, together with its sparse coverage at 40°N, is the most likely reason it does not track the LFE counts.
