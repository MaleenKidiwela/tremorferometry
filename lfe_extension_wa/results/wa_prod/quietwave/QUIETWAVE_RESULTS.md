# Quiet-day waveform test of the WA production detections (NW, core network) — 2026-10-06

## Plain summary

- **The quiet-day NW set at P ≥ 0.85 (999 detections) is mostly real.** Outside the 12-s template window (where the
  P wave and the coda sit, and where the detector never looked), these windows match an independent tremor-day
  stack of the same family and station about as well as tremor-day detections at the same P do: **0.98 of the
  tremor-day level [90 % range 0.87–1.07]**. Scrambled-control, reversed-control and random-time windows score ≈ 0 on
  the same test.
- **But 70 % of that set is "quiet" only because no tremor catalogue exists for those days**: 574 are before
  2009-08-06 (outside the WECC ETS windows) and 121 are after 2026-05-27, when the PNSN catalogue on disk ends
  (all 121 on 2026-06-28 to 07-01 and 09-03 to 09-05). On those days the detections look 100 % real (1.04
  [0.92–1.12]); they are probably tremor the day classes miss, not inter-ETS repeats.
- **On days that really are quiet (catalogue era, 2009-08-06 to 2026-05-27) the picture is weaker**: P ≥ 0.75 →
  0.43 of tremor-day level [0.30–0.55]; P ≥ 0.85 → 0.80 [0.59–0.96]; P ≥ 0.9 → 0.90 [0.73–1.19]. As a rough false
  share: ~55 % at 0.75, ~20 % (4–41 %) at 0.85, ~10 % (0–27 %) at 0.9.
- **What the two controls catch:**
  - **Scrambled control:** in most cases one site has a real signal and the other sites line up by chance. Its
    dominant site carries ~75 % of the sum (random times: 47 %), and that site's own unshifted time falls within
    ±2 s of a data candidate in 76–81 % of cases (random times: 5–6 %). That time matches a P ≥ 0.5 detection of the same
    family in 16 % (P ≥ 0.75) to 28 % (P ≥ 0.85) of cases (random times: 0 %). Two or more sites each landing on a
    P ≥ 0.5 detection: 0.8 % (chance level 0.3 %). The windows have normal energy and only about a tenth of real-LFE
    evidence outside the template window. So these are mostly noise sums carrying one site's real LFE.
  - **Reversed control:** large non-LFE transients. Correlation with the forward template is 0.00, the same as
    random times. Coherence outside the window is 0. Window/pre-window energy has a median of 1.26–2.0 (data
    1.04–1.08, random 1.00), envelope peak/median is 12–16 (data about 6), and 9–10 % contain a spike. The stacks
    show the template window sitting just in front of the onset of a large arrival. These are not catalogued local
    earthquakes (2007–09 ANSS M ≥ 1: 0.6 % vs 0.1–0.5 % for random times). They sit near low-P data candidates (any
    NW candidate within ±30 s: 86–90 % vs 39–41 % random), but only 2–4 % sit near a P ≥ 0.85 one.
- **Which false share is closer?** For the whole quiet set at 0.85 the waveforms give about 2 % [−7 %, 13 %], so
  **the scrambled 5 % is closer than the reversed 16 %**. On the truly quiet subset the waveforms give about 20 %, which
  falls between the scrambled control (8 %) and the reversed control (26 %), and the range covers both. At 0.75 on
  that subset (~57 %), the reversed control (48 %) is closer than the scrambled one (30 %). Neither control is a copy of
  the data's false alarms. Treat them as a bracket: scrambled as the low end, reversed as the high end.
- **Recommendation (not a decision):** for NW counts on catalogue-quiet days, use **P ≥ 0.85 with a ~10–25 %
  false-share bracket**, or **P ≥ 0.9** if inter-ETS rates need ≲ 10 % (quiet-day families: 58 at 0.9 vs 94 at 0.85; catalogue era 36 vs 57). Do not use
  0.75 on quiet days. Separately, the day classes should mark days with no tremor catalogue as "unknown" rather than "quiet".
  This applies to pre-2009-08 days outside the WECC windows and to everything after 2026-05-27.

Caveats: the real-share numbers use tremor-day detections at the same P as "100 % real". That holds if quiet-day
LFEs are as strong as tremor-day LFEs of the same family. Weaker inter-ETS events would push the estimate down, so
the false shares above may be too high. The 90 % ranges come from a day-block bootstrap and do not cover that
assumption. The catalogue-era 0.85 subset is small (304 data candidates, of which 239 can be compared with the
tremor-day stack).

## What was done

- **Exact replay of the production matched filter.** Inputs: candidates from `conf/P/prod/P_*.parquet` (net core,
  reg NW, `role` not sealed), with day classes from `wa_prod_post.day_class` (quiet = no PNSN tremor within 60 km;
  before 2009-08-06, outside the WECC windows).
  - Data and random times: each station window starts at `i + round(offset_s·20)`, the template moveout
    (`wa_prod_mf.run_day`).
  - Reversed control: same position, compared with the reversed template.
  - Scrambled control: `((i − shift_site) % NV) + d`, where `shift_site = run_<date>.json shifts['core|c1|<site>']`
    (`wa_prod_common.site_shift`, one per site and date).
  - Stations: the day's core members after the drop-dead rule. Weights: production site merge B.
  - Waveforms: the 20 Hz 1–8 Hz caches the matched filter read (`find_cache`).
  - Check: the recomputed network cc equals the stored `net_cc` for all 28,502 candidates (max |Δ| 2.5e-4;
    `tables/replication_check.json`).
- **Sets:**
  - Quiet days: all P ≥ 0.75 rows (data 2,956, scrambled 634, reversed 1,012). Subsets are taken at 0.85 and 0.9.
  - Tremor reference (tremor30): every control row (scrambled 257, reversed 89) plus data on 200 sampled days
    (10,299).
  - Random times: one random time of the same family and day for every data row.
  - The tremor-day reference stack (tremor data P ≥ 0.85, even days of the month) is independent of every set
    scored against it. Tremor rows are scored on odd days only.
- Scripts: `scripts/qw_extract.py`, `qw_windows.py` (24 CPU processes, ~13 min), `qw_analyze.py`, `qw_extra.py`.
  No GPU, downloads, deletions or git were used. The production queue and its markers were not touched.

## Test 1 — stacks (`tables/t1_stacks_summary.csv`, per family × station: `t1_stacks_per_family_station.csv`)

Families with ≥ 5 events in the group. "Netw" means coef-weighted over the family's stations, then the median
over families. SNR is the stack rms in [0, 12 s] / rms in [−20, −8 s]. "Outside" is the stack correlation with the
independent tremor-day stack, using [−8, 0) s and [12, 24) s only.

| quiet days, all eras | families | stack cc fwd (netw) | stack SNR (netw) | 5-event stack cc / SNR | outside cc (median) |
|---|---|---|---|---|---|
| data P ≥ 0.75 | 88 | 0.36 | 1.51 | 0.28 / 1.25 | 0.028 |
| scrambled P ≥ 0.75 | 32 | 0.33 | 1.22 | 0.28 / 1.12 | 0.000 |
| reversed P ≥ 0.75 | 57 | 0.00 | 5.14 | 0.00 / 4.97 | −0.001 |
| random (data's family-days) | 88 | −0.00 | 1.02 | 0.00 / 1.02 | −0.002 |
| data P ≥ 0.85 | 46 | 0.42 | 1.62 | 0.34 / 1.35 | 0.054 |
| scrambled P ≥ 0.85 | 4 | 0.36 | 1.20 | 0.32 / 1.17 | 0.008 |
| reversed P ≥ 0.85 | 6 | −0.01 | 3.82 | −0.01 / 3.73 | −0.002 |
| data P ≥ 0.9 | 19 | 0.44 | 1.67 | 0.39 / 1.52 | 0.076 |
| tremor-day data P ≥ 0.85 (odd days) | 71 | 0.43 | 1.56 | 0.33 / 1.24 | 0.060 |

Reading: the correlation inside the template window **cannot separate data from the scrambled control** (0.36 vs
0.33). Every candidate was selected for matching there. The energy and the out-of-window match do separate them:
quiet data look like tremor-day data, and the scrambled control sits close to random times. Reversed-control
stacks are energetic but contain no forward-template signal. Example families (template plus the four stacks, in
units of each stack's pre-window noise): `figures/example_family_{3,36,83,101}.png`. In family 3 the
scrambled-control stack shows a clear S arrival at the B013/B943 site, which is the "one real site" pattern.

## Test 2 — single windows (`tables/t2_single_windows.csv`, per candidate: `t2_network_per_candidate.parquet`)

| quiet, P ≥ 0.85 | station cc fwd median [p25, p75] | network cc fwd (weighted) median | window/pre energy median | share with window/pre > 1.2 | env peak/median | spike share |
|---|---|---|---|---|---|---|
| data (999) | 0.052 [−0.008, 0.156] | 0.209 | 1.063 | 0.131 | 6.0 | 0.004 |
| scrambled (53) | 0.035 [−0.009, 0.126] | 0.187 | 1.012 | 0.000 | 5.2 | 0.000 |
| reversed (163) | −0.001 [−0.039, 0.036] | −0.003 (reversed template: 0.175) | 1.533 | 0.595 | 14.4 | 0.092 |
| random (999) | 0.000 [−0.042, 0.042] | −0.001 | 1.001 | 0.015 | 5.3 | 0.001 |
| tremor-day data (5,503) | 0.056 | 0.191 | 1.062 | 0.148 | 7.3 | 0.003 |

`figures/single_window_netcc_quiet.png` (network cc) and `figures/single_window_energy_ratio.png` (energy) show
the same pattern at 0.75 and 0.9.

**Selection-free real share** (`tables/t4_selection_free_summary.csv`). This is each window's out-of-window match
to the tremor-day stack, averaged over the set and divided by the family-matched tremor-day mean at the same P
(1 = as real as tremor-day detections). 90 % day-bootstrap ranges are in brackets.

| quiet days | P ≥ 0.75 | P ≥ 0.85 | P ≥ 0.9 |
|---|---|---|---|
| data, all | 0.81 [0.71–0.90] | 0.98 [0.87–1.07] | 1.05 [0.95–1.16] |
| data, catalogue era (truly quiet) | 0.43 [0.30–0.55] | 0.80 [0.59–0.96] | 0.90 [0.73–1.19] |
| data, no catalogue (pre-2009-08 / post-2026-05) | 1.01 [0.93–1.09] | 1.04 [0.92–1.12] | 1.11 [1.00–1.24] |
| scrambled, all | 0.10 [0.04–0.17] | 0.18 [−0.01–0.37] | (2 rows) |
| reversed, all | 0.02 [−0.03–0.06] | −0.01 [−0.08–0.09] | −0.01 |
| random, all | 0.00 | −0.01 | −0.04 |

The energy version (`real_frac_esnr`) agrees for data (0.99 / 1.02 / 1.09 overall; 0.84 / 0.82 / 0.98 in the
catalogue era). Transient-flagged data candidates (envelope peak/median > 12 or a spike) still carry LFE evidence:
0.74 vs 0.84 for the rest in the catalogue era at 0.85 (`tables/t4_transient_split.csv`).

## Test 3 — timing and earthquakes (`tables/t3_timing_summary.csv`, `t3_earthquake_coincidence_2007_2009.csv`, `t3_day_activity.csv`)

**Scrambled control.** "Dominant site" = the site with the largest weighted cc. Its real time is
`(i − shift) % NV`, checked against data candidates within ±2 s.

| quiet, P ≥ 0.85 (53) | scrambled | random times, same shifts applied |
|---|---|---|
| dominant site on a same-family data detection with P ≥ 0.5 | 0.28 (0.16 at P ≥ 0.75) | 0.00 |
| dominant site on any NW detection with P ≥ 0.5 | 0.34 (0.21) | 0.01 |
| dominant site on any NW candidate (≥ 8 MAD, any P) | 0.81 (0.76) | 0.06 |
| ≥ 2 sites each on a NW detection with P ≥ 0.5 | 0.00 (0.008) | 0.006 (0.003) |
| ≥ 2 sites each on any NW candidate | 0.45 (0.49) | 0.13 |
| dominant-site share of the sum (median) | 0.74 | 0.47 |

**Reversed control (quiet).**

| | P ≥ 0.75 | P ≥ 0.85 | random |
|---|---|---|---|
| within ±2 s of any NW candidate | 0.13 | 0.11 | 0.06 |
| within ±30 s of any NW candidate | 0.86 | 0.90 | 0.39–0.40 |
| within ±30 s of a NW detection with P ≥ 0.5 | 0.34 | 0.32 | 0.04 |
| within ±30 s of a NW detection with P ≥ 0.85 | 0.02 | 0.04 | 0.01 |
| ANSS M ≥ 1 event in [t − 60, t + 30] s (2007–09 only) | 0.006 | 0.000 | 0.001–0.003 |

**Day activity.** Quiet-day data detections cluster on active days: the median day has 23 NW P ≥ 0.85
detections, and 72 % of them fall on days with ≥ 10. The controls fire on inactive days (median 0; 1–4 % on days
with ≥ 10). If the data detections were noise, the controls should fire on the same days. They do not.

## Files

- `QUIETWAVE_RESULTS.md` (this file); `tables/` (all numbers); `figures/` (4 example families, 2 distributions).
- `rows.parquet` (every scored candidate, with `row_id`); `nw_core_candidates.parquet` and
  `core_data_all_regions.parquet` (extracts from the P tables); `work/` (per-day windows, float16, 3.4 GB;
  `work_smoke/` = 1 smoke day). The whole folder is 6.5 GB. It can all be regenerated from the scripts, so the
  extracts and `work/` can be removed once you OK it.
