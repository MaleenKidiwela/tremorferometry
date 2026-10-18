# Q5 plan — are the deep between-ETS LFE bursts small slow slips?

Written 2026-10-07 2:25 pm PT, **before any Q5 number was computed**. Script `svi_science/scripts/q5_bursts.py`
(CPU, < 10 GB, no downloads). Outputs here; figures `svi_science/figures/q5_*.png`.

## Inputs (frozen, re-used as they are)
- **SVI** (primary, 2014-01-01 … 2026-05-27): the 76 BE events `results/s7/between_events.csv` (P ≥ 0.75, counted once,
  outside the 79 episodes ± 2 d, no PNSN epicentre ≤ 30 km of the family that day) and the 688 BE_broad events
  `results/s7/between_broad_events.parquet` (outside episodes ± 2 d, any local tremor). Family positions
  `results/s0/families.csv`. Eligible family-days rebuilt exactly as `s7_between.py` (famday_detectability prod
  observable ∩ exposure_svi; BE pool = outside pad and n30 = 0; BE_broad pool = outside pad).
- **WA-NW** (primary E3, 2015-01-01 … 2026-05-27): the 220 BE events `results/a_wa_between/between_events.csv` and the
  934 BE_broad events (`between_broad_events.parquet`, series data). Positions `families_wa.csv`. Eligible pools rebuilt
  exactly as `a_wa_between.py` from `exposure_famday.parquet` (BE pool day_class quiet, outside pad; broad pool day_class
  ≠ unknown, outside pad). **Secondary**: E1/E2 events `between_events_E1E2.csv` (172) — used only in tests 2–4 and as
  a secondary column of test 1 (PNSN catalogue covers them).
- PNSN tremor `catalogs/pnsn_tremor_cascadia_full.csv` (ends 2026-05-27).
- GNSS: the S7/A estimator (Gualandi X.mat daily series; offset = median of days d+1…d+5 minus median of d−5…d−1,
  projected on each station's own ETS azimuth from S1 Tier A, mean over SVI-14 or WA-12 stations).
  Green's functions `lfe_extension/results/slip_compare/slip_cache.npz` (G, cutde half-space, Slab2 mesh), μ = 40 GPa
  (the slip_build value).

## Burst definition (same as S7 / A)
- Burst = events chained by gaps < 24 h. "Big burst" = ≥ 3 events. "Multi-family burst" = ≥ 3 events from ≥ 2 families.
- Counts: SVI 40 bursts (7 big); WA E3 62 (14 big). Unit of inference for every test is the **burst**, because events
  inside a burst are not independent; event-level numbers are columns.

## Test 1 — coincident PNSN tremor by distance
- For each event: number of PNSN epicentres at 0–30, 30–60, 60–100 km from its family, (a) same UTC day, (b) ± 6 h.
- BE events have no tremor ≤ 30 km that day **by construction**, so for BE the bins 30–60 and 60–100 km are tested; the
  0–30 km bin is tested on BE_broad (between episodes, local tremor allowed) against its own pool.
- A burst is "with tremor in bin b" if any of its events is.
- Scrambled-time control ("random days"): each burst is moved to a random eligible day of its first event's family
  (same pool rule as its set), keeping its number of days; 1000 draws. Event level: each event moved to a random
  eligible day of its family (the S7/A rule).
- **Pass rule (per region, per bin):** burst share above the 99th percentile of the random-day control (p < 0.01) ⇒
  "tied to tremor in that bin". Event level is a column. SVI repeats the WA result if the 30–60 or 60–100 km bin passes
  at burst level in SVI.
- Coverage column: share of major-ETS family-days with ≥ 1 epicentre ≤ 30 km, per region (how well PNSN sees each
  region's tremor). Not a test.

## Test 2 — migration inside bursts
- Bursts with ≥ 4 events whose families span ≥ 2 positions ≥ 5 km apart (else position cannot change).
- Theil–Sen slope of along-strike km vs time (h) and of down-dip km vs time; speed in km/h.
- Control: event times permuted inside the burst (positions kept), 5000 draws; two-sided p on |slope|.
- **Pass rule (per region):** "bursts migrate" if the number of bursts with p < 0.05 (either axis, Bonferroni ×2 → p <
  0.025 per axis) exceeds the 99th percentile of Binomial(n_bursts, 0.05). Speeds of passing bursts are then compared
  with: ETS fronts ~5–15 km/day (0.2–0.6 km/h), rapid tremor reversals ~10–20× faster (≈ 2–10 km/h, Houston et al. 2011),
  along-dip tremor streaks tens of km/h (Ghosh et al. 2010).
- Column: the same slope on PNSN epicentres within 100 km of the burst centroid, burst start − 6 h … end + 6 h, when
  ≥ 5 epicentres (time-permutation control as above).

## Test 3 — duration and size (descriptive, no pass rule)
- Per burst: duration (h), n families, max family separation (km), extent of coincident PNSN tremor (along-strike span,
  same window as test 2).
- Compared with the region's **minor episodes** (non-major in episodes_svi_final / episodes_nw_final, inside the study
  period): span (days), along-strike span of P ≥ 0.75 counted-once LFE families, number of PNSN epicentres.
- Moment–duration scaling for slow earthquakes (Ide et al. 2007: M0 ≈ T × 10^12.5±0.5 N·m/s) gives the Mw a burst of
  that duration would have if it were a slow slip; used in test 5.

## Test 4 — same patches, regular recurrence?
- Per source group: bursts it takes part in, distinct years, share of all bursts.
- Inter-burst intervals for (a) all bursts in the region and (b) the bursts of the top group (SVI G10, WA G019 — named
  now from S7/A, not chosen after looking at intervals): coefficient of variation (CV) of intervals.
- Control: each burst moved to a random eligible day of its first family (1000 draws), CV recomputed.
- **Pass rule:** "regular recurrence" if CV < 1st percentile of the control; "clustered" if CV > 99th percentile; else
  "consistent with random timing". Median interval reported against the ETS repeat time.

## Test 5 — GNSS on the biggest bursts
- Day set: the big bursts (≥ 3 events) of each region; day = the UTC day holding most of the burst's events.
  WA secondary set: add the E1/E2 big bursts.
- Observed: mean projected offset (mm). Control: 1000 sets of season-matched (± 30 day-of-year) allowed days (outside
  all episodes ± 8 d, as S7/A), one per burst. z = (obs − control mean)/σ. "Slip-like" if z > 3.
- Single-burst offsets for the two biggest bursts per region (z against single random days).
- **Bound** (always stated): upper bound = max(obs − mean, 0) + 3σ mm, per region and for both regions combined
  (inverse-variance weighting of the per-burst moment estimates).
- **Mm → Mw:** for each burst, pure-thrust uniform slip on all mesh patches within R = 10 km (5 and 20 km columns) of
  its family centroid; predicted mean projected offset at the same stations per N·m; M0 bound = bound / mean predicted
  offset per N·m; Mw = (2/3)(log10 M0 − 9.1).
- Compared with the Mw the Ide scaling predicts for the bursts' durations: if the bound is above that Mw, the GNSS
  stack **cannot** test the slow-slip idea, and we say so.

## What would answer the question
- **For "small slow slips":** tremor tied at burst level (test 1, both regions); migration at slow-slip-like speeds
  (test 2); durations/extents on the small end of the minor-episode distribution (test 3); repeated deep patch (test 4);
  GNSS positive (test 5, unlikely by design).
- **Against:** no tremor link, no migration, single-family point-like bursts lasting minutes, GNSS bound below the
  scaling prediction.
- Expected outcome: GNSS is a bound only; the tremor and migration tests carry the answer.

Language: "scrambled control" / "random-day control", never "null". Times in PT.
