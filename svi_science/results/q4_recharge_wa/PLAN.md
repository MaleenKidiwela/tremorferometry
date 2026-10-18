# Q4 — Washington-NW confirmatory test of the slow deep "recharge" (PLAN, frozen before any Q4 number)

Written 2026-10-07 ~2:30 pm PT, before any Q4 statistic was computed. Only file schemas, family positions and the
episode lists were looked at while writing it (no WA rate-vs-time-since-ETS numbers).
Script `svi_science/scripts/q4_recharge_wa.py`; outputs here; figure `svi_science/figures/q4_recharge_curves.png`.
"Time-reversed control" = templates reversed in time; "scrambled control" = station data shifted to other times of day.

## Question
SVI (P3 step 4): deep activity between big ETS is low right after a big ETS and builds for > 1 year (big-ETS clock,
D30 0.69, time-shift p 0.002, τ ≈ 322 d (85–400), fit at its 400-d limit; the 39-major clock fails, p 0.077; the two
halves disagree). Does Washington-NW (every day 2007-01-01 … 2026-05-27) show the same slow deep rebuild, and on the same
time scale?

## Data and counting (same as P3 step 4 / A, adapted to WA)
- **Event:** core-network NW candidate (`wa_prod/conf/P/prod/P_*.parquet`, net core, reg NW) with
  `P_own_WAVE_PNET ≥ 0.75` (= catalogue_v2 headline own_pnet ≥ 0.75), **counted once by dedup variant D**
  (`post_v2/dedup/flags/`). Series: data, time-reversed control (rev), scrambled control (c1); all run identically.
- **Families / positions:** `svi_science/results/a_wa_between/families_wa.csv` (copies take their source's primary
  position, `families_v1_3.csv`).
- **Deep (headline):** position-based — the deepest third (Slab2 depth at the primary epicentre ≥ its 2/3 quantile) of
  the A/B-located families. NOT Hyp2000 depths (RESULTS_A caveat).
  Columns: deepest third by straight-frame down-dip km (A/B); deepest third by Slab2 depth over all families incl. C.
- **Exposure:** scanned family-days of the deep families (`post_v2/dedup/rates_daily_v2_dedup.parquet`, NW family rows,
  zero rows included). Rate = counted-once deep events per family-day.
- **Activity counted:** ALL deep LFEs on days outside every MAJOR episode (69 majors of
  `episode_rates/episodes_nw_final.csv`), whatever the clock — so minor-episode and local-tremor-day events count, ETS
  events never do. Days inside majors are dropped from both events and exposure.
  - Note found while reading P3's code (before any WA number): P3's big-ETS-clock row dropped only the days of the BIG
    episodes, so deep events inside non-big majors were counted at later t. Q4's headline drops all majors (the stated
    definition). The P3-as-coded variant (drop only the reference episodes) is a column for WA, and SVI is re-run both
    ways (column "SVI re-check", P3 data and code copied unchanged otherwise).
- **Period / eras:** E1 2007-01-01 … 2011-12-31, E2 2012-01-01 … 2014-12-31, E3 2015-01-01 … 2026-05-27 (episode list
  ends with PNSN coverage). Days are assigned to eras by date; cycles keep their identity across an era boundary.
  Before 2009-08 only the three WECC windows are listed (all big); unlisted minor episodes then simply count as activity.

## Clocks and episode ends (as P3)
- **Big-ETS clock (headline):** NW episodes with `big == True` in `b_ets_record/episodes_b.csv` (primary two-class split
  among majors + WECC windows; 31 episodes 2007–2025). Columns: a2 (one per year) list; all 69 majors as clock.
- t = days since the reference episode ended; days 1 … 400 (TMAX = 400 as P3); bins 0–7, 7–14, 14–30, 30–60, 60–90,
  90–150, 150–240, 240–400.
- **Ends:** D1 = episode end (headline); D2 = end + 7 d (t' = t − 7, t ≥ 8); D3 = the family's own end (last day inside
  the reference episode with ≥ 1 PNSN epicentre within 30 km of the family, `wa_prod/post/tremor_family_day.parquet`;
  episode end where no such day / before 2009-08), family-day table, no surrogate (as P3).

## Statistics
- **D30** = 1 − rate(t 1–30) / rate(t 31–400). Positive = quieter in the first month.
- **β** = slope of a Poisson fit rate ∝ t^β over t = 1 … 400 (exposure offset). Positive = rises with time since the ETS.
  Column β_late = same over t = 31 … 400 (is the rise still going after the first month?).
- **Fit** rate = r∞ (1 − a e^(−t/τ)) (P3 code; τ bounded 1–400 d); τ 95 % range from 1,000 bootstraps over cycles.
  τ ≥ 395 d is marked "at the limit" (no levelling-off within 400 d).
- **Time-shift control (p):** each cycle's daily series (events and exposure together) rotated by a random lag within the
  cycle, 1,000×, recompute D30 and β; p = (#surrogate ≥ observed + 1)/1001 (P3 surrogate). Column: p from a cycle
  permutation (each cycle's daily series replaced by another cycle's of ≥ the same length, truncated) is NOT run — the
  rotation is the only pre-registered control.
- **Time-reversed / scrambled controls:** same deep set, same tables; events and D30 reported. If a control exceeds 10 %
  of the data events it is said plainly.
- **Era-pooled shape normalised per era:** each era's exposure multiplied by that era's mean rate (over the usable days of
  that table), so the pooled rate is "relative to its own era"; same fit / D30 / β / surrogates / bootstrap. Reported,
  not part of the rule.
- **PNSN deep tremor (detection-independent):** `catalogs/pnsn_tremor_cascadia_full.csv`, epicentres within 15 km of a
  headline deep family (P3 used 15 km of deep families), per day (n = 1 per covered day), 2009-08-06 … 2026-05-27, same
  clocks/ends, majors' days dropped; per era and pooled.
- **Fixed 9-borehole run:** if `wa_fixed10/post_v2_core9/catalogue_v2.parquet` is complete when the main run finishes, a
  column with its data rows (P_own ≥ 0.75, NW, headline deep families via family_uid; raw counting, since no dedup
  flags exist for it; exposure = deep families × days). If not complete, noted and not used.

## Pass rule (fixed now)
**Recharge SUPPORTED in Washington** if ALL of:
1. E3, big-ETS clock, D1, headline deep set: D30 > 0 with time-shift p < 0.01, AND β > 0 with time-shift p < 0.01;
2. E1 and E2 (same clock/definitions): D30 > 0 and β > 0 (sign only);
3. PNSN deep tremor (big clock, D1, 2009-08-06 … 2026-05-27 pooled): D30 > 0 and β > 0 (sign only).

**NOT SUPPORTED** if E3 D30 ≤ 0, or E3 β ≤ 0, or E3 p(D30) ≥ 0.05.
Anything else = "same direction, not confirmed", naming the part that misses.

**Time scale vs SVI:** τ (E3 and era-pooled) with its bootstrap 95 % range. "Same time scale" if the WA range overlaps
SVI's 85–400 d; "different" if not. If both are at the 400-d limit the answer is "both still rising at 400 d" (no
measured τ on either side).

## Columns (never change the verdict)
All-majors clock; a2 clock; D2; D3; P3-as-coded exclusion; down-dip / incl.-C deep sets; P ≥ 0.85; source-group-once and
raw counting; live-core-days only (days with ≥ 50 % of the era's median live core count, `a_wa_between/live_core_by_day.csv`);
SVI re-check (P3 big clock with all 39 majors dropped).
