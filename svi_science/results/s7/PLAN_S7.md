# S7 — LFEs between ETS episodes: PLAN (frozen before any S7 number)

Written 2026-10-07 ~3:30 pm PT, before any S7 statistic was computed. (Only the size of the candidate pool was looked
at while choosing the definition: P ≥ 0.75 de-duplicated data, outside episodes ± 2 d = 721, of which 109 have no PNSN
epicentre within 30 km of the family that day; scrambled control 11 and 3.)
Script `svi_science/scripts/s7_between.py`; outputs `svi_science/results/s7/`, figures `svi_science/figures/s7_*.png`.

## Unit and period
- Event = one row of `s0/analysis_table.parquet` (run prod) with series data, P_same ≥ 0.75, keep_xfam8 (counted once
  across neighbouring families). Scrambled control = the same rows of series c1 (stations shifted to different days),
  through every step. P ≥ 0.5 = column only (its control is ~half of the quiet-day data).
- Period 2014-01-01 … 2026-05-27 (UTC days). The PNSN tremor catalogue ends 2026-05-27, so after that date neither the
  episodes nor the tremor test are defined; later events are counted in a separate line and not used.
- Dates/joins are UTC days (as the catalogue and PNSN exposure); times shown to the reader are Pacific time.

## Definition of "between-episode" (fixed here, used everywhere)
**BE (headline):** the event's UTC day is outside every one of the 79 SVI episode windows widened by 2 days
([start − 2 d, end + 2 d], `episodes_svi_final.csv`, majors and minors), AND no PNSN epicentre lies within 30 km of the
event's family on that day (`n_trem30 == 0`).
Columns (reported side by side, never substituted): BE_broad = outside episodes ± 2 d with any local tremor;
`n_tremor_day == 0` (no PNSN epicentre in the whole family box that day); `quiet` flag (box-quiet and ≥ 3 d from any box
tremor day); outside the 79 episodes without the 2-day pad.
Eligible family-days (the exposure for every rate and shuffle): family × UTC day in the period, family observable that
day (`famday_detectability`, run prod), satisfying the same BE condition for that family.

## Tests and pass/fail rules
1. **Inventory.** Counts data vs control per definition; per family and per group (map, depth), per year (timeline).
   *Clustered in time*: number of distinct UTC days holding BE events vs a family-preserving day shuffle (each event moved
   to a random eligible day of its own family, 1000×). Clustered if observed distinct days < 1st percentile (p < 0.01).
   Bursts = events chained by gaps < 24 h; reported: number, sizes, share of events in bursts of ≥ 3.
2. **Next ETS.** τ = days from the event to the start of the next MAJOR episode (39 majors; events after the last major
   start are censored). Superposed-epoch rate = events per eligible family-day in τ bins 0–7, 7–14, 14–30, 30–60, 60–120,
   120–240, > 240 d. Statistic F30 = share of events with τ < 30 d. *Ramps up before ETS* if F30 > 99th percentile of the
   family-preserving day shuffle AND the excess holds for events ≥ 30 d after the previous major's end (so it is not the
   tail of the last ETS). The same for days since the previous major end (tail test), reported alongside.
   *Sits where the next ETS starts*: events with τ ≤ 60 d before one of the 10 SVI majors with an S3 LFE start position
   (`lfe_start10`); statistic = median |along-strike(family) − start|; control = events re-paired with a random other
   one of those 10 episodes (5000×). Pass if p < 0.01 with ≥ 20 events; fewer than 20 events = "too few to test".
   Column: the same with a PNSN-tremor start (median along-strike of the first 10 % of box epicentres) for all 39 majors.
3. **Small slip / tremor elsewhere.** Share of BE events whose family has ≥ 1 PNSN epicentre 30–100 km away on the same
   UTC day, vs the same share over the family-preserving shuffle; tied to distant tremor if above the 99th percentile.
   Column: tremor 30–100 km within ± 6 h of the event.
   **GNSS stack** (Merlin S2 step 7): clusters of BE days (dates < 3 d apart merged; cluster day = day with most events).
   Offset = M5[d + 1] − M5[d − 5] (5-day medians after minus before) at the 14 SVI-box stations (S1 set), projected on
   each station's own ETS azimuth (S1), station mean; north subset as a column. Days whose [d − 6, d + 6] touches an
   episode window ± 2 d are dropped (same for controls). Control: 1000 sets, each cluster day replaced by a random allowed
   day within ± 30 days of year. Slip-like signal only if the mean offset is positive (ETS direction) and > 3σ of the
   control; otherwise the 3σ value is the bound (mm per ~6 days), compared with ~5 mm for a big ETS. Reference column:
   the same estimator on the first day of the 40 minor episodes (can the stack see small slip at all?).
4. **Down-dip.** ETS share per family = its event count inside the 39 majors. Statistic = mean Bostock depth (and mean
   down-dip km) of BE events vs 10,000 draws of the same number of events from the families' ETS shares. *Deeper than
   ETS* if observed > 99th percentile. Column: Spearman of per-family (BE + 0.5)/(ETS + 0.5) vs depth over families with
   ≥ 10 ETS events; per-group table.

## Checks carried (not filters)
- Fixed-13-station run (prod_S13, to 2021-03-31): inventory and the depth test repeated.
- Scrambled control through tests 1, 2 and 4 (counts too small for statistics; listed).
- n_sta_day of BE events vs eligible days (is BE activity just better-detected days?).
- Detection-problem flags of S0 carried in the event table.
- No Merlin review of this plan (caller asked to run directly); results provisional until a control or NW repeats them.
