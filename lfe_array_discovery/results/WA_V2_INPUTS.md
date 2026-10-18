# Inputs for `freeze_wa_v2` — station-hour health census and hour-matched quiet sets

*Written 2026-09-20 by agent I. **This file is an INPUT TABLE, not a freeze.** It contains no
Washington result and no engine constant. `freeze_wa_v2` is written by whoever opens WA_GO; this
document exists so that nothing below has to be re-derived then. `results/freeze_wa_v1.json` and
`scripts/x_step5_queue.sh` were NOT touched, and `results/WA_GO` was NOT created.*

Authorisation: Merlin's Step-1 audit (`results/merlin_review_step1_2026-09-20.md`), Q3 "Legitimate
WA work that does not wait on the north control … (b) the WA **station-health census** and
hour-matched quiet set (CPU, prepares v2)", plus the defect list F1 (no station-health screen) and
F2 (quiet hours {0, 23}, no earthquake screen).

---

## 1. Station-hour health census (defect F1)

Rule, applied verbatim to EVERY station-hour of the caches (tremor and quiet alike), computed on
the frozen 20 Hz 1–8 Hz cache so the Welch band *is* the passband:

| test | veto when |
|---|---|
| flat | any component has zero rms, or > 50 % of the hour equals its own median (gap-filled constant) |
| H/V | horizontal/vertical rms ratio outside [1/30, 30] |
| tone | 1–8 Hz Welch PSD peak/median > 1e4 on any component |
| short | fewer than half the hour's samples present |

Script `scripts/b_wa_health.py`. Outputs: `results/wa_station_hour_veto_<year>.csv` (one row per
station-hour with reason codes), `results/wa_station_veto_summary_<year>.csv`,
`results/wa_array_veto_summary_<year>.csv`.

### 1.1 N-Olympic 2010 (cache tranches T2A–T2D + the T2B quiet days)

147 stations, 9,289 station-days, **222,936 station-hours**.

| array | station-hours | % vetoed | % flat | % H/V | % tone | % short | stations | stations > 30 % |
|---|---|---|---|---|---|---|---|---|
| UW | 3,936 | **4.37** | 0.71 | 3.73 | 0.00 | 0.58 | 3 | 0 |
| XG.DR | 18,936 | 2.79 | 0.01 | 2.79 | 0.00 | 0.00 | 11 | 1 |
| PB | 32,400 | 0.79 | 0.03 | 0.03 | 0.56 | 0.19 | 19 | 0 |
| XU.C | 1,584 | 0.44 | 0.00 | 0.06 | 0.00 | 0.38 | 1 | 0 |
| XG.PA | 17,088 | 0.23 | 0.08 | 0.08 | 0.05 | 0.11 | 12 | 0 |
| XG.TB | 20,208 | 0.12 | 0.00 | 0.00 | 0.00 | 0.12 | 14 | 0 |
| CN | 3,504 | 0.11 | 0.09 | 0.06 | 0.03 | 0.00 | 2 | 0 |
| XU.BS | 17,448 | 0.05 | 0.00 | 0.01 | 0.00 | 0.05 | 14 | 0 |
| XU.W | 7,008 | 0.04 | 0.00 | 0.04 | 0.00 | 0.00 | 4 | 0 |
| XG.BH | 16,728 | 0.02 | 0.01 | 0.00 | 0.00 | 0.02 | 11 | 0 |
| XU.PL | 10,512 | 0.01 | 0.00 | 0.01 | 0.00 | 0.00 | 6 | 0 |
| XG.CL | 29,376 | 0.01 | 0.00 | 0.01 | 0.00 | 0.00 | 20 | 0 |
| XG.GC | 21,768 | 0.00 | — | — | — | — | 13 | 0 |
| XG.SL | 1,008 | 0.00 | — | — | — | — | 3 | 0 |
| XG.LC | 21,432 | 0.00 | — | — | — | — | 14 | 0 |

**No array exceeds 20 % vetoed station-hours** (worst: UW 4.4 %), so Merlin's conditional — "the
site-probe admission was blind to it and the WA coverage numbers in `P_STEP2_SITES.md` are
overstated" — **does not fire for 2010**. The WA 2010 cache carries no SHVB-class station.

Worst individual stations: **XG.DR09 30.1 % (all H/V)** — the only station over 30 %;
PB.B019 10.3 % (all tone); UW.GNW 8.3 % (H/V, and GNW is a dropped/withheld site anyway);
UW.DOSE 3.3 % (short); PB.B007 3.0 %.

**Coverage recomputed anyway** (`scripts/b_wa_coverage.py`, dropping stations > 30 % vetoed, i.e.
XG.DR09 alone), with the same code that produced the published numbers:

| variant | sub-region | cells | admitted sites | % cells ≥5 sites ≤60 km | % ≤30 km | median sites |
|---|---|---|---|---|---|---|
| as published | N-Olympic | 130 | 20 | 94.6 | 50.8 | 10 |
| health-adjusted | N-Olympic | 130 | 20 | **94.6** | 50.8 | 10 |
| as published | S-Olympic | 84 | 20 | 83.3 | 13.1 | 7 |
| health-adjusted | S-Olympic | 84 | 20 | **83.3** | 13.1 | 7 |

and on `freeze_wa_v1`'s own 145 active cells with its 11-site table: 11/11 sites survive, coverage
**97.2 % → 97.2 %**, median 8 sites. (The freeze records 97.93 % for the same quantity; my
recompute of the *published* number gives 97.2 %, a 0.7-point difference that comes from the cell
depth column used for the hypocentral distance — worth reconciling when v2 is written, but it does
not change any decision.)

**Consequence for v2:** adopt the station-hour veto table as a per-hour mask (`W_HEALTH_VETO`), not
as a station-level exclusion — the census says the WA problem is hour-shaped, not station-shaped.

### 1.2 Puget / S-WA 2008 (cache tranche T4)

83 stations, 5,690 station-days, **136,560 station-hours**.

| array | station-hours | % vetoed | % flat | % H/V | % tone | % short | stations | stations > 30 % |
|---|---|---|---|---|---|---|---|---|
| UW (GNW only) | 1,800 | **12.22** | 2.56 | 8.72 | 1.28 | 1.83 | 1 | 0 |
| ZK | 480 | 7.50 | 0.00 | 0.00 | 0.00 | 7.50 | 2 | 0 |
| YW | 9,648 | 6.99 | 0.28 | 5.67 | 0.87 | 0.27 | 6 | **1** |
| XU.BD | 10,920 | 1.01 | 0.35 | 0.04 | 0.05 | 0.61 | 6 | 0 |
| XU.BS | 10,632 | 0.89 | 0.24 | 0.19 | 0.08 | 0.57 | 6 | 0 |
| XU.PL | 9,888 | 0.86 | 0.50 | 0.03 | 0.01 | 0.34 | 6 | 0 |
| XU.E | 8,424 | 0.77 | 0.46 | 0.12 | 0.17 | 0.26 | 5 | 0 |
| XU.N | 20,928 | 0.65 | 0.16 | 0.09 | 0.03 | 0.44 | 13 | 0 |
| PB | 23,736 | 0.62 | 0.28 | 0.10 | 0.00 | 0.35 | 16 | 0 |
| XU.S | 21,936 | 0.43 | 0.14 | 0.01 | 0.04 | 0.25 | 12 | 0 |
| XU.W | 16,488 | 0.13 | 0.02 | 0.04 | 0.04 | 0.03 | 9 | 0 |
| XU.C | 1,680 | 0.00 | — | — | — | — | 1 | 0 |

**Again no array exceeds 20 %** (worst: UW = GNW alone at 12.2 %, and GNW already fails the 2008
site probe at 1.7 dB). One station over 30 %: **YW.FACB 35.3 % (29.6 % H/V, 4.4 % tone)** — the
closest thing in the WA caches to the SHVB failure mode, and YW is a 6-station array of which only
2 were admitted anyway.

Coverage recomputed dropping YW.FACB:

| variant | sub-region | cells | admitted sites | % cells ≥5 sites ≤60 km | % ≤30 km | median sites |
|---|---|---|---|---|---|---|
| as published | Puget-SWA | 146 | 37 | 84.2 | 44.5 | 9 |
| health-adjusted | Puget-SWA | 146 | 36 | **84.2** | 44.5 | 9 |
| as published | S-Olympic | 84 | 37 | 94.0 | 14.3 | 7 |
| health-adjusted | S-Olympic | 84 | 36 | **90.5** | 13.1 | 7 |

So the only coverage number that moves anywhere in Washington is **S-Olympic 2008: 94.0 % → 90.5 %**
(one station, −3.5 points). Everything else is unchanged.

---

## 2. Hour-matched, earthquake-screened quiet set (defect F2)

`freeze_wa_v1` inherits `QUIET_HOURS_PER_DAY = 2` → `linspace(0, 23, 2)` → clock hours **{0, 23}
only**, against candidate hours that span all 24 clock hours, and applies no earthquake screen (the
defect that let one regional event set T = 12.50 on the 2005 bed single-handed).

Replacement rule (`scripts/b_wa_quiet.py`), as Merlin prescribed for 1c′:

1. 40 hours on the 20 quiet days of the window;
2. every clock hour 0–23 appears at least once; the remaining 16 slots are allocated by largest
   remainder to match the **candidate-hour clock histogram**;
3. hours spread over the quiet days (≈2 per day, never two identical (day, hour) slots);
4. **earthquake screen** (catalog query only, no waveform download): ComCat M ≥ 2.5 within 300 km
   of the sub-region centre excludes the hour ± 10 min; M ≥ 5.5 anywhere excludes the origin hour
   and the two following hours (teleseismic body- and surface-wave arrivals);
5. hours where more than half of the bed's cached station-hours are vetoed by §1 are skipped.

### 2.1 N-Olympic 2010 → `results/wa_quiet_hours_2010.csv` (+ `.json`, `_skipped.csv`)

* candidate hours: **454** on 29 tremor days (the same PNSN rule `freeze_wa_v1` used — reproduced,
  not re-derived); clock histogram 13–25 per hour, i.e. genuinely all-day.
* quiet days: the 20 of `freeze_wa_v1` (2010-06-14..06-27, 2010-07-14..07-19).
* earthquakes found in the quiet span: 3 local M ≥ 2.5 within 300 km, 43 global M ≥ 5.5 →
  **113 (day, hour) slots excluded**; 38 placement attempts were rejected (earthquake or no data).
* result: **40 hours, all 24 clock hours present, 1–3 hours on each of the 20 quiet days**, maximum
  health-veto fraction in any selected hour 0.03.
* quota per clock hour: 2 for hours 0–13 and 22–23, 1 for hours 14–21 (the candidate histogram is
  mildly night-heavy).

### 2.2 Puget / S-WA 2008 → `results/wa_quiet_hours_2008.csv`

* candidate hours: **113** hot hours on the 13 array-energy tremor days (there is no tremor
  catalogue before 2010-01-06; the flag and its frozen p95 = 2.9127 are reproduced from
  `results/hour_energy_2008_raw.csv` + `results/probe_days_2008.json`).
* quiet days: the 20 literature-quiet days of `results/hour_energy_flag_2008.csv`
  (2008-02-01..03-31 and 2008-07-05..07-31).
* earthquakes in that 6-month span: 44 local M ≥ 2.5 within 300 km of (46.99, −122.55), 254 global
  M ≥ 5.5 → **728 (day, hour) slots excluded**, 33 placement attempts rejected.
* result: **40 hours, all 24 clock hours present, 1–3 per quiet day**.

⚠ **Finding that v2 must not ignore: the 2008 candidate-hour set is diurnal.** Its clock histogram
is 12–13 hours each at 15–20 UTC and 0–2 at 01–12 UTC — i.e. ~86 % of the "tremor" hours fall in
07–16 local PDT. Tremor is not diurnal; cultural noise is. Two consequences: (i) the hour-matched
quiet set is *essential* for the 2008 pass (the {0, 23} quiet hours of the v1 rule would have
sampled precisely the hours the candidate set does not contain, making the threshold meaningless);
(ii) the 2008 array-energy hour flag itself should be audited (e.g. weekday/weekend and per-station
contrast) before any Puget 2008 discovery run, because it may be selecting daytime anthropogenic
energy rather than tremor. Recorded here, not acted on — no constant was changed.

---

## 3. What `freeze_wa_v2` should take from here

1. `quiet_hours_csv = results/wa_quiet_hours_2010.csv` (and the 2008 file for the Puget pass) —
   read verbatim through `w_common.QUIET_HOURS_CSV` (`W_QUIET_HOURS_CSV`), replacing
   `quiet_hours_per_day = 2`; the WA smoke threshold T = 6.50 is void and must be re-derived on
   this list, with the mandatory leave-one-hour-out column.
2. `health_veto_csv = results/wa_station_hour_veto_2010.csv` (`W_HEALTH_VETO`) as a per-station-hour
   mask in every sum, quiet and tremor alike.
3. The census result itself: **no WA array is >20 % degenerate**, so the coverage figures of
   `P_STEP2_SITES.md` and `freeze_wa_v1` stand as published (recompute logged in
   `results/wa_coverage_health_adjusted_2010.csv`).
4. Still owed by v2 and NOT prepared here (out of this agent's scope): the window-transient veto,
   the NaN/zero-median fix flags, the pair-seed rule, the restriction of the 11-site spread subset
   to the autocorrelation stage only (F5), and `check_constants` against the 1c′ freeze.

Scripts: `scripts/b_wa_health.py`, `scripts/b_wa_quiet.py`, `scripts/b_wa_coverage.py`.
