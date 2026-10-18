# Q4 results — does Washington-NW show SVI's slow deep rebuild after a big ETS?

Run 2026-10-07, ~2:29–2:44 pm PT. The rules were frozen first in `PLAN.md`. Scripts:
- `svi_science/scripts/q4_recharge_wa.py` (build + run, ~13 min CPU)
- `q4_core9.py` (fixed 9-borehole column)
- `q4_posthoc.py` (NOT pre-registered, descriptive)
- `q4_figures.py`

Logs are in `q4.log`. The verdict is in `q4_summary.json`. There was no Merlin review, so the result is provisional.
"Time-reversed control" = templates reversed in time. "Scrambled control" = station data shifted to other times of day.

## Plain answer
- **No. Washington does not confirm the slow deep rebuild.** The pre-registered verdict is **NOT SUPPORTED**: in E3
  (2015–26), the time-shift control gives p = 0.09 for the first-month deficit (the rule needs < 0.01, and ≥ 0.05 fails).
- **What does repeat: the first month after a big ETS is almost empty of deep LFEs.**
  - This holds in every era: 0.1× the later rate in E1, E2 and E3.
  - It also holds in PNSN deep tremor, which does not depend on our detector.
- **Why the LFE test still fails:** Washington's deep activity between ETS comes in a few huge bursts.
  - Four days hold 48 % of all E3 events: 2018-11-19 (185 events), 2017-11-11 (177), 2021-04-17 (107) and
    2022-04-23 (65). All fall 154–228 d after a big ETS.
  - So shifting days within a cycle easily makes a deficit this big (median 0.45, 99th percentile 0.95, observed 0.90).
  - Counting days with any deep LFE instead of events (not pre-registered), the WA first-month deficit is gone:
    E3 0.24, p 0.36.
- **"Rises for more than a year" is not seen in E3.**
  - The E3 rate peaks 150–240 d after the ETS (29 per 1,000 family-days), then falls at 240–400 d (4.5).
  - Only E1 and E2 keep rising past 240 d, with 5 and 7 cycles.
- **Time scale: it does not match SVI in any way that means something.**
  - The LFE fits overlap SVI's range: WA E3 τ = 358 d (140–400), all eras 245 d (112–400), SVI 85–400. But each upper
    end sits at the 400-d fit limit, and the WA curve is a burst bump rather than a smooth rise.
  - The detection-independent PNSN deep tremor recovers in about 1.5 months: τ = 48 d (31–103), a much shorter time
    scale than a year.
- **SVI re-check:** SVI's own pass survives the stricter definition (all majors dropped):
  - first-month deficit p = 0.007;
  - τ = 171 d (46–400).
  - Counting days instead of events, SVI's deficit stays (p 0.007), but the recovery shortens to τ ≈ 53 d (23–204).
  - So in both regions, the "> 1 year" part rests on event volume in bursts. The robust part is a quiet first month or two.

## Definitions (frozen, `PLAN.md`)
- **Events:** core-net NW, `P_own_WAVE_PNET ≥ 0.75`, counted once by dedup variant D. There are 26,818 counted-once
  data events from 2007-01-01 to 2026-05-27.
- **Deep:** position-based. The deepest third of A/B-located families by Slab2 depth at the primary epicentre
  (≥ 44.5 km): 39 families, 27 source groups. Down-dip-third and incl.-C sets are columns.
- **Activity:** all deep LFEs on days outside every one of the 69 NW majors (minor-episode and local-tremor-day events
  included). Majors' days are dropped from both counts and exposure.
- **Exposure:** scanned family-days from `rates_daily_v2_dedup`.
- **Clocks:** the headline clock is the 31 NW big ETS of `b_ets_record` (primary list). Columns use all 69 majors and
  the a2 list (19). t = days since the end, out to 400 d.
- **Ends:** D1 regional end (headline); D2 end + 7 d; D3 the family's own last tremor day within 30 km.
- **Eras:** E1 2007–11, E2 2012–14, E3 2015–2026-05-27. Each era is tested on its own. "All eras" = each era's exposure
  scaled by its own mean rate.
- **Note found before running (from P3's code):** P3's big-clock row dropped only the big episodes, so deep events
  inside non-big majors counted at later t. Q4 drops all majors. The P3-as-coded version is a column, and SVI was
  re-run both ways.

## Numbers (big-ETS clock, D1, headline deep set; `recharge_fits_wa.csv`)
Rates are per 1,000 family-days. "Control 99th" = 99th percentile of the time-shift control (each cycle rotated, 1,000×).

| set | events | cycles | 0–30 d / 31–400 d | D30 (95 % range) | control median / 99th | p (D30) | β (p) | τ d (95 % range) |
|---|---|---|---|---|---|---|---|---|
| **E3 2015–26 (confirmatory)** | 1,099 | 21 | 0.96 / 9.9 | **0.90** (0.80–0.94) | 0.45 / 0.95 | **0.090** | 1.03 (0.076) | 358 (140–400*) |
| E1 2007–11 | 886 | 5 | 1.5 / 16.2 | 0.90 (0.76–0.99) | 0.23 / 0.97 | 0.073 | 0.54 (0.040) | 100 (49–233) |
| E2 2012–14 | 373 | 7 | 0.5 / 13.4 | 0.96 (0.59–0.99) | 0.89 / 0.99 | 0.198 | 1.91 (0.039) | 400* (111–400) |
| all eras, normalised per era | 2,358 | 31 | 0.10 / 1.15 × mean | 0.92 (0.86–0.95) | 0.36 / 0.89 | 0.005 | 0.87 (0.005) | 245 (112–400*) |
| PNSN deep tremor 2009-08…2026-05 (per day) | 8,772 | 29 | 531 / 2,000 | 0.73 (0.62–0.84) | 0.16 / 0.56 | 0.001 | 0.36 (0.001) | **48 (31–103)** |
| PNSN deep tremor E3 only | 5,777 | 21 | 566 / 1,959 | 0.71 | 0.15 / 0.60 | 0.001 | 0.36 (0.001) | 48 (32–81) |
| time-reversed control (all eras) | 7 | 31 | — | — | — | — | — | — |
| scrambled control (all eras) | 15 | 31 | — | — | — | — | — | — |

\* at the fit's 400-d limit.

**Rule check:**
- E3: D30 > 0 but p = 0.090, and β > 0 but p = 0.076 → fails.
- E1/E2 same sign: yes.
- PNSN same sign: yes.

Since p(D30) ≥ 0.05, the verdict is **NOT SUPPORTED**. Both controls are < 1 % of the data.

**Binned shape, rate relative to its own mean** (figure panel a):

| t (d) | 0–7 | 7–14 | 14–30 | 30–60 | 60–90 | 90–150 | 150–240 | 240–400 |
|---|---|---|---|---|---|---|---|---|
| E1 | 0.05 | 0.10 | 0.13 | 0.13 | 0.86 | 1.07 | 1.21 | 1.40 |
| E2 | 0.11 | 0.05 | 0 | 0.17 | 0 | 0.31 | 1.87 | 2.50 |
| E3 | 0.04 | 0.05 | 0.19 | 0.08 | 0.19 | 0.49 | **3.44** | 0.54 |
| PNSN deep tremor | 0.04 | 0.59 | 0.29 | 0.57 | 1.11 | 1.42 | 1.19 | 1.09 |
| SVI (all majors dropped) | 0.07 | 0.04 | 0 | 0.32 | 0.72 | 0.91 | 1.37 | 1.46 |

## Columns (none changes the verdict; E3 p for D30)
- **Episode ends:**
  - D2 end + 7 d: 0.066 (all eras 0.003).
  - D3 family's own end: D30 0.89, τ 400* (no control, as in P3).
- **Deep sets and cut-off:**
  - Down-dip third: 0.093.
  - Incl. C families: 0.085.
  - P ≥ 0.85: 0.095 (all eras 0.002).
- **Counting:**
  - Source-group once: 0.076.
  - Raw: 0.096.
  - Live-core days only: 0.090.
- **Clocks:**
  - a2 list (one per year): 0.187.
  - P3-as-coded (non-big majors' events kept): 0.071, all eras 0.002.
  - **All 69 majors: E3 D30 −0.03 (p 0.50); all eras 0.20 (p 0.21). No deficit at all** (panel b), as in SVI's
    39-major clock.
- **Fixed 9-borehole run** (`recharge_fits_core9.csv`): its catalogue finished at ~2:25 pm PT, so it was used.
  - It is WAVE-only and raw-counted, because no dedup flags exist for it.
  - E3: D30 0.35 (p 0.12). E1 0.73 (p 0.08). E2 0.40 (p 0.20). All eras 0.42 (p 0.013).
  - Same direction, weaker.
  - Its controls at 0.75 are large: E3 usable 3,633 data vs 990 time-reversed vs 2,001 scrambled. So it is a noisy
    column.
- **SVI re-check (P3 data):**
  - As coded, it reproduces P3: D30 0.685, p 0.002, τ 322 (86–400).
  - With all majors dropped: D30 0.977, p 0.007, τ 171 (46–400).
  - Halves: 2014–19 p 0.051, 2020–26 p 0.034.
- **Post-hoc, NOT pre-registered** (`recharge_fits_posthoc_daycount.csv`, panel d): each day counted active or not.

  | region | D30 | p | τ (d) |
  |---|---|---|---|
  | WA E1 | 0.23 | 0.45 | — |
  | WA E2 | 0.56 | 0.36 | — |
  | WA E3 | 0.24 | 0.36 | — |
  | WA all eras | 0.30 | 0.15 | — |
  | SVI | 0.86 | 0.007 | 53 (23–204) |

  The busiest E3 days are in `busiest_days_E3.csv`.

## What this means
- **Washington agrees with SVI on one thing:** right after a big ETS, the deep zone goes quiet.
  - Deep LFEs are at ~0.1× for the first 1–2 months in all three eras.
  - Deep PNSN tremor is quiet for ~2 weeks and below normal for ~2 months.
- **It does not show a slow, year-long rebuild.**
  - In E3, deep LFE output peaks 5–8 months after a big ETS, in a few large bursts, and then drops.
  - The tremor recovers with τ ≈ 1.5 months.
- **The clock matters.** Measured from any major ETS, there is no deficit at all. This is the same as in SVI.
- **Caveats:**
  - Big-ETS cycles in WA are short: 31 in 19 years, and only E1/E2 have many gaps over 240 d.
  - The 400-d window is mostly unsampled in E3.
  - Bursts make every event-count test low-powered. A day-count or burst-count rule would be the fairer pre-registered
    test next time.

## Files
- **Plan and verdict:** `PLAN.md`, `q4_summary.json`
- **Tables:**
  - `recharge_fits_wa.csv`, `recharge_curves_wa.csv` (every set, all eras)
  - `recharge_fits_core9.csv`, `recharge_curves_core9.csv`
  - `recharge_fits_posthoc_daycount.csv`, `recharge_curves_posthoc_daycount.csv`
  - `busiest_days_E3.csv`
  - `events_p075.parquet` (event cache)
- **Logs:** `q4.log`, `q4_run.out`, `q4_core9.out`
- **Figure:** `svi_science/figures/q4_recharge_curves.png`. Panels:
  - a) WA big clock by era + all eras + SVI (as coded and all majors dropped)
  - b) all-majors clock
  - c) PNSN deep tremor
  - d) post-hoc day count, WA vs SVI
