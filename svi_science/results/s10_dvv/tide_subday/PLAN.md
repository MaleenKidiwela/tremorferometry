# S10 follow-up — tidal-phase dv/v (A) and sub-daily dv/v inside big episodes (B): rules fixed before any number

Written Sun 2026-10-11, ~10:20 am PT, before any tidal or sub-daily stack was measured. Approved by the user the same
day. Base: the S10 pilot (`../pilot/PLAN.md`, Revisions 1–2) — its rows, windows and estimator are reused unchanged.
eps = −dv/v (positive eps = slower). All plotted values are dv/v = −eps. "Scrambled control" = detector rows whose
station data were shifted to another time of the same day; here always the pilot's noise version `c2` (Revision 1: the
c1 rows read at their detector frame time). Times in PT. Order: family 10 first, then 87 and 17, then the pooled-family
result. Nothing below is changed after a number is seen; additions go at the end with time and reason.

## 0. What is reused unchanged (from `scripts/s10_measure.py` and the pilot tables)
- Rows: `../pilot/tables/rows.parquet` (counted-once data repeats with P ≥ 0.5 on active days; era rows carry the era P)
  and the pilot's c1 rows; `row_station.parquet` for station/epoch membership.
- Windows: `../pilot/win/<NET.STA>.npz` (data and site-shifted c1) and `<NET.STA>__c1f.npz` (the c2 noise control):
  2–8 Hz, 40 Hz, despike MAD 8, per-window RMS normalisation, raw (no SVD).
- Per (family, station, response epoch) series, exactly as `run_series`: P-weighted grand stack G0 → S pick → cc
  realign of every data and c2 window (±0.5 s on S−1…S+2 s) → G1 → S re-pick tS. Functions imported from
  `s10_measure.py` (`Ref`, `realign`, `pick_S`, `wmean`, `stretch_trace`), not copied. Z component.
- Estimator: S-anchored stretch of the reference, ±1 % grid, 0.01 % step, parabolic refinement, windows τ = 2–4 s
  (headline) and 3–5 s. |eps| ≥ 0.99 % = no peak in range → not a measurement (pilot Revision 2), dropped from every
  average.
- Stations: the 7 core stations (B011, B926, B009, B010, PGC E1–E3, VGZ E1–E4, SNB E1–E5) and the POLARIS + LZB
  stations for 2003–2006, with the pilot's epochs (KELB and KLNB kept separate). Never reference across an epoch.
- Stack weights: P (P ≥ 0.5). Averages of eps over cells, days, series, stations and families: CC²-weighted, CC = the
  peak correlation of that measurement with its reference.

## A. Tidal-phase dv/v

### A1. Tidal quantities per repeat
- Stress model and hourly tables from S5 unchanged: `results/s5/tidal_stress_bostock_hourly.parquet` for times before
  2014-01-01 and `results/s5_setup/tidal_stress_svi_hourly.parquet` from 2014-01-01 (fid = family number). Continuity
  check before use: in 2014-01 the two tables must agree (max |Δ UDSS| < 1 Pa) for families 10, 87, 17.
- UDSS = `tau_updip_Pa` (body + ocean, + = up-dip shear on the hanging wall = thrust-promoting); σn = `sigma_n_Pa`
  (tension +). Values at the repeat time by linear interpolation of the hourly series; rates = hourly centred difference,
  interpolated the same way (Pa/h).
- **Headline phase θ = M2 phase of UDSS**, computed with `s5_lib` (harmonic fit per era table, θ = 0 at peak positive
  UDSS of M2, θ = −90° at peak rising UDSS rate). Columns: M2 phase of σn.
- Time used = the catalogue row time t (detector time; within seconds of the origin time — 0.01° of M2 phase).
- **Bins.** Phase: 8 bins of 45° centred at 0°, 45°, …, 315°. Value: 8 quantile bins of UDSS per family (edges from that
  family's repeats). Columns with the same binning: UDSS rate, σn, σn rate (quantile bins), M2 σn phase.
- **Fake cycles (control ii):** phase = 2π·frac((t − 2003-01-01 00:00 UTC)/P) for P = 12.0 h and 13.5 h, same 8 bins.
  Flag fixed now: 12.0 h equals the S2 tidal period and the second harmonic of the solar day (time-of-day detection
  effects), so it is not tide-free; a third fake, **P = 10.5 h**, is added as a column for that reason.

### A2. Designs (all three computed; D1 is the headline the pass rule is judged on)
- **D1 — within-day differential (main).** For every active day d of a series and bin b with ≥ 3 data windows: cell
  stack S_db (P-weighted) measured against that day's all-bin stack R_d (P-weighted, all its windows in the series).
  Per bin: eps_b = CC²-weighted mean over days. Every bin draws from the same days, so day-to-day and episode-to-episode
  changes cancel. Column: leave-bin-out day reference (R_d without the cell's windows), observed only.
- **D2 — day-matched pooled.** Per series and bin, the P-weighted stack of all bin-b windows over all days is measured
  against Σ_d W_db·R_d / Σ_d W_db (W_db = ΣP of bin b on day d): the same days with the same day weights. One
  measurement per (series, bin); bins with < 3 windows skipped. Higher signal-to-noise than D1, same day control.
- **D3 — simple pooled (as asked).** Per series and bin, the P-weighted stack of all bin-b windows vs the epoch grand
  stack G1 (all windows, not leave-bin-out; the injection test measures the dilution this causes).
- Aggregation: series → station (all epochs of the station, CC²-weighted per bin) → pooled over all stations → pooled
  over the three families. Sites for the consistency rule: B009, B010, B011 and PGC = one site (PGC); every other
  station is its own site.

### A3. Fits
- Phase variables: eps_b = c + p·cos θ_b + q·sin θ_b, weighted least squares, weights = Σ CC² of bin b. Amplitude
  A = √(p² + q²), phase φ = atan2(q, p). Reported as dv/v: amplitude A and the phase of maximum dv/v (φ + 180°).
- Value variables: weighted linear slope of eps_b on the bin's mean value (kPa or kPa/h), → dv/v per kPa = −slope.

### A4. Controls (pre-registered)
- **(i) random phase shift, 1,000 draws (seed 20261011):** every repeat's θ gets an independent random offset U(0, 360°);
  for value variables every repeat's time is moved by an independent offset U(0, 12.42 h) and the value is read at the
  moved time. Same binning, same design, same fit → distribution of A (phase) or slope (value) per station, pooled and
  pooled-family. Column: a common random offset per (family, day) — keeps each day's within-day clustering.
  p = (1 + #{control ≥ observed}) / 1001.
- **(ii) fake cycles** 12.0 h and 13.5 h (and 10.5 h, column): same designs and fits, each against its own (i)
  distribution. If a fake cycle passes the A5 rule, the tidal result is reported as "not distinguishable from a fake
  cycle".
- **(iii) scrambled control:** the c2 windows (same days, same realignment) through D1 (cells vs the c2 day's all-bin
  stack, ≥ 3 windows) and D3 (c2 bin stacks vs the data grand stack), binned by their own times' tidal phase; their own
  (i) distribution (200 draws).
- **(iv) population check per bin** (per family, all series pooled): mean P; mean within-day log10 S amplitude
  (S-window RMS × stored raw RMS, minus that station-day's mean); mean realignment lag (ms) and realignment cc with G0
  (proxies for sub-location / waveform population; no per-repeat locations exist); S6 matched-filter size dM where the
  repeat has one (2014+). Each gets the same cos fit and a (i) p-value. A population effect is flagged if p < 0.01.
- **Injection (calibration, family 10, all series):** every data window stretched about S by eps_i = a·cos θ_i, a =
  0.1 % (and 0.03 %); all three designs re-run. Recovery = fitted A projected on the injected phase / a. The bound in A6
  is divided by the recovery.

### A5. Pass rule (per family and pooled-family, per design; judged on D1)
All three: (1) pooled A (2–4 s) above the 99th percentile of control (i); (2) at least 2 stations at different sites
whose own vector, projected on the pooled phase, is above the 95th percentile of their own control projection; (3) the
pooled 3–5 s vector projected on the 2–4 s phase is positive. Fake cycles judged by the same rule.

### A6. Otherwise: bound
- 99 % upper bound on the tidal dv/v amplitude: (A_obs + 2.58·σ_c) / recovery, σ_c = SD of the control (i) cos
  component. Per kPa: (|slope| + 2.58·SD of control slopes) / recovery from the UDSS value bins; also the amplitude bound
  divided by the M2 UDSS amplitude at the family.
- Literature comparison from memory, labelled as such: tidal dv/v of ~10⁻⁴–10⁻³ (0.01–0.1 %) measured at shallow depth
  (< 1 km) for tidal strains ~10⁻⁸.

## B. Sub-daily dv/v inside the biggest episodes

### B1. Stacks
- Episodes (r23 `episodes_2003_2026.csv`, start − 2 d … end + 2 d): 2004-07, 2005-09, 2003-02, 2013-09, 2008-05,
  2011-08, 2024-08, 2023-04, 2019-08. Families 10, then 87, 17.
- Bins: 6 h (headline) and 12 h (column), aligned to 00:00 UTC (= 5 pm PDT / 4 pm PST the day before). Only
  (family, series, bin) with ≥ 20 data windows.
- Per bin: P-weighted stack vs the epoch grand stack G1 **leaving that bin's windows out** (the within-epoch reference),
  2–4 s and 3–5 s; railed values dropped. c2 scrambled-control stack of the same bin vs the same reference where it has
  ≥ 3 windows (column).

### B2. Front arrival at the family
- S3 PNSN tremor-front time at the family where it exists (`results/s5/step3_front_arrivals.csv`, the S3b method:
  2019-08, 2023-04, 2024-08). Otherwise the family's own onset: 10th-percentile time of its counted-once high-P repeats
  (P ≥ 0.75; era 0.90) in the episode. The pilot's front date is drawn as well.

### B3. Plot
- One figure per family, one subplot per episode: dv/v (per series minus its CC²-weighted episode mean) vs time in PT;
  points coloured by CC, marker = station; front arrival marked; no error bars.

### B4. Test for a change through the episode
- Per (family, series, episode): Δ_front = CC²-mean(bins after the front) − CC²-mean(bins before), needs ≥ 2 bins on
  each side; Δ_third = mean(last third) − mean(first third) of the time-ordered bins, needs ≥ 3 bins (third size =
  ⌊k/3⌋ ≥ 1).
- Pooled Δ = weighted mean over (series, episode), weight = harmonic mean of the two sides' Σ CC².
- Time-shuffled control: bin time labels permuted within each (series, episode), 2,000 draws, same statistic → two-sided
  p, SD.
- Pass: p < 0.01 in 2–4 s, same sign in 3–5 s, same sign at ≥ 2 stations at different sites. Otherwise bound
  |Δ| + 2.58·SD_shuffle.
- The same for families 87 and 17 and pooled over families.

## Outputs
`tables/` (per-repeat tidal table, cell/bin tables, fits, control distributions, injection, population, sub-daily bins,
tests), `figures/`, `RESULTS.md` (plain words first, then numbers and controls). Scripts `svi_science/scripts/
s10_tide_*.py`; queue `s10_tide_queue.sh` with DONE markers. CPU only, ≤ 6 workers, < 60 GB RAM; the multi-family
scale-up's files are not touched.

## Addendum 1 (2026-10-11, ~10:25 am PT — after the B run, before any post-hoc number)
- B4 before/after-front test as written has **0 usable (series, episode) units**: bins with ≥ 20 repeats almost never
  fall before the front (by construction the family's onset is the 10th-percentile repeat time, and ≥ 20-repeat bins
  come with the activity). The first-vs-last-third test is unaffected (57 units at 6 h, 72 at 12 h, pooled families).
- Post-hoc column (labelled as such, not a pass rule): the same Δ_front with ≥ 1 bin on each side, same shuffle.
