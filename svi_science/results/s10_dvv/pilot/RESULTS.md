# S10 pilot — daily active-day LFE stacks for coda dv/v: results in plain words

Written Sun 2026-10-11, ~9:15 am PT. Rules: `PLAN.md` (written first; Revisions 1–2 added before any gate verdict was
read). Scripts `svi_science/scripts/s10_*.py` (queue `s10_queue.sh`). Tables `tables/`, figures `figures/`.
ε = −dv/v (positive ε = slower). "Scrambled control" = the detector run with each station site's data shifted to
another time of the same day. Everything here is provisional: one pilot, three families.

## What was run
- Families 10, 87, 17; active days (Σ P of counted-once P ≥ 0.5 repeats ≥ 3): 120, 108, 96 days, 2003-01 → 2026-09.
  8,444 data repeats and 2,949 scrambled-control rows on those days.
- Stations: B011, B926, B009, B010, PGC (3 response epochs), VGZ (4), SNB (5), and for 2003–2006 the 12 POLARIS
  stations + LZB. 83 (family, station, epoch) series with ≥ 3 days; 48 of them at the 7 core stations.
- Windows: 2014+ repeats from the cut store (all 2,665 matched); everything else from the raw day files at the same
  offsets. The two readers agree to cc ≥ 0.99998 on 297 window-components (`tables/verify_cut_vs_raw.csv`).
- One stack per (family, station, day), P-weighted, Z component; S re-picked and every window cc-realigned on the
  epoch grand stack; stretch about S against the leave-day-out epoch stack in 2–4 s and 3–5 s; 200-draw bootstrap
  error per day; pairwise doublet graph per series; per-episode stacks; step injection; precision vs n.

## Answers

**1. Do daily active-day stacks give usable dv/v?** Not as a daily series. A typical day stack holds 14–17 repeats and
its ε has a bootstrap error of **0.39 %** (median over 1,844 core-station days; best stations PGC and B010 0.30 %,
B011 0.34 %, B926/B009 0.45 %, VGZ/SNB 0.5 %). 19 % of the core days (4–9 % at PGC/B010/B011, 27–37 % at
B926/SNB/VGZ) give no correlation peak inside ±1 % at all and are dropped as non-measurements. Expected
signals (≤ 0.05–0.1 %) are 4–8 × below one day's error.

**2. At what precision?** Measured with disjoint half-stacks in the busiest days (`figures/s10_precision_vs_n.png`,
`tables/precision.csv`):
- PB boreholes (2–4 s): σ ≈ 0.65 / 0.44 / 0.29 / 0.20 % at n = 10 / 20 / 50 / 100 repeats (σ ∝ n^−0.5). That gives
  **n ≈ 400 repeats for 0.1 % and ≈ 1,600–1,800 for 0.05 %** (3–5 s: ≈ 200–350 and 500–1,200).
- Surface broadband (PGC/VGZ/SNB) and POLARIS: σ ≈ 0.37 % at n = 100 and 0.27 % at n = 200, falling only as
  n^−0.3 (the small-n end is capped by the ±1 % search). Extrapolated, 0.1 % needs thousands of repeats; this is
  not reliable beyond the measured range.
- Episode stacks (all repeats of one episode) at the 5 best stations: median n 37, σ 0.26 %; 0.17 % for n
  100–300 and 0.10 % for n > 300. Only the largest episodes reach 0.1 %.
- So 0.05 % per family-station is out of reach on any one day or episode; it needs pooling over many families,
  stations or episodes, and that pooling would have to be shown to be honest.

**3. Do the honesty gates pass?** No series passes all five (0 of 48 core, 0 of 83 total; `tables/gates.csv`).
- (a) Causality (coda / pre-P noise RMS > 1.5): **passes** 44 of 47 core series (coda 2–15 × the noise). The
  exceptions are VGZ for family 87 in three epochs (0.5–1.4): no usable coda there.
- (b) Scrambled control: the noise control (Revision 1: same scrambled rows read at the detector frame time) is
  consistent with zero in 15 of 23 core series that have ≥ 3 control days. The mirror rule (|mirror| < |coda| on
  ≥ 70 % of days) **fails everywhere** (median 15 %). The reason is that the pre-P window holds no coherent
  signal, so its "stretch" lands at random inside ±1 %. The rule as written cannot pass here, and failing it does
  not show a coda artefact. The site-shifted control from the plan is not a noise control at the PGC site. The
  scrambled detections are carried by that one site (station z ≈ 6 at B011/PGC vs ≈ 0.5 elsewhere), so its
  windows there are real LFEs (control-stack cc with the coda reference 0.65–0.69 at PGC/B010/B011, 0.13–0.17 at
  B926/VGZ/SNB).
- (c) P ≥ 0.75-only vs P-weighted within errors: passes 11 of 48 core series (median 84 % of days within 2σ,
  bar 90 %).
- (d) Lapse test (3–5 s should follow 2–4 s): **fails** 43 of 48 (median slope 0.40, sign agreement 43 %, bar 0.6 /
  70 %). The day-to-day ε differences behave as window-specific noise, not as a velocity change.
- (e) Step injection ±0.1 % / ±0.3 % into the second half: passes 28 of 42 core series. All four boreholes pass
  except B010 family 17 (leakage 0.36). PGC E1/E3 also pass. The failures are short or noisy SNB/VGZ/PGC-E2
  epochs (daily recovery 0.6–1.35, leakage up to 0.9). The graph estimator recovers better than the daily one
  (39 of 42).
- (f) Metadata overlay: drawn on every panel (`figures/s10_daily_fam*.png`). Each epoch is its own zero, so no
  instrument step can enter a series. KELB → KLNB passed the gain-normalised check for family 87 only (cc 0.79 vs
  criterion 0.73, Δε −0.10 ± 0.25 %, coda/S ratio 1.06) and was merged there (`tables/daily_kelbm.parquet`). For
  family 17 it did not merge (stretch at the ±1 % edge). Family 10 has no KELB-b days.

**4. Is there any change above the errors and the scrambled control?** Nothing that meets the PLAN §7 bar
(> 3σ, outside the control, reproduced at a second station or window).
- Within an episode (`figures/s10_within_episode.png`, `tables/within_episode.csv`): ε minus the episode mean
  against days since the front reached the family. Pooled over 1,386 station-days, the slope is
  −0.0004 ± 0.0037 %/day and day 0–3 minus later is 0.000 %. The largest bin offsets (family 87 days 3–6
  −0.13 %, family 17 day 3 +0.15 %) are ≤ 1.5σ.
- Long-term (`tables/longterm.csv`): per-series trends are 0.00 ± 0.01 %/yr on the boreholes over 16–20 yr. 7 of 39
  core series have |z| > 2 (≈ 2 expected), but they disagree between families at the same station: B011 −0.018 %/yr
  for families 10 and 17, +0.003 for 87; B926 −0.025 for 87, 0.00 for 10 and 17. That is not one medium change.
  Episode-level trends: 4 of 37 beyond 2σ.
- Same-day agreement between stations (`tables/core_agreement.csv`) is mostly within errors (χ²_red ≈ 0.4–1.8),
  with weak correlations. The one repeated positive pair is PGC–B010 (ρ 0.44–0.69 in 4 family-epoch pairs), but
  B011, which sits at PGC itself, shows none. So it is not a medium signal; a likely cause is that every station
  stacks the same LFE subset each day.

## POLARIS addition (2003–2006)
- Active days with a stack per station (families 10 / 17 / 87): MGCB 14/20/16, TWBB 14/20/16, SILB 14/21/16,
  SSIB 14/20/16, TWKB 14/20/16, TSJB 14/20/0, KLNB 14/16/11, LZB 20/22/0, SHDB 5/11/11, SHVB 8/9/10, KHVB 8/9/0,
  TWGB 5/13/0, KELB 6/5/5. The family has 20 / 23 / 20 active days in 2003–2006.
- **They add no new days.** Every one of those days already has a PGC (and mostly VGZ/SNB) stack. What they add is
  more stations on the same days.
- Precision is like the CN stations (median day σ 0.2–0.8 %, median 0.42 %); causality passes on 32 of 35 series
  (SHVB family 87 fails; TWBB and TWGB family 17 have S−P too long for the mirror window).
- Agreement with PGC/VGZ/SNB on the same days: 90 of 95 station pairs agree within errors (χ² p > 0.01). The
  weighted correlations are scattered around zero (−0.9 … +0.8), and the points do not lie on the 1:1 line
  (`figures/s10_polaris_vs_pgc.png`). They agree because both sides are noise of the size of their error bars, not
  because they track a common change.

## Verdict (provisional)
- The P-weighted daily unit works mechanically: coda is causal and the estimator recovers injected steps on the
  boreholes and PGC. But one day's stack is a ~0.3–0.5 % measurement. The window-lapse test fails, so day-to-day ε
  is noise, not velocity.
- This is a precision limit. It is not a sign of a broken method.
- Usable precision needs ~400 (2–4 s) or ~200–350 (3–5 s) repeats per stack on a borehole for 0.1 %, i.e.
  big-episode stacks or pooling across families. That is the regime where Merlin's per-episode design was aimed.
- No within-episode or long-term change was found above the errors and the scrambled control. The lapse,
  P-subset and step-injection gates also bound what these data can say.
- Not done in this pilot: horizontals (N/E ε are stored as columns but not gated) and other families.

## SVD (added Sun 2026-10-11, ~9:45 am PT; rules `svd/SVD.md`, written before any number)
Plain words first. The original pipeline's SVD-Wiener filter (from `scripts/dvv_roll30cal.py`, unchanged) makes the
day stacks look cleaner and halves the error bars. But it does this by **shrinking every measured change, real or
injected, to about half its size**. A known step we inject comes back at ~0.45 of its true size, and about 0.4 of
it leaks into the days before the step. In units of a true velocity change the SVD daily series is therefore
*less* precise than the raw one, not more. This is a provisional pilot result for three families and five stations.

- **Does SVD shrink the per-day error? Yes on paper.** Median day σ (2–4 s) 0.35 % raw → 0.21 % SVD (ratio 0.56 per
  day; 3–5 s the same, 0.55). The σ includes the filter: the SVD was re-run inside every bootstrap draw. A joint
  bootstrap (all days at once) gives the same ratio (0.55). The independent half-split check agrees with these σ
  (empirical / predicted 0.84 SVD, 0.90 raw; no series above 1.11). cc with the reference rises from 0.71 to 0.93,
  and days with no peak inside ±1 % drop from 13 % to 4 %. Effective rank: median 12 singular values above the
  noise level (1–17), out of a median 71 days per series; Σ of the Wiener weights 21.
- **Does it leak? Yes, in every series.** Pre-registered step test (+0.1 % and +0.3 %, three step dates between
  episodes per series, both windows): SVD recovers 0.10–0.70 of the step (median 0.44), raw 0.87–1.50 (median
  1.02). Measured against a fixed reference, 0.33–0.50 of the step shows up in the days *before* it (raw: 0 by
  construction). **0 of 21 SVD series pass; 18 of 21 raw series pass** (the raw fails are B010 fam 10 and two short PGC
  E2 series). A one-episode or one-day transient of +0.3 % also comes back at 0.45–0.47 (raw 1.02–1.03). The raw
  and SVD ε of the same days correlate at 0.90 with slope 0.46. The filter scales the measurement down; it does not
  remove noise that is separate from the signal.
- **Net effect.** Error in units of a true change (σ / step recovery, a post-hoc column) is median 0.35 % raw and
  0.48 % SVD: worse in all 21 series (×1.42). By the rule fixed in SVD.md ("shrinks" needs σ ratio < 0.8 AND
  half-split ≤ 1.5 AND injection pass), **0 of 21 series** have an SVD gain.
- **Does the lapse test now pass? No.** It passes in 8 of 21 series with SVD vs 7 of 21 raw (median 3–5 s vs 2–4 s
  slope 0.54 vs 0.45, sign agreement 73 % vs 66 %). The filter shrinks both windows by the same factor, so it cannot
  turn noise into a consistent change with lapse time.
- **Scrambled control (noise windows) through the same path:** consistent with zero about as often as raw. Pass
  12/16 raw, 11/16 SVD own-matrix, 12/18 SVD inserted into the data matrix. Most control days still give no
  peak inside ±1 % (median ~70 % of days), so these counts rest on a few days per series. The control cc with the reference rises
  from 0.13 raw to 0.21 when the noise stack is filtered inside the data matrix. The filter pulls noise slightly
  toward the reference shape, but nowhere near the data cc (0.93).
- **Is it the one parameter?** Checked post hoc. Estimating the noise level from the pre-P noise window instead of
  the original's [−4, −2) s before S (which holds P coda here) lowers the noise level ~15 % and raises the rank, but
  recovery only goes from 0.48 to 0.51 and leakage from 0.34 to 0.32 (20 series, 0 pass). The shrinkage is not set by
  that choice.
- **Why it differs from the earlier "SVD innocent" benchmark (likely, not tested here):** the original filter ran on
  30-day rolling stacks. Each of those rows already averages ≥ 5 days, and neighbouring rows share most of their
  content. One active-day stack of ~15 repeats is far noisier, so the weak stretch pattern falls under the Wiener
  noise level and gets damped.
- **Reading.** For the daily unit the filter buys nothing usable: the smaller error bars are paid for by an
  equal-or-larger loss of signal, plus leakage across the step. The raw daily stacks (or bigger stacks: episodes,
  pooled families) remain the honest path. Not tried: SVD on episode stacks or rolling stacks of this record, other
  Wiener weights, horizontals.
- Figures: `svd/figures/s10_svd_stacks_fam{10,17,87}.png` (packed day rows, raw | SVD, B011 and PGC; the raw side
  here is realigned, unlike `figures/s10_daily_stacks_fam10.png`), `svd/figures/s10_svd_daily_fam{10,17,87}.png`
  (daily dv/v = −ε vs calendar time, raw and SVD, ± σ, ETS shading). Tables `svd/tables/`.
