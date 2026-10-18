# S12-LVL, LFE side — PLAN (rules written before any measurement)

Written Tue 2026-10-13, ~10:20 am PT, after the forward calibration (task 1) and before any task-2..5 number.
Goal: find, in our LFE data, the ETS-related Vs increase inside the low-velocity layer (LVL) that Gosselin et al.
(2020, Sci. Adv. 6:eaay5174) see in receiver functions at CN.PGC, SNB and VGZ (+0.05–0.1 km/s across the LVL top in the
15 days after the ETS midpoint, then decaying; layer ~3.5 km). The user is clear the signal exists: every observable below
gets a forward-calibrated target, and a non-detection is reported as (measured bound) / (target) with the next adaptation.

## 1. Forward calibration (task 1) — DONE, numbers in RESULTS.md §1
Model: Vs inside a 3.5 km layer whose top passes through the LFE hypocentre, dipping 10.2° toward N52E (plane fit to
Slab2 under the SVI families; families sit 3.7 km above Slab2 on median), LVL Vs 3.0 km/s, change **δ = +2 %** nominal
(range 1–3 %). Script `scripts/s12_forward.py` (imports `s11_kernel_mc.py`, not edited). Tables `tables/fwd_*`.

| observable | target for δ = +2 % (median pair; range over pairs / kernels / δ 1–3 %) |
|---|---|
| coda dv/v, S-anchored 3C stretch, 2–5 s | **+0.35 to +0.44 %** (kernels ell20/ell40/two3_20; pairs q10–q90 0.15–0.85 %) |
| coda dv/v per lapse window 1–2 / 2–3 / 3–4 / 4–6 s | +0.28–0.38 / +0.34–0.45 / +0.32–0.43 / +0.34–0.40 % |
| LVL-base S reflection (SrS) lag after S | lag 1.4–2.3 s (median 1.8 s); **change −59 ms** (−30 to −90 ms) |
| oceanic-Moho S reflection (7 km below top) | lag 1.5–4.1 s (median 3.0 s); change ≈ −59 ms (same LVL legs) |
| LVL-base PrP after P (Z) | lag 0.05–1.0 s; change −1.7 s × δVp/Vp |
| P/S amplitude ratio, source inside the LVL | +5 % (Vp fixed) / +2.5 % (δVp = δVs/2) / 0 (Vp/Vs constant or source above the LVL) |
| S−P of direct waves | ≈ 0 (direct rays leave upward; LVL not crossed) |
| coda decay | no velocity target; sensitivity d lnE(τ) = −ω (time in LVL) δQ⁻¹, time in LVL ≈ 0.2 τ (0.47 s at 2–2.5 s, 1.1 s at 6 s) |

Key consequence (written before measuring): a Gosselin-sized change is **large** for the S-anchored coda (~0.4 %), ~25×
the S11 round-7 within-episode bound (lag 0–1 vs ≥ 3: < 0.016 %, 2σ). So, if the change is local to the families, it
cannot switch on between lag 0–1 and lag ≥ 3 days. The open places are (i) the contrast between ETS LFEs and
**between-ETS** LFEs (S11's ΔD compares episodes with neighbours and cancels any change shared by all majors) and
(ii) pre-front LFE days, (iii) phases that see the layer directly (the LVL-base reflection).

## 2. LVL-base reflection search (task 2) — `scripts/s12_reflect.py`
Data: S11 grand stacks (P-weighted, 3C, 2–8 Hz, store + store_ext), every causal series (causal25_3c > 1.5), epochs of
one (family, station) combined. Observables: residual log-envelope (horizontal energy for SrS, vertical for SrP; about P
for PrP) after a quadratic coda trend fit 0.6–6.5 s; and the S-wavelet matched filter |cc| on the horizontals.
Moveout stack over all pairs for a grid of layer models H 1.5–6 km × Vs 2.6–3.8 km/s (flat-layer ray theory).
**Scrambled control:** predicted lags permuted among the pairs of the same station (keeps each station's own
reverberations and lag distribution), 1000 draws, maximum over the model grid (scan-corrected p).
**Detection rule (all of):** scan-corrected p < 0.01; the best model's peak > 0 at ≥ 2 of the 4 sites (PGC site, B926,
SNB, VGZ); same sign in odd- and even-numbered families. Record sections (envelopes sorted by distance with predicted
curves) are drawn whatever the result.
If detected: delay of the phase vs direct S per epoch class (task-3 classes; within-ETS pre/post front, post-ETS vs
between-ETS, per episode where n allows); detection of a change = |Δ| ≥ 3σ, front-scrambled p < 0.01, sign negative
(faster). Target −59 ms. Independently of detection, the windowed SrS shift (±0.3 s around the nominal predicted lag,
lag search ±0.15 s) is carried as a column in task 3 (its target depends on the reflection's share of the window energy;
not a headline).

## 3. Superposed-epoch coda dv/v and decay (task 3) — `scripts/s12_epoch.py`
Data: S11 round-3 cc-weighted day stacks (store_ccw, + _ext in SVI), causal series; fronts = S11 round-7 headline front
dates (S3 onset, else own detection) per (family, major). Not duplicating round 7: it measured 2–5 s per day
(lag 0–1 vs ≥ 3). Here: lapse-resolved windows, decay, and the **between-ETS baseline**.
Classes per family: pre (−10..−1 d), 0, 1, 2, 3, 4–6, 7–10, 11–20, tail 21–60 d after the family's front; **betw** =
> 60 d after the previous and > 30 d before the next major front; Gosselin frame: −15..−1, 0..15, 16..45 d from the
major's regional midpoint (median front). Aggregate **0..20** = front lags 0–20 d (classes 0 through 11–20).
Per series: class stack (all majors superposed) vs the series' grand stack; 1–2, 2–3, 3–4, 4–6, 2–5 s joint 3C stretch
(bank ±3 %); contrast = class − betw. Decay: per-day noise-corrected ln E per window and slope 1–6 s, class = weighted
day mean. Pool: inverse-variance over series (series σ from 200 bootstrap draws over majors and inter-ETS intervals);
pooled σ from the same joint draws for all series.
Controls: front-scrambled (front of each (family, major) replaced by start(major) + rel(family, other major)), 200 draws;
time-of-day scrambled control stacks (Sc) through the same classes.
**Detection rule (headline = 0..20 − betw, 2–5 s and per window):** |contrast| ≥ 3σ; positive (speed-up); same sign at
≥ 3 of 4 sites; time-of-day control within 2σ of 0; lapse pattern not against the forward one (no window of opposite
sign at ≥ 2σ). For pre vs 0..20: same rule plus front-scrambled p < 0.01. Any detection is checked for the seasonal
cycle (classes' calendar months) and for source composition (bursts vs majors: the 'betw' days are bursts) before it
is called a velocity change.
Non-detection: 2σ bound / forward target, and the next adaptation.
Order: SVI first, then WA-NW (same code, S11_REGION=WA_NW; WA has no LVL forward run — SVI kernels are used as the
target scale, labelled).

## 4. Source-side Vp/Vs proxy (task 4) — same script, same classes and controls
ln(A_P/A_S) on Z and on H (projection on the reference wavelets in P −0.15..+0.6 s and S −0.15..+0.6 s) and δ(S−P).
Target: +5 % / +2.5 % / 0 (table above). Rule: as §3. S−P is a geometry check (target ≈ 0); S11 ballistic round 3 already
bounds it (arrival vs between +0.008 ± 0.016 %).

## 5. Shear-wave splitting (task 5) — only if time allows; direct S does not cross the LVL (target ≈ 0) unless the
reflection is found, so splitting would be measured on the reflected phase.

## Common
- Resources ≤ 12 CPU, < 80 GB. Long jobs in `scripts/s12_queue.sh` (nohup, DONE markers), block 15 in
  `scripts/resume_after_reboot.sh`.
- Plots: line graphs, dv/v points coloured by weight / CC, no error bars on dv/v.
- Any rule change after seeing data is logged in NOTES.md as an amendment with time and reason.

## 6. Coordinator tasks A–D (added Tue 2026-10-13 ~10:25 am PT, before any A–D number; user: "we are trying to actually
find something new" — follow the fault valve in time and space with sources inside the fault zone)
Measurement `scripts/s12_novel_measure.py` (one pass per region over store_ccw), analysis `scripts/s12_novel.py`.
Observables: coda dv/v 2–5 s (best precision) and 4–6 s (largest LVL share relative to the 25–45 km layer, 0.45 vs 0.30
at 1–2 s), per-day ln-energy decay slope (Qc proxy), and for C the source scalars. Forward target for coda: a_LVL × δ,
a_LVL from §1 (SVI series-specific; WA uses the SVI median, labelled).

**A. Propagating onset.** Daily front-lag classes D−10..D20 (all majors superposed, per series, joint bootstrap). Model:
dv/v(lag) = c + A·f(lag), f = share of the family's LVL coda sensitivity already behind the front, with the front moving
at v = 10 km/d and the along-strike sensitivity Gaussian with σ = 7 km (f(0) = 0.5, f(1 d) = 0.92; v = 5 km/d column).
A is the LVL-change signal in dv/v units; target A = a_LVL × 2 % ≈ +0.45 % (2–5 s). Also A with the regional-midpoint
alignment (non-propagating alternative; Gosselin frame classes). Control: front-scrambled daily classes (200 draws), A
refitted per draw. Detection: |A| ≥ 3σ (bootstrap), front-scrambled p < 0.01, A > 0, A(local) > A(regional).
**B. Recharge between ETS.** Day-level dv/v (2–5, 4–6 s) and decay slope of between-ETS days (and all days ≥ 60 d after
a front) regressed within series on time since the previous major front (per year), with an annual cycle (cos, sin) and
series offsets; pooled WLS, cluster bootstrap over inter-ETS intervals (200). Prediction (fault valve): dv/v falls as Pf
rebuilds (negative slope). Target: if the Gosselin post-ETS +2 % decays to the inter-ETS level, the inter-ETS recharge
is a fall of a_LVL × 2 % ≈ 0.45 % over the cycle (~14 mo SVI) → slope ≈ −0.4 %/yr if linear, or concentrated in the first
months (exponential, τ 30–60 d). Detection: |slope| ≥ 3σ, negative, same sign in ≥ 3 of 4 SVI sites / both WA halves
(odd/even families), annual cycle included; column without the annual cycle.
**C. Source properties across the front.** Per-day S-window spectral ratio ln(E 5–8 / E 2–4 Hz), S envelope width,
ln(P/S) Z (RMS), within-series contrasts by front lag (pre, 0, 1, 2–3, 4–20) vs between-ETS, adjusted for day stack
size (ln events per day, ln weight) by within-series regression. S6 per-detection size across the front is quoted (S6
test 1). Expected size: a Pf drop of 1–10 MPa against an effective normal stress of a few MPa could change the stress drop
by tens of % → corner frequency by ~10–30 % (∝ Δσ^(1/3)) → spectral ratio change of order 0.1–0.3 in ln units (labelled
order-of-magnitude, not a calibrated target). Detection: |Δ| ≥ 3σ (joint bootstrap), same sign at ≥ 3 sites, not
explained by the size adjustment.
**D. Scaling with episode size.** Per major: pooled (over series) contrast of its 0..20 stack vs the between-ETS stack;
WLS on episode duration Du (headline, round-6 sizes), Fs and log n̄ (stack length) as columns. Fault valve: bigger ETS →
more drainage → larger positive dv/v. Target: 0.45 % × relative size spread (if the change scales linearly with
duration, the slope per 10 d ≈ 0.45 % × 10 d / mean Du). Detection: slope ≥ 3σ (bootstrap over majors), positive, robust
to log n̄.

## 7. Confirmation of the WA-NW recharge candidate (written Tue 2026-10-13 ~10:30 am PT, before any §7 number)
Candidate (§6 B + post-hoc checks): WA-NW between-ETS dv/v 2–5 s falls −0.33 ± 0.14 %/yr with time since the previous
front (11 families with round-7 fronts; scrambled p 0.004). Scripts `scripts/s12_b_measure.py` (per-day measurement of
EVERY store_ccw series; families without a round-7 front take the major's median front), `scripts/s12_recharge2.py`.
Model as §6 B: within-series, cc²-weighted, per-site annual cycle; σ = cluster bootstrap over inter-ETS intervals (200).

**Forward lapse patterns (computed before the data, from §1 kernels, slope per window relative to 2–5 s):** a change in
the LVL is flat with lapse (a_LVL 1–2/2–3/3–4/4–6 s = 0.19/0.22/0.22/0.20, ell40; ratio 4–6 / 1–2 ≈ 1.05–1.2); a shallow
0–5 km change grows (a_s 0.17/0.26/0.30/0.32, ratio 1.6–2.1); a 25–45 km change falls (a_d 0.65/0.54/0.51/0.45, 0.7).
So "the slope grows with later lapse" would point to a SHALLOW cause, not the LVL. The test fits slope_w = k·pattern_w for
each pattern and reports χ² per pattern; it is descriptive (classifies depth), not a pass/fail item.

**Confirmed only if all of C1–C4 hold:**
- C1 families: the slope is negative at ≥ 2σ in the families NOT among the 11 of §6 (independent set), OR ≥ 65 % of the
  families with ≥ 15 between-ETS day-measurements have a negative slope (binomial p < 0.05).
- C2 source / family mix: (a) the fit restricted to families with between-ETS days both < 6 and ≥ 6 months after the ETS
  keeps the same sign and ≥ 50 % of the slope; (b) adding the day's S−P (P-window lag after S alignment) as a covariate
  keeps ≥ 50 % of the slope; (c) the S−P slope with time since ETS and the family-position (along-strike, lat/lon)
  drift are reported.
- C3 stations: same sign in boreholes (PB) and in surface stations, each fitted alone.
- C4 season: with the annual cycle per STATION plus a semi-annual term the slope stays ≤ −2σ, and a season-matched
  scrambled control (time since ETS shuffled within series and calendar quarter, 500 draws) gives p < 0.01.
Descriptive: per-station slopes; lapse pattern (above); shape — exponential A·exp(−t/τ) (τ grid 20–540 d) vs linear,
ΔAIC and τ with an interval bootstrap; days 21–60 d after the front added as a column (they lie on the early curve).
**Replication.** WA-SW: counted before measuring (catalogue_v2, region SW, ΣP ≥ 3 family-days): 0 between-ETS family-days
at P ≥ 0.5, 0.6 or 0.75 → **no replication possible with existing data** (logged). SVI by era: the §6 B fit split
2003–2013 (backward catalogue) / 2014–2026. A replication = negative slope at ≥ 2σ in that era.
**Tie-in.** Per inter-ETS interval, the early level (60–180 d, within-series residual after season) against the size of
the preceding major (Du headline; T, Fs columns; round-6 sizes): a positive slope at ≥ 2σ is "consistent with the
round-7 size hint".
**1–3 Hz reflection restack** (`scripts/s12_lowband.py`): SVI near stations, per (family, station) up to 150 highest-P
data repeats re-read from the existing raw day files, band 1–3 Hz, cc realignment on the S pulse, stacked; then the §2
scan, scrambled control and injection unchanged (envelope smoothing scaled to the band). Detection rule as §2. Target:
reflection ~1.8 s after S, −59 ms for +2 %; a 1 km gradational base still reflects at 1–3 Hz (λ_S ≈ 1.5 km).
