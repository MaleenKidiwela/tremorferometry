# S15 RESULTS — VLFE timing against the local LFE sequence (Q1) and tidal modulation (Q2)

Written 2026-10-17 ~3:05 pm PDT, updated ~3:25 pm after Merlin's audit. Rules: `PLAN.md` (pre-registered ~2:45 pm, before
any number; Amendment 1 ~2:50 pm adds a fairer re-anchored control before its numbers existed; Amendment 2 ~2:55 pm adds
robustness columns for Q2; Amendment 3 ~3:15 pm = Merlin's audit tests, written before their numbers). Log:
`NOTES.md`. Scripts: `scripts/s15_q1.py`, `s15_q1b.py`, `s15_detect.py`, `s15_q2.py`, `s15_q2b.py`, `s15_figs.py`
(each runs in < 20 s). Times PT. "Random-draw control" = the same statistic computed on LFE times drawn from the same
families and episode.

## Plain answer — Q4 (where can we see VLFEs?), slip-rate column, largest events (added 2026-10-17 ~5:45 pm PDT, PLAN Amendment 5)
- **The catalogue sees VLFEs along a ~50 km stretch, not only at C0/C1, and VLFEs occur all along it.** Green's-function
  thrust sources injected through the production detector and battery are recovered at Mw 3.8 as well at four SVI nodes
  48.34–48.52 N (0.42–0.55 in 2012, 0.61–0.76 in 2024) as at the C0/C1 nodes (0.38–0.60 / 0.72–0.77); 100–150 catalogued
  VLFEs have their concurrent LFEs within 16 km of each of those nodes. So "only C0/C1 host strong VLFEs" was a reading of the
  templates' names, not of the data.
- **Outside that stretch the catalogue is blind** (recovery < 0.2, Merlin's "cannot speak" outcome): the NW end of the SVI
  cloud (48.70 N: 0.04 / 0.17) and the southern Olympic/Puget node (47.80 N: 0.00 / 0.01). Under the NE Olympics (48.16 N,
  123.18 W) it is partial (0.04 in 2012, 0.37 in 2024). No node shows recovery ≥ 0.5 with zero real events, so there is no
  sign of physical confinement within what we can see; the Slab2 check is not triggered.
- VLFE locations sit a median 8.8 km NW of their own concurrent LFEs in every episode and both front directions — the known
  location bias, not a real offset.
- Slip-rate column: VLFE days are ordinary slip-rate days for the LFEs around them (median percentile ~50), 0–5 days later
  than those LFEs relative to the regional peak (inside the ±5 d resolution).
- Largest events: every episode tops out at Mw 3.4–3.9, with ≤ ~5 completeness-corrected events ≥ Mw 3.6; big episodes make
  more small VLFEs, not bigger ones. No b-value is measurable.

## Plain answer — Q3 (slip partition) and Q5 (migration rate), added ~4:00 pm PDT (PLAN Amendment 4, Merlin design `MERLIN_PARTITION_MIGRATION_DESIGN.md`)
- **Overlap check (Q3 step 0): the terms are additive.** Re-stacking S13 without LFEs within ±60 s of a catalogued VLFE
  keeps the per-LFE 20–50 s moment at 0.85 (2012-08), 0.93 (2018-06) and 0.92 (2024-08) of the full value (rule ±20 %).
  The VLFE-adjacent LFEs (36–68 per episode) carry 2.6–6.6× more moment each, so 9–17 % of the "continuous" term sits in
  VLFE windows. The partition is LFE : continuous : discrete : GNSS, with that overlap as a stated correction.
- **Partition:** per episode the LFE self-moment is ~1/10 of the continuous slow moment (ratio 8–15 in every episode) and the
  catalogued VLFEs add 0.5–3× the continuous term after ÷ completeness (0.3–2.4× raw); all seismic slow moment together is
  0.1–0.8 % of the GNSS Tier-B moment. So, in units of the LFEs' own moment: LFEs 1, continuous ≈ 10, strong VLFEs ≈ 5–25 (completeness-corrected),
  GNSS Tier B ≈ 4,000–7,000 (main 2018–2024 episodes).
- **Size test: not claimable.** VLFEs per zone LFE (÷ compl_3.4) do not follow episode size (ρ −0.14 vs regional slip,
  n 7), and they track the LFE network's station count (ρ −0.90) — fewer stations, fewer LFEs, more VLFEs per LFE — so by
  the pre-set rule the size test is uninterpretable.
- **Q5, VLFE productivity vs local passage speed: flat.** Against ±24 h matched pseudo-VLFEs, VLFEs per LFE are the same in
  fast, medium and slow local passage (slow/fast 0.91, percentile 14; 4 of 8 episodes one way). Against whole-episode draws
  VLFEs avoid the fastest-passage places (1.59, percentile 100), but that compares places, not timing, and is a column.
- **But VLFEs favour days when the front is advancing fast.** 28 % of VLFEs fall in the top tercile of daily front speed
  (> 10.5 km/d) vs 22 % of matched LFE draws (percentile 100); detection is lower on those days, not higher.
- **The RTR excess was a timing effect.** Against ±24 h matched draws the 2010–13 RTR fraction is 25.6 % vs 24.8 %
  (percentile 56): VLFEs are late in the sequence, RTRs are too, and nothing is left once time is matched.
- **Tides by speed:** declustered VLFE M2 a 0.71 / 0.74 / 1.00 (fast / medium / slow local passage; n 63 / 48 / 53) — a 0.29
  rise, just under the claimable 0.3; the matched LFEs rise the same way (0.43 / 0.50 / 0.63). The Poisson speed × tide
  model shows no trend for VLFEs (A 1.02 / 0.84 / 1.16 by daily speed, SE 0.15–0.19).
- What it means: the 20–50 s slow moment is dominated by a continuous part ~10× the LFEs' own moment, with strong VLFEs
  a comparable discrete add-on; none of it is more than ~1 % of the geodetic slip. Where and how fast the front passes does
  not set how many VLFEs a patch makes; VLFEs ride the times of fast forward advance and the late, tide-sensitive phase.

## Plain answer (updated ~3:25 pm PDT after Merlin's audit, `MERLIN_AUDIT_Q1Q2.md`, and PLAN Amendment 3)
- **Q1 (timing): PASS-late.** VLFEs fire later in each patch's own LFE sequence than a random LFE of the same families
  does. Median time since local LFE onset: **1.79 d for VLFEs vs 1.16 d for the random-draw control** (99th percentile
  of 1,000 control medians 1.35 d; observed above all 1,000). The sign holds in **9 of 9** big episodes and in both
  front directions (NW-ward +0.53 d, SE-ward +0.74 d). Mean rank of a VLFE in its local sequence 0.62 (0.5 = random LFE).
  Defining concurrency with P ≥ 0.5 LFEs gives the same verdict (2.51 vs 1.88 d, p99 2.13, 9/9, both directions).
  - Caveat (Merlin's "would change my mind" item): VLFEs with a single concurrent family in NW-ward fronts (52–59 events,
    ~20 % of the NW-ward set) are **not** late, and stay flat when concurrency is enlarged with P ≥ 0.5 LFEs (−0.2 d,
    percentile 18). The NW-ward lag is carried by the multi-family events.
- **Q2 (tides): stronger than the LFEs.** Headline, declustered (first VLFE per hour): **VLFE M2 a = 0.74** (von Mises
  κ 0.80, peak/trough rate ratio 5.0, n 223) **vs 0.45 for matched LFE draws** (99th percentile 0.64); both peak at
  −42° to −47°, between peak up-dip shear rate (−90°) and peak shear (0°). The honest LFE baseline, with the VLFE's own
  ±1 h removed from the draws, is 0.30 (all events) to 0.36 (declustered). All-event values (a 1.05, κ 1.23, peak/trough
  12, Rayleigh z 141) are inflated by bursts: the tide changes how many VLFEs a burst has, not just whether bursts occur.
  Detection does not follow the M2 phase (injections: b ≤ 0.012).
  - **New finding:** VLFE tidal modulation grows with time since local onset — a 0.91 / 1.11 / 1.37 for Δt < 1 / 1–3 /
    ≥ 3 d, with θ0 moving −53° → −48° → −16° (from peak stress rate toward peak stress). It does not grow with the LFE
    burst around the VLFE (a 1.17 / 1.10 / 0.97 by ±10-min LFE count tercile; 0.92 without LFEs) or with size (1.12 below
    vs 0.98 above Mw_scan 3.26): the tide changes VLFE count, not VLFE size.
  - **Gate test (Poisson regression, hourly, 12 episodes, 495 VLFEs): ACCEPT, but at the threshold.** With the local
    LFE-rate term in, VLFE M2 log-rate amplitude 1.02 (bootstrap 0.74–1.31) vs 0.51 (0.40–0.64) for the zone LFEs:
    ratio **2.01 (1.32–2.78)**; rule ≥ 2 and CI above 1. Without 2016–17 the ratio is 2.00 (intermediate by 0.004).
    VLFE count rises super-linearly with local LFE rate (exponent β 1.36), and β × LFE amplitude = 0.69 lies inside the
    bootstrap range of the VLFE amplitude (difference −0.05 to +0.63), so a convex sampling of slip rate is not excluded
    by this test. Against all pool LFEs, or using only hours with zone LFEs, the VLFE excess is clear (ratio 3.0 and 6.1,
    difference CIs above 0). O1 cross-check: not detected for either (VLFE a 0.15, LFE 0.16, random p99 0.19) — below
    the 0.33 that exponential scaling would predict for VLFEs.
- **What it means physically:** VLFEs are a late-phase, strongly tide-gated part of local slow slip; they are not a
  proportional sample of LFE activity. They come from the same slip (same tidal phase), avoid the first day of local
  activity, and their tidal response strengthens and shifts toward peak stress as the local sequence matures — the same
  path S5 found for LFEs, but about twice as strong. Whether that excess is a separate tidal gate or a strongly convex
  response to slip rate is not settled: the regression gives ~2× and β 1.36 can account for part of it, while the lack
  of any rise with burst intensity argues against a pure convex response.

## Scope (stated once)
The closed catalogue is used as it is (540 ETS VLFEs, 2010–2026; nothing re-vetted). Its limits: per-episode
completeness (C0/C1 injections, about 0.14–0.89 at Mw 3.8 depending on episode) and ±16 km locations (with an ~11.5 km
NW bias). Q1 avoids the location limit by anchoring each VLFE to its own concurrent LFEs; Q2 is insensitive to it.
Q1 answers "when within the local LFE sequence", at the day scale, for VLFEs that have LFEs within ±10 min (422 of the
431 with an eligible control). It does not say where the VLFE source is relative to the front at better than family scale.

## Q1 — timing
Raw figures first: `figures/q1_day_scale.png` (Δt histogram vs control; rank histogram; the 1,000 control medians;
per-episode differences), `figures/q1_variants_direction.png`, `figures/q1_minutes.png`, `figures/q1_detection_vs_local_day.png`.
Tables: `tables/q1_events.csv` (per VLFE, all columns), `q1_day_summary.csv` / `q1b_day_summary.csv`, `q1_day_per_episode.csv` /
`q1b_day_per_episode.csv`, `q1_day_per_direction.csv` / `q1b_day_per_direction.csv`, `q1_other_scales.json`, `q1b_other.json`,
`q1_episode_directions.csv`, `detect_recovery_vs_local_day.csv`.

### Check of the inputs
- The LFE pool (SVI P ≥ 0.75 + WA-NW selected) reproduces the catalogue's `n_lfe_10m` for all 540 events: 431 with
  ≥ 1 LFE within ±10 min, 109 without.
- Front directions: s3 where available; otherwise the sign of the LFE onset progression in our pool (2011-08 NW-ward,
  ρ +0.90; 2012-08 SE-ward, −0.93; 2013-09 NW-ward, +0.75; 2025-09b1 NW-ward, +0.88). Of the 9 big episodes, 7 are NW-ward
  and 2 SE-ward (2012-08, 2018-06).

### Day scale (Merlin's single best first test)
| control | n | VLFE median Δt (d) | control median (d) | control 1st / 99th pct | percentile | big-9 same sign | NW / SE same sign | verdict |
|---|---|---|---|---|---|---|---|---|
| **re-anchored** (Amendment 1; verdict read here) | 422 | **1.79** | 1.16 | 0.99 / 1.35 | 100 | **9/9** | yes / yes | **PASS-late** |
| fixed-weight (pre-registered §3.2) | 431 | 1.77 | 1.15 | 0.96 / 1.36 | 100 | 9/9 | yes / yes | PASS-late |
| fixed-weight, without 2016–17 | 429 | 1.78 | 1.14 | 0.96 / 1.36 | 100 | 9/9 | yes / yes | PASS-late |
| single concurrent family only (no intersection effect) | 65–73 | 1.30–1.46 | 0.93–1.13 | — / 1.37–1.56 | 95 | 7–8/9 | yes / yes | below the 99 % bar |

- The re-anchored control places pseudo-VLFEs at random LFE times of the same families (each with ≥ 1 other LFE within
  ±10 min) and recomputes their concurrent set and onset exactly as for a VLFE. It moved the control median by only 0.01 d,
  so the bias found in the audit (Amendment 1) is small; pseudo-VLFEs have 4.0 concurrent families on average, VLFEs 3.8.
- Per big episode (VLFE median − control median, re-anchored): 2011-08 +2.04 d, 2012-08 +0.30, 2013-09 +1.20, 2018-06 +0.89,
  2019-08 +0.55, 2020-11 +0.62, 2023-04 +0.87, 2024-08 +0.28, 2025-09b1 +1.13 (n 4 with an eligible control; +0.16 with the
  fixed control on n 11). Individually above the 99th percentile: 2011-08, 2018-06, 2023-04; the rest 89–97.
- Rank: mean 0.62, median 0.67; 70 % of VLFEs fall in the later half of their local sequence (figure panel b).
- Robust onset (5th-percentile LFE time instead of the first LFE): 1.50 vs 0.93 d, still above all control medians.
- In plain words: the VLFE comes about 0.6 day (≈ 3 km of front travel at 5 km/day) later in the local sequence than an
  average LFE of the same families. Most of the shift is a shortage of VLFEs in the first day of local activity (panel b).

### Report-back conditions (columns, §0)
- **Detection vs day in the local sequence: flat.** Recovery of C0/C1 injections, days 0–2 vs ≥ 2 after local onset:
  Mw 3.6 0.61 vs 0.57 (+7 %), Mw 3.8 0.74 vs 0.70 (+5 %); big 9 episodes Mw 3.8 0.81 vs 0.80 (+1 %). Not near the 30 %
  line, and where it differs it favours early detection, so it cannot make VLFEs look late.
- **Concurrent-LFE vs VLFE-location variants.** Pooled signs agree (+0.63 vs +0.26 d, location variant at percentile
  98.6). By direction the location variant gives SE-ward +1.39 d and NW-ward −0.18 d. This is the predicted signature of
  the ~11.5 km NW location bias (a NW-placed VLFE meets its neighbourhood's onset early in a NW-ward front and late in a
  SE-ward one); the average of the two directions, where the bias cancels to first order, is +0.61 d — the same as the
  concurrent-LFE answer.
- **Both front directions** show the lag with the concurrent-LFE anchor (above all 1,000 control medians in each).

### The 109 VLFEs without LFEs
- 83 have families with LFEs within 20–30 km in the episode. Location-anchored Δt minus control: NW-ward −0.70 d
  (percentile 0.3), SE-ward +0.53 d (percentile 92). **The same-sign check fails**: their apparent timing flips with front
  direction, as the location bias predicts. The bias-cancelled average is about −0.1 d, i.e. no sign of lateness for these
  events; with ±16 km locations their timing cannot be read at the day scale.

### Minutes scale (LFE counts around the VLFE time, 431 events)
| bin (min) | [−10, −2) | [−2, 0) | [0, +2) | [+2, +10] |
|---|---|---|---|---|
| LFEs | 992 | 346 | 294 | 912 |
| per minute | 124 | 173 | 147 | 114 |

- VLFEs sit in LFE bursts: the rate within ±1 min is ~1.3× the 2–10-min level, peaking within ±20 s of the VLFE time.
  Around random LFEs of the same families the peak is broader and taller (1.9×).
- Slightly more LFEs in the 2 minutes before the VLFE than after: after/before 0.85 (binomial p 0.04; outer bins 0.92,
  p 0.07). Around random LFEs the same ratio is 1.05 (1st–99th 0.85–1.26), so the VLFE value is at the control's 1st
  percentile. In the figure the LFE rate builds over ~1 min into the VLFE and falls to background within ~20 s after it.
  This is a weak (p ≈ 0.04) asymmetry, worth stating as a hint only. The VLFE time is the matched-filter time; the
  LFE–VLFE peak lies within ±20 s of it, so the 2-min bins are not affected by the alignment.

### Hours scale (RTRs 2010–13, Bostock-density list, 30 windows)
- Inside an RTR (time in the window, along-strike within ±20 km): **24 % of VLFEs (31/129) vs 15 % for the random-draw
  control** (99th percentile 22 %; percentile 99.85). With re-anchored pseudo-VLFEs: 25.6 % (30/117) vs 15.9 % (99th 23.9 %).
- Carried by 2012-08: 44 % vs 23 %. 2011-08 12–14 % vs 11 %; 2010-08b1 1 of 8–9 vs 8–11 %; 2013-09 0 of 18 vs 4–5 %.
  So the RTR excess is real in the pooled test but rests mainly on one episode.
- S3 day-scale returns (SVI bursts, 2015–2024): 4 of 289 VLFEs (1.4 %) vs 1.7 % control — no excess. In the two episodes Merlin named:
  2018-06 3 of 83 inside returns (control 4.3 %), 2023-04 0 of 58 (control 1.2 %).

## Q2 — tides
Raw figures first: `figures/q2_tides.png` (M2 and K1 phase histograms relative to random times, VLFEs vs matched LFEs;
rate vs stress decile; per-episode modulation), `figures/q2_detection_vs_tide.png`. Tables: `tables/q2_summary.csv` (all sets ×
τ_updip / τ / Coulomb 0.4 / σn × M2 / K1), `q2_events.csv`, `q2_per_episode.csv`, `q2_stress_deciles.csv`,
`q2_detection_vs_tide.csv`, `q2b_clusters.json`, `q2_info.json`.

### Headline (M2, tidal up-dip shear, S5 phase convention)
| set | n | VLFE a | θ0 | Rayleigh z | matched LFEs a (1st–99th) | LFE θ0 | verdict |
|---|---|---|---|---|---|---|---|
| all with a control | 514 | **1.05** | −42° | 141 | **0.47** (0.36–0.60) | −42° | **stronger** |
| with LFEs | 431 | 1.08 | −39° | 126 | 0.49 (0.36–0.63) | −45° | stronger |
| no LFEs | 83 | 0.92 | −59° | 18 | 0.47 (0.18–0.72) | −22° | stronger |
| without 2016–17 | 505 | 1.05 | −42° | 138 | 0.49 (0.37–0.61) | −42° | stronger |
| declustered (first per 1 h) | 223 | 0.74 | −47° | 31 | 0.45 (0.28–0.64) | −43° | stronger |
| declustered (first per 6 h) | 137 | 0.76 | −59° | — | 0.42 (p99 0.66) | — | stronger |
| 2010–13 | 129 | 1.01 | −41° | 33 | 0.44 (0.18–0.72) | −59° | stronger |
| 2014–26 | 385 | 1.06 | −42° | 108 | 0.50 (0.37–0.64) | −37° | stronger |

- a = 2R is the rate modulation (rate ∝ 1 + a cos(θ − θ0)); a > 1 means the rate nearly vanishes at the opposite phase.
  In the histogram, VLFE rate is 2.3–2.5× the random-time rate at −90° to −30° and 0.1–0.25× at +90° to +180°.
- Cluster bootstrap of the VLFE value (resampling 1-h / 6-h clusters): 95 % range 0.87–1.22 / 0.90–1.20, entirely above the
  LFE 99th percentile 0.60.
- Random-time control (uniform within ±24 h of each VLFE): for all 514, mean a 0.08, 99th percentile 0.21 (0.21–0.46 in the
  smaller subsets); every set is above its 99th percentile.
- Same answer for total shear in the convergence direction (τ: a 1.05, −41°) and Coulomb μ′ 0.4 (a 1.01, −73°); σn alone
  peaks at −132° (clamping side), also 2× the LFEs — the components are correlated.
- **Per big episode** (VLFE a vs matched-LFE median): 2011-08 1.36 vs 0.40, 2012-08 0.91 vs 0.57, 2013-09 0.74 vs 0.45,
  2018-06 1.51 vs 0.53, 2019-08 1.19 vs 0.97, 2020-11 1.08 vs 0.76, 2023-04 1.07 vs 0.40, 2024-08 0.96 vs 0.50,
  2025-09b1 1.80 vs 0.81. VLFEs exceed the LFE median in 9/9 and the LFE 99th percentile in 5/9 (the others at percentile 88–99). Phases agree
  within ~30° in 7/9 (2012-08 and 2023-04 differ by 31° and 32°).
- **Stress level** (total tidal up-dip shear, deciles at random times): VLFE rate 2.1–2.2× random in the top two deciles
  and 0.23–0.31× in the bottom three; matched LFEs 1.5× and 0.52–0.60×. Mean stress at event times: +0.71 σ (VLFE) vs
  +0.42 σ (LFE) vs +0.07 σ (random); with the median family std of ~1.0 kPa this is Δτ ≈ 630 Pa vs 340 Pa, an implied
  exponential sensitivity aσ ≈ 1.6 kPa for VLFEs vs ≈ 2.9 kPa for the LFEs around them (S5's whole-catalogue LFE value with
  τ: 4.7 kPa; the LFEs near VLFEs are themselves late in the episode, where S5 found the tidal response largest).
- **K1**: VLFE a 0.21 vs LFE 0.26 — same. K1 is locked to the solar day, so it is mixed with time-of-day effects (VLFEs
  and LFEs both have a time-of-day signal, R 0.19 and 0.22; injection recovery varies by 2–5 % with K1 phase).

### Checks on the tidal result
- **Detection does not follow M2.** Injected C0/C1 events (17 runs, 27,200 per magnitude): recovery modulation with M2
  phase b = 0.004 (Mw 3.6) and 0.012 (Mw 3.8), 95 % bound ≤ 0.03, against a VLFE modulation of 1.05. Coastal long-period
  noise following the ocean tide therefore cannot produce the VLFE signal.
- The matched LFE draws are within ±24 h of each VLFE (488 events), ±48 h (4) or the whole episode (22); 26 VLFEs had no
  LFE pool and are left out of the comparison (they carry phases in `q2_events.csv`).
- The VLFE stress point is the SVI family nearest the VLFE's concurrent-LFE position (median 3 km away); location error does
  not matter for the phase (Merlin: ≤ 14 % of the shear std across the 10 families nearest C0/C1).

## Amendment 5 — Q4 detection map, slip-rate column, largest events (`scripts/s15_q4_inj.py`, `s15_q4_queue.sh`, `s15_q4_summary.py`, `s15_q4b.py`, `s15_figs3.py`)
Raw figure first: `figures/q4_detection_map.png`. Tables: `tables/q4_recovery.csv` (per node, year, Mw, period, Wilson 95 %),
`data/q4/` (sources, checks, injection lists and scored rows), `q4_centroids.csv`, `q5orig_sliprate.csv`, `q_bsubstitute.csv`.

### Q4 part 1 — recovery of Green's-function thrust sources (300 injections per node × Mw × year; ETS + between, paired times)
| node | where | 2012 Mw 3.6 / 3.8 | 2024 Mw 3.6 / 3.8 | ratio to GF C0/C1 (Mw 3.8, 2012 / 2024) | real VLFEs anchored ≤ 16 km |
|---|---|---|---|---|---|
| G0 (77) | C0 node, 48.16 N 123.58 W | 0.21 / 0.38 | 0.61 / 0.72 | reference | 133 |
| G1 (94) | C1 node, 48.25 N 123.58 W | 0.42 / 0.60 | 0.70 / 0.77 | reference | 141 |
| N1 (112) | SVI, 48.34 N 123.45 W | 0.29 / 0.50 | 0.53 / 0.61 | 1.01 / 0.81 | 101 |
| N2 (109) | SVI, 48.34 N 123.85 W | 0.31 / 0.55 | 0.67 / 0.76 | 1.13 / 1.01 | 108 |
| N3 (126) | SVI, 48.43 N 123.85 W | 0.31 / 0.54 | 0.67 / 0.75 | 1.11 / 1.00 | 150 |
| N4 (143) | SVI, 48.52 N 123.85 W | 0.18 / 0.42 | 0.55 / 0.69 | 0.86 / 0.93 | 110 |
| N5 (176) | SVI NW end, 48.70 N 123.99 W | 0.01 / 0.04 | 0.06 / 0.17 | 0.08 / 0.22 | 0 |
| W1 (14) | S Olympics/Puget, 47.80 N 122.90 W | 0.01 / 0.00 | 0.01 / 0.01 | 0.00 / 0.01 | 3 |
| W2 (80) | NE Olympics, 48.16 N 123.18 W | 0.02 / 0.04 | 0.19 / 0.37 | 0.09 / 0.49 | 68 |
| empirical C0 / C1 (whole year) | — | 0.72, 0.84 / 0.79, 0.89 | 0.59, 0.70 / 0.74, 0.76 | — | — |
- **Check of the GF sources.** On the near vertical channels that dominate the matched filter the GF waveforms match the
  empirical C0/C1 stacks (cc 0.72–0.90 at ±1 s lag); horizontals and distant channels fit poorly (median zero-lag cc over
  all channels 0.15–0.23). In 2012 the GF references are recovered about half as often as the empirical C0/C1 sources
  (0.38–0.60 vs 0.79–0.89 at Mw 3.8); in 2024 they match (0.72–0.77 vs 0.74–0.76). So absolute GF recoveries are a lower
  bound, especially in 2012; the ratios to the GF C0/C1 references are the like-for-like numbers.
- **Reading by Merlin's outcomes.** Visible (recovery ≥ 0.5, ratio ≈ 1): the SVI band 48.15–48.52 N between 123.45 and
  123.85 W — real VLFEs occur throughout it. Blind (< 0.2): 48.70 N on the NW end and 47.80 N in the south. Partial: the
  NE-Olympic node in 2024 (0.37), where 68 VLFEs have concurrent LFEs within 16 km but only 8 are located there (location
  bias pushes them NW). No visible node lacks real events, so nothing here argues for physical confinement; equally, the
  catalogue cannot say whether VLFEs happen beyond the visible band. Finding 13's "a few patches" should read "the band of
  the SVI LFE cloud the network can see, 48.15–48.5 N".
- Only 2012 and 2024 were injected; the visible band will be narrower in years with fewer stations (2015–2018, 2021–22).

### Q4 part 2 — VLFE vs concurrent-LFE along-strike centroids (`tables/q4_centroids.csv`, `q4_centroids_all.json`)
| episode | dir | n | VLFE along (median, km) | concurrent-LFE along | difference (per-event median) |
|---|---|---|---|---|---|
| 2010-08b1 | NW | 8 | −10.1 | −45.5 | +27.8 |
| 2011-08 | NW | 38 | −20.2 | −44.1 | +21.1 |
| 2012-08 | SE | 54 | −21.7 | −32.2 | +9.2 |
| 2013-09 | NW | 18 | +7.1 | −17.9 | +12.3 |
| 2018-06 | SE | 73 | +1.6 | −5.3 | +3.8 |
| 2019-08 | NW | 29 | −35.8 | −41.5 | +7.3 |
| 2020-11 | NW | 23 | −28.0 | −33.7 | +2.9 |
| 2023-04 | NW | 48 | +1.6 | −36.4 | +28.0 |
| 2024-08 | NW | 118 | −7.7 | −11.0 | +4.7 |
| 2025-09b1 | NW | 11 | +7.9 | +8.7 | +2.2 |
- VLFE locations sit a median 8.8 km NW (+ along strike) of their own concurrent LFEs (431 events): 10.1 km in NW-ward
  fronts, 6.2 km in SE-ward fronts. The offset has the same sign in every episode and both front directions, so it is the
  known location bias (~11.5 km NW), not a physical displacement of the VLFE source. Anchored to their LFEs, the VLFEs span
  the same along-strike range as the LFE activity of each episode (10–90 %: −53 to +8 km), moving with the episode
  (−45 km in 2010–11, −5 to −11 km in 2018 and 2024).

### Original Q5 (column) — VLFEs in the regional slip-rate history (`tables/q5orig_sliprate.csv`; resolution ±5 d)
- Percentile of the regional slip rate on VLFE days within its distribution on zone-LFE days, median per episode:
  2010-08b1 45, 2011-08 78, 2012-08 20, 2013-09 42, 2017-03 18, 2018-06 76, 2019-08 64, 2019-11 80, 2020-11 54, 2023-04 52,
  2024-08 34, 2025-09b1 39. Mean over episodes ~52: VLFE days are typical slip-rate days for the LFEs around them.
- Median day offset from the regional slip-rate peak: VLFEs −4 to +3 d in 10 of 12 episodes (2012-08 +10 and 2025-09b1 +9,
  where the peak sits early in the window), zone LFEs −6 to 0 d in the same 10. Per episode, VLFEs fall 0–5 d later than
  the zone LFEs in 11 of 12 (2017-03 1 d earlier), in line with Q1, but this is inside the ±5 d resolution of the slip-rate history. Not a headline.

### b-value substitute (`tables/q_bsubstitute.csv`)
| episode | n | largest Mw_scan | n ≥ 3.6 | compl_3.6 | n ≥ 3.6 ÷ compl (95 %) | n ≥ 3.8 |
|---|---|---|---|---|---|---|
| 2010-08b1 | 9 | 3.93 | 3 | 0.64 | 4.7 (1.0–13.8) | 1 |
| 2011-08 | 41 | 3.53 | 0 | 0.72 | 0 (0–5.1) | 0 |
| 2012-08 | 57 | 3.54 | 0 | 0.88 | 0 (0–4.2) | 0 |
| 2013-09 | 19 | 3.44 | 0 | 0.80 | 0 (0–4.6) | 0 |
| 2018-06 | 83 | 3.61 | 3 | 0.70 | 4.3 (0.9–12.6) | 0 |
| 2019-08 | 31 | 3.55 | 0 | 0.71 | 0 (0–5.2) | 0 |
| 2020-11 | 28 | 3.65 | 1 | 0.75 | 1.3 (0–7.5) | 0 |
| 2023-04 | 58 | 3.77 | 1 | 0.75 | 1.3 (0–7.4) | 0 |
| 2024-08 | 136 | 3.88 | 4 | 0.89 | 4.5 (1.2–11.5) | 1 |
| 2025-09b1 | 21 | 3.44 | 0 | 0.34 | 0 (0–11) | 0 |
- Every episode tops out at Mw 3.4–3.9 (scan Mw, ×2 in moment), and no episode has more than ~5 completeness-corrected
  events ≥ Mw 3.6. The large episodes differ in how many Mw 3.1–3.5 events they make, not in their largest event; with
  0–4 events ≥ 3.6 per episode, a b-value is not measurable.

## Q3 — slip partition (PLAN Amendment 4; `scripts/s15_q3_rows.py`, `s15_q3_restack.sh`, `s15_q3.py`)
Raw figure first: `figures/q3_partition.png`. Tables: `tables/q3_restack_result.csv`, `q3_restack_rows.csv`, `q3_partition.csv`
(all columns, raw and corrected), `q3_size_test.json`. Re-stack outputs: `s13_vlfe/tables/s15ex_svi_*`, `s15in_svi_*`
(same stacker, stations and data as the S13 table run of each episode; pooled unit only).

### Step 0 — overlap re-stack
| episode | full per-LFE M0 | without ±60 s of a VLFE (n) | ratio | only the ±60 s LFEs (n) | ratio to "without" | share of continuous in VLFE windows |
|---|---|---|---|---|---|---|
| 2012-08 | 4.46e12 | 3.77e12 (2,796), z 5.1 | **0.85** | 2.50e13 (68), z 9.5 | 6.6 | 17 % |
| 2018-06 | 6.01e12 | 5.61e12 (1,462), z 8.3 | **0.93** | 2.97e13 (36), z 6.3 | 5.3 | 9 % |
| 2024-08 | 5.01e12 | 4.60e12 (1,166), z 8.0 | **0.92** | 1.18e13 (59), z 4.8 | 2.6 | 13 % |
- All three within ±20 % → **additive** by the pre-set rule. The complement stacks well below the ~1e14 per window a full
  VLFE in every window would give (VLFEs are not exactly aligned with the LFE windows), so the "share in VLFE windows" is the
  overlap that a partition should subtract (column), not the full VLFE moment.

### Per-episode partition (selected columns; N m)
| episode | VLFEs | SVI LFEs | LFE self | continuous | discrete raw | discrete ÷ c3.4 | c3.4 | GNSS Tier B | regional slip (mm) | cont / LFE | disc / cont (corr) | VLFE per 100 zone LFEs (corr) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2010-08b1 | 9 | 1,474 | 8.2e14 | 7.2e15 | 2.6e15 | 5.7e15 | 0.45 | — | 14.1 | 8.7 | 0.80 | 5.8 |
| 2011-08 | 41 | 1,445 | 7.0e14 | 7.7e15 | 4.0e15 | 7.4e15 | 0.53 | — | 8.3 | 10.9 | 0.97 | 20.5 |
| 2012-08 | 57 | 2,864 | 1.6e15 | 1.3e16 | 5.7e15 | 9.4e15 | 0.61 | — | 14.9 | 8.3 | 0.73 | 13.4 |
| 2013-09 | 19 | 1,359 | 7.1e14 | 7.4e15 | 2.3e15 | 3.9e15 | 0.59 | — | 14.9 | 10.3 | 0.53 | 7.4 |
| 2018-06 | 83 (17 relay self) | 1,498 | 8.3e14 | 9.0e15 | 1.05e16 | 1.7e16 | 0.49 | 5.6e18 | 18.5 | 10.8 | 1.89 | 20.5 |
| 2019-08 | 31 | 662 | 3.3e14 | 4.3e15 | 3.6e15 | 5.9e15 | 0.60 | 1.3e18 | 3.7 | 13.0 | 1.39 | 22.3 |
| 2019-11 | 5 | 613 | 3.0e14 | 3.5e15 | 1.1e15 | 2.1e15 | 0.53 | 1.2e18 | 4.0 | 11.7 | 0.61 | 2.6 |
| 2020-11 | 28 | 617 | 3.1e14 | 3.6e15 | 3.4e15 | 7.2e15 | 0.46 | 1.5e18 | 4.5 | 11.8 | 2.01 | 21.6 |
| 2023-04 (flag) | 58 | 881 | 4.9e14 | 4.4e15 | 6.9e15 | 2.2e16 | 0.31 | 4.4e18 | 16.3 | 9.1 | 5.0 | 53.7 |
| 2024-08 | 136 | 1,225 | 6.8e14 | 6.1e15 | 1.47e16 | 1.86e16 | 0.79 | 4.7e18 | 15.4 | 9.1 | 3.03 | 47.5 |
| 2025-09b1 (flag) | 21 | 170 | 1.2e14 | 6.3e14 | 2.3e15 | 1.1e16 | 0.21 | 4.3e18 | 13.8 | 5.5 | 17 | 223 |
| 2015-12 (flag) | 3 | 1,138 | 5.9e14 | 6.2e15 | 5.2e14 | 1.1e16 | 0.05 | 5.3e18 | 15.8 | 10.5 | 1.77 | 19.8 |
| 2021-01 | 3 | 394 | 2.3e14 | 2.8e15 | 5.8e14 | 1.4e15 | 0.43 | 3.8e18 | 11.2 | 12.0 | 0.50 | 3.7 |
| 2021-10 | 1 | 344 | 1.6e14 | 2.5e15 | 1.5e14 | 2.1e15 | 0.07 | 1.3e18 | 1.0 | 15.1 | 0.83 | 13.3 |
| 2022-03 | 0 | 381 | 2.2e14 | 3.1e15 | 0 | 0 | 0.14 | 1.5e18 | 4.7 | 14.0 | 0 | — (0 zone LFEs) |
- LFE self-moment uses each family's `mag_med` (5.6e11 N m per LFE on average). Continuous / LFE self is 8–15 in every
  episode (S13's ~10–12× accompaniment). Discrete ÷ compl_3.4 is an under-correction (the median VLFE, Mw 3.26, is below
  the lowest injected magnitude); the Mw ≥ 3.5 bin-corrected floor is only 0–25 % of it (`disc_m0_ge35_bincorr`).
- Seismic slow moment (LFE + continuous + corrected discrete) is 0.1–0.8 % of GNSS Tier-B moment (2014+). Tier B is
  regional (×~1.4 low, 34 km half-potency radius), so this compares a patch-scale seismic sum with island-scale slip.
- Cross-era trends in disc/cont (0.5–1.0 in 2010–13 vs 1.4–3.0 in 2018–24) are not admissible: VLFE completeness and the
  LFE catalogue (backward vs forward, 15 vs 13 stations/day) both change between eras.

### Size test (matched-completeness subset: 2011-08, 2012-08, 2013-09, 2018-06, 2019-08, 2020-11, 2024-08)
| quantity | ρ vs regional slip | ρ vs Tier A D (n 4) | ρ vs compl_3.4 | ρ vs n_sta_day | max / min |
|---|---|---|---|---|---|
| VLFE per 100 zone LFEs ÷ compl_3.4 | −0.14 | −0.40 | 0.18 | **−0.90** | 6.5 |
| VLFE per 100 zone LFEs, raw | −0.07 | 0.00 | 0.32 | −0.70 | 8.7 |
| disc / cont ÷ compl_3.4 | 0.11 | 0.40 | −0.07 | −0.68 | 5.8 |
- **Verdict: uninterpretable** by the pre-set rule (|ρ| with n_sta_day ≥ 0.6). The size correlation itself is flat. The
  station-count link is the LFE denominator: years with 13 LFE stations per day (2019–24) have fewer LFEs per patch than
  years with 15, so "VLFEs per LFE" rises. A per-LFE productivity across eras therefore needs an LFE completeness model the
  project does not yet have (`s0/PLAN_S0.md` forbids the Pdet column as a correction).

## Q5 — migration rate (PLAN Amendment 4; `scripts/s15_q5.py`)
Raw figure first: `figures/q5_migration.png` (pace distributions of VLFEs, matched pseudo-VLFEs and all LFEs; VLFE-per-LFE by
tercile; tidal a by tercile). Tables: `tables/q5_summary.json`, `q5_per_episode.csv`, `q5_episode_speeds.csv`, `q5_events.csv`,
`q5_tide_by_tercile.csv`, `q5_poisson_speed_tide.csv`, `q5_detection_by_speed.csv`.

- Local passage pace = Theil–Sen of family onsets (S3 rule) vs along-strike within ±10 km, folded to the front direction.
  Tercile edges from all pool LFEs: pace ≤ 0.055 d/km (faster than 18 km/d, or simultaneous/reversed; 21 % of LFEs and
  18 % of VLFEs have pace ≤ 0) / 0.055–0.129 / > 0.129 d/km (slower than 7.8 km/d). Daily front speed terciles: < 5.3,
  5.3–10.5, > 10.5 km/d.
- Per-episode main speeds from our pool agree with s3 in 6 of 8 shared episodes; 2019-08 and 2023-04 disagree in sign (as in
  Q1's direction table), which makes their folded local paces small — a caveat for those two.

| test | VLFE shares fast / med / slow | control shares | slow ÷ fast (per LFE) | control 1st–99th | percentile | verdict |
|---|---|---|---|---|---|---|
| **local pace vs ±24 h matched pseudo-VLFEs (n 419)** | 0.31 / 0.36 / 0.33 | 0.29 / 0.37 / 0.34 | **0.91** | 0.82–1.23 | 14 | **flat** (4 of 8 episodes; NW 0.92, SE 0.70) |
| local pace vs whole-episode draws (column) | 0.31 / 0.36 / 0.33 | 0.41 / 0.31 / 0.28 | 1.59 | 0.76–1.30 | 100 | VLFEs avoid fast-passage places |
| daily front speed vs ±24 h draws (n 360) | 0.28 / 0.37 / 0.35 | 0.22 / 0.39 / 0.36 | 0.75 | 0.82–1.22 | 0 | more VLFEs on fast-advance days |
- Fast-forward class (top daily-speed tercile): 28.3 % of VLFEs vs 22.3 % of matched draws (99th percentile 25.8 %).
- RTRs 2010–13 against ±24 h matched draws: 25.6 % vs 24.8 % (99th 31.6 %, percentile 56). Against whole-episode draws it
  was 25.6 % vs 15.9 % (Amendment 1). So the RTR excess comes from VLFE timing (late), not from reversals as such.
- VLFE size does not change with speed: median Mw_scan 3.28 / 3.24 / 3.26 (fast / medium / slow).
- Tides: declustered VLFE M2 a 0.71 / 0.74 / 1.00 vs matched LFEs 0.43 / 0.50 / 0.63 (fast / medium / slow local pace;
  ±0.15 per tercile). The 0.29 VLFE rise is just under the claimable 0.3 and the LFEs rise alike, so it is the S5
  maturity effect, not a separate speed control. Poisson model (hourly, daily-speed tercile × tide, onset bins and LFE rate
  in): VLFE A 1.02 / 0.84 / 1.16 (slowest / medium / fastest; SE 0.16 / 0.15 / 0.19), LFE A 0.43 / 0.47 / 0.59 — no trend.
- Detection vs speed (column, 30 % flag): recovery fast / slow local pace 0.57 / 0.63 (Mw 3.6), 0.69 / 0.78 (Mw 3.8),
  outer difference 11–12 %; slow / fast daily speed 0.69 / 0.55 (Mw 3.6), 0.79 / 0.68 (Mw 3.8), −21 % / −14 %. No 30 % flag.
  Recovery is lower on fast-advance days, so the fast-day VLFE excess is not a detection effect (it is, if anything,
  understated by ~15–20 %).

## Amendment 3 — Merlin audit tests (run ~3:15–3:22 pm PDT; `scripts/s15_q2c.py`, P ≥ 0.5 runs via `S15_PMIN=0.5`)
Tables: `tables/q2c_gate.csv`, `q2c_columns.json`, `q2c_a_by_bins.csv`, `q2c_episodes.csv`, `q2c_hourly.parquet`,
`tables/p050/` (all Q1 tables with P ≥ 0.5 concurrency).

### Q1 columns
| subset | n | VLFE median Δt | control median | control p99 | percentile |
|---|---|---|---|---|---|
| single concurrent family, NW-ward (P ≥ 0.75) | 52 | 1.15 d | 1.10 d | 1.58 d | 60 |
| single concurrent family, SE-ward (P ≥ 0.75) | 13 | 2.35 d | 1.19 d | 2.21 d | 99.5 |
| same NW-ward events, concurrency from P ≥ 0.5 (now 2.3 families each) | 59 | 1.42 d | 1.62 d | 2.12 d | 18 |
| same SE-ward events, P ≥ 0.5 | 14 | 4.03 d | 3.36 d | 5.22 d | 82 |
| all with LFEs, P ≥ 0.5 (re-anchored) | 431 | 2.51 d | 1.88 d | 2.13 d | 100 (9/9 episodes; NW +0.49, SE +0.87 d) |

- The full P ≥ 0.5 test passes the same rule; onsets move earlier with the looser catalogue, so all Δt grow.
- The NW-ward single-family events show no lag at either threshold. Merlin flagged exactly this as the result that would
  reopen a direction-dependent anchor effect, so the NW-ward part of the lag should be read as a property of the
  multi-family (busier) VLFEs, not of all VLFEs.

### Q2 columns
- κ (von Mises) and peak/trough rate ratio exp(2κ): all events κ 1.23, ratio 12; declustered κ 0.80, ratio 5.0.
- Matched-LFE baseline without the VLFE's own ±1 h: a 0.30 (p99 0.41) for all events, 0.36 (p99 0.54) declustered —
  against 0.47 / 0.44 with it. VLFE values 1.05 / 0.74 are above both.
- VLFE a by bins (`q2c_a_by_bins.csv`):

  | by | bins | a | θ0 |
  |---|---|---|---|
  | days since local onset | < 1 / 1–3 / ≥ 3 d (n 125 / 197 / 109) | 0.91 / 1.11 / 1.37 | −53° / −48° / −16° |
  | ±10-min LFE count | ≤ 3 / 4–6 / ≥ 7 (n 179 / 117 / 135); no LFEs (83) | 1.17 / 1.10 / 0.97; 0.92 | −32° / −51° / −40°; −59° |
  | size | Mw_scan < 3.26 / ≥ 3.26 (n 263 / 251) | 1.12 / 0.98 | −46° / −37° |

### Gate: Poisson regression of hourly counts (PLAN Amendment 3.1)
12 episodes with ≥ 5 VLFEs (495 VLFEs, 9,000 hours); zone = pool families within 20 km of each episode's VLFE centroid
(19–45 families); onset bins pre / 0–1 / 1–3 / 3–6 / ≥ 6 d; episode offsets; 500 day-block bootstrap resamples.

| model | VLFE A (95 %) | θ0 | β | LFE A (95 %) | LFE θ0 | ratio (95 %) | β·A_LFE | A_V − β·A_LFE (95 %) | outcome |
|---|---|---|---|---|---|---|---|---|---|
| **headline: zone-LFE rate term** | **1.02 (0.74–1.31)** | −34° | 1.36 | **0.51 (0.40–0.64)** | −63° | **2.01 (1.32–2.78)** | 0.69 | −0.05 to +0.63 | **ACCEPT (at threshold)** |
| no LFE-rate term | 1.25 (1.00–1.55) | −42° | — | 0.51 | −63° | 2.47 (1.85–3.34) | — | — | ACCEPT |
| all pool LFEs as rate term and LFE model | 0.97 (0.75–1.27) | −44° | 1.25 | 0.32 (0.26–0.38) | −39° | 3.05 (2.15–4.27) | 0.40 | +0.31 to +0.90 | ACCEPT |
| hours with ≥ 1 zone LFE (1,016 h) | 1.19 (0.87–1.65) | −45° | 0.63 | 0.20 (0.12–0.29) | −75° | 6.09 (3.48–10.6) | 0.12 | +0.70 to +1.52 | ACCEPT |
| without 2016–17 episodes | 1.01 (0.75–1.28) | −34° | 1.36 | 0.51 (0.39–0.65) | −63° | 2.00 (1.33–2.85) | 0.69 | −0.04 to +0.66 | intermediate |
| O1 (zone rate term) | 0.14 (0.04–0.45) | −30° | 1.46 | 0.06 (0.02–0.21) | −115° | 2.4 (0.55–13.7) | 0.09 | −0.14 to +0.35 | uninformative (low power) |

- A is the amplitude in log rate (rate ∝ exp(A cos(θ − θ0))), comparable to κ, not to a.
- Reading: the pre-set acceptance is met by 0.008 in ratio. The headline model cannot separate "~2× more tide-sensitive"
  from "convex response to a modulated slip rate" (β·A_LFE inside the range). The variants where the LFE rate is better
  measured relative to its own tide (all pool LFEs; non-empty hours) favour a genuine VLFE excess. The zone-LFE rate is
  sparse (median 0 LFEs per hour; 11 % of hours non-empty), which limits how much of the tide the rate term can absorb.
- O1 event-level cross-check (all events): VLFE a 0.15 (θ0 −95°), matched LFEs 0.16 (p99 0.27), random p99 0.19 — not
  detected; declustered 0.25 vs LFE 0.13 (random p99 0.30). Exponential scaling from M2 predicts ~0.33 for VLFEs, so O1
  does not confirm the M2 excess (it does not contradict it strongly either, given random p99 0.19–0.30).
- K1 (Merlin's note): "same" is uninterpretable because the solar-day signal (R 0.19) is as large as the K1 signal; the
  exponential model would predict VLFE K1 a ≈ 0.46 against 0.21 observed.

## What the two results say together
- VLFEs are rare at the very start of local LFE activity and favour the following days, when S5 found LFEs most tide-sensitive
  and when Bostock's RTRs run (the RTR excess in 2012-08 fits this). Even against LFEs of the same families within ±24 h,
  VLFEs are about twice as tidally modulated, so the tidal result is not just a consequence of the late timing.
- Both point the same way: VLFEs are a late-phase, strongly tide-gated part of local slow slip; they are not a
  proportional sample of LFE activity. A "VLFE = random LFE" picture (Gomberg 2016 sum of LFEs with VLFE times drawn at
  random from LFE times) does not fit at the day scale or in tidal sensitivity. The tie between Q1 and Q2 is direct: VLFE
  tidal modulation rises with days since local onset (a 0.91 → 1.37, θ0 −53° → −16°), as S5 found for LFEs.
- Superseded by Merlin's audit: the earlier phrase "large excursions of an already-modulated slip rate" is not supported as
  stated — VLFE tidal concentration does not rise with burst intensity and the densest LFE bursts are the least tidal.
  The regression leaves a strongly convex response to slip rate (β 1.36) as a partial alternative.
- Not shown here: whether VLFEs are the sum of LFEs in moment (Q3), and anything about VLFEs outside the C0/C1/Strait zone
  the catalogue can see.

## Files
- `PLAN.md`, `NOTES.md`, `README.md`, `MERLIN_SCIENCE_PLAN.md`.
- `MERLIN_AUDIT_Q1Q2.md` (audit); `scripts/s15_q2c.py` (Amendment 3: gate, O1, columns).
- Amendment 5: `scripts/s15_q4_inj.py`, `s15_q4_queue.sh`, `s15_q4_summary.py`, `s15_q4b.py`, `s15_figs3.py`; `data/q4/`; figure `q4_detection_map.png`.
- `MERLIN_PARTITION_MIGRATION_DESIGN.md`; `scripts/s15_q3_rows.py`, `s15_q3_restack.sh`, `s15_q3.py` (Q3), `s15_q5.py` (Q5), `s15_figs2.py`;
  figures `q3_partition.png`, `q5_migration.png`.
- `scripts/s15_common.py` (loaders; `S15_PMIN` env for the P ≥ 0.5 column), `s15_q1.py` (day/minutes/hours, fixed control), `s15_q1b.py` (re-anchored control),
  `s15_detect.py` (recovery vs local day), `s15_q2.py` (tides), `s15_q2b.py` (tidal checks), `s15_figs.py`.
- `tables/`: as listed in each section; draw arrays `q1_*.npy`, `q1b_*.npy`, `q2_phases.npz`.
- `figures/`: `q1_day_scale.png`, `q1_variants_direction.png`, `q1_minutes.png`, `q1_detection_vs_local_day.png`,
  `q2_tides.png`, `q2_detection_vs_tide.png`.
