# C — LFE activity at a patch after the slip front passes (SVI + Washington-NW): results

Rules were frozen first in `PLAN_C.md` (2026-10-07 ~1:25 pm PT). Run 2026-10-07 ~1:25–1:35 pm PT. Scripts:
`svi_science/scripts/c_decay.py` (pre-registered), `c_decay_post.py` and `c_decay_check.py` (checks added afterwards, marked
**post-hoc**), `c_figures.py`. Headline = P ≥ 0.75, counted once. SVI 10 big episodes (7,753 LFEs), NW 14 big episodes
(9,295 LFEs). Patch = SVI group (12) or NW 10-km along-strike bin × down-dip third. Front arrival = PNSN tremor
(10th-percentile time per 5-km bin), independent of our detections.

## 1. Plain-words answer

- **A patch switches on when the front arrives, stays active for about 2–3 days, and then stops. It does not
  decay like aftershocks.** In single patch-episodes the "plateau, then a sharp stop" model wins in most cases where the data
  can tell (SVI 23 of 31, NW 17 of 28, SVI fixed stations 13 of 16). Plateau length is 2.6–3.0 days after the busiest
  half-day. Stacked over many patches, the curve *looks* Omori-like or exponential, but that comes from stacking: patches
  with different on-times are added together.
- **Duration does not clearly depend on depth.** In SVI the deepest patches are shorter (2.3 vs 2.8–3.0 days; fixed
  stations 1.3 vs 3.2 days), but this fails the pre-registered test (p 0.09; fixed stations p 0.12). In Washington-NW
  there is no trend (p 0.80). Duration does not follow episode size (tremor count or GNSS offset). **Post-hoc:** it does
  follow front speed. In Washington-NW, slower fronts keep a patch active longer (ρ −0.91, p 0.0001). Duration × speed
  is about the same in every episode: an active band roughly **15 km** wide (interquartile range 10–18 km; SVI 13 km)
  travelling with the front.
- **Inside those days, LFEs come in bursts lasting hours, not at random.** The coefficient of variation of the time
  between events is 2.5–2.8. Poisson would give 1.0, and shuffling events within each day gives 1.6–1.7. About 70 % of
  events sit in bursts, against 25 % after shuffling within each day. This holds on the fixed-station run, after merging
  events less than 60 s apart, and is much weaker in the scrambled control (stations shifted to different days). Most of
  the extra clustering is at the scale of hours: shuffling within each hour keeps nearly all of it. Burstiness does not
  depend on depth.
- **Burst sizes are neither a clean power law nor Poisson-like.** Bursts are much larger than shuffling allows: 16–22 % of
  bursts have ≥ 5 events, against < 1 % after shuffling. The bulk fits an exponential better than a power law (SVI
  p 0.0009; NW undecided, exponential after merging within 10 s), with a somewhat heavy tail above ~10 events.
- **Re-activation later in the same episode is rare (about 6 % of patch-episodes) and not tied to depth.** When it
  happens, several patches restart on the same day, i.e. a second front, not a local restart. SVI 2018-06 restarted on
  Jul 8 PT (groups 9 and 11; the back-running front from S3). Washington-NW 2021-09 restarted on Oct 20–24 PT
  (4 patches, PNSN tremor confirms).

Figures: `figures/c_rate_vs_front.png` (stacked rate vs time since front per depth class; single-patch shape; duration
vs down-dip), `figures/c_iet.png` (inter-event times vs controls, incl. the scrambled control), `figures/c_burst_sizes.png`.

## 2. Numbers with controls

### Q1 shape (Rule S) — FAIL as a stack shape; cell answer = plateau-then-stop
| run | stack winner (ΔAIC) | decided cells: plateau / Omori / exp | undecided | plateau T (median) | real plateaus (T > 2w, T > 0.5 d) |
|---|---|---|---|---|---|
| SVI headline | Omori (103) | 23 / 6 / 2 | 21 of 52 | 2.6 d, w 0.02 d | 18 of 23 |
| SVI fixed-S13 | Omori (68) | 13 / 3 / 0 | 14 of 30 | 2.8 d | 12 of 13 |
| NW headline | exponential (119) | 17 / 9 / 2 | 28 of 56 | 3.0 d, w 0.11 d | 15 of 17 |

The stack winner differs from the cell winner, so under Rule S the stack shape is a product of stacking. Omori p values
of the few Omori-winning cells are 1.0–4.8, mostly > 2: steeper than aftershock decay. Stacked lags: median 2.7 d (SVI),
1.8 d (NW) after tremor arrival, 90 % of events within 7.7 / 8.4 days. 4 % (SVI) and 13 % (NW) of LFEs come before the
local tremor arrival; in SVI deep families this rises to 10 % (cf. S3b, not a general pattern). The median single-patch
profile (busiest PT day = 1) goes 0.05–0.5 the day before, then 0.2–0.6 the day after → 0.1–0.2 → ≈ 0 by day +4 (fig. panel b). Day
binning and a busiest-day anchor make this look more peaked than the plateau fits.

### Q1 duration vs depth (Rule D) — FAIL in both regions
| run | cells | D80 (30-event thinned) shallow / mid / deep, days | within-episode ρ (dd) | p | tail T90 ρ, p (column) |
|---|---|---|---|---|---|
| SVI headline | 72 | 2.77 / 2.83 / 2.27 | −0.23 | 0.089 | −0.35, 0.010 |
| SVI P ≥ 0.5 | 90 | 3.10 / 2.98 / 2.26 | −0.22 | 0.054 | −0.30, 0.008 |
| SVI fixed-S13 | 37 | 3.21 / 2.54 / 1.34 | −0.31 | 0.12 | −0.43, 0.023 |
| NW headline | 115 | 3.21 / 2.83 / 2.90 | +0.03 | 0.80 | −0.05, 0.66 |
| post-hoc single sources: SVI families | 65 | 2.50 / 2.92 / 1.72 (6 deep) | −0.17 | 0.28 | |
| post-hoc single sources: NW sources | 85 | 2.08 / 2.67 / 2.03 | +0.00 | 0.98 | |

Reading: a consistent *hint* that the deepest SVI patches are active for a shorter time (same sign in all four SVI
variants), not a result.

### Q1 duration vs episode size (Rule Z) — FAIL
Episode median D80 vs PNSN tremor count: SVI ρ +0.30 (p 0.41), NW +0.43 (p 0.12). Vs GNSS offset: SVI +0.32 (p 0.36), NW
+0.47 (p 0.11). Vs LFE count: +0.25 / +0.51 (p 0.48 / 0.07).
**Post-hoc:** vs |front speed| (from S3): NW ρ −0.91, p 0.0001 (−0.93 for migrating episodes only); SVI −0.35, p 0.34.
Band width = D80 × speed: SVI median 13.4 km (IQR 9.8–16.5), NW 14.6 km (IQR 10.1–17.7). Patch-crossing time is not the
explanation: SVI patches need a median 1.8 d for the front to cross them (D80 2.9 d); NW cells 0.3 d (D80 3.1 d); D80
exceeds 2× the crossing time in 44 % (SVI) and 81 % (NW) of cells.

### Q2 burst structure (Rules B1–B3)
| run | cells | median CV | Poisson (uniform) | shuffle within day: median CV (q99) | within hour (q99) | events in bursts: obs / day-shuffle / hour-shuffle |
|---|---|---|---|---|---|---|
| SVI headline | 72 | **2.76** | 0.97 | 1.72 (1.85) | 2.68 (2.72) | 70 % / 25 % / 52 % |
| SVI fixed-S13 | 37 | **2.59** | 0.98 | 1.61 (1.76) | 2.48 (2.53) | 71 % / 26 % / 51 % |
| NW headline | 115 | **2.47** | 0.96 | 1.59 (1.67) | 2.41 (2.44) | 74 % / 26 % / 58 % |
| SVI scrambled control (P ≥ 0.5, cells ≥ 10) | 33 | 1.18 | 0.84 | 1.04 (1.14) | 1.16 (1.18) | 42 % / 18 % / 36 % |
| post-hoc: SVI merged < 60 s | 69 | 2.48 | | 1.60 (1.72) | 2.42 (2.45) | 63 % / 24 % / 47 % |
| post-hoc: NW merged < 60 s | 95 | 2.19 | | 1.46 (1.54) | 2.17 (2.20), p 0.27 | 64 % / 25 % / 51 % |

- B1 "clustered beyond the daily rate": **PASS** (p 0.002 for SVI, NW and fixed-S13). "Clustered below 1 h": formally passes
  (p 0.002), but the excess is small and disappears in NW once events < 60 s apart are merged. **The bursts are hours long.**
- B2 scrambled control: CV ratio to the within-day shuffle is 1.12, against 1.53 for the data (outside the 20 % band), so
  the clustering is **not** called detection-driven. The scrambled control has 4 % of gaps < 30 s; the data have 13 %
  (SVI) and 17 % (NW), 67–78 % of them from the same family/source. These may be real repeats seconds apart or residual
  double detections; merging them does not change B1.
- B3 depth: FAIL. Within-episode ρ(CV, dd) +0.08 (p 0.55) SVI, +0.14 (p 0.24) NW; burst fraction p 0.22 / 0.57. Median
  CV shallow / mid / deep: SVI 2.64 / 2.96 / 2.78; NW 2.19 / 2.30 / 3.00 (NW deep higher, not significant within episode).
- Burst sizes (sizes ≥ 2, gap < 0.1 × mean): SVI n = 1,650, mean 3.2, max 19, P(≥ 5) 16 %, Vuong z −3.3 → exponential
  (p 0.0009); fixed-S13 z −2.6 (p 0.011, undecided); NW n = 1,543, mean 3.7, max 24, P(≥ 5) 22 %, z −0.4 → undecided
  (merged < 10 s: z −2.6, p 0.0096 → exponential). Within-day shuffle: mean 2.2–2.3, max 7–9, P(≥ 5) 0.8–0.9 %. In the
  figure the observed tail lies between the exponential and power-law fits.

### Q3 re-activation (Rule R) — rare, FAIL for depth
| run | cells ≥ 10 | re-activating | expected if the patch followed the whole region's timing | shallow / mid / deep | ρ (dd), p |
|---|---|---|---|---|---|
| SVI headline | 90 | 5 (5.6 %) | 4.2 % | 3 % / 3 % / 12 % | +0.09, 0.50 |
| SVI fixed-S13 | 50 | 5 (10 %) | 1.6 % | 12 % / 5 % / 18 % | −0.01, 0.97 |
| NW headline | 191 | 11 (5.8 %) | 4.2 % | 7 % / 2 % / 8 % | +0.03, 0.71 |

Returns come a median of 8–9 days after the patch's busiest day. 40–45 % of them fall within ±1 day of another patch's
return (SVI 2018-06: groups 9 and 11 on Jul 8 PT; NW 2021-09: four patches on Oct 20–24 PT; NW 2015-12: two deep bins on
Jan 13–14 PT). PNSN tremor in the same along-strike range also returns in 1 of 5 SVI cases and 7 of 11 NW cases. Deep SVI
group 10 re-activates in 2015-12 and 2019-11 (also on fixed-S13). Too few cases to call this a depth effect.

## 3. Comparison with prior work (honest)

- **Omori-like decay:** we see none at the patch level. Local LFE activity is controlled by the front (on for ~2–3 days,
  then off), not by aftershock-like relaxation. Any study that stacks rate against time since front arrival will see an
  Omori-like or exponential tail produced by mixing patches with different on-times. Our rate-vs-time stacks warn against
  reading that tail as a decay law.
- **Local duration of a few days and an active band ~15 km behind the front** agree with tremor-based descriptions of
  Cascadia ETS, where tremor persists for a few days behind the main along-strike front (Houston et al. 2011; Rubin &
  Armbruster 2013; Wech & Creager). New here: the same picture from LFE families at patch level over 2015–2026. The
  duration ∝ 1/speed relation (NW only, post-hoc) implies a band of about constant width, which is a testable claim.
- **Wech & Creager (2011) down-dip continuum; Bostock et al. (2015) recurrence by depth in SVI:** the strong depth
  contrasts there are *between* episodes (deeper = more frequent, smaller events; cf. our S7: between-episode LFEs are deep).
  *Within* the passage of a big ETS front we find no robust depth dependence of duration, burstiness or re-activation, only
  an SVI hint that the deepest patches switch off sooner. Both fit together: depth matters for how often a patch slips, not
  for how it behaves while the big front crosses it.
- **Frank et al. (Guerrero) LFE bursts:** we also find LFEs strongly clustered in bursts lasting hours inside the active period.
  Frank et al. report burst behaviour changing along dip; within Cascadia ETS fronts we do not see that (B3 fails in both
  regions).
- **Thomas et al. 2016/2018; Shelly & Johnson 2011 (Parkfield):** families there range from near-continuous to episodic.
  During ETS our patches all behave like the episodic type: a burst of a few days with hour-scale clustering inside.

## 4. Caveats
- Detection completeness differs by depth and over time. All depth tests compare cells within the same episode, and
  durations are thinned to 30 events; the fixed-13-station run (6 episodes) gives the same signs. These steps reduce the
  problem but do not remove it (deep SVI families have fewer carrier stations).
- NW down-dip positions of single sources are uncertain (S3). NW depth classes come from source dd_km thirds.
- The plateau fits hit the sharpness bound (w = 0.02 d), so the stop is sharper than we can resolve, not exactly 0.02 d.
- The post-hoc checks (speed relation, single sources, merging < 60 s) were not pre-registered and are provisional.
- No Merlin review was done before running (caller asked for direct work).
