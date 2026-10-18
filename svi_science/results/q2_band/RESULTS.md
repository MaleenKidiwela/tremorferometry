# Q2: is the active band behind the slip front the same in every ETS, 2003–2026? Results

The tests and pass rules were frozen first in `PLAN.md` (2026-10-07, ~2:25 pm PT). The main run went from 2:28 to 2:58 pm PT. Post-checks
that were not pre-registered ran at ~3:00 pm PT and are marked **post-hoc** below.

Scripts:
- `svi_science/scripts/q2_band.py` (pre-registered)
- `q2_post.py` (post-hoc)
- `q2_figures.py`

Outputs are in `results/q2_band/`: `episodes.csv`, `units.csv`, `controls.csv`, `thinning*.csv`, `s13.csv`, `tremor.csv`,
`result.json` and `result_post.json`. No Merlin review (the caller asked for direct work), so treat this as provisional.

**Data.** All four datasets use the same code and the same families or sources:

| dataset | catalogue | big episodes | LFEs |
|---|---|---|---|
| SVI 2003–13 | Bostock | 10 | 256,920 |
| SVI 2014–26 | ours | 10 | 7,753 |
| Washington-NW 2007–14 | ours | 10 | 10,597 |
| Washington-NW 2015–26 | ours | 14 | 9,295 |

The headline set is P ≥ 0.75, counted once.

**Band width is measured two ways:**
- **From duration:** |front speed| × the median active time of single sources in the episode. Active time is D80, the 10–90 % time,
  taken over 30-event subsamples so that it does not depend on the number of events.
- **Direct:** the along-strike extent (10–90 %) of LFEs inside 12-h windows, also from 15-event subsamples.

Front speed is a robust line fit of along-strike position against time, from LFEs, in all eras. Before 2009 there is no tremor
catalogue, so LFEs are the only source in the early years.

## 1. Plain-words answer

- **Yes: the band has the same width in every episode, over 23 years, in both regions.**
  - In southern Vancouver Island (SVI), the band from duration is about **15 km** wide.
    - Bostock's 2003–13 catalogue: 14.7 km in each episode (interquartile range 14.2–15.3, scatter 8 %).
    - Ours, 2014–26: 15.6 km.
  - In Washington-NW it is about **11–12 km** (2007–14: 10.6 km; 2015–26: 11.8 km).
  - No episode differs from the others by more than its measurement error.
  - There is no trend with year.
  - Over ~20 years the band changed by less than about +7 to +9 km and by less than −2 to −3 km (95 % ranges).
- **The structure is stationary: slower fronts are active longer at each source.**
  - Duration × speed stays fixed: a source switches on when the front arrives and stays on for about 2–3 days.
  - Duration against speed: ρ −0.90 in Washington-NW (both eras) and −0.83 in Bostock's SVI.
- **Bostock's catalogue and ours give the same band once detector recall is matched.**
  - Thinning Bostock's catalogue down to our number of LFEs leaves the median width unchanged.
  - It does raise the episode-to-episode scatter from 8 % to about 33 %, which is the scatter we see in our own catalogues.
  - So our noisier per-episode widths come from having fewer detections, not from real changes.
- **The band does not depend on episode size or migration direction.**
- **Depth: in SVI the deep (down-dip) edge of the band is narrower; in Washington-NW it is not.**
  - In Bostock's SVI catalogue, deep sources stay active for 2.3 days against 3.5 days for shallow ones (p 0.0004).
  - That gives about 11 km of band width at depth against 16 km up-dip.
  - Our SVI catalogue shows the same sign but is not significant, so by the frozen rule this stays provisional.
  - Washington-NW shows nothing in either era.
- **SVI vs Washington: possibly about 3 km wider in SVI, not robust.**
  - SVI sources stay active longer (3.0 vs 2.1 days) while fronts move at the same speed (4.9 vs 4.8 km/d).
  - The difference disappears in the common era and in the patch and arrival-time versions.
  - The direct measure points the other way.

Figures:
- `figures/q2_band_width_by_year.png`: band width per episode, 2003–2026, by region, with the scrambled control and thinned Bostock.
- `figures/q2_duration_vs_speed.png`: duration vs speed with constant-width lines; band width by down-dip third; source duration vs
  down-dip position.
- `figures/q2_thinning.png`: Bostock full catalogue vs thinned to our count.

## 2. Numbers with controls

### Width per era (single-front episodes: front p < 0.01 and |v| ≥ 2 km/d)
| era | n ep | width from duration, median (IQR) | direct, median | \|v\| | source D80 | episode-to-episode scatter (CV) | Cochran Q p (excess scatter τ) |
|---|---|---|---|---|---|---|---|
| SVI Bostock 2003–13 | 8 | **14.7 km** (14.2–15.3) | 18.6 km | 4.8 km/d | 3.1 d | 8 % | 0.65 (0 %) |
| SVI ours 2014–26 | 9 (8 with a source ≥ 30) | **15.6 km** (13.7–19.3) | 14.8 km | 5.6 km/d | 2.7 d | 40 % | 0.24 (16 %) |
| WA-NW 2007–14 | 9 | **10.6 km** (9.3–12.1) | 25.6 km | 4.3 km/d | 2.4 d | 25 % | 0.96 (0 %) |
| WA-NW 2015–26 | 11 | **11.8 km** (10.2–14.2) | 24.1 km | 6.2 km/d | 2.0 d | 30 % | 0.17 (20 %) |

Region pooled, Cochran's Q:
- SVI: p 0.52 (τ 0)
- NW: p 0.62 (τ 0)

So no episode is wider or narrower than its error allows (T2).

Episodes left out of the width-from-duration tests because the front does not migrate:
- SVI-B: 2007-01, 2013-09
- SVI-O: 2021-10
- NW-E: 2012-04
- NW-O: 2017-03, 2018-05, 2020-10

SVI-O 2022-03 has no source with ≥ 30 events. With a 15-event floor it gives 12.3 km.

### T1, scrambled control (source positions shuffled among events within each episode, front speed kept): PASS everywhere
| era | width from duration: observed / control median | direct: observed / control median |
|---|---|---|
| SVI-B | 14.7 / 51.2 km | 18.6 / 52.7 km |
| SVI-O | 15.6 / 39.8 km | 14.8 / 43.2 km |
| NW-E | 10.6 / 45.8 km | 25.6 / 55.2 km |
| NW-O | 11.8 / 50.6 km | 24.1 / 53.9 km |

All are below the 1st percentile of 200 control sets (p 0.005, the smallest p possible with 200 sets). The band is a real, narrow
feature. It is not "everything is active for the whole episode".

### T3 trend with year, and T4 eras: no change
- **Width from duration, trend:**
  - SVI 2003–26: +0.08 km/yr, p 0.27. Change over 21.5 yr between −3.1 and +9.1 km.
  - NW 2007–26: +0.12 km/yr, p 0.26. Change between −2.1 and +6.7 km.
  - **Post-hoc** 15-event floor: SVI −4.4 to +2.7 km (p 0.51); NW −1.7 to +7.8 km (p 0.18).
  - Speed and source duration show no trend either (p 0.39–0.66).
- **Eras, width from duration:**
  - SVI Bostock vs ours: 14.7 vs 15.6 km, difference −0.9 km (95 % −6.9 to +1.7), p 0.57.
  - Same comparison with Bostock thinned, largest events kept: 15.7 vs 15.6 km, p 0.96.
  - **Post-hoc** 15-event floor: 14.1 vs 13.5 km (−1.8 to +2.0), p 0.67. Random thinning p 0.24; largest-M thinning p 0.09.
  - Washington 2007–14 vs 2015–26: 10.6 vs 11.8 km, −1.2 km (−5.7 to +2.0), p 0.32.
- **Direct width in SVI:** narrower in our era (18.6 vs 14.8 km, p 0.007; trend p 0.008). This is not a real change. Thinned
  Bostock gives +2.1 km (p 0.08, random) or −1.6 km (p 0.52, largest-M), so the sign does not hold and by the frozen rule no era
  difference is called. The denser catalogue simply catches more distant events active at the same time.
  - Washington direct width: 25.6 vs 24.1 km, p 0.88.

### T6 thinning: Bostock thinned per episode to 639 LFEs (the median of our SVI big episodes)
- Median ratio thinned / full:
  - front speed: 1.01 (random), 1.14 (largest-M)
  - width from duration: 1.06 (largest-M)
  - direct width: 0.94 (random), 0.71 (largest-M)
- Random thinning leaves **no** source with ≥ 30 events, because the LFEs are spread over 115 families. So the pre-registered
  30-event width is undefined for random thinning.
- **Post-hoc**, 15-event floor: ratio 1.03 (random) and 1.17 (largest-M). Patch version: 1.07 / 0.91.
- A single thinned draw has an episode-to-episode CV of 0.33 (10–90 %: 0.22–0.52). This matches our own catalogues (0.25–0.40).
- Our lower recall therefore keeps the median width and only adds scatter. The Bostock era's 8 % scatter is the better estimate of
  how stable the band really is.

### T5 regions (SVI vs Washington-NW)
- Width from duration, all eras: SVI 15.0 vs NW 11.2 km, difference +3.9 (2.1–5.6), **p 0.004**.
  - Post-hoc 15-event floor: 13.7 vs 11.3, p 0.002.
- The cause is longer source durations in SVI (3.0 vs 2.1 d, p 0.015). Speeds are the same (4.9 vs 4.8 km/d).
- It is **not robust**:
  - common era 2014–26: p 0.09
  - patch version: p 0.23
  - arrival-time speed: p 0.18
  - the direct width goes the other way (SVI 16.7 vs NW 24.6 km, p 0.017; 14.8 vs 24.1 km in 2014–26, p 0.005)
- The direct width depends on the source layout: SVI families sit in a compact strip, while Washington-NW sources are spread
  ~120 km along strike and several patches can be active at once. Call this a possible difference of ~2–4 km, not a result.
- Sources in Washington-NW are source groups (`sources_nw.csv`), not single families. This is a difference in method between the
  regions; it does not affect the comparisons between eras.

### T7 depth (within episode, source D80 vs down-dip coordinate; permutation)
| era | ρ | p | width shallow / mid / deep (km) | D80 shallow / mid / deep (d) |
|---|---|---|---|---|
| SVI Bostock (1,021 sources-episodes) | **−0.13** | **0.0004** | 16.4 / 14.8 / 11.4 | 3.5 / 3.2 / 2.3 |
| SVI ours (63) | −0.17 | 0.27 | 14.9 / 15.1 / 5.6 | 2.7 / 3.0 / 1.0 (few deep sources) |
| WA-NW 2007–14 (109) | −0.06 | 0.68 | 10.1 / 10.0 / 11.7 | 2.4 / 2.1 / 2.9 |
| WA-NW 2015–26 (81) | −0.13 | 0.35 | 16.6 / 13.9 / 12.2 | 2.5 / 2.1 / 2.0 |

The frozen rule needs p < 0.01 in both eras of a region, so it fails. SVI is **provisional**: the deep edge is narrower, strongly
in Bostock and with the same sign in ours. This also matches C (SVI deepest patches shorter in all four variants). Washington-NW
shows no depth dependence.

### T8 episode size and T9 direction: no dependence
- **Width from duration vs size**, permutation within era:
  - SVI: ρ +0.03 against along-strike extent, −0.29 against episode duration, −0.23 against PNSN tremor count, +0.01 against LFE
    rank; all p ≥ 0.29.
  - NW: ρ −0.07 / +0.09 / +0.21 / +0.36 for the same four; all p ≥ 0.12.
- **Direct width in NW does grow with episode extent** (ρ 0.60, p 0.0015) and with LFE rank (0.58, p 0.005). The direct measure
  sees how much of the margin is active at once, not the band behind the front.
- **Duration vs speed:**
  - NW ρ −0.90, p 0.0001 (2007–14: −0.77; 2015–26: −0.87)
  - SVI Bostock −0.83, p 0.010
  - SVI ours −0.10, p 0.82 (8 episodes, few sources each)
  - SVI all, 15-event floor: −0.59, p 0.013
- **Direction (T9):**
  - NW: SE-ward 11.4 vs NW-ward 11.0 km, p 1.0.
  - SVI (NaN-safe, post-hoc): 13.9 vs 15.2 km, p 0.18.
  - 15-event floor: p 0.15 (SVI) and 0.35 (NW).

### T10 fixed-S13 stations (SVI, 6 episodes 2015-12 … 2021-01): PASS by rule (4 of 6 agree within error)
- Ratios S13 / headline for width from duration: 0.99, 1.62, 1.25, 1.44, 0.68, and undefined for 2021-01 (no S13 source ≥ 30).
  2018-06 disagrees.
- Direct width ratio: 0.99.

### T11 tremor check (PNSN, 25 single-front episodes from 2010)
- Tremor front speeds agree with LFE front speeds: ρ 0.81.
- The tremor band from duration is wider: 22–28 km, a median of 1.7× the LFE value. It uses 5-km bins plus epicentre errors.
- Tremor width does not follow the LFE width episode by episode (ρ −0.16, p 0.46).
- Tremor width also shows no trend (SVI p 0.63, NW p 0.73) and no era difference (SVI p 1.0, NW p 0.27).
- So tremor confirms the stationarity but cannot resolve the 10–15 km band itself.

### Columns
- Patch width, as in C: SVI 17.9 km, NW 16.6 km. This matches C's 13–15 km plus the time the front needs to cross a patch.
  - Its SVI era difference (18.8 vs 14.4 km, p 0.015; random thinning p 0.006; largest-M p 0.24) comes from patch size. Bostock
    lights all families in a group, while our detections come from fewer families. Not called.
- Arrival-time speed × D80: SVI 19.6 km, NW 16.0 km, no trend (p 0.53 / 0.86).
- The distance behind the fitted front line (26 / 45 km) is dominated by the front not moving at a constant speed. Not used.

## 3. Comparison with prior work
- **Houston et al. 2011; Rubin & Armbruster 2013; Wech & Creager:** tremor stays active for a few days behind the main front, a band
  of order 10–20 km. We find the same with LFE families, and show for the first time that it is **the same in every big ETS over
  2003–2026**, in two regions and with two independent detectors.
- **Slower fronts, longer local activity:** C found this in NW 2015–26 (post-hoc). It now holds in NW 2007–14 and in Bostock's SVI
  2003–13, so it is a property of the front, not of one dataset.
- **Bostock et al. 2015 (SVI) recurrence by depth:** deep families repeat more often in smaller events. Our narrower deep edge in SVI
  (provisional) fits that picture: deep sources may release their share in a shorter burst while the front passes.

## 4. Caveats
- 8–11 single-front episodes per era. "Same" means no difference at this precision; the 95 % ranges above are the claim.
- Our catalogues have few sources with ≥ 30 events per episode (SVI 0–16, NW 2–27). Per-episode widths are therefore noisy, and two
  of our SVI episodes rest on one or no source. Bostock's era (100+ sources per episode) carries the precision.
- Width from duration is |v| × median D80, so it assumes one front speed per episode. Episodes with no single front are excluded.
- The definitions of a "source" differ between regions (SVI single family; NW source group). This affects the region comparison
  only.
- Post-hoc fixes:
  - The thinned-Bostock columns had been joined by episode name and leaked onto WA-NW rows with the same name. They were blanked
    in `episodes.csv`; no pre-registered statistic used those rows.
  - Several statistics were made NaN-safe.
