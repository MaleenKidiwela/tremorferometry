# Q2 — Is the active band behind the slip front the same in every ETS, 2003–2026? (SVI + Washington-NW)

Pre-registered 2026-10-07 ~2:25 pm PT, before any number of this study was computed. Known beforehand (from C, post-hoc):
NW 2015–26 duration × speed ≈ 15 km, SVI ≈ 13 km. No Merlin review (caller asked for direct work). Script
`svi_science/scripts/q2_band.py` (+ `q2_figures.py`). Times reported in PT. CPU only, < 30 GB.

## Data and eras (headline = P ≥ 0.75, counted once — same as C)
| region / era | events | episodes |
|---|---|---|
| SVI-B 2003–2013 (Bostock) | `p3_bostock/bostock_events_mapped.parquet` (fid → our family, our along/down-dip frame) | the 10 episodes of `p3_bostock/episodes_bostock_final.csv`, events in [start − 2 d, end + 3 d) |
| SVI-O 2014–2026 (ours) | `s3/events_svi.parquet`, series data, P ≥ 0.75, counted_once | the 10 S3 big episodes |
| NW-E 2007–2014 (ours, older) | `catalogue_v2.parquet` NW, selected_v2, P_own_pnet ≥ 0.75, counted once (dedup D flags), family → source (`source_groups_E3`, as s3_build) → `s3/sources_nw.csv` position | B primary-big NW episodes before 2015 that also pass the S3 "big" rule (≥ 300 counted-once LFEs), windows w0–w1 of `b_ets_record/episodes_b.csv` |
| NW-O 2015–2026 (ours) | `s3/events_nw.parquet` as in C | the 14 S3 big NW episodes |

Unit = single source (SVI family; NW source) — primary, because a multi-km patch adds the front's crossing time to the
duration. Column: C's patches (SVI group; NW 10-km bin × down-dip third).

## Measurements per episode (identical code in every era; LFEs primary, PNSN tremor as a check)
1. **Front speed v** = Theil–Sen slope of along-strike position vs time over the episode's events (S3 function: ≤ 1,500
   random events, pairs at the same position ignored). Time-shuffle p (positions permuted, 500) and day-block bootstrap CI
   (200). Column: arrival-based speed = Theil–Sen of 5-km bin centre vs bin arrival (10th-percentile LFE time, bins ≥ 10 events).
   **Single-front episode** = p_speed < 0.01 and |v| ≥ 2 km/d. Width W_dur is only tested on single-front episodes.
2. **Unit active duration** D80_n30 = median over 200 random 30-event subsamples of (t90 − t10) — count-free; units with
   ≥ 30 events in the episode. Column floor 15 events (D80_n15).
3. **Band width, from duration: W_dur = |v| × median unit D80_n30** (primary). Uncertainty: σ²(log W) = σ²(log v) from the
   bootstrap + σ²(log median D80) from resampling units (500).
4. **Band width, direct: W_snap** = along-strike extent of simultaneous activity: non-overlapping 12-h PT windows with
   ≥ 15 events; extent = median over 100 random 15-event subsamples of (q90 − q10 of along-strike position); W_snap =
   median over windows. Count-free by construction. (Includes up to v × 0.5 d of migration inside a window; same rule
   everywhere.) Available for all episodes, migrating or not.
5. Column W_front = q10–q90 of distance behind the fitted front line (position − v·t residual), single-front episodes.

## Tests and pass rules
- **T1 band is real (control):** within each episode, permute the unit labels (positions) among events (keeps the episode's
  time history and every unit's count), 200 sets. The band must WIDEN: PASS per region-era if the median W_dur and median
  W_snap over episodes are below the 1st percentile of the same median over control sets (p < 0.01).
- **T2 same in every episode (heterogeneity):** Cochran's Q on log W_dur with per-episode σ, per region-era and per region
  pooled. "Same in every episode within measurement error" if Q p ≥ 0.01; else report the excess scatter τ (between-episode
  SD in %). Also report CV.
- **T3 stationary over the years:** per region, Theil–Sen slope of W (W_dur and W_snap) vs onset year, permutation p
  (10,000). "Changes with time" if p < 0.01. Report the 95 % CI of the change over the span (km) as the bound.
- **T4 era comparison:** SVI-B vs SVI-O; NW-E vs NW-O: difference of medians, Mann–Whitney p and bootstrap 95 % CI.
  "Eras differ" if p < 0.01 AND the same sign with Bostock thinned (T6, both thinning versions). "Same within X km" stated
  from the CI.
- **T5 regions:** SVI vs NW (all eras pooled; then common era 2014–26): Mann–Whitney, p < 0.01 = differ.
- **T6 thinning (detector recall):** Bostock thinned, per episode, to the median count of our SVI big episodes, (a) at random,
  (b) keeping his largest-magnitude events (a realistic recall model); 50 random draws for (a). Recompute v, D80_n30, W_dur,
  W_snap. Sensitivity = median ratio thinned / full. T4 must hold with the thinned values for any era call.
- **T7 depth:** within-episode Spearman ρ(unit D80_n30, unit down-dip coordinate), averaged over episodes with ≥ 3 units;
  dd labels permuted within episode (5,000) → two-sided p. "Band width depends on depth" if p < 0.01 in a region in BOTH eras
  (same sign). Also W_dur per down-dip third (|v| × median D80 of the units in that third).
- **T8 episode size:** Spearman ρ(W, size) with permutation stratified by era (10,000). Sizes: (i) along-strike extent of the
  episode's LFEs (q05–q95; comparable across eras), (ii) episode duration (t10–t90 of LFEs), (iii) PNSN tremor count n_trem30
  (2009-08 on), (iv) LFE count ranked within era. p < 0.01 = depends on size; 0.01–0.05 hint.
- **T9 direction:** SE-ward vs NW-ward (sign of v), single-front episodes, per region: Mann–Whitney with era-stratified
  permutation. p < 0.01 = depends on direction.
- **T10 fixed stations (SVI):** S13 run (6 episodes 2015-12 … 2021-01): W_dur, W_snap vs headline, paired; must agree within
  the per-episode σ in ≥ 4 of 6.
- **T11 tremor check (where PNSN exists: NW 2009-08 on, SVI 2010-08 on):** same v / D80 / W from tremor epicentres within
  10 km of our sources (as s3_build), unit = 5-km along-strike bin. Report Spearman ρ(W_LFE, W_tremor) across episodes and the
  median ratio. Tremor is a check, not a gate (epicentre errors widen W_snap).

## Not claimed whatever the numbers
"Same" means "no difference detectable at this precision"; the CI is the claim. Absolute width is bounded below by the
source spacing (SVI families cluster in groups) and above by the 12-h window migration in W_snap; comparisons between eras
use the same families and the same code, so these cancel. Before 2009-08 there is no PNSN tremor; front speed and duration
come from LFEs only in all eras.
