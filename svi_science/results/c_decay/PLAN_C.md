# C — LFE activity at a patch after the slip front passes (SVI and Washington-NW)

Pre-registered 2026-10-07 ~1:25 pm PT, before any number of this study was computed. (No Merlin review: caller asked
for direct work; rules below are frozen.) Script: `svi_science/scripts/c_decay.py`. Times reported in PT.

## Data (headline and columns)
- **SVI headline**: `results/s3/events_svi.parquet`, series data, P_same ≥ 0.75, counted_once (keep_xfam8), episode in
  the 10 big SVI episodes of `results/s3/episodes_kinematics.csv`. Columns: P ≥ 0.5 counted once; raw P ≥ 0.75.
- **SVI scrambled control** (stations shifted to different days): series c1, P ≥ 0.5, counted once (P ≥ 0.75 has only
  49 rows in big episodes — too few). Used pooled only (Q2).
- **SVI fixed stations**: `results/s0/analysis_table_S13.parquet` (run prod_S13, data, P ≥ 0.75, keep_xfam8), assigned
  to the same padded episode windows (s3_build rule). Covers to 2021-03-31 → 6 big episodes (2015-12 … 2021-01).
- **NW headline**: `results/s3/events_nw.parquet`, data, P_same ≥ 0.75, selected_v2, counted_once (dedup D), episode in
  the 14 big NW episodes. Position = source (`sources_nw.csv`). No scrambled-control rows exist for NW → time-shuffle
  controls only. Quality-C sources flagged as a column (depth uncertain).

## Patches and depth classes
- SVI patch = the 12 pre-registered S0 groups. NW patch = 10-km along-strike bin (as S3 §7) × down-dip third (below).
- Depth class: SVI = family down-dip coordinate thirds over the 115 families (same cut as S3b); NW = source down-dip
  coordinate (dd_km) thirds over the 100 sources. A patch's class/dd = event-weighted mean of its members' dd.
- Cell = patch × big episode. Analysis floors: Q1 duration / Q2 ≥ 30 headline events; Q3 ≥ 10 events.

## Local front arrival (detection-independent, as S3b)
Per episode, per 5-km along-strike bin: arrival = 10th-percentile time of PNSN tremor epicentres in the bin (≥ 20
epicentres). Lag of an LFE = time − arrival(bin of its along-strike position). Events in bins with no arrival: dropped from
lag work only (count reported). Column: LFE-own arrival = 5th percentile of the cell's own event times.

## Q1 — rate vs time since front arrival: shape, duration, depth, episode size
- **Stack**: per region × depth class, histogram of lags (0.25 d bins, −10 … +30 d), event-weighted (column:
  each episode weighted equally).
- **Shape fit**: lags t' = lag − t_pk in [0, 25 d], t_pk = centre of the highest 0.5-d bin in [−2, +10] d. Three
  2-parameter densities, maximum likelihood on the truncated window:
  M1 Omori-like (t'+c)^−p; M2 exponential with floor exp(−t'/τ)+b; M3 plateau-then-stop 1/(1+exp((t'−T)/w)).
  Same fits per cell (≥ 30 events after the cell's own peak bin); cell winner = lowest AIC with ΔAIC ≥ 2, else undecided.
- **Rule S (shape)**: a region's shape is called X if X wins the stack by ΔAIC ≥ 10 AND X is the most common decided winner
  among cells (≥ 40 % of decided cells and ≥ 1.5× the next). If the stack winner and the cell winner differ, the stack shape
  is reported as a stacking product, the cell shape is the answer.
- **Duration** per cell: D80 = t90 − t10 of the cell's event times; count-free version D80_n30 = median over 200
  random 30-event subsamples (headline). Column: tail T90 = t90 − cell arrival (median arrival of its events).
- **Rule D (depth)**: within-episode Spearman ρ(D80_n30, cell dd), averaged over episodes with ≥ 3 cells; control =
  dd labels permuted among cells within each episode (5,000) → two-sided p. "Duration depends on depth" if p < 0.01 in a
  region AND (SVI) the fixed-S13 run gives the same sign. Passing in one region only = provisional.
- **Rule Z (episode size)**: Spearman ρ(episode median D80_n30, episode size); size = PNSN tremor count in the box
  (detection-independent; columns: headline LFE count, S1 GNSS offset D (SVI) / D_WA (NW)). Permutation 10,000:
  p < 0.01 = depends on size; 0.01 ≤ p < 0.05 = hint only (10 / 14 episodes, weak power).

## Q2 — burst structure inside the active period
- Per cell (≥ 30 events): inter-event times (IET) between consecutive events, excluding gaps > 3 d (those are Q3).
  CV = std/mean (Poisson = 1). Burst = run of events linked by IET < g; headline g = 0.1 × the cell's mean IET
  (columns g = 10 min, 60 min). Fraction of events in bursts (size ≥ 2).
- **Controls (time shuffles within each cell)**: S_uni = same count uniform over [first, last] (Poisson reference);
  S_day = each event moved to a random time within its own PT day (keeps daily counts); S_hour = within its own hour
  (keeps hourly counts → tests clustering below 1 h). 500 sets each.
- **Rule B1**: "clustered beyond the daily rate" if the median CV over cells exceeds the 99th percentile of the same median
  over the 500 S_day sets (p < 0.01), in a region; SVI must also pass on fixed-S13. Same with S_hour → "clustered at < 1 h".
- **Rule B2 (scrambled control)**: SVI c1 (P ≥ 0.5) per-cell CV vs its own S_day; if the scrambled control is clustered
  beyond S_day as strongly as the data (ratio CV/CV_Sday within 20 % of the data's), the data clustering is called
  possibly detection-driven. If c1 has < 5 cells with ≥ 30 events, use cells ≥ 10 and say so.
- **Rule B3 (depth)**: CV and burst fraction vs cell dd with the Rule-D within-episode permutation, p < 0.01.
- **Burst sizes**: pooled per region (and per depth class), sizes s ≥ 2; discrete power law (zeta, s ≥ 2) vs discrete
  exponential (geometric, s ≥ 2), Vuong normalized log-likelihood ratio: power law if R > 0 with p < 0.01, exponential if
  R < 0 with p < 0.01, else undecided. Compare tail P(s ≥ 5) with the S_day sets.

## Q3 — re-activation later in the same episode
- Per cell (≥ 10 events): re-activation = (S3 rule at patch level) a PT day with ≥ 3 events, then ≥ 3 consecutive days with
  0 events, then a day with ≥ 3 events. Fraction of cells re-activating per depth class; delay of the return.
- **Chance control**: per cell, n_cell times drawn from all headline events of the same episode and region (the
  patch follows the region-wide time history), 500 sets → expected re-activation fraction.
- **Rule R (depth)**: statistic = Spearman ρ(reactivated 0/1, cell dd), control = dd labels permuted among cells within
  episode AND within count tercile (5,000). p < 0.01 = depends on depth; SVI also checked on fixed-S13.
- Column: whether PNSN tremor in the patch's along-strike range (± 5 km) also shows the return (≥ 3 epicentres on a return
  day after ≥ 3 days with ≤ 1), and whether returns in an episode fall on the same days across patches (front-driven).

## Not claimed whatever the numbers
Absolute per-patch completeness differs by depth (deep families have fewer carriers); within-episode comparisons and
count thinning reduce but do not remove this. Any depth result passing only in one region stays provisional.
