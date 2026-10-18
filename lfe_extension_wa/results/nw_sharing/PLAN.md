# NW event sharing between families — PLAN (frozen before any pair count)

Written 2026-10-07 ~12:10 pm PT, before any pair count was made. Script: `lfe_extension_wa/scripts/nw_sharing.py`
(sub-commands). Outputs only in this folder. Nothing is deleted; no production file is edited. Later changes go in
"Revisions" with the time and reason.

Question: do some of the 132 NW headline families (`wa_prod/setup/families.csv`, `reg == "NW"`) keep picking up the same
events on the permanent network, so that they are not distinct sources there? What does that do to the NW family count
and to the per-family / per-episode rate results?

## 1. Pair table (all 132 x 131 / 2 NW pairs, per era E1 / E2 / E3)
- Repeats of family i = headline data candidates: core network, series data, P_own_WAVE_PNET >= 0.75
  (`conf/P/prod/P_*.parquet`, the rows `rates_daily_v2` counts). N_i(era) = their number.
- Same-event pair count n(i<-j) = number of i's headline repeats that have at least one candidate of family j (core, data,
  any P) with
  - |t_i - t_j| <= 2 s, and
  - the same predicted arrival at their best shared era-core station:
    |(t_i - t_j) - Doff(i, j)| <= 0.5 s (Doff = offset difference of the two templates at the shared station with the
    largest min(template weight); `wa_post_dedup.share_mats`, imported read-only),
  - and station sharing |A ∩ B| / min(|A|, |B|) >= 0.5 over the era-core template stations.
  This is the variant D pair condition (any distance + arrival check) without the "higher net_cc wins" ordering, so it
  is symmetric in meaning (both directions are counted).
- Columns kept for every pair: n without the station-sharing condition; n with only the +-2 s window (no arrival
  check); n with the partner also a headline repeat; the directional D-flag winner count from `post_v2/dedup/flags/`.
- Sharing fraction: s(i, j) = n(small <- big) / N_small, small = the family with fewer headline repeats in that era.
  Both directional fractions are columns.
- Chance level c(i, j): the same count with family j's candidates shifted by +-60, 120, 300, 600 s (8 shifts, mean),
  divided by N_small.
- Other columns: epicentral distance (lat_use / lon_use), depth difference, shared era-core template stations (list,
  overlap), weighted template-station Jaccard, catalogue v1.2 run of origin, `union_nn_km`, `union_flag_lt2km`,
  `dup_group`, location quality, erh / erz / n_picks, Royer-truth match (one-to-one Hungarian on the digitized NW
  Hyp2000 templates, matched if <= 5 km; distance a column).
- Template similarity: the NW union and catalogue v1.2 tested only pairs <= 10 km apart (NW union) or cross-catalogue
  pairs (v1.2), so no NW-NW similarity exists for most pairs. It is computed here with the NW union rule unchanged
  (array-era 2-8 Hz stacks `hypo_28/stacks_loc.npz`; >= 2 common stations with cc >= 0.5 and lag spread <= 0.15 s =
  "duplicate"; spread > 0.3 s = "neighbour") for every linked pair and for a control set of unlinked pairs at matched
  distance. Array-era co-detection: share of the smaller family's original discovery detections (`det_source`) with a
  partner detection at the pair's modal offset +-0.25 s, minus the same with the partner shifted by +-60 / 120 s.

## 2. Link rule and source groups (PRE-REGISTERED)
- Primary era: E3 (2015-01-01 .. 2026-09-15, fixed E3 core network; the headline era).
- **Link i-j iff s(i, j) >= 0.30 AND s(i, j) >= 5 x c(i, j) AND n(small <- big) >= 10.**
- Source groups = connected components of the link graph (singletons are groups of one).
- Reported: number of groups (= families effectively distinct on the permanent network), group sizes, members; per group
  the fraction of member pairs directly linked (chaining check); groups' max member separation.
- Sensitivity (columns only, not the answer): link threshold 0.20 and 0.50; E1 and E2 graphs with the same rule;
  cross-era persistence of E3 links.
- Does grouping follow: location (all member pairs <= 10 km), catalogue v1.2 dup_group / union_flag_lt2km, run of
  origin (pairs from the same run vs different runs, against the all-pair base rate), station sharing (share of the
  linked pairs whose top-weight shared station is B013 or B943, against all pairs).

## 3. Location error or distinct sources?
- For linked pairs: separation d vs sigma = sqrt(erh_i^2 + erh_j^2) (array-era Hyp2000). "Within location error" if
  d <= 2 sigma. Depth difference vs sqrt(erz_i^2 + erz_j^2) the same way. Each family's n_picks and quality.
- Same numbers for all unlinked pairs at the same distances (base rate).
- Array-era waveform test and co-detection (above): a linked pair that is a "duplicate" or co-detects in the array era
  is one source mislocated; a linked pair that is a "neighbour" with no array-era co-detection is two sources the
  permanent network cannot separate.

## 4. Science re-run at the source-group level (sensitivity row, k-free)
- Episodes, windows, exposure, the share test, the neighbourhood rule (40 km, else region), bootstrap, BH, labels and the
  residual-correlation regression are the episode_rates functions (imported read-only; `episode_rates.py` is not edited,
  its default behaviour unchanged): NW E3 MAJOR episodes of `episode_rates/episodes_nw_final.csv`, not excluded.
- Group unit: window = union of the members' episode windows; position = mean of member positions; count N = group
  repeats at P >= 0.75 minus the time-reversed control, where member candidates of the same group are counted once
  (a member candidate is dropped if a higher-net_cc candidate of another member lies within +-2 s with the arrival
  check of section 1; same for the control series); D = number of distinct exposure days of the members; X = mean of the
  members' X (secondary column only); carrier Jaccard = weighted Jaccard on the members' max template weights.
- Acceptance check first: the same code with singleton groups must reproduce the family-level V0 result (89 tested,
  58 over-dispersed, median chi2/df 4.03).
- Report next to the family-level (V0) and the de-duplicated (`svi_science/results/s0b/before_after.csv`, NW E3 at
  0.75) rows: tested, over-dispersed / tested, median chi2/df, distance and carrier coefficients with permutation p,
  residual correlation by distance bin.

## Revisions
(none yet)

### Revision 1 — 2026-10-07 ~12:20 pm PT (after the array-era tests ran, before the episode re-run)
PLAN sec. 3 named the two array-era tests but gave no number for "co-detects". Fixed now, from the control set (174
unlinked pairs at matched distance; only 4 of them, 2.3 %, reach 0.1): a linked pair is **array-confirmed same source**
if array-era co-detection excess >= 0.1 OR the waveform test says "duplicate" on >= 5 stations (the NW union rule asks
only >= 2 stations, and 9 of 174 control pairs pass on 2-3 stations at cc ~ 0.5, so 2 stations is too weak). Otherwise the
linked pair is **distinct on the array** (two sources the permanent network cannot separate). This is a diagnostic
classification chosen after seeing the linked-pair results; it is not part of the pre-registered link rule. The episode
re-run (sec. 4) gets one extra sensitivity row: groups built from array-confirmed links only.
