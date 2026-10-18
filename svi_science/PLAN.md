# Southern Vancouver Island LFE science plan (2026-10-07)

Goal: use the 12-year southern Vancouver Island (SVI) LFE catalogue (115 Bostock families, 2014-01 → 2026-09, per-repeat
confidence), its rates, the per-episode results and the modelled GNSS slow slip to find patterns not described before.
Studies run one at a time, in the order below. Each study gets its own `PLAN` section frozen before any number is
computed, a Merlin review before it runs, and a results section after. Washington-NW can repeat any study that works.

## What we have (inputs)

| Input | Path | Notes |
|---|---|---|
| Per-candidate catalogue | `lfe_extension/results/refine/conf/P_by_candidate.parquet` | 1.81 M rows; time, family, P_same, scrambled control (series c1), station set S13/S14/S15, n_sta_day, quiet flag |
| Per family-day counts | `lfe_extension_wa/results/episode_rates/counts_svi_fid_date.parquet` | P ≥ 0.5 / 0.75 / 0.85 / 0.9, data and scrambled control; run prod and fixed-13-station prod_S13 (to 2021) |
| Tremor exposure | `lfe_extension_wa/results/episode_rates/exposure_svi.parquet` | PNSN epicentres ≤ 30 km per family-day |
| Episodes | `lfe_extension_wa/results/episode_rates/episodes_svi_final.csv` | 39 major episodes 2014–2026 (frozen rule) |
| Episode consistency | `lfe_extension_wa/results/episode_rates/RESULTS.md`, `robustness/` | ~1/3 of families over-dispersed; residuals shared < 20 km |
| Monthly rates | `lfe_extension/results/refine/rates/` | Bostock-equivalent rates, k factors, step checks |
| Modelled slip | `lfe_extension/results/slip_compare/slip_svi_family_daily.parquet` | daily, thrust-positive, relative to trend; **resolves only a ~35 km blob** (RESOLUTION.md) |
| GNSS time series | `/home/jovyan/cascadia-sse-strain/data/2026-08-18/X.mat` | 14 stations in the SVI box |
| Cut repeat windows | `lfe_extension/results/cut/svi_ext/` | −15…+30 s, all stations; needed for LFE size and dv/v |
| Family locations | `lfe_catalog_repair/results/e5r_families_2012.csv` | Bostock hypocentres, mag_med |
| Known detection problems | 2017-03 (B010 out, VGZ swap), 2025-10 (VGZ 2 of 30 days), PGC changes 2015-01-12 / 2017-08-10 | carried as flags everywhere |

## Rules for every study
- Headline cut-off P ≥ 0.75; 0.5 and 0.85 as columns. The scrambled control goes through every step and gets the same table.
- Rank, don't discard: no family or episode is dropped; flags and covariates instead (station set, n_sta_day, detection-problem episodes).
- Detectability is the main confounder. Every result is checked against the fixed-13-station run (to 2021) and against n_sta_day.
- Pass/fail rules written before computing. Surprising results are provisional until a control or a second region repeats them.
- GNSS slip is used only at regional scale (whole SVI, or north vs south at most), never per family or per depth.
- Results go to `svi_science/results/<study>/`, notes to `svi_science/NOTES.md`, and the decision tree page (new section C).

## Order and studies

### S0. Shared foundations (half a day)
Build one frozen analysis table: every candidate with time, family, P, series, station set, n_sta_day, episode id,
time since the episode's start at that family, and a detectability index per family-day. Also family metadata
(location, depth, along-strike / down-dip coordinates, carrier stations). Everything below reads only this table.

### S1. Episode size from GNSS (Merlin Tier A + Tier B gate; half a day)
- Question: how big was each ETS under SVI, measured without the slip model, and can the slip model be trusted at regional scale?
- Method: station offsets across each episode at the 14 SVI-box stations (5-day medians, each station's own ETS direction,
  station bootstrap); random-window null. Then validate the slip model: resolution trace, north/south spike test,
  modelled vs observed offsets at ALBH / PGC5 / SC04, per-episode Mw vs published values.
- Pass: ≥ 15 of the major episodes above 3σ. Fail (< 8): the GNSS side of S2 is dropped; only the detection curve stays.
- Output: one geodetic size per episode, plus a north/south slip share if the gate passes.

### S2. LFEs vs slip size (needs S1)
- Question: does the number of LFEs in an ETS grow in proportion to the slip, faster, or slower? Does the ratio drift over 12 years?
- Method: log N = a + b·log D across major episodes; ratios N/D, tremor/D, N/tremor tested for trends by permutation;
  north/south LFE share vs north/south slip share. LFE moment (from family magnitudes) vs geodetic moment as a column.
- New if: b differs from 1 by more than 2σ (tremor is reported as roughly linear), or LFEs per mm of slip drift while
  tremor per mm does not (p < 0.01), or the along-strike shares match (ρ > 0.5, p < 0.01 against shuffled and rotated controls).

### S3. Slip fronts per episode (LFE-only)
- Question: where does each ETS start, which way and how fast does it move, where does it stop, and does it come back
  (secondary fronts, reversals)? Do these repeat from episode to episode?
- Method: per family per episode, first and last activation times (robust: time to 10 % / 90 % of the family's episode count);
  fit a front (start point, velocity) per episode; look for back-propagating bursts at family resolution; compare
  start points and directions across ~39 episodes. GNSS timing (regional, ~1 week) only as a cross-check.
- Controls: times shuffled within each episode; detection sensitivity map (a family seen by fewer stations activates "later").
- New if: start points or directions cluster (repeatable nucleation sites), or front speed relates to episode size (S1).

### S4. Episode fingerprints and the slip budget
- Question: is the "each ETS lights different patches" pattern random, or does a patch that was quiet in one ETS make up
  for it in the next (a slip-budget memory), or are some patches consistently early / late / strong?
- Method: per episode a share vector over families (from the existing consistency tables); correlation of residual shares
  between consecutive episodes (negative = catch-up, positive = persistent patches); clustering of episodes by fingerprint;
  relation of fingerprints to S3 start points and to the north/south slip share (S1).
- Controls: episode order shuffled; scrambled control fingerprints; detection-problem episodes flagged.
- New if: consecutive-episode residuals anti-correlate (p < 0.01): evidence that patches balance slip over time.

### S5. Tidal sensitivity per family over 12 years (needs a tide setup)
- Question: how strongly does each family respond to tidal stress, and does that change within an episode, between
  episodes, over the years, or with depth? (Tidal sensitivity is a proxy for friction and fluid pressure at the patch.)
- Setup: tidal stress on the fault plane at each family (solid-earth tide plus ocean loading; needs a tide package such
  as SPOTL or pyTMD with an ocean model; not installed yet). Shear stress in the slip direction, normal stress, Coulomb.
- Method: tidal phase of every LFE vs the phase distribution of times when the family could have fired (active windows);
  modulation amplitude per family, per episode, and early vs late in the episode at that family.
- Controls: scrambled-control detections; random times inside active windows; shuffled episodes.
- New if: sensitivity changes systematically over the 12 years or between episodes, or rises late in episodes per family.

### S6. LFE size per detection (uses the cut windows)
- Question: do LFEs get bigger or smaller with episode size, time within the episode, tidal phase, or over the years?
- Method: amplitude of each detection relative to its family template at the best stations → relative magnitude;
  size distributions per family per episode. Station gain changes (sensor swaps) handled per station epoch.
- Controls: scrambled-control windows give the noise-amplitude floor; check against instrument epochs.

### S7. Between-episode LFEs (Merlin step 7)
- Question: are the 262 quiet-day LFEs (P ≥ 0.75) tied to small slip transients, to tides, or to nearby earthquakes?
- Method: GNSS displacement rate stacked on quiet-day LFE days vs 1,000 season-matched random sets; clustering in time and space.
- Expectation: probably null at this GNSS resolution; reported as a bound.

### S8. Family-to-family triggering
- Question: within an episode, does family A firing tend to come just before family B, consistently across episodes?
- Pitfall: neighbouring families detect the same event (co-detection). Duplicate events must be merged first.
- Method: lag histograms of de-duplicated event times between family pairs; consistency of preferred direction across episodes.

### S9. Triggering by passing earthquake waves
- Question: do large distant earthquakes (M ≥ 7, 2014–2026) raise LFE rates, per family, depending on where we are in the ETS cycle?
- Method: rate in the hours after surface waves arrive vs matched windows; dynamic stress from station peak velocity.

### S10. dv/v from per-episode repeat stacks (end goal; separate phase)
- Question: does seismic velocity at depth change during or after ETS (e.g. fluid pressure)?
- Method: per family per episode stacks from the cut windows, coda stretching with the project's existing dv/v machinery
  and the fixed-station logic from the WA run.

## Order and timing
S0 → S1 → S3 and S4 (LFE-only) → S2 → tide setup → S5 → S6 → S7 → S8 → S9 → S10.
Each study: plan section frozen, Merlin review, run, results, decision tree update, user check-in before the next.

---
## Revision 1 (2026-10-07, after Merlin's review — supersedes the order and units above where they differ)

**Size check (Merlin, from the files):** at P ≥ 0.75 the production catalogue holds 13,466 repeats (scrambled control 83)
in 12 years; median family 73 events; only 390 family × major-episode cells with ≥ 10 events (111 with ≥ 30); 11 episodes
with ≥ 300 events. Per-family-per-episode designs are underpowered as first written.

**Changes:**
1. Analysis unit: de-duplicated repeats (keep_xfam8: the strongest family within 8 s, so one event is not counted by
   several neighbouring families) pooled into ~8–12 spatial family groups. Per-episode work only on the 11 episodes with
   ≥ 300 events. Pre-registered floors: ≥ 30 per cell, ≥ 300 per episode, ≥ 100 per group-episode. Both de-duplicated and
   raw variants carried as columns. Never use P_any or ΣP as a count (control mass too high); 0.75 headline, 0.5 column.
2. S0b (new, first): re-run the episode-consistency check (lfe_extension_wa/results/episode_rates) with de-duplication.
   The existing share tables did not apply keep_xfam8, so the "neighbours < 3 km move together (0.73)" result may be
   partly the same events counted by neighbouring families. 3,088 of 13,466 P ≥ 0.75 repeats are duplicates.
3. New order: S0 (+ S0b) → S1 → S2 → S7 → S3 (with S8 merged in as lag structure inside fronts) → S4 → S6 → S5 → S9 → S10.
4. S3: onset delays depend on how many events a family has; control by thinning busy families to quiet families' counts
   and by comparing LFE onset with PNSN tremor onset within 10 km (detectability-free). Frame as repeatability across the
   11 big episodes (fronts 2003–2013 already studied: Rubin & Armbruster 2013, Peng & Rubin).
5. S4: consecutive-episode anti-correlation is built in by the method (split episodes < 30 d apart, shares sum to one,
   residuals about a mean). Use only pairs ≥ 90 d apart, leave-pair-out means, episode-order shuffle; primary statistic =
   variance ratio of each group's cumulative residual (anti-persistence), not lag-1 correlation.
6. S5 tides: solid-earth tide alone is not acceptable at 30–45 km under Vancouver Island; ocean loading (Pacific shelf
   and inland waters) is the same size with a different phase → needs a tide package plus an ocean model (TPXO/FES).
   Bursts break the usual phase statistics → time-shift surrogates by whole days. Separate M2 from K1/O1; the scrambled
   control's phase histogram must be flat. First reproduce Royer et al. (2015) on the 2003–2013 windows
   (cut/svi_ext_bostock2013). Pooled design only (per-family-per-episode impossible).
7. S6: no control windows were cut (data series only); the noise floor comes from the −15…−5 s part of each window.
   Report sizes only above a per-family-month completeness amplitude; test amplitude vs net_mad (winner's curse).
   Gain steps: NLLB BH→HH 2017-02-23 (1.77×), PGC 2015-01-12.
8. S7: define "between-episode" once (282 at no tremor that day, 118 quiet-flag, 1,037 outside all 79 episodes); 0.75 only.
9. S9: only ~5–8 distant earthquakes shook Vancouver Island enough → one pooled test, a bound.
10. S10: event population changes between episodes can mimic dv/v → mirror control; separate phase.
11. Added: SVI vs Washington-NW cross-segment comparison (lead/lag of shared episodes, leakage check for S1);
    decay of LFE rate after onset per group (Omori-like); participation map of the 40 minor episodes by depth.

Most likely to be publishable (Merlin): S2 (12-year LFE–slip scaling, if S1 passes), S3 (repeatable nucleation sites),
S7 (which families fire between episodes and timing to the next ETS). Likely nulls: S5 trend, S9, S10.

---
## Part 2 (2026-10-07 ~1:30 pm PT, user): exhaust BOTH extended catalogues (S. Vancouver Island + Washington) before dv/v
Every item ends with a plain-words verdict, notes, tree update. Replication in the second region is the priority.
A. Washington-NW repeat of S7: between-episode LFEs (deep? bursty? quiet month after each ETS?) — NW headline counted once by source (catalogue v1.3 / nw_sharing groups, dedup variant D).
B. ETS over the long record, both regions: inter-event intervals, durations, sizes (LFE count, tremor count, Tier-A GNSS where available), regularity and drift 2007/2014–2026, seasonal phase of ETS onsets and of between-episode rates (Rayleigh / permutation vs month), cross-border lead/lag (which segment starts first; delay).
C. Activity at a patch after the front passes, both regions: decay shape (Omori-like vs exponential, time constant per group/depth), burst size distribution, inter-event time statistics.
D. Washington fixed-station cross-era rates (task 3): did NW rates change 2007–2026?
E. Tides (S5) in both regions once the tide-stress setup is validated.
F. S6 LFE size (SVI; running).
G. Triggering by nearby intraslab and distant large earthquakes (bound).

---
## Part 3 (2026-10-07 ~1:30 pm PT, after Merlin's big-picture review — MERLIN_BIG_PICTURE_2026-10-07.md): the paper question
Is the deep plate interface's within-cycle behaviour (recharge time after ETS, tidal-sensitivity rise after the front,
post-front decay) the same in every cycle over 23 years (Bostock 2003–2013 + ours 2014–2026), in two segments?
Steps (cheapest invalidation first): 1 Bostock catalogue → family ids; 2 Bostock-era episodes; 3 S7 replication on Bostock
(+ WA item A); 4 recharge curve 2003–2026; 5 tides with Royer 2015 positive control; 6 strainmeter stack on deep-burst days;
7 fixed-station ratios (bound); 8 front-arrest test; 9 dv/v bound last. Detectability-free observables only for evolution claims.
