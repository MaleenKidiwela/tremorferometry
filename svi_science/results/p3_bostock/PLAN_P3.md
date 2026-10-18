# P3 — Bostock 2003–2013 baseline, quiet month after ETS, recharge curve: PLAN (frozen before any P3 test)

Written 2026-10-07 ~1:45 pm PT, before any P3 test statistic was computed. Spec: `svi_science/PLAN.md` Part 3 steps 1–4,
`svi_science/MERLIN_BIG_PICTURE_2026-10-07.md` Prescription 1–4. Scripts `svi_science/scripts/p3_*.py`; outputs here;
figures `svi_science/figures/p3_*.png`. Any later change is appended under "Revisions" with time and reason.

## What was looked at before writing this (inventory only, no test)
- Bostock catalogue (`lfe_catalog_repair/results/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv`): 269,423 rows,
  130 distinct (lat, lon, dep) keys; all rows match an `e5r_families_2012.csv` key exactly; all 115 production families
  (`svi_science/results/s0/families.csv`, fid = bkey) receive events (256,920 rows).
- Its events fall on only **230 UTC days**, in 10 runs of 11–31 days: 2003-02-26…03-10, 2004-07-02…27, 2005-09-03…25,
  2007-01-24…02-03, 2008-05-10…06-02, 2009-05-10…06-07, 2010-08-16…09-08, 2011-08-19…09-09, 2012-08-30…09-25,
  2013-09-13…10-13. Not one event outside them. Per-year rows: 2003 32,978; 2004 59,659; 2005 57,504; 2007 3,147;
  2008 14,224; 2009 22,112; 2010 18,114; 2011 18,663; 2012 28,149; 2013 14,873.
- Bostock et al. (2015, JGR, `notes/JGR Solid Earth - 2015 - Bostock - Magnitudes…pdf`, §2): "a catalog of LFE detections
  from southern Vancouver Island assembled from 10 major ETS episodes between 2003 and 2013"; 2003–2005 are the most
  complete years (dense temporary stations); templates were then applied to 2007–2013 episodes with subsets of common
  stations. **Detection ran on the ETS windows only.**
- Consequence fixed here (not a result): no between-episode day and no post-episode month was ever scanned in Bostock's
  catalogue. The S7 between-episode tests, the post-ETS quiet month and the recharge curve **cannot be measured on
  Bostock 2003–2013**. Absence of events outside the windows is "not scanned", never "quiet". This is stated as
  "not testable", not as pass or fail (Merlin's kill condition cannot be triggered or cleared by this catalogue).

## Step 1 — Bostock catalogue → our 115 family ids (`p3_bostock_map.py`)
- Key: exact (lat, lon, dep) rounded to 4 decimals → `e5r_families_2012.csv` bkey (= `svi_ref_common.bostock()` rule).
  Times = catalogue UTC (origin offsets of a few seconds are irrelevant at day level; not applied).
- PASS if all of: (a) ≥ 110 of 115 production fids receive ≥ 1 event; (b) the mapped per-year totals over the 130 keys equal
  the raw per-year row counts exactly; (c) per family per year counts equal `svi_ref_common.bostock()` for 2010, 2011, 2012
  (2013 is that module's sealed year; read here directly from the CSV, no fitting happens in P3); (d) 2012 per-family
  counts equal `e5r_families_2012.csv` n_all. Else STOP and report (headline reverts to 12 years).
- Outputs: `bostock_events_mapped.parquet` (fid, time_utc, date, dep, mag, group, downdip_km), `step1_checks.json`,
  per-year × group table.

## Step 2 — Bostock-era episodes from Bostock's own family exposure (`p3_bostock_episodes.py`)
- Exposure (Bostock era): family f exposed on UTC day d if it has ≥ 1 Bostock event that day (115 families). This replaces
  PNSN tremor (none before 2009-08; SVI PNSN exposure starts 2010-01). It is an LFE exposure, so thresholds do not mean
  the same as for tremor; reported as such.
- Rule = the final ladder rule of `episode_rates/PLAN.md` Revision 1: qualifying day = ≥ 5 % of the 115 families exposed;
  runs bridged across gaps ≤ 5 d; keep if span ≥ 5 d; MAJOR if family-exposure-days ≥ 0.15 × 115 × 14 = 241.5.
  Column: the frozen first rung (10 %, gap 3, MAJOR 0.25 → 402.5).
- PASS if each of the 10 targets (2003-02, 2004-07, 2005-09, 2007-01, 2008-05, 2009-05, 2010-08, 2011-08, 2012-09,
  2013-09) overlaps a MAJOR episode (start..end within ± 15 d of that month). Splits (> 1 MAJOR per window) reported.
  Because the catalogue only contains these windows, recovery is near-automatic; the step's value is the episode
  start/end dates and the check that the frozen rule does not break or split them.
- Cross-check 2010–2013: overlap with the PNSN-based SVI list of `b_ets_record/episodes_b.csv` (2010-08b1, 2011-08,
  2012-08, 2013-09): start/end differences in days.

## Step 3 — Quiet month after a major ETS ("S7 replication")
- **Bostock 2003–2013: not testable** (above). Stated; no time-shift surrogate is computed because there are no
  between-episode days to shift into. No scrambled control exists for Bostock's catalogue (stated).
- Instead, three episode-end definitions on the two catalogues that scanned every day (`p3_quiet_month.py`):
  - **D1 regional end** (S7/A as run): τ = days since the previous MAJOR's end.
  - **D2 end + 7 d**: τ' = τ − 7; events and eligible days with τ < 7 are dropped; F30 = share with 0 ≤ τ' < 30.
    Why: the episode rule bridges gaps ≤ 5 d, so the first ~5 days after an "end" are, by construction, days with little
    tremor; D2 asks whether the deficit lasts beyond that built-in quiet stretch.
  - **D3 the family's own end** (third definition, justified): τ = days since the last day inside that MAJOR's span on
    which the family itself had tremor within 30 km (`n30 ≥ 1`); if it had none in that MAJOR, the regional end. Why:
    the regional end is when tremor stops anywhere in the box, which can be days after slip stopped at a given patch;
    the "re-arming" clock of a patch should start when that patch stopped. It is the most local, physically motivated end
    available without using the LFE detector itself (so not circular).
  - Statistic and control exactly as S7/A: F30 = share of BE events with τ < 30 vs the family-preserving day shuffle
    (each event moved to a random eligible day of its own family, 1000×), p_low = share of shuffles ≤ observed.
    BE_broad (any local tremor) as a column.
  - Data: SVI 2014-01-01…2026-05-27 (S7 tables: P ≥ 0.75, counted once keep_xfam8, 39 majors); WA-NW E3 2015-01-01…
    2026-05-27 (A tables: P ≥ 0.75, counted once D, 69 majors). WA E1/E2 as columns.
  - Rule per dataset: **"holds under all three ends"** if F30 < shuffle median and p_low < 0.01 for D1, D2 and D3;
    **"holds only at the regional end"** if D1 passes and D2 or D3 does not; **"does not hold"** if D1 fails.
  - Status: SVI is the discovery dataset, so SVI rows are a robustness check, not a confirmation. WA's confirmatory
    verdict is A's (D1, E3: failed); D2/D3 for WA are descriptive and cannot turn that into a pass.
- Three-dataset replication status of the quiet month is then written as: SVI (discovery) / Bostock (not testable) /
  WA (A's confirmatory verdict + D2/D3 description).

## Step 4 — Recharge curve: all deep activity vs time since the last major ended (`p3_recharge.py`)
- Deep = families of S0 groups whose centroid down-dip coordinate is ≥ +10 km (groups 1, 4, 10, 12; 37 families).
  Columns: families with Bostock depth ≥ 40 km (22); group 10 alone; and the up-dip groups (centroid ≤ −10 km:
  3, 6, 8, 11) as the comparison curve.
- Events (ours): `s0/analysis_table.parquet`, run prod, data, P_same ≥ 0.75, keep_xfam8; **every** such event outside the
  reference majors is used — between-episode, minor-episode, local-tremor days alike. Exposure = observable deep
  family-days (`famday_detectability`, run prod). Scrambled control (c1) through the same steps (column). Period
  2014-01-01…2026-05-27.
- Reference lists ("major"): primary = the 39 SVI majors (`episodes_svi_final.csv`, as S7); column = the B "big" list
  (`b_ets_record/episodes_b.csv`, region SVI, big == True, 2014+, 14 episodes). Days inside a reference episode are
  excluded; t = days since the previous reference end (D1; D2 and D3 as columns), up to the next reference start.
- Curve: rate per 1,000 deep family-days in t bins 0–7, 7–14, 14–30, 30–60, 60–90, 90–150, 150–240, 240–400 d,
  normalised by the fitted plateau r∞ (shape only).
- Fit: rate(t) = r∞ · (1 − a · exp(−t / τ)), Poisson maximum likelihood on daily deep counts with exposure;
  a ∈ [−5, 1] (a < 0 = an after-ETS excess that decays), τ ∈ [1, 400] d. Recovery time constant = τ.
  95 % CI from 1000 bootstraps over cycles (each cycle = one reference end to the next reference start, resampled
  with replacement). Deficit statistic D30 = 1 − rate(0–30 d) / rate(30 d … next start).
- Control for "is there a recovery at all": time-shift surrogates — within every cycle the daily count and exposure
  series are circularly shifted by a random amount (1000×), keeping each cycle's total and its bursts but breaking the
  link to the previous end. **Recovery present** if observed D30 > 99th percentile of the surrogates (p < 0.01) and the
  bootstrap 95 % CI of a excludes 0. Otherwise "no measurable recovery", with the 99th-percentile D30 as the bound.
- **Era comparison (pre-registered):** the 2003–2013 LFE curve cannot be built (Bostock scanned ETS windows only).
  What is compared instead, same rules:
  - (E-a) ours split in two halves, 2014-01…2019-12 vs 2020-01…2026-05 (same detector; station changes flagged).
  - (E-b) deep PNSN tremor (detectability-free relative to our detector), 2010-01-22…2013-12-31 vs 2014-01-01…
    2026-05-27: epicentres within 15 km (epicentral) of a deep family; daily counts; reference = PNSN-rule majors
    (`episodes_b.csv` SVI, major == True for 2010–2013, `episodes_svi_final.csv` majors 2014+). Exposure = 1 per day.
  - Rule: **"different"** if the 95 % bootstrap CIs of τ do not overlap AND the bootstrap distribution of the D30
    difference excludes 0 at p < 0.01 (two-sided); otherwise **"no detectable difference"**, reported as "a D30
    difference larger than X (95 % CI half-width) would have been seen".
- Bostock-era substitute that IS measurable (descriptive, labelled as such, shape only): per major episode the deep share
  = deep-group events / all 115-family events, Bostock 2003–2013 (10 episodes) and ours 2014–2026 (B big list; 39 majors
  as column). Reported: within-era spread (CV), within-era trend with year and dependence on the interval since the
  previous big ETS (Spearman, permutation p, 10,000×; pass bar p < 0.01). Absolute deep shares are NOT compared across
  eras (recall vs depth differs between detectors; our recall ~0.18–0.20 of his). Caveat carried: Bostock's network
  changed between episodes (2003–2005 densest).

## Figures
- `p3_bostock_timeline.png`: Bostock daily events 2003–2013 by depth class with the Step-2 episodes, plus the scanned/
  not-scanned shading; and ours 2014–2026 for context (same depth classes).
- `p3_recharge_curves.png`: normalised deep-rate vs time since major end; ours (39 majors; big list), halves, PNSN deep
  tremor by era; surrogate band.
- `p3_quiet_month.png`: F30 vs shuffle for D1/D2/D3, SVI and WA.

## Rules carried
P ≥ 0.75 headline; rank, don't discard (no family or episode dropped; flags as columns); "scrambled control", never
"null"; dates are UTC days, times shown in PT; RAM < 30 GB, CPU only.

## Revisions
- 2026-10-07 ~1:50 pm PT (bookkeeping, no rule change; noticed in the run log after the first step-4 line had printed — the code used the group list, so only this text was wrong): the deep groups 1, 4, 10, 12
  hold **32** families (5 + 11 + 6 + 10), not 37 as written above (arithmetic slip). Up-dip groups 3, 6, 8, 11 hold 41.
- Same time (implementation detail, no rule change): the big list used as a clock includes 2013-09 (PNSN rule, B) as the
  first reset for 2014 days; for the 39-major clock the first reset is the 2014-06 major's end (days before it unused).
  For the PNSN-tremor era split, days (not cycles) are assigned to eras by date; a cycle crossing 2014-01-01 is cut there.
