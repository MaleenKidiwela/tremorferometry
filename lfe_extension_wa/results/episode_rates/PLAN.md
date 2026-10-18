# Per-family x per-ETS-episode rate consistency — PLAN (frozen before any LFE count)

Written 2026-10-06 ~8:20 pm PT, before any per-episode LFE count was made. Spec: `MERLIN_DESIGN_2026-10-06.md` (this folder).
Script: `lfe_extension_wa/scripts/episode_rates.py` (one file, sub-commands). Outputs in this folder.
Any later change to a rule below is appended in the "Revisions" section with the reason and the time, never edited in place.

Regions:
- **A. WA-NW** — the 132 NW families of `lfe_extension_wa/results/wa_prod/setup/families.csv` (`reg == "NW"`), headline label
  P_own_WAVE_PNET, fixed core network of the era (production run `conf/P/prod/`).
- **B. SVI** — the Bostock (2015) families with production candidates in `lfe_extension/results/refine/conf/P_by_candidate.parquet`
  (fam = bkey; positions = lat/lon of `lfe_catalog_repair/results/e5r_families_2012.csv`), 2014-01-01 .. 2026-09-15.

## 1. Exposure (tremor only)
- Exposure of family f on UTC day d: n_trem30(f, d) = number of PNSN tremor epicentres within 30 km (epicentral, `wa_common.km`)
  of the family position that day. Exposed = n_trem30 >= 1.
- WA: `wa_prod/post/tremor_family_day.parquet` column n30 (checked against `rates_daily_v2.parquet` n_trem30).
  Catalogue covered 2009-08-06 .. 2026-05-27 (`post_v2/catalogue_span.json`).
- SVI: built here from `catalogs/pnsn_tremor_cascadia_full.csv` (covered 2010-01-06 .. 2026-05-27) and the Bostock positions.
- Days after 2026-05-27 have no exposure: no episode can be defined there (stated, not filled).

## 2. Episodes (fixed before counting)
- Regional episode = maximal run of days on which >= 10 % of the region's families are exposed; runs separated by <= 3
  non-qualifying days are bridged (merged); keep if span (first to last qualifying day, inclusive) >= 5 d.
- MAJOR if family-exposure-days (sum over the span of the number of exposed families) >= 0.25 x N_fam x 14
  (NW N_fam = 132 -> 462; SVI N_fam = number of families in the P table, 115 expected -> 402.5).
- WA before 2009-08-06: the WECC windows of `wa_prod_post.WECC` (2007-01-12..02-02, 2008-04-30..06-02, 2009-04-28..06-03;
  the brief calls this wa_common.WECC — it lives in wa_prod_post.py) are MAJOR episodes by definition; exposure unknown
  (X = NaN, D = window length, every family's window = the whole WECC window).
- Per-family window (f, e) = union of [d - 2, d + 2] over f's exposure days d inside the episode span. A family with no
  exposure day in the episode has an empty window (N = 0, X = 0, D = 0) and is listed with `exposed = False`.
- Robustness column only: N_span = count over the whole episode span +- 2 d (same for every family).

### Acceptance check (before counting)
- "Regional entries" of the repo's reference episode tables:
  - NW: `data/extension_stations/tremor_episodes.csv`, region NW_WA: every literature/scan-days entry and every PNSN entry with
    n_strong_days >= 14 (ETS-scale; the 14 matches the 14-day scale of the MAJOR rule; small 2-6-day entries cannot be MAJOR
    by construction and are not ETS). Covered years: 2007 .. 2026-05.
  - SVI: `data/extension_stations/tremor_episodes.csv` has no SVI region. Reference = the band 48.5 N entries of the repo's other
    `fault_tomography/inversion/res_catalog_sse/tremor_episodes.csv` (all 'seg N'; onset .. peak), years 2014 .. 2026-05.
- An entry is "in the MAJOR list" if its date span overlaps a MAJOR episode span.
- If any entry fails: revise BEFORE counting, in this order, stopping at the first that passes, and record it below:
  (i) bridge gaps <= 5 d; (ii) qualifying-day fraction 10 % -> 7.5 % -> 5 %; (iii) MAJOR factor 0.25 -> 0.2 -> 0.15.
  If still failing, keep the last rule and list the failing entries as a deviation (they are not added by hand).
- MAJOR episodes that are not in the reference list are kept (reported as "not in reference").

## 3. Per-(family, episode) table (identical columns both regions)
Counts are summed over the family window days.
- WA (NW): label P_own_WAVE_PNET, core network, counts at P >= 0.5 / 0.75 / 0.85 / 0.9 for data, time-reversed (rev) and
  scrambled (c1) series, recomputed from `conf/P/prod/P_*.parquet` (identical to `rates_daily_v2` n0xx own_pnet — checked).
  N = n075_data - n075_rev (headline, main control). N_scr = n075_data - n075_c1. N050/N085/N090 = data - rev at that cut.
  sumP_net = sumP_data - sumP_rev (column).
- SVI: P_same >= 0.75 etc.; N = n_data - n_c1 (scrambled control is SVI's only control). N_scr = N (same).
  Two series: **S13 fixed** (primary: run prod_S13 for dates < 2021-04-01, run prod for later dates, whose effective set must be
  S13 — checked) and **eff** (production, effective set per day). The SVI ΣP-excess columns of famday.parquet are NOT used.
- X = sum n_trem30 over window; D = exposure days; active days = window days with >= 1 data detection at 0.75;
  onset lag = first active day - first exposure day; offset lag = last active day - last exposure day (days, signed).
- live-carrier days: WA = sum over window days of |family template stations (tpl_library) ∩ that day's loaded core set|
  (`mf/prod/run_<date>.json` members.core); SVI = same with the 2012 usable-template station set
  (`e5r_template_qc_2012.csv`, usable & template) ∩ the day's loaded stations (`refine/feat/<run>/<date>.json`).
- station-change flag: WA = any join/leave/node change of a family carrier in `post_v2/core_set_changes_v2.csv` within
  [window start - 45 d, window end + 45 d] (era-start rows excluded); SVI = any `svi_ref_common.CHANGES` station change
  touching a family carrier (or a station-set change) in the same +-45 d range; the changes are listed in a column.
- In-sample (WA): on V2008 (2008-05-01..06-18) and V2010 (2010-07-14..09-18) days, detections in build half (half_of == 0)
  and guard (-1) are dropped and the test half x F is reported (F = 86400 / (6 x 7080) = 2.034, verified empirically from
  the test-half share on out-of-sample days), flag `insample_halves`. Labels: 2007-01 "sealed", 2009-05 "semi",
  fit months (wa_prod_common.FIT_MONTHS) "fit", 2011-04 "excluded" (counted, flagged, never tested), 2012-09/10 "t2012
  templates" for families with any core candidate tpl_insample == True in that episode (expected 43-44 NW fids; checked
  against tpl_library tag t2012new). Episodes overlapping an excluded month are listed but not tested.
- SVI in-sample: none expected (truth/fit days are 2010-2013; production starts 2014) — confirmed from
  `refine/day_table.csv` roles and stated.
- Cross-era bridging column (WA only): k_v2 (E3 own_pnet, corrected `post_v2/rates_k_factors_v2.csv`, cut 075) x N; E1/E2 PNET k
  unavailable (NaN), so E1/E2 cross-era columns are net counts only.

## 4. Pre-registered sanity checks (run before step 5; reported verbatim)
- (a) WA build/test ratio = (n075 data - rev in build half) / (same, test half), pooled over out-of-sample MAJOR episodes
  (not V2008/V2010, not sealed/semi/excluded): PASS if pooled ratio in [0.85, 1.15] and every such episode with >= 400
  half-counted detections has its 95 % Poisson interval overlapping [0.85, 1.15]. In-sample episodes reported (expect ~2x).
- (b) main control <= 2 % of data (count at 0.75, regional, per MAJOR episode), both regions, every MAJOR episode
  (WA rev; SVI c1). WA c1 share reported as a column. FAIL lists the episodes.
- (c) Spearman(log regional N, log regional X) >= 0.5 across MAJOR episodes within one station era (WA E3 2015-2026; SVI S13
  fixed overall, and eff per set S15/S14/S13 where >= 4 episodes). If (c) fails for a region's primary series: STOP before
  step 5 for that region and report (detectability dominates).
- Named cases reported: SVI 2017-03, SVI 2025-10, WA 2024-08 (N, X, D, N/D, N/X vs era median).
- Merlin's NW out-of-sample per-episode reference ("net per family-tremor-day") compared with N / sum D.

## 5. Consistency (k-free, within one station era)
- Eras: WA primary E3 (2015-01-01 .. 2026-09-15); E1 (2007-11) and E2 (2012-14) as separate secondary blocks. Episodes are
  assigned to the era of their start date; an episode crossing an era boundary is listed but tested in neither.
  SVI primary = S13 fixed series, all MAJOR episodes 2014-2026; SVI eff series = secondary.
- Test set for family f: MAJOR episodes of the era where f is exposed (D_fe >= 1), not excluded, and (WA) not in-sample
  (in-sample E1/E2 episodes enter the E1/E2 secondary tests with the x2.034 test-half counts, flagged).
- Counts used in the tests: max(round(N), 0).
- Neighbourhood of f: families within 40 km (f included); if < 5 other families, the whole region.
  T_e = sum over the neighbourhood of N_ge (each over its own window).
- Primary multinomial share test: p_f = sum_e N_fe / sum_e T_e; expected N_fe = p_f T_e; chi2 = sum_e (N_fe - p_f T_e)^2 /
  (T_e p_f (1 - p_f)); df = n_e - 1 (episodes with T_e > 0); report chi2/df and a parametric-bootstrap p (2,000 draws,
  N*_fe ~ Binomial(T_e, p_f), p_f re-estimated per draw); BH q across the region's tested families (column).
- Secondary Poisson dispersion: rho_f = sum_e N_fe / sum_e X_fe; disp = sum_e (N_fe - rho_f X_fe)^2 / (rho_f X_fe) / (n_e - 1),
  episodes with X_fe > 0.
- No test when sum_e N_fe < 20 or n_e < 3 (row still listed, `tested = False`, reason column).
- Labels: "consistent" if q >= 0.05; "over-dispersed" if q < 0.05; additional label "exposure-explained" if disp < 2 while
  chi2/df > 3.
- Same tables run on the controls (WA rev and c1 counts, SVI c1) where sum >= 20.

### Artefact attribution
- Pearson residual r_fe = (N_fe - E_fe) / sqrt(T_e p_f (1 - p_f)).
- Pair correlation of residuals for family pairs with >= 5 shared tested episodes; regress corr on log10(distance km) and
  carrier-Jaccard (WA: core-station template sets of the era; SVI: 2012 usable-template sets ∩ S13); permutation p (1,000
  permutations of the family -> (position, carrier set) assignment) for each coefficient and for the Spearman of corr vs
  each predictor.
- Per (family, episode) labels, only for |r_fe| >= 2 (else "within noise"):
  - "station-coincident" if the station-change flag is set AND the mean residual of the other tested families that share a
    changed carrier has the same sign as r_fe with |mean| >= 1;
  - "spatially coherent" if the mean residual of the other tested families within 20 km (>= 3 of them) has the same sign
    with |mean| >= 1;
  - both -> "spatially coherent + station-coincident"; neither -> "isolated".
- SVI fixed-S13 ratio per (family, episode) before 2021-04: R = N_S13 / N_eff at P_same >= 0.75 (N_eff >= 5); family median
  and robust sigma (1.4826 MAD; if 0, Poisson sigma of the ratio); "set-affected" if |R - median| > 2 sigma.
- WA fixed-10 subset re-run: only if station-coincident steps show across eras; it needs a new matched-filter pass, so it
  is out of scope here (would be proposed, not run).

## 6. Pass/fail bars summary
(a) 0.85-1.15; (b) <= 2 %; (c) Spearman >= 0.5 (else stop step 5 for that region); tests need sum N >= 20 and >= 3 episodes;
labels per q < 0.05, disp < 2 & chi2/df > 3, |r| >= 2.

## Revisions
(none yet)

### Revision 1 — acceptance check (2026-10-06 ~8:25 pm PT, BEFORE any LFE count)
Frozen rule (10 %, gap 3, >= 5 d, MAJOR 0.25): NW 17/18 reference entries in MAJOR, SVI 12/13 (`acceptance_*_frozen.csv`).
Failing: NW literature entry 2009-08-05..08-25 (the PNSN 2009 pull starts 2009-08-06 and is sparse: only 2009-08-09 (68
families) and 08-16 (113) qualify; best run 08-16..08-21, 223 family-exposure-days < 462); SVI band-48.5 entry
2020-09-13..09-29 (a 4-day burst 2020-09-12..15 with 80-106 families exposed, shorter than the 5-day minimum; nothing after).
Ladder run in the frozen order (`acceptance_*_ladder.csv` holds the last rung; log lines in `episode_rates.log`):
  (i) gap 5: NW 17/18, SVI 12/13; (ii) +frac 7.5 %: 17/18, 12/13; frac 5 %: 17/18, 12/13;
  (iii) +MAJOR 0.2: 17/18, 12/13; MAJOR 0.15: NW 18/18 (passes), SVI 12/13 (2020-09 still fails: the minimum length, not
  a threshold, excludes it).
Applied as pre-registered: final rule for both regions = qualifying fraction 5 %, gaps <= 5 d bridged, span >= 5 d,
MAJOR factor 0.15 (NW 0.15 x 132 x 14 = 277.2; SVI 0.15 x 115 x 14 = 241.5). SVI 2020-09 is listed as a failing entry
(deviation; not added by hand). Because the last rung is much looser than the design's MAJOR, every episode row also
carries `in_frozen_major` (overlaps a MAJOR episode of the frozen 0.25 rule), and the consistency results are repeated on
that subset as a robustness line. The primary result uses the final (pre-registered ladder) rule.
