# SVI backward extension 2003-01 → 2013-12 — production plan

Approved by the user 2026-10-08 (GATES.md Part C, `svi_science/results/backward_gates/GATES.md`). Launched Thu 2026-10-08
1:09 am PT. "Scrambled control" = the same detector on the same day with each station site's data shifted to a
different time of that day; it runs every day in P5.

**Progress, one command:** `bash /home/jovyan/tremorferometry/svi_science/scripts/bw_prod_status.sh`

Markers and logs: this folder (`P*_DONE` / `P*_FAILED`, `*_PID`, `logs/`, `tables/`). Raw data:
`lfe_extension/traces_backward/<NET.STA>/<YYYY>/<DOY>.mseed`; caches `lfe_extension/cache20hz_backward/`. Station-days
already in the forward caches (`cache20hz_r1`, `cache20hz`, `lfe_catalog_repair/cache20hz`, `lfe_picker/cache20hz`, all
the same chain) are not downloaded again and are read from there. Nothing is deleted by any stage.

## Status (updated Sun 2026-10-11 12:15 am PT; P1–P5 rows as of their last update)

| # | stage | script | status |
|---|---|---|---|
| P1 | download | `bw_prod_p1a_v2.sh` (NRCan, 7 lanes; replaces `bw_prod_p1a.sh`) + `bw_prod_p1b.sh` (EarthScope) | queue B **DONE** (all years). Queue A: 2003 **DONE**; restarted 7:52 am PT as v2 (**PID 703856**): KELB → ETS windows ± 2 months → rest of 2004–2009 |
| P2 | 20 Hz cache | `bw_prod_p2_v2.sh` (**PID 703857**; replaces `bw_prod_p2.sh`): KELB, then each window, then each year | 2003, 2010–2013 **DONE**; rest waiting on P1A markers |
| P3 | templates | `bw_prod_p3.sh` + `bw_prod_kelb_templates.py` | **DONE**: G3 era set frozen; SNB-CMG40T era template NOT adopted (held-out recall 0.397 vs 0.399 and 0.345 vs 0.350 — no gain; `templates/p3_config.json`). Era set with KELB: `templates/kelb/` (used by every era run from now on) |
| P4 | quiet rule + thresholds | `bw_p4b_*.py`, `bw_p4b_q.sh` (P4-2b, `p4b/P4B.md`) | **P4-2b computed Sat 10-10 8:25 pm PT**: B2/B3 pass (8.75); CN4/B1 bracketed (TRIM fails); era 8.25 (cross-check pass); ladder queue PID 966687 (below) |
| P5 | matched filter | `bw_prod_p5_v2.sh` (**PID 703858**; replaces `bw_prod_p5.sh`) | **DONE**: frozen 2003, 2010, 2011, 2012, 2013; era 2003 (old site set). Next: `p5e_2003_kelb`, then each window, then each full year |
| P6 | confidence + one look | `bw_p6_fit.py`, `bw_p6_apply.py`, `bw_p6_calib.py`, `bw_p6_onelook.py` (`p6/P6.md`) | **DONE Sat 10-10 11:57 pm PT**. P for all 1.6 M candidates. Control P > 0.5: CN4 1.4 %, B1 2.2 %, B2 0.8 %, B3 1.7 %, B4 1.2 % (no refit triggered); **era 13.5 % → flagged, P\* 0.90**. One look opened 11:56 pm PT: recall at P ≥ 0.75 frozen CN4 2003 0.009, B1 2005 0.013, B2 2007 0.082, B3 2008 0.138 / 2009 0.057; era 2003 0.138, 2005 0.214. No set passes all SCOPE §5 bars → G5: ΣP bracket |
| P7 | controls + catalogue | `bw_p7_*.py` (`p7/P7.md`) | **DONE Sun 10-11 12:00 am PT**. Ladder K to S15: B2 0.68, B1 2.01, CN4 1.77, era 0.95 (era bracket [0.66, 1.59]; ETS k 3–10× quiet k). Steps (b) FAIL (placebo flags 35 %), fairness (c) FAIL (frozen), (d) mostly FAIL, (e) PASS (1.04); withheld S-energy PASS (not B1: too few); invariance FAIL except B3; G6 = control_subtracted_bracket everywhere. Catalogue `data/catalogue/svi_catalogue_v1_4_backward.parquet` + `svi_rates_v1_4_backward_*`. Summary `RESULTS_BACKWARD.md` |

### Changes of 2026-10-09 (user-approved)

1. **PO.KELB added (2003-01-01 → 2004-06-03).** KELB is the KLNB site before its rename (same coordinates to 0.0005°,
   same Guralp CMG3ESP_100sec / Nanometrics Trident). Metadata difference: overall gain 6.17e7 (2002-09-12 → 2003-07-31)
   and 3.07e8 (→ 2004-06-03) vs KLNB 3.97e8. KELB uses the KLNB era templates, placements and weights; only the
   `rel_amp` feature sees the gain (amplitude reference = KLNB's × (3.07/3.97)²; before 2003-07-31 KELB `rel_amp` is
   ~25× too small — noted, not corrected).
2. **Which era run counts:** for **2003** use **`p5e_2003_kelb`** (all 2003 era days, with KELB); `p5e_2003` (before KELB)
   is kept on disk and superseded. For 2004–2006 the era runs `p5e_<Y>` include KELB from the start. Adding KELB
   lengthens the scrambled-control site list, so the per-day site shifts of era runs differ from those of the old
   `p5e_2003` (statistically equivalent; the frozen runs keep the forward site list).
3. **NRCan order:** for the years not yet downloaded, the Bostock ETS windows ± 2 months come first
   (W2004 2004-05-02→09-27, W2005 2005-07-03→11-25, W2007 2006-11-24→2007-04-03, W2008 2008-03-10→08-02,
   W2009 2009-03-10→08-07), then the rest year by year. P2 and P5 follow per window (markers `P1A_W*`, `P2_W*`,
   `P5W_*`); P5 uses the same run names (`p5f_<Y>`, `p5e_<Y>`) for window and full-year passes, so days are never run twice.
4. **Repair passes no longer re-ask days with no file:** a day whose first pass failed all 20 tries on all three channels is a
   real gap (if the data existed, the chance of 60 failures at NRCan's ~1/6 success rate is ~1e-5). In 2003 the old
   repair passes spent ~14 h on such days, mostly **NLLB, missing for all of 2003** (363 days), plus 15–16 days each of VGZ
   and LZB. The 2004 pass stopped (by PID) at day 729 of 3,871; finished days are kept and skipped.

## P4 — the decision that is needed (not decided here)

The forward thresholds (T_core 8.75 / T_ext 8.25 and their per-epoch re-derivations) come from quiet-day false alarms,
with "quiet" = no PNSN tremor in the family box that day and ≥ 3 days from any. PNSN tremor starts 2009-08 (local
catalogue 2010-01), so before 2010 another quiet rule or another threshold rule is needed. Options:

1. **Catalogue distance:** quiet = ≥ 30 days from any Bostock ETS window. Simple and independent of our detector, but
   his catalogue covers only the ten big ETS windows, so minor inter-ETS tremor counts as "quiet" → thresholds come out
   higher (conservative).
2. **Array-energy flag** (`lfe_array_discovery` `w_energy_hours.py`, Merlin's Q3 fallback): hourly tremor-band energy on
   our own stations. Data-driven, but its sensitivity changes with the station set, i.e. exactly at the era boundaries.
3. **Both:** quiet = option 1 AND not flagged by option 2. Strictest; fewest quiet days (needs ≥ 200 quiet family-days
   per network, SCOPE G4).
4. **No quiet days for the threshold at all:** set T per network era so that its scrambled-control rate per observable
   family-day equals the forward E0 rate at 8.75 / 8.25 (the matching used in G1/G3); quiet days are then needed only
   for the between-episode analyses, where option 1 or 3 is used.
5. Wech/WTD tremor catalogue 2005–2014: would be the natural source, but was not fetchable from here (2026-09).

**Recommendation (for the user / Merlin to confirm):** option 4 for the thresholds — it uses the same rule in every
era, needs no external catalogue, and the anchor check is that 2010–2013 (B4 = S15) reproduces the forward 8.75 / 8.25
within one grid step — plus option 3 as the quiet-day flag for the between-episode science, with PNSN replacing option 2
from 2010. P5 stores every candidate ≥ 8 MAD, so no matched-filter work has to wait for this decision.

## ETA summary (PT, estimated Fri 2026-10-09 7:55 am from the measured NRCan rate ≈ 6 station-days/min)

- KELB downloaded + cached ~10 am Fri; `p5e_2003_kelb` done ~noon Fri.
- ETS windows (≈ 4.5k NRCan station-days) downloaded ~11 pm Fri; windows through P5 ~early Sat 10-10.
- All P1 (≈ 13.5k NRCan station-days left) ~**Sun 2026-10-11 late evening**; P5 complete ~**Mon 2026-10-12 morning**.
- P6/P7 after the P4 decision, ~1 day. Storage now: raw 0.83 TB, caches 0.55 TB; 58 TB free.

## P4-2b — Merlin review 2 (2026-10-11; adopted, user said "keep going")
Root cause of "2003–05 thresholds ~10× too strict": 25 PGC burst days (18 in Feb 2015, after the 2015-01-12 PGC change) carry 85 % of the CN4 scrambled-control candidates in the P4-2 reference days. The pooled-mean rule is hijacked by them. With burst-robust statistics every frozen set lands at 8.75 ± 0.15 (median-matched CN4 8.80, B1 8.70, B2 8.85, B3 8.80). Real 2003–05 CN4 control at 8.75 is 0.076–0.106, the same as S15. The 2005-11 "step" is the same artefact.
Second flaw (structural): in the frozen detector the merged PGC site (PGC+B009/B010/B011) weighs ~2× any other site and is shifted as one unit in the control. So for CN4/B1/B2/B3 the scrambled control is a one-site-signal model, and data − control is biased low (more so on tremor days).
Prescription:
1. P4-2b: same 509 days. (i) per-day median control rate; (ii) pooled mean after dropping the top-5 % S15fix control days. Reference S15fix at 8.75. Pass: every frozen set within [8.50, 9.00] by both → T_core 8.75 all frozen sets, T_ext 8.00. Statistics disagreeing > 0.5 → bracket.
2. Era detector: same statistics on its env-quiet days (quiet_env3_CN4), matched to the S15fix levels (expected ~8.25). Cross-check frozen CN4 on the same days → 8.75 ± 0.25, else bracket. Use the control, never data FA.
3. Burst-day column (control ≥ 5× yearly median); flag, never drop.
4. Single-site column n_z2 ≥ 2. B2/B3 between rates = bracket [all, n_z2 ≥ 2]. Per-site-cap test logged as a decision point after P6, on the 2005-11..2006-06 overlap.
5. Completeness ladder: old set run 60 d after each boundary (2005-07-01, 2005-11-02, 2007-11-13, 2009-08-27); k from env-quiet ΣP-excess with bootstrap CI; k_CN4/era etc. from the overlap, ETS vs between separately; agree within ×2, else bracket.
6. ΣP science threshold-free: keep all ≥ 8 MAD; rate = (ΣP_data − ΣP_control)/observable family-days; T only for the headline count; never T before P.
7. Step checks before any 2003–13 vs 2014–26 claim: (a) 1–2 pass; (b) no step > 2σ after k in ±60 d; (c) slope of excess on control ≈ 0 (env-quiet days); (d) Bostock ETS recall at P ≥ 0.75 agrees with ladder k within ×1.5; (e) 2010–13 between rates on env-quiet days (N ≤ 12).
Keep the pooled-mean values beside the robust ones; this is a recorded protocol change.

### P4-2b — status and adopted values (Sat 2026-10-10, 8:25 pm PT; details `p4b/P4B.md`, rules written before computing)

| detector | T_core MED | T_core TRIM | T_core POOLED (P4-2, kept) | outcome |
|---|---|---|---|---|
| frozen CN4 (to 2005-06) | 8.79 | 9.82 | 10.76 | TRIM outside [8.50, 9.00], MED−TRIM 1.03 > 0.5 → **bracket [8.79, 9.82]** |
| frozen B1 (2005-07..10) | 8.69 | 9.56 | 10.38 | same → **bracket [8.69, 9.56]** |
| frozen B2 (2005-11..2007-11) | 8.84 | 8.90 | 8.91 | pass → **8.75** |
| frozen B3 (2007-11..2009-08) | 8.76 | 8.77 | 8.79 | pass → **8.75** |
| S15 (from 2009-08) | 8.72 | 8.75 | 8.75 | reference → **8.75** |
| era (2003-02..2006-06), env-quiet days | 8.25 | 8.24 | 8.07 | frozen-CN4 cross-check 8.57 / 8.57 (within 8.75 ± 0.25) → **8.25** |

- T_ext 8.00 for every set (as prescribed; era T_ext also 8.00 by both statistics).
- Rule 1 ("every frozen set within [8.50, 9.00] by both") **fails** on CN4/B1 TRIM. The prescribed trim drops the top-5 %
  **S15fix** control days, which contain only 15 of CN4's 25 burst days. The burst-day claim itself is verified: 85.4 % of
  CN4 control ≥ 8.75 falls on 25 days, 18 of them in Feb 2015 and 23 after 2015-01-12, and the PGC site is the top site for
  87 % of them. Diagnostic only, not adopted: trimming by each set's own top-5 % gives CN4 9.03, B1 8.88. CN4/B1 carry
  the bracket until the user / Merlin decide. T matters only for the headline count (rule 6).
- Flag columns: sidecars `p4b/sidecar/<run>_flags.parquet` (`burst_day`, `multi_site_z2`, keyed run/date/row) for every
  finished P5 run. p5f_2008/2009 and the ladders follow from the queue.
- Ladder: `lad_CN4_20050701` done. `lad_B1_20051102` was running at 8:25 pm PT. `lad_B2_20071113` 2008 part waits for
  `P2_2008_DONE` (~Sun 10-11) and `lad_B3_20090827` waits for `P2_2009_DONE` (~Sun 10-11 late / Mon 10-12 early). The k script
  `bw_p4b_kladder.py` is ready and waits for P from P6 (smoke-tested only).

### Decision-point log
- **D-P6-1 (open, 2026-10-11): era P recalibration.** The era model (fitted on the 2004 ETS build blocks, 76 % positive)
  gives the scrambled control P > 0.5 in 13.5 % of candidates and ranks no better than net_mad on 2003/2005. Pre-registered
  outcome: flagged, P\* 0.90 for singles, ΣP excess for rates. Option: refit with quiet-day control rows as negatives (new
  pre-registered step). User / Merlin to decide.
- **D-P7-1 (open, 2026-10-11): step-check statistic.** The day-bootstrap z flags 35 % of no-change month starts (placebo);
  a 2003–13 vs 2014–26 between-episode claim needs a test calibrated on the placebo distribution.
- **D-P4b-1 (open; P now exists, still not run — needs a capped-weight GPU matched-filter run): per-site-cap test.** In the frozen CN4/B1/B2/B3 detectors the merged PGC site weighs ~2× any
  other site and is shifted as one unit in the scrambled control. Test after P6, on the 2005-11..2006-06 era/frozen
  overlap: cap any one site's share of the network sum and compare the ΣP excess and k_era/B2 with and without the cap.
  Decide then whether the cap enters the frozen B2/B3 between-episode rates. Until then, B2/B3 between rates are the bracket
  [all, n_z2 ≥ 2].
- **D-P4b-2 (open, now): CN4/B1 T_core.** The P4-2b TRIM rule failed for CN4/B1, so they are bracketed ([8.79, 9.82] and
  [8.69, 9.56]). Options: keep the bracket, adopt MED (8.75), or re-register TRIM with the set's own top-5 % days
  (diagnostic 9.03 / 8.88). User / Merlin to choose. Nothing downstream needs it before the headline count.


### P4-2b decision (main session, 2026-10-10 8:24 pm PT)
CN4/B1 keep the pre-registered BRACKET (MED vs TRIM differ > 0.5): headline counts for 2003-01..2005-10 frozen detector reported at both ends. Reason: the rule was fixed before computing; TRIM fails only because the S15fix-chosen trim misses 10 of CN4's 25 burst days, which favours MED, but choosing MED after seeing this would be a second post-hoc change. Cost is small: rates are ΣP-based (threshold-free), and the era detector (T_core 8.25, cross-check PASS) is the primary detector for 2003-02..2006-06. B2/B3 = 8.75, era = 8.25, T_ext 8.00.
