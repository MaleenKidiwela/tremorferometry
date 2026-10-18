# SVI backward extension 2003-01 → 2013-12 — go/no-go gates

Scope and plan: `svi_science/results/backward_scope/SCOPE.md`. Scripts: `svi_science/scripts/bw_gate_*.py`.
"Scrambled control" = the same detector run on the same day, with each station site's data shifted to a
different time of that day.

## Verdicts in plain words (written 2026-10-07, 10:15 pm PT; rules in Part A were fixed before computing)

- **G0 download pilot — PASS.** NRCan does deliver the pre-2010 CN data. Every CN station-day of the pilot
  (264/264) and of the G2 windows arrived complete (3 components), but only after many tries: a single NRCan
  request succeeds 16 % of the time, so a station-day costs ~19 requests and ~70 s per lane. A second pass that
  re-asks only for the missing channel recovered all 42 incomplete days. The detector and its scrambled
  control ran on all 66 pilot days with 0 station-day errors; on 2008 ETS days it finds ~3× more candidates
  than its scrambled control, on quiet days ~2×.
- **G1 reduced networks — PASS.** On the 2011–2012 truth days, the 4 CN stations alone keep 77–85 % of the
  full network's recall at the same scrambled-control rate (pass bar 40 %); 112 of the 116 templated families
  stay on ≥ 2 sites; 2010 (out of sample) 68–80 %. The 2005-11 network (13 stations) keeps 83–85 % (bar 80 %).
  The cost of fewer stations shows up as lower precision (CN4 0.30–0.36 vs 0.38–0.48 at matched control rate)
  and in the confidence model: its scrambled control gets P > 0.5 in 4.9 % (2011–12) and 6.4 % (2010) of
  candidates, above the SCOPE bar of 3 % (full set with the same code: 3.1 %). So CN4-era between-episode rates
  would need a stricter P cut or the bracket treatment.
- **G2 old VGZ / SNB instruments — PASS for 2006-04 → 2011, MARGINAL before.** VGZ (CMG3T, 2007–09) and SNB
  (GD2_v12 2006–08; GD1 2008–11) transfer as well as the never-changed stations (2008 and 2010 windows).
  For 2003–2006, the test first had to be repaired: **Bostock's 2004 and 2005 event times are truncated to
  whole seconds** (every `sec` value is an integer in those two years only), which smears every stack, also at
  the never-changed PGC. With his events re-timed onto our axis by the other CN stations, VGZ (CMG3ESP/GD1,
  2003–07) passes; SNB (CMG40T, before 2006-04-20) misses the pre-registered bar (cc ratio 0.95 vs 1.03 for the
  no-change stations; median cc 0.69 vs 0.73 in 2012). The pre-registered rebuild check cannot decide it: the
  never-changed NLLB fails that check too. So SNB before 2006-04-20 is a small, possible loss, not a broken
  template; it is tested as an A/B in production stage S3 (below).
- **G3 POLARIS era templates — PASS, by a large margin.** Templates for 11 POLARIS stations + LZB, built from
  the 2004 ETS build blocks (1,105 usable family-station templates, 125 families), raise held-out recall of
  Bostock's 2004 events from 0.035–0.046 (CN4 frozen) to 0.42–0.43 at an equal or lower scrambled-control
  rate (×9–12; 90 % interval of the gain +0.34…+0.44), precision 0.56–0.67 → 0.73–0.74, rediscovered families
  57–67 → 118. Two cautions: these templates were scored on the same episode they were built from (other
  blocks), and Bostock's catalogue was itself made on these stations, so part of the gain is "seeing what he
  saw". The CN4-only recall in 2004 (0.035–0.046) is far below its 2010–12 value (~0.13–0.16): his 2004
  catalogue is much denser, so CN4 alone would badly under-count 2003–2005 relative to his catalogue.
- **Go.** All four gates pass (G2 with one marginal station-epoch). Production should run the era-augmented
  detector for 2003-02 → 2006-06 and the frozen network after, with the plan in Part C. Nothing large has been
  downloaded and no production stage has been started.

## Part A — pass rules, written BEFORE anything was computed (2026-10-07, 7:20 pm PT)

Everything below reuses the forward SVI detector unchanged (`lfe_extension/scripts/svi_ref_mf.py`: frozen
2012 templates, weights, moveout, mergeB sites, ≥ 8 × MAD candidates, one scrambled control per day). New
station sets and the pre-2010 station epochs are added from a wrapper (`bw_gate_mf.py`), never by editing the
forward scripts. 2013 truth is NOT used anywhere (sealed; its one look is spent).

### G1 — can the 4 CN stations (2003-01 → 2005-06) and the 2005-11 network detect the families? (no download)

Runs: the forward truth days already on disk (`lfe_extension/results/refine/day_table.csv`), role `fit`
(2011 inside Bostock's window, 2012 test blocks only) and role `s9` (2010, reported).
- **CN4** = PGC, VGZ, NLLB, SNB (2011, 2012, 2010).
- **B2** = the 2005-11 set = CN4 + B001, B003, B004, B005, B006, B007 + B009, B010, B011 (13; 2011, 2012).
- **B1** = CN4 + Olympic PB (B001, B003–B007) (10; 2011, 2012), extra, for 2005-07..10.
- Reference = the forward full-set runs on the same days (S15 for 2011/2012: `r3t2011_nc`, `r3t2012_nc`;
  S9 for 2010: `r3t2010_nc`), read-only from `lfe_extension/results/refine/feat/`.

Scoring (exactly `svi_ref_r3.py`): one-to-one recall of Bostock's events (own family, ± 2 s) by all candidates
≥ T, where T for the reduced set is chosen so that its scrambled-control rate (control candidates per
observable family-day) MATCHES the full set's rate at the forward operating thresholds 8.25 and 8.75 MAD.
Precision = fraction of data candidates ≥ T that match one of his events of the same family (± 2 s).
Families "detectable on ≥ 2 sites" = families observable (≥ 2 live sites) on ≥ 90 % of the scored days.
"Rediscovered" family = ≥ 5 of his events matched at the matched rate (r3 definition).

**G1 pass rule, CN4** (decides whether 2003-01 → 2005-10 can be used between episodes):
1. pooled 2011 + 2012 recall of CN4 ≥ 0.40 × the full set's recall, at the matched scrambled-control rate, at
   BOTH operating rates (8.25 and 8.75); AND
2. ≥ 100 of 115 families are observable on ≥ 2 sites on ≥ 90 % of the scored days, AND the number of
   rediscovered families is ≥ 50 % of the full set's.
- Fail → 2003-01 → 2005-10 is ETS-window-only (and the era-augmented LZB/POLARIS variant is the fallback).
- Reported beside it, not part of the pass (the SCOPE's model bar, for the later confidence stage): an
  out-of-fold WAVE confidence model fitted on the CN4 rows exactly like the forward one (2011 + 2012 test
  blocks, 5 folds grouped by day, HistGradientBoosting + isotonic): AUC(same), calibration slope,
  matched/ΣP, control P > 0.5. SCOPE bar: AUC(same) ≥ 0.80 and control P > 0.5 ≤ 3 %. The same code is run
  on the full-set rows as the reference. 2010 is scored once by the CN4 model, out of sample.

**G1 pass rule, B2** (decides tier 2, 2005-11 → 2009): pooled 2011 + 2012 recall ≥ 0.80 × the full set's at
the matched rate, at both operating rates (SCOPE S1 bar).

### G2 — do the 2012 templates still fit VGZ and SNB on their older instruments? (≤ 20 GB download)

Track A-2 step-4 protocol (`lfe_extension/scripts/a4_epoch_stack.py` functions, unchanged: cut 3-C 20 Hz
1–8 Hz windows at the predicted S, rms-normalise, cc-align, stack; signed 3-C cc against the 2012 stack over
S−1..S+7 s maximised over ± 2 s → cc and S shift dS; registration = median leave-one-out per-event cc).
One change, stated up front: the reference events are **Bostock's own catalogue events** of that ETS (truth
time + frozen 2012 per-family offset), not withheld-station detections. Bostock's catalogue exists for these
windows and is independent of the station under test, so it is the cleaner reference.
- New leg: the old-epoch ETS window. Control leg: the same station on 2012 Bostock events in the TEST blocks
  only (the 2012 templates were built on the build blocks), each family capped at the new leg's event
  count (matched N). Families need ≥ 20 events; ≤ 400 per family (largest magnitude first).
- Epochs and windows:
  | station-epoch | sensor / digitiser | test window(s) |
  |---|---|---|
  | VGZ BH 2007-02-28 → 2009-07-30 | CMG3T / GD2_v21 | 2008-05-10 → 06-02 ETS |
  | VGZ BH 2003-02-27 → 2007-02-28 | CMG3ESP / GD1 | 2004-07-02 → 07-27, 2005-09-03 → 09-25 |
  | SNB BH 2008-12-11 → 2011-03-17 | CMG3ESP / GD1 | 2010-08-16 → 09-08 (on disk) |
  | SNB BH 2006-04-20 → 2008-12-11 | CMG3ESP / GD2_v12 | 2008-05-10 → 06-02 |
  | SNB BH → 2006-04-20 | CMG40T / GD1 | 2004-07, 2005-09 |
- No-change controls (same instrument as in 2012, same windows, same procedure): PGC (one epoch 2000–2015)
  and NLLB (CMG3ESP/GD2_v12 in both 2003–2007-01 and 2007-02 → 2017) in every window; PB B009–B011, B926 and
  B001/B003–B007 where they existed.
- Per station-epoch: cc ratio = median(cc new) / median(cc control leg); registration ratio = median over
  families of registration new / registration control leg; dS median and MAD; polarity flips (cc ≤ −0.6).
  The no-change reference = the same two ratios pooled over the no-change stations of that window.

**G2 pass rule** per (station, epoch), applied to each window and to the windows pooled:
cc ratio ≥ no-change cc ratio − 0.05 AND registration ratio ≥ no-change registration ratio − 0.05 AND
|median dS| ≤ 0.10 s AND polarity flips ≤ 10 % of families. (The 0.05 tolerance is the spread of the three
no-change-like transfers in A-2: 0.84–0.90.)
- Where an epoch fails: rebuild its templates from Bostock's events of that epoch in the BUILD (even) 2-h
  blocks, and validate on the TEST (odd) blocks: cc of the rebuilt template with the test-block stack, and
  test-block registration, compared with the same numbers for the 2012 templates on the 2012 build/test split.
  Rebuild validates if both are ≥ the 2012 values − 0.05. Gate passes for that epoch if transfer passes or the
  rebuild validates.

### G0 — can the pre-2010 data be downloaded and run end to end? (pilot)

Pilot = all template stations alive at the time, **2008-05-01 → 2008-06-04** (covers the whole 2008 ETS
05-10 → 06-02) and **2008-08-01 → 08-31** (quiet). CN from NRCan
(`https://www.earthquakescanada.nrcan.gc.ca/fdsnws`), one request per channel, ≥ 10 tries per channel with
backoff; PB from EarthScope. The G2 windows 2004-07 and 2005-09 go through the same downloader and are also
reported. Downloads go to `lfe_extension/traces_backward/<NET.STA>/<year>/<doy>.mseed`, 20 Hz caches to
`lfe_extension/cache20hz_backward/` (forward chain `a4_cache20.py`).

**G0 pass rule:** for each CN station, 3-component days with ≥ 50 % coverage ≥ 95 % of its metadata-epoch days
in the pilot (pooled over CN), AND the cache + detector (+ scrambled control) runs on every pilot day with
≤ 3 % station-day errors (the forward guard). Reported, not part of the pass: candidates ≥ 8.75 MAD per
day, data vs scrambled control, ETS days vs quiet days (B3 = forward S14 set). The 2008 detector recall is
deliberately NOT scored here: 2003–2009 truth is kept unspent for the one look per network model (SCOPE R7).
If CN < 80 %: stop and ask about NRCan bulk / AutoDRM (SCOPE S2).

### G3 — POLARIS era templates (added by the coordinator 2026-10-07; rules written 8:58 pm PT, before any G3 data)

Bostock built his 2003–2013 catalogue on POLARIS + LZB; the frozen templates cover only stations running in
2010–13, so 2003–2006 would otherwise run on CN4 alone. G1/G2 are not redone.
- Data: 2004-07-01 → 07-28 (the 2004 ETS, 55k catalogue events) + 2004-06-04 → 06-17 (14 quiet days before the
  episode; after the KELB → KLNB rename of 2004-06-03, so one station code). Stations: C8 SHVB, MGCB, TWBB,
  SHDB, KHVB (EarthScope first, NRCan fallback); PO TWKB, TSJB, KLNB, SILB, SSIB, TWGB (NRCan, ≥ 10 tries per
  channel); CN LZB (NRCan). CN4 quiet days added. Budget ≤ 20 GB on top of G2.
- Era templates: for every (family, new station), stack at Bostock's 2004-07 event times in the BUILD (even)
  2-h blocks only, exactly the forward e5r recipe (predicted S, ± 20 s, rms-normalised, cc-aligned, S pick,
  12-s template from S − 1 s, 3-C L2-normalised; usable if stack S-SNR ≥ same-N random-time stack + 3 dB;
  ≥ 30 build events; weight 10^(SNR/20) median-normalised, capped at 3; amplitude reference = median 12-s
  window energy on the quiet days). Time axis: Bostock's 2004 times sit ~0.3–0.6 s off the frozen axis
  (`bw_gate_g2_lagdiag.py`), so each family's event times are first shifted by the median per-event lag of
  the never-changed PGC frozen template at PGC (build-block events). Each new station is its own site.
- Detector runs (forward `svi_ref_mf.py` via the wrapper, same 40 days, same scrambled-control site list
  for both): **(a) CN4 frozen templates**; **(b) CN4 + POLARIS/LZB era templates**.
- Scoring on the TEST (odd) blocks of the 2004-07 episode only: one-to-one recall of Bostock's events (own
  family, ± 2 s), precision, at T_a = 8.25 / 8.75 for (a) and the T for (b) that matches (a)'s scrambled-control
  rate per observable family-day; quiet-day data and control rates reported. Per family: which stations
  carry it (frozen and era).
- Note on truth use: this spends the 2004-07 test-block look for the CN4 (B0) detector; the later B0
  one-look (SCOPE S7) must use 2003 + 2005-09 and report 2004 as already seen.

**G3 pass rule:** the era-augmented detector (b) raises held-out recall at the matched scrambled-control rate,
at both operating rates: recall(b) / recall(a) ≥ 1.10 AND the day-bootstrap 90 % interval of recall(b) −
recall(a) is above 0 AND the number of rediscovered families (≥ 5 matched) does not fall. Pass → production
uses the era-augmented detector for 2003 → 2006 (while POLARIS runs) and logs the detector change at each era
boundary (`net_set`, `tpl_epoch` columns). Fail → CN4 frozen templates in 2003–2005-10.


## Part B — numbers

### G0 (`tables/g0_decision.json`, `g0_download_by_station.csv`, `g0_requests.csv`, `g0_detector_days.csv`, `figures/g0_pilot.png`)

| item | value |
|---|---|
| CN pilot station-days, 3-C ≥ 50 % coverage / metadata-epoch days | **264 / 264 (100 %)**; G2 windows 2004-07, 2005-09: 212 / 212 |
| NRCan success per single request | 0.16 (≈ 19 requests and ≈ 70 s per station-day per lane, 8 lanes) |
| channels still missing after 20 tries (first pass) → recovered by the repair pass | 43 → 42 of 42 days complete |
| PB (EarthScope) pilot quiet month | 293 / 310 3-C days; B005 has 17 real gaps (204 at EarthScope) |
| detector + scrambled control on pilot days | 66 / 66 days, 901 station-days ok, 0 errors (guard 0.00 %) |
| candidates ≥ 8.75 MAD per day, data / control | ETS 273 / 91; quiet 25 / 13 |
| cache chain identity | re-cached PB.B001 2008-05-15 is sample-identical to the forward cache (max diff 0) |
| sizes | raw 13.1 GB (G0 pilot 8.7 + G2 windows 4.4), 20 Hz caches ~21 GB |

### G1 (`tables/g1_recall.csv`, `g1_decision.json`, `g1_model.csv`, `g1_observable.csv`, `g1_precision_same_any.csv`, `figures/g1_recall.png`)

Recall of Bostock's events, reduced set vs full set at the full set's scrambled-control rate (2011 fit days + 2012 test blocks):

| set | T_ref | 2010 (vs S9) | 2011 | 2012 | pooled 2011–12 ratio | bar | precision full → set (pooled) |
|---|---|---|---|---|---|---|---|
| CN4 | 8.25 | 0.205 → 0.139 | 0.200 → 0.157 | 0.184 → 0.137 | **0.77** | 0.40 | 0.38 → 0.30 |
| CN4 | 8.75 | 0.166 → 0.134 | 0.161 → 0.139 | 0.146 → 0.121 | **0.85** | 0.40 | 0.48 → 0.36 |
| B2 (2005-11 set) | 8.25 | – | 0.200 → 0.174 | 0.184 → 0.149 | **0.85** | 0.80 | 0.38 → 0.35 |
| B2 | 8.75 | – | 0.161 → 0.138 | 0.146 → 0.116 | **0.83** | 0.80 | 0.48 → 0.45 |
| B1 (CN4 + Olympic PB) | 8.25 / 8.75 | – | | | 0.71 / 0.83 | 0.40 | |

- CN4: families observable on ≥ 2 sites on ≥ 90 % of days 112 / 116 templated; rediscovered families 157 vs 175
  (8.25) and 152 vs 168 (8.75), summed over the two years (≥ 50 % bar met).
- At 8.25 the CN4 matched threshold sits at the 8.0 floor of the candidate list, so CN4 is compared at a LOWER
  control rate than the full set (0.93–1.43 vs 1.05–1.68 per family-day): its recall ratio there is conservative.
- Confidence model (same forward recipe, out of fold; `g1_model.csv`): AUC(same) S15 0.869, CN4 0.870 (2010 out
  of sample 0.835), B2 0.854; calibration slope 1.01–1.04; control P > 0.5: S15 3.1 %, **CN4 4.9 % (2010: 6.4 %)**,
  B2 2.6 %, B1 5.8 %. SCOPE bar (≤ 3 %): B2 passes, CN4 and B1 do not.
- Precision drop of CN4 is not only neighbour-family confusion: any-family precision also falls 0.67–0.70 → 0.54–0.55.

### G2 (`tables/g2_*.csv`, `g2w_*.csv`, `g2r_*.csv`, `g3_delta_pgc_check_*.csv`; `figures/g2_transfer.png`, `g2w_transfer.png`, `g2r_transfer.png`)

Three reference-event versions were run (only the first was pre-registered; see Deviations):

| station-epoch | window | Bostock times (pre-reg.) cc / reg ratio vs no-change | withheld detections | Bostock re-timed (best powered) |
|---|---|---|---|---|
| VGZ CMG3T/GD2_v21 2007-02 → 2009-07 | 2008-05 | 0.997 / 1.27 vs 0.919 / 1.17 **PASS** | 1.17 / 1.28 vs 1.19 / 1.13 PASS | 1.07 / 1.35 vs 1.11 / 1.15 **PASS** |
| SNB CMG3ESP/GD2_v12 2006-04 → 2008-12 | 2008-05 | 1.020 / 1.18 vs 0.919 / 1.17 **PASS** | 1.31 / 1.47 vs 1.19 / 1.13 PASS | 1.058 / 1.20 vs 1.111 / 1.15 fail by 0.003 |
| SNB CMG3ESP/GD1 2008-12 → 2011-03 | 2010-08 | 0.986 / 1.06 vs 1.023 / 1.08 **PASS** | – | – |
| VGZ CMG3ESP/GD1 2003-02 → 2007-02 | 2004-07 + 2005-09 | 0.38 / 0.65 vs 0.42 / 0.74, dS +0.35 s FAIL (no-change stations fail alike) | 1.06 / 1.14 vs 1.14 / 1.46 FAIL (n 31) | 0.995 / 1.07 vs 1.034 / 1.02 **PASS** (median cc 0.575 vs 0.578 in 2012) |
| SNB CMG40T/GD1 → 2006-04 | 2004-07 + 2005-09 | 0.30 / 0.82 vs 0.42 / 0.74 FAIL (same) | 1.64 / 1.41 vs 1.14 / 1.46 PASS (n 57) | **0.951** / 1.09 vs 1.034 / 1.02 **FAIL** (median cc 0.69 vs 0.73) |

- dS = 0.00 s and polarity flips ≤ 5 % in every re-timed and withheld test; no timing shift from any old instrument.
- Why the pre-registered version fails in 2004/2005 for every station: at PGC, the frozen template sits on average
  +0.45 s after Bostock's 2004 times (family-to-family MAD 0.08 s) with 0.30 s scatter inside a family; in 2008 and
  2012 the same measurement gives +0.00 s and 0.05 s (`g3_delta_pgc_check_*.csv`). The catalogue's 2004 and 2005
  seconds are all integers (2003 and 2007–2013 have 0.05-s resolution): truncation gives exactly a +0.5 s mean and
  ±0.5 s spread. Matching with ± 2 s is unaffected; stacking is not.
- Rebuild check (pre-registered; build-block stack vs test-block stack, against the 2012 split; `g2r_rebuild_summary.csv`):
  fails for VGZ, SNB **and the never-changed NLLB** (0.29–0.33 vs 0.50) while PGC passes, so it measures the smaller,
  noisier 2004/05 events, not the instrument; it cannot decide SNB-40T.

### G3 (`tables/g3_recall.csv`, `g3_decision.json`, `g3_quiet.csv`, `g3_family_recall.csv`, `g3_carriers.csv`, `g3_template_qc.csv`; `figures/g3_recall.png`)

| | (a) CN4 frozen | (b) CN4 + POLARIS/LZB era |
|---|---|---|
| T (b matched to a's control rate) | 8.25 / 8.75 | 8.00 / 8.05 (8.00 = list floor; b's control rate 0.10 vs a's 0.22) |
| held-out recall, 27,555 test-block events | 0.046 / 0.035 | **0.427 / 0.422** (×9.2 / ×12.2; 90 % CI of gain +0.34…+0.44) |
| precision (own family ± 2 s) | 0.56 / 0.67 | 0.73 / 0.74 |
| rediscovered families (≥ 5) | 67 / 57 | 118 / 118 |
| quiet days 2004-06-04..17, data / control per family-day | 0.29 / 0.12 (8.25), 0.11 / 0.05 (8.75) | 0.12 / 0.04, 0.10 / 0.04 |
| families on ≥ 2 sites; median carriers | 112; 4 | 126; 13 |

- Era templates usable per station (of 126 families with ≥ 30 build events): KLNB 113, SSIB 108, TWKB 107, SILB 103,
  MGCB 101, LZB 99, SHDB 98, SHVB 83, TWBB 80, TSJB 73, KHVB 71, TWGB 69. Per-family carriers: `g3_carriers.csv`
  (`frozen_cn4`, `era` columns).
- Data: 9.7 GB raw (12 era stations + CN4 quiet days); every NRCan day complete after the repair pass except SHDB
  2004-06-06 (204 on all 60 tries: a real gap).

### Deviations from Part A (all decided after a result, stated here)

1. **G2 reference events.** The pre-registered Bostock-time version could not test 2004/2005 because his times are
   truncated there (above). Two added versions: Track A-2's original withheld-station detections (few families per
   station, noisy) and Bostock's events re-timed onto the frozen axis by the other CN stations (`bw_gate_retime.py`,
   ± 0.6 s around his time + 0.5 s). The re-timed version is the best powered and is the one the verdict uses.
2. **G3 time axis.** Part A said "shift each family by the PGC lag"; the per-event re-timing above replaced it
   (truncation is per event, not per family). The per-family PGC lag is still written (`g3_delta_pgc.csv`).
3. **G1 family count.** The first G1 observability count looped over keys 0–114; the family table has 130 keys (116
   templated). Fixed before the verdict (112 / 116); the recall numbers were never affected.
4. **G0** includes a channel-repair pass (re-ask only for the missing channel) that was not in the plan; 2008 truth
   was not scored, as planned.
5. Housekeeping: when the first download queue was stopped (by PID) to shorten time-outs, `rm -f` was run on
   `traces_backward/*/*/*.tmp` (incomplete writes of that killed run, if any). Nothing else was deleted.

## Part C — recommended production plan (not started; for the user to approve)

Stages write DONE/FAILED markers as the forward queues; the scrambled control runs every day.

| # | stage | what | ETA (PT, if started Thu 2026-10-08 9 am) |
|---|---|---|---|
| P1 | download | CN4 2003-01 → 2009-12 from NRCan (≈ 9.7k station-days still missing), POLARIS PO (6 stations, ≈ 5.8k) and LZB 2003 → 2006-06 from NRCan, C8 (5 stations) and PB from EarthScope; 8 NRCan lanes + repair pass; CN/UW 2010–2013 from EarthScope. Raw ≈ 0.8 TB (plan A 0.48 TB + POLARIS/LZB ≈ 0.3 TB) | NRCan-bound: ≈ 17k station-days at ≈ 7/min ≈ 40 h → **Sat 10-10 ~1 am**; EarthScope parts ≈ 4 h in parallel |
| P2 | 20 Hz cache | forward chain (`bw_gate_cache.py`), ≈ 1.3 TB | +3 h after each block of P1 (runs behind it) |
| P3 | templates | era templates as G3 for the POLARIS/LZB stations; SNB-40T A/B (2012 vs re-timed era template) on the 2005-09 ETS, decided by the G3 rule; template epoch log | ≈ 2 h, Sat 10-10 morning |
| P4 | quiet rule + thresholds | pre-2010 quiet-day flag (decide with Merlin: ≥ 30 d from any Bostock window, or array energy); T_core/T_ext per network era by the forward false-alarm rule | ≈ 2 h + Merlin |
| P5 | matched filter | 2003-01 → 2013-12, every day, data + scrambled control: frozen network where it existed, **era-augmented detector 2003-02 → 2006-06**, both variants on the overlap months so each era boundary has a measured step; ≈ 4,018 days ≈ 6–8 GPU-h | Sat 10-10 evening |
| P6 | confidence + one look | B0–B4 models as SCOPE §5; era-augmented model fitted on 2004 build blocks only; one look each at 2003, 2005-09, 2007, 2008, 2009 (2004-07 test blocks already seen in G3) | Sun 10-11 |
| P7 | controls + catalogue | scrambled control, withheld S-energy, station-set invariance, step checks at every POLARIS station loss (2004-08 KHVB, 2004-09 SHVB, 2005-01 SHDB, 2005-05 TWGB, 2006-05/06 the rest) | Sun 10-11 – Mon 10-12 |

Storage: raw ≈ 0.8 TB + caches ≈ 1.3 TB, peak ≈ 2.1 TB of 58 TB free. Main risks: NRCan throughput (the only
long pole; more lanes untested), the CN4-era confidence model's control rate (4.9–6.4 % at P > 0.5), and the
truncated 2004/2005 catalogue times (harmless for ± 2 s matching, must be re-timed for any stacking).
