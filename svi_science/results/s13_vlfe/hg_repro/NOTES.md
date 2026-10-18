# hg_repro NOTES (PT times)

## 2026-10-15
- ~2:00 am PT: read the paper in full (12 pp.), RESULTS.md scan pilot sections, MERLIN_SCAN_NOGO_REVIEW.md, hg2018_table1.csv.
  Station sets per template are NOT published; Fig. 1 shows 9 stations (CN.NLLB, CN.VGZ, UW.LRIV, UW.DOSE, UW.GNW, UW.SP2,
  TA.A04D, TA.B05D, TA.D03D). Fig. 2 shows T1 uses SP2, T2 B05D, T4 A04D.
- Paper inconsistencies noted: inter-ETS = 6 Sep 2011–30 Aug 2012 (§2.1) vs 1 Oct 2011–30 Aug 2012 (§3.3); 18-month run
  1 Jul 2011–1 Dec 2012 (§2.1) vs 1 Jun 2011–31 Dec 2012 (§3.3); noise windows "during the 2011 ETS" (§3.3) run over
  1 Jun 2011–1 Dec 2013 (§2.3). Text 1,394 detections in the 2011 ETS run vs ~776 implied by Fig. 6 for the same window.
- Digitised Fig. 6 and Fig. 7 (300-dpi renders in data/). Fig. 7 digitiser first version picked up the legend line
  (~0.44–0.8 jumps) -> fixed with curve tracking (±12 px), max step now ≤ 0.03.
- On-disk raw (TRACE_ROOTS) for 1 Jun 2011–31 Dec 2012: NLLB 576/580 d, VGZ 580, LRIV 580, DOSE 573, GNW 579, SP2 410
  (missing Jun–Jul, Sep–Dec 2011), A04D 395, D03D 382 (missing Jun–Jul, Sep–Dec 2011), B05D 78. 2013 (to 1 Dec): all ≥ 325 d
  except B05D 50. EarthScope: TA has LH (1 sps); UW.SP2 BH only; TA.A04D not in ../data/resp/TA_bb.xml.
- 2:15 am PT: prep test days TA.B05D 2012-061 and UW.SP2 2011-278 were made before their neighbour days were fetched (edge NaN) -> moved (not deleted) to data/disp1sps_stale/; they are re-made in the main prep.
- 2:22 am PT: launched scripts/hg_queue.sh (prep TA -> wait fetch -> prep UW -> cont -> mf real -> noise -> analyze); markers PREP_DONE CONT_DONE MF_REAL_DONE MF_NOISE_DONE ANALYZE_DONE. Far refs built (REFS_DONE; 580 days each, 2011-06..2012-12).
- ~2:35 am PT: launched scripts/hg_queue2.sh (waits for main queue; GF caches -> controls -> stack-and-invert); markers GF_DONE CONTROLS_DONE STACKINV_DONE. Fetch finished: 3.20 GB transferred (TA LH 1,148 days, SP2 167 + DOSE 7 BH days decimated; 39 days no data at EarthScope), 1.1 GB stored in data/raw.
- ~2:53 am PT BUG found and fixed: den (90-s window norm) used one global cumsum over 914 days; the 0.82-0.93 m glitches
  at TA.A04D and TA.D03D on 2012-06-07 destroyed float64 precision afterwards, so ~70 % of their later windows were set
  dead (A04D Z live 0.55, D03D Z 0.54). Fix: cumsums restarted every hour with local demean (hg_cont.py --den-only);
  live now 0.97-0.997. Queues stopped (PIDs killed); stale real_T*/noise_*/noise_times moved (not deleted) to
  data/mf_stale_den_v1/ with the old den files; real runs + noise re-run. Controls' window-RMS cumsum now caps |x| at
  50x the 99.9th percentile. Also: V-snr (template SNR >= 3) leaves 0 channels for T1, T2, T5, 1 for T3, 3 for T4 — the
  variant is empty for most templates (recorded; it stays a column).
- ~3:05 am PT: real runs redone with the fixed den (MF_REAL_DONE 03:04); T1 24 channels (D03D dead at its origin), others 27.
  Noise stage running (4 × 25 noise templates), then analyze; queue 2 waits on it (GF -> controls -> stack-and-invert for
  thrA and the best K1 threshold). Agent stopping here per the >20-min rule; resume = read RESULTS.md status + tables.
- ~3:10–4:15 am PT: noise stage, analyze, GF, controls, stack-and-invert finished (all markers). τ = 10 s (template mean VR
  0.28 vs 0.24/0.26/0.22 for 2/5/20 s). best_threshold.json: thr_rate079 (K1 4.4/day vs 8.5).
- ~4:20 am PT (coordinator request): added scripts/hg_lfe.py (LFE coincidence ±60 s, SVI P≥0.5 + WA-NW selected_v2, vs
  random times and scrambled-control detections; stack self-similarity) and scripts/hg_fig_selbias.py
  (figures/fig2_selection_bias.png). Post-hoc analyses, not pre-registered; reported as such.
- ~4:40 am PT RESULTS.md written. Headlines: literal threshold → 7 detections in 18 months (the paper's literal recipe is
  inconsistent with its own 0.6–0.8/day false rate); rate-calibrated threshold → 4.4/day (paper 8.5), ETS/inter 2.0
  (paper 1.3); T1 per-period pattern reproduces. Survival vs our scrambled controls: ETS 0.70/0.74, inter-ETS 0.45; vs
  H&G-style noise templates: 0.31/0.39 and 0.01. Noise-template stacks resemble their selecting window at CC 0.60–0.79
  (selection bias) → their Fig.-2/stack validation is not discriminating. Whole-space CMT cannot reproduce their stack
  inversions (template VR 0.13–0.47, T1/T3/T5 49–104 km off). ETS detections carry LFEs 2.5–3.8× random (T2/T4 4.6–5.5×);
  inter-ETS detections at the random rate (0.4 % vs 0.5 %).
- ~10:50–11:00 am PT (R2, user: "can't we just reverse engineer to reproduce what they got"): tried to fetch the
  supporting information (Table S1, Text S1). Wiley/AGU pages and every downloadSupplement URL give HTTP 403 (Cloudflare
  challenge); the open copies (ETH, sonar.ch) have no supplement; figshare nothing. Not accessible -> continue with the
  digitised targets; no detection-by-detection match possible. Log notes/papers/HG2019_SI/SOURCE_LOG.md. The AGU FM15
  abstract of the same work says they then used 6×MAD and removed ANSS earthquake windows -> added variant V8 (masks).
- ~11:00 am PT: per-template split of the Fig-6 total cannot be recovered from Figs 6+7 (non-negative fit 574/725/1436/
  213/1411, bootstrap ranges ×2 or more; equal split fits almost as well) -> primary per-template targets = equal totals
  × Fig-7 fractions (tables/r2_targets.csv). (Caught the pandas datetime64[us] trap in the first fit attempt: all bins
  came out equal; cast to ns.)
- ~11:00 am PT: PLAN.md Amendment R2 written (targets, calibration = 1 Jul–6 Sep 2011, held out = inter-ETS, 2012 ETS,
  post, noise false rate; pass = pooled inter-ETS and 2012-ETS counts within ±25 %; variants V0–V8; choice by pre-ETS count,
  tie-break daily r). Then launched scripts/r2_queue1.sh (velocity arrays + variants V0–V6 for T1–T5, 5 workers;
  marker R2_VAR_DONE). New scripts: hg_r2.py, hg_r2_calib.py, hg_r2_heldout.py, hg_r2_stage2.py, r2_queue2.sh.
- ~11:10 am PT DISCLOSURE (test run): to test hg_r2_calib/heldout/stage2 before chaining them, I ran them in the
  scratchpad on fake variant files built from the R1 V0 peaks. That printed held-out ratios for V0 and V8 on real data
  (C1 thresholds: V0 inter-ETS 0.77, 2012 ETS 1.32; V8 0.98 / 1.52 of the paper). The variant choice is made by the code
  with the rule fixed in Amendment R2 (pre-ETS count, tie-break daily r); nothing in the rule, targets or variant list was
  changed after this. The test also caught a bug: the daily-r reference was flat (pandas date_range is datetime64[us] ->
  interp out of range); fixed before any real calibration run. Stage-2 test on R1 data reproduced the R1 counts (inter-ETS
  1,353 real, station-day-shift 643, polarity 548), so the stage-2 code path matches R1.
- ~11:12 am PT: launched scripts/r2_queue2.sh (waits for R2_VAR_DONE -> calib/choice -> [V7] -> held-out -> MF at chosen
  setting if not V0/V8 -> stage 2). Markers R2_CALIB_DONE R2_HELDOUT_DONE R2_MF_DONE R2_DONE; logs logs/r2_*.log.
- ~11:12 am PT: the user downloaded the supplement (Table S1: 4,355 detections, 1 Jul 2011–29 Nov 2012; SI text). Stopped
  the R2 queues by PID at ~11:14 am PT (stage 1 was mid-way, no variant file written, no R2 number exists). Facts from
  Table S1 itself: per-template totals 167/2,219/1,375/127/467 (T2+T3 = 82 %); within-template spacing down to 22–24 s
  (T4 ≥ 490 s); their detections of the template events sit 9 s (T1, T4, T5) / 12 s (T2) before the Table-1 origin
  (T3 none) -> time stamp ≈ window start ≈ origin − 9 s. 472 times have fractional seconds.
- ~11:15 am PT: PLAN Amendment R2a written (Table-S1 targets replace the equal split; choice by recall of Table-S1
  calibration-period detections within ±30 s; new variants V5b window origin−9 s, V6b declustering 20 s; detection match
  ±5/±30 s; tests on Table S1's own detections). Restarted scripts/r2_queue1.sh (variants incl. V5b, V6b) at 11:15 am PT.
  Note: in killing the old queue I also matched my own shell with pgrep -f (harmless); the queue scripts now use exact
  patterns. New scripts hg_r2_match.py, hg_r2_s1.py; tables/hg2019_tableS1.csv, tables/r2_targets_S1.csv.
- ~11:18 am PT: blind code test of the whole R2a chain on fake variant files (stdout discarded, only errors checked), so
  that no held-out number is seen before the automatic choice.
- 2026-10-15 11:26 PT: R2a calibration choice (1 Jul-6 Sep 2011 only, Table-S1 targets, before any held-out count): V4 (vel,0,sum,all,90), recall of Table-S1 calibration detections within 30 s 0.101, pre-ETS |ln| 0.007, daily r 0.34; groups beating V0: {'other': 'V4'}; V7 needed: False. Table tables/r2_calibration.csv.
- ~11:26–11:40 am PT: R2a calibration chose V4 (velocity): recall of Table-S1 calibration detections within ±30 s only
  0.101 (all variants 0.009–0.101); held-out FAIL (inter-ETS 2,094 vs 2,923; 2012 ETS 530 vs 419); detection match all
  periods 7 % (±30 s). Coordinator/user: treat the low overlap as our bug until shown otherwise -> PLAN Amendment R2b
  (~11:40 am PT, before any R2b number): self-detection check, stacks at Table-S1 times vs controls incl. a
  selection-matched noise-template control, overlap vs rank / tolerance / lag histogram, self-consistency of our own
  variants, an 11-point processing grid on the calibration period. scripts/hg_r2b.py, scripts/r2b_queue.sh launched
  ~11:42 am PT (waits for R2_MF_DONE); marker R2B_DONE.
- ~11:48 am PT: PLAN Amendment R2c (moment budget of Table-S1 detections) written before computing; scripts/hg_r2c.py.
  First run (11:50): the pre-registered plain SUM is dominated by a few huge windows (T2 inter-ETS sum 6.9e18 N m with
  median scale only 0.25; bootstrap 3.5e18–1.0e19). POST-HOC (flagged): added a robust total = N × 5–95 % trimmed mean
  of the scales, for data and both noise references; rerun ~11:51 (logs/r2c_v2.log). Both are reported.
- ~12:08 pm PT: R2b finished (self, overlap, stacks, grid: best G9 20-s declustering 0.111 vs V4 0.101; G0/G10 code checks
  reproduce V4/V0 exactly; T4/T5 recovered 0–1 % at every grid point). R2c rerun with robust totals done. RESULTS.md
  section "Reverse-engineered reproduction" written. Verdict: both explanations, in different parts of the catalogue —
  a shared strong core (their top 10 % matched 52–77 %), a noise-dominated near-threshold tail (our own variants agree
  only 24–67 %; noise templates fire as often as real templates at Table S1's count level), and a T4/T5
  station/processing mismatch (their T4/T5 between-ETS windows stack to CC 0.03–0.06 with our templates).
- ~12:28 pm PT (user: "don't you have at least a few real-looking VLFEs in between?"): PLAN Amendment R3 written before
  any R3 number (candidate set, per-event tests a–g with fixed pass rules, selection-matched held-out-station aggregate,
  classification). New scripts/hg_r3.py (rows / events YEAR / aux / agg), scripts/hg_r3_summary.py, queues
  r3_queue.sh, r3_events.sh, r3_ets_queue.sh; outputs data/r3/, tables/r3_*, figures/r3_*. Merlin not consulted (flagged).
- 12:31 pm PT: first launch used 6 event workers + aux + agg (8 processes > 6 rule) → killed within seconds, restarted
  with 4 event workers. (My pkill -f pattern also killed my own shell; harmless, now kill by PID.) Events stage first
  failed on pickling pandas namedtuples → SimpleNamespace rows; re-run 12:33 pm PT. Markers R3_EVENTS_DONE, R3_AUX_DONE,
  R3_AGG_DONE (all ~12:35 pm PT).
- Candidate set: 2,886 Table-S1 detections in 8 Sep 2011 – 28 Aug 2012; 148 matched ±30 s after the 9/12 s offset;
  130 distinct (T2 100, T3 28, T1 1, T5 1).
- Result: 0/130 real-looking (top 30: 0); scrambled control 0/650; selection-matched noise 1/1,300. 101/130 inside
  catalogued-earthquake windows (78 % vs 52 % random; surface waves 72 % vs 47 %, body waves 22 % vs 8 %, M ≥ 3 ≤ 2,500 km
  12 % vs 4 %). Held-out-station G −0.018 (top 30) / −0.003 (all), 27–28th percentile of selection-matched noise.
- ~12:36 pm PT POST-HOC positive control (not in R3; same code, R3_TAG=ets): Table-S1 detections inside the 2011 and 2012
  ETS windows (self-detections removed): 115 distinct candidates → 8 real-looking (top 30: 4), all T2, Mw 3.1–3.5, ≤ 6 km
  from a family, all with LFEs; includes the scan-pilot-3 event 2012-09-19 19:18 PDT. Scrambled 2/575, noise 1/1,150.
  Between vs ETS: Fisher p = 0.002. So the battery has power; the between-ETS zero is not a too-strict-rules artefact.
- Post-hoc power note for G: template events at the held-out stations have window/pre-event RMS 0.4–1.7; T2 vs T4
  held-out CC 0.36; ETS top-30 G 0.045 (top of random, 63rd pct of selection-matched) → G only sensitive to
  template-sized repeats.
- Near-misses (fail one test, post-hoc look): cands 7, 49, 116, 128 — three far-loud long-period wavefields, one
  single-station PFB pulse; waveforms figures/r3_nearmiss_wave_cand*.png.
- RESULTS.md section "Individual between-ETS candidates" written ~12:45 pm PT.
