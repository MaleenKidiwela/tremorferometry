Written at launch (2026-10-03), before any production result existed. Each item: what differs from WA_PROD_PLAN.md + the
user's decisions of 2026-10-03, and why.

1. **D1, A2011 leakage: rebuild instead of flag.** The step-4b templates (`tpl_step4b_a`) were built from A2011 array
   events. They are replaced by `v10b` = the same frozen recipe on the V2010 BUILD-half array detections at the same nodes
   (7,636 templates, 41 nodes, 200 families; step4b_a had 9,246). A2011 is then out-of-sample for every station (plan D1
   recommendation). The V2010 build half stays in-sample.
2. **Fixed rate networks by a rule written before any run** (`setup/core_sets.json`, `core_set_table.csv`): permanent
   station, raw on >= 80 % of the era's days, a W4a template source (own node with >= 30 families, or a W2-approved
   carry-over) on >= 95 % of its raw days, LRIV out after 2019-05-30 (plan D10 / SVI policy), never XU.W060. Headline E3
   2015-2026 = 19 stations (B001, B003-B007, B009-B011, B013, B014, B017, B018, B019, B926, B941, B943, GMW, GNW), not the
   21 of the plan: PGC, VGZ, DOSE and HDW have no usable W4a template for part of 2015-2026 (their control-level epochs are
   re-templated, user decision), so they cannot be in a FIXED set; B019, B941 and GMW are now on disk (tiers B/C). E1
   2007-2011 = 12 stations, E2 2012-2014 = 37.
3. **Re-templated epochs and other new nodes enter only the rich pass** (W4b after the confidence model; T2B rule P_own_WAVE
   >= 0.3 from the network without the node, families <= 50 km, the node's first two WA ETS else its first 24 months, with
   the DUC_FIN per-channel scrambled-control check: N >= 20, >= 10 day-shifted controls, coherence above all). Node changes
   of non-core stations have no W2 swap check, so they are re-templated rather than carried over. The sealed months are
   never used to build a template.
4. **B014 2023 check** is the W2 monthly stack-test record before / after its 2023-08 .. 2024-08 outage (no metadata change
   exists): ok share 0.30 before, 0.50 after -> no change detected, old template kept (rule in `setup/b014_check.json`).
5. **Health exclusions**: only W2 dead months and dead days (969 station-days). W2 'flat' windows (>= 1 % repeated samples)
   are NOT exclusions: they flag years of normal PB borehole data. Flat / dead channels are caught per window by the
   frozen energy-floor guard and the drop-dead rule (> 5 % of the day guarded -> station absent that day). Non-core
   stations have no W2 health table; the same per-day guard applies to them.
6. **Scrambled-control shifts** are deterministic per (date, site) (the SVI_REF fix) instead of wa_mf's sequential draws,
   so a site's shift does not depend on which other sites are present. Data rows reproduce wa_mf.py exactly (smoke: 10,698 /
   10,698 detections of the go/no-go core-21 run on 2010-08-20, max |d net_mad| = 0); control rows differ by construction.
7. **Truth runs**: networks E3 + all with 2 scrambled controls + the time-reversed control (`truth`), E1 + E2 with 1 scrambled
   + time-reversed (`truth_e12`) -- GPU memory. XU / XG station-days before 2012 are never used (they made the truth).
8. **W8 PNET gate** uses the frozen fine-tune `ft1.pt`, PNET block on the `truth` run (E3 + all). Rule as the user fixed it
   (dAUC own >= 0.01 with day-bootstrap CI > 0, slope 0.8-1.2); the plan's extra condition on the never-trained stations
   (B017, B018, GNW, HDW) is reported as a column, not applied.
9. **Fit rows** are a uniform random subsample (<= 600,000 data rows of the fit days; base rate kept) to bound compute.
10. **Refit risk 'older model'** = FULL features fitted on the V2008 + V2010 test halves only (the truth the earlier WA
    confidence used); the earlier model itself was never saved as a file.
11. **W9 member probability** is a simplified column: softmax(beta x net_mad) over the <= 3-km neighbours' candidates within
    +-2 s (rows with P_nbr3_WAVE >= 0.3), beta by maximum likelihood on the truth fit rows -- not the step-2 station-level
    attribution, which needs a separate group matched filter.
12. **W11**: k factors for the headline labels (own per family; nbr3 per Olympic A position group); the nbr3 count of
    NW / SW / Arm B families is a column without k. The daily-stack activity column (stack_activity.py, a separate GPU
    engine) and the SVI withheld-station inter-ETS validation were NOT run: inter-ETS rates are a control-subtracted bracket.
13. **W12 windows**: S times at non-template stations from a straight-ray harmonic-mean S speed through the USGS CVM v1.7
    1-D profile (royer_NW_WA / royer_SW_WA) + the family's median template residual. 'Every recording station' is limited
    to stations <= 150 km from a Washington family (the California stations of lfe_south/traces are otherwise in the roots);
    stations with no known coordinates are listed in `cut_roots/skipped_stations.csv`. Pre-registered size limit: the
    pilot month extrapolated to all raw must be < 3.0 TB, else the full cut stops for a user decision.
14. **Span** stays 2007-01-01 .. 2026-09-15 (the W1 / W2 span), although tier A now holds data to 2026-09-30.
15. **Rich pass network** = every station with a W4a or W4b template; it runs only on days with a W4b node live.

16. (main session, 2026-10-04) The first fit crashed: fam_code (264 families) exceeds the 255-category limit of
    HistGradientBoosting. fam_code is now a plain numeric feature (only in the FULL / NO_DENSITY sensitivity models;
    the headline WAVE model never uses it); reg_code stays categorical. Queue relaunched; finished stages resume.
16b/main 10-04 09:10Z: production MF stalled 05:00Z→09:10Z waiting for cache/done_p2_2013-09; that month has no phase-2 (traces_dl2) jobs, so the cache stage never wrote its marker. Markers written for all months with 0 phase-2 jobs. Code fix for future runs: write done_p2 for empty months.
17. (main session, 2026-10-06) w4b_targets crashed 10-05 19:43Z: cache manifest_p1_2013-09.csv has only 'empty' rows and was written without a 'node' column. node_days() now skips such files (backup .bak_pre_nodecol); relaunched, finished stages resume. ~11 h idle before it was noticed.
18. (main session, 2026-10-06) User decisions after the cut pilot and the quiet-day waveform test: window floor P_max >= 0.2 (conf/floor_user.json overrides the automatic 0.1; pilot at 0.1 kept in cut_pilot_floor01/); NW quiet-day cut-off 0.85; days without a tremor catalogue labelled 'unknown'. Queue relaunched for pilot -> cut -> verify -> stack check -> delete list -> report.
19. (main session, 2026-10-06) Floor 0.2 pilot extrapolated 3.63 TB (> 3 TB). User: floor 0.3 (lowest headline cut-off is 0.5). Window selection now uses P_sel = max(P_max_WAVE, P_own_WAVE_PNET) because the headline score is the PNET model (backup .bak_pre_pnetsel). Pilot month 2015-10 at 0.3: 170,050 candidates. Relaunched.
20. (post-processing v2, 2026-10-06; `scripts/wa_prod_post_v2.py`, outputs `post_v2/`, v1 `post/ cutoff_check/ rates/` untouched) User decisions: day classes tremor (PNSN tremor ≤ 30 km; before 2009-08-06 the WECC windows) / quiet (catalogue covers the day, 2009-08-06..2026-05-27, nothing ≤ 30 km) / unknown (no catalogue). Headline single repeats and rates: P_own_WAVE_PNET ≥ 0.75 for NW and Arm B on every class, ≥ 0.5 for SW on tremor days only; Olympic A group rates on P_nbr3_WAVE ≥ 0.75 (applied on every class; quiet / unknown are mostly false by the scrambled control). WAVE model and the earlier class rules (tremor 0.75, SW 0.5, NW quiet 0.85) are kept as secondary columns. PNET k factors cannot be measured for E1 / E2 (the truth_e12 run has no PNET block), so E1 / E2 headline rates are net counts without k. The 'quiet' class now includes days with tremor 30–60 km away (v1 used 60 km); the 60-km subset is the diagnostic class quiet_beyond60.
21. (main session, 2026-10-06 18:41Z) Pod restarted 17:41Z during the full cut (memory ~176 GB+ with 10 workers on 2007-01). Relaunched with WA_PROD_CUT_WORKERS=6; no other change.
22. (main session, 2026-10-06 20:11Z) Second pod death 19:37Z (anon 181 GB, 6 workers, 2007-01 at day ~20). Cause: cut_windows.cut used
    pool.imap, which reads ahead without limit; the single gzip h5 writer is slower than the workers, so the raw-day backlog grew
    until OOM (the 1-month pilot was bounded by its 30 jobs, the 7,198-day cut is not; worker count does not bound it). Fix:
    ordered bounded_imap with at most 2 x workers days in flight (output identical, same order). Relaunched, 6 workers.
23. (main session, 2026-10-06 22:14Z) Cut stopped after 14 months (2007-01..2008-02 done) to remove the next OOM risks: (a) the month index
    was held as Python tuples until month end (peak 83 GB at 25 M windows in 2007-01; 2009-05 / 2010-09 have ~4x the candidates, ~100 M
    windows) -> now built per day and streamed to index_<m>.parquet (one row group per day), same columns; (b) fix_offset_source read
    whole month files -> now row group by row group; (c) stackcheck loaded the whole store -> now a filtered read (same selection);
    (d) CW.verify loaded the whole index + all segments -> now reads only the n sampled rows (load_index_sample) and their months'
    segments (different random rows than before). mem_watchdog: 5-s polls, kill at 175 GB (was 60 s / 185 GB).
    Backups: cut_windows.py.bak_pre_streamidx, wa_prod_cut.py.bak_pre_streamidx, scripts/mem_watchdog.sh.bak_60s.
    Memory test before relaunch: months 2008-02 (compare with production output) + 2010-09 (busiest) in results/cut/memtest_1006/.
24. (user, 2026-10-06 ~22:30Z) Keep window floor P_sel >= 0.3, native sample rate, all components; let the full cut finish at full size (~1-2 TB projected).
25. (main session, 2026-10-06 23:15Z) memtest_1006: streamed index reproduces production 2008-02 exactly (3,821,369 rows, same columns
    / dtypes / values; h5 datasets identical). But 2010-09 with 6 workers and 12 days in flight peaked at 153 GB (busy days ~6 GB of
    waveforms each). In-flight default lowered to workers + 1 (env CUT_INFLIGHT); production cut to run with 4 workers (writer-bound,
    no throughput loss). Re-test: first 12 days of 2010-09, 4 workers.
26. (main session, 2026-10-06 23:57Z) Re-test passed: first 12 days of 2010-09 (76.1 M windows, 13.5 GB), 4 workers, in-flight 5 -> peak
    anon 104 GB, ~90 GB sustained. Production queue relaunched with WA_PROD_CUT_WORKERS=4; resumes at 2008-03 (14 months done).
    (The queue cleared the WA_PROD_FAILED marker of the deliberate stop itself on start.)
27. (user, 2026-10-07 ~02:00Z) User stopped the window cut ('window cutting is a waste of time'). All WA_PROD processes, the cut, the day-file test and the report_v2 watcher stopped. On disk: 16 cut months (2007-01..2008-04) in results/cut/wa_ext + partial 2008-05 .tmp; test dirs memtest_1006, memtest2_dayfiles. Uncommitted day-file mode in cut_windows.py (CUT_DAYFILES, backup .bak_pre_dayfiles) untested. Nothing deleted.
28. (user, 2026-10-07) Deleted all WA cut windows (user OK): cut/wa_ext (110 GB), cut/memtest_1006 (20 GB), cut/memtest2_dayfiles, wa_prod/cut_pilot{,_floor01,_floor02} (43 GB). Small summaries (PILOT.json, VERIFY, cut_summary, done_*.json, offsets) kept in wa_prod/cut_abandoned_summaries_2026-10-07/. Raw continuous data untouched; next steps read the raw waveforms directly.
29. (main session, 2026-10-06 ~8 pm PT) k-factor fix (Merlin finding, verified): k_factors_v2 counted truth events from every build/test half while the detection side (truth_rows row_role == "fit") holds only the test half on V2008 / V2010 days (build half = insample, edges = guard). Verified: 2008-05 fit rows half 1 only; truth events halves 0/1/−1 = 25,530/24,444/895 → k ~2× too high on those days. Fix: truth events restricted to the (date, half) blocks present in the fit rows (keeps 24,444 of 50,869 in 2008-05, 30,662 of 65,041 in 2010-08, all 32,820 in 2011-08). v2 rates and tables re-run; old outputs in post_v2/pre_kfix_2026-10-06/, script backup wa_prod_post_v2.py.bak_pre_kfix. The superseded v1 script (wa_prod_post.py l.232–234) has the same pattern; v1 rates in rates/ are left as they are (superseded, now also known biased).
30. (main session, 2026-10-06 ~9:10 pm PT) Zero-day fix (found by the per-episode check, verified): rates_daily_v2 held only unit-days with ≥ 1 candidate (any series), so scanned days with no candidate were missing. NW E3 family-tremor-days: 10,691 of 57,831 (18.5 %) missing, and every one of them had 16–18 live template stations in that day's core run → real zeros. Fix: scanned_unit_days() (a unit is scanned on a day if ≥ 1 of its template nodes in setup/tpl_library.csv is live in mf/prod/run_<date>.json members.core; Olympic A group = any member) adds zero rows before day classes and k are applied. k factors were already computed over all scanned days (unchanged). Rates and tables re-run; previous outputs in post_v2/pre_zerodays_2026-10-06/, script backup wa_prod_post_v2.py.bak_pre_zerodays.
31. (main session, 2026-10-07 ~10:45 am PT) "Counted once" columns, user-approved; flags only, nothing deleted, raw counts stay the headline. Neighbouring families detect the same event, and the WA rates had no cross-family de-duplication (Olympic A groups merge only their own members at ±2 s). New script `scripts/wa_post_dedup.py` (`wa_prod_post_v2.py` not edited; the fixed-10 queue uses it). It writes new files only in `post_v2/dedup/` (README_dedup.md). Rule = the SVI keep_xfam8 rule: per day × series (data, time-reversed and scrambled control alike), a candidate is flagged when a station-sharing family's candidate within ±2 s has a higher raw net_cc. Variants: A any distance (upper bound); B ≤ 10 km; C = B and D = A with the two detections' best-station arrivals agreeing within 0.5 s. At the headline cut, A flags 31 % of NW repeats in every era. 82–86 % of those pairs are arrival-consistent, against about 5 % for time-shifted chance pairs, and most same-event pairs are 10–40 km apart. So B flags only 3–7 % and D 27–28 %. NW E3 tremor net per family-day: raw 0.282, B 0.261, D 0.207, A 0.195. k refitted on de-duplicated truth + detections (truth events are double counted too: A 55,231 of 198,804), giving rates 3.70 / 3.55 / 2.98 / 2.70. The existing k must not be applied to de-duplicated counts. The raw columns reproduce rates_daily_v2 exactly.
