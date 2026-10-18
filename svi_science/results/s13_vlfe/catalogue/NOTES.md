# catalogue NOTES — 2003–2026 strong-VLFE catalogue, S Vancouver Island (stage 1)

## 2026-10-15 (Thu) ~1:00–1:20 pm PT — set-up, PLAN_CAT, e′, E2 launch
- Read: MERLIN_CAT_DESIGN.md, hg_repro RESULTS (all), PLAN R2–R3, r3 tables, hg_r3.py, scan_columns2.py, PLAN_SCAN S5.
- 14 events identified from `../hg_repro/tables/r3_ets_events.csv` with ratio_far (n_far ≥ 2, ≤ 0.551) as (a): 8 vetted
  (a ∧ b ∧ c ∧ e ∧ f) = the R3 eight; 6 near-misses fail only (e), each by a head2 window of a M 4.6–5.0 teleseism
  8,000–11,000 km away (Merlin F2 confirmed).
- Disclosure (before PLAN_CAT): to size the battery compute I printed raw V4 detection COUNTS per MF level
  (5.5/6/7/8/10 MAD) for T1–T3, the 16 scrambled-control catalogues and 25 noise templates over 2011-06-01–2013-01-01
  (e.g. T2 at 5.5 MAD: real 2,373, control mean 601). No battery, no ETS split, no rate per year.
- PLAN_CAT.md written ~1:10 pm PT, before any battery/e′/E3 number. Amendment ~1:15 pm PT (before any E3 number):
  complete linkage instead of single linkage for location clusters (10-km scan grid would chain).
- Step-5 metadata: the EarthScope fdsnws/availability service (and its extent endpoint) answers "This service has been
  retired" (also at service.earthscope.org). Switched to the MUSTANG percent_availability metric (daily QC metadata,
  no waveforms). NCEDC needs full ISO times (HTTP 400 fixed). NRCan holds no CN LHZ channels in 250–800 km; CN
  broadband BHZ/HHZ listed as a separate class (Merlin's fallback class).
- e′ validation (`scripts/cat_eprime.py`, ~1:09 pm PT): ETS 8 → 14 (the 6 near-misses, nothing else), between ETS 0/130,
  scrambled 2/650, selection-matched noise 1/1,300 (unchanged) → **e′ accepted** (rule: ≤ 14, ≤ 2/650, 0).
  ETS-window scrambled 2 → 5/575 (column).
- E2 (`scripts/cat_e2.py`, queue `scripts/e2_queue.sh`, marker `E2_DONE` 1:13 pm PT): 115,209 battery rows (81,469 unique
  times); expensive columns only for the 4,302 times passing (c) and (e) or the e′ time part (PLAN_CAT §3 shortcut).
  Runtime ~2 min with 8 workers.

## ~1:13–1:20 pm PT — E2 results (tables/e2_separation.csv)
- Accepted tiers exist. Lowest MF level that passes: 5.5 MAD.
  - e′ + two-template agreement: ETS 18, between ETS 1, scrambled controls 0 (mean 0, max 0), noise templates 0.
  - (e) + two-template: ETS 12, between 0, controls 0.
  - e′ single template, 6 MAD + S5: ETS 21, between 2, controls 0.13/yr (max 1), noise templates 0.62 per template-yr.
- Fails: e′ single template at 5.5 MAD without S5 (ETS 33, between 12, controls 1.48/yr > 1).
- The one between-ETS e′ two-template event: 2012-02-13 05:44:40 PST, T3 (agreement 6.05 MAD), z_best 8.3, ratio_far
  0.40 (5 refs), S5 pass, 2 SVI LFEs within ±10 min, near SVI family 79 (2.4 km); fails (e) by a M 4.8 head2 window.
- Per-test pass fractions are per row and include rows not given the expensive columns (they count as failing a/b/f).

## ~1:15–1:25 pm PT — E3/E4 launch, far-reference metadata
- Station data: `scripts/cat_cont.py` built 0.02–0.05 Hz velocity for CN.PGC (555 of 914 frame days), CN.SNB (510),
  CN.PFB (334), CN.YOUB (135) in the hg_repro frame (`data/cont/`). Stations: 13 distinct (Merlin's "16" = 9 + 7 with
  VGZ, NLLB, LRIV counted twice).
- Bank (`scripts/cat_e3.py bank`, 1:15:54 pm PT): 14 events → 6 complete-linkage clusters (C0: 1, 90, 3, 6, 107;
  C1: 11, 62, 36, 96; C2: 32, 77; C3: 0; C4: 15; C5: 42 near-miss); 6 leave-one-out stacks (vetted events 0 and 15 are
  single-event clusters, so other templates must find them); H&G T2 and T4 re-cut on 39 channels; 25 noise templates.
- MF queue `scripts/e3_queue.sh` launched 1:16 pm PT (156 jobs: 32 with 16 scrambled controls, 124 real-only; 8 workers,
  ~43 GB anon). Code check on the first outputs: C1 finds its own members at 15–16 MAD and the other vetted events at
  6–8 MAD; the H&G T2 re-cut finds all 8 vetted events at 5.8–9.5 MAD (before thresholds; not a result).
  **Marker `catalogue/E3_DONE`; resume = re-run `scripts/e3_queue.sh` (finished MF jobs are skipped), then write the E3
  section of RESULTS.md from `tables/e3_summary.csv`, `e3_refind.csv`, `e3_thresholds.csv`.**
- Far-reference metadata done (`FARREF_META_DONE`): see RESULTS §3. Bug fixed before reporting: metadata days per year
  were the longest single epoch, not the union of epochs (US.NEW 2011 read 243 instead of 365).

## ~1:20–1:25 pm PT — far-reference download (user approval relayed by the coordinator ~1:20 pm PT)
- Set (from `tables/farref_per_year.csv`): IU.COR, US.NEW, US.HAWA (to 2024), US.MSO, US.WVOR — LHZ only, 2003-01-01 –
  2026-10-15; ≥ 4 per year in every year (≥ 3 rule met); azimuth spread 89–177° (no LHZ exists north or west).
- No year has < 2 LHZ far stations → no other data class fetched.
- `../scripts/scan_fetch_ref_all.py` (nohup, launched 1:24 pm PT): EarthScope dataselect, monthly requests, 2 in parallel,
  retries 10–40 s, default user agent, no personal data; day files in `../data/raw_ref_all/NET.STA/YYYY/DOY.mseed`
  (raw_ref2011 layout, readable by the scan_prep_v3/scan_prep_ref path when used); StationXML `../data/resp_scan/*.all.xml`;
  COR/NEW/HAWA 2011–12 already on disk are skipped; stops if > 10 GB. Log `../logs/ref_all_download.log`; marker
  `../REF_ALL_DONE`; resume = re-run the script (month markers). First months: COR 2003 ~6 MB/month (≈ 70 MB/yr).

## ~1:16–1:35 pm PT — E3/E4 results (marker E3_DONE 1:33 pm PT; ~17 min)
- Thresholds at 0.5/template-day: all13 6.98/6.60 MAD (2011/2012), hg9 6.86/6.71, cn4 6.08/5.84, cn4gnw 6.29/6.05.
- LOO re-find at 0.5: all13 6/8 (fail ≥ 7), hg9 7/8, cn4 2/8, cn4gnw 1/8. cn4 scrambled 1.16/yr (fail ≤ 1).
  → by the pre-set rule 2003–09 = Mw ≥ 4.0 tier; uniform catalogue from 2010. Scope: test events Mw 3.1–3.5.
- all13 ETS 49 / between 1 / scrambled 0.13/yr.
- E2 waveform figures (`scripts/cat_wavefig.py`, 12 between-ETS events): the 2012-02-13 two-template event shows no pulse
  and far references carry coherent long-period wavetrains → most likely noise passing the battery.
- RESULTS §4–6 written ~1:40 pm PT.

## ~1:40–2:00 pm PT — Amendment P (pilot 2014 + early-era injection), set-up
- Coordinator relayed: proceed with Merlin step 6 (pilot 2014) and a 4-station cross-template injection (2011–12).
  PLAN_CAT Amendment P written ~1:45 pm PT before any 2014 or injection number.
- Far-reference download finished 1:45 pm PT: 5.53 GB (5.2 GiB apparent; `du` shows ~20 GB allocated on this
  filesystem's large blocks), 39,290 day files; no-data months: HAWA 30 (ends 2024-09), WVOR 9, MSO 1. Marker
  `../REF_ALL_DONE` (its time stamp was written in UTC by mistake — "20:45" = 1:45 pm PT; script fixed).
- Shared-code edits (additive; 2011/2012 behaviour unchanged): `../scripts/scan_common.py` gains a SCAN_YEAR=2014 branch
  and year-prefixed tables for years other than 2011/2012 (`tables/scan2014_*`; pre-edit copy
  `data/scan_common_before_2014branch.py`); `scan_refs.py` also reads `resp_scan/NET.STA.all.xml`; `scan_prep_ref.py`
  also reads `data/raw_ref_all`. ComCat regional M ≥ 3 within 2,500 km for 2014 fetched (event metadata).
- Scan year 2014 queue `../scripts/scan2014_queue.sh` launched 1:40 pm PT (prep 8 pilot + 6 H&G/regional + 3 far
  stations, noise, masks, bp, refqc, GPU detect real with no polarity draws); marker `../data/scan_2014/SCAN2014_DONE`.
- Injection machinery `scripts/cat_inj.py`. Code checks before use:
  - local scan statistic vs stored 2011 zmax at 40 random seconds: median |Δz| 6e-7, same node 39/39 (one all-NaN time);
  - local MF vs stored E3 cn4 values at event 0: C0 6.26 / C2 4.93 / HG_T4 5.68 = stored.
- Sources: C0 (5 events, Mw_ref 3.43 from mean M0 of 4 members with Mw; PFB missing in C0 → no PFB signal injected)
  and C1 (4 events, Mw_ref 3.21).
- Smoke test (2 injections, not results): MF z of cross-templates on 4 stations stays ~4–5 at Mw 3.4 and 4.2 (CC saturates
  at the waveform similarity between clusters), while the scan z grows ×18. The early-era scan (no LRIV) placed one C1
  injection on the grid's south edge.
- Master queue `scripts/p_queue.sh` (waits for SCAN2014_DONE; B: `cat_inj.py early 2011/2012` + summary, marker
  `INJ_EARLY_DONE`; A: `cat_p14.py cont bank mf detect battery inj_list inj summary`, marker `P14_DONE`).

## ~2:10–3:00 pm PT — B results, queue crash, A04D fix
- p_queue crashed 2:11 pm PT in `cat_inj.py early_summary` (a string column, mf_best_tpl, entered the groupby mean). The
  1,000 early injections had finished without errors (2011: 490, 2012: 510). Fixed (bool columns only; medians of MF and
  scan z added), summary re-run ~2:50 pm PT → RESULTS §7. Mw_c(early) = none ≤ 4.2.
- Queue rewritten: every stage is a `run` with STEP markers in `data/p14/` and a FAIL line in `logs/p_queue.log` on any
  non-zero exit; a missing scan-2014 marker with no running queue also writes FAIL.
- Pilot cont build showed TA.A04D 2014 with 0 days: prep had logged "noresp" for all 365 days (no response in the scan
  lookup order). Fix: copied `hg_repro/data/resp/TA.A04D.xml` (BH epochs cover 2014) to `../data/resp_scan/TA.A04D.all.xml`
  and re-ran `scan_prep_ref.py TA.A04D` for 2014 (365 ok). The A04D-less pilot outputs were set aside, not deleted:
  `data/p14/bank_v1_noA04D.npz`, `mf_v1_noA04D/`, `cont/*TA.A04D.v1_empty.npy`, `STEP_p14_{cont,bank}_DONE.v1_stale`.
  Queue relaunched ~2:53 pm PT.
- 2014 prep failures (no data / errors), for the record: B05D 323 days, LZB 197, YOUB 120 (ends 2014-09-03), PFB 23,
  SP2 10, VGZ 7, COR 7, NLLB 5, SNB 5, NEW 2.

## ~2:53–3:15 pm PT — pilot 2014 results (P14_DONE 3:03 pm PT)
- NO-GO: truth 0/8, between-ETS 3 vs scrambled mean 0.44 (rule ≤ 1.44), scrambled 0.79/yr (ok), completeness 81.5 % at
  Mw 3.8 (ok). Thresholds from the 2014 noise templates are in `tables/p14_noise_threshold.csv`.
- Diagnosis (post hoc, columns only): the scan statistic at the T6–T13 origins is 1.5–3.6. The raw VGZ/LRIV records show no
  pulse above background, and PGC was very noisy in Dec 2014 (99th pct 700–1,000 nm). The December scan network was
  PGC N/Z, VGZ E/Z, LRIV and PFB Z on some days (SNB month-excluded, YOUB ended, LZB excluded).
- RESULTS §8 written ~3:15 pm PT.

## ~3:15–3:40 pm PT — Merlin review 2 + addendum: Q0, Q2, Q1′ launch
- PLAN_CAT Amendment Q (~3:21 pm PT) and Q1′ (~3:36 pm PT) written before their numbers.
- Q0: hg_r3.lfe_lists clip fixed (2013-01-01 → 2027-01-01; copy `data/hg_r3_before_lfe_fix.py`). US.MSO/US.WVOR added to
  `../scripts/scan_refs.py` REFS and `../tables/scan2_ref_stations.csv` (rows appended; pre-edit copy
  `data/scan2_ref_stations_before_Q0.csv`; 2011/2012 unaffected, no bp files there). 2014: MSO/WVOR prepped, bp, new
  `scan_2014/templates_ref.npz` (the 2011 symlink kept as `templates_ref.npz.symlink_2011`), refqc
  (`scripts/q0_refs.sh`). cat_common FAR5 now includes MSO/WVOR (no effect on 2011–12).
- 2014 pilot re-run with validated bank + 5 refs + no between gate: ETS 0 / between 3 / other 1, controls 0.22/yr,
  completeness 0.80 at Mw 3.8, truth 0/8 → NO-GO (truth only). PGC statement in RESULTS §8 retracted (checked: normal noise).
- Q2 triage: all 3 between events (all 2014-01-23 10:59–11:55 PST) "likely noise" (vr_thrust 0.07–0.11, ratio_GNW 0.81–1.07).
- Sensor check: CAFE BS11/W020/W040/PL11/W060 = CMG40T-1 (1 s) in 2011-06..2012-05 (StationXML), so Merlin's dense sets
  with CAFE cannot carry 2011 templates; SQM/HDW/OTR not broadband in 2011–12. Sets re-defined by rule (Q1′).
- Shared-code additions for later pilots: scan_common SCAN_YEAR 2013/2018 branches; `../tables/sc_gain.csv` + 12 rows
  (g = 1, Z only, for the post-2017 PGC/VGZ/SNB/NLLB HH epochs; pre-edit copy `data/sc_gain_before_Q4.csv`);
  `../scripts/scan_year_queue.sh YEAR`; ComCat regional M ≥ 3 lists for 2013 and 2018. Pilot code `scripts/cat_pilot.py`
  (PILOT=2013|2018; 2008 to be added per the addendum) written, not yet run.
- Q1 queue `scripts/q1_queue.sh` launched 3:38 pm PT (cont for 8 stations, bank, MF 6 sets × 39 templates plain +
  weighted, detect, battery 2011/2012, summary); marker `Q1_DONE`; FAIL lines in `logs/q_main.log`.

## ~5:40–6:02 pm PT — results of Q1 / pilots (Q), Merlin review 3, Amendment R
- Q1_DONE 5:0x pm, PILOT2013_DONE, PILOT2008_DONE 5:38 pm PT (numbers in RESULTS §10–11). Weighted ≥ plain everywhere;
  dense2011 8/8 (scrambled 0.26/yr), hg9 8/8, era2016 6/8, era2007 2/8.
- 2008 pilot non-informative: scan 2008 had only PGC N/Z (VGZ, SNB month-excluded all year by the 2011-referenced noise rule;
  PFB/YOUB no 2008 data on disk).
- Amendment R written ~5:55 pm PT before any 2018 number; the Q verdicts were kept (`tables/pilot2013_*_Q.*`, `pilot2008_verdict_Q.json`).
- 2013 R re-summary on n = 200 (`scripts/cat_r.py`): T2b 48 vs 21.0 (p < 0.001) and T3′ 3/3 reproduce Merlin. A first run
  gave T2b = 0 because family ids were str vs int; fixed before reporting. Between-ETS 4 → Tier A′ 0.
- Queues: `scripts/r_queue.sh` (waits SCAN2018_DONE; 2013 n = 800 → R2013_DONE; pilot 2018 → PILOT2018_DONE),
  `scripts/v2_queue.sh` (after PILOT2018_DONE; scan v2 on 2013 → SCANV2_2013_DONE); scan year 2018 building
  (`../scripts/scan_year_queue.sh 2018`). FAIL lines → `logs/q_main.log`.
- Shared-code edits for v2 (env-gated, defaults unchanged): scan_detect.py (SCAN_GRID, SCAN_GAIN_FALLBACK), scan_columns2.py
  (SCAN_GRID), hg_r3.py ev_init (SCAN_RUN, SCAN_BP).

## ~6:29 pm PT — 2011 carrier fetch (user approval relayed by the coordinator, stated as ~7:20 pm PT)
- `scripts/fetch_2011carriers.py` (nohup): UW.SP2, TA.A04D, TA.D03D, UW.WISH, US.NLWA, every 2011–2012 day without any file
  in the trace roots, hg_repro raw, disp1sps_v3 or hg_repro disp1sps: NLWA 392, WISH 305, D03D 166, A04D 151, SP2 137 days
  (1,151). TA: LH 1 sps as delivered (~0.45 MB/day). UW/US BH: decimated in memory to 1 sps, integer counts, STEIM2.
  Output `../data/raw/NET.STA/YYYY/DOY.mseed` (a trace root); log `logs/fetch_2011carriers.log`; marker `FETCH2011_DONE`;
  2-GB stop. Note for production: SP2/A04D/D03D 2011-06 → 2013 days already exist in `hg_repro/data/raw` + `hg_repro/data/disp1sps`
  (fetched for hg_repro), which is not a scan trace root. The production prep must read them (or the v3 day files must be
  built from them) before 2011–12 production.

## ~7:40–7:50 pm PT — Amendment S: 2018 scan v2, carriers, production queue
- Coordinator relay: scan v2 on 2013 accepted (ETS 18; completeness at Mw 3.8 0.853 / 0.814; controls 0.54/yr; Tier A′ 0);
  FETCH2011_DONE (coordinator fixed a mixed 40/50-Hz day in fetch_2011carriers.py; 0.219 GB stored).
- Amendment S written ~7:45 pm PT before its numbers. Pilot 2018 under R (v1 scan) passed T1′/T2′/T3′/controls with
  completeness 0.27/0.22 at Mw 3.8 (recorded in S).
- `scripts/v2_year.sh YEAR` (generalised v2 queue) launched for 2018 at 7:45 pm PT → `SCANV2_2018_DONE`.
- Production mode added to `cat_pilot.py` (env PROD_YEAR; pilot behaviour and RNG unchanged without it; copies of the
  pre-edit scripts in `data/cat_pilot_before_prod.py`, `data/cat_r_before_prod.py`), `cat_prod.py` (stations rule,
  per-episode table, events, known-answer ledger, summary), `carriers_prep.sh` (links hg_repro v3-recipe day files of
  SP2/A04D/D03D into disp1sps_v3, preps the fetched raw days, rebuilds 2011/2012 reference bp arrays into `bp2_ref_v2`).
- Shared code: `../scripts/scan_common.py` generic branch for 2003–2026 production years; `hg_r3.py` and `cat_scanv2.py`
  read SCAN_REFBP (defaults unchanged). ComCat regional lists for 2010, 2015–2017 fetched (event metadata).
- `scripts/prod_queue.sh` launched ~7:49 pm PT: waits for SCANV2_2018_DONE → carriers_prep (CARRIERS_DONE) → years 2011, 2012,
  2013, 2018, 2014, 2015, 2016, 2017, 2010 (scan years 2010/2015–2017 built first) → `PROD_<Y>_DONE`. FAIL lines in
  `logs/q_main.log`. Resume: re-run `scripts/prod_queue.sh`.

## ~8:45–8:56 pm PT — production fix, restart, relay queue (Amendment T)
- 2011 production failed at "cont" (8:02 pm PT, KeyError 2011: pilot CFG dict indexed before the production override). Fixed
  (`.get`), leap years now 366 days, majors flagged; per-year prep of carriers missing from v3 (2014 WISH/NLWA; 2010 and
  2015–2017 after their scan years). Station sets 2011–2014 / 2018 printed in `tables/prod_stations.csv`
  (2011: 13; 2012: 14; 2013: 14; 2014: 12 before WISH/NLWA; 2018: 10).
- Amendment T written ~8:52 pm PT (production fix, 2011 known-answer check, relay rule).
- `scripts/prod_queue.sh` restarted 8:51 pm PT (finished stages skipped; 2011 cont/bank done 8:53 pm). `scripts/cat_relay.py`
  + `scripts/relay_queue.sh` (waits for PROD_2010/2016/2017/2018_DONE; markers RELAY_<Y>_DONE) launched 8:56 pm PT.
- ~9:30 pm PT — PROD_2011 checked: ledger 7/7, 2011-08 episode 41 tier events, controls 0.69/yr, completeness 0.79/0.71 at
  Mw 3.8 (n = 800), between 3 → Tier A′ 0. Queue left running (2012 in progress). RESULTS §13.

## 2026-10-16 (Fri) ~3:10–3:25 am PT — production 2010–2018 done, summary
- PROD_2010…2018_DONE (last ~3:09 am PT); relay_queue started 3:09 am PT (2018 bank built: see `tables/prod2018_relay_*`;
  markers RELAY_2018/2017/2016_DONE). No FAIL after the 8:02 pm fix.
- `scripts/cat_prodsum.py` → `tables/prod_years.csv`, `prod_episodes_all.csv`, `prod_ledger_all.csv`, `prod_between_all.csv`,
  `figures/prod_episodes.png`. RESULTS §14 written.
- Headline: 163 tier events in 51 episodes; 141 in the 6 patch-active episodes (126/129 LFE-carrying); Tier A′ 0 in all 9
  years; ledger 49/49; completeness good 2011–14, poor 2015–18; scrambled > 1/yr 2014–17.

## ~3:20–3:35 am PT — relay results, Amendment U (relay re-run) + V (plan), flagged episodes
- Relay verdicts (coordinator relay + tables): 2018 accepted (15/17, 0.34/yr), 2017 rejected (0/3), 2016 none.
- Amendment U (~3:22 am PT) before any re-run number; relay-cluster injection sources (`cat_inj.py relay_sources`;
  `data/inj/sources_pre_relay.npz` kept); PROD_RELAY mode in cat_pilot (copies `data/*_before_prodr.py`);
  `scripts/prodr_queue.sh` launched 3:23 am PT (markers PRODR_<Y>_DONE).
- Flagged episodes (`scripts/cat_flag.py`): 2010-08b1 "no LFEs" traced to a second lfe_lists clip (lower bound 2011-05-31) →
  fixed, relabel queued (`scripts/relabel_queue.sh`, RELABEL_DONE after PRODR_2016_DONE). 2017-03 events look like real
  patch VLFEs (clean near pulses, far refs quiet); Q2's GNW clause rejects them (GNW noisy in 2017).
- Amendment V written (plan only): forward relay chain 2019–2026 and backward relay 2007–2009 via 2010-08b1.

## 2026-10-16 ~5:55–6:15 am PT — relay re-run + relabel results, final summary
- PRODR_2018/2017/2016_DONE, RELABEL_DONE (last ~6:00 am PT); no FAIL.
- `scripts/cat_prodsum.py final` → `tables/prod_*_final.csv`, `figures/prod_episodes_final.png` (2016–2018 from prodr);
  non-final tables refreshed after the relabel. RESULTS §16 written, §14 marked partly superseded.
- 2018 relay re-run: C0/C1 completeness 0.28 → 0.56 at Mw 3.8; controls 0.67/yr; 2018-06 65 tier (17 relay self, 48 new).
  2016–17 relay re-run: scrambled 3.7/2.65 per yr (fail). 2010 T2a 0 → 0.80, 2010-08b1 8/8 with LFEs.

## 2026-10-16 ~6:05–6:20 am PT — Amendment W (control failure 2016–17, 2016-02-25 candidate)
- Amendment W (~6:05 am PT) before diagnostics; W1 diagnostics (`scripts/cat_diag.py`): LRIV broken 2016–2017 (420–480 nm/s);
  control tier detections cluster on LRIV-broken days; sign-leakage not the cause; W2 candidates not adoptable as written
  (glitch clause fails everywhere; year-own threshold not specific) → W2b (RMS-only station-month QC), disclosed as
  post-diagnostic (~6:08 am PT). QC mode in cat_pilot (env PROD_QC; masks data, controls, noise templates, injections).
  `scripts/w_queue.sh` launched ~6:15 am PT.
- W3: 2016-02-25 12:11 PST candidate = noise in raw waveforms (coherent wavefield everywhere, far refs 42–62 nm); passes Q2
  only because broken LRIV sets the near maximum. `scripts/cat_flag.py` now takes a time prefix and a run name.
- W4 (forward relay 2019–2026) waits for the W re-runs; 2007–2009 held.
- ~7:45 am PT — w_queue bank step for prodrq2016 stuck 60 min (infinite noise-window search under QC masks). Killed PID 374447,
  fixed (QC-excluded channels not required; mask cache; try cap), queue restarted 7:43 am PT (bank now < 1 s).
  2013 W-check: tier 18 → 19, completeness at Mw 3.8 0.831 → 0.831 / 0.826 → 0.779, controls 0.725 → 0.725 → W2b adopted.

## ~9:00–9:20 am PT — forward relay chain 2019–2026 launched (Amendment X)
- Worker counts now env-configurable (CAT_NW in cat_pilot/cat_e2/cat_inj/cat_scanv2; scan_year_queue W=${CAT_NW}).
  w_queue restarted 9:14 am PT with CAT_NW=4 (its prodrq2017 inj_list stage re-runs; finished stages skipped).
- Multi-year relay banks in cat_pilot (PROD_BANK_YEARS, PROD_NAME; prefixes R<y>, members <y>:), cat_relay/cat_inj/cat_prod
  generalised to run names; OTR/HDW added as relay-capable stations; per-year W2b exclusions (`qc_exclusions.py YEAR`);
  ComCat 2019–2026 regional lists fetched. Pre-edit copies `data/*_before_fwd.py`.
- `scripts/fwd_queue.sh` (CAT_NW=4) launched 9:17 am PT; state files `data/fwd_bank_years.txt`, `data/fwd_fails.txt`.
- ~11:15–11:35 am PT — post-W2b diagnostics (W_DONE 11:10 am PT): W2b zeroing verified (data, noise templates); no other
  broken channel; 2011-bank 2017 excess not significant (p ≈ 0.3); relay excess scales with bank size (2017: same 21 bank-pair
  rows + 21 relay rows). RESULTS §18. Fix proposal to Merlin (not applied). Forward chain untouched: prodf2019 in inj
  (inj_list took 38 min).

## 2026-10-16 ~11:20–11:40 am PT — Merlin review 4 → Amendment Y
- Forward chain held (fwd_queue.sh killed after RELAYF_2019; scan-year 2020 build left running). 2019 bank-dependent outputs moved
  to `data/prodf2019_preW2bBank/` + table copies `*_preW2bBank*`; markers renamed `.preW2bBank`; fwd_fails reset to 0.
- Amendment Y (~11:28 am PT). `scripts/y_queue.sh` (2018 W2b + relay re-cut; Y_2018W2B_DONE), `scripts/y_checks.sh`
  (Y3/Y4; Y_CHECKS_DONE) launched. RELAY2018_W2B env switch for the 2018 relay paths in cat_pilot.
- Y5: 2017-03 PNSN tremor within 30 km of C0/C1 2,506 (as in patch majors) vs 42 SVI LFEs → 2017 T2 untestable. Y6: LRIV Z
  fault red (×435 at 200 s, ×79 at 33 s, ×5 at 5 s) → sensor fault.
- ~12:06 pm PT — user decision recorded (Amendment Y7): noise-template control is primary; scrambled = column; relay acceptance = LOO ≥ 85 % only. y_checks updated (noise-template raw + Tier A′ rates).
- ~12:10 pm PT — 2018 W2b relay bank accepted (18 vetted, LOO 18/18, scrambled 0.22/yr, real ETS 86, between 1;
  noise-template tier rate 0). rsources failed on a cont path (QC/relay runs reuse data/prod<Y>/cont) → fixed in cat_inj;
  y_queue resumed (rsources done, prodrq2018 running), y_checks and y_ranked relaunched. Amendment Y8 (ranked catalogue;
  `scripts/cat_ranked.py`, `y_ranked.sh` → Y_RANKED_DONE) pre-registered ~12:07 pm PT; score fit running.
- ~12:15 pm PT — Y8 score fit done (`data/ranked/score_model.json`): 32 positives (known VLFEs found at ≥ 5.5 MAD) vs 140,368
  noise-template detections 2011–2013; 5-fold CV AUC 0.9997. Caveat for the write-up: the positives were selected by the
  battery and two-template agreement, so separation on battery columns is partly built in; with 32 positives some
  coefficients are unstable (e.g. far_missing +2.0, zmad −0.5). The curves (real vs noise-template per catalogue) are the
  product; the score is only an ordering.
- ~1:10 pm PT — Y_2018W2B_DONE: 2018 relay (W2b) accepted 18/18; prodrq2018 T1 83, T2 pass, T3 2/3, scrambled 0.34/yr, noise 0, completeness 0.68/0.60 @Mw3.8. Forward chain restarted (bank 2011 + W2b-2018 relay, capped; CAT_NW=6); fwd_queue patched (RELAY2018_W2B=1, cap).

## 2026-10-16 1:40 pm PT — Y8b launched (fix of broken Y8 ranked score)
- Y8 score broken (missing battery columns raised score). Fix pre-registered as PLAN_CAT Y8b; scripts/cat_ranked2.py +
  scripts/y_ranked2.sh; cat_e2.run_battery got gate_all=True. Smoke test on 12 rows of 2018 passed (scratch only).
- ~17.6k agreement rows (real 8,974; noise templates 1,844; scrambled 4 draws ~6.8k) through the scan-v2 battery, 2 workers.
- Marker Y_RANKED2_DONE; per-year data/ranked/R2BAT_<Y>_DONE. Resume: `nohup bash scripts/y_ranked2.sh &` (skips done years).

## 2026-10-16 ~2:05 pm PT — Y7 rule wired into the verdicts (coordinator request)
- cat_r.py controls.pass = noise-template rate <= 1 per template-yr; scrambled kept as columns (ctrl_rate_between_per_yr,
  scrambled_le_1_per_yr). cat_pilot.py summary items.controls same rule (takes effect for later summary steps).
- fwd_queue.sh: relay acceptance already read cat_relay's Y7 verdict (LOO >= 85 % only); the 2-consecutive-fail rule uses
  completeness + acceptance, never the scrambled gate. Only the header comment was changed (same byte length; the queue runs).
- prodf2019 R verdict re-derived (old copy tables/prodf2019_R_verdict_v2_preY7rule.json): T1 36 pass, T2 pass, T3 3/3,
  noise-template 0/yr, scrambled 1.97/yr (column), between raw 0 / Tier A' 0, GO = True. T1's lambda still uses the
  scrambled ETS-window mean (Amendment R; unchanged).

## 2026-10-16 ~2:40 pm PT — Y8c + RESULTS §20 draft
- Y8c (PLAN_CAT, pre-registered before its numbers, flagged as a post-data change for Merlin): CV-L2 monotone logistic,
  logit as score; CV still prefers the weakest L2 (separable), but the logit spreads. Outputs catalogue_ranked3_*, ranked3_*.
- §20: the 4–8× between-ETS real excess is on any footing (noise ×11, scrambled A1/A2). Single-template raw rates real ≈ noise;
  agreement is 2–6× easier for the 11 similar real templates. Top 5 between tierA rows = 4 LRIV single-station transients +
  1 PGC spike; top days = strong 20–50 s wavefield (storm metric at 5–10 s misses 2013-01-18). 0 of 80 tierA rows pass battery′.
- Bug caught: waves() shadowed its `pre` argument, so top-day figures got numeric names; renamed to r20_topday_wave_*.png.

## 2026-10-16 ~3:00 pm PT — Amendment Y9 (Merlin review 5)
- Pre-registered Y9 (PLAN_CAT) before numbers. Steps 1, 2, 3, 5 done (Y9_STEP1/2/3/5_DONE); step 4 matched noise banks
  running (scripts/y9_nb.sh, per-stage markers data/y9nb/<run><Y>/STEP_*_DONE, final Y9_STEP4_DONE).
- Q2: injections 93–100 % pass vr ≥ 0.12, but known VLFEs Q2 3/26 (GNW clause 11/26, vr 16/26). Ladder level 5 between:
  39 real (29 days) vs noise×11 1.8, A1 2.5, A2 4.0; ETS 146. Eye labels of the 39: wavefield 22, eq-wave 8, spike 4,
  candidates 5 (2016-06-13 07:41 PDT clearest; 2015-04-13 06:06 PDT regional-looking). Emergence threshold 2.04 (85 % of
  known); 2014-01-23 rows 4.5 % above.
- Mechanical checks failed and reported: 24 head2 rows in the top 50 between; known VLFEs at max level 21/26.
- Bug note: loso fill changed 0 → −1 in cat_ranked2.feats (Y9).

## 2026-10-16 ~3:25 pm PT — Amendment Y10 + forward chain restarted; relay2018 injection diagnosis
- FWD_DONE (3:10 pm) came from the old stop rule (2019, 2020 "failed" on combined completeness). Y10: no stop rule, flag only.
  FWD_DONE renamed FWD_DONE.stopRule2020; fwd_queue.sh restarted 3:12 pm PT (CAT_NW=6) at scan year 2021, bank 2011 ∪ 2018-W2b relay.
- Relay2018 injections at Mw 3.8: MF agreement is fine (0.96–1.00 in 2019/2020), the loss is battery item f only (pass_f 0.86
  in 2018, 0.33 in 2019, 0.00 in 2020) with scan z_best 560–4,950 (C0/C1: 34–77). Cause: the relay scan-side injection
  sources (`data/inj/sources.npz`, keys RC*|sc|*) were stacked from `scan_2018/bp_2018v2`, which has no W2b mask; UW.LRIV
  is W2b-excluded all of 2018 (Z all months) and all 18 relay events are in July 2018, so the sources carry LRIV at
  4.9–13.5 µm (C0: 0.11 µm), and UW.WISH at 1.0–2.8 µm (C0: 0.08 µm). Injected into a healthy 2019/2020 LRIV, one station
  carries the scan stack, so leave-one-station-out fails. An injection-source bug, not blind templates. Proposed fix (not
  applied; report first): rebuild relay |sc| (and check |mf|) sources without W2b-excluded channels and with a WISH QC check,
  then re-run inj + summary for prodrq2018, prodf2019, prodf2020 and the chain years.

## 2026-10-16 ~3:40 pm PT — Y10b relay-source fix
- cat_inj.run_relay_sources masks W2b-excluded channels per member month; 2018 RC sources rebuilt (32 LRIV keys gone, rest identical).
- scripts/y10b_queue.sh queued (waits for Y9_STEP4_DONE): inj/summary/R/psum for prodrq2018, prodf2019, prodf2020; marker Y10B_DONE.
- WISH July 2018 passes W2b (0.56–0.73× ref); its big source amplitudes are horizontal long-period noise; reported only.

## 2026-10-16 ~4:00 pm PT — prodf2021 bank KeyError fixed
- run_bank weights: r.get(ch, NaN) (UW.OTR, new in 2021, has no data in the 2011 template windows) -> weight NaN -> channel excluded for that template, same as the bank .get for missing template channels. The chain had stopped (exit 1, nothing skipped); restarted 3:58 pm PT, bank done, mf running.

## 2026-10-16 ~5:05 pm PT — 2021 weak-year diagnosis (existing files only)
- Station set 2021 vs 2020: CN.SNB gone (83 days in 2021, none after; trace archives too), UW.OTR added (no 2011-template data).
- W2b 2021: CN.VGZ Z excluded Jul–Dec (MF-band noise 3.6×, 8.4, 11.2, 11.4, 9.6, 8.9× ref; E rises to ~1.0×, N fine);
  UW.SP2 E all year, N Jun–Dec (as 2020). All others normal (0.4–2.3× seasonal).
- 2011-bank weight shares (2020 / 2021): LRIV 0.16/0.18, VGZ 0.13/0.15, PGC 0.12/0.13, DOSE 0.11/0.12, WISH 0.10/0.12,
  SNB 0.10/–, GNW 0.09/0.10, NLLB 0.09/0.10, SP2 0.09/0.10, OTR 0/0. From July 2021, VGZ Z (the patch-nearest vertical) and
  SNB are both gone, and PFB/YOUB/LZB have been absent since 2019.
- Patch LFEs (svi_c010_075) per episode window: 2021-01 1,486; 2021-05 547; 2021-07 162; 2021-08 85; 2021-10 2,239 (2020-11:
  1,842). The 2021-10 ETS did reach the patches, so the low 2021-10 completeness (0.05) is instrumental.
- VGZ beyond 2021 (raw traces, Z 20–50 s, 12 days/yr, VGZ/PGC ratio): 2020 2.0, 2021 2.8, 2022 2.5, 2023 3.0, 2024 2.8,
  2025 1.0, 2026 1.0, so it is likely degraded 2022–2024 and recovered (or re-instrumented) 2025–26; W2b decides monthly.
  SNB has no data after March 2021.

## 2026-10-16 ~5:20 pm PT — (c) why 2021-10 completeness is 0.05 (prodf2021 injections, Mw 3.8; diagnostic only)
- Battery item f (leave-one-station-out on the scan) is the collapse: C0 pass_f Jan 0.44, May 0.88, Jul 0.00, Aug 0.03,
  Oct 0.01; C1 ~1.0 until Aug, then 0.40 in Oct; all relay sources ≤ 0.06 from July. Items a, c, e′ stay 0.94–1.00, b 0.8–0.95.
- Cause: the scan's near network in 2021 is only NLLB, PGC, SNB, VGZ, LRIV (`tables/scan2021_noise_month.csv`); SNB ends in
  March and the scan noise rule excludes VGZ Jul–Dec, so from July the scan rests on PGC + LRIV (+ the distant NLLB) and
  dropping either one kills the scan z (C0 z_best ~32 but z_loso < 0.5 z_best).
- Second, smaller loss: MF agreement. Median best/second MF z fall from ~8.0/7.5 (Jan–May) to 6.4/5.7 (Oct) for C0 with
  24 live channels in Oct (27 before), so C0 agreement drops 1.00 → 0.61 and some relays 0.08–0.5. The 5.5 MAD threshold is
  per template-year, not per day; the loss is lower SNR from losing VGZ Z (and SNB), not a moving threshold.

## 2026-10-16 ~5:08 pm PT — Y9 step 4 verdict (accepted banks) and Y9_DONE
- Accepted 14/24 banks: level-5 between 2.58 (1.62–3.91) per bank-yr → 14.3 (9.0–21.7) scaled to 2010–18 vs real 39
  (27.7–53.3); 2011–13 real 11 (6.1/yr) vs expected 4.6, p = 0.008. Outcome intermediate (neither A ≥ 20 nor B ≤ 5).
  All-banks column 46.2. Agreement ratio 0.10–44.6 follows the base window (2012 N02 base = 26.6k detections/yr).
- 11 real 2011–13 level-5 rows: 7 wavefield, 2 eq-wave, 2 spike, 0 candidates. cat_y9nb.summary now uses accepted banks;
  old verdict kept as y9nb_verdict_v1_allbanks.json. RESULTS §21 completed; Y9_DONE written.

## 2026-10-16 ~5:18 pm PT — Y9b candidates, Y11 f′ re-score, station query (RESULTS §22 draft)
- Candidates: 2015-04-13 = regular earthquake (HF burst at −100 s; pre-registered HF window missed it, disclosed);
  2016-06-13 07:41 = VLFE-like; other three unclear.
- f′: controls rise ≤ 1, C0 completeness +0.38 in low-f months (2021: 0.15 → 0.64). Numeric adopt rule met; held for step 1
  and because FDSN shows healthy broadbands near the patches (BUCKS 24 km from 2023-09, SYMB 32 km since 2016; not local).
- Correction: svi_c010_075 = LFE candidates; true patch LFE counts 2021-01 0, 2021-10 55.
- Compute queue scripts/y11_queue.sh (after Y10B_DONE): inj38 → ladder 2019-21 → Y9b banks; markers Y11_INJ_DONE,
  Y11_LADDER_DONE, Y9B_DONE.

## 2026-10-16 ~5:40 pm PT — Y12 started (stations first)
- Fetch running (scripts/fetch_y12.py; SYMB 2021→, BUCKS 2023-09→; ~0.45 MB/day stored, ~24 MB/day transferred; ETA
  FETCH_Y12_DONE ~6:35 pm PT). Responses for SYMB, BUCKS, RNWY merged into data/resp (originals *.preY12).
- Code hooks (env-only, chain unaffected): Y12_EXTRA (cat_scanv2 extra scan stations), INJ_SOURCES (cat_inj sources file),
  LADDER_RUNS (ladder run override). cat_y12.py (qc, GF sources, mkrun, rescore); y12_pilot.sh (2021 as prodfy2021).
- Master queue scripts/yq_master.sh (≤ 2 workers, after Y10B_DONE): Y11 inj38 → Y12 pilot 2021 → Y11 ladder 2019–21
  (2021 = prodfy) → Y9b banks. Markers Y11_INJ_DONE, Y12_PILOT2021_DONE, Y11_LADDER_DONE, Y9B_DONE.

## 2026-10-16 ~7:45 pm PT — Y12 pilot evaluated (RESULTS §23); prodfy 2022–26 started; pre-Y12 chain stopped
- QC pass 2021: SYMB, OSQM, RNWY, SLDQ. C0 2021-10 f 0.006 → 0.36 (×1; ×0.5 0.07, ×2 0.49); f′ ~0.47 regardless.
  Real: T1 2 → 5, T3 2/3, 2021-10 one tier event; no relay bank (2 vetted). f′ not applied (f pass 0.55 at ×1).
- Y11 step 1: C0 drops LRIV 99.4 %, relay drops PGC 99.0 % (accept).
- Pre-Y12 fwd_queue killed 7:37 pm PT after prodf2023 inj_list; prodfy_queue.sh started (RELAY_RUN_PREFIX=prodfy, MF seeded
  from prodf2022/2023). Marker per year PRODFY_<Y>_DONE, FWDY_DONE.

## 2026-10-16 ~8:08 pm PT — Y11_DONE (RESULTS §24)
- 2019-21 level-5 between: 14 real vs noise 0, A1 1.5, A2 3.0, matched-bank placeholder 4.0; eye: 9 wavefield, 1 eq-wave, 1 unclear, 3 candidates (2020-06-11 13:42, 2020-09-11 05:00, 2020-10-25 19:25 PDT). 2021 per-episode table written; none bound-only on the Y12 network.

## 2026-10-16 ~8:15 pm PT — Y13 queued; prodfy 2022 done
- Y13 pre-registered (PLAN_CAT); scripts/cat_y13.py + y13_queue.sh (after Y9B_DONE, 2 workers; 2021+ years wait for prodfy MF/battery). Marker Y13_DONE.
- prodfy2022 done 8:11 pm PT: Y12 stations OSQM, RNWY, SLDQ; CN.SYMB failed 2022 QC only on the teleseism ratio (2.25 > 2, M7.6 Mexico 2022-09-19 at 3,812 km; level ratio 0.99) -- rule kept, reported. No 2022 relay bank.
- ~8:22 pm PT: Y12b (teleseism ratio = column only; PLAN_CAT). prodfy2022 SYMB re-run queued in lane2 (after Y9B_DONE, before Y13); first run kept as no-SYMB variant.
- ~8:35 pm PT: RESULTS §25. Y9b verdict "A" was a pooled-row bug (2010-2021 live time); corrected 2010-18 control 17.2 (13.3-21.9) vs real 39 -> intermediate; real/control 2.8x in 2011-17. cat_y9nb summary fixed (no re-run).
- ~8:45 pm PT: Aug-2022 blindness diagnosed (tables/y12_2022_episode_diag.csv): f ratio clause fails seasonally (C0/C1 Mw 3.8 max station share 0.62-0.69 Jan-May -> 0.84-0.90 Jul-Oct, LRIV dropped; z_loso 18-29 so f-prime passes 1.00 every month); VGZ Z data only 12/31 days in Aug 2022 (gaps Jul-Oct); agreement 0.70-0.75; e-prime fine. Patch LFE075: 0/0/0/1/0/0; tremor <=30 km: Feb 2, Mar 139, May 5, Aug 218, Aug-b1 0, Dec 16.
- ~11:00 pm PT: Y14 pre-registered; f-prime in cat_common (apply_rules + rescore_inj; pass_f_preY14 kept); y14_rescore.sh running (marker Y14_RESCORE_DONE ~11:20 pm); cat_y14.py prints the per-year f vs f-prime table afterwards.
- ~11:15 pm PT: RESULTS §26 table written (y14_summary.csv). No GO changes; completeness up (2019-23 +0.11..+0.37); relay LOO under f 2019 0.72, 2020 0.36, 2023 0.56 (not re-run).
- ~4:40 am PT: RESULTS §27 (Y13 verification). C0 46 between / 220 ETS vs g1 63-130 / 57-101: g1 not ETS-modulated; top-33 by eye mostly earthquake waves (regional, far refs first); no tremor/ETS-edge clustering. Y13 "A" not supported -> isolated, unconfirmed.
