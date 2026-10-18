# MERLIN — review after the 2014 NO-GO (2026-10-15 ~3:30 pm PT; saved by main session, condensed)
Verdict: method viable; 2014 failed as a TEST YEAR, not the detector. (1) The 2014-11 episode never reached the template patches: P≥0.75
LFEs within 10 km of clusters C0/C1: 2014-11 8/3 vs 2011-08 121/83, 2012-08 157/205, 2013-09 92/20, 2015-12 80/61, 2018-06 74/80,
2019-08 85/30, 2020-11 66/45, 2021-10 28/27, 2023-04 81/52, 2024-08 66/67, 2025-09b1 11/9. (2) "4-station early era" premise wrong:
7–10 near broadbands on disk in 2003–05 (CN PGC/VGZ/SNB/NLLB/LZB + C8 SHVB/KHVB/MGCB/TWBB/SHDB) and 2007–09 (+UW SQM/HDW/GNW, XU BS11/W040).
(3) Plain channel sum diluted by horizontals ~20× noisier than verticals (adding GNW lowered cn4 2/8→1/8).
Route: patch catalogue with era-matched empirical banks, 2010–2026 (2007– if era2007 passes); pilots 2013-09 then 2018-06; 2003–06 bootstrap later, gated.
Flaws: P2 GO rule gated on between-ETS count (must be a result, not a gate); 2014 bank included H&G T6–T13 which are not observable
(behave like noise templates) → recount without them: between 2, ETS 0, other 1; BUG hg_r3.lfe_lists() clips LFE/tremor lists to
< 2013-01-01 → all 2014 LFE labels void; PGC "700–1,000 nm" claim inconsistent with noise tables (PGC Z 21–50 nm/day) → resolve/retract;
2014 used 3 far refs (MSO/WVOR unprepped); E3 "ETS 45–49" includes self-detections (honest benchmark: 15 tier events in 2012-09 with 2011 templates).
Prescription:
0 fixes (~1 h): lfe_lists frame-year; prep 5 far refs per pilot year; validated bank = T1–T5 + C0–C5 (+T2/T4 re-cuts), others columns;
  GO rule: between-ETS count reported, never a gate. Accept: 2014 re-summary ETS 0 / between 2 / other 1.
1 station-era table + 2011–12 LOO re-test with sets era2007 {PGC,VGZ,SNB,NLLB,GNW,SQM,HDW}, era2016 {PGC,VGZ,LRIV,SNB,PFB,NLLB,DOSE,GNW,SP2,HDW,SQM,OTR},
  max2011; pre-registered weighted sum (per-channel weight = template RMS / pre-origin 10-min RMS). Accept: weighted ≥ plain; hg9 ≥ 7/8;
  era2016 ≥ 7/8 with scrambled ≤ 1/yr; era2007 ≥ 5/8 with scrambled ≤ 1/yr (→ 2007–09 joins; 2008-05 pilot).
2 triage the 3 2014 between-ETS events (pre-stated class: real-looking iff vr_thrust ≥ 0.12 AND no comparable far wavetrain AND
  ratio_GNW ≤ 0.551 if north of 48.1 N; expectation: all "likely noise"); vr_thrust ≥ 0.12 as a column for the next pilot.
3 Pilot 2013-09 (same network as E3): truth T1 ETS tier ≥ 5 & Poisson p < 0.01; T2 ≥ 70 % with LFEs and ≥ 3× random; T3 near-station pulse;
  controls 8+8 scrambled + noise-pair ≤ 1/yr; completeness ≥ 40 % Mw 3.4, ≥ 70 % 3.6, ≥ 80 % 3.8 (C0/C1 injections ETS and between).
  Fail with good completeness = finding (VLFEs did not repeat).
4 Pilot 2018-06 (era2016 bank), same rules.
5 Full run 2010–2026 (2007– if era2007 passed), ~1.5–2 days, per-year Mw_c, per-episode counts/moment, patch-recurrence table,
  between-ETS Poisson bounds, known-answer ledger.
6 2003–06 bootstrap, gated rehearsal on 2012-09 with 4-station cascade; else bound only.
Honest claims if pilots fail: 18 strong SVI ETS VLFEs 2011–12 at LFE patches (0 in 26 controls); non-repeat in 2013/2018 at stated
completeness; no real between-ETS VLFE 2011–12 + Poisson bounds; H&G 2014 templates not observable; 2014-11 did not activate the patches.

## ADDENDUM (Merlin, density-aware revision, ~3:45 pm PT) — supersedes the station sets/pilots above where different
Broadband density on disk (sensor-checked, ≥ 30 s flat): ≥150 real days within 160 km: 2003 ~9, 2004 ~17, 2005 ~13, 2006 6, 2007 15,
2008 15, 2009 10, 2010 13, 2011 14, 2012–2014 19–20, 2015–2020 15–18, 2021 21, 2022 27, 2023 31, 2024–26 33–34 (VLFE runs used 5–13).
In ETS windows: 2011 template week 24 stations; 2012 22; 2013-09 24; 2018-06 16; 2007-01 24; 2008-05 26; 2009-05 12–13; 2010-08 26;
2004-07 16 (POLARIS, all north of strait); 2005-09 12; 2003-02 6–8. CAFE BS11/W020/W040/PL11/W060 recorded every 2011 template day →
2011 templates can be cut on 11–13 stations that exist in 2007–2010. Patch activity early: 2003-02 130/222; 2004-07 627/1325; 2005-09
460/634; 2007-01 43/7; 2008-05 147/87; 2009-05 68/51; 2010-03 0/0; 2010-08 109/58.
What density buys: decay-ratio separation 2σ (k=3 ring) → 3.5σ (k=10) → 5σ (k=20); LOO z scales ~√N_ch with a noise-weighted sum
(predicted ~8 at 33 ch, 11–13 at 70 ch vs saturating 4.6–5.4 at 12 ch); Mw_c −0.14 to −0.15; 2007–2010 rescued by CAFE/TA (11–13 carriers).
Revised plan: 0 fixes; 1 dense bank + weighted sum on 2011–12 (sets dense2011/12, era2007, era2010, era2016, hg9, cn4; plain vs weighted;
  accept weighted ≥ plain, hg9 ≥ 7/8, dense ≥ 8/8 scrambled ≤ 1/yr, era2007 ≥ 5/8 scrambled ≤ 1/yr → 2007–2010 enters);
2 triage 2014 events; 3 pilot 2013-09 (dense, 24 stations); 4 pilot 2008-05 (CAFE era, 26 stations; backup 2010-08) instead of 2018-06;
5 full run 2007–2026; 6 2003–06 gated scan-only bootstrap rehearsal, then POLARIS-era templates.
Optional downloads (need user OK): (i) TA 2006–2008 LH 1-sps (A04A, B04A, B05A, C03A, C04A, C05A, D03A, D04A) ≈ 3 GB;
(ii) southern 2003–2009 ETS + control-month windows (UW.OPC, SQM, GNW 2003–06, US.OCWA) ≈ 15 GB.
