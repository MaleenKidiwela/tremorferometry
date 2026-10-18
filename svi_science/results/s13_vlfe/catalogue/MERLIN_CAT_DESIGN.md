# MERLIN — design for the 2003–2026 strong-VLFE catalogue (2026-10-15 ~1 pm PT; saved by main session, condensed)
Verdict: viable, but NOT with the synthetic-grid scan as detector (2012 scan-only pool through the R3 battery: real 76/yr vs controls
94–96 between ETS; no separation at any z). What produced 8/115 vs 0/130 was the empirical-template matched filter (H&G T2, 27 ch,
90-s real waveform) as selector + scan/far-refs as vetter. Detector = empirical template bank (H&G T2 + stacks of our vetted events +
H&G 2014-12 templates); vetting = scan + far references + refined earthquake test; thresholds from controls per template-year
(never from Table S1 counts). Pre-2010 needs a restricted-station test on 2011–12 first.
Flaws: F1 scan as detector fatal; F2 test (e) over-broad (6 of 14 ETS real-looking lost to head2 windows of distant M4.6–5.0
teleseisms) → e'; F3 H&G station set complete only 2011–14 (6 stations 2012–26; 2–3 before 2009); F4 far LHZ refs: HAWA 2005-03,
NEW 2006-01, OMAK 2008-11, I02D 2010-08, COR pre-2009 unverified; GNW location-confounded → use ratio_far (n_far ≥ 2) as rule
(reproduces 8/115 & 0/130), GNW a column; F5 V4 thresholds fitted to Table S1 → control-based thresholds for new years; F6 MF
completeness unmeasured → cross-template injection of scaled real waveforms (leave-one-out); F7 battery passes ~0.35 % of scrambled
candidates; full V4 ~2,000 between-ETS dets/yr → ~7 false/yr alone → need a surrogate tier (two-template agreement / higher MF);
F8 live-time denominators (e' removes ~50 %); F9 PNSN tremor only 2010+, sc_gain only 2011 epochs → per-epoch gains or g=1;
F10 the 8 are a post-hoc positive control (independent support = LFE coincidence + near-station decay).
Staged plan:
1 PLAN_CAT.md pre-registered (tiers, thresholds, battery a–f + e', completeness, live time, pass rules).
2 E2 (2011–12, ~1 h): battery on every V4 T1/T2/T3 detection vs 8+8 control catalogues and selection-matched noise, vs MF level,
  with (e)/e', ± S5 columns, two-template-agreement tier. Accept: a tier with ETS ≥ 8, between-ETS ≤ max(control)+2, control ≤ 1/yr.
3 e' (fail if bw5/regional/reg3, head2 of M ≥ 5.5, or any head2 with ≥ 2 far refs excited). Accept: ETS 8 → ≤ 14, scrambled ≤ 2/650,
  between-ETS stays 0.
4 E3/E4 (2011–12, ~2 h): empirical templates from the 8 vetted + 6 near-miss events on all 16 stations, clustered by location; MF on
  station sets {all}, {H&G 9}, {PGC,VGZ,NLLB,SNB}, {+GNW}; controls; threshold = noise-template rate 0.5/template/day; battery; LOO.
  Accept: all-station ≥ 7/8 re-found by a template not containing them; 4-station ≥ 5/8 with control ≤ 1/yr.
5 Far references: FDSN metadata query for LHZ 2003–2026 at 250–800 km (IU, US, UW, TA incl. 2006–08, BK, CC, CN, II); ≥ 3 per year;
  ~57 MB/station-year compressed → ≈ 7 GB for 5 stations × 23 yr (NEEDS USER OK). If 2003–06 < 2 LHZ stations: ask for a new class
  (candidate-window BHZ ±1.5 h, CN.BBB/PNT, ~1 GB/station-year).
6 Pilot 2014 (truth: H&G 2014 T6–T13): GO ≥ 5/8 re-found, between-ETS ≤ control + 1 (controls ≤ 1/yr), ≥ 80 % at Mw 3.8.
7 Pilot 2008 (only if E3 4-station passed); else 2003–09 = Mw ≥ 4.0 tier + bound; uniform catalogue starts 2010.
8 Full run 2003–2026 (~2–3 days wall; ≤ 8 workers); tiers A (controls ≤ 0.3/yr + battery + two-template agreement), B (≤ 1/yr + battery).
9 Outputs: event table (UTC/PT, template, MF, location, Mw, ratio_far, e' flag, LFE/tremor status, tier, waveform fig), per-ETS counts
  /moment/completeness/LFE-locked moment/GNSS, per-year between-ETS count and Poisson bound above uniform Mw_c (~3.6–3.8 for 2010+,
  ~4.0 earlier), known-answer ledger.
10 WA-NW afterwards.
