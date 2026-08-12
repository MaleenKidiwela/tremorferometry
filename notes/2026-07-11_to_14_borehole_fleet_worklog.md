# Borehole LFE-coda dv/v fleet — work log 2026-07-11 → 2026-07-14

Narrative + decision record for the multi-day push that produced the **27-station borehole fleet**.
The authoritative pipeline + live status ledger is `notes/FINAL_PIPELINE.md`; this is the journey, the
decisions, and the test results behind it. Advisor throughout: **Merlin** (many consults) — several of its
catches are logged below.

## OUTCOME
**27 strong stations (≥20 causality-certified) · 2,068 certified LFE families · + 3 flagged** (kept,
excluded from inversion). The full 45-station PB borehole pool is exhausted; the rest are closed with
documented reasons. Next phase = broadbands (`notes/2026-07-14_broadband_dvv.md`), then the 4-D inversion.

Strong: B011(119) B035(113) B926(110) B943(109) B928(106) B024(100) B001(98) B935(98) B005(97) B036(96)
B012(87) B039(86) B010(85) B003(83) B040(80) B018(78) B031(73) B004(64) B033(55) B927(55) B017(51) B023(49)
B032(45) B013(42) B201(39) B941(26) B022(24). Flagged: B045(4) B009(10) B027(13).

## TIMELINE
- **07-11 (crash recovery):** server crashed mid-run. Reconstructed state from disk; resumed. Reworked the
  orchestration to a **rolling buffer** (delete each station's traces after finalize) + held downloads
  ("downloads way ahead"). Reinstated then investigated the **battery gate**.
- **07-12:** built the autonomous **continue-through-29 driver** (throttled download feeder + serial
  processor + supervisor). Fixed the **300-cap** and **adaptive-floor** bugs. Drove the fleet toward 29.
- **07-13:** fleet at ~25 strong; low-cal tail stations starved (B045/B030/B932). Merlin audits caught real
  bugs (below). Froze the **inclusion rule before the inversion**. Endgame reframed count→quality.
- **07-14:** southern **Ducellier catalog** → **B935 rescued** (98, catalog-confirmed); **B017 recovered**
  (51, a healthy station that had fallen out). Fleet = 27 strong. Mendocino declared **source-free**.

## KEY DECISIONS
1. **Battery gate: REINSTATED then DROPPED.** Ran the Tier-1 stack-vs-random battery on the causality-
   certified subset; it **can't separate real from ringing on cap-off data** (B009 5/10 FAIL, B010 67/85
   FAIL, 0 TRUSTED anywhere; ringing coda_sigma ≥ certified). Causality is the certification; battery is
   informational only. (`battery_gate.sh`/`battery_gate_dvv.py` kept but not called.)
2. **300 = DENSIFY BUDGET, not a result cap.** Certification is per-family, so dv/v uses ALL causality-
   certified families (never truncate to 300). Nearly recapped B005 to top-300 (would have discarded 37
   reliable families) — user caught it: "isn't the goal dv/v from reliable LFEs?" Correct.
3. **Adaptive candidate floor 0.5 → 0.2.** The 0.5 floor starved low-calibration stations to 0 families
   (B045 got 1046 cand, B030 got 34) — the B011-trained picker scores real LFEs low off-B011. The 300 cap
   bounds densify cost regardless. (⚠ first attempt only edited the docstring — Merlin caught it — line 18
   was still `max(0.5,…)`; fixed for real, verified 1046→15000.)
4. **Fleet inclusion rule FROZEN (rule-before-result, before any inversion):** a station enters the
   inversion + pooled stats only if **z_certified ≥ 20 AND causality survival ≥ 15%**; flagged stations
   stay on disk + map labeled, excluded from analysis. Prevents post-hoc fleet tuning.
5. **Composition by QUALITY, not count.** Refused GOLD-1/2 padding stations (B028 is 0.4 km from B027;
   B204 ~20 km from B201 — no new interface coverage). Finalize at the honest count; don't force 29/30.
6. **Rolling buffer + gap audit + map auto-refresh** run per finalize; interior-year <200-day gaps flagged.

## TEST RESULTS / FINDINGS
- **GOLD does NOT predict causality-certified count.** B045 (June GOLD 52) recovered to 15k candidates →
  172 families → 19 LFE-label → **4 causality-certified** (median causality 0.79 = noise-dominated).
  Lowering the floor recovers *candidates*, not *certified families*.
- **B005 "over-selection" was valid, not waste.** 544 spatial bins → 564 templates → 97 certified; capping
  to 300 would drop 37 real families (P barely predicts certification — flat ~15–21% across P quartiles).
- **B935 (98 certified)** decisively overturned its "junk GOLD-6" refusal — it's one of the 4 boreholes that
  *generated* the Ducellier catalog (records 34/66 families). 3 of those 4 (B039/B040/B935) are strong here.
- **Mendocino source-free (catalog-confirmed).** B045/B049/B932 appear in 0 Ducellier families (nearest
  81–109 km); an independent purpose-built catalog confirms no LFEs within range → geometry, not detector
  weakness. Citable network-design limit.
- **2017 gap = a DOWNLOAD truncation** (data only fetched to day 59), not a processing gap; the old "fix"
  re-densified the same 59 days. B011 truly fixed (re-download → 356 days); B926 2017 still pending.

## BUGS CAUGHT (mostly by Merlin — reason the user wants it consulted for decisions)
- Floor "fix" that only edited the docstring (line 18 unchanged) — would have failed every recovery.
- recap300 that would delete 37 reliable causality-certified families to enforce a books-uniform 300.
- The 300-cap not actually enforced (bin-coverage overrode it → B005 564, B011 397, B035 324, B943 394).
- Map silently 2 stations behind (new prepped stations lacked coords in candidate_stations_post2020.csv;
  added IRIS coords for all 23).

## DEFERRED BACKLOG (before/with the inversion)
B926 2017 gap fix; B012 2025 gap; B004 battery backfill (traces on disk, free); B006/B007 mirror-correction
validation triplet (co-located with strong B005). Then the **4-D δβ/β inversion** over the 27 strong stations.
