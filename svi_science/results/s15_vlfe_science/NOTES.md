# S15 NOTES (log)

- 2026-10-17 ~2:15 pm PDT: read MERLIN_SCIENCE_PLAN.md, catalogue RESULTS §14–§27, findings page finding 13. Located the
  540 = sum of `<run>_ets_events_v2.csv` over the 17 runs of the findings-page table (verified row counts).
- ~2:40 pm: coordinator relayed the user's instruction: use the catalogue as closed, state limits once, detection-recovery
  check is a column, not an abandon gate. Put into PLAN §0.
- ~2:45 pm: PLAN.md written (pre-registration) before any number.
- ~2:46 pm: first Q1 run (s15_q1.py, 11 s). Loader checks: recomputed ±600 s counts reproduce n_lfe_10m exactly (431 / 109).
  Episode-window dates are UTC dates (one ETS event at 2022-11-30 03:37 UTC fell outside under a PT reading). pandas
  datetime64[us] trap hit once in the episode loader (fixed with .dt.as_unit('ns')). 15 Bostock-only families lacked
  positions in s0/families.csv → taken from lfe_catalog_repair/results/e5r_families_2012.csv.
  Directions: own LFE-onset sign agrees with s3 on 6/8 shared episodes (disagrees 2019-08, 2023-04, both |ρ| ≤ 0.10);
  s3 used where present. 2011-08 NW (ρ +0.90), 2012-08 SE (−0.93), 2013-09 NW (+0.75), 2025-09b1 NW (+0.88).
- ~2:50 pm: AUDIT of run A: the fixed-weight control is biased toward "late" (VLFE time sits after the latest onset of its
  concurrent families; the control draws from the union). PLAN Amendment 1 adds a re-anchored control (pseudo-VLFEs at
  random LFE times with their own concurrent sets). Run-A numbers kept as reported.
- ~2:53 pm: Q2 first run (s15_q2.py, 14 s): M2 τ_updip a = 1.05 (n 514 with a control), θ0 −42°, matched LFEs a 0.47 (p99 0.60), θ0 −42° → 'stronger' by the §3.14 rule. Amendment 2 adds robustness checks (clusters, detection vs tidal phase, per-episode). 26 VLFEs have no LFE pool (no control); 22 needed the whole-episode fallback.
- ~2:52 pm: re-anchored control (s15_q1b.py): control median 1.16 d vs fixed 1.145 d — the audited bias is small.
  9 VLFEs (7 in 2025-09b1) have no eligible pseudo-VLFE in their pool (all pool LFEs isolated) and drop out of the
  re-anchored test (422 left). Verdict PASS-late under both controls.
- ~2:53 pm: detection column (s15_detect.py): recovery days 0–2 vs ≥ 2 after local onset +5–7 % (big-9 +1 %), flat.
- ~2:58 pm: Q2 checks (s15_q2b.py): detection vs M2 phase flat (b ≤ 0.012, 95 % ≤ 0.03); cluster bootstrap 0.87–1.22;
  declustered 1 h / 6 h a 0.74 / 0.76 vs LFE p99 0.64 / 0.66; VLFE > LFE median in 9/9 big episodes.
- ~3:00 pm: figures (s15_figs.py) rendered and checked by eye; label overlaps fixed once.
- ~3:05 pm: RESULTS.md written. Self-check corrections before hand-back: "6/9 episodes above the LFE 99th pct" → 5/9;
  random-time a "≤ 0.1" → mean 0.08, p99 0.21; the no-2016–17 row is the fixed-weight control; S3 return sentence rephrased
  to what was computed.
- Not done (not asked): findings page / tree update.
- ~3:12 pm: Merlin audit received (MERLIN_AUDIT_Q1Q2.md): Q1/Q2 stand; gate regression + O1 + columns pre-registered as PLAN Amendment 3 (~3:15 pm) before computing.
- ~3:16 pm: P ≥ 0.5 Q1 column (S15_PMIN=0.5 → tables/p050/): PASS-late (2.51 vs 1.88 d). NW-ward single-family events
  stay flat (1.42 vs 1.62 d, pct 18) — Merlin's "would change my mind" item partly met; reported as a caveat.
- ~3:22 pm: s15_q2c.py (4 min; 500-resample day-block bootstrap ×6 models): gate ACCEPT at threshold (ratio 2.01,
  CI 1.32–2.78; without 2016–17 2.00 → intermediate); β 1.36, β·A_L 0.69 inside the VLFE range. O1 not detected.
  New: VLFE a rises with Δt (0.91/1.11/1.37), not with burst count or size. LFE baseline without own ±1 h 0.30.
  First attempt crashed on object-dtype epoch from Q.build (cast to float).
- ~3:25 pm: RESULTS.md updated (plain answer rewritten, Amendment 3 section added, interpretation reworded per Merlin).
- ~3:52 pm: Merlin Q3/Q5 design received; pre-registration pasted into PLAN as Amendment 4 before any number.
- ~3:51 pm: Q3 step 0 re-stacks launched (scripts/s15_q3_restack.sh, 3 parallel, marker Q3_RESTACK_DONE; rows ex/in: 2012-08 2796/68, 2018-06 1462/36, 2024-08 1166/59 pooled LFEs).
- ~3:55 pm: re-stacks done in 4 min (ex/full 0.85 / 0.93 / 0.92 → additive). s15_q3.py: partition table; size test
  uninterpretable (ρ −0.90 with n_sta_day; size ρ −0.14). 2022-03 has 0 zone SVI LFEs (centroid = C0/C1 midpoint; 2022 LFE gap).
  The 2018 relay self-detections are prod2018_relay_events.csv (the design's 'prodrq2018_relay_events.csv' does not exist).
- ~3:56 pm: s15_q5.py (11 s): local-pace productivity flat (0.91, pct 14); daily fast-advance excess (28 % vs 22 %, pct 100);
  RTR excess gone against ±24 h draws (25.6 vs 24.8 %); tide by pace 0.71/0.74/1.00 (LFE rises too). Detection no 30 % flag.
  Caveat: our pool's main-speed sign disagrees with s3 for 2019-08, 2023-04.
- ~4:00 pm: RESULTS.md updated (Q3/Q5 plain answer + sections), figures q3_partition.png, q5_migration.png.
- ~5:23 pm: Amendment 5 pre-registered (Q4 GF injections at 7 nodes + C0/C1 GF references, 2012 and 2024; slip-rate column; b-value substitute).
- ~5:30 pm: GF sources built (data/q4/sources_<Y>.npz). Check vs empirical C0/C1 stacks: zero-lag cc median only 0.15–0.23
  over all channels, but the near verticals that dominate the MF match at lag ±1 s (LRIV.Z 0.87–0.90, VGZ.Z 0.82–0.90,
  PGC.Z 0.72–0.80); horizontals and far channels fit poorly. Amplitude GF/empirical 0.8–0.9 (C0), ~0.3 (C1). So GF sources
  are a rougher stand-in than the empirical stacks; the GF C0/C1 references (G0/G1) are the like-for-like baseline.
- ~5:35 pm: Q4 injections done (2012 5 min, 2024 4 min; 5,400 each, 0 row errors; marker S15_Q4_DONE). SVI nodes
  48.34–48.52 N recovered like the GF C0/C1 references (ratio 0.8–1.1); NW end 48.70 N and S-Olympic 47.80 N < 0.2;
  NE-Olympic 0.04 (2012) / 0.37 (2024). GF references recover ~half of empirical C0/C1 in 2012, equal in 2024.
  s15_q4b.py: centroids (VLFE 8.8 km NW of own LFEs, both directions = bias), slip-rate column, b-substitute.
- ~5:45 pm: RESULTS.md updated (Amendment 5 plain answer + section); figure q4_detection_map.png.
