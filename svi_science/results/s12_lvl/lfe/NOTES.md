# S12-LVL, LFE side — decision log (times PT)

## Request (Tue 2026-10-13 ~9:45 am PT, via main session)
Find in our LFE data the ETS-related LVL velocity change of Gosselin et al. 2020 (PGC/SNB/VGZ: Vs contrast at the LVL top
+0.05–0.1 km/s in the 15 d after the ETS midpoint, decaying; layer ~3.5 km). The user is clear the signal exists: design
observables sensitive to it, then search; do not conclude "nothing" from a measurement that cannot see it. Tasks: (1)
forward calibration first; (2) LVL-base reflection (main idea); (3) superposed-epoch coda dv/v + Qc per lapse window
aligned on the front, with a between-ETS baseline; (4) P/S and S−P across the front; (5) splitting if time.
Existing waveforms/stacks only; ask before any download.

## Why S11 could not see it (from the request, kept as the design premise)
- ΔD compares episodes with neighbours: a change repeated after every ETS cancels.
- The coda averages 25–45 km; direct rays from LFEs go up and miss the LVL below the source.

## Log
(Times before 10:25 am were first written from memory and were ~10 min late; corrected from the log files: forward MC
9:52, coef 9:58, PLAN §1–5 ~10:02, queue 10:07, PLAN §6 ~10:13, queue2 10:18, all PT.)
- 9:50 am — Read: Gosselin text; S11 REASONING/NOTES; round5, redpy, ballistic round3, round7 RESULTS; S10 pilot tables.
  Data found: S11 day-stack stores (`s11_depth/tables/store[_ext]`, P-weighted, 3C, 40 Hz, 2–8 Hz, tau −14..+10 s about S,
  every window RMS-normalised over the 24 s, so amplitude ratios are valid only within one component); round-3
  cc-weighted day stores (`round3/{svi,wa_nw}/tables/store_ccw[_ext]`); round-7 front dates per (family, major)
  (`round7/{svi,wa_nw}/tables/r7_front_times.csv`, set = headline).
- 9:52 am — Geometry: plane fit to Slab2 under the SVI families: dip 10.2° toward N52E, residual 0.6 km; families sit
  3.7 km above Slab2 on median (q10–q90 0.2–5.4 km). LVL top placed at the hypocentre (nominal), with tallies for a 2 / 5 km
  layer, the source 1 km inside, and the top 1 km below the source.
- 9:55 am — DECISION F1: LVL Vs 3.0 km/s, Vp/Vs 2.3 (CVM v1.7 south_VI has no LVL: 4.25 km/s at 34.5–42 km). Why: the
  LVL is ~25–30 % slow in published RF models; Gosselin's +0.05–0.1 km/s is then +1.5–3 %; nominal δ = +2 %.
- 9:56 am — Forward MC (s12_forward.py mc, 791 pairs × 4 variants × 8000 particles, 70 s). Check: the variant with CVM
  speed inside the LVL reproduces S11's a_d for 2–5 s 3C (0.548 vs 0.547 in kernel.log). Unit sum of the CVM categories
  is 1.000 for that variant (exact) and 0.94–0.95 with the slow LVL (the LVL's extra slowness is carried by the LVL tally).
- 10:05 am — RESULT F1: a_LVL (2–5 s, 3C) median 0.18–0.22 over kernels (pairs q10–q90 0.07–0.43). So a +2 % LVL change
  reads **+0.35–0.44 %** in S-anchored coda dv/v. Per window 1–2/2–3/3–4/4–6 s: 0.14–0.19/0.17–0.22/0.16–0.22/0.17–0.20.
  The LVL is large for the S-pinned kernel because its excess residence sits right at the source and the layer is slow.
  Ray theory: LVL-base SrS 1.4–2.3 s after S (median 1.8 s) with 2.95 s spent in the LVL → −59 ms for +2 %.
  Moho SrS (7 km) median 3.0 s after S. PrS arrives near S (−5.3..+1.7 s), SrP 3–12 s before S (P coda), PrP 0.05–1 s
  after P. Near-station epicentral distances are 30–55 km (median 42 km), S incidence at the source ~60°.
- 10:10 am — Consequence written into PLAN before measuring: round-7's within-episode bound (< 0.016 %, lag 0–1 vs ≥ 3)
  is ~25× below the coda target, so a local change cannot switch on inside the LFE-active days. The searchable contrasts
  are ETS vs between-ETS LFEs, pre-front LFEs, and the reflection.
- 10:12 am — PLAN.md written (rules for tasks 2–5) before any task-2..5 number. Code tests on 6/36 series and a 20-draw
  reflection test were run only to debug; their numbers are not read as results (test markers `*_test_DONE`; I removed
  the 20-draw REFLECT_DONE marker I had just made so the full 1000-draw run is not skipped).
- 10:13 am — DECISION E1 (classes): between-ETS = > 60 d after the previous and > 30 d before the next major front of the
  family. Why: Gosselin's change is significant in 15–60 d bins and decays; 60 d keeps the baseline in the decayed state.
  Columns betw_early (60–240 d) / betw_late (> 240 d) show any slow recovery.
- 10:13 am — DECISION E2 (measurement): superpose all majors into one class stack per series, then stretch (not per day):
  1-s lapse windows are too short for single-day stacks (S10: day σ 0.3–0.5 % in 2 s windows). Bootstrap over majors and
  inter-ETS intervals with the same draws for every series (joint), so co-located instruments are not counted as
  independent.
- 10:14 am — Queue launched (`scripts/s12_queue.sh`; block 15 added to `scripts/resume_after_reboot.sh`).
- 10:08 am — AMENDMENT R1 (post-hoc, after seeing SrP p 0.001 at the grid edge): the SrP/PrP scrambled control was
  invalid by construction (lags expressed after S, so P-locked energy was never permuted). Re-anchored on P, lags < 0.6 s
  after P excluded; rerun → SrP p 0.68. Logged in RESULTS §2.
- 10:10 am — Injection threshold for the SrS scan: detection (p < 0.01) from 0.08 × direct S; sharp LVL base predicts ~0.13.
- 10:11 am — Coordinator message: priority is something new; add A (propagating onset), B (recharge), C (source
  properties), D (episode size) after tasks 1–3. PLAN §6 written ~10:13 am before any A–D number.
- 10:12 am — Footprint MC (LVL change only within 15 / 25 km of PGC/SNB/VGZ), nominal coefficients reproduced exactly;
  pooling gains delta_foot15/25 columns (re-ran POOL for both regions; values of the other columns unchanged).
- 10:21 am — A–D done (queue2). WA B −0.33 ± 0.14 %/yr (2.3σ) — the only candidate.
- 10:23 am — AMENDMENT B1 (post-hoc checks, labelled): scrambled control (time-since-ETS shuffled within series, 500),
  leave-one-site-out, calendar trend, binned curve (`scripts/s12_recharge_check.py`). WA: scrambled p 0.004, LOSO −0.41..−0.26,
  trend −0.37; SVI −0.006 (p 0.47). The 3σ bootstrap rule is still not met → candidate.
- Housekeeping: I removed only my own code-test markers (`REFLECT_DONE` of the 20-draw test, `NOVEL_SVI_test_DONE`) and
  two `fwd_*_test` files from a 6-pair test, minutes after creating them. No result, figure or data of others touched.
- Not done: task 5 (splitting) — direct S does not cross the LVL and the reflection was not found; C time-of-day control.
- ~10:30 am — Coordinator: confirm the WA recharge candidate with existing data only. PLAN §7 written first (C1–C4 rules,
  forward lapse patterns: LVL flat, shallow grows, 25–45 km falls → "grows with lapse" would mean shallow). WA-SW counted
  before measuring: 0 between-ETS family-days at any P threshold → no replication possible.
- 10:34 am — `s12_b_measure.py` WA-NW: 860 series, 78 families (12 with fronts), 3374 between-ETS days (vs 738 in §6).
- 10:41 am — RESULT §7 WA: not confirmed. All families −0.055 ± 0.053 %/yr; the 39 families outside the §6 set +0.008 ±
  0.059; C4 (season) fails. The −0.33 belonged to the 11 front families. Tie-in to ETS size: none.
- 10:42 am — queue3 continues: SVI b-measure + era split, 1–3 Hz low-band build (raw day files, 7 workers, ~1–2 h), scan,
  injection, figures. Resume block 15b added.
- 10:52 am — 11-family check: not compact (p 0.73), slope extreme among random sets (p 0.009), WA ~99 km from PGC, deeper (42.8 vs 39.4 km) but other deep families show nothing (−0.03 ± 0.09) → scattered, most likely noise.
- 2026-10-14 2:50 pm PT — Cause of the 10-13 1–3 Hz scan crash (rc 1, queue.log 17:50 UTC): 10 of the 380 low-band pairs
  have NaN S−P (no sp in offsets_sp.csv) → `resid_logenv_P` took `max()` of an empty P window. The injection then failed
  only because `refl_traces_lb.npz` was never written. FIX: `resid_logenv_P` returns NaN for a non-finite or out-of-range
  P time (those pairs get no P-side residual; S-side phases unchanged; 2–8 Hz run unaffected, all sp finite there).
- 2:52–2:53 pm PT — Rerun scan + injection only (`scripts/s12_queue4.sh`, nohup, OMP 1; build reused; markers
  REFLECT_LB_DONE, REFLECT_INJECT_LB_DONE, QUEUE4_DONE). RESULT: no reflection (SrS LVL base scan p 0.19, nominal p 0.90);
  threshold 0.12 (p 0.02) – 0.18 (p 0.002) × direct S vs 0.08 at 2–8 Hz; sharp-base target ~0.13 → 1–3 Hz is LESS
  sensitive. PrP nominal p 0.001 is an artifact (only 8 pairs have a defined metric at the nominal lag; < 10 → NaN) —
  scan p 0.34 is the valid number. Deviation: nmax 100 used (PLAN said up to 150). RESULTS §7b.
- 2:55 pm PT — SVI era split written up (RESULTS §7c) from `tables/svi/r2_verdict.json`. POST HOC (labelled):
  `scripts/s12_era_check.py` adds a per-era scrambled control (300 draws: p 0.74 / 0.66) and the B009 instrument-veto
  refit (28 days dropped; 2003–2013 −0.0037 → −0.0046 %/yr). The interval bootstrap with 7 intervals per era moves σ by
  ~1.5× between runs (0.017 vs 0.027; 0.013 vs 0.020); both reported.
