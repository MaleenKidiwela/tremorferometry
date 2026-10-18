# Merlin audit of E5 step 5 (Sept-2012 southern-VI replication) — 2026-09-11 ~13:15 UTC (Fable 5.1; read-only; verbatim summary)

**Verdict.** The pre-registered replication FAILED and the agent reported it honestly; nothing was tuned. The agent's diagnosis (quiet-day
excess = noise coherent across the PGC site, invisible to the per-station rolled control) is right in direction — strengthened by two
independent tests — but its mechanism proof is incomplete (the SITE-rolled control on the same network was never run; the one-per-site run
is confounded because it also halves recall), and the proposed fix (site-merge) must be judged on recall at matched false-alarm rate,
because the cluster carries real signal. The decision that matters most: STOP treating the sensitivity gap as a template-N problem. Own-
family recall is flat in template event count above ~50 build events and is set by distance to the single good site; no on-disk lever
reaches Bostock-class recall on this geometry. Fix the FA model (one diagnostic pass, one design pass), re-scope the Thomas science to
what ~25-35% family-dependent completeness supports, and gate the post-2019 download on a validated chance control, not on recall.

## Stands / falls
Stands: download, split, offsets, from-scratch templates (q-control 0.78), per-station registration at the model, own-family recall 0.234
(chance ~0), FAIL on the FA bars, daily MAD uncontaminated by signal (ETS/quiet MAD ratio p10/50/90 0.93/1.02/1.13), main-front speed by
first-burst, re-rupture interval distribution.
Falls / relabel: "catalog-definition MET" cells (chance-inflated); "one-per-site settles it"; "merge every co-located cluster" as a fleet
rule (one site only); "FA wrong, NOT sensitivity short" (both apply); secondary-front count 13 vs 20 as a completeness statement (metric
not density-monotonic: full 20, random 35% -> 26, Mc1.8 -> 14); uniformly-decimated reference; the 14-Sept front "recovered" (gate-edge
coin flip, rho 0.5996 vs 0.6).

## Flaws
F1 MAJOR Mechanism inferred not measured: e4_template_test.py:522 draws shift_s per station; :559 rolls each station separately -> the
   control destroys within-site AND between-site coherence. One-per-site removes B009/B010/B011/B005-B007 entirely; on doy 243-244 own
   recall falls 22/123 -> 10/123 (detections 668 -> 124). Cluster stations carry real signal (registration 1.66-1.87 MAD).
F2 MAJOR Noise-vs-LFE argued by absence. Positive tests (Merlin ran on e5r_detections_e5r_snr.parquet, L=8, >=8 MAD): 44.5% of quiet-day
   detections have a DIFFERENT family's detection within +-2 s (54% within +-5 s), 6-20x Poisson; those pairs' inter-family distance median
   10.8 km, p90 42 km, 21% beyond 30 km vs ETS-day own-matched 2.7 km / 2% -> a single-site transient maps onto every family at similar
   hypocentral distance from PGC (annulus). Night/day contrast 30-60x (04-11 UTC 321-1023 vs 13-35 in 14-22 UTC) vs ~5x for matched
   ETS-day detections. Excess correlates with PGC-cluster weight (rho 0.52) not LRIV-cluster (rho -0.13). Positive stack test still needed.
F3 MAJOR Scoring: e5r_score.py:207-216 sets met_cat from raw cat recall while chance_cat = 0.064; chance-corrected: overall 0.30 (>=0.25),
   M1.8-2.0 0.33 < 0.35 (flips to NOT MET), M>=2.0 0.40. The "rolled <=5%" bar uses the per-station roll, which understates chance ~4.6x on
   quiet days; honest ETS-day floor fraction ~ quiet rate/ETS rate ~ 4.47/12.4 ~ 36% at 8 MAD. Own-family recall (chance 0.0002) is the only
   chance-clean sensitivity number — primary everywhere.
F4 MAJOR Generalisation: only the PGC site (PGC + B009/B010/B011, 5-22 m, coastal) shows excess; the LRIV cluster (293-303 m, inland)
   none. Rule that generalises: always use a SITE-rolled chance control (site = stations within ~1 km); merge only where it shows
   within-site excess.
F5 MAJOR Framing: sensitivity IS short, for geometry not N. Event-weighted own recall by build-event count: 0.10 (<=50), 0.23, 0.21, 0.27,
   0.25 (200-300) — flat above 50; by distance from PGC 0.27-0.36 within 30 km, 0.17-0.18 beyond 40 km; by # of >=7-dB stations within
   40 km: 0.17/0.19/0.28/0.32 for 0/1/2/3. Family median magnitude correlates with distance from PGC (rho 0.69) -> recall looks
   magnitude-flat. PGC sub-stack S-SNR saturates at ~12 dB by N~100; far boreholes still rising at 200. N binds only for 30 families with
   <50 build events (6% of test events). "Recall of his catalog" is a lower bound on relative sensitivity (his marginal events selected on
   HIS noise) but 0.23 is far from "equivalent".
F6 MAJOR Secondary-front COUNT not density-monotonic (rho gate on all events in a 3-h window; thinning raises rho variance). "13 vs 20"
   uninterpretable; metric mis-specification.
F7 Uniform random decimation is the wrong reference: recall is family-dependent (p10/50/90 = 0/0.17/0.42; 22 templated families <=0.05
   hold 1,951 of 13,381 test events) and the catalog carries a ~40% ETS-day floor -> per-family Bernoulli thinning + uniform-time floor,
   >=20 seeds. "Re-rupture 0.25 within 0.04 of 0.29" untested until then.
F8 MINOR Cross-family declustering by max network MAD (e5r_score.py:120) biases toward high-floor PGC-dominated families; template-less
   families get 0.136 "recall" by aliasing. Assign by raw weighted cc among families sharing >=50% of stations.
F9 MINOR Far-station S picks at +1.7..+3.3 s sit at the +-4 s search edge (e5r_templates.py:91) — some clipped. EQ veto radius unverified;
   nocturnal pattern rules out EQ dominance regardless.

## Answers
A Noise, high confidence (F2 + PNSN zero on quiet days while our catalog runs 100-1,000/day). Missing: positive stack test + site-rolled
  control. B Honest precision: 8 MAD ~20-25% overall / 40-55% ETS days (floor ~350/day, non-stationary: 09-26/27 1,861/1,551 with zero
  Bostock); 9 MAD ~45-55% / 65-75%. Use >=9 MAD for anything quantitative NOW. Future chance control: roll SITES (common shift per site),
  keep station-rolled as lower bound. C Levers at matched FA: (1) site FA fix decides 8 vs 9 MAD (+40-60% relative); (2) path-corrected S /
  wider pick search at far boreholes ~+5%; (3) drop <5 dB stations +-5%; (4) 2010+2012 templates for 7 persistent stations +0-10%;
  (5) template length nil at matched FA; (6) more stations excluded/irrelevant. CEILING of all on-disk levers ~0.30-0.35 own recall.
  Bostock-class (>=0.5) needs sites within ~30 km of each family — geometry, not software. D 14-Sept front = coin flip; report rho +
  speed with CI in Thomas's windows. >=52% completeness not reachable on 15 stations. E F3, F6-F8; 2-h split sound for this test.

## Prescription (no downloads)
1 SITE-rolled chance control on the unchanged network (one GPU pass ~7 min): stations within 1 km share one shift ({PGC,B009,B010,B011},
  {LRIV,B005,B006,B007}); write both controls; store per-site contribution to the sum at each detection. Mechanism confirmed if quiet-day
  site-rolled rate within 1.3x of data (~4.5) and PGC-site share at quiet detections >> its weight share. Failure: site-rolled ~ station-rolled
  -> excess is BETWEEN-site coherent; do not merge; 1b decides.
1b POSITIVE stack test: re-run the 5 worst quiet-FA families (bkeys 103, 32, 75, 128 + next) on quiet days with VGZ and SNB EXCLUDED; stack
  VGZ/SNB 3-C at the >=8-MAD detection times with the family moveout (no refinement); S-SNR vs same-N random stack. Noise: <= null + 1 dB
  at N>=300. >= null + 3 dB with the arrival at predicted S -> real inter-ETS LFEs (quiet-day bar void; a finding).
2 SITE-MERGED network, two pre-registered variants, one GPU pass each: within-site w-weighted average of member cc series; site weight
  (A) max member w, (B) sqrt(sum w²). Chance control = site-rolled. Operating T = smallest of {7.5, 8, 8.5, 9} meeting quiet <=1/template/day
  AND site-rolled/data <=1.3 AND ETS-day site-rolled fraction <=5%. Sensitivity bars unchanged (own >=0.25 overall, >=0.20 at M1.6-1.8; the
  M1.8-2.0 >=0.35 / M>=2.0 >=0.60 bars stay on record but encode a slope the data refute). Adopt the variant with higher own recall at its
  compliant T; "merge did not destroy sensitivity" = own >=0.20. Failure: neither reaches quiet <=1 below 9 MAD -> cap PGC site <=25% of the
  sum or add a two-site coincidence gate.
3 Recall levers at the step-2 operating point (real if >=+0.02 own recall): (a) pick search +-8 s + per-station median residual correction
  for LRIV/B001-B007; (b) network restricted to stations >=5 dB; (c) 2010+2012 combined templates for the 7 persistent stations (align 2010
  stack to 2012 by cc; never touch the 2012 test half).
4 Fix scoring/reference (CPU): chance-corrected catalog recall; per-family Bernoulli-thinned Bostock reference + uniform floor, 20 seeds;
  replace secondary-front count by speed + rho + CI in Thomas's named windows and a density-robust statistic on completeness-matched counts.
5 Decision rules: (a) "Bostock-equivalent" only if own recall >=0.5 at compliant FA — NOT reachable on this geometry; honest product = "a
  detector with a validated chance model at 25-35% family-dependent completeness" — say so. (b) Post-2019 download (VI 2024-08 first;
  Puget needs T3 discovery first) only when step 1 confirms the mechanism, step 2 passes its FA bars on BOTH the quiet and post-SSE
  (09-26..10-05) windows, and step 4's re-scoped metrics are pre-registered. (c) RE-SCOPE NOW: reachable at ~30% = main-front speed
  (first-burst), re-rupture interval distribution, first-burst migration per family, per-family rates with injected-template completeness
  weights; NOT reachable = secondary-front counts, re-rupture family fraction, tides, rate decay. Cross-cycle comparisons only on
  completeness-matched quantities.

## What would change Merlin's mind
Step 1 site-rolled ~ station-rolled AND 1b S-SNR >= null + 3 dB at excluded stations -> quiet-day detections are real inter-ETS LFEs at
hundreds/day, FA bars void, recall side a near-pass. Step 3c >= +0.05 -> N was binding; 2011/2013 ETS-peak downloads (~5-7 GB each) become
the lever. Neither likely (<10%).
