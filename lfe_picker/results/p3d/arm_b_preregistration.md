# Step 7 — Arm B (seed-and-grow) PRE-REGISTRATION. WRITE-ONLY; nothing below has been run.

Written 2026-09-19 because Step 3 did **not** license Step 6 on the V1 network
(predicted δ_full = 0.57 σ against the 1.5 σ bar; `results/p3d/step3_verdict.json`).
Merlin's prior for this arm is ~25 %, against ≤10 % for the pairwise chain.
This file is the protocol. It is frozen by being written before any of it is executed; its sha256 is
recorded in `results/p3d/freeze.json` under `arm_b`.

## Why a different arm
The pairwise chain asks every pair of windows to link on its own. Its statistic is the site sum, and
Step 3 measured that the site sum is dominated by ONE site (PGCq alone δ = 2.44 σ; every added site
dilutes it). Bostock's own construction does not ask that: it starts from a LOUD window, uses it as
a template, and lets the matched filter over the whole record supply the N that a single pair can
never have. The quantity that has to be large is then not the per-pair tail probability but the
number of events a single loud seed can pull in — which grows with the record length, not with δ⁴.

## Protocol (pre-registered, in order)

1. **Seeds.** On PNSN tremor hours inside the family box, over the pooled TEST episodes, compute the
   per-site incoherent S-energy in 12-s windows at 6-s hop on the dense ≤60 km site set
   (V1: PGCq, VGZ, SYMB, BUCKS; V2 adds LHLYB, JRDRB, ALBHB, DAOB if they deliver). Take the top
   **1 %** of windows by the *median* per-site z (median, not sum, so a single loud site cannot
   nominate a seed on its own — the Step-3 dilution result says the sum would be PGCq alone).
   Withheld sites (B926q, GOBB, LOPEZ, TURTL) are NOT in the seed statistic.
2. **Matched filter.** Each seed's 12-s multi-station window becomes a template; run it across the
   pooled test episodes with the frozen chain and the frozen two-site coincidence gate. Hits above
   the **1e-4** quantile of the site-sum null (index × 1,800 s shift, 5 seeds) are kept.
3. **Stack and iterate twice.** `align_stack` (2 iterations, ±1.5 s), re-pick S on the aligned
   stack, rebuild the 12-s template from S−2 s, re-run. Two iterations, then stop — iterating to
   convergence is a forking path and is not permitted.
4. **De-duplication.** Cluster third-iteration templates by multi-station stack cc ≥ each station's
   cross-family p95 (`results/gates/g_crossfam_cc.csv`); one candidate per cluster; a candidate
   whose PGC stack cc exceeds 0.524 against any of the 130 accepted families' 2012 PGC templates is
   an ALIAS of that family, not a discovery, and is reported as such.
5. **Acceptance.** A cluster is a candidate family only if
   (a) its third-iteration stack passes **Gate 1 at ≥2 WITHHELD sites** (`g_gates.evaluate`,
       `cfg.withheld = B926q, GOBB, LOPEZ, TURTL`), and
   (b) its member count exceeds **N-shift max + 3** with N-shift mean ≤10 % of the true count
       (the P3 § 6.3 claim rule), and
   (c) Gates 0 and 2 pass and Gate 3 passes or is untestable, and
   (d) it re-fires in an episode it was not discovered in (cross-episode reproducibility).
6. **Kill test FIRST.** Before any of the above is believed, run steps 1–5 on the **KEPT** families:
   seeds restricted to windows containing a known hp event must rediscover those families at
   ≥0.5 × the truth template's own per-family recall (the relative clause of `q_rediscover.py`).
   If the construction cannot rediscover families it is *given*, it is not run on the unknown.
7. **Scoring.** The 10 isolated hidden families are the truth, with the SAME caveat Merlin recorded:
   only 12, 29, 62 (V1) gain a genuinely new near site and only ≤6 of the 10 have ≥20 hp truth
   events in 2023–26, so "≥5/10" is unreachable by construction. The pre-registered bar is on the
   SCOREABLE isolated stratum (hidden families with ≥20 hp truth events in the test episodes):
   **PROCEED if ≥2/3 are rediscovered**; unscoreable families are reported by name with their Δsites.

## What would make this worth running
Step 3 says the single-site statistic at PGCq is already 2.44 σ, and SYMB alone 1.08 σ. Arm B does
not need a large per-pair δ; it needs one loud window per family per episode and enough record to
matched-filter. The measured hp populations (S18 alone: 2,192 hp events over 115 families) say that
condition is met for the large families and fails for exactly the isolated hidden ones that the
score depends on. **The honest expectation is therefore: Arm B will rediscover LARGE families and
will not settle the isolated-hidden score.** That is worth knowing, and it is a different claim from
"discovery works"; it must be written that way if it is run.
