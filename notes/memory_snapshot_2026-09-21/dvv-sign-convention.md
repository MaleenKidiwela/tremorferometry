---
name: dvv-sign-convention
description: "The stored dv/v column is eps = -(dv/v); map fixed 2026-08-20, inversion code fixed 2026-08-22 (results on disk still old)"
metadata: 
  node_type: memory
  type: project
  originSessionId: b2cdc227-d9f7-47e6-b5ae-77db95eb27af
  modified: 2026-08-22T23:27:41.791Z
---

The `dvv` column in every `data/daily_dvv_*.csv` is the STRETCH eps = dt/t, which equals **-(dv/v)**.
A velocity INCREASE returns a NEGATIVE number. Proven by forward test (Merlin, 2026-08-20): synthesise a
physical dv/v of +1.00%, pass it through the repo's own `stretch_dvv`, get -0.99% back.

Root cause: `src/tremorferometry/dvv.py` had a "Convention:" note contradicting its own forward model six
lines above; `dvv_roll30cal.py:22` documented "dv/v = -eps" but never applied it. Survived three audits.

**Fixed 2026-08-20 at ONE boundary only** — `build_borehole_dvv_map.py` has `DVV_SIGN = -1.0` applied at the
display pivots. Estimator and CSVs deliberately left as eps (bit-identical, provenance). **Never negate in
two places — a second flip re-inverts it.**

**Inversion CODE fixed 2026-08-22** (independently re-verified: forward test of `stretch_dvv` on a
physically-constructed speedup returned -0.0099 for +1% dv/v; `G.npz` inspected — 100% negative entries,
rows sum to -1 = stores -K). `invert_dvv_4d.py` no longer flips G (`eps = G_stored·m` as stored → m = true
δβ/β); `SITES`, `idx_data`, `closure_mean` negated at save so every saved/printed product is in TRUE
velocity-change sign. **`inversion_4d.npz` + figures on disk are still the OLD sign-inverted run — not
re-run yet (user said wait).** Old note for reference: the false premise (column is dv/v) is stated in
`invert_epoch.py:6`. Unaffected: all magnitudes,
nulls, checkerboards, resolution, the stability bound (sign-symmetric). Affected: ledger R6's pre-registered
"sign CORRECT" — a bookkeeping question, needs re-scoring. The round-2 synthetic injection CANNOT detect
this (plants in model space, forward-models through the same G — self-consistent by construction).

**FULL-CHAIN AUDIT 2026-08-22 (evening, daily pipeline): ALL HOPS CONSISTENT.** Verified at code level:
`stretch_dvv` samples cur at t(1+eps) vs ref at t → returned dvv = eps = dt/t = −(dv/v) (root ✓);
`dvv_roll30cal.py` writes eps unchanged ✓; assembler ×100 only ✓; `build_G_captured.py` stores G = −K
(rows sum −1) ✓; `invert_dvv_4d_daily.py` uses G as stored (no flip), so MODEL = TRUE δβ/β, s & c in
eps units, negated ONCE at save (SITES = −S, daily_offsets c_true_dvv = −c) ✓; transmission T,
deconvolution A/agg, and all movies/figures operate on true-sign products, no further negation ✓.
Empirical: +1% injection → T median +0.69 (0.05% of entries negative = near-zero-illumination ringing,
|T|~0.01, floored at 0.1 before any division — benign). One DOC fix: methods artifact eq (6) symbol
K→G_c. The annual `inversion_4d.npz` on disk REMAINS old-sign.

**User's methodological caution (2026-08-20):** do NOT use "ETS should cause a slowdown" as an acceptance
test for the sign. Tremor is a proxy for slow slip, not slip, and whether slip produces a measurable velocity
change is the hypothesis under test — validating tooling against it bakes the hypothesis in. The forward test
settles the convention on its own. Planned instead: model fault slip from GNSS and compare per-location slip
to dv/v (parked as of 2026-08-20).
