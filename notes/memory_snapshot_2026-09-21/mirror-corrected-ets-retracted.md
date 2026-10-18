---
name: mirror-corrected-ets-retracted
description: "USER RETRACTED the mirror-corrected +0.043% ETS value (2026-08-23): do not use it as a signal-class anchor"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: eeb752df-b1e8-4172-9ab5-f9abf3317b7f
  modified: 2026-08-23T04:41:11.196Z
---

2026-08-23, user: "forget the mirror corrected value, its not good." The mirror-corrected ETS-onset
residual (+0.043%, p=0.019, fig mirror_corrected_dvv_ets.png — cited in [[dvv-sign-convention]]-era
notes and the FINAL_PIPELINE memory as THE validated signal class) is NO LONGER to be used as a
reference signal size.

**Why:** user judged the mirror-corrected value unreliable (reason not elaborated; the supervisor had
already paused fleet-wide mirror CLEANING on 2026-07-22 — mirror-as-control/certification remains OK).

**How to apply:** (1) never quote 0.043% (or the 0.05% rounded bar) as "the validated signal class";
(2) the daily-product claim bar (transmission_deconv mask at 0.05% agg sigma) is a PLACEHOLDER
convention pending the empirical per-claim nulls (Merlin round-3 item B: circular-shift nulls through
the full pipeline supersede any fixed bar); tile noise RANKING stands, the pass/fail line is soft;
(3) mirror-as-CONTROL (causality certification, coda-vs-mirror ratio) is unaffected; (4) before any
mirror-v2 noise-reduction run, re-confirm with the user — their confidence in mirror CORRECTION is low.
