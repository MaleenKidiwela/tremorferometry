---
name: rank-dont-discard
description: User 2026-09-23 — never filter candidate families with cuts; keep every row, turn every criterion into a column, and show real-vs-null counts at every level. The 85%-of-known-families test is a VETO on proposed cuts, never a licence to apply one, because undiscovered families need not resemble the known population
metadata:
  type: feedback
---

User (2026-09-23): "this is a dumb rule, you are discarding them, what if there are undiscovered families" — objecting to
"anything we add has to pass 85 % of the 58 rediscovered families or it is deleted".

**Why they are right:** calibrating acceptance criteria against a known population (Bostock's northern families, one array,
one band, one size distribution) can only reject criteria that are too strict FOR THAT POPULATION. It gives no protection to a
real family that is smaller, deeper, noisier, or seen by fewer stations — which is exactly what an undiscovered family may be.
Passing the 85 % test is not evidence a criterion is valid.

**How to apply:**
1. Do not filter. Every grown family keeps its row in the output permanently.
2. Every criterion becomes a COLUMN (measured value), never a gate.
3. Run the identical chain on a site-shift/whole-chain null and produce the SAME table for it.
4. Evidence = the real-vs-null count as a function of each criterion (a separation curve), plus a ranking — not a single
   pre-chosen threshold.
5. The ≥85 %-of-known-families test survives ONLY as a veto: a proposed cut that fails it may not be used at all. It never
   licenses a cut.
6. Only the acceptance steps the papers actually perform (dedupe by shared detections or waveform correlation; location) may
   reduce a candidate list, and even they are reported with counts before and after.
See [[justify-criteria-on-known-families]] (superseded in part by this), [[reproduce-known-answer-first]].
