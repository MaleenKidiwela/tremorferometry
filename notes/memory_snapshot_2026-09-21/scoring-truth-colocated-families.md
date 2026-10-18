---
name: scoring-truth-colocated-families
description: 2026-09-20 lesson — per-family one-to-one purity/recall against hp-only single-family labels is NOT a noise measure when Bostock families sit 0.1–1.6 km apart and 97.5 % of hp events are co-detected by a neighbour's template; use a ≤3 km cluster truth at low p_real before calling a template impure
metadata:
  type: feedback
---

P4-0b (lfe_picker/results/P4_0B_RESULTS.md) reported "purity collapse as templates grow" as the mechanism for 0/18 families assembled; Merlin's close-out showed two-thirds of the "impure" hits were the same family's frozen detections below the p_real 0.8 cut or a neighbour family 0.1–1.6 km away (97.5 % of hp events have n_other5 > 0). Re-scored with a 3-km cluster truth the relative purity of G44 went 0.70 → 0.85.

**Why:** an incomplete single-label truth makes purity fall exactly as a template gains sensitivity; the anti-correlation carries no information about noise admission. Same trap for the front estimator: per-family "first activation" is not independent across co-located families.
**How to apply:** before interpreting any per-family purity/recall on the Bostock catalog, (1) build the truth from all non-hidden families within ≤3 km of the target at low p_real, deduplicated at 2 s; (2) score the reference template against the same truth; (3) for front/timing estimators, de-duplicate events with n_other5 or cluster families at ≤3 km first. Also: a growth-threshold floor that the FA rule never reaches (T=8.0 gave 0.012/template-day vs 0.1 allowed) is a convention, not a measured limit — report it as such. See [[dont-over-conclude]], [[lfe-extension-study-2026-09-16]].
