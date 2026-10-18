---
name: reproduce-known-answer-first
description: User 2026-09-20 — when a published result exists on data we hold (Bostock's 130 families on POLARIS), a failed rediscovery means OUR CODE is wrong; debug against the known answer instead of declaring the method/archive non-reproducible or "closing the study"
metadata:
  type: feedback
---

User (2026-09-20): "do you realize that you are saying none of the work by Bostock is reproducible using their own stations? ... if it fails, Bostock families were proven that their methodology worked" and "your code is wrong, if you can't reproduce Bostock, so you are writing code or waiting".

**Why:** a positive control with a published answer is a debugging target, not a hypothesis test; pre-registration/STOP rules are for discovery claims, and applying them to an implementation phase produced a day of circular control re-specifications and a proposed "study closes" verdict that was logically invalid.
**How to apply:** (1) when known truth exists, trace where the chain loses known cases stage by stage and fix it; keep a small sealed holdout for the final blind check only; (2) never relay "not reproducible from the archive / method doesn't work" when the original authors demonstrably succeeded on the same data; (3) Merlin (and freezes) only at the final blind check during reproduction work; (4) be writing code, not waiting — launch the debugging agent immediately. See [[dont-over-conclude]], [[use-merlin-for-decisions]] (scope it to claims, not debugging).
