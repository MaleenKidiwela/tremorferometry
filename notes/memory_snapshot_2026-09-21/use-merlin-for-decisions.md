---
name: use-merlin-for-decisions
description: "Consult the Merlin agent BEFORE any method/scientific decision or when unsure — don't decide unilaterally"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: de9e1cec-3d07-493e-8e0a-d25339535f7a
---

The user wants the **merlin** agent (senior scientific advisor, Fable 5) consulted BEFORE making any
substantial method/scientific/statistical decision, or whenever a question or judgment call arises —
not decided unilaterally.

**Why:** In the 29-borehole dv/v run I made several method changes on my own (lowered the adaptive
candidate floor 0.5→0.2, added the 300-family densify cap, wrote a recap300 that would have discarded
reliable causality-certified families) — and the recap300 one was wrong; the user had to catch it
("isn't the goal dv/v from reliable LFEs?"). Merlin exists precisely to vet these before they cost time/
resources or corrupt a result.

**How to apply:** When choosing between methodologies, changing a pipeline parameter, when a result looks
suspicious, or facing any decision with a scientific/statistical dimension, spawn the `merlin` agent with
the specific claims/decisions to audit and act on its prescription. Merlin critiques/prescribes; it never
writes code. Related: [[project-mandate-2026-06-09]], the logic-supervisor advisor.
