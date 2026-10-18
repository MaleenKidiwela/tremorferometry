---
name: minimize-fable-delegate-opus
description: "User 2026-09-16 — Fable (main session) usage is near its limit; do coordination only, delegate all coding/analysis/monitoring to Opus agents, keep replies short"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: d6a923e6-f9f1-4ee1-8d37-955b6b2e0f35
  modified: 2026-09-16T23:15:33.860Z
---

User (2026-09-16, late in the lfe_extension session): "Please use opus agents and you as orchestrator because we are approaching
limit on fable."

**Why:** Fable quota is the scarce resource; Opus agent tokens are not. Long Fable analysis turns (like the dv/v-vs-Bostock identity
test I ran inline) burn the budget.
**How to apply:** main session = plan, brief, arbitrate, log to NOTES, answer the user briefly. Any computation, file inspection
beyond a few lines, monitoring loops, or write-ups → an Opus agent (subagent_type general-purpose, model opus). Prefer one
long-running agent that polls and assembles (e.g. Track D) over repeated Fable check-ins. Reinforces [[orchestrate-with-opus-subagents]]
and [[be-concise]].

**Addendum 2026-09-20:** resuming an Opus agent re-reads its whole transcript — agent E reached ~500k tokens and each resume cost that again. Rule: once an agent's transcript passes ~200k tokens, spawn a FRESH agent for the next task with file pointers (results files, freeze json, scripts) instead of resuming; reserve resumes for short follow-ups on small-context agents.
