---
name: orchestrate-with-opus-subagents
description: "User wants the main (Fable) session to ORCHESTRATE — delegate coding and simple tasks to Opus sub-agents, keep each study in its own folder, document everything"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 03499166-a677-42d8-9397-e2cf5e28621f
  modified: 2026-09-10T03:57:38.498Z
---

User directive (2026-09-10): "please use sub agents that are not fable (use opus) to do simple tasks, or coding";
"keep everything documented, also do it in a seperate folder"; "please act as the orchestrator".

**Why:** Fable time is expensive; coding is routine. The user wants the main session to plan, write specs, review and
integrate, while Opus agents write/run scripts. Studies must be self-contained and traceable.

**How to apply:** For any new study: create `<study>/` with README.md (inputs, conventions), PLAN.md (pre-registered
questions/quantities/controls/steps), NOTES.md (append-only dated ledger), scripts/, results/, figures/. Write a precise
spec (inputs with exact paths, definitions, outputs, doc requirements) and launch `Agent(subagent_type="general-purpose",
model="opus")`; run independent steps as parallel agents. Review agent output, append to NOTES.md, report concisely.
Never let agents modify files outside the study folder. Related: [[be-concise]], [[dont-over-conclude]].
