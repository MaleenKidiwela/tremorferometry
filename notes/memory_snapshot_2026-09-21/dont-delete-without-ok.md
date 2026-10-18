---
name: dont-delete-without-ok
description: Never delete figures/outputs without EXPLICIT confirmation; a question about deletion is not approval
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 5cecf374-37e4-4bcc-aea9-79b2b6cdc569
---

Do NOT delete figures, results, or any outputs unless the user EXPLICITLY says to delete them.
A question like "did you delete the previous PGC figures?" is the user CHECKING, **not** an
instruction to delete — ask before acting.

**Why:** I deleted 8 old PGC figures after misreading such a question as approval; the user said
"Don't delete" and I had to restore them from git. The user wants old outputs KEPT for reference/
comparison even when re-doing a station on the clean pipeline.

**How to apply:** deleting raw WAVEFORMS after processing is fine when the user asks for it (re-
downloadable) — but figures/dv-v/results are deliverables: keep them unless told otherwise. When
in doubt, ask. (Git tracks `figures/` so committed ones are recoverable via `git checkout -- <f>`,
but don't rely on that.) See [[proactive-note-taking]].
