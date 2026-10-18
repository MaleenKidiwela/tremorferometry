---
name: agents-must-not-idle-wait
description: "User 2026-09-19 — Opus agents left waiting on long GPU jobs (wait loops, interim reports, drafting while waiting) burned half the weekly Opus quota; never leave an agent alive while a shell job runs"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: d6a923e6-f9f1-4ee1-8d37-955b6b2e0f35
  modified: 2026-09-19T22:13:54.272Z
---

User (2026-09-19): "it wasted 1/2 of my weekly usage" — after the P3D agent (750k+ tokens) sat in wait loops across ~9 hours of
GPU passes, sent interim status turns, and drafted write-ups while waiting; two other finished agents were still alive in
wait loops days later.

**Why:** every agent wake-up (wait-loop expiry, interim report, drafting) is a paid Opus call; a GPU pass needs none of them.
**How to apply:** (1) an agent LAUNCHES long jobs as a nohup shell queue and ENDS ITS TURN immediately; (2) the main session
watches for completion with a shell Monitor/until-loop (cheap), then resumes the agent ONCE for scoring/write-up; (3) brief
agents explicitly: "do not wait, poll, or draft while jobs run — exit after launching"; (4) TaskStop every agent as soon as its
report is in; (5) fewer, longer-scoped agents beat many; (6) no interim status turns ever. Reinforces
[[minimize-fable-delegate-opus]] — quota is the scarce resource on both models.

**Addendum 2026-09-20 — detach queues properly.** Two nohup queues launched by an Opus agent (Step 1b/1c, lfe_array_discovery) were SIGTERM-killed (-15) the moment the agent's session ended, because they were children of the agent's background Bash task. Agent D's queue, launched from a foreground Bash call, survived. Rule: launch every long queue as `setsid nohup bash queue.sh >> log 2>&1 < /dev/null &` from a FOREGROUND Bash call (never inside a run_in_background task), and have the orchestrator verify with `pgrep -af queue.sh` after the agent exits. Queue scripts must be resumable per step (guard each stage by its output file) so a relaunch does not redo hours of work — 1b had no guards and re-ran from scratch.
