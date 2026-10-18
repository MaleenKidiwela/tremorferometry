---
name: crash-recovery-context
description: "This project's sessions crash repeatedly; how to recover state fast"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: be83f309-c83d-4336-b8a2-99a689098c27
---

The user's environment crashes repeatedly ("crashed again"). On restart they ask to "review memory and files" to reconstruct what was in progress.

**Why:** much of the work is done in an interactive Python REPL/IPython (no `.ipynb` files in the repo), so the in-flight commands are lost on crash — only on-disk artifacts survive.

**How to apply:** to recover, check (1) memory for the active project, (2) `git status` for uncommitted/untracked files, (3) newest files by mtime (`find -mmin` / `ls -lat`) vs. last commit time to see what was built after the last commit, (4) whether scripts have produced their outputs yet. Keep [[nisqually-precursor-project]] updated as work progresses so recovery is instant. Working dir is `/home/jovyan/tremorferometry`.
