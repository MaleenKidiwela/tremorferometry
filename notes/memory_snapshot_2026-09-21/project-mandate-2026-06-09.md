---
name: project-mandate-2026-06-09
description: "Standing autonomous mandate + hard rules — goal is the best subsurface monitoring system; flag don't delete; one station's waveforms at a time; watch memory"
metadata: 
  node_type: memory
  type: project
  originSessionId: 18c98477-b5c6-4af0-9788-f240c207a136
---

Granted 2026-06-09: standing permission to autonomously do whatever advances the MAIN GOAL — a **crisp 4D
velocity field of the Cascadia interface; "the best subsurface monitoring system there is."** Long-term
stations are done; post-2020-only stations are now in scope for spatially dense recent epochs
(`data/candidate_stations_post2020.csv`, 248 candidates ≥5yr within 60km of the dense band).

**Hard rules (user, same day):**
1. **NEVER delete/overwrite existing products** — bad/gappy families get FLAGGED (`data/family_quality_flags.csv`,
   filter on `quality` col), never removed. Future families: free rein.
2. **Few stations at a time; delete each station's waveforms before the next download** (storage discipline).
3. **Monitor pod memory** (anon, not memory.current) so it doesn't crash — `scripts/mem_watchdog.sh` runs in
   background (WARN ≥90G, kill largest python ≥118G, logs to logs/mem_watchdog.log); restart it after a crash.

**Why:** user wants sustained autonomous progress between check-ins; crashes and storage bloat have burned us.
**How to apply:** before any deletion ask "is this a product or a re-downloadable input?"; keep the watchdog
alive; when batching stations, serialize downloads. See [[margin-station-status]], [[dont-delete-without-ok]].
