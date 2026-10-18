---
name: success-criteria
description: "What \"best subsurface monitoring system\" means for THIS project — judged by the geophysics literature, not generic/biology/cybersecurity \"monitoring\""
metadata: 
  node_type: memory
  type: project
  originSessionId: 18c98477-b5c6-4af0-9788-f240c207a136
---

User reframe (2026-06-09): "best subsurface monitoring system" is NOT a product or generic monitoring
pipeline ("this isnt biology or cyber security"). The comparison class is **geophysical time-lapse
velocity imaging**: ambient-noise dv/v (Brenguier et al.), coda-wave interferometry, repeating-earthquake
doublets, InSAR/GPS on the surface side. "Best" = beat that literature on concrete, falsifiable axes:

1. **DEPTH** — the novel axis. Nearly all dv/v monitoring images the shallow crust (surface src/rcv);
   ours is anchored on the deep megathrust (~30-40 km) via on-fault LFE sources. This is the real claim.
2. **SENSITIVITY** — smallest real δβ/β resolvable per cell (now ~0.1-0.15%; temporal-methods work aims lower).
3. **SPATIAL RESOLUTION** — fault cells resolved by ≥3 stations (post-2020 densification lever).
4. **TEMPORAL RESOLUTION & HONESTY** — detect a sudden change WHEN it happens, no smearing, no faked
   precursor (the differential step-injection benchmark, `scripts/bench_step_recovery.py`, is the referee).
   **TARGET CADENCE = per-DAY (aspirational; user 2026-06-10), per-month is only the current SNR-limited
   compromise, NOT the goal.** Judge metrics at the finest honestly-resolvable cadence; push toward daily.
   We are not yet stable at per-day — closing that gap is the objective. Implication: denser family coverage
   (more genuine families, even firing-day-only ones) + per-day estimators (pairwise/doublet, no 30-d smoothing)
   are what BUY daily resolution → strengthens the family-expansion case [[family-cwi-predictor]],
   [[dvv-temporal-estimators]]. Goal restated in `.claude/agents/logic-supervisor.md` (the Fable advisor).

**Why:** keeps us from over-hyping/over-engineering ops and anchors decisions to the science bar.
**How to apply:** justify each method/station/feature by which axis it advances vs the dv/v literature;
report results in those terms. See [[project-mandate-2026-06-09]], [[fault-tomography-next-phase]],
[[coda-window-2to4-correction]]. Don't borrow "monitoring system" buzz from other fields.
