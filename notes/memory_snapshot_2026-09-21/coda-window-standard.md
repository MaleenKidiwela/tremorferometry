---
name: coda-window-standard
description: "Coda-window choice for dv/v — 1-4s is the margin standard; old PGC 0-2s was wrong (ballistic, not coda); 1-3 vs 1-4 standardization is open"
metadata: 
  node_type: memory
  type: reference
  originSessionId: 400c26f9-17ff-46fc-914c-9c6a7671906d
---

The coda window must start **after** the direct arrival so it samples multiply-scattered later energy (where dv/v is stable and coherent), not the ballistic/early pulse.

- **Margin standard = 1–4 s** (HDW, CPW, B018, new PGC default to it).
- **GNW original** and the **old PGC** correct-window plot were **1–3 s** — also a valid coda window; 1–3 s and 1–4 s give essentially the same structure (annual cycle + station-specific multi-year lows), which is good evidence the dv/v is real coda behavior, not a window artifact.
- **0–2 s is WRONG** — it starts at the direct arrival and samples the ballistic/early part, so it is not proper coda interferometry. The standalone old `figures/smoke_dvv_PGC_79.png` (0–2 s) is spiky/unstable for this reason; the 1–3 s top panel of `smoke_dvv_PGC_79_window_compare.png` is the correct old PGC.

**PGC validation (2026-05-30):** new clean PGC re-run at 1–3 s (~96.7k measurements) reproduces the old 1–3 s signal — both show a 2012–2014 low + recovery — confirming the clean pipeline (resample_poly + coverage families + per-era for the 2017 BHZ→HHZ change) recovers the real signal, not an artifact.

**OPEN DECISION:** the margin is mixed — GNW + old PGC on 1–3 s, HDW/CPW/B018/new PGC on 1–4 s. For a single cross-station comparison the whole margin should be re-cut to one window (1–4 s is the recent default). Not yet done. Related: [[margin-station-status]], [[nisqually-precursor-project]].
