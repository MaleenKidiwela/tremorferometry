# RESUME — Royer & Bostock (2014) Washington reproduction (written 2026-09-26 ~16:40Z)

Read this first after a reset. Detailed log: `notes/2026-09-25_Notes.md` (sections A–T).
Agent log: `lfe_array_discovery/results/ab_royer/AB_ROYER_LOG.md`. Plan: `notes/MARGIN_CATALOG_PLAN_2026-09-24.md`.

## User rules for this work (also in memory)
- Follow Royer's method and GET FAMILIES; no side tests/extras unless asked. New locatable families count too.
- Finish NW first, then SW, then NVI; extension to 2026 later (Phase 3).
- Say "scrambled control", never "null". Don't gate on ETS/tremor windows (labels only).
- Agents: launch-and-exit, TaskStop on hand-back. No git commits/pushes. Use the GPU.

## Where it stands
- Chain: `lfe_array_discovery/scripts/ab_*.py`; queue `scripts/ab_royer_queue.sh` (v4: NW grow→locate, SW grow→locate,
  AB_ROYER_WA_DONE, then NVI). Main queue is **PAUSED** (marker results/AB_ROYER_FAILED.v5_paused_for_NW_first).
- Autocorrelation reproduces Royer: NW 4603 (vs 4915), SW 3217 (vs 3306), NVI 2039 (vs 4267, open).
- NW growth (8 iterations, then continued to 14): 169 families, template SNR ~3.9. Located 45 (Hyp2000, GSC 1-D, Vp/Vs 1.73, ratio 1.5).
- **Eye review (2–8 Hz restack, results/ab_royer/NW/hypo_28/review/eye_class.csv): 33 real LFE families, 136 noise.**
  The 33: all located, depth median 35 km (22 inside 31–46), 21/122 Royer templates matched within 5 km (p=0.017).
  Noise mostly from station W060 (in 126/136 noise families' growth stations) + FACB; ringing at the B013/W030/B943 site; GNW.
- Truth: Royer template locations digitized from the vector figures → `results/ab_royer/truth/` (scorer scripts/ab_truth_score.py).
- Slab-fixed depth (data/slab_models/): no gain — the limit is picks. SW slab follow-up waits for SW hypo (pid 2245051).

## Running now
1. **NW station-QC growth rerun** (W060 + FACB removed from the pool): `scripts/nw_qc_queue.sh`, started 16:19Z,
   outputs NW/grow_qc + NW/hypo_qc, markers results/NW_QC_{DONE,FAILED,HEARTBEAT}. ~20 s/day, ~1 h per iteration, up to 14 iterations.
2. Agent writing `scripts/ab_grow_fast.py` (CPU prefetch overlapping the GPU, equivalence-tested) for SW and later runs.

## Next
- When NW_QC_DONE: 2–8 Hz restack (scripts/ab_restack28.py) → eye review → count real families → locate → score vs Royer.
- Eye-pick the real families to tighten locations (erh median 5.4 km now).
- Then SW with the same station QC + ab_grow_fast; then NVI; then extension to 2026 (78 velocity stations in the box in 2025, 52 of them started 2012 or later).
- Open: all-days pass code exists (ab_alldays*.py) but is stopped; the station-extension block at the end of ab_royer_queue.sh must not run unless asked.
