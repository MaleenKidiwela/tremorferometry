# lfe_picker

Array-level LFE detector for a margin-wide catalog of located families.

Read in this order: `PLAN.md` (design, stages, pre-registered bars, Merlin's amendments) →
`results/merlin_review_plan_2026-09-17.md` (the review that re-scoped the phase) →
`results/P0O_RESULTS.md` (oracle + referee freeze: verdict **DISCOVERY-ONLY**) →
`results/P1_RESULTS.md` (two learned single-station baselines) →
`results/P2_RESULTS.md` (learned CF + explicit Slab2 back-projection) → `NOTES.md` (dated ledger).

**Referee discipline** is enforced by `results/referee_protocol.json` and logged in
`results/referee_ledger.csv`: 2013 is the FINAL referee (one look, spent in P0-O and closed), 2011
gets exactly one evaluation per stage gate (P0-O, P1b and P2 have each spent one), 2010 and 2012 are
the development years where every choice was made, 2008/2009 are training-only.
30 hidden families and 3 held-out stations (`results/holdouts.json`) are excluded from every training
set and reported only as test strata.

Scripts are prefixed by stage: `o_*`, `o1..o3_*` = P0-O (read-only inputs to later stages);
`p_*`, `p1b_*`, `p2_*` = P1/P2. Inputs read-only: `../lfe_extension` (catalog v1, traces, frozen
chain), `../lfe_catalog_repair` (references, Lin's code), `../codex-LFE` (metadata).

Large artefacts: `results/mm/` is the memory-mapped ZNE day store (91 GB — rebuild with
`scripts/p_build_mm.py`), `results/p2/cf_*/` the per-station-day characteristic functions.
