---
name: lfe-catalog-goal-2026-09-10
description: CURRENT GOAL (user, 2026-09-10) — Thomas et al. 2025 (Sci. Adv. adz2832) at PNW scale over the full record — a reliable LOCATABLE LFE catalog (DL picker -> association -> location) to resolve fronts/secondary fronts/re-rupture/healing/tides for EVERY SSE
metadata:
  type: project
---

**Since 2026-09-10 the project goal is: do Thomas et al. 2025 (Sci. Adv. adz2832, 'Rapid fault healing from cementation controls the
dynamics of deep slow slip and tremor') for ALL of PNW over the ENTIRE record.** They used the Bostock family-located LFE catalog + PNSN
tremor for ONE SSE (Sept 2012, S. Vancouver Island): main front ~10 km/day, secondary fronts 10-40x faster (2-6 h) re-rupturing segments
within days, per-family rate decay, tidal sensitivity growing through the SSE, secondary-front stress drops ~8 kPa. We want that for every
SSE, every segment, 2005/2009-2026. The reliable, locatable LFE catalog is the MEANS. dv/v is closed ([[sse-model-vs-dvv-2026-09-03]]);
`tremor_slip_geometry/` is the PNSN-tremor precursor of the same dynamics questions, not a separate goal.
Catalog must deliver: event times + locations (or family locations + event times where sparse), family assignment, an amplitude->moment
module (Bostock-2015-style), documented time-varying completeness, Slab2 geometry per source. Goal statement is at the top of
`lfe_catalog_repair/PLAN.md`.

**Why:** the 1.5 TB matched-filter detections are template-shaped noise with a small true fraction (~108 false hits/template/day at
cc>=0.80; individual LFEs are ~1.8 MAD on a single station) — NOT an event catalog ([[mf-detection-catalogs-exist]]). The user pushed
back on the template-matching limit: a deep-learning picker using the multi-second 3-C signature is not bound by it (Lin 2024 catalogued
>1 M LFEs on Vancouver Island that way). Depth/location is the point — "locatable".

**How to apply:** workspace `lfe_catalog_repair/` (README/PLAN/NOTES/RESULTS). Method = replicate the METHOD not the weights: SeisBench
LFEDetect (`cascadia` vs joint `jcms`) -> PyOcto (N>=3 on VI, N>=4 elsewhere) -> NonLinLoc + casc1.6 (S-weighted; family/cluster centroids
only), fine-tune before trusting anything south of Puget. Stages: E1 resources (done), E2 pilot on VI ground truth (Lin/Bostock 2010,
VI 2017, off-truth WISH/YACT/B943), F1 fleet inventory + association regions + cost, then E2->fleet decision (streamed download ->
pick -> delete, ~20 TB, weeks). Pass bar: family hypocentres to ~2-3 km with fleet geometry; Lin (residual<0.5, N>=4) and Bostock as
separate calibration targets. Orchestration per [[orchestrate-with-opus-subagents]]. Expect Oregon/N-California to yield an activity
time series + template matching against Ducellier/Shelly families, not hypocentres, unless fine-tuning works.

**STATUS 2026-09-11 (E4, lfe_catalog_repair/RESULTS.md):** ML picker (LFEDetect) -> PyOcto = 1-3% of Bostock's events, chance-limited
-> FAMILY SEED only; denoising HURTS; fine-tuning unresolved. **TEMPLATE ROUTE OPEN on our station class**: Bostock's own events register at
4.4-6.2 MAD (network median, 6 CN ring stations) once the SeisBench stacks are read as Z,N,E (stored Z12, 1=N — a sorted() bug had rotated
them). All 14 POLARIS inner stations closed 2004-07 -> Bostock's 2010/2012 detections were made on ~our 7-station ring; implied threshold
3.6 MAD. Threshold binds, not geometry: at 4 MAD Puget/N-Oregon reach ~0.9 of Bostock density from 2023. Required design: SNR-weighted
network stacking, 8-12 s templates, network-sum threshold per template per day. Thomas harness (`scripts/e4_thomas_harness.py`,
catalog-agnostic): main front reproduced; SECONDARY FRONTS not yet algorithmically reproducible even on Bostock's catalog; tide = placeholder.
The 0.078 coherent s/n anchor is RETIRED. Merlin audits: results/merlin_audit_e2/e3_2026-09-11.md, merlin_vet_e4_redesign_2026-09-11.md.
