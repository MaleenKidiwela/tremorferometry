---
name: broadband-dvv-validated-vs-borehole
description: PGC broadband dv/v VALIDATED against co-located B011 borehole (mirror-null method); justifies the 239-station broadband/land fleet rollout
metadata: 
  node_type: memory
  type: project
  originSessionId: de9e1cec-3d07-493e-8e0a-d25339535f7a
---

**2026-07-15: broadband LFE-coda dv/v is validated against the borehole ground truth.** PGC (surface broadband, two eras: BHZ 2010-2017 + HHZ 2018-2026, per-era refs) vs B011 (borehole, ~290 m, same LFE family catalog). Both eras PASS; the broadband reproduces the borehole dv/v at station AND per-family level → the 239-station broadband/land fleet rollout is justified, and PGC is a bona-fide 16-yr broadband dv/v station.

**Reusable METHOD lesson (co-located validation):** at short separation the shared noise field makes a raw dv/v correlation ~guaranteed high even for pure-noise coda (the project's "raw dv/v = 65% artifact" finding). The valid test = MIRROR-CORRECTED, deseasoned residuals vs a **PGC-mirror-vs-B011-coda NULL** (same triggers/days/noise, zero LFE energy). Merlin caught this before I coded the wrong (raw-corr) test — [[use-merlin-for-decisions]].

**Numbers:** raw coda r 0.61/0.81 → mirror null N1 **0.04/0.04** → corrected S1 **0.38/0.43** (CIs exclude 0; slope 0.62/0.60 in [0.3,3]). Per-family matched-vs-mismatched Δ 0.31/0.38, **p=0.0005** both eras. Moving-block bootstrap (block=60d, ~55-65 eff. samples/era).

**Seam caveat (Aug-2017):** cross-seam dv/v LEVEL is UNMEASURABLE by construction (per-era refs); the swap is harmless to within-era dv/v (both eras pass independently). A real site step could NOT be tested cleanly because B011's MIRROR product has a 0-day gap May-Nov 2017 (mirror npz predates the coda 2017 gap-fix). Deseasoned-coda fallback (-0.34%, p=0.055) is marginal + mirror-uncorrected → NO step claimed. TODO if seam ever matters: rebuild B011 mirror to fill 2017.

**2026-07-15 (eve) DECISION — BROADBAND-ONLY fleet (user dropped SHORT-PERIOD EHZ).** Fleet = the 178 BH/HH stations in `data/broadband_fleet_order.csv` (region-stratified), running via `scripts/fleet_broadband_all.sh` (sequential, GPU-serialized, rolling-buffer). EHZ-only stations (~61; UW 25 of 92) are OUT. GOOD broadband stations so far: **PGC 205** (validated vs B011), **SHB 113** (56 BHZ+57 HHZ, 22-yr), **CLRS 62** (validated vs B926) — all on `figures/borehole_dvv_map.html` (33 stations). Driver traps fixed: download ALWAYS resumes (partial fragments were skipped→0 candidates); `select_families_coverage` 4th arg `rank` (fleet) vs `label` (borehole) — the fleet picker's absolute P doesn't transfer across instrument classes so select by RANK, causality is the sole verdict (Merlin caught my wrong "EHZ needs new picker" conclusion — it was a mis-reinstated P≥0.5 gate). `score_family_stacks_picker` fixed for Z-only + 40 Hz. Full recipe/decisions: `notes/FLEET_BROADBAND_PIPELINE.md`.

**2026-07-15 also: DISCOVERY path VALIDATED (fleet GO).** CLRS (broadband, no borehole in discovery) → 62 causality-certified families, 30% survival; they fire at real B926 LFE times (shift-null 4–∞×) and sit 92% within 10 km of B926's certified families. Full recipe proven: download → 40 Hz → PNSN detect → picker `--target-fs 40` → adaptive top-30k → cluster **`--min-years 1`** → densify → causality. FLEET RULES (bake into driver, these are traps): (1) `score_candidates --target-fs 40` (picker trained at 40 Hz); (2) `discover_gpu --min-years 1 --min-family-members 3 --cc-threshold 0.80 --fs 40` (default min-years 3 gave 48 vs 199 families — a silent 4× yield bug Merlin caught); (3) EVERY coincidence claim needs a +120 s shift-null (raw coincidence measures the tremor catalog, not LFE identity); (4) frozen inclusion bar ≥20 certified AND ≥15% survival — do NOT bend per-station. Fleet order: `data/broadband_fleet_order.csv` (239 stations, region-stratified, SHB #1).

Scripts: `scripts/b011_pgc_validation.py`, `scripts/pgc_mirror_build.sh`, `scripts/clrs_two_arm.sh`, `scripts/clrs_new_cohort.sh`. Results: `data/b011_pgc_validation.json`. Full write-up: `notes/2026-07-14_broadband_dvv.md` (2026-07-15 entries). B011 cert list = `data/b011_fwd_vs_rev_coda.csv` (ratio>1.5 = 119 fam), NOT a causality_cert.csv. Related: [[coda-window-2to4-correction]], [[lfe-feature-fingerprint]], [[use-merlin-for-decisions]].
