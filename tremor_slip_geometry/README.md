# tremor_slip_geometry — PNSN tremor vs Gualandi (2025) slip: relative geometry through 2010-2026

Self-contained study folder (user 2026-09-10: "do it in a separate folder", "keep everything documented", "no dv/v yet").
- PLAN.md       pre-registered questions, quantities, controls, steps S1-S6 (dv/v bridge deferred)
- NOTES.md      running ledger: every run, number, issue, decision (append-only, dated)
- scripts/      one script per step (s1_coregistration.py, s2_depth_occupancy.py, s3_efficiency_maps.py, s4_fronts.py,
                s5_site_recurrence.py, s6_model_checks.py); each has a docstring stating inputs, definitions, outputs
- results/      CSV tables written by the scripts (s<step>_*.csv)
- figures/      PNG figures (s<step>_*.png)

Inputs (read-only, outside this folder):
- catalogs/pnsn_tremor_cascadia_full.csv                          PNSN tremor master (time, lat, lon; depth is a placeholder)
- fault_tomography/inversion/res_catalog_sse/slip_patches_5mm.{csv,npz}          unbridged slip-patch catalog + labels (T x 742), inc30_mm
- fault_tomography/inversion/res_catalog_sse/slip_patches_5mm_mirror.{csv,npz}   back-slip control catalog
- fault_tomography/inversion/res_catalog_sse/tremor_vs_slip_patches.csv          per-patch tremor confirmation (confirmed, seg, catalog)
- fault_tomography/inversion/res_catalog_sse/components_at_cells_trailing30.npz  A_raw (T x 13 vbICA amplitudes), days
- /home/jovyan/cascadia-sse-strain  (mesh data/cascadia_mesh.txt; sse_strain.load_fault_ascii / load_solution /
  component_slip_patterns -> Ld; sse_strain.megathrust.slab_depth_interpolator)   thrust slip per triangle = -(Ld @ A_raw.T) [m]
Conventions: thrust-positive slip; slip start/end are trailing-30-d threshold crossings (slip within (t-30, t] > 5 mm);
local along-strike / down-dip directions from the Slab2 depth gradient; segments north >= 46N, central 44-46N, south < 44N.
Execution: coding by Opus sub-agents from specs written by the main session; main session reviews, integrates, documents.
