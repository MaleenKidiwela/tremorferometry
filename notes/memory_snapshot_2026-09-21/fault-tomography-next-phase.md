---
name: fault-tomography-next-phase
description: The NEXT phase after station dv/v runs — invert all LFE-coda dv/v into a 4-D δβ/β map on the Cascadia plate interface
metadata: 
  node_type: memory
  type: project
  originSessionId: 400c26f9-17ff-46fc-914c-9c6a7671906d
---

**After the per-station dv/v runs, the project's payoff phase is FAULT-PLANE TOMOGRAPHY.** Workspace: `fault_tomography/` (PLAN.md + theory/01_kernel.md + theory/02_fault_tomography.md — the two theory notes the user supplied 2026-06-06, persisted there because uploads are ephemeral).

**Goal:** invert every (LFE-family × station × lapse-window × epoch) coda observable into a **4-D movie of shear-velocity change δβ/β(x₀,y₀,T) on the plate interface** + a per-family **source-migration map**.

**The physics (why it's novel):** deep repeating LFE source (~30 km, on-fault) + surface receiver → coda kernel = **deep lobe pinned in the fault zone** (common-mode across all stations sharing that family) + surface lobe + tube. Common-mode deep lobe is what separates fault-zone change from shallow site change. Half-space diffusion kernel = one 1-D quadrature per model point (Note1 eq 13-14). **KEY INSIGHT (opposite of ambient-noise):** EARLIEST usable S-coda gives the SHARPEST fault image (source already at depth). Two orthogonal channels: coherent stretch ∝t → velocity; incoherent variance flat-in-t → source migration (δ=√3·β·σ_τ). Inverse: d=G_f m_f+G_s m_s+G_p m_p; bipartite (family×station) structure makes fault-vs-site identifiable; fault-surface Laplacian reg; LSQR.

**What we already have that plugs in:** `data/long_window_daily_<STA>.npz` (coda stacks = raw material, KEEP FOREVER), `data/daily_dvv_<STA>_coda_1to4.csv` (=stretch d_abi; patch id `<lat>_<lon>__c<n>` embeds fault location), cc_max (decorrelation seed), Slab2 (`data/cas_slab2_input_04-18.csv`) for the fault mesh. All discovery is PNSN-catalog-driven → families are SHARED physical sources across stations → the bipartite matrix is latent, just needs assembly.

**Implications for ongoing station runs:** (1) keep every long_window npz; (2) **azimuthal coverage per fault patch resolves the deep lobe** → station CLUSTERS over shared LFE patches > isolated new patches; bias remaining picks toward new azimuths on already-imaged patches; (3) early coda matters (stacks start −3s ✓).

**First thing to compute (A1 feasibility):** cross-station fault-patch overlap — how many 0.05° fault cells are imaged by ≥3 stations with azimuth spread. Computable now from the dv/v csvs (match on `<lat>_<lon>` patch prefix). Decides if the common-mode deep lobe is resolvable. See [[margin-station-status]].

**UPDATE 2026-06-09 — BUILT & VALIDATED (multi-window + near-surface).** Joint inversion of 3 coda windows
(1-3/2-4/3-5 s, deseasoned 30-cal-day data) with a SHARED fault field + per-(station,window) site terms:
`fault_tomography/inversion/invert_multiwin.py`. **Checkerboard 0.88** (vs 0.77 single-window — depth leverage
is real; fault = across-window common mode). Synthetic spike test (`synth_resolution.py`): a 0.3% fault change
is recoverable 77-101% margin-wide, tight where ≥4 stations / smeared at 1-2 (noise floor ~0.1%). **KEY: the
site terms ARE the near-surface field** — the seasonal lives there (B036 0.32% annual), the deep fault stays
clean (raw≈deseasoned) → one inversion = TWO products (deep fault + shallow/site). Resolution bottleneck =
fault/site separation = stations-per-cell, so MORE SOUTHERN STATIONS (not more windows) is the coverage fix.
PARKED: null test (southern hot real vs artifact); SVD strictly-causal rerun (SVD bleeds future + smears sharp
drops). Full arc in `notes/2026-06-09_Notes.md`. See [[coda-window-2to4-correction]].

**UPDATE 2026-06-09b — kernel physics, cross-station null, 3D-volume verdict.**
(1) **Kernel = two-lobed.** The single-scatter sensitivity is a 3D prolate ellipsoid (foci=source,receiver);
3D MASS profile is U-shaped: ~32% at the SOURCE lobe (fault), ~30% at the RECEIVER lobe (shallow-under-station),
only ~12% in the mid-crust. So ~⅓ of coda dv/v is fault-sensitive, ~⅓ shallow — surface signal is HALF the
measurement, not a contaminant; fault/surface separation is mandatory. (`smoke_kernel_crosssection.png`,
`smoke_kernel_massprofile.png`, `volume3d_checkerboard.py`/`synth_resolution.py` for the math.)
(2) **Cross-station null (model-free stability).** Same source @48.4N seen by 6 stations → dv/v UNCORRELATED
(mean pairwise −0.03 deseasoned). No common deep change ≳0.1–0.3% — the stable-interface result without any
inversion. Bipartite is real at the LOCATION level (family IDs are per-station; 58 exact 0.05° sources on ≥3
stations, 196 cells on ≥3). (`smoke_sharedsource_48.4N.png`.)
(3) **3D VOLUME: fine NO, coarse YES.** Free 3D volume doesn't resolve at fine scale (checkerboard 0.11 — depth is
the unconstrained per-datum direction) but DOES at coarse scale (1.8°/60km → 0.86) INCLUDING shallow-vs-deep
separation. So a COARSE 3D (~150–200 km horizontal, 3 depth zones 10/26/42 km, no lumped site terms — shallow layer
= near-surface resolved) is meaningful. `invert_volume3d.py` → `fault_volume3d.npz` + `fault_movie_volume3d.mp4`;
`scripts/build_volume3d_html.py` → `cascadia_volume3d.html` (interactive: pick depth zone + slide time). Always
resolution-test before rendering a volume. Full arc: `notes/2026-06-09_Notes.md` (now §§1–9).
