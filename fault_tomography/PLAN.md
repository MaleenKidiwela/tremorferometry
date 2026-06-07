# Fault-Plane δβ/β Tomography from LFE Coda — PLAN

**Workspace for the NEXT phase**: turn the per-station dv/v curves (the station-running phase)
into a **4-D tomographic movie of shear-velocity change δβ/β(x₀,y₀,T) on the Cascadia plate
interface**, plus a **map of where the LFE sources migrate**. Foundation theory lives in
`theory/01_kernel.md` (the volumetric sensitivity kernel) and `theory/02_fault_tomography.md`
(assembling all pairs into a fault-plane inverse problem).

Status: PLANNING. Station-running (the data-collection phase) is ~25 stations done and continuing
in parallel; nothing here blocks on it, but it gets sharper as more stations land.

---

## 0. The deliverable, in one line
A movie of δβ/β on the plate interface through ETS cycles + a per-family source-migration map,
inverted from the (LFE-family × station × lapse-window × epoch) coda observables we already produce.

---

## 1. The method in operational terms (what the two notes actually require)

**Geometry that makes this new.** Deep repeating LFE source (~30 km, *on the fault*) + surface
receiver. The coda sensitivity kernel becomes: a **deep lobe pinned in the fault zone at the
source** + a **surface lobe at the station** + a **connecting tube**. The deep lobe is
**common-mode across every station that records the same family** — that shared deep lobe is the
entire lever that separates a fault-zone change from a shallow site change.

**Two observables per (family a, station b, lapse window i, epoch T):**
- **(a) stretch δτ/τ → velocity change.** Coherent, grows ∝ lapse time, common-mode across stations
  → fault-zone δβ/β. **We already measure exactly this** (our `dvv`).
- **(b) decorrelation 1−CC / residual variance σ_τ → scattering change OR source migration.** We
  have `cc_max` per measurement as the seed.

**Velocity vs source migration are orthogonal fingerprints** (Note1 §8): velocity = coherent stretch
∝ t; source displacement δ = zero-mean but variance σ_τ² = δ²/(3β²), *flat* in t. Fit each coda
window for both → they don't alias. δ = √3·β·σ_τ.

**Inverse problem** (Note2 §5): `d = G_f m_f + G_s m_s + G_p m_p + ε`
- `m_f` = fault δβ/β map (target), `m_s` = per-station site terms, `m_p` = coarse off-fault path field.
- The (family × station) **bipartite structure** makes fault-vs-site identifiable: a true fault change
  is the part of family a's signal *coherent across all stations recording it*; a site change is what's
  common across all families recorded at station b.
- Regularize with a **fault-surface Laplacian** (anisotropic: smoother along-strike than down-dip).
  Solve LSQR/CGLS; inspect resolution matrix + checkerboard tests.

**THE KEY INSIGHT that flips our intuition** (Note1 §6, Note2 §6): because the source is *already at
depth*, the **earliest usable S-coda gives the SHARPEST fault image** (tightest deep lobe). This is
the OPPOSITE of ambient-noise / surface–surface monitoring, where you push to *late* coda for depth.
→ We should measure a **ladder of lapse windows** and expect the early ones to carry the fault signal.
Q_i only sets how late we can go before coda hits noise; it does **not** change kernel shape.

---

## 2. What we ALREADY have (plugs straight in)

| Asset | File / form | Role in the inversion |
|---|---|---|
| Per-(family,station,day) **coda stacks** | `data/long_window_daily_<STA>.npz` (13 s, −3..+10 s @ 40 Hz; arrays `stacks`,`patches`,`dates`,`n_det`,`t`) | **Raw material.** Covers early→late coda → multi-lapse-window ready. KEEP FOREVER. |
| Per-(family,station,day) **stretch** | `data/daily_dvv_<STA>_coda_1to4.csv` (`patch,date,dvv,cc_max,dvv_err,n_det`) | = d_abi for the 1–4 s window. `patch` = `<lat>_<lon>__c<n>` → **fault location is embedded**. |
| **Decorrelation seed** | `cc_max` column | Channel (b) start. |
| **Source/family locations** | `data/<sta>_pnsn_families_100km.summary.csv` (`lat,lon,snr`) + Slab2 depth | X_a = X(ξ_a) on the fault. |
| **Fault geometry** | `data/cas_slab2_input_04-18.csv` (Slab2 Cascadia) + `data/station_slab2_depth.csv` | Build the dipping interface mesh Σ. |
| **Station locations** | `scripts/plot_pb_path_map.py` COORDS + original 5 | R_b. |
| **Coverage / array** | ~25 stations, tens of families each, 41.5–49.2°N | Azimuthal looks per fault patch. |

**Crucial latent structure:** all discovery is driven by the SAME master catalog
(`catalogs/pnsn_tremor_cascadia_full.csv`), so a given LFE patch is the *same physical source* at
every station that records it. Same 0.05° fault grid cell (`<lat>_<lon>` prefix of the patch id) =
same patch. → The family×station bipartite matrix is already implicit in our data; we just assemble it.

---

## 3. Gaps to build (G#)

- **G1 — Bipartite assembly.** Pivot per-station dv/v into a `(fault_patch × station × lapse × epoch)`
  tensor; match families across stations by fault grid cell. *Buildable NOW.*
- **G2 — Multi-lapse-window stretch.** Re-measure dv/v on the existing stacks for a window ladder
  (e.g. 0.5–1.5, 1–2, 2–4, 4–7, 7–10 s). No new downloads. *Buildable NOW from `.npz` stacks.*
- **G3 — Medium parameters.** D = βℓ*/3 and Q_i from coda-envelope decay of the stacks (per region).
  *Buildable NOW.*
- **G4 — Fault mesh + family embedding.** Slab2 → Σ mesh + place every family. *Buildable NOW.*
- **G5 — Kernel module.** Half-space diffusion kernel (Note1 eq 13–14, one 1-D quadrature per model
  point) + single-scatter ellipsoid (eq 15) for early windows. Unit tests: ∫K dV = 1; homogeneous
  perturbation → δτ/t = −δβ/β. *Buildable NOW.*
- **G6 — Two-channel measurement.** Separate coherent stretch (∝ t) from incoherent variance (flat) per
  (family,station). *Needs G2.*
- **G7 — Inversion.** Assemble G = [G_f | G_s | G_p], weight C_d⁻¹ scaled by captured-sensitivity η_abi,
  fault-Laplacian regularization, LSQR; resolution matrix + checkerboard. *Needs G1–G5.*
- **G8 — 4-D.** Per-epoch (≈30/60-day) inversion with temporal smoothing → the movie + source migration.
  *Needs G7.*

---

## 4. Ordered work plan

**Phase A — buildable NOW, no new data (do in parallel with finishing the station queue):**
- **A1 Feasibility check** *(do first)*: cross-station fault-patch overlap — how many fault cells are
  imaged by ≥3 stations, and with what azimuthal spread? This decides whether the common-mode deep lobe
  is actually resolvable. If overlap is thin, it reshapes which stations we still run.
- **A2** Assemble the `(patch × station × lapse × epoch)` data tensor from existing dv/v + stacks (G1).
- **A3** Multi-lapse-window dv/v re-measurement on the `.npz` stacks (G2).
- **A4** Medium params D, Q from coda envelopes (G3).
- **A5** Fault mesh + embed families from Slab2 (G4).

**Phase B** — Kernel module + unit tests (G5).
**Phase C** — Forward operator G (fault + site + path blocks, η-weighting); checkerboard/resolution tests (G7 setup).
**Phase D** — Static (single-epoch) inversion → first δβ/β fault map; validate against a full-wavefield synthetic with a planted on-fault patch.
**Phase E** — 4-D epoch inversion → movie; per-family source-migration channel (G6, G8).

---

## 5. Implications for the STATION-RUNNING happening now (avoid rework)

1. **KEEP every `long_window_daily_<STA>.npz`** — it's the raw material for every lapse window we'll
   ever want. (We already keep these; raw traces are fine to delete.) ✓
2. **Azimuthal coverage per fault patch is what resolves the deep lobe.** A cluster of stations over a
   shared, rich LFE patch is worth MORE than an isolated new patch. → bias remaining station picks
   toward adding **new azimuths to already-imaged patches**, not only new ground.
3. **Early coda matters.** Our stacks start pre-S (−3 s) so early coda is captured ✓. We may add
   narrower/earlier windows in A3.
4. **Same master catalog → shared physical families** → the bipartite structure is free; just assemble it.

---

## 6. Open decisions to settle before Phase C
- Coda-window ladder: how early do we trust the S-coda, and how many windows?
- Reference handling (all-time vs per-era) consistent across windows and with the kernel's lapse-time weighting.
- Fault-mesh resolution: patch size vs family spacing vs deep-lobe radius √(2D·t_min) (Note2 eq 10).
- Regularization: anisotropy of the fault Laplacian (strike vs dip).
- Diffusion vs RTE/single-scatter kernel per lapse window (early windows need the ellipsoid/RTE form).

---

## 7. Folder layout (this workspace)
```
fault_tomography/
  PLAN.md                  <- this file
  theory/
    01_kernel.md           <- Note1: volumetric sensitivity kernel (deep source / surface receiver)
    02_fault_tomography.md <- Note2: all-pairs fault-plane inverse problem
  (to come) data_assembly/ kernels/ inversion/ figures/
```
