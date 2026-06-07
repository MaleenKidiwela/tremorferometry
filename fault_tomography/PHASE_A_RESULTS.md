# Phase A — Results log

## A1: fault-patch × station overlap feasibility (2026-06-06)

**Question:** does our data have the (fault-patch × station) bipartite redundancy the deep-lobe
inversion needs? The common-mode deep lobe is only resolvable where a patch is imaged by ≥3 stations
spanning a decent azimuth range (theory/01 §6, theory/02 §6).

**Method:** `data_assembly/a1_patch_overlap.py` — every dv/v patch id `<lat>_<lon>__c<n>` carries its
fault-cell location; count distinct stations per fault cell at several patch sizes. Read-only on the 21
stations that currently have `daily_dvv_<STA>_coda_1to4.csv` (PB boreholes + HDW; the GNW/CPW/PGC/NLLB
`_perera` files are not yet folded in → real overlap is a bit higher).

**Result — resolution vs. overlap tradeoff:**

| patch size | #patches | ≥2 sta | ≥3 sta | ≥4 sta | max sta |
|---|---|---|---|---|---|
| 0.05° (~6 km)  | 1074 | 215 | 30  | 2  | 5 |
| **0.10° (~11 km)** | 570 | 276 | **117** | 34 | 5 |
| 0.15° (~17 km) | 372 | 246 | 131 | 50 | 8 |
| 0.20° (~22 km) | 231 | 163 | 113 | 54 | — |
| 0.30° (~33 km) | 126 | 97  | 71  | 51 | — |

**Verdict: FEASIBLE at ~10–17 km fault-patch resolution.** The native 0.05° grid over-fragments
(finer than both the deep-lobe radius and the cross-station sharing scale), giving a misleadingly sparse
30 patches. At **0.10–0.15° (11–17 km)** — which matches the physical deep-lobe radius √(2D·t_min) for
crustal D and few-second lapse — we get **~120–130 fault patches imaged by ≥3 stations** (34–54 by ≥4).
That is a viable tomographic target, and it's a floor: it grows when we (a) fold in GNW/CPW/PGC/NLLB,
(b) finish the remaining stations, (c) optionally do a deliberate *shared-patch* densification.

**Decisions this settles / informs:**
- **Fault-mesh resolution → ~0.10–0.15° (11–17 km)** along the interface (was an open decision in PLAN §6).
- **Selection-strategy note for any future densification:** the per-station coverage selection optimized
  each station's *own* azimuthal spread, NOT cross-station patch sharing. To strengthen the deep-lobe
  inversion we could add a pass that deliberately images a common set of high-value patches from every
  station that can see them (boosts ≥4-station coverage where it matters).
- The ~117 well-covered patches (0.10°) define where the fault map will be *data-resolved*; the
  fault-surface Laplacian (theory/02 §5) fills the rest via regularization.

Figure: `figures/a1_patch_overlap.png` (0.10° patches colored by #stations; black rings = ≥3-station).

**Next (A2):** assemble the full `(patch × station × lapse × epoch)` data tensor on a 0.10–0.15° mesh,
folding in the `_perera` stations; then A3 multi-lapse-window dv/v on the stacks, A4 medium params, A5 mesh.

## A1b: co-located stations — use the borehole, drop the redundant surface partner (2026-06-06, user)

For the bipartite inversion, two **co-located** stations view the same fault patches from nearly the same
azimuth → the second adds NO azimuthal diversity, only noise and over-weighting. So for each co-located
pair, keep ONE station (prefer the clean single-era borehole) in the tomography set:
- **CPW ↔ B018** → use **B018** (clean borehole). DROP CPW from the tomography. CPW's only unique offering
  is its pre-2007 years, but pre-2007 has too few stations for any spatial inversion anyway → no loss.
  (CPW remains a valid STANDALONE dv/v curve — its clean 2001–2011 decade — just not a tomography input.)
- **PGC ↔ B011** (≈49 km, near-co-located): decide per-pair — PGC has the longer/cleaner broadband record,
  B011 is the borehole. Likely keep PGC for its length; revisit at A2.
- **GNW, NLLB** have no truly co-located borehole → keep (GNW especially, for its 1995–2026 length).
Action for A2: when building the data tensor, deduplicate co-located stations by patch+azimuth before
counting "≥3-station" coverage, so the redundancy metric isn't inflated by co-located pairs.

## BOREHOLE-FIRST BUILD (2026-06-06) — machinery validated end-to-end

Decision (user): build the tomography on the 22 PB boreholes first (uniform single-era, all-time-ref,
no per-era/co-location complications), then add the heterogeneous stations.

- **A2** (`data_assembly/a2_borehole_tensor.py`): 8.16M dv/v rows → **603 fault cells (0.10°), 123 with
  ≥3 stations**; monthly (cell×station×epoch) tensor 2005–2026. *Catch:* Slab2 INPUT DB mixes constraint
  types (bathymetry 1–5 km, deep tomography 60–440 km); naive interp put 27 tremor cells <15 km. Fixed by
  using only interface-tracking constraints (drop BA/TO, keep 10–70 km band, clip 12–55) → cells 14–55 km,
  median ~30 (right for Cascadia LFEs). Depth is APPROXIMATE (smooth Slab2 *output* grid firewall-blocked).
- **B kernel** (`kernels/kernel.py`): half-space diffusion kernel by 1-D quadrature. *Catch:* theory note
  eq.14 normaliser drops a (πa)^1.5 factor → ∫K dV≈16000 instead of 1; restored from Chapman–Kolmogorov →
  ∫K dV=1.03, two-lobed (deep@source + surface@station). Also added the **single-scattering ellipsoid
  kernel** (eq.15).
- **KEY PHYSICS:** our **1–4 s coda is EARLY coda** → must use the single-scatter ellipsoid kernel, NOT
  diffusion (which smears over √(2Dt)~40 km). Switching kernels took captured on-fault sensitivity
  **η 0.04 → 0.98** and checkerboard recovery **0.18 → 0.77**.
- **C forward** (`inversion/forward.py`): G_f = −K(ellipsoid) per (cell,station) pair over the 603 nodes;
  1064 data pairs. Medium (A4, literature first pass): β=3.5 km/s, ℓ*=40 km → D=47 km²/s.
- **D validation** (`figures/checkerboard_test.png`): planted ±1% checkerboard → recovered **corr 0.77 at
  the 123 ≥3-station cells** with LSQR + graph-Laplacian regularisation. **The end-to-end machinery works
  and resolves the interface where the boreholes have ≥3-station coverage.**

Remaining for a real result: add per-station **site terms G_s** (fault-vs-site separation — the seasonal
signal must land on sites, not the fault), then invert REAL epochs → first δβ/β map → 4-D loop. Then fold
in the non-borehole + COLT/COR/GNW stations (northern azimuths + pre-2007 reach).
