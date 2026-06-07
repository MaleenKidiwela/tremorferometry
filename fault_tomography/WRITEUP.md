# Time-Lapse Shear-Velocity Tomography of the Cascadia Plate Interface from LFE-Coda Interferometry

**Working methods/approach write-up.** Companion theory in `theory/01_kernel.md` (volumetric sensitivity
kernel) and `theory/02_fault_tomography.md` (fault-plane inverse problem); feasibility in
`PHASE_A_RESULTS.md`; engineering plan in `PLAN.md`. Draft — numbers current as of 2026-06-06.

---

## 1. Summary

We use low-frequency earthquakes (LFEs) on the Cascadia subduction interface as **fixed, repeating,
deep seismic sources** and measure the time-varying travel-time change (dv/v) of their **scattered S-wave
coda** at a margin-wide network of surface and borehole stations. Because each LFE family is a point
source sitting *on the fault* at ~30–40 km depth, the coda sensitivity kernel for any surface receiver has
a sensitivity lobe **pinned in the fault zone at the source** — and that lobe is *common to every station
that records the same family*. Inverting the full set of (LFE-family × station × lapse-window × epoch) coda
observations therefore yields a **4-D image of shear-velocity change δβ/β(x₀, y₀, t) on the plate
interface itself**, plus a by-product map of where the LFE sources migrate. To our knowledge this is the
first attempt to turn LFE coda into a spatially-resolved, time-lapse velocity image *of the megathrust*.

## 2. Why this geometry is new

Classical coda-wave interferometry and ambient-noise monitoring assume a **surface source and surface
receiver**, and treat the coda as a mix of a 2-D surface-wave field (early coda) and a 3-D body-wave field
(late coda). The resulting depth sensitivity has two lobes beneath each station and decays with depth — it
is hard to attribute a measured dv/v to a specific depth, let alone to the fault.

Our geometry breaks three of those assumptions, each of which *helps*:

1. **The source is deep and repeating.** An LFE family is an effectively fixed point source you can
   re-record for decades — the ideal repeatable probe.
2. **A deep source radiates no surface waves.** At depth ≫ wavelength the source-side field is pure 3-D
   S-body-wave intensity; the only place a surface-wave channel can enter is at the free surface near the
   receiver. The source leg is single-mode and clean.
3. **The receiver is at the surface**, carrying the free-surface amplification and being the *only* leg
   where body↔surface coupling occurs.

The consequence (derived in `theory/01_kernel.md`) is a kernel with **one lobe pinned in the fault zone at
the source, one lobe at the surface station, and a connecting tube.** Two properties of the deep lobe make
the network inversion possible:

- **Common-mode.** Every station that records a given family shares the *same* deep lobe, so a real
  fault-zone velocity change appears as a dv/v signal **coherent across all of that family's stations**,
  while a shallow change near one station appears only in that station's data. This coherent-vs-local
  contrast is the core discriminant.
- **Early-coda = sharpest fault image.** Because the source is *already at depth*, the earliest usable
  S-coda (least-scattered, most direct paths) weights the deep lobe most tightly; later coda diffuses and
  shifts weight toward the shallow receiver side. This is the **opposite** of ambient-noise monitoring,
  where one pushes to *late* coda to reach depth. Scanning lapse-time windows thus gives a depth-resolving
  family of kernels, and the earliest windows carry the fault signal.

## 3. Data

**Stations (27 processed, 40.5°N–49.2°N along the margin):** 22 PBO/NOTA borehole seismometers (network
PB, EHZ 100 Hz, single stable instrument era → all-time reference) spanning the whole transect, plus 5
long-running surface/broadband stations (UW short-period GNW/HDW/CPW, CN broadband PGC/NLLB; multi-era →
per-era referencing). Borehole stations are quiet and stable — most yield coda cross-correlation cc ≈
0.97–0.99. The surface stations add **multi-decade reach** (GNW/HDW back to 1995/1999, pre-dating the
borehole network) and independent azimuths in the north.

**Sources:** LFE families from the PNSN tremor catalogue, located on the interface (depths from Slab2,
Hayes et al. 2018; ~30–40 km, deepening downdip). Each family is detected and stacked per station; per
station we keep a coverage-balanced + top-10%-SNR selection of ~30–117 families.

**Observable (per family a, station b, day):** the coda stretch δτ/τ — the relative travel-time change of
the 13-s long-window stack (−3..+10 s about the LFE), measured by coda stretching in a 1–4 s lapse window
against an all-time (single-era) or per-era (multi-era) mean reference. We also retain the maximum
correlation cc as the seed for the decorrelation/source-migration channel. Each family's stacks form a
daily time series spanning the station's record.

**Bipartite structure (the key data property):** because discovery is driven by one master catalogue, a
given LFE patch is the *same physical source* at every station that records it. The patch identifiers embed
the source's fault grid cell, so the dataset is naturally a **(fault-patch × station × lapse-window ×
epoch)** tensor — a row effect per fault patch and a column effect per station — which is what makes the
fault-vs-site separation identifiable (§5).

**Feasibility (Phase A1).** The deep-lobe inversion needs fault patches imaged by ≥3 stations spanning a
range of azimuths. Counting cross-station overlap on the real data: at the native 0.05° (~6 km) grid only
30 patches reach ≥3 stations, but that grid over-fragments below the physical resolution. At a
physically-appropriate **0.10–0.15° (11–17 km) fault patch** — matching the diffusive deep-lobe radius
√(2D·t) — **~117–131 patches are imaged by ≥3 stations** (34–54 by ≥4). This sets the fault-mesh
resolution and confirms the inversion is feasible; coverage is densest (and deepest in time) in the north,
borehole-only and post-2007 in the south.

## 4. Forward model

For a relative velocity perturbation δβ/β localized near model point **x**, the coda travel-time change
recorded from source **s** at receiver **r** in lapse window t is (Pacheco & Snieder 2005)

  δτ(t)/t = −∫_V K(x; s, r, t) · (δβ/β)(x) d³x,  with ∫_V K dV = 1,

where the sensitivity kernel is the occupation-time (Markov-bridge) product of intensity propagators,

  K(x; s, r, t) = [1 / (t · P(s,r,t))] ∫₀ᵗ P(s,x,u) P(x,r,t−u) du.

Using the half-space diffusion propagator with an image source for the reflecting free surface gives an
explicit kernel evaluated by **one 1-D numerical quadrature per model point** (`theory/01_kernel.md`,
eq. 13–14) — no full-wavefield simulation. For the earliest windows we use the single-scattering kernel,
a prolate ellipsoid with foci at the source and receiver. Intrinsic attenuation Q cancels in the kernel
*shape* (it only sets how late in lapse time the coda stays above noise), so the kernel needs only the
diffusivity D = βℓ*/3 and the geometry. The diffusion form is most trustworthy beyond ~one transport
mean free path of each sensor — which, for a ~30 km source and crustal ℓ*, is **exactly the fault-zone
region we care about**; within ℓ* of the surface station we splice in a radiative-transfer kernel.

## 5. Inverse problem

Collapsing the volumetric field onto the interface Σ (the velocity changes we seek — fluids, pore
pressure, damage/healing — are fault-localized) turns K into an effective 2-D fault-plane kernel κ(ξ) =
K·h(ξ), with a **captured-sensitivity fraction η ≤ 1** per datum that we use as a data weight (low-η pairs
are off-fault-dominated and down-weighted). Discretizing the interface on the 0.10–0.15° mesh, the full
forward map is

  d = G_f m_f + G_s m_s + G_p m_p + ε,

where **m_f = δβ/β on the fault (target)**, m_s = per-station shallow *site* terms, m_p = a coarse
off-fault *path* field. The design matrix G_f has a **near-diagonal self-imaging core** — each family most
strongly illuminates its own patch — plus controlled neighbour leakage. Identifiability comes from the
bipartite structure: a true fault change is the part of a family's signal coherent across *all* its
stations; a site change is what is common across *all* families at a station. We solve the regularized
least-squares / Bayesian system (LSQR) with a fault-surface Laplacian (optionally anisotropic, smoother
along strike than dip) and report the resolution matrix plus checkerboard tests to map where the
family+station layout actually resolves the interface.

**Two orthogonal channels.** The same coda yields two non-aliasing observables (Snieder 2006): a
**coherent stretch growing ∝ lapse time → velocity change** (the inversion above), and an **incoherent,
lapse-time-flat variance → source relocation**, giving each family's migration distance δ = √3·β·σ_τ. We
estimate the velocity field and the per-family source migration separately, then optionally jointly.

**4-D.** Every observable is a time series, so we invert per epoch (≈30–60 day windows, with temporal
smoothing) to produce **δβ/β(x₀, y₀, t)** — a movie of interface velocity change through episodic
tremor-and-slip (ETS) cycles — alongside a map of where the LFE sources themselves migrate.

## 6. Why a single station is not enough (a worked caution)

A single-station dv/v change is fundamentally ambiguous: it could be a fault-zone change, a shallow site
change, or data degradation, and one station cannot separate them. Our station B201 illustrates this — it
shows a multi-year dv/v decrease, but the change is **~4× larger than at its nearest neighbour 19 km away
and absent further out**, i.e. not network-coherent, so it is most likely local rather than a megathrust
signal. The common-mode-across-stations test that flags this is precisely the discriminant the fault-plane
inversion formalizes. The whole method exists to convert such ambiguous single-station curves into a
proper fault-vs-site-vs-noise attribution.

## 7. Status and plan

- **Phase A (data + feasibility) — in progress.** 27 stations processed into per-family-day dv/v;
  cross-station overlap confirms feasibility at ~11–17 km resolution (§3). Remaining A-steps: assemble the
  full (patch × station × lapse × epoch) tensor; re-measure a *ladder* of lapse windows from the existing
  stacks (early windows for the fault, per §2); estimate D and Q from coda-envelope decay; build the Slab2
  fault mesh and embed the families.
- **Phase B–E.** Implement and unit-test the kernel module (eq. 13–15); assemble G with site/path nuisance
  blocks and η-weighting; run checkerboard/resolution tests; invert a single epoch for a first δβ/β map and
  validate against a full-wavefield synthetic with a planted on-fault patch; then iterate over epochs for
  the 4-D map and the source-migration channel.

## 8. Significance

If it works, this delivers the first **spatially-resolved, time-lapse image of shear-velocity change on the
Cascadia plate interface**, driven by sources that sit on the fault itself. Targets include imaging
velocity changes that accompany ETS / slow slip, testing for pre-seismic changes around large events (the
project grew from a Nisqually-2001 precursor question on the 30-year GNW record), and mapping where LFE
sources migrate — a second fault-process observable obtained for free from the same coda.

---

*References as in `theory/01_kernel.md` and `theory/02_fault_tomography.md` (Pacheco & Snieder 2005, 2006;
Snieder 2006; Snieder & Vrijlandt 2005; Margerin et al. 2016, 2019; Obermann et al. 2013b, 2016; Paasschens
1997; Barajas et al. 2022; Aster et al. 2018; Tarantola 2005; Hayes et al. 2018 [Slab2]).*
