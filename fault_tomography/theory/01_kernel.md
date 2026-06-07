# Coda-Wave Sensitivity Kernel for a Deep Repeating Source and a Surface Receiver

**S-wave coda monitoring of LFE families (Cascadia, source depth ≈ 30 km) — forward and inverse problem**

> Persisted from upload (Note1). Companion: `02_fault_tomography.md`. Operational summary in `../PLAN.md`.

---

## 0. What is actually new here

Classical coda-wave interferometry (CWI) and passive image interferometry build their travel-time
sensitivity kernels for a **surface source / surface receiver** geometry, treating the coda as a
*composite* of a 2-D surface-wave field (early coda) and a 3-D body-wave field (late coda) (Obermann
et al. 2013b, 2016; Margerin et al. 2016). That composite is why the standard kernel has two strong
lobes *beneath each station* and weak sensitivity in between, decaying with depth.

This geometry breaks three of those assumptions:

1. The source is **deep** (z_s = d ≈ 30 km), inside the fault zone, and is a *repeating* LFE family —
   an effectively fixed point source you can re-record indefinitely.
2. **A deep source does not radiate surface waves.** At d ≫ λ the source-side field is *pure 3-D
   body-wave (S) intensity*. There is no 2-D component on the source leg ("only S in the coda"). The
   only place a 2-D channel enters is the free surface near the receiver.
3. The **receiver is at the free surface**, so the receiver leg carries free-surface amplification and
   is the *only* leg where body↔surface coupling can occur.

Consequence: a kernel with one lobe **pinned in the fault zone** (at the source) and one lobe at the
surface station, joined by a sensitivity tube. Because every surface station shares the *same* deep
source, the deep lobe is **common-mode across the whole network** — exactly what lets you separate a
fault-zone change from a shallow site change.

---

## 1. Assumptions and notation

| Symbol | Meaning |
|---|---|
| s = (0,0,d) | source position (LFE family centroid), depth d |
| r = (r_h, 0) | receiver position, on the free surface z=0 |
| x = (x_h, z) | scatterer / model point, z ≥ 0 downward |
| β | shear-wave (S) speed; single propagating mode |
| ℓ | scattering mean free path; τ = ℓ/β |
| ℓ* | transport mean free path, ℓ* = ℓ/(1−⟨cosθ⟩) |
| D | energy diffusivity, D = βℓ*/3 (3-D) |
| σ_a = 1/τ_i | intrinsic absorption rate (τ_i = Q_i/ω) |
| t | lapse time in the coda (from LFE origin time) |
| P(x₁,x₂,t) | mean intensity (energy-density) propagator |
| K(x; s,r,t) | travel-time sensitivity kernel |
| δβ/β | relative shear-velocity perturbation field (target) |
| δ | source-location displacement between two LFE realizations |

Working assumptions: single S mode, isotropic background, statistically homogeneous scattering (RTE
applies; diffusion at t ≫ τ); free surface = perfectly reflecting for S energy (Neumann); weak localized
perturbations → first-order (Born) sensitivity; coda windows long vs period, short vs medium-change time.

---

## 2. The intensity propagator

**Diffusion limit (t ≫ τ):**
∂P/∂t = D∇²P − σ_a P + δ(x−x₀)δ(t),  D = βℓ*/3.   (1)

Infinite-medium Green's function:
P_∞(x,x₀,t) = (4πDt)^(−3/2) · exp[−|x−x₀|²/(4Dt)] · e^(−σ_a t),  t>0.   (2)

**Early coda — keep the ballistic term (Paasschens 1997):**
P_RTE(R,t) = [e^(−βt/ℓ)/(4πR²)]·δ(t−R/β)  +  P_sc(R,t)·Θ(t−R/β),  R=|x−x₀|.   (3)
The halo P_sc interpolates single-scattering → Gaussian (2). The factor e^(−βt/ℓ) is the exponential
mean-free-path decay. Use (3) in early S-coda; (2) is its long-lapse limit.

**Half-space via images (reflecting free surface):**
P_½(x,x₀,t) = g(x,x₀,t) + g(x,x₀',t),  with x₀'=(x_{0h},−z₀), g=(4πDt)^(−3/2)exp[−|x−x₀|²/4Dt].  (4)
- Source leg (x₀ = s = (0,0,d)): two-term (real at d + image at −d); **no surface-wave term**.
- Receiver leg (z₀=0): the two images coincide → **factor-2 free-surface amplification** (5).
- Direct source→receiver: P_½(s,r,t) = 2(4πDt)^(−3/2) exp[−(r_h²+d²)/(4Dt)] e^(−σ_a t).   (6)

---

## 3. The travel-time sensitivity kernel

A trajectory's δT from a localized δβ/β near x: δT = −∫_path (δβ/β)(dℓ/β)  (7). So each unit of *time*
the path spends near x contributes −(δβ/β) to fractional travel-time change. Expected occupation time
(Markov bridge): ρ(x) = [1/P(s,r,t)]∫₀ᵗ P(s,x,u)P(x,r,t−u)du   (8).

**Forward relation + kernel (Pacheco & Snieder 2005):**
δτ(t)/t = −∫_V K(x;s,r,t)·(δβ/β) d³x.   (9)
K(x;s,r,t) = [∫₀ᵗ P(s,x,u)P(x,r,t−u)du] / [∫_V d³x' ∫₀ᵗ P(s,x',u)P(x',r,t−u)du].   (10)

Chapman–Kolmogorov collapses the denominator to t·P(s,r,t):
K = [1/(t·P(s,r,t))]∫₀ᵗ P(s,x,u)P(x,r,t−u)du,  ∫_V K d³x = 1.   (11)
Consequences: (i) homogeneous perturbation → δτ/t = −δβ/β (the stretch); (ii) **intrinsic absorption
cancels** in K's *shape* (Q_i only sets usable lapse range).

**Explicit half-space diffusion kernel** (a = 4D, ρ_s=|x_h|, ρ_r=|x_h−r_h|):
K(x;s,r,t) = (1/N(t)) ∫₀ᵗ du/[u(t−u)]^(3/2) · [e^(−(ρ_s²+(z−d)²)/(au)) + e^(−(ρ_s²+(z+d)²)/(au))] · e^(−(ρ_r²+z²)/(a(t−u))).   (13)
N(t) = e^(−(r_h²+d²)/(at)) / √t.   (14)
**One 1-D numerical integral per model point — no wavefield simulation.**

---

## 4. Single-scattering kernel — the isochron ellipsoid

Early S-coda (Born, Pacheco & Snieder 2006):
K_ss(x) ∝ [e^(−(R_s+R_r)/ℓ)/(R_s²R_r²)]·δ(t−(R_s+R_r)/β),  R_s=|x−s|, R_r=|x−r|.   (15)
δ-function = prolate ellipsoid with foci at source & receiver (plus mirror at image source). Kernel
evolves: thin ellipsoidal shell (early) → two-lobe diffusive cloud (late). For a deep source the
early-coda ellipsoid threads the full crustal column including the fault zone.

---

## 5. Validity caveat

Diffusion kernel (13) is exact only beyond ~one ℓ* from both source and receiver. Within ℓ* of a sensor
the flux is directional → need the specific-intensity (RTE) kernel (Margerin et al. 2016) or Monte-Carlo
RTE. **Favorable for us:** the deep lobe (fault zone, ~30 km) is plausibly near/beyond ℓ* from the
surface receiver → diffusion kernel is most trustworthy *exactly where we want the answer*. Splice an
RTE/specific-intensity kernel within ℓ* of the surface station (the only receiver-side body↔surface
coupling region). Source leg has **no 2-D term** → the only splice is at the receiver.

---

## 6. Geometry & why the deep lobe is special

K is large where *both* legs are large → **deep lobe** at the source (inside the LFE fault zone; absent
in surface–surface geometry), **surface lobe** at the station (factor-2), **connecting tube** (early-coda
ellipsoid = sharp limit).
- **Lapse-time depth migration:** early t weights least-scattered/direct paths → deep tube + fault zone;
  growing t diffuses → shallow near-receiver sensitivity gains weight. Scanning windows t₁<t₂<… gives a
  depth-resolving family of kernels.
- **Common-mode deep lobe (network trick):** all stations share source s → all kernels have the deep
  lobe in the same place → a real fault change is a **coherent δv/v across all station pairs & lapse
  times**; a shallow change near one station appears only in that station's pairs.

---

## 7. Forward problem (two observables)

(a) **Stretch** (velocity channel): (δτ/τ)_{a,i} = −∫_V K(x;s,r_a,t_i)(δβ/β) d³x.   (16)
(b) **Decorrelation** (scatterer-change channel): 1−CC_{a,i} ≈ ∫_V K^dec(x;s,r_a,t_i) δΣ(x) d³x,
   K^dec ∝ [∫₀ᵗ P(s,x,u)P(x,r,t−u)du]/P(s,r,t).   (17)

---

## 8. Velocity change vs. source-location change

- **Velocity change** → coherent mean shift; δτ **grows ∝ lapse time**, same sign across window (stretch).
- **Source displacement δ** → path-length change only on the source leg; mean vanishes, variance nonzero
  and lapse-time-independent: ⟨δτ⟩=0, σ_τ² = δ²/(3β²) (single S mode).   (18)
- **Discriminant:** fit each window for a stretch (coherent, ∝t) AND residual scatter (incoherent, flat).
  Common-mode growing stretch → fault-zone velocity change; flat decorrelation/variance with ⟨δτ⟩≈0 →
  source relocation of magnitude δ = √3·β·σ_τ. They do not alias.

---

## 9. The inverse problem

Discretize half-space into cells V_j. d = G m + ε,  G_{(ai),j} = −K(x_j;s,r_a,t_i)V_j,  m_j = (δβ/β)(x_j).
(19). Tikhonov: m̂ = argmin ‖Gm−d‖²_{C_d⁻¹} + λ²‖Lm‖²  (20). Bayesian/Gaussian solution + resolution R (21).
Solve LSQR/CGLS; inspect R + checkerboard. **Joint source location:** decoupled (stretch→velocity;
variance→δ) or fully joint [G | H] with H = source-displacement partials. **Design:** multiple lapse
windows → depth; azimuthal station coverage → resolves deep lobe; multiple families → 4-D sparse
tomography; Q_i sets usable lapse range, not kernel shape.

---

## 10. Minimal computational recipe

1. Estimate β, ℓ* (→ D = βℓ*/3), Q_i from LFE coda-envelope decay.
2. Per (station a, lapse window t_i): evaluate half-space kernel (13)–(14) on the 3-D grid by 1-D
   quadrature; within ℓ* of a sensor use RTE/specific-intensity (or MC-RTE); use the ellipsoidal
   single-scatter kernel (15) for earliest windows.
3. Measure per pair/window the stretch (δτ/τ) and the residual decorrelation/variance.
4. Build G (19), invert (20)–(21) for fault-zone δβ/β; estimate δ = √3·β·σ_τ from the incoherent channel.
5. Validate against full-wavefield simulations with a planted deep velocity perturbation.

---

## References (key)
Pacheco & Snieder 2005 (JASA), 2006 (GJI); Snieder 2006 (PAGEOPH); Snieder & Vrijlandt 2005 (JGR);
Margerin et al. 2016, 2019 (GJI); Obermann et al. 2013b, 2016 (GJI); Paasschens 1997 (PRE); Barajas et
al. 2022 (GJI); Planès et al. 2014; Rossetto et al. 2011; Sens-Schönfelder & Wegler 2006; Poupinet et
al. 1984; Robinson et al. 2011; Singh et al. 2019. (Full DOIs in the original upload.)
