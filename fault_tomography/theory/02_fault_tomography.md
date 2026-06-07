# Mapping δβ(x₀,y₀) on the Fault Plane from Many LFE-Family / Station Coda Pairs

**Companion to `01_kernel.md`.** That doc gives the volumetric kernel K(x; s,r,t); this one assembles
all (LFE family, station, lapse window) observations into a single inverse problem whose model is the
relative shear-velocity change *as a function of position on the fault plane*, m(ξ) = δβ/β(x₀,y₀).

> Persisted from upload (Note2). Operational summary in `../PLAN.md`.

---

## 0. The idea in one paragraph

Each LFE family is a fixed repeating source *on the fault* at fault coordinate ξ_a = (x₀a, y₀a). The
travel-time kernel for that family at any surface station has its **deep lobe pinned at ξ_a** — strongest
sensitivity to the velocity change of the family's *own* fault patch. So family a "self-images" patch ξ_a,
with weaker leakage to neighbors (deep-lobe width) plus off-fault sensitivity along the up-going tube and
under the station. Tile the fault with the families you have, add every station and lapse window as a
redundant look → fault-plane tomography with a near-diagonal core (family ↔ its patch) + a controllable
off-fault nuisance space that network redundancy suppresses.

---

## 1. Indices, model, the one modeling choice

a=1..A family (source at X_a=X(ξ_a), fault coord ξ_a); b=1..B station at R_b; i=1..I lapse window t_i;
ξ=(x₀,y₀) on the dipping fault surface Σ (embedding X(ξ)); k=1..K fault-mesh node (family locations ⊂ nodes);
h(ξ) fault-zone thickness; A_k patch area; **m_k = (δβ/β)(ξ_k) = target**.

**Modeling choice:** velocity changes are *localized to the fault zone* → the volumetric field collapses
onto Σ. Justified physically (fault zones concentrate fluids/pore-pressure/damage that move β) and
geometrically (deep lobe pinned to Σ by the on-fault source). **Cost:** off-fault changes (shallow site
hydrology, bulk crust) project through the tube + surface lobe → carried as nuisance (§4) or they bias m̂.

---

## 2. 3-D kernel → fault-plane kernel

**Thin-zone collapse:** if δβ/β is confined to thickness h(ξ) about Σ,
δτ/t |_{a,b,i} = −∫_Σ κ_{abi}(ξ)·(δβ/β)(ξ) dA + e^off_{abi},  with **κ_{abi}(ξ) = K(X(ξ); X_a, R_b, t_i)·h(ξ)** [area⁻¹].   (1)

**Captured-sensitivity fraction:** η_{abi} ≡ ∫_Σ κ_{abi}(ξ) dA ≤ 1   (2) = fraction of the coda's
sensitivity actually on the fault for that pair. Largest at short lapse times (tight deep lobe), shrinks
as the field diffuses off-fault and toward the station. **Report η_{abi}** — it's the honest measure of
how fault-confined each datum is; low-η pairs are off-fault-dominated → down-weight them.

Deep-lobe peak width (fault-plane resolution length) ≈ √(2D·t_i) in the diffusive regime (or ~ℓ*).

---

## 3. Discretized forward operator over all pairs

d_{abi} = (δτ/τ)_{abi} = Σ_k G_{(abi),k} m_k + e^off_{abi} + ε_{abi},
**G_{(abi),k} = −κ_{abi}(ξ_k)·A_k = −K(X_k; X_a, R_b, t_i)·h_k·A_k.**   (3)
Stack all valid triples → d = G m + e^off + ε   (4).

**Near-diagonal core:** because κ_{abi} peaks at ξ_a, rows of family a load most onto node k(a)
(co-located with the source). Ordering family nodes first → G ≈ (block-)diagonal self-imaging part +
neighbor leakage:
d_{abi} ≈ −κ^self_{abi} A_{k(a)} m_{k(a)}  −Σ_{k∈N(a)} κ_{abi}(ξ_k) A_k m_k  + e^off + ε.   (5)
→ well-conditioned where families are dense: every family-hosting patch is directly illuminated by its
own source.

---

## 4. Off-fault nuisance & the bipartite (two-way) structure

e^off_{abi} = a **site term c_b** (shallow change near station b, surface lobe, depends on b not source) +
a **path/bulk term p_{abi}** (smooth off-fault crust along the a→b tube).
**Augmented model: d = G_f m_f + G_s m_s + G_p m_p + ε**   (6), m_f=fault map (target),
m_s={c_b} per-station site, m_p=coarse 3-D off-fault field.

**Why the network makes it identifiable:** to leading order d_{abi} ≈ α_{abi} m_{k(a)} + β_{abi} c_b + path
(7) — additive in a source-side (fault, function of a) + receiver-side (site, function of b) effect: a
**bipartite/two-way structure** (row effect per fault patch, column effect per station).
1. **Common-mode over stations isolates the fault:** ⟨d⟩_a = avg over (b,i) ≈ ⟨α⟩_a m_{k(a)} + ⟨βc⟩_a; if
   site changes are zero-mean across the network, m̂_{k(a)} ≈ −⟨d⟩_a/(⟨κ^self⟩_a A_{k(a)}). A *true* fault
   change is the part of family a's signal **coherent across all its stations**.
2. **Common-mode over families isolates the site** (average over all families at station b → c_b).
> If site changes are not zero-mean (regional shallow transient), the marginal is biased → keep m_s explicit.

---

## 5. The overall inverse problem

min over [m_f, m_s, m_p] of:
‖W^{1/2}(G_f m_f + G_s m_s + G_p m_p − d)‖² + λ_f‖L_Σ m_f‖² + λ_s‖m_s‖² + λ_p‖L_3 m_p‖²   (8)
- W = C_d⁻¹ data weights, **scaled by η_{abi}** (eq 2) so low-fault-sensitivity pairs count less.
- **L_Σ = surface Laplacian on the fault mesh, optionally anisotropic** (smoother along strike than dip) —
  couples family-less patches to neighbors so m_f is defined on the full mesh.
- λ_s damps site terms; L_3 smooths the path field.

Bayesian/Gaussian form with C_m = blkdiag(C_f,C_s,C_p), G = [G_f G_s G_p]:
m̂ = (GᵀC_d⁻¹G + C_m⁻¹)⁻¹ GᵀC_d⁻¹ d,  Ĉ_m̂ = (GᵀC_d⁻¹G + C_m⁻¹)⁻¹,  R = Ĉ_m̂ GᵀC_d⁻¹G.   (9)
Solve LSQR/CGLS. Fault block of R → spatial resolution; cross-blocks of Ĉ → fault↔site↔path trade-off.

---

## 6. Resolution & experiment design on the fault

- **Resolution length:** L_res(ξ) ≈ max( √(2D·t_min)  [deep-lobe radius],  Δ_LFE(ξ)  [family spacing] )   (10).
  Two levers: denser families tighten Δ_LFE; **earlier lapse windows tighten the deep lobe** — OPPOSITE of
  surface–surface monitoring (source already at depth → earliest usable S-coda = sharpest fault image). Q_i
  only sets how early-to-late before coda hits noise.
- **Azimuthal station coverage per family** decides whether the on-fault patch separates from the up-going
  tube. One station can't tell "change at ξ_a" from "change anywhere on my tube"; an azimuthal ring pins it.
- **Point-spread / checkerboard tests** through the real G → where the layout resolves the fault vs where
  off-fault leakage (null space) dominates.

---

## 7. Time-lapse (4-D) & the source-migration channel

**4-D map:** every observable is a time series → d_{abi}(T) = Σ_k G_{(abi),k} m_k(T) + …, m_k(T)=(δβ/β)(ξ_k,T).
Invert each epoch (or jointly with temporal smoothing ‖∂_T m_f‖²) → **δβ/β(x₀,y₀,T) across the slab — a
movie of shear-velocity change on the plate interface through an ETS cycle.** ← THE DELIVERABLE.

**Source migration (second-moment channel):** a family displacement δ_a → no coherent stretch but variance
σ²_{τ,a} = δ_a²/(3β²), common across its stations. Estimate per-family δ_a from the incoherent (variance)
channel; orthogonal to m_f (first moment ∝t vs second moment flat) → no aliasing. Gives a map of **where
the LFE sources themselves migrate** — a second fault-process observable for free.

---

## 8. Algorithm

1. **Geometry & medium:** place families on fault mesh X(ξ_a); estimate β, ℓ* (D=βℓ*/3), Q_i, thickness h(ξ).
2. **Kernels:** per (a,b,i) compute volumetric kernel (01_kernel eq 13), restrict to Σ, ×h_k → κ_{abi}(ξ_k);
   record η_{abi}. Ellipsoid/RTE kernel for earliest windows & within ℓ* of stations.
3. **Data:** per (a,b,i,T) measure coherent stretch (δτ/τ) → d, and residual variance σ_τ → migration channel.
4. **Assemble** G = [G_f G_s G_p] (eq 3,6); weight C_d⁻¹ scaled by η_{abi}.
5. **Quick look:** network-common marginals (§4) → regularization-free first map + fault-vs-site sanity check.
6. **Invert** (8)/(9) with fault-surface Laplacian; tune λ_f,λ_s,λ_p by L-curve / cross-validation.
7. **Resolve & validate:** resolution matrix, point-spread, full-wavefield synthetic with a planted on-fault patch.
8. **Iterate over epochs T** → 4-D map; co-estimate per-family migration δ_a from the variance channel.

---

## References (key)
Pacheco & Snieder 2005; Snieder 2006; Snieder & Vrijlandt 2005; Margerin et al. 2016; Obermann et al. 2016;
Brenguier et al. 2008 (Parkfield postseismic); Aster, Borchers & Thurber 2018; Tarantola 2005. (DOIs in original.)
