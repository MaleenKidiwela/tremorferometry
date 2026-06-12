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

## NULL TEST: Is the southern "hot zone" temporal RMS real or a coverage artifact? (2026-06-10)

**Question:** The real multi-window inversion shows elevated temporal RMS in the sparse south
(42–44°N: median 0.09–0.15%, 1–2 stations; 40–42°N: median 0.29%, 1–2 stations).  Does this
reflect genuine fault velocity variability, or is it noise amplified through a sparse design matrix?

**Method:**  `fault_tomography/inversion/null_test.py`
- Planted fault field: ZERO everywhere, all months.
- Forward-modelled every real (cell, station, window, month) datum using the real kernels and
  real coverage pattern.
- Added synthetic noise drawn per datum from the REAL per-(station, window) residual std, estimated
  by pooling residuals from the real inversion across all solved months
  (noise_std range: 0.020–0.909%, median 0.098%).
- Inverted with identical joint LSQR + Laplacian + site-term machinery (λ_f=0.4, λ_s=0.05).
- Repeated 25 realisations (different RNG seeds).

**Noise std used (per station|window pair):**
- Estimated from real inversion residuals (dd − Ĝf mf − Ĝs site), pooled across all ok months.
- Min: 0.020%, Median: 0.098%, Max: 0.909%.

**Results — null temporal RMS vs observed, by latitude band and station coverage:**

| Lat band | nsta | n cells | obs median (%) | null 95th pct (%) | null max (%) | obs pctile in null |
|---|---|---|---|---|---|---|
| 40-42 | 1 | 67 | 0.294 | 0.135 | 0.152 | 100th |
| 40-42 | 2 | 28 | 0.289 | 0.133 | 0.146 | 100th |
| 40-42 | 3 | 11 | 0.289 | 0.132 | 0.140 | 100th |
| 42-44 | 1 | 59 | 0.131 | 0.105 | 0.131 | 100th |
| 42-44 | 2 | 43 | 0.145 | 0.097 | 0.118 | 100th |
| 42-44 | 3 | 25 | 0.093 | 0.088 | 0.111 | 100th |
| 42-44 | >=4 | 6 | 0.052 | 0.053 | 0.060 | 100th |
| 44-46 | 1 | 47 | 0.050 | 0.079 | 0.099 | 100th |
| 44-46 | 2 | 22 | 0.057 | 0.073 | 0.083 | 82th |
| 44-46 | 3 | 19 | 0.050 | 0.067 | 0.074 | 56th |
| 44-46 | >=4 | 22 | 0.049 | 0.058 | 0.070 | 100th |
| 46-48 | 1 | 39 | 0.030 | 0.041 | 0.047 | 64th |
| 46-48 | 2 | 34 | 0.033 | 0.037 | 0.045 | 100th |
| 46-48 | 3 | 25 | 0.034 | 0.038 | 0.040 | 92th |
| 46-48 | >=4 | 32 | 0.028 | 0.035 | 0.039 | 52th |
| 48-50 | 1 | 74 | 0.053 | 0.094 | 0.118 | 100th |
| 48-50 | 2 | 48 | 0.052 | 0.073 | 0.115 | 78th |
| 48-50 | 3 | 38 | 0.048 | 0.065 | 0.088 | 100th |
| 48-50 | >=4 | 18 | 0.044 | 0.051 | 0.057 | 100th |

**Verdict:**

- **40-42 1-sta (n=67):** obs median 0.294% sits at 100th pctile of the null — SIGNIFICANT (exceeds null 95th = 0.135%). Likely real signal or unmodelled systematic.
- **40-42 2-sta (n=28):** obs median 0.289% sits at 100th pctile of the null — SIGNIFICANT (exceeds null 95th = 0.133%). Likely real signal or unmodelled systematic.
- **42-44 1-sta (n=59):** obs median 0.131% sits at 100th pctile of the null — SIGNIFICANT (exceeds null 95th = 0.105%). Likely real signal or unmodelled systematic.
- **42-44 2-sta (n=43):** obs median 0.145% sits at 100th pctile of the null — SIGNIFICANT (exceeds null 95th = 0.097%). Likely real signal or unmodelled systematic.

**Overall conclusion:**
The null test directly answers the open question from the synthetic resolution test (§4, 2026-06-09).
Cells where the observed temporal RMS lies within the null distribution (obs_pctile ≤ 95th) are
consistent with pure noise amplification through the sparse design matrix and should NOT be
interpreted as real fault velocity change. Cells/classes where the observed RMS significantly
exceeds the null distribution warrant further investigation.

Figure: `figures/null_test_southern.png`
Result archive: `fault_tomography/inversion/null_test_results.npz`

## NULL TEST: Is the southern "hot zone" temporal RMS real or a coverage artifact? (2026-06-10)

**Question:** The real multi-window inversion shows elevated temporal RMS in the sparse south
(42–44°N: median 0.09–0.15%, 1–2 stations; 40–42°N: median 0.29%, 1–2 stations).  Does this
reflect genuine fault velocity variability, or is it noise amplified through a sparse design matrix?

**Method:**  `fault_tomography/inversion/null_test.py`
- Planted fault field: ZERO everywhere, all months.
- Forward-modelled every real (cell, station, window, month) datum using the real kernels and
  real coverage pattern.
- Added synthetic noise drawn per datum from the REAL per-(station, window) residual std, estimated
  by pooling residuals from the real inversion across all solved months
  (noise_std range: 0.036–0.700%, median 0.193%).
- Inverted with identical joint LSQR + Laplacian + site-term machinery (λ_f=0.4, λ_s=0.05).
- Repeated 25 realisations (different RNG seeds).

**Noise std used (per station|window pair):**
- Estimated from real inversion residuals (dd − Ĝf mf − Ĝs site), pooled across all ok months.
- Min: 0.036%, Median: 0.193%, Max: 0.700%.

**Results — null temporal RMS vs observed, by latitude band and station coverage:**

| Lat band | nsta | n cells | obs median (%) | null 95th pct (%) | null max (%) | obs pctile in null |
|---|---|---|---|---|---|---|
| 40-42 | 1 | 57 | 0.190 | 0.290 | 0.323 | 40th |
| 40-42 | 2 | 34 | 0.198 | 0.289 | 0.324 | 48th |
| 40-42 | 3 | 20 | 0.191 | 0.288 | 0.308 | 48th |
| 40-42 | >=4 | 14 | 0.194 | 0.290 | 0.304 | 44th |
| 42-44 | 1 | 59 | 0.132 | 0.218 | 0.265 | 40th |
| 42-44 | 2 | 43 | 0.141 | 0.218 | 0.255 | 48th |
| 42-44 | 3 | 25 | 0.113 | 0.198 | 0.247 | 40th |
| 42-44 | >=4 | 6 | 0.095 | 0.131 | 0.164 | 60th |
| 44-46 | 1 | 47 | 0.063 | 0.078 | 0.104 | 100th |
| 44-46 | 2 | 22 | 0.069 | 0.069 | 0.091 | 98th |
| 44-46 | 3 | 19 | 0.059 | 0.090 | 0.112 | 64th |
| 44-46 | >=4 | 22 | 0.058 | 0.069 | 0.091 | 98th |
| 46-48 | 1 | 37 | 0.054 | 0.066 | 0.075 | 100th |
| 46-48 | 2 | 25 | 0.055 | 0.065 | 0.074 | 100th |
| 46-48 | 3 | 21 | 0.051 | 0.060 | 0.069 | 92th |
| 46-48 | >=4 | 60 | 0.055 | 0.059 | 0.068 | 96th |
| 48-50 | 1 | 60 | 0.097 | 0.160 | 0.196 | 100th |
| 48-50 | 2 | 49 | 0.088 | 0.131 | 0.191 | 100th |
| 48-50 | 3 | 34 | 0.082 | 0.121 | 0.158 | 100th |
| 48-50 | >=4 | 46 | 0.069 | 0.083 | 0.097 | 100th |

**Verdict:**

- **40-42 1-sta (n=57):** obs median 0.190% sits at 40th pctile of the null (null 95th = 0.290%, null max = 0.323%). CANNOT BE DISTINGUISHED from noise amplification — likely a COVERAGE ARTIFACT.
- **40-42 2-sta (n=34):** obs median 0.198% sits at 48th pctile of the null (null 95th = 0.289%, null max = 0.324%). CANNOT BE DISTINGUISHED from noise amplification — likely a COVERAGE ARTIFACT.
- **42-44 1-sta (n=59):** obs median 0.132% sits at 40th pctile of the null (null 95th = 0.218%, null max = 0.265%). CANNOT BE DISTINGUISHED from noise amplification — likely a COVERAGE ARTIFACT.
- **42-44 2-sta (n=43):** obs median 0.141% sits at 48th pctile of the null (null 95th = 0.218%, null max = 0.255%). CANNOT BE DISTINGUISHED from noise amplification — likely a COVERAGE ARTIFACT.

**Overall conclusion:**
The null test directly answers the open question from the synthetic resolution test (§4, 2026-06-09).
Cells where the observed temporal RMS lies within the null distribution (obs_pctile ≤ 95th) are
consistent with pure noise amplification through the sparse design matrix and should NOT be
interpreted as real fault velocity change. Cells/classes where the observed RMS significantly
exceeds the null distribution warrant further investigation.

Figure: `figures/null_test_southern.png`
Result archive: `fault_tomography/inversion/null_test_results_calT35.npz`
