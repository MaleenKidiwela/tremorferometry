# Cascadia LFE-coda dv/v — End-of-Run Methodology & Results (2026-07-21)

Detailed methodology and coverage/results summary for the full campaign: from the borehole array to the
178-station broadband fleet, through to the deep-resolution characterization. Companion figures:
`fault_tomography/inversion/res_catalog/atlas.png` (sensitivity atlas) and `checkerboard_epochs.png`
(multi-scale checkerboard). Interactive dv/v map: `figures/borehole_dvv_map.html` (164 stations).

---

## 1. Objective

Build the highest-resolution 4-D shear-velocity-change (δβ/β) monitor of the Cascadia plate interface from
LFE coda-wave interferometry, and honestly characterize **how well the deep interface is resolved, and how that
resolution evolves as station coverage grows over 2009–2026.** The novel axis is DEPTH (the deep megathrust),
where the ambient-noise / doublet literature is weak.

---

## 2. Data & instruments — a two-tier network

**Everything processed at 40 Hz** (native-40 channels used directly; higher rates decimated via anti-aliased
`resample_poly`, ratio per-station from `Fraction/gcd`). No station skipped for sample rate.

- **Tier 1 — borehole array (30 stations).** PB-network Cascadia borehole seismometers (B0xx/B9xx), the
  validated backbone. 27 are inversion-grade (≥20 causality-certified families); 2,095 certified families total.
  Top: B011 119, B035 113, B926 110, B943 109, B928 106, B024 100.
- **Tier 2 — broadband/land fleet (178 stations processed).** Region-stratified over the whole margin
  (39.7–50.4°N), high-gain broadband seismometers only (BH?/HH?, 2nd SEED char `H`). **124 INCLUDE / 55 FLAG**
  (~69% include). Networks: UW, CN, CC, NC, BK, others; N. California nets (BK/NC) via NCEDC, rest via IRIS/
  EarthScope. Top: OTR 155, YELM 155, WISH 153, DOORS 146, MARQ 141, LTC 139, JEDS 137, LDH 136.
- **Anchors (3).** PGC, SHB, CLRS — co-located / discovery validation stations (multi-era).

**Certified-family inventory: 12,703 causality-certified ("good") families across 187 stations** (boreholes
2,095 · broadband+anchors 10,608). Certification = coda 2–4 s RMS / pre-arrival mirror RMS > 1.5 (see §3).

---

## 3. Per-station pipeline (the validated recipe)

Driver: `scripts/fleet_station.sh` (fleet), settled borehole recipe in `notes/FINAL_PIPELINE.md`. Stages:

1. **Download** (`download_broadband.py`) — always resumable (never skip on a partial fragment); `.nd` nodata
   sentinels; BK/NC → NCEDC provider.
2. **Band pick** (`pick_band.py`) — one vertical band spanning the most of the 2010–2026 tremor era; requires a
   high-gain seismometer (2nd char `H`, never `HN`=accelerometer); span aggregated PER BAND across epochs.
3. **Detect** — PNSN-tremor-driven candidate detection at 40 Hz (`discover_nllb_pnsn_driven.py`), bbox ±0.9°/
   ±1.35° around the station. Guard: <100 candidates → FLAG (no tremor/data).
4. **Score** — LFE-picker probability per candidate (`score_candidates.py --target-fs 40`, MANDATORY 40 Hz
   resample — the picker was trained at 40 Hz).
5. **Adaptive top-N** (`adaptive_cand_threshold.py`) — rank-based top-30k, floor 0.05 (garbage bound only).
6. **Cluster** (`discover_gpu.py`, GPU) — matched-filter families, cc≥0.80, ≥3 members, **--min-years 1**
   (default 3 is a 4×-yield trap). Guard: <50 families → FLAG (skip densify).
7. **Coverage-select the 300 densify budget by family SNR** (`select_families_coverage.py 'snr'`). SNR is
   predictive of causality survival (AUC 0.68); the family-stack picker P(LFE) is ANTI-predictive (AUC 0.40)
   and was dropped. Densify is a BUDGET (≤300); dv/v uses ALL certified families.
8. **Densify** — forward matched-filter only (`densify_gnw_gpu.py`, 2–8 Hz, threshold 0.8), 40 Hz.
9. **Daily stacks** (`build_long_window_3comp.py`) — per family, per day: normalized average of that day's
   matched detections (cc≥0.80), **≥20 detections/day** or the day is dropped (SNR bar). Durable product:
   `data/long_window_daily_<STA>_Z.npz` (195 survive on disk, 1.3 GB — NEVER deleted by the trace janitor).
10. **dv/v** (`dvv_roll30cal.py`) — 30-calendar-day TRAILING rolling stack, SVD-Wiener filter, coda stretch in
    the **2–4 s** window (off the pinned direct-S; 1–4 s dilutes the stretch ~50× toward zero), `--origin-anchor`
    for true scale, 1-day step, ≥5 daily stacks in the trailing window or no point (gap-honest).
11. **Causality certification** (`finalize_causality.py`) — a family is reliable iff RMS(coda 2–4 s) /
    RMS(mirror −2..0 s) > 1.5. The mirror (pre-arrival, acausal) window doubles as the noise estimate.

**Inclusion gate (frozen 2026-07-13, rule-before-result):** INCLUDE iff ≥20 certified families AND survival
(certified ÷ densified) ≥15%. Else FLAG (kept, completed, banked — not discarded).

---

## 4. Key methodological decisions

- **40 Hz everywhere** — picker trained at 40 Hz; wrong-Nyquist features otherwise.
- **Causality (mirror) is the SOLE LFE verdict** — the coda/pre-arrival ratio certifies reliable families and the
  mirror simultaneously estimates the noise. Reverse-densify retired; the mirror does both jobs.
- **Family selection by SNR, never picker p_lfe** — the family-stack absolute-P gate is anti-predictive.
- **Survival = certified ÷ DENSIFIED**, not clustered (a denominator bug had wrongly failed BBO at 116 certified).
- **Rolling-buffer storage** — raw traces freed after each station DONE; daily-stack npz + dv/v CSVs + cert files
  are the durable products and are never deleted.
- **Merlin (Fable-5 advisor) consulted before every method/parameter decision** — caught the min-years trap, the
  survival denominator, the +120 s shift-null requirement, and (below) the entire resolution-analysis framing.

**The 8 pipeline bugs caught & fixed during the run** (regression tells in `notes/FLEET_BROADBAND_PIPELINE.md`
§5): (1) selection by p_lfe → SNR; (2) survival denominator; (3) download skipping fragments; (4) BK/NC not at
IRIS → NCEDC; (5) pick_band scored per-epoch → aggregate per band; (6) accelerometer HN mis-pick → require 2nd
char H (KSXB recovered 2→78 certified); (7) GPU not serialized; (8) appB family-stack picker gate dropped.

---

## 5. Coverage & FLAG taxonomy

**INCLUDE (map-grade): 164 stations** on the interactive map (≥20 certified). FLAGs are kept, completed, and
banked — flagged on family COUNT/SURVIVAL, never on family RELIABILITY, so their certified families remain valid.

FLAG diagnosis precedence (before calling any FLAG "genuine weak"): band (accelerometer?) → provider (NCEDC?) →
download completeness → only then genuine-weak. FLAG classes and the banked lists:

- **CODA-LESS (new sub-type, confirmed 3×: LRIV, PR01, PR03).** Families cluster and produce dv/v, but the
  causality distribution has NO high tail — `max(causality) < 1.4` across ALL families (vs 3–13 for real
  stations). Emergent/microseism-dominated sites (N. Olympic coast; Cascade-volcano "PR" array near Rainier/
  St. Helens). Diagnostic is the tail SHAPE, not detection/family count (KSXB has 10× LRIV's detections and
  certifies 78). Per-station, NOT per-array (PR02 INCLUDEd at 131).
- **GENUINE-WEAK / <50 families** (PCOR, WEAV, DOSE, TKEY, UNFR, KAUT, …) — real data, too few coherent repeaters.
- **DATA-SPARSE at IRIS** — per-station, NOT a network rule: most CN is well-served (TOFB 89, CBB 84, SNB 115,
  MGRB 117, TXDB 127); only individual sparse stations need NRCan → `data/cn_sparse_nrcan.txt` (MYRA, SHPB, BOIB).
- **SHORT-RECORD** (ONRC: 3 day-files at IRIS) — recent installs, too little record for multi-year dv/v.
- **TIER-2 (count-OK / survival-short, ≥20 cert but <15% survival)** → `data/tier2_survival_borderline.txt`:
  BHW 21@7%, KRP 38@13%, BHAM 33@14%, MOGU 21@11%, PIT6 26@9%. Reliable families, noisy pool. Inversion-inclusion
  decided later on evidence.
- **COUNT-SHORT / SURVIVAL-OK (<20 cert but survival ≥15%)** → `data/count_short_survival_ok.txt`: TRIPT 17@20%,
  CRBN 14@23%. Family-supply-limited; the mirror of tier-2. Also inversion-inclusion candidates.
- **BUG-FLAG re-run queue** → `data/rerun_queue.txt` (6): FISH, JCC, GOBB, SYMB, KMR, KHMB. (JCC re-flagged 0
  candidates even via NCEDC — needs a closer look.)

---

## 6. Validation (before the fleet was trusted)

- **B011 ↔ PGC (borehole vs co-located broadband, ~290 m apart).** Mirror-null method (raw correlation is invalid
  at 290 m — shared noise field guarantees high correlation even for pure noise; Merlin). Both instrument eras
  PASS; per-family matched-vs-mismatched p=0.0005. Broadband reproduces the borehole. `scripts/b011_pgc_validation.py`.
- **CLRS discovery path (no borehole answer key).** Two-arm + adaptive top-N; +120 s shift-null confirmed the
  signal is LFE identity, not tremor intensity. 62 certified families = the same LFEs a co-located borehole sees.
  → the picker-discovery path is valid where no borehole exists. This is what licenses the fleet.

---

## 7. Deep-resolution characterization (the resolution-through-time study)

**Goal:** how well is the DEEP interface resolved, and how does that improve as coverage grows? Merlin-vetted
design; NO velocity inversion is run — only resolution/information characterization. Scripts:
`assemble_res_catalog.py` → `build_G_captured.py` → `sensitivity_atlas.py` (+ `gates_AB.py`,
`checkerboard_epochs.py`).

### 7.1 Model & forward operator
- Model: δβ/β on 0.10° interface cells at their Slab2 depth (872 cells, 456 deep >30 km). Combined network:
  ALL 12,703 certified families from ALL 187 stations → **6,726 (cell, station) data** (families collapse to
  cells) over **173 distinct sites** (co-located stations clustered <2 km).
- Kernel: single-scattering early-coda ellipsoid (`kernels/kernel.py::kernel_singlescatter`, 2–4 s window,
  β=3.5, ℓ*=40 km, free-surface image). δτ/τ = −∫K·δβ/β (rows sign-flipped).

### 7.2 The captured-fraction correction (retiring "0.77")
The kernel is VOLUMETRIC; an interface (a 2-D surface) captures only a small share of a 3-D coda shell. Measured
**on-interface captured fraction: median 0.047** (±1.5 km slab; grid-convergence PASS, base 0.051 vs refined
0.050). The historic checkerboard used the unit-sum kernel, which renormalizes every ray to full interface
sensitivity — inflating interface amplitudes ~12× and making weak rays look fully informative. **Fix (Merlin,
adopted as the production operator, not a test-only patch): multiply each row by its captured fraction f**
(calibrated operator Gc = unit-sum × f). This is simultaneously the physical calibration and the F5 down-
weighting. **The earlier 2-D INTERFACE 0.77 is retired** as optimistic (stale 28-station/1–4 s assembly, no
causality/cc filter, no site terms, unit-sum inflation). NOTE: the 3-D VOLUME 0.86 is NOT retired — the volume
model captures the full sensitivity so it reproduces honestly on the new network (0.85 @ 250 km, §7.5b); the
inflation was specific to projecting a volumetric kernel onto a thin interface.

### 7.3 Two gates before any verdict
The initial per-month resolution matrix showed χ²/N flooring at ~8 (flat in λ) — the interface under-fits.
Merlin flagged this could be a removable artifact. Two gates settled it:
- **Gate A (mirror confound).** On the 7 stations with a mirror product, raw vs mirror-corrected post-site
  residual: median 0.331% → 0.316% (~5% change). **The family-specific floor is REAL, not the removable mirror
  artifact** (the mirror noise is station-common, already absorbed by the site terms).
- **Gate B (statics).** Post-site residual is **99% time-varying, 1% static** (so epoch-difference maps get no
  free pass), lag-1 month autocorr 0.60 → only ~3 independent samples/yr (annual stacking gain ×1.7, not ×√12).

**Verdict (experiment-backed, CORRECTED — see §7.5 for the correction).** After per-station SITE terms absorb the
shallow common-mode (~53% of variance), the family-specific residual is ~0.25–0.34% at MONTHLY cadence, so a
single month does not give a clean cell map. BUT — an initial reading of this as "a margin-wide INDEX, not a map"
was OVER-PESSIMISTIC: it came from a λ-selection bug (the atlas PSF-floor λ wrongly transplanted into the
checkerboard → over-smoothing). Corrected with oracle λ (§7.5): the deep interface is a RESOLVABLE MAP at
~70–140 km scale, sharpening from ~180 km as the fleet grew, with the geometry itself resolving to the 44 km cell
scale (the limit is SNR, not geometry). The "index/bound" statement holds only at the strictest DETECTABILITY
tier (real signals at monthly cadence); it is NOT the resolution limit.

### 7.4 The sensitivity/information atlas through time
Operator: calibrated fault + per-station site terms; λ_f frozen by the PSF-floor criterion (smallest λ s.t.
median deep PSF ≥ 30 km = kernel footprint; χ²/N is λ-insensitive because unmodeled variance dominates, so at
fine scales the reported resolution is regularization-set BY CONSTRUCTION — the data-driven content is in the
noise/precision products). Realistic per-datum σ = time-varying residual (0.25%); measurement-floor σ = 0.13%.

- **(a) Deep-index precision vs year** — posterior σ of the deep-cell area-weighted mean δβ/β: **0.64% (2009,
  13 sites) → 0.12% (2021+, ~145 sites)**, ~5× improvement, plateauing once the fleet saturated coverage.
- **(b) Resolved independent deep modes (SVD of the whitened, site-projected deep operator, sv>1)** — **~3
  (2010–2020, borehole-dominated) → ~15–25 (2021+, fleet era)**, growing with per-month active coverage (~85 →
  ~145 sites), fluctuating month-to-month with ETS-driven family activity. These modes ARE the large-scale
  patterns (margin-wide mean, N–S gradient, …).
- **(c) Per-cell 2σ minimum-detectable-amplitude maps** per epoch — the map-shaped, signal-agnostic product
  ("in year T the network could have detected a deep transient of X% at this cell/region").

### 7.5 Multi-scale interface checkerboard — CORRECTED (oracle λ, three noise tiers)
Merlin audited the checkerboard code and found the atlas PSF-floor λ=2.069 was wrongly reused in the checkerboard
(different whitening) → over-smoothing → artificially low recovery (the earlier "0.06/0.69" and the frozen-λ table
below it). CORRECTED (`checkerboard_corrected.py`): per (epoch, scale, noise tier), sweep λ and report the best
achievable (a resolution test reports the estimator class's ceiling; the frozen production λ lives ONLY in the
atlas). Deep well-covered cells. Finest resolved scale (corr>0.7):

| epoch | sites | noise-free (geometry) | σ=0.13% (measurement) | σ=0.25% (total residual) |
|-------|-------|-----------------------|-----------------------|--------------------------|
| 2011  | 51    | 44 km (corr ~0.99)    | 180 km                | none (250 km=0.61)       |
| 2016  | 62    | 44 km                 | 100 km                | none                     |
| 2021  | 136   | 44 km                 | **70 km**             | 180 km                   |
| 2025  | 149   | 44 km                 | **70 km**             | 180 km                   |
| all   | 188   | 44 km (~1.00)         | **70 km** (44 km=0.68)| 140 km                   |

KEY: (i) **the geometry resolves 44 km fully (noise-free ~1.0 every epoch) — no geometric obstruction**; the
collapse is 100% SNR. (ii) The fleet SHARPENED resolution from ~180 km (2011) to ~70 km (2021+) at measurement
noise. (iii) A ±1% interface checker produces ~0.035–0.05% data-space signal — matching the observed +0.043% ETS
transient, confirming the calibration reproduces the real signal class. (iv) Noise reduction (fleet-wide mirror
correction, longer averaging) buys resolution directly toward 44 km — the roadmap. The old frozen-λ numbers (35 km
0.05 / 250 km 0.72) are SUPERSEDED.

### 7.5b INTERFACE vs VOLUME parameterization — which "0.77/0.86" survives
Comparison on the SAME new network (deep well-covered cells), OLD testing regime vs the HONEST corrected one
(oracle λ, §7.5):

| 2-D interface, deep | OLD (unit-sum, no-site, inv-crime noise, λ=0.3) | HONEST (calibrated + site + real noise, oracle λ) |
|---|---|---|
| 44 km  | 0.92 | ~0.43 (σ=0.25%) / ~0.66 (σ=0.13%) |
| 250 km | 0.96 | ~0.73 (σ=0.25%) / ~0.84 (σ=0.13%) |

Under OLD assumptions the new (bigger) network recovers 0.92/0.96 — HIGHER than the old 0.77/0.86, so the network
genuinely improved. The honest recovery is lower mainly because of **capture-weighting** (a thin interface sees
~5% of the volumetric coda sensitivity) plus realistic noise — but it is NOT a collapse (0.43–0.66 at 44 km,
0.73–0.84 at 250 km). (An earlier ablation reported 0.01/0.06 at the honest tier; those came from the same
frozen-λ over-smoothing bug and are SUPERSEDED — the buggy `checkerboard_ablation.csv` was removed.) **Only the
old 2-D INTERFACE 0.77 is retired as inflated.**

The **3-D VOLUME checkerboard** on the new network (`volume3d_checkerboard_new.py`; normalized kernel over the
crustal volume, site terms, realistic noise) — where the sensitivity actually lives, so no capture dilution —
**RECOVERS WELL and reproduces the old ~0.86 honestly:** overall 0.70 (70 km) / 0.81 (140 km) / **0.85 (250 km)**,
holding at depth (32–48 km: 0.63–0.85). Depth layers trade off (single 2–4 s window → weak vertical separation),
so "by-depth" curves track each other. **The old VOLUME finding stands; it was the interface projection that was
ill-posed.**

**Synthesis (the honest resolution statement, CORRECTED).** The coda sensitivity is volumetric. Broad/VOLUMETRIC
deep δβ/β is resolvable at large scale (~140–250 km, incl. at depth; corr 0.81–0.85). For the INTERFACE
parameterization (the harder, more specific hypothesis): the network GEOMETRY resolves down to the 44 km cell
scale (noise-free corr ~1.0) — there is NO geometric obstruction; resolution is SNR-limited. At measurement noise
the finest resolved interface scale is ~70 km (recent), sharpened from ~180 km by the fleet; at the total
family-specific residual (0.25%, the detectability tier) it is ~140–180 km. So the deep interface is a RESOLVABLE
MAP at ~70–140 km, improving with coverage, and heading toward 44 km as the noise floor drops (fleet-wide mirror
correction, longer averaging). Consistent with the atlas: mode count 3→24 and index precision 0.64%→0.12% (both
λ-independent, stand). The "margin-wide index/bound" reading applies ONLY at the strictest single-month
detectability tier — it is a detectability statement, NOT the resolution limit. (Correction history: an initial
frozen-λ checkerboard bug depressed interface recovery and led to a premature "index, not map" verdict; Merlin's
audit found the λ transplant, and the oracle-λ rerun gives the numbers above.)

**Caveats stated (Merlin ruling 4c wording):** the coda kernel is volumetric; an interface parameterization
captures a median ~5% of shell mass per datum. The operator is calibrated for this (rows × f), so amplitudes are
physically scaled under the interface-confinement assumption; the near-receiver 7–14 km shell is family-
independent and absorbed by site terms; mid-crustal volumetric changes are not modeled and contaminate fault
cells (bounded by the leakage probe — future work). The earlier 0.77 checkerboard is superseded.

---

## 7.6 The 4-D inversion (v1, raw tensor) — `invert_dvv_4d.py`
Merlin-approved recipe: per-pair harmonic DESEASON → per-pair DEMEAN (anomaly maps; absolute level
unidentifiable) → stack. Products: (1) annual anomaly movie 2010–2026 (`inversion_annual.png`); (2) ETS-phase
COMPOSITE difference (ETS months − inter-ETS, tremor-catalog tercile split) as the primary science figure
(`inversion_ets_index.png`); (3) resolved-mode time series (from the atlas). λ RE-PICKED by a 100 km checkerboard
oracle → **λ_f=0.433** (recovery 0.93), NOT the over-smoothing 2.069. Site terms per station, 0.10° grid, whitened.

**Honesty gates (all PASS):**
- Recent deep index data-space **0.0039%** < the 0.02% ETS-null network bound.
- **Scrambled-year NULL**: max deep |anomaly| 0.35% vs real annual max 0.17% → **every annual anomaly is WITHIN
  the null** → no resolvable secular deep change.
- ETS-composite **closure**: forward-predicted data std 0.012% < observed 0.043% (not inflating noise).
- ETS-composite **fault variance-reduction beyond site terms: 2.3%** (marginal real fault info).
- **~129/456 deep cells** resolved to 2σ_m ≤ 0.5% (recent).

**Result:** the deep megathrust interface is **velocity-stable within resolution — a null-gated, now
spatially-resolved BOUND** (~28% of the deep interface constrained), consistent with the ETS-null history.
Caveats: v1 on the raw tensor (mirror-correction v2 would lower the noise floor and sharpen resolution — the
geometry already resolves 44 km); annual maps contain NO ETS transients by construction; amplitudes are
model-space (~20× data-space); anomaly maps only (no absolute level). Files: `inversion_4d.npz`,
`inversion_annual.png`, `inversion_ets_index.png`.

## 8. Bottom line

- **157+ stations with certified deep dv/v** (30 boreholes + 124 broadband + 3 anchors), **12,703 certified
  families** — the largest LFE-coda dv/v network assembled on Cascadia, validated borehole↔broadband and via the
  no-borehole discovery path.
- **The deep megathrust dv/v resolves as ~large-scale (≥~200 km) regional structure + a margin-wide index**, not
  a fine 4-D map. The broadband fleet took the deep-index precision from 0.64% → 0.12% (~5×) and the number of
  above-noise deep modes from ~3 → ~24, but the robustly-resolved scale stayed ~200 km. This is the honest,
  quantified answer to "spatial resolution of the deep dv/v through time," and it reinforces the deep-megathrust-
  is-velocity-stable (bound) result rather than claiming a deep 4-D signal.

## 9. Future work
- Fleet-wide mirror pass (stacks survive on disk) → cleaner family-specific noise; re-test whether annual-cadence
  regional difference maps become viable.
- Two-depth shallow-leakage probe (shallow 3–12 km vs mid-crust 15–30 km) to quantify off-interface contamination.
- Drain the re-run queue (6) and NRCan retry for CN-sparse (3).
- Decide tier-2 / count-short inversion inclusion on robustness evidence.
- The 4-D δβ/β inversion itself (boreholes + broadband), reported at the resolution this study establishes:
  regional/large-scale, with per-epoch detectability from the atlas.
