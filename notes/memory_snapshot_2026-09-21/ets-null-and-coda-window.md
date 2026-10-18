---
name: ets-null-and-coda-window
description: ETS produces no station-consistent deep-interface dv/v change (<~0.003%); LFE coda usable window is 1-4s only (no late coda)
metadata: 
  node_type: memory
  type: project
  originSessionId: 400c26f9-17ff-46fc-914c-9c6a7671906d
---

**⚠ UPDATE 2026-06-10b (southern hot RESOLVED):** the "southern hot zone" (40-42N, flagged 2.2x-null at 28 sta)
was a SINGLE-STATION ARTIFACT (B935 site-term leakage). Adding the Mendocino subarray (B045/B932/B049) → 35-station
true-scale null test: 40-42N observed RMS now sits at the 40-48th percentile of the pure-noise null (consistent
with noise, NOT real). Interface is velocity-stable MARGIN-WIDE; detection floor is coverage-set (~0.06% north,
~0.29% south). Lesson: denser coverage KILLS single-station artifacts. (notes/2026-06-10 §13; null_test_calT35.)

**⚠ UPDATE 2026-06-10:** finding #2 (stable interface) SURVIVES the true-scale fix [[dvv-temporal-estimators]] —
true-scale multi-window inversion network fault index std ~unchanged (0.030%→0.033%). BUT absolute bounds here
are in the OLD 0.44×-distorted scale; multiply ~1.6–2.3× for honest dv/v (ETS bound ~0.003%→~0.005–0.007% true).
Also #1 below is SUPERSEDED on window: usable band is **2–4 s**, not 1–4 s (1–4 s diluted by pinned S, see
[[coda-window-2to4-correction]]); the S-arrival/no-late-coda physics still holds. Conclusion intact, numbers rescaled.

**Two durable findings from the ETS analysis (2026-06-06). Full log: `fault_tomography/ETS_ANALYSIS.md`.**

**1. LFE CODA WINDOW = 1–4 s ONLY (hard constraint, empirically verified).** Stacked coda envelope of clean
borehole stations (`figures/smoke_coda_envelope_check.png`): the **S arrives at t = 1.02 s** (stack t=0 is a
reference ~1 s ahead of the S, NOT the S itself); the **1–4 s window has S/N ≈ 11**; the **4–8 s window is
S/N ≈ 0.7 — BELOW the noise floor.** Weak LFE coda decays into noise by ~4 s. So **there is no usable late
coda** — classical "late coda for sensitivity" CWI is INFEASIBLE with LFEs. 1–4 s is the only usable window,
which also means dv/v precision is coda-length-limited. (Also confirms 0–1 s = pre-S/ballistic, the reason
the old 0–2 s window was wrong; see [[coda-window-standard]].)

**2. ETS produces NO station-consistent deep-interface dv/v change.** Superposed-epoch of dv/v on Cascadia
ETS onsets (ETS timing self-derived from the PNSN tremor catalog, recurrence ~12.9 mo north = real).
Battery of tests all NULL, either sign: per-segment (N/C/S) monthly, north weekly, detrend-robust, active-
family (detection-weighted), pooled-37-events sign-agnostic, and the decisive **station-by-station agreement**
test (% of 22 stations positive sits at ~50% at every lag = random). **Upper bound < ~0.002–0.003%** on the
deep-interface velocity change during ETS. **This smallness is PHYSICALLY EXPECTED, not a non-detection:** dv/v
is driven by crack opening/closing; at LFE depth (30–45 km) the ~1 GPa confining pressure clamps cracks shut, so
dv/dstress collapses with depth — the deep interface is intrinsically velocity-stress-INSENSITIVE, so a tiny
change is predicted regardless of the ETS forcing. The contrast with shallow ambient-noise ~0.1% IS the
depth-dependence of crack compliance. Reframed: *the deep megathrust is velocity-stable through ETS, as expected
where high confining pressure suppresses the crack compliance that drives shallow dv/v* — a deep-depth
constraint, not a failure. (Slow slip's low stress drop also contributes, but the depth/pressure effect dominates.)

**METHOD LESSONS (reusable):** (1) **station-by-station agreement >> region-pooling** — averaging stations
into a regional index FIRST manufactures spurious hints dominated by 1–2 stations; demand the stations agree.
(2) **Test BOTH signs** — dilatant hardening/drainage can RAISE velocity during slip; don't assume a drop.
(3) **Verify the coda window empirically** (envelope + S/N), don't assume. See [[fault-tomography-next-phase]].
