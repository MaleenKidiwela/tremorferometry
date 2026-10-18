---
name: coda-window-2to4-correction
description: "The 1-4s coda window DILUTES dv/v ~50x toward zero (pinned direct-S dominates); the correct band is 2-4s, measured via causal rolling stack + SVD-Wiener"
metadata: 
  node_type: memory
  type: project
  originSessionId: 400c26f9-17ff-46fc-914c-9c6a7671906d
---

**MAJOR METHOD CORRECTION (2026-06-07, user-driven). The 1–4 s coda window we used for the entire
margin SUPPRESSES dv/v ~50× toward zero. The correct band is ~2–4 s.** This qualifies/overturns the
small-dv/v story (incl. the ETS null) — those were measured in a diluted window, not a stable fault.

**The dilution, quantified (B927 family c660):**
- dv/v(1–4 s) std 0.016% vs dv/v(2–4 s) std 0.173% → **10.8× by scatter.** Regression
  dvv(1–4)=**0.018**·dvv(2–4) → 1–4 s retains only **2% of the real coda signal** → **54× suppression.**
  Peak-to-peak 0.085% (1–4) vs 1.18% (2–4).
- **Mechanism (the "why", airtight):** in the 1–4 s window **98% of the energy is the direct S (1–2 s
  slice); the actual coda (2–4 s) is 2%.** The stretch is a correlation (energy-weighted dot product),
  and the S is **PINNED at zero shift by the matched-filter template alignment** (it's the source's locked
  direct arrival). So the cc-vs-ε peak is dragged to ε≈0 by the dominant motionless S. Prediction
  0.98·0 + 0.02·(coda) = 2% retention = the measured 0.018 slope **exactly.** Flatness = measuring the
  pinned S, NOT a stable fault.

**The correct band = 2–4 s.** Two constraints fix it: (lower) get OFF the pinned direct-S (S at ~1.0 s,
verified envelope peak 0.97 s) → start ≥~2 s; (upper) stay INSIDE the coherent coda → end ≤~4 s. The
coherence floor is real: even with max stacking + SVD-Wiener, 3–5 s cc 0.63 / 17% usable, 4–6 s cc 0.51 /
3% — the coda is physically DECORRELATED past ~4 s (a coherence statement, not amplitude; supersedes the
old envelope-amplitude argument in [[ets-null-and-coda-window]]). 2–5 s also works (traces same excursions)
but cc 0.76 < 2–4's 0.81 — 2–4 is cleanest. Bounds are family/station-dependent; re-verify per station.

**THE METHOD (what makes off-S coda measurable; single daily off-S coda only coheres ~2% of days):**
per family — (1) **30-day CAUSAL rolling WAVEFORM stack** (trailing: each point = mean of the 30 daily
stacks up to & including that date, none after; count-based 30 stacks ≈ 30 days; cumsum); (2) **SVD-Wiener
filter** the rolling-stack matrix (row-normalise, random-matrix noise floor σ_n·(√ne+√ns), Wiener gain
sv²/(sv²+floor²), keep components above floor); (3) **stretch 2–4 s** vs the filtered all-time reference.
The rolling stack is the dominant lever (2%→50% coherent days); SVD-Wiener adds polish (50%→60%) without
changing curve shape (denoise, not fabricate). **NEVER median dv/v values — the stacking IS the smoothing.**

**⚠ ROLLING MUST BE CALENDAR-BASED, NOT COUNT-BASED (2026-06-08 fix).** The original
`dvv_2to4_roll30svd.py` used a 30-*stack* trailing window (cumsum), which on gappy stations stitches
across data gaps: PGC's "30-day" windows had a *median 77-day* span, 81% >45 days, reaching back up to
285 days to collect 30 stacks (clean stations like B927/GNW were fine, median 29 d). Fix =
`scripts/dvv_roll30cal.py` (+ `run_cal_allstations.sh`): TRUE 30-CALENDAR-DAY trailing window
(searchsorted, min 5 stacks, else the date gets no point). Clean stations unchanged (100% kept, median 30
stacks); gappy ones honest (PGC 92% kept, median 18 stacks, cc 0.70→0.64 — the old high cc was the stitching
artifact). Output `data/daily_dvv_{S}_{1to3,2to4,3to5}_cal.csv` (windows 1-3, 2-4, 3-5). The `_cal` set
supersedes the `_roll` set for the map/analysis.

**Rollout:** `scripts/dvv_2to4_roll30svd.py` (parallel across families) + `scripts/run_2to4_roll_allstations.sh`
→ `data/daily_dvv_{STATION}_2to4_roll.csv` for the 28 finalized stations (GNW uses GNW_e23 npz, PGC uses
PGC_despike). See also [[coda-window-standard]] (now superseded: 1–4 s was diluting every station).

**FULL TOMO REDONE IN 2–4 s (2026-06-07) — and the ETS null SURVIVES the correction.** Deseasoned 2–4 s
(per-family annual+semiannual removal; `data/daily_dvv_{S}_2to4_des.csv`) → tensor `a2_2to4des.py`
(prefix `tomo24des_*`) → 4-D inversion `fault_4d_tomo24des.npz` (kernel ALSO moved to 2–4 s via
`CODA_W1/W2` env in forward.py — checkerboard 0.77, unchanged because the deep self-imaging lobe is anchored
at the source cell regardless of window). **Kernel-selection CV (`kernel_cv.py`): all kernels indistinguishable
(~36–39% held-out VR; site terms soak the signal so the fault kernel is weakly constrained) — single-scatter
2–4 s is justified physically, not by a data win.** **Robust 4-D result (≥4-station cells only): network
interface δβ/β std 0.019%, FLAT through ETS (corr w/ tremor 0.03–0.08) → ETS null CONFIRMED in the un-diluted
window.** So the 50× dilution was NOT hiding an ETS signal; the deep fault really is velocity-stable.
Honest bound loosened 0.003% (flattered) → ~0.02% (real network std), still ≪ shallow ~0.1%. The raw movie's
±0.5% swings were UNDER-CONSTRAINED southern cells (B036/B935/B039) echoing one noisy station — collapse to
±0.08% with ≥4 stations. Robust cells are MOSTLY NORTHERN → this is a northern-interface statement; southern
Cascadia is UNRESOLVED (needs more stations, not better processing). Also: full-GNW (1995+) 2–4 s shows NO
Nisqually-2001 precursor (−0.034%, within 0.27% noise; sparse early coverage) → strengthens
[[nisqually-drop-artifact-check]] (founding drop = response-removal artifact). Movie `figures/fault_movie_2to4.mp4`.
