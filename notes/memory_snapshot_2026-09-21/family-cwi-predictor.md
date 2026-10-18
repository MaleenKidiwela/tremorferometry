---
name: family-cwi-predictor
description: Template waveform SHAPE predicts whether an LFE family is CWI-suitable (GOOD) vs gappy/noise (BAD) — free pre-densify screen; spectral centroid + kurtosis rule
metadata: 
  node_type: memory
  type: project
  originSessionId: 18c98477-b5c6-4af0-9788-f240c207a136
---

Established 2026-06-10 (Opus agent + independently verified). We CAN predict, before any densify, which
LFE families give clean continuous coda (GOOD) vs gappy/bursty (BAD), from the **discovery-stage stacked
template's WAVEFORM SHAPE**.

- Cheap detection-statistics (n_det, span, active-months, monthly-CV, SNR) do NOT predict at family level —
  per-station AUC ~0.5 (pooled signal is a pure station confound). Catalog recurrence/location/depth = chance.
- **Template shape DOES:** GroupKFold-by-station (whole stations held out) AUC 0.92 combined; 0.89
  within-station deconfounded (label-shuffle control → 0.74 floor). Top features: spectral centroid, spectral
  entropy, kurtosis, crest factor — all template-intrinsic.
- Physics: GOOD = low spectral centroid (~3.4 Hz), low kurtosis (~1.2), emergent/coherent = real LFE;
  BAD = high centroid (~4.6 Hz), spiky (kurt ~7.9), noise-like contaminant cluster or micro-EQ. So it predicts
  "clean LFE template vs contaminant" — exactly the triage we want for CWI.

**OPERATIONAL RULE (free, pre-densify, verified):** skip families with **spectral centroid > 4.3 Hz AND
kurtosis > 4** → flags 11% of pool, 97% precision on BAD, removes 47% of BADs, loses 0.4% of GOODs. Aggressive
(HistGB P_good<0.3): skip 21%, remove 73% BADs, lose 4.7% GOODs. Add as a PRE-SCREEN in family selection
(MARGIN_WORKFLOW step 2f) BEFORE densify → saves compute + improves dv/v. Largely supersedes the sampled-day
densify probe for triage (reserve sampling for borderline 0.3≤P_good≤0.6).

**REFINEMENT (2026-06-10): continuity is ONE axis = template spectral quality.** Among genuine LFEs (pass the
screen), continuity (active-day fraction) is partially predictable (grpCV R²=0.31) — STILL driven by spectral
centroid (within-sta r=−0.53), NOT catalog recurrence/location/depth (r≈0, the "busy patch → steadier" hypothesis
fails). So there is NO separate repeat-rate signature: the centroid/kurtosis screen already captures most
predictable CWI-suitability; residual burst-vs-steady among clean templates is set by the slow-slip cycle (needs
direct sampled-densify measurement, not predictable).

Files: `scripts/family_predictor.py`, `notes/FAMILY_PREDICTOR.md`, `data/family_predictor_full_features.csv`
(+ perfeature_auc, importance), `figures/family_predictor_auc.png`. Labels from `data/family_quality_flags.csv`.
Caveats: NLLB has no .npz template (can't screen); feature direction flips at 2/17 stations → use the GLOBAL
rule, not per-station. sklearn 1.9.0 user-installed. See [[family-selection-rule]], [[dvv-temporal-estimators]].
