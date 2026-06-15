# LFE Feature-Fingerprint Exploration — toward a probabilistic LFE picker

*Workspace created 2026-06-13. Side exploration off the main tremorferometry tree.
Self-contained; read this first.*

## 0. GOAL (the north star)
**Find credible waveform features / fingerprints that identify LFEs**, well enough to eventually
train an **EQTransformer-style probabilistic LFE picker** (per-window/per-sample P(LFE)). Document
everything as we go. This is the supervised, learned answer to the project's central problem: the
discovery+densify detector (2-s, cc>=0.8, Z-only, envelope-peak seeding) **admits noise** (~90/family/day),
and simple thresholds can't separate emergent LFEs from impulsive/cultural/EQ contaminants.

Success = (1) a feature space where known LFEs separate from known non-LFEs; (2) a map of where each
densified family (GOLD/TRUSTED/UNDET/FAIL) lands in it; (3) a learned embedding + classifier that
outputs a calibrated LFE probability and beats the current detector on the trust-battery referee.

## 1. WHY features (the physics we exploit)
LFEs are **emergent, low-frequency-dominant (~1-8 Hz, sharp falloff >8-10 Hz), low-SNR, S-dominant,
no surface waves**. Contaminants differ in physically measurable ways:
- regular micro-EQs: impulsive onset, energy to 20-40+ Hz, clear P
- blasts: Rg/surface waves, ripple-fired spectral combs, weekday/midday timing
- cultural/anthropogenic: narrowband/harmonic or broadband, diurnal
- glitches: spikes, non-physical envelopes
The project already proved simple spectral shape separates the *spiky* class (centroid+kurtosis, AUC 0.92)
but NOT the *emergent-noise* class — so we need richer features + (eventually) learned embeddings, and
ideally polarization / multi-station info (the true discriminator).

## 2. LABELS (what makes this supervised)
- **Positives (LFE):** Lin (2023) southern-VI catalog (1.06M dets, 2005-2017) -> windows at VI stations.
  Curate clean positives: N>=6 stations, residual<~0.8 s, depth 25-45 km (interface).
  Co-located borehole **B011** (Pat Bay, EH? 100 Hz; same instrument family as the WA targets) is the
  bridge that avoids VI-broadband -> WA-borehole domain shift.
- **Hard negatives:** trust-battery **FAIL** families; **reversed-template** fakes; fingerprinted
  CULTURAL/BLAST/NATURAL families (data/family_fingerprint.csv); regional **ANSS/PNSN earthquakes**.
- **Unlabeled / to-be-placed:** every densified family's detection windows, by station & grade.

## 3. WINDOWS
- Detection time t_d (mf_<sta>_all.csv) ~= S-arrival (template anchored near S). Use a **long window**
  `[t_d - 10 s, t_d + 30 s]` (40 s) to capture pre-event noise (SNR), onset, S, and coda.
- Native fs (borehole EH? 100 Hz; broadband BH? 40 Hz). HF features (>20 Hz) borehole-only -> flagged.
- QC: drop windows with gaps/zero-fill or <90% expected length.
- Sample ~150 detections/family (high-cc + random mix) for per-family stats; keep per-detection rows for
  within-family variability ("how much do the features change").

## 4. FEATURE CATALOG (Stage-1 hand-crafted; feature_defs.py)
Spectral (Welch PSD): centroid, spread, flatness, rolloff85, peak-freq, band fractions
{1-2,2-4,4-8,8-16,16-Nyq}, HF ratio (E>8 / E[2-8]), log-log spectral slope, spectral entropy.
Envelope/temporal (Hilbert, 2-8 Hz): crest factor, env kurtosis, env skew, rise time (10->90%),
duration (>half-max), zero-crossing rate, SNR (signal/pre-noise).
Polarization (3-comp, when horizontals present): H/V energy ratio (LFEs S-dominant -> high),
rectilinearity, planarity.
(All physically motivated; interpretable; cheap; the baseline the learned embedding must beat.)

## 5. STAGED ROADMAP
- **Stage 1 — hand-crafted features** (extract_features.py): per-window features -> per-family fingerprints.
  Deliver: feat_<sta>.parquet (per-detection) + feat_<sta>_family.csv (median/IQR). Start B927.
- **Stage 2 — unsupervised structure** (cluster_explore.py): standardize -> PCA/UMAP -> HDBSCAN/KMeans.
  Map GOLD vs FAIL vs Lin-positive vs ANSS-EQ; quantify feature variability & redundancy; do FAIL
  families land apart from GOLD? Where do densified families sit?
- **Stage 3 — learned fingerprints**: spectrogram **autoencoder** (CNN, torch) -> latent vectors as
  fingerprints; cluster latents; compare to hand-crafted; self-supervised/contrastive optional.
- **Stage 4 — the picker**: supervised model (CNN/transformer on waveform or spectrogram) trained on
  Lin positives + hard negatives -> calibrated P(LFE) per window; then per-sample (EQTransformer-style)
  for picking. Validate on held-out + the trust-battery referee + reversed-template nulls.

## 6. VALIDATION / ANTI-CIRCULARITY
- Lin is **external** ground truth (independent of our detector). Use time-shifted & reversed-template nulls.
- Never train on features derived from the same cc/template that defined the family (leakage).
- Referee every cleaning/relabel with the trust battery (coda-sigma, day/night, reversed-fake max).
- Hold out stations (train VI/borehole, test WA borehole) to measure transfer honestly.

## 7. KNOWN RISKS (we are our own skeptic; Fable offline)
- **Domain shift** VI-broadband -> WA-borehole (mitigated by borehole positives at B011/B004).
- **Emergent-noise wall**: single-station spectral features may not separate emergent LFEs from emergent
  ambient noise; the real discriminator is multi-station moveout. Z-only features = a ceiling; document it.
- **Label noise**: Lin is itself ML-derived (has false positives) -> curate by N/residual/depth.
- **Component availability**: WA targets are EHZ-only on disk; 3-comp features need horizontal re-download.

## 8. STATUS LOG
- 2026-06-13: workspace created; PLAN written. Stage-1 extractor built. B927 horizontals (EH1/EH2)
  downloading (append to existing EHZ). Demo Stage-1 run on B927 (vertical) launched. torch+umap installing.
  Context from main tree: notes/2026-06-13_Notes.md (envelope-peak critique, Lin-as-ground-truth, B011 bridge,
  lag-histogram = episode co-activity not event-identity).
- 2026-06-13 (later): **RESULT 1 DONE — credible LFE fingerprint** (LIN-vs-RAND AUC 0.95-0.98 at PGC+B011,
  SNR-INDEPENDENT, physical: spectral slope/polarization/2-8Hz band/HF-depletion). See RESULTS.md.
  **RESULT 2 — B011 trust-families = coherent low-frequency HF-depleted repeaters, spectrally LFE-consistent**
  (not noise, not EQ) — reached only after catching 3 of my own artifacts (alignment confound, classifier
  extrapolation, jitter-degraded LIN stack reference). torch/umap/hdbscan installed. All Lin-region 3-comp
  downloads complete (PGC/B011/B926). KEY LESSON: align positives precisely before stacking; inspect raw
  feature distributions not just classifier P. NEXT: aligned-fair family verdict + stack polarization; ANSS
  EQ hard negatives; PGC/B926 replicate; Stage-3 autoencoder; Stage-4 picker.
- 2026-06-14: **R3-R10 complete — RESULTS.md (findings) + DECISIONS.md (rationale) + notes/2026-06-14 are the
  authoritative current record** (this PLAN's earlier sections are the original roadmap). R3 EQ hard-neg
  (LFE-vs-EQ 0.97); R4 per-window picker v0 (LFE AUC 0.975); R5 PGC transfer; R6 autoencoder (latent 0.965);
  R7 per-sample U-Net (0.85, single-station-capped); R8 continuous activity detector (ETS 0.94/quiet 0.05);
  R9 tremor-window picker (EQ/BLAST 0.97, in-tremor-noise 0.81); R10 wired into real GPU discovery (filter at
  FAMILY level; B011 contamination is in-tremor-noise not EQ/blast → picker bites harder at EQ/blast stations).
