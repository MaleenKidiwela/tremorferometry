# LFE Feature-Fingerprint — RESULTS (living doc, started 2026-06-13)

Goal: credible features/fingerprints that identify LFEs → toward an EQTransformer-style
probabilistic LFE picker. Method/plan in PLAN.md. This file = findings + the skeptical trail.

---

## RESULT 1 — A credible, physical, SNR-independent LFE fingerprint  ✅ (the core deliverable)

Trained on ground truth: curated **Lin (2023)** LFEs (`LIN`, N≥6, residual<1.2, within 35 km of station)
vs **random-time** windows (`RAND`) at the *same* station. 40 s windows ([−10,+30] s about the S anchor),
hand-crafted spectral + envelope + 3-component polarization features (`feature_defs.py`).

| station | LIN-vs-RAND CV-AUC (logistic / RF) | AUC w/o amplitude features | snr Cohen's d |
|---|---|---|---|
| **CN.PGC** (broadband) | 0.95 / 0.978 | **0.978** | −0.09 |
| **PB.B011** (borehole) | 0.978 / 0.983 | **0.983** | +0.02 |

**The fingerprint (LIN vs RAND, by effect size):** steep negative **spectral slope**, high **polarization**
(planarity, rectilinearity, H/V), energy concentrated in **2–8 Hz** band fractions, **depleted HF**
(low `hf_ratio`, low `bf_16up`), low spectral centroid/rolloff. All physically the LFE signature.

**Key skeptic control — it is NOT just "signal vs silence":** removing every amplitude/SNR feature leaves
the AUC *unchanged* (0.978/0.983), and `snr` has ~zero effect size. The discrimination is **waveform
character** (spectral shape + polarization), not loudness. This is the foundation a picker needs.

Replicated at a broadband (PGC) and a borehole (B011) → not an instrument artifact.

---

## RESULT 2 — Where do the trust-battery families land? (B011) — a 3-artifact skeptical saga

Question: do B011's densified, trust-graded families (38 GOLD…) look like LFEs in this space?

**2a. Per-window scoring = noise-diluted (as expected).** Individual densified detection windows score
P(LFE)≈0.10 for ALL grades. This is the project's core truth in feature space: the lists are ~90/day and
mostly noise-matches, so random detections are mostly noise. `mean = 0.05·0.9 + 0.95·0.05 ≈ 0.09` ✓.
Per-window is too blunt — must test the COHERENT signal (stacks), like the trust battery does.

**2b. Stack test, naive → looked alarming, was an ARTIFACT.** Coherent family stacks scored P(LFE)≈0.48
(full model) — intermediate between LIN(0.95) and RAND(0.04). With the alignment-robust model (spectral +
polarization only) they *dropped* to 0.165, seemingly "not LFE". **Three artifacts, each caught:**
  1. **Alignment confound** — family windows are matched-filter-aligned (sharp stacks); my LIN windows are
     anchored at a fixed OT+11 s with real S-time jitter (smeared stacks). This unfairly inflates the
     onset/envelope features for families. → dropped envelope features (robust model).
  2. **Classifier extrapolation** — the robust score *fell*, not rose. Reason: the family stacks lie
     *beyond* the LIN training cloud on the low-frequency side, so a model trained to separate LIN(centroid
     7 Hz) from RAND(8.7 Hz) mis-scores centroid=4 Hz. The classifier score on stacks is unreliable here.
  3. **Degraded positive reference** — the LIN *stacks* are smeared by anchor jitter (centroid 7.2 Hz is
     too high for clean LFEs), so they are a poor positive for stack-comparison.

**2c. Raw spectra settle it (alignment-robust ground truth).** Mean stack spectral features:

| class | sp_slope | centroid | bf_2_4 | hf_ratio | rolloff85 |
|---|---|---|---|---|---|
| **FAM (families)** | −1.76 | **4.0 Hz** | **0.79** | **0.05** | **4.0 Hz** |
| LIN (real LFEs)    | −1.78 | 7.2 Hz | 0.38 | 0.42 | 12.3 Hz |
| RAND (noise)       | −1.02 | 8.7 Hz | 0.23 | 0.92 | 22.5 Hz |

On **every** physical low-frequency metric the family stacks are the **most** LFE-like — sharp ~3 Hz peak,
79% of energy in 2–4 Hz, near-zero HF, steepest falloff (`figures/spectra_compare_b011.png`).

**CORRECTED CONCLUSION:** B011's trust-battery families are **coherent, sharply low-frequency (2–5 Hz),
HF-depleted repeating signals — spectrally LFE-consistent.** They are NOT noise (far from RAND) and NOT
impulsive earthquakes (those would be HF-rich). The "0.165 → not LFE" was an artifact; had I stopped at the
classifier number I'd have reported the opposite of the truth. (The visual "impulsive onset" in
`stack_examples_b011.png` is the matched-filter alignment sharpening a *low-frequency* arrival, not HF.)

**Open caveat:** "coherent low-frequency repeater, S-band, HF-depleted" is the LFE signature but not a
location/depth proof. A fully airtight per-family verdict needs (i) cross-correlation-ALIGNED LIN stacks for
a fair comparison, and (ii) polarization of the family stacks (only Z was saved this run).

---

## RESULT 3 — The fingerprint REJECTS EARTHQUAKES, not just silence  ✅ (picker-critical)

Added a hard-negative class: 3,811 ANSS/USGS regional **earthquakes** (M≥1, ≤145 km, anchored at estimated
S-arrival) at B011. Pairwise RandomForest CV-AUC:

| pair | AUC |
|---|---|
| LIN vs RAND (noise) | 0.983 |
| **LIN vs EQ (earthquake)** | **0.973** (0.972 without amplitude features) |
| EQ vs RAND | 0.922 |

LFE-vs-EQ physical contrast (class means): LFEs are **lower-frequency** (centroid 7.3 vs 10.2 Hz, hf_ratio
0.53 vs 2.0, steeper slope), **more S-polarized** (hv_ratio **3.69** — highest of all classes, vs EQ 2.70,
RAND 1.43), and **emergent** (lower env kurtosis/crest) vs the impulsive EQs. `figures/multiclass_pca_b011.png`
shows LIN as a tight distinct cluster, EQ spread to one side, RAND broad in the middle.

**Implication:** the feature set supports a genuine **3-class LFE / EQ / noise** discriminator at AUC ~0.97-0.98
on every LFE-vs-X pair — the foundation for an EQTransformer-style picker (which is fundamentally per-sample
multi-class). SNR-independent throughout.

## RESULT 4 — Picker v0: a calibrated per-window LFE/EQ/NOISE model  ✅ (the deliverable's first form)

`build_picker_v0.py` trains a 3-class RandomForest (isotonic-calibrated) on B011's LFE+EQ+NOISE windows.
5-fold cross-validated:

| class | one-vs-rest CV-AUC | precision | recall | f1 |
|---|---|---|---|---|
| **LFE** | **0.975** | 0.884 | 0.901 | **0.892** |
| EQ | 0.923 | 0.825 | 0.709 | 0.763 |
| NOISE | 0.951 | 0.867 | 0.916 | 0.891 |

Overall accuracy 0.863. Top feature importances: **hv_ratio (0.156)**, env_skew, sp_slope, bf_4_8, env_crest,
bf_2_4 — polarization + spectral shape + emergent-onset, as expected. Saved `models/picker_v0_b011.joblib`
(+ `.json`, scaler + feature list + classes bundled). EQ recall is the weakest (0.71 — small EQs blur into
noise); LFE detection is the strong suit. This is the per-window precursor to the per-sample EQTransformer-style
picker (Stage 4 proper). `figures/picker_v0_plfe_b011.png` = P(LFE) histogram by true class.

## RESULT 5 — Replicates across instrument types (transfer)  ✅

PGC (CN broadband, BH? 40 Hz) vs B011 (PB borehole, EH? 100 Hz):

| pair | B011 borehole | PGC broadband |
|---|---|---|
| LIN vs RAND | 0.983 | 0.978 |
| LIN vs EQ | 0.973 | 0.970 |
| picker v0 LFE one-vs-rest AUC | 0.975 | 0.972 |

Same discriminators (polarization, spectral slope, emergent onset). So the fingerprint is **not a
single-instrument artifact**. CAVEAT: *absolute* feature values differ by instrument (PGC broadband carries
strong <2 Hz microseism the borehole attenuates → PGC centroid/band-fractions shift), so a cross-station
picker needs **per-station training or instrument-response normalization**; within-station separation is
robust either way. (PGC LFE f1 0.76 < B011 0.89 only because PGC has 940 LIN positives vs B011's 4,516.)

## RESULT 6 — Unsupervised spectrogram AUTOENCODER rediscovers the LFE fingerprint  ✅

`extract_spectrograms.py` → fixed 40×64 log-spectrograms (0–25 Hz) for LIN/EQ/RAND + family windows;
`train_autoencoder.py` → conv autoencoder (latent 32, 40 epochs, torch CPU), purely unsupervised.
**A linear probe on the learned latent separates LIN vs (EQ+RAND) at AUC 0.965** — matching the supervised
hand-crafted features. So LFE-discriminative structure emerges from the data *without labels*. UMAP of the
latent (`figures/ae_umap_b011.png`) shows LIN as a distinct cluster vs EQ/noise; individual family windows
spread broadly (noise-diluted, as in RESULT 2 — it's the stacks that are LFE-consistent). This validates that
a learned model (the basis of an EQTransformer-style picker) will find the same signal, and gives reusable
embeddings (`data/ae_latent_b011.npz`).

## RESULT 7 — EQTransformer-style per-sample picker (Stage 4 v1): works, but single-station is limited

Built the full per-sample pipeline: `build_picker_dataset.py` (3,509 labeled 41 s / 3-comp / 50 Hz segments
with per-sample [LFE, EQ, NOISE] Gaussian targets at Lin-LFE and ANSS-EQ arrivals) → `train_phasenet.py`
(1D U-Net, PhaseNet/EQTransformer-style, softmax per sample, early stopping).

**Held-out per-sample AUC: LFE 0.848, EQ 0.853.** The P(LFE) trace genuinely tracks LFE presence
(`figures/phasenet_examples_b011.png` — the predicted stream rises around true LFEs). BUT **pick-level
precision is poor** (precision 0.03, recall 0.47 at peak>0.3): broad, noisy bumps with many false peaks.

**Why — two honest, structural reasons (not just tuning):**
1. **Imprecise labels.** Lin gives origin times; my per-sample target sits at OT+11 s (S estimate) ± several
   s of real jitter → I had to widen the label to σ=1.5 s → the model can only produce broad bumps → poor
   localization. A production picker needs **precise arrival times** (matched-filter the Lin catalog at the
   station, or use family matched-filter detection times).
2. **Single-station + emergent signal = the wall.** LFEs are emergent and low-SNR; the field detects them by
   **network matched-filter / summation**, not single-station per-sample picking, precisely because one
   station can't localize an emergent arrival. This is the same "emergent-noise wall" flagged at the project
   level — a single-station per-sample picker has a real ceiling.

**Verdict:** the architecture works and produces an LFE-responsive probability stream (AUC 0.85), but the
**per-WINDOW classifier (RESULT 4, AUC 0.975) is the strong, deployable detector now**; the per-sample picker
needs (a) precise matched-filter labels and (b) ideally **multi-station input** to reach production quality.
Model saved `models/phasenet_b011.pt`.

**Follow-up test (precise labels) — confirms the limit is STRUCTURAL, not tuning.** I built cross-correlation
pick refinement (`refine_lin_picks.py`): a global LFE reference aligned only at **mean CC 0.22** (lag std
2.2 s) — i.e. a *single* LFE template doesn't fit because **LFEs are patch-specific** (different families =
different waveforms) and individually emergent/noisy. Retraining the per-sample picker on these refined labels
left AUC ~0.83 (unchanged). So the ceiling is the **single-station + emergent-signal wall**, not label jitter.
Truly precise single-station picks would require per-family templates (= the existing family detector,
noise-diluted); the real lever is multi-station network input.

## RESULT 8 — Deployed continuous scan: the per-window picker is a strong LFE-ACTIVITY detector  ✅

`scan_continuous.py` slides the v0 model across whole days. Auto-picked busy vs quiet:

| day | Lin LFEs | mean P(LFE) |
|---|---|---|
| 2007-029 (ETS) | 245 | **0.941** |
| 2012-172 (quiet) | 0 | **0.047** |

**20× day-level contrast** (`figures/scan_b011.png`): P(LFE) is pinned high through the ETS day and flat near
zero on the quiet day; EQ bumps appear separately. So as a deployed detector it is an excellent **LFE-activity
/ tremor** discriminator. Caveat (honest): on the ETS day it's high *most* of the day, and the *hourly*
correlation with the sparse discrete Lin catalog is weak (r=0.15) — because tremor is near-continuous LFE
activity, so the picker tracks the *activity envelope*, not individual cataloged picks. (First-glance "it fires
everywhere = broken" was wrong; the quiet-day 0.047 proves it discriminates — it's the base-rate that's high
on an ETS day because the day really is LFE-rich.) Discrete-event timing is the per-sample/network job (R7).

## RESULT 9 — THE GOAL: pick LFEs inside PNSN tremor windows, reject EQ/blast/cultural  ✅ (with one hard residual)

`build_tremor_picker.py` — classes LFE (Lin) / EQ (ANSS earthquake) / BLAST (ANSS explosion, 537) / TNOISE
(random in-tremor-window times not near an LFE = the ambient/cultural the discovery detector grabs).

| LFE vs … | AUC |
|---|---|
| **EQ (earthquake)** | **0.968** |
| **BLAST (explosion)** | **0.968** |
| **TNOISE (in-tremor ambient/cultural)** | **0.809** ← the hard one |
| pooled contaminants | 0.865 |

**Key finding:** EQs and blasts — the impulsive, HF-rich confusers — are **cleanly rejected (0.97)**. The hard
residual is **in-tremor ambient/cultural noise (0.81)**, because it shares the LFE low-frequency band (it's the
emergent-noise wall again). Base-rate-aware operating point (recall 0.8): the picker **enriches LFEs ~3×**
(pool 5%→15%, 10%→28% purity) and removes EQ/blast, but does not alone produce a pure LFE set — the TNOISE
residual remains. Model saved `models/tremor_picker_b011.joblib`. Figure `figures/tremor_picker_b011.png`.

**Why this still solves the discovery problem:** discovery doesn't need per-candidate purity — it **clusters**
candidates (cc≥0.8, ≥3 members/≥3 yr). Incoherent TNOISE does NOT cluster into persistent families, so
clustering removes it for free. What clustering does NOT remove is **repeating EQ multiplets and quarry-blast
sequences** — which DO form spurious "families" (the original contamination!). This classifier kills exactly
those (0.97). So: **tremor-picker (rejects EQ/blast) + clustering (rejects incoherent TNOISE) = clean families.**

## RESULT 10 — Wired into real GPU discovery (filtered vs unfiltered) on B011

Scored 415k B011 candidates (2010-13) with the tremor-picker, ran `discover_gpu` on the full pool vs
the filtered pool. Classified every resulting family by member coincidence (±6 s) with Lin LFEs / ANSS
earthquakes / blasts.

**Finding A — filter at the FAMILY level, not the candidate level.** Per-candidate filtering (P≥0.4)
fragments families (members lost → fall below ≥3-member/≥3-yr) → kept only **19%** of Lin-confirmed LFE
families for a 26%→38% purity gain (bad trade). Clustering first, then scoring each family by mean member
P(LFE), is strictly better (same purity at higher retention):
| family-P(LFE) thr | families | LFE retained | purity | contaminant removed |
|---|---|---|---|---|
| 0.2 | 2874 | **72%** | 32% | 46% |
| 0.3 | 1745 | 49% | 35% | 69% |
| 0.4 | 940 | 29% | 38% | 84% |

So a *gentle* family-level filter (thr~0.2-0.3) removes ~half the contaminant families while keeping
most real LFE families — a usable post-discovery quality filter.

**Finding B — at B011 the contamination is NOT EQ/blast.** Of 4,875 baseline families, only **5 were EQ,
0 BLAST**; the bulk (3,616) is **UNK** (coherent repeaters not in Lin = in-tremor cultural/ambient AND
genuine LFE families Lin doesn't catalog). The picker's *strength* is rejecting EQ/blast (0.97) — but
B011 barely has those; its contaminant is the *hard* in-tremor/uncatalogued class (0.81). So the picker
would help **much more at an EQ/blast-contaminated station** (project's anthropogenic B018/B033/CPW/GNW or
natural-seismicity B935/B017/B022). Also: Lin-confirmation **underestimates** true purity (many UNK are
real LFE families outside Lin's catalog), so the real purity gain is better than the table's % implies.

**Verdict on "does it make better seeds?":** Yes — best as a **family-level** filter; it cleanly removes
EQ/blast multiplet seeds (the classic contamination) and gently trims low-P(LFE) families. Its leverage is
**station-dependent**: large where EQ/blast contamination dominates, modest at LFE-pure sites like B011.

## RESULT 11 — Longer templates do NOT clean up clustering (a tested correction)

Hypothesis (mine): a longer / 3-comp template raises N_eff → cleaner clustering, fewer spurious families.
TESTED the length half on B011 baseline candidates (same pool, vary template length + cc):

| config | families | LFE-confirmed | purity |
|---|---|---|---|
| 2 s, cc 0.80 (baseline) | 4,875 | 1,254 | **25.7%** |
| 4 s, cc 0.80 | 37 | 3 | 8.1% |
| 4 s, cc 0.65 | 2,973 | 721 | 24.3% |
| 6 s, cc 0.65 | 254 | 45 | 17.7% |

**Result: no length beats 2 s for clustering purity.** Why my N_eff argument failed *for clustering*:
- The N_eff/longer-template benefit applies to **matched-filter DETECTION** (matching a known template to
  continuous data — densify), where added *coherent* samples raise significance.
- For **CLUSTERING candidates into families**, the added samples are LFE **coda, which decorrelates between
  occurrences**. So cc≥0.8 over a long window is unmeetable for genuine emergent LFEs (4 s → families collapse
  to 37), and the few that survive are **coherent IMPULSIVE contaminants** (machinery / repeating EQs) → 8%
  purity (the "high-cc selects against LFEs" trap again). Lowering cc to compensate just restores ~the 2 s
  purity with fewer families.
- **Implication: ~2 s is approximately matched to the *coherent duration* of an LFE waveform** (the coda
  decorrelates beyond it). Longer adds incoherent samples → helps neither clustering nor (genuine-signal)
  detection. So 2 s is a defensible choice for the *clustering* step; the real lever for cleaner seeds is the
  **picker** (reject EQ/blast) + per-detection LFE character, NOT template length.
- **Untested remaining lever: 3-component** templates add *independent coherent* info (LFE S is strong on
  horizontals, and the S-burst is coherent unlike the coda), so it could help where length doesn't — but it
  needs a discover_gpu code change (Z-only → 3-comp CC). Not yet run.

## METHODOLOGICAL LESSONS (carry forward to the picker)
- **Positives must be precisely aligned.** Anchor-jitter smearing degrades stacked LFE references and
  inflates their HF. For Stage-4, align positives by cross-correlation / matched filter, or use Lin arrival
  info — do not rely on a fixed OT+offset for stacking.
- **Classifiers extrapolate badly.** Always inspect raw feature distributions, not just P; the family stacks
  were *more* extreme than the positive cloud.
- **Per-window ≠ family.** Noise dilution (~90/day) means individual densified windows can't judge a family;
  use coherent stacks.
- Single-window LFE-vs-RAND is solid (jitter just shifts the LFE within a 40 s window; doesn't cancel it).
  Only *stacking* is alignment-sensitive.

## FILES (complete inventory — updated through R10)
**Docs:** PLAN.md (goal/roadmap), RESULTS.md (this, findings R1-R10), DECISIONS.md (rationale + rejected paths).

**Code (lfe_features/*.py):**
- `feature_defs.py` — the 23 hand-crafted features (spectral / envelope / 3-comp polarization).
- `extract_features.py` — features at densified-family detection windows (R1).
- `extract_lin_positives.py` — LIN (Lin LFE) + RAND ground-truth windows (R1).
- `extract_eq_negatives.py` — ANSS events → EQ or BLAST windows (--eventtype/--label) (R3/R9).
- `extract_tremor_noise.py` — in-tremor ambient (TNOISE) windows (R9).
- `analyze_lfe_features.py` — LIN-vs-RAND AUC + SNR control + GOLD/FAIL placement (R1/R2).
- `analyze_multiclass.py` — LFE/EQ/noise multiclass + PCA (R3/R5).
- `cluster_explore.py` — PCA/separation/feature-redundancy explorer (R2 framework).
- `family_stack_test.py` — coherent family stacks, full+robust scoring, saves stacks npz (R2).
- `plot_stack_examples.py` — stack waveform+spectrum panels (R2).
- `build_picker_v0.py` — calibrated 3-class per-window LFE/EQ/NOISE picker (R4).
- `extract_spectrograms.py` + `train_autoencoder.py` — spectrogram conv-AE, learned latents (R6).
- `build_picker_dataset.py` + `train_phasenet.py` — per-sample U-Net picker dataset + training (R7).
- `refine_lin_picks.py` — cross-correlation pick refinement (R7 follow-up; CC 0.22 → rejected).
- `scan_continuous.py` — deploy v0 as a continuous LFE-activity scanner (R8).
- `build_tremor_picker.py` — LFE vs EQ/BLAST/TNOISE picker for tremor windows (R9).
- `score_candidates.py` — tag discovery candidates with P(LFE) → baseline/filtered parquets (R10).
- `discovery_filter_demo.py` — CPU envelope-peak+cluster proof (hit scale wall — needs GPU) (R10).
- `compare_discovery.py` — classify GPU-discovered families by ground truth, baseline vs filtered (R10).
- `fam_level_filter.py` — family-level P(LFE) filter sweep (R10, the right filter placement).

**Models (lfe_features/models/):** `picker_v0_{b011,pgc}.joblib(.json)` (R4 per-window),
`phasenet_b011.pt` (R7 U-Net), `tremor_picker_b011.joblib(.json)` (R9 tremor picker).

**Key data (lfe_features/data/):** `feat_ref_{b011,pgc}.parquet` (LIN/RAND), `feat_{b011,b927}.parquet`
(densified), `feat_{eq,blast,tnoise}_b011.parquet` (negatives), `spec_b011.npz`+`ae_latent_b011.npz` (R6),
`picker_ds_b011.npz` (R7), `stacks_b011.npz`+`stack_plfe_b011.csv` (R2), `lin_precise_b011.csv` (R7),
`fam_level_b011.csv`/`disc_filtered_classified_b011.csv` (R10). Discovery I/O in main `data/`:
`b011_cand_{baseline,filtered}.parquet`, `b011_disc_{baseline,filtered}.{npz,summary.csv,members.parquet}`,
`{eq,blast}_catalog_b011.csv`.

**Figures (lfe_features/figures/):** `spectra_compare_b011.png` (key, R2), `multiclass_pca_{b011,pgc}.png`
(R3/R5), `picker_v0_plfe_{b011,pgc}.png` (R4), `ae_umap_b011.png` (R6), `phasenet_examples_b011.png` (R7),
`scan_b011.png` (R8), `tremor_picker_b011.png` (R9), plus `stack_examples_b011.png`/`stack_plfe_b011.png`/
`plfe_by_grade_b011.png`/`featspace_pca_b927.png` (R2).

**Cross-refs:** main-tree notes/2026-06-13_Notes.md + notes/2026-06-14_Notes.md; memory lfe-feature-fingerprint.

## NEXT (refined after R1-R7)
DONE: R3 hard negatives (EQ), R4 per-window picker, R5 PGC transfer, R6 autoencoder, R7 per-sample U-Net v1.
Highest-leverage remaining, in order:
1. **Precise LFE labels** — matched-filter the curated Lin catalog at the station (or use family matched-filter
   detection times) to get arrival times to ±0.1 s. This directly lifts the R7 per-sample picker (its ceiling
   is label jitter) AND fixes the R2 alignment confound (aligned LIN stacks → airtight family verdict).
2. **Multi-station input** — the field's SNR comes from the network; a 2-3 station per-sample picker should
   beat single-station (addresses the emergent-signal wall). Re-add network coincidence (project Phase C).
3. Family-stack polarization (only Z saved in R2) for a complete family verdict.
4. Replicate the EQ/picker stack at B926 (Lin region) and transfer-test to WA borehole B927 (horizontals on disk).
5. Productionize: deploy the R4 per-window picker as a scanner now (it's the strong detector); compare its
   detections to the trust-battery families and to a held-out Lin window.
