---
name: lfe-feature-fingerprint
description: "LFE feature-fingerprinting toward an EQTransformer-style LFE picker — workspace lfe_features/, credible fingerprint found, trust-families are LFE-consistent"
metadata: 
  node_type: memory
  type: project
  originSessionId: c508ae52-cfc6-40fd-99ec-a76553165d3a
---

New project thread (2026-06-13, user-set goal): find credible waveform features/fingerprints that identify
LFEs → build an **EQTransformer-style probabilistic LFE picker**. Self-contained workspace **`lfe_features/`**
(PLAN.md, RESULTS.md, 7 scripts). Uses **Lin (2023) catalog as external ground-truth POSITIVES** and the
Lin-region stations **PGC (broadband) / B011, B926 (boreholes)** — B011 is the key bridge (graded, borehole,
co-located with PGC in the densest part of Lin's catalog; same instrument family as the WA targets).

**RESULT 1 (solid): a credible, physical, SNR-INDEPENDENT LFE fingerprint.** Curated Lin LFEs (LIN: N≥6,
residual<1.2, ≤35 km) vs random-time windows (RAND), 40 s windows, hand-crafted spectral+envelope+3-comp
features → LIN-vs-RAND CV-AUC **0.95–0.98 at both PGC and B011**. AUC is UNCHANGED when all amplitude/SNR
features are dropped (snr Cohen's d ≈ 0) → it's waveform CHARACTER not loudness. Fingerprint = steep negative
spectral slope, high polarization (H/V, planarity, rectilinearity), energy in 2–8 Hz, HF-depleted.

**RESULT 2 (verified): B011 trust-battery families are coherent, sharply low-frequency (2–5 Hz), HF-depleted
repeaters — spectrally LFE-CONSISTENT** (not noise, not impulsive-EQ). Reached only after catching 3 of my
own artifacts; the naive classifier "0.165 → not LFE" was WRONG.

**RESULTS 3–6 (same day):** R3 added ANSS regional EARTHQUAKE hard-negatives → **LIN-vs-EQ AUC 0.973**
(fingerprint rejects EQs; LFEs are lower-freq + more S-polarized hv_ratio 3.69>2.70 + emergent). R4 = **picker
v0**, a calibrated 3-class LFE/EQ/NOISE RandomForest, **LFE one-vs-rest AUC 0.975, f1 0.89**, top feature
hv_ratio — saved `lfe_features/models/picker_v0_{b011,pgc}.joblib`. R5 = **transfers** to PGC broadband
(LIN-RAND 0.978, LIN-EQ 0.970); but absolute feature values differ by instrument (broadband microseism) →
per-station train or response-normalize for cross-station. R6 = **unsupervised spectrogram autoencoder**
(torch) latent linear-probe LIN-vs-(EQ+RAND) **AUC 0.965** → LFE structure emerges WITHOUT labels.

**Why:** Don't repeat these. (1) **Align positives precisely BEFORE stacking** — Lin windows anchored at a
fixed OT+offset have S-time jitter that smears stacks and inflates their HF → degrades the LFE reference
(LIN stack centroid 7 Hz, too high). (2) **Never trust a classifier P at the edge of the training cloud** —
family stacks were *more* low-freq than the LIN cloud → extrapolation scored them low; the RAW spectra
(figures/spectra_compare_b011.png) told the truth. (3) **Per-window can't judge a family** (~90/day noise
dilution) → use coherent stacks (the trust-battery logic).

**RESULT 7 (Stage-4 capstone): EQTransformer-style per-sample U-Net picker** (build_picker_dataset.py +
train_phasenet.py, 1D U-Net, per-sample [LFE,EQ,NOISE] softmax). Held-out per-sample **LFE AUC 0.848**, EQ
0.853; the P(LFE) trace tracks LFEs BUT pick-level precision is poor (broad/noisy). Two STRUCTURAL limits:
(1) imprecise labels — Lin gives origin times, my per-sample target is OT+11s ± several s jitter → wide labels
→ poor localization → need PRECISE matched-filter picks; (2) single-station + emergent LFE = the wall (the
field uses NETWORK matched-filter). So the per-WINDOW picker v0 (AUC 0.975) is the deployable detector NOW;
the per-sample picker needs precise labels + multi-station input. Model models/phasenet_b011.pt.

**RESULT 8: deployed v0 per-window picker as a continuous scanner** (scan_continuous.py) → ETS day mean
P(LFE) **0.941** vs quiet day **0.047** (20× contrast) = strong LFE-ACTIVITY/tremor detector
(figures/scan_b011.png); tracks the tremor envelope not sparse discrete catalog picks (hourly r=0.15).
**Precise-label test (R7 follow-up) FAILED**: a global LFE template aligns at only CC 0.22 → LFEs are
patch-specific + emergent → no single-station precise labels → per-sample picker stuck ~0.83. CONCLUSION: the
real lever is MULTI-STATION (network), as the field does; single-station per-sample LFE picking is structurally
capped. Deliverables complete: credible fingerprint + per-window picker (deployed) DONE.

**RESULT 9-10 (the actual goal — pick LFEs in PNSN tremor windows, reject EQ/blast/cultural):** Added BLAST
(ANSS explosions) + TNOISE (in-tremor ambient) classes. tremor_picker (build_tremor_picker.py): **LFE-vs-EQ
0.97, LFE-vs-BLAST 0.97, LFE-vs-TNOISE 0.81** (the hard residual = in-tremor low-freq ambient, shares LFE
band). Wired into real GPU discovery (score_candidates.py → discover_gpu on filtered vs baseline, B011
2010-13). **Filter at the FAMILY level, NOT candidate level** (candidate-level fragments families: kept only
19% of LFE families). Family-level (cluster → score family by mean member P(LFE)): thr 0.2 keeps 72% LFE
families / removes 46% contaminants; thr 0.3 keeps 49% / removes 69%. **B011's contamination is NOT EQ/blast**
(only 5 EQ/0 blast families of 4875) — it's the hard in-tremor/uncatalogued class → picker helps only modestly
HERE; it would help MUCH more at EQ/blast-contaminated stations (anthropogenic B018/B033/CPW/GNW; natural-
seismicity B935/B017/B022). Lin-confirmation UNDERESTIMATES purity (many UNK families = real LFEs outside Lin).
Picker model: models/tremor_picker_b011.joblib. ALSO: continuous per-window picker (R8) = strong activity
detector (ETS 0.94 vs quiet 0.05); per-sample U-Net (R7) single-station-capped.

## STATE 2026-06-18 (B011 full pipeline + categorization) — see notes/2026-06-18_Notes.md
**THE 23 FEATURES** (lfe_features/feature_defs.py, FEATURE_ORDER): spectral(13)= sp_slope, sp_centroid, sp_spread,
sp_peakfreq, sp_rolloff85, sp_flatness, sp_entropy, bf_1_2/bf_2_4/bf_4_8/bf_8_16/bf_16up (band-energy fractions),
hf_ratio; envelope/onset(7)= env_crest, env_kurt, env_skew, rise_time, duration, zcr, snr; polarization(3,3-comp)=
hv_ratio, rectilinearity, planarity. Heavy lifters: hv_ratio(top), sp_slope, 2-8Hz band fractions, hf_ratio,
sp_centroid, env_skew/crest/kurt. SNR-INDEPENDENT (dropping snr/env_crest leaves AUC unchanged). LFE = low-freq
(centroid ~3Hz), steep neg slope, HF-depleted, S-dominant (high H/V), emergent.

**THE PICKER *IS* THE CATEGORIZER** (model lfe_features/models/tremor_picker_b011.joblib, 4-class LFE/EQ/BLAST/
TNOISE; polarization is INSIDE it). Lin-catalog coincidence = EXTERNAL validation only, coverage-limited (no
Olympic side) — do NOT categorize by it. The spiky-shape screen (centroid>4.3 & kurt>4) only catches EQ/blast,
NOT the low-freq in-tremor-noise contaminant. To categorize a FAMILY: run the picker on its coherent 3-comp STACK
(seed-member stacks of 3-9 work for a first pass; DENSIFIED stacks of 1000s are cleaner).

**B011 DISCOVERY THRESHOLDS** (candidates = PNSN-tremor-window envelope-peaks 2010-13, scored by picker, filtered,
clustered cc>=0.8/>=3memb/>=3yr): P>=0.4 -> 610 families; **P>=0.5 -> 314**; P>=0.7 -> 52. Seed clusters are tiny
(3-9 members) -> need densify for clean stacks/dv/v.

**HEADLINE RESULT: of the 314 P>=0.5 templates, picker says 274 REAL LFE / 35 earthquake-multiplets / 5 in-tremor-
noise; 126 high-confidence (family-stack P(LFE)>=0.7), 19 at >=0.9.** (data/family_picker_p50.csv). First-pass from
3-9-member seed stacks; densify (mf_b011p50_, running) gives the sharper re-categorization.

**DENSIFY + dv/v DONE for the 52 (P>=0.7):** densify scripts/densify_gnw_gpu.py (--fs 100, --despike-mad 8,
out-prefix mf_b011p70_) -> mf_b011p70_all.csv (17M det) -> build_long_window_resp.py (--no-deconv --fs 40
--min-det 20) -> long_window_daily_B011p70.npz (183,876 daily stacks) -> **dvv_2to4_roll30svd.py** (30-day rolling
+ SVD, 2-4s; the CORRECT method, "stacking IS the smoothing, NO median-of-dvv") -> daily_dvv_B011p70_2to4_roll.csv
(mean cc 0.97). Result: clean ANNUAL SEASONAL ~±0.2% (shallow thermo-hydrological), expected. DESEASONING: harmonic
(annual+semiannual) fit LEFT ~1/3 of seasonal (only 94% var); use **day-of-year CLIMATOLOGY subtraction** (100% var
removed) — residual = coherent multi-year drift (2010-11 bump, std 0.05%), likely interannual shallow climate.
PLOT RULE: honest non-clipped y-axis; dvv column is a FRACTION (×100 for %); avoid 60-d median (user banned).

**3-COMPONENT CHECK:** one Lin-confirmed family (c86) stacked VERTICAL-dominant (H/V 0.29) not S-dominant — flagged
the "continuous + vertical-dominant" suspect class (daily-repeater grid: most families EPISODIC at ETS, a handful
c86/c177/c9/c88/c24 fire continuously = suspect). cc>=0.95 detections are 834x ETS-episodic (genuine); the ~90/day
all-cc floor is the detector noise (continuity in raw lists is NOT a red flag — expected).

**MAPPING:** no cartopy on pod; fetch coastline geojson (raw.githubusercontent nvkelso/natural-earth-vector
ne_10m_coastline.geojson). cas_slab2_input is scattered EQ data NOT a slab surface (don't contour it). Families
sit in the VI ETS tremor band, median ~80km W/NW of B011 (locations are PNSN tremor-bin tags, not relocated).

## How to apply / NEXT
Resume from `lfe_features/RESULTS.md` §NEXT + notes/2026-06-18. Open: let densify(314, mf_b011p50_) finish ->
re-categorize 314 from clean stacks + dv/v on the real subset. Bigger levers: run picker filter at an
EQ/blast-contaminated station (B018/B935) where 0.97 EQ/blast rejection bites; MULTI-STATION per-sample picker
(=project Phase C). Relates to [[detection-cleaning-circular-t2a-gap]] (supervised Phase-C path),
[[no-claude-coauthor-trailer]]. GIT: lfe_features/ work is on main (authored Maleen Kidiwela <seismic@uw.edu>,
NO Claude trailer); large data/models gitignored.
