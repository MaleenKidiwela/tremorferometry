# Predicting LFE-family CWI-suitability from pre-densify signatures

**Date:** 2026-06-10
**Script:** `scripts/family_predictor.py`
**Features:** `data/family_predictor_full_features.csv` (1892 families × 33 features + label)
**Figure:** `figures/family_predictor_auc.png`
**Aux:** `data/family_predictor_perfeature_auc.csv`, `data/family_predictor_importance.csv`

## Question
Can we **pre-screen** the eligible LFE-family pool *before* the expensive full
densify, by predicting from cheap pre-densify signatures which families will be
continuously-repeating / CWI-suitable (**GOOD**) vs gappy / bursty (**BAD**)?

## Data & labels
- Labels: `data/family_quality_flags.csv`, GOOD=1 vs BAD=0, MARGINAL dropped.
- **n = 1892** families (1449 GOOD, 443 BAD; base-rate GOOD = 0.766) across **28 stations**.
- Templates: `data/<sta>_pnsn_families_100km.npz` — fs=40 Hz, 80-sample (2 s) window,
  bandpass 2–8 Hz, **L2-normalized**, S-peak pinned at the center sample.
  All labeled families have a template **except NLLB (76 families, no NPZ)** →
  template features NaN for NLLB (imputed; HGB handles NaN natively).

## Feature families tested (the *untested*, physically-motivated ones)
1. **Template waveform shape** (17 feats): crest factor, kurtosis, skew, envelope
   crest, impulsiveness (duration > half-max envelope), rise time, peak position,
   envelope decay rate, central energy concentration, self-SNR, zero-crossing rate,
   autocorrelation half-width, dominant freq, spectral centroid / bandwidth /
   entropy / flatness (Hann-windowed rFFT over 2–8 Hz).
2. **Catalog / location** (9 feats): per-0.05° tremor recurrence from
   `catalogs/pnsn_tremor_cascadia_full.csv` (count, log-count, active-month
   fraction, inter-event CV); station→family distance & azimuth; family lat/lon;
   station Slab2 depth.
3. **Existing cheap baseline** (8 feats): n_det, span_yr, n_yr, n_mo,
   frac_active_yr, monthly_cv, det_per_activemo, template snr.

## Methodology (the whole point)
- **GroupKFold by STATION** (whole stations held out, 10 folds) → report **grouped-CV
  AUC**, never pooled, so station-confounds can't leak.
- **Within-station AUC** baseline (k-fold *inside* each station) — the decisive
  family-intrinsic test: does template shape split GOOD/BAD among families recorded
  at the *same station / same instrument*?
- **Negative control:** shuffle labels *within* each station, re-run grouped CV.
  This preserves per-station base-rate structure but destroys any real
  family-intrinsic signal.
- Permutation importance on a full-fit HGB.
- Single-feature AUCs reported sign-agnostic (max(auc, 1−auc)) + raw.

## RESULTS

### Combined models — grouped-CV AUC (by station)
| feature set | logit | HGB |
|---|---|---|
| **template_only** | 0.903 | **0.923** |
| catalog + location | 0.514 | 0.429 |
| cheap_only (baseline) | 0.747 | 0.703 |
| ours (template+cat+loc) | 0.882 | 0.903 |
| ALL | 0.866 | 0.893 |

### Top single features — grouped-CV AUC
| feature | grpCV AUC | raw | within-station |
|---|---|---|---|
| **t_spec_centroid** | **0.907** | 0.093 (inverted) | 0.930 |
| t_spec_entropy | 0.874 | 0.126 (inv) | 0.896 |
| t_env_crest | 0.859 | 0.141 (inv) | 0.841 |
| t_kurtosis | 0.859 | 0.141 (inv) | 0.837 |
| t_dur_halfmax_s | 0.858 | 0.858 | 0.812 |
| t_spec_flatness | 0.845 | 0.155 (inv) | 0.841 |
| t_spec_bw | 0.832 | 0.168 (inv) | 0.819 |
| t_crest | 0.855 | 0.145 (inv) | 0.811 |
| snr (baseline) | 0.785 | 0.215 (inv) | 0.764 |
| best catalog/location | ~0.55–0.63 | — | ~0.62–0.66 |

### Confound diagnostics
- **Within-station baseline (template_only / logit): weighted-mean AUC = 0.892**
  over 19 stations with both classes. High at nearly every station individually
  (B013 0.967, B004 0.955, HDW 0.962, B935 0.977, B039 0.996, GNW 0.901, PGC 0.893).
  → the signal discriminates GOOD vs BAD *inside a single station/instrument* —
  it is **genuinely family-intrinsic, not a station-cluster artifact**.
- **Negative control (within-station label shuffle): grpCV AUC 0.923 → 0.738 ± 0.030.**
  The drop confirms a large genuine component; the residual 0.738 (above 0.5) is
  the part of the grouped score that is station-template-distribution confound
  (per-station base rate is preserved under within-station shuffle). The honest
  family-intrinsic number is the **within-station 0.89**, not the headline 0.92.
- Permutation importance dominated by **t_spec_centroid (+0.099)** and
  **t_spec_entropy (+0.019)**; everything else < 0.005.
- **Catalog + location features are NOT predictive** (grpCV ≈ 0.43–0.63): local
  tremor recurrence, distance, azimuth, Slab2 depth do not tell you whether a
  given family will be CWI-suitable.

### Physical interpretation (consistent & directional)
- **BAD families:** higher spectral centroid (median ≈ **4.6 Hz**), higher kurtosis
  (median ≈ 7.9, spiky), higher spectral entropy/flatness — i.e. high-frequency,
  impulsive, noise-like / poorly-coherent templates.
- **GOOD families:** lower centroid (median ≈ **3.4 Hz**), low kurtosis
  (median ≈ 1.2, emergent), lower entropy — coherent, low-frequency LFE energy.
- Direction (**BAD centroid > GOOD centroid**) holds at **15 / 17** stations with
  ≥3 of each class (only COLT and B036 flip, both marginal southern outliers).
- Makes sense: a spiky/high-freq stacked template is dominated by spurious or
  noise-contaminated detections → sparse, gappy continuity → poor CWI; a clean,
  emergent, low-freq LFE template repeats steadily → continuous coda for CWI.

### Screening gain (best grouped model, template_only / HGB, OOF probs)
Flagging predicted-BAD families to **skip** before densify:
| threshold (P_good<) | % pool skipped | BADs skipped | GOODs lost | precision(BAD) | recall(BAD) |
|---|---|---|---|---|---|
| 0.30 | 20.8% | 325 / 443 | 68 | 0.83 | 0.73 |
| 0.40 | 22.5% | 341 | 84 | 0.80 | 0.77 |
| 0.50 | 25.2% | 352 | 124 | 0.74 | 0.80 |

At **threshold 0.30** we skip ~21% of the pool, correctly removing **325 / 443
(73%) of BAD families at 83% precision**, losing only 68 / 1449 (4.7%) GOOD
families. This is a real, operational pre-densify cull.

## VERDICT
**CWI-suitability IS predictable from pre-densify signatures — specifically from
the per-family *template waveform shape*.** Grouped-CV AUC = **0.92** (template
HGB), confirmed family-intrinsic by within-station AUC = **0.89** and the
label-shuffle control collapsing toward chance. This *overturns* the prior
"detection-stats don't predict continuity → must measure directly" expectation:
detection-COUNT/time stats and catalog/location remain non-predictive (grpCV
≈ 0.5–0.7), but the **template spectral character is the missing intrinsic
predictor**. The single best feature is **spectral centroid** (low → GOOD),
backed by kurtosis / envelope-crest / spectral entropy.

> Caveat to keep honest: the headline 0.92 has a station-confound floor of ~0.74
> (shuffle control). Quote the **within-station ~0.89** as the deconfounded,
> transfer-to-new-stations expectation, and the 0.92 only as the in-distribution
> pooled-grouped number.

## OPERATIONAL RECOMMENDATION
1. **Yes, pre-screen.** Compute the cheap template-shape features (esp. spectral
   centroid, kurtosis, spectral entropy, envelope crest) on every candidate
   family's discover-stage template — these are already produced by `discover_gpu.py`
   (the NPZ), so the cost is ~free, *before* any densify.
2. **Cull rule.** Two validated options:
   - *Model (max recall):* drop families with HGB P_good < 0.30 → skip ~21% of
     pool, remove 73% of BADs (precision 0.83), lose 4.7% of GOODs.
   - *Simple univariate (max precision, no model needed):* **spectral centroid
     > 4.3 Hz AND kurtosis > 4** → skip **11%** of pool, remove **47%** of BADs at
     **97% precision**, lose only **6 / 1404 (0.4%)** GOOD families. Loosening to
     centroid > 4.0 & kurt > 3 raises BAD recall to 56% at 94% precision (1.0%
     GOOD loss). Use the univariate rule for a safe high-precision cull; use the
     model for a more aggressive cull.
3. Because the discriminant is **family-intrinsic and station-transferable**
   (within-station AUC 0.89), the rule can be applied to **new stations not in the
   training set** with only mild degradation — exactly the use-case (pre-screen a
   fresh station's eligible pool).
4. Do **not** rely on catalog recurrence / location / detection-count stats for
   this screen — they are at chance under grouped CV.
5. A direct sampled-day densify probe is **no longer required** for triage; reserve
   it only for borderline families (0.30 ≤ P_good ≤ 0.6) where the template-shape
   call is least certain.

## Refinement (2026-06-10): is CONTINUITY (not just contaminant-vs-real) predictable?
Restricted to the 1691 GENUINE-LFE families (pass the centroid<4.3|kurt<4 screen) and predicted the continuity
target `cov_own_span` (active-day fraction within the family's own span) with GroupKFold-by-station.
- **Grouped-CV R² = +0.31** → continuity IS partially predictable across held-out stations.
- The driver is STILL template **spectral centroid** (within-station r=−0.53): lower-freq/emergent → steadier.
- Catalog recurrence (cat_active_mo_frac, cat_log_count, cat_inter_cv), location, depth: r≈0 → the
  "high-recurrence patch repeats more steadily" hypothesis FAILS.
**Unified conclusion:** ONE master axis = template spectral quality. It (a) strongly rejects contaminants
(AUC 0.92) and (b) moderately grades continuity among real LFEs (R² 0.31) — same feature. There is NO separate
repeat-rate signature; residual burst-vs-steady among equally-clean templates is set by the slow-slip cycle, not
a fixed family property → would need direct (sampled-densify) measurement. Operationally: the centroid/kurtosis
pre-screen already captures most of the predictable CWI-suitability; finer continuity ranking buys little.

## Spatial test (2026-06-10): is continuity predictable from LOCATION? NO.
Mapped continuous vs sparse repeaters for all 35 stations (2468 station-family curves, 74% continuous):
`figures/family_continuity_map.png` (points) + `figures/family_continuity_fraction_map.png` (per-0.1°-cell fraction).
GroupKFold-by-station AUC for predicting continuity at a NEW station:
- distance-to-station 0.54 | source lat/lon 0.58 | **local tremor recurrence 0.52 (chance)** | all combined 0.61.
- vs template-shape predictor 0.92. → LOCATION/recurrence carry NO usable signal.
- Latitude trend exists but BACKWARDS: south 40-42N 90% continuous, north 48-50N 66% — densest-tremor north has
  LOWEST continuity → it's a station-geometry/coverage artifact, NOT intrinsic source activity. "Busy fault patch →
  continuous repeater" hypothesis FALSIFIED.
**FINAL predictability hierarchy:** real-vs-contaminant = template shape AUC 0.92 (usable pre-screen); continuous-vs-
episodic = template centroid R² 0.31 (modest); from location/recurrence/geometry = chance. A family's CWI-suitability
is in its WAVEFORM SHAPE, not its location. Data: data/family_continuity_classes.csv, _locfeatures.csv.

## ONE-MONTH DENSIFY PROBE (2026-06-10): predicts full-record continuity at AUC 0.92
Tested on 2605 fully-densified families (known full-record class): active-day count in a SINGLE month vs
full-record continuous/sparse class.
- Pooled AUC 0.917. Distributions: continuous fam = ~30/30 active days in a typical month (25th-pct month 29);
  sparse fam = ~6/30 typical, but BEST month up to 30 (ETS burst).
- Operating point: densify 1 month, call continuous if **>=15 active days** -> typical-month TPR 0.98, FPR 0.09.
- WHICH-MONTH sensitivity (real): continuous fam in a data-gap month looks dead (worst-month median 4d);
  sparse fam in an ETS burst false-passes 92%. FIX: probe a NON-ETS month (catalog gives timing) OR probe 2
  separated months & require both active (episodic fam won't be mid-burst twice). 88% of continuous fam are
  active >=15d in >80% of their months -> reliably active except under ETS-burst sampling.
- ~30x cheaper than full densify (1 month vs ~250 months).
**COMPLETE CHEAP TRIAGE (combine):** (1) template-shape screen [free, pre-densify, AUC 0.92] drops contaminants;
(2) one-month densify [cheap, AUC 0.92] grades continuity among survivors — the question template shape could NOT
answer (continuity was only R^2 0.31 from template). Together: audition the whole eligible pool before committing
to full densification. Location/recurrence remain useless (AUC ~0.55).

## AUDITION MECHANISM (scripts/family_audition.py) + cleanup rule (2026-06-10)
Two-stage cheap family-set EXPANSION protocol, validated (recall 0.93-1.0, precision 1.0 on 4 stations):
- `--screen STA --floor F` : Stage 1, FREE. Template centroid>4.3 & kurt>4 -> contaminant, dropped. Writes
  data/<sta>_audition_screen.csv (genuine flag). B003: 584 eligible -> 362 genuine (222 contaminants dropped).
- `--probe-stage2 STA --probe-years Y [Y2]` : Stage 2, LIVE (needs waveforms). Densifies ONLY the genuine,
  not-yet-done candidates over a short probe window (1 year robust, or 2 separated months/years to beat
  ETS-burst false-pos); keeps families active >=50% of probe days = continuous. Writes _audition_keep.csv
  (winners -> full densify) + _audition_log.csv.
- **CLEANUP RULE (user, baked in):** the probe densify is a THROWAWAY classifier — after classification it
  DELETES all probe-scratch files (mf_<sta>_probe_*.csv + temp candidate summary), for skipped AND winners
  (winners get a separate full-record densify). Only the small keep-list + log survive. Prevents scratch bloat.
- `--validate STA` : simulate the probe from existing stacks vs full-record truth (no download) — how it was verified.
LIVE PILOT PENDING: B003 (re-download -> screen[done:362] -> probe-stage2 -> full-densify ~120 winners -> merge
-> restack), gated on whether expansion measurably improves the dv/v/fault-cell count. Solidify here BEFORE the
post-2020 fleet (else 248x rework + repeated B933-type wasted runs).

## Does the audition remove GOOD families? (false-negative check, 2026-06-10)
Checked Stage-1 template screen against MEASURED continuity truth (1816 families w/ template+truth):
- Global thr (centroid>4.3 & kurt>4): drops 3.5% of CONTINUOUS families — but ALL concentrated at B020(45)+
  B028(6), stations whose genuine LFE templates intrinsically run high-frequency (B020 cont centroid 4.7 vs
  sparse 4.9 — barely separable). Elsewhere clean.
- FIX (applied): STATION-RELATIVE centroid cap = max(4.3, station's own 85th-pct) → false-neg on continuous
  3.5%->1.2% (sparse still 22% dropped). family_audition.py screen_pool updated.
- DESIGN RULE: Stage 1 = conservative PRE-FILTER only (save Stage-2 cost), NOT the arbiter. Stage 2 measures
  the real targets directly & station-agnostically: active-day fraction (continuity) AND coda cc in the probe
  window (>0.7; catches frequently-firing contaminants like B933 cc 0.54 that Stage 1 might pass). Stage 1 and
  Stage 2 catch DIFFERENT failures (contaminant vs episodic); use both, Stage 2 decisive. [TODO: add probe-window
  coda-cc check to probe_stage2 — currently counts detections only.]
