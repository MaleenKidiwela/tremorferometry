# The Family Trust Test — deciding which LFE families to believe, per station

*Design doc, 2026-06-10. Motivated by the B003/B018 contamination investigation (notes/2026-06-10 §16-18):
detection lists can be heavily noise-contaminated while the stacked dv/v is real (B018), or genuinely
contaminated (B003 expansion). Detection-count statistics alone cannot tell these apart — only
product-level waveform tests can. This battery decides, per (station, family): TRUSTED / SITE-CARRIER /
CONTAMINATED / UNDETERMINED.*

---

## 1. The problem it solves

Our matched-filter detector (2-s Z-only templates, fixed cc≥0.8 ≈ 3σ) admits two kinds of detections:
real repeats of a deep LFE source, and ambient/cultural noise windows that happen to correlate with the
template. A family's detection list is therefore a MIX with unknown proportions. The dv/v is measured on
daily STACKS of those windows, where real repeating signal reinforces coherently and noise averages toward
the station's ambient-noise autocorrelation (ACF). So a family's dv/v can be:
- carried by real LFE coda → genuine deep-fault-sensitive measurement (keep for tomography);
- carried by the noise ACF → a real but SHALLOW signal (single-station autocorr interferometry of the
  near-surface) mislabeled as a deep source (keep only as a site channel);
- garbage (no coherent content at all) → exclude.

Three hard lessons shape the design:
1. **Detection-list statistics cannot condemn the product.** B018's lists look awful (day/night 2.5, ~90
   det/day, cc hugging 0.83, tremor-corr ~0.05 — partly because the top-100/day cap censors everything),
   yet its stacks contain unambiguous real repeating coda. Count stats are FLAGS, never verdicts.
2. **Anything inside the selection window is coherent by construction.** Detections are chosen to match the
   template at −1..+1 s, so coherence there proves nothing — even pure noise matches stack into the template
   shape there. Validity must be tested OUTSIDE the selection support (the 2–4 s coda), where 2–8 Hz noise
   decorrelates in ~0.2 s and selection cannot manufacture coherence.
3. **Every statistic needs an empirical null.** Thresholds tuned by eye drift; the same pipeline run on
   windows that CANNOT contain real repeats (random times; time-reversed templates) gives each station its
   own honest null distribution.

## 2. The battery

### Tier 0 — free flags (existing csv/npz, minutes per station; flags only, never verdicts)
- **Template shape**: spectral centroid + kurtosis of the discovery template. Spiky/high-frequency
  (kurt>4 & centroid above the station-health-gated cut) = the impulsive non-LFE class (B933). Validated
  AUC 0.92 — but it does NOT catch the noise-ACF class (emergent templates look "GOOD").
- **cc-distribution shape**: noise-matched lists hug the threshold (median ≈0.83–0.87, exponential tail —
  the analytic prediction for cc≥0.8 at N_eff≈24); real repeats add a separated high-cc mode.
- **Cap saturation** (fraction of family-days at the top-100 cap) and **day/night ratio** (solar-diurnal
  detections = cultural noise admixture). High values mean "the detection list is polluted", NOT "the dv/v
  is wrong".

### Tier 1 — the decisive waveform tests (needs the station's raw waveforms once)

**T1a. Stack-vs-random coda ratio (the core test).**
For each family: stack N≈300 randomly sampled detections (bandpass 2–8 Hz, normalized); separately stack
N random times from the same record; repeat the random stack K≈10 times to get a null distribution.
Score = coda-window (2–4 s, OUTSIDE the selection support) amplitude of the detection stack, in σ-units of
the random-stack null. Also recorded: the S-window ratio (information, not verdict — selection inflates it).
- Real repeating family: coda ≫ null (B018 anchor case: 3.2× amplitude).
- Noise-matched family: coda ≈ the ACF tail ≈ null level.
Why it works: selection conditions only on −1..+1 s; coherent structure seconds later can only come from a
genuinely repeating source-side wavefield.

**T1b. Source-independence split (day/night stack agreement).**
Split the family's detections into local-day and local-night, stack each, compare the two stacks' coda
(2–4 s) cross-correlation, and (if both support a dv/v) the agreement of day-only vs night-only dv/v series.
The earth does not know the time of day: medium-carried signal → the arms MATCH; noise-source-carried
signal → day and night noise differ → the arms diverge. (Orthogonal split variants: odd/even weeks.)

**T1c. Time-reversed-template null families (the built-in calibration).**
Run the SAME matched filter on ~30 sampled days using the station's templates TIME-REVERSED (identical
spectrum and N_eff; cannot match real repeats). The resulting "families" are guaranteed-fake; push them
through T1a/T1b. They define the station's empirical fail distribution — every threshold in the battery is
set from these, not by hand. (Also yields a per-family false-detection-rate estimate: reversed/forward rate.)

### Tier 2 — dv/v-level attribution (needs the dv/v products + one template-free ACF run)

**T2a. Shallow-share slope.** Compute the station's template-free Z-autocorrelation dv/v
(`scripts/autocorr_dvv.py` — no templates, no catalog). Regress each family's dv/v on it (deseasoned,
slope not correlation — correlation is uninformative because the stable deep fault makes even genuine
families correlate with the site signal).
- slope ≈ ⅓ → the kernel-predicted receiver-lobe share of a GENUINE two-lobed family (B018: 0.34);
- slope ≈ 1 → the family IS the shallow monitor (SITE-CARRIER);
- calibrate the dividing line per station with the T1c null families' slopes.

**T2b. Rolling-stack coda coherence.** The family's existing cc_max in the 2–4 s product. A collapse
(≈0.5–0.6 while the station's good families sit ≥0.8) = no usable coda (the B933 failure), regardless of
everything else.

## 3. Scoring → tiers

Per family, in order:
1. T1a coda ratio below the reversed-template null (or T2b collapse) → **CONTAMINATED** (exclude from all
   velocity products; detection list may still interest the noise channel).
2. T1a passes but T2a slope ≈ 1 and T1b arms disagree → **SITE-CARRIER** (real shallow signal; keep as an
   explicitly-labeled near-surface/site channel; NEVER attach a deep-fault kernel to it).
3. T1a passes, T1b arms agree, T2a slope ≈ ⅓ → **TRUSTED** (fault-tomography eligible; weight by its
   T1a σ-score and detection independence rather than raw counts).
4. Anything mixed → **UNDETERMINED** (usable for station-level products with down-weighting; re-test with
   more samples before tomography use).
All scores are saved continuous (σ-units, slopes, agreement cc), so the inversion can weight instead of
binarize. Nothing is deleted — class labels go in `data/family_trust_<sta>.csv`.

## 4. Validation before rollout (non-negotiable)
Run the full battery first on **B018** (waveforms on disk; ground truth established): its coverage-selected
families must score TRUSTED; its time-reversed nulls must score CONTAMINATED; disagreements get inspected
individually before the battery is believed. Second target: **B003** (re-download), where we EXPECT the
expansion families to fail and the verdict on the original 81 is genuinely open — both outcomes informative.
Only then station-by-station rollout (one re-download at a time, janitor-cleaned).

## 5. Costs
- Tier 0: free, minutes/station.
- Tier 1: the station's waveforms on disk (re-download ~30–50 G where janitored, a few hours, resumable);
  the tests themselves read ~300 windows × (families + nulls) ≈ 1–2 h/station single-pass (group window
  reads by day-file); reversed-template null ≈ minutes of GPU on 30 days.
- Tier 2: ACF dv/v ≈ 20 min/station (already built); regressions free.

## 6. What the battery does NOT decide
- It does not prove a TRUSTED family's source is ON the plate interface (location comes from the catalog
  seeding); cross-station detection coincidence remains the gold-standard source test for shared patches.
- It does not rescue the detection lists for rate-based science (the cap censoring stands).
- It does not measure deep velocity CHANGE — it only certifies which families are valid instruments. The
  deep signal itself remains whatever the certified families collectively say (currently: stable interface).
