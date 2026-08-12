# Broadband / land-station dv/v — work log (started 2026-07-14)

**GOAL:** extend the LFE-coda dv/v pipeline from PB boreholes (27-strong fleet DONE) to BROADBAND / land
stations, producing RELIABLE dv/v. Pilot = **PGC** (CN.PGC, BHZ 40 Hz, ~290 m from borehole B011). The
fingerprint/picker must be VALIDATED and TRAINED for broadbands (user requirement) — not transferred from
the borehole picker. Advisor: Merlin (consulted 2026-07-14).

## KEY CONSTRAINTS / DECISIONS
- **Boreholes are done; do NOT touch the borehole pipeline or the B011 picker.** Broadbands get a SEPARATE
  fs-portable feature set + a fresh broadband picker.
- **Train a FRESH broadband picker** (Merlin), not a recalibrated borehole one — surface microseism <2 Hz and
  HF attenuation shift absolute features. Run the 3-way comparison (raw / recalibrated / fresh) to quantify.
- **Ground truth = co-located B011→PGC labels (primary).** B011's causality-certified families give precise,
  abundant detection times that label PGC's surface record — no circularity (labels from borehole waveforms,
  features from PGC). Lin catalog = held-out 2nd source. Ducellier catalog = out-of-region VALIDATION gate
  only (its cc~0.12 network times are seconds-jittery; not for training).
- **fs-portability:** PGC pilot done at 40 Hz native (internally consistent). BEFORE rolling to HH-100 Hz
  broadbands, fix feature_defs (features integrate to Nyquist → fs-dependent): common band ~2–18 Hz, slope
  2–15 Hz, hf/bf on 8–16(18) Hz, decimate all to 40 Hz. (Merlin flaw #3.)
- **Scalable strategy:** any broadband within ~20–30 km of a certified PB borehole can be labeled this same
  co-located way — that, not Ducellier, is the fleet-wide broadband training approach.

## TESTS & RESULTS
### T1 — Stack test: does the LFE coda SURVIVE at the surface? (Merlin step 1, DECISIVE) — ✅ PASS
`scripts/pgc_stack_test.py`. Stack PGC BHZ at B011 certified-family high-cc detection times (2010–2017),
measure causality ratio RMS(coda 2–4 s)/RMS(mirror −2..0 s).
- ⚠ **Method bug caught before reporting:** first run gave 0% — artifact of anchoring coda/mirror on the
  *found peak*, so the mirror window swallowed the direct-S energy. Fixed to FIXED-ORIGIN windows (t=0 =
  detection origin, S ~ +1 s), matching `finalize_causality`.
- **RESULT (corrected): 70/119 certified families = 59% show a causal coda at PGC's surface** (ratio p50
  1.65, p75 2.11, p90 2.83, max 4.36). Arrival lag +1.12 s = matches B011's S-anchor → co-location physically
  confirmed. **Verdict: broadband coda dv/v is VIABLE** (≥30% bar cleared, better than borehole ~1/3 rate).
- LESSON: causality-ratio tests must use fixed-origin windows; peak-anchoring is biased by the direct S.
- Prior: RESULTS.md R5 already had LFE-vs-RAND AUC 0.978 at PGC (window detection works); T1 answers the
  open coda/stack question that dv/v actually needs.

### T2 — Broadband picker training (IN PROGRESS)
⚠ **Label bug caught on first extraction (fixed):** (1) positives at cc≥0.85 are noise-dominated at the
window level (matched-filter fires ~76% of all minutes across families) → borderline labels; (2) TNOISE was
placed in the rare no-trigger GAPS = quiet times ≈ RAND (the EASY negative), not the hard in-tremor rumble,
and only 503 survived. **Fix:** positives = cc≥0.92 (high-confidence, cleaner individual windows); TNOISE =
random times inside PNSN tremor windows near B011 (34,110 tremor detections 2010–2017 = genuine in-tremor
noise), ≥60 s from a positive. Re-running.
Building training set via `lfe_features/extract_colocated_bb.py` → `feat_coloc_pgc.parquet`:
- COLOC positives (B011 cert high-cc → PGC windows) + TNOISE (in-tremor noise, the hard negative that
  decides deployability). Existing on disk: feat_ref_pgc (Lin+RAND), feat_eq_pgc (EQ).
- Next: train fresh RF+isotonic picker, grouped CV (by family AND year), SNR/amplitude-drop control.
- Acceptance (Merlin): vs RAND ≥0.95, vs TNOISE ≥0.75, amplitude-control AUC drop ≤0.02.
  Failure: vs TNOISE <0.65 → picker can't score during tremor, broadband discovery noise-flooded → stop.

**★ FINDING — co-located B011 positives are TOO WEAK per-window for training.** Trained on COLOC (B011
cc≥0.92 → PGC) the picker got only **0.52 vs RAND** (coin-flip). Diagnostic: COLOC positive features look
almost like RAND noise (2–4 Hz energy 0.003 vs Lin 0.010; hf_ratio 1.25 vs Lin 0.72). Co-located detections
are precise in TIME but many are WEAK LFEs invisible in a single surface window — they only emerge STACKED
(hence T1 passed at 59% but per-window fails). **So co-location validates the coda/STACK, not per-window
detection.** For the picker use STRONG per-window-visible LFEs = the Lin network catalog. (Fleet implication:
broadbands without a Lin-type catalog need a STRENGTH filter on the co-located labels, not raw cc≥0.92.)

**★ RESULT — fresh PGC broadband picker (Lin positives, GroupKFold by year), `train_broadband_picker.py`
→ models/picker_broadband_pgc.joblib:** vs RAND **0.932** · vs EQ **0.896** · **vs TNOISE 0.769** (hard
in-tremor test — PASSES ≥0.75, well above the 0.65 fail line) · amplitude-control drop **0.002** (fingerprint
is shape, not loudness). vs RAND (0.932) just under Merlin's 0.95 because grouped-by-YEAR CV is stricter than
R5 and only 940 Lin positives. **VERDICT: broadband LFE detection VIABLE at PGC.** Next: full PGC pipeline
run (discover with this picker → densify → causality → dv/v).

### ★★ MERLIN RESULTS-VALIDATION AUDIT (2026-07-14) — caught real issues; NO-GO on full run yet
- **FATAL (fixed): the grouped-by-YEAR CV was fake.** `train_broadband_picker.py prep()` assigned RANDOM
  years (inputs lacked a `year` col) → the CV was a plain random split, so 0.932/0.769 are the LENIENT
  numbers and vs RAND MISSES the 0.95 bar under the easy protocol. Fixed prep() to derive year from `time`
  (unit='s'!) and error on missing — but see next.
- **MAJOR: positives are ONE ETS episode.** Lin curated positives near PGC = 940, ~98% Aug–Sep 2010 (PGC
  waveforms on disk start 2010; the Lin near-PGC catalog is mostly 2005–2011, 5,633 in 2005). So cross-
  episode generalization (what production 2018–2026 needs) is UNTESTED, and year-grouped CV is impossible
  until we **download PGC BHZ 2005–2009** (~10k positives across ≥5 episodes). This is THE picker fix.
- **RESULT 1 (stack test) independently VERIFIED sound** (59%, tpk 1.03–1.22 s). But it's 8-YEAR stacks;
  dv/v uses **30-CAL-DAY** stacks → the decisive unchecked question (running: `pgc_30day_stack_test.py`).
- **DECISION A (Lin positives) sound for PGC**, but RESULT 2's "weak-but-real" interpretation is ASSERTED,
  not proven (competing: most B011 cc≥0.92 are FALSE noise triggers). Matters for the fleet strategy; needs
  the coincidence/stack-convergence diagnostic. Fleet claim (worklog "any broadband w/in 20–30 km of a
  borehole labeled this way") is FALSIFIED per-window — co-location is a STACK/coda transplant tool, not a
  per-window labeling tool. Broadbands outside Lin coverage → picker-transfer + causality (borehole precedent).
- **MINOR:** amplitude-control is vacuous (all features already scale-invariant — cite "shape by design",
  not as an empirical control). Document TNOISE exclusion rule (it IS reproducible via extract_tremor_noise.py).
- **Reassurance:** the 30-station borehole fleet ran on B011's picker (vs TNOISE 0.809) with causality as the
  referee; a PGC picker at 0.769 is the same class. Picker quality is NOT the risk — validation honesty +
  the 30-day-stack question are. **Merlin: NO-GO on the full pipeline until the 30-day test + honest re-report
  + 2005–2009 download clear.**

### T3 — 30-DAY rolling-stack coda test (Merlin's DECISIVE dv/v gate) — ✅ PASS (strong)
`scripts/pgc_30day_stack_test.py`. For the 70 T1 families, bin B011 cc≥0.85 times into 30-cal-day blocks,
stack PGC per block, per block compute causality ratio + coda-cc vs the family's all-time coda.
- **RESULT: 43 families have ≥50% of their 30-day blocks passing (ratio>1.5 AND coda_cc≥0.6).** 6142 blocks,
  58% pass; median block ratio 1.61, coda_cc **0.86**, n 300. Acceptance was ≥20 → **PGC dv/v VIABLE
  YEAR-ROUND** (>2× the bar; coda is strong AND temporally stable at the 30-day resolution dv/v uses).
- This answers the question T1 (8-yr stacks) couldn't: dv/v lives on 30-day stacks, and they carry a usable,
  stable coda at PGC's surface for 43 families → clears the ≥20-certified inclusion rule handily.

## STATUS: physics GREEN (T1+T3). Remaining before production PGC dv/v (Merlin):
1. Honest picker re-report + **download PGC BHZ 2005–2009** → true cross-episode GroupKFold CV (the picker is
   a candidate pre-ranker; causality is the dv/v referee, so this is about validation honesty not blocking).
2. Download PGC 2018–2026 (production dv/v span; on-disk is only 2010–2017).
3. Full PGC pipeline run: discover (broadband picker) → densify → 30-day stacks → causality → dv/v → mirror
   clean → metadata QC (2013 glitch) → ≥20-certified inclusion. Then generalize to other broadbands.

### ★★ PGC TWO-INSTRUMENT DECISION (Merlin 2026-07-14) — decimate, but TWO SEGMENTS
PGC archive reality (probed IRIS): downloadable waveforms only **BHZ 40 Hz 2010–Aug 2017**, then **HHZ 100 Hz
Aug 2017–2026** (sensor swap, NO overlap; nothing pre-2010). Options weighed: (A) 2nd 100 Hz picker, (B)
decimate HHZ→40 Hz + reuse the 40 Hz picker.
- **VERDICT = B (decimate to 40 Hz, ONE picker).** Everything the method uses is ≤16 Hz (mostly 2–8 Hz),
  well under 40 Hz Nyquist; polyphase 100→40 is transparent + is ALREADY the fleet standard (`densify_gnw_gpu.py`,
  `build_long_window_3comp.py`). **Option A is UNEXECUTABLE** — no strong LFE catalog covers 2017–2026 to
  train a 100 Hz picker (Lin is 2005–2011 only; co-located labels are per-window-useless).
- **BUT PGC = TWO dv/v SEGMENTS, offset UNJOINED.** No BHZ/HHZ overlap → the Aug-2017 join is empirically
  uncalibratable; an instrument step and a velocity step at the same date are unfalsifiably confounded.
  ⚠ BUG TRAP: `dvv_roll30cal.py` line 53 builds the stretch reference as the ALL-TIME mean → run across the
  swap it manufactures a fake step. **Reference must be split per era (run dv/v on per-era npz).** Mirror-
  correct per era too.
- **Decisive instrumental-vs-real test:** B011 (290 m away, unchanged instrument) must show the Aug-2017 step
  too if it's a REAL velocity change; if PGC eras disagree and B011 is smooth → offset is instrumental.
- **★ RESPONSE CHECK PASSED (Merlin step 1):** BHZ & HHZ both flat-to-velocity across 2–8 Hz (norm |resp| at
  2/4/8 Hz = 1.000/0.997/0.985 BHZ, 0.985/0.992/1.000 HHZ) → **raw-bandpass-per-era valid, NO deconvolution**
  (deconvolution is itself a suspected artifact source — Nisqually memory). Eras are shape-comparable.
- **Single-2010-episode picker caveat is now PERMANENT** (no pre-2010 data → Merlin's 2005–2009 CV plan is
  dead). Lean on causality (picker-independent referee); validate HHZ era via B011 co-located T1/T3 stack
  tests (non-circular). PGC picker = candidate pre-ranker only.
- **PLAN:** response ✅ → download HHZ 2017–2026 (RUNNING, decimate 100→40) → T1+T3 on HHZ era (decisive gate
  for segment 2; ≥30% causal / ≥20 families 30-day) → one densify whole record w/ BHZ-era templates → causality
  certify PER ERA → dv/v PER-ERA references → report 2 segments + B011 instrumental-vs-real figure.

### T4 — HHZ-era (2018-2026) segment-2 physics gate — ✅ PASS (stronger than BHZ era)
`scripts/pgc_hhz_stack_test.py` (decimates 100->40). B011 certified cc>=0.85 times 2018-2026, stack PGC HHZ.
- **T1: 104/119 families (87%) causal ratio>1.5** (BHZ era 59%). **T3: 54 families >=50% of 30-day blocks
  pass** (BHZ era 43). Median block ratio 1.51, coda_cc 0.84. Both >> bars.
- **VERDICT: PGC = TWO viable dv/v segments, 16-year record** (2010-2017 BHZ: 43 reliable families; 2018-2026
  HHZ: 54). Decimation 100->40 clean; newer HHZ era is actually cleaner (better modern data). Download of
  2017-2026 HHZ (2931 files) done.
- NOTE (self-critique): let the HHZ download sit idle after it finished before running this — should chain
  next step to completion, not wait.

## STATUS: BOTH PGC segments physics-GREEN. Next → produce the actual dv/v time series (per-era reference,
reliable families only, mirror-corrected), then the B011 instrumental-vs-real comparison across Aug-2017.

## SEQUENCE (Merlin steps)
1. ✅ Stack test (coda survives — 59%).  2. Fix feature band (deferred to fleet rollout).  3. Train fresh
PGC picker on co-located labels + validate vs RAND/EQ/TNOISE.  4. 3-way picker comparison.  5. Held-out 2nd
CN broadband (Lin-labeled).  6. Ducellier gate.  7. Full PGC pipeline run (discover→densify→causality→dv/v),
all standard gates + ≥20-certified inclusion rule. Then generalize to all broadbands.

## PGC facts
BHZ/BHN/BHE 40 Hz, 2009–2017 on disk (46 GB, 2586 days); production dv/v will need a 2018–2026 download.
Known PGC 2013 glitch year — apply despike-mad 8, check traces/day. Co-located borehole B011: 119 certified
families, 2010–2026 mf detections = the labels.

## 2026-07-15 — PGC dv/v COMPLETE (two-segment deliverable)
Pipeline `scripts/pgc_dvv_pipeline.sh` finished. PGC (~290 m from borehole B011) built its own daily
stacks at B011 reliable-family detection times, causality-certified AT PGC, coda stretch 2-4s, per-era refs.

| Segment | Instr | Years | Reliable families (caus-cert) | dv/v std |
|---|---|---|---|---|
| Era 1 | BHZ 40Hz | 2010-2017 | 161 | 0.481% |
| Era 2 | HHZ 100->40Hz | 2018-2026 | 201 | 0.370% |
| Total | | 16 yr | 362 | |

Decisions/results:
- PER-ERA references held (no forced jump across Aug-2017 sensor swap) — Merlin's critical requirement met.
- HHZ era cleaner (std 0.37 vs 0.48, 201 vs 161) -> decimate-100->40 path is sound; two unjoined segments.
- Confirms stack->causality->per-era-dv/v recipe works on SURFACE BROADBAND, not just boreholes.
- Both eras clear the >=20-family inclusion bar comfortably. PGC = viable broadband station.
Outputs: data/daily_dvv_PGC{bhz,hhz}_Z_2to4.csv ; data/pgc{bhz,hhz}_causality_cert.csv ; data/pgc{bhz,hhz}_3comp_summary.json
Next (Phase 0): (2) B011 instrumental-vs-real check across Aug-2017; (3) validate PICKER-DISCOVERY path on CLRS (co-located cross-check).

## 2026-07-15 — B011 (borehole) ↔ PGC (broadband) CO-LOCATED VALIDATION  ★ PROGRAM VALIDATED
The load-bearing test: does the broadband dv/v reproduce the trusted co-located borehole? If the one
co-located pair agrees, the whole broadband program is trustworthy. PGC (surface BB) and B011 (borehole)
are ~290 m apart, share the SAME LFE family catalog (PGC built at B011's detection times).

### Merlin caught a fatal flaw in my first plan (consulted BEFORE coding — memory rule paid off)
My original plan = correlate raw PGC vs B011 dv/v. Merlin: INVALID. At 290 m they share the same noise
field (microseism/weather), so raw dv/v corr is ~guaranteed high EVEN IF PGC's coda were pure noise
(this project's own "raw dv/v = 65% artifact" finding). Verdict must live in MIRROR-CORRECTED, deseasoned
residuals, referenced against a PGC-mirror-vs-B011-coda NULL (same triggers/days/noise, zero LFE energy).
=> had to BUILD PGC mirror dv/v first (scripts/pgc_mirror_build.sh: build_mirror_npz + dvv_roll30cal per era).
Merlin also fixed my wrong B011 cert path -> data/b011_fwd_vs_rev_coda.csv (ratio>1.5 = 119 fam).

### Method (scripts/b011_pgc_validation.py, replicates cm()/des()/beta from plot_mirror_corrected_dvv.py)
Per era: F_e = PGC-era-cert ∩ B011-cert. Build cm() cross-family median (cc>=0.6, >=3 fam/day, 15d roll)
for PGC/B011 × coda/mirror. R = des(coda) − beta·des(mirror). S1 = Spearman(R_PGC,R_B011). N1 = mirror null
= Spearman(des(PGC mirror), R_B011). Moving-block bootstrap (block=60d, ~55-65 effective samples/era).
PASS(e): S1 CI>0 AND (S1-N1) CI>0 AND slope in [0.3,3]; INCONCLUSIVE if S1<0.2. Program VALIDATED iff both PASS.

### RESULTS — both eras PASS, strongly
| metric | BHZ 2010-2017 | HHZ 2018-2026 |
| raw coda r (confounded)   | 0.606 | 0.808 |
| mirror null N1 (~0 target) | 0.038 | 0.041 |
| corrected S1 (real signal) | 0.379 [0.178,0.550] | 0.434 [0.269,0.591] |
| S1 - N1 (excludes 0)       | 0.341 [0.043,0.630] | 0.393 [0.098,0.663] |
| slope OLS (in [0.3,3])     | 0.62 [0.34,0.90] | 0.60 [0.41,0.80] |
| families / common days     | 94 / 2546 | 114 / 3020 |
| VERDICT | PASS | PASS |

KEY: raw 0.6-0.8 collapses to N1≈0.04 under the mirror null -> the agreement is NOT shared weather; the
~0.4 that survives is genuine LFE-coda dv/v the broadband reproduces from the borehole. The raw→corrected
drop (~0.6→0.38) quantifies how much naive broadband-borehole "agreement" is artifact — a number to reuse
for every future broadband station.

### PER-FAMILY (matched vs mismatched) — the tomography-critical result
BHZ: matched med r=0.466 (n=94) vs mismatched 0.156 (n=1879), Δ=0.310, p=0.0005
HHZ: matched med r=0.552 (n=114) vs mismatched 0.176 (n=2280), Δ=0.376, p=0.0005
=> broadband carries PER-FAMILY velocity information (not just station-wide common mode) -> usable in the
4-D inversion, not just as a station-median series.

### SEAM (Aug-2017 sensor swap) — BUG CAUGHT + honest verdict = UNMEASURABLE (no step claimed)
BUG (self-caught): first run printed "step=+nan%, p=0.006, REAL step". Spurious — obs came out NaN and NaN
comparisons collapse the null p to its floor (1/180). Fixed: guard NaN -> report UNMEASURABLE; added
diagnostics + a deseasoned-coda fallback.
Cause: B011's MIRROR product has 0 days in May-Nov 2017 (coda has 209) — the mirror npz predates the 2017
coda gap-fix -> corrected residual empty across the seam. Real data gap, not a logic error.
Fallback (deseasoned coda, NOT mirror-corrected): step=-0.345%, p=0.055 — MARGINAL and artifact-prone
(mirror-uncorrected -> could be the 65% noise). => NOT claiming a real site step.
HONEST STATEMENT: (a) PGC's cross-seam dv/v LEVEL is unmeasurable by construction (per-era refs); (b) the
swap is validated as HARMLESS to within-era dv/v (both eras PASS independently); (c) a real site step at
Aug-2017 cannot be established with clean (mirror-corrected) methods given B011's mirror 2017 gap.
DATA-QUALITY TODO: rebuild B011 mirror to fill the 2017 gap if the seam step ever becomes important.

### IMPLICATION
Broadband LFE-coda dv/v is VALIDATED against the borehole ground truth (station-level AND per-family, both
instrument eras). This justifies the full broadband/land fleet rollout (239 stations) toward the spatial-
resolution goal. PGC is a bona-fide 16-yr broadband dv/v station.
Outputs: data/b011_pgc_validation.json ; logs/b011_pgc_validation.log ; scripts/b011_pgc_validation.py

## 2026-07-15 — 40 Hz consistency fix (score_candidates) — FLEET-CRITICAL
User caught: I launched CLRS detection at --fs 100 (native HHZ). The whole pipeline standard is 40 Hz, and
crucially the BROADBAND PICKER was trained at 40 Hz (extract_colocated_bb.py FS=40.0, skips non-40Hz traces).
Bug found: lfe_features/score_candidates.py used the trace's NATIVE fs (line 40) -> on 100 Hz HHZ it would
extract features at 100 Hz (Nyquist 50 vs 20) = wrong spectral slope/HF-depletion/2-8Hz-energy -> garbage
P(LFE) from a 40Hz-trained picker. discover_gpu/densify/stacks/dvv were already 40 Hz; discover_nllb resamples
to --fs internally; only the scorer didn't.
FIX: added --target-fs to score_candidates.py; resamples each day (resample_poly) to 40 Hz BEFORE feature
extraction. Fleet rule: broadband scoring ALWAYS --target-fs 40 (boreholes omit = native, back-compat).
CLRS redone: detection --fs 40, picker=picker_broadband_pgc copied to tremor_picker_clrs, score --target-fs 40.

## 2026-07-15 — CLRS picker-DISCOVERY validation (fleet-gating) — Merlin two-arm test
Validating the discovery path (finding LFE families with NO borehole answer key) on CLRS (CN broadband, ~40 m
from borehole B926 = the answer key). Pipeline: PNSN-driven detection -> broadband picker P(LFE) -> adaptive
top-N -> GPU cluster -> densify -> causality.

### Merlin caught my flawed metric (shift-null)
My "AUC 0.66, monotonic 8-10x enrichment" evidence for picker transfer was UNSOUND. Merlin ran a +120s TIME-
SHIFT null -> AUC UNCHANGED (0.658). The picker score tracks TREMOR-EPISODE INTENSITY, not per-event LFE
identity. RULE now: every coincidence claim in this project needs a shifted-time null or it's just measuring
the tremor catalog. Also: AUC 0.66 is NOT proof of degradation (in-domain 0.93-0.98 was Lin-vs-RAND, a
different protocol). The coincidence-AUC has no authority over the fleet decision -> run end-to-end instead.
Also caught a live trap: clrs_cand_filtered.parquet was EMPTY (scored --thr 0.5 but max P=0.33); adaptive
top-N is rank-based, re-run needed.

### Two-arm decisive test (picker vs picker-blind control)
Adaptive top-30k (thr 0.201, 31786 cand). Control arm = top-30k by SNR (picker-blind, 1% overlap). Both
through discover_gpu @40Hz:
| arm | families | members | B926-coincidence real | +120s null |
| picker  | 48 | 155 | 5.00% | 0.00% |
| control | 92 | 676 | 1.59% | 0.23% (7x) |
=> PICKER ADDS PURITY not count: picker families are REAL LFEs (5% coincide w/ confident B926 vs 0% shift-null
= clean pass of the mandatory shift-null); control finds MORE clusters but ~3x LESS pure (many non-LFE
repeaters). Picker families span 48.2-49.0N (med 48.75) = B926's zone (48.82N). Picker arm PASSES Merlin
step-2 (>=30 families, coincidence >=2x null).

### DECISIVE GATE (running): densify 48 picker families -> stacks -> dv/v -> causality
GATE = fleet inclusion rule: >=20 causality-certified AND >=15% survival (NOT B926-parity; B926 itself = 110
certified from 530 discovered). If pass -> discovery path validated, launch fleet. If fail (B045 signature:
many families, single-digit certified) -> densify control to localize blame (picker vs site).
Scripts: scripts/clrs_two_arm.sh, scripts/clrs_densify_causality.sh. Data: clrs_{picker,control}_families.*

## 2026-07-15 — CLRS yield: PARAMETER-DRIFT BUG (min-years), not broadband physics (Merlin)
Gate result was 16/48 certified (33% survival) -> under the >=20 bar. Merlin found the cause: my
scripts/clrs_two_arm.sh called discover_gpu WITHOUT --min-years, inheriting the DEFAULT 3, while ALL fleet
boreholes (incl. B926's 530-family benchmark) used --min-years 1. So I compared CLRS at 3x-stricter clustering.
FIX: re-cluster same picker top-30k at --min-years 1 -> 199 families (vs 48). Bar stays frozen at >=20/>=15%
(Merlin rejected lowering it = rule-after-result, and rejected widening candidates below P=0.2 = noise-buying).
Purity shift-null on the 151 NEW (1-2yr) families: coincidence 4.44% vs +120s null 1.11% = 4.0x -> real LFEs.
Densifying the 151 new SEPARATELY (per-cohort survival tracking; junk-guard: new survival <10% => junk).
FLEET RULE: discover_gpu MUST be called with --min-years 1 --min-family-members 3 --cc-threshold 0.80 --fs 40
(bake into the fleet driver; the default min-years 3 is a trap).

## 2026-07-15 — CLRS DISCOVERY PATH VALIDATED ★ FLEET GO
After the min-years fix: densified the 151 NEW families separately (junk-guard).
| cohort | certified | survival |
| original (>=3yr) | 16/48 | 33% |
| new (1-2yr)      | 46/151 | 30% |
| TOTAL            | 62 | ~30% |
New-cohort survival 30% ≈ original 33% => NOT junk (Merlin junk-guard <10% cleanly cleared). 62 >> 20 bar.

IDENTITY CHECK (Merlin step 4) — CLRS certified families ARE the same LFEs B926 sees:
- TIME: family member times coincide with real B926 cc>=0.92 detections at 4-inf x the +120s shift-null.
- SPACE: 92% of CLRS certified families within 10 km of a B926 certified family (100% within 20 km, median 5.5 km);
  exact 0.05-deg bin overlap 27% is just grid quantization.
CONCLUSION: discovery path works at a no-borehole broadband. FLEET IS GO. CLRS = 62-family broadband dv/v station.
