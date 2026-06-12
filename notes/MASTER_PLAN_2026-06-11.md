# MASTER EXECUTION PLAN — family trust, product recovery, and the path back to the 4-D field
*Written 2026-06-10/11 as a self-contained handoff (audience: a fresh Opus session). Read this FIRST, then
the pointed files. The project state changed drastically on 2026-06-10 — do not trust conclusions in older
notes without checking §2 below for whether they were revised.*

---

## 0. CONTEXT — what this project is

Cascadia "tremorferometry": LFE families (catalog-seeded repeating low-frequency earthquakes ON the plate
interface, ~30–45 km deep) are matched-filter-detected at ~35 stations; daily coda stacks → stretch dv/v
(2–4 s lapse, 30-cal-day causal rolling, 1-day step) → multi-window joint inversion → δβ/β on the interface.
GOAL: a crisp, honest, time-resolved (per-month now, **per-day aspiration**) image of megathrust velocity
change — judged against the dv/v literature on DEPTH (the novel axis), sensitivity, resolution, honesty
(memory: success-criteria). Environment & pipeline mechanics: notes/MARGIN_WORKFLOW.md (env paths, CPU cap
32, ONE GPU job at a time, resample_poly only, etc.). User rules: NEVER delete products (flag only); few
stations' waveforms at a time (janitor auto-cleans); Fable agents are ADVISORS only (read-only, never write
code); the logic-supervisor agent (.claude/agents/logic-supervisor.md) should adversarially review every
major conclusion BEFORE you act on it — it caught multiple critical errors.

## 1. THE CENTRAL DISCOVERY OF 2026-06-10 (read carefully — it reframes everything)

**The detector admits noise.** Detection = matched filter, 2-s Z-only templates, FIXED cc≥0.8 ≈ only ~3σ
(N_eff ≈ 2·B·T ≈ 24). Weak/emergent (low-SNR) templates therefore match pervasive ambient/cultural noise
~50–100×/day (analytically predicted ≈ observed). Consequences, all VERIFIED:
- Detection lists at many stations are a MIX of real LFE repeats + noise matches. Margin diurnal screen
  (solar day/night detection ratio; real LFEs are lunar-tidal/episodic, not solar): 10/37 stations >1.5.
- Detection-COUNT statistics (rates, day/night, tremor-corr, continuity) are BROKEN as family verdicts:
  the top-100/day cap censors bursts, and stacking suppresses noise admixture. **Lists can look awful while
  the stacked product is real — B018 proved this.**
- Stacking noise-matched windows = accidental single-station ambient-noise AUTOCORRELATION interferometry:
  a REAL but SHALLOW (near-surface) dv/v signal mislabeled as deep-fault coda. Not "random garbage" —
  mislabeled depth. The B003 EXPANSION (346 low-SNR families added 06-10) is this: verified by the
  SNR-rate cliff (det-rate 95/day below template-SNR 15 vs 2/day above; Spearman −0.46) + solar diurnality
  2.65 + failed dt-vs-lapse dilation test (the "seasonal" there is partly composition artifact).
- **BUT the B018 investigation (notes/2026-06-10 §18) shows the ORIGINAL coverage-selected families carry
  REAL repeating coda**: stack of 300 random detections vs 300 random-time windows → S-window 10.8×, COD A
  (2–4 s, OUTSIDE the 2-s selection support, where selection cannot fake coherence) 3.2× the random null;
  family-dv/v-on-ACF regression slope 0.34 ≈ the kernel's ⅓ receiver-lobe prediction for GENUINE families
  (contamination predicts ~1); coherent non-seasonal drops (2016-09..2017-03 min −0.38%, 98% of 56 families)
  absent from the shallow ACF. **B018's dv/v is REAL.**
- Template SNR is a CREST measure → an SNR≥15 floor is an ANTI-LFE filter (selects spiky non-LFE class;
  98% template-shape-BAD; would delete B001's whole pool). NEVER use raw template SNR as a genuineness cut.
  (cc of a detection proves "correlates with a 2-s wavelet", NOT "is an LFE" — noise matches sit at cc 0.85.)
- DEEP-SIGNAL STATUS: the dv/v is genuinely deep-SENSITIVE (⅓ source-lobe mass; B018 slope confirms) but the
  deep fault is velocity-STABLE — cross-station shared-source coherence ≈ 0 (505 pairs), ETS response ≈ 0.
  Everything that visibly wiggles in any curve is shallow/site/seasonal/noise. The deep result is a BOUND.
  The B018-type coherent drops are non-ACF medium signals of UNDETERMINED depth — candidate, not claim.

## 2. STATUS LEDGER — what stands, what fell, what's pending (with evidence files)

STANDS (verified, safe to build on):
- ETS null / stable interface (cross-station, station-by-station agreement; does NOT depend on contaminated
  families; artifact channel would bias toward false POSITIVE, got null) — fault_tomography/ETS_ANALYSIS.md,
  memory ets-null-and-coda-window (with 06-10 updates).
- Shared-source correlation null (505 pairs, excess −0.006) — data/sharedsource_correlations.csv.
- B933 drop (3 independent lines) — notes/2026-06-10 §11.
- Two-estimator strategy + differential step-injection benchmark (boxcar = background, pairwise = events,
  TV disqualified) — notes/DVV_TEMPORAL_METHODS.md, memory dvv-temporal-estimators.
- Stretch-anchor physics: the fixed point is the PINNED S at t=+1.0 s (envelope-verified), NOT t=0, NOT
  sample 0. dvv_roll30cal.py --origin-anchor now anchors at t_s=1.0. Old products under-read scale
  window-dependently (sample-0: 0.30/0.44/0.56×; t=0-anchor "_calT": 0.50/0.67/0.75× — i.e. the EXISTING
  _calT generation is STILL ~1.5:1 inter-window inconsistent and needs ONE more regeneration with the fixed
  code). Relative/temporal structure unaffected throughout.
- B018 dv/v REAL (the §1 evidence) — figures/B018_stack_vs_random.png, B018_autocorr_vs_family.png.
- Cross-family correlogram: within-station spatial structure is WHITE at monthly cadence (corr flat 0.315
  →0.308 over 1–130 km; station-common mode 0.31) → within-station spatial detail in the monthly fault map
  is noise-fitting; honest spatial resolution = station spacing — figures/crossfamily_correlogram.png.
- Family-shape predictor for the SPIKY contaminant class (AUC 0.92, B933-class) — but it does NOT catch the
  noise-ACF class (emergent templates score GOOD), and its GOOD labels are partly built on "continuity"
  which the noise class manufactures → labels need re-derivation after Tier-1 truth exists.

FELL / REVOKED:
- B003 expansion "gate PASSED" (+87/+17 cells, RMS drop) — the added families are the noise-ACF class;
  the gains were contamination bookkeeping. fault_4d_multiwin_calT35_b003exp.npz = do not use.
- "cc 0.85 proves real LFEs" (notes §16 lesson) — WRONG, corrected in §18.
- "SNR≥15 = real" (my inversion of it) — WRONG, see §1.
- Stage-2 continuity probe as a PURITY gate — adverse selection (noise fires daily; real LFEs are episodic).
  Continuity stays a DISPLAY-quality criterion only.
- The 28-station southern-hot attribution nuance: 35-sta true-scale null test says 40-42°N consistent with
  noise; but scale-fix and station-set changed together (confound; Fable FLAW-4) — the cheap 28-sta-on-_calT
  rerun to separate them was never run. The southern bound stands; the "subarray killed it" story is
  plausible-not-proven.

PENDING / IN FLIGHT (check before redoing):
- Tier-0 margin-wide family audit: scripts/family_trust_tier0.py → data/family_trust_provisional.csv
  (launched 06-10 ~23:00; check logs/family_trust_tier0.log for TIER0_DONE).
- B018 day/night split-stack dv/v test: scripts/b018_daynight_test.sh → logs/b018_daynight.log
  ("DAYNIGHT TEST DONE"); arms must MATCH for medium-carried signal. Calibrates Tier-1b.
- B018 template-free ACF dv/v exists: data/autocorr_dvv_B018.csv (scripts/autocorr_dvv.py — works, but NO
  spectral whitening yet; its own seasonal phase/amplitude not trustworthy until whitened; its NON-seasonal
  part is what the 0.34 slope rests on).
- Continuous-only dv/v map with expanded B003 (fault_tomography/cascadia_dvv_map_continuous.html) — built
  pre-revocation; B003's layer shows contaminated expansion families. Rebuild after trust labels exist.

## 3. THE PROGRAM — phases, gates, exact steps

### PHASE A — Family Trust Test (the current centerpiece; design: notes/FAMILY_TRUST_TEST.md)
A1. Collect Tier-0 output (running): data/family_trust_provisional.csv. Per-station summary → triage
    ranking. REMEMBER: flags ≠ verdicts (B018 lesson). Deliver the table + a margin map of CLEANISH fraction.
A2. Build Tier-1 battery script (scripts/family_trust_tier1.py):
    - T1a stack-vs-random: per family, stack N≈300 random detections (bandpass 2–8 Hz, normalize by 2–4 s
      window std — see the working prototype inline in the 06-10 transcript / figures/B018_stack_vs_random.png);
      K≈10 random-time stacks → null σ; score = coda(2–4 s) amplitude in null-σ units. Efficient impl: group
      window reads by day-file, one pass per station.
    - T1b day/night split: per family, day-only vs night-only stacks → coda cc between arms (and, station-
      level, the already-built b018_daynight_test.sh pattern for dv/v-series agreement).
    - T1c reversed-template null: time-flip ALL the station's templates, run densify_launcher.py on ~30
      sampled days (spread across seasons/years), --start-year/--end-year per chunk or filtered day list;
      push the fake families through T1a/T1b → the station's empirical FAIL distribution; thresholds from
      these, never hand-set. (Also gives per-family false-rate = reversed/forward rate.)
    - T2a shallow-share slope: scripts/autocorr_dvv.py per station (ADD spectral whitening per-hour before
      autocorrelating — Fable flagged the un-whitened seasonal as source-contaminated), then per-family
      deseasoned regression slope of family dv/v on ACF dv/v. Genuine ≈ ⅓, site-carrier ≈ 1, calibrate
      cut with T1c nulls.
    - Scoring → TRUSTED / SITE-CARRIER / CONTAMINATED / UNDETERMINED with continuous scores
      (data/family_trust_<sta>.csv). NOTHING deleted.
A3. VALIDATION GATE (non-negotiable, B018 first — waveforms ON DISK):
    - B018's 56 coverage families: expect mostly TRUSTED (its c107 anchor passed at 3.2×).
    - B018 reversed-template nulls: must score CONTAMINATED.
    - Any real-family CONTAMINATED verdict → inspect individually BEFORE believing the battery.
    - Have the logic-supervisor (Fable) adversarially review the battery results before rollout.
A4. Rollout, one station at a time (priority from A1 triage; re-download → run battery → janitor cleans):
    suggested order: B003 (the open question: expansion expected CONTAMINATED, original 81 genuinely
    uncertain), then the worst Tier-0 stations (likely B004, COR, GNW, B033...), then the rest. ~3–5 h/station
    wall-clock, mostly download. Capture per-station: trust table + the three calibration figures.

### PHASE B — Rebuild certified products (after ≥ the first few stations are labeled)
B1. Regenerate true-scale dv/v with the CORRECTED S-anchor (t_s=1.0; code already fixed in dvv_roll30cal.py):
    a NEW generation suffix (suggest `_calS`) for ALL stations × 3 windows; keep _cal/_calT untouched.
B2. Rebuild per-station dv/v + maps using ONLY TRUSTED (+ optionally UNDETERMINED down-weighted) families;
    SITE-CARRIER families become an explicit per-station near-surface channel (useful: it is an independent
    site-term constraint for the inversion — Fable's FLAW-8 fix — and pairs with the whitened ACF product).
B3. Rerun the multi-window inversion on certified _calS data, with TWO upgrades that are already specified:
    (i) row-weighting by measurement count/SE (unweighted rows let one noisy family = 900-day family);
    (ii) per-(cell,month) resolved masks saved in the npz (per-month resolution accounting — the movie must
    not display Laplacian interpolation as data; Fable FLAW-2). Then rerun null_test.py (SFX/TAG/REAL_NPZ
    env-parameterized) INCLUDING the 28-station-on-_calS variant to settle the southern attribution confound.
B4. Re-derive the headline numbers on certified families ONLY: ETS null bound, network fault index, site
    terms. Expect them to hold (they never rested on the contaminated class) — but verify, don't assume.

### PHASE C — Detector rebuild (the root-cause fix; spec from Fable, 06-10)
The fixed cc=0.8 / 2-s / Z-only detector is mathematically unfixable by thresholding (needed threshold >1).
Rebuild: 3-component (requires re-downloading horizontals — download_station.py already has --channel),
template support [−1,+2] s (S-window only, DISJOINT from the 2–4 s measurement band — no selection-circularity),
per-template per-day MAD-adaptive threshold (~8–9×MAD), network coincidence where station pairs <~50 km.
Pilot ONE on-band borehole (B001 or PGC area) end-to-end; compare against Tier-1-certified families; only
then consider fleet re-detection. This phase is big — do not start it before Phases A/B are banked.

### PHASE D — Back to the 4-D field and the per-day goal
- Spatial honesty: present the field at station-spacing resolution (correlogram result); fault map = bound.
- Per-day: pool certified families per station (genuine-events/day budget — compute it per station from
  certified detections; ~100/day at good stations) via the pairwise estimator + Kalman layer; cross-station
  coherence at daily cadence is the deep-transient alarm. 4-D injection-recovery harness (space-time
  checkerboard through real monthly sampling) replaces the static checkerboard as the resolution metric.
- Family EXPANSION (the original goal that started 06-10): only after the trust battery gates it — audition
  pool → template-shape screen → Tier-1 battery on probe → full densify of certified additions only.

## 4. PITFALLS / LESSONS (cost a full day — do not repeat)
1. VERIFY before interpreting; the lead flip-flopped 4× on B003/B018 by asserting mechanisms ahead of data.
   Pattern that worked: user skepticism → decisive cheap test → verdict. Use the logic-supervisor early.
2. Correlation is usually the wrong statistic (use regression SLOPE for mixes; deep-stable makes even
   genuine families correlate with site signals).
3. Anything inside the selection window is coherent by construction; test OUTSIDE it.
4. Detection-count statistics cannot condemn stack products (cap censoring + stack noise-suppression).
5. "Continuity" selects FOR the noise class. "High template SNR" selects AGAINST LFEs. Use template SHAPE
   (spiky class) + Tier-1 waveform tests.
6. Plotting: the 30-d/1-d rolling IS the smoothing — never smooth on top (user explicitly banned 60-d
   medians etc.); plot data as-is; y-limits honest (no clipping); the established style = faint per-family
   + bold cross-family median (see smoke_dvv_*.png).
7. pandas traps that bit repeatedly: .corr/.kurt/.cov are DataFrame METHODS — name columns differently;
   Date.now-type nondeterminism banned in workflows; check `_calT`-vs-`_cal` scale (×2.26) when comparing.
8. Memory/notes hygiene: every conclusion that later fell was already written to notes/memory — when a
   verdict changes, UPDATE the old entry (notes/2026-06-10 §16 has the correction pattern).
9. Background jobs: launch tracked; only ONE GPU job; janitor auto-deletes waveforms once stacks exist
   (pause it if waveforms must persist, as for Tier-1 testing).
10. Fable agents: advisory ONLY (read-only). They run their own checks and have been right 3-of-3 times
    when the lead was wrong — but they also over-condemned once (the "all contaminated" cascade needed the
    B018 stack test to scope correctly). Adversarial review ≠ oracle; decisive tests rule.

## 5. KEY FILE INVENTORY (for this program)
- Design doc: notes/FAMILY_TRUST_TEST.md. This plan: notes/MASTER_PLAN_2026-06-11.md.
- Day log with all 06-10 findings: notes/2026-06-10_Notes.md (§14-18 are the critical ones).
- Tier-0 output: data/family_trust_provisional.csv (+ logs/family_trust_tier0.log).
- B018 evidence: figures/B018_stack_vs_random.png, figures/B018_autocorr_vs_family.png,
  data/autocorr_dvv_B018.csv, data/b018_family_audit.csv, scripts/b018_daynight_test.sh (+ its outputs
  data/daily_dvv_B018_2to4_{dayonly,nightonly}.csv, data/long_window_daily_B018_{dayonly,nightonly}.npz).
- Tools: scripts/autocorr_dvv.py (template-free ACF dv/v; ADD WHITENING), scripts/family_audition.py
  (Stage-1 screen + probe machinery), scripts/family_trust_tier0.py, scripts/dvv_roll30cal.py (S-anchor
  fixed, t_s=1.0), scripts/bench_step_recovery.py + scripts/dvv_pairwise.py (estimator suite),
  fault_tomography/inversion/{invert_multiwin,null_test}.py (env-parameterized SFX/OUT/TAG/REAL_NPZ).
- Agent definitions: .claude/agents/logic-supervisor.md (Fable advisor), .claude/agents/station-dvv-pipeline.md.
- Memory: ~/.claude/.../memory/MEMORY.md index — read success-criteria, dvv-temporal-estimators,
  family-cwi-predictor, ets-null-and-coda-window, project-mandate-2026-06-09, margin-station-status.
  ⚠ memory family-cwi-predictor + margin-station-status contain pre-06-10-evening claims that §2 above
  revises (cc-discriminator lesson, B003 expansion success) — update them as Phase A lands.

## 6. IMMEDIATE NEXT ACTIONS (in order)
1. Check the two in-flight runs (Tier-0 sweep; B018 day/night) — collect, summarize, append to notes.
2. Build Tier-1 (A2), validate on B018 (A3), Fable review.
3. B003 re-download + full battery (A4 start) — settles the open "original 81" question.
4. Phase B regeneration (S-anchor `_calS`) once ≥5-10 stations are labeled.
5. Keep the user in the loop at each gate; they drive priorities. Concise reporting; no over-smoothing;
   honest axes; flag-don't-delete.

---
## 7. UPDATE — 2026-06-11 PM (supersedes the matching items in §2/§6; full log notes/2026-06-11_Notes.md §13)

**PHASE A COMPLETE — trust-battery rollout done, 10 stations.** Two gates standardized: GOLD = coda_sigma >
station's reversed-template fake-MAX (strict tail); TRUSTED = > fake-95th. Totals 341 GOLD / 433 TRUSTED /
114 FAIL over 781 scored families. Rich tier B022(67 GOLD)/B040(62)/B001(60)/B018(40)/B201(23); solid
B928(33)/B014(19); marginal B003(16)/B017(15)/B935(6). Elevated-floor stations B201/B022/B014 (fakes reach
14-16sigma) have a coherent persistent site signal; per-station reversed calibration auto-compensates.
Files: data/family_trust_tier1_<STA>.csv + data/family_trust_master.csv. Memory: [[trust-battery-rollout-2026-06-11]].

**CORRECTION 1 — detection-level cleaning is CIRCULAR (discarded).** Selecting detections by coda-match to the
all-time-mean coda, then measuring stretch-deviation from that same mean, suppresses the signal it measures;
scatter/seasonal fell for ALL families incl. the FAIL control. cc_max->0.98 is a selection artifact. Do NOT use.

**CORRECTION 2 — T2a is a SITE filter, NOT a depth certificate (Fable, big reframe).** T2a (regress family dv/v on
station Z-ACF dv/v) measures only whether a family's dv/v TIME-SERIES tracks the shallow signal = contamination
filter. It does NOT certify deep KERNEL (geometry) or deep CHANGE (temporal). Do not read slope~1/3 as deep-proof.
T2a-FIXED (block-bootstrap CIs, scripts/exp_t2a_fixed.py) on B018: 21/39 GOLD SHALLOW-COUPLED (candidate
site-carriers), 15 DECOUPLED (deep-or-noise), 3 ANOMALOUS. **T1 sigma is ORTHOGONAL to shallow-coupling
(Spearman -0.08): cannot sigma-weight out site-carriers -> the inversion input MUST be T2a-filtered, not just
trust-ranked.** Memory: [[detection-cleaning-circular-t2a-gap]].

**DEEP STATUS reaffirmed:** deep velocity CHANGE = upper bound ~0.02% (positive confirmation not possible with
current data; inter-station coherence is confounded by regional shallow climate AND already null = 505-pair
shared-source). Deep SENSITIVITY (geometry) IS positively confirmable via dt-vs-lapse DILATION + multi-window
kernel -> requires the `_calS` per-window true-scale regen (still pending; `_calT` ~1.5:1 inter-window inconsistent).

**REVISED NEXT ACTIONS (await user steer at the fork):**
1. (recommended next) `_calS` per-window true-scale regen -> lapse-dilation deep-SENSITIVITY test (the honest
   positive-depth path; Phase B1). 
2. Apply T2a site filter across stations before any inversion (needs each station's Z-ACF dv/v; only B018+B935
   computed -> others need re-download). Inversion uses DECOUPLED, not all-GOLD.
3. Disk: 13GB leftovers (ERW/SNB/LRIV/LZB) + B935 41GB (poor pilot, dead ambient) pending user OK to clear.
