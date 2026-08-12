# Borehole LFE-coda dv/v — new tests & findings (2026-07-10)

Consolidated write-up of the methodological findings made while building the 29-station borehole
3-component pipeline. Context: LFEs are detected by matched-filter "densify" with the cap OFF
(`--top-n 0`), then clustered into families; we produce coda-wave dv/v per family and stack the
reliable ones. Everything below came out of getting that to actually work on B926 and B011.

---

## 1. The certification cascade — how we tell a real LFE family from matched-filter ringing

The central problem: cap-off densify yields hundreds of candidate "families," but most are the
matched filter correlating with noise (ringing), not real LFEs. Three successive tests:

### 1a. Trust battery (250-sample stack-vs-random) — FAILS on cap-off data
The old gate stacked ~250 sampled detections per family and compared the coda amplitude to a
random-time null, calibrated by reversed-template fakes. On the current cap-off family set (B926,
300 families) it **false-negatived nearly everything**: 154 FAIL / 142 UNDET / 4 TRUSTED, with real
families' coda-σ median **−0.09 ≈ fakes' −0.08**. A 250-detection sample against a noise-dominated
stream *is* mostly noise, so real and fake are indistinguishable. (A stale B011 battery on old,
strong hand-picked P>0.9 families looked fine — 47 TRUSTED, real 13.7 ≫ fake 1.8 — proving the
battery works on strong families but collapses on the weak cap-off coverage set we need.) → retired.

### 1b. The detection *stream* is noise-dominated — the foundational finding
Reversed-template densify (time-flipped templates match no real LFE) produced **886 M** detections
for B011 vs **887 M** forward → **rev/fwd = 1.00**. The raw detection *rate* carries ≈ zero LFE
information. You cannot gate families on detection count; only the **stacked coda coherence**
separates real from ringing. This is why the whole certification is grand-stack based.

### 1c. Reverse densify (fwd-vs-reversed grand-stack coda ratio) — the working gate
Per family, `ratio = RMS(fwd grand coda 2–4 s) / RMS(rev grand coda 2–4 s)`, keep ratio > 1.5.
Consistent across stations: **B926 110/300 reliable (37%), B011 119/397 (30%)**; both median
ratio ~1.1 → the median family is ringing. ~1/3 real is a robust property of cap-off borehole
densify. This became the certification.

### 1d. Causality replaces reverse densify — the big efficiency finding
`causality = RMS(fwd coda 2–4 s) / RMS(fwd MIRROR −2..0 s)` on the *same forward stack* — the
pre-arrival mirror window (before the LFE) is just background noise, i.e. the noise floor. Result:
- corr(causality, full fwd/rev ratio) = **1.00** on both stations.
- `caus>1.5` reproduces the reversed-densify reliable set at **96% exact-family overlap**; the few
  differences are threshold jitter (ratio ≈ 1.5) or "mirror-hot" families causality *conservatively*
  rejects.
- Merlin verified it's a **mechanical identity**, per-family: `fwd_mirror / rev_coda` = 1.005,
  log-correlation 0.99, ±5% at 5–95%. Noise autocorrelation is symmetric in lag, so the window
  before the arrival equals the noise floor after it.
- ⚠ Corrections: the reversed floor is NOT a global scalar (CoV 0.4, family-specific — causality
  works because it uses the *per-family* mirror); and "fwd/rev AND causality = two independent nulls
  that agreed" is **retracted** — they are the *same* statistic measured twice (corr 1.00). The
  certification is singly-confirmed, still sound.

**Consequence:** full reverse densify (≈50% of per-station GPU) retired for certification →
causality (>1.5). Keep only a ~2–3% **stratified sampled reverse** per station, for the jobs
causality can't do (see §4 caveat + fake-rate control), plus a per-station identity gate
(`median(mirror/reverse) ∈ [0.9,1.1]`, else fall back to full reverse for that station).

---

## 2. The 3-component (horizontal) investigation

### 2a. Horizontal weakness is PURITY, not alignment
Initial (wrong) diagnosis: horizontals stacked at Z-detection times fail because the S arrival is
mis-timed on horizontals. **Wrong** — S is simultaneous across a sensor's 3 components, and a
constant lag is translation-invariant to stacking. The real cause is that the cap-off stream is
noise-dominated; on the *certified* families + episode days the horizontals do show structure.
⚠ **Per-day lag search is FATAL**: a 1-sample (25 ms) shift injects ~1.25% fake dv/v via the
S-anchor — 5× the real Z signal. Never align-then-stack.

### 2b. H1 dead, H2 coda real
- **H1**: anti-causal (causality median 0.77, 0/110 pass). It carries family-specific energy
  (split-half 0.91) but timed *before* the pick, so the 2–4 s window misses it. Not noise, but not
  usable in the standard window.
- **H2**: the coda *waveform* is real — family-specific and reproducible (base-rate 21% on certified
  vs 2% on non-certified families; split-half within 0.96 / cross 0.00; not narrowband; independent
  of Z, coda cc 0.35).

---

## 3. Reliable families ≈ a third; the reversed floor is nearly free from the forward stack
Quantified certified sets: B926 110/300, B011 119/397. And §1d means we get this from the forward
stack alone — the reversed densify was re-measuring, at 50% of the GPU, a noise floor the forward
mirror window already contains.

---

## 4. THE DEEPEST LESSON — reproducibility ≠ real (stacked-noise-field artifact)

Testing whether H2 *dv/v* (not just its coda waveform) is a usable measurement exposed a subtle,
important trap:
- H2 dv/v **self-reproduces** across disjoint family halves at r = 0.95.
- **But so does H1** (r = 0.83–0.99) — the channel we *certified as pure noise-floor*.
- Therefore split-half family reproducibility does **NOT** certify a dv/v. Family splits do not split
  the contaminating **days**: stacking the same days' noise field across two family halves
  manufactures a highly reproducible *fake* common-mode dv/v.
- The quoted r(Z,H2)dv/v = 0.25 / −0.03 are near-meaningless anyway (two 15-day rolling-median
  series → huge autocorrelation, tiny effective N, no null; |r| ≈ 0.27 arises by chance).

**Correction to an earlier claim:** "H2 PROVEN REAL" established the H2 *coda waveform* is real; it
did **not** establish the H2 *dv/v* is real. Those were conflated. The H2 dv/v is **artifact-
dominated** (worse than noise-dominated).

**Decision:** horizontal coda dv/v is a **waste at boreholes → Z-only fleet-wide.** Keep the free
3-comp stacking + archive H stacks (preserves future direct-S shear-wave splitting / anisotropy and
the H1 anti-causal-window analysis); produce no horizontal dv/v/gate/cert.

**The decisive test (running now): reversed-stack dv/v control.** Run the identical stretch on the
REVERSED (noise-triggered, no real coda) daily stacks; correlate the reversed common-mode with the
forward common-mode, raw and deseasoned (remove trend + annual + semiannual).
- **Applies to Z too, and gates the whole fleet:** deseasoned r(fwd Z, rev Z) ≤ 0.2 → Z payload is
  real coda-borne signal. ≥ 0.4 → STOP: Z is the same stacked-noise artifact and every dv/v
  conclusion (including the seasonal structure) is suspect.
- Preliminary check is encouraging: the certified-Z common-mode is nearly orthogonal to the H1
  noise-mode (deseasoned r +0.10 / +0.02). **[RESULT PENDING — update here when the control finishes.]**

---

## 5. Pipeline / efficiency findings (engineering, but load-bearing)

- **Densify is I/O/CPU-bound, GPU sits at 0% util** during it → run **forward + reversed densify
  concurrently** (both fit the 46 GB L40S) to fill the idle GPU; ~2× per-station throughput, no
  method change.
- **GPU idles during the long CPU stacks tail** (a big station's stacks take hours because it reads
  ~2 B detection rows) → **pipeline GPU densify of station N+1 to overlap the CPU stacks of station
  N**; only the GPU-densify stages are serialized.
- **Gap fills go reliable-only**: after certification is frozen, a re-densify (e.g. the 2017 fill)
  runs only the reliable families, forward only (no reversed, no re-cert) — ~6× cheaper. (The
  reliable-only shortcut applies ONLY to re-densify/gap-fills, never the first pass — you must
  densify all candidates to certify them.)
- **Data QC catch:** B926 & B011 had 2017 truncated to 59 days (day-59 download cut-off) — a mid-
  record gap that faked a discontinuity in the dv/v; the other stations were fine. Always audit
  per-year day counts.
- **★ Fixed candidate P-threshold STARVES sparse/low-calibration stations (B004 caught).** The picker is
  trained on B011, so its P(LFE) drifts LOW at other stations; a fixed P≥0.7 candidate cut over-filters.
  B004: only 2,439 candidates at P≥0.7 → 41 families → ~0 causality-reliable (a false "done"). At P≥0.6
  it has 12.7k (like B926's 14k). B004 is also genuinely tremor-poorer (~12k PNSN tremor within 0.5° vs
  ~75k at B926/B011). **Fix = ADAPTIVE candidate threshold: pick P to yield ~15k candidates, floor 0.5 /
  cap 0.7 — NOT a fixed 0.7.** Safe because causality certifies at the END (a wide early net just costs
  densify; a too-tight early cut is unrecoverable without a full re-run). B004 re-thresholded to P≥0.55
  (22k) and re-queued; starved outputs in data/b004_p70_starved_backup/. Also: FLAG any station finalizing
  with <~20 certified families instead of silently marking it done. (B927/B928 got 300 — not all stations
  are affected, only sparse/low-P ones.)

---

## Status ledger (as of this write-up)
- Certification: causality (>1.5) primary; full reverse retired; ~2–3% sampled reverse kept for
  fake-rate control + per-station identity gate + (former) H2 gate.
- Horizontals: dv/v DROPPED fleet-wide; Z-only. 3-comp stacking kept + archived.
- ⚠ **BLOCKING:** reversed-stack dv/v control on Z must clear (deseasoned r ≤ 0.2) before fleet
  rollout — RUNNING.
- Done: B926, B011 (certified Z dv/v). In flight: B001 finalize, B004/B927/B928 concurrent batch.
- Endgame unchanged: after all 29 stations' certified Z dv/v → multi-window 4-D δβ/β inversion.
