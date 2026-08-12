# RESUME 2026-07-10 — Borehole 3-component dv/v pipeline (READ FIRST after reset)

## GOAL
3-component (Z/EH1/EH2) coda dv/v for **each family** of **each Cascadia PB borehole (29 stations)**,
from **real LFEs**. Inversion (4-D fault field) = FINAL step after ALL 29 stations' dv/v.
Station list = the 29 trust-battery boreholes: `results/family_verdicts/ALL_families_labeled.csv`.

## ★ CRITICAL FIXES THIS SESSION (do not regress)
1. **Densify truncation bug (FATAL, fixed):** `densify_gnw_gpu.py` had `N_FIX=3_456_016` hardcoded
   for 40 Hz → truncated 100 Hz days to first 9.6 h. Now `N_FIX=int(round(fs))*86400+16` (fs-aware).
2. **Run everything at 40 Hz** (Option B, Merlin+user). Densify/stacks/dv/v at **fs=40**; discovery+picker
   stay 100 Hz (fingerprint trained at 100 Hz — do NOT move it). 40 Hz makes N_FIX correct, matches the
   28-station fleet the inversion assembles with, 2.5× faster. Resample template npz 100→40 once
   (`resample_poly 2/5`, m 200→80): `data/<sta>_disc_p70_2010_2026_m3_40hz.npz`.
3. **Transient-Tk GPU fix:** densify recomputes template FFT per TB=128 batch (not resident) → any family
   count fits (was OOM at >~600 templates). GPU ~15-19 GB.
4. **Memory cap = 187 GB cgroup** (`free` lies/shows 1.5 TB). Watch `anon` in /sys/fs/cgroup/memory.stat,
   NOT total (cache). OOM watchdog: `scripts/oom_watch.sh` (alerts anon>155). 3-ch stacker streams per-year.
5. **CUDA-12 libs:** prepend `/opt/conda/lib/python3.13/site-packages/nvidia/*/lib` to LD_LIBRARY_PATH +
   `CUDA_PATH=/opt/conda/targets/x86_64-linux` for every GPU step. (memory: cupy-cuda12-libpath)
6. **Kill jobs by** `ps -eo pid,comm,cmd | awk '$2~/python/ && /densify_gnw/'` — pgrep/pkill -f MATCH YOUR
   OWN shell (the pattern is in the cmdline) and leave orphans. Orphan spawn-workers pile up; clean with
   `ps -eo pid,ppid,cmd | awk '$2==1 && /multiprocessing.spawn/{print $1}' | xargs kill -9`.

## PER-STATION PIPELINE (borehole-first; land stations deferred)
1. **Download** 3-comp EH1/EH2/EHZ 2007-2026 (`download_borehole_3comp.py`) — NET. Keep **6 boreholes**
   buffered; delete a station's traces only AFTER its trace-dependent steps (densify+reversed+stacks+trust
   nulls) done, then download next (rolling buffer). dv/v+inversion use saved stacks, no traces.
2. **Candidate detection** (`discover_nllb_pnsn_driven.py --candidates-only`, fs=100, master catalog
   `pnsn_tremor_cascadia_full.csv`, bbox=±0.9lat/±1.35lon) — CPU. Many stations already have
   `<sta>_pnsn_candidates_100km.parquet` from prior sessions.
3. **Picker score** (`score_candidates.py --thr 0.7`, needs `tremor_picker_<sta>.joblib` = COPY of
   `tremor_picker_b011.joblib`; fingerprint transfers, only P-calibration drifts) — CPU. keep P≥0.7.
4. **Cluster** (`discover_gpu.py --fs 100 --cc-threshold 0.80 --min-family-members 3 --min-years 1`) — GPU (~15s).
5. **App B** (`score_family_stacks_picker.py`) → keep **pred==LFE label** (NOT P>0.9 — cross-station P runs
   low; label is calibration-robust) — CPU.
6. **Coverage select ~300/station** (`select_families_coverage.py <disc_prefix> <picker_csv> 300`): dedup
   near-dup templates (cc≥0.9) → keep every 0.05° bin (coverage) → super-high P>0.95 GUARANTEED → fill by
   P(LFE) to ~300. Writes `<disc>_sel300.summary.csv`. (B011→397, B926→300.)
7. **Densify** (`densify_gnw_gpu.py --fs 40 --top-n 0 (CAP OFF) --despike-mad 8 --min-snr 0`, templates=
   40hz npz, summary=sel300) → `mf_<sta>p90f40_*.csv` — GPU (~1.5-2h/station).
8. **Reversed densify** (SAME but templates time-flipped: `..._40hz_rev.npz`) → `mf_<sta>p90f40rev_*.csv`
   = noise-match floor — GPU. **FINDING: fwd≈rev in COUNT (noise-dominated stream!) → cannot gate on count,
   must gate on STACK coherence (trust battery).** User chose FULL reversed densify (not sampled).
9. **3-ch stacks** (`build_long_window_3comp.py --fs 40 --mf-csv-glob 'data/mf_<sta>p90f40_*.csv'`) →
   `long_window_daily_<sta>p90f40_{Z,H1,H2}.npz` — CPU.
10. **Trust battery** (`family_trust_tier1.py --mf 'mf_<sta>p90f40_*.csv' --rev-mf 'mf_..rev_*.csv'` — I
    edited --mf to accept globs) → `family_trust_tier1_<sta>.csv` (GOLD/TRUSTED/UNDET/FAIL) — CPU. Needs traces.
11. **dv/v** (`dvv_roll30cal.py --npz ..._{Z,H1,H2}.npz --window 2 4 --origin-anchor` [MUST use
    origin-anchor, t_s=1.0] → `daily_dvv_<sta>p90f40_<C>_2to4.csv`) — CPU. Run per component.
    ⚠ PENDING: gate dv/v by trust verdicts + reversed-noise before it's "certified".

RESOURCES: GPU = serial bottleneck (cluster+densify+reversed). Everything else CPU (parallelize: discover
next stations on CPU while GPU densifies) or NET. GPU=L40S 46GB. CPU cap 32 (workers ≤ ~30 total).

## CURRENT STATE (2026-07-10 ~07:00 UTC)
- **B926 DONE through dv/v** (300 families): fwd+rev densify done, 3-ch stacks done, 3-comp dv/v done
  (`daily_dvv_B926p90f40_{Z,H1,H2}_2to4.csv`, 1.35M rows each). **UNGATED** — Z coherent ±0.5%, H1/H2 noisy
  (rail-hit on noise days). Trust battery RUNNING (`family_trust_tier1_B926.csv` pending). Figures:
  b926_3comp_dvv, b926_family_3comp_dvv, b926_300_traces, b926_family_dailystack (family c13 = clean
  16-yr coherent LFE), b926_sel300_map.
- **B011:** 397 families selected (`b011_disc_p70_2010_2026_m3_sel300.summary.csv`), forward densify
  RUNNING (chained via `scripts/chain_b011_densify.sh`, ~1/17). Then needs reversed+stacks+trust+dv/v.
- **B001, B004:** CPU discovery prep running/queued (`scripts/discover_prep_queue.sh` — B004/B001/B003/
  B927/B928; candidates+picker-score; clustering pending GPU).
- **Waveforms on disk (6):** B001 B004 B011 B926 B927 B928 (+PGC). Downloads STOPPED (buffer full);
  resume one-at-a-time as stations clear (`download_borehole_3comp.py`).
- **27 stations still need full discovery→dv/v.** GPU-serial ≈ multi-day.

## IMMEDIATE NEXT
1. B926 trust battery done → **gate B926 dv/v** (drop FAIL families + noise days where fwd stack coda σ
   ≯ reversed) → re-plot certified 3-comp dv/v (should tighten H1/H2).
2. B011 densify done → reversed densify → stacks → trust → dv/v.
3. Cluster (GPU) the CPU-prepped stations (B004/B001...) → App B → select → densify. Rolling buffer.
4. After ALL 29 dv/v: multi-window inversion (`fault_tomography/inversion/invert_multiwin.py`).

## ★★ HORIZONTAL dv/v DROPPED — Z-ONLY FLEET (Merlin, 2026-07-10) + ⚠ Z CONTROL PENDING
- **Horizontal (H1/H2) coda dv/v is a WASTE → Z-only fleet-wide.** Keep 3-comp STACKING (free, one pass)
  and ARCHIVE H stacks, but produce NO horizontal dv/v / gate / cert / finalize branch.
- ⚠ **KEY CORRECTION (I conflated two things):** the H2 coda *WAVEFORM* is real (base-rate 21%v2%,
  split-half 0.96 — that stands), BUT the H2 *dv/v* is a **stacked-noise-field ARTIFACT**, not signal.
  Proof: the dv/v-family split-half test self-reproduces at r=0.95 for H2 — **but ALSO 0.83-0.99 for H1**,
  the channel we certified as pure noise. So "reproducible ≠ real": stacking the SAME days' noise field
  across disjoint family halves manufactures a reproducible fake common-mode. Family-split does NOT split
  the contaminating DAYS. The r(Z,H2)dv/v=0.25/-0.03 I quoted are near-meaningless (two 15-day rolling
  medians, huge autocorrelation, no null; |r|~0.27 arises by chance).
- ⚠⚠ **DECISIVE control now RUNNING (`scripts/rev_dvv_control.py`):** run the SAME stretch on the REVERSED
  (noise-triggered, no real coda) daily stacks, correlate rev common-mode vs fwd. **Z acceptance:
  deseasoned r(fwdZ,revZ) ≤ 0.2 → Z payload real; ≥ 0.4 → STOP FLEET, Z is the same artifact & all dv/v
  conclusions (incl seasonal structure) suspect.** Merlin's prelim check is ENCOURAGING (certified-Z
  common-mode ~orthogonal to H1 noise-mode: deseasoned r +0.10/+0.02), but the rev-Z control is the proper
  cert. This was the open contamination item in this file — MUST clear before fleet rollout.
- Deferred horizontal science (from archived stacks, later): direct-S shear-wave splitting per era;
  H1's anti-causal family-specific energy via a per-component window.

## ★★ CAUSALITY REPLACES FULL REVERSE DENSIFY (Merlin-verified, 2026-07-10) — big fleet speedup
- **Reliable families = causality (fwd coda 2-4s RMS / fwd MIRROR -2..0s RMS) > 1.5** on the forward grand
  stack. This reproduces the full fwd-vs-reversed reliable set at **99% agreement** (B926/B011), corr 1.00.
- **WHY (mechanical identity, Merlin verified per-family):** the pre-arrival mirror window IS the noise
  floor that reverse densify measures — fm/rc median 1.005-1.009, log-corr 0.99, ±5% at 5-95%. Noise
  autocorrelation is symmetric in lag → pre-window ≈ post-window floor. So reverse densify is REDUNDANT
  with the mirror for CERTIFICATION.
- ⚠ **Two corrections (my earlier claims were WRONG):** (1) the reversed floor is NOT a "global scalar" —
  CoV 0.46 = ×3 spread, it's family-specific; causality works because it uses the PER-FAMILY mirror. (2)
  "fwd-vs-rev AND causality = two independent nulls that agreed" is RETRACTED — they are the SAME statistic
  measured twice (corr 1.00). Certification is SINGLY-confirmed (still sound, just not doubly).
- **★ FLEET PROTOCOL (Merlin):** (a) certify Z on **causality>1.5** (primary, no full reverse). (b) Keep a
  **K≈150-day STRATIFIED sampled reverse per station** (~10/yr + ~50 episode days, ~2.7% of full reverse)
  — NOT for certification but for 3 jobs causality can't do: the **fake-RATE λ(t) control** (amplitude null
  ≠ rate null; still a MUST per the 07-10 audit), **per-station identity verification**, and the **H2 gate**
  (episode-day reversed H2 stacks). (c) **Per-station identity gate:** from the sampled reverse check
  median(fm/rc)∈[0.9,1.1] AND Spearman(fm,rc)≥0.95 → certify on causality; else run FULL reverse for THAT
  station (mirror≠floor there, e.g. P-coda in -2..0s at different slab geometry, or Z anti-causal energy).
  (d) FIRST rollout station: K=365 for a tighter out-of-sample identity check, then drop to K=150.
  (e) H2 gate KEEPS fwd/rev on horizontals (H1 proved horizontal mirror carries anti-causal energy).
  (f) Log mirror-hot families (caus<1.5 but sampled fwd/rev>1.5) as REJECT-flag; >10% mirror-hot → stop.
- **Current batch (B004/B927/B928) keeps full reverse** (already running concurrent = ~free wall-clock, and
  gives identity validation on 5 stations); switch to this cheap protocol AFTER the batch. Experiment:
  `scripts/cheap_cert_experiment.py`.

## ★ TRUST BATTERY vs REVERSE DENSIFY (B926+B011 evidence, 2026-07-10)
- **Detection stream is noise-dominated (crisp proof):** B011 forward 887M detections vs REVERSED 886M →
  **rev/fwd = 1.00**. Time-flipped template (matches no real LFE) triggers as often as the real one →
  raw detection RATE carries ≈0 LFE info; NEVER gate on count. Only grand-stack coda coherence separates.
- **Reverse-densify cert is reliable + consistent:** B926 110/300 (37%) real, B011 119/397 (30%) real; both
  median fwd/rev ratio ~1.1 → the MEDIAN family is ringing (symmetric coda). ~1/3 real is a robust property
  of cap-off borehole densify. `data/<sta>_fwd_vs_rev_coda.csv`, ratio>1.5 = certified.
- **Trust battery (250-sample) is underpowered on cap-off sets:** B926 (current 300 fams) → 154 FAIL/142
  UNDET/4 TRUSTED, real coda_σ −0.09 ≈ fake −0.08 (false-negatives all). A STALE B011 battery (2026-06-11,
  71 OLD strong P>0.9 fams, 0 overlap w/ current) looked fine (47 TRUSTED, real 13.67≫fake 1.77) → the
  battery works on strong hand-picked fams but COLLAPSES on the weak cap-off coverage set. Verdict depends
  on family set, not station. RETIRED as the gate; certify on fwd-vs-rev grand stacks (null = full record,
  ~1e8 detections, vs the battery's 250) cross-checked by causality. Battery = fast screen only.

## ★ B926 DEEP-DIVE FINDINGS (station 1 — apply to all)
- **Trust battery (family_trust_tier1.py) is UNDER-POWERED for cap-off data** — its 250-detection
  sample ≈ pure noise (stream is noise-dominated, fwd densify count ≈ reversed count) → false-FAILs
  everything (B926: 154 FAIL/142 UNDET/4 TRUSTED/1 GOLD; F coda_sigma ≈ R). DO NOT trust its verdict.
- **CERTIFY on GRAND stacks** (Merlin): the null at median 1.15 measured two EQUIVALENT ways (NOT
  independent — see the CAUSALITY section; corr 1.00, they are the same coda/noise-floor statistic) —
  (a) causality: grand-stack coda(2-4s) RMS / mirror(-2..0s) RMS [PRIMARY, no reverse needed];
  (b) fwd-vs-reversed: fwd grand coda / rev grand coda (needs reversed daily stacks). Keep
  families with ratio > 1.5. B926: **110/300 real-coda** (rest = matched-filter RINGING, symmetric coda).
  Saved data/b926_fwd_vs_rev_coda.csv.
- **DAY-GATE (n_fwd≥2·n_rev) is TOO AGGRESSIVE** — fwd≈rev means only ~2% of days survive → median gets
  NOISIER. DON'T over-gate. **The MEDIAN over certified families is the robust estimator** (ringing
  averages to ~0). B926 certified Z dv/v median std 0.23% = usable. (per-day counts are in npz n_det.)
- **★★ 3-COMPONENT VERDICT (Merlin-guided, decisively tested on B926 — supersedes my earlier wrong
  "alignment" story).** My alignment diagnosis was PHYSICALLY WRONG: S is simultaneous across a sensor's
  3 components; a constant lag is translation-invariant to stacking. Horizontals' weakness is a PURITY
  problem, not alignment. ⚠ **Per-day lag search = FATAL** (a 1-sample/25ms shift injects ~1.25% fake dv/v
  via the S-anchor = 5× the Z std). NEVER align-then-stack.
  Purity-stratified H test (zero GPU, existing fwd+rev daily-stack npz, restrict to 110 certified families
  + PNSN episode days; `scripts/purity_h_test_b926.py`, `purity_h_decisive_b926.py`):
    · **Z: 97/110 pass** (caus median 2.11, fwd/rev 2.14) — positive control OK.
    · **H1: DEAD** — caus median 0.77 (ANTI-causal), 0 pass. Not recoverable at boreholes.
    · **H2: 23/110 pass** (caus>1.5 & fwd/rev>1.5). PROVEN REAL, not artifact, by two decisive tests:
      base-rate **21% certified vs 2% non-certified** (gate is informative, not a correlated-noise tail);
      split-half within-vs-cross **gap 0.96 ≈ Z's 0.98** (family-specific reproducible coda; a shared
      artifact waveform would give within≈cross). Not narrowband (Guard C). r(Z,H2)dv/v=0.25 → H2 is
      NOT Z-bleed (independent) but noisy.
  **Gated 3-comp dv/v** (`gate_3comp_dvv_b926.py` → figures/b926_certified_3comp_gated_dvv.png):
    Z 110 fam **std 0.22%** (science-grade); H2 23 fam **std 0.59%** (2.7× noisier, weakly tracks Z).
  **★ POLICY (Merlin):** fleet default **Z-only is the science payload**; the stacker ALREADY writes
  H1/H2 stacks at Z-times for FREE, so the H-gate is a free CPU by-product per station — each station
  EARNS an H2 second product only by passing the gate on its own data (no 3× GPU; reject that framing).
  H2 = a NOISY consistency check on Z, NOT a new science axis; the 4-D inversion does NOT need it. H1 is
  not recoverable at boreholes (anti-causal). Deliverable per station = **Z certified dv/v + H2 dv/v on
  the families that pass the free gate; H1 reported as failing the physical-reality gate.**
  · Asymmetry checks (`purity_h_asym_b926.py`): certified-family azimuth span 146-330° (range 185°) →
    radiation-pattern does NOT explain H1-dead; H1/H2 raw noise floor 1.31× → H1 is NOT just a noisier
    axis; H2-vs-Z coda |cc| median **0.35** → H2 is INDEPENDENT of Z (not bleed), a genuine independent
    horizontal measurement. ⓘ H1 nuance: split-half gave H1 within 0.91/cross 0.00 = it DOES carry
    family-specific energy, but ANTI-CAUSALLY timed (peaks in mirror -2..0 s, before the Z pick) → the
    standard 2-4 s window misses it. H1 is real signal in the wrong window vs the Z detection, NOT noise
    — a future shifted-window (per-component pick) candidate, not usable in the standard pipeline now.
  · Reusable scripts: `certify_fwd_vs_rev.py <TAG>` (Z cert), `finalize_3comp_dvv.py <STA> <TAG>`
    (cert + H-gate + gated Z+H2 dv/v + figure + <sta>_3comp_summary.json). Per-station driver +
    5-station batch (B001/B004/B003/B927/B928) delegated to a background rollout agent; status in
    logs/rollout_status.log. B011 auto-finalizes via chain_b011_certify.sh.

## OPEN CERTIFICATION ITEMS (Merlin audit, notes/2026-07-10_Notes.md)
- Gate dv/v by trust battery + reversed-noise BEFORE trusting (cap-off stream is winter-modulated noise;
  seasonal noise can fake seasonal dv/v). NOT done yet — dv/v so far is ungated/preliminary.
- Replace weak Lin-coincidence test with **B011↔PGC cross-instrument coincidence (±2 s)** (co-located,
  still not run) — the real external-source check.
- 3-comp agreement ≠ deep; depth needs T2a decoupling + lapse-dilation (later).
