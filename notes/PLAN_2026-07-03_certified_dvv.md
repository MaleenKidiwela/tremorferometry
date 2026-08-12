# PLAN 2026-07-03 — Certified-LFE dv/v: certification protocol + depth-resolved inversion
*Merlin-audited (Fable 5, adversarial review of the naive plan "GOLD∩picker∩T2a → _calS regen →
multi-window inversion"). Supersedes the Phase-B section of notes/MASTER_PLAN_2026-06-11.md as the
execution plan; the MASTER_PLAN status ledger (§2, §7) still governs what stands/fell.*

**User goal:** dv/v separated by depth (shallow / [intermediate] / deep megathrust) via the
multi-window inversion — but FIRST a defensible way to confirm the input families are REAL LFEs.

---

## 0. FLAWS FOUND IN THE NAIVE PLAN (severity-ordered; each shaped a step below)

1. **FATAL — B011-p50 "274 real LFE families" is CIRCULAR + domain-shifted.** Candidates were
   pre-filtered P(LFE)≥0.5 by `tremor_picker_b011.joblib`, clustered, then family stacks scored by
   the *same model* — and the model was trained on single windows but applied to stacks (the exact
   R2 extrapolation failure lfe_features/RESULTS.md documents). The 274 = a candidate list, NOT a
   certification.
2. **FATAL as fleet recipe — GOLD∩picker∩T2a is undefined at ~26/29 stations** (T2a at 3 stations;
   picker validated at B011+PGC only, needs per-station retraining per R5). Intersection silently
   degrades to GOLD-only, which is proven insufficient for the deep channel (T1 σ ⊥ shallow-coupling).
3. **FATAL for 3-layer framing — invert_multiwin.py has NO intermediate layer** (model = fault cells
   + per-(station,window) scalar site terms). Overlapping 1-3/2-4/3-5 s windows share kernel lobes →
   dt-vs-lapse supports a 2-parameter split (intercept=shallow, slope=distributed dilation).
   An intermediate layer without proof = regularization artifact.
4. **MAJOR — _calS validated at ONE station×window** (B928 2to4; slope 1.491 vs _calT ≈ predicted
   1/0.67 ✓ encouraging). No 1to3/3to5 _calS exists → "inter-window inconsistency fixed" is untested.
5. **MAJOR — B3 inversion upgrades never implemented**: monthly aggregation unweighted (no n/SE),
   LSQR rows uniform, no per-(cell,month) resolution masks in npz.
6. **MAJOR — checkerboard is an inverse crime** (same G forward+inverse, deterministic sin() noise).
   Use the null_test.py residual-noise pattern for any resolution claim.
7. **MAJOR conceptual — T1 and the picker SHARE their dangerous failure mode**: a persistent coherent
   low-frequency cultural/site repeater passes T1 (elevated-floor stations show fakes at 14-16σ) and
   is the picker's weak class (in-tremor-ambient AUC 0.81). Their agreement ≠ certification. Only
   EXTERNAL-source evidence certifies against it → Gate S below (the missing layer).
8. **MAJOR preventive — do NOT substitute T2a with family-derived common-mode** at T2a-less stations:
   deep signal is common-mode too → regressing it out biases the deep bound toward null (flattering).
9. MINOR: cc_max>0.7 row-selection stretch-bias (quantify once); TNOISE label impurity; Lin train/test
   event leakage 2010-13 (quote holdout numbers); regenerate the Jun-10 B928 calS file (predates the
   final anchor commit) rather than trust it.

---

## PART 1 — CERTIFICATION PROTOCOL ("am I relying on real LFE families?")

Three gates with DISTINCT roles (picker demoted to seed-filter + raw-feature character table;
continuous scores used for weighting within the certified set, never as a gate substitute):
- **Gate R** (real coherent repeater) = T1 trust battery, reversed-fake-calibrated — EXISTS, 29 stations.
- **Gate S** (EXTERNAL source, not station-local) = NEW, the missing layer: cross-instrument
  coincidence + cultural-periodicity nulls + out-of-sample episodicity.
- **Gate D** (depth-use eligibility) = T2a-DECOUPLED, deep channel only; SHALLOW-COUPLED families are
  RETAINED as shallow-channel data.

### Steps (decisive + cheap first)
- **P1.1 (free) Freeze the certified-input definition** in a spec note: certified-deep =
  R-pass (GOLD, or TRUSTED+S-pass) ∧ S-pass ∧ D-DECOUPLED (D evaluated only at T2a-capable stations;
  others feed shallow field + site terms only). Spec names exact file/column per gate.
- **P1.2 (cheap CPU, DECISIVE) B011↔PGC cross-instrument coincidence.** Co-located (~150 m),
  different networks/instruments; detection lists ON DISK (data/mf_b011p50_*.csv, data/mf_pgc_*.csv,
  old mf_b011_*). Per family: fraction of detections with a PGC detection within ±2 s; null =
  within-day circularly-shifted times; negative control = reversed-fake families. Earth-source must
  appear on both; station-local noise must not.
  GATE: real ≫ shuffled ≫ fakes. **If old B011 GOLD families sit at chance → STOP EVERYTHING**
  (T1 passes station-local artifacts wholesale; certification concept fails; report to user).
- **P1.3 (cheap-med, waveforms on disk) Picker-independent T1 on the 314 p50 families** via
  scripts/family_trust_tier1.py + fresh reversed-template nulls for the p50 set. This REPLACES the
  circular "274". GATE: picker-LFE families pass ≥ old-B011 GOLD rate (38/71). FAIL (<~40% beat
  fake-max) → p50 pilot dead as inversion input (and that's an important finding).
- **P1.4 (cheap) Raw-feature character table** on densified 3-comp stacks (sp_slope, centroid,
  bf_2_4, hf_ratio, rolloff85 + stack polarization — R2's open caveat) vs three REFERENCE
  distributions: CC-ALIGNED Lin stacks (fixes R2 alignment confound), reversed-fake stacks, the 35
  EQ-multiplet stacks. Thresholds from references, not a model. Replaces p_lfe as character verdict.
- **P1.5 (cheap, list-based) Cultural nulls the battery lacks:** (i) weekday/weekend rate ratio
  (7-day periodicity = anthropogenic, day/night misses it); (ii) OUT-OF-SAMPLE episodicity — p50
  seeded from 2010-13 tremor windows but densify covers 2007-17 → 2007-09/2014-17 detections
  concentrating in PNSN ETS episodes is a genuine held-out test. Mind cap-censoring (flags → inspect,
  split-stacks where censored).
- **P1.6 T2a wherever waveforms exist** (B011, PGC, B926, B927 now), whitened autocorr_dvv +
  exp_t2a_fixed. Verify the existing B018/B927/B941 T2a used the WHITENED ACF (check dates vs
  whitening commit). Institutional rule: every future station download runs autocorr+T2a before
  janitor. NEVER family-common-mode as substitute (Flaw 8).
- **P1.7 (free) Cross-layer audit table** T1×character×coincidence at B011: DISAGREEMENT cells get
  individual inspection; agreement rate is NOT certification (Flaw 7) — say so in the note.

**Sufficiency after Gate S added:** R kills incoherent noise; S kills station-local coherent signal
(the class R+picker both miss); character kills EQ/blast multiplets; D keeps site-carriers out of the
deep channel. Remaining honest caveat: absolute source depth/location still rests on catalog seeding +
cross-station coincidence — keep that caveat in every deliverable.

---

## PART 2 — DEPTH-RESOLVED dv/v PROGRAM (ordered stages, each with a falsifiable gate)

- **S0 (cheap CPU, FIRST — can invalidate the deliverable) Identifiability test.** Real
  (family,station,window) geometry + kernels → forward operators for 2-layer (site+fault) vs 3-layer
  (site+mid-path+fault). Inject synthetic intermediate-only anomaly; measure leakage; inspect SVD of
  stacked multi-window G. Evaluate BOTH the current 1-3/2-4/3-5 ladder AND adding a 4-7 s window
  (stacks run −3..+10 s → free; the only real lapse-lever extension).
  GATE: intermediate recoverable, <30% leakage, resolved → keep 3 layers. EXPECTED: FAIL → deliverable
  = shallow + deep-bound, declared BEFORE compute is spent. (If even 2-layer shows heavy site↔fault
  leakage → design rethink, report first.)
- **S1 (cheap) _calS validation, one station, ALL windows.** Regenerate B011-or-B928 1to3/2to4/3to5
  (+4to7 if adopted) with dvv_roll30cal.py --origin-anchor (t_s=1.0). Differential step-injection
  (bench_step_recovery.py): GATE = all windows read injected amplitude within ~10% (vs 1.5:1 now);
  dt-vs-lapse of injected arm = pure slope through S-anchor, ~zero intercept.
  FAIL → anchor physics still wrong → HALT everything downstream.
- **S2 (medium, CPU-days) Fleet _calS regen**, all stations × windows on existing npz. Regen all
  families, filter at assembly. First reconcile canonical deseason (untracked
  lfe_features/deseason_dvv.py vs established _des pipeline) before any _calS_des.
- **S3 (medium) B011-p50 END-TO-END PILOT** — the right pilot (waveforms on disk + co-located
  independent instrument as control). Part-1-certified families ONLY → daily stacks (3-comp if
  feasible; polarization = free QC) → _calS dv/v on the window ladder → **LAPSE-DILATION test**
  (per family/epoch fit dt vs lapse-center: slope = distributed/deep-sensitive dilation, intercept =
  shallow static). This is the POSITIVE deep-sensitivity demonstration that makes the deep bound
  publishable rather than vacuous.
  GATES: (i) injection-recovery on this exact dataset (uniform stretch → slope only; early-window-only
  stretch → intercept only); (ii) seasonal must load on intercept/shallow channel; (iii) day/night
  dv/v arm agreement (b018_daynight pattern).
  FAIL signature: seasonal loads on SLOPE → decomposition broken OR real distributed seasonal —
  either way STOP and report; no fleet inversion with unexplained slope-seasonal.
- **S4 (cheap-med) Inversion upgrades, THEN fleet runs.** invert_multiwin.py: (i) carry n + scatter
  through monthly aggregation, weight LSQR rows 1/SE; (ii) save per-(cell,month) resolution masks;
  (iii) replace inverse-crime checkerboard with null_test.py residual-noise pattern + cross-kernel
  run (forward diffusion / invert single-scatter) to bound kernel error. Run on _calS_des certified
  input: deep channel = T2a-capable stations only; site-carrier families explicitly feed the
  site/shallow channel. Include the deferred 28-station-consistent-scale rerun (southern-attribution
  confound, MASTER_PLAN §2 FELL last item).
  GATES: null-test zero-field RMS vs observed per lat band; injected deep step recovered ≥~70%;
  injected shallow step leaks <~10% into fault field. Gate failure → claim DROPPED, not massaged.
- **S5 (cheap) Deliverables, honestly scoped:** margin-wide shallow/site dv/v field at
  station-spacing resolution; deep channel = time series + bound from certified, T2a-filtered,
  scale-consistent input + lapse-dilation sensitivity demonstration attached; intermediate layer only
  if S0 passed. Every map frame masked by S4 resolution masks.

---

## CUT / DEFER
- CUT intermediate layer from the promise unless S0 passes (expected: it won't).
- CUT classifier-P-on-stacks as certification anywhere. Picker survives as (a) family-level seed
  filter at discovery (thr~0.2-0.3) and (b) raw-feature character table.
- DEFER: per-sample U-Net (R7), multi-station picker, AE latents, 3-comp discovery templates —
  none gates Parts 1/2. DEFER Phase C detector rebuild until Parts 1-2 banked.
- DEFER fleet-wide picker-seeded re-discovery: only if B011-p50 beats catalog-seeded families on
  certified yield AND dv/v noise. No dedicated T2a re-download campaign — piggyback.

## DECISIVE-TEST LADDER (what to run first, all cheap, this week)
1. P1.2 B011↔PGC coincidence (free, lists on disk) — can falsify the whole certification concept.
2. S0 identifiability (CPU) — decides 2-layer vs 3-layer before any regen.
3. P1.3 T1-on-p50 (waveforms on disk) — the non-circular replacement for "274".
4. S1 _calS all-window injection at one station — go/no-go for the fleet regen.

## WHAT WOULD CHANGE THE VERDICT (Merlin)
- S0 shows resolvable intermediate → restore 3-layer deliverable.
- GOLD families fail PGC coincidence → T1 certifies station-local signal → priority flips to Phase C
  detector rebuild immediately.
- _calS windows still inconsistent after t_s=1.0 → stretch-anchor physics misdiagnosed → multi-window
  program halts regardless of certification quality.

*Merlin session: agentId a400fb391c20c3c45 (resumable). Key files he verified: invert_multiwin.py
(lines 53/63/86-101/124-127/145-150), dvv_roll30cal.py (19-29), daily_dvv_B928_2to4_calS.csv
(slope 1.491 vs _calT), mf_b011p50_*/mf_pgc_* on disk, ALL_families_labeled.csv.*
