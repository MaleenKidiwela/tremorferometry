# Track A (method) — results, 2026-09-16

Executes `results/merlin_review_plan_2026-09-16.md` steps 0–3, plus step 3b (user amendment). Everything
written here lives in `lfe_extension/`; `../lfe_catalog_repair` and `../codex-LFE` were read-only inputs and
every script used was COPIED into `scripts/` and edited there. The measurement chain is unchanged from the
validated E4–E6 path (20 Hz, 1–8 Hz Butterworth-4 zero-phase, existing component handling) — Merlin F6.

**Verdict per step**

| step | verdict |
|---|---|
| 0 — chain frozen, E6 re-scored, matched-FA ROC | **ACCEPT** (reproduces Merlin's table to 0.0005) |
| 1 — mergeB, 12 s, HF veto, FA-only thresholds | **ACCEPT** (recall at FA 1.0 = 0.2022 ∈ [0.18, 0.24]; site-rolled/data 0.90–1.03) |
| 2 — frozen cross-year transfer | **PASS both directions** (2012 → 2010 ratio 1.07 / 1.09; reverse 0.96 / 0.945 on held-out data) |
| 3 — template iteration | **REJECTED** (held-out 2010 recall at matched FA falls 0.012–0.017; adoption needed ≥ +0.01) |
| 3b — are quiet-day detections real? | **decided**: ≥ 86 % noise at T_ext; ~22 % real repeats at T_core |

Steps 4–6 were not started, per the brief.

---

## STEP 0 — measurement chain frozen, E6 re-scored, matched-FA ROC — **ACCEPT**

### (a) Merlin F4 confirmed: the on-disk E6 score JSONs were mislabeled
`e6_score_e6_{sites,mergeA,mergeB}.json` were written without `--e6`, so the key `L{L}_T{int(T)}` collided:
their `T7` is really **T = 7.5**, their `T8` is really **T = 8.5**, and **T = 8.0 was overwritten and absent**.
Re-scored with `--e6` into `results/e6_rescore/` (`scripts/a0_rescore_e6.sh`). Example — mergeB own recall /
quiet FA: old "T8" 0.1448 / 0.317 (= T 8.5); true T 8.0 is **0.1886 / 0.967**, which is the number Merlin
quoted from the FA tables.

Two things the spec did not anticipate:
* **9 scoreable passes, not 11.** `e6_lev_a_mv` and `e6_lev_c_mv` are build-half moveout passes run without
  `--budget`/`--catalog-out`; they have no budget CSV and no detections parquet.
* **The lever QC CSVs are partial and crash the scorer.** `e6_template_qc_2012_a.csv` holds only the 7 rebuilt
  stations, `e6_template_qc_2012_c.csv` only the 7 persistent ones (and has no `usable` column). The correct
  per-family station set was taken from the template npz actually used
  (`results/e6_rescore/qc/qc_from_npz_lev_{a,c}.csv`).

### (b) One-to-one matching (F11) — the delta, reported once
`results/e6_rescore/one2one_delta.csv`, 8 variants × 5 thresholds:
**own-family recall is identical to four decimals everywhere (Δ = 0.0000)**. Within-family declustering at the
template length already makes the detection side one-to-one, and no two Bostock events of one family share a
detection. The any-family CATALOG recall falls by 0.042–0.147 (mergeB at T = 8: 0.295 → 0.176) and its 600-s
chance falls with it (0.041 → 0.032). One-to-one greedy-by-time is the default for all Track A scoring from
here.

### (c) Matched-false-alarm ROC — reproduces Merlin's table
`scripts/a0_matched_fa_roc.py` → `results/e6_matched_fa_roc.csv` (grid: `results/e6_fa_recall_grid.csv`).
T*(FA) log-interpolated in FA on the {7, 7.5, 8, 8.5, 9, 10} budget grid; own recall interpolated linearly in
T (Merlin's recipe) and also recomputed exactly at T* (agrees to ≤ 0.002).

| variant | recall @ FA 1.0 | @ 0.5 | @ 0.1 | ETS own precision @ FA 1.0 | Merlin's value |
|---|---|---|---|---|---|
| sites (no merge) | 0.1653 | 0.1396 | 0.0923 | 0.421 | 0.165 / 0.140 / 0.092 · 0.42 |
| mergeA | 0.1711 | 0.1445 | 0.1051 | 0.473 | 0.171 / 0.144 / 0.105 · 0.47 |
| **mergeB** | **0.1901** | **0.1626** | **0.1169** | 0.452 | 0.190 / 0.163 / 0.117 · 0.45 |
| mergeB + cap PGC 25 % | 0.1669 | 0.1446 | 0.0988 | 0.473 | 0.167 / 0.145 / 0.099 · 0.47 |
| mergeB + 2-site coincidence | 0.1808 | 0.1563 | 0.1128 | 0.486 | 0.181 / 0.156 / 0.113 · 0.49 |
| lev_a (wide far-borehole pick) | 0.1866 | 0.1606 | 0.1156 | 0.451 | 0.187 / 0.161 / 0.116 · 0.45 |
| lev_b (≥ 5 dB stations) | 0.1993 | 0.1694 | 0.1246 | 0.441 | 0.199 / 0.169 / 0.125 · 0.44 |
| lev_c (2010 + 2012 templates) | 0.1989 | 0.1719 | 0.1248 | 0.456 | 0.199 / 0.172 / 0.125 · 0.46 |

Maximum deviation **0.0005** (bar ± 0.005) → **ACCEPT**; no scorer audit needed.
Figure: `figures/a0_matched_fa_roc.png`.

### (d) `--min-sites` (F8)
`e4_template_test.py --min-sites N` skips a family-day unless ≥ N distinct SITES have a usable template and
data; `n_sites` is now recorded on every budget row, netfam row and detection. Because the gate only decides
whether a family-day is scored (the surviving family's network sum is unchanged), it is exactly equivalent to
a post-hoc filter on `n_sites`, so the frozen point below is computed without the gate and its cost reported
separately: on 2012 / 15 stations **51 of 5,919 family-days (0.9 %)** are single-site, carrying 36 of 1,574
quiet-day detections at T = 8 and 1.4 % of all detections; exactly 1 of 116 families is single-site throughout.

---

## STEP 1 — variant, template length, HF veto, frozen operating point — **ACCEPT**

### (a) 12-s templates adopted
One GPU pass (`scripts/a1_run_mergeB_L12.sh`, 8.5 min on the H100), identical to the E6 mergeB pass except
`--lengths 12`. Own recall at matched quiet-day FA 1.0: **0.2022 (12 s) vs 0.1901 (8 s), Δ = +0.0121** ≥ the
pre-registered +0.01 → **adopt 12 s**. (FA 0.5: 0.1747 vs 0.1626; FA 0.1: 0.1225 vs 0.1169.)

### (b) HF-ratio veto (Merlin F14, flagged SPECULATIVE) — **adopted; the largest single gain**
`scripts/a1_hf_veto.py`: E(10–19 Hz)/E(2–8 Hz) at CN.PGC over the scored window, from the NATIVE-rate traces
(the 20 Hz cache cannot carry the band; PGC is BH 40 Hz in 2012 → Nyquist 20 Hz, so the band stops at 19 Hz).
PB.B011 (EH 100 Hz, literal 10–20 Hz) computed alongside as a band-choice control and agrees.

| | own-matched detections | quiet-day detections | separation |
|---|---|---|---|
| PGC median HF ratio (8 s) | 0.247 | 4.547 | 18× |
| B011 median HF ratio (8 s) | 0.163 | 0.999 | 6× |

Cut = 99th percentile over BUILD-half own-matched detections (3.511 at 8 s, 3.638 at 12 s); the test half never
touched it. **Quiet-day FA falls 61–65 % (8 s) / 58–60 % (12 s) at a test-half own-recall cost of 0.3–0.7 %
relative** — bar was ≥ 30 % at ≤ 2 %, and this also clears Merlin's "would change my mind" line (FA cut ≥ 50 %).
The quiet-day detections are HF-rich transients (uncatalogued micro-earthquakes / site glitches), not
sub-threshold LFEs. Figure: `figures/a1_hf_ratio.png`.

### (c) Thresholds by FALSE ALARMS ONLY
`scripts/a1_operating_point.py`, 0.25-MAD grid (exact, not interpolated: the detections parquet reproduces the
budget counts at every tabulated T). Quiet days = 2012 doy 228–242, verified **0 PNSN tremor detections in the
family box**.

| configuration | T_ext (quiet ≤ 1.0 AND post-SSE ≤ 1.5) | own recall | ETS own prec. | T_core (quiet ≤ 0.1) | own recall | ETS own prec. |
|---|---|---|---|---|---|---|
| mergeB 8 s, no veto | 8.50 (0.317 / 0.989) | 0.1448 | 0.574 | 9.00 (0.084) | 0.1125 | 0.669 |
| mergeB 8 s + HF veto | 8.00 (0.371 / 1.414) | 0.1874 | 0.467 | 8.75 (0.063) | 0.1267 | 0.628 |
| mergeB 12 s, no veto | 8.50 (0.386 / 1.179) | 0.1644 | 0.531 | 9.25 (0.078) | 0.1156 | 0.668 |
| **mergeB 12 s + HF veto (FROZEN)** | **8.25** (0.277 / 1.042) | **0.1841** | 0.488 | **8.75** (0.097) | **0.1463** | 0.596 |
| **what the veto buys** | −0.25 MAD | **+0.020** | | −0.50 MAD | **+0.031** | |

Recall by magnitude at T_ext: M < 1.6 0.161, M 1.6–1.8 0.183, M 1.8–2.0 0.187, M ≥ 2.0 0.258 — still nearly
magnitude-flat, as E5r found. **The post-SSE bar binds, not the quiet-day bar**: at T = 8.25 the unvetoed 12-s
run has quiet FA 0.555 (passes ≤ 1.0) but post-SSE FA 1.516, one part in 300 over the 1.5 bar. Track C should
expect the operating point to be set by residual post-SSE activity, not by the quiet background.
Figure: `figures/a1_fa_vs_T.png`.

### Accept bars
* own recall at matched FA 1.0 = **0.2022 ∈ [0.18, 0.24]** ✔
* quiet-day site-rolled / data = **0.903–1.033** over the whole threshold grid, inside [0.77, 1.3] ✔
  *Caveat:* the site-rolled control's detection TIMES are not stored, so the HF veto cannot be applied to the
  control; this bar is evaluated on the unvetoed detector.

### FROZEN OPERATING POINT — `results/frozen_operating_point.json`
variant **mergeB** (within-site w-weighted mean cc for {PGC, B009, B010, B011} and {LRIV, B005, B006, B007},
site weight √Σw²; no PGC cap, no coincidence gate) · **12-s templates** · **T_ext = 8.25**, **T_core = 8.75**
(MAD of the family-day network statistic) · **HF veto ON**, PGC 10–19 Hz / 2–8 Hz, cut **3.6384** ·
weights `e5r_weights_snr_2012.csv` (sha256[:16] a214c27d096cac97) · templates `e5r_templates_2012.npz`
(64b45a5bf6a34c1a) · moveout `e4_moveout_e5r_mv_build.csv` (78da98f68fbff432) · one-to-one matching ·
≥ 2-site rule implemented and its cost quantified. Every input hash is recorded in the JSON.

---

## STEP 2 — frozen cross-year transfer — **PASS (the make-or-break test)**

Track B had already delivered **PB.B009 / B010 for 2010 doy 222–256**, so the 2010 network is *identical* to the
2012 ablation network: nine stations / six sites (PGC + B009 + B010 + B011 = one site, VGZ, SNB, NLLB, B926,
LRIV). 2010 ground truth = Bostock 2010 events re-keyed so `fam_id = bkey` + the **90 accepted** offsets of
`e3_bostock_family_offsets.csv` → 15,931 events.

| | T_ext = 8.25 | T_core = 8.75 |
|---|---|---|
| own recall, 2012 ablation (test half) | 0.1983 | 0.1589 |
| own recall, 2010 frozen transfer | 0.2117 | 0.1724 |
| **ratio (bar 0.7–1.3)** | **1.067 PASS** | **1.085 PASS** |
| quiet FA 2012 / 2010 (strict-quiet) | 0.773 / 1.158 | 0.260 / 0.613 |
| FA ratio (bar 0.5–2.0) | 1.50 PASS | 2.36 not met |
| own precision 2012 / 2010 | 0.353 / 0.309 | 0.471 / 0.402 |

Same picture with the frozen HF veto applied to both years (recall ratios 1.066 / 1.086). **No time-base or
moveout problem exists** — 2012 templates, weights, moveout and MAD thresholds transfer unchanged.
Figure: `figures/a2_transfer_recall_by_mag.png`.

**The T_core FA ratio miss is an artefact of the 2010 quiet window, not of transfer.** Of Merlin's candidate
days (222–226, 250–256) only {223, 224, 226, 252, 254, 255, 256} are PNSN- and Bostock-quiet — three pre-SSE
and four post-SSE (doy 251 alone carries 91 tremor detections and 21 Bostock events). Split:

| window | quiet FA 2010 @ T_ext | ratio | @ T_core | ratio |
|---|---|---|---|---|
| pre-SSE (223, 224, 226) | 0.233 | **0.30** | 0.078 | **0.30** |
| post-SSE (252, 254, 255, 256) | 1.851 | 2.40 | 1.015 | 3.90 |

2010 pre-SSE is ~3× *quieter* than 2012 pre-SSE; the entire excess is post-SSE, the same behaviour that sets
T_ext in 2012.

**Reverse direction (2010-built templates → 2012).** No 2010 from-scratch templates existed, so they were built
with the same machinery (`e5r_templates.py` parameterised by `--year / --gt / --stations / --cache-dirs /
--out-dir`): **716 usable (family, station) templates over 90 families**, q-control cc to the SeisBench stacks
0.79 (PGC) / 0.79 (VGZ) / 0.87 (SNB) — as good as the 2012 build's 0.78.

| comparison | T_ext | T_core |
|---|---|---|
| all 2012 events, all families | 0.175 vs 0.262 → 0.67 | 0.142 vs 0.214 → 0.66 |
| all 2012 events, the 90 common families | 0.202 vs 0.269 → 0.75 | 0.164 vs 0.221 → 0.74 |
| **2012 TEST half, the 90 common families** | **0.2025 vs 0.2109 → 0.96** | **0.160 vs 0.170 → 0.945** |

The first two rows are contaminated (for the 2012-built templates the 2012 build half is training data). On
the held-out half and the families both sets cover the two directions agree within 5 % → **reverse transfer
PASS**. The gap in the "all families" row is pure coverage: only 90 of 116 families have an accepted 2010
offset and ≥ 30 2010 events.

---

## STEP 3b (user amendment) — are the quiet-day detections false alarms, or real inter-ETS repeats?

Quiet-day detections are an upper bound on the FA rate, not a measurement, because families also fire between
ETS episodes and Bostock's catalog (no events outside ETS windows) cannot say otherwise. Test: one GPU pass of
the frozen configuration with **CN.VGZ, CN.SNB and PB.B926 excluded from the network sum** (113 of 116
templated families have a usable template at one of the three), then stack that withheld station at the
family's quiet-day detection times on the predicted S and compare with 120 matched random-time stacks.
ETS-day own-matched detections are the positive control.
(`scripts/a3b_run_withheld.sh`, `scripts/a3b_quietday_stack_test.py` → `results/a3b_quietday_stack_test.csv`,
figure `figures/a3b_quietday_stack_test.png`.)

| | T_ext = 8.25 | T_core = 8.75 |
|---|---|---|
| families with ≥ 20 quiet-day detections (testable) | 42 | 25 |
| PASS / FAIL / inconclusive | **0 / 33 / 9** | **1 / 16 / 8** |
| median quiet-day S-SNR excess over null | +0.13 dB | +0.46 dB |
| median cc of the stack to the family template (null) | 0.168 (0.157) | 0.167 (0.170) |
| **ETS positive control** | **76 / 81 register**, median **+7.57 dB**, cc **0.735** | 70 / 75 |

Power is ample: the control's N (median 62) is comparable to the quiet-day N (median 52) and the null scatter
is 1.2 dB, so a real repeat population would be seen.

**Quantitative bound.** Calibrating each family's coherent signal-to-noise from its own ETS control
(excess = 10 log₁₀(1 + p²·N·(a/σ)²), p = 1 for the control) and solving for the real-repeat fraction p of the
quiet-day detections (`results/a3b_realfraction_bound.json`):

| threshold | families | pooled quiet excess | control | **p (real fraction)** | 95 % upper bound |
|---|---|---|---|---|---|
| T_ext = 8.25 | 36 | +0.104 ± 0.215 dB | +7.69 dB | **0.057** | 0.140 |
| T_core = 8.75 | 21 | +0.776 ± 0.296 dB | +7.58 dB | **0.225** | 0.325 |

**Verdict.** At **T_ext the quiet-day rate is a genuine false-alarm measurement** — ≥ 86 % of those detections
are noise — so the FA-based threshold rule stands as written and is not materially conservative. At **T_core
about 22 % (≤ 33 %) of the surviving quiet-day detections are real inter-ETS repeats**, so the true T_core FA
is ≈ 0.78× the measured one (0.097 → ≈ 0.075 per template-day): T_core is mildly conservative. Only one family
(bkey 75, 21 detections) passes individually at T_core and none at T_ext, so per-family promotion of inter-ETS
detections to events is **not** yet justified — the population estimate is the usable result, and it says the
inter-ETS LFE rate in these families is small but not zero. Caveats: 71 of 113 families have < 20 quiet-day
detections and are untestable; the 12-station withheld network is noisier on quiet days than the frozen
15-station one, so the bound errs conservatively.

---

## STEP 3 — template iteration (Merlin's protocol)

### Selection
Two GPU passes of the frozen configuration with the two-site coincidence gate (`--coinc-gate 2,2.5`), catalog
from T_core = 8.75, then tremor days only, the frozen HF veto and cross-family assignment:

| selection | detections | → tremor days & T_core | → HF veto | → cross-family | families | median events/family | ≥ 30 | ≥ 50 |
|---|---|---|---|---|---|---|---|---|
| gated, 15 stations | 10,135 | 8,977 | 8,736 | **6,468** | 115 | 42 | 68 | 49 |
| leave-one-site-out (PGC site removed from the sum) | 2,218 | 2,166 | 2,147 | **1,882** | 105 | **11** | 20 | 8 |

The PGC site is the only one over Merlin's 40 % weight-share trigger (it exceeds it on 53 % of family-days;
every other site's median share is ≤ 0.14), and its leave-one-site-out selection lands well under Merlin's
fallback threshold of 50 events per family per SSE. **Fallback taken: PGC-site templates are re-stacked only
for the 19 families with ≥ 30 LOSO events; the other families keep the frozen PGC-site templates**, because
re-stacking that site on events it selected itself is the self-reinforcement the protocol exists to prevent.

### Stacking
`e5r_templates.py` unchanged (rms-normalised windows, linear stack, cc alignment ±1.5 s, 2 iterations, one
shift per station for all three components), `--min-events 30`. Merged over the frozen set: **1,219 pairs —
426 frozen pairs replaced, 75 genuinely NEW (family, station) templates, 718 frozen kept.**

| | far stations (B001–B007, LRIV) | near / high-SNR (PGC site, B926, NLLB) | VGZ / SNB |
|---|---|---|---|
| median ΔS-SNR (iterated − frozen) | **+0.6 … +1.9 dB** | **−1.1 … −3.8 dB** | +0.05 / +0.18 dB |
| new templates added | 65 of 75 | 2 | 4 |

Iteration does the TRANSFER job Merlin predicted (F5), not the sensitivity job: it adds templates at the far
stations Bostock's catalog never constrained, and *costs* S-SNR at the near stations because the selected sets
are smaller (median 42 events, 47 for the PGC site) than Bostock's 115–160. Stability is good: |ΔS pick| median
0.05 s (62 of 426 move > 0.5 s, mostly B004), differential moveout changes 0 of 426 beyond 0.3 s, max
cross-family cc falls at every station (PGC 0.70 → 0.65, SNB 0.65 → 0.57), cc to the SeisBench stacks
0.797 → 0.785. Figure `figures/a3_iteration.png`.

### The six rejection tests — and a calibration problem the spec did not anticipate
Applied literally, the ABSOLUTE clauses of tests (a) "< 0.6", (d) "> 0.7" and (e) "> 50 % of the 1–8 Hz energy
in one 1-Hz bin" fire at essentially the same rate on the **FROZEN** templates (13 / 320 / 443 of 1,144 pairs)
as on the iterated ones (32 / 282 / 433 of 1,219): they describe this template class (LFE spectra peak at
1–3 Hz; families a few km apart share waveforms), not damage done by iterating. The script therefore reports
both the literal verdict and a **baseline-conditioned** one that keeps every relative clause as written and
applies each absolute clause only where the frozen template satisfied it.

| test | literal | baseline-conditioned |
|---|---|---|
| (a) q-control cc vs SeisBench | 32 | 19 |
| (c) S pick > 0.5 s | 62 | 62 |
| (c) differential moveout > 0.3 s | 0 | 0 |
| (d) cross-family collapse | 282 | 33 |
| (e) spectral narrowing | 433 | 28 |
| (f) S-SNR falls > 1 dB | 177 | 177 |
| **any (of 1,219 pairs)** | 709 | **261** |
| **families accepted** | **2 / 116** | **52 / 116** |

51 of the 64 rejected families are rejected by (f) alone — the small-N problem above, not drift.

### Rejection test (b) — held-out-year recall at matched false-alarm rate: **the iteration is REJECTED**
The 2012 test half is no longer clean (the iterated templates were built from detections on both halves), so
the decisive test is the held-out year. 2010, nine persistent stations, frozen vs iterated templates compared
at the SAME 2010 quiet-day FA (`scripts/a3_score_heldout.py` → `results/a3_heldout_matched_fa.csv`):

| matched 2010 quiet FA | T frozen | recall frozen | T iterated | recall iterated | **Δ** | Merlin's rule |
|---|---|---|---|---|---|---|
| 1.00 | 8.36 | 0.2021 | 8.36 | 0.1846 | **−0.0174** | reject (> 0.01 drop) |
| 0.50 | 8.91 | 0.1608 | 8.91 | 0.1466 | **−0.0141** | reject |
| 0.25 | 9.54 | 0.1257 | 9.54 | 0.1132 | **−0.0125** | reject |

Restricted to the **52 families the rejection tests accepted**, the iteration is exactly neutral
(Δ = −0.0003 / −0.0002 / −0.0004): no gain, no loss (`a3_heldout_matched_fa_acceptedfams.csv`). The same sign
appears even on the CONTAMINATED 2012 test half (frozen 0.1847 → iterated 0.1696 at T_ext, quiet-day
detections ×1.14), so this is not a held-out-year artefact.

Convergence (Merlin's stop rule cc(iter k, k−1) ≥ 0.98): cc between the iterated and frozen 12-s templates
over ±1 s is median **0.70** and **0 of 426 pairs reach 0.98** (`results/a3_cc_iter_vs_frozen.csv`) — the
templates changed a lot; they just did not change for the better.

### STEP 3 VERDICT — **iteration NOT adopted; the frozen Bostock-seeded templates stand**
Merlin's adoption rule (≥ +0.01 held-out recall at matched FA) is not met anywhere, and the whole-catalog
result triggers his explicit rejection clause. No second re-stack was run: it would be built on an already
rejected first iteration, and the expected ceiling was ≤ +0.02.

**Why it fails, in one line:** the PGC site carries > 40 % of most families' weight, and its leave-one-site-out
selection yields only ~11 events per family — so the stations that dominate the network statistic get *worse*
templates (−1 to −4 dB S-SNR), while the gains are confined to far stations that carry little weight.

**Recommended next experiment (not run):** a *targeted* iteration that re-stacks ONLY the far stations
(B001–B007, LRIV — the ones that gain 0.6–1.9 dB and account for 65 of the 75 new templates) and keeps the
frozen templates everywhere else. That isolates the part of the iteration that demonstrably works from the
part that does the damage, and needs one template merge plus two GPU passes. It is also the natural route to
Merlin's real target for A2 — coverage for the 14 template-less families and for new instrument epochs.
