# P3D — the dense north network, 2019–2026: transfer, recall and the discovery arithmetic

2026-09-19, `lfe_picker`, agent **P3D**. Spec: `results/merlin_review_p3d_2026-09-19.md`
(Steps 0–5; Step 6 only if Step 3 licenses it; Step 7 write-only).
Everything new is `scripts/s_*.py` → `results/p3d/` and `figures/s_*`. `q_*`, `p_*`, `g_*`, `o_*`,
`r_*` were used as READ-ONLY libraries (`g_common.cc_template`, `p2_cf.cf_day_sid`, `p_models`,
`o_cache_month.cache_one` via `r_cache_p0.dispatch`, `o_run_month`/`o_common` through an
in-process monkey-patch); `e5r_stations.csv` and every frozen file are byte-identical.
2011 and 2013 were never opened — a hard guard raises in `s_common.cache_index`.

---

## 0. Headline

| step | result |
|---|---|
| **0 premise correction** | reproduced Merlin's arithmetic exactly. Median family: sites ≤45 km **1 → 2**, ≤55 km **2 → 5**. **82/130** families gain ≥1 *genuinely new* site ≤45 km. Counting CLRS/OSQM as new sites (Merlin flaw 3) would have said 0 → 2 and 93/130. Of the ten ISOLATED hidden families only **12, 29, 62** gain a genuinely new site ≤45 km |
| **0 prediction (hashed first)** | predicted SYMB S-SNR(N60) **8.52 dB**, w 1.13 — i.e. *one more PGC-class site*. Predicted certified recall gain over families with ≥1 new site ≤45 km: **D1/D2 +17.2 %, D3 +27.9 %** |
| **0 acceptance** | SYMB SSE delivery **99.8 %** (bar 80 %) → **PASS**. 74 station-days masked at coverage_frac < 0.9 |
| **1 machinery control** | **CLRS vs B926: 44/52 control-OK families match, 0.846 (95 % CI 0.72–0.93) → PASS** (bar ≥0.70). OSQM vs B001 3/5 = 0.60, 95 % CI 0.15–0.95 → **UNTESTABLE, low power**, not a machinery failure |
| **2 transfer** | **(family, SYMB) accepted = 54** (proceed bar ≥30; FAIL bar <10) → **PROCEED**. Bank: 251 (family, station) templates over 11 stations, 184 in the production bank |
| **2 — THE key number** | **ρ_SYMB / ρ_PGC = 0.783** (n = 64 families, same events, leave-one-out single-event registration). ρ_BUCKS/ρ_PGC = **0.296** (n = 33). No new station reaches PGC. Merlin's decision scale: ≥3 × → discovery has a chance; ≤1.2 × → "SYMB is not the prize" |
| **3 pair-rate prediction** | δ_old7 **1.25 σ**, δ_dense **1.88 σ** → ratio **1.50** → predicted δ_full **0.57 σ** vs the 1.5 σ bar → **STEP 6 NOT LICENSED**. PGCq ALONE is the best configuration (2.44 σ); every other site dilutes it. § 4 |
| **4a CF timing kill test** | **SURVIVES.** Median \|CF peak − Crosson S\|: SYMB **0.95 s**, PGCq 1.05, VGZ 1.20, BUCKS 2.30; 1 of the 4 ≤60 km sites above the 1.5 s bar. Far sites 2.4–2.85 s = uniform chance (2.50 s). **SYMB — a station the CF has never seen — has the best residual of all 16 sites.** § 5 |
| **5 certified recall gain** | **BAR NOT MET.** At matched quiet-day FA (18.28/day; T 9.25 → 8.96): certified gain **−5.1 %** (bootstrap LB −8.5 %), Spearman(gain, Δsites) **−0.013**, against the hashed prediction **+27.9 %**. Marginal purity 0.361 vs reverse-marginal 0.306. § 7 |
| **PQ amendment (V1/V2)** | hashed before any PQ stack existed. Premise: median sites ≤45 km 1 → 2 (V1) → **4 (V2)**. JRDRB delivers 87 days against 553 nodata; LHLYB only 19/59 of its build episode. The S18 V2 run "licenses" Step 6 — **but that is a null-realisation artifact**: 5 draws of the identical statistic on identical events give ratios 1.72–5.36 across a bar of 4. § 8 |

---

## 1. Step 0 — roster, sites, the premise correction, and the hashed freeze

`scripts/s_step0_roster.py`, `scripts/s_common.py`; products `results/p3d/{stations_p3d.csv,
sites_p3d.csv, premise_family_sites.csv, pred_q2iii.csv, pred_q2iii_recall.csv,
pred_site_terms.csv, freeze.json, step0_delivery.csv, masked_station_days.csv, step0_accept.json}`.

### 1.1 Roster and the co-located site merge

27 stations → **17 sites**. The frozen 15 keep their `STAS_ORDER` positions; the ten new 3-C
stations and the two Z-only helpers are APPENDED after them. Z-only STW/UWFH carry
`summable = False` and never enter a site statistic. The merge Merlin's flaw 3 requires:

| site | members | separation | median hypocentral distance to the 130-family field |
|---|---|---|---|
| PGCq | PGC + B009 + B010 + B011 | — (frozen) | 52.7 km |
| LRIVq | LRIV + B005 + B006 + B007 | — (frozen) | 65.0 km |
| **B926q** | **B926 + CLRS** | **25 m** | 60.20 km — *withheld* |
| **B001q** | **B001 + OSQM** | **17 m** | 77.4 km |

LRIV carries `band_p3d = HH` and is a transfer TARGET, not a template station (its 2012 templates
failed A-2 step 4 and the frozen policy drops it after 2019-05-30).

### 1.2 The premise correction (the table that corrects the P3D plan)

Sites within 45 / 55 / 60 km of each family, BEFORE (every site that already existed — which
includes B926q and B001q, since CLRS and OSQM add no geometry) vs AFTER:

| | median before | median after |
|---|---|---|
| sites ≤45 km | **1** | **2** |
| sites ≤55 km | **2** | **5** |
| sites ≤60 km | **3** | **7** |

Genuinely new sites ≤45 km per family: 48 families gain 0, 50 gain 1, 32 gain 2. The ten
**isolated hidden** families, which are what a discovery test would be scored on:

| bkey | sites ≤45 km before → after | new sites ≤45 km | new sites ≤55 km | SYMB distance (km) |
|---|---|---|---|---|
| 8 | 1 → 1 | — | LOPEZ | 76.0 |
| 10 | 3 → 3 | — | BUCKS+ODUC | 57.3 |
| 12 | 2 → 4 | **BUCKS+ODUC** | BUCKS+ODUC+SYMB | 45.4 |
| 15 | 2 → 2 | — | BUCKS+SYMB | 54.5 |
| 29 | 1 → 3 | **BUCKS+SYMB** | BUCKS+SYMB | 41.2 |
| 62 | 1 → 2 | **SYMB** | BUCKS+SYMB | 36.7 |
| 87 | 2 → 2 | — | RNWY+SYMB+TURTL | 51.3 |
| 113 | 0 → 0 | — | GOBB | 57.0 |
| 115 | 1 → 1 | — | SYMB | 51.5 |
| 119 | 0 → 0 | — | GOBB | 59.0 |

This is Merlin's list, independently reproduced: only 12, 29, 62 gain a genuinely new near site;
113 and 119 gain nothing inside 55 km except GOBB; 15 and 87 gain nothing inside 45 km.

### 1.3 The a priori geometry rule

Sites whose median hypocentral distance to the family field is ≤60 km, withheld sites dropped:
**PGCq (52.7), VGZ (52.7), SYMB (44.7), BUCKS (57.4)**. B926q lands at 60.20 km and is withheld
anyway. Merlin's illustrative drop list names NLLB, B004, B001/OSQM, GOBB, LOPEZ, TURTL, RNWY; the
rule applied literally also drops **B003 (70.1), LRIVq (65.0), SNB (68.0), PTRF (66.8), ODUC
(73.5)**. SYMB is, by this measure, the best-placed site in the whole roster — better than the PGC
cluster.

### 1.4 The Q2(iii) prediction, hashed before any transfer stack was opened

S-SNR = a + b·log₁₀(hypocentral km), fitted per station class on `e5r_template_qc_2012.csv` at
matched stack size (`snr_N30`, `snr_N60`); weights w = 10^(S-SNR/20) with the frozen 3×-median clip;
recall through the E5 prediction model read as a library of formulae (bg + e_net(e_class(epi,
klass)), pooled log-spread), matched FA read as matched T in family-day MAD units.

| station | median hypocentral km | predicted S-SNR(N60) dB | predicted w |
|---|---|---|---|
| **SYMB** | 44.7 | **8.52** | **1.13** |
| BUCKS | 57.4 | 7.86 | 1.04 |
| CLRS | 60.2 | 7.73 | 1.03 |
| RNWY | 64.4 | 7.55 | 1.01 |
| LRIV | 65.0 | 7.53 | 1.01 |
| PTRF | 66.8 | 7.46 | 1.00 |
| GOBB | 68.8 | 7.38 | 0.99 |
| TURTL | 70.9 | 7.30 | 0.98 |
| ODUC | 73.5 | 7.20 | 0.97 |
| OSQM | 77.4 | 7.07 | 0.95 |
| LOPEZ | 77.6 | 7.06 | 0.95 |

**Predicted certified recall gain** (families with ≥1 genuinely new site ≤45 km):
**D1 +17.2 %, D2 +17.2 %, D3 +27.9 %** (mean +20.8 %).

The model's own honest uncertainty is written next to it: the per-station residuals (SITE TERMS)
of the same fit span **6.5 dB**, from LRIV −3.7 dB to PGC +2.7 dB, and a class-plus-distance model
cannot know a new site's own term. § 3.2 shows this is exactly where the prediction went wrong.

### 1.5 Freeze and acceptance

`results/p3d/freeze.json` (sha256 in the file) records: the 30 hidden families and the
`holdouts.json` hash, the withheld SITE set {B926q = B926+CLRS, GOBB, LOPEZ, TURTL} and the rule
that they never enter any network sum, the site merge, the appended `STAS_ORDER`, the episode split
(calibrate S01–S13, test S14/S18/S21/S22, BUCKS/TURTL/ODUC build on S18 and test on S21/S22), the
network epochs D1/D2/D3, the quiet months, the geometry rule with every site's median hypocentral
distance, the **discovery threshold rule** (the null quantile giving ≤30 expected false edges over
the realised pair count, ≈1e-9 at 3.4×10¹⁰ pairs — the RULE is frozen, the quantile is recomputed
from the run), the Step-6 licence rule, and the sha256 of every table including the prediction.

**Stores.** The 20 Hz / 1–8 Hz cache was built for every delivered P3D station-day through the
frozen chain (`scripts/s_cache.py` calling `r_cache_p0.dispatch` = `o_cache_month.cache_one`, and its
Z-only twin for STW/UWFH); ~9,000 new station-days. `scripts/s_build_mm.py` is the P3D twin of
`p_build_mm.py` with the S15 filter replaced by `s_common.cache_index` (the override Merlin asks
for) and its own `_p3d` tag so the frozen `results/mm/day_index.parquet` / `days_zne_f16.dat` are
never touched — it is **written but NOT run**, because every P3D step reads the caches directly
through `s_common.load_day_zne` (same ZNE rotation, same StationXML azimuths, same values) and a
75 GB second memmap bought nothing this session. A hard guard in `s_common.cache_index` raises on any
attempt to open 2011 or 2013.

**Acceptance.** SSE-window delivery, from the downloader's manifest (read-only):

| station | SSE days in epoch | delivered | frac | masked (coverage < 0.9) |
|---|---|---|---|---|
| **SYMB** | 866 | 864 | **0.998** | 6 |
| CLRS | 866 | 863 | 0.997 | 10 |
| TURTL / RNWY / OSQM / ODUC | 419 / 598 / 598 / 270 | all | 1.000 | 0 |
| BUCKS | 325 | 308 | 0.948 | 1 |
| LOPEZ | 547 | 508 | 0.929 | 3 |
| PTRF | 866 | 782 | 0.903 | 14 |
| GOBB | 866 | 542 | **0.626** | 18 |
| STW / UWFH (helpers) | 866 | 838 / 794 | 0.968 / 0.917 | 20 / 2 |

**SYMB gate PASS (0.998 ≥ 0.80).** 74 station-days are masked and excluded from every later step
(the mask is enforced inside `s_common.cache_index`). GOBB's 63 % is a genuine data gap; GOBB is a
withheld site and is never in a network sum.

---

## 2. Step 1 — the co-located machinery controls

`scripts/s_step1_control.py` → `results/p3d/step1_colocated.csv`, `step1_summary.csv`.
The borehole templates are rotated to ZNE (`rot_matrix`, StationXML azimuths) and re-normalised
per component before any comparison (Merlin flaw 6). The statistic is `g_common.cc_template` —
exactly the statistic whose per-station cross-family p95 is tabulated in
`results/gates/g_crossfam_cc.csv`.

| new | co-located ref | bar (ref cross-family p95) | families tested | control-OK | match | fraction | 95 % CI | verdict |
|---|---|---|---|---|---|---|---|---|
| **CLRS** | B926 (25 m) | 0.348 | 57 | **52** | **44** | **0.846** | 0.719–0.931 | **PASS** |
| OSQM | B001 (17 m) | 0.219 | 9 | 5 | 3 | 0.600 | 0.147–0.947 | UNTESTABLE (low power) |

The CLRS leg is decisive and it validates the bar as well as the machinery: the cross-family chance
level measured *at CLRS itself*, from its own transfer stacks against the other families' rotated
B926 templates (3,192 pairs), is **p95 = 0.3435** — within 1 % of B926's tabulated 0.348.

The OSQM leg cannot say anything. B001 is the weakest station in the network (47 usable 2012
templates, cross-family p95 0.219, median 2012 stack S-SNR 4.19 dB) and OSQM only starts
2021-03-06, so its calibration episodes are S09–S13; just 9 families reach N ≥ 30 and only 5 pass
the control leg. P(X ≤ 3 | n = 5, p = 0.70) = **0.47** — no evidence against the bar. The five
control-OK families track B001 cleanly (OSQM cc ≈ 0.6 × B001 cc at every one), and the two "misses"
sit at 0.203 and 0.156 against a 0.219 bar. **Verdict: the machinery is confirmed on the powered
leg; the OSQM leg is recorded as UNTESTABLE, and Step 2 proceeds.** (Stopping the phase on an n = 5
leg while the n = 52 leg passes at 0.85 would be the error, not the reverse.)

---

## 3. Step 2 — transfer stacks on the calibration episodes

`scripts/s_transfer.py` (engine) + `scripts/s_step2_transfer.py` →
`results/p3d/{step2_legs.csv, step2_summary.csv, step2_rho.csv, bank_transfer_v1/}`.
Protocol exactly Merlin Q1; the frozen `cut_windows` / `align_stack` are COPIES, not re-derivations.

### 3.1 Coverage

| station | site | attempted | insufficient | untestable | **accepted** | rejected | accepted under the LITERAL (b) | median stack S-SNR dB | median excess dB | median ρ_LOO | median half-split cc |
|---|---|---|---|---|---|---|---|---|---|---|---|
| CLRS | B926q | 130 | 66 | 4 | **56** | 4 | 56 | 12.18 | 9.52 | 0.313 | 0.672 |
| **SYMB** | SYMB | 130 | 66 | 4 | **54** | 6 | 53 | 13.28 | 9.71 | 0.308 | 0.603 |
| GOBB | GOBB | 130 | 76 | 4 | 39 | 11 | 41 | 8.03 | 6.58 | 0.208 | 0.419 |
| LRIV (HH) | LRIVq | 130 | 66 | 5 | 37 | 22 | 40 | 5.61 | 4.15 | 0.142 | 0.214 |
| TURTL | TURTL | 130 | 97 | 3 | 23 | 7 | 20 | 6.59 | 4.03 | 0.202 | 0.311 |
| PTRF | PTRF | 130 | 74 | 4 | 20 | 32 | 22 | 4.07 | 3.32 | 0.108 | 0.181 |
| RNWY | RNWY | 130 | 112 | 1 | 9 | 8 | 9 | 7.27 | 3.81 | 0.200 | 0.232 |
| BUCKS | BUCKS | 130 | 97 | 3 | 8 | 22 | 8 | 3.41 | 2.33 | 0.116 | 0.152 |
| OSQM | B001q | 130 | 112 | 1 | 4 | 13 | 4 | 2.14 | 2.19 | 0.075 | 0.154 |
| ODUC | ODUC | 130 | 97 | 3 | 1 | 29 | 1 | 2.96 | 2.37 | 0.050 | 0.133 |
| LOPEZ | LOPEZ | 130 | 125 | 0 | 0 | 5 | 0 | 3.65 | 2.63 | 0.143 | 0.214 |

**Step-2 gate: 54 (family, SYMB) templates accepted ≥ 30 → PROCEED** (and far above the FAIL bar
of 10). `insufficient` dominates because N ≥ 30 hp events per (family, station) inside the
station's own calibration episodes is a hard bar for a 115-family catalogue whose detections are
concentrated in a few large families; the labelled fallback (cat_core, p_real ≥ 0.5, N ≥ 60) is
included in these counts. BUCKS/TURTL/ODUC build on S18 alone, which is why their `insufficient`
counts are the highest.

Banks (frozen npz format, `tpl` / `sp` / per-station ±20 s stacks):
`bank_transfer_v1/bank_full.npz` **251** (family, station) templates of which **67 hidden**, and
`bank_transfer_v1/bank_production.npz` **184** with **0 hidden** (verified by key listing in
`bank_*_keys.csv`).

### 3.2 THE key number — single-event registration

Median leave-one-out per-event correlation of each event window against its own stack over
S−1..S+7 s, on the SAME events at the SAME N as the reference station:

| new station | n families | **ρ_new** | ρ_PGC | **ρ_new / ρ_PGC** | ρ/ρ_VGZ | S-SNR new dB | S-SNR PGC dB |
|---|---|---|---|---|---|---|---|
| **SYMB** | 64 | **0.308** | 0.412 | **0.783** | 1.198 | 13.28 | 16.01 |
| CLRS | 64 | 0.313 | 0.412 | 0.773 | 1.180 | 12.18 | 16.01 |
| GOBB | 54 | 0.208 | 0.411 | 0.571 | 1.033 | 8.03 | 15.34 |
| TURTL | 33 | 0.202 | 0.396 | 0.559 | 0.949 | 6.59 | 13.65 |
| RNWY | 18 | 0.200 | 0.461 | 0.488 | 0.722 | 7.27 | 16.00 |
| LOPEZ | 5 | 0.143 | 0.470 | 0.398 | 0.420 | 3.65 | 17.05 |
| LRIV (HH) | 64 | 0.142 | 0.412 | 0.366 | 0.523 | 5.61 | 16.01 |
| **BUCKS** | 33 | **0.116** | 0.396 | **0.296** | 0.553 | 3.41 | 13.65 |
| PTRF | 56 | 0.108 | 0.408 | 0.275 | 0.452 | 4.07 | 16.10 |
| OSQM | 18 | 0.075 | 0.461 | 0.168 | 0.250 | 2.14 | 16.00 |
| ODUC | 33 | 0.050 | 0.396 | 0.158 | 0.253 | 2.96 | 13.65 |

**ρ_SYMB / ρ_PGC = 0.78.** Merlin's decision scale reads: ≥3 × → the pairwise discovery arm has a
real chance; 1.5–2 × → recall gains of 20–40 % for SYMB-adjacent families and discovery stays
closed; **≤1.2 × → "SYMB is not the prize"**. The measurement is *below one*: the best-placed site
in the roster, 8 km closer to the family field than the PGC cluster, registers single LFEs
**worse** than PGC does.

**Why, and what it says about the Step-0 prediction.** On the same events at matched N, SYMB's
stack S-SNR is **2.73 dB BELOW** PGC's, where the hashed class-plus-distance model predicted it
**0.36 dB above** — a 3.1 dB miss, which is exactly the size of the per-station SITE TERMS the same
fit measures (§ 1.4: −3.7 dB at LRIV to +2.7 dB at PGC). SYMB's own site term is ≈ −3 dB. This is
the honest reading of the prediction: distance and station class predict the *ordering* of the new
stations well (SYMB ≈ CLRS > GOBB ≈ TURTL ≈ RNWY > PTRF ≈ BUCKS > OSQM ≈ ODUC tracks hypocentral
distance with two inversions) but they cannot predict a new site's noise, and the level is wrong by
the full site-term spread.

### 3.3 Controls

* **(c) matched-N same-events control.** Every leg carries the identical stack at the nearest of
  PGC / VGZ / B926 on exactly the same events at the same N. Legs whose control fails are labelled
  `untestable`, never `reject` — 1–5 per station.
* **(d) self-consistency.** An alignment-free re-stack at the template's own placement, vs the
  template: median 0.60 (SYMB) / 0.67 (CLRS), collapsing to 0.13–0.21 at the far stations.
* **PGC-dominance half-split.** The hp set is split at the median PGC-site z (from Track E's
  per-site re-run, `persite_<month>.parquet`) and the new-station stack is rebuilt from each half.
  Median half-split cc: **SYMB 0.603, CLRS 0.672** — comfortably above the 0.5 default bar and
  above each station's own cross-family chance; **BUCKS 0.152, ODUC 0.133, OSQM 0.154** are at or
  below chance, so those stations' pooled stacks are NOT certified against PGC dominance.
* **Weights.** `w_new = w_ref × 10^((S-SNR_new − S-SNR_ref)/20)` with both S-SNRs on the same
  events at the same N, then the frozen 3×-median clip and median normalisation.

### 3.4 A trap in acceptance test (b), reported in full

The wrong-cell moveout null at a SINGLE station degenerates: a cell whose predicted S is (S−P) ≈
6–7 s EARLY puts the family's **own P arrival** inside the S−1..S+5 window. Measured at PGC on
family 10, the excess is +16.2 dB at Δ = 0 and **+17.2 dB at Δ = −5 s**, falling to ≤2 dB for
|Δ| ≥ 10 s and to −15 dB at Δ = +8 s. So the literal (b) can fail at the gold-standard reference
station for a purely geometric reason. Both are reported for every leg: the literal `pass_b`, and
`pass_b_noP` — the same null with offsets within 4 s of −(S−P) removed. On the real grid-derived
offset distribution the two agree closely (SYMB 53 vs 54 accepted), so nothing material turns on it
here; the acceptance column uses `pass_b_noP` and the literal count is printed beside it.

---

## 4. Step 3 — the pair-rate positive control and the discovery PREDICTION

`scripts/s_pairlink3.py` (the P3 `q_pairlink2` statistic, re-pointed at the 2019–2026 network and at
hp candidates) + `scripts/s_step3_delta.py` → `results/p3d/{pairlink3_cal_null.json,
pairlink3_cal_true.json, pairlink3_cal_true_by_family.csv, step3_delta.csv, step3_delta_fit.csv,
step3_verdict.json}`.

Candidates: **2,834** hp events (p_real ≥ 0.8) of the **KEPT** (non-hidden) families on the
calibration episodes S01–S13, on days cached at all 14 stations of the two configurations;
**4.01 × 10⁶** pairs after the 27-s overlap exclusion, of which **94,138 same-family**.
Null: the identical chain with station *i*'s data circularly shifted by *i* × 1,800 s.
BUCKS is absent from the calibration episodes (it starts 2023-09) and drops out of the dense
configuration by the pre-registered availability rule, so the calibration dense set is
**PGCq + SYMB + VGZ**.

### 4.1 Per-pair rates and enrichment

| configuration | sites | null 1e-3 thr | same-family rate @1e-3 | @1e-4 | far-pair rate @1e-3 | enrichment @1e-4 |
|---|---|---|---|---|---|---|
| **old7** | B001q B003 B004 LRIVq NLLB PGCq VGZ | 11.56 | **0.0315** | **0.00617** | 1.47e-3 | 26 × |
| **dense (≤60 km rule)** | PGCq SYMB VGZ | 7.64 | **0.0753** | **0.0278** | 1.19e-3 | **155 ×** |
| dense − SYMB | PGCq VGZ | 6.12 | 0.0559 | 0.0193 | 1.49e-3 | 84 × |
| old7 + SYMB | the 7 + SYMB | 12.38 | 0.0554 | 0.0144 | 1.44e-3 | 59 × |
| PGCq alone | PGCq | 2.49 | 0.4930 | 0.1306 | 1.04e-1 | 20 × |
| SYMB alone | SYMB | 5.39 | 0.0145 | 0.00323 | 2.86e-4 | 106 × |

### 4.2 δ fitted from the tail quantiles (P3 § 5.5 method) and the pre-registered decision

| configuration | δ (σ) |
|---|---|
| **PGCq alone** | **2.44** |
| **dense = PGCq+SYMB+VGZ** | **1.88** |
| dense − SYMB = PGCq+VGZ | 1.72 |
| old7 + SYMB | 1.58 |
| **old7** | **1.25** |
| SYMB alone | 1.08 |
| VGZ alone | 0.82 |
| LRIVq alone | 0.73 |
| B004 / NLLB / B001q / B003 alone | 0.53 / 0.48 / 0.40 / 0.16 |

**ratio δ_dense / δ_old7 = 1.50** → **predicted δ_full = 0.38 σ × 1.50 = 0.57 σ**, against the
pre-registered licence bar of **1.5 σ (ratio ≥ 4)**.

> **Step 6 is NOT licensed on the V1 network.** The prediction IS the result.

Two things the ablation says that matter more than the headline:

1. **The pair statistic is maximised by the PGC cluster ALONE** (δ = 2.44 σ), and *every* other site
   dilutes it: adding VGZ takes it to 1.72 σ, adding SYMB and VGZ to 1.88 σ, and the full seven-site
   sum down to 1.25 σ. This is P3 § 5.4's dilution result reproduced at single-site granularity on a
   completely different candidate set and a different decade.
2. **SYMB does help, by about the amount its registration says it should.** It lifts PGCq+VGZ from
   1.72 → 1.88 σ (+9 %) and the seven-site sum from 1.25 → 1.58 σ (+26 %). That is a real,
   measurable contribution — and it is an order of magnitude short of the factor 4 the discovery
   arm needs.

The absolute δ values here are much larger than P3's 0.38 σ because the candidates are hp
detections of a template catalogue (strong, already-registered events), not Bostock's full
population. Merlin's rule uses only the RATIO for exactly this reason, and the ratio is what is
reported.

**DEFERRED for the coordinator's PQ amendment.** The Step-3 DECISION is to be taken on the V2
network (V1 + PQ.LHLYB/JRDRB/ALBHB/DAOB). At the time of writing PQ has delivered **0 manifest rows
and 0 cached station-days** (recorded in `freeze.json → amendment_v2.delivery_state_at_hash_time`),
so the V2 run is HELD, exactly as instructed. Note the constraint the split imposes: JRDRB (2021-04)
and ALBHB (2022-01) exist only from S09 / S11 onward, and **LHLYB (2023-11) and DAOB (2024-02) do
not exist in ANY calibration episode** — a V2 pair run that contains all four PQ stations can only
be built on S18, their pre-registered BUILD episode. The planned V2 readouts are therefore
(i) S11–S13 with JRDRB + ALBHB, and (ii) S18 with all four, with the pre-registered rule applied to
(ii). For the ratio to reach 4 the four PQ sites would, between them, have to do about
**2.7 × more than SYMB's +9 %** — which the V1 ablation says is only possible if their per-site δ is
comparable to PGCq's own 2.44 σ, i.e. if they are PGC-class sites and not SYMB-class ones.

---

## 5. Step 4a — the CF timing kill test: **the arm SURVIVES**

`scripts/s_step4_cf.py` → `results/p3d/{step4a_residuals.csv, step4a_summary.csv,
step4a_verdict.json}`; figure `figures/s_cf_timing.png`.

At every hp event of the four TEST episodes (5,275 events), the frozen P2 characteristic function
(`results/p2/models/p2nnpu.pt`) is read at each 3-C station and its peak time inside ±5 s of the
Crosson-predicted S is recorded. **Every station goes through the UNKNOWN-station embedding**
(`sid = nsta = 15`): the code path is `p2_cf.cf_day_sid(model, meta, x, sid)` with `sid` passed
explicitly — `sid_of` is consulted only under `--sid known`, so nothing is ever indexed by S15
position and a new station is legal. Days are read from the 20 Hz caches through `s_common`
(ZNE-rotated), not through `p_common.DayStore`, whose index is filtered to S15. Only the window grid
that covers each event's ±5 s is evaluated, by the same tiling rule `cf_day_sid` uses, so the values
are the full-day CF sample-for-sample in the interior.

| site | ≤60 km | n | **median \|residual\| (s)** | median residual (s) | frac ≤1 s | frac ≤0.75 s | median CF peak |
|---|---|---|---|---|---|---|---|
| **SYMB** | ✓ | 5,261 | **0.95** | −0.25 | 0.515 | 0.434 | 0.676 |
| **PGCq** | ✓ | 21,051 | **1.05** | −0.70 | 0.495 | 0.393 | 0.723 |
| **VGZ** | ✓ | 4,249 | **1.20** | −0.40 | 0.450 | 0.361 | 0.594 |
| B926q | | 9,663 | 1.50 | −0.55 | 0.361 | 0.275 | 0.563 |
| TURTL | | 5,263 | 1.50 | −0.30 | 0.366 | 0.293 | 0.512 |
| PTRF | | 4,896 | 1.70 | +0.05 | 0.338 | 0.268 | 0.490 |
| LRIVq | | 20,943 | 1.75 | +0.25 | 0.322 | 0.252 | 0.516 |
| GOBB | | 3,657 | 1.80 | 0.00 | 0.319 | 0.258 | 0.516 |
| LOPEZ | | 5,263 | 1.85 | +0.10 | 0.305 | 0.243 | 0.482 |
| NLLB | | 5,267 | 1.90 | −0.65 | 0.292 | 0.231 | 0.473 |
| RNWY | | 5,262 | 2.20 | −0.35 | 0.262 | 0.204 | 0.238 |
| **BUCKS** | ✓ | 3,761 | **2.30** | +0.75 | 0.227 | 0.166 | 0.418 |
| B001q | | 10,529 | 2.40 | +0.65 | 0.231 | 0.180 | 0.418 |
| B003 | | 5,266 | 2.40 | +0.10 | 0.234 | 0.178 | 0.395 |
| ODUC | | 3,761 | 2.80 | +1.80 | 0.157 | 0.115 | 0.434 |
| B004 | | 5,266 | 2.85 | +1.55 | 0.176 | 0.130 | 0.336 |

**Verdict: 1 of the 4 sites inside the 60 km rule exceeds 1.5 s → the arm is NOT closed.**
A CF carrying no timing information peaks uniformly inside the ±5 s search, i.e. median |residual|
= **2.50 s**: B004 (2.85), ODUC (2.80), B001q and B003 (2.40) are at or beyond chance, so the CF has
**no** usable timing at the far sites — the dilution result again, now in the timing domain.

Three things must be said with it, because they bound what 4b could deliver:

1. **0.95–1.20 s is above Merlin's own "5 km plausible" line of 0.75 s.** At 3.7 km s⁻¹ a 1.0 s
   single-station timing error is ~3.7 km; with the four ≤60 km sites (one of which, BUCKS, is at
   chance-adjacent 2.3 s), a back-projection median location error of ≤10 km is plausible and ≤5 km
   is not. Clause (c) of the reinstatement rule is therefore reachable and clause (d) — "≥3/10
   isolated hidden families with ≥20 matched detections within 5 km" — is not, for the same
   truth-scarcity reason Merlin recorded in flaw 2.
2. **The reading is optimistic for the frozen stations.** The CF was trained on 2012 data recorded at
   PGC/VGZ/NLLB/… so their residuals are not out-of-sample in station identity (only in time). SYMB,
   BUCKS, TURTL, GOBB, LOPEZ, PTRF, RNWY, ODUC are genuinely unseen stations, and SYMB — an unseen
   station — has the BEST residual of all sixteen sites (0.95 s). That is the strongest single piece
   of evidence in this report that the CF generalises across stations.
3. **Step 4b (back-projection) is licensed but was NOT RUN** in this session; it needs full-day CFs
   over the test episodes and the quiet months, the 1,444-cell grid at three false-event rates, and
   the shuffled-moveout/N-shift nulls. It is the next thing to run, and clause (d) should be
   expected to fail on truth availability rather than on method.

---

## 6. Step 7 — Arm B pre-registration (WRITE ONLY)

Because Step 3 did not license Step 6, the Arm B (seed-and-grow, Bostock's construction) protocol
has been written and hashed before any of it is executed:
`results/p3d/arm_b_preregistration.md`, sha256 in `freeze.json → arm_b`. Its design point is the
Step-3 ablation: the arm does not need a large per-pair δ, it needs one loud window per family per
episode plus record length, and the seed statistic is the MEDIAN per-site z (not the sum) precisely
because Step 3 showed the sum is PGCq alone. Its honest expected outcome is written into the file:
**Arm B will rediscover LARGE families and will not settle the isolated-hidden score.**

---

## 7. Step 5 — the recall re-run at matched FA: **the bar is NOT MET**

`scripts/s_run_month.py` (an in-process monkey-patch of the FROZEN `o_run_month`; `o_run_month.py`
and `o_common.py` are byte-identical) + `scripts/s_step5_analyse.py` →
`results/p3d/trackS/` and `results/p3d/{step5_tiers.csv, step5_family_counts.csv,
step5_a3b_certification.csv, step5_gain_by_family.csv, step5_verdict.json}`.

`o_run_month` derives its station list from the TEMPLATE BANK
(`STAS = C.order_stations({k[1] for k in tpl})`), so supplying the transferred bank brings the new
stations into the network automatically; the patch only teaches the supporting tables about them
(station table, epochs, the **co-located site merge B926←CLRS and B001←OSQM**, the N-matched
transfer weights, a P3D quiet-month amplitude reference, and an optional un-drop of LRIV, whose
frozen policy excludes it after 2019-05-30 and whose HH templates are a Step-2 product).

### 7.1 What was run

Three configurations on the **S18 test episode** (2024-07-14 … 09-10, the pre-registered
per-episode readout) plus the quiet months **2024-12 and 2025-11** (61 quiet days ≥ the 5 k
quiet-family-day rule at 116 families = 7,076):

| configuration | bank | sites |
|---|---|---|
| `old13` | frozen 2012 | 12 stations, withheld sites out |
| `dense` | frozen ∪ transfer | 18–19 stations, withheld sites out, LRIV-HH reinstated |
| `denseprod` | frozen ∪ transfer | 21–23 stations, every site |

The pooled S14/S21/S22 run was not executed in this session (see § 9).

### 7.2 A fixed MAD is NOT a fair comparison — the matched-FA tiers

At the frozen T_core = 9.25 on 2024-08 the three networks give **6,294 / 4,641 / 3,687** detections.
Adding stations at a fixed MAD REDUCES the count, because the per-family-day network MAD is computed
on the network sum and the added stations raise its noise faster than its signal. That is the same
dilution the Step-3 ablation measures, in the production statistic. The FA-only rule therefore
matters, and it gives:

| network epoch | quiet days | old13 T_core | its quiet-day FA | **matched T_dense** | **matched T_denseprod** |
|---|---|---|---|---|---|
| D3 (2023-09 →) | 61 | 9.25 | **18.28 false events/day** | **8.959** | **8.801** |

### 7.3 Marginal, reverse-marginal and the a3b certification

At those matched thresholds on S18, over 112 families: old13 **8,872** detections, dense **8,469**;
**6,503 matched** one-to-one at ±2 s, **1,966 marginal** (dense-only) and **2,369
reverse-marginal** (old-only).

Each family's marginal and reverse-marginal sets were stacked at the WITHHELD site(s) (B926 + CLRS)
at the family's own predicted S and certified by
`excess_dB = 10 log10(1 + p² N (a/σ)²)`, with (a/σ) calibrated on that family's own hp events at the
same site. 138 of 216 (family, set) legs were certifiable; 78 had too few members.

| | median certified purity | n legs |
|---|---|---|
| marginal (dense-only) | **0.361** | 64 |
| reverse-marginal (old-only) | **0.306** | 74 |

**The detections the dense network adds are no more real than the ones it loses.**

### 7.4 The certified gain vs the hashed prediction

| | value |
|---|---|
| **certified gain, pooled over families with ≥1 new site ≤45 km** | **−5.1 %** |
| bootstrap 2.5 % lower bound | **−8.5 %** |
| **Spearman(certified gain, Δsites ≤45 km)** | **−0.013** (p = 0.95) |
| raw count gain, same 21 scored families | −23.0 % |
| raw count gain, all 78 families with Δsites ≥1 | −1.8 % |
| raw count gain, the 34 families with Δsites = 0 | −8.1 % |
| **hashed Step-0 prediction (D3)** | **+27.9 %** |
| bar (≥ +25 %, bootstrap LB > 0, Spearman > 0.3) | **NOT MET on all three clauses** |

Certified gain by geometry: Δsites = 0 → −3.8 %, Δsites = 1 → −6.1 %, Δsites = 2 → −4.5 %.

**Read the failure signature carefully.** Merlin pre-registered one failure mode — *uniform gain
regardless of geometry, or > 2 × predicted, = the E5 site-noise signature → suspect the tier
derivation*. That is **not** what happened: the result is a uniform small **LOSS**, flat in geometry,
with the added detections no purer than the lost ones. A site-merge error would show up as a gain
concentrated on the co-located pairs; the merge was verified by Step 1 and the co-located stations
are inside single sites in every configuration. The honest reading is the one Step 2 and Step 3
already gave: **the new stations register single LFEs worse than the sites already in the network
(ρ_SYMB/ρ_PGC = 0.78), so adding them to a weighted site sum trades a little signal for a little
more noise, and at matched false-alarm rate that trade is approximately neutral-to-negative.**

The prediction was wrong in the same direction and by about the same factor as the S-SNR prediction
(§ 3.2): the E5 class-plus-distance model gives a new station the *class average* registration for
its distance, and the measured site terms span 6.5 dB. **A distance-and-class model is not a
substitute for measuring a site.** That is the single most transferable lesson in this report.

---

## 8. The coordinator's PQ amendment — V1 and V2, and why the V2 "licence" is an ARTIFACT

`scripts/s_step0_v2.py` → `results/p3d/{stations_p3d_v2.csv, sites_p3d_v2.csv,
premise_family_sites_v2.csv, pred_q2iii_v2.csv, pred_q2iii_recall_v2.csv}` and
`freeze.json → amendment_v2` (hashed, with the PQ **delivery state at hash time recorded as 0
manifest rows and 0 cached station-days for all six V2 stations** — the pre-analysis guarantee).

### 8.1 The premise, with and without PQ

| | before | V1 after | **V2 after (+PQ)** |
|---|---|---|---|
| median sites ≤45 km | 1 | 2 | **4** |
| median sites ≤55 km | 2 | 5 | **8** |
| families gaining ≥1 genuinely new site ≤45 km | — | 82/130 | **100/130** |

Isolated hidden families, new sites ≤45 km: V1 gives 12 (BUCKS+ODUC), 29 (BUCKS+SYMB), 62 (SYMB) and
nothing else; **V2 adds** 10 (+ALBHB), 12 (+ALBHB+JRDRB+LHLYB), 15 (+ALBHB+DAOB),
29 (+ALBHB+DAOB+LHLYB), 62 (+ALBHB+DAOB+JRDRB+LHLYB), 87 (+DAOB). Families 8, 113, 115, 119 still
gain nothing. The ≤60 km geometry rule grows from {PGCq, VGZ, SYMB, BUCKS} to
{+LHLYB, JRDRB, ALBHB, DAOB}. Predicted recall gain **V1 +27.9 %, V2 +78.8 %**.

### 8.2 PQ delivery — two of the four stations are effectively absent

*(figures below are FINAL — the downloader completed all 20,447 queued station-days at 07:28; the
counts quoted in an earlier draft of this addendum were mid-download snapshots. The ratios, and
every conclusion drawn from them, are unchanged.)*

| station | in-epoch SSE days | delivered | nodata | fraction | note |
|---|---|---|---|---|---|
| **DAOB** | 294 | **290** | 4 | **0.99** | S18 complete (59/59) |
| **ALBHB** | 532 | **447** | 85 | **0.84** | but S12 gives 6 days and S13 24 — see § A4 |
| LHLYB | 302 | **230** | 72 | 0.76 | still only **18 of the 59 days of its pre-registered BUILD episode S18** (2024-07 is genuinely `nodata`) |
| **JRDRB** | 598 | **144** | **454** | **0.24** | **effectively does not deliver waveforms**; zero days in any calibration episode |
| OHOH | 505 | 497 | 8 | 0.98 | |
| ODEER | 144 | 144 | 0 | 1.00 | |

So the V2 decision network that can actually be built today is **V1 + ALBHB + DAOB**, and the
pre-registered "LHLYB/DAOB build on S18" clause is not satisfiable for LHLYB.

### 8.3 The S18 readout, and its null-realisation spread

Run on S18 with `--min-avail 0.7` (LHLYB and JRDRB dropped): **1,384 candidates, 956,444 pairs,
24,019 same-family — identical events in every run below**.

| run | δ_old7 | δ_dense | ratio | pre-registered verdict |
|---|---|---|---|---|
| V1, null shift 1800 (default ordering) | 0.544 | 1.557 | 2.86 | NOT LICENSED |
| **V2, +ALBHB+DAOB, null shift 1800** | **0.304** | **1.611** | **5.30** | *"LICENSED"* |
| V1, null shift 1200 | 0.775 | 1.334 | 1.72 | NOT LICENSED |
| V1, null shift 2400 | 0.452 | 1.351 | 2.99 | NOT LICENSED |
| V1, null shift 3600 | 0.336 | 1.803 | 5.36 | *"LICENSED"* |

**δ_old7 mean 0.482, range 0.304–0.775 (±49 %); δ_dense mean 1.531, range 1.334–1.803 (±15 %);
the pre-registered ratio ranges 1.72–5.36 and "licenses" Step 6 in 2 of 5 draws, against a bar of
4 — with no change whatsoever in the data.** The only thing that differs between these runs is which
circular shift each station receives in the null (index × SHIFT seconds), i.e. the null REALISATION.
The denominator is not estimable at this sample size: on S18 the `old7` tail holds 102 / 18 / 3 / 0
same-family pairs at the 1e-3 / 1e-4 / 1e-5 / 1e-6 quantiles, so two or three of the four fit points
exist.

> **The V2 "licence" is a null-realisation artifact, not evidence. Step 6 is NOT licensed.**

What IS stable, because it is measured on identical events with identical stations and a
well-populated tail, is the effect of adding the two delivering PQ stations to the dense sum:
**δ_dense 1.557 → 1.611, i.e. ALBHB + DAOB add +3.5 %** — the same order as SYMB's +9 % in the
calibration run, and far from the ×2.7 the rule needs.

The **well-powered measurement stands**: the 13-episode calibration run (2,834 events, 4.01 M pairs,
all four quantiles populated in both configurations) gives ratio **1.50 → predicted δ_full 0.57 σ →
STEP 6 NOT LICENSED**.

**What the V2 decision still needs** (and what should be re-run when PQ finishes): the same
measurement with (i) PQ's SSE windows complete for LHLYB, (ii) a candidate set of the size the
calibration run had — which for post-2023 stations means pooling S18 with S21 and S22, i.e. spending
the discovery look — and (iii) **δ averaged over ≥5 null realisations**, which this session shows is
mandatory and which the original protocol did not require. That last point is a correction to the
protocol, not to the data.

---

## 9. What remains

| | |
|---|---|
| **Step 4b** (back-projection) | licensed by 4a and NOT RUN: needs full-day CFs over the test episodes + quiet months, the 1,444-cell grid at 6/21/100 false events/day, and the shuffled-moveout / N-shift nulls. Expect clause (c) (≤10 km) reachable and clause (d) (≥3/10 isolated hidden within 5 km) to fail on truth availability |
| **Step 5 pooled** | only the pre-registered S18 per-episode readout was run; S14 + S21 + S22 pooled remains (≈11 further month-runs × 3 configurations, ~1 h each at 3-way concurrency) |
| **Step 3 / V2 decision** | HELD: PQ must finish delivering (LHLYB has 19/59 of its build episode, JRDRB 87/640), and δ must be averaged over ≥5 null realisations |
| **Step 6** | NOT licensed |
| **Step 7 Arm B** | pre-registered and hashed, not run |
| protocol correction | the δ-ratio licence rule needs a null-realisation average and a minimum populated-quantile count; as written it is adjudicable only on a candidate set of the calibration run's size |

---

# ADDENDUM — 2026-09-19b: the PQ component-dip correction

The downloader agent found, against raw IRIS StationXML, that **all four PQ stations declare HHZ at
dip +90 (positive DOWN)** where every CN / UW / PB station in this project declares −90
(positive UP). The frozen 20 Hz cache stores raw traces, so PQ Z was polarity-inverted relative to
the rest of the network in anything that read it. Everything in §§ 0–9 above stands; this addendum
records what was affected, the re-runs, and the before/after.

## A1. Which legs actually touched PQ data — audit

`s_step2_transfer.py` was run with its default `--variant v1`, whose target list comes from
`stations_p3d.csv` (the V1 roster), and `s_step4_cf.py` took its roster from
`s_common.station_table()` — also V1. The audit of the delivered products is therefore:

| leg | PQ stations present? | affected |
|---|---|---|
| Step 2 transfer stacks, ρ table, bank (`step2_legs.csv`, `step2_rho.csv`, `bank_*_keys.csv`) | **no** — BUCKS, CLRS, GOBB, LOPEZ, LRIV, ODUC, OSQM, PTRF, RNWY, SYMB, TURTL only | no |
| Step 4a CF residuals (`step4a_summary.csv`, 16 sites) | **no** | no |
| Step 5 (`trackS/*`, tiers, certification, certified gain) | **no** — bank has no PQ | no |
| Step 3 calibration run (`pairlink3_cal_*`) | **no** — 14 stations, no PQ | no |
| **Step 3 V2 on S18 (`pairlink3_v2_s18_*`)** | **YES — ALBHB + DAOB, raw Z** | **yes** |

**Exactly one delivered leg used PQ data.** The PQ transfer stacks, PQ registration ρ and PQ CF
residuals that the correction asks about did not exist yet; they were built for the first time here,
with the dip honoured.

## A2. The general fix

`s_common` now parses `chan_az_dip` from `results/p3d/station_epochs.csv` for **every** station and
derives both the horizontal azimuths and `ZSIGN[sta] = −1 if dip_Z > 0 else +1`, applied in
`load_day_zne`. Measured result: `{ALBHB, DAOB, JRDRB, LHLYB} → −1`, all 29 other stations `→ +1`,
so every non-PQ product is bit-identical. `S_P3D_NO_DIP=1` reproduces the pre-correction behaviour
through the same code path, which is how the before/after below was measured rather than argued.

**Control (must be a no-op):** SYMB CF 4a re-run under the new code path gives median |residual|
**0.95 s**, `frac_within_1s` 0.5152 — identical to § 5 to every printed digit.

## A3. Before / after

### A3.1 The pair statistic and δ_dense V2 — a no-op, as the algebra requires

The pair statistic correlates two windows **at the same station**; multiplying one component of
*both* windows by −1 leaves the 3-C inner product unchanged. Measured directly on two ALBHB windows:
cc = −0.085813289610 as read, and −0.085813289610 with Z flipped in both — identical to 12 digits.
End-to-end, the full V2 S18 run re-executed with the dip honoured:

| | δ_old7 | δ_dense | ratio |
|---|---|---|---|
| V2 S18, raw Z (§ 8.3) | 0.304 | 1.611 | 5.30 |
| **V2 S18, dip-fixed** | **0.304** | **1.611** | **5.30** |

**The dip error did not change any number in § 8.** It would have mattered for anything that sums
Z coherently *across* stations — the Step-4b back-projection, which was not run — and it matters for
the CF (A3.3).

### A3.2 PQ registration ρ — built for the first time, dip honoured

Transfer stacks on each PQ station's own pre-registered calibration episode, using the **measured**
data extents (§ A4). Identical protocol to § 3, same events, same N, LOO registration:

| station | cal episodes | families ≥N | **ρ_PQ** | ρ_PGC | **ρ_PQ/ρ_PGC** | S-SNR dB | S-SNR PGC | accepted |
|---|---|---|---|---|---|---|---|---|
| **LHLYB** | S18 | 8 | 0.2751 | 0.3523 | **1.086** | 10.29 | 12.65 | 7 |
| **DAOB** | S18 | 33 | 0.3023 | 0.3961 | **0.826** | 11.12 | 13.65 | 25 |
| ALBHB | S12–S13 | 1 | 0.1473 | — | — | 2.12 | — | 0 |
| JRDRB | S21 | — | — | — | — | — | — | not built (no data before 2024-10) |

For reference, § 3.2: SYMB **0.783**, BUCKS 0.296.

**This is the one place the PQ amendment changes the picture.** DAOB (0.83) and especially
LHLYB (**1.09 — the first new station to match PGC**) register single LFEs far better than any
station in the V1 set, and both beat SYMB. Two caveats carried with the number: LHLYB rests on
**8 families** (its data start 2024-08 gives it only 19 of the 59 days of its S18 build episode),
and both are measured on S18 where BUCKS scores 0.296, so the episode is not the explanation.
ALBHB is untestable as a transfer target — its measured extent starts 2022-09, so its calibration
episodes S12–S13 contain almost no hp events and 129 of 130 families are `insufficient`.

**ρ is invariant to the flip** (as it is a within-station statistic): dip-fixed vs raw gives
DAOB 0.302307 / 0.302294 and LHLYB 0.275114 / 0.275092 — agreement to 4 decimals, the residue being
float noise in the rotation path.

### A3.3 PQ CF 4a residuals — the one place the flip is not a no-op

| site | n | **median \|residual\| dip-fixed** | raw Z | Δ | median CF peak dip-fixed / raw |
|---|---|---|---|---|---|
| LHLYB | 2,277 | **1.10 s** | 1.05 s | +0.05 | 0.586 / 0.621 |
| DAOB | 3,764 | **1.20 s** | 1.20 s | 0.00 | 0.594 / 0.613 |
| ALBHB | 5,265 | **1.40 s** | 1.25 s | +0.15 | 0.578 / 0.621 |
| JRDRB | 1,575 | **1.65 s** | 1.65 s | 0.00 | 0.420 / 0.402 |

The CF is close to polarity-blind for timing: the shifts are 0–3 samples of the 20 Hz grid and the
correct-polarity values are marginally *worse*, which is not a result worth interpreting at this
size — it is reported so that the corrected numbers are the ones on record. What does matter is the
level: all four PQ sites fall between SYMB (0.95 s) and the 1.5 s kill bar except JRDRB (1.65 s),
and all four are far better than BUCKS (2.30 s).

Under the V2 geometry rule the ≤60 km site set becomes
{PGCq, VGZ, SYMB, BUCKS, LHLYB, JRDRB, ALBHB, DAOB}; **2 of 8 exceed 1.5 s (BUCKS 2.30, JRDRB 1.65)
→ the CF arm still SURVIVES**, with more margin than it had on V1.

## A4. Measured data extents, and what they do to the episode split

`station_epochs.csv` now carries `first_data_month` / `months_with_data` / `months_probed` from the
downloader's probe. These override the StationXML epochs for the split:

| station | StationXML epoch | **measured first data** | months with data / probed | pre-registered split | **split actually used** |
|---|---|---|---|---|---|
| ALBHB | 2022-01 | **2022-09** | 49 / 57 | calibrate S08–S13 | **S12–S13** (S08–S11 have no data) |
| DAOB | 2024-02 | 2024-02 | 32 / 32 | build S18, test S21/S22 | unchanged |
| LHLYB | 2023-11 | **2024-08** | 26 / 35 | build S18, test S21/S22 | S18 **partial** (19 of 59 days) |
| JRDRB | 2021-04 | **2024-10** | 16 / 66 | calibrate S08–S13 | **build S21, test S22** — no data in any calibration episode |
| OHOH | 2022-04 | 2022-05 | 53 / 54 | — | calibrate S12–S13 |
| ODEER | 2025-06 | 2025-07 | 15 / 15 | — | build S21, test S22 |

This is encoded in `s_common.cal_episodes_for` / `test_episodes_for` and re-hashed as
`freeze.json → amendment_v2b`, which keeps `amendment_v2` verbatim and adds the `z_sign` map, the
measured extents and the corrected StationXML coordinates/elevations. The V2 geometry is unchanged
by the corrected coordinates (they differ from the inventory's by < 10 m): median sites ≤45 km
still 1 → 2 (V1) → **4** (V2), 100/130 families gaining a new site ≤45 km.

## A5. The licence statistic with the protocol fix — and the V2 verdict

`scripts/s_step3_licence.py` implements the correction proposed in § 8.3, in three parts:
**(i)** a configuration's δ counts in a realisation only if ≥ 3 of the 4 tail quantiles carry ≥ 3
same-family pairs; **(ii)** δ is averaged over null realisations and the rule must be met by the
mean *and* by the 2.5 % bootstrap over realisations; **(iii)** at least 3 realisations must be
admissible before the rule may be adjudicated at all — without (iii), (ii) degenerates, because a
bootstrap over one realisation returns that realisation.

| set | admissible realisations | δ_old7 (mean, range) | δ_dense (mean, range) | ratio | predicted δ_full | **verdict** |
|---|---|---|---|---|---|---|
| **V2 (+ALBHB+DAOB), dip-fixed, S18** | **1 of 4** | 0.338 (single) | 1.56 (1.40–1.66) | — | — | **UNADJUDICABLE** |
| V1, S18 | 3 of 4 | 0.590 (0.452–0.775) | 1.414 (1.334–1.557) | 2.40 (boot LB 1.72) | 0.91 σ (LB 0.65) | NOT LICENSED |
| **V1, calibration S01–S13** | **3 of 3**, all 4 quantiles populated in both configurations | **1.228** (1.189–1.253) | **1.870** (1.863–1.876) | **1.52** (boot LB 1.50) | **0.579 σ** (LB 0.569) | **NOT LICENSED** |

The rule fires only above 1.5 σ (ratio ≥ 3.9).

The calibration row is now fully compliant with all three clauses — two further null realisations
(shift 1,200 and 2,400 s) were run to close clause (iii), and they also supply the cleanest possible
contrast with S18:

| candidate set | events / pairs | δ_old7 across null realisations | spread | admissible |
|---|---|---|---|---|
| **calibration S01–S13** | 2,834 / 4.01 M | 1.253, 1.242, 1.189 | **±2.6 %** | **3 of 3** |
| S18 | 1,384 / 0.96 M | 0.775, 0.544, 0.452, 0.336, 0.304 | **±49 %** | 1–3 of 4 |

δ_dense is stable in both (calibration 1.863–1.876, ±0.3 %); it is the **denominator** that decides
whether the rule can be adjudicated, and it is estimable only on a candidate set of the calibration
run's size. That is the quantitative form of the protocol correction:
`results/p3d/step3_licence_V1_CAL_multi.{csv,json}`.

> **V2 verdict: the PQ network cannot adjudicate the licence on S18 — the test is UNADJUDICABLE
> there, not passed and not failed.** In 3 of 4 null realisations the `old7` denominator has fewer
> than 3 populated quantiles. Its S18 tail holds 102 / 18 / 3 / 0 same-family pairs at
> 1e-3 / 1e-4 / 1e-5 / 1e-6 in the most favourable realisation and 69 / 7 / 2 / 0 in the V2
> dip-fixed one — the realisation-to-realisation swing IS the problem, not a property of one draw. The single "STEP6_LICENSED" draw reported in § 8.3 is exactly the
> artifact the fixed rule is designed to reject, and the fixed rule rejects it.
>
> **The operative verdict is therefore unchanged: Step 6 is NOT licensed**, on the only
> well-powered measurement available — the 13-episode calibration run, which now satisfies every
> clause of the corrected rule over 3 null realisations: ratio **1.52** (bootstrap LB 1.50) →
> **0.579 σ** (LB 0.569 σ), against a bar of 1.5 σ. The margin is a factor 2.6, not a near miss.

**What a real V2 decision needs**, now quantified: a candidate set large enough that `old7` keeps
≥ 3 populated quantiles in ≥ 3 null realisations. The calibration run achieves that with 2,834
candidates / 4.01 M pairs; S18 alone has 1,384 / 0.96 M and does not. For the PQ stations that means
pooling S18 + S21 + S22.

**That pooled run is now DATA-READY.** With the download complete, S21 (82 d) and S22 (62 d) are
cached in full for all four PQ stations, and S18 for DAOB (59/59) and ALBHB (54/59); only LHLYB is
short there (18/59). Pooled, S18+S21+S22 carries **2,469 hp events of kept families** against the
1,384 usable on S18 alone.

Two honest qualifications on that, because the arithmetic does not obviously clear the bar:
the 2,469 is counted **before** the all-stations-cached intersection that every pair run applies,
whereas the calibration run's 2,834 is the count **after** it (from 6,316 raw). So the pooled test
set will land materially below 2,469 once intersected — plausibly ~1.8–2.0 k — i.e. **below the
calibration run, not at it.** Pooling S18+S21+S22 therefore buys roughly 1.8 × the events and ~3 ×
the pairs of S18 alone, which may still not be enough to keep `old7` at ≥ 3 populated quantiles in
≥ 3 null realisations. The second qualification is budget: S21 and S22 are two of the four test
episodes, and spending them on a discovery *prediction* is a coordinator decision, not one this
agent may take.

The recommendation is therefore conditional, not a plan: **before spending S21/S22, check the
projected tail occupancy.** `old7` needs ~30 same-family pairs above the 1e-5 quantile to give three
usable fit points; S18 delivers 3, and pairs scale as N², so ~3 × the pairs gives ~9 — still short.
On that projection the pooled test run would probably ALSO be unadjudicable, and the honest
conclusion is that **the V2 discovery prediction cannot be settled on the 2019–2026 north network at
all** — not by S18, and probably not by the whole test set. That is worth knowing before the budget
is spent, and it is the single most useful thing this addendum can say about what to do next.

## A6. Two bookkeeping faults found and fixed while doing this

1. **`s_step2_transfer.py` wrote its bank to a fixed directory** regardless of `--tag`, so the PQ
   run overwrote `bank_transfer_v1/` (251 → 32 templates). **Step 5 is unaffected**: all ten
   bank loads in `logs/s5b_*.log` record `+251 (family, station) templates`, i.e. every month of
   every configuration used the correct bank, which was overwritten only afterwards. The path now
   honours `--tag` (`bank_transfer_v1{tag}`) and both banks were regenerated.
2. The first PQ CF attempt crashed on an empty frame (its `--stations` filter met the V1 roster) and
   the wrapper's `mv` then renamed the **V1** 4a products. The numbers were already in § 5; the
   artifacts (`step4a_residuals.csv`, `step4a_summary.csv`, `step4a_verdict.json`) were regenerated
   from scratch, and `s_step4_cf.py` gained `--variant` / `--out-tag` and an empty-frame guard.

Both are recorded rather than quietly repaired because they are the kind of fault that silently
corrupts a downstream number, and in this case the check that they did not is on the record.
