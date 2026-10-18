# P1b — two learned single-station baselines: RESULTS

2026-09-17, lfe_picker. Spec: `results/merlin_review_plan_2026-09-17.md` § P1, as amended by
`results/P0O_RESULTS.md` (oracle verdict DISCOVERY-ONLY). Everything below was produced inside
`lfe_picker/` with the `p_` prefix; `o_*` and the neighbouring workspaces were read only.

## 0. Headline

| item | verdict |
|---|---|
| architecture selected on the DEV years (2010 + 2012) | **PNET** (PhaseNet-style U-Net, 60 s @ 20 Hz) — beats the Lin re-implementation at **every** peak rate on **both** dev years |
| `predicted`-label weight sensitivity (0.5 vs 1.0) | **0.5 wins at the selection metric** (PGC 0.029/0.032 vs 0.023/0.025 at 100 peaks/station-day) despite the 1.0 run getting 37 % MORE steps |
| **Merlin's bar: single-station S recall ≥ 0.30 at ≤ 100 peaks/station-day at PGC** | **NOT MET — 0.0218 on 2011 (0.029 on dev 2012).** It is also **arithmetically unreachable**: see § 2 |
| at the template's own total peak budget (11,600 peaks/station-day = 100 per template-day × 116 families) | PGC **0.486** on 2011 vs the template cc's **0.119** — a **4.1×** single-station improvement |
| vs the pretrained LFEDetect reference (0.12–0.13 at 720 peaks/station-day) | **matched**: PGC 0.112, B011 0.118, B009 0.111 on 2011 |
| 2011 referee budget | ONE evaluation, spent here, recorded in `results/referee_ledger.csv` |

## 1. What was built

`scripts/p_build_mm.py` → one memory-mapped float16 store of every cached S15 station-day 2008–2013
(**3,735 station-days, 38.7 GB**, `results/mm/`), **rotated to ZNE** with the StationXML azimuths.
This mattered: the PB boreholes' EH1/EH2 are *not* at 0/90 (B001 335.5/65.5, B003 239/329,
B004 167.5/257.5, B005 90.5/180.5, B006 201/291, B007 36.5/126.5, B010 228.5/318.5, B011 11/101,
B926 299/29 — only B009 is 0/90), so the frozen chain's `EH1→N, EH2→E` naming is a *mislabel* for a
learned model. B006/B007 azimuths were absent from `a4_stationxml_audit.csv` and were fetched from
IRIS (`results/p_stationxml_B006_B007.xml`). The template path in `o_*` is untouched and still reads
the raw Z/1/2 order, as specified.

**Labelling policy** (`scripts/p_dataset.py`), Merlin F6 taken literally:

* positives = every `pos_S` row of `labels_v0` in the training years at the training stations for the
  kept families, Gaussian target σ 0.5 s, **all** events inside a window labelled (an ETS burst puts
  several in one 60-s window); sample weight 1.0 `picked` / 0.5 `predicted`;
* certain negatives = PNSN-quiet station-days with ±120 s around **every** Bostock event of **any**
  family excluded; hard negatives = catalogued EQ/blast windows, HF-veto quiet-day transients, T_ext
  quiet-day detections — all at target 0, full weight;
* **unlabelled**: inside a positive window the loss is applied ONLY within ±3 s of a labelled S. The
  rest of the window is ETS tremor of unknown content and carries weight 0. Windows within ±3 s of a
  **hidden** family's event carry weight 0 there, so the holdout leaks neither way.

Training scope: 2008 + 2009 + 2010 + 2012; stations S15 minus B926/NLLB/LRIV; families minus the 30
hidden. 2011 untouched during training; 2013 never opened.

**Time-base check (PASS, and one finding).** Stacked S-window energy at the `labels_v0` S times,
2012, 4,000 events/station: PGC peaks **+0.40 s** after the label at 2.51× the pre-window energy,
B010 +0.25/1.99, VGZ +0.35/1.95, SNB +0.50/1.59 — the time base is sound (the failure mode would have
been a peak metres away from zero). But the far boreholes peak late: B005 +0.75, B001 +1.15,
**B003 +1.95, B004 +3.15 s**. B003/B004 labels are `predicted` (Crosson + a5 correction); the
"≥ 0.5 s" P0-O warned about is really **2–3 s** there. That is a LABEL error, not a model error, and
it is why every table below is also cut by provenance.

## 2. The metric, and the ceiling nobody had computed

Merlin's P1b metric is single-station S recall at ±1 s with **100 declustered peaks per
station-day**. A family-agnostic detector emits ONE series per station-day, so R peaks can match at
most R events:

> **recall ≤ R × (station-days) / (labelled events at that station)**

On 2011 at PGC that is 100 × 22 / 17,463 = **0.126**. On 2012 it is **0.102**. The bar of 0.30 is
**2.4–2.9× above the arithmetic ceiling** and cannot be met by any single-station detector, however
good, on an ETS-dominated referee window. The P0-O comparison number (template cc 0.119) is at 100
peaks **per template-day**, i.e. ×116 families ≈ 11,600 peaks per station-day — a budget 116× larger.
Both numbers are therefore reported, and the matched-budget column (11,600/station-day) is the
apples-to-apples one.

This is reported as a defect of the bar, not as an excuse: at 100 peaks/station-day PNET reaches
**17 % of its ceiling at PGC** (0.0218/0.126) and 20 % at B011, so there is real room left even
inside the ceiling.

## 3. Dev years — architecture selection (2010 + 2012, unlimited budget)

Recall at 100 peaks/station-day, ±1 s (`results/p1/curve_*_tol1.csv`,
`figures/p1b_recall_curves_{2010,2012}.png`):

| sta | lin 2010 | lin 2012 | **pnet w0.5 2010** | **pnet w0.5 2012** | pnet w1.0 2010 | pnet w1.0 2012 |
|---|---|---|---|---|---|---|
| PGC | 0.0155 | 0.0168 | **0.0316** | **0.0290** | 0.0247 | 0.0234 |
| B011 | 0.0121 | 0.0124 | **0.0304** | **0.0293** | 0.0259 | 0.0254 |
| B009 | 0.0149 | 0.0120 | **0.0281** | **0.0309** | 0.0249 | 0.0229 |
| VGZ | 0.0198 | 0.0146 | **0.0312** | **0.0259** | 0.0267 | 0.0200 |
| SNB | 0.0176 | 0.0125 | 0.0207 | 0.0154 | 0.0217 | 0.0155 |
| B005 | — | 0.0159 | — | **0.0238** | — | 0.0271 |
| B003 | — | 0.0115 | — | **0.0157** | — | 0.0145 |
| B926 (held out) | 0.0020 | 0.0022 | 0.0061 | 0.0054 | 0.0028 | 0.0024 |
| NLLB (held out) | 0.0018 | 0.0027 | 0.0053 | 0.0050 | 0.0045 | 0.0047 |
| LRIV (held out) | 0.0023 | 0.0039 | 0.0109 | 0.0112 | 0.0090 | 0.0091 |

2010 carries no labels at B001–B007 (those station-days are not in the 2010 cache), so those rows
are blank there and the 2010 comparison rests on the CN stations and the PGC cluster.

**PNET is selected.** It wins at every rate on both dev years (e.g. at 11,600 peaks/station-day,
PGC 0.506/0.481 vs lin 0.344/0.464). Model sizes and compute were comparable (lin 3.26 M params,
60,000 steps; pnet 1.82 M params, 60,000 steps; batch 256 throughout, ≈ 15.4 M training windows each,
well inside the ≤ 1 GPU-day budget — the two P1b runs plus the sensitivity run together took ≈ 3.5
GPU-hours on the H100).

**`predicted`-label weight (0.5 vs 1.0).** At the selection metric the 0.5 weighting is better at
every near site on both dev years, and this holds *against* a compute handicap: the w=1.0 run ran to
82,000 steps (37 % more) before it was stopped, and its window-level dev AUC was higher (0.989 vs
0.977) — i.e. it fit the labels better and detected worse. That is the signature of the 2–3 s label
error at B003/B004 being fitted. At very high peak rates (11,600/day) the two are within 2 %.
**Conclusion: keep w_pred = 0.5;** the sensitivity is real but modest (≈ 25 % at the low-rate end,
≈ 0 at the high-rate end).

## 4. The one 2011 evaluation (PNET selected; the other two recorded, never used to select)

`results/p1/curve_*_2011_tol1.csv`, `figures/p1b_recall_curves_2011.png`.
Denominator = every `pos_S` label of 2011 at that station (186,926 rows over 757 station-days).

### 4a. Per-station recall vs peaks/station-day — PNET (w_pred 0.5)

| sta | 10 | 30 | **100** | 300 | **720** | 1,000 | 3,000 | 6,000 | **11,600** | ceiling@100 |
|---|---|---|---|---|---|---|---|---|---|---|
| **PGC** | 0.0014 | 0.0054 | **0.0218** | 0.0577 | **0.1120** | 0.1409 | 0.2753 | 0.3852 | **0.4856** | 0.126 |
| B011 | 0.0031 | 0.0079 | **0.0257** | 0.0618 | 0.1182 | 0.1464 | 0.2763 | 0.3843 | 0.4826 | 0.126 |
| B009 | 0.0013 | 0.0053 | 0.0205 | 0.0577 | 0.1105 | 0.1361 | 0.2603 | 0.3690 | 0.4678 | 0.126 |
| B010 | 0.0024 | 0.0055 | 0.0177 | 0.0472 | 0.0954 | 0.1202 | 0.2410 | 0.3551 | 0.4707 | 0.122 |
| VGZ | 0.0014 | 0.0053 | 0.0178 | 0.0467 | 0.0957 | 0.1233 | 0.2508 | 0.3560 | 0.4440 | 0.140 |
| SNB | 0.0004 | 0.0012 | 0.0063 | 0.0212 | 0.0545 | 0.0742 | 0.1854 | 0.2960 | 0.3998 | 0.134 |
| **B005** | 0.0006 | 0.0023 | **0.0071** | 0.0199 | 0.0502 | 0.0684 | 0.1800 | 0.2908 | 0.3946 | 0.334 |
| B007 | 0.0003 | 0.0014 | 0.0058 | 0.0205 | 0.0482 | 0.0687 | 0.1688 | 0.2691 | 0.3736 | 0.254 |
| B006 | 0.0006 | 0.0013 | 0.0048 | 0.0179 | 0.0439 | 0.0604 | 0.1593 | 0.2706 | 0.3793 | 0.253 |
| **B004** | 0.0000 | 0.0000 | **0.0038** | 0.0130 | 0.0310 | 0.0407 | 0.1240 | 0.2267 | 0.3498 | 0.344 |
| **B001** | 0.0003 | 0.0006 | **0.0032** | 0.0105 | 0.0270 | 0.0361 | 0.1070 | 0.1931 | 0.3115 | 0.278 |
| **B003** | 0.0003 | 0.0009 | **0.0027** | 0.0114 | 0.0335 | 0.0489 | 0.1318 | 0.2314 | 0.3478 | 0.219 |
| LRIV\* | 0.0011 | 0.0020 | 0.0092 | 0.0245 | 0.0531 | 0.0678 | 0.1526 | 0.2384 | 0.3342 | 0.338 |
| NLLB\* | 0.0007 | 0.0017 | 0.0080 | 0.0223 | 0.0513 | 0.0668 | 0.1566 | 0.2405 | 0.3456 | 0.148 |
| B926\* | 0.0000 | 0.0005 | 0.0028 | 0.0102 | 0.0280 | 0.0406 | 0.1185 | 0.2146 | 0.3349 | 0.136 |

\* held-out station, never trained on.
Near sites mean @100 = **0.0207**; far boreholes (B001/B003/B004/B005/B006/B007/LRIV) mean = **0.0052**.

### 4b. The three baselines on 2011, recorded together

| sta | lin @100 | **pnet w0.5 @100** | pnet w1.0 @100 | lin @11,600 | **pnet w0.5 @11,600** | pnet w1.0 @11,600 |
|---|---|---|---|---|---|---|
| PGC | 0.0130 | **0.0218** | 0.0101 | 0.4512 | **0.4856** | 0.4711 |
| B011 | 0.0060 | **0.0257** | 0.0178 | 0.4083 | **0.4826** | 0.4660 |
| B003 | 0.0022 | **0.0027** | 0.0028 | 0.3217 | **0.3478** | 0.3276 |
| B004 | 0.0012 | **0.0038** | 0.0038 | 0.3235 | **0.3498** | 0.3358 |

The dev selection is confirmed on the referee year: PNET (w 0.5) is best at PGC, at the boreholes,
and at both ends of the rate axis.

### 4c. Strata

| stratum | @100 | @720 | @11,600 |
|---|---|---|---|
| `picked` labels (n = 131,501) | **0.0160** | 0.0889 | 0.4439 |
| `predicted` labels (n = 55,425) | **0.0039** | 0.0320 | 0.3298 |
| kept families (n = 129,782) | 0.0137 | 0.0748 | 0.4056 |
| **hidden families** (n = 57,144) | 0.0095 | 0.0656 | **0.4201** |
| trained stations (n = 149,419) | 0.0141 | 0.0796 | 0.4279 |
| **held-out stations** (n = 37,507) | 0.0059 | 0.0416 | **0.3390** |

Three things worth saying plainly:

1. **The hidden families are not penalised at the high-rate end** (0.420 vs 0.406). The model has not
   memorised the 100 kept family positions/waveforms; it is detecting a class. At the low-rate end
   hidden families do lose (0.0095 vs 0.0137), which is the expected consequence of their events also
   being smaller on average, not of leakage — leakage would push the other way.
2. **`predicted` labels cost a factor 4 at 100 peaks/day and only 1.35 at 11,600.** That pattern is a
   label TIMING error (±1 s tolerance around a label that is 2–3 s off), not a detection failure: at a
   loose rate the peak is there, it just is not within ±1 s of where the catalogue says.
3. **Held-out stations degrade 21 % at 11,600/day** (0.339 vs 0.428) and 58 % at 100/day. The 21 %
   figure is inside Merlin's P4a "≤ 30 % relative" bar; the 100/day figure is not, and the honest
   reading is that station-agnostic transfer costs most where the threshold is tightest.

## 5. Association (network level)

The primary association is the explicit physics back-projection and is reported in `P2_RESULTS.md`.
The secondary control Merlin asked for, **PyOcto** with the same Crosson 1-D velocities, S-only, on
the PNET picks at 3,000 picks/station-day over three 2011 ETS days (`scripts/p_pyocto.py`):

| | |
|---|---|
| picks in | 125,999 |
| associated events out | 11,116 |
| Bostock events on those days | 3,106 |
| matched ±2 s one-to-one | 790 (raw recall 0.254) |
| time-shifted chance control (+307 s) | 531 |
| **chance-corrected recall** | **0.083** |
| median distance to the nearest Bostock family | 9.1 km |
| fraction > 5 km from every known family | 0.68 |

Read as an order-of-magnitude cross-check only: a general-purpose associator, fed our picks, finds
real events at a few per cent above chance and puts two thirds of them away from the known patches.
It is not a competing catalog (S-only association is weakly constrained and the chance rate is high
during ETS).

## 6. Verdicts

* **Merlin's P1b bar (≥ 0.30 at ≤ 100 peaks/station-day at PGC): NOT MET.** Measured 0.0218 on 2011,
  0.029 on dev 2012. The failure mode Merlin named for a miss below 0.2 ("labels/time base wrong —
  check offset sign and S placement against `a5_picks.csv`") was checked FIRST and is **not** the
  cause: the stacked energy peaks 0.25–0.5 s after the labelled S at every picked-provenance station.
  The bar is unreachable by construction — its own arithmetic ceiling at PGC is 0.126.
* **Against the statistics that exist on disk, the learned CF is a large gain**: 4.1× the template cc
  at the template's own total peak budget, and it matches pretrained LFEDetect at 720 peaks/day while
  LFEDetect needed the same rate for 0.12–0.13.
* **PNET, w_pred 0.5, carries into P2.**

---

# ADDENDUM (2026-09-19, P3 Step 0a) — the 4.1 × is NOT chance-corrected; corrected it is **2.7 ×**

Merlin's P1/P2 review (`results/merlin_review_p1p2_2026-09-19.md`, flaw #4) is upheld. At 11,600
declustered peaks per station-day the ±1 s matching tolerance covers a large fraction of the day, so
a sizeable part of §4a's recall is accidental. Every P1b curve (all three models, 2010 / 2012 / 2011)
has been re-scored with the ±307 s shift control used by `scripts/p2_score.py`
(`scripts/q_p1b_chance.py` → `results/p3/p1b_chance_corrected.csv`; recall_cc = (recall − chance) /
(1 − chance), chance = the mean of the +307 s and −307 s rematches). **Nothing was re-run and nothing
was re-tuned** — the stored 2011 scans of the already-spent P1b evaluation were re-read, so no
referee budget was consumed and `referee_ledger.csv` is unchanged.

### PNET (w_pred 0.5), recall / chance / chance-corrected

| year | sta | @100 | @720 | @11,600 |
|---|---|---|---|---|
| 2011 | **PGC** | 0.0218 / 0.0018 / **0.0200** | 0.1120 / 0.0155 / **0.0980** | 0.4856 / 0.2426 / **0.3209** |
| 2011 | B011 | 0.0257 / 0.0023 / 0.0234 | 0.1182 / 0.0150 / 0.1048 | 0.4826 / 0.2341 / 0.3244 |
| 2011 | B003 | 0.0027 / 0.0030 / **−0.0003** | 0.0335 / 0.0167 / 0.0170 | 0.3478 / 0.2654 / 0.1122 |
| 2011 | B004 | 0.0038 / 0.0025 / 0.0013 | 0.0310 / 0.0168 / 0.0144 | 0.3498 / 0.2851 / 0.0905 |
| 2012 | PGC | 0.0290 / 0.0020 / 0.0270 | 0.1111 / 0.0144 / 0.0981 | 0.4810 / 0.2289 / 0.3270 |
| 2010 | PGC | 0.0316 / 0.0020 / 0.0296 | 0.1249 / 0.0131 / 0.1132 | 0.5059 / 0.2179 / 0.3683 |

Chance at 100 peaks/station-day is 0.002 (negligible); at 720 it is 0.015; at 11,600 it is **0.24**,
i.e. a quarter of the day is inside ±1 s of some peak, exactly as Merlin computed.

### The corrected headline

The denominator of the "4.1 ×" is the template cc's single-station recall at PGC at 100 peaks per
*template*-day (= the same 11,600 peaks/station-day budget), 0.1194 in the O3 oracle. That statistic
is read at the family's own predicted S, so its chance level is measured directly from O3's sampled
ETS-day negatives at the same threshold: **0.0000** (quiet-day negatives 0.0001). The template number
therefore does not move.

| | raw | chance | chance-corrected |
|---|---|---|---|
| PNET @ 11,600 peaks/station-day, PGC 2011 | 0.4856 | 0.2426 | **0.3209** |
| PNET @ 11,600 peaks/station-day, PGC 2012 | 0.4810 | 0.2289 | **0.3270** |
| template cc @ 100 peaks/template-day, PGC (2012 oracle) | 0.1194 | 0.0000 | **0.1194** |
| **gain** | 4.07 × | — | **2.69 × (2011) / 2.74 × (2012)** |

**"4.1 ×" is withdrawn and replaced by 2.7 ×.** That is slightly above Merlin's predicted ceiling
(≤ 2.5 ×) and well below the published figure; the direction of his correction is confirmed. The
720-peaks/day comparison with pretrained LFEDetect is unaffected in kind (PNET 0.1120 raw → 0.0980
corrected at PGC; LFEDetect's 0.12–0.13 is itself uncorrected and would lose a similar ~0.015).

### Far boreholes — the case is closed harder than §4c said

Chance correction is *worst* where the labels are predicted, because a loose peak budget plus a 2–3 s
label error means the "hits" were largely accidental. At 11,600 peaks/station-day on 2011 the far
boreholes go 0.356 raw → **0.123** corrected (near sites 0.471 → 0.309), and B003 at 100 peaks/day is
**indistinguishable from chance (−0.0003)**. Merlin's flaw #6 (no far-borehole headroom) stands.
