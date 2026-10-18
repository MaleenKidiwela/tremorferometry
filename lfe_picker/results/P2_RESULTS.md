# P2 — learned characteristic function + explicit Slab2 back-projection: RESULTS

2026-09-17, lfe_picker. Spec: `results/merlin_review_plan_2026-09-17.md` § P2, under the P0-O verdict
**DISCOVERY-ONLY** — *nothing here is a completeness or sensitivity claim*. The 2011 evaluation was
run ONCE, in a single pass covering all three configurations, and is recorded in
`results/referee_ledger.csv`. 2013 was never opened.

## 0. Headline

| item | result |
|---|---|
| **network recall at matched quiet-day FA (2011)** | **0.0363 at 6 false events/day** and **0.0674 at 21/day**, against the frozen template's **0.162 / 0.200** at the same rates → **0.22× / 0.34×** (chance-corrected 0.0242 / 0.0419 → 0.15× / 0.21×) |
| **pre-registered moveout-shuffle null (bar ≤ 0.10)** | **FAIL — 0.90 at 6/day, 0.92 at 21/day.** Diagnosed, not excused: see § 4 |
| supplementary **fixed-cell** moveout null | **PASS, decisively** — at Bostock event times the true table gives 4.75 site-MAD units vs 2.65 shuffled and 0.15 at random times; 6.9× enrichment above the shuffled p99 |
| **fraction of detections > 5 km from every known family** | **0.61** (> 10 km: 0.38) — the discovery-relevant number |
| location error vs the matched family | median **25.6 km** (p90 70 km) — no useful per-event location |
| **calibration within ±0.1** | **NO** — max \|p_real − observed\| = 0.29; ≤ 0.11 over the four bins holding 96 % of detections |
| SSL pretraining gate | **FAILED (AUC 0.757 < 0.85) → pretraining dropped** |

## 1. What was built

**CF.** The P1b-selected PhaseNet-style U-Net (60 s, 20 Hz, 3-C ZNE, peak-normalised, station-agnostic
weights + a station embedding dropped 30 % of the time), warm-started from `pnet_w05` and retrained
for 100,000 steps × batch 256 with two additions:

* **real-waveform transplants** (`scripts/p2_train.py`, `Transplanter`): a loud Bostock event's
  *actual* 3-C windows, taken at its own S time at up to 12 stations on its own day, tapered, scaled
  by an amplitude annealed 1.0 → 0.2, and added into a PNSN-quiet recipient day at the moveout of a
  **random Slab2 cell**. 6,135,239 transplanted windows = **24.0 % of all training samples** (Merlin's
  cap is 30 %); the donor waveform is real, only its position on the interface — and hence its
  inter-station timing — is new, so the CF and the back-projection are exercised away from the 130
  known patches without any position prior being taught.
* **nnPU on the unlabelled ETS-day tremor**, `π = 0.10` (§ 3). This was the single largest lever.

**Back-projection** (`scripts/p2_backproj.py`), no learned location head anywhere: 1,444 Slab2 cells
of 3 km over 48.0–49.0 N / −124.5..−123.0 E, depth fixed on Slab2, Crosson 1-D S times + the a5
station corrections, per-station log-odds MAD-normalised over the station-day, merged into 9 sites by
the mean of MAD units, summed at the predicted S for every cell on a 0.25 s origin-time lattice,
per-cell per-day MAD normalisation, detections = declustered peaks of max_c Z(c, t), threshold set by
**quiet-day false events/day only**.

**Grid honesty note.** Slab2 is a median **7.5 km shallower** than Bostock's own family depths
(p10 −17.3, p90 +6.6 km). A family and its nearest grid cell therefore differ by a median RMS of
**0.50 s** in differential travel time (`results/p2/family_on_grid.csv`). That is a floor on the
location accuracy below, and it is the spec's grid, not a choice made here.

**CF-lag calibration (dev only).** `scripts/p_cflag.py` measures, on 2012, where the trained CF peaks
relative to `t_origin + tt(family's own cell, station)`. After P2 training the lags are −0.20 … +0.30 s
at all 15 stations and the lag against the *label* is ≈ 0.00 s everywhere — including B003/B004, whose
raw S-energy stacks were 2–3 s late under the `predicted` labels (§ P1 § 1). The per-site differential
medians (±0.26 s, `results/p2/site_delay.json`) are added to the travel times. They change nothing
material: the three normalisation/delay variants tested on dev agree to < 10 % at every FA rate.

## 2. The 2011 evaluation (one pass, three configurations)

`bash scripts/p2_run_year.sh results/p2/models/p2nnpu.pt p2nnpu 2011 233,238,239` produced the CF for
all 757 station-days, the 15-station back-projection, its moveout-shuffle null, and the 12-station run
(B926/NLLB/LRIV withheld) used for the independent calibration check. 16 ETS days (18,373 Bostock
events), 8 quiet days.

| quiet-day FA | z | ETS detections | recall | chance (±307 s shift) | chance-corrected | template at the same FA | ratio |
|---|---|---|---|---|---|---|---|
| **6 /day** (= template T_core) | 6.55 | 3,194 | **0.0363** | 0.0124 | 0.0242 | **0.1619** | **0.22×** |
| **21 /day** (= template T_ext) | 5.89 | 7,076 | **0.0674** | 0.0266 | 0.0419 | **0.2002** | **0.34×** |
| 100 /day | 4.92 | 21,921 | 0.1729 | 0.0841 | 0.0970 | — | — |

The gap is not in the threshold, it is in **specificity**. At 6 quiet-day false events/day the frozen
template emits **368 detections per ETS day** (a 61× activity contrast); P2 emits 200/day (a 33×
contrast) but only 21 % of those match a Bostock event. Per Merlin's F2/F8 framing this is exactly
what the P0-O oracle licensed: a discovery engine, not a completeness engine.

### Strata (at 6 /day)

| stratum | n events | recall |
|---|---|---|
| kept families | 13,128 | 0.0341 |
| **hidden families** | 5,245 | **0.0418** |
| families with `picked` S labels | 17,605 | 0.0368 |
| families with only `predicted` S labels | 768 | 0.0247 |
| M ≤ 1.4 | 1,290 | 0.0147 |
| 1.4 < M ≤ 1.8 | 12,104 | 0.0320 |
| 1.8 < M ≤ 2.2 | 4,904 | 0.0514 |
| M > 2.2 | 70 | **0.1286** |
| detections with 8 sites / 9 sites | 1,850 / 1,344 | — |
| families with ≥ 1 event recovered | **91 / 130** | median per-family recall 0.017 |

Two things this settles. (i) **The 30 hidden families are recovered slightly BETTER than the kept
ones** (0.0418 vs 0.0341): the detector has learned an LFE class, not the 100 training families.
(ii) The recall rises monotonically and steeply with magnitude (8.7× from the smallest to the largest
bin) — a chance process would be flat in magnitude, so the detections are real events.

## 3. The nnPU decision (made on dev, before the 2011 look)

Merlin F6 allows ETS-day non-event windows to be *either* excluded *or* used with nnPU and a class
prior. Both were built and compared on 2012:

| dev 2012, 15 stations | `p2core` (unlabelled EXCLUDED) | **`p2nnpu`** (π = 0.10) |
|---|---|---|
| recall at 6 /day | 0.0022 | **0.0428** |
| chance-corrected at 6 /day | 0.0019 | **0.0302** |
| recall at 21 /day | 0.0151 | **0.1022** |
| ETS/quiet detection contrast at z = 6 | 0.9× | **29.7×** |
| moveout-shuffle ratio at 6 /day | 0.56 | 0.87 |

`p2core` failed in a specific, diagnosable way: with no pressure anywhere inside tremor, its CF
saturated at ≈ 0.20 across the whole of an ETS day and rose only to ≈ 0.27 at a catalogued LFE — a
**1.3× contrast with a 2.0–2.35 s FWHM**. Its per-cell per-day normalisation was then dominated by
tremor, so its highest network z values came from *quiet* days (glitches and regional earthquakes),
and at the pre-registered FA it produced fewer detections per ETS day than per quiet day. nnPU fixes
exactly that. π = 0.10 was set a priori, not tuned: catalogued Bostock events occupy ≈ 3 % of an ETS
station-day within ±3 s, and the catalogue is believed ~1/3 complete.

`p2core` is retained as the ablation; it was never run on 2011.

## 4. The moveout-shuffle null — a FAIL, and what it actually measures

Pre-registered bar: shuffled-moveout detections ≤ 10 % of true-moveout detections on ETS days.
Measured on 2011: **0.901 at 6/day, 0.924 at 21/day, 0.925 at 100/day** (and 0.95–1.06 when the null
is given its own quiet-day FA calibration). On dev the same numbers are 0.87–0.95. **FAIL.**

The shuffle is real — permuting the site → travel-time assignment changes the differential moveout by
an RMS of 10.2 s against a true differential spread of 6.2 s. So the failure is not a broken null.

What it means is visible in the threshold dependence (dev, `p2core`): the ratio falls 0.98 → 0.13 as
z goes 3.0 → 6.5. The detection statistic is **max over 1,444 cells**. With that many cells a wrong
travel-time table still finds *some* cell whose alignment happens to match a tremor burst, so a
count-based comparison of two maxima is a weak test of the geometry — it mostly measures how often
tremor is loud, which both tables see.

The supplementary **fixed-cell** null (`scripts/p2_fixedcell_null.py`) removes the maximisation: at
each of the 18,373 Bostock events on 2011 ETS days, read the site-summed CF at the moveout of that
event's **own** family cell, under (a) the true table, (b) the shuffled table, (c) the true table at
random times on the same day.

| | median | 75 % | max |
|---|---|---|---|
| **true moveout at the event** | **4.75** | 8.30 | 35.1 |
| shuffled moveout at the event | 2.65 | 5.52 | 23.5 |
| true moveout at random times | 0.15 | 2.81 | 31.4 |

Fraction of true-moveout values above the **shuffled p99**: **0.069 (6.9× chance)**; above the
random-time p99: **0.136 (13.6× chance)**. The physics carries real information; the pre-registered
null statistic is the wrong instrument for a max-over-cells detector, and a fair version of it should
be written before P3.

## 5. Location, and the centroid-collapse diagnostic

| at 6 /day | value |
|---|---|
| location error of matched detections vs their Bostock family (median) | **25.6 km** (p90 70.3) |
| distance from every detection's best cell to the nearest of the 130 families (median) | 6.9 km |
| **fraction > 5 km from every known family** | **0.609** |
| fraction > 10 km from every known family | 0.380 |

There is **no centroid collapse** — the opposite. Detections do not pile onto the 130 catalogued
positions (as the learned-location-head pathology of E3 would have produced); they spread over the
grid, with a concentration along the 48.4–48.6 N band where the family cluster sits
(`figures/p_p2_map_p2nnpu_2011_{233,238,239}.png`). But a 25.6 km median location error means the
per-event cell is **not a usable location**; the 0.61 "away from every family" figure must therefore
be read as *"the back-projection's argmax is far from the known patches"*, **not** as *"61 % of these
are new sources"*. With this location accuracy the honest use of P2 output is as a time-domain
proposal generator whose locations are refined later by stacking, exactly the role Merlin's F10
assigned it.

Contributors to the 25.6 km: the 0.50 s family-to-cell travel-time mismatch from the Slab2 depth
offset (§ 1), the 2 s FWHM of the CF peak, and the near-degenerate moveout geometry that § 4 exposes.

## 6. Calibration against withheld-station realness

Isotonic `p_real` fitted on the **12-station dev run** (2012, B926/NLLB/LRIV excluded from the
back-projection) and applied to the 12-station 2011 run, so the three withheld stations are a genuinely
independent test. "Withheld" below = the fraction of a detection's withheld stations whose S-window
energy excess (E(S−1..S+5) / E(20 s earlier), that station-day's MAD units, at the predicted S of the
detection's own cell) exceeds that station-day's own 200-sample random-time p95. Under "all false" that
rate is 0.05 by construction.

| p_real bin | n | mean p_real | Bostock match rate | withheld-station excess rate | \|dev\| |
|---|---|---|---|---|---|
| 0.1–0.2 | 921 | 0.196 | 0.148 | 0.121 | 0.048 |
| 0.2–0.3 | 937 | 0.277 | 0.179 | 0.121 | 0.098 |
| 0.3–0.4 | 558 | 0.348 | 0.235 | 0.163 | 0.113 |
| 0.4–0.5 | 268 | 0.437 | 0.306 | 0.189 | 0.131 |
| 0.5–0.6 | 82 | 0.581 | 0.329 | 0.207 | 0.252 |
| 0.7–0.8 | 32 | 0.726 | 0.438 | 0.302 | 0.288 |

**Verdict: NOT within ±0.1** (max 0.288). Two qualifications, both stated rather than used as an
escape: (i) the four lowest bins, holding 2,684 of 2,798 detections, are within 0.11; (ii) the Bostock
match rate is a **lower bound** on realness, so a calibrated `p_real` *should* sit above it — the
withheld-station column is the independent estimate, and it rises monotonically 0.121 → 0.302 against
a 0.05 chance floor, i.e. a **2.4× to 6.0× enrichment** that tracks `p_real` in the right direction.
What fails is the absolute level at high `p_real`, which was fitted on a year with ~20 % more events
per ETS day than 2011 (2012: 1,388 Bostock events per ETS day; 2011: 1,148).

## 7. Pretraining gate — FAILED, pretraining dropped

`scripts/p_ssl_probe.py`. Pretext = cross-station envelope reconstruction, trained on **2014–2016
only**, with the target station ≥ 15 km from every input station (so the PGC and LRIV intra-site
clusters, Merlin's F4, cannot be what is learned). Linear probe = logistic regression on the frozen
encoder's pooled features, 2012 Bostock event windows vs in-tremor non-event windows of the **same**
station-days, grouped 4-fold by station-day.

**AUC 0.757** (folds 0.766 / 0.749 / 0.765 / 0.747) against the 0.85 bar → **pretraining is dropped**,
and nothing in P2 uses it. This is F4 confirmed rather than assumed.

## 8. Files

Models `results/p2/models/{p2core,p2nnpu,ssl_pretext}.pt`; 2011 detections
`results/p2/det_p2nnpu_2011{,_shuf,_s12}.parquet` and the scored table
`results/p2/det_scored_p2nnpu_2011_s12.parquet`; per-cell probability maps for the 2011 SSE days
`results/p2/maps_p2nnpu_2011_{233,238,239}.npz`; scores `results/p2/score_p2nnpu_2011*.json`,
`calib_p2nnpu_2011_s12.csv`, `fixedcell_null_2011.csv`, `ssl_probe.json`; figures
`figures/p_p2_map_p2nnpu_2011_*.png`, `figures/p_p2_diag_p2nnpu_2011_s12.png`.

## 9. What P2 is and is not

* It **is** a source-agnostic proposal generator that, at the template's own quiet-day false-alarm
  rate, finds a fifth of what the templates find while knowing nothing about the 130 families —
  including the 30 it was never shown, at the same rate as the ones it was.
* It is **not** a completeness engine (P0-O settled that before a GPU was touched), and it is **not**
  a locator: 25.6 km median error.
* The pre-registered moveout-shuffle null **fails**, and the honest reading is that the statistic is
  wrong for a max-over-cells detector, not that the geometry is uninformative — the fixed-cell version
  of the same null passes at 6.9× chance. **A corrected null must be pre-registered before P3**, and
  no discovery claim should be made until it is.

---

# ADDENDUM (2026-09-19, P3 Step 0b/0c) — **P2 back-projection CLOSED as locator/discovery engine**

Merlin's P1/P2 review (`results/merlin_review_p1p2_2026-09-19.md`) reversed §4 and §5 of this
document, on numbers recomputed from P2's own stored outputs. Both recomputations are reproduced
here (`scripts/q_step0_bc.py`, outputs in `results/p3/`). Nothing was re-run, nothing re-tuned, no
referee budget consumed, `referee_ledger.csv` unchanged.

## 1. The corrected event-level null — the shuffled detector keeps 60 % of the recall

§4 compared *counts* of detections and called the pre-registered null "the wrong instrument". The
statistic that is right for a back-projection detector is the shuffled run's **own chance-corrected
Bostock recall**, against the true run's, at the *same* thresholds
(`results/p3/p2_corrected_event_null.csv`):

| year | quiet-day FA | true recall (cc) | moveout-shuffled recall (cc) | **ratio** | (count ratio, §4) |
|---|---|---|---|---|---|
| 2011 | 6 /day | 0.0242 | **0.0145** | **0.600** | 0.901 |
| 2011 | 21 /day | 0.0419 | 0.0258 | 0.616 | 0.924 |
| 2011 | 100 /day | 0.0970 | 0.0686 | 0.707 | 0.925 |
| 2012 (dev) | 6 /day | 0.0302 | 0.0186 | 0.616 | 0.865 |
| 2012 (dev) | 21 /day | 0.0633 | 0.0438 | 0.692 | 0.910 |

A travel-time table that is *wrong* recovers **60–70 % of the events the right one recovers**. The
pre-registered bar was ≤ 0.10. The null was not the wrong instrument: it measured the detector as
used, and it was right. Of the true run's 0.0242 chance-corrected recall at 6/day, at most ~40 % is
attributable to the geometry, and 0.0242 is itself 0.15 × the frozen template's 0.1619. The
fixed-cell null in §4 is not withdrawn — the physics does carry information when the maximisation
over 1,444 cells is removed — but it does not license the max-over-cells detector.

## 2. The "61 % > 5 km from every known family" equals the uniform-over-grid expectation

`results/p3/p2_uniform_grid_expectation.json`. Of the 1,444 Slab2 cells, **64.5 %** are more than
5 km from every one of the 130 families and 42.2 % are more than 10 km (median 7.7 km). The observed
detection figures are **60.9 %** (2011) and **59.1 %** (2012) > 5 km, and 38.0 % / 37.5 % > 10 km, on
argmax cells spread over 1,229 (2011) and 1,310 (2012) of the 1,444 cells. The observed numbers are
*below* the uniform expectation. **The 61 % carries no information about new sources** and is
withdrawn as a discovery-relevant statistic; §5's 25.6 km median location error is the operative
number.

## 3. The loudest ETS days — the per-cell per-day MAD normalisation collapses

`results/p3/loud_day_collapse.csv`. Detections per day at the 6/day operating point, against the
Bostock event count of the same day and the frozen template's T_core catalog:

| year | doy | Bostock events | **P2 detections** | template (T_core) |
|---|---|---|---|---|
| 2011 | 239 | 2,262 | **10** | 402 |
| 2011 | 240 | 1,225 | **10** | 201 |
| 2011 | 241 | 1,447 | **5** | 195 |
| 2011 | 242 | 1,169 | **10** | 208 |
| 2011 | all other ETS days (mean) | 1,023 | 263 | 241 |
| 2012 | 253 | 1,117 | **1** | 183 |
| 2012 | 255 | 1,843 | **17** | 502 |
| 2012 | 257 | 1,326 | 13 | 505 |
| 2012 | 260 | 1,801 | 17 | 875 |
| 2012 | 267 | 526 | 549 | 165 |
| 2012 | all other ETS days (mean) | 1,428 | 201 | 522 |

On the four loudest 2011 days — which hold 6,103 of the year's 18,373 ETS-day events — P2 emits
**5–10 detections per day**, a 26–50 × collapse against its own ETS-day average, while the template
detector is essentially flat (135–178 per 1,000 events against 237 on the other days). The same
signature appears on 2012 doy 253, 255, 257 and 260. The cause is the two-sided per-cell per-day MAD:
on a day where tremor is on for most of the day the MAD is set *by the signal*, so the threshold
tracks the activity it is supposed to detect. `lowhalf` normalisation existed
(`results/p2/det_probe2012lh.parquet`, `det_dev_lowhalf_yes.parquet`) and was not used in the
production runs. The reported recall of 0.036 at 6/day is therefore **depressed**, and any future use
of this CF must normalise against a quiet reference or the lower half of the distribution.

## 4. What stands, what falls

* **Falls:** P2 as a locator; P2 as an event-level discovery engine; the 61 % figure; the count-based
  reading of the moveout null.
* **Stands:** the learned CF detects LFEs *as a class* (hidden families recovered at the same rate as
  trained ones, 0.0418 vs 0.0341); recall rises 8.7 × with magnitude, so the detections are real
  events; the fixed-cell null (physics carries information when the cell is fixed); SSL gate FAILED.
* **Surviving roles for the CF (each with a kill test, Merlin P3 Step 2):** a trigger-ablation arm
  for a repeater search, and a candidate per-detection realness feature for Track E.
