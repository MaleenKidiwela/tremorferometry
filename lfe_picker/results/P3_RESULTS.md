# P3 — closure bookkeeping and the template-free repeater search: RESULTS (Steps 0, 1a, and the CF completeness kill test)

2026-09-19, `lfe_picker`. Spec: `results/merlin_review_p1p2_2026-09-19.md` § "P3 prescription",
Steps 0–2. Everything below was produced inside `lfe_picker/` with the `q_` prefix; `o_*`, `p_*`,
`g_*` and the neighbouring workspaces were read only (the one change needed to a frozen script — a
station override — was made as an in-process monkey-patch in `scripts/q_run_month.py`, documented
there; `o_run_month.py`, `o_common.py` and the frozen operating-point json are byte-identical).

**Referee discipline.** 2012 (dev, unlimited) carries every measurement in §2–§4. 2013 was not
opened. 2011 was **not re-run**: Step 0 re-scores the *stored products* of the P1b and P2 evaluations
that were already spent and recorded on 2026-09-17, which selects and tunes nothing;
`results/referee_ledger.csv` is unchanged and no new referee budget was consumed.

---

## 0. Headline

| step | result |
|---|---|
| **0a** chance-corrected P1b | the **"4.1 ×" is withdrawn → 2.69 × (2011) / 2.74 × (2012)**. Chance at 11,600 peaks/station-day is 0.24. Addendum in `P1_RESULTS.md` |
| **0b** corrected event-level null | the moveout-shuffled detector keeps **0.60** of the true detector's chance-corrected recall (bar ≤ 0.10); "61 % > 5 km from every family" = the uniform-grid expectation **64.5 %**. Addendum in `P2_RESULTS.md`, headed "P2 back-projection CLOSED as locator/discovery engine" |
| **0c** loud-day collapse | on 2011 doy 239–242 P2 emits **5–10 detections/day** against 1,169–2,262 Bostock events and 195–402 template detections — a 26–50 × collapse; the template is flat |
| **0d** relative-clause scorer | `scripts/q_rediscover.py` + the truth-template reference table (§1), measured on 2012 with SNB + B926 withheld |
| **1a** pair-linkability kill test | **STOP. 0 of the 10 isolated hidden families reach the bar** (≥ 30 same-family pairs above the null's 1e-6 quantile); all ten score **0**, and the best scores 3 even at a 10 × looser threshold. The linkage signal is real (9–38 × enrichment) but ~10³ too small |
| **Step 2 (CF completeness kill test)** | AUC **0.689** in the 8–9 net-MAD bin → Merlin's "in between — report"; **0.755** in the 7–8 bin, but that is an upper bound (the CF trained on 2012) |
| Steps 1b, 2 (trigger ablation) | **NOT RUN** — Step 1a's STOP condition fired |
| **1a POOLED** (2010 + 2012, coordinator-approved follow-up) | still **0 of 10** at the bar (best family 3). Pooling doubles the same-family pairs and raises linkable pairs 11 → 49, a 2.2 × gain in per-pair rate — real, and ~30 × short. § 5 |
| **site-subset ablation** | the far sites DILUTE: dropping the 3 farthest sites multiplies linkable same-family pairs by **4–7 ×** at working thresholds and triples the enrichment (32 × → 68–113 ×). Still 0/10. § 5.4 |
| **Step 6 sizing** | Puget/Olympic: 111 roster stations within 80 km of the Slab2 30–45 km band, **median 30 per cell**, 241 SSE days 2020–2026 → **1.1 TB** pull. Holdout must become gates + station holdout + pre-registered nulls + cross-episode reproducibility. § 6 |

---

## 1. Step 0 — closure bookkeeping

### 1a. Chance-correcting P1b (`scripts/q_p1b_chance.py` → `results/p3/p1b_chance_corrected.csv`)

Every P1b curve (3 models × 2010/2012/2011 × 15 stations × 9 peak rates, 1,053 rows) re-scored with
the ±307 s shift control of `scripts/p2_score.py`. Full tables and discussion:
`P1_RESULTS.md` § ADDENDUM. The headline:

| | raw | chance | corrected |
|---|---|---|---|
| PNET @ 11,600 peaks/station-day, PGC 2011 | 0.4856 | 0.2426 | **0.3209** |
| template cc @ 100 peaks/template-day, PGC | 0.1194 | 0.0000 | **0.1194** |
| **gain** | 4.07 × | | **2.69 ×** |

The template statistic is read at the family's own predicted S, and its chance level measured on O3's
sampled ETS-day negatives at the same threshold is 0.0000 — only the numerator moves. Far boreholes
go 0.356 → **0.123** corrected; B003 at 100 peaks/station-day is indistinguishable from chance.

### 1b. The corrected event-level null, and the uniform-grid expectation (`scripts/q_step0_bc.py`)

`results/p3/p2_corrected_event_null.csv`, `p2_uniform_grid_expectation.json`. Reproduced exactly:
shuffled/true chance-corrected recall ratio **0.600** at 6/day on 2011 (0.616 / 0.707 at 21 / 100 per
day; 0.616–0.745 on dev 2012), against a count ratio of 0.90; uniform-over-grid expectation for
"> 5 km from every family" **64.5 %** against an observed 60.9 % (2011) / 59.1 % (2012). Both are
written into `P2_RESULTS.md` under the heading Merlin asked for.

### 1c. The loud-day collapse (`results/p3/loud_day_collapse.csv`)

| year | doy | Bostock | **P2 @ 6/day** | template T_core |
|---|---|---|---|---|
| 2011 | 239 | 2,262 | **10** | 402 |
| 2011 | 240 | 1,225 | **10** | 201 |
| 2011 | 241 | 1,447 | **5** | 195 |
| 2011 | 242 | 1,169 | **10** | 208 |
| 2011 | other ETS days (mean) | 1,023 | 263 | 241 |
| 2012 | 253 | 1,117 | **1** | 183 |
| 2012 | 255 | 1,843 | **17** | 502 |
| 2012 | 267 | 526 | 549 | 165 |
| 2012 | other ETS days (mean) | 1,428 | 201 | 522 |

The two-sided per-cell per-day MAD is set by the signal on days when tremor is on for most of the
day. 2012 doy 257 and 260 show the same signature (13 and 17 detections at 1,326 and 1,801 events);
doy 267 is the opposite case (few events, so the MAD is a noise estimate again). `lowhalf` exists in
`p2_backproj.py` and was not used in the production runs.

### 1d. The relative-clause rediscovery scorer and its reference table

`scripts/q_rediscover.py` is `g_rediscover.py` with clauses (b) and (c) made RELATIVE:

| clause | old (g_rediscover) | new (q_rediscover) |
|---|---|---|
| (a) centroid | ≤ 5 km | ≤ 5 km (unchanged, absolute) |
| (b) one-to-one recall | ≥ 0.30 | **≥ 0.5 × recall_ref(family)** |
| (c) member purity | ≥ 0.50 | **≥ 0.8 × purity_ref(family)** |
| (d) kept-explained | ≤ 0.50 | ≤ 0.50 (unchanged, absolute) |
| — | — | families with `recall_ref = 0` are reported as **NO-REFERENCE**, never as a failure |

`recall_ref` / `purity_ref` come from a **frozen-template run with SNB + B926 withheld**, on the same
year and window: `scripts/q_run_month.py 2012-0{8,9},2012-10 --drop-stations SNB,B926`
(948 of 1,144 (family, station) templates survive; 13 stations; 47,160 detections ≥ 8 MAD, 10,030
at T_core). Reference table: `results/p3/truth_template_reference.csv` (T_core) and
`truth_template_reference_ext.csv` (T_ext).

**The truth-template per-family reference (T_core, window 2012-08-15 … 2012-10-05, SNB + B926 withheld)**

| stratum | families | Bostock events | detections | **median recall_ref** | **median purity_ref** |
|---|---|---|---|---|---|
| hidden, isolated NN ≥ 5 km | 10 | 2,001 | 771 | **0.196** | **0.567** |
| hidden, mid 2.5–5 km | 10 | 2,561 | 911 | 0.114 | 0.438 |
| hidden, close < 2.5 km | 10 | 2,831 | 794 | 0.128 | 0.456 |
| kept | 100 | 20,756 | 7,459 | 0.072 | 0.482 |

Per isolated hidden family (the stratum the P3 bar is written on):

| bkey | events | det | recall_ref | purity_ref | → needs recall ≥ | → needs purity ≥ |
|---|---|---|---|---|---|---|
| 12 | 391 | 132 | 0.069 | 0.205 | 0.035 | 0.164 |
| 62 | 308 | 37 | 0.071 | 0.595 | 0.036 | 0.476 |
| 113 | 113 | 19 | 0.106 | 0.632 | 0.053 | 0.505 |
| 115 | 124 | 28 | 0.185 | 0.821 | 0.093 | 0.657 |
| 29 | 121 | 81 | 0.190 | 0.284 | 0.095 | 0.227 |
| 15 | 148 | 49 | 0.203 | 0.612 | 0.101 | 0.490 |
| 8 | 243 | 132 | 0.222 | 0.409 | 0.111 | 0.327 |
| 87 | 201 | 104 | 0.279 | 0.538 | 0.139 | 0.431 |
| 119 | 81 | 51 | 0.284 | 0.451 | 0.142 | 0.361 |
| 10 | 271 | 138 | 0.332 | 0.652 | 0.166 | 0.522 |

Two of the thirty hidden families (58 and 125, both in the < 2.5 km stratum) have `recall_ref = 0` on
this network — the truth templates themselves do not detect them — and are NO-REFERENCE. The old
absolute clauses demanded ≥ 0.30 recall and ≥ 0.50 purity, i.e. **more than the truth templates
achieve on 8 of these 10 families**, which is exactly Merlin's flaw #5 and the reason
`g_rediscover.py`'s own calibration run scored the template detector at 9/30.

---

## 2. Step 1a — the pair-linkability kill test: **STOP**

`scripts/q_pairlink.py`; 75 s of H100 for both passes (the whole test is 3.8 × 10⁸ pairs × 13
stations × 121 lags ≈ 1.7 × 10¹⁵ flops, expressed as fp16 matmuls).

### 2.1 What was run

* **Candidates**: every Bostock 2012 ETS-day event (doy 249–268), 27,758 rows, kept + hidden, family
  labels used only afterwards for scoring; 27,681 survive the requirement that all 13 stations have
  valid data (99.7 %).
* **Network**: S15 minus SNB and B926 → 13 stations, 7 sites (PGC-cluster and LRIV-cluster merged,
  divisor 4 each).
* **Windows**: 24 s, 3-C, 20 Hz, starting at `t_origin + c_sta`; `c_sta` = Crosson S from the grid
  centre (48.5 N, −123.75 E, Slab2 21.9 km) + the a5 station correction, one constant per station
  (PGC 10.17 s … NLLB 24.11 s). `t_origin` = Bostock DB time + frozen 2012 family offset.
* **Statistic**: unit-norm 3-C cc for every pair at every station over lags ±3 s (121 lags);
  MAD-normalised per station per (day of *i*, day of *j*) block; site mean; **sum over the 7 sites**;
  **max over lag**.
* **Null**: identical, with each station's record circularly shifted by (station index × 1,800 s).
  Validity is decided on the *true* window in both modes, so null and true see the same event set.
* **Threshold**: the null's 1e-6 quantile = **23.32** (null max 29.76).

**One correction that had to be made before any of this meant anything.** The first null run put
7,472 pairs above 100 with a clean gap below — all of them pairs of candidates **less than a second
apart** (Bostock lists many simultaneous events under different families during a burst). Their
24-s windows overlap almost completely, so their cc is trivially ≈ 1 in the true run *and* in the
null, and they set the whole tail. Pairs closer in time than window + lag (27 s) are therefore
excluded from both runs (`--min-dt`, recorded in the output json); 57,460 of 3.82 × 10⁸ pairs are
removed and the null's 1e-6 quantile falls from 172.06 to 23.32.

### 2.2 The readout

**Same-family pairs above the null's 1e-6 quantile, per isolated hidden family** (bar: ≥ 30 for
≥ 5 of 10; STOP if < 3):

| bkey | 8 | 10 | 12 | 15 | 29 | 62 | 87 | 113 | 115 | 119 |
|---|---|---|---|---|---|---|---|---|---|---|
| same-family pairs available | 28,608 | 35,847 | 75,311 | 10,808 | 7,228 | 45,296 | 19,568 | 6,036 | 7,348 | 2,890 |
| **above 1e-6 threshold** | **0** | **0** | **0** | **0** | **0** | **0** | **0** | **0** | **0** | **0** |
| above 1e-5 (19.96) | 2 | 0 | 0 | 0 | 0 | 3 | 0 | 0 | 0 | 0 |
| above 1e-4 (15.66) | 15 | 19 | 0 | 0 | 0 | 19 | 9 | 1 | 0 | 0 |
| above 1e-3 (11.68) | 130 | 370 | 236 | 14 | 15 | 229 | 154 | 13 | 30 | 15 |

**0 of 10 at the bar → STOP** (`< 3`). Merlin's Step 6 is triggered.

Whole-catalogue totals (3,986,404 same-family pairs, 9,773,831 cross-family pairs < 5 km apart,
369,287,294 cross-family pairs > 5 km apart):

| threshold (null quantile) | same-family above | rate | alias (< 5 km) rate | far cross rate | **enrichment same / far** |
|---|---|---|---|---|---|
| 1e-3 (11.68) | 22,411 | 5.6 × 10⁻³ | 3.7 × 10⁻³ | 6.5 × 10⁻⁴ | **8.7 ×** |
| 1e-4 (15.66) | 950 | 2.4 × 10⁻⁴ | 1.3 × 10⁻⁴ | 1.9 × 10⁻⁵ | **12.6 ×** |
| 1e-5 (19.96) | 78 | 2.0 × 10⁻⁵ | 5.2 × 10⁻⁶ | 8.9 × 10⁻⁷ | **22 ×** |
| **1e-6 (23.32)** | **11** | 2.8 × 10⁻⁶ | 1.0 × 10⁻⁷ | 7.3 × 10⁻⁸ | **38 ×** |

The linkage signal is **real** — same-family pairs are 9–38 × enriched in the tail, monotonically in
threshold, and the enrichment rises with magnitude (at the 1e-3 quantile the rate is 5.3 × 10⁻³ for
pairs whose weaker member is M ≤ 1.4 and 1.7 × 10⁻² for M > 2.0,
`results/p3/pairlink_true_by_magnitude.csv`). It is simply ~10³ too small to build a graph from:
**0.56 % of same-family pairs reach even the 1e-3 quantile, where the far-cross-family false rate
already produces 240,000 spurious edges against 22,000 true ones — an 11 : 1 false-to-true edge
ratio before clustering has begun.**

Figures: `figures/q_pairlink_distributions.png` (true vs null, and true split by family relation),
`figures/q_pairlink_by_family.png` (per-family counts vs the bar, and the enrichment).

### 2.3 Four checks before accepting the STOP

**(i) The null is if anything too strict, and the verdict survives a matched null.** The null's tail
is *heavier* than the true run's (left panel of `q_pairlink_distributions.png`) — the circular shift
moves some windows into quiet parts of the day, where two windows correlate better than two tremor
windows do. Redoing the test with an *in-distribution* null — the 3.69 × 10⁸ cross-family pairs more
than 5 km apart, which see exactly the same data and the same normalisation — gives a 1e-6 threshold
of 19.83 instead of 23.32, 86 same-family pairs above it (2.2 × 10⁻⁵, 21.6 × enrichment), and a
per-family maximum of 3. **Still 0 of 10 at the bar.**

**(ii) The window placement is not what fails** (`scripts/q_pairdiag.py`,
`results/p3/pairlink_window_diag.csv`; doy 251, the 12 most prolific families, 117,973 same-family
and 1,112,216 cross-family pairs). Enrichment of same-family pairs at the cross-family 0.9999
quantile:

| configuration | enrichment |
|---|---|
| A — 24 s at the grid-centre constant (**what the kill test uses**) | 1.9 × |
| B — 24 s at the family's OWN predicted S (oracle placement) | 2.8 × |
| C — 8 s at the family's own S (oracle placement + oracle duration) | **6.0 ×** |
| D — 4 s at the family's own S | 4.2 × |

Giving the search the answer — the family's own cell and an LFE-length window — buys about 3 ×.
The test is short by ~10³.

**(iii) The machinery is sound; it is the pair statistic that is weak** (`scripts/q_stackdiag.py`,
`results/p3/pairlink_stack_audit.csv`). Same events, same windows, same normalisation, same site sum,
same max-over-lag, two statistics:

| statistic | median | baseline | excess |
|---|---|---|---|
| event vs event, same family | 6.42 | 6.10 (cross-family) | **+0.32** |
| event vs its family's own cc-aligned stack (leave-one-out) | 10.72 | 6.10 | **+4.6** |

The data, the ZNE rotation, the frozen time base and the code all deliver the expected matched-filter
excess; the pair statistic is **14 × weaker in excess**, which is the textbook √(signal fraction)
penalty for correlating one noisy event against another instead of against a stack. This is the
physical reason for the STOP, and it is not fixable by tuning: at the observed same-family rate above
a matched 1e-6 threshold (2.2 × 10⁻⁵), reaching 30 linkable pairs needs **≈ 1.4 million same-family
pairs, i.e. ~1,670 events of one family in one ETS**. The most prolific family in the 2012 ETS has
**561**, the median has 190.

**(iv) The alias channel is not what saved it either.** Cross-family pairs less than 5 km apart link
at 65 % of the same-family rate at the 1e-3 quantile and 53 % at 1e-4 — so even if the threshold were
lowered until clusters formed, a majority of the linked pairs would be the wrong family. Alias rate
at the pre-registered threshold: **1.0 × 10⁻⁷** (1 pair in 9.77 M).

### 2.4 Verdict

> **Template-free single-ETS repeater discovery is NOT viable on the S15 permanent network.** With
> Bostock's own event times as the candidate set — the best candidate list any trigger could ever
> produce — and with the 7-site network sum, single events of the same family do not correlate with
> each other above a multiple-comparison threshold. Merlin's **Step 6** applies: re-scope P5 to the
> dense post-2019 regions (Puget 2023+: ~26 stations within 60 km of a cell, i.e. ~4 × the sites and
> a ~2 × lower per-site threshold) and/or accumulate over multiple ETS episodes before stacking.

Steps 1b and the Step-2 trigger ablation were therefore **not run**: every candidate set they compare
(P2 detections, single-station CF peaks, uniform tremor-hour windows) is a *subset or superset of*
the Bostock times already used here, and all of them feed the same pair-linking stage that has just
failed by three orders of magnitude on the strongest possible input. Running them could only produce
a worse number.

---

## 3. Step 2 — the CF completeness kill test (run: it is independent of the repeater chain)

`scripts/q_cf_completeness.py` → `results/p3/cf_completeness.{csv,json}`,
`cf_completeness_detections.parquet`.

Candidates = every stored 2012 ETS-day detection of the frozen template detector
(`../lfe_extension/results/e5r_detections_a1_mergeB_L12.parquet`, L = 12, net_mad ≥ 7), 60,671 with
≥ 5 sites of CF. Feature = the learned CF (log-odds, **lowhalf**-normalised per station-day, averaged
within the PGC and LRIV clusters, summed over the 9 sites) read at each station's predicted S for
that detection's own family. Truth = the detection's own-family one-to-one Bostock match
(`match_own`, Track D's ±2 s greedy matcher; overall match rate 0.182).

| net_mad bin | n | matched | unmatched | **AUC** | median CF (matched) | median CF (unmatched) |
|---|---|---|---|---|---|---|
| 7–8 (below T_core) | 42,143 | 3,308 | 38,835 | **0.755** | 6.90 | 1.66 |
| **8–9** | 10,375 | 2,596 | 7,779 | **0.689** | 7.18 | 3.35 |
| ≥ 9 | 8,153 | 5,126 | 3,027 | **0.645** | 8.76 | 5.75 |
| all | 60,671 | 11,030 | 49,641 | 0.762 | 7.78 | 2.13 |

**Verdict on Merlin's rule: 0.689 in the 8–9 bin → "in between — report".** The thread is neither
closed nor cleared. Three things must be said with the number:

1. This is an **upper bound**: the CF was trained on 2012 Bostock positives, so it has seen these
   very events. Any retrain decision needs an out-of-sample repeat (2010 is the free dev year).
2. The AUC *rises* as the threshold falls — 0.755 in the 7–8 band, which is precisely the band a
   completeness gain would have to come from — but that band is also where `match_own` is least
   reliable as truth, because Bostock's catalogue is incomplete and an unmatched detection there is
   not necessarily false. That bias pushes the AUC **down**, while the training overlap pushes it
   **up**; the two are not separable from this run.
3. The "all" row (0.762) is inflated by the correlation between net_mad and realness across bins and
   should not be quoted.

Recommendation (not executed, per the spec): **do not retrain yet.** Repeat this exact measurement on
2010, where the CF saw no labels, and use Track E's realness framework rather than `match_own` as the
target. If the 2010 8–9-bin AUC also lands near 0.7, the CF is worth ~one extra tier of template
threshold and belongs in Track E as a feature — not as a detector.

The other two Step-2 arms (trigger ablation over P2 detections / single-station CF peaks / uniform
tremor-hour windows) are **not run**, for the reason in § 2.4.

---

## 4. Files

| file | what |
|---|---|
| `scripts/q_p1b_chance.py` | Step 0a — ±307 s chance correction of every stored P1b scan |
| `scripts/q_step0_bc.py` | Steps 0b/0c — corrected event-level null, uniform-grid expectation, loud-day table |
| `scripts/q_rediscover.py` | Step 0d — relative-clause P3 scorer + `truth_reference()` |
| `scripts/q_run_month.py` | frozen matched filter with a `--drop-stations` / `--bank` override (monkey-patch; `o_*` untouched) |
| `scripts/q_pairlink.py` | Step 1a — the GPU pair-linkability test and its circular-shift null |
| `scripts/q_pairdiag.py`, `q_stackdiag.py` | Step 1a audits — window ablation, pair-vs-stack |
| `scripts/q_cf_completeness.py` | Step 2 — the CF completeness kill test |
| `scripts/q_figs.py` | the two Step-1a figures |
| `results/p3/p1b_chance_corrected.csv` | 1,053 rows: raw / chance / corrected recall, 3 models × 3 years × 15 stations × 9 rates |
| `results/p3/p2_corrected_event_null.csv`, `p2_uniform_grid_expectation.json`, `loud_day_collapse.csv` | Step 0b/0c products |
| `results/p3/truth_template_reference{,_ext}.csv` | the truth-template per-family reference (SNB + B926 withheld) |
| `results/p3/trackQ/` | the withheld-station frozen-template run on 2012 (detections, family-days, stacks, run json) |
| `results/p3/pairlink_{null,true}*.{json,npz,csv,parquet}` | Step 1a: histograms, per-family counts at four thresholds, magnitude table, the 458 stored edges, event table |
| `results/p3/pairlink_{window_diag,stack_audit}.csv` | the two audits |
| `results/p3/cf_completeness.{csv,json}` | the CF completeness AUC table |
| `figures/q_pairlink_distributions.png`, `q_pairlink_by_family.png` | Step 1a figures |
| `scripts/q_pairlink2.py`, `q_cache2010.py`, `q_pool_report.py`, `q_figs2.py` | § 5: the pooled 2010 + 2012 repeat, the 2010 borehole cache build, the report assembly, the figure |
| `results/p3/pairlink2_{pool,y2012}_{null,true}*`, `pool_*.csv/json` | § 5 products (thresholds, per-family counts per site configuration, ablation, extrapolation) |
| `figures/q_pairlink_pooled.png` | § 5 figure |

`results/p3/pairlink_{smoke,smoketrue,nulldiag}*` are the development smoke runs (2 and 4 days) and
the diagnostic that located the time-overlap contamination; kept, not deleted, and not part of any
result.

---

## 5. Step 1a, POOLED — 2010 + 2012 together, and the site-subset ablation

`scripts/q_pairlink2.py` (the Step-1a recipe generalised to several years and to arbitrary site
subsets), `scripts/q_pool_report.py`, `scripts/q_figs2.py`; figure
`figures/q_pairlink_pooled.png`; products `results/p3/pairlink2_*`, `pool_*.csv/json`.
Dev years only — 2011 and 2013 were not touched, `referee_ledger.csv` is unchanged.

### 5.1 What had to be built first

The picker's memory-mapped day store has **no B001–B007 for 2010**, so the 13-station network cannot
be run on that year out of the box. The raw day traces do exist read-only in
`../lfe_extension/traces` (doy 222–256), so `scripts/q_cache2010.py` built the missing 209
station-days of 20 Hz cache by calling `o_cache_month.cache_one` as a library — the frozen chain,
unmodified, writing into the picker's own `cache20hz/`. All 209 built clean (100 Hz native,
full days, component order `EHZ,EH1,EH2`, 100 % non-zero).

`q_pairlink2.py` then reads the **20 Hz caches** rather than the mm store, in both years. Those
caches hold raw counts (CN.PGC 2012 peaks at 1.6 × 10⁷), which **overflow float16**; the first run
returned NaN thresholds for that reason. Each window is now rms-scaled before the float16 cast —
free, because every downstream normalisation is per row — and a finiteness assertion guards the cut.

**Control.** Running the pooled code on 2012 alone reproduces the original mm-store result of § 2
essentially exactly: same-family pairs 3,994,687 vs 3,986,404, and pairs above the 1e-6 quantile
**11 vs 11**, above 1e-5 **78 vs 78**, above 1e-4 972 vs 950. The cache path, the ZNE-vs-raw
component order and the mm path agree, as they must (a common orthogonal rotation leaves the 3-C
inner product of two windows at the same station unchanged).

### 5.2 The pooled run

39,563 usable candidates (**11,853 from 2010** + **27,710 from 2012**; 2010 contributes 16 of its 17
ETS days, 2012 all 20), **782,521,738 pairs** after the 27-s overlap exclusion, of which 8,118,843
same-family, 18,411,064 cross-family < 5 km apart and 755,991,831 cross-family > 5 km apart.
Pooled null 1e-6 quantile **22.59** (2012 alone: 23.31). Cross-year pairs are what does the work:
pooling **2.03 ×** the events gives **2.03 ×** the same-family pairs.

**Same-family pairs above threshold, per isolated hidden family** (bar ≥ 30 for ≥ 5 of 10):

| bkey | events 2010 | events 2012 | pooled | same-family pairs | @1e-4 | @1e-5 | **@1e-6** | (@1e-6, 2012 only) |
|---|---|---|---|---|---|---|---|---|
| 8 | 109 | 240 | 349 | 60,620 | 32 | 6 | **0** | 0 |
| 10 | 267 | 269 | 536 | 143,004 | 163 | 11 | **1** | 0 |
| 12 | 94 | 389 | 483 | 116,226 | 6 | 1 | **1** | 0 |
| 15 | 188 | 148 | 336 | 56,115 | 0 | 0 | **0** | 0 |
| 29 | 23 | 121 | 144 | 10,262 | 0 | 0 | **0** | 0 |
| 62 | 18 | 302 | 320 | 50,881 | 28 | 5 | **0** | 0 |
| 87 | 221 | 199 | 420 | 87,715 | 69 | 6 | **3** | 0 |
| 113 | 91 | 111 | 202 | 20,190 | 0 | 0 | **0** | 0 |
| 115 | 136 | 122 | 258 | 33,062 | 11 | 1 | **0** | 0 |
| 119 | 23 | 77 | 100 | 4,905 | 2 | 0 | **0** | 0 |

**0 of 10 reach the bar; the best family reaches 3.** The STOP verdict of § 2.4 stands.

The lever does work, in the direction and by roughly the size expected — it is simply far too small:

| | 2012 only | **pooled 2010 + 2012** | gain |
|---|---|---|---|
| same-family pairs | 3,994,687 | 8,118,843 | 2.03 × |
| same-family pairs above the 1e-6 quantile | 11 | **49** | **4.5 ×** |
| per-pair rate above the 1e-6 quantile | 2.75 × 10⁻⁶ | **6.04 × 10⁻⁶** | 2.19 × |
| best isolated hidden family | 0 | **3** | — |

Enrichment over far cross-family pairs, pooled: **8.4 × / 16.8 × / 25.7 × / 32.1 ×** at the
1e-3 / 1e-4 / 1e-5 / 1e-6 quantiles (2012 alone: 8.7 / 12.6 / 22 / 38 ×). The alias rate is
4.9 × 10⁻³ / 2.3 × 10⁻⁴ / 1.1 × 10⁻⁵ / 1.3 × 10⁻⁶ at the same four thresholds — still 55–67 % of the
same-family rate at the two loose thresholds, so the wrong-family channel does not go away with
pooling either.

### 5.3 The arithmetic

At the pooled isolated-stratum rate above the 1e-6 quantile (5 pairs in 582,980 = 8.6 × 10⁻⁶), a
family needs **3.5 million same-family pairs ≈ 2,650 events** to put 30 pairs above threshold. The
median isolated hidden family has **328 pooled events over two ETS episodes (~164 per episode)** and
the most prolific family in the pooled set has 759. So the requirement is

> **≈ 16 ETS episodes** for a median isolated family, **≈ 10** for the most prolific, at the present
> 7-site network — against the 2 episodes pooled here.

That is the quantitative form of the § 2.4 verdict, and it is also the first number that makes the
Step-6 branch look like a real option rather than a hope: a 2020–2026 Puget window contains ~8–10
ETS episodes.

### 5.4 Do the far sites dilute? **Yes — clearly.**

Every per-site lag-resolved z was kept, so 14 site-subset configurations were scored from the same
matmuls, each with its own null-derived thresholds. Each family is scored in the subset chosen for
it. `n4_*` = the 4 sites nearest that family among the 7 discovery sites; `o4_*` = the coordinator's
literal list, PGC-cluster + B926 + VGZ + the family's nearest far borehole. **Both use the family's
position, which a template-free search does not know, and `o4_*` puts the withheld station B926 back
in — they are dilution diagnostics, not candidate detectors.**

Isolated hidden families, same-family pairs above threshold (582,980 pairs available):

| configuration | @1e-3 | @1e-4 | @1e-5 | @1e-6 | families at the bar |
|---|---|---|---|---|---|
| **7 discovery sites (protocol)** | 4,542 | 311 | 30 | 5 | **0/10** |
| 4 nearest of those 7 (oracle site choice) | 7,708 | **1,306** | **221** | 4 | **0/10** |
| PGCq + B926 + VGZ + nearest far borehole (oracle + withheld) | **9,435** | 1,264 | 153 | 8 | **0/10** |

Whole catalogue, same-family pairs above threshold (8,118,843 available): 58,988 / 3,506 / 270 / 49
at 7 sites, against 56,159 / 6,457 / **1,069** / 22 for the nearest-4 and 81,511 / 9,005 / **1,038** /
47 for the oracle-4. Per-configuration enrichment at the 1e-6 quantile reaches **113 ×**
(`n4_B001+LRIVq+PGCq+VGZ`, 38 families) and **68 ×** (`o4_B926+LRIVq+PGCq+VGZ`, 54 families) against
**32 ×** for the protocol 7.

So at every working threshold (1e-3 to 1e-5) **dropping the three far, low-SNR sites multiplies the
number of linkable same-family pairs by 4–7 × and roughly triples the enrichment.** The far sites
(NLLB at 82 km, B001 at 68, B004 at 61 from the family field) add noise to the network sum faster
than they add signal. At the extreme 1e-6 quantile the counts converge again, because a 4-site sum
has a relatively heavier tail and its own threshold rises to meet it. Nothing reaches the bar in any
configuration.

This is the most useful thing the ablation says about Step 6: **what the re-scope needs is stations
CLOSE to the source, not merely more of them.**

### 5.5 What √N_sites would buy

Fitting the same-family shift against the null from the four tail quantiles gives
**δ = 0.38 σ** at 7 sites (σ = the spread of the null network statistic). The network sum grows like
N_sites in signal and √N_sites in noise, so δ ∝ √N_sites:

| sites | δ (σ) | rate above the 1e-6 quantile | same-family pairs for 30 | events of one family |
|---|---|---|---|---|
| 7 (now) | 0.38 | 6.2 × 10⁻⁶ | 4.9 M | **3,120** |
| 10 | 0.46 | 8.7 × 10⁻⁶ | 3.5 M | 2,630 |
| 14 | 0.54 | 1.3 × 10⁻⁵ | 2.4 M | 2,180 |
| 20 | 0.65 | 2.0 × 10⁻⁵ | 1.5 M | **1,730** |
| 26 | 0.74 | 3.0 × 10⁻⁵ | 1.0 M | 1,430 |

Going 7 → 20 sites cuts the events-per-family requirement by only **1.8 ×** (3,120 → 1,730). Site
count alone therefore does not close the gap; combined with episode pooling it brings the
requirement to **~7–11 ETS episodes** for a prolific family, which is what a 2020–2026 Puget window
offers — and that is before any credit for the closer stations, which § 5.4 says is the bigger term.

---

## 6. Step 6 — what a post-2019 Puget / Oregon re-scope would actually need

Merlin's Step 6 is the branch Step 1a selects. This section sizes it; nothing here has been run.

### 6.1 The network that replaces S15

Roster = `../data/candidate_stations_post2020.csv` (277 post-2020 candidate stations, margin-wide).
Target = the Slab2 **30–45 km depth band** (the LFE depth range), thinned to a ~5 km grid, and the
stations within **80 km** of it:

| segment | band cells | roster stations in segment | **stations within 80 km of the band** | stations within 80 km **per cell** (p10 / median / p90) |
|---|---|---|---|---|
| Puget / Olympic, 47.0–48.5 N | 545 | 122 | **111** | 9 / **30** / 35 |
| Oregon, 43.5–46.0 N | 123 | 83 | **67** | 0 / **12** / 21 |
| combined, 43.5–48.5 N | 721 | 180 | **171** | 7 / **28** / 35 |

For comparison, the S15 network that failed Step 1a offers **7 sites** (9 before withholding SNB and
B926) at 26–82 km from the family field, with only the PGC cluster inside 30 km. Puget/Olympic is
the segment that actually changes the problem; **Oregon at a median of 12 stations per cell is not
obviously better than what we already have** and should not be pulled first.

### 6.2 The data pull

SSE windows = days with ≥ 50 PNSN tremor detections inside the segment, 2020–2026
(`../catalogs/pnsn_tremor_cascadia_full.csv`): **241 days** (Puget/Olympic), **193 days** (Oregon),
**421 days** in the union. Median day file on disk for a 3-C station-day is **20.5 MB**; the frozen
20 Hz float32 cache costs another **20.7 MB** per station-day.

| scope | station-days | raw mseed | + 20 Hz cache | total |
|---|---|---|---|---|
| Puget / Olympic only | 111 × 241 = 26,751 | 548 GB | 554 GB | **1.10 TB** |
| Oregon only | 67 × 193 = 12,931 | 265 GB | 268 GB | 0.53 TB |
| both, on the union of days | 171 × 421 = 71,991 | 1.48 TB | 1.49 TB | **2.97 TB** |

Disk is not the constraint (72 TB free). The constraint is FDSN wall time: P0-O pulled 1,566
station-days in ~6 minutes at ≤ 8 concurrent requests (~4 station-days/s), so 27,000 station-days is
a few hours and 72,000 is most of a day, with retries. **Recommendation: pull Puget/Olympic only
(1.1 TB, ~26,751 station-days) and only for the SSE windows.** Continuous data is not needed until a
family bank exists.

### 6.3 The holdout design, with no Bostock truth

South and east of Bostock's field there is no event catalog, so the P0-O design — 30 hidden families
scored by one-to-one recall and purity against a truth catalog — **cannot be reproduced**. What is
left, in descending order of strength:

1. **Gates 0–3 are the only judge** (`results/GATES_RESULTS.md`). They are detector-agnostic and
   already calibrated: Gate 1 power 0.77 / 0.91 / 0.92 at N = 20 / 50 / 100 with false pass 1/150
   over three negative mechanisms; Gate 0 alias false-pass 12 % under the adopted rule; Gate 2
   essentially untestable below N ≈ 150–200; Gate 3 needs ≥ 2 × 30 d of calendar — which a 2020–2026
   span **has** and a single ETS does not, so Gate 3 becomes usable for the first time.
2. **Station holdout replaces family holdout.** Withhold ≥ 2 sites from discovery for the whole run
   and require Gate 1 on them. With 30 stations per cell this costs almost nothing, and it is the
   only leakage-proof holdout available without truth.
3. **Pre-registered nulls become the primary evidence, not a supporting check** (Merlin Step 3):
   N-shift (the whole chain on station-shifted data, 5 seeds), N-quiet (the chain on quiet days),
   N-moveout. Claim rule for any gate-PASS family: true count ≥ N-shift max + 3 **and** N-shift mean
   ≤ 10 % of true count; otherwise "candidate", never "discovery".
4. **Cross-episode reproducibility as an internal holdout**: discover on episodes A, require the
   family to re-fire in episodes B. 2020–2026 gives ~8–10 independent ETS episodes per segment, so
   this is a genuine out-of-sample test in time and is the closest available substitute for a
   held-out truth catalog.
5. **Literature seeds, if the existence test passes first.** Ducellier 2004–2011 (66 templates /
   37 located families) covers the Olympic segment and Shelly 2018–2024 (27 families, 61,441
   detections) the far south; Shelly's catalog is on disk. **But P0-O § O4 is still DEFERRED**: none
   of B045, B935, RVIT, PETY exists in any trace tree here, so "do the roster stations see these
   families at all" has never been tested. That test must pass before any seed is used as truth and
   before any southern recall number is quoted.

### 6.4 What 7 → ~20 sites is worth, and what it is not

The network statistic is a sum of per-site MAD-z values, so coherent signal grows like N_sites and
the noise like √N_sites: the per-pair network SNR gains **√(20/7) = 1.69 ×**, equivalently the cc a
single site must show to clear a fixed network threshold falls by the same 1.69 ×. The quantitative
projection from the pooled measurement is in § 5.4. The honest summary is that **site count alone
does not close the gap**: the gap Step 1a measures is ~10³ in tail probability, and √N over a
realistic 7 → 20–26 improves the same-family/null separation by well under an order of magnitude in
rate. The parts of the Puget re-scope that could matter more than the count are (a) **distance** —
several stations within 20–30 km of every cell instead of one cluster at 26 km and five sites at
45–82 km, which raises per-site SNR directly rather than as √N, and (b) **episode count** — ~8–10
ETS episodes instead of one, which is quadratic in pairs (§ 5). Those two, together, are the case
for the re-scope; the site count on its own is not.
