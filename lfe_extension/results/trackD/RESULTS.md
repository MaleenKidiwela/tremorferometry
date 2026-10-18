# Track D — the product: completeness model, catalog v1.1, and the science QC

Generated 2026-09-17 11:32 UTC by `scripts/d4_results.py` from the files in `results/trackD/`. Executes steps 6–8 of `results/merlin_review_plan_2026-09-16.md`. Inputs: Track C's monthly outputs (read-only), Track A's frozen operating point and referee runs, Track A-2's relocations and epoch policy.

## 0.0 Catalog v1.1 — what changed since v1, and why

Two independent changes, both upstream of Track D, both from Track C:

**(1) The cc-floor bug.** Track C found and fixed a floor bug in the per-site cross-correlation and re-ran **10 months**, reporting 486 artifact rows removed (`results/trackC/rerun_ccfloor_log.csv`). Read off that log directly, across those ten months: rows with |cc| > 1 **233 → 0**, rows with net_mad > 50 **205 → 0**, Track-E-flagged rows **179 → 0**, and recorded detections 23,160 → 21,851 (1309 fewer; the maximum |cc| anywhere goes from 3.99e+11 to 0.38). This is a DATA change: those months' detection files are different, so every Track D product built on them is rebuilt here.

**(2) The tiers were re-derived on the corrected data** by the same frozen false-alarm-only rule. The old table is kept as `tiers_by_epoch.pre_ccfloor.csv`:

| epoch | T_core_v1 | T_ext_v1 | T_core_v1_1 | T_ext_v1_1 |
|---|---|---|---|---|
| E0 | 8.75 | 8.25 | 8.75 | 8.25 |
| E1 | 8.25 | 8.0 | 9.0 | 8.0 |
| E2 | 8.75 | 8.0 | 9.0 | 8.0 |
| E3 | 9.25 | 8.0 | 9.25 | 8.0 |

E0 and E3 are unchanged; **E1's T_core moves 8.25 → 9.00 and E2's 8.75 → 9.00**. That is the dominant effect on the catalog — a threshold change, not a data change.

**Counts, per epoch and per tier:**

| epoch_id | core_v1 | ext_v1 | core_v1_1 | ext_v1_1 | d_core | d_ext |
|---|---|---|---|---|---|---|
| E0 | 19154 | 40222 | 18958 | 39841 | -196 | -381 |
| E1 | 26076 | 38358 | 10660 | 38334 | -15416 | -24 |
| E2 | 15649 | 38120 | 12307 | 38093 | -3342 | -27 |
| E3 | 32412 | 128115 | 32405 | 128107 | -7 | -8 |
| ALL | 93291 | 244815 | 74330 | 244375 | -18961 | -440 |

`cat_core` falls 93,291 → 74,330 (−20.3 %), almost all of it in E1 and E2 and almost all of it the threshold change; `cat_ext` is essentially unchanged (244,815 → 244,375, −0.18 %) because T_ext moved only for epochs that were already at 8.00. The E0 and E3 deltas (-196 and -7 core events) are the cc-floor fix alone, since their thresholds did not move — that is the clean measure of the bug's size on the catalog: well under 1 %.

**What this does NOT change.** The 2012 referee is Track A data and is untouched by both, so every referee-based result is bit-identical: the completeness-model validation bar, the two pre-registered Thomas bars, the thinned-Bostock reference, and the frozen α of the §3b-2 front estimator (no redesign was done — the frozen α = 0.05 was applied as-is).

**Verdicts: none flipped.** The per-verdict comparison is in §0 and the sections below; the changes are in the completeness NUMBERS (E1 T_core recall 0.184 → 0.129, E2 0.146 → 0.129) and in the floor λ (E1 and E2 quiet-day rates fall 1.5–4× with the higher T_core), not in any pass/fail.

Version v1 is preserved in full: its catalog, family table and family-day table under their own names, and its completeness / floor / Thomas outputs in `results/trackD/v1_snapshot/`.

## 0. Headline

* **Catalog**: 74,330 `cat_core` events and 244,375 `cat_ext` events over 115 families, 2014-01-01 .. 2026-09-15, from 153 finished Track C months (the complete 2014-01..2026-09 record).
* **Completeness-model validation bar (±0.05 on the held-out 2010 year): MET** — largest |predicted − measured| own recall over all magnitude bins and both tiers = 0.0236.
* **Thomas bar — main-front first-burst speed, cat_core: MET** (ours -5.275 vs Bostock 2012 -4.695, 12.4 % vs a 20 % tolerance).
* **Thomas bar — re-rupture interval median, cat_core: NOT MET** (ours 1.051 vs Bostock 2012 0.7039, 49.4 % vs a 5 % tolerance).
* **Thomas bar — main-front first-burst speed, cat_ext: MET** (ours -5.228 vs Bostock 2012 -4.695, 11.4 % vs a 20 % tolerance).
* **Thomas bar — re-rupture interval median, cat_ext: NOT MET** (ours 0.8663 vs Bostock 2012 0.7039, 23.1 % vs a 5 % tolerance).
* **Main fronts recovered (§3b-2, post-registration estimator, approved 2026-09-17)**: 27 of 31 SSE windows fitted and **12 significant** at `cat_core` (the pre-registered estimator managed 3), 6 propagating SE and 6 NW, |speed| 1.3–8.6 km/day.

### Network epochs and the tier thresholds in force

**Track C is the authority on both.** Source: `/home/jovyan/tremorferometry/lfe_extension/results/trackC/tiers_by_epoch.csv (Track C, FA-only rule; 4 epochs)`. Track C defines the network epochs and re-derives T_core and T_ext for each by the frozen false-alarm-only rule (T_core: quiet-day FA ≤ 0.1 per template-day; T_ext: quiet-day ≤ 1.0 AND post-SSE ≤ 1.5), applying them as a re-flag of the `cat_core` / `cat_ext` columns. Every Track D step — the epoch ablation, the Mc fit, the floor λ, the catalog — reads the same epochs and the same numbers through `scripts/d_tiers.py`, so the completeness model can never disagree with the flags. A row whose threshold could not be derived yet (too few quiet family-days) falls back to the frozen referee point 8.75 / 8.25 and is marked so.

| epoch | t0 | t1 | station_set | n_stations | T_core | T_ext | source | status |
|---|---|---|---|---|---|---|---|---|
| E0 | 2014-01-01 | 2017-02-28 | S15 | 15 | 8.75 | 8.25 | Track C, derived | derived (derived 9.00/8.00, within 0.25 -> frozen anchor kept) |
| E1 | 2017-03-01 | 2018-12-31 | S15 | 15 | 9.0 | 8.0 | Track C, derived | derived |
| E2 | 2019-01-01 | 2021-03-31 | S14 | 14 | 9.0 | 8.0 | Track C, derived | derived |
| E3 | 2021-04-01 | 2099-12-31 | S13 | 13 | 9.25 | 8.0 | Track C, derived | derived |

**An independent confirmation.** Track D's 2012 station-set ablation says the 13-station network (S13) needs **T = 9.33** to reproduce the frozen quiet-day false-alarm rate at T_core. Track C's FA-only re-derivation, on the production record and with no reference to that ablation, lands on **T = 9.25** for epoch E3. Two different data sets and two different procedures agree to 0.08 MAD — the strongest single check that the completeness model and the flags describe the same detector.

Track D measures per-family recall by ablating the **station set**, which is the thing that actually changes the detector; two of Track C's epochs (E0 and E1) share the 15-station set, so the ablation axis is S15 / S14 / S13 and `d_tiers.py` maps each epoch onto it. `completeness_family_epoch.csv` is the map already expanded onto Track C's epoch axis; `completeness_family_stationset.csv` is the measurement itself.

## 1. Completeness model (Merlin Step 6)

### 1a. Epoch / station-set ablation

The three station sets the production run realises (see the epoch table above for the mapping onto Track C's network epochs):

| station set | stations | what changed |
|---|---|---|
| S15 | 15 | the frozen 2012 network (Track C epochs E0, E1) |
| S14 | 14 | UW.LRIV dropped from 2019-05-30, A-2 epoch-identity FAIL (Track C epoch E2) |
| S13 | 13 | LRIV dropped + CN.SNB ends 2021-03-28 (Track C epoch E3) |

**CN.PFB (ends 2018-12-20) carries no 2012 template and is not in the network at all**, so that boundary has no effect on this detector and needed no ablation. Four GPU passes (`scripts/d1_run_ablation.sh`) re-ran the frozen configuration with the S14 and S13 station sets on the 2012 test half and on the held-out 2010 year; the frozen HF veto is applied to all of them. `P9` is the 9-station persistent set (= the 2010 network), Track A step 2's control.

| station_set | ref_year | tier | n_fam | n_events | recall@frozen T | quiet FA @frozen T | T @matched FA | recall@matched FA | recall_fam_median | recall_fam_p10 | recall_fam_p90 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| S13_2010 | 2010 | T_core | 90.0 | 15931.0 | 0.1771 | 0.3559 | 8.9231 | 0.1645 | 0.1631 | 0.0179 | 0.3168 |
| S14_2010 | 2010 | T_core | 90.0 | 15931.0 | 0.1738 | 0.2968 | 8.7636 | 0.173 | 0.1494 | 0.0171 | 0.3206 |
| S15_2010 | 2010 | T_core | 90.0 | 15931.0 | 0.1718 | 0.2919 | 8.75 | 0.1718 | 0.1436 | 0.0171 | 0.3199 |
| S13_2010 | 2010 | T_ext | 90.0 | 15931.0 | 0.2201 | 0.7106 | 8.4216 | 0.2043 | 0.2056 | 0.0317 | 0.3938 |
| S14_2010 | 2010 | T_ext | 90.0 | 15931.0 | 0.214 | 0.5653 | 8.2628 | 0.213 | 0.197 | 0.0306 | 0.3853 |
| S15_2010 | 2010 | T_ext | 90.0 | 15931.0 | 0.2107 | 0.5554 | 8.25 | 0.2107 | 0.1879 | 0.0306 | 0.3843 |
| P9 | 2012 | T_core | 130.0 | 13381.0 | 0.1582 | 0.1143 | 8.8202 | 0.1532 | 0.0947 | 0.0 | 0.3098 |
| S13 | 2012 | T_core | 130.0 | 13381.0 | 0.1526 | 0.2901 | 9.3272 | 0.1185 | 0.0994 | 0.0 | 0.301 |
| S14 | 2012 | T_core | 130.0 | 13381.0 | 0.146 | 0.0977 | 8.7528 | 0.1457 | 0.0931 | 0.0 | 0.2992 |
| S15 | 2012 | T_core | 130.0 | 13381.0 | 0.1463 | 0.0971 | 8.75 | 0.1463 | 0.0909 | 0.0 | 0.2979 |
| P9 | 2012 | T_ext | 130.0 | 13381.0 | 0.1977 | 0.3258 | 8.3223 | 0.1911 | 0.1383 | 0.0 | 0.3534 |
| S13 | 2012 | T_ext | 130.0 | 13381.0 | 0.1909 | 0.7443 | 8.7763 | 0.1505 | 0.1381 | 0.0 | 0.3472 |
| S14 | 2012 | T_ext | 130.0 | 13381.0 | 0.1841 | 0.279 | 8.2539 | 0.1837 | 0.1354 | 0.0 | 0.3534 |
| S15 | 2012 | T_ext | 130.0 | 13381.0 | 0.1841 | 0.2766 | 8.25 | 0.1841 | 0.1334 | 0.0 | 0.3526 |

**Reading (Merlin F3).** At a FIXED MAD threshold a smaller network is not a fairer detector, it is a lower one. Dropping LRIV (S15 → S14) costs nothing on either axis — recall and false alarms are unchanged — which independently supports Track A-2's decision to drop it. Losing SNB in 2021-03 (S14 → S13) is the real break: at the frozen thresholds that network sits at ~3× the quiet-day false-alarm rate the operating point was frozen at, and once that is matched it loses ~19 % of its recall. The held-out 2010 year reproduces the direction:

| tier | ablation | ratio_2012 | ratio_2010 | ratio_2012_from15 | cross_year_delta |
|---|---|---|---|---|---|
| T_ext | drop LRIV | 0.9312 | 1.0155 | 1.0 | 0.0843 |
| T_ext | drop LRIV+SNB | 0.966 | 1.0444 | 1.0374 | 0.0784 |
| T_core | drop LRIV | 0.9225 | 1.0117 | 0.998 | 0.0892 |
| T_core | drop LRIV+SNB | 0.9646 | 1.0311 | 1.0434 | 0.0665 |

**The DIRECTION transfers, the MAGNITUDE does not.** The E3 false-alarm penalty is 3.0x on 2012 but only 1.2x on 2010 at T_core — the 2010 'strict quiet' set is seven days, four of them post-SSE, and it already sits at 3x the 2012 quiet-day rate, so it has little headroom. **The ablation is therefore the MECHANISM, not the calibration**: the number attached to each event is λ measured directly from the PRODUCTION quiet days of that epoch (section 1c), not the ablation ratio.

Per-family recall for every epoch and both tiers, with Wilson intervals and the matched-FA value, is `results/trackD/completeness_family_epoch.csv`; the per-epoch thresholds are `completeness_matched_fa.csv`.

### 1b. Daily sensitivity — Mc(family, day)

Model (full statement in `scripts/d1_mc_model.py`):

```
P(detect | family f, day d, magnitude M) = Phi( b_f + alpha*(M - 1.8) + c1*u(f,d) + c2*v(f,d) )
Mc(f,d,T) = 1.8 - (b_f + c1*u + c2*v) / alpha
```

`b_f` is a ridge-shrunk per-family intercept fitted separately for each epoch's station set (so it carries component 1 as well); `u` = log10(family-day MAD scale / its epoch median); `v` = the SNR-weighted log10 noise-amplitude factor built from the daily median 8-s window energy, normalised per station per INSTRUMENT epoch so the 2017 CNSN gain changes cancel.

| tier | station_set | b0 | alpha | c1 | c2 | lam |
|---|---|---|---|---|---|---|
| T_core | P9 | -1.1893 | 1.2493 | 0.2618 | -0.1736 | 1.0 |
| T_core | S13 | -1.2225 | 1.2073 | 0.4915 | -0.3679 | 1.0 |
| T_core | S14 | -1.2438 | 1.2234 | 1.5146 | -0.2425 | 1.0 |
| T_core | S15 | -1.1625 | 1.2013 | 1.3771 | -0.1886 | 1.0 |
| T_core | S15_2010 | -1.0627 | 1.7034 | -4.9765 | -0.0062 | 1.0 |
| T_ext | P9 | -0.9702 | 1.13 | 1.4408 | -0.1295 | 1.0 |
| T_ext | S13 | -0.8835 | 1.0654 | 0.7508 | -0.1374 | 1.0 |
| T_ext | S14 | -0.9592 | 1.0487 | 1.9277 | -0.0823 | 1.0 |
| T_ext | S15 | -1.0381 | 1.1297 | 2.0047 | -0.1006 | 1.0 |
| T_ext | S15_2010 | -0.8424 | 1.642 | -3.3875 | 0.0203 | 1.0 |

**Validation — the pre-registered bar.** Fit on the 2012 ablation whose station set equals the 2010 network (9 stations), predicted on the HELD-OUT 2010 year:

| tier | bin | n | measured | predicted | delta | within_0.05 |
|---|---|---|---|---|---|---|
| T_ext | M<1.6 | 3382 | 0.1517 | 0.1753 | 0.0236 | True |
| T_ext | M1.6-1.8 | 6315 | 0.2067 | 0.193 | -0.0137 | True |
| T_ext | M1.8-2.0 | 4788 | 0.2362 | 0.2293 | -0.0069 | True |
| T_ext | M>=2.0 | 1446 | 0.2822 | 0.2869 | 0.0048 | True |
| T_ext | all | 15931 | 0.2107 | 0.2086 | -0.0021 | True |
| T_ext | M<1.6 (2012 in-sample) | 3981 | 0.1723 | 0.1704 | -0.0019 | True |
| T_ext | M1.6-1.8 (2012 in-sample) | 4769 | 0.1948 | 0.1896 | -0.0052 | True |
| T_ext | M1.8-2.0 (2012 in-sample) | 3462 | 0.2085 | 0.2138 | 0.0053 | True |
| T_ext | M>=2.0 (2012 in-sample) | 1169 | 0.2635 | 0.2662 | 0.0027 | True |
| T_ext | all (2012 in-sample) | 13381 | 0.1977 | 0.1968 | -0.0008 | True |
| T_core | M<1.6 | 3382 | 0.11 | 0.1187 | 0.0087 | True |
| T_core | M1.6-1.8 | 6315 | 0.1465 | 0.1386 | -0.0079 | True |
| T_core | M1.8-2.0 | 4788 | 0.1773 | 0.1697 | -0.0077 | True |
| T_core | M>=2.0 | 1446 | 0.2158 | 0.2248 | 0.0091 | True |
| T_core | all | 15931 | 0.1543 | 0.1515 | -0.0027 | True |
| T_core | M<1.6 (2012 in-sample) | 3981 | 0.1188 | 0.1145 | -0.0043 | True |
| T_core | M1.6-1.8 (2012 in-sample) | 4769 | 0.138 | 0.1333 | -0.0046 | True |
| T_core | M1.8-2.0 (2012 in-sample) | 3462 | 0.1459 | 0.1564 | 0.0105 | True |
| T_core | M>=2.0 (2012 in-sample) | 1169 | 0.2164 | 0.2089 | -0.0076 | True |
| T_core | all (2012 in-sample) | 13381 | 0.1412 | 0.1403 | -0.0009 | True |

Production Mc: median **2.60** (T_ext) / **2.76** (T_core), IQR 2.42–3.09. E5's independent prediction model gave Mc(8 MAD) ≈ 2.25–2.29 for this network; ours sits above it because the tiers are 8.25 / 8.75 MAD.

**Caveat — extrapolation.** The referee windows are late-summer ETS windows, so the day terms were calibrated over a narrow range: only 52 % of production family-days fall inside the referee 1–99 % range (column `day_terms_in_calib_range`). Because |c2| ≤ 0.32, a day term 0.3 log-units outside the range moves Mc by ≤ 0.1 magnitude, so the extrapolation is small but it is not zero.

**Cross-check against the independent E4/E5 measurements.** In this parameterisation alpha = (amplitude slope in log10 units per magnitude) / (spread of log10 of the network statistic at fixed magnitude). E4 measured a within-family amplitude slope of 10^(0.38 dM) but noted that the measured within-family RECALL slope is about half of it (alpha_amp ~ 0.2, E5's 'central case'), and E5 measured the spread sigma_log = 0.175 pooled. Those give 0.20 / 0.175 = 1.14 -- against our fitted **1.20**. The sensitivity model therefore reproduces two independently measured quantities it was never shown.

Another caveat: `c1` (the MAD-scale term) comes out POSITIVE on 2012 — on busy days the cc series' own MAD tracks activity, so `u` is partly an activity covariate rather than a pure noise term — and its sign flips in a 2010-only fit. Only the 2012 coefficients are exported, and the out-of-sample test above is what licenses them.

### 1c. False-alarm floor λ(family, tier, epoch, season)

λ = (catalog detections on quiet family-days / quiet observable family-days) × (1 − p_real), with p_real = 0.225 at T_core and 0.057 at T_ext from Track A step 3b (the withheld-station stack test), gamma-Poisson shrunk (k = 10 template-days) toward the epoch-season pool. A quiet family-day is Track C's `quiet_flag`: no PNSN tremor in the family box that day AND ≥ 3 days from any day that has one. Cells with no quiet day at all fall back family-epoch → epoch-season → epoch → all, recorded per row in `lambda_src`.

| tier | epoch_id | season | n_fam | quiet_family_days | quiet_det | rate_pooled | lambda_pooled |
|---|---|---|---|---|---|---|---|
| T_core | E0 | DJF | 115 | 16401.0 | 1807.0 | 0.11018 | 0.08539 |
| T_core | E0 | JJA | 115 | 13219.0 | 613.0 | 0.04637 | 0.03594 |
| T_core | E0 | MAM | 115 | 14357.0 | 1327.0 | 0.09243 | 0.07163 |
| T_core | E0 | SON | 115 | 14490.0 | 903.0 | 0.06232 | 0.0483 |
| T_core | E1 | DJF | 115 | 6202.0 | 78.0 | 0.01258 | 0.00975 |
| T_core | E1 | JJA | 115 | 6548.0 | 312.0 | 0.04765 | 0.03693 |
| T_core | E1 | MAM | 115 | 9651.0 | 449.0 | 0.04652 | 0.03606 |
| T_core | E1 | SON | 115 | 9085.0 | 794.0 | 0.0874 | 0.06773 |
| T_core | E2 | DJF | 115 | 12442.0 | 185.0 | 0.01487 | 0.01152 |
| T_core | E2 | JJA | 115 | 8050.0 | 501.0 | 0.06224 | 0.04823 |
| T_core | E2 | MAM | 115 | 8229.0 | 324.0 | 0.03937 | 0.03051 |
| T_core | E2 | SON | 115 | 9660.0 | 426.0 | 0.0441 | 0.03418 |
| T_core | E3 | DJF | 114 | 21294.0 | 636.0 | 0.02987 | 0.02315 |
| T_core | E3 | JJA | 114 | 30470.0 | 2004.0 | 0.06577 | 0.05097 |
| T_core | E3 | MAM | 114 | 22506.0 | 1010.0 | 0.04488 | 0.03478 |
| T_core | E3 | SON | 114 | 24204.0 | 995.0 | 0.04111 | 0.03186 |
| T_ext | E0 | DJF | 115 | 16401.0 | 4571.0 | 0.2787 | 0.26282 |
| T_ext | E0 | JJA | 115 | 13219.0 | 1741.0 | 0.1317 | 0.1242 |
| T_ext | E0 | MAM | 115 | 14357.0 | 3236.0 | 0.2254 | 0.21255 |
| T_ext | E0 | SON | 115 | 14490.0 | 2448.0 | 0.16894 | 0.15931 |
| T_ext | E1 | DJF | 115 | 6202.0 | 623.0 | 0.10045 | 0.09473 |
| T_ext | E1 | JJA | 115 | 6548.0 | 2545.0 | 0.38867 | 0.36651 |
| T_ext | E1 | MAM | 115 | 9651.0 | 3285.0 | 0.34038 | 0.32098 |
| T_ext | E1 | SON | 115 | 9085.0 | 3511.0 | 0.38646 | 0.36443 |
| T_ext | E2 | DJF | 115 | 12442.0 | 1289.0 | 0.1036 | 0.0977 |
| T_ext | E2 | JJA | 115 | 8050.0 | 2632.0 | 0.32696 | 0.30832 |
| T_ext | E2 | MAM | 115 | 8229.0 | 1780.0 | 0.21631 | 0.20398 |
| T_ext | E2 | SON | 115 | 9660.0 | 2354.0 | 0.24369 | 0.2298 |
| T_ext | E3 | DJF | 114 | 21294.0 | 3699.0 | 0.17371 | 0.16381 |
| T_ext | E3 | JJA | 114 | 30470.0 | 12668.0 | 0.41575 | 0.39206 |
| T_ext | E3 | MAM | 114 | 22506.0 | 6823.0 | 0.30316 | 0.28588 |
| T_ext | E3 | SON | 114 | 24204.0 | 6777.0 | 0.28 | 0.26404 |

**Seasonal-stability bar** (Merlin: winter quiet-day FA within 2x summer), evaluated per epoch — and the answer is not a single number:

| tier | epoch_id | DJF_rate | JJA_rate | DJF_over_JJA | bar_met |
|---|---|---|---|---|---|
| T_core | E0 | 0.1102 | 0.0464 | 2.3759 | False |
| T_core | E1 | 0.0126 | 0.0476 | 0.2639 | False |
| T_core | E2 | 0.0149 | 0.0622 | 0.2389 | False |
| T_core | E3 | 0.0299 | 0.0658 | 0.4541 | False |
| T_ext | E0 | 0.2787 | 0.1317 | 2.1161 | False |
| T_ext | E1 | 0.1005 | 0.3887 | 0.2585 | False |
| T_ext | E2 | 0.1036 | 0.327 | 0.3169 | False |
| T_ext | E3 | 0.1737 | 0.4158 | 0.4178 | False |

**8 of 8 epoch-tier cells miss the ≤ 2x band, and they miss it in BOTH directions**: winter is ~2.4x summer in E0 but summer is ~2–4x winter in E1, E2 and E3. The bar is reported as NOT met and nothing is tuned to it; what the full record shows is that the quiet-day floor is seasonally UNSTABLE rather than seasonally BIASED, which is why λ is carried per (family, tier, epoch, season) on every event rather than as one number.


### 1d. What drives the seasonal instability of the floor

The ratios are reported as measured and nothing is tuned to them. Two candidate causes were separated (`scripts/d1d_season_diagnostic.py`, `scripts/d1d_season_stack_test.py`; figure `figures/d1d_season_diagnostic.png`): a NOISE / MAD-adaptivity effect, or genuinely more inter-ETS LFE activity in the busier season.

**Scope note.** The withheld-station test in (iv) below was designed and run when only 2014–2015 was on disk, where winter was the busy season; its six months are all in epoch **E0** (DJF 2015-02 + 2014-02 + 2014-12 vs JJA 2014-08 + 2015-08 + 2014-07, matched quiet-day counts). It therefore answers the question for E0's WINTER excess specifically. The full record then showed the sign of the seasonal imbalance flips in E1–E3, which is itself the main finding of (ii) and does not depend on (iv).

**(i) The detector's threshold does not adapt to the ambient level.** The score is in units of each family-day's own network-MAD, and that MAD is formed from an amplitude-NORMALISED cross-correlation, so it barely moves when the ground gets louder. Quiet-day detections per template-day against threshold:

| epoch_id | season | 8.0 | 8.25 | 8.5 | 8.75 | 9.0 | 9.5 | 10.0 | 11.0 |
|---|---|---|---|---|---|---|---|---|---|
| E0 | DJF | 0.5483 | 0.3336 | 0.2093 | 0.1333 | 0.0858 | 0.0365 | 0.0175 | 0.0055 |
| E0 | JJA | 0.2773 | 0.1637 | 0.1002 | 0.0624 | 0.0396 | 0.018 | 0.0093 | 0.0026 |
| E0 | MAM | 0.4499 | 0.2805 | 0.1826 | 0.1262 | 0.0934 | 0.054 | 0.0359 | 0.0161 |
| E0 | SON | 0.3609 | 0.2159 | 0.136 | 0.0892 | 0.0606 | 0.03 | 0.0173 | 0.0064 |
| E1 | DJF | 0.1482 | 0.0921 | 0.0606 | 0.041 | 0.0282 | 0.0158 | 0.0098 | 0.0047 |
| E1 | JJA | 0.4959 | 0.2946 | 0.1762 | 0.1156 | 0.0776 | 0.0397 | 0.0211 | 0.0104 |
| E1 | MAM | 0.391 | 0.2312 | 0.1394 | 0.0901 | 0.0604 | 0.0271 | 0.0135 | 0.0052 |
| E1 | SON | 0.6404 | 0.4316 | 0.3024 | 0.2184 | 0.1586 | 0.0935 | 0.0566 | 0.0171 |
| E2 | DJF | 0.1573 | 0.1011 | 0.0664 | 0.0473 | 0.0339 | 0.0202 | 0.0129 | 0.0052 |
| E2 | JJA | 0.4891 | 0.3294 | 0.2289 | 0.1667 | 0.1267 | 0.0729 | 0.044 | 0.0159 |
| E2 | MAM | 0.3064 | 0.2011 | 0.1382 | 0.1 | 0.0689 | 0.0361 | 0.0208 | 0.0062 |
| E2 | SON | 0.3573 | 0.2394 | 0.1622 | 0.1148 | 0.0835 | 0.0484 | 0.0297 | 0.0124 |
| E3 | DJF | 0.2931 | 0.2046 | 0.1477 | 0.1063 | 0.079 | 0.0483 | 0.0304 | 0.0116 |
| E3 | JJA | 0.6366 | 0.4358 | 0.3082 | 0.2223 | 0.1657 | 0.0969 | 0.0609 | 0.0268 |
| E3 | MAM | 0.4677 | 0.3152 | 0.2214 | 0.1616 | 0.1202 | 0.0708 | 0.044 | 0.0172 |
| E3 | SON | 0.4481 | 0.306 | 0.2157 | 0.1535 | 0.1124 | 0.0625 | 0.0358 | 0.0123 |

Tail decay (dex per MAD over 8.25–10; more negative = a lighter tail):

| epoch_id | season | tail_slope_dex_per_MAD | n_quiet_days |
|---|---|---|---|
| E0 | DJF | -0.734 | 16401 |
| E0 | JJA | -0.712 | 13219 |
| E0 | MAM | -0.504 | 14357 |
| E0 | SON | -0.626 | 14490 |
| E1 | DJF | -0.553 | 6202 |
| E1 | JJA | -0.644 | 6548 |
| E1 | MAM | -0.7 | 9651 |
| E1 | SON | -0.5 | 9085 |
| E2 | DJF | -0.503 | 12442 |
| E2 | JJA | -0.494 | 8050 |
| E2 | MAM | -0.565 | 8229 |
| E2 | SON | -0.512 | 9660 |
| E3 | DJF | -0.471 | 21294 |
| E3 | JJA | -0.486 | 30470 |
| E3 | MAM | -0.485 | 22506 |
| E3 | SON | -0.53 | 24204 |

The winter tail is **not** heavier — its decay (-0.734 dex/MAD) is as steep as or steeper than summer's (-0.712). The winter excess is a near-uniform elevation of the whole tail, not an extra population of large excursions.

**(ii) The floor is DECOUPLED from the ambient level — the strongest result of this section.** Winter is the loudest season in every epoch, yet it carries the highest quiet-day rate in only ONE of them, and the family-day MAD scale barely moves at all:

| epoch_id | season | mad_scale_median | n_sites_median | n_quiet_days | e8_rel_median | e8_rel_p90 | n_station_days | u_median | v_median | Mc_core_median |
|---|---|---|---|---|---|---|---|---|---|---|
| E0 | DJF | 0.0166 | 6.0 | 16401 | 1.6113 | 84.8232 | 5174 | 0.0187 | 0.0215 | 2.6843 |
| E0 | JJA | 0.0164 | 6.0 | 13219 | 1.1127 | 72.0773 | 4423 | 0.0069 | -0.0096 | 2.6983 |
| E0 | MAM | 0.0168 | 6.0 | 14357 | 1.1885 | 58.1461 | 4451 | 0.0179 | 0.0015 | 2.6853 |
| E0 | SON | 0.0162 | 6.0 | 14490 | 1.153 | 73.3187 | 4245 | 0.0055 | -0.0138 | 2.7061 |
| E1 | DJF | 0.0172 | 6.0 | 6202 | 1.192 | 2.4425 | 1848 | 0.0285 | 0.0857 | 2.7444 |
| E1 | JJA | 0.0161 | 6.0 | 6548 | 0.9738 | 2.9934 | 2744 | 0.0027 | -0.0089 | 2.7643 |
| E1 | MAM | 0.0168 | 6.0 | 9651 | 1.1505 | 4.8182 | 2821 | 0.0162 | 0.0138 | 2.7549 |
| E1 | SON | 0.0161 | 6.0 | 9085 | 0.8545 | 2.0798 | 2883 | 0.0038 | 0.0796 | 2.7849 |
| E2 | DJF | 0.0173 | 6.0 | 12442 | 1.0953 | 2.3838 | 3537 | 0.0255 | 0.0679 | 2.7571 |
| E2 | JJA | 0.0167 | 6.0 | 8050 | 0.839 | 1.968 | 2632 | 0.0132 | -0.0096 | 2.7613 |
| E2 | MAM | 0.0168 | 6.0 | 8229 | 0.841 | 1.7739 | 3045 | 0.0121 | 0.0077 | 2.7639 |
| E2 | SON | 0.0167 | 6.0 | 9660 | 0.7958 | 2.2933 | 2655 | 0.0125 | -0.0309 | 2.7568 |
| E3 | DJF | 0.0169 | 5.0 | 21294 | 0.9625 | 2.0258 | 6093 | 0.0296 | 0.0365 | 2.7609 |
| E3 | JJA | 0.0163 | 5.0 | 30470 | 0.8925 | 2.1896 | 7499 | 0.0104 | 0.0638 | 2.8427 |
| E3 | MAM | 0.0167 | 5.0 | 22506 | 0.9594 | 2.259 | 7091 | 0.0211 | 0.0435 | 2.8019 |
| E3 | SON | 0.0162 | 5.0 | 24204 | 0.866 | 2.4055 | 6390 | 0.0092 | 0.0506 | 2.8064 |

Across the four epochs the DJF/JJA ratio of relative daily noise energy is E0 1.45, E1 1.22, E2 1.31, E3 1.08 — winter is louder every time — while the MAD scale stays within a few per cent of 0.0165 in every epoch and season. A detector scored in units of an amplitude-NORMALISED cross-correlation's own MAD is therefore nearly blind to the ambient level: **a louder season does not raise the bar, and a louder season is not what raises the floor.** What raises the floor is the population of COHERENT transients in that epoch-season, which the daily median energy does not measure — and which has no fixed seasonal sign.

(For the first epoch specifically, E0:) relative daily noise energy DJF 1.61 vs JJA 1.11 (a factor 1.45 in energy, 1.20 in amplitude), against a MAD scale that moves by only +1.7 %. So a louder winter does NOT raise the bar — the false-alarm floor is free to rise with the number of coherent winter transients. This is the mechanism; it is a property of the frozen MAD-relative operating point, not a fault in the month.

**(iii) Within a season, the noisiest days are not the ones that fire.** Slope of the quiet-day detection count against the daily noise term `v`:

| epoch_id | season | tier | n_family_days | rate | slope_per_dex_v | corr |
|---|---|---|---|---|---|---|
| E0 | DJF | T_ext | 16401 | 0.3336 | -1.6778 | -0.1383 |
| E0 | DJF | T_core | 16401 | 0.1333 | -0.7526 | -0.1202 |
| E0 | JJA | T_ext | 13219 | 0.1637 | 0.4258 | 0.0364 |
| E0 | JJA | T_core | 13219 | 0.0624 | 0.161 | 0.0274 |
| E0 | MAM | T_ext | 14357 | 0.2805 | -0.9801 | -0.0296 |
| E0 | MAM | T_core | 14357 | 0.1262 | -0.4977 | -0.0212 |
| E0 | SON | T_ext | 14490 | 0.2159 | -0.2786 | -0.0291 |
| E0 | SON | T_core | 14490 | 0.0892 | -0.2054 | -0.0421 |
| E1 | DJF | T_ext | 6202 | 0.0921 | -0.3948 | -0.1144 |
| E1 | DJF | T_core | 6202 | 0.041 | -0.1915 | -0.0949 |
| E1 | JJA | T_ext | 6548 | 0.2946 | -1.7862 | -0.1197 |
| E1 | JJA | T_core | 6548 | 0.1156 | -0.9616 | -0.1328 |
| E1 | MAM | T_ext | 9651 | 0.2312 | -0.7394 | -0.0378 |
| E1 | MAM | T_core | 9651 | 0.0901 | -0.3168 | -0.0321 |
| E1 | SON | T_ext | 9085 | 0.4316 | -0.2135 | -0.0132 |
| E1 | SON | T_core | 9085 | 0.2184 | -0.1498 | -0.0141 |
| E2 | DJF | T_ext | 12442 | 0.1011 | -0.1461 | -0.0696 |
| E2 | DJF | T_core | 12442 | 0.0473 | -0.0802 | -0.0631 |
| E2 | JJA | T_ext | 8050 | 0.3294 | -0.3685 | -0.0365 |
| E2 | JJA | T_core | 8050 | 0.1667 | -0.1855 | -0.0337 |
| E2 | MAM | T_ext | 8229 | 0.2011 | -0.1725 | -0.0192 |
| E2 | MAM | T_core | 8229 | 0.1 | -0.144 | -0.0286 |
| E2 | SON | T_ext | 9660 | 0.2394 | -0.3226 | -0.0374 |
| E2 | SON | T_core | 9660 | 0.1148 | -0.2019 | -0.0418 |
| E3 | DJF | T_ext | 21294 | 0.2046 | -0.2039 | -0.021 |
| E3 | DJF | T_core | 21294 | 0.1063 | -0.1323 | -0.0199 |
| E3 | JJA | T_ext | 30470 | 0.4358 | -0.1819 | -0.031 |
| E3 | JJA | T_core | 30470 | 0.2223 | -0.106 | -0.0278 |
| E3 | MAM | T_ext | 22506 | 0.3152 | -0.2329 | -0.0451 |
| E3 | MAM | T_core | 22506 | 0.1616 | -0.128 | -0.0425 |
| E3 | SON | T_ext | 24204 | 0.306 | 0.0615 | 0.0174 |
| E3 | SON | T_core | 24204 | 0.1535 | 0.0354 | 0.0171 |

Inside winter the slope is NEGATIVE — on the loudest days the correlation coefficient itself is depressed and fewer detections clear the threshold. The seasonal excess is therefore not carried by the loudest days; it is carried by the winter population of coherent transients (storm microseism bursts and uncatalogued micro-earthquakes) that the daily median energy does not measure.

**(iv) The decisive test — is the winter excess real LFE activity?** Track A step 3b's withheld-station test, repeated per season on production months with CN.VGZ, CN.SNB and PB.B926 excluded from the network sum (`scripts/d1d_run_withheld.sh`, six GPU passes), the frozen HF veto applied, quiet family-days taken from Track C's own `quiet_flag`, and months chosen so the two seasons have MATCHED quiet-day counts (DJF 2015-02 + 2014-02 + 2014-12 vs JJA 2014-08 + 2015-08 + 2014-07).

Production's quiet-day rate is ~4x lower than 2012's, so no single family has enough quiet-day detections for step 3b's coherent per-family stack (2–7 here against 20+ there). The primary statistic is therefore a POOLED INCOHERENT one: for every quiet-day detection, r = E(S ± 2 s) / E(pre-arrival 4 s) at the withheld station, each window normalised by its own noise, averaged in log10 over thousands of detections. **Two nulls**, and the second is the one that matters:

* **time null** — matched random times on the same quiet days, same offsets: *is there anything there at all?*
* **moveout-shuffled null** — the SAME detection times with ANOTHER family's S offset at that station (offsets differing by > 5 s): *is the energy at the FAMILY-SPECIFIC predicted arrival?* A broadband regional transient that made the network fire raises both arms; a real LFE raises only the real one. Step 3b did not have this control; at production's purity it is essential.

| season | tier | n_quiet_family_days | n_det | excess_db | excess_db_se | z | p_null | mo_mean_log10r | excess_vs_mo_db | p_mo |
|---|---|---|---|---|---|---|---|---|---|---|
| DJF | T_core | 6293 | 1940 | 0.9178 | 0.1482 | 2.8185 | 0.0488 | 0.0236 | 0.5864 | 0.0476 |
| DJF | T_ext | 6293 | 4559 | 0.4489 | 0.0967 | 2.131 | 0.0488 | 0.028 | 0.1727 | 0.0476 |
| JJA | T_core | 6325 | 1808 | 1.5922 | 0.2628 | 13.5628 | 0.0244 | 0.0854 | 0.737 | 0.0476 |
| JJA | T_ext | 6325 | 3659 | 1.0326 | 0.1373 | 14.672 | 0.0244 | 0.0521 | 0.5066 | 0.0476 |

(seasonal contrast not rendered: float() argument must be a string or a real number, not 'Series')


## 2. Catalog v1 (Merlin Step 7)

`results/trackD/lfe_extension_catalog_v1_1.parquet` (+ `.csv`): **{meta['n_events']:,} rows** (74,330 `cat_core`, 244,375 `cat_ext`; the tiers overlap), 115 families, 2014-01-01 .. 2026-09-15. 0 duplicates removed; month ownership is half-open (a detection belongs to the month of its own UTC date) and every month's count agrees with `months_done.csv` (0 mismatches).

**The tiers are not nested.** 3694 rows are `cat_core` but not `cat_ext`: T_core > T_ext, so every core detection clears the ext threshold, but the frozen cross-family assignment is evaluated separately on each tier's own set and a detection can lose to a station-sharing neighbour at T_ext that is itself below T_core. Both flags are stored; use the one the analysis needs, never their union as a single tier.

**TIME BASIS.** `time` is the origin-axis time the network sum was formed on = **Bostock's database clock for that family plus the family's fitted offset** (`e5r_family_offsets_2012.csv`, 3.25–21.25 s, median 12.25 s). It is **NOT an origin time**: the datum is Bostock's and the offset was fitted, so only DIFFERENCES carry physical meaning — within a family exactly, and between families to the accuracy of the offsets.

Events per year:

| year | cat_core | cat_ext | families |
|---|---|---|---|
| 2014 | 4511 | 9057 | 115 |
| 2015 | 7542 | 15195 | 115 |
| 2016 | 6382 | 14252 | 115 |
| 2017 | 4411 | 17874 | 115 |
| 2018 | 6772 | 21797 | 115 |
| 2019 | 6347 | 19539 | 115 |
| 2020 | 4777 | 15077 | 115 |
| 2021 | 5988 | 21749 | 115 |
| 2022 | 5963 | 24478 | 114 |
| 2023 | 5491 | 23323 | 114 |
| 2024 | 7354 | 27834 | 114 |
| 2025 | 6318 | 23553 | 114 |
| 2026 | 2474 | 10647 | 114 |

Events and floor by epoch:

| epoch_id | cat_core | cat_ext | families | median_lambda_core | median_Mc_core |
|---|---|---|---|---|---|
| E0 | 18958 | 39841 | 115 | 0.0972 | 2.5597 |
| E1 | 10660 | 38334 | 115 | 0.0286 | 2.6667 |
| E2 | 12307 | 38093 | 115 | 0.0284 | 2.6232 |
| E3 | 32405 | 128107 | 114 | 0.0499 | 2.4714 |

How each epoch's floor was estimated (number of (family, season) cells per fallback level, T_core) — `family-epoch-season` is the direct measurement, the rest are fallbacks for cells with no quiet observable family-day:

| epoch_id | epoch-season | family-epoch-season |
|---|---|---|
| E0 | 4 | 460 |
| E1 | 4 | 460 |
| E2 | 4 | 460 |
| E3 | 8 | 456 |

`family_table_v1.csv`: 130 families, **14 template-less** (no usable 2012 template at any station — they can never be detected and are carried explicitly), 115 with at least one `cat_core` event; 10 carry our own accepted hypocentre (σ_h ≈ 2.0 km from the matched synthetic floor), 88 carry a reported-only one.

`family_day_table_v1.parquet`: one row per (family, day) with `observable`, `n_sites`, `Mc_*`, `Pdet_*`, `lambda_*` and the day terms.

## 3. Science QC — the Thomas harness per SSE (Merlin Step 8)

**SSE windows.** PNSN tremor in the family box, ≥ 50 epicentres/day, clustered with gaps ≤ 20 d, span ≥ 10 d, padded ±3 d → **31 windows** 2014–2026, 15 of them with ≥ 10 active days ('major'). A second set with a 5-day clustering gap (`sse_windows_sub.csv`) is run as a reported sensitivity. Windows, with the harness result for each, are `results/trackD/thomas_per_sse.csv`; the figures are `figures/d3_sse_<id>_<tier>.png` and `figures/d3_thomas_summary_<tier>.png`.

Three legs are run per window: `cat_core`, `cat_ext`, and `cat_core_hp` = `cat_core` restricted to ≥ 5 contributing SITES — a detector-internal purity leg, reported and never barred, because the harness's leading edge is the 10th percentile of event position and at ~50 % precision the background can pin it (the harness then reports 'no departure').

### 3a. The two pre-registered bars, on the 2012 referee

| bar | tier | ours | bostock | rel_diff | tol | met |
|---|---|---|---|---|---|---|
| main-front first-burst speed | cat_core | -5.2754 | -4.6949 | 0.1237 | 0.2 | True |
| re-rupture interval median | cat_core | 1.0513 | 0.7039 | 0.4936 | 0.05 | False |
| main-front first-burst speed | cat_ext | -5.2283 | -4.6949 | 0.1136 | 0.2 | True |
| re-rupture interval median | cat_ext | 0.8663 | 0.7039 | 0.2307 | 0.05 | False |

**The re-rupture bar failure is the decimation bias, not a catalog fault.** Bostock's 2012 catalog thinned per family at OUR measured recall plus OUR measured quiet-day floor (`d3_thinned_reference.py`, 20 seeds) is the fair reference:

| tier | metric | ours | thinned_mean | thinned_p2_5 | thinned_p97_5 | bostock_full | inside_thinned_95 |
|---|---|---|---|---|---|---|---|
| T_core | front_firstburst_kmday |  | -4.9603 | -5.2733 | -4.6715 | -4.6949 | False |
| T_core | reactivation_interval_d_p50 |  | 0.8567 | 0.674 | 1.0585 | 0.7039 | False |
| T_core | n_secondary |  | 20.55 | 15.475 | 28.575 | 20.0 | False |
| T_core | n_secondary_null_circ |  | 13.88 | 12.6238 | 15.5888 | 2.25 | False |
| T_core | frac_fam_reactivated |  | 0.1571 | 0.1259 | 0.1808 | 0.4923 | False |
| T_ext | front_firstburst_kmday |  | -4.9895 | -5.3148 | -4.8022 | -4.6949 | False |
| T_ext | reactivation_interval_d_p50 |  | 0.915 | 0.6931 | 1.1389 | 0.7039 | False |
| T_ext | n_secondary |  | 12.75 | 7.0 | 19.0 | 20.0 | False |
| T_ext | n_secondary_null_circ |  | 9.68 | 8.2538 | 10.7362 | 2.25 | False |
| T_ext | frac_fam_reactivated |  | 0.1901 | 0.127 | 0.2468 | 0.4923 | False |

### 3b. How often the harness's estimators actually work at this completeness

| tier | n_windows | first_burst_defined | first_burst_plausible | lead_track_defined | secondary_above_null_p95 | n_with_reactivation | reactivation_median_over_sses_d |
|---|---|---|---|---|---|---|---|
| cat_core | 31 | 7 | 3 | 15 | 16 | 21 | 1.02 |
| cat_core_hp | 31 | 5 | 3 | 8 | 11 | 16 | 0.989 |
| cat_ext | 31 | 13 | 6 | 23 | 4 | 29 | 0.935 |

**This is the honest limit of the science QC.** The main-front first-burst estimator returns a physically plausible southward speed (−3 to −15 km/day) in only **3 of 31** windows at `cat_core`, 4 at `cat_core_hp` and 6 at `cat_ext`; in most windows the harness reports 'no departure' or a degenerate propagation window. The cause is visible in the space-time figures: its leading edge is the 10th percentile of event position, and at ~50 % precision the false-alarm background — which is spread along the whole strike — pins that percentile before the front arrives. Where the estimator does work the speeds are cat_core: -7.84;-5.82;-6.05, cat_core_hp: -8.08;-6.71;-6.25, cat_ext: -8.55;-3.28;-4.88;-8.58;-7.21;-4.56 km/day, bracketing Bostock's 2012 value of −4.70 and Thomas's 7.3 km/day reference. The re-activation interval, which needs no front, is available in most windows and its median across SSEs (~0.9–1.0 d) sits where the per-family-thinned reference says it should.

A front estimator matched to this catalog's purity — a first-activation fit restricted to families whose burst exceeds their own λ — was named here as the obvious next step. It was subsequently approved and built; see §3b-2.

### 3b-2. A purity-matched main-front estimator (POST-REGISTRATION, approved 2026-09-17)

**This section is not part of the pre-registered plan.** The estimator below was designed after §3b showed the harness's front estimator failing on this catalog; it was approved by the coordinator on 2026-09-17, designed on the 2012 referee ALONE, frozen, and only then run on the extension years. Code: `scripts/d3b2_front_v2.py`; figure `figures/d3b2_front_v2_summary.png`; results `thomas_front_v2.csv`.

**Definition.** Two changes to the harness's first-burst front, both aimed at the diagnosed failure and nothing else:

1. **Purity gate.** A family activates on the first day whose event count for that family is inconsistent with that family's OWN measured false-alarm floor: `n_f(d) ≥ N_MIN` AND `P(N ≥ n_f(d) | Poisson(λ_f)) < α / n_days` (Bonferroni over the window's days), with λ from §1c for that (family, tier, epoch, season). The activation time is the start of the first burst on that day (harness definition: gap < 600 s, ≥ 4 events).

2. **Window from the activations.** The propagation window is the 5th–95th percentile of the activation times, not the 10th-percentile leading-edge track — so background families cannot set it.

Everything else is the harness's (Theil-Sen, ≥ 5 families, 300-draw family bootstrap). One along-strike frame, computed once from all 130 Bostock hypocentres, is used for every window so speeds are comparable between SSEs.

**Free parameters and how they were fixed.** `N_MIN = 4` and `MIN_FAM = 5` are inherited from the harness; the 5–95 % window is fixed a priori. The ONLY fitted parameter is α, chosen from the pre-declared grid {0.05, 0.01, 0.001} as the first value putting ALL FOUR 2012 referee catalogs within 20 % of Bostock's 2012 first-burst speed (−4.695 km/day):

| alpha | leg | n_events | n_fam_activated | n_fit | win_days | speed_kmday | ci_lo | ci_hi | rel_diff_vs_bostock2012 | within_20pct |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.05 | bostock2012 | 28149.0 | 129.0 | 115 | 13.401 | -5.45 | -5.882 | -5.013 | 0.161 | True |
| 0.05 | ours2012_core | 7329.0 | 106.0 | 94 | 13.596 | -5.033 | -5.357 | -4.671 | 0.072 | True |
| 0.05 | ours2012_ext | 11405.0 | 110.0 | 98 | 13.73 | -4.615 | -5.024 | -4.087 | 0.017 | True |
| 0.05 | bostock_thinned_to_our_completeness |  |  | 20 |  | -4.844 | -4.933 | -4.73 | 0.032 | True |
| 0.01 | bostock2012 | 28149.0 | 129.0 | 115 | 13.401 | -5.45 | -5.898 | -4.973 | 0.161 | True |
| 0.01 | ours2012_core | 7329.0 | 105.0 | 93 | 13.42 | -5.115 | -5.424 | -4.765 | 0.09 | True |
| 0.01 | ours2012_ext | 11405.0 | 110.0 | 98 | 13.726 | -4.787 | -5.186 | -4.338 | 0.02 | True |
| 0.01 | bostock_thinned_to_our_completeness |  |  | 20 |  | -4.798 | -4.909 | -4.71 | 0.022 | True |
| 0.001 | bostock2012 | 28149.0 | 129.0 | 115 | 13.401 | -5.45 | -5.89 | -5.031 | 0.161 | True |
| 0.001 | ours2012_core | 7329.0 | 105.0 | 93 | 13.42 | -5.125 | -5.434 | -4.755 | 0.092 | True |
| 0.001 | ours2012_ext | 11405.0 | 110.0 | 98 | 13.566 | -4.922 | -5.265 | -4.522 | 0.048 | True |
| 0.001 | bostock_thinned_to_our_completeness |  |  | 20 |  | -4.798 | -4.907 | -4.71 | 0.022 | True |

All three α values pass on all four legs — the estimator is insensitive to it (the four referee speeds move by ≤ 0.3 km/day across the whole grid) — so the rule selects **α = 0.05** and it is frozen there (`thomas_front_v2_design.json`). The extension years were not looked at until that point.

**Nulls.** The harness's circular and day-shuffle nulls are run unchanged, and a third is added because one of them is degenerate here: the day-shuffle permutes family labels among the SAME day's events, so every family's DAILY count is exactly preserved and a day-resolution activation statistic is nearly invariant under it. The evidence is in the table below — `null_dayshuf_mean` reproduces the observed speed to a few per cent in almost every window. It is reported for continuity and is NOT a significance test of this estimator. The added **position-shuffle** null keeps every activation time and permutes the along-strike coordinate across families, which asks directly whether activation ORDER relates to POSITION.

**Plausibility flag** (pre-declared): ≥ 5 families fitted, a bootstrap CI excluding zero, and |speed| ≤ 20 km/day (4× Thomas's 7.3 km/day reference).

| tier | n_windows | fitted | plausible | plausible_and_xshuf_sig | n_SE | n_NW | median_abs_speed | v1_plausible_for_comparison |
|---|---|---|---|---|---|---|---|---|
| cat_core | 31 | 27 | 13 | 12 | 6 | 6 | 5.554 | 3 |
| cat_ext | 31 | 30 | 11 | 8 | 5 | 3 | 4.39 | 6 |

**Yield: 27 of 31 windows fitted and 11 significant at `cat_core`, against 3 of 31 for the pre-registered estimator** — the diagnosis in §3b was right and the fix works. Full per-SSE table (speed, CI, direction, propagation window, nulls, plausibility) is `results/trackD/thomas_front_v2.csv`; per-family activations are `thomas_front_v2_activations.csv`.

**cat_core**:

| sse_id | t0 | major | n_events | n_fam_activated | n_fit | win_start | win_days | speed_kmday | ci_lo | ci_hi | direction | rho | null_circ_p | null_xshuf_p | null_dayshuf_mean | plausible |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SSE201406 | 2014-06-10 | False | 537 | 16 | 14.0 | 2014-06-14 11:21:00 | 20.589 | 1.065 | -5.295 | 6.54 | NW | 0.218 | 0.439 | 0.317 | 1.331 | False |
| SSE201411 | 2014-11-14 | True | 1094 | 37 | 33.0 | 2014-11-22 14:14:00 | 20.653 | -0.983 | -2.712 | 3.211 | SE (Thomas's southward) | -0.293 | 0.244 | 0.268 | -1.041 | False |
| SSE201507 | 2015-07-05 | False | 399 | 7 | 5.0 | 2015-07-22 10:06:00 | 5.599 | 0.229 | -0.859 | 9.103 | NW | 0.3 | 0.368 | 1.0 | 0.349 | False |
| SSE201510 | 2015-10-05 | False | 854 | 11 | 9.0 | 2015-10-06 17:22:00 | 11.883 | -8.018 | -50.256 | 5.677 | SE (Thomas's southward) | -0.583 | 0.049 | 0.244 | -3.976 | False |
| SSE201512 | 2015-12-16 | True | 2947 | 75 | 67.0 | 2015-12-19 00:50:00 | 13.425 | -8.551 | -9.893 | -7.553 | SE (Thomas's southward) | -0.853 | 0.024 | 0.024 | -8.136 | True |
| SSE201608 | 2016-08-03 | False | 1464 | 15 | 13.0 | 2016-08-04 00:05:00 | 28.925 | -2.229 | -3.238 | 0.79 | SE (Thomas's southward) | -0.346 | 0.024 | 0.024 | -2.372 | False |
| SSE201703 | 2017-03-15 | True | 1294 | 56 | 50.0 | 2017-03-22 02:22:00 | 8.852 | 8.477 | 7.114 | 9.605 | NW | 0.852 | 0.024 | 0.024 | 8.407 | True |
| SSE201709 | 2017-09-17 | True | 1048 | 13 | 11.0 | 2017-09-21 03:59:00 | 24.788 | 5.806 | -0.056 | 15.182 | NW | 0.382 | 0.024 | 0.024 | 5.417 | False |
| SSE201805 | 2018-04-30 | False | 421 | 10 | 8.0 | 2018-05-03 01:59:00 | 14.082 | -5.229 | -5.913 | -0.28 | SE (Thomas's southward) | -0.81 | 0.098 | 0.024 | -4.556 | True |
| SSE201806 | 2018-06-17 | True | 3618 | 74 | 66.0 | 2018-06-26 20:19:00 | 12.867 | -5.681 | -6.258 | -5.006 | SE (Thomas's southward) | -0.902 | 0.024 | 0.024 | -5.896 | True |
| SSE201904 | 2019-04-19 | False | 531 | 9 | 7.0 | 2019-04-22 07:29:00 | 13.137 | -0.018 | -13.464 | 5.331 | SE (Thomas's southward) | -0.179 | 0.878 | 1.0 | -2.457 | False |
| SSE201908 | 2019-08-08 | True | 2227 | 44 | 38.0 | 2019-08-12 00:05:00 | 11.318 | 4.39 | 3.529 | 4.986 | NW | 0.763 | 0.024 | 0.024 | 4.359 | True |
| SSE201911 | 2019-11-05 | True | 1903 | 49 | 43.0 | 2019-11-08 19:23:00 | 16.392 | -2.396 | -5.021 | -1.513 | SE (Thomas's southward) | -0.844 | 0.024 | 0.024 | -2.122 | True |
| SSE202010 | 2020-10-12 | True | 1513 | 37 | 33.0 | 2020-11-02 10:40:00 | 10.889 | 7.323 | 6.148 | 8.935 | NW | 0.926 | 0.024 | 0.024 | 7.153 | True |
| SSE202101 | 2021-01-27 | True | 993 | 47 | 41.0 | 2021-02-03 06:28:00 | 7.175 | -6.703 | -7.494 | -5.229 | SE (Thomas's southward) | -0.916 | 0.024 | 0.024 | -5.89 | True |
| SSE202106 | 2021-06-23 | False | 207 | 7 | 5.0 | 2021-06-28 04:42:00 | 18.173 | -7.836 | -38.58 | -5.31 | SE (Thomas's southward) | -0.9 | 0.22 | 0.22 | -7.017 | True |
| SSE202110 | 2021-09-30 | True | 2010 | 45 | 39.0 | 2021-10-04 03:47:00 | 24.391 | -0.767 | -1.768 | 0.865 | SE (Thomas's southward) | -0.207 | 0.073 | 0.195 | -0.372 | False |
| SSE202203 | 2022-03-14 | True | 1600 | 52 | 46.0 | 2022-03-23 17:45:00 | 13.593 | -4.235 | -5.413 | -2.941 | SE (Thomas's southward) | -0.509 | 0.024 | 0.024 | -3.923 | True |
| SSE202207 | 2022-07-14 | True | 1721 | 23 | 19.0 | 2022-07-17 12:02:00 | 24.936 | 0.061 | -2.196 | 4.894 | NW | -0.135 | 0.829 | 0.902 | 0.28 | False |
| SSE202212 | 2022-12-01 | False | 526 | 7 | 5.0 | 2022-12-08 03:55:00 | 7.95 | -2.64 | -21.897 | 0.865 | SE (Thomas's southward) | -0.7 | 0.122 | 0.902 | -3.388 | False |
| SSE202304 | 2023-04-02 | True | 2640 | 69 | 61.0 | 2023-04-05 09:55:00 | 15.422 | 5.602 | 4.542 | 6.805 | NW | 0.806 | 0.024 | 0.024 | 5.542 | True |
| SSE202306 | 2023-06-18 | False | 522 | 9 | 7.0 | 2023-06-20 00:54:00 | 17.403 | 0.245 | -27.249 | 3.943 | NW | 0.071 | 0.561 | 0.927 | 1.293 | False |
| SSE202308 | 2023-08-13 | False | 919 | 13 | 11.0 | 2023-08-30 06:53:00 | 28.764 | -2.149 | -3.785 | 0.296 | SE (Thomas's southward) | -0.6 | 0.049 | 0.024 | -1.859 | False |
| SSE202401 | 2024-01-26 | False | 504 | 10 | 8.0 | 2024-01-27 09:57:00 | 5.146 | 24.664 | 6.358 | 36.216 | NW | 0.738 | 0.024 | 0.024 | 22.338 | False |
| SSE202407 | 2024-07-16 | True | 4742 | 88 | 78.0 | 2024-07-19 02:20:00 | 40.249 | 5.507 | 4.248 | 6.129 | NW | 0.737 | 0.024 | 0.024 | 5.38 | True |
| SSE202508 | 2025-08-12 | True | 4401 | 75 | 67.0 | 2025-08-18 03:23:00 | 54.495 | 1.303 | 0.224 | 5.046 | NW | 0.247 | 0.024 | 0.024 | 1.253 | True |
| SSE202604 | 2026-04-03 | False | 1246 | 21 | 19.0 | 2026-04-06 22:28:00 | 50.447 | -1.395 | -1.806 | 0.361 | SE (Thomas's southward) | -0.504 | 0.049 | 0.024 | -1.085 | False |

**cat_ext**:

| sse_id | t0 | major | n_events | n_fam_activated | n_fit | win_start | win_days | speed_kmday | ci_lo | ci_hi | direction | rho | null_circ_p | null_xshuf_p | null_dayshuf_mean | plausible |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SSE201406 | 2014-06-10 | False | 1127 | 21 | 19.0 | 2014-06-11 00:44:00 | 23.293 | 2.126 | -0.152 | 6.525 | NW | 0.371 | 0.073 | 0.049 | 2.042 | False |
| SSE201411 | 2014-11-14 | True | 1824 | 43 | 37.0 | 2014-11-22 14:05:00 | 20.664 | -2.115 | -2.986 | 2.348 | SE (Thomas's southward) | -0.335 | 0.024 | 0.024 | -1.785 | False |
| SSE201507 | 2015-07-05 | False | 738 | 12 | 10.0 | 2015-07-21 07:20:00 | 8.107 | 0.376 | -5.728 | 11.227 | NW | 0.382 | 0.415 | 0.927 | -0.774 | False |
| SSE201510 | 2015-10-05 | False | 1444 | 17 | 15.0 | 2015-10-07 00:03:00 | 11.303 | -1.151 | -6.24 | 3.684 | SE (Thomas's southward) | -0.225 | 0.61 | 0.512 | -1.145 | False |
| SSE201512 | 2015-12-16 | True | 4348 | 82 | 72.0 | 2015-12-19 06:17:00 | 13.134 | -7.943 | -9.126 | -6.54 | SE (Thomas's southward) | -0.779 | 0.024 | 0.024 | -7.865 | True |
| SSE201608 | 2016-08-03 | False | 3027 | 25 | 21.0 | 2016-08-04 04:51:00 | 37.876 | -1.043 | -3.249 | -0.104 | SE (Thomas's southward) | -0.492 | 0.146 | 0.098 | -1.012 | True |
| SSE201703 | 2017-03-15 | True | 3225 | 77 | 69.0 | 2017-03-15 02:56:00 | 15.904 | 7.051 | 4.919 | 8.347 | NW | 0.607 | 0.024 | 0.024 | 7.22 | True |
| SSE201709 | 2017-09-17 | True | 2937 | 22 | 18.0 | 2017-09-21 00:35:00 | 24.883 | -0.055 | -2.035 | 2.924 | SE (Thomas's southward) | -0.028 | 0.829 | 0.902 | -0.233 | False |
| SSE201805 | 2018-04-30 | False | 1701 | 33 | 29.0 | 2018-05-02 17:06:00 | 13.82 | -3.168 | -4.443 | -1.953 | SE (Thomas's southward) | -0.682 | 0.024 | 0.024 | -2.931 | True |
| SSE201806 | 2018-06-17 | True | 8578 | 97 | 87.0 | 2018-06-20 07:46:00 | 18.699 | -4.557 | -5.289 | -3.518 | SE (Thomas's southward) | -0.695 | 0.024 | 0.024 | -4.758 | True |
| SSE201904 | 2019-04-19 | False | 1400 | 12 | 10.0 | 2019-04-20 19:28:00 | 13.156 | 0.286 | -23.178 | 10.497 | NW | -0.127 | 0.415 | 0.854 | -1.185 | False |
| SSE201905 | 2019-05-28 | False | 679 | 8 | 6.0 | 2019-06-06 15:10:00 | 2.861 | 28.816 | -91.323 | 142.803 | NW | 0.657 | 0.024 | 0.902 | 19.186 | False |
| SSE201908 | 2019-08-08 | True | 5121 | 64 | 56.0 | 2019-08-12 00:01:00 | 11.159 | 2.223 | 0.169 | 3.608 | NW | 0.17 | 0.024 | 0.073 | 2.089 | True |
| SSE201911 | 2019-11-05 | True | 4718 | 66 | 58.0 | 2019-11-08 09:02:00 | 16.86 | -2.035 | -3.14 | -1.39 | SE (Thomas's southward) | -0.681 | 0.024 | 0.024 | -2.089 | True |
| SSE202002 | 2020-02-01 | False | 642 | 14 | 12.0 | 2020-02-03 01:40:00 | 17.076 | 1.14 | -4.063 | 4.356 | NW | 0.224 | 0.512 | 0.293 | -0.235 | False |
| SSE202010 | 2020-10-12 | True | 3649 | 57 | 51.0 | 2020-10-15 02:01:00 | 30.354 | 3.804 | 0.469 | 5.874 | NW | 0.356 | 0.024 | 0.024 | 4.134 | True |
| SSE202101 | 2021-01-27 | True | 2534 | 63 | 55.0 | 2021-01-31 23:43:00 | 9.191 | -6.119 | -6.788 | -5.111 | SE (Thomas's southward) | -0.798 | 0.024 | 0.024 | -5.917 | True |
| SSE202106 | 2021-06-23 | False | 975 | 14 | 12.0 | 2021-06-25 01:55:00 | 22.013 | -1.129 | -5.574 | 5.039 | SE (Thomas's southward) | -0.385 | 0.537 | 0.171 | -1.279 | False |
| SSE202110 | 2021-09-30 | True | 6102 | 73 | 65.0 | 2021-10-01 21:33:00 | 26.202 | -0.69 | -1.443 | 0.124 | SE (Thomas's southward) | -0.228 | 0.049 | 0.098 | -0.757 | False |
| SSE202203 | 2022-03-14 | True | 5722 | 82 | 72.0 | 2022-03-23 02:28:00 | 14.248 | 0.085 | -1.111 | 1.36 | NW | 0.111 | 0.854 | 0.927 | 0.086 | False |
| SSE202207 | 2022-07-14 | True | 5815 | 54 | 48.0 | 2022-07-17 12:57:00 | 31.289 | -0.914 | -1.674 | 0.124 | SE (Thomas's southward) | -0.268 | 0.024 | 0.049 | -0.854 | False |
| SSE202212 | 2022-12-01 | False | 1654 | 25 | 21.0 | 2022-12-05 00:53:00 | 10.536 | 3.258 | 0.28 | 4.867 | NW | 0.473 | 0.024 | 0.073 | 3.582 | True |
| SSE202304 | 2023-04-02 | True | 7444 | 98 | 88.0 | 2023-04-05 00:33:00 | 15.325 | 4.222 | 3.459 | 5.162 | NW | 0.731 | 0.024 | 0.024 | 4.212 | True |
| SSE202306 | 2023-06-18 | False | 2030 | 33 | 29.0 | 2023-06-20 01:13:00 | 16.593 | 0.303 | -3.597 | 5.282 | NW | -0.017 | 0.585 | 0.756 | -0.326 | False |
| SSE202308 | 2023-08-13 | False | 4365 | 41 | 37.0 | 2023-08-16 03:40:00 | 37.019 | -0.495 | -1.694 | 2.166 | SE (Thomas's southward) | -0.146 | 0.22 | 0.585 | -0.363 | False |
| SSE202401 | 2024-01-26 | False | 1478 | 17 | 15.0 | 2024-01-26 10:08:00 | 18.844 | 4.432 | -0.351 | 19.199 | NW | 0.6 | 0.024 | 0.073 | 7.597 | False |
| SSE202407 | 2024-07-16 | True | 14946 | 103 | 91.0 | 2024-07-18 01:07:00 | 40.234 | 0.464 | -0.012 | 0.992 | NW | 0.18 | 0.049 | 0.024 | 0.591 | False |
| SSE202503 | 2025-03-20 | False | 825 | 7 | 5.0 | 2025-03-26 05:18:00 | 17.613 | -7.597 | -30.413 | 22.726 | SE (Thomas's southward) | -0.9 | 0.231 | 0.122 | -9.486 | False |
| SSE202508 | 2025-08-12 | True | 14015 | 94 | 84.0 | 2025-08-16 23:08:00 | 54.688 | 0.028 | -0.32 | 0.549 | NW | -0.021 | 0.78 | 0.756 | 0.018 | False |
| SSE202604 | 2026-04-03 | False | 4242 | 38 | 34.0 | 2026-04-06 00:46:00 | 50.319 | -0.336 | -1.918 | 3.326 | SE (Thomas's southward) | 0.026 | 0.268 | 0.561 | -0.261 | False |

**What it says, stated conservatively.** Of the 11 significant `cat_core` fronts, 6 propagate SE (Thomas's 'southward', the sense of the 2012 episode) and 6 propagate NW; |speed| spans 1.3–8.6 km/day with a median of 5.6, bracketing Bostock's 2012 −4.70 and Thomas's 7.3 km/day reference. Bidirectional along-strike propagation is expected in Cascadia, so a mix of senses is not by itself suspicious — but note that the four windows with the largest |speed| are also the ones whose SSE window is long enough (> 25 days) to contain more than one episode, and this estimator fits ONE front per window. Splitting multi-episode windows is the next refinement and is not done here.

### 3c. Per-SSE results (pre-registered estimator)

**cat_core** (31 windows with ≥ 50 events):

| sse_id | t0 | t1 | major | n_events | win_days | front_firstburst_kmday | front_lead_kmday | n_secondary | n_secondary_null_circ | n_secondary_null_dayshuf | reactivation_interval_d_p50 | frac_fam_reactivated |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SSE201406 | 2014-06-10 | 2014-07-07 | False | 537 | 0.0 |  |  | 2.0 | 0.15 | 0.6 | 2.491 | 0.014 |
| SSE201411 | 2014-11-14 | 2014-12-13 | True | 1094 | 7.329 |  | 14.409 | 10.0 | 2.0 | 2.9 | 1.001 | 0.026 |
| SSE201507 | 2015-07-05 | 2015-07-31 | False | 399 | 0.0 |  |  | 1.0 | 0.05 | 0.0 |  | 0.0 |
| SSE201510 | 2015-10-05 | 2015-10-21 | False | 854 | 11.137 |  | 0.0 | 0.0 | 1.15 | 0.4 | 2.741 | 0.033 |
| SSE201512 | 2015-12-16 | 2016-01-09 | True | 2947 | 14.503 | -7.841 | -6.653 | 11.0 | 7.3 | 2.4 | 1.02 | 0.12 |
| SSE201608 | 2016-08-03 | 2016-09-16 | False | 1464 | 21.369 |  | 1.416 | 1.0 | 0.3 | 0.15 | 1.497 | 0.037 |
| SSE201702 | 2017-02-10 | 2017-02-27 | False | 142 | 0.0 |  |  | 0.0 | 0.0 | 0.05 |  | 0.0 |
| SSE201703 | 2017-03-15 | 2017-04-06 | True | 1294 | 0.0 |  |  | 11.0 | 3.95 | 3.0 | 0.661 | 0.082 |
| SSE201709 | 2017-09-17 | 2017-10-27 | True | 1048 | 6.154 |  |  | 0.0 | 0.2 | 0.05 |  | 0.0 |
| SSE201805 | 2018-04-30 | 2018-05-20 | False | 421 | 0.0 |  |  | 1.0 | 0.0 | 0.2 |  | 0.0 |
| SSE201806 | 2018-06-17 | 2018-07-18 | True | 3618 | 12.789 | -5.823 | -5.543 | 20.0 | 11.35 | 2.15 | 0.965 | 0.102 |
| SSE201904 | 2019-04-19 | 2019-05-07 | False | 531 | 0.0 |  |  | 0.0 | 0.1 | 0.05 |  | 0.0 |
| SSE201905 | 2019-05-28 | 2019-06-19 | False | 165 |  |  |  | 0.0 | 0.0 | 0.35 |  | 0.0 |
| SSE201908 | 2019-08-08 | 2019-08-27 | True | 2227 | 2.594 |  | -5.863 | 16.0 | 12.7 | 1.7 | 0.558 | 0.09 |
| SSE201911 | 2019-11-05 | 2019-11-28 | True | 1903 | 11.593 | 0.038 | -0.301 | 16.0 | 10.25 | 2.15 | 0.953 | 0.09 |
| SSE202002 | 2020-02-01 | 2020-02-23 | False | 240 | 1.536 |  |  | 0.0 | 0.05 | 0.05 | 1.641 | 0.042 |
| SSE202010 | 2020-10-12 | 2020-11-19 | True | 1513 | 0.0 |  |  | 12.0 | 5.4 | 1.85 | 0.797 | 0.088 |
| SSE202101 | 2021-01-27 | 2021-02-18 | True | 993 | 5.008 | -6.048 | -5.922 | 12.0 | 4.3 | 1.75 | 1.03 | 0.065 |
| SSE202106 | 2021-06-23 | 2021-07-17 | False | 207 | 0.959 |  | -11.157 | 1.0 | 0.0 | 0.1 |  | 0.0 |
| SSE202110 | 2021-09-30 | 2021-11-01 | True | 2010 |  |  |  | 9.0 | 4.8 | 2.4 | 0.921 | 0.047 |
| SSE202203 | 2022-03-14 | 2022-04-10 | True | 1600 | 0.0 |  |  | 7.0 | 6.45 | 1.5 | 1.042 | 0.087 |
| SSE202207 | 2022-07-14 | 2022-09-04 | True | 1721 | 17.557 |  | -4.151 | 3.0 | 0.15 | 0.55 | 1.2 | 0.027 |
| SSE202212 | 2022-12-01 | 2022-12-18 | False | 526 | 0.0 |  |  | 1.0 | 0.05 | 0.0 |  | 0.0 |
| SSE202304 | 2023-04-02 | 2023-04-27 | True | 2640 | 5.462 | 9.809 | -2.62 | 13.0 | 8.9 | 2.85 | 1.308 | 0.123 |
| SSE202306 | 2023-06-18 | 2023-07-10 | False | 522 |  |  |  | 0.0 | 0.0 | 0.0 |  | 0.0 |
| SSE202308 | 2023-08-13 | 2023-10-07 | False | 919 | 21.79 |  | -2.386 | 2.0 | 0.2 | 0.1 | 0.838 | 0.044 |
| SSE202401 | 2024-01-26 | 2024-02-22 | False | 504 | 0.0 |  |  | 1.0 | 0.0 | 0.2 |  | 0.0 |
| SSE202407 | 2024-07-16 | 2024-09-08 | True | 4742 | 19.111 |  | -2.033 | 12.0 | 15.15 | 2.35 | 1.024 | 0.22 |
| SSE202503 | 2025-03-20 | 2025-04-14 | False | 204 | 0.0 |  |  | 0.0 | 0.0 | 0.0 | 0.54 | 0.037 |
| SSE202508 | 2025-08-12 | 2025-10-22 | True | 4401 | 40.065 | -1.959 | -0.513 | 17.0 | 10.15 | 3.45 | 0.609 | 0.139 |
| SSE202604 | 2026-04-03 | 2026-05-30 | False | 1246 | 37.597 | -1.778 | -0.627 | 3.0 | 0.15 | 0.25 | 2.775 | 0.016 |

**cat_ext** (31 windows with ≥ 50 events):

| sse_id | t0 | t1 | major | n_events | win_days | front_firstburst_kmday | front_lead_kmday | n_secondary | n_secondary_null_circ | n_secondary_null_dayshuf | reactivation_interval_d_p50 | frac_fam_reactivated |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SSE201406 | 2014-06-10 | 2014-07-07 | False | 1127 | 2.263 |  | -11.42 | 0.0 | 1.1 | 0.7 | 0.942 | 0.01 |
| SSE201411 | 2014-11-14 | 2014-12-13 | True | 1824 | 7.709 |  | -4.313 | 11.0 | 5.65 | 2.05 | 0.883 | 0.039 |
| SSE201507 | 2015-07-05 | 2015-07-31 | False | 738 | 1.401 |  |  | 0.0 | 0.15 | 0.05 |  | 0.0 |
| SSE201510 | 2015-10-05 | 2015-10-21 | False | 1444 | 11.204 | -2.002 | 1.27 | 1.0 | 2.1 | 0.45 | 0.926 | 0.036 |
| SSE201512 | 2015-12-16 | 2016-01-09 | True | 4348 | 14.508 | -8.554 | -5.58 | 11.0 | 6.7 | 2.45 | 1.045 | 0.115 |
| SSE201608 | 2016-08-03 | 2016-09-16 | False | 3027 | 24.092 | -3.278 | 0.198 | 1.0 | 3.1 | 0.8 | 0.991 | 0.047 |
| SSE201702 | 2017-02-10 | 2017-02-27 | False | 349 | 0.0 |  |  | 1.0 | 0.05 | 0.3 |  | 0.0 |
| SSE201703 | 2017-03-15 | 2017-04-06 | True | 3225 | 6.61 |  | -5.074 | 9.0 | 4.4 | 2.1 | 0.86 | 0.096 |
| SSE201709 | 2017-09-17 | 2017-10-27 | True | 2937 | 13.243 |  | -2.537 | 1.0 | 3.0 | 0.25 | 0.582 | 0.009 |
| SSE201805 | 2018-04-30 | 2018-05-20 | False | 1701 | 11.195 | -0.111 | -1.409 | 1.0 | 2.65 | 0.8 | 0.928 | 0.009 |
| SSE201806 | 2018-06-17 | 2018-07-18 | True | 8578 | 21.521 | -4.877 | -1.779 | 5.0 | 4.8 | 2.1 | 1.094 | 0.148 |
| SSE201904 | 2019-04-19 | 2019-05-07 | False | 1400 | 1.346 |  | -19.96 | 1.0 | 1.15 | 0.3 | 0.902 | 0.014 |
| SSE201905 | 2019-05-28 | 2019-06-19 | False | 679 | 0.0 |  |  | 1.0 | 0.95 | 0.4 | 0.655 | 0.011 |
| SSE201908 | 2019-08-08 | 2019-08-27 | True | 5121 | 2.61 | -8.582 | -5.286 | 7.0 | 8.35 | 0.75 | 0.994 | 0.108 |
| SSE201911 | 2019-11-05 | 2019-11-28 | True | 4718 | 11.75 | -0.055 | -1.109 | 11.0 | 12.6 | 0.75 | 0.762 | 0.134 |
| SSE202002 | 2020-02-01 | 2020-02-23 | False | 642 | 15.561 |  | -1.245 | 0.0 | 0.75 | 0.2 | 4.167 | 0.024 |
| SSE202010 | 2020-10-12 | 2020-11-19 | True | 3649 | 22.58 |  | -0.267 | 9.0 | 17.2 | 1.6 | 0.77 | 0.09 |
| SSE202101 | 2021-01-27 | 2021-02-18 | True | 2534 | 4.25 | -7.208 | -5.055 | 9.0 | 10.5 | 0.9 | 1.029 | 0.083 |
| SSE202106 | 2021-06-23 | 2021-07-17 | False | 975 | 1.684 |  | -9.873 | 2.0 | 1.0 | 0.35 | 1.855 | 0.03 |
| SSE202110 | 2021-09-30 | 2021-11-01 | True | 6102 | 0.0 |  |  | 5.0 | 7.3 | 1.05 | 0.894 | 0.071 |
| SSE202203 | 2022-03-14 | 2022-04-10 | True | 5722 | 1.52 |  | -0.357 | 6.0 | 4.4 | 0.9 | 0.567 | 0.132 |
| SSE202207 | 2022-07-14 | 2022-09-04 | True | 5815 | 16.955 | -4.561 | -0.28 | 3.0 | 6.6 | 1.2 | 0.995 | 0.053 |
| SSE202212 | 2022-12-01 | 2022-12-18 | False | 1654 |  |  |  | 1.0 | 2.1 | 0.6 | 1.277 | 0.009 |
| SSE202304 | 2023-04-02 | 2023-04-27 | True | 7444 | 4.467 | 1.622 | -0.896 | 10.0 | 3.95 | 1.65 | 0.796 | 0.184 |
| SSE202306 | 2023-06-18 | 2023-07-10 | False | 2030 | 8.274 |  | -0.314 | 0.0 | 4.15 | 0.45 | 0.935 | 0.009 |
| SSE202308 | 2023-08-13 | 2023-10-07 | False | 4365 | 8.94 | 2.808 | -0.038 | 1.0 | 6.1 | 0.9 | 1.639 | 0.07 |
| SSE202401 | 2024-01-26 | 2024-02-22 | False | 1478 | 0.0 |  |  | 0.0 | 0.8 | 0.15 | 1.775 | 0.028 |
| SSE202407 | 2024-07-16 | 2024-09-08 | True | 14946 | 26.797 | -1.464 | -0.238 | 7.0 | 8.9 | 0.75 | 0.871 | 0.254 |
| SSE202503 | 2025-03-20 | 2025-04-14 | False | 825 | 0.0 |  |  | 1.0 | 0.55 | 0.35 | 2.12 | 0.02 |
| SSE202508 | 2025-08-12 | 2025-10-22 | True | 14015 | 37.786 | -1.961 | -0.326 | 6.0 | 13.25 | 1.25 | 0.602 | 0.158 |
| SSE202604 | 2026-04-03 | 2026-05-30 | False | 4242 | 0.0 |  |  | 2.0 | 5.55 | 1.0 | 1.26 | 0.062 |

**cat_core_hp** (31 windows with ≥ 50 events):

| sse_id | t0 | t1 | major | n_events | win_days | front_firstburst_kmday | front_lead_kmday | n_secondary | n_secondary_null_circ | n_secondary_null_dayshuf | reactivation_interval_d_p50 | frac_fam_reactivated |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SSE201406 | 2014-06-10 | 2014-07-07 | False | 382 | 0.0 |  |  | 0.0 | 0.1 | 0.05 |  | 0.0 |
| SSE201411 | 2014-11-14 | 2014-12-13 | True | 861 | 0.0 |  |  | 8.0 | 0.8 | 2.2 | 1.001 | 0.036 |
| SSE201507 | 2015-07-05 | 2015-07-31 | False | 287 |  |  |  | 0.0 | 0.0 | 0.0 |  | 0.0 |
| SSE201510 | 2015-10-05 | 2015-10-21 | False | 407 | 0.0 |  |  | 0.0 | 0.15 | 0.0 | 0.926 | 0.021 |
| SSE201512 | 2015-12-16 | 2016-01-09 | True | 2667 | 14.507 | -8.083 | -6.593 | 10.0 | 7.85 | 2.7 | 1.012 | 0.133 |
| SSE201608 | 2016-08-03 | 2016-09-16 | False | 642 | 21.462 |  | -0.443 | 0.0 | 0.1 | 0.0 | 2.691 | 0.032 |
| SSE201702 | 2017-02-10 | 2017-02-27 | False | 81 |  |  |  |  |  |  |  | 0.0 |
| SSE201703 | 2017-03-15 | 2017-04-06 | True | 1182 | 0.0 |  |  | 14.0 | 3.05 | 2.75 | 0.661 | 0.098 |
| SSE201709 | 2017-09-17 | 2017-10-27 | True | 864 | 0.0 |  |  | 0.0 | 0.0 | 0.2 |  | 0.0 |
| SSE201805 | 2018-04-30 | 2018-05-20 | False | 346 | 0.0 |  |  | 0.0 | 0.0 | 0.25 |  | 0.0 |
| SSE201806 | 2018-06-17 | 2018-07-18 | True | 3159 | 11.545 | -6.709 | -5.948 | 22.0 | 11.0 | 3.4 | 0.912 | 0.099 |
| SSE201904 | 2019-04-19 | 2019-05-07 | False | 422 |  |  |  | 0.0 | 0.0 | 0.0 |  | 0.0 |
| SSE201905 | 2019-05-28 | 2019-06-19 | False | 144 |  |  |  | 0.0 | 0.0 | 0.15 |  | 0.0 |
| SSE201908 | 2019-08-08 | 2019-08-27 | True | 1553 | 0.0 |  |  | 17.0 | 7.85 | 1.75 | 0.714 | 0.109 |
| SSE201911 | 2019-11-05 | 2019-11-28 | True | 1680 | 14.27 | -2.073 | -3.483 | 11.0 | 9.9 | 2.1 | 0.754 | 0.082 |
| SSE202002 | 2020-02-01 | 2020-02-23 | False | 214 | 0.0 |  |  | 0.0 | 0.0 | 0.0 | 1.641 | 0.059 |
| SSE202010 | 2020-10-12 | 2020-11-19 | True | 1292 | 0.0 |  |  | 12.0 | 3.2 | 1.9 | 0.797 | 0.111 |
| SSE202101 | 2021-01-27 | 2021-02-18 | True | 887 | 5.002 | -6.25 | -5.584 | 11.0 | 3.25 | 2.25 | 1.037 | 0.065 |
| SSE202106 | 2021-06-23 | 2021-07-17 | False | 59 |  |  |  |  |  |  |  |  |
| SSE202110 | 2021-09-30 | 2021-11-01 | True | 915 | 9.169 |  | -2.566 | 8.0 | 1.1 | 1.55 | 0.977 | 0.037 |
| SSE202203 | 2022-03-14 | 2022-04-10 | True | 932 | 0.308 |  |  | 8.0 | 2.35 | 2.35 | 1.031 | 0.115 |
| SSE202207 | 2022-07-14 | 2022-09-04 | True | 147 | 0.0 |  |  | 0.0 | 0.0 | 0.1 |  | 0.0 |
| SSE202212 | 2022-12-01 | 2022-12-18 | False | 52 |  |  |  |  |  |  |  | 0.0 |
| SSE202304 | 2023-04-02 | 2023-04-27 | True | 1815 | 5.459 | 8.399 | -2.599 | 17.0 | 7.1 | 3.8 | 1.97 | 0.137 |
| SSE202306 | 2023-06-18 | 2023-07-10 | False | 71 |  |  |  |  |  |  |  | 0.0 |
| SSE202308 | 2023-08-13 | 2023-10-07 | False | 208 |  |  |  | 0.0 | 0.0 | 0.0 |  | 0.0 |
| SSE202401 | 2024-01-26 | 2024-02-22 | False | 86 |  |  |  |  |  |  |  | 0.0 |
| SSE202407 | 2024-07-16 | 2024-09-08 | True | 2407 | 0.0 |  |  | 21.0 | 7.35 | 2.95 | 1.01 | 0.173 |
| SSE202503 | 2025-03-20 | 2025-04-14 | False | 79 |  |  |  |  |  |  |  | 0.0 |
| SSE202508 | 2025-08-12 | 2025-10-22 | True | 1581 | 24.916 |  | 0.0 | 8.0 | 0.9 | 2.6 | 0.597 | 0.127 |
| SSE202604 | 2026-04-03 | 2026-05-30 | False | 333 | 0.0 |  |  | 0.0 | 0.0 | 0.0 |  | 0.0 |

Catalog size inside every SSE window (including windows too small to run the harness):

| sse_id | cat_core | cat_core_hp | cat_ext |
|---|---|---|---|
| SSE201406 | 537.0 | 382.0 | 1127.0 |
| SSE201411 | 1094.0 | 861.0 | 1824.0 |
| SSE201507 | 399.0 | 287.0 | 738.0 |
| SSE201510 | 854.0 | 407.0 | 1444.0 |
| SSE201512 | 2947.0 | 2667.0 | 4348.0 |
| SSE201608 | 1464.0 | 642.0 | 3027.0 |
| SSE201702 | 142.0 | 81.0 | 349.0 |
| SSE201703 | 1294.0 | 1182.0 | 3225.0 |
| SSE201709 | 1048.0 | 864.0 | 2937.0 |
| SSE201805 | 421.0 | 346.0 | 1701.0 |
| SSE201806 | 3618.0 | 3159.0 | 8578.0 |
| SSE201904 | 531.0 | 422.0 | 1400.0 |
| SSE201905 | 165.0 | 144.0 | 679.0 |
| SSE201908 | 2227.0 | 1553.0 | 5121.0 |
| SSE201911 | 1903.0 | 1680.0 | 4718.0 |
| SSE202002 | 240.0 | 214.0 | 642.0 |
| SSE202010 | 1513.0 | 1292.0 | 3649.0 |
| SSE202101 | 993.0 | 887.0 | 2534.0 |
| SSE202106 | 207.0 | 59.0 | 975.0 |
| SSE202110 | 2010.0 | 915.0 | 6102.0 |
| SSE202203 | 1600.0 | 932.0 | 5722.0 |
| SSE202207 | 1721.0 | 147.0 | 5815.0 |
| SSE202212 | 526.0 | 52.0 | 1654.0 |
| SSE202304 | 2640.0 | 1815.0 | 7444.0 |
| SSE202306 | 522.0 | 71.0 | 2030.0 |
| SSE202308 | 919.0 | 208.0 | 4365.0 |
| SSE202401 | 504.0 | 86.0 | 1478.0 |
| SSE202407 | 4742.0 | 2407.0 | 14946.0 |
| SSE202503 | 204.0 | 79.0 | 825.0 |
| SSE202508 | 4401.0 | 1581.0 | 14015.0 |
| SSE202604 | 1246.0 | 333.0 | 4242.0 |

## 4. Known limitations

1. **Completeness is ~15 % (T_core) / ~18 % (T_ext) of Bostock-class events** on the 15-station set at the referee thresholds, and it varies by epoch because the thresholds do: per epoch, at each epoch's own T, E0 14.6 %; E1 12.9 %; E2 12.9 %; E3 12.2 % at T_core (`completeness_by_epoch_summary.csv`). This catalog is a decimated view of the LFE population, not an equivalent of Bostock's. Every decimation-sensitive statistic must be compared against the per-family-thinned reference, not against his raw catalog.
2. **Completeness is strongly family-dependent**: per-family own recall on the S15 set at T_core has p10/p50/p90 = 0.000 / 0.091 / 0.298, and 14 of 130 families have no usable template at all and are NEVER detected. Any spatial statistic (front position, along-strike coverage) inherits this and must use the per-family recall column.
3. **Time basis**: Bostock's DB clock plus a fitted per-family offset — not an origin time (section 2).
4. **UW.LRIV is dropped after 2019-05-30** by Track A-2's epoch-identity result. The ablation above shows the loss is free (recall and false alarms both unchanged), and every family-day records it in `excluded_sta` / `n_sta_excluded`.
5. **CN.SNB ends 2021-03-28** and that IS costly. At the frozen threshold the 13-station network sits at ~3× the quiet-day false-alarm rate it was frozen at; Track C's per-epoch re-derivation responds by raising E3's T_core to 9.25, which restores the false-alarm budget and pays for it in recall (E3 T_core completeness 0.122 against E0's 0.146, −16 %). Either way the post-2021 record is the weakest part of the catalog; the per-event λ and Mc columns carry it.
6. **The post-2017 HF veto is provisional.** The cut is 3.6384 on the PGC BH epochs and 3.5129 on the HH epochs (Track A-2 RULE-M); the HH cut's 95 % bootstrap CI contains the frozen value, but it is calibrated from one SSE per epoch and the SSE-to-SSE spread inside an epoch is ±35 %. `hfr`, `hf_fs` and `hf_cut` are stored on every row so any re-calibration is a re-flag, not a re-run (`scripts/c_reflag_hf.py`).
7. **Quiet-day detections are not all false alarms.** Track A step 3b measured ~22 % real inter-ETS repeats at T_core (≤ 14 % at T_ext); λ is corrected by that fraction and `lambda_fa_low` carries the 95 %-upper-bound version. Only one family passes the withheld-station test individually, so no inter-ETS detection is promoted to a confirmed event.
8. **Magnitudes.** Our detections carry `rel_amp` (an SNR-weighted cc·√energy ratio in units of each station's 2012 quiet-day window rms), not a magnitude. Mc is expressed on Bostock's magnitude scale through the referee fit; it is a completeness statement about Bostock-class events, not a magnitude for our own events.
9. **The secondary-front count is uninformative at this completeness.** In the thinned reference the circular null is itself 10–14 detections against Bostock's 2.25 — the estimator's null scales with catalog density, so the count is reported and never interpreted.
10. **The quiet-day floor is seasonally unstable, and not in a fixed direction.** All eight epoch-tier cells miss Merlin's ≤ 2× winter/summer band, but E0 misses it high (winter 2.4× summer) and E1–E3 miss it low (summer 2–4× winter). Section 1d shows the floor is decoupled from the ambient level — winter is the loudest season in every epoch — so no seasonal correction is applied; λ is carried per (family, tier, epoch, season) instead.
11. **The PRE-REGISTERED main-front estimator mostly does not work at this purity** — a plausible southward speed in only 3 of 31 SSE windows at `cat_core` — because its leading edge is the 10th percentile of event position and the along-strike false-alarm background pins it. The purity-matched replacement of §3b-2 (post-registration, approved 2026-09-17) fits 27 of 31 and finds 11 significant, but it is a NEW estimator: its α was fixed on the 2012 referee, and the pre-registered bars in §3a are still the bars. Its own remaining limitation is that it fits ONE front per SSE window, and the four largest |speeds| come from windows long enough (> 25 days) to hold more than one episode.

