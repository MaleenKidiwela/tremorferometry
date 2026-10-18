# Does the Washington five-clause screen survive contact with Bostock's own families?

Agent Q, 2026-09-21.  Workspace `results/q_screenval/`; queue `scripts/q_sv_queue.sh`
(markers `results/Q_SV_DONE|FAILED`); code `scripts/q_sv_data.py`, `q_sv_screen.py`,
`q_sv_tables.py`, `q_sv_wa.py`.

**The objection.** The five-clause screen applied to every Washington candidate
(`scripts/n_screen.py`, `results/WA_2008_RESULTS.md` sec.4) contains criteria that Bostock and
Chestler never used.  Arguing about them is pointless; the question is empirical and is now
measured.  Every clause, **exactly as coded**, has been applied to the grown templates of the
**58 Bostock development families that the same recipe rediscovered blind** on the 2005 POLARIS
bed, to the **64 templates matching the 18 located sealed families**, to the 18 candidate-new
templates, and to the **1,201 unlocatable/unmatched templates** (the noise population).  Those
58 families are not candidates: they are Bostock's catalogue, matched at his own &ge;30 %
one-to-one &plusmn;2 s bar (`results/REPRO_BOSTOCK.md` sec.8-sec.11, `NORTH_DISCOVERY_SUMMARY.md`).

> **Headline.  The screen as coded keeps 2 of those 58 families (3 %).  Without clause (e) it
> keeps 7 (12 %).  A screen that rejects 97 % of Bostock's own catalogue cannot be used to
> conclude anything from "0 candidates" in Washington.**

---

## 1. What was measured, and on what

| | |
|---|---|
| bed | 2005 POLARIS, 15 stations, 24 ETS days 2005-08-29&hellip;09-21 (the K_REPRO3 run) |
| templates | `results/repro/grow15/tpl_it0.npz`, **1,460** grown groups |
| ETS-day detections | `results/repro/grow15/det_it0.parquet`, 64,959 on 24 days |
| kept/hidden match table | `results/repro/classify15.csv` (177 kept-matching templates &rarr; 58 families; 64 hidden-matching &rarr; 18; 18 candidate-new; 1,201 unlocatable) |
| locations | `results/repro/loc15_it0.csv` &mdash; `k_locate.py`, grid 47.90&ndash;49.28 N / &minus;124.90&hellip;&minus;122.82 / 18&ndash;56 km, 50-realisation shuffled-station control, Slab2 via `p_common.slab2_interp` |
| **NEW, clause (a)** | the same templates re-fired unchanged at 8 &times; MAD over the **20 literature-quiet 2005 days** (2005-06-25&hellip;07-14, `k_common.QUIET_DAYS`; 7&ndash;12 stations cached per day) &rarr; `results/q_screenval/det_quiet15.parquet`, 43,224 detections |
| **NEW, clause (e)** | first-half and second-half re-stacks of every template on its own detections (`k_grow.build_templates` verbatim, the same rule `n_relocate.py --half` uses) and re-location on the same grid &rarr; `reloc_h1.csv`, `reloc_h2.csv` |
| clause (b) reference | 2005 has no tremor catalogue, so the tremor-hour set is the hour-binned **Bostock LFE catalogue** with **the template's own family left out**; both an *any-LFE-in-the-hour* set (base rate 0.77) and a **base-rate-matched** set (top 23 % of hours, base 0.23, matching WA 2008's audited base rate of 0.24) are reported |

Unit of account.  The pre-stated rule speaks of *families*, so a clause passes a family if **any**
of that family's templates passes it &mdash; the same convention under which the family counts as
rediscovered.  Per-template fractions are given beside it.

### 1.1 The quiet-day false-alarm rate, north vs Washington

This is the number the 2008 write-up used to diagnose the Washington templates as "noise
repeaters, not tremor-locked":

| bed | quiet-day detections per template-day | tremor-day rate | ratio |
|---|---|---|---|
| **2005 POLARIS (58 Bostock families rediscovered)** | **1.480** (20 quiet days) | **1.854** (24 ETS days) | **1.25** |
| Puget/S-Olympic 2008 (agent N, real chain) | 1.525 (19 quiet days) | 1.541 (20 tremor days) | 1.01 |
| 2010 templates re-fired in 2011 (`WA_2011_REFIRE_RESULTS.md`) | 1.89 | 0.66 | 0.35 |

The network-level statistic that condemned the 2008 bed puts the bed on which the recipe
demonstrably works only 1.25&times; away from it.  Only the 2011 re-fire (ratio 0.35, i.e.
*anti*-correlated with tremor) is qualitatively different.  Caveat: the 2005 quiet days carry
7&ndash;12 cached stations against 9&ndash;15 on the ETS days, and a thinner network makes the
8 &times; MAD bar noisier &mdash; the same confound the Washington comparison has.

---

## 2. Every clause, applied as coded

| clause (as coded) | TRUE families (58) | their 177 templates | hidden families (18) | their 64 templates | 18 candidate-new | noise population (1,201) |
|---|---|---|---|---|---|---|
| (a) quiet-day FA &le; 0.1/tpl-day AND tremor rate &ge; 3&times; it | **49/58 = 84%** | 140/177 = 79% | 17/18 = 94% | 46/64 = 72% | 14/18 | 662/1201 = 55% |
| (b) tremor co-timing > same-day base rate (any-LFE hours, base 0.77) | **33/58 = 57%** | 59/177 = 33% | 14/18 = 78% | 41/64 = 64% | 9/18 | 659/1201 = 55% |
| (c1) located on &ge; 5 stations | **50/58 = 86%** | 133/177 = 75% | 16/18 = 89% | 47/64 = 73% | 18/18 | 399/1201 = 33% |
| (c2) shuffled-station control passes | **31/58 = 53%** | 52/177 = 29% | 11/18 = 61% | 23/64 = 36% | 18/18 | 56/1201 = 5% |
| (c3) INTERIOR to the search grid | **22/58 = 38%** | 31/177 = 18% | 8/18 = 44% | 15/64 = 23% | 6/18 | 52/1201 = 4% |
| (c4) &#124;&Delta;z&#124; &le; 8 km of Slab2 | **28/58 = 48%** | 54/177 = 31% | 8/18 = 44% | 17/64 = 27% | 18/18 | 113/1201 = 9% |
| (d) &ge; 3 distinct tremor days | **52/58 = 90%** | 127/177 = 72% | 14/18 = 78% | 47/64 = 73% | 18/18 | 958/1201 = 80% |
| (e) 2nd-half re-stack relocates within 10 km | **13/58 = 22%** | 17/177 = 10% | 8/18 = 44% | 12/64 = 19% | 2/18 | 22/1201 = 2% |

Per-family rows for all 58: sec.2.3 and `results/q_screenval/per_family_58.csv`.
Per-template table with every column: `results/q_screenval/screen_north.csv`.

### 2.1 Cumulative, in `n_screen.py`'s own order

| cumulative, n_screen's order | TRUE families left | their templates | hidden families | their templates | candidate-new | noise templates |
|---|---|---|---|---|---|---|
| (c1) located on &ge; 5 stations | **50/58** | 133 | 16/18 | 47 | 18 | 399 |
| (c2) shuffled-station control passes | **31/58** | 52 | 11/18 | 23 | 18 | 56 |
| (c3) INTERIOR to the search grid | **17/58** | 19 | 8/18 | 14 | 6 | 14 |
| (c4) &#124;&Delta;z&#124; &le; 8 km of Slab2 | **11/58** | 11 | 5/18 | 10 | 6 | 0 |
| (d) &ge; 3 distinct tremor days | **11/58** | 11 | 5/18 | 10 | 6 | 0 |
| (a) quiet-day FA &le; 0.1/tpl-day AND tremor rate &ge; 3&times; it | **10/58** | 10 | 5/18 | 9 | 6 | 0 |
| (b) tremor co-timing > same-day base rate (any-LFE hours, base 0.77) | **7/58** | 7 | 4/18 | 8 | 2 | 0 |
| (e) 2nd-half re-stack relocates within 10 km | **2/58** | 2 | 2/18 | 5 | 0 | 0 |

### 2.2 Whole screens

| screen | TRUE families recovered | their templates | hidden families | candidate-new | noise templates |
|---|---|---|---|---|---|
| the screen as coded, all five clauses | **2/58 = 3%** | 2 | 2/18 | 0/18 | 0/1201 = 0.0% |
| as coded without (e) | **7/58 = 12%** | 7 | 4/18 | 2/18 | 0/1201 = 0.0% |
| as coded without (c3) and (e) | **10/58 = 17%** | 13 | 4/18 | 7/18 | 0/1201 = 0.0% |
| justified clause only: (c1) | **50/58 = 86%** | 133 | 16/18 | 18/18 | 399/1201 = 33.2% |
| (c1)+(c4&prime;) only | **43/58 = 74%** | 106 | 14/18 | 18/18 | 273/1201 = 22.7% |
| (c1)+(c4&prime;)+(d) | **38/58 = 66%** | 83 | 10/18 | 18/18 | 196/1201 = 16.3% |
| **recommended tier A** = (c1)+(c4&prime;)+(d)+(a&prime;) | **33/58 = 57%** | 69 | 9/18 | 15/18 | 132/1201 = 11.0% |
| **recommended tier B** = tier A + control (c2) | **21/58 = 36%** | 33 | 7/18 | 15/18 | 16/1201 = 1.3% |
| **recommended tier C** = tier B + co-timing (b&prime;) | **13/58 = 22%** | 16 | 6/18 | 10/18 | 10/1201 = 0.8% |

The last three rows are the recommendation of sec.5.

### 2.3 Per-family table (58 rediscovered Bostock families)

Y = at least one of that family's templates passes; `best |dz|` = the smallest |&Delta;z| among its
located templates.

| bkey | templates | max n_det | max days | c_a | c_b | c_c1 | c_c2 | c_c3 | c_c4 | c_d | c_e | best &#124;&Delta;z&#124; | clauses passed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 1 | 9 | 1 | Y | &middot; | &middot; | &middot; | &middot; | &middot; | &middot; | &middot; | &mdash; | 1/8 |
| 3 | 5 | 58 | 9 | Y | Y | Y | Y | Y | Y | Y | Y | 0.6 | 8/8 |
| 5 | 5 | 31 | 8 | Y | Y | Y | Y | Y | Y | Y | Y | 0.1 | 8/8 |
| 6 | 5 | 19 | 7 | Y | Y | Y | Y | &middot; | Y | Y | &middot; | 0.5 | 6/8 |
| 9 | 1 | 11 | 7 | &middot; | Y | Y | &middot; | Y | &middot; | Y | &middot; | 12.1 | 4/8 |
| 13 | 2 | 9 | 5 | &middot; | &middot; | Y | &middot; | &middot; | Y | Y | &middot; | 1.9 | 3/8 |
| 14 | 3 | 10 | 3 | Y | Y | Y | &middot; | &middot; | &middot; | Y | &middot; | 29.3 | 4/8 |
| 16 | 3 | 44 | 7 | Y | Y | Y | Y | Y | Y | Y | Y | 1.3 | 8/8 |
| 20 | 3 | 11 | 4 | Y | Y | Y | &middot; | &middot; | Y | Y | &middot; | 6.7 | 5/8 |
| 21 | 1 | 7 | 3 | Y | &middot; | Y | &middot; | &middot; | &middot; | Y | &middot; | 8.0 | 3/8 |
| 23 | 6 | 26 | 4 | Y | &middot; | Y | Y | &middot; | &middot; | Y | &middot; | 8.2 | 4/8 |
| 30 | 1 | 38 | 3 | Y | &middot; | Y | Y | &middot; | &middot; | Y | &middot; | 32.6 | 4/8 |
| 33 | 3 | 54 | 5 | &middot; | Y | Y | Y | Y | Y | Y | &middot; | 7.2 | 6/8 |
| 37 | 3 | 14 | 4 | Y | &middot; | Y | Y | &middot; | Y | Y | &middot; | 0.6 | 5/8 |
| 38 | 2 | 18 | 4 | Y | &middot; | Y | &middot; | Y | Y | Y | &middot; | 0.5 | 5/8 |
| 39 | 1 | 7 | 3 | &middot; | &middot; | Y | &middot; | &middot; | &middot; | Y | &middot; | 8.4 | 2/8 |
| 41 | 1 | 32 | 8 | &middot; | Y | Y | Y | &middot; | &middot; | Y | &middot; | 22.2 | 4/8 |
| 42 | 1 | 5 | 3 | Y | Y | &middot; | &middot; | &middot; | &middot; | Y | &middot; | &mdash; | 3/8 |
| 44 | 8 | 11 | 5 | Y | Y | Y | Y | Y | Y | Y | &middot; | 2.7 | 7/8 |
| 46 | 3 | 12 | 5 | Y | &middot; | Y | &middot; | &middot; | &middot; | Y | &middot; | 29.2 | 3/8 |
| 47 | 5 | 36 | 6 | Y | Y | Y | Y | Y | Y | Y | Y | 0.4 | 8/8 |
| 49 | 1 | 9 | 2 | &middot; | &middot; | Y | Y | &middot; | &middot; | &middot; | &middot; | 18.1 | 2/8 |
| 52 | 1 | 13 | 4 | Y | &middot; | Y | Y | &middot; | &middot; | Y | &middot; | 17.6 | 4/8 |
| 53 | 1 | 6 | 1 | Y | &middot; | &middot; | &middot; | &middot; | &middot; | &middot; | &middot; | &mdash; | 1/8 |
| 54 | 7 | 53 | 7 | Y | Y | Y | Y | Y | Y | Y | Y | 0.4 | 8/8 |
| 55 | 4 | 37 | 6 | Y | Y | Y | &middot; | &middot; | &middot; | Y | &middot; | 18.9 | 4/8 |
| 56 | 8 | 66 | 19 | Y | Y | Y | Y | Y | Y | Y | Y | 4.0 | 8/8 |
| 61 | 2 | 10 | 3 | Y | &middot; | Y | &middot; | &middot; | &middot; | Y | &middot; | 10.6 | 3/8 |
| 64 | 1 | 4 | 1 | Y | &middot; | Y | Y | &middot; | &middot; | &middot; | &middot; | 15.8 | 3/8 |
| 65 | 6 | 34 | 6 | Y | Y | Y | Y | Y | &middot; | Y | Y | 15.6 | 7/8 |
| 66 | 5 | 19 | 8 | Y | Y | Y | Y | &middot; | Y | Y | &middot; | 1.7 | 6/8 |
| 67 | 1 | 12 | 4 | Y | &middot; | Y | &middot; | &middot; | &middot; | Y | &middot; | 39.9 | 3/8 |
| 68 | 3 | 16 | 7 | Y | Y | Y | &middot; | Y | Y | Y | &middot; | 3.9 | 6/8 |
| 71 | 2 | 5 | 3 | Y | Y | Y | &middot; | &middot; | &middot; | Y | &middot; | 16.3 | 4/8 |
| 72 | 5 | 34 | 6 | Y | Y | Y | Y | Y | Y | Y | Y | 2.0 | 8/8 |
| 73 | 2 | 14 | 3 | Y | &middot; | Y | Y | &middot; | &middot; | Y | &middot; | 35.0 | 4/8 |
| 74 | 1 | 6 | 3 | Y | &middot; | Y | &middot; | &middot; | &middot; | Y | &middot; | 10.0 | 3/8 |
| 75 | 1 | 5 | 1 | Y | &middot; | &middot; | &middot; | &middot; | &middot; | &middot; | &middot; | &mdash; | 1/8 |
| 76 | 1 | 15 | 5 | Y | &middot; | Y | &middot; | &middot; | &middot; | Y | &middot; | 9.3 | 3/8 |
| 81 | 3 | 14 | 4 | Y | &middot; | Y | &middot; | &middot; | &middot; | Y | Y | 14.1 | 4/8 |
| 84 | 1 | 3 | 1 | &middot; | Y | &middot; | &middot; | &middot; | &middot; | &middot; | &middot; | &mdash; | 1/8 |
| 86 | 7 | 35 | 4 | Y | Y | Y | Y | &middot; | Y | Y | Y | 1.4 | 7/8 |
| 91 | 2 | 13 | 5 | Y | Y | Y | Y | Y | Y | Y | &middot; | 0.5 | 7/8 |
| 92 | 3 | 17 | 3 | Y | Y | Y | Y | Y | Y | Y | &middot; | 0.7 | 7/8 |
| 96 | 11 | 61 | 5 | Y | Y | Y | Y | Y | Y | Y | Y | 1.0 | 8/8 |
| 97 | 2 | 15 | 3 | Y | &middot; | Y | Y | &middot; | Y | Y | &middot; | 2.5 | 5/8 |
| 101 | 2 | 54 | 5 | Y | Y | Y | Y | Y | Y | Y | &middot; | 7.5 | 7/8 |
| 104 | 3 | 9 | 3 | Y | Y | Y | Y | &middot; | Y | Y | &middot; | 0.4 | 6/8 |
| 105 | 3 | 23 | 4 | Y | &middot; | Y | Y | Y | Y | Y | Y | 2.5 | 7/8 |
| 106 | 2 | 15 | 5 | Y | Y | &middot; | &middot; | &middot; | &middot; | Y | &middot; | &mdash; | 3/8 |
| 107 | 8 | 38 | 7 | Y | Y | Y | Y | Y | Y | Y | Y | 3.9 | 8/8 |
| 108 | 1 | 8 | 3 | &middot; | Y | &middot; | &middot; | &middot; | &middot; | Y | &middot; | &mdash; | 2/8 |
| 111 | 6 | 54 | 10 | Y | Y | Y | Y | Y | Y | Y | &middot; | 4.4 | 7/8 |
| 114 | 3 | 17 | 3 | &middot; | Y | Y | &middot; | &middot; | Y | Y | &middot; | 3.8 | 4/8 |
| 116 | 2 | 11 | 3 | Y | Y | Y | Y | Y | Y | Y | &middot; | 5.5 | 7/8 |
| 117 | 1 | 8 | 4 | Y | &middot; | Y | &middot; | &middot; | &middot; | Y | &middot; | 13.7 | 3/8 |
| 127 | 1 | 9 | 3 | Y | &middot; | Y | &middot; | Y | &middot; | Y | &middot; | 24.2 | 4/8 |
| 128 | 1 | 10 | 3 | Y | &middot; | &middot; | &middot; | &middot; | &middot; | Y | &middot; | &mdash; | 2/8 |

---

## 3. Verdicts by the pre-stated rule

Rule, fixed before the numbers: *a clause that rejects &gt; 15 % of the true Bostock families is
UNJUSTIFIED as coded; a clause that passes &ge; 85 % of true families **and** rejects &ge; 50 % of
the noise population is JUSTIFIED; anything in between is reported with a data-driven
reformulation.*

| clause | verdict by the pre-stated rule |
|---|---|
| (a) quiet-day FA &le; 0.1/tpl-day AND tremor rate &ge; 3&times; it | **UNJUSTIFIED as coded** (rejects 16% of the true families) |
| (b) tremor co-timing > same-day base rate (any-LFE hours, base 0.77) | **UNJUSTIFIED as coded** (rejects 43% of the true families) |
| (c1) located on &ge; 5 stations | **JUSTIFIED** (86% pass, rejects 67% of noise) |
| (c2) shuffled-station control passes | **UNJUSTIFIED as coded** (rejects 47% of the true families) |
| (c3) INTERIOR to the search grid | **UNJUSTIFIED as coded** (rejects 62% of the true families) |
| (c4) &#124;&Delta;z&#124; &le; 8 km of Slab2 | **UNJUSTIFIED as coded** (rejects 52% of the true families) |
| (d) &ge; 3 distinct tremor days | passes 90% but rejects only 20% of the noise population &mdash; **retain, no power** |
| (e) 2nd-half re-stack relocates within 10 km | **UNJUSTIFIED as coded** (rejects 78% of the true families) |

**One clause of the five survives: (c1), "located on &ge; 5 stations".**  Seven do not.

### 3.1 The reformulations, measured on the same populations

| reformulated clause | TRUE families (58) | their 177 templates | hidden families (18) | 18 candidate-new | noise population (1,201) |
|---|---|---|---|---|---|
| (a&prime;) quiet-day FA &le; 0.25/tpl-day AND tremor rate &ge; 3&times; it | **50/58 = 86%** | 145/177 = 82% | 18/18 = 100% | 15/18 | 778/1201 = 65% |
| (a&Prime;) tremor- vs non-tremor-day rate, day-label shuffle p &lt; 0.05 | **13/58 = 22%** | 17/177 = 10% | 5/18 = 28% | 4/18 | 397/1201 = 33% |
| (b&prime;) co-timing > base, base-rate-matched tremor hours (base 0.23) | **47/58 = 81%** | 112/177 = 63% | 15/18 = 83% | 11/18 | 360/1201 = 30% |
| (b&Prime;) same, at p &lt; 0.05 | **32/58 = 55%** | 59/177 = 33% | 11/18 = 61% | 4/18 | 126/1201 = 10% |
| (c2&prime;) control loosened to rms &lt; shuffled median | **43/58 = 74%** | 112/177 = 63% | 15/18 = 83% | 18/18 | 262/1201 = 22% |
| (c4&prime;) &#124;&Delta;z&#124; &le; 19 km of Slab2 | **43/58 = 74%** | 106/177 = 60% | 14/18 = 78% | 18/18 | 273/1201 = 23% |
| (e&prime;) FIRST-half vs SECOND-half re-stack within 10 km | **11/58 = 19%** | 15/177 = 8% | 2/18 = 11% | 0/18 | 2/1201 = 0% |

**(a) &mdash; the fixed 0.1/day cap and the 3&times; ratio.**  Sweeping the cap:

| clause (a) cap (ratio 3&times; kept throughout) | TRUE families passing | noise rejected |
|---|---|---|
| &le; 0.1/tpl-day | 49/58 = 84% | 45% |
| &le; 0.25/tpl-day | 50/58 = 86% | 35% |
| &le; 0.5/tpl-day | 51/58 = 88% | 30% |
| &le; 1.0/tpl-day | 51/58 = 88% | 28% |
| no cap, ratio only | 51/58 = 88% | 28% |

The cap is very nearly right: at 0.1/day it rejects 16 % of the true families (one point over the
bar), at 0.25/day 14 %, and loosening it further buys nothing because the residual losses are the
3&times;-ratio leg, not the cap.  But the clause only ever rejects 28&ndash;45 % of the noise
population, so it can never reach JUSTIFIED &mdash; it is a cheap sanity clause, not a
discriminator.  **The data-driven replacement proposed in the brief &mdash; "the tremor-day rate
exceeds the non-tremor-day rate at p &lt; 0.05 under a day-label shuffle" &mdash; is much worse:
it keeps 13 of 58 true families (22 %) while rejecting only 33 % of the noise.**  With 24 tremor
days, 20 quiet days and a median of ~10 detections per true-family template, the permutation test
has almost no power; the fixed threshold beats it on both axes.  Recommendation: keep the
threshold form, loosen the cap to 0.25/template-day, and do **not** adopt the shuffle.  For
scale: the only published quiet-period false-alarm number from this group (Thomas &amp; Bostock
2015, sec.4) is **2&ndash;3 false detections per day for their high-quality templates**, i.e.
20&ndash;30&times; our cap, and they treat it as acceptable.

**(b) &mdash; tremor co-timing.**  As coded the clause compares the observed in-tremor-hour
fraction with a same-day uniform redraw.  Its verdict depends entirely on the base rate of the
tremor-hour set: with an *any-LFE* hour set the 2005 base rate is 0.77 (mean 0.92 over each
template's own days) and the clause keeps 57 % of the true families; with a set matched to WA
2008's own base rate (0.23) it keeps 81 % and rejects 70 % of the noise.  Washington's audited
hour set already has base 0.24, so **(b&prime;) is the correct 2005 analogue of the Washington
clause and (b) as coded here understates it.**  81 % is still below the 85 % bar.  Adding a
p &lt; 0.05 requirement to (b) is a bad trade: 55 % of true families, 90 % noise rejection.

**(c1) &mdash; located on &ge; 5 stations.**  The only JUSTIFIED clause: 86 % of the true families,
67 % of the noise population rejected.  Note this is bed-dependent: on a 15-station bed &ge;5 is a
real bar; on agent L's 96-station 2010 bed 1,498 of 1,499 templates clear it and the clause does
nothing (sec.5.2).

**(c2) &mdash; the shuffled-station control.**  The single strongest discriminator in the whole
screen (rejects 95 % of the noise population) and the single most expensive: it rejects 47 % of
the true families as coded, and **even at its loosest useful setting** (observed rms below the
*median*, not the 5th percentile, of the 50 shuffles) it still rejects 26 %.  No threshold makes
it a justified gate.  Recommendation: **demote it from a gate to a tier.**

**(c3) &mdash; interior to the search grid.**  Rejects 62 % of the true families.  This is not a
property of a family, it is a property of **our grid**: of 597 located templates, 318 rail on the
18 km depth floor (58 of the 133 kept-family templates do), 99 on the 56 km floor, 122 on a
latitude face and 110 on a longitude face.  Bostock's families sit at 25&ndash;45 km, so a
depth-floor rail means "the moveout did not constrain depth", not "this is not a family".  The
data-driven form is not a looser clause but a **grid-adequacy diagnostic**: report the rail
fraction, and widen the grid (or report depth as unconstrained) until the fraction of *known*
families railing is under 15 %.  Keeping it as a gate discards two thirds of a real catalogue.

**(c4) &mdash; |&Delta;z| &le; 8 km of Slab2.**  The distribution of |&Delta;z| for the true
families tells the width directly.  Per-template median 9.7 km, 85th percentile 25.3 km; taking
the best-located template per family, median 6.1 km, 85th percentile 18.6 km:

| &#124;&Delta;z&#124; width | of the 50 LOCATED true families | of all 58 | noise population rejected |
|---|---|---|---|
| 5 km | 24/50 located = 48% | 24/58 = 41% | 93% |
| 8 km | 28/50 located = 56% | 28/58 = 48% | 91% |
| 10 km | 33/50 located = 66% | 33/58 = 57% | 88% |
| 12 km | 34/50 located = 68% | 34/58 = 59% | 87% |
| 15 km | 37/50 located = 74% | 37/58 = 64% | 83% |
| 18 km | 41/50 located = 82% | 41/58 = 71% | 78% |
| 19 km | 43/50 located = 86% | 43/58 = 74% | 77% |
| 20 km | 43/50 located = 86% | 43/58 = 74% | 76% |
| 25 km | 45/50 located = 90% | 45/58 = 78% | 73% |
| 30 km | 47/50 located = 94% | 47/58 = 81% | 70% |

**19 km is the width at which 85 % of the *located* true families survive**, and at that width the
clause still rejects 77 % of the noise population, so `(c4&prime;)` is justified *conditional on a
location existing*.  8 km is roughly a factor of two too tight for this 1-D locator on a
15-station bed; it is not a statement about LFEs, it is the locator's depth error.

**(d) &mdash; &ge; 3 distinct tremor days.**

| clause (d) | TRUE families passing | noise rejected |
|---|---|---|
| &ge; 2 days | 53/58 = 91% | 7% |
| &ge; 3 days | 52/58 = 90% | 20% |
| &ge; 4 days | 35/58 = 60% | 34% |
| &ge; 5 days | 25/58 = 43% | 44% |
| &ge; 6 days | 16/58 = 28% | 52% |

&ge; 3 days is the right knee (90 % of true families) but it is nearly powerless: 80 % of the noise
population also fires on &ge; 3 days.  Retain it as a cheap definitional clause ("a family
repeats"), expect nothing from it, and do not tighten it &mdash; &ge; 4 days already costs 30 % of
the true families for 14 points of noise rejection.

**(e) &mdash; the second-half re-stack.**  This is the clause that produced "0 candidates" in
Washington, and on the north bed it is catastrophic.  Re-stacking a **true Bostock family's**
template from the second half of its own detections and re-locating it moves the hypocentre by a
median of **35.5 km** (quartiles 16.3 / 35.5 / 60.7 km, n = 120 of the 177 kept-family templates
where both locations exist).  Only 13 of 58 families (22 %) reproduce within 10 km.  The purer
form &mdash; first half against second half, which removes the difference between a node-time
stack and a detection stack &mdash; is worse still: median 33.9 km, 11 of 58 families (19 %).
Even a 50 km bar keeps only 62 % of them.  **Location stability under re-stacking is not a
property that Bostock's own families have on this bed**; the 1-D locator plus a 12-s stack of ~20
detections simply does not reproduce to 10 km.  The clause must be dropped.  What the split-half
displacement *is* good for is an honest **location-uncertainty estimate** to report beside each
candidate.

---

## 4. What the literature actually required

**Access note.** Bostock, Royer, Hearn &amp; Peacock (2012), G-cubed 13, doi 10.1029/2012GC004391
is genuinely closed (Unpaywall `is_oa: false, oa_status: closed, oa_locations: []`; OpenAlex and
Semantic Scholar agree; Wiley 403; Royer's UBC thesis, the obvious long restatement, is behind a
UBC bot block).  Everything below is quoted from sources that were actually read, and each is
labelled as the 2012 paper itself (never), a later paper by the same group describing the 2012
catalogue, or a sister application of the same chain.

**Bostock et al. (2015), JGR 120, 6329** (PDF, `pages.uoregon.edu/amthomas/pubs/bostock_jgr_2015.pdf`).
The catalogue is inherited, not re-screened: "This catalog comprises those templates from *Bostock
et al.* [2012] and *Royer and Bostock* [2013], as well as new additions for a total of 130 distinct
source regions on the plate boundary."  The **only culling stated in the paper** is the
declustering rule: the catalogue "has been culled such that only the detection with the strongest
correlation in any given 4 s interval &hellip; is retained, resulting in a total of 269,586
detections across all years."  The two threshold rules in sec.2 are **per channel**, for the
magnitude work, not family acceptance: a station/channel is admitted "if (i) the correlation
coefficient of the template channel with the final template stack exceeds 0.4 and (ii) the energy
within the detection window exceeds 70 % of that in the total (6 s long) channel trace", and an
individual detection's channel is used if its cc with the template channel exceeds 0.25.  Arrival
clarity appears only as a *property* of the result: "templates reveal particularly simple
waveforms: dipolar pulses in particle velocity and single pulses in particle displacement; hence,
the direct wave is easily identified and consistent across these stations."  Tremor association is
structural &mdash; "Tremor data from these episodes enabled the construction of an initial set of
high-SNR LFE templates."

**Royer, Thomas &amp; Bostock (2015), JGR 120, doi 10.1002/2014JB011430** is the one reachable
source that quantifies a criterion applied to the 2012 catalogue: it "was culled to **93
geographically independent LFE families** with a total of 148,847 detections from ETS episodes
between 2003 and 2013" (and the northern Washington catalogue "was regrouped into 100
geographically independent LFE families").  130 templates &rarr; 93 independent families, so the
*dedupe* was real and non-trivial &mdash; and it is the **only** family-acceptance rule positively
attested for Bostock 2012.  The threshold is not given; the paper defers method to the 2012 paper.

**Sister applications of the same chain** supply the thresholds the 2012 paper withholds.  Plourde,
Bostock, Audet &amp; Thomas (2015), GRL 42, 4849: "Final templates were tested for uniqueness using
two distinct methods &hellip; the templates were considered to originate from the same family if
&gt; 10 % of detections were shared &hellip; A correlation coefficient &gt; 0.3 was considered to
suggest templates originated from the same family"; and a *de facto* pickability gate, "A total of
59 initial templates from Cascadia produced waveforms sufficiently clear to be picked.  These 59
templates were found to represent 34 distinct LFE families."  Thomas &amp; Bostock (2015),
*Tectonophysics*, is the most explicit acceptance recipe in the reachable literature: "We require a
minimum of **100 detections** in the 3-day window to stack an LFE template &hellip; We eliminate
duplicate families by ensuring that no family shares more than 10 % of detection times, defined as
a difference between detection times of [&lt;] 5 s, with any other family.  Finally, the remaining
LFE families are located using HYPOINVERSE."  Two of its statements bear directly on our clauses.
**Location was not a gate**: "Of these 18, we were unable to locate five families because the
templates were not of sufficient quality to accurately pick arrivals" &mdash; and those five were
kept.  And the group's own **quiet-period false-alarm estimate**: "Assuming, as a worst case
scenario, that the detection rate prior to the arrival of the SSE in Oregon is representative of
the false detection rate, relatively high-quality templates &hellip; have **2&ndash;3 false
detections per day**, while over 700 events over a 3-day period were used in template
construction."  **Their own high-quality templates sit 20&ndash;30&times; above our 0.1/day cap,
and they treat that as acceptable.**  Finally, Frank et al. (2014) place the 2012 catalogue
explicitly in the *automatic* camp &mdash; templates found "with automatic algorithms applied to
previously identified tremor [&hellip; Bostock et al., 2012 &hellip;]" &mdash; not the
visual-inspection lineage, and Rubin &amp; Armbruster (2013) note that 2012 itself acknowledged
that "some time windows pass the detection threshold for more than one LFE family", i.e. families
were not required to be mutually exclusive at the detection level.

**Chestler (2017), UW dissertation sec.2.1&ndash;2.2** (open access; the same text as Chestler &amp;
Creager 2017, JGR 122).  Seeds are peaks "larger than 18 times the median absolute deviation" on
each Ghosh et al. (2012) array beam, boxcar-ed at &plusmn;1 s and stacked across arrays.  A
candidate is verified by autocorrelating 15-s 2&ndash;8 Hz 3-C windows on "a set of 3-stations from
each array", combining them with "an nth-root stack where n = 2", and looking "for times when the
stacked autocorrelation function exceeds 45 times the median absolute deviation"; "up to 150 of the
time windows that correlate best with the initial detection are stacked to create a template
waveform".  Then the one explicitly qualitative gate in the published Cascadia literature: **"If
the template has clear P- and S-wave arrivals, we strengthen it by autocorrelating the template
waveform through time to look for additional correlating time windows and re-stack to create
template waveforms with higher SNRs"**, iterated three times, with membership set by "a stricter
correlation threshold that varies from 0.05 to 0.20 depending on the number of stations that are
active at the time of the LFE".  Families are located by grid search on the Preston et al. (2003)
3-D model from P and S differential times, refined by Waldhauser &amp; Ellsworth double difference.
Duplicates are removed by R_com: "When two LFE families have R_com &gt; 0.5, the family with the
lower median correlation coefficient is removed.  This is nearly always the family with fewer
detections.  This reduced the number of LFE families from 57 to 43."  **No quiet-day false-alarm
rate, no minimum day count, no Slab2 agreement, no interior-to-grid clause, no split-half
location-stability test.**  ETS association is structural (the AofA beams are tremor beams), never
a statistic.

### 4.1 Clause by clause: codification or addition?

| our clause | published counterpart | verdict |
|---|---|---|
| (c1) located on &ge; 5 stations | **partial codification.** Both catalogues locate their families (Bostock 2012 via Royer et al. 2015; Chestler 3-D grid search + double difference) &mdash; but Thomas &amp; Bostock (2015) explicitly **kept 5 of 18 families they could not locate.** Our shuffled-station control has no published counterpart. | keep as a gate; keep unlocated templates on a reported side list |
| (c2) shuffled-station control | **addition.** No published location control of any kind. | demote to a tier |
| (b) tremor co-timing above a base rate | **weak codification.** Both groups get tremor association for free by searching tremor/ETS data; neither tests it statistically. | keep as a tier, base-rate-matched |
| (a) quiet-day false-alarm rate | **addition, and far stricter than the group's own practice.** The only published quiet-period estimate is Thomas &amp; Bostock (2015)'s 2&ndash;3 false detections per day for *high-quality* templates, used post hoc and not as a gate. | keep, loosened, as a sanity clause; never as the deciding clause |
| (d) &ge; 3 distinct tremor days | **addition.** No day-count rule exists anywhere; the nearest relative is Thomas &amp; Bostock's &ge; 100 detections in a 3-day window. | keep, no power |
| (c4) &#124;&Delta;z&#124; &le; 8 km of Slab2 | **addition.** Interface proximity is always reported as a *result*, never used as a filter; Slab2 postdates both catalogues. | reformulate to 19 km |
| (c3) interior to the search grid | **addition**, and a property of our grid, not of a family. | drop; report the rail fraction |
| (e) second-half re-stack within 10 km | **addition.** Nothing resembling a split-half or reproducibility test appears in any source read; "iterate until convergence" is an alignment criterion, not an acceptance test. | drop |
| &mdash; | **published, not implemented by us:** the *dedupe* (the one rule positively attested for Bostock 2012 &mdash; 130 templates culled to 93 geographically independent families; Plourde &gt;10 % shared detections or cc &gt; 0.3; Chestler R_com &gt; 0.5, 57&rarr;43), and Chestler's "clear P- and S-wave arrivals". | see sec.6 |

Six of our eight coded conditions have no counterpart in any published Cascadia LFE catalogue,
and the two most destructive of them &mdash; (c3) and (e) &mdash; are inventions of this project.
Conversely, the one criterion that **is** positively attested for Bostock 2012, a duplicate merge,
is the one we do not apply (`REPRO_BOSTOCK.md` sec.7 measured that neither published uniqueness
test separates these families, so the omission is itself evidence-based &mdash; but it means our
"58 families / 177 templates" is a template count, not an independent-family count).

---

## 5. The recommended screen, and what Washington looks like under it

### 5.1 The screen

**Gates (tier A)** &mdash; all four required:

1. **(c1)** located on &ge; 5 stations from the template's own stack picks.  *JUSTIFIED: 86 % of
   the true families, rejects 67 % of the noise population.*
2. **(c4&prime;)** |&Delta;z| &le; **19 km** of Slab2, applied to located templates.  *86 % of the
   located true families, rejects 77 % of the noise population.*  The width is the 85th percentile
   of |&Delta;z| for known families under this locator and must be re-measured for any other
   locator or bed.
3. **(d)** &ge; **3** distinct tremor days.  *90 %, almost no power &mdash; a definitional clause.*
4. **(a&prime;)** quiet-day false alarms &le; **0.25** per template-day **and** the tremor-day rate
   &ge; 3 &times; that.  *86 %, 35 % noise rejection.*  Not the day-label shuffle, which was
   measured to be far worse.

**Tier B** = tier A + the shuffled-station control.  **Tier C** = tier B + tremor co-timing above a
**base-rate-matched** tremor-hour set.  Both are reported as tiers because each rejects &gt; 15 %
of Bostock's own families; neither may be used to turn a positive count into zero.

**Dropped:** (c3) interior-to-grid &rarr; replaced by reporting the fraction of located templates
railing on each grid face, with the grid widened until known families stop railing.  (e)
second-half stability &rarr; replaced by reporting the split-half location displacement as the
candidate's location uncertainty.

On the north bed this screen recovers **33 of the 58 true families at tier A, 21 at tier B, 13 at
tier C**, against **11.0 % / 1.3 % / 0.8 %** of the 1,201-template noise population &mdash; versus
**2 of 58** and 0 % for the screen as coded.

### 5.2 Washington 2008 (agent N) re-screened

`results/n_puget2008/2008{,_null}/screen_final.csv` + `grow/det_it0.parquet` + `det_fa.parquet`,
unchanged; re-screened by `scripts/q_sv_wa.py` &rarr; `results/q_screenval/wa2008_*_rescreen.csv`.

| | real (160 templates) | site-shift null (36) |
|---|---|---|
| **screen as coded (a)&ndash;(e)** | **0** | **0** |
| (c1) located &ge; 5 stations | 156 | 36 |
| (c4&prime;) &#124;&Delta;z&#124; &le; 19 km | 62 | 15 |
| (d) &ge; 3 days | 110 | 17 |
| (a&prime;) quiet FA &le; 0.25/tpl-day, &ge; 3&times; | 70 | 9 |
| **recommended tier A** | **17** | **2** |
| **tier B** (+ control) | **4** | **0** |
| **tier C** (+ co-timing) | **4** | **0** |

The 2008 bed is not empty under a screen calibrated on known families: **17 real against 2 null at
tier A (10.6 % of real templates against 5.6 % of null templates), and 4 real against 0 null at
tier B/C.**  Four templates is not a discovery &mdash; the null has only 36 templates, so 0 there
is weak evidence &mdash; but it is a list to re-fire on 2007, which "0 candidates" was not.  The
pair-stage diagnosis of `WA_2008_RESULTS.md` sec.17 (tremor/quiet pair ratio 0.89) is untouched by
this and remains the real problem with that bed.

### 5.3 Washington 2010 (agent L) re-screened

`results/l_wa2010{,_null}/verify_it0.csv`.  Clause (a) cannot be evaluated: the 2010 grown
templates were never re-fired on quiet days, so (a&prime;) is treated as passing and tier A is
only (c1)+(c4&prime;)+(d).

| | real (1,499 templates) | site-shift null (1,500) | real/null |
|---|---|---|---|
| **L's screen as published** (located, control, &#124;&Delta;z&#124; &le; 8, &ge; 3 days, co-timing) | **50** | **17** | 2.9 |
| (c1) located &ge; 5 stations | 1498 | 1500 | 1.0 |
| (c4&prime;) &#124;&Delta;z&#124; &le; 19 km | 891 | 876 | 1.0 |
| (d) &ge; 3 days | 1376 | 1269 | 1.1 |
| **recommended tier A** (no quiet-day data) | **825** | **740** | 1.1 |
| **tier B** (+ control) | **234** | **59** | 4.0 |
| **tier C** (+ co-timing) | **127** | **27** | 4.7 |

Two things follow.  (i) On a 96-station bed the location gates are free &mdash; 1,498 of 1,499
templates locate on &ge; 5 stations and |&Delta;z| &le; 19 km passes 59 % of the null &mdash; so
tier A is meaningless there and **the recommended thresholds must be re-derived per bed**, which
is exactly what sec.5.1 says about (c4&prime;).  (ii) The real-vs-null **contrast improves**:
2.9&times; under L's screen, 4.0&times; at tier B, 4.7&times; at tier C, with 127 real against 27
null.  The 2010 bed's problem is not the screen, it is that the null produces hundreds of
screen-passing templates at all &mdash; the pair-stage failure diagnosed in PLAN.md sec.18.
**Clause (a) is the one clause that would help there and it has never been run on that bed: the
2010 templates should be re-fired on the 40 hour-matched quiet hours.**

---

## 6. Caveats, and what this does not say

* The "noise population" is 1,201 templates that do not locate or do not match a catalogued
  family.  It is *not* certified noise: some of it is real LFE signal from families Bostock never
  listed, and `REPRO_BOSTOCK.md` sec.11.3 argues seven of the candidate-new templates are probably
  unlisted neighbours.  Noise-rejection fractions are therefore upper bounds on the true
  specificity.
* The 2005 quiet days carry fewer stations (7&ndash;12) than the ETS days (9&ndash;15).  A thinner
  network makes the per-day 8 &times; MAD bar noisier, which inflates the quiet-day false-alarm
  rate.  The Washington comparison has the same confound in the same direction.
* Clause (b) has no true tremor catalogue for 2005; the substitute is Bostock's own LFE catalogue,
  leave-one-family-out.  A family's own detections cannot inflate its score, but the hour set is
  still built from the same physical process the templates detect.
* The 85 % / 50 % rule is a choice, fixed in advance.  A clause that keeps 60 % of true families
  is not wrong in itself &mdash; a high-purity catalogue may be worth it &mdash; but it must not
  then be used to argue that a bed has nothing on it.
* Two published practices this project does **not** implement are Chestler's "clear P- and S-wave
  arrivals" gate and the family dedupe (Bostock 2012's cull to 93 geographically independent
  families, Plourde's &gt;10 %-shared-detections / cc &gt; 0.3, Chestler's R_com &gt; 0.5).  The
  arrival-clarity gate is the only qualitative screen either group applied and is worth codifying
  (e.g. a stack-envelope P and S peak test).  The dedupe was measured on this bed
  (`REPRO_BOSTOCK.md` sec.7) and neither published test separates these families &mdash; so our
  counts are template counts, and duplicate candidates must be resolved by location.

## 7. Files

| what | where |
|---|---|
| per-template screen, all clauses and reformulations | `results/q_screenval/screen_north.csv` |
| per-family table, 58 rediscovered Bostock families | `results/q_screenval/per_family_58.csv` |
| markdown fragments behind every table here | `results/q_screenval/fragments.json` |
| quiet-day re-fire of the 2005 templates | `results/q_screenval/det_quiet15.parquet`, `quiet_day_stations.csv` |
| first/second-half re-stacks and re-locations | `results/q_screenval/reloc_h1.csv`, `reloc_h2.csv` (+ `_tpl.npz`) |
| Washington re-screens | `results/q_screenval/wa2008_{real,null}_rescreen.csv`, `wa2010_{real,null}_rescreen.csv`, `wa_rescreen.json` |
| code | `scripts/q_sv_data.py`, `q_sv_screen.py`, `q_sv_tables.py`, `q_sv_wa.py`, queue `scripts/q_sv_queue.sh` |
| logs | `logs/q_sv_queue.log`, `q_sv_quiet.log`, `q_sv_h1.log`, `q_sv_h2.log` |
