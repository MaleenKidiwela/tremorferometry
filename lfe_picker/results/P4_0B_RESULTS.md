# P4-0b — the TRUE Arm B kill test on the PQ-era north

Run 2026-09-20. Pre-registration `results/p4/p4_0b_freeze.json` (written and hashed before any
pass, sha256 `28fe6d5d02a8cc1c…`). Spec: `results/merlin_review_p4_2026-09-20.md` §P4-0b and
P4-2 steps 2–6; construction: `results/p3d/arm_b_preregistration.md` steps 2–5, executed here for
the first time. Ablation pre-registration `results/p4/p4_0b_abl_prereg.json`.

---

## §0 Verdict (frozen, pre-registered)

> **ARM (i) — GIVEN seeds: 0 of 18 eligible kept families assembled. Bar ≥ 2/3 (≥ 12). → FAIL.**
> **ARM (ii) — TEMPLATE-FREE seeds: 1 kept family rediscovered (family 44, by template `F187`,
> 34 one-to-one ±2 s matches = 36 % of its 95 hp events, 0 aliases, null max 2). Bar ≥ 3. → FAIL.**
> **Pre-registered consequence: STOP P4-2 — do not run Puget discovery.**
> No family was NO-REFERENCE: all 18 were testable, so the 0/18 is a measured zero, not an
> untestable one. The five pass-0 survivors all produced coherent iterated stacks (per-site S-SNR
> median 11.95 dB, alignment-free S-excess ≥ +3 dB at 10/10 sites for all five) — the construction
> ran and produced real waveforms; it did not reach the pre-registered bar.

### Documented deviations from Merlin's letter (all in the freeze, all decided before the run)

1. **Region-centre depth.** Slab2 at the box centre (48.50 N, −123.75 E) returns **21.9 km**, about
   13 km above the Bostock family field. The moveout depth used is the **median family depth
   35.27 km**; lat/lon are Merlin's `p_common.GRID_LAT/GRID_LON`. Using Slab2 would have mis-cut
   every template by ≈ 3 s in the common term.
2. **`net_mad` conversion.** Merlin's equal-weight site sum (per-site mean of member-station
   lowhalf-MAD z, clipped at 8) is converted to `net_mad` by the frozen runner's per-(template,
   day) median/MAD normalisation. Without it one fixed per-template threshold is not comparable
   across days carrying 5, 6 or 7 sites, and per-day availability is deliberately **not**
   intersected. Consequence to keep in mind: the z-clip bounds the sum at 8 × n_sites, so `net_mad`
   saturates near 20 on a 6-site day.
3. **Arm-(i) reference = truth-stack, not transferred template.** Merlin's wording is "the family's
   own transferred template". The P3D transfer bank on disk covers only 26–184 (family, station)
   pairs and **none at PGC/VGZ/ALBHB**, so a 7-site transferred reference does not exist for these
   18 families. The reference is each family's own **truth-stack** 12-s template (`align_stack`
   over *all* its hp events at each site, cut from S−2 s) grown through the identical machinery.
   It uses more events than any transferred template would, so the relative bar is **conservative**;
   absolute recall/purity are reported alongside.
4. **Arm-(ii) one-offset matching.** A template-free seed's origin axis is the region-centre
   moveout, so a true family's hits arrive displaced by **one** constant `b = t_seed − τ₀`,
   identical for every event (derivation in the freeze). `b` is estimated as a single nuisance
   parameter per (template, family) over ±30 s and the search is **charged to a 20-realisation null**
   in which the template's hits are redrawn uniformly inside their own days. Strict `b = 0` numbers
   are reported too — for `F187`/family 44 they are 33 vs 34 matches, so the offset bought almost
   nothing.

---

## §1 Setting

| item | value |
|---|---|
| sites (7) | PGCq (PGC+B009+B010+B011), VGZ, SYMB, BUCKS, LHLYB, ALBHB, DAOB |
| never in any sum | B926, NLLB, LRIV (holdouts.json) + CLRS, GOBB, LOPEZ, TURTL (P3D withheld sites) |
| episodes | S18 (2024-07-14…09-10), S21 (2025-08-04…10-24), S22 (2026-04-01…06-01) |
| episode days | **203 of 203 usable** (≥ 3 sites); site counts **5 sites: 21 d, 6 sites: 53 d, 7 sites: 129 d** — days are per-station, *not* intersected |
| quiet days | **20**, 2024-09-06 … 2026-06-24, all 7 sites, no PNSN tremor within 150 km of the centre |
| hidden / sealed | 30 hidden families never seeded, scored or matched; 2011 and 2013 sealed; no Puget byte read or written; no FDSN request |

**Eligible kept families (18)** — non-hidden, ≥ 50 hp events (`realness_v1_1`, `p_real ≥ 0.8`) on
usable episode days. bkey (hp count): 13 (51), 21 (60), 24 (53), 32 (68), 36 (59), 39 (55),
44 (95), 52 (79), 56 (96), 61 (130), 75 (50), 81 (72), 86 (65), 93 (52), 98 (115), 110 (62),
114 (215), 116 (55). Bar (i) is therefore ≥ 12/18.

**Threshold rule (frozen FA-only, `c_tiers_by_epoch.py`).** Smallest `T` on a 0.25 grid at or above
the 8.0 floor with ≤ 0.1 false hits per template-day on the 20 quiet days. No template saturated.
`T_seed` actually used:

| bank | n | median | range | fraction at the 8.0 floor |
|---|---|---|---|---|
| given pass-0 | 54 | 8.00 | 8.00–11.25 | 0.963 |
| given N-moveout | 54 | 8.00 | 8.00–8.00 | 1.000 |
| given it1 / it2 | 5 / 5 | 8.00 / 8.00 | 8.00–13.50 / 8.00–15.75 | 0.80 / 0.80 |
| truth-stack reference | 18 | 8.00 | 8.00–8.00 | 1.000 |
| free pass-0 | 1500 | 8.00 | 8.00–20.25 | 0.992 |
| free N-moveout | 1500 | 8.00 | 8.00–20.25 | 0.993 |
| free it1 / it2 | 2 / 2 | 8.00 / 8.50 | 8.00–8.00 / 8.00–9.00 | 1.00 / 0.50 |

The quiet-day false-alarm rate at `net_mad ≥ 8` is 0 for almost every template, so the FA-only rule
returns the floor — the growth threshold is set by the recording floor, not by measured noise.

**Compute.** Frozen queue **3.04 h** wall (pid 1040534), ablation **0.64 h** (pid 1046043), both on
one shared H100. Measured growth rate **≈ 940 template-station-days/s** (free pass-0: 1,500
templates × 203 days × ≈ 9.5 stations in 51.4 min). Longest single steps: free pass-0 growth
51.4 min, free N-moveout growth 46.8 min, given pass-0 growth 14.3 min, reference bank build 8.4 min.

---

## §2 Arm (i) — given seeds: mechanics

54 seeds (3 hp events per eligible family, distinct days, distinct episodes preferred, ranked by the
frozen catalogue's own `net_mad`), 20-s templates cut from predicted S − 6 s at each family's own
Bostock cell.

### Pass-0 survival: 54 → 5

| clause | templates passing it (of 54) | templates failing it |
|---|---|---|
| ≥ 20 independent hits | **5** | 49 |
| ≥ 3 distinct days | 36 | 18 |
| ≥ 2 episodes | 39 | 15 |
| ≤ 20 % of hits within ±10 min of the seed | 32 | 22 |
| N-moveout ratio ≥ 3 | 23 | 31 |
| **all five (survivors)** | **5** | 49 |

The independent-hit floor is the sole binding clause: applying it first takes 54 → 5, and the other
four clauses then remove nothing. Pass-0 independent-hit counts are **median 6, p75 11.5, max 34**
over the whole 203-day record.

**The five survivors** (all from distinct families — 32, 44, 61, 75, 110):

| tid | indep | days | episodes | frac within ±10 min of seed | N-moveout hits | ratio |
|---|---|---|---|---|---|---|
| G32_1 | 29 | 26 | 2 | 0.034 | 3 | 9.7 |
| G44_0 | 34 | 8 | 2 | 0.029 | 1 | 34.0 |
| G61_0 | 22 | 5 | 2 | 0.091 | 2 | 11.0 |
| G75_1 | 32 | 11 | 2 | 0.063 | 1 | 32.0 |
| G110_1 | 23 | 7 | 3 | 0.043 | 7 | 3.3 |

**Burst convergence: absent.** The fraction of hits within ±10 min of the seed is 0.011–0.065 per
family (median 0.027); **0 of 5 families exceed 50 %**. Merlin's "report immediately" failure
signature did not fire — the decluster and independence guards work.

### Iteration trajectories (exactly two iterations, no third)

`align_stack` (2 iterations, ±1.5 s) → re-pick S at the highest-S-SNR site → one global shift → 12-s
templates from S−2 s. Global picks moved by −1.35…+0.50 s at it1 and −1.65…+0.25 s at it2 (i.e. the
assumed cell was already close). Stack quality rose monotonically: per-site S-SNR median
**9.85 → 11.95 dB**, S-excess median **9.64 → 12.15 dB**, all 10 stations registering. The cheap
alignment-free screen passes for **5/5** templates at **10/10 sites** — which is exactly why Merlin
demoted it from evidence to screen (§F3).

Hit counts also rose: G44_0 34 → 103 → 114, G75_1 32 → 131 → 147, G61_0 22 → 67 → 92, G110_1
23 → 79 → 92, G32_1 29 → 37 → 35.

### Scoring against the truth-stack reference (`score_given.csv`)

| bkey | n_hp | ref indep | recall_ref | purity_ref | need recall (0.5×) | need purity (0.8×) | grown tid | grown hits | recall | purity | assembled |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 13 | 51 | 55 | 0.627 | 0.582 | 0.314 | 0.465 | — | 0 | 0 | — | no |
| 21 | 60 | 74 | 0.550 | 0.446 | 0.275 | 0.357 | — | 0 | 0 | — | no |
| 24 | 53 | 45 | 0.623 | 0.733 | 0.311 | 0.587 | — | 0 | 0 | — | no |
| **32** | 68 | 17 | 0.206 | 0.824 | 0.103 | 0.659 | G32_1 | 35 | 0.015 | 0.029 | no |
| 36 | 59 | 47 | 0.458 | 0.574 | 0.229 | 0.460 | — | 0 | 0 | — | no |
| 39 | 55 | 11 | 0.109 | 0.545 | 0.055 | 0.436 | — | 0 | 0 | — | no |
| **44** | 95 | 101 | 0.453 | 0.426 | 0.226 | **0.341** | G44_0 | 114 | **0.358** | **0.298** | **no** |
| 52 | 79 | 42 | 0.342 | 0.643 | 0.171 | 0.514 | — | 0 | 0 | — | no |
| 56 | 96 | 100 | 0.417 | 0.400 | 0.208 | 0.320 | — | 0 | 0 | — | no |
| **61** | 130 | 66 | 0.254 | 0.500 | 0.127 | 0.400 | G61_0 | 92 | 0.092 | 0.130 | no |
| **75** | 50 | 116 | 0.600 | 0.259 | 0.300 | **0.207** | G75_1 | 147 | **0.360** | **0.122** | **no** |
| 81 | 72 | 70 | 0.431 | 0.443 | 0.215 | 0.354 | — | 0 | 0 | — | no |
| 86 | 65 | 52 | 0.446 | 0.558 | 0.223 | 0.446 | — | 0 | 0 | — | no |
| 93 | 52 | 104 | 0.615 | 0.308 | 0.308 | 0.246 | — | 0 | 0 | — | no |
| 98 | 115 | 119 | 0.496 | 0.479 | 0.248 | 0.383 | — | 0 | 0 | — | no |
| **110** | 62 | 98 | 0.726 | 0.459 | 0.363 | 0.367 | G110_1 | 92 | 0.274 | 0.185 | no |
| 114 | 215 | 80 | 0.233 | 0.625 | 0.116 | 0.500 | — | 0 | 0 | — | no |
| 116 | 55 | 70 | 0.455 | 0.357 | 0.227 | 0.286 | — | 0 | 0 | — | no |

**The two near-passes, explicitly.**

* **Family 44 / G44_0.** Recall **0.358** against a reference recall of 0.453 and a required 0.226 —
  it clears the recall clause by 58 %. It fails on **purity 0.298 against the required 0.341**
  (reference purity 0.426). It is short by 0.043 in purity and by nothing else.
* **Family 75 / G75_1.** Recall **0.360** against a reference recall of 0.600 and a required 0.300 —
  clears recall. Fails on **purity 0.122 against the required 0.207** (reference purity 0.259).

Both near-passes fail the same way: after two iterations the template detects a fair share of its
family's events and a larger share of things that are not that family's events. G32_1 is the
opposite case — it passed every survival clause but has purity 0.029 against its nominal family, so
it was never that family's template at all.

---

## §3 Sensitivity ablation — **NOT THE FROZEN RESULT**

Pre-registered separately (`results/p4/p4_0b_abl_prereg.json`; all 10 frozen `v_` script hashes
verified unchanged; outputs isolated in `results/p4/p4_0b_abl/`). The frozen verdict of §0 stands
and nothing here modifies it. The question: the ≥ 20-independent-hit floor was set for a multi-year
record and applied to 203 days, where it removed 49 of 54 seeds before any template could iterate.
Seeds, pass-0 templates, pass-0 and N-moveout hit tables, the T rule, the two-iteration count, the
dedupe rule, the relative clauses and the truth-stack reference are all the frozen ones.

| | floor | days | pass-0 survivors | families reaching iteration | **assembled at it2** | near-misses |
|---|---|---|---|---|---|---|
| frozen | ≥ 20 | ≥ 3 | 5 | 5 | **0 / 18** | 2 |
| **A1** | ≥ 8 | ≥ 3 | 18 | **11** | **1 / 18** (family 110) | 5 |
| **A2** | ≥ 5 | ≥ 2 | 21 | **12** | **1 / 18** (family 110) | 6 |

Lowering the floor roughly doubles the families that reach iteration (5 → 11 → 12) and moves the
score from 0/18 to 1/18. The one family that assembles is **110**, and it assembles through a
*different* seed than the frozen survivor (G110_0, recall 0.387 ≥ 0.363 needed, purity 0.453 ≥
0.367 needed) — the frozen run's G110_1 reached only purity 0.185. A1 near-misses (recall clause
met, purity clause failed): families 21, 44, 75, 81, 116. A2 adds 61. Every single near-miss fails
on purity, none on recall.

### A3 — purity versus hit count (diagnostics only, no GPU)

Fraction of a pass-0 template's independent hits that match its **own** family's hp events
(one-to-one ±2 s), over all 54 pass-0 templates:

| independent hits | templates | median own-family purity |
|---|---|---|
| 0–4 | 21 | **1.000** |
| 5–7 | 9 | 0.714 |
| 8–12 | 13 | 0.500 |
| 13–19 | 6 | 0.528 |
| ≥ 20 | 5 | **0.364** |

Spearman(purity, n_indep) = **−0.68**. And the per-iteration trajectory of the five frozen
survivors shows the same trade *within* each template — recall rises, purity falls, monotonically:

| tid | indep 0→1→2 | purity 0→1→2 | recall 0→1→2 |
|---|---|---|---|
| G44_0 | 34 → 103 → 114 | 0.471 → 0.311 → 0.298 | 0.168 → 0.337 → 0.358 |
| G75_1 | 32 → 131 → 147 | 0.281 → 0.145 → 0.122 | 0.180 → 0.380 → 0.360 |
| G110_1 | 23 → 79 → 92 | 0.391 → 0.215 → 0.185 | 0.145 → 0.274 → 0.274 |
| G61_0 | 22 → 67 → 92 | 0.364 → 0.149 → 0.130 | 0.062 → 0.077 → 0.092 |
| G32_1 | 29 → 37 → 35 | 0.034 → 0.027 → 0.029 | 0.015 → 0.015 → 0.015 |

**What it means.** The seeds that are pure are the ones that are starved, and they do not become
pure families after two iterations: every template that grows enough hits to clear a survival floor
has already admitted a majority of non-family hits, and iterating on that hit set makes the
admission worse, not better. The purity a template ends with is set at pass 0 and degrades with
growth; it is not recovered by stacking.

### A4 — the reference at three thresholds

The frozen reference recall spans 0.11–0.73, so the reference itself was checked at the frozen
operating point's tiers as well as at the 8.0 floor the run used:

| threshold | median reference recall | median reference purity | median independent hits |
|---|---|---|---|
| T = 8.00 (frozen P4-0b FA-only value) | 0.454 | 0.489 | 70 |
| T_ext = 8.25 | 0.448 | 0.515 | 60 |
| T_core = 8.75 | 0.398 | 0.565 | 48.5 |

Raising the threshold trades reference recall for reference purity as expected and moves both by
under 0.08 — the 0.11–0.73 spread is a property of the *families* (family 39 recall 0.109, family 32
0.206, family 110 0.726 at T = 8.0), not an artefact of the 8.0 floor. Since the relative clauses
scale with the reference, a stricter reference would have *lowered* the recall bar and *raised* the
purity bar — the direction the near-passes were already failing in.

### Why the floor was not the cause

The floor determined how many templates were *eligible* to be scored, not whether a scored template
met the bar. Dropping it from 20 to 8 and then to 5 more than doubled the families reaching
iteration (5 → 11 → 12) but moved the score only from 0/18 to 1/18, and every added family failed on
the purity clause, never on recall (A1 near-misses 21/44/75/81/116; A2 adds 61). A3 supplies the
mechanism: purity and hit count are anti-correlated at −0.68 across pass-0 templates and fall
monotonically within every survivor across both iterations. A lower floor therefore admits *more
templates of lower purity*, which is the axis the bar is failing on. The 0/18 is a purity result,
not a floor result — the floor changes 0/18 to 1/18 and no lower floor can change it further,
because below 5 independent hits the N-moveout ratio can no longer be evaluated at all.

---

## §4 Arm (ii) — template-free seeds

### Seed statistics (the corrected rule worked)

| stage | count |
|---|---|
| PNSN tremor hours (≥ 5 epicentres within 50 km of the centre, episode days) | 883 |
| candidate 12-s windows at 6-s hop | 529,800 (529,695 scored) |
| in band (≥ 3 sites at z ≥ 3 **and** 4 ≤ z_med ≤ 25) | **48,730** |
| after ±60 s earthquake veto (`trackC/eq_catalog.csv`, 26,465 events) | 48,220 |
| after ≤ 3 seeds/hour and ≥ 120 s separation | **1,590** |
| after HF ratio E(10–19)/E(2–8) < 3.5 at LHLYB and SYMB (native 100 Hz) | 1,590 of 1,590 — **no window vetoed** |
| budget, ordered by z_med | **1,500** (z_med 4.54–25.00) |

**The cap worked.** Max **3 seeds per hour** (the cap), **top-10-hours fraction 2.0 %** against
Merlin's >20 % failure signature, **609 distinct hours** over **59 distinct days**, top-5 days
16.1 %. This is the direct repair of P3D's F2 (456 seeds in one hour, top-10 hours holding 21 % of
9,874 seeds). The HF veto never fired — every in-band candidate had HF ratio ≈ 0.1–0.3, two orders
below the 3.5 cut, so on this network the HF clause is inert.

### Survival: 1,500 → 2

| clause | templates passing (of 1,500) |
|---|---|
| ≥ 20 independent hits | **19** |
| ≥ 3 distinct days | 165 |
| ≥ 2 episodes | 224 |
| N-moveout ratio ≥ 3 + near-seed clause, with the above | **2** |

| tid | seed day / episode | seed z_med | indep | days | eps | near-seed | N-moveout | ratio |
|---|---|---|---|---|---|---|---|---|
| F187 | 2025-10-05 20 h, S21 | 22.91 | 22 | 6 | 2 | 0.091 | 6 | 3.67 |
| F1184 | 2024-08-25 23 h, S18 | 6.54 | 20 | 17 | 3 | 0.200 | 5 | 4.00 |

Both iterated twice, both survived, both distinct under dedupe (0 merged pairs), both pass the
cheap screen — F187 at 9/9 registering sites (S-excess +5.5 to +14.9 dB), F1184 at 9/10 (+1.6 dB at
BUCKS, +7.7 to +19.1 dB elsewhere). After two iterations F187 has 95 independent hits on 7 days /
2 episodes; F1184 has 24 on 18 days / 3 episodes.

### The one rediscovery

**F187 → family 44.** 34 one-to-one ±2 s matches of its 95 independent hits to family 44's 95 hp
events = **35.8 % of the family** (bar: ≥ 30 % and ≥ 20 matches). Fitted constant offset
**b = −0.65 s**; under **strict b = 0** it is **33 matches / 34.7 %**, so the one nuisance parameter
bought a single match. The 20-realisation within-day null gives **mean 0.60, max 2** matches —
34 against a null max of 2. **0 aliases**: its next-best kept family is 36 with 9 matches (15.2 %,
below both clauses) and then 32 with 4. F187's seed is the loudest class of in-band window
(z_med 22.9, near the top of the 4–25 band).

**F1184 matched no kept family.** Its best score against any of the 18 is **1 match** (family 44,
1.1 %), with null max 1 — i.e. indistinguishable from chance. It has 24 independent hits spread over
**18 distinct days and 3 episodes**, registers at all 10 stations with +7.7 to +19.1 dB
alignment-free S excess at 9 of them, and passed the N-moveout scramble at ratio 4.0. It is a
coherent, repeating, well-registered source that is **not** one of the 18 eligible kept families.
**The 30 hidden families were NOT checked against it** — by the governance rule they are never
matched — so what F1184 is remains open and is not claimed here in either direction.

---

## §5 Seed-z band of true hp events

`results/p4/p4_0b/hp_seed_z.json` — the median-over-sites z of the 12-s **region-moveout** energy
window at every one of the **2,469** hp events of kept families in the pooled episodes, all scored
on ≥ 3 sites. This is the number that would fix a Puget seed rule.

| percentile | 5 | 10 | 25 | **50** | 75 | 90 | 95 |
|---|---|---|---|---|---|---|---|
| z_med | −0.51 | −0.23 | 0.53 | **2.32** | 6.09 | 12.99 | 20.28 |

* **64.6 %** of true hp events have z_med **below 4** — below Merlin's lower band edge.
* **3.6 %** are above 25 — above the upper edge.
* Only **31.8 %** fall inside the pre-registered **4 ≤ z_med ≤ 25** band *and* carry ≥ 3 sites at
  z ≥ 3.
* Only **50.8 %** have ≥ 3 sites individually at z ≥ 3 at all.

**Caveat.** The frozen 20 Hz cache is band-passed **1–8 Hz**, so this "2–8 Hz" energy is measured in
the cache's 1–8 Hz band. The extra 1–2 Hz octave adds microseism power to both the window and its
background, which will bias z_med low relative to a native 2–8 Hz measurement. The direction of the
bias is known, the size is not measured here.

**Implication.** On this network the 4–25 seed band discards roughly two-thirds of real LFEs before
any detector runs, and the ≥ 3-sites-at-z≥3 clause alone discards half. A Puget seed rule inheriting
these numbers unchanged would be searching the top third of the LFE amplitude distribution. Merlin
pre-registered exactly this measurement as the thing P4-0b should return before P4-2's band is
frozen ("if median-z of true LFE windows on 3+ sites is typically < 4, the seed rule must lower its
floor and the seed budget must rise").

---

## §6 What this shows and what it does not

**It shows**, on 7 PQ-era sites (PGCq, VGZ, SYMB, BUCKS, LHLYB, ALBHB, DAOB) and 203 pooled S18 +
S21 + S22 days, with per-day availability and the frozen measurement chain:

1. Bostock's construction — a **single-event seed → matched-filter growth → two iterated stacks** —
   does **not** assemble known Bostock families to within 0.5× recall and 0.8× purity of those
   families' own truth stacks: 0 of 18 under the pre-registered clauses, 1 of 18 when the survival
   floor is lowered to 8 or 5 independent hits.
2. The failure is a **purity** failure, not a recall failure and not a threshold, seed-count or
   survival-floor failure: recall clauses are met by 5–6 families in the ablation, purity clauses by
   one; purity and hit count are anti-correlated (Spearman −0.68) and purity falls monotonically
   across both iterations for every survivor.
3. A template-free pass over 1,500 correctly-capped seeds **rediscovers one** kept family
   (family 44, 34 matches vs null max 2, 0 aliases) and returns one further coherent repeating
   source that matches no eligible kept family.
4. The corrected seed cap works (max 3/hour, top-10-hours 2 %, 609 hours over 59 days), and the true
   hp-event seed-z distribution sits mostly **below** the band that was pre-registered for Puget.
5. The burst-convergence failure mode Merlin asked to be reported immediately **did not occur**
   (0 of 5 families, hit fractions within ±10 min of the seed 0.011–0.065).

**It does not show** — each of the following is untested here, and this report advocates none of
them:

* that the construction fails on a **multi-year record**. The survival floor of 20 independent hits
  was written for one; 203 days gave a median of 6.
* that it fails with **more than two iterations**. Exactly two were run, by pre-registration, and the
  trajectory over those two is recall-up / purity-down, which does not extrapolate.
* that it fails with a **different growth statistic** — a different per-site clip, weighted rather
  than equal-weight site summing, a coincidence clause stronger than 2 sites at z ≥ 3, a
  cross-correlation threshold set on tremor days rather than on quiet days, or a per-site moveout
  re-pick at iteration instead of one global shift.
* that it fails with a **different template length, pre-roll, or cell**. One combination (20 s from
  S−6 s at pass 0, 12 s from S−2 s after) was run.
* anything about **Puget**. No Puget data was read. The result gates P4-2 by the pre-registered
  rule; it measures nothing about Puget's network or sources.
* what **F1184** is. Hidden families were not matched against it, by governance.

---

## §7 Files and reproduction

**Frozen run.** Pre-registration `results/p4/p4_0b_freeze.json` (created 2026-09-20T04:13:52Z,
sha256 `28fe6d5d02a8cc1c4f39e32b751f7740a8e4864ab552f504bc66012c6c31b229`), hashing the ten scripts
of the run:

| script | sha256 (first 16) |
|---|---|
| `scripts/v_common.py` | `345e91ed89e5ebd3` |
| `scripts/v_engine.py` | `59d5d60d7ca3effe` |
| `scripts/v_seeds.py` | `bb2a1902f389118a` |
| `scripts/v_bank.py` | `0b69e01a5fc3ca95` |
| `scripts/v_run.py` | `0449a64d4a7d2bd2` |
| `scripts/v_iterate.py` | `1d5f74033fd15c5d` |
| `scripts/v_dedupe.py` | `760e558903948679` |
| `scripts/v_score.py` | `317d57691c3d337b` |
| `scripts/v_freeze.py` | `e943d33b6857edb1` |
| `scripts/v_p4_0b_queue.sh` | `57c17b7c861aecab` |

Outputs `results/p4/p4_0b/` — `P4_0B_SCORE.json`, `score_given.csv`, `score_free.csv`,
`seeds_given.{parquet,json}`, `seeds_free.{parquet,json}`, `hp_seed_z.{parquet,json}`,
`bank_*.npz`, `run_*_{thresholds.csv,hits.parquet,indep.parquet,survival.csv}`,
`it{1,2}_given_diag.csv`, `dedupe_given.csv`, `screen_{given,free}.csv`,
`availability_episode_days.csv`, `quiet_days.csv`. Marker `results/p4/P4_0B_DONE`, log
`logs/p4_0b_queue.log`. Re-run: `nohup bash scripts/v_p4_0b_queue.sh > logs/p4_0b_queue.log 2>&1 &`
— every step is guarded by its output file and resumes.

**Ablation.** Pre-registration `results/p4/p4_0b_abl_prereg.json` (created 2026-09-20T05:05:31Z;
records that all ten frozen hashes above still matched at launch), new scripts `scripts/v_abl.py`
and `scripts/v_p4_0b_abl_queue.sh`. Outputs `results/p4/p4_0b_abl/` —
`P4_0B_ABL_SUMMARY.json`, `A1_/A2_score.json`, `A1_/A2_score_by_iteration.csv`,
`A3_{diagnostics.json,template_purity.csv,purity_by_hitcount.csv,frozen_survivor_trajectory.csv}`,
`A4_{reference.json,reference_by_threshold.csv}`. Marker `results/p4/P4_0B_ABL_DONE`, log
`logs/p4_0b_abl_queue.log`. `v_abl.iterate_abl` is a verbatim copy of `v_iterate.iterate` with the
hard-coded 10-hit stacking floor exposed as a parameter; `v_iterate.py` itself was not edited.

**Read-only inputs.** `results/holdouts.json`, `results/p3d/{freeze.json,sse_windows.csv,
arm_b_preregistration.md,stations_p3d_v2.csv,station_epochs.csv}`,
`../lfe_extension/results/{frozen_operating_point.json,a5_station_corrections.csv,
trackE/realness_v1_1.parquet,trackC/eq_catalog.csv}`,
`../lfe_catalog_repair/results/e5r_families_2012.csv`,
`catalogs/pnsn_tremor_cascadia_full.csv`, the 20 Hz caches and `traces/` for the native-rate HF read.

---

# Corrections after Merlin close-out (2026-09-20)

Source: `results/merlin_review_p4_closeout_2026-09-20.md`. **The frozen verdict of §0 stands and is
over-determined**: arm (i) 0/18 and arm (ii) 1 family under the pre-registered clauses, hence
STOP P4-2. Nothing below changes it. What changes are two *mechanism* statements in §3 and §6, which
were wrong and are retracted here. The sections above are left unedited; this is the correction of
record. Two CPU-only measurements were run for it (`results/p4/p4_0b_rescore/`, scripts
`v_rescore.py`, `v_f1184_gate1.py`); both are labelled **post-hoc, not the frozen verdict**.

## C-1 (F1) — RETRACTION: the "purity collapse" mechanism

**§3 and §6 claimed the arm-(i) failure was a purity failure intrinsic to the construction —
"growth admits non-family hits as recall rises", "purity is set at pass 0 and degrades with growth".
That claim is withdrawn.** It was an artefact of the scoring truth, not a property of the detector.

The frozen scoring counted a hit correct only if it one-to-one matched a **p_real ≥ 0.8 detection
labelled to one single family**. Measured on the usable episode days:

* **97.5 %** of the 2,466 kept-family hp events carry `n_other5 > 0` — they are co-detected by at
  least one *other* family's template within 5 s.
* The 18 eligible families are not isolated sources. Nearest non-hidden neighbour, hypocentral:
  52 → **0.12 km** (51), 75 → **0.17 km** (74), 21 → **0.45 km** (20), 36 → 0.92 km (35),
  32 → 1.37 km (36), 56 → 1.50 km (55), 114 ↔ 116 → 1.63 km. Family 44's 0.79 km neighbour (43) and
  family 61's 1.43 km neighbour (60) are **hidden**, so they cannot enter any truth here.

A per-family one-to-one purity against a single-label, p_real ≥ 0.8 truth therefore does not measure
noise admission: a hit on a real event of the neighbouring family 1 km away, or on the same family's
own detection at p_real = 0.5, scores as impure. **The A3 anti-correlation (purity 1.00 at 0–4 hits
→ 0.36 at ≥ 20, Spearman −0.68) is exactly what an incomplete truth produces as a template becomes
more sensitive, and carries no information about noise.** §3's paragraph "Why the floor was not the
cause" is correct in its conclusion (the floor is not what decides the bar) but wrong in its stated
reason; the sentence "the 0/18 is a purity result" is withdrawn.

### Post-hoc re-score with a cluster truth — *not the frozen verdict*

`scripts/v_rescore.py` → `results/p4/p4_0b_rescore/{rescore_cluster_truth.csv,rescore_summary.json}`.
Truth for family X = frozen-catalogue detections (`realness_v1_1`) on the usable episode days of
**every non-hidden family within 3.0 km hypocentral of X**, deduplicated at 2 s; primary at
p_real ≥ 0.2, secondary at any p_real. Templates, hits, reference and the 0.5× / 0.8× relative
clauses are unchanged. Hidden families are excluded from the truth and from all reporting.

| | truth p_real ≥ 0.2 | truth any p_real |
|---|---|---|
| **frozen** it2 templates | **1 / 18** (family 44) | **2 / 18** (families 44, 75) |
| **A2 ablation** it2 templates | **5 / 18** (21, 44, 56, 75, 110) | **6 / 18** (21, 44, 56, 75, 110, 116) |

Merlin's acceptance was "≤ 6/18 confirms closure; ≥ 12/18 is a failure signature". **The maximum is
6/18. The failure signature did not fire; closure is confirmed.** His predicted frozen count (2/18)
and his measured G44_0 relative purity (0.85) are reproduced exactly: under any-p cluster truth
G44_0 has purity 0.667 against REF44's 0.782 → relative **0.852**, and it assembles. Family 75
assembles at 0.806 relative purity. The five frozen survivors under any-p cluster truth:

| family | cluster (non-hidden ≤ 3 km) | truth events | tid | hits | recall / need | purity / need | assembled |
|---|---|---|---|---|---|---|---|
| 32 | 3 families | 3,530 | G32_1 | 35 | 0.007 / 0.002 | 0.686 / 0.753 | no |
| 44 | 1 (43 is hidden) | 721 | G44_0 | 114 | 0.105 / 0.055 | **0.667 / 0.626** | **yes** |
| 61 | 3 | 2,689 | G61_0 | 92 | 0.017 / 0.011 | 0.511 / 0.715 | no |
| 75 | 6 | 2,491 | G75_1 | 147 | 0.038 / 0.018 | **0.639 / 0.634** | **yes** |
| 110 | 1 | 501 | G110_1 | 92 | 0.108 / 0.082 | 0.587 / 0.669 | no |

The corrected reading of arm (i): **the construction assembles 1–2 of 18 families under the frozen
survival floor and 5–6 of 18 when the floor is lowered, on a 203-day record.** That is a
power-and-seed result, not an intrinsic incapacity, and it is still far below the 12/18 bar.

## C-2 (F4) — the 8.0 growth floor was never exercised

The FA-only rule never bound. `run_given0_quiethits.parquet` holds **13 hits at `net_mad ≥ 8.0` over
54 templates × 20 quiet days = 1,080 template-days → 0.0120 false hits per template-day, against an
allowance of 0.100** — eight times below. `v_engine.run_day` applies `max(thr, floor)`, so nothing
below 8.0 was ever recorded and the rule could not return a lower T. The pass-0 starvation that
removed 13 of 18 families before iteration (median 6 independent hits) is therefore **partly a floor
convention, not a measured noise limit**. §1's observation that "the growth threshold is set by the
recording floor, not by measured noise" was correct but was not carried into §3/§6 as a limitation:
the ablation moved only the *survival* floor, never this one, so "the floor was not the cause" was
tested for one floor and asserted for the other.

## C-3 (F3) — RETRACTION: arm (ii) failed at the SEED stage, and §5's implication is withdrawn

§5 concluded that "the 4–25 seed band discards roughly two-thirds of real LFEs" and implied the
remedy was to lower the floor and raise the budget. **That implication is withdrawn.** Measured
against the frozen catalogue on the usable episode days (`seeds_free.parquet`, ±2 s, one-to-one not
required):

| seeds | within ±2 s of a **kept hp event** | within ±2 s of **any frozen detection** |
|---|---|---|
| all 1,500 | **2.20 %** | **12.60 %** |
| same-day random chance (30 reps) | 0.21 % ± 0.13 | 2.82 % ± 0.38 |
| enrichment | **10.7×** | **4.5×** |

| z_med band | n | any detection | kept hp |
|---|---|---|---|
| 4–6 | 234 | 9.4 % | 1.7 % |
| 6–8 | 252 | 13.5 % | 2.8 % |
| 8–12 | 356 | 12.9 % | 2.5 % |
| 12–25 | 658 | 13.2 % | 2.0 % |

The seed statistic **is** enriched for catalogued activity — 4.5–10.7× chance, so it is not picking
noise at random — but in absolute terms **~88 % of seeds are not catalogued at all, only 2.2 % are
kept LFEs, and the yield is flat in z_med**. A flat yield means no choice of band, floor or budget
improves it: moving the band cannot raise a fraction that does not vary with the band. Merlin
reports 9.7 % / ~7 % chance for the any-detection figure against my 12.6 % / 2.8 %; the difference is
in the denominator convention (his chance level is higher), and both readings agree on the two facts
that matter — the absolute LFE fraction is a few per cent, and it is flat across z.

**Corrected statement:** there is no energy-based seed rule for LFEs on this network. Arm (ii)'s
single rediscovery is consistent with arm (i)'s per-true-seed yield (~9 % pass-0 survival) times a
~2–5 % LFE fraction among seeds. For Puget, which has **only** template-free seeds, this is the
decisive negative, and it is independent of P4-0a's geometry result.

§5's percentile table itself stands as measured, with F9's caveat sharpened: `hp_seed_z` is computed
at the **region-centre** moveout and in the **1–8 Hz** cache, so "64.6 % of true LFEs below z_med 4"
is the operational property of the free-arm seed rule, not the LFE amplitude distribution.

## C-4 (F5) — the truth-stack reference is circular

`v_bank.build_ref` stacks each family's reference template from the same hp events that
`v_score._hp_by_family` then scores against, on the same days. The reference is trained on the test
set. Deviation 3 called this "conservative"; it is (it inflates `recall_ref`/`purity_ref`, making the
relative bar harder), but it is **not fair**, and "conservative" understated the problem. The
cluster re-score in C-1 keeps the same reference and still moves the count to 1–2/18 frozen and
5–6/18 under A2, so the circularity is not what decides the outcome — but no future run should reuse
this construction. An independent (transferred) reference is the correct comparator.

## C-5 (F6) — the N-moveout null is a window-exclusion null

`v_bank.build_nmoveout` re-cuts the seed with a wrong cell's per-station timing. A pass-0 template
cut from a **real event** carries the true inter-station timing inside its samples regardless of the
cell assumed for the cut, so the network sum still coheres at the true origin — unless the shift
(median 8 s, range 4–15 s, `bank_given_nm_cells.csv`) pushes the arrival outside the 20-s window at
enough stations. It is therefore a *window-exclusion* null, not a moveout scramble. It worked
operationally (17 of the 19 free templates with ≥ 20 independent hits failed it; all five given
survivors passed) and should be documented and cited as such, not as a test of moveout coherence.

## C-6 — credit where it is due: F187 → family 44

The one positive datum of the whole exercise, and it should be stated plainly. `F187`'s seed is a
12-s tremor-hour window **184 s from the nearest kept hp event** (and 1.0 s from a low-`p_real`
frozen detection). Growth from that non-catalogued window, two iterations, and dedupe produced a
template whose 95 independent hits match **34** of family 44's 95 hp events one-to-one at ±2 s —
35.8 %, against a 20-realisation within-day null of mean 0.60 / max 2, with **0 aliases** (next-best
family 36 at 9 matches). Its hits carry a median of **21 PNSN epicentres within 60 km ±30 min**.
**A template-free seed converged on a real, known family.** That the construction can do this at all
is the reason the north negative is read as power-and-seed-limited rather than as incapacity.

(One caution on how tremor co-timing is quoted: F187's raw co-timing fraction is 1.000, but its
same-day random base rate is 0.993 — all seven of its days are tremor days, so the fraction is
uninformative. The discriminating statistic is the epicentre count: median 21 for F187 versus
median 0 for F1184.)

## C-7 — F1184 Gate 1 at the never-summed sites: **FAIL** (descriptive)

`scripts/v_f1184_gate1.py` → `results/p4/p4_0b_rescore/F1184_gate1.{json,csv}`. Alignment-free mean
stack of F1184's 24 independent hits at the bank's region-centre moveout (cumulative iteration shift
−4.45 s), carried to the seven sites that were never in any network sum; `s_common.s_excess_db`
maximised over a ±1.5 s constant offset; 200 same-day random-time nulls charged the identical search.

| site | n stacked | S excess (dB) | null median | null p99 | above p99 |
|---|---|---|---|---|---|
| CLRS | 24 | **−1.82** | +0.53 | +2.71 | no |
| GOBB | 19 | **−2.42** | +0.45 | +2.59 | no |
| LOPEZ | 24 | **−2.60** | +0.66 | +2.68 | no |
| TURTL | 24 | **−1.02** | +0.72 | +3.22 | no |
| LRIV | 24 | **−1.58** | +0.37 | +2.29 | no |
| NLLB | 24 | **−0.67** | +0.47 | +2.08 | no |
| B926 | 18 | **−1.17** | +0.38 | +1.94 | no |

**0 of 7 above null p99 → Gate 1 FAILS** (bar: ≥ 2). Every observed excess is *negative* while every
null median is positive. Tremor co-timing **0.500** against a same-day base rate of **0.522**
(p05 0.458, p95 0.625) — indistinguishable — with a median of **0** PNSN epicentres within 60 km
±30 min. F1184 therefore shows no arrival at any site outside the detection network, is not
tremor-locked, and matched no eligible kept family (best 1 match, null max 1). The +7.7 to +19.1 dB
excess it showed at the *detection* sites in §4 is what Merlin's F3 predicted: selection bias, since
the hits were chosen by correlation with the template at exactly those sites. §4's description of
F1184 as "a coherent, repeating, well-registered source" should be read as **a repeating waveform
coherent across the detection network only**; on this evidence it is not an LFE family. As stated in
§4 and preserved here, the 30 hidden families were never matched against it.

## C-8 (F10) — VERIFIED as a mechanism, immaterial in this run

Merlin's speculative F10: `v_engine.DayData.station_z` computes the median and lower-half MAD of the
cc series over **all** samples, gaps included (where cc ≈ 0), which on a gappy station-day could
deflate the MAD and inflate z. Measured (`results/p4/p4_0b_rescore/F10_gap_mad_check.json`, CPU,
one template, exact `station_z` arithmetic reproduced in numpy):

| station-day | zero-sample fraction | MAD all samples | MAD valid only | z inflation |
|---|---|---|---|---|
| DAOB 2024-08-12 (the gappiest in the run) | 4.54 % | 0.0667 | 0.0707 | **+5.9 %** |
| B009 2024-08-15 | 0.89 % | 0.0486 | 0.0492 | **+1.2 %** |
| DAOB 2024-08-12 with a synthetic 30 % gap | 30 % | 0.0321 | 0.0722 | **+125 %** |
| B009 2024-08-15 with a synthetic 30 % gap | 30 % | 0.0248 | 0.0481 | **+94 %** |

**The mechanism is real** — a 30 % gap deflates the lowhalf MAD by a factor ≈ 2 and inflates z by
~100 %. **It is immaterial here**: across 631 station-days scanned over the first 70 episode days,
only **0.16 %** have more than 1 % zero samples and **none** exceed 10 %; the worst day in the entire
run costs +5.9 % in z at one station. Gap samples also cannot themselves become hits — the per-sample
`valid` mask zeroes their contribution to the site mean and to the ≥ 2-sites-at-z ≥ 3 coincidence
count — so the only leak is the scale. **Answer: confirmed as a mechanism, bounded at ≈ 6 % in z on
the single worst station-day, no effect on any result reported here.** It should nonetheless be
fixed (restrict the median/MAD to valid samples) before the code is reused on gappier data.

## C-9 — corrected summary of what P4-0b establishes

Replacing the first two bullets of §6:

1. On 7 PQ-era sites and 203 pooled episode days, the construction assembles **1–2 of 18** known
   families under the pre-registered survival floor and **5–6 of 18** with the floor at 5–8
   independent hits, when scored against a 3-km cluster truth that reflects how the frozen catalogue
   actually labels these events. Against a bar of 12/18 this is a clear negative.
2. The binding limits are **record length and seeds**, not an intrinsic purity ceiling: 13 of 18
   families never produced a template with ≥ 20 independent hits in 203 days, at a growth threshold
   the false-alarm rule never exercised (C-2); and the template-free seed rule selects catalogued
   LFEs at 2–5 % regardless of z-band (C-3).

Everything else in §6 — including the full list of what is *not* shown — stands as written.

---

## Standing conclusion (Merlin, close-out review Q4, verbatim)

> "We tested whether new LFE families can be found without a pre-existing catalogue — first in the
> north, where Bostock's families are truth, then as a plan for Puget. They cannot, with the tools we
> have. Pairwise linkage (P3D) found no clusters. The real seed-and-grow construction, run on the
> north with 7 good sites over 203 episode days, assembled at most 2–3 of 18 known families from
> single true events (13 never produced a template that survived the first pass), and a
> template-free seed rule picked windows that were catalogued LFEs only 2–7 % of the time,
> rediscovering one family. Puget has 6–10 usable sites per sub-region at 50–65 km from the tremor
> field, covering at most a fifth of the active cells with three sites, about 67 build-period tremor
> days, and no seeds. The Puget discovery line is closed. The north negative is a power-and-seed
> result, not proof that the construction is intrinsically incapable — the earlier 'purity collapse'
> reading was a scoring artefact and is withdrawn. The only route to LFE families that has worked in
> this project remains growth from known events or stacks (the transferred-template catalogue)."
