# Reproducing Bostock's LFE families from the public archive — north (southern Vancouver Island)

Agent K, 2026-09-20. Full working ledger: `results/REPRO_BOSTOCK.md`. Method = Bostock et al.
(2012) / Brown et al. (2008) implemented literally from the published descriptions (6-s windows,
1–8 Hz, per-component cc summed over all channels, threshold = median + 8 × MAD, 4-s declustering,
components linked into groups, stacked with Bungum–Husebye realignment, 12-s station-tapered
templates, network matched filter at 8 × MAD, iterate). No thresholds calibrated on quiet days, no
survival gate, no scramble ratio, no coincidence clauses — none of those appear in the papers.

## The headline

| | result |
|---|---|
| **KEPT families (development truth), 2005 POLARIS bed, 15 stations** | **58 of 96** rediscovered at Bostock's ≥30 % one-to-one ±2 s criterion (63 of 100 scoring on raw catalogue times) |
| **HIDDEN families (sealed until the single blind run)** | **22 of 30** — 6 of the 10 "isolated" stratum |
| Ceiling of this bed (templates built from the catalogue itself, identical filter) | 89 of 96 kept, **29 of 30 hidden** |
| So, as a fraction of what this bed can possibly yield | kept **58/89 = 65 %**, hidden **22/29 = 76 %** |
| Detection recall per family | 2–4 % of catalogued events (ceiling on the same bed: 31 %) |
| **Permanent-station bed (2012, 15 permanent stations, same code)** | **3 of 77** — ceiling on that bed 33 of 77 |

So the array-era discovery **is** reproducible from the public archive: the engine, seeded only by
its own autocorrelation, recovers roughly two thirds of Bostock's families on a single 24-day ETS
episode, and the blind number (73 %) matches the development number (60 %).

## What had gone wrong before

Five defects, all ours, none of them a threshold:

1. **The 2005 catalogue times are 100 % whole-second** (2012: 2.6 %). Every "truth" stack the study
   had built was incoherent, which made a sound bed look unusable. Fixed by matched-filtering each
   event against an iterated stack of the family's loudest events: stack coherence at disjoint
   verification sites 1.26 → 3.76 (shifted-catalogue null 1.26).
2. **The same jitter one stage later:** the autocorrelation's 0.5-s row grid leaves ±0.25 s on every
   node time, whose plain-stack coherent gain |sinc(f·0.5 s)| is zero at 2, 4, 6 and 8 Hz. This is
   exactly what Bostock 2015 §2's Bungum–Husebye realignment exists to remove; it had never been
   implemented.
3. **Invented statistics.** The frozen engine replaced the published plain channel sum with a
   per-station lower-half-MAD z, clipped at 8, site-averaged and row-renormalised, then added
   `n_site_z3 ≥ 2`, `n_sites ≥ 4`, a quiet-day false-alarm threshold and a ≥20-hit/≥3-day/scramble
   survival gate. None is in the papers. Removing them turned a statistic that had produced 0
   usable seeds into one whose strongest 100 detections are 71 % true.
4. **Gap handling.** Zero-energy windows normalise to unit-norm noise and correlate at ~1; seven of
   576 hours produced 1.16 M of 1.32 M "detections", and `1/sqrt(0)` in the matched filter made
   every gap edge a detection.
5. **The merge.** Plourde's >10 %-shared-detection test is specified at *sample* precision; at the
   ±2 s we had used it chained neighbouring families and destroyed 20 real ones (50 → 9). Measured
   against truth, *neither* published uniqueness test separates families here (same-family stack cc
   0.072 vs different-family 0.060), because Bostock's families sit 0.1–1.6 km apart and co-detect.
   No merge is applied; duplicates must be resolved by location.

## Was "template-free discovery is impossible on the permanent network" true?

Largely yes, and now with numbers. Same code, same thresholds, 2012 permanent bed (15 permanent
stations, 52 days): the autocorrelation runs at **1.7× chance** (0.05 % true pairs vs 0.03 %)
against **11×** on the array bed (12.85 % vs 1.17 %); only 4 groups reach purity ≥0.8 against 205;
**3 of 77** families are rediscovered. The growth machinery is not the problem — fed catalogue-built
templates the same bed yields 33 of 77 families at 33 % recall. It is the pair stage that dies
without the array.

## Three beds, side by side

| bed | stations | days | true autocorrelation pairs | chance | **families discovered** | ceiling |
|---|---|---|---|---|---|---|
| 2005 POLARIS array | 15 | 24 | **12.85 %** | 1.17 % | **58 of 96** | 89 of 96 |
| 2012 permanent | 15 | 52 | **0.05 %** | 0.03 % | **3 of 77** | 33 of 77 |
| 2010 permanent | 8 | 34 | (re-running after two cache defects) | | | |

## Candidate new families

Of 1,457 grown templates: 177 match a kept family (58 distinct), 64 match a hidden family (18
distinct), 1,201 do not locate, and **18 pass the candidate rule** (located on ≥5 stations,
shuffled-pick control passed, within 8 km of Slab2, repeating on ≥3 days). Ten of the eighteen rail
at 18 km, the shallowest node of the search grid — their depth is unconstrained, and four of those
also sit 20–80 km from any catalogued family, where a persistent local noise source is the likelier
explanation. **Seven survive a conservative reading** (interior depth 26–38 km, clear control
margin, 3–5 days): grp 534, 1215, 912, 889, 247, 745, 822. Six of the seven lie **1.3–8.1 km from a
catalogued family** and are most plausibly unlisted neighbours inside Bostock's own field; only
**grp 1215 (47.90 N, −123.66, 32 km depth, 24.6 km from the nearest kept family)** is spatially
distinct. Table: `results/repro/candidate_new_families.csv`.

**Nothing here is claimed as a discovery.** The confirmation test is a re-fire: POLARIS-BC ran
2004-04 → 2005-10 and the 2004-07 ETS episode (59,659 catalogued events over 18 busy days, ~900
station-days, downloadable) has never been touched. A candidate becomes a family if its unchanged
2005 template re-fires there — ≥20 independent detections on ≥3 days at the same 8 × MAD rule, with
a location reproducing to within the 2005 uncertainty. That has not been run.

## In plain terms

Bostock's 130 families **are** recoverable from the public archive with his own published recipe.
Starting from nothing but the raw 2005 waveforms, the chain finds its own seeds by correlating
tremor against itself, and recovers **58 of the 96 development families and 22 of the 30 families
that were sealed away before any of this was written** — 65 % and 76 % of what the same 15 stations
can yield even when handed the answer. The earlier conclusion that this was impossible came from
five defects of ours, not from the data: whole-second catalogue times that made every truth stack
incoherent, the same sub-sample jitter again at the stacking stage, a home-made detection statistic
in place of the published one, unhandled data gaps, and a duplicate-merging rule applied 80× too
loosely.

What the **permanent network alone** can and cannot do is now a number rather than an opinion. It
can *run* families: given a correct template it recovers 33 of 77 families on the 2012 bed at a
third of their events. It cannot *find* them: its autocorrelation sits at 1.7× chance instead of
11×, only 4 clean groups ever form, and 3 families are discovered. Template-free discovery needs
the array — which is exactly why every published Cascadia catalogue was built on one.

## Addendum (orchestrator, 2026-09-21T00:13:44Z): 2010 permanent bed
2010 PERMANENT bed (8 stations with coordinates, 34 days; K_REPRO6, logs/k_repro6_queue.log): 1,181 groups, 19 at purity ≥0.5 covering 15 kept families, 7 at ≥0.8; autocorrelation discovery **5 of 60** kept families in window; catalogue-built ceiling **50 of 60**. Same pattern as 2012 (3/77 vs 33/77): the permanent network runs known families but does not find them.
