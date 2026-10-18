# Reproducing Bostock et al. (2012) on the 2005 POLARIS bed — dated ledger

Agent K. Development truth = the **100 KEPT families**. The 30 hidden families
(`../lfe_picker/results/holdouts.json`) are untouched until the single blind run.
No Merlin, no freezes, no gates. Every entry = what was measured, what was found wrong, what
changed.

---

## 2026-09-20 — §0. The published recipe, beside what our engine did

Sources actually read (PDF or publisher text; Wiley/Elsevier full text is 403 from here, so the
Brown-2008 and Bostock-2012 parameters marked (†) come from the abstracts plus the verbatim
restatements in the open-access secondary literature listed below).

| parameter | published | our frozen engine (`w_autocorr.py`, `w_common.py`) |
|---|---|---|
| band | **1–8 Hz** (Bostock 2015 §2, templates "filtered between 1–8 Hz") | 1–8 Hz ✔ |
| sample rate | not stated | 20 Hz |
| autocorrelation window | **6 s** (†Royer & Bostock 2014 restating Bostock 2012: "autocorrelations of 6-second long windows") | **10 s** ✘ |
| channels | **each component** of each station, summed (†); Bostock 2015 uses P on Z, S on N+E | 3-C mean per station, then site mean ✘ |
| stations | **7** stations, southern Vancouver Island (†Bostock 2012 abstract) | 12 "sites" after a 1-km merge, 11 clean |
| network function | **plain sum/mean of the per-channel cc** (†; Shelly 2007 convention) | per-station lower-half-MAD **z**, clipped at 8, site-averaged, then the row re-normalised by its own median / lower-half MAD ✘ |
| detection threshold | **8 × MAD** of that summed function (Shelly 2007; restated verbatim in Lin et al. 2024 Seismica §1) | FA-only threshold derived from 20 **quiet days**, engineering floor 5.0; 1c′ landed on **T = 8.0–12.5** in the z-statistic ✘ |
| extra coincidence clauses | **none stated** | `n_site_z3 ≥ 2` and `n_sites ≥ 4` ✘ (invented) |
| declustering | strongest detection in a **4 s** interval (Bostock 2015 §2 culling rule) | 10 s |
| data span per autocorrelation | Brown 2008 tested on **1 hour**; Bostock 2012 on 2004 + 2005 tremor | 1 hour ✔ |
| grouping | "**waveform cluster analysis**" of the 5,775 detections → **140 families** (†Bostock 2012 abstract) | ≥5 members in one hour **and** ≥2 ten-minute bins (1c), or ≥2-node components (1c′) ✘ |
| template build | iterative shifting against the running composite stack until convergence (Bostock 2015 §2, after Bungum & Husebye 1971); Plourde 2015 (same group): initial 60 s template, later 12 s with a station-dependent taper; phase-weighted stacks early, linear stacks late | `align_stack` ±1.5 s, 2 passes; 20 s pass-0 / 12 s iterated template cut at an energy-onset argmax searched over −4…+16 s |
| matched-filter threshold | 8 × MAD, same rule | FA-calibrated `FA_PER_TPL_DAY = 0.1` per template-day ✘ |
| iterations | "until convergence"; Plourde: several, phase-weighted then linear | 2 |
| survival / QC | none stated; uniqueness tested afterwards by shared-detection fraction >10 % or inter-template cc >0.3 (Plourde 2015) | **≥20 independent hits, ≥3 distinct days, scramble ratio ≥3, `frac_near_seed` ≤0.20** ✘ (all invented) |
| seeds | autodetection **plus visual identification** (Plourde 2015 §2) | purely automatic |

Sources: Bostock, Thomas, Savard, Chuang & Rubin (2015) *JGR* 120, 6329 §2 (PDF,
`pages.uoregon.edu/amthomas/pubs/bostock_jgr_2015.pdf`); Plourde, Bostock, Audet & Thomas (2015)
*GRL* 42, 4849 §2 (same host); Lin, Thomas, Bachelot, Toomey, Searcy & Melgar (2024) *Seismica*
2.4 §1 (open access) for the "eight times the median absolute deviation" rule; Bostock, Royer,
Hearn & Peacock (2012) *G-cubed* 13, 2012GC004391 abstract (7 stations, 5,775 LFEs, 140 families,
waveform cluster analysis); Royer & Bostock (2014) *EPSL* 402, 247 as quoted for the 6-s window.
**Assumptions recorded where a parameter is not stated:** sample rate = our 20 Hz cache; MAD taken
over the whole block's correlation values (not per row); minimum |lag| = one window length (6 s);
row step 0.5 s.

---

## 2026-09-20 — §1. The 2005 catalogue times are whole-second, and that alone made every truth
stack incoherent

`bostock2015_sloweq_FULL_2003-2013.csv`: **100.0 %** of the 57,504 2005 rows have integer `sec`
(2012: 2.6 %).

`scripts/k_s1_loud.py` — for the 20 most productive kept families in the 24-day bed
(2005-08-29…09-21, 16,900 events): loudness = network-median 3-s S-energy excess over the local
±60 s tremor background at the 4 nearest detection sites (top decile = LOUD); then each event is
matched-filtered against an iterated stack of the family's loudest events (Bostock 2015 §2's own
Bungum–Husebye alignment), with the **alignment stations and the verification stations disjoint**
so the coherence check is not the selection statistic.

| stack coherence (1 = incoherent) | catalogue times | corrected times |
|---|---|---|
| alignment sites (selection-biased) | 1.19 | 19.3 |
| **verification sites (disjoint, honest)** | **1.26** | **3.76** |
| verification sites, +300 s shifted-catalogue null | 1.03 | **1.26** |

Corrections: |dt| median 0.55 s, sd 0.76 s. Network cc of a loud event against its family stack
0.21 (median), non-loud 0.13.

**Consequence.** Merlin's F3 was right that the truth stacks were incoherent, and the repair is the
one above. With corrected times the *predicted* S is fine after all: the empirical envelope-peak
offset from `fam_tt` (`scripts/k_s2a_stacks.py`) is **0.13–1.15 s at every surface site**
(TWKB 0.55, MGCB 0.40, GLBC 0.45, PGC 0.13; peak heights 8–62 MAD), not the 2.8 s median that
STEP1CP_RESULTS §a reported from incoherent stacks. Only the four boreholes (55–60 km) stay
incoherent (coh 1.10–1.18) — they do not see these events.

---

## 2026-09-20 — §2. What the pair statistic can and cannot see (measured, not assumed)

`scripts/k_s2c_sweep.py`, loud pairs of the 6 most productive kept families in their busiest
hours, at corrected times, cut at the **empirical** S, 20 Hz 1–8 Hz, per-station 3-C cc averaged
over 9–11 stations, against random window pairs in the same hour:

| window | matched-filter cc (template vs event) | pair cc (event vs event) |
|---|---|---|
| 1.5 s | 0.224 | 0.051 |
| 6 s | 0.184 | 0.027 |
| 10 s | 0.165 | 0.021 |

pair cc ≈ (template cc)² exactly as expected when both events are noise-dominated: the pair
statistic is the *square* of the matched filter's. That is the physics, and it is why the loud
tail — not the median event — is what autocorrelation sees.

---

## 2026-09-20 — §3. The literal recipe, run

`scripts/k_auto.py` implements the table in §0 with every invented clause removed: 6-s windows,
0.5-s row step, full-sample lags, per-component normalised cc, **plain mean over all 33 channels**,
threshold = median + 8 × MAD of the hour's own correlation values, |lag| ≥ 6 s, 4-s declustering.

* 8 × MAD is empirically **5.15 σ** of the network function (so the rule is not a high bar, it is
  a *measured* one).
* Two hours of 2005-09-13 (10h, 11h), 11 clean sites: **5,957 detections in 8 s of GPU time.**

Scored (`scripts/k_score_pairs.py`) against the corrected-time catalogue of only the 20 families
that have corrections so far — a detection is a true pair of family k if its lag reproduces a
catalogued inter-event interval of k to ±0.6 s and its first window sits within ±8 s of that
event:

| subset | matched | families | chance (+300 s) |
|---|---|---|---|
| all 5,957 | **35.0 %** | 6 | 5.4 % |
| top 2,000 by cc | 43.0 % | 6 | 2.1 % |
| top 500 | 50.0 % | 5 | 0.0 % |
| top 100 | **78.0 %** | 2 | 0.0 % |

**This is the stage the study had recorded as a null.** The frozen engine's per-station z /
site-merge / clip / `n_site_z3 ≥ 2` / `n_sites ≥ 4` / quiet-day FA threshold, applied to 10-s
windows, turned a statistic that is 78 % pure in its top 100 into 0 usable seeds.

Next: corrected times for all 100 kept families, the same autocorrelation over all 24 ETS days,
then grouping → stacking → iterated matched filter, scored by Bostock's own criterion.

---

## 2026-09-20 — §4. Full episode, pair stage (24 ETS days, 11 clean sites)

576 station-hours, 4 s of GPU each. Two data-quality defects found and fixed in `k_auto.py`
(both are gap handling, not thresholds):

1. **Degenerate windows.** A window inside a gap has ~zero energy; normalised, it is a unit-norm
   noise vector that correlates ~1 with another such window. 7 of 576 hours had their block sd
   blown to 27–473 and produced **1,159,000 of 1,322,000 "detections"**. Degenerate windows are
   now excluded on both axes (`emed`/`good` mask).
2. **`cc` is not comparable across blocks.** Each hour has its own noise level, so ranking seeds by
   raw `cc` puts the noisiest hours on top — the global top 100 by `cc` were 0 % true. Every
   detection now carries `nmad = (cc − median)/MAD` of its own block.

| clean detections (480 hours, 143,441 pairs) | matched | families | chance (+300 s) |
|---|---|---|---|
| all | 8.41 % | **88 of 96** | 0.92 % |
| top 50,000 by nmad | 13.9 % | 77 | 1.00 % |
| top 10,000 by nmad | 24.3 % | 54 | 1.66 % |
| top 2,000 by nmad | 34.1 % | 29 | 3.45 % |
| top 500 by nmad | 45.8 % | 10 | 4.20 % |
| **top 100 by nmad** | **71.0 %** | 9 | 0.00 % |

8 × MAD is **5.13–5.24 σ** of the network function on every hour measured.
7 stations (the 7 nearest, Bostock's count) vs 11 on the same two days: 11 sites win —
`results/repro/auto7` vs `auto11_2d` (logs `k_score_auto7.log`, `k_score_auto11_2d.log`).

## 2026-09-20 — §5. Grouping (`k_group.py`)

Their "waveform cluster analysis" read literally: windows (merged at the 4-s declustering interval)
are nodes, detections are edges, a candidate family is a connected component. No minimum size is
stated; we keep components with ≥3 nodes.

| budget (strongest detections by nmad) | nodes | components | ≥3 members | purity ≥0.5 | families covered |
|---|---|---|---|---|---|
| 5,000 | 289 | — | 37 | 16 | 12 |
| 20,000 | 2,845 | 1,069 | 235 | 70 | 41 |
| 60,000 | 8,762 | 2,860 | 701 | 202 | 57 |
| **143,441 (all above 8×MAD)** | 19,245 | 5,771 | **1,397** | **337** | **75** |

No giant-component collapse (largest 116 nodes). Note that groups matching one of the 30 **hidden**
families are scored 0 here by construction — they are never inspected.

## 2026-09-20 — §6. Stack → template → matched filter → iterate (`k_grow.py`)

Two more defects of the same family, found by running it:

3. **Matched-filter denominator in gaps** — `1/sqrt(0)` gave cc ≈ 10¹³ and every "detection" on a
   gappy day was a gap edge (798 → 226 detections for the same 19 templates once fixed).
4. **The 0.5-s row grid leaves ±0.25 s of jitter on every node time**, whose coherent gain for a
   plain stack is |sinc(f·0.5 s)| — **zero at 2, 4, 6 and 8 Hz**. This is the same defect as the
   whole-second catalogue times in §1, one stage later, and it is exactly what Bostock 2015 §2's
   Bungum–Husebye step exists to remove. Implemented: round 0 plain stack, then each member
   re-shifted on the **network** sum of its per-station cc against the current stack, 2 rounds.

Chain as run: long 44-s stack from node − 12 s at every station → 12-s template cut 2 s before the
**stack envelope peak at that station** (Plourde's station-dependent taper) → network mean of the
per-station normalised 3-C cc over each ETS day → median + 8 × MAD → strongest per 4 s → restack
the top 200 detections → repeat.

**Positive control on clean seeds.** 25 groups at purity ≥0.8 with ≥8 nodes (19 distinct kept
families), grown over all 24 ETS days, scored by Bostock's criterion (≥30 % of the grown family's
detections match one family's catalogued events one-to-one at ±2 s):

| iteration | grown families | meeting the ≥30 % bar | **distinct kept families rediscovered** |
|---|---|---|---|
| 0 | 25 | 14 | **14 of 19 seeded** |
| 1 | 24 | 16 | 14 |
| 2 | 24 | 14 | 13 |

So the growth stage is not the blocker: given a clean seed it reproduces the family. Running now
over **all 1,397 groups** (`scripts/k_repro2_queue.sh`, markers `K_REPRO2_DONE|FAILED`) with
dedupe by Plourde's >10 %-shared-detections rule and full per-stage loss accounting.

---

## 2026-09-20 — §7. Full run, and where the families are lost

`k_repro2_queue.sh` grew **all 1,397 groups** through 4 passes (1,255 templates at pass 0; detections
per template rose 13 → 201 as the stacks enriched). Scored with Bostock's criterion (≥30 % of a
grown family's detections match one kept family's catalogued events one-to-one at ±2 s):

| | templates | ≥30 % bar | **distinct kept families** |
|---|---|---|---|
| iteration 0, no merge | 1,250 | 137 | **50 of 96** |
| iteration 1 / 2 / 3, no merge | 1,207 / 1,139 / 1,010 | 136 / 123 / 114 | 47 / 47 / 40 |
| iteration 0, merged at ±2 s shared times | 204 | 37 | 28 |
| iteration 0, merged at **sample precision** | 255 | 41 | 30 |
| iteration 3, merged at sample precision | 238 | 19 | 18 |

**The merge was the dominant loss (50 → 30), and it is not recoverable by tightening it.** Plourde
et al. (2015) share detection times "to the nearest sample (Δt = 0.025 s)"; we had used ±2 s, an
80× looser test — fixed. But the decisive measurement is that **neither published uniqueness test
separates our families at all**:

| labelled edge (both ends pass the ≥30 % bar) | n | template stack cc | shared-detection fraction |
|---|---|---|---|
| SAME kept family | 23 | med 0.072 (p10 0.059, p90 0.102) | med 0.143 |
| DIFFERENT kept families | 12 | med 0.060 (p10 0.049, p90 0.085) | med 0.167 |

At the cross-family p95 (cc 0.088) only 26 % of same-family edges survive; the share test is worse
(13 %). Plourde applied cc > 0.3 to converged high-SNR templates and to families far enough apart
to be separable; ours top out at cc 0.242, and Bostock's southern-Vancouver-Island families sit
0.1–1.6 km apart and co-detect. **Decision: no merge is applied.** The count reported is the number
of DISTINCT kept families represented among the grown templates; the duplication (1,250 templates
for ~50 + hidden families) is a real limitation of the deliverable catalogue, to be resolved by
location, not by waveform similarity. `scripts/k_dedupe2.py` implements both tests and the
calibration; a third bug is recorded in it (an `np.roll` by up to 640 samples on a 121-sample lag
axis made every template cc NaN).

**Per-stage loss (96 kept families with corrected times):**

| stage | families surviving |
|---|---|
| seen by the pair stage (≥1 true autocorrelation pair) | **88** |
| ≥1 group at purity ≥0.5 | 75 |
| ≥1 group at purity ≥0.8 | 58 |
| grown and passing Bostock's ≥30 % bar | **50** |

Of the 75 seeded families, 47 were rediscovered. Seed size is the discriminator:

| seed group size (nodes) | rediscovered / seeded |
|---|---|
| ≤4 | 26 / 47 |
| 5–6 | 5 / 10 |
| **7–10** | **11 / 12** |
| >10 | 5 / 6 |

Failed families are also smaller in the catalogue (median 297 events vs 516 for rediscovered ones).
Iterating does NOT help (50 → 47 → 47 → 40): as a template enriches it also absorbs its
neighbours' events, so the one-to-one purity falls. Pass 0 is the best pass.

**Station count is a real lever** (2005-09-13 + 09-16, same code, same threshold):

| stations | detections | true pairs | families | top-2,000 purity |
|---|---|---|---|---|
| 11 clean | 24,344 | 23.8 % | 42 | 76 % |
| **15 (adds JRBC, LZB, VGZ, SNB)** | 25,362 | **30.5 %** | **45** | 92 % |

Running now (`scripts/k_repro3_queue.sh`, markers `K_REPRO3_DONE|FAILED`): the whole episode at 15
stations, regroup, grow every group, plus (i) a looser 6×MAD FIRST pass with 8×MAD afterwards, and
(ii) the **catalogue-built ceiling** — templates made from the corrected-time catalogued events of
every kept family, through the identical matched filter — which is the reference the recall must be
judged against.

---

## 2026-09-20 — §8. 15 stations, the ceiling, and the two-stage threshold

Whole 2005 episode, 15 stations (the 11 clean + JRBC, LZB, VGZ, SNB — the latter four were
withheld only by the old blind-**Gates** design, which does not govern this reproduce loop and
never withheld *families*), identical literal chain.

| stage | 11 stations | **15 stations** |
|---|---|---|
| autocorrelation detections above 8×MAD | 143,441 | 124,496 |
| true same-family pairs | 8.4 % (chance 0.92 %) | **12.9 % (chance 1.17 %)** |
| kept families seen at the pair stage | 88 | **92 of 96** |
| groups ≥3 members | 1,397 | 1,514 |
| groups at purity ≥0.5 / families covered | 337 / 75 | 393 / **77** |
| grown templates passing Bostock's ≥30 % bar | 137 | **177** |
| **distinct kept families rediscovered** | 50 | **58 of 96** |

Per-stage loss, 15 stations: **92 seen → 77 with a pure group → 60 at purity ≥0.8 → 58 grown.**
15 families have true autocorrelation pairs but never form a pure group; 23 form one and are lost
at growth (all small: median seed 4 nodes).

**The ceiling** (`k_ceiling.py`): templates built directly from up to 200 corrected-time catalogued
events per kept family, spread over the window, pushed through the **identical** matched filter at
8 × MAD on the same 15 stations and the same 24 days:

| | templates passing ≥30 % | distinct families | recall_all | recall_60 | recall_loud | purity |
|---|---|---|---|---|---|---|
| ceiling, pass 0 | 93 / 96 | **89** | 0.310 | 0.222 | 0.479 | 0.915 |
| ceiling, pass 1 | 89 / 96 | 85 | 0.360 | 0.259 | 0.542 | 0.846 |
| **grown from autocorrelation seeds** | 177 / 1,457 | **58** | 0.020 | 0.007 | 0.041 | 0.556 |

So with **perfect** seeding this network recovers 89 of 96 kept families and only **31 % of each
family's catalogued events** (22 % of the 60-s-independent subset, 48 % of the loud subset) at
Bostock's own threshold. Recall is therefore a station-count/SNR fact of this 15-station bed, not a
failure of the chain — **the ≥30 % one-to-one purity bar, not recall, is the right success
measure**, and the discovery figure to quote is 58 of the ceiling's 89 (65 %).

**A looser first pass makes it worse, not better.** Running pass 0 at 6 × MAD and 8 × MAD
afterwards (the "gather then tighten" idea) gives **0** families at pass 0 and 13 at pass 1, versus
58 with 8 × MAD throughout: the extra members are false and the stack never converges. Bostock's
single threshold is right. (`results/repro/grow15L`, logs `k_score15L_it*.log`.)

---

## 2026-09-20 — §9. THE BLIND CHECK (30 sealed families), and a scorer defect found first

**A void run, named.** The first blind attempt reported 0 of 30. It was scoring the 2005 templates
against **2010** events: `k_repro4_queue.sh`'s `run_bed` helper `export`s `K_YEAR/K_DAY0/K_DAY1`
and `unset`s them at the end, but the 2010 permanent bed aborted early (`KeyError: 'PFB'` —
PFB/YOUB/OCP/STW are absent from the station table), so the `unset` never ran and the blind stage
inherited `K_YEAR=2010`. The hidden truth file spans t = 1.2819e9–1.2840e9 (2010-08-16…09-08)
against a detection axis of 1.1253e9–1.1273e9 (2005). Every one of 1,457 templates returned
`n_match = 0` by construction. A second, smaller defect was found in the same audit: `k_s1_loud`
built sub-second corrections for only **11 of 30** hidden families (its reference bootstrap needs
≥10 loud events after cutting), and the scorer silently scored only those 11 — the log's "of 11".

**Neither defect touched the engine.** The corrections are not needed for this scoring at all:
the one-to-one tolerance is ±2 s and whole-second quantisation is ≤0.5 s. Control: scoring the
**kept** families with raw whole-second times gives **63 of 100** (vs 58 of 96 with corrected
times) — the scorer works on raw times.

**The blind result, on the unchanged K_REPRO3 templates:**

| | | |
|---|---|---|
| **HIDDEN families rediscovered** | **22 of 30** | 97 of 1,457 templates pass the ≥30 % bar |
| **ISOLATED stratum** | **6 of 10** | families 8, 10, 15, 62, 87, 115 |
| per-family medians | purity 0.655 | recall_all 0.024, recall_60 0.015 |
| kept, for comparison | 58 of 96 (corrected) / 63 of 100 (raw) | |

73 % of the sealed families are recovered by seeds the engine found for itself, against 60–63 % of
the development families. The blind number is consistent with — slightly better than — the
development number, which is what an honest blind check should look like. This is the first and
only blind measurement; nothing was retuned after seeing it.

## 2026-09-20 — §10. Was "template-free discovery is impossible on the permanent network" true?

The identical literal chain, same code, no parameter changed, on the **2012 permanent-station bed**
(`lfe_catalog_repair/cache20hz`, 52 days 2012-08-15…10-05, 15 permanent stations: PGC, VGZ, SNB,
NLLB, LRIV, B001, B003, B004, B005, B006, B007, B009, B010, B011, B926; 77 of the 100 kept
families have events in that window):

| | 2005 POLARIS array bed (15 stations) | **2012 permanent bed (15 stations)** |
|---|---|---|
| autocorrelation detections above 8×MAD | 124,496 | 2,966,785 |
| **true same-family pairs** | **12.85 %** (chance 1.17 %) | **0.05 %** (chance 0.03 %) |
| families seen at the pair stage | 92 of 96 | 58 of 77 (at 1.7× chance) |
| groups ≥3 members | 1,514 | 835 |
| groups at purity ≥0.5 / families | 393 / **77** | 17 / **13** |
| groups at purity ≥0.8 | 205 | **4** |
| **distinct families rediscovered (≥30 %)** | **58 of 96** | **3 of 77** |
| ceiling (catalogue-built templates, same filter) | 89 of 96 | **33 of 77** |

**The answer is: essentially yes, for the permanent network — and the reason is the pair stage, not
the machinery.** On the permanent bed the growth stage still works when handed a correct template
(the ceiling recovers 33 of 77 families at recall_all 0.33, recall_loud 0.61), but the
autocorrelation is at **1.7× chance** (0.05 % vs 0.03 %) instead of 11× (12.85 % vs 1.17 %), so
almost no clean group ever forms: 4 groups at purity ≥0.8 against 205 on the array bed. Discovery
collapses by a factor of ~20 (3 vs 58) while the ceiling only halves (33 vs 89). The earlier
project verdict was right about the permanent network and wrong only in being stated as a fact
about the method: with the array, the same code finds 58 of 96 kept and 22 of 30 sealed families.
(The 2010 permanent bed aborted on a station-table gap and is being re-run with the 8 stations that
have coordinates.)

---

## 2026-09-20 — §11. The blind denominator, the classification, and the 2010 bed

### 11.1 The hidden ceiling — 22/30 has the same shape as 58/96

Catalogue-built templates for the 30 sealed families (up to 200 corrected-time events each,
spread over the window) through the identical 15-station matched filter at 8 × MAD
(`results/repro/ceiling_hid2005/`, score `ceiling_hid2005_score.csv`, log `k_ceiling_hid_score.log`):

| stratum | discovered (autocorrelation seeds only) | ceiling (catalogue-built templates) | fraction of ceiling |
|---|---|---|---|
| **KEPT (development)** | 58 of 96 | 89 of 96 | **65 %** |
| **HIDDEN (sealed, blind)** | **22 of 30** | **29 of 30** | **76 %** |

Ceiling medians for the hidden set: purity 0.910, recall_all 0.266, recall_60 0.221,
recall_loud 0.429. Blind result with *corrected* hidden times as a confirmation:
**22 of 30**, 94 templates, purity 0.627 (`logs/k_blind_corr.log`) — identical to the raw-time
number reported in §9, so the blind figure does not depend on the time correction.

### 11.2 Classification of all 1,457 grown 2005 templates

`scripts/k_classify.py` → `results/repro/classify15.csv` (locations from
`results/repro/loc15_it0.csv`; hidden labels from the **correct** `score_hidden2005_raw.csv`;
stack cc against both ceilings, `ceiling15b/tpl_it0.npz` and `ceiling_hid2005/tpl_it0.npz`):

| class | templates | distinct |
|---|---|---|
| matches a KEPT family (≥30 %) | 177 | 58 families |
| **matches a HIDDEN family (≥30 %)** | **64** | **18 families** |
| candidate new family (located, shuffled-pick control passed, \|Δz\| ≤ 8 km of Slab2, ≥3 days) | **18** | 18 |
| unlocatable (<5 stations with a coherent stack, or control failed, or off-interface) | 1,201 | — |

(The hidden-class count of 18 is smaller than the 22 of §9 because classification also requires a
location on ≥5 stations; 4 of the 22 blind-matched families are carried only by templates that do
not locate.)

### 11.3 The candidate new families — and why most of them are not new

`results/repro/candidate_new_families.csv`, all 18, with the nearest catalogued family:

| grp | n_det | days | n_sta | lat | lon | dep | Slab2 | Δz | rms | rms_shuf p05 | nearest kept (km) | nearest hidden (km) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 400 | 2505 | 24 | 6 | 47.90 | −124.82 | 22.0 | 24.3 | −2.3 | 1.41 | 3.59 | 22 (80.5) | 18 (78.0) |
| 843 | 67 | 17 | 7 | 48.14 | −124.90 | 18.0 | 24.3 | −6.3 | 3.77 | 4.18 | 93 (65.7) | 18 (71.7) |
| 80 | 66 | 17 | 5 | 48.08 | −124.86 | 18.0 | 22.8 | −4.8 | 6.49 | 6.71 | 93 (69.3) | 18 (71.0) |
| 128 | 26 | 12 | 5 | 47.90 | −122.82 | 18.0 | 16.1 | 1.9 | 7.17 | 7.68 | 16 (47.0) | 8 (40.2) |
| 1424 | 18 | 10 | 6 | 48.02 | −122.82 | 18.0 | 16.0 | 2.0 | 0.96 | 3.21 | 16 (35.5) | 8 (29.5) |
| 556 | 17 | 4 | 5 | 48.38 | −123.58 | 18.0 | 21.2 | −3.2 | 0.79 | 1.04 | 38 (4.3) | 29 (4.0) |
| 295 | 13 | 5 | 9 | 48.44 | −123.46 | 18.0 | 14.1 | 3.9 | 1.84 | 3.42 | 35 (5.8) | 69 (9.2) |
| 1248 | 12 | 5 | 11 | 48.68 | −123.94 | 18.0 | 19.3 | −1.3 | 3.41 | 3.76 | 103 (4.3) | 95 (16.6) |
| 1375 | 11 | 9 | 5 | 48.23 | −122.82 | 18.0 | 15.8 | 2.2 | 4.67 | 5.31 | 16 (21.3) | 8 (20.7) |
| **534** | 11 | 4 | 8 | 48.20 | −123.78 | 38.0 | 31.5 | 6.5 | 0.57 | 1.02 | 13 (8.1) | 1 (2.3) |
| **1215** | 10 | 4 | 5 | 47.90 | −123.66 | 32.0 | 35.9 | −3.9 | 0.59 | 0.95 | 0 (24.6) | 2 (31.3) |
| **912** | 10 | 4 | 6 | 48.50 | −124.18 | 28.0 | 29.7 | −1.7 | 0.84 | 1.34 | 84 (5.8) | 58 (8.4) |
| **889** | 10 | 4 | 5 | 48.59 | −123.90 | 32.0 | 27.3 | 4.7 | 0.29 | 0.71 | 91 (1.3) | 82 (9.8) |
| 527 | 9 | 3 | 9 | 48.29 | −123.46 | 18.0 | 21.9 | −3.9 | 1.61 | 2.71 | 20 (2.2) | 15 (3.7) |
| **247** | 9 | 5 | 7 | 48.23 | −123.14 | 28.0 | 28.0 | 0.0 | 0.54 | 1.88 | 16 (6.3) | 8 (3.7) |
| **745** | 8 | 4 | 6 | 48.59 | −123.58 | 26.0 | 29.7 | −3.7 | 0.62 | 1.36 | 101 (4.0) | 78 (7.3) |
| **822** | 8 | 3 | 5 | 48.44 | −123.82 | 30.0 | 36.6 | −6.6 | 0.10 | 0.47 | 55 (2.3) | 60 (4.5) |
| 1268 | 6 | 3 | 6 | 48.50 | −124.06 | 18.0 | 18.4 | −0.4 | 0.42 | 1.09 | 57 (4.0) | 58 (2.3) |

**Read this conservatively. Ten of the eighteen rail at 18.0 km, the shallowest node of the search
grid**, which means their depth is unconstrained, not shallow; they should be discounted, and four
of them (400, 843, 80, 128/1375) also sit 20–80 km from any catalogued family, off the family field
entirely, where a persistent local noise source is the likelier explanation (grp 400 fires on all
24 days with 2,505 detections — that is a site, not an LFE family). **Seven have interior depths
(26–38 km), lie within 8 km of Slab2, pass the shuffled-pick control with a clear margin, and repeat
on 3–5 days: 534, 1215, 912, 889, 247, 745, 822.** Six of those seven lie **1.3–8.1 km from a
catalogued family**, i.e. they are most plausibly unlisted neighbours inside Bostock's own field
rather than a new region; only 1215 (24.6 km from the nearest kept, 31.3 km from the nearest
hidden) is spatially distinct. None correlates with any kept or hidden ceiling stack above
cc = 0.15, but that test has no power here (§7: same-family stack cc is only 0.072).

**Honest count: 7 candidate new northern families, of which 6 are probable unlisted neighbours and
1 (grp 1215, 47.90 N / −123.66, 32 km) is spatially distinct.** None is claimed; the confirmation
is a re-fire test (§11.5).

### 11.4 Three beds, one table

| bed | stations | days | detections | true pairs | chance | pure groups (≥0.5) / families | **discovered** | ceiling |
|---|---|---|---|---|---|---|---|---|
| **2005 POLARIS array** | 15 | 24 | 124,496 | **12.85 %** | 1.17 % | 393 / 77 | **58 of 96** | 89 of 96 |
| **2012 permanent** | 15 | 52 | 2,966,785 | **0.05 %** | 0.03 % | 17 / 13 | **3 of 77** | 33 of 77 |
| 2010 permanent | 8 | 34 | — | — | — | — | — | — |

The 2010 bed crashed twice: first on `KeyError: 'PFB'` (PFB/YOUB/OCP/STW have no coordinates in the
station table), then because a 2010 station-day is cached with a single channel and a `(1, N)`
array reached the autocorrelation. `k_common.load_day` now rejects non-3-C station-days and the bed
is re-running on the 8 stations that have coordinates (`scripts/k_repro6_queue.sh`, markers
`K_REPRO6_DONE|FAILED`). Its row is left blank rather than guessed.

### 11.5 How a candidate would be confirmed

POLARIS-BC ran 2004-04 → 2005-10, and the 2004-07 ETS episode (59,659 catalogued events over 18
busy days) is downloadable — about 900 station-days on the same stations. A candidate is confirmed
if its 2005 template, unchanged, re-fires in 2004-07: ≥20 independent detections on ≥3 days at the
same 8 × MAD rule, with a location that reproduces to within the 2005 uncertainty. That is a test
on data the template has never seen, in a different episode, and it is the only evidence that would
turn "candidate" into "new family". It has not been run.

### §11.6 (orchestrator, 2026-09-21T00:13:44Z) — 2010 permanent bed
2010 PERMANENT bed (8 stations with coordinates, 34 days; K_REPRO6, logs/k_repro6_queue.log): 1,181 groups, 19 at purity ≥0.5 covering 15 kept families, 7 at ≥0.8; autocorrelation discovery **5 of 60** kept families in window; catalogue-built ceiling **50 of 60**. Same pattern as 2012 (3/77 vs 33/77): the permanent network runs known families but does not find them.
