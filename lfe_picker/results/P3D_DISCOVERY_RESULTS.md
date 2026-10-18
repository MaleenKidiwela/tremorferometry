# P3D — the discovery look: template-free repeater search on the PQ-era dense network

2026-09-19, `lfe_picker`, agent P3D. Authority: coordinator decision to spend the discovery look on
pooled **S18 + S21 + S22** with the V2 (PQ-era) network, per Merlin P3D § Q3 / Step 6 and the P3
Step-3 nulls, as **one look** recorded in `results/referee_ledger.csv`.

Pre-registration: **`results/p3d/discovery_freeze.json`**, written and hashed **before any pass**
(`31da03e3fdc5113e53ce297e9efb1712d74ca4e0e3be153a11ca749890bfdc00`). Every bar, clause, null and
threshold rule below is read from that file at run time; none is re-decided here.

---

## 1. What was frozen, and the two forks the coordinator settled

### 1.1 Networks

| | id | rule | sites | days | episodes | windows | pairs |
|---|---|---|---|---|---|---|---|
| **PRIMARY** | NET-B | the V2 ≤60 km site set **minus JRDRB** | 7 | 129 | S18 17 + S21 50 + S22 62 | 381,000 | 7.26e10 |
| SECONDARY | NET-D | PGC cluster + the two best-registering others by Step-2 ρ | 3 | 162 | S18 18 + S21 82 + S22 62 | 762,600 | 2.91e11 |

NET-B sites: **PGCq, VGZ, SYMB, BUCKS, LHLYB, ALBHB, DAOB**. Withheld and never summed:
**B926q (B926+CLRS), GOBB, LOPEZ, TURTL**. Co-located pairs merged; PQ HHZ dip +90 honoured.

NET-D is named by the ρ rule as **PGCq + LHLYB (ρ/ρ_PGC = 1.086) + DAOB (0.826)**. *Deviation
recorded in the freeze*: the instruction said "over all 203 days", but LHLYB is cached on 162 of the
203, so NET-D realises **162 days**. The two 203-day alternatives (PGCq+BUCKS+DAOB, PGCq+SYMB+DAOB)
would violate the ρ rule by substituting a worse-registering site, so the ρ rule was honoured and
the shortfall stated rather than the rule bent.

### 1.2 Why S18 forced a choice at all

Requiring all 8 V2 ≤60 km sites leaves **112 of 203 pooled days and zero of S18's 59** — LHLYB
delivers 18/59 S18 days (2024-07 is genuinely `nodata`) and JRDRB 0/59. Dropping JRDRB, the station
already documented as non-delivering (24 % of in-epoch SSE days, 0 calibration days, worst 4a
residual), restores an S18 readout at the cost of one site.

### 1.3 Candidates, threshold, truth

* **Candidates** — uniform 12-s windows at 6-s hop on PNSN tremor hours within 50 km of the
  130-family field. **The PNSN catalog ends 2026-05-27 23:52**, so S22 contributes tremor hours from
  **57 of its 62 padded days**; the last 5 days are cached but carry no tremor-hour label and yield
  no candidates. Tremor hours: S18 658, S21 838, S22 177.
* **Threshold** — the null quantile giving **≤30 expected false edges** over the realised pair count
  (the rule is frozen; the quantile is recomputed from the run).
* **Truth T1** — hp events (`p_real ≥ 0.8`) of the 30 hidden families in the three episodes:
  **18 of 30 scoreable** (≥20 T1 events), **6 of 10 in the isolated stratum**.
* **Bars** — PRIMARY: ≥2/3 of the scoreable isolated stratum (bkeys **8, 10, 15, 87, 115, 119**).
  Two SECONDARY readouts, both explicitly *never bars*: (i) identity-only on the PQ-rich isolated
  families **12, 29, 62**; (ii) a Δsites-weighted rediscovery fraction over all 10.

### 1.4 The structural fact, fixed in advance

Recorded in the freeze before any pass, so the reading could not be chosen afterwards:
the **scoreable** isolated families gain **0, 1, 2, 1, 0, 0** new PQ sites within 45 km, while the
three isolated families that gain **5 each** (12, 29, 62) carry only **8, 8 and 4** T1 events and are
unscoreable. Medians: scoreable 44.5 T1 / 0.5 new PQ sites; unscoreable 6.0 T1 / 5.0 new PQ sites.
Spearman(n_T1, Δsites) = −0.17 (p = 0.63, n = 10) — **not significant**; the structural fact is what
is reported, not a correlation.

Consequence, also fixed in advance: **a PRIMARY failure is not evidence that the PQ geometry fails
to help**, because the primary bar is carried by families PQ barely touches. A primary failure with
a positive SECONDARY-(i) readout is the signature of geometry helping where it exists but truth
being absent; a primary failure with a null SECONDARY-(i) is the stronger negative.

---

## 2. The engine, and the two things that had to be right before spending the look

`scripts/s_discover_pairs.py` — the P3 Step-1a statistic (per-site MAD-z of the 24-s 3-C normalised
cross-correlation, max over ±3 s lag, summed over the configuration's sites, 27-s overlap exclusion,
index × shift circular-shift null) applied to **uniform tremor-hour windows** instead of catalogue
detections, emitting the **edge list** the discovery chain needs.

**(a) The fused accumulator.** The inherited code held one `[n_i, n_j·121]` fp32 array *per site* and
summed them — 7–8× the memory traffic for no arithmetic. Fusing to a single accumulator was
validated against the per-site path on random data: **max |difference| 0.000e+00**, argmax lag
identical on **100 %** of pairs. Not "close" — bit-identical.

**(b) j-blocking, after an OOM.** The first full pass died with CUDA OOM: some tremor days carry
~11 k windows, and the j-side unfold costs n_j·121·1440 fp16 *per station* (~3.9 GB each) plus an
11 GB fp32 output tile. Days are now split into sub-blocks of ≤3,000 windows. Regression-checked on
a 30 k-window probe:

| | pre-blocking | post-blocking | |
|---|---|---|---|
| pairs | 446,870,460 | 446,870,460 | same |
| threshold | 29.682254791259766 | 29.682254791259766 | same |
| quantile | 6.71335e-08 | 6.71335e-08 | same |
| runtime | 162 s | **67 s** | 2.4× faster |

Measured full-scale rate on NET-B: **2.8e7 pairs/s**, i.e. **~45 min per 7.26e10-pair pass** — the
30 k probe's tiny GEMMs had implied 4 h and were a floor, not a projection.

**(c) The scorer's matcher, proven before the data landed.** The one-to-one ±2 s matcher that decides
every rediscovery clause was exercised against known truth (bkey 10, 161 T1 events in the pooled
episodes): the truth against itself matches **161/161**; jittered by ±1 s, still **161**; shifted by
+30 s, **4** (the residual rate LFE bursts produce by accident); a *different* family's events, **0**.
So a rediscovery claim cannot come from a mis-wired matcher, and the accidental-match floor is
measured rather than assumed.

---

## 3. The three arms — inputs

| arm | trigger | candidates | threshold rule |
|---|---|---|---|
| **A pairwise** | none (uniform grid) | 381,000 windows | ≤30 expected false edges |
| **B seed-and-grow** | top 1 % by **median** per-site S-energy z | **3,809 seeds** of 380,943 scored windows | 1e-4 site-sum quantile |
| **C CF-triggered** | frozen P2 CF (`p2nnpu`), unknown-station embedding, dip honoured | 282,019 peaks (720/station-day at 7 sites, declustered 6 s) → **145,661** merged candidates | ≤30 expected false edges |

ARM B's seed statistic is the **median** over sites, not the sum — pre-registered, because Step 3
showed the site sum is dominated by the PGC cluster and a sum would let one loud site nominate
every seed. All three arms run on **identical data and identical nulls**, as required.

---

## 4. Results — all three arms, one table

| | ARM A pairwise | ARM B seed-and-grow | ARM C CF-triggered |
|---|---|---|---|
| candidates | 379,686 windows | 3,809 seeds × 379,686 | 144,893 CF candidates |
| pairs | **7.21e10** | 9.15e8 | 1.05e10 |
| threshold rule | ≤30 expected false edges | fixed 1e-4 quantile | ≤30 expected false edges |
| threshold | **64.795** (q = 4.16e-10) | **17.084** (q = 1e-4, rank 91,545) | **36.048** (q = 2.86e-9, rank 30) |
| **edges above threshold** | **455** | 327,763 | **873** |
| components ≥20 members | 2 (largest 132) | 3 (largest 45) | 2 (largest 125) |
| **components passing S excess ≥+3 dB at ≥3 sites** | **0 / 2** | **0 / 3** | **0 / 2** |
| PRIMARY: scoreable isolated rediscovered | **0 / 6** | **0 / 6** | **0 / 6** |
| SECONDARY (i) identity-only on 12, 29, 62 | none | none | none |
| SECONDARY (ii) Δsites-weighted fraction | 0.00 | 0.00 | 0.00 |

**No arm produced a single cluster that survives the chain.** The best site excess reached by any
component, in any arm, was:

| arm | cluster | members | span | per-site S excess (dB) |
|---|---|---|---|---|
| A | 13840 | 132 | 392 d | ALBHB 1.48, BUCKS 2.33, DAOB 0.54, LHLYB 2.69, PGCq 1.82, **SYMB 2.80**, VGZ 0.54 |
| B | 254611 | 32 | 1.7 d | ALBHB −1.46, BUCKS 0.47, DAOB 0.33, LHLYB 1.44, **PGCq 2.60**, SYMB 1.02, VGZ −0.30 |
| **C** | **10640** | **62** | 8.7 d | **ALBHB 3.27**, BUCKS 0.15, DAOB 0.41, **LHLYB 3.61**, PGCq 0.04, SYMB 0.86, VGZ −0.05 |

ARM C's cluster 10640 is the closest approach in the whole study: **two** sites (ALBHB, LHLYB —
both PQ) clear +3 dB, against a clause requiring three. Note which sites they are: the two PQ
stations nearest the family field. Cluster 13840 (A) and 3739 (C) span 392 days — the full pooled
window — so they are not bursts; they are the kind of long-lived, low-excess association a
site-coherent noise process produces.

**ARM C retention.** The pre-registered rule keeps ARM C only if it rediscovers ≥2 more scoreable
isolated families than ARM A. It rediscovered **0 vs 0** → **ARM C is NOT retained.**

**N-shift claim-rule nulls: not triggered.** The 5-seed N-shift null exists to separate "discovery"
from "candidate" for gate-PASS clusters matching no Bostock family. **No cluster in any arm reached
the gates**, so there was nothing to adjudicate and the 5 × 45 min was not spent. The conditional is
recorded, not skipped.

## 5. The kill test — is the construction broken, or is it starved?

A null result is only evidence if the machinery can link families it is *given*. Each arm's own
edge list was asked: of the windows coinciding (±2 s) with a **KEPT** family's hp event, do their
edges land on another hp event of the **same** family?

| arm | windows on a kept-family event | edges with both ends labelled | same-family | chance | enrichment | P(≥ observed) |
|---|---|---|---|---|---|---|
| A | 507 / 379,686 (0.13 %) | 0 | 0 | 0.0285 | — | — |
| **B** | 507 / 379,686 | **15** | **5** | 0.0285 | **11.7 ×** | **4.4e-5** |
| C | 225 / 144,893 | 1 | 1 | 0.0285 | 35 × (n = 1) | — |

**The construction is not broken — it is starved.** ARM B's 5 same-family edges out of 15 is
11.7 × chance at p = 4.4e-5, the same sign and the same order as P3's 8–38 × enrichment on
catalogue candidates. The linkage physics is real and reproduces. What is missing is *quantity*:
5 same-family edges across 4 families, where building one family needs tens.

## 6. The finding that decides how the null must be read

The candidate grid can only represent a fraction of the truth, and that fraction is **not uniform
across the families the bar is made of**. Of hp events in the pooled episodes, the fraction lying
within ±2 s of any candidate window is **23.2 %** (kept), **23.5 %** (all hidden), **18.4 %**
(isolated hidden). Per family of the PRIMARY bar:

| bkey | T1 events | representable | where its truth lives |
|---|---|---|---|
| **8** | 27 | **0 (0 %)** | S18 17, S21 10 |
| **10** | 161 | **0 (0 %)** | S18 142, S21 19 |
| **15** | 46 | **0 (0 %)** | S18 45, S21 1 |
| 87 | 101 | 65 (64 %) | S21 65, S22 24, S18 12 |
| 115 | 43 | 33 (77 %) | S21 31, S18 12 |
| 119 | 35 | 22 (63 %) | S21 20, S22 8, S18 7 |

**Three of the six families in the PRIMARY bar have ZERO representable truth**, because their events
live in S18 and **NET-B keeps only 17 of S18's 59 days** (LHLYB and ALBHB are not cached on the
other 42). So the PRIMARY bar — ≥2/3 of six, i.e. ≥4 — was **capped at 3/6 before any detector ran**.

> **The PRIMARY bar was unreachable by construction on NET-B. Its failure is therefore not evidence
> about the detector, the network's geometry, or the PQ stations.**

This was not knowable at branch-decision time: it required the realised candidate grid. Measured
across the three branches that were on the table:

| configuration | days | S18 days | representable T1 per PRIMARY family | families reachable |
|---|---|---|---|---|
| **NET-B (chosen PRIMARY)** | 129 | 17 | 0, 0, 0, 65, 33, 22 | **3 / 6** |
| NET-D (ρ rule) | 162 | 18 | 10, 19, 0, 73, 43, 29 | 5 / 6 |
| max-coverage 3 (PGCq+BUCKS+DAOB) | 201 | 57 | 27, 161, 45, 85, 43, 31 | **6 / 6** |

The branch chosen on geometry and registration quality was the **worst of the three** on the axis
that turned out to decide the test, and the branch rejected for violating the ρ rule was the only
one with full truth coverage. That is the single most useful thing this look produced.

## 7. Verdicts

* **PRIMARY (NET-B, all three arms): NOT MET — 0/6 — and uninterpretable**, because only 3 of the 6
  bar families had any representable truth.
* **SECONDARY (i)**, identity-only on the PQ-rich isolated families 12, 29, 62: **none rediscovered**
  in any arm. Descriptive only, as pre-registered — with 8, 8 and 4 truth events it never had power.
* **SECONDARY (ii)**, Δsites-weighted fraction over the 10 isolated families: **0.00** in all arms.
* **ARM C**: not retained (0 vs 0 against a ≥2 margin).
* **ARM B kill test**: the construction links kept families at **11.7 × chance (p = 4.4e-5)** but
  produces only 5 such edges — real physics, ~10²–10³ short, exactly P3's conclusion on a different
  decade, a different network and template-free candidates.
* **N-shift claim-rule nulls**: not triggered (no gate-PASS cluster).
* **NET-D**: **not run — engine.** At 2.91e11 pairs and the measured 2.8e7 pairs/s it is ~2.9 h per
  pass, above the ≤2 GPU-h/pass condition. Given § 6 it is also the configuration that *should* have
  been primary (5/6 reachable vs 3/6), so this is the first thing to run if the look is revisited.

**Reading fixed in advance (freeze § c).** A primary failure with a null SECONDARY-(i) was
pre-registered as "the stronger negative". It is, however, subordinate to § 6: on NET-B the primary
test could not have passed, so the honest summary is **not** "the dense PQ network cannot discover
families" but **"this configuration could not test that, and the one that could was not run."**

---

# EXPLORATORY / POST-HOC RE-RUN — 2026-09-19

**Status: EXPLORATORY, NOT A LICENCE.** `discovery_freeze.json → amendment_exploratory_2026-09-19`
(sha256 `c9034088…`), written before any of these passes. **The ledger row P3D-DISC remains SPENT
and CLOSED** — nothing below reopens the discovery look. Reason for the re-run: the pre-registered
PRIMARY network NET-B left 3 of the 6 scoreable-isolated families with **zero** representable truth,
so the pre-registered look could not test the question it was built for. Every other rule —
candidate definition, threshold rule, nulls, gates, scorer, both SECONDARY readouts, the kill test —
is identical; only the station set changes.

## E1. Reachable truth per configuration — *stated before the results*

| config | sites | stations | days | S18 / S21 / S22 | tremor h | windows | pairs | **PRIMARY families reachable** | isolated truth representable | all hidden truth |
|---|---|---|---|---|---|---|---|---|---|---|
| **NET-B** (pre-registered) | PGCq VGZ SYMB BUCKS LHLYB ALBHB DAOB | 10 | 129 | 17 / 50 / 62 | 635 | 381,000 | 7.26e10 | **3 / 6** | 121/434 (**28 %**) | 445/1303 (34 %) |
| **NET-D** (ρ rule) | PGCq LHLYB DAOB | 6 | 162 | 18 / 82 / 62 | 1,271 | 762,600 | 2.91e11 | **5 / 6** | 183/434 (**42 %**) | 738/1303 (57 %) |
| **NET-M** (max-coverage) | PGCq BUCKS DAOB | 6 | 201 | 57 / 82 / 62 | 1,646 | 987,600 | 4.88e11 | **6 / 6** | 412/434 (**95 %**) | 1264/1303 (**97 %**) |

Per PRIMARY family (representable T1 / total T1):

| config | 8 | 10 | 15 | 87 | 115 | 119 |
|---|---|---|---|---|---|---|
| NET-B | **0**/27 | **0**/161 | **0**/46 | 65/101 | 33/43 | 22/35 |
| NET-D | 10/27 | 19/161 | **0**/46 | 73/101 | 43/43 | 29/35 |
| **NET-M** | **27**/27 | **161**/161 | **45**/46 | 85/101 | 43/43 | 31/35 |

NET-M removes the defect that made the pre-registered look untestable: **95 %** of isolated-hidden
truth is representable, against 28 % on NET-B.

## E2. Results — every configuration and arm

| config | arm | pairs | threshold | edges | components ≥20 | **pass S excess** | PRIMARY | SEC (i) | SEC (ii) |
|---|---|---|---|---|---|---|---|---|---|
| NET-B | A | 7.21e10 | 64.795 | 455 | 2 | 0 | 0/6 | 0 | 0.00 |
| NET-B | B | 9.15e8 | 17.084 | 327,763 | 3 | 0 | 0/6 | 0 | 0.00 |
| NET-B | C | 1.05e10 | 36.048 | 873 | 2 | 0 | 0/6 | 0 | 0.00 |
| **NET-M** | A | **4.86e11** | 33.79 | 1,469 | 5 | **0** | 0/6 | 0 | 0.00 |
| **NET-M** | B | 4.54e9 | 8.633 | 1,593,057 | **percolated** | — | — | 0 | — |
| **NET-M** | C | 1.91e10 | 23.349 | 310 | 2 | **0** | 0/6 | 0 | 0.00 |
| NET-M | B, z≥16 *sensitivity* | 4.54e9 | 16.0 | 2,429 | 9 | 0 | 0/6 | 0 | 0.00 |
| **NET-D** | A | 2.90e11 | 26.122 | 2,167 | 7 | **1** | 0/6 | 0 | 0.00 |

**Rediscovery by stratum: zero, in every configuration and every arm.** No isolated-hidden family
(scoreable or not), no mid or close-stratum hidden family, was rediscovered. SECONDARY (i) on the
PQ-rich families 12, 29, 62: none. SECONDARY (ii) Δsites-weighted fraction: 0.00 throughout.

**NET-M ARM B percolated.** At the pre-registered 1e-4 quantile (threshold 8.633) the graph fuses
into a single 404,165-window component — 41 % of all candidates — and nothing falls in the
admissible [20, 20000] size range. That is a *result*: on a 3-site network the site sum has a
heavier tail and the pre-registered quantile is too permissive. A descriptive sensitivity at z ≥ 16
de-percolates it into 9 components, none of which passes the S-excess clause.

**A structural note on the S-excess clause.** It requires ≥ +3 dB at **≥3 sites**. NET-M and NET-D
have exactly 3 sites, so on them the clause demands **unanimity** — strictly harder than on the
7-site NET-B. Several clusters cleared +3 dB at one or two sites (NET-D 269205: LHLYB 6.28;
NET-D 5507: LHLYB 9.77; NET-D 26183: DAOB 4.55, PGCq 4.00; NET-M-B-sens 81255: DAOB 8.73,
PGCq 4.27) and were rejected on the third.

## E3. The one cluster that passed — NET-D 51599 — and why it is not a source

| | |
|---|---|
| members | 124 windows, span 591 d (S21 101, S18 13, S22 10) |
| **independent time episodes (>27 s apart)** | **33** — 91 of 123 gaps are ≤27 s at the 6-s hop |
| S excess in the chain | DAOB 6.87, LHLYB 3.83, PGCq 3.60 dB — all three sites, clause passed |
| best Bostock match | bkey 25, **1** one-to-one match — effectively unmatched |
| **Gate 1 at the withheld sites** (a3b, alignment-free, on the 33 independent episodes) | **0 / 5 pass**: B926 1.13 (null p99 1.98), CLRS 0.48 (2.20), GOBB −0.15 (2.82), LOPEZ 1.97 (2.90), TURTL 0.62 (2.92) |
| **Gate 0 identity** | PGC stack top \|cc\| vs the 130 2012 PGC templates = **0.231** (bkey 59), against the 0.524 chance level → not a known family |
| PGC stack on the 33 independent episodes | S-SNR **1.70 dB**, S excess **0.59 dB** — no S arrival |

The member count is inflated ~4× by window adjacency: the 27-s overlap exclusion forbids a *pair*
within 27 s, but two adjacent windows can still join the same component transitively through a
shared distant partner. The chain's S excess was measured on all 124 windows **after
`align_stack`**; the alignment-free test at the withheld sites, on the 33 genuinely independent
episodes, shows nothing at any site, and PGC — a site *in* the detection network — shows no S
arrival at all. **51599 is an alignment-and-adjacency artefact, not a source.**

**Consequently no gate-PASS cluster matching no Bostock family exists in any configuration, and the
5-seed N-shift claim-rule nulls were not triggered.** Nothing is reported as "discovery" or as
"candidate, nulls pending".

## E4. The kill test — where the exploratory run actually changed something

| config | arm | windows on a kept-family event | edges both-labelled | same-family | enrichment | verdict |
|---|---|---|---|---|---|---|
| NET-B | B | 507 / 379,686 | 15 | 5 | 11.7 × | starved |
| NET-M | A | 1,431 / 986,125 | 1 | 0 | — | no power |
| **NET-M** | **B** | **1,431 / 986,125** | **476** | **196 (41.2 %)** | **14.5 ×** | **LINKS KEPT FAMILIES** |
| NET-M | C | 328 / 144,893 | 0 | 0 | — | no power |
| NET-D | A | 891 / 761,984 | 0 | 0 | — | no power |

This is the substantive gain from the re-run. On NET-M the seed-and-grow construction links kept
families **196 times** at **14.5 × chance** (41.2 % of both-labelled edges same-family, against a
2.85 % chance rate) — where on NET-B it managed 5. Truth representability went 28 % → 95 % and the
kill test went from "starved" to "works". **The linkage physics is confirmed on a network that can
actually see the truth; what still fails is the step from linked pairs to a clustered, gated
family** — percolation at the pre-registered threshold, and no surviving S excess at three sites.

## E5. Verdicts

* **Rediscovery: 0 in every stratum, every configuration, every arm.** The PRIMARY bar is NOT MET on
  NET-M and NET-D as well — and on those two, unlike NET-B, it **could** have been met (6/6 and 5/6
  families reachable). That makes these the informative negatives the pre-registered run was not.
* **No gate-PASS cluster anywhere**; the single S-excess survivor (NET-D 51599) fails Gate 1 at all
  five withheld sites and Gate 0 identity. N-shift nulls not triggered.
* **The construction is confirmed to work and confirmed to be insufficient**: 14.5 × same-family
  enrichment on NET-M, and still no family assembled.
* **ARM C** is not retained under the pre-registered rule in any configuration (0 vs 0 everywhere).
* Two engine faults were found and fixed while scoring, both of which had turned a *result* into a
  crash: `components()` died when the graph percolated (one giant component above the size cap), and
  the chain wrote a column-less CSV when nothing survived. Both now emit typed, empty products and
  report the component-size distribution.

**What this run licenses.** Nothing, by construction — it is post-hoc. What it *establishes* is
narrower and more useful than the pre-registered run could have: on a network where 95 % of the
hidden truth is representable, template-free pair linkage recovers same-family structure at 14.5 ×
chance and still cannot assemble a family. The failure is now located at the **clustering and
gating** step, not at linkage and not at truth coverage — which is a different, testable target for
any future pre-registered attempt.
