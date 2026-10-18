# P0-O — oracle and referee freeze: RESULTS

2026-09-17, lfe_picker. Spec: `results/merlin_review_plan_2026-09-17.md` § "Prescription — P0-O" plus
the P0 label / holdout definitions. Everything below was produced inside `lfe_picker/`; `../lfe_extension`,
`../lfe_catalog_repair` and `../codex-LFE` were read only. Scripts are copies with an `o_` / `o1..o3`
prefix; nothing outside this folder was edited.

## 0. Headline

| item | verdict |
|---|---|
| **O1 referee freeze** | done — `results/referee_protocol.json`, `results/holdouts.json`, both hashed |
| **O1 2008/2009 download** | finished — 851 station-days (10 PB boreholes); **CN / C8 / UW are NOT available pre-2010** |
| **O2 P1a baseline** | **ACCEPT** — 2011 own recall 1.11× (T_core) / 1.09× (T_ext) the 2012 S15 reference, inside 0.7-1.3× |
| **O3 headroom oracle** | **DISCOVERY-ONLY** — energy-network recall / template recall = **0.00** at both pre-registered rates (bar: < 1.25×) |
| **O4 southern truth test** | DEFERRED, nothing queued — `results/o4_deferred.md` |
| **P0 label table** | done — `results/labels_v0.parquet`, 1,479,426 rows, 1,011,174 positive (event, station) S labels |

---

## O1. Referee freeze

`results/referee_protocol.json` — the protocol, with every input hashed.

| role | year | Bostock window | data on disk | budget |
|---|---|---|---|---|
| **FINAL REFEREE** | 2013 | 2013-09-13 .. 10-13 (14,873 events) | 2013-08-30 .. 10-28, 60 d × 19 stations | ONE look at the end of the phase — **spent once here (O2)** and recorded; not to be used again |
| **STAGE GATE** | 2011 | all 2011 events (18,663) | 2011-08-05 .. 09-25, 52 d × 19 stations | exactly one evaluation per stage gate — **one spent here (O2)** |
| development | 2010 | — | 2010-08-10 .. 09-13 | unlimited |
| development | 2012 | — | 2012-08-15 .. 10-05 | unlimited; all fitting and the oracle live here |
| training breadth | 2008 | — | 2008-05-01 .. 06-10 (downloaded) | training only |
| training breadth | 2009 | — | 2009-05-01 .. 06-15 (downloaded) | training only |

**Time base.** Per-family offsets for EVERY year = the frozen 2012 offsets
(`e5r_family_offsets_2012.csv`, 130/130 accepted); `t_origin = Bostock DB time + offset_s(bkey)`.
Never refit on 2010 / 2011 / 2013. Family map = exact (lat, lon, dep) match to
`e5r_families_2012.csv`; 0 Bostock rows unmapped over 2003-2013.

**sha256 (first 16 hex; full values in the json).**

| input | sha256_16 |
|---|---|
| `bostock2015_sloweq_FULL_2003-2013.csv` | `04cb6d76703e9fdc` |
| `e5r_family_offsets_2012.csv` | `abfdd198b5ee7a47` |
| `e5r_templates_2012.npz` | `64b45a5bf6a34c1a` |
| `frozen_operating_point.json` | `5dc7ea70640decfe` |
| `e5r_stations.csv` | `50a3999c6bbe8e4e` |
| `e5r_families_2012.csv` | `2bc48a64302cab50` |
| `e5r_weights_snr_2012.csv` | `a214c27d096cac97` |
| `e4_moveout_e5r_mv_build.csv` | `78da98f68fbff432` |

The templates / weights / moveout / stations / families hashes reproduce the `sha256_16` recorded
inside the frozen operating point — the referee is pinned to the same bytes Track C ran on.

### Holdouts (`results/holdouts.json`, seed 20260917)

30 hidden families, 10 per nearest-neighbour stratum (epicentral NN distance between Bostock family
hypocentres), drawn by `numpy.default_rng(20260917)` from families with ≥ 50 Bostock events in 2012
**and** ≥ 20 in 2011 (eligibility uses no property of our detector):

| stratum | eligible | picked | families (bkey) |
|---|---|---|---|
| isolated, NN ≥ 5 km | 19 | 10 | 8, 10, 12, 15, 29, 62, 87, 113, 115, 119 |
| NN 2.5-5 km | 30 | 10 | 1, 2, 4, 7, 18, 28, 59, 69, 88, 112 |
| close, NN < 2.5 km | 58 | 10 | 17, 26, 27, 43, 58, 60, 78, 82, 95, 125 |

NN distribution over the 130 families: median 2.34 km, p10 0.00, p90 5.92, 56 % below 3 km — the
stratification is necessary, not decorative. sha256 of the hidden list:
`4ac40235db8fee0c` (field `sha256_of_hidden_family_list`).
Held-out stations: **B926, NLLB, LRIV** — excluded from all training and all labels, evaluated only as
unseen stations (Merlin P4a).

### 2008 / 2009 download — finished, and a real limitation

`scripts/o_download_plan.py` (FDSN metadata 2008-01-01..2010-01-01) + `scripts/o_download_worker.py`,
≤ 8 concurrent requests, output in `lfe_picker/traces` and `lfe_picker/results/download_manifest.csv`.
1,566 station-days queued, **851 data / 715 nodata, 17.0 GB**, finished in ~6 min.

| network | result |
|---|---|
| PB boreholes (B001, B003, B004, B005, B006, B007, B009, B010, B011, B926) | **all 87 days per station** (B003 80, B005 75 — the rest `nodata` at the source) |
| CN (PGC, VGZ, SNB, NLLB, LZB, PFB, YOUB), C8 SOKB | **0 days** — IRIS returns HTTP 204 for every 2008/2009 day |
| UW LRIV | not installed until 2009-08-27 → `no_epoch` in both windows |

IRIS has the CN *metadata* for 2008/2009 but not the waveforms, and the NRCAN FDSN endpoint
(`earthquakescanada.nrcan.gc.ca/fdsnws`, and three alternates) is unreachable from this pod. **So the
2008/2009 training breadth is borehole-only** — 10 PB stations, no CN surface station. That matters
for P1/P2: a model trained with 2008/2009 in the mix sees no pre-2010 CN example, and CN is where most
of the template information lives (O3 below). Recorded, not worked around.

---

## O2. P1a baseline — the frozen template detector on 2011 and 2013

Cache built with the copied `o_cache_month.py` from `../lfe_extension/traces` (1,989 station-days,
2011-08 .. 2013-10); detector = `o_run_month.py`, an unmodified copy of `c_run_month.py` pointed at the
picker's cache and output directory. Because `epoch_of_month` has no network epoch before 2014, the
runs used the **frozen anchor T_core 8.75 / T_ext 8.25** and the frozen BH HF cut 3.6384 (PGC is BH in
both years, the epoch the cut was calibrated on). 15 stations (S15), 116 families with templates.

Scoring = Track D's recipe verbatim (`d1_score_epochs.one2one`): own family, ±2 s, one-to-one greedy by
time, HF veto applied, recall read off at the frozen thresholds. Denominator = **all** Bostock events of
the year inside the data window (2011: 18,663; 2013: 14,873; every event of both years falls inside the
referee sub-window, so the two scopes coincide).

| year | tier | events | recalled | **own recall** | 95 % CI | 2012 S15 reference | **ratio** | accept 0.7-1.3× |
|---|---|---|---|---|---|---|---|---|
| 2011 | T_core 8.75 | 18,663 | 3,022 | **0.1619** | 0.1567-0.1673 | 0.1463 | **1.107** | **PASS** |
| 2011 | T_ext 8.25 | 18,663 | 3,736 | **0.2002** | 0.1945-0.2060 | 0.1841 | **1.088** | **PASS** |
| 2013 | T_core 8.75 | 14,873 | 2,448 | **0.1646** | 0.1587-0.1706 | 0.1463 | 1.125 | (referee — recorded only) |
| 2013 | T_ext 8.25 | 14,873 | 2,952 | **0.1985** | 0.1921-0.2050 | 0.1841 | 1.078 | (referee — recorded only) |

No time-base or moveout problem: the failure mode the spec warned about (< 0.5×) is nowhere near, and
the 2011 and 2013 values agree with each other to 2 %. The slight excess over 2012 is expected — the
2012 reference is measured on the test half only, and the 2011/2013 windows are more ETS-dominated.

### False alarms (quiet family-day = no PNSN tremor in the family box and ≥ 3 days from any tremor day)

| year | tier | quiet family-days | FA / template-day (raw) | FA / template-day (catalog) | **network-wide false events/day** raw / catalog |
|---|---|---|---|---|---|
| 2011 | T_core | 575 | 0.0522 | 0.0487 | **6.0 / 5.6** |
| 2011 | T_ext | 575 | 0.1826 | 0.1600 | 21.0 / 18.4 |
| 2013 | T_core | 1,380 | 0.1072 | 0.1014 | **12.3 / 11.7** |
| 2013 | T_ext | 1,380 | 0.2848 | 0.2601 | 32.8 / 29.9 |

Network-wide rate = FA per template-day × the median number of observable families per day (115).
The 2013 T_core value, **12.3 events/day**, lands on Merlin's pre-registered 12/day almost exactly, which
independently validates that number. The pre-registered **116/day for T_ext does not reproduce**: the
frozen T_ext runs at ~30 events/day on 2013 quiet days (and ~32 on 2012). O3 therefore evaluates at both
12/day and 116/day as instructed, but 116/day should be read as "a much looser point on the curve",
not as the frozen T_ext.

*Only allowed look at 2013.* The 2013 numbers above are recorded and are not used to select, tune or
compare anything downstream.

Files: `results/trackO/{detections,family_days,stacks,run}_*.{parquet,npz,json}`,
`results/o2_p1a_{2011,2013}.json`, `o2_recall_{year}.csv`, `o2_recall_family_{year}.csv`,
`o2_fa_{year}.csv`; figure `figures/o2_p1a_baseline.png`.

---

## O3. Per-station headroom oracle (2012 development year) — **DISCOVERY-ONLY**

Two statistics on exactly the same station-day series, the same cache, the same placement, the same
site merge (`scripts/o3_oracle.py`, 364 s on the H100):

* **ENERGY** — the E4-style incoherent 3-C S-window energy excess, `log10 E(S−1..S+5 s) / E(the same
  6-s window 20 s earlier)`, 1-8 Hz at 20 Hz, MAD-normalised over the valid samples of that station-day;
* **CC** — the 3-C normalised correlation with that family's 2012 template, MAD-normalised per
  (family, station, day) — the statistic the frozen detector already uses.

Both are read at the same predicted template start, merged within the PGC cluster {PGC,B009,B010,B011}
and the LRIV cluster {LRIV,B005,B006,B007} by the mean of their MAD units, and summed over the 9 sites.
Positives: Bostock's 2012 **test-half, non-guard** events (13,381; 13,032 evaluable — the rest fall
within 2 s of a day edge or an invalid window). Negatives: 154,714 sampled origin times — quiet-day
(doy 228-242) and ETS-day ≥ 120 s from any Bostock event of any family. False alarms: declustered peaks
of the network series over 1,613 quiet family-days, converted to network-wide false events/day by × 115.

### The decision table

| network-wide quiet-day false events/day | ENERGY recall | CC recall (same combination) | frozen template recall | **energy / template** |
|---|---|---|---|---|
| **12/day** (pre-registered, = the frozen T_core) | **0.000** | 0.074 | **0.148** | **0.00** |
| **116/day** (pre-registered) | **0.000** | 0.117 | **0.241** | **0.00** |
| 1,000/day | 0.020 | 0.188 | beyond the recorded range | — |
| 10,000/day | 0.174 | 0.324 | beyond the recorded range | — |
| 100,000/day | 0.452 | 0.564 | beyond the recorded range | — |

**Pre-registered rule: energy-network recall < 1.25 × template ⇒ P2 is DISCOVERY-ONLY.** Measured ratio
**0.00** at both rates. The verdict is not marginal.

Two more readings of the same result:

* **Rank separation, not thresholds.** Network AUC of ENERGY against the negatives is 0.785 (quiet) /
  0.787 (ETS); CC is 0.995 / 0.996. So the energy statistic *does* carry real signal — positives sit at
  a median of 5.3 summed MAD units against ~0.0 for random times — it is simply far weaker, and its
  quiet-day tail is heavy (declustered network peaks: p99 14, p99.99 93, max 350, against a positive
  p90 of 15.5). Glitches and teleseisms, not LFEs, own the top of the energy distribution.
* **The price of matching the template.** To reach the frozen T_core recall of 0.148 the ENERGY network
  must run at **7,519 network-wide false events/day** (626× the template's 12/day); the CC network with
  the same naive combination needs 371/day (31×). The frozen detector's SNR weights and per-family-day
  network-MAD normalisation are worth a factor ~31 in false alarms by themselves.

### Single-station ROC (`figures/o3_station_roc.png`, `results/o3_station_roc.csv`)

| | median AUC vs quiet | median AUC vs ETS |
|---|---|---|
| ENERGY | 0.644 | 0.646 |
| CC | 0.697 | 0.696 |

The medians hide the structure that actually matters for a learned per-station model:

| station group | AUC energy | AUC template cc |
|---|---|---|
| near sites: PGC, B009, B010, B011, B926, NLLB | 0.62-0.67 | **0.72-0.87** |
| far boreholes: B001, B003, B004, B005, B006, B007, LRIV | 0.60-0.65 | **0.52-0.56** |
| CN surface: SNB, VGZ | 0.64-0.67 | 0.70 |

**At the seven far boreholes the source-agnostic energy statistic BEATS the family template** (≈ 0.64
against ≈ 0.53, barely above chance). That is where a learned per-station characteristic function has
something to add: not more recall at the near sites, but usable evidence at stations the template
cannot read. It is also exactly where the labels are weakest (§ P0: those stations' S times are
PREDICTED, not picked).

### The target a learned per-station model must beat

Single-station recall at 100 declustered peaks per day (energy: per station-day; cc: per template-day):

| sta | recall @100 peaks, ENERGY | recall @100 peaks, CC |
|---|---|---|
| PGC | 0.001 | **0.119** |
| B011 | 0.003 | 0.120 |
| B010 | 0.006 | 0.118 |
| B009 | 0.007 | 0.110 |
| B926 | 0.003 | 0.105 |
| VGZ | **0.075** | 0.049 |
| SNB | **0.032** | 0.034 |
| NLLB | 0.011 | 0.043 |
| B005 | **0.022** | 0.005 |
| B001 / B003 / B004 / B006 / B007 / LRIV | 0.001-0.007 | 0.005-0.010 |

Merlin's P1b bar is **single-station S recall ≥ 0.30 at ≤ 100 peaks/station-day at PGC**. The template
cc, the strongest single-station statistic we have there, reaches **0.119**; the energy excess reaches
0.001. A learned model must therefore beat 0.12 at PGC to be worth anything, and 0.30 is a 2.5× jump
over the best incoherent-plus-template statistic on disk. That bar is not obviously reachable, and P1b
should be read as a test of it rather than an expectation.

Files: `results/o3_oracle_decision.json`, `o3_station_roc.csv`, `o3_network_curve.csv`,
`o3_rate_axis.csv`, `o3_template_curve.csv`, `results/oracle/o3_*_oracle2012.*`;
figures `figures/o3_network_roc.png`, `figures/o3_station_roc.png`.

---

## O4. Southern truth-existence test — DEFERRED, nothing queued

`results/o4_deferred.md`. Confirmed by inspection: none of B045, B935, RVIT, PETY exists in any trace
tree on disk. Shelly's catalog (27 families, 61,441 detections, 2018-2024) IS on disk. Until the test
passes, no southern recall number may be quoted.

---

## P0. Label table — `results/labels_v0.parquet` (1,479,426 rows)

| kind | rows | what it is |
|---|---|---|
| `pos_S` | **1,011,174** | one per (Bostock event, station with data and a 2012 template); `t` = absolute S arrival |
| `pos_P` | 433,229 | only where a stack P pick exists (`a5_picks.t_P` finite) |
| `neg_certain_day` | 704 | PNSN-quiet station-days; rule attached: exclude ±120 s around every Bostock event |
| `hard_neg_eq` | 14,818 | catalogued earthquakes / blasts (Track C `eq_catalog.csv`) × stations with data |
| `hard_neg_hfveto` | 16,427 | HF-veto quiet-day transients of the frozen detector |
| `hard_neg_text` | 1,399 | T_ext quiet-day detections of the frozen detector |
| `unlabelled_ets_day` | 1,675 | ETS station-days — non-event windows here are UNLABELLED, never negatives (Merlin F6) |

**S time rule.** `DB time + frozen 2012 family offset + [ a5 stack S pick, else Crosson ttS ] + a5
per-station correction`; σ = 0.5 s on every row; `provenance` ∈ {picked, predicted}. (The a5 table
stores `t_S` with the station correction removed for locating, so the correction is added back here to
give the expected *observed* arrival.) P rows carry the stack pick only.

| by provenance | rows |
|---|---|
| picked | 697,325 (69 %) |
| predicted | 313,849 (31 %) |

Per station (picked / predicted) — the asymmetry Merlin flagged is confirmed and is now a column, not a
footnote:

| sta | picked | predicted | | sta | picked | predicted |
|---|---|---|---|---|---|---|
| B009 | 101,143 | 8,464 | | B001 | 9,779 | 35,442 |
| B010 | 101,680 | 9,746 | | B003 | 11,353 | 43,473 |
| B011 | 95,675 | 13,054 | | **B004** | **0** | **38,925** |
| B926 | 83,949 | 17,508 | | B005 | 11,923 | 26,113 |
| PGC | 69,540 | 4,866 | | B006 | 18,851 | 32,047 |
| SNB | 56,689 | 14,190 | | B007 | 21,090 | 27,814 |
| VGZ | 55,030 | 10,949 | | LRIV | 11,904 | 15,352 |
| NLLB | 48,719 | 15,906 | | | | |

B004 has **no** stack S pick for any family: every one of its 38,925 positive labels is a predicted time.

| year | stations | Bostock events | pos_S | quiet days | ETS days | PNSN coverage |
|---|---|---|---|---|---|---|
| 2008 | 10 (PB only) | 14,224 | 91,451 | — | 22 | **no** |
| 2009 | 10 (PB only) | 22,112 | 148,667 | — | 24 | **no** |
| 2010 | 14 | 18,114 | 137,376 | 1 | 17 | yes |
| 2011 | 18 | 18,663 | 186,926 | 8 | 16 | yes |
| 2012 | 15 | 28,149 | 297,834 | 19 | 20 | yes |
| 2013 | 18 | 14,873 | 148,920 | 16 | 22 | yes |

**A trap caught while building this.** The PNSN tremor catalog begins 2010-01-22 in the family box. A
"no tremor detection that day" rule therefore declares *every* 2008 and 2009 day quiet — including the
ETS days themselves — which would have poisoned the certain-negative class with ~800 station-days of
pure ETS (Merlin's F6, by a different route). The builder now refuses to declare any certain negative
in a period the tremor catalog does not cover, and additionally requires zero Bostock events that day.
2008/2009 contribute positives and unlabelled ETS days only.

**Holdouts are flagged, not dropped**: 421,496 rows carry `hidden_family = True` and 250,431 carry
`held_out_station = True`. Every P1-P3 training set must filter on those two columns.

Summary json: `results/labels_v0_summary.json`.

---

## What this means for the phase

1. The referee is sealed and the time base is frozen; 2013 has been looked at once and is now closed.
2. The frozen template detector reproduces itself on two unseen years — the referee setup is sound, and
   P1a is **0.162 / 0.200** (2011, T_core / T_ext) at 6 and 21 network-wide false events per day.
3. **P2 is discovery-only.** No sensitivity or completeness claim may be attached to the learned model
   on the north permanent network; the production detector for any accepted family stays the frozen
   template machinery. The learned model's job is to propose sources the templates cannot see.
4. The one place the oracle says a learned per-station model has genuine headroom is the **far
   boreholes**, where the template cc is at chance (AUC 0.52-0.56) and a source-agnostic transient
   statistic already reaches 0.64 — and that is also where the labels are predicted rather than picked,
   so P1b must be scored separately by label provenance.
5. 2008/2009 breadth is borehole-only; CN pre-2010 is not obtainable from any reachable FDSN endpoint.
