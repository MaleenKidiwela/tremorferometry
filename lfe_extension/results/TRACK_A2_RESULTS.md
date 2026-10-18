# TRACK A-2 RESULTS — Merlin steps 4 (epoch transfer / identity) and 5 (relocation with controls)

Agent: Track A-2 engineer, 2026-09-16.  Everything new is prefixed `a4_` (epochs) or `a5_` (relocation).
`../lfe_catalog_repair` and `../codex-LFE` were read ONLY; every script used is a copy in `scripts/`.
Nothing belonging to the other two agents (`a0_`–`a3_`, `e5r_*`, `e6_*`, `c_*`/`trackC`) was written or
edited.  Frozen configuration throughout (`results/a1_frozen_a1_mergeB_L12.json`): mergeB site merge,
12-s templates, 2012 templates / weights / moveout, T_core = 8.75 MAD, T_ext = 8.25, 20 Hz / 1–8 Hz
measurement chain, two-site coincidence gate (2 sites each ≥ 2.5 MAD of their own site series).

---

## STEP 4 — epoch transfer and identity

### Protocol actually run (Merlin's, plus two controls he did not specify)

For each (station, boundary) the reference events are the detections of a matched-filter pass with
**that station EXCLUDED from the network sum**, at T_core = 8.75 on the tremor days of the first cached
SSE of the new epoch; the station's own 3-C data is then cut at the predicted S arrival
(t_detection + ttS, Crosson 1-D to the Bostock hypocentre), rms-normalised, cc-aligned (one shift per
station, all three components — `e5r_templates.cut_windows` / `align_stack`, unchanged) and stacked;
the new stack is compared with the 2012 from-scratch stack of the same (family, station) over
S−1..S+7 s by a signed 3-C cc maximised over ±2 s, whose lag is the S-pick shift dS.

Two controls were added because the pre-registered bars are absolute and the two sides are not
otherwise measured identically:

* **matched-N previous-epoch control leg** — the same withheld-station selection and stack run on the
  2012 ETS tremor days (doy 247–270), which lie in the PREVIOUS epoch for all five boundaries, with each
  family capped at the number of events the new leg used (cc between two stacks rises with N on either
  side, and the 2018/2021 SSE windows give fewer detections per family than the 2012 ETS).  This is the
  *floor* of the procedure: what cc, dS and registration come out when the instrument has NOT changed.
* **B926 negative control** — B926 is PB EH 100 Hz throughout with no epoch boundary, so the identical
  procedure across the same 2012 → 2018 gap measures what the YEAR costs by itself.

`registration` = median LEAVE-ONE-OUT per-event cc of the individual event windows against their own
stack over S−1..S+7 (LOO removes the trivial N-dependence); the pre-registered bar is
registration ≥ 0.8 × the previous-epoch (control-leg) value for the same family and station.

### Windows and networks

| boundary | station | new-epoch window | network (target excluded) |
|---|---|---|---|
| 2017-08-10 BH→HH | PGC | 2018 doy 171–196 (SSE 2018-06-20..07-15) | 14 stations; **B009/B010/B011 kept** — the natural control Merlin specifies (same site, unchanged instruments) |
| 2017-03-03 sensor swap | VGZ | 2018 doy 171–196 | 14 stations |
| 2017-02-23 BH→HH | NLLB | 2018 doy 171–196 | 14 stations |
| 2017-06-08 sensor swap | SNB | 2018 doy 171–196 | 14 stations |
| 2019-05-30 BH→HH | LRIV | 2021 doy 276–302 (SSE 2021-10-03..10-29) | 13 stations, B005/B006/B007 kept (same site) |
| 2017-03-03 (secondary) | VGZ17 | 2017 doy 77–93 (SSE 2017-03-18..04-03) | only PGC, SNB, LRIV, B011 are on disk for 2017 |
| none (negative control) | B926 | 2018 doy 171–196 | 14 stations |

NLLB's boundary (2017-02-23) precedes the cached 2017 SSE but NLLB has **no 2017 data on disk at all**,
so it was tested in the 2018-06 SSE; SNB's boundary (2017-06-08) post-dates the 2017 SSE, likewise 2018.

### StationXML audit (`results/a4_stationxml_audit.csv`)

Fetched independently from EarthScope, level=response, all epochs.  **Azimuths and dips are identical
across every boundary** (0/90 horizontals, −90 vertical); what changes is the sensor and the gain:

| station | boundary | sensor before → after | sensitivity (counts / m s⁻¹) |
|---|---|---|---|
| PGC | 2017-08-10 | Guralp CMG3T_360sec → Nanometrics Trillium 120 PH | 5.091e9 → 4.796e8 (**÷10.6**) |
| VGZ | 2017-03-03 | Guralp CMG3ESP_60sec → Nanometrics Trillium 120 PH | 6.788e9 → 4.796e8 (**÷14.2**) |
| NLLB | 2017-02-23 | Guralp CMG3ESP_60sec → Nanometrics Trillium 120 PH | 8.486e8 → 4.796e8 (÷1.8) |
| SNB | 2017-06-08 | Guralp CMG3ESP_60sec → Nanometrics Trillium 120 PH | 8.485e8 → 4.796e8 (÷1.8) |
| LRIV | 2019-05-30 | Velocity transducer, BH → HH | 7.556e8 → 5.037e8 (÷1.5) |

This is one CNSN fleet upgrade: PGC, VGZ, NLLB and SNB all end on the *same* Trillium 120 PH response.
Our chain does not remove the response, but cc is amplitude-normalised, so a gain change alone cannot
move cc; and PGC (÷10.6) transfers essentially perfectly, which rules the gain out as the mechanism for
any station that does not.

### Results (`results/a4_epoch_identity.csv`, `results/a4_epoch_summary.csv`, `figures/a4_epoch_identity_<STA>.png`)

n = families stacked on the new leg; "ctl-OK" = families whose matched-N 2012 control leg itself passes
|cc| ≥ 0.6 and |dS| ≤ 0.5 s (i.e. families where the test can say anything at all).

| station | boundary | n | ctl-OK | ACCEPT | ACCEPT / ctl-OK | cc median new → ctl | cc ratio | dS median (MAD) | registration ratio | polarity flips |
|---|---|---|---|---|---|---|---|---|---|---|
| **PGC** | 2017-08-10 BH→HH | 66 | 61 | 46 | **0.74** | 0.805 ← 0.893 | 0.90 | **+0.00 s (0.00)** | **0.94** | 0 |
| **SNB** | 2017-06-08 sensor | 74 | 56 | 37 | **0.63** | 0.622 ← 0.737 | 0.84 | +0.00 s (0.05) | 0.89 | 3 |
| **NLLB** | 2017-02-23 BH→HH | 56 | 36 | 22 | **0.50** | 0.619 ← 0.709 | 0.87 | +0.00 s (0.05) | 0.88 | 1 |
| **VGZ** | 2017-03-03 sensor | 65 | 43 | 7 | 0.14 | 0.478 ← 0.662 | 0.72 | +0.00 s (0.05) | 0.71 | 0 |
| **LRIV** | 2019-05-30 BH→HH | 13 | 5 | 2 | 0.40 | **−0.216** ← 0.542 | — | +0.05 s (0.35) | 0.75 | 0 |
| VGZ (2017 SSE, 4-station net) | 2017-03-03 | 8 | 7 | 3 | 0.43 | 0.522 ← 0.696 | 0.75 | +0.00 s (0.00) | 0.60 | 0 |
| **B926 — NEGATIVE CONTROL, no boundary** | 2012→2018 only | 56 | 20 | 3 | **0.05** | 0.404 ← 0.510 | **0.79** | +0.00 s (0.05) | **0.66** | 1 |

**The negative control is the headline.** B926 has no instrument change of any kind, yet across the same
2012 → 2018 gap and the same procedure its cc median falls 0.510 → 0.404 (ratio 0.79), its registration
ratio is 0.66, and only **5 %** of its control-OK families satisfy the pre-registered bars.  So the
absolute bars |cc| ≥ 0.6 and registration ≥ 0.8 × previous epoch are failed by a station whose
instrument did not change: as written they measure the change in the EVENT POPULATION and the
detector's per-event purity between a 2012 ETS and a 2018 SSE, not instrument identity.  Read against
that control:

* **PGC (0.90 / 0.94), NLLB (0.87 / 0.88) and SNB (0.84 / 0.89) all transfer BETTER than the
  no-change control.** There is no detectable instrument-induced identity loss at those three
  boundaries.  PGC is the strongest case and it is also the one with the natural control Merlin asked
  for (B009/B010/B011 unchanged in the network): its six best families overlay the 2012 BH templates at
  cc 0.88–0.93 (`figures/a4_epoch_identity_PGC.png`).
* **VGZ (0.72 / 0.71) is statistically indistinguishable from the no-change control (0.79 / 0.66).**
  Its 86 % "failure" rate against the absolute bars is therefore NOT evidence of a metadata problem.
* **LRIV is the one genuine failure.** Its cc median goes NEGATIVE (−0.22) while the control leg on the
  same 13 families is +0.54, and its dS MAD (0.35 s) is 7× every other station's.  Only 13 of 130
  families have a usable LRIV template at all, and LRIV is the lowest-registration station tested.
* **dS is the cleanest result of step 4: the median S-pick shift is +0.00 s at EVERY boundary**, with a
  MAD of 0.00–0.05 s (LRIV 0.35 s).  No timing shift was introduced by any of the five instrument
  changes.  A gain change of ÷10.6 (PGC) or ÷14.2 (VGZ) leaves the waveform shape and timing intact.
* Polarity: 5 flips over 338 family-stations, **one of them at the no-change control B926** → these are
  noise at the cc ≈ −0.6 boundary, not real polarity reversals.  No station shows a systematic flip.

**Chance level of the |cc| ≥ 0.6 bar (extra control, Merlin step-3 criterion (d)).** Correlating each new
stack against a DIFFERENT family's 2012 stack at the same station: |cc| p95 = **0.67 at PGC**, 0.52 NLLB,
0.50 SNB, 0.42 VGZ, 0.42 LRIV, 0.32 B926; the fraction of cross-family pairs reaching 0.6 is 10.5 % at
PGC and ≤ 3 % elsewhere.  **At PGC the 0.6 bar has a ~10 % false-accept rate from neighbouring families
alone** — neighbouring families genuinely look alike at a single station.  Epoch identity should be
judged against this per-station chance level, not against a universal 0.6.

**StationXML audit** was triggered by the > 50 % rule at VGZ, LRIV, VGZ17 — *and at the no-change
negative control B926*, which is itself the proof that the rule as written fires on the year effect.
The audit found **no metadata fault at any boundary**: azimuths 0/90, dips 0/−90, unchanged everywhere
(table above).  What changed is one CNSN sensor upgrade (Guralp → Nanometrics Trillium 120 PH) with a
gain drop of ÷1.8 to ÷14.2, and PGC — the largest gain drop — transfers best.

### What Track C should use

1. **Per-epoch template sets are written and ready**: `results/a4_templates_<STA>_<epoch>.npz`, in the
   `--extra-tpl` format (`tpl` {(family, station): {12.0: (3, 240) L2-normalised}}, `sp` placements,
   plus the full ±20 s stacks), polarity already applied: **PGC/2018_HH 46 templates, SNB/2018_HH 37,
   NLLB/2018_HH 22, VGZ/2018_HH 7, LRIV/2021_HH 2, B926/2018_EH 3.**
2. **Do not re-stack PGC, NLLB or SNB just because of the 2017 boundary** — the 2012 templates transfer
   across it better than a no-change station transfers across the same years.  Re-stacking per epoch is
   still worth doing for SENSITIVITY (the new-epoch stacks are built from the new epoch's own events),
   but it is not required for identity.
3. **LRIV must not carry 2012 templates past 2019-05-30.** Either re-stack it per epoch from
   withheld-station detections (only ~13 families have enough events) or drop LRIV from the network for
   the post-2019 epochs and absorb the loss in the completeness model.  LRIV is in any case one of the
   two lowest-yield stations for relocation as well (S-pick acceptance 0.41).
4. **Re-derive the epoch bars empirically.** The pre-registered |cc| ≥ 0.6 / registration ≥ 0.8 bars
   should be replaced by *relative* bars: cc ratio and registration ratio against a same-gap no-change
   station (B926 here; B001/B003/B004/B009/B010/B011 are equally eligible, they never change), and by
   the per-station cross-family chance level instead of the universal 0.6.
5. **Two operational facts for the production loop**: snap CN sampling rates to the nearest integer
   before merging (41 of 504 station-days in 2021 are otherwise lost), and run at most 2–3 concurrent
   `e4_template_test.py` passes on the shared H100.


---

## STEP 5 — relocation with controls

### (a) Synthetic floor FIRST — the 5 km bar IS testable on this geometry

`scripts/a5_synth_floor.py`, machinery in `scripts/a5_locate.py` (a memory-lean copy of
`d3_stage4_locate.locate_from_tables`, **verified identical to the original on 12 random
family × prior cases**, so the estimator is D3's, not a new one).  Nine SITES — PGC cluster = ONE
station, LRIV cluster = ONE station (Merlin F10(i)) — Bostock's 130 hypocentres forward-modelled through
Crosson 1-D, 0.3 s S / 0.2 s P Gaussian pick noise plus a per-site static drawn U(−0.5, +0.5) s
(applied to S and, scaled by 1/1.78, to P), 50 realisations per family, grid 124 × 124 km at 1 km and
10–60 km depth at 0.5 km.

| configuration | median Δh | p90 Δh | median \|Δz\| | at bounds |
|---|---|---|---|---|
| noise-free (grid + table resolution) | 0.41 km | 0.56 | 0.27 km | 0 % |
| **S+P, free depth, with noise** | **1.43 km** | 2.67 | 1.70 km | 0 % |
| S+P, Slab2 prior | 1.44 km | 2.76 | 1.92 km | 0 % |
| S-only, free depth | 1.52 km | 3.00 | 2.16 km | 0 % |
| matched to the pick sets actually obtained, all 88 located families | 2.73 km | 6.37 | 2.27 km | 0 % |
| matched, restricted to the 14 ACCEPT-geometry families | **2.03 km** | 4.17 | 2.29 km | 0 % |

**The synthetic median horizontal error is 1.4–2.0 km, well below 5 km → the bar is testable.**
(E3's "S-only solutions collapse to the network centroid" pathology was a 3-station problem; with 9
sites and azimuthal gap ≤ 180° the geometry is not the limitation.)

### (b) Station corrections from intraslab earthquakes, not from Bostock

USGS, M ≥ 2, depth 25–50 km, 48.0–49.0 N / −124.5..−123.0 E, 2010–2026: 42 events, of which **7 have
≥ 14 of our stations on disk** (2011-08-06, 2011-09-06, 2012-08-29, 2014-08-01, 2014-12-10, 2015-02-07,
2015-12-17; M 2.08–3.30, depth 26–47 km).  Picked on our own traces at native rate, 2–15 Hz,
envelope-onset + AIC, accepted at peak/noise-rms ≥ 4 → **97 P and 98 S picks**.

Median residual vs Crosson (s), = the station correction applied in (c):

| | PGC | VGZ | SNB | NLLB | LRIV | B926 | B001 | B003 | B004 | B005 | B006 | B007 | B009 | B010 | B011 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S | +0.33 | −0.05 | +0.67 | +0.57 | +0.57 | +0.01 | +0.31 | +1.13 | +0.17 | +0.48 | +0.49 | +0.46 | +0.11 | +0.11 | −0.05 |
| P | +0.10 | +0.05 | +0.94 | +0.43 | +0.14 | +0.30 | +0.17 | +0.87 | +1.03 | +0.09 | +0.09 | +0.09 | +0.27 | +0.22 | +0.13 |
| n(S) | 7 | 7 | 4 | 6 | 7 | 5 | 7 | 7 | 6 | 7 | 7 | 7 | 7 | 7 | 7 |

**These are 2–3× SMALLER than the +1.7 to +3.25 s that E3 reported for the same stations from Bostock's
hypocentres** — i.e. most of that apparent "station residual" was family mislocation in the reference,
exactly the circularity Merlin's F10(ii) warned about.  CAVEAT (`figures/a5_eq_picks.png`): the **P**
residual grows systematically with epicentral distance (≈0 s below 40 km, +0.5..+2.5 s beyond 60 km) —
a Crosson-model bias that no constant per-station term can absorb.  The **S** residuals show no such
trend.  This is why the S-only sensitivity run below matters.

### (c) Relocation

P on Z and S on the horizontals of OUR own 2012 from-scratch per-(family, station) stacks, which are
already ±20 s long inside `e5r_templates_2012.npz` (P sits 4–8 s before S, well inside the window, so no
longer stacks had to be rebuilt).  Envelope-onset + AIC, noise window ending 3 s before the predicted P,
SNR ≥ 3.  Site collapse = median over the cluster members of (member pick − that member's earthquake
correction).  Run A free depth, Run B Slab2 prior (σ = 5 km); cost function, σ_S = 0.30 s, σ_S−P = 0.50 s
unchanged from D3.

Pick yield: 1,144 stacks → 768 S and 481 P picks → 446 site-S picks → **88 families with ≥ 3 picked
sites; 14 with ≥ 4 picked sites within 60 km AND azimuthal gap ≤ 180°** — the pre-registered ACCEPT
subset.  The binding constraint is the pick yield at the far sites, not the network: S-pick acceptance
is 0.86–0.92 at the PGC cluster, 0.72–0.81 at NLLB/SNB/VGZ, 0.77 at B926, but 0.17–0.41 at
B001/B003/B005–7/LRIV and **0.00 at B004**.

| run | n | median Δh | p90 | median Δz | at bounds | shuffled-null Δh | obs/null | spread strike / dip |
|---|---|---|---|---|---|---|---|---|
| **A free depth (primary)** | 14 | **3.61 km** | 5.76 | −2.20 (\|Δz\| 2.93) | 0 % | 29.4 km | **0.12** | 0.95 / 0.97 |
| **B Slab2 prior** | 14 | **3.45 km** | 5.76 | −3.36 | 0 % | 28.9 km | 0.12 | 0.95 / 0.98 |
| A, no station corrections | 14 | 4.88 km | 6.82 | −4.61 | 0 % | 27.9 km | 0.17 | 0.84 / 0.93 |
| A, S-only (no S−P) | 14 | 2.88 km | — | −2.02 | 0 % | 29.4 km | 0.10 | 0.97 / 0.99 |
| (matched synthetic floor) | 14 | 2.03 km | 4.17 | — | 0 % | — | — | — |
| all 88 located families (reported, not barred) | 88 | 5.26 km | — | −3.91 | 0 % | — | — | — |

**Pre-registered bars, ACCEPT subset (14 families):**

| bar | value | verdict |
|---|---|---|
| horizontal median ≤ 5 km to OWN Bostock family | 3.61 km (A) / 3.45 km (B) | **MET** |
| observed median < 0.5 × shuffled-pick null (100 realisations) | 0.12 × | **MET** |
| spread ≥ 0.7 × Bostock's along strike and dip | 0.95 / 0.97 | **MET** |
| ≤ 5 % of depths at the grid bounds | 0 % | **MET** |
| depth — reported, no bar | median Δz −2.20 km (free) / −3.36 km (Slab2), median \|Δz\| 2.93 km | reported |

10 of the 14 accepted families are individually within 5 km; 42 of all 88 located families are.
The observed 3.6 km is 1.8× the matched synthetic floor, so most of the residual is pick noise, and the
result is robust to both the station corrections (4.88 km without them, still under the bar) and to
dropping P altogether (2.88 km, i.e. *better* — consistent with the distance-dependent P bias in (b)).

Pick-gate sensitivity (pre-declared as a sensitivity, the primary stays SNR ≥ 3): at SNR ≥ 2.5 the
yield rises to 97 located families and **21** in the accept subset, with median Δh 4.02 km (Run A) /
3.88 km (Run B), obs/null 0.13, spread 1.04 / 0.99, 0 % at bounds — every bar still met, 50 % more
families.  That is the cheapest lever if Track C needs more located families.

**The remaining 74 families are reported separately** in `results/a5_relocation.csv`
(`verdict = reported_only(geometry)`): they have fewer than 4 picked sites within 60 km or a gap > 180°,
almost always because the far boreholes give no S pick on their stacks.

---

## Deliverables

| file | contents |
|---|---|
| `results/a4_epoch_identity.csv` | per family-station: station, boundary, epoch, family, n_events, cc, dS, reg_ratio, verdict, polarity_flip, control-leg values, stack S-SNR |
| `results/a4_epoch_summary.csv` | one row per (station, boundary): counts, medians, ratios, cross-family chance level, `audit_stationxml` flag |
| `results/a4_templates_<STA>_<epoch>.npz` | accepted per-epoch templates, `--extra-tpl` format, polarity applied (PGC 46, SNB 37, NLLB 22, VGZ 7, LRIV 2, B926 3) |
| `results/a4_stacks_<STA>_<epoch>_<leg>.npz` | every new/control stack, including the rejected families |
| `results/a4_stationxml_audit.csv` | response epochs, azimuths, dips, sensitivities for the six stations |
| `figures/a4_epoch_identity_<STA>.png` | 6-family before/after overlays + cc, dS and registration panels |
| `results/a5_relocation.csv` | per family: Run A and Run B lat/lon/z, Bostock lat/lon/z, Δh, Δz, n_sites, gap, verdict |
| `results/a5_controls.json` (+ `_nocorr`, `_sonly`, `_snr25`) | all four controls and the pre-registered verdicts |
| `results/a5_synth_floor.{csv,json}` (+ `_matched`) | synthetic floor, per realisation and summarised |
| `results/a5_station_corrections.csv`, `a5_eq_picks.csv`, `a5_eq_catalog_raw.csv` | earthquake station terms and the picks behind them |
| `results/a5_picks.csv`, `a5_picks_site.csv`, `a5_geometry.csv` | family-stack picks, site-collapsed picks, per-family geometry |
| `figures/a5_synthetic_floor.png`, `a5_relocation_map.png`, `a5_eq_picks.png` (+ variants) | the step-5 figures |

Scripts (all in `scripts/`, all copies — nothing outside `lfe_extension/` was written):
`a4_cache20.py`, `a4_cache_windows.sh`, `a4_run_selection.sh`, `a4_epoch_stack.py`, `a4_run_stacks.sh`,
`a4_run_matched.sh`, `a4_finish2.sh`, `a4_epoch_report.py`, `a5_locate.py`, `a5_synth_floor.py`,
`a5_eq_corrections.py`, `a5_relocate.py`.

## Caveats and known limits

1. The step-4 absolute bars are not usable as written — see the B926 negative control.  Every verdict
   above is stated relative to it.
2. The VGZ 2017 re-test has only 8 families because just 5 stations have 2017 data on disk; when the
   downloader reaches 2017 this should be redone with the full network.
3. Step-4 event selection did not apply the frozen HF-ratio veto: its cut (3.6384) was calibrated on
   PGC BH 10–19 Hz in 2012 and PGC is HH 100 Hz after 2017-08-10, so transferring the cut across the
   band change is itself an untested epoch question.  On SSE tremor days at T_core the veto mainly
   removes quiet-day glitches, so its absence should bias the reference sets only slightly.
4. The relocation ACCEPT subset is 14 of 130 families.  The limit is the S-pick yield on the far
   borehole stacks (0.00 at B004), not the network geometry — the synthetic floor says the geometry
   would deliver 1.4 km if every site were picked.
5. The earthquake P residuals carry a distance-dependent Crosson-model bias that a constant station
   term cannot absorb; the S-only run is given for exactly that reason and is the more trustworthy of
   the two (2.88 km).
6. Bostock's family depths are the reference for Δz, and they are themselves model-dependent; the
   −2 to −3 km offset should not be read as an absolute depth correction.

---

## ADDENDUM (Track C hand-off) — per-epoch HF-veto cuts, and the epoch-boundary list

### 1. The HF-ratio cut: the shift is a POPULATION effect, not an instrument effect

Three independent diagnostics, in the order they kill the candidate explanations:

**(a) Instrument response (`a4_hf_response.py`, `results/a4_hf_response_scan.csv`).** Evaluating
|Resp(f)|² from StationXML for every CN.PGC epoch, the band-averaged response predicts the measured
ratio to shift by only **×1.16** across 2017-08-10 in the 10–19 Hz band (×1.06 in 10–16).  Track C
measured **×2.41**.  So the BH anti-alias roll-off accounts for at most a fifth of the shift.

**(b) Band redefinition does NOT fix it** (`a4_hf_bands.py`, `results/a4_hf_band_scan.csv`; the HF ratio
recomputed from the native traces on all 27 days of 2015-12 and all 21 of 2018-06 — it reproduces Track
C's medians, 0.199 and 0.480, exactly).  HH/BH shift of the median by band: **10–19 ×2.41, 10–16 ×1.86,
10–14 ×2.06, 8–16 ×1.59, 12–16 ×1.86.**  No band is epoch-invariant, so the coordinator's candidate 2 is
ruled out.

**(c) The distribution moves with the MONTH, not the sensor.** With the *same* BH instrument the median
of all detections is 1.019 on 2012 doy 228–279 and 0.199 on 2015-12 — a 5× spread with no boundary
crossed.  Restricted to a matched LFE-like core (net_mad ≥ T_core, n_sites ≥ 4, PNSN tremor days) the
medians are **0.256 (2012 BH) / 0.169 (2015-12 BH) / 0.247 (2018-06 HH)**: the BH-to-BH spread is as
large as any BH-to-HH difference.  What changed in 2018-06 is the *fraction of contaminants* in the
detection list, which is exactly what the veto is meant to remove — the frozen cut vetoes 12.0 % of all
2018-06 detections but only **0.5 % of the LFE-like core**.

### The rule, and why not the one proposed

| rule | cut on 2012 | cut on 2015-12 (same instrument) | stable? |
|---|---|---|---|
| referee p99 (own-matched) — the frozen rule | 3.6384 | (no referee) | — |
| **RULE-M** = F × median(core), F = 14.22 fixed on 2012 | 3.6384 | 2.40 [2.36, 2.47] | yes (±35 % SSE-to-SSE) |
| RULE-Q = quantile q* = 99.196 % of the core | 3.5977 | 1.31 | no |
| **RULE-T** = p99 of the core (the proposed rule) | 3.2456 | **1.2077** | **NO — 2.7× between two BH months** |

RULE-T is a tail statistic of a mixture whose mixing fraction changes month to month, so it is not
stable even within one instrument epoch; it is reported in `a4_hf_veto_cuts.csv` as `diagnostic` and
**not recommended**.  RULE-M keys on the *median of the LFE-like core*, which the diagnostics above show
to be the epoch-stable quantity, and it reproduces the referee cut on 2012 by construction.
Validated on 2012 (`results/a4_hf_rule_validation.csv`), RULE-M, RULE-Q and the frozen cut all meet the
original adoption bars identically: quiet-day FA reduction **0.59** (T_ext) / 0.57 (T_core) at a
relative test-half recall cost of **0.003** (bars: ≥ 0.30 and ≤ 0.02).

### Deliverable — `results/a4_hf_veto_cuts.csv`

| station_epoch | band | fs | cut | status | method | n_calib |
|---|---|---|---|---|---|---|
| PGC_BH 2010-08-01..2017-08-10 | 10–19 | 40 | **3.6384** | recommended | referee p99 (own-matched, 2012 build half) = the frozen cut | 6549 |
| PGC_HH 2017-08-10..2026 | 10–19 | 100 | **3.5129** | recommended | RULE-M: F = 14.22 × core median of 2018-06 | 1104 |
| (+6 diagnostic rows: RULE-Q and RULE-T per leg) | | | | diagnostic | | |

**The practical conclusion is reassuring: the HH cut the rule asks for, 3.51, has a 95 % bootstrap CI of
[3.34, 3.83], which CONTAINS the frozen 3.6384.** The frozen cut was never badly wrong after
2017-08-10 — the 6× rise in the vetoed fraction is the veto correctly removing a contaminant population
that is genuinely larger in 2018-06, not a mis-calibration.  Track C can adopt 3.5129 for the HH epoch
(or keep 3.6384 and lose nothing measurable); what it should NOT do is treat the post-2017 `hf_pass`
column as suspect.  Ratio distributions and veto fractions per epoch: `results/a4_hf_epoch_profile.csv`,
`figures/a4_hf_epoch_cuts.png`.

**Caveat / follow-up**: RULE-M is calibrated per epoch on that epoch's first cached SSE; the SSE-to-SSE
spread inside one instrument epoch is ±35 %, larger than any instrument step measured here.  If Track C
wants tighter control the right granularity is **per SSE**, not per epoch — `a4_hf_epoch_cuts.py` will
produce that from any month's detections parquet.  Also note that no separate cut is needed for the
PGC BH sub-boundary 2015-01-12 (CMG3T_NSN → CMG3T_360sec): its predicted response effect is 4.6 % and it
is far inside the SSE-to-SSE spread.

### 2. Epoch-boundary list — CONFIRMED, all five are real response changes

From the independent StationXML pull (`results/a4_stationxml_audit.csv`); azimuths (0/90) and dips
(0/−90) are unchanged at every one of them, so none is an orientation problem:

| boundary | band change? | sensor change | gain change | own template epoch? |
|---|---|---|---|---|
| PGC 2017-08-10 | BH 40 → HH 100 | Guralp CMG3T_360sec → Nanometrics Trillium 120 PH | ÷10.6 | **yes** |
| NLLB 2017-02-23 | BH 40 → HH 100 | Guralp CMG3ESP_60sec → Trillium 120 PH | ÷1.8 | **yes** |
| SNB 2017-06-08 | **no** (HH→HH) | Guralp CMG3ESP_60sec/GD2_v21 → Trillium 120 PH | ÷1.8 | **yes** |
| VGZ 2017-03-03 | **no** (HH→HH) | Guralp CMG3ESP_60sec/GD2_v21 → Trillium 120 PH | ÷14.2 | **yes** |
| LRIV 2019-05-30 | BH 40 → HH 100 | velocity transducer, re-rated | ÷1.5 | **yes** |

SNB 2017-06-08 and VGZ 2017-03-03 are absent from `results/station_epochs.csv` only because that table
keys on the seismometer BAND; both are genuine response epochs in StationXML (different sensor, different
overall sensitivity) and both must be treated as template epochs.  The whole set is one CNSN upgrade —
PGC, VGZ, NLLB and SNB all end on the *same* Trillium 120 PH response (4.79595e8 counts m⁻¹ s).

### 3. Per-epoch template sets — proposed, and the coordinator's decision

| sta | epoch_start | npz | templates = frozen-2012 + new-epoch |
|---|---|---|---|
| PGC | 2017-08-10 | `results/a4_templates_merged_PGC_2017-08-10.npz` | 111 = 65 + **46** |
| NLLB | 2017-02-23 | `results/a4_templates_merged_NLLB_2017-02-23.npz` | 88 = 66 + **22** |
| VGZ | 2017-03-03 | `results/a4_templates_merged_VGZ_2017-03-03.npz` | 97 = 90 + **7** |
| SNB | 2017-06-08 | `results/a4_templates_merged_SNB_2017-06-08.npz` | 101 = 64 + **37** |

These are MERGED sets (`a4_make_epoch_sets.py`): every 2012 template of that station, replaced family by
family by the new-epoch re-stack wherever the epoch-identity test ACCEPTED it (polarity applied; the 8-s
template of a replaced entry is the re-normalised first 8 s of its own 12-s template).  Handing over the
identity-filtered `a4_templates_<STA>_<epoch>.npz` alone would have deleted two thirds of PGC's family
coverage; the merged set strictly dominates the status quo — full coverage, epoch-matched waveforms
wherever the identity evidence supports them.

**LRIV was deliberately left out of the proposal** (the file `a4_templates_merged_LRIV_2019-05-30.npz`
exists, and `a4_make_epoch_sets.py --with-lriv` adds the row). Its epoch identity FAILED — cc median
−0.22 against its own 2012 stacks, worse than the no-change control — so relabelling its post-2019 rows
`per_epoch:<file>` would have asserted a transfer that was not verified.  **B926 is also not listed**:
it has no epoch boundary, it was only the negative control, and its 3-family "accepted" set would
cripple a station that needs no per-epoch treatment at all.

**Coordinator decision (2026-09-16), recorded here as the operative state.** The proposal above is
preserved at `results/trackC/templates_by_epoch.a4_proposed.csv`; the operative
`results/trackC/templates_by_epoch.csv` now carries the coordinator's own schema and decision:
**keep the frozen 2012 templates for PGC, NLLB, SNB and VGZ** (label `2012_final` — consistent with
this report's recommendation 2, since all four transfer at or above the no-change control) and
**DROP LRIV from the network after 2019-05-30** (`drop=1`), with the loss recorded per family-day for
the completeness model.  The merged sets remain on disk if a future pass wants the epoch-matched
waveforms for sensitivity rather than identity.
