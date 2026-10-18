# RESULTS — lfe_catalog_repair (2026-09-10; D1 + D3 pilot complete)

## Verdict
The project's LFE families cannot yield an event catalog by matched filtering with the current template design, from stored hits
(D1) or from re-detected waveforms (D3), single-station or network-stacked over 9 stations. This is a measured limit, not a tuning
problem.

## The numbers that decide it
- A median Lin/Bostock LFE at PGC has amplitude s/n = 0.078 of the ambient 1-8 Hz noise (130-event refined-pick stack).
- Single-channel cc noise for a 2-s template: sigma 0.168 (theory 1/sqrt(2BL) = 0.189). The old fixed cut cc >= 0.80 is ~6 sigma
  -> 100-260 false alarms per template per day; the ETS excess is +25 % spread over all cc (D1).
- One LFE registers at 1.8 MAD on a 3-C station and 3.5 MAD on the 9-station / 21-channel network sum -> undetectable at any
  threshold. Network stacking behaves as theory (noise / 3.3 for 21 channels) and sigma falls as 1/sqrt(template length).
- After MAD thresholding the false-alarm load falls 4,000x (258 -> 0.06 det/template/day) but the survivors still show ETS/quiet
  >= 5 for 0-0.5 % of templates, tremor rho ~ 0, and Lin recall/precision equal to a +600 s shifted control.
- Stacking stored detection times at neighbouring stations gives no S arrival above a random-time stack (waveform-level proof
  that the stored times are template-shaped noise).
- Required for detection at k = 6 MAD: ~20 three-component stations with 6-s network-consistent templates; k = 9: ~47 stations.
- Location geometry is NOT the obstacle: with 9 stations and 0.2-s picks, family hypocentres to 0.6-0.9 km horizontal, 0.7-1.7 km
  depth are reachable (synthetic test, Crosson 1976 1-D model). No real family reached 4 usable stations.

## Route — SUPERSEDED 2026-09-10 (see PLAN.md 'CORRECTED PLAN'): the network-wide catalog is built with a deep-learning picker
## (Lin et al. 2024 replication), not template matching. The text below documents the earlier reasoning.
## Route (revised after re-reading lfe_features/ R1-R11, the June 2026 picker work)
Both threads converge: no single-station route yields individual LFEs — the per-sample picker hit the same structural wall (R7:
pick precision 0.03, precise labels changed nothing) that D3 quantifies (1.8 MAD per 3-C station). And R11 showed longer templates
do not help because the LFE waveform is coherent only over ~2 s. Therefore:
1. FLEET-WIDE LFE-ACTIVITY catalog (not events): deploy the existing per-window picker (R4 LFE AUC 0.975; R9 rejects EQ/blast at
   0.97; retrain with TNOISE + autoencoder latents to push the 0.81 in-tremor residual) as a sliding scan over 3-C continuous data ->
   per-station hourly P(LFE) activity with EQ/blast rejection. A PNSN-independent seismic tremor/LFE-activity measure margin-wide =
   the third measure for the tremor-slip scaling paper. Cost: the same streaming waveform pass as any fleet job; compute light.
2. EVENT-level catalog with hypocentres ONLY where the network allows it (Puget / southern Vancouver Island, >= 20 3-C stations):
   network matched filter seeded by Bostock/Lin families (2-s templates, all channels summed with the family moveout), family
   hypocentres with the stage-4 machinery (geometry gives ~1 km), event-level relative locations for >= 4-station detections.
3. Do NOT pursue single-station per-sample picking or longer templates further; use Bostock/Lin catalogs as ground truth and, in
   the north, as the event catalog.
## Flags for the main project
- The densify stacks used for dv/v are built from these detections: the true-LFE fraction per stack is small and the false alarms
  are template-shaped (they correlate with the template by selection). The dv/v pipeline's own validations (causality gate, PGC vs
  B011, mirror null) stand on their own evidence, but the "detections" must not be read as LFE occurrence times.
- Two data QC holes found: 58 % of WISH's and 27 % of YACT's raw re-detections were zero-amplitude degenerate windows; 12.9 % of
  PGC days are component-short and carried 97.7 % of its raw detections. Worth checking in the dv/v stacking inputs.
- Bug found in an existing file: the Lin parquet is datetime64[us]; .astype("int64")/1e9 yields epoch 1.27e6 (also affects
  lin_precise_b011.csv joins).
Files: 30 scripts, ~70 results, figures d3_*; NOTES.md (949 lines, dated per stage).

---

# F1 — fleet inventory & association regions (2026-09-10, metadata only; report `results/f1_fleet_inventory.md`)
Deliverable: `results/f1_fleet_inventory.md`; tables `results/f1_{stations,station_years,cell_density,regions,region_stations,cost_estimate,groundtruth}.csv`;
figure `figures/f1_fleet_inventory.png`; verification `scripts/f1_verify.py` (31/31 cross-table checks PASS). No waveforms downloaded.

## Inventory (f1_stations.csv, f1_station_years.csv)
- **691 stations = 278 project + 413 non-project FDSN** with a 3-C >=20 Hz triplet within 150 km of an ETS source cell.
  Project classes: 181 broadband, 59 short-period, 36 PB borehole short-period, 2 broadband-borehole.
- **57 project stations never have 3-C >=20 Hz** (32 NC + 25 UW, all Z-only) -> unusable for LFEDetect's ZNE input.
- Project 3-C stations per year: **70 (2009) -> 116 (2018) -> 216 (2021) -> 208 (2026)**; **2,167.6 effective station-years**
  (4,603.1 incl. non-project). ZNE 4,425 vs Z12 893 station-years. 77 % of fleet 3-C station-years are 100-Hz HH.

## Geometry (f1_cell_density.csv, f1_regions.csv; 137 cells, 0.25 deg, Slab2 30-45 km, 40-50 N)
- **At a 100 km radius the PLAN's locatable fraction is saturated: >=3 and >=4 stations in 100 % of cells in 100 % of
  region-years, both rosters.** Station count within 100 km discriminates nothing in Cascadia. All reporting therefore moved
  to a **60 km** radius (proxy anchored on D3's s/n = 0.078 at ~40-50 km) plus a broadband-only split.
- Fraction of cells with **>=4 3-C stations <=60 km** (fleet): 2009 = 0.19-0.87 by region -> **1.00 everywhere from 2021**.
  **>=6 <=60 km**, 2012 -> 2024: R1 0.23->0.95, R2 0.76->1.00, R3 0.43->1.00, R4 0.18->1.00, R5 0.07->1.00, R6 0.00->1.00,
  R7 0.07->0.93. First year >=90 % of cells clear the >=6/60 km bar: R2 2017, R6 2019, R1/R3/R4/R5 2020, R7 2021 —
  but 2009/2010 for R2/R3 if the 413 non-project stations are included.
- Median azimuthal gap at 60 km (2019-26): R4 81 < R3 89 < R5 98 < R2 100 < R6 107 < R1 118 < **R7 134 deg** (worst).

## Cost for the fleet picker run (f1_cost_estimate.csv; measured inputs, not guesses)
- bytes/sample measured from real day files via miniSEED record headers: **1.416 B/sample at 40 Hz, 0.932 at 100 Hz**
  (fitted b(fs) = c0 + c1/fs per station-year; b(20 Hz) = 2.22 extrapolated).
- GPU: **52.6 station-days/min per queue = 6.93 GPU-min/station-year** (E2 measured, L40S, 4 concurrent queues x 2 torch
  threads; snapshot 38 runs / 1,962 station-days at 22:31 — E2 is still appending). **1,867 S picks/station-day** measured.
- **Project fleet: 17.6 TB native download, 9.1 TB (STEIM2) / 16.4 TB (float32) stored at 20 Hz, 250 GPU queue-hours =
  ~63 h wall-clock on one L40S, 1.48e9 S picks.** Project+extra: 36.8 TB, 532 queue-hours (~133 h), 3.14e9 picks.
  Association ~123 CPU-h + location ~8 CPU-h (ASSUMED rates). **The download dominates the GPU by ~5x**; ~15/6/3 days at a
  sustained 10/25/50 MB/s (no download rate was measured). `/home/jovyan` has ~75 TB free — the 80 GB streaming constraint
  in the D3 plan was a self-imposed study-folder budget.

## Ground truth (f1_groundtruth.csv)
- **Lin 2024 exists only at 48.07-49.32 N** (R1 675,615 / R2 382,499 events); its own quality cut (residual<0.5 AND N>=4)
  keeps **9.8 %** (104,071). Bostock 269,423 events / **122 distinct hypocentres for 130 published families**, 48.12-48.93 N.
  SeisBench LFE stacks = **129 families, 21 stations, P AND S picks** (the only fine-tuning substrate), all R1+R2.
  Ducellier 34 located families and Shelly 27 families in R7 (Shelly's release has no coordinates; region from the paper).
- **R3/R4/R5/R6 have zero LFE ground truth — 83 of 137 source cells (61 %).**

## Verdict per region (plain language; cites Lin 2024 §4.3 and Muenchmeyer et al. 2024, see the E1 NOTES entry)
- **Event catalog with hypocentres: R2 (S Vancouver Island / JdF) over the whole span; R1 from ~2018-2020.** R2 is the only
  region with both pre-2016 geometry and the ground truth to measure precision.
- **Conditional (fine-tune + validate out-of-region first): R3 Puget/Olympic and R4 N Oregon from 2019** — best azimuthal
  gaps in the margin and broadband-dominated (so the small-aperture-borehole failure Muenchmeyer hit at San Andreas does not
  apply), but **no in-region ground truth** and network-limited before 2019.
- **Activity time series only: R5, R6** (bar crossed only in 2019-20; mid-record completeness jump would poison rate work)
  **and R7** (worst geometry, but the only southern region with public family ground truth -> per Muenchmeyer, template
  matching against the Shelly/Ducellier families beats the DL picker there).
- No margin-wide temporally-uniform event catalog is supportable. Honest fleet product = R1+R2 (+R3/R4 from 2019) events
  **plus** a margin-wide activity series carrying the per-region-year completeness proxy from `f1_regions.csv`.

## Caveats
Metadata only: availability is an **upper bound** (no gap/noise/glitch check; this project measured 12.9 % component-short
PGC days, 58 % WISH / 27 % YACT zero-amplitude windows). `frac_ge*` is a **geometry proxy, not detectability**. Instrument
class is per-station, not per-station-year; azimuthal gap ignores depth; the 7 regions overlap (150 km membership) so member
counts are not additive; 2026 is truncated at 2026-09-10; b(20 Hz) and the assoc/locate CPU rates are extrapolations /
assumptions; the grid covers only Slab2 30-45 km.

---

# E2 — LFEDetect -> PyOcto -> NonLinLoc pilot (2026-09-10; picks, association and locations complete)
All numbers below are regenerated by `scripts/e2_summary.py` into `results/e2_summary.txt`, which
names the source file for each table. Wording is provisional: this is one pilot region, two windows,
pretrained weights, no fine-tuning.

## What was run
SeisBench **LFEDetect** (Münchmeyer et al. 2024) weights `cascadia` and `jcms` over the 3-C day
files already on disk: CN.PGC 2010-05-01..11-01 (184 d), a 10-station Vancouver Island network for
the 34 densest Lin days (2010-08-10..09-14), an 8-station network for the 2017 ETS
(2017-02-20..04-21), and UW.WISH / PB.B943 / UW.YACT for 2017-2018 (730 d each) — **3,162
station-days per weight set, 6,324 station-day runs over the two weight sets, 9.18 M picks**
(`results/e2_picks_*.parquet`, inventory
`e2_picks_inventory.csv`). Then **PyOcto** association (Crosson 1-D, 36 configurations: 2 sets x 2
models x pick threshold 0.2/0.3/0.5 x >= 3/4/5 stations) and **NonLinLoc 6.04** (Crosson 1-D LAYER
cards, S 0.4 s / P 1.5 s; casc1.6 was NOT used — it is a 17.7 GB ASCII block and Münchmeyer's
released workflow is 1-D).

## The numbers that decide it
- **The picker does see Lin/Bostock LFEs.** At PGC the S picks pile up on the predicted arrival with
  a **20-50 sigma** excess over the chance background, at **-0.12 s** relative to Lin OT + ttS
  (Crosson 1-D) in every distance bin — a sub-second time base. Both catalogues agree.
- **But single-station recall is 12-13 % and precision 0.6-4 %.** Chance-free peak excess
  (+-1 s, thr 0.1): recall **0.130** vs the Lin high-quality subset (residual<0.5, N>=4, n=10,695)
  and **0.123** vs Bostock within 40 km; precision 0.008 / 0.006. In a fixed +-3 s window at thr
  0.3: recall 0.058 (null 0.015), precision 0.023 (null 0.006). Lin et al. report 90 % recall /
  95 % precision for their own (associated) catalogue.
- **Association is what makes it behave like a catalogue.** vi2010 cascadia, thr 0.3: 280 events/day
  at >= 3 stations, 100 at >= 4, 30 at >= 5; **ETS/quiet 9.6 / 17.3 / 34.2** and
  **Spearman(daily events, daily PNSN tremor in the box) +0.93 / +0.88 / +0.87**. vi2017 (no ground
  truth): 126 / 38 / 13 events/day, ETS/quiet 12.0 / 18.0 / 20.0, rho +0.78 / +0.76 / +0.74. The raw
  single-station pick rate only reaches ETS/quiet 1.3-2.7 (PGC) and 1.1-2.1 (off-VI).
- **Most associated events are unverified.** Precision vs Bostock (+-2 s, chance-corrected):
  0.05 (>=3 stations, thr 0.2) -> 0.13 (>=4, thr 0.3) -> **0.22** (>=5, thr 0.3), with recall
  falling 0.063 -> 0.024 -> 0.012. 65 % of >=3-station events are bare triples. Because Bostock's
  catalogue is itself incomplete (18,114 events in the window against Lin's 90,523) this precision
  is a LOWER bound; the same figure against Lin's high-quality subset scaled by that subset's share
  of the Lin catalogue would allow ~0.2-0.4. Either way the majority of events is unconfirmed.
- **Event locations are not usable; the cloud is.** 9,515 events located in 9 min. Per event:
  **25.3 km** median horizontal offset from the time-matched Bostock event (IQR 11-47), 45.3 km from
  the time-matched Lin event, formal ERH 49.8 km. The **epicentre cloud** sits on the LFE belt: the
  median distance to the nearest of the 130 Bostock family hypocentres is **6.2 km against a
  45.9 km uniform-random null**. Location quality is network-limited, not method-limited: distance
  to the nearest family falls 6.6 -> 3.2 km and ERH 57 -> 9 km as the station count goes 3 -> 7.
  > **RETRACTED 2026-09-11 (E3 step 1, Merlin F4).** The uniform-over-the-bounding-box null is
  > wrong: S-only 3-station solutions collapse toward the network centroid, which is where the
  > Bostock families are. Locating the SHIFTED-PICK chance events through the identical setup puts
  > them **8.1 km** from the nearest family against the observed 6.2 km (ratio 0.76) with the same
  > depth pathology. There is no "7x concentration on the LFE belt". And "network-limited, not
  > method-limited" is unsupported: E3 step 2 shows the N>=4 geometry resolves 3.7 km with 0.3 s
  > picks, so the 25 km is the ASSOCIATION.
- **Depth is not resolved from S-only picks**: median 5.8 km where Slab2 is at 28.8 km, 47 % of
  events pinned above 5 km, only 17 % in the 25-45 km interface band, ERZ 20 km, and the offset from
  the matched Bostock depth is a **-24 km bias**, not scatter.
- **The picker run is CPU-bound, not GPU-bound, and scales nearly linearly with concurrency.**
  Profiling one queue (`results/e2_profile.csv`): the **GPU forward is 10-12 % of the wall time**;
  SeisBench's own 1-8 Hz bandpass + decimation to 20 Hz is **33 % at 40 Hz and 75 % at 100 Hz**
  (and F1 says 77 % of fleet 3-C station-years are 100-Hz HH), obspy read 9-48 %, pick
  post-processing 5-7 %. The L40S still earns its place (a CPU-only forward is 2.2-2.3 s/day
  against 0.12-0.13 s) but it is idle most of the time. Scaling (`results/e2_scaling.csv`): the
  per-queue rate falls only 4-5 % from 4 to 12 concurrent queues, GPU utilisation never exceeds
  32 %, and 12 queues use 24 of the pod's 32 allowed CPUs and 10.2 GB of the 46 GB card. Thread
  limits are irrelevant (95.5 vs 90.8 station-days/min per queue with OMP unset vs 2).
  **Fleet projection per weight set, 100-Hz rates:** project fleet (2,168 station-years =
  791,174 station-days) **70.6 h at 4 queues, 24.6 h at 12**; project+extra (4,603 station-years)
  149.9 h at 4, 52.2 h at 12. A preprocessing worker pool feeding one GPU process would raise the
  32-CPU ceiling to ~1,100-1,700 station-days/min = **8-12 h for the project fleet** (a further
  ~2-3x, estimate only); moving the bandpass+decimation onto the GPU would remove the CPU ceiling
  altogether. Either way **the 17.6 TB download dominates** (F1: ~3-15 days at 50-10 MB/s).
  The pilot's own aggregate (47.4 station-days/min per queue over 6,324 station-day runs) was measured
  on a pod shared with the association, location and off-GT jobs and mixes 40-/100-Hz stations;
  use the scaling table for planning.
- **Two data-quality findings.** (i) The `cascadia` weights fire on the **first window of a day
  file** (0.39-0.67 picks per station-day in the first 5 s, 5-10x the expected rate; `jcms` 10-30x
  cleaner). Because every station is cut on the same UTC boundary this manufactured a spurious
  zero-lag WISH x B943 coincidence (5.9x background at a lag the geometry forbids); it is removed by
  dropping the first/last 60 s of each day, and it does not affect daily rates or the associated
  events. (ii) The **degenerate-window QC hole that poisoned D3 does not transfer to the picker**:
  0.10-0.24 % of cascadia picks (0.03-0.04 % of jcms) fall in gap-affected 60-s windows against a
  0.75 % baseline, and 0 % in fully empty windows — LFEDetect does not fire on zero-filled data.
- **A catalogue-time-base error worth flagging.** The Slow-Earthquake-DB `Bostock2015-LFE` time is
  **not an origin time on our time base**: our picks peak at **+12.4 s** after the predicted arrival
  built from it, and the offset is more nearly constant WITHOUT the travel time (+24-26 s from the
  observed S arrival at PGC) than with it. An independent Lin-vs-Bostock comparison shows the same
  sign and size (+10..+16 s). The previous calibration (C = +0.086 s from the median offset to the
  nearest Lin event) was chance-dominated: 16,026 of 18,114 Bostock events have a Lin neighbour
  within +-20 s by chance. All Bostock-referenced numbers here use the measured offset and are
  lower bounds; the convention should be settled with the catalogue's authors.

## Verdict
> **PARTLY OVERTURNED 2026-09-11 by the E3 controls — read the "E3" section at the end of this file
> first.** The two things this verdict PASSES on (the ETS/quiet 5-34x contrast and rho ~0.9) are
> reproduced by a shifted-pick chance null and carry no information about whether the events are
> real; the E2 operating point (thr 0.2-0.3, N>=4) is 33-59 % chance. The things it FAILS on are
> mostly benchmark artifacts: the Lin 90/95 % comparison is between two different quantities, and
> the 25 km hypocentre error is the association step, not the picker or the velocity model.

**The pilot passes on activity, throughput and structure, and fails on per-event quality.**
Association gives a tremor-tracking event series (ETS/quiet 5-34x, rho ~0.9) at a fleet cost of
25-70 h on this one GPU — but with the pretrained weights on this geometry most events are unconfirmed (precision
0.13-0.22 against Bostock, a lower bound),
recall is a few per cent of the reference catalogues, and hypocentres are 25 km with no depth. That
is 1-2 orders below Lin et al. 2024 and below the PLAN's GOAL (Thomas et al. 2025 at PNW scale
needs event-level times AND few-km locations, family assignment, moment, stable completeness).
**Recommendation: fine-tune LFEDetect on the substrate already on disk (1,817 Bostock stacks with
P and S picks for 129/130 families, Shelly MTJ, Ducellier, plus our hard negatives and the
first-window artifact) before committing the 17.6 TB fleet download; operating point N >= 4 at pick
threshold 0.2-0.3; keep both weight sets until fine-tuning decides; do not publish S-only depths.**

Files: scripts/e2_*.py + e2_run_*.sh; results/e2_* (~90 files, index in results/e2_summary.txt);
figures/e2_*.png (5); ledger NOTES.md '2026-09-10 E2 pilot (resumed after the 21:03 pod crash)'.

---

# E3 (2026-09-11) — the on-disk controls. The E2 scorecard inverts again.

**Nothing was downloaded except one paper PDF. No E2 output was deleted or overwritten.**
Ledger with every run, number and bug: `NOTES.md`, entry "2026-09-11 E3 steps 0-6".

## STEP 1 VERDICT (stated first, because everything else depends on it)

**The E2 activity PASSes are RETRACTED. At the operating point E2 recommended, a large fraction of
the associated events is chance coincidence, and the ETS/quiet contrast that was the main evidence
for the catalogue being real is reproduced, ratio for ratio, by pure chance.**

A per-station constant time shift (|shift| uniform in 300-3,000 s, independent per station) leaves
every station's pick count and daily rate — and therefore the ETS modulation — exactly intact while
destroying real inter-station coincidence. Anything PyOcto still finds is combinatorial. Five seeds
per configuration; a second variant with per-station-per-hour shifts (Münchmeyer's own control
design) changes the answer by 4-14 % relative.

| set | thr | N>= | observed ev/day | null ev/day (± s.d.) | **null / observed** | branch |
|---|---|---|---|---|---|---|
| vi2010 | 0.2 | 3 | 669.5 | 616.0 ± 3.7 | **0.92** | chance-dominated |
| vi2010 | **0.2** | **4** | 279.1 | 163.2 ± 3.5 | **0.59** | **chance-dominated** |
| vi2010 | 0.2 | 5 | 99.1 | 22.6 ± 0.6 | 0.23 | unresolved |
| vi2010 | 0.3 | 3 | 271.9 | 199.1 ± 2.3 | **0.73** | chance-dominated |
| vi2010 | **0.3** | **4** | 96.1 | 31.8 ± 1.2 | **0.33** | **unresolved** |
| vi2010 | 0.3 | 5 | 27.7 | 2.43 ± 0.35 | **0.088** | PASS |
| vi2010 | 0.5 | 4 | 15.0 | 1.55 ± 0.17 | **0.103** | PASS |
| vi2010 | 0.5 | 5 | 3.50 | 0.059 ± 0.029 | 0.017 | PASS |
| vi2017 | 0.3 | 3 / 4 / 5 | 126.4 / 36.2 / 9.45 | 88.8 / 10.5 / 0.70 | 0.70 / 0.29 / 0.074 | — |
| vi2010 jcms | 0.3 | 3 / 4 / 5 | 114.6 / 33.6 / 8.97 | 66.6 / 6.69 / 0.28 | 0.58 / 0.20 / 0.031 | — |

The audit's rule is ">=50 % at N>=4 -> chance-dominated; <=15 % -> PASS stands; between ->
unresolved, report the ratio". **E2's two recommended operating points land at 0.59 and 0.33**, i.e.
one chance-dominated and one unresolved. Every N>=3 configuration is 0.42-0.92: the N>=3 catalogue
is chance. **There IS a clean operating point — pick threshold 0.5 with N>=4 (15 events/day, 10 %
chance) or threshold 0.3 with N>=5 (28 events/day, 8.8 %)** — which is also where Münchmeyer et al.'s
own scrambled-pick control sits (5-10 %).

**ETS vs quiet, separately.** Merlin's F10 had to be fixed first: with the 30-day rolling PNSN
tremor statistic built on the FULL record (not inside the 34-day window with `min_periods=5`),
**all 34 days of vi2010 are ETS days and none is quiet** — E2's "ETS/quiet 34x" was entirely a
window-edge artifact and the ratio is undefined there. vi2017 has 35 ETS and 4 quiet days:

| N>= | obs ETS | null ETS | obs quiet | null quiet | obs ETS/quiet | **null ETS/quiet** |
|---|---|---|---|---|---|---|
| 3 | 196.0 | 139.6 | 36.0 | 25.8 | 5.4× | **5.4×** |
| 4 | 57.9 | 17.5 | 5.50 | 1.40 | 10.5× | **12.5×** |
| 5 | 15.3 | 1.17 | 0.75 | 0.10 | 20.4× | **11.7×** |

A per-station activity contrast of r becomes r^N in an N-fold chance catalogue, so the ETS/quiet
contrast and the ρ ≈ 0.9 tremor correlation of the associated events carry **no information** about
whether the events are real. Corroboration of the mechanism: the fitted chance exponent
(d log null-rate / d log network-pick-rate) is **2.16 / 3.11 / 3.97 for N>=3/4/5**, i.e. exactly N.

**The location "cloud" statement is also retracted (F4).** Locating the null events through the
identical NonLinLoc setup: the observed epicentres sit a median **6.2 km** from the nearest of the
130 Bostock family hypocentres — and the **chance events sit 8.1 km away, a ratio of 0.76**. The
"45.9 km uniform-random null → 7.4× concentration on the LFE belt" was an artifact of drawing the
null uniformly over a bounding box; S-only 3-station solutions collapse toward the network centroid,
which is where the families are. The null's depth pathology is the same (median 4.5 km, 57 % above
5 km depth).

Figure `figures/e3_null_rates.png`; tables `results/e3_nullevents_rates.csv`, `e3_nullhr_rates.csv`.

## STEP 2 — where E2's 25 km actually comes from: the association, not the picker

The prescribed method could not be used: **every SeisBench stack trace is aligned on its own P
arrival at sample 800**, so the dataset carries no cross-station time base and only the S−P
differential time survives. Four configurations replaced it (Merlin re-specified after being shown
the evidence). Figure `figures/e3_relocation.png`.

* **A — real stack picks, S−P grid search, same Crosson 1-D model as E2** (672 both-picked traces,
  SHVB excluded for bad metadata coordinates): median **3.1 km** horizontal, **−3.6 km** depth for
  families with ≥6 stations and gap <200°. **Passes the audit's 5 km bar: the model and the method
  are not the problem.** The depth bias matches Crosson over-predicting S−P by ~0.3 s.
* **C — synthetic S-only through the exact NLL control**, times read from NonLinLoc's own grids so
  the model error is identically zero; all 42 three-to-six-station subsets, 10,836 locations per
  noise level. At **N=3 the answer is independent of pick noise** (17.7 km / −27.6 km at σ=0.05 s;
  18.6 km / −28.0 km with heavy tails; a third of solutions pinned to the grid top) — **E2's −27 km
  depth bias and shallow pile-up are 100 % geometry.** But at **N≥4 the geometry is fine**: 3.7 km
  at the measured 0.3 s pick noise, against E2's 25 km.
* **D — real picks, known events, no PyOcto**: 6.97 km / −30.5 km at N=3, against E2's 33 km. Real
  per-station S-pick residual s.d. is **0.40–0.56 s**, not the 1–3 s assumed.

**The ladder C → D → E2 puts the dominant term in the association step**: same picker, same data,
same locator, and the only difference between D (7 km) and E2 (33 km) is PyOcto. That is the same
conclusion step 1 reaches from the rate side, by independent evidence.

## STEP 3 — recall is SNR-limited, but the injection test cannot measure the ceiling
> **SWAP-AFFECTED (E4 step 0).** The injected waveform was inserted with its horizontals
> swapped (`e3_inject.py:132`, fixed), so the picker saw the LFE 90 deg rotated on the
> horizontals. The injected recall curve below is therefore a lower bound of unknown size, and
> the "11x shortfall against the real curve" is partly this bug and partly the x-axis anchor
> (Merlin F4: the 0.078 anchor is a COHERENT stack amplitude, not a single-event incoherent
> s/n). Re-run with the fix and a directly measured incoherent s/n: **E4 step 3**.

1,800 injections per station (Bostock stacks at the host station, 6 amplitude levels, real quiet
days, both weight sets, un-injected baseline measured and subtracted). Figure
`figures/e3_injection.png`.
**Neither branch of the audit fires: the injected curve lies far BELOW the real one** — at matched
s/n (rms over 6 s / rms of the 1–8 Hz noise, the definition behind D3's 0.078) the real Bostock
recall is ~11× the injected recall. Matching the curves implies a real LFE behaves like s/n ≈ 0.22–0.25,
about 3× the D3 anchor.
The obvious explanation was tested and **rejected**: splitting the real population by local event
density gives isolated LFEs (≤1 neighbour within ±30 s) recall 0.128/0.171/0.197 versus dense
(≥10) 0.100/0.114/0.279 — burst context does not explain the gap. What remains is that the 0.078
anchor is low and/or that a family **stack** is a pessimistic proxy for a single event (stacking
jittered events smears the waveform). **Consequence: the injection test cannot be used to claim the
picker is at its SNR ceiling, and any fine-tune must be judged on real single-event positives, not
on the stacks** (which are in-sample for the `cascadia` weights anyway).
**False alarms and the lever.** F on real quiet days at PGC (cascadia): 720 / 247 / 113 / 32 S picks
per station-day at thr 0.1 / 0.2 / 0.3 / 0.5. For the N>=4 chance rate to fall to ≤10 % of the real
rate, **F must fall by 2.3× at thr 0.2, 1.7× at thr 0.3, and is already there at thr 0.5.** The
required factor is ~2 and raising the pick threshold already delivers it.

## STEP 3b — denoising before picking: DeepDenoiser HURTS, spectral Wiener is neutral
> **PARTLY SWAP-AFFECTED (E4 step 0).** The recall and false-alarm rows inherit the injection
> component swap (`e3_denoise.py:194, 317`, fixed). The **fidelity finding stands** -- the
> rms-kept / waveform-correlation / spectral-ratio measurements compare each arm against the
> same injected signal component by component, so they are orientation-independent:
> DeepDenoiser still deletes 88-93 % of the LFE. The conclusion "denoising is not a lever" is
> unchanged.

Three arms, identical days, injections, band-pass and sampling rate (100 Hz throughout).
DeepDenoiser `original` was already cached offline. Figure `figures/e3_denoise.png`.

| | rms of the LFE kept | waveform corr. | spec. ratio 1–4 Hz | F (PGC, thr 0.3) | real-day recall vs Bostock (thr 0.1) | vs Lin-hq |
|---|---|---|---|---|---|---|
| raw | 1.00 | 1.00 | 1.00 | 113 /day | **0.134** | **0.111** |
| **DeepDenoiser** | **0.07–0.12** | **0.44–0.60** | **0.04–0.06** | **956 /day** | **0.034** | **0.010** |
| spectral Wiener | 0.92–0.94 | 0.96–0.99 | 0.92 | 120 /day | 0.131 | 0.111 |

**DeepDenoiser removes 88–93 % of the LFE amplitude and strips it hardest exactly in the 1–4 Hz band
where LFE energy lives, while multiplying the false-alarm rate by 8.5×** — worse on both axes at
once. It is trained on impulsive local earthquakes and treats an emergent low-SNR LFE as noise. The
classical spectral-Wiener baseline preserves the waveform and lands within a few per cent of raw on
every metric. **Denoising is not a lever for this problem**; a useful denoiser would have to be
trained on LFEs, and step 3 says the limit is signal amplitude rather than removable additive noise.

## STEP 4 — **RETRACTED 2026-09-11 (E4 step 0): COMPONENT SWAP.** ~~Bostock templates as a classifier:
## a weak excess, and D3 stands~~

> **RETRACTION.** Every number in this section was computed with the two horizontal components
> of the DATA swapped relative to the templates. The SeisBench stacks are stored `Z12` with
> **1 = N, 2 = E**; `e3_template_classifier.py` built the data matrix as `[Z] + sorted(non-Z)`,
> which on CN stations is BHZ, BHE, BHN = **Z, E, N**. Independently verified here
> (`scripts/e4_polarization_check.py`, `results/e4_polarization_check.csv`): on 327 stacks with
> P snr > 6 dB and rectilinearity > 0.7 the horizontal P principal axis agrees with
> `path_backazimuth_deg` to a median **7.9 deg** under 1 = N (90 % within 25 deg) against
> **37.3 deg** under 1 = E (36 %). So the templates' N was dotted with the data's E throughout.
> The conclusions "Bostock's own events register at 0.60 MAD with their own templates" and
> "D3 stands" are therefore **UNTESTED, not established** (Merlin re-audit F1, F2, F3). The
> placement was also network-coherence-only with one common lag and an 8-s template; the
> corrected test is **E4 step 1** below. Fixed at `scripts/e3_template_classifier.py:190`
> (old behaviour reproducible with `E3_LEGACY_CHAN_ORDER=1`). The `mad_median` column of
> `results/e3_classifier_summary.csv` is separately buggy (normalised by the best-of-129
> background) -- see `results/e3_classifier_summary.WARNING_mad_median.txt` and the corrected
> copy `results/e4_e3_classifier_summary_relabelled.csv`.


Network cc of the 3-C data at the 6 overlap stations against all 129 family templates, one common
origin-time lag over ±3 s (grouped `conv1d`, GPU). Only 362 of 774 (family, station) templates
exist, so the effective network is **2.8 stations**. Figure `figures/e3_classifier.png`.
* **Associated events show a real but small excess**: 3.7 % exceed the chance-event 99th percentile
  against 1.1 % by construction for the null and 1.8 % for random times; 23.3 % reach 6 MAD of their
  own template's background against 14.4 % for random times. Roughly **one associated event in ten**
  carries a family match above background — a hybrid detect-then-classify product is weakly viable.
* **Bostock's own catalogued events register at only 0.60 MAD with their own family templates**
  (p90 3.66; 2.3 % above 6 MAD). The "≥6–8 MAD → D3's 3.5 MAD was a template-quality artifact"
  branch does not fire. Scaled to Bostock's 21 stations this is ~1.0–1.5 MAD — **D3's retirement of
**Placement check and the honest caveat.** The best-fitting common lag is nearly UNIFORM over the
+-3 s search in every set (fraction with |lag| < 0.5 s: assoc 0.155, bostock 0.157, null 0.138,
random 0.143, against 0.167 for a uniform distribution; peak/median bin 1.2-1.4 everywhere). So no
set shows a coherent origin-time solution -- which is itself the finding, but it means this test
cannot by itself separate "the templates genuinely do not register on 2.8 stations" from "a
residual error in the template placement". Two things argue for the former: the own-template cc is
shifted the RIGHT WAY (+0.65 MAD, median 0.1164 against the per-family random background 0.1059 --
a grossly wrong placement would give zero shift), and Muenchmeyer et al. report the same null
result for Cascadia in their Figure S10. A stronger placement test (per-station single-template cc
at Bostock times against random times, which does not require network coherence) was NOT run.

  single-template matched filtering stands**. *(RETRACTED with the rest of step 4. The clause
  "and it agrees with Münchmeyer's own Figure S10" is also WITHDRAWN on its own merits, per
  Merlin F2: S10 reports that their UNCATALOGED DL detections correlate with family stacks like
  noise. It says nothing about CATALOGUED Bostock events failing against their own templates --
  those exceeded Bostock's network-cc threshold by construction, so a null there is an
  implementation red flag, not corroboration.)*

## STEP 5 — the Bostock time base is a per-FAMILY property, and it is now measured

Per-family peak of (PGC S pick − DB time − ttS): **90 of 100 families accepted, median +11.75 s,
s.d. 4.14 s, range +2.50 to +20.50 s.** The offsets do vary by more than 5 s across families, but
the variation is real, not noise: **between-family sum-of-squares 0.565 of the total versus
between-station 0.009** (per-station offsets are all +11.25 to +13.00 s), and the PGC offsets
cross-validate against ≥3 other stations in an independent window at **r = 0.959, residual s.d.
1.19 s**. F7's apparent station dependence was the residual of a family mixture.
Re-scoring with the table raises the chance-corrected Bostock precision by 22–41 % relative
(N≥4: 0.129 → 0.158; N≥5: 0.218 → 0.307) and sharpens the residual peak to 3.2–6.8× background.
**The convention itself is still unidentified** — no single-station reference-arrival model reduces
the scatter (best candidate 3.80 s versus 4.14 s raw) — so it needs the paper or the DB maintainers;
Lin-hq stays primary. The chance-dominated `C` estimator in `e2_groundtruth.py` (16,026 of 18,114
Bostock events have a Lin event within 20 s *by chance*) has been removed from the script; the CSV
it wrote is left untouched with a sibling `*.WARNING_chance_dominated.txt` beside it. Figure `figures/e3_family_offsets.png`.

## STEP 0 — Münchmeyer et al. (2024): we are where the picker's authors are

`results/e3_munchmeyer_numbers.md`. Their Cascadia DL catalogue is **8,369 events over 59 days
(141.8/day) against 148,650 Bostock events over the identical windows (2,519.5/day) — a ratio of
0.056, i.e. ~18× smaller**. Our vi2010 N≥4 thr-0.3 run gives 100 events/day against 2,663
Bostock/day: **the same place on the same curve.** Their headline "60 % of our LFEs have a reference
LFE within 10 s" has a chance floor we measured at **28.3 %**, so it is 2.1× chance. They **do**
report a chance-association control (per-station-per-hour shifts → 5–10 % of the catalogue), which
is the direct precedent for step 1 and agrees with the audit's ≤15 % bar. Their PyOcto parameters
are in a supplement that is not in the arXiv package (404) — **an open item, not a guess.**

## STEP 6 — documentation corrections

**Off-ground-truth (F6).** Merlin's numbers reproduce exactly. PNSN tremor days in 2017–18 within
40/60/80/100 km: **WISH 26/62/104/134, YACT 19/44/72/94, B943 109/146/173/210.** WISH and YACT have
almost no tremor near them and cannot test generalisation. Inside the 2017 source zone the SURFACE
stations respond — ETS/quiet **SYMB 2.26, PGC 2.25, VGZ 2.02, CLRS 1.99, SNB 1.81** against the
co-located borehole **B011 1.78** — and in the direct 0.6 km surface/borehole control PGC vs B011
recall is 0.102 vs 0.125 (thr 0.1) and 0.033 vs 0.059 (thr 0.3), i.e. **1.2–1.8×, not categorical**;
the best single station of all is the surface station **SNB (0.167)**. *"Surface flat, borehole
responds" is distance to the source, not instrument type, and must not seed a borehole-only design.*
**Lin attribution (F2).** Lin's 90 %/95 % are held-out balanced-test-window numbers with an equal
noise class, not continuous-record catalogue numbers; his published catalogue is a 15 s N≥3
coincidence catalogue of which his own quality cut keeps 9.8 %. E2's "FAIL vs Lin's 90 %/95 %"
compared two different quantities and is withdrawn.
**F1 readiness at 40 km (F8).** `results/f1_readiness_40km.csv`, recomputed from F1's own tables
(method validated by reproducing F1's stored n60 exactly in 75–85 % of cell-years, within ±1 in
92–95 %). The measured usable radius is ~40 km (single-station recall 16 % <25 km, 11 % 25–40 km,
3.2 % 40–60 km). Pooled over all 137 cells, project fleet: **fraction with ≥4 3-C stations
0.34 at 40 km vs 0.92 at 60 km in 2018, and 0.78 vs 1.00 in 2024; with ≥6, 0.08 vs 0.56 in 2018 and
0.34 vs 0.985 in 2024.** F1's report is otherwise unchanged.

## What E3 changes

1. **The operating point moves** to thr 0.5 / N≥4 or thr 0.3 / N≥5 (8.8–10 % chance, 15–28 ev/day).
   The E2 recommendation (thr 0.2–0.3 / N≥4) is 33–59 % chance.
2. **Do not fine-tune and do not start a fleet download on this evidence.** Step 3 shows the
   injection test cannot establish an SNR ceiling and that a fine-tune must be judged on real
   single-event positives; step 3b shows denoising is not the lever; step 1 shows the needed
   false-alarm reduction is only ~1.7–2.3× and the threshold already delivers it.
3. **Event hypocentres are an association problem, not a picker problem** (step 2's C→D→E2 ladder),
   and event-level depth from S-only picks is unrecoverable at N=3 by geometry alone.
4. **Prescription step 6 (network back-projection of the S-probability traces) is now triggered** by
   step 1's verdict. It was conditional and was not run here.

# E4 (2026-09-11) — the components were swapped; with them fixed, Bostock's own templates DO register

**Nothing was downloaded except the PMC text of Thomas et al. 2025 and FDSN station METADATA.
No E2/E3 output was deleted or overwritten.** Ledger with every run, number, bug and fix:
`NOTES.md`, entry "2026-09-11 E4". Merlin's mid-run vet
(`results/merlin_vet_e4_redesign_2026-09-11.md`) replaced the step-1 binary bar with a sensitivity
model, the hold-out with alternate 2-h blocks, and added the decimation and prediction legs; both
the original bars and the model are reported.

## STEP 1 VERDICT (stated first) — the template route is OPEN on our station class

**E3 step 4's "Bostock's own events register at 0.60 MAD with their own templates" was an artefact
of the component swap and a bad placement. With Z,N,E fixed, a per-station empirical moveout and
per-template backgrounds, Bostock's catalogued LFEs register at a network median of 5.3-6.2 MAD on
SIX 2010 stations, against a random-time background of 2.73 MAD — an excess of +2.6 to +3.4 MAD —
and 12-26 % of them exceed an 8-MAD network threshold at which the per-station time-shifted control
fires on 0.00-0.02 % of the same events.**

2010 outer ring PGC/LZB/VGZ/SNB/YOUB/PFB (C8.SOKB has NO data on disk: 34 zero-byte `.mseed.nd`
files), 34 days, 15,596 catalogued events, 88 families with a template and an accepted time-base
offset, median **3** of the 6 stations carrying a template for a given family (9 channels).
Figure `figures/e4_step1_registration.png`; tables `results/e4_step1_{registration,network,channels,magslope}.csv`.

| template | net MAD median (M1.8-2.1) | excess over random | p90 | >= 6 MAD | **>= 8 MAD** | shifted control >= 8 | detections/template/day at 8 MAD (median over family-days; RENAMED E5 step 0 -- these are TOTAL declustered detections, not false alarms; quiet-day FA is in the E5 section) |
|---|---|---|---|---|---|---|---|
| 2 s | 4.41 | +1.66 | 6.77 | 0.181 | 0.037 | 0.0000 | 0 |
| 4 s | 5.30 | +2.57 | 8.24 | 0.369 | **0.117** | 0.0000 | 1 |
| 6 s | 5.65 | +2.92 | 8.97 | 0.440 | **0.168** | 0.0002 | 1 |
| 8 s | 5.89 | +3.17 | 9.51 | 0.481 | **0.206** | 0.0002 | 2 |
| 12 s | 6.15 | +3.42 | 10.27 | 0.525 | **0.256** | 0.0002 | 4 |

**Against the original brief's acceptance** ("median >= 3 MAD for M >= 1.8 with 2-4 s templates AND
recall >= 30 % at 8 MAD with chance <= 10 %"): the MAD half passes at every length >= 4 s (5.3-6.2
raw, +2.6-3.4 over background), the chance half passes by three orders of magnitude (0.02 % against
a 10 % bar), and **the recall half does not** — 12 % at 4 s, 21 % at 8 s, 26 % at 12 s against 30 %.
This is the *intermediate* case the brief asked to be reported with numbers. It is not the
FAILURE branch: failure was "median < 1.5 MAD after the fix and the empirical moveout".

**Against Merlin's sensitivity model** (his prediction for a correct implementation was a median
1.5-2.5 MAD per station for M >= 1.8 at 4-8 s): measured per-station excess at M1.8-2.1, 8-s
template, with the prediction `q x r1 / MAD_scale` (q = 0.89, r1 = 10^(S_snr_db/20)/sqrt(N_stack)):

| station | epi (km) | r1 | MAD scale | measured excess | predicted | **measured/predicted** | Z-only excess |
|---|---|---|---|---|---|---|---|
| LZB | 23.8 | 0.062 | 0.050 | 2.34 | 1.11 | **2.11** | 0.74 |
| PGC | 37.5 | 0.052 | 0.046 | 2.15 | 1.00 | **2.16** | 0.44 |
| VGZ | 36.1 | 0.045 | 0.048 | 1.15 | 0.84 | 1.37 | 0.13 |
| PFB | 41.9 | 0.033 | 0.067 | 0.68 | 0.43 | 1.58 | -0.39 |
| YOUB | 52.0 | 0.058 | 0.054 | 1.87 | 0.95 | **1.96** | 0.66 |
| SNB | 55.5 | 0.040 | 0.064 | 0.98 | 0.56 | 1.75 | 0.15 |

Measured registration is **1.4-2.2x the stack-SNR prediction** (median 1.9x). That is Merlin's
">= 2x -> the stack-SNR convention understates r1" signature, and step 3 confirms it from a
completely independent measurement: the per-event INCOHERENT s/n at PGC, measured directly with no
stacking, is **0.57-1.02**, not the 0.078 coherent-stack anchor the project has used since D3.

**Three structural facts the model needs, all measured here.**
1. **The LFE is on the HORIZONTALS.** Single-channel (Z) excess against its own random-time
   background is 0.13-0.46 MAD while the 3-C mean is 0.79-1.57 — a gain of **3.4-6.1x, not
   sqrt(3) = 1.73**. PFB's Z-only excess is *negative*. This is exactly the quantity the component
   swap destroyed.
2. **Longer templates keep winning.** The excess rises monotonically 2 -> 12 s at every station and
   on the network. The LFE signal is coherent over at least 12 s, so it is not a 2-s S packet; that
   is also why a 1-s-resolution phase picker under-performs a matched filter here.
3. **The MAD-vs-magnitude slope is FLAT, and it was measured, not assumed.** Within family, the
   excess rises by only 1.28 MAD per magnitude unit on a base of 1.45 MAD, i.e. an amplitude
   scaling of **10^(0.38 dM)** where moment scaling would give 10^(1.50 dM). Between families the
   slope is *negative* (distant families need bigger events to be catalogued). The independent
   incoherent-s/n measurement gives the same answer: 0.57 at M1.45 to 1.02 at M2.00 = 10^(0.46 dM).
   Bostock's catalogue magnitudes compress the true amplitude range by ~3-4x in log terms.

**Empirical moveout matters.** The median best single-station lag is **+0.41 s (IQR +0.08..+0.64,
s.d. 0.39 s)** relative to Crosson ttP + the stack's own S-P — exactly Merlin's F3 prediction. Pass 2
(the moveout applied) lifts the network excess at 8 s from +2.41 to +3.34 MAD and the 8-MAD recall
from 14.0 % to 22.3 %.
**The moveout is not overfitting.** In step 2's leg A the moveout was measured on the BUILD half of
the alternate-2-h-block split and scored only on the TEST half: the network excess at 8 s,
M1.8-2.1, is **2.83 MAD with 16.8 % at 8 MAD**, against 3.17 / 20.6 % when the moveout is measured
on the whole window. Same class, fully held out.

## STEP 2 — the q-control passes; the borehole transfer works; the network sum does not

Hold-out = alternate 2-h blocks with a 60-s edge guard (the Aug-build/Sep-test split of the brief is
impossible: the 2010 ETS is over by 1 September). Each leg measured its empirical moveout on the
BUILD half and was scored on the TEST half. Figure `figures/e4_step2_transfer.png`.

**(a) q-control — our stacking reproduces Bostock's templates.** Re-stacking each ring station's own
2010 data at the build-half Bostock times (per-family offset + Bostock's own S-P) for 82 families
with >= 40 build-half events (median **90** events per stack) gives
**cc(our 2010 re-stack, SeisBench stack) = 0.73-0.86 by station** (PGC 0.79, LZB 0.83, VGZ 0.73,
SNB 0.84, YOUB 0.86, PFB 0.78), against Merlin's >= 0.6 bar and his "< 0.3 -> the offset sign or the
S-P application is wrong, stop" failure branch. **The E3 step-5 time base and the placement are
correct.** Merlin's literal bar (N >= 200) cannot be evaluated -- the build half gives at most 167
events per family -- but at **N >= 150** (32 pairs) the median cc is **0.84** and 84 % exceed 0.6.
His stop-line (cc < 0.3 at large N) fires on **five pairs, and all five are the SAME family (15) at
five different stations**, with stack S-SNR 1.0-5.1 dB and an S pick 0.9-2.2 s off the prediction:
a single bad family, not a global sign error. Family 15 should be blocklisted. The held-out registration ratio (our 90-event template / the published stack) is
0.45-0.90 by station against a predicted q-ratio of 0.34-0.57 — **our short-stack templates are
BETTER than q(N) = 1/sqrt(1+1/(r1^2 N)) predicts**, again because r1 is understated.

**(b) transfer to stations Bostock never had.** Stacks built purely from 2010 data at Bostock times
(Crosson-predicted S, +-20 s, iteratively cc-aligned, S picked on the stack envelope), N ~ 90:

| station | class | epi (km) | stack S-SNR (dB) | same-N random null (dB) | held-out excess (MAD, 8 s) | vs the co-located PGC 2010 template |
|---|---|---|---|---|---|---|
| B011 | borehole, 0.3 km from PGC | 40.9 | **8.60** | -0.03 | **1.19** | **0.85** |
| B926 | borehole | 46.9 | **9.66** | 0.19 | **1.87** | 1.33 |
| LRIV | surface bb | 51.9 | 2.76 | 0.19 | 0.39 | 0.28 |
| STW | surface, **Z only** | 41.4 | 0.06 | -0.08 | 0.04 | 0.03 |
| (ring reference) PGC / LZB / YOUB | | | 11.2 / 11.8 / 10.2 | 0.08 / -0.32 / -0.10 | | |

**B011 passes Merlin's bar (>= 0.7x the co-located PGC 2010-only re-stack) at 0.85, and B926 exceeds
it.** ~90 events are enough to build a usable borehole template; LRIV (52 km) and the single-channel
STW are not usable. **How many events a template needs:** S-SNR rises from ~4 dB at N ~ 40 to
~10-12 dB at N ~ 200 on the good stations; the transfer is limited by distance and channel count,
not by the number of events available.

**(c) the network gain is where it fails.** Held-out network excess at 8 s, all magnitudes:
ring with published stacks **2.91 MAD**; ring with our 2010 re-stacks 2.81; **ring + B011 + B926 +
LRIV + STW 2.64** — *adding four stations makes the network slightly WORSE*. The cause is arithmetic,
not physics: an equal-weight channel mean has `e_net = (sum_s e_s)/sqrt(n)`, so a station helps only
if its own excess exceeds about half the current mean. B011 (1.19) and B926 (1.87) clear that bar;
LRIV (0.39) and STW (0.04) do not and dilute the sum. The SNR-weighted optimum
`sqrt(sum_s e_s^2)` would have risen from 3.70 to 4.31 MAD (+17 %) by adding them.
**Measured network / equal-weight prediction = 0.86 (ring) and 0.82 (ring+transfer)** -- within
Merlin's "30 % of sqrt(sum e_s^2)/1.38" band in both cases, so the combination is behaving as
designed; the design is the problem. **An SNR-weighted network stack is now a required change.**

## STEP 3 — the injection is FIXED and RE-ANCHORED, and it now lies ABOVE the real curve

`scripts/e4_inject.py`. Three fixes: components; the x-axis anchor; and own N-event stacks as the
injected waveform. Figure `figures/e4_step3_injection.png`.

**The anchor, measured directly and with a control.** Per-event incoherent s/n at PGC =
sqrt(rms^2(S-1..S+5 s) / rms^2(the same window 20 s earlier) - baseline), 1-8 Hz, no stacking, for
the 7,176 Bostock events with epi_PGC < 40 km. Control: the identical statistic at 4,000 random
times on the same days gives a power ratio of **1.077**, and at the Bostock times shifted by
300-3000 s, 1.102. Against that baseline the measured s/n is **0.57 (M1.45), 0.70 (M1.55),
0.78 (M1.70), 0.96 (M1.87), 1.02 (M2.00), 1.46 (M>=2.1)** -- **7-13x the 0.078 coherent anchor.**
Isolated events (<= 1 other Bostock event within +-30 s) give 0.59-0.97, essentially the same as
dense ones (0.82), so this is NOT burst contamination.

**Injected recall at PGC (cascadia, thr 0.1, chance-corrected against the same un-injected day):**

| s/n_rms | SeisBench family stack | our own stack, all events | our own stack, N=20 |
|---|---|---|---|
| 0.2 | 0.133 | 0.117 | 0.067 |
| 0.3 | 0.450 | 0.350 | 0.117 |
| 0.5 | 0.900 | 0.667 | 0.450 |
| 1.0 | 0.983 | 0.967 | 0.817 |

Real recall at PGC is 0.109-0.256 over M1.44-1.90. Placed at the MEASURED anchor (s/n 0.57-1.02),
the real events should be recovered 67-98 % of the time; they are recovered 11-26 %.
**The acceptance ("injected and re-anchored real curves agree within 2x") FAILS, in the opposite
direction from E3** — where E3 had the injection 11x BELOW the real curve, it is now ~4-6x ABOVE it.
The curves cross at s/n ~ 0.2-0.3, so the picker behaves as if a real LFE were a s/n 0.2-0.3
*stack-shaped* signal: 3-4x the D3 anchor, 3x below the measured window energy.
~~**Per the brief's failure branch, this REOPENS "fine-tune on template-aligned real positives with a
family-held-out split by `source_family_id`" as an option. No fine-tune was run.**~~ **WITHDRAWN (E5 step 0, 2026-09-11; results/merlin_audit_e4_2026-09-11.md F7):** the 0.57-1.02 anchor is the INCOHERENT S-window power excess while the matched filter sees the coherent part (~0.10-0.20 in cc, 1.3-2.5x the 0.078); anchored on coherent amplitude the injected recall 0.13-0.45 matches the real 0.11-0.26 within 2x, so the E3-audit acceptance PASSES and fine-tuning is CLOSED (deferred indefinitely). The physical
reading, consistent with step 1, is that a real single LFE has the same 6-s window energy as a stack
but spread over >= 8-12 s and variable in shape, so a 1-s-resolution picker trained on aligned
stacks under-reads it while a long matched filter does not.

## STEP 4 — tighter association tolerances halve the chance fraction

`scripts/e4_run_assoc_null.sh`, `e4_assoc_null_score.py`; 32 PyOcto runs at pick_match_tolerance
1.2 s, edt_pick_std 0.5 s, zlim 20-50 km, `n_picks` run DIRECTLY at 4 and 5.
Figure `figures/e4_step4_assoc_null.png`.

| thr | N | obs/day | null/day | **null/obs** | hourly-shift variant | E3 (tol 2.5 s, std 1.5 s, z 0-70) |
|---|---|---|---|---|---|---|
| 0.3 | 4 | 71.5 | 12.54 +- 0.44 | **0.175** | 0.200 | 0.330 |
| 0.3 | 5 | 19.0 | 0.49 +- 0.08 | **0.026** | 0.038 | 0.088 |
| 0.5 | 4 | 11.3 | 0.52 +- 0.06 | **0.046** | 0.052 | 0.103 |
| 0.5 | 5 | 2.41 | 0.03 | **0.012** | 0.024 | 0.017 |

The acceptance (null/obs <= 0.10 at N >= 4 thr 0.3 with obs >= 60/day) is **not** met — 0.175 at
71.5 events/day — but the chance fraction is halved at every operating point while 74 % of the
events are kept, and thr 0.3/N>=5 (0.026, 19/day) and thr 0.5/N>=4 (0.046, 11/day) are now
**cleaner than Muenchmeyer's own 5-10 % scrambled-catalogue standard.** Bostock recall
(per-family offsets, +-2 s, chance-corrected) at thr 0.3/N>=4: 0.7 % / 2.9 % / 4.5 % / **9.4 %**
by magnitude bin, 3.4 % overall. **The DL catalogue remains a family seed, not a Thomas-grade
catalogue.** EQ/blast screen: `data/eq_catalog_{pgc,b011}.csv` hold 44 regular earthquakes in the
vi2010 window; 0.45-1.2 % of our associated events fall within 60 s of one (chance 0.18 %) —
present but negligible.

## STEP 5 — the Thomas harness runs on Bostock's own 2012 catalogue; two metrics reproduce, one does not

`scripts/e4_thomas_harness.py` takes any catalogue with (time, family_id, lat, lon, depth,
magnitude) and runs unchanged on ours later. Thomas's own methods are qualitative — the paper gives
no detection algorithm, binning, completeness cut, decay form or tidal method — so the operational
criteria are OURS and are written to `results/e4_thomas_criteria_bostock2012.json`.
Figure `figures/e4_step5_thomas.png`; summary `results/e4_thomas_summary_bostock2012.csv`.
Geometry: along-strike azimuth **-40.3 deg (N320E)** from the Cascadia mesh; family along-strike
span 112.5 km, down-dip 62.1 km. 28,149 events, 130 families, 2012-08-30..09-25.

| magnitude cut | events kept | families >= 50 ev | **main-front speed (km/day, southward)** | secondary fronts (vs null) | re-rupture within 5 d | Schuster R, 1st / 3rd of the SSE |
|---|---|---|---|---|---|---|
| none | 100 % | 117 | **3.3** (lead edge 3.4) | 0 (0.8) | 0.934 | 0.140 / 0.251 |
| M >= 1.7 | 51.7 % | 72 | **5.0** (5.0) | 2 (4.8) | 0.909 | 0.048 / 0.352 |
| M >= 1.8 | 33.5 % | 55 | **4.9** (5.1) | 3 (5.4) | 0.888 | 0.081 / 0.394 |
| M >= 1.9 | 18.2 % | 40 | **5.1** (5.2) | 0 (5.2) | 0.834 | 0.133 / 0.432 |
| M >= 2.0 | 8.3 % | 13 | **5.2** (5.3) | 0 (0.4) | 0.717 | 0.058 / 0.517 |
| M >= 2.1 | 3.3 % | 2 | **5.3** (5.3) | 0 (0.0) | 0.627 | 0.096 / 0.543 |

* ~~**Main front: REPRODUCED, and it needs almost no completeness.**~~ **WITHDRAWN (E5 step 0, 2026-09-11; results/merlin_audit_e4_2026-09-11.md F5):** the 3.3-5.3 km/day is window-dependent (days 0-27 -3.3; propagation phase days 6-21 -4.5..-5.2) and Thomas's reference is 80 km in 11 d = 7.3 km/day; re-estimated with a first-BURST front in E5 step 4. Original text: The 6-h activity-centroid and
  leading-edge tracks give 3.3-5.3 km/day SOUTH-EASTWARD along strike — Thomas's ~5 km/day, and
  within a factor 1.0-1.5 at every completeness down to **3 % of Bostock's events**. (The
  first-activation regression is the weak estimator: -1.2 km/day on the full catalogue, because with
  28,149 events every family activates on day 1; it converges to -5.6 km/day as Mc rises.)
* ~~**Tidal sensitivity growing through the SSE: REPRODUCED qualitatively, and it survives to 3 %
  completeness**~~ **WITHDRAWN (E5 step 0, 2026-09-11; results/merlin_audit_e4_2026-09-11.md F6):** a clustering artefact -- Schuster assumes independent events, but the 28,149 events sit in 2,047 bursts; the matched-n per-family test gives R_early 0.68 > R_late 0.42 (the opposite). Burst-level R with a time-shift null in E5 step 4. Original text: — the Schuster R of the first third of the episode is 0.05-0.14 against 0.25-0.54
  for the last third at every cut. **The tidal phase is a PLACEHOLDER**: no Cascadia interface
  tidal-stress model and no tide package is on disk (the only tidal code,
  `catalogs/ducellier_catalog/src/`, reads a California-only `.mat`), so the phase comes from a
  degree-2 equilibrium solid-Earth tide computed here — a proxy for the vertical tide, not for shear
  stress on the interface. Nothing quantitative should be taken from it.
* ~~**Secondary fronts: NOT REPRODUCED, and the bottleneck is the ANALYSIS, not the catalogue.**~~ **WITHDRAWN (E5 step 0, 2026-09-11; results/merlin_audit_e4_2026-09-11.md F4):** the detector was gated to 1-10 km behind a main front taken from the first-EVENT Theil-Sen fit (-1.18 km/day, rho -0.12), which never leaves the origin while the real front runs +35 -> -45 km; 0 vs 0.8 is uninformative and the bottleneck was NOT established either way. Redone with a leading-edge-tracked front in E5 step 4. Original text:
  With our operational criteria (>= 3 families bursting in a 6-h window, along-strike span >= 5 km,
  |Spearman rho| >= 0.7 between position and burst-start time, speed 10-40x the main front or
  50-400 km/day absolute, duration 2-6 h, centroid 1-10 km behind the main front, and a rate
  anomaly >= 3x the families' own background) we detect **0-3 per episode, never above the
  circular-shift null (0-5.4)**, on Bostock's *complete* catalogue. Thomas found four by eye on the
  same data. Either his criteria are tighter than the paper states or a different statistic is
  needed; **this must be settled on his catalogue before any download is justified by it.**
* **Re-rupture: robust to ~18 % completeness.** The fraction of bursts re-activating a family within
  5 days falls 0.93 -> 0.83 down to Mc 1.9 and then to 0.63 at Mc 2.1, while the median inter-burst
  interval inflates 1.8 h -> 3.9 h -> 6.5 h. Usable to Mc ~1.9, biased beyond.
* **Per-family rate decay: NOT usable with this estimator.** The Omori p median swings 0.04 -> 0.22
  -> 0.70 -> 0.48 -> 0.63 -> 0.73 across the decimation and only 13-37 of 117 families ever yield a
  fit. Flagged, not reported as a result.

## STEP 6 — the station metadata kills one of the two candidate downloads

`results/e4_station_availability.csv`, `e4_station_epochs_full.csv`. Unrestricted FDSN channel
metadata, EarthScope + NRCan, with a network and coordinate filter (station CODES collide: TWGB and
TWKB also exist in Taiwan's TW, LCBC in Colombia's CM, MGCB in Chile's ZR, OCP in CI and UU, B011 in
SN and ZF — without the filter POLARIS stations are reported as "open").

* **Of the 21 SeisBench template stations, exactly SEVEN were recording 3-C at >= 40 Hz in
  Aug-Sep 2010 and in Aug-Oct 2012: LZB, PFB, PGC, SNB, SOKB, VGZ, YOUB** — that is the outer ring
  we already have on disk, plus SOKB.
* **The 14 inner POLARIS/C8 stations closed between 2004-06 and 2007-02** (KELB 2004-06-03,
  TWGB 2005-05-10, GOWB 2005-08-19, SHVB 2005-10-13, GLBC 2005-10-16, LCBC 2005-10-12,
  TWKB 2006-05-14, SSIB 2006-06-16, TSJB 2006-06-19, SILB 2006-07-04, JRBC 2007-01-15,
  TWBB 2007-02-14, MGCB 2007-02-19, KHVB 2004-08-05). **A POLARIS Aug-Sep 2010 pull is IMPOSSIBLE**
  — those instruments did not exist then. Their data cover ~2002-2007, i.e. Muenchmeyer's 2003-2005
  windows, not ours. **This retires Merlin's first-download option (a) outright**, and it also means
  Bostock's 2010 and 2012 detections were made with essentially the same seven-station geometry we
  are testing on.
* **The southern-VI ring has since eroded**: LZB closed 2014-06-17, YOUB 2014-09-03, SOKB 2015-01-18,
  PFB 2018-12-20, SNB 2021-03-28. Only PGC and VGZ of the 21 template stations are still running.
  The current replacement set is PGC, VGZ, NLLB, CLRS, SYMB, PTRF (CN), B011, B926 (PB) and
  LRIV, OCP, STW (UW, the last two Z-only for LFE purposes; their post-2016 3-C channels are
  strong-motion `EN?`).
* **Volumes (nothing downloaded).** 0.42 MB per channel-hour at 100 Hz, measured from the files on
  disk. (a) the 7 template stations for 2010-08-01..09-30 (61 d): **12.6 GB native, 5.0 GB at
  40 Hz** — and 6 of the 7 are already here. (b) southern VI for 2012-08-15..10-15 (62 d), the 7
  template stations plus B011, B926, NLLB, LRIV = 11 stations: **20.1 GB native, 8.1 GB at 40 Hz.**

## STEP 7 — the prediction table, and the threshold is the fleet decision

`scripts/e4_prediction_table.py`, `results/e4_prediction_table.csv`, `e4_prediction_sanity.md`.
Figure `figures/e4_step7_prediction.png`. Model, all inputs measured in E4: per-station 3-C excess
**2.15 MAD (< 40 km) / 0.83 MAD (40-60 km)** for surface broadband and **1.53 MAD** for boreholes
(a power law in distance does NOT fit these points — at 36-42 km the measured excess runs 0.68
(PFB), 1.15 (VGZ), 1.19 (B011), 2.15 (PGC): site dominates distance over this range, and a fitted
power law failed the sanity check below by a whole magnitude); amplitude scaling 10^(0.38 dM);
network degradation 0.86; q(N=500)/q(published) = 0.77/0.89; k-th nearest station at
`nearest_km x sqrt(k)`.

**Sanity check on Bostock's own 21-station array** (median family distance 31 km): the model
predicts a network excess of **7.0 MAD** at M1.8-2.1, i.e. a completeness magnitude of **M2.10 at an
8-MAD threshold** — above his catalogue's median M1.70. ~~Inverting instead: if his catalogue is
complete to its 10th percentile M1.43, **his operating threshold in these units is ~3.6 MAD** —
which is D3's 3.5 MAD, arrived at from a completely different calculation. So the model is not
broken; **an 8-MAD threshold is roughly one magnitude unit more conservative than Bostock's own
practice**, and the fleet question is what false-alarm rate is affordable, not whether the signal is
there.~~ **WITHDRAWN (E5 step 0, 2026-09-11; results/merlin_audit_e4_2026-09-11.md F2/F3):** the
model compared the EXCESS statistic to the threshold, but detection is on the RAW statistic
(raw = excess + 2.73 background), so every Mc in this step is ~0.48 mag too high at 8 MAD (1.3 at
4 MAD) and the 7.0-MAD / M2.10 sanity number above inherits the same error. Corrected (E5 step 2)
the model reproduces Bostock's modal M~1.7 at 8 MAD with no threshold inversion; the "3.6 MAD =
D3's 3.5" confirmation compared two different quantities (a detection threshold vs a predicted
network signal), and a literal 3.6-MAD threshold would make his catalogue >90 % noise. The table
below is kept for the record only; its 4-MAD columns are FA-PROHIBITED (E5 step 1 quiet-day FA):

| region | year | stations within 60 km | Mc at 8 MAD | Mc at 4 MAD | fraction of Bostock density at 4 MAD | events/family/SSE at 4 MAD (p50) |
|---|---|---|---|---|---|---|
| S. Vanc. Isl. / JdF | 2019 | 11 | 2.62 | 1.84 | 0.27 | 51 |
| S. Vanc. Isl. / JdF | 2023 | 17 | 2.42 | 1.64 | 0.63 | 120 |
| Puget / Olympic | 2019 | 14 | 2.69 | 1.91 | 0.17 | 32 |
| Puget / Olympic | 2023 | 26 | 2.19 | 1.40 | 0.92 | 177 |
| N. Oregon | 2023 | 20 | 2.21 | 1.43 | 0.91 | 174 |
| C. Oregon | 2023 | 12 | 2.46 | 1.68 | 0.54 | 103 |
| N. Vanc. Isl. | 2023 | 9 | 2.76 | 1.98 | 0.09 | 17 |
| NorCal / MTJ | 2025 | 14 | 2.43 | 1.64 | 0.61 | 117 |

~~**Crossed with the decimation leg:** the main front needs >= 3 % of Bostock's density and re-rupture
needs >= 18 %. At a 4-MAD threshold, Puget/Olympic and N. Oregon from 2023 (0.91-0.92) and southern
VI from 2023 (0.63) clear both; at 8 MAD none does. **The threshold, and therefore the network false
alarm rate, is the binding constraint — not the station geometry.** On the 2010 ring the measured
false-alarm rate is 82 per template per day at 6 MAD and 2 per template per day at 8 MAD, with the
shifted control at 0; that gap is what an SNR-weighted network stack and a longer template have to
close.~~ **WITHDRAWN (E5 step 0, 2026-09-11; results/merlin_audit_e4_2026-09-11.md F1):** '82 / 2 per template per day' are TOTAL declustered DETECTIONS (median over family-days), most of them real LFEs on ETS days; the quiet-day false-alarm rate is measured in E5 step 1 (see the E5 section: FA budget table) and the only FA-viable threshold is ~8 MAD, so the '4-MAD-class threshold suffices in three regions' conclusion is VOID. The decimation requirements themselves (main front 3 %, re-rupture 18 %) were computed with the estimators withdrawn above and are re-derived in E5 step 4.

## What E4 changes

1. **E3 step 4 is retracted and reversed.** Template matching with the public stacks DOES register
   Bostock's events on a six-station 2010 network (+2.6 to +3.4 MAD, 12-26 % at 8 MAD, chance 0.02 %).
   "D3 stands" is withdrawn. Back-projection (E2-audit step 6 / E3 prescription 7) is **not**
   triggered.
2. ~~**The project's s/n anchor is wrong by 7-13x.**~~ **WITHDRAWN (E5 step 0, 2026-09-11; results/merlin_audit_e4_2026-09-11.md F7): coherent and incoherent amplitude are two different quantities; the coherent anchor is 1.3-2.5x the 0.078, not 7-13x.** The measured per-event incoherent s/n at PGC is
   0.57-1.02, not 0.078; and the magnitude scaling of Bostock's catalogue is 10^(0.4 dM), not
   10^(1.5 dM). Every sensitivity statement conditioned on 0.078 needs re-deriving.
3. **The network combination must change** from an equal-weight channel mean to an SNR-weighted
   stack, and templates should be 8-12 s, not 2-4 s.
4. **Transfer to boreholes works at ~90 events** (B011 0.85x the co-located surface template);
   transfer to distant or single-channel stations does not.
5. **A POLARIS 2010 download is impossible** (stations closed 2004-2007), which removes one of the
   two candidate downloads outright. ~~The Sept-2012 southern-VI pull is 8.1 GB at 40 Hz for 11
   stations.~~ SUPERSEDED (E5): the 2012 replication set is the user's 15-station set (PLAN E5); volume in the E5 section.
6. ~~**The secondary-front metric is blocked on the analysis, not the data**, and must be settled on
   Bostock's own catalogue first.~~ **WITHDRAWN (E5 step 0, 2026-09-11; results/merlin_audit_e4_2026-09-11.md F4): the E4 detector was disabled by its own front model; neither 'blocked on the analysis' nor its opposite was established. Redone in E5 step 4.**

# E5 (2026-09-11) — the false-alarm budget, the corrected prediction, the weighted network, the fixed Thomas harness

**Agent: Fable 5.1 (the Opus launch died on its quota; identical spec). Nothing downloaded except the PMC text
of Thomas et al. 2025. No E2-E4 output deleted or overwritten; new files are `e5_`-prefixed; the correlator's
new behaviour sits behind `--budget` / `--weights`.** Ledger: `NOTES.md` "2026-09-11 E5 steps 0-4" (every run,
number, bug and code line). Specification: `results/merlin_audit_e4_2026-09-11.md` (F1-F12). Station sets per the
user: the 2012 replication set is PGC, VGZ, NLLB, SNB, LRIV, B001, B003-B007, B009-B011, B926 (15; LZB, YOUB, PFB and
the C8 array excluded); on the 2010 data the OPERATING SUBSET (PGC, VGZ, SNB + transfers B011, B926, LRIV, NLLB) is
primary, the full E4 ring + transfers and the bare E4 six-station ring are reference columns.

## STEP 0 — the record is corrected
The E4 "false alarms per template per day" were TOTAL declustered detections (median over family-days); every
sentence built on them and on the 4-MAD column is struck in RESULTS E4 and NOTES E4 with a one-line reason and a
pointer to the audit (3.6-MAD inversion, "4-MAD-class threshold suffices", "main front reproduced at 3.3 km/day",
"tidal growth reproduced", "secondary fronts not reproduced", "fine-tuning reopened", "11 stations / 8.1 GB");
PLAN GOAL now reads "~7 km/day (80 km in 11 days; the paper quotes ~5-10)".

## STEP 1 — the detection budget: 8 MAD holds on quiet days; precision is a station-count problem
`e4_template_test.py --budget` (every family with a template scored on every day; declustered detections,
rolled control, +-2 s matches to the same / any Bostock family, test-block split, extras near a neighbouring
family), `scripts/e5_budget_table.py`; templates = Bostock's published stacks on the ring + BUILD-half transfers
for B011/B926/LRIV/NLLB (82 families, 4/8/12 s), moveout from the build half. Quiet day = <= 50 Bostock events
(doy 228, 229, 249-251); 19 ETS days. Figure `figures/e5_fa_budget.png`; tables `results/e5_budget_table_*.csv`,
`e5_budget_bynsta_*.csv`, `e5_budget_bynsta_precision_*.csv`.
**Bug fixed first:** the rolled control did not roll each station's gap mask with its series, so a 0.2 % gap at
LRIV on doy 251 gave 50.8 "chance" detections/template/day; fixed (`vmask_s`), all runs re-done; the quiet-day
rolled rate now sits just below the quiet-day data rate (0.64 vs 0.76), as a chance control must.

| network (8-s template) | T (MAD) | quiet-day FA /template/day, mean (p90) | ETS-day total (rolled) | precision own / any family, +-2 s | recall of Bostock events, +-2 s (test half) |
|---|---|---|---|---|---|
| operating subset, 7 stations, 127 families | 7 | 8.5 (18.0) | 11.8 (7.8) | 0.14 / 0.28 | 0.264 (0.206) |
| | **8** | **0.76 (2.0)** | 2.72 (1.09) | **0.37 / 0.54** | **0.154 (0.117)** |
| | 10 | 0.014 (0) | 0.48 (0.07) | 0.65 / 0.78 | 0.048 (0.035) |
| full ring + transfers, 10 stations | 7 | 5.3 (10.0) | 10.6 (5.0) | 0.24 / 0.37 | 0.408 (0.347) |
| | **8** | **0.42 (1.0)** | 3.05 (0.49) | **0.56 / 0.69** | **0.268 (0.220)** |
| | 10 | 0.005 (0) | 0.76 (0.02) | 0.84 / 0.90 | 0.099 (0.080) |
| E4 six-station ring, published stacks | 8 | 0.51 (1.0) | 3.00 (0.93) | 0.42 / 0.60 | 0.195 |
| | 10 | 0.008 (0) | 0.59 (0.05) | 0.73 / 0.85 | 0.067 |

* **Quiet-day FA passes at 8 MAD (0.42-0.76 <= 1) and at 7 MAD (5.3-8.5 <= 10) for 8-s templates**; 12-s templates
  fail the 8-MAD bar (1.26 on the subset, 0.62 on the ring) — F12: length buys recall (0.154 -> 0.192) at the FA
  price, so the FA-matched pairs are 8 s @ 8 MAD or 12 s @ ~9 MAD. 5-6 MAD: 650-700 and 70-96 FA/template/day
  (> 95 % noise) — 4 MAD is FA-prohibited with room.
* **ETS-day precision FAILS the 85 % bar on the own-family definition (37-56 % at 8 MAD)**, but the failure lives in
  the families with 1-3 templates on the network (precision 0-11 %, 2.3 quiet-day FA at one station). Families with
  >= 4 (subset) / >= 5 (ring) stations: **58-75 % / 67-78 % own-family, 66-82 % / 77-86 % any-family at 8 MAD; 90-97 %
  at 10 MAD**. Every 2012 family has >= 7 of the 15 stations within 60 km (median 11).
* **The extras are LFE-related, not noise:** 69 % of the unmatched 8-MAD detections on ETS days fall within 60 s of a
  NEIGHBOURING family's catalogued event vs 38 % for random times (43 vs 38 % at the noise-dominated 5 MAD) —
  sub-threshold LFEs Bostock's seven-station detector did not catalogue and/or cross-family aliases. The quiet-day
  rate (~25 % of the ETS-day detections) is the noise floor.
* Recall by magnitude is FLAT (subset: 0.125 / 0.150 / 0.162 / 0.225 for M<1.5 / 1.5-1.8 / 1.8-2.1 / >= 2.1); by
  stations carrying a template (ring + transfers): 1-4 stations 0.11-0.27, 5-7 0.23-0.27, 8 stations 0.37.
* **Operating threshold: 8 MAD with 8-s templates** (the "quiet FA >= 5 -> 9-10 MAD" signature did not fire; 10 MAD
  is reported alongside as the precision-driven alternative — recall drops to a third). Carried into steps 2-3.

## STEP 2 — the corrected prediction model, and what it says about the 15-station 2012 set
`scripts/e5_prediction_model.py`; `results/e5_prediction_{spread,insample,arrays,2012_families,2012_recall,
regionyear,measured_curve,measured_curve_within}.csv`, `e5_prediction_model.md`; figure `figures/e5_recall_curves.png`.
Fix (F2): raw = excess + 2.73 is compared with T; recall(T, M) = P(raw >= T) with the MEASURED spread of the raw
statistic at fixed magnitude (2010 ring, 8 s, M1.8-2.1: p10/50/90 = 3.39 / 5.89 / 9.51 -> pooled sigma_log 0.175,
within-family 0.140). Per-station registration at published quality: surface < 40 km 2.15 MAD, 40-60 km 0.83,
borehole 2.12 (B011 1.65, B926 2.60 after dividing the ~90-event transfers by the measured q-ratio 0.72); template
quality q(N) calibrated on that ratio (r1_eff 0.088: q(50/100/200) = 0.53/0.66/0.78); 12 s = 1.10 x 8 s; network
equal-weight g (sum e)/sqrt(n) and SNR-weighted g sqrt(sum e^2) with g = 0.86 measured (1/1.38 as the conservative
bracket); magnitude alpha 0.38 (E4 within-family) with 0.2 / 0 / 0.46 / 0.58 as brackets.
* **In-sample PASS** (each 2010 family's actual station set): predicted raw median 5.47 vs measured 5.89 (bar +-0.5);
  frac >= 8 0.165 vs 0.206 (bar +-0.05) with the pooled spread — the within-family spread gives 0.114 because the
  recall tail is between-family registration variance. Excess under-predicted by 1.15x (labelled bracket only).
* **Bostock's 21-station array: Mc(8) = 1.75 with all stations within 60 km (1.98 with the templated stations
  only), inside the 1.6-1.9 target.** The 7-station 2010-13 array (ring + SOKB) predicts Mc(8) = 2.54 (all within
  60 km) / 2.77 (templated) against his 2010 modal M 1.75: **our implementation registers ~1.7-2.0x below Bostock's
  on the same seven stations (+0.8-1.0 mag)** — his 1,400-2,000-event iteratively refined templates and combination
  vs our published-stack pipeline. That deficit is the honest calibration of everything below.
* **Measured recall is magnitude-FLAT when pooled** (0.18-0.25 from M1.35 to M2.25 at 8 MAD; the distant families
  carry the largest catalogued magnitudes) and rises only mildly within family (0.157 -> 0.272 over +-0.35 mag,
  half the alpha = 0.38 slope). The per-bin 2012 numbers are therefore model-shaped; alpha 0.2 is the central case;
  **the pooled overall recall is the number the download should be held to.**
* **Prediction for the 15-station 2012 set** (130 families, per-family n 50/192/392, build half n/2 -> 126 templated;
  stations within 60 km median 11, p10 7; q(build) median 0.66; Mc(8) median 2.25-2.29):

| set | template | weights | alpha | T | recall all | M<1.6 | M1.6-1.8 | M1.8-2.0 | M>=2.0 |
|---|---|---|---|---|---|---|---|---|---|
| 15 stations | 8 s | equal | 0.38 | 8 | 0.267 | 0.175 | 0.259 | 0.340 | 0.414 |
| 15 stations | 8 s | SNR | 0.38 | 8 | **0.279** | 0.184 | 0.273 | 0.355 | 0.429 |
| 15 stations | 8 s | SNR | **0.2** | 8 | **0.321** | 0.259 | 0.323 | 0.370 | 0.395 |
| 15 stations | 8 s | SNR | 0 (flat) | 8 | 0.378 | 0.367 | 0.385 | 0.388 | 0.360 |
| 15 stations | 8 s | SNR, g 1/1.38 | 0.38 | 8 | 0.209 | 0.135 | 0.202 | 0.269 | 0.335 |
| 15 stations | 12 s | equal / SNR | 0.38 | 8 | 0.310 / 0.324 | 0.207 / 0.218 | 0.303 / 0.318 | 0.391 / 0.408 | 0.469 / 0.483 |
| 15 stations | 8 s | SNR | 0.38 | **10** | 0.093 | 0.040 | 0.082 | 0.135 | 0.206 |
| 18 (+LZB/YOUB/PFB) | 8 s | equal / SNR | 0.38 | 8 | 0.298 / 0.316 | 0.193 / 0.206 | 0.289 / 0.307 | 0.380 / 0.401 | 0.477 / 0.496 |
| 18 (+LZB/YOUB/PFB) | 12 s | equal / SNR | 0.38 | 8 | 0.345 / 0.364 | 0.228 / 0.242 | 0.336 / 0.356 | 0.436 / 0.457 | 0.536 / 0.556 |

  **At 8 MAD the 15-station set is predicted to recall 27-32 % of Bostock's 28,149 events (8 s; 31-38 % at 12 s at
  the FA price), 21 % under Merlin's conservative gain, 9 % at 10 MAD.** Merlin's download condition (>= 20 % at
  8 MAD) is met by every 8-MAD row; the pre-registered per-bin bars (>= 60 % at M >= 2.0, >= 35 % at 1.8-2.0,
  >= 20 % at 1.6-1.8) are NOT all predicted to be met (M >= 2.0: 0.40-0.43 at 8 s, 0.48 at 12 s) and rest on a
  magnitude slope the 2010 data only half support. **Cost of excluding LZB/YOUB/PFB: -0.03 to -0.04 overall
  (0.28 -> 0.32 at 8 s), -0.07 at M >= 2.0; families with >= 50 % predicted recall 9 -> 20 (8 s), 30 -> 41 (12 s).**
* Region-year at 8 MAD (`e5_prediction_regionyear.csv`, q at N = 96, alpha 0.38): Mc(8) equal / SNR and the recalled
  fraction of the 2012 magnitude distribution — Puget 2023+ 1.84 / 1.79 (0.43-0.46), N-Oregon 2023+ 1.89 / 1.86
  (0.40-0.42), S-VI 2023 2.04 / 2.02 (0.33-0.34) and 2025 2.03 / 1.97 (0.34-0.37), C-Oregon 2023 2.15 / 2.13 (0.29),
  NorCal 2025 2.15 / 2.10 (0.29-0.31), N-VI 2.36 / 2.33 (0.21-0.22). At 10 MAD every region loses half.
* **Download volume (no fetch): 15 stations x 3 channels x 51 days (2012-08-15..10-05) at 40 Hz = 11.2 GB**
  (F1's measured 1.416 B/sample STEIM2 at 40 Hz; 9.3 GB at E4's 1.17 B/sample); 17.0 GB at the native rates (12
  stations at 100 Hz).

## STEP 3 — the SNR-weighted network, held out: +22-28 % excess, +25-33 % recall at matched false-alarm rate
`e4_template_test.py --weights results/e5_weights_snr.csv --half test` (alternate 2-h blocks, 60-s guard; templates,
moveout AND weights from the BUILD half only; weights = 10^(build-half stack S-SNR dB/20), capped at 3x the median,
median-normalised, 820 (family, station) pairs; the rolled control is recomputed WITH the weights);
`scripts/e5_weighting_table.py`; figure `figures/e5_weighting.png`; families 1 and 15 excluded (they register at
0.5-1.5 MAD excess and 0-3 % recall on every network; reported separately in `e5_weighting_table.csv`).

| network, 8-s template, test half | weights | held-out excess over random (M1.8-2.1 / all) | frac >= 8 MAD, M1.8-2.1 (+-1 s / +-2 s) | rolled >= 8 at events | quiet-day FA at 8 MAD | ETS rolled /template/day | FA-matched threshold | recall at matched FA, M1.8-2.1 |
|---|---|---|---|---|---|---|---|---|
| ring + transfers (10) | equal | 3.34 / 3.29 | 0.242 / 0.246 | 0.000 | 0.42 | 0.51 | 8.0 | 0.246 |
| ring + transfers (10) | **SNR** | **4.07 / 4.07** | **0.341 / 0.349** | 0.000 | 0.68 | 0.74 | **8.27** | **0.308** |
| operating subset (7) | equal | 2.35 / 2.28 | 0.129 / 0.133 | 0.000 | 0.76 | 1.07 | 8.0 | 0.133 |
| operating subset (7) | **SNR** | **2.94 / 2.96** | **0.195 / 0.203** | 0.000 | 1.36 | 1.47 | **8.32** | **0.177** |

* **The weights are not noise:** build-half stack S-SNR predicts the test-half single-station registration with
  Spearman rho = 0.82 over 565 (family, station) pairs (within-family median 0.86) — the failure branch does not fire.
* **Acceptance on the 10-station network: held-out excess 4.07 >= 3.3 PASS (E4 leg C had 2.64), frac >= 8 0.341 >=
  0.22 PASS, rolled >= 8 at the events unchanged at 0.000 PASS — but the quiet-day FA at 8 MAD is WORSE (0.42 -> 0.68;
  0.76 -> 1.36 on the subset) and the rolled RATE rises 40-46 %:** concentrating the sum on the strong channels makes
  its tails heavier at a fixed MAD multiple. At the FA-matched threshold (8.3 MAD) the honest gain is **+25 % recall
  on the ring + transfers (0.246 -> 0.308 at M1.8-2.1; 0.231 -> 0.301 all magnitudes; +27 % at 12 s) and +33 % on the
  operating subset (0.133 -> 0.177)**, and the 10-station weighted network still clears both bars there (excess 4.07,
  recall 0.31). On the 7-station operating subset the bars are not reached (2.94 MAD, 0.195 at 8 MAD; 0.177 at
  matched FA). Merlin's F11 ("right in direction but small, +17 %") is confirmed in size; station count remains the
  bigger lever (7 -> 10 stations: 0.133 -> 0.246 equal-weight).
* Test-block precision at 8 MAD: ring equal 0.53 own / 0.67 any -> weighted 0.47 / 0.62 (same at matched FA within
  noise); at 10 MAD 0.81-0.83 / 0.88-0.90 both.

## STEP 4 — the Thomas harness with fixed estimators: front ~5 km/day, secondary fronts real but the strict null is ~2/episode, no tidal signal, re-rupture in half the families
`scripts/e5_thomas_harness.py`; `results/e5_thomas_summary_bostock2012.csv` (+ per-cut detail files); figure
`figures/e5_thomas_harness.png`. PMC text (fetched): Thomas gives "~5 km/day" (Fig. 2A caption), "80 km in 11 days
or ~10 km/day", secondary fronts 1-10 km from the main front, 10-40x faster, 2-6 h, "four are marked with arrows"
(UNDATED in the text), Fig. 2C "on 14 September ... ~12 km/hour", family 141 "six secondary fronts between 12 and 16
September"; no algorithm, no tidal model. The validation targets are therefore the 14-Sept front and family 141's six.

| Mc (events kept) | first-burst front, km/day (CI) | leading edge, km/day | secondary fronts | circular null p95 | day-shuffle null p95 | 14-Sept anchor | fam-141 fronts 12-16 Sept (core / span) | tide R null-p (all; thirds) | families re-activated (>= 12 h) | interval median (d) |
|---|---|---|---|---|---|---|---|---|---|---|
| none (100 %) | **-4.69 (-5.00..-4.38)**, rho -0.92 | -4.83 | **20** | 4.1 | 3.1 | rho 0.5996 (gate 0.6) | 2 / 5 | 0.07; 0.65 / 0.07 / 0.48 | 49 % | 0.70 |
| 1.7 (52 %) | -5.10 (-5.34..-4.83) | -5.06 | 15 | 9.1 | 4.1 | found, 8.5 km/h | 2 / 4 | 0.055; 0.94 / 0.07 / 0.22 | 35 % | 0.82 |
| 1.8 (33 %) | -5.02 (-5.35..-4.59) | -5.17 | 14 | 14.1 | 4.0 | found, 9.8 km/h | 1 / 1 | 0.03; 0.99 / 0.035 / 0.22 | 33 % | 0.93 |
| 1.9 (18 %) | -5.19 (-5.52..-4.69) | -5.18 | 14 | 18.0 | 5.0 | no | 0 | 0.045; 0.89 / 0.025 / 0.24 | 24 % | 0.94 |
| 2.0 (8 %) | -4.70 (-6.59..-2.25) | -5.54 | 6 | 12.1 | 2.0 | no | 0 | 0.055; 0.94 / 0.41 / 0.17 | 13 % | 0.63 |
| 2.1 (3 %) | -8.67 (-140..-0.2) | -4.19 | 5 | 5.0 | 1.05 | no | 0 | 0.47; 0.87 / 0.83 / 0.24 | 6 % | 1.15 |

* **(a) Main front = first BURST per family inside a leading-edge-derived propagation window (09-07 10:26 ..
  09-21 09:51, 14.0 d; leading edge +22 -> -55.5 km, 73 km travelled): -4.7 km/day (CI -5.0..-4.4, rho -0.92, 98
  families), leading edge -4.8, centroid -4.4; on Thomas's own dates (09-06 -> 09-17) the leading edge moves 50 km
  (full) to 62-67 km (Mc >= 1.7) in 11 d = 4.6-6.1 km/day.** That is the paper's "~5 km/day", at the lower edge of
  the 5-10 accept band and 30 % below the 7.3 km/day implied by "80 km in 11 days". Robust to Mc 2.0 (8 % of the
  events; CI +-0.3 km/day to 18 %), breaks at Mc 2.1.
* **(b) Secondary fronts (leading-edge-tracked front, 3-h / 30-min windows, event-level, >= 8 events / >= 4
  families / >= 5 km, between-family Theil-Sen 1-10 km/h with bootstrap CI excluding 0, |rho| >= 0.6): 20 at full
  completeness against nulls of 2.3 +- 1.5 (circular per-family shift) and 1.8 +- 1.2 (within-day family shuffle).**
  They are compact fronts inside zone-wide background — core extent (events within 5 km of the fit) median 15 km
  (5-30), core fraction 0.81, 3-34 families, speeds 1.5-9.2 km/h = 8-45x the main front, 0.6-17 km behind the
  leading edge — i.e. Thomas-class objects. Validation: six detections fall in 12-16 Sept, two with family 141 in
  their core (five within their span) against his six on that family; **the 14-Sept ~12 km/h front is in the data
  (09:34-11:34 UTC, 7.1-8.5 km/h, 690-830 events from 63-69 families) but sits AT the rho gate at full completeness
  (0.5996 vs 0.6; not tuned) and is detected at Mc 1.7 (8.5 km/h) and 1.8 (9.8 km/h).** So: the detector finds
  Thomas-class fronts well above both nulls, recovers one text anchor at Mc 1.7-1.8 (borderline at full
  completeness) and family 141 partly (2 of 6 by core); the "<= 1 null per episode" bar is not met (the detector's
  own null level is ~2/episode). Not the criteria-mismatch failure branch, not a clean pass. The circular null
  GROWS as the catalogue thins (9 -> 18 at Mc 1.7-1.9), so the count is above the strict null only at >= 52 % of
  Bostock's events (Mc <= 1.7) while it stays above the day-shuffle null at every cut.
* **(c) Tides (PLACEHOLDER equilibrium tide, burst-level, time-shift null of 200 draws; Schuster dropped): no
  significant phase clustering at any stage** — R 0.28 (null p 0.07), thirds p 0.65 / 0.07 / 0.48; the E4
  "growth" is gone. The last third is never significant. No claim; the external dependency is a Cascadia
  interface tidal-stress model.
* **(d) Re-rupture = a burst after >= 12 h of family quiescence: 86 re-activations in 64 of 130 families (49 %),
  4 % of bursts; interval p10/50/90 = 0.53 / 0.70 / 1.84 d, 99 % within 5 d.** The family fraction is completeness-
  sensitive (49 -> 35 -> 24 % at Mc 1.7 / 1.9); the interval median is stable to Mc 1.9.
* **(e) Completeness each metric needs (fraction of Bostock's 2012 events): main-front speed 8 %; secondary-front
  count above the strict null >= 52 % (above the weak null 3 %; the 14-Sept anchor needs 33-52 %); re-rupture
  family fraction >= 52 %, interval median >= 18 %; tides not testable without a stress model.**
* **Crossed with step 2's 15-station prediction (27-32 % recall at 8 MAD, ~35-45 % at M >= 1.8): the 2012 replication
  can reach the main-front speed and the re-rupture interval distribution; the secondary-front count above the
  strict null and the re-rupture family fraction need >= 52 % and are out of reach at 8 MAD (a weighted 12-s
  catalogue at 31-38 % gets to the day-shuffle-null level only); tides are not testable.** Merlin's "what would
  change my mind" clause fires: secondary fronts need >= 35-50 % of Bostock density.

## What E5 changes, and the decisions a human should make
1. **The operating point is 8 MAD with 8-s templates** (quiet-day FA 0.4-0.8 per template per day on the 2010
   networks; 12 s needs ~9 MAD). 4-6 MAD are FA-prohibited (70-700/template/day). The E4 threshold argument is void.
2. **Precision is set by the number of stations carrying a template, not by the threshold:** >= 5-station families
   give 67-78 % own-family / 77-86 % any-family precision at 8 MAD and 90-97 % at 10 MAD; every 2012 family has
   >= 7 of the 15 stations within 60 km. A "report only families with >= 4-5 templates" rule (or a per-station
   coincidence gate) should be part of the step-5 design — a decision, not made here.
3. **The corrected model is in-sample-honest (5.47 vs 5.89 MAD; 0.165 vs 0.206 at 8 MAD) and reproduces Bostock's
   21-station array (Mc(8) 1.75)** but says our pipeline registers ~1.7-2.0x below his on the same seven stations.
4. **15-station 2012 prediction at 8 MAD: 27-32 % overall recall (8 s; 31-38 % at 12 s), 21 % under the conservative
   gain, ~40-45 % at M >= 2.0; 9 % at 10 MAD.** Merlin's >= 20 % download condition is met; the pre-registered
   per-bin bars (60 / 35 / 20 %) are not all predicted to be met and rest on a magnitude slope the 2010 data only half
   support. Excluding LZB/YOUB/PFB costs 3-4 points of recall (7 at M >= 2.0). **Volume: 11.2 GB at 40 Hz** (17.0 GB
   native) for 2012-08-15..10-05.
5. **SNR weighting is worth +25-33 % recall at matched false-alarm rate** and should be in the step-5 matcher, with
   the threshold set on quiet-day FA rather than on a fixed MAD multiple.
6. **Thomas metrics the 2012 replication can reach at 27-32 % recall: the main-front speed (needs 8 %; ours gives
   ~5 km/day on Bostock's catalogue, the paper's own number, 30 % below the 7.3 of "80 km in 11 days") and the
   re-rupture interval distribution (18 %). Out of reach: the secondary-front count above the strict null and the
   re-rupture family fraction (both need >= 52 % of Bostock's events); tides untestable without an interface stress
   model.** Merlin's clause "secondary fronts need >= 35 % of Bostock density -> replan the science before the fleet"
   applies: the 2012 download would validate the detector and the front speed, not the secondary-front science.
7. **Step 5 (the Sept-2012 download) is NOT launched** — it needs the user's explicit go on: the 15-station set at
   11.2 GB; the per-bin pass bars (replace with the pooled >= 25 % or accept the model-slope caveat); the >= 5-station
   reporting rule; whether a 12-s / 9-MAD or an 8-s / 8-MAD operating pair is pre-registered.

# E5 step 5 (2026-09-11) — the 2012 replication: ONE approved download, from-scratch templates, Bostock's 28,149 events as the referee

**Agent: Fable 5.1 (Opus quota exhausted). The only waveform fetch of the session is STEP A below (user "Go" ~06:05 UTC); nothing else was
downloaded. Nothing of E2-E5 was deleted or overwritten; new files are `e5r_`-prefixed plus the catalog `results/e5_catalog_2012_sVI*.parquet`;
`e4_template_test.py` changed only behind new flags. Ledger: NOTES.md "2026-09-11 E5 step 5 (2012 replication)". Pre-registration:
`results/merlin_audit_e4_2026-09-11.md` prescription 5 -- every bar below is reported as met / not met at the pre-registered operating point
(8 MAD, SNR-weighted network, 8-s templates; 12 s alongside); nothing was tuned to a bar. 9 and 10 MAD are DIAGNOSTIC columns.**

## PASS / FAIL — pre-registered bars, Bostock's TEST-half events (13,381; alternate 2-h blocks, 60-s guard), 15 stations, 2012-08-15..10-05

Two recall definitions are reported because the audit's wording admits both: **own** = a detection of the event's OWN family within +-2 s
(like-for-like with the E5 prediction, which modelled the family's own network statistic); **catalog** = a kept detection of ANY family
within +-2 s after cross-family declustering (what "our catalog vs his events" literally asks), with its chance level from the same catalog
shifted by 600 s. The own-family chance (rolled network at the events) is 0.000 at every threshold.

| bar (audit §5) | 8 s / 8 MAD, SNR-weighted (PRE-REGISTERED) | 12 s / 8 MAD | 8 s / **9 MAD** (smallest T meeting the quiet-day bar) | 8 s / 8 MAD equal-weight (reference) |
|---|---|---|---|---|
| recall overall >= 25 % | own **0.234 NOT MET**; catalog 0.346 (chance 0.064) MET | own **0.252 MET**; cat 0.365 (0.074) | own 0.146 NOT MET; cat 0.226 (0.028) NOT MET | own 0.155; cat 0.259 (0.038) |
| recall M >= 2.0 >= 60 % (n = 1,169) | own 0.286 / cat 0.442 **NOT MET** | 0.314 / 0.471 NOT MET | 0.194 / 0.317 NOT MET | 0.197 / 0.324 |
| recall M 1.8-2.0 >= 35 % (n = 3,462) | own 0.234 NOT MET / cat 0.371 MET | 0.259 NOT MET / 0.395 MET | 0.143 / 0.236 NOT MET | 0.159 / 0.280 |
| recall M 1.6-1.8 >= 20 % (n = 4,769) | own 0.236 **MET** / cat 0.338 MET | 0.250 / 0.354 MET | 0.143 NOT MET / 0.216 MET | 0.157 / 0.255 |
| (M < 1.6, n = 3,981; no bar) | 0.218 / 0.307 | 0.230 / 0.323 | 0.139 / 0.208 | 0.138 / 0.220 |
| rolled control <= 5 % of ETS-day (>= 100 Bostock events) detections | **8.3 % NOT MET** | 7.5 % NOT MET | 3.3 % MET | 9.2 % NOT MET |
| pre-SSE quiet days (08-15..08-29) <= 1 detection / template / day | **4.47 (p90 10) NOT MET** (rolled 0.96) | 4.72 NOT MET | 0.63 MET (rolled 0.10) | 1.14 NOT MET (rolled 0.39) |
| extras >= 90 % on >= 100-event days, with lower median amplitude | **38 % NOT MET**; amplitude 0.19 vs 0.65 matched (lower: yes) | 42 % NOT MET; 0.17 vs 0.53 | 46 % NOT MET; 0.24 vs 0.73 | 43 % NOT MET; 0.20 vs 0.56 |
| catalog size / matched to any Bostock event (+-2 s) / to own family | 30,201 / 7,009 (23 %) / 5,100 (17 %) | 33,968 / 7,824 / 5,739 | 8,911 / 4,335 (49 %) / 3,430 (38 %) | 14,400 / 5,160 (36 %) / 3,969 |

**Verdict in one sentence: at the pre-registered point the from-scratch 15-station matcher recovers 23 % (own family) / 35 % (catalog) of
Bostock's test-half events with an own-family chance of 0.000 -- 0.84x the E5 prediction overall and magnitude-FLAT (0.66x at M >= 1.8) --
but every false-alarm bar fails, and the diagnosis below shows why: the pre-SSE "quiet-day" detections are site-coherent noise amplified by
the two co-located station clusters (PGC/B009/B010/B011 and LRIV/B005/B006/B007), which the per-station rolled control cannot see; at the
FA-compliant 9 MAD the recall halves (0.146 / 0.226).** The failure signature is the audit's "chance > 15 % -> FA model wrong" branch (catalog
precision 23 % at 8 MAD), not the "sensitivity short" branch (per-station registration is at the E5 model, below); earthquake/glitch
contamination is negligible (62 of 30,201 detections within 60 s of a catalogued earthquake, chance 0.1 %).

## Measured vs predicted recall (own-family, +-2 s; E5 step-2 model for the actual 15-station set, g 0.86; `results/e5r_score_e5r_snr.json`)

| bin | Bostock test events | measured 8 s / 8 MAD | predicted SNR alpha 0.38 | ratio | predicted SNR alpha 0.2 | ratio | measured 12 s / 8 MAD (pred. 0.38) | ratio | measured 8 s / 9 MAD | ratio vs 0.38 |
|---|---|---|---|---|---|---|---|---|---|---|
| all | 13,381 | **0.234** | 0.279 | **0.84** | 0.321 | 0.73 | 0.252 (0.324) | 0.78 | 0.146 | 0.52 |
| M < 1.6 | 3,981 | 0.218 | 0.184 | 1.18 | 0.259 | 0.84 | 0.230 (0.218) | 1.05 | 0.139 | 0.75 |
| M 1.6-1.8 | 4,769 | 0.236 | 0.273 | 0.87 | 0.323 | 0.73 | 0.250 (0.318) | 0.79 | 0.143 | 0.52 |
| M 1.8-2.0 | 3,462 | 0.234 | 0.355 | 0.66 | 0.370 | 0.63 | 0.259 (0.408) | 0.63 | 0.143 | 0.40 |
| M >= 2.0 | 1,169 | 0.286 | 0.429 | 0.67 | 0.395 | 0.72 | 0.314 (0.483) | 0.65 | 0.194 | 0.45 |

* The pooled recall lands at 0.84x the alpha-0.38 prediction and 0.73x the alpha-0.2 "central" one (E5 said the pooled number is the one to
  hold the download to; 27-32 % predicted, 23.4 % measured at 8 s, 25.2 % at 12 s). The magnitude dependence is even flatter than E5's
  half-slope warning: 0.218 -> 0.286 from M < 1.6 to M >= 2.0 (a factor 1.3 where alpha 0.38 predicts 2.3), so the M >= 2.0 bar (60 %) was
  never within reach of this detector class on this geometry; the E5 prediction over-shoots at M >= 1.8 by 1.5x and under-shoots below M 1.6.
* Recall by the number of stations carrying a usable template (own / catalog): 1-4 stations (92 events) 0.076 / 0.141; 5-7 (867) 0.171 /
  0.278; 8-10 (3,074) 0.274 / 0.392; 11-15 (7,132) 0.240 / 0.349; families without a template (280 events, 14 families) 0 / 0.136.
  Per family (`e5r_score_e5r_snr_recall_family.csv`): own recall p10/50/90 = 0 / 0.168 / 0.419; 5 families above 0.5, 36 at or below 0.05
  (14 of them template-less); Spearman with station count 0.42, with build-half event count 0.45, with median magnitude -0.15 -> **strongly
  family-dependent completeness** (the audit's third signature), driven by template coverage, not by magnitude.
* Per-station registration at the test events (M1.8-2.1, 8 s, excess over the random-time median, `e5r_score_e5r_snr_station_registration.csv`):
  PGC 2.11 (E5 model 2.15), B010 1.87, B009 1.66, B011 1.66 (model 1.65), B926 1.66 (model 2.60), NLLB 1.12 (86 km), LRIV 1.04, VGZ 1.01 (1.15),
  B007 1.01, B004 1.00, B003 0.93, B006 0.93, B005 0.78, SNB 0.69 (0.98), B001 0.63. **The physics per station is at the model** (near
  stations within 5 %, far ones above the > 60 km extrapolation); the shortfall is in the network combination and the family coverage.
* +-5 s diagnostic: own 0.265 / catalog 0.424 (8 s); build-half moveout correction nil (median 0.00 s).

## The false-alarm diagnosis (why the FA bars fail; nothing re-tuned)

1. **Bostock's 2010-2013 catalog has NO events outside the ETS windows** (2012: none before 08-30 or after 09-25; 2011 Aug-Sep only; 2013
   Sep-Oct only), so "Bostock <= 50 events/day" cannot certify a day as LFE-free -- the quiet-day bar rested on an untested assumption.
2. The 7,272 quiet-day detections (8 s, >= 8 MAD) are day-clustered (08-22..25 and 08-27/28 carry 500-1,700 each; 08-21 zero), 68 % fall in
   04-11 UTC (night: lowest noise, hence lowest daily MAD), 0 PNSN tremor detections in the box on any quiet day, 0.2 % within 60 s of a
   catalogued earthquake, and their time clustering (74-90 % with a neighbour within 60 s vs 44-78 % Poisson) is the same as on ETS days.
3. The per-family quiet-day rate follows the weight fraction carried by the PGC/B009/B010/B011 cluster: Spearman 0.52; 9.5 / template / day
   for families with > 70 % cluster weight vs 0.1-0.3 for < 35 % (`NOTES`); the three worst families (bkeys 103, 32, 128) carry 69-83 %.
4. **The one-station-per-site test settles it** (`scripts/e5r_run_onepersite.sh`; same templates, weights, moveout and code; PGC, VGZ, NLLB,
   SNB, LRIV, B001, B003, B004, B926; quiet window): quiet-day detections / rolled control = **1.04 at 8 MAD (0.85 vs 0.81), 1.00 at 7, 1.3 at 9-10**,
   against **4.63 for the full 15-station network (4.47 vs 0.96)**. Co-located sensors see the same site noise; it adds coherently in the network
   sum, while rolling each station's series independently destroys it -> the rolled control under-states the chance level of this geometry by
   ~4.6x on quiet days and ~3x after the SSE (09-26..10-05: 3.2x the rolled rate, 34 PNSN tremor detections in the box in ten days).
   Consequence for the ETS days too: ours / rolled = 6.7 over 09-06..09-24, ours / Bostock = 0.56, matched / Bostock = 0.25.
5. What the extras are, then: at 8 MAD 77 % of the catalog is unmatched; 38 % of the extras fall on >= 100-event days (16-17 of 52 days), their
   median amplitude is 0.19 vs 0.65 for matched events (as the bar requires), 25 % are within 60 s of another family's Bostock event vs 17 %
   for random times (E5's 2010 ring: 69 vs 38 %). A minority are sub-threshold LFEs; the majority is the site-coherent floor.
   **Design consequence (a decision, not made here): co-located stations must be merged at the site level (one 3-C series per site, or a
   site-mean cc) before the network sum, and the chance control must roll SITES, not stations.**

## Thomas metrics — our catalog vs Bostock's, same estimators (`scripts/e5_thomas_harness.py`, crash guard only; `results/e5r_thomas_comparison.csv`)

Our catalog = 8 s / 8 MAD SNR-weighted (30,201 events incl. the quiet-day floor; magnitude proxy uncalibrated -> no Mc cut); reference columns =
Bostock's 28,149 and Bostock randomly decimated to our catalog-level recall (0.346 -> 9,789 events; `e5r_bostock2012_random_decimated.csv`).

| metric | Bostock 2012 | Bostock decimated to 34.6 % | **ours 8 s / 8 MAD** | ours 12 s / 8 MAD | bar | met? |
|---|---|---|---|---|---|---|
| main front, first-burst Theil-Sen (km/day, CI) | -4.69 (-5.00..-4.38), rho -0.92, 98 fam | -4.73 (-5.08..-4.37) | **-4.48 (-4.94..-3.98), rho -0.80, 93 fam** | -4.71 (-5.17..-4.27) | within +-20 % of his | **MET (0.95x; 1.00x at 12 s)** |
| main front, 6-h leading edge (km/day) | -4.83 (-5.17..-4.52); 72.8 km travelled | -4.87 | -0.49 (-1.23..+0.02); 18.7 km | -0.50 | -- | not reproduced (see below) |
| propagation window (from the leading-edge track) | 09-07 10:26 .. 09-21 09:51 (14.0 d) | 09-07 .. 09-21 | 08-30 13:25 .. 09-21 09:09 (21.8 d) | same | | |
| secondary fronts (count; circular null mean / p95; day-shuffle null) | 20; 2.25 / 4.05; 1.75 | 26; 4.35 / 7.05; 2.5 | **13; 0.10 / 1.0; 1.05** | 7; 0.0 / 0.0; 0.7 | count within +-1 | **NOT MET** (13 vs 20) |
| secondary-front speed median (km/h) / core extent (km) / behind edge (km) | 4.5 / 14.6 / 6.3 | 4.3 / 12.7 / 6.0 | 6.4 / 17.7 / 10.3 | 6.3 / 21.1 / 9.1 | | Thomas-class objects |
| 14-Sept ~12 km/h anchor found / fam-141 fronts 12-16 Sept (span / core) | no (rho 0.5996 at gate) / 5 / 2 | no / 6 / 5 | yes (1) / 3 / 1 | no / 2 / 1 | | |
| bursts / re-activations (>= 12 h quiet) | 2,047 / 86 | 761 / 51 | 1,179 / 40 | 1,301 / 26 | | |
| **re-rupture: fraction of families re-activated** | **0.49** | 0.29 | **0.25** | 0.22 | within +-0.10 of his | **NOT MET** (-0.24; within 0.04 of the decimated reference) |
| re-rupture interval p10 / p50 / p90 (d) | 0.53 / 0.70 / 1.84 | 0.53 / 0.86 / 1.87 | 0.53 / 0.73 / 1.72 | 0.53 / 0.72 / 2.12 | | median within 5 % |
| tides (PLACEHOLDER equilibrium tide; burst-level R, time-shift null p) | 0.28, p 0.07 | 0.28, p 0.12 | 0.23, p 0.64 | 0.22, p 0.57 | none | no signal anywhere (as E5) |

* **The main-front speed IS reproduced by the first-burst estimator (-4.48 vs -4.69 km/day, 5 %; -4.71 at 12 s)** -- the burst definition
  (>= 4 events within 10-min gaps) filters the incoherent floor, so the estimator E5 rated as needing only 8 % completeness survives a
  23 %-precision catalog. **The leading-edge estimator does NOT**: our catalog has >= 10 events per 6-h bin from 08-15 on (the site-coherent
  floor spread over the whole 110-km zone), so the 10th-percentile track sits near the southern end from the start, the propagation window
  opens on 08-30 and the "speed" is -0.5 km/day. This is the front estimator's known sensitivity to false alarms, not new physics.
* **Secondary fronts: 13 vs his 20** (12-s catalog 7) -- both far above our own nulls (0.1 / 1.05), with Thomas-class speeds (6.4 km/h median) and
  extents; his randomly decimated catalog gives 26. The count bar (+-1) is not met, exactly as E5's decimation leg predicted for a ~30 %
  catalog ("out of reach at 8 MAD"); with the extra 77 % floor the detector's null falls to ~0 while the count falls, so no claim either way.
* **Re-rupture family fraction 0.25 vs his 0.49** -- not met against his value, but within 0.04 of his own catalog decimated to our recall
  (0.29): the shortfall is completeness (E5: needs >= 52 %), not a different process. The interval distribution reproduces (median 0.73 vs
  0.70 d, p10 0.53 both, 100 % within 5 d), as E5 predicted (needs 18 %).
* Per-family rate decay: not evaluated (E5 flagged the Omori estimator as unusable even on his catalog). Tides: placeholder, no claim.

## The product — `results/e5_catalog_2012_sVI.parquet` (8 s / 8 MAD, 30,201 rows) and `e5_catalog_2012_sVI_L12.parquet` (12 s, 33,968 rows)

| column | meaning |
|---|---|
| `time` | detection time (datetime64[ns], UTC) = the ORIGIN time of the network statistic's peak, i.e. Bostock's per-family time base (DB time + 2012 offset) |
| `family_id` | Bostock's 2012 family (bkey 0-129 of `results/e5r_families_2012.csv`); `sb_family_id` = the SeisBench id it maps to (NaN for the unmatched family; three pairs share one) |
| `lat, lon, depth` | the family hypocentre INHERITED from Bostock's catalog (S-only, family-level; no event-level location) |
| `rel_amp`, `magnitude` | relative amplitude = SNR-weighted mean over contributing stations of cc3_s x sqrt(E3_s / E3_quiet_s) (template-projected amplitude in units of each station's quiet-day window rms); `magnitude` = log10(rel_amp), UNCALIBRATED (no moment scaling, no site/distance correction) |
| `n_stations` | stations with a usable template contributing to the family's network that day (3-15) |
| `network_mad`, `cc_best` | the detection statistic (network sum in MAD units of that family-day) and the raw SNR-weighted network cc |
| `template_len_s` | 8 or 12 |
| `matched_bostock_id`, `matched_bostock_mag`, `matched_same_family` | row index of `results/e5r_gt_bostock2012.parquet` of the nearest Bostock event within +-2 s (NaN if none), its magnitude, and whether it is the same family |
| `doy` | 2012 day of year |

Caveats: (i) 77 % of the 8-MAD rows are unmatched and, per the diagnosis, mostly site-coherent noise -- use `network_mad >= 9` (8,911 rows,
49 % matched) or the `matched_*` columns for anything quantitative; (ii) locations are Bostock's family hypocentres; (iii) the amplitude proxy
is uncalibrated; (iv) cross-family declustering keeps the highest network MAD within +-2 s, so near-duplicate families (bkeys 51/52, 74/75,
104/105) alias onto one another; (v) 2012-08-15..10-05 only.

## Other deliverables
* Download: `results/e5r_download_table.csv` -- 15.69 GB, 774 station-days, 6 no-data (PB.B010 doy 228-233); VGZ/SNB are HH? 100 Hz in 2012.
* Templates: `results/e5r_template_qc_2012.csv`, `e5r_templates_2012.npz`, `e5r_weights_snr_2012.csv` -- 1,144 usable (family, station)
  templates over 116 families (stations per family p10/50/90 = 6/10/14); q-control cc vs the SeisBench stacks 0.78 at PGC/VGZ/SNB.
* Time base: `results/e5r_family_offsets_2012.csv` -- 2012 - 2010 drift median +0.50 s, MAD 0.74 s (63 families).
* Figures: `figures/e5r_download_coverage.png`, `e5r_template_snr_vs_N.png`, `e5r_offsets.png`, `e5r_quietday_fa.png`, `e5r_recall_by_mag.png`,
  `e5r_daily_counts.png`, `e5r_front_tracks.png`, `e5r_secondary_fronts.png`.
* Scores: `results/e5r_score_e5r_{snr,eq}.json` + `_fa/_recall_mag/_recall_nsta/_recall_family/_daily/_extras/_station_registration.csv`;
  Thomas: `results/e5_thomas_summary_{ours2012_L8,ours2012_L12,bostock2012_dec}.csv` (+ per-run detail files), `e5r_thomas_comparison.csv`.

## What E5 step 5 changes, and the decisions a human should make
1. **The template route works on the 2012 geometry at the level E5 predicted for the signal (0.84x pooled; per-station registration at the
   model) but NOT at the level it predicted for the noise: the 2010-calibrated FA model (0.4-0.8 quiet-day detections / template / day at 8 MAD)
   does not transfer to a network with two co-located four-sensor clusters (4.5 / template / day; one-per-site 0.85).** The fix is structural
   (site-level merge before the network sum; site-rolled chance control), not a threshold.
2. At the FA-compliant 9 MAD the catalog is 49 % precise but recalls 15 % (own) / 23 % (catalog) of Bostock's events -- below every recall bar.
   Whether a site-merged 15-station network recovers the 8-MAD recall at the 9-MAD false-alarm rate is the single most valuable next test and
   needs NO download (same cached data; one GPU pass).
3. Recall is magnitude-flat (1.3x from M < 1.6 to M >= 2.0), so "complete above Mc" is the wrong framing for this detector class; completeness
   is FAMILY-dependent (coverage), which is what a Thomas-style analysis must carry as a per-family weight.
4. Thomas metrics: main-front speed reproduced (first-burst, 5 %); secondary-front count and re-rupture family fraction NOT reproduced, at the
   completeness E5 said they need (>= 52 %) -- the re-rupture fraction matches his catalog decimated to our recall within 0.04.
5. Per the pre-registration this is a FAIL on 6 of 8 bars at 8 MAD (own-family definition; 4 of 8 on the catalog definition) and on 7 of 8 at
   9 MAD; the two design lessons (co-location; magnitude-flat, family-dependent completeness) are the return on the 15.7 GB.
