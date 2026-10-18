# F1 — fleet inventory and association-region design for a fleet-wide LFE picker run

**Run** 2026-09-10 22:31-22:35 (resumed after the 21:03 pod crash), `scripts/f1_fleet_inventory.py`,
log `logs/f1_fleet_inventory.log`; verification `scripts/f1_verify.py` -> `logs/f1_verify.log`
(**31/31 cross-table consistency checks PASS**).
**No waveforms were downloaded.** Every availability number comes from FDSN *station/channel
metadata* (one `level=channel` request per network, 14 networks, cached as StationXML) plus the
Slab2 mesh and the catalogs already on disk.

> **Read every availability and cost number as provisional and as an UPPER BOUND.** Metadata says a
> channel *existed*; it does not say the day is usable. This project's own measurements on the same
> stations: 12.9 % of PGC days are component-short, 58 % of WISH and 27 % of YACT re-detection
> windows were zero-amplitude, and glitchy station-years can lose ~90 % of days to QC
> (NOTES 2026-09-10 D3; memory `glitch-days-false-dvv-drop`). A data-quality pass would only ever
> *reduce* these counts.

Outputs: `results/f1_{stations,station_years,cell_density,regions,region_stations,cost_estimate,groundtruth}.csv`,
figure `figures/f1_fleet_inventory.png`.

---

## (i) Station table

**Master table** `f1_stations.csv` — 691 stations: **278 "project"** (every station this project has
ever used: `data/broadband_fleet_order.csv`, `data/candidate_stations_post2020.csv`, the 3 anchors)
and **413 "extra"** = non-project FDSN stations that (a) have a 3-C >=20 Hz triplet at some point and
(b) lie within 150 km of at least one ETS source cell. All 278 project stations returned FDSN
metadata (0 missing).

| class (project roster) | N | | class (extra roster) | N |
|---|---|---|---|---|
| broadband (BH/HH, surface) | 181 | | broadband | 387 |
| short-period (EH/SH, surface) | 59 | | broadband_borehole | 4 |
| borehole short-period (PB, depth >=20 m) | 36 | | borehole short-period | 13 |
| broadband_borehole | 2 | | short-period | 9 |
| **total** | **278** | | **total** | **413** |

Networks — project: UW 92, NC 48, PB 38, UO 34, BK 24, CN 23, CC 18, US 1.
Extra: UW 129, CC 68, BK 41, UO 39, TA 34, PO 25, NC 21, CN 20, C8 14, PB 9, OO 9, NV 3, IU 1.
Roles: 237 fleet-broadband, 38 borehole, 3 anchors. Of the 278, **190 already have a dv/v product,
187 a causality certificate, 185 an inversion pair** (`has_dvv`, `has_cert`, `in_inversion`).

**57 project stations never have a 3-C triplet at >=20 Hz** (32 NC + 25 UW; all of them are
`z_only_ever` = a vertical-only short-period record). They are unusable for LFEDetect (ZNE input,
20 Hz) although they carry Z-only data the dv/v pipeline used. List in the run log.

**3-C >=20 Hz availability per year** (`f1_station_years.csv`; one row per station-year, 691 x 18 =
12,438 rows; `days` = exact channel-epoch overlap with the calendar year; **2026 is truncated at
2026-09-10**, hence 252 days):

| year | project stations 3-C | eff. station-years | ZNE | Z12 | borehole | extra stations 3-C |
|---|---|---|---|---|---|---|
| 2009 | 70 | 66.2 | 32 | 38 | 37 | 102 |
| 2010 | 73 | 67.8 | 35 | 38 | 37 | 122 |
| 2012 | 77 | 74.2 | 38 | 39 | 37 | 130 |
| 2014 | 81 | 77.3 | 42 | 39 | 37 | 133 |
| 2016 | 94 | 79.9 | 55 | 39 | 37 | 140 |
| 2018 | 116 | 100.5 | 77 | 39 | 37 | 117 |
| 2019 | 140 | 120.8 | 101 | 39 | 37 | 148 |
| 2020 | 172 | 139.8 | 133 | 39 | 37 | 131 |
| 2021 | 216 | 183.8 | 177 | 39 | 37 | 127 |
| 2023 | 212 | 201.5 | 173 | 39 | 37 | 218 |
| 2025 | 209 | 198.7 | 170 | 39 | 37 | 279 |
| 2026 | 208 | 137.5 | 169 | 39 | 38 | 287 |

Full per-year table in `logs/f1_verify.log` §(i). Totals: **2,167.6 effective station-years**
(project) and **4,603.1** (project + extra). Orientation codes across all 3-C station-years:
ZNE 4,425, **Z12 893** (the PB boreholes and a few others: horizontals of unknown azimuth — fine for
a detector and for S-wave arrival times, a caveat for anything using polarisation).
Sample rates of the fleet's 3-C station-years: **100 Hz 1,882 / 40 Hz 374 / 50 Hz 181** (band mix
HH 1,182, EH 700, BH 555) — the fleet is predominantly 100-Hz HH, which sets the download volume.
The jump from 70 (2009) to ~210 (2021+) project stations is real station installation, dominated by
the post-2020 candidate list; those candidates have *not* been processed by the dv/v pipeline.

---

## (ii) Source-cell station density

Source grid: 0.25 deg cells whose Slab2 interface depth is **30-45 km** (the ETS/LFE band),
40.0-50.0 N -> **137 cells** (R1 22, R2 17, R3 23, R4 22, R5 27, R6 11, R7 15). Per cell-year and per
roster, `f1_cell_density.csv` gives station counts at **60 / 100 / 150 km**, the broadband-only
counts (`n60_bb`, `n100_bb`), nearest station, azimuthal gap at both 60 and 100 km, and the borehole
fraction.

**Median** number of project-fleet 3-C stations per cell (`f1_regions.csv` `med_n60/100/150`;
the per-region mean differs from the median by up to ~2.5 stations in either direction — both are
printed in `logs/f1_verify.log` §(ii)):

| region | n60 2009 | n60 2018 | n60 2024 | n100 2009 | n100 2018 | n100 2024 | n150 2009 | n150 2024 |
|---|---|---|---|---|---|---|---|---|
| R1 N Vancouver Island | 4.0 | 7.0 | 7.5 | 9.0 | 14.0 | 18.0 | 17.5 | 31.5 |
| R2 S Vancouver Island / JdF | 7.0 | 11.0 | 14.0 | 15.0 | 21.0 | 32.0 | 22.0 | 57.0 |
| R3 Puget / Olympic | 5.0 | 5.0 | 16.0 | 13.0 | 18.0 | 42.0 | 23.0 | 82.0 |
| R4 N Oregon | 4.0 | 6.0 | 13.0 | 10.0 | 16.0 | 31.5 | 16.0 | 57.5 |
| R5 C Oregon | 3.0 | 5.0 | 9.0 | 6.0 | 12.0 | 19.0 | 11.0 | 33.0 |
| R6 S Oregon | 3.0 | 5.0 | 8.0 | 6.0 | 11.0 | 19.0 | 9.0 | 30.0 |
| R7 NorCal / MTJ | 3.0 | 5.0 | 7.0 | 9.0 | 13.0 | 19.0 | 15.0 | 39.0 |

(The complete 18-year tables for all three radii, median *and* mean, are in `logs/f1_verify.log`
§(ii). Maps for 2012 / 2018 / 2024 are panels b-d of the figure.)

**The single most important result of F1:** at a **100 km** radius, *every* cell has >=4 3-C stations
in *every* year for *both* rosters (`frac_ge3` = `frac_ge4` = 1.00 throughout, 252/252 region-years).
Station *count within 100 km* is therefore not a discriminating quantity anywhere in Cascadia —
which also means the PLAN's "locatable fraction (>=3 / >=4 / >=6 stations)" as literally specified is
saturated and uninformative. The tables below are therefore reported at **60 km**, the radius at
which an LFE is plausibly above the noise (D3 measured a median Lin/Bostock LFE at PGC, ~40-50 km
from its family, at amplitude s/n = 0.078), *and* split by instrument class.
**60 km is a geometry proxy, not a detectability measurement** — F1 cannot measure detectability from
metadata.

---

## (iii) Association regions and locatable fraction

Seven along-strike regions (lat bands), `f1_regions.csv` (region x year x roster) and
`f1_region_stations.csv` (membership = within 150 km of any cell of the region; regions therefore
overlap and member counts are **not** additive).

Fraction of a region's cells with **>=4 3-C stations within 60 km** (project fleet):

| year | R1 | R2 | R3 | R4 | R5 | R6 | R7 |
|---|---|---|---|---|---|---|---|
| 2009 | 0.68 | 0.82 | 0.87 | 0.55 | 0.19 | 0.36 | 0.33 |
| 2012 | 0.68 | 0.82 | 0.87 | 0.55 | 0.59 | 0.36 | 0.33 |
| 2016 | 0.73 | 1.00 | 0.87 | 0.86 | 0.85 | 1.00 | 0.40 |
| 2018 | 0.95 | 1.00 | 0.87 | 0.91 | 0.89 | 1.00 | 0.87 |
| 2021 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| 2026 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.93 |

Fraction with **>=6 within 60 km** (the E1 "replication-grade" bar) and with **>=4 BROADBAND within
60 km**:

| criterion / year | R1 | R2 | R3 | R4 | R5 | R6 | R7 |
|---|---|---|---|---|---|---|---|
| >=6 <=60 km, 2012 | 0.23 | 0.76 | 0.43 | 0.18 | 0.07 | 0.00 | 0.07 |
| >=6 <=60 km, 2018 | 0.64 | 1.00 | 0.43 | 0.68 | 0.44 | 0.27 | 0.33 |
| >=6 <=60 km, 2024 | 0.95 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.93 |
| >=4 broadband <=60 km, 2012 | 0.00 | 0.24 | 0.04 | 0.00 | 0.00 | 0.00 | 0.13 |
| >=4 broadband <=60 km, 2018 | 0.77 | 0.76 | 0.04 | 0.73 | 0.33 | 0.18 | 0.80 |
| >=4 broadband <=60 km, 2024 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| >=3 <=60 km, 2009 | 0.95 | 0.94 | 1.00 | 0.86 | 0.52 | 0.73 | 0.67 |

Region summary, median over the dense era **2019-2026**, project fleet (all columns in
`f1_regions.csv`):

| region | cells | member stations (fleet / of which broadband) | med n60 | med n60 broadband | med n100 | med nearest km | med az gap <=60 km | med az gap <=100 km | borehole frac <=60 km | frac >=4 @60 | frac >=6 @60 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| R1 N Vancouver Island (48.7-50.0 N) | 22 | 64 / 45 | 7.5 | 6.0 | 18.0 | 18.5 | 118 | 66 | 0.30 | 1.00 | 0.95 |
| R2 S Vancouver Island / JdF (47.9-48.7 N) | 17 | 94 / 65 | 14.0 | 8.0 | 32.0 | 13.0 | 100 | 52 | 0.38 | 1.00 | 1.00 |
| R3 Puget / Olympic (46.6-47.9 N) | 23 | 139 / 97 | 16.0 | 12.0 | 42.0 | 14.1 | 89 | 47 | 0.21 | 1.00 | 1.00 |
| R4 N Oregon (45.2-46.6 N) | 22 | 109 / 80 | 13.0 | 11.0 | 31.5 | 12.1 | 81 | 53 | 0.12 | 1.00 | 1.00 |
| R5 C Oregon (43.6-45.2 N) | 27 | 75 / 51 | 9.0 | 7.0 | 18.5 | 15.7 | 98 | 80 | 0.22 | 1.00 | 1.00 |
| R6 S Oregon (42.0-43.6 N) | 11 | 61 / 39 | 7.5 | 5.5 | 18.5 | 11.1 | 107 | 65 | 0.24 | 1.00 | 1.00 |
| R7 NorCal / MTJ (40.0-42.0 N) | 15 | 93 / 53 | 7.0 | 7.0 | 19.0 | 25.6 | 134 | 70 | 0.11 | 0.97 | 0.93 |

Adding every non-project FDSN station in the region ("all" roster) lifts med n60 substantially
(R3 16 -> 25, R4 13 -> 18.5, R2 14 -> 16.5) and closes the R1 gap (118 -> 104 deg), and pushes the
>=6-at-60-km bar back in time by up to a decade in the north (R2 satisfies it from 2009, R3 from
2010) — the 413 extra stations are the cheapest available improvement to association geometry.

**First year each criterion is met, project fleet** (`logs/f1_verify.log` §(vi)):

| region | >=4 @60 km in *all* cells | >=6 @60 km in >=90 % | >=4 broadband @60 km in >=90 % | >=6 @60 km in >=90 %, all-FDSN roster |
|---|---|---|---|---|
| R1 | 2020 | 2020 | 2019 | 2018 |
| R2 | 2016 | 2017 | 2019 | **2009** |
| R3 | 2019 | 2020 | 2021 | **2010** |
| R4 | 2019 | 2020 | 2020 | 2019 |
| R5 | 2019 | 2020 | 2020 | 2020 |
| R6 | 2016 | 2019 | 2020 | 2016 |
| R7 | 2021 | 2021 | 2019 | 2021 |

The azimuthal-gap condition costs little in the dense era but matters in R1 and R7:
`frac_ge4_60_gap180` (>=4 stations within 60 km *and* max azimuthal gap <=180 deg) is 0.91 in R1 and
0.87-1.00 in R7 for 2021+, against 1.00 elsewhere; before 2016 it falls to 0.00-0.33 in R5/R6/R7.

---

## (iv) Data volume and compute for a fleet-wide picker run

`f1_cost_estimate.csv` (per year x scope, plus TOTAL/ASSOC/LOCATE rows). Measured inputs, not guesses:

* **bytes/sample** from the study folder's own day files, sample counts read from the miniSEED record
  headers (600 files): **1.416 B/sample at 40 Hz, 0.932 B/sample at 100 Hz** (STEIM2 encodes
  differences, so the rate matters); the cost table evaluates a fitted b(fs) = c0 + c1/fs per
  station-year. b(20 Hz) = 2.22 B/sample is an **extrapolation** — treat the decimated-store figure
  as indicative and the float32 figure (4 B/sample exactly) as the safe bound.
* **GPU throughput** from the E2 pilot's own runs, `results/e2_throughput.csv`, snapshot at
  22:31 = 38 runs / 1,962 station-days on the **L40S**: **52.6 station-days/min aggregated per queue**
  (per-run range 32-100), i.e. **6.93 GPU-min per station-year**. E2 ran **4 concurrent queues x 2
  torch threads**, so this per-queue rate already includes 4-way device sharing and the
  `gpu_pick_hours` column is **queue-hours**: divide by the number of concurrent queues for
  wall-clock. *`e2_throughput.csv` is still being appended to by the concurrently running E2 pilot,
  so this figure will drift; the table's `note` column records the snapshot it used.*
* **pick rate** measured over the same runs: median **1,867 S picks per station-day** (raw LFEDetect
  output at the E2 threshold — overwhelmingly false picks; association is what filters them).

| scope | eff. station-years | download, native bands | download, lowest adequate band | 20 Hz store (STEIM2, extrap.) | 20 Hz store (float32) | GPU picker | wall-clock, 1 L40S, 4 queues | S picks |
|---|---|---|---|---|---|---|---|---|
| **project fleet** (221 stations) | 2,167.6 | **17.6 TB** | 17.5 TB | 9.1 TB | 16.4 TB | **250 queue-hours** | **~63 h** | 1.48e9 |
| project + extra FDSN (634 stations) | 4,603.1 | **36.8 TB** | 36.4 TB | 19.4 TB | 34.8 TB | 532 queue-hours | ~133 h | 3.14e9 |

Association and location (order of magnitude, arithmetic in the CSV `note` column; the rates
themselves are **assumptions**, not measurements): PyOcto at 2e5 picks/CPU-min -> **~123 CPU-hours**
(fleet) / ~261 (fleet+extra), embarrassingly parallel by region x day; NonLinLoc at 30 ms/event on a
precomputed travel-time grid, assuming 1 associated event per 1,500 raw picks -> ~1.0e6 events and
**~8 CPU-hours** (fleet) / ~17 (fleet+extra).

Consequences worth deciding on:

1. **The download dominates, not the GPU.** ~63 wall-clock GPU-hours against 17.6 TB over the wire.
   At a sustained aggregate 10 / 25 / 50 MB/s the fleet download is ~15 / ~6 / ~3 days
   (arithmetic on the 17.6 TB; **no download rate was measured** — the E2/D3 download logs carry no
   timing). Streaming per station (download -> decimate to 20 Hz -> pick -> delete) is what the D3
   plan assumed and remains the right shape.
2. **The 100-Hz HH fleet is what costs 17.6 TB.** Only 20 station-years of the fleet have a cheaper
   band that still yields 3-C >=20 Hz (`fs_lo < fs`), so requesting the lowest adequate band saves
   only ~0.3 % — the volume is intrinsic.
3. **Disk is not the binding constraint it was assumed to be:** `/home/jovyan` has ~75 TB free, so
   even 9-16 TB of decimated 20 Hz data is storable; the 80 GB cap that shaped D3 was a
   self-imposed study-folder budget. Whether to spend shared-filesystem space is a human decision.

---

## (v) Ground truth and fine-tuning coverage per region

`f1_groundtruth.csv`; provenance and licences in `results/e1_lin_resources.md`. All of these were
already local (E1); nothing was re-downloaded for F1.

| region | Lin 2024 events (all / after the paper's own residual<0.5 & N>=4 cut) | Bostock 2015 events / families | SeisBench LFE stacks families (**with P+S picks**) | Ducellier southern-Cascadia families | Shelly 2025 MTJ families | PNSN tremor |
|---|---|---|---|---|---|---|
| R1 N Vancouver Island | 675,615 / 66,539 | 31,624 / 24 | 24 | 0 | 0 | 107,762 |
| R2 S Vancouver Island / JdF | 382,499 / 37,532 | 237,799 / 98 | 105 | 0 | 0 | 90,458 |
| R3 Puget / Olympic | 0 | 0 | 0 | 0 | 0 | 74,940 |
| R4 N Oregon | 0 | 0 | 0 | 0 | 0 | 32,342 |
| R5 C Oregon | 0 | 0 | 0 | 0 | 0 | 45,046 |
| R6 S Oregon | 0 | 0 | 0 | 0 | 0 | 143,374 |
| R7 NorCal / MTJ | 0 | 0 | 0 | 34 | 27 | 239,647 |

Detail and caveats:

* **Lin 2024** (Zenodo 10016020): 1,058,114 events, but they occupy **only lat 48.07-49.32 N,
  lon -124.57..-122.94** — i.e. R1+R2 and nothing else. Its own quality cut (residual < 0.5 s **and**
  N >= 4 stations) keeps **104,071 events = 9.8 %**. No per-station picks were released.
* **Bostock 2015** (Slow Earthquake Database, **no redistribution**): 269,423 events, lat
  48.12-48.93 N; 122 distinct (lat,lon) hypocentres against 130 published families (some families
  share rounded coordinates). Origin times only, no picks, no family label.
* **SeisBench `LFEStacksCascadiaBostock2015`** (CC-BY-4.0, on disk): 1,817 three-component 20 Hz
  stacks, **129 families at 21 stations, with `trace_P_arrival_sample` and `trace_S_arrival_sample`**
  — the only per-station-pick substrate we have, and the fine-tuning material. Split train 1,071 /
  dev 184 / test 562. Entirely R1+R2.
* **Ducellier & Creager 2022** (`catalogs/ducellier_catalog/`, local, 2.5 GB): 37 located families
  from the Plourde 2015 template list, lat 39.45-41.72 N -> 34 inside R7, 3 south of 40 N (outside
  the source grid). ~478k LFEs, 66 templates.
* **Shelly 2025 MTJ** (USGS, CC0): 61,441 detections / 27 families 2018-2024. **The release carries
  no coordinates**, so its region assignment (R7) is taken from the paper (Mendocino/Gorda), not
  from data. Note that `catalogs/mtj_ducellier/LFEcat_MTJ_2018-2024.csv` in the main tree is a saved
  404 HTML page, not data (E1; left in place).
* **PNSN tremor** (748,002 detections, all regions): a weak validator everywhere and **not** LFE
  ground truth — tremor epicentres, and our own candidate windows are already PNSN-driven, so it is
  not independent of the existing pipeline.

**R3, R4, R5, R6 have no LFE ground truth of any kind** — not one located family, not one pick.
That is 83 of the 137 source cells (61 %).

---

## (vi) Recommendation: which regions can support an event catalog with hypocentres

Plain language, and deliberately conservative. Two published warnings govern this, both from E1 —
cited, not re-derived:

* **Lin et al. 2024 §4.3**: applied to PB stations from central to southern Cascadia, their picker's
  peaks align with SSE timing only for the northern B001-B013, and "the overall results are noisy and
  may include many false detections"; the 95 % precision is conditional on Vancouver-Island-like
  noise and paths, and they recommend **transfer learning** elsewhere. Their own catalog collapses in
  early 2014 when only three stations (NLLB, PFB, PGC) were available — three stations is the floor
  at which >=3-station association barely functions.
* **Münchmeyer et al. 2024**: template matching "allows ... to identify LFEs with fewer stations than
  deep learning"; their DL Cascadia catalog was "substantially smaller" than the reference catalog;
  and they had to **exclude San Andreas entirely** because a small-aperture borehole array gave
  unusable locations — a direct warning wherever the local geometry is PB-borehole-dominated.

Given those, and F1's geometry:

**Event catalog with hypocentres — go.**
* **R2 (S Vancouver Island / Juan de Fuca)** — the only region that is unambiguously ready across the
  whole 2009-2026 span: >=6 3-C stations within 60 km for 76 % of cells already in 2012 (100 % from
  2018; 100 % from 2009 if non-project stations are included), the smallest azimuthal gap (52 deg at
  100 km, 100 deg at 60 km), the nearest median station (13 km) — and it holds essentially all the
  ground truth (Bostock 98 families, 105 of the 129 pick-bearing stacks, 382k Lin events). This is
  the replication region: validate here or nowhere.
* **R1 (N Vancouver Island)** — go from ~2018-2020 (>=6 at 60 km reaches 0.64 in 2018, 0.95 in 2020).
  Ground truth exists but is thinner (24 families, 24 stacks) and its azimuthal gap at 60 km is the
  second worst (118 deg), so expect worse depth control than R2. Adding non-project stations moves
  its readiness to 2018.

**Event catalog with hypocentres — conditional, fine-tune and validate first (no local ground truth).**
* **R3 (Puget / Olympic)** and **R4 (N Oregon)** — geometry is *excellent* from 2019-2021 (R3 reaches
  med n60 = 16, gap 89 deg; R4 med n60 = 13, gap 81 deg, the best gaps in the margin) and both are
  broadband-dominated, so Münchmeyer's borehole-aperture problem does not apply. But there is **zero
  local LFE ground truth**, so precision cannot be measured in-region; Lin's own transfer-learning
  warning applies directly. Route: fine-tune on the 129 pick-bearing stacks + our hard negatives,
  then validate against borehole/surface consistency and tremor timing, and publish hypocentres only
  if a held-out check in R2 survives the same fine-tuning. Before ~2019 R3's `>=6 at 60 km` is only
  0.43 and R4's 0.18-0.68 — so any R3/R4 event catalog is **network-limited before 2019** even in the
  best case.

**Activity time series only (plus template matching where families exist).**
* **R5 (C Oregon)** and **R6 (S Oregon)** — no ground truth, and the network only crosses the >=6-at-
  60-km bar in 2019-2020 (0.07 and 0.00 respectively in 2012). A 2009-2026 event catalog here would
  have a completeness that changes by an order of magnitude mid-record, which is exactly the artifact
  that makes rate/scaling analyses meaningless. Report a per-station activity time series with an
  explicit completeness curve.
* **R7 (NorCal / MTJ)** — the worst geometry in the margin (median azimuthal gap 134 deg at 60 km,
  nearest station 25.6 km, `>=6 at 60 km` only 0.93 even in 2024 and 0.07 in 2012), *but* the only
  southern region with public family-level ground truth (Shelly 27 families CC0, Ducellier 34 located
  families local). Per Münchmeyer, **template matching against those families is the better tool
  here than the DL picker**; the DL product for R7 should be an activity time series, with located
  events only for the Shelly/Ducellier families themselves.

**Cross-cutting recommendation.** A single margin-wide, temporally uniform event catalog is not
supportable: the fleet's own geometry crosses the replication bar in different years in each region
(2009-2016 in R2, 2019-2021 elsewhere). The honest fleet product is (a) an **event catalog with
hypocentres for R1+R2 over the full span, and for R3+R4 from 2019 subject to fine-tuning
validation**, plus (b) a **margin-wide activity time series** with the per-region-year completeness
proxy from `f1_regions.csv` attached to it. Including the 413 non-project FDSN stations is the
cheapest single improvement available (it advances R2/R3 readiness by ~8-10 years) and costs 2.1x the
download and 2.1x the GPU time.

---

## Caveats (all of them)

1. **Metadata-only.** No gap, noise, glitch or response check. Availability is an upper bound (see
   the box at the top).
2. **"Locatable" is a geometry proxy.** `frac_ge*` counts stations inside a radius; it does not model
   detection at all. The 60 km radius is a design choice anchored on D3's s/n = 0.078 at ~40-50 km
   and on Lin's station-family distances, not a measured detection threshold. At 100 km the criterion
   is saturated (always 1.00) and therefore useless.
3. **Instrument class is per station, not per station-year** (`klass` from the whole channel record),
   so a station counted as "broadband" in `n60_bb` may have contributed a non-BH/HH triplet in a
   particular year.
4. **Azimuthal gap** is computed from epicentral azimuths of stations inside the radius only, ignores
   depth and ignores which stations actually produce a usable pick.
5. **Region membership uses a 150 km radius**, so the seven regions overlap: member-station counts
   across regions are not additive (1,505 region-station rows over 691 stations).
6. **2026 is truncated at the run date** (2026-09-10, 252 days) for open-ended channel epochs; every
   other year is a full calendar year.
7. **b(20 Hz) = 2.22 B/sample is extrapolated** from the 40 and 100 Hz measurements; the float32
   figure is exact arithmetic and bounds it.
8. **The GPU rate is a live snapshot** of a file the E2 pilot is still writing (38 runs at 22:31) and
   was measured with 4 queues sharing one L40S; it is a per-queue rate, so `gpu_pick_hours` are
   queue-hours.
9. **Association/location CPU costs are assumptions** (PyOcto picks/CPU-min, NonLinLoc ms/event,
   picks-per-event), unlike the download and GPU figures which are measured.
10. **No download rate was measured**, so the wall-clock download estimate is parametric only.
11. The **source grid is restricted to Slab2 30-45 km depth**; LFEs outside that band (shallower
    up-dip or in the Gorda plate) are not represented by any cell.
12. **Shelly's region assignment comes from the paper**, not from coordinates (the release has none);
    **Bostock gives 122 distinct hypocentres for 130 published families**.
13. Lin's catalog is used as epicentres/centroids only: 34 % of its depths sit in the deepest grid
    bin and its single-event location scatter is ~22 km (E1).
