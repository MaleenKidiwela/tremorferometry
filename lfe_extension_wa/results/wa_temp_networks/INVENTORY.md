# Temporary seismic networks near the Washington LFE families, 2012-01 → 2026-09

Written 2026-10-02 by `scripts/wa_temp_inv.py` (read-only queries; nothing deleted). Families = the v1.1
representatives of `data/catalogue/families_v1_2.csv` (NW 132, SW 79, Olympic arm A 75, arm B 43; 329 positions,
`pos_lat/lon`, else `approx_lat/lon`). Distances are epicentral, station to the nearest family.

Files (all in this folder):

- `by_network.csv`: one row per network (the table below).
- `inventory.csv`: one row per station (4,932 stations ≤ 75 km): group used, components, rate, sensor, live days,
  sample results, size estimate, nearest 3 families, status, `select` reason.
- `channels.csv`: every FDSN channel epoch; `samples.csv`: every 10-minute sample request and its answer.
- `selected_for_download.csv`: the 630 stations passed to the download queue (`scripts/wa_dl2.py`, tier T), with rank.
- `sparse_families.csv`: families with < 3 usable stations within 30 km. **There are none**: every family has
  ≥ 3 permanent bridge-node stations within 30 km, so the "nearest station for sparse families" rule added nothing.

## Method

1. FDSN station service, box 45.7–49.7 °N, 125.9–120.9 °W, 2012-01-01 → 2026-09-30, at EarthScope, NCEDC and
   NRCan. Temporary = FDSN code starting with X, Y, Z or a digit, plus TA. NCEDC serves no temporary network here.
   NRCan serves C2, QW and QX in addition to CN/PO/C8 (QW and QX are Q-prefix permanent codes, kept as a check).
   Channel level, then kept stations ≤ 75 km from a family.
2. Data check: one 10-minute sample per station-year (the `data/extension_stations/BRIDGE_PLAN.md` convention),
   vertical channel of the preferred seismometer group (H, then L or P geophone, then N accelerometer). Within 40 km,
   up to 4 more times (10/30/70/90 % of the epoch-year) if the middle is empty. A network group with ≥ 10 empty
   stations and none with data gets one try per station after that.
3. Size: median sample bytes × 144 × components, × live days × (years with data / years sampled).
4. Download rule: status `ok` (public, seismometer, data served) and ≤ 30 km from a family; plus the nearest
   temporary station(s) of a family with < 3 usable stations within 30 km. Usable = permanent bridge-node stations
   live in 2012–2026 + selected temporary stations. No family qualifies.

Limits: a year passes if any sample returns data, so outages shorter than a year are not seen. The FDSN availability
services (EarthScope fdsnws/availability, station `includeavailability`, ph5ws) are retired (HTTP 410 / "retired"),
so these samples are the only data check. Node networks whose dataselect returns 204 may be archived in PH5 only.

## Networks found (≤ 75 km of a family)

Columns: groups and components of the station group used; nearest station (family region); stations ≤ 10 / ≤ 30 km;
sample station-years with data / sampled; size of all stations; selected stations / size.

| net | experiment | span (in 2012–26) | stations | channels used | dmin km (sta, region) / ≤10 / ≤30 | samples with data | est. all (GB) | selected / GB | status |
|---|---|---|---|---|---|---|---|---|---|
| XU | Earthscope integrated investigations of Cascadia subduction (the array the NW families were found on) | 2012-01 → 2012-05 | 20 | SH 50 Hz 3-C (2 EH 100 Hz), CMG-40T-1 | 0.0 (W020, Olympic A) / 20 / 20 | 16/20 | 30 | 16 / 30 | ok 16; 4 empty |
| TA | USArray Transportable Array | 2012-01 → 2016-06 | 7 | BH 40 Hz 3-C, Trillium 240 | 6.5 (D03D, Olympic A) / 3 / 4 | 31/31 | 101 | 4 / 51 | ok |
| X4 | Monitoring active faults for tectonic mapping in WA (CVP, LR, MA lines) | 2016-06 → 2020-05 | 28 | HH 100 Hz 3-C, CMG-3T | 2.9 (LR12, SW) / 7 / 22 | 63/65 | 352 | 22 / 280 | ok |
| XD | iMUSH, Mount St. Helens: 41 broadband + 3,076 1-C nodes | 2014-05 → 2016-12 (nodes 4–7 days, 2014-07/08) | 3,117 | HH 100 Hz 3-C CMG-3T; nodes DPZ 250 Hz GS-11V | 0.2 (node 40098, SW) / 113 / 580 | 2,922/3,199 | 801 | 549 / 134 (6 BB 51, 543 nodes 83) | ok 2,922; 195 empty |
| ZZ | Cascadia ship-to-shore, onshore (CW line) | 2012-06 → 2012-07 | 15 | HH 100 Hz 3-C, L28 | 0.6 (CW12, SW) / 5 / 12 | 15/15 | 21 | 12 / 16 | ok |
| YG | Cascadia ship-to-shore, onshore (CO line) | 2012-07 (13 days) | 6 | EH 100 Hz 3-C, L28 | 9.8 (CO03, SW) / 1 / 3 | 6/6 | 2 | 3 / 0.9 | ok |
| 7D | Cascadia Initiative OBS (FN, J57 sites; seafloor) | 2012-01 → 2014-07 | 40 | HH 125 Hz 3-C Trillium Compact (5 BH 50 Hz); 5 pressure-only | 2.8 (FN05A, NW) / 7 / 26 | 62/63 | 340 | 21 / 200 | ok 34; 5 pressure only; 1 empty |
| X6 | JdF crust hydration OBS (OC1206A; seafloor) | 2012-06 (12 days) | 6 | EL 200 Hz 3-C L28LB | 27.4 (27, NW) / 0 / 2 | 6/6 | 17 | 2 / 5.5 | ok |
| YR | Open-access OBS controlled-source dataset (seafloor) | 2021-06/07 (10 days) | 9 | EL 200 Hz 3-C | 25.9 (S114, NW) / 0 / 1 | 9/9 | 13 | 1 / 0.5 | ok |
| Z4 | Large-N onshore complement to the 2020 Langseth survey (nodes) | 2021-05 → 06 | 26 | DPZ 500 Hz 1-C nodes | 5.7 (8003, NW) / 7 / 9 | 0/26 | – | 0 | metadata only, no data served |
| Z6 | Seattle Basin edge (nodes) | 2019-07 → 08 | 100 | DP 500 Hz 3-C nodes | 10.7 (87, Olympic B) / 0 / 100 | 0/100 | – | 0 | metadata only, no data served |
| 4C | GeoGirls, Mount St. Helens | 2015-08 (1 day) | 5 | DP 250 Hz 1-C | 29.1 (DAN, SW) / 0 / 1 | 0/5 | – | 0 | metadata only, no data served |
| 1D | MSH node array | 2014 | 629 | DPZ 250 Hz nodes | 60.1 / 0 / 0 | 0/629 | – | 0 | > 30 km; no data served |
| 2N, 4K, YI | Rainier 2020, Rainier 2023, MSH 2017 node arrays | 2020, 2023, 2017 (5–6 weeks) | 95, 50, 120 | DP 250/500 Hz 3-C nodes | 33.1, 33.2, 67.3 / 0 / 0 | 0 | – | 0 | > 30 km; no data served |
| 3H, 6N | Rainier 2022, Rainier 2021 node arrays | 2022-08, 2021-07 → 09 | 48, 56 | DP 500 Hz 3-C nodes | 33.2 / 0 / 0 | 48/48, 52/56 | 255, 507 | 0 | > 30 km (data served) |
| 2R, 9J, X5 | MSH 2024, MSH 2021, MSH 2022 node arrays | 2024, 2021, 2022 | 73, 33, 68 | DP 250/500 Hz 3-C nodes | 65–66 / 0 / 0 | all | 603, 337, 499 | 0 | > 30 km (data served) |
| YH | Play Fairway, St. Helens shear zone | 2016-06 → 12 | 20 | HH 100 Hz 3-C KS2000 | 52.0 / 0 / 0 | 20/20 | 233 | 0 | > 30 km |
| Z2 | Sediment-transport geophones | 2024-09 → 2026-04 | 5 | EH 125 Hz 3-C L-22 | 33.5 / 0 / 0 | 10/12 | 24 | 0 | > 30 km |
| QW, QX | Canadian EEW (CNEEWS) and extended hosts (NRCan) | 2020 →, 2023-11 → | 35, 39 | CN/HN accelerometers (Titan, Fortimus, TSMA3) | 5.6 / 5 / 14; 7.5 / 3 / 18 | 0/185; 1/115 | – | 0 | accelerometer only (and not served by NRCan dataselect) |
| 1H | MOCHA magnetotellurics | 2013 → 2018 | 23 | LF/LQ | 2.7 | – | – | 0 | not seismic |
| C2 | Canadian geomagnetic observatory (Victoria) | – | 1 | magnetometer | 9.0 | – | – | 0 | not seismic |
| ZD, ZE | Nisqually River bedload (infrasound, nodes) | 2023 → 2025 | 14, 32 | HD infrasound; GP 1000 Hz | 38.5, 35.6 | – | – | 0 | excluded, restricted (closed at provider) |
| **Z5** | Rainier NP wide-aperture node array 2025 | 2025-06 → 08 | 212 | DP 250 Hz nodes | 18.4 / 0 / 8 | not queried | – | 0 | **excluded, restricted** (not public; the user handles it separately) |

Not found in 2012-01 → 2026-09 within 75 km (searched, no stations returned): XA, XG, XO, XQ, XT, YW, ZG, ZO, ZU, 1E.
(XG, YW, XQ and ZU are pre-2012 arrays.) The UW network's own short deployments (BES, BST, RAD nodes 2017–2020) are
under the permanent code UW; the ≤ 10 km ones are download tier D2 (`scripts/wa_dl2.py`).
Permanent networks in the box that are not part of this inventory: AM (Raspberry Shake), C0, C8, CC, CN, GD, GM, GS,
NP, NV (ONC seafloor), OW, PB, PQ, RE, SB, SY, UO, US, UW.

## Download (queue `scripts/wa_dl2.py`, tier T, after tiers C and D)

Selected 630 stations (31,353 station-days, ≈ 0.80 TB at metadata live days; ≈ 0.72 TB scaled by the sample hit
rate), in this rank order (within a rank, nearest station first):

1. XU 16 stations, 2012-01 → 05 (the discovery array, ties 2012 to the array-era templates) ≈ 30 GB
2. TA D03D, D04D, D04E, A04D, 2012 → 2016 (ETS 2012, 2013, 2014, 2016) ≈ 51 GB
3. X4 22 stations, 2016-06 → 2020-05 (ETS 2017, 2018, 2019; 15 near SW families) ≈ 280 GB
4. XD broadband KRES, MA05, ME02, MG02, MB05, MH01, 2014 → 2016 ≈ 51 GB
5. ZZ 12 + YG 3 (summer 2012, no ETS) ≈ 17 GB
6. 7D 21 OBS (2012 → 2014, seafloor noise) ≈ 200 GB
7. X6 2 + YR 1 OBS (10–12 days each) ≈ 6 GB
8. XD 543 1-C nodes (4–7 days, 2014-07/08, no ETS) ≈ 83 GB

The queue's 2.0 TB cap (C 0.61 + D 1.07 TB first) leaves ≈ 0.32 TB, so ranks 1–2 and most of rank 3 fit; ranks
3 (rest) to 8 and tier D2 (≈ 0.07 TB) wait for a cap decision. Every row stays in `results/wa_dl2/plan.csv`;
raising `--cap-tb` and relaunching `scripts/wa_dl2_queue.sh` resumes exactly where it stopped.

XU.W060 is in the selection: it is the bursty station that took over NW family growth (memory
`noisy-station-hijacks-growth`); downloading it is harmless, keep it out of any growth pool.
