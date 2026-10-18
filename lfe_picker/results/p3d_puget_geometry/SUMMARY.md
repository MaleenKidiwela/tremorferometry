# T-study: hypocentral network density over the Cascadia LFE source zones (read-only geometry)

`t_` prefix; everything written here or as `figures/t_*.png`. **No waveforms were downloaded.**
Station metadata comes from the IRIS/EarthScope FDSN **station** service (`level=channel`) and from
files already in the repo. Two passes are recorded: **Pass A** (the original brief — project roster,
Puget vs Oregon vs North) and **Pass B** (the revised brief — *every* station 2019+ in the northern
box, with a ranked list of stations that would help discovery). Pass B supersedes Pass A wherever
they overlap; Pass A still carries the Oregon answer, which is outside Pass B's box.

## Method (both passes)

* **Distances are HYPOCENTRAL**: `sqrt(epicentral^2 + source depth^2)`. A 30-45 km deep source is
  never closer than its own depth, so "sites within 45 km" is a much harsher test than map distance.
* **Sites, not stations**: stations within **1 km** are merged (single-linkage) into one site, because
  co-located sensors are not independent (this is the P3D `SITE_CLUSTERS` rule generalised —
  PGC+B009+B010+B011, LRIV+B005+B006+B007, B926+CLRS, B001+OSQM all collapse automatically).
* **Alive in a year** = a 3-C epoch overlapping that calendar year by >= 90 days.
* **Cells**: Slab2 interface points on a ~5 km grid with interface depth **28-45 km**, depth from the
  same `LinearNDInterpolator(+Nearest fallback, clip 12-55 km)` recipe as
  `fault_tomography/inversion/assemble_res_catalog.py` on `data/cas_slab2_input_04-18.csv`
  (etype BA/TO dropped). The file is Slab2's scattered *control-point input*, not the gridded model,
  so map contours drawn from it are indicative only (smoothed, hull-masked, in the figure).
* **Tremor weight**: PNSN `catalogs/pnsn_tremor_cascadia_full.csv` epicentres within 20 km of a cell.
  `tremor_active` = >= 100 epicentres (Pass B counts 2019+ only; the catalog ends 2026-05-27).
* **Azimuthal gap** = largest gap between the azimuths of the sites within 60 km (360 deg if < 2 sites).

## Pass B -- the full 2019+ inventory, box 46.5-49.5 N / -125.5..-121.5 E

IRIS `network=*`, `level=channel`, 2019-01-01 -> 2026-09-19, any network. NCEDC was queried too and
adds nothing in this box (it is a California archive). **All 237 detection-capable stations are
`open`** (no restricted metadata).

| | count |
|---|---|
| stations with any 3-C set >= 20 Hz | **546** |
| ... of which real 3-C **seismometers** (non-citizen) = *detection-capable* | **237** |
| ... surface broadband / borehole / surface short-period | 203 / 19 / 15 |
| accelerometer-only 3-C (strong-motion; NP, GS, urban UW) | 303 (flagged, excluded from the statistics) |
| Z-only stations | 9 (`z_only_2019plus.csv`) |
| in the 183-station project roster | 80 of 237 |
| detection-capable **sites** after the 1 km merge | **207** (74 alive 2019, 95 in 2021, 136 in 2023, **158 in 2025**) |

Networks contributing 3-C seismometers: UW 132, CC 34, PB 19, CN 18, GS 9, YR 8, X4 6, Z2 5, **PQ 4**,
C0 1, US 1. **PQ = "Public Safety Geoscience Program Canadian Research Network" (NRCan)** — four open
100 Hz broadbands on southern Vancouver Island (ALBHB, DAOB, JRDRB, LHLYB, 2021-2024) that are in no
project list and sit directly over the Bostock family field. X4/YR/Z2 are temporary experiments
(flagged `temporary_net`).

### (1) The 130 Bostock families (`family_geometry_by_year.csv`)

Sites within 45 km hypocentral, over the 130 families:

| network | year | median | mean | p90 | median <=55 km | nearest (km) | gap (deg) | >=2 sites | >=3 sites |
|---|---|---|---|---|---|---|---|---|---|
| S15 | any | 1 | 0.82 | 2 | 2 | 44.5 | 182-192 | 20% | 2% |
| S15+SYMB+BUCKS | 2019/2021 | 1 | 1.32 | 2 | 3 | 40.2 | 146-157 | 37% | 8% |
| S15+SYMB+BUCKS | 2023/2025 | 1 | 1.52 | 3 | 3 | 38.8 | 143 | 45% | 19% |
| **all open 3-C** | 2019 | 1 | 1.42 | 3 | 4 | 39.4 | 125 | 42% | 14% |
| **all open 3-C** | 2021 | 2 | 1.85 | 3 | 5 | 38.9 | 107 | 63% | 28% |
| **all open 3-C** | 2023 | 3 | 2.46 | 4 | 6 | 37.1 | 86 | 72% | 52% |
| **all open 3-C** | **2025** | **4** | **3.31** | 6 | 8 | 36.9 | 85 | **75%** | **73%** |

**The headline: the network that would support northern LFE discovery already exists and is open —
we are simply not using it.** Adding only the four PQ stations to S15+SYMB+BUCKS does nearly all of it:

| network | median sites <=45 km | >=2 | >=3 |
|---|---|---|---|
| S15 | 1 | 20% | 2% |
| S15+SYMB+BUCKS | 1 | 45% | 19% |
| **S15+SYMB+BUCKS + 4 PQ** | **3** | **75%** | **72%** |
| all 207 sites (any epoch) | 4 | 75% | 73% |

16 families still have **no** site within 45 km even with everything (bkeys 16, 17, 83, 88, 110, 112,
113, 114, 116, 118, 119, 123, 124, 127, 128, 129) — they are the 40-47 km deep families under the
Strait of Georgia, where the depth alone eats the budget; 33 families have < 2.

### (1b) Against Bostock's own 2003-2005 template array

The 21 stations Bostock used to build the family templates (15 POLARIS/C8 temporaries that closed
2004-2007 + 6 permanent CN: PGC, VGZ, SNB, LZB, PFB, YOUB), taken from
`lfe_catalog_repair/results/external/other/seisbench_lfestacks_cascadia/metadata.csv`, merged into
sites the same way. Per Bostock family, sites within 45 km hypocentral
(`bostock_array_comparison.csv`, `bostock_array_comparison_summary.csv`):

| network | sites | median <=45 km | mean | p90 | >=1 | >=2 | >=3 | median <=55 km | median nearest (km) |
|---|---|---|---|---|---|---|---|---|---|
| **Bostock 2003-2005 array** | 21 | **5** | 5.78 | 11 | 84% | 78% | **77%** | 13 | 37.5 |
| S15 (2025-alive) | 8 | 1 | 0.82 | 2 | 60% | 20% | 2% | 2 | 44.5 |
| S15+SYMB+BUCKS | 10 | 1 | 1.52 | 3 | 87% | 45% | 19% | 3 | 38.8 |
| **S15+SYMB+BUCKS+PQ4** | 14 | **3** | 3.12 | 6 | 87% | 75% | **72%** | 7 | 37.1 |

Bostock could see these families because he had, typically, **five** sites inside 45 km. Our working
set has one. Adding the four open NRCan PQ broadbands restores the fraction of families with >= 3
close sites to essentially the 2003-2005 level (72% vs 77%) with 14 sites instead of 21 — that is the
cheapest available fix for northern discovery, and it needs no new hardware, only downloads.

### (2) Interface cells (`cell_geometry_by_year.csv`), 1453 cells, 987 tremor-active

Tremor-active cells, sites <= 45 km:

| network | year | median | mean | p90 | median <=55 | nearest | gap | >=3 sites |
|---|---|---|---|---|---|---|---|---|
| S15(+SYMB+BUCKS) | any | 0 | 0.3-0.5 | 1-2 | 0 | 61 | 360 | 0-4% |
| all open 3-C | 2019 | 1 | 1.67 | 3 | 3 | 41.3 | 153 | 17% |
| all open 3-C | 2021 | 2 | 2.20 | 5 | 4 | 40.4 | 136 | 31% |
| all open 3-C | 2023 | 2 | 3.15 | 7 | 7 | 39.2 | 108 | 46% |
| all open 3-C | **2025** | **3** | 3.43 | 7.4 | 8 | 39.1 | 100 | **50%** |

By 0.5 deg latitude band (all open 3-C, tremor-active cells, median sites <= 45 km, 2019 -> 2025):

| band | 2019 | 2025 | gap 2025 | >=3 sites 2025 |
|---|---|---|---|---|
| 46.5-47.0 (Chehalis / south Olympia) | 3 | **5** | 94 | 75% |
| 47.0-47.5 (Shelton - Tacoma - Hood Canal) | 1 | **4** | **68** | 68% |
| 47.5-48.0 (Olympics interior / Seattle) | 1 | 2 | 105 | 45% |
| 48.0-48.5 (north Olympics, the family field) | 1 | **4** | 72 | 76% |
| 48.5-49.0 (Salish Sea / Victoria) | 1 | 2 | 137 | 24% |
| 49.0-49.5 (south BC) | 1 | 1 | 186 | 6% |

### (3) Ranked stations outside S15+SYMB/BUCKS (`station_ranking.csv`, 148 candidate sites, 2025-alive)

Greedy = the site that lifts the most units to **>= 3 sites within 45 km**, added one at a time on top
of S15+SYMB+BUCKS (north: families; Puget: tremor-weighted cells 46.5-48.0 N).

**North (Bostock families)** — the entire useful set is four NRCan PQ broadbands plus PTRF:

| # | site | net | lat, lon | families <=45 km | <=55 km | greedy gain | 3-C since |
|---|---|---|---|---|---|---|---|
| 1 | **JRDRB** | PQ | 48.433, -124.054 | 53 | 75 | +20 families to >=3 | 2021-04 |
| 2 | **LHLYB** | PQ | 48.536, -123.730 | 66 | 104 | +30 | 2023-11 |
| 3 | **DAOB** | PQ | 48.520, -123.419 | 36 | 90 | +12 | 2024-02 |
| 4 | **ALBHB** | PQ | 48.390, -123.488 | 53 | 92 | +6 | 2022-01 |
| 5 | PTRF | CN | 48.544, -124.413 | 13 | 38 | +2 | 2019 (roster) |
| - | ODUC / RNWY | UW | 47.97 / 48.61 | 7 / 2 | 21 / 33 | small | 2024 / 2021 |

**Olympics / Puget (tremor-weighted cells)** — a deep UW bench, all built 2021-2025:

| # | site | net | lat, lon | cells <=45 km | tremor-weighted | 3-C since |
|---|---|---|---|---|---|---|
| 1 | ODUC | UW | 47.970, -123.872 | 81 | 3.5 | 2024-05 |
| 2 | OHOH | UW | 47.860, -123.932 | 73 | 3.9 | 2022-04 |
| 3 | ODEER | UW | 47.949, -123.265 | 44 | **14.3** | 2025-06 |
| 4 | B013+B943 | PB | 47.813, -122.911 | 21 | 10.1 | (borehole, long record) |
| 5 | EQUIL | UW | 47.79, -122.86 | 22 | 8.7 | 2022 |
| 6 | DOSE | UW | 47.717, -122.972 | 18 | 9.5 | 2019 |
| 7 | SKOKO | UW | 47.409, -123.367 | 67 | 12.0 | 2022-07 |
| 8 | EAGLE | UW | 47.258, -122.999 | 38 | 13.8 | 2021-11 |
| 9 | PAN4H | UW | 47.176, -123.270 | 65 | 10.3 | 2021-05 |
| 10 | WINDI | UW | 47.314, -123.332 | 66 | 10.7 | 2021-11 |
| also | B018 (PB borehole, 46.98/-123.02, 60 cells, tw 13.8, full record), WYNO, DOORS, GNW, HDW, TKEY | | | | | |

## Pass A -- project roster only, Puget vs Oregon vs North (original brief)

Files `cells.csv`, `north_families.csv`, `summary.csv`, `sites_all.csv`, `sites_project.csv`,
`stations_epochs.csv`, figure `figures/t_puget_geometry.png`. Station pool = `f1_stations.csv`
(a superset of the 183 roster and the 277 post-2020 candidates) in 42.8-49.8 N; two pools are
carried: `project` (roster+candidates, 154 sites) and `all` (+ FDSN extras, 351 sites).

Median sites <= 45 km over **tremor-active** cells, 2024:

| region | pool=project | pool=all | p90 (all) | gap (all) |
|---|---|---|---|---|
| North 48-49 N | 2 (mean 2.0) | 2 (mean 2.7) | 6 | 113 |
| Puget 47-48 N | 2 (mean 2.5) | **3 (mean 4.0)** | 9 | 92 |
| Oregon 43.5-46 N | 2 (mean 2.0) | 3 (mean 2.9) | 6 | 113 |

**Plain answer to the original question: yes, but in one place — the Puget band 47.0-47.5 N.**
There the 2024 all-network median is **5 sites within 45 km** (mean 5.1, p90 10, azimuthal gap 68 deg),
against the northern discovery network's 1-2. Oregon's only comparable patch is **45.5-46.0 N**
(median 4, gap 68 deg, i.e. the Portland - Mt Hood corridor); south of 44.5 N Oregon is *worse* than
the north (median 1-2, gap 155-180 deg with the project roster). Puget's density is new: the band's
median went 1 (2020) -> 3 (2023-2026), built by ~10 UW 3-C stations installed 2020-2022 (WINDI,
SKOKO, WYNO, ANGUS, PAN4H, HERD, OHOH, MONTE, DOORS, EAGLE); Oregon's rise (mean 2.0 -> 2.9) is the
UO buildout (DEAN, TDWTR, VITI, LOGZ, OLBLU, GROND).

## Caveats

* Geometry only. Nothing here says a station is *quiet* enough to detect LFEs — hypocentral distance
  and azimuthal gap are necessary, not sufficient. Site noise is the untested half.
* Pass A's North family reference is epoch-blind for the `P3Dfrozen15` variant (all 15 given a
  2003-2026 interval); its `S15` variant uses the real epochs, which is why S15 collapses to the
  frozen-15 numbers from 2021 (SOKB/LZB/PFB/YOUB died 2014-2018, SNB 2021-03-28).
* `summary.csv` year 0 = "ignore epochs" (the discovery-era geometry).
* Strong-motion (accelerometer-only) stations are inventoried and counted per cell
  (`n45_strongmotion`) but never treated as detection-capable.
* Tremor weighting uses epicentres, not depths, so it smears the down-dip edge by a few km.

## Reproduce

```
python results/p3d_puget_geometry/t_fetch_epochs.py   # Pass A station epochs (IRIS, level=channel)
python results/p3d_puget_geometry/t_geometry.py       # Pass A cells/summary
python results/p3d_puget_geometry/t_figure.py         # figures/t_puget_geometry.png
python results/p3d_puget_geometry/t_inventory.py      # Pass B: all networks 2019+ in the box
python results/p3d_puget_geometry/t_geometry2.py      # Pass B: families, cells, ranking
python results/p3d_puget_geometry/t_figure2.py        # figures/t_stations_2019plus_map.png
```

`figures/t_stations_2019plus_map.png` panel (b) also carries the Bostock 2003-2005 array (filled
orange = POLARIS temporaries, open orange = the 6 permanent CN), labelled, so the old and new
geometries can be read against each other directly.
