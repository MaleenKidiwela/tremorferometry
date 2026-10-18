# Southern truth inventory — what is on disk, verified row by row

Agent T, 2026-09-23 16:10 UTC. Every number below was read from the files named beside it; nothing is quoted from a paper abstract.

## 1. Shelly et al. (2025) — Mendocino / southern edge of Cascadia

* Catalogue `lfe_catalog_repair/results/external/other/LFEcat_MTJ_2018-2024_shelly2025.csv` — **61,441 detections, 27 families**, 2018-01-01 .. 2024-11-30 (USGS data release 10.5066/P1TCKK7G, CC0).
* Locations from the release metadata `Shelly2025_metadata.xml`: latitude 40.027..40.088 N, longitude -123.713..-123.532, depth 22.2..28.6 km. The 27 families occupy a patch about 7 km N-S by 15 km E-W — one tight NE-dipping zone.
* Detections per family: min 594, median 1,558, max 7,867. Every family spans essentially the whole 2018-2024 record.
* **No templates are released** — only times and locations. Our templates must be built from his times on our own stations (stage 2).
* `catalogs/mtj_ducellier/LFEcat_MTJ_2018-2024.csv` is a saved 404 HTML page (1205 bytes), not data — confirmed again here; the working copy is the one above.

## 2. Ducellier & Creager (2022) — southern Cascadia, eight years

* `catalogs/ducellier_catalog/` (the author's own Zenodo repository, complete).
* **66 families**, one directory each under `data/Ducellier/catalogs/<family>/`, with two pandas frames (`year,month,day,hour,minute,second,cc,nchannel`):
  * `catalog_2004_2011.pkl` — the PERMANENT-network catalogue: **477,715 detections**, 2004-01-01 .. 2011-12-31.
  * `catalog_2007_2009.pkl` — the FAME temporary-array catalogue: **107,924 detections**, 2007-07 .. 2009-12.
* Locations from `data/Plourde_2015/templates_list.txt` (Plourde et al. 2015 templates, which Ducellier re-used): latitude 39.45..41.72 N, longitude -123.47..-122.43, depth 15.0..46.8 km, with published uncertainties eH 0.8-5.2 km, eZ 0.3-6.4 km. The 66 template families collapse to 34 numbered Plourde families plus MF1/MF2/BC3.
* **Templates ARE released**: `src/templates_v1/<family>/<STA>_<CHA>.pkl` (obspy Trace + orientation) — 66 families, 20 Hz, response-removed to velocity, 1.5-9 Hz, 10-45 s long (constant within a family), median 16 channels per family, 47 distinct stations. `src/templates/` is the FAME set and `src/templates_2007_2009/` the intermediate set.
* His station tables are on disk too: `data/Ducellier/stations_permanent.txt` (57 rows / 47 stations, his selected permanent bed), `stations_permanent_all.txt` (83, the pool), `stations_temporary.txt` (33, FAME), and the per-network channel/epoch tables `data/networks/{BK,NC,PB,FAME}.txt`.
* **37 of 57 rows (65 %) are SINGLE-COMPONENT VERTICALS** (NC EHZ/SHZ: GBB GHM GRO GSN GTC GVV KBN KBS KCS KFP KIP KKP KOM KPP KRK LAM LBF LBK LBP LDB LGB LGP LPG LRB LRR LSF LSH LSR LTC LVR LWH KTR). His detector correlates PER CHANNEL and averages, so a vertical-only station is a full member of his network. Any loader of ours that requires 3-C would silently run a thinned network — the southern chain here therefore accepts 1-C and 2-C station-days by construction.
* His published parameters, read from `src/find_all_LFEs_parallel.py` __main__ and `src/get_data.py`: band 1.5-9.0 Hz zero-phase, 20 Hz (dt 0.05), response removed to VEL with pre_filt (0.2, 0.5, 10, 15) and water_level 80, per-channel normalised cc averaged over channels with a complete gapless block, threshold **8 x MAD** of the block, declustering `floor(1/(dt*freq0))` = 20 samples = 1 s, 12-hour blocks.

## 3. Other southern files checked

* `catalogs/ducellier_b935/` — 34 per-family CSVs with a single `time` column (45,844 rows total): a derived extract of the same catalogue restricted to families that B935 sees. No picks, no locations. Kept as a cross-check only.
* `catalogs/mtj_ducellier/` — the 404 page described above.
* Waveform trees: no southern station-day existed anywhere on disk before this run (`lfe_picker/results/o4_deferred.md` §'Why it is deferred' re-verified).

## 4. Station geometry we derived (FDSN, not assumed)

* `results/fdsn_channels_south.csv` — 30,501 channel epochs from EarthScope + NCEDC for NC, BK, PB, UW, CC, US, UO (and the XQ/FAME temporaries, kept for reference) in the box 38.6-42.4 N, -124.8 .. -121.6, 2004-2026.
* Stations alive on 2026-06-01: 204 (195 of them 3-C). Stations alive in BOTH 2009 and 2026: 59, of which 50 are 3-C today — so unlike the north (POLARIS decommissioned), the southern permanent bed can be reproduced *and* extended on the same instruments.
* Stage-1 rosters (`results/stage1_rosters.csv`): permanent stations within 80 km hypocentral, alive on the sampled detection days AND alive today, nearest 10 kept as the download budget. Shelly: 23 distinct stations, 10 per family, median hypocentral distance 38 km. Ducellier: 35 distinct stations, 7 per family (min 3), median 63 km.
* Ducellier reproduction beds: **bed A** = his exact set (47 stations, 37 Z-only rows); **bed B** = bed A plus 25 stations within 80 km of a family that exist in the window but were not in his list. Both are run and reported separately.

## 5. What stage 1 downloads

* Existence test: **12,606 station-days** (1301 distinct days, 54 stations) — the busiest days of each family, enough for >= 150 of its own detections (median 406 per family). Estimated 100-160 GB.
* Ducellier reproduction: **2,800 station-days** over 2008-04-01 .. 2008-05-10, the richest window in his catalogue (64 of 66 families with >= 20 events in 30 days). Estimated ~70 GB.
* Nothing else is fetched until a stage needs it.

