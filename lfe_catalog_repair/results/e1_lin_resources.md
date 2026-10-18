# E1 — Resources for replicating the Lin et al. (2024) deep-learning LFE catalog network-wide

Compiled 2026-09-10. All downloads under `lfe_catalog_repair/results/external/`.
Sizes, SHA256 and licences are recorded per item in §7.

---

## 0. Headline answers

| Question | Answer |
|---|---|
| Are Lin's **trained weights** public? | **No.** The Zenodo/GitHub release is scripts only — no `.tf`/`.h5`/checkpoint files, no weight download link in the paper's availability statement. |
| Is Lin's **code** public? | **Yes** (Zenodo 10.5281/zenodo.10076720 + GitHub `jiunting/Cas_LFE`, CC-BY-4.0 on Zenodo, no LICENSE in the repo). Downloaded. |
| Is Lin's **catalog** public? | **Yes** — Zenodo 10.5281/zenodo.10016020, one 101 MB CSV, 1,058,114 located LFEs, CC-BY-4.0. Downloaded (no 504/403 this time; `curl -L` to the Zenodo API `/content` endpoint works). |
| Does the catalog include station picks? | **No.** Only located events (window start, origin time, lon/lat/depth, residual, N stations). The per-station pick files (`Detections_S_*/cut_daily_*.csv`) were never released. |
| Is there a **usable pretrained LFE picker** for Cascadia? | **Yes, from the other paper**: SeisBench `LFEDetect` (Münchmeyer et al. 2024) ships pretrained weights including a `cascadia` model trained on the same Bostock (2015) catalog, plus joint multi-region models. Downloaded (8 weight sets, ~1.1 MB each). |
| Is the **training data** obtainable? | **Yes, fully.** (a) Bostock (2015) origin times + family hypocentres + magnitudes from the Slow Earthquake Database — **downloaded complete, 269,423 events / 130 families, an exact match to the published count**; (b) per-family, per-station **P and S picks** from the SeisBench `LFEStacksCascadiaBostock2015` dataset (CC-BY-4.0) — the public equivalent of Lin's private `sav_family_phases.npy`; the (lat,lon,depth) join between (a) and (b) is verified for 128/130 families. Both downloaded. |
| Realistic completeness outside the VI/Puget network? | **Poor without retraining/transfer learning.** Lin's own §4.3 test on PB stations from central to southern Cascadia: "the overall results are noisy and may include many false detections". Münchmeyer: template matching "allows … to identify LFEs with fewer stations than deep learning", and the DL Cascadia catalog was "substantially smaller" than the reference. See §6. |
| Recommendation | **Do not try to reproduce Lin's weights.** Run SeisBench `LFEDetect` (cascadia + joint weights) as the picker, PyOcto for association, NonLinLoc (already installed) for location; fine-tune on Cascadia-wide data using the Bostock picks + our own hard negatives. See §8. |

---

## 1. Lin et al. (2024) — paper, code, model design

**Citation.** Lin, J.-T., Thomas, A. M., Bachelot, L., Toomey, D. R., Searcy, J., & Melgar, D. (2024).
*Detection of Hidden Low-Frequency Earthquakes in Southern Vancouver Island with Deep Learning.*
**Seismica 2(4)**. https://doi.org/10.26443/seismica.v2i4.1134
Corresponding author Jiun-Ting Lin (LLNL, lin51@llnl.gov). LLNL-JRNL-855845.
Funding NSF 1848302, USGS NEHRP G22AP00299.

**Code.** Lin, J.-T. (2023b). *jiunting/Cas_LFE: First release of Cascadia LFE detection scripts.*
Zenodo https://doi.org/10.5281/zenodo.10076720 (v1.0.1, 218 KB zip, CC-BY-4.0).
Live repo: https://github.com/jiunting/Cas_LFE (default branch `main`, HEAD `ea03213`, no LICENSE file,
1 star; last push 2024-03-28). Both downloaded; the GitHub HEAD adds one file over the Zenodo zip
(`daily2input_S_20240326.py`, the PB-network inference driver).
Availability statement in the paper: *"The codes for LFE detection were taken from Thomas et al. (2021)
and are available at https://zenodo.org/records/10076720"*. **No model weights are named anywhere.**

**Architecture** (`unet_tools.py::make_large_unet` / `make_large_unet_drop`; SI Figure S1):
- 1-D **U-Net**, TensorFlow/Keras, taken verbatim from Thomas et al. (2021).
- Encoder Conv1D(32f,k=21) → MaxPool → Conv1D(64f,k=15) → MaxPool → Conv1D(128f,k=11) → MaxPool →
  Flatten → Dense(188) → reshape; decoder mirrors with UpSampling1D + skip concatenations;
  output Conv1D(1, k=21, **sigmoid**) → Flatten. Loss binary cross-entropy, Adam lr=1e-4.
- Filter count scaled by `fac` ∈ {0.5, 1, 2} ("size 0.5/1/2"). **Chosen: size 2** (SI Fig S1). Optional
  dropout variant. Batch size 32, **30 epochs**, ~700,000 waveforms per epoch.
- **Two separate models**, one for P and one for S (`-pors` flag). Only the **S** model is used for
  location; the P model's picks were discarded (σ = 4.2 s, frequently latching onto S).

**Input** (`daily2input_*.py::ZEN2inp`, paper §2.2):
- 3-component, **100 Hz**, **15 s windows** (1500 samples) → input tensor (1500, **6**).
- Pre-processing: `detrend('linear')` → `taper(0.02)` → **`highpass 1.0 Hz`** → interpolate to 100 Hz.
  **No instrument-response removal** — but amplitudes are normalised per window by the max over all
  three components before the log transform, so gain sensitivity is largely removed.
- Feature map per component: `x → [ln(|x|+ε), sgn(x)]` with ε = 1e-6 (Thomas et al. 2021 "log features").
  6 channels = 3 components × (log-amplitude, sign).
- Inference: non-overlapping 15 s windows, decision **threshold y ≥ 0.1**.

**Training labels** (§2.1):
- Source catalog = **Bostock et al. (2015)**, "downloaded from the slow earthquake database (Kano et al. 2018)".
- **Per-arrival picks**, not origin times: for each catalogued *arrival-time pick* they cut a 30 s window
  centred on the pick; the target is a **Gaussian of σ = 0.4 s** centred on the P or S arrival.
- Noise class: an equal number of windows taken before the P arrival, at random times ≥ 180 s from any
  known LFE; target identically zero.
- ~**500,000 waveforms** each for P and for S, plus a similar amount of noise.
  **269,422 events** total; **25 % of events held out** for testing, split **by event ID** (not by trace)
  to prevent leakage.
- Data generator augmentation: random 15 s sub-window from the 30 s trace, so the pick lands anywhere.

**Stations.** 14 training/testing stations, networks **PO** (POLARIS, all HH) and **CN** (mixed BH/EH/HH):
`LZB, YOUB, KLNB, TWKB, SNB, SILB, TSJB, TWGB, PFB, GOWB, NLLB, PGC, SSIB, VGZ`
(explicit list in `cat2nlloc.py:44`; channel map for CN in `daily2input_S_20220919.py:150`).
10 **unseen** test stations = the **C8** network: `BPCB, GLBC, JRBC, LCBC, MGCB, PHYB, SHDB, SHVB, SOKB, TWBB`
(`daily2input_S_20220919.py:456`). Wider-region test (SI Fig S9) used **PB** boreholes B004→B036.

**Association** (§2.3): purely temporal — a detection is any 15 s window with y ≥ 0.1; an event requires
**≥ 3 stations detecting in the same 15 s window**. (`find_template.py` additionally offers a 30 km
centroid-distance filter and N_min = 5 for template building.) No declustering, no PyOcto/GaMMA.
`detc2PhaseLink.py` / `PhaseLink_classification.py` exist but are exploratory; `cat2nlloc.py` /
`plot_NLLoc_result.py` are a NonLinLoc side-branch — **the published catalog used neither**.

**Location** (§2.3): direct **grid search** over M = 1,024,800 nodes — 1 km spacing, 120 (lon) × 140 (lat)
× 0–60 km depth, centred on (−123.75, 48.7). Travel times by the **Toomey et al. (1994)** ray/finite-
difference scheme via the **Stingray** MATLAB package (`Locating/stingray_utils.py`,
`Locating/make_input_files.py` writes `srControl/srGeometry/srStation/srModel/srElevation` .mat files).
Origin time eliminated by demeaning the residual vector (their Eq. 2), so only S arrivals are needed;
minimum-misfit node wins (Eq. 3).

**Velocity model**: **Savard, Bostock & Christensen (2018)**, G-Cubed 19(6),
https://doi.org/10.1029/2017GC007417 — read from a local file `Savard_VpVs.txt`
(`longitude latitude depth Vp Vs DWS`), interpolated to the 1 km Stingray grid, with topography from
`svi_topo.nc` to convert the sea-level datum to elevation datum.
**Not distributed** — not in the repo, not in the paper's data statement, no supplement file at the
publisher (Wiley returns 403 to automated fetches). Would have to be requested from the authors, or
extracted from Savard's UBC thesis (https://open.library.ubc.ca/cIRcle/collections/ubctheses/24/items/1.0371609).
See §5 for substitutes.

**Reported performance.** At threshold 0.1: accuracy 92 % (S) / 90 % (P); precision ~95 % both;
recall ~90 % both. ROC AUC: S 0.97 (all test data), 0.98 (M > 2.2 only); P 0.92 / 0.96.
S arrival-time misfit −0.2 ± 3.6 s at threshold 0.1, −0.17 ± 2.6 s at 0.5. P misfit σ = 4.2 s (unusable).
Single-event location scatter **~22 km**; family-centroid averaging → 10 km (N = 10), 8 km (N = 30);
**depth essentially unresolved** (S-only).

---

## 2. The Lin LFE catalog (downloaded)

`results/external/lin_catalog/EQloc_001_0.1_3_S.csv` — 106,324,644 B, md5 matches Zenodo, CC-BY-4.0.
Citation: Lin, J.-T. (2023a). *A Deep Learning-Based Low-Frequency Earthquake Catalog in Southern
Vancouver Island.* Zenodo. https://doi.org/10.5281/zenodo.10016020 (published 2023-10-18).
Filename decodes as `EQloc_<thresh 0.01? see note>_<0.1 = y-threshold>_<3 = N_min>_S`.

Columns and content (1,058,114 rows, verified):

| column | meaning |
|---|---|
| `starttime` | start of the 15 s association window (UTC, always :00/:15/:30/:45) |
| `OT` | grid-searched origin time |
| `lon`,`lat` | best-fit grid node (−124.575…−122.939 E, 48.068…49.320 N) |
| `depth` | **negative-down**, 0 to −60 km, 1 km grid |
| `residual` | mean demeaned travel-time misfit, s |
| `dt` | `OT − starttime` (−33.2 … +11.2 s) |
| `N` | number of stations in the association |

Time span **2005-01-01 to 2017-02-21** (the paper says the model was run on 2005–2018; detections
collapse after early 2014 because only NLLB, PFB, PGC remained). Yearly counts:
2005 196,398 · 2006 66,953 · 2007 86,737 · 2008 131,521 · 2009 138,434 · 2010 138,959 · 2011 95,421 ·
2012 94,533 · 2013 78,089 · 2014 19,394 · 2015 8,191 · 2016 3,390 · 2017 94.

**Quality caveats measured here (important before using it as ground truth):**
- N = 3 for **46 %** of events; N ≥ 4 for 54 %, N ≥ 5 for 30.5 %.
- `residual < 0.5 s` (the paper's "high-quality" cut, used for its Fig. 10) keeps only **41.3 %** (436,687);
  `< 1.0 s` keeps 53.4 %.
- **Depths pile up on the grid boundaries**: 34 % of events fall in the deepest bin (50–60 km) and 15 %
  in 0–10 km, i.e. a third of the catalog is pinned against the 60 km bottom plane. Depth carries no
  information; treat the catalog as epicentres only, and preferably as family/cluster centroids.

---

## 3. Bostock LFE catalog and per-station arrival information

### 3a. Slow Earthquake Database (event times, family hypocentres, magnitudes) — downloaded
The old URL `www-solid.eps.s.u-tokyo.ac.jp/~sloweq/` is **retired** ("no longer maintained"); the live
database is **https://eqdb.sakura.ne.jp/**. Cascadia holdings: `Bostock2015-LFE`, `Ducellier2022a-LFE`,
`Ducellier2022b-LFE`, `Shelly2025-LFE`, `PNSN-Tremor`.

There is no documented API, but the site's own AJAX endpoint works from `curl`. **Recipe** (recorded so
it can be repeated):

```bash
FN=sloweq.$RANDOM.csv
curl -s -X POST https://eqdb.sakura.ne.jp/php/sqlIntrTrem.php \
  -H 'Referer: https://eqdb.sakura.ne.jp/' \
  --data-urlencode 'sdate=2008-01-01' --data-urlencode 'span=0366' \
  --data-urlencode 'tbldb=data' \
  --data-urlencode 'catalogs[]=Bostock2015-LFE' \
  --data-urlencode "output=/home/eqdb/www/download/$FN" \
  --data-urlencode 'format=year,month,day,hour,min,sec,lat,lon,dep,mag,catalog,ref,udate' \
  --data-urlencode 'policies[]=<the policy text>'
sleep 10 && curl -L "https://eqdb.sakura.ne.jp/download/$FN" -o out.csv
```
Traps: the `catalogs` value must be the **button id** (`Bostock2015-LFE`), not the `value` attribute
(`Bostock2015`) — with the latter the query silently returns 0 rows. Omit `output` to get JSON directly.
Requests whose result set is large enough to time out return **HTTP 500 / a 404 on the download file**;
break the request into shorter `span`s (5–10 days works everywhere).

**Downloaded and COMPLETE:** `results/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv` —
**269,423 rows, an exact match to the published event count** of Bostock et al. (2015), and **130 unique
`(lat, lon, dep)` triples = the 130 Bostock families** (depths 28.79–46.6 km, lat 48.12–48.93,
lon −124.34…−123.10; magnitudes 0.62–2.58; 343–4,521 events per family, median 1,943).
Built from the yearly CSV chunks (`bostock/chunks/`, 152,260 events) plus the July-2004 and Sept-2005
ETS episodes recovered in 5-day JSON windows (`bostock2015_sloweq_2004-2005_json_partial.csv`,
117,163 events) — those two windows are the ones whose result sets are large enough to 500/404 the CSV
route. Per-year: 2003 32,978 · 2004 59,659 · 2005 57,504 · 2007 3,147 · 2008 14,224 · 2009 22,112 ·
2010 18,114 · 2011 18,663 · 2012 28,149 · 2013 14,873. 2006 and 2014 are genuinely empty (no ETS
episode in the catalog).

**The DB record carries no family label, no station picks and no arrival times** — only the origin time,
the repeated family hypocentre, and the magnitude. Family identity is recovered by grouping on the
`(lat, lon, dep)` triple, and **the join to the SeisBench pick metadata (§3b) is verified**: at 2-decimal
rounding, 128 of the 130 sloweq families match 128 of the 129 SeisBench families on the
`(lat, lon, depth)` key; the remaining two need a nearest-neighbour match. This closes the loop —
origin times (here) + per-station P/S picks (§3b) = Lin's `total_mag_detect_0000_cull_NEW.txt` +
`sav_family_phases.npy`, fully reconstructed from public sources.

**Licence:** the Slow Earthquake Database user policy is restrictive — *"Refer designated references …
Acknowledge our project. Not for commercial use. **No redistribution.**"* Fine for internal use and for
citing; **do not republish these files.**

### 3b. The public equivalent of Lin's `sav_family_phases.npy` — downloaded, CC-BY-4.0
Lin's training pipeline needs, per family, the **station-by-station P and S arrival times** plus the
family hypocentre and a `catShift`. Those live in two files that were never released
(`total_mag_detect_0000_cull_NEW.txt`, `sav_family_phases.npy`).

**They are recoverable from SeisBench**: `seisbench.data.LFEStacksCascadiaBostock2015` —
*"Low-frequency earthquake stacks underneath Vancouver Island … based on the catalog by Bostock et al
(2015). Compiled to SeisBench format by Münchmeyer et al (2024)", licence **CC BY 4.0***.
Downloaded directly from the SeisBench mirror (no library needed) to
`results/external/other/seisbench_lfestacks_cascadia/{metadata.csv, waveforms.hdf5}`.

Contents (verified): **1,817 three-component stacks, 20 Hz**, covering **129 of the 130 families** at
**21 stations** (`GLBC GOWB JRBC KELB KHVB LCBC LZB MGCB PFB PGC SHVB SILB SNB SOKB SSIB TSJB TWBB TWGB
TWKB VGZ YOUB`, network PO; HHZ 1230 / BHZ 335 as the naming channel), median 14 stations per family
(6–21), median 874 events per stack (5–3693), median S SNR 5.45 dB.
Per-trace fields include `source_family_id`, `source_latitude/longitude/depth_km`,
`trace_S_arrival_sample` (1453 traces), `trace_P_arrival_sample` (1067), `trace_P_predicted_arrival_sample`
(750), `path_P/S_residual_s`, `source_horizontal/vertical_uncertainty_km`, `trace_P/S_snr_db`, and a
train/dev/test `split`.

⇒ Join the family hypocentres in §3a with `source_*` here to attach the family ID, and the training
windows Lin used are then reproducible: `arrival = origin_time + (stack pick − stack origin)`, cut from
FDSN data.

### 3c. Bostock's own distributions
Bostock et al. (2012, G-Cubed, 10.1029/2012GC004391) and Bostock et al. (2015, JGR,
10.1002/2015JB012195) do not release machine-readable family/arrival files; the Savard & Bostock (2015)
BSSA electronic supplement (bssa_105-4/2014301-esupp) contains **figures only** — checked, no catalog,
no templates, no picks. Wiley blocks automated fetching (403), so the JGR/G-Cubed supporting-information
items were not enumerated; if per-station template arrival tables are wanted beyond §3b, ask
M. Bostock (UBC EOAS) directly.

---

## 4. Other public Cascadia LFE catalogs, pickers and tremor tools

| Resource | What it is | Access / status |
|---|---|---|
| **Münchmeyer et al. (2024)** — *Deep learning detects uncataloged low-frequency earthquakes across regions*, Seismica 3(1), https://doi.org/10.26443/seismica.v3i1.1185 (preprint arXiv:2311.13971) | PhaseNet-based LFE **phase picker** trained on LFE *stacks* + INSTANCE noise, for Cascadia (Bostock 2015), San Andreas (Shelly 2017), Guerrero (Frank 2014), Nankai (JMA); workflow = SeisBench picker → **PyOcto** associator → **NonLinLoc** with a 1-D model. | **Weights public via SeisBench** (`seisbench.models.LFEDetect`). All 8 weight sets downloaded: `cascadia, mexico, san_andreas, nankai, cms, jcms, jcs, jms` (~1.1 MB each). Datasets CC-BY-4.0. |
| **Shelly et al. (2025)** — *A northeast-dipping zone of low-frequency earthquakes at the southern edge of Cascadia subduction*, GRL 52, e2025GL116116 | **61,441 detections, 27 LFE families**, 2018-01-01→2024-11-30, Mendocino/Gorda edge. Threshold 9 × daily MAD of the channel-summed CC, mean CC ≥ 0.12. | USGS data release **10.5066/P1TCKK7G**, **CC0**. Downloaded (`other/LFEcat_MTJ_2018-2024_shelly2025.csv`, 4.31 MB, 61,440 rows: time, s_of_day, ccsum, thresh, cc_div_thr, meancc, nchan, family ID, idnum). **⚠ The local copy at `catalogs/mtj_ducellier/LFEcat_MTJ_2018-2024.csv` (1,205 B) is a saved 404 HTML page, not data.** |
| **Ducellier & Creager (2022)** — *An eight-year-long LFE catalog for Southern Cascadia*, JGR 10.1029/2021JB022986 | ~478 k LFEs, southern Cascadia (FAME + permanent networks), 66 families; per-family pandas frames (`year…second, cc, nchannel`) + **templates** + relative travel times. | Zenodo 10.5281/zenodo.6374239 (`catalog-v3.zip`, 587 MB, "other-open"). **Already local** and complete: `catalogs/ducellier_catalog/` (2.5 GB clone) — `data/Ducellier/catalogs/<family>/catalog_2007_2009.pkl` + `catalog_2004_2011.pkl` for 66 families, `src/templates*/`, `src/timearrival/*.pkl`, `src/traveltimes/`, plus MATLAB `locate_family.m`. Also in the Slow Earthquake DB as `Ducellier2022a-LFE` / `Ducellier2022b-LFE`. Not re-downloaded. |
| Local `catalogs/ducellier_b935/` | 34 per-family CSVs, single `time` column (detection times) — a derived extract, no picks. | already local |
| **Ducellier & Creager (2022)** — *Depth and thickness of tectonic tremor in the NE Olympic Peninsula*, JGR 10.1029/2021JB022708 | **Tremor** depths, not an LFE catalog. Code: Zenodo 10.5281/zenodo.5725990 (`timelags-v2.0.zip`, 23.7 MB). | public, not downloaded (tremor-only) |
| **Thomas et al. (2021)** — *Identification of LFEs on the San Andreas fault with deep learning*, GRL 10.1029/2021GL093157 | **The architecture Lin reused**, with **trained Parkfield weights**. `amtseismo/Parkfield-LFE-CNN` (Zenodo 10.5281/zenodo.4777132): TF checkpoints for size 0.5/1/2 × σ = 0.05/0.1/0.2/0.4 × drop/no-drop, plus `apply_cnn.py`, `lfe_unet_tools.py`, and 500-sample LFE/noise training h5s. | Downloaded the two configs matching Lin's choice (`large_2.0_unet_lfe_std_0.4.tf`, `drop_large_2.0_…`) + scripts + both h5. **No LICENSE file in the repo** — cite Thomas et al. (2021); check with the author before redistributing. |
| **Wech tremor detector (WECC)** | The PNSN tremor catalog behind `catalogs/pnsn_tremor_*.csv`. | **Detection/location code is NOT public.** PNSN's GitHub org publishes only `tremor_api` and `tremor-client` (the web API/UI). The catalog itself is served by that API (already used in this project) and is mirrored in the Slow Earthquake DB as `PNSN-Tremor`. |
| **SeisBench pretrained models with LFE classes** | Only `LFEDetect` (above). All other SeisBench models (PhaseNet, EQTransformer, GPD, EQCCT, …) are earthquake pickers with no LFE class. Münchmeyer's workflow deliberately runs an ordinary PhaseNet alongside and *removes* events detected as both. |

---

## 5. Tools and velocity models

| Tool | Status on this pod |
|---|---|
| **NonLinLoc** | **already installed**: `/home/jovyan/.local/bin/{NLLoc,Vel2Grid,Grid2Time,…}` |
| **seisbench** | installable, `pip install --user seisbench` resolves to **0.12.5**, dry-run clean (torch 2.12.0+cpu, obspy 1.5.0 already present). Not installed yet. |
| **pyocto** | installable, **0.2.0** on PyPI, dry-run clean. Source https://github.com/yetinam/pyocto (MIT). |
| **GaMMA** | **not on PyPI under a usable name** (`gamma` 0.0.1 is an unrelated squat). Install from source: `pip install git+https://github.com/AI4EPS/GaMMA` (MIT, active). PyOcto is the better fit here — it is what Münchmeyer's LFE workflow uses. |
| Python grid-search alternative | `pykonal` 0.4.1 on PyPI (Eikonal solver) if a NonLinLoc-free path is wanted; Lin's own approach (demeaned S-only travel-time misfit on a precomputed grid) is ~40 lines of numpy once a travel-time table exists. |
| torch / tensorflow | torch 2.12.0 **CPU-only** build is installed; **tensorflow is NOT installed** — required to *run* Lin's or Thomas's code as-is. A GPU is available on this pod (see PLAN.md), so a CUDA torch build + SeisBench is the lower-friction route than adding TF. |

**Velocity models**

- **Savard et al. (2018)** — what Lin actually used. Not publicly distributed (see §1). Request from authors.
- **Stephenson et al. (2017) casc1.6** — USGS OFR 2017–1152, data release **10.5066/F7NS0SWM**.
  `casc1.6_velmdl.zip`, 1.43 GB, whole margin (40.2–50 N, 121–129 W, to 59.4 km); it holds a single
  ASCII file `casc1.6_velmdl.txt` that **expands to 17.7 GB** — unzip selectively / stream it.
  **Downloaded** to `results/external/velocity/`. (The existing
  `/home/jovyan/cascadia-sse-strain/data/cascadia_casc16_profile.csv` is a 64-row 1-D extract of this.)
- **CVM v1.7 (2025)** — Stephenson/Moschetti et al., USGS OFR 2025–1045, data release **10.5066/P14HJ3IC**,
  CC0. netCDF4 in three layers: `CVM17_L01.nc` 6.18 GB (0–1.2 km, 200 m), `CVM17_L2.nc` 3.62 GB
  (1.5–9.9 km, 300 m), **`CVM17_L3.nc` 763 MB (10.8–59.4 km, 900 m — the LFE depth range)**.
  ⚠ **Not scriptable**: the large files sit behind ScienceBase's JS file manager; both
  `.../manager/download/<cuid>` and `.../catalog/file/get/...?f=__s3__` return the 4 KB SPA shell.
  Only the FGDC metadata came down (`velocity/CVM17_metadata.xml`). Fetch these three by hand from
  https://www.sciencebase.gov/catalog/item/65b40c6ad34e36a390458d76 if v1.7 is wanted; casc1.6 is
  scriptable and adequate for a 1-D or smoothed 3-D LFE grid.
- For a 1-D starting model, Münchmeyer's workflow used a plain 1-D model with NonLinLoc; the project
  already has `Crosson (1976)` in use (RESULTS.md stage-4) and Slab2 depths in
  `data/station_slab2_depth.csv`.

---

## 6. Expected completeness outside the Vancouver Island / Puget network

This is the crux, and both papers are explicit.

**Lin et al. §4.3 (verbatim):** *"We test the model performance in a broader region using the Plate
Boundary Observatory stations (PB) … spanning from central to southern Cascadia. The comparison with the
SSEs … reveals that although the peak LFE detections align with the timing of large SSE occurrences,
specifically for the northern stations close to southern Vancouver Island (i.e. B001-B013), the overall
results are noisy and may include many false detections (Supplementary Figure S9). We note that even
though the model's statistics indicate a low false detection rate, this assessment is based on testing
data recorded in southern Vancouver Island similar to the training data … noisier performance is
expected due to different noise characteristics and source-station paths … **Transfer learning … with
additional training data in this region may improve the detections**."*
And §4.1: *"It remains to be seen how far beyond the study region the CNN can reliably identify LFEs;
this likely depends on the spatial variability of the LFE source and high-frequency noise."*

**What the numbers imply for station density**
- The < 1 % false-event figure is *derived*, not measured: Lin's §3.2 argument is "assuming there is
  only one LFE in the 15 s time window and all the detections are made independently", combined with the
  95 % precision on VI test data. It is therefore conditional on **≥ 3 stations in the same 15 s window**
  *and* on VI-like noise, and it compounds: whatever the per-station false-detection fraction q is, the
  false-event rate goes as ~q^N. With q = 0.05 (VI) and N = 3 that is 1.3e-4; if q degrades to 0.2 in an
  unfamiliar noise regime it becomes 8e-3 — ~60× worse. That is exactly the "noisy … many false
  detections" PB result, and it is why the association rule, not the picker threshold, is the knob to
  tighten off-VI (require N ≥ 4 and travel-time consistency, not just coincidence).
- Lin's own catalog collapses when the network thins: detections fall off a cliff in early 2014 *"due to
  a lack of stations in our data set, with only three stations (i.e. NLLB, PFB, PGC) available"* — with
  3 stations the association is only just satisfiable and the catalog is ~1/10 of its 2009–2013 rate.
  **Three stations within ~50 km is the hard floor; ~6+ is where the catalog becomes stable.**
- Location, not detection, is the harder constraint: 22 km single-event scatter with 4–10 stations of
  VI geometry, and depth unresolvable from S alone. Münchmeyer's San Andreas case had to be **excluded
  from analysis** because a small-aperture borehole array + high pick uncertainty gave unusable
  locations — a direct warning for any Cascadia segment covered mainly by clustered PB boreholes.
- Münchmeyer's conclusion, on the sensitivity trade-off: template matching *"leads to a more sensitive
  model, as evidenced by higher event counts. Furthermore, it allows template matching to identify LFEs
  with fewer stations than deep learning"*, and their DL Cascadia catalog was *"substantially smaller"*
  than the reference catalog for the same period. **A DL catalog is more diverse (source-agnostic, finds
  new families) but not more complete per family.**
- Encouraging counter-evidence: the model generalises to unseen **stations** within the region (Lin's 10
  C8 stations, SI Fig S3, patterns consistent with Bostock), and across **regions** at the picker level
  (Münchmeyer: AUC 0.97–1.00 above 0 dB SNR, ≥ 0.88 even at −10 dB; leave-one-region-out catalogs
  correlate 0.70–0.90 with references). The joint multi-region model beat the Cascadia-only model.

**Realistic expectation for this project**
- **Southern Vancouver Island / Strait of Juan de Fuca / north Olympic (CN+PO+C8+PB, ≥ 8 stations
  within 60 km):** replication-grade. Expect ~Lin-level event counts, epicentres to ~10 km on
  cluster centroids, no depth.
- **Puget Lowland / central Washington (dense PNSN + PB):** workable, but the picker must be
  fine-tuned — noise character (urban, borehole EH vs surface BH/HH) differs from PO/CN.
- **Oregon and northern California (sparse, PB-borehole-dominated, long inter-station distances):**
  expect an **activity time series, not an event catalog**. This is precisely where Lin's PB test was
  "noisy". Ducellier's southern-Cascadia template catalog and Shelly's MTJ catalog are the better
  ground truth / product there.
- The project's own D3 measurement (RESULTS.md: median Lin LFE at PGC has amplitude s/n 0.078 of the
  1–8 Hz ambient; 1.8 MAD per 3-C station under template matching) is **consistent** with a DL picker
  succeeding where a matched filter fails — the picker exploits the multi-second 3-C waveform *shape*,
  not a single correlation scalar — but it does not exempt us from the ≥ 3-station association rule.
  Association is what buys the < 1 % false-event rate, and association needs stations.

---

## 7. Download manifest

All paths relative to `lfe_catalog_repair/results/external/`.

| path | bytes | SHA256 (or md5 where the source publishes md5) | source URL | licence |
|---|---|---|---|---|
| `lin_catalog/EQloc_001_0.1_3_S.csv` | 106,324,644 | `702dec6c700713cf7153b59b3c259d2ab501fcb656603cc7dc2fcf77e71db8ab` (md5 `5e4bb932f36736b2c1fe67b7d2a07e28` = Zenodo) | https://zenodo.org/api/records/10016020/files/EQloc_001_0.1_3_S.csv/content | CC-BY-4.0 |
| `lin_catalog/record_10016020.json` | 4,296 | — | https://zenodo.org/api/records/10016020 | — |
| `lin_cas_lfe_code/Cas_LFE-v1.0.1.zip` | 223,536 | `c1bf7a55974b67c49273dd20b6216a9b049520ca6ac9b9ea7c005ae57e87e027` | https://zenodo.org/api/records/10076720/files/jiunting/Cas_LFE-v1.0.1.zip/content | CC-BY-4.0 |
| `lin_cas_lfe_code/zenodo_unzip/` | — | (unpacked, 56 .py files) | — | CC-BY-4.0 |
| `lin_cas_lfe_code/github_main/` | ~4.5 MB | git HEAD `ea03213` | https://github.com/jiunting/Cas_LFE | no LICENSE file |
| `lin_cas_lfe_code/Lin2024_Seismica_1134.pdf` | 8,876,570 | `c9fb66ef2f59a2fc9cefaaf1652762ac797276640d604530a5323ee8fcd4b967` | https://seismica.library.mcgill.ca/article/download/1134/1681/11898 | CC-BY-4.0 (Seismica) |
| `lin_cas_lfe_code/Lin2024_supplementary.docx` | 8,666,999 | `45ceb40bace844d2fe40ba4e1619c02930819f91cb91cd3d1c32975dbe75d556` | https://seismica.library.mcgill.ca/article/view/1134/1674 | CC-BY-4.0 |
| `bostock/bostock2015_sloweq_FULL_2003-2013.csv` **(use this one)** | 25,315,472 | `04cb6d76703e9fdcfe5ddfd4274bea32479240e8a9f0e2a1e858ee892d7d8679` | eqdb.sakura.ne.jp POST/JSON recipe (§3a), merged | Slow Earthquake DB policy — **no redistribution** |
| `bostock/bostock2015_sloweq_2003-2013.csv` | 14,345,945 | `bc5d878c1762c3ce5ca27d04c38ab58857e183b450ccbe69fa4d3ac30239e31b` | same (CSV route, yearly) | same |
| `bostock/bostock2015_sloweq_2004-2005_json_partial.csv` | 11,364,893 | `ca80216026ab72cde397f597983476211479cac615501be267787aafa64182d8` | same (JSON route, 5-day windows) | same |
| `bostock/chunks/bostock_*.csv` | — | — | same | same |
| `other/seisbench_lfestacks_cascadia/metadata.csv` | 559,383 | `442adeaf4db40fbda993a9766aca0745cbcff034dc34c04406c2be9e0c5947d0` | https://seisbench.gfz.de/mirror/datasets/lfestackscascadiabostock2015/metadata.csv | CC-BY-4.0 |
| `other/seisbench_lfestacks_cascadia/waveforms.hdf5` | 39,708,932 | `ccfa11fb273a157c0876dad8eea7c7be93db59e5e01a52e8791fc86f6a1cbbfb` | …/waveforms.hdf5 | CC-BY-4.0 |
| `other/seisbench_lfedetect/{cascadia,cms,jcms,jcs,jms,mexico,nankai,san_andreas}.{pt,json}` | 8 × ~1.11 MB + 8 × ~0.8 KB | see `sha256sum` in dir | https://seisbench.gfz.de/mirror/models/v3/lfedetect/<name>.{pt,json}.v1 (primary root: https://hifis-storage.desy.de/Helmholtz/HelmholtzAI/SeisBench/models/v3/lfedetect/) | SeisBench (LGPL-3.0 code); cite Münchmeyer et al. 2024 |
| `other/LFEcat_MTJ_2018-2024_shelly2025.csv` | 4,310,270 | `503bc5545b19b4a16107d9dadde43816022913ed88a328234f3993923a51790e` | https://www.sciencebase.gov/catalog/item/68278bb3d4be02693eeab6f4 (doi 10.5066/P1TCKK7G) | **CC0** |
| `other/Shelly2025_metadata.xml` | 11,855 | `bd6c796b816f0a80e793945301d01826106e9e4f93d32af131b75ba293d250bc` | same | CC0 |
| `other/parkfield_lfe_cnn/models_v2/large_2.0_unet_lfe_std_0.4.tf.data-00000-of-00001` | 43,511,943 | `83d2239d1da9b1c2b88378869f2efdbed9ea1b21fc498ae897905e1128c8dfb8` | https://raw.githubusercontent.com/amtseismo/Parkfield-LFE-CNN/main/… | no LICENSE; cite Thomas et al. 2021 (Zenodo 10.5281/zenodo.4777132) |
| `other/parkfield_lfe_cnn/models_v2/large_2.0_unet_lfe_std_0.4.tf.index` | 2,374 | `a0c6930a8a00d2e04a51ecc55b630af600f12059e1dbe9e2a51e7d2e81e7c64d` | same | same |
| `other/parkfield_lfe_cnn/models_v2/drop_large_2.0_unet_lfe_std_0.4.tf.data-00000-of-00001` | 43,511,943 | `52ac0e9cff4b9bde75fa1580170cd74d47df081f6c48da29c051d56bd4855688` | same | same |
| `other/parkfield_lfe_cnn/models_v2/drop_large_2.0_unet_lfe_std_0.4.tf.index` | 2,374 | `1c5d6fb3aec5ce8c0189ea88270add28120f942f56c74735825fc1994e0292fe` | same | same |
| `other/parkfield_lfe_cnn/parkfield_lfe_data_500.h5` | 36,014,048 | `bf81c73db383fcb22667491deffa2127aaa544cff03b0ae06cdde83217398bb5` | same | same |
| `other/parkfield_lfe_cnn/parkfield_noise_data_500.h5` | 36,014,048 | `4a081c6bad3f86adc7c3a20b7209cd4e65f3da4f0b576fd0e214b8f8b75586e6` | same | same |
| `other/parkfield_lfe_cnn/{README.md,apply_cnn.py,figure4.py,lfe_unet_logfeatures.py,lfe_unet_tools.py}` | 2,771 / 4,892 / 4,648 / 7,415 / 19,598 | see `sha256sum` in dir | same | same |
| `velocity/casc1.6_velmdl.zip` | 1,430,939,131 | `32fdddc43fbf15ff9736789b86b2f100f1c82b8439a747b09ecbc164cdda7dec` | https://www.sciencebase.gov/catalog/file/get/59f1e68be4b0220bbd9dd4b4?f=__disk__d0/92/1e/… (doi 10.5066/F7NS0SWM) | USGS public domain |
| `velocity/casc1.6_metadata.xml` | 6,761 | `7f81d4963b8b440f416d0218be19bc5519c15a4e1aa194e9d40f0a3001997e65` | same item | USGS public domain |
| `velocity/CVM17_metadata.xml` | 16,838 | `326b03badc2c792729fff9cc55d5ff8b5a23fee0da6cf98b27289481d6f2ec8d` | https://www.sciencebase.gov/catalog/item/65b40c6ad34e36a390458d76 (doi 10.5066/P14HJ3IC) | CC0 |

**Not downloaded, deliberately:** CVM17 L01/L2/L3 netCDFs (10.6 GB, manual download only — see §5);
Ducellier `catalog-v3.zip` (587 MB — already local and complete at `catalogs/ducellier_catalog/`);
SeisBench LFE stacks for Mexico/San Andreas (only needed if training the joint model from scratch);
`timelags-v2.0.zip` (tremor, out of scope).

---

## 8. Recommendation

**Replicate the *method*, not the *weights*. Use Münchmeyer's released model as the picker and Lin's
paper as the association/location design.**

1. **Picker — use SeisBench `LFEDetect`, not a Lin reimplementation.**
   Weights are public and already on disk. `cascadia.json` gives the exact contract: 3-component **ZNE**,
   **20 Hz**, **1200 samples = 60 s**, bandpass **1–8 Hz zerophase**, `norm="peak"`, classes `PSN`,
   sigmoid output, `grouping="instrument"`, default `overlap=600`. Start with `cascadia`, and benchmark
   `jcms` (joint Nankai+Cascadia+Mexico+SanAndreas) against it — Münchmeyer found the **joint model
   outperformed the Cascadia-only one**, including for Cascadia. This costs days, not weeks: no TF
   install, no 269 k-window training-set build, a CUDA torch + `pip install --user seisbench` and the
   fleet streaming pass we already have.
   Reproducing Lin's weights instead means re-deriving `sav_family_phases.npy` (§3b makes this possible)
   *and* re-downloading ~1 M 30 s training windows *and* installing TensorFlow to run 2020-era Keras
   checkpoint code — for a model that, on its own author's assessment, is *worse* outside VI.
2. **Keep Lin's model only as the calibration target.** Its catalog (§2) is the ground truth for the
   VI segment: for 2005–2013 compare our detections window-by-window against `EQloc_001_0.1_3_S.csv`,
   restricted to `residual < 0.5 s` and `N ≥ 4`, and against the Bostock origin times (§3a) which are
   the *independent* truth. Report recall against Bostock, and agreement-with-Lin separately.
3. **Association: PyOcto, with Lin's 15 s / N ≥ 3 rule as the sanity floor.** PyOcto is what the
   released picker was designed around, is MIT, installs cleanly, and unlike Lin's pure coincidence rule
   it enforces travel-time consistency — which is exactly the discriminant that keeps false-positive
   rates from exploding where the noise character is unfamiliar. Report N-station distributions; require
   N ≥ 4 for anything called an event outside VI.
4. **Location: NonLinLoc (already installed) with casc1.6, S-weighted.** Use Lin's demeaned-residual
   trick if P picks stay unusable. Publish **cluster/family centroids, never single-event depths** — a
   third of Lin's own catalog is pinned to his grid's bottom plane (§2), and both papers say depth is
   unresolved. If per-family hypocentres are wanted at VI-like precision, the §3b stacks give
   ready-made P *and* S picks at 21 stations for 129 families.
5. **Fine-tune before trusting anything south of Puget.** The training substrate now exists and is
   permissively licensed: §3b stacks (129 Cascadia families, CC-BY-4.0) + §3a origin times + Shelly's
   CC0 MTJ families (27 families, 2018–2024 — the *only* public southern-Cascadia LFE ground truth with
   family labels) + Ducellier's 66 southern-Cascadia families already on disk + **our own hard negatives**
   (the lfe_features R9 EQ/blast set and the in-tremor ambient windows that scored 0.81, plus the
   zero-amplitude degenerate windows found in WISH/YACT). Münchmeyer's noise+stack synthesis recipe
   (combine 1–3 stacks with real noise at controlled SNR) is the cheap way to make labelled data for any
   station whose noise we can sample — it needs no new picks.
6. **Set expectations by segment now, not after the run** (§6): replication-grade only where ≥ 6
   3-component stations sit within ~60 km; three stations is the floor at which Lin's own catalog
   collapsed; Oregon/northern California should be scoped as an **activity time series** plus template
   matching against Ducellier/Shelly families, not as an event catalog with hypocentres.
7. **One correctness fix noted in passing:** `catalogs/mtj_ducellier/LFEcat_MTJ_2018-2024.csv` is a
   1,205-byte saved 404 page, not the Shelly catalog. The real file (61,440 rows) is now at
   `results/external/other/LFEcat_MTJ_2018-2024_shelly2025.csv`. Left in place — outside this task's
   write scope.
