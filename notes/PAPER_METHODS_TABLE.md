# Published LFE methods vs. what our code actually does

Built 2026-09-24 by agent Z from the SIX full-text PDFs on disk plus Ducellier & Creager's LaTeX
source. **Every parameter below is quoted verbatim from the paper.** Where a paper does not state
something, the cell says `SILENT` — it is *not* filled in from a neighbouring paper. (Filling
Washington cells from Bostock 2012 is exactly the error that produced the 6-s / 51-station
Washington run.)

Sources read:
- `notes/1-s2.0-S0012821X13004780-main.pdf` + `notes/royer2014_supp.pdf` — **Royer & Bostock 2014**, EPSL 402, 247–256
- `notes/Geochem Geophys Geosyst - 2012 - Bostock - ...pdf` — **Bostock et al. 2012**, G-cubed 13, Q11007
- `notes/JGR Solid Earth - 2017 - Chestler - ...pdf` — **Chestler & Creager 2017**, JGR 122, 3099–3114
- `notes/1-s2.0-S0040195115003807-main.pdf` (= `codex-LFE/cache/papers/thomas_bostock_2015_central_oregon.pdf`, accepted-MS twin) — **Thomas & Bostock 2015**, Tectonophysics 658, 111–116
- `notes/JGR Solid Earth - 2015 - Bostock - Magnitudes and moment-duration scaling...pdf` — **Bostock et al. 2015**, JGR 120, 6329–6350
- `notes/Geophysical Research Letters - 2015 - Plourde - ...pdf` — **Plourde et al. 2015**, GRL 42, 4849–4855
- `notes/Geophysical Research Letters - 2025 - Shelly - A Northeast-Dipping Zone of Low-Frequency Earthquakes at the Southern Edge.pdf` — **Shelly, Goldberg, Wech & Thomas 2025**, GRL 52, e2025GL116116 (added 2026-09-25 by agent AA; ⚠ its **Supporting Information S1 carries the "full methods" and is NOT on disk**, so every parameter below comes from the main text, the figure captions or the acknowledgements, and the rest is SILENT)
- `notes/Geochem Geophys Geosyst - 2014 - Sweet - A family of repeating low-frequency earthquakes at the downdip edge of tremor and.pdf` — **Sweet, Creager & Houston 2014**, G-cubed 15, 3713–3721 (added 2026-09-25 by agent AA; the fuller statement of the Array-of-Arrays front end that Chestler & Creager 2017 compresses into one paragraph)
- `catalogs/ducellier_catalog/paper/agujournaltemplate.tex`, `answer_reviewers.tex`, `supplement.pdf`, and his released code `catalogs/ducellier_catalog/src/find_all_LFEs_parallel.py` — **Ducellier & Creager 2022**

---

## 0. THE HEADLINE: who located what, and how many

| Paper | Families / templates produced | **Did they locate them?** | With what |
|---|---|---|---|
| **Bostock et al. 2012** (SVI) | **140** families ("grouped into 140 families using a hierarchical cluster analysis", §2.1 ¶6) | **YES — all 140.** "we may use P- and S-traveltime picks to generate a representative hypocentral location for each LFE family… we will consider each template as representing a distinct LFE location" (§3 ¶8) | **hypoinverse** [Klein 2002]; ~16 P-times + 22 S-times per template |
| **Royer & Bostock 2014** (NVI/NW/SW) | initial templates **224 NW / 70 SW / 54 NVI** (§3.2); **final templates 122 NW / 54 SW / 47 NVI** (§3.3) | **YES — twice.** "We perform two locations, one using a standard linearized inversion (Hyp2000, Klein, 2002) and one using the double difference location algorithm (hypoDD, Waldhauser, 2001)" (§5) | **Hyp2000 then hypoDD.** hypoDD retains "91%, 92%, 77% and 72% of the original events … for the NVI, SVI, NW and SW datasets" (§5) |
| **Chestler & Creager 2017** (Olympics) | **57 initially → 43 final** ("This reduced the number of LFE families from 57 to 43", §2.2); 32,387 unique LFEs | **YES — all 43.** "To locate LFE families we first produce P and S wave differential times… We then perform a grid search using the 3D velocity model from Preston et al. [2003] to locate each LFE family" (§2.1) | **3-D grid search (Preston et al. 2003) + double-difference** [Waldhauser & Ellsworth 2000] |
| **Thomas & Bostock 2015** (central OR) | **18** families, ">15,000 individual LFE detections" (§3) | **PARTLY — 13 of 18.** "Of these 18, we were unable to locate five families because the templates were not of sufficient quality to accurately pick arrivals" (§3). Of the 13, **only 4 have free depths**: "most locations have poor depth control so we fix their depths to the plate interface. Of the remaining four events … we find depths ranging from 31 to 34 km" (§3) | **HYPOINVERSE** [Klein 2002], run with **two** velocity models (PNSN "O0"; Williams et al. 2011 local 1-D) |
| **Bostock et al. 2015** (SVI magnitudes) | **130** templates ("We employ 130 low-frequency earthquake (LFE) templates", Abstract; "a total of 130 distinct source regions on the plate boundary", §2); 269,586 detections | **INHERITED, not re-located.** The 130 are "those templates from Bostock et al. [2012] and Royer and Bostock [2013], as well as new additions" (§2). This paper states no new location step — **SILENT** on locating | n/a (uses "the 1-D model used for routine earthquake location by the Geological Survey of Canada" only for template alignment, §2 step 1) |
| **Plourde et al. 2015** (S. Cascadia) | **59 templates → 34 distinct families** ("These 59 templates were found to represent 34 distinct LFE families", §3); 9,900 LFEs, mean 290/family, max 870 | **YES — all 34.** "Hypocentral depths range from 28 km to 47 km" (§3); every family is labelled with its depth in Fig. 1 | **SILENT on the software.** Only: "For hypocentral location, we used a one-dimensional velocity model modified from Verdonck and Zandt [1994]" (§2) |
| **Ducellier & Creager 2022** (S. Cascadia, 8 yr) | uses **Plourde's *unreduced* 66 families** ("For each of 66 LFE families they provided hypocentral locations…"); "~\citeA{PLO_2015} have later reduced the number of LFE families to 37 … but we chose to use the initial detections"; **6 dropped → 60 carried** ("We eliminate six families for which too few permanent stations have good templates") | **NO — adopted Plourde's locations.** Reviewer reply: "We use the same locations as P[l]ourde et al. (2015). We have not tried to relocate the families using the new templates, except for family A" | family A only, relocated by S−P and P−P times → "30 to 50 kilometers east … depth is 20-25 kilometers" |

| **Sweet, Creager & Houston 2014** (W. Washington, downdip edge) | **1** family — "We create and analyze a 5 year catalog of **9000** repeats of a **single LFE family**" (§1); final catalogue **8942** events (§2) | **YES.** "we employed a **double-difference relocation method** [Waldhauser and Ellsworth, 2000] to find precise relative locations for the **700 most highly correlating events**" (§4); absolute location from "differential S and P times between stations, and S minus P times at individual stations" | **hypoDD** on ~**4 million** autocorrelation differential times (1.4 M P, 2.6 M S), **73 %** variance reduction, formal errors **100 m** |
| **Shelly et al. 2025** (Mendocino triple junction) | **27** families — "we developed **27 separate waveform templates (families)**"; the number is "**largely arbitrary** … chosen as sufficient to collectively detect activity during the 2018–2024 episodes of tremor identified in this area by PNSN" (§2) | **YES — all 27.** "we determined initial locations for each family using a **grid search for minimum mean residual** … We refined these locations with **hypoDD**" (§2) | **grid search in the layered P-wave model of Yoon and Shelly (2024)**, then **hypoDD** [Waldhauser & Ellsworth 2000] on "differential times derived **both from phase pick times and from waveform cross-correlation**", relocated jointly with regional seismicity 2015-01-01→2024-10-11 from **NCEDC** picks |

> **Read this before any Washington run.** Royer & Bostock located **122 NW + 54 SW + 47 NVI**
> templates with Hyp2000 *and* relocated them with hypoDD (Table S1 parameters). The orchestrator's
> prior belief that they did not locate families is wrong.

---

## 1. Parameter-by-parameter table

Columns: R&B14 = Royer & Bostock 2014 · B12 = Bostock 2012 · C&C17 = Chestler & Creager 2017 ·
S25 = Shelly et al. 2025 · SW14 = Sweet, Creager & Houston 2014 ·
T&B15 = Thomas & Bostock 2015 · B15 = Bostock 2015 · P15 = Plourde 2015 · D&C22 = Ducellier &
Creager 2022 · **OURS** = this repo, with file:line.

### 1.1 Stations

| | Value (verbatim) |
|---|---|
| **R&B14 — autocorrelation** | NW: "we applied network autocorrelation to 4 tremor episodes … using the **7 three-component stations B001, SQM, BS11, W020, W040, GNW and PL11**". SW: "we performed network autocorrelation using 2 tremor episodes … and employed the **7 three-component stations N050, S030, W070, N060, PL11, S040 and N070**". NVI: "applied to two tremor episodes … for **different combinations of 7 of the 10 three-component stations VI10, VW03, VI11, VW02, VI52, VI53, VI05, VI08, VI06, VI04**" (§3.1) |
| **R&B14 — network cross-correlation** | "choosing **subsets of 7 to 10 high SNR, three-component stations available within a template** and scanning through all available tremor dates" (§3.3). Selection criterion = high SNR + availability |
| **R&B14 — data pool** | NW "comprises **25 stations** skirting the eastern flanks of the Olympic Peninsula"; SW "includes **25 stations** running approximately east–west to the southwest of Puget Sound" (§2.1); NVI = "**27 broadband seismometers** deployed along two mutually perpendicular arms" (§2.2); SVI = "a line of some **10 stations** … but also incorporated an additional **>20 stations**" (§2.2) |
| **B12 — autocorrelation ("anchor")** | "seismograms recorded during the 2004 (July 4–26) and 2005 (September 3–25) ETS episodes from POLARIS-BC stations **TWBB, TSJB, TWKB, MGCB, KELB**, and permanent stations **PFB, LZB and PGC**"; the text then calls these "the **7** three-component anchor stations listed above" (§2 ¶4, §2.1 ¶5). ⚠ **The paper lists 8 codes but says 7 — an internal inconsistency in the paper itself.** "5–10 km station spacing" |
| **B12 — network cross-correlation** | "we extend our data set to other stations shown in Figure 1 on southern Vancouver Island, the northern Olympic Peninsula and the Gulf and San Juan Islands" (§2 ¶4). Scanning uses "**subset templates comprising 6 or 7 stations** … Selection of stations for the subset template is based on the SNR level of the template waveforms and, of course, the availability of station data" (§2.2 ¶7) |
| **C&C17** | "the AofA experiment, which consisted of **eight 1 km aperture arrays, each with 10–20 three-component stations** … deployed from June 2009 to September 2010 and August 2011"; autocorrelation "using **a set of three stations from each array**"; "Templates also include stations from the **CAFE** experiment located near the AofA" (§2.1). Individual station codes: **SILENT** |
| **T&B15** | Networks: "Global Seismic Network, the Plate Boundary Observatory, the Central Oregon Locked Zone Array (COLZA) …, the Flexarray Along Cascadia Experiment, and the Pacific Northwest Seismic Network"; "Data from **all high gain stations operating in west central Oregon** were downloaded from the IRIS DMC" (§2). Candidate-template detection used exactly **3**: "We applied the method of Savard and Bostock (2014) to stations **ALP0, ALS0, and BWO0** … due to their generally high signal-to-noise ratios and proximity to the 2009 SSE" (§2). Fig. 1 labels 18 codes: COR, B030, B031, B032, TAKO, TOLO, EUO, ALP0, ALS0, BWO0, MAP0, TOL0, YAC1, FACK, FACL, FACN, FACO, FACR (networks XA, IU, PB, UW, YW) |
| **B15** | "a template is constructed from all possible station/channels that have recorded signal in any one of the total 10 ETS episodes" (§2 step 5). Named codes: **SILENT in text** (Fig. 1 map only) |
| **P15** | "Data were collected for **25 stations from the Mendocino Experiment**, as well as **7 permanent stations from the Northern California Seismic Network**, for the period of 21 March to 30 April 2008" (§2). Individual codes: **SILENT** |
| **D&C22** | FAME (temporary) + BK + NC + PB. "We used both one-component and three-components seismic stations. Depending on availability, we used channels **BHZ, EHZ, HHZ, or SHZ** … We restricted ourselves to seismic stations **less than 100 kilometers away from the epicenter** of an LFE family". Full list in Tables S1/S2. Released roster = `data/Ducellier/stations_permanent.txt`, **57 stations**; per-family station lists in `families_permanent.txt` (**65 families**) |
| **SW14** | "the **1 km aperture Big Skidder (BS) array of six 3-component stations** … This array is part of the **Cascadia Arrays For EarthScope (CAFE)** experiment", correlated "for **all available channels**" and stacked "across all **18 channels**" (§2). Location added "data from the BS array and from **three isolated stations** with noise levels small enough to record these tiny LFEs" (§4) — Figure 1 labels them **W020, B013, FACB**. Amplitudes used "**six 3-component stations** within the 1 km aperture BS array", "**12 total horizontal channels**" (§5) |
| **S25** | Acknowledgements only: "Stations used in this study are operated by the U.S. Geological Survey (**NC** network, doi: 10.7914/SN/NC), UC Berkeley (**BK** network, doi: 10.7932/BDSN), and the Earthscope Consortium (**PB** network …)". **Codes, count, and 3-C-vs-vertical: SILENT** — Figure 1 draws them as triangles, Figure 2's caption says only "Station and channel codes are indicated at right"; the full station table is in Supporting Information S1, **not on disk**. One channel rule IS stated: "Channels added in the final stacking (**only used in stacking and not for event detection**) are label at right in red" (Fig. 2b) — his detection channel set is a strict subset of his stacking set, by an unstated rule. Future reach: "the stacked waveforms with the **expanded station set** used in this study for hypocentral locations might facilitate monitoring back to at least **2009**, when most of the **PB** network borehole seismometers … became operational" (§4) |
| **OURS — north (SVI)** | 13 `DET_SITES` / **11 `CLEAN11`** = GLBC, TWKB, MGCB, TSJB, TWBB, PGC, SILB, SSIB, LCBC, B003, B005 — `lfe_array_discovery/scripts/k_common.py:47-49`; queues run a 15-station `S15` set, `k_repro3_queue.sh:13` |
| **OURS — Washington** | **52 probe-admitted stations** (`lfe_array_discovery/results/site_probe_2008.csv`, `admit==True`), selected by `n_common.py:171-208` from a quiet-day/tremor contrast probe. Set = PB.B001/B003/B005/B006/B007/B013/B018/B201/B941/B943 + 40 XU/YW codes |
| **OURS — south** | Ducellier's own roster, read from his `stations_permanent.txt` (`lfe_south/scripts/t7_duc_mf.py`, `t6_duc_cache.py`) |
| **OURS — Shelly path** | **23 stations, 270 templates over 27 families** (`lfe_south/results/templates_south_build.csv`): BK.BLAS/BRIC/HALS/JEWT/PRDS/RVIT/SPRL/USAL, NC.KBN/KBSB/KBU/KCR/KIP/KLR/KMR/KPR/KSM/KTJ, PB.B047/B049/B933/B934/B935 — our own choice by hypocentral distance, **not his** (his list is SILENT). Per family the median is **10** stations / 16–20 channels |
| **OURS — Olympic arrays** | Beams built per 1-km array, `lfe_array_discovery/scripts/o_beam.py:49-58` |

### 1.2 Time windows, blocks, and what is correlated with what

| | Value (verbatim) |
|---|---|
| **R&B14** | "All data were divided into **24 hour-long segments**" (§3). Autocorrelation: "Each of these data sets was independently analyzed on an **hour-by-hour basis**. Each hour-long segment was divided into **15 s windows lagged by 0.5 s** to produce a total of **7170 individual windows per hour**. **Each window was correlated with all other windows in the same 1-hour record for an individual (station-component) channel**" (§3.1) |
| **B12** | "All data are organized into **24-hour long files**" (§2 ¶4). "we divide each hour-long segment of data into **6-second long windows with 0.5 s overlap** resulting in a total of **7188 windows per hour**. **Each window is correlated with all other windows within the same hour-long record**"; "Since we restrict our search for repeat events to hour-long segments, we are not able to detect events that occur less than twice per hour" (§2.1 ¶5) |
| **C&C17** | Ghosh beams: "For each array and **every 1 min time window** … **every 50% overlapping 1 min time window**". LFE verification: "we **autocorrelate 15 s windows** of three-component data … **centered on each peak through time**" (§2.1). Block length for the MAD: **SILENT** |
| **T&B15** | "autocorrelating records on each station channel for a subset of quality stations over a specified time period, **in our case 1 h**" (§2). **Window length: SILENT** — "we experimented with window durations and adjusted the timing of windows to capture move out across the network", no value given |
| **B15** | Template waveforms "**6 s in length**"; a "**2.0 s duration 'detection' window (including 0.25 s half-cosine taper at either edge)**" (§2 steps 1, 3) |
| **P15** | "**windows of 60 s** around these LFEs were used as initial templates"; "we used initial hypocentral estimates to create a **station-dependent taper window** — this allowed us to use a **short window of 12 s**" (§2) |
| **D&C22** | Template build: "we downloaded an **80-second-long** seismic waveform starting 10 seconds before the LFE detection time … cut the first and last 10 seconds of data to obtain a **one-minute-long template**". Scan: "**we download one hour of seismic data**. Then for each station and each channel, we cross-correlate the one-hour long signal with the one-minute-long template". Shorter final templates: "we thus reduced the length of the templates to **25 to 40 seconds** … **average duration equal to 30 seconds**". Reviewer reply confirms the block: "We analyzed the seismic data **one hour a time**. For each hour of data, we assumed that there was an LFE each time the average cross-correlation over all stations and channels is higher than 8 times the median absolute deviation. **The threshold is thus different for each hour of data**" |
| **SW14** | "A **15 s window** containing the P and S waves for the reference event was autocorrelated with **12 h of data** for all available channels" (§2). Detection stage: "a stack of fully normalized **running autocorrelations** between each stacked template and its corresponding continuous **5 year long seismogram**", with the MAD computed "**For each day**" (§2). Hop/lag step: **SILENT** |
| **S25** | Template/candidate window: "**Original candidate LFE waveform (6 s)** is highlighted in blue" (Fig. 2a). **Hop and block length: SILENT.** The MAD span is stated only through the threshold: "**9 times daily median absolute deviation**" (Fig. 4a) → **per day** |
| **OURS — autocorrelation** | **`WIN_S = 6.0`** s, **`HOP_S = 0.5`** s, `MIN_LAG_S = 6.0` s — `k_auto.py:42-44` (and identically `lfe_south/scripts/ks_auto.py:42-44`, `o_auto.py` via `KA.*`). Blocks: hour loop, `k_auto.py:174-205`; `--block` default `"hour"` at `k_auto.py:150` |
| **OURS — templates** | Long stack window `LONG_PRE=12.0`, `LONG_LEN=44.0`; template cut **`TPL_LEN = 12.0` s**, `TPL_PRE = 2.0` s — `k_grow.py:31-34`; south `t9_templates.py:29-30` (`--pre 2.0 --len 12.0`) |
| **OURS — Ducellier replication** | **`BLOCK_H = 12`** hours — `lfe_south/scripts/t7_duc_mf.py:26` |
| **OURS — Shelly path** | template **12 s** with **2 s pre-arrival** (`t9_templates.py:29-30`), i.e. 2× Shelly's 6 s; MAD block = the **whole day** (`t11_mf.py:91`), which happens to match his "daily MAD" |
| **OURS — Olympic beams** | `WIN_N = 1200` (60 s @ 20 Hz), `HOP_N = 600` (50 %) — `o_beam.py:49-50` |

### 1.3 Band, sample rate, response removal

| | Value (verbatim) |
|---|---|
| **R&B14** | "band-pass filtered between **1 and 8 Hz**, **resampled to 40 sps**" (§3). Response removal: **SILENT** |
| **B12** | "filtered between **1 and 10 Hz** and **resampled to 40 sps**" (§2 ¶4). Response removal: **SILENT** |
| **C&C17** | Ghosh beams: "filtered vertical-component seismograms between **5 and 9 Hz**". LFE autocorrelation: "**band-passed between 2 and 8 Hz**" (§2.1). Moments: "we **deconvolve our data to obtain displacement** and filter it between **2 and 8 Hz**" (§2.3) — i.e. response removed *for moments only*. Sample rate: **SILENT** |
| **T&B15** | "Data were demeaned, detrended, **band-pass filtered between 1 and 10 Hz, resampled to 40 sps**, and stored in 24 h data files" (§2). Response removal: **SILENT** |
| **B15** | "Template waveforms, filtered between **1–8 Hz**" (§2 step 1). Sample rate: **SILENT** |
| **P15** | "Data were **filtered between 1.5 and 9 Hz**" (§2). Sample rate implied by "**Δt = 0.025 s**" (§2) = **40 sps**. Response removal: **SILENT** |
| **D&C22** | "we detrended the data, tapered the first and last 5 seconds of the data with a Hann window, **removed the instrument response**, **bandpassed filter between 1.5 and 9 Hz**, **resampled the data to 20 Hz**" (§Data) |
| **SW14** | Detection band: **SILENT**. Amplitude measurement: "All stations had instrument responses **deconvolved to displacement**, and were band-pass filtered from **2 to 8 Hz**" (§5). Sample rate: **SILENT** |
| **S25** | "Waveforms are filtered **2−9 Hz**" (Fig. 2 caption). **Sample rate: SILENT. Response removal: SILENT** |
| **OURS — north/WA** | **1–8 Hz, 20 Hz, response NOT removed** — `k_common.py:28` (`FS = 20.0`), `k_common.py:141` ("the frozen 20 Hz 1-8 Hz cache"), chain described at `d_equiv.py:5-8` |
| **OURS — south** | **1.5–9 Hz, 20 Hz, response REMOVED** (`output="VEL"`, `pre_filt=(0.2,0.5,10,15)`, `water_level=80`) — `lfe_south/scripts/t6_duc_cache.py:25,81-82`; `t_common.py:5` (`FS=20.0`) |
| **OURS — Shelly path** | **1–8 Hz, 20 Hz, response NOT removed** — the stage-1 cache `t3_cache.py:17` (`bp(...,f1=1.0,f2=8.0)`), `t_common.py:5`. So our band is **1–8 Hz against his 2–9 Hz** and our rate 20 Hz against **SILENT** |
| **OURS — Olympic beams** | steering band `BAND_STEER = (5.0, 8.0)` — `o_beam.py:52`; beam formed over the cached 1–8 Hz (documented deviation from Ghosh's 5–9, `o_beam.py:18-20`) |

### 1.4 Components and how they are combined

| | Value (verbatim) |
|---|---|
| **R&B14** | Autocorrelation on **3-C**: "**The resulting time-series from all stations and all components corresponding to the same hour were stacked** to create 'network autocorrelation' records" (§3.1). Cross-correlation: "the **summed** network cross-correlation coefficient" (§3.3). Stacking normalisation: "waveforms (**normalized to unit maximum amplitude across three components**)" (§3.3) |
| **B12** | "the results are stored as individual station/component autocorrelation records. **Autocorrelation records from all stations and components corresponding to the same hour are summed** to produce a network autocorrelation record" (§2.1 ¶5) |
| **C&C17** | Beams: **vertical only** ("filtered vertical-component seismograms"). Autocorrelation: **3-C** ("15 s windows of three-component data"); "We **stack the autocorrelation functions for all channels and stations using an nth-root stack where n = 2**" (§2.1) |
| **T&B15** | "autocorrelating records on **each station channel**"; "The resulting time series are then **stacked** to create a network autocorrelation function" (§2). Fig. 2 shows "Horizontal (top row) and vertical (bottom row) component seismograms" → 3-C |
| **B15** | "**P (Z component) and S (N, E components)**"; "note, P and S waves are treated independently" (§2 steps 1-2) |
| **P15** | **SILENT** on component handling |
| **D&C22** | Both 1-C and 3-C stations; "if we **stack the cross-correlation signals for all the channels and all the stations**"; his code takes the **mean** over channels (`find_all_LFEs_parallel.py:321` `meancc`) |
| **SW14** | 3-C throughout: "The resulting **fully normalized autocorrelation functions—one for each channel**—were then combined using a **third-root stack** to identify time windows that correlated best across all **18 channels**" (§2). Detection stage stacks "fully normalized running autocorrelations" (stack type not restated). Amplitudes on **horizontal** channels only, keeping "only channels where the time of the peak amplitude matched to within **0.2 s** of each other" (§5) |
| **S25** | Both phases are used: "LFE analysis provides much higher precision in source location by utilizing **both P- and S-phase arrival times**" (§1); picks are made for "**P (orange) and S (purple)**" on the stacks (Fig. 2b). The statistic is a network "**correlation sum**" (Fig. 4a) — whether it is summed over channels or stations, and any weighting, is **SILENT** |
| **OURS — autocorrelation** | 3-C; **plain MEAN over all channels**: `acc /= float(nch)` — `k_auto.py:100` (nch = 3 × n_stations, `k_auto.py:83`) |
| **OURS — matched filter** | per-station 3-C normalised cc, **network mean** — `k_grow.py:196` (`acc = acc / cnt`) |
| **OURS — south** | mean over channels present in the block — `t7_duc_mf.py:125` |
| **OURS — Shelly path** | per-channel normalised cc **summed at a common lag** — `t11_mf.py:87`. ⚠ **This was a defect, not a choice**: each station's template is cut at its own offset `tt_s + tau_s` (`t9_templates.py:93`), which spans **1.49–5.85 s** within one family, so summing at a common lag stacked stations several coda cycles out of phase and put our detection times on no published time axis. Fixed in the copy `t11_mf_aa.py` (shift each channel by `tt+tau-pre` first); after the shift our peak lands at his published time with **median dt 0.00 s** |

### 1.5 Detection statistic and threshold

| | Value (verbatim) |
|---|---|
| **R&B14 — autocorrelation** | "When the network autocorrelation exceeded **8 times the median absolute deviation** (Shelly et al., 2006), we registered the corresponding window pairs as LFE detections" (§3.1). MAD span: the hour (records are per-hour) |
| **R&B14 — cross-correlation** | "when the **summed network cross-correlation coefficient exceeded 8 times the median absolute deviation, now for a 24-hour period**, we logged detections" (§3.3) — **the MAD span changes from 1 h to 24 h between stages** |
| **B12** | "Detections are logged by searching for windows when the network autocorrelation exceeds a threshold of **8 times the median absolute deviation**" (§2.1 ¶5); scanning: "logging times where this quantity again exceeds **8 times the median absolute deviation**" (§2.2 ¶7) |
| **C&C17** | **Two different multiples, on two different statistics.** Beam peaks: "For each array, j, we select all the peaks, l, that are **larger than 18 times the median absolute deviation (MAD)**". Autocorrelation: "look for times when the **stacked autocorrelation function exceeds 45 times the MAD**" (§2.1). Stack type = **2nd-root**, not linear |
| **T&B15** | "When the network autocorrelation function exceeds a given threshold (typically **eight times the median absolute deviation, 8*MAD**)"; network CC: "times when the network cross-correlation **exceeds eight times the median absolute deviation**" (§2) |
| **B15** | Not a detector. Channel admission: "the **correlation coefficient of the template channel with the final template stack exceeds 0.4**" and "the **energy within the detection window exceeds 70%** of that in the total (6 s long) channel trace"; per-detection: "correlation coefficient with the corresponding template channel **greater than 0.25**" (§2 steps 4-5) |
| **P15** | **SILENT.** No MAD multiple, no threshold value is stated anywhere in the paper |
| **D&C22** | "Whenever the value of the average cross-correlation is higher than a threshold (we chose a threshold equal to **eight times the median absolute deviation (MAD) of the stacked cross-correlation** as was done by Shelly et al. 2007)". **Second-stage filter:** "we kept only LFE detections for which the **product of the stacked cross-correlation value by the number of channels recording at that time is higher than a threshold** … different for each LFE family, but … constant over time for a given LFE family", chosen so that "**half the LFE detections above the threshold in the FAME catalog are also in the network catalog**" (and vice versa). Effect: "Filtering the catalog has removed about **84 %** of the LFEs" (FAME) and "about **63 %**" (networks). His code: `MAD = np.median(np.abs(meancc - np.mean(meancc)))`, `threshold = 8.0`, per **1-hour** block |
| **SW14** | "For each day, we calculate the **median of the absolute deviation (MAD)** of the stacked autocorrelation function. We chose to define a detection as a time when the **stacked autocorrelation function exceeded 10\*MAD**, producing a typical threshold in which the **mean correlation typically exceeded 0.2–0.3**" (§2). Template-finding stage used the **third-root** stack; MAD span = **the day** |
| **S25** | "a **correlation sum detection threshold of 9 times daily median absolute deviation**" (Fig. 4a) — **9×MAD, MAD over the DAY**. No second-stage filter is described. Detection-rate caveat: "Minor variations in detection rate over time may relate to **network performance** (Figure S1)" |
| **OURS — autocorrelation** | `NMAD = 8.0`; `thr = median + 8 × MAD` of the **whole block** (hour), computed from a subsampled pass-1 — `k_auto.py:45,88-108`. Linear stack (mean), no n-th root |
| **OURS — matched filter** | `thr = med + 8 × MAD` **per group, per day** — `k_grow.py:37,203-205` |
| **OURS — Ducellier replication** | `THR = 8.0`, `mad = median(|mcc − mean(mcc)|)` over a **12-hour** block — `t7_duc_mf.py:26,126` |
| **OURS — Shelly path** | `thr`: `nm = (c − median(c))/MAD(c)` over the **day**, detections at `nm ≥ 8` — `t11_mf.py:91-93`. Two mismatches with the paper: the multiple is **8, his is 9**; and the cc itself was **not bounded** — `ncc()` divides by `sd = sqrt(max(M·S2 − S1·S1, 1e-30))` from cumulative sums, and on the zero-filled gaps the 20 Hz cache writes (`t3_cache.py:31-32`) the true window variance is 0, the clamp gives `sd = 1e-15` and the "normalised" cc reaches **2.4 × 10¹⁵**. Measured on family 23, 2022-05-11: BK.RVIT 33 % constant samples, NC.KIP 50 %; summed-cc MAD **25.7** (a 20-channel normalised sum cannot exceed 20); **16.6 %** of samples above 8×MAD; **128,922 detections that one day** for 27 families against Shelly's **71**. Fixed in `t11_mf_aa.py` (a window is live only if it contains a sample-to-sample change and its std exceeds 1e-4 of the day's median): **21 detections** the same day, 11 of them his |
| **OURS — Olympic** | inherits `KA.WIN_S / KA.NMAD` (6 s, 8×MAD, linear) — `o_auto.py:237`. **No 18×MAD beam-peak stage, no 45×MAD, no 2nd-root stack** |

### 1.6 Declustering / culling

| | Value (verbatim) |
|---|---|
| **R&B14** | "This catalogue was then **culled to exclude overlapping detections and retain only those with high SNR**" (§3.1). No numeric rule — the SNR bar is **SILENT** |
| **B12** | **SILENT** (no declustering rule stated) |
| **C&C17** | **SILENT** at detection; deduplication is at family level (§1.8 below) |
| **T&B15** | "Detections are then **culled to eliminate duplicates and retain only high-quality detections** (Royer and Bostock, 2014)" (§2). Numeric rule: **SILENT** |
| **B15** | "culled such that **only the detection with the strongest correlation in any given 4 s interval** (approximate maximum duration of the S wave coda) is retained" (§2) |
| **P15** | **SILENT** |
| **D&C22** | "we kept only LFEs that are separated by at least **one second** and, when two LFEs are separated by less than one second, we keep only the one with the higher value of the stacked cross-correlation". His code: `maxdiff = floor(1/(dt*freq0))` with `freq0 = 1.0` → 20 samples @ 20 Hz = **1 s**, chaining consecutive runs |
| **SW14** | "Because the S minus P time is about 6 s, and there are often spurious detections when the **P wave of the template aligns with the S wave** of an LFE on the seismogram, we require the **time between adjacent LFE detections to be greater than 6 s**" (§2). Then three named post-hoc culls on 9600 detections: "(1) identify and remove events that occur **within 30 s from the time of any nearby earthquakes in the PNSN catalog (60 detections)**; (2) remove **large amplitude events that occur during known periods of logging near our stations (30 days in 2007)**; (3) remove **isolated events that did not occur as part of larger LFE swarms** and which were visually determined to be non-LFEs (**180 detections**)" → **8942** (§2) |
| **S25** | "when an overall catalog is compiled, we **retain only the detection with highest correlation sum within each 4 s window**" (§3) — a **fixed 4 s window** rule, and it doubles as the cross-family dedupe |
| **OURS — north/WA** | `DECL_S = 4.0` s max-pool — `k_auto.py:46,114-122`; `k_grow.py:38,206-211` (cited to Bostock 2015) |
| **OURS — south** | 1 s run-chaining, faithful to his `clean_LFEs` — `t7_duc_mf.py:26,44-52,92` |
| **OURS — Shelly path** | **4 s greedy exclusion radius** (`t11_mf.py:29-35,39`), which is **not** his rule: his is the maximum inside each **fixed 4 s window** (§3). `t11_mf_aa.py --decl bin` implements the fixed-window form |

### 1.7 Grouping detections into families

| | Value (verbatim) |
|---|---|
| **R&B14** | "we proceeded to **combine all channels corresponding to a given detection into a single 'super' trace and cross-correlated all such traces against one another** to determine maximum correlation coefficients and corresponding lags. The correlation coefficients were used to **populate a similarity matrix employed within a hierarchical cluster analysis**, allowing the detection waveforms to be grouped into clusters" (§3.2). **Linkage and cut threshold: SILENT. Minimum cluster size: SILENT** |
| **B12** | "These initial detections were **grouped into 140 families using a hierarchical cluster analysis [Rowe et al., 2002] based on correlation coefficients of composite wave trains constructed by lining up all station/component records** corresponding to an individual detection window" (§2.1 ¶6). **Linkage / cut / min size: SILENT** |
| **C&C17** | "we **group the original LFE and all the correlating time windows into an LFE family, imposing a stricter correlation threshold that varies from 0.05 to 0.20 depending on the number of stations** that are active at the time of the LFE (**0.05, 0.1, 0.14, 0.15, and 0.2 for less than 6 stations, 7–9 stations, 10–12 stations, 13–15 stations, and more than 15 stations**, respectively)" (§2.1). No similarity matrix / linkage — family = the set that correlates with the template |
| **T&B15** | No clustering step at all. "We require a **minimum of 100 detections in the 3-day window** to stack an LFE template" (§2) |
| **B15** | n/a (inherits families) |
| **P15** | No clustering step; candidates were found by "a combination of **autodetection methods [Savard and Bostock, 2015] and visual identification**" (§2) |
| **D&C22** | n/a (inherits Plourde's 66 families) |
| **SW14** | No clustering step — one family, found by eye: "We first identified this LFE through **visual inspection** of the data. We chose **one of the highest signal-to-noise instances** of this LFE and used that reference event to build up a template" (§2) |
| **S25** | No clustering step. "We formed **initial waveform templates by manual selection of waveforms during times of PNSN-detected tremor** in this area … In many cases … **multiple nearby time windows were tried as templates; we retained the one that produced the strongest detections at other times**" (§2) |
| **OURS** | **Connected components of a union-find graph**, NOT hierarchical clustering: nodes = detection window times merged within `MERGE_S = 4.0` s, edges = autocorrelation detection pairs — `k_group.py:28,51-65,81-99`. `--min-members` default **3** (`k_group.py:72`); `--budget` **20000** strongest detections carried (`k_group.py:71`). The docstring is explicit that this is our reading, not the paper's algorithm (`k_group.py:4-11`) |

### 1.8 Template construction, iteration, convergence

| | Value (verbatim) |
|---|---|
| **R&B14** | "Waveforms for all channels … within a given cluster were **shifted and stacked** to produce an initial LFE template" (§3.2). "**New templates were formed by stacking waveforms (normalized to unit maximum amplitude across three components) for all newly registered detections. After several iterations of network cross-correlation and stacking, we obtained final sets of templates**" (§3.3). **Number of iterations: "several" — no number. Convergence rule: SILENT.** "Each template possesses **100's to 1000's of independent detections**" |
| **B12** | "we adopted an **iterative procedure** of alternately i) performing LFE detection by network cross correlation with LFE templates, and ii) waveform stacking to update LFE templates. **After 5 such iterations neither template SNR nor numbers of detections continued to show significant improvement**"; "Most templates comprise between **1000–5000 independent detections** made over the 8 ETS episodes" (§2.2 ¶7) |
| **C&C17** | "**up to 150 of the time windows that correlate best with the initial detection are stacked** to create a template waveform. We found that **stacking more than 150 windows does not result in templates with better SNRs** … **After three iterations**, we group …" (§2.1) |
| **T&B15** | "we use a technique known as **phase-weighted stacking on the intermediate iterations** … In order to refine the templates, we **iteratively apply cross-correlation and stacking until no additional detections are registered (generally 10–15 iterations)**. If additional stations begin to show coherent energy, those stations are included in network correlation upon subsequent iterations" (§2) |
| **B15** | Alignment: "**Fine alignment is achieved through an approach [Bungum and Husebye, 1971] wherein individual template waveforms are iteratively shifted using correlation with the composite waveform stack from a previous iteration until convergence is achieved**" (§2 step 2) |
| **P15** | "**Phase-weighted stacking was employed until the final iteration** of network cross correlation and stacking … **in this study we found linear stacking to generate more detections in later iterations**, once the template had a sufficiently high SNR" (§2). Iteration count: **SILENT** |
| **D&C22** | "we took the **150 LFE detection times with the best cross-correlation value**, we downloaded one minute of seismic data around each detection, and **linearly stacked the waveforms**"; initial templates: "we then **linearly stacked all the waveforms after normalizing each waveform with the root mean square (RMS)**". Iterations: **effectively 2 rounds** (FAME templates → FAME catalog → permanent-network templates) |
| **SW14** | "Templates were built by **linearly stacking the 80 best correlating windows, forcing identical time offsets for each station/channel. This process was repeated three times** to further improve the quality of the templates" (§2) |
| **S25** | "a **multi-step, iterative procedure of cross-correlation of candidate templates, stacking the detections with the highest correlation sum across the network, and picking phase arrivals (P and S) on the stacked waveforms**. The stacked waveforms are then used to construct new templates with higher signal-to-noise ratios, with the process repeating", "based on the approach developed by **Shelly and Hardebeck (2010)**" (§2); "**Iterative detection and stacking (2–3 iterations)** greatly improved the signal-to-noise ratio of the waveform templates" (§2). **Number of windows stacked per iteration, alignment rule, and convergence test: SILENT** ("full methods are available in Supporting Information S1", which is not on disk) |
| **OURS** | `MAX_STACK = 200` best-cc members per iteration (`k_grow.py:39,276`); alignment = Bungum–Husebye re-shift, `ALIGN_LAG = 1.0` s, `ALIGN_ITERS = 2` / `align_rounds=2` (`k_grow.py:35-36,64`); outer `--iters` default **3** (`k_grow.py:234`); **linear stack only, no phase-weighted stacking anywhere**; no convergence test — a fixed iteration count. Station admitted if stack envelope peak ≥ `MIN_PEAK_MAD = 3.0` (`k_grow.py:40,124-130`), family kept if ≥ 4 stations (`k_grow.py:131`) |

### 1.9 Dedupe / uniqueness between families

| | Value (verbatim) |
|---|---|
| **R&B14** | **SILENT** |
| **B12** | Acknowledged but not enforced: "**Although some of the 140 composite waveform templates exhibit overlap in contributing detections, we will consider each template as representing a distinct LFE location**" (§3 ¶8) |
| **C&C17** | "We assess the similarity of two LFE families by the ratio, **Rcom**, of the number of common detections between the two families to the number of LFE detections in the family with the least number of detections. **When two LFE families have Rcom > 0.5, the family with the lower median correlation coefficient is removed** … This reduced the number of LFE families from 57 to 43. After this step, **if there are still LFEs that are included in multiple families, we keep the LFE in the family that it correlates best with and remove it from the others**" (§2.2) |
| **T&B15** | "We **eliminate duplicate families by ensuring that no family shares more than 10 % of detection times, defined as a difference between detection times of <5 s**, with any other family" (§2) |
| **B15** | n/a |
| **P15** | Two independent tests: (1) "the number of shared detections times was computed, **allowing for a lag of up to 25 s** … the templates were considered to originate from the same family if **>10 % of detections were shared**"; (2) "cross-correlation coefficients were computed for each template pair, using only stations for which a phase was picked; again, this allowed for a lag of up to 25 s. **A correlation coefficient >0.3** was considered to suggest templates originated from the same family" (§2) |
| **D&C22** | Did not dedupe — used Plourde's **unreduced 66**, explicitly rejecting his later reduction to 37 |
| **SW14** | n/a (one family). Reported as an observation instead: "we find **no other detectable LFEs within 10 km** of this LFE family" (§6) |
| **S25** | No family is removed — all 27 are kept and the dedupe is at the detection level: "Due to their close locations, the densely distributed families often contain **overlapping detections**, where multiple waveform templates detect the same event in the data. In this case, when an overall catalog is compiled, we **retain only the detection with highest correlation sum within each 4 s window**" (§3) |
| **OURS** | `k_dedupe.py` / `k_dedupe2.py` exist but are **not in the k_auto→k_group→k_grow chain** run by the repro queues. The literal chain carries every connected component ≥3 with no cross-family uniqueness rule |

### 1.10 Location

| | Value (verbatim) |
|---|---|
| **R&B14** | "**Hyp2000** (Klein, 2002)" and "**hypoDD** (Waldhauser, 2001)", "both distributed by the United States Geological Survey" (§5). Inputs: "**a combination of ordinary phase picks from our LFE catalogue and differential travel times from phase correlation of P- and S-waves**". hypoDD parameters in **Table S1**: `IDAT = 3` (catalog + cross-correlation), **5 parameter blocks of `NITER = 8` = 40 iterations total**, `WTCCP/WTCCS` 0.01/0.01 → 1/1, `WTCTP/WTCTS` 1/1 → 0.01/0.01, `WRCC/WDCC` −9 → 6/8, `WRCT/WDCT` −9 → 5/2 → 0.5/2, `DAMP` 50, 50, 40, 40, 40. Velocity model: **named only as "the 1-D, fixed VP/VS model used for LFE location"** (§7) — the model itself is **SILENT**. Uncertainties: "Mean, nominal uncertainties in horizontal and vertical locations as supplied by Hyp2000 … are **1.8 km, 2.16 km (NVI), 1.1 km, 1.9 km (SVI), 2.0 km, 1.3 km (NW), and 3.7 km, 3.5 km (SW)**" |
| **B12** | "These times were inverted for location using the **hypoinverse software** distributed by the United States Geological Survey [Klein, 2002]"; "**averages of approximately 16 P-times and 22 S-times per template**" (§3 ¶8). Velocity model: **SILENT** |
| **C&C17** | "we first produce **P and S wave differential times by cross-correlating stacked templates across different stations** and **S minus P times by cross-correlating S (on horizontal channels) versus P (on vertical channels) at the same stations**. We then perform a **grid search using the 3D velocity model from Preston et al. [2003]** … We **convert this P wave model to an S wave model assuming the earth is a Poisson solid** [i.e. VP/VS = √3 ≈ 1.732]. Finally, we improve relative locations between nearby families using the **double difference method of Waldhauser and Ellsworth [2000]** applied to differential P and S times determined by the autocorrelation of template waveforms" (§2.1) |
| **T&B15** | "the remaining LFE families are located using **HYPOINVERSE** (Klein, 2002) with **two different velocity models. The first is the O0 model, a regional velocity model employed by the PNSN for earthquake location, and the second is a local 1D model employed by Williams et al. (2011) derived from the 2D model of Gerdom et al. (2000)**" (§2). Depths mostly fixed: "most locations have poor depth control so we **fix their depths to the plate interface**" (§3) |
| **B15** | No location. 1-D model used only for initial alignment: "**the 1-D model used for routine earthquake location by the Geological Survey of Canada**" (§2 step 1) |
| **P15** | "we used a **one-dimensional velocity model modified from Verdonck and Zandt [1994]**" (§2). **Software: SILENT. Absolute vs differential: SILENT** |
| **D&C22** | Adopted Plourde's. Family-A relocation only: "We locate the family using **S minus P times for each station, or using P-P times between stations, or using a combination of both**" (reviewer reply) |
| **SW14** | "we **cross correlated stacked templates at different stations to measure differential S and P times between stations, and S minus P times at individual stations**. Using data from the BS array and from **three isolated stations** …, we employed a **double-difference relocation method [Waldhauser and Ellsworth, 2000]** to find precise relative locations for the **700 most highly correlating events** … nearly **4 million autocorrelation-obtained differential times (1.4 million for P waves, 2.6 million for S waves)** … a **variance reduction of 73 %**" (§4). "Individual LFE locations have **formal errors of 100 m**". Velocity model for the absolute location: **SILENT**; depths are reported only against interface models — "LFE depths are **within 1 km of two different plate interface models [Preston et al., 2003; McCrory et al., 2012]** and exhibit a similar dip angle to the east", while "The plate interface model of Audet et al. [2010] is about **5 km shallower** than the LFEs" (§4). **Absolute depth in km and the family's latitude/longitude: SILENT in the text** (Figure 1 map only) |
| **S25** | "we determined initial locations for each family using a **grid search for minimum mean residual using the layered P-wave velocity model from Yoon and Shelly (2024)**, following the approach of **Shelly and Hardebeck (2010)**. We refined these locations with **hypoDD (Waldhauser & Ellsworth, 2000)**, using **differential times derived both from phase pick times and from waveform cross-correlation**. We relocated the LFE families **together with nearby seismicity from 01 January 2015 to 11 October 2024**, using phase pick times obtained from the **NCEDC**" (§2). No 3-D model, deliberately: "LFE source locations could benefit from use of a three-dimensional (3D) velocity model; however, based on experience elsewhere, **the differences in locations would likely be minor, particularly for relative locations**" (§3). VP/VS or the S model: **SILENT** |
| **OURS** | Grid search on a **Crosson 1-D** model with **VP/VS = 1.78** — `uu_common.py:8,177,220`; grids `uu_locate.py:37-41` (`BASE`/`WIDE`/`WA_BASE`/`WA_WIDE`, 0.03° lat × 0.04° lon × 2 km depth); observables = envelope-max S times / AIC onsets / station-pair XC differentials / back-projection (`uu_locate.py:5-8`); station terms leave-one-family-out (`uu_locate.py:165`). Relative pass = our **own damped-LSQR double-difference**, not hypoDD: `LINK_KM=25`, `MAX_NEIGH=12`, `MIN_SHARED=4`, `WIN_S=6.0`, `MAXLAG_S=2.5`, `CC_MIN=0.4`, `ITERS=4`, `DAMP_REL=1.0`, `DAMP_ABS=0.05` — `uu_dd.py:30-39` |

### 1.11 Falsifiable reported results

| | Value (verbatim) |
|---|---|
| **R&B14** | Initial autocorrelation detections: "**4915, 3306 and 4267 initial detections** … from the NW, SW and NVI datasets". Initial templates: "**224, 70 and 54**". Final: "**122, 54 and 47**". Depths: "Hypocentres beneath SVI range between **29 and 40 km** depth whereas those below Washington state fall between **31 and 46 km** depth (NW) and **32 and 49 km** depth (SW)… NVI … **33 and 39 km**". Epicentres "bounded by the **25 and 38 km** slab depth contours for the Audet et al. (2010) model and the **32 and 45 km** contours for the McCrory et al. (2012) model". Tremor episodes used — NW: **Jan 17–30 2007, May 01–23 2008, Aug 07–24 2010, Aug 04–23 2011**; SW: **Jan 14–Feb 01 2007, May 01–25 2008**; NVI: **Sep 05–11 2006, Jun 13–18 2007** |
| **B12** | "an initial suite of **5775 LFEs** detected in 2004 and 2005 at a select set of 7 stations is grouped into **140 families**" (Abstract); "For the **46 days of tremor** analyzed at the 7 stations for 2004 and 2005, we logged **5775 initial LFE detections**"; "The LFEs define a relatively tight surface **between 27 and 37 km depth**" (§3 ¶9); ETS windows **2004 July 4–26** and **2005 Sept 3–25**; scanning years **2003, 2007–2011** |
| **C&C17** | "Initially, we found **57 LFE families** … each with **250–1700 detections** … There are **34,479 LFEs** in total"; final "**32,387 unique LFEs sorted into 43 families**"; moments "**M0 ≈ 1.4 × 10^10 to 1.9 × 10^12 N m (Mw = 0.7–2.1)**", computed for **34,264** LFEs; characteristic moments per family "**0.5 to 2.7 × 10^11 N m (Mw = 1.1–1.5)**"; "characteristic moments for LFE families are on average **40 % smaller during inter-ETS** time periods" |
| **T&B15** | "**18 families** and over **15,000** individual LFE detections" (§3; the Abstract says "nearly 16 thousand events"); "**13 locatable** LFE families" (Fig. 3 caption); four free depths "**31 to 34 km**", and "using the O0 model … the depths … **increases by an average of 4 km**"; **"propagates south at a rate of 5 km/day"**; slip front first seen in the northernmost families "**on August 30th**"; scan span "**46 days of continuous data from August 3 through September 17, 2009**"; min 100 detections in the 3-day window |
| **B15** | "**130** LFE templates … **over 269,000 independent detections** … between 2003 and 2013"; exactly "**269,586 detections**"; "**10 major ETS episodes**"; "b values **≥ ~5**", best power-law "**β = 5.19 ± 0.09**" at "M̌0 = 2.49 × 10^12 N m (M̌W = 2.20)" → "b value ≈ 6.3"; "typically **half the total moment release occurs within the first 12–24 h** of LFE activity" |
| **P15** | "**59 templates** … **34 distinct LFE families**"; "**9900 LFEs**, or a **mean of 290 per family**, with **870** as the largest number of detections in a single family"; "Cascadia LFEs extend from **40°N to 41.8°N**"; "Hypocentral depths range from **28 km to 47 km**"; "**The main front of tremor epicenters migrates along strike at 9 km d⁻¹**"; "one instance of **rapid tremor reversal**, observed to propagate in the opposite direction at **10–20 km h⁻¹**"; 3 southern families on Maacama (39.50 N/123.36 W/**15 km**; 39.56 N/123.33 W/**21 km**) and Bucknell Creek (39.45 N/123.10 W/**16 km**), with **869, 302, 325** repeats over 1 Mar – 30 Apr 2008 |
| **D&C22** | Recovery against Plourde: "**For 19 families, we recovered all the LFE detections** initially present … **For 61 families, we recovered more than 90 %** of the initial LFE detections. We also **added 74 % more LFEs**"; filtering removes "**about 84 %**" (FAME) / "**about 63 %**" (networks); "**Only nine families have four or less seismic stations with good templates**"; "We **eliminate six families**"; family A "40.09N and 34.9km depth"; catalog span **2004–2011**, stopped because "In November 2011, several one-component stations stopped recording" |
| **SW14** | "**The 9000 individual LFEs** from this repeating family cluster into **198 swarms** that recur roughly every week" (Abstract); catalogue **9600 → 8942** events for "**October 2006 to September 2011**" (§2). Swarm = "at least **four detections** such that the largest gap between detections is **less than 3 h**" → **4–268 LFEs per swarm (median 29)**, **150** isolated detections, durations **10 min–12 h (median 1 h)**, "**average recurrence interval of 8 days**" with inter-swarm times "as short as 3 h or as long as 4 weeks" (§3). Time-predictability: cumulative swarm moment vs **time until the next swarm** correlates **0.52 ± 0.03** over 40 swarm definitions, "**p values ranging from 10⁻⁹ to 10⁻¹³**", against **0.16 ± 0.03** for time SINCE the last swarm (§3). Locations: patch "**approximately 2 km long and 500 m wide**", elongated along plate convergence, vertical spread "**<300 m**", core "roughly circular with a **radius of 300 m**" holding **46 %** of the 700 events (§4). Amplitudes "**0.2 to 18 nm**" peak-to-peak displacement; moments "**2.9×10⁹ N m (Mw0.3) to 2.3×10¹¹ N m (Mw1.5)**", 5-year total "**2.6×10¹⁴ N m (Mw3.6)**"; amplitude distribution **exponential** with "**A₀ ≈ 4.5 nm**" and "c₃ … about **10⁴**", power-law fit would need "**b-value (slope) of about 4**"; complete above **Mw 1** (§5). Migration within swarms: "**five where locations migrated SW (updip), two … NE (downdip), and two … no evidence**", "usually had speeds of **1 km/h**" over **1–2 km** (§4). ETS: "During each of the **five ETS events between 2007 and 2011** …, we detect LFE activity within **1–3 days following vigorous updip tremor activity**", and ETS swarms are "**not systematically bigger or different than the other 50 swarms** that occur over the course of each ETS cycle" (§6). LFE locking efficiency "reaches only about **20 %**" in two **50 m radius** areas (§6); stress drops "**20 and 1 kPa**" for radii 100 and 300 m (§6) |
| **S25** | "**27** stacked LFE waveform templates that we use to detect events from **2018 to 2024**" (Abstract); Figure 4a span "**01 January 2018–30 November 2024**". Geometry: "the **27 LFE family hypocenters form an arcuate alignment at 22–29 km depth, extending nearly 20 km along a strike of ∼290° and dipping toward the north-northeast at approximately 45°**" (§3); "the majority of the identified LFE families locate **near the downdip edge and in the middle along strike**, with more distributed sources on the updip, western, and eastern edges" (§3); the zone sits "**∼15 km south of the southern edge of Gorda slab seismicity**" (Abstract) and "offset **50–100 km west of the main tremor band**" (Key Points). Timing: "small bursts of activity every few days in each family, with a **typical recurrence time of only about 2 days**" (§3), restated as "typically every **2–3 days** in each family" (§4); "**migration velocities of ∼40–80 km/hr**" (§3); "**Notably absent from our data set are any families that exhibit large, infrequent bursts of activity associated with episodic slow slip**" (§3). Geodesy: "this initial attempt was **unsuccessful** in resolving any deformation related to the LFE activity" (§2). Why 2018: "The analysis presented here goes back to **2018, due to station availability**" (§4). Released catalogue as it sits on our disk: **27 families / 61,441 detections**, 2018-01-01→2024-11-30, depths **22.2–28.6 km** (`lfe_south/results/SOUTH_TRUTH_INVENTORY.md` §1) |

### 1.12 Data availability

| | Value (verbatim) |
|---|---|
| **R&B14** | Waveforms: **no statement**. Catalogue release: **SILENT** (Appendix A supplementary material contains only Table S1 + Figs S1–S4) |
| **B12** | "We gratefully acknowledge the Geological Survey of Canada and Incorporated Research Institutions for Seismology for ready access to high quality waveform data." Catalogue release: **SILENT** |
| **C&C17** | "**All the waveforms from the Array of Arrays and CAFE experiments can be found at the IRIS Data Management Center.**" Catalogue: only "Table S1 in the supporting information" (characteristic moments per family). Full detection catalogue: **SILENT** |
| **T&B15** | "The facilities of IRIS Data Services … were used for access to waveforms". Catalogue release: **SILENT** |
| **B15** | "We gratefully acknowledge the IRIS Data Management Center and the Geological Survey of Canada for allowing access to the seismological waveform data". Catalogue release: **SILENT** |
| **P15** | "Data were accessed through the Incorporated Research Institutions for Seismology Data Management Center, **doi:10.7914/SN/XQ_2007**". Catalogue release: **SILENT** |
| **D&C22** | "The seismic recordings from the FAME experiment … can be downloaded from the IRIS website [doi:10.7914/sn/xq_2007]. Waveform data from the permanent networks were accessed through the NCEDC, doi:10.7932/NCEDC. … **The Python scripts used to analyze the data and make the figures can be found on the first author's Github account, accessible through Zenodo** [ariane_ducellier_2022_6374239]." — **the only one of the seven with a released code+catalogue**, and it is on disk at `catalogs/ducellier_catalog/` |
| **SW14** | **No data availability section.** Supporting Information is "**Readme, Figures S1–S10**" only; the paper states no catalogue or code release — **SILENT**. (The BS array is CAFE, whose waveforms are at the IRIS DMC per C&C17's statement.) |
| **S25** | "**The full catalog of LFE detections is available in Shelly (2025)**" — a USGS data release, on our disk as `lfe_catalog_repair/results/external/other/LFEcat_MTJ_2018-2024_shelly2025.csv` (27 families, 61,441 rows, CC0). "Waveform, phase, and event catalog data were obtained from the **NCEDC**, ncedc.org, doi: 10.7932/NCEDC and the **EarthScope SAGE** Data Management Center"; GNSS from **GAGE**; the P-wave model "was obtained from Benz (2024)" (doi:10.5066/P14QGYIY). **Templates: not released. Code: SILENT** |

---

## 2. MISMATCH LIST, ranked by likely impact

### M1 — Washington: window length is 6 s, Royer used **15 s** (and 40 sps, not 20)
- Paper: R&B14 §3.1 — "**15 s windows lagged by 0.5 s** … 7170 individual windows per hour"; §3 — "resampled to **40 sps**".
- Ours: `lfe_array_discovery/scripts/k_auto.py:42` `WIN_S = 6.0`; `k_common.py:28` `FS = 20.0`.
- Why it matters most: a 15-s window spans the whole S coda and the 6-s window does not; at 20 Hz
  the 1–8 Hz band is at 2.5× Nyquist with a filter roll-off right at the band edge. The `6.0`
  constant is documented in `k_auto.py:6-8` as coming from *Bostock 2012*, which is a **different
  region and a different paper**. This is the transplant that caused the Washington failure.
- Change: `k_auto.py:42` (and the identical `lfe_south/scripts/ks_auto.py:42`, plus `o_auto.py`
  which imports `KA.WIN_S`). Sample rate: `k_common.py:28` and the cache chain (`d_equiv.py:5-8`).

### M2 — Washington: wrong stations, and 3 of Royer's 7 NW anchors are not even admitted
- Paper: R&B14 §3.1 NW = **B001, SQM, BS11, W020, W040, GNW, PL11**; SW = **N050, S030, W070, N060, PL11, S040, N070**.
- Ours: the 52-station probe-admitted set (`lfe_array_discovery/results/site_probe_2008.csv`,
  selected at `n_common.py:171-208`). Of Royer's NW seven, **only B001, BS11, PL11 and W040 are in
  our set — SQM, GNW and W020 are absent** (GNW was explicitly rejected by our 2008 site probe,
  see `n_report.py:200`). Of the SW seven, **N070 is absent**; the other six are present.
- Why it matters: Royer's autocorrelation is a *sum over 7 channels × 3 components*; substituting
  45 other stations changes the statistic's noise floor and its MAD, so the 8×MAD bar is no longer
  the published bar. Running 51–52 stations is not "more of the same experiment".
- Change: pass an explicit `--stas` for the NW/SW runs (`k_auto.py:148`), do not use `admitted()`.

### M3 — Grouping: connected components vs **hierarchical clustering on a similarity matrix**
- Papers: B12 §2.1 ¶6 "**hierarchical cluster analysis [Rowe et al., 2002] based on correlation
  coefficients of composite wave trains**"; R&B14 §3.2 "**populate a similarity matrix employed
  within a hierarchical cluster analysis**".
- Ours: union-find connected components on the detection-pair graph — `k_group.py:51-65,96-99`,
  with node merge `MERGE_S = 4.0` (`k_group.py:28`) and `--min-members 3` (`k_group.py:72`).
- Why it matters: connected components **chain** — A~B and B~C put A and C in one family even when
  cc(A,C)≈0. Hierarchical clustering with a cut does not. This systematically under-produces
  families (few giant components instead of 140/224) and pollutes every template built from them.
  Our own docstring flags this as an interpretation (`k_group.py:4-11`), but the papers are explicit.
- Change: `k_group.py:96-107` — build the composite "super trace" per detection, cross-correlate
  all pairs, cluster the similarity matrix. Linkage and cut are **SILENT in both papers**, so they
  must be swept, not guessed.

### M4 — Olympics: Chestler's thresholds and stack type are not implemented at all
- Paper: C&C17 §2.1 — beam peaks at "**18 times the MAD**"; autocorrelation of "**15 s windows … 2 and 8 Hz**"; "**nth-root stack where n = 2**"; "**exceeds 45 times the MAD**"; template from "**up to 150**" windows; "**After three iterations**"; family membership at a station-count-dependent cc of **0.05/0.1/0.14/0.15/0.2**.
- Ours: `o_auto.py:237` inherits `KA.WIN_S=6.0`, `KA.NMAD=8.0`, linear stack, band 1–8 Hz;
  `k_grow.py:39` `MAX_STACK=200`; no beam-peak/boxcar stage; no 2nd-root; no station-count cc ladder.
- Why it matters: 8×MAD on a linear stack and 45×MAD on a 2nd-root stack are not comparable
  statistics. Our detector is far looser at the autocorrelation stage and has no beam-peak
  pre-selection at all, so it is not scoring the same experiment against Chestler's 43 families.
- Change: `o_auto.py` (window/band/nmad/root), a new beam-peak stage against `o_beam.py` output,
  `k_grow.py:39` 200→150, `k_grow.py:234` iters 3 (already matches Chestler).

### M4b — Olympics: **Sweet et al. 2014 is the fuller statement of that front end** (new, 2026-09-25)
- Chestler & Creager compress the array stage into one paragraph; Sweet, Creager & Houston 2014
  (same group, same CAFE/AofA data, one family) spell the correlation machinery out:
  **15 s window containing P and S**, autocorrelated against **12 h** of data on **all 18 channels**
  of the 6-station, 1-km-aperture **BS array**; per-channel **fully normalised** autocorrelations
  combined with a **THIRD-root stack** (C&C17 uses a 2nd-root stack — the two papers differ);
  template = **linear stack of the 80 best-correlating windows, forcing identical time offsets for
  each station/channel**, **repeated three times**; detection = **10×MAD of the stacked
  autocorrelation, MAD per DAY**, which lands at "mean correlation typically exceeded **0.2–0.3**";
  minimum separation between detections **6 s, justified as the S−P time** (a P-on-S alias rule);
  then three named culls (±30 s of a PNSN catalogue earthquake, logging days, visually non-LFE
  isolated events).
- Ours (`o_auto.py:237` inheriting `KA.WIN_S=6.0`, `KA.NMAD=8.0`, linear stack, `MAX_STACK=200`,
  `DECL_S=4.0`) matches none of it: **6 s vs 15 s**, **linear vs third-root**, **8×MAD vs 10×MAD**,
  **200 vs 80** stacked windows, **4 s vs 6 s** minimum separation, and no forced identical offsets.
- Why it matters: Sweet's 6-s minimum separation is derived from the **S−P time**, i.e. it exists to
  kill the P-aligned-on-S alias. Our 4 s neither kills that alias nor is justified anywhere. And an
  **n-th-root stack is a different statistic** — its MAD multiple cannot be carried across.
- Change: `o_auto.py` (window 15 s, nth-root, nmad), `k_grow.py:39` 200→80 for the Olympic path,
  `k_auto.py:46` 4 s→S−P-derived, and add the three post-hoc culls as REPORTED columns.

### M5 — Ducellier replication: our MAD block is 12 h, his is 1 h
- Paper: D&C22 reviewer reply — "We analyzed the seismic data **one hour a time**. For each hour of
  data … **The threshold is thus different for each hour of data**"; §Method — "we download **one
  hour** of seismic data". His released code confirms it: `find_all_LFEs_parallel.py:220,227-229`
  loops `for hour in range(0, nhour)` with `nhour = ceil((t2-t1)/3600)` and computes
  `MAD = np.median(np.abs(meancc - np.mean(meancc)))` inside that loop (line 321). The 12-hour
  figure in his `__main__` is only how he *chunked the job*, not the MAD span.
- Ours: `lfe_south/scripts/t7_duc_mf.py:26` `BLOCK_H = 12`, MAD at `t7_duc_mf.py:126` over the
  whole 12-h block.
- Why it matters: a 12-h MAD is dominated by whatever the quietest 6 hours look like, so the 8×MAD
  bar moves; also our "complete gapless block only" channel rule (`t7_duc_mf.py:113`) drops a
  channel for all 12 h if it has one gap, where he would only lose the affected hour. Both push
  recall down — which is consistent with our 0.836 against his catalogue.
- Change: `t7_duc_mf.py:26` → `BLOCK_H = 1`, and re-derive the gapless rule per hour.

### M6 — Declustering: 4 s is Bostock 2015's *magnitude-catalogue cull*, not a detector rule
- Papers: B15 §2 — "culled such that only the detection with the strongest correlation in any
  given 4 s interval … is retained" (applied *after* detection, to the 269,586-event magnitude
  catalogue). R&B14 §3.1 says only "**culled to exclude overlapping detections and retain only
  those with high SNR**" — numerically **SILENT**. B12 and P15 are **SILENT**.
- Ours: `DECL_S = 4.0` applied inside the detector in both stages — `k_auto.py:46`, `k_grow.py:38`.
- Why it matters: 4 s at a 0.5 s hop suppresses 8 neighbouring rows; in a tremor burst where LFEs
  genuinely recur faster than 4 s (exactly Plourde's and Thomas's regime) it deletes real events.
  Also the SNR half of Royer's cull is not implemented at all, so we keep low-SNR detections he
  discarded before clustering.
- Change: `k_auto.py:46`, `k_grow.py:38`; add an SNR cull at `k_auto.py:134-137`.

### M7 — Iteration and stacking: no phase-weighted stacking, fixed iteration count, 200 not 150
- Papers: T&B15 §2 — "**phase-weighted stacking on the intermediate iterations**", "**until no
  additional detections are registered (generally 10–15 iterations)**"; P15 §2 — "**Phase-weighted
  stacking was employed until the final iteration**"; B12 §2.2 — "**After 5 such iterations**";
  C&C17 §2.1 — "**up to 150**", "**After three iterations**"; D&C22 — "**the 150 LFE detection
  times with the best cross-correlation value**"; R&B14 §3.3 — "**several iterations**" (no number).
- Ours: linear stack only; `--iters` default **3** (`k_grow.py:234`); `MAX_STACK = 200`
  (`k_grow.py:39`); no convergence test.
- Why it matters: the five papers use five different iteration counts and two different stack
  types. A single hard-coded 3 cannot reproduce Thomas (10–15, PWS) or Bostock 2012 (5).
- Change: `k_grow.py:39,234` and add a PWS option to `build_templates` (`k_grow.py:50-136`).

### M8 — Dedupe: no cross-family uniqueness rule in the running chain
- Papers: C&C17 `Rcom > 0.5` + best-correlation reassignment; T&B15 ">10 % of detection times
  within <5 s"; P15 ">10 % shared with ±25 s lag" **or** "cc > 0.3". B12 explicitly declines to
  dedupe; R&B14 **SILENT**.
- Ours: `k_dedupe.py` / `k_dedupe2.py` are not wired into `k_repro*_queue.sh`.
- Why it matters: family counts are the headline falsifiable number in three of the seven papers
  (57→43, 59→34, 18 after a 10 % rule). Without the rule our count is not comparable.

### M9 — Band and response: south is 1.5–9 Hz with response removed, north/WA is 1–8 Hz raw counts
- Papers: B12 **1–10 Hz**; R&B14 **1–8 Hz**; T&B15 **1–10 Hz**; C&C17 **5–9 Hz** (beam) / **2–8 Hz**
  (autocorrelation); P15 and D&C22 **1.5–9 Hz**. Response removal is stated **only** by D&C22 (and
  by C&C17 for moments).
- Ours: north/WA 1–8 Hz, no response (`d_equiv.py:5-8`); south 1.5–9 Hz with response
  (`t6_duc_cache.py:25,81-82`).
- Impact: our north band is Royer's, not Bostock 2012's (1–10). Scoring against Bostock's 140
  families at 1–8 Hz is a 2-Hz-narrower detector than his.

### M10 — Location: neither Hyp2000 nor hypoDD nor Preston 3-D is used anywhere
- Papers: R&B14 Hyp2000 + hypoDD (Table S1: IDAT 3, 5×8 = 40 iterations, DAMP 50/50/40/40/40);
  B12 hypoinverse; T&B15 HYPOINVERSE with O0 *or* Williams 1-D; C&C17 3-D Preston grid search +
  Waldhauser & Ellsworth DD, VP/VS = √3; P15 1-D Verdonck & Zandt (software SILENT).
- Ours: Crosson 1-D grid search, **VP/VS = 1.78** (`uu_common.py:177,220`), plus a home-grown
  damped-LSQR DD (`uu_dd.py:30-39`).
- Impact: lower than M1–M5 for *finding* families, but it is decisive for scoring depth ranges
  (R&B14's 31–46 km NW / 32–49 km SW; C&C17's interface; P15's 28–47 km). Our VP/VS 1.78 vs
  Chestler's Poisson 1.732 alone shifts depths by ~1 km per 10 km of S−P.

### M11 — Mendocino: our Shelly run used an **unbounded cc**, no moveout, 8×MAD, 1–8 Hz and 12 s
- Paper (S25): band "**2−9 Hz**" (Fig. 2); candidate window "**6 s**" (Fig. 2a); threshold "**9 times
  daily median absolute deviation**" on a network "**correlation sum**" (Fig. 4a); dedupe "retain only
  the detection with **highest correlation sum within each 4 s window**" (§3); templates from
  "**2–3 iterations**" of detect-stack-repick (§2); **P and S** picked on the stacks (§2).
- Ours (`lfe_south/scripts/t11_mf.py`, `t9_templates.py`): **1–8 Hz** cache, **12 s** template with
  2 s pre, **8×MAD**, **4 s greedy** decluster, **one** stacking pass at his published times, no
  picks — and two outright defects: (a) `ncc()`'s variance clamp turns the cache's **zero-filled
  gaps** into cc values of 10¹⁰–10¹⁵, which is what produced **12,848,734 detections against his
  61,441**; (b) the per-station cc series were **summed without the per-station moveout**
  `tt_s + tau_s` (spread 1.49–5.85 s within one family), so true events did not stack and our times
  sat on no published axis.
- A third fault was in the SCORING, not the detector: `t8_score.py` put **all 61,441** of his events
  in the denominator while the filter ran only on the **385 held-out days**, where only **5,021** of
  his events live — the recall **cap** was **0.082**, so "recall 0.026" was 0.32 of the achievable
  maximum, not 0.026 of it.
- Change: `t11_mf_aa.py` (live-window mask + moveout shift + `--tmin 9` + fixed-4-s-window dedupe +
  `--hp2` for 2–9 Hz), `t9_templates_aa.py --len 6.0` for his window length, `aa_restack.py` for his
  2–3 iterations, `t8_score_aa.py --days-mode source|fam` for the denominator. Measured effect on
  2022-05-11: **128,922 → 21** detections for 27 families, of which **11 match his 71**.

### M12 — What we would still not be doing, even after M11: his template iteration and his stations
- Paper (S25 §2): templates are built by **iterating** detection → stack the highest-correlation-sum
  detections → **re-pick P and S** → re-detect, "**2–3 iterations**", with extra channels added at the
  final stacking step and "**only used in stacking and not for event detection**".
- Ours: a single RMS-normalised stack at his published times, every (family, station) template
  averaged with **unit weight** including ones whose own stack SNR is 0.45 (measured across the 270
  templates: stack S-window-to-pre-window SNR p10 **0.84**, median **1.50**, p90 **3.52**).
- Why it matters: his statistic is a sum over channels that each carry the family; ours averages in
  channels that do not, which lowers the mean-cc signal faster than the MAD — the same arithmetic
  that made Ducellier bed B worse than bed A (`SOUTH_RESULTS.md` §0). The paper's rule for which
  channels enter the detection sum is **SILENT**, so it must be swept and reported, never guessed.

---

## 3. What reproducing each paper actually requires

### Bostock et al. 2012 — southern Vancouver Island (our 58/96 blind)
- **Stations:** exactly the anchors `TWBB, TSJB, TWKB, MGCB, KELB, PFB, LZB, PGC` for the
  autocorrelation (the paper calls them "7" but lists 8 — run both). Scanning with 6–7-station
  subsets chosen by template SNR.
- **Window:** 6 s, 0.5 s overlap, hour blocks, 7188 windows/hour. **1–10 Hz, 40 sps.**
- **Grouping:** hierarchical clustering (Rowe et al. 2002) on the composite-wave-train similarity
  matrix. Linkage/cut SILENT → sweep.
- **Threshold:** 8×MAD, both stages. **No declustering rule is published.**
- **Iterations:** 5, judged by "neither template SNR nor numbers of detections continued to show
  significant improvement".
- **Location:** hypoinverse, ~16 P + 22 S per template. No dedupe (deliberately).
- **Score against:** 5775 initial detections from 46 tremor days in 2004+2005 at the anchors;
  140 families; 1000–5000 detections per template over 8 ETS episodes; depths 27–37 km.

### Royer & Bostock 2014 — NW and SW Washington (the one we got wrong)
- **Stations:** NW autocorrelation = `B001, SQM, BS11, W020, W040, GNW, PL11`. SW = `N050, S030,
  W070, N060, PL11, S040, N070`. NVI = 7 of the 10 `VI10, VW03, VI11, VW02, VI52, VI53, VI05,
  VI08, VI06, VI04`. Cross-correlation = subsets of **7 to 10** high-SNR 3-C stations per template.
  Do **not** substitute probe-admitted stations.
- **Window:** **15 s, lagged 0.5 s, 7170 windows/hour, hour blocks. 1–8 Hz, 40 sps.**
- **Dates:** NW 2007-01-17→30, 2008-05-01→23, 2010-08-07→24, 2011-08-04→23. SW 2007-01-14→02-01,
  2008-05-01→25. NVI 2006-09-05→11, 2007-06-13→18.
- **Grouping:** super-trace cross-correlation → similarity matrix → hierarchical clustering.
- **Threshold:** 8×MAD over the **hour** for autocorrelation; 8×MAD over the **24-hour period** for
  network cross-correlation. Cull overlapping + low-SNR detections.
- **Location:** Hyp2000, then hypoDD with Table S1 (IDAT 3; NITER 8×5 = 40; WTCCP/WTCCS
  0.01→1→1→1→1 / WTCTP/WTCTS 1→1→0.01→0.01→0.01; WRCC/WDCC −9,−9,6/8,6/8,6/8; WRCT/WDCT
  −9,−9/8,5/2,5/2,0.5/2; DAMP 50,50,40,40,40).
- **Score against:** initial detections **4915 / 3306 / 4267** (NW/SW/NVI); initial templates
  **224 / 70 / 54**; final templates **122 / 54 / 47**; hypoDD retention **77 % / 72 % / 91 %**;
  depths **31–46 km (NW)**, **32–49 km (SW)**, **33–39 km (NVI)**; Hyp2000 σ 2.0/1.3 km (NW),
  3.7/3.5 km (SW).

### Chestler & Creager 2017 — Olympic Peninsula (our 61/69 spatial match)
- **Stations:** the 8 AofA 1-km arrays (10–20 3-C stations each), 2009-06→2010-09 and 2011-08,
  plus nearby CAFE stations in the templates. Array codes are SILENT in the paper.
- **Pipeline:** Ghosh multibeam (5–9 Hz, vertical, 1-min 50 %-overlap, slowness grid) →
  peaks > **18×MAD** per array → 2-s boxcar stack across arrays → candidate LFEs →
  **15-s 3-C autocorrelation at 2–8 Hz using 3 stations per array**, **2nd-root stack**,
  **> 45×MAD** → template from **≤150** best windows → **3 iterations** → family membership at cc
  **0.05/0.1/0.14/0.15/0.2** for <6 / 7–9 / 10–12 / 13–15 / >15 stations.
- **Dedupe:** `Rcom > 0.5` removes the lower-median-cc family (57→43), then unique-assign each LFE.
- **Location:** differential P and S times across stations + S−P at each station → grid search in
  the **Preston et al. 2003 3-D** model, VP/VS = √3 → Waldhauser & Ellsworth DD.
- **Score against:** **57 → 43** families; **34,479 → 32,387** LFEs; 250–1700 detections/family;
  Mw 0.7–2.1; characteristic moments 0.5–2.7 × 10^11 N m; inter-ETS Mc 40 % smaller.

### Thomas & Bostock 2015 — central Oregon
- **Key negative result to honour:** "**Network autocorrelation … failed to identify a single LFE
  template**" on 3 days (2009-09-03→05) — because station spacing is ~30 km. Do not expect the
  autocorrelation chain to work here.
- **Stations:** candidate templates from exactly **ALP0, ALS0, BWO0** via the Savard & Bostock
  cross-station method; network correlation then over the 18 stations of Fig. 1 (XA/IU/PB/UW/YW),
  adding stations as coherent energy appears.
- **Window:** SILENT for the correlation window. **1–10 Hz, 40 sps, 24-h files.**
- **Threshold:** 8×MAD. **≥100 detections in the 3-day window** to stack a template.
- **Iterations:** phase-weighted stacking on intermediate iterations, linear at the end, "until no
  additional detections are registered (generally **10–15 iterations**)".
- **Dedupe:** no family shares **>10 %** of detection times (match window **<5 s**).
- **Scan:** 46 days, 2009-08-03 → 2009-09-17.
- **Location:** HYPOINVERSE with the PNSN **O0** model and the **Williams et al. 2011** local 1-D;
  depths fixed to the interface except 4 families.
- **Score against:** **18** families, **>15,000** detections, **13** locatable, **4** free depths at
  **31–34 km** (O0 deeper by ~4 km), slip front reaching the northernmost families **2009-08-30**
  and propagating south at **5 km/day**.

### Bostock et al. 2015 — SVI magnitudes (not a discovery paper)
- **Inputs:** the 130 templates from B12 + R&B14 + new additions; 10 ETS episodes 2003–2013.
- **Windows:** templates 1–8 Hz, 6 s; detection window 2.0 s with 0.25 s half-cosine tapers,
  centred on the envelope maximum of the final stack.
- **Alignment:** 1-D GSC routine-location model for the initial shift, then Bungum & Husebye
  iterated re-shift to convergence, P and S treated independently.
- **Channel gates:** template-vs-stack cc > **0.4** AND ≥ **70 %** of the 6-s trace energy inside
  the 2-s window; per-detection cc > **0.25**.
- **Cull:** strongest correlation per **4 s**.
- **Score against:** **269,586** detections; b ≥ ~5; β = 5.19 ± 0.09 at M̌0 = 2.49 × 10^12 N m;
  half the moment in the first 12–24 h of an ETS.

### Plourde et al. 2015 — southern Cascadia (the seed set behind Ducellier)
- **Stations:** 25 Mendocino Experiment (XQ, doi:10.7914/SN/XQ_2007) + 7 NCSN permanent,
  2008-03-21 → 2008-04-30. Codes SILENT.
- **Band:** 1.5–9 Hz; Δt = 0.025 s (40 sps).
- **Candidates:** Savard & Bostock autodetection **plus visual identification** — there is no
  autocorrelation stage to reproduce. 60-s windows as initial templates.
- **Iteration:** network cross-correlation + stacking; phase-weighted until the final iteration,
  then linear; station-dependent taper allows a **12-s** working window. Iteration count SILENT.
- **Threshold:** **SILENT** — the paper never states a MAD multiple. Do not assume 8.
- **Uniqueness:** shared-detection >10 % (lag ≤25 s) **or** template-pair cc >0.3 → same family.
- **Location:** 1-D model modified from Verdonck & Zandt 1994; software SILENT.
- **Score against:** **59 templates → 34 families**; **9900** LFEs; mean **290**/family, max **870**;
  40°N–41.8°N; depths **28–47 km**; main front **9 km/day**; one RTR at **10–20 km/h**; three
  strike-slip families at 15/21 km (Maacama) and 16 km (Bucknell Creek) with 869/302/325 repeats.

### Sweet, Creager & Houston 2014 — one family at the downdip edge, western Washington
- **Stations:** the **6 three-component stations of the 1-km-aperture Big Skidder (BS) CAFE array**
  (18 channels) for detection; add **W020, B013, FACB** (Fig. 1) for location.
- **Discovery:** visual identification of one high-SNR instance — there is no blind stage to
  reproduce. Reference event → **15 s window containing P and S** → autocorrelate against **12 h** of
  data on all 18 channels → **third-root stack** of the per-channel fully normalised functions.
- **Template:** **linear stack of the 80 best-correlating windows, identical time offsets forced per
  station/channel, three iterations**.
- **Detection:** stack of fully normalised running autocorrelations over the **5-year** record,
  **10×MAD with the MAD computed per day** (mean cc ≈ 0.2–0.3 at that bar), **minimum 6 s between
  adjacent detections because S−P ≈ 6 s**.
- **Culls (report each, with counts):** −60 within 30 s of a PNSN earthquake; −(30 logging days in
  2007); −180 isolated visually-non-LFE events → **9600 → 8942**.
- **Location:** template-vs-template cross-correlation for inter-station differential P and S times
  plus S−P at each station → **hypoDD** on ~4 M differential times for the **700** best events;
  expect **73 %** variance reduction and **100 m** formal errors.
- **Score against:** **1** family, **8942** detections over **2006-10 → 2011-09**; **198** swarms;
  recurrence **8 days**; swarm moment vs time-to-next **r = 0.52 ± 0.03** (p 10⁻⁹–10⁻¹³); patch
  **2 km × 500 m**, **<300 m** thick, **300 m** core with 46 % of events; amplitudes **0.2–18 nm**;
  **Mw 0.3–1.5**, total **Mw 3.6**; exponential amplitudes **A₀ ≈ 4.5 nm**; in-swarm migration
  **1 km/h**; LFE activity **1–3 days after** updip ETS tremor in all five 2007–2011 episodes.
- **Falsifiable behaviour that needs no catalogue:** the moment-vs-time-to-next-swarm correlation,
  the 8-day recurrence, the ≤20 % locking efficiency, and the 1 km/h in-swarm migration are all
  testable on any family we grow in that area.

### Shelly et al. 2025 — Mendocino triple junction (our 12.8 M-detection failure)
- **Stations:** **NC + BK + PB**; individual codes, count and component mix are **SILENT** (the table
  is in Supporting Information S1, which is not on disk). Our own 23-station set stands in.
- **Band / window:** **2–9 Hz**; candidate window **6 s**; sample rate **SILENT**.
- **Templates:** manual candidate selection during PNSN-detected tremor, several nearby windows
  tried and the strongest kept; then **2–3 iterations** of detect → stack the highest-correlation-sum
  detections → re-pick **P and S** on the stack. Extra channels may join the **final stack only** and
  are "not [used] for event detection".
- **Detection:** network **correlation sum ≥ 9 × the daily MAD**; then **keep only the strongest
  detection in each fixed 4 s window**, which is also the only dedupe between families.
- **Location:** grid search for minimum mean residual in the **Yoon & Shelly (2024)** layered P model
  → **hypoDD** with catalogue-pick *and* cross-correlation differential times, jointly with regional
  seismicity 2015-01-01 → 2024-10-11 (**NCEDC** picks). No 3-D model, deliberately.
- **Score against:** **27** families, **61,441** detections, **2018-01-01 → 2024-11-30**; depths
  **22–29 km**; an arcuate alignment **~20 km long**, strike **~290°**, dipping **NNE at ~45°**,
  **~15 km south** of the southern edge of Gorda slab seismicity and **50–100 km west** of the main
  tremor band; per-family recurrence **~2 days** (2–3 days) with **no** ETS-like large bursts;
  intra-burst migration **40–80 km/h**.
- **What his paper does NOT give us:** templates (not released), station codes, sample rate, the
  channel-selection rule for the detection sum, and the iteration's stack size. Those are SILENT and
  must be swept as columns.

### Ducellier & Creager 2022 — southern Cascadia, 8 years (our 0.836 recall)
- **Families:** Plourde's **unreduced 66** (not 34, not 37); **60** survive to the permanent-network
  catalogue. Locations adopted from Plourde unchanged (family A excepted).
- **Stations:** FAME + BK + NC + PB; BHZ/EHZ/HHZ/SHZ; 1-C and 3-C; **< 100 km from the family
  epicentre**. Released roster `data/Ducellier/stations_permanent.txt` (57) and per-family lists in
  `families_permanent.txt` (65 rows, durations 30–40 s).
- **Preprocessing:** detrend → 5-s Hann taper → **remove response** → **1.5–9 Hz** → **20 Hz**.
- **Templates:** FAME stage = **60 s**; permanent stage = best **150** detections, linear stack,
  template shortened to **25–40 s (mean 30 s)**.
- **Detection:** per-channel normalised cc → **mean over channels with a complete record** →
  **8×MAD per ONE-HOUR block** → declusters runs within **1 s** keeping the max.
- **Second filter:** keep only detections where `cc × n_channels > per-family threshold`, the
  threshold set so half the survivors are cross-confirmed by the other catalogue. Removes **84 %**
  (FAME) / **63 %** (networks).
- **Score against:** 19 families with 100 % recovery of Plourde's detections; **61 of 66** families
  with >90 % recovery; +74 % added LFEs; 9 families with ≤4 good permanent-station templates;
  6 families eliminated; catalogue 2004–2011.
