# Method from the papers — Thomas & Bostock (2015) central Oregon

Sources (all on disk):
- **T&B** = Thomas & Bostock 2015, Tectonophysics 658:111, `notes/1-s2.0-S0040195115003807-main.pdf`
- **S&B** = Savard & Bostock **2015**, BSSA 105(4):2128, `notes/oregon_sources/Savard_Bostock_2015_BSSA_cross-station_LFE_detection.pdf`.
  T&B cite it as "Savard and Bostock (2014)", which is the AGU Fall Meeting abstract (T&B reference list). The BSSA paper is the full method.
- **CODE** = Savard's original MATLAB, `notes/oregon_sources/CS-quake-detection/matlab_original/` (github.com/savardge/CS-quake-detection; README: "code from the cross-station correlation detection algorithm ... published in BSSA 2015"). Read-only. The version on GitHub is dated after the paper (its date list is set to Aug 2014 CVI/SeaJade), but it has an `'OR'` region branch and an Oregon station file `Dsta_OR.mat` (42 stations incl. ALP0, BWO0, **no ALS0**).
- **S&P** = Schimmel & Paulssen 1997 GJI (phase-weighted stack), `notes/oregon_sources/Schimmel_Paulssen_1997_GJI_phase-weighted_stacks.pdf`.
- Rubin & Armbruster 2013 (G3), cited by T&B for the idea only, is **not** on disk (Wiley blocked the download; OpenAlex lists it as bronze open access). Armbruster et al. 2014 JGR is on disk. Neither is needed for the recipe, since T&B say they applied S&B's method.

Rule: where CODE is more specific than S&B, **follow CODE** and log the difference (§4). Every choice the papers leave open is listed in §5 and becomes a CLI column. Nothing is tuned to reach "18 families".

## 1. Data (T&B §2)
1. "Data from all high gain stations operating in west central Oregon were downloaded from the IRIS DMC. Data were demeaned, detrended, band-pass filtered between **1 and 10 Hz**, resampled to **40 sps**, and stored in **24 h** data files." → `scripts/or_cache40.py`.
   - T&B don't mention instrument correction. S&B §"Search Window Admission" says "corrected for instrument response". CODE reads pre-built binaries, so we can't tell. We apply **no** response correction. Every later step normalises amplitude, so correction only changes the relative spectral shape between instrument types (a column, §5).
2. Stations: T&B Fig. 1 lists COR, B030–B032, TAKO, TOLO, EUO, ALP0, ALS0, BWO0, MAP0, TOL0, YAC1, FACK, FACL, FACN, FACO, FACR, plus B027/B028/HEBO/... in the text/figures. Our 30 are in `notes/thomas2015_candidate_stations.csv`. FACR is not in our list, and FACM has no data.

## 2. Cross-station detector (S&B §Method + CODE) → `scripts/xs_detect.py`
T&B §2: "We applied the method of Savard and Bostock (2014) to stations **ALP0, ALS0, and BWO0** … due to their generally high signal-to-noise ratios and proximity to the 2009 SSE", on the 3 days **09/03–09/05**.

| step | S&B text | CODE (followed) |
|---|---|---|
| 1 window | "scanning window length to be **18 s with 17.5 s overlap**" (SVI; chosen from network max S-delay vs LFE duration) | `tlen=18`, `tinc=0.5` (main_scc.m:7-8); 7200 windows/h |
| 2 channels | horizontals only; "north and east components are treated as independent traces" | N and E rows (main_scc.m:327) |
| 3 normalise | "subtract the mean and normalize to the maximum absolute amplitude of the two horizontal channels" | per station **per hour** by max(|N|,|E|) (main_scc.m:285-287), then per window detrend, demean, unit max of N,E (cleanwf.m:57-68) |
| 4 spectral admission T1 | drop traces "with more than T1 % of their total energy at frequencies lower than 2 Hz or higher than **5 Hz**", T1 = 70 % | **frequency-weighted** energy Σ S·f: below 2 Hz > 70 % **or** above **6 Hz** > 70 % (cleanwf.m:95,139) |
| 5 spectral admission T2 | energy within **±0.5 Hz** of the dominant frequency > T2 = 70 % | ±**5 bins** (±0.28 Hz for 18 s) of the PSD maximum > **70** % (cleanwf.m:82-89; parameters.txt "Latest" says 85) |
| 6 envelope SNR | — (not in paper) | drop station if Σenv/max(env) < **60** on N or E (cleanwf.m:49-52; `Tlow`) |
| 7 station floor | "a minimum of **four** different stations" | ≥4 at main_scc.m:348, cleanwf.m:163, lfedetect.m:76, lfedetect.m:237 |
| 8 lag limits | precritical-incidence lag window per station pair; same-station N/E pair **0.3 s** | gettwin.m: twin=1.5·D/3.6 km/s; for D<25 km, (25/cos45°−√((25−D)²+25²))/3.6; mccc2 searches ±twin/2 → same-station ±**0.15 s** |
| 9 envelope stage | cc of Hilbert envelopes; keep traces with cc > **Cenv** "typically within 0.4–0.5" | `ThreshCC=0.43` (lfedetect.m:28); a channel is kept if it is in ANY pair > 0.43 |
| 10 align + t0 | MCCC (VanDecar & Crosson 1990) on the kept envelopes, stack, t0 = stack max, mask ±**4 s** | mccc2 → shift.m (FFT shift) → mean envelope max → ±4 s (lfedetect.m:132-139) |
| 11 duplicates | — | windows overlapping a detected window by >4 s: keep the one with the larger envelope-stack max (lfedetect.m:151-206); a detection is written when it leaves a 36-window (18 s) memory |
| 12 waveform stage | largest subset S3 with max σ_i < **σmax = 0.3 s**; Cmin from 1.0 down "by a constant decrement"; lags "within **±2 s**" | ccsearch.m: decrement **0.01**, `twin=2` → lags ±**1 s**; σ from mccc2 eq. 3 with n−2; sigrthresh=0.3 (lfedetect.m:24) |
| 13 times | t_S = t0 + δt1 + δt2, mean when both components | ts = tabs − tdel − mean(tdel0) (lfedetect.m:259-262); hyp2000 weight = fix(σ/0.1) |
| 14 location | HYPOINVERSE, depth fixed 32 km, keep epicentres with horizontal error < 3 km | writeHYPO.m (a fake P 3 s before the earliest S, weight 3); the location/erh filter isn't in CODE |
| 15 P / hypocentre | P on Z by correlating with S windows; mean cc > 0.4, one cc > 0.5, σ < 0.1 s; erh, erz < 3 km | detectP*.m (not ported: T&B do not use detection hypocentres) |

## 3. Templates, growth, dedupe, location (T&B §2–3) → `xs_grow.py`, `xs_scan46.py`, `xs_locate.py`
1. Candidate templates: S&B Discussion: "Using our cross-station detection **time windows as initial templates** in the iterative network-matched detection or stacking method of Bostock et al. (2012)". → each detection is one candidate. Its template is its own 18 s window on the detecting stations (one absolute window, so moveout is kept).
2. "we use network cross-correlation on **3 days** of data to register additional detections defined as times when the network cross-correlation exceeds **eight times the median absolute deviation**" → `ab_grow.matched_filter` (imported, not edited), 8×MAD per 24 h.
3. "Detections are then stacked to create a new template with better SNR (Gibbons and Ringdal, 2006)."
4. "we use … **phase-weighted stacking on the intermediate iterations**" (S&P 1997; Thurber et al. 2014). PWS = linear stack × |mean unit phasor|^ν; S&P use **ν = 2**.
5. "iteratively apply cross-correlation and stacking **until no additional detections are registered (generally 10–15 iterations)**."
6. "If additional stations begin to show coherent energy, those stations are included in network correlation upon subsequent iterations." Fig. 2 shows COR joining after iteration 1.
7. "We require a minimum of **100 detections in the 3-day window** to stack an LFE template."
8. "LFE templates that make the cut are cross-correlated through **46 days** … August 3 through September 17, 2009."
9. "We eliminate duplicate families by ensuring that **no family shares more than 10 % of detection times**, defined as a difference between detection times of **< 5 s**, with any other family."
10. HYPOINVERSE with **two models**: PNSN **O0** and the **Williams et al. (2011)** local 1-D model (from Gerdom et al. 2000). "most locations have poor depth control so we **fix their depths to the plate interface**"; the 4 free depths are 31–34 km (Williams), +4 km on average with O0. 5 of 18 were not locatable because their "templates were not of sufficient quality to accurately pick arrivals".

## 4. Code vs paper (S&B) differences — CODE followed by default
| # | paper | code | column |
|---|---|---|---|
| D1 | T1 band edges 2 / **5** Hz, plain energy | 2 / **6** Hz, **f-weighted** energy | `TF_HIGH`, weighting |
| D2 | T2 ±0.5 Hz | ±5 bins (±0.28 Hz) | T2 width |
| D3 | T2 = 70 | cleanwf 70; parameters.txt "Latest" 85 | T2 |
| D4 | — | envelope ratio Tlow = 60 | on/off |
| D5 | waveform lag ±2 s | ±1 s (`twin=2` halved in mccc2) | CS_TWIN |
| D6 | same-station lag 0.3 s | ±0.15 s | SAME_STA_TWIN |
| D7 | "largest subset S3" | **Q1**: ccsearch never accepts the set containing **every** channel (ccsearch.m:47) | `--allow-full-set` |
| D8 | — | **Q2**: the double-count check runs *before* ccsearch, so a stronger-envelope window that then fails ccsearch still erases the earlier detection | counted (`q2_lost`) |
| D9 | band 1–8 Hz (SVI) | — | T&B use 1–10 Hz; T1's 6 Hz edge was tuned for 1–8 Hz data |
| D10 | response-corrected | unknown | none applied |
| D11 | history on absolute time | CODE mixes hour-relative `tt` across hours (bug) | fixed here |
| D12 | — | the last 18 s of memory at day end are never written | flushed here |

## 5. What the papers leave SILENT (defaults are marked; sweep them later as columns, never tune to 18)
Detector (T&B → S&B on only 3 stations):
- **S1 station floor.** S&B/CODE need ≥4 stations, which three stations can never meet. The only value possible is **3 (all three)** → `--min-sta 3`. T&B don't say how they adapted it.
- **S2 location of 3-station detections.** Three S times with depth fixed give an exactly determined epicentre, so the 3 km erh filter has no meaning. **Not applied**; each detection is a candidate.
- S3 thresholds for Oregon: Cenv (0.4–0.5 "depending on noise level"), T1, T2, σmax. Default: CODE values.
- S4 1–10 Hz band with a T1 edge tuned for 1–8 Hz data (D9).
- S5 instrument correction (D10).

Growth:
- **S6 template window length**: "we experimented with window durations". Default = the **18 s** detection window.
- S7 components in the network correlation: default **Z,N,E** (Fig. 2 shows horizontals and vertical).
- S8 network function: sum vs channel-mean (G&R). Default **mean** (ab_grow's).
- S9 per-sample station floor for the MAD pool: default **3**.
- S10 declustering within a template: default one detection per template length.
- S11 **station-admission rule** for "coherent energy". Default: odd/even half-stack correlation in the template window > that station's own random-time 99th percentile (200 draws × 50 windows per half). The statistic is a column for every pool station at every iteration.
- S12 PWS exponent ν (default 2, S&P) and which iterations count as "intermediate". Default: PWS for every non-final iteration, and the **final template is the linear stack**.
- S13 iteration cap: "generally 10–15". Default **15**, with a converged flag per template.
- S14 convergence per template: "no additional detections" is read as n(it) ≤ n(it−1).
- S15 whether duplicate merging happens *during* growth. T&B place it after the 46-day scan. Default: **after** (`--merge-each-iter` = ablation).
- S16 "eliminate" duplicates vs merge them. Default: **eliminate** (keep the family with more 46-day detections), with `--mode merge` as the column.

Location:
- **S17 the O0 and Williams (2011) layer tables** are NOT on disk and are NEEDED. GSC and Crosson models are available as placeholders.
- S18 picking method (manual in T&B). Default: ab_hypo zero-phase picks (Hilbert-rotated extremum, ratio R ≥ 1.5).
- S19 "poor depth control" criterion and which interface model was used for fixing. Default: write both the free and the fixed solution; the interface is McCrory 2012 (Fig. 1 contours are McCrory 2004).
- S20 Vp/Vs: default 1.73.
