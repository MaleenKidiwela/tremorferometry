# Washington N-Olympic 2010 -- the arrays run AS ARRAYS, scored against Chestler & Creager

Agent O, 2026-09-21.  PLAN.md sec.17: the goal is to **reproduce the published Washington LFE families**, so the 2010 N-Olympic null of `WA_2010_RESULTS.md` / `WA_2011_REFIRE_RESULTS.md` is treated as OUR reproduction failure, to be debugged the way the north was.  The prime suspect was the site treatment: 84 of the 96 admitted stations sit in ten clusters <= 1.7 km across, and the 2010 run gave every one of them an independent vote in the network correlation sum.  Chestler & Creager (2017), Ghosh et al. (2012) and Sweet et al. (2014) used the arrays AS ARRAYS.

## 0. What the papers actually did with the arrays (fetched and quoted)

**Chestler & Creager (2017), *JGR* 122, 3099 (10.1002/2016JB013717)** -- the 43-family Olympic
catalogue.  The methods chapter is reproduced verbatim in Chestler's UW dissertation *Using
Low-Frequency Earthquakes to Investigate Slow Slip* (sec. 2.1, open access, fetched 2026-09-21):

> "We use high signal-to-noise (SNR) level waveforms, calculated by Ghosh et al. (2012) using
> multibeam-backprojection (MBBP), to create a catalog of low-frequency earthquakes.  The data
> were collected during the Array of Arrays (AofA) experiment, which consisted of eight 1-km
> aperture arrays, each with 10-20 3-component stations ... For each array and every 1-minute
> time window Ghosh et al. (2012) filtered vertical-component seismograms between 5-9 Hz.  Then
> they apply a **delay-and-sum beamforming** approach to stack data and find the **vector
> slowness** corresponding to the direction of a plane wave that produces the stack with the most
> power.  Stacking seismograms from each window at the optimal slowness produces a continuous
> time series with excellent SNR.  This is repeated for each array and every 50-percent
> overlapping one-minute time window."

> "We use these high SNR waveforms to search for new LFE families.  For each array, j, we select
> all the peaks that are larger than **18 times the median absolute deviation** ... we create a
> **boxcar-like time series** for each array [+-1 s] ... Finally, we stack the time series for
> all arrays and use the peaks in the resulting function as potential LFEs."

> "To verify that potential LFEs repeat in time, using a set of **3 stations from each array**,
> we autocorrelate **15-second windows** of 3-component data, band-passed between **2 and 8 Hz**,
> centered on each peak through time.  We stack the autocorrelation functions for all channels
> and stations using an **nth-root stack where n=2** and look for times when the stacked
> autocorrelation function exceeds **45 times the median absolute deviation**.  For each
> potential LFE, up to **150** of the time windows that correlate best with the initial detection
> are stacked to create a template waveform ... After **three iterations**, we group the original
> LFE and all the correlating time windows into an LFE family, imposing a stricter correlation
> threshold that varies from 0.05 to 0.20 depending on the number of stations that are active."

> "To locate LFE families we first produce P- and S-wave differential times by cross correlating
> stacked templates across different stations and S minus P times by cross correlating S- (on
> horizontal channels) versus P- (on vertical channels) waves at the same stations.  We then
> perform a grid search using the **3D velocity model from Preston et al. (2003)** ... Finally we
> improve relative locations between nearby families using the **double difference** method of
> Waldhauser and Ellsworth (2000)."

> "We assess the similarity of two LFE families by the ratio R_com of the number of common
> detections ... When two LFE families have **R_com > 0.5**, the family with the lower median
> correlation coefficient is removed ... This reduced the number of LFE families from **57 to
> 43**."

So the published Olympic chain is, in order: (i) a **slowness-grid delay-and-sum BEAM per 1-km
array** (Ghosh et al. 2012, vertical component, 5-9 Hz, 1-minute windows, 50 % overlap);
(ii) 18 x MAD peaks on each beam, combined across arrays by a +-1 s boxcar stack, as the
CANDIDATE LFE times; (iii) a 15-s, 2-8 Hz, 3-component matched filter seeded by those windows on
**3 stations per array**, network-combined by a **2nd-root stack** at **45 x MAD**; (iv) templates
from the best 150 windows, 3 iterations; (v) R_com > 0.5 de-duplication; (vi) location by
differential P/S times in a 3-D model, then double difference.

**What this study can and cannot copy.**  Steps (i)-(iv) are copied here (sec. 1 and sec. 6).
Step (vi) is NOT: this project has no 3-D Olympic model and no P picks on LFE stacks, so the
location stage stays the study's own -- stack-and-pick S on the Crosson 1-D model with intraslab
station terms, a 50-realisation shuffled-station control and the wide-grid `edge` flag
(`n_locate.py`).  Every family here is therefore located less precisely than Chestler's, which is
why the Chestler check (sec. 5) is run at 5 km and 10 km, not at 1 km.


### 0b. The sites, defined by geometry

Single-link clustering of the 96 admitted stations at 1.5 km, partitioned by network (a borehole and a surface sensor 100 m apart are not one array).  Every (cluster, network) sub-group with >= 2 members is ONE site; singletons stay single stations.  **96 stations -> 22 sites (10 arrays covering 84 stations + 12 single stations).**

| site | net | n | aperture_km | lat | lon | members |
|---|---|---|---|---|---|---|
| xBH | XG | 9 | 1.06 | 48 | -123 | BH01,BH03,BH04,BH05,BH06,BH07,BH08,BH09,BH10 |
| xBS | XU | 9 | 1.22 | 48 | -123 | BS01,BS02,BS03,BS04,BS06,BS22,BS25,BS26,BS27 |
| xDR | XG | 10 | 1.17 | 48 | -123 | DR01,DR02,DR03,DR04,DR05,DR06,DR07,DR08,DR10,DR12 |
| xGC | XG | 13 | 1.38 | 47.9 | -123 | GC02,GC03,GC04,GC05,GC06,GC07,GC08,GC09,GC10,GC11,GC12,GC13,GC14 |
| xLC | XG | 13 | 1.37 | 48.1 | -123 | LC02,LC03,LC04,LC05,LC06,LC07,LC08,LC09,LC10,LC11,LC12,LC13,LC14 |
| xPA | XG | 7 | 1.32 | 48.1 | -123 | PA02,PA03,PA04,PA06,PA07,PA09,PA12 |
| xPBB005 | PB | 3 | 0.281 | 48.1 | -124 | B005,B006,B007 |
| xPBB013 | PB | 2 | 0.044 | 47.8 | -123 | B013,B943 |
| xPL | XU | 6 | 1.22 | 47.5 | -123 | PL01,PL02,PL03,PL04,PL05,PL11 |
| xTB | XG | 12 | 1.68 | 48 | -123 | TB01,TB04,TB05,TB06,TB07,TB08,TB09,TB10,TB11,TB12,TB13,TB14 |


Single stations: B001, B003, B014, B017, B018, GNW, LRIV, VGZ, W020, W030, W040, W060


### 0c. Two controls that changed the beam, both measured on this bed

**(i) The slowness search is not identifiable at this aperture.**  Ghosh's beam picks, per 1-minute window, the slowness that maximises 5-8 Hz stack power.  Run here, it raises the median semblance by 1.3-2.8x over the zero-slowness stack -- but it raises it by the SAME factor on the hour-matched QUIET hours, where there is no tremor to beam at.  The gain is the bias of taking a maximum over the slowness grid, not signal: a 1-km array at 6.5 Hz has a main lobe ~1/(fD) = 0.15 s/km wide, i.e. about seven resolution cells in the whole search.  Measured (median over 4 tremor and 4 quiet hours per array):

| site | M | floor | tremor_zero | tremor_best | tremor_s | quiet_zero | quiet_best | quiet_s |
|---|---|---|---|---|---|---|---|---|
| xBH | 9 | 0.111 | 0.111 | 0.202 | 0.19 | 0.106 | 0.15 | 0.29 |
| xBS | 9 | 0.111 | 0.131 | 0.233 | 0.13 | 0.166 | 0.249 | 0.18 |
| xDR | 10 | 0.1 | 0.11 | 0.207 | 0.13 | 0.107 | 0.163 | 0.2 |
| xGC | 13 | 0.077 | 0.073 | 0.11 | 0.3 | 0.095 | 0.133 | 0.25 |
| xLC | 13 | 0.077 | 0.088 | 0.125 | 0.11 | 0.107 | 0.169 | 0.29 |
| xPA | 7 | 0.143 | 0.219 | 0.272 | 0.07 | 0.143 | 0.345 | 0.21 |
| xPBB005 | 3 | 0.333 | 0.271 | 0.495 | 0.37 | 0.249 | 0.478 | 0.37 |
| xPL | 6 | 0.167 | 0.172 | 0.238 | 0.3 | 0.155 | 0.226 | 0.35 |
| xTB | 12 | 0.083 | 0.096 | 0.147 | 0.13 | 0.124 | 0.164 | 0.24 |


The zero-slowness semblance sits at the incoherent floor 1/M in BOTH populations.  A per-window slowness would therefore inject a random +-0.1-0.3 s jitter into the array trace -- fatal for a matched filter, which needs the beam waveform to be the same function of time every day.  **Adopted: the zero-slowness delay-and-sum beam** (the equal-weight mean of the rms-normalised members), which is the same estimator at the only slowness the data support.


**(ii) How coherent IS a 1-km array in the LFE band?**  Catalogued intraslab earthquakes (25-46 km deep, inside the box) are real deep sources at the LFE incidence.  Coherence N x <stack^2>/<trace^2> of their S window across each array (1 = incoherent, M = perfect), with the pre-event noise window as the control:

| band | M | coh_S | coh_noise | n | perfect | incoherent_floor |
|---|---|---|---|---|---|---|
| 1-2 | 7.5 | 1.93 | 1.21 | 70 | 7.5 | 1 |
| 2-3 | 7.5 | 1.71 | 0.961 | 70 | 7.5 | 1 |
| 3-4 | 7.5 | 1.46 | 0.931 | 70 | 7.5 | 1 |
| 4-6 | 7.5 | 1.24 | 0.902 | 70 | 7.5 | 1 |
| 6-8 | 7.5 | 1.06 | 0.876 | 70 | 7.5 | 1 |
| 1-4 | 7.5 | 1.59 | 1.11 | 70 | 7.5 | 1 |
| 1-8 | 7.5 | 1.28 | 1.05 | 70 | 7.5 | 1 |


So the arrays are only ~25 % coherent at 1-2 Hz and essentially incoherent above 4 Hz.  **This is the measurement that matters for the site treatment**: beam-stacking an array destroys most of the per-member information about the signal, while the 2010 run's 13-votes-per-cluster keeps it but also lets the cluster's SHARED noise vote 13 times.  Per-member statics from the same earthquakes do not cross-validate (leave-one-event-out coherence is not improved), so no statics are applied.


| site | n_events | n_members_with_static | static_sd_ms | coh_plain | coh_static_LOO | coh_noise_plain | coh_noise_static_LOO |
|---|---|---|---|---|---|---|---|
| xBH | 2 | 0 |  | 0.891 | 1.06 | 1.28 | 0.933 |
| xBS | 40 | 8 | 7.93 | 1.49 | 1.46 | 1.32 | 1.32 |
| xDR | 14 | 10 | 54.8 | 1.29 | 1.14 | 1.07 | 1.1 |
| xGC | 2 | 0 |  | 1.16 | 0.943 | 1.11 | 1.11 |
| xLC | 17 | 10 | 53.6 | 1.16 | 0.899 | 0.982 | 0.94 |
| xPA | 17 | 7 | 95.1 | 1.32 | 1.01 | 1.06 | 1.01 |
| xPBB005 | 0 | 0 |  |  |  |  |  |
| xPBB013 | 0 | 0 |  |  |  |  |  |
| xPL | 28 | 6 | 10.8 | 1.43 | 1.38 | 1.18 | 1.17 |
| xTB | 23 | 11 | 40.3 | 1.11 | 1.12 | 1.01 | 1.07 |


### 0d. The site treatment, chosen by measurement (not by assumption)

The identical literal autocorrelation (`k_auto`, 6 s, 1-8 Hz, 8 x MAD of the block, |lag| >= 6 s, 4 s decluster) run three ways on the same tremor and quiet hours:

| mode | n_blocks_tremor | det_per_hour_tremor | traces_tremor | thr_sd_tremor | n_blocks_quiet | det_per_hour_quiet | traces_quiet | thr_sd_quiet | tremor_quiet_ratio |
|---|---|---|---|---|---|---|---|---|---|
| member | 12 | 4350 | 96 | 4.917 | 12 | 3140 | 78 | 5.221 | 1.385 |
| beam | 12 | 130.2 | 21 | 5.168 | 12 | 131.4 | 21 | 5.248 | 0.9905 |
| siteavg | 12 | 8633 | 96 | 5.113 | 12 | 545.9 | 78 | 5.218 | 15.81 |


**Chosen: `siteavg`** -- highest tremor/quiet pair-count ratio at 8 x MAD on 12 + 12 hours restricted to the member set present in both populations; the 2010 run's whole-bed value was 4.2x and the north's truth-based value 11x


## 1. Pair stage

```json
{
 "n_blocks_tremor": 454,
 "n_blocks_tremor_clean": 454,
 "n_blocks_tremor_rejected": 0,
 "n_det_tremor": 849471,
 "nmad_in_sd_median": 5.109324251917286,
 "sta_per_block_median": 96.0,
 "det_per_hour_tremor": 1871.0814977973569,
 "nmad": {
  "p50": 8.5495023727417,
  "p90": 9.69809627532959,
  "p99": 12.365299224853516,
  "max": 44.70164108276367
 },
 "top100_nmad_min": 35.52286911010742,
 "top100_hours": 2,
 "top100_days": 2,
 "top500_nmad_min": 20.084474563598633,
 "top500_hours": 8,
 "top500_days": 7,
 "top2000_nmad_min": 16.555816650390625,
 "top2000_hours": 17,
 "top2000_days": 14,
 "top10000_nmad_min": 12.004586219787598,
 "top10000_hours": 80,
 "top10000_days": 24,
 "top50000_nmad_min": 10.086075782775879,
 "top50000_hours": 254,
 "top50000_days": 27,
 "n_blocks_quiet": 40,
 "n_blocks_quiet_clean": 40,
 "n_det_quiet": 9100,
 "det_per_hour_quiet": 227.5,
 "tremor_quiet_ratio": 8.224534056252118,
 "sta_per_block_median_quiet": 78.0,
 "nmad_in_sd_median_quiet": 5.214930701266695,
 "tremor_quiet_ratio_above_nmad2000": 4405286343.612334,
 "pair_lag_in_tremor_hour": 1.0,
 "pair_lag_in_tremor_hour_note": "1.0 by construction: an autocorrelation block is one hour"
}
```


## 2. Screen -- real (`2010`)

| clause | name | alone | cumulative |
|---|---|---|---|
| c_c1 | (c) located on >= 5 stations | 971 | 971 |
| c_c2 | (c) shuffled-station control passes | 186 | 186 |
| c_c3 | (c) INTERIOR to the wide grid | 113 | 37 |
| c_c4 | (c) |dz| <= 8 km of Slab2 | 222 | 14 |
| c_d | (d) >= 3 distinct tremor days | 1326 | 12 |
| c_a | (a) quiet-day FA <= 0.1/tpl-day and >= 3x lower than tremor | 535 | 3 |
| c_b | (b) tremor co-timing > same-day base rate | 1020 | 3 |


quiet-day network FA 6.00989898989899/template-day; tremor-day rate 5.94452571693951/template-day; pass (a)-(d) **3**
; pass (e) stability **0**; CANDIDATES **0**



## 2. Screen -- site-shift null (`2010_null`)

| clause | name | alone | cumulative |
|---|---|---|---|
| c_c1 | (c) located on >= 5 stations | 761 | 761 |
| c_c2 | (c) shuffled-station control passes | 48 | 48 |
| c_c3 | (c) INTERIOR to the wide grid | 67 | 4 |
| c_c4 | (c) |dz| <= 8 km of Slab2 | 167 | 2 |
| c_d | (d) >= 3 distinct tremor days | 1219 | 2 |
| c_a | (a) quiet-day FA <= 0.1/tpl-day and >= 3x lower than tremor | 514 | 1 |
| c_b | (b) tremor co-timing > same-day base rate | 1149 | 1 |


quiet-day network FA 0.33495962314939437/template-day; tremor-day rate 0.5934932937299856/template-day; pass (a)-(d) **1**
; pass (e) stability **0**; CANDIDATES **0**



## 5. The Chestler check

**A. screened families -> nearest published Olympic coordinate**

```json
{
 "n_screened_located": 0,
 "median_km_3d": null,
 "min_km_3d": null,
 "random_median_km": null,
 "random_p05_km": null,
 "frac_random_at_least_as_close": null,
 "n_within_5km": 0,
 "n_within_10km": 0,
 "n_within_20km": 0
}
```


**B. where are THEIR families in OUR output** (any interior-located grown template within 10 km of each of the 69 derived coordinates; the site-shift null answers the same question and is the false-hit rate)

```json
{
 "n_chestler_coords": 69,
 "n_interior_located_real": 113,
 "n_interior_located_null": 67,
 "n_coords_hit_10km_real": 18,
 "n_coords_hit_5km_real": 5,
 "n_coords_hit_10km_null": 14,
 "n_coords_hit_5km_null": 0,
 "median_nearest_km_real": 13.571941731615428,
 "median_nearest_km_null": 14.914688279822334
}
```


| row | lat | lon | dep | nearest_grp | nearest_km | hit_10km | n_within_10km |
|---|---|---|---|---|---|---|---|
| 22 | 47.9 | -123 | 34.2 | 770 | 3.36 | True | 2 |
| 7 | 47.8 | -123 | 27 | 612 | 4.22 | True | 3 |
| 5 | 47.8 | -123 | 27 | 612 | 4.64 | True | 3 |
| 12 | 47.9 | -123 | 35.1 | 770 | 4.88 | True | 2 |
| 24 | 48 | -123 | 35.4 | 770 | 4.99 | True | 2 |
| 4 | 47.8 | -123 | 35 | 524 | 5.61 | True | 2 |
| 27 | 48 | -123 | 36.7 | 544 | 5.85 | True | 2 |
| 56 | 48 | -123 | 39 | 544 | 6.15 | True | 1 |
| 32 | 48 | -123 | 36.3 | 544 | 6.28 | True | 3 |
| 25 | 48 | -123 | 36.8 | 544 | 6.47 | True | 3 |
| 28 | 48 | -123 | 37 | 625 | 6.67 | True | 2 |
| 54 | 48.2 | -123 | 38 | 885 | 6.69 | True | 1 |
| 26 | 48 | -123 | 38.3 | 544 | 7.19 | True | 1 |
| 30 | 48 | -123 | 40 | 544 | 7.65 | True | 1 |
| 3 | 47.7 | -123 | 24 | 276 | 7.72 | True | 2 |
| 20 | 47.9 | -123 | 37.7 | 524 | 9.13 | True | 1 |
| 51 | 47.7 | -123 | 19 | 276 | 9.43 | True | 2 |
| 47 | 47.9 | -123 | 37 | 524 | 9.58 | True | 1 |
| 66 | 47.9 | -123 | 45 | 155 | 10.2 | False | 0 |
| 48 | 48 | -123 | 40 | 561 | 10.8 | False | 0 |
| 17 | 48.2 | -123 | 41 | 885 | 10.9 | False | 0 |
| 10 | 47.9 | -123 | 38.2 | 524 | 11.2 | False | 0 |
| 33 | 47.9 | -123 | 38 | 524 | 11.3 | False | 0 |
| 9 | 47.9 | -123 | 38.1 | 524 | 11.4 | False | 0 |
| 6 | 47.9 | -123 | 37 | 524 | 11.4 | False | 0 |
| 11 | 47.9 | -123 | 43 | 155 | 11.5 | False | 0 |
| 19 | 47.9 | -123 | 38.8 | 524 | 11.5 | False | 0 |
| 18 | 47.9 | -123 | 38.1 | 524 | 12 | False | 0 |
| 52 | 47.9 | -123 | 38 | 561 | 12.1 | False | 0 |
| 45 | 48 | -123 | 40 | 561 | 12.3 | False | 0 |
| 46 | 48 | -123 | 40 | 561 | 12.3 | False | 0 |
| 13 | 48 | -123 | 41 | 561 | 12.7 | False | 0 |
| 65 | 47.9 | -123 | 45 | 155 | 13 | False | 0 |
| 8 | 47.9 | -123 | 41 | 561 | 13.4 | False | 0 |
| 14 | 47.9 | -123 | 40 | 561 | 13.6 | False | 0 |
| 49 | 47.9 | -123 | 39 | 561 | 13.8 | False | 0 |
| 39 | 48.1 | -123 | 41.7 | 544 | 14 | False | 0 |
| 34 | 48.1 | -123 | 42.3 | 544 | 14.1 | False | 0 |
| 58 | 47.9 | -123 | 40 | 524 | 14.3 | False | 0 |
| 62 | 47.9 | -123 | 40 | 524 | 14.3 | False | 0 |
| 38 | 48.1 | -123 | 41.8 | 544 | 14.7 | False | 0 |
| 60 | 48 | -123 | 41 | 561 | 14.9 | False | 0 |
| 21 | 48 | -123 | 40.5 | 561 | 15 | False | 0 |
| 53 | 48 | -123 | 42 | 561 | 15.2 | False | 0 |
| 42 | 48.1 | -123 | 42 | 544 | 15.2 | False | 0 |
| 2 | 47.9 | -123 | 42 | 524 | 15.7 | False | 0 |
| 29 | 48 | -123 | 42.6 | 561 | 16.3 | False | 0 |
| 41 | 48.2 | -123 | 45.2 | 1095 | 16.3 | False | 0 |
| 61 | 48 | -123 | 42 | 561 | 17.1 | False | 0 |
| 15 | 48 | -123 | 42 | 561 | 17.3 | False | 0 |
| 44 | 48.1 | -123 | 41.9 | 544 | 17.4 | False | 0 |
| 36 | 48.1 | -123 | 42 | 544 | 17.4 | False | 0 |
| 40 | 48.2 | -123 | 44.8 | 1095 | 17.7 | False | 0 |
| 37 | 48.3 | -123 | 41 | 1095 | 18 | False | 0 |
| 55 | 48 | -123 | 43 | 561 | 18.2 | False | 0 |
| 50 | 48 | -123 | 42 | 561 | 18.2 | False | 0 |
| 63 | 48 | -123 | 43 | 524 | 18.2 | False | 0 |
| 59 | 48 | -123 | 40 | 561 | 18.7 | False | 0 |
| 43 | 48.1 | -123 | 41.8 | 889 | 18.9 | False | 0 |
| 23 | 48 | -123 | 41.9 | 561 | 19 | False | 0 |
| 31 | 48.1 | -123 | 43 | 561 | 21 | False | 0 |
| 67 | 48 | -123 | 42 | 524 | 21 | False | 0 |
| 16 | 48 | -123 | 40 | 561 | 21.2 | False | 0 |
| 35 | 48.3 | -123 | 49 | 1095 | 21.4 | False | 0 |
| 57 | 48 | -123 | 41 | 561 | 22.1 | False | 0 |
| 64 | 48.3 | -123 | 50 | 1095 | 23.1 | False | 0 |
| 1 | 48.1 | -123 | 44 | 561 | 23.6 | False | 0 |
| 68 | 48.1 | -123 | 44 | 561 | 26.7 | False | 0 |
| 0 | 48.1 | -123 | 42 | 485 | 29.3 | False | 0 |


## 6. Confirmation on the sealed 2011 window

Rule, fixed before the 2011 numbers existed and identical to `WA_2011_REFIRE_RESULTS` sec.5: >= 3 distinct 2011 tremor days, quiet-day FA <= 0.1/template-day on the 7 strictly-quiet days, and an INTERIOR 2011 location within 10 km (3-D) of the 2010 location.  The null templates go through byte-identical code: their confirmation count is the false-confirmation rate.


**real**: ```{"tag": "real", "n_refired": 3, "n_tremor_days_2007": 26, "n_quiet_days_2007": 7, "n_det_2007": 25, "n_det_quiet_2007": 2, "rate_tremor_per_tpl_day": 0.32051282051282054, "rate_quiet_per_tpl_day": 0.09523809523809523, "n_ge3_days": 3, "n_fa_ok": 1, "n_loc_ok": 0, "n_CONFIRMED": 0}```


**null**: ```{"tag": "null", "n_refired": 1, "n_tremor_days_2007": 26, "n_quiet_days_2007": 7, "n_det_2007": 2, "n_det_quiet_2007": 2, "rate_tremor_per_tpl_day": 0.07692307692307693, "rate_quiet_per_tpl_day": 0.2857142857142857, "n_ge3_days": 0, "n_fa_ok": 0, "n_loc_ok": 0, "n_CONFIRMED": 0}```


real:
| grp | lat_2008 | lon_2008 | dep_2008 | n_det_2007 | n_days_2007 | n_det_quiet_2007 | fa_2007_per_tpl_day | c_days | c_fa | lat_2007 | lon_2007 | dep_2007 | edge_2007 | ctrl_pass_2007 | move_km | c_loc | CONFIRMED |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 67 | 47.5 | -122 | 23 | 17 | 9 | 1 | 0.143 | True | False | 47.5 | -123 | 59 | True | False | 40.6 | False | False |
| 328 | 47.3 | -123 | 31 | 4 | 3 | 0 | 0 | True | True | 47.3 | -123 | 23 | False | True | 30.5 | False | False |
| 434 | 48.2 | -124 | 29 | 4 | 3 | 1 | 0.143 | True | False | 47.8 | -123 | 39 | False | True | 60.6 | False | False |


null:
| grp | lat_2008 | lon_2008 | dep_2008 | n_det_2007 | n_days_2007 | n_det_quiet_2007 | fa_2007_per_tpl_day | c_days | c_fa | lat_2007 | lon_2007 | dep_2007 | edge_2007 | ctrl_pass_2007 | move_km | c_loc | CONFIRMED |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 950 | 47.3 | -123 | 29 | 2 | 1 | 2 | 0.286 | False | False |  |  |  |  | False |  | False | False |


## 7. Arm B -- Chestler's OWN detection chain, run literally

```json
{
 "n_seeds": 756,
 "n_keep": 739,
 "iter0": {
  "n_det": 54417,
  "n_templates": 756
 },
 "iter1": {
  "n_det": 143469,
  "n_templates": 756
 },
 "iter2": {
  "n_det": 208058,
  "n_templates": 756
 },
 "iter3": {
  "n_det": 197566,
  "n_templates": 756
 },
 "dedupe": {
  "iteration": 3,
  "n_families_before": 756,
  "n_families_after": 739,
  "rcom": 0.5,
  "n_detections": 197566
 }
}
```


## 8. Files

| what | path |
|---|---|
| site definition | /home/jovyan/tremorferometry/lfe_array_discovery/results/o_oly2010/sites.csv |
| array coherence (intraslab EQ) | /home/jovyan/tremorferometry/lfe_array_discovery/results/o_oly2010/eq_statics_events.csv, eq_statics_validation.csv, array_coherence_bands.csv |
| beam / slowness controls | /home/jovyan/tremorferometry/lfe_array_discovery/results/o_oly2010/smoke_beam.json |
| site-treatment smoke | /home/jovyan/tremorferometry/lfe_array_discovery/results/o_oly2010/mode_smoke.json |
| beam traces | /home/jovyan/tremorferometry/lfe_array_discovery/beamcache/<SITE>/<year>/<doy>.npz |
| real chain | /home/jovyan/tremorferometry/lfe_array_discovery/results/o_oly2010/2010/ |
| site-shift null chain | /home/jovyan/tremorferometry/lfe_array_discovery/results/o_oly2010/2010_null/ |
| Chestler check | /home/jovyan/tremorferometry/lfe_array_discovery/results/o_oly2010/chestler_check*.csv|json |
| 2011 confirmation | /home/jovyan/tremorferometry/lfe_array_discovery/results/o_oly2010/2011/ |
| Arm B (Chestler chain) | /home/jovyan/tremorferometry/lfe_array_discovery/results/o_oly2010/chestler_arm/ |
