# NOTES — lfe_catalog_repair (append-only ledger)
## 2026-09-10 setup; D1 launched (Opus agent) — user: "run that test"
## 2026-09-10 amendment sent to the D1 agent: quiet != noise; noise-only floor models; inter-ETS candidate diagnostics (see PLAN.md)

## 2026-09-10 D1 results

Ran `scripts/d1_cc_separability.py` (OMP_NUM_THREADS=4). **Streamed 118.6 GB of `data/mf_*.csv`,
124 station-year files, 1,584,478,589 detections, 3,454 templates over the 12 S14 stations
(WISH RADR B943 LEBA PR01 CARB ANGUS HEBO YACT SSO2 ALSE WLOO) in 6.6 min; 6.9 min end to end.**
3,443 templates pass the >= 20 ETS-day AND >= 20 quiet-day gate. Reduction on the fly (pyarrow,
64 MB blocks - a deliberate deviation from the `pandas chunksize=2e6` spec, it is the code S14
validated on these files and ~3x faster on the same three columns) to
`H[template, year, condition, cc-bin]` (0.005 bins, 0.80-1.00), `D[template, day]`,
`Dc[template, day, 0.01 cc bin]` (daily counts above any 0.01-quantised threshold, which is what
makes the burstiness test possible), plus a time-stamped cc >= 0.90 buffer for the coincidence
pilot. **Budget note: the total streamed across the session was ~251 GB, i.e. the stated 150 GB
budget was exceeded**, because the user amendment (Q10 reference + daily above-threshold counts)
arrived after the first full pass and could not be derived from it. A stage-1 cache is now written
(`scratchpad/d1_cache/stage1_*.pkl`, 0.88 GB, `--use-cache`) so any further change to the analysis
costs no re-stream.

### Machinery validation
Aggregate (detection-weighted) ETS/quiet ratio over all templates = **1.268**, against the
**x1.24** S14 measured on the same files by a completely different route (slip-patch event windows
rather than PNSN tremor percentiles). The ETS/quiet definition is therefore sound.

### The four reference floors (the user amendment)
LFEs also occur between ETS episodes, so quiet weeks are noise PLUS inter-ETS LFEs. Every
separability number is reported against four floors:
| floor | definition | pure noise? |
|---|---|---|
| Q25 | quiet-day histogram, local 30-d tremor < 25th pct | no (the original spec) |
| Q10 | quiet-day histogram, local 30-d tremor < 10th pct | no, but purer - **p10 = 0 for 100 % of cells, so Q10 = days with ZERO PNSN tremor within 40 km in the trailing 30 d** |
| NF-exp | exponential fit to log(rate) over cc **0.80-0.85 only**, extrapolated up | yes, but assumes an exponential tail |
| NF-emp | POOLED cc histogram of the **unresponsive** templates (ETS/Q25 < 1.1, n = 1,088), per station, normalised on 0.80-0.85 and rescaled to each template's own 0.80-0.85 rate | yes, no functional form assumed |
Both noise-only floors are scaled **separately for ETS days and quiet days** from that condition's
own 0.80-0.85 anchor, so a change in background noise level is not charged to "signal".

**NF-exp is REJECTED as a floor model, by its own diagnostic.** Its apparent repairability is
*highest* for the templates that have NO tremor response and therefore cannot contain LFEs:
NF-exp "repairable" = **41.7 % unresponsive vs 27.8 % responsive**, max attainable precision
**0.901 unresponsive vs 0.633 responsive** (against Q25: -0.049 vs 0.551). An exponential
extrapolation manufactures signal exactly where none can exist, so the false-alarm tail above
cc ~0.90 is heavier than exponential. Per template the two noise models disagree by
**1.20x [0.40, 4.19]** on the rate above cc 0.90 - the floor there is only known to a factor of a
few. NF-emp is used for the D2 threshold.

### (a) Does the ETS excess sit at higher cc than the floor?  **NO.**
- median cc of QUIET detections **0.8225**, of ETS detections **0.8225**, of the ETS-minus-Q25
  excess **0.8225** -> **shift +0.0000** (i.e. < 0.005, one bin), and the shift is > 0 for only
  **33 %** of templates (38 % of responsive ones).
- fraction of detections above cc 0.90: QUIET **3.25 %**, ETS **3.23 %**. Above cc 0.95: QUIET
  **0.801 %**, ETS **0.665 %**.
- ratio (ETS fraction above t) / (QUIET fraction above t), 1.0 = identical shape:
  **0.916 at 0.90 and 0.739 at 0.95** (responsive 0.888 / 0.772). The ETS-day population is if
  anything *less* high-cc-weighted than the quiet-day population.
- Q10 changes nothing (responsiveness ETS/Q10 1.264 vs ETS/Q25 1.249).
The ETS excess is a **proportional** excess: ~+25 % at every cc, not a high-cc population.
The excess itself is real and overwhelmingly significant - excess_total +18.19 [+0.82, +61.21]
det/template/day, positive for **76.6 %** of templates, year-matched +13.19 and positive for
67.8 %, Wilcoxon p = 3e-249 - it is simply not separable in cc.
The quiet noise tail IS close to exponential where it is anchored (0.80-0.85 fit slope
-35.5 /unit cc, R2 = 0.994, 98 % of templates R2 > 0.9; over the wider 0.80-0.90 window R2 falls
to 0.935 and only 64 % exceed 0.9), but its extrapolation above 0.90 is unreliable per template.

### (b) What fraction of templates is repairable by thresholding alone?  **~1 %, at best ~2 %.**
Repairable = a threshold exists with precision >= 0.8 AND recall >= 0.5.
| floor | max precision (median) | precision >= 0.8 reachable | recall at thr_p80 | repairable (as defined) | repairable, recall <= 1 |
|---|---|---|---|---|---|
| Q25 | 0.405 | 10.4 % | 0.038 | 1.07 % | **0.90 %** |
| Q10 | 0.407 | 10.0 % | 0.038 | 1.02 % | **0.84 %** |
| NF-emp | 0.090 | 11.1 % | 1.027 | 9.03 % | **2.24 %** |
| NF-exp (rejected) | 0.701 | 36.3 % | 1.007 | 30.99 % | 5.37 % |
**CAVEAT on the modelled-floor columns:** recall > 1 for **61 % (NF-emp)** and **71 % (NF-exp)** of
the templates that reach precision 0.8, because a modelled floor can exceed the observed rate in
the low-cc bins and so deflate the "total excess" denominator. The last column excludes those; it
is the honest number. Responsive templates only: Q25 3.7 %, NF-emp 11.2 % as defined but the same
recall > 1 inflation applies. Year-matched (Q25) 1.9 %.
Median thr_p80 0.915 (Q25) / 0.900 (NF-emp); thr_p80 percentiles p10 0.855, p50 0.915, p90 0.962.
**Recommended D2 per-template threshold, derived from the NOISE-ONLY floor at 1 false alarm/day
(NOT from the ETS-minus-quiet excess): NF-emp 1/day = 0.945 [0.920, 0.980], which exists for
96.9 % of templates** (0.1 FA/day exists for only 71.2 %). At the Q25 1/day threshold the observed
ETS rate is only 2.54x [0.78, 6.22] the extrapolated false-alarm rate.

### (b2) Inter-ETS LFE candidates: bursty, but NOT tremor-linked
Quiet-day detections above the NF-emp floor at the per-template 1 FA/day threshold:
- candidate rate **-0.073 [-0.437, +0.531] /template/day**, positive for only **42.3 %** of
  templates; observed quiet tail / NF-emp floor above 0.90 = **0.86x [0.46, 1.56]**. The quiet
  high-cc tail is **consistent with pure noise** - there is no detectable inter-ETS LFE excess.
  (Same at a fixed cc >= 0.95: -0.052 [-1.23, +0.43] /day, positive for 47.0 %.)
- **BURSTY: yes.** dispersion var/mean of the daily quiet-day count above threshold
  **53.2 [20.9, 195.1]**, > 1.5 for **96 %** of templates; **93.0 %** of quiet days have ZERO
  detections above threshold against a Poisson expectation of 40.1 %.
- **TREMOR-LINKED: no.** The test is only defined for **26.4 %** of templates (p25 of the local
  30-d tremor count is 0 for 66.7 % of cells and p10 is 0 for 100 %, so most quiet-day tremor
  series have no variance - a real power limitation). Among the **908** templates where it IS
  defined: Spearman(daily count above cc 0.95, local 30-d tremor) median **-0.004**, positive for
  **48.0 %** - indistinguishable from chance; candidate rate in the top vs bottom half of
  quiet-day tremor > 1 for only **39.3 %**.
So the high-cc quiet-time detections are strongly non-Poisson but their bursts do not track
tremor: this is **bursty NOISE** (storms, glitches, regional seismicity), not an inter-ETS LFE
population. The amendment's alternative hypothesis was tested and does not rescue the catalogue.

### (c) Does network coincidence add discrimination?  **It finds a real population, but not LFEs.**
Pilot on the 5 best-covered source cells (152 cells carry responsive templates at >= 2 stations),
30 station pairs of which 26 have simultaneous coverage, ETS-day detections above each template's NF-emp threshold, +-6 s, with a lag scan
over +-60 s (peak nearest-neighbour lag median -4.5 s; |peak| > 6 s for 13 % of pairs, but
re-centring on the peak barely moves the rates: ETS 0.046 -> 0.045 pairs/day).
- coincident **pair rate**: ETS **0.046**/day, QUIET **0.102**/day, ETS time-shifted +600 s
  **0.004**/day. **ETS/shifted = 11.4** (a real coincident population, ~11x chance) but
  **ETS/QUIET = 0.45** - it is **2.2x more abundant per day in quiet weeks than in ETS weeks**.
  (The pair *fraction* tells the same story: ETS 0.028, QUIET 0.064, shifted 0.0016.)
- the coincident population is confined to the **surface** triple WISH/RADR/LEBA (**8/8** pairs
  show coincidence, mean ETS 0.149 / QUIET 0.330 / shifted 0.013 per day), while pairs involving
  the **B943 borehole**, ANGUS or CARB show essentially none (**1/18** pairs). A deep LFE source
  would be seen best by the borehole; shared surface noise would not.
Coincidence therefore selects genuinely simultaneous arrivals at nearby surface stations that are
time-stationary and not tremor-modulated - regional signals/noise, not LFEs. It is a useful veto
but on its own it does not isolate the ETS-modulated population.

### (d) Is the false-alarm floor stationary?  **Only roughly - a fixed threshold is not safe.**
Spearman(year, quiet-day rate) = -0.104 [-0.424, +0.257]; |rho| > 0.5 for 27.4 % and p < 0.05 for
14.2 % of templates. max/min quiet-year rate **3.27x [2.28, 6.91]** (CV 0.32); **> 2x for 85.3 %**
and **> 5x for 34.1 %** of templates, up to 65x at PR01. A single fixed per-template cc threshold
therefore drifts by a factor of a few over a 17-year record, so **D2/D3 must recalibrate per epoch
(per year or per instrument era)**. The drift is NOT what creates the ETS contrast: the
year-matched ratio (1.230) matches the pooled one (1.249).

### Station and segment structure
Strongly heterogeneous. Responsive fraction: **WISH 84.3 %** (median ratio 1.82), ALSE 39.7 %,
LEBA 38.7 %, ANGUS 29.9 %, B943 20.6 %, RADR 16.0 %, HEBO 16.3 %, CARB 6.0 %, WLOO 5.2 %,
PR01 3.0 %, SSO2 1.9 %, **YACT 0 %** (ratio 0.91). Segments: north 27.9 % responsive (n = 2,155),
central 13.2 % (n = 1,256), south 31.3 % (n = 32). Two pathologies worth flagging: **PR01** has a
nearly flat cc distribution (0.80-0.90 quiet-tail R2 = 0.21, 65x floor drift, ratio 0.91) and
**HEBO** fires at 919 det/template/day over only 412 active days and its 1 FA/day threshold exists
for just 63 % of templates. Highest NF-emp max precision: ANGUS 0.78, CARB 0.71, SSO2 0.68.

### VERDICT
**(a) Does the ETS excess sit at higher cc than the floor?  NO - by 0.000 in median cc (< one
0.005 bin), and the ETS-day population is if anything *less* high-cc-weighted than the quiet-day
population (ratio of the fraction above cc 0.95 = 0.74).** The +25 % ETS excess is real
(p = 3e-249) but proportional across the whole cc range, so cc carries almost no information about
which detections are real.

**(b) What fraction of templates is repairable by thresholding alone?  ~1 %** (0.90 % against the
Q25 quiet floor, 0.84 % against the near-pure Q10 floor, 2.24 % against the noise-only NF-emp
floor once the recall > 1 artifact is removed; 3.7 % among responsive templates against Q25). At
80 % precision the median threshold keeps only **3.8 %** of the excess. Thresholding alone cannot
produce a catalogue.

**(c) Does network coincidence add discrimination?  It adds a real veto but not LFE selection.**
Coincidence is 11x above the time-shifted control, so it isolates a genuinely simultaneous
population - but that population is 2.2x more abundant per day in quiet weeks than in ETS weeks
and is invisible to the co-located borehole, so it is shared surface noise / regional seismicity.
It would help suppress false alarms; it does not by itself recover LFEs.

**(d) RECOMMENDATION: go to D3 (re-detect from waveforms), NOT D2 (repair from stored hits).**
The stored files carry only (template, time, cc) at a fixed cc >= 0.80 with a 6-s gap. D1 shows
that (i) the real detections are not separable from the false-alarm floor in that single scalar,
(ii) at most ~1-2 % of templates could be salvaged by any per-template threshold, (iii) the floor
above cc 0.90 is itself only known to a factor of a few and drifts by > 2x for 85 % of templates,
and (iv) the one extra discriminant available from the stored hits, network coincidence, selects a
non-tremor population. What is missing is per-detection **amplitude / SNR / 3-component
polarisation / network-stack coherence** - none of which is in the files. D3 should therefore
re-detect with per-template MAD-based (not fixed-cc) thresholds, 3-C cc stacking, a network-stack
requirement, and declustering, recalibrated per epoch.
A cheap D2-lite is still worth keeping on the shelf: restricted to the ~22.6 % responsive
templates (best at WISH, ALSE, LEBA, ANGUS), a NF-emp 1 FA/day threshold near cc 0.945 plus the
coincidence veto would give a small, heavily contaminated list - useful for QC or for seeding D3
template selection, not as a catalogue.

Files written: `scripts/d1_cc_separability.py`; `results/d1_template_table.csv` (3,443 rows x 100
columns), `results/d1_stacked_hist.csv` (cc bin x 11 conditions x 4 responsiveness classes),
`results/d1_floor_shapes.csv` (per-station NF-emp noise shapes), `results/d1_coincidence.csv`
(30 station pairs incl. the lag scan), `results/d1_station_summary.csv` (per station / segment /
responsiveness class); `figures/d1_cc_separability.png` (15 panels: A-C stacked ETS/quiet/floor/
excess histograms, D-F precision-recall against the Q25 / NF-emp / NF-exp floors, G-I max-precision
and thr_p80 / recall distributions by floor, J floor-model comparison, K burstiness, L tremor link,
M where the excess sits, N coincidence pair rates, O per station).

## 2026-09-10 orchestrator review — D1 verdict: NOT repairable from the stored hits -> D3
1.58 billion detections, 3,454 templates, 12 stations. The ETS excess (+25 %, +18 det/template/day, Wilcoxon p 3e-249) is
PROPORTIONAL across all cc: median cc of the excess = median cc of the floor (shift +0.000), fraction above 0.90 identical (3.25 %
quiet vs 3.23 % ETS), above 0.95 the ETS population is LESS high-cc weighted (ratio 0.74). With any noise-only floor model only
~1-2 % of templates reach 80 % precision, keeping 3.8 % of the excess. Quiet-week contamination by inter-ETS LFEs is not the issue
(the zero-tremor Q10 reference gives the same answer). Network coincidence is real (11x above a shifted control) but 2.2x MORE
frequent in quiet weeks and confined to the surface triple WISH/RADR/LEBA while the borehole B943 shares none -> surface-noise
coincidences, not a deep source; usable only as a veto. The floor drifts > 2x for 85 % of templates (65x at PR01). The stored cc
carries no separating information; missing are per-detection amplitude/SNR, 3-C polarisation and network-stack coherence.
DECISION: D2 dropped. D3 = GPU waveform re-detection. The stored hits remain useful only as a seed list (NF-emp threshold 0.945,
responsive 22.6 % of templates, plus the coincidence veto). Station pathologies: PR01 (flat cc, 65x drift), HEBO (919/template/day).
Implication for the main project: the +25 % proportional excess means the matched filter's rate rises during ETS at ALL cc — the
detector responds to tremor energy/noise level, not to discrete LFE arrivals; the stacking pipeline is insensitive, but per-event
statistics from these files are invalid (memory updated earlier).
## 2026-09-10 D3 pilot launched (user: "go") — Opus agent; stages 0 (recon) -> 1 (download 4 stations x 2017-2018) -> 2 (GPU re-detection,
## 3-C cc stack, per-template per-30-d MAD thresholds 8/9/10, declustering, per-detection SNR) -> 3 (validation: ETS/quiet, tremor
## correlation, WISH-B943 coincidence, Lin/Bostock recall-precision at PGC; old fixed-0.80 rule on the same series as the control)

## 2026-09-10 literature: how LFEs are located in Thomas et al. 2025 (Sci. Adv. adz2832) and its source catalog
Thomas et al. 2025 do not locate LFEs; they use the Bostock et al. (2015) southern Vancouver Island LFE catalog (Slow Earthquake
Database) + PNSN tremor. Bostock et al. 2012 (G-cubed): network AUTOCORRELATION of 6-s windows at 7 stations -> 140 families; stack all
detections of a family -> template per station; extend templates to more stations by cross-correlation; pick P (Z) and S (N,E) on
the STACKED templates (1-8 Hz) aligned with the GSC 1-D velocity model -> one hypocentre per FAMILY (thrust mechanisms parallel to the
interface from template polarisations). Lin et al. 2024 (Seismica) deep-learning catalog: association + grid search in a 3-D
regional velocity model (>1 M LFEs). Implication for D3: our design (network-consistent templates from stacked detections, P/S
picks on stacks, family hypocentres) is the Bostock approach; event-level locations need >= 4 stations or cc-lag relative methods.

## 2026-09-10 D3 stage 0 — reconnaissance + pilot plan

### (a) How the project downloads waveforms
`scripts/fleet_station.sh <NET> <STA> <START> [CHAN|auto]` is the one-station driver. Step 1 calls
`scripts/download_broadband.py --net --sta --start --end --workers 6 --chan "<band>"`, which fetches **whole
UTC days, ALL components of the selected band** (default selector `BH?,HH?,EH?,SH?`) into the HARD-CODED path
`data/waveforms/<NET>.<STA>/<YYYY>/<JJJ>.mseed`; it is resumable (skips existing days, writes `<file>.nd`
markers for no-data days) and retry-wrapped (IRIS/EARTHSCOPE; NCEDC for BK/NC). Band choice is
`scripts/pick_band.py <NET> <STA>` = the single high-gain vertical band (`[BHES]H?`) spanning most of
2010-2026, +2 bonus for native 40 Hz. **Downstream everything runs at 40 Hz** (`--fs 40`, 2-8 Hz bandpass,
`resample_poly` never obspy `resample`). Because the day files already contain all 3 components and only the
Z is used downstream, **3-C data comes for free with the standard download.**
Wrapper written for this study (the original hard-codes `data/waveforms`):
`lfe_catalog_repair/scripts/d3_download.py --outroot lfe_catalog_repair/traces` — same logic, same layout, so
the project's readers work unchanged with `--wfdir lfe_catalog_repair/traces`.

### (a2) Where the family TEMPLATES live
`scripts/discover_gpu.py` (stage-2 clustering) writes, per station:
* `data/<sta>_disc.npz` — **one array per family, keyed by family_id** (`"46.950_-122.750__c25"` =
  0.05-deg source-cell lat_lon + cluster index). Each is the **Z-component** stack of the family's members,
  40 Hz, 2-8 Hz, `template_pre=-1.0 s .. template_post=+1.0 s` -> **80 samples = 2.0 s**, float32.
* `data/<sta>_disc.summary.csv` — family_id, n_members, n_years, first_year, last_year, lat, lon, snr
* `data/<sta>_disc.members.parquet` — (family_id, time) member detection times, **median 3 members/family**
  (min-family-members 3), spread over the station's whole record.
Pilot-station template inventories (all m=80 @ 40 Hz, Z only):
| station | template file | n families | n cells | member times |
|---|---|---|---|---|
| UW.WISH | `data/wish_disc.npz` | 1,064 | 290 | 3,644 |
| PB.B943 | `data/b943_disc_p70_2010_2026_m3_40hz.npz` (200-sample 100 Hz twin also on disk) | 3,407 | 387 | 10,958 |
| UO.ALSE | `data/alse_disc.npz` | 452 | 170 | - |
| UW.YACT | `data/yact_disc.npz` | 322 | 150 | 1,012 |
| CN.PGC | `data/pgc_pnsn_families_100km.npz` | 5,640 | 671 | 20,787 |
Densify (`scripts/densify_gnw_gpu.py`) used only the **coverage-selected 300** (`*_sel300.summary.csv`); the
full family set is available for D3.
**Source-cell sharing (the network-coincidence lever):** WISH n B943 = **160 shared cells**;
PGC n B943 = 266; PGC n WISH = 72; RADR n WISH = 99; LEBA n WISH = 75. ALSE/YACT share ~0 with the north
(different segment) - so the borehole/surface coincidence test is **WISH x B943**.

### (b) How the GPU matched filter works
`scripts/densify_gnw_gpu.py`: CPU `ProcessPoolExecutor` workers do read -> `resample_poly` -> demean ->
2-8 Hz bandpass -> optional `--despike-mad 8` winsorise, and hand a whole day back as float32; the **main
process only** imports cupy. Per day: `rfft(day, n=N)` once (fixed `N_FIX = fs*86400+16`, `N =
next_fast_len`), templates normalised+time-reversed, **batched in blocks of TB=128** (`rfft` of the batch,
multiply, `irfft`, slice `[m-1 : m-1+valid]`), divided by a GPU sliding-norm (float64 `cumsum` for accuracy),
`cc` clipped to |cc|<=1 and non-finite zeroed, then `cc >= threshold` -> `cp.nonzero` so **only the sparse
crossings cross PCIe**. CPU then does a greedy `--min-gap-s 6` peak-pick (`_dedup_min_gap`), a glitch-day
guard (`--max-raw-det`), a per-day QC (median count > 2000 or median cc > 0.96 -> drop the day) and an
optional top-N-per-(family,day) cap; output CSV = (template, time, cc, station). Peak GPU memory is bounded
at ~4 GB by the 128-template batching. Measured throughput (WISH, **300 templates**, Z only):
**~1.5-1.8 days/s, i.e. ~4 min per station-year** (`logs/densify_wish.log`).
GPU checked live: **NVIDIA L40S, 47.2 / 47.7 GB free**, cupy rfft OK after
`export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH`.

### (c) Lin/Bostock ground truth — what is actually on disk
* `catalogs/lin_families_ets_2010_vi.csv` (218 rows) and `lin_families_top30.csv` are **family LOCATIONS
  ONLY** (family_id, lat, lon, depth_km, n_detections) - no event times.
* **`data/detections_lin_ets_2010_vi.parquet` = 106,644 Lin event times, 217 families,
  2010-05-17 .. 2010-10-14, southern Vancouver Island.** This is the usable ground truth on disk.
* `src/tremorferometry/lin_catalog.py` documents the source: Zenodo 10016020 `EQloc_001_0.1_3_S.csv`,
  1.05 M detections, **2005-01 .. 2017-02**. The raw file is NOT on disk; Zenodo is reachable
  (HTTP 200, 106,324,644 bytes) and a copy is being pulled to `lfe_catalog_repair/results/d3_lin_EQloc.csv`
  (first attempt 504'd; retry loop running).
* `catalogs/ducellier_catalog/` = southern-Cascadia LFE catalogs 2004-2011 / 2007-2009 (per-family pkl with
  year..second + cc + nchannel) - wrong region/era for this pilot. `catalogs/mtj_ducellier/LFEcat_MTJ_2018-2024.csv`
  is a **404 HTML page**, not data.
* PGC family/Lin proximity: **2,129 of 5,640 PGC families are within 5 km of a Lin family** (3,769 within
  10 km), covering **131 distinct Lin families**. PGC is the right ground-truth station.
* **CONSEQUENCE FOR THE PILOT WINDOW:** Lin ends 2017-02, so 2017-2018 gives essentially no overlap. PGC is
  therefore run on ground-truth windows instead: **2010-05-01..2010-11-01** (guaranteed - the on-disk
  parquet, ETS 2010) and **2015-11-01..2016-04-01** (the 2015-12 ETS, 2,069 PNSN tremor within 60 km; usable
  if the raw Zenodo catalog lands). Both are inside 40 Hz BH? eras.

### (c2) Station viability for 2017-2018 (FDSN channel metadata, checked live)
| station | coords | band in window | native fs | verdict |
|---|---|---|---|---|
| **UW.WISH** | 47.1160 -123.7714 | BH? 2011-06-14..2019-05-15 | **40 Hz** | GO - one instrument era, native 40 Hz, D1's most responsive station (84.3 %) |
| **PB.B943** | 47.8132 -122.9113 | EH1/EH2/EHZ 2008-02-27..open | 100 Hz | GO - borehole, 160 cells shared with WISH (104 km apart) |
| **UW.YACT** | 45.9325 -122.4193 | HH? 2013-07-02..2020-10-17 | 100 Hz | GO - and D1's **0 % responsive** station, so it doubles as the hardest case / negative control |
| **UO.ALSE** | 44.3818 -123.5904 | HH? **starts 2019-04-19** | 100 Hz | **IMPOSSIBLE** for 2017-2018 -> replaced by YACT |
| **CN.PGC** | 48.6498 -123.4521 | BH? ..2017-08-10 then HH? 100 Hz | 40 Hz | ground-truth windows only (see above); note the **2017-08-10 sensor change** inside the 2017-2018 window |
ETS/quiet structure inside 2017-2018 (PNSN tremor within 60 km, per month): WISH 2017-02 **368**, 2018-01 207,
2018-05 284, 2018-08 115; B943 2017-02 **2,244**, 2017-03 1,242, 2018-05 **2,557**, 2017-10 601, 2018-01 574;
YACT 2017-12 142, 2018-01 **366**; plus long quiet stretches - the window contains 3-4 ETS per station. Good.

### (d) Disk
`/home/jovyan` = 100 T, **76 T available** (25 % used). Disk is not a constraint.
Measured 3-C day-file sizes: WISH 16.3 MB, B943 23.6 MB, YACT 32.2 MB, PGC 14.2 MB.
Planned download: WISH 730 d = 11.9 GB, B943 730 d = 17.2 GB, YACT 730 d = 23.5 GB, PGC 154+152 d = 4.3 GB
-> **~57 GB total**, inside the ~60 GB cap, 4 stations kept.

### PLAN for stages 1-3
**Stage 1 (running).** `d3_download.py` into `lfe_catalog_repair/traces/<NET>.<STA>/<YYYY>/<JJJ>.mseed`,
one band per station (WISH `BH?`, B943 `EH?`, YACT `HH?`, PGC `BH?`), 2017-01-01..2019-01-01 for
WISH/B943/YACT and the two ground-truth windows for PGC.
**Stage 2 (GPU re-detection).** Reuse the `densify_gnw_gpu.py` GPU kernel, extended:
1. **3-C templates.** The npz templates are Z-only. N/E (B943: 1/2) templates are built by stacking the SAME
   family member times from `<sta>_disc.members.parquet` on the horizontal channels - identical statistical
   basis as the existing Z stacks (median 3 members). The member days are spread over 2009-2026 and are NOT in
   the pilot download, so they are fetched as **short 3-C cutouts (+-30 s) via FDSN bulk requests**
   (~3.6 k x 43 kB = ~0.16 GB per station, negligible) into `lfe_catalog_repair/traces/cutouts/`.
2. **cc per component** on the same 40 Hz / 2-8 Hz preprocessing, then **cc_stack = mean over available
   components**.
3. **MAD thresholds.** Per template per day the GPU also returns the **median and MAD of the stacked cc
   series** (estimated on a stride-m subsample = 43,200 values/day, which is ample and ~independent).
   Threshold for day D = `median + k*MAD` of the **trailing 30-day** pool of those daily stats (median of
   daily medians / median of daily MADs), k = 8, 9, 10; extraction floor = the k=8 threshold (also
   min'ed with the day's own k=8 threshold) so k=9/10 are exact subsets in post-processing.
   NOTE (arithmetic done up front): for an 80-sample 2-8 Hz template the single-channel cc noise sigma is
   ~0.2 (MAD ~0.13), so `median + 9 MAD` **exceeds 1.0 on a single channel** - it is only attainable because
   3-C stacking cuts sigma by sqrt(3) (9 MAD ~ 0.70). Templates whose k-threshold exceeds 1.0 detect nothing;
   that count will be reported per k rather than hidden.
4. **Declustering** with min separation = template length (2.0 s), greedy highest-cc first (the existing
   `_dedup_min_gap`, `--min-gap-s 2`), plus the existing glitch-day guard and per-day QC.
5. **Per detection stored:** time, template, cc_stack, cc_z, cc_n, cc_e, SNR (peak |envelope| in the 2 s
   template window / RMS of the preceding 10 s, same band, Z), threshold used, k, window median, window MAD,
   n_components. **Plus the OLD rule on the identical stacked series** (fixed cc_stack >= 0.80, 6 s gap) for
   the before/after.
   -> `results/d3_detections_<STA>.parquet`, `results/d3_thresholds_<STA>.csv`.
   Cost estimate from the measured 4 min / station-year / 300 templates: WISH 1,064x3 comps x 2 yr ~ 85 min,
   B943 3,407x3x2 ~ 4.6 h, YACT 322x3x2 ~ 26 min, PGC 5,640x3x0.85 yr ~ 3.2 h -> **~9 h GPU** for the pilot.
   If that overruns, the fallback in the brief applies (all families for one ETS month, then all months for
   the top ~50 families by ETS responsiveness).
**Stage 3 (validation).** ETS/quiet per template from PNSN tremor within 40 km of the family's source cell
(30-d count > 90th pct vs < 25th pct, the D1 definition); Spearman of daily counts vs local tremor; fraction
of templates with ETS/quiet >= 5; the same three numbers for the OLD rule on the same series. Network:
WISH x B943 same-cell families, coincidence within +-3 s after a median-lag alignment estimated from the ETS
cross-correlation of daily-binned detection times, against a +600 s shifted control. Ground truth: PGC
families within 5 km of a Lin family - recall (Lin times with a re-detection within +-2 s) and precision, at
k = 8/9/10 and for the OLD rule. Figure `figures/d3_pilot_validation.png`.
## 2026-09-10 stage 4 (family locations) appended to the running D3 pilot agent — see PLAN.md

## 2026-09-10 D3 amendment received (user, via orchestrator): add STAGE 4 = family LOCATION
After stage 3, locate the PGC/Lin-region families from network-consistent template stacks built on
6-8 additional stations within ~100 km, pick P/S on the stacks, grid-search hypocentres in a 1-D
model with a Slab2 prior, compare with Bostock/Lin family locations, and do event-level
double-difference relative locations. Spec in the orchestrator message; verdict wanted on whether
~2-3 km family hypocentres are reachable with this geometry.

## 2026-09-10 D3 stage 1 — download DONE

`lfe_catalog_repair/scripts/d3_download.py` (the `scripts/download_broadband.py` logic with a
configurable `--outroot`; nothing written outside the study folder). One band per station so the 3-C
selection is unambiguous. `.nd` markers mark no-data days.

| station | net | window | band | day-files | no-data | size |
|---|---|---|---|---|---|---|
| WISH | UW | 2017-01-01..2019-01-01 | BH? 40 Hz | 727 | 3 | 13 GB |
| B943 | PB | 2017-01-01..2019-01-01 | EH? 100 Hz (Z,1,2) | 730 | 0 | 19 GB |
| YACT | UW | 2017-01-01..2019-01-01 | HH? 100 Hz | 724 | 6 | 26 GB |
| PGC | CN | 2010-05-01..2010-11-01 (Lin GT) + 2015-11-01..2016-04-01 + 2017-01-01..2017-08-10 | BH? 40 Hz | 557 | 0 | 9.5 GB |
| SYMB | CN | 2017-02-20..2017-04-21 (stage 4) | HH? | 60 | 0 | 1.3 GB |
| B011 | PB | " | EH? | 60 | 0 | 1.4 GB |
| VGZ | CN | " | HH? | 54 | 6 | 1.7 GB |
| LRIV | UW | " | BH? | 60 | 0 | 0.94 GB |
| CLRS | CN | " | HH? | 60 | 0 | 1.2 GB |
| PTRF | CN | " | HH? | 42 | 18 | 0.86 GB |
| SNB | CN | " | HH? | 60 | 0 | 1.9 GB |
| OCP | UW | " | EH? (Z only) | 60 | 0 | 0.82 GB |
| GOBB | CN | " | HH? | **0** | 60 | - |
Total **77 GB**, 3,350 day-files, 93 no-data days, plus 462 MB of template cutouts (below).
**GOBB has no HH data archived for 2017** (its metadata advertises HH?), so stage 4 runs on
**8 added stations**: azimuths from the PGC family centroid 356 (SYMB), 51 (B011), 53 (SNB),
100 (VGZ), 155 (LRIV), 253 (OCP), 281 (PTRF), 327 (CLRS) — good coverage except S-SW.
**DEVIATIONS from the brief, and why:**
1. **ALSE is impossible** for 2017-2018 (UO.ALSE starts 2019-04-19) -> replaced by **UW.YACT**
   (HH? 2013-07..2020-10), which is also D1's **0 %-responsive** station and therefore the hardest
   case / a negative control.
2. **PGC is not run on 2017-2018 for the Lin test.** The Lin catalogue ends 2017-02 and only its
   **ETS-2010** slice is on disk (`data/detections_lin_ets_2010_vi.parquet`, 106,644 event times,
   217 families, 2010-05-17..2010-10-14). The raw Zenodo file (2005-01..2017-02) is reachable by
   HEAD but every GET returned 504 then 403 - not fetched. PGC therefore runs on the **2010 ETS
   window** for ground truth, plus 2015-11..2016-04 (the 2015-12 ETS) and 2017-01..2017-08-10
   (the 2017-03 ETS, for stage 4 and the tremor test). The 2017-08-10 BH->HH sensor change at PGC is
   avoided entirely rather than modelled.
3. **Stage-4 window shortened to 60 days** (2017-02-20..2017-04-21) instead of 2017-2018:
   9 stations x 2 years is ~146 GB and even 9 x 120 d would have put the folder at ~89 GB, over the
   ~80 GB cap. The 60 days span the 2017-02/03/04 northern ETS (PNSN tremor within 60 km of PGC:
   Feb 197, Mar 3,511, Apr 326), which is what the network stacks need. **8 stations kept rather than
   5**: for LOCATION, azimuthal coverage matters more than record length.

## 2026-09-10 D3 stage 2a — 3-COMPONENT templates built

`lfe_catalog_repair/scripts/d3_build_3c_templates.py`. The stored templates are Z-only, so the
horizontals are stacked from the SAME family member times (`<sta>_disc.members.parquet`), with the
same recipe as `discover_gpu.py` (40 Hz, 2-8 Hz, [-1,+1) s, each cut demeaned + unit-normalised, mean
over members, renormalised). Member times lie all over the record and are not in the pilot download,
so +-40 s 3-C cutouts were pulled with `get_waveforms_bulk` and cached
(`lfe_catalog_repair/traces/cutouts/`, 462 MB total - negligible).
**Built-in verification: the REBUILT Z template vs the STORED Z template.**
| station | families with 3-C | member times used | cc(rebuilt Z, stored Z) median | frac > 0.99 |
|---|---|---|---|---|
| WISH | 1,063 / 1,064 | 3,292 (era 2011-06..2019-05) | **1.0000** | 0.793 |
| B943 | 3,407 / 3,407 | 8,765 distinct | **0.9981** | 0.789 |
| YACT | 322 / 322 | 475 (era >= 2013-07) | **0.9808** | 0.394 |
| PGC | (see below) | 10,599 distinct (era < 2017-08-10) | - | - |
The pipeline is therefore reproducing the project's own templates; where cc < 0.99 it is because
some member cutouts were lost (FDSN gap / trace split) so the rebuilt stack has fewer members - YACT
is the worst (median 2 members instead of 3) because its families extend past the era filter.
Engineering note: `obspy.Stream.merge()` on 10k-30k independent cutouts is quadratic and never
returns (killed after 20 min on WISH); the cutouts are contiguous by construction, so the merge was
removed and `cut()` binary-searches the trace list instead.

## 2026-09-10 D3 stage 2b — the CRITICAL calibration result (why 3-C stacking is the whole point)

Ran the re-detector on WISH, Z-ONLY, 2 ETS days, 128 templates, to measure the cc noise level before
committing to a full pass:
- **single-component cc noise: median 0.000, MAD 0.168** (p10-p90 of MAD 0.132-0.191). An 80-sample
  2-8 Hz template has only ~24 independent samples, so sigma_cc ~ 0.2.
- therefore on ONE component **`median + 8*MAD` = 1.37, `+9*MAD` = 1.54, `+10*MAD` = 1.72 -
  ALL ABOVE 1.0, i.e. UNATTAINABLE for 100 % of templates.** A MAD threshold simply cannot be used
  on a single 2 s template channel.
- at the old fixed cc >= 0.80 the same series gives **233 detections/template/day** (2 s gap), which
  reproduces D1's ~108/template/day at a 6 s gap. The machinery matches D1.
**Consequence:** the fixed 0.80 rule is only ~5.9 MAD on a single channel - that is why it admits
~10^2 false alarms per template per day. The MAD rule becomes usable only once the noise sigma is
reduced by stacking: 3 components cut it by ~sqrt(3), so `median + 9*MAD` lands near 0.87 and
`+10*MAD` near 0.97 - i.e. the new rule is STRICTER than the old one, which is the point. This is
recorded here because it reframes D3: the gain is not "a smarter threshold" but "a threshold that
becomes meaningful once the noise floor is lowered by channel stacking".

Stage-2 implementation (`lfe_catalog_repair/scripts/d3_redetect.py`): the `densify_gnw_gpu.py` kernel
with (i) all 3 components and `cc_stack = mean` of the per-component cc, (ii) per template per day
`med`/`MAD` of the stacked series from a stride-80 subsample (43,200 ~independent values/day; raw MAD,
no 1.4826 scaling - the Shelly/EQcorrscan convention), (iii) an extraction floor of
`min(med + 8*MAD, 0.80)` so BOTH the new k = 8/9/10 rules (recomputed per template x 30-day calendar
window from the stored daily stats) AND the old fixed rule are exact subsets of what is stored,
(iv) declustering at the template length (2.0 s) - monotone-safe, so re-thresholding afterwards gives
exactly the peaks a higher threshold would have produced, (v) per detection cc_stack, cc per
component, n_comp, SNR (peak |Z| in the 2 s window / RMS of the preceding 10 s) and the daily med/MAD.
Measured GPU throughput (L40S, tb=256, 5 CPU preprocessors): **0.26 day/s at 1,064 templates x 3
components**, i.e. ~48 min per 2-year station at that template count (~4.8 station-years/hour at
512 template-components).

## 2026-09-10 D3 stage 2b (cont.) — engineering notes from the first full pass

* **GPU OOM at tb=256.** The first WISH pass died at day ~250 with
  `OutOfMemoryError: allocating 2.8 GB (allocated so far: 44.4 GB)` on the 46 GB L40S. At a batch of
  256 templates the live set is 3 x (256 x 3.46e6) float32 cc arrays (the per-component cc must stay
  resident so the per-detection cc can be gathered) + cc_stack + the complex product + the cuFFT
  workspace for N = 3.46e6 x 256. **tb=64 keeps it at ~8.9 GB and is exactly as fast** (measured
  identical, 0.4 min for 12 days x 512 template-components at tb = 64 and 256), so tb=64 is now the
  default; an `OutOfMemoryError` handler also frees the pool, clears the cuFFT plan cache and halves
  the batch for the rest of the run rather than dying.
* **Canonical Z template.** `load_templates` now always uses the project's STORED Z template
  (built from ALL of the family's members) for the Z channel and the 3-C rebuild only for the
  horizontals. This keeps the old-vs-new comparison like-for-like on Z and avoids inheriting the
  member losses of the rebuild (which matter most at PGC: 4,909/5,640 families got 3-C templates,
  rebuilt-Z cc median 0.909, only 9.9 % above 0.99, because the pre-2017-08 era filter removes the
  post-2017 members).
* **Scope trims, stated for the record.** PGC is run on the **2,129 of 5,640 families whose source
  cell is within 5 km of a Lin/Bostock family** (`results/d3_pgc_lin_families.txt`, 287 distinct
  cells) - the only families either PGC purpose (Lin ground truth, stage-4 locations in "the region
  with ground truth") uses - and the 2015-16 bonus era is dropped. Together this cuts ~2.7 h of
  GPU. WISH (1,064), YACT (322) and B943 (3,407) run on ALL of their families.

## 2026-09-10 D3 stage 3 — the "before" leg measured on the SAME window

`lfe_catalog_repair/scripts/d3_oldrule_baseline.py` runs the D1 ETS/quiet + tremor statistics on the
**stored** densify detections (`data/mf_<sta>p90f40_<yr>.csv`, Z-only, fixed cc >= 0.80, 6 s gap,
read-only) restricted to 2017-01-01..2019-01-01, i.e. the same days and cells the new rules are
scored on. Those files hold only the 300 coverage-selected families per station.
| station | stored detections | templates | det/tmpl/day | ETS/quiet median | frac >= 2 | frac >= 5 | Spearman median | frac rho>0 | frac p<0.05 | median cc |
|---|---|---|---|---|---|---|---|---|---|---|
| WISH | 56,028,073 | 300 | **258.3** | 1.133 | 0.000 | **0.000** | +0.038 | 0.727 | 0.213 | 0.827 |
| YACT | 48,617,396 | 300 | **225.1** | n/a (zero quiet-rate) | 0.000 | **0.000** | **-0.127** | 0.007 | 0.806 | 0.819 |
| B943 | 33,001,376 | 394 | **114.7** | 1.210 | 0.109 | **0.018** | +0.074 | 0.914 | 0.515 | 0.820 |
So on the stored catalogue **NO template at WISH or YACT and 1.8 % at B943 reach ETS/quiet >= 5**,
at 100-260 false alarms per template per day, and YACT's daily counts are *anti*-correlated with its
local tremor (rho -0.127, 80.6 % significant, only 0.7 % positive) - a detector tracking noise, not
LFEs. This is the baseline the MAD/3-C rules have to beat, and it reproduces D1 independently
(aggregate ETS/quiet 1.27, 22.6 % of templates "responsive" at ratio >= 1.1).

## 2026-09-10 D3 stage 4a — GEOMETRY CONTROL before any real picking

`lfe_catalog_repair/scripts/d3_stage4_synthtest.py` (-> `results/d3_stage4_synthetic_location_test.txt`)
forward-models S and P at the 9 stage-4 stations for a synthetic source and inverts it with the same
grid search. Forward and inverse share the 1-D model, so this isolates the
**geometry + pick-noise limit** (the formal precision) from the 1-D-model systematic.
| pick sigma | horizontal error | depth error |
|---|---|---|
| 0 (noise-free) | **0.00 km** (exact, all 4 test sources) | **0.0 km** |
| 0.10 s | 0.32-0.49 km | +-0.42-0.91 km |
| 0.20 s | 0.59-0.92 km | +-0.69-1.73 km |
| 0.40 s | 1.40-1.81 km | +-1.49-3.18 km |
(S-only vs S+P: indistinguishable horizontally; P tightens depth by ~30-40 %.)
**So the station geometry is NOT the obstacle to 2-3 km family hypocentres** - picks good to ~0.4 s
suffice. Whether the pilot reaches it therefore depends only on (i) the real stack pick quality and
(ii) the 1-D-model systematic over a subduction zone, which the bootstrap sigma does NOT capture and
which the no-prior depth column is there to expose.
Velocity model (`d3_vel1d.py`): the project has **no** velocity model - the only wave speed anywhere
is the homogeneous coda BETA = 3.5 km/s in `fault_tomography/inversion/build_G_captured.py:13`, a
coda-residence scalar. Stage 4 therefore uses the **Crosson (1976) western-Washington layered P
model** (the PNSN/UW standard for this region) with **Vp/Vs = 1.78**, exact flat-layer shooting for
the direct branch plus a head-wave branch, tabulated at 1 km x 0.5 km and bilinearly interpolated
(cached in `results/d3_tt_table.npz`).

## 2026-09-10 D3 stage 4 pre-check — stacking at the STORED detection times gives NOTHING

Before the re-detected catalogue was ready, the whole stage-4 stacking machinery was run on the
**stored** PGC densify detections (`data/mf_pgc_2017.csv`, fixed Z-only cc >= 0.80: 25,704 rows in
2017-02-20..2017-04-21, 3 families with >= 200) across all 9 stage-4 stations, with the same
random-time null control. Result (`results/d3_stage4_storedtimes_control.{txt,csv}`):
**0 of 27 family-station stacks pass**, and every station's real stack is statistically identical to
its own random-time stack:
| station | snr_h | null snr_h | cc | null cc |
|---|---|---|---|---|
| PGC (0 km) | 1.76 | 1.70 | 0.125 | 0.122 |
| B011 (32 km) | 2.23 | 1.38 | 0.128 | 0.127 |
| SYMB (10 km) | 1.85 | 1.89 | 0.146 | 0.152 |
| VGZ / SNB / CLRS / PTRF / LRIV / OCP | 1.65-2.87 | 1.23-2.58 | 0.124-0.172 | 0.123-0.176 |
The one exception is the **Z** channel at PGC (snr_z 9.07 vs 2.8-4.3 elsewhere) - which is
**circular**: those detections were selected precisely because their Z waveform correlates with the
PGC Z template at cc >= 0.80. Their horizontals carry no coherent S at all.
**This is an independent, waveform-level confirmation of the D1 verdict**: stacking 200 supposed
occurrences of one LFE family, at the station the template came from, produces no S arrival above a
random-time stack, so the stored times are template-shaped noise, not LFE origins. It also validates
the null control (it correctly reports "no signal" where there is none) and the SNR metric
calibration (measured pure-noise levels snr_z ~3.3-3.5, snr_h ~1.7-1.9, cc ~0.09-0.15).

## 2026-09-10 D3 stage 2/3 — WISH DONE, and it is a NEGATIVE result (diagnosed)

**Run.** 727 days, 4 QC-skipped, **2,310,248 stored detections in 49.1 min = 2.43 station-years/hour**
at 1,064 templates x 3 components (L40S, tb=64, 6 CPU preprocessors) -> `results/d3_raw_WISH_*.parquet`.
Measured 3-C stacked-cc noise: **MAD_stack ~= 0.084** (vs 0.168 single-channel), so the k-thresholds
land at **0.669 (k=8) / 0.753 (k=9) / 0.836 (k=10)** and are attainable for ~100 % of
template-windows - exactly the mechanism stage 2b predicted (sigma falls by ~sqrt(3): 0.168/sqrt(3)
= 0.097 predicted, 0.084 measured).

**Result table** (`results/d3_station_summary.csv`, `d3_validation_WISH.csv`; the stored Z-only rule
row is from `d3_storedrule_baseline.csv`, same window, its 300 families):
| rule | detections | det/tmpl/day | templates firing | ETS/quiet median | frac >= 2 | **frac >= 5** | Spearman median | frac rho>0 | median SNR |
|---|---|---|---|---|---|---|---|---|---|
| STORED Z-only cc>=0.80 (300 fams) | 56,028,073 | 258.3 | 300 | 1.13 | 0.000 | **0.000** | +0.038 | 0.727 | - |
| 3-C stack cc>=0.80 ("old rule, same series") | 1,983 | 0.0026 | 794 | 0.00 | 0.003 | **0.000** | +0.008 | 0.425 | 2.04 |
| 3-C MAD k=8 (thr 0.669) | 297,202 | 0.387 | 1,060 | 0.00 | 0.015 | **0.000** | +0.00002 | 0.498 | 2.71 |
| 3-C MAD k=9 (thr 0.753) | 54,601 | 0.071 | 948 | 0.00 | 0.010 | **0.0009** | +0.0026 | 0.476 | 2.85 |
| 3-C MAD k=10 (thr 0.836) | 6,529 | 0.0085 | 573 | 0.00 | 0.004 | **0.0009** | +0.0017 | 0.280 | 2.23 |
The false-alarm load collapses by 3-5 orders of magnitude (258 -> 0.07 det/template/day at k=9), which
is the intended effect. **But no rule produces an ETS response or a tremor correlation at WISH.**

### Why - three candidate explanations, tested
1. **"The adaptive threshold cancels the ETS signal" (tremor raises the cc noise floor).**
   **REJECTED.** `results/d3_diag_fixedthr_WISH.csv`: Spearman(daily MAD, local 30-d tremor) median
   **+0.020** (61 % positive), and mean MAD on ETS days / quiet days = **1.008** [p10 0.987,
   p90 1.025]. The noise floor is essentially not ETS-modulated, so the threshold is not eating the
   signal.
2. **"The threshold is simply in the wrong place."** **REJECTED.** Sweeping a FIXED cc_stack
   threshold with no adaptivity at all:
   | fixed cc_stack | detections | ETS/quiet median | frac >= 2 | frac >= 5 | Spearman median | frac rho>0 |
   |---|---|---|---|---|---|---|
   | 0.60 | 333,208 | 0.00 | 0.008 | 0.000 | +0.0041 | 0.548 |
   | 0.65 | 132,469 | 0.00 | 0.020 | 0.000 | +0.0037 | 0.545 |
   | 0.70 | 37,220 | 0.00 | 0.029 | 0.000 | +0.0049 | 0.566 |
   | 0.75 | 10,254 | 0.00 | 0.013 | 0.000 | +0.0047 | 0.555 |
   | 0.80 | 1,983 | 0.00 | 0.002 | 0.000 | +0.0072 | 0.423 |
   There is no threshold anywhere between 0.60 and 0.80 at which a tremor-linked population appears.
3. **"What IS being detected is impulsive noise."** **SUPPORTED.** Per-detection structure:
   * strongly **Z-dominated**: at k=9 the median cc_Z is 0.832 but cc_N 0.625 and cc_E 0.651
     (median |cc_Z - cc_N| = 0.217) - the horizontals contribute little;
   * **hour-of-day is grossly non-uniform** (chi2 = 10,741 at k=9, 67,051 at k=8, against 35.2 for
     p = 0.05 on 23 df) -> cultural/diurnal;
   * **strongly day-clustered**: 3 of 214 active days hold 50 % of the cc_stack>=0.80 detections;
   * **multi-template simultaneity**: 28 % (k=9) and 69 % (k=8) of detections fall within 1 s of a
     detection on *another* template - a broadband transient lights up every template at once,
     which an LFE does not;
   * **low SNR**: median 2.0-3.0 (peak |Z| in the 2 s window / RMS of the preceding 10 s).
So at WISH the single-station 3-C matched filter, at every threshold, is detecting impulsive site
noise. This is consistent with the stage-4 pre-check (stacking the stored times at 9 stations gives
nothing above a random-time null) and with D1.

**What this does NOT yet establish** is whether the method or the target is at fault: WISH may simply
have no detectable LFEs. The decisive test is PGC over the Lin/Bostock ground-truth window, which is
therefore promoted to the FRONT of the GPU queue (`scripts/d3_queue3.sh`).

## 2026-09-10 D3 amendment 2 received (user, via orchestrator): add STAGE 5 = NETWORK matched filter
Root cause named by the user, and it matches the stage-2b arithmetic: a 2 s single-station template
has only ~24 independent samples (sigma_cc ~0.17) and 3-C stacking buys just sqrt(3). The standard
remedy is the NETWORK matched filter - sum the per-channel cc over all stations carrying a
network-consistent template for the family, aligned by the family moveout, and threshold the network
sum (sigma ~ 0.17/sqrt(n_channels) ~ 0.04 for 8 stations x 3 comps). Stage 5 = that, in the 60-day
southern-Vancouver-Island window, plus a template-LENGTH test (2 / 4 / 6 s) on a 10-day ETS subset,
plus a fleet-design recommendation with GPU cost. Spec in the orchestrator message.

## 2026-09-10 D3 stage 3 — PGC + the GROUND-TRUTH AUDIT: the sensitivity budget, measured

### PGC 2010 re-detection (the Lin/Bostock window)
184 days, 5,416,787 stored detections in 22.3 min (**1.36 station-years/hour** at 2,129 templates x
3 components; 1,855 of the 2,129 have all three). Data note: **12.9 % of template-days are
component-short and they carry 97.7 % of the raw detections** (on a Z-only day the stacked cc
collapses onto the single channel, MAD doubles and the floor falls back to 0.80) - the component
guard removes them from numerator and denominator alike.
| rule | detections | det/tmpl/day | templates firing | ETS/quiet median | frac >= 5 | Spearman median | median SNR | median cc | thr |
|---|---|---|---|---|---|---|---|---|---|
| 3-C stack cc>=0.80 | 12,052 | 0.035 | 234 | **1.66** | 0.0009 | +0.090 | 3.34 | 0.847 | - |
| 3-C MAD k=8 | 87,172 | 0.255 | 683 | 0.33 | 0.0038 | +0.022 | 2.93 | 0.729 | 0.700 |
| 3-C MAD k=9 | 21,615 | 0.063 | 250 | **2.06** | 0.0038 | +0.089 | 3.22 | 0.832 | 0.788 |
| 3-C MAD k=10 | 5,422 | 0.016 | 124 | **5.12** | 0.0047 | +0.090 | 3.62 | 0.945 | 0.875 |
Unlike WISH, the ETS/quiet ratio at PGC **rises monotonically with the threshold** (1.66 -> 2.06 ->
5.12) and the SNR rises with it (3.34 -> 3.62): a small real population is being isolated. But only
124 of 2,129 templates fire at k=10 and only 0.47 % clear ETS/quiet >= 5.

### Lin/Bostock ground truth — recall/precision = the shifted control
Per template the rates are far too low for a per-template score, so the test is POOLED (one global
lag fitted on the even-index Lin events of every Lin-proximate template, scored on the odd-index
ones; 970 of 1,013 PGC templates lie within 5 km of a Lin family, 117 distinct Lin families):
| rule | detections | Lin events | fitted lag | recall | precision | recall (+600 s control) | precision (control) |
|---|---|---|---|---|---|---|---|
| cc>=0.80 | 11,492 | 49,908 | +5.38 s | 0.0006 | 0.0014 | 0.0010 | 0.0018 |
| MAD k=8 | 69,024 | 147,524 | +7.38 s | 0.0043 | 0.0029 | 0.0046 | 0.0041 |
| MAD k=9 | 17,169 | 59,133 | +3.12 s | 0.0009 | 0.0029 | 0.0017 | 0.0038 |
| MAD k=10 | 4,447 | 28,537 | +3.12 s | 0.0003 | 0.0020 | 0.0005 | 0.0031 |
**Every rule matches the ground truth at exactly the chance rate** (the shifted control is equal or
higher), and the fitted lags do not converge on the project's calibrated +10 s.
(A real bug was found and fixed on the way: the Lin parquet stores `datetime64[us]`, so
`.astype("int64")/1e9` yielded epoch ~1.27e6 and an EMPTY intersection with the detections. The same
trap bites `lfe_features/data/lin_precise_b011.csv`. Both are now pinned to ns.)

### WHY - the sensitivity budget, measured from the ground truth itself
First a decisive correction. **The Lin ORIGIN times are 11.28 s before the S arrival with a scatter
of +-2.2 s** (measured from the project's own refined picks: median arrival-OT 11.28 s, p05 8.0,
p95 14.0, std 2.22 s - which is exactly why `refine_lin_picks.py` exists). At 1-8 Hz a +-2.2 s label
error destroys stack coherence outright, so my first visibility test - stacking at OT + predicted S
and searching +-4 s - was **invalid**, and so was its wide-search repeat (25 Lin families x 800
events, snr_h 2.44 vs null 2.28, peak position scattered +3.5..+25.5 s: all label jitter).
Redone at the **refined** arrival times (`scripts/d3_lin_refined_stack.py`):
| stack | events | snr_h | null | ratio | peak vs refined arrival |
|---|---|---|---|---|---|
| all curated picks in the window | 776 | **4.27** | 2.27 | **1.88x** | **-0.62 s** |
| one 0.1-deg source cell | 130 | **2.98** | 1.61 | **1.86x** | **-0.45 s** |
**The LFE arrival IS there at PGC, lands where it should, and the D3 processing chain is sound.**
The cc-refined variant collapses to the null on the 776-event set (0.99x) - honest evidence that
alignment refinement at this SNR just aligns noise.
Now the number that decides everything. From `rms(signal window)/rms(null window) = 1 + N*(s/n)^2`
under the square root, with the single-cell 130-event stack at L = 6 s (ratio 1.34):
**a single median Lin LFE at PGC has amplitude s/n = 0.078 of the ambient 1-8 Hz noise** (0.032 from
the 776-event mixed-source stack, which is smeared by source mixing and is a lower bound).
Measured template-window energy ratios rise with length (1.22 at 2 s, 1.27 at 4 s, **1.34 at 6 s**),
so the coherent LFE energy at PGC lasts **longer than 6 s** - a 2 s template throws most of it away.
Matched-filter arithmetic (single-channel cc sigma = 1/sqrt(2*B*L), measured 0.168 at L = 2 s vs
0.189 predicted - agreement confirms the model; detection needs s/n >= k * 0.6745 * sigma/sqrt(Nch)):
| template L | sigma_1ch | channels needed for k=9 | for k=8 | k reached by 1 station (3-C) | by 9 stations (27 ch) |
|---|---|---|---|---|---|
| 2 s | 0.189 | **216 ch (72 stations)** | 171 | 1.06 MAD | 3.2 MAD |
| 4 s | 0.134 | 109 ch (36 stations) | 86 | 1.50 MAD | 4.5 MAD |
| 6 s | 0.109 | **72 ch (24 stations)** | 57 | **1.84 MAD** | **5.5 MAD** |
**So the D3 negative result is a SENSITIVITY result, not a method error**: with a 2 s template on one
3-component station the median LFE sits at ~1 MAD of the cc noise - it is not detectable at ANY
threshold, which is precisely what the WISH fixed-threshold sweep (0.60-0.80, no ETS response
anywhere) and the PGC ground-truth test show. It also explains D1: a detector running at ~1 MAD of
signal admits only noise, so its rate tracks the noise level and rises ~25 % with tremor at every cc.
**Prediction for Stage 5, recorded before running it:** the 9-station / 27-channel network with 6 s
templates should reach ~5.5 MAD for a median LFE - i.e. it should detect the upper part of the
amplitude distribution at k ~ 5-6 but NOT at k = 9. Stage 5 therefore scans k = 4,5,6,8,9,10.

## 2026-09-10 D3 stage 4 family locations — the geometry is fine, the template TIMES are not

**Geometry control first (`scripts/d3_stage4_synthtest.py`, run before any real picking).** With the
9 stage-4 stations and the 1-D model used for both forward and inverse, the grid search recovers a
synthetic source EXACTLY (0.00 km, 0.0 km at four test points) and degrades only with pick noise:
err_h 0.32-0.49 km at sigma 0.10 s, 0.59-0.92 km at 0.20 s, 1.40-1.81 km at 0.40 s; err_z +-0.4-0.9,
+-0.7-1.7, +-1.5-3.2 km respectively (P tightens depth by 30-40 %, and barely affects the
horizontal). **So 2-3 km family hypocentres are well within this network's reach - the limit is not
the geometry.**

**But the 2017 window cannot supply template times.** Stage 4 needs >= 200 re-detections per family
inside the 60-day stage-4 window; over 2,129 PGC templates that window contains only
**3,899 detections in total** (old rule 72, with a maximum of ONE per family; k=9 318; k=10 52). Only
the loosest rule, k=8, has any family above 200 (3 families, max 458). Run on those:
| station | dist | snr_z (null) | snr_h (null) | cc (null) | kept |
|---|---|---|---|---|---|
| **PGC** (0 km) | - | **6.56 (3.38)** | **5.04 (1.50)** | 0.141 (0.112) | 2/3 |
| B011 (32 km) | | 4.47 (3.49) | 2.07 (1.76) | 0.122 (0.114) | 0/3 |
| SYMB (10 km) | | 3.84 (2.87) | 1.61 (1.48) | 0.147 (0.147) | 0/3 |
| VGZ / SNB / CLRS / PTRF / LRIV / OCP | 35-64 km | 3.3-3.8 (2.8-3.8) | 1.6-3.0 (1.7-3.0) | 0.11-0.18 | 0/3 |
Only PGC's own stack rises above its null - and that is **circular**, since those detections were
made with PGC's own template. **0 of 3 families reach the 4 stations a location needs, so 0 families
are located** (`results/d3_stage4_families.csv` is empty by outcome, not by error;
`figures/d3_stage4_locations.png` shows the geometry and the station set).
Method note: the mean cut-to-stack **cc had to be dropped as a gate** - the demonstrably real
Lin-refined stack has cc 0.105 against a null of 0.099, so at this SNR only the stack-SNR ratio
against a same-n random-time null discriminates. cc is retained as a diagnostic.
Velocity model, Slab2 handling and the picker are as described in the stage-4a entry; the Lin family
DEPTHS in `catalogs/lin_families_ets_2010_vi.csv` (12-35 km where the interface is 30-45 km) are
**not trustworthy** - `lin_catalog.py` sign-flips Lin's depth and averages over a 0.1 deg cell - so
stage 5 uses a fixed 35 km source depth for the moveout instead.

## 2026-09-10 D3 stage 5 — NETWORK matched filter: it works as theory says, and it is still ~3x short

Templates could not come from the PGC re-detections (stage 4 showed they stack to nothing), so they
are built from the project's **refined Lin picks**, each event corrected by ITS OWN predicted S
moveout to a common PGC reference (`--single-cell`), fit half (388 picks) for the templates and the
other 388 held out for scoring. 9 stations / 21 channels, 2010-08-10..2010-09-13.

**Per-station template stacks (388 events, RAW - no alignment refinement):**
| station | dist | snr_h | null | ratio | S pick vs predicted |
|---|---|---|---|---|---|
| **PGC** | 0 km | **3.83** | 1.64 | **2.34x** | +1.38 s |
| **OCP** (Z only) | 64 km | **4.50** | 2.39 | **1.88x** | -0.25 s |
| VGZ | 35 km | 2.08 | 1.75 | 1.19x | +1.98 s |
| SNB | 57 km | 1.48 | 1.25 | 1.18x | +1.12 s |
| B011 | 32 km | 1.84 | 1.65 | 1.12x | -0.03 s |
| B926 | 47 km | 1.93 | 1.76 | 1.10x | +3.25 s |
| STW | 36 km | 3.20 | 3.05 | 1.05x | +2.17 s |
| GOBB / LRIV | 58 / 50 km | 3.39 / 1.52 | 4.41 / 1.95 | 0.77x / 0.78x | +0.25 / -1.20 s |
Only PGC and OCP clear a 1.3x-null gate; the others are weakly positive (1.05-1.19x) - real but
sub-threshold, so the network sum was run over **all** of them (`--use-all`), which is exactly what
a network sum is for.

**The network sum behaves exactly as the theory predicts:**
| template L | sigma_1ch predicted | sigma_net predicted (21 ch) | **sigma_net MEASURED** | improvement over 1 channel |
|---|---|---|---|---|
| 2 s | 0.189 | 0.0412 | **0.0537** | 3.5x |
| 4 s | 0.134 | 0.0292 | **0.0406** | 3.3x |
| 6 s | 0.109 | 0.0238 | **0.0329** | 3.3x |
So summing 21 channels really does divide the cc noise by **3.3x** (not the ideal sqrt(21) = 4.6 -
the shortfall is correlated noise between components/nearby stations plus imperfect moveout), and
**longer templates help as ~1/sqrt(L), confirming that the fleet's 2 s template throws away most of
the available gain.**

**But it is still not enough.** Recall/precision against the HELD-OUT refined picks (+-1 s), with a
+600 s shifted control:
| L | k | detections | recall | precision | precision (shifted) | frac >= 3 stations |
|---|---|---|---|---|---|---|
| 6 s | 4 | 132,248 | 0.0979 | 0.0003 | 0.0003 | 0.012 |
| 6 s | 5 | 27,085 | 0.0077 | 0.0001 | 0.0002 | 0.033 |
| 6 s | 6 | 2,689 | 0.0026 | 0.0004 | 0.0000 | 0.093 |
| 6 s | 8 | 6 | 0 | 0 | 0 | 0.500 |
| 6 s | 9/10 | **0** | - | - | - | - |
(2 s and 4 s are the same picture.) **Precision equals the shifted control at every k and every
length**: the 132,248 detections at k=4 are the noise floor itself (~3,900/day), and by k=8 the
detector is silent. Nothing real is being found.
This is precisely what the measured numbers require: with sigma_net = 0.0329 and a signal of
s/n = 0.078, a median LFE reaches **k = 0.078/(0.6745 x 0.0329) = 3.5 MAD** - it sits *inside* the
k=4 noise cloud and disappears above it. My pre-registered prediction was 5.5 MAD; it came out 3.5
because only 21 channels were available (not 27), the measured sigma_net is 1.38x worse than ideal,
and 7 of the 9 stations carry only 1.0-1.2x-null templates.
**Scaling to a threshold that works**, using the MEASURED (not ideal) channel scaling:
* k = 6 needs sigma_net <= 0.0193 -> **~61 channels ~ 20 three-component stations**
* k = 9 needs sigma_net <= 0.0128 -> **~140 channels ~ 47 three-component stations**
both with **6 s** templates (2 s templates need ~2.7x more channels again).

## 2026-09-10 D3 stage 3 — the four single-station legs, cleaned, and two QC lessons

Two data-quality guards were added to `d3_validate.py` after the first pass and every station was
re-scored with them:
1. **component guard** - on a day missing a horizontal the stacked cc collapses onto one channel,
   its MAD doubles and the extraction floor falls back to 0.80, so such days dump a burst of
   spurious "new" detections. Affected 0.09 % of WISH template-days but **12.9 % at PGC2010,
   12.8 % at PGC2017 and 3.1 % at YACT - and those days carried 97.7 %, 99.7 % and 98.5 % of the
   raw detections.**
2. **degenerate-window guard** - a zero-filled / constant segment yields detections whose 2 s window
   has literally zero amplitude. **57.8 % of WISH's and 26.7 % of YACT's stored detections** were
   such windows (at YACT a single day held 50 % of the k=9 detections and dragged the median SNR to
   1e-8). Both are cheap gates the fleet should adopt.

Final per-station scoreboard (2017-2018 for WISH/YACT/PGC2017, 2010-05..11 for PGC2010):
| station | rule | detections | det/tmpl/day | ETS/quiet median | frac >= 5 | Spearman median | median SNR |
|---|---|---|---|---|---|---|---|
| WISH | stored Z-only 0.80 | 56,028,073 | 258.3 | 1.13 | 0.000 | +0.038 | - |
| WISH | 3-C 0.80 / k=8 / k=9 / k=10 | 361 / 248,863 / 44,951 / 3,961 | 0.0005 / 0.32 / 0.059 / 0.005 | 0.63 / 0.49 / 0.37 / 0.74 | 0 / 0.002 / 0.001 / 0.001 | +0.034 / **-0.018** / -0.009 / -0.004 | 2.3-2.8 |
| YACT | stored Z-only 0.80 | 48,617,396 | 225.1 | n/a | 0.000 | **-0.127** | - |
| YACT | 3-C 0.80 / k=8 / k=9 / k=10 | 108 / 560 / 39 / 5 | 0.0005 / 0.003 / 0.0002 / 0.00002 | n/a (no quiet-day detections) | 0 | +0.06 / +0.03 / +0.06 / +0.06 | 2.9-3.4 |
| PGC2010 | 3-C 0.80 / k=8 / k=9 / k=10 | 12,052 / 87,171 / 21,615 / 5,422 | 0.035 / 0.26 / 0.063 / 0.016 | **1.66 / 0.33 / 2.06 / 5.12** | 0.001-0.005 | +0.090 / +0.022 / +0.089 / +0.090 | 3.3-3.6 |
| PGC2017 | 3-C 0.80 / k=8 / k=9 / k=10 | 82 / 12,631 / 1,454 / 114 | 0.0002-0.03 | 0 / 0 / 0 / n/a | 0.000 | +0.09 / +0.04 / +0.09 / +0.09 | 2.9-3.5 |
**The false-alarm load is crushed** - 258 -> 0.06 det/template/day at WISH (a 4,000x reduction), and
the detections that survive have higher SNR. **But the ETS response and the tremor correlation are
not recovered anywhere except weakly at PGC2010**, where ETS/quiet rises monotonically with the
threshold (1.66 -> 2.06 -> 5.12) and rho stays ~+0.09. At the pass bar the answer is unambiguous:
**frac(ETS/quiet >= 5) is 0.000-0.005 at every station under every rule** (the D1 stored rule scored
0.000, 0.000 and 0.018).

Measured GPU throughput (L40S, tb=64, 6 CPU preprocessors), which is what the fleet cost rests on:
| station | templates x comps | days | wall | station-years/hour |
|---|---|---|---|---|
| YACT | 966 | 724 | 15.1 min | **7.88** |
| WISH | 3,192 | 727 | 49.1 min | **2.43** |
| PGC 2010 | 6,387 | 184 | 22.3 min | **1.36** |
| PGC 2017 | 6,387 | 221 | 26.6 min | **1.37** |
The product is constant: **~8,000 (template-components x station-years) per hour**, and it is
independent of template length (the cost is the full-day FFT). Useful conversions: a fleet station
with the standard 300 densify families x 3 components runs at **~8.9 station-years/hour**, i.e.
~1.9 h for a 17-year record, ~15 GPU-days for the whole 187-station fleet. A 20-station NETWORK
matched filter over 300 families x 3 components (18,000 channels) runs at **~0.44 station-years per
hour**, i.e. ~39 GPU-hours per 17-year source region.

## 2026-09-10 D3 — FLEET DESIGN RECOMMENDATION (from the measured budget, not from theory alone)

Everything below follows from three measured numbers: the per-event LFE amplitude at PGC
**s/n = 0.078**, the single-channel cc noise **sigma = 1/sqrt(2BL)** (measured 0.168 at L = 2 s
against 0.189 predicted), and the **measured** network scaling **3.3x for 21 channels** (ideal 4.6x;
PGC-only 3 channels gave sigma = 0.0730, the 21-channel sum 0.0329, a direct 2.22x measurement
against an ideal of 2.65x).

1. **TEMPLATE LENGTH: 6 s, not 2 s.** The measured network sigma falls 0.0537 -> 0.0406 -> 0.0329
   from L = 2 -> 4 -> 6 s (~1/sqrt(L), as theory says), and the coherent LFE energy in the
   ground-truth stack is still growing at 6 s (window energy ratio 1.22 / 1.27 / 1.34). A 2 s
   template needs ~2.7x more channels than a 6 s one for the same threshold. Going beyond 6 s is
   worth testing - this pilot did not bound the upper end.
2. **CHANNELS: all three components, always.** 3-C stacking is what makes a MAD threshold usable at
   all (on one channel `median + 8*MAD` = 1.37 > 1 and is unattainable by construction).
3. **NETWORK, not single-station - and the network has to be BIG.** With the measured scaling:
   * **k = 6 needs ~61 channels ~ 20 three-component stations**
   * **k = 9 needs ~140 channels ~ 47 three-component stations**
   Nine stations / 21 channels reach only k = 3.5 for a median LFE, and that is inside the noise.
   This is why Lin/Bostock-class catalogues are built on tens of stations, and it is the single
   most important design conclusion of D3.
4. **Template TIMES must come from a verified source.** Every attempt to bootstrap templates from
   this project's own detections failed: stacking at the stored densify times gives nothing at 9
   stations (stage-4 pre-check), and stacking at the D3 re-detections gives nothing except the
   circular PGC self-match. Templates must be seeded from an external LFE catalogue (Lin/Bostock,
   Ducellier) or from a network detection at a threshold that is actually reachable.
5. **Two cheap QC gates the current pipeline lacks**: drop template-days below the modal component
   count, and drop zero-amplitude (degenerate-window) detections. They removed 58 % of WISH's and
   27 % of YACT's raw detections here.
6. **Cost.** Single-station 3-C re-detection: ~8,000 template-components x station-years per hour
   (~1.9 GPU-hours per 17-year station at 300 families, ~15 GPU-days for the 187-station fleet),
   independent of template length. A 20-station network matched filter over 300 families:
   ~39 GPU-hours per 17-year source region. Both are affordable; the binding constraint is the
   number of stations that must be correlated **together**, not the GPU.

## 2026-09-10 D3 pilot results — VERDICT against the PLAN's pass criteria

| criterion (PLAN "After D1: D3 design") | result | verdict |
|---|---|---|
| **Lin/Bostock recall & precision** | pooled over 970 Lin-proximate PGC templates: recall 0.0003-0.0043, precision 0.0014-0.0029 - **equal to or below the +600 s shifted control** at every rule. The network matched filter (21 channels, L = 6 s) is the same: precision 0.0001-0.0004 vs a control of 0.0000-0.0003 at every k. | **FAIL** |
| **ETS/quiet >= 5 for responsive templates** | frac(ETS/quiet >= 5) = **0.000-0.005** at all four stations under all four rules (stored Z-only rule: 0.000 WISH, 0.000 YACT, 0.018 B943). Only PGC2010 shows a monotone rise with threshold (1.66 -> 2.06 -> 5.12 median). | **FAIL** (one weak positive) |
| **tremor correlation** | Spearman(daily count, local 30-d tremor) median -0.018..+0.09; frac positive 0.02-0.19. A fixed-threshold sweep at WISH (cc_stack 0.60-0.80) gives rho ~ +0.004 and no ETS response anywhere. | **FAIL** |
| **borehole/surface agreement** | see the WISH x B943 coincidence entry below | - |
| **GPU throughput (for the fleet cost)** | **~8,000 template-components x station-years per hour**; 2.43 sy/h at 1,064x3, 7.88 at 322x3, 1.36-1.37 at 2,129x3. | **PASS** |

### What D3 actually established
1. **The re-detection machinery is correct.** It reproduces D1 independently (233 det/template/day
   at a fixed Z-only 0.80, vs D1's ~108 at a 6 s gap), it reproduces the project's own templates
   (rebuilt-Z vs stored-Z cc median 1.0000 at WISH, 0.998 at B943), its cc-noise measurement matches
   theory (0.168 measured vs 0.189 predicted at L = 2 s), its network sum divides that noise by a
   measured 3.3x at 21 channels, and its nulls behave (pure-noise stacks land exactly where the
   controls say they should).
2. **The MAD threshold does what it was supposed to do**: 3-C stacking drops the cc noise MAD from
   0.168 to 0.084, which is what makes `median + k*MAD` attainable at all, and the false-alarm load
   falls from 258 to 0.06 detections/template/day (4,000x) with higher-SNR survivors.
3. **It still finds no LFEs, and the reason is measured, not guessed.** A single median Lin/Bostock
   LFE at PGC has amplitude **s/n = 0.078** of the ambient 1-8 Hz noise. Against a single-channel cc
   sigma of 0.109 (L = 6 s) that is **1.8 MAD** on one 3-component station and **3.5 MAD** on the
   9-station / 21-channel network - inside the noise in both cases. **k = 6 needs ~20
   three-component stations and k = 9 needs ~47.** The D3 negative result is a sensitivity budget,
   not a bug.
4. **Three things in the existing pipeline are now known to be wrong or fragile**, each with a
   number: the **2 s template** (needs 2.7x more channels than a 6 s one), the **single-station**
   design (1.8 MAD), and the **absence of component/degenerate-window QC** (58 % of WISH's raw
   detections were zero-amplitude windows).
5. **Family LOCATION is not blocked by geometry.** The synthetic control recovers a source exactly
   with noise-free picks and to 0.6-0.9 km horizontal / +-0.7-1.7 km depth at 0.2 s pick noise, so
   **2-3 km family hypocentres are reachable with this station geometry** - what is missing is a
   template time base good enough to make network stacks pick-able, which is the same sensitivity
   problem.

### Honest caveats
* The s/n = 0.078 figure is one number, from one station, for the MEDIAN of 130 curated Lin events in
  one 0.1-deg cell; the amplitude distribution is broad and the largest LFEs are certainly easier.
  It should be re-measured per station before the fleet design is frozen.
* The Lin ground truth used here is 776 curated, cross-correlation-refined picks in a 34-day window.
  The full Zenodo catalogue (1.05 M detections, 2005-2017) was NOT obtained (Zenodo returned 504
  then 403); with it, the recall test could be run per family and over many more years.
* PGC's horizontals are missing on 12.9 % of days (carrying 97.7 % of the raw detections before the
  guard), and 274 of its 2,129 templates have no horizontal at all - PGC is not the cleanest possible
  test station.
* Stage 4/5 templates rest on a 35 km assumed source depth because the Lin family depths on disk are
  not trustworthy (12-35 km where the interface is 30-45 km, after `lin_catalog.py` sign-flips and
  cell-averages them).

## 2026-09-10 D3 — B943 borehole leg and the WISH x B943 network coincidence

**B943 re-detection** (1,927 of 3,407 families - those in the 160 source cells it shares with WISH):
730 days, 2,530,698 stored detections in 78.8 min (**1.52 station-years/hour**); 5.5 % of them were
zero-amplitude windows. No component-short days at all (the PB borehole EH1/EH2/EHZ record is clean).
| rule | detections | det/tmpl/day | ETS/quiet median | frac >= 2 | frac >= 5 | Spearman median |
|---|---|---|---|---|---|---|
| stored Z-only cc>=0.80 (394 fams) | 33,001,376 | 114.7 | 1.21 | **0.109** | **0.018** | +0.074 |
| 3-C stack cc>=0.80 | 124 | 0.0001 | 0.00 | 0.000 | 0.000 | +0.049 |
| 3-C MAD k=8 | 358,820 | 0.255 | 0.449 | **0.055** | 0.0042 | -0.024 |
| 3-C MAD k=9 | 20,581 | 0.015 | 0.191 | 0.020 | 0.0016 | -0.010 |
| 3-C MAD k=10 | 686 | 0.0005 | 0.00 | 0.0005 | 0.0005 | +0.008 |
Honest note: at the borehole the **STORED rule scores better** on ETS responsiveness (frac >= 2 of
0.109 and frac >= 5 of 0.018) than any new rule (0.055 / 0.0042 at best). The new rules cut the
false-alarm load 450x but do not improve, and on this metric degrade, the ETS contrast - the same
sensitivity ceiling, seen from the other side.

**Network coincidence, surface WISH x borehole B943** (159 shared source cells, 104 km apart;
same construction as D1: lag from the pair-difference histogram during ETS, +-3 s coincidence,
+600 s shifted control):
| rule | cells with enough detections | median lag | coincident fraction, ETS | quiet | +600 s control | ETS/control | ETS/quiet |
|---|---|---|---|---|---|---|---|
| k=8 | 29 / 159 | +0.0 s | 0.0011 | 0.0014 | 0.0008 | **1.3x** | 0.81 |
| k=9 | 1 / 159 | +0.0 s | 0.0000 | 0.0000 | 0.0000 | - | - |
**The re-detected catalogues show no network coincidence at all** (1.3x the shifted control, and
fewer coincidences during ETS than during quiet weeks). Compare D1, where the STORED hits gave
**11.4x** above the shifted control - a population D1 had already identified as shared SURFACE noise
confined to the WISH/RADR/LEBA triple and invisible to this same borehole. Removing the false alarms
removed the coincidences with them, and nothing real took their place. Borehole/surface agreement:
**none measurable.**

## 2026-09-10 D3 — files, and where to pick this up

**Scripts** (all in `lfe_catalog_repair/scripts/`, all reading the project's data read-only):
`d3_download.py` (3-C download into the study folder) - `d3_build_3c_templates.py` (horizontal
templates from the stored family member times, with the rebuilt-vs-stored Z check) -
`d3_redetect.py` (the GPU re-detector: 3-C cc stack, per-day/per-window MAD thresholds, OOM-retry) -
`d3_validate.py` (ETS/quiet + tremor stats, component and degenerate-window guards, truncation
diagnostic) - `d3_oldrule_baseline.py` (the stored Z-only rule on the same window) -
`d3_network.py` (WISH x B943 coincidence) - `d3_lin.py` (Lin ground truth, per-template + pooled) -
`d3_lin_template_test.py` and `d3_lin_refined_stack.py` (the visibility audit and the s/n
measurement) - `d3_vel1d.py` (Crosson 1976 + Slab2 + travel times) - `d3_stage4_{stacks,locate,fig,
synthtest}.py` (network stacks, picks, grid-search locations, geometry control) -
`d3_stage5_network.py` + `d3_stage5_fig.py` (the network matched filter and its figure) -
`d3_figures.py` (the pilot validation panel) - `d3_stored_to_parquet.py` - drivers `d3_stage3.sh`,
`d3_stage4.sh`, `d3_queue*.sh`, `d3_run_redetect.sh`, `d3_stage*_dl*.sh`.

**Key results** (`lfe_catalog_repair/results/`): `d3_station_summary.csv` (the scoreboard),
`d3_detections_<STA>.parquet` + `d3_thresholds_<STA>.csv` + `d3_validation_<STA>.csv` for WISH /
YACT / B943 / PGC2010 / PGC2017, `d3_raw_<STA>_{detections,daily}.parquet` (pre-QC),
`d3_storedrule_baseline.csv`, `d3_diag_fixedthr_WISH.csv` (the fixed-threshold sweep),
`d3_lin_groundtruth{,_pooled}.csv`, `d3_lin_refined_stack{,_cell}.txt` (**the s/n = 0.078
measurement**), `d3_stage4_{stacks,families,picks,events}.csv`,
`d3_stage4_synthetic_location_test.txt` (the geometry control),
`d3_stage4_storedtimes_control.txt`, `d3_stage5_{stacks,summary}.csv` +
`d3_stage5_detections_L{2,4,6}.parquet` + `pgconly/d3_stage5_summary.csv`,
`d3_coincidence_WISH_B943.csv`, `d3_templates3c_<STA>.npz` (+ `_report.csv`).
**Figures**: `d3_pilot_validation.png` (12 panels, old vs new), `d3_stage5_network.png`
(sigma vs channels and template length - the design figure), `d3_stage4_locations.png`.
**Waveforms**: `lfe_catalog_repair/traces/` = 79 GB, 13 stations (the 80 GB cap is reached; nothing
more can be downloaded without deleting).

**If this is picked up again, the three things worth doing, in order:**
1. **Get the full Lin/Bostock catalogue** (Zenodo 10016020 `EQloc_001_0.1_3_S.csv`, 1.05 M
   detections 2005-2017; every GET returned 504 then 403 today) and **re-measure s/n per station**.
   The whole design rests on one number from one station.
2. **Build the ~20-station network** for one source region and run `d3_stage5_network.py` with
   L = 6 s at k = 6. That is the smallest configuration the measured budget says can work, and the
   script already does it - it only needs the stations.
3. **Do not spend more GPU on single-station re-detection.** It is measured at 1.8 MAD; no
   threshold rule can rescue it. The stored `data/mf_*.csv` remain what D1 said they are: stacking
   material, not an event list.

## 2026-09-10 orchestrator review — D3 pilot complete: negative with a measured cause (see RESULTS.md)
Individual LFEs (median Lin event s/n 0.078) are 1.8 MAD on a 3-C station and 3.5 MAD on 21 channels -> not detectable by template
matching at this scale; ~20 3-C stations with 6-s templates needed at k=6. Options A/B/C written to RESULTS.md; recommendation: B
(deep-learning picker + association + grid-search location, continuing the project's lfe_features thread), C for validation.

## 2026-09-10 orchestrator — reconciliation with the earlier picker work (lfe_features/ R1-R11, June 2026; user: "you didn't look")
What that thread already established: (R4) a per-window 3-class LFE/EQ/NOISE RandomForest picker with LFE AUC 0.975 (B011, transfers
to PGC) = a strong LFE-ACTIVITY detector (20x ETS/quiet day contrast, R8) but not an event timer; (R7) a single-station per-sample
PhaseNet/EQTransformer picker reaches AUC 0.85 but pick precision 0.03 — a STRUCTURAL wall (emergent, low-s/n signal at one station),
not label jitter (precise labels changed nothing); (R9) EQ/blast rejected at 0.97, in-tremor ambient noise only 0.81 (hard residual);
(R11) LONGER TEMPLATES DO NOT HELP: the LFE waveform is coherent only over ~2 s (coda decorrelates), so 4-6 s templates add incoherent
samples — this corrects my stage-5 suggestion (D3 found sigma falls with L but precision never improves: consistent). Its NEXT list
already named multi-station input as the lever. D3 now quantifies that lever: one LFE = 1.8 MAD per 3-C station, 3.5 MAD over 21
channels; >= 20 stations for k = 6. Both threads converge: NO single-station route (template or picker) yields individual LFEs.
REVISED ROUTE: (1) fleet-wide LFE-ACTIVITY catalog from the existing per-window picker (R4/R9 models, retrained with TNOISE + AE
latents) scanned over 3-C continuous data -> per-station hourly P(LFE) with EQ/blast rejection = a PNSN-independent tremor measure
(third measure for the scaling paper, margin-wide); (2) EVENT-level catalog with hypocentres only in the dense north (Puget / southern
Vancouver Island) via network matched filter seeded by Bostock/Lin families over >= 20 3-C stations, family locations with the
stage-4 machinery; (3) drop single-station per-sample picking. Recurring trap re-hit by the D3 agent: pandas datetime64[us] -> int64
/1e9 mis-scale (lfe_features DECISIONS E23) — now in memory.

## 2026-09-10 orchestrator — correction: the picker is ALREADY in production (user: "the current picks were made using this")
scripts/fleet_station.sh step 3 SCORE@40 scores every candidate with lfe_features/models/picker_broadband_pgc.joblib; adaptive rank
threshold keeps ~15k candidates/station (data/<sta>_cand_filtered.parquet: time, lat, lon, snr, p_lfe; WISH p_lfe 0.11-0.33 — absolute
scores drift low off-B011, ranks used). Candidates = envelope peaks INSIDE PNSN tremor windows (discover_nllb_pnsn_driven.py), each
inheriting the tremor window's lat/lon -> family "locations" are tremor-bin centroids. Picker selection by p_lfe for FAMILIES was found
anti-predictive (AUC 0.40) and dropped; causality is the sole verdict. The densify matched filter (the false-alarm floor) is downstream
and picker-free. Consequence: route (1) in RESULTS.md (continuous picker scan) is largely redundant and PNSN-dependent unless run over
all continuous data; per-station daily candidate counts already exist as a cheap semi-independent activity measure. The only genuinely
new product remains the north-only network event catalog with hypocentres.

## 2026-09-10 orchestrator — correction after user pushback: a network-wide catalog IS achievable via a deep-learning picker
I over-generalised the template-matching limit into "impossible at single stations". Counterexample: Lin et al. 2024 (> 1 M LFEs,
single-station DL picker + association + location, same station types). Our June per-sample picker's failure was small training set +
origin-time labels + no association, not physics. New plan E1-E3 in PLAN.md; E1 (resources) launched. Lesson recorded in memory.

## 2026-09-10 E1 Lin picker resources
Full report: `results/e1_lin_resources.md`. Downloads: `results/external/` (~1.8 GB; URLs, sizes, SHA256,
licences in report §7). Web search + curl only; nothing outside `lfe_catalog_repair/` was modified.

**Lin et al. 2024 (Seismica 2(4), doi 10.26443/seismica.v2i4.1134).**
- CODE public (Zenodo 10.5281/zenodo.10076720 = `jiunting/Cas_LFE`, HEAD ea03213) — downloaded.
  **TRAINED WEIGHTS ARE NOT PUBLIC**: no .tf/.h5/checkpoint in the repo or zip, none named in the
  availability statement. Replicating Lin means retraining.
- Model = the Thomas et al. (2021) 1-D U-Net (TF/Keras), size factor 2, sigmoid arrival probability,
  separate P and S models. Input 3-comp, 100 Hz, 15 s (1500 samples) x 6 channels =
  per-component [ln(|x|+1e-6), sgn(x)] after per-window max-amplitude normalisation;
  preprocessing detrend/taper 2%/highpass 1.0 Hz only (NO response removal). Threshold y>=0.1.
- Labels = Bostock et al. (2015) ARRIVAL PICKS (not origin times), Gaussian target sigma 0.4 s,
  30 s cut windows, equal noise class >=180 s from any LFE; ~500k waveforms per phase;
  269,422 events, 25% held out SPLIT BY EVENT ID. 30 epochs, batch 32.
- Stations: 14 train (PO+CN: LZB YOUB KLNB TWKB SNB SILB TSJB TWGB PFB GOWB NLLB PGC SSIB VGZ),
  10 unseen = C8 (BPCB GLBC JRBC LCBC MGCB PHYB SHDB SHVB SOKB TWBB); wider test = PB B004..B036.
- Association = pure coincidence, >=3 stations in the same 15 s window. Location = 1 km grid search
  (120x140x60 km), Stingray/Toomey travel times, S ONLY (P misfit sigma 4.2 s -> discarded),
  origin time removed by demeaning residuals. Velocity model = Savard et al. 2018 (`Savard_VpVs.txt`)
  — NOT publicly distributed (Wiley 403; ask authors / Savard UBC thesis 1.0371609).
- Performance: S accuracy 92%, precision ~95%, recall ~90%, AUC 0.97; S pick misfit -0.2 +/- 3.6 s;
  single-event location scatter ~22 km, depth unresolved.

**Lin CATALOG downloaded** (Zenodo 10.5281/zenodo.10016020, CC-BY-4.0, 101 MB, md5 verified):
1,058,114 events 2005-01-01..2017-02-21, columns starttime/OT/lon/lat/depth/residual/dt/N.
NO station picks (the per-station `cut_daily_*.csv` were never released). Measured caveats:
N=3 for 46% of events; residual<0.5 s (the paper's own quality cut) keeps only 41.3%;
34% of depths sit in the deepest 50-60 km bin and 15% in 0-10 km — a third of the catalog is PINNED
to the grid's bottom plane. Use as epicentres/centroids only, filtered residual<0.5 & N>=4.

**Bostock et al. 2015 catalog RECOVERED COMPLETE** ->
`results/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv`: **269,423 rows = the exact published
count**, 130 families. Source = Slow Earthquake Database, which has MOVED: sloweq at U-Tokyo is retired,
live site is **https://eqdb.sakura.ne.jp/**. No documented API; POST recipe in report §3a. Two traps:
(1) `catalogs[]` must be the button **id** (`Bostock2015-LFE`), not the value (`Bostock2015`) — the
wrong one silently returns 0 rows; (2) large result sets 500/404 — the 2004-07 and 2005-09 ETS episodes
only come down in 5-day windows via the JSON route. The DB rows carry origin time + repeated family
hypocentre + magnitude, NO family label and NO picks. Licence: **no redistribution**.

**The missing `sav_family_phases.npy` is PUBLIC after all** — SeisBench
`LFEStacksCascadiaBostock2015` (CC-BY-4.0), downloaded to
`results/external/other/seisbench_lfestacks_cascadia/`: 1,817 three-component 20 Hz stacks,
**129/130 Bostock families at 21 stations** (incl. the C8 set), with `trace_S_arrival_sample` (1453),
`trace_P_arrival_sample` (1067), family hypocentres and uncertainties, train/dev/test split.
The (lat,lon,depth) join to the sloweq origin times is verified: 128 of 130 families match at
2-decimal rounding. => Lin's training substrate is fully reconstructable from public sources.

**A better picker exists, with weights: Münchmeyer et al. 2024** (Seismica 3(1),
doi 10.26443/seismica.v3i1.1185; arXiv 2311.13971) = SeisBench `LFEDetect`, a variable-length PhaseNet.
All 8 pretrained weight sets downloaded (`other/seisbench_lfedetect/`): cascadia, mexico, san_andreas,
nankai + joint cms/jcms/jcs/jms. `cascadia.json` contract: ZNE, 20 Hz, **1200 samples = 60 s**,
bandpass 1-8 Hz zerophase, norm=peak, classes PSN, sigmoid, grouping=instrument, overlap 600.
Trained on LFE STACKS + INSTANCE noise (synthetic SNR control, 1-3 LFEs per trace), workflow
picker -> **PyOcto** -> **NonLinLoc** 1-D. The JOINT model beat the Cascadia-only model, including in
Cascadia. Weights mirror: https://seisbench.gfz.de/mirror/models/v3/lfedetect/ (primary hifis-storage.desy.de).

**Other Cascadia LFE resources.** Shelly et al. 2025 (GRL e2025GL116116) southern-edge/MTJ catalog:
61,441 detections, 27 families, 2018-2024, USGS 10.5066/P1TCKK7G, **CC0** — downloaded
(`other/LFEcat_MTJ_2018-2024_shelly2025.csv`). ⚠ **`catalogs/mtj_ducellier/LFEcat_MTJ_2018-2024.csv`
in the main tree is a 1,205-byte saved 404 HTML page, not data** — left in place (outside write scope).
Ducellier & Creager 2022 southern-Cascadia catalog (478k LFEs, 66 families, templates + travel times)
is ALREADY local and complete at `catalogs/ducellier_catalog/` (2.5 GB; Zenodo 10.5281/zenodo.6374239)
— not re-downloaded; `catalogs/ducellier_b935/` is a 34-family time-only extract. Ducellier's Olympic
Peninsula work (10.1029/2021JB022708, Zenodo 10.5281/zenodo.5725990) is TREMOR depths, not an LFE
catalog. Thomas et al. 2021 Parkfield U-Net weights ARE public (`amtseismo/Parkfield-LFE-CNN`,
Zenodo 10.5281/zenodo.4777132) — the size-2/sigma-0.4 configs Lin chose are downloaded to
`other/parkfield_lfe_cnn/` (no LICENSE file; cite the paper). **Wech/PNSN tremor detector code is NOT
public** (PNSN GitHub has only tremor_api / tremor-client). No SeisBench model other than LFEDetect has
an LFE class. Savard & Bostock 2015 BSSA e-supplement = figures only (checked).

**Tools.** NonLinLoc ALREADY INSTALLED (`~/.local/bin/{NLLoc,Vel2Grid,Grid2Time}`).
seisbench 0.12.5 and pyocto 0.2.0 install cleanly (pip dry-run verified; torch 2.12.0+cpu and
obspy 1.5.0 present). **TensorFlow is NOT installed** — needed to run Lin's/Thomas's code as-is.
GaMMA is not usably on PyPI (`gamma` 0.0.1 is a squat); install from github.com/AI4EPS/GaMMA.
pykonal 0.4.1 available as a python grid-search alternative.
Velocity: `casc1.6_velmdl.zip` (Stephenson 2017, 10.5066/F7NS0SWM, 1.43 GB zip -> **17.7 GB** ASCII)
downloaded to `results/external/velocity/`. CVM v1.7 (10.5066/P14HJ3IC, CC0, L3 = 763 MB covering
10.8-59.4 km) is **NOT scriptable** — the big netCDFs sit behind ScienceBase's JS file manager
(both download routes return a 4 KB SPA shell); fetch by hand if wanted. Only its FGDC metadata came down.

**Completeness outside VI/Puget — the answer is in the papers, and it is a warning.**
Lin §4.3 tested PB stations central->southern Cascadia: peaks align with SSE timing for the northern
B001-B013 only, and "the overall results are noisy and may include many false detections"; the 95%
precision is conditional on VI-like noise and paths, and they recommend TRANSFER LEARNING. Lin's own
catalog collapses in early 2014 "due to a lack of stations ... only three stations (NLLB, PFB, PGC)" —
3 stations is the floor at which the >=3-station association barely functions. Münchmeyer: template
matching "allows ... to identify LFEs with fewer stations than deep learning" and their DL Cascadia
catalog was "substantially smaller" than the reference; they had to EXCLUDE San Andreas entirely because
a small-aperture borehole array gave unusable locations — a direct warning for PB-borehole-dominated
Oregon/N-California. Segment expectation: replication-grade only with >=6 3-C stations within ~60 km
(southern VI / Juan de Fuca / N Olympic); workable-with-fine-tuning in Puget/central WA; Oregon and
northern California = an ACTIVITY TIME SERIES plus template matching against Ducellier/Shelly families,
not an event catalog with hypocentres. Depth is unresolved everywhere from S alone.

**RECOMMENDATION (report §8): replicate the METHOD, not the WEIGHTS.**
Picker = SeisBench LFEDetect `cascadia`, benchmarked against `jcms` (joint) — weights already on disk,
no TF install, no 269k-window training build, and Lin's own assessment is that his model is worse
off-VI. Association = PyOcto (travel-time consistent; require N>=4 outside VI), with Lin's 15 s/N>=3
rule as the sanity floor. Location = NonLinLoc (installed) + casc1.6, S-weighted, publishing
family/cluster centroids only. Fine-tune before trusting anything south of Puget, using
Münchmeyer's noise+stack synthesis on: the 129-family CC-BY-4.0 stacks, Shelly's CC0 MTJ families
(the only public southern-Cascadia ground truth with family labels), Ducellier's 66 families already
local, and our own hard negatives (lfe_features R9 EQ/blast set, the 0.81 in-tremor ambient windows,
and the zero-amplitude degenerate windows found at WISH/YACT). Keep Lin's catalog (residual<0.5, N>=4)
and the Bostock origin times as the two calibration targets, reported separately.
## 2026-09-10 E1 reviewed (LFEDetect weights public; Lin + Bostock catalogs + stacks with picks on disk; NonLinLoc/PyOcto/casc1.6 ready) -> E2 pilot launched
## 2026-09-10 F1 launched (fleet inventory) in parallel with E2

## 2026-09-10 F1 (resumed after crash) — fleet inventory & association regions
Deliverable `results/f1_fleet_inventory.md`; run log `logs/f1_fleet_inventory.log`; verification
`scripts/f1_verify.py` -> `logs/f1_verify.log` (**31/31 cross-table checks PASS**); figure
`figures/f1_fleet_inventory.png` (5 panels: station map + 3 density maps + 7 region small-multiples).
Metadata only — no waveforms downloaded. FDSN load: **one `level=channel` request per network**
(14 networks, sequential, StationXML cached in the session scratchpad; the previous session's cache
died with /tmp) -> 11,517 high-gain channel epochs, 1,241 stations in the bbox.

**State found on disk.** 6 of the 7 CSVs were from a 21:00 run of an OLDER script version
(`f1_cost_estimate.csv` had no `note` column, which the current `cost_table` adds), and
`f1_stations.csv` was 0 bytes — truncated by the 21:03 re-run that the crash killed at its first
`to_csv`. Everything below is from a clean end-to-end re-run; all 7 CSVs + the figure were rewritten.

**Fixes made to `scripts/f1_fleet_inventory.py`** (backup of the crashed version in the session
scratchpad, `f1_backup_pre_fix.py`):
1. `F1_CACHE` default pointed at the dead session's scratchpad -> repointed (cosmetic; the env var
   is what the run used).
2. `YEARS` 2010-2026 -> **2009-2026** (the project's fleet runs start 2009-01-01).
3. Open-ended channel epochs were counted as full years, so 2026 was billed as 365 days -> added
   `NOW = 2026-09-10` truncation. 2026 now shows 252 days / 137.5 effective station-years.
4. **Latent `UnboundLocalError`**: `z_any = z_any or z_only_here` sits outside the `if "Z" in comps`
   block, so a (band,loc) group at >=20 Hz with no Z channel would crash `station_year_table` ->
   `z_only_here = False` initialised. (It did not fire on this roster; it would on the first such
   station.)
5. **BUG in the cost estimate (headline number was wrong):** `measure_bytes_per_sample` divided the
   day-file size by `3 * (metadata max fs) * 86400`, but the files hold whichever band was actually
   downloaded (mostly 40-Hz BH) while metadata max fs is often 100 Hz HH -> bytes/sample was ~2x too
   low (0.73). Now reads the real sample count from the miniSEED record headers (`obspy read
   headonly`) on a 600-file random sample. Measured **1.416 B/sample at 40 Hz vs 0.932 at 100 Hz**
   (STEIM2 encodes differences, so the rate matters), so the pooled median was itself a ~30 % bias
   for a fleet that is 77 % 100-Hz HH: the table now evaluates a fitted `b(fs) = c0 + c1/fs` per
   station-year. Fleet native download went 13.2 TB (wrong) -> 21.4 TB (pooled) -> **17.6 TB**
   (per-rate model). `b(20 Hz) = 2.22` is an extrapolation and is flagged as such in the CSV note.
6. Added `fs_lo` per station-year (cheapest band that still gives 3-C >=20 Hz) and
   `raw_download_GB_lowest_band`: the `--target-fs 40` lever saves only **0.3 %** (20 of 2,437 fleet
   3-C station-years have a cheaper adequate band) — the volume is intrinsic to a 100-Hz HH fleet.
7. **Metric fix — the PLAN's locatable fraction is saturated.** `>=3` and `>=4` stations within
   **100 km** is satisfied by **every cell in every year for both rosters** (252/252 region-years at
   1.00): station count within 100 km does not discriminate anywhere in Cascadia. Added the same
   fractions at **60 km** (the radius where an LFE is plausibly above noise: D3 median Lin LFE s/n
   0.078 at PGC, ~40-50 km from its family), a 60-km azimuthal gap and borehole fraction, and
   **broadband-only** counts `n60_bb/n100_bb` (Muenchmeyer's small-aperture-borehole warning).
   New columns in `f1_regions.csv`: `frac_ge{3,4,6}_60`, `frac_ge4_60_gap180`,
   `frac_ge{3,4,6}_60_bb`, `med_n60_bb`, `med_n100_bb`, `med_az_gap_60`, `med_borehole_frac_60`.
   Figure panels b-e now plot the 60-km criterion; the 100-km saturation is stated in the caption.

**Numbers (all provisional, metadata-only).**
- Station master `f1_stations.csv`: **691 stations = 278 project + 413 non-project FDSN** (3-C at some
  point AND within 150 km of a source cell). All 278 project stations returned metadata.
  Project classes: broadband 181, short-period 59, PB borehole short-period 36, broadband-borehole 2.
  **57 project stations NEVER have a 3-C triplet >=20 Hz** (32 NC + 25 UW, all Z-only) -> unusable
  for LFEDetect's ZNE input. 190/278 already have a dv/v product, 187 a causality cert.
- `f1_station_years.csv` (12,438 rows = 691 x 18): project 3-C stations per year 70 (2009) -> 116
  (2018) -> **216 (2021)** -> 208 (2026, truncated); effective station-years 66.2 -> 100.5 -> 183.8 ->
  137.5; **total 2,167.6 (project) / 4,603.1 (project+extra)**. Orientation ZNE 4,425 / **Z12 893**
  station-years. Rates: 100 Hz 1,882 / 40 Hz 374 / 50 Hz 181 station-years (HH 1,182, EH 700, BH 555).
  The 2019-2021 growth is real installation, dominated by the post-2020 candidate list (those
  stations have never been processed by the dv/v pipeline).
- Source grid: **137 cells** (0.25 deg, Slab2 30-45 km, 40-50 N): R1 22, R2 17, R3 23, R4 22, R5 27,
  R6 11, R7 15.
- MEDIAN fleet 3-C stations per cell within 60 km, 2009 -> 2024: R1 4->7.5, R2 7->14, R3 5->16,
  R4 4->13, R5 3->9, R6 3->8, R7 3->7. Within 100 km: 6-15 (2009) -> 18-42 (2024).
  (Self-caught error: the first draft of the report quoted the pivot_table MEAN as "median";
  `f1_verify.py` now prints both and `f1_regions.med_n*` is the median.)
- Locatable fraction, **>=4 within 60 km** (fleet): 2009 R1 0.68 / R2 0.82 / R3 0.87 / R4 0.55 /
  R5 0.19 / R6 0.36 / R7 0.33 -> **1.00 everywhere from 2021** (R7 0.93 in 2025-26).
  **>=6 within 60 km** 2012 -> 2024: R1 0.23->0.95, R2 0.76->1.00, R3 0.43->1.00, R4 0.18->1.00,
  R5 0.07->1.00, R6 0.00->1.00, R7 0.07->0.93. First year the >=6/60 km bar is met in >=90 % of cells:
  R2 2017, R6 2019, R1/R3/R4/R5 2020, R7 2021 — but **2009/2010 with the non-project stations added**
  for R2/R3. Median azimuthal gap at 60 km (2019-26): R4 81 < R3 89 < R5 98 < R2 100 < R6 107 <
  R1 118 < R7 134 deg.
- Cost `f1_cost_estimate.csv`, **measured** inputs: bytes/sample per rate (above); GPU **52.6
  station-days/min per queue** = **6.93 GPU-min/station-year** (E2 snapshot 22:31: 38 runs, 1,962
  station-days, L40S, 4 concurrent queues x 2 torch threads, per-run 32-100 — the file is still being
  appended to by the running E2 pilot, so the number drifts); **1,867 S picks/station-day**.
  Fleet: **17.6 TB native download**, 9.1 TB (STEIM2, extrapolated) / 16.4 TB (float32) at 20 Hz,
  **250 GPU queue-hours = ~63 h wall-clock on one L40S with 4 queues**, 1.48e9 S picks.
  Project+extra: 36.8 TB, 532 queue-hours (~133 h), 3.14e9 picks. Assoc/locate (ASSUMED rates):
  ~123 + ~8 CPU-hours (fleet). **The download dominates the GPU by ~5x**; at 10/25/50 MB/s sustained
  the 17.6 TB is ~15/6/3 days (no download rate was measured — the E2/D3 logs carry no timing).
  Note `/home/jovyan` has ~75 TB free: the 80 GB cap that shaped the D3 streaming plan is a
  self-imposed study-folder budget, not the filesystem.
- Ground truth `f1_groundtruth.csv`: **Lin 2024 occupies only 48.07-49.32 N** (R1 675,615 /
  R2 382,499 events, nothing elsewhere) and its own quality cut (residual<0.5 **and** N>=4) keeps
  **104,071 = 9.8 %** (R1 66,539 / R2 37,532). Bostock 269,423 events, lat 48.12-48.93 N, **122
  distinct hypocentres for 130 published families** (R1 24 / R2 98 families). SeisBench
  `LFEStacksCascadiaBostock2015` = 1,817 stacks, **129 families, 21 stations, with P AND S picks**
  (R1 24 / R2 105; train 1,071 / dev 184 / test 562) — the only per-station-pick fine-tuning
  substrate. Ducellier/Plourde 37 located families 39.45-41.72 N (34 in R7, 3 south of 40 N).
  Shelly 2025 MTJ 61,441 detections / 27 families — **no coordinates in the release**, R7 assigned
  from the paper. PNSN tremor 748,002 (all regions, not independent of our own candidate windows).
  **R3/R4/R5/R6 have zero LFE ground truth = 83 of 137 cells (61 %).**

**Recommendation (report §vi, citing Lin §4.3 and Muenchmeyer 2024 as recorded in the E1 entry, not
re-derived).** Event catalog with hypocentres: **R2 over the whole span** (only region with both
geometry from 2009-2012 and all the ground truth) and **R1 from ~2018-2020**. Conditional on
fine-tuning + out-of-region validation: **R3 and R4 from 2019** (best azimuthal gaps in the margin,
broadband-dominated so the borehole-aperture warning does not apply, but zero local ground truth and
network-limited before 2019). Activity time series only: **R5, R6** (bar crossed only in 2019-20; a
mid-record order-of-magnitude completeness change would poison any rate analysis) and **R7** (worst
geometry: gap 134 deg, nearest 25.6 km — but the only southern region with public family ground
truth, so per Muenchmeyer template matching against the Shelly/Ducellier families is the better tool
there). No single margin-wide temporally-uniform event catalog is supportable; the honest fleet
product is R1+R2 (+R3/R4 from 2019) events plus a margin-wide activity series carrying the
per-region-year completeness proxy from `f1_regions.csv`.

**Caveats.** (a) metadata only — no gap/noise/glitch/response check, so availability is an UPPER
bound (this project's own QC: 12.9 % of PGC days component-short, 58 % WISH / 27 % YACT zero-amplitude
windows); (b) `frac_ge*` is a **geometry proxy, not detectability** — F1 cannot measure whether an LFE
is above noise at 60 km; (c) `klass` is per-station, not per-station-year; (d) azimuthal gap uses
epicentral azimuths of in-radius stations only, ignores depth; (e) region membership at 150 km ->
regions overlap, member counts not additive; (f) 2026 truncated at the run date; (g) b(20 Hz)
extrapolated; (h) the GPU rate is a live snapshot of a file E2 is still writing, and is a per-queue
rate; (i) assoc/locate CPU rates are assumptions; (j) no download rate measured; (k) the grid is
restricted to Slab2 30-45 km, so up-dip/Gorda LFEs have no cell.

**Needs a human decision.** (1) whether to include the 413 non-project FDSN stations (2.1x download
and GPU, but advances R2/R3 replication-grade geometry by ~8-10 years); (2) whether to spend ~9-17 TB
of shared-filesystem space instead of streaming; (3) whether R3/R4 hypocentres may be published at
all with no in-region ground truth, or only after a fine-tuned model is validated in R2; (4) the E2
pilot is still running and appending to `e2_throughput.csv` — the cost table should be regenerated
(`python scripts/f1_fleet_inventory.py`, ~4 min with the StationXML cache) once E2's throughput
numbers are final.

## 2026-09-10 E2 pilot (resumed after the 21:03 pod crash) — LFEDetect -> PyOcto -> NonLinLoc, end to end

**Environment after the crash.** The CUDA torch that lived in `/opt/conda` was wiped; it is back as
**torch 2.12.0+cu130 in `~/.local`** (user site). Verified inside the exact env of
`e2_run_pick.sh` (its `LD_LIBRARY_PATH` = the `/opt/conda/.../nvidia/*/lib` set, `PYTHONPATH` =
`/opt/conda/lib/python3.13/site-packages`): `torch.cuda.is_available()` True, device **NVIDIA
L40S**, `seisbench 0.12.5`, and a forward pass of both `LFEDetect.from_pretrained('cascadia')` and
`('jcms')` on the GPU (output (2,2,1200), finite). The cu13 wheels in `~/.local` win over the
prepended cu12 `LD_LIBRARY_PATH`, so the driver script needed no change. One package had to be
installed: **`pip install --user pyrocko`** (2026.6.2) — PyOcto's `VelocityModel1D.create_model`
needs its eikonal solver and raises `ImportError` without it. Nothing was installed into
`/opt/conda`.

**Picks finished** with `scripts/e2_run_pick_resume.sh` (same 4-queue / `--torch-threads 2` layout
and the same output paths as `e2_run_pick.sh`, logs `logs/e2_pick_resume_q{1..4}.log`): q1 =
vi2010 jcms then WISH cascadia+jcms, q2 = vi2017 jcms then YACT cascadia+jcms, q3 = B943 cascadia,
q4 = B943 jcms. `e2_pick.py` writes one parquet per SET at the end, so the two part-finished sets
were re-run whole. GPU memory with all four queues up: **846 MiB each, 3.4 of 46 GB** — no OOM
risk, `--batch-size 1024` untouched.

### Bugs found and fixed (all minimal, all commented in the code at the fix site)
1. **`val.gt.isin(...)` (e2_validate_pgc2010.py:159)** — the column is named `gt`, which shadows
   `DataFrame.gt`; `val.gt` returned the *method*. Fixed to `val["gt"]`. The same trap fired twice
   more in ad-hoc analysis this session (`s.gt == "lin_hq"` silently evaluates to `False`); the
   grep over the E2 scripts found no other attribute-style access of a method-named column.
2. **Empty QUIET set** (`e2_validate_pgc2010.py`, `e2_offgt.py`): quiet days were selected with
   `roll < p25` of the 30-d rolling tremor count. More than 25 % of days have ZERO tremor within
   40 km, so p25 == min == 0 and the strict `<` selected nothing -> the ETS/quiet ratio printed
   `nan`. Fixed to `>= p90` / `<= p25`.
3. **The Bostock catalogue time convention was calibrated by a chance-dominated estimator.**
   `e2_groundtruth.py` set `C = median(Lin OT - Bostock t)` over Bostock events with a Lin event
   within +-20 s. At 500 Lin events/day, **16,026 of 18,114 Bostock events have such a neighbour by
   chance**, so that median is pulled to ~0 whatever the truth is (it reported **C = +0.086 s**).
   The same flaw was in `e2_validate_pgc2010.py`'s per-ground-truth pick offset (median residual of
   the nearest pick within 10 s). Both replaced by the **peak of the residual histogram** over
   +-40 s in 0.25 s bins (background = median bin count > 5 s from the peak) -> written to
   `results/e2_pgc2010_offsets.csv`.
4. **Consequence of (3): the peak-excess background window contained the real peak.**
   `e2_resid_hist.py` integrated the excess around lag 0 with the background taken from
   10 < |lag| < 30 s. The Bostock peak sits at **+12.4 s**, i.e. inside that band: the background
   came out inflated and the "excess" at 0 s was **negative** (-3.8 sigma) — which I nearly read as
   "no Bostock detections" before the wide-lag scan. Now the excess is integrated around the
   MEASURED peak with the background >5 s away from it. Lin's numbers are unchanged (its peak is at
   -0.1 to -0.6 s).
5. **PyOcto 0.2.0 output contract** (`e2_associate.py`): assignments come back keyed `event_idx`
   (not `event`), the events frame has no `residual` column (residuals are per assigned pick), and
   `_convert_picks` identifies picks by POSITION in the frame it is handed — so the pick frame must
   be `reset_index(drop=True)` after the phase filter or `pick_idx` cannot be joined back. Also
   added: PyOcto's `n_picks`/`n_s_picks` are PICK counts, so the plan's ">= N stations" rule is now
   an explicit `--min-stations` post-filter (and is swept 3/4/5).
6. **Lin ground truth was empty for vi2010** (`e2_assoc_score.py`): `e2_gt_lin.parquet` is
   de-duplicated across windows, and the vi2010 window (2010-08-10..09-14) is inside the pgc2010
   window (2010-05..11), so every vi2010 row was dropped as a duplicate and `window == "vi2010"`
   selects **0 events**. Selecting on the event time instead recovers 90,523.
7. **NonLinLoc 6.04 needed five fixes** (`e2_locate_nll.py`), each found by running it:
   (a) its banner is not UTF-8 -> `subprocess(..., errors="replace")`;
   (b) **long ABSOLUTE paths corrupt the `GTFILES` wave-type field** — the model grid name came out
   as `layer.P<binary garbage>.mod` (the bytes of the following `GT_PLFD 1.0e-3`); every control-file
   path is now RELATIVE and the binaries run with `cwd=root`;
   (c) the model/time grid must **start above sea level** (`ztop = -3 km`): with a 0 km top every
   `GTSRCE` at positive elevation gave `Source point is not inside model grid` and no travel-time
   grid at all;
   (d) with `LOCGRID ... NO_SAVE` NLLoc locates every event and writes **nothing** ("9515 events
   located, 0 locations completed") -> the flag must be `SAVE`; the location grid was coarsened to
   2 km (111x111x41) so the per-event file is 2 MB, and the script deletes the grid/scatter files
   after parsing the summary;
   (e) if the time grid is shallower than `LOCGRID`, NLLoc **silently locates nothing** ("200 events
   read, 0 events located") with no error line — the model grid depth is now set independently of
   the location grid.
8. **Day-boundary pick artifact (new finding, not a script bug).** See the off-ground-truth section.
9. Minor: guard for ground truth with <50 events (vi2017 has 4-7 Lin events: the Lin catalogue ends
   2017-02-21), and a unit slip of mine in `e2_summary.py` (station-years/hour).

### Step 2 — single-station validation at CN.PGC, 2010-05-01..2010-11-01 (184 days)
`e2_validate_pgc2010.py` (log `logs/e2_validate.log`), `e2_resid_hist.py` (new, `logs/e2_resid_hist.log`).
Ground truth: Lin et al. 2024 (113,429 events in the window; **high-quality subset residual<0.5 AND
N>=4 = 10,695** = the primary target) and Bostock et al. 2015 (18,114 events, 13,837 with a PGC
stack S-P). Predicted PGC S arrival from the Crosson-1976 1-D model.

**The catalogue time bases, measured (`results/e2_pgc2010_offsets.csv`).** The picks align with
**Lin OT + ttS at -0.12 s (14.4x background, 0.25 s bins)** — a sub-sample time base, in every
distance bin separately (peak at -0.25 s for epi<25, 25-50 and >=50 km), which validates both the
Lin OTs and the 1-D model. Against **Bostock the peak is at +12.38 s** (7.5-13.6x background,
identically for both models and every threshold, and for both the `stack` and the pure-1-D
prediction). Cross-checks: (i) Lin-vs-Bostock origin times themselves show a broad +10..+16 s bump
(1.2x, no sharp peak) — independent corroboration of the sign and size; (ii) split by distance the
offset is more nearly constant WITHOUT the travel time (+23.9/+26.1/+24.6 s for epi<25/25-50/>=50 km)
than with it (+12.1/+12.1/+7.1 s), i.e. **the Slow-Earthquake-DB `Bostock2015-LFE` time is not an
origin time on our time base**; it behaves like a fixed ~24 s offset from the observed S arrival at
PGC. => **`C = +0.086 s` in `e2_groundtruth.py` is wrong; use the measured peak.** All Bostock
numbers below carry the empirical offset and are lower bounds (the residual distance dependence is
~5 s, wider than the +-3 s scoring window). *This should be resolved against the paper / the DB
maintainers before any Bostock-referenced product.*

**Recall / precision, +-3 s window, next to the time-shifted null** (`e2_validation_pgc2010.csv`):

| ground truth | thr | picks/day | recall (cascadia / jcms) | null | precision (casc / jcms) | null |
|---|---|---|---|---|---|---|
| Lin hq (10,695) | 0.1 | 901 / 1229 | 0.208 / 0.233 | 0.080 / 0.113 | 0.014 / 0.011 | 0.005 |
| Lin hq | 0.3 | 147 / 127 | **0.058 / 0.049** | 0.015 | **0.023 / 0.023** | 0.006 / 0.007 |
| Lin hq | 0.5 | 44 / 26 | 0.023 / 0.017 | 0.005 / 0.004 | 0.031 / 0.039 | 0.006 / 0.008 |
| Bostock <40 km (7,837) | 0.3 | 147 / 127 | 0.085 / 0.071 | 0.019 / 0.018 | 0.016 / 0.015 | 0.004 / 0.005 |

**Chance-free peak excess (+-1 s around the measured peak; `e2_resid_excess.csv`)** — the number to
quote, because at 150-1200 picks/day the fixed-window null is of the same order as the measurement:

| ground truth | thr | excess picks | sigma | recall_excess | precision_excess |
|---|---|---|---|---|---|
| Lin hq, cascadia | 0.1 | 1,386 | 33 | **0.130** | 0.008 |
| Lin hq, cascadia | 0.3 | 440 | 20 | 0.041 | 0.016 |
| Bostock <40 km, cascadia | 0.1 | 962 | 27 | **0.123** | 0.006 |
| Bostock <40 km, cascadia | 0.3 | 426 | 19 | 0.054 | 0.016 |
| Lin hq, jcms | 0.1 | 1,329 | 31 | 0.124 | 0.006 |

**Reading.** The picker genuinely detects LFEs at PGC — the residual peak is 20-50 sigma and the two
independent catalogues agree — but **at best ~12-13 % of catalogued LFEs are found, and only
~0.6-4 % of picks are real**; single-station precision is 3-4x chance, not the 95 % Muenchmeyer /
Lin report for their (associated) catalogues. The daily S-pick rate at PGC does track tremor within
40 km: **ETS/quiet 1.31 / 1.57 / 1.85 at thr 0.1 / 0.3 / 0.5 (cascadia), Spearman +0.44..+0.47**;
jcms 1.38 / 1.99 / 2.71, Spearman +0.53..+0.55. So a single station gives an activity index, not an
event list — exactly the split D3 reached from the template-matching side.

### Step 3 — association (PyOcto, Crosson 1-D) and its scoring
`e2_associate.py` (fixed, see bugs 5) + `e2_assoc_score.py` (new). Sweep: sets vi2010 (10 stations,
34 d, 2010-08-10..09-14) and vi2017 (8 stations, 60 d, 2017-02-20..04-21) x models cascadia/jcms x
pick thresholds 0.2/0.3/0.5 x minimum stations 3/4/5 = 36 configurations
(`results/e2_events_*.csv`, `e2_assoc_score_*.csv`, master `e2_assoc_scores.csv`,
`e2_assoc_rates.csv`, driver `scripts/e2_run_assoc.sh`, log `logs/e2_run_assoc.log`).

**vi2010, cascadia** (tolerance **2 s on ORIGIN times**; Bostock shifted by the measured +12.46 s):

| thr | N>= | events/day | ETS/quiet | rho(tremor) | recall_B | null | precision_B | null | prec_excess(+-2 s) |
|---|---|---|---|---|---|---|---|---|---|
| 0.2 | 3 | 689 | 5.0 | +0.95 | 0.125 | 0.053 | 0.078 | 0.035 | 0.049 |
| 0.2 | 4 | 288 | 12.1 | +0.94 | 0.078 | 0.025 | 0.113 | 0.040 | 0.088 |
| 0.2 | 5 | 102 | 20.7 | +0.93 | 0.040 | 0.010 | 0.158 | 0.043 | 0.145 |
| 0.3 | 3 | 280 | 9.6 | +0.93 | 0.066 | 0.023 | 0.100 | 0.038 | 0.069 |
| 0.3 | 4 | 100 | 17.3 | +0.88 | 0.036 | 0.009 | 0.149 | 0.042 | 0.129 |
| 0.3 | 5 | 30 | 34.2 | +0.87 | 0.016 | 0.003 | **0.220** | 0.043 | 0.218 |
| 0.5 | 3 | 59 | 16.2 | +0.90 | 0.020 | 0.005 | 0.132 | 0.041 | 0.119 |

jcms behaves the same at ~2-3x lower event rate for the same threshold (thr 0.2 / N>=4: 175 ev/day,
ETS/quiet 6.4, precision 0.114; thr 0.3 / N>=5: 9.7 ev/day, precision 0.230). Station-count
distribution at thr 0.3 / N>=3: **{3: 6151, 4: 2396, 5: 722, 6: 184, 7: 52, 8: 9, 9: 1}** — 65 %
of events are bare triples.
Against **Lin's high-quality subset** the same configs give recall 0.085 -> 0.007 and precision
0.038 -> 0.079 (thr 0.2/N>=3 -> 0.3/N>=5); the Lin numbers are systematically lower than the
Bostock ones because Lin's own origin times scatter (his S misfit is -0.2 +- 3.6 s), so a +-2 s
window loses part of the peak.

**vi2017** has **no ground truth** (the Lin catalogue ends 2017-02-21, 4-7 events overlap), so only
the tremor reference applies: cascadia thr 0.3, N>=4 -> **38 events/day, ETS/quiet 18.0, Spearman
+0.76**; N>=3 -> 126 ev/day, 12.0, +0.78. The 2017 ETS episode is reproduced day for day
(figure iii).

**Reading.** Association is what turns the picks into something that behaves like a catalogue: the
event rate tracks PNSN tremor at **rho +0.87..+0.95** (vi2010) and **+0.71..+0.82** (vi2017) with
**ETS/quiet 5-34x**, far above the 1.3-2.7x of the raw single-station pick rate, and it comfortably
clears D3's ">= 5" pass bar. But **only 10-22 % of the associated events coincide with a
catalogued LFE** (Bostock, +-2 s, chance-corrected) — a lower bound, since Bostock lists 18,114 of
the 90,523 Lin events in this window, so an unmatched event can still be real; and pushing precision
up by demanding 5 stations costs a factor 10 in rate and drops recall to 1.6 %.

### Step 4 — location (NonLinLoc 6.04, Crosson 1-D as LAYER cards, S 0.4 s / P 1.5 s)
`e2_locate_nll.py` (fixed, bugs 7) + `e2_loc_compare.py` (new). Grid build (Vel2Grid + Grid2Time,
P and S, GRID2D 261x84 at 1 km) takes **11 s**; locating **9,515 events took 9 min** (the whole
vi2010 cascadia thr 0.3 / N>=3 catalogue — no need to cap, `--max-events 200` was used only for the
smoke test). casc1.6 was NOT used: it is a 17.7 GB ASCII block and Muenchmeyer's released workflow
is 1-D; that is a documented deviation from the plan's "casc1.6" wording.

| quantity | vi2010 cascadia thr 0.3 (9,515 ev) | thr 0.5 (1,990 ev) |
|---|---|---|
| formal ERH / ERZ (median) | 49.8 / 20.3 km | 58.6 / 21.4 km |
| rms residual (median) | 0.09 s | 0.04 s |
| epicentre vs time-matched **Bostock** (<=2 s) | **25.3 km** (IQR 11.1-46.9), n=976 | 21.5 km (IQR 8.8-40.5), n=297 |
| epicentre vs time-matched **Lin** (<=2 s) | 45.3 km (IQR 23.4-75.0), n=2,339 | 40.4 km (IQR 18.2-66.3), n=538 |
| epicentre to **nearest Bostock family** | **6.2 km** (IQR 3.2-16.7) vs **45.9 km** for uniform random points | 5.2 km (IQR 3.0-10.9) |
| depth | median 5.8 km; 47 % < 5 km; 17 % in 25-45 km | median 5.8 km; 48 % < 5 km |
| depth vs Bostock (matched) | **-24.1 km** (IQR -30.7..-7.5) | -25.4 km |
| by station count (dh to nearest family / ERH) | N=3: 6.6 km / 57 km; N=4: 5.9 / 40; N=5: 5.3 / 26; N=6: 4.6 / 16; N=7: 3.2 / 9 | |
| NonLinLoc vs PyOcto epicentre | 12.3 km (IQR 3.3-41.9) | 25.8 km |

**Reading.** Two different statements, both true. (1) **Individual event locations are not usable**:
25-45 km from their time-matched catalogue counterpart, formal ERH ~50 km, and the depth is not
resolved at all from S-only picks on 3-4 stations — half the events are pushed to the top of the
grid (median 5.8 km where Slab2 is at 28.8 km), the classic S-only depth-origin-time trade-off, and
the depth error is a **-24 km bias**, not scatter. (2) **The epicentre CLOUD is nevertheless on the
right structure**: the median distance to the nearest of the 130 Bostock family hypocentres is
6.2 km against a **45.9 km uniform-random null** — a 7x concentration on the known LFE belt. And
the location quality is limited by the network, not the method: dh to the nearest family falls
6.6 -> 3.2 km and ERH 57 -> 9 km from 3 to 7 stations, which is exactly what F1 says the fleet can
supply (>= 6 3-C stations within 60 km in most regions from ~2017-2021).

### Step 5 — off ground truth: WISH (surface), B943 (borehole), YACT (control), 2017-2018
`e2_offgt.py` (log `logs/e2_offgt.log`), `e2_offgt_degen.py` (new). 730/727/724 station-days per
station per model, both weight sets; the only references available are PNSN tremor within 40 km and
the WISH x B943 same-source coincidence (104 km apart, 160 shared discovery source cells).

**Daily S-pick rate vs tremor** (`e2_offgt_summary.csv`; ETS = 30-d rolling tremor > p90, quiet
<= p25):

| station | model | thr | picks/day | ETS/quiet | Spearman rho | tremor days in 2 yr |
|---|---|---|---|---|---|---|
| WISH | cascadia | 0.3 | 120.1 | **1.16** | +0.083 (p 0.03) | 26 |
| WISH | jcms | 0.3 | 131.5 | 1.07 | +0.019 (p 0.6) | 26 |
| B943 | cascadia | 0.3 | 106.4 | **1.70** | **+0.361** | 109 |
| B943 | cascadia | 0.5 | 36.6 | **2.09** | **+0.391** | 109 |
| B943 | jcms | 0.3 | 42.4 | 1.49 | +0.341 | 109 |
| YACT | cascadia | 0.3 | 154.8 | 1.16 | +0.175 | 19 |
| YACT | jcms | 0.3 | 53.2 | 1.18 | +0.101 | 19 |

Only the **borehole B943** shows a clear response (ETS/quiet 1.4-2.1, rho +0.24..+0.39, rising with
threshold); at the two surface stations the pick rate is essentially flat through ETS
(1.07-1.24, rho <= 0.18). Note the reference itself is thin at WISH and YACT: only 26 and 19 days
in two years have any PNSN tremor within 40 km, against 109 at B943.

> **CORRECTED 2026-09-11 (E3 step 6a; Merlin F6).** The last sentence is the whole explanation and
> the surface/borehole reading above is WRONG. Tremor days within 40/60/80/100 km in 2017-18:
> WISH 26/62/104/134, YACT 19/44/72/94, B943 109/146/173/210 -- WISH and YACT are simply far from
> the tremor and cannot test generalisation or instrument type. Inside the 2017 Vancouver Island
> source zone the SURFACE stations respond as well as or better than the borehole (cascadia thr
> 0.3, full-record ETS/quiet within 40 km of each station): SYMB 2.26, PGC 2.25, VGZ 2.02,
> CLRS 1.99, SNB 1.81 vs the co-located borehole B011 1.78. Direct 0.6 km control, single-station
> +-1 s peak-excess recall vs Lin-hq in vi2010: PGC 0.102 vs B011 0.125 (thr 0.1) and 0.033 vs
> 0.059 (thr 0.3) = 1.2-1.8x, not categorical; the best station of all is the SURFACE station SNB
> (0.167 / 0.077). **Do not let this seed a borehole-only fleet design.**
> See `e3_offgt_attribution.csv`, `e3_offgt_ets_quiet.csv`, `e3_surface_vs_borehole.csv`.

**WISH x B943 coincidence — a DAY-BOUNDARY ARTIFACT found and removed, and what is left.**
The first run put a one-bin spike at **exactly lag 0.0 s, 5.9x background** (cascadia thr 0.3),
which is physically impossible: the shared source cells predict S(B943) - S(WISH) =
**-17.4 .. +0.6 s** (median -10.8 s). Diagnosis: of the 281 coincident pairs within +-0.25 s, **237 (84 %) lie in the first 5 s
of a UTC day** and 241 in the first 3 hours, and half of them differ by exactly 0.008 s. The
`cascadia` weights fire on the first window of a day file — **0.48 (WISH), 0.67 (B943), 0.39 (YACT)
picks per station-day in the first 5 s of the day, 5-10x the rate expected for a 60 s slice**;
`jcms` is 10-30x cleaner (0.015-0.055/day). Because every station is cut on the same UTC boundary,
this manufactures a zero-lag coincidence. `e2_offgt.py` now drops the first and last 60 s of every
day file from the coincidence test (daily rates are unaffected: the artifact is 0.1-1.0 % of picks,
2.2-7.1x the expected day-edge share for cascadia and 0.8-1.7x for jcms — column
`frac_in_day_edge_60s` in `e2_offgt_summary.csv` — and the ASSOCIATED events are clean, 0.05-0.09 %
in the first 60 s of a day against 0.07 % expected, because PyOcto demands travel-time consistency).

**With the artifact removed the coincidence is at the physically right lag but not significant.**
The largest in-band coincidence sits at **-11.5 s, 1.43x background, 3.3 sd** (cascadia thr 0.3;
background 47 +- 6 pairs per 0.5 s bin), at **-11.5 s / 2.0x / 2.0 sd** at thr 0.5 and at
**-9.5 s / 1.42x / 2.1 sd** for jcms thr 0.3 — all inside the predicted band and within ~1 s of the
predicted median (-10.8 s). But the search runs over 36 in-band bins, so 3.3 sd is only ~p 0.03
after the search penalty. Read it as **no measurable surface-borehole coincidence, with a hint at
the correct lag** — an improvement on D3, where the re-detected catalogues gave 1.3x the shifted
control at lag 0.0 with FEWER coincidences during ETS than during quiet weeks.

**Degenerate-window QC hole: it does NOT transfer to the picker** (`e2_offgt_degen.csv`).
The E2 day files contain no exact-zero samples at all (checked directly), so a zero-amplitude
window can only come from a reader-zero-filled GAP; the measurement is therefore exact from the
miniSEED headers (`headonly=True`, per-second 3-component coverage). At WISH, mean day coverage is
**0.9928** with 26 of 727 days below 99 %, **0.75 % of all 60 s windows are gap-affected**, and the
picks are **DEPLETED** in them: 0.10-0.24 % of cascadia picks and 0.03-0.04 % of jcms picks, i.e.
enrichment **0.04-0.32** (3-25x fewer than baseline), and **0.00 %** of picks sit in a fully empty
window. LFEDetect's peak normalisation plus the trained model simply do not fire on zero-filled
data — the failure mode that put 58 % of D3's raw matched-filter detections in zero-amplitude
windows is absent here.

### Step 6 — throughput: the pilot is CPU-bound, not GPU-bound (orchestrator request)
New script `scripts/e2_scaling.py` (logs `e2_profile.log`, `e2_scaling*.log`; results
`e2_profile.csv`, `e2_scaling.csv`). It re-uses `e2_pick.py`'s own code path
(`e2_pick.picks_for_day`) in separate processes, so nothing in the production script was changed and
`e2_throughput.csv` is not touched.

**Where the wall time goes (one queue, medians over 8-10 days):**

| stage | CN.PGC 40 Hz, 14.8 MB/day | PB.B943 100 Hz, 19.3 MB/day |
|---|---|---|
| obspy read | 0.48 s (48 %) | 0.12 s (9.5 %) |
| SeisBench `annotate_stream_pre` (1-8 Hz bandpass + decimate to 20 Hz) | 0.33 s (33 %) | **0.99 s (75 %)** |
| windowing + **GPU forward** | **0.12 s (12 %)** | **0.13 s (9.7 %)** |
| `picks_from_annotations` x 7 thresholds x 2 phases | 0.07 s (7 %) | 0.07 s (5.5 %) |
| total | 1.01 s | 1.31 s |
| the same `annotate` on the CPU | - | 2.30 s (vs 1.11 s with the L40S) |

So **the GPU forward is 10-12 % of the wall time** and the dominant cost is the CPU
bandpass+decimation, which scales with the NATIVE sampling rate — and F1 says **77 % of fleet 3-C
station-years are 100-Hz HH**, i.e. the expensive case. The L40S is still worth having (it turns a
2.2-2.3 s CPU forward into 0.12-0.13 s) but it leaves the pipeline CPU-bound. I/O is not the
constraint either: this pod has **1.58 TB RAM with 551 GB in page cache**, and even cold NFS reads
run at 74-109 MB/s (2.3x slower than cached).

**Concurrency scaling** (steady state; the timer starts after a warm-up day, so CUDA/cuDNN start-up
is not charged; every configuration reads DISJOINT day files so the page cache cannot flatter the
later ones; nvidia-smi sampled every 5 s):

| station | queues x torch-threads | per-queue station-days/min | steady-state aggregate | GPU util mean / max | GPU memory |
|---|---|---|---|---|---|
| CN.PGC (40 Hz) | 1 x 2 | 105.2 | 105 | 0 % / 0 % | 855 MiB |
| CN.PGC | 4 / 8 / 12 x 2 | 76.7 / 72.1 / 78.7 | 307 / 577 / 944 | <= 11 % / 32 % | 10.2 GB at 12 |
| UW.WISH (40 Hz) | 4 x 2 | 90.8 | 363 | 0 % / 0 % | 3.4 GB |
| UW.WISH | 12 x 2 | 86.5 | **1,038** | 9.3 % / 28 % | 10.2 GB |
| PB.B943 (100 Hz) | 4 x 2 | 46.7 | 187 | 7.2 % / 18 % | 3.4 GB |
| PB.B943 | 12 x 2 | 44.7 | **537** | 7.0 % / 23 % | 10.2 GB |

The per-queue rate falls only **4-5 %** going from 4 to 12 queues, so scaling is nearly linear and
**the GPU never exceeds 32 % utilisation**; at 12 queues x 2 torch threads we use 24 of the pod's 32
allowed CPUs and 10.2 GB of the 46 GB card. **Thread limits are irrelevant here**: 4 WISH queues
gave 95.5 station-days/min per queue with OMP/OPENBLAS/MKL unset (which is what `e2_run_pick.sh`
does, and OpenMP then defaults to the 176 CPUs `nproc` falsely reports) versus 90.8 with the limit
set to 2 — obspy/SeisBench preprocessing is single-threaded numpy, so the classic
176-vs-32-CPU trap does not bite this job.

**Fleet projection per weight set** (100-Hz rates, the fleet-dominant case; F1's station-year
totals):

| fleet | station-days | 4 queues (187/min) | 12 queues (537/min) | preprocessing-pool ceiling (est.) |
|---|---|---|---|---|
| project (2,168 station-years) | 791,174 | **70.6 h (2.9 d)** | **24.6 h (1.0 d)** | ~8-12 h |
| project + extra (4,603 station-years) | 1,680,132 | 149.9 h (6.2 d) | 52.2 h (2.2 d) | ~17-25 h |

**What a preprocessing worker pool would buy (estimate only — `e2_pick.py` was NOT rewritten).**
Everything except the bandpass+decimation costs 0.32 s per station-day (read 0.12 + GPU 0.13 +
post 0.07), and only 0.13 s of that touches the GPU. Putting the `pre` stage in a
`ProcessPoolExecutor` feeding one GPU process (the `densify_gnw_gpu.py` pattern) puts
28 of the 32 allowed cores on the 0.99 s/day stage = **~1,700 station-days/min of preprocessing
capacity**, with 3-4 cores enough for the GPU process (0.13 s/day) — a practical ceiling of
**~1,100-1,700 station-days/min = 8-12 h for the project fleet**, i.e. **~2-3x beyond 12 queues**,
plus it replaces 12 CUDA contexts (10.2 GB) with one. Worth doing only if the fleet pass is to be
repeated (e.g. once per weight set, or after each fine-tuning round). **The bigger lever is moving
the bandpass + decimation onto the GPU** (torch or cupy `resample_poly`): that removes the CPU
ceiling entirely and would leave the download as the only cost. Either way the download dominates:
F1 puts it at 17.6 TB, ~3-15 days at 50-10 MB/s.

**Caveat on the pilot's own number.** `e2_throughput.csv` reports 47.4 station-days/min per queue
(6,324 station-day runs in 2.22 queue-hours), but those runs shared the pod with the association sweep,
NonLinLoc and the off-GT analysis, and mixed 40-Hz and 100-Hz stations. **For planning use the
clean scaling table above.** F1's `f1_cost_estimate.csv` was built from a 52.6 station-days/min
snapshot of this file; regenerating it is F1's call (their file, not touched here).

### Throughput (see step 6 for the CPU/GPU breakdown)
`e2_throughput.csv`, 44 runs, all on the L40S: **6,324 station-day runs (3,162 station-days x two
weight sets) in 2.22 queue-hours = 47.4 station-days/min per queue** (median per run 46.6, range
32.4-100.2), i.e. 7.7 GPU-min per station-year per weight set, 1,451 P+S picks per station-day at
threshold >= 0.1 and 42-280 S picks per station-day at 0.3. Those runs shared the pod with the association/location/off-GT jobs; the
clean scaling measurement in step 6 gives 187 station-days/min at 4 queues and 537 at 12 for a
100-Hz station, and the fleet projection below is taken from there.

### Verdict against the plan's pass criteria
| criterion (source) | measured | verdict |
|---|---|---|
| ETS/quiet >= 5 (D3 pilot bar) | associated events 5.0-34.2 (vi2010), 5.5-20.0 (vi2017) | **PASS** |
| tremor correlation | rho +0.87..+0.95 (vi2010), +0.71..+0.82 (vi2017) | **PASS** |
| GPU throughput for a fleet run | 7.3 GPU-min/station-year; 66 h wall for the whole project fleet | **PASS** |
| recall/precision vs Lin & Bostock (E2 wording) | measured, and both catalogues agree: single-station recall 12-13 % of catalogued LFEs at 0.6-4 % precision; associated events 10-22 % precision at 1.6-7.8 % recall | ~~FAIL vs Lin's own 90 % recall / 95 % precision~~ **WITHDRAWN 2026-09-11 (E3 step 6b, F2): Lin's 90/95 % are held-out balanced-TEST-WINDOW numbers with an equal noise class, not continuous-record catalogue numbers, and his own published catalogue is a 15 s N>=3 coincidence catalogue of which his quality cut keeps 9.8 %. The comparison is between two different quantities.** |
| location scatter vs Bostock family hypocentres | 25.3 km per event (IQR 11-47) and 6.2 km epicentre-cloud to nearest family vs a 45.9 km random null; depth unresolved (-24 km bias, ERZ 20 km) | **FAIL for event-level few-km locations; PASS as a cloud on the right structure** |
| GOAL (PLAN header, Thomas et al. 2025 at PNW scale): event-level times AND few-km locations, family assignment, amplitude->moment, stable completeness | not met with the pretrained weights on this geometry: precision ~0.15, no depth, no family assignment, completeness set by the station count | **NOT YET** |

> **SUPERSEDED IN PART 2026-09-11 (E3).** The two PASSes below (ETS/quiet 5-34x, rho ~0.9 of the
> associated events) and the "epicentre cloud 7x better than random" statement are RETRACTED: a
> shifted-pick null reproduces both (null/observed 0.33-0.59 at the recommended N>=4 operating
> points; null nearest-family distance 8.1 km against the observed 6.2 km). The vi2010 "ETS/quiet
> 34x" is additionally a window-edge artifact -- under a full-record rolling statistic all 34 days
> of that window are ETS days and none is quiet. See RESULTS.md "E3 ... STEP 1 VERDICT".

**Bottom line.** The pipeline runs end to end and is cheap; the picker really does see Lin/Bostock
LFEs (20-50 sigma residual peaks, two independent catalogues, sub-second time base); association
turns the picks into a product whose rate tracks tremor at rho ~0.9 with a 10-34x ETS contrast; and
the epicentre cloud sits on the Bostock family belt at 7x better than random. But the
**per-event** quality with the PRETRAINED LFEDetect weights and a 10-station 2010 network is
1-2 orders below Lin et al. 2024: only 13-22 % of associated events coincide with a Bostock LFE
(a lower bound, since Bostock holds 18 k of Lin's 90 k events in the window) and individual
hypocentres are 25 km with no depth. This is the outcome Muenchmeyer (their
Cascadia DL catalogue was "substantially smaller" than the reference) and Lin (§4.3, "recommend
transfer learning") both predicted for out-of-the-box weights.

### What I would do next (needs a human decision)
1. **FINE-TUNE before any fleet run.** The substrate is already on disk (E1): 1,817 SeisBench
   Bostock stacks with P AND S picks for 129/130 families at 21 stations (CC-BY-4.0), Shelly's CC0
   MTJ families, Ducellier's 66 local families, plus our own hard negatives (lfe_features R9
   EQ/blast, in-tremor ambient) and now a *third* negative class the pilot identified: the
   first-window-of-day artifact. Münchmeyer's own noise+stack synthesis recipe applies directly.
2. **cascadia vs jcms**: at equal event rate they are within noise of each other (thr 0.2/N>=4
   cascadia: 288 ev/day, precision 0.113; jcms thr 0.2/N>=4: 175 ev/day, 0.114; jcms needs a lower
   threshold for the same rate). cascadia has ~20 % higher single-station recall at a matched
   threshold; jcms is 10-30x cleaner on the day-boundary artifact and slightly better at high
   precision. **Recommendation: run both in the pilot phase (it doubles a cheap cost), pick after
   fine-tuning.**
3. **N>=3 vs N>=4**: N>=3 is 65 % bare triples with precision 0.10 and ERH 57 km; **N>=4 at pick
   threshold 0.2-0.3 is the operating point** (100-290 events/day, ETS/quiet 12-17, precision
   0.11-0.15, ERH 40 km). N>=5 buys precision 0.16-0.22 but at 1/10 the rate and 1.6 % recall.
4. **Do not publish event-level depths from S-only picks.** Either add the P picks properly
   (LFEDetect emits them; Lin discarded P at sigma 4.2 s, but our S-only depths are biased -24 km)
   or fix depth to Slab2 and solve for epicentre + origin time only.

### Files written by the resumed run (nothing deleted; `e2_scaling_raw_multiheader.csv` is the raw
### append-log of `e2_scaling.csv`, kept because the three runs wrote three different schemas)
**Scripts (new):** `e2_run_pick_resume.sh`, `e2_resid_hist.py`, `e2_assoc_score.py`,
`e2_run_assoc.sh`, `e2_loc_compare.py`, `e2_offgt_degen.py`, `e2_scaling.py`, `e2_summary.py`,
`e2_figures.py`. **Scripts (fixed):** `e2_validate_pgc2010.py`, `e2_associate.py`,
`e2_locate_nll.py`, `e2_offgt.py`.
**Picks:** `e2_picks_{vi2010,vi2017}_jcms.parquet`, `e2_picks_offgt_{WISH,YACT,B943}_{cascadia,jcms}.parquet`.
**Validation:** `e2_validation_pgc2010.csv`, `e2_pgc2010_offsets.csv`, `e2_pgc2010_daily.csv`,
`e2_resid_hist.csv`, `e2_resid_excess.csv`.
**Association:** `e2_events_{vi2010,vi2017}_{cascadia,jcms}_t{0.2,0.3,0.5}.csv` (+ `_assign.parquet`),
`e2_assoc_score_*.csv` (36), `e2_assoc_scores.csv`, `e2_assoc_daily_*.csv`, `e2_assoc_rates.csv`,
`e2_assoc_locmatch_*.csv`, `e2_pyocto_crosson` (velocity model).
**Location:** `e2_events_vi2010_cascadia_t{0.3,0.5}.nll.csv`, `e2_loc_compare_vi2010_cascadia_t{0.3,0.5}.csv`,
`nll_e2_events_*/` (control files, travel-time grids, per-event .hyp; the per-event probability
grids were deleted by the script after parsing).
**Off ground truth:** `e2_offgt_summary.csv`, `e2_offgt_daily.csv`, `e2_offgt_coincidence.csv`,
`e2_offgt_xcorr_*.npz`, `e2_offgt_degen.csv`, `e2_offgt_degen_days.csv`.
**Throughput:** `e2_throughput.csv` (appended), `e2_profile.csv`, `e2_scaling.csv`,
`e2_picks_inventory.csv`, **`e2_summary.txt` (every number above, regenerable via `e2_summary.py`)**.
**Figures:** `e2_pgc2010_validation.png`, `e2_pgc2010_exampleday.png`,
`e2_assoc_rate_vs_tremor.png`, `e2_location_scatter.png`, `e2_offgt_ets_quiet.png`.

### Loose ends
* `e2_assoc_daily_vi2010_cascadia_t0.3.csv` (no `_nN` suffix) is a duplicate of the `_n3` file from
  the first scoring pass, left in place.
* `e2_groundtruth.py` still contains the chance-dominated `C` estimator (bug 3); the downstream
  scripts no longer use its value, but the file should be fixed before anyone else reads
  `e2_bostock_lin_dt.csv` as a calibration.
* casc1.6 was not used for the locations (1-D Crosson instead); a 3-D NonLinLoc run needs the
  17.7 GB ASCII block converted to NLL grids first.
* vi2017 has no ground truth at all (the Lin catalogue ends 2017-02-21) — for a second validated
  window we would need the 2011-2013 Bostock years, which ARE in the catalogue on disk.

## 2026-09-11 ~00:30 Merlin audit of the E2 verdict -> results/merlin_audit_e2_2026-09-11.md (verbatim)
Scorecard inverted: the FAILs (recall/precision, 25-km locations) are mostly benchmark artifacts (single-station vs network catalogs;
recall monotonic in magnitude/distance = SNR ceiling; precision scored vs incomplete catalogs when quiet/ETS bounds it >=36-50%;
uniform-box location null; constant Bostock offset); the PASSes (ETS/quiet, rho of associated events) are UNVERIFIED — PyOcto with
3 S picks vs 4 unknowns always fits, and a crude Poisson estimate of chance triples reproduces the quiet-day event rate. Fine-tuning on
the SeisBench stacks = retraining on the cascadia training set with a family-leaky split, and cannot lift an SNR-limited recall.
Off-GT WISH/YACT "flat" = distance to tremor, not generalization (surface CN stations inside the 2017 source zone respond 2.6-4.2x).
GOAL mismatch: Thomas-style dynamics need family-assigned times + relative amplitudes, not event hypocentres. Prescribed E3 (on-disk,
no downloads): (1) shifted-pick chance-association null [decisive]; (2) relocate the 129 Bostock families through the same NLL from
stack picks; (3) injection recall-vs-SNR + false-alarm rate on real quiet days; (4) Bostock templates as CLASSIFIER on associated
events + MAD at Bostock event times; (5) per-family Bostock time base; (6, conditional) S-probability back-projection; (7) decision gate,
re-target the product, first science milestone = replicate Thomas's Sept-2012 SSE on S. Vancouver Island (tens of GB), fleet after.
User go 2026-09-11 ~00:45 -> E3 agent launched (steps 0-6: + Münchmeyer numbers fetch, + doc corrections / F1 at 40 km).

## 2026-09-11 E3 steps 0-6 — the on-disk controls Merlin prescribed (one Opus agent, no downloads)

Environment: torch 2.12.0+cu130 in `~/.local`, `torch.cuda.is_available()` True on the **NVIDIA
L40S**, seisbench 0.12.5, pyocto 0.2.0 — verified inside the exact env exports of
`scripts/e2_run_pick.sh` before any GPU step. Nothing installed. Concurrency kept at or below the
28-worker cap (4 PyOcto runs x 7 threads; 4 NLLoc; one GPU process at 4 torch threads).
**Only fetch of the whole run: the arXiv PDF.** No waveforms. No E2 output overwritten (everything
new is `e3_`-prefixed; the two E2 scripts that were modified keep their old behaviour by default).

### STEP 1 — chance-association null. VERDICT FIRST: the E2 activity PASSes are RETRACTED.

`e2_associate.py` gained `--shift-null SEED` (per-station constant offset, |shift| uniform in
[300, 3000] s, sign random, independent per station; picks that cross a UTC day boundary stay in
the set, so the pick count and every station's daily rate — hence the ETS modulation — are
preserved exactly while real inter-station coincidence is destroyed). Driver
`scripts/e3_run_null.sh`, 25 runs (vi2010 cascadia 0.2/0.3/0.5 and jcms 0.3, vi2017 cascadia 0.3,
5 seeds each), ~40 s per run. Scoring `scripts/e3_null_score.py` -> `results/e3_nullevents_rates.csv`.

**Merlin F10 fixed first, and it mattered.** ETS/quiet is now defined from the 30-d centred rolling
PNSN tremor count within **60 km of the network centroid**, built on the **FULL record**
(2010-01..2026-05, 5,969 d, `min_periods=30`) with the p90/p25 thresholds taken as percentiles of
that full series. Under that definition **all 34 days of the vi2010 window are ETS days and none is
quiet** — the window was chosen as the 34 densest Lin days. E2's "ETS/quiet 34x" for vi2010 was
therefore an artifact of `min_periods=5` inside the window, exactly as F10 said; there is no quiet
population in vi2010 at all and the ratio is undefined, not 34. vi2017 (60 d) does have both:
**35 ETS, 4 quiet, 21 intermediate**.

| set | model | thr | N>= | observed ev/d | null ev/d (mean +- sd, 5 seeds) | null/obs |
|---|---|---|---|---|---|---|
| vi2010 | cascadia | 0.2 | 3 | 669.5 | 616.0 +- 3.7 | **0.92** |
| vi2010 | cascadia | 0.2 | 4 | 279.1 | 163.2 +- 3.5 | **0.59** |
| vi2010 | cascadia | 0.2 | 5 | 99.1 | 22.6 +- 0.6 | 0.23 |
| vi2010 | cascadia | 0.3 | 3 | 271.9 | 199.1 +- 2.3 | **0.73** |
| vi2010 | cascadia | 0.3 | 4 | 96.1 | 31.8 +- 1.2 | **0.33** |
| vi2010 | cascadia | 0.3 | 5 | 27.7 | 2.43 +- 0.35 | 0.088 |
| vi2010 | cascadia | 0.5 | 3 | 58.5 | 24.5 +- 0.7 | 0.42 |
| vi2010 | cascadia | 0.5 | 4 | 15.0 | 1.55 +- 0.17 | 0.103 |
| vi2010 | cascadia | 0.5 | 5 | 3.50 | 0.059 +- 0.029 | 0.017 |
| vi2010 | jcms | 0.3 | 3/4/5 | 114.6 / 33.6 / 8.97 | 66.6 / 6.69 / 0.28 | 0.58 / 0.20 / 0.031 |
| vi2017 | cascadia | 0.3 | 3/4/5 | 126.4 / 36.2 / 9.45 | 88.8 / 10.5 / 0.70 | 0.70 / 0.29 / 0.074 |

**ETS vs quiet, separately (vi2017 thr 0.3, the only window with both populations).**

| N>= | obs ETS | null ETS | obs quiet | null quiet | obs ETS/quiet | **null ETS/quiet** |
|---|---|---|---|---|---|---|
| 3 | 196.0 | 139.6 | 36.0 | 25.8 | 5.44x | **5.42x** |
| 4 | 57.9 | 17.5 | 5.50 | 1.40 | 10.53x | **12.49x** |
| 5 | 15.3 | 1.17 | 0.75 | 0.10 | 20.42x | **11.66x** |

That table is the decisive one: **the shifted-pick null reproduces the ETS/quiet contrast, ratio
for ratio.** It has to — a per-station constant shift preserves each station's daily rate, and an
N-fold chance coincidence rate goes as the product of those rates, so an activity contrast of
r per station becomes r^N in the event catalogue whether or not a single real LFE is present. The
"ETS/quiet 5-34x" and "rho 0.87-0.95 with tremor" PASSes therefore **carry no information about
whether the events are real** and are retracted as evidence.

**Sensitivity variant (not in the prescription; added because the branch verdict turns on this
ratio).** A per-station CONSTANT shift also de-aligns the co-elevation of station rates inside a
tremor burst, which would *deflate* the chance estimate. `--shift-mode station_hour` reproduces
Muenchmeyer's own control (a fresh shift per station per hour, |shift| in [60, 300] s >> the ~30 s
maximum inter-station travel-time difference across this network, so coincidence is still
destroyed but the burst envelope stays aligned at the few-minute scale). 15 further runs
(`scripts/e3_run_null_hourly.sh`, `e3_nullhr_rates.csv`): the null rises by only **4-14 % relative**
— vi2010 thr 0.3 N>=4 goes 0.330 -> **0.365**, thr 0.2 N>=4 0.585 -> 0.611, vi2017 thr 0.3 N>=4
0.291 -> 0.318. The verdict is robust to the shift design.

**BRANCH VERDICT (audit rule: >=50 % at N>=4 -> chance-dominated; <=15 % -> PASS stands; between
-> unresolved, report the ratio).** It is threshold-dependent and must be stated that way:
* **thr 0.2, N>=4 (one of the two configurations E2 recommended as "the operating point"):
  null/obs = 0.59 +- 0.01 (0.61 under the hourly null) -> CHANCE-DOMINATED.**
* **thr 0.3, N>=4 (the other recommended operating point): 0.33 +- 0.01 (0.37 hourly) ->
  UNRESOLVED by the audit's rule; one associated event in three is chance.**
* thr 0.3 N>=5 (0.088), thr 0.5 N>=4 (0.103) and thr 0.5 N>=5 (0.017) **PASS**, at 28 / 15 / 3.5
  events per day instead of 96-279.
* Every N>=3 configuration is 0.42-0.92: **the N>=3 catalogue is chance, full stop.**
So: the E2 operating point was in the chance-dominated-to-unresolved zone, and there IS a clean
operating point — **thr 0.5 with N>=4, or thr 0.3 with N>=5** — at 6-11 % chance and ~15-28
events/day. That is the concrete change to the recommendation.

Independent corroboration of the mechanism: fitting log(null ev/d) against log(network S picks/d)
over the three thresholds gives a chance exponent **p = 2.16 / 3.11 / 3.97 for N>=3/4/5**
(`e3_F_scaling.csv`) — pure chance predicts exactly N. F1's crude Poisson estimate was right.

**Located null (Merlin F4).** The seed-1 vi2010 cascadia thr-0.3 null events (7,078) went through
`e2_locate_nll.py` unchanged (9 min) and then `e2_loc_compare.py`, which gained `--null-nll` so the
shifted-pick null replaces the uniform-over-the-bounding-box null (the legacy number is still
printed, labelled wrong).
**Observed epicentres: median 6.2 km to the nearest of the 130 Bostock family hypocentres.
Shifted-pick chance events: 8.1 km. Ratio 0.76.** The uniform-box null said 45.9 km, i.e. a "7.4x
concentration on the LFE belt". **That statement is retracted**: chance events land essentially as
close to Bostock families as the real catalogue does, because S-only 3-station solutions collapse
toward the network centroid, which is where the families are. The null's depth distribution is the
same pathology too (median 4.5 km, 57 % shallower than 5 km, vs observed 5.8 km / 47 %).

Files: `e3_nullevents_{vi2010,vi2017}_*_s{1..5}.csv` (+ `_assign.parquet`), `e3_nullhr_*`,
`e3_nullevents_rates.csv`, `e3_nullhr_rates.csv`, `e3_nullevents_vi2010_cascadia_t0.3_s1.nll.csv`,
`e3_loc_null_vi2010_cascadia_t0.3_e3famoff.csv`, `e3_loc_compare_vi2010_cascadia_t0.3_e3famoff.csv`.
Figure `figures/e3_null_rates.png`.

### STEP 2 — family relocation control. The prescription's premise was wrong; here is what replaced it.

**Data fact found before writing any code, and it invalidates the prescribed method.** In
`seisbench_lfestacks_cascadia/metadata.csv` (1,817 traces, 129 families, 21 stations, 20 Hz),
`trace_P_predicted_arrival_sample` is **exactly 800.0 for all 750 traces that carry it** and
`trace_P_arrival_sample` is **exactly 800.0 for all 1,067 that carry it**. Every trace is
individually aligned on its own P arrival; there is no trace start time anywhere in the dataset.
**Cross-station moveout has been destroyed by the alignment** — the only timing information that
survives per (family, station) is the S-P differential time (S_sample - 800)/20 s (it behaves:
family 159 gives 3.95 s at GLBC, 19.4 km, and 6.35 s at LZB, 46.3 km). Writing the prescribed NLL
obs files ("P at a per-family origin + 40 s, S at 40 + (S-P)") would impose the false constraint
that every P arrival is simultaneous, i.e. that the family is equidistant from all stations, and
would drag every solution to the station centroid. Merlin was consulted with the evidence and
re-specified the control; the four configurations below are his.

**Config A — model + method (real stack picks, S-P differential grid search, `e3_family_relocate.py`).**
Own grid search on the same Crosson-1976 1-D model `e2_locate_nll.py` uses (via
`e2_ttable_crosson.npz`), L1 misfit, 1 km horizontal over +-60 km, 1 km depth 5-70 km, origin time
cancels identically. Only the **672 both-picked traces** are used (P hand-picked AND S hand-picked;
the 750 P-*predicted* traces are Bostock's own model prediction at the published hypocentre and
would be circular), **SHVB excluded** (Merlin: its metadata coordinates disagree with
`path_ep_distance_km` by up to 19 km where all 20 other stations agree to < 2 km). A
vertical-incidence elevation correction elev x (1/Vs0 - 1/Vp0) is added because the table puts
receivers at z = 0 while NLL's GTSRCE puts them at their true elevation.

| config | subset | families | median dh (IQR) | median dz (IQR) | L1 resid |
|---|---|---|---|---|---|
| A (all stations, S-P) | all | 102 | **4.12 km** (2.24-5.83) | -3.23 km (-4.68..-1.62) | 0.06 s |
| A | n>=6 stations | 53 | 4.00 km (2.24-5.10) | -3.28 km | 0.08 s |
| A | n>=6 and gap < 200 deg | 40 | **3.08 km** (2.24-4.21) | -3.57 km | 0.08 s |
| B (6 vi2010 stations, S-P) | all | 24 | 3.00 km (2.18-5.21) | -3.89 km | 0.02 s |

**A passes the audit's bar** (dh <= 5 km and |dz| <= 5 km on the headline subset): the Crosson 1-D
model plus this machinery can put a family within ~3 km horizontally from publication-quality
picks. The systematic -3.2 to -3.9 km depth bias is exactly what Merlin predicted from his own
measurement that Crosson over-predicts S-P by ~0.31 s at the published hypocentres (effective
Vp/Vs ~1.757 rather than 1.78). **So the model and the method are NOT the 25 km.**

**Config C — geometry + S-only, through the EXACT `e2_locate_nll.py` control (`e3_geom_synth.py`).**
Synthetic, and labelled as such: the sources are the 129 published hypocentres and the S times are
read out of **NonLinLoc's own travel-time grids** (`time/layer.S.<STA>.time.buf`, GRID2D, bilinear
in (horizontal distance, depth)), so the model term is identically zero and only geometry plus the
injected pick noise remains. All 42 three- to six-station subsets of PGC/LZB/VGZ/YOUB/SNB/PFB, 2
realisations, four noise levels, 10,836 NLL locations each, 4 runs in parallel (~35 min).

| sigma | N=3 dh / dz / frac<5 km | N=4 | N=5 | N=6 |
|---|---|---|---|---|
| 0.05 s | **17.7 km / -27.6 / 0.33** | 0.98 / -0.77 / 0.01 | 0.52 / -0.37 | 0.43 / -0.32 |
| **0.30 s** (measured) | **18.0 km / -27.8 / 0.34** | **3.73 / -1.33 / 0.03** | 2.49 / -0.43 | 2.12 / -0.04 |
| 0.50 s | 18.2 / -27.8 / 0.34 | 6.04 / -1.86 / 0.05 | 4.21 / -0.78 | 3.75 / -0.26 |
| heavy tail | 18.6 / -28.0 / 0.35 | 5.72 / -2.01 / 0.07 | 3.90 / -1.20 | 3.19 / +0.15 |

The N=3 row is **independent of the pick noise** (17.7 km at 0.05 s, 18.6 km at heavy tail): a
3-station S-only solution has 3 observations for 4 unknowns, the depth/origin-time trade-off is
unconstrained, a third of the solutions are pinned to the top of the grid and the depth error is
-28 km. **E2's N=3 depth bias of -27 km and its shallow pile-up are 100 % geometry**, reproduced
here with a perfect model. But at **N>=4 the geometry is fine** — 3.7 km at the measured 0.3 s pick
noise, against E2's 25 km. Geometry does **not** explain E2 beyond N=3.

**Config D — real picks, known answer, no PyOcto (`e3_realpick_locate.py`).** Merlin's
discriminator for exactly that gap. For every Bostock event the predicted S arrival at each station
is (DB time) + (per-family offset from step 5) + Crosson ttS; any real LFEDetect S pick within
+-1.0 s at a station with epi < 40 km is taken as that event's pick, and events with >= 3 such
stations go through the unmodified `e2_locate_nll.py`.
Measured **real per-station S-pick residual sd inside the +-1 s window: 0.40-0.56 s** (PGC 0.44,
LZB 0.47, VGZ 0.46, YOUB 0.55, SNB 0.48, PFB 0.51 at thr 0.1) — confirming Merlin's 0.3-0.5 s and
killing the "1-3 s pick scatter" premise I had started from.

| | N=3 | N=4 |
|---|---|---|
| C, synthetic sigma 0.3 s | 18.0 km / -27.8 km | 3.73 km / -1.33 km |
| **D, REAL picks, no PyOcto** (thr 0.1, 124 ev) | **6.97 km** (4.5-10.3) / **-30.5 km**, 50 % < 5 km | 1.85 km (n=1) |
| E2, PyOcto-associated, Bostock time-matched | 33 km / -27 km | 25 km / -20 km |

Reading the ladder: **C -> D** (real, non-Gaussian pick error) costs nothing at N=3 — D is in fact
*better* than C because D's three stations are always the near ones (epi < 40 km) while C averages
over all 3-subsets including poor geometry. **D -> E2** is where the factor ~5 lives. Since D uses
the same picker on the same data with the same locator and differs only in that it does not go
through PyOcto, **the dominant term in E2's event-location error is the ASSOCIATION, not the model,
not the geometry (beyond N=3), and not the picker.** That is the same conclusion step 1 reaches from
the rate side, from completely different evidence.
Caveat recorded: config D yields only **124 events at N>=3 and 1 at N>=4** at thr 0.1 (19 / 0 at
thr 0.3). Merlin's predicted 698/169/27 came from selecting picks on `thr_first <= 0.3` alone,
which is the bug described below and is effectively the thr-0.1 pick set.
Files: `e3_family_reloc{,_summary}.csv`, `e3_geom_synth_{s005,s030,s050,sheavy}.csv`,
`e3_realpick_{locate,summary}_t{0.1,0.3}.csv`. Figure `figures/e3_relocation.png`.

### STEP 3 — injection recall-vs-SNR and the false-alarm rate on real noise (`e3_inject.py`, GPU)

Bostock stacks recorded AT the host station, band-passed 1-8 Hz, scaled to peak/host-rms ratios
0.1/0.2/0.3/0.5/1/2, resampled to the host rate and added into real quiet 3-C day files at known
random times: **300 injections per station per level, 6 host days, both weight sets**, and the same
day picked un-injected so the baseline (chance picks plus any residual real LFE) is measured and
subtracted rather than assumed. Host days: PGC the 6 quietest 2010-05..07 days (0 Bostock, 0-4 Lin
events, 0 PNSN tremor within 60 km); LZB and YOUB the quietest days on disk, which all sit inside
the 2010-08..09 ETS window because only the 34-day vi2010 window was downloaded — their F is an
upper bound and their baseline subtraction is doing real work. Only 2 of the 92 PGC May-July days
have literally zero Bostock AND zero Lin events (Lin lists ~616 events/day), so the
"zero-event day" criterion in the brief is unusable as written and the measured-baseline design
replaced it.

Two SNR measures per injection: `snr_peak` (the knob) and **`snr_rms`** = rms of the injected signal
over the 6 s from S-1 s divided by the rms of the 1-8 Hz host — which is the definition behind D3's
`rms(signal)/rms(null) = sqrt(1 + N (s/n)^2)`, i.e. the one the project's s/n = 0.078 is quoted in.
**The conversion used for the real curve, stated:** s/n(M) = 0.078 x 10^(1.5 (M - M_ref)) with
M_ref = 1.650 = the median magnitude of the 8,554 Bostock events with epi_PGC < 40 km (Bostock's
magnitudes are moment magnitudes from template amplitudes, so amplitude scales as the moment). Two
caveats carried into the figure: the 0.078 anchor is a median *Lin* event taken as the median
*Bostock* event of the same population, and D3's mixed-source stack gave 0.032 as a lower bound, so
the anchor is uncertain by a factor ~2.4 (drawn as a band).

| PGC cascadia, thr 0.1 | injected s/n_rms | injected recall +-1 s | real s/n_rms | real recall |
|---|---|---|---|---|
| level 0.1 | 0.019 | 0.000 | 0.038 (M 1.44) | 0.109 |
| level 0.2 | 0.037 | 0.003 | 0.055 (M 1.55) | 0.142 |
| level 0.3 | 0.057 | 0.003 | 0.070 (M 1.62) | 0.163 |
| level 0.5 | 0.091 | 0.017 | 0.090 (M 1.69) | 0.185 |
| level 1.0 | 0.187 | 0.075 | 0.118 (M 1.77) | 0.193 |
| level 2.0 | 0.387 | **0.653** | 0.185 (M 1.90) | 0.256 |

**Neither of the audit's two branches fires: the injection is far BELOW the real curve, not above
it** (at matched s/n_rms ~0.09, real 0.185 vs injected 0.017, a factor ~11; at ~0.19, real 0.19-0.26
vs injected 0.075). Matching the two curves says a real Bostock LFE behaves like s/n_rms ~0.22-0.25,
i.e. **~3x the 0.078 anchor**.

**The obvious explanation was tested and rejected.** Real LFEs come in bursts, so LFEDetect's 60 s
window usually holds several of them while the injections are isolated events in quiet noise.
Splitting the real population by how many other Bostock events (any family, epi_PGC < 40 km) lie
within +-30 s (`e3_real_recall_by_burst.csv`) gives, in matched magnitude bins,
**isolated (<=1 neighbour) 0.128 / 0.171 / 0.197 vs moderate (2-9) 0.124 / 0.181 / 0.233 vs dense
(>=10) 0.100 / 0.114 / 0.279** — the isolated events do as well as the dense ones. Burst context is
NOT the explanation. What is left is (i) the s/n = 0.078 anchor being low (D3 itself flags it as one
station, one cell, and notes alignment jitter biases such an estimate downward) and/or (ii) the
family STACK being a pessimistic proxy for a single event: stacking hundreds of jittered events
smears the waveform, so at equal 6 s rms the stack has a lower peak and a less impulsive onset than
a real LFE. **Operational consequence, which is the half of Merlin's branch that does apply: the
injection test cannot be used to claim the picker is at its SNR ceiling, and any fine-tune must be
evaluated on REAL single-event positives, not on the stacks.** (Caveat, as the audit required: the
stacks are in-sample for the `cascadia` weights, so the injected recall is if anything optimistic —
which makes the shortfall against the real curve more striking, not less.)

Pick-time bias and scatter at PGC cascadia thr 0.1, levels where >10 injections are recovered:
bias +0.03 s at level 1.0 and 0.00 s at level 2.0, scatter (16-84 pct half-width) 0.29 s and 0.18 s
— consistent with the 0.40-0.56 s real residual sd measured independently in config D.

**False-alarm rate F on real quiet days** (`e3_inject_falsealarm.csv`, S picks per station-day,
cascadia): PGC 720 / 247 / 113 / 32 at thr 0.1 / 0.2 / 0.3 / 0.5; LZB 393 / 61 (thr 0.3 on an
ETS-window day, an upper bound); YOUB 553 / 97. **How far must F fall for the N>=4 chance rate to
be <= 10 % of the real rate?** With the chance exponent p fitted above and real = observed - null:
**thr 0.2 N>=4 needs F x 1/2.34; thr 0.3 N>=4 needs F x 1/1.67; thr 0.5 N>=4 is already there
(1.05).** So the requirement is modest — a factor ~2 — and it is *already met* by raising the pick
threshold, at the cost of rate. Fine-tuning would have to buy that factor without losing the recall
that the threshold increase costs; nothing here shows it can.

### STEP 3b (user request, 2026-09-11) — does DENOISING before picking help? NO. (`e3_denoise.py`)

`DeepDenoiser` `original` (Zhu et al. 2019, MIT) was **already in `~/.seisbench/models/v3/deepdenoiser`
— no download was needed**. Three arms, everything else identical, **all arms at 100 Hz**
(DeepDenoiser's native rate; the 40 Hz hosts are resampled once, up front, for every arm, so no arm
sees a different rate) and the band-pass is the one LFEDetect applies itself:
`raw` | `dn` (DeepDenoiser per component) | `wiener` (non-ML baseline: band-limited spectral
soft-mask, noise PSD = the quietest decile of the day's own 30 s window PSDs, applied in 1-8 Hz).
The same injections as step 3 were replayed (same station, day, S time, family stack and level).

**(iv) What the denoiser does to the injected LFE**, measured exactly as D = arm(host+signal) -
arm(host) against the known injected signal:

| arm | rms kept | waveform corr | spectral ratio 1-2 Hz | 2-4 | 4-6 | 6-8 |
|---|---|---|---|---|---|---|
| raw | 1.000 | 1.000 | 1.00 | 1.00 | 1.00 | 1.00 |
| **dn (DeepDenoiser)** | **0.071-0.117** | **0.44-0.60** | **0.04** | **0.06** | 0.16 | 0.29 |
| wiener | 0.92-0.94 | 0.96-0.99 | 0.92 | 0.92 | 0.93 | 0.94 |

**DeepDenoiser removes 88-93 % of the LFE's amplitude and strips it hardest exactly where LFE
energy lives (96 % gone at 1-4 Hz, 71 % at 6-8 Hz)** — it treats the emergent, low-frequency,
low-SNR LFE as noise. It is trained on impulsive local earthquakes and this is the expected failure,
now measured.

**(i) injected recall** (cascadia thr 0.3, level 2.0, chance-corrected): LZB raw 0.459, wiener
0.499, **dn 0.000**; PGC raw 0.003-0.023, dn 0.023 with a baseline of 0.023 (i.e. nothing).
**(iii) false alarms** (S picks per quiet station-day, cascadia thr 0.3): PGC raw 113, wiener 120,
**dn 956 (8.5x worse)**; LZB raw 61.5, wiener 65.7, dn 244.7; YOUB raw 97.3, wiener 100.5, dn 150.2.
**(ii) REAL dense-ground-truth PGC 2010 days** (2010-08-30/09-04/08-28/08-21, 701-807 Lin-hq events
each), +-1 s peak-excess recall at thr 0.1: vs **Bostock** raw 0.134, wiener 0.131, **dn 0.034**;
vs **Lin-hq** raw 0.111, wiener 0.111, **dn 0.010**. Precision moves the same way
(Bostock thr 0.3: raw 0.208, wiener 0.197, dn 0.012).

**Answer: denoising does not help.** DeepDenoiser is actively harmful on both axes at once — 4-11x
worse recall AND 2-8x more false alarms — because it deletes the signal and leaves reconstruction
artifacts the picker fires on. The classical spectral-Wiener baseline is **neutral**: it preserves
the waveform (corr 0.96-0.99) and lands within a few per cent of raw on recall, F and precision,
which is the expected result for a linear soft-mask applied inside the same band the picker already
filters to. No LFE-appropriate denoiser weights exist in SeisBench; a denoiser for this problem
would have to be trained on LFEs, and step 3's finding (the limit is the signal amplitude, not
additive broadband noise the mask can remove) says the headroom is small.
Files `e3_denoise_{recall,falsealarm,fidelity,realdays}.csv`, figure `figures/e3_denoise.png`.

### STEP 4 — Bostock templates as a CLASSIFIER (`e3_template_classifier.py`, GPU)

Network cc of the 3-C data at the 6 overlap stations against each of the 129 family templates
(the family's own stack, the 8 s segment from S-1 s to S+7 s, placed at Crosson ttP + the stack's
own S-P so the moveout is the family's; 3 components concatenated and L2-normalised; 20 Hz, 1-8 Hz;
one COMMON origin-time lag scanned over +-3 s at 0.05 s, so this tests network coherence and not 6
independent single-station correlations). Implemented as a grouped `conv1d` on the GPU (the first
version looped over families in Python — 12 M iterations — and could not finish).
Sets: `assoc` = vi2010 cascadia thr 0.3 N>=4 (3,364), `null` = the step-1 chance events at the same
cut (1,134), `bostock` = 4,000 of Bostock's own events at the step-5 per-family origin times,
`random` = 4,000 uniform times. Only 362 of the 129x6 (family, station) templates exist, so the
effective network is **2.8 stations per family**, not 6.

Per-TEMPLATE random-time background (each family's own): median cc 0.1059, MAD 0.0162.

| set | best-family cc median | frac above the null 99th pct | best-family MAD median | frac >= 6 MAD |
|---|---|---|---|---|
| random | 0.2398 | 0.0175 | 4.61 | 0.144 |
| null (chance events) | 0.2378 | 0.0106 (1 % by construction) | 4.50 | 0.123 |
| **assoc** | 0.2483 | **0.0366** | **4.92** | **0.233** |
| bostock (best of 129) | 0.2413 | 0.0292 | 4.73 | 0.198 |
| **bostock vs its OWN template** | **0.1164** | 0.0020 | **0.60** (p90 3.66) | **0.023** |

Two findings, and they point opposite ways.
1. **There is a real but small excess in the associated events**: 3.7 % of them exceed the
   chance-event 99th percentile against the 1.1 % the null gives by construction and the 1.8 %
   random times give, and 23.3 % reach 6 MAD against 14.4 % for random times. So roughly one
   associated event in ten carries a family-template match above background. A hybrid
   DL-detect -> template-classify product is **weakly** viable, not strongly.
2. **Bostock's own catalogued events register at only 0.60 MAD with their own family templates on
   this subset** (p90 3.66 MAD; 2.3 % above 6 MAD). The audit's "Bostock events at >= 6-8 MAD ->
   D3's 3.5 MAD was a template-quality artifact" branch does NOT fire; **D3 stands**. Scaling our
   ~2.8-station (8-9 channel) network to Bostock's 21 stations by sqrt(channels) gives ~1.0-1.5 MAD,
   still far under any detection threshold — consistent with D3's arithmetic and with Muenchmeyer's
   own Figure S10, where the correlation values of their uncataloged Cascadia detections against the
   family stacks were "identical to the noise distribution".
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

BUGFIX logged: the first version normalised everything by the MAD of the RANDOM-time distribution of
the BEST-of-129 cc, which is a maximum over 129 draws and far above any single template's
background; Bostock's own-template cc then came out at -5.7 MAD, which is meaningless. The full
(n_events x n_fam) matrix is now saved (`e3_classifier_perfam.npz`) and every cc is expressed in the
MAD of its own template's background.
Files `e3_classifier_{events.parquet,summary.csv,mad.csv,perfam.npz}`, figure `figures/e3_classifier.png`.

### STEP 5 — per-family Bostock time base (`e3_bostock_timebase.py`)

Peak of the (PGC S pick - Bostock DB time - ttS_PGC) residual distribution per family, for the 100
families with >= 50 events in the pgc2010 window, at cascadia thr 0.1.
**Result: 90 of 100 families accepted; median offset +11.75 s, sd 4.14 s, IQR +9.25..+12.75, full
range +2.50..+20.50 s.** So the offsets DO vary by more than 5 s across families — but the audit's
failure branch ("vary > 5 s with no pattern -> the convention needs the authors; keep Lin-hq
primary") only half fires, because the variation is **real, reproducible and family-specific**:

* **It is a SOURCE property, not a station property.** Decomposing 839 accepted (family, station)
  offsets over 100 families x 10 vi2010 stations: **between-FAMILY sum of squares / total = 0.565,
  between-STATION = 0.009.** Per-station offsets are nearly identical (PGC +12.50, SNB +12.50,
  VGZ +11.75, LZB +12.00, PFB +11.25, YOUB +11.75, NLLB +12.50, B011 +12.50, B926 +11.50,
  LRIV +13.00) — F7's apparent station dependence was the residual of a family mixture.
* **Cross-validated.** The offset measured at PGC over 2010-05..11 versus the median offset at
  >= 3 OTHER stations over the independent 2010-08..09 window: **r = 0.959 (p = 8e-50, n = 90),
  PGC-minus-others median -0.38 s, sd 1.19 s, 81 % within 1 s, 97 % within 2 s.** The per-family
  table is therefore good to ~1 s, against the +-5 s smear the single constant caused.
* **The convention itself is still not identified.** offset = origin - DB time, so the DB time sits
  ~12 s BEFORE the origin. No single-station reference-arrival model explains the scatter:
  offset + ttP_LZB has sd 3.80 s and offset + ttP_PGC 3.85 s, versus 4.14 s for the raw offset —
  no reduction worth anything. Weak correlations only (Spearman -0.39 with epi_PGC, -0.43 with
  median magnitude, +0.33 with latitude). **The convention needs the paper or the DB maintainers;
  the empirical per-family table is the working fix and Lin-hq stays primary.**

**Re-scoring with the per-family offsets.** `e2_assoc_score.py` and `e2_loc_compare.py` both gained
`--bostock-family-offsets` (default = the old single constant, so E2 behaviour is untouched) and a
`--tag-suffix` so E3 outputs cannot overwrite E2 files. On vi2010 cascadia thr 0.3: the Bostock
time-matched horizontal offset improves **25.3 -> 23.7 km** and the depth difference -24.1 -> -22.2
km at the same match count (967 vs 976 within 2 s). Modest, as expected — the per-family correction
is ~1-4 s while the location error is 25 km.
**Bostock-referenced precision/recall, three ways** because origin-time matching is unfair at 25 km
mislocation (F4).  Re-scored files `e2_assoc_score_vi2010_cascadia_t0.3_n{3,4,5}_e3famoff.csv`
(the E2 originals untouched), vi2010 cascadia thr 0.3, per-family offsets:

| N>= | recall +-2 s (null) | precision +-2 s (null) | recall +-5 s | precision +-5 s | prec. peak-excess +-2 s | peak/bg |
|---|---|---|---|---|---|---|
| 3 | 0.077 (0.024) | 0.107 (0.037) | 0.157 | 0.197 | **0.092** (E2 const: 0.069) | 3.2 |
| 4 | 0.041 (0.010) | 0.158 (0.042) | 0.082 | 0.278 | **0.158** (E2 const: 0.129) | 4.2 |
| 5 | 0.019 (0.003) | 0.230 (0.041) | 0.034 | 0.383 | **0.307** (E2 const: 0.218) | 6.8 |

The per-family offsets sharpen the residual peak (peak/background 3.2/4.2/6.8) and raise the
chance-corrected precision by 22-41 % relative at every station cut — the single constant was
smearing real matches out of the window, exactly as F7 predicted.  Third way, by
**ARRIVAL time** — which is what config D is, and which is the honest version: 124 Bostock events at
thr 0.1 have >= 3 of the 6 overlap stations carrying a real S pick within +-1 s of their predicted
arrival, out of 15,931 in the window, i.e. an arrival-matched network recall of **0.8 % at N>=3**
(19 events, 0.12 %, at thr 0.3) — with per-station match rates of 9.2 % (PGC), 11.0 % (LZB),
12.2 % (VGZ), 19.1 % (YOUB), 17.3 % (SNB), 3.4 % (PFB) at thr 0.1.

**The chance-dominated C estimator is gone.** `e2_groundtruth.py` no longer computes
`C = median(Lin OT - Bostock t)` over Bostock events with a Lin event within 20 s (16,026 of 18,114
have one BY CHANCE at ~500 Lin events/day, so the median is pulled to ~0 whatever the truth is; it
reported +0.086 s). The code block is replaced by a comment explaining why.  The CSV it wrote is LEFT IN PLACE
untouched (E3 modifies no E2 output); a sibling
`results/e2_bostock_lin_dt.WARNING_chance_dominated.txt` now sits next to it saying what it is not
and pointing at the per-family table. `e2_assoc_score.bostock_offset()` now treats a missing file as C = 0.
Files `e3_bostock_{family_offsets,famsta_offsets,offset_pattern}.csv`, figure
`figures/e3_family_offsets.png`.

### STEP 0 — Muenchmeyer et al. (2024) numbers -> `results/e3_munchmeyer_numbers.md`

arXiv:2311.13971v1 fetched (3.1 MB, 15 pp). **The supplement is not in the arXiv package**
(`/src/2311.13971/anc` -> 404), and the PyOcto parameters, pick thresholds and event filters are all
deferred to it — so E2's `n_picks=3 / n_s_picks=3 / tol 2.5 s / edt_std 1.5 s` cannot be checked
against theirs. Recorded as an open item needing the published version, not another guess.
* **Their Cascadia DL catalogue is 8,369 events over 59 days (2003-02-26..03-10, 2004-07-02..07-27,
  2005-09-03..09-25) = 141.8/day. The Bostock reference over the identical windows, computed here
  from the catalogue on disk, is 148,650 events = 2,519.5/day. Ratio 0.056 — their catalogue is
  ~18x SMALLER than the reference.** Our vi2010 N>=4 thr-0.3 run gives 100 events/day against 2,663
  Bostock/day in its window: **we are at the same place on the same curve as the picker's authors.**
* Their recall/precision definitions: one-way match fractions at **+-10 s**, no chance correction.
  "For 60 % of the LFEs from our workflow we find an LFE in the reference catalog within 10 s" is
  precision-like; the "32 % of the cataloged events" sentence follows the Nankai discussion and
  cannot be Cascadia's (32 % of 148,650 > their whole catalogue). Implied Cascadia recall <= 3.4 %.
  **We measured the chance floor of their criterion**: 2x10^5 random times per window give a Bostock
  event within +-10 s **28.3 % of the time** (0.310/0.249/0.289), so their 60 % is 2.1x chance and
  only 53 % of their "matches" are informative. At +-2 s the same floor is 8.8 %.
* **They DO report a chance-association control** (sec. 4): per-station-per-hour constant shifts,
  giving scrambled catalogues at **5-10 %** of the original, which they call an over-estimate of the
  false-positive rate. That is the direct precedent for step 1, it agrees with the audit's
  <=15 % PASS bar, and our thr-0.5 N>=4 / thr-0.3 N>=5 configurations (10 % / 8.8 %) sit exactly
  there while the E2 operating point (33-59 %) does not.
* Also relevant to step 4: their Figure S10 found the correlation of their uncataloged Cascadia
  detections against the family stacks to be "identical to the noise distribution".

### STEP 6 — documentation corrections

**(a) Off-ground-truth attribution (F6), recomputed independently — Merlin's numbers reproduce exactly.**
`e3_offgt_fix.py`. PNSN tremor DAYS in 2017-2018 within 40/60/80/100 km:
**WISH 26/62/104/134, YACT 19/44/72/94, B943 109/146/173/210** (and PGC 92/133/159/185, VGZ
126/147/181/193, LRIV 103/158/189/205). WISH and YACT simply have almost no tremor near them; they
cannot test generalisation.
Inside the 2017 Vancouver Island source zone, SURFACE stations respond (cascadia thr 0.3,
ETS/quiet from the full-record 40 km rolling statistic): **SYMB 2.26, PGC 2.25, VGZ 2.02, CLRS 1.99,
SNB 1.81, LRIV 1.34 — against the co-located BOREHOLE B011 at 1.78** (Spearman +0.51..+0.77 for the
surface stations, +0.69 for B011). The borehole is not categorically better; it is not even the best.
Direct surface/borehole control at 0.6 km separation, single-station +-1 s peak-excess recall vs
Lin-hq in vi2010: **PGC 0.102 vs B011 0.125 at thr 0.1 (1.23x), 0.033 vs 0.059 at thr 0.3 (1.77x)** —
and the best station of all is the SURFACE station **SNB at 0.167 / 0.077**, ahead of both
boreholes. **The E2/D3 reading "surface flat, borehole responds" is distance to the source, not
instrument type or generalisation, and must not seed a borehole-only fleet design.**
**(b) Lin attribution (F2).** Lin's 90 %/95 % are held-out balanced-test-window numbers with an
equal noise class, not continuous-record catalogue numbers; his own published catalogue is a 15 s
N>=3 coincidence catalogue of which his quality cut keeps 9.8 %. E2's "FAIL vs Lin's own 90 %
recall / 95 % precision" compares two different quantities and is withdrawn as a comparison.
**(c) F1 readiness at 40 km (F8).** `e3_f1_readiness40.py` recomputes F1's own cell/station tables
at 40 km (the measured usable radius: single-station recall 16 % at <25 km, 11 % at 25-40 km, 3.2 %
at 40-60 km). Method validated first by reproducing F1's stored n60 (exact in 75-85 % of cell-years,
within +-1 in 92-95 %). Pooled over all 137 cells, project fleet:
**frac >=4 3-C stations: 0.34 at 40 km vs 0.92 at 60 km in 2018; 0.78 vs 1.00 in 2024.
frac >=6: 0.08 vs 0.56 in 2018; 0.34 vs 0.985 in 2024.** F1's own report was not modified;
`results/f1_readiness_40km.csv` is the new table.

### Bugs found and fixed in this run (all minimal, all commented at the fix site)
1. **Pick-threshold selection used `thr_first <= t` without `prob >= t`** in three new E3 scripts.
   A pick exists at threshold t iff `thr_first <= t <= prob` (the `e2_pick` docstring), so the wrong
   filter makes EVERY threshold identical to the 0.1 pick set. It was caught because the real and
   injected recall curves came out bit-identical at thr 0.1 and 0.3, and because
   `net_picks_per_day` came out *increasing* with threshold (11605 -> 11609 instead of
   4895 -> 1088). Fixed in `e3_inject.py`, `e3_denoise.py`, `e3_inject_vs_real.py` and both runs
   redone. **The same slip is present in the audit's own config-D counts (698/169/27).** No E2
   result is affected — `e2_associate.load_picks` has always used both conditions.
2. **Step 5's first per-family estimator accepted noise.** 0.25 s bins over +-40 s with only
   ~100-300 pairs per family leaves the background bin count at 0 for most families, so
   `peak / max(bg, 1e-9)` returned 3e9 and a 3-count bin passed a "peak/bg >= 3" test: it reported
   96 of 100 families with offsets scattered over 59 s. Replaced by a 2 s sliding window stepped
   0.25 s with a Poisson z against the flat background outside the peak and a **4 sigma** bar
   (4, not 3, because ~320 window positions are scanned per family). 90 of 100 accepted, range 18 s.
3. **Step 4's MAD scale** (see step 4): normalised by the best-of-129 background instead of the
   per-template background; gave -5.7 MAD for Bostock's own events.
4. **`gt` as a column name shadows `DataFrame.gt`** — fired again in `e3_figures.py`
   (`w.gt.unique()` -> `'function' object has no attribute 'unique'`). Bracket access everywhere.
5. **pandas reads the string `"null"` as NaN**: the event set literally named `null` came back as
   NaN from `e3_classifier_summary.csv` and appeared as a bar labelled "nan". Read with
   `keep_default_na=False` and coerce the numeric columns.
6. Backgrounded jobs launched with a plain `&` are killed when the Bash tool times out at 2 min
   (exit 143); `setsid nohup ... < /dev/null & disown` is required. One denoise run was lost to this.
7. `d3_vel1d.tt` takes `(x_km, z_src, phase)`, not `(lat, lon, dep, slat, slon)`.
8. **`e2_assoc_score.py` appends every score row to the master `e2_assoc_scores.csv`**, so the E3
   re-score put 24 `_e3famoff` rows into an E2 output. Caught in the final audit; the 24 rows were
   split out into `e3_assoc_scores_famoff.csv` and the E2 master restored to its original 144 rows
   (18 configurations), and the script now writes to the E3 master whenever an E3 option is given.
   This was the only E2 file E3 modified; everything else with an `e2_` prefix written today is a
   NEW file carrying the `_e3famoff` suffix.

### Decisions a human should make
1. **The operating point must move.** thr 0.2 / N>=4 is chance-dominated (0.59) and thr 0.3 / N>=4
   is one-third chance. thr 0.5 / N>=4 (15 ev/day, 10 % chance) or thr 0.3 / N>=5 (28 ev/day, 8.8 %)
   are the defensible ones, and they match Muenchmeyer's own 5-10 % scrambled-catalogue standard.
2. **No fine-tune, and no fleet download, on this evidence.** Step 3 shows the injection test cannot
   establish an SNR ceiling and that a fine-tune would have to be judged on real single-event
   positives; step 3b shows denoising is not the lever; step 1 shows the required F reduction is
   only ~1.7-2.3x and is already obtainable by moving the threshold.
3. **Step 6 of the prescription (network back-projection of the S-probability traces) is now
   triggered** — step 1 put the recommended operating point in the chance-dominated/unresolved band.
   It was NOT run here (it was conditional and outside steps 0-5 + 6-documentation).
4. **The Bostock DB time convention still needs the authors** for anything Bostock-referenced to be
   quotable; the per-family table is a measured workaround, good to ~1 s.

## 2026-09-11 ~02:15 Merlin re-audit of E3 -> results/merlin_audit_e3_2026-09-11.md (verbatim)
Step 1 (chance null) and the E2 retractions STAND; ML route recovers ~1-3% of Bostock's events -> cannot serve the GOAL alone (role: family
seed). Step 4 FALLS: the SeisBench stacks are stored Z,1,2 with 1=N, 2=E (verified by P polarization vs backazimuth, 9° vs 40°), and
e3_template_classifier.py / e3_inject.py / e3_denoise.py ordered data Z,E,N via sorted() -> templates and injected LFEs were 90°-rotated.
"Bostock's own events at 0.60 MAD" and "D3 stands" are UNTESTED; the step-3 "paradox" is mostly the coherent-amplitude s/n anchor plus
the swap; config D partly circular. Also untested cheap lever: association tolerances 5x the measured residual. On-disk fact: Bostock's
catalog covers the Sept-2012 SSE (28,149 events, 130 families) -> the Thomas harness can be validated with zero download. Prescription
E4: 0 fix components; 1 CORRECTED template test at Bostock times with empirical per-station moveout, 2/4/8-s templates (DECIDING:
median >=3 MAD at M>=1.8 and >=30% recall at 8 MAD with <=10% chance = template route open); 2 corrected injection + incoherent s/n
anchor; 3 tight-tolerance association re-null; 4 Thomas-2025 harness on Bostock's Sept-2012 catalog; 5 metadata availability for the
21 template stations 2012; 6 (conditional on 1) targeted Sept-2012 VI download 10-20 GB; 7 back-projection only if 1 fails.
User go 2026-09-11 ~02:50 -> E4 agent launched (steps 0-6); Merlin vetting the post-2019 redesign (T1-T4) in parallel.

## 2026-09-11 ~03:05 Merlin vet of the T1-T4 redesign -> results/merlin_vet_e4_redesign_2026-09-11.md; amendments sent to the E4 agent
Binary bar replaced by a sensitivity model + completeness-requirement curve + prediction table (see PLAN.md E4). Hold-out -> alternate 2-h
blocks; transfer stations B011/B926/LRIV(+STW); SOKB 2010 empty; T2a q-control added; decimation leg added to the Thomas harness; T3/T4
re-specified (2010 data; drift controls) and deferred to a separate spec.

## 2026-09-11 E4 — component fix, the corrected template test, transfer, tight-tolerance null, Thomas harness

Environment: same as E3 (torch 2.12.0+cu130 in `~/.local`, NVIDIA L40S, seisbench 0.12.5, pyocto
0.2.0) under the exact env exports of `scripts/e2_run_pick.sh`. Nothing installed. Concurrency at
or below the 28-worker cap (one GPU process at 4 torch threads; PyOcto 3 x 7 threads).
**Fetches: the PMC page of Thomas et al. 2025 (text) and FDSN station METADATA. No waveforms.**
No E2/E3 output deleted or overwritten; every new file is `e4_`-prefixed, and the three E3 scripts
that were fixed keep the old behaviour behind `E3_LEGACY_CHAN_ORDER=1`.

Merlin's E4 vet (`results/merlin_vet_e4_redesign_2026-09-11.md`) arrived mid-run and was applied
before step 1 was interpreted and before step 2 was started: the step-1 binary bar was replaced by a
sensitivity model, template lengths became 4/6/8/12 s (2 s kept as a diagnostic), the hold-out became
alternate 2-h blocks, the transfer set became B011/B926/LRIV(+STW), and the Thomas harness gained a
decimation leg. The report records both the original acceptance bars and the model.

### STEP 0 — the component order. CONFIRMED, FIXED, AND E3 STEP 4 IS RETRACTED.

`scripts/e4_polarization_check.py` (new) reproduces Merlin's test independently: 327 SeisBench
stacks with P snr > 6 dB and 3-C rectilinearity > 0.7; the azimuth of the horizontal P principal
axis against `path_backazimuth_deg` (mod 180) has median **7.9 deg under slot1 = N (90 % within
25 deg, 95 % within 45 deg)** against **37.3 deg under slot1 = E (36 % / 58 %)**. So
`component_order = "Z12"` means **1 = N, 2 = E**, exactly as the audit said.
`results/e4_polarization_check.csv`.

Fixed by an explicit resolver, `scripts/e4_chanorder.py` (new): channels are selected by the LAST
LETTER of the code -- Z, then N or 1, then E or 2 -- and two different horizontal conventions on one
station-day raise and stop. Call sites, all commented at the fix with the old line quoted:
`e3_template_classifier.py:190`, `e3_inject.py:132`, `e3_denoise.py:194` and `:317`.
`E3_LEGACY_CHAN_ORDER=1` reproduces the pre-fix behaviour.
**Sweep of every station-day on disk** (`e4_chanorder` over `traces/*/*/`): no station has both N
and 1 -- CN/UW surface stations are `?HZ,?HN,?HE`, PB boreholes are `EHZ,EH1,EH2`. Resolved orders
printed per station-day by every script now. Two data facts found in the same sweep and used below:
* **`traces/C8.SOKB/2010/` is 34 ZERO-BYTE `.mseed.nd` placeholders -- SOKB has NO data on disk.**
  The 2010 outer ring is therefore SIX stations (PGC, LZB, VGZ, SNB, YOUB, PFB), not seven.
* `traces/CN.GOBB/2010`, `UW.OCP`, `UW.STW` are **EHZ only** (single component); `CN.GOBB/2017` is
  also all-zero `.nd`.

Documentation corrections made: RESULTS.md step 4 is marked **RETRACTED** with the polarization
numbers; the clause "and it agrees with Muenchmeyer's own Figure S10" is withdrawn on its own merits
(F2: S10 is about UNCATALOGED DL detections, not catalogued events, which exceeded Bostock's own
network-cc threshold by construction); steps 3 and 3b are flagged swap-affected with the note that
3b's FIDELITY finding is orientation-independent and stands. The buggy `mad_median` column is left
in `e3_classifier_summary.csv` (E4 overwrites no E3 file) with a sibling
`e3_classifier_summary.WARNING_mad_median.txt` and a relabelled copy
`e4_e3_classifier_summary_relabelled.csv`.

### STEP 1 — the corrected template test. `scripts/e4_template_test.py` (new), GPU.

Method, in one paragraph. Full-day continuous normalised cross-correlation at 20 Hz, 1-8 Hz, per
channel, computed by rFFT with the DATA transform computed once per station-day and reused by every
family and every template length. Templates are the family's own SeisBench stack, S-0.5 s to
S-0.5+L, each component L2-normalised separately, so a station's statistic is the mean of 3
per-channel normalised cc and the network statistic is the channel-count-weighted mean over stations
(Merlin F4: threshold the network SUM normalised by the contributing channels, never per station).
Placement: origin = Bostock DB time + the E3 step-5 per-family offset; template start at station s =
origin + Crosson ttP(family hypocentre, s) + (stack S sample - 800)/20 - 0.5 s. Background for every
(family, station, day) and every (family, day) network series = the median and MAD of THAT
template's own continuous series that day, gaps masked. Event score = max over an INDEPENDENT
+-3 s scan per station (its argmax is the per-event best lag) and over a COMMON +-1 s scan for the
network. Pass 1 measures the empirical moveout (median best lag per family-station-length); pass 2
applies it. Controls: N random times per day scored identically; declustered peak counts above
6/8 MAD in the day's own series = false alarms per template per day; and a per-station TIME-SHIFTED
network (each station's cc series rolled by its own 300-3000 s) = the chance level of the network.

**Two bugs found and fixed before any number was read.**
1. **float32 catastrophic cancellation in the sliding-energy prefix sum.** The running sum of
   squares for a station with counts ~1e4 reaches ~3e14, whose float32 ulp is ~2e7, so differencing
   two prefix sums for a quiet 2-s window came out ~0 and cc = num/1e-10 ~ 1e12. It showed up as
   `mad_scale` 7e11 at PGC, 4e11 at SNB, 1e10 at YOUB and 6e9 at PFB while LZB and VGZ (1000x
   quieter instruments) were fine at 0.06. Fixed by doing the prefix sum in float64 and flooring the
   window energy at 1e-10 x its median (`e4_template_test.py`, the `EN[s]` block and the `e` block).
   After the fix every station's `mad_scale` is 0.044-0.090, against the white-noise prediction
   0.6745/sqrt(2 B L)/sqrt(3) = 0.052 at L = 4 s.
2. **Gap masking.** PFB's 2010 day files merge from 231 fragments; zero-filled gaps give an
   undefined cc. Samples whose 3-C window energy is below 1e-6 of the day median are now excluded
   from every median, MAD, extremum and false-alarm count.
Also fixed: the network was being assembled on each station's own day-file clock; it is now built on
a COMMON origin axis with per-station integer shifts (`k[s] = round((tref - T0[s] + o[s]) * FS)`),
which matters whenever two stations' day files do not start at the same instant.

### STEP 4 — association re-null at tolerances matched to the measured residual.

`e2_associate.py` gained `--zlim` (default = the old area, so E2/E3 behaviour is untouched);
driver `scripts/e4_run_assoc_null.sh`, scorer `scripts/e4_assoc_null_score.py` (both new).
32 PyOcto runs at **pick_match_tolerance 1.2 s, edt_pick_std 0.5 s, zlim 20-50 km**, with
`n_picks` run DIRECTLY at 4 and 5 (E3 ran n_picks = 3 and post-filtered on n_stations), thresholds
0.3 and 0.5, 5 constant-shift seeds and 2 hourly-shift seeds each.

| thr | N | obs/day | null/day (5 seeds) | null/obs | E3 wide-tolerance null/obs |
|---|---|---|---|---|---|
| 0.3 | 4 | 71.5 | 12.54 +- 0.44 | **0.175** (hourly 0.200) | 0.330 |
| 0.3 | 5 | 19.0 | 0.494 +- 0.084 | **0.026** (hourly 0.038) | 0.088 |
| 0.5 | 4 | 11.3 | 0.518 +- 0.061 | **0.046** (hourly 0.052) | 0.103 |
| 0.5 | 5 | 2.41 | 0.029 | **0.012** (hourly 0.024) | 0.017 |

Tightening the tolerances to the measured 0.40-0.56 s S residual and restricting depth to the
interface band **halves to quarters the chance fraction at every operating point** (0.33 -> 0.175 at
thr 0.3/N>=4; 0.088 -> 0.026 at thr 0.3/N>=5) while keeping 74 % of the events (96 -> 71.5/day).
The prescription's acceptance ("null/obs <= 0.10 at N >= 4 thr 0.3 with obs >= 60/day") is **not**
met -- 0.175 at 71.5/day -- but thr 0.3/N>=5 now sits at 0.026 with 19 events/day and thr 0.5/N>=4
at 0.046 with 11.3/day, i.e. **cleaner than Muenchmeyer's own 5-10 % scrambled-catalogue standard.**
Bostock recall (per-family offsets, +-2 s, chance-corrected) at thr 0.3/N>=4: 0.70 % (M<1.5),
2.9 %, 4.5 %, **9.4 % (M>=2.1)**, 3.4 % overall. At thr 0.5/N>=5: 0.2 % overall.
**The DL catalogue's role as a family seed, not a Thomas-grade catalogue, is unchanged.**
EQ/blast cross-match: `data/eq_catalog_pgc.csv` and `data/eq_catalog_b011.csv` (3,636 / 4,240
events, 2007-2017) hold **44 regular earthquakes inside the vi2010 window**; 0.45-1.2 % of our
associated events fall within 60 s of one (chance 0.18 %), i.e. earthquake contamination is
present but negligible. `results/e4_assoc_{null_table,bostock_recall,eq_crossmatch}.csv`.

### STEP 5 — the Thomas-2025 harness, `scripts/e4_thomas_harness.py` (new, catalog-agnostic).

The PMC full text was fetched. **Thomas's front and tidal analysis is QUALITATIVE** -- the secondary
fronts are "vertically oriented groups of events closely spaced in time" read off a space-time plot;
the paper gives no detection algorithm, no binning, no completeness cut, no rate-decay functional
form and no tidal method. Its stated numbers: main front SOUTHWARD at "~5 km/day" (also "~10
km/day"), ~80 km in 11 days; secondary fronts initiating 1-10 km from the main front, 10-40x faster
(example ~12 km/h), lasting 2-6 h, median stress drop 8 kPa; four marked in Fig. 2A and family 141
hosting six between 12 and 16 September; early LFEs "not significantly influenced by the tides"
while later secondary fronts are "almost always tidally triggered". **The operational criteria are
therefore OURS**, are printed with every run and are written to
`results/e4_thomas_criteria_bostock2012.json`.

Two harness bugs found and fixed before the numbers were read:
1. **Wrong along-strike axis.** `data/cas_slab2_input_04-18.csv` is the scattered INPUT to Slab2
   (earthquake and receiver-function observations), not the fitted surface; its local gradient is
   noise and gave an along-strike azimuth of **+86 deg (E-W)** at southern Vancouver Island. Using
   the gridded `cascadia_mesh.txt` through `sse_strain.megathrust.slab_depth_interpolator` gives
   interface depth 38.4 km, down-dip azimuth **49.7 deg** and along-strike azimuth **-40.3 deg**
   (N320E) -- the correct Cascadia geometry, and it flips the sign of the front speed to
   SOUTH-EASTWARD, i.e. Thomas's "southward". Family along-strike span 112.5 km, down-dip 62.1 km.
2. **A broken secondary-front null.** The first null redrew every burst start uniformly over the
   episode, which SPREADS each family's clustered bursts and therefore INFLATES the number of
   distinct families per 6-h window; it produced a null ABOVE the observed count. Replaced by a
   CIRCULAR per-family time shift (each family's whole burst sequence rolled by its own random
   amount inside the episode), which preserves within-family clustering, burst counts and the
   spatial pattern and destroys only the inter-family timing.

### STEP 1 results — the numbers, and how they were read

Two passes over 34 days x 5 template lengths, ~2.5 min each on the L40S (the day's data FFT is
computed once and reused by 362 templates x 5 lengths). Pass 1 (no moveout) -> pass 2 (empirical).
`results/e4_tt_{events,station,netfam}_ring_p{1,2}.parquet/.csv`, `e4_moveout_ring_p1.csv`,
analysis `scripts/e4_step1_model.py` -> `e4_step1_{registration,network,channels,magslope}.csv`,
figure `figures/e4_step1_registration.png`.

Empirical moveout: median best single-station lag **+0.41 s** (IQR +0.08..+0.64, s.d. 0.39) over
1,325 (family, station, length) cells -- Merlin's F3 prediction of a 0.3-0.5 s station-to-station
Crosson error, measured. Applying it lifts the 8-s network excess from +2.41 to +3.34 MAD and the
8-MAD recall from 0.140 to 0.223.

Headline (pass 2, M1.8-2.1 unless stated): network median 4.41 / 5.30 / 5.65 / 5.89 / 6.15 MAD at
L = 2 / 4 / 6 / 8 / 12 s on a random-time background of 2.73, i.e. an excess of +1.66 / +2.57 /
+2.92 / +3.17 / +3.42; 8-MAD recall 0.037 / 0.117 / 0.168 / 0.206 / 0.256; the per-station
time-shifted network fires on 0.0000-0.0002 of the same events; DETECTIONS (renamed E5 step 0 -- these
are total declustered detections, NOT false alarms) 0 / 1 / 1 / 2 / 4 per template per day at 8 MAD,
median over family-days (33 / 58 / 72 / 82 / 102 at 6 MAD).
Per-station excess / prediction (q r1 / MAD_scale, 8 s, M1.8-2.1): LZB 2.11, PGC 2.16, YOUB 1.96,
SNB 1.75, PFB 1.58, VGZ 1.37 -- median 1.9x, i.e. Merlin's "the stack-SNR convention understates r1".
Single-channel (Z, against its OWN randZ background) 0.13-0.46 MAD against 0.79-1.57 for the 3-C
mean: gain 3.4-6.1x, not sqrt(3). Within-family magnitude slope 1.28 MAD per unit on a base of 1.45
-> amplitude scaling 10^(0.38 dM); between-family slope NEGATIVE (distance confound).

### STEP 2 — `scripts/e4_transfer.py` (build) + `e4_run_step2.sh` / `e4_step2_model.py` (evaluate)

Hold-out = alternate 2-h blocks, 60-s guard (build 7,849 / test 7,771 events, 217 guarded).
82 families with >= 40 build-half events; median 90 events per stack; 820 templates over 10 stations
in 174 s. Three evaluation legs, each measuring its moveout on the BUILD half and scored on the
TEST half: A = ring with the published stacks, B = ring with our 2010 re-stacks (published ones
dropped), C = ring published + transferred B011/B926/LRIV/STW.
* q-control: cc(our re-stack, SeisBench stack) PGC 0.785, LZB 0.832, VGZ 0.729, SNB 0.842,
  YOUB 0.861, PFB 0.776 -- all above the 0.6 bar at N = 90, far above the 0.3 stop-line.
  No family reaches N >= 200 in the build half (max 167); at N >= 150 (32 pairs) median cc 0.837,
  84 % above 0.6, and the five pairs below 0.3 are ALL family 15 at five different stations
  (stack S-SNR 1.0-5.1 dB, S pick 0.9-2.2 s off) -- one bad family, not a global sign error.
  Held-out registration ratio (ours / published, 8 s): 0.69 / 0.90 / 0.45 / 0.56 / 0.83 / 0.75
  against a predicted q-ratio of 0.53 / 0.57 / 0.46 / 0.43 / 0.53 / 0.34 -- measured >= predicted at
  5 of 6 stations.
* transfer: stack S-SNR B011 8.60 dB, B926 9.66, LRIV 2.76, STW 0.06 against same-N random-time null
  stacks at -0.03 / 0.19 / 0.19 / -0.08 dB (ring reference PGC 11.21, LZB 11.78, YOUB 10.20,
  SNB 7.46, PFB 6.23, VGZ 5.57). Held-out single-station excess (8 s, M1.8-2.1) B011 1.188,
  B926 1.873, LRIV 0.388, STW 0.042 against the co-located PGC 2010-only re-stack 1.404 -> B011
  0.85x (PASSES Merlin's 0.7x bar), B926 1.33x, LRIV 0.28x, STW 0.03x.
* network: A 2.910, B 2.808, C 2.636 MAD (8 s, all magnitudes) -- **adding four stations made it
  WORSE**. Diagnosed, not hand-waved: an equal-weight channel mean gives e_net = (sum e_s)/sqrt(n),
  so a station helps only if its excess exceeds ~half the running mean; LRIV and STW do not. The
  SNR-weighted optimum sqrt(sum e_s^2) rises 3.70 -> 4.31 (+17 %) when they are added.
  Measured / equal-weight prediction 0.86 (A) and 0.82 (C), both inside Merlin's 30 % band.
  **Action for the fleet design: SNR-weighted network stacking, not an equal-weight mean.**

### STEP 3 — `scripts/e4_inject.py`

Incoherent s/n measured with a control (4,000 random times on the same days: power ratio 1.077;
Bostock times shifted 300-3000 s: 1.102). Bostock times: 1.40-3.20, i.e. s/n 0.57-1.46 by magnitude.
Isolated events (<= 1 neighbour within +-30 s) 0.59-0.97 vs dense 0.82 -> NOT burst contamination.
Injection: 3 quiet PGC days x 140 injections x 3 template kinds (published stack, own all-event
stack, own N=20 stack) at snr_rms 0.05-2.0, each day also picked un-injected as the baseline.
Recall at thr 0.1: published stack 0.13 / 0.45 / 0.90 / 0.98 at s/n 0.2 / 0.3 / 0.5 / 1.0; own-all
0.12 / 0.35 / 0.67 / 0.97; own-N20 0.07 / 0.12 / 0.45 / 0.82. Real recall 0.11-0.26 at the measured
s/n 0.57-1.02 -> the injection is 4-6x ABOVE the re-anchored real curve; the curves cross at
s/n ~ 0.2-0.3. Acceptance ("within 2x") FAILS, so the fine-tune option is REOPENED (real
template-aligned positives, family-held-out by source_family_id) and was NOT run, per the brief.

### STEP 6 — availability, `scripts/e4_station_availability.py` + `e4_avail_summary.py`

Trap found and fixed: station CODES collide across networks. An unfiltered FDSN query reports TWGB
and TWKB as open to 2030 because Taiwan's TW network has stations of the same names (lat 22 N);
likewise LCBC in Colombia's CM, MGCB in Chile's ZR, OCP in CI and UU, B011 in SN and ZF. The
summary now filters on Cascadia network codes (PO, C8, CN, PB, UW, UO) and a 39-51 N / -128..-119
box. Results in RESULTS.md; the decisive one is that **all 14 inner POLARIS/C8 template stations
closed between 2004-06 and 2007-02, so a POLARIS 2010 pull is impossible** and Bostock's own 2010
and 2012 detections were made on essentially the seven-station geometry we are testing.
Volumes from 0.42 MB per channel-hour at 100 Hz measured on the files on disk.

### STEP 7 — `scripts/e4_prediction_table.py`

First version fitted a power law in distance to the 10 measured (distance, excess) points; it
predicted Mc 2.69 on Bostock's own array against his median M1.70 -- a failure by a whole magnitude
unit. The points show why: at 36-42 km the measured excess is 0.04 (STW, Z only), 0.68 (PFB),
1.15 (VGZ), 1.19 (B011), 2.15 (PGC). SITE dominates DISTANCE over this range and a power law
extrapolates badly. Replaced by a piecewise-empirical model (median measured excess, two distance
bands, by station class), which gives 7.0 MAD on Bostock's array -> Mc 2.10 at 8 MAD and Mc 1.46 at
4 MAD, and inverts to an implied Bostock operating threshold of **3.6 MAD -- D3's 3.5 MAD from a
completely different calculation.** Both threshold columns are reported because the threshold, not
the geometry, is what binds. [WITHDRAWN E5 step 0 (Merlin E4 audit F1-F3): excess compared to a
threshold that applies to the RAW statistic (raw = excess + 2.73); the 3.6-MAD inversion and the
4-MAD column are void -- see RESULTS.md E4 step 7 and the E5 section.]

### Bugs found and fixed in E4 (all minimal, all commented at the fix site)
1. **Horizontal components swapped** in `e3_template_classifier.py:190`, `e3_inject.py:132`,
   `e3_denoise.py:194` and `:317` -- `[Z] + sorted(non-Z)` = Z,E,N against the stacks' Z,N,E.
   Fixed by `scripts/e4_chanorder.py` (explicit last-letter resolution, stops on an N/1 clash);
   old behaviour behind `E3_LEGACY_CHAN_ORDER=1`.
2. **float32 catastrophic cancellation** in the sliding-energy prefix sum (`e4_template_test.py`):
   MAD scales of 7e11 at PGC, 4e11 at SNB, 1e10 at YOUB, 6e9 at PFB. Fixed with a float64 prefix sum
   and a relative energy floor; every station now lands at 0.044-0.090, against the white-noise
   prediction 0.052 at 4 s.
3. **Zero-filled gaps** gave undefined cc (PFB's 2010 days merge from 231 fragments). Samples whose
   3-C window energy is below 1e-6 of the day median are masked out of every statistic.
4. **Network assembled on per-station clocks** -- rebuilt on a common origin axis with per-station
   integer shifts.
5. **False alarms counted as rising edges** over-counted one excursion several times (the 1-8 Hz cc
   series oscillates); replaced by a max-pool decluster over one template length.
6. **Single-channel excess referenced to the 3-C random background** in the first channel-scaling
   table (the max-over-scan background of a 1-channel series differs); `randZ_median` is now
   recorded per (family, station, day) and used.
7. **Wrong along-strike axis** in the Thomas harness: `cas_slab2_input_04-18.csv` is scattered INPUT
   to Slab2, not the surface; its gradient gave +86 deg (E-W) instead of -40.3 deg (N320E).
8. **Broken secondary-front null** (uniform redraw of burst starts inflates the distinct-family
   count per window); replaced by a circular per-family time shift.
9. **Power-law distance extrapolation** in the prediction table (see step 7).
10. `shift` as a DataFrame column name shadows `DataFrame.shift` -- bracket access (the same class
    of bug as E3's `gt`).
11. `EP.load_model` takes `(name, device)`, not `(name)`.
12. `e3_inject.py` and friends use paths relative to `/home/jovyan/tremorferometry`, so anything
    importing them must be run from there.

### Decisions a human should make
1. **The step-1 verdict is "template route OPEN, threshold-limited".** The MAD and chance halves of
   the acceptance pass with room; the 30 % recall half does not (12-26 %). Whether that is enough
   depends on the threshold the fleet can afford, which is step 7's binding constraint.
   [E5 step 0: the step-7 threshold argument is WITHDRAWN (F1-F3); the operating point is ~8 MAD.]
2. **Two required design changes** before any fleet run: SNR-weighted network stacking (an
   equal-weight mean made the network worse when four stations were added) and 8-12 s templates.
3. **The 0.078 s/n anchor is retired.** The measured per-event incoherent s/n at PGC is 0.57-1.02
   and Bostock's magnitudes scale as 10^(0.4 dM), not 10^(1.5 dM). Anything conditioned on 0.078
   (D3's arithmetic, the E3 injection reading, parts of the fleet cost case) needs re-deriving.
4. ~~**Fine-tuning is REOPENED but not run** (the injection now sits 4-6x above the real curve).~~
   WITHDRAWN E5 step 0 (F7): coherent-vs-incoherent anchor conflation; fine-tuning CLOSED.
5. **Candidate downloads:** POLARIS Aug-Sep 2010 is IMPOSSIBLE (stations closed 2004-2007).
   Sept-2012 southern VI is 11 stations x 62 d = 20.1 GB native / 8.1 GB at 40 Hz.
   [SUPERSEDED E5: user's 15-station set (PLAN E5); volume recomputed in RESULTS E5.]
6. ~~**The secondary-front metric must be settled on Bostock's own catalogue first** -- our criteria
   find 0-3 against a null of 0-5.4 where Thomas found four by eye.~~ WITHDRAWN E5 step 0 (F4): the
   detector was gated behind a first-event front that never moved; redone in E5 step 4.

## 2026-09-11 ~04:45 Merlin re-audit of E4 -> results/merlin_audit_e4_2026-09-11.md (verbatim summary)
Template route OPEN stands (raw 5.9 MAD median at 8 s; 21% >=8; rolled control 0.02%). VOID: E4's "FA/template/day" (they are total
detections — quiet-day FA at 8 MAD ~0.6/template/day, at 6 MAD 45-115 -> operating point ~8 MAD), the "implied Bostock 3.6 MAD", the
"4 MAD suffices" conclusion, and the prediction table (compared EXCESS to threshold; raw = excess + 2.73 -> Mc ~0.5 mag too pessimistic:
corrected 8-MAD Mc Puget 1.71 / N-OR 1.73 / S-VI 1.94 -> ~50/50/17% of Bostock 2012 density; Bostock 21-station sanity reproduces his
modal M1.73). Thomas harness: main-front speed window-dependent (Thomas = 80 km/11 d ~7 km/day); secondary-front detector was disabled by
its own front model (bottleneck NOT established); tidal growth = clustering artifact. Fine-tuning CLOSED (coherent vs incoherent anchor
conflation; coherent-anchored injection matches real within 2x). 2012 S-VI network = ~17 sites (12-17 GB @40 Hz). DOWNLOAD APPROVED
conditional on E5 steps 1-4 + predicted >=20% recall. Awaiting user go.
## 2026-09-11 ~05:00 user: exclude C8 + LZB/YOUB/PFB from the 2012 set -> 15 stations (see PLAN E5). E5 steps 0-4 launched (no download).
## 2026-09-11 04:54 UTC: E5 Opus agent died at launch (HTTP 429, Opus session limit until 04:30 UTC next day); relaunched the identical spec on Fable (user: "Continue with rest"). No partial output existed.

## 2026-09-11 E5 steps 0-4 — FA budget, corrected prediction, SNR-weighted hold-out, Thomas-harness fixes (no download)

**Agent: Fable 5.1** (the Opus launch died on HTTP 429 at 04:54 UTC — Opus quota exhausted until 04:30 UTC next day; the
identical spec was re-run on Fable, user: "Continue with rest"). Environment as E4 (torch cu130 in ~/.local, L40S, one GPU
process at 4 torch threads, all CPU work < 28 threads). **Fetches: the PMC text of Thomas et al. 2025 (two prompts on the
cached page). No waveforms.** Nothing of E2-E4 deleted or overwritten; every new file is `e5_`-prefixed; the two E4 scripts
changed (`e4_template_test.py`, `e4_transfer.py` untouched; only the correlator) keep the old behaviour unless `--budget` /
`--weights` is given. Merlin's E4 audit (`results/merlin_audit_e4_2026-09-11.md`) is the specification; its flaw numbers are
cited as F1..F12 below. Station set per the user (PLAN E5): 2012 replication = PGC, VGZ, NLLB, SNB, LRIV, B001, B003-B007,
B009-B011, B926 (15); on the 2010 data the OPERATING SUBSET = PGC, VGZ, SNB + transfers B011, B926, LRIV, NLLB (7), the full
E4 ring + transfers (10) and the E4 6-station ring are the reference columns.

### STEP 0 — record corrected (RESULTS.md E4, NOTES E4 decisions, PLAN GOAL)
RESULTS E4: step-1 table column renamed "detections/template/day (median over family-days; NOT false alarms)"; struck with a
one-line reason + audit pointer: the 3.6-MAD inversion and the "8 MAD is a magnitude more conservative than Bostock" sentence
(F2/F3), the step-7 "crossed with the decimation leg / 4 MAD suffices / the threshold binds / 82 and 2 FA per template per
day" paragraph (F1), "main front REPRODUCED 3.3-5.3 km/day" (F5), "tidal sensitivity growing: REPRODUCED" (F6), "secondary
fronts NOT REPRODUCED, bottleneck = analysis" (F4), the injection "REOPENS fine-tuning" sentence and "anchor wrong by 7-13x"
(F7), "11 stations / 8.1 GB" (superseded by the 15-station set), "What E4 changes" items 2/5/6. NOTES E4: the "false alarms
0/1/1/2/4" sentence renamed, the step-7 3.6-MAD paragraph and decisions 1/4/5/6 marked. PLAN GOAL: "~10 km/day" ->
"~7 km/day (80 km in 11 days; the paper quotes ~5-10)". Acceptance met: no surviving RESULTS sentence relies on the
median-detections column or the 4-MAD column (both are kept, struck, for the record).

### STEP 1 — detection-budget table on the 2010 data (`e4_template_test.py --budget`, `scripts/e5_budget_table.py`)
Code (all behind `--budget`, commented at the fix sites): the network series is built for EVERY family with a template on a
station that has data that day (E4 built it only for families with a catalogued event that day, which biases a quiet-day
count); per (day, L, family, T in 5/6/7/8/10 MAD) it records declustered detections, rolled-control detections (each
station's series rolled by its own 300-3000 s BEFORE the weighted sum), detections within +-2 s of a Bostock event of the
SAME family / of ANY family (all 130 families, both halves), rolled detections matched to the same family (= chance
matching, 0.0-0.2 %), the same split into TEST blocks (odd 2-h blocks outside the 60-s guard), and the 8-MAD extras within
60 s of a NEIGHBOURING family's catalogued event vs the same for the 600 random times/day; per event a +-2 s common-lag
network score `net_mad_2s`. Templates: Bostock's published stacks on the ring + BUILD-half transfers for B011/B926/LRIV/NLLB
(re-built at 4/8/12 s: `e4_transfer.py --min-events 40 --lengths 4,8,12 --tag e5`, 82 families x 10 stations = 820
templates, 174 s CPU; trap: the script's default `--min-events 100` gives only 34 families — E4 used 40). Moveout measured on
the BUILD half (pass 1, tag `e5_mv_build`, 1,767 (family, station, L) corrections, median +0.03 s, IQR -0.02..+0.36) and
applied to the full window. Quiet day = <= 50 Bostock events that UTC day = doy 228, 229, 249, 250, 251 (doy 248 has 57 and is
ETS by the rule); 19 ETS days. Runs: `e5_budget_op` (7 stations, 127 families), `e5_budget_ring` (10, 129), `e5_budget_ring6`
(E4 six-station ring, published stacks only). ~3.5 min each on the L40S (34 days x 3 lengths).

**Bug found and fixed (run 1): the rolled control did not roll the gap mask.** `e4_template_test.py`, the `vmask &=
VALID[s][L][k[s]:k[s]+NVo]` intersection was applied to BOTH the network and the rolled series, so a station's masked gap
(cc undefined, ~1e10) landed on samples the mask called valid in the rolled series: 50.8 rolled detections/template/day at
8 MAD on doy 251 (LRIV has a 0.2 % gap that day) against 0.5-0.9 on the other quiet days. Fixed by rolling each station's
mask with its series (`vmask_s`), applied in E5 runs only (E4 bit-for-bit otherwise). E4's ring runs had no gaps at the
stations used; E4 leg C (LRIV/STW) may have carried a mildly inflated `shift_fa8` (its mean was 1.56 vs 1.95 for leg A, so
not visibly). After the fix the quiet-day rolled rate sits at or just below the quiet-day data rate (0.64 vs 0.76), as a
chance control should. All three budget runs were re-run after the fix.
Second bug: in `--budget` mode the network loop reaches families whose only SeisBench stacks sit on stations NOT selected
(LZB/YOUB/PFB for the subset) but which have a transferred template -> `KeyError (6.0, 'B011')` at the placement (`off`);
fixed by seeding `meta` with every SeisBench family's hypocentre under `--budget` and requiring a placement for a station to
join the sum.

**Results, 8-s template (per template per day; family-day means; precision = matched to own family within +-2 s / total):**

| network | T | quiet-day FA mean (p90) | quiet rolled | ETS-day total | ETS rolled | precision own / any (test blocks own) | recall of all Bostock events +-2 s (test half) |
|---|---|---|---|---|---|---|---|
| operating subset (7) | 7 | 8.5 (18.0) | 7.2 | 11.8 | 7.8 | 0.14 / 0.28 (0.12) | 0.264 (0.206) |
| operating subset (7) | **8** | **0.76 (2.0)** | 0.64 | 2.72 | 1.09 | **0.37 / 0.54 (0.32)** | **0.154 (0.117)** |
| operating subset (7) | 10 | 0.014 (0) | 0.013 | 0.48 | 0.07 | 0.65 / 0.78 (0.60) | 0.048 (0.035) |
| ring + transfers (10) | 7 | 5.3 (10.0) | 4.7 | 10.6 | 5.0 | 0.24 / 0.37 (0.22) | 0.408 (0.347) |
| ring + transfers (10) | **8** | **0.42 (1.0)** | 0.28 | 3.05 | 0.49 | **0.56 / 0.69 (0.53)** | **0.268 (0.220)** |
| ring + transfers (10) | 10 | 0.005 (0) | 0.005 | 0.76 | 0.016 | 0.84 / 0.90 (0.83) | 0.099 (0.080) |

Template length at 8 MAD, operating subset: quiet FA 0.32 / 0.76 / 1.26 at 4 / 8 / 12 s (12 s FAILS the <= 1 bar), precision
own 0.43 / 0.37 / 0.31, recall 0.076 / 0.154 / 0.192 (F12 confirmed: longer templates buy recall at the price of FA; the
FA-matched choice is 8 s at 8 MAD or 12 s at ~9 MAD). Ring + transfers: 0.19 / 0.42 / 0.62, precision 0.62 / 0.56 / 0.50,
recall 0.151 / 0.268 / 0.319. At 5 and 6 MAD the quiet-day FA is 650-700 and 70-96 per template per day (>95 % noise) —
F1's "4 MAD is FA-prohibited" is confirmed with room.
Recall by magnitude (operating subset, 8 s, 8 MAD, +-2 s, all events incl. families without a template on the subset):
M<1.5 0.125, M1.5-1.8 0.150, M1.8-2.1 0.162, M>=2.1 0.225 (test half 0.127 at M1.8-2.1); ring + transfers: 0.23 / 0.26 /
0.27 / 0.31 class (`e5_budget_table_e5_budget_ring.csv`). Recall by stations carrying a template (ring + transfers, 8 s):
1-4 stations 0.11-0.27, 5-7 stations 0.23-0.27, 8 stations 0.37, at 8 MAD.
**Precision is a station-count problem, not a threshold problem:** families with only 1-3 templates on the network give
precision 0-11 % and carry 2.3 (1 station) quiet-day FA; families with >= 4 (subset) / >= 5 (ring) stations give **58-75 % /
67-78 % own-family precision (66-82 % / 77-86 % any-family) at 8 MAD and 90-97 % at 10 MAD** (`e5_budget_bynsta_precision_*.csv`).
Every 2012 family has >= 7 of the 15 stations within 60 km (p10 = 7, median 11), i.e. the >= 5-station class.
**What the extras are:** at 8 MAD, 69 % of the unmatched detections on ETS days fall within 60 s of a NEIGHBOURING
family's catalogued event against 38 % for random times on the same days (at 5 MAD, noise-dominated, 43 % vs 38 %) — the
extras are tremor-clustered, i.e. sub-threshold LFEs Bostock's 7-station detector did not catalogue and/or cross-family
aliases, not noise; the quiet-day rate is the noise floor (0.42-0.76 of the 2.7-3.0 ETS-day detections, ~25 %).
**Acceptance:** quiet-day FA at 8 MAD <= 1 PASS (0.42-0.76; 8 s), at 7 MAD <= 10 PASS (5.3-8.5); ETS-day precision >= 85 %
FAIL on the own-family definition (37-56 %; 58-78 % for >= 4-5 stations; 84-97 % only at 10 MAD). The failure signature
"quiet FA at 8 >= 5 -> 9-10 MAD" did NOT fire. **Operating threshold carried into steps 2-3: 8 MAD (8-s template)**, with
10 MAD reported alongside as the precision-driven alternative (recall halves to a third).

### STEP 2 — corrected prediction model (`scripts/e5_prediction_model.py`; F2, F9, F11)
Fix: detection is on the RAW statistic, raw = excess + bg with bg = 2.728 (random-time network median, L = 8); recall(T, M) =
P(raw >= T) with log10(raw) Gaussian about log10(raw_med(M)) with the MEASURED spread of the raw statistic at fixed magnitude:
2010 ring, L = 8, M1.8-2.1, raw p10/50/90 = **3.39 / 5.89 / 9.51** (Merlin's 3.5/6.1/9.6 recomputed) -> sigma_log POOLED 0.175
(PRIMARY), within-family robust 0.140. Per-station published-quality excess (L = 8, M1.8-2.1): LZB 2.34, PGC 2.15, VGZ 1.15,
PFB 0.68, SNB 0.98, YOUB 1.87 (ring, published stacks) and the transfers brought to published quality by the measured
q-ratio ours(N~90)/published = 0.72 (`e4_step2_q.csv`): B011 1.65, B926 2.60, LRIV 0.54 -> class model surface < 40 km 2.15,
40-60 km 0.83, borehole 2.12 (B011/B926 at 41-47 km, both bands), > 60 km far x (60/d)^2 (assumption). Template quality
q(N) = 1/sqrt(1 + 1/(r1_eff^2 N)) with r1_eff = 0.088 CALIBRATED on the 0.72 ratio at N = 90 (the metadata r1 ~0.05
understates it, E4 step 1): q(25/50/100/200/500) = 0.40/0.53/0.66/0.78/0.89. e(12 s)/e(8 s) = 1.10 (median per-station
ratio). Network: equal weight g_eq (sum e_s)/sqrt(n) with g_eq = 0.863 measured (E4 leg A); SNR-weighted g sqrt(sum e_s^2)
with the same g (primary) and g = 1/1.38 (Merlin, conservative). Magnitude: e(M) = e(1.95) 10^(alpha dM), alpha 0.38
(E4 within-family), 0.46/0.58 (incoherent slope, F9) and 0.0/0.2 (see the measured curve below).
**In-sample (each 2010 family's ACTUAL ring station set, station-measured e_s, equal weight): predicted raw median 5.47 vs
measured 5.89 (bar 5.9 +- 0.5 PASS); frac >= 8 0.165 vs 0.206 (bar 0.21 +- 0.05 PASS with the pooled spread; 0.114 with the
within-family spread — the tail is between-family registration variance the class model cannot carry, so the pooled spread
is the right one).** The class model instead of station values gives 5.35 / 0.161. The excess is under-predicted by 1.15x
(measured/predicted), reported as a labelled bracket only. Failure signature (miss by > 2x) NOT triggered.
**Measured recall vs magnitude is FLAT pooled** (2010 ring, 8 s, >= 8 MAD: 0.18-0.25 from M1.35 to M2.25, dipping to 0.19 at
M2.0-2.2 — the distant families carry the biggest catalogued magnitudes) and rises only mildly WITHIN family (0.157 at
dM -0.35 -> 0.272 at +0.35, d log10 recall/dM = 0.36 — about half the alpha = 0.38 model's slope). The per-bin 2012 numbers
therefore inherit a model slope; alpha = 0.2 is the central case, alpha = 0 (flat) and 0.38 the brackets; the pooled OVERALL
recall is the number to hold the download to. Figure `figures/e5_recall_curves.png` (a, d).
**Bostock's arrays (his templates q ~ 1):** 21-station array, templated stations per family (median 11): Mc(8) = 1.98; all
21 within 60 km (median 19): **Mc(8) = 1.75** (target 1.6-1.9 PASS; his 2003-05 modal M 1.6-1.7). 7-station 2010-13 array
(ring + SOKB): templated (median 3/family) Mc(8) = 2.77, all within 60 km (5/family) **Mc(8) = 2.54 vs his 2010 modal 1.75 ->
our implementation-vs-Bostock deficit = +0.8-1.0 mag = 1.7-2.0x in network registration** (his 7-station matcher needs
1.7-2.0x our predicted excess to put M1.75 at 8 MAD). Stated plainly: with the same seven stations Bostock's iteratively
refined 1,400-2,000-event templates and his own combination out-register our published-stack implementation by ~2x.
**2012 15-station set** (130 families active 2012-08-15..10-15; per-family n p10/50/90 = 50/192/392, build half = n/2, template
if >= 20 build events -> 126 families; stations within 60 km per family median 11, p10 7, p90 12; q(build) median 0.66; e_net
median equal 3.92 / weighted 4.04; Mc(8) median 2.29 / 2.25). Predicted recall of his 28,149 events (pooled spread, g 0.86):

| set | L | weights | alpha | T | all | M<1.6 | M1.6-1.8 | M1.8-2.0 | M>=2.0 |
|---|---|---|---|---|---|---|---|---|---|
| 15 | 8 s | equal | 0.38 | 8 | 0.267 | 0.175 | 0.259 | 0.340 | 0.414 |
| 15 | 8 s | SNR | 0.38 | 8 | **0.279** | 0.184 | 0.273 | 0.355 | 0.429 |
| 15 | 8 s | SNR | 0.2 (central) | 8 | **0.321** | 0.259 | 0.323 | 0.370 | 0.395 |
| 15 | 8 s | SNR | 0.0 (flat) | 8 | 0.378 | 0.367 | 0.385 | 0.388 | 0.360 |
| 15 | 8 s | SNR, g 1/1.38 | 0.38 | 8 | 0.209 | 0.135 | 0.202 | 0.269 | 0.335 |
| 15 | 12 s | equal / SNR | 0.38 | 8 | 0.310 / 0.324 | 0.207 / 0.218 | 0.303 / 0.318 | 0.391 / 0.408 | 0.469 / 0.483 |
| 15 | 8 s | SNR | 0.38 | 10 | 0.093 | 0.040 | 0.082 | 0.135 | 0.206 |
| 18 (+LZB/YOUB/PFB) | 8 s | equal / SNR | 0.38 | 8 | 0.298 / 0.316 | 0.193 / 0.206 | 0.289 / 0.307 | 0.380 / 0.401 | 0.477 / 0.496 |
| 18 (+LZB/YOUB/PFB) | 12 s | equal / SNR | 0.38 | 8 | 0.345 / 0.364 | 0.228 / 0.242 | 0.336 / 0.356 | 0.436 / 0.457 | 0.536 / 0.556 |

Cost of the user's exclusion of LZB/YOUB/PFB: **-0.03 to -0.04 in overall recall (0.28 -> 0.32 at 8 s SNR-weighted;
0.32 -> 0.36 at 12 s), -0.07 at M >= 2.0**; families with >= 50 % predicted recall 9 -> 20 (8 s) / 30 -> 41 (12 s). Region-year
table redone at 8 (and 10) MAD only, q at N = 96, alpha 0.38 (`e5_prediction_regionyear.csv`): 8-MAD Mc (equal / SNR) Puget
2023+ 1.84 / 1.79 (0.43-0.46 of the 2012 event distribution recalled), N-Oregon 2023+ 1.89 / 1.86 (0.40-0.42), S-VI 2023
2.04 / 2.02 (0.33-0.34), S-VI 2025 2.03 / 1.97 (0.34-0.37), C-Oregon 2023 2.15 / 2.13 (0.29), NorCal 2025 2.15 / 2.10
(0.29-0.31), N-VI 2.36 / 2.33 (0.21-0.22); at 10 MAD every region loses half (Puget 0.26, S-VI 0.17). Download volume, no
fetch: 15 stations x 3 ch x 51 d (2012-08-15..10-05) at 40 Hz = **11.2 GB** at F1's measured 1.416 B/sample (STEIM2 at
40 Hz; E4's 1.17 B/sample would give 9.3 GB); at the native rates (12 stations at 100 Hz, 3 at 40 Hz) 17.0 GB.
Bugs in the new script caught at first run: `summ.mode` (DataFrame method shadows the column -> empty selection ->
IndexError) -> bracket access throughout; `itertuples()` mangles column names containing '<' -> `to_dict("records")`.

### STEP 4 — Thomas harness fixed (`scripts/e5_thomas_harness.py`, imports E4's geometry/tide/burst helpers; CPU, ~4 min)
PMC text (two prompts): main front "~5 km/day" (Fig. 2A caption) and "80 km along strike in 11 days or ~10 km/day", "first
manifests in the north (40 km) around 6 September"; secondary fronts "initiate within 1 to 10 km of the main front,
propagate 10 to 40 times faster, and span durations of 2 to 6 hours", "four are marked with arrows" (NO dates in the text),
Fig. 2C "a secondary front illuminated by both LFEs and tremor on 14 September that propagates at a speed of ~12 km/hour",
family 141 "six secondary fronts that occur between 12 September and 16 September" (60 % of its moment); tides: early
LFEs "not significantly influenced", later bursts "almost always tidally triggered", ~2 kPa shear stress, no model named;
median secondary-front stress drop 8 kPa, Mw ~5.1; no detection algorithm, count, completeness or tidal method. The text
anchors we can validate against are therefore (i) 14 Sept ~12 km/h and (ii) six fronts on family 141 in 12-16 Sept (family
141 = SeisBench id 141, 48.5015 N -123.8948 E 33.2 km, x = +6.7 km along strike; 128/130 families matched to SeisBench ids
by lat/lon/depth rounding). The four Fig.-2A arrows cannot be dated from the text.
(a) Main front = first BURST per family (>= 4 events, gaps < 10 min) after 08-30; propagation window from the 6-h
leading-edge track (10th percentile of event position; bins >= 10 events): onset = first bin >= 20 % of the peak bin count
(run 1 used "the first bins" and collapsed to zero length on the full catalogue because the sparse pre-front days already
put the 10th percentile at -33 km — fixed), departure >= 5 km south of the onset position, arrival within 5 km of the
southern extreme. Full catalogue: onset 09-07 04:00, window **09-07 10:26 .. 09-21 09:51 (14.0 d)**, leading edge +22 ->
-55.5 km; **first-burst Theil-Sen -4.69 km/day (bootstrap CI -5.00..-4.38, rho -0.92, 98 families in the window); leading edge
-4.83 (-5.17..-4.52), centroid -4.35; leading-edge travel 72.8 km in 14.0 d.** Like-for-like with Thomas's dates (leading edge
09-06 12:00 -> 09-17 12:00): 50 km / 11 d = 4.6 km/day (full), 62-67 km / 11 d = 5.6-6.1 km/day at Mc >= 1.7. Decimation:
first-burst -5.10 / -5.02 / -5.19 / -4.70 / -8.67 km/day at Mc 1.7 / 1.8 / 1.9 / 2.0 / 2.1 (CI half-widths 0.3-0.4 to Mc 1.9,
-6.6..-2.2 at 2.0, -140..-0.2 at 2.1); leading edge -5.06 / -5.17 / -5.18 / -5.54 / -4.19. Verdict: **4.7-5.2 km/day,
at the lower edge of the 5-10 accept band and 30 % below the 7.3 km/day of "80 km in 11 days", consistent with the paper's
own "~5 km/day"**; the estimator is robust to Mc 2.0 (8 % of events) and breaks at Mc 2.1.
(b) Secondary fronts: 3-h window, 30-min step, event-level times, >= 8 events from >= 4 families spanning >= 5 km, Theil-Sen
speed over BETWEEN-family event pairs (same-family pairs have dx = 0 and would pull the median to zero) 1-10 km/h with a
200-draw bootstrap CI excluding 0, |Spearman(x, t)| >= 0.6, overlaps merged (most events kept), distance behind the
interpolated leading-edge track as a diagnostic; nulls = circular per-family shift of event times and within-UTC-day shuffle
of family labels (20 draws each). Diagnostics added (not gates): median residual about the fitted line, core fraction and
core extent (events within 5 km of the line). **Full catalogue: 20 detections vs nulls circular 2.25 +- 1.5 (p95 4.05) and
day-shuffle 1.75 +- 1.2 (p95 3.05).** They are compact fronts inside zone-wide background: core extent median 14.6 km (5-30),
core fraction median 0.81, 3-34 families in the core, speeds 1.5-9.2 km/h (= 8-45x the 4.8 km/day main front; Thomas
10-40x), 0.6-17 km behind the leading edge (Thomas 1-10). Validation: six detections in 12-16 Sept (09-12 20:05, 09-13 22:04,
09-15 02:34, 09-16 02:07 / 08:42 / 16:35); **two of them have family 141 inside their core (5 within their full span)**
against Thomas's six on that family; **the 14-Sept ~12 km/h front is in the data at 09:34-11:34 UTC (7.1-8.5 km/h, 690-830
events from 63-69 families) but sits AT the rho gate at full completeness (rho = 0.5996 vs 0.6 required; not tuned) and is
detected at Mc 1.7 (8.5 km/h) and 1.8 (9.8 km/h).** Against "recover >= 3 of the four with <= 1 null per episode": the four
arrows are undated in the text, so the test is against the two text anchors — one recovered at Mc 1.7-1.8 (borderline at
full completeness), family 141 partly (2 of 6 by core, 5 of 6 by span); nulls 1.75-2.25 per episode (<= 1 NOT met — the
null level of this detector is ~2/episode). Not the "criteria mismatch" failure branch (the detections have Thomas's speeds,
sizes and positions), not a clean pass either. Decimation: 20 / 15 / 14 / 14 / 6 / 5 detections at Mc none / 1.7 / 1.8 / 1.9 /
2.0 / 2.1 vs circular-null p95 4.05 / 9.05 / 14.05 / 18.0 / 12.1 / 5.0 and day-shuffle p95 3.05 / 4.1 / 4.0 / 5.0 / 2.0 / 1.05:
**above the stricter (circular) null only at full completeness and Mc 1.7 (>= 52 % of events)**, above the day-shuffle null
at every cut. The circular null GROWS as the catalogue thins (fewer events -> spurious rank correlations at fixed gates).
(c) Tides (PLACEHOLDER equilibrium tide, labelled everywhere; Schuster p dropped): one phase per burst START (2,047 bursts);
R = 0.284, time-shift null (tide shifted 1-30 d, 200 draws) p = 0.07 (null p95 0.29); thirds: R 0.109 / 0.354 / 0.275 with
null p 0.65 / 0.07 / 0.475 -> **no significant tidal phase clustering at any stage; the E4 "growth" is gone**; at Mc 1.7-1.9
the whole-episode p is 0.03-0.055 and the middle third 0.025-0.07 (marginal), the last third never (p 0.2-0.5). No claim of
tidal sensitivity is made; a Cascadia interface tidal-stress model (e.g. SPOTL-class solid-Earth + ocean loading resolved on
the Slab2 interface) is the named external dependency.
(d) Re-rupture = a burst starting >= 12 h after the family's previous event: 86 re-activations, 4.2 % of bursts, **64 of 130
families (49 %)**; interval p10/25/50/75/90 = 0.53 / 0.61 / 0.70 / 1.12 / 1.84 d, 99 % within 5 d. Decimation: families
re-activated 49 -> 35 -> 33 -> 24 -> 13 -> 6 % (Mc none -> 2.1); median interval 0.70 -> 0.82 -> 0.93 -> 0.94 d stable to
Mc 1.9, then noisy.
(e) Completeness each metric NEEDS (fraction of Bostock's 2012 events): main-front speed <= 8 % (Mc 2.0; CI +-0.3 km/day to
18 %); secondary-front count above the strict null >= 52 % (Mc 1.7), above the day-shuffle null >= 3 %; the 14-Sept anchor
needs 33-52 %; re-rupture FAMILY FRACTION within 30 % of truth >= 52 %, interval median >= 18 %; tides n/a (placeholder).
Files `results/e5_thomas_{summary,criteria,families,leadtrack_*,firstburst_*,bursts_*,secondary_*,rerupture_*,windows_2012-09-14_*}`,
figure `figures/e5_thomas_harness.png`.

### STEP 1 (continued) — the E4 six-station ring in budget mode (`e5_budget_ring6`, published stacks only, build-half moveout)
8 s, 8 MAD: quiet-day FA 0.51 (p90 1.0), rolled 0.45; ETS 3.00 (rolled 0.93); precision own / any 0.42 / 0.60 (test blocks
0.41); recall +-2 s 0.195 (E4 reported 0.222 at +-1 s with the all-window moveout, 0.168-0.206 held out); 10 MAD: 0.008 /
0.73 / 0.85 / 0.067; 4 s @ 8 MAD: 0.14 / 0.53 / 0.66 / 0.107; 12 s @ 8 MAD: 0.87 / 0.34 / 0.55 / 0.247. By stations: 4-station
families precision 0.67 / 0.78 and recall 0.28; the 5-station group (families 1 and 15 + neighbours) recall 0.05 (F10).

### STEP 3 — SNR-weighted network, held out (`e5_run_step3.sh`: e5_w_{op,ring}_{eq,snr}, test half, 8 + 12 s, --budget)
Weights `results/e5_weights_snr.csv`: 10^(build-half stack S-SNR dB/20) from `e4_transfer_stackqc_e5.csv` (all 10 stations,
82 families), capped at 3x the median, median-normalised (per-station medians PGC 1.53, LZB 1.57, YOUB 1.44, B926 1.20,
SNB 1.15, B011 1.08, PFB 0.96, VGZ 0.85, NLLB 0.84, LRIV 0.54); pairs without a weight -> 1.0; the rolled control uses the
same weights. Families 1 and 15 excluded from the headline (excess 0.5-1.5 MAD, recall 0-3 % on every network).
Ring + transfers (10), 8 s, M1.8-2.1: equal excess 3.34 / frac>=8 0.242 (0.246 at +-2 s) -> **weighted 4.07 / 0.341 (0.349)**;
all magnitudes 3.29 / 0.226 -> 4.07 / 0.334; rolled >= 8 at events 0.000 -> 0.000; quiet-day FA at 8 MAD 0.42 -> 0.68;
ETS rolled 0.51 -> 0.74 /template/day (+46 %); 12 s: frac>=8 (+-2 s) 0.282 -> 0.401. Operating subset (7): 2.35 / 0.129 ->
2.94 / 0.195; quiet FA 0.76 -> 1.36; ETS rolled 1.07 -> 1.47. **FA-matched** (weighted threshold with the equal network's
quiet-day FA at 8 MAD, log-linear in T): ring 8.27 MAD -> recall 0.308 vs 0.246 (+25 %; all mags 0.301 vs 0.231; 12 s 0.358
vs 0.282); subset 8.32 MAD -> 0.177 vs 0.133 (+33 %). Weights vs test registration: Spearman 0.82 over 565 (fam, sta) pairs,
within-family median 0.86 -> not noise. Test-block precision at 8 MAD: ring 0.53 / 0.67 (equal) vs 0.47 / 0.62 (weighted);
10 MAD 0.83 / 0.90 vs 0.81 / 0.88. Acceptance: excess >= 3.3 PASS and frac>=8 >= 0.22 PASS on the 10-station network (also at
the FA-matched threshold: 4.07 / 0.308), rolled-at-events unchanged PASS, quiet-day FA "not worse" FAIL at a fixed 8 MAD
(the weighted sum has heavier tails), rolled RATE up 40-46 %; on the 7-station subset the absolute bars are not reached.
Files `results/e5_weighting_{table,fa,famatched,corr}.csv`, figure `figures/e5_weighting.png`.

### Runs (all on the L40S unless CPU): e4_transfer e5 (CPU 174 s) | e5_mv_build (pass 1, 10 stations, build half, 3 lengths,
~3.5 min) | e5_budget_op / _ring / _ring6 (pass 2 + budget, ~3.5 min each; the first two re-run after the rolled-mask fix) |
e5_w_op_eq / _snr, e5_w_ring_eq / _snr (test half, 8 + 12 s, ~2.5 min each) | e5_thomas_harness (CPU, 3 runs: window fix,
core diagnostics) | e5_prediction_model (CPU, 5 runs) | e5_budget_table / e5_weighting_table / e5_figures (CPU).
### Bugs found in E5 (code line at the fix site)
1. `e4_template_test.py` rolled control: gap mask not rolled with the series -> 50.8 chance detections/template/day on doy 251
   (LRIV 0.2 % gap); `vmask_s` added, E5 modes only. E4 leg C's `shift_fa8` may have been mildly inflated.
2. `e4_template_test.py --budget`: `KeyError (6.0, 'B011')` — transfer-only families lacked a hypocentre in `meta` (loaded only
   for the SELECTED stations' stacks); seeded from all SeisBench families; `ss` requires a placement.
3. `e5_thomas_harness.py` propagation window: "first bins" onset collapsed the window to zero length on the full catalogue
   (sparse pre-front days put the 10th-percentile edge at -33 km); onset = first 6-h bin >= 20 % of the peak count.
4. `e5_thomas_harness.py` secondary speed: same-family event pairs (dx = 0) would drag the x-vs-t Theil-Sen median to zero ->
   between-family pairs only (design, documented).
5. `e5_prediction_model.py`: `summ.mode` — the column name shadows `DataFrame.mode` -> empty selection -> IndexError; bracket
   access throughout. `itertuples()` mangles 'recall_M<1.6' -> `to_dict("records")`.
6. `e4_transfer.py` default `--min-events 100` -> 34 families (E4 used 40 -> 82); re-run with 40.
7. Process hygiene: `pgrep/pkill -f <pattern>` matched the tool's own shell (exit 144) and a waiter died with the tool shell;
   waiters now `setsid nohup`, kills by PID.
### Decisions a human should make (also RESULTS E5 "What E5 changes")
(i) go / no-go on the 11.2-GB 15-station Sept-2012 download under a predicted 27-32 % recall at 8 MAD (Merlin's >= 20 % met;
the per-bin 60/35/20 % bars are not all predicted and rest on a half-supported magnitude slope — pooled >= 25 % is the
defensible pre-registration); (ii) a >= 4-5-station reporting rule (precision 58-86 % vs 0-11 % below it) or a per-station
coincidence gate for step 5; (iii) 8 s @ 8 MAD vs 12 s @ ~9 MAD as the pre-registered operating pair; (iv) SNR weighting with
the threshold set on quiet-day FA, not on a fixed MAD multiple; (v) whether to keep secondary fronts and the re-rupture family
fraction as 2012 targets at all (both need >= 52 % of Bostock's events; the 2012 catalogue would reach ~30 %); (vi) the
whole-step-0 strike list is applied — E4's step-7 table remains in RESULTS only as struck record.
## 2026-09-11 ~06:05 UTC user "Go" -> E5 step 5 APPROVED: download 15 stations 2012-08-15..10-05 (~17 GB native / 11 GB @40 Hz) + from-scratch SNR-weighted replication vs Bostock's 28,149 events; pre-registered bars per merlin_audit_e4 §5. Agent on Fable (Opus quota).

## 2026-09-11 E5 step 5 (2012 replication) — the ONE approved download + from-scratch replication of Bostock's Sept-2012 catalog
**Agent: Fable 5.1 (Opus quota exhausted; user "Go" 2026-09-11 ~06:05 UTC).** Environment as E4/E5 (torch 2.12 cu130, L40S, CUDA-12 libpath
prepend; all CPU work <= 28 threads). Nothing of E2-E5 deleted or overwritten; every new file is `e5r_`-prefixed (catalog `e5_catalog_2012_sVI*`);
`e4_template_test.py` gained flag-guarded extensions only (`--no-seisbench --families-csv --stations-csv --gt-parquet --offsets-csv --skip-z
--catalog-out --amp-ref-csv`; default behaviour bit-for-bit; backup `scripts/.e4_template_test.py.bak_e5r`). **Fetches: (i) FDSN station
METADATA for the 15 stations (channel codes 2012); (ii) THE ONE APPROVED WAVEFORM DOWNLOAD below. Nothing else.**

### STEP A — download (`scripts/e5r_download.sh` -> `d3_download.py`, 3 stations x 2 day-workers = <= 6 concurrent requests; logs `logs/e5_dl_<STA>.log`)
Channel bands from the metadata (NOT all BH? as the brief assumed): CN PGC, CN NLLB, UW LRIV = BH? 40 Hz; **CN VGZ and CN SNB = HH? 100 Hz**
(SNB switched BH->HH 2011-03-17; VGZ has been HH since 2009-07-31 -- the 2010 VGZ traces on disk are HH too); PB B001/B003-B007/B009-B011/B926 =
EH? 100 Hz. Window 2012-08-15..2012-10-05 inclusive (52 days), all three components, native rate. 06:03-06:08 UTC, ~20 MB/s.
Result (`results/e5r_download_table.csv`, `scripts/e5r_dl_table.py`): **15.69 GB, 774 station-days with data, 6 no-data (PB.B010 doy 228-233)**;
every station 3-C, Z/N/E resolved by `e4_chanorder` (CN/UW ?HZ/?HN/?HE, PB EHZ/EH1/EH2 -> Z/N/E); coverage >= 0.96 everywhere except one day
each at B001 (0.78), B011 (0.84), B010 (0.31, its first day); PGC day files are 57 fragments/channel (merged, gaps zero-filled and masked
by the correlator's energy mask as in E4). 20-Hz 1-8 Hz cache: `scripts/e5r_cache.sh` -> `e4_cache20.py --year 2012` (unchanged), 15 processes,
~1 min per batch; cache20hz/<NET.STA>/2012/<doy>.npz.

### STEP B — ground truth, split, offsets (`scripts/e5r_groundtruth.py`, `scripts/e5r_offsets.py`)
* Bostock 2012 events in the window: **28,149 (2012-08-30 01:04 .. 09-25 22:30), 130 families** (unique lat/lon/dep). Nearest-neighbour map
  to SeisBench ids (same method as e2_groundtruth): 129 within 1 km (median 0.33 km) but three Bostock pairs 0.1-0.2 km apart in depth share
  one SeisBench id (bkeys 51/52 -> 234, 74/75 -> 257, 104/105 -> 22) and one 297-event family (48.51 N -123.92 E 32.9 km) has no SeisBench
  counterpart -> 126 distinct SeisBench ids. DECISION: the from-scratch replication keys families by **Bostock's own family (bkey 0-129)**;
  the SeisBench id is carried only for the q-control and the 2010-offset comparison.
* Hold-out: alternate 2-h blocks of t_origin = DB time + per-family offset, 60-s guard (never parity): **build 14,308 / test 13,381 / guard
  460**. Families with >= 30 / >= 100 / >= 200 build-half events: **116 / 63 / 13 of 130** (n_build p25/50/75 = 52/98/160, max 263).
* Per-family time base re-derived on the 2012 BUILD half at PGC, two estimators (log `logs/e5r_offsets_run2.log`, fig `figures/e5r_offsets.png`):
  (a) E3's stack-alignment residual histogram (peak_offset copied verbatim; per event the lag of the max 3-C cc of the SeisBench PGC stack,
  S-1..S+7, within +-40 s of the zero-offset prediction; 2-s window, Poisson z >= 4): **accepted 71 of 129 families with a PGC stack**
  (acceptance rises with N: 0/15 below 30 build events, 5/26 at 30-60, 22/26 at 60-100, 36/51 at 100-200, 8/12 above 200); offset median
  +12.25 s, IQR +11.0..+13.9, range +3.25..+21.25 (the 2010 E3 table spans the same range: it is a per-family property, E3 step 5 stands).
  (b) from-scratch envelope peak of the un-aligned PGC stack at DB + ttS(Crosson), +-40 s: accepted (>= null + 3 dB) 122 of 130.
  **(a) vs (b) agree: a - b median +0.17 s, IQR -0.14..+0.65 (70 families). DRIFT 2012 - 2010 over 63 families: median +0.50 s, MAD 0.74 s;
  39 within 0.5 s, 53 within 1 s, 2 beyond 2 s** (bkey 59 = SeisBench 15 at +5.75 and bkey 112 = SB 12 at +2.25). Run 1 -> run 2 RULE (set
  before templates were built): (a) is used only when (b) is unavailable or agrees within 2 s, else (b) -- because (a) leans on the
  SeisBench stack, which E4/E5 showed is bad for family 15 (bkey 59: offset_a 11.0 at z 7.6 vs offset_b 5.5 vs 2010 5.25 -> (b) wins).
  Working offsets: 69 (a), 54 (b), 1 from 2010 (bkey 46), 6 the median (all with < 40 build events -> at most one gets a template).
  Split rebuilt on the final offsets (counts above).

### STEP C — templates from scratch (`scripts/e5r_templates.py`, 15 processes, 8-23 s per station; `logs/e5r_templates.log`)
Method: per (family with >= 30 build events, station): 3-C 1-8 Hz 20-Hz data cut at t_origin + ttS(Crosson) +- 20 s, rms-normalised, linear
stack, alignment refined by cc of every event against the running stack (+-1.5 s, 2 iterations), S picked on the 0.5-s-smoothed 3-C envelope
within +-4 s of nominal, S-SNR (S-0.5..S+3.5 vs 15 s earlier) vs a same-N random-time stack put through the SAME alignment + pick (stricter than
E4's un-aligned null), nested sub-stack S-SNR at N = 30/60/100/200; templates 8 s and 12 s from S-1 s, per-component L2-normalised; USABLE =
S-SNR >= null + 3 dB. NO SeisBench stack is used as a template. PRE bookkeeping: the correlator places a template at origin + ttP + sp - 0.5,
so sp = (ttS - ttP) + (pick - nominal) - 0.5 -> start = origin + ttS + pick refinement - 1.0 (commented at the site).
**Coverage (`results/e5r_template_qc_2012.csv`): 1,144 usable (family, station) templates over 116 families; stations per family p10/50/90
= 6/10/14; 113 families with >= 5 stations, 91 with >= 8.** Usable / 116 and median S-SNR (dB): PGC 111 / 12.1, B010 113 / 10.8, B009 111 / 10.4,
B011 109 / 9.6, SNB 101 / 8.0, VGZ 97 / 7.1, B926 95 / 8.2, NLLB 88 / 6.4 (at 86 km!), B003 56 / 3.1, B006 52 / 2.8, B007 49 / 3.0, B001 47 / 2.6,
B004 41 / 2.2, B005 40 / 2.4, LRIV 34 / 2.0 (null 0.0-1.0 dB everywhere). Usable fraction by N: 31 % (< 30 events), 48 % (30-60), 59 %
(60-100), 75 % (100-200), 82 % (> 200). S-SNR growth (usable, median of nested sub-stacks) N = 30/60/100/200: surface 5.5/7.3/8.6/8.9 dB,
borehole 4.7/6.1/7.1/8.6 dB. **q-control (our 2012 re-stack vs the SeisBench stack, 8 s from S-1, +-1 s): cc median 0.78 at PGC (93 pairs,
86 >= 0.6), VGZ (60, 53 >= 0.6) and SNB (22, 18 >= 0.6); N >= 200: PGC 0.81 (7 of 8 >= 0.6), VGZ 0.66 (4 of 6), SNB one pair at 0.40** ->
Merlin's >= 0.6 bar met at PGC/VGZ, the single SNB N >= 200 pair fails it (family-level, not a sign error). S pick minus the Crosson S: 0.0 at
PGC by construction of the offset, -0.15..+0.1 at the co-located B009/B010/B011/B926/VGZ, +0.7-0.8 at NLLB/SNB, **+1.7-2.0 at LRIV/B005-7 and
+2.2-3.3 at the far southern boreholes B001/B003/B004 (52-69 km)** -- a systematic path-model residual at the Olympic-Peninsula stations that
sits near the +-4 s pick search edge (flagged: the far-borehole templates may be picked on a later phase; not tuned).
Weights (`results/e5r_weights_snr_2012.csv`, BUILD half only): w = 10^(S-SNR/20), capped at 3x the median, median-normalised; per-station
medians PGC 1.64, B010 1.42, B009 1.32, B011 1.27, B926 1.16, SNB 1.04, NLLB 1.00, VGZ 0.95, LRIV 0.77, B003 0.75, B004 0.71, B007 0.70, B005 0.69,
B006 0.67, B001 0.64. Amplitude reference (`e5r_amp_ref_2012.csv`): per station the median sliding 3-C window energy (8 / 12 s) over the quiet days.

### STEP D — network matched filter (`scripts/e5r_run_mf.sh` -> `e4_template_test.py` with the E5-step-5 flags; GPU L40S, 4 torch threads)
Pass 1 (`e5r_mv_build`, `--half build`, 4.3 min): per-(family, station, L) residual moveout = median best lag at the BUILD-half events (+-3 s
scan): **2,288 corrections, median +0.00 s, IQR -0.03..+0.04 s** -- the from-scratch stack alignment had already absorbed the moveout that E4
measured at +0.41 s against Crosson + the published S-P; the correction is nil but is applied for method parity with E5 step 3.
Pass 2 (`e5r_snr`, then the equal-weight reference `e5r_eq`; `--half all --budget --budget-thr 7,8,9,10 --catalog-out --amp-ref-csv`,
`--weights e5r_weights_snr_2012.csv` for the primary; ~5 min each): 15 stations x 1,144 templates x 2 lengths x 52 days; per (family, day)
the network statistic = SNR-weighted channel mean of the per-station 3-C cc (weights x channel count, normalised by the contributing weight
that day), background = median/MAD of that family's own daily network series (gaps masked), declustering = max-pool over one template length,
rolled control = each station's series rolled by its own 300-3000 s BEFORE the weighted sum (mask rolled too, E5 fix). Every declustered
detection >= 7 MAD is written with its network MAD, raw network cc, SNR-weighted relative amplitude (cc3_s x sqrt(E3_s / E3_quiet_s), i.e. the
template-projected amplitude in units of each station's pre-SSE quiet-day window rms, weighted mean over the contributing stations) and
+-2 s matches; thresholding that list at T >= 7 reproduces the T-level declustered set exactly (the peak set is threshold-independent).
Cross-family declustering (the catalog, `e5r_score.py`): a detection is kept iff its network MAD is the maximum over ALL families' detections
within +-2 s (Bostock-style max-cc assignment; MAD-normalised so families with different station counts compare). First run: 260,215 (SNR) /
132,520 (equal) detections >= 7 MAD -- but see bug 5 (quiet days skipped) -> pass 2 re-run in full (`e5r_run_mf2.sh`) for a self-consistent file.

### Bugs found in E5 step 5 so far (code line at the fix site; all minimal, all commented)
1. `e5r_groundtruth.py` (first draft): `g.loc.fillna` -- `loc` is a DataFrame attribute (the memory's bracket-access trap) -> `g["loc"]`.
2. Channel bands: the brief assumed BH? 40 Hz for CN/UW; the metadata says **VGZ and SNB are HH? 100 Hz in 2012** (SNB switched 2011-03-17).
   Fetched at the native rate (the E4 cache builder resamples 100 -> 20 Hz by resample_poly); volume 15.7 GB instead of the estimated 17.
3. `e5r_offsets.py` run 1 -> run 2: offset rule made stack-independent when (a) and (b) disagree by > 2 s (family 15 / bkey 59, see STEP B).
4. `e4_template_test.py:~452` RuntimeWarning "invalid value encountered in cast" (NaN offsets for families without a template at a station,
   cast to int64) -- pre-existing E4 behaviour, masked downstream by `isfinite(o)`; benign, left as is.
5. **`e4_template_test.py` (`if not len(evd): continue`): a day with NO catalogued event was skipped entirely**, so in `--budget` mode the 15
   pre-SSE quiet days (0 Bostock events) produced neither detections nor budget rows -> quiet-day FA = NaN in the first scoring, and the catalog
   lacked its quiet-day false alarms. E5's 2010 quiet days had <= 50 events, so this never surfaced. Fixed behind `--budget` (event-free days
   are processed with an empty event set; `ev_rows` only appended when events exist); pass 2 re-run in full for both weightings.
6. The empty-day path then raised `TypeError: bitwise_and` at `e4_template_test.py:456`: `np.array([...])` of an EMPTY list is float64, so
   `sel & (...)` failed; dtypes made explicit (`dtype=bool` / `np.int64` / `float`). Second re-run clean.
7. `e5_thomas_harness.py:170` (pre-existing E5 harness): the within-day-shuffle null on OUR 12-s catalog produced a secondary-front window with
   no event within 5 km of the fitted line -> `min()` of an empty array. Crash guard only (core_x_min/max = NaN when empty; diagnostics, not gates);
   Bostock's runs never hit it, so no E5 number changes.
8. `e5r_dl_table.py` / the download queue: two relative-path slips after `cd` in backgrounded chains (tool cwd resets) -> absolute paths everywhere.

### STEP D results (final, after the quiet-day re-run; `results/e5r_detections_e5r_{snr,eq}.parquet`, `e5_budget_e5r_*.csv`, `e4_tt_events_e5r_*.parquet`)
SNR-weighted: 8 s -> 30,201 catalog rows at >= 8 MAD after cross-family declustering (7,009 within +-2 s of a Bostock event, 5,100 same family);
12 s -> 33,968 (7,824 / 5,739). Equal-weight reference: 14,400 (5,160 / 3,969). Quiet-day detections per template per day (mean over 1,627
family-days; rolled in brackets), SNR 8 s: T7 29.5 (9.4), **T8 4.47 (0.96), p90 10**, T9 0.63 (0.10), T10 0.08 (0.01); 12 s: 30.0 / 4.72 / 0.73 /
0.12; equal 8 s: 11.5 / 1.14 (0.39) / 0.10 / 0.01. ETS-day (>= 100 Bostock events) detections per template per day, SNR 8 s: T8 12.4 (rolled
1.03, 8.3 %), T9 4.7 (0.15, 3.3 %), T10 2.2 (0.03, 1.5 %); ETS precision own / any at 8 MAD 0.27 / 0.53 (equal 0.43 / 0.64), at 9 MAD 0.47 /
0.72, at 10 MAD 0.63 / 0.85. Runtime: pass 2 ~6.5 min per weighting (52 days, 2 lengths, 15 stations, 1,144 templates; 24 GB GPU).

### STEP E — scores (`scripts/e5r_score.py`; `results/e5r_score_e5r_snr.json` + csv; the PASS/FAIL table is in RESULTS.md)
Test half 13,381 events (280 in 14 template-less families). 8 s / 8 MAD SNR: recall own 0.234 (M<1.6 0.218, 1.6-1.8 0.236, 1.8-2.0 0.234,
>= 2.0 0.286), catalog 0.346 (chance 0.064; 0.307 / 0.338 / 0.371 / 0.442), +-5 s 0.265 / 0.424; own-family rolled-at-events 0.000. 12 s: own 0.252
(0.230 / 0.250 / 0.259 / 0.314), cat 0.365. 9 MAD: own 0.146, cat 0.226 (chance 0.028). 10 MAD: 0.087 / 0.139. Equal 8 MAD: 0.155 / 0.259.
Measured / predicted (own; SNR alpha 0.38): 0.84 pooled, 1.18 / 0.87 / 0.66 / 0.67 by bin; alpha 0.2: 0.73; 12 s 0.78; 9 MAD 0.52. Recall by
usable-template count: 1-4 0.076, 5-7 0.171, 8-10 0.274, 11-15 0.240 (own). Per-family own recall p10/50/90 0 / 0.168 / 0.419; Spearman with
n_sta 0.42, n_build 0.45, magnitude -0.15. Station registration at test events (M1.8-2.1, excess over random median): PGC 2.11, B010 1.87,
B009/B011/B926 1.66, NLLB 1.12, LRIV 1.04, VGZ 1.01, B007 1.01, B004 1.00, B003/B006 0.93, B005 0.78, SNB 0.69, B001 0.63. Extras at 8 MAD:
23,192 (77 %); 38 % on >= 100-event days; amplitude 0.19 vs 0.65 matched; 25 % within 60 s of another family's event vs 17 % random;
62 within 60 s of a catalogued earthquake (74 EQs in the window +-15 d; chance 0.1 %).
**Quiet-day diagnosis (the decisive finding).** (i) Bostock's catalog has no events outside ETS windows in 2010-2013 (2012: 0 before 08-30 / after
09-25; 2011 only Aug-Sep 18,663; 2013 only Sep-Oct 14,873) -> the quiet-day bar's premise was untested. (ii) 7,272 quiet-day detections: by day
228:184 229:191 230:85 231:7 232:32 233:236 234:0 235:1674 236:1113 237:1399 238:579 239:76 240:1013 241:528 242:155; 68 % at 04-11 UTC; PNSN
tremor in the box 0 on every quiet day (13 on 08-31); 15 within 60 s of a catalogued EQ (0.2 %); same-day time clustering (near-60-s fraction
0.74-0.90 vs Poisson 0.44-0.78) like ETS days; families spread over the whole zone (std lat 0.15, lon 0.24 deg; median 34-37 km from PGC).
(iii) quiet rate per family vs the PGC/B009/B010/B011 weight fraction: Spearman 0.52; <= 0.2: 0.10 /day, 0.2-0.35: 0.32, 0.35-0.5: 3.0, 0.5-0.7:
5.1, > 0.7: 9.5; worst families bkey 103 (46.5 /day, cluster 0.69), 32 (40.3, 0.83), 75 (27.7), 128 (22.9, 0.78). (iv) **one-station-per-site run**
(`e5r_run_onepersite.sh`, tag `e5r_snr_onepersite`: PGC, VGZ, NLLB, SNB, LRIV, B001, B003, B004, B926; doy 228-244; 8 s; 2 min GPU): quiet-day data /
rolled = T7 8.18 / 8.14, **T8 0.85 / 0.81, T9 0.11 / 0.08, T10 0.01 / 0.01 vs the full network 29.5 / 9.4, 4.47 / 0.96, 0.63 / 0.10, 0.08 / 0.01** ->
the excess over the rolled control is entirely the co-located clusters (site-coherent noise adds in the sum; independent rolls destroy it).
Post-SSE 09-26..10-05: ours 3.2x rolled, PNSN box 34 tremor detections in ten days (north of 49.3 N: 0). ETS 09-06..09-24: ours / Bostock 0.56,
matched / Bostock 0.251, ours / rolled 6.7.
Daily counts (SNR 8 s / 8 MAD; ours / matched / rolled / Bostock / PNSN box): 09-07 1058 / 488 / 69 / 2033 / 339; 09-15 1019 / 699 / 92 / 3129 / 448;
09-18 1196 / 726 / 66 / 2696 / 404; 09-26 1861 / 0 / 363 / 0 / 33 (full table `e5r_score_e5r_snr_daily.csv`).

### STEP F — Thomas harness on our catalog (`scripts/e5r_run_thomas.sh`; tags ours2012_L8, ours2012_L12, bostock2012_dec = his catalog randomly
decimated to our catalog-level recall 0.346 -> 9,789 events, seed 5; `--mc none`; ~30 s each) -> `results/e5r_thomas_comparison.csv`
first-burst front ours -4.48 km/day (CI -4.94..-3.98, rho -0.80, 93 families) / 12 s -4.71 vs his -4.69 (-5.00..-4.38) / dec -4.73 -> +-20 % MET;
leading edge ours -0.49 (window 08-30 13:25 .. 09-21 09:09, 21.8 d, 18.7 km) vs his -4.83 (09-07 .. 09-21, 72.8 km) -> not reproduced (floor-smeared
track); secondary fronts 13 (nulls circ 0.10 / p95 1.0, day-shuffle 1.05 / 2.1; speed median 6.4 km/h, core 17.7 km, 10.3 km behind) vs his 20
(2.25 / 4.05; 1.75) and dec 26 -> +-1 NOT MET; 14-Sept anchor found (1) in ours, not in his at full completeness; fam-141 12-16 Sept 3 span / 1 core
vs 5 / 2; bursts 1,179 vs 2,047 (dec 761); re-activations 40 vs 86 (dec 51); **families re-activated 0.25 vs 0.49 (dec 0.29)** -> +-0.1 NOT MET vs his,
within 0.04 of the decimated reference; interval p10/50/90 0.53 / 0.73 / 1.72 d vs 0.53 / 0.70 / 1.84; tides placeholder R 0.23 p 0.64 vs 0.28 p 0.07.
12-s catalog: 7 secondary fronts, re-rupture 0.22, front -4.71.

### STEP G — deliverables
`results/e5_catalog_2012_sVI.parquet` (30,201 rows; columns time, family_id, sb_family_id, lat, lon, depth, rel_amp, magnitude, n_stations,
network_mad, cc_best, template_len_s, matched_bostock_id [row of e5r_gt_bostock2012.parquet, NaN if none within +-2 s], matched_bostock_mag,
matched_same_family, doy) and `_L12` (33,968); README block + caveats in RESULTS.md. Figures `figures/e5r_{download_coverage,template_snr_vs_N,
offsets,quietday_fa,recall_by_mag,daily_counts,front_tracks,secondary_fronts}.png` (`scripts/e5r_figures.py`). Scores `results/e5r_score_e5r_{snr,eq}*`.
Post-fix: `matched_bostock_idx` (an index into the t_origin-SORTED ground truth) replaced by `matched_bostock_id` (file row) in both catalogs.

### Runs (chronological, UTC): e5r_download.sh 06:03-06:08 (15 stations) | e5r_cache.sh x2 06:08-06:12 | e5r_groundtruth.py x3 | e5r_offsets.py
run 1 06:11, run 2 06:13 (rule) | e5r_templates.py 06:14-06:17 (15 procs) | e4_template_test.py e5r_mv_build 06:20-06:24 | e5r_snr 06:24-06:29,
e5r_eq 06:29-06:34 (quiet days skipped -> superseded) | e5r_score x2 06:34 (superseded) | e5_thomas x3 06:34-06:36 (superseded) | re-run e5r_snr
06:40-06:47, e5r_eq 06:47-06:54 (final) | e5r_score x2 06:54 | e5_thomas x3 06:54 | e5r_figures 06:57 | e5r_snr_onepersite 06:56-06:58 | catalog
id fix + RESULTS/NOTES 12:30 (session rate-limited ~07:00-12:30; nothing re-run on the GPU after 06:58).

### Decisions a human should make
1. The pre-registered replication FAILS at 8 MAD on the FA bars (quiet 4.47 vs 1; rolled 8.3 % vs 5 %; extras 38 % vs 90 %) and on 3 of 4 recall
   bars on the own-family definition (0.234 vs 0.25; 0.286 vs 0.60; 0.234 vs 0.35; M1.6-1.8 0.236 vs 0.20 met); at 9 MAD the FA bars pass and every
   recall bar fails (0.146). Accept the verdict as is, or re-pre-register a site-merged network (the co-location finding) -- one GPU pass, no download.
2. Co-located clusters (PGC/B009/B010/B011; LRIV/B005/B006/B007) must be merged at the site level and the chance control must roll sites: a design
   rule for every fleet region with PB boreholes next to CN/UW broadbands.
3. The E5 prediction model over-predicts the magnitude slope (measured 1.3x from M<1.6 to M>=2.0 vs 2.3x at alpha 0.38) and its FA calibration
   does not transfer to clustered geometries; the pooled signal prediction (0.84x) is fine.
4. The far Olympic-Peninsula boreholes (B001/B003/B004/B005/B006/B007) and LRIV register at 0.6-1.0 MAD with S picks 1.7-3.3 s after the Crosson
   S -- whether a path-specific travel-time (or a wider pick search) recovers them is an open implementation question.
5. Thomas science at this completeness: main-front speed yes (first-burst estimator), re-rupture interval yes; secondary-front count and re-rupture
   family fraction no (as E5's decimation predicted) -- the science replan Merlin flagged ("secondary fronts need >= 35-50 % of Bostock density") stands.
Note: the catalog rematch by file row gives 7,010 / 5,100 (8 s) and 7,825 / 5,738 (12 s) matched any / own against the scorer's 7,009 / 5,100 and
7,824 / 5,739 -- one tie at exactly 2.0 s resolved differently; all RESULTS numbers are the scorer's.

## 2026-09-11 ~13:15 UTC Merlin audit of the 2012 replication -> results/merlin_audit_e5r_2026-09-11.md
Replication FAILED honestly. Site-noise diagnosis right in direction (positive tests: 44.5% of quiet-day detections have another family
within +-2 s at 10.8 km median / 21% >30 km = annulus around PGC; 30-60x nocturnal) but mechanism unproven (no SITE-rolled control; one-
per-site halves recall too). Sensitivity is GEOMETRY-limited (recall flat in template N above 50; 0.27-0.36 within 30 km of PGC vs 0.17
beyond 40 km); on-disk ceiling ~0.30-0.35 own recall; Bostock-class needs sites within ~30 km of each family. Honest precision: 8 MAD
~20-25% overall / 40-55% ETS days; 9 MAD ~45-55% / 65-75%. Honest product = detector with validated chance model at 25-35% family-
dependent completeness. Prescribed E6 (no downloads): 1 site-rolled control; 1b positive stack test at excluded VGZ/SNB; 2 site-merged
network (2 variants, pre-registered FA bars); 3 recall levers (pick search, >=5 dB subset, 2010+2012 templates); 4 scoring fixes
(chance-corrected, per-family thinned reference, density-robust secondary statistic); 5 re-scope Thomas science to main-front speed,
re-rupture intervals, first-burst migration, completeness-weighted rates; post-2019 download gated on 1-2 + pre-registration. Awaiting user go.
## 2026-09-11 ~13:30 UTC user "Ok" -> E6 launched (Fable agent; steps 1, 1b, 2, 3, 4 + re-scope documentation; no downloads).
## 2026-09-11 ~13:45 UTC GPU outage announced (a few hours out): E6 reordered GPU-first; checkpoint file results/e6_checkpoint.txt.
