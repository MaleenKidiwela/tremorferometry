## GOAL (user, 2026-09-10 22:50) — Thomas et al. 2025 (Sci. Adv. adz2832) at PNW scale, full record
Do what Thomas et al. did for ONE SSE (Sept 2012, southern Vancouver Island, Bostock family-located LFE catalog + PNSN tremor: main
front ~7 km/day (80 km in 11 days; the paper quotes ~5-10; corrected E5 step 0); secondary fronts 10-40x faster, 2-6 h, re-rupturing segments within days; per-family LFE rate decay; tidal
sensitivity evolving through the SSE; secondary-front stress drops ~8 kPa median) for EVERY SSE, over ALL of Cascadia/PNW
(BC-WA-OR-N.CA), over the ENTIRE available record (2005/2009-2026). The catalog is the means; the space-time dynamics of slow slip
(fronts, secondary fronts/RTRs, re-rupture, healing, tidal modulation) across cycles and segments is the science.
Catalog requirements this imposes (beyond "picks"): (1) event-level origin times (s) AND locations (few km) or, where the network
is too sparse, family-level locations with event-level times (Bostock style); (2) family/cluster assignment so per-family rates and
re-rupture can be measured; (3) an amplitude -> moment/magnitude module (Bostock 2015 style template-amplitude scaling) for
moment release and stress drop; (4) completeness that is STABLE across time and network changes (or a documented, time-varying
magnitude of completeness) so cycles can be compared; (5) Slab2/interface geometry at each source for tidal stress and
along-strike/down-dip framing. Stage plan below (E1 done -> E2 pilot -> F1 inventory -> fleet decision) stands; the analysis
endpoint above decides the pass criteria. The tremor_slip_geometry study is the PNSN-tremor precursor of the same questions.

# PLAN — LFE catalog repair (2026-09-10)
D1 separability diagnostic (this round): for ~12 stations already streamed in S14, per template: cc distributions (0.80-1.00, 0.005
   bins, per-day rates) in ETS weeks vs quiet weeks (ETS/quiet defined from PNSN tremor within 40 km of the template's source cell:
   30-d count > 90th pct vs < 25th pct); excess = ETS rate - quiet rate per bin; where does the excess sit? Precision/recall vs
   threshold; noise-tail exponential fit from quiet weeks; per-template verdicts; network-coincidence pilot on a few source cells.
D2 (if D1 positive) full streaming pass: per-template thresholds, declustering, network coincidence -> lfe_catalog_repair/results/
   lfe_events.parquet; validation vs PNSN tremor and Lin/Bostock families.
D3 (if D1 negative) design of a waveform re-detection (MAD thresholds, 3-C cc stacking, network stacking) and its cost.
Amendment (user, 2026-09-10): QUIET weeks are noise + inter-ETS LFEs, not pure noise. D1 therefore also uses (i) a < 10th-percentile
tremor reference, (ii) a NOISE-ONLY floor model (exponential tail fitted on cc 0.80-0.85 only; and the pooled cc-histogram shape of
unresponsive templates scaled to each template's low-cc rate), (iii) an inter-ETS LFE candidate count (quiet-time detections above the
floor; burstiness and tremor correlation as the reality check). D2 thresholds must come from the noise-only floor (target false-alarm
rate), never from ETS-minus-quiet, so inter-ETS LFEs are kept.
Compute note (user, 2026-09-10): a GPU is available. D2 parsing is I/O-bound (use cuDF if installed for CSV parsing/histograms). D3
re-detection should run on the GPU via the existing scripts/discover_gpu.py matched-filter path (~200x CPU), with per-template MAD
thresholds and 3-C cc stacking added; apply the CUDA-12 LD_LIBRARY_PATH fix after any pod restart (memory: cupy-cuda12-libpath).

## After D1 (2026-09-10): D2 dropped; D3 design
D3 pilot (proposed): re-download 2 years (one full ETS cycle each) for 4-5 stations — WISH + B943 (north, surface + borehole, same
cell families), one central-Oregon station (YACT or ALSE), one Vancouver Island station covering Lin/Bostock families (ground truth).
Re-detect on GPU (scripts/discover_gpu.py path): 3-C cc stacking per template; per-template threshold = median + 9 MAD of the
daily cc series (recomputed per 30-d window to follow the floor drift); network stack of cc across stations for same-cell families;
declustering (template duration); per-detection amplitude/SNR, polarisation, and station count stored. Validation: Lin/Bostock family
occurrence times; ETS/quiet ratio >= 5; tremor correlation; borehole/surface agreement. If the pilot passes -> full fleet re-detection.

## D3 stage 4 — family LOCATION pilot (user 2026-09-10: "continue with the catalog refinement and location")
Region with ground truth: southern Vancouver Island / northern Puget (PGC + Lin/Bostock families). 6-8 extra fleet stations within
~100 km with 2017-2018 data (B011, CN/PB). Network-consistent templates = stacks of the extra stations' 3-C data at the PGC detection
times, aligned by Slab2-predicted S times in a 1-D model then refined by cross-correlation; P (Z) and S (N/E) picks on stacks;
grid-search family hypocentres (0.5 km) with Slab2 as weak prior; bootstrap uncertainties; comparison with Lin/Bostock locations and
the nominal discovery cell; event-level relative locations from cc-lag differences for >= 4-station events (Bostock 2012 approach).
Pass criterion: family hypocentres to ~2-3 km with the fleet geometry; then fleet-wide.

## D3 stage 5 — NETWORK matched filter + template length (orchestrator, 2026-09-10, from the stage-2b calibration)
Root cause of the false-alarm floor: a 2-s single-station template has ~24 independent samples (sigma_cc ~ 0.17), so cc >= 0.80 is
only ~6 MAD; 3-C stacking buys sqrt(3). Remedy = network matched filter (sum of per-channel cc over all stations with a network-
consistent template, aligned by the family moveout; sigma ~ 0.04 for 8 stations x 3 comps) and longer templates (4 s, 6 s around S).
Stage 5 compares single-station 3-C vs network-stacked detections on the same 60-day window (ETS/quiet, tremor correlation,
>= 3-station fraction, Lin recall/precision where possible) and template lengths; delivers the FLEET design recommendation.
Fleet cost sketch from the pilot: GPU ~48 min per 2-year station at ~1,000 templates x 3-C (~5 station-years/h); 3-C waveforms
~13 GB per station per 2 yr (40 Hz) -> ~20 TB for 195 stations x 17 yr: must be streamed per station (download -> detect -> delete),
~7 h/station serial -> ~2 months serial, ~2-4 weeks with download/GPU pipelining. Network detection multiplies template-channels
but shares the same waveform pass. Decide after stage 5.

## CORRECTED PLAN (2026-09-10, after user pushback) — deep-learning picker catalog for the whole network
The template-matching limit (1.8 MAD per station) does NOT apply to a deep-learning picker that uses the multi-second, 3-C LFE signature:
Lin et al. 2024 catalogued > 1 M LFEs on southern Vancouver Island with a single-station picker + association + 3-D grid-search location
on the same station types. Our June picker failed for fixable reasons (3.5k training windows, origin-time labels, no association).
E1 (running): obtain Lin's model/code/training data, the Bostock catalog with arrival picks, tools (seisbench, pyocto/GaMMA, NonLinLoc),
   a 3-D velocity model.
E2: run the picker fleet-wide on streamed 3-C data (GPU), associate (PyOcto/GaMMA), locate (grid search / NonLinLoc, Slab2 prior).
E3: validate vs Lin (VI), tremor + borehole-surface consistency elsewhere; completeness map by region.
Fallback: retrain the same architecture on Bostock's precise picks + our June hard negatives (EQ/blast/in-tremor noise).
Timeline ~3-4 weeks pipelined; the download (20 TB streamed) is the cost; picking is fast on the GPU.

## E2 pilot (2026-09-10) — LFEDetect picker -> PyOcto -> NonLinLoc on data already on disk
Picker: SeisBench LFEDetect (Münchmeyer et al. 2024; weights `cascadia` and joint `jcms`; ZNE 20 Hz, 60-s windows, 1-8 Hz). Data on disk:
PGC 2010-05..2010-11 (Lin/Bostock ground truth), 9 VI stations 2017-02-20..04-21, WISH/B943/YACT 2017-18. Add 4-5 VI stations for
2010-05..2010-11 (within cap) so association can be tested against Lin/Bostock. Steps: picks -> PyOcto (N >= 3 VI, >= 4 elsewhere) ->
NonLinLoc with casc1.6 (S-weighted) -> recall/precision vs Lin and Bostock, location scatter vs Bostock family hypocentres, GPU throughput.

## F1 (2026-09-10) — fleet inventory & association regions (no waveforms): station table with 3-C availability per year (FDSN metadata),
source-cell station density (60/100/150 km) per year, 5-7 association regions, locatable fraction per region-year, data/GPU cost estimate,
ground-truth/fine-tuning coverage per region. Runs in parallel with the E2 pilot.

## E3 (2026-09-11, Merlin audit of E2 — results/merlin_audit_e2_2026-09-11.md) — on-disk controls BEFORE any fleet download
Status: steps 0-6 DONE 2026-09-11 ~01:55 (RESULTS.md E3) — step-1 verdict: E2 PASSes RETRACTED (chance-dominated at N>=4 thr 0.2/0.3); step 6 TRIGGERED, awaiting Merlin re-audit (~02:00) + user go. Steps 1-5 first (one Opus agent, no downloads); 6 conditional on 1; 7 = decision gate.
1 shifted-pick PyOcto null (is the associated catalog real or combinatorial?)  2 relocate Bostock families through the same NLL
(picker vs geometry vs S-only vs model)  3 injection recall-vs-SNR + false-alarm rate on real quiet days (can fine-tuning buy anything?)
4 Bostock templates as classifier on associated events; MADs of Bostock's own events with his own templates on the 6-station subset
5 per-family Bostock time base + re-score  6 (if 1 shows chance dominance) network back-projection of LFEDetect S-probability traces
7 decision gate: product = event times + family label + relative amplitude + completeness curve; hypocentres by-product; fine-tune only
if 3 shows the false-alarm floor binding (real Cascadia noise negatives, day-edge class, family-held-out split); FIRST SCIENCE MILESTONE
= replicate Thomas 2025 Sept-2012 SSE on southern Vancouver Island with our catalog (targeted VI download 2012-08..10, tens of GB);
fleet decision after. Also: F1 readiness at 40 km; correct off-GT and Lin-attribution text in RESULTS/NOTES.

## E4 (2026-09-11, Merlin re-audit of E3 — results/merlin_audit_e3_2026-09-11.md) — corrected template test + Thomas harness, on disk
Status: STEPS 0-7 DONE 2026-09-11 ~02:40 (RESULTS.md E4). VERDICT: with the components fixed and an empirical
moveout, Bostock's own templates DO register on the six-station 2010 ring (+2.6..+3.4 MAD over the random-time
background, 12-26 % at 8 MAD, shifted-control chance 0.02 %) -> E3 step 4 RETRACTED, back-projection NOT triggered,
template route OPEN but threshold-limited. Also: the 0.078 s/n anchor is wrong by 7-13x; the network must be
SNR-weighted (equal-weight mean got WORSE when 4 stations were added); a POLARIS 2010 download is IMPOSSIBLE
(stations closed 2004-2007). Original status line: LAUNCHED 2026-09-11 ~02:50 (user "Ok"); Opus agent on steps 0-6 (outer-ring template test first, POLARIS only as control); Merlin vetting the post-2019 redesign in parallel. Component-order bug (stacks Z,N,E; scripts Z,E,N) voids E3 step 4 and taints steps 3/3b.
0 fix component order  1 corrected Bostock-template test at Bostock event times, empirical per-station moveout, 2/4/8-s templates
(DECIDING: template route open if median >=3 MAD at M>=1.8 and >=30% recall at 8 MAD with <=10% chance)  2 corrected injection +
incoherent s/n re-anchor (closes or reopens fine-tuning)  3 association re-null with tolerances matched to the 0.4-0.6 s residual
4 Thomas-2025 analysis harness built and validated on Bostock's own Sept-2012 catalog (zero download)  5 metadata availability of the
21 template stations for 2012  6 CONDITIONAL on 1: targeted Sept-2012 southern-VI download (10-20 GB), our catalog, Thomas metrics
7 back-projection only if 1 fails. Also correct RESULTS.md step 4 / S10 sentence, flag 3/3b, fix mad_median.
REDESIGN (user objection 2026-09-11 ~02:40: "POLARIS isn't going to work for our available stations after 2019"): validate on the station
class and spacing we will USE (CN outer ring + PB boreholes, all still operating) in 2010-2013 where Bostock/Lin ground truth exists;
POLARIS = positive control only. Tests: T1 corrected template test on the 2010 outer ring + boreholes; T2 template TRANSFER to stations
Bostock never had (B011, B926, NLLB, LRIV, OCP; build Aug 2010, test Sep 2010) = the mechanism for every post-2019 station; T3 family
DISCOVERY with no prior catalog on 2017 VI data (autocorrelation / DL-seeded), checked against Bostock's family map; T4 CHAINING across
station turnover (2010 -> 2017). Post-2019 pipeline per region: discover families on the dense network in ETS windows -> per-station
templates by transfer -> network matched filter with per-template MAD thresholds -> family locations from stacked P/S picks -> relative
amplitudes -> per-family per-year completeness curve. Dense catalog = 2019-2026 cycles where density is Bostock-class (Puget 12, S-VI 9,
N-OR 9 within 40 km; C/S-OR + MTJ 6); activity index elsewhere/earlier. T3/T4 wait for Merlin's vet (running) before being handed over.
MERLIN VET of the redesign (2026-09-11 ~03:05, results/merlin_vet_e4_redesign_2026-09-11.md) — APPLIED to the running E4 agent:
the fleet decision is QUANTITATIVE, not pass/fail: (i) T1-T2 -> a calibrated SENSITIVITY MODEL (registration in MAD vs magnitude, distance,
station class, template length, channel count; compared to r1·q·A(M)/sigma predicted from the stack metadata; expected ~1.5-2.5 MAD at
M>=1.8 on the 6-station ring); (ii) Thomas harness -> COMPLETENESS-REQUIREMENT curve by decimating Bostock's Sept-2012 catalog (Mc 1.7-2.1);
(iii) PREDICTION TABLE per region-year (8-MAD completeness magnitude, fraction of Bostock density, per-family counts per SSE) crossed with
(ii) = the download/fleet decision document. Design fixes: hold-out = alternate 2-h blocks (Aug/Sep split impossible); transfer stations
B011/B926/LRIV(+STW Z); SOKB 2010 has no data (ring = 6); T2a q-control (re-stack own 2010 data at Bostock times, cc>=0.6 vs SeisBench
stack); network-SUM threshold per template per day normalised by channel count (never per-station); T3 discovery on 2010 not 2017 (DL-seeded
clustering by differential S-arrival first; 12,703 families as seeds only after a viability cc test; controls: split-half, shifted null,
tremor-free window); T4 chaining with leave-one-station-out refinement + drift controls. Scoping: Bostock-class density = Puget 2021+,
S-VI 2023+ (not 2019+). First-download rule: POLARIS only if real events >=2x below prediction AND injection at prediction; Sept-2012 VI
only if predicted Mc<=2.0 AND fronts survive at that completeness (cheaper first: ETS-peak days 2009/11/12/13 to grow templates); 2024-08
VI window only after chaining passes.

## E5 (2026-09-11, Merlin re-audit of E4 — results/merlin_audit_e4_2026-09-11.md) — FA budget, corrected prediction, harness fixes, then 2012
Status: steps 0-4 DONE 2026-09-11 ~05:50 on Fable (Opus quota exhausted; RESULTS.md E5, NOTES 'E5 steps 0-4'). Operating point 8 MAD / 8 s (quiet-day FA 0.4-0.8 per template per day; precision 58-86 % for >= 5-station families, 37-56 % overall); corrected model in-sample PASS, Bostock 21-station Mc(8) 1.75, our pipeline ~1.7-2x below his on 7 stations; 15-station 2012 prediction 27-32 % recall at 8 MAD (11.2 GB at 40 Hz), LZB/YOUB/PFB exclusion costs 3-4 points; SNR weighting +25-33 % recall at matched FA; harness: front ~5 km/day (accept-band edge), 20 secondary fronts vs null ~2 (14-Sept anchor at the rho gate), no tidal signal, re-rupture in 49 % of families; secondary fronts and the re-rupture family fraction need >= 52 % of Bostock's events. STEP 5 NOT LAUNCHED — awaits the user's explicit go. 0 correct the record (FA -> "detections"; strike 3.6/4-MAD, "front 3.3", "tidal growth", "secondary
not reproduced"; GOAL "~7 km/day")  1 detection-budget table on the 2010 ring (thresholds 5-10; quiet vs ETS days; matched vs rolled)
2 fix prediction model (raw vs T; recall curves; SNR-weighted; predict the actual 2012 set; region-year at 8 MAD)  3 SNR-weighted network
held-out on 2010 (blocklist fams 1, 15)  4 Thomas harness fixes (first-burst front; leading-edge-tracked secondary detector validated vs
Thomas's four; burst-level tides with time-shift null; re-rupture by burst) + decimation  5 CONDITIONAL DOWNLOAD Sept-2012 S-VI, ~17
stations, 40 Hz, 12-17 GB, from-scratch replication of Bostock's 28,149 events with pre-registered pass numbers (recall >=25% overall,
>=60% at M>=2.0; rolled <=5%; quiet <=1/template/day)  6 T3 discovery on 2010 in parallel; T4 after; DL fine-tuning deferred indefinitely.
USER DECISION 2026-09-11 ~05:00: the 2012 replication station set EXCLUDES the C8 array (BPCB, SOKB, PA01-03) and the closed CN stations
LZB, YOUB, PFB. SET = PGC, VGZ, NLLB, SNB (CN); LRIV (UW); B001, B003, B004, B005, B006, B007, B009, B010, B011, B926 (PB) = 15 stations,
all operating through >=2020 (14 through 2026) and all former dv/v fleet stations (B006/B007 have no family products). E5 steps 0-4
launched on disk; the download (step 5) requires a separate explicit user go after the 15-station prediction is shown.
STEP 5 STATUS (2026-09-11, Fable agent; RESULTS.md 'E5 step 5', NOTES 'E5 step 5 (2012 replication)'): download DONE (15.7 GB, 15 stations,
774 station-days; VGZ/SNB are HH? 100 Hz); from-scratch templates 1,144 (family, station); catalog results/e5_catalog_2012_sVI.parquet.
Pre-registered verdict at 8 MAD / 8 s SNR-weighted: recall own-family 0.234 (cat 0.346; chance 0.000 / 0.064) = 0.84x predicted, magnitude-FLAT;
FA bars FAIL (quiet-day 4.47 / template / day, rolled 8.3 %, extras 38 % on ETS days) -- cause = site-coherent noise via the two co-located
4-sensor clusters (one-station-per-site run: quiet detections = rolled control); at the FA-compliant 9 MAD recall 0.146 / 0.226. Thomas: main-front
speed reproduced (first-burst -4.48 vs -4.69 km/day), secondary-front count (13 vs 20) and re-rupture family fraction (0.25 vs 0.49; 0.29 for his
catalog decimated to our recall) NOT reproduced, as E5's decimation predicted. Next (human decision): site-merged network, one GPU pass, no download.

## E6 (2026-09-11, Merlin audit of the 2012 replication — results/merlin_audit_e5r_2026-09-11.md) — FA mechanism, site-merge, re-scope
Status: LAUNCHED 2026-09-11 ~13:30 UTC (user "Ok"; Fable agent). No downloads. 1 site-rolled chance control (mechanism)  1b positive stack test at excluded VGZ/SNB
(noise vs real inter-ETS LFEs)  2 site-merged network, variants A/B, pre-registered FA bars (quiet <=1/template/day, site-rolled/data
<=1.3, ETS site-rolled <=5%), sensitivity bars unchanged  3 recall levers (pick search +-8 s + station residual; >=5 dB subset; 2010+2012
templates)  4 scoring/reference fixes  5 decision rules: "Bostock-equivalent" not reachable on this geometry (needs sites within ~30 km of
each family); honest product = validated-chance detector at 25-35% completeness; Thomas science RE-SCOPED to main-front speed, re-rupture
intervals, first-burst migration, completeness-weighted rates; post-2019 VI 2024-08 download only after 1-2 pass + pre-registration.
