# Merlin re-audit of E3 — 2026-09-11 ~02:15 (Fable 5.1; read-only advisor; verbatim)

**Verdict.** Step 1 is sound and its retraction of E2's PASSes stands; the ML route, at any chance-clean operating point, recovers
~1-3% of Bostock's events and cannot serve the GOAL. But E3's most consequential negative — step 4's "Bostock's own events register
at 0.60 MAD with their own templates, D3 stands" — is UNFOUNDED: the classifier (and the injection, and the denoiser test) put the
stacks' horizontal components in the wrong order. The stacks are stored Z-1-2 with 1 = N and 2 = E (verified from P-wave
polarization: 222 stacks with P snr > 6 dB and rectilinearity > 0.7, horizontal P azimuth agrees with path_backazimuth_deg to a
median 9° (91% within 25°) if slot 1 = N, vs 40° (34%) if slot 1 = E), while e3_template_classifier.py line 183, e3_inject.py lines
120-123 and e3_denoise.py lines 187-188/307-308 order the data Z, E, N via sorted(). The single decision that matters is therefore not
fleet-vs-no-fleet but whether Bostock-style network template matching registers his events on the 2010 data on disk once the
components are fixed and the moveout is empirical — that on-disk test decides whether the Thomas product is reachable on Vancouver
Island at all, and it costs a GPU-hour, not a download.

## What stands, what falls
E2 stands: picker sees LFEs (PGC S-pick excess 20-50 sigma, -0.12 s); single-station recall SNR-limited (0.11->0.26 over M 1.44->1.90);
throughput. E2 falls: ETS/quiet and rho PASSes; 7.4x cloud; "network-limited not method-limited"; operating point thr 0.2-0.3/N>=4;
fine-tune-on-stacks; "FAIL vs Lin 90/95%".
E3 stands: step 1 (bias caveats F7); step 2 configs A and C; step 5 per-family offsets (r=0.959, SS_family 0.565 vs SS_station 0.009);
step 0 Münchmeyer numbers; step 3b fidelity finding (DeepDenoiser deletes 88-93% of the LFE — orientation-independent); step 6 docs.
E3 falls: step 4 in its entirety ("D3 stands" -> "untested"); step 3's interpretation (the "paradox" is mostly an x-axis anchor problem
plus the swap); step 2's "25 km is the association" softened to "association plus pick outliers, on N=3 only"; step 3b recall/FA rows
inherit the swap (conclusion unchanged).

## Flaws, severity-ordered
F1 — Horizontal components swapped in classifier, injection and denoiser tests. FATAL for step 4; MAJOR for step 3/3b. waveforms.hdf5
  data_format/component_order = Z12; LFEStacksCascadiaBostock2015 default component_order="Z12"; scripts build [Z]+sorted(non-Z) ->
  BHZ,BHE,BHN = Z,E,N. Template N dotted with data E; injected LFE seen 90°-rotated by the picker (cascadia.json.v1 is ZNE by channel
  code). The 0.60 MAD, the weak excess, "D3 stands" and the 11x injection shortfall are all measured with a rotated signal.
F2 — "matches Münchmeyer Fig. S10" is a misreading. S10 says their UNCATALOGED DL detections correlate with family stacks like noise;
  it says nothing about catalogued Bostock events failing with their own templates — those exceeded his network-cc threshold by
  construction. Step 4's result contradicts the construction of the catalog it tests -> implementation red flag, not corroboration.
F3 — Step 4's placement cannot deliver network coherence even without F1: template position = origin + Crosson ttP + stack S-P, ONE
  common lag over +-3 s; station-to-station Crosson ttP error ~0.3-0.5 s decorrelates all but one station at 2-4 Hz. Per-station test
  "NOT run". 8-s template (S-1..S+7) dilutes a signal coherent over ~2 s. Near-uniform best-lag is exactly what F1+F3 predict.
F4 — Step 3's "paradox" is mostly the real curve's x-axis: the 0.078 anchor = rms(stack)/rms(null-stack) on a 130-event stack = COHERENT
  amplitude; jitter/variability push it below the single-event amplitude. With F1, "real LFE behaves like s/n 0.22-0.25" is the expected
  direction. "Injection cannot establish an SNR ceiling; fine-tune on real positives" is premature. F-reduction arithmetic fine.
F5 — Config D partly circular (picks within +-1 s of truth-predicted arrival removes association error AND every outlier > 1 s); D at
  N>=4 is n=1; "N>=4 gives 3.7 km" rests on synthetic C. Under-read number: of 15,931 Bostock events, 123 have >=3 and ONE has >=4 near
  stations with a real pick within +-1 s at thr 0.1 — the picker essentially never registers a catalogued LFE at four stations of this
  network. That, not PyOcto, is why N>=4 Bostock recall excess is 1.8-2.9%.
F6 — Untested cheap lever: association tolerances 5x wider than the measured residual (pick_match_tolerance 2.5 s, edt_pick_std 1.5 s,
  zlim 0-70 km vs S residual sd 0.40-0.56 s). Chance ~ tolerance^(N-1) x fittable volume; ~1.2 s and 20-50 km depth should cut N>=4
  chance several-fold at ~5% real-pick loss.
F7 — Step-1 null biases both directions, net robust: constant shifts de-align within-burst co-elevation (null LOW; hourly variant +4-14%);
  pick consumption by real events (null HIGH). Neither moves 0.59 into <=0.15. N>=4/5 numbers come from n_picks=3 runs filtered on
  n_stations; exponent fit on 3 points; N>=5 ETS/quiet rests on 3 quiet-day events.
F8 — Bookkeeping: e3_classifier_summary.csv still carries buggy mad_median (-5.70); step 5's "independent window" lies inside the PGC
  window (cross-validation is across stations, not time — still the right test); clean catalog not screened vs EQ/blast catalogs.

## Answers
A. Accept step 1 and the retractions. Clean operating points are >=90% non-chance but "real" != "Bostock LFE": chance-corrected Bostock
  recall 0.9%/1.6% (N>=5, +-2/+-5 s), 1.8%/2.9% (N>=4). Order 1e-2 of Bostock.
B. Config D not clean (F5). Still: tighten tolerances to the measured residual, restrict depth to the interface band, do NOT require P
  (sigma_P 4.2 s), run PyOcto at n_picks = intended N, re-null. Improves purity for the DL catalog's only viable role (family seed).
C. Paradox explanation ranking: (1) anchor is coherent amplitude; (2) F1 rotation; (3) stack smoothing; (4) chance-bias excluded by the
  peak-excess estimator. Decide by: re-run injection with Z,N,E; measure per-event INCOHERENT s/n directly (rms² ratio S window vs 20 s
  earlier, per magnitude bin, no stacking); inject own N-event stacks (5/20/100/all) at matched snr_rms. Fine-tuning closed until then;
  if reopened: family-held-out split by source_family_id; positives = single-event windows at Bostock times aligned by (fixed) cc.
D. Reconciliation candidates: F1, F3, 2.8 vs ~7 stations x 3 comps, Bostock's iteratively refined templates vs compiled stacks. Settling
  test = corrected step 4 on 2010 data (Prescription 1), scaled by sqrt(channels) to his network.
E. Thomas product needs hundreds of events/family and times to ~10 min; ML route gives 0.1-0.2 events/family/day. Ranking: (ii)/(iv)
  network template matching FIRST — the only detector class that produced 2,500/day here, gives family label + relative amplitude for
  free, validated on disk in 2010 vs Bostock's 18,114 events; (i) back-projection second; (iii) sparse ML catalog REJECTED for the GOAL.
  PNW-wide, the DL clean catalog seeds families where none exist. ON-DISK FACT: Bostock's catalog covers the Sept-2012 SSE — 28,149
  events, 130 families, peak 3,129/day (2012-08-15..10-15) -> the Thomas analysis harness can be built and validated with ZERO download.

## Prescription (cheapest invalidation first; no downloads before step 6)
0. Fix component order in e3_template_classifier.py (183), e3_inject.py (120-123), e3_denoise.py (187-188, 307-308): select channels
   explicitly Z, then N/1, then E/2; print resolved order per station-day. Acceptance: Z,N,E for every CN station; polarization check
   reproduces ~9° under 1=N. Stop if any station has both N and 1 channels.
1. CORRECTED STEP 4 on vi2010 — THE DECIDING EXPERIMENT. (a) per-station per-own-family cc at Bostock times (step-5 origin + Crosson
   ttP + stack S-P), independent lag scan +-3 s per station, vs random times -> per-station MAD by magnitude bin (M<1.5, 1.5-1.8,
   1.8-2.1, >=2.1) and median best lag per (family, station) = empirical moveout. (b) network cc with empirical moveout, common +-1 s
   lag; template lengths 2 s (S-0.5..S+1.5), 4 s, 8 s. (c) median/p90 MAD at Bostock events by magnitude, fraction >=6/>=8 MAD, FA per
   template per day at 6/8 MAD from random times and per-station shifted control. Acceptance ("template route open"): median >=3 MAD
   for M>=1.8 on the 6-station subset with 2-4 s templates, AND recall of Bostock M>=1.8 at an 8-MAD network threshold >=30% with
   shifted-control chance <=10%. Failure: median <1.5 MAD after fix + empirical moveout -> public stacks do not register on this
   subset; go to step 5, rank back-projection up.
2. Corrected injection + re-anchored real curve: e3_inject.py at PGC with Z,N,E; per-event incoherent s/n from the 8,554 Bostock
   events with epi_PGC<40 km; own N-event stacks (5/20/100/all) at matched snr_rms. Acceptance: curves agree within 2x -> SNR-limited,
   fine-tuning CLOSED. Failure: injection >> real after both fixes -> fine-tune on template-aligned real positives, family-held-out.
3. Association re-null with tight tolerances: pick_match_tolerance 1.2 s, edt_pick_std 0.5 s, zlim 20-50 km, n_picks = 4 and 5
   directly, thr 0.3/0.5, 5 constant-shift + 2 hourly seeds; Bostock recall by magnitude with per-family offsets; EQ/blast cross-match.
   Acceptance: null/obs <=0.10 at N>=4 thr 0.3 with obs >=60/day. Seed role only.
4. THOMAS-2025 HARNESS ON BOSTOCK'S CATALOG (on disk, Sept 2012): per-family rates, along-strike main-front position vs time,
   secondary-front detection by Thomas's criteria, re-rupture, tidal-phase statistics (needs a tidal stress model — check disk; a small
   dependency, not waveforms). Acceptance: reproduce Thomas's main-front speed, secondary-front count/speed class, tidal-sensitivity
   evolution within their uncertainties. Failure: cannot reproduce from Bostock's own catalog -> the ANALYSIS is the bottleneck.
5. Metadata-only availability query (NRCan + IRIS) for the 21 template stations, 3-C >=40 Hz, 2010-08..09 and 2012-08-15..10-15
   (stack metadata station_network_code unreliable — labels CN stations "PO"). No waveforms.
6. CONDITIONAL on step 1 passing: targeted Sept-2012 VI download (CN PGC/LZB/VGZ/YOUB/SNB/PFB/NLLB, PB B011/B926 + any operating
   template station; ~1 GB/station at 40 Hz -> 10-20 GB). Templates for stations without stacks by stacking 2010 on-disk data at Bostock
   times; corrected network matcher; score vs Bostock's 28,149 events (recall by magnitude, shifted control); run the step-4 harness on
   OUR catalog. Acceptance: >=50% recall of Bostock M>=1.6 with <=10% chance and Thomas metrics reproduced. Failure: recall <20% at
   M>=1.8 -> network-limited in 2012; re-plan around 2018+ geometry.
7. Back-projection (E2-audit step 6) ONLY if step 1 fails: LFEDetect_S streams for vi2010, 5-km Slab2 grid, moveout-stack, per-node
   threshold from shifted control; <2x gain over thr 0.3/N>=5 = failure.
Also: correct RESULTS.md step 4 (retract "D3 stands" + the S10 sentence), flag step 3/3b as swap-affected, fix mad_median column.

## What would change Merlin's mind
With Z,N,E fixed, empirical per-station moveout and 2-4 s templates, if Bostock's M>=1.8 events still register below ~1.5 MAD median on
the six 2010 stations -> template matching with the public stacks is dead here, D3 re-instated for this subset, DL route (back-projection
first) the only on-disk path, GOAL out of reach at PNW scale. If the corrected injection lands on the re-anchored real curve ->
fine-tuning closed for good.
