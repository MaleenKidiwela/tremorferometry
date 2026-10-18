# Merlin audit of the E2 pilot verdict — 2026-09-11 ~00:30 (Fable 5.1; read-only advisor; verbatim)

**Verdict.** The E2 pilot is honest engineering with a mis-read scorecard: the two headline FAILs (recall/precision, event
locations) are mostly benchmark artifacts, and the two headline PASSes (ETS/quiet 5-34x, rho 0.9 of *associated* events) are
unverified and, on the numbers in the files, plausibly dominated by chance association of S-only picks. The single decision that
matters most is not "fine-tune or not" but whether the PyOcto event catalog is real or combinatorial — a shifted-pick null that
takes minutes on data already on disk. Second, the pipeline is optimizing the wrong product for the stated GOAL: Thomas-style
dynamics need family-assigned event times and relative amplitudes, not 25-km event hypocentres; the recommended fine-tune (on the
same stacks the `cascadia` weights were trained on) cannot fix an amplitude-limited recall and should not precede the on-disk
controls below.

## Flaws, severity-ordered
**F1. No chance-association null; the activity PASS is unproven and likely inflated — FATAL for the PASS claims (until tested).**
e2_associate.py: n_picks=3, n_s_picks=3, n_p_picks=0, pick_match_tolerance 2.5 s, velocity tol 2.0 s. Three S picks vs four unknowns
are always fittable. The only nulls in e2_assoc_score.py shift the GROUND TRUTH (chance-MATCH null), not the picks. Per-station S
rates at thr 0.3: 118-404/day (network ~2,500/day; 1,250-1,500 quiet, 3,500-4,950 ETS). Crude Poisson chance triples: ~45/day quiet
(observed N>=3 ~50/day); ~5.5/day N>=4 (observed 11/day) -> quiet-day catalog consistent with pure chance; ETS-day catalog bounded
by chance at roughly the observed level (850 vs 480 at N>=3; 330 vs 196 at N>=4). log(events/day) vs log(network picks/day) slope
2.3/3.4/3.7 for N>=3/4/5 (pure chance predicts ~N) — suggestive not decisive. rho 0.93 of associated events adds nothing beyond
single-station pick rates (rho 0.4-0.9 per station). Precision vs Bostock 0.15 (0.27 at 5 s) / vs Lin-hq 0.06 is what a
chance-dominated catalog looks like. Everything (operating point, ETS/quiet PASS, rho PASS, cloud statement) rests on this.

**F2. The recall FAIL compares single-station continuous recall to Lin's balanced-test-set numbers; recall is amplitude-limited — MAJOR.**
Lin's 90%/95% are held-out test-window numbers with an equal noise class, not catalog numbers; Lin's catalog is a 15-s N>=3
coincidence catalog of which his own quality cut keeps 9.8%. PGC recall vs Bostock (<40 km, thr 0.1) by magnitude: 0.056 (M<1.5),
0.098, 0.142, 0.179, 0.241 (M 1.95-2.1), 0.281 (M>=2.1); by distance 0.162 (<25 km), 0.110 (25-40), 0.032 (40-60). Median LFE
s/n 0.078 at PGC -> even M2.1 is s/n ~0.4. Per-station recall vs Lin-hq (thr 0.1): SNB 0.17, YOUB 0.16, B011 0.14, PFB 0.13,
PGC 0.12, LZB 0.08, NLLB 0.07, B926 0.06, VGZ 0.02, LRIV (87 km) 0.01. No re-weighting lifts a picker above the SNR of the data.
The +-1 s excess method itself is sound (0.5/1/2 s capture 76/93/100% of the peak; no P-as-S secondary peak).

**F3. The precision FAIL (0.6-4%) is an artifact of scoring against incomplete catalogs; quiet-day rates bound precision at >=24-63%
on ETS days — MAJOR.** 1 - quiet/ETS = 0.24/0.36/0.46 (cascadia thr 0.1/0.3/0.5), 0.28/0.50/0.63 (jcms). The picker is not the
precision problem at single-station level; the false-alarm FLOOR (~110-140 picks/day at thr 0.3) is, because it feeds F1.

**F4. The location "cloud" null is wrong; "network-limited not method-limited" is unsupported — MAJOR.** e2_loc_compare.py draws the
null uniformly over the bounding box of located epicentres; S-only 3-station solutions collapse toward the network centroid, where
Bostock's families are. 6.6->3.2 km with 3->7 stations is also what chance events do. OT peaks of associated events are broad
(2.5-5.5x bg vs 14-20x for arrival peaks) -> a +-2 s OT tolerance is unfair with 25-km mislocation; match on arrival times or >=5 s.

**F5. The fine-tune recommendation re-trains on the training set, with a leaky split, to fix a problem it cannot fix — MAJOR.**
cascadia.json: "trained on LFE stacks from Cascadia (Bostock et al. 2015)" = the 1,817 proposed stacks. metadata.csv split is by
trace: all 129 families and 21 stations appear in >1 split -> family-leaky validation. Recall is SNR-limited (F2). Genuinely new
information: real Cascadia noise negatives (INSTANCE noise is Italian), the day-edge artifact class, real single-event positives.
Those target the false-alarm floor — right lever for F1 — but only after an injection/FA test, with a family-held-out split.

**F6. Off-GT "surface flat / borehole responds" is distance-to-source, not generalization/instrument — MAJOR (mis-attribution).**
PNSN tremor days 2017-18 within 40/60/80/100 km: WISH 26/62/104/134, YACT 19/44/72/94, B943 109/146/173/210. Inside the 2017 VI
source zone (thr 0.3, tremor within 40 km): surface SYMB 4.17, PGC 3.13, CLRS 2.84, VGZ 2.56, SNB 2.13 ETS/quiet (rho 0.5-0.77) vs
co-located borehole B011 2.32. PGC vs B011 recall vs Lin-hq: 0.115 vs 0.137 (thr 0.1), 0.037 vs 0.064 (thr 0.3) = 1.2-1.7x, not
categorical. WISH/YACT cannot test generalization; correct NOTES before it seeds a borehole-only design.

**F7. Bostock time-base offset is station/distance-dependent; a single +12.4 s smears Bostock scores — MAJOR for Bostock numbers.**
Per-station peaks vi2010: +12.1..+13.1 s at 8 stations, PFB +9.4, LRIV +7.9; PGC +12.4 (<40 km) but +6.6/+9.1 (40-60 km). Behaves
like a reference-arrival-based time. Fix: per-family (or per family-station) offsets from residual histograms. e2_groundtruth.py
lines 230-247 still write the chance-dominated C.

**F8. F1's 60-km detectability radius overstates usable geometry ~2x — MINOR-to-MAJOR.** Measured single-station recall 16/11/3% at
<25/25-40/40-60 km -> effective radius ~40 km; recompute frac_ge4/ge6 at 40 km. F1 recommends template matching in R7 while the
project retired it — cannot both stand (see step 4).

**F9. VGZ Lin-referenced recall anomalous (0.017, peak -2.4 s; vs Bostock 0.069 at +12.1 s) — MINOR.** Exclude VGZ from Lin scoring
until resolved.

**F10. Minor.** ETS/quiet inside the 34-d vi2010 window uses a centred 30-d rolling sum with min_periods=5 -> "quiet" = window edges,
the 34x is edge-driven. Fixed-window offsets from cascadia thr 0.3 applied to both models (harmless). E2/F1 GPU-hours differ by snapshot.

## Answers A-F
A. Mostly benchmark artifact with a real SNR limit underneath (F2, F3, F4). Münchmeyer's own Cascadia numbers: SPECULATIVE (PDF not on
disk) — recollection: their DL Cascadia catalog was a small fraction of Bostock's, framed as uncatalogued LFEs. Executor: fetch arXiv
2311.13971, record Cascadia event count vs reference, PyOcto params, any chance-association control.
B. Unseparated in E2. Controls on disk: (i) locate the 129 Bostock families from the SeisBench stack P+S picks through the identical
NLL setup with all 21 stations (model+method error; stack rms 0.09 s), then vi2010-overlap stations S-only at 3/4/5 stations
(geometry+S-only; 83/30/6 families have >=3/4/5 such stations: PGC 107, LZB 100, VGZ 67, YOUB 42, SNB 26, PFB 20); (ii) locate the
step-1 null events identically; (iii) dh for time-matched N>=5 vs unmatched. Depth = S-only 3-station degeneracy + grid-top pinning
(LOCGRID from 0 km) — expected.
C. Not fine-tuning first. Minimal deciding experiment = shifted-pick PyOcto null (step 1). Then family-relocation control (2), then
injection recall-vs-SNR + false-alarm-on-real-noise (3) — the only test that says whether fine-tuning has anything to buy. Fleet
as-is -> event rate ~F^N in the per-station false-alarm rate = a tremor proxy, not an LFE catalog. Alternative: network
back-projection of the continuous S-probability traces (step 6).
D. Distance, not generalization (F6). Formal: recompute off-GT ETS/quiet with tremor within 60-100 km (expect flat); report PGC vs
B011 as the surface/borehole control.
E. Thomas 2025 used Bostock FAMILY hypocentres + event-level template times + template-scaled amplitudes. Fronts need family spacing
~5 km and times to ~10 min; tides need times to ~30 min and hundreds of events/family; stress drops need relative moments. Event-level
few-km hypocentres NOT required; family assignment, completeness stability, relative amplitude ARE. E2 optimizes the wrong thing.
Re-target: DL detection -> family classification by network cross-correlation against family templates at the DL-predicted time ->
per-family times + cc-scaled amplitudes (a classifier on ~100-300 candidates/day, not a detector on 86,400 s/day — D3's
false-alarm-floor argument does not apply). Where no families exist (R3-R6), bootstrap families from DL-associated events by
stacking (Bostock 2012's autocorrelation step, DL-seeded).
F. See F1-F10; most consequential: no association null (F1), precision vs incomplete catalogs (F3), uniform-box location null (F4),
constant Bostock offset (F7), F1's 60-km radius (F8).

## Prescription (no downloads; all inputs on disk)
1. **Chance-association null.** Add to e2_associate.py a per-station constant time offset (300-3,000 s, independent per station,
   preserves daily rates and ETS modulation, destroys real inter-station coincidence). vi2010 cascadia thr 0.2/0.3/0.5, vi2017
   cascadia 0.3, 5 realizations each; existing N>=3/4/5 post-filter + e2_assoc_score rates; feed thr-0.3 null events through
   e2_locate_nll + e2_loc_compare for the null nearest-family distance and depth histogram. Acceptance: observed vs null events/day
   per (thr, N), ETS and quiet days separately, null mean+spread; null nearest-family median. Branch: null >=50% of observed at
   N>=4 -> chance-dominated: retract PASSes and cloud statement, no fine-tune, go to steps 5 and 6. Null <=15% -> PASS stands;
   precision 0.15-0.27 then reflects catalog incompleteness/out-of-family events; steps 2-4 refine the product.
2. **Family relocation control.** NLL obs files from seisbench_lfestacks_cascadia/metadata.csv (trace_P/S_arrival_sample, 20 Hz)
   for 129 families; locate through the exact e2_locate_nll control in (a) all 21 stations P+S, (b) vi2010-overlap stations P+S,
   (c) same stations S-only incl. every 3- and 4-station subset; compare to published family hypocentres. Interpretation: (a)
   within 3-5 km / depth 5 km -> model+method fine; (c) tens of km + depth collapse -> E2's 25 km is geometry+S-only; (a) >10 km
   -> Crosson-1D/NLL setup wrong first.
3. **Injection recall-vs-SNR + false-alarm rate on real noise.** Bostock stack waveforms (PGC, YOUB, LZB families; PGC stacks for
   B011) scaled to 0.1/0.2/0.3/0.5/1/2 x host 1-8 Hz RMS, injected into real quiet days (zero Bostock events, e.g. PGC May-Jul
   2010), hundreds per station per SNR; both weight sets; recall within +-1 s, bias/scatter, pick rate on un-injected windows
   (F per day per threshold). Overlay with the recall-vs-magnitude curve via s/n 0.078. Interpretation: injection matches real curve
   -> SNR-limited, fine-tune can only lower F (compute how far F must fall for N>=4 chance to be negligible from step 1). Injection
   far above real curve -> stacks vs real single-event variability -> fine-tune must use REAL single-event positives. Caveat: stacks
   are in-sample for cascadia -> optimistic.
4. **Template-as-classifier on associated events.** For vi2010 N>=4 thr-0.3 events (and step-1 null events): network cc of 3-C data
   at the 6 overlap stations vs each of 129 family templates (stack P/S moveout), +-3 s around the predicted arrival; best family,
   its cc, cc distribution. Same at Bostock event times (step-5 offsets) -> how many MADs Bostock's own events register with his own
   templates on this subset. Acceptance: fraction of events whose best cc exceeds the null-event 99th percentile; MAD distribution at
   Bostock times. Clear excess -> hybrid product viable; Bostock events at >=6-8 MAD -> D3's "3.5 MAD" was a template-quality
   artifact, Bostock route open in R1/R2; no excess and ~3 MAD -> D3 stands, family assignment needs more stations.
5. **Per-family Bostock time base.** For families with >=50 events in pgc2010: peak of (PGC S pick - Bostock time - ttS) -> offset
   table; re-run assoc_score and loc_compare with it; fix/delete the C estimator in e2_groundtruth.py; report at +-2 s, +-5 s and
   arrival-time matching. Failure: offsets vary >5 s with no pattern -> convention needs the authors; keep Lin-hq primary.
6. **(Conditional on step 1 chance dominance) Network back-projection of S-probability traces.** Save the LFEDetect_S annotation
   stream (20 Hz float16) per station-day for vi2010; 5-km Slab2-constrained grid; moveout-shift and sum (Crosson 1-D); per-node
   threshold from quiet-day statistics or per-station-shifted null; score vs Bostock/Lin-hq as in step 1 at matched events/day.
7. **Decision gate + re-targeting (no fleet download before this).** Product = (i) event times, (ii) family/cluster label, (iii)
   relative amplitude, (iv) documented completeness curve; hypocentres a by-product. Fine-tune only if step 3 shows F binding, with
   real Cascadia noise negatives, the day-edge class, family-held-out split. FIRST SCIENCE MILESTONE: replicate Thomas 2025's
   Sept 2012 SSE on southern Vancouver Island with our own catalog (Bostock's catalog + Thomas's results = empirical pass/fail) —
   a targeted download of the VI network for 2012-08..10 (tens of GB), not 17.6 TB; the fleet decision follows. Also recompute F1
   readiness at 40 km (F8), correct the off-GT text (F6) and the Lin attribution (F2) in RESULTS/NOTES.

## What would change Merlin's mind
Shifted-pick null <=10% of observed N>=4 ETS-day rate AND null events not concentrated near Bostock families -> activity PASS and
cloud statement genuine; priorities become steps 2-5 with a stronger case for a regional R2 run. Injection recall at s/n 0.3-0.5 far
above the real curve -> waveform variability not SNR is the ceiling -> fine-tune on real positives moves to the top. Bostock's own
events at ~3 MAD with his own templates at these six stations -> D3's retirement stands; family route needs a denser network (~2018+).
