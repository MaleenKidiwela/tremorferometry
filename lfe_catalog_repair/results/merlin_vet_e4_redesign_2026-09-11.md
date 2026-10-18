# Merlin vet of the post-2019 validation redesign (T1-T4) — 2026-09-11 ~03:05 (Fable 5.1; read-only; verbatim summary)

**Verdict.** Principle right, redesign salvageable, but it inherits the wrong DECISION STRUCTURE from the E4 step-1 bar: it asks a binary
"does the method register on our station class?" when the fleet decision is quantitative — what completeness magnitude does the
post-2019 network deliver, and what completeness does each Thomas metric need? Both answerable with zero download. Convert Tests 1-2 into
a calibrated per-station SENSITIVITY MODEL (registration vs distance, station class, template event count, template length, channel
count) and convert the Thomas harness into a COMPLETENESS-REQUIREMENT curve by decimating Bostock's own Sept-2012 catalog; downloads are
then justified by prediction-meets-requirement, not by pass/fail on a 6-station subset. Three designs do not work as written: the Aug/Sep
hold-out (ETS over by Sept 1), the "within 10 km of a Bostock family" discovery criterion (always true), and "per-template per-station MAD
thresholds" (reinstates the single-station false-alarm floor D3 measured).

## Flaws
F1 The E4 step-1 acceptance bar is mis-calibrated (MAJOR/FATAL for branch logic). From stack metadata (trace_S_snr_db, trace_stack_count)
   the implied single-event amplitude SNR r1 = 10^(snr/20)/sqrt(N): PGC 0.065, LZB 0.068, YOUB 0.058, VGZ 0.048, SNB 0.043, PFB 0.034
   (agrees with D3's 0.078 at PGC). Dead POLARIS inner stations were 2x better: KELB 0.141, SOKB 0.134, KHVB 0.096 — Bostock's advantage
   was proximity/site. cc noise per channel 1/sqrt(2BL) (0.134 at 4 s, 0.109 at 6 s); network gain 1.38x short of ideal; only 362 of
   129x6 (family, ring-station) templates exist (2.8 stations/family). Prediction for a CORRECT implementation on the 6-station ring:
   median ~1.5-2.5 MAD for M>=1.8 at 4-8 s templates. The bar "median >=3 / failure <1.5 -> back-projection" can fire failure on a
   working method. Replace before interpreting.
F2 Aug-build/Sep-test impossible (MAJOR): 2010 window peaks 2,551 events on 08-25, 57 on 09-05; top-15 families have 235-327 events in
   Aug and 1-35 in Sep; only 6 families >=50 in both. 35% of events have a same-family neighbour within 30 s -> never split by parity.
F3 Test 3 criterion uninformative and 2017 is the wrong data (MAJOR): ~1 family per 30-40 km² -> a random point is within 10 km of 8-10
   families. 2017 has NO event-level truth (Lin ends 2017-02-21, Bostock 2013). Test discovery on 2010 (18,114 Bostock events).
F4 "Per-template per-station MAD thresholds" reinstates the single-station FA floor (MAJOR): threshold on the NETWORK SUM per template
   per day, normalised by the number of contributing channels that day (d3_stage5_network.py ~376-419 pattern).
F5 "sqrt(added channels)" is the wrong expectation (MAJOR for T2): measured 3.3x vs ideal 4.6x at 21 channels; transferred-template
   quality q = 1/sqrt(1+1/(r1² N)) -> 0.51 at N=138, 0.68 at N=342, 0.89 at N=1,500. Criterion = prediction consistency from per-station
   measurements.
F6 "2019-2026 = Bostock-class" is not what F1 says (MAJOR scoping): Puget frac_ge6_40 0.43 (2019), 0.17 (2020), 0.70 (2021), 0.96 (2022),
   1.0 (2023+); S-VI 6-7 within 40 km through 2021, 9 from 2023 (10 UIL stations from 2024). Bostock-class cycles: Puget 2021-02,
   2022-03/04, 2023-04, 2024-08, 2025-10; S-VI 2023+. Sensitivity is set by station-source DISTANCE distribution, not counts.
F7 POLARIS trigger: only if real events register >=2x below prediction AND injection recovers at prediction.
F8 Transfer geometry: families within 40 km — B011 73, B926 52, LRIV 31, STW 50 (EHZ only), NLLB 6, OCP 0 -> drop NLLB/OCP/GOBB. Key by
   net.sta. C8.SOKB 2010 = zero-byte .nd files (no data) -> the ring is 6 stations. SNB is BH 40 Hz (2010) / HH 100 Hz (2017).
F9 Thomas harness: no Thomas PDF on disk; only tidal code is Ducellier's California-only .mat -> tidal leg needs an external
   Cascadia interface tidal-stress model; fronts/re-rupture/rate-decay need no tides.
F10 Do NOT inherit D3's weak stacking as a prediction for T2: D3 stacked Lin picks by 0.1° CELL (several families, ~1 s moveout smear).
   Bostock per-family times have 0.01-s precision; per-family constant offsets are harmless to stacking.
F11 Test 4 has no truth in 2017 -> self-consistency only.

## Prescription
0 Redefine T1-T2 deliverable as a SENSITIVITY MODEL: per (family, station) excess-in-MAD at Bostock times over random-time background,
  per magnitude bin, per template length (4, 6, 8, 12 s), per station with median family distance; fit registration e(M, d, class) and
  the MAD-vs-M slope (measure it). E4 step-1 acceptance becomes: per-station registration of M>=1.8 at PGC/LZB/VGZ/YOUB within a factor 2
  of r1·q·A(M)/sigma_station predicted from stack metadata (q=0.89, sigma = 1/sqrt(2BL)/sqrt(3)). Signatures: <0.5x prediction at all four
  WITH injection at prediction -> real-event variability, POLARIS justified; <0.5x AND injection short -> implementation bug; >=2x ->
  the stack-SNR convention understates r1, rebuild the model from measurements.
1 T2a q-control where truth exists: re-stack each ring station's own 2010 data at Bostock times (per-family offset + Bostock S-P, ±20 s,
  S pick on the stack) on the BUILD half for the 77 families with >=100 events; (i) cc(2010 re-stack, SeisBench stack) >=0.6 for N>=200;
  (ii) held-out registration ratio 2010-template/SeisBench-template = measured q, within ±0.15 of 1/sqrt(1+1/(r1²N)). cc<0.3 for N>=200 ->
  offsets/moveout wrong, stop.
2 HOLD-OUT: alternate 2-hour blocks build/test per family over the whole window with 60-s guards at block edges (secondary: first 60% /
  last 40% by time). Never by event parity.
3 T2b transfer to B011, B926, LRIV (+STW Z-only): stack build-half events at predicted S (±20 s), pick S on the stack, refine moveout by cc;
  6 and 8 s templates. Event counts: N >= 1/r1² for q>=0.7 (~400 at r1 0.05), N >= 1/(3r1²) for q>=0.5 (~130) -> ~35 usable, ~77 marginal
  families. Pass per station: held-out single-station registration >=0.7x the co-located PGC 2010-only re-stack (B011 is 0.3 km from PGC
  = the borehole-vs-surface transfer control); stack S-SNR >=1.5x same-N random null for N>=200. Network: ring+transfer held-out
  registration > ring-only and within 30% of sqrt(sum e_s²)/1.38; shifted-control FA at 8 MAD unchanged within Poisson. Failure: network
  <= ring-only -> moveout refinement fits noise (check by refining on random times).
4 THOMAS HARNESS + DECIMATION (parallel, zero download): on 2012-08-15..10-15 (28,149 events, 130 families; per-family p10/50/90 =
  49/192/392); decimate by Mc 1.7/1.8/1.9/2.0/2.1 (kept 100/33.5/18.2/8.3/3.3%; families >=50 events 130->55->40->13->2); re-run every
  metric -> the completeness each metric NEEDS. Fetch Sci. Adv. adz2832 text for the exact front/tidal criteria before coding.
5 PREDICTION TABLE: per region-year, stations within 60 km with distance + class; r_s from the measured registration-vs-distance curve
  (broadband vs borehole), q from achievable chain-start N, sigma from template length, gain with the 1.38 degradation -> 8-MAD
  completeness magnitude, implied fraction of Bostock density and per-family counts per SSE. Cross with step 4. Sanity: the model on
  Bostock's 21-station array must predict his median M1.73 detectable.
6 T3 DISCOVERY on 2010 (10 3-C stations, 08-19..09-04): (a) DL-seeded — cluster clean LFEDetect events (thr 0.5/N>=4) by differential
  S-arrival pattern (not NLL location), stack per cluster, iterate MF -> re-stack twice; (b) our 12,703 families as seeds only after a
  viability cc test vs SeisBench PGC stacks / transferred B011/B926 stacks (one Bostock family at cc >> others = real seed); (c) blind
  autocorrelation last (single-pair cc ~ r1² ~ 0.003; only M>~2.3 pairs visible, 41 events). Controls: split-half (odd/even active days,
  inter-half stack cc >=0.6), per-station random-shift null, tremor-free window (May-Jul 2010 PGC ~0). Pass: >=50% of discovered families
  match ONE Bostock family by stack-vs-stack cc (best >=2x second) and held-out event-time recall above shifted control; shuffled null
  <=20% as many families. Only then rehearse on 2017.
7 T4 CHAINING 2010->2017 with drift controls: detect 2017-03-13..04-09 with 2010-built templates on persisting stations (PGC, VGZ, SNB
  BH->HH, B011, LRIV); stack CLRS/SYMB/PTRF at those times; leave-one-station-out refinement; refine only from network detections >=8 MAD;
  cc(2017 re-stack, SeisBench) >=0.6; Jaccard(2010-chained vs SeisBench-direct detection lists) >=0.7; shifted-control FA flat across
  generations; 2017 family locations within ~5 km of Bostock's. Failure: cc falls while counts rise -> self-reinforcement.
8 FIRST-DOWNLOAD RULE: (a) POLARIS Aug-Sep 2010 (~10 GB) only under the step-0 signature "real events >=2x below prediction AND injection
  at prediction". (b) Sept 2012 southern VI (CN ring + B011/B926 + LRIV/STW; 10-20 GB native, ~5 GB at 40 Hz) when step 5 predicts
  Mc <= 2.0 on that set AND step 4 shows main+secondary fronts survive at that completeness; cheaper first: ETS-peak days (~20/episode)
  for 2009/2011/2012/2013 at persistent stations (~2-3 GB/episode) to push chain-start templates to N >= 1,000 (q ~0.9). (c) Post-2019 VI:
  2024-08-12..09-01 (largest post-2019 VI episode), stations within 60 km (~15 GB at 40 Hz) only after step 7 passes and (b) reproduces
  the reachable Thomas metrics.

## Missing from the post-2019 pipeline (E)
network-sum thresholding with daily channel normalisation; EQ/blast veto (NRCan/PNSN catalogs + a regular-earthquake picker veto);
response removal for the amplitude module only (2-8 Hz; calibrate on the 2010-2013 Bostock overlap); completeness by per-family per-year
injection of the family's own template; S-only family locations with P from stacks where Z SNR permits + Slab2 prior; Cascadia interface
tidal-stress model; cross-family declustering (assign to max-cc family).

## What would change Merlin's mind
Ring registration >=2x the r1 prediction (or a much steeper MAD-vs-M slope) -> the ring reaches Mc ~1.7-1.8 and the binary framing would
have sufficed. If step 4 shows Thomas's headline metrics need >=50% of Bostock's density, no station set in f1_readiness_40km.csv reaches
it -> replan the science before any download.
