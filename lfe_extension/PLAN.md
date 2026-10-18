# lfe_extension — extend the Bostock (2015) LFE family catalog to 2026 on the permanent network

Started 2026-09-16. User decision: extend Bostock's 130 southern-Vancouver-Island families forward from
2013 to 2026, **continuous** (not ETS windows only), on the long-running stations, Bostock's own way
(template matching + template iteration), with our own family locations as a check. Discovery of new
families is NOT in scope here.

Findings that fix the design (see `lfe_catalog_repair/RESULTS.md` E4–E6 and `codex-LFE/docs/`):
- Bostock's 2007–2013 detections used only the permanent ring (PGC, VGZ, SNB, LZB, PFB, YOUB); the 15
  POLARIS/C8 temporaries were gone by 2007. Implied threshold ~3.6 MAD after years of template iteration.
- Our from-scratch 2012 templates (Bostock event times only, 15 stations) reach cc 0.78 to his stacks and
  recover 23 % of his events at 8 MAD (25 % precision) / 15 % at 9 MAD (~50 % precision). False alarms are
  site-coherent noise from the co-located borehole clusters (PGC/B009/B010/B011; LRIV/B005/B006/B007), not
  random: one-station-per-site brings quiet-day FA to the rolled control. E6 site-merge A halves FA at
  8 MAD (0.43/template/day) at a recall cost (0.19). No template iteration has ever been run.
- The dv/v families are NOT Bostock families (waveform identity test 2026-09-16, scratchpad); they cannot
  seed anything. Seeds = Bostock's 130 families (times + SeisBench stacks P/S).
- Codex audit corrections adopted: one-to-one matching everywhere, rotate boreholes to ZNE from
  StationXML azimuths, 35-day divisor bug, roster = 183 stations of figures/borehole_dvv_map.html.

## Stations (permanent, over or near the families)
Roster: CN PGC, VGZ, SNB (→2021-03); PB B001, B003, B004, B005, B006, B007, B009, B010, B011, B926.
Off-roster permanents that carried Bostock/our templates (add, they cost nothing): CN NLLB, LZB (→2014-06),
YOUB (→2014-09), PFB (→2018-12); UW LRIV; C8 SOKB (→2015-01).
Channel band per station per epoch from StationXML (PGC BH→HH change; SNB BH→HH 2011-03-17; PB EH? 100 Hz).

## Data
On disk (`lfe_catalog_repair/traces`): 2010 ETS (Aug 10–Sep 13, 13 stns), 2012 ETS (Aug 15–Oct 5, 15 stns),
2017 ETS (Feb 20–Apr 21, 10 stns), PGC May–Nov 2010.
To download (Track B): 2011 ETS (2011-08-05..09-25) and 2013 ETS (2013-08-30..10-28) for all stations
(overlap years → template iteration + held-out validation); then **continuous 2014-01-01 → 2026-09-15**,
oldest first, month by month, all stations in parallel. ~1.4 TB at 20 MB/station-day; disk has 74 TB.

## Tracks (run concurrently — never wait on the download)
- **Track B (download, background, starts now):** `scripts/download_plan.py` builds the per-station
  epoch/band table (reuse codex-LFE metadata, read-only) and the ordered queue; `scripts/download_run.sh`
  streams it with ≤8 concurrent FDSN requests (IRIS for CN/PB/UW/C8), resumable, `.nd` no-data markers,
  `results/download_manifest.csv` (station, date, bytes, status) updated per day-file so Track C can follow.
- **Track A (method, on cached 2010+2012 data, GPU):**
  A1 site-level merge before the network sum (one 3-C series or site-mean cc per site), site-rolled chance
     control (roll SITES not stations). Operating point chosen on 2012 quiet days + Bostock recall, then frozen.
  A2 template iteration, Bostock style: detect at high confidence (≥9 MAD, ≥2 sites), re-stack per station
     (SNR-weighted, sinc-aligned), re-run; ≤3 iterations; guard against drift with a fixed q-control
     (cc to the SeisBench stack must not fall) and against self-reinforcement by validating on the
     HELD-OUT overlap year (templates from 2012 → test 2010, and vice versa; later 2011/2013 when down).
  A3 family location from OUR stacks: pick P (Z) and S (H) on the iterated per-station stacks, locate with
     Crosson 1-D (+ Slab2 prior as a separate run), compare with Bostock's hypocentres (bar: ≤ 5 km median).
  A4 per-family completeness model: recall vs Bostock by family and by network epoch (which stations alive).
- **Track C (production, follows the manifest):** per (month, station-set) block once all stations for
  that month are on disk: preprocess (40 Hz, 3-C, ZNE, 2–8 Hz), network matched filter with the frozen
  templates of the current epoch, detections + daily per-station stacks saved; re-stack templates at
  each instrument change / lost station, with the q-control; keep raw data (no deletion needed).
- **Track D (product):** events (time, family, network MAD, rel. amplitude, n_sites), family table
  (Bostock hypocentre + OUR hypocentre + uncertainty), per-family per-epoch completeness, Thomas harness
  run per SSE (`lfe_catalog_repair/scripts/e5_thomas_harness.py`) as the science QC.

## Bars (pre-registered; report met / not met, never tune to them)
- FA: quiet-day (no PNSN tremor in box) detections ≤ 1 / template / day AND ≤ 5 % of ETS-day count, with
  the SITE-rolled control; earthquake veto ≤ 0.5 % of detections within 60 s of a catalogued EQ.
- Recall vs Bostock on the held-out overlap year (own family, ±2 s, one-to-one): report; target ≥ 30 % at
  the FA-compliant point (E5 measured 15 % at 9 MAD without iteration/merge).
- Location: our family hypocentres within 5 km (median) of Bostock's.
- Science QC: main-front speed within 20 % of Bostock's on 2012 (already met at 5 %); secondary-front
  count and re-rupture fraction reported vs his catalog decimated to our recall.

## Governance
Merlin reviews this plan before Track A starts (download is method-agnostic and starts immediately).
Opus agents code; main session specs/reviews. Everything in this folder: README / PLAN / NOTES (dated
ledger) / scripts / results / figures / logs. Read the parent folders read-only; copy, don't import.

## Amendments after Merlin review (2026-09-16; full text `results/merlin_review_plan_2026-09-16.md`)
- PRODUCT = TWO-TIER catalog: T_core (quiet-day FA ≤ 0.1 / template / day) everywhere; T_ext (quiet FA ≤ 1.0,
  post-SSE FA ≤ 1.5) only inside tremor windows or with its false-alarm floor λ attached to every event.
- Thresholds set by FALSE-ALARM criteria ONLY; recall is reported at the resulting T, never used to pick T.
- Site merge = E6 variant B (within-site weighted-mean cc, site weight √Σw²), settled by the matched-FA ROC table.
- Honest recall expectation 0.18–0.24 at FA 1.0 on 15 stations (2012), lower post-2021; iteration is for
  TRANSFER (epochs, template-less families), expected ≤ +0.02 recall.
- Bars replaced: ETS-day site-rolled ≤5 % (unfalsifiable) → quiet-day site-rolled/data ∈ [0.77, 1.3] + reported
  one-to-one precision; recall ≥30 % → within 0.7–1.3× the station-set ablation prediction; location 5 km →
  horizontal median ≤5 km with shuffled-pick null and synthetic floor, depth reported.
- New bars: frozen cross-year transfer (2012→2010, 2010→2012); epoch identity |cc| ≥ 0.6; completeness model
  validated on held-out 2010 within ±0.05; winter quiet-day FA within 2× summer.
- Measurement chain frozen on the validated E4–E6 path (20 Hz, 1–8 Hz, existing component handling); ≥2-site
  observability rule; one-to-one matching.
- Track B queue re-ordered: one SSE window per instrument epoch first (2015-12-05..2016-01-15, 2018-06-10..07-25,
  2021-10-01..11-05) + one winter quiet month (2016-12), then oldest-first fill.
- (user, 2026-09-16) Quiet-day detections are an UPPER BOUND on false alarms, not a measurement: families fire between
  ETS episodes. Step 3b (required before Track C): withheld-station stack test on quiet-day detections per family
  (stack at a station excluded from detection; S arrival at the predicted moveout vs a random-time stack). Families that
  pass get their inter-ETS detections promoted to events with their own completeness/floor; the rest stay as λ.
