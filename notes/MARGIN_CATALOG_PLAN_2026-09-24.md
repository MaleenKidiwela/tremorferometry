# MARGIN-WIDE LFE CATALOGUE — THE PLAN
Written 2026-09-24 after the user pointed out that this phase was being run reactively.
Supersedes ad-hoc sequencing. Update this file when a decision point resolves; do not let it go stale.

## THE GOAL
A located LFE family catalogue spanning Cascadia, with detections through 2026, good enough to
resolve slip fronts, secondary fronts, re-rupture and healing in EVERY episode — not just the north.
Everything below is instrumental to that.

## WHY IT IS HARD, IN ONE LINE
Dense temporary arrays FIND families; permanent networks only RUN them (measured: permanent beds
discover 3/77 and 5/60 while catalogue templates recover 33/77 and 50/60). Every published Cascadia
catalogue was built on an array that no longer exists. So: reproduce each published set on its own
array-era data, then carry it forward on permanent stations.

## PHASE 1 — PROVE THE ENGINE ON EVERY PUBLISHED CATALOGUE  (the validation phase)
| # | Region | Catalogue | Testable against | Status |
|---|---|---|---|---|
| 1 | Vancouver Is. 48–49N | Bostock 2012/15, 130 fam | times + locations | **DONE** 58/96 blind, 22/30 sealed |
| 2 | S. Cascadia 39.5–41.7N | Ducellier 2022, 66 fam | times + locations + templates | **DONE** recall 0.836, 63/66 |
| 3 | Olympics 47.8–48.3N | Chestler 2017, ~43 fam | 69 locations only | **DONE (spatial)** 61/69 vs null 3/69 |
| 4 | Mendocino 40N | Shelly 2025, 27 fam | times + locations | **DIAGNOSED + RE-RUNNING** (aa_shelly) — the 12.8 M detections were OUR bug, see below |
| 5 | Puget/S-WA 46–47.5N | Royer & Bostock 2014 | **EXISTENCE + 29–46 km band + region** (see below) | RUNNING (ab_royer, pid 1904278 (v2), launched 2026-09-25 09:12Z) — paper-literal rebuild (40 sps, 15 s, their 14 stns, 6 episodes); y_wa_finish died silently 09-24 after beam_2008, superseded |
| 5b | N. Vancouver Is. ~50N | Royer & Bostock 2014 NVI, 54->47 fam | existence + count/depth reference | ADDED 2026-09-25 to ab_royer (POLARIS PO at EarthScope, 2006-09 + 2007-06); was omitted |
| 6 | C. Oregon ~44N | Thomas & Bostock 2015, 18 fam | count + depth + migration rate | QUEUED, not started |
### A published POSITIVE RESULT is itself a testable claim (user, 2026-09-24)
Earlier framing — "Royer & Bostock released nothing, so their Washington families cannot be tested" —
is WRONG and is retracted. They reported that families EXIST in Washington on the interface at 29–46 km.
That is falsifiable on our side without their catalogue: run their data and either we find locatable,
repeating, interface-depth families above the null, or we have failed to reproduce them. Their event
list would only tell us WHICH families we recovered, not WHETHER any are there.
Consequently the two previous Puget zeros (WA_2008_RESULTS.md, WA_2008Q_RESULTS.md) are BUG REPORTS
against our code, not evidence about Washington — and a named candidate bug already exists (the array
treatment that rescued the Olympic bed was never applied to the CAFE beds).
The same logic applies to Thomas & Bostock's 18 Oregon families.
**Scoring rule for any region whose authors released no catalogue:** the test is existence —
do we produce located families, above the whole-chain null, inside the published region and depth band?
Their stated numbers (family count, depth band, migration rate) are REFERENCE LINES, never gates.

**Phase 1 exit test:** every region either reproduced, or its failure diagnosed to a named defect in
our code (never to "not reproducible" — the authors succeeded on this data).

## PHASE 2 — MAKE LOCATIONS CATALOGUE-GRADE  (blocks Phase 3 quality, not Phase 3 start)
The rebuilt locator (L9_xcdiff) is 4.0 km against Bostock's published coordinates in the NORTH but
carries a 35 km median bootstrap radius in WASHINGTON. Until that is understood, WA families are a
population-level match, not a named family list.
- 2a. Diagnose the north/WA gap: station geometry, aperture, depth-grid rails (WA rail 0.53).
- 2b. Reduce it or characterise it honestly per region. Per-family uncertainty is a catalogue column.
- 2c. Decision point: if WA cannot go below ~10 km, the WA product is "families with 20–35 km
      uncertainty", which is still usable for FRONT tracking (the science goal) but not for
      structure. Say so explicitly rather than implying precision we do not have.

## PHASE 3 — BUILD THE CATALOGUE  (per region, same machinery each time)
For each region that clears Phase 1:
  a. templates from each PERMANENT station's own stack at the array-era detection times
  b. frozen extension chain forward to 2026 (FA-only tiers per network epoch)
  c. per-family detection catalogue + completeness by station-set ablation + per-family location+σ
- North: **DONE** (130 families to 2026, lfe_extension v1.1)
- S. Cascadia: next after Shelly scores — Ducellier's 66 are reproduced, so this is mechanical
- Olympics: after Phase 2a tells us what uncertainty to publish
- Puget/S-WA, Oregon: only if Phase 1 clears them

## PHASE 4 — THE SCIENCE THE CATALOGUE IS FOR
Slip-front tracking per episode, margin-wide; secondary fronts; re-rupture; healing; tides.
The parked front-estimator v2 (lfe_extension/results/trackD §3b-2) is the first tool, and it already
has one owed fix: de-duplicate co-detected events across families ≤3 km before fitting a front.


### Royer extension run (user-approved 2026-09-25)
After the literal reproduction (AB_ROYER_DONE), the same queue runs an EXTENSION: literal detections +
templates unchanged; extra public stations (UW/PB/TA/XU/CN within 75/100/150 km, swept) added ONLY at
growth and location (Hyp2000+hypoDD, Table S1). Separate outputs results/ab_royer/ext/, markers
AB_ROYER_EXT_{DONE,FAILED}; null gets identical treatment. Purpose: better depth/locations (feeds D3),
larger families, independent confirmation, and which permanent stations carry each family (Phase 3).


### ROYER MANDATE (user, 2026-09-25): do not stop until the Royer method yields locatable families
Clarified ~19:50Z: NEW locatable families from their method count as success; matching Royer's templates is a second column, not a gate.
Scoring truth: their Fig.4 / Fig.S1 / cross-sections are vector graphics -> digitized template
locations (results/ab_royer/truth/royer2014_templates_digitized.csv, scorer scripts/ab_truth_score.py).
TARGET ("good portion"): in each of NW, SW, NVI — >=50 % of their templates matched one-to-one by a
located family of ours within 5 km epicentral (and far above the rotated/shifted null), depths inside
their band, family count within a factor ~2 of 122/54/47. Loop: run -> score per stage vs their
4915/3306/4267 -> 224/70/54 -> 122/54/47 -> diagnose the first stage that departs as OUR bug -> fix ->
re-run. Log every cycle in notes/YYYY-MM-DD_Notes.md and results/ab_royer/AB_ROYER_LOG.md.

#### Every avenue (user 2026-09-25: "if they did, there's no way we can't") — status tracked here
A1 Literal recipe chain ab_royer (NW/SW/NVI), per-stage counts vs paper — RUNNING (agent AB)
A2 Digitized truth from their vector figures (Fig.4, S1, sections, mech maps) + one-to-one scorer — RUNNING
A3 Literature/thesis hunt: Royer UBC thesis (cIRcle), follow-ups reusing the templates (Royer, Thomas &
   Bostock 2015 tides etc.), data repos, Rowe 2002, Bostock 2012 velocity model, pre-2009 tremor lists — RUNNING
A4 Location-driven positive control: at each digitized Royer location, predict P/S moveouts, beam/stack
   tremor-day data to build a template, matched-filter -> do repeating families exist there on OUR data?
   (separates "our data/stations" from "our detection chain") — after A2
A5 Synthetic injection: repeating synthetic LFEs at Royer depths into real CAFE noise at realistic SNR;
   check each stage (autocorr, cluster, grow, locate) recovers them — unit test of our chain — after AB relaunch
A6 Positive controls already in hand: Sweet 2014 family (10 km from Royer's nearest), our Olympic 2010 run
   and Chestler 61/69 overlap the NW area — compare their families against A2 truth
A7 Stage-by-stage diffs vs paper counts (4915/3306/4267 -> 224/70/54 -> 122/54/47); first departing stage = bug
A8 Parameter readings where paper is silent: swept as columns, never tuned to the counts

## DECISION POINTS (what result changes what)
- **D1. WA mode test (imminent, y_wa_finish stage 2).** If site-averaging lifts the 2008 pair ratio
  from ~0.9 the way it lifted 2010 (1.39→15.81), the two previous Puget zeros were configuration and
  Royer & Bostock's families were never searched for. If it does NOT lift, the Puget bed differs from
  the Olympic bed in a way that must then be named (station quality? aperture? no families there?).
- **D2. Shelly score — RESOLVED as a code fault (2026-09-25, agent AA).** The difference from
  Ducellier was NOT stations or era: (a) `t11_mf.py`'s cc normalisation is unbounded wherever the
  20 Hz cache zero-fills a gap (sd clamped to 1e-15 -> cc up to 2.4e15; 128,922 detections on one
  day against his 71), and Ducellier's path is immune because his cache carries a `valid` mask;
  (b) the network cc was summed without the per-station moveout `tt+tau` (1.49-5.85 s spread inside
  one family), and Ducellier's templates are all cut at the same absolute catalogue time so he needs
  no shift; (c) the recall denominator used all 61,441 of his events against 385 days holding only
  5,021 (cap 0.082). Fixed in COPIES (t11_mf_aa.py etc.); re-run in flight as
  `lfe_south/scripts/aa_shelly_queue.sh`. His published parameters are now in PAPER_METHODS_TABLE
  (2-9 Hz, 6 s, 9xMAD daily, fixed-4-s dedupe, 2-3 template iterations).
- **D3. Phase 2a.** Whether WA locations can reach ~10 km decides whether the WA product is a named
  family catalogue or a front-tracking population.
- **D4. Oregon.** Only after D1 resolves, so we are not running two unproven beds at once.

## WHAT IS RUNNING RIGHT NOW
**2026-09-26: superseded by notes/RESUME_2026-09-26_royer.md (NW station-QC rerun; main queue paused).** Older lines below are historical.

- `lfe_south/scripts/t_south_queue.sh` pid 1509294 — now at the extension stage (do not kill).
  NOTE its extension detections use the unfixed `t11_mf.py`, so §5 of SOUTH_RESULTS is inflated by
  the same normalisation bug and is owed a re-run.
- `lfe_south/scripts/aa_shelly_queue.sh` — the fixed Shelly re-run (rows A/B/C/D + ablation N),
  outputs under `lfe_south/results/aa_shelly/`.
- `lfe_array_discovery/scripts/ab_royer_queue.sh` pid 1904278 (v2) — Royer & Bostock 2014 literal rebuild (D1 now answered by this). Recipe: results/ab_royer/METHOD_FROM_PAPER.md; markers results/AB_ROYER_{HEARTBEAT,AUTOCORR_DONE,DONE,FAILED}. y_wa_finish died 2026-09-24 ~16:30 without marker (superseded).
- Watch armed: notes/queue_watch.log trips if either dies or finishes.

## RULES THAT GOVERN ALL OF IT (see memory/)
Published recipe literally; RANK DO NOT DISCARD (no gates, every criterion a column, null gets the
identical table, evidence is the separation curve); a failed reproduction is our bug; agents launch
detached queues and exit; never pkill -f; document as you go.

## WHAT I WILL NOT DO
Start a new region while a decision point on a running one is unresolved. That is what went wrong
today (Oregon launched mid-Washington, then stopped).
