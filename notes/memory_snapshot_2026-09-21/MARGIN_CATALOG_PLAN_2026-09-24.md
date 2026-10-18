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
| 4 | Mendocino 40N | Shelly 2025, 27 fam | times + locations | RUNNING (t_south) |
| 5 | Puget/S-WA 46–47.5N | Royer & Bostock 2014 | **EXISTENCE + 29–46 km band + region** (see below) | RUNNING (y_wa_finish) — first proper run |
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

## DECISION POINTS (what result changes what)
- **D1. WA mode test (imminent, y_wa_finish stage 2).** If site-averaging lifts the 2008 pair ratio
  from ~0.9 the way it lifted 2010 (1.39→15.81), the two previous Puget zeros were configuration and
  Royer & Bostock's families were never searched for. If it does NOT lift, the Puget bed differs from
  the Olympic bed in a way that must then be named (station quality? aperture? no families there?).
- **D2. Shelly score (imminent).** Second independent confirmation that the chain transfers. If it
  fails where Ducellier succeeded, the difference is Shelly's stations/era and must be diagnosed.
- **D3. Phase 2a.** Whether WA locations can reach ~10 km decides whether the WA product is a named
  family catalogue or a front-tracking population.
- **D4. Oregon.** Only after D1 resolves, so we are not running two unproven beds at once.

## WHAT IS RUNNING RIGHT NOW
- `lfe_south/scripts/t_south_queue.sh` pid 1509294 — Shelly reproduction, then extension.
- `lfe_array_discovery/scripts/y_wa_finish_queue.sh` — 2008+2007 CAFE as arrays, new locator (D1).
- Watch armed: notes/queue_watch.log trips if either dies or finishes.

## RULES THAT GOVERN ALL OF IT (see memory/)
Published recipe literally; RANK DO NOT DISCARD (no gates, every criterion a column, null gets the
identical table, evidence is the separation curve); a failed reproduction is our bug; agents launch
detached queues and exit; never pkill -f; document as you go.

## WHAT I WILL NOT DO
Start a new region while a decision point on a running one is unresolved. That is what went wrong
today (Oregon launched mid-Washington, then stopped).
