# lfe_south — a LOCATED LFE catalogue south of Bostock's field

Agent T, started 2026-09-23. Goal (user): extend the located-LFE margin SOUTH, by the route that
is already proven — permanent networks RUN families they are handed, they do not FIND them
(north: 3/77 in 2012 and 5/60 in 2010 discovered, vs 33/77 and 50/60 from catalogue-built
templates). Every southern family below is already found and located by someone else.

## The two southern truths (details: results/SOUTH_TRUTH_INVENTORY.md)
* **Ducellier & Creager (2022)** — 66 families, 39.45–41.72 N, 15–46.8 km deep; 477,715 detections
  2004–2011 on the PERMANENT networks plus 107,924 on the FAME array 2007–2009. **Templates,
  station lists and every parameter are on disk.** 65 % of his permanent bed is vertical-only.
* **Shelly et al. (2025)** — 27 families in one tight zone at 40.03–40.09 N, 22–29 km deep;
  61,441 detections 2018-01 .. 2024-11. Times and locations only, no templates.

## Stages (the queue `scripts/t_south_queue.sh` runs them in this order)
1. **Existence test** (Merlin's o4 spec) — stack the modern permanent stations at each family's
   published times, tau searched ±30 s, against 200 random-time nulls. Reported as a
   DISTRIBUTION, ranked; nothing filtered.
2. **Reproduction of Ducellier** — HIS templates, HIS stations (bed A = his exact set including
   every vertical-only station; bed B = his set plus stations installed since), HIS parameters
   (1.5–9 Hz, 20 Hz, response-removed, per-channel cc averaged over complete blocks, 8×MAD,
   1-s decluster, 12-h blocks), scored one-to-one at ±2 s against his own catalogue with three
   denominators (all / 60-s-independent / loud decile) and the precision that Bostock's 30 %
   criterion reads.
3. **Templates from our own stacks**, per station, build half / test half.
4. **Reproduction of Shelly** — our templates from HIS times, scored on days the template never
   saw.
5. **Blind discovery** on the southern permanent bed with the reproduced Bostock recipe
   (ks_*.py, ported unchanged but for the 1-C/2-C loader fidelity fix) — the southern analogue of
   58/96, with the section-2 reproduction as its ceiling.
6. **Extension to 2026** — continuous record, matched filter, FA-ONLY tiers per network epoch by
   the frozen c_tiers rule, per-family per-year counts, station-set ablation.
7. `results/SOUTH_RESULTS.md` — reproduction first, then existence, then extension.

## Binding rules
* Published recipes, unedited. The only edits are fidelity edits forced by the southern network
  (accepting 1-C/2-C stations, rational resampling), each recorded in NOTES.md.
* **RANK, DO NOT DISCARD** (user 2026-09-23). No family is ever dropped for failing a clause.
  Every criterion is a COLUMN with its measured value. Merlin's "≥2 stations above the null p99"
  and "≥10 of 27" and Bostock's "≥30 %" are REPORTED REFERENCE LINES on a distribution.
  Every family with any station above its own null median is carried forward; the extension is
  the test.
* Launch-and-exit: one detached flocked queue, per-stage guards, TERM/INT/HUP traps.
