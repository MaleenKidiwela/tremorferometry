# Method selection: bounded comparison before scaling

Latest experiment: [training-derived physical-site weighting](SITE_WEIGHTING_EXPERIMENT.md)
did not improve calibration recovery under the existing control gate (29 versus
31 matched references). Equal weighting remains the baseline. This rules out one
proposed improvement on this panel, not the feasibility of an LFE catalog.

Current follow-up: [completed miss diagnosis](MISS_DIAGNOSIS.md) points to limited
four-site pick support and low template-match scores. Prioritize a controlled
sensitivity experiment and waveform review before further association tuning or
GPU scaling. The method choice remains open.

The decision is which detection method to carry forward. Bulk waveform downloads,
margin-wide picking, and a production catalog remain deferred. The current CPU
experiment compares roster-restricted template matching with corrected neural
picking and association; it does not commit the project to neural association.

September 12 update: the [revised interleaved comparison](INTERLEAVED_METHOD_RESULTS.md)
is complete. It includes parallel CPU inference with verified borehole orientation,
new observed/control association, and shared-reference scoring with explicit
prediction coverage and model-dependent training clocks. Keep matched filtering
as the known-repeat baseline; diagnose missed events and review candidate identity
before choosing a full-margin method. The sections below preserve the initial
reuse experiment and its design history.

## Existing evidence and the missing comparison

Claude's corrected injections and continuous model runs establish useful starting
points. His best association experiment used ten stations, only five of which
belong to the user's map roster. Our previous matched-filter experiments used
2012 dates and different station subsets. Comparing their headline counts or
recalls to the 2010 association run would not rank the methods fairly.

The current northern development panel is August 10–September 13, 2010, with
CN.PGC, CN.SNB, CN.VGZ, PB.B011, and PB.B926. All five have saved Cascadia-model
picks. Use this already-examined panel to debug and compare implementations;
do not call it an untouched test.

## First association comparison

`codex_lfe.reuse_association` consumes saved picks, station coordinates/elevations,
and the existing Crosson travel-time table. It restricts station identity before
association and discards peaks within one second of UTC day edges. The reported
event interval excludes 3,120 seconds at each end for observed and shifted runs.
This bounds outer-edge effects of shifts up to 3,000 seconds plus the 120-second
association overlap; internal gaps and detailed prediction coverage remain
limitations. Input `t_start`/`t_end` columns describe individual pick excursions,
not continuous inference coverage.

Predeclared development grid: thresholds 0.3 and 0.5; at least four or five
distinct S-picking stations, with optional P picks; three constant-station and
three station-hour shift realizations per setting. Seeds are 1, 2, 3. Save all
settings, including zero-event results. No threshold is selected on counts alone.

Retain the legacy spatial search and 1.2-second association tolerance to isolate
roster feasibility. The depth search is 20–50 km, so apparent slab-depth locations
are not independent validation. The passed `edt_pick_std=0.5` parameter is
inactive because exponential EDT is disabled, as in the original runner. This
was confirmed in the installed PyOcto source and [official documentation](https://pyocto.readthedocs.io/en/stable/pages/associator.html).
Consequently the earlier improvement cannot be attributed to changing that
standard deviation. PyOcto locations remain preliminary.

### Physical-site correction

The first run exposed a limitation of station-wise controls: CN.PGC and PB.B011
are approximately 290 m apart. They must not be treated as independent noise
sites. The second development run groups stations within 1 km, shifts the two
nearby instruments together in both control modes, and requires S support at
all four distinct physical-site groups. One kilometre is an explicit conservative
grouping rule, not proof that more distant stations have independent noise.
This change was motivated by geometry after the first run; v2 is another
development diagnostic, not a preregistered held-out result.

The distinct-site rule is a post-association filter and can miss alternative
associations that a rerun on an individually selected instrument subset might
find. Saved pick identities and assignments make this inspectable. Counts from
the original station-wise controls remain archived and must not be substituted
for the grouped controls when making a method decision.

The local miniSEED header audit found 171 of 175 possible station-day files:
PGC has all 35 days, the other four lack September 13. Their header streams
contain 89 gap/overlap records in total. All present files have three channel
names, but header inspection does not establish sample validity, orientation,
or continuous inference coverage. No raw waveform downloads were needed for
these association experiments.

Reproduce the second run from `codex-LFE`:

```bash
python -m codex_lfe.reuse_association \
  --output results/map_network_association_sites_development_v2 \
  --colocated-km 1 --minimum-sites 4 --threads 8
```

The output directory contains input/source hashes, a copy of the runner, retained
input picks, grouped control offsets, event/assignment tables, individual run
counts, and the comparison. Completed runs are reused only with matching
configuration/source hashes and verified event/assignment output hashes.

## Arms to compare on shared data

| Arm | Source of candidate events | What would justify using it |
| --- | --- | --- |
| Known-family matched filtering | Verified training-only published-family templates | Good recovery of repeats at acceptable measured contamination; this arm alone cannot establish new-family discovery. |
| LFE neural picking plus physical association | Published-model P/S proposals, associated before review | Independently supported events/families beyond the template inventory, useful recovery, tolerable association errors. |
| Neural discovery plus template extension | Verified neural seeds produce training-only templates | Additional independently supported weak repeats at comparable contamination, sufficient to justify the added complexity. |

Use identical station sets, available waveform intervals, and final-event scoring.
Keep source discovery distinct from known-family recovery. Reconcile event-origin
and template-arrival clocks on development data before scoring. A broader origin
search and a different regional velocity model are sensitivity checks, not a way
to select whichever map looks most plausible.

## Evidence required for a choice

- One-to-one reference agreement with source/timing uncertainty and exposure
  exclusions reported. Reference nonmatches are not automatically false events.
- Multistation waveform review with uncertain cases retained; separate phase
  timing, ordinary-earthquake confusion, duplicates, and wrong associations.
- Withheld allowed-site evidence where enough independent sites remain. Requiring
  all five northern sites leaves none for that check; a four-site discovery run
  must withhold the fifth before association, not remove it retrospectively.
- Both shifted controls through the full pipeline, with matched time support,
  counts and variability. Control ratios are not proven false-discovery rates.
- Whole-date and source-family holdouts for the final decision benchmark, plus
  pretrained-weight exposure checks. Existing development results cannot be
  relabeled independent by retrospectively splitting them.
- A modest central/southern transfer panel and a borehole/surface comparison
  before choosing a margin-wide method. Download only the missing benchmark
  waveforms after dates and stations are specified; bulk acquisition comes later.

No numerical precision/completeness requirement has yet been agreed for the
scientific catalog. Report tradeoff curves rather than inventing an acceptance
bar or declaring a winner from an arbitrary control-ratio cutoff. Choose the
simplest method supported by the shared benchmark. The hybrid is a hypothesis,
not the predetermined outcome.
