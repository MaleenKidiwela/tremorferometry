# NW dominance-cap repeat on two more episodes: rules fixed BEFORE any run (2026-10-03 01:43Z)

Question: should NW production cap any site's share of a family's site weight at 0.3?
On V2010 (results/wa_pre/cap/) the cap raised NW own precision 0.318 -> 0.365 and NW family-windows above both
scrambled controls 81 -> 111, while own recall fell 0.124 -> 0.107. That was one 39-day episode. This test repeats it
on two other episodes, NW (Royer_NW_WA, 132 families) only.

## Episodes
- **V2008** (2008-05-01 .. 05-26 episode days + the 3 "post" days 2008-06-16 .. 18, i.e. the whole V2008 window as in
  the go/no-go V2008 supplement). Scored: **test half only** (templates come from the build half); HF cuts calibrated on
  the V2008 build half (`gonogo_score.py --wins V2008`, exactly as `results/wa_gonogo/score_core21_v2008`).
- **A2011** (2011-07-30 .. 08-27). Scored: **all of A2011** (an episode the templates never saw). Scored together with
  V2010 (`gonogo_score.py --wins V2010,A2011`, exactly as the go/no-go `score_core21`), so HF cuts come from the V2010
  build half and A2011 uses V2010's per-epoch tiers (not used in the headline). The V2010 rows of these runs are a
  reproducibility check, not new evidence.
- Not used: 2007-01 and 2009-05 (sealed one-look), 2011-04 (excluded).

## Runs (all on the same 21-core network, `results/wa_gonogo/setup`, 20 usable stations / 13 sites; PB.B926 has no
build-half template and contributes nothing)
- Matched filter `scripts/pre_mf_cap.py` (through `scripts/pre_run.py`, download roots cut out), unchanged, with the
  V2010 cap run's arguments: `--tpl build --all-templates --drop-dead 0.05 --n-control 2` (same seeds -> same scrambled
  control shifts per day), weight modes:
  - **uncapped** (`--weight-mode frozen`; must reproduce `results/wa_gonogo/mf_core21` / `mf_core21_v2008`),
  - **capped** (`--weight-mode cap --cap 0.3`, water-filling),
  - **dominant-site-out** (`--weight-mode lodo`), for reference.
- HF column `wa_hf.py --tpl build --all-templates` on each forward run; scorer `gonogo_score.py --tplmode all`
  unchanged.
- **Time-reversed-template control** (reference only, not part of the rule): the same three weight modes with
  `--tpl buildrev --no-control` (every template channel reversed, same offsets and weights; templates
  `results/wa_pre/ctl/res/templates/tpl_buildrev.npz`). Reported per variant: >= 8.75 x MAD candidates per observable
  scanned episode family-day, own-family "precision" of those candidates on the scored half (its chance level), and
  the number of family-windows where data >= 2x the reversed control on label days.

## Metrics (NW, from `summary.csv` of each score, row window x Royer_NW_WA; frozen core tier 8.75 x MAD)
- NW own precision = `prec_own_cat_core` (catalogue detections matched one-to-one +-2 s to own-family truth).
- NW family-windows above both scrambled controls = `sepboth2_core_data`: family-windows with
  (data + 1) / (max(control 1, control 2) + 1) >= 2 on the family's label days (raw >= 8.75 counts). Control-vs-control
  (`sep2_core_ctl`) is reported as the false-pass level.
- NW own recall = `recall_own_cat_core`.

## Acceptance rule (fixed by the main session)
Adopt the 0.3 cap for **NW production only** if, on **BOTH** V2008 and A2011, capped vs uncapped:
1. NW own precision rises by >= 0.02 (absolute: prec_cap - prec_uncapped >= 0.020), AND
2. NW family-windows above both scrambled controls do not fall (sepboth2_cap >= sepboth2_uncapped), AND
3. NW own recall falls by no more than 20 % relative (recall_cap >= 0.80 x recall_uncapped).
Otherwise NW stays uncapped. Dominant-site-out and the time-reversed control do not enter the decision. No other
criterion, threshold or episode is added after seeing results.

## Data completeness (a stage that fails stops the test; no verdict is given from it)
- Before the runs: on-disk core station-days on the V2008 and A2011 days (raw in a non-download root or cache) whose
  20 Hz cache is missing or unreadable must be <= 3 %. Counted on the 20 stations that can contribute; the 21-station
  count (with PB.B926, raw on disk but never cached and no template) is reported beside it. Precedent: the go/no-go V2008
  supplement, 27 / 539 = 5.0 % on 21 stations, 26 of them PB.B926; 1 / 513 = 0.19 % on the 20.
- After each matched-filter run: cached core station-days not loaded for any reason other than the frozen > 5 %-guarded
  drop must be <= 3 %.
- HF stage: station-day jobs that return no ratios must be <= 3 % of jobs.
- Reproducibility checks (stop and investigate if they fail): the uncapped reruns equal the go/no-go matched-filter
  output (every candidate, net_mad); the capped / dominant-site-out V2010 rows reproduce `results/wa_pre/cap/cap_comparison.csv`
  (NW own precision 0.3653 / 0.6765, family-windows 111 / 128).
