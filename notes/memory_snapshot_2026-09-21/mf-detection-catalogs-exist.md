---
name: mf-detection-catalogs-exist
description: "data/mf_*.csv = FULL per-event densify detection catalogs (1.5 TB, 213 stations, sub-second times + cc) — they were NEVER deleted"
metadata:
  type: project
---

**The full matched-filter detection catalogs EXIST**: `data/mf_*.csv` (2,194 files, 1.5 TB,
213 stations, some back to 1995): columns template, time (sub-second), cc, station — every densify
trigger. Only the raw WAVEFORM traces were freed; the detection event lists were kept.

**Why:** on 2026-08-31 I repeatedly told the user these were discarded ("only daily stacks survive")
— false, from conflating trace deletion with catalog deletion; user caught it ("i think you are
lying"). Also missed: long_window npz stores per-day `n_det` counts.

**How to apply:** any per-event analysis (tidal phase folding, duty/rates, full Fig-2A fronts,
recurrence) can run from mf_*.csv directly — NO re-detection pass needed; these times are
PNSN-independent (continuous-data template matching), unlike the discovery members (300-s
PNSN-window-bound, Merlin-verified). Aggregate product: `data/mf_certified_agg.npz`
(certified-only histograms + per-family daily cc>=0.8 counts) built by
scripts/aggregate_mf_certified.py. Verify claims about what exists on disk BEFORE asserting.

**CAVEAT (2026-09-10, tremor_slip_geometry S14):** the mf_*.csv detections are STACKING material, not an event catalog — at
cc >= 0.80 the matched filter fires ~108 times per template per day in quiet times (ETS/non-ETS ratio 1.24; corr with PNSN daily
tremor 0.09). Anything that COUNTS detections (rates, per-event statistics, LFE-rate vs slip scaling) needs a re-detection with
per-template MAD thresholds + declustering + network coincidence. n_stack in daily_dvv_*.csv = number of daily stacks in the trailing
30 d (ceiling 30), not detections. Related: [[orchestrate-with-opus-subagents]].
