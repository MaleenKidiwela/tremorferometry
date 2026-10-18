# SUPERSEDED (2026-09-30): the Ducellier2022 files in this folder

`Ducellier2022_<month>.parquet`, `*_coverage.csv`, `plan_Ducellier2022_*`, and in `results/`:
`catalog_Ducellier2022.parquet`, `catalog_Ducellier2022_per_family_year.csv`, `tiers_Ducellier2022.csv`
were made by `t16_extend.sh` with the UNFIXED `t11_mf.py` (zero-filled gaps -> unbounded normalised cc;
all channels summed at a common lag) on OUR OWN stacked templates at the 12 nearest modern stations
(193,445,746 raw detections, 27.2 M "core").  They are kept on disk but must not be used.

Replacement: `results/DUC_EXT_FIXED_RESULTS.md` and `results/duc_fix/` (queue
`scripts/duc_ext_fixed_queue.sh`): his templates, his stations (re-templated where the instrument
changed), his hourly 8 x MAD and cc x nchannel rule, the sh_v2 fixes.
