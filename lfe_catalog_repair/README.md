# lfe_catalog_repair — can the densify matched-filter detections be turned into an LFE event catalog?

Background (tremor_slip_geometry S14, 2026-09-10): data/mf_<sta>p90f40_<year>.csv (template,time,cc,station; 1.5 TB, 2194
station-year files) are the densify detections at a FIXED threshold cc >= 0.80 (6-s gap). They fire ~108 times per template per
day at rest; the ETS/non-ETS ratio is only ~1.24 and the correlation with PNSN daily tremor 0.09 -> as stored they are stacking
material, not an event list. Real LFEs are a ~20 % excess on a constant false-alarm floor.
Question D1: do the real detections (the ETS excess) sit at distinctly higher cc than the floor? If yes -> per-template
threshold + declustering + network coincidence recovers a catalog from the stored hits (no waveforms). If no -> re-detection from
waveforms with per-template noise thresholds is required.
Folder: PLAN.md, NOTES.md (ledger), scripts/, results/, figures/. Inputs read-only: data/mf_*.csv, catalogs/pnsn_tremor_cascadia_full.csv,
tremor_slip_geometry/results/s1_events.csv, tremor_slip_geometry/results/s14_coverage.csv (which station-years were streamed).
