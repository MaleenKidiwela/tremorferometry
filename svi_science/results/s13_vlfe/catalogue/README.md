# catalogue — 2003–2026 strong-VLFE catalogue, southern Vancouver Island

Goal: a catalogue of strong very-low-frequency earthquakes (20–50 s, Mw ~3–4) on the SVI plate interface, 2003–2026,
built as Merlin designed it (`MERLIN_CAT_DESIGN.md`): empirical-template matched filter as detector, the R3 per-event
battery (far/near ratio, plane wave, location, earthquake windows, single station) as vetter, thresholds from noise
templates per template-year, scrambled controls everywhere. Rules: `PLAN_CAT.md` (pre-registered). Log: `NOTES.md`.
Answers: `RESULTS.md`. Parent work: `../hg_repro/` (Hutchison & Ghosh 2019 reproduction, R2/R3), `../` (S13 scan).

Stage 1 (2026-10-15, data on disk 2011–12): Merlin steps 1–4 and the metadata half of step 5.

Scripts (`scripts/`):
- `cat_common.py`: paths, pass rules a–f, e′, S5, scan-node table (location, footprint, family distance), head2 magnitudes.
- `cat_eprime.py`: e′ validation on the R3 sets → `tables/eprime_*`, `figures/eprime_ratio_vs_z.png`.
- `cat_e2.py` (+ `e2_queue.sh`, marker `E2_DONE`): battery on every V4 T1–T3 detection ≥ 5.5 MAD, the 16 scrambled
  controls and 10 noise templates → `data/e2/`, `tables/e2_*`, `figures/e2_*`.
- `cat_cont.py`: continuous 0.02–0.05 Hz velocity for CN.PGC/SNB/PFB/YOUB in the hg_repro frame → `data/cont/`.
- `cat_e3.py` (+ `e3_queue.sh`, marker `E3_DONE`): empirical template bank, matched filter on four station sets with
  controls and noise templates, thresholds, leave-one-out re-find, battery → `data/e3/`, `tables/e3_*`, `figures/e3_*`.
- `cat_farref_meta.py` (marker `FARREF_META_DONE`): FDSN station metadata + MUSTANG daily availability for far LHZ
  references 250–800 km, CN broadband class → `tables/farref_*`. No waveforms.
