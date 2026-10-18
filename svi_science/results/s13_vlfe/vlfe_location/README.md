# vlfe_location — how well are the catalogue VLFEs located?

User 2026-10-16: "these [VLFEs] are not located right". This folder checks the location of the 2010–2018 catalogue
VLFEs (`../catalogue/`, read only) against the LFEs that fire at the same time, measures the resolution the scan grid
actually achieves, tests whether locations copy the detecting template, and tries cheaper better Green's functions.

- `PLAN.md` — pre-registered tests 1–4 (2026-10-16 ~10:02 am PT).
- `NOTES.md` — work log + amendment A1.
- `RESULTS.md` — plain answers with numbers.
- `scripts/` — t1_lfe_centroids.py (Test 1), t2_zsurface.py YEAR (Test 2 + data windows), t23_summary.py (Tests 2–3 +
  figures), t4_relocate.py (Test 4).
- `tables/`, `figures/`, `data/` (z surfaces, data windows, template banks).
