# P0-O item 4 — southern truth-existence test: DEFERRED (2026-09-17)

**Status: not run, nothing queued.** The test cannot be run today and nothing was downloaded for it,
by instruction. This note records exactly what it needs so it can be started the moment the southern
data lands, and what must NOT be claimed until it passes.

## Why it is deferred
No 2019+ southern station-day exists anywhere on disk. The only trace trees in the project are
`../lfe_extension/traces` (19 stations, all southern-Vancouver-Island / Olympic: C8.SOKB, CN.LZB,
CN.NLLB, CN.PFB, CN.PGC, CN.SNB, CN.VGZ, CN.YOUB, PB.B001/B003/B004/B005/B006/B007/B009/B010/B011/B926,
UW.LRIV) and `../lfe_catalog_repair/traces` (the same plus CN.CLRS, CN.GOBB, CN.PTRF, CN.SYMB, PB.B943,
UW.OCP, UW.STW, UW.WISH, UW.YACT). **None of B045, B935, RVIT or PETY is present**, and the picker's own
`traces/` holds only the 2008/2009 borehole windows fetched for P0-O item 1.

## What the test is (Merlin, P0-O step 4)
For each of David Shelly's 27 Mendocino-triple-junction LFE families, stack the roster stations at his
detection times with a per-family time offset searched over ±30 s, and compare against the same search
on 200 random-time null stacks.

* **Pass** = at least 10 of the 27 families show an S-window energy excess above the null p99 at ≥ 2
  stations.
* **Fail** ⇒ P4-south is reduced to *injection sensitivity + activity series only*, and **no southern
  recall number may be quoted**, because there would be no denominator: the families would not be
  known to register on our stations at all.

## What it needs
| input | where | status |
|---|---|---|
| Shelly 2025 detection times, 27 families | `../lfe_catalog_repair/results/external/other/LFEcat_MTJ_2018-2024_shelly2025.csv` (61,441 rows, `ID`/`idnum` = family, 2018-01-01 .. 2024-12-31; time from `%year,month,day,hour,minute,second`) | **on disk** |
| Shelly metadata | `../lfe_catalog_repair/results/external/other/Shelly2025_metadata.xml` | on disk |
| stations | **PB.B045, PB.B935, UW.RVIT (or NC), NC.PETY** — plus, if they can be had, PB.B049 / PB.B932 / PB.B933 as the nearer Mendocino boreholes already scouted by the dv/v work | **NOT on disk** |
| continuous data | the SSE windows of 2019-2024 that overlap Shelly's detections, ~3-C at native rate, then the frozen 20 Hz / 1-8 Hz chain | **to download** |
| station metadata | FDSN `level=channel` for the four stations, 2018-2026, band preference HH > BH > EH | to fetch |

Estimated volume: 4 stations × the Shelly-active windows. A single 60-day window is ~5 GB; the whole
2018-2024 overlap at 4 stations is ~50-60 GB — trivial against the 72 TB free, so the constraint is
only that the download has not been authorised, not disk.

## How to run it when the data lands
1. `scripts/o_download_plan.py` with `ROSTER = [("PB","B045"),("PB","B935"),("UW","RVIT"),("NC","PETY")]`,
   `WINDOW0/WINDOW1 = 2018-01-01 / 2025-01-01` and one priority per Shelly-active window; then
   `scripts/o_download_worker.py --workers 8`. Output stays in `lfe_picker/traces` and
   `lfe_picker/results/download_manifest.csv` (never `../lfe_extension`).
2. `scripts/o_cache_month.py` for each month (the frozen chain, unchanged).
3. A new `scripts/o4_shelly_stack.py`: for each family, for each station, stack 3-C windows at
   `t_detection + tau` over a ±30 s grid of `tau` (0.05 s steps is enough at 20 Hz / 1-8 Hz), take the
   statistic `10*log10( E(S-1..S+5 s) / E(S-12..S-6 s) )` used by Gate 1, and repeat the identical
   search on 200 random-time stacks of the same size. Report per (family, station): best tau, excess
   in dB, null p99, pass/fail; then the family-level rule above.

## Until it passes
* No southern recall, completeness or Mc number may be quoted for Shelly's families.
* P4-south is "injection sensitivity with the north-frozen mapping, labelled **uncalibrated**".
* The same existence test is owed for Ducellier's 2004-2011 families before any claim there; Merlin's
  F14 notes that his templates live on FAME temporaries plus KHBB/YBH/WDC/KRMB/KSXB/B039, none of which
  is on our roster either.
