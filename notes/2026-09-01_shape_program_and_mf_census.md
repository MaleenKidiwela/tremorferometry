# 2026-09-01 — MF catalog census + waveform-shape program (design→pilots→maps→audit)

## MAJOR DATA-INVENTORY CORRECTION (memory: mf-detection-catalogs-exist.md)
`data/mf_*.csv` = FULL per-event densify detections (1.5 TB, 2,194 files, 213 stations, sub-second
times + cc) — NEVER deleted; I had falsely claimed only day-stacks survived (user caught it).
Also: long_window npz carry per-day `n_det` counts. Consequences: per-event science (tides, rates,
recurrence, full Fig-2A) needs NO re-detection; Merlin's tidal-salvage path is executable from disk.

## CENSUS (scripts/aggregate_mf_certified.py -> data/mf_certified_agg.npz)
Certified-only (minus 229-family blocklist, date vetoes): 12,112 (station|family) series,
**15.2 BILLION detections at cc>=0.8** (3.46B at cc>=0.85), zero file errors. ~200 triggers/day/family
=> trigger-floor dominated; the full-record Fig-2A (figures/lfe_fronts_full_mf.png) is a featureless
wash — EMPIRICAL CONFIRMATION of Merlin's threshold warning; fake-template duty calibration now
non-optional. Rate-ANOMALY view (figures/lfe_rate_anomaly.png) shows structure but the blue->red
2019-21 step is ROSTER GROWTH (strips sum stations) — per-(station,family) normalization needed
(daily_cc80 matrix in the npz supports it; NOT yet built).

## SHAPE PROGRAM (user-driven: "what does waveform shape vary across space tell us")
- Merlin DESIGN ruling: A3 (coda/direct) + A1 (early-coda decay) well-posed WITHOUT response deconv
  (ratios/log-slopes cancel station response); A2 absolute duration DEAD on this archive (2-8 Hz
  prefilter destroyed sub-2Hz corners) -> A2r = co-located in-band spectral contrasts; decay = 
  "effective early-coda decay" (lapse/tS ~1-1.6, NOT intrinsic Qc), fit +4..+9 s, mirror-floor
  truncation (3x), NO Wiener, horizontals primary; everything through two-way FE design.
- PILOTS (scripts/pilot_shape_metrics.py): reliability decay +0.90 / cdr +1.00 / cent +1.00 (PASS);
  variance: station FE 0.59-0.78, residual 1.9-9x above split-half noise; co-located centroid diffs
  36-50x floor. GRADUATED per Merlin's pre-registered bar.
- STAGE 2 (scripts/shape_stage2_maps.py -> res_catalog/shape_maps.npz, figures/shape_maps.png):
  midpoint-cell FE maps: cdr R2 0.82 checker-recovery 0.78 (coherent: central coda-rich vs south
  direct-dominated); decay 0.57 (sparse, n=1873); A2r spread 0.26 Hz everywhere.
- Merlin RESULTS audit: **C1 "medium map" FAILS as stated** (two live confounds: no source FE — and
  C3 proves source variation large; cdr direct-window = template segment -> per-(station,family)
  selection bias, structural); survives only as "organized structure exists". C2 parked (noise-
  geography selection). **C3 (co-located source differences) = the presentable finding** pending
  checks. Upgrade path: (1) within-family cross-station discriminator, (2) selection-immune
  coda-only ratio E[+6.5..+9]/E[+2.5..+5], (3) stage-3 source-FE solve + 2-D checker + lambda
  oracle, (4) C3 checks, (5) C2 later.
- STEP-1 DISCRIMINATOR RUN (342 contrasts): slope +1.35 vs map prediction (medium~+1, source~0),
  corr +0.25, permutation p=0.0005 -> PARTIAL PROMOTION of C1-as-medium (path-dependence real but
  loose). C3 check (i): corr(|diff|, stack days) = -0.11, robust at high-n (0.30 Hz median) — PASS.
  C3 check (ii): ordering replication 28/41 across stations — middling.
- REMAINING before findings: step-2 selection-immune metric map (THE gate for all cdr products);
  stage-3 source-FE solve; then Merlin re-ruling. NOTHING here is a presentable finding yet except
  C3 is close.

## Standing figure inventory from today
lfe_fronts_full_mf.png (trigger-floor wash), lfe_fronts_cc85.png, lfe_rate_anomaly.png
(roster-confounded), shape_maps.png (provisional), + earlier lfe_fronts_allyears.png (discovery
seeds, episode-bound per Merlin's 300-s window code-check), lfe_tidal_modulation.png (DEAD —
inherited PNSN phase), sse2012_fronts.png.

## FINAL VALIDATION ROUND (Merlin ruling + cor-discriminator, 2026-09-01 late)
- STEP-2 verdict was ambiguous (corr cdr-vs-cor maps +0.35); Merlin ruling: SOUTHERN
  "direct-dominated" feature RETRACTED (selection artifact + weakest resolution — phrase as "not
  robust to selection-immune measurement"); decay-map near-orthogonality carries NO physics content
  until split-half MAP reliability is measured; slope>1 in discriminators = expected Laplacian
  shrinkage (same physics as B4 transmission), corr 0.25-0.27 => map explains ~6% of contrast
  variance => claims limited to FIELD EXISTENCE + ONE REGIONAL SIGN, never cell structure.
- **COR-DISCRIMINATOR (the gate): PASS — slope +1.34, corr +0.27, permutation p = 0.0002** (bar
  p<0.01 positive). The "all-contaminated" failure configuration did NOT occur.
- STANDING RESULT (moderate confidence, Merlin-conditioned): **"A path-borne coda-enrichment field
  exists along the margin, elevated in the central section (44-48N)"** — with caveats: no cell-level
  structure claimable; no amplitudes (no lambda oracle; maps shrunk); south unresolved-and-retracted.
  Plus C3: co-located source-character variation ~0.26-0.33 Hz, 30-50x floor, robust to stack size
  (check i PASS); ordering replication 28/41 (check ii middling).
- REMAINING (Merlin d-list): split-half maps + attenuation-corrected inter-map corr (settles the
  three-orthogonal-maps question + whether decay/C2 has ANY map information); ringing variant
  ([3..5] or 4-8 Hz denominator); STAGE-3 source-FE solve + 2-D checker + lambda oracle (cor as
  primary input, cdr demoted to diagnostic) before anything beyond the narrowed claim.

## RESUME POINTER (session restart safe, 2026-09-02)
Nothing running. All products current on disk. Priority queue on resume:
1. Split-network coherence test on the 47.2N deep decline (cheapest decisive test).
2. End-to-end c-channel amplitude calibration (gates ALL slow-trend claims; Merlin B4 step-4 spec
   in notes/2026-08-22 design note — must run THROUGH the deconvolution pipeline).
3. Empirical per-claim nulls + split-half tau (supersede the placeholder 0.05% bar).
4. Fake-template duty calibration (gates rate/duty findings; 15.2B census made it non-optional).
5. Shape stage-3 (source-FE solve, cor as primary metric); movies + methods-artifact refresh.
User decisions open: mirror v2 (assumed dead, unconfirmed); June fault_4d*.npz legacy deletion.
Read first: this note + 2026-08-31 veto ledger + memory MEMORY.md.
