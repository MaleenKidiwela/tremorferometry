# LFE-feature exploration — DECISIONS LOG (rationale + rejected paths)

Companion to RESULTS.md (findings) and PLAN.md (goal/roadmap). This records *why* each choice was made
and what was tried and rejected, so the reasoning is reconstructable. Chronological within sections.

## A. Framing / strategy
1. **Pivoted from Lin time-coincidence → feature-space fingerprinting.** The B011-vs-Lin coincidence test
   was confounded (family coords are PNSN tremor-bin centroids not LFE relocations; lists noise-diluted;
   tremor = episode-level smear). Coincidence can't cleanly confirm "is this an LFE" → switched to
   classifying the *waveform features*.
2. **Lin (2023) as ground-truth positives; Lin-region stations PGC/B011/B926; B011 primary.** B011 is a
   PB borehole in the dense part of Lin's catalog, co-located with PGC, already trust-graded, and the SAME
   instrument family (EH? borehole) as the WA deployment targets → dissolves the VI→WA domain-shift worry.
3. **Deleted all stored waveforms except B927; re-downloaded PGC/B011/B926 in 3-component, Lin window
   (through 2017).** User decisions (AskUserQuestion): 3-comp for max optionality (polarization + future
   multi-station), Lin-window because that's where the labels are.
4. **Keep BOTH products:** the general continuous detector (R8) AND the tremor-window picker (R9) — per
   user ("general continuous is still good, don't scrap it"). They answer different questions.

## B. Windows & features
5. **40 s windows [t−10, t+30], event/S at +10 s.** "Long enough" to hold pre-event noise (for SNR),
   emergent onset, S, and coda; the +10 s offset anchors all classes consistently.
6. **23 hand-crafted features** (spectral shape + envelope/onset + 3-comp polarization), RandomForest.
   Chosen over deep features first because interpretable, physical, cheap, and the project already proved
   spectral shape separates classes (centroid/kurtosis AUC 0.92). The autoencoder (R6) is the learned check.
7. **Dropped Lin's depth filter; curate with N≥6, residual<1.2 (≤35 km).** Lin depth is bimodal-unreliable
   (quantiles −60/−59/−30/−3/0) and over-cut positives 11.6k→1.9k. Do NOT trust Lin depth.
8. **SNR/amplitude control** (re-score dropping all amplitude features) baked into every analysis — to prove
   discrimination is waveform *character*, not loudness. (It held: AUC unchanged, snr d≈0.)
9. **Negative classes built in stages:** RAND (quiet noise) → +EQ (ANSS earthquakes) → +BLAST (ANSS
   explosions) → +TNOISE (random in-tremor times not near an LFE = the real in-tremor confuser). Driven by
   the goal sharpening to "reject EQ/blast/cultural inside tremor windows."
10. **Event anchoring:** EQ/BLAST at estimated S = OT + hypocentral-dist/3.6 km·s⁻¹ (distances vary, so per-
    event); LIN at OT+11 s; TNOISE/candidates at the detection time itself.

## C. Skeptical corrections (artifacts caught — see RESULTS R2, R7, R8, R10)
11. **Family-stack verdict, alignment confound:** family stacks are matched-filter-aligned (sharp), my LIN
    stacks are OT-anchored with S-time jitter (smeared) → unfair on onset/envelope features → re-scored on
    spectral+polarization (alignment-robust) only.
12. **Classifier extrapolation:** the robust score *fell* because family stacks lie *beyond* the LIN cloud on
    the low-freq side → trusted RAW spectral features over the classifier P (the spectra showed families are
    *more* LFE-like). Lesson: never trust a classifier P at the edge of its training distribution.
13. **Degraded positive reference:** LIN *stacks* are HF-inflated by anchor jitter → poor positive; precise
    alignment is needed for stack comparisons.
14. **Per-sample picker (R7):** widened LFE label σ to match true arrival uncertainty + early stopping;
    **precise-label refinement REJECTED** — a global LFE template aligns at only CC 0.22 (LFEs patch-specific
    + emergent), retrain unchanged → concluded the limit is structural (single-station + emergent), not labels.
15. **Continuous scanner (R8):** first read "P(LFE)=0.94 everywhere = broken" was WRONG; the quiet-day 0.05
    proved it discriminates (ETS days really are LFE-rich). It's an activity detector, not an event timer.
16. **Discovery filter placement (R10):** per-candidate filtering fragments families (kept only 19% of LFE
    families) → switched to FAMILY-level filtering (cluster first, score family by mean member P(LFE)).

## D. Discovery integration (R10)
17. **Scoped to 2010–2013** (Lin-dense, multi-ETS, tractable: 415k candidates scorable in ~3 min).
18. **Faithful discovery params:** fs=100 (borehole), cc≥0.8, complete-linkage, ≥3 members, **min-years=3**,
    max-bin 2000 — identical to the production `discover_gpu` so the comparison is real.
19. **Ground-truth family classification** by member coincidence (±6 s) with Lin/EQ/blast; explicitly noted
    Lin-confirmation UNDERESTIMATES purity (UNK families include real LFEs Lin doesn't catalog at B011's edge).
20. **Baseline = scored candidates (valid features), Filtered = P(LFE)≥thr subset** — controlled comparison
    (same pool, only the filter differs).

## E. Tooling / environment (ephemeral pod)
21. Reinstalled obspy, pyarrow, torch (CPU), umap-learn, hdbscan, cupy-cuda12x (`--user`, persists on NFS);
    cupy needs `CUDA_PATH=/opt/conda/targets/x86_64-linux`.
22. **Perf:** picker dataset build was bandpass+resampling whole days (~55 min projected) → slice-then-filter
    short windows (~75 s). ~100× speedup.
23. **pandas traps fixed:** `to_datetime` returns datetime64[us] here → force ns before int64 (else 1000×
    time mis-scale); `.to_numpy()`/de-arrow before sklearn (pyarrow-backed cols break indexing).
24. **bracket-trap (memory):** `pgrep -f build_picker_dataset` matched my own shell command → false "still
    alive"; the Python proc was already dead.

## F. Open decisions / not yet done (for the next session)
- Test the tremor-picker filter at an **EQ/blast-contaminated station** (B018 anthropogenic / B935 natural-
  seismicity) where its 0.97 EQ/blast rejection actually bites — B011 has almost no EQ/blast families so the
  win is muted there.
- **AE-latent concat** to attack the hard in-tremor-noise residual (0.81) — needs an AE *retrain* (only
  latents were saved, not weights) on spectrograms incl. TNOISE/BLAST.
- **Multi-station** per-sample picker (the structural lever for R7) — project Phase C.
- Replicate at B926; per-family-aligned family verdict (R2 finish).
