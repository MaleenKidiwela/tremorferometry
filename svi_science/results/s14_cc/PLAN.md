# S14 — does the stretching CC show changes that dv/v does not? (rules fixed first)

Written Wed 2026-10-14, 1:59–2:20 pm PT, before any CC-vs-time number of this study was computed.
Question (user, 2026-10-14): "sometimes CC shows things that dv/v doesn't show — have you looked into it?" Clarified: the CC
is **the CC of the LFE dv/v stretching measurement itself** (the colour of the points on our dv/v plots). So the primary
observable is that existing CC; nothing new is defined for the headline. Wording: "scrambled control" = labels or dates
shuffled within a family-station series (keeping the season); "c2 control" = the detector's noise windows (station data
shifted to another time of the same day) through the identical path. Times PT. dv/v = −ε, in %.

## 0. Status
- 2026-10-14 ~3:05 pm PT — pilot, full run (SVI + WA-NW) and post hoc done; `RESULTS.md`. Amendment 1 in NOTES.

## 1. What CC already exists (inventory, 1:40–1:58 pm PT)
| level | file | content | reference |
|---|---|---|---|
| episode, SVI, 76 fam × 7 stns (PGC, B009, B010, B011, B926, VGZ, SNB), 8,946 cells | `s11_depth/tables/episode_raw.parquet` | ε and CC of the joint 3C S-anchored stretch at 2–5 s (headline), 3–6, 2–3.5, 3.5–6 s; Z-only and H-only 2–5 s; 2-parameter (shift ±50 ms × stretch) CC; c2-control CC per window; n, ΣP, n_days, **pre_rms_ratio** (pre-P RMS / 2–5 s RMS, Z, of the episode stack) | leave-episode-out within response epoch |
| episode, WA-NW, 78 fam × 75 stns, 18,428 cells | `s11_depth/wa_nw/tables/episode_raw.parquet` | same columns | same |
| episode, SVI, Z only (S10 multi) | `s10_dvv/multi/tables/episode_stacks.parquet` | cc24, cc35 (Z, 2–4 / 3–5 s) + c2 | leave-episode-out |
| day, SVI 64 fam / WA 12 fam (round-7 headline series) | `s11_depth/round7/{svi,wa_nw}/tables/r7_front_days_{near,far}.parquet` | day ε and CC (3C 2–5 s), weight ws, data + c2 rows; front times `r7_front_times.csv` | series grand stack |
| day, SVI fams 10/17/87, 20 stns (S10 pilot) | `s10_dvv/pilot/tables/daily.parquet` | cc24/cc35 (Z) + c2 | epoch grand stack |
| day, fam 10 at PGC (the user's figure series) | `s10_dvv/pilot/tables/fam10_PGC_3c25_daily.csv`, `fam10_PGC_episode_dvv.csv` | raw and SVD 3C 2–5 s / Z 2–4 s ε and CC | epoch grand stack |
| day stacks (to recompute anything) | `s11_depth/tables/store/<NET.STA>.npz`, `wa_nw/tables/store/`, `round3/*/tables/store_ccw/` | P-weighted 3C day stacks St, weights Ws, grand stack G1, tS, c2 day stacks | — |
Not existing: CC of the direct-S window; CC between successive episodes; per-day pre-arrival noise. These are computed from the
stores (§4) as secondary columns.

## 2. Expected size (rough, stated before measuring)
- **RF anchor** (rf_gosselin §1–11): PGC ΔCC_SV −0.12 (15 d post vs pre), about half from a pre-ETS high (+0.04) and half a
  post-ETS low (−0.06), lowest 6–15 d after the GPS midpoint, back by ~24–27 d. That is a single converted phase that
  crosses the LVL once.
- **LFE coda sees the LVL only partly.** From the S11 kernels (ell 40 km, `kernel_mc.npz`), coda paths in 2–5 s past S spend
  2.8–3.1 s inside the 10-km source sphere against 2.4–2.7 s for the direct S: ~0.1–0.4 s extra, i.e. 1–3 % of the lapse
  time, of which a 3.5-km layer below the LFEs holds roughly a quarter to a half.
- (i) **Velocity change in the LVL only** (Gosselin dVs +0.05–0.1 km/s ≈ +2–4 %): mean apparent stretch ≈ 0.03 × 0.1 s /
  3.5 s ≈ 0.05–0.1 % (already near our dv/v bounds 0.02–0.05 %, so mostly excluded as a stretch); the decorrelating part is
  the path-to-path spread of LVL residence (σ ≈ 0.1–0.2 s → σ_δt ≈ 3–6 ms) → 1−CC ≈ 2π² f² σ_δt² ≈ **0.003–0.01** at 4 Hz.
- (ii) **Scatterer / impedance change in the LVL** (fluid pressure): the S reflection off the LVL bottom (3.5 km below the
  LFE, two-way ≈ 2.5–3 s past S) lands in the 2–3.5 s window, amplitude ≈ 0.1–0.2 of direct S, i.e. a 10–30 % share of
  that window's energy. If that phase changes like the RF (ΔCC −0.12), the window CC changes by ≈ share × 0.12 ≈
  **0.01–0.04**, largest at 2–3.5 s, smaller at 3.5–6 s, ~0 in the direct-S window.
- Stated expectation: **|ΔCC| ≈ 0.005–0.04 (episode or post-ETS time scale), largest in 2–3.5 s, sign negative after
  the front / after ETS.** In Fisher z at CC 0.5–0.8 this is Δz ≈ 0.007–0.1. A timing glitch of δt gives 1−CC ≈ 2π²f²δt²
  = 0.03 at 10 ms, 4 Hz — the same size; hence the clock checks below.

## 3. Confounds and how each is handled (fixed)
- **Stack size and noise.** CC rises with n and falls with noise. Value used: Fisher z = atanh(clip(CC, −0.99, 0.99)).
  Covariates: log n (episode; day level: log ws), log pre_rms_ratio (episode; day level: recomputed per day stack), log n_days.
  Three forms of every comparison:
  - **REG (headline):** ANCOVA, z = α_series + β·covariates + γ·class (+ season terms where stated), OLS, series fixed
    effects; γ is the class contrast adjusted for n and noise. β common per region.
  - **RAW:** the same without covariates (shown, not used for the verdict).
  - **MATCHED:** within series, each cell of one class paired to the nearest cell of the other class with |Δlog n| ≤ 0.3
    and |Δlog r| ≤ 0.3 (no reuse); mean paired Δz.
- **Season.** Every scrambled control shuffles labels only among cells within ±45 days of year in the same series.
- **Instrument / clock.** Flags: PGC 2019-11, 2021-01 episodes; B010 2009-03, 2022-03, 2024-08 (ballistic RESULTS); B013
  after 2015-01-01 (WA). Diagnostic column Δcc_shift = cc2p_25 − cc3c_25 (CC recovered by a time shift); cells with
  |shift2p_25| ≥ 20 ms flagged. Every verdict is re-run with flagged cells dropped.
- **Activity noise (tremor in the window).** The c2-control CC column goes through the same comparison; a contrast that the
  c2 column also shows at > half the size is called activity noise, not a medium change.
- **dv/v leaking into CC.** Stretch is fitted, so a pure stretch does not lower CC; but railed or mis-fitted ε could.
  Column: |ε anomaly| as an extra covariate (REG+ε); verdict must hold in it.

## 4. Steps
**Step A (first, raw):** side-by-side time series of the existing dv/v and its CC (same series, raw values, points coloured
by CC, no error bars, ETS majors shaded): (1) family 10 at PGC daily (the user's plot series, raw and SVD); (2) families 10, 87,
17 at PGC, B009, B010, B011, episode level 2–5 s; (3) pooled SVI and WA-NW (median over series of series-demeaned values per
episode). Then the list "CC changes while dv/v stays flat" (§5 C7). This figure goes at the top of RESULTS.md.
**Step B (store columns, secondary):** from the stores, rebuild each episode stack (Σ Ws·St over the episode's days) and its
leave-episode-out reference; check it reproduces cc3c_25 (rule: median |Δcc| ≤ 0.005, else stop and fix); compute CC of the
direct-S window (−0.3…+1.5 s past S, same joint 3C stretch), CC between successive episode stacks (2–5 s, after stretch),
and per-day pre_rms_ratio for the day-level tests.
**Step C:** the comparisons C1–C7 on the pilot; **Step D:** if the pilot is sane (§6), all SVI and WA-NW series.
Script names `svi_science/scripts/s14_*.py`; ≤ 4 workers, OMP_NUM_THREADS=1, RAM < 60 GB; nohup + DONE markers in this folder.

## 5. Comparisons (fixed; windows 2–5 s headline, also 2–3.5, 3.5–6, 3–6, Z 2–5, H 2–5, direct S)
- **C1 ETS vs between:** major-episode cells vs burst cells (REG γ = major − burst).
- **C2 RF-like pre vs post ETS:** burst cells whose mid date is 1–45 d after a major's end ("post") vs 1–60 d before a
  major's start ("pre"); γ = post − pre. RF predicts negative. Also post-only vs all other bursts.
- **C3 within episode (day level, front-linked):** r7 day rows of majors, lag = day − front date at that family
  (`r7_front_times.csv`, headline set), bins ≤ −1, 0, 1, 2, 3, 4–6, ≥ 7; REG with series×major fixed effects,
  covariates log ws and log day-noise; contrasts **late (≥ 3) − early (0–1)** (RF-like: lowest late) and **pre-front
  (≤ −1) − early**. Scrambled control: each family's relative front offsets deranged among its majors, 500 draws.
- **C4 episode to episode:** (a) cross-family coherence of REG residuals (no class terms): M = median over families of
  the correlation of a family's residual series with the mean of the other families' (same station, common episodes ≥ 8);
  (b) co-located: correlation of PGC vs B009/B010/B011 residuals for the same family-episode, and of the site means
  across episodes; (c) successive-episode CC (Step B) vs gap and vs class of the later episode. Scrambled control: episode
  labels shuffled within series within ±45 d of year, 1000 draws.
- **C5 seasonal:** annual harmonic of REG residuals pooled per region; scrambled = per-series random circular day-of-year
  shift, 1000 draws. Compared with the dv/v season (0.07–0.08 %, peak ~21 Aug–3 Sep).
- **C6 long-term:** Theil–Sen slope of residual vs year within each epoch (column; reference drift caveat).
- **C7 CC-only list:** a cell is "CC low, dv/v flat" if its REG residual (no class terms) < −3 robust σ of the series
  (MAD × 1.4826) and |ε − series median| < sig3c_25. An episode is a **CC-only episode** at a station if ≥ 3 families flag
  it there; it is **site-common** if the co-located instruments also average < −1 robust σ. Listed with n, noise, Δcc_shift,
  glitch flag. Expected number from the scrambled control (labels shuffled; same rule).

## 6. Pass rules (fixed)
A comparison (C1, C2, C3, C5) is **"CC change detected (provisional)"** only if all hold:
1. scrambled-control p < 0.01 two-sided (pilot) / < 0.003 (full run) for REG γ;
2. MATCHED γ has the same sign and differs from REG by < 2σ;
3. co-located: same sign at ≥ 3 of the 4 PGC-site instruments (SVI); WA: same sign at both multi-instrument sites
   (B005/B006/B007/LRIV and B013/B943) where measurable;
4. odd- vs even-numbered families: same sign;
5. holds (rule 1) with glitch-flagged cells and |shift2p| ≥ 20 ms cells dropped, and in REG+ε;
6. c2-control column does not show the same-sign contrast at > half the size.
Otherwise: no detected change; report the 2σ bound on ΔCC (and Δz) with its scope — the window (lapse → the volume the
kernel samples), the time scale (episode / days since front / weeks after ETS), and the comparison.
C4 coherence: M above the 99th percentile of the scrambled control = "episode CC anomalies are shared"; then the co-located
correlation says whether it is site-common (instrument or very shallow) or path-specific.

**Pilot sanity (go to full run):** Step B reproduction rule passes; β(log n) > 0 and β(log noise) < 0 in REG; scrambled γ
centred on 0 (|median| < 0.2 σ). Pilot = families 10, 87, 17 at PGC, B009, B010, B011 (co-located: an independent-instrument
check). Full run = all SVI series (7 stations) and all WA-NW series.

## 7. Plots
Raw series first. CC and dv/v series: points coloured by CC, no error bars, line graphs allowed. σ in tables/text only.
