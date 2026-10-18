# RESULTS — rf_gosselin (core stages; written 2026-10-13 ~10:40 am PT)

Stages done: RF_CORE, REPRO (1993–2019, 0.05–0.5 Hz), REPRO07 (0.05–0.7 Hz), OOS (2019-06..2026), FULL (1993–2026).
Pending in the queue: recovery curves (REC_REPRO, REC_FULL), the extension stations (SVI NLLB/YOUB/MGB/GOBB; WA GNW/SP2/
LON/OFR/WISH/NLWA), and slip-front timing (SLIP_SVI, SLIP_WA). Each section is added here when its DONE marker appears.

Main figure: `figures/repro/fig2_CN.PGC.png` (the paper's Fig. 2 layout). Tables: `tables/<tag>/summary.csv`,
`tables/<tag>/<STA>_result.json`. Δ = mean(post 15 d) − mean(pre 15 d) around the ALBH GPS midpoints. The scrambled
control shuffles RF dates among the events 20,000 times. "pct" = fraction of scrambled draws at least as extreme in
the expected direction (one-sided).

## 1. Reproduction, 1993-07 … 2019-06, 21 ETS (the paper's count), 0.05–0.5 Hz

| station | RFs | ΔCC_SV 15 d | Welch p | pct | ΔCC_SV 150 d | ΔCC_SH 15 d | ΔdVs top 15 d (km/s) | pct | ΔdVs bottom (km/s) | pct |
|---|---|---|---|---|---|---|---|---|---|---|
| **PGC** | 1224 | **−0.122** (−15 % of median 0.79) | **0.002** | **0.002** | −0.036 | −0.075 (p 0.25) | +0.52 | 0.18 | +0.18 | 0.66 |
| SNB | 797 | −0.059 | 0.52 | 0.27 | +0.009 | −0.258 (p 0.03) | −0.25 | 0.75 | −0.34 | 0.18 |
| VGZ | 836 | −0.146 (−30 % of 0.49) | 0.13 | 0.09 | −0.012 | −0.050 | +0.03 | 0.47 | −0.43 | 0.04 |
| LZB | 510 | +0.025 | 0.75 | — | −0.077 | −0.319 (p 0.01) | −0.20 | — | −0.54 | — |
| PFB | 561 | −0.094 | 0.39 | — | −0.040 | −0.108 | +0.25 | — | −0.16 | — |

What this shows:
- **PGC reproduces the paper's main result.** After ETS, SV receiver functions match the long-term harmonic model worse.
  The drop is −0.12 in the first 15 days, beyond the scrambled control at the 99.8 % level. It stays outside the 95 % band
  for 15–60-day bins and shrinks with bin width (−0.122, −0.070, −0.060, −0.056, −0.031, −0.029, −0.036 for 15…150 d),
  which is the same decay shape as the paper's Fig. 2E.
- **VGZ** has the same sign and a larger drop (−0.15), but it is not significant (9th percentile, p = 0.13; only 52 RFs in the bins).
  Its LVL-bottom contrast falls (−0.43 km/s, 4th percentile), the "bottom reduction" the paper reports.
- **SNB**: same sign, small (−0.06), not significant. Its dVs-top has the wrong sign.
- LZB and PFB show no significant SV change, as in the paper.
- **dVs size is not reproduced.** PGC's LVL-top dVs rises (+0.52 km/s, right sign) but sits inside the scrambled band (82nd
  percentile). Our single-RF dVs scatter is about ±1.5 km/s, so the paper's +0.05–0.1 km/s cannot be resolved with our
  amplitude conversion. Our conversion probably differs from theirs (slab-normal correction, depth picks).
- Their SH was flat. Ours is also not significant at PGC (−0.075, p = 0.25). It is significant at SNB and LZB, which the
  paper did not see.

**Score (PLAN §3, fixed before results): partly reproduced.** ΔCC_SV has the right sign and size at all three of PGC,
SNB and VGZ, but is significant only at PGC. The dVs-top size test fails everywhere (PGC is 5× too large and not significant).

## 2. 0.7 Hz variant (1993–2019)
- PGC: ΔCC_SV −0.125 (p 0.002, pct 0.002). It decays the same way, so the higher band confirms PGC.
- VGZ: −0.140 (pct 0.08). SNB: −0.061 (pct 0.22).
- dVs top: PGC +0.45 (pct 0.13), VGZ +0.15. VGZ bottom −0.57 (pct 0.007). PFB top +0.87 (pct 0.01; unexpected, single station).

## 3. Out of sample, 2019-06 … 2026 (9 ETS; the harmonic model is built from this period alone)
- **PGC: ΔCC_SV −0.094 (pct 0.04, p 0.10), 363 RFs.** The 30–60-day bins are outside the 95 % band (−0.100, −0.075, −0.067).
  **The signal replicates independently**, slightly smaller.
- VGZ: −0.31 (pct 0.03), but only 12/5 RFs in the bins.
- SNB: closed 2021-03 (36 RFs), no test.

## 4. Full record, 1993–2026 (30 ETS)
- **PGC: ΔCC_SV −0.115, p 0.001; no scrambled draw out of 20,000 is as low.** It decays to −0.030 at 150 d.
  dVs top +0.56 (pct 0.035, now one-sided significant).
- VGZ: −0.146 (pct 0.056); bottom −0.46 (pct 0.014).
- SNB: −0.03 (no change).

## 5. Timing (provisional; from the test run of the slip-front code, SVI 2010–2025)
Tremor within 30 km of PGC peaks at a median of +2.5 days after the GPS midpoint (VGZ +0.1 d). The first tremor arrives
about 10–15 days before. So the post-ETS RF change starts as the local slip front passes the station.

## Pending (queue, see logs/queue.log)
REC_REPRO, REC_FULL (recovery curves, 3-day steps); RF_EXT, EXT_SVI, EXT_WA; SLIP_SVI, SLIP_WA (these rewrite §5).

---
# Added 2026-10-14 (PT): recovery, extension and slip-front stages

The queue finished 2026-10-14 1:13 pm PT (QUEUE_DONE); every stage in the "Pending" list above is done. §5 above is
superseded by §9. Recovery curves (`scripts/rf05_recovery.py`): 15-day (and 7-day) sliding windows, centres every 3 days
(PLAN E4 said 7-day steps) from −60 to +120 d; anomaly = window mean − mean of all RFs of the run; grey band = 95 % of a
scrambled control (RF values shuffled among event dates, 5,000 draws). With 61 centres, about 3 centres per curve fall
outside a 95 % band by chance, so single isolated excursions mean little; runs of adjacent centres matter.

## 6. Recovery curves, core stations (REC_REPRO 1993–2019, REC_FULL 1993–2026)
Tables `tables/{repro,full}/recovery_<STA>.csv`, figures `figures/{repro,full}/recovery_<STA>.png`. ΔCC_SV, 15-day window.

| station, run | 15 d before mid (anomaly) | lowest after mid | at centre | outside 95 % band (centres, d) | back to 0 at | 7-day window lowest |
|---|---|---|---|---|---|---|
| **PGC full** | **+0.044** | **−0.061** | **+9 d** | **below: +9, +12, +15**; above: −18, −12, −9 | **+27 d** | −0.080 at +6 d |
| PGC repro | +0.061 | −0.060 | +9 d | below: +9; above: −15, −12, −9 | +21 d | −0.121 at +6 d |
| VGZ full | +0.060 | −0.128 | +9 d | none (band −0.148) | +33 d | −0.174 at +6 d |
| SNB full | +0.086 | −0.089 | +24 d | none | +33 d | — |
| LZB (= repro) | +0.141 | −0.007 | +48 d | curve sits above the band from −60 to +18 d | — | — |
| PFB (= repro) | +0.110 | −0.098 | +24 d | none | +42 d | — |

(LZB and PFB stop in 2014 and 2018, so their repro and full curves are the same RFs.)

Other quantities at PGC (full): CC_SH dips briefly at +3 d (−0.059, inside the band). dVs top: nothing outside the band near
the ETS (+0.24 km/s at +9 d, band ±0.43). **dVs bottom drops at +21 to +30 d** (−0.73 km/s at +27 d, below the band at 4 adjacent
centres; repro only at +33 d). VGZ dVs bottom is below the band at +6 and +18 d.

What this shows:
- **The PGC change is short.** CC_SV starts falling at the GPS midpoint, is lowest 6–15 days after it (−0.06 against the
  long-term mean; −0.08 with 7-day windows), and is back at the long-term level by window centre +24 to +27 days. So the
  change lasts about **3 weeks** and has recovered within about **a month**. This matches the paper's "about two weeks"
  and the bin-width decay in §1.
- **About half of the 15-day Δ is a high before the ETS.** At PGC, CC_SV in the 2–3 weeks before the midpoint is +0.044
  above the long-term mean (outside the band at −18, −12, −9 d) and −0.05 to −0.06 below it in the 2 weeks after. The pre-vs-post
  Δ of −0.12 (§1) is the sum of the two. VGZ has the same shape (+0.06 before, −0.10 after), not significant.
- The LVL-bottom contrast at PGC falls 3–4 weeks after the midpoint, later than the CC drop. Its size (−0.7 km/s) has the same
  scale problem as dVs top (§1), so read it as a timing, not a size.
- SNB, LZB, PFB: no post-ETS drop outside the band. LZB's curve is offset upward over −60…+18 d (not an ETS drop; unexplained).
- Recovery time scale for the LFE search: **onset at the midpoint (within ~3 days), strongest days +3 to +15, gone by ~+25 d.**

## 7. Extension, southern Vancouver Island (EXT_SVI + REC_EXT_SVI; 1993–2026, ALBH GPS midpoints, 30 ETS)
Tables `tables/ext_svi/`, figures `figures/ext_svi/`. Δ = mean(post 15 d) − mean(pre 15 d); pct = one-sided scrambled-control
percentile in the expected direction (CC down, dVs top up); 20,000 draws.

| station (years) | RFs | n pre/post | ΔCC_SV 15 d | Welch p | pct | ΔCC_SV 30 d (pct) | ΔCC_SV 150 d | ΔCC_SH 15 d (p) | ΔdVs top (pct) | tremor ≤30 km / episode |
|---|---|---|---|---|---|---|---|---|---|---|
| **NLLB** (2004–26) | 1128 | 45/37 | **−0.089** | 0.046 | **0.040** | −0.075 (0.021) | −0.009 | −0.094 (0.24) | +0.10 (0.72) | 84 |
| **YOUB** (2003–14) | 436 | 15/13 | **−0.205** | 0.088 | **0.026** | −0.073 (0.16) | −0.010 | −0.143 (0.32) | +0.46 (0.88) | 1767 |
| GOBB (2017–26) | 303 | 15/7 | −0.109 | 0.41 | 0.15 | −0.050 (0.29) | +0.006 | −0.359 (0.12) | −0.11 (0.39) | 323 |
| MGB (2005–17) | 600 | 29/23 | +0.013 | 0.86 | 0.59 | +0.005 (0.55) | −0.005 | +0.009 (0.92) | −0.09 (0.39) | 974 |

Recovery curves (15-day window): NLLB's curve is **above** the band before the ETS (−21 to −9 d) and only −0.03 at its lowest
(+30 d, inside the band); YOUB is lowest at +6 d (−0.14, inside its wide band, 12 RFs), back to 0 by +21 d; GOBB and MGB show
nothing outside the band near the ETS (MGB one centre at −3 d).

What this shows:
- **Two of four Vancouver Island stations show the same sign as PGC, beyond the 95 % scrambled control: NLLB −0.089 (pct 0.040)
  and YOUB −0.205 (pct 0.026).** GOBB has the same sign but is not significant (only 7 RFs in the post bin). MGB shows no change.
- At NLLB the change holds in 30-day bins (pct 0.021) and fades by 120–150 d, the same decay shape as PGC. But its recovery
  curve says the Δ comes from **higher CC before the ETS** (+0.07), not a post-ETS drop (+0.01). At YOUB it is mostly a
  post-ETS drop (+0.06 before, −0.14 after) but rests on 15 + 13 RFs.
- Scored with the PLAN §3 station rule as a yardstick (not pre-set for extension stations): NLLB passes R1–R5 for CC and R1–R2
  for dVs top (+0.10 km/s, not significant). YOUB passes R1 and R3 but its size (−0.205) is just over the R2 upper limit (0.20)
  and its dVs top (+0.46) is outside R2.
- Scope of the no-change: MGB (west coast, 2005–2017, 23 post-ETS RFs) shows no 15-day pre/post CC_SV change around ALBH midpoints.

## 8. Extension, Washington (EXT_WA + REC_EXT_WA; 2007–2026, SEAT GPS midpoints, 27 ETS)
Tables `tables/ext_wa/`, figures `figures/ext_wa/`. Same columns as §7.

| station | RFs | n pre/post | ΔCC_SV 15 d | Welch p | pct (CC down) | ΔCC_SV 30 / 150 d | ΔCC_SH 15 d (p) | ΔdVs top (pct) | tremor ≤30 km / episode |
|---|---|---|---|---|---|---|---|---|---|
| **LON** (Mt Rainier) | 1236 | 62/66 | **−0.118** | **0.003** | **0.003** | −0.093 / −0.034 (pct 0.016) | −0.054 (0.08) | pick failed | 0 |
| GNW (Kitsap) | 1351 | 72/78 | **+0.075** (opposite sign) | 0.09 | 0.97 (0.034 for CC up) | +0.015 / −0.024 | −0.008 (0.88) | −0.42 (opposite sign, 0.02 for down) | 905 |
| SP2 (Seattle) | 1131 | 60/56 | +0.015 | 0.68 | 0.65 | −0.003 / −0.022 | −0.057 (0.33) | pick failed | 42 |
| OFR (Olympic coast, 2007–13) | 316 | 19/20 | −0.016 | 0.82 | 0.42 | −0.036 / −0.004 | −0.006 (0.94) | pick failed | 0 |
| WISH (Grays Harbor) | 1140 | 60/68 | +0.001 | 0.96 | 0.52 | −0.006 / +0.012 | −0.093 (0.28) | pick failed | 1 |
| NLWA (SW Olympics) | 1320 | 66/77 | −0.010 | 0.79 | 0.39 | −0.012 / +0.005 | −0.139 (0.055) | pick failed | 1 |

Recovery curves: LON is above the band before the ETS (−15 to −6 d, +0.07) and only −0.03 after (inside the band). GNW is below
the band at +24/+27 d and +54 to +66 d (−0.088 at +54 d), a late dip, not the PGC shape. SP2, OFR, WISH, NLWA: nothing outside
the band near the ETS.

What this shows:
- **One Washington station, LON, shows the PGC sign at high significance: ΔCC_SV −0.118, pct 0.003, p 0.003**, and it decays with
  bin width. But (a) LON is ~50 km east of the tremor band (no PNSN tremor within 30 km in any of 24 episodes), its slab is at
  ~72 km, below the 0–60 km depth range we use; (b) its Δ is mostly a pre-ETS high (+0.07 before vs −0.03 after); (c) unlike
  PGC it stays outside the band even with 150-day bins (−0.034, pct 0.016), which looks like a slow or seasonal pattern rather
  than a 2–3-week change. WA midpoints cluster in May and August (13 of 27) and LON CC is +0.04 in July and −0.03 in Sep–Oct, so
  a season effect is possible; the scrambled control does not keep the season. Not tested.
- **The one WA station that sits on the tremor, GNW, shows the opposite sign:** CC_SV rises (+0.075, 3.4 % tail for a rise) and its
  LVL-top dVs falls (−0.42 km/s, pct 0.02), at 72/78 RFs in the bins.
- Scope of the no-change: SP2, OFR, WISH and NLWA show no 15-day pre/post CC_SV change around the SEAT-GPS midpoints, 2007–2026.
  OFR, WISH and NLWA have almost no tremor within 30 km (≤ 1 event per episode), so for them this is not a test of the LVL
  above tremor; SP2 is at the edge (42 tremor events per episode).
- Across the 10 extension stations, 3 pass the one-sided 5 % level in the PGC direction (NLLB, YOUB, LON; ~0.5 expected by
  chance, binomial P(≥3) ≈ 1 % if independent) and 1 passes it in the opposite direction (GNW).

## 9. Timing: GPS midpoint vs local slip front (SLIP_SVI, SLIP_WA; supersedes §5)
`tables/slipfront_SVI_ALBH.csv`, `tables/slipfront_WA_SEAT.csv`: per station and episode (2009-06 onward, PNSN tremor within
30 km of the station; PLAN §5 said 50 km), times relative to the GPS midpoint of the 2nd percentile ("first"), median and 98th
percentile ("last") of local tremor. Medians over episodes:

| station | episodes with tremor | first | median | last | spread of median (IQR) | tremor / episode |
|---|---|---|---|---|---|---|
| PGC | 17/17 | −14.5 d | **+2.5 d** | +6.1 d | 3.6 d | 911 |
| VGZ | 17/17 | −9.7 | +0.1 | +6.5 | 4.2 | 1410 |
| SNB | 16/17 | −14.1 | +1.3 | +5.9 | 8.4 | 111 |
| LZB | 17/17 | −6.5 | +3.6 | +9.9 | 4.9 | 2029 |
| PFB | 15/17 | −3.1 | +3.2 | +11.3 | 13.0 | 981 |
| NLLB | 16/17 | −16.7 | +1.3 | +9.6 | 12.3 | 84 |
| YOUB | 16/17 | −6.5 | +2.2 | +10.2 | 16.0 | 1767 |
| MGB | 16/17 | −8.3 | +0.5 | +5.2 | 20.2 | 974 |
| GOBB | 17/17 | −17.3 | −1.8 | +6.6 | 14.1 | 323 |
| GNW (WA) | 24/24 | −15.5 | +0.9 | +9.9 | — | 905 |
| SP2 (WA) | 23/24 | −21.6 | +0.2 | +22.2 | — | 42 |
| LON, OFR, WISH, NLWA | 0–1 /24 | no local tremor | | | | 0–1 |

Alignment check (quick re-binning of the existing per-RF tables, `scripts/rf06_align_check.py` → `tables/align_check_15d.csv`;
same RFs, only the episodes with a local tremor time, 2010–2025; scrambled control 4,000 draws). ΔCC_SV, 15-day bins:

| station | episodes | aligned on GPS midpoint | on median local tremor | on first local tremor |
|---|---|---|---|---|
| **PGC** | 17 | **−0.133 (pct 0.0003)** | **−0.128 (pct 0.001)** | −0.013 (pct 0.38) |
| VGZ | 17 | −0.088 (0.21) | +0.069 (0.71) | +0.095 (0.82) |
| NLLB | 16 | −0.100 (0.041) | +0.022 (0.63) | +0.098 (0.95) |
| YOUB | 16 | −0.085 (0.29) | +0.188 (0.90) | −0.054 (0.35) |
| GOBB | 17 | −0.109 (0.15) | +0.095 (0.80) | +0.015 (0.55) |
| MGB | 16 | +0.143 (0.97) | −0.028 (0.38) | +0.046 (0.73) |
| SNB | 16 | +0.063 (0.68) | −0.086 (0.26) | +0.098 (0.77) |
| GNW | 24 | +0.041 (0.84) | −0.015 (0.36) | +0.024 (0.71) |
| SP2 | 23 | +0.024 (0.71) | −0.031 (0.23) | +0.040 (0.86) |

What this shows:
- **At PGC the change lines up with the main tremor passage under the station, not with the first tremor nearby.** Aligned on
  the midpoint or on the median local tremor time the Δ is the same (−0.133 vs −0.128, both beyond the scrambled control);
  aligned on the first local tremor (about 2 weeks earlier) it disappears (−0.013).
- **Midpoint vs local front cannot be told apart at PGC**: local tremor peaks only +2.5 days after the ALBH midpoint (IQR 3.6 d),
  less than the 15-day bins and the RF sampling (~1 RF per week). Neither alignment is better.
- At the other stations the test has little power: few RFs per episode, and per-episode tremor times scatter by 8–20 days
  (IQR), so re-aligning moves RFs across the bin edge. NLLB's midpoint-aligned Δ does not survive local-tremor alignment, but
  NLLB has sparse tremor (84 per episode, IQR 12 d), so its local time is poorly defined.

## 10. Checks and things that look broken
- **US.NLWA (re-downloaded 10-14 after the '00' location-code fix) looks sensible**: 2,312 of 2,527 events with data, 1,343
  RFs above 5 dB (58 %, same pass rate as GNW and SP2), 1,320 in 2007–2026 (~66 per year; 2022 only 22), median CC_SV 0.77,
  coherent harmonics. Its CC result (no change) is usable.
- **dVs is not usable at SP2, LON, NLWA, OFR, WISH**: the LVL picks failed. SP2 top = bottom = 58 km (window edge); LON 57/59 km
  (slab ~72 km, outside the 0–60 km range); NLWA top at 15 km (window edge, slab ~28 km); OFR top at 39.25 km (window edge);
  WISH top has a positive amplitude (no negative pulse found). Their dVs columns (e.g. NLWA dVs top +0.19 km/s at every bin
  width, pct 0.999 at 150 d) are not LVL measurements. CC_SV/CC_SH are unaffected (0–60 km).
- Thin bins: GOBB 15/7, YOUB 15/13, OFR 19/20 RFs; their recovery curves have empty 7-day windows (source of the
  "All-NaN slice" warnings in the REC_EXT logs; harmless).
- The scrambled control shuffles RF dates without keeping the calendar month. ETS midpoints cluster by season (SVI Aug–Sep,
  WA May and Aug) and CC_SV has month-to-month swings of ±0.04–0.06 at PGC, NLLB, LON, GNW. A season-preserving control
  (shuffle within calendar month) is the obvious next check for every result above. Not run.

## 11. Answers in short (as of 2026-10-14)
1. **Recovery at PGC:** CC_SV drops from the GPS midpoint, is lowest at +6 to +15 d (−0.06 vs long-term mean; −0.12 pre-vs-post
   including a +0.04 pre-ETS high), and is back to the long-term level by +24 to +27 d, i.e. **~3 weeks, recovered within a month**.
2. **Other stations:** on Vancouver Island, **NLLB (−0.089, pct 0.040) and YOUB (−0.205, pct 0.026)** show the same sign
   beyond the scrambled control (NLLB driven by a pre-ETS high; YOUB on 28 RFs); GOBB same sign, not significant; MGB none.
   In Washington, **LON (−0.118, pct 0.003)** has the same sign but sits off the tremor and looks slow/seasonal; **GNW, the only
   WA station over the tremor, shows the opposite sign** (+0.075); SP2, OFR, WISH, NLWA none (2007–2026, 15-day pre/post CC_SV).
3. **Timing:** at PGC the change goes with the main tremor passage (= midpoint +2.5 d), not with the first local tremor
   (~2 weeks earlier). The midpoint and the local front are too close at PGC to say which one it follows better.

---
# Added 2026-10-14 ~3:05 pm PT: same-month scrambled control, seasonal regression, pre-ETS window test

## 12. Same-month scrambled control and seasonal regression (PLAN Amendment 2, written before computing)
Script `scripts/rf07_samemonth.py` (existing per-RF tables only). Tables `tables/samemonth/samemonth_{bins,curve,season,lead}.csv`,
figures `figures/samemonth/recovery_<tag>_<STA>.png` (left: recovery curve with the old grey band and the same-month purple band;
right: monthly CC_SV and the seasonal fit, ETS months in blue). Same-month control: the ETS dates (not the RF values) are moved to
other years, same calendar date ±15 d, > 60 d from any real midpoint; 20,000 draws (5,000 for the curves). The observed Δ
reproduce the stored rf04 values to 1e-16. "SM median" = the median Δ of the same-month draws, i.e. what season and RF sampling
alone produce. pct = one-sided, CC down.

| station (run) | ΔCC_SV 15 d | SM median | SM σ | pct same-month | pct old | 30 d (pct SM) | 150 d (pct SM) | season-predicted Δ15 | season-removed Δ15 (pct SM) |
|---|---|---|---|---|---|---|---|---|---|
| **PGC (repro 1993–2019)** | **−0.122** | −0.001 | 0.041 | **0.003** | 0.002 | −0.070 (0.006) | −0.036 (0.004) | −0.008 | −0.114 (0.005) |
| **PGC (full 1993–2026)** | **−0.115** | −0.005 | 0.034 | **0.0004** | 0.0003 | −0.073 (0.002) | −0.030 (0.006) | −0.004 | −0.111 (0.0008) |
| PGC (2019–26 only) | −0.094 | +0.008 | 0.055 | 0.032 | 0.041 | −0.100 (0.006) | −0.024 (0.31) | −0.006 | −0.088 (0.037) |
| SNB (repro / full) | −0.059 / −0.029 | −0.002 | 0.075 | 0.23 / 0.36 | 0.27 / 0.38 | −0.120 (0.016) / −0.082 (0.071) | +0.009 / +0.004 | −0.006 / −0.001 | −0.053 / −0.027 |
| VGZ (repro / full) | −0.146 / −0.146 | −0.016 / −0.008 | 0.12 / 0.09 | 0.13 / 0.067 | 0.09 / 0.056 | −0.139 (0.056) / −0.134 (0.023) | −0.012 / −0.022 | +0.002 / −0.002 | −0.148 / −0.145 (0.067) |
| **NLLB** | **−0.089** | +0.004 | 0.053 | **0.044** | 0.040 | −0.075 (0.019) | −0.010 (0.41) | −0.012 | −0.077 (0.058) |
| **YOUB** | **−0.205** | −0.011 | 0.108 | **0.037** | 0.026 | −0.073 (0.11) | −0.010 (0.019*) | +0.011 | −0.216 (0.034) |
| GOBB | −0.109 | −0.007 | 0.086 | 0.11 | 0.15 | −0.050 (0.29) | +0.006 | −0.002 | −0.107 (0.11) |
| MGB | +0.013 | −0.008 | 0.067 | 0.63 | 0.59 | +0.005 | −0.005 | +0.004 | +0.009 |
| GNW (opposite sign) | +0.075 | +0.002 | 0.043 | 0.044 for CC up (two-sided 0.088) | 0.034 for up | +0.015 | −0.024 (0.045) | −0.005 | +0.080 (0.038 for up) |
| SP2 | +0.015 | −0.004 | 0.036 | 0.70 | 0.65 | −0.003 | −0.022 (0.17) | −0.002 | +0.017 |
| **LON** | **−0.118** | +0.006 | 0.042 | **0.0016** | 0.003 | −0.093 (0.0008) | **−0.034 (0.015)** | **+0.000** | **−0.118 (0.0014)** |
| OFR | −0.016 | +0.046 | 0.065 | 0.19 | 0.42 | −0.036 (0.025*) | −0.004 (0.0008*) | +0.020 | −0.036 (0.12) |
| WISH | +0.001 | −0.002 | 0.025 | 0.55 | 0.52 | −0.006 | +0.012 (0.013 for up) | −0.002 | +0.003 |
| NLWA | −0.010 | −0.001 | 0.034 | 0.39 | 0.39 | −0.012 | +0.005 | −0.002 | −0.007 |

\* pct driven by a shifted same-month median on a short record (YOUB 150 d median +0.060; OFR median +0.04–0.05 on 7 years
and 280 RFs): the control is thin there, not a result.

Seasonal regression (annual + semiannual, fitted on the RFs > 60 d from any midpoint; `samemonth_season.csv`):

| station | annual amp. | semiannual amp. | F-test p (4 terms) | CC_SV peak (day of year) |
|---|---|---|---|---|
| PGC (full) | 0.029 | 0.020 | 0.004 | ~211 (end July; high Aug–Sep, low Feb and Dec) |
| NLLB | 0.044 | 0.021 | 0.009 | ~203 |
| YOUB | 0.058 | 0.051 | 0.017 | ~271 |
| GNW | 0.050 | 0.008 | 0.010 | ~141 |
| NLWA | 0.031 | 0.021 | 0.016 | ~204 |
| VGZ, SNB, SP2, WISH, OFR | 0.02–0.04 | 0.004–0.044 | 0.06–0.67 | — |
| LON, GOBB, MGB | 0.005–0.018 | 0.005–0.009 | 0.89–0.99 | — |

What this shows:
- **PGC still passes (rule S1).** Its 15-day drop (−0.122 in 1993–2019, −0.115 in 1993–2026) is beyond 99.7 % and 99.96 % of the
  same-month draws. Season and RF sampling alone give ~0 (draw median −0.001 / −0.005). The 30–150-day bins also stay beyond the
  control. The 2019–26 out-of-sample run holds too (pct 0.032).
- **PGC's post-ETS dip survives; its pre-ETS high does not stand out.** Against the same-month band the recovery curve is below at
  +9 to +18 d (lowest −0.061 at +9 d). The +0.04–0.05 high at −18 to −9 d is now inside the band (old band: outside).
- **LON is not seasonal (rule S2: "surviving").** −0.118, same-month pct 0.0016. The seasonal model predicts +0.000 of it. With the
  season removed it is still −0.118 (pct 0.0014). LON's CC has no annual or semiannual cycle (amplitude 0.005, F p 0.99). The
  150-day persistence survives as well (−0.034, pct 0.015). On the same-month curve LON is above the band at almost every centre
  from −60 to −3 d (+0.03 to +0.10), then drops to −0.03 at the midpoint. In the same weeks the same-month draws dip (median down to
  −0.03), so season would push the other way. So LON's pattern is a ~2-month high before the WA midpoint that ends at it, and season
  does not explain it. What it is remains open. LON has no tremor within 30 km and its slab (~72 km) lies below the 0–60 km
  depth range used.
- **Other stations (rule S3).** NLLB (0.044) and YOUB (0.037) keep their reading, NLLB just inside the 5 % level. VGZ (0.067), GOBB,
  SNB, MGB, SP2, WISH and NLWA are unchanged: not significant. GNW's opposite-sign rise is at 4.4 % one-sided (8.8 % two-sided), so
  under the pre-set two-sided rule it does not hold. At 15 d, 4 of 13 stations pass in the PGC direction (PGC, NLLB, YOUB, LON). About
  0.7 would pass by chance.
- **The seasonal regression agrees.** CC_SV does have a seasonal cycle at 5 stations (annual 0.03–0.06, F p ≤ 0.02). At PGC the cycle
  is high in August–September, when SVI ETS cluster. But the cycle changes little across a 30-day pre/post pair: the season-predicted
  15-day Δ is −0.004 at PGC (3 % of the observed), −0.012 at NLLB (13 %), +0.011 at YOUB (opposite sign) and 0.000 at LON. A fit on
  all RFs gives the same (|Δ| ≤ 0.013).
- **Scope.** This tests the 15–150-day pre/post CC_SV contrast around GPS midpoints against a control that keeps the calendar date
  (±15 d), on the same RFs. It does not test other non-ETS causes tied to the same dates.

## 13. Pre-ETS window on stations not used to see it (acf_pgc PLAN §10, Amendment C2, pre-registered)
H = mean CC_SV in days [−21, 0) minus mean in [−60, −25), same-month control (20,000 draws), pct = one-sided for H > 0.

| station | role | H | SM median | SM σ | pct (same-month) | RFs in the two windows |
|---|---|---|---|---|---|---|
| GOBB | primary | +0.011 | +0.001 | 0.057 | 0.43 | 21 / 22 |
| MGB | primary | −0.054 | +0.012 | 0.051 | 0.89 | 43 / 59 |
| GNW | primary | −0.070 | +0.020 | 0.034 | 0.99 | 102 / 151 |
| SP2 | primary | +0.031 | +0.003 | 0.028 | 0.15 | 87 / 142 |
| OFR / WISH / NLWA | secondary (no local tremor) | +0.056 / −0.009 / −0.020 | | | 0.07 / 0.78 / 0.92 | |
| PGC (full; where the lead was seen) | reference | +0.024 | +0.010 | 0.026 | 0.29 | 74 / 118 |
| NLLB, SNB repro, LON (seen) | reference | +0.122 / +0.098 / +0.018 | | | 0.001 / 0.043 / 0.081 | |

**The pre-ETS high is not confirmed on the stations not used to see it.** H > 0 at 2 of 4, and the Stouffer-combined p is 0.90.
The pass rule needed p ≤ 0.05 and 3 of 4 → **fail**. GNW goes the other way (−0.070), and its 15-day rise in §8 points the same way.
In this form (pre-ETS weeks vs 4–8 weeks before) PGC's own H is only +0.024 (pct 0.29). The PGC "high" in §6 is relative to
the long-term mean, and the weeks −60…−25 are high too. NLLB is the one station with a strong H (+0.12, pct 0.001), but NLLB was
among those used to see the lead. Combined reading with the ACF test: `../acf_pgc/RESULTS.md`, "Lead test".
