# P4-0 / P4-0a — Puget–Olympic station shortlist: the site proxy CALIBRATES, and all three sub-regions STOP

*2026-09-20. Everything here is `u_`-prefixed: `scripts/u_*.py`, `results/p4/`, `figures/u_*.png`.
Nothing outside `lfe_picker/` was written; `../lfe_extension/traces` and `../catalogs/` were read only.*

## 0. Headline

1. **The proxy calibrates.** On the 14 north stations with a measured single-event registration
   ρ (P3D §3.2 + §A3.2), the pre-registered model ρ/ρ_PGC ~ (distance, quiet-night noise vs PGC,
   tremor/quiet contrast) gives **leave-one-out Spearman +0.688 (p = 0.007, n = 14), all three signs
   as expected, LOO median |error| 0.166 in ρ/ρ_PGC** — inside Merlin's ≥0.6 / ≤0.25 acceptance.
   So Puget selection is **calibrated**, not uncalibrated.
2. **But the predictive content is mostly geometry.** Distance alone reaches LOO Spearman +0.732;
   the fitted contrast coefficient is ≈0 (+0.0001 per dB). The noise term contributes
   −0.0126 per dB (10 dB of extra site noise ⇒ −0.13 in ρ/ρ_PGC). The contrast is a *univariately*
   good proxy (Spearman +0.556 in the 1–8 Hz cache band, **+0.754 in the native 2–8 Hz band**) but
   it is nearly collinear with distance, so it adds little once distance is in. Reported plainly
   because the headline number alone would overstate what the site measurement buys.
3. **All three sub-regions STOP**, under the pre-registered rule *and* under every relaxation
   tested. The binding clause is not the noise bars — it is
   *"≥4 admitted sites within 45 km hypocentral of ≥50 % of the sub-region's tremor-active cells"*.
   Admitted sites: R1 3, R2 3, R3 3 (specified far bar) or 5 / 4 / 6 (intent far bar). The fraction
   of cells with ≥4 admitted sites inside 45 km is **0.00 everywhere**, and 0.00–0.07 even if the
   clause is loosened to ≥3 sites.
4. **The geometry itself is not the problem; the alignment between geometry and quality is.**
   With ALL 118 candidate sites the clause is satisfiable — R1 0.55, R2 0.60, R3 0.34 of cells have
   ≥4 sites within 45 km. The sites that supply that coverage are the noisy UW lowland broadbands
   that fail the tremor-contrast bar; the sites that pass the bars are the quiet boreholes and
   CN/PQ stations sitting 45–70 km away. **Puget has the network Bostock had, but not at the
   stations that register tremor.**

**Verdict: STOP P4** under the pre-registered P4-0a rule. What would change it is stated in §5.

## 1. What was run, and what it cost

| stage | script | station-days | downloaded |
|---|---|---|---|
| P4-0 probe (139 Puget candidates × 8 days) | `u_select`,`u_download`,`u_cache`,`u_stats`,`u_inventory` | 969 planned, 928 with statistics | 21.8 GiB (95 linked from `../lfe_extension/traces`) |
| P4-0a north calibration (15 stations × 8 days) | `u_north`, `u_hf_north` | 120 planned, 116 delivered | 3.1 GiB (19 linked) |
| P4-0a far-tremor control (134 stations × 4 days) | `u_fartremor` | 498 planned, 468 delivered | 11.1 GiB (48 linked) |
| P4-0a native-rate 10–19 Hz control | `u_fartremor` (HF pass) | 928 + 116 station-days re-read at native rate | — |

**36 GiB total, ~25 min wall.** Everything went through the same chain: `r_download_day.fetch_day`
→ `o_cache_month.cache_one` (frozen 20 Hz / 1–8 Hz; a documented rational-decimation twin for the
22 candidates at 50/125 Hz) → `u_stats.one` (velocity PSD via StationXML evalresp; incoherent 6-s
window power). Channel dip was read from StationXML for every channel and applied on read
(`zsign = −1` for the 4 stations declaring the vertical at dip +90).

## 2. Calibration on the north (Merlin F5)

Probe days, chosen inside S18/S21/S22 with ≥100 PNSN epicentres within 40 km of the Bostock family
centroid (48.50 N, −123.78, 35.3 km), and matched quiet days with **no** epicentre within 150 km on
the day or ±1 d. The same 8 days for all 15 stations, so the PGC ratio is exact.

| tremor day | n(40 km) | matched quiet day |
|---|---|---|
| 2024-08-21 | 502 | 2024-08-06 |
| 2024-08-25 | 478 | 2024-09-07 |
| 2025-10-06 | 474 | 2025-10-25 |
| 2025-10-07 | 465 | 2025-10-26 |

(S22 contributed no day with ≥100 epicentres within 40 km of the centroid, so the design is 4+4 —
the same shape as the Puget probe, not the 6+6 Merlin sketched.)

### 2.1 The table

`noise` = median 2–8 Hz velocity PSD on the QUIET days in dB relative to CN.PGC on the same days;
`contrast` = 10 log10 (tremor / quiet) of the median 6-s window power; `dist` = mean hypocentral
distance to the 130 Bostock family locations.

| station | ρ/ρ_PGC | dist km | noise Z dB | noise H dB | noise 3-C dB | contrast 3-C dB (1–8 Hz) | contrast 3-C dB (native 2–8 Hz) |
|---|---|---|---|---|---|---|---|
| PQ.LHLYB | 1.086 | 46.4 | −3.6 | −4.7 | −4.8 | +6.5 | +18.4 |
| CN.PGC | 1.000 | 53.2 | 0.0 | 0.0 | 0.0 | +5.5 | +13.9 |
| PQ.DAOB | 0.826 | 51.6 | −1.6 | −0.4 | −1.3 | +3.6 | +14.1 |
| CN.SYMB | 0.783 | 46.8 | −8.4 | −7.3 | −7.6 | +8.2 | +21.7 |
| CN.CLRS | 0.773 | 63.1 | −8.8 | −9.9 | −9.9 | +4.8 | +14.5 |
| CN.GOBB | 0.571 | 68.7 | +4.7 | +4.1 | +3.5 | −2.0 | +3.8 |
| UW.TURTL | 0.559 | 70.8 | −7.3 | −8.6 | −8.5 | −3.1 | +8.8 |
| UW.RNWY | 0.488 | 64.8 | +14.4 | +14.4 | +14.3 | +5.7 | +7.0 |
| UW.LOPEZ | 0.398 | 76.7 | −4.2 | −4.7 | −4.9 | −8.4 | +5.2 |
| UW.LRIV | 0.366 | 68.1 | −2.0 | +1.4 | +1.5 | +4.7 | +13.7 |
| UW.BUCKS | 0.296 | 59.4 | +19.4 | +16.4 | +16.7 | −1.4 | +3.2 |
| CN.PTRF | 0.275 | 66.6 | +3.8 | +2.1 | +2.0 | −2.0 | +3.3 |
| UW.OSQM | 0.168 | 80.2 | −4.4 | +1.9 | +1.9 | −2.5 | +6.6 |
| UW.ODUC | 0.158 | 74.6 | −5.3 | −2.3 | −2.4 | +1.2 | +4.5 |
| *CN.VGZ (derived ρ, held out)* | *0.590* | *55.3* | *+3.8* | *+2.7* | *+2.6* | *+0.2* | — |

VGZ is not tabulated against PGC in §3.2; its ρ/ρ_PGC ≈ 0.59 is recovered from the ρ/ρ_VGZ column
and is **held out of the fit**. As an external check the fitted model predicts 0.704 for it
(error 0.11, inside the LOO median error).

### 2.2 Spearman, signs, and leave-one-out

| proxy (n = 14) | Spearman vs ρ/ρ_PGC | p | expected sign? |
|---|---|---|---|
| mean hypocentral distance to the family field | **−0.745** | 0.002 | yes (negative) |
| quiet-night noise vs PGC, 3-C | −0.385 | 0.175 | yes (negative) |
| quiet-night noise vs PGC, horizontals | −0.385 | 0.175 | yes |
| quiet-night noise vs PGC, vertical | −0.077 | 0.794 | yes, but ≈0 |
| tremor/quiet contrast, 1–8 Hz cache band | +0.556 | 0.039 | yes (positive) |
| **tremor/quiet contrast, native 2–8 Hz** | **+0.754** | 0.002 | yes (positive) |
| distance to the NEAREST family | +0.134 | 0.648 | no — not usable |

| model (LOO) | Spearman | median &#124;err&#124; | β |
|---|---|---|---|
| **prescribed** dist + noise_3C + contrast_3C | **+0.688** (p 0.007) | **0.166** | 1.943, −0.0218, −0.0126, +0.0001 |
| dist + noise + contrast (native 2–8 Hz) | +0.684 | 0.187 | 1.839, −0.0207, −0.0115, +0.0030 |
| dist + noise_H + contrast_H | +0.710 | 0.169 | — |
| dist + contrast | +0.701 | 0.196 | — |
| **dist only** | **+0.732** | 0.181 | — |
| contrast only | +0.385 | 0.234 | — |
| noise only | +0.099 | 0.240 | — |

**ACCEPT** on the prescribed model (≥0.6, expected signs, ≤0.25). The honest caveat is in §0.2:
`dist only` is as good, so the *site measurement* is worth ~nothing beyond geometry in this fit,
and the ranking below is, in practice, a geometry ranking with a small noise correction.
Two further caveats: ρ is measured against a single northern source field, and the calibration set
spans only 46–80 km, so extrapolation beyond ~85 km is unsupported. `figures/u_rho_calibration.png`.

## 3. The far-tremor control (Merlin's ≤ +1 dB bar)

### 3.1 The specified construction is geometrically impossible

Reusing north tremor days on which Puget had no tremor cannot be done: the Bostock family field
(48.50 N) sits **~25 km north of the top of the P4 band**, so every day with ≥100 north epicentres
has tremor within 50 km of a Puget active cell. Checked day by day over 2023-01-01 → 2026-05-27:
**zero** qualifying days (also zero at a 100 km buffer).

Substituted, at Merlin's own 100–200 km scale: days with ≥100 PNSN epicentres **120–350 km** from
the nearest Puget active cell and **not one** within 120 km — both azimuths, so a station cannot
pass by being blind to one direction.

| far-tremor day | n | median distance | source |
|---|---|---|---|
| 2024-06-24 | 277 | 171 km | north Vancouver Island (49.8 N) |
| 2025-06-08 | 244 | 170 km | north Vancouver Island |
| 2024-06-11 | 185 | 212 km | Oregon (45.9 N) |
| 2025-01-25 | 258 | 285 km | Oregon (45.6 N) |

### 3.2 The absolute +1 dB bar is a factor 3 inside its own noise

The same statistic measured on the 4 *quiet* days scatters by a **median 3.4 dB** (max deviation
from the quiet median; p90 6.3 dB; median max/min spread 5.1 dB) across 110 stations. A far-tremor
"contrast" of +1 dB is therefore not distinguishable from ordinary day-to-day variability, and the
far days are also seasonally mismatched to the quiet days. Both readings are carried:

* **specified (primary):** far contrast ≤ +1 dB — 69 / 134 stations pass;
* **intent:** local − far ≥ 3 dB (*"is the contrast as large for distant tremor as for local?"*) —
  30 / 134 pass, and it is the reading that isolates the real failures: **PB.B003 far +21.4 dB**
  (a glitch day: its 2025-05-16 window power is 10³ above its own normal), **UW.LOPEZ far 2.5 dB
  vs local 0.0**, **UW.OSQM far 1.1 vs local 3.3**.

Both verdicts are reported throughout; neither changes the STOP.

### 3.3 The 10–19 Hz control

Measured at the **native** rate (the 20 Hz cache is band-limited to 1–8 Hz), four 1.5-h blocks per
day. **106 of 133** stations have a smaller contrast in 10–19 Hz than in the native 2–8 Hz band, as
a low-frequency source requires. The same pass exposed that the cache-band statistic is really a
1–8 Hz number diluted by the 1–2 Hz microseism: native 2–8 Hz contrasts run 3–8 dB larger
(SYMB 8.2 → 21.7 dB, LRIV 4.7 → 13.7 dB), and correlate with ρ far better (+0.754 vs +0.556).
**The per-day ≥ +3 dB admission bar therefore uses the native 2–8 Hz statistic**, which is the band
Merlin specified; the cache-band contrast remains the calibration variable because that is the band
the north fit was built in.

## 4. The Puget rank and the sub-region decision

Predicted ρ/ρ_PGC applies the north model per sub-region, with the distance term built the same way
as in the north (tremor-weighted mean hypocentral distance to that sub-region's active cells; the
model is linear, so this equals the mean of the per-cell predictions).

### 4.1 Admission funnel (134 probed stations)

| bar | passing |
|---|---|
| predicted ρ/ρ_PGC ≥ 0.5 | 35 |
| tremor contrast ≥ +3 dB on ≥3 of 4 tremor days (native 2–8 Hz) | 31 |
| *(same bar on the diluted 1–8 Hz cache band, for reference)* | *20* |
| far-tremor ≤ +1 dB (specified) | 69 |
| local − far ≥ 3 dB (intent) | 30 |
| 10–19 Hz contrast < 2–8 Hz contrast | 106 |
| delivery ≥ 90 % of the days in the station's own epoch, coverage ≥ 0.9 | 117 |
| had ≥3 tremor days in its epoch at all | 116 |
| **ALL (specified far bar)** | **12 stations / 9 unique** |
| **ALL (intent far bar)** | **23 stations / 15 unique** |

### 4.2 Sub-regions

| rule | sub-region | admitted sites | withheld | cells with ≥4 admitted sites ≤45 km | ≥3 ≤45 km | ≥4 ≤60 km | median gap | verdict |
|---|---|---|---|---|---|---|---|---|
| specified | R1 47.0–47.5 | 2 (+1 withheld) | HERD | 0.00 | 0.00 | 0.00 | 360° | **STOP** |
| specified | R3 47.5–48.0 | 2 (+1) | SKOKO | 0.00 | 0.00 | 0.00 | 360° | **STOP** |
| specified | R2 48.0–48.3 | 2 (+1) | SYMB | 0.00 | 0.01 | 0.00 | 308° | **STOP** |
| intent | R1 | 3 (+2) | HDW, KTSAP | 0.00 | 0.07 | 0.36 | 297° | **STOP** |
| intent | R3 | 4 (+2) | HDW, LRIV-cluster | 0.01 | 0.05 | 0.44 | 309° | **STOP** |
| intent | R2 | 2 (+2) | PGC-cluster, SYMB | 0.00 | 0.00 | 0.06 | 360° | **STOP** |

Admitted / withheld stations (intent far bar, the more permissive of the two):
R1 SKOKO, PAN4H, KTSAP(w), HDW(w), HERD · R2 B005/B006/B007, B001, LRIV, SYMB(w), B009/B010/B011(w)
· R3 B943, B001, HDW(w), B005/B006/B007(w), KTSAP, SKOKO, LRIV(w).
Note the 1-km site merge is what hurts R2: **B005+B006+B007+LRIV are ONE site**, and so are
PGC+B009+B010+B011 and OSQM+B001 — six admitted *stations* collapse to two *sites*.

Every long-record PB borehole in the box is carried as **backfill only** where it fails a bar:
B001, B003, B004, B013, B014, B017, B018, B941, B943 (B005/B006/B007 are admitted; B941/B014 fail
predicted ρ; B003 fails the far control on a glitch day).

### 4.3 Why — and what would change it

The clause is unattainable at this *number* of sites, not at this *quality*:

| site set | R1 | R3 | R2 | (fraction of cells with ≥4 sites ≤45 km) |
|---|---|---|---|---|
| all 118 candidate sites | 0.55 | 0.34 | 0.60 | the geometry ceiling |
| 106 eligible probed sites | 0.53 | 0.31 | 0.57 | |
| top-12 by predicted ρ | 0.26 | 0.10 | 0.15 | |
| top-20 by predicted ρ | 0.40 | 0.20 | 0.17 | |
| top-30 by predicted ρ | 0.51 | 0.31 | 0.54 | |
| admitted (intent bar) | 0.00 | 0.01 | 0.00 | |

So ~30 sites per sub-region are needed to satisfy the clause, and the bars admit 2–6. The
geometrically decisive sites and the bar they fail (R1, ranked by cells within 45 km):
WINDI (62 cells) — contrast 2/4 days; PAN4H (57) — far +2.8 dB; **SKOKO (56) — admitted**;
ANGUS (46) — contrast 1/4; **HERD (45) — admitted**; WYNO (42, contrast +11.7 dB, the strongest in
R1) — delivery (2 of its 8 days returned no data); EAGLE (36) / DOORS (36) — contrast; MONTE, B017,
WISH — predicted ρ. In R2: BUCKS (71 cells) — contrast −4.8 dB; B003 (48) — far control;
SLDQ (48) — predicted ρ 0.48, just under the bar.

Three things would change the verdict, in order of cost:
1. **Re-read the clause at 60 km** (R1 0.36, R3 0.44 under the intent bar) — still below 0.5, so
   this alone does not rescue it;
2. **More probe days.** 4 tremor days makes "≥3 of 4" a 75 % clause on a noisy per-day statistic;
   Merlin specified 6+6 (≥4 of 6 = 67 %). WINDI, DOORS, ANGUS and EAGLE sit at 1–2 of 4. A 6+6
   re-probe of the ~25 geometrically decisive sites is ~3 GiB and ~10 min and is the cheapest
   test of whether the contrast bar or the station is the problem;
3. **Accept the pivot Merlin himself named** — growth from *transferred* northern-family-class
   templates placed at Puget cells (moveout-only seeds), which does not need a discovery-grade
   local network.

## 5. The two judgement calls (both turn out to be moot)

* **PSD bar, power vs amplitude.** Resolved as instructed: the **power reading (≤ +4.8 dB) is
  primary**; 80 / 134 stations pass it, 97 / 134 pass the amplitude reading (≤ +9.5 dB). It changes
  nothing: **every** admitted station passes the power reading (9/9 specified, 15/15 intent), so the
  PSD clause is not binding anywhere in v1. It is reported as a column, not applied as a separate
  bar, because the same noise measurement already enters predicted ρ.
* **`chain = rational` stations (50/125 Hz, rational rather than integer decimation).** **None of
  them is admitted** under either rule, so no frozen-chain equivalence check is owed before any
  template product. Should a later re-probe admit one, the check is still required and the
  `chain` column flags them.

## 6. Faults found and fixed while doing this

1. **`loc` parsed as a float.** `pd.read_csv(..., na_values=[""])` on `candidates.csv` turned
   location code `"00"` into `0.0`, so the response lookup missed for 11 stations × 4 far days and
   their PSD/energy came back in **raw counts** — producing nonsense "far contrasts" of 165 dB.
   Fixed with `dtype={"loc": str}`; all 468 far station-days now carry an evalresp response. The
   P4-0 probe statistics were never affected (0 `no_response` rows there).
2. **Delivery denominator.** Delivery was being computed over all 8 probe days, which penalises a
   station for days before it was installed. Now `n_days_ok / n_days_planned` (the days the station
   had an epoch for): 89 → 117 of 134 pass.
3. **Contrast band.** The "2–8 Hz" contrast taken from the 20 Hz cache is really 1–8 Hz and is
   diluted; the admission bar now uses the native-rate 2–8 Hz statistic (§3.3).

## 7. Files

`results/p4/`: `north_stations.csv`, `north_probe_days.csv`, `north_station_day_stats.csv`,
`north_hf_station_day_stats.csv`, `north_calibration.csv`, `rho_calibration.json`,
`far_tremor_days.csv`, `far_station_day_stats.csv`, `hf_station_day_stats.csv`,
`shortlist_v1.csv` (specified far bar), `shortlist_v1_margin.csv` (intent far bar),
`subregion_summary.csv`, `subregion_summary_margin.csv`, plus the P4-0 products
(`candidates.csv`, `probe_days.csv`, `station_probe_days.csv`, `noise_inventory.csv`,
`shortlist_v0.csv`, `station_rank_all.csv`, `shortlist_cell_geometry.csv`).
`figures/`: `u_noise_inventory.png`, `u_rho_calibration.png`, `u_shortlist_v1.png`.
`scripts/`: `u_common.py`, `u_select.py`, `u_download.py`, `u_cache.py`, `u_stats.py`,
`u_inventory.py`, `u_north.py`, `u_calibrate.py`, `u_fartremor.py`, `u_hf_north.py`,
`u_shortlist_v1.py`, `u_p4_0_queue.sh`.

---

# ADDENDUM — 2026-09-20: the 6+6 re-probe at Merlin's original clause. STOP stands.

The first pass evaluated the contrast bar on **4 tremor days at "≥3 of 4"**; Merlin specified
**6 tremor + 6 matched quiet days at "≥4 of 6"**. That shortcut was ours, so the clause was
restored and re-run. **Pre-registered before the first byte was fetched**
(`results/p4/reprobe_prereg.json`, 2026-09-20T04:15:49Z: clause, day-selection rule, station list,
and SHA-256 of 11 input files + 7 scripts).

## A.1 What was added

Two extra tremor days chosen on **PNSN epicentre count alone** (≥100 within 40 km of a 47.0–48.3 N
active cell, ≥3 days from any day already used, one from 2023-03/04 and one from 2024-08 / 2025-09-10),
and two matched quiet days by the unchanged quiet rule (no epicentre within 75 km on the day or ±3 d):

| new tremor day | epicentres ≤40 km | matched quiet day |
|---|---|---|
| 2023-04-11 | 415 | 2023-03-19 |
| 2024-08-17 | 671 | 2024-08-03 |

Station set as specified — every candidate **site** contributing to ≥10 cells of any sub-region at
≤45 km hypocentral (50 sites), plus every site then scoring 1–2 of 4 on the native contrast
(62 sites): **101 stations / 84 sites**, larger than the ~25–35 estimate. 338 station-days queued,
**323 delivered (10 GiB), 2.5 min wall**. Far-tremor days unchanged. Calibration unchanged.
A station that was not re-probed keeps its 4-day measurement and is judged on "≥3 of 4"; the
`contrast_need` column records which clause each row was judged under.

## A.2 What changed — the bars

| bar (134 probed stations) | 4-day pass | 6-day pass |
|---|---|---|
| tremor contrast (native 2–8 Hz) | 31 | **59** |
| far-tremor ≤ +1 dB (specified) | 69 | **50** |
| local − far ≥ 3 dB (intent) | 30 | **41** |
| 10–19 Hz < 2–8 Hz | 106 | 109 |
| delivery ≥ 90 % | 117 | 121 |
| predicted ρ ≥ 0.5 | 35 | 35 (unchanged by construction) |

**The 4-day/≥3-of-4 reading was indeed over-strict: the contrast bar's pass count nearly doubled.**

**And the re-probe settles the far-bar question empirically.** Six stations that passed the
absolute ≤ +1 dB bar on 4 quiet days fail it on 6 — B005 (1.58 dB), B006 (1.42), LRIV (1.16),
SKOKO (1.75), HDW (1.12), HERD (2.58) — while their **local − far margins are 6.5 to 10.4 dB**,
i.e. every one of them is an unambiguous local-tremor responder. Merely adding two quiet days to
the baseline flipped them. That is direct evidence that the absolute +1 dB clause measures baseline
noise, not far-tremor response (§3.2: the quiet-day-to-quiet-day scatter of the same statistic is
median 3.4 dB). Under the intent bar **nothing was lost and 12 stations were gained**.

## A.3 The four verdicts side by side

| | rule | R1 47.0–47.5 | R3 47.5–48.0 | R2 48.0–48.3 |
|---|---|---|---|---|
| **4-day** | specified far bar | 3 sites (2+1) · n45≥4 **0.00** · STOP | 3 (2+1) · **0.00** · STOP | 3 (2+1) · **0.00** · STOP |
| **4-day** | intent far bar | 5 (3+2) · **0.00** (≥3: 0.07) · STOP | 6 (4+2) · **0.00** (≥3: 0.01) · STOP | 4 (2+2) · **0.00** · STOP |
| **6-day** | specified far bar | **0 sites** · 0.00 · STOP | 2 (2+0) · 0.00 · STOP | 5 (3+2) · 0.00 (≥3: 0.01) · STOP |
| **6-day** | intent far bar | **8 (6+2)** · 0.024 (≥3: **0.218**) · STOP | **10 (8+2)** · 0.083 (≥3: 0.122) · STOP | **7 (5+2)** · 0.000 (≥3: 0.014) · STOP |

(`n45≥4` = fraction of the sub-region's tremor-active cells with ≥4 admitted, non-withheld sites
inside 45 km hypocentral; the clause needs ≥0.50. Median azimuthal gap improves to 233° in R1 and
291° in R3 under the best configuration, still far from closed.)

**Admitted / withheld (6-day, intent far bar — the most permissive defensible configuration):**
* **R1** WINDI, SKOKO, GNW, DOORS, HDW, HERD + withheld PAN4H, KTSAP
* **R3** B943, B001, OSQM, GNW, KTSAP, SKOKO, COYL, DOSE, EQUIL + withheld HDW, LRIV-cluster
* **R2** B006/B007/B005-cluster, B001, LRIV, OSQM, DAOB, SYMB, B009/B010/B011 + withheld ALBHB, VGZ

Stations that flipped **in** with the longer probe (all under the intent bar; none flipped out):
R1 **WINDI, GNW, DOORS**; R3 **COYL, DOSE, EQUIL, GNW, OSQM**; R2 **ALBHB, DAOB, OSQM, VGZ**.
Under the specified bar the flips are the far-bar artefact described in A.2:
R1 lost HDW, HERD, SKOKO (→ 0 sites); R3 lost B005, B006, HDW, LRIV, SKOKO and gained EQUIL;
R2 lost B005, B006, LRIV and gained DAOB, VGZ.

## A.4 Conclusion

**STOP stands for all three sub-regions, at Merlin's original clause, under both far-bar readings.**
The re-probe did what it was supposed to do — it showed the 4-day contrast clause was too strict
(31 → 59 stations) and nearly doubled the admitted sets under the intent bar — and the sub-region
rule still fails by an order of magnitude: the best configuration reaches **0.02–0.08** of cells
with ≥4 admitted sites inside 45 km against a bar of **0.50**, and **0.22** even at ≥3 sites.
This is the §4.3 arithmetic, now confirmed with twice the evidence per station: with 6–10 admitted
sites per sub-region the clause is unreachable, because it needs ~30. The limit is the number of
Puget sites that both sit near the interface and register tremor, not the length of the probe.

Files: `reprobe_prereg.json`, `probe_days_v2.csv`, `probe_queue_v2.csv`, `station_day_stats_v2.csv`,
`hf_station_day_stats_v2.csv`, `reprobe_download_manifest.csv`;
`shortlist_v1.csv` / `shortlist_v1_margin.csv` are now the **6-day** products
(`shortlist_v1_6day*.csv`), with the first pass preserved as `shortlist_v1_4day*.csv`;
`subregion_summary{,_margin,_4day,_4day_margin,_6day,_6day_margin}.csv`;
`figures/u_shortlist_v1_6day.png`. Script `scripts/u_reprobe.py`.

---

# Corrections after Merlin close-out (2026-09-20)

Source: `results/merlin_review_p4_closeout_2026-09-20.md` (F2, F7, F8, Q1). The **STOP verdict for
all three sub-regions stands and is over-determined** — nothing below reverses it. The original text
above is left unedited; these are corrections to its *reasoning*, not to its conclusion.

### C-1 (F2) — the "≥ 50 % of cells" clause was unattainable, so it could not have been a
### discovery-feasibility gate

§0.3 treats "≥ 4 admitted sites within 45 km hypocentral of ≥ 50 % of the sub-region's
tremor-active cells" as the binding clause. Re-measured from `results/p4/shortlist_cell_geometry.csv`,
set `all_candidates` — **every one of the 118 candidate sites admitted, no quality bar applied at
all**:

| sub-region | active cells | frac with ≥ 4 sites ≤ 45 km | frac with ≥ 3 |
|---|---|---|---|
| R1 (47.0–47.5 N) | 165 | **0.545** | 0.685 |
| R3 (47.5–48.0 N) | 156 | **0.340** | 0.449 |
| R2 (48.0–48.3 N) | 142 | **0.599** | 0.768 |

**R3 could never have passed the 0.50 bar under any selection**, and R1/R2 clear it only if
essentially nothing is rejected. The clause was a bound-coverage criterion applied as a
discovery-feasibility gate; that was a specification error (Merlin's, acknowledged in his F2), and
§0.3's framing of it as "the binding clause" overstates what failing it means. §0.4's own sentence —
"with ALL 118 candidate sites the clause is satisfiable — R1 0.55, R2 0.60, R3 0.34" — already
contains the correction; it should have been read as evidence that the bar was mis-set, not only as
evidence about station quality.

**The geometry finding that stands** is the *admitted-set* coverage, which needs no arbitrary bar to
interpret (`results/p4/subregion_summary_6day_margin.csv`, intent far bar, 6-day probe):

| sub-region | admitted sites (summing + withheld) | frac cells ≥ 4 sites ≤ 45 km | frac cells ≥ 3 | median n45 | median azimuthal gap |
|---|---|---|---|---|---|
| R1 | 8 (6 + 2) | 0.024 | 0.218 | 1.0 | 233° |
| R3 | 10 (8 + 2) | 0.083 | 0.122 | 0.0 | 291° |
| R2 | 7 (5 + 2) | 0.000 | 0.014 | 0.0 | 330° |

That is the durable statement: **6–10 usable sites per sub-region, sitting 50–65 km from the tremor
field, covering at most ~22 % of active cells at three sites and ~8 % at four**, with median
azimuthal gaps of 233–330°. Bostock's families were built with typically five sites ≤ 45 km.

### C-2 (Q1) — the far-tremor bar

**The intent substitution (local − far ≥ 3 dB) is the operative bar**, not the absolute ≤ +1 dB
clause. Merlin accepts it: the absolute clause sits inside the 3.4 dB quiet-day scatter of the same
statistic and flipped six unambiguous local responders (B005, B006, LRIV, SKOKO, HDW, HERD; margins
6.5–10.4 dB) when two quiet days were added. Both readings were carried in the appendix and neither
changes the outcome; the report should have led with the differential reading. All admitted-site
counts quoted in C-1 are under the intent bar.

### C-3 (F8) — "the proxy calibrates" should not have been the headline

The LOO Spearman +0.688 is real, but §0.2 already records that **distance alone reaches +0.732** and
that the fitted tremor-contrast coefficient is ≈ 0. The shortlist is therefore **a geometry ranking
with a small site-noise correction** (−0.0126 in ρ/ρ_PGC per dB), not a calibrated multi-term site
proxy. It is acceptable as used — the contrast bar is applied separately as a hard filter — but the
ordering of §0.1 and §0.2 inverts the honest emphasis.

### C-4 (F7) — withheld selection

Withholding by ρ rank put **PAN4H** (third-best geometry in R1, 57 cells) into R1's withheld pair.
Withheld sites should be chosen to minimise coverage loss, not by rank. Moot: no sub-region
proceeds.

### C-5 (Q1) — a cell-restricted attempt is legitimate in principle, and moot in fact

Merlin notes that the natural unit for a discovery attempt is a contiguous cluster of tremor-active
cells with ≥ 4 admitted sites ≤ 45 km plus 2 withheld ≤ 60 km, and that R1 would offer 4 such cells
(36 at ≥ 3 sites) and R3 13 (19 at ≥ 3). That variant is not pursued: P4-0b (below, and
`results/P4_0B_RESULTS.md`) shows the construction has **no template-free entry point** — the seed
rule selects catalogued LFEs at a rate flat in z — and Puget's BUILD period is ~67 tremor days
against the 203 episode days on which the north arm already starved. Legitimate, moot.
