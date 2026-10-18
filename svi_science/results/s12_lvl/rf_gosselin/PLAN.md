# S12 / rf_gosselin — reproduce Gosselin et al. 2020 (Sci. Adv. 6:eaay5174), then extend

Written 2026-10-13 (PT) before any data were touched. Amendments go at the bottom, dated, never edits above.

## 0. Why
Gosselin et al. (2020) found, with receiver functions (RFs), that the shear-wave velocity of the Cascadia
low-velocity layer (LVL, ~35 km under southern Vancouver Island) rises by ~0.05–0.1 km/s for about two weeks
after each ETS midpoint (fault-valve: pore pressure drops after slip). The user is clear the signal is real.
We (1) reproduce it on their stations and period, (2) extend it (2019–2026 out-of-sample, full 1993–2026,
WA stations, finer time resolution) and (3) turn it into a calibrated target (window, sign, size) for the LFE
coda dv/v search.

## 1. Their recipe (from Materials and Methods) and our implementation
Their group's code is `rfpy` (Audet); we installed rfpy 0.1.2 and copy its conventions where the paper is silent.

| Step | Paper | Ours |
|---|---|---|
| Stations | CN PGC, SNB, VGZ, LZB, PFB | same (IRIS/EarthScope FDSN copy of CN data) |
| Events | M > 5.8, 30–90°, 1993–2019, SNR > 5 dB on vertical | USGS ComCat M ≥ 5.8, 30–90°, 1993–2026; SNR = rfpy `calc_snr` (0.05–1 Hz, 30 s after vs before P) > 5 dB |
| Window | 120 s P coda | download P−150 s … P+150 s; rfpy layout (noise = P−120…P−5, signal = P−5…P+110) |
| Pre-processing | decimate to 5 Hz, low-pass 1 Hz | low-pass 1 Hz (2-pole, zero-phase), resample 5 Hz; response removed to velocity (instrument changes over 26 yr) |
| Rotation | Z-R-T then up-going P, SV, SH (Bostock 1998) | rfpy `rotate(align='PVH')`, surface vp 6.0, vs 3.5 km/s |
| Deconvolution | Wiener spectral, noise from pre-event (Audet 2010) | rfpy Wiener: S·P*/(P·P* + N·N*) |
| RF filter | 0.05–0.5 Hz; variant up to 0.7 Hz | both: 0.05–0.5 Hz (main) and 0.05–0.7 Hz (variant) |
| Depth | time→depth with a southern Vancouver Island profile (Ramachandran 2006), 0–60 km | 1-D Ps moveout with a layered SVI average (Vp 6.0/6.4/6.8/7.6 at 0/10/25/40 km, Vp/Vs 1.75), dz = 0.25 km |
| Harmonics | k = 0,1,2 jointly on SV and SH (Eq. 1), SVD, all data | same, per station, all RFs of the analysed period |
| CC | observed vs harmonic-predicted RF, per event, SV and SH separately; outliers MAD > 2.5 removed | same (0–60 km) |
| dVs | Ps amplitude of LVL top and bottom → Δβ/β with the Aki–Richards linear P→S coefficient at slab-normal incidence (slab strike/dip at station), β = 4.6 km/s; cubic-in-time trend removed; MAD 2.5 outliers | same; LVL top/bottom depths picked once per station from the k = 0 harmonic (negative then positive pulse in the slab window), slab geometry from McCrory 2012 / Slab2 in repo |
| ETS times | ALBH GPS east, NAM08, detrend, 14-day running mean, 5th-order polynomial per episode, midpoint = 2nd-derivative zero | UNR IGS20 ALBH tenv3 (1994–2026), same steps; cross-checked against repo lists (tremor-based SVI episodes `b_ets_record/episodes_b.csv`, `s1/episode_geodetic.csv`, `data/ets_cycles_pgc.csv`) |
| Bins | pre = [mid−W, mid), post = [mid, mid+W), W = 15, 30, 45, 60, 90, 120, 150 d; two-sample t test | same |
| Bootstrap | shuffle RF dates without replacement, 20,000 draws, 95 % envelope | same ("scrambled control": the RF dates are shuffled among the events) |

WA ETS times: same GPS recipe on a WA station (SC04 / P436 Olympic; SEAT for Puget), cross-checked with the
repo NW tremor list (`episodes_b.csv` region NW).

## 2. Data
- Event list: `data/rf_teleseis/events_M58_1993_2026.csv`.
- Waveforms: `data/rf_teleseis/raw/<NET>.<STA>/<event_id>.mseed` + `raw/<NET>.<STA>.xml` (response). IRIS FDSN
  bulk requests, ≤ 3 in parallel, cached and resumable (per-chunk DONE markers), under nohup; restart = block 14 of
  `scripts/resume_after_reboot.sh`.
- Station set. Core (paper): CN PGC, SNB, VGZ, LZB, PFB.
  Extension (long broadband records over the ETS zone; PB boreholes excluded as unsuitable):
  SVI: CN NLLB (2003–), CN YOUB (2003–14), CN MGB (2005–17), CN GOBB (2017–).
  WA: UW GNW (1995–), UW SP2 (1997–), UW LON (1992–), UW OFR (2001–13, Olympics), UW WISH (2004–), US NLWA (2006–, Olympics).

## 3. How the reproduction is scored (fixed now, before results)
Main configuration: 1993–2019 (their period), 0.05–0.5 Hz, 15-day bins.
For each of PGC, SNB, VGZ we compute ΔCC_SV = mean CC(post) − mean CC(pre) and ΔdVs_top, ΔdVs_bot.
- **R1 sign**: ΔCC_SV < 0 and ΔdVs_top > 0.
- **R2 size**: |ΔCC_SV| between 0.025 and 0.20 (they report 0.05–0.10; we accept half to double);
  ΔdVs_top between 0.025 and 0.2 km/s (they report 0.05–0.1).
- **R3 significance**: ΔCC_SV beyond the 95 % scrambled-control level (one-sided, 20,000 date shuffles) or t-test p < 0.05.
- **R4 shape**: |Δ| at 15 days larger than at 150 days (decay with bin width) for ΔCC_SV.
- **R5 SH**: ΔCC_SH not significant (their result).
A station "reproduces" if R1 + R2 + R3 hold for ΔCC_SV and R1 + R2 for ΔdVs_top. **Reproduced** = at least 2 of
the 3 stations, PGC among them. **Partly reproduced** = signs right at ≥ 2 of 3 stations but size or significance short.
**Not reproduced** = otherwise. LZB and PFB are expected to show nothing (no score). ΔdVs_bot < 0 is reported, not scored
(they themselves had low confidence).
Note on ambiguity: the paper says "mean CC drops ~5–10 %"; we read this as an absolute CC difference of 0.05–0.10 and
also report the relative change (ΔCC / mean CC). Both are shown.

If not reproduced, per the project rule ("reproduce the known answer first"), we debug OUR chain against the paper
(RF section fig S1-like, harmonic components fig S2-like, CC series fig S3-like, ETS dates table S2-like) before
saying anything about the signal. Things to vary in that order: ETS midpoints (tremor- vs GPS-based), the SNR/outlier
rules, the depth conversion, the deconvolution noise window, the station's channel epochs.

## 4. Extensions (after scoring)
E1 2019-10 … 2026 alone (out-of-sample; harmonic model from all data, scrambled control from the same subset).
E2 full 1993–2026. E3 WA stations with WA ETS midpoints. E4 recovery curve: post-ETS ΔCC and ΔdVs in sliding
windows (7-day steps, 7- and 15-day widths) from −60 to +120 days. E5 0.7 Hz variant.

## 5. Output for the LFE search
`target.json`: per station — window of the change relative to the ETS midpoint (and to the local slip-front arrival,
from tremor: first tremor within 50 km of the station), sign, and size as dVs and as % of LVL Vs (LVL Vs taken
from Audet et al. 2009 / Bostock: ~3.0–3.3 km/s; stated), with the scrambled-control percentile.

## 6. Rules
Plots: shaded bootstrap bands like the paper (no error bars on dv/v elsewhere). "Scrambled control", never "null".
Times in PT. CPU ≤ 8, RAM < 40 GB. Kill by PID only; delete nothing. Do not touch s11_depth or the eq-fetch folders.

## Amendment 1 (2026-10-13 ~10:10 am PT, before any result)
- Instrument response is not deconvolved; sensitivity + orientation only (the RF is a ratio of components of one sensor).
- CN waveforms come from NRCAN/CNDC FDSN (EarthScope's CN copy lacks most pre-2010 data).
- E1 (2019-06..2026): the harmonic model is built from that period's RFs alone, so the test is fully independent of
  the reproduction period.
- ETS midpoints entering a run are those inside the run's period (1994-07..2018-07 = 21 episodes for the reproduction).
- WA (E3) uses 2007-2026 with the NW tremor list refined on SEAT GPS; pre-2007 WA timing is not built yet.
- dVs uses the RF value averaged over ±0.5 km at the LVL top/bottom depths picked once per station (not a per-RF
  extremum, which is biased by noise).

## Amendment 2 (Wed 2026-10-14 ~2:45 pm PT) — same-month scrambled control + seasonal regression. Written BEFORE computing.
Why: the scrambled control above shuffles CC values among RF dates and ignores the calendar. ETS midpoints cluster by
season (SVI Aug–Sep, WA May and Aug) and CC_SV swings by month (RESULTS §8, §10), so a seasonal CC pattern could pass as
an ETS effect. Script `scripts/rf07_samemonth.py` (new; reads only the existing `tables/<tag>/<STA>_events.csv` and
`tables/ets_mid_*.csv`, no new RFs). Quantity: `cc_sv_clean` (the same column as §1–9).
**Runs**: core PGC, SNB, VGZ in `repro` (1993–2019) and `full` (1993–2026), PGC also `oos2019`; extension NLLB, YOUB,
GOBB, MGB (`ext_svi`, ALBH midpoints) and GNW, SP2, LON, OFR, WISH, NLWA (`ext_wa`, SEAT midpoints, 2007–2026).
**Same-month scrambled control** (20,000 draws for the bins, 5,000 for the recovery curve): each real midpoint inside the
station's RF span is moved to a random other year (any year with RFs) at the same calendar date ± U{−15…+15} days; a
candidate is rejected if it lies < 60 d from ANY real midpoint of the list, outside the station's first–last RF date, or
< 30 d from a date already drawn in the same draw (up to 200 tries, else that episode is dropped from the draw). RF values
stay on their own dates; only the ETS dates move. Statistic recomputed exactly as rf04/rf05 (nearest-midpoint
assignment; pre [−W, 0), post [0, W); recovery = 15-d window mean minus mean of all RFs, centres −60…+120 every 3 d).
Reported: obs Δ, the same-month median Δ (= what season + sampling alone give), robust σ, one-sided pct (CC down),
and the old shuffle pct beside it. Obs Δ must reproduce the stored rf04 values (check).
**Seasonal regression**: CC_SV = a + b₁cos ωt + b₂sin ωt + c₁cos 2ωt + c₂sin 2ωt (ω = 2π/365.25 d), least squares,
fitted on the RFs > 60 d from any real midpoint ("quiet" RFs, so the ETS signal cannot enter the fit); also fitted on
all RFs (reported). Report annual and semiannual amplitudes, F-test p of the 4 seasonal terms vs a constant (quiet RFs),
the **season-predicted Δ** = Δ of the fitted model evaluated at the RF dates in the real pre/post bins, and the
**season-removed Δ** (CC minus the model) scored with both controls.
**Pass / reading rules (fixed now):**
- S1 "PGC still passes": PGC ΔCC_SV (15 d) has same-month one-sided pct ≤ 0.05 in BOTH `repro` and `full`.
- S2 LON: called **seasonal** if its 15-d same-month pct > 0.05 OR the season-predicted Δ is ≥ 50 % of the observed Δ.
  Called **surviving** if its same-month pct ≤ 0.05 AND its season-removed Δ has same-month pct ≤ 0.05.
  The 150-d persistence survives only if the 150-d Δ has same-month pct ≤ 0.05.
- S3 every other station: keeps its §1/§7/§8 reading if the same-month pct ≤ 0.05 (in the expected direction; GNW's
  opposite-sign result is scored with the two-sided pct). Recovery curves: list the centres outside the same-month 95 % band.
- 14 station-runs × several bins are tested; ~1 in 20 cells passes by chance; this is stated with the results.
- Also in rf07 (pre-registered in `../acf_pgc/PLAN.md` §10, Amendment C2, same time): the pre-ETS window
  H = mean CC_SV [−21, 0) − mean [−60, −25) d with the same-month control, scored on GOBB, MGB, GNW, SP2 (stations whose
  pre-window was never computed); OFR, WISH, NLWA reported; PGC as reference.
