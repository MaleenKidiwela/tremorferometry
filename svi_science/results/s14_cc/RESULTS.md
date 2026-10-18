# S14 RESULTS — does the stretching CC show things that dv/v does not?

Written Wed 2026-10-14, ~2:50–3:05 pm PT. Rules: `PLAN.md` (written 1:59–2:20 pm PT, before any CC-vs-time number);
decision log `NOTES.md` (Amendment 1 after the pilot). Scripts `svi_science/scripts/s14_{common,raw,store,tests,figures,posthoc,posthoc_fig2}.py`.
"CC" = the CC of our existing LFE dv/v stretching measurement (the colour on the dv/v plots). Contrasts are in Fisher z
(Δz); at CC ≈ 0.5, ΔCC ≈ 0.75 Δz; at CC ≈ 0.3, ΔCC ≈ 0.9 Δz. "Scrambled control" = labels or dates shuffled within each
family-station series, keeping the season. "c2" = the detector's noise windows run through the same path. dv/v in %.
Everything is provisional.

## Raw first: dv/v and its CC, same series, side by side
![raw pilot](figures/figA2_pilot_episode_dvv_cc_raw.png)
`figures/figA2_pilot_episode_dvv_cc_raw.png`: families 10, 87, 17 at PGC, B009, B010, B011 (co-located), per episode,
2–5 s, leave-episode-out reference. Also `figA1_fam10_PGC_dvv_cc_raw.png` (family 10 at PGC, daily raw + SVD: the series
on our usual plot) and `figA3_pooled_dvv_cc_raw.png` (median over all series per episode, both regions).

What the raw plots show:
- **Yes, CC moves a lot while dv/v does not.** Episode CC runs 0.1–0.95 within one series; family 10's daily CC spans
  0.2–0.9 inside a single episode. dv/v stays within its usual ±0.2–0.4 % scatter.
- **The CC drops are shared by the four co-located instruments** (e.g. family 87 in 2007, 2013-04, 2018, 2024-08;
  family 10 in 2022-08). They are the same events stacked at four instruments.
- **Pooled, the CC curve follows the number of stacked events** (blue dotted line in A3): episodes with fewer repeats
  have lower CC. The c2 noise-window CC (red) stays flat.
- Exceptions visible by eye: **B009 2007–2009** (family 87: CC 0.08–0.1 at B009 while PGC/B010/B011 are 0.3–0.8).

**Raw "CC low while dv/v flat" list** (series CC < −3 robust σ and |dv/v anomaly| < its bootstrap σ;
`tables/rawlist_*.csv`):

| | cells | CC low | CC low & dv/v flat | episodes where ≥ 3 families do it at one station |
|---|---|---|---|---|
| SVI | 6,273 | 185 | 68 (1.1 %) | 5 (each 3 of 12–56 families measured) |
| WA-NW | 13,736 | 212 | 62 (0.5 %) | 2 |

In those cells the stacks are smaller than usual for their series (median ln n −1.0 SVI, −0.4 WA, i.e. 1.5–3× fewer
events) and noisier (ln noise ratio +0.3). One SVI case (B010, 2009-03 burst) is a known station-timing glitch.

## Why: stack size and noise set most of the CC
![confound](figures/figB1_confound_SVI.png)
- Within a series, Fisher z of the episode CC rises with ln n and falls with the noise ratio (pre-P RMS / coda RMS).
  REG model (series fixed effects + ln n + ln noise + ln n_days), 2–5 s: **explains 47 % of the CC variance in SVI
  (pilot 63 %) and 24 % in WA**. Slopes: +0.13 / +0.11 per ln n; −0.22 / −0.20 per ln noise (SVI / WA).
- The store rebuild reproduces every table CC exactly (median |Δcc| 0.00000; SVI 8,937, WA 18,428 cells).

## Pre-registered comparisons (CC after stack size and noise; PLAN §5–6)
Full run: SVI 76 families × 7 stations (8,893 cells), WA-NW 78 families × 75 stations (18,360 cells); day level from
the same stores (leave-day-out reference). Pilot (families 10/87/17 at the PGC site, 434 cells) gives the same reading
and is in `tables/tests_SVI_pilot.csv`. Contrast figure: `figures/figB2_contrasts_SVI_pilot_SVI_WA_NW.png`.

| comparison (2–5 s unless noted) | SVI Δz ± σ, scrambled p | WA-NW Δz ± σ, scrambled p | pass rules |
|---|---|---|---|
| C1 ETS major − between-ETS burst | +0.003 ± 0.005, p 0.54 | −0.012 ± 0.004, p 0.013 (3–6 s −0.013, p 0.001) | **no**. WA 3–6 s fails: glitch-dropped p 0.66; c2 noise windows show the same −0.012 (activity noise); B013/B943 disagree |
| C2 RF-like: bursts 1–45 d after a major − bursts 1–60 d before | **+0.053 ± 0.020**, p 0.013 (higher after, opposite to the RF) | −0.029 ± 0.014, p 0.10 (RF direction) | **no**. SVI: flagged cells dropped → −0.016, p 0.83 |
| C2b post-ETS bursts − all other bursts | +0.049 ± 0.014, p 0.005 | −0.017 ± 0.010, p 0.11 | no (p > 0.003; signs differ by region) |
| C3 within episode: days ≥ 3 after the local front − days 0–1 | −0.014 vs scrambled −0.016, p 0.81 | +0.006 vs +0.006, p 0.99 | **no front-locked change** |
| C3 days before the front − days 0–1 | +0.001, p 0.23 (direct S −0.042, p 0.002) | −0.014, p 0.49 | no; the SVI direct-S value is also in the c2 noise windows (−0.040) → activity noise |
| C5 seasonal (annual harmonic) | 0.001 z, p 0.84 | **0.0055 z, p 0.004**, peak ~8 Jan (2–3.5 s and direct S also p 0.004) | no (rule p < 0.003); WA = a lead of ~0.005 in CC |
| C6 long-term trend | +0.002 z/yr (41 % of series negative) | −0.001 z/yr (54 %) | none |
| C4c successive episodes, CC vs time gap | +0.012 z per ln-day | +0.002 | no progressive decorrelation over 2003–2026 |

Same reading in every window (direct S, 2–3.5, 2–5, 3.5–6, 3–6 s; Z and H). The 2–3.5 s window (where an LVL-bottom
reflection would land) shows nothing that 2–5 s does not.

**Expected size vs what we can see** (PLAN §2: ΔCC 0.005–0.04, i.e. Δz ≈ 0.007–0.05 at CC 0.5):
- C1 (σ 0.004–0.005 z) would see the upper two-thirds of that range; nothing is there.
- C2 is the RF-like test. Its σ is 0.014–0.020 z, so the 2σ bound is ≈ 0.03–0.04 z (≈ 0.025–0.03 in CC). That rules out
  only the top of the expected range (ΔCC ≈ −0.04). SVI even reads the other way (+0.05, which goes away when the
  flagged cells are dropped).
- C3 (scrambled σ 0.010–0.015 z) rules out a front-locked drop larger than ≈ 0.02–0.03 z.

**Scope of the no-change statements:**
- Medium: the volume the coda window samples, 2–5 s past S. In the S11 kernels this is mostly the path and ~10 km
  around the source, with only ~0.1–0.4 s of extra time near the source.
- Time scales: one episode against the others; days since the local slip front; 1–45 days after an ETS against 1–60 days before.
- Comparison: always within one family-station series.
- This is not the RF's 15-day post-ETS window at a single conversion point. LFE repeats exist only during episodes and
  in bursts, so the post-ETS weeks are sampled only by 240 SVI / 795 WA burst cells.

**C7, CC-only cells after stack size and noise.** Low tail (residual < −3σ): 192 SVI / 306 WA; mirror high tail
(> +3σ): 136 / 279. CC-low with dv/v flat: 79 / 109; CC-high with dv/v flat: 74 / 97. So the low tail is about as heavy
as the high tail. Station-episodes where ≥ 3 families are low and dv/v is flat: 5 SVI / 3 WA. In none of them is the
co-located instrument also low. `tables/cconly_*.csv`.

## What is left: a station-episode CC term that dv/v does not have (C4 + POST HOC)
- **C4 (pre-registered):** after stack size and noise, CC anomalies are shared by the families at one station in one
  episode. M (median correlation of a family with the mean of the others) = **0.28 SVI / 0.24 WA**, against a scrambled
  99th percentile of 0.13 / 0.10 (p ≤ 0.003, the smallest possible with 300 draws). Same in every window (direct S 0.35 / 0.33).
- Co-located instruments share part of it:
  - PGC site: r 0.21–0.41 between PGC, B009, B010 and B011 for the same family-episode (`figB4_coloc_SVI.png`);
  - B013–B943: r 0.53;
  - LRIV group: r 0.15–0.31.
- **POST HOC (after reading C4; `tables/posthoc_*.json`, `figures/figC1_*`, `figC2_*`):**
  - Size: the station-episode term has SD ≈ **0.07 in z (≈ 0.05 in CC)**. It is 27 % (SVI) / 17 % (WA) of the CC
    variance left after stack size and noise.
  - dv/v has no such term: station-episode SD **0.00 % (SVI) / 0.07 % (WA)**. The CC and dv/v station-episode terms
    correlate at only r 0.18 / 0.05. **This is the "CC shows things dv/v doesn't" in numbers.**
  - Not regional:
    - the part common to all stations explains 8 % (SVI) / 2 % (WA) of the variance;
    - it does not follow the pooled dv/v (r 0.07 / −0.03);
    - site-to-site correlation of the station-episode means has a median of 0.12.
  - Not ETS-linked (C1–C3) and not a whole-trace time shift: allowing a ±50 ms shift (cc2p) keeps the term (0.068 → 0.066 z).
    In WA it follows the c2 noise-window CC somewhat (r 0.30). In SVI it does not (r −0.08).
  - The biggest cases are single instruments, and the co-located instrument decides it:

    | station, episodes | CC (Δz, mean over families) | dv/v there | co-located instrument | reading |
    |---|---|---|---|---|
    | **B013, 2012-08 … 2013-05** (4 episodes, 6–35 families) | −0.34 to −0.46 | flat (−0.04 to +0.22 %) | B943 normal (−0.12 to +0.17) | B013 instrument |
    | **B009, 2007-01 … 2009-03** (5 episodes) | −0.34 to −0.52 | also low in bursts (−0.4 to −0.6 %) | PGC/B010/B011 normal | B009 instrument |
    | B013 + B943, 2026-09 burst (4 families) | −0.27 / −0.37 | −0.32 / −0.26 % | both low | small stack, recent; open |

## In plain words
- **What the data show:** the stretching CC does carry things that dv/v does not. Most of it is how many events went into
  the stack and how noisy it was: about half the CC variance in SVI and a quarter in WA. After removing that, a
  station-by-episode part remains (≈ 0.05 in CC). It is shared by all families at that station in that episode, and
  dv/v has no matching part.
- **What it is tied to:** single instruments, not the deep region. The clearest cases are B013 in 2012–2013 (CC down,
  dv/v flat, B943 at the same site normal) and B009 in 2007–2009. It does not line up with ETS, the slip front, the weeks
  after ETS, or the season (SVI; WA a 0.005 lead in January).
- **The RF-type post-ETS drop is not seen in LFE coda CC.** The test can only exclude a drop of about 0.03 CC or more
  (top of the expected range). SVI reads slightly higher after ETS (not robust). WA reads slightly lower (−0.03 z, p 0.10).
- **Useful by-product:** a station-episode CC drop with dv/v flat, and the co-located instrument normal, is a cheap flag
  for instrument trouble. B013 2012-08 … 2013-05 is a new one, not in our glitch list. It should be added to the vetoes
  (needs a look at B013 metadata). **Done 10-14, see the next section.**

## Instrument check and vetoes (Wed 2026-10-14, ~2:40–2:58 pm PT; `tables/instrument_meta/`, `tables/veto_impact_*`)
- **Metadata.** StationXML shows one channel epoch for each station, with no documented change in sensor, gain,
  orientation or sample rate. B013 has one epoch since 2007-01-06 and B009 since 2005-11-02. The faults show up only in
  EarthScope's daily QC metrics (MUSTANG sample RMS per component).
- **B013 → veto 2012-07-13 … 2013-07-24.**
  - From 2012-07-13 both horizontals are 4–60× noisier than co-sited B943 (EH2 up to ~2,800 counts vs a normal 25–50).
    The vertical is normal.
  - The fault ends with a data outage on 2013-07-10 … 07-23, most likely a service visit.
  - This matches the CC drop (2012-08 … B2013-05-24) and CC is normal again at B2013-08-15.
  - After the outage, B013's EH1 runs at about half its old level relative to B943, and stays there. This is kept,
    because the CC is normal. It is the likely undocumented cause of the ballistic B013 horizontal-S step (round 3).
- **B009 → veto 2007-01-08 … 2009-10-01.**
  - On 2007-01-08, B009's EH1 output dropped ~35 % against its own EHZ and EH2 overnight. There were no gaps, no timing
    flags and no metadata change.
  - The winter EHZ/EH1 ratio was 0.65–0.85 before and from 2009-11 on, and 1.2–1.5 in between. It recovered over
    May–October 2009, so the end date is uncertain within that span.
  - This matches the CC drop (2007-01 … B2009-03-09) and CC is normal again at B2009-10-10.
- **How it is applied** (flag, nothing deleted):
  - the windows are listed in `svi_science/instrument_vetoes.csv`;
  - S11 `load_cells` drops those cells (default on; `S11_APPLY_VETO=0` gives the old behaviour);
  - the S14 `glitch` flag now includes them;
  - the tables on disk are all pre-veto.
- **Impact on S11** (linear refits of the existing cells): the vetoed cells carry 0.1 % (SVI) and 0.2 % (WA) of the weight.

  | | before | after veto |
  |---|---|---|
  | SVI Q3 composite (15 majors) | −0.000 ± 0.040 % | −0.002 ± 0.040 % |
  | SVI Q1 / Q2 | 0.175 ± 0.052 / +0.125 ± 0.111 % | 0.177 ± 0.052 / +0.109 ± 0.111 % |
  | WA Q3 composite (32 majors) | +0.012 ± 0.021 % | +0.008 ± 0.020 % |
  | WA Q1 / Q2 | 0.164 ± 0.071 / −0.003 ± 0.089 % | 0.162 ± 0.072 / −0.005 ± 0.089 % |
  | WA with the B013 series also split at 2013-07-24 | | Q2 +0.052 ± 0.093 % (0.6σ shift), Q1 0.168 |

  No reading changes. The gates and injections were not rerun.

## Files
- Tables:
  - `tests_{SVI,WA_NW,SVI_pilot}.csv` (C1–C3, all forms and rule columns);
  - `tests_extra_*.json` (C4–C7, REG fit);
  - `c3_profile_*.csv`;
  - `store_{episode,day}_*.parquet` (rebuilt episode CC, direct-S CC, successive-episode CC, day CC with
    leave-day-out reference, day noise);
  - `store_check_*.txt`;
  - `rawlist_*`, `cconly_*`, `posthoc_*`.
- Figures: `figA1`–`figA3` (raw), `figB1` (confound), `figB2` (contrasts vs scrambled band), `figB3` (day profile vs front),
  `figB4` (co-located), `figC1`/`figC2` (post hoc episode-common / station-episode).
