# Q6 — What stops an ETS front? And is the Strait of Juan de Fuca a boundary?

Written 2026-10-07 ~5:02 pm PT. Rules: `PLAN.md` (pre-registered 4:55 pm PT; changes and one post-hoc check in its §9).
Scripts `svi_science/scripts/q6_{build,tests,lfe,posthoc,figures}.py`. Primary data = PNSN tremor 2010-01 → 2026-05,
the whole margin from N. California to central Vancouver Island on one along-margin axis (0 = mid-Strait, + toward NW).

## Plain answer

- **Fronts tend to stop where the fault slipped recently.** As a tremor front moves along the margin, the chance that
  it stops at the next 5-km stretch is about **10 %** if that stretch had tremor within the last month. It drops to
  about **7 %** if the stretch had been quiet for more than about 3 months. In other words, a front runs farther into
  fault that has been quiet longer. The reverse idea, that fronts stop where the fault has been quiet longest, is wrong
  in the tremor: the sign is the opposite every time.
  The effect is real against the scrambled control (episodes moved to each other's dates), with p < 0.001. It holds
  in every version: since the last big ETS, in Vancouver Island + Washington alone, in 2015–2026, and when recent
  activity is ignored. But it is **modest**. Recent slip is one push among many and not a wall: most fronts stop at
  stretches that were not recently active.
  - By the frozen rule this is **PARTIAL**, not SUPPORTED. The family-resolved SVI LFE check has the wrong raw sign.
    That check has almost no power: only 6 stop events inside the SVI family box, because fronts usually run out of
    the box. Its own scrambled control is even more negative.
- **The Strait of Juan de Fuca is NOT a special place for starting or stopping ETS.**
  - Starts: 8 of 133 big episodes start within ±15 km of mid-Strait. Starting points spread evenly along the margin
    would put about 10 there; starting points that follow where tremor occurs would put about 14 there.
  - Stops: 29 front stops fall in the Strait zone, against 25 from where tremor occurs (p 0.25) and 23 when each run
    is shifted randomly by up to 100 km (p 0.08).
  - Fronts that reach the Strait cross it a bit less often: 38 % against 54 % elsewhere (p 0.05). But the low-crossing
    stretch runs from the Strait to beyond Lake Cowichan (about −15 to +100 km). It is not a narrow line at the Strait.
  - The Strait zone is not a gap in tremor detection: it has more epicentres than the stretches next to it.
- **S3's "5 of 24 starts near the Strait" came from the edges of the LFE boxes.** The Washington LFE family box ends at
  about −20 to −30 km on the margin axis, and the SVI box begins just north of the Strait. So fronts entering either box
  "start" there. Merging the overlapping SVI and Washington episodes into one along-margin record leaves 2 of 18 starts
  in the Strait zone, against about 2 for starts spread evenly and 3 for starts that follow where the families are.

## Numbers (with controls)

### Q6.1 Front arrest (PNSN tremor, 475 big episodes, 9,494 walk steps, 875 stops)

Model: front continues into the next 5-km bin (1) or stops (0) ~ log10(days since that bin last had tremor) + distance
from start + episode size + the bin's usual tremor rate. β = logit change per 1 SD of log time. Referee = episode-order
shuffle (2,000 draws; each episode's footprint moved to another episode's date; keeps each bin's usual rate).

| version | β real | scrambled control mean (5–95 %) | one-sided p | rows / stops |
|---|---|---|---|---|
| **stop hazard, any earlier episode (PRIMARY)** | **+0.135** | −0.039 (−0.100…+0.020) | **0.0005** | 9,494 / 875 |
| stop hazard, since last big ETS | +0.168 | −0.019 | 0.0005 | 9,377 / 870 |
| stop hazard, Vancouver Is. + WA only (−200…+200 km) | +0.215 | −0.014 | 0.001 | 2,782 / 244 |
| stop hazard, ignoring episodes that ended < 7 d before | +0.114 | −0.043 | 0.0005 | 9,487 / 875 |
| same, VI + WA only | +0.123 | −0.014 | 0.020 | 2,782 / 244 |
| stop hazard, censored rows kept | +0.131 | −0.038 | 0.0005 | 9,632 / 888 |
| stop hazard, 2010–2014 / 2015–2026 | +0.091 / +0.173 | −0.030 / −0.037 | 0.043 / 0.0005 | 2,323 / 7,171 |
| enter vs not (±150 km of start), any / big | +0.135 / +0.256 | −0.049 / −0.058 | 0.0005 / 0.0005 | 25,551 |
| enter vs not, VI + WA, any / big | +0.268 / +0.482 | −0.018 / −0.014 | 0.001 / 0.0005 | 7,216 |
| *post-hoc*: only rows with ≥ 14 d / ≥ 60 d since last tremor | +0.166 / +0.199 | −0.029 / +0.002 | 0.001 / 0.006 | 7,151 / 3,143 |

The post-hoc rows rule out a possible artefact of how episodes are built. Two bursts closer than 2 days and 15 km are
merged into one episode, so a front can never be seen entering a bin that another burst used within about 2 days.
Above 14 days, the stop chance still falls from 10.3 % (14–30 d) to 9.7 %, 8.1 % and 7.0 % (> 99 d). The scrambled
control stays flat at 8.3–10.5 %.

Figure: `figures/q6_arrest_vs_tsince.png`. The stop chance by quintile of time since last tremor is 10.6, 10.3, 9.7,
8.6 and 6.8 %. The scrambled control stays flat at about 9.1–9.5 % (5–95 % band 8.2–10.5 %).

**Pass rule (PLAN §3).** PNSN primary p < 0.01 ✔. VI+WA β > 0 ✔. SVI LFE β > 0 ✘ (−0.35) → **PARTIAL**.

**LFE family-resolved checks.**

| catalogue | walk β (p) | enter β (p) | stops | note |
|---|---|---|---|---|
| SVI ours 2014+, any episode | −0.35 (0.61; control mean −0.46) | +0.30 (0.22) | 6 | fronts usually leave the ~100-km box, so there are almost no stops |
| SVI ours, since last big | +0.17 (0.47) | +0.95 (0.08) | 6 | |
| Bostock 2003–2013 | — | — | 0 | all 10 of his ETS cover every bin of his families, so there is no arrest information |
| joint SVI 2003–2026 (big only) | +0.31 (0.27) | +0.19 (0.08) | 6 | |
| Washington-NW E3 (column; lumpy positions) | **+0.85 (0.017)** | **+0.73 (0.015)** | 13 | history = majors only |

The Washington LFE check agrees with the tremor (p ≈ 0.02). The SVI LFE check is too small to say anything.

### Q6.2 Strait of Juan de Fuca (Strait zone = ±15 km of mid-Strait; Victoria +9, Port Angeles −17 km)

| test (Vancouver Is. + WA, −200…+200 km) | in Strait zone | uniform | matched to tremor | other control | verdict |
|---|---|---|---|---|---|
| starts of big episodes (n 133) | 8 | 10.0 (p 0.79) | 13.9 (p 0.97) | own-footprint 13.5 (p 1.0) | **no clustering** |
| starts, ±10 / ±25 km | 2 / 17 | 6.7 / 16.6 | 9.3 / 22.7 | | none |
| stops of big episodes (n 242) | 29 | 18.2 (p 0.010) | 25.3 (p 0.25) | rigid shift 22.8 (p 0.08) | **not beyond controls** |
| all episodes ≥ 30 epicentres: starts (n 455) | 52 | 34.1 (p 0.002) | 47.6 (p 0.27) | | explained by the high tremor density there |
| all ≥ 30: stops (n 866) | 104 | 65.0 (p 0.0001) | 90.6 (p 0.08) | shift 74.7 (p 0.0005) | column only: more small-episode stops than the shift control, not beyond the tremor-matched control |
| whole margin, big starts (n 460) | 8 | 11.2 | 16.9 | | no |

- **Crossing.** 15 of the 40 fronts that reached the Strait crossed it (38 %). Across the other 30-km stretches the
  figure is 210 of 386 (54 %); Fisher p 0.047. The Strait ranks 4th lowest of 75 overlapping stretches. But the dip is
  as deep from +15 to +45 km (38 %) and at +75…+105 km (39 %). Fronts stall along the whole Victoria → Lake Cowichan
  stretch, not just at the Strait.
- **Detectability.** Epicentres per 30-km stretch: Strait 26,641; next stretches 23,991 (SE) and 23,577 (NW). The
  Strait is not a detection hole.
- **LFE check** (margin axis, n 18 after merging overlapping SVI + WA big episodes): 2 starts in the Strait zone,
  against 1.8 for starts spread evenly (p 0.56), 3.0 for starts that follow where families are (p 0.82) and 4.5 for the
  episode's own events (p 1.0). Unmerged, as S3 counted them (n 24): 4 against 2.4 / 4.0 / 2.7 (p 0.22 / 0.58 / 0.05).
  The Washington episodes "near the Strait" start at −35 to −40 km on the margin axis, at the north edge of the
  Washington family box.

## Figures

- `figures/q6_fronts_margin.png` — every big tremor episode 2010–2026 along the whole margin (top) and zoomed on
  Vancouver Island + Washington (bottom). The grey bar is the front run, the black dot is the start, the blue ▼ / orange
  ▲ are the NW / SE stops, and the Strait is shaded.
- `figures/q6_arrest_vs_tsince.png` — stop chance against time since last tremor (real against the scrambled control),
  and β for every version against the control band.
- `figures/q6_strait.png` — starts and stops along −200…+200 km with the uniform and tremor-matched expectations, and
  the crossing rate per 30-km stretch.

## Caveats

- "Front" here means the contiguous along-margin run of 5-km bins with ≥ 3 epicentres, walking out from where the episode
  started (one empty bin is bridged). Episodes are space-time clusters of PNSN tremor: 2 days, 15 km. A big episode is
  one with ≥ 206 epicentres (data-driven split), giving 475 margin-wide.
- "Last active" counts any tremor cluster, including small inter-ETS bursts. The median time since last tremor at a
  front step is about 1 month, so in this tremor record "recently slipped" mostly means slipped in the last weeks to
  months. The "since last big ETS" version gives the same sign.
- The arrest effect is a tendency, not a rule. The stop chance per 5-km step goes from about 10 % to about 7 %. The
  family-resolved SVI check cannot confirm it (6 stops); Washington LFEs agree (p 0.02).
- Every result comes from a single along-margin coordinate; down-dip structure is ignored.
