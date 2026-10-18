# Early-2006 bursts and the 2006-03-28/29 excess: what fires, who carries it, verdicts

Written Sun 2026-10-11 (PT). Follow-up to `../FIRST_LOOK.md` ("Things to look at"). **PROVISIONAL**, like the first look:
P is not applied and the era threshold (9.0) is a placeholder. No Bostock events were read. Only family positions
(`data/catalogue/families_v1_3.csv`, SVI:k = detector family k) were used. 15 of the 130 detector families are not in
v1.3, so their positions come from `e5r_families_2012.csv`; this matters only for 2006-06-18.
Counts are candidates with net_mad ≥ 9.0 (frozen B2 T_core, and era at the same provisional value). "Scrambled control" =
the same detector on the same day, with each station site shifted to another time of that day. 2006 is before DST,
so **PST = UTC − 8**.

## In plain words

- **The Jan / Apr / May / Jun 2006 "noise bursts" are not noise. They are small tremor bursts between ETS episodes.**
  Each one fires a tight group of 2–6 neighbouring families. The envelope tremor detector finds hours of coherent tremor
  on the same days and in the same area. On Jan 16 the record sections show LFE pulses arriving at every station.
- **Why the frozen control tracked the data:** the frozen B2 detector has only four sites: B004, NLLB, VGZ and the "PGC site".
  The PGC site is four co-located stations (PGC, B009, B010, B011) that the scrambled control shifts *together*.
  That one site carries almost the whole network sum: dom_top 0.87–0.91, z ≈ 6 at all four stations, z < 1 at the
  other sites. So its tremor pulses also pass in the control. This is the W060 lesson in another form. Here the
  dominant site is recording real signal, not noise, but the B2 frozen control is still not a fair noise measure on
  tremor days. The era detector adds independent POLARIS sites (KLNB, SILB, SSIB, z 2.4–3.8 in data vs 0.3–0.7 in
  control) and separates data from control 5–160×.
- **2006-03-28/29 (and 27, 30): undecided, leaning towards weak real LFE activity.** It is much weaker and more scattered
  than the other bursts. It is not a one-station artefact, and the SNB gap does not explain it. The case for and against
  is below.

## Verdict per burst

| burst (UTC) | families that fire (data − control) | where (approx.) | envelope flag N (CN4 / all) | era data / ctrl | verdict |
|---|---|---|---|---|---|
| 2006-01-15/16 | 114, 116, 112, 119, 110 (100 % of excess in top 6) | 48.76–48.86 N, 123.58–123.72 W, 42–46 km: NW and down-dip edge of the family area (along-strike +16…+19 km, down-dip +24…+37 km) | 216 / 249, 173 / 196; located 48.80 N 123.70 W | 514 / 90, 260 / 71 | **likely real tremor burst** |
| 2006-04-24 | 87, 88 (98 %) | 48.56 N, 123.30–123.38 W, 41–45 km (down-dip, central) | 207 / 287; 48.55 N 123.35 W | 249 / 14 | **likely real** |
| 2006-05-01 | 17, 16, 30, 8 (99 %) | 48.21–48.35 N, 123.10–123.20 W, 40–42 km (SE end, Juan de Fuca Strait) | 118 / 284; 48.30 N 123.15 W | 319 / 2 | **likely real** |
| 2006-05-03 | 10, 15, 20, 21, 24 (99 %) | 48.22–48.33 N, 123.41–123.46 W, 34 km (S, up-dip of 05-01) | 72 / 181; 48.23–48.30 N 123.35–123.45 W | 350 / 10 | **likely real** |
| 2006-05-11 | 116, 112, 119, 110 | same NW cluster as January | 93 / 98; 48.70–48.80 N 123.55–123.75 W | 89 / 40 | **likely real (small)** |
| 2006-05-25/26 | 112, 114, 116, 110, 109, 119 (94 %) | same NW cluster as January | 32 / 22, 247 / 243; 48.80 N 123.70 W | 156 / 33, 693 / 97 | **likely real** |
| 2006-06-18 | 125, 126, 128, 120, 122 (4 of 5 not in v1.3) | 48.87–48.91 N, 123.88–123.93 W, 40–41 km (furthest NW) | 277 / 201; 48.90 N 123.95 W | 274 / 18 (frozen only 17 / 7) | **likely real** |
| **2006-03-28/29** (and 03-27, 03-30) | 85 / 77 families (frozen); top 6 = only 35–40 % of excess: 75, 77, 72, 73, 74, 81, 89, 90 | frozen top families 48.53–48.58 N, 123.84–124.20 W, 33–35 km; the rest scattered along about 50 km of strike | **CN4 0 / 8** (only 3 CN4 stations); **all 33 / 59**, located 48.05–48.40 N, 123.3–122.8 W (SE edge) | 21 / 1, 23 / 1 | **undecided, leaning weak real** |

Timing within the day:
- Jan 15: frozen detections at 07–10 and 14–24 UTC. Envelope tremor windows nearly all day.
- Jan 16: 00–07 UTC (Jan 15 16:00–23:00 PST).
- May 26: 00–14 UTC.
- Jun 18: 04–23 UTC.
- 03-28: **03:00–13:45 UTC (Mar 27 19:00 – Mar 28 05:45 PST)**.
- 03-29: **07:00–15:00 UTC (23:00–07:00 PST)**.
- 03-30: 01–11 UTC.

## 2006-03-28/29 in detail

**For real (weak LFE / tremor activity):**
- Frozen data is 3× its control (731/243, 732/278). The era detector, with 15 independent stations, gives 21/1 and 23/1:
  small numbers, but a clean separation.
- The excess is not carried by the PGC site alone. NLLB, a separate site 65 km north, has median z 3.0 in data vs 1.25 in
  control (on 03-29: VGZ 2.0 vs 1.0). The scrambled control breaks this alignment, so it is cross-site coherence.
- In 20-min close-ups (`fig_zoom_2006-03-28_0900.png`, `fig_zoom_2006-03-29_0930.png`), 10–30-s emergent pulses appear
  at the same moment on 8–10 stations: the PGC site, TWKB, MGCB, LZB, SILB, VGZ, KLNB. Examples are 09:12 and 09:15 UTC
  on 03-28, and the 09:34 UTC burst on 03-29. The quiet reference day 03-20 shows no such common pulses
  (`fig_zoom_2006-03-20_0818.png`). Jan 16 shows the same pattern, but stronger (`fig_zoom_2006-01-16_0130.png`).
- The night-time (04–13 UTC) 2–8 Hz noise floor on 03-28/29 is 1.3–1.4× the reference nights at PGC and B010, 1.7× at
  SILB and 1.4× at SSIB (`station_rms_hourly_abs.csv`). That is what weak tremor does, but storms do it too:
  03-22/23 is similar, with no detections.
- The all-station envelope variant counts 03-28 and 03-29 as active (33 and 59 windows ≥ k 16). The windows are
  emergent (impulsiveness 1.7–1.9, vs ~2 on the clear tremor days and > 8 for earthquakes) and fall in the same night
  hours as the detections.
- There is no earthquake in the repo catalogues that explains it. The ANSS hypocentre files
  (`data/slab_models/mccrory2012/ANSS_hypos_*`) have only an M1.3 at 49.15 N 122.69 W on 03-27 15:46 UTC. The one local
  event visible in the record section (all stations, ~13:43 UTC on 03-28) produces only a few detections.

**Against / unresolved:**
- The main envelope flag (CN4) is not active: 0 and 8 windows. SNB has no data from 03-26 to 03-29, so CN4 ran on
  3 stations.
- No compact source. 85 families fire, spread over about 50 km along strike (10–90 %: −21…+31 km) and about 23 km
  down-dip, with the top 6 at only 35–40 % of the excess. Every other burst puts 94–100 % in 2–6 neighbouring families.
  The envelope locations (SE edge, Juan de Fuca Strait) do not match the frozen top families (48.53 N 123.9 W). So the
  **area is not determined**. "Southern Vancouver Island / Juan de Fuca Strait" is as close as the evidence allows.
- The signal is invisible in 10-s envelopes (`fig_recsec_2006-03-28.png`). Detections happen only while the PGC site is
  quiet: they stop when PGC/B009/B010/B011 noise jumps at about 15:30 UTC (07:30 PST, working hours). So the PGC site
  still sets when the detector can see it (top station in 40 % of detections).
- No high-frequency depletion. hfr_rel is 1.2–1.7 for both data and control. On the clear tremor days it is 0.1–0.4
  and on the 2007-01 ETS reference 0.74. That is expected if the tremor is too weak to dominate the 2–8 Hz energy,
  but it is not positive evidence.

**What would settle it:** apply the confidence model P to these candidates; stack the 03-28/29 detections of families
72–77 and compare the stack with their templates (moveout across the era stations); run the CN4 envelope detector on the
3 stations with SNB replaced by LZB, as a check only.

## Station metadata near these dates (from the run json files and `bw_gate_common.py`)

- **SNB:** no data 2006-03-26 → 03-29 (back 03-30), exactly over the March burst. 03-26 has no burst (data 4 / control 4),
  so losing SNB does not create detections by itself. Sensor change 2006-04-20 (CMG40T → CMG3ESP): its 2–8 Hz level
  afterwards is about 0.3× the March value.
- **B004:** absent 2006-02-04. Very bursty in January (night 2–8 Hz floor 23–516× its March level), but its z is ≈ 0.5–0.9
  and it never carries detections.
- **VGZ:** no change (next epoch 2007-02-28).
- **PB Olympic boreholes enter the runs:** B005 04-12, B001 04-13, B006/B007 04-21, B003 04-27. All after the
  March burst; they contribute z < 1 on the May days.
- **POLARIS / C8 ends:** MGCB after 05-13, TWKB after 05-14, SSIB 06-16, TSJB 06-19, TWBB 06-20, KLNB 06-22; SILB 07-04;
  LZB capped 07-01. The era detector has no output after 06-30.

## Known events

- Repo catalogues: the PNSN tremor catalogue starts 2009, and no Kao/GSC tremor catalogue is in the repo. The ANSS
  hypocentre files have an **M3.3 intraslab event at 48.57 N 123.53 W, 44 km, 2006-01-15 12:29 UTC**, under the January
  burst area. The envelope tremor ran all day on Jan 15 (from 00 UTC), before this event, so the burst is not its aftershock
  sequence. There is an M2.9 at 48.59 N 124.32 W, 40 km, on 2006-05-11 00:51 UTC, the day of the small May 11 burst.
- A quick web search found no published report of tremor on these exact 2006 dates (PNSN's old ETS pages and the PGC
  TAMS system are mentioned, not checked further).

## Consequence for the first look

The "2006 noise months" line in `FIRST_LOOK.md` should be read as **inter-ETS tremor bursts that the frozen B2 detector
sees through one site**. Their control rate (0.19–0.27 per family-day) is a property of the four-site B2 geometry during
tremor. It is not station noise. This bears on the P4 control-rate target for B2: tremor days inflate the frozen control.

## Files

- Scripts:
  - `svi_science/scripts/bw_burst2006_feat.py`: summaries, station shares, hourly counts, family counts, envelope windows,
    `fig_hourly.png`, `fig_families.png`.
  - `bw_burst2006_recsec.py`: `fig_recsec_<date>.png`, 2–8 Hz envelope sections plus 90-s N-component zooms.
  - `bw_burst2006_zoom.py`: `fig_zoom_*.png`.
  - `bw_burst2006_noise.py`: `station_rms_hourly_abs.csv`.
  - All are read-only on `mf/feat`, the 20 Hz caches and `p4/envflag/win`; single process.
- Tables: `day_summary.csv`, `station_share.csv` (per-station median z, fraction top, fraction z ≥ 3, data vs control),
  `hourly.csv`, `fam_counts.csv`, `env_windows.csv`, `station_env_hourly.csv`, `station_rms_hourly_abs.csv`.
