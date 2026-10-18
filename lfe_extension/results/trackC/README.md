# Track C — production catalog (2014–2026), outputs and how to run it

Everything here is produced by `scripts/c_*.py` on the **frozen measurement chain**
(`results/frozen_operating_point.json`, Merlin F6 + Track A step 1):
20 Hz · 1–8 Hz Butterworth-4 zero-phase · component order Z / N|1 / E|2 resolved, **not** rotated
(borehole EH1/EH2 → N/E) · 12-s templates · site-merge **B** (`{PGC,B009,B010,B011}` and
`{LRIV,B005,B006,B007}`, site weight √Σw²) · SNR weights `e5r_weights_snr_2012.csv` ·
per (family, day) network MAD · **T_ext = 8.25**, **T_core = 8.75** · HF-ratio veto at the PGC site
(E(10–19 Hz)/E(2–8 Hz) over the scored window, cut **3.6384**).

**Verified**: `c_run_month.py` reproduces the Track A frozen run (`a1_mergeB_L12`) **bit for bit** —
2012 doy 228–232 (quiet) 253/253 detections and doy 248–252 (ETS) 3426/3426 detections, with
Δ net_mad = Δ net_cc = Δ amplitude = **0.0** — with `--min-sites 0 --cc-device {gpu|cpu}` (the CPU
offload path is numerically identical, it only changes where the float32 series are held).

> **Template status (Track A closed 2026-09-16).** Template iteration was REJECTED on the held-out
> year, so the from-scratch 2012 set (`e5r_templates_2012.npz`, 12 s) is **final** for every
> station-epoch whose instruments match 2012. Every detection carries `epoch_template`:
> `2012_final` when all contributing stations are on such an epoch, `2012_provisional` when at least
> one is past its post-2012 instrument boundary (NLLB 2017-02-23, VGZ 2017-03-03, SNB 2017-06-08,
> PGC 2017-08-10, LRIV 2019-05-30) and is therefore being matched with a template built on different
> instruments, and `per_epoch:<file>` once Track A's A-2 agent drops in
> `results/a4_templates_<STA>_<YYYY-MM-DD>.npz`. `n_sta_provisional` counts the provisional stations
> behind each detection. **No code change is needed when those files appear** — see
> "Per-epoch templates" below.

---

## Pipeline

| script | what it does |
|---|---|
| `scripts/c_month_ready.py` | is a month runnable? which stations contribute? (reads `results/download_manifest.csv` + `results/station_epochs.csv`) |
| `scripts/c_cache_month.py` | 20 Hz cache for a month → `cache20hz/<NET.STA>/<YYYY>/<DOY>.npz` + `daily_noise_<YYYY-MM>.parquet` |
| `scripts/c_run_month.py` | the network matched filter for a month → `detections_`, `family_days_`, `stacks_`, `run_` |
| `scripts/c_qc_month.py` | monthly QC → `qc_<YYYY-MM>.csv`, appends `qc_summary.csv`, `figures/trackC_monthly_qc.png` |
| `scripts/c_eq_catalog.py` | one-off: the regional earthquake catalog for the EQ flag → `eq_catalog.csv` |
| `scripts/c_watch.py` + `scripts/c_run_all.sh` | the resumable watcher that chains all of the above, month by month |

```bash
nohup bash scripts/c_run_all.sh > /dev/null 2>&1 &   # start / restart the production watcher
tail -f logs/trackC_watch.log                         # what it is doing
cat results/trackC/months_done.csv                    # what is finished
kill $(cat logs/trackC_watch.pid)                     # stop — safe at any point, fully resumable

bash scripts/c_run_all.sh --once                      # one month then stop
bash scripts/c_run_all.sh --months 2015-12,2018-06    # a named list, in order
bash scripts/c_run_all.sh --allow-partial             # dry runs: run months whose download is incomplete
bash scripts/c_run_all.sh --from-ym 2014-01           # skip the 2010-2013 referee windows
touch results/trackC/STOP                             # graceful stop: finish this month, then exit
```

The watcher holds an `flock` on `results/trackC/.watch.lock`; a second one exits immediately.
**Epoch ends are EXCLUSIVE**, matching `scripts/download_plan.py`: a station-day ON an epoch end date
was queued `no_epoch` and never fetched, so treating it as inclusive made five months wait forever on
a station-day that will never arrive (2014-06 CN.LZB 06-17, 2014-09 CN.YOUB 09-03, 2015-01 C8.SOKB
01-18, 2017-02 CN.NLLB 02-21, 2018-12 CN.PFB 12-20 — all runnable since 2026-09-17). On a band-switch
date the two epochs abut, so the station stays alive through the new epoch and nothing is lost.
It loops: oldest **runnable**, unprocessed month → cache → run → QC → append to `months_done.csv`;
if nothing is runnable it sleeps 10 min and looks again. A month is **runnable** when every
station-day of every station alive in that month has a finished manifest row
(`data|nodata|error|present_elsewhere|no_epoch`); a station **contributes** on a day when its row is
`data|present_elsewhere` with `coverage_frac ≥ 0.5`.

The H100 is shared with other agents: `c_run_month.py` waits at the start of each day until
`--gpu-need-gb` (6 GB) is free, chooses per day whether to hold the per-station cc series on the GPU
or in host RAM (`--cc-device auto`, threshold `--gpu-free-gb` 20 GB), halves the family chunk on a
CUDA OOM and falls back to the host if even the allocation fails. None of this changes a number.

---

## Files

### `detections_<YYYY-MM>.parquet` — one row per declustered network detection
Every detection **at or above `record_thr` = 8.0 MAD** is written, so both tiers are a *filter*, never
a re-run. Rows are never deleted: the HF veto, the tiers and the cross-family assignment are columns.

| column | meaning |
|---|---|
| `t`, `time` | detection time: epoch seconds, and the same as datetime64. This is the ORIGIN-axis time used by the network sum (station template start = `t + ttP + (S−P) − 0.5 s + moveout`). |
| `fam` | Bostock family (`bkey`, 0–129) |
| `date`, `doy` | UTC day of the detection |
| `net_mad` | network score in MAD of THAT family-day's own network series |
| `net_cc` | raw weighted network cc at the peak (the quantity the cross-family rule compares) |
| `mad_scale`, `med_scale` | the family-day MAD and median the score was formed from |
| `tier_ext`, `tier_core` | `net_mad ≥ 8.25` / `≥ 8.75` |
| `n_sites`, `n_sta` | distinct SITES and stations contributing that family-day (≥ 2 sites always, by the observability rule) |
| `rel_amp` | SNR-weighted mean of cc·√(E3/E3_ref) — units of each station's 2012 quiet-day window rms |
| `con_pgc`, `wsh_pgc` | the PGC site's contribution to the network sum and its share of the total weight (dominance diagnostic) |
| `hfr` | E(10–19 Hz)/E(2–8 Hz) at CN.PGC over the scored window, from the NATIVE-rate trace; NaN when PGC has no data |
| `hf_pass` | `not (hfr defined and hfr > 3.6384)` — an undefined ratio is never vetoed (same rule as the calibration) |
| `hf_veto_adopted` | `True` — the frozen point adopted the veto, so the product tiers are `tier_* & hf_pass` |
| `tremor_flag`, `n_tremor_day` | a PNSN tremor detection inside the family box within ±2 h / the count in the box that UTC day |
| `eq_flag`, `eq_dt` | a catalogued earthquake within 60 s / the distance in seconds to the nearest one |
| `keep_xfam_ext`, `keep_xfam_core` | **FROZEN cross-family assignment** (Merlin Step 7) evaluated on the T_ext and the T_core set: among families sharing ≥ 50 % of their (2012) stations, a detection is dropped if a station-sharing family has a higher RAW weighted cc within ±2 s |
| `keep_mad_ext`, `keep_mad_core` | the older max-network-MAD assignment, stored for comparison only |
| `epoch` | the network configuration of that day (`STA:band|…`) |
| `quiet_flag` | the family-day was quiet (no PNSN tremor in the box that day and ≥ 3 days from any tremor day) — Track A step 3b found ~22 % of quiet-day T_core detections are REAL inter-ETS repeats (≤ 14 % at T_ext), so this is the label the per-family per-season floor λ is computed from, **not** a false-alarm label |
| `cat_ext`, `cat_core` | **catalog membership**: `tier_* AND hf_pass AND keep_xfam_*`. Vetoed and cross-family-suppressed rows stay in the file (every column preserved) so λ can be re-derived; they are simply not members |
| `epoch_template` | `2012_final` \| `2012_provisional` \| `per_epoch:<file>` — see the note at the top |
| `n_sta_provisional` | contributing stations whose template set is labelled `2012_provisional` (0 under the current policy) |
| `n_sta_post_boundary` | contributing stations being matched ACROSS a post-2012 instrument change (the real diagnostic) |
| `n_sta_excluded` | stations the epoch policy dropped from this family-day (LRIV after 2019-05-30) |
| `hf_cut` | the HF cut actually applied to this row (frozen 3.6384, or the per-epoch cut when `results/a4_hf_veto_cuts.csv` exists) |
| `hf_fs` | the HF station's NATIVE sampling rate that day (40 Hz PGC BH, 100 Hz PGC HH). **The frozen cut 3.6384 was calibrated on PGC BH 40 Hz and is NOT epoch-invariant** — see the caveat below |

**The catalog** = `cat_core` everywhere, plus `cat_ext` inside tremor windows or with its false-alarm
floor λ attached (PLAN.md amendment; Merlin F1).

### CAVEAT — the frozen HF cut is epoch-dependent
The veto cut 3.6384 is the 99th percentile of the build-half own-matched ratio at **CN.PGC BH, 40 Hz**,
whose anti-alias filter rolls off inside the 10–19 Hz band. On **PGC HH, 100 Hz** (from 2017-08-10) the
same band is not rolled off, so the whole ratio distribution shifts up: measured here, the median HF
ratio of the detections goes 0.199 (2015-12, BH) → 0.480 (2018-06, HH) and the fraction above the frozen
cut goes 2.0 % → 12.0 %. Applying the frozen cut unchanged therefore vetoes ~6× more detections after
2017-08-10 than before it. `hf_fs` is stored on every row and `qc_summary.csv` tracks the shift
(`hfr_median`, `hfr_frac_above_cut`, `flag_hfr_shift`), but **the cut itself needs a per-epoch
re-calibration from Track A** (same rule, own-matched 99th percentile, per instrument epoch) before the
post-2017 catalog is used. Until then, treat `hf_pass` after 2017-08-10 as conservative.

### TIERS PER NETWORK EPOCH — `results/trackC/tiers_by_epoch.csv`
T_core and T_ext are FALSE-ALARM thresholds, and the false-alarm rate at a fixed MAD is a property of
the NETWORK, not of the year: a smaller network is noisier, so holding 8.25 / 8.75 fixed across the
record would mean a catalog whose purity silently changes with the station list. The two tiers are
therefore re-derived **per network epoch by the same FA-only rule that defined them**
(`scripts/c_tiers_by_epoch.py`):

* `T_core` = smallest T on a 0.25-MAD grid at or above the 8.0 recording floor with quiet-day
  FA ≤ 0.1 per template-day; `T_ext` = smallest T with quiet-day FA ≤ 1.0 **and** post-SSE FA ≤ 1.5.
* quiet family-day = observable, no PNSN tremor in the family box that day, ≥ 3 days from any tremor
  day (the stored `quiet_flag`); post-SSE family-day = observable, 1–10 days after the last tremor day
  of an episode (an episode ends when ≥ 3 tremor-free days follow) — the production analogue of Track
  A's 2012 doy 270–279 window.
* FA(T) counts detections with `net_mad ≥ T` **and** `hf_pass` on those family-days (the frozen point
  was derived with the veto), divided by the number of those family-days. No cross-family assignment
  (the rate is per template per day) and no control correction — the raw rate, as in Track A.
* Epochs: **E0** 2014-01..2017-02 (15 stations) · **E1** 2017-03..2018-12 (post CN sensor upgrade) ·
  **E2** 2019-01..2021-03 (no PFB, LRIV dropped from 2019-05-30) · **E3** 2021-04.. (no SNB).
* **E0 anchor rule**: if E0's own quiet days reproduce the frozen 8.25 / 8.75 within 0.25 MAD, the
  frozen pair stays. They do (E0 derives 9.00 / 8.00 on 36,865 quiet family-days), so **E0 = 8.75 /
  8.25**, unchanged.
* A row already marked `derived` is FROZEN — the tiers do not drift as months arrive. Pass
  `--rederive` to recompute one deliberately.

Because `net_mad` is stored, a tier change is post-processing, never a re-run:

```bash
python scripts/c_tiers_by_epoch.py            # (re-)derive; writes tiers_by_epoch.csv
python scripts/c_reflag_tiers.py --all --dry-run
python scripts/c_reflag_tiers.py --all        # rewrite tier_*, keep_xfam_*, cat_*, family-day counts
```
`c_reflag_tiers.py` recomputes the cross-family assignment as well, because that rule is evaluated ON
the thresholded set. `net_mad`, the templates and the HF cut are never touched. Every detection carries
`T_core_used`, `T_ext_used` and `tier_epoch`, and the runner reads `tiers_by_epoch.csv` for new months.

The watcher does this automatically: after every month it re-runs `c_tiers_by_epoch.py` (which can only
ADD an epoch, since derived rows are frozen), and the first time an epoch becomes derivable it runs
`c_reflag_tiers.py --all` once, so an epoch never ends up carrying two different thresholds across the
record. `logs/trackC_tiers.log` records both.

**Status 2026-09-17**: only E0 has enough processed quiet days to derive. E1 and E3 each have a single
processed month and it is an SSE month with ZERO quiet family-days; E2 has no processed month at all.
All three sit on the E0 anchor until their quiet months are processed, at which point
`c_tiers_by_epoch.py` will derive them (they are not frozen) and `c_reflag_tiers.py --all` must be run
once to apply them across the record.

### Epoch template POLICY — `results/trackC/templates_by_epoch.csv`
This file, written from Track A-2's step-4 epoch-transfer results, is the authority on what each
station-epoch uses. Columns (`sta`, `epoch_start` required):

| column | meaning |
|---|---|
| `sta`, `epoch_start` | the station and the EFFECTIVE epoch start the row applies to (`YYYY-MM-DD`; the later of the StationXML epoch start and the station's post-2012 instrument boundary) |
| `npz` | a per-epoch template file to use instead of the frozen 2012 set; **empty means "keep the 2012 set for this epoch"** and that decision is explicit, not a default |
| `label` | the `epoch_template` value to stamp (e.g. `2012_final`) |
| `drop` | `1` = **exclude this station from the network for this epoch entirely** |
| `note` | why |

Current contents (A-2 step 4, judged against the B926 no-change negative control, not against the
unusable absolute bars):
* **PGC 2017-08-10, NLLB 2017-02-23, SNB 2017-06-08, VGZ 2017-03-03 → `2012_final`, keep the 2012
  templates.** All four transfer across their 2017 sensor upgrade at least as well as B926 — a station
  with no instrument change at all — transfers across the same 2012→2018 gap (cc ratios 0.90 / 0.87 /
  0.84 / 0.72 vs the control's 0.79; registration 0.94 / 0.88 / 0.89 / 0.71 vs 0.66), the median S-pick
  shift is +0.00 s at every boundary, and there are no real polarity flips.
* **LRIV after 2019-05-30 → DROPPED from the network.** Its stack cc against the 2012 template goes
  NEGATIVE (−0.216) where the matched-N control leg on the same families is +0.542, and its dS MAD is
  0.35 s, 7× every other station. Every family-day records the loss in `excluded_sta` /
  `n_sta_excluded` so the completeness model can absorb it.

Because of this policy `epoch_template` is `2012_final` everywhere; the diagnostic that a station is
being matched across an instrument change is kept separately as `n_sta_post_boundary`.

> **Only one agent may write this file.** On 2026-09-16 the Track A-2 agent overwrote it with a version
> pointing every 2017 boundary at `results/a4_templates_merged_<STA>_<date>.npz` and containing **no
> LRIV row**, under which LRIV would have kept carrying 2012 templates past 2019-05-30. That version is
> preserved verbatim as `templates_by_epoch.a4_proposed.csv`; the per-epoch stacks it names are a
> SENSITIVITY option that the coordinator has not adopted. If the policy is ever changed, re-run every
> affected month — the template set changes the detections themselves, unlike the HF cut.

### Per-epoch templates (drop-in, no code change)
`c_run_month.py` builds ONE template set per (station, effective epoch) and looks for, in order: the
`npz` column of `templates_by_epoch.csv` above, then `results/a4_templates_<STA>_<YYYY-MM-DD>.npz`, then
`results/a4_templates_<STA>_<YYYYMMDD>.npz`. If one is found it is used for that station in that epoch
(its own `tpl` and `sp`, its own placement) and the label becomes the file's basename. A policy row with
an empty `npz` pins the frozen 2012 set for that epoch and suppresses the filename search, so the A-2
files `results/a4_templates_{PGC,NLLB,SNB,VGZ}_2018_HH.npz` are deliberately NOT used (their epoch tags
do not match the date convention either). A station-epoch with no template set at all is a hard error —
the runner refuses to run a station silently blind.

### Per-epoch HF-veto cuts — `results/a4_hf_veto_cuts.csv` (pending from A-2)
This file now exists. The runner reads it and stamps `hf_cut` and `hf_cut_src` on every detection; the
frozen 3.6384 remains the fallback. Two schemas are accepted — A-2's
(`station_epoch` like `PGC_BH_2010-08-01..2017-08-10`, `band`, `fs`, `cut`, `method`, …) and the generic
(`sta`, `fs`, `epoch_start`, `epoch_end`, `cut`) — with three rules:
* rows whose `band` is not the **10–19 Hz** band the stored `hfr` was computed in are **dropped**; a cut
  calibrated in another band would be silently wrong;
* the MOST CONSTRAINED matching row wins;
* among equally constrained rows **the anchor wins**: a row whose `method` says it is the frozen referee
  cut cannot be displaced by a later recalibration of the same epoch, so the operating point the catalog
  was frozen on never moves. Anything the anchor rule does not settle is reported inside `hf_cut_src`.

In force today: **BH epochs 3.6384 (frozen, unchanged), HH epochs (after 2017-08-10) 3.5129.** Months already processed do **not** need the matched filter re-run —
`hfr` and `hf_fs` are stored, so:

```bash
python scripts/c_reflag_hf.py --all --dry-run   # what would change
python scripts/c_reflag_hf.py --all             # re-flag hf_pass / cat_ext / cat_core + regenerate QC
```
The previous values are kept as `hf_pass_prev` / `hf_cut_prev`, a summary goes to
`results/trackC/hf_reflag_log.csv`, and nothing else in the file is touched.

### `family_days_<YYYY-MM>.parquet` — one row per (family, day)
`fam`, `date`, `doy`, `observable`, `n_sites`, `n_sta`, `sites` (which ones), `mad_scale`,
`med_scale`, `n_valid_s` (seconds of valid network series), `n_det_rec` / `n_det_ext` / `n_det_core`,
`n_tremor_day`, `quiet_flag` (no tremor in the box that day AND ≥ 3 days from any tremor day),
`epoch`, `excluded_sta` (stations the epoch policy dropped that would otherwise have had a template —
the completeness input for the LRIV drop). `observable = False` means fewer than 2 sites had a usable template and data — the
family-day is **not scored**, which is what the completeness model needs (Merlin F8).

### `stacks_<YYYY-MM>.npz` — per (family, station, day) stacks of the T_core detections
`fam`, `sta`, `doy`, `n` (events in the stack), `stack` (n_entries, 3, 320) float32 — the linear mean
of rms-normalised 16-s windows starting 2 s before the template start at that station, 20 Hz, same
component order as the cache. Only non-empty combinations are stored. For re-stacking and for the
epoch-identity test (Merlin Step 4).

### ROBUST ENERGY FLOOR in the cc normalisation (bug fix, 2026-09-17)
Found by Track E (`results/trackE/bad_persite_detections.csv`). The matched filter normalises each
correlation by the sliding window's own energy, and the old guards were a denominator clamp at
1e-10 x the all-channel median plus `VALID = (3-C window energy > 1e-6 x its median) AND (every channel
energy > 0)`. Neither catches the real failure mode: gaps are merged with `fill_value=0`, so a dead
channel is **exactly zero** for most of the day (UW.LRIV 2019-02-06: one channel zero on 1,725,189 of
1,727,760 windows). The 3-C SUM is then carried by the two healthy channels and passes both tests,
while the windows that STRADDLE the edge of the dead interval have a tiny-but-POSITIVE energy on that
channel — so `min > 0` admits them, the denominator is not clamped, and
`cc = numerator / sqrt(tiny)` blows up to ~1e11.

The guard is now a floor **relative to the station-day median, applied per channel**, and the cc is
**set to zero** at masked windows rather than merely masked, so a blown-up value can never leak into
the site mean, the network sum, or the declustering max-pool:

```
med3    = median(e3)                     # 3-C window energy: the ROBUST SCALE (positive while ANY
                                         # channel is alive; median(ew) is NOT usable - on UW.LRIV
                                         # 2018-03-02 it is exactly 0 and every relative threshold
                                         # collapsed to 0 on the very station that needed the guard)
ref_ch  = max(med_ch, 1e-2 * med3)       # a dead/near-dead channel gets the fallback reference
bad     = any_c( ew[c] < 1e-3 * ref_ch[c] )  or  e3 < 1e-3 * med3  or  min_c ew[c] <= 0
if med3 <= 0: bad = all                  # every channel dead -> the station-day is unusable
bad     = dilate(bad, +-n)               # window j covers [j, j+n): every window sharing a bad SAMPLE
cc[:, bad] = 0 ;  VALID = ~bad
```
Both the fallback reference and the dilation are load-bearing, and each was found by a failed attempt:
* without the dilation the fix BACKFIRES — zeroing only the core bad windows removes the huge value
  that was suppressing its neighbours in the declustering max-pool, and contaminated neighbours the
  plain energy test does not reach become detections instead (2018-03-02 went from 28 impossible rows
  to 35, with `net_cc` up to 1.7e14);
* without the `med3` scale the guard is a no-op wherever one channel is dead for a large fraction of
  the day, because `median(ew)` is then 0.
`--energy-floor` (default 1e-3) exposes the fraction; `n_windows_masked_by_energy_floor` is recorded
per station in `run_<YYYY-MM>.json`.

**Effect on the frozen 2012 run**: the detection SET and every raw `net_cc` are unchanged bit for bit
(253/253 on doy 228–232, 3426/3426 on doy 248–252, max |Δnet_cc| = **0.0**). `net_mad` moves by at most
4.7e-4 (quiet) / 2.4e-3 (ETS), and only because two 2012 station-days contain a real data gap
(PB.B001 doy 230, CN.NLLB doy 250) whose gap-edge windows (6,050 and 3,293 after dilation) are now
excluded from the network median/MAD estimate — i.e. the change is exactly the windows the old floor
should have excluded. No detection anywhere else in 2012 is touched.

### CN sampling-rate snap (operational fix, 2026-09-16)
CN day-files come back from IRIS as many fragments whose per-fragment sampling rate is recorded as
99.99999 / 100.00001 Hz rather than 100 Hz. When two fragments of the SAME channel disagree, obspy
refuses to merge and the whole station-day is lost. `c_cache_month.py` snaps every fragment rate to the
nearest integer when it is within 1e-4 relative (< 0.01 Hz at 100 Hz, far below one sample over a day)
and counts the snaps in `n_fs_snapped`; a rate further out than that still raises, so a genuinely
different rate is never silently glued together. Verified: a day that needs no snap comes out
byte-identical to the pre-fix cache.

### `daily_noise_<YYYY-MM>.parquet` — per station-day
`netsta`, `year`, `doy`, `date`, `status`, `order`, `fs_native`, `n_fs_snapped`, `e8_Z/N/E` (the daily **median 8-s
window energy per component**, Merlin Step 6 component 2), `frac_nonzero`, `cov_manifest`, `secs`.

### `qc_<YYYY-MM>.csv` / `qc_summary.csv` / `figures/trackC_monthly_qc.png`
Per SITE: quiet-day detections per template-day at both tiers, that site's **running floor λ** (the
median of its earlier months) and `flag_2lambda` (Merlin: flag a site whose quiet-day rate exceeds
2 λ). Per network: unobservable fraction, EQ-flag fraction, HF-ratio median / p90 / fraction above
the cut and its shift against the running median (the response-change alarm), detections per tier,
station-day coverage and per-station noise level.

### `run_<YYYY-MM>.json`, `months_done.csv`, `eq_catalog.csv`
Run metadata (configuration, counts, per-day timings, runtime), the watcher's state table, and the
regional earthquake catalog (PNSN `eq_catalog_pgc.csv` 2009–2017 + USGS FDSN 2009–2026, de-duplicated
within 3 s / 20 km, 26,465 events).

---

## Switching from provisional to final templates

`c_run_month.py` resolves its template set in this order:
1. `--templates /path/to/set.npz` on the command line;
2. a line in `results/TRACK_A_RESULTS.md` of exactly the form

   ```
   FINAL_TEMPLATES: results/a3_templates_final.npz
   ```
   (absolute, or relative to `lfe_extension/`);
3. otherwise `results/frozen_operating_point.json → inputs.templates.path`
   (= `e5r_templates_2012.npz`, tagged `2012_provisional`).

The tag written into every row is the npz basename, so a re-run is self-identifying. **Re-running is
the whole record again**: delete (or move) `results/trackC/months_done.csv` and the month files, then
restart the watcher. Nothing is ever deleted automatically.

Per-epoch templates: when Track A supplies one set per instrument epoch, the switch needs a mapping
file `results/trackC/templates_by_epoch.csv` (`epoch_id, npz`) — `c_run_month.py` currently applies
ONE set to every epoch and stamps `epoch_template` accordingly; the per-epoch loader is the one code
change still owed to Track A's hand-over.
