# STAGE 0 -- audit of the WASHINGTON chain for the zero-fill / no-validity-mask bug class

Agent AB, 2026-09-25.  Static reading of the Washington scripts, plus the cache manifests' headers
(`results/cache_T*.parquet`), plus a **read-only** measurement of the cached waveforms themselves
(`scripts/ab_gapscan.py` -> `results/ab_royer/gapscan_cache20.csv`, `gapscan_hours.csv`).
Nothing in the detection chain was executed and no data product of the chain exists.

**Disclosure.** `ab_gapscan.py` was launched at 03:33Z before the instruction to write code without
running it arrived.  `kill -TERM` was issued against the pid returned by the shell, but that pid was
the backgrounding subshell and the python process survived as an orphan; it ran to completion at
03:41Z.  It is read-only (it opens cached `.npz` files and writes two CSVs), it downloaded nothing and
used no GPU.  Its output is used below **because it corrects the header-based reading in one
important respect**, and the correction is flagged where it applies.

The trigger: in `lfe_south/scripts/t11_mf.py` the normalised cross-correlation divides by
`sd = sqrt(max(M*S2 - S1*S1, 1e-30))` while `t3_cache.py:31-32` zero-fills gaps with no validity
mask, so a constant run gives `sd ~ 1e-15` and a "normalised" cc of 2.4e15 -- 128,922 detections in
one day against Shelly's 71.  Independently, agent G found a station-hour more than half flat giving
`ew.median() = 0` -> `cc = inf` -> `z = NaN` in the north engine.

---

## VERDICT IN ONE PARAGRAPH

**The southern bug itself -- an unbounded cc produced by a clamped denominator -- is NOT present in
the Washington chain.** Where the south *floors* the denominator (`max(var, 1e-30)`, so a numerator
that is only FFT round-off divides by 1e-15 and explodes), Washington sets the denominator to
`+inf` for a rejected window (`torch.where(good, ew.sqrt(), inf)`, `k_auto.py:76`;
`torch.where(e > 1e-4*emed, e.sqrt(), inf)`, `k_grow.py:167-169`), which forces `cc = 0`.  That is a
validity mask expressed through the denominator, and every cc in the Washington engines is
Cauchy-Schwarz bounded by 1.  **But the same bug CLASS is present in three weaker forms** -- a
guard that is RELATIVE to the median and therefore degenerates on a mostly-dead channel; masked
samples entering the median/MAD pool; and a network divisor that is a constant channel count rather
than the per-sample count of channels actually contributing.  All three are real and all three are
quantified below.  Their measured exposure on the recorded Washington runs is small at day resolution
(0.62 % of the 2008 Puget station-days, 0.04 % of the 2010 Olympic ones), but the per-HOUR measurement
on Royer's own bed -- the resolution at which `k_auto` actually takes its MAD -- finds **5x** that
exposure and finds it on **four stations the day-level test calls completely clean**, because a single
dead component inside a live three-component station is invisible to every all-component zero test this
project has used.  **No recorded Washington number is invalidated by this bug class; in particular it
does not explain the 2008 zeros, whose named cause is elsewhere (the 6 s window, the 20 Hz rate and the
wrong station set) and whose direction is opposite.**  Two further defects surfaced on the way and are
worth as much as the bug itself: the 20 Hz cache does keep **exact** 50,000-second zero runs, and the
manifest column written to detect them (`frac_nonzero`) is simply **wrong** -- 1.000000 on a file that
is 41.9 % non-zero.

---

## 1. What the cache actually stores, and why its own gap indicator is blind

`lfe_picker/scripts/o_cache_month.py:70` -- `st.merge(method=1, fill_value=0)` -- then band-passes
at the native rate and decimates.  No mask is written to the npz.

**Exact zero runs DO survive the band-pass and the decimation.**  Measured on the file itself,
`cache20hz/UW.GNW/2010/224.npz` (2010-08-12, a Royer NW10 station-day) is **41.93 % non-zero** on all
three components, exactly matching the download manifest's `coverage_frac` of 0.418, and its longest
all-zero run is **50,174 s**.  So the southern pathology's raw ingredient -- a long exactly-constant
run reaching the correlation -- is present in this cache too; only the guard differs.
*(This corrects a wrong inference I first drew from the manifest header, that filtering smeared the
zeros away.  It does not.)*

**A second, independent defect: the manifest's own gap column is wrong.**  For that same station-day
`results/cache_T*.parquet` records `frac_nonzero = 1.000000` while the file is 0.419.  So the column
built to detect this (`o_cache_month.py:53,92`) does not describe the file on disk and cannot be used
as a gap indicator anywhere.  Cache-wide it reports only 10 of 23,744 station-days below 1.0 and
**none** below 0.5, which is false.  Two columns of the same manifest *are* reliable and did flag the
day: `coverage_frac` (0.418, carried from the download manifest) and `e8_Z/N/E` (0, 0, 0).  No
downstream script reads either.

Royer's own episode station-days, measured on the waveforms (`gapscan_cache20.csv`, 756 of the 833
present in the 20 Hz cache):

| quantity | value |
|---|---|
| station-days with any all-3-component zero sample | **6 (0.8 %)** -- all UW.GNW |
| station-days more than half zero | **1** (UW.GNW 2010-08-12, 58.1 %) |
| station-days with any zero-energy 15 s window | **11** |
| station-days tripping the degenerate guard (`emed == 0` on some component) | **1** |
| **station-hours** tripping the degenerate guard | **58 of 18,138 (0.32 %)** |

## 2. The three defects, file and line

### D1 -- the energy guard is RELATIVE to the median, so it degenerates on a mostly-dead channel
`k_auto.py:69-70`
```
emed = ew.median(dim=1, keepdim=True).values      # per component, over the hour
good = ew > 1e-4 * emed                           # (3, NV)
```
`k_grow.py:166-169` is the same statement over the day, and `o_auto.py:66-76` is a copy of
`k_auto.py`.  When a component is more than half dead, `emed == 0` exactly and the guard becomes
`ew > 0` -- i.e. **no guard at all**: every gap-smear window is admitted.  The cc stays bounded, but
two windows that each straddle a gap EDGE both contain the band-pass impulse response, which is a
clean 1-8 Hz transient, so their normalised cc is genuinely high and they pair as a "detection".
The block-level pre-check `if float(ew.median()) <= 0: continue` (`k_auto.py:66`) does not catch this
because it is the median over all three components flattened, not per component.

**Measured incidence.**  Two views, and the finer one matters more because `k_auto` takes its MAD over
the **hour**, not the day.

*Per station-day, from the manifests' `e8_Z/N/E` (the median 6 s window energy per component, i.e.
`emed` at the rate and window those runs used):* 81 station-days cache-wide carry a component with
`e8 == 0` exactly.

| run | cached station-days in its admitted roster | with a zero-median component | stations |
|---|---|---|---|
| 2008 Puget (52-station probe roster) | 3,867 | **24 (0.62 %)** | YW.FACB (19), PB.B018, XU.N050, XU.BS11 |
| 2010 Olympic (97-station roster) | 6,817 | **3 (0.04 %)** | UW.GNW, UW.LRIV, XG.PA02 |
| Royer NW+SW episodes (the 7+7 named stations) | 756 | 1 (UW.GNW 2010-08-12) | UW.GNW |

*Per station-HOUR, measured directly on Royer's own station-days (`gapscan_hours.csv`):* **58 of
18,138 hours (0.32 %)** trip `emed == 0` on at least one component -- an order of magnitude more
exposure than the day view shows, and **on four more stations**:

| subarray episode | station | hours tripping the guard | day-level `frac_zero` |
|---|---|---|---|
| NW10 | GNW | 14 (2010-08-12) | 0.581 |
| NW08 | PL11 | 11 (2008-05-09) | **0.000** |
| SW08 | PL11 | 11 (2008-05-09) | **0.000** |
| SW08 | S040 | 11 (2008-05-09) | **0.000** |
| SW08 | N060 | 9 (2008-05-09) | **0.000** |
| NW07 / NW08 | GNW | 1 + 1 | 0.022 / 0.032 |

**The PL11 / S040 / N060 rows are the important ones and they are invisible to every test this project
has ever applied.**  Their day has no all-three-component zero sample at all (`frac_zero = 0`), and no
15 s window of zero total energy (`frac_win_zero = 0`), yet for 9-11 hours of it **one or two
individual components are dead for more than half the hour**, so that component's `emed` is exactly 0
and its guard becomes `ew > 0`.  The reject fraction the guard *would* have applied on those hours
(`frac_win_tiny`) runs 0.19-0.55, so between a fifth and a half of that channel's windows should have
been excluded and none were.  2008-05-09 is exactly the day whose `coverage_frac` is 0.75 (PL11, S040)
and 0.84 (N060) -- the gap is real, it is per-component, and only a per-component test sees it.

### D2 -- masked samples enter the median/MAD pool
`k_auto.py:92-100` accumulates the whole `acc` block into `acc_stats` and takes the median and MAD of
everything, including the (row, column) pairs whose denominator was `inf` and whose value is
therefore an exact 0.  `k_grow.py:202-204` does the same per group over the day.  A day or hour with
a substantial gap therefore has a median and a MAD pulled toward zero by a spike of exact zeros, so
`med + 8*MAD` is a **lower** bar than 8 x MAD of the live data would be.  Direction: **more**
detections, all of them in the live part of the record, i.e. a false-alarm inflation that is
invisible in the log (the log prints the threshold, not the pool composition).

### D3 -- the network divisor is a constant channel count, not the per-sample count
`k_auto.py:90` `acc /= float(nch)` with `nch` fixed for the hour; `k_grow.py:200`
`acc = acc / cnt.clamp_min(1.0)` with `cnt` incremented **once per station** (`k_grow.py:199`,
a scalar per group) rather than per sample.  So at a time when 3 of 7 stations are inside a gap, the
network mean is still divided by 7.  Direction: **dilution**, i.e. **missed** detections, and it is
time-varying within a day, which means the statistic is not stationary across the block whose MAD is
being taken.  This is the same arithmetic that made Ducellier bed B worse than bed A.

### D4 -- `k_common.load_day_native` zero-fills with no guard at all
`k_common.py:214-222` builds `x = np.zeros(n)` and inserts each trace, then band-passes; there is no
mask and no downstream energy guard on this path.  Used only by `k_s2_pair.py`, `k_s2c_sweep.py` and
`k_s2d_weights.py` (the `--chain native` sweeps), not by any run whose numbers are on record.

### `o_beam.py` -- the array-beam front end has the relative-guard defect too
`o_beam.py:126` admits a 10-minute member window unless `x.std() <= 0`, then normalises it with
`Z.append(x / x.std())`.  A member whose 10 minutes are mostly gap has a tiny std, is amplified to
unit std, and enters the delay-and-sum beam **with equal weight**.  Only exactly-flat is rejected.
`o_beam.py:272` (`(seg**2).mean() > 0`) is the same exact-zero test.  There is no MAD threshold in
`o_beam.py`, so this contaminates the array trace rather than a detection count.

### `uu_*.py` -- clean for the adopted locator
`uu_stacks.py:60,75,81`, `uu_picks.py:82-90`, `uu_dd.py:47,54` all guard with `if d <= 0: skip` or
`np.where(den > 0, num/den, 0)` and operate on STACKS, never on raw day data, so they never see a
gap directly.  One floor-instead-of-mask remains: `uu_picks.py:68-70` floors the AIC variances at
`1e-30` inside a logarithm, so `log(1e-30) = -69` makes a flat interval the global minimum and pins
the onset there.  It is inert for the adopted locator: the AIC onset (`A_onset`) was measured and
**discarded** (`WA_LOCATION_REBUILD.md` sec.2 and sec.4: epi 16.94 km vs 4.00 km), and `L9_xcdiff`
uses `obs='xc'` with `hcol='h_envH'`.

## 3. What this invalidates -- and what it does not

| number on record | touched? | reasoning |
|---|---|---|
| **2008 Puget zeros** (`WA_2008_RESULTS.md`, `WA_2008Q_RESULTS.md`) | **NO -- and the bug cannot explain them** | D1/D2 push the detection count **up**; the 2008 result was zero families.  Its named cause stands: the 6 s window, 20 Hz rate and 52 substituted stations of `PAPER_METHODS_TABLE.md` M1/M2.  0.62 % of that bed's station-days carry the degenerate guard, which can add false alarms but cannot remove families. |
| **2010 Olympic 61-of-69 coordinate match** | **NO** (day-level exposure 3 in 6,817 = 0.04 %; per-hour unmeasured) | the whole-chain site-shift null received the identical treatment and hit 3 of 69 at 10 km against the real run's 61.  A gap artefact appears in both arms, so it cannot manufacture that separation.  The residual uncertainty is that per-hour exposure on this roster was not measured (sec.4). |
| **280 dedupe-rich 2010 templates** (null maximum 16) | **NO** | same run, same 0.04 % exposure, same null.  An artefact shared by both arms cannot produce a real-vs-null separation of 280 versus 16. |
| **`WA_LOCATION_REBUILD.md` locator calibration (4.0 km vs 27.6 km)** | **NO** | it reads stacks, not day data, and the one floor-instead-of-mask defect is in a variant that was measured and discarded. |
| **anything computed with `k_common.load_day_native`** | not on record | the `--chain native` sweeps only. |

**Therefore: the Washington numbers on record are NOT retracted by this audit, and no re-run of them
is owed on account of this bug class.** What is owed, and what this bed does, is not to carry the
defect forward: the Royer chain reads a cache that stores a real per-sample validity mask
(`ab_cache40.py`), admits a window only if every one of its samples is valid, keeps masked pairs out
of the median/MAD pool, and divides the network sum by the **per-sample** count of contributing
channels (`ab_common.window_stats`, `ab_auto.autocorr_hour`, `ab_grow.matched_filter`).  No relative
energy guard appears anywhere in the `ab_` chain.

## 4. The measurement

`scripts/ab_gapscan.py` (also stage 0 of `scripts/ab_royer_queue.sh`, where it is guarded by its own
output so it will be skipped) opens the cached station-days of Royer's 14 stations over the six
episode windows and writes, per station-day and per hour: `frac_zero`, the longest constant run in
seconds, the number of 15 s windows of exactly zero energy, `emed_zero_comp` (how many components trip
the degenerate guard), and `frac_win_tiny` (the fraction of windows the `1e-4 * emed` guard would
reject).  Outputs, both present: `results/ab_royer/gapscan_cache20.csv` (834 rows) and
`gapscan_hours.csv` (18,138 rows).  It reads the 20 Hz cache and nothing else.

**It only covers Royer's own station-days.**  The per-hour incidence on the 2008 Puget and 2010 Olympic
rosters -- where the numbers on record live -- is still known only at day resolution from `e8_*`.
Given that the hour view found 5x the exposure of the day view on this bed, and found it on four
stations the day view called clean, the honest statement about those rosters is: **day-level exposure
is 0.62 % and 0.04 %, and the per-hour figure is unmeasured and is probably several times larger.**
Extending `ab_gapscan.py` to those two rosters is a one-argument change and is the first thing to run
if anyone wants to tighten section 3's verdict beyond "the null absorbs it".

## 5. Provenance

Read: `scripts/k_common.py`, `k_auto.py`, `k_grow.py`, `o_auto.py`, `o_beam.py`, `uu_common.py`,
`uu_stacks.py`, `uu_picks.py`, `uu_locate.py`, `uu_dd.py`, `d_cache.py`, `d_equiv.py`,
`lfe_picker/scripts/o_cache_month.py`, `lfe_south/scripts/t3_cache.py`, `t11_mf.py`.
Header inspection: `results/cache_T2A..T5.parquet`, `results/site_probe_2008.csv`,
`results/site_probe_2010.csv`, `results/fdsn_stations_40-50N_2003-2014.txt`.
Waveform measurement (read-only): `results/ab_royer/gapscan_cache20.csv`, `gapscan_hours.csv`,
`logs/ab_gapscan.log`, and a direct read of `cache20hz/UW.GNW/2010/224.npz`.
