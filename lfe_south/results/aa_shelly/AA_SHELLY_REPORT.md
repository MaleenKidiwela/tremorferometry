# Shelly (2025) reproduction — the 12.8 M detections, diagnosed and fixed

Agent AA, 2026-09-25.  Companion: `notes/PAPER_METHODS_TABLE.md` §S25 rows (his actual parameters, read from the paper) and `notes/CATALOGUE_SUMMARY.md` §10.

## 0. The fault, measured

`t11_mf.py` normalises the cc by `sd = sqrt(max(M*S2 - S1*S1, 1e-30))` from cumulative sums.  The 20 Hz cache zero-fills gaps (`t3_cache.py:31-32`, obspy `merge(fill_value=0)` + `trim(pad=True, fill_value=0)`) and stores no validity mask, so a gap is an exactly constant run: the true window variance is 0, the cumsum difference goes negative from round-off, the clamp gives sd = 1e-15, and the 'normalised' cc reaches 1e15.

| measurement | value |
|---|---|
| BK.RVIT HHZ, family 23, 2022-05-11: constant-sample fraction | 0.333 |
| same channel, max abs 'normalised' cc | 2.4e15 |
| NC.KIP EHZ same day: constant fraction / max cc | 0.501 / 2.8e10 |
| channel-summed cc series that day: MAD (cannot exceed nchan=20) | 25.7 |
| fraction of samples above 8xMAD | 0.166 |
| detections that ONE day, 27 families, t11_mf.py | 128,922 |
| detections that ONE day, 27 families, fixed | 21 |
| Shelly's own detections that day, 27 families | 71 |

Second fault: `t11_mf.py` sums each station's cc at a common lag, but every station's template is cut at its own offset tt_s+tau_s (`t9_templates.py:93`), which for family 1 spans 1.49-5.85 s.  The stack was therefore incoherent and its time axis was not Shelly's.  After shifting each channel by (tt_s+tau_s-pre) the peak at his published times lands at dt median 0.00 s.

Third fault, in the SCORING: `t8_score.py` used all 61,441 published events as the denominator while the filter ran only on the held-out days, where only 5,021 of his events live — the recall cap was 0.082, so the published '0.026' was 0.32 of the achievable maximum.  `t8_score_aa.py --days-mode source|fam` fixes the denominator.


## Row A — 12 s templates, cache band 1.5-9 Hz, moveout-aligned

* family-days processed **10395**, channels/family-day median **19**, flat-masked fraction median **0.000**
* nmad of the daily statistic: p50 **0.00**, p99 **3.47**, p99.9 **4.64**, day max median **10.5** (medians over family-days)
* detections per family-day at each threshold (mean / median):

| threshold xMAD | mean det/family-day | median | max |
|---|---|---|---|
| 6 | 355.2 | 66 | 13080 |
| 7 | 246.0 | 8 | 10837 |
| 8 | 207.0 | 3 | 10127 |
| 9 | 182.8 | 1 | 10127 |
| 10 | 166.1 | 1 | 10126 |
| 12 | 146.1 | 0 | 10126 |
| 15 | 131.3 | 0 | 10126 |

* **quiet-day false-alarm rate** (family-days where Shelly has NO detection, n=3684): mean **345.32** det/family-day at 9xMAD, vs **360.58** on his active family-days (n=6711); his own rate on active family-days is **4.01**
  * at 6xMAD: quiet **345.32** vs active **360.58** det/family-day
  * at 7xMAD: quiet **240.55** vs active **249.05** det/family-day
  * at 8xMAD: quiet **202.27** vs active **209.67** det/family-day
  * at 9xMAD: quiet **178.53** vs active **185.20** det/family-day
  * at 10xMAD: quiet **162.09** vs active **168.30** det/family-day
  * at 12xMAD: quiet **143.38** vs active **147.60** det/family-day
  * at 15xMAD: quiet **130.28** vs active **131.79** det/family-day

* **separation curve** — our statistic at his 26895 published times on these days: nmad p10 **3.79**, p50 **6.30**, p90 **11.56** (best within +/-2 s), against a noise p99.9 of **4.64**; timing offset dt median **0.00** s
  * fraction of his events above 6xMAD: **0.544**
  * fraction of his events above 7xMAD: **0.408**
  * fraction of his events above 8xMAD: **0.298**
  * fraction of his events above 9xMAD: **0.218**
  * fraction of his events above 10xMAD: **0.161**
  * fraction of his events above 12xMAD: **0.087**

| tmin | his events | ours | matched | pooled recall | precision | fam with >=1 match | fam >=30 % precision |
|---|---|---|---|---|---|---|---|
| 6.0 | 26895 | 3692010 | 14400 | 0.535 | 0.0039 | 27/27 | 0 |
| 7.0 | 26895 | 2557525 | 10838 | 0.403 | 0.0042 | 27/27 | 0 |
| 8.0 | 26895 | 2152215 | 7936 | 0.295 | 0.0037 | 27/27 | 0 |
| 9.0 | 26895 | 1900564 | 5799 | 0.216 | 0.0031 | 27/27 | 0 |
| 10.0 | 26895 | 1726559 | 4291 | 0.160 | 0.0025 | 27/27 | 0 |
| 12.0 | 26895 | 1518761 | 2308 | 0.086 | 0.0015 | 27/27 | 0 |
| 15.0 | 26895 | 1364369 | 999 | 0.037 | 0.0007 | 25/27 | 0 |
| nan | 26895 | 3692010 | 14400 | 0.535 | 0.0039 | 27/27 | 0 |

Per-family rows at every threshold are in `mfA_score.csv` — no family is dropped.

## Row B — 12 s templates, 2-9 Hz (Shelly's band), moveout-aligned

* family-days processed **10395**, channels/family-day median **19**, flat-masked fraction median **0.000**
* nmad of the daily statistic: p50 **0.00**, p99 **3.51**, p99.9 **4.75**, day max median **11.9** (medians over family-days)
* detections per family-day at each threshold (mean / median):

| threshold xMAD | mean det/family-day | median | max |
|---|---|---|---|
| 6 | 416.1 | 119 | 12924 |
| 7 | 256.7 | 21 | 10736 |
| 8 | 206.1 | 7 | 10130 |
| 9 | 180.8 | 3 | 10129 |
| 10 | 165.1 | 2 | 10129 |
| 12 | 146.0 | 0 | 10128 |
| 15 | 129.9 | 0 | 10127 |

* **quiet-day false-alarm rate** (family-days where Shelly has NO detection, n=3684): mean **413.46** det/family-day at 9xMAD, vs **417.50** on his active family-days (n=6711); his own rate on active family-days is **4.01**
  * at 6xMAD: quiet **413.46** vs active **417.50** det/family-day
  * at 7xMAD: quiet **258.43** vs active **255.69** det/family-day
  * at 8xMAD: quiet **209.12** vs active **204.51** det/family-day
  * at 9xMAD: quiet **184.22** vs active **178.88** det/family-day
  * at 10xMAD: quiet **169.14** vs active **162.87** det/family-day
  * at 12xMAD: quiet **152.16** vs active **142.59** det/family-day
  * at 15xMAD: quiet **139.16** vs active **124.82** det/family-day

* **separation curve** — our statistic at his 26895 published times on these days: nmad p10 **4.68**, p50 **7.85**, p90 **13.62** (best within +/-2 s), against a noise p99.9 of **4.75**; timing offset dt median **0.00** s
  * fraction of his events above 6xMAD: **0.740**
  * fraction of his events above 7xMAD: **0.608**
  * fraction of his events above 8xMAD: **0.480**
  * fraction of his events above 9xMAD: **0.368**
  * fraction of his events above 10xMAD: **0.279**
  * fraction of his events above 12xMAD: **0.156**

| tmin | his events | ours | matched | pooled recall | precision | fam with >=1 match | fam >=30 % precision |
|---|---|---|---|---|---|---|---|
| 6.0 | 26895 | 4325023 | 19614 | 0.729 | 0.0045 | 27/27 | 0 |
| 7.0 | 26895 | 2668002 | 16164 | 0.601 | 0.0061 | 27/27 | 0 |
| 8.0 | 26895 | 2142882 | 12803 | 0.476 | 0.0060 | 27/27 | 0 |
| 9.0 | 26895 | 1879120 | 9833 | 0.366 | 0.0052 | 27/27 | 0 |
| 10.0 | 26895 | 1716118 | 7440 | 0.277 | 0.0043 | 27/27 | 0 |
| 12.0 | 26895 | 1517487 | 4153 | 0.154 | 0.0027 | 27/27 | 0 |
| 15.0 | 26895 | 1350316 | 1815 | 0.067 | 0.0013 | 27/27 | 0 |
| nan | 26895 | 4325023 | 19614 | 0.729 | 0.0045 | 27/27 | 0 |

Per-family rows at every threshold are in `mfB_score.csv` — no family is dropped.

## Row C — 6 s templates (Shelly's length), 2-9 Hz, moveout-aligned

* family-days processed **10395**, channels/family-day median **19**, flat-masked fraction median **0.000**
* nmad of the daily statistic: p50 **0.00**, p99 **3.50**, p99.9 **4.72**, day max median **10.8** (medians over family-days)
* detections per family-day at each threshold (mean / median):

| threshold xMAD | mean det/family-day | median | max |
|---|---|---|---|
| 6 | 371.9 | 102 | 12984 |
| 7 | 242.4 | 16 | 10971 |
| 8 | 203.2 | 4 | 10130 |
| 9 | 181.5 | 2 | 10130 |
| 10 | 167.0 | 1 | 10130 |
| 12 | 148.2 | 0 | 10129 |
| 15 | 131.9 | 0 | 10128 |

* **quiet-day false-alarm rate** (family-days where Shelly has NO detection, n=3684): mean **366.29** det/family-day at 9xMAD, vs **374.96** on his active family-days (n=6711); his own rate on active family-days is **4.01**
  * at 6xMAD: quiet **366.29** vs active **374.96** det/family-day
  * at 7xMAD: quiet **245.09** vs active **241.00** det/family-day
  * at 8xMAD: quiet **207.75** vs active **200.62** det/family-day
  * at 9xMAD: quiet **186.56** vs active **178.75** det/family-day
  * at 10xMAD: quiet **172.46** vs active **163.93** det/family-day
  * at 12xMAD: quiet **155.14** vs active **144.33** det/family-day
  * at 15xMAD: quiet **141.49** vs active **126.60** det/family-day

* **separation curve** — our statistic at his 26895 published times on these days: nmad p10 **4.54**, p50 **7.41**, p90 **12.33** (best within +/-2 s), against a noise p99.9 of **4.72**; timing offset dt median **0.00** s
  * fraction of his events above 6xMAD: **0.711**
  * fraction of his events above 7xMAD: **0.562**
  * fraction of his events above 8xMAD: **0.421**
  * fraction of his events above 9xMAD: **0.307**
  * fraction of his events above 10xMAD: **0.220**
  * fraction of his events above 12xMAD: **0.112**

| tmin | his events | ours | matched | pooled recall | precision | fam with >=1 match | fam >=30 % precision |
|---|---|---|---|---|---|---|---|
| 6.0 | 26895 | 3865746 | 18870 | 0.702 | 0.0049 | 27/27 | 0 |
| 7.0 | 26895 | 2520265 | 14957 | 0.556 | 0.0059 | 27/27 | 0 |
| 8.0 | 26895 | 2111752 | 11219 | 0.417 | 0.0053 | 27/27 | 0 |
| 9.0 | 26895 | 1886890 | 8186 | 0.304 | 0.0043 | 27/27 | 0 |
| 10.0 | 26895 | 1735493 | 5858 | 0.218 | 0.0034 | 27/27 | 0 |
| 12.0 | 26895 | 1540140 | 2963 | 0.110 | 0.0019 | 27/27 | 0 |
| 15.0 | 26895 | 1370883 | 1107 | 0.041 | 0.0008 | 26/27 | 0 |
| nan | 26895 | 3865746 | 18870 | 0.702 | 0.0049 | 27/27 | 0 |

Per-family rows at every threshold are in `mfC_score.csv` — no family is dropped.

## Row N — ABLATION: normalisation fixed but NOT moveout-aligned (t11_mf.py stacking)

* family-days processed **378**, channels/family-day median **12**, flat-masked fraction median **0.000**
* nmad of the daily statistic: p50 **0.00**, p99 **3.50**, p99.9 **4.72**, day max median **8.3** (medians over family-days)
* detections per family-day at each threshold (mean / median):

| threshold xMAD | mean det/family-day | median | max |
|---|---|---|---|
| 6 | 216.5 | 92 | 3214 |
| 7 | 40.9 | 8 | 1359 |
| 8 | 7.2 | 1 | 372 |
| 9 | 1.3 | 0 | 69 |
| 10 | 0.4 | 0 | 9 |
| 12 | 0.3 | 0 | 5 |
| 15 | 0.3 | 0 | 4 |

* **quiet-day false-alarm rate** (family-days where Shelly has NO detection, n=152): mean **181.57** det/family-day at 9xMAD, vs **240.01** on his active family-days (n=226); his own rate on active family-days is **2.70**
  * at 6xMAD: quiet **181.57** vs active **240.01** det/family-day
  * at 7xMAD: quiet **31.29** vs active **47.35** det/family-day
  * at 8xMAD: quiet **5.44** vs active **8.38** det/family-day
  * at 9xMAD: quiet **1.18** vs active **1.32** det/family-day
  * at 10xMAD: quiet **0.48** vs active **0.33** det/family-day
  * at 12xMAD: quiet **0.43** vs active **0.19** det/family-day
  * at 15xMAD: quiet **0.42** vs active **0.19** det/family-day

* **separation curve** — our statistic at his 611 published times on these days: nmad p10 **3.36**, p50 **4.92**, p90 **7.30** (best within +/-2 s), against a noise p99.9 of **4.72**; timing offset dt median **0.30** s
  * fraction of his events above 6xMAD: **0.237**
  * fraction of his events above 7xMAD: **0.124**
  * fraction of his events above 8xMAD: **0.046**
  * fraction of his events above 9xMAD: **0.018**
  * fraction of his events above 10xMAD: **0.008**
  * fraction of his events above 12xMAD: **0.000**

| tmin | his events | ours | matched | pooled recall | precision | fam with >=1 match | fam >=30 % precision |
|---|---|---|---|---|---|---|---|
| 6.0 | 611 | 81842 | 133 | 0.218 | 0.0016 | 20/27 | 0 |
| 7.0 | 611 | 15458 | 71 | 0.116 | 0.0046 | 16/27 | 0 |
| 8.0 | 611 | 2722 | 27 | 0.044 | 0.0099 | 13/27 | 0 |
| 9.0 | 611 | 477 | 11 | 0.018 | 0.0231 | 6/27 | 2 |
| 10.0 | 611 | 148 | 5 | 0.008 | 0.0338 | 4/27 | 1 |
| 12.0 | 611 | 109 | 0 | 0.000 | 0.0000 | 0/27 | 0 |
| 15.0 | 611 | 108 | 0 | 0.000 | 0.0000 | 0/27 | 0 |
| nan | 611 | 81842 | 133 | 0.218 | 0.0016 | 20/27 | 0 |

Per-family rows at every threshold are in `mfN_score.csv` — no family is dropped.

## Row D — row A templates re-stacked on our own build-day detections (Shelly's iteration)

* family-days processed **10395**, channels/family-day median **19**, flat-masked fraction median **0.000**
* nmad of the daily statistic: p50 **0.00**, p99 **3.48**, p99.9 **4.64**, day max median **7.8** (medians over family-days)
* detections per family-day at each threshold (mean / median):

| threshold xMAD | mean det/family-day | median | max |
|---|---|---|---|
| 6 | 274.3 | 48 | 13128 |
| 7 | 192.8 | 3 | 10527 |
| 8 | 165.0 | 0 | 10127 |
| 9 | 147.9 | 0 | 10127 |
| 10 | 136.3 | 0 | 10127 |
| 12 | 121.6 | 0 | 10126 |
| 15 | 106.2 | 0 | 10120 |

* **quiet-day false-alarm rate** (family-days where Shelly has NO detection, n=3684): mean **275.12** det/family-day at 9xMAD, vs **273.91** on his active family-days (n=6711); his own rate on active family-days is **4.01**
  * at 6xMAD: quiet **275.12** vs active **273.91** det/family-day
  * at 7xMAD: quiet **194.63** vs active **191.78** det/family-day
  * at 8xMAD: quiet **165.94** vs active **164.46** det/family-day
  * at 9xMAD: quiet **147.43** vs active **148.08** det/family-day
  * at 10xMAD: quiet **134.87** vs active **137.04** det/family-day
  * at 12xMAD: quiet **119.50** vs active **122.76** det/family-day
  * at 15xMAD: quiet **105.74** vs active **106.51** det/family-day

* **separation curve** — our statistic at his 26895 published times on these days: nmad p10 **2.13**, p50 **2.94**, p90 **4.06** (best within +/-2 s), against a noise p99.9 of **4.64**; timing offset dt median **0.00** s
  * fraction of his events above 6xMAD: **0.012**
  * fraction of his events above 7xMAD: **0.010**
  * fraction of his events above 8xMAD: **0.009**
  * fraction of his events above 9xMAD: **0.009**
  * fraction of his events above 10xMAD: **0.008**
  * fraction of his events above 12xMAD: **0.007**

| tmin | his events | ours | matched | pooled recall | precision | fam with >=1 match | fam >=30 % precision |
|---|---|---|---|---|---|---|---|
| 6.0 | 26895 | 2851762 | 276 | 0.010 | 0.0001 | 27/27 | 0 |
| 7.0 | 26895 | 2004025 | 239 | 0.009 | 0.0001 | 26/27 | 0 |
| 8.0 | 26895 | 1714991 | 211 | 0.008 | 0.0001 | 25/27 | 0 |
| 9.0 | 26895 | 1536917 | 197 | 0.007 | 0.0001 | 25/27 | 0 |
| 10.0 | 26895 | 1416542 | 181 | 0.007 | 0.0001 | 24/27 | 0 |
| 12.0 | 26895 | 1264114 | 165 | 0.006 | 0.0001 | 22/27 | 0 |
| 15.0 | 26895 | 1104287 | 146 | 0.005 | 0.0001 | 22/27 | 0 |
| nan | 26895 | 2851762 | 276 | 0.010 | 0.0001 | 27/27 | 0 |

Per-family rows at every threshold are in `mfD_score.csv` — no family is dropped.
