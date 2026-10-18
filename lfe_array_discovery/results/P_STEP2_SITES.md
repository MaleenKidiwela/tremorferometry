# P_STEP2_SITES — delivery, array-era site table, site probe, coverage and withheld sets
Agent D, 2026-09-20. Covers Merlin Step 2 (metadata/site table), Step 4 (array-era site probe)
and the Step 0 inputs (active-cell coverage, withheld sets), plus the rational-chain
equivalence check and the sealed-set verification. Every number below is computed from
`results/download_manifest.csv`, `results/channels_2006_2012.csv` and `results/station_flags.csv`;
scripts `scripts/d_{report,sites,probe,coverage,equiv}.py`.

## 0. Headline findings (read these first)
1. **ZK is vertical-only in the DATA.** Its StationXML declares EH1/EH2/EHZ, but 1,604 of 1,642
   delivered ZK station-days contain a single channel (EHZ); only 27 carry 3 channels. Merlin's
   Step-2 failure signature therefore fires — for **ZK, not for XG**. The 2008 dense array
   contributes a Z-only site; all 3-C counts below are DELIVERED 3-C, never declared 3-C.
2. **XG = 7 sites at 1 km, 50 Hz, rational chain — confirmed**, and XG is genuinely 3-C in the
   data (9,588 of 9,588 delivered station-days carry 3 channels for the 102 3-C stations).
3. **The XG SL array died 2010-06-29**, before the 2010-07-20..09-10 window. 8 named arrays →
   7 alive in the ETS window; SL only appears in the June quiet days.
4. **XG.CL is admitted 0/16** by the probe (median tremor contrast 1.9 dB): the easternmost XG
   array does not register N-Olympic tremor. Six XG arrays clear Merlin's ≥5-admitted-members
   bar (BH 9, DR 11, GC 13, LC 13, PA 7, TB 12) — exactly his acceptance threshold.
5. **Merlin's withheld rule is infeasible as written.** Nowhere do 3 permanent sites sit within
   60 km of ≥80 % of a sub-region's cells (N-Olympic: 2; everywhere else: 1). A 60 km radius
   cannot cover an ~150 km-long sub-region. Step 0 has to choose between visibility and
   coverage; both options are costed in section 5.
6. Rational twin vs frozen chain: **cc ≥ 0.99985 on every XG component, 1.000000 on PB** →
   equivalent at the 0.999 bar. Sealed set: **5,062/5,062 files hashed, zero cache entries and
   zero unsealed manifest rows inside 2011-07-20..09-10 → SEALED OK.**

## 1. Delivery by tranche
(T6 = ±3 min event windows, so 'station-days' are event-windows and there is no plan count;
T7 = the extension tail, still streaming.)

```
tranche  planned  attempted  data  nodata  error  linked    gib  frozen  rational  n_stations  delivery_frac  median_cov
    T2A   1748.0       1748  1636     112      0     287  25.78    1271       477          33         0.9359         1.0
    T2B    640.0        640   592      48      0       0   9.49     460       180          32         0.9250         1.0
    T2C   5651.0       5651  5616      35      0       0  50.27       0      5651          91         0.9938         1.0
    T2D   1606.0       1606  1523      83      0       0  22.56       0      1606          25         0.9483         1.0
     T3   5062.0       5062  4453     609      0     317  55.25    1325      3737         122         0.8797         1.0
     T4   8218.0       8218  7380     838      0     234 115.26    3968      4250         183         0.8980         1.0
     T5   5019.0       5019  4530     489      0       0  78.73     653      4366          68         0.9026         1.0
     T6      NaN      12951  9143    3808      0       0   0.65    5097      7854         384         0.7060         NaN
     T7 113529.0        717   566     151      0       0   9.42     662        55          13         0.7894         1.0
```

Per-network delivery fraction:

```
tranche    T2A    T2B    T2C    T2D     T3     T4     T5     T6     T7
net                                                                   
CC         NaN    NaN    NaN    NaN    NaN    NaN    NaN  0.706    NaN
CN       1.000  1.000    NaN    NaN  1.000    NaN    NaN  0.412  0.000
PB       0.979  0.963    NaN    NaN  0.961  0.752  0.850  0.871  0.960
TA         NaN    NaN    NaN    NaN    NaN    NaN    NaN  0.974    NaN
US         NaN    NaN    NaN    NaN    NaN    NaN    NaN  0.930    NaN
UW       0.844  1.000    NaN    NaN  1.000  0.994  0.993  0.874  0.827
XG       0.863  0.788  0.994    NaN  0.805    NaN    NaN  0.464    NaN
XU         NaN    NaN    NaN  0.948  0.955  0.935  0.907  0.691    NaN
YW         NaN    NaN    NaN    NaN    NaN  0.872    NaN  0.904    NaN
Z3         NaN    NaN    NaN    NaN    NaN    NaN    NaN  0.426    NaN
ZK         NaN    NaN    NaN    NaN    NaN  0.914    NaN  0.448    NaN
```

Delivered channel counts (the 3-C audit that produced finding 1):

```
nc    1.0  2.0    3.0
net                  
CC      2    2    494
CN      0    0    527
PB      0    0   8776
TA      0    0    643
US      0    0    120
UW    506   78   2401
XG      0    0   9588
XU      1   17  13328
YW      0    0    704
Z3      0    0     52
ZK   1604   11     27
```

**Stations below 70 % delivery in their window: 57 station-tranche entries.** Full list in
`results/delivery_low_stations.csv`; the ones that matter:

```
tranche net   sta  planned  data  delivery_frac
    T2A  XG  SL05       53     0          0.000
    T2A  UW  DOSE       52    19          0.365
    T2A  PB  B203       53    33          0.623
    T2B  XG  SL05       20     0          0.000
    T2B  PB  B202       20     6          0.300
    T2B  XG  LC02       20     6          0.300
    T2C  XG  BH11       59    27          0.458
    T2D  XU  BS11       73    42          0.575
    T2D  XU  BS23       73    42          0.575
     T4  PB  B202       77     0          0.000
     T4  PB  B203       77     0          0.000
     T4  PB  B204       77     0          0.000
     T4  XU  N100       77     0          0.000
     T4  ZK H06BL       50     1          0.020
     T4  ZK  K06L       50     1          0.020
     T4  XU  N080       77     4          0.052
     T4  YW  FACF       77    20          0.260
     T4  XU  PL03       77    27          0.351
     T4  XU  N010       77    28          0.364
     T4  ZK  I05B       19     8          0.421
     T4  ZK  I04B        2     1          0.500
     T4  ZK  J04B        2     1          0.500
     T4  ZK  K04B        2     1          0.500
     T4  ZK  L04B        2     1          0.500
     T4  PB  B020       68    38          0.559
     T4  ZK  J08L       26    16          0.615
     T4  ZK   J08       19    12          0.632
     T4  XU  E010       77    49          0.636
     T4  PB  B023       35    23          0.657
     T4  ZK I05BL       26    18          0.692
     T5  XU  W030       74    13          0.176
     T5  XU  PL03       74    21          0.284
     T5  XU  W090       74    23          0.311
     T5  XU  W060       74    25          0.338
     T5  PB  B013       61    37          0.607
     T5  XU  E010       74    45          0.608
     T5  PB  B003       74    46          0.622
     T5  XU  S060       74    49          0.662
```

* `XG.SL05` 0/53 in T2A is not a fetch failure — the whole SL array stopped on 2010-06-29
  (finding 3). The GC+SL site still delivers 53/53 days through its 13 GC members.
* `PB.B202/B203/B204` 0/77 in T4: no PB data at Mount St Helens in 2008 for those boreholes.
* The ZK entries are the array's real duty cycle: ZK delivered only 2008-03-05, 03-11 and
  2008-05-05..05-22 — i.e. **the ETS itself and almost nothing else** (section 3).

## 2. Step 2 — the array-era site table
`results/sites_array_era_final.csv` (281 sites, 1 km single-linkage over every box station with a
2006–2012 epoch): site_id, members, members_3c, lat, lon, Slab2 depth below the site, fs,
n_components, sensor + sensor class, chain, flags `ehz_only`/`rational`/`short_period`, and for
each array-era window the delivered `days_any`, `days_3c`, `members_3c` and
`members_declared3c_zonly_data`.

**3-C SITES per window per network — this replaces PLAN.md §1's station counts:**

```
window       w2007_01  w2008_05  w2010_08  w2011_08
net3c                                              
CN                  0         0         2         2
PB                  4        11        14        14
PB,TA,UW            1         1         1         1
PB,XU               1         1         1         1
TA,UW,XU,Z3         0         1         1         1
UW                  1         1         1         1
XG                  0         0         7         7
XU                 41        39         3         3
XU,Z3               1         1         1         1
XU,Z3,ZK            1         1         1         1
YW                  0         6         0         0
```

Totals and sub-regions (3-C sites):

| set | 2007-01 | 2008-05 | 2010-08 | 2011-08 |
|---|---|---|---|---|
| whole box | 50 | 62 | 32 | 32 |
| N-Olympic 47.8–48.3 | 9 | 11 | 14 | 14 |
| S-Olympic 47.3–47.8 | 6 | 9 | 6 | 6 |
| Puget/S-WA 46.2–47.3 | 35 | 40 | 7 | 7 |

Delivered 3-C STATIONS per window: {'w2007_01': 67, 'w2008_05': 82, 'w2010_08': 131, 'w2011_08': 120}. Declared-3-C but Z-only in the archive: {'w2007_01': 0, 'w2008_05': 82, 'w2010_08': 0, 'w2011_08': 0}
(all 82 of them ZK, in 2008-05).

**XG confirmation.** 102 of XG's 188 box stations are 3-C (the other 86 are EHZ-only and are
excluded from every sum). All 102 are **SH?/EH? at 50 Hz → rational chain**. At the project's
1 km merge they form **7 sites** (GC and SL chain together through a <1 km link); by deployment
name there are 8 arrays, of which SL is dead in the ETS window:

```
site_id       lat         lon  n_members_3c fs    chain  w2010_08_members_3c  w2010_08_days_3c  w2011_08_members_3c
   S242 47.935244 -123.049838            19 50 rational                   13                53                   13
   S246 47.973112 -123.138540            14 50 rational                   12                53                   11
   S248 48.005519 -123.085057            11 50 rational                   10                53                   10
   S249 48.005657 -123.313276            11 50 rational                   11                53                   10
   S256 48.054878 -123.464704            13 50 rational                   10                53                    9
   S257 48.055400 -123.210930            14 50 rational                   13                53                   10
   S260 48.068733 -122.970047            20 50 rational                   16                53                   13
```

ZK's 97 stations are one site (S245, 47.957/−122.928, shared with XU/Z3/XG members); in 2008-05
that site has 8 delivered-3-C members and 82 declared-3-C-but-Z-only ones.
## 3. Step 4 — array-era site probe (P4-0 statistic, native rate)
Statistic reused verbatim from `lfe_picker/scripts/u_fartremor.hf_one`: native sample rate,
4 disjoint 1.5 h blocks (00/06/12/18 UTC), zero-phase Butterworth-4 in 2–8 Hz and in the 10–19 Hz
control, median over the day of the mean square in 6 s windows, 3-C median. The contrast is a
per-station ratio of the same statistic on two day sets, so the instrument response and the
digitiser gain cancel exactly — response correction cannot change a single dB here, and the
channel epoch is the same on both day sets for every station used.

### 3.1 N-Olympic 2010 (`results/site_probe_2010.csv`)
Tremor days (PNSN ≥100 epicentres within 40 km of the XG array centroid 48.014/-123.176): **2010-08-15, 2010-08-16, 2010-08-18, 2010-08-19, 2010-08-20, 2010-08-21** — the expected
2010-08-16..08-28 band, shifted two days earlier by the epicentre count.
Quiet days (6 best-delivered from the frozen T2B catalog-quiet set): 2010-07-14, 2010-07-15, 2010-07-16, 2010-07-17, 2010-07-18, 2010-07-19.
Far days (S-WA-only tremor): 2010-07-25, 2010-09-03, 2010-09-04, 2010-09-05, 2010-09-06, 2010-09-07.
Admission: ≥4/6 days ≥+3 dB AND local−far ≥+3 dB → **97 of 131 stations admitted**.

```
array net  n_probed  n_admitted  median_contrast_db  median_local_minus_far_db  median_hf_db  flag_lt5
   CN  CN         2           1                3.93                      -2.66          1.00      True
   PB  PB        19          10                3.19                       5.22          0.54     False
   UW  UW         3           2                6.53                       8.14          1.19      True
XG.BH  XG        10           9               13.64                      13.42          6.66     False
XG.CL  XG        16           0                1.91                       6.83         -1.05      True
XG.DR  XG        11          11               15.43                      15.29          4.49     False
XG.GC  XG        13          13               13.07                      12.86          3.93     False
XG.LC  XG        13          13               11.11                      10.60          1.26     False
XG.PA  XG        10           7                5.62                       9.93          3.89     False
XG.TB  XG        12          12               12.48                      12.24          1.91     False
   XU  XU        22          19               10.75                      11.34          3.64     False
```

* **XG.CL fails outright (0/16 admitted, median contrast 1.9 dB)** although its local−far is
  +6.8 dB — CL sits at −122.97, the far eastern end, and simply does not see this source region
  at 2–8 Hz. Flag it: an array with <5 admitted members.
* **XG.SL does not appear at all** — dead since 2010-06-29.
* Six XG arrays have ≥5 admitted members (BH 9, DR 11, GC 13, LC 13, PA 7, TB 12), which is
  exactly Merlin's Step-4 acceptance bar ("≥6 XG sites with ≥5 admitted members each").
* **CN.PFB is rejected on local−far (−10.6 dB)**: the 'far' S-WA days 2010-09-03..07 are louder at
  PFB than the N-Olympic days, i.e. that episode was not far from Vancouver Island. The far
  clause is doing its job but its day set is only 'far' for stations south of the Strait.
* **UW.DOSE is UNTESTABLE** — it was installed 2010-07-21, after every quiet day in the frozen
  T2B set, so it has no baseline (n_quiet = 0).
* PB: 10 of 19 admitted. Admitted with 6/6 days: B001 (10.5 dB), B014 (9.3), B013 (6.0), B943 (6.5).
  The Mount St Helens boreholes (B201–B204) and the south-coast ones (B019/B020/B023) are
  correctly rejected (contrast ≤ 0.5 dB).

### 3.2 Puget/S-Olympic 2008 (`results/site_probe_2008.csv`)
No tremor catalog exists before 2010-01-06, so the tremor days come from the **array's own hourly
energy flag** (Merlin Q3 fallback), computed on the 20 Hz cache: per station the hourly median 6 s
window power in the cache band, divided by that station's median quiet-day hourly value; the
network median of that ratio per hour; an hour is hot above the p95 of the same series on the
literature-quiet days (p95 = 2.91); a day is a tremor day at ≥6 hot hours.
**It is an energy flag, not a catalog** (`results/hour_energy_flag_2008.csv`).

The flag is a strong independent check on itself: it fires on 2008-05-05..05-23 plus 05-29 and
06-10 and nowhere else — i.e. it reconstructs the published 2008-05 Washington ETS from the
waveforms alone, with 0–2 hot hours on every literature-quiet day except 2008-03-11 (7 hot hours,
which is probably real tremor on a day the literature calls quiet — excluded from the baseline).

Tremor days used: **2008-05-06, 2008-05-07, 2008-05-15, 2008-05-16, 2008-05-21, 2008-05-22**. Quiet days: 2008-02-14, 2008-03-24, 2008-03-31, 2008-07-05, 2008-07-19, 2008-07-28
plus the in-window fallback set 2008-03-05, 2008-04-20, 2008-04-26, 2008-04-27, 2008-05-25, 2008-06-01, 2008-06-08.
The fallback was forced: ZK, Z3 and YW have little or no data in the 2008-02/03 and 2008-07
literature-quiet pools, so the 6 quietest in-window days by the flag (0 hot hours, lowest
network-median ratio) are used as their baseline; every station records which baseline it used.
**There are no far-tremor days in 2008** (no catalog) → the local−far clause is UNTESTABLE and
admission rests on the contrast clause alone: **52 of 82 stations with a usable baseline admitted**.

```
array net  n_probed  n_admitted  median_contrast_db  median_hf_db  flag_lt5
   PB  PB        16          10                5.10          1.95     False
   UW  UW         1           0                1.67          4.66      True
   XU  XU        57          40                5.35          3.38     False
   YW  YW         6           2                2.50          1.39      True
   ZK  ZK         2           0                 NaN           NaN      True
```

* **ZK: 471 of 479 attempted station-days fail on channel order — one channel (EHZ) only.**
  Only 2 ZK stations have 3 components and neither clears the bar. The 2008 dense array gives
  a Z-only site, exactly Merlin's Step-2 failure signature. Flag: <5 admitted members.
* **UW: only GNW is 3-C in 2008 and it fails (1.7 dB)** — HDW is EHZ/EH2 until 2022-09, so the
  2008 window has no admitted UW station at all.
* YW 2/6, XU 40/57, PB 10/16. XU carries the 2008 window.

## 4. Coverage with admitted sites
Active cells (Merlin's definition): Slab2 interface points 28–45 km inside the box with ≥30 PNSN
2010–2011 epicentres within 8 km → 403 active cells in the box. **There is no 2008 tremor catalog,
so the Puget/S-Olympic 2008 rows use the SAME 2010–11 cells** — stated, not hidden; if the 2008
ETS lit up different cells, the 2008 coverage numbers are for the wrong target.

```
        setting subregion  n_cells  n_sites_available  n_sites_admitted  pct_cells_ge5_sites_60km_all  pct_cells_ge5_sites_60km  pct_cells_ge5_sites_30km  median_sites_60km  n_withheld_candidates
 N-Olympic 2010 N-Olympic      130                 32                20                          98.5                      94.6                      50.8               10.0                      9
 N-Olympic 2010 S-Olympic       84                 32                20                          91.7                      83.3                      13.1                7.0                      9
Puget-SOly 2008 Puget-SWA      146                 62                37                          99.3                      84.2                      44.5                9.0                      8
Puget-SOly 2008 S-Olympic       84                 62                37                          98.8                      94.0                      14.3                7.0                      8
```

`pct_cells_ge5_sites_60km_all` uses every delivered-3-C site; the other two columns use ADMITTED
sites only. N-Olympic 2010 keeps **94.6 % of cells with ≥5 admitted sites within 60 km hypocentral**
(median 10 sites) — Merlin's 100 % becomes 94.6 % once the probe and the delivered-3-C audit are
applied, which is still the best geometry this project has had. The 30 km epicentral figure
(50.8 %) is well below his 75 % estimate.

## 5. Withheld sets for the Step 0 freeze
**Merlin's rule (3 permanent sites within 60 km of ≥80 % of the sub-region's cells) is infeasible
everywhere.** Sites meeting the 80 % clause: N-Olympic 2, S-Olympic 1, Puget/S-WA 1. A 60 km
radius cannot reach 80 % of a sub-region that is ~150 km long, so the clause cannot be satisfied
by any three stations. Full ranked candidate table in `results/withheld_candidates.csv`:

```
        setting subregion  meets_80pct site_id                                    site        nets                         members  frac_cells_60km  pct_60km_without_it  loss_pct  needed_by_extension
 N-Olympic 2010 N-Olympic         True    S254                                    B001          PB                         PB.B001            0.869                 94.6       0.0                 True
 N-Olympic 2010 N-Olympic        False    S229                          B013+B943+W030       PB,XU         PB.B013;PB.B943;XU.W030            0.415                 94.6       0.0                False
 N-Olympic 2010 N-Olympic        False    S170                                     GNW          UW                          UW.GNW            0.085                 94.6       0.0                False
 N-Olympic 2010 N-Olympic        False    S158                                    B014          PB                         PB.B014            0.046                 94.6       0.0                 True
 N-Olympic 2010 N-Olympic        False    S067                                    B018          PB                         PB.B018            0.000                 94.6       0.0                False
 N-Olympic 2010 N-Olympic        False    S073                                    B017          PB                         PB.B017            0.000                 94.6       0.0                 True
 N-Olympic 2010 N-Olympic        False    S275                                     VGZ          CN                          CN.VGZ            0.500                 91.5       3.1                 True
 N-Olympic 2010 N-Olympic         True    S258                B005+B006+B007+B04A+LRIV    PB,TA,UW PB.B005;PB.B006;PB.B007;UW.LRIV            0.877                 89.2       5.4                 True
 N-Olympic 2010 N-Olympic        False    S259                                    B003          PB                         PB.B003            0.362                 89.2       5.4                 True
 N-Olympic 2010 S-Olympic        False    S254                                    B001          PB                         PB.B001            0.262                 83.3       0.0                 True
 N-Olympic 2010 S-Olympic        False    S258                B005+B006+B007+B04A+LRIV    PB,TA,UW PB.B005;PB.B006;PB.B007;UW.LRIV            0.190                 83.3       0.0                 True
 N-Olympic 2010 S-Olympic        False    S259                                    B003          PB                         PB.B003            0.000                 83.3       0.0                 True
 N-Olympic 2010 S-Olympic        False    S275                                     VGZ          CN                          CN.VGZ            0.000                 83.3       0.0                 True
 N-Olympic 2010 S-Olympic        False    S229                          B013+B943+W030       PB,XU         PB.B013;PB.B943;XU.W030            0.607                 82.1       1.2                False
 N-Olympic 2010 S-Olympic        False    S073                                    B017          PB                         PB.B017            0.048                 82.1       1.2                 True
 N-Olympic 2010 S-Olympic        False    S158                                    B014          PB                         PB.B014            0.357                 73.8       9.5                 True
 N-Olympic 2010 S-Olympic        False    S067                                    B018          PB                         PB.B018            0.179                 71.4      11.9                False
 N-Olympic 2010 S-Olympic         True    S170                                     GNW          UW                          UW.GNW            0.857                 61.9      21.4                False
Puget-SOly 2008 Puget-SWA        False    S026                                    B201          PB                         PB.B201            0.014                 84.2       0.0                False
Puget-SOly 2008 Puget-SWA        False    S212 C001+C002+C003+C004+C005+C04A+C04A+DOSE TA,UW,XU,Z3                         XU.C04A            0.000                 84.2       0.0                False
Puget-SOly 2008 Puget-SWA        False    S229                          B013+B943+W030       PB,XU                 PB.B013;PB.B943            0.000                 84.2       0.0                False
Puget-SOly 2008 Puget-SWA        False    S254                                    B001          PB                         PB.B001            0.000                 84.2       0.0                 True
Puget-SOly 2008 Puget-SWA        False    S258                B005+B006+B007+B04A+LRIV    PB,TA,UW         PB.B005;PB.B006;PB.B007            0.000                 84.2       0.0                 True
Puget-SOly 2008 Puget-SWA        False    S259                                    B003          PB                         PB.B003            0.000                 84.2       0.0                 True
Puget-SOly 2008 Puget-SWA        False    S070                                    B941          PB                         PB.B941            0.192                 83.6       0.7                False
Puget-SOly 2008 Puget-SWA         True    S067                                    B018          PB                         PB.B018            0.856                 81.5       2.7                False
Puget-SOly 2008 S-Olympic        False    S254                                    B001          PB                         PB.B001            0.262                 94.0       0.0                 True
Puget-SOly 2008 S-Olympic        False    S258                B005+B006+B007+B04A+LRIV    PB,TA,UW         PB.B005;PB.B006;PB.B007            0.190                 94.0       0.0                 True
Puget-SOly 2008 S-Olympic        False    S067                                    B018          PB                         PB.B018            0.179                 94.0       0.0                False
Puget-SOly 2008 S-Olympic        False    S026                                    B201          PB                         PB.B201            0.000                 94.0       0.0                False
Puget-SOly 2008 S-Olympic        False    S259                                    B003          PB                         PB.B003            0.000                 94.0       0.0                 True
Puget-SOly 2008 S-Olympic        False    S070                                    B941          PB                         PB.B941            0.012                 92.9       1.2                False
Puget-SOly 2008 S-Olympic        False    S229                          B013+B943+W030       PB,XU                 PB.B013;PB.B943            0.607                 84.5       9.5                False
Puget-SOly 2008 S-Olympic         True    S212 C001+C002+C003+C004+C005+C04A+C04A+DOSE TA,UW,XU,Z3                         XU.C04A            0.821                 82.1      11.9                False
```

Two costed options (`results/withheld_proposal.json`):

**Option A — keep every extension permanent in the discovery sum** (Merlin's stated preference).
* N-Olympic 2010|N-Olympic: withhold ['PB.B013;PB.B943;XU.W030', 'UW.GNW', 'PB.B018'], cell visibility [0.415, 0.085, 0.0];
  coverage ≥5 sites/60 km 94.6 % → 94.6 % (30 km: 50.8 %).
* N-Olympic 2010|S-Olympic: withhold ['PB.B013;PB.B943;XU.W030', 'PB.B018', 'UW.GNW'], cell visibility [0.607, 0.179, 0.857];
  coverage ≥5 sites/60 km 83.3 % → 51.2 % (30 km: 10.7 %).
* Puget-SOly 2008|Puget-SWA: withhold ['PB.B201', 'XU.C04A', 'PB.B018'], cell visibility [0.014, 0.0, 0.856];
  coverage ≥5 sites/60 km 84.2 % → 81.5 % (30 km: 39.7 %).
* Puget-SOly 2008|S-Olympic: withhold ['PB.B018', 'PB.B941', 'XU.C04A'], cell visibility [0.179, 0.012, 0.821];
  coverage ≥5 sites/60 km 94.0 % → 81.0 % (30 km: 9.5 %).
  For N-Olympic this costs **nothing** in coverage (94.6 → 94.6 %) but two of the three withheld
  sites see almost no cells (GNW 8.5 %, B018 0 %): **Gate 1 at B018 is vacuous for N-Olympic**, and
  only the B013+B943 cluster (41.5 %) is a real test. Option A cannot supply Merlin's "Gate 1 PASS
  at ≥2 withheld sites" for N-Olympic.

**Option B — allow an extension permanent to be withheld from the DISCOVERY sum** (it can still
carry its own templates in the extension, since it is only barred from the discovery stacks).
* N-Olympic 2010|N-Olympic: withhold ['PB.B001', 'CN.VGZ', 'PB.B005;PB.B006;PB.B007;UW.LRIV'], visibility [0.869, 0.5, 0.877] (sum 2.246), extension-needed [True, True, True]; coverage 94.6 % → 83.1 % (30 km: 43.1 %).
* N-Olympic 2010|S-Olympic: withhold ['PB.B001', 'PB.B005;PB.B006;PB.B007;UW.LRIV', 'PB.B013;PB.B943;XU.W030'], visibility [0.262, 0.19, 0.607] (sum 1.059), extension-needed [True, True, False]; coverage 83.3 % → 82.1 % (30 km: 10.7 %).
* Puget-SOly 2008|Puget-SWA: withhold ['PB.B201', 'XU.C04A', 'PB.B018'], visibility [0.014, 0.0, 0.856] (sum 0.87), extension-needed [False, False, False]; coverage 84.2 % → 81.5 % (30 km: 39.7 %).
* Puget-SOly 2008|S-Olympic: withhold ['PB.B001', 'PB.B018', 'PB.B013;PB.B943'], visibility [0.262, 0.179, 0.607] (sum 1.048), extension-needed [True, False, False]; coverage 94.0 % → 82.1 % (30 km: 9.5 %).
  For N-Olympic this buys two genuinely-testing withheld sites (B001 86.9 %, B005-cluster 87.7 %,
  VGZ 50.0 %) for an 11.5-point coverage cost (94.6 → 83.1 %). **Recommended for N-Olympic**, since
  a Gate-1 clause that cannot fail is worth less than 11 points of redundancy.

**Permanents that must stay available to the extension** (they are the whole post-2011 network):
PB B001 B003 B004 B005 B006 B007 B014 B017, UW LRIV, UW DOSE, CN VGZ, CN PFB. Being withheld from
a discovery sum does not remove a station from the extension; it only bars it from the stacks that
produced the family it is later asked to confirm.

**Puget/S-Olympic 2008 has essentially no withheld capacity**: B018 (85.6 % visibility) is the only
admitted permanent site that sees the Puget cells; B941 sees 19.2 %, everything else 0–1.4 %, and
no UW station is admitted. Gate 1 in the 2008 window will rest on B018 alone.

## 6. Rational-chain equivalence (P4-0 §5)
3 XG + 3 PB stations × 2 days; the rational twin's 20 Hz output against the frozen integer chain
applied to a 100 Hz resample of the same trace; interior samples only (5 s guard).

```
 netsta    chain   min_cc   med_cc  max_diff_over_rms  med_rms_ratio    verdict
PB.B001   frozen 1.000000 1.000000           0.000000       1.000000 equivalent
PB.B003   frozen 1.000000 1.000000           0.000000       1.000000 equivalent
PB.B005   frozen 1.000000 1.000000           0.000000       1.000000 equivalent
XG.BH01 rational 0.999895 0.999924           0.567428       0.995638 equivalent
XG.CL13 rational 0.999921 0.999953           0.447639       0.997837 equivalent
XG.TB01 rational 0.999849 0.999961           0.375706       0.996001 equivalent
```

Verdict: **equivalent** — every XG component has cc ≥ 0.99985 (bar 0.999) and PB is bit-identical
(the rational twin degenerates to the frozen chain at 100 Hz, as designed). One caveat for the
freeze: the *peak* sample difference on XG reaches 0.38–0.57 × RMS and the RMS ratio is 0.996,
i.e. the two anti-alias paths differ by ~0.4 % in amplitude and disagree on isolated samples near
the 8 Hz corner. That is far below the cc bar but it is not zero, so XG template products should
carry `chain=rational` in the freeze and any XG-vs-PB amplitude comparison should not be trusted
below the 1 % level.

## 7. Sealed set
`results/sealed_check.txt`: 5,062 files on disk, 5,062 hashed, 0 unhashed, 0 hashed-but-missing;
**0 cache20hz entries inside 2011-07-20..09-10** and **0 manifest rows for those days written
outside `traces_sealed/`** (the extension queue has the sealed window removed by construction).
**VERDICT: SEALED OK.** The only reads of the sealed tree so far are `os.walk` + sha256 in
`d_seal.py` and the day-count inventory used for the 2011-08 column of the site table — no
waveform sample has been loaded.

## 8. What this changes in the plan
* PLAN §1's station table is superseded by §2 above (sites, delivered-3-C).
* PLAN §10 "XG 8 sites" → **7 sites, 6 usable arrays in the 2010 ETS window** (SL dead, CL fails the probe).
* PLAN §10 "2008-05 XU 44 + ZK 1 + YW 6 + PB 16" → the ZK site is **Z-only in the data**; the 2008
  window's 3-C power is XU (40 admitted) + PB (10) and nothing else.
* Merlin's withheld rule needs the explicit relaxation in §5 before Step 0 can be frozen.
* The 2008 Gate-1 design rests on a single withheld site (B018) — worth telling Merlin before Step 8.