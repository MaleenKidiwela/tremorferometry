# Royer & Bostock (2014): the 25-station arrays and the scan dates, reconstructed

2026-09-25. Sources: sec.2.1 of the paper, Fig.1D (rendered at 600 dpi from
`notes/1-s2.0-S0012821X13004780-main.pdf`), `results/fdsn_stations_40-50N_2003-2014.txt`, and the
PNSN tremor pages. Code: `scripts/ab_common.py` (`POOL`, `EPOCHS`, `AZ12_EXTRA`) and
`scripts/ab_tremor_days.py`.

## 1. Station pools (sec.2.1 and Fig.1D)

The legend of Fig.1 colours the stations: NW = red, SW = green, shared NW/SW = purple.

| colour | stations read off Fig.1D |
|---|---|
| red (NW) | STW B007 OPC B04A LRIV B001 BLN BS11 W020 B943 B013 FACB W030 W040 DOSE C04A HDW GNW (18) |
| purple (shared) | SQM PL11 W060 SMW N050 N060 N080 (7) |
| green (SW) | NLWA N015 N020 N030 WISH S015 N040 D03A S020 B017 S030 W070 W080 S040 S050 B018 N070 S055 B019 E03A W090 (21) |

* **NW pool** = red + purple = **25 labels**. This matches the paper's "25 stations".
* **SW pool** = green + purple = **28 labels**, against the paper's "25". This is an unresolved
  ambiguity. Several labels are successive stations at the same site, and those sites
  operated in different years:
  * D03A (TA, until 2008-02) / WISH (UW), same coordinates;
  * B04A (TA, until 2008-02) / B007 (PB) / LRIV (UW, from 2009-08);
  * C04A (TA to 2008-02, then XU) / DOSE (UW, from 2010-07);
  * B013 / B943 (PB) / W030 (XU).

  No label is removed. Whether a station operated on a given day is decided by the FDSN archive:
  an absent station-day returns `nodata`, which is recorded in `fetch_manifest.csv`.
* Networks and epochs come from the FDSN station table. C04A's network changes with the date:
  TA until 2008-02-18, XU from 2008-03-05.
* For B001, B007, B013, B943, B017, B018 and B019, the borehole 1/2 azimuths come from FDSN
  StationXML (queried 2026-09-25).
* Vertical-only station-days (UW short-period EHZ stations such as STW, BLN, HDW and SMW, and SQM
  after 2011-04) are cached with zero horizontals. The paper chooses "three-component stations"
  for the growth subset (sec.3.3), so these are never picked for a growth template. They do
  contribute their vertical channel to the sec.3.2 cluster stacks ("all channels available") and
  P picks to location.
* **Growth** (sec.3.3): each template uses the **7-10 highest-SNR three-component stations** of
  its subarray pool. SNR is the envelope peak over the median of the template stack. The stations
  are ranked, with no SNR cut.

## 2. Scan dates (sec.3.3 "all available tremor dates")

The paper gives no dates for the growth scan. The rule, applied in order:

1. Every sec.3.1 autocorrelation day is a scan day.
2. **From 2009-08-06**, when a public machine-readable tremor location catalogue exists (PNSN/Wech;
   `catalogs/pnsn_tremor_cascadia_full.csv` plus the 2009 API pull in
   `results/ab_royer/tremor/pnsn_tremor_2009.csv`): a UTC day is a tremor day when it holds
   **at least 24 PNSN epicentres (2.5-min windows, i.e. at least 1 hour of tremor) within 25 km of
   any station of the subarray's pool**.
   * NW is scanned from 2009-08-06 to 2011-12-31.
   * SW has no data after 2008 (sec.2.1).
   * How the day count depends on the threshold is in `tremor/rule_sensitivity.csv`: NW 2009-08..2011 has
     155 / 113 / 86 / 68 days at >=1 / 12 / 24 / 48 epicentres.
3. **Before 2009-08-06** no machine-readable location catalogue is public. The PNSN archive
   (`assets.pnsn.org/oldets`) holds only daily map images. The scan days are therefore the PNSN
   published northern-Washington ETS windows (Wech's WECC; `oldets/compare_maps/Density.html`):
   * winter 2007: Jan 12 - Feb 2;
   * spring 2008: Apr 30 - Jun 2;
   * spring 2009: Apr 28 - Jun 3 (NW only; the PNSN 2009 ETS page reports tremor from about
     Apr 27 in southern Puget Sound, passing under the northern Olympic coast about May 10).

   These windows cover the whole N-Washington / S-Vancouver-Island ETS, so they cannot be
   restricted to the subarray region. They over-include days, which costs only compute.

**Result:** NW has **180 scan days**: 2007: 22, 2008: 34, 2009: 47, 2010: 32, 2011: 45.
SW has **56 scan days**: 2007: 22, 2008: 34.

**Where the paper cannot be followed exactly:**
* Inter-ETS tremor days before 2009-08 cannot be identified because no public catalogue exists.
  "All available tremor dates" in 2007-2009 is therefore limited to the ETS windows.
* The exact dates Royer scanned are not recoverable.
