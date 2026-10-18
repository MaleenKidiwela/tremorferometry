#!/usr/bin/env python
"""Royer & Bostock (2014) Washington reproduction -- shared definitions.  Prefix `ab_`.

Everything in this module is quoted from the paper (`notes/1-s2.0-S0012821X13004780-main.pdf`) or
measured from disk.  Nothing is inherited from Bostock et al. (2012): the 6 s / 20 Hz / 51-station
Washington configuration that produced the earlier Puget zeros came from transplanting that paper's
parameters into this region, and every value below supersedes it.

Paper, verbatim
---------------
sec.3   "All data were divided into 24 hour-long segments, band-pass filtered between 1 and 8 Hz,
        resampled to 40 sps"
sec.3.1 NW: "we applied network autocorrelation to 4 tremor episodes (January 17-30, 2007,
        May 01-23, 2008, August 07-24, 2010 and August 04-23, 2011) using the 7 three-component
        stations B001, SQM, BS11, W020, W040, GNW and PL11."
        SW: "we performed network autocorrelation using 2 tremor episodes (January 14-February 01,
        2007, May 01-25, 2008) and employed the 7 three-component stations N050, S030, W070, N060,
        PL11, S040 and N070."
        "Each hour-long segment was divided into 15 s windows lagged by 0.5 s to produce a total of
        7170 individual windows per hour.  Each window was correlated with all other windows in the
        same 1-hour record for an individual (station-component) channel.  The resulting time-series
        from all stations and all components corresponding to the same hour were stacked to create
        'network autocorrelation' records.  When the network autocorrelation exceeded 8 times the
        median absolute deviation ... we registered the corresponding window pairs as LFE
        detections.  This catalogue was then culled to exclude overlapping detections and retain
        only those with high SNR.  After processing all data, we had selected 4915, 3306 and 4267
        initial detections ... from the NW, SW and NVI datasets"
sec.3.2 "combine all channels corresponding to a given detection into a single 'super' trace and
        cross-correlated all such traces against one another to determine maximum correlation
        coefficients and corresponding lags.  The correlation coefficients were used to populate a
        similarity matrix employed within a hierarchical cluster analysis ... resulted in 224, 70
        and 54 initial templates"
sec.3.3 "choosing subsets of 7 to 10 high SNR, three-component stations available within a template
        and scanning through all available tremor dates.  As before, when the summed network
        cross-correlation coefficient exceeded 8 times the median absolute deviation, now for a
        24-hour period, we logged detections.  New templates were formed by stacking waveforms
        (normalized to unit maximum amplitude across three components) for all newly registered
        detections.  After several iterations ... 122, 54 and 47"
sec.5   Hyp2000 then hypoDD (Table S1: IDAT 3, 5 x NITER 8 = 40 iterations, DAMP 50,50,40,40,40);
        "Hypocentres beneath SVI range between 29 and 40 km depth whereas those below Washington
        state fall between 31 and 46 km depth (NW) and 32 and 49 km depth (SW)";
        Hyp2000 mean nominal uncertainty 2.0/1.3 km (NW), 3.7/3.5 km (SW); hypoDD retains 77 % (NW)
        and 72 % (SW).

What is SILENT in the paper (and therefore a swept COLUMN here, never a guessed constant)
-----------------------------------------------------------------------------------------
* the linkage and the cut of the hierarchical cluster analysis (sec.3.2);
* the numeric definition of "cull overlapping detections and retain only those with high SNR";
* the minimum cluster size for an initial template;
* the number of iterations ("several") and the convergence test;
* which 7-10 stations are the "high SNR" subset for a given template.

THE VALIDITY MASK
-----------------
The 20 Hz cache used by every earlier Washington run is built with obspy
`merge(method=1, fill_value=0)` (lfe_picker/scripts/o_cache_month.py:70) and stores no mask, so a
data gap arrives at the correlation as a run of (after band-passing, near-)zero samples.  The
engines defend against it with a RELATIVE guard, `ew > 1e-4 * median(ew)`, which degenerates to
`ew > 0` on any channel that is more than half dead -- 81 station-days of this cache have a
component whose median 6 s window energy is exactly 0.  This cache carries `m`, a per-sample bool
validity mask, and every statistic below is computed over valid samples only.
"""
from __future__ import annotations

import glob
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
ROOT = os.path.dirname(STUDY)
RES = f"{STUDY}/results/ab_royer"
LOGS = f"{STUDY}/logs"
CACHE40 = f"{STUDY}/cache40hz"
TRACES = f"{STUDY}/traces"
sys.path.insert(0, HERE)
for d in (RES, LOGS, CACHE40):
    os.makedirs(d, exist_ok=True)

# ----------------------------------------------------------------- paper parameters ------------
FS = 40.0                  # sec.3  "resampled to 40 sps"
BAND = (1.0, 8.0)          # sec.3  "band-pass filtered between 1 and 8 Hz"
WIN_S = 15.0               # sec.3.1 "15 s windows"
HOP_S = 0.5                # sec.3.1 "lagged by 0.5 s"
NWIN_HOUR = 7170           # sec.3.1 "a total of 7170 individual windows per hour"
NMAD = 8.0                 # sec.3.1 and sec.3.3 "8 times the median absolute deviation"
BLOCK_AUTO = "hour"        # sec.3.1 "on an hour-by-hour basis"
BLOCK_GROW = "day"         # sec.3.3 "now for a 24-hour period"
SUBSET_MIN, SUBSET_MAX = 7, 10      # sec.3.3 "subsets of 7 to 10 high SNR, three-component stations"

SUBARRAY = {
    "NW": dict(
        stas=["B001", "SQM", "BS11", "W020", "W040", "GNW", "PL11"],
        episodes=[("2007-01-17", "2007-01-30"), ("2008-05-01", "2008-05-23"),
                  ("2010-08-07", "2010-08-24"), ("2011-08-04", "2011-08-23")],
        n_det_pub=4915, n_tpl0_pub=224, n_tpl_pub=122,
        depth_pub=(31.0, 46.0), sig_h_pub=2.0, sig_z_pub=1.3, dd_keep_pub=0.77),
    "SW": dict(
        stas=["N050", "S030", "W070", "N060", "PL11", "S040", "N070"],
        episodes=[("2007-01-14", "2007-02-01"), ("2008-05-01", "2008-05-25")],
        n_det_pub=3306, n_tpl0_pub=70, n_tpl_pub=54,
        depth_pub=(32.0, 49.0), sig_h_pub=3.7, sig_z_pub=3.5, dd_keep_pub=0.72),
    # sec.2.2 / 3.1: northern Vancouver Island, POLARIS portable deployment (Audet et al. 2008),
    # autocorrelation "for different combinations of 7 of the 10 three-component stations VI10,
    # VW03, VI11, VW02, VI52, VI53, VI05, VI08, VI06, VI04" on Sept 05-11 2006 and June 13-18
    # 2007.  Literal reading implemented in ab_auto.py: per HOUR, the first 7 of these 10 (in the
    # paper's listed order) that recorded that hour.  Depths "between 33 and 39 km" (sec.5,
    # profiles I, J).
    "NVI": dict(
        stas=["VI10", "VW03", "VI11", "VW02", "VI52", "VI53", "VI05", "VI08", "VI06", "VI04"],
        episodes=[("2006-09-05", "2006-09-11"), ("2007-06-13", "2007-06-18")],
        n_det_pub=4267, n_tpl0_pub=54, n_tpl_pub=47,
        depth_pub=(33.0, 39.0), sig_h_pub=1.8, sig_z_pub=2.16, dd_keep_pub=0.91),
}
N_AUTO = 7          # sec.3.1: every subarray autocorrelates 7 three-component stations

# Network codes and coordinates from results/fdsn_stations_40-50N_2003-2014.txt (FDSN station
# service, 40-50 N, 2003-2014).  XU = the Flexible Array CAFE experiment; PB = PBO borehole;
# UW = PNSN permanent.  XU.N0xx / S0xx / W07x end 2008-12-31, which is exactly why the paper says
# "only data from the 2007 and 2008 tremor episodes were available for the SW subarray".
NET = {"B001": "PB", "SQM": "UW", "BS11": "XU", "W020": "XU", "W040": "XU", "GNW": "UW",
       "PL11": "XU", "N050": "XU", "S030": "XU", "W070": "XU", "N060": "XU", "S040": "XU",
       "N070": "XU"}
COORD = {                       # lat, lon, elev_m
    "B001": (48.043072, -123.131409, 237.0),
    "SQM": (48.073681, -123.048058, 45.0),
    "BS11": (47.956799, -122.925201, 585.7),
    "W020": (47.873100, -122.945602, 349.6),
    "W040": (47.601002, -123.121300, 336.7),
    "GNW": (47.564130, -122.824980, 220.0),
    "PL11": (47.454300, -123.153999, 373.8),
    "N050": (47.224201, -123.395599, 122.8),
    "S030": (47.029099, -123.296501, 112.7),
    "W070": (46.924500, -123.227898, 408.5),
    "N060": (47.150501, -123.164497, 150.7),
    "S040": (47.040100, -123.026604, 170.5),
    "N070": (47.051800, -122.715302, 1.8),
}


# ----------------------------------------------------------------- the 25-station pools --------
# sec.2.1: NW "comprises 25 stations skirting the eastern flanks of the Olympic Peninsula to the west
# of Puget Sound", SW "includes 25 stations running approximately east-west to the southwest of Puget
# Sound"; the stations are those drawn in Fig.1D (red = NW, green = SW, purple = shared NW/SW).
# Read off Fig.1D at 600 dpi (results/ab_royer/POOL_RECONSTRUCTION.md):
#   red    STW B007 OPC B04A LRIV B001 BLN BS11 W020 B943 B013 FACB W030 W040 DOSE C04A HDW GNW
#   purple SQM PL11 W060 SMW N050 N060 N080
#   green  NLWA N015 N020 N030 WISH S015 N040 D03A S020 B017 S030 W070 W080 S040 S050 B018 N070
#          S055 B019 E03A W090
# NW = red + purple = 25 labels (exactly the paper's count).  SW = green + purple = 28 labels for the
# paper's "25": the ambiguity is documented; co-located successive stations (D03A/WISH, B04A/B007/
# LRIV, C04A/DOSE, B013/B943/W030) are separate labels in the figure, and which ones operated in an
# episode is decided by the data, not by us.
# NVI: the 27 stations of Fig.1B (POLARIS PO arms + CNSN CN / C8 stations), read at 600 dpi:
#   MAYB VI30 VI50 VI52 VI53 VI10 VW03 VW02 VI12 VI11 TASB VW01 VI08 WOSB VI07 VI09 VI06 VI54 VI05
#   VI55 VI04 VI03 VI56 VI02 NCRB VI01 VI57   (27, exactly the paper's count; VI31/32/51 exist in
#   the PO inventory but are not drawn).  MAYB, WOSB, NCRB are CN EHZ-only (vertical) stations whose
#   data are at the NRCan FDSN; TASB is C8 (HH?, EarthScope).
POOL = {
    "NVI": ["MAYB", "VI30", "VI50", "VI52", "VI53", "VI10", "VW03", "VW02", "VI12", "VI11",
            "TASB", "VW01", "VI08", "WOSB", "VI07", "VI09", "VI06", "VI54", "VI05", "VI55",
            "VI04", "VI03", "VI56", "VI02", "NCRB", "VI01", "VI57"],
    "NW": ["STW", "B007", "OPC", "B04A", "LRIV", "B001", "BLN", "BS11", "W020", "B943", "B013",
           "FACB", "W030", "W040", "DOSE", "C04A", "HDW", "GNW",
           "SQM", "PL11", "W060", "SMW", "N050", "N060", "N080"],
    "SW": ["NLWA", "N015", "N020", "N030", "WISH", "S015", "N040", "D03A", "S020", "B017", "S030",
           "W070", "W080", "S040", "S050", "B018", "N070", "S055", "B019", "E03A", "W090",
           "SQM", "PL11", "W060", "SMW", "N050", "N060", "N080"],
}
_NVI_COORD = {  # FDSN (EarthScope / NRCan) station service, 2006-2007 epoch
    "VI01": (50.3548, -125.9197, 13.0), "VI02": (50.3313, -126.071, 524.0),
    "VI03": (50.3427, -126.1807, 295.0), "VI04": (50.3133, -126.2479, 245.0),
    "VI05": (50.2863, -126.3714, 512.0), "VI06": (50.2674, -126.4528, 333.0),
    "VI07": (50.2481, -126.543, 393.0), "VI08": (50.2026, -126.6011, 198.0),
    "VI09": (50.1674, -126.5749, 913.0), "VI10": (50.2242, -126.9406, 156.0),
    "VI11": (49.9894, -126.8568, 156.0), "VI12": (50.0537, -127.1186, 9.0),
    "VI30": (50.5587, -127.003, 9.0), "VI50": (50.5489, -126.8241, 33.0),
    "VI52": (50.3898, -126.7654, 635.0), "VI53": (50.3116, -126.7571, 292.0),
    "VI54": (50.1214, -126.4415, 311.0), "VI55": (49.9981, -126.319, 396.0),
    "VI56": (49.9253, -126.21, 319.0), "VI57": (49.8394, -126.1009, 169.0),
    "VW01": (50.1651, -126.6607, 407.0), "VW02": (50.1421, -126.8086, 448.0),
    "VW03": (50.0838, -126.9228, 152.0), "TASB": (49.9266, -126.6513, 18.0),
    "MAYB": (50.4025, -127.1737, 1418.0), "NCRB": (50.4035, -126.0553, 1293.0),
    "WOSB": (50.1608, -126.5704, 961.0)}
_NVI_NET = {s: "PO" for s in _NVI_COORD}
_NVI_NET.update(TASB="C8", MAYB="CN", NCRB="CN", WOSB="CN")

# network EPOCHS (FDSN station service, results/fdsn_stations_40-50N_2003-2014.txt).  C04A was TA
# until 2008-02-18 and XU from 2008-03-05; the others have one network in 2007-2011.
EPOCHS = {"C04A": [("TA", "2005-10-22", "2008-02-18"), ("XU", "2008-03-05", "2010-12-31")]}
_POOL_NET = {"STW": "UW", "B007": "PB", "OPC": "UW", "B04A": "TA", "LRIV": "UW", "BLN": "UW",
             "B943": "PB", "B013": "PB", "FACB": "YW", "W030": "XU", "DOSE": "UW", "C04A": "XU",
             "HDW": "UW", "W060": "XU", "SMW": "UW", "N080": "XU", "NLWA": "US", "N015": "XU",
             "N020": "XU", "N030": "XU", "WISH": "UW", "S015": "XU", "N040": "XU", "D03A": "TA",
             "S020": "XU", "B017": "PB", "W080": "XU", "S050": "XU", "B018": "PB", "S055": "XU",
             "B019": "PB", "E03A": "TA", "W090": "XU"}
_POOL_COORD = {
    "STW": (48.151, -123.671, 0.0), "B007": (48.058, -123.504, 0.0), "OPC": (48.100, -123.413, 0.0),
    "B04A": (48.057, -123.504, 0.0), "LRIV": (48.057, -123.504, 0.0), "BLN": (48.007, -122.973, 0.0),
    "B943": (47.813, -122.911, 0.0), "B013": (47.813, -122.911, 0.0), "FACB": (47.887, -122.646, 0.0),
    "W030": (47.813, -122.911, 0.0), "DOSE": (47.717, -122.972, 0.0), "C04A": (47.717, -122.972, 0.0),
    "HDW": (47.649, -123.053, 0.0), "W060": (47.320, -123.312, 0.0), "SMW": (47.319, -123.344, 0.0),
    "N080": (47.165, -122.552, 0.0), "NLWA": (47.392, -123.869, 0.0), "N015": (47.110, -124.033, 0.0),
    "N020": (47.124, -123.971, 0.0), "N030": (47.190, -123.859, 0.0), "WISH": (47.116, -123.771, 0.0),
    "S015": (47.117, -123.658, 0.0), "N040": (47.194, -123.608, 0.0), "D03A": (47.116, -123.771, 0.0),
    "S020": (47.057, -123.550, 0.0), "B017": (46.996, -123.558, 0.0), "W080": (46.716, -123.285, 0.0),
    "S050": (47.007, -122.871, 0.0), "B018": (46.980, -123.020, 0.0), "S055": (46.839, -122.731, 0.0),
    "B019": (46.653, -123.652, 0.0), "E03A": (46.546, -123.563, 0.0), "W090": (46.597, -123.424, 0.0)}
# borehole 1/2 azimuths from FDSN StationXML (queried 2026-09-25; constant over 2006-2012)
AZ12_EXTRA = {"B013": (31.5, 301.5), "B943": (301.5, 31.5), "B017": (206.0, 296.0),
              "B018": (21.5, 111.5), "B019": (305.0, 35.0), "B007": (36.5, 126.5),
              "B001": (335.5, 65.5)}


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


for _s, _n in _POOL_NET.items():
    NET.setdefault(_s, _n)
for _s, _n in _NVI_NET.items():
    NET.setdefault(_s, _n)
    COORD.setdefault(_s, _NVI_COORD[_s])
for _s, _c in _POOL_COORD.items():
    COORD.setdefault(_s, _c)


def _fdsn_coords():
    """Exact lat/lon/elev for pool stations from the FDSN station table (primary network)."""
    p = f"{STUDY}/results/fdsn_stations_40-50N_2003-2014.txt"
    if not os.path.exists(p):
        return
    F = pd.read_csv(p, sep="|")
    F.columns = [c.strip("# ").strip() for c in F.columns]
    for s_ in _POOL_COORD:
        r = F[(F.Station == s_) & (F.Network == NET[s_])]
        if len(r):
            r = r.iloc[-1]
            COORD[s_] = (float(r.Latitude), float(r.Longitude), float(r.Elevation))


_fdsn_coords()

# ----------------------------------------------------------------- MENDO (2026-09-29) ----------
# Blind Royer-chain discovery at the Mendocino triple junction (Shelly et al. 2025 region).  Additive:
# nothing above changes for NW/SW/NVI.  Stations = Shelly's 17 S1 detection stations + the location-stage
# stations that have raw data in lfe_south/traces (KBS, KKP, KCPB, KCPC, B045, B046 have none).
# Coordinates / networks: lfe_south/results/fdsn_channels_south.csv (latest epoch operating >= 2018).
# Raw miniSEED: lfe_south/traces/<NET.STA>/<yyyy>/<doy>.mseed (same layout as TRACES).  Cache: 40 sps
# rebuilt from raw by ab_cache40.py (the lfe_south 20 sps cache is not used: FS is 40 everywhere here).
# Days: results/ab_royer/tremor/scan_days_MENDO.csv (572 = union of lfe_south
# results/aa_shelly/stage2_days_{build,test}.json, Shelly2025 keys); autocorrelation days
# results/ab_royer/tremor/auto_days_MENDO.csv (written by scripts/mendo_setup.py).
SOUTH = f"{ROOT}/lfe_south"
MENDO_S1 = ["BLAS", "BRIC", "GHOP", "HALS", "JEWT", "PRDS", "SPRL", "RVIT", "HUNT",
            "KMPB", "KMR", "KPR", "KBU", "B049", "B932", "B934", "B933"]
MENDO_LOC = ["KBSB", "KIP", "KBN", "KFP", "KCR", "KSM", "B047", "B935"]
MENDO_BAD = ["BRIC", "RVIT", "KCR", "KSM"]      # dead channels / huge glitches (lfe_south SH_V2 step 1)
RAW_ROOT = {}


def _mendo_setup():
    p = f"{SOUTH}/results/fdsn_channels_south.csv"
    if not os.path.exists(p):
        return None
    F = pd.read_csv(p, dtype=str)
    F = F[F.cha.str[:2].isin(["HH", "BH", "EH"]) & (F.end.str[:4].astype(int) >= 2018)]
    for s_ in MENDO_S1 + MENDO_LOC:
        g = F[F.sta == s_].sort_values("start")
        if not len(g):
            continue
        r = g.iloc[-1]
        NET[s_] = str(r.net)
        COORD[s_] = (float(r.lat), float(r.lon), float(r.elev))
        RAW_ROOT[s_] = f"{SOUTH}/traces"
    s3 = [s_ for s_ in MENDO_S1 if s_ in COORD and s_ not in ("KPR", "KBU")]   # the 3-C S1 stations
    la0 = float(np.mean([COORD[s_][0] for s_ in MENDO_S1 if s_ in COORD]))
    lo0 = float(np.mean([COORD[s_][1] for s_ in MENDO_S1 if s_ in COORD]))
    d = {s_: np.hypot(111.19 * (COORD[s_][0] - la0),
                      111.19 * np.cos(np.radians(la0)) * (COORD[s_][1] - lo0)) for s_ in s3}
    # autocorrelation preference order (NVI rule in ab_auto: per hour the first 7 three-component
    # stations present): nearest the S1-array centroid first, the known-bad BRIC / RVIT last
    # + stations flagged by the self-repeat screen (scripts/mendo_station_screen.py) go after those
    scr = os.environ.get("AB_MENDO_SCREEN", f"{RES}/tremor/auto_screen_MENDO.csv")
    flg = set()
    if os.path.exists(scr):
        _S = pd.read_csv(scr)
        flg = set(_S[_S.flag.astype(bool)].sta)
    auto = sorted(s3, key=lambda s_: (s_ in flg, s_ in MENDO_BAD, d[s_]))
    pool = auto + [s_ for s_ in MENDO_S1 + MENDO_LOC if s_ in COORD and s_ not in auto]
    ad = f"{RES}/tremor/auto_days_MENDO.csv"
    eps = [(x, x) for x in pd.read_csv(ad).date.astype(str)] if os.path.exists(ad) else []
    SUBARRAY["MENDO"] = dict(stas=auto, episodes=eps, n_det_pub=0, n_tpl0_pub=0, n_tpl_pub=27,
                             depth_pub=(22.2, 28.6), sig_h_pub=np.nan, sig_z_pub=np.nan,
                             dd_keep_pub=np.nan)
    POOL["MENDO"] = pool
    return (la0, lo0)


MENDO_CENTROID = _mendo_setup()
MENDO_SET = set(POOL.get("MENDO", []))

# ----------------------------------------------------------------- EXTENSION pool --------------
# results/ab_royer/ext/pool_<SUB>.csv (ab_ext_pool.py): every public station of the allowed networks
# within EXT_RMAX km of the subarray, with its operating epochs.  Loaded here so NET/COORD know
# the extra stations; a code already used by a literal-pool station keeps the literal entry.
EXT_DIR = f"{RES}/ext"
EXT = {}
for _sub in ("NW", "SW", "NVI"):
    _p = f"{EXT_DIR}/pool_{_sub}.csv"
    if os.path.exists(_p):
        _E = pd.read_csv(_p)
        EXT[_sub] = _E
        for _r in _E.itertuples():
            if _r.sta not in NET:
                NET[_r.sta] = _r.net
                COORD[_r.sta] = (float(_r.lat), float(_r.lon), float(_r.elev))


def ext_stations(sub, radius_km):
    E = EXT.get(sub)
    if E is None:
        raise SystemExit(f"extension pool for {sub} not built (ab_ext_pool.py)")
    return sorted(set(E[E.dist_km <= radius_km].sta) | set(POOL[sub]))


def ext_operating(sta, date, sub):
    """True if the extension table says `sta` operated on `date` (literal-pool stations: True)."""
    E = EXT.get(sub)
    if E is None or sta in POOL[sub]:
        return True
    r = E[E.sta == sta]
    return bool(((r.start <= date) & (r.end >= date)).any())


def pool_from_arg(arg, sub):
    """'' -> the literal 25-station pool; 'ext:<R>' -> the extension pool within R km."""
    if not arg:
        return list(POOL[sub])
    if arg.startswith("ext:"):
        return ext_stations(sub, float(arg.split(":")[1]))
    return arg.split(",")


def netsta(sta):
    """Directory name (fixed primary network per station)."""
    return f"{NET[sta]}.{sta}"


def net_for(sta, date):
    """The network code the station had ON `date` (for the FDSN request)."""
    for sub, E in EXT.items():
        if any(sta in P_ for P_ in POOL.values()):
            break
        r = E[(E.sta == sta) & (E.start <= date) & (E.end >= date)]
        if len(r):
            return str(r.iloc[0].net)
    for net, a, b in EPOCHS.get(sta, []):
        if a <= date <= b:
            return net
    return NET[sta]


def days(sub):
    """Every date of every episode of one subarray, as YYYY-MM-DD."""
    out = []
    for a, b in SUBARRAY[sub]["episodes"]:
        out += [x.strftime("%Y-%m-%d") for x in pd.date_range(a, b)]
    return sorted(set(out))


def scan_days(sub):
    """sec.3.3 "all available tremor dates" -- results/ab_royer/tremor/scan_days_<sub>.csv
    (ab_tremor_days.py); falls back to the sec.3.1 days if the list has not been built."""
    p = f"{RES}/tremor/scan_days_{sub}.csv"
    if not os.path.exists(p):
        return days(sub)
    return sorted(pd.read_csv(p).date.astype(str))


def station_days(sub):
    """Every station-day the chain reads: the 7 autocorrelation stations on the sec.3.1 days, and
    the 25-station pool (sec.2.1, Fig.1D) on every scan day (sec.3.3)."""
    out = {(s, d) for d in days(sub) for s in SUBARRAY[sub]["stas"]}
    out |= {(s, d) for d in scan_days(sub) for s in POOL[sub]}
    return sorted(out, key=lambda x: (x[1], x[0]))


def station_table(stas=None):
    """DataFrame indexed by station code with lat/lon/elev -- the shape uu_common expects."""
    stas = list(COORD) if stas is None else list(stas)
    T = pd.DataFrame([dict(sta=s, net=NET[s], lat=COORD[s][0], lon=COORD[s][1],
                           elev=COORD[s][2]) for s in stas]).set_index("sta")
    return T


def epi_km(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = [np.asarray(v, float) for v in (la1, lo1, la2, lo2)]
    return np.sqrt((111.19 * (la2 - la1)) ** 2
                   + (111.19 * np.cos(np.radians(0.5 * (la1 + la2))) * (lo2 - lo1)) ** 2)


# ----------------------------------------------------------------- site-shift null -------------
# WHOLE-CHAIN NULL.  AB_NULL=<k> (k >= 1) makes every stage below read a DIFFERENT DAY for each
# station: each station's day series is circularly shifted WITHIN ITS OWN EPISODE by a
# station-specific whole-day offset.  Members of a site move together, so site-coherent noise, the
# diurnal cycle and the station's own gap pattern are preserved exactly, and only true inter-station
# coherence -- the thing an LFE family is -- is destroyed.  Every stage inherits it from this module,
# which is what makes the null a whole-chain null rather than a detector-only one.
NULL_K = int(os.environ.get("AB_NULL", "0"))


def _episode_days():
    """Null shift groups: the scan days of each (subarray, year), which contain the sec.3.1 days.
    (Before the scan-day list existed these were the sec.3.1 episodes; extra scan days outside
    them would otherwise NOT be shifted under the null.)"""
    out = []
    for sub in SUBARRAY:
        sd = pd.Series(scan_days(sub))
        for _, g in sd.groupby(sd.str[:4]):
            out.append((sub, list(g)))
    return out


_EPI_DAYS = _episode_days()
_MENDO_SITE = {}


def _mendo_site(sta):
    """Site index of a MENDO pool station: single-linkage groups of stations < 3 km apart
    (ab_mute.SITE_KM), numbered in POOL['MENDO'] order."""
    if not _MENDO_SITE:
        P = POOL["MENDO"]
        lab = list(range(len(P)))
        for i, a in enumerate(P):
            for j, b in enumerate(P[:i]):
                if epi_km(*COORD[a][:2], *COORD[b][:2]) < 3.0:
                    old, new = lab[i], lab[j]
                    lab = [new if x == old else x for x in lab]
        u = sorted(set(lab), key=lab.index)
        _MENDO_SITE.update({s_: u.index(x) for s_, x in zip(P, lab)})
    return _MENDO_SITE.get(sta, len(_MENDO_SITE) + sorted(COORD).index(sta))


def null_remap(sta, date, k=None):
    """The date this station actually reads under null realisation k (identity when k == 0)."""
    k = NULL_K if k is None else k
    if not k:
        return date
    for gsub, days_ in _EPI_DAYS:
        if date in days_:
            import zlib
            n = len(days_)
            # AUDIT FIX A12: DISTINCT shifts (a seeded permutation of 1..n-1 per realisation and
            # group) -- the old per-station hash collided, and two stations sharing a shift keep
            # their TRUE mutual coherence inside the null.  Stations are indexed inside the group's
            # own subarray pool (autocorrelation stations first), so the 7 autocorrelation stations
            # never collide and the 25-station pool collides only when the group has < 25 days.
            order = list(SUBARRAY[gsub]["stas"]) + [x for x in POOL[gsub]
                                                    if x not in SUBARRAY[gsub]["stas"]]
            order += [x for x in sorted(COORD) if x not in order]
            rng = np.random.default_rng(zlib.crc32(f"{k}|{gsub}|{days_[0]}".encode()))
            sh = rng.permutation(np.arange(1, n))
            if gsub == "MENDO":
                # MENDO scrambled control: stations < 3 km apart are ONE site and shift together
                # (as armb_null_shift.py / w_null_siteshift.py); sites indexed in pool order
                h = int(sh[_mendo_site(sta) % sh.size])
            else:
                h = int(sh[order.index(sta) % sh.size])
            return days_[(days_.index(date) + h) % n]
    return date


# ----------------------------------------------------------------- 40 Hz cache -----------------
def cache_path(sta, date):
    ts = pd.Timestamp(date)
    return f"{CACHE40}/{netsta(sta)}/{ts.year}/{int(ts.dayofyear):03d}.npz"


def raw_path(sta, date):
    ts = pd.Timestamp(date)
    return f"{RAW_ROOT.get(sta, TRACES)}/{netsta(sta)}/{ts.year}/{int(ts.dayofyear):03d}.mseed"


def load_day(sta, date):
    """(t0, X (3,N) float32 ZNE at 40 sps 1-8 Hz, m (N,) bool validity) or None.

    `m` is False wherever the sample is inside a data gap or inside the band-pass transient that a
    gap edge injects.  Nothing downstream may use a sample whose mask is False.

    Under a site-shift null the WAVEFORM comes from `null_remap(sta, date)` while `t0` is rewritten
    to the REQUESTED date, so every station is still presented on one common absolute clock and only
    the inter-station coherence has been destroyed."""
    src = null_remap(sta, date)
    p = cache_path(sta, src)
    if not os.path.exists(p):
        return None
    z = np.load(p)
    X = z["x"]
    if X.shape[0] != 3:
        return None
    m = z["m"].astype(bool) if "m" in z.files else np.ones(X.shape[1], bool)
    t0 = float(z["t0"])
    X = np.ascontiguousarray(X, dtype=np.float32)
    if sta in MENDO_SET:                       # 2026-09-29 MENDO dead-channel / stuck-span guard
        v = dead_guard(X, m)
        if v is None:
            return None
        X, m = v
    mdir = os.environ.get("AB_MUTE_DIR", "")
    if mdir:                                   # 2026-09-26 burst muting (scripts/ab_mute.py)
        X, m = apply_mute(mdir, sta, src, t0, X, m)
    if src != date:
        t0 += (pd.Timestamp(date) - pd.Timestamp(src)).total_seconds()
    return t0, X, m


DEAD_WIN_S = 2.0
DEAD_ABS = 1e-6          # counts^2 per 2 s window: a station whose loudest channel median is below
                         # this carries only float residue (real integer-count data read >= ~1)


def dead_guard(X, m):
    """MENDO port of lfe_south SH_V2 step 1 (sh_v2_common.prep_day) to the 40 sps cache, applied in
    load_day BEFORE the burst mute so every stage (energy for the mute, autocorrelation, clustering,
    growth, restacks) sees the same cleaned day.
      * 2 s window energies on valid samples; a channel whose PLAIN daily median is <= GUARD_STA_REL
        (1e-4) x the station's loudest channel median is DEAD for the day (BK.BRIC / BK.RVIT HHE,
        NC.KCR / KSM EHZ residue 1e-17..1e-13) -> set to exactly 0 (every stage then drops it);
      * all channels dead, or loudest median < DEAD_ABS -> station absent that day (None);
      * stuck / flat spans on a LIVE channel: 2 s windows with energy < GUARD_REL (1e-3) x that
        channel's median are marked invalid in m (padded one window each side), exactly like gaps.
    The station-level network floor (GUARD_NET_REL) stays in day_refs (growth / clustering)."""
    n = int(round(DEAD_WIN_S * FS))
    N = X.shape[1] // n * n
    if N == 0:
        return None
    x2 = (X[:, :N].astype(np.float64) ** 2).reshape(3, -1, n).sum(-1)
    mv = np.asarray(m[:N], bool).reshape(-1, n).all(-1)
    if mv.sum() < 10:
        return None
    pm = np.array([np.median(x2[c, mv & (x2[c] > 0)]) if (mv & (x2[c] > 0)).any() else 0.0
                   for c in range(3)])
    top = float(pm.max())
    if top < DEAD_ABS:
        return None
    dead = pm <= GUARD_STA_REL * top
    X = X.copy() if dead.any() else X
    X[dead] = 0.0
    live = ~dead
    low = (x2[live] < GUARD_REL * pm[live][:, None]).any(0) & mv
    if low.any():
        low = low | np.r_[low[1:], False] | np.r_[False, low[:-1]]
        bad = np.repeat(low, n)
        m = np.asarray(m, bool).copy()
        m[:N] &= ~bad
        X = X.copy()
        X[:, :N][:, bad] = 0.0
    return X, m


def apply_mute(mdir, sta, date, t0, X, m):
    """2026-09-26 burst muting + line notching (scripts/ab_mute.py).  The mute mask is ANDed into the
    validity mask: muted samples become exactly like gap samples (m False, x 0 as a gap is; never
    zeroed without masking).  Then the persistent narrow spectral lines listed for this station-day
    are notched (zero-phase, 0.06 Hz wide) over each valid run >= 60 s.  Mask files are on a 2 s grid
    from midnight; no mask file = nothing muted / notched."""
    ts = pd.Timestamp(date)
    p = f"{mdir}/mask/{netsta(sta)}/{ts.strftime('%Y-%m-%d')}.npz"
    if not os.path.exists(p):
        return X, m
    z = np.load(p)
    mu = z["mute"].astype(bool)
    lines = z["lines"] if "lines" in z.files else np.zeros((0, 3))
    if mu.any():
        nw = int(round(float(z["win"]) * FS))
        i0 = int(round((ts.timestamp() - t0) * FS))           # sample index of midnight
        k = np.flatnonzero(mu)
        a = np.clip(i0 + k * nw, 0, X.shape[1])
        b = np.clip(i0 + (k + 1) * nw, 0, X.shape[1])
        d = np.zeros(X.shape[1] + 1, np.int32)
        np.add.at(d, a, 1)
        np.add.at(d, b, -1)
        ms = np.cumsum(d[:-1]) > 0
        if ms.any():
            m = m & ~ms
            X = X.copy()
            X[:, ms] = 0.0
    if len(lines):
        import ab_mute
        X = ab_mute.notch(X, m, lines)
    return X, m


class DayCache:
    """{sta: (t0, X, m)} for ONE date.  Built once per date and reused."""

    def __init__(self, stas):
        self.stas = list(stas)
        self.date = None
        self.D = {}

    def day(self, date):
        if date != self.date:
            self.D = {}
            for s in self.stas:
                v = load_day(s, date)
                if v is not None:
                    self.D[s] = v
            self.date = date
        return self.D

    def hour(self, date, hh):
        """{sta: (X (3,nh), m (nh,))} for one whole hour on a common absolute time base."""
        D = self.day(date)
        t_h = pd.Timestamp(f"{date} {hh:02d}:00:00").timestamp()
        nh = int(round(3600 * FS))
        out = {}
        for s, (t0, X, m) in D.items():
            i0 = int(round((t_h - t0) * FS))
            if i0 < 0 or i0 + nh > X.shape[1]:
                continue
            out[s] = (X[:, i0:i0 + nh], m[i0:i0 + nh])
        return out

    def cut(self, sta, t_abs, len_s, date=None):
        """(3, n) float64 segment, or None if it is outside the day OR touches an invalid sample."""
        date = (pd.Timestamp(float(t_abs), unit="s").strftime("%Y-%m-%d") if date is None
                else date)
        v = self.day(date).get(sta)
        if v is None:
            return None
        t0, X, m = v
        n = int(round(len_s * FS))
        i0 = int(round((t_abs - t0) * FS))
        if i0 < 0 or i0 + n > X.shape[1]:
            return None
        if not m[i0:i0 + n].all():
            return None
        seg = X[:, i0:i0 + n].astype(np.float64)
        return seg if np.isfinite(seg).all() else None


# ----------------------------------------------------------------- window statistics -----------
def window_stats(x, m, n, hop):
    """Per-window energy and validity on the 0.5 s row grid, mask-aware.

    x (C, N) float, m (N,) bool, n samples per window, hop samples between window starts.
    Returns (ew (C, NW) float64 window energies with invalid samples contributing 0,
             vw (NW,) bool  True only where EVERY sample of the window is valid).
    No relative-median guard is used anywhere: validity is the mask, not an energy ratio."""
    x = np.asarray(x, float)
    C, N = x.shape
    NW = (N - n) // hop + 1
    mm = np.asarray(m, bool)
    xs = np.where(mm[None, :], x, 0.0)
    c = np.concatenate([np.zeros((C, 1)), np.cumsum(xs * xs, axis=1)], 1)
    starts = np.arange(NW) * hop
    ew = c[:, starts + n] - c[:, starts]
    cv = np.concatenate([[0], np.cumsum(mm.astype(np.int64))])
    vw = (cv[starts + n] - cv[starts]) == n
    return ew, vw


# ----------------------------------------------------------------- 2026-09-25 gap guard --------
# DEFECT (NW growth v4b): `e > 0` let near-zero-energy channel-windows into the per-channel cc
# denominator -> network cc up to 106 on 7 days.  Two causes seen on 2008-07-05:
#   (a) the day-long cumsum in window_stats loses all precision after a huge glitch (FACB N:
#       |x| 7e6 -> 79 % of later windows read e = 0 or 3e-7 instead of ~1e3);
#   (b) dead channels carrying only tiny non-zero residue (SMW Z: e ~ 1e-36..1e-5).
# FIX: exact (block-local) window energies + a RELATIVE guard.  A channel-window is valid only if
#   mask-valid, AND <= GUARD_MAX_FLAT of its samples are exactly zero / constant (repeat of the
#   previous sample), AND e > GUARD_REL x ref_c.  ref_c (channel_refs) = median energy of the
#   channel's TRULY non-zero windows that day (not flat, and above GUARD_STA_REL x the station's
#   loudest component level -- so a mostly-dead channel's median is not dragged to 0 or to its
#   residue); a channel with < GUARD_MIN_LIVE_S of such windows is dead; a station whose level is
#   < GUARD_NET_REL x the day's median station level is dead (day_refs; unit-free).
GUARD_REL = 1e-3
GUARD_STA_REL = 1e-4
GUARD_MAX_FLAT = 0.01
GUARD_NET_REL = 1e-6        # station floor vs the day's median station level (unit-free)
GUARD_MIN_LIVE_S = 600.0     # a channel needs >= 10 min of live windows that day to define ref_c


def _box_energy(x2, n, starts):
    """Exact sliding sums of x2 (C, N) over [s, s+n) for s in starts, without a day-long cumsum:
    block-local suffix + prefix sums (no cancellation; error relative to one block only)."""
    C, N = x2.shape
    K = -(-N // n) + 1
    pad = np.zeros((C, K * n))
    pad[:, :N] = x2
    B = pad.reshape(C, K, n)
    pre = np.concatenate([np.zeros((C, K, 1)), np.cumsum(B, axis=2)], 2)          # pre[k, r]
    suf = np.concatenate([np.cumsum(B[:, :, ::-1], axis=2)[:, :, ::-1],
                          np.zeros((C, K, 1))], 2)                                 # suf[k, r]
    k = starts // n
    r = starts % n
    return suf[:, k, r] + pre[:, k + 1, r]


def _flat_samples(x, m):
    """(C, N) bool: sample is exactly zero, repeats the previous sample, or is mask-invalid."""
    flat = (x == 0)
    flat[:, 1:] |= (x[:, 1:] == x[:, :-1])
    flat |= ~np.asarray(m, bool)[None, :]
    return flat


def channel_refs(x, m, n, rel_sta=GUARD_STA_REL, max_flat=GUARD_MAX_FLAT, sub=20):
    """(C,) per-channel reference window energy for one station-day (0 = dead channel), computed
    on a `sub`-sample grid of windows (medians only).

    "Truly non-zero" windows: live (flat fraction <= max_flat) AND e > rel_sta x the station's
    loudest component median -- a dead channel carrying filter residue (e ~ 1e-17 on FACB N
    2008-07-05) is not flat, so the median must also exclude sub-residue windows or it collapses
    to the residue level.  A component with no such windows (>= GUARD_MIN_LIVE_S) is dead."""
    x = np.asarray(x, np.float64)
    C, N = x.shape
    if N < n:
        return np.zeros(C)
    starts = np.arange(0, N - n + 1, sub)
    e = _box_energy(x * x, n, starts)
    fl = _flat_samples(x, m).astype(np.int32)
    cf = np.concatenate([np.zeros((C, 1), np.int64), np.cumsum(fl, 1)], 1)
    live = ((cf[:, starts + n] - cf[:, starts]) / n) <= max_flat
    min_live = max(int(GUARD_MIN_LIVE_S * FS / sub), 1)
    med = np.array([float(np.median(e[c, live[c]])) if live[c].sum() >= min_live else 0.0
                    for c in range(C)])
    top = float(med.max()) if med.size else 0.0
    ref = np.zeros(C)
    if top <= 0:
        return ref
    for c in range(C):
        good = live[c] & (e[c] > rel_sta * top)
        if good.sum() >= min_live:
            ref[c] = float(np.median(e[c, good]))
    return ref


def day_refs(D, n, net_rel=GUARD_NET_REL):
    """{sta: channel_refs} for one date's {sta: (t0, X, m)}, plus the NETWORK floor: a station
    whose loudest channel ref is < net_rel x the median over stations is dead (all-channel
    residue, e.g. SMW 2008-07-05 at ~1e-20 counts), whatever its internal ratios."""
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max(min(len(D), 16), 1)) as ex:
        R = dict(zip(D, ex.map(lambda s: channel_refs(D[s][1], D[s][2], n), D)))
    tops = np.array([r.max() for r in R.values() if r.max() > 0])
    if tops.size:
        floor = net_rel * float(np.median(tops))
        for s in R:
            if R[s].max() < floor:
                R[s] = np.zeros_like(R[s])
    return R


def window_stats_guarded(x, m, n, hop, rel=GUARD_REL, max_flat=GUARD_MAX_FLAT, ref=None):
    """Like window_stats, but (ew exact energies (C, NW), okc (C, NW) bool per-channel validity
    under the relative gap guard above).  okc already includes the mask test."""
    x = np.asarray(x, np.float64)
    mm = np.asarray(m, bool)
    C, N = x.shape
    NW = (N - n) // hop + 1
    starts = np.arange(NW) * hop
    xs = np.where(mm[None, :], x, 0.0)
    ew = _box_energy(xs * xs, n, starts)
    cv = np.concatenate([[0], np.cumsum(mm.astype(np.int64))])
    vw = (cv[starts + n] - cv[starts]) == n
    fl = _flat_samples(x, mm).astype(np.int32)
    cf = np.concatenate([np.zeros((C, 1), np.int64), np.cumsum(fl, 1)], 1)
    ff = (cf[:, starts + n] - cf[:, starts]) / n
    if ref is None:
        ref = channel_refs(x, mm, n, max_flat=max_flat)
    okc = vw[None, :] & (ff <= max_flat) & (ref[:, None] > 0) & (ew > rel * ref[:, None])
    return ew, okc


def segment_channel_ok(seg, ref, n_ref, rel=GUARD_REL, max_flat=GUARD_MAX_FLAT):
    """(C,) bool: the gap guard for one already-cut, mask-valid segment (C, n), given the
    station-day channel refs (computed for n_ref-sample windows; the segment energy is scaled to
    that length).  Used by the stacking / super-trace code."""
    seg = np.asarray(seg, np.float64)
    e = (seg ** 2).sum(1) * (float(n_ref) / seg.shape[1])
    fl = _flat_samples(seg, np.ones(seg.shape[1], bool)).mean(1)
    return np.isfinite(seg).all(1) & (fl <= max_flat) & (ref > 0) & (e > rel * ref)


class RefCache:
    """{(sta, date): channel refs} for a DayCache (window length n samples), built per date with
    the network floor of day_refs."""

    def __init__(self, dc, n):
        self.dc, self.n, self.R = dc, int(n), {}

    def get(self, sta, date):
        if date not in self.R:
            self.R = {date: day_refs(self.dc.day(date), self.n)}   # one date resident
        return self.R[date].get(sta, np.zeros(3))


def valid_runs(m, min_len):
    """[(i0, i1), ...] maximal runs of True in m with i1 - i0 >= min_len."""
    mm = np.asarray(m, bool).astype(np.int8)
    d = np.diff(np.concatenate([[0], mm, [0]]))
    a = np.flatnonzero(d == 1)
    b = np.flatnonzero(d == -1)
    return [(int(i), int(j)) for i, j in zip(a, b) if j - i >= min_len]


def mad(v):
    v = np.asarray(v, float)
    med = float(np.median(v))
    return med, float(np.median(np.abs(v - med)))


def snr_window(x, m=None):
    """Royer sec.3.1's cull needs an SNR per detection; the paper is SILENT on its definition.

    Reported (never thresholded inside the detector): the ratio of the peak 1 s rms inside the
    window to the median 1 s rms of the window, summed over 3 components -- a scale-free,
    band-limited impulsiveness that does not depend on absolute gain."""
    x = np.asarray(x, float)
    n1 = max(int(round(1.0 * FS)), 1)
    e = (x ** 2).sum(0)
    if m is not None:
        e = np.where(np.asarray(m, bool), e, np.nan)
    k = np.convolve(np.nan_to_num(e), np.ones(n1) / n1, mode="valid")
    if k.size == 0 or np.median(k) <= 0:
        return np.nan
    return float(np.sqrt(k.max() / np.median(k)))


BAND_PREF = ["BH", "HH", "EH", "SH"]


def pick_band(st):
    """AUDIT FIX A8: keep ONE three-component (location, band) of an obspy Stream.

    FDSN answers to `BH?,HH?,EH?,SH?` hold two bands on instrument-transition days (UW.SQM
    2011-04..07 BH+EH, XU.BS11 2008 BH/HH) and the cache refuses mixed rates.  Preference: a
    (loc, band) with a Z and two horizontals, band order BH > HH > EH > SH, empty location first.
    Returns (stream, tag); if no 3-C set exists the stream is returned unchanged (the cache then
    records `not three components`)."""
    groups = {}
    for tr in st:
        ch = tr.stats.channel
        groups.setdefault((tr.stats.location, ch[:2]), set()).add(ch[-1])
    ok = [(k, v) for k, v in groups.items()
          if "Z" in v and len(v & set("NE12")) >= 2]
    if not ok:
        return st, "no3c"
    ok.sort(key=lambda kv: (BAND_PREF.index(kv[0][1]) if kv[0][1] in BAND_PREF else 99,
                            kv[0][0] != ""))
    loc, band = ok[0][0]
    return st.select(location=loc, channel=band + "?"), f"{loc}.{band}"


def az12(sta, date, chans):
    """(az1, az2) of the two horizontal channels `chans` on `date`: the fixed tables first, else the
    FDSN StationXML (cached per station under results/ab_royer/az12/)."""
    if sta in AZ12_EXTRA:
        return AZ12_EXTRA[sta]
    import json
    d = f"{RES}/az12"
    os.makedirs(d, exist_ok=True)
    p = f"{d}/{net_for(sta, date)}.{sta}.json"
    if not os.path.exists(p):
        from obspy.clients.fdsn import Client
        inv = Client("EARTHSCOPE", timeout=120).get_stations(
            network=net_for(sta, date), station=sta, level="channel")
        rows = [dict(ch=c.code, loc=c.location_code, az=float(c.azimuth or 0.0),
                     start=str(c.start_date)[:10], end=str(c.end_date)[:10])
                for n in inv for st_ in n for c in st_]
        tmp = f"{p}.tmp{os.getpid()}"
        json.dump(rows, open(tmp, "w"))
        os.replace(tmp, p)
    rows = json.load(open(p))
    out = []
    for ch in chans:
        hit = [r for r in rows if r["ch"] == ch and r["start"] <= date
               and (r["end"] in ("None", "") or r["end"] >= date)]
        if not hit:
            return None
        out.append(hit[0]["az"])
    return tuple(out)


def torch_dev():
    import torch
    return "cuda" if torch.cuda.is_available() else "cpu"
