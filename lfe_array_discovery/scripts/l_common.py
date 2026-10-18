#!/usr/bin/env python
"""Washington N-Olympic 2010 bed adapter for agent K's LITERAL Bostock chain.  Prefix `l_`.

NOTHING in `scripts/k_*.py` is edited.  Every `k_` module reads its bed facts through
`k_common` (`K.load_day`, `K.cache_index`, `K.station_table`) at CALL time, so this module
rebinds those IN PROCESS -- exactly the way `x_common.activate()` rebinds `w_common` -- and the
k_ stages then run unchanged on the WA data.  The recipe constants (6 s windows, 1-8 Hz, plain
per-component mean, 8 x MAD, 4 s decluster, Bungum-Husebye realignment, 12 s templates, no merge)
live in k_auto / k_group / k_grow and are NOT touched here.

What the WA bed genuinely needs, and nothing else:

  1. A STATION TABLE.  `w_common.station_table()` holds only the 45 north-era stations, which is
     what aborted the 2010 permanent bed with `KeyError: 'PFB'` (ledger sec.9).  Built here from
     `results/channels_2006_2012.csv` for the probe-admitted WA stations, net-qualified.
  2. BOREHOLE AZIMUTHS.  `k_common.load_day` rejects any cache file whose channel order is not
     Z,N,E unless the station is in `s_common.AZ12`; the PB boreholes are cached (EHZ, EH1, EH2)
     and would be silently dropped.  AZ12 is extended from the StationXML azimuths in
     `results/channels_2006_2012.csv`, so they are properly ROTATED, not merely accepted.
  3. A NET-QUALIFIED CACHE INDEX.  `k_common.cache_index()` keys on the bare station code, and the
     WA cache contains C8.PA02..PA12 beside XG.PA02..PA12.  The index is pre-filled here from the
     explicit (net, sta) list so no code can be resolved to the wrong network, and only this
     study's cache root is used.
  4. The STATION-HOUR HEALTH VETO (`results/wa_station_hour_veto_2010.csv`) as a per-hour mask,
     plus the station-level exclusion of any station > 30 % vetoed hours (XG.DR09 alone).
  5. The whole-chain SITE-SHIFT NULL: a station-specific whole-day offset applied on read, with
     the nominal day's `t0` restamped so every downstream stage is byte-identical.

Environment
  L_NULL_SHIFT=1   activate the day-shift null (offsets from L_NULL_SEED, default 11)
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
RES = f"{STUDY}/results"
CACHE = f"{STUDY}/cache20hz"
ROOT = "/home/jovyan/tremorferometry"
PICKER = f"{ROOT}/lfe_picker"
PNSN = f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv"
SLAB2 = f"{ROOT}/data/cas_slab2_input_04-18.csv"
CHESTLER = (f"{ROOT}/codex-LFE/results/olympic_reference_discovery_v1/"
            "derived_coordinate_rows_not_a_catalog.csv")
sys.path.insert(0, HERE)
sys.path.insert(0, f"{PICKER}/scripts")

# ---- the bed (results/WA_V2_INPUTS.md sec.1.1, results/P_STEP2_SITES.md) -----------------------
YEAR = 2010
DAY0, DAY1 = "2010-07-20", "2010-09-10"
SUBREGION = "N-Olympic"
LAT0, LAT1, LON0, LON1 = 47.8, 48.3, -123.7, -122.7
PROBE = f"{RES}/site_probe_2010.csv"
VETO = f"{RES}/wa_station_hour_veto_2010.csv"
QUIET = f"{RES}/wa_quiet_hours_2010.csv"
HOURS = f"{RES}/l_wa2010/candidate_hours.csv"
VETO_STATION_MAX = 0.30          # WA_V2_INPUTS sec.3: drop a station > 30 % vetoed hours
OUT = os.environ.get("L_OUT", f"{RES}/l_wa2010")

_S = {}


def log(*a):
    import time as _t
    print(f"[{_t.strftime('%H:%M:%S')}]", *a, flush=True)


# ------------------------------------------------------------------ the admitted station set ---
def admitted():
    """DataFrame(net, sta, netsta, veto_frac) of the probe-admitted WA 2010 stations that survive
    the station-level health bar.  No station is withheld (the north's one real lever is station
    count: ledger sec.7/8, 11 -> 15 stations raised true pairs 8.4 % -> 12.9 %)."""
    if "A" not in _S:
        A = pd.read_csv(PROBE)
        A = A[A.admit].copy()
        V = pd.read_csv(VETO, usecols=["netsta", "veto"])
        f = V.groupby("netsta").veto.mean()
        A["veto_frac"] = [float(f.get(n, 0.0)) for n in A.netsta]
        A["health_drop"] = A.veto_frac > VETO_STATION_MAX
        _S["A"] = A.sort_values("netsta").reset_index(drop=True)
    return _S["A"]


def stations():
    A = admitted()
    return list(A[~A.health_drop].sta)


def netsta():
    A = admitted()
    return {r.sta: r.netsta for r in A[~A.health_drop].itertuples()}


# ------------------------------------------------------------------ station table ---------------
def station_table():
    """One row per admitted station: lat, lon, net, netsta, fs_native.  Indexed by bare `sta`,
    because every k_/s_ helper addresses stations by bare code -- but resolved through the
    explicit (net, sta) pair, so C8.PA02 can never stand in for XG.PA02."""
    if "ST" not in _S:
        C = pd.read_csv(f"{RES}/channels_2006_2012.csv").dropna(subset=["lat", "lon"])
        ns = netsta()
        C["netsta"] = C.net + "." + C.sta
        C = C[C.netsta.isin(set(ns.values()))]
        g = (C.groupby(["net", "sta"], as_index=False)
               .agg(lat=("lat", "median"), lon=("lon", "median"), elev_m=("elev", "median"),
                    fs_native=("fs", "max"), band=("band", "first")))
        g["netsta"] = g.net + "." + g.sta
        miss = sorted(set(ns) - set(g.sta))
        if miss:
            raise SystemExit(f"no channel metadata for {miss}")
        _S["ST"] = g.set_index("sta")
    return _S["ST"]


def az12_wa():
    """{sta: (az_N-equivalent, az_E-equivalent)} for every admitted station whose cache is not
    already Z,N,E -- i.e. the PB boreholes cached as (EHZ, EH1, EH2)."""
    C = pd.read_csv(f"{RES}/channels_2006_2012.csv")
    ns = netsta()
    C["netsta"] = C.net + "." + C.sta
    C = C[C.netsta.isin(set(ns.values()))]
    out = {}
    for s in ns:
        c = C[(C.sta == s) & (C.chan.str.len() == 3)]
        a1 = c[c.chan.str.endswith("1")]
        a2 = c[c.chan.str.endswith("2")]
        ah = a1[a1.chan.str[:2].isin(["EH", "BH", "HH", "SH"])]
        bh = a2[a2.chan.str[:2].isin(["EH", "BH", "HH", "SH"])]
        if len(ah) and len(bh):
            out[s] = (float(ah.azimuth.iloc[0]), float(bh.azimuth.iloc[0]))
    return out


# ------------------------------------------------------------------ cache index -----------------
def cache_index():
    """{(sta, year, doy): path} built from the EXPLICIT (net, sta) list and this study's cache
    root only.  Replaces `k_common.cache_index`, whose bare-code glob is ambiguous in WA."""
    idx = {}
    for s, ns in netsta().items():
        for p in glob.glob(f"{CACHE}/{ns}/*/*.npz"):
            parts = p.split(os.sep)
            idx[(s, int(parts[-2]), int(parts[-1][:3]))] = p
    return idx


# ------------------------------------------------------------------ health veto -----------------
def health_veto():
    """{(date, hour): set of vetoed bare station codes} over the admitted set."""
    if "V" not in _S:
        V = pd.read_csv(VETO, usecols=["netsta", "date", "hour", "veto"])
        ns = set(netsta().values())
        V = V[V.veto & V.netsta.isin(ns)]
        d = {}
        for r in V.itertuples():
            d.setdefault((r.date, int(r.hour)), set()).add(r.netsta.split(".")[1])
        _S["V"] = d
    return _S["V"]


# ------------------------------------------------------------------ hours -----------------------
def bed_days():
    return [d.strftime("%Y-%m-%d") for d in pd.date_range(DAY0, DAY1)]


def candidate_hours():
    """PNSN tremor hours: >= 5 epicentres within 50 km of an active cell of the sub-region
    (w_common.AC_HOUR_MIN_EPI / AC_HOUR_EPI_KM, reproduced through the frozen WA adapter)."""
    if os.path.exists(HOURS):
        return pd.read_csv(HOURS)
    import x_common as X
    W = X.activate()
    H = W.candidate_hours()
    os.makedirs(os.path.dirname(HOURS), exist_ok=True)
    H.to_csv(HOURS, index=False)
    return H


def quiet_hours():
    """The 40 hour-matched, earthquake-screened quiet hours of results/WA_V2_INPUTS.md sec.2.1."""
    return pd.read_csv(QUIET)[["date", "hour"]]


def tremor_hour_set():
    H = candidate_hours()
    return {(r.date, int(r.hour)) for r in H.itertuples()}


def pnsn(pad_days=3):
    """PNSN tremor epicentres inside the bed window (+- pad)."""
    T = pd.read_csv(PNSN, usecols=["time", "lat", "lon"])
    T["time"] = pd.to_datetime(T.time, errors="coerce").astype("datetime64[ns]")
    T = T.dropna(subset=["time"])
    t0 = pd.Timestamp(DAY0) - pd.Timedelta(days=pad_days)
    t1 = pd.Timestamp(DAY1) + pd.Timedelta(days=pad_days + 1)
    T = T[(T.time >= t0) & (T.time < t1)].reset_index(drop=True)
    T["t"] = T.time.astype("int64") / 1e9
    return T


def active_cells():
    import x_common as X
    return X.active_cells(sub=True)


def station_corrections():
    """Intraslab-earthquake S station terms (results/x_station_corrections.csv, x_corrections.py).
    A station without a usable term gets 0.0 and is FLAGGED -- never silently imputed."""
    p = f"{RES}/x_station_corrections.csv"
    corr, flagged = {}, []
    if os.path.exists(p):
        C = pd.read_csv(p).set_index("sta")
        for s in stations():
            if s in C.index and bool(C.loc[s, "usable"]):
                v = float(C.loc[s, "corr_S"])
                corr[s] = v if np.isfinite(v) else 0.0
            else:
                corr[s] = 0.0
                flagged.append(s)
    else:
        corr = {s: 0.0 for s in stations()}
        flagged = list(stations())
    return corr, flagged


# ------------------------------------------------------------------ the null --------------------
def null_offsets(stas, seed=None):
    """Station-specific WHOLE-DAY offsets (never 0) over the bed window."""
    seed = int(os.environ.get("L_NULL_SEED", "11")) if seed is None else seed
    days = bed_days()
    rng = np.random.default_rng(seed)
    return {s: int(rng.integers(1, len(days))) for s in sorted(stas)}


# ------------------------------------------------------------------ activation ------------------
def activate(null=None):
    """Rebind k_common onto the WA bed.  Idempotent.  Returns the k_common module."""
    import k_common as K
    if _S.get("ACT"):
        return K
    import s_common as S
    ST = station_table()
    K._C["ST"] = ST
    S.station_table = lambda: ST
    S.AZ12.update(az12_wa())
    for s in list(S.AZ12):
        S.ZSIGN.setdefault(s, 1.0)
    K._IDX.clear()
    K._IDX.update(cache_index())
    K.BED_YEAR, K.DAY0, K.DAY1 = YEAR, DAY0, DAY1
    K.bed_days = bed_days
    K.CLEAN11 = stations()
    K.DET_SITES = stations()

    null = (os.environ.get("L_NULL_SHIFT", "") == "1") if null is None else null
    if null:
        days = bed_days()
        di = {d: i for i, d in enumerate(days)}
        off = null_offsets(stations())
        _S["OFF"] = off
        base = K.load_day

        def _shifted(sta, date):
            """The station's data from day d+off[sta], restamped with day d's t0 -- so every
            downstream stage (autocorrelation, grouping, stacking, matched filter, location) is
            byte-identical code on a network whose inter-station coherence has been destroyed."""
            if date not in di or sta not in off:
                return base(sta, date)
            src = days[(di[date] + off[sta]) % len(days)]
            v = base(sta, src)
            if v is None:
                return None
            return float(pd.Timestamp(date).timestamp()), v[1]

        K.load_day = _shifted
        log(f"[null] day-shift active, offsets {json.dumps({k: off[k] for k in list(off)[:6]})} ...")
    _S["ACT"] = True
    return K


def bed_summary():
    A = admitted()
    keep = A[~A.health_drop]
    return dict(subregion=SUBREGION, day0=DAY0, day1=DAY1,
                n_admitted=int(len(A)), n_health_dropped=int(A.health_drop.sum()),
                health_dropped=list(A[A.health_drop].netsta),
                n_stations=int(len(keep)),
                by_net={k: int(v) for k, v in keep.groupby("net").size().items()},
                stations={k: sorted(v) for k, v in
                          keep.groupby("net").sta.apply(list).items()})
