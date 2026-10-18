#!/usr/bin/env python
"""SEALED 2011 re-fire bed adapter (2011-07-20 .. 2011-09-10, N-Olympic).  Prefix `m_`.

This is the CONFIRMATION step for the WA 2010 candidate families.  Nothing is tuned: the
templates, the matched filter, the 8 x MAD-of-the-block rule, the stacking/alignment code and the
locator are agent K's, reached through `l_common.activate()` exactly as the 2010 run reached them.
The only things this module changes are BED facts:

  * the window is 2011-07-20 .. 2011-09-10 (the sealed tranche T3);
  * the station universe is the 2010 probe-admitted, 2010-health-kept set INTERSECTED with the
    stations that exist in 2011 and pass the SAME station-level health bar (> 30 % vetoed
    station-hours) on the 2011 census -- a station absent in 2011 is simply dropped, never
    substituted;
  * the tremor hours/days and the quiet days are recomputed for 2011 with the IDENTICAL rules
    (>= 5 PNSN epicentres within 50 km of an active cell; quiet = no PNSN tremor within 100 km of
    the region centre and >= 3 days from any such day).

Verified: the same code reproduces the 2010 numbers exactly (454 candidate hours on 29 days).
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
RES = f"{STUDY}/results"
CACHE = f"{STUDY}/cache20hz"
sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import l_common as L                                                        # noqa: E402

YEAR = 2011
DAY0, DAY1 = "2011-07-20", "2011-09-10"
TRANCHE = "T3"
VETO = f"{RES}/wa_station_hour_veto_2011.csv"
VETO_2010 = f"{RES}/wa_station_hour_veto_2010.csv"
OUT = os.environ.get("M_OUT", f"{RES}/m_2011")

# ---- the rule constants, copied from w_common (NOT re-tuned) -----------------------------------
AC_HOUR_MIN_EPI = 5          # candidate hour: >= 5 PNSN epicentres ...
AC_HOUR_EPI_KM = 50.0        # ... within 50 km of any ACTIVE CELL of the sub-region
QUIET_TREMOR_KM = 100.0      # quiet day: no PNSN tremor within 100 km of the region centre ...
QUIET_MIN_DAYS = 3           # ... and >= 3 days from any such tremor day
GRID_LAT0, GRID_LON0 = 48.00, -123.20
VETO_STATION_MAX = L.VETO_STATION_MAX

_S = {}
log = L.log


def bed_days():
    return [d.strftime("%Y-%m-%d") for d in pd.date_range(DAY0, DAY1)]


# ------------------------------------------------------------------ stations -------------------
def avail():
    """The sealed 2011 tranche as the download manifest records it (status == data)."""
    M = pd.read_csv(f"{RES}/download_manifest.csv", keep_default_na=False)
    T = M[(M.tranche == TRANCHE) & (M.status == "data")]
    T = T.drop_duplicates(["net", "sta", "day"], keep="last").copy()
    T["netsta"] = T.net + "." + T.sta
    return T


def _veto_frac(path):
    if not os.path.exists(path):
        return {}
    V = pd.read_csv(path, usecols=["netsta", "veto"])
    return V.groupby("netsta").veto.mean().to_dict()


def station_frame():
    """One row per 2010 probe-admitted station with both years' health fractions and the
    2011 availability, so the drop reason of every station is on the record."""
    if "SF" in _S:
        return _S["SF"]
    A = pd.read_csv(L.PROBE)
    A = A[A.admit].copy()
    f10, f11 = _veto_frac(VETO_2010), _veto_frac(VETO)
    have = set(avail().netsta)
    A["veto_frac_2010"] = [float(f10.get(n, 0.0)) for n in A.netsta]
    A["veto_frac_2011"] = [float(f11.get(n, np.nan)) for n in A.netsta]
    A["in_2011"] = [n in have for n in A.netsta]
    A["drop_2010_health"] = A.veto_frac_2010 > VETO_STATION_MAX
    A["drop_2011_health"] = A.veto_frac_2011.fillna(0.0) > VETO_STATION_MAX
    A["keep_2011"] = A.in_2011 & ~A.drop_2010_health & ~A.drop_2011_health
    _S["SF"] = A.sort_values("netsta").reset_index(drop=True)
    return _S["SF"]


def stations():
    """The 2011 station universe: 2010-admitted AND 2010-health-kept AND present in 2011 AND
    2011-health-kept.  Bare station codes, as every k_ helper addresses them."""
    A = station_frame()
    return sorted(A[A.keep_2011].sta)


def stations_2010():
    A = station_frame()
    return sorted(A[~A.drop_2010_health].sta)


# ------------------------------------------------------------------ tremor / quiet -------------
def pnsn(pad_days=3):
    T = pd.read_csv(L.PNSN, usecols=["time", "lat", "lon"])
    T["time"] = pd.to_datetime(T.time, errors="coerce").astype("datetime64[ns]")
    T = T.dropna(subset=["time"])
    t0 = pd.Timestamp(DAY0) - pd.Timedelta(days=pad_days)
    t1 = pd.Timestamp(DAY1) + pd.Timedelta(days=pad_days + 1)
    T = T[(T.time >= t0) & (T.time < t1)].reset_index(drop=True)
    T["t"] = T.time.astype("int64") / 1e9
    return T


def candidate_hours(day0=None, day1=None):
    """w_common.candidate_hours' rule, evaluated for an arbitrary window: hours with
    >= 5 PNSN epicentres within 50 km of an active cell of the sub-region."""
    import s_common as S
    day0, day1 = day0 or DAY0, day1 or DAY1
    AC = L.active_cells()
    T = pd.read_csv(L.PNSN, usecols=["time", "lat", "lon"])
    T["time"] = pd.to_datetime(T.time, errors="coerce").astype("datetime64[ns]")
    T = T.dropna(subset=["time"])
    days = {d.strftime("%Y-%m-%d") for d in pd.date_range(day0, day1)}
    T = T[T.time.dt.strftime("%Y-%m-%d").isin(days)]
    if not len(T):
        return pd.DataFrame(columns=["date", "hour", "n_epi"])
    tla, tlo = T.lat.to_numpy(), T.lon.to_numpy()
    cla, clo = AC.lat.to_numpy(), AC.lon.to_numpy()
    near = np.zeros(len(T), bool)
    for i in range(0, len(cla), 200):
        sl = slice(i, min(i + 200, len(cla)))
        d = S.epi_km(cla[sl][:, None], clo[sl][:, None], tla[None, :], tlo[None, :])
        near |= (d <= AC_HOUR_EPI_KM).any(0)
    T = T[near]
    key = pd.DataFrame(dict(date=T.time.dt.strftime("%Y-%m-%d").to_numpy(),
                            hour=T.time.dt.hour.to_numpy()))
    H = key.groupby(["date", "hour"], as_index=False).size().rename(columns={"size": "n_epi"})
    return H[H.n_epi >= AC_HOUR_MIN_EPI].sort_values(["date", "hour"]).reset_index(drop=True)


def tremor_dates():
    return sorted(set(candidate_hours().date))


def quiet_frame():
    """Every bed day with its distance (days) to the nearest PNSN tremor day within
    QUIET_TREMOR_KM of the region centre."""
    import s_common as S
    T = pd.read_csv(L.PNSN, usecols=["time", "lat", "lon"])
    T["time"] = pd.to_datetime(T.time, errors="coerce").astype("datetime64[ns]")
    T = T.dropna(subset=["time"])
    d = S.epi_km(T.lat.to_numpy(), T.lon.to_numpy(), GRID_LAT0, GRID_LON0)
    td = sorted({t.date() for t in T.time[d <= QUIET_TREMOR_KM]})
    tds = np.array([pd.Timestamp(x).value for x in td], dtype="int64")
    trem = set(tremor_dates())
    rows = []
    for dd in bed_days():
        t = pd.Timestamp(dd)
        gap = float(np.abs(tds - t.value).min() / 86400e9) if len(tds) else 1e9
        rows.append(dict(date=dd, days_from_tremor=gap, quiet=bool(gap >= QUIET_MIN_DAYS),
                         tremor_day=bool(dd in trem)))
    return pd.DataFrame(rows)


def quiet_dates():
    Q = quiet_frame()
    return sorted(Q[Q.quiet].date)


def nontremor_dates():
    Q = quiet_frame()
    return sorted(Q[~Q.tremor_day].date)


# ------------------------------------------------------------------ activation -----------------
def activate():
    """Rebind k_common onto the 2011 bed: l_common's station table / borehole azimuths /
    net-qualified cache index (they are year-agnostic), then the 2011 window and station set."""
    K = L.activate(null=False)     # the 2011 re-fire is never day-shifted
    K.BED_YEAR, K.DAY0, K.DAY1 = YEAR, DAY0, DAY1
    K.bed_days = bed_days
    K.CLEAN11 = stations()
    K.DET_SITES = K.CLEAN11
    return K


def bed_summary():
    A = station_frame()
    H = candidate_hours()
    Q = quiet_frame()
    return dict(window=[DAY0, DAY1], year=YEAR, tranche=TRANCHE,
                n_admitted_2010=int(len(A)),
                n_kept_2010=int((~A.drop_2010_health).sum()),
                n_present_2011=int((A.in_2011 & ~A.drop_2010_health).sum()),
                n_absent_2011=int(((~A.in_2011) & ~A.drop_2010_health).sum()),
                absent_2011=sorted(A[(~A.in_2011) & (~A.drop_2010_health)].netsta),
                n_drop_2011_health=int((A.in_2011 & ~A.drop_2010_health
                                        & A.drop_2011_health).sum()),
                drop_2011_health=sorted(A[A.in_2011 & ~A.drop_2010_health
                                          & A.drop_2011_health].netsta),
                n_stations=int(A.keep_2011.sum()),
                by_net={k: int(v) for k, v in
                        A[A.keep_2011].groupby("net").size().items()},
                stations=sorted(A[A.keep_2011].netsta),
                n_candidate_hours=int(len(H)), n_tremor_days=int(H.date.nunique()),
                tremor_dates=tremor_dates(),
                n_quiet_days=int(Q.quiet.sum()), quiet_dates=quiet_dates(),
                n_nontremor_days=int((~Q.tremor_day).sum()))
