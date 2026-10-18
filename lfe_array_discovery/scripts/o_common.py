#!/usr/bin/env python
"""ARRAY-AWARE N-Olympic bed adapter (2010 discovery, 2011 confirmation).  Prefix `o_`.

Why this exists (PLAN.md sec.17).  The 2010 N-Olympic run treated all 96 three-component
stations as INDEPENDENT sites, but 64 of them sit in six XG clusters <= 1.4 km across and a
further 15 in two XU (CAFE) clusters of the same size.  Summing the per-station correlation over
the members of a 13-element cluster gives that one patch of ground 13 votes in the network
function, and everything those 13 seismometers share -- cultural noise, a local glitch, wind on
one ridge -- is coherent across all of them.  Chestler & Creager (2017), Ghosh et al. (2012) and
Sweet et al. (2014) never did that: they BEAMFORMED each 1-km array into ONE high-SNR trace and
worked with the arrays as arrays.  This module rebuilds the bed on that footing.

SITES, defined by geometry and nothing else
  Single-link clustering of the 96 admitted stations at 1.5 km, then partitioned by NETWORK
  (a borehole and a surface sensor 100 m apart are not one array).  Every (cluster, network)
  sub-group with >= 2 members becomes ONE BEAM site; singletons stay single stations.  The
  partition is frozen to results/o_oly2010/sites.csv on first use and read from there after.

Everything else is the l_common / n_common pattern: NOTHING in scripts/k_*.py is edited, the
k_ modules read their bed facts through k_common at CALL time, and this module rebinds them in
process.  `o_run.py` additionally registers this module as `n_common` in sys.modules so agent N's
five-clause screen (n_screen.py), wide-grid locator (n_locate.py), re-stack/re-locate
(n_relocate.py), re-fire (n_refire.py) and confirmation (n_confirm.py) run here UNCHANGED.

Environment
  O_BED=2010|2011     discovery bed / sealed confirmation bed (default 2010)
  O_NULL_SHIFT=1      whole-chain site-shift null (offsets from O_NULL_SEED, default 11)
  O_OUT               output directory
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
BEAM = f"{STUDY}/beamcache"
ROOT = "/home/jovyan/tremorferometry"
PICKER = f"{ROOT}/lfe_picker"
PNSN = f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv"
CHESTLER = (f"{ROOT}/codex-LFE/results/olympic_reference_discovery_v1/"
            "derived_coordinate_rows_not_a_catalog.csv")
sys.path.insert(0, HERE)
sys.path.insert(0, f"{PICKER}/scripts")

BED = os.environ.get("O_BED", "2010")
OUTROOT = os.environ.get("O_OUTROOT", f"{RES}/o_oly2010")
SITES_CSV = f"{OUTROOT}/sites.csv"
STATICS = f"{OUTROOT}/statics.json"
BEAMLOG = f"{OUTROOT}/beam_hours_{{bed}}.csv"

CLUSTER_KM = 1.5          # single-link distance that defines "one site"
MIN_BEAM_MEMBERS = 3      # an array-site hour with fewer live members is vetoed
VETO_STATION_MAX = 0.30   # WA_V2_INPUTS sec.3, unchanged

if BED == "2010":
    YEAR, DAY0, DAY1 = 2010, "2010-07-20", "2010-09-10"
    HOURS_CSV = f"{RES}/l_wa2010/candidate_hours.csv"
    QUIET_HOURS_CSV = f"{RES}/wa_quiet_hours_2010.csv"
    VETO_YEARS = (2010,)
elif BED == "2011":
    YEAR, DAY0, DAY1 = 2011, "2011-07-20", "2011-09-10"
    HOURS_CSV = f"{RES}/m_2011/candidate_hours.csv"
    QUIET_HOURS_CSV = ""
    VETO_YEARS = (2011,)
else:
    raise SystemExit(f"O_BED must be 2010 or 2011, not {BED!r}")

OUT = os.environ.get("O_OUT", f"{OUTROOT}/{BED}")
_S = {}


def log(*a):
    import time as _t
    print(f"[{_t.strftime('%H:%M:%S')}]", *a, flush=True)


# =============================================================== the 96 member stations ========
def members_2010():
    """The 96 probe-admitted, health-kept stations of the 2010 run -- the member pool."""
    if "M96" not in _S:
        import l_common as L
        A = L.admitted()
        _S["M96"] = A[~A.health_drop].reset_index(drop=True)
    return _S["M96"]


def member_table():
    """lat/lon/net for every member station (all years)."""
    if "MT" not in _S:
        C = pd.read_csv(f"{RES}/channels_2006_2012.csv").dropna(subset=["lat", "lon"])
        C["netsta"] = C.net + "." + C.sta
        keep = set(members_2010().netsta)
        C = C[C.netsta.isin(keep)]
        g = (C.groupby(["net", "sta"], as_index=False)
               .agg(lat=("lat", "median"), lon=("lon", "median"), elev_m=("elev", "median"),
                    fs_native=("fs", "max"), band=("band", "first")))
        g["netsta"] = g.net + "." + g.sta
        _S["MT"] = g.set_index("sta")
    return _S["MT"]


def member_cache_index():
    """{(member_sta, year, doy): path} over the raw 20 Hz station cache."""
    if "MIDX" not in _S:
        idx = {}
        for r in members_2010().itertuples():
            for p in glob.glob(f"{CACHE}/{r.netsta}/*/*.npz"):
                parts = p.split(os.sep)
                idx[(r.sta, int(parts[-2]), int(parts[-1][:3]))] = p
        _S["MIDX"] = idx
    return _S["MIDX"]


def member_az12():
    """{member: (az1, az2)} for members cached as (?HZ, ?H1, ?H2) -- the PB boreholes."""
    if "AZ" not in _S:
        import l_common as L
        _S["AZ"] = L.az12_wa()
    return _S["AZ"]


# =============================================================== site definition ================
def _build_sites():
    import scipy.cluster.hierarchy as h
    from scipy.spatial.distance import squareform
    M = members_2010()
    T = member_table()
    stas = [s for s in M.sta if s in T.index]
    la = np.radians([float(T.lat[s]) for s in stas])
    lo = np.radians([float(T.lon[s]) for s in stas])
    R = 6371.0
    P = np.vstack([R * np.cos(la) * np.cos(lo), R * np.cos(la) * np.sin(lo),
                   R * np.sin(la)]).T
    D = np.sqrt(((P[:, None, :] - P[None, :, :]) ** 2).sum(-1))
    lab = h.fcluster(h.linkage(squareform(D, checks=False), "single"), CLUSTER_KM, "distance")
    rows = []
    for c in sorted(set(lab)):
        idx = [i for i in range(len(stas)) if lab[i] == c]
        for net in sorted({T.net[stas[i]] for i in idx}):
            sub = [i for i in idx if T.net[stas[i]] == net]
            mem = sorted(stas[i] for i in sub)
            ap = float(D[np.ix_(sub, sub)].max()) if len(sub) > 1 else 0.0
            if len(mem) >= 2:
                code = "x" + _code(net, mem)
                kind = "ARRAY"
            else:
                code = mem[0]
                kind = "SINGLE"
            rows.append(dict(site=code, kind=kind, net=net, n_members=len(mem),
                             aperture_km=round(ap, 3),
                             lat=float(np.mean([T.lat[s] for s in mem])),
                             lon=float(np.mean([T.lon[s] for s in mem])),
                             members=",".join(mem)))
    S = pd.DataFrame(rows).sort_values(["kind", "site"]).reset_index(drop=True)
    if S.site.duplicated().any():
        raise SystemExit(f"duplicate site codes: {sorted(S.site[S.site.duplicated()])}")
    os.makedirs(OUTROOT, exist_ok=True)
    S.to_csv(SITES_CSV, index=False)
    return S


def _code(net, mem):
    """A short, collision-free code for an array: the common station prefix, else net+first."""
    p = os.path.commonprefix(mem)
    p = "".join(ch for ch in p if ch.isalpha())
    return (p if len(p) >= 2 else f"{net}{mem[0]}")[:6]


def sites_all():
    """The FROZEN site table (all sites defined on the 2010 member pool)."""
    if "SA" not in _S:
        _S["SA"] = (pd.read_csv(SITES_CSV) if os.path.exists(SITES_CSV) else _build_sites())
    return _S["SA"]


def _present(sta):
    """Does this member have cached data inside this bed's window?"""
    idx = member_cache_index()
    for d in window_days():
        ts = pd.Timestamp(d)
        if (sta, ts.year, ts.dayofyear) in idx:
            return True
    return False


def sites():
    """Sites usable in THIS bed: members restricted to those with data in the window; an ARRAY
    keeps >= 2 surviving members, a SINGLE its one member.  No station is ever substituted."""
    if "S" in _S:
        return _S["S"]
    S = sites_all().copy()
    keep = []
    for r in S.itertuples():
        mem = [m for m in str(r.members).split(",") if _present(m)]
        if (r.kind == "ARRAY" and len(mem) >= 2) or (r.kind == "SINGLE" and len(mem) == 1):
            keep.append(dict(site=r.site, kind=r.kind, net=r.net, n_members=len(mem),
                             n_members_2010=int(r.n_members), aperture_km=float(r.aperture_km),
                             lat=float(r.lat), lon=float(r.lon), members=",".join(mem)))
    _S["S"] = pd.DataFrame(keep).sort_values("site").reset_index(drop=True)
    return _S["S"]


def site_members():
    return {r.site: str(r.members).split(",") for r in sites().itertuples()}


def arrays():
    S = sites()
    return {r.site: str(r.members).split(",") for r in S[S.kind == "ARRAY"].itertuples()}


def stations():
    """The k_ chain's 'stations' ARE the sites."""
    return list(sites().site)


def netsta():
    return {r.site: (r.net + "." + r.site) for r in sites().itertuples()}


# =============================================================== k_common bed facts =============
def station_table():
    if "ST" not in _S:
        S = sites().copy()
        S["netsta"] = S.net + "." + S.site
        S["elev_m"] = 0.0
        S["fs_native"] = 20.0
        S["band"] = "beam"
        _S["ST"] = S.set_index("site")[["net", "lat", "lon", "elev_m", "fs_native", "band",
                                        "netsta"]]
    return _S["ST"]


def az12_wa():
    """Only SINGLE sites can need rotation; a beam is written as Z,N,E already."""
    AZ = member_az12()
    return {r.site: AZ[r.site] for r in sites().itertuples()
            if r.kind == "SINGLE" and r.site in AZ}


def beam_path(site, date):
    ts = pd.Timestamp(date)
    return f"{BEAM}/{site}/{ts.year}/{ts.dayofyear:03d}.npz"


def cache_index():
    """{(site, year, doy): path} -- the BEAM cache for arrays, the station cache for singles."""
    if "IDX" in _S:
        return _S["IDX"]
    idx = {}
    S = sites()
    mi = member_cache_index()
    for r in S.itertuples():
        if r.kind == "ARRAY":
            for p in glob.glob(f"{BEAM}/{r.site}/*/*.npz"):
                parts = p.split(os.sep)
                idx[(r.site, int(parts[-2]), int(parts[-1][:3]))] = p
        else:
            for (sta, y, d), p in mi.items():
                if sta == r.site:
                    idx[(r.site, y, d)] = p
    _S["IDX"] = idx
    return idx


# =============================================================== days and hours =================
def window_days():
    return [d.strftime("%Y-%m-%d") for d in pd.date_range(DAY0, DAY1)]


def candidate_hours():
    H = pd.read_csv(HOURS_CSV)
    return H[["date", "hour"]].copy()


def tremor_hour_set():
    return {(str(r.date), int(r.hour)) for r in candidate_hours().itertuples()}


def tremor_dates():
    return sorted(set(candidate_hours().date))


def quiet_days():
    """2010: the 20 hour-matched, earthquake-screened quiet days of WA_V2_INPUTS sec.2.1 (whole
    calendar days OUTSIDE the window, so a template-level false-alarm RATE is well defined).
    2011: the 7 strictly-quiet days of the sealed window (results/m_2011/quiet_days.csv) -- the
    sealed set contains no more, exactly as agent M recorded (WA_2011_REFIRE_RESULTS sec.1)."""
    if "QD" in _S:
        return _S["QD"]
    if BED == "2010":
        q = sorted(set(pd.read_csv(QUIET_HOURS_CSV).date))
    else:
        Q = pd.read_csv(f"{RES}/m_2011/quiet_days.csv")
        q = sorted(Q[Q.quiet.astype(bool)].date)
    _S["QD"] = q
    return q


def quiet_hours():
    """The pair-stage quiet CONTRAST.  2010: the 40 hour-matched quiet hours.  2011: not used."""
    if BED == "2010":
        return pd.read_csv(QUIET_HOURS_CSV)[["date", "hour"]]
    return pd.DataFrame([dict(date=d, hour=h) for d in quiet_days() for h in range(0, 24, 4)])


def bed_days():
    """Every day this bed reads: tremor days + quiet days (the beam cache is built for these)."""
    return sorted(set(tremor_dates()) | set(quiet_days()))


def beam_days():
    return bed_days()


# =============================================================== health veto ====================
def _member_veto():
    """{(date, hour): set of vetoed MEMBER station codes}."""
    if "MV" in _S:
        return _S["MV"]
    ns = set(members_2010().netsta)
    d = {}
    for y in VETO_YEARS:
        p = f"{RES}/wa_station_hour_veto_{y}.csv"
        if not os.path.exists(p):
            continue
        V = pd.read_csv(p, usecols=["netsta", "date", "hour", "veto"])
        V = V[V.veto & V.netsta.isin(ns)]
        for r in V.itertuples():
            d.setdefault((r.date, int(r.hour)), set()).add(r.netsta.split(".")[1])
    _S["MV"] = d
    return d


def health_veto():
    """{(date, hour): set of vetoed SITE codes}.  A SINGLE site inherits its member's veto; an
    ARRAY site is vetoed only when the beam of that hour was built from < MIN_BEAM_MEMBERS live
    members (recorded by o_beam.py in beam_hours_<bed>.csv), which is the array analogue."""
    if "V" in _S:
        return _S["V"]
    MV = _member_veto()
    S = sites()
    single = {r.site for r in S.itertuples() if r.kind == "SINGLE"}
    out = {}
    for k, bad in MV.items():
        s = bad & single
        if s:
            out[k] = set(s)
    p = BEAMLOG.format(bed=BED)
    if os.path.exists(p):
        B = pd.read_csv(p)
        for r in B[B.n_members < MIN_BEAM_MEMBERS].itertuples():
            out.setdefault((str(r.date), int(r.hour)), set()).add(str(r.site))
    _S["V"] = out
    return out


# =============================================================== geometry / corrections ========
def active_cells():
    import x_common as X
    return X.active_cells(sub=True)


def near_stations(km=60.0):
    import s_common as S
    ST = station_table()
    AC = active_cells()
    return sorted(s for s in stations()
                  if float(np.min(S.epi_km(float(ST.lat[s]), float(ST.lon[s]),
                                           AC.lat.to_numpy(), AC.lon.to_numpy()))) <= km)


def station_corrections():
    """Intraslab-earthquake S station terms.  A SINGLE site takes its member's term; an ARRAY
    site takes the MEDIAN of its members' usable terms (the beam is the members' mean, so the
    mean site term is the right first-order correction) and is FLAGGED if none of its members
    has one.  Nothing is imputed."""
    p = f"{RES}/x_station_corrections.csv"
    corr, flagged = {}, []
    C = pd.read_csv(p).set_index("sta") if os.path.exists(p) else pd.DataFrame()
    SM = site_members()
    for s in stations():
        vals = []
        for m in SM[s]:
            if len(C) and m in C.index and bool(C.loc[m, "usable"]):
                v = float(C.loc[m, "corr_S"])
                if np.isfinite(v):
                    vals.append(v)
        if vals:
            corr[s] = float(np.median(vals))
        else:
            corr[s] = 0.0
            flagged.append(s)
    return corr, flagged


# =============================================================== the null =======================
def null_offsets(stas, seed=None):
    seed = int(os.environ.get("O_NULL_SEED", "11")) if seed is None else seed
    rng = np.random.default_rng(seed)
    n = max(len(tremor_dates()), 2)
    return {s: int(rng.integers(1, n)) for s in sorted(stas)}


# =============================================================== activation =====================
def activate(null=None):
    """Rebind k_common onto the ARRAY bed.  Idempotent.  Returns the k_common module."""
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

    null = (os.environ.get("O_NULL_SHIFT", "") == "1") if null is None else null
    if null:
        td, qd = tremor_dates(), quiet_days()
        part = {d: ("t", i) for i, d in enumerate(td)}
        part.update({d: ("q", i) for i, d in enumerate(qd)})
        lists = {"t": td, "q": qd}
        off = null_offsets(stations())
        _S["OFF"] = off
        base = K.load_day

        def _shifted(site, date):
            """This site's data from another day of the SAME partition, restamped with day d's
            t0.  Inter-SITE coherence is destroyed; the intra-array beam, the noise, the glitches
            and the geometry are untouched, and every downstream stage is byte-identical code."""
            if date not in part or site not in off:
                return base(site, date)
            kind, i = part[date]
            Lst = lists[kind]
            if len(Lst) < 2:
                return base(site, date)
            v = base(site, Lst[(i + off[site]) % len(Lst)])
            if v is None:
                return None
            return float(pd.Timestamp(date).timestamp()), v[1]

        K.load_day = _shifted
        log(f"[null] site-shift active (within-partition), "
            f"{json.dumps({k: off[k] for k in list(off)[:6]})} ...")
    _S["ACT"] = True
    return K


def bed_summary():
    S = sites()
    A = S[S.kind == "ARRAY"]
    corr, flagged = station_corrections()
    H = candidate_hours()
    return dict(
        bed=BED, day0=DAY0, day1=DAY1, subregion="N-Olympic",
        n_member_stations=int(len(members_2010())),
        n_sites=int(len(S)), n_array_sites=int(len(A)),
        n_single_sites=int((S.kind == "SINGLE").sum()),
        n_members_in_arrays=int(A.n_members.sum()),
        cluster_km=CLUSTER_KM,
        arrays=[dict(site=r.site, net=r.net, n=int(r.n_members),
                     aperture_km=float(r.aperture_km), lat=round(float(r.lat), 4),
                     lon=round(float(r.lon), 4), members=r.members)
                for r in A.itertuples()],
        singles=sorted(S[S.kind == "SINGLE"].site),
        n_tremor_hours=int(len(H)), n_tremor_days=int(H.date.nunique()),
        n_quiet_days=len(quiet_days()), quiet_days=quiet_days(),
        n_intraslab_corr=int(len(corr) - len(flagged)), n_corr_flagged=int(len(flagged)),
        corr_flagged=sorted(flagged))
