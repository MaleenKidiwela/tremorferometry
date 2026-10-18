#!/usr/bin/env python
"""Build every download queue (tranches 2-5, 7) + results/freeze_windows.json.

One CSV per sub-tranche in results/queues/, columns
    tranche,net,sta,loc,band,fs,chain,three_c,date,outroot,group
Station/day eligibility comes ONLY from results/station_flags.csv + channels_all.csv (tranche 1):
a station-day is queued when a 3-C velocity band (>= 20 Hz) is declared on that day.  EHZ-only
stations are skipped and listed in results/skipped_ehz_only.csv - the single exception is the
handful of EXPLICITLY NAMED permanent stations (UW HDW is EHZ/EH2 until 2022-09) which are queued
with three_c=False so they can never enter a 3-C sum by accident.

Quiet days: 2010 from the PNSN catalog (no epicentre within 100 km of any tremor-active Slab2 cell
on the day or +-3 d); 2008 and 2007 have NO tremor catalog (PNSN starts 2010-01-06) so their quiet
days are LITERATURE-quiet: >= 30 d from the published ETS window, recorded as such.
"""
from __future__ import annotations
import os
import sys
import json
import argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d_common as D                                                   # noqa: E402

COLS = ["tranche", "net", "sta", "loc", "band", "fs", "chain", "three_c", "alts", "date",
        "outroot", "group"]
# published ETS windows used for the LITERATURE-quiet rule (no catalog before 2010-01-06)
ETS_LIT = {"2008": ("2008-05-01", "2008-06-05"), "2007": ("2007-01-10", "2007-02-20")}
QUIET_MIN_GAP_D = 30


def days(a, b):
    return [str(d.date()) for d in pd.date_range(a, b, freq="D")]


# ------------------------------------------------------------------ per-day band resolution
def resolve_days(C, net, sta, dlist, allow_z_only=False):
    """-> dict day -> [(loc, band, fs, chain, three_c), ...] in preference order.

    The list matters: many XG/XU stations declare BOTH EH? and SH? for the same epoch while only
    one of them actually has data (XG.CL13 2010-08-20 returns nodata on EH, data on SH), so the
    fetcher must be able to fall back to the next declared band instead of writing a `.nd`.
    """
    g = C[(C.net == net) & (C.sta == sta)]
    g = g[g.chan.str[-1].isin(list("ZNE12")) & g.band.isin(D.BAND_PREF) & (g.fs >= D.FS_MIN)]
    if not len(g):
        return {}
    dd = np.array(dlist)
    best = {}
    for (loc, band), gg in g.groupby(["loc", "band"]):
        bpref = D.BAND_PREF.index(band)
        lpref = D.LOCPREF.index(loc) if loc in D.LOCPREF else 99

        def cov(chars):
            m = np.zeros(len(dd), bool)
            fsd = np.zeros(len(dd))
            for r in gg[gg.chan.str[-1].isin(list(chars))].itertuples():
                k = (dd >= r.start) & (dd <= r.end)
                m |= k
                fsd = np.where(k, np.maximum(fsd, r.fs), fsd)
            return m, fsd
        mz, fz = cov("Z")
        mn, fn = cov("N1")
        me, fe = cov("E2")
        three = mz & mn & me
        fsd = np.maximum(np.maximum(fz, fn), fe)
        ok = three | (mz if allow_z_only else np.zeros(len(dd), bool))
        for i in np.where(ok)[0]:
            key = (0 if three[i] else 1, bpref, lpref)
            best.setdefault(dd[i], []).append((key, loc, band, float(fsd[i]), bool(three[i])))
    return {d: [(v[1], v[2], v[3], D.chain_for(v[3]), v[4]) for v in sorted(vs)]
            for d, vs in best.items()}


def rows_for(C, stations, dlist, tranche, group, outroot, named=frozenset()):
    """stations = [(net, sta), ...]; returns (rows, skipped) with skipped = no usable band."""
    rows, skipped = [], []
    for net, sta in stations:
        r = resolve_days(C, net, sta, dlist, allow_z_only=((net, sta) in named))
        if not r:
            rz = resolve_days(C, net, sta, dlist, allow_z_only=True)
            skipped.append(dict(tranche=tranche, group=group, net=net, sta=sta,
                                reason=("vertical_only_in_window" if rz
                                        else "no_epoch_in_window")))
            continue
        if all(not v[0][4] for v in r.values()) and (net, sta) not in named:
            skipped.append(dict(tranche=tranche, group=group, net=net, sta=sta,
                                reason="vertical_only_in_window"))
            continue
        for d in dlist:
            vs = r.get(d)
            if not vs:
                continue
            v = vs[0]
            alts = ";".join(f"{a[0]}|{a[1]}|{a[2]:g}" for a in vs[1:])
            rows.append(dict(tranche=tranche, net=net, sta=sta, loc=v[0], band=v[1], fs=v[2],
                             chain=v[3], three_c=v[4], alts=alts, date=d, outroot=outroot,
                             group=group))
    return rows, skipped


# ------------------------------------------------------------------ quiet days
def quiet_days_2010(n_want=20):
    """No PNSN epicentre within 100 km of ANY tremor-active cell on the day or +-3 d."""
    from scipy.spatial import cKDTree
    G = D.active_cells()
    A = G[G.tremor_active]
    T = D.tremor_catalog()
    T = T[(T.time >= "2010-04-01") & (T.time < "2010-12-01")]
    tree = cKDTree(D.proj(A.lat.values, A.lon.values))
    dmin = tree.query(D.proj(T.lat.values, T.lon.values))[0]
    T = T.assign(dmin=dmin, day=T.time.dt.floor("D"))
    pools = [("2010-06-01", "2010-07-31"), ("2010-05-01", "2010-05-31"),
             ("2010-10-01", "2010-10-31")]
    ets = set(days(*D.WINDOWS["t2_nolympic2010"]))
    allday = pd.date_range("2010-04-01", "2010-11-30", freq="D")
    n = (T[T.dmin <= D.QUIET_RADIUS_KM].groupby("day").size()
         .reindex(allday, fill_value=0))
    qok = n.rolling(2 * D.QUIET_PAD_D + 1, center=True,
                    min_periods=2 * D.QUIET_PAD_D + 1).max() == 0
    out, used_pools = [], []
    for a, b in pools:
        if len(out) >= n_want:
            break
        pool = [str(d.date()) for d in qok[qok].index
                if a <= str(d.date()) <= b and str(d.date()) not in ets and str(d.date()) not in out]
        if pool:
            used_pools.append(f"{a}..{b} ({len(pool)} quiet)")
        take = pool[:max(0, n_want - len(out))] if len(out) + len(pool) > n_want else pool
        # spread evenly through the pool instead of taking the first ones
        if len(pool) > (n_want - len(out)) > 0:
            idx = np.round(np.linspace(0, len(pool) - 1, n_want - len(out))).astype(int)
            take = [pool[i] for i in sorted(set(idx))]
        out += take
    meta = dict(rule=f"no PNSN epicentre within {D.QUIET_RADIUS_KM:g} km of any tremor-active "
                     f"Slab2 cell on the day or +-{D.QUIET_PAD_D} d",
                n_active_cells=int(A.shape[0]), pools_used=used_pools,
                n_found=len(out), n_wanted=n_want, catalog_quiet=True)
    return sorted(out), meta


def quiet_days_literature(key, pools, n_want):
    """No tremor catalog before 2010-01-06 -> >= 30 d from the published ETS window."""
    e0, e1 = ETS_LIT[key]
    lo = str((pd.Timestamp(e0) - pd.Timedelta(days=QUIET_MIN_GAP_D)).date())
    hi = str((pd.Timestamp(e1) + pd.Timedelta(days=QUIET_MIN_GAP_D)).date())
    out = []
    per = int(np.ceil(n_want / len(pools)))
    for a, b in pools:
        pool = [d for d in days(a, b) if d <= lo or d >= hi]
        if not pool:
            continue
        k = min(per, len(pool), n_want - len(out))
        idx = np.round(np.linspace(0, len(pool) - 1, k)).astype(int)
        out += [pool[i] for i in sorted(set(idx))]
    meta = dict(rule=f">= {QUIET_MIN_GAP_D} d from the published {key} ETS window {e0}..{e1}; "
                     f"LITERATURE-quiet, NOT catalog-quiet (PNSN starts 2010-01-06)",
                pools=[f"{a}..{b}" for a, b in pools], n_found=len(out), n_wanted=n_want,
                catalog_quiet=False)
    return sorted(out), meta


# ------------------------------------------------------------------ station sets
def sets(F):
    """Station lists from the tranche-1 flags (3-C only unless explicitly named)."""
    F = F.copy()
    F["key"] = F.net + "." + F.sta
    three = F[F.three_c]
    S = {}
    S["PB"] = [(r.net, r.sta) for r in three[three.net == "PB"].itertuples()]
    S["UW_named"] = [("UW", s) for s in D.UW_ARRAY_ERA]
    S["CN_named"] = [("CN", s) for s in D.CN_ARRAY_ERA]
    S["XU"] = [(r.net, r.sta) for r in three[three.net == "XU"].itertuples()]
    S["ZK"] = [(r.net, r.sta) for r in three[three.net == "ZK"].itertuples()]
    S["Z3"] = [(r.net, r.sta) for r in three[three.net == "Z3"].itertuples()]
    S["YW"] = [(r.net, r.sta) for r in three[three.net == "YW"].itertuples()]
    xg = three[three.net == "XG"].copy()
    xg["array"] = xg.sta.str[:2]                      # BH CL DR GC LC PA SL TB = the 8 arrays
    xg["span"] = (pd.to_datetime(xg.end) - pd.to_datetime(xg.start)).dt.days
    thin = xg.sort_values(["array", "span", "sta"], ascending=[True, False, True]) \
             .groupby("array").head(1)
    S["XG_thin"] = [(r.net, r.sta) for r in thin.itertuples()]
    S["XG_rest"] = [(r.net, r.sta) for r in xg.itertuples()
                    if (r.net, r.sta) not in set(S["XG_thin"])]
    S["XG_all"] = S["XG_thin"] + S["XG_rest"]
    return S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-quiet-2010", type=int, default=20)
    a = ap.parse_args()
    F = D.station_flags()
    C = pd.read_csv(os.path.join(D.RES, "channels_all.csv"), keep_default_na=False,
                    dtype={"loc": str})
    C["loc"] = C["loc"].astype(str).replace("nan", "")
    C["fs"] = pd.to_numeric(C.fs, errors="coerce").fillna(0.0)
    S = sets(F)
    NAMED = set(S["UW_named"]) | set(S["CN_named"])          # may be queued vertical-only
    q2010, m2010 = quiet_days_2010(a.n_quiet_2010)
    q2008, m2008 = quiet_days_literature("2008", [("2008-02-01", "2008-03-31"),
                                                  ("2008-07-01", "2008-07-31")], 20)
    q2007, m2007 = quiet_days_literature("2007", [("2006-11-01", "2006-11-30"),
                                                  ("2007-04-01", "2007-04-30")], 15)
    d2010 = days(*D.WINDOWS["t2_nolympic2010"])
    d2011 = days(*D.WINDOWS["t3_sealed2011"])
    d2008 = days(*D.WINDOWS["t4_puget2008"])
    d2007 = days(*D.WINDOWS["t5_cafe2007"])

    thin_set = S["XG_thin"] + S["PB"] + S["UW_named"] + S["CN_named"]
    plan = [
        ("T2A", "thin",     thin_set,        d2010,          D.TRACES),
        ("T2B", "quiet",    thin_set,        q2010,          D.TRACES),
        ("T2C", "xg_rest",  S["XG_rest"],    d2010 + q2010,  D.TRACES),
        ("T2D", "xu",       S["XU"],         d2010 + q2010,  D.TRACES),
        ("T3",  "sealed",   S["XG_all"] + S["XU"] + S["PB"] + S["UW_named"] + S["CN_named"],
                                             d2011,          D.TRACES_SEALED),
        ("T4",  "puget2008", S["XU"] + S["ZK"] + S["PB"] + [("UW", "HDW"), ("UW", "GNW")]
                             + S["YW"],      d2008 + q2008,  D.TRACES),
        ("T5",  "cafe2007", S["XU"] + S["PB"] + [("UW", "HDW"), ("UW", "GNW")],
                                             d2007 + q2007,  D.TRACES),
        ("T7",  "extension", D.EXT_ROSTER,   None,           D.TRACES),
    ]
    summary, skipped_all = [], []
    for tr, group, stns, dlist, outroot in plan:
        if tr == "T7":
            allday = days(*D.WINDOWS["t7_extension"])
            s0, s1 = D.WINDOWS["t3_sealed2011"]
            # the sealed re-fire window must exist ONLY in traces_sealed/ - otherwise the
            # permanents' 2011-08 days would sit unsealed in traces/ and the seal is worthless.
            allday = [d for d in allday if not (s0 <= d <= s1)]
            ets = [d for d in allday if d[5:7] in ("08", "09")]
            dlist = ets + [d for d in allday if d[5:7] not in ("08", "09")]
        stns = sorted(set(stns))
        rows, skipped = rows_for(C, stns, dlist, tr, group, outroot, named=NAMED)
        skipped_all += skipped
        Q = pd.DataFrame(rows, columns=COLS)
        if tr == "T7":                                        # ETS months first, then oldest-first
            Q["_k"] = np.where(Q.date.str[5:7].isin(["08", "09"]), 0, 1)
            Q = Q.sort_values(["_k", "date", "net", "sta"]).drop(columns="_k")
        else:
            Q = Q.sort_values(["date", "net", "sta"])
        p = os.path.join(D.QUEUES, f"q_{tr}_{group}.csv")
        Q.to_csv(p, index=False)
        ns = Q[["net", "sta"]].drop_duplicates()
        summary.append(dict(tranche=tr, group=group, queue=os.path.basename(p),
                            n_stations=len(ns), n_stations_3c=int(Q.three_c.sum() > 0 and
                                                                  Q[Q.three_c][["net", "sta"]]
                                                                  .drop_duplicates().shape[0]),
                            n_days=len(set(Q.date)), n_station_days=len(Q),
                            n_rational=int((Q.chain == "rational").sum()),
                            outroot=os.path.basename(outroot)))
        D.log(f"[queue] {tr:4s} {group:10s} {len(ns):4d} stations x {len(set(Q.date)):5d} days "
              f"= {len(Q):7d} station-days ({int((Q.chain=='rational').sum())} rational) -> {p}")
    SU = pd.DataFrame(summary)
    SU.to_csv(os.path.join(D.RES, "queue_summary.csv"), index=False)
    EH = F[F.ehz_only][["net", "sta", "lat", "lon", "band", "fs", "n_comp", "chans", "sensor",
                        "start", "end"]].copy()
    EH["tranche"], EH["group"], EH["reason"] = "ALL", "flags", "ehz_only_never_summed"
    SK = pd.concat([EH, pd.DataFrame(skipped_all)], ignore_index=True)
    SK.to_csv(os.path.join(D.RES, "skipped_ehz_only.csv"), index=False)
    D.log(f"[queue] {len(SK)} station-tranche entries skipped (EHZ-only / no band) "
          f"-> results/skipped_ehz_only.csv")
    # ---------------------------------------------------------------- freeze
    fz = dict(
        study="lfe_array_discovery", role="agent D download+metadata queue",
        written_utc=pd.Timestamp.utcnow().isoformat(),
        box=dict(lat=[D.LAT0, D.LAT1], lon=[D.LON0, D.LON1]),
        networks=D.NETS, band_pref=D.BAND_PREF, fs_cache=D.FS_CACHE,
        site_merge_km=D.SITE_MERGE_KM, max_concurrent_fdsn=D.MAX_CONCURRENT_FDSN,
        windows={k: list(v) for k, v in D.WINDOWS.items()},
        quiet_days=dict(y2010=dict(days=q2010, **m2010), y2008=dict(days=q2008, **m2008),
                        y2007=dict(days=q2007, **m2007)),
        station_sets={k: [f"{n}.{s}" for n, s in v] for k, v in S.items()},
        sealed_window_excluded_from_extension=list(D.WINDOWS["t3_sealed2011"]),
        named_vertical_only_allowed=[f"{n}.{s}" for n, s in sorted(NAMED)],
        extension_roster=[f"{n}.{s}" for n, s in D.EXT_ROSTER],
        tranches=summary,
        cache_chain=dict(frozen="lfe_picker/scripts/o_cache_month.py:cache_one",
                         rational="lfe_picker/scripts/u_cache.py:cache_one_rational",
                         note="every manifest/cache row carries chain=frozen|rational"),
        sha256=D.script_hashes())
    D.jdump(fz, os.path.join(D.RES, "freeze_windows.json"))
    D.log(f"[queue] freeze -> results/freeze_windows.json "
          f"(2010 quiet {len(q2010)}/{a.n_quiet_2010}, 2008 {len(q2008)}/20, 2007 {len(q2007)}/15)")
    print(SU.to_string(index=False))


if __name__ == "__main__":
    main()
