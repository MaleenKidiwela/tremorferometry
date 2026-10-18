#!/usr/bin/env python
"""Second production download queue for the WA extension (user-approved 2026-10-02): tiers C + D of
results/WA_PROD_PLAN.md §4 and the useful TEMPORARY-network data of results/wa_temp_networks/INVENTORY.md.
Runs beside scripts/wa_dl.py (tiers A + B, results/wa_dl/); the two never fetch the same station-day (see `skip`).

  tier C   SW carriers 2012-01-01 .. 2026-09-30, while live: UW.CPW RRHS PCFR SMW VVHS RER, PB.B019 B941
  tier D   UW 3-component stations whose first metadata node starts >= 2017-01-01, <= 10 km from a WA family
           (v1.1 representatives, pos_lat/lon else approx_lat/lon) and live >= 365 days in the span (the plan's 22)
  tier D2  the same rule for the short UW deployments (live < 365 days: BES*, BST*, SDSN/SDSS ...)
  tier T   temporary networks: stations <= 30 km from a family, plus the nearest temporary station(s) of every
           family that has < 3 usable stations within 30 km (rule in results/wa_temp_networks/INVENTORY.md);
           Z5 never (restricted); closed/non-seismic/accelerometer-only networks never.

Conventions = scripts/wa_dl.py / wa_fetch.py: live, band and location codes from the bridge-node tables
(tiers C, D, D2) or from the FDSN channel epochs in results/wa_temp_networks/channels.csv (tier T: the preferred
seismometer group of the epoch, all its components); provider EARTHSCOPE (CN falls back to NRCan); one UTC day per
request; tmp + rename writes; `<DOY>.mseed.nd` = verified empty.

Output root (a NEW trace root, NOT yet in wa_common.TRACE_ROOTS; append it after the GO/NO-GO run in
results/wa_gonogo/ has finished, it only holds days missing from every other root):
  lfe_extension_wa/traces_dl2/<NET.STA>/<YYYY>/<DOY>.mseed
A day already raw in another root is not planned (TRACE_ROOTS) or is SYMLINKED (other projects' roots, status
`linked`). A day in results/wa_dl/plan.csv (tiers A + B) is never planned here.

  python scripts/wa_dl2.py plan        -> results/wa_dl2/plan.csv, plan_summary.csv
  python scripts/wa_dl2.py run [--workers 8 --rate 3 --cap-tb 2.0]
  python scripts/wa_dl2.py summary     -> results/wa_dl2/completeness_by_station.csv
Manifest (append-only, resumable): results/wa_dl2/manifest.csv
  tier, netsta, date, status (data|nodata|linked|skipped_ab|error), bytes, source, note, utc
"""
from __future__ import annotations

import argparse
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, FIRST_COMPLETED, wait

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W                                                     # noqa: E402

OUT = f"{W.RES}/wa_dl2"
DL = f"{W.STUDY}/traces_dl2"
LOG = f"{W.STUDY}/logs/wa_dl2.log"
MAN = f"{OUT}/manifest.csv"
PLAN = f"{OUT}/plan.csv"
AB_PLAN, AB_MAN = f"{W.RES}/wa_dl/plan.csv", f"{W.RES}/wa_dl/manifest.csv"
TINV = f"{W.RES}/wa_temp_networks"
D0, D1 = "2012-01-01", "2026-09-30"
TIER_C = ["UW.CPW", "UW.RRHS", "UW.PCFR", "UW.SMW", "UW.VVHS", "PB.B019", "PB.B941", "UW.RER"]
ORDER = {"C": 0, "D": 1, "T": 2, "D2": 3}           # under the 2.0 TB cap: C, D, then T (ranked), D2 last
OTHER_RAW = [f"{W.ROOT}/{r}" for r in ["lfe_picker/traces", "lfe_south/traces", "lfe_array_discovery/traces_sealed",
                                         "lfe_array_discovery/traces_eq", "data/waveforms"]]
NRCAN = "https://earthquakescanada.nrcan.gc.ca"
os.makedirs(OUT, exist_ok=True)


def log(*a):
    W.log(LOG, *a)


def families():
    F = pd.read_csv(f"{W.ROOT}/data/catalogue/families_v1_2.csv", low_memory=False)
    F = F[F.region.isin(W.REGIONS) & (F.representative_v1_1 == True)].copy()
    F["flat"] = F.pos_lat.fillna(F.approx_lat)
    F["flon"] = F.pos_lon.fillna(F.approx_lon)
    return F[F.flat.notna()].reset_index(drop=True)


def dmin(F, la, lo):
    return float(np.nanmin(np.hypot((F.flon.values - lo) * 111.195 * np.cos(np.radians((F.flat.values + la) / 2)),
                                    (F.flat.values - la) * 111.195)))


def nodes():
    fs = [f"{W.EXT}/{r}_bridge_nodes_{k}.csv" for r in ["NW_WA", "SW_WA", "OLY2010"] for k in ["strict", "lenient"]]
    N = pd.concat([pd.read_csv(f, keep_default_na=False, dtype={"loc": str}) for f in fs]).drop_duplicates("node")
    N["netsta"] = N.net + "." + N.sta
    return N


def _ok(p):
    try:
        return os.path.getsize(p) > 0
    except OSError:
        return False


def dl_path(ns, d):
    ts = pd.Timestamp(d)
    return f"{DL}/{ns}/{ts.year}/{int(ts.dayofyear):03d}.mseed"


def tier_d(N, F):
    """tier D / D2 station lists from the bridge nodes (UW, 3-C, first node >= 2017, <= 10 km)."""
    G = N.groupby("netsta").agg(lat=("lat", "first"), lon=("lon", "first"), first=("start", "min"),
                                ncomp=("n_comp", "max")).reset_index()
    G = G[G.netsta.str.startswith("UW.") & (G["first"] >= "2017-01-01") & (G.ncomp >= 3)].copy()
    G["d"] = [dmin(F, a, b) for a, b in zip(G.lat, G.lon)]
    G = G[G.d <= 10.0]
    days = pd.date_range(D0, D1).strftime("%Y-%m-%d")
    live = {}
    for ns in G.netsta:
        n = N[N.netsta == ns]
        live[ns] = sum(bool(((n.start <= d) & (n.end > d)).any()) for d in days)
    G["live"] = G.netsta.map(live)
    return sorted(G[G.live >= 365].netsta), sorted(G[(G.live > 0) & (G.live < 365)].netsta), G


def temp_rows():
    """tier T station-day rows from the inventory selection (results/wa_temp_networks/selected_for_download.csv)."""
    f = f"{TINV}/selected_for_download.csv"
    if not os.path.exists(f):
        log("[plan] no temporary selection file yet; tier T empty")
        return []
    S = pd.read_csv(f, keep_default_na=False, dtype={"loc": str})
    C = pd.read_csv(f"{TINV}/channels.csv", keep_default_na=False, dtype={"loc": str})
    out = []
    for r in S.itertuples():
        c = C[(C.net == r.net) & (C.sta == r.sta) & (C.cha.str[:2] == r.group) & (C["loc"] == r.loc)]
        for ts in pd.date_range(max(D0, r.start), min(D1, r.end[:10])):
            d = ts.strftime("%Y-%m-%d")
            dn = (ts + pd.Timedelta(days=1)).strftime("%Y-%m-%d")
            if not ((c.start < dn) & (c.end > d)).any():
                continue
            out.append(dict(tier="T", netsta=f"{r.net}.{r.sta}", date=d, bands=r.group, locs=r.loc,
                            est_mb=float(r.mb_day), provider=r.provider, prio=int(r.prio)))
    return out


def plan():
    N, F = nodes(), families()
    C = pd.read_csv(f"{W.RES}/wa_prod_plan/download_candidates_2012_2026.csv").set_index("ns")
    D, D2, G = tier_d(N, F)
    G.to_csv(f"{OUT}/tier_d_stations.csv", index=False)
    log(f"[plan] tier D {len(D)}: {' '.join(D)}")
    log(f"[plan] tier D2 {len(D2)}: {' '.join(D2)}")
    ab = set()
    if os.path.exists(AB_PLAN):
        A = pd.read_csv(AB_PLAN, usecols=["netsta", "date"])
        ab = set(zip(A.netsta, A.date))
    days = pd.date_range(D0, D1)
    rows = []
    for tier, sts in [("C", TIER_C), ("D", D), ("D2", D2)]:
        for ns in sts:
            n = N[N.netsta == ns]
            for ts in days:
                d = ts.strftime("%Y-%m-%d")
                live = n[(n.start <= d) & (n.end > d)]
                if len(live):
                    rows.append(dict(tier=tier, netsta=ns, date=d, bands=",".join(sorted(set(live.band))),
                                     locs=",".join(sorted(set(live["loc"]))), est_mb=float(C.mbday.get(ns, 15.0)),
                                     provider="EARTHSCOPE", prio=0))
    rows += temp_rows()
    P = pd.DataFrame(rows)
    keep, link, nd = [], [], []
    for ns, d in zip(P.netsta, P.date):
        ts = pd.Timestamp(d)
        rel = f"{ns}/{ts.year}/{int(ts.dayofyear):03d}.mseed"
        on = next((f"{r}/{rel}" for r in W.TRACE_ROOTS if _ok(f"{r}/{rel}")), "")
        keep.append(not on and (ns, d) not in ab)
        link.append(next((os.path.realpath(f"{r}/{rel}") for r in OTHER_RAW if _ok(f"{r}/{rel}")), ""))
        nd.append(next((f"{r}/{rel}.nd" for r in W.TRACE_ROOTS if os.path.exists(f"{r}/{rel}.nd")), ""))
    P["link_src"], P["known_nd"] = link, nd
    n_ab = sum((s, d) in ab for s, d in zip(P.netsta, P.date))
    P = P[keep].reset_index(drop=True)
    P.to_csv(PLAN, index=False)
    S = P.groupby(["tier", "netsta"]).agg(days=("date", "size"), first=("date", "min"), last=("date", "max"),
                                         linkable=("link_src", lambda v: int((v != "").sum())),
                                         est_gb=("est_mb", lambda v: v.sum() / 1e3)).reset_index()
    S.to_csv(f"{OUT}/plan_summary.csv", index=False)
    T = S.groupby("tier").agg(stations=("netsta", "size"), days=("days", "sum"), est_tb=("est_gb", "sum"))
    T["est_tb"] = T.est_tb / 1e3
    log(f"[plan] {len(P)} station-days on {P.netsta.nunique()} stations (dropped {n_ab} already in the A+B plan); "
        f"by tier: {T.round(3).to_dict('index')}; linkable {int((P.link_src != '').sum())}; "
        f"est total {S.est_gb.sum() / 1e3:.3f} TB")
    print(T.to_string())
    return P


_CL = {}


def client(prov):
    from obspy.clients.fdsn import Client
    k = (prov, threading.get_ident())
    if k not in _CL:
        _CL[k] = Client(prov, timeout=240)
    return _CL[k]


class Rate:
    """Global limiter: at most `rps` request starts per second across all threads."""

    def __init__(self, rps):
        self.dt, self.t, self.lk = 1.0 / rps, 0.0, threading.Lock()

    def wait(self):
        with self.lk:
            now = time.time()
            t = max(now, self.t + self.dt)
            self.t = t
        if t > now:
            time.sleep(t - now)


def _nodata(msg):
    return ("No data" in msg or "204" in msg or "NoData" in msg or "FDSNNoDataException" in msg)


class ABGuard:
    """re-reads the A+B manifest every few minutes: a station-day it already holds is never fetched here."""

    def __init__(self):
        self.t, self.s, self.lk = 0.0, set(), threading.Lock()

    def has(self, ns, d):
        with self.lk:
            if time.time() - self.t > 300:
                try:
                    M = pd.read_csv(AB_MAN, usecols=["netsta", "date", "status"], keep_default_na=False)
                    M = M[M.status.isin(["data", "linked"])]
                    self.s = set(zip(M.netsta, M.date))
                except Exception:
                    pass
                self.t = time.time()
        return (ns, d) in self.s


def fetch(r, rate, guard):
    from obspy import UTCDateTime
    ns, d = r.netsta, r.date
    p = dl_path(ns, d)
    base = dict(tier=r.tier, netsta=ns, date=d)
    if _ok(p):
        return dict(base, status="data", bytes=os.path.getsize(p), source="resumed", note="")
    if os.path.exists(p + ".nd"):
        return dict(base, status="nodata", bytes=0, source="resumed", note="")
    if guard.has(ns, d):
        return dict(base, status="skipped_ab", bytes=0, source="wa_dl", note="held by the A+B queue")
    if isinstance(r.link_src, str) and r.link_src:
        os.makedirs(os.path.dirname(p), exist_ok=True)
        if os.path.lexists(p):
            os.remove(p)
        os.symlink(r.link_src, p)
        return dict(base, status="linked", bytes=0, source=r.link_src, note="symlink, no download")
    net, sta = ns.split(".")
    t0 = UTCDateTime(d)
    cha = ",".join(b + "?" for b in r.bands.split(","))
    loc = ",".join(x if x else "--" for x in str(r.locs).split(",")) if isinstance(r.locs, str) else "--"
    p0 = r.provider if isinstance(r.provider, str) and r.provider else "EARTHSCOPE"
    provs = [p0] + ([NRCAN] if net in ("CN", "C8") and p0 != NRCAN else [])
    last = ""
    for prov in provs:
        for k in range(4):
            try:
                rate.wait()
                st = client(prov).get_waveforms(net, sta, loc, cha, t0, t0 + 86400)
                if len(st) == 0:
                    break
                os.makedirs(os.path.dirname(p), exist_ok=True)
                st.write(p + ".tmp", format="MSEED")
                os.replace(p + ".tmp", p)
                return dict(base, status="data", bytes=os.path.getsize(p), source=prov,
                            note=",".join(sorted({tr.id for tr in st})))
            except Exception as e:
                last = f"{type(e).__name__}:{str(e)[:120]}".replace("\n", " ")
                if _nodata(last) or (prov != "EARTHSCOPE" and "404" in last):
                    break
                time.sleep(10 * (k + 1) * (3 if ("429" in last or "503" in last) else 1))
    if last and not (_nodata(last) or "404" in last):
        return dict(base, status="error", bytes=0, source=",".join(provs), note=last)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p + ".nd", "w").close()
    return dict(base, status="nodata", bytes=0, source=",".join(provs), note=last[:80])


def _append(rows):
    if rows:
        D = pd.DataFrame(rows)
        D["utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        D.to_csv(MAN, mode="a", header=not os.path.exists(MAN), index=False)


def manifest():
    if not os.path.exists(MAN):
        return pd.DataFrame(columns=["tier", "netsta", "date", "status", "bytes", "source", "note", "utc"])
    return pd.read_csv(MAN, keep_default_na=False).drop_duplicates(["netsta", "date"], keep="last")


def dl_bytes():
    tot = 0
    for dp, _, fs in os.walk(DL):
        for f in fs:
            fp = os.path.join(dp, f)
            if not os.path.islink(fp):
                try:
                    tot += os.path.getsize(fp)
                except OSError:
                    pass
    return tot


def heartbeat(msg):
    with open(f"{OUT}/WA_DL2_HEARTBEAT", "w") as fh:
        fh.write(f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {msg}\n")


def run(a):
    P = pd.read_csv(PLAN, keep_default_na=False, dtype={"locs": str})
    if a.stations:
        P = P[P.netsta.isin(a.stations.split(","))]
    if a.days:
        P = P[P.date.isin(a.days.split(","))]
    if a.tiers:
        P = P[P.tier.isin(a.tiers.split(","))]
    M = manifest()
    fin_st = ["data", "nodata", "linked", "skipped_ab"]
    done = set(zip(M[M.status.isin(fin_st)].netsta, M[M.status.isin(fin_st)].date))
    if a.retry_errors_only:
        err = set(zip(M[M.status == "error"].netsta, M[M.status == "error"].date))
        P = P[[(s, d) in err for s, d in zip(P.netsta, P.date)]]
    P = P[[(s, d) not in done for s, d in zip(P.netsta, P.date)]]
    if "prio" not in P:
        P["prio"] = 0
    P["prio"] = pd.to_numeric(P.prio, errors="coerce").fillna(0)
    P = P.assign(_t=P.tier.map(ORDER)).sort_values(["_t", "prio", "date", "netsta"])  # tier, T rank, oldest first
    if a.only:
        P = P.head(a.only)
    cap = a.cap_tb * 1e12
    have = dl_bytes()
    log(f"[run] {len(P)} station-days to do (tiers {sorted(P.tier.unique())}), workers {a.workers}, "
        f"rate {a.rate}/s, traces_dl2 now {have / 1e9:.1f} GB, cap {a.cap_tb} TB")
    rate, guard = Rate(a.rate), ABGuard()
    rows, n, gb_new, t0 = [], 0, 0.0, time.time()
    it = iter(P.itertuples())
    stop = None
    with ThreadPoolExecutor(a.workers) as ex:
        pend = set()
        while True:
            while stop is None and len(pend) < 2 * a.workers:
                r = next(it, None)
                if r is None:
                    break
                pend.add(ex.submit(fetch, r, rate, guard))
            if not pend:
                break
            fin, pend = wait(pend, return_when=FIRST_COMPLETED)
            for f in fin:
                x = f.result()
                rows.append(x)
                n += 1
                gb_new += x["bytes"] / 1e9
            if have + gb_new * 1e9 > cap and stop is None:
                stop = f"storage cap {a.cap_tb} TB reached: traces_dl2 {(have + gb_new * 1e9) / 1e12:.3f} TB"
                log("[run] STOP " + stop)
            if len(rows) >= 50:
                _append(rows)
                rows = []
            if n % 200 == 0:
                el = time.time() - t0
                msg = (f"{n}/{len(P)}  +{gb_new:.1f} GB  {n / max(el, 1):.2f} sd/s  "
                       f"ETA {(len(P) - n) / max(n / max(el, 1), 1e-6) / 3600:.1f} h")
                log("[run] " + msg)
                heartbeat(msg)
    _append(rows)
    log(f"[run] finished {n} station-days, +{gb_new:.1f} GB in {(time.time() - t0) / 3600:.2f} h")
    heartbeat(f"finished {n}")
    summary()
    if stop:
        with open(f"{OUT}/STOPPED_STORAGE", "w") as fh:
            fh.write(stop + "\n")
        sys.exit(3)


def summary():
    P = pd.read_csv(PLAN, keep_default_na=False)
    M = manifest()
    J = P.merge(M[["netsta", "date", "status", "bytes"]], on=["netsta", "date"], how="left")
    J["status"] = J.status.replace("", pd.NA).fillna("pending")
    J["bytes"] = pd.to_numeric(J.bytes, errors="coerce").fillna(0)
    out = []
    for (tier, ns), g in J.groupby(["tier", "netsta"]):
        c = g.status.value_counts().to_dict()
        have = c.get("data", 0) + c.get("linked", 0)
        out.append(dict(tier=tier, netsta=ns, planned=len(g), data=c.get("data", 0), linked=c.get("linked", 0),
                        nodata=c.get("nodata", 0), skipped_ab=c.get("skipped_ab", 0), error=c.get("error", 0),
                        pending=c.get("pending", 0), gb=round(g.bytes.sum() / 1e9, 2),
                        got_frac=round(have / max(len(g), 1), 4), first=g.date.min(), last=g.date.max()))
    S = pd.DataFrame(out)
    S.to_csv(f"{OUT}/completeness_by_station.csv", index=False)
    log(f"[summary] {J.status.value_counts().to_dict()}  {J.bytes.sum() / 1e9:.1f} GB downloaded")
    return S


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["plan", "run", "summary"])
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--rate", type=float, default=3.0, help="max request starts per second (all threads)")
    ap.add_argument("--cap-tb", type=float, default=2.0)
    ap.add_argument("--only", type=int, default=0)
    ap.add_argument("--days", default="")
    ap.add_argument("--stations", default="")
    ap.add_argument("--tiers", default="")
    ap.add_argument("--retry-errors-only", action="store_true")
    a = ap.parse_args()
    {"plan": lambda: plan(), "run": lambda: run(a), "summary": lambda: print(summary().to_string(index=False))}[a.cmd]()
