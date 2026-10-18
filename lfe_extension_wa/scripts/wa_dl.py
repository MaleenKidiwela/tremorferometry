#!/usr/bin/env python
"""Production downloads for the WA extension, tiers A + B of results/WA_PROD_PLAN.md §4 (user-approved 2026-10-02).

  tier A  gaps of the 21 core stations 2012-01-01 .. 2026-09-30
          CN.PGC CN.VGZ PB.B001 B003-B007 B009-B011 B926 B013 B014 B017 B018 B943 UW.LRIV DOSE GNW HDW
  tier B  NW/Olympic UW short-period + US.NLWA while live, same span:
          UW.BLN SQM OSD SVOH GMW OTR BHW HTW STW SP2 WISH STOR, US.NLWA
  Tiers C (SW carriers) and D (new 3-C UW) are NOT approved and are not in TIERS.

Same conventions as scripts/wa_fetch.py: live = a node of data/extension_stations/{NW_WA,SW_WA,OLY2010}_bridge_nodes_
{strict,lenient}.csv with start <= day < end (this reproduces the plan's `live` counts in
results/wa_prod_plan/download_candidates_2012_2026.csv); bands = the bands of the live node(s); location = the live
nodes' location codes (avoids pulling duplicate co-located sensors); provider EARTHSCOPE, CN falls back to NRCan;
one UTC day per request, tmp + rename writes, `<DOY>.mseed.nd` = verified empty.

Output root (a NEW trace root, appended to wa_common.TRACE_ROOTS):
  lfe_extension_wa/traces_dl/<NET.STA>/<YYYY>/<DOY>.mseed
A day already present as raw in a root the WA chain does not read (lfe_picker/traces, lfe_south/traces,
lfe_array_discovery/traces_sealed, data/waveforms) is SYMLINKED into traces_dl (status `linked`), not re-fetched.
A day already present in a TRACE_ROOT is not planned at all.

  python scripts/wa_dl.py plan                    -> results/wa_dl/plan.csv, plan_summary.csv
  python scripts/wa_dl.py run [--workers 8 --rate 3 --cap-tb 0.8 --only N]
  python scripts/wa_dl.py summary                 -> results/wa_dl/completeness_by_station.csv
Manifest (append-only, resumable): results/wa_dl/manifest.csv
  tier, netsta, date, status (data|nodata|linked|error), bytes, source, note, utc
"""
from __future__ import annotations

import argparse
import os
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, FIRST_COMPLETED, wait

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W                                                     # noqa: E402

OUT = f"{W.RES}/wa_dl"
DL = f"{W.STUDY}/traces_dl"
LOG = f"{W.STUDY}/logs/wa_dl.log"
MAN = f"{OUT}/manifest.csv"
PLAN = f"{OUT}/plan.csv"
D0, D1 = "2012-01-01", "2026-09-30"
TIERS = {
    "A": ["CN.PGC", "CN.VGZ", "PB.B001", "PB.B003", "PB.B004", "PB.B005", "PB.B006", "PB.B007", "PB.B009",
          "PB.B010", "PB.B011", "PB.B926", "UW.LRIV", "PB.B013", "PB.B014", "PB.B017", "PB.B018", "PB.B943",
          "UW.DOSE", "UW.GNW", "UW.HDW"],
    "B": ["UW.BLN", "UW.SQM", "US.NLWA", "UW.OSD", "UW.SVOH", "UW.GMW", "UW.OTR", "UW.BHW", "UW.HTW", "UW.STW",
          "UW.SP2", "UW.WISH", "UW.STOR"],
}
OTHER_RAW = [f"{W.ROOT}/{r}" for r in ["lfe_picker/traces", "lfe_south/traces", "lfe_array_discovery/traces_sealed",
                                         "data/waveforms"]]
NRCAN = "https://earthquakescanada.nrcan.gc.ca"
os.makedirs(OUT, exist_ok=True)


def log(*a):
    W.log(LOG, *a)


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


def plan():
    N = nodes()
    C = pd.read_csv(f"{W.RES}/wa_prod_plan/download_candidates_2012_2026.csv").set_index("ns")
    days = pd.date_range(D0, D1)
    rows = []
    for tier, sts in TIERS.items():
        for ns in sts:
            n = N[N.netsta == ns]
            for ts in days:
                d = ts.strftime("%Y-%m-%d")
                live = n[(n.start <= d) & (n.end > d)]
                if not len(live):
                    continue
                rel = f"{ns}/{ts.year}/{int(ts.dayofyear):03d}.mseed"
                on = next((f"{r}/{rel}" for r in W.TRACE_ROOTS if _ok(f"{r}/{rel}")), "")
                if on and not on.startswith(DL):
                    continue                                            # chain already reads it
                nd = next((f"{r}/{rel}.nd" for r in W.TRACE_ROOTS if r != DL and os.path.exists(f"{r}/{rel}.nd")), "")
                other = next((os.path.realpath(f"{r}/{rel}") for r in OTHER_RAW if _ok(f"{r}/{rel}")), "")
                rows.append(dict(tier=tier, netsta=ns, date=d, bands=",".join(sorted(set(live.band))),
                                 locs=",".join(sorted(set(live["loc"]))), link_src=other, known_nd=nd,
                                 est_mb=float(C.mbday.get(ns, 15.0))))
    P = pd.DataFrame(rows)
    P.to_csv(PLAN, index=False)
    S = P.groupby(["tier", "netsta"]).agg(days=("date", "size"), first=("date", "min"), last=("date", "max"),
                                         linkable=("link_src", lambda v: int((v != "").sum())),
                                         known_nd=("known_nd", lambda v: int((v != "").sum())),
                                         est_gb=("est_mb", lambda v: v.sum() / 1e3)).reset_index()
    S.to_csv(f"{OUT}/plan_summary.csv", index=False)
    log(f"[plan] {len(P)} station-days on {P.netsta.nunique()} stations; tier A {int((P.tier == 'A').sum())}, "
        f"tier B {int((P.tier == 'B').sum())}; linkable {int((P.link_src != '').sum())}, known-nodata "
        f"{int((P.known_nd != '').sum())}; est {S.est_gb.sum():.0f} GB (fetch part "
        f"{P[(P.link_src == '') & (P.known_nd == '')].est_mb.sum() / 1e6:.3f} TB)")
    print(S.to_string(index=False))


_CL = {}


def client(prov):
    from obspy.clients.fdsn import Client
    k = (prov, threading.get_ident())
    if k not in _CL:
        _CL[k] = Client(prov, timeout=180)
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


def fetch(r, rate):
    """One station-day. Mirrors wa_fetch.fetch, writing under traces_dl."""
    from obspy import UTCDateTime
    ns, d = r.netsta, r.date
    p = dl_path(ns, d)
    base = dict(tier=r.tier, netsta=ns, date=d)
    if _ok(p):
        return dict(base, status="data", bytes=os.path.getsize(p), source="resumed", note="")
    if os.path.exists(p + ".nd"):
        return dict(base, status="nodata", bytes=0, source="resumed", note="")
    if isinstance(r.link_src, str) and r.link_src:
        os.makedirs(os.path.dirname(p), exist_ok=True)
        if os.path.lexists(p):
            os.remove(p)                                                # a dangling link left by a crash
        os.symlink(r.link_src, p)
        return dict(base, status="linked", bytes=0, source=r.link_src, note="symlink, no download")
    # an earlier `.nd` marker in another root (column known_nd) is NOT trusted: it may come from a different band
    # request; the day is re-queried (a no-data answer is cheap)
    net, sta = ns.split(".")
    t0 = UTCDateTime(d)
    cha = ",".join(b + "?" for b in r.bands.split(","))
    loc = ",".join(x if x else "--" for x in str(r.locs).split(",")) if isinstance(r.locs, str) else "--"
    provs = ["EARTHSCOPE"] + ([NRCAN] if net in ("CN", "C8") else [])
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


def run(a):
    P = pd.read_csv(PLAN, keep_default_na=False)
    if a.days:
        P = P[P.date.isin(a.days.split(","))]
    if a.stations:
        P = P[P.netsta.isin(a.stations.split(","))]
    M = manifest()
    done = set(zip(M[M.status.isin(["data", "nodata", "linked", "nodata_known"])].netsta,
                   M[M.status.isin(["data", "nodata", "linked", "nodata_known"])].date))
    if a.retry_errors_only:
        err = set(zip(M[M.status == "error"].netsta, M[M.status == "error"].date))
        P = P[[(s, d) in err for s, d in zip(P.netsta, P.date)]]
    P = P[[(s, d) not in done for s, d in zip(P.netsta, P.date)]]
    P = P.assign(_t=P.tier.map({"A": 0, "B": 1})).sort_values(["_t", "date", "netsta"])   # A first, oldest first
    if a.only:
        P = P.head(a.only)
    cap = a.cap_tb * 1e12
    have = dl_bytes()
    log(f"[run] {len(P)} station-days to do (tiers {sorted(P.tier.unique())}), workers {a.workers}, "
        f"rate {a.rate}/s, traces_dl now {have / 1e9:.1f} GB, cap {a.cap_tb} TB")
    rate = Rate(a.rate)
    rows, n, gb_new, t0 = [], 0, 0.0, time.time()
    tlast = t0
    it = iter(P.itertuples())
    stop = None
    with ThreadPoolExecutor(a.workers) as ex:
        pend = set()
        while True:
            while stop is None and len(pend) < 2 * a.workers:
                r = next(it, None)
                if r is None:
                    break
                pend.add(ex.submit(fetch, r, rate))
            if not pend:
                break
            fin, pend = wait(pend, return_when=FIRST_COMPLETED)
            for f in fin:
                x = f.result()
                rows.append(x)
                n += 1
                gb_new += x["bytes"] / 1e9
            if have + gb_new * 1e9 > cap and stop is None:
                stop = f"storage cap {a.cap_tb} TB reached: traces_dl {(have + gb_new * 1e9) / 1e12:.3f} TB"
                log("[run] STOP " + stop)
            if len(rows) >= 50:
                _append(rows)
                rows = []
            if time.time() - tlast > 60:
                tlast = time.time()
                el = time.time() - t0
                log(f"[run] {n}/{len(P)}  +{gb_new:.1f} GB  {n / max(el, 1):.2f} sd/s  "
                    f"ETA {(len(P) - n) / max(n / max(el, 1), 1e-6) / 3600:.1f} h")
                with open(f"{OUT}/PROGRESS", "w") as fh:
                    fh.write(f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {n}/{len(P)} "
                             f"+{gb_new:.1f}GB\n")
    _append(rows)
    log(f"[run] finished {n} station-days, +{gb_new:.1f} GB in {(time.time() - t0) / 3600:.2f} h")
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
    N = nodes()
    days = pd.date_range(D0, D1).strftime("%Y-%m-%d")
    out = []
    for (tier, ns), g in J.groupby(["tier", "netsta"]):
        n = N[N.netsta == ns]
        live = sum(bool(((n.start <= d) & (n.end > d)).any()) for d in days)
        c = g.status.value_counts().to_dict()
        have = c.get("data", 0) + c.get("linked", 0)
        out.append(dict(tier=tier, netsta=ns, live_days=live, already_on_disk=live - len(g), planned=len(g),
                        data=c.get("data", 0), linked=c.get("linked", 0), nodata=c.get("nodata", 0),
                        nodata_known=c.get("nodata_known", 0), error=c.get("error", 0),
                        pending=c.get("pending", 0), gb=round(g.bytes.sum() / 1e9, 2),
                        complete_frac=round((live - len(g) + have) / max(live, 1), 4),
                        first=g.date.min(), last=g.date.max()))
    S = pd.DataFrame(out)
    S.to_csv(f"{OUT}/completeness_by_station.csv", index=False)
    log(f"[summary] {J.status.value_counts().to_dict()}  {J.bytes.sum() / 1e9:.1f} GB downloaded")
    return S


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["plan", "run", "summary"])
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--rate", type=float, default=3.0, help="max request starts per second (all threads)")
    ap.add_argument("--cap-tb", type=float, default=0.8)
    ap.add_argument("--only", type=int, default=0)
    ap.add_argument("--days", default="")
    ap.add_argument("--stations", default="")
    ap.add_argument("--retry-errors-only", action="store_true")
    a = ap.parse_args()
    {"plan": lambda: plan(), "run": lambda: run(a), "summary": lambda: print(summary().to_string(index=False))}[a.cmd]()
