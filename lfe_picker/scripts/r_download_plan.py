#!/usr/bin/env python
"""P3D step 1 — station/band epoch table + SSE windows + ordered download queue.

Adapted from scripts/o_download_plan.py (itself from ../lfe_extension/scripts/download_plan.py).
Roster = the ten close 3-C permanent stations that post-date Bostock, plus two Z-only helpers.

Outputs (all under results/p3d/):
  station_epochs.csv   net,sta,band,fs,loc,start,end,lat,lon,elev,n_comp,chan_az_dip,mode
                       chan_az_dip = "HHZ:0.0/-90.0;HHN:0.0/0.0;HHE:90.0/0.0"  (ZNE rotation)
  sse_windows.csv      id,core_start,core_end,start,end,n_days,n_tremor  (PNSN northern episodes)
  download_queue.csv   net,sta,band,loc,fs,date,doy,priority,status,path
  download_manifest.csv (header only if absent; the worker appends)
  download_plan_log.txt

Priorities
  P0  SSE windows 2019-2026 (PNSN northern episodes, padded +-5 d)
  P1  one quiet month per year (fewest tremor epicentres in the box; 0 where possible)
  P2  everything else 2019-01-01..2026-09-15, oldest first

Usage: python scripts/r_download_plan.py [--no-fdsn]
"""
import os
import csv
import json
import argparse
from datetime import date, datetime, timedelta, timezone

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(HERE, "results", "p3d")
TRACES = os.path.join(HERE, "traces")
CATALOG = os.path.abspath(os.path.join(HERE, "..", "catalogs", "pnsn_tremor_cascadia_full.csv"))

# ---- roster ------------------------------------------------------------------
ROSTER_3C = [
    ("CN", "SYMB"), ("CN", "CLRS"), ("CN", "PTRF"), ("CN", "GOBB"),
    ("UW", "RNWY"), ("UW", "BUCKS"), ("UW", "TURTL"), ("UW", "LOPEZ"),
    ("UW", "OSQM"), ("UW", "ODUC"),
    # added 2026-09-19 (coordinator): southern Vancouver Island PQ broadbands within 45 km of
    # 36-66 Bostock families each, plus two Olympic UW broadbands.  All six PROBED first
    # (r_probe_delivery.py / r_probe_extent.py) - metadata open is not delivery.
    ("PQ", "LHLYB"), ("PQ", "JRDRB"), ("PQ", "ALBHB"), ("PQ", "DAOB"),
    ("UW", "OHOH"), ("UW", "ODEER"),
]
ROSTER_Z = [("UW", "STW"), ("UW", "UWFH")]          # Z-only helpers, EHZ
ROSTER = ROSTER_3C + ROSTER_Z
ZONLY = {s for _, s in ROSTER_Z}
BAND_PREF = ["HH", "BH", "EH"]                      # seismometers only
BAND_PREF_Z = ["EH", "HH", "BH"]                    # helpers: EHZ is the prescribed stream

META0 = date(2014, 1, 1)                            # metadata query window
META1 = date(2026, 9, 19)
DL0 = date(2019, 1, 1)                              # download window (inclusive)
DL1 = date(2026, 9, 15)

# PNSN northern-episode box
BOX = dict(lat0=48.0, lat1=49.2, lon0=-124.6, lon1=-122.8)
SSE_MIN_PER_DAY = 50        # a "tremor day"
SSE_MAX_GAP = 20            # days; cluster tremor days with gaps <= this
SSE_MIN_DAYS = 10           # minimum span of an episode (core, in days)
SSE_PAD = 5                 # days padded on each side


# ---------------------------------------------------------------- metadata
def _fetch_inventory(net, stas, log):
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime
    last = None
    for prov in ["IRIS"]:
        try:
            c = Client(prov, timeout=180)
            inv = c.get_stations(network=net, station=",".join(stas), level="channel",
                                 starttime=UTCDateTime(str(META0)),
                                 endtime=UTCDateTime(str(META1)))
            log.append(f"{net}: metadata from {prov} ({sum(len(n) for n in inv)} stations)")
            return inv
        except Exception as e:
            last = e
            log.append(f"{net}: metadata FAILED at {prov}: {type(e).__name__}: {e}")
    raise RuntimeError(f"no metadata for {net}: {last}")


def _as_date(t, default):
    if t is None:
        return default
    return t.datetime.date() if hasattr(t, "datetime") else t


def station_epochs(inv, log):
    """One row per contiguous (loc, band) interval that is the preferred band at that time
    and carries the required number of components (3, or 1 for the Z-only helpers)."""
    rows = []
    for netobj in inv:
        for sta in netobj:
            zonly = sta.code in ZONLY
            pref = BAND_PREF_Z if zonly else BAND_PREF
            need = 1 if zonly else 3
            per = {}
            for ch in sta:
                band = ch.code[:2]
                if band not in pref:
                    continue
                if zonly and ch.code[-1] != "Z":
                    continue
                s = _as_date(ch.start_date, META0)
                e = _as_date(ch.end_date, META1)
                per.setdefault((ch.location_code, band), []).append(
                    (s, e, ch.code, float(ch.sample_rate),
                     float(ch.azimuth if ch.azimuth is not None else 0.0),
                     float(ch.dip if ch.dip is not None else 0.0)))
            if not per:
                log.append(f"  {netobj.code}.{sta.code}: NO usable band in {pref}")
                continue
            bounds = {META0, META1}
            for v in per.values():
                for s, e, _, _, _, _ in v:
                    bounds.add(max(s, META0))
                    bounds.add(min(e, META1))
            bnd = sorted(b for b in bounds if META0 <= b <= META1)
            segs = []
            for a, b in zip(bnd[:-1], bnd[1:]):
                if a >= b:
                    continue
                mid = a + (b - a) / 2
                avail = {}
                for (loc, band), v in per.items():
                    live = [(c, fs, az, dip) for s, e, c, fs, az, dip in v if s <= mid < e]
                    comps = {c[-1] for c, _, _, _ in live}
                    if len(comps) >= need:
                        # one entry per channel code, deterministic order Z, N|1, E|2
                        byc = {}
                        for c, fs, az, dip in live:
                            byc[c] = (fs, az, dip)
                        avail[(loc, band)] = byc
                if not avail:
                    continue
                key = sorted(avail, key=lambda k: (pref.index(k[1]), k[0]))[0]
                byc = avail[key]
                fs = max(v[0] for v in byc.values())
                segs.append((a, b, key[0], key[1], fs, len(byc), byc))
            for a, b, loc, band, fs, ncomp, byc in segs:
                azd = ";".join(f"{c}:{byc[c][1]:.1f}/{byc[c][2]:.1f}"
                               for c in sorted(byc, key=lambda c: "ZNE12".find(c[-1])))
                if rows and rows[-1]["sta"] == sta.code and rows[-1]["net"] == netobj.code \
                        and rows[-1]["band"] == band and rows[-1]["loc"] == loc \
                        and rows[-1]["_end"] == a and rows[-1]["chan_az_dip"] == azd:
                    rows[-1]["_end"] = b
                    continue
                rows.append(dict(net=netobj.code, sta=sta.code, loc=loc, band=band, fs=fs,
                                 _start=a, _end=b, n_comp=ncomp, lat=sta.latitude,
                                 lon=sta.longitude, elev=sta.elevation,
                                 chan_az_dip=azd, mode="Z" if zonly else "3C"))
            log.append(f"  {netobj.code}.{sta.code}: " + "; ".join(
                f"{r['band']}{r['fs']:.0f}Hz loc'{r['loc']}' n={r['n_comp']} "
                f"{r['_start']}..{r['_end']}"
                for r in rows if r["sta"] == sta.code and r["net"] == netobj.code))
    return rows


EP_COLS = ["net", "sta", "band", "fs", "loc", "start", "end", "lat", "lon", "elev",
           "n_comp", "chan_az_dip", "mode", "first_data_month", "last_data_month",
           "months_with_data", "months_probed"]


def data_extent():
    """Empirical month-by-month delivery from r_probe_extent.py, if it has been run.

    StationXML epochs lie (PQ.JRDRB opens 2021-04-01 and the archive is empty until 2024-10),
    so the epoch table carries the MEASURED extent alongside the declared one.  The queue is
    still built from the declared epoch - delivery is recorded empirically in the manifest,
    which is the GOBB rule - but downstream code can trust these columns.
    """
    p = os.path.join(RESULTS, "probe_month_extent.csv")
    out = {}
    if not os.path.exists(p):
        return out
    with open(p) as f:
        for r in csv.DictReader(f):
            k = (r["net"], r["sta"])
            o = out.setdefault(k, dict(first="", last="", n=0, tot=0))
            o["tot"] += 1
            if int(r["days_with_data"]) > 0:
                o["n"] += 1
                o["first"] = o["first"] or r["month"]
                o["last"] = r["month"]
    return out


def write_epochs(rows):
    ext = data_extent()
    p = os.path.join(RESULTS, "station_epochs.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(EP_COLS)
        for r in sorted(rows, key=lambda r: (r["net"], r["sta"], r["_start"])):
            e = ext.get((r["net"], r["sta"]), {})
            w.writerow([r["net"], r["sta"], r["band"], f"{r['fs']:g}", r["loc"],
                        r["_start"], r["_end"], f"{r['lat']:.5f}", f"{r['lon']:.5f}",
                        f"{r['elev']:.1f}", r["n_comp"], r["chan_az_dip"], r["mode"],
                        e.get("first", ""), e.get("last", ""), e.get("n", ""),
                        e.get("tot", "")])
    return p


def read_epochs():
    rows = []
    with open(os.path.join(RESULTS, "station_epochs.csv")) as f:
        for r in csv.DictReader(f):
            r["_start"] = date.fromisoformat(r["start"])
            r["_end"] = date.fromisoformat(r["end"])
            r["fs"] = float(r["fs"])
            rows.append(r)
    return rows


# ---------------------------------------------------------------- SSE windows
def tremor_days(log):
    """{date: n} epicentres inside the northern box, whole catalog."""
    import pandas as pd
    d = pd.read_csv(CATALOG, usecols=["time", "lat", "lon"])
    d["time"] = pd.to_datetime(d["time"])
    b = d[(d.lat >= BOX["lat0"]) & (d.lat <= BOX["lat1"]) &
          (d.lon >= BOX["lon0"]) & (d.lon <= BOX["lon1"])]
    cnt = b.groupby(b.time.dt.date).size()
    log.append(f"PNSN catalog {CATALOG}: {len(d)} rows, {len(b)} inside the box "
               f"({BOX['lat0']}-{BOX['lat1']}N, {BOX['lon0']}..{BOX['lon1']}E); "
               f"catalog spans {d.time.min().date()}..{d.time.max().date()}")
    return {k: int(v) for k, v in cnt.items()}


def sse_windows(cnt, log):
    days = sorted(d for d, n in cnt.items() if n >= SSE_MIN_PER_DAY and DL0 <= d <= DL1)
    clusters, cur = [], []
    for d in days:
        if cur and (d - cur[-1]).days > SSE_MAX_GAP:
            clusters.append(cur)
            cur = []
        cur.append(d)
    if cur:
        clusters.append(cur)
    out = []
    for c in clusters:
        span = (c[-1] - c[0]).days + 1
        if span < SSE_MIN_DAYS:
            continue
        s = max(DL0, c[0] - timedelta(days=SSE_PAD))
        e = min(DL1, c[-1] + timedelta(days=SSE_PAD))
        ntr = sum(cnt.get(c[0] + timedelta(days=k), 0) for k in range(span))
        out.append(dict(id=f"S{len(out)+1:02d}", core_start=c[0].isoformat(),
                        core_end=c[-1].isoformat(), start=s.isoformat(), end=e.isoformat(),
                        n_tremor_days=len(c), span_days=span,
                        padded_days=(e - s).days + 1, n_tremor=ntr))
    p = os.path.join(RESULTS, "sse_windows.csv")
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "core_start", "core_end", "start", "end",
                                          "n_tremor_days", "span_days", "padded_days",
                                          "n_tremor"])
        w.writeheader()
        w.writerows(out)
    log.append(f"SSE windows -> {p}: {len(out)} episodes, "
               f"{sum(r['padded_days'] for r in out)} padded days")
    for r in out:
        log.append(f"  {r['id']} {r['start']}..{r['end']} (core {r['core_start']}.."
                   f"{r['core_end']}, {r['span_days']} d, {r['n_tremor']} epicentres)")
    return out


def quiet_months(cnt, log):
    """One quiet calendar month per year 2019..2026 - fewest epicentres in the box."""
    per = {}
    d = DL0
    while d <= DL1:
        per.setdefault(d.strftime("%Y-%m"), 0)
        per[d.strftime("%Y-%m")] += cnt.get(d, 0)
        d += timedelta(days=1)
    chosen = []
    for y in range(DL0.year, DL1.year + 1):
        cand = sorted(((v, k) for k, v in per.items() if k.startswith(str(y))))
        if not cand:
            continue
        chosen.append(cand[0][1])
        log.append(f"  quiet month {y}: {cand[0][1]} ({cand[0][0]} epicentres in the box)")
    return set(chosen)


# ---------------------------------------------------------------- queue
def build_queue(rows, sse, quiet, log):
    byst = {}
    for r in rows:
        byst.setdefault((r["net"], r["sta"]), []).append(r)
    p0 = set()
    for w in sse:
        d = date.fromisoformat(w["start"])
        e = date.fromisoformat(w["end"])
        while d <= e:
            p0.add(d)
            d += timedelta(days=1)

    def prio(d):
        if d in p0:
            return "P0"
        if d.strftime("%Y-%m") in quiet:
            return "P1"
        return "P2"

    q = []
    d = DL0
    while d <= DL1:
        pr = prio(d)
        for net, sta in ROSTER:
            eps = [e for e in byst.get((net, sta), []) if e["_start"] <= d < e["_end"]]
            if not eps:
                q.append(dict(net=net, sta=sta, band="", loc="", fs="", date=d.isoformat(),
                              doy=f"{d.timetuple().tm_yday:03d}", priority=pr,
                              status="no_epoch", path=""))
                continue
            e = eps[0]
            path = os.path.join(TRACES, f"{net}.{sta}", f"{d.year}",
                                f"{d.timetuple().tm_yday:03d}.mseed")
            q.append(dict(net=net, sta=sta, band=e["band"], loc=e["loc"],
                          fs=f"{float(e['fs']):g}", date=d.isoformat(),
                          doy=f"{d.timetuple().tm_yday:03d}", priority=pr,
                          status="pending", path=""))
        d += timedelta(days=1)
    order = {"P0": 0, "P1": 1, "P2": 2}
    q.sort(key=lambda r: (order[r["priority"]], r["date"], r["net"], r["sta"]))
    p = os.path.join(RESULTS, "download_queue.csv")
    tmp = p + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["net", "sta", "band", "loc", "fs", "date", "doy",
                                          "priority", "status", "path"])
        w.writeheader()
        w.writerows(q)
    os.replace(tmp, p)
    summ = {}
    for r in q:
        summ.setdefault(r["priority"], {}).setdefault(r["status"], 0)
        summ[r["priority"]][r["status"]] += 1
    log.append("queue -> " + p)
    for pr in sorted(summ):
        log.append(f"  {pr}: " + ", ".join(f"{k}={v}" for k, v in sorted(summ[pr].items())))
    bysta = {}
    for r in q:
        if r["status"] == "pending":
            bysta.setdefault((r["net"], r["sta"]), {}).setdefault(r["priority"], 0)
            bysta[(r["net"], r["sta"])][r["priority"]] += 1
    for k in sorted(bysta):
        log.append(f"    {k[0]}.{k[1]}: " + ", ".join(
            f"{p2}={n}" for p2, n in sorted(bysta[k].items())))
    return q, summ


MANIFEST_COLS = ["net", "sta", "band", "date", "doy", "status", "bytes", "n_channels",
                 "fs", "coverage_frac", "n_frag", "updated_utc", "note"]


def seed_manifest(q, log):
    """The worker only walks `pending` rows, so the out-of-epoch days are recorded here."""
    p = os.path.join(RESULTS, "download_manifest.csv")
    have = set()
    if os.path.exists(p):
        with open(p) as f:
            for r in csv.DictReader(f):
                have.add((r["sta"], r["date"]))
    else:
        with open(p, "w", newline="") as f:
            csv.writer(f).writerow(MANIFEST_COLS)
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    n = 0
    with open(p, "a", newline="") as f:
        w = csv.writer(f)
        for r in q:
            if r["status"] != "no_epoch" or (r["sta"], r["date"]) in have:
                continue
            w.writerow([r["net"], r["sta"], "", r["date"], r["doy"], "no_epoch",
                        0, 0, "", "", "", now, "outside station epoch"])
            n += 1
    log.append(f"manifest seeded with {n} no_epoch rows -> {p}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-fdsn", action="store_true")
    a = ap.parse_args()
    os.makedirs(RESULTS, exist_ok=True)
    log = [f"P3D download plan built {datetime.now(timezone.utc).isoformat(timespec='seconds')}"]
    if a.no_fdsn:
        rows = read_epochs()
        log.append(f"epochs reused ({len(rows)} rows)")
    else:
        rows = []
        bynet = {}
        for net, sta in ROSTER:
            bynet.setdefault(net, []).append(sta)
        for net, stas in bynet.items():
            rows += station_epochs(_fetch_inventory(net, stas, log), log)
        log.append("epochs -> " + write_epochs(rows))
        missing = {(n, s) for n, s in ROSTER} - {(r["net"], r["sta"]) for r in rows}
        if missing:
            log.append(f"WARNING no epoch found for {sorted(missing)}")
    cnt = tremor_days(log)
    sse = sse_windows(cnt, log)
    quiet = quiet_months(cnt, log)
    q, summ = build_queue(rows, sse, quiet, log)
    seed_manifest(q, log)
    with open(os.path.join(RESULTS, "download_plan_log.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    print("\n".join(log))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
