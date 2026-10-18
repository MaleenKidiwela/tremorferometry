#!/usr/bin/env python
"""Track B step 1 - station/band epoch table + ordered download queue.

Queries FDSN station metadata (level=channel) for the lfe_extension roster, picks ONE
seismometer band per epoch (HH > BH > EH; accelerometers HN/EN and LH/SN ignored), writes

  results/station_epochs.csv    net,sta,loc,band,fs,start,end,lat,lon,n_comp
  results/epoch_crosscheck.csv  comparison against the read-only codex-LFE triplet table
  results/download_queue.csv    net,sta,band,loc,fs,date,doy,priority,status,path
  results/download_manifest.csv seeded with the present_elsewhere rows (runner appends)

Priorities
  P0 2010-08-10..2010-09-13  PB B001..B010 only
  P1 2011-08-05..2011-09-25  all stations
  P2 2013-08-30..2013-10-28  all stations
  P3 2014-01-01..2026-09-15  all stations, oldest month first

Days whose file already exists in this tree or in ../lfe_catalog_repair/traces (same band)
are marked present_elsewhere; days outside any epoch are marked no_epoch.

Usage: python scripts/download_plan.py [--no-fdsn]   (--no-fdsn reuses station_epochs.csv)
"""
import os
import csv
import glob
import json
import argparse
from datetime import date, datetime, timedelta, timezone

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(HERE, "results")
TRACES = os.path.join(HERE, "traces")
# The 2008/2009 windows do not overlap ANY tree on disk (the earliest existing day is
# 2010-08-10), so there is nothing to index and the foreign trees are left alone.
OTHER_TREES = []
CODEX_EPOCHS = os.path.abspath(os.path.join(
    HERE, "..", "codex-LFE", "results", "full_history_work_plan_v1",
    "nonoverlapping_triplet_epochs.csv"))

ROSTER = [
    ("CN", "PGC"), ("CN", "VGZ"), ("CN", "SNB"), ("CN", "NLLB"), ("CN", "LZB"),
    ("CN", "YOUB"), ("CN", "PFB"),
    ("UW", "LRIV"),
    ("C8", "SOKB"),
    ("PB", "B001"), ("PB", "B003"), ("PB", "B004"), ("PB", "B005"), ("PB", "B006"),
    ("PB", "B007"), ("PB", "B009"), ("PB", "B010"), ("PB", "B011"), ("PB", "B926"),
]
BAND_PREF = ["HH", "BH", "EH"]          # seismometers only, in preference order
# lfe_picker P0-O item 1: ONLY the two pre-2010 ETS training windows.  Everything from 2010
# onward is already on disk in ../lfe_extension/traces and is NEVER re-fetched here.
WINDOW0 = date(2008, 1, 1)
WINDOW1 = date(2010, 1, 1)              # exclusive for metadata query
P0_STAS = set()
PRIORITIES = [
    ("O2008", date(2008, 5, 1), date(2008, 6, 10)),
    ("O2009", date(2009, 5, 1), date(2009, 6, 15)),
]
NRCAN = "https://earthquakescanada.nrcan.gc.ca/fdsnws"


# ---------------------------------------------------------------- metadata
def _fetch_inventory(net, stas, log):
    """level=channel inventory for one network; IRIS first, NRCan fallback for CN."""
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime
    provs = ["IRIS"] + ([NRCAN] if net == "CN" else [])
    last = None
    for prov in provs:
        try:
            c = Client(prov, timeout=180)
            inv = c.get_stations(network=net, station=",".join(stas), level="channel",
                                 starttime=UTCDateTime(str(WINDOW0)),
                                 endtime=UTCDateTime(str(WINDOW1)))
            log.append(f"{net}: metadata from {prov} ({sum(len(n) for n in inv)} stations)")
            return inv
        except Exception as e:                                   # pragma: no cover
            last = e
            log.append(f"{net}: metadata FAILED at {prov}: {type(e).__name__}: {e}")
    raise RuntimeError(f"no metadata for {net}: {last}")


def _as_date(t, default):
    if t is None:
        return default
    d = t.datetime.date() if hasattr(t, "datetime") else t
    return d


def station_epochs(inv, log):
    """-> list of dicts: one row per contiguous (loc, band) interval that is the preferred
    band available at that time and has 3 components."""
    rows = []
    for netobj in inv:
        for sta in netobj:
            # gather channel epochs per (loc, band), seismometers only
            per = {}
            for ch in sta:
                band = ch.code[:2]
                if band not in BAND_PREF:
                    continue
                s = _as_date(ch.start_date, WINDOW0)
                e = _as_date(ch.end_date, WINDOW1)
                per.setdefault((ch.location_code, band), []).append(
                    (s, e, ch.code[-1], float(ch.sample_rate)))
            if not per:
                continue
            # boundary times
            bounds = {WINDOW0, WINDOW1}
            for v in per.values():
                for s, e, _, _ in v:
                    bounds.add(max(s, WINDOW0))
                    bounds.add(min(e, WINDOW1))
            bnd = sorted(b for b in bounds if WINDOW0 <= b <= WINDOW1)
            segs = []
            for a, b in zip(bnd[:-1], bnd[1:]):
                if a >= b:
                    continue
                mid = a + (b - a) / 2
                avail = {}
                for (loc, band), v in per.items():
                    comps, fss = set(), set()
                    for s, e, comp, fs in v:
                        if s <= mid < e:
                            comps.add(comp)
                            fss.add(fs)
                    if len(comps) >= 3:
                        avail[(loc, band)] = (len(comps), max(fss))
                if not avail:
                    continue
                # preference: band order, then lowest location code
                key = sorted(avail, key=lambda k: (BAND_PREF.index(k[1]), k[0]))[0]
                segs.append((a, b, key[0], key[1], avail[key][1], avail[key][0]))
            # merge contiguous identical segments
            for seg in segs:
                if rows and rows[-1]["sta"] == sta.code and rows[-1]["net"] == netobj.code \
                        and rows[-1]["band"] == seg[3] and rows[-1]["loc"] == seg[2] \
                        and rows[-1]["_end"] == seg[0]:
                    rows[-1]["_end"] = seg[1]
                    continue
                rows.append(dict(net=netobj.code, sta=sta.code, loc=seg[2], band=seg[3],
                                 fs=seg[4], _start=seg[0], _end=seg[1], n_comp=seg[5],
                                 lat=sta.latitude, lon=sta.longitude))
            log.append(f"  {netobj.code}.{sta.code}: " + "; ".join(
                f"{r['band']}{r['fs']:.0f}Hz {r['_start']}..{r['_end']}"
                for r in rows if r["sta"] == sta.code and r["net"] == netobj.code))
    return rows


def write_epochs(rows):
    p = os.path.join(RESULTS, "station_epochs.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["net", "sta", "band", "fs", "loc", "start", "end", "lat", "lon", "n_comp"])
        for r in sorted(rows, key=lambda r: (r["net"], r["sta"], r["_start"])):
            w.writerow([r["net"], r["sta"], r["band"], f"{r['fs']:g}", r["loc"],
                        r["_start"], r["_end"], f"{r['lat']:.5f}", f"{r['lon']:.5f}",
                        r["n_comp"]])
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


# ---------------------------------------------------------------- cross-check
def crosscheck(rows, log):
    """Interval-by-interval comparison with the codex-LFE triplet table, clipped to our
    metadata window (the codex rows start long before 2010-08-01)."""
    out = []
    if not os.path.exists(CODEX_EPOCHS):
        log.append(f"cross-check SKIPPED: {CODEX_EPOCHS} missing")
        return out
    ours = {}
    for r in rows:
        ours.setdefault((r["net"], r["sta"]), []).append(r)
    with open(CODEX_EPOCHS) as f:
        for c in csv.DictReader(f):
            key = (c["network"], c["station"])
            if key not in ours:
                continue
            cs = max(date.fromisoformat(c["start_utc"][:10]), WINDOW0)
            ce = min(date.fromisoformat(c["end_utc_exclusive"][:10]), WINDOW1)
            if cs >= ce:
                continue                      # codex epoch entirely before our window
            cb = sorted({ch.split(".")[-1][:2] for ch in c["channels"].split(",")})
            mine = [r for r in ours[key] if r["_start"] < ce and r["_end"] > cs]
            my_bands = sorted({r["band"] for r in mine})
            # overlapping-interval band agreement
            ok = bool(mine) and all(
                (set(cb) == {r["band"]}) for r in mine
                if max(r["_start"], cs) < min(r["_end"], ce))
            out.append(dict(net=key[0], sta=key[1], codex_band=",".join(cb),
                            codex_start=cs.isoformat(), codex_end=ce.isoformat(),
                            ours_bands=",".join(my_bands),
                            ours_start=min((r["_start"] for r in mine), default=""),
                            ours_end=max((r["_end"] for r in mine), default=""),
                            agree="yes" if ok else "no",
                            codex_channels=c["channels"]))
    p = os.path.join(RESULTS, "epoch_crosscheck.csv")
    cols = ["net", "sta", "codex_band", "codex_start", "codex_end", "ours_bands",
            "ours_start", "ours_end", "agree", "codex_channels"]
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(out)
    bad = [o for o in out if o["agree"] == "no"]
    log.append(f"cross-check: {len(out)} codex epochs overlap our window, {len(bad)} disagree "
               f"-> {p}")
    for o in bad:
        log.append(f"    DISAGREE {o['net']}.{o['sta']} codex {o['codex_band']} "
                   f"{o['codex_start']}..{o['codex_end']} vs ours {o['ours_bands']}")
    missing = sorted({(r['net'], r['sta']) for r in rows} -
                     {(o['net'], o['sta']) for o in out})
    log.append(f"    not in codex table (no cross-check possible): {missing}")
    return out


# ---------------------------------------------------------------- existing files
def index_existing(log):
    """{(net.sta, year, doy): (path, band)} over the READ-ONLY foreign trees only.

    Days already in our own traces/ are NOT marked present_elsewhere: they either already
    carry a manifest row (the runner skips them) or fetch_day verifies the file on disk and
    writes a real `data` row with bytes/coverage - which is what Track C wants to read.
    """
    from obspy import read
    idx = {}
    for root in OTHER_TREES:
        if not os.path.isdir(root):
            continue
        n = 0
        for p in glob.glob(os.path.join(root, "*", "*", "*.mseed")):
            parts = p.split(os.sep)
            stak, yr, fn = parts[-3], parts[-2], parts[-1]
            if "." not in stak or not yr.isdigit():
                continue
            doy = fn.split(".")[0]
            if not doy.isdigit():
                continue
            key = (stak, int(yr), int(doy))
            if key in idx:
                continue
            idx[key] = (p, None)
            n += 1
        log.append(f"existing .mseed under {root}: {n}")
    # band per (station, year) - read one file per station-year (bands change rarely,
    # and never inside a year for this roster except at a documented switch date)
    byyear = {}
    for (stak, yr, doy), (p, _) in idx.items():
        byyear.setdefault((stak, yr), []).append((doy, p))
    bands = {}
    for k, v in byyear.items():
        v.sort()
        b = set()
        for _, p in [v[0], v[-1]]:
            try:
                st = read(p, headonly=True)
                b |= {t.stats.channel[:2] for t in st}
            except Exception:
                pass
        bands[k] = b
    for key in idx:
        idx[key] = (idx[key][0], bands.get((key[0], key[1]), set()))
    return idx


# ---------------------------------------------------------------- queue
def build_queue(rows, existing, log):
    byst = {}
    for r in rows:
        byst.setdefault((r["net"], r["sta"]), []).append(r)
    q = []
    seen = set()
    for prio, d0, d1 in PRIORITIES:
        d = d0
        while d <= d1:
            for net, sta in ROSTER:
                if (sta, d) in seen:
                    continue                 # already queued at a higher priority
                seen.add((sta, d))
                if prio == "P0" and not (net == "PB" and sta in P0_STAS):
                    continue
                eps = [e for e in byst.get((net, sta), []) if e["_start"] <= d < e["_end"]]
                if not eps:
                    q.append(dict(net=net, sta=sta, band="", loc="", fs="", date=d.isoformat(),
                                  doy=f"{d.timetuple().tm_yday:03d}", priority=prio,
                                  status="no_epoch", path=""))
                    continue
                e = eps[0]
                key = (f"{net}.{sta}", d.year, d.timetuple().tm_yday)
                st, path = "pending", ""
                if key in existing:
                    p, bands = existing[key]
                    if not bands or e["band"] in bands:
                        st, path = "present_elsewhere", p
                q.append(dict(net=net, sta=sta, band=e["band"], loc=e["loc"],
                              fs=f"{e['fs']:g}", date=d.isoformat(),
                              doy=f"{d.timetuple().tm_yday:03d}", priority=prio,
                              status=st, path=path))
            d += timedelta(days=1)
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
    for prio in sorted(summ):
        log.append(f"  {prio}: " + ", ".join(f"{k}={v}" for k, v in sorted(summ[prio].items())))
    return q, summ


MANIFEST_COLS = ["net", "sta", "band", "date", "doy", "status", "bytes", "n_channels",
                 "fs", "coverage_frac", "n_frag", "updated_utc", "note"]


def seed_manifest(q, log):
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
            if r["status"] != "present_elsewhere" or (r["sta"], r["date"]) in have:
                continue
            w.writerow([r["net"], r["sta"], r["band"], r["date"], r["doy"],
                        "present_elsewhere", "", "", r["fs"], "", "", now, r["path"]])
            n += 1
    log.append(f"manifest seeded with {n} present_elsewhere rows -> {p}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-fdsn", action="store_true",
                    help="reuse results/station_epochs.csv instead of querying FDSN")
    a = ap.parse_args()
    os.makedirs(RESULTS, exist_ok=True)
    log = []
    if a.no_fdsn:
        rows = read_epochs()
        log.append(f"epochs reused from results/station_epochs.csv ({len(rows)} rows)")
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
    crosscheck(rows, log)
    existing = index_existing(log)
    q, summ = build_queue(rows, existing, log)
    seed_manifest(q, log)
    with open(os.path.join(RESULTS, "download_plan_log.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    print("\n".join(log))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
