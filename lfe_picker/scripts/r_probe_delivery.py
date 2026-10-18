#!/usr/bin/env python
"""P3D step 1b - PROBE whether a candidate station ACTUALLY delivers waveforms.

Metadata "open" is not delivery (the CN.GOBB lesson: wide-open StationXML, whole SSE episodes
absent from the archive).  For every station this fetches level=channel metadata, then asks each
provider in turn for THREE SSE days per station-year and records exactly what came back.

Providers tried in order: IRIS, then the NRCan FDSN endpoint.  NOTE (2026-09-19): NRCan is
reachable at https://www.earthquakescanada.nrcan.gc.ca/fdsnws - the bare
earthquakescanada.nrcan.gc.ca host answers 301 and obspy does not follow the redirect, which is
why NRCan was recorded as "unreachable from this pod" in the earlier work.

  python scripts/r_probe_delivery.py --stations PQ.LHLYB,PQ.JRDRB,... [--per-year 3]

Writes results/p3d/probe_delivery.csv (one row per station-day probed) and
       results/p3d/probe_summary.csv (one row per station) and prints the verdict table.
Probe waveforms are written into traces/ like any other day, so nothing is fetched twice.
"""
import os
import csv
import sys
import argparse
from datetime import date, datetime, timezone

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))
RESULTS = os.path.join(HERE, "results", "p3d")
TRACES = os.path.join(HERE, "traces")
IRIS = "IRIS"
NRCAN = "https://www.earthquakescanada.nrcan.gc.ca/fdsnws"
PROVIDERS = [("IRIS", IRIS), ("NRCAN", NRCAN)]


def nrcan_query(net, sta, band, d, hours=1):
    """Ask NRCan dataselect over PLAIN HTTP and report what it says.

    obspy cannot be pointed at NRCan (its service discovery raises FDSNNoServiceException even
    though /fdsnws/station/1/version and /fdsnws/dataselect/1/query both answer 200), and
    fetch_day would short-circuit on the `.nd` marker the IRIS attempt just wrote - so the
    fallback has to be an independent request.  Returns (status, http_code, n_bytes).
    """
    import urllib.request
    import urllib.error
    url = (f"{NRCAN}/dataselect/1/query?net={net}&sta={sta}&loc=*&cha={band}?"
           f"&start={d.isoformat()}T00:00:00&end={d.isoformat()}T{hours:02d}:00:00")
    try:
        with urllib.request.urlopen(url, timeout=180) as r:
            body = r.read()
        return (("data" if len(body) > 512 else "nodata"), r.status, len(body))
    except urllib.error.HTTPError as e:
        return (("nodata" if e.code in (204, 404) else "error"), e.code, 0)
    except Exception as e:                                          # network / DNS / TLS
        return ("error", f"{type(e).__name__}", 0)
BAND_PREF = ["HH", "BH", "EH"]
META0, META1 = date(2019, 1, 1), date(2026, 9, 19)


def fetch_meta(net, stas, log):
    """-> {(net,sta): dict(lat, lon, elev, band, fs, loc, start, end, chan_az_dip, prov)}"""
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime
    out = {}
    for pname, purl in PROVIDERS:
        todo = [s for s in stas if (net, s) not in out]
        if not todo:
            break
        try:
            inv = Client(purl, timeout=180).get_stations(
                network=net, station=",".join(todo), level="channel",
                starttime=UTCDateTime(str(META0)), endtime=UTCDateTime(str(META1)))
        except Exception as e:
            log.append(f"  metadata {net} @ {pname}: FAIL {type(e).__name__}: {str(e)[:120]}")
            continue
        for n in inv:
            for s in n:
                per = {}
                for ch in s:
                    if ch.code[:2] not in BAND_PREF:
                        continue
                    per.setdefault((ch.location_code, ch.code[:2]), {})[ch.code] = (
                        float(ch.sample_rate),
                        float(ch.azimuth if ch.azimuth is not None else 0.0),
                        float(ch.dip if ch.dip is not None else 0.0),
                        (ch.start_date.datetime.date() if ch.start_date else META0),
                        (ch.end_date.datetime.date() if ch.end_date else META1))
                cand = {k: v for k, v in per.items() if len(v) >= 3}
                if not cand:
                    continue
                key = sorted(cand, key=lambda k: (BAND_PREF.index(k[1]), k[0]))[0]
                byc = cand[key]
                st = min(v[3] for v in byc.values())
                en = max(v[4] for v in byc.values())
                out[(n.code, s.code)] = dict(
                    net=n.code, sta=s.code, lat=s.latitude, lon=s.longitude,
                    elev=s.elevation, loc=key[0], band=key[1],
                    fs=max(v[0] for v in byc.values()),
                    start=max(st, META0), end=min(en, META1),
                    n_comp=len(byc), meta_prov=pname,
                    chan_az_dip=";".join(
                        f"{c}:{byc[c][1]:.1f}/{byc[c][2]:.1f}"
                        for c in sorted(byc, key=lambda c: "ZNE12".find(c[-1]))))
                log.append(f"  metadata {n.code}.{s.code} @ {pname}: {key[1]} "
                           f"{max(v[0] for v in byc.values()):g} Hz loc'{key[0]}' "
                           f"{max(st, META0)}..{min(en, META1)}")
    return out


def sse_days():
    ws = []
    with open(os.path.join(RESULTS, "sse_windows.csv")) as f:
        for r in csv.DictReader(f):
            ws.append((r["id"], date.fromisoformat(r["core_start"]),
                       date.fromisoformat(r["core_end"])))
    return ws


def pick_days(meta, wins, per_year):
    """3 SSE days per station-year, spread over the episodes that fall in that year."""
    out = []
    for y in range(meta["start"].year, min(meta["end"].year, 2026) + 1):
        cand = []
        for wid, a, b in wins:
            if a.year != y:
                continue
            span = (b - a).days
            for frac in (0.25, 0.5, 0.75):
                from datetime import timedelta
                d = a + timedelta(days=int(round(frac * span)))
                if meta["start"] <= d < meta["end"]:
                    cand.append((wid, d))
        seen, keep = set(), []
        # spread across DIFFERENT episodes first
        for wid, d in cand:
            if wid in seen:
                continue
            seen.add(wid)
            keep.append((wid, d))
        for wid, d in cand:
            if len(keep) >= per_year:
                break
            if (wid, d) not in keep:
                keep.append((wid, d))
        out += [(y, wid, d) for wid, d in keep[:per_year]]
    return out


def probe(meta, days, log):
    from r_download_day import fetch_day
    rows = []
    for y, wid, d in days:
        got = None
        rec = fetch_day(meta["net"], meta["sta"], meta["band"], d, meta["loc"], TRACES, IRIS)
        rec.update(provider="IRIS", year=y, sse_id=wid)
        rows.append(dict(rec))
        if rec["status"] == "data":
            got = "IRIS"
        else:
            # IRIS empty -> ask NRCan INDEPENDENTLY (see nrcan_query for why not fetch_day)
            st, code, nb = nrcan_query(meta["net"], meta["sta"], meta["band"], d)
            rows.append(dict(net=meta["net"], sta=meta["sta"], band=meta["band"],
                             date=d.isoformat(), doy=f"{d.timetuple().tm_yday:03d}",
                             year=y, sse_id=wid, provider="NRCAN", status=st, bytes=nb,
                             n_channels="", fs="", coverage_frac="", n_frag="",
                             note=f"http {code} (1 h probe)"))
            if st == "data":
                got = "NRCAN"
        log.append(f"    {meta['net']}.{meta['sta']} {d} ({wid}): "
                   f"{'DATA from ' + got if got else 'no data from any provider'}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", required=True, help="NET.STA,NET.STA,...")
    ap.add_argument("--per-year", type=int, default=3)
    a = ap.parse_args()
    os.makedirs(RESULTS, exist_ok=True)
    log = [f"probe {datetime.now(timezone.utc).isoformat(timespec='seconds')}"]
    bynet = {}
    for k in a.stations.split(","):
        n, s = k.split(".")
        bynet.setdefault(n, []).append(s)
    meta = {}
    for n, ss in bynet.items():
        meta.update(fetch_meta(n, ss, log))
    allrows, summ = [], []
    wins = sse_days()
    for k in a.stations.split(","):
        n, s = k.split(".")
        m = meta.get((n, s))
        if m is None:
            log.append(f"  {k}: NO 3-C metadata at any provider")
            summ.append(dict(net=n, sta=s, meta="none", delivers="NO_METADATA"))
            continue
        days = pick_days(m, wins, a.per_year)
        log.append(f"  {k}: probing {len(days)} SSE days "
                   f"({m['start']}..{m['end']}, {m['band']} {m['fs']:g} Hz)")
        rows = probe(m, days, log)
        allrows += rows
        ok = [r for r in rows if r["status"] == "data"]
        provs = sorted({r["provider"] for r in ok})
        cov = [float(r["coverage_frac"]) for r in ok if r["coverage_frac"]]
        yrs = sorted({r["year"] for r in rows})
        yrs_ok = sorted({r["year"] for r in ok})
        summ.append(dict(
            net=n, sta=s, lat=round(m["lat"], 4), lon=round(m["lon"], 4),
            band=m["band"], fs=f"{m['fs']:g}", loc=m["loc"], n_comp=m["n_comp"],
            epoch_start=m["start"], epoch_end=m["end"], meta_prov=m["meta_prov"],
            chan_az_dip=m["chan_az_dip"],
            probed=len(days), delivered=len(ok),
            delivers=("YES" if len(ok) == len(days) else
                      "PARTIAL" if ok else "NO"),
            from_provider=",".join(provs) if provs else "",
            years_probed=",".join(str(x) for x in yrs),
            years_no_data=",".join(str(x) for x in yrs if x not in yrs_ok),
            cov_min=round(min(cov), 4) if cov else "",
            mb_per_day=(round(sum(int(r["bytes"]) for r in ok) / 1e6 / len(ok), 1)
                        if ok else "")))
    if allrows:
        cols = ["net", "sta", "band", "date", "doy", "year", "sse_id", "provider", "status",
                "bytes", "n_channels", "fs", "coverage_frac", "n_frag", "note"]
        with open(os.path.join(RESULTS, "probe_delivery.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(allrows)
    with open(os.path.join(RESULTS, "probe_summary.csv"), "w", newline="") as f:
        cols = sorted({c for r in summ for c in r})
        order = ["net", "sta", "delivers", "from_provider", "delivered", "probed",
                 "years_no_data", "band", "fs", "loc", "n_comp", "epoch_start", "epoch_end",
                 "meta_prov", "lat", "lon", "cov_min", "mb_per_day", "chan_az_dip"]
        cols = [c for c in order if c in cols] + [c for c in cols if c not in order]
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(summ)
    print("\n".join(log))
    print()
    print("  {:<12}{:<9}{:<14}{:>10}{:<14}{:>9}".format(
        "station", "delivers", "from", "days ok", "  years empty", "MB/day"))
    for r in summ:
        print("  {:<12}{:<9}{:<14}{:>10}{:<14}{:>9}".format(
            f"{r['net']}.{r['sta']}", r["delivers"], r.get("from_provider", ""),
            f"{r.get('delivered','?')}/{r.get('probed','?')}",
            "  " + (r.get("years_no_data") or "-"), str(r.get("mb_per_day", ""))))
    print(f"\n-> {RESULTS}/probe_summary.csv, probe_delivery.csv")


if __name__ == "__main__":
    main()
