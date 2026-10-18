#!/usr/bin/env python
"""Track B step 2 - fetch ONE station-day of raw 3-C miniSEED and verify it.

Adapted from ../lfe_catalog_repair/scripts/d3_download.py (retry-wrapped FDSN, `.nd`
no-data markers, resumable). NATIVE sample rate: no decimation, no filtering, no response
removal, no merging of gaps - the bytes that come back are what is written.

Layout: <outroot>/<NET>.<STA>/<YYYY>/<DOY>.mseed   (DOY zero-padded to 3)
        <...>.mseed.nd  = verified empty day (no data at the provider)

fetch_day() returns a dict with status/bytes/n_channels/fs/coverage_frac/n_frag/note that
the runner appends to results/download_manifest.csv.

Standalone:  python scripts/download_day.py --net CN --sta PGC --band BH --date 2014-03-01
"""
import os
import time
import argparse
from datetime import date

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_OUTROOT = os.path.join(HERE, "traces")
NRCAN = "https://earthquakescanada.nrcan.gc.ca/fdsnws"
_CLIENTS = {}


def get_client(prov="IRIS", timeout=120):
    """One obspy Client per (provider, thread)."""
    import threading
    from obspy.clients.fdsn import Client
    key = (prov, threading.get_ident())
    c = _CLIENTS.get(key)
    if c is None:
        c = Client(prov, timeout=timeout)
        _CLIENTS[key] = c
    return c


def day_path(outroot, net, sta, d):
    return os.path.join(outroot, f"{net}.{sta}", f"{d.year}",
                        f"{d.timetuple().tm_yday:03d}.mseed")


def verify(path, band):
    """Read headers only -> (n_channels, fs, coverage_frac, n_frag, note)."""
    from obspy import read
    st = read(path, headonly=True)
    chans = {}
    for tr in st:
        chans.setdefault(tr.stats.channel, []).append(tr)
    cov, frags, fss = [], [], set()
    for ch, trs in chans.items():
        fs = float(trs[0].stats.sampling_rate)
        fss.add(fs)
        secs = sum(tr.stats.npts / fs for tr in trs)
        cov.append(min(secs / 86400.0, 1.0))
        frags.append(len(trs))
    note = ""
    if len(fss) > 1:
        note = "mixed_fs:" + ",".join(f"{f:g}" for f in sorted(fss))
    if any(not c.startswith(band) for c in chans):
        note = (note + ";" if note else "") + "chan:" + ",".join(sorted(chans))
    return (len(chans), max(fss) if fss else 0.0,
            min(cov) if cov else 0.0, max(frags) if frags else 0, note)


def fetch_day(net, sta, band, d, loc="", outroot=DEFAULT_OUTROOT, prov="IRIS",
              tries=3, timeout=120):
    """Download one station-day. Returns a manifest dict. Never raises."""
    from obspy import UTCDateTime
    p = day_path(outroot, net, sta, d)
    rec = dict(net=net, sta=sta, band=band, date=d.isoformat(),
               doy=f"{d.timetuple().tm_yday:03d}", status="", bytes="", n_channels="",
               fs="", coverage_frac="", n_frag="", note="")
    if os.path.exists(p):                                    # already on disk -> verify
        try:
            n, fs, cov, nf, note = verify(p, band)
            rec.update(status="data", bytes=os.path.getsize(p), n_channels=n, fs=f"{fs:g}",
                       coverage_frac=f"{cov:.4f}", n_frag=nf, note=(note + ";resumed").strip(";"))
        except Exception as e:
            rec.update(status="error", note=f"unreadable_existing:{type(e).__name__}")
        return rec
    if os.path.exists(p + ".nd"):
        rec.update(status="nodata", bytes=0, note="resumed")
        return rec
    t0 = UTCDateTime(d.isoformat())
    last = ""
    for k in range(tries):
        try:
            st = get_client(prov, timeout).get_waveforms(
                net, sta, loc if loc else "*", band + "?", t0, t0 + 86400)
            if loc:
                st = st.select(location=loc)
            st = st.select(channel=band + "?")
            if len(st) == 0:
                raise_nodata = True
            else:
                raise_nodata = False
            if raise_nodata:
                os.makedirs(os.path.dirname(p), exist_ok=True)
                open(p + ".nd", "w").close()
                rec.update(status="nodata", bytes=0)
                return rec
            os.makedirs(os.path.dirname(p), exist_ok=True)
            st.write(p + ".tmp", format="MSEED")
            n, fs, cov, nf, note = verify(p + ".tmp", band)
            os.replace(p + ".tmp", p)
            rec.update(status="data", bytes=os.path.getsize(p), n_channels=n, fs=f"{fs:g}",
                       coverage_frac=f"{cov:.4f}", n_frag=nf, note=note)
            return rec
        except Exception as e:
            msg = f"{type(e).__name__}:{str(e)[:120]}"
            last = msg
            if "No data" in str(e) or "204" in str(e) or "FDSNNoDataException" in msg:
                os.makedirs(os.path.dirname(p), exist_ok=True)
                open(p + ".nd", "w").close()
                rec.update(status="nodata", bytes=0)
                return rec
            for junk in (p + ".tmp",):
                if os.path.exists(junk):
                    try:
                        os.remove(junk)
                    except OSError:
                        pass
            time.sleep(3 * (k + 1) ** 2)                      # 3, 12, 27 s backoff
    rec.update(status="error", note=last)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True)
    ap.add_argument("--sta", required=True)
    ap.add_argument("--band", required=True, help="HH / BH / EH (band code, no component)")
    ap.add_argument("--date", required=True, help="YYYY-MM-DD")
    ap.add_argument("--loc", default="")
    ap.add_argument("--outroot", default=DEFAULT_OUTROOT)
    ap.add_argument("--prov", default="IRIS")
    a = ap.parse_args()
    print(fetch_day(a.net, a.sta, a.band, date.fromisoformat(a.date), a.loc,
                    os.path.abspath(a.outroot), a.prov))


if __name__ == "__main__":
    main()
