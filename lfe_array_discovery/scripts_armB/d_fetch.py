#!/usr/bin/env python
"""Run one download queue with the project's FDSN day fetcher (REUSED, not rewritten).

  * `r_download_day.fetch_day` verbatim (retry-wrapped FDSN, `.nd` no-data markers, native sample
    rate, nothing filtered/merged, resumable) behind a <= 8-thread pool;
  * every queued station-day that already exists in a READ-ONLY sibling trace tree
    (../lfe_extension/traces, ../lfe_picker/traces, ../lfe_catalog_repair/traces) is SYMLINKED
    instead of being re-fetched, provided its channels match the band we want;
  * appends to results/download_manifest.csv (tranche, net, sta, day, status, bytes, chain,
    linked, + loc/band/fs/coverage_frac/n_channels/n_frag/note).

  python scripts/d_fetch.py results/queues/q_T2A_thin.csv [--workers 8] [--limit N]
"""
from __future__ import annotations
import os
import sys
import csv
import time
import argparse
import threading
from datetime import date, datetime, timezone
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "..", "lfe_picker", "scripts"))
import d_common as D                                                   # noqa: E402
from r_download_day import fetch_day                                   # noqa: E402

COLS = ["tranche", "net", "sta", "day", "status", "bytes", "chain", "linked", "loc", "band",
        "fs", "three_c", "n_channels", "coverage_frac", "n_frag", "outroot", "group",
        "updated_utc", "note"]
_lock = threading.Lock()
_n = {"done": 0, "data": 0, "nodata": 0, "error": 0, "linked": 0, "bytes": 0}
_t0 = time.time()


def day_path(outroot, net, sta, d):
    return os.path.join(outroot, f"{net}.{sta}", f"{d.year}",
                        f"{d.timetuple().tm_yday:03d}.mseed")


def link_existing(net, sta, band, d, outroot):
    """Symlink a matching day from a read-only sibling tree (never copied, never modified)."""
    ns = f"{net}.{sta}"
    doy = d.timetuple().tm_yday
    dst = day_path(outroot, net, sta, d)
    if os.path.exists(dst) or os.path.exists(dst + ".nd"):
        return "present"
    for root in D.LINK_ROOTS:
        src = os.path.join(root, ns, str(d.year), f"{doy:03d}.mseed")
        if not os.path.exists(src):
            continue
        try:
            from obspy import read
            chans = {tr.stats.channel for tr in read(src, headonly=True)}
        except Exception:
            continue
        if not any(c.startswith(band) for c in chans):
            continue
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        try:
            os.symlink(os.path.realpath(src), dst)
        except FileExistsError:
            pass
        return "linked"
    return ""


def candidates(row):
    """[(loc, band, fs), ...] in preference order: the chosen band then the declared alternates."""
    out = [((row.get("loc") or "").strip(), row["band"], row["fs"])]
    for a in (row.get("alts") or "").split(";"):
        if not a.strip():
            continue
        p = a.split("|")
        if len(p) == 3:
            out.append((p[0].strip(), p[1], p[2]))
    return out


def one(row):
    net, sta = row["net"], row["sta"]
    outroot = row["outroot"]
    d = date.fromisoformat(row["date"])
    cands = candidates(row)
    linked = ""
    try:
        linked = link_existing(net, sta, cands[0][1], d, outroot)
    except Exception as e:
        linked = f"link_err:{type(e).__name__}"
    tried, rec, loc, band, fs = [], None, cands[0][0], cands[0][1], cands[0][2]
    for k, (loc, band, fs) in enumerate(cands):
        rec = fetch_day(net, sta, band, d, loc, outroot, prov="IRIS", tries=3, timeout=180)
        tried.append(band)
        if rec["status"] == "data" or k == len(cands) - 1:
            break
        # nodata on a DECLARED band -> drop the marker and try the next declared band
        nd = day_path(outroot, net, sta, d) + ".nd"
        if os.path.exists(nd):
            try:
                os.remove(nd)
            except OSError:
                pass
    note = rec["note"]
    if len(tried) > 1:
        note = (note + ";" if note else "") + "bands_tried:" + ",".join(tried)
    try:
        chain = D.chain_for(float(fs))
    except (TypeError, ValueError):
        chain = row["chain"]
    out = dict(tranche=row["tranche"], net=net, sta=sta, day=row["date"], status=rec["status"],
               bytes=rec["bytes"], chain=chain, linked=linked, loc=loc, band=band,
               fs=rec["fs"] or fs, three_c=row["three_c"], n_channels=rec["n_channels"],
               coverage_frac=rec["coverage_frac"], n_frag=rec["n_frag"], outroot=outroot,
               group=row["group"],
               updated_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
               note=note)
    with _lock:
        _n["done"] += 1
        _n[out["status"]] = _n.get(out["status"], 0) + 1
        if linked == "linked":
            _n["linked"] += 1
        try:
            _n["bytes"] += int(out["bytes"] or 0)
        except (TypeError, ValueError):
            pass
        with open(D.MANIFEST, "a", newline="") as f:
            csv.DictWriter(f, COLS, extrasaction="ignore").writerow(out)
        if _n["done"] % 200 == 0:
            el = time.time() - _t0
            D.log(f"  {_n['done']} done  data={_n.get('data',0)} nodata={_n.get('nodata',0)} "
                  f"err={_n.get('error',0)} linked={_n['linked']}  {_n['bytes']/2**30:.1f} GiB  "
                  f"{el/60:.1f} min  ({_n['done']/max(el,1):.2f} sta-day/s)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("queue")
    ap.add_argument("--workers", type=int, default=D.MAX_CONCURRENT_FDSN)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    with open(a.queue) as f:
        rows = list(csv.DictReader(f))
    done = set()
    if os.path.exists(D.MANIFEST):
        with open(D.MANIFEST) as f:
            for r in csv.DictReader(f):
                if r.get("status") in ("data", "nodata"):
                    done.add((r["net"], r["sta"], r["day"], r["outroot"]))
    else:
        with open(D.MANIFEST, "w", newline="") as f:
            csv.DictWriter(f, COLS).writeheader()
    todo = [r for r in rows
            if (r["net"], r["sta"], r["date"], r["outroot"]) not in done]
    if a.limit:
        todo = todo[:a.limit]
    D.log(f"[fetch] {os.path.basename(a.queue)}: {len(todo)} station-days to fetch "
          f"({len(rows)-len(todo)} already terminal), {min(a.workers, D.MAX_CONCURRENT_FDSN)} "
          f"concurrent FDSN requests")
    with ThreadPoolExecutor(max_workers=min(a.workers, D.MAX_CONCURRENT_FDSN)) as ex:
        list(ex.map(one, todo))
    el = (time.time() - _t0) / 60
    D.log(f"[fetch] {os.path.basename(a.queue)} finished: {_n} in {el:.1f} min")
    return 0 if _n.get("error", 0) <= 0.25 * max(len(todo), 1) else 1


if __name__ == "__main__":
    sys.exit(main())
