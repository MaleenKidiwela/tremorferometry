#!/usr/bin/env python
"""OPEN THE 2011 SEAL, once, and record it.  Prefix `m_`.

  1. verify every entry of `results/sealed_2011_manifest.sha256` (symlinks are followed -- 317 of
     the 5,062 entries are links into lfe_extension/traces, created by the downloader);
  2. write `results/SEALED_2011_OPENED_<UTC timestamp>` with the reason and the verification
     result;
  3. append the same line to `results/SEALED_2011_DO_NOT_READ` so the marker and the seal note
     can never disagree.

Nothing under traces_sealed/ is modified.  The seal is opened exactly once: if an
`SEALED_2011_OPENED_*` marker already exists this script re-verifies and exits 0 without writing a
second marker.

  python scripts/m_seal_open.py [--workers 16]
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
RES = f"{STUDY}/results"
SEALED = f"{STUDY}/traces_sealed"
MANIFEST = f"{RES}/sealed_2011_manifest.sha256"
NOTE = f"{RES}/SEALED_2011_DO_NOT_READ"
REASON = ("confirmation re-fire of the WA 2010 candidate templates on the sealed 2011 window; "
          "nothing tuned -- the templates, the 8xMAD matched filter, the stacking and the "
          "locator are the unchanged 2010 code (agent M, PLAN.md sec.15 / "
          "results/WA_2010_RESULTS.md addendum).")


def sha256_file(p, buf=1 << 20):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while True:
            b = f.read(buf)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=16)
    a = ap.parse_args()
    t0 = time.time()
    rows = []
    with open(MANIFEST) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            h, rel = line.split("  ", 1)
            rows.append((h, rel))
    print(f"[seal] {len(rows)} manifest entries", flush=True)

    def check(hr):
        h, rel = hr
        p = os.path.join(SEALED, rel)
        if not os.path.exists(p):
            return ("missing", rel)
        try:
            return ("ok", rel) if sha256_file(p) == h else ("mismatch", rel)
        except Exception as e:                                              # noqa: BLE001
            return (f"error:{type(e).__name__}", rel)

    n_ok, bad = 0, []
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        for i, (st, rel) in enumerate(ex.map(check, rows, chunksize=8), 1):
            if st == "ok":
                n_ok += 1
            else:
                bad.append((st, rel))
            if i % 1000 == 0:
                print(f"   {i}/{len(rows)}  {time.time()-t0:.0f}s  bad={len(bad)}", flush=True)
    dt = time.time() - t0
    print(f"[seal] verified {n_ok}/{len(rows)} in {dt:.0f}s; {len(bad)} not ok", flush=True)
    for st, rel in bad[:20]:
        print(f"   BAD {st} {rel}")
    if bad:
        raise SystemExit(f"seal verification FAILED on {len(bad)} files")

    prev = sorted(glob.glob(f"{RES}/SEALED_2011_OPENED_*"))
    ts = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    if prev:
        print(f"[seal] already opened: {os.path.basename(prev[0])} -- integrity re-verified, "
              "no second marker written")
        return
    line = (f"SEAL OPENED {ts}\nby: agent M ({os.path.basename(sys.argv[0])}, pid {os.getpid()})\n"
            f"reason: {REASON}\n"
            f"integrity: {n_ok}/{len(rows)} manifest entries verified sha256-clean in {dt:.0f}s "
            f"(317 of them are symlinks into lfe_extension/traces and were followed)\n"
            f"scope: traces_sealed/ is READ ONLY from here on; it is cached to cache20hz/<NET.STA>"
            "/2011/<DOY>.npz with the frozen chain and never modified.\n")
    with open(f"{RES}/SEALED_2011_OPENED_{ts.replace(':', '')}", "w") as f:
        f.write(line)
    with open(NOTE, "a") as f:
        f.write("\n" + line)
    print(line)


if __name__ == "__main__":
    main()
