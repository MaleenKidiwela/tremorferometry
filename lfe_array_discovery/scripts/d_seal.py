#!/usr/bin/env python
"""Seal the 2011-07-20..2011-09-10 re-fire set: hash it and mark it DO-NOT-READ.

Writes results/sealed_2011_manifest.sha256 (sha256 <two spaces> relative path, sorted) over every
file under traces_sealed/, and results/SEALED_2011_DO_NOT_READ.  Nothing in traces_sealed/ is ever
cached, read or scored until Step 7 of the plan opens it once.
"""
from __future__ import annotations
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d_common as D                                                   # noqa: E402


def main():
    files = []
    for dp, _, fns in os.walk(D.TRACES_SEALED):
        for fn in fns:
            p = os.path.join(dp, fn)
            files.append((os.path.relpath(p, D.TRACES_SEALED), p))
    files.sort()
    out = os.path.join(D.RES, "sealed_2011_manifest.sha256")
    nb = 0
    with open(out, "w") as f:
        for rel, p in files:
            nb += os.path.getsize(p)
            f.write(f"{D.sha256_file(p)}  {rel}\n")
    with open(os.path.join(D.RES, "SEALED_2011_DO_NOT_READ"), "w") as f:
        f.write("SEALED SET - 2011-07-20..2011-09-10 N-Olympic re-fire test.\n"
                "traces_sealed/ must NOT be read, cached, plotted or scored until Step 7 of\n"
                "PLAN.md opens it ONCE with the frozen W1 templates.  Integrity: compare against\n"
                "results/sealed_2011_manifest.sha256 before opening.\n"
                f"{len(files)} files, {nb/2**30:.2f} GiB.\n")
    D.log(f"[seal] {len(files)} files, {nb/2**30:.2f} GiB hashed -> {out}")


if __name__ == "__main__":
    main()
