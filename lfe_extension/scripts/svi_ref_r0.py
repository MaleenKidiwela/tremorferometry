#!/usr/bin/env python
"""SVI_REF R0 -- stop the idle Track C watcher and freeze catalogue v1.1 as the baseline.

  check     the watcher (scripts/c_watch.py, PID 629057, idle 14 days) must be gone: results/trackC/STOP was touched
            (2026-10-01T22:28:14Z, it exits at its next poll); waits up to --wait-s, never kills anything; writes OUT/r0/watcher.json
  snapshot  sha256 of every v1.1 output (results/trackC, results/trackD, results/trackE, the frozen operating point,
            data/catalogue/families_v1_1.csv) -> OUT/r0/v11_hashes.csv
  verify    re-hash and compare (run at the end of the queue): OUT/r0/v11_verify.json; exit 7 if anything changed
"""
from __future__ import annotations
import os, sys, json, glob, time, argparse, subprocess
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svi_ref_common as R

O0 = f"{R.OUT}/r0"
PID = 629057


def alive(pid):
    try:
        out = subprocess.run(["ps", "-o", "pid=,cmd=", "-p", str(pid)], capture_output=True, text=True).stdout.strip()
        return out
    except Exception:
        return ""


def v11_files():
    fs = []
    for d in (f"{R.RES}/trackC", f"{R.RES}/trackD", f"{R.RES}/trackE"):
        fs += [p for p in glob.glob(f"{d}/**/*", recursive=True) if os.path.isfile(p)]
    fs += [f"{R.RES}/frozen_operating_point.json", f"{R.ROOT}/data/catalogue/families_v1_1.csv"]
    return sorted(p for p in fs if os.path.exists(p) and not p.endswith("/STOP"))


def hashes(fs):
    def one(p):
        st = os.stat(p)
        return dict(path=p, bytes=st.st_size, mtime=st.st_mtime, sha256=R.sha256(p))
    with ThreadPoolExecutor(8) as ex:
        return pd.DataFrame(list(ex.map(one, fs)))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("what", choices=["check", "snapshot", "verify"])
    ap.add_argument("--wait-s", type=float, default=1200.0)
    a = ap.parse_args(); R.ensure(O0)
    if a.what == "check":
        t0 = time.time(); s = alive(PID)
        while s and time.time() - t0 < a.wait_s:
            time.sleep(30); s = alive(PID)
        log_tail = ""
        p = f"{R.EXT}/logs/trackC_watch.log"
        if os.path.exists(p):
            log_tail = open(p).read().splitlines()[-1:]
        r = dict(pid=PID, still_running=bool(s), ps=s, stop_file_present=os.path.exists(f"{R.RES}/trackC/STOP"),
                 watcher_log_last=log_tail, checked_utc=pd.Timestamp.now("UTC").isoformat(),
                 note="STOP touched 2026-10-01T22:28:14Z; nothing was killed")
        json.dump(r, open(f"{O0}/watcher.json", "w"), indent=1)
        R.log(f"[r0] watcher {r}")
    elif a.what == "snapshot":
        p = f"{O0}/v11_hashes.csv"
        if os.path.exists(p):
            R.log(f"[r0] {p} exists (frozen); not rewritten"); return
        H = hashes(v11_files()); H.to_csv(p, index=False)
        R.log(f"[r0] v1.1 snapshot: {len(H)} files, {H.bytes.sum() / 1e9:.2f} GB -> {p}")
    else:
        H0 = pd.read_csv(f"{O0}/v11_hashes.csv")
        H1 = hashes([p for p in H0.path if os.path.exists(p)])
        J = H0.merge(H1, on="path", how="left", suffixes=("", "_now"))
        changed = J[(J.sha256 != J.sha256_now)]
        new = sorted(set(v11_files()) - set(H0.path))
        r = dict(n_files=len(H0), n_changed_or_missing=len(changed), changed=changed.path.tolist()[:50], new_files=new[:50],
                 verified_utc=pd.Timestamp.now("UTC").isoformat())
        json.dump(r, open(f"{O0}/v11_verify.json", "w"), indent=1)
        R.log(f"[r0] verify: {len(changed)} of {len(H0)} v1.1 files changed or missing; {len(new)} new files")
        if len(changed):
            sys.exit(7)


if __name__ == "__main__":
    main()
