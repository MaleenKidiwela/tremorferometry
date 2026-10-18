#!/usr/bin/env python
"""DUC_FIN step 5 -- DELETE-CANDIDATE list (nothing is deleted) for the Ducellier stations, 2008-01-01 .. 2012-07-31.

Kinds: raw miniSEED (traces/, traces_duc/), cache_duc (his chain, 20 sps), cache20hz (t3 chain) of his 47 stations, and
results/duc_fix/prep/<date>.npz (one file per day, all stations).  A file is a CANDIDATE only when
  * every DUC_FIN stage that reads it has finished (feature MF, HF, PNET, cut, fit: marker files in OUT), and
  * (raw) the cut read it without error (results/cut/duc/stationdays status ok), and
  * no other job still needs that station-day:
      - sh_back (Shelly backward, HALTED before BK2): S6 stations NC.KMPB NC.KMR PB.B049 PB.B932 PB.B933 PB.B934,
        2008-01-01 .. 2017-12-31 (raw + cache20hz)   [of these only PB.B933 is a Ducellier station]
      - DUC_FIX full-record extension (scripts/duc_ext_fixed_queue.sh, months 2012-08 .. 2026-09 not scanned): its final
        re-templating of the 2011-10..12 new epochs uses 24-month windows (to 2013-10..12) and then RE-SCANS every
        month from the epoch start (pass 2) -> raw + cache_duc + prep of 2011-10-01 .. 2012-07-31 are kept (new pass-2
        candidates would need new windows cut from raw)
      - the known-answer window 2008-04-01 .. 2008-05-12 (duc_fix stage 1, D2011): cache_duc + prep kept (raw listed)
      - Mendocino discovery (lfe_array_discovery, 572 Shelly days 2018-2024) and the ~600 sampled Shelly days 2018-2023:
        checked programmatically for overlap with 2008-2012 (none expected)
Regenerate column: how the file could be rebuilt.  Writes CUT_OUT/../DELETE_CANDIDATES_ducellier.csv,
KEEP_ducellier.csv, DELETE_CANDIDATES_ducellier.txt (summary) and the consumer scan (other scripts that read these
folders) into the summary.
"""
import os, sys, json, glob, subprocess
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import duc_fin_common as F
import duc_fix_common as DC
from t_common import ROOT, WS, TR, log

S6 = ["NC.KMPB", "NC.KMR", "PB.B049", "PB.B932", "PB.B933", "PB.B934"]
KEEP_RESCAN = ("2011-10-01", "2012-07-31")
KEEP_KA = ("2008-04-01", "2008-05-12")
STAGES = ["MF_DONE", "HF_DONE", "PNET_DONE", "CUT_DONE", "FIT_DONE"]


def day_of(y, fn):
    return str((pd.Timestamp(f"{y}-01-01") + pd.Timedelta(days=int(fn[:3]) - 1)).date())


def files_of(root, ns, kind, ext):
    out = []
    for y in range(2008, 2013):
        d = f"{root}/{ns}/{y}"
        if not os.path.isdir(d): continue
        for fn in os.listdir(d):
            if not fn.endswith(ext) or ".tmp" in fn: continue
            day = day_of(y, fn)
            if F.T0 <= day < F.T1: out.append((f"{d}/{fn}", ns, day, kind))
    return out


def other_days():
    """days of 2008-2012 that the Mendocino discovery / Shelly sampled days use (expected: none)."""
    s = set()
    p = f"{ROOT}/lfe_array_discovery/results/ab_royer/tremor/scan_days_MENDO.csv"
    if os.path.exists(p):
        s |= {str(x)[:10] for x in pd.read_csv(p).iloc[:, 0]}
    for h in ("build", "test"):
        q = f"{WS}/results/aa_shelly/stage2_days_{h}.json"
        if os.path.exists(q):
            for k, v in json.load(open(q)).items():
                if k.split("|")[0] == "Shelly2025": s |= set(v)      # the Ducellier2022 keys belong to the finished AA study
    return {d for d in s if F.T0 <= d < F.T1}


def consumer_scan():
    cmd = (f"grep -rlE 'lfe_south/(traces|cache_duc|cache20hz)|traces_duc|cache_duc' {ROOT} --include=*.py --include=*.sh "
           f"--exclude-dir=.git 2>/dev/null")
    try: out = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=600).stdout.split()
    except Exception as e: out = [f"scan failed: {e}"]
    return sorted(out)


def main():
    out = F.CUT_OUT; parent = os.path.dirname(os.path.abspath(out))
    done = {s: os.path.exists(f"{F.OUT}/{s}") for s in STAGES}
    stas = F.stations()
    rows = []
    for ns in stas:
        rows += files_of(TR, ns, "raw_traces", ".mseed") + files_of(DC.TR2, ns, "raw_traces_duc", ".mseed")
        rows += files_of(f"{WS}/cache_duc", ns, "cache_duc", ".npz") + files_of(f"{WS}/cache20hz", ns, "cache20hz", ".npz")
    for p in sorted(glob.glob(f"{F.FIX}/prep/*.npz")):
        d = os.path.basename(p)[:10]
        if F.T0 <= d < F.T1: rows.append((p, "ALL", d, "prep_duc_fix"))
    A = pd.DataFrame(rows, columns=["path", "netsta", "day", "kind"])
    A["size_bytes"] = [os.path.getsize(p) for p in A.path]
    SD = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{out}/stationdays_*.parquet"))], ignore_index=True) \
        if glob.glob(f"{out}/stationdays_*.parquet") else pd.DataFrame(columns=["path", "status", "n_win"])
    cutok = dict(zip(SD.path, SD.status)); nwin = dict(zip(SD.path, SD.n_win))
    od = other_days()
    why = []
    for r in A.itertuples():
        w = []
        if not all(done.values()): w.append("DUC_FIN stages not all done: " + ",".join(k for k, v in done.items() if not v))
        if r.kind.startswith("raw"):
            st = cutok.get(r.path)
            if st is None: w.append("raw file not read by the DUC_FIN cut")
            elif st != "ok": w.append(f"cut status {st}")
        if r.netsta in S6 and r.day < "2018-01-01" and r.kind in ("raw_traces", "raw_traces_duc", "cache20hz"):
            w.append("sh_back BK2 (S6 2008-2017; queue halted, not finished)")
        if KEEP_RESCAN[0] <= r.day <= KEEP_RESCAN[1] and r.kind != "cache20hz":
            w.append("DUC_FIX full-record extension: final 24-month re-templating of the 2011-10..12 epochs + pass-2 re-scan of 2011-10..2012-07")
        if KEEP_KA[0] <= r.day <= KEEP_KA[1] and r.kind in ("cache_duc", "prep_duc_fix"):
            w.append("known-answer window (duc_fix stage 1)")
        if r.day in od: w.append("Mendocino discovery / Shelly sampled day")
        why.append("; ".join(w))
    A["keep_reason"] = why
    regen = {"raw_traces": "FDSN re-download (NCEDC / IRIS): python scripts/duc_fix_data.py --month <M> --stage fetch",
             "raw_traces_duc": "FDSN re-download (NCEDC / IRIS): python scripts/duc_fix_data.py --month <M> --stage fetch",
             "cache_duc": "python scripts/duc_fix_data.py --month <M> --stage cache (needs the raw day file)",
             "cache20hz": "python scripts/t3_cache.py <plan> (needs the raw day file)",
             "prep_duc_fix": "python scripts/duc_fix_prep.py --t0 <D> --t1 <D+1> (needs cache_duc of that day)"}
    A["regenerate"] = A.kind.map(regen)
    A["windows_cut"] = [int(nwin.get(p, 0)) if k.startswith("raw") else np.nan for p, k in zip(A.path, A.kind)]
    A["size_GB"] = (A.size_bytes / 1e9).round(4)
    C = A[A.keep_reason == ""].drop(columns=["keep_reason"]); K = A[A.keep_reason != ""]
    C.to_csv(f"{parent}/DELETE_CANDIDATES_ducellier.csv", index=False); K.to_csv(f"{parent}/KEEP_ducellier.csv", index=False)
    by = A.assign(cand=A.keep_reason == "").groupby(["kind", "cand"]).size_bytes.agg(["size", "sum"]).reset_index()
    by["TB"] = (by["sum"] / 1e12).round(3)
    kr = K.keep_reason.str.split("; ").explode().value_counts()
    kst = K.groupby("keep_reason").agg(files=("path", "size"), TB=("size_bytes", lambda x: round(x.sum() / 1e12, 3)),
                                       stations=("netsta", lambda x: ",".join(sorted(set(x))[:60])),
                                       first_day=("day", "min"), last_day=("day", "max")).reset_index()
    kst.to_csv(f"{parent}/KEEP_ducellier_summary.csv", index=False)
    cons = consumer_scan()
    txt = [f"DELETE-CANDIDATE list for the Ducellier stations 2008-01-01..2012-07-31 (NOTHING DELETED) -- {pd.Timestamp.now('UTC'):%Y-%m-%d %H:%M} UTC",
           f"candidates: {len(C)} files, {C.size_bytes.sum() / 1e12:.3f} TB -> {parent}/DELETE_CANDIDATES_ducellier.csv",
           f"kept: {len(K)} files, {K.size_bytes.sum() / 1e12:.3f} TB -> KEEP_ducellier.csv (summary KEEP_ducellier_summary.csv)",
           f"DUC_FIN stage markers: {done}", f"days of 2008-2012 used by Mendocino discovery / Shelly sampled days: {len(od)}",
           "by kind (cand = candidate):", by.to_string(index=False), "keep reasons (files):", kr.to_string(),
           "scripts that reference lfe_south raw / cache folders (for the user to check):"] + [f"  {c}" for c in cons]
    open(f"{parent}/DELETE_CANDIDATES_ducellier.txt", "w").write("\n".join(txt) + "\n")
    log("\n".join(txt[:6]) + "\n" + by.to_string(index=False))


if __name__ == "__main__":
    main()
