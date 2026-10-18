#!/usr/bin/env python
"""SVI_REF R9 -- the delete-CANDIDATE list (NOTHING is deleted; per the user, nothing goes before the Vancouver work is finished
and the user approves an exact list; this is input to a later JOINT SVI + Washington decision).

Rules (Merlin):
  * the 18 stations shared with the Washington extension (lfe_extension_wa/results/step4b/setup/stations.csv; WA reads
    lfe_extension/traces through wa_common.TRACE_ROOTS and the SVI 20 Hz cache) are HELD -- raw AND cache -- until Washington
    has scanned and cut;
  * only CN.NLLB (not in the WA roster) can be a candidate, and only for station-days that (i) the R8 cut read OK and whose
    cut verified 20/20, (ii) R2 measured (raw health row), and (iii) no symlink of lfe_array_discovery (T7 tranche) or
    lfe_picker points at; its 20 Hz cache is listed separately (regenerable only while the raw exists);
  * truth-year raw (2010-2013) is kept; results/, the cut store, the logs are kept.
Outputs OUT/r9/: DELETE_CANDIDATES_svi.csv (one row per file), KEEP_svi_summary.csv (per station-year, with the reason),
DELETE_CANDIDATES_svi.txt (summary), symlink_readers.csv.

  python scripts/svi_ref_delete.py
"""
from __future__ import annotations
import os, sys, json, glob, subprocess, warnings
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svi_ref_common as R
import c_common as C

O9 = f"{R.OUT}/r9"
CAND = "CN.NLLB"


def symlink_readers():
    rows = []
    for proj in ("lfe_array_discovery", "lfe_picker", "lfe_extension_wa"):
        root = f"{R.ROOT}/{proj}"
        if not os.path.isdir(root): continue
        try:
            out = subprocess.run(["find", root, "-type", "l", "-lname", f"*lfe_extension/traces/*"], capture_output=True, text=True,
                                 timeout=1800).stdout.split()
        except Exception as e:
            rows.append(dict(project=proj, netsta="?", n_links=-1, note=f"scan failed {e}")); continue
        s = pd.Series([os.readlink(p).split("/lfe_extension/traces/")[1].split("/")[0] for p in out]) if out else pd.Series(dtype=str)
        for ns, n in s.value_counts().items():
            rows.append(dict(project=proj, netsta=ns, n_links=int(n), note=""))
    return pd.DataFrame(rows)


def main():
    warnings.filterwarnings("ignore"); R.ensure(O9)
    wa = pd.read_csv(f"{R.ROOT}/lfe_extension_wa/results/step4b/setup/stations.csv")
    wa_set = set(wa.netsta)
    L = symlink_readers(); L.to_csv(f"{O9}/symlink_readers.csv", index=False)
    linked = {(r.netsta) for r in L.itertuples() if r.n_links > 0}
    cut_ok = os.path.exists(f"{R.CUT}/VERIFY.txt") and open(f"{R.CUT}/VERIFY.txt").read().startswith("PASS")
    SD = pd.concat([pd.read_parquet(p) for p in glob.glob(f"{R.CUT}/stationdays_*.parquet")], ignore_index=True) if cut_ok else pd.DataFrame(columns=["path", "status"])
    okpaths = set(SD[SD.status == "ok"].path) if len(SD) else set()
    H = pd.read_parquet(f"{R.OUT}/r2/r2_health_stationday.parquet") if os.path.exists(f"{R.OUT}/r2/r2_health_stationday.parquet") else pd.DataFrame(columns=["sta", "date"])
    r2days = {(r.sta, r.date) for r in H.itertuples()}
    cand, keep = [], []
    for ns in sorted(os.listdir(C.TRACES)):
        for y in sorted(os.listdir(f"{C.TRACES}/{ns}")):
            d = f"{C.TRACES}/{ns}/{y}"
            fs = [f for f in os.listdir(d) if f.endswith(".mseed")]
            nb = sum(os.path.getsize(f"{d}/{f}") for f in fs)
            reasons = []
            if int(y) <= 2013: reasons.append("truth-year raw (2010-2013): kept")
            if ns in wa_set: reasons.append("shared with the WA extension (step4b roster): hold raw and cache until WA has scanned and cut")
            if ns in linked: reasons.append(f"symlinked by {','.join(L[(L.netsta == ns) & (L.n_links > 0)].project)}")
            if ns != CAND and not reasons: reasons.append("not the NLLB candidate")
            if reasons:
                keep.append(dict(kind="raw", netsta=ns, year=y, n_files=len(fs), GB=round(nb / 1e9, 3), reason="; ".join(reasons))); continue
            for f in fs:
                p = f"{d}/{f}"; day = str((pd.Timestamp(f"{y}-01-01") + pd.Timedelta(days=int(f[:3]) - 1)).date())
                why = []
                if p not in okpaths: why.append("not fully cut + verified")
                if ("NLLB", day) not in r2days: why.append("no R2 health row")
                if why: keep.append(dict(kind="raw", netsta=ns, year=y, n_files=1, GB=round(os.path.getsize(p) / 1e9, 4), reason="; ".join(why)))
                else: cand.append(dict(kind="raw", path=p, netsta=ns, day=day, bytes=os.path.getsize(p),
                                       reason="NLLB (not in the WA roster, no symlink reader), cut read OK + verified 20/20, R2 measured"))
    for root, kind in ((C.CACHE, "cache20hz"), (R.CACHE_R1, "cache20hz_r1")):
        if not os.path.isdir(root): continue
        for ns in sorted(os.listdir(root)):
            for y in sorted(os.listdir(f"{root}/{ns}")):
                d = f"{root}/{ns}/{y}"; fs = [f for f in os.listdir(d) if f.endswith(".npz")]
                nb = sum(os.path.getsize(f"{d}/{f}") for f in fs)
                if ns == CAND and int(y) >= 2014:
                    for f in fs:
                        cand.append(dict(kind=kind, path=f"{d}/{f}", netsta=ns, day=f"{y}:{f[:3]}", bytes=os.path.getsize(f"{d}/{f}"),
                                         reason="NLLB 20 Hz cache (regenerable from raw only while the raw exists; SVI_REF stages done)"))
                else:
                    keep.append(dict(kind=kind, netsta=ns, year=y, n_files=len(fs), GB=round(nb / 1e9, 3),
                                     reason="shared with WA (cache) or truth-year cache: held" if ns in wa_set or int(y) <= 2013 else "held"))
    Cd = pd.DataFrame(cand); K = pd.DataFrame(keep)
    Cd.to_csv(f"{O9}/DELETE_CANDIDATES_svi.csv", index=False)
    Ks = K.groupby(["kind", "netsta", "year", "reason"]).agg(n_files=("n_files", "sum"), GB=("GB", "sum")).reset_index() if len(K) else K
    Ks.to_csv(f"{O9}/KEEP_svi_summary.csv", index=False)
    tb = Cd.bytes.sum() / 1e12 if len(Cd) else 0.0
    msg = (f"SVI_REF R9 delete-CANDIDATE list (NOTHING DELETED; for a later joint SVI + Washington decision by the user):\n"
           f"  candidates {len(Cd)} files, {tb:.3f} TB ({Cd.groupby('kind').bytes.sum().div(1e9).round(1).to_dict() if len(Cd) else {}} GB)\n"
           f"  held {Ks.n_files.sum() if len(Ks) else 0} files, {Ks.GB.sum() / 1e3 if len(Ks) else 0:.3f} TB; reasons "
           f"{Ks.groupby('reason').GB.sum().div(1e3).round(3).to_dict() if len(Ks) else {}}\n"
           f"  cut verified: {cut_ok}; symlink readers: {L[L.n_links > 0].groupby('project').n_links.sum().to_dict() if len(L) else {}}\n")
    open(f"{O9}/DELETE_CANDIDATES_svi.txt", "w").write(msg)
    R.log("[r9] " + msg)


if __name__ == "__main__":
    main()
