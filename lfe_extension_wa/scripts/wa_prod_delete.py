#!/usr/bin/env python
"""WA PRODUCTION W13: the delete-CANDIDATE list, joint with SVI (lfe_extension/results/refine/r9/).  NOTHING is deleted.
The list goes to the user only after the cut-off table and the rates exist (the queue runs this last).

Candidates (one row per file):
  raw_svi_shared  lfe_extension/traces 2014-2026 of the stations SVI held for Washington (SVI R9 'shared with the WA
                  extension'): the WA cut read the file OK and the WA cut verified 20/20
  raw_wa_dl       lfe_extension_wa/traces_dl, traces_dl2 (WA's own downloads; real files only, symlinks listed apart):
                  the WA cut read the file OK and verified 20/20
  cache_*         cache20hz_prod (this run), cache20hz (WA trial), lfe_extension/cache20hz(_r1) (SVI): regenerable only
                  while the raw exists -- listed so the user can decide caches and raw together
Held (summary rows, with the reason): truth raw 2006-2011 (incl. the sealed months 2007-01 and 2009-05), 2012-2013 raw kept
unless WA downloaded it, files any other project's symlink reads, every file the cut did not read OK, and everything of
lfe_array_discovery / lfe_picker / lfe_catalog_repair (LISTED ONLY: those projects' decision).  results/, the cut store
and logs are never candidates.  Storage policy: finish south before north (memory), so this list is input to a later joint
decision.

  python scripts/wa_prod_delete.py
"""
from __future__ import annotations
import os, sys, json, glob, subprocess, warnings
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_prod_common as P
W = P.W

OD = f"{P.OUT}/delete"
SVI_TR = f"{P.ROOT}/lfe_extension/traces"


def readers(targets):
    """every symlink in the projects that points into one of the target roots (WA's own cut_roots views excluded)."""
    rows = []
    for proj in ("lfe_array_discovery", "lfe_picker", "lfe_extension", "lfe_extension_wa", "lfe_south", "lfe_catalog_repair"):
        root = f"{P.ROOT}/{proj}"
        if not os.path.isdir(root):
            continue
        expr = []
        for tg in targets:
            expr += ["-lname", f"*{tg}/*", "-o"]
        cmd = ["find", root, "-path", f"{P.OUT}/cut_roots", "-prune", "-o", "-type", "l", "("] + expr[:-1] + [")", "-print"]
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, timeout=7200).stdout.split()
        except Exception as e:
            rows.append(dict(project=proj, path="?", link=f"scan failed {e}")); continue
        for p in out:
            try:
                rows.append(dict(project=proj, path=p, link=os.readlink(p)))
            except OSError:
                pass
    return pd.DataFrame(rows, columns=["project", "path", "link"])


def main():
    warnings.filterwarnings("ignore"); P.ensure(OD)
    cut_ok = os.path.exists(f"{P.CUT}/VERIFY.txt") and open(f"{P.CUT}/VERIFY.txt").read().startswith("PASS")
    SD = pd.concat([pd.read_parquet(p) for p in glob.glob(f"{P.CUT}/stationdays_*.parquet")], ignore_index=True) if cut_ok else pd.DataFrame(columns=["path", "status"])
    okp = set(os.path.realpath(p) for p in SD[SD.status == "ok"].path) if len(SD) else set()
    R = readers([SVI_TR, f"{P.STUDY}/traces_dl", f"{P.STUDY}/traces_dl2"]) if not P.SMOKE else pd.DataFrame(columns=["project", "path", "link"])
    R.to_csv(f"{OD}/symlink_readers.csv", index=False)
    linked = set(os.path.realpath(x) if os.path.exists(x) else x for x in R.link if isinstance(x, str))
    svi_keep = pd.read_csv(f"{P.ROOT}/lfe_extension/results/refine/r9/KEEP_svi_summary.csv") if os.path.exists(
        f"{P.ROOT}/lfe_extension/results/refine/r9/KEEP_svi_summary.csv") else pd.DataFrame(columns=["netsta", "reason"])
    shared = set(svi_keep[svi_keep.reason.astype(str).str.contains("shared with the WA")].netsta)
    cand, keep = [], []

    def scan(root, kind, rule):
        if not os.path.isdir(root):
            return
        for ns in sorted(os.listdir(root))[: (2 if P.SMOKE else None)]:
            d0 = f"{root}/{ns}"
            if not os.path.isdir(d0):
                continue
            for y in sorted(os.listdir(d0)):
                if not y.isdigit():
                    continue
                for e in os.scandir(f"{d0}/{y}"):
                    if not e.name.endswith(".mseed"):
                        continue
                    p = e.path; rp = os.path.realpath(p); is_link = os.path.islink(p)
                    try:
                        b = 0 if is_link else e.stat().st_size
                    except OSError:
                        b = 0
                    why = rule(ns, int(y), p, rp, is_link)
                    if why is None:
                        cand.append(dict(kind=kind, path=p, netsta=ns, year=int(y), bytes=b,
                                         reason="cut read OK + verified 20/20; no other reader"))
                    else:
                        keep.append(dict(kind=kind, netsta=ns, year=int(y), bytes=b, reason=why))

    def rule_raw(ns, y, p, rp, is_link):
        if y <= 2011:
            return "truth-era raw 2006-2011 (incl. sealed months): kept"
        if is_link:
            return "symlink (frees nothing; the target belongs to another project)"
        if not cut_ok:
            return "cut not verified"
        if rp not in okp:
            return "not read OK by the WA cut"
        if p in linked or rp in linked:
            return "read by another project's symlink"
        return None

    def rule_svi(ns, y, p, rp, is_link):
        if ns not in shared:
            return "not an SVI station held for WA (SVI's own decision, R9)"
        if y < 2014:
            return "pre-2014 SVI raw: SVI's truth years / not part of the joint list"
        return rule_raw(ns, y, p, rp, is_link)

    scan(SVI_TR, "raw_svi_shared", rule_svi)
    scan(f"{P.STUDY}/traces_dl", "raw_wa_dl", rule_raw)
    scan(f"{P.STUDY}/traces_dl2", "raw_wa_dl2", rule_raw)
    for root, kind in ((P.CACHE_NEW, "cache_wa_prod"), (f"{P.STUDY}/cache20hz", "cache_wa_trial"),
                       (f"{P.ROOT}/lfe_extension/cache20hz", "cache_svi"), (f"{P.ROOT}/lfe_extension/cache20hz_r1", "cache_svi_r1")):
        if not os.path.isdir(root):
            continue
        for ns in sorted(os.listdir(root)):
            for y in sorted(os.listdir(f"{root}/{ns}")) if os.path.isdir(f"{root}/{ns}") else []:
                fs = [e for e in os.scandir(f"{root}/{ns}/{y}") if e.name.endswith(".npz")]
                nb = sum(e.stat().st_size for e in fs)
                if y.isdigit() and int(y) <= 2011:
                    keep.append(dict(kind=kind, netsta=ns, year=int(y), bytes=nb, reason="truth-era cache: kept"))
                else:
                    cand.append(dict(kind=kind, path=f"{root}/{ns}/{y}/", netsta=ns, year=int(y) if y.isdigit() else -1, bytes=nb,
                                     reason=f"20 Hz cache, regenerable from raw only while the raw exists ({len(fs)} files)"))
    for proj in ("lfe_array_discovery/traces", "lfe_picker/traces", "lfe_catalog_repair/traces"):
        root = f"{P.ROOT}/{proj}"
        if os.path.isdir(root):
            for ns in sorted(os.listdir(root)):
                keep.append(dict(kind="other_project_raw", netsta=ns, year=-1, bytes=np.nan,
                                 reason=f"{proj}: listed only -- that project's decision"))
    C = pd.DataFrame(cand); K = pd.DataFrame(keep)
    C.to_csv(f"{OD}/DELETE_CANDIDATES_wa.csv", index=False)
    Ks = K.groupby(["kind", "netsta", "year", "reason"], dropna=False).bytes.sum().reset_index() if len(K) else K
    Ks.to_csv(f"{OD}/KEEP_wa_summary.csv", index=False)
    svi_c = f"{P.ROOT}/lfe_extension/results/refine/r9/DELETE_CANDIDATES_svi.csv"
    sv = pd.read_csv(svi_c) if os.path.exists(svi_c) else pd.DataFrame(columns=["kind", "bytes"])
    msg = ("WA PRODUCTION delete-CANDIDATE list (NOTHING DELETED; joint SVI + Washington decision by the user):\n"
           f"  WA candidates {len(C)} entries, {C.bytes.sum() / 1e12 if len(C) else 0:.3f} TB by kind "
           f"{(C.groupby('kind').bytes.sum() / 1e9).round(1).to_dict() if len(C) else {}} GB\n"
           f"  SVI R9 candidates {len(sv)} files, {sv.bytes.sum() / 1e12 if len(sv) else 0:.3f} TB (lfe_extension/results/refine/r9/)\n"
           f"  held {Ks.bytes.sum() / 1e12 if len(Ks) else 0:.3f} TB (+ other projects' raw, listed only); reasons "
           f"{(Ks.groupby('reason').bytes.sum() / 1e12).round(3).to_dict() if len(Ks) else {}}\n"
           f"  cut verified: {cut_ok}; symlink readers found: {len(R)}\n")
    open(f"{OD}/DELETE_CANDIDATES_wa.txt", "w").write(msg)
    P.log("[delete] " + msg)


if __name__ == "__main__":
    main()
