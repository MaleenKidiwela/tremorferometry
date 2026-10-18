#!/usr/bin/env python
"""WA pre-production W0 + W1 (2026-10-02).

W0  freeze: sha256 of every existing WA input the production run will use (scripts, templates, setups, confidence
    and PNET models, HF cuts, node / tremor-day / episode tables, PNSN catalogue, family catalogue, truth fractions)
    -> results/wa_pre/w0/frozen_inputs.csv (+ `verify` re-hashes and lists any change).
W1  completeness of the 21 core stations 2007-01-01 .. 2026-09-15 (production span):
      expected   = station-day inside a metadata epoch of the station (node tables, = FDSN StationXML epochs)
      raw        = a non-empty day file in any raw root on disk (download roots traces_dl*/ and the sealed
                   lfe_array_discovery/traces_sealed are NEVER read); prod_visible = in wa_common.TRACE_ROOTS
                   (minus downloads), i.e. what the WA cache builder would read
      readable   = obspy headonly read OK with >= 1 HH/BH/EH trace; cov = recorded fraction of the day of the
                   cache band's best channel (from the record headers)
      cache      = WA-chain 20 Hz cache file (lfe_extension_wa/cache20hz) present; readable = np.load of x, m, live
                   with the right shapes (the zip CRC is checked on read); svi_cache = an SVI-chain cache exists
                   (other format, not usable without an adapter -- information only)
    per station-month table + failing months (> 3 % of expected station-days missing or unreadable).
    Stage guard: > 3 % of PRESENT raw / cache files erroring -> W1_FAILED (and exit 5).
Day roles: per family-month truth_array_frac (results/wa_gonogo/truth_array_frac.csv, reused) + role.

  python scripts/pre_w01.py freeze | verify | complete --workers 24 | roles
"""
from __future__ import annotations
import os, sys, glob, hashlib, argparse, time, warnings
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pre_common as P
W = P.W

OD0 = f"{P.PRE}/w0"
OD1 = f"{P.PRE}/w1"


def frozen_list():
    S = P.STUDY; R = P.ROOT
    pats = [f"{S}/scripts/wa_*.py", f"{S}/scripts/step*.py", f"{S}/scripts/confidence_*.py", f"{S}/scripts/pnet_*.py",
            f"{S}/scripts/stack_activity*.py", f"{S}/scripts/t2b.py", f"{S}/scripts/gonogo_*.py", f"{S}/scripts/merlin_*.py",
            f"{S}/scripts/*.sh",
            f"{S}/results/templates/*", f"{S}/results/setup/*", f"{S}/results/step4b/setup/*",
            f"{S}/results/merlin_tests/t1_setup/*", f"{S}/results/wa_gonogo/setup/*",
            f"{S}/results/hf_cuts.csv", f"{S}/results/wa_gonogo/score_core21/hf_cuts.csv",
            f"{S}/results/confidence/model_info.json", f"{S}/results/confidence/candidates_confidence.parquet",
            f"{S}/results/pnet/models/*", f"{S}/results/step1_prob/model_info.json", f"{S}/results/step1_prob/candidates.parquet",
            f"{S}/results/wa_gonogo/truth_array_frac.csv", f"{S}/results/wa_gonogo/truth_array_frac_months.csv",
            f"{S}/results/wa_prod_plan/inventory_by_station.csv", f"{S}/results/WA_PROD_PLAN.md",
            f"{W.EXT}/*_bridge_nodes_strict.csv", f"{W.EXT}/*_tremor_days.csv", f"{W.EXT}/tremor_episodes.csv",
            f"{R}/catalogs/pnsn_tremor_cascadia_full.csv", f"{R}/data/catalogue/families_v1_1.csv",
            f"{R}/data/catalogue/families_v1_2.csv", f"{R}/lfe_south/scripts/cutoff_check.py",
            f"{R}/lfe_south/scripts/cut_windows.py", f"{R}/lfe_extension/scripts/svi_ref_fit.py"]
    fs = []
    for p in pats:
        fs += [f for f in sorted(glob.glob(p)) if os.path.isfile(f)]
    return sorted(set(fs))


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def _hrow(p):
    st = os.stat(p)
    return dict(path=os.path.relpath(p, P.ROOT), bytes=st.st_size, mtime=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(st.st_mtime)),
                sha256=sha(p))


def freeze(a):
    os.makedirs(OD0, exist_ok=True)
    out = f"{OD0}/frozen_inputs.csv"
    if os.path.exists(out):
        print("exists, not overwritten:", out); return
    fs = frozen_list()
    with ProcessPoolExecutor(8) as ex:
        rows = list(ex.map(_hrow, fs, chunksize=4))
    D = pd.DataFrame(rows)
    D.to_csv(out, index=False)
    P.log(f"[w0] {len(D)} files, {D.bytes.sum()/1e9:.2f} GB hashed -> {out}")


def verify(a):
    D = pd.read_csv(f"{OD0}/frozen_inputs.csv")
    rows = []
    for r in D.itertuples():
        p = f"{P.ROOT}/{r.path}"
        s = sha(p) if os.path.exists(p) else "MISSING"
        rows.append(dict(path=r.path, same=s == r.sha256, now=s))
    V = pd.DataFrame(rows)
    V.to_csv(f"{OD0}/verify_{time.strftime('%Y%m%dT%H%M')}.csv", index=False)
    print(f"[w0 verify] {len(V)} files, changed {int((~V.same).sum())}")
    P.log(f"[w0 verify] {len(V)} files, changed {int((~V.same).sum())}: {V[~V.same].path.tolist()[:10]}")


# ------------------------------------------------------------------------------------------------ W1
BANDS = ["HH", "BH", "EH"]


def day_job(args):
    """FAST check (2026-10-02 22:10Z revision: a full headonly read of every file ran at ~140 MB/s, ~6 h for 2.9 TB).
    raw: file present, non-empty, FIRST record parses (obspy get_record_information: channel, rate, record length),
    file size is a whole number of records.  The FULL read of every raw day file is done by W2 raw (pre_w2.py) for all
    station-days not already fully read by SVI R2; W1's final readability combines both (w1_readability_final.csv).
    cache: npz zip directory + the small members m / live read (zip CRC checked) + x array header (shape, dtype)."""
    ns, ds = args
    warnings.filterwarnings("ignore")
    r = dict(netsta=ns, date=ds)
    files = P.raw_files(ns, ds)
    r["raw"] = bool(files)
    r["prod_visible"] = P.prod_raw_path(ns, ds) is not None
    r["raw_path"] = files[0] if files else ""
    r["n_raw_roots"] = len(files)
    r["readable"] = False; r["err"] = ""
    if files:
        from obspy.io.mseed.util import get_record_information
        try:
            rp = os.path.realpath(files[0]); sz = os.path.getsize(rp)
            info = get_record_information(rp)
            rl = int(info.get("record_length", 0))
            r["rec_len"] = rl; r["fs_first"] = float(info.get("samp_rate", np.nan)); r["chan_first"] = str(info.get("channel", ""))
            r["bytes"] = sz
            r["whole_records"] = bool(rl > 0 and sz % rl == 0)
            r["readable"] = bool(rl > 0 and r["chan_first"][:2] in BANDS)
            if not r["readable"]:
                r["err"] = f"first record channel {r['chan_first']} rec_len {rl}"
        except Exception as e:
            r["err"] = f"{type(e).__name__}: {str(e)[:60]}"
    ts = pd.Timestamp(ds); tail = f"{ns}/{ts.year}/{int(ts.dayofyear):03d}.npz"
    c = f"{P.WA_CACHE}/{tail}"
    r["cache"] = os.path.exists(c)
    r["cache_ok"] = False
    if r["cache"]:
        try:
            import zipfile
            with zipfile.ZipFile(c) as zf:
                names = set(zf.namelist())
                with zf.open("x.npy") as fx:
                    ver = np.lib.format.read_magic(fx)
                    shp, _, dt = (np.lib.format.read_array_header_1_0 if ver == (1, 0) else np.lib.format.read_array_header_2_0)(fx)
                m = np.load(zf.open("m.npy")); lv = np.load(zf.open("live.npy"))
                node = str(np.load(zf.open("node.npy"))) if "node.npy" in names else ""
            r["cache_ok"] = bool(tuple(shp) == (3, W.NDAY) and m.shape == (W.NDAY,) and lv.shape == (3,))
            r["cache_node"] = node; r["cache_valid"] = float(m.mean())
            if not r["cache_ok"]:
                r["cache_err"] = f"shape {shp}"
        except Exception as e:
            r["cache_err"] = f"{type(e).__name__}: {str(e)[:60]}"
    r["svi_cache"] = any(os.path.exists(f"{s_}/{tail}") for s_ in P.SVI_CACHES)
    return r


def complete(a):
    os.makedirs(OD1, exist_ok=True)
    out = f"{OD1}/w1_station_days.parquet"
    if os.path.exists(out):
        print("exists:", out); X = pd.read_parquet(out)
    else:
        N = P.nodes(); core = P.core_list()
        jobs = []
        allday = [d.strftime("%Y-%m-%d") for d in pd.date_range(P.T0, P.T1, freq="D")]
        EXP = {}
        for ns in core:
            e = P.expected_days(ns, N); EXP[ns] = e
            jobs += [(ns, d) for d in allday]          # every day: also catches raw OUTSIDE the metadata epochs
        P.log(f"[w1] {len(jobs)} station-days to check, {a.workers} workers")
        t0 = time.time()
        with ProcessPoolExecutor(a.workers) as ex:
            rows = list(ex.map(day_job, jobs, chunksize=64))
        X = pd.DataFrame(rows)
        X["expected"] = [d in EXP[ns] for ns, d in zip(X.netsta, X.date)]
        X.to_parquet(out, index=False)
        P.log(f"[w1] done {time.time()-t0:.0f}s")
    X["month"] = X.date.str[:7]
    X["cov50"] = X.readable
    E = X[X.expected]
    M = E.groupby(["netsta", "month"]).agg(expected=("date", "size"), raw=("raw", "sum"), prod_visible=("prod_visible", "sum"),
                                           readable=("readable", "sum"), cov50=("cov50", "sum"),
                                           raw_unreadable=("raw", lambda v: int((v & ~E.loc[v.index, "readable"]).sum())),
                                           cache=("cache", "sum"), cache_ok=("cache_ok", "sum"),
                                           cache_err=("cache", lambda v: int((v & ~E.loc[v.index, "cache_ok"]).sum())),
                                           svi_cache=("svi_cache", "sum")).reset_index()
    M["missing"] = M.expected - M.readable
    M["missing_frac"] = M.missing / M.expected
    M["not_prod_visible"] = M.readable - E.groupby(["netsta", "month"]).apply(lambda g: int((g.readable & g.prod_visible).sum())).values
    M["fail_3pct"] = M.missing_frac > 0.03
    M.to_csv(f"{OD1}/w1_station_month.csv", index=False)
    F = M[M.fail_3pct].copy()
    F["why"] = np.where(F.raw < F.expected, "no raw file on disk", "unreadable raw")
    F.to_csv(f"{OD1}/w1_failing_months.csv", index=False)
    # outside metadata epochs but raw present (information)
    O = X[~X.expected & X.raw].groupby(["netsta", X.date.str[:4]]).size().rename("n").reset_index()
    O.to_csv(f"{OD1}/w1_raw_outside_epochs.csv", index=False)
    # per station summary
    S = M.groupby("netsta").agg(expected=("expected", "sum"), readable=("readable", "sum"), raw_unreadable=("raw_unreadable", "sum"),
                                cache=("cache", "sum"), cache_err=("cache_err", "sum"), svi_cache=("svi_cache", "sum"),
                                n_months=("month", "size"), n_fail_months=("fail_3pct", "sum")).reset_index()
    S["missing_frac"] = 1 - S.readable / S.expected
    S.to_csv(f"{OD1}/w1_station_summary.csv", index=False)
    # era summary
    M["era"] = np.where(M.month < "2012-01", "2007-2011", np.where(M.month < "2014-01", "2012-2013", "2014-2026"))
    ER = M.groupby("era")[["expected", "readable", "raw_unreadable", "cache", "cache_err"]].sum()
    ER["missing_frac"] = 1 - ER.readable / ER.expected
    ER.to_csv(f"{OD1}/w1_era_summary.csv")
    pres = int(E.raw.sum()); bad = int((E.raw & ~E.readable).sum())
    cpres = int(E.cache.sum()); cbad = int((E.cache & ~E.cache_ok).sum())
    msg = (f"[w1] expected {len(E)}; raw present {pres}, unreadable {bad} ({bad/max(pres,1):.3%}); WA cache present {cpres}, "
           f"error {cbad} ({cbad/max(cpres,1):.3%}); missing {int(M.missing.sum())} ({M.missing.sum()/len(E):.2%}); "
           f"failing station-months {int(M.fail_3pct.sum())} of {len(M)}")
    P.log(msg); print(msg); print(ER.to_string()); print(S.to_string())
    if bad / max(pres, 1) > 0.03 or cbad / max(cpres, 1) > 0.03:
        P.marker("W1_FAILED", msg); sys.exit(5)


# ------------------------------------------------------------------------------------------------ roles
def roles(a):
    os.makedirs(OD1, exist_ok=True)
    T = pd.read_csv(f"{P.STUDY}/results/wa_gonogo/truth_array_frac.csv")

    def role(m):
        if m == "2007-01":
            return "sealed_one_look_headline"
        if m == "2009-05":
            return "sealed_one_look_semi_independent"
        if m in P.EXCLUDED_MONTHS:
            return "excluded_no_array_stations"
        return "fit"
    T["role"] = T.month.map(role)
    note = {"2008-05": "build half (V2008) = in-sample (template days, never fitted or scored); test half = fit",
            "2010-08": "build half (V2010) = in-sample; test half = fit",
            "2011-08": "A2011: fit; in-sample for the stations whose templates come from tpl_step4b_a (D1)",
            "2011-04": "truth made only by permanent stations (0 array stations): not fit truth",
            "2009-05": "sealed, one look; ~half of the truth stations are permanent (semi-independent)",
            "2007-01": "sealed, one look, headline (array-majority truth)"}
    T["note"] = T.month.map(note).fillna("")
    out = f"{OD1}/family_month_roles.csv"
    if not os.path.exists(out):
        T.to_csv(out, index=False)
    # Olympic A / Arm B truth: 2010-07..09 only (no truth_array_frac computed; reported out-of-fold only)
    tr = pd.read_parquet(f"{P.STUDY}/results/merlin_tests/t1_setup/truth.parquet")
    F = pd.read_csv(f"{P.STUDY}/results/merlin_tests/t1_setup/families.csv")
    tr = tr.merge(F[["fid", "region"]], on="fid", how="left") if "region" not in tr else tr
    tr["month"] = tr.date.str[:7]
    O = tr[tr.region.isin(["Olympic_2010", "Olympic_2010_armB"])].groupby(["fid", "region", "month"]).size().rename("n_truth").reset_index()
    O["role"] = "fit_out_of_fold_only (2010 is the only truth; no sealed check)"
    O.to_csv(f"{OD1}/family_month_roles_olympic.csv", index=False)
    # day roles: every day 2007-01-01 .. 2026-09-15
    scan = pd.read_csv(f"{P.STUDY}/results/merlin_tests/t1_setup/scan_days.csv")
    tm = set(T.month)
    rows = []
    for d in pd.date_range(P.T0, P.T1, freq="D"):
        ds = d.strftime("%Y-%m-%d"); m = ds[:7]
        if m in tm or m in ("2010-07", "2010-09"):
            r = role(m) if m in tm else "fit"
        else:
            r = "other"
        rows.append(dict(date=ds, month=m, role=r))
    D = pd.DataFrame(rows)
    sd = scan.groupby("date").det_source.apply(lambda v: ";".join(sorted(set(x.split("ab_royer/")[-1].rsplit("/", 1)[0] for x in v))))
    D["truth_scanned_by"] = D.date.map(sd).fillna("")
    D["half_rule"] = np.where(D.month.isin(["2008-05", "2010-08"]), "build half in-sample / test half fit", "")
    D.loc[D.date.isin(["2010-07-14", "2010-07-15", "2010-07-16", "2010-07-17", "2010-07-18"]), "half_rule"] = "V2010 quiet days: build half in-sample / test half fit"
    out = f"{OD1}/day_roles.csv"
    if not os.path.exists(out):
        D.to_csv(out, index=False)
    print(T.groupby(["month", "role"]).agg(n_fam=("fid", "size"), n_truth=("n_truth", "sum"),
                                         med_array_frac=("truth_array_frac", "median")).to_string())
    print(D.role.value_counts().to_string())
    P.log(f"[roles] {len(T)} family-months, {len(O)} Olympic family-months, {len(D)} days: {D.role.value_counts().to_dict()}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["freeze", "verify", "complete", "roles"])
    ap.add_argument("--workers", type=int, default=24)
    a = ap.parse_args()
    dict(freeze=freeze, verify=verify, complete=complete, roles=roles)[a.what](a)
