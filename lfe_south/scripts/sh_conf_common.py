#!/usr/bin/env python
"""SH_CONF (2026-10-01) -- shared definitions for the per-repeat confidence of Shelly (2025)'s 27 Mendocino families.

Two candidate runs (every >= 8 x MAD 4-s-bin maximum, nothing cut):
  v3b  the SH_V3B it0 run on the 385 Shelly2025 TEST days 2018-2024 (results/sh_v3b/mf_it0_test.parquet; v3b prep)
  fwd  the forward run 2023-12-01 .. 2026-09-28 (results/sh_ext2026/candidates_raw8x.parquet; segments cont / A / NK / B / C,
       each with its own template set; sh_ext2026 prep)
Each run also gets ONE site-rolled scrambled control (series 'c1'): the same detector with every <= 3-km site's
correlation series rolled by its own random shift of the same day (sh_conf_mf.py).

Day roles (leakage, results/SH_CONF_RESULTS.md):
  build_any  the day is a v3b template-build day of ANY family (stage2_days_build.json)  -> in-sample, never fitted/evaluated
  fit        v3b run, 2018-2023, not build_any                                      -> fit (5-fold grouped-by-day out-of-fold)
  ho2024     v3b run, 2024, not build_any                                           -> independent evaluation 1
  fwd_ovl    fwd run, day < 2024-12-01, not build_any, not a fit day                 -> independent evaluation 2 (his catalogue exists)
  fwd_ext    fwd run, day >= 2024-12-01                                             -> application only
"""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t_common import ROOT, WS, RES, CACHE, FS, TR, npz_path, day_path, log

OUT = os.environ.get("SH_CONF_DIR", f"{RES}/sh_conf")
V3B = f"{RES}/sh_v3b"; EXT = f"{RES}/sh_ext2026"
SRC = "Shelly2025"; FAMS = [str(i) for i in range(1, 28)]
STAS = ["BK.BLAS", "BK.BRIC", "BK.GHOP", "BK.HALS", "BK.JEWT", "BK.PRDS", "BK.SPRL", "BK.RVIT", "BK.HUNT",
        "NC.KMPB", "NC.KMR", "NC.KPR", "NC.KBU", "PB.B049", "PB.B932", "PB.B934", "PB.B933"]
T_HIS0 = pd.Timestamp("2018-01-01"); T_HIS1 = pd.Timestamp("2024-12-01")
KMR_SW1 = pd.Timestamp("2025-09-17"); KMR_SW2 = pd.Timestamp("2026-03-12")
TMIN = 8.0            # candidate tier (x daily MAD)
GAP_S = 4.0           # the MF's fixed 4-s bins
NRAND = 1000          # random origins per family-day (site noise level)
MATCH_S = 2.0
PNET_HO = ["BK.JEWT", "NC.KBU", "PB.B933"]     # stations never used in the PNET fine-tune
ERR_MAX = 0.03        # a stage fails if > 3 % of the station-days it read errored
TPL = {"v3b": f"{V3B}/tpl_it0", "cont": f"{V3B}/tpl_it0", "A": f"{V3B}/tpl_it0", "NK": f"{EXT}/tpl_noKMR",
       "B": f"{EXT}/tpl_B", "C": f"{EXT}/tpl_C"}
PREPDIR = {"v3b": f"{V3B}/prep", "fwd": f"{EXT}/prep"}


def read_days(path, source=SRC):
    dd = json.load(open(path)); s = set()
    for k, v in dd.items():
        if k.split("|")[0] == source:
            s |= set(v)
    return sorted(s)


def fam_days(path, source=SRC):
    dd = json.load(open(path))
    return {k.split("|")[1]: set(v) for k, v in dd.items() if k.split("|")[0] == source}


def day_table():
    """one row per (run, day): seg, template set, prep dir, day flags."""
    test = read_days(f"{RES}/aa_shelly/stage2_days_test.json")
    BF = fam_days(f"{RES}/aa_shelly/stage2_days_build.json")
    build_any = set().union(*BF.values())
    rows = [dict(run="v3b", date=d, seg="v3b") for d in test]
    segs = {}
    for seg in ("cont", "A", "B", "C"):
        for d in read_days(f"{EXT}/days/{seg}.json"):
            segs[d] = seg
    for d in (KMR_SW1, KMR_SW2):
        segs[str(d.date())] = "NK"
    proc = pd.read_csv(f"{EXT}/days_processed.csv").day.astype(str).str[:10]
    for d in sorted(proc):
        rows.append(dict(run="fwd", date=d, seg=segs.get(d, "?")))
    D = pd.DataFrame(rows)
    if (D.seg == "?").any():
        raise SystemExit(f"forward days without a segment: {D[D.seg == '?'].date.tolist()[:5]}")
    D["tpl"] = D.seg.map(TPL); D["prep"] = D.run.map(PREPDIR)
    D["build_any"] = D.date.isin(build_any)
    D["year"] = D.date.str[:4].astype(int)
    D["his_cov"] = pd.to_datetime(D.date) < T_HIS1
    fit_days = set(D[(D.run == "v3b") & (D.year <= 2023) & ~D.build_any].date)
    D["role"] = "fwd_ext"
    D.loc[(D.run == "v3b") & D.build_any, "role"] = "insample"
    D.loc[(D.run == "v3b") & (D.year <= 2023) & ~D.build_any, "role"] = "fit"
    D.loc[(D.run == "v3b") & (D.year == 2024) & ~D.build_any, "role"] = "ho2024"
    ov = (D.run == "fwd") & D.his_cov
    D.loc[ov & D.build_any, "role"] = "insample"
    D.loc[ov & ~D.build_any & D.date.isin(fit_days), "role"] = "fwd_fitday"
    D.loc[ov & ~D.build_any & ~D.date.isin(fit_days), "role"] = "fwd_ovl"
    return D


def build_fam_days():
    return fam_days(f"{RES}/aa_shelly/stage2_days_build.json")


def load_tpl(tag):
    M = pd.read_csv(tag + ".csv"); M = M[M.source == SRC].copy(); M["fam"] = M.fam.astype(str)
    Z = np.load(tag + ".npz"); T = {k: Z[k].astype(np.float64) for k in M.key if k in Z}
    return M, T


def sites():
    """netsta -> site id (<= 3 km complete-chain clustering, sh_v2_common.neighbour_sites)."""
    import sh_v2_common as C2
    C = C2.coords()
    sid, _ = C2.neighbour_sites(STAS, C)
    return {s: int(k) for s, k in zip(STAS, sid)}


def families():
    F = pd.read_csv(f"{RES}/families_south.csv")
    F = F[F.source == SRC].copy(); F["fam"] = F.fam.astype(int).astype(str)
    return F


def fam_groups(km=3.0):
    """fam -> list of families within <= km epicentral distance (incl. itself)"""
    import sh_v2_common as C2
    F = families()
    la, lo = F.lat.to_numpy(float), F.lon.to_numpy(float)
    d = C2.epi_km(la[:, None], lo[:, None], la[None], lo[None])
    fl = F.fam.tolist()
    return {f: [fl[j] for j in np.nonzero(d[i] <= km)[0]] for i, f in enumerate(fl)}


def his_catalogue():
    from sh_v3_common import his_catalogue as h
    S = h(); S["t_epoch"] = S.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
    return S


def candidates(run, seg=None):
    if run == "v3b":
        D = pd.read_parquet(f"{V3B}/mf_it0_test.parquet")
        D = D[D.net_mad >= TMIN].copy(); D["seg"] = "v3b"
    else:
        D = pd.read_parquet(f"{EXT}/candidates_raw8x.parquet")
        D["t_epoch"] = D.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
    D["fam"] = D.fam.astype(str)
    D["date"] = D.t.dt.floor("D").dt.strftime("%Y-%m-%d")
    if seg is not None:
        D = D[D.seg == seg]
    return D.sort_values(["date", "fam", "t_epoch"]).reset_index(drop=True)


def site_shifts(date, site_ids, ns_samples):
    """deterministic per-site roll (samples) for the scrambled control of `date`: distinct shifts in
    [10 min, 1 day - 10 min], >= 2 min apart."""
    rng = np.random.default_rng(20261001 + int(date.replace("-", "")))
    lo, hi = int(600 * FS), int(ns_samples - 600 * FS)
    out = {}
    u = sorted(set(site_ids))
    while True:
        v = rng.integers(lo, hi, len(u))
        if len(u) < 2 or np.min(np.diff(np.sort(v))) >= 120 * FS:
            break
    for k, s in zip(u, v):
        out[int(k)] = int(s)
    return out


def rand_origins(date, n_samples, n=NRAND):
    rng = np.random.default_rng(1000 + int(date.replace("-", "")))
    return rng.integers(0, n_samples, n)


def s_index(i, sh, shift, nw, pre_n):
    """station data-array index of the S pick for published-axis origins i (data: shift=0, invalid outside the
    array -> -1; control: rolled modulo nw exactly as the scrambled-control MF built it)."""
    i = np.asarray(i, np.int64)
    if shift:
        return ((i + sh + shift) % nw) + pre_n
    s = i + sh
    return np.where((s >= 0) & (s < nw), s + pre_n, -1)


def write_status(stage, tag, rows):
    os.makedirs(f"{OUT}/status", exist_ok=True)
    p = f"{OUT}/status/{stage}_{tag}.csv"
    pd.DataFrame(rows, columns=["run", "date", "netsta", "status", "msg"]).to_csv(p, index=False)


def guard(stage):
    """fail (exit 5) when > ERR_MAX of the station-days read by `stage` errored."""
    import glob
    fs = sorted(glob.glob(f"{OUT}/status/{stage}_*.csv"))
    if not fs:
        log(f"[guard:{stage}] no status files"); sys.exit(5)
    S = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    S = S.drop_duplicates(["run", "date", "netsta"], keep="last")
    n_ok = int((S.status == "ok").sum()); n_err = int((S.status == "error").sum())
    frac = n_err / max(n_ok + n_err, 1)
    msg = (f"[guard:{stage}] station-days read ok {n_ok}, error {n_err} ({frac:.2%}); other {S.status.value_counts().to_dict()}")
    log(msg)
    open(f"{OUT}/status/GUARD_{stage}.txt", "w").write(msg + "\n")
    if frac > ERR_MAX:
        log(f"[guard:{stage}] FAIL: error fraction {frac:.2%} > {ERR_MAX:.0%}; examples:\n"
            + S[S.status == "error"].head(10).to_string())
        sys.exit(5)


def mfcheck():
    """every run-day has a feature file; the re-run reproduces every stored candidate (net_mad / ccsum)."""
    D = pd.read_csv(f"{OUT}/day_table.csv", dtype={"date": str})
    if len(sys.argv) > 2: D = D[D.date.isin(sys.argv[2].split(","))]
    miss, bad, n_st, n_dat = [], [], 0, 0
    for r in D.itertuples():
        p = f"{OUT}/feat/{r.run}/{r.date}.json"
        if not os.path.exists(p): miss.append((r.run, r.date)); continue
        m = json.load(open(p)); n_st += m["n_cand_stored"]; n_dat += m["n_data"]
        if m["max_err_mad"] > 1e-3 or m["max_err_ccsum"] > 1e-4 or m["n_data"] != m["n_cand_stored"]:
            bad.append((r.run, r.date, m["max_err_mad"], m["n_data"], m["n_cand_stored"]))
    msg = (f"[mfcheck] run-days {len(D)}, missing {len(miss)}, not reproduced {len(bad)}; candidates found {n_dat} of {n_st} stored")
    log(msg + (f"; e.g. missing {miss[:5]} bad {bad[:5]}" if (miss or bad) else ""))
    open(f"{OUT}/status/MFCHECK.txt", "w").write(msg + "\n")
    if len(miss) + len(bad) > ERR_MAX * len(D): sys.exit(7)


if __name__ == "__main__":
    if sys.argv[1] == "mfcheck":
        mfcheck()
    elif sys.argv[1] == "guard":
        guard(sys.argv[2])
    elif sys.argv[1] == "days":
        D = day_table(); os.makedirs(OUT, exist_ok=True); D.to_csv(f"{OUT}/day_table.csv", index=False)
        log(f"[days] {D.groupby(['run', 'role']).size().to_dict()}")
