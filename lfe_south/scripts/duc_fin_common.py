#!/usr/bin/env python
"""DUC_FIN (2026-10-01) -- shared definitions for finishing the Ducellier & Creager (2022) extension on the data
already on disk (2008-01 .. 2012-07): 2011 diagnosis, re-templating of the new 2011-10..12 epochs, the feature +
scrambled-control matched filter, repeat windows for dv/v, the per-repeat confidence and the delete-candidate list.

Everything is written under OUT = $DUC_FIN_DIR (default results/duc_fin); smoke tests point DUC_FIN_DIR at a
scratch folder.  Nothing under results/duc_fix/ is modified: duc_fix's mf1/, prep/, station_epochs.csv and
cache_duc/ are READ only.

Day roles (fixed BEFORE any fit, see day_table()):
  ho_month    6 held-out months spread over 2008-2010 (HO_MONTHS) -- never fitted, independent check 1
  pnet_train  every 10th remaining 2008-2010 day -- the PNET fine-tune windows come only from these days, so they are
              in-sample for PNET: never fitted, never evaluated (P is still written)
  fit         the other 2008-2010 days (5 outer folds grouped by day, out-of-fold P)
  ind2011a    2011-01..04 (his catalogue complete) -- independent check 2
  gap2011     2011-05..12 (his catalogue lost stations; results/duc_fin/D2011.md) -- evaluated separately, not fitted
  ext2012     2012-01..07 (no catalogue) -- application only
"""
import os, sys, json, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t_common import ROOT, WS, RES, TR, FS, log, day_path
import duc_fix_common as DC

OUT = os.environ.get("DUC_FIN_DIR", f"{RES}/duc_fin")
FIX = f"{RES}/duc_fix"
T0, T1 = "2008-01-01", "2012-08-01"                       # [T0, T1)
MONTHS = [str(m) for m in pd.period_range("2008-01", "2012-07", freq="M")]
HO_MONTHS = ["2008-02", "2008-09", "2009-03", "2009-10", "2010-05", "2010-09"]   # fixed 2026-10-01 before fitting
PNET_HO = ["PB.B935", "NC.KHMB", "NC.LBP"]                 # never used in the PNET fine-tune (LBP = 1-component)
ERR_MAX = 0.03
NRAND = 1000
MATCH_S = 2.0
TMIN = 8.0
CUT_OUT = os.environ.get("DUC_FIN_CUT_DIR", f"{RES}/cut/duc")


def day_table():
    days = pd.date_range(T0, T1, freq="D", inclusive="left")
    D = pd.DataFrame(dict(date=[str(d.date()) for d in days]))
    D["month"] = D.date.str[:7]; D["year"] = D.date.str[:4].astype(int)
    D["role"] = "ext2012"
    D.loc[D.year <= 2010, "role"] = "fit"
    D.loc[D.month.isin(HO_MONTHS), "role"] = "ho_month"
    cand = D.index[(D.role == "fit")]
    D.loc[cand[5::10], "role"] = "pnet_train"
    D.loc[(D.year == 2011) & (D.month <= "2011-04"), "role"] = "ind2011a"
    D.loc[(D.year == 2011) & (D.month >= "2011-05"), "role"] = "gap2011"
    return D


def netof():
    H = DC.his_rows(); return dict(zip(H.sta, H.net))


def stations():
    fp = DC.families(); used = set(",".join(fp.stations).split(",")); no = netof()
    return sorted(f"{no[s]}.{s}" for s in used)


def fam_locs():
    """his family locations (Plourde 2015 templates_list.txt: name, family, lat, lon, depth)."""
    p = f"{DC.DUC}/data/Plourde_2015/templates_list.txt"
    rows = []
    for ln in open(p).read().splitlines()[1:]:
        q = ln.split()
        if len(q) >= 5: rows.append((q[0], float(q[2]), float(q[3]), float(q[4])))
    return pd.DataFrame(rows, columns=["fam", "lat", "lon", "depth"])


def sites():
    """netsta -> site id: <= 3 km clustering (duc_fix_clean.neighbour_sites).  His 47 stations: no two within 3 km."""
    import duc_fix_clean as CL
    C = CL.coords(); st = stations()
    sid, _ = CL.neighbour_sites(st, C)
    return {s: int(k) for s, k in zip(st, sid)}


def fam_groups(km=3.0):
    import duc_fix_clean as CL
    F = fam_locs(); la, lo = F.lat.to_numpy(), F.lon.to_numpy()
    d = CL.epi_km(la[:, None], lo[:, None], la[None], lo[None])
    fl = F.fam.tolist()
    return {f: [fl[j] for j in np.nonzero(d[i] <= km)[0]] for i, f in enumerate(fl)}


def site_shifts(date, site_ids, ns_samples):
    """deterministic per-site roll (samples) for the scrambled control of `date`: distinct shifts in
    [10 min, 1 day - 10 min], >= 2 min apart (sh_conf_common.site_shifts, own seed)."""
    rng = np.random.default_rng(20261002 + int(date.replace("-", "")))
    lo, hi = int(600 * FS), int(ns_samples - 600 * FS)
    u = sorted(set(int(s) for s in site_ids))
    while True:
        v = rng.integers(lo, hi, len(u))
        if len(u) < 2 or np.min(np.diff(np.sort(v))) >= 120 * FS:
            break
    return {k: int(s) for k, s in zip(u, v)}


def rand_origins(date, n_samples, n=NRAND):
    rng = np.random.default_rng(2000 + int(date.replace("-", "")))
    return rng.integers(0, n_samples, n)


# ------------------------------------------------------------------------------------------- his catalogue
def his_catalogue():
    """his 'net' catalogue (every 8 x MAD detection, all 66 families) with his catalogue rule flag hp
    (cc x nch >= threshold_perm, threshold 0 = eliminated family -> never in his catalogue)."""
    H = DC.his_catalogue(); thr = DC.thresholds()
    H["thr"] = H.fam.map(thr).fillna(0.0)
    H["hp"] = (H.thr > 0) & (H.cc * H.nch >= H.thr)
    H["t_epoch"] = H.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
    return H


def match_one_to_one(ta, fa, tb, fb, tol=MATCH_S, same_fam=True):
    """greedy one-to-one matching (closest pairs first) of events a (times ta, families fa) to events b.
    Returns (ia, ib, dt) arrays of matched index pairs."""
    ta = np.asarray(ta, float); tb = np.asarray(tb, float)
    ob = np.argsort(tb); tbs = tb[ob]
    lo = np.searchsorted(tbs, ta - tol, "left"); hi = np.searchsorted(tbs, ta + tol, "right")
    n = hi - lo
    ia = np.repeat(np.arange(len(ta)), n)
    if not len(ia): return np.zeros(0, int), np.zeros(0, int), np.zeros(0)
    off = np.arange(len(ia)) - np.repeat(np.cumsum(n) - n, n)
    jb = ob[np.repeat(lo, n) + off]
    dt = np.abs(tb[jb] - ta[ia])
    if same_fam:
        ok = np.asarray(fa)[ia] == np.asarray(fb)[jb]
        ia, jb, dt = ia[ok], jb[ok], dt[ok]
    o = np.argsort(dt, kind="stable"); ia, jb, dt = ia[o], jb[o], dt[o]
    ua = np.zeros(len(ta), bool); ub = np.zeros(len(tb), bool); keep = np.zeros(len(ia), bool)
    for k in range(len(ia)):
        a, b = ia[k], jb[k]
        if ua[a] or ub[b]: continue
        ua[a] = ub[b] = True; keep[k] = True
    return ia[keep], jb[keep], dt[keep]


def nearest(tq, tref, fq=None, fref=None):
    """distance (s) from each query to the nearest reference event (same family if fq/fref given)."""
    tq = np.asarray(tq, float); out = np.full(len(tq), np.inf)
    if fq is None:
        r = np.sort(np.asarray(tref, float))
        if not len(r): return out
        j = np.clip(np.searchsorted(r, tq), 1, max(len(r) - 1, 1))
        return np.minimum(np.abs(r[j - 1] - tq), np.abs(r[np.minimum(j, len(r) - 1)] - tq))
    fq = np.asarray(fq); fref = np.asarray(fref); tref = np.asarray(tref, float)
    for f in np.unique(fq):
        m = fq == f; r = np.sort(tref[fref == f])
        if not len(r): continue
        j = np.clip(np.searchsorted(r, tq[m]), 1, max(len(r) - 1, 1))
        out[m] = np.minimum(np.abs(r[j - 1] - tq[m]), np.abs(r[np.minimum(j, len(r) - 1)] - tq[m]))
    return out


# ------------------------------------------------------------------------------------------- arrival offsets
def crh(path):
    v, z = [], []
    for ln in open(path).read().splitlines()[1:]:
        if not ln.strip(): continue
        v.append(float(ln[:5])); z.append(float(ln[5:]))
    return np.array(v), np.array(z)


def ts_direct(depth, xkm, v, ztop):
    """direct-S travel time (s) from a source at `depth` km to a surface receiver at epicentral distance xkm
    in the 1-D layered model (layer tops ztop, velocities v): ray shooting by bisection on the ray parameter."""
    th = []
    for i in range(len(v)):
        a = ztop[i]; b = ztop[i + 1] if i + 1 < len(v) else 1e9
        if a >= depth: break
        th.append((min(b, depth) - a, v[i]))
    h = np.array([t[0] for t in th]); vv = np.array([t[1] for t in th])
    if xkm < 1e-3: return float(np.sum(h / vv))
    pmax = 1.0 / vv.max() * (1 - 1e-9)
    lo, hi = 0.0, pmax
    for _ in range(80):
        p = 0.5 * (lo + hi)
        X = np.sum(h * p * vv / np.sqrt(1 - (p * vv) ** 2))
        if X < xkm: lo = p
        else: hi = p
    p = 0.5 * (lo + hi)
    return float(np.sum(h / (vv * np.sqrt(1 - (p * vv) ** 2))))


def arrival_offsets(all_stations=None):
    """S-arrival offset (s after the detection time = his template first sample) for every (family, station).
    Family template stations: the peak of HIS template's energy envelope (1-s smoothing of x^2, each channel
    normalised, summed over the station's channels) -> source 'template'.  Every other station: direct-S time in the USGS CVM v1.7
    mendo_MTJ 1-D S model (data/velocity_models/usgs_cascadia_v1.7/profiles/mendo_MTJ_S.crh) from his family location
    (Plourde 2015 list) + the family's median (template peak - predicted) correction -> source 'moveout_cvm17'.
    Returns DataFrame fam, netsta, offset_s, offset_source, pred_s, resid_s (template stations)."""
    import duc_fix_clean as CL
    v, z = crh(f"{ROOT}/data/velocity_models/usgs_cascadia_v1.7/profiles/mendo_MTJ_S.crh")
    T = DC.his_templates(); no = netof(); C = CL.coords(); F = fam_locs().set_index("fam")
    stas = sorted(set(all_stations or []) | set(stations()))
    stas = [s for s in stas if s in C]
    rows = []
    for f in sorted(T):
        if f not in F.index: continue
        la, lo, dep = F.loc[f, ["lat", "lon", "depth"]]
        en = {}
        for (sta, cha), tp in T[f].items():
            e2 = tp * tp / max(float(np.mean(tp * tp)), 1e-30)
            en.setdefault(sta, []).append(np.convolve(e2, np.ones(20) / 20, "same"))
        pred = {ns: ts_direct(dep, float(CL.epi_km(la, lo, *C[ns])), v, z) for ns in stas}
        res = []
        tpl = {}
        for sta, es in en.items():
            ns = f"{no.get(sta, 'NC')}.{sta}"
            tpl[ns] = float(np.argmax(np.sum(es, 0))) / FS
            if ns in pred: res.append(tpl[ns] - pred[ns])
        corr = float(np.median(res)) if res else 0.0
        for ns in stas:
            if ns in tpl:
                rows.append((f, ns, tpl[ns], "template", pred[ns] + corr, tpl[ns] - pred[ns] - corr))
            else:
                rows.append((f, ns, pred[ns] + corr, "moveout_cvm17", pred[ns] + corr, np.nan))
    return pd.DataFrame(rows, columns=["fam", "netsta", "offset_s", "offset_source", "pred_s", "resid_s"])


# ------------------------------------------------------------------------------------------- status + guard
def write_status(stage, tag, rows):
    os.makedirs(f"{OUT}/status", exist_ok=True)
    pd.DataFrame(rows, columns=["date", "netsta", "status", "msg"]).to_csv(f"{OUT}/status/{stage}_{tag}.csv", index=False)


def guard(stage, out=None):
    """exit 5 when > ERR_MAX of the station-days that stage was expected to read errored or were silently missing
    (status 'error' or 'missing'); legitimate absences (no data at the data centre, no epoch) are 'nodata'."""
    out = out or OUT
    fs = sorted(glob.glob(f"{out}/status/{stage}_*.csv"))
    if not fs:
        log(f"[guard:{stage}] no status files"); sys.exit(5)
    S = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True).drop_duplicates(["date", "netsta"], keep="last")
    bad = S.status.isin(["error", "missing"])
    n_ok = int((S.status == "ok").sum()); n_bad = int(bad.sum())
    frac = n_bad / max(n_ok + n_bad, 1)
    msg = (f"[guard:{stage}] station-days ok {n_ok}, error/missing {n_bad} ({frac:.2%}); all {S.status.value_counts().to_dict()}")
    log(msg)
    os.makedirs(f"{out}/status", exist_ok=True)
    open(f"{out}/status/GUARD_{stage}.txt", "w").write(msg + "\n")
    if frac > ERR_MAX:
        log(f"[guard:{stage}] FAIL: {frac:.2%} > {ERR_MAX:.0%}; examples:\n" + S[bad].head(15).to_string())
        sys.exit(5)


def plan_status():
    """(sta, day) -> duc_fix plan cache status for 2008-01..2012-07 (built/cached = data on disk)."""
    out = {}
    for m in MONTHS:
        p = f"{FIX}/data/plan_{m}.csv"
        if not os.path.exists(p): continue
        P = pd.read_csv(p, dtype={"loc": str}, usecols=["net", "sta", "day", "cache"])
        for r in P.itertuples():
            out[(f"{r.net}.{r.sta}", r.day)] = str(r.cache).split(":")[0]
    return out


if __name__ == "__main__":
    if sys.argv[1] == "days":
        os.makedirs(OUT, exist_ok=True); D = day_table(); D.to_csv(f"{OUT}/day_table.csv", index=False)
        log(f"[days] {D.role.value_counts().to_dict()}")
    elif sys.argv[1] == "guard":
        guard(sys.argv[2])
    elif sys.argv[1] == "offsets":
        os.makedirs(OUT, exist_ok=True); O = arrival_offsets(); O.to_csv(f"{OUT}/arrival_offsets.csv", index=False)
        r = O.resid_s.dropna()
        log(f"[offsets] {len(O)} family-stations ({(O.offset_source == 'template').sum()} template); template-station "
            f"residual vs CVM1.7 moveout: median |r| {r.abs().median():.2f} s, p90 {r.abs().quantile(0.9):.2f} s")
