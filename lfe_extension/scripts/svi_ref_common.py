#!/usr/bin/env python
"""SVI_REF (2026-10-02) -- shared definitions for the southern-Vancouver-Island refinement queue
(results/SVI_REFINE_PLAN.md as amended by Merlin's review, coordinator brief 2026-10-02).

Nothing of catalogue v1.1 is ever overwritten: every output goes under OUT = $SVI_REF_DIR
(default results/refine/) or the cut store CUT = $SVI_REF_CUT_DIR (default results/cut/svi_ext/).
Inputs read only: results/trackC/*, results/trackE/*, results/frozen_operating_point.json, the
20 Hz caches (lfe_extension, lfe_catalog_repair, lfe_picker) and the raw miniSEED.

"Scrambled control" = the same detector run on the same day with each station SITE's data shifted to a
different time of that day (one shift per site), so nothing lines up across sites except by chance.

Day roles (frozen in OUT/day_table.csv BEFORE any fit; see day_table()):
  prod          2014-01-01 .. 2026-09-15 production days (application only)
  fit           2011 truth days inside Bostock's catalogue window, and 2012 truth days (the row-level
                block role then splits 2012 into test blocks = fit, build blocks = insample, edges = guard)
  sealed        2013 truth days inside his window -- the one-time check (user decision 2026-10-02);
                scored ONCE by svi_ref_fit.py open_sealed, never fitted or tuned on
  s9            2010 truth days (9-station network): its own category, reported, never fitted
  nocov         truth-year days on disk outside his catalogue window (no labels; scrambled control only)
"""
from __future__ import annotations
import os, sys, json, glob, time, hashlib
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import c_common as C

ROOT = "/home/jovyan/tremorferometry"
EXT = C.EXT
RES = C.RES
REP = C.REP
PICK = f"{ROOT}/lfe_picker"
SOUTH = f"{ROOT}/lfe_south"
OUT = os.environ.get("SVI_REF_DIR", f"{RES}/refine")
CUT = os.environ.get("SVI_REF_CUT_DIR", f"{RES}/cut/svi_ext")
SEALED_YEAR = int(os.environ.get("SVI_REF_SEALED_YEAR", "2013"))
CACHE_R1 = os.environ.get("SVI_REF_CACHE_R1", f"{EXT}/cache20hz_r1")     # R1 rebuilds land here (never over cache20hz)
FS = C.FS
ERR_MAX = 0.03
NRAND = 1000
MATCH_S = 2.0
REC_THR = 8.0
PROD_T0, PROD_T1 = "2014-01-01", "2026-09-16"                   # [T0, T1): v1.1 runs 2014-01-01 .. 2026-09-15
MONTHS = [str(m) for m in pd.period_range("2014-01", "2026-09", freq="M")]
BOST = f"{REP}/results/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv"
FAMS12 = f"{REP}/results/e5r_families_2012.csv"
OFFS12 = f"{REP}/results/e5r_family_offsets_2012.csv"

S15 = list(C.STAS_ORDER)
S14 = [s for s in S15 if s != "LRIV"]
S13 = [s for s in S14 if s != "SNB"]
S9 = ["PGC", "VGZ", "NLLB", "SNB", "LRIV", "B009", "B010", "B011", "B926"]
SSETS = dict(S15=S15, S14=S14, S13=S13, S9=S9)
SITE_LIST = sorted(set(C.site_map(S15).values()))                 # fixed order -> a site's shift never depends on the set

TRUTH = {2010: dict(t0="2010-08-10", t1="2010-09-14", caches=[f"{REP}/cache20hz", f"{EXT}/cache20hz"], ssets=["S9"]),
         2011: dict(t0="2011-08-05", t1="2011-09-26", caches=[f"{PICK}/cache20hz"], ssets=["S15", "S14", "S13"]),
         2012: dict(t0="2012-08-15", t1="2012-10-06", caches=[f"{REP}/cache20hz"], ssets=["S15", "S14", "S13"]),
         2013: dict(t0="2013-08-30", t1="2013-10-29", caches=[f"{PICK}/cache20hz"], ssets=["S15", "S14", "S13"])}
BLOCK_S, GUARD_S = 7200.0, 60.0                                    # e5r_groundtruth.py: alternate 2-h blocks, 60-s guard

# station / instrument changes for the step checks (plan R7) and the CNSN swap audit (Merlin R2)
CHANGES = [("2014-06-17", "LZB ends"), ("2014-09-03", "YOUB ends"), ("2015-01-12", "PGC response sub-change"),
           ("2017-02-23", "NLLB BH->HH"), ("2017-03-03", "VGZ sensor swap"), ("2017-06-08", "SNB sensor swap"),
           ("2017-08-10", "PGC BH->HH"), ("2018-12-20", "PFB ends"), ("2019-05-30", "LRIV BH->HH, dropped"),
           ("2021-03-28", "SNB ends")]
SWAPS = [("PGC", "2015-01-12"), ("NLLB", "2017-02-23"), ("VGZ", "2017-03-03"), ("SNB", "2017-06-08"),
         ("PGC", "2017-08-10"), ("LRIV", "2019-05-30")]


def log(m):
    print(f"[{time.strftime('%H:%M:%S', time.gmtime())}] {m}", flush=True)


def ensure(*dirs):
    for d in dirs:
        os.makedirs(d, exist_ok=True)


def sha256(path, nbytes=None):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(1 << 22)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def to_epoch(ts):
    """datetime-like -> epoch seconds (cast to [ns] first: the pandas datetime64[us] trap)."""
    return pd.to_datetime(ts).values.astype("datetime64[ns]").astype(np.int64) / 1e9


# ------------------------------------------------------------------------------------------- truth
_GT = {}


def bostock(years=(2010, 2011, 2012, 2013), allow_sealed=False):
    """Bostock (2015) events mapped to bkey by EXACT (lat, lon, dep) match to e5r_families_2012.csv, time base =
    the FROZEN 2012 per-family offsets (t_origin = DB time + offset_s) -- lfe_picker o2_score.gt_year verbatim.
    The sealed year is withheld unless allow_sealed (only svi_ref_fit.open_sealed passes it)."""
    key = (tuple(years), allow_sealed)
    if key in _GT:
        return _GT[key].copy()
    # 2013 is lfe_picker's sealed referee: loadable ONLY by the one-time open stage of the real run
    # (SEALED_YEAR == 2013 and allow_sealed); a smoke run with a stand-in sealed year can never read it.
    years = [y for y in years if (allow_sealed or y != SEALED_YEAR) and (y != 2013 or (allow_sealed and SEALED_YEAR == 2013))]
    FAM = pd.read_csv(FAMS12); OFF = pd.read_csv(OFFS12)
    omap = {int(r.bkey): float(r.offset_s) for r in OFF.itertuples() if bool(r.accept)}
    B = pd.read_csv(BOST)
    B = B[B.year.isin(years)]
    B["time"] = pd.to_datetime(dict(year=B.year, month=B.month, day=B.day, hour=B.hour, minute=B["min"],
                                    second=B.sec.astype(int)), errors="coerce") + pd.to_timedelta(B.sec % 1.0, unit="s")
    k = {(round(r.lat, 4), round(r.lon, 4), round(r.dep, 4)): int(r.bkey) for r in FAM.itertuples()}
    B["fam"] = [k.get((round(a, 4), round(b, 4), round(c, 4)), -1) for a, b, c in zip(B.lat, B.lon, B.dep)]
    B = B[B.fam >= 0].copy()
    B["t_origin"] = to_epoch(B.time) + B.fam.map(omap)
    B = B.dropna(subset=["t_origin"]).sort_values("t_origin").reset_index(drop=True)
    B["date"] = pd.to_datetime(B.t_origin, unit="s").dt.strftime("%Y-%m-%d")
    B["year"] = B.date.str[:4].astype(int)
    B = B[["t_origin", "fam", "date", "year", "mag", "lat", "lon", "dep"]]
    _GT[key] = B
    return B.copy()


def coverage(year):
    """[first, last] UTC date of Bostock's events of that year inside the data window (his catalogue window)."""
    if year == 2013:
        return ("2013-09-13", "2013-10-13")         # lfe_picker/results/referee_protocol.json bostock_window: no 2013 row read
    B = bostock([year], allow_sealed=True)
    w = TRUTH[year]
    B = B[(B.date >= w["t0"]) & (B.date < w["t1"])]
    return (B.date.min(), B.date.max()) if len(B) else (None, None)


def block_role(t):
    """2012 row-level role from the origin-axis time: test block (odd) / build block (even) / guard edge."""
    t = np.asarray(t, float)
    blk = np.floor(t / BLOCK_S).astype(np.int64); pos = t - blk * BLOCK_S
    r = np.where(blk % 2 == 0, "insample", "fit").astype(object)
    r[(pos < GUARD_S) | (pos > BLOCK_S - GUARD_S)] = "guard"
    return r


def day_table():
    rows = [(d, "prod", "prod") for d in (str(x.date()) for x in pd.date_range(PROD_T0, PROD_T1, freq="D", inclusive="left"))]
    for y, w in TRUTH.items():
        a, b = coverage(y)
        for x in pd.date_range(w["t0"], w["t1"], freq="D", inclusive="left"):
            d = str(x.date())
            inside = a is not None and a <= d <= b
            if not inside: role = "nocov"
            elif y == SEALED_YEAR: role = "sealed"
            elif y == 2010: role = "s9"
            elif y == 2013: role = "unused2013"            # smoke runs only: 2013 stays sealed
            else: role = "fit"
            rows.append((d, f"truth{y}", role))
    D = pd.DataFrame(rows, columns=["date", "kind", "role"])
    D["year"] = D.date.str[:4].astype(int); D["month"] = D.date.str[:7]
    return D


def load_day_table():
    p = f"{OUT}/day_table.csv"
    if not os.path.exists(p):
        raise SystemExit(f"[common] {p} missing: run `svi_ref_common.py days` first (frozen before any fit)")
    return pd.read_csv(p, dtype={"date": str})


# ------------------------------------------------------------------------------------------- scrambled control
def site_shifts(date, nv):
    """deterministic per-site shift (samples) for the scrambled control of `date`: one per site of SITE_LIST (fixed
    order, so a site's shift does not depend on the station set), in [10 min, nv - 10 min], >= 2 min apart."""
    rng = np.random.default_rng(20261002 + int(date.replace("-", "")))
    lo, hi = int(600 * FS), int(nv - 600 * FS)
    while True:
        v = rng.integers(lo, hi, len(SITE_LIST))
        if np.min(np.diff(np.sort(v))) >= 120 * FS:
            break
    return {q: int(s) for q, s in zip(SITE_LIST, v)}


def rand_idx(date, fam, valid_idx, n=NRAND):
    rng = np.random.default_rng(3000 + int(date.replace("-", "")) * 131 + int(fam))
    if len(valid_idx) == 0:
        return np.zeros(0, np.int64)
    return np.sort(rng.choice(valid_idx, min(n, len(valid_idx)), replace=len(valid_idx) < n))


# ------------------------------------------------------------------------------------------- status + guard
def write_status(stage, tag, rows):
    ensure(f"{OUT}/status")
    pd.DataFrame(rows, columns=["date", "netsta", "status", "msg"]).to_csv(f"{OUT}/status/{stage}_{tag}.csv", index=False)


def guard(stage, pattern=None, out=None):
    """exit 5 when > ERR_MAX of the station-days the stage was expected to read errored or are silently missing
    (status 'error' or 'missing'); legitimate absences (no data at the data centre / no epoch) are 'nodata'."""
    out = out or OUT
    pat = pattern or f"{stage}_*"
    fs = sorted(glob.glob(f"{out}/status/{pat}.csv"))
    if not fs:
        log(f"[guard:{stage}] no status files ({out}/status/{pat}.csv)"); sys.exit(5)
    S = pd.concat([pd.read_csv(f, dtype=str) for f in fs], ignore_index=True).fillna("")
    S = S.drop_duplicates(["date", "netsta"], keep="last")
    bad = S.status.isin(["error", "missing"])
    n_ok = int((S.status == "ok").sum()); n_bad = int(bad.sum())
    frac = n_bad / max(n_ok + n_bad, 1)
    msg = f"[guard:{stage}] station-days ok {n_ok}, error/missing {n_bad} ({frac:.2%}); all {S.status.value_counts().to_dict()}"
    log(msg)
    ensure(f"{out}/status")
    open(f"{out}/status/GUARD_{stage}.txt", "w").write(msg + "\n")
    if frac > ERR_MAX:
        log(f"[guard:{stage}] FAIL: {frac:.2%} > {ERR_MAX:.0%}; examples:\n" + S[bad].head(20).to_string())
        sys.exit(5)
    return frac


# ------------------------------------------------------------------------------------------- matching
def nearest_dt(tq, fq, tref, fref, same=True):
    """|dt| to the nearest reference event (same family if same)."""
    tq = np.asarray(tq, float); out = np.full(len(tq), np.inf)
    if not len(tref):
        return out
    if not same:
        r = np.sort(np.asarray(tref, float))
        j = np.clip(np.searchsorted(r, tq), 1, max(len(r) - 1, 1))
        return np.minimum(np.abs(r[j - 1] - tq), np.abs(r[np.minimum(j, len(r) - 1)] - tq))
    fq = np.asarray(fq); fref = np.asarray(fref); tref = np.asarray(tref, float)
    for f in np.unique(fq):
        m = fq == f; r = np.sort(tref[fref == f])
        if not len(r):
            continue
        j = np.clip(np.searchsorted(r, tq[m]), 1, max(len(r) - 1, 1))
        out[m] = np.minimum(np.abs(r[j - 1] - tq[m]), np.abs(r[np.minimum(j, len(r) - 1)] - tq[m]))
    return out


def one2one_recall(t_det, f_det, t_ev, f_ev, tol=MATCH_S):
    """lfe_picker o2_score.one2one per family (own family, +-2 s, greedy by time): -> matched mask over the EVENTS."""
    t_det = np.asarray(t_det, float); f_det = np.asarray(f_det); t_ev = np.asarray(t_ev, float); f_ev = np.asarray(f_ev)
    hit = np.zeros(len(t_ev), bool)
    for f in np.unique(f_ev):
        ie = np.flatnonzero(f_ev == f); b = np.sort(t_det[f_det == f])
        if not len(b):
            continue
        o = ie[np.argsort(t_ev[ie], kind="stable")]; j = 0; nb = len(b)
        for i in o:
            ti = t_ev[i]
            while j < nb and b[j] < ti - tol:
                j += 1
            if j < nb and abs(b[j] - ti) <= tol:
                hit[i] = True; j += 1
    return hit


def netsta_of():
    ep = C.epochs()
    return {s: ep[ep.sta == s].netsta.iloc[0] for s in sorted(ep.sta.unique())}


def daylists():
    """OUT/days/*.json for the queue (smoke runs restrict them through SVI_REF_* environment lists)."""
    od = f"{OUT}/days"; ensure(od)
    env = lambda k: [x for x in os.environ.get(k, "").split(",") if x]
    DT = load_day_table()
    prod = env("SVI_REF_PROD_DAYS") or DT[DT.kind == "prod"].date.tolist()
    json.dump(prod, open(f"{od}/prod.json", "w"))
    first = []
    for y in TRUTH:
        d = env(f"SVI_REF_TRUTH_DAYS_{y}") or DT[DT.kind == f"truth{y}"].date.tolist()
        json.dump(d, open(f"{od}/truth_{y}.json", "w")); first += d
    fm = json.load(open(f"{OUT}/r3/flagged_months.json")) if os.path.exists(f"{OUT}/r3/flagged_months.json") else []
    r3m = env("SVI_REF_R3M_DAYS") or [d for d in DT[DT.kind == "prod"].date if d[:7] in set(fm)]
    json.dump(r3m, open(f"{od}/r3m.json", "w")); first += r3m
    s13 = env("SVI_REF_S13_DAYS") or [d for d in prod if d < "2021-04-01"]
    json.dump(s13, open(f"{od}/prodS13.json", "w"))
    json.dump(sorted(set(first)), open(f"{od}/prep_first.json", "w"))
    json.dump(sorted(set(prod) | set(s13) - set(first)), open(f"{od}/prep_rest.json", "w"))
    log(f"[days] prod {len(prod)}, truth {[(y, len(json.load(open(f'{od}/truth_{y}.json')))) for y in TRUTH]}, r3m {len(r3m)}, "
        f"S13 {len(s13)}, prep first {len(set(first))}")


def busiest_month():
    best = (None, -1)
    for m in MONTHS:
        p = f"{C.TRC}/run_{m}.json"
        if os.path.exists(p):
            n = json.load(open(p)).get("n_detections", 0)
            if n > best[1]: best = (m, n)
    return best[0]


if __name__ == "__main__":
    if sys.argv[1] == "days":
        ensure(OUT)
        p = f"{OUT}/day_table.csv"
        if os.path.exists(p) and "--force" not in sys.argv:
            log(f"[days] {p} exists (frozen); not rewritten"); sys.exit(0)
        D = day_table(); D.to_csv(p, index=False)
        json.dump(dict(created_utc=pd.Timestamp.now("UTC").isoformat(), sha256=sha256(p), sealed_year=SEALED_YEAR,
                       coverage={y: coverage(y) for y in TRUTH}, roles=D.groupby(["kind", "role"]).size().astype(int)
                       .reset_index().astype(str).values.tolist()), open(f"{OUT}/day_table_meta.json", "w"), indent=1)
        log(f"[days] {D.groupby(['kind', 'role']).size().to_dict()}; coverage {[(y, coverage(y)) for y in TRUTH]}")
    elif sys.argv[1] == "daylists":
        daylists()
    elif sys.argv[1] == "busiest":
        print(os.environ.get("SVI_REF_PILOT_MONTH") or busiest_month())
    elif sys.argv[1] == "guard":
        guard(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
