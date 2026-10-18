#!/usr/bin/env python
"""Shared definitions for the Washington PRODUCTION run (results/WA_PROD_PLAN.md W0-W13 as amended by Merlin's review and
the user's decisions of 2026-10-03).  Every output goes under OUT = $WA_PROD_DIR (default results/wa_prod/), the cut store
under CUT = $WA_PROD_CUT_DIR (default results/cut/wa_ext/), new 20 Hz caches under CACHE_NEW = $WA_PROD_CACHE (default
lfe_extension_wa/cache20hz_prod/).  Nothing that exists is overwritten; nothing is deleted; no git.

"Scrambled control" (series c1, c2) = the same detector with each station SITE's cc series shifted to a different time of
the same day (one shift per site), so nothing lines up across sites except by chance.  "Time-reversed control" (series
rev) = the same detector with every template channel reversed in time (same offsets, weights, channels, guard): the main
false-alarm reference (user 2026-10-03).

Networks (one matched-filter pass computes several network sums from the same station cc series):
  core   the FIXED rate network of the era (E1 2007-2011, E2 2012-2014, E3 2015-2026-09-15 = the headline set);
         a station is in an era set by the rule in wa_prod_setup.core_sets (written before any run)
  E1/E2/E3  the three era sets, all computed on the truth days (k factors and the confidence fit for every set)
  all    every station with a W4a template (frozen T1 / V2010-build rebuild / carry-over approved by W2) that day
  rich   'all' + the W4b templates (re-templated control-level epochs + T2B new nodes); rich pass only
"""
from __future__ import annotations
import os, sys, json, glob, time, hashlib, zlib
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wa_common as W

ROOT = W.ROOT
STUDY = W.STUDY
RES = W.RES
OUT = os.environ.get("WA_PROD_DIR", f"{RES}/wa_prod")
CUT = os.environ.get("WA_PROD_CUT_DIR", f"{RES}/cut/wa_ext")
CACHE_NEW = os.environ.get("WA_PROD_CACHE", f"{STUDY}/cache20hz_prod")
CACHE_READ = [CACHE_NEW, W.CACHE, f"{RES}/wa_pre/w2/cache20hz"]          # first hit wins (all built by wa_cache.cache_one)
SMOKE = os.environ.get("WA_PROD_SMOKE", "0") == "1"
LOGF = os.environ.get("WA_PROD_LOG", f"{STUDY}/logs/wa_prod.log")   # override: fixed-10 re-run (2026-10-07)

T0, T1 = "2007-01-01", "2026-09-15"                                       # production span, inclusive
ERAS = {"E1": ("2007-01-01", "2011-12-31"), "E2": ("2012-01-01", "2014-12-31"), "E3": ("2015-01-01", "2026-09-15")}
SEALED_MONTHS = {"2007-01": "sealed_headline", "2009-05": "sealed_semi"}
EXCLUDED_MONTHS = {"2011-04"}
FIT_MONTHS = {"2007-10", "2008-05", "2009-02", "2009-04", "2010-03", "2010-07", "2010-08", "2010-09", "2011-08", "2011-12"}
TRUTH_MONTHS = sorted(set(SEALED_MONTHS) | EXCLUDED_MONTHS | FIT_MONTHS)
WINDOW_BUILD = {"V2008": ("2008-05-01", "2008-06-18"), "V2010": ("2010-07-14", "2010-09-18")}   # build halves in-sample

# raw roots: wa_common.TRACE_ROOTS (incl. traces_dl) + traces_dl2 once WA_DL2_DONE exists + the roots wa_cache did not read
DL2 = f"{STUDY}/traces_dl2"
DL2_DONE = f"{RES}/wa_dl2/WA_DL2_DONE"
DL2_STOP = f"{RES}/wa_dl2/STOPPED_STORAGE"
EXTRA_ROOTS = [f"{ROOT}/lfe_picker/traces", f"{ROOT}/lfe_south/traces", f"{ROOT}/data/waveforms"]
NEVER_ROOT = ("traces_sealed",)
BANNED_NETS = {"7D", "X6", "YR", "Z5"}                                    # user 2026-10-03 (+ XD nodes, below)
NO_TEMPLATE = {"XU.W060"}                                                 # bursty: never in any template pool
REGIONS = {"Royer_NW_WA": "NW", "Royer_SW_WA": "SW", "Olympic_2010": "OlyA", "Olympic_2010_armB": "ArmB"}
TREM_REG = {"Royer_NW_WA": "NW_WA", "Royer_SW_WA": "SW_WA", "Olympic_2010": "OLY2010", "Olympic_2010_armB": "OLY2010"}
TPL_TAGS_W4A = ["build", "v10b", "t2012new", "step4b_t2012", "step4b_t2013", "step4b_t2014"]   # first listed wins
INSAMPLE_TAGS = {"t2012new": "T2012", "step4b_t2012": "T2012", "step4b_t2013": "T2013", "step4b_t2014": "T2014"}
P_CUT = 0.75
CUTS = [0.5, 0.6, 0.7, 0.75, 0.8, 0.9]
ERR_MAX = 0.03
NRAND = 1000


def dl2_ready():
    return os.path.exists(DL2_DONE) or os.path.exists(DL2_STOP)


def prod_roots():
    r = [x for x in W.TRACE_ROOTS if not any(n in x for n in NEVER_ROOT)]
    if dl2_ready() and DL2 not in r:
        r.append(DL2)
    for x in EXTRA_ROOTS:
        if x not in r:
            r.append(x)
    return r


def patch_wa():
    """make the frozen wa_* modules read the production roots and the production cache list (never writes elsewhere)."""
    W.TRACE_ROOTS[:] = prod_roots()

    def cache_path(netsta, date, _dirs=tuple(CACHE_READ)):
        ts = pd.Timestamp(date)
        tail = f"{netsta}/{ts.year}/{int(ts.dayofyear):03d}.npz"
        for c in _dirs:
            p = f"{c}/{tail}"
            if os.path.exists(p):
                return p
        return f"{_dirs[0]}/{tail}"
    W.cache_path = cache_path
    W.CACHE = CACHE_NEW
    return W


def new_cache_path(netsta, date):
    ts = pd.Timestamp(date)
    return f"{CACHE_NEW}/{netsta}/{ts.year}/{int(ts.dayofyear):03d}.npz"


def find_cache(netsta, date):
    ts = pd.Timestamp(date)
    tail = f"{netsta}/{ts.year}/{int(ts.dayofyear):03d}.npz"
    for c in CACHE_READ:
        p = f"{c}/{tail}"
        if os.path.exists(p):
            return p
    return None


def raw_path(netsta, date, roots=None):
    ts = pd.Timestamp(date)
    for r in (roots or prod_roots()):
        p = f"{r}/{netsta}/{ts.year}/{int(ts.dayofyear):03d}.mseed"
        if os.path.exists(p) and os.path.getsize(p) > 0:
            if any(n in os.path.realpath(p) for n in NEVER_ROOT):
                continue
            return p
    return None


def log(*a):
    os.makedirs(os.path.dirname(LOGF), exist_ok=True)
    W.log(LOGF if not SMOKE else None, *a)


def ensure(*d):
    for x in d:
        os.makedirs(x, exist_ok=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(1 << 22)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def days(t0=T0, t1=T1):
    return [str(d.date()) for d in pd.date_range(t0, t1, freq="D")]


def era_of(date):
    for k, (a, b) in ERAS.items():
        if a <= date <= b:
            return k
    return None


def month_role(m):
    if m in SEALED_MONTHS:
        return SEALED_MONTHS[m]
    if m in EXCLUDED_MONTHS:
        return "excluded"
    if m in FIT_MONTHS:
        return "fit"
    return "prod"


def site_shift(date, site, nv, k=1):
    """scrambled-control shift (samples) of a SITE on a day: deterministic in (date, site, control k), independent of
    which other sites are present (SVI_REF fix), uniform in [300, 3000] s (the frozen WA ROLL_S range)."""
    s = zlib.crc32(f"{date}|{site}|{k}".encode())
    rng = np.random.default_rng(s)
    return int(rng.uniform(*W.ROLL_S) * W.FS) % nv


# ------------------------------------------------------------------------------------------------ status + guard
def write_status(stage, tag, rows):
    ensure(f"{OUT}/status")
    pd.DataFrame(rows, columns=["date", "netsta", "status", "msg"]).to_csv(f"{OUT}/status/{stage}_{tag}.csv", index=False)


def guard(stage, pattern=None, fail=True):
    """FAIL (exit 5) when > 3 % of the station-days a stage was expected to read errored or are silently missing.
    Status values: ok / error / missing (counted) ; nodata / excluded / notemplate / dead (not counted: legitimate)."""
    pat = pattern or f"{stage}_*"
    fs = sorted(glob.glob(f"{OUT}/status/{pat}.csv"))
    if not fs:
        log(f"[guard:{stage}] no status files {OUT}/status/{pat}.csv")
        if fail:
            sys.exit(5)
        return 1.0
    S = pd.concat([pd.read_csv(f, dtype=str) for f in fs], ignore_index=True).fillna("")
    S = S.drop_duplicates(["date", "netsta"], keep="last")
    bad = S.status.isin(["error", "missing"])
    n_ok = int((S.status == "ok").sum()); n_bad = int(bad.sum())
    frac = n_bad / max(n_ok + n_bad, 1)
    msg = f"[guard:{stage}] station-days ok {n_ok}, error/missing {n_bad} ({frac:.2%}); all {S.status.value_counts().to_dict()}"
    log(msg)
    ensure(f"{OUT}/status")
    open(f"{OUT}/status/GUARD_{stage}.txt", "w").write(msg + "\n")
    if frac > ERR_MAX and fail:
        log(f"[guard:{stage}] FAIL {frac:.2%} > {ERR_MAX:.0%}:\n" + S[bad].head(30).to_string())
        sys.exit(5)
    return frac


# ------------------------------------------------------------------------------------------------ tables from setup
def families():
    return pd.read_csv(f"{OUT}/setup/families.csv")


def stations():
    return pd.read_csv(f"{OUT}/setup/stations.csv", keep_default_na=False)


def nbr_map():
    """{fid: [fids <= 3 km, itself included]}; a family without a position is its own only neighbour."""
    F = families()
    out = {}
    for r in F.itertuples():
        v = [int(float(x)) for x in str(r.nbr3).split(";") if x not in ("", "nan")]
        out[int(r.fid)] = v or [int(r.fid)]
    return out


def core_sets():
    return json.load(open(f"{OUT}/setup/core_sets.json"))


def day_table():
    return pd.read_csv(f"{OUT}/setup/day_table.csv", dtype=str)


def tremor_days():
    out = {}
    for rg in ["NW_WA", "SW_WA", "OLY2010"]:
        t = pd.read_csv(f"{W.EXT}/{rg}_tremor_days.csv")
        out[rg] = dict(zip(pd.to_datetime(t.day).dt.strftime("%Y-%m-%d"), t.n_windows.astype(int)))
    return out


def to_epoch(ts):
    return pd.to_datetime(ts).values.astype("datetime64[ns]").astype(np.int64) / 1e9


def marker(name, text=""):
    ensure(OUT)
    open(f"{OUT}/{name}", "w").write(time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + " " + text + "\n")


def done(name):
    return os.path.exists(f"{OUT}/{name}")
