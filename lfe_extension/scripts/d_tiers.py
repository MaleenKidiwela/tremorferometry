#!/usr/bin/env python
"""Track D — the epoch axis and the per-epoch tier thresholds.

**Track C is the authority on both.**  `results/trackC/tiers_by_epoch.csv` defines the NETWORK EPOCHS
and, for each, re-derives T_core and T_ext by the frozen false-alarm-only rule (T_core: quiet-day FA
<= 0.1 per template-day; T_ext: quiet-day <= 1.0 AND post-SSE <= 1.5), applied as a re-flag of the
`cat_core` / `cat_ext` columns.  Every Track D step reads the epochs and the thresholds through this one
loader, so the completeness model, the floor lambda and the catalog can never disagree with the flags.
A row whose threshold is blank (not enough quiet family-days yet) falls back to the frozen referee point.

Track C's epochs are instrument/roster epochs; Track D's per-family recall is measured by ablating the
STATION SET on the referee years, and two of Track C's epochs share a station set.  The map is:

| Track C epoch | months | station set | referee ablation run |
|---|---|---|---|
| E0 | 2014-01..2017-02 | 15 stations | `a1_mergeB_L12` (the frozen 2012 run) |
| E1 | 2017-03..2018-12 | 15 stations (same set; CN sensors upgraded, 2012 templates kept per A-2 step 4) | `a1_mergeB_L12` |
| E2 | 2019-01..2021-03 | 14 from 2019-05-30 (LRIV dropped); the first five months still have 15 | `d1_abl2012_E2` |
| E3 | 2021-04.. | 13 (LRIV dropped + SNB ends) | `d1_abl2012_E3` |

If the file is absent the fallback axis is the station-set axis itself (S15 / S14 / S13)."""
from __future__ import annotations
import os
import numpy as np, pandas as pd

EXT = "/home/jovyan/tremorferometry/lfe_extension"
TIERS_CSV = f"{EXT}/results/trackC/tiers_by_epoch.csv"
FROZEN = dict(T_ext=8.25, T_core=8.75)

# station sets (the thing the ablation measures) and the referee runs that realise them
STATION_SETS = {
    "S15": dict(n=15, stations="PGC,VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926",
                run2012="a1_mergeB_L12", run2010="a2_xfer2010",
                note="the frozen 15-station network"),
    "S14": dict(n=14, stations="PGC,VGZ,NLLB,SNB,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926",
                run2012="d1_abl2012_E2", run2010="d1_abl2010_E2",
                note="LRIV dropped from 2019-05-30 (A-2 epoch-identity FAIL)"),
    "S13": dict(n=13, stations="PGC,VGZ,NLLB,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926",
                run2012="d1_abl2012_E3", run2010="d1_abl2010_E3",
                note="LRIV dropped + CN.SNB ends 2021-03-28"),
}
# Track C epoch -> station set.  Anything unknown falls back to S15 unless its months start after the
# SNB end, in which case S13 (so a re-cut epoch table cannot silently mislabel the network).
EPOCH_SET = {"E0": "S15", "E1": "S15", "E2": "S14", "E3": "S13"}
FALLBACK_EPOCHS = [("S15", "2014-01-01", "2019-05-29"), ("S14", "2019-05-30", "2021-03-27"),
                   ("S13", "2021-03-28", "2026-12-31")]


def _pick(cols, *names):
    low = {str(c).lower(): c for c in cols}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    return None


def _months_to_dates(m):
    """'2014-01..2017-02' -> ('2014-01-01', '2017-02-28')."""
    try:
        a, b = str(m).split("..")
        t0 = pd.Period(a.strip(), freq="M").start_time
        t1 = pd.Period(b.strip(), freq="M").end_time
        return str(t0.date()), str(t1.date())
    except Exception:
        return None


def load(path=TIERS_CSV):
    """-> (EPOCHS, TIERS, source)

    EPOCHS  list of dicts: epoch_id, t0, t1, station_set, n_stations, note
    TIERS   epoch_id -> {'T_core': x, 'T_ext': y, 'derived': bool}
    """
    if not os.path.exists(path):
        EP = [dict(epoch_id=e, t0=a, t1=b, station_set=e, n_stations=STATION_SETS[e]["n"],
                   note=STATION_SETS[e]["note"]) for e, a, b in FALLBACK_EPOCHS]
        TI = {e["epoch_id"]: dict(FROZEN, derived=False) for e in EP}
        return EP, TI, (f"{path} NOT PRESENT -> the station-set axis S15/S14/S13 with the frozen "
                        f"referee point 8.25 / 8.75")
    R = pd.read_csv(path)
    ce = _pick(R.columns, "epoch", "epoch_id", "name")
    cm = _pick(R.columns, "months", "month_range", "range")
    cc = _pick(R.columns, "T_core", "t_core", "tcore")
    cx = _pick(R.columns, "T_ext", "t_ext", "text")
    cn = _pick(R.columns, "note")
    cs = _pick(R.columns, "status")
    if ce is None:
        EP = [dict(epoch_id=e, t0=a, t1=b, station_set=e, n_stations=STATION_SETS[e]["n"],
                   note=STATION_SETS[e]["note"]) for e, a, b in FALLBACK_EPOCHS]
        return EP, {e["epoch_id"]: dict(FROZEN, derived=False) for e in EP}, \
            f"UNREADABLE {path} (columns {list(R.columns)}) -> station-set axis, frozen thresholds"
    EP, TI = [], {}
    for r in R.itertuples():
        e = str(getattr(r, ce))
        rng = _months_to_dates(getattr(r, cm)) if cm else None
        t0, t1 = rng if rng else ("1900-01-01", "1900-01-01")
        ss = EPOCH_SET.get(e)
        if ss is None:                       # unknown epoch label: infer from its start date
            ss = "S13" if t0 >= "2021-03-28" else ("S14" if t0 >= "2019-05-30" else "S15")
        EP.append(dict(epoch_id=e, t0=t0, t1=t1, station_set=ss, n_stations=STATION_SETS[ss]["n"],
                       note=(str(getattr(r, cn)) if cn else "") + f" | station set {ss}: "
                            + STATION_SETS[ss]["note"]))
        vc = getattr(r, cc) if cc else np.nan
        vx = getattr(r, cx) if cx else np.nan
        ok = np.isfinite(pd.to_numeric(vc, errors="coerce")) and np.isfinite(pd.to_numeric(vx, errors="coerce"))
        st = str(getattr(r, cs)) if cs else ""
        TI[e] = dict(T_core=float(vc) if ok else FROZEN["T_core"],
                     T_ext=float(vx) if ok else FROZEN["T_ext"],
                     derived=bool(ok and st.startswith("derived")), status=st)
    EP.sort(key=lambda d: d["t0"])
    return EP, TI, f"{path} (Track C, FA-only rule; {len(EP)} epochs)"


def epoch_of(dates, EPOCHS=None):
    EPOCHS = EPOCHS if EPOCHS is not None else load()[0]
    d = np.asarray(dates).astype(str)
    out = np.full(len(d), "", dtype=object)
    for e in EPOCHS:
        out[(d >= e["t0"]) & (d <= e["t1"])] = e["epoch_id"]
    return out


def set_of(epoch_id, EPOCHS=None):
    EPOCHS = EPOCHS if EPOCHS is not None else load()[0]
    for e in EPOCHS:
        if e["epoch_id"] == epoch_id:
            return e["station_set"]
    return "S15"


def thresholds_for_set(ss, EPOCHS=None, TIERS=None):
    """the thresholds to score a referee ablation of station set `ss` with: the tightest (largest T_core)
    of the epochs that use that set, so the referee number is never more optimistic than production."""
    EPOCHS = EPOCHS if EPOCHS is not None else load()[0]
    TIERS = TIERS if TIERS is not None else load()[1]
    cand = [TIERS[e["epoch_id"]] for e in EPOCHS if e["station_set"] == ss and e["epoch_id"] in TIERS]
    if not cand:
        return dict(FROZEN)
    return dict(T_core=max(c["T_core"] for c in cand), T_ext=max(c["T_ext"] for c in cand))


if __name__ == "__main__":
    EP, TI, src = load()
    print(f"[tiers] source: {src}")
    for e in EP:
        t = TI[e["epoch_id"]]
        print(f"  {e['epoch_id']}: {e['t0']} .. {e['t1']}  set {e['station_set']} "
              f"({e['n_stations']} stations)  T_core {t['T_core']}  T_ext {t['T_ext']}"
              f"  {'derived' if t.get('derived') else 'frozen fallback'}  [{t.get('status','')}]")
    for ss in ("S15", "S14", "S13"):
        print(f"  station set {ss}: scored at {thresholds_for_set(ss, EP, TI)}")
