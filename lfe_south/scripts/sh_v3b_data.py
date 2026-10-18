#!/usr/bin/env python
"""SH_V3B (2026-09-29): data + prep for the 4 S1 stations v3 lacked (BK.GHOP, NC.KMPB, PB.B932, BK.HUNT).

  plan : results/sh_v3b/download_plan.csv = (net, sta, day) for every Shelly2025 build + test day
         (results/aa_shelly/stage2_days_{build,test}.json, 2018-01-01..2024-12-01) inside each station's
         FDSN operating period (results/fdsn_channels_south.csv, HH/EH/BH channels).  Download and cache
         then run with the UNCHANGED t2_fetch.py / t3_cache.py (native-rate miniSEED -> 1-8 Hz, 20 sps
         cache20hz/<NET.STA>/<yyyy>/<doy>.npz, HH > BH > EH > SH band preference, no response removal --
         exactly how the existing cache was built).
  prep : results/sh_v3b/prep/<date>.npz.  prep_day (sh_v2_common, unchanged) is run on the v2 prep
         station set (aa_shelly/templates6_build, 23 stations) + the 4 new ones into results/sh_v3b/prep_full/;
         the merged file keeps the v2 prep entries (mute / lines / ref) of the stations v2 already prepared
         BYTE-IDENTICAL (so v3 -> v3b changes only the station set) and adds the new stations' entries
         from the 27-station run (their burst ratios see all neighbours, as v2 would have).
  python scripts/sh_v3b_data.py plan | prep [--workers 12] [--limit N]
"""
import os, sys, json, argparse, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sh_v2_common as C2
from t_common import RES, log

NEW = ["BK.GHOP", "NC.KMPB", "PB.B932", "BK.HUNT"]
V3B = f"{RES}/sh_v3b"; PD = os.environ.get("SH_V3B_PREPROOT", V3B)
FULL = f"{PD}/prep_full"; OUT = f"{PD}/prep"; V2P = f"{RES}/sh_v2/prep"
DAYFILES = [f"{RES}/aa_shelly/stage2_days_build.json", f"{RES}/aa_shelly/stage2_days_test.json"]
T0, T1 = pd.Timestamp("2018-01-01"), pd.Timestamp("2024-12-01")


def all_days():
    days = set()
    for p in DAYFILES:
        for k, v in json.load(open(p)).items():
            if k.split("|")[0] == "Shelly2025": days |= {pd.Timestamp(x) for x in v}
    return sorted(d for d in days if T0 <= d < T1)


def plan():
    F = pd.read_csv(f"{RES}/fdsn_channels_south.csv", dtype=str)
    F = F[F.cha.str[:2].isin(["HH", "BH", "EH", "SH"])]
    rows = []
    for ns in NEW:
        net, sta = ns.split(".")
        g = F[(F.net == net) & (F.sta == sta)]
        per = [(pd.Timestamp(a), pd.Timestamp(b[:10]) if not b.startswith("3000") else pd.Timestamp("2100-01-01"))
               for a, b in zip(g.start, g.end)]
        for d in all_days():
            if any(a <= d < b for a, b in per): rows.append(dict(net=net, sta=sta, day=str(d.date())))
    P = pd.DataFrame(rows); os.makedirs(V3B, exist_ok=True); P.to_csv(f"{V3B}/download_plan.csv", index=False)
    log(f"[v3b plan] {len(P)} station-days: {P.groupby('sta').size().to_dict()}")


_G = {}
def _job(d):
    warnings.filterwarnings("ignore")
    d = pd.Timestamp(d); out = f"{OUT}/{d.date()}.npz"
    if os.path.exists(out): return []
    rows = C2.prep_day(d, _G["stas"], _G["C"], _G["nref"])          # -> FULL/<date>.npz (27 stations)
    zf = dict(np.load(f"{FULL}/{d.date()}.npz"))
    p2 = f"{V2P}/{d.date()}.npz"
    z2 = dict(np.load(p2)) if os.path.exists(p2) else {"stas": np.array([], dtype=str), "nref": np.int64(_G["nref"])}
    s2 = [str(s) for s in z2["stas"]]
    add = [str(s) for s in zf["stas"] if str(s) in NEW]
    save = {k: v for k, v in z2.items() if k not in ("stas", "nref")}
    for ns in add:
        for q in ("mute", "lines", "ref"): save[f"{ns}|{q}"] = zf[f"{ns}|{q}"]
    tmp = out + ".tmp.npz"
    np.savez_compressed(tmp, stas=np.array(s2 + add), nref=z2.get("nref", _G["nref"]), **save)
    os.replace(tmp, out)
    return [r for r in rows if r["netsta"] in NEW]


def prep(workers, limit, only=""):
    from concurrent.futures import ProcessPoolExecutor
    os.makedirs(FULL, exist_ok=True); os.makedirs(OUT, exist_ok=True)
    C2.PREP = FULL                                                   # prep_day writes here (read by workers via fork)
    M = pd.read_csv(f"{RES}/aa_shelly/templates6_build.csv"); M = M[M.source == "Shelly2025"]
    _G["stas"] = sorted(set(M.netsta.unique()) | set(NEW)); _G["C"] = C2.coords(); _G["nref"] = int(round(6.0 * C2.FS))
    miss = [s for s in NEW if s not in _G["C"]]
    if miss: raise SystemExit(f"no coordinates for {miss}")
    days = [str(d.date()) for d in all_days()]
    if limit: days = days[:limit]
    if only: days = [d for d in days if d in only.split(",")]
    log(f"[v3b prep] {len(days)} days, {len(_G['stas'])} stations (new {NEW}), workers={workers}")
    rows = []
    with ProcessPoolExecutor(workers) as ex:
        for i, r in enumerate(ex.map(_job, days)):
            rows += r
            if i % 50 == 0: log(f"  {i + 1}/{len(days)}")
    old = f"{PD}/prep_summary_new.csv"
    df = pd.DataFrame(rows)
    if os.path.exists(old): df = pd.concat([pd.read_csv(old), df]).drop_duplicates(["day", "netsta"], keep="last")
    df.to_csv(old, index=False)
    if len(df):
        log(f"[v3b prep] new-station days {df.groupby('netsta').size().to_dict()}; muted frac mean "
            f"{df.groupby('netsta').muted_frac.mean().round(4).to_dict()}; station-dead {df.groupby('netsta').station_dead.sum().to_dict()}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what", choices=["plan", "prep"])
    ap.add_argument("--workers", type=int, default=12); ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--only", default="", help="comma list of dates (smoke)")
    a = ap.parse_args()
    plan() if a.what == "plan" else prep(a.workers, a.limit, a.only)
