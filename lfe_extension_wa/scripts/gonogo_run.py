#!/usr/bin/env python
"""GO/NO-GO Step B runner (Merlin 2026-10-02): the T1 all-templates matched filter with the network RESTRICTED to
the 21 continuous on-disk core stations (results/wa_prod_plan/inventory_by_station.csv, raw 2014-26 >= 3,900 days).

Nothing in the frozen chain is changed: this wrapper only points wa_common.RES at results/wa_gonogo (templates are
read through the symlink results/wa_gonogo/templates -> ../templates, i.e. the SAME build-half T1 templates) and
WA_SETUP at results/wa_gonogo/setup (T1 setup with stations.csv cut to the core).  Outputs never touch existing
results.

  python scripts/gonogo_run.py setup
  python scripts/gonogo_run.py guard            # station-day completeness of the core on V2010 + A2011 (pre-run)
  python scripts/gonogo_run.py mf  -- --tag core21 --windows V2010,A2011 --tpl build --n-control 2 --drop-dead 0.05 --all-templates
  python scripts/gonogo_run.py hf  -- --mf core21 --tpl build --workers 8 --all-templates
  python scripts/gonogo_run.py postguard core21 # every on-disk core station-day loaded (or dropped as > 5 % guarded)
"""
import os, sys, json, glob, shutil
HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
GO = f"{STUDY}/results/wa_gonogo"
os.environ["WA_SETUP"] = f"{GO}/setup"
sys.path.insert(0, HERE)
import numpy as np
import pandas as pd
import wa_common as W

W.RES = GO
W.SETUP = f"{GO}/setup"
# lfe_array_discovery/traces_sealed is EXCLUDED: it is that project's sealed holdout (never read, never counted)
RAW_ROOTS = ["lfe_array_discovery/traces", "lfe_extension/traces",
             "lfe_catalog_repair/traces", "lfe_extension_wa/traces", "lfe_picker/traces", "lfe_south/traces",
             "data/waveforms"]


def core_list():
    inv = pd.read_csv(f"{STUDY}/results/wa_prod_plan/inventory_by_station.csv")
    c = inv[inv.raw_c14_26 >= 3900].netsta.tolist()
    assert len(c) == 21, c
    return c


def raw_any(ns, d):
    ts = pd.Timestamp(d)
    for r in RAW_ROOTS:
        p = f"{W.ROOT}/{r}/{ns}/{ts.year}/{int(ts.dayofyear):03d}.mseed"
        if os.path.exists(p) and os.path.getsize(os.path.realpath(p)) > 0:
            return True
    return False


def setup():
    os.makedirs(W.SETUP, exist_ok=True)
    src = f"{STUDY}/results/merlin_tests/t1_setup"
    for f in ["families.csv", "truth.parquet", "scan_days.csv"]:
        shutil.copy(f"{src}/{f}", f"{W.SETUP}/{f}")
    S = pd.read_csv(f"{src}/stations.csv")
    core = core_list()
    missing = sorted(set(core) - set(S.netsta))
    S[S.netsta.isin(core)].to_csv(f"{W.SETUP}/stations.csv", index=False)
    if not os.path.exists(f"{GO}/templates"):
        os.symlink("../templates", f"{GO}/templates")
    json.dump(dict(core=core, not_in_trial_roster_no_build_template=missing,
                   n_core_used=int(S.netsta.isin(core).sum()), n_sites=int(S[S.netsta.isin(core)].site.nunique())),
              open(f"{W.SETUP}/core.json", "w"), indent=1)
    print(open(f"{W.SETUP}/core.json").read())


def guard(wins=("V2010", "A2011"), out="guard_station_days.csv"):
    core = core_list()
    rows = []
    for w in wins:
        for d, k in W.window_days(w):
            for ns in core:
                c = W.cache_path(ns, d)
                ok = False; err = ""
                if os.path.exists(c):
                    try:
                        z = np.load(c); _ = z["x"].shape; ok = True
                    except Exception as e:
                        err = repr(e)[:80]
                rows.append(dict(window=w, date=d, kind=k, netsta=ns, raw_on_disk=raw_any(ns, d),
                                 cache=os.path.exists(c), cache_ok=ok, err=err))
    X = pd.DataFrame(rows)
    X.to_csv(f"{GO}/{out}", index=False)
    ond = X[X.raw_on_disk | X.cache]
    bad = ond[~ond.cache_ok]
    fr = len(bad) / max(len(ond), 1)
    print(f"core station-days on the {X.date.nunique()} days {list(wins)}: {len(X)}; on disk (raw or cache) {len(ond)}; "
          f"missing/error cache {len(bad)} ({fr:.2%})")
    print(X.groupby(["netsta", "window"]).cache_ok.sum().unstack().to_string())
    if fr > 0.03:
        open(f"{GO}/GUARD_FAILED", "w").write(f"{len(bad)} of {len(ond)} on-disk core station-days missing/error\n")
        sys.exit(1)


def postguard(tag, gfile="guard_station_days.csv"):
    S = pd.read_csv(f"{W.SETUP}/stations.csv")
    G = pd.read_csv(f"{GO}/{gfile}")
    G = G[G.netsta.isin(S.netsta) & G.cache_ok]
    rows = []
    for p in sorted(glob.glob(f"{GO}/mf_{tag}/run_*.json")):
        r = json.load(open(p))
        for ns in S.netsta:
            rows.append(dict(date=r["date"], netsta=ns, loaded=ns in r["nodes"], dead=ns in r.get("dead", {})))
    L = pd.DataFrame(rows).merge(G[["date", "netsta"]], on=["date", "netsta"], how="inner")
    lost = L[~L.loaded]
    fr = len(lost) / max(len(L), 1)
    print(f"[postguard {tag}] cached core station-days {len(L)}; loaded {int(L.loaded.sum())}; not loaded {len(lost)} "
          f"({fr:.2%}); of the not-loaded, dropped by the T1 rule (> 5 % of the day guarded): {int(L.dead.sum())}; "
          f"not loaded for any other reason (error / missing): {int((~L.loaded & ~L.dead).sum())}")
    fr = int((~L.loaded & ~L.dead).sum()) / max(len(L), 1)
    L.to_csv(f"{GO}/mf_{tag}/postguard.csv", index=False)
    if fr > 0.03:
        open(f"{GO}/GUARD_FAILED", "w").write(f"postguard {tag}: {len(lost)} of {len(L)} not loaded\n")
        sys.exit(1)


if __name__ == "__main__":
    cmd = sys.argv[1]
    rest = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else sys.argv[2:]
    if cmd == "setup":
        setup()
    elif cmd == "guard":
        guard(*([tuple(rest[0].split(",")), rest[1]] if rest else []))
    elif cmd == "postguard":
        postguard(*rest)
    elif cmd in ("mf", "hf"):
        mod = __import__("wa_mf" if cmd == "mf" else "wa_hf")
        sys.argv = [f"wa_{cmd}.py"] + rest
        mod.main()
