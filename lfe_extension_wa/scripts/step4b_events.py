#!/usr/bin/env python
"""Step 4b: templating inputs for the bridge rule.

--mode array     source (a): the families' OWN array-era detections, on days the added station epochs are
                 live: the V2010 BUILD half (the test half stays untouched for validation) + every A2011
                 detection (2011-08 ETS, NW runs; never scored).  Nodes = every epoch of the added
                 stations + every epoch of a trial station that has no build stack (e.g. GNW BH40,
                 SQM EH100, WISH BH40, NLWA 2011, SP2 2011), as seen in the 20 Hz cache.
                 -> results/step4b/events_array.parquet, results/step4b/nodes_array.csv
--mode episode   source (b), for episode E (T2012 / T2013 / T2014), after its pass 1:
                 source families = the families that SEPARATE: >= 3 sites and data/control >= 2 at the
                 frozen T_core (tier counts per family-day, the validation_by_family.csv definition) in
                 the trial's validation table (V2008/V2010), in the step-4b V2010 validation, or in any
                 step-4b trial pass so far (tremor days of the window; the bridge rule iterates);
                 events = their catalogue detections (cat_core: T_core, HF pass, cross-family kept) on the
                 region's tremor days of E, accumulated with the previous episodes' events;
                 nodes = every epoch live in E with no BUILD stack (added stations + new epochs).
                 -> results/step4b/events_<E>.parquet, nodes_<E>.csv, sources_<E>.csv
"""
from __future__ import annotations
import os, sys, glob, json, argparse
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

PFX = os.environ.get("STEP4B_PREFIX", "step4b_")   # smoke: step4b_smoke_
OUT = f"{W.RES}/step4b" + ("" if PFX == "step4b_" else "/smoke")
os.makedirs(OUT, exist_ok=True)
LOG = f"{W.STUDY}/logs/step4b_events.log"
TREM = {"Royer_NW_WA": "NW_WA", "Royer_SW_WA": "SW_WA", "Olympic_2010": "OLY2010", "Olympic_2010_armB": "OLY2010"}
EPISODES = ["T2012", "T2013", "T2014"]


def tremor_days():
    out = {}
    for rg, tr in TREM.items():
        t = pd.read_csv(f"{W.EXT}/{tr}_tremor_days.csv")
        out[rg] = set(pd.to_datetime(t[t.columns[0]]).dt.strftime("%Y-%m-%d")[t[t.columns[1]] > 0])
    return out


def build_nodes():
    return set(pd.read_csv(f"{W.RES}/templates/tpl_qc_build.csv").node)


def sep_from_validation(p, windows=None):
    V = pd.read_csv(p)
    if windows:
        V = V[V.window.isin(windows)]
    V = V[V.n_days_obs > 0]
    ok = (V.med_sites >= 3) & (V.data_per_famday_core >= 2 * V.ctl_per_famday_core) & (V.data_per_famday_core > 0)
    return set(V[ok].fid.astype(int))


def sep_from_trial(mf):
    """families with median >= 3 observable sites and data/control >= 2 (frozen T_core tier counts) on the
    window's tremor days"""
    FD = pd.concat([pd.read_parquet(p) for p in glob.glob(f"{W.RES}/mf_{mf}/famday_*.parquet")])
    F = W.load_families()
    td = tremor_days()
    reg = dict(zip(F.fid, F.region))
    FD = FD[FD.observable & FD.fid.isin(reg)]
    FD = FD[[d in td[reg[f]] for f, d in zip(FD.fid, FD.date)]]
    g = FD.groupby("fid").agg(ms=("n_sites", "median"), nd=("n_core_data", "sum"), nc=("n_core_control", "sum"))
    return set(g[(g.ms >= 3) & (g.nd > 0) & (g.nd >= 2 * g.nc)].index.astype(int))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["array", "episode"], required=True)
    ap.add_argument("--episode", default="")
    a = ap.parse_args()
    S = W.load_stations()
    new_ns = set(S[S.why.astype(str).str.startswith("step4b")].netsta)
    bn = build_nodes()
    M = pd.read_csv(f"{W.RES}/cache_manifest.csv")
    M = M[M.status.isin(["ok", "cached"])]
    if a.mode == "array":
        T = pd.read_parquet(f"{W.SETUP}/truth.parquet")
        E = T[((T.window == "V2010") & (T.half == 0)) | (T.window == "A2011")][["fid", "t", "window"]]
        E.to_parquet(f"{OUT}/events_array.parquet", index=False)
        days = {d for w in ["V2010", "A2011"] for d, _ in W.window_days(w)}
        m = M[M.date.isin(days) & M.netsta.isin(set(S.netsta))]
        N = m[m.netsta.isin(new_ns) | ~m.node.isin(bn)].groupby(["netsta", "node"]).size().rename("n_days").reset_index()
        N.to_csv(f"{OUT}/nodes_array.csv", index=False)
        W.log(LOG, f"[events:array] {len(E)} events ({E.window.value_counts().to_dict()}), {E.fid.nunique()} families; "
              f"{len(N)} nodes on {N.netsta.nunique()} stations")
        return
    ep = a.episode
    k = EPISODES.index(ep)
    src = sep_from_validation(f"{W.RES}/score_val/validation_by_family.csv")
    p = f"{W.RES}/score_{PFX}v2010/validation_by_family.csv"
    if os.path.exists(p):
        src |= sep_from_validation(p, ["V2010"])
    F = W.load_families()
    src &= set(F.fid)
    trial_src = set()
    for e in EPISODES[:k + 1]:
        for mf in [f"{PFX}{e.lower()}p1", f"{PFX}{e.lower()}p2"]:
            if glob.glob(f"{W.RES}/mf_{mf}/famday_*.parquet"):
                trial_src |= sep_from_trial(mf)
    fam = sorted(src | trial_src)
    pd.DataFrame(dict(fid=fam, from_validation=[f in src for f in fam], from_trial=[f in trial_src for f in fam])
                 ).merge(F[["fid", "family_uid", "region"]], on="fid").to_csv(f"{OUT}/sources_{ep}.csv", index=False)
    td = tremor_days()
    reg = dict(zip(F.fid, F.region))
    D = pd.read_parquet(f"{W.RES}/score_{PFX}{ep.lower()}p1/detections_flagged.parquet")
    D = D[D.cat_core & D.fid.isin(fam)]
    D = D[[d in td[reg[f]] for f, d in zip(D.fid, D.date)]]
    E = D[["fid", "t"]].assign(window=ep)
    if k > 0:
        prev = f"{OUT}/events_{EPISODES[k - 1]}.parquet"
        if os.path.exists(prev):
            E = pd.concat([pd.read_parquet(prev), E], ignore_index=True)
        else:
            W.log(LOG, f"[events:{ep}] WARNING no previous-episode events {prev}")
    E.to_parquet(f"{OUT}/events_{ep}.parquet", index=False)
    nodes = {}
    for pj in glob.glob(f"{W.RES}/mf_{PFX}{ep.lower()}p1/run_*.json"):
        r = json.load(open(pj))
        for ns, nd in r["nodes"].items():
            if ns in new_ns or nd not in bn:
                nodes[(ns, nd)] = nodes.get((ns, nd), 0) + 1
    N = pd.DataFrame([dict(netsta=ns, node=nd, n_days=n) for (ns, nd), n in sorted(nodes.items())],
                     columns=["netsta", "node", "n_days"])
    N.to_csv(f"{OUT}/nodes_{ep}.csv", index=False)
    W.log(LOG, f"[events:{ep}] source families {len(fam)} (validation {len(src)}, trial {len(trial_src)}); events "
          f"{len(E)} ({E.window.value_counts().to_dict()}); {len(N)} nodes without a build stack")


if __name__ == "__main__":
    main()
