#!/usr/bin/env python
"""Trial re-templating input: which station EPOCHS (bridge nodes) are new in the trial window, and the
pass-1 permanent detections to template them from (the bridge rule: a new epoch is templated at an
episode where already-templated sites detect the family).

new node = a node live on trial days whose (station, node) has no build-set stack at all (a station on
the same node as in the validation windows keeps its template: no rebuild).
events   = pass-1 catalogue detections (cat_core: T_core, HF pass, cross-family kept, >= 2 sites).

  python scripts/wa_newnodes.py --mf t2012p1
  -> results/trial/new_nodes.csv, results/trial/pass1_events.parquet
"""
from __future__ import annotations
import os, sys, glob, json, argparse
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mf", default="t2012p1")
    ap.add_argument("--tpl", default="build")
    a = ap.parse_args()
    od = f"{W.RES}/trial"
    os.makedirs(od, exist_ok=True)
    Q = pd.read_csv(f"{W.RES}/templates/tpl_qc_{a.tpl}.csv")
    have = set(zip(Q.netsta, Q.node))
    prior = Q.groupby("netsta").node.apply(lambda v: ";".join(sorted(set(v)))).to_dict()
    live = {}
    for p in glob.glob(f"{W.RES}/mf_{a.mf}/run_*.json"):
        r = json.load(open(p))
        for ns, nd in r["nodes"].items():
            live.setdefault((ns, nd), 0)
            live[(ns, nd)] += 1
    rows = [dict(netsta=ns, node=nd, n_days=n, prior_nodes=prior.get(ns, ""))
            for (ns, nd), n in sorted(live.items()) if (ns, nd) not in have]
    N = pd.DataFrame(rows, columns=["netsta", "node", "n_days", "prior_nodes"])
    N.to_csv(f"{od}/new_nodes.csv", index=False)
    D = pd.read_parquet(f"{W.RES}/score_{a.mf}/detections_flagged.parquet")
    E = D[D.cat_core][["fid", "t"]]
    E.to_parquet(f"{od}/pass1_events.parquet", index=False)
    W.log(f"{W.STUDY}/logs/newnodes.log", f"[newnodes] {len(N)} new nodes: {list(N.node)}; "
          f"{len(E)} pass-1 cat_core events, {E.fid.nunique()} families")


if __name__ == "__main__":
    main()
