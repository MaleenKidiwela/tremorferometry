#!/usr/bin/env python
"""Split the Puget/S-Olympic group list into growth chunks and merge the chunk outputs.  Prefix `n_`.

`k_grow.matched_filter` allocates one (n_templates, 1.728e6) float32 accumulator per day, which
is 10 GB for 1,500 templates; the GPU accumulator is one day-length row per template, so growth is run in
chunks of the SAME sorted order k_grow itself uses (n_nodes, nmad_max descending) and the outputs
concatenated.  Group ids are global, so nothing is renumbered and the merged files are what a
single k_grow call would have produced.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("split")
    s.add_argument("--nodes", required=True)
    s.add_argument("--summary", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--size", type=int, default=1000)
    s.add_argument("--max-groups", type=int, default=0)
    m = sub.add_parser("merge")
    m.add_argument("--dirs", required=True)
    m.add_argument("--out", required=True)
    a = ap.parse_args()

    if a.cmd == "split":
        os.makedirs(a.out, exist_ok=True)
        NO = pd.read_csv(a.nodes)
        G = pd.read_csv(a.summary).sort_values(["n_nodes", "nmad_max"], ascending=False)
        n_all = len(G)
        if a.max_groups and n_all > a.max_groups:
            G = G.head(a.max_groups)
        k = 0
        for i in range(0, len(G), a.size):
            g = G.iloc[i:i + a.size]
            g.to_csv(f"{a.out}/chunk{k:02d}_summary.csv", index=False)
            NO[NO.grp.isin(set(g.grp))].to_csv(f"{a.out}/chunk{k:02d}_nodes.csv", index=False)
            k += 1
        json.dump(dict(n_groups_total=int(n_all), n_groups_grown=int(len(G)),
                       max_groups=int(a.max_groups), n_chunks=k, size=a.size),
                  open(f"{a.out}/_split.json", "w"), indent=1)
        print(f"{n_all} groups, growing {len(G)} in {k} chunk(s) of {a.size}")
        print(k)
    else:
        os.makedirs(a.out, exist_ok=True)
        dirs = [d for d in a.dirs.split(",") if d]
        D = [pd.read_parquet(f"{d}/det_it0.parquet") for d in dirs
             if os.path.exists(f"{d}/det_it0.parquet")]
        if not D:
            raise SystemExit("no chunk detections to merge")
        pd.concat(D, ignore_index=True).to_parquet(f"{a.out}/det_it0.parquet", index=False)
        arr, meta = {}, {}
        for d in dirs:
            p = f"{d}/tpl_it0.npz"
            if not os.path.exists(p):
                continue
            z = np.load(p, allow_pickle=True)
            meta.update(json.loads(str(z["meta"])))
            for k in z.files:
                if k != "meta":
                    arr[k] = z[k]
        np.savez_compressed(f"{a.out}/tpl_it0.npz", **arr,
                            meta=np.array(json.dumps(meta)))
        print(f"merged {len(dirs)} chunks -> {len(meta)} template-stations, "
              f"{sum(len(x) for x in D)} detections")


if __name__ == "__main__":
    main()
