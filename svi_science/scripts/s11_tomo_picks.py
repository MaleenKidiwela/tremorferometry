#!/usr/bin/env python
"""S11 tomo pilot: build per-family P/S pick tables with the FROZEN Royer picker (a4f_hypo.make_picks).

Sets (tomo/PLAN.md §0):
  A  NW run-1 picks.csv as on disk (+ eye_class column); also recomputed from its stacks (sanity check)
  B  catalogue Royer_NW_WA families (v1_3), picks recomputed on each family's hypo_28 stack_file
  C  catalogue Royer_SW_WA families, same
Output: tomo/tables/picks_{A,B,C}.csv with columns set, fam, grp, sta, tP, rP, tS, rS, src, + family metadata.
"""
import os
import sys

import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
LFE = f"{ROOT}/lfe_array_discovery/scripts"
sys.path.insert(0, LFE)
os.chdir(f"{ROOT}/lfe_array_discovery")
import a4f_hypo as H  # noqa: E402

OUT = f"{ROOT}/svi_science/results/s11_depth/tomo/tables"
AB = f"{ROOT}/lfe_array_discovery/results/ab_royer"


def load_stacks(npz, grps=None):
    z = np.load(npz)
    ST = {}
    for k in z.files:
        g, s = k.split("_", 1)
        if grps is not None and int(g) not in grps:
            continue
        ST.setdefault(int(g), {})[s] = z[k].astype(np.float64)
    return ST


def main():
    # ---- set A
    A = pd.read_csv(f"{AB}/NW/hypo_28/picks.csv")
    E = pd.read_csv(f"{AB}/NW/hypo_28/review/eye_class.csv")
    A["eye_class"] = A.grp.map(E.set_index("grp")["class"])
    L = pd.read_csv(f"{AB}/NW/hypo_28/loc_headline.csv").set_index("grp")
    L33 = pd.read_csv(f"{AB}/NW/hypo_28/review/lfe33_locations.csv").set_index("grp")
    for c in ("lat", "lon", "dep"):
        A[f"hyp_{c}"] = A.grp.map(L[c])
    for c, d in (("lat", "lat"), ("lon", "lon"), ("depth_km", "dep")):
        A[f"eye_{d}"] = A.grp.map(L33[c])
    A["set"], A["fam"], A["src"] = "A", A.grp.astype(str), "NW/hypo_28/picks.csv"
    # sanity: recompute
    ST = load_stacks(f"{AB}/NW/hypo_28/stacks_loc.npz")
    R = H.make_picks(ST)
    m = A.merge(R, on=["grp", "sta"], suffixes=("", "_re"))
    dd = {c: float(np.nanmax(np.abs(m[c] - m[f"{c}_re"]))) for c in ("tP", "rP", "tS", "rS")}
    nan_mis = int(((m.tP.isna()) != (m.tP_re.isna())).sum() + ((m.tS.isna()) != (m.tS_re.isna())).sum())
    print(f"[A] {len(A)} rows, recompute matched {len(m)}; max |diff| {dd}; NaN mismatches {nan_mis}", flush=True)
    A.to_csv(f"{OUT}/picks_A.csv", index=False)

    # ---- sets B, C
    F = pd.read_csv(f"{ROOT}/data/catalogue/families_v1_3.csv", low_memory=False)
    for tag, reg in (("B", "Royer_NW_WA"), ("C", "Royer_SW_WA")):
        G = F[F.region == reg].copy()
        rows = []
        for sf, g in G.groupby("stack_file"):
            path = f"{ROOT}/{sf}"
            grps = set(int(x) for x in g.family_id)
            ST = load_stacks(path, grps)
            P = H.make_picks(ST)
            miss = grps - set(P.grp.unique()) if len(P) else grps
            if miss:
                print(f"[{tag}] {sf}: {len(miss)} families without stacks: {sorted(miss)[:10]}", flush=True)
            P["src"] = sf
            P = P.merge(g[["family_id", "family_uid", "lat", "lon", "depth_km", "erh", "erz", "n_detections",
                           "era_start", "era_end", "location_quality", "real_criterion"]]
                        .assign(grp=g.family_id.astype(int)), on="grp", how="left")
            rows.append(P)
        P = pd.concat(rows, ignore_index=True)
        P["set"], P["fam"] = tag, P.family_uid
        P.to_csv(f"{OUT}/picks_{tag}.csv", index=False)
        print(f"[{tag}] {P.fam.nunique()} families, {len(P)} station rows; P {P.tP.notna().sum()}, "
              f"S {P.tS.notna().sum()}; rP>=1.2&rS>=1.2 rows {((P.rP >= 1.2) & (P.rS >= 1.2)).sum()}", flush=True)


if __name__ == "__main__":
    main()
