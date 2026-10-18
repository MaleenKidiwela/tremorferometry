#!/usr/bin/env python
"""SH_KMRFIX task 2 (2026-10-01, CPU): forward catalogue v2 = segments cont + A (17-station it0 run, unchanged) up to
2025-09-16, then the 16-station no-KMR run mf_NK for 2025-09-17 .. 2026-09-28.  Same Text S1 rules as the original
(sh_v3_common.apply_rules, k = 9, 4-s declustering over all families, applied to the merged set as sh_ext2026.catalog).
Checks: cat_NK_B / cat_NK_C reproduced from mf_NK; the v2 rows before 2025-09-17 equal the original catalogue rows.
Every >= 8x candidate of v2 gets det_id_v2 and is mapped to the nearest cut-window candidate (det_id): exact (cont/A,
NK switch days), else same family within 2 s, else any family within 2 s, else unmapped.
Outputs results/sh_kmrfix/{catalogue_S1_9x_v2.parquet,candidates8x_v2.parquet,t2_checks.json}."""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sh_v3_common import apply_rules
WS = "/home/jovyan/tremorferometry/lfe_south"; RES = f"{WS}/results"; EXT = f"{RES}/sh_ext2026"; CUT = f"{RES}/cut/sh_ext"
OUT = f"{RES}/sh_kmrfix"; os.makedirs(OUT, exist_ok=True)
T0, T_END = pd.Timestamp("2023-12-01"), pd.Timestamp("2026-09-29"); SW1 = pd.Timestamp("2025-09-17")
COLS = ["t", "fam", "ccsum", "meancc", "nlive", "net_mad", "nchan", "nsta", "seg"]


def days(seg):
    d = json.load(open(f"{EXT}/days/{seg}.json")); v = sorted(set(sum(d.values(), [])) if isinstance(d, dict) else d)
    return set(pd.to_datetime(v))


def keys(D):
    return D.fam.astype(str) + "|" + np.round(D.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9 * 20).astype(np.int64).astype(str)


def main():
    parts = []
    for seg, f in [("cont", "mf_cont"), ("A", "mf_A"), ("NK", "mf_NK")]:
        D = pd.read_parquet(f"{EXT}/{f}.parquet"); D["fam"] = D.fam.astype(str); D["seg"] = seg
        ds = days(seg); n0 = len(D); D = D[D.t.dt.floor("D").isin(ds)]
        print(f"[cat] {seg}: {len(D)} of {n0} rows on its {len(ds)} days")
        parts.append(D)
    D = pd.concat(parts, ignore_index=True).sort_values("t").reset_index(drop=True)
    D = D[(D.t >= T0) & (D.t < T_END)].reset_index(drop=True)
    dd = D.groupby([D.t.dt.floor("D"), "seg"]).size().reset_index().groupby("t").size()
    assert (dd == 1).all(), "a day in two segments"
    CAT = apply_rules(D, 9)[COLS].reset_index(drop=True)
    CAT.to_parquet(f"{OUT}/catalogue_S1_9x_v2.parquet", index=False)
    chk = {}
    # check 1: cat_NK_B / C reproduced
    NK = parts[2]
    for ep in ("B", "C"):
        E = apply_rules(NK[NK.t.dt.floor("D").isin(days(ep))].reset_index(drop=True), 9)
        S = pd.read_parquet(f"{EXT}/cat_NK_{ep}.parquet"); S["fam"] = S.fam.astype(str)
        same = len(set(keys(E)) & set(keys(S)))
        inv2 = CAT[CAT.t.dt.floor("D").isin(days(ep))]
        chk[f"cat_NK_{ep}"] = dict(stored=len(S), recomputed=len(E), identical_rows=same, in_v2_catalogue_on_those_days=len(inv2),
                                   v2_rows_identical=len(set(keys(inv2)) & set(keys(S))))
    # check 2: v2 before the swap == original catalogue
    O = pd.read_parquet(f"{EXT}/cat_S1_9x_internal.parquet"); O["fam"] = O.fam.astype(str)
    a = CAT[CAT.t < SW1]; b = O[O.t < SW1]
    chk["pre_swap"] = dict(v2=len(a), original=len(b), identical=len(set(keys(a)) & set(keys(b))))
    b2 = O[O.t >= SW1]; a2 = CAT[CAT.t >= SW1]
    chk["post_swap"] = dict(v2=len(a2), original_v1=len(b2), common_rows=len(set(keys(a2)) & set(keys(b2))),
                            v2_over_v1=len(a2) / len(b2))
    chk["total"] = dict(v2=len(CAT), v1=len(O))
    # >= 8x candidates + det_id_v2 + mapping to the cut windows
    R = D[D.net_mad >= 8.0][COLS].reset_index(drop=True)
    R["t_epoch"] = R.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
    R["det_id_v2"] = np.arange(len(R), dtype=np.int64)
    R["in_S1_9x_v2"] = keys(R).isin(set(keys(CAT))).values
    Dc = pd.read_parquet(f"{CUT}/detections.parquet", columns=["det_id", "fam", "t_epoch", "seg"]); Dc["fam"] = Dc.fam.astype(str)
    Dc = Dc.sort_values("t_epoch").reset_index(drop=True)
    R["det_id"] = -1; R["map_dt_s"] = np.nan; R["map_type"] = "none"; R["map_fam"] = ""
    ta = Dc.t_epoch.to_numpy()
    # same family
    for f, g in R.groupby("fam"):
        c = Dc[Dc.fam == f]; tc = c.t_epoch.to_numpy()
        if not len(tc): continue
        j = np.clip(np.searchsorted(tc, g.t_epoch.to_numpy()), 1, len(tc) - 1)
        d0 = g.t_epoch.to_numpy() - tc[j - 1]; d1 = g.t_epoch.to_numpy() - tc[j]
        k = np.where(np.abs(d1) < np.abs(d0), j, j - 1); dt = np.where(np.abs(d1) < np.abs(d0), d1, d0)
        m = np.abs(dt) <= 2.0
        idx = g.index[m]
        R.loc[idx, "det_id"] = c.det_id.to_numpy()[k[m]]; R.loc[idx, "map_dt_s"] = dt[m]
        R.loc[idx, "map_type"] = np.where(np.abs(dt[m]) < 0.026, "exact", "same_fam"); R.loc[idx, "map_fam"] = f
    # any family for the rest
    rest = R.index[R.map_type == "none"]
    if len(rest):
        tq = R.loc[rest, "t_epoch"].to_numpy()
        j = np.clip(np.searchsorted(ta, tq), 1, len(ta) - 1)
        d0 = tq - ta[j - 1]; d1 = tq - ta[j]
        k = np.where(np.abs(d1) < np.abs(d0), j, j - 1); dt = np.where(np.abs(d1) < np.abs(d0), d1, d0)
        m = np.abs(dt) <= 2.0
        R.loc[rest[m], "det_id"] = Dc.det_id.to_numpy()[k[m]]; R.loc[rest[m], "map_dt_s"] = dt[m]
        R.loc[rest[m], "map_type"] = "other_fam"; R.loc[rest[m], "map_fam"] = Dc.fam.to_numpy()[k[m]]
    R.to_parquet(f"{OUT}/candidates8x_v2.parquet", index=False)
    mt = {}
    for name, sub in [("all_8x", R), ("S1_9x", R[R.in_S1_9x_v2]), ("NK_all_8x", R[R.seg == "NK"]), ("NK_S1_9x", R[(R.seg == "NK") & R.in_S1_9x_v2])]:
        vc = sub.map_type.value_counts().to_dict(); mt[name] = dict(n=len(sub), **{k: int(v) for k, v in vc.items()},
                                                                    frac_mapped_same_fam=float(sub.map_type.isin(["exact", "same_fam"]).mean()),
                                                                    abs_dt_med_same_fam=float(sub[sub.map_type == "same_fam"].map_dt_s.abs().median()),
                                                                    abs_dt_p90_same_fam=float(sub[sub.map_type == "same_fam"].map_dt_s.abs().quantile(0.9)))
    chk["mapping"] = mt
    json.dump(chk, open(f"{OUT}/t2_checks.json", "w"), indent=1, default=float)
    print(json.dumps(chk, indent=1, default=float))


if __name__ == "__main__":
    main()
