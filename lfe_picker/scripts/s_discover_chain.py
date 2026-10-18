#!/usr/bin/env python
"""P3D DISCOVERY -- the chain from a pair-edge list to gated, located candidate families.

Merlin P3D Q3: kNN -> connected components -> align_stack -> S excess >= +3 dB at >= 3 sites ->
g_locate 2-D interface -> 12-s templates.  Gates 0-3 and the scorer are applied afterwards by
`s_discover_score.py`; this script produces the proposals.

Inputs : results/p3d/disc_<tag>_edges.parquet, disc_<tag>_windows.parquet, disc_<tag>.json
Outputs: results/p3d/disc_<tag>_clusters.csv    one row per surviving component
         results/p3d/disc_<tag>_proposals.parquet   g_gates-ready (pid, t, lat, lon)
         results/p3d/disc_<tag>_bank.npz        12-s templates in the frozen format
"""
from __future__ import annotations
import os, sys, json, time, argparse, warnings
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                       # noqa: E402
import s_transfer as T                                                     # noqa: E402
import s_discover_pairs as DP                                              # noqa: E402

OUT = S.P3D


def components(E, n_nodes, knn=10, min_members=20, max_members=20000):
    """kNN sparsification (each node keeps its strongest `knn` edges) then connected components."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    if not len(E):
        return np.full(n_nodes, -1), pd.DataFrame()
    keep = np.zeros(len(E), bool)
    for col in ("i", "j"):
        o = np.lexsort((-E.z.values, E[col].values))
        seen = {}
        for idx in o:
            k = E[col].values[idx]
            c = seen.get(k, 0)
            if c < knn:
                keep[idx] = True
                seen[k] = c + 1
    Ek = E[keep]
    g = coo_matrix((np.ones(len(Ek)), (Ek.i.values, Ek.j.values)), shape=(n_nodes, n_nodes))
    g = g + g.T
    n_c, lab = connected_components(g, directed=False)
    sz = np.bincount(lab, minlength=n_c)
    good = np.flatnonzero((sz >= min_members) & (sz <= max_members))
    rows = [dict(cid=int(c), n_members=int(sz[c])) for c in good]
    # A loose threshold can fuse everything into one giant component; that is a RESULT (the graph
    # has percolated and the threshold is too permissive for this network), not an error, so the
    # size distribution is reported and an empty, well-typed table is returned.
    big = int((sz > max_members).sum())
    print(f"  [components] {n_c} components; sizes: max {int(sz.max()) if n_c else 0}, "
          f"{int((sz >= min_members).sum())} with >= {min_members} members, "
          f"{big} ABOVE the {max_members} cap (percolated)", flush=True)
    if not rows:
        return lab, pd.DataFrame(columns=["cid", "n_members"]).astype({"cid": int, "n_members": int})
    return lab, pd.DataFrame(rows).sort_values("n_members", ascending=False).reset_index(drop=True)


def site_stacks(dc, sites_st, bkeyless_times, cst, ST, max_ev=400, rng=None):
    """Per-SITE aligned stack at the pair statistic's own placement (grid-centre moveout).
    Returns {site: dict(stack, s_idx, excess_db, snr_db, n)}."""
    out = {}
    t_ev = np.asarray(bkeyless_times, float)
    if rng is not None and len(t_ev) > max_ev:
        t_ev = np.sort(rng.choice(t_ev, max_ev, replace=False))
    for site, stas in sites_st.items():
        best = None
        for s in stas:
            Wd, keep = S.cut_fast(dc, s, t_ev + cst[s] - T.WIDE, T.NWIDE, norm=False)
            if Wd is None or len(Wd) < 5:
                continue
            W = T.win_at(Wd, 0)
            if W is None:
                continue
            st_un = W.mean(0)
            st_al, W_al, sh = S.align_stack(W, iters=2, lag=T.OFF_MAX)
            s_idx = S.pick_s(st_al, T.NOM, search_s=4.0)
            rec = dict(site=site, sta=s, n=len(W), s_idx=s_idx,
                       excess_db=T.excess_offsearch(st_un),
                       snr_db=S.s_snr_db(st_al, s_idx),
                       rho=T.loo_registration(W_al, st_al, s_idx),
                       stack=st_al.astype(np.float32),
                       dt_s=(s_idx - T.NOM) / S.FS + cst[s])
            if best is None or (np.nan_to_num(rec["excess_db"], nan=-99) >
                                np.nan_to_num(best["excess_db"], nan=-99)):
                best = rec
        if best is not None:
            out[site] = best
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--knn", type=int, default=10)
    ap.add_argument("--min-members", type=int, default=20)
    ap.add_argument("--min-excess-db", type=float, default=3.0)
    ap.add_argument("--min-sites", type=int, default=3)
    ap.add_argument("--max-clusters", type=int, default=400)
    ap.add_argument("--max-events", type=int, default=400)
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    t0 = time.time()
    meta = json.load(open(f"{OUT}/disc_{a.tag}.json"))
    Wdf = pd.read_parquet(f"{OUT}/disc_{a.tag}_windows.parquet")
    E = pd.read_parquet(f"{OUT}/disc_{a.tag}_edges.parquet")
    print(f"[chain] {len(Wdf)} windows, {len(E)} edges", flush=True)
    lab, comp = components(E, len(Wdf), a.knn, a.min_members)
    print(f"[chain] {len(comp)} components with >= {a.min_members} members "
          f"(largest {comp.n_members.max() if len(comp) else 0})", flush=True)
    if not len(comp):
        pd.DataFrame(columns=["cid", "n_members", "n_sites_stacked", "n_sites_excess",
                              "excess_by_site", "pass_excess", "t_first", "t_last",
                              "span_days"]).to_csv(f"{OUT}/disc_{a.tag}_clusters.csv", index=False)
        pd.DataFrame(columns=["cid", "t"]).to_parquet(f"{OUT}/disc_{a.tag}_members.parquet",
                                                      index=False)
        json.dump(dict(tag=a.tag, n_components=0, note="no component reached the member floor"),
                  open(f"{OUT}/disc_{a.tag}_chain.json", "w"), indent=1)
        print("[chain] nothing to carry forward")
        return
    ST = pd.read_csv(f"{S.P3D}/stations_p3d_v2.csv").set_index("sta")
    stas = meta["stations"]
    site_of = {s: ST.site[s] for s in stas}
    sites_st = {}
    for s in stas:
        sites_st.setdefault(site_of[s], []).append(s)
    cst, ctr, cdep = DP.station_constants(stas, ST)
    dc = S.DayCache(maxn=700, years=list(range(2019, 2027)), stas=stas)
    import g_locate as GL
    import g_common as GC
    site_ll = {q: (float(np.mean([ST.lat[s] for s in v])), float(np.mean([ST.lon[s] for s in v])))
               for q, v in sites_st.items()}
    grid = GL.Grid(centre=ctr, sites=sorted(site_ll), site_ll=site_ll) \
        if "centre" in GL.Grid.__init__.__code__.co_varnames else None
    rng = np.random.default_rng(20260919)
    rows, props, tpl = [], [], {}
    for r in comp.head(a.max_clusters).itertuples():
        mem = np.flatnonzero(lab == r.cid)
        tt = Wdf.t.values[mem]
        stk = site_stacks(dc, sites_st, tt, cst, ST, a.max_events, rng)
        ok = {q: v for q, v in stk.items()
              if np.isfinite(v["excess_db"]) and v["excess_db"] >= a.min_excess_db}
        rec = dict(cid=int(r.cid), n_members=int(r.n_members), n_sites_stacked=len(stk),
                   n_sites_excess=len(ok),
                   excess_by_site=";".join(f"{q}:{v['excess_db']:.2f}" for q, v in sorted(stk.items())),
                   pass_excess=bool(len(ok) >= a.min_sites),
                   t_first=float(tt.min()), t_last=float(tt.max()),
                   span_days=float((tt.max() - tt.min()) / 86400.0))
        if rec["pass_excess"]:
            tS = {q: stk[q]["dt_s"] for q in ok}
            rec["tS"] = json.dumps({k: round(v, 3) for k, v in tS.items()})
            props.append(dict(cid=int(r.cid), n=len(mem)))
            for q in ok:
                v = stk[q]
                s_idx = v["s_idx"]
                aa = s_idx - int(round(T.TPL_PRE * S.FS))
                if 0 <= aa and aa + T.NTPL <= v["stack"].shape[1]:
                    z = v["stack"][:, aa:aa + T.NTPL].astype(np.float64)
                    n = np.linalg.norm(z, axis=1, keepdims=True)
                    tpl[(float(r.cid), v["sta"])] = {12.0: (z / np.where(n > 0, n, 1)).astype(np.float32)}
        rows.append(rec)
        if len(rows) % 25 == 0:
            print(f"  [chain] {len(rows)}/{min(len(comp), a.max_clusters)}  {time.time()-t0:.0f}s",
                  flush=True)
    C = pd.DataFrame(rows)
    C.to_csv(f"{OUT}/disc_{a.tag}_clusters.csv", index=False)
    np.savez_compressed(f"{OUT}/disc_{a.tag}_bank.npz", tpl=np.array(tpl, dtype=object),
                        sp=np.array({k: 0.0 for k in tpl}, dtype=object))
    # member times for the surviving clusters, for the gates / scorer
    keep = set(C[C.pass_excess].cid)
    M = pd.DataFrame(dict(cid=lab[np.isin(lab, list(keep))] if keep else [],
                          t=Wdf.t.values[np.isin(lab, list(keep))] if keep else []))
    M.to_parquet(f"{OUT}/disc_{a.tag}_members.parquet", index=False)
    json.dump(dict(tag=a.tag, n_components=int(len(comp)), n_scored=int(len(C)),
                   n_pass_excess=int(C.pass_excess.sum()), min_excess_db=a.min_excess_db,
                   min_sites=a.min_sites, knn=a.knn, min_members=a.min_members,
                   runtime_s=time.time() - t0),
              open(f"{OUT}/disc_{a.tag}_chain.json", "w"), indent=1, default=str)
    print(f"[chain] {int(C.pass_excess.sum())}/{len(C)} components pass the S-excess clause "
          f"(>= {a.min_excess_db} dB at >= {a.min_sites} sites)  {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
