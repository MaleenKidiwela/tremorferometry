#!/usr/bin/env python
"""s10_multi_inventory -- S10 multi step 1 (multi/PLAN.md §1-2): select families (>= 10 active days over >= 10 years),
rows on active days (data P >= 0.5 and scrambled-control c1 rows), expected-S times at the 7 stations (data and c2 =
c1 row read at its frame time t + offset), cut-store match, pilot-window reuse keys, shared episode units.
Read-only on every input; writes svi_science/results/s10_dvv/multi/tables."""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_multi_common as M

M.ensure()
log = S.Log(f"{M.LOGD}/inventory.log")
RES = f"{S.ROOT}/svi_science/results"


def load_rows():
    fc = ["run", "date", "series", "k", "fam", "P_same", "keep_xfam8", "time_utc"]
    A = pd.read_parquet(f"{RES}/s0/analysis_table.parquet", columns=fc)
    A = A[(A.run.astype(str) == "prod") & A.keep_xfam8].copy()
    A = A.assign(run="prod", series=A.series.astype(str), P=A.P_same,
                 t=A.time_utc.astype("datetime64[ns]").astype("int64") / 1e9,
                 detector="forward", src="fw")[["run", "date", "series", "k", "fam", "P", "t", "detector", "src"]]
    B = pd.read_parquet(f"{S.ROOT}/data/catalogue/svi_catalogue_v1_4_backward.parquet",
                        columns=["run", "date", "series", "k", "fam", "P", "t", "detector", "primary", "keep_xfam8"])
    B = B[B.primary & (B.keep_xfam8 == 1)].copy()
    B = B.assign(src="bw", series=B.series.astype(str))[["run", "date", "series", "k", "fam", "P", "t", "detector", "src"]]
    R = pd.concat([B, A], ignore_index=True)
    R["fam"] = R.fam.astype(int); R["k"] = R.k.astype(np.int64); R["run"] = R.run.astype(str); R["date"] = R.date.astype(str)
    return R


def offsets_table(fams):
    os.environ.setdefault("BW_ERA_TPL", f"{RES}/backward/templates/kelb/g3_era_templates.npz")
    os.environ.setdefault("BW_EPOCHS", "1")
    import bw_gate_common as G
    R, C = G.patch_forward(backward=True)
    fz = C.frozen(); L = float(fz["template_length_s"]); tp, _ = C.resolve_templates(fz, "")
    tpl, sp, F, wgt, mv, amp = C.load_frozen_tables(fz, tp, L)
    off = C.placement(tpl, sp, F, mv, C.station_table(), L)
    O = pd.read_csv(f"{S.CUT}/offsets_used.csv")
    rows = []
    for f in fams:
        for ns in M.STATIONS:
            sta = ns.split(".")[1]; k = (float(f), sta)
            o_tpl = off.get(k, np.nan) + 0.5
            ou = O[(O.fam == f) & (O.netsta == ns)]
            o_cut = float(ou.offset_s.iloc[0]) if len(ou) else np.nan
            rows.append(dict(fam=f, netsta=ns, offset_cut=o_cut, offset_cut_src=ou.offset_source.iloc[0] if len(ou) else "",
                             offset_tpl=o_tpl, offset=o_cut if np.isfinite(o_cut) else o_tpl, sp=float(sp.get(k, np.nan))))
    T = pd.DataFrame(rows)
    T.to_csv(f"{M.TB}/offsets_sp.csv", index=False)
    # pilot consistency
    P0 = pd.read_csv(f"{M.PILOT}/tables/offsets_sp.csv")
    cmp = T.merge(P0, on=["fam", "netsta"], suffixes=("", "_pilot"))
    log(f"offsets: {len(T)} fam-stations, missing offset {T.offset.isna().sum()}, missing sp {T.sp.isna().sum()}; "
        f"vs pilot max |d offset| {np.nanmax(np.abs(cmp.offset - cmp.offset_pilot)):.2e}, |d sp| {np.nanmax(np.abs(cmp.sp - cmp.sp_pilot)):.2e}")
    return T


def episode_units(dates_nonmajor):
    p = pd.read_csv(f"{RES}/p3_bostock/episodes_bostock_final.csv")
    p = p[p.start < "2010-01-01"][["episode", "start", "end"]].assign(src="p3_bostock")
    b = pd.read_csv(f"{RES}/b_ets_record/episodes_b.csv")
    b = b[(b.region == "SVI") & b.major.astype(bool)][["episode", "start", "end", "big"]].assign(src="episodes_b")
    E = pd.concat([p, b], ignore_index=True)
    E["big"] = E.big.fillna(True).astype(bool)          # all Bostock 2003-09 episodes are big ETS
    E["cls"] = "major"
    E["w0"] = (pd.to_datetime(E.start) - pd.Timedelta(days=2)).dt.strftime("%Y-%m-%d")
    E["w1"] = (pd.to_datetime(E.end) + pd.Timedelta(days=2)).dt.strftime("%Y-%m-%d")
    return E.sort_values("start").reset_index(drop=True)


def assign_units(dates, E):
    d = pd.Series(sorted(set(dates)))
    lab = pd.Series("", index=d.values, dtype=object)
    for r in E.itertuples():
        m = (d >= r.w0) & (d <= r.w1)
        lab[d[m].values] = np.where(lab[d[m].values] == "", r.episode, lab[d[m].values])
    rest = sorted(lab.index[lab == ""])
    bursts = []
    prev = None; cur = None
    for x in rest:
        tx = pd.Timestamp(x)
        if prev is None or (tx - prev).days > 7:
            cur = f"B{x}"; bursts.append(dict(episode=cur, start=x, end=x, cls="burst"))
        bursts[-1]["end"] = x
        lab[x] = cur; prev = tx
    B = pd.DataFrame(bursts)
    U = pd.concat([E[["episode", "start", "end", "cls", "big", "src", "w0", "w1"]], B.assign(big=False, src="bursts", w0=B.start, w1=B.end)],
                  ignore_index=True)
    return lab, U.sort_values("start").reset_index(drop=True)


def main():
    R = load_rows()
    d = R[(R.series == "data") & (R.P >= 0.5)]
    sp = d.groupby(["fam", "date"]).P.sum()
    act = sp[sp >= 3].rename("sumP").reset_index()
    act["yr"] = act.date.str[:4]
    g = act.groupby("fam").agg(n_days=("date", "size"), n_years=("yr", "nunique"))
    sel = g[(g.n_days >= M.MIN_DAYS) & (g.n_years >= M.MIN_YEARS)].reset_index()
    sel.to_csv(f"{M.TB}/families.csv", index=False)
    fams = sorted(sel.fam.astype(int))
    log(f"families selected: {len(fams)} (>= {M.MIN_DAYS} active days over >= {M.MIN_YEARS} years); family-days {sel.n_days.sum()}")
    act = act[act.fam.isin(fams)].drop(columns="yr")
    act["n_data_rows"] = d.groupby(["fam", "date"]).size().reindex(pd.MultiIndex.from_frame(act[["fam", "date"]])).to_numpy()
    O = offsets_table(fams)
    key = set(zip(act.fam, act.date))
    on = np.array([(f, dd) in key for f, dd in zip(R.fam, R.date)])
    X = R[on & (((R.series == "data") & (R.P >= 0.5)) | (R.series == "c1"))].copy().reset_index(drop=True)
    X["pthr"] = np.where(X.detector == "era", 0.90, 0.75)
    X["hiP"] = X.P >= X.pthr
    X["row_id"] = np.arange(len(X))
    # cut-store match (data rows)
    D = pd.read_parquet(f"{S.CUT}/detections.parquet", columns=["det_id", "fam", "t"])
    D["fam"] = D.fam.astype(int); D = D[D.fam.isin(fams)].sort_values("t")
    q = X[X.series == "data"][["row_id", "fam", "t"]].copy()
    q["tt"] = pd.to_datetime(q.t, unit="s").astype("datetime64[ns]"); q = q.sort_values("tt")
    m = pd.merge_asof(q, D.rename(columns={"t": "tt_cut"}), left_on="tt", right_on="tt_cut", by="fam",
                      tolerance=pd.Timedelta("0.1s"), direction="nearest")
    X["det_id"] = X.row_id.map(dict(zip(m.row_id, m.det_id))).fillna(-1).astype(np.int64)
    log(f"rows {len(X)}: data {(X.series == 'data').sum()}, c1 {(X.series == 'c1').sum()}; data in cut store "
        f"{(X.det_id >= 0).sum()} (2014+: {((X.det_id >= 0) & (X.date >= '2014')).sum()} of {((X.series == 'data') & (X.date >= '2014')).sum()})")
    # pilot reuse key
    P0 = pd.read_parquet(f"{M.PILOT}/tables/rows.parquet", columns=["run", "date", "series", "k", "fam", "row_id"])
    P0 = P0.assign(run=P0.run.astype(str), date=P0.date.astype(str), series=P0.series.astype(str)).rename(columns={"row_id": "pilot_row_id"})
    X = X.merge(P0, on=["run", "date", "series", "k", "fam"], how="left")
    X["pilot_row_id"] = X.pilot_row_id.fillna(-1).astype(np.int64)
    pf = X.fam.isin(M.PILOT_FAMS)
    log(f"pilot families: {pf.sum()} rows, matched to pilot rows {(X.pilot_row_id[pf] >= 0).sum()} (pilot table has {len(P0)})")
    # episode units
    E = episode_units(None)
    lab, U = assign_units(act.date, E)
    act["episode"] = act.date.map(lab)
    act = act.merge(U[["episode", "cls"]], on="episode", how="left")
    U.to_csv(f"{M.TB}/episode_units.csv", index=False)
    act.to_csv(f"{M.TB}/active_days.csv", index=False)
    log(f"episode units: majors {(U.cls == 'major').sum()} ({act[act.cls == 'major'].episode.nunique()} with active days), "
        f"bursts {(U.cls == 'burst').sum()}; family-days in majors {(act.cls == 'major').sum()}, in bursts {(act.cls == 'burst').sum()}")
    # per-station expected S (data and c2 alike: t + offset)
    off = O.set_index(["fam", "netsta"]).offset
    rows = []
    for ns in M.STATIONS:
        ep = S.epoch_of(ns, X.date)
        o = np.array([off.get((f, ns), np.nan) for f in X.fam])
        tS = X.t.to_numpy() + o
        keep = (ep != "") & np.isfinite(tS)
        rows.append(pd.DataFrame(dict(row_id=X.row_id.to_numpy()[keep], netsta=ns, epoch=ep[keep], tS=tS[keep],
                                      date=X.date.to_numpy()[keep])))
    RS = pd.concat(rows, ignore_index=True)
    X.to_parquet(f"{M.TB}/rows.parquet", index=False)
    RS.to_parquet(f"{M.TB}/row_station.parquet", index=False)
    inv = RS.merge(X[["row_id", "series"]], on="row_id").groupby(["netsta", "epoch", "series"]).size().unstack(fill_value=0)
    inv.to_csv(f"{M.TB}/inventory_rows_by_station_epoch.csv")
    log("row-stations by station/epoch:\n" + inv.to_string())
    open(f"{M.OUT}/INVENTORY_DONE", "w").write("done\n")


if __name__ == "__main__":
    main()
