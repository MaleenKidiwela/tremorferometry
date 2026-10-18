#!/usr/bin/env python
"""s10_inventory -- S10 pilot step 1 (PLAN §1-3, §9): rows, active days, per-station expected-S times for data and
scrambled-control rows, window sources. Read-only on every input; writes svi_science/results/s10_dvv/pilot/tables.

    BW_ERA_TPL=.../kelb/g3_era_templates.npz BW_EPOCHS=1 python s10_inventory.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S

S.ensure()
log = S.Log(f"{S.LOGD}/inventory.log")


def offsets_table():
    os.environ.setdefault("BW_ERA_TPL", f"{S.ROOT}/svi_science/results/backward/templates/kelb/g3_era_templates.npz")
    os.environ.setdefault("BW_EPOCHS", "1")
    import bw_gate_common as G
    R, C = G.patch_forward(backward=True)
    fz = C.frozen(); L = float(fz["template_length_s"]); tp, _ = C.resolve_templates(fz, "")
    tpl, sp, F, wgt, mv, amp = C.load_frozen_tables(fz, tp, L)
    off = C.placement(tpl, sp, F, mv, C.station_table(), L)
    O = pd.read_csv(f"{S.CUT}/offsets_used.csv")
    O = O[O.fam.isin(S.FAMS)]
    rows = []
    allsta = S.CORE + S.POLARIS
    for f in S.FAMS:
        for ns in allsta:
            sta = ns.split(".")[1]
            k = (float(f), sta)
            o_era = off.get(k, np.nan) + 0.5
            ou = O[(O.fam == f) & (O.netsta == ns)]
            o_cut = float(ou.offset_s.iloc[0]) if len(ou) else np.nan
            src = ou.offset_source.iloc[0] if len(ou) else ""
            use = o_era if (ns in S.POLARIS or not np.isfinite(o_cut)) else o_cut
            rows.append(dict(fam=f, netsta=ns, offset_cut=o_cut, offset_cut_src=src, offset_tpl=o_era, offset=use,
                             sp=float(sp.get(k, np.nan))))
    T = pd.DataFrame(rows)
    T.to_csv(f"{S.TB}/offsets_sp.csv", index=False)
    log("offsets (cut vs template placement, core stations): max |diff| "
        f"{np.nanmax(np.abs((T.offset_cut - T.offset_tpl)[T.offset_cut_src == 'template'])):.4f} s")
    return T


def load_rows():
    fc = ["run", "date", "series", "k", "fam", "P_same", "keep_xfam8", "time_utc"]
    A = pd.read_parquet(f"{S.ROOT}/svi_science/results/s0/analysis_table.parquet", columns=fc,
                        filters=[("fam", "in", list(S.FAMS))])
    A = A[(A.run.astype(str) == "prod") & A.keep_xfam8].copy()
    A = A.assign(run="prod", series=A.series.astype(str), P=A.P_same, t=A.time_utc.astype("datetime64[ns]").astype("int64") / 1e9,
                 detector="forward", src="fw")[["run", "date", "series", "k", "fam", "P", "t", "detector", "src"]]
    B = pd.read_parquet(f"{S.ROOT}/data/catalogue/svi_catalogue_v1_4_backward.parquet",
                        columns=["run", "date", "series", "k", "fam", "P", "t", "detector", "primary", "keep_xfam8"],
                        filters=[("fam", "in", list(S.FAMS))])
    B = B[B.primary & (B.keep_xfam8 == 1)].copy()
    B = B.assign(src="bw", series=B.series.astype(str))[["run", "date", "series", "k", "fam", "P", "t", "detector", "src"]]
    R = pd.concat([B, A], ignore_index=True)
    R["fam"] = R.fam.astype(int); R["k"] = R.k.astype(np.int64)
    return R


def main():
    O = offsets_table()
    R = load_rows()
    d = R[(R.series == "data") & (R.P >= 0.5)]
    sp = d.groupby(["fam", "date"]).P.sum()
    act = sp[sp >= 3].rename("sumP").reset_index()
    act["n_data_rows"] = d.groupby(["fam", "date"]).size().reindex(pd.MultiIndex.from_frame(act[["fam", "date"]])).to_numpy()
    log("active days per family: " + act.groupby("fam").size().to_dict().__repr__())
    key = set(zip(act.fam, act.date))
    on = np.array([(f, dd) in key for f, dd in zip(R.fam, R.date)])
    X = R[on & (((R.series == "data") & (R.P >= 0.5)) | (R.series == "c1"))].copy().reset_index(drop=True)
    X["pthr"] = np.where(X.detector == "era", 0.90, 0.75)
    X["hiP"] = X.P >= X.pthr
    # scrambled-control frame info
    X["i"] = -1; X["tref"] = np.nan; X["NVo"] = -1
    X["shifts"] = ""
    c1 = X.series == "c1"
    for (run, date), g in X[c1].groupby(["run", "date"]):
        base = f"{S.FEAT_FW}/{run}" if run == "prod" else f"{S.FEAT_BW}/{run}"
        Fp = pd.read_parquet(f"{base}/{date}.parquet", columns=["k", "fam", "series", "t", "i"]).set_index("k")
        fd = pd.read_parquet(f"{base}/{date}_fd.parquet")
        fd = fd[fd.series == "c1"].set_index("fam")
        sh = json.load(open(f"{base}/{date}.json"))["shifts"]
        ff = Fp.reindex(g.k)
        assert (ff.fam.to_numpy() == g.fam.to_numpy()).all() and (np.abs(ff.t.to_numpy() - g.t.to_numpy()) < 1e-3).all(), (run, date)
        X.loc[g.index, "i"] = ff.i.to_numpy()
        X.loc[g.index, "tref"] = fd.tref.reindex(g.fam).to_numpy()
        X.loc[g.index, "NVo"] = fd.NVo.reindex(g.fam).to_numpy()
        X.loc[g.index, "shifts"] = json.dumps(sh)
    X["row_id"] = np.arange(len(X))
    # cut-store match for data rows (fam, time within 0.1 s)
    D = pd.read_parquet(f"{S.CUT}/detections.parquet", columns=["det_id", "fam", "t"])
    D["fam"] = D.fam.astype(int); D = D[D.fam.isin(S.FAMS)].sort_values("t")
    q = X[X.series == "data"][["row_id", "fam", "t"]].copy()
    q["tt"] = pd.to_datetime(q.t, unit="s").astype("datetime64[ns]"); q = q.sort_values("tt")
    m = pd.merge_asof(q, D.rename(columns={"t": "tt_cut"}), left_on="tt", right_on="tt_cut", by="fam",
                      tolerance=pd.Timedelta("0.1s"), direction="nearest")
    X["det_id"] = X.row_id.map(dict(zip(m.row_id, m.det_id))).fillna(-1).astype(np.int64)
    log(f"data rows matched to the cut store: {(X.det_id >= 0).sum()} of {(X.series == 'data').sum()} "
        f"(2014+: {((X.det_id >= 0) & (X.date >= '2014')).sum()} of {((X.series == 'data') & (X.date >= '2014')).sum()})")
    # per-station expected S times
    off = O.set_index(["fam", "netsta"])
    rows = []
    for ns in S.CORE + S.POLARIS:
        sta = ns.split(".")[1]; site = S.SITE.get(sta, sta)
        ep = S.epoch_of(ns, X.date)
        ok = ep != ""
        if not ok.any():
            continue
        o = np.array([off.offset.get((f, ns), np.nan) for f in X.fam])
        tS = X.t.to_numpy() + o
        for j in np.flatnonzero(c1.to_numpy() & ok):
            r = X.iloc[j]
            sh = json.loads(r.shifts) if r.shifts else {}
            if site not in sh or r.NVo <= 0:
                tS[j] = np.nan; continue
            s_ = int(sh[site]) % int(r.NVo)
            tS[j] = r.t + (((int(r.i) + s_) % int(r.NVo)) - int(r.i)) / S.FS_DET + o[j]
        keep = ok & np.isfinite(tS)
        rows.append(pd.DataFrame(dict(row_id=X.row_id.to_numpy()[keep], netsta=ns, epoch=ep[keep], tS=tS[keep],
                                      date=X.date.to_numpy()[keep])))
    RS = pd.concat(rows, ignore_index=True)
    X.drop(columns=["shifts"]).to_parquet(f"{S.TB}/rows.parquet", index=False)
    RS.to_parquet(f"{S.TB}/row_station.parquet", index=False)
    act.to_csv(f"{S.TB}/active_days.csv", index=False)
    inv = RS.merge(X[["row_id", "fam", "series"]], on="row_id").groupby(["fam", "netsta", "epoch", "series"]).size().unstack(fill_value=0)
    inv.to_csv(f"{S.TB}/inventory_rows_by_station_epoch.csv")
    log(f"rows {len(X)} (data {(X.series == 'data').sum()}, c1 {c1.sum()}), row-stations {len(RS)}")
    log("row-stations by family/station/epoch:\n" + inv.to_string())


if __name__ == "__main__":
    main()
