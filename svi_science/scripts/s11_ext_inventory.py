#!/usr/bin/env python
"""s11_ext_inventory -- S11 Revision 2 (PLAN §4): inventory for the SVI families NOT in the S10-multi set (115 - 76),
with the identical logic of s10_multi_inventory.py (same row sources, P >= 0.5 data + c1 rows on active days with
sum P >= 3, same template offsets, cut-store match, per-station expected S and response epochs). The family-selection
rule (>= 10 active days over >= 10 years) is NOT applied: every remaining SVI family with an active day is taken.
Dates map to the existing S10-multi episode units; days outside every unit form new bursts by the same 7-day rule
(labels 'X<first date>'). Read-only on all inputs; writes svi_science/results/s11_depth/tables/ext/.
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s11_common as C

RES = f"{S.ROOT}/svi_science/results"
EXT = f"{C.TB}/ext"
ROW0 = 10_000_000


def load_rows():        # identical to s10_multi_inventory.load_rows
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


def offsets_table(fams):  # identical computation to s10_multi_inventory.offsets_table, written to tables/ext
    os.environ.setdefault("BW_ERA_TPL", f"{RES}/backward/templates/kelb/g3_era_templates.npz")
    os.environ.setdefault("BW_EPOCHS", "1")
    import bw_gate_common as G
    R, Cm = G.patch_forward(backward=True)
    fz = Cm.frozen(); L = float(fz["template_length_s"]); tp, _ = Cm.resolve_templates(fz, "")
    tpl, sp, F, wgt, mv, amp = Cm.load_frozen_tables(fz, tp, L)
    off = Cm.placement(tpl, sp, F, mv, Cm.station_table(), L)
    O = pd.read_csv(f"{S.CUT}/offsets_used.csv")
    rows = []
    for f in fams:
        for ns in C.STATIONS:
            sta = ns.split(".")[1]; k = (float(f), sta)
            o_tpl = off.get(k, np.nan) + 0.5
            ou = O[(O.fam == f) & (O.netsta == ns)]
            o_cut = float(ou.offset_s.iloc[0]) if len(ou) else np.nan
            rows.append(dict(fam=f, netsta=ns, offset_cut=o_cut, offset_cut_src=ou.offset_source.iloc[0] if len(ou) else "",
                             offset_tpl=o_tpl, offset=o_cut if np.isfinite(o_cut) else o_tpl, sp=float(sp.get(k, np.nan))))
    T = pd.DataFrame(rows); T.to_csv(f"{EXT}/offsets_sp.csv", index=False)
    return T


def main():
    os.makedirs(EXT, exist_ok=True)
    log = C.Log(f"{C.LOGD}/ext_inventory.log")
    cat = pd.read_csv(f"{S.ROOT}/data/catalogue/families_v1_3.csv")
    svi = set(pd.to_numeric(cat[cat.region == "S_Vancouver_Is"].family_id, errors="coerce").dropna().astype(int))
    have = set(pd.read_csv(f"{C.MTB}/families.csv").fam)
    want = sorted(svi - have)
    log(f"SVI catalogue families {len(svi)}, in S10 multi {len(have & svi)}, to add {len(want)}: {want}")
    R = load_rows()
    d = R[(R.series == "data") & (R.P >= 0.5)]
    sp = d.groupby(["fam", "date"]).P.sum()
    act = sp[sp >= 3].rename("sumP").reset_index()
    act = act[act.fam.isin(want)].copy()
    g = act.assign(yr=act.date.str[:4]).groupby("fam").agg(n_days=("date", "size"), n_years=("yr", "nunique")).reset_index()
    g.to_csv(f"{EXT}/families.csv", index=False)
    fams = sorted(g.fam.astype(int))
    log(f"families with >= 1 active day: {len(fams)} ({len(want) - len(fams)} have none: {sorted(set(want) - set(fams))}); family-days {len(act)}")
    act["n_data_rows"] = d.groupby(["fam", "date"]).size().reindex(pd.MultiIndex.from_frame(act[["fam", "date"]])).to_numpy()
    O = offsets_table(fams)
    key = set(zip(act.fam, act.date))
    on = np.array([(f, dd) in key for f, dd in zip(R.fam, R.date)])
    X = R[on & (((R.series == "data") & (R.P >= 0.5)) | (R.series == "c1"))].copy().reset_index(drop=True)
    X["pthr"] = np.where(X.detector == "era", 0.90, 0.75)
    X["hiP"] = X.P >= X.pthr
    X["row_id"] = ROW0 + np.arange(len(X))
    D = pd.read_parquet(f"{S.CUT}/detections.parquet", columns=["det_id", "fam", "t"])
    D["fam"] = D.fam.astype(int); D = D[D.fam.isin(fams)].sort_values("t")
    q = X[X.series == "data"][["row_id", "fam", "t"]].copy()
    q["tt"] = pd.to_datetime(q.t, unit="s").astype("datetime64[ns]"); q = q.sort_values("tt")
    m = pd.merge_asof(q, D.rename(columns={"t": "tt_cut"}), left_on="tt", right_on="tt_cut", by="fam",
                      tolerance=pd.Timedelta("0.1s"), direction="nearest")
    X["det_id"] = X.row_id.map(dict(zip(m.row_id, m.det_id))).fillna(-1).astype(np.int64)
    X["pilot_row_id"] = -1
    log(f"rows {len(X)}: data {(X.series == 'data').sum()}, c1 {(X.series == 'c1').sum()}; data in cut store {(X.det_id >= 0).sum()}")
    # episode units: existing map first, then existing unit windows, then new bursts (same 7-day rule)
    A0 = pd.read_csv(f"{C.MTB}/active_days.csv"); U0 = pd.read_csv(f"{C.MTB}/episode_units.csv")
    lab = dict(zip(A0.date, A0.episode))
    dates = sorted(set(act.date))
    for x in dates:
        if x in lab:
            continue
        hit = U0[(U0.w0 <= x) & (U0.w1 >= x)]
        if len(hit):
            lab[x] = hit.episode.iloc[0]
    rest = [x for x in dates if x not in lab]
    newu = []; prev = None
    for x in rest:
        tx = pd.Timestamp(x)
        if prev is None or (tx - prev).days > 7:
            newu.append(dict(episode=f"X{x}", start=x, end=x, cls="burst", big=False, src="s11_ext_bursts", w0=x, w1=x))
        newu[-1]["end"] = x; newu[-1]["w1"] = x; lab[x] = newu[-1]["episode"]; prev = tx
    NU = pd.DataFrame(newu, columns=["episode", "start", "end", "cls", "big", "src", "w0", "w1"])
    NU.to_csv(f"{EXT}/episode_units_new.csv", index=False)
    act["episode"] = act.date.map(lab)
    cls = dict(zip(U0.episode, U0.cls)); cls.update(dict(zip(NU.episode, NU.cls)))
    act["cls"] = act.episode.map(cls)
    act.to_csv(f"{EXT}/active_days.csv", index=False)
    log(f"ext family-days {len(act)}: in existing units {act.episode.isin(set(U0.episode)).sum()}, in majors "
        f"{(act.cls == 'major').sum()}; new burst units {len(NU)}")
    off = O.set_index(["fam", "netsta"]).offset
    rows = []
    for ns in C.STATIONS:
        ep = S.epoch_of(ns, X.date)
        o = np.array([off.get((f, ns), np.nan) for f in X.fam])
        tS = X.t.to_numpy() + o
        keep = (ep != "") & np.isfinite(tS)
        rows.append(pd.DataFrame(dict(row_id=X.row_id.to_numpy()[keep], netsta=ns, epoch=ep[keep], tS=tS[keep],
                                      date=X.date.to_numpy()[keep])))
    RS = pd.concat(rows, ignore_index=True)
    X.to_parquet(f"{EXT}/rows.parquet", index=False); RS.to_parquet(f"{EXT}/row_station.parquet", index=False)
    inv = RS.merge(X[["row_id", "series"]], on="row_id").groupby(["netsta", "series"]).size().unstack(fill_value=0)
    log("row-stations:\n" + inv.to_string())
    open(f"{C.OUT}/EXT_INVENTORY_DONE", "w").write("done\n")


if __name__ == "__main__":
    main()
