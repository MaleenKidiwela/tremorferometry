#!/usr/bin/env python
"""Score one stack_activity run: windows of W days, stack correlation, z vs random times, p vs control.

For W in {1, 5, 10, 30} (calendar blocks of W days from the first day of each window label; a continuous
run uses its first day): each (family, station-epoch) stack of the window = sum of its daily stacks.
  rho       stack-template correlation, 2-8 Hz: mean over the template's live channels of
            <S_c, tau_c> / (|S_c| |tau_c|)
  u         noise-unit form: mean_c <S_c, tau_c/|tau_c|> / sqrt(e_c / L)  (e = summed window energy)
  z         (rho - mean of the 6 random-times stacks of that window) / pooled SD of the random-times
            stacks (RMS over all windows of that station-template and W) — for data and each control
  A         family activity = Stouffer sum over stations  sum z / sqrt(n)  (all templates; A_clear = clear
            templates only), computed identically for data and control1..3
  p_ctl     empirical p of A against the scrambled-control pool of that family and W (control1..3, all
            windows); each control is scored against the other two controls the same way
Validation mode: truth = the family's own array-era detections in the window, test half when the run is
--scored test, counted on days its source run scanned; W > 1 windows contain only scanned days (W = 1
rows keep unscanned days with truth = NaN). Tremor days (PNSN region tables) are a LABEL column only.

  python scripts/stack_activity_score.py --run valA --mode validation --truth results/setup/truth.parquet \
         --scan results/setup/scan_days.csv
"""
from __future__ import annotations
import os, sys, glob, json, argparse, time
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W
from stack_activity import load_templates, SERIES, OUT

WS = [1, 5, 10, 30]
TREM = {"Royer_NW_WA": "NW_WA", "Royer_SW_WA": "SW_WA", "Olympic_2010": "OLY2010", "Olympic_2010_armB": "OLY2010"}
LOG = f"{W.STUDY}/logs/stack_activity.log"


def label_of(ds):
    for w, segs in W.WINDOWS.items():
        for a, b, _ in segs:
            if a <= ds <= b:
                return w
    return None


def tremor_days():
    out = {}
    for rg, tr in TREM.items():
        t = pd.read_csv(f"{W.EXT}/{tr}_tremor_days.csv")
        out[rg] = set(pd.to_datetime(t[t.columns[0]]).dt.strftime("%Y-%m-%d")[t[t.columns[1]] > 0])
    return out


def rho_u(S, e, T):
    """S (N,10,3,L) f32, e (N,10,3), T (N,3,L) unit per live channel -> rho, u (N,10)"""
    live = (np.abs(T).sum(-1) > 0)                                  # (N,3)
    dot = np.einsum("nscl,ncl->nsc", S, T)
    nS = np.sqrt((S ** 2).sum(-1))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(nS > 0, dot / nS, np.nan)
        u = np.where(e > 0, dot / np.sqrt(e / S.shape[-1]), np.nan)
    lv = live[:, None, :]
    r = np.where(lv, r, np.nan); u = np.where(lv, u, np.nan)
    with np.errstate(invalid="ignore"):
        return np.nanmean(r, -1), np.nanmean(u, -1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--mode", choices=["validation", "trial"], required=True)
    ap.add_argument("--truth", default="")
    ap.add_argument("--scan", default="")
    ap.add_argument("--families", default=f"{W.RES}/setup/families.csv")
    ap.add_argument("--ws", default="1,5,10,30")
    a = ap.parse_args()
    t0 = time.time()
    rd = f"{OUT}/run_{a.run}"
    cfg = json.load(open(f"{rd}/config.json"))
    od = f"{OUT}/score_{a.run}"
    os.makedirs(od, exist_ok=True)
    F = pd.read_csv(a.families)
    REG = dict(zip(F.fid, F.region))
    Q, _, T28 = load_templates(cfg["tpl"].split(","))
    files = sorted(glob.glob(f"{rd}/day_*.npz"))
    days = [os.path.basename(p)[4:14] for p in files]
    cont = bool(cfg.get("date_range"))
    origin = {}
    for ds in days:
        lab = "cont" if cont else label_of(ds)
        origin.setdefault(lab, ds)
    scored_test = cfg["scored"] == "test"
    # scanned days per family (validation)
    scan = None
    if a.mode == "validation":
        SC = pd.read_csv(a.scan)
        ds_by_src = SC.groupby("det_source").date.apply(set).to_dict()
        scan = {r.fid: ds_by_src.get(r.det_source, set()) for r in F.itertuples()}
        TR = pd.read_parquet(a.truth)
        if scored_test:
            TR = TR[TR.half == 1]
        TRc = TR.groupby(["fid", "date"]).size()
        TRt = {f: np.sort(g.t.to_numpy()) for f, g in TR.groupby("fid")}
    TD = tremor_days()
    ws = [int(x) for x in a.ws.split(",")]
    rows = []                   # station-window rows
    # accumulators for W > 1: key (W, lab, block) -> dict rowkey -> [S, e, cnt, ndays]
    acc = {w: {} for w in ws if w > 1}
    cur = {w: None for w in ws if w > 1}
    days_in = {}

    def flush(w, key):
        A = acc[w].pop(key, None)
        if not A:
            return
        keys = list(A)
        S = np.stack([A[k][0] for k in keys]); e = np.stack([A[k][1] for k in keys])
        c = np.stack([A[k][2] for k in keys])
        T = np.stack([T28[k] for k in keys])
        r, u = rho_u(S, e, T)
        for i, k in enumerate(keys):
            rows.append((w, key[1], key[2], k[0], k[1], k[2], A[k][3], c[i], r[i], u[i]))

    for ip, (p, ds) in enumerate(zip(files, days)):
        z = np.load(p)
        fid = z["fid"]; ns = z["netsta"]; node = z["node"]
        S = z["S"].astype(np.float32); e = z["e"]; c = z["cnt"]
        keys = [(int(f), s, nd) for f, s, nd in zip(fid, ns, node)]
        lab = "cont" if cont else label_of(ds)
        dnum = (pd.Timestamp(ds) - pd.Timestamp(origin[lab])).days
        if 1 in ws and len(keys):
            T = np.stack([T28[k] for k in keys])
            r, u = rho_u(S, e, T)
            for i, k in enumerate(keys):
                rows.append((1, lab, dnum, k[0], k[1], k[2], 1, c[i], r[i], u[i]))
        days_in.setdefault((1, lab, dnum), []).append(ds)
        for w in acc:
            key = (w, lab, dnum // w)
            if cur[w] is not None and cur[w] != key:
                flush(w, cur[w])
            cur[w] = key
            days_in.setdefault(key, []).append(ds)
            B = acc[w].setdefault(key, {})
            for i, k in enumerate(keys):
                if scan is not None and ds not in scan.get(k[0], ()):
                    continue
                if k in B:
                    B[k][0] += S[i]; B[k][1] += e[i]; B[k][2] += c[i]; B[k][3] += 1
                else:
                    B[k] = [S[i].copy(), e[i].astype(np.float64), c[i].astype(np.int64), 1]
        if (ip + 1) % 10 == 0:
            W.log(LOG, f"[score:{a.run}] {ip+1}/{len(files)} days {time.time()-t0:.0f}s")
    for w in acc:
        if cur[w] is not None:
            flush(w, cur[w])
    ns_ = len(SERIES)
    R = pd.DataFrame(dict(W=[x[0] for x in rows], label=[x[1] for x in rows], block=[x[2] for x in rows],
                          fid=[x[3] for x in rows], netsta=[x[4] for x in rows], node=[x[5] for x in rows],
                          n_days_stacked=[x[6] for x in rows]))
    C = np.stack([x[7] for x in rows]); RH = np.stack([x[8] for x in rows]); U = np.stack([x[9] for x in rows])
    for j, s in enumerate(SERIES):
        if not s.startswith("random") or s == "random1":
            R[f"cnt_{s}"] = C[:, j]
        R[f"rho_{s}"] = RH[:, j]
    for j, s in enumerate(SERIES[:4]):
        R[f"u_{s}"] = U[:, j]
    rnd = RH[:, 4:]
    R["rand_mean"] = np.nanmean(rnd, 1)
    R["rand_var"] = np.nanvar(rnd, 1, ddof=1)
    R["u_rand_mean"] = np.nanmean(U[:, 4:], 1)
    QI = Q.set_index(["fid", "netsta", "node"])
    mi = pd.MultiIndex.from_arrays([R.fid, R.netsta, R.node])
    R["clear"] = QI.clear.reindex(mi).fillna(False).to_numpy().astype(bool)
    # network weight of the template (wa_mf rule: n_channels x SNR weight) for the weighted Stouffer A_w
    R["wt"] = (QI.w * QI.live.astype(str).str.count("1")).reindex(mi).to_numpy()
    grp = R.groupby(["W", "fid", "netsta", "node"]).rand_var
    R["rand_sd_pool"] = np.sqrt(grp.transform("mean"))
    for s in SERIES[:4]:
        R[f"z_{s}"] = (R[f"rho_{s}"] - R.rand_mean) / R.rand_sd_pool
    # window metadata
    R["start"] = [min(days_in[(w, l, b)]) for w, l, b in zip(R.W, R.label, R.block)]
    R["end"] = [max(days_in[(w, l, b)]) for w, l, b in zip(R.W, R.label, R.block)]
    R.to_parquet(f"{od}/station_window.parquet", index=False)
    # ---------------- family windows
    ok = R.cnt_data > 0
    out = []
    for (w, lab, b, f), g in R[ok].groupby(["W", "label", "block", "fid"]):
        row = dict(W=w, label=lab, block=b, fid=f, region=REG.get(f), start=g.start.iloc[0], end=g.end.iloc[0],
                   n_sta=len(g), n_sta_clear=int(g.clear.sum()), n_days_stacked=int(g.n_days_stacked.max()),
                   cnt_data=int(g.cnt_data.sum()))
        for s in SERIES[:4]:
            zz = g[f"z_{s}"].to_numpy(); m = np.isfinite(zz)
            row[f"A_{s}"] = zz[m].sum() / np.sqrt(m.sum()) if m.any() else np.nan
            mc = m & g.clear.to_numpy()
            row[f"A_clear_{s}"] = zz[mc].sum() / np.sqrt(mc.sum()) if mc.any() else np.nan
            wt = g.wt.to_numpy()
            mw = m & np.isfinite(wt)
            row[f"A_w_{s}"] = (wt[mw] * zz[mw]).sum() / np.sqrt((wt[mw] ** 2).sum()) if mw.any() else np.nan
            row[f"sumrho_{s}"] = float(np.nansum(g[f"rho_{s}"]))
        row["sumrho_rand_mean"] = float(np.nansum(g.rand_mean))
        out.append(row)
    FW = pd.DataFrame(out)
    # p vs the scrambled-control pool (per run, W, family); each control vs the other two
    for pre in ["A", "A_w", "A_clear"]:
        for s in SERIES[:4]:
            FW[f"p_{pre}_{s}"] = np.nan
        for (w, f), g in FW.groupby(["W", "fid"]):
            ctl = {s: g[f"{pre}_{s}"].to_numpy() for s in SERIES[1:4]}
            pool = np.concatenate([v for v in ctl.values()]); pool = pool[np.isfinite(pool)]
            x = g[f"{pre}_data"].to_numpy()
            FW.loc[g.index, f"p_{pre}_data"] = [(1 + (pool >= v).sum()) / (1 + len(pool)) if np.isfinite(v) and len(pool) else np.nan for v in x]
            for s in SERIES[1:4]:
                pl = np.concatenate([ctl[t] for t in SERIES[1:4] if t != s]); pl = pl[np.isfinite(pl)]
                FW.loc[g.index, f"p_{pre}_{s}"] = [(1 + (pl >= v).sum()) / (1 + len(pl)) if np.isfinite(v) and len(pl) else np.nan for v in ctl[s]]
    # labels: truth (validation), tremor days
    dl = {k: v for k, v in days_in.items()}
    FW["days"] = [";".join(sorted(dl[(w, l, b)])) for w, l, b in zip(FW.W, FW.label, FW.block)]
    FW["n_tremor_days"] = [len(set(d.split(";")) & TD.get(r, set())) for d, r in zip(FW.days, FW.region)]
    FW["n_days"] = FW.days.str.count(";") + 1
    if a.mode == "validation":
        nt, nsc = [], []
        for d, f in zip(FW.days, FW.fid):
            dd = [x for x in d.split(";") if x in scan.get(f, ())]
            nsc.append(len(dd))
            nt.append(int(sum(TRc.get((f, x), 0) for x in dd)) if dd else np.nan)
        FW["n_scanned_days"] = nsc
        FW["n_truth"] = nt
        # diagnostic: fraction of the full-network top-K candidates within +-2 s of a truth detection
        cands = []
        for ds in days:
            p = f"{rd}/cand_{ds}.parquet"
            if os.path.exists(p):
                cands.append(pd.read_parquet(p).assign(date=ds))
        if cands:
            CA = pd.concat(cands, ignore_index=True)
            CA["hit"] = False
            for f, g in CA.groupby("fid"):
                tt = TRt.get(f)
                if tt is None or not len(tt):
                    continue
                i = np.searchsorted(tt, g.t.to_numpy())
                d1 = np.abs(g.t.to_numpy() - tt[np.clip(i - 1, 0, len(tt) - 1)])
                d2 = np.abs(g.t.to_numpy() - tt[np.clip(i, 0, len(tt) - 1)])
                CA.loc[g.index, "hit"] = np.minimum(d1, d2) <= W.MATCH_S
            H = CA.groupby(["fid", "date", "series"]).hit.mean().unstack("series")
            m1 = FW.W == 1
            for s in H.columns:
                FW.loc[m1, f"cand_truth_frac_{s}"] = [H[s].get((f, d), np.nan) for f, d in zip(FW[m1].fid, FW[m1].days)]
    FW.to_parquet(f"{od}/family_window.parquet", index=False)
    FW.drop(columns=["days"]).to_csv(f"{od}/family_window.csv", index=False)
    W.log(LOG, f"[score:{a.run}] {len(R)} station-window rows, {len(FW)} family-window rows, {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
