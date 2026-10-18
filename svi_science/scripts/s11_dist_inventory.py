#!/usr/bin/env python
"""s11_dist_inventory -- S11 distant-station lever, step 1 (distant/PLAN.md §2-§3): far pairs (80 < hypo <= 150 km),
candidate stations with raw files on disk, predicted S from 1-D travel times anchored on the near template offsets,
sample-date stacks -> S pick, pick-SNR, causality, scrambled-control pick-SNR, cc20, usable pairs, expected usable cells.
STOP marker (DIST_STOP) when expected usable cells < 100. Read-only on every input; writes distant/<region>/tables.
    S11_REGION=SVI   python s11_dist_inventory.py --workers 10
    S11_REGION=WA_NW python s11_dist_inventory.py --workers 10
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s11_common as C
import s11_dist_common as D
D.set_grid(D.INV_GRID)
import s10_windows as SW
from s10_measure import env, wmean

PER_YEAR = 4
SEARCH = 6.0
NDRAW, NSUB = 20, 20
MIN_CELLS = 100
_RTS = None
_AV = {}


class DistRaw(SW.Raw):
    """s10_windows.Raw on the distant roots; previous day read only when the span starts < 120 s after midnight
    (equivalent: a day file starts at 00:00, so a later span lies inside it)."""
    def __init__(self, netsta):
        self.ns = netsta; self.cache = {}

    def day(self, d):
        if d in self.cache:
            return self.cache[d]
        import obspy
        p = _AV.get((self.ns, d.strftime("%Y-%m-%d")))
        if p is None:
            p = D.raw_path(self.ns, d, _RTS)
        out = None
        if p is not None:
            try:
                st = obspy.read(p)
                bands = {}
                for tr in st:
                    c = tr.stats.channel
                    if c[-1] in SW.CORD:
                        bands.setdefault(c[:2], set()).add(SW.CORD[c[-1]])
                full = [b for b, cs in bands.items() if 0 in cs]
                if full:
                    band = sorted(full, key=lambda b: (-len(bands[b]), {"HH": 0, "EH": 1, "BH": 2}.get(b, 3)))[0]
                    st = st.select(channel=f"{band}?"); comps = {}
                    for c in set(tr.stats.channel for tr in st):
                        s2 = st.select(channel=c).copy()
                        try:
                            s2.merge(method=0, fill_value=0)
                        except Exception:
                            continue
                        tr = s2[0]
                        comps[SW.CORD[c[-1]]] = (tr.data.astype(np.float64), float(tr.stats.sampling_rate), float(tr.stats.starttime.timestamp))
                    out = comps
            except Exception:
                out = None
        self.cache[d] = out
        while len(self.cache) > 4:
            self.cache.pop(next(iter(self.cache)))
        return out

    def window(self, tS):
        a, b = tS - S.CUT_PRE, tS + S.CUT_POST
        res = [None, None, None]
        d0 = pd.Timestamp(a, unit="s").normalize(); d1 = pd.Timestamp(b, unit="s").normalize()
        days = {d0, d1}
        if a - d0.timestamp() < 120:
            days.add(d0 - pd.Timedelta(days=1))
        for c in range(3):
            for d in sorted(days):
                Dd = self.day(d)
                if not Dd or c not in Dd:
                    continue
                x, fs, t0 = Dd[c]
                i0 = int(np.floor((a - t0) * fs)); i1 = int(np.ceil((b - t0) * fs)) + 1
                if i0 >= 0 and i1 <= x.size:
                    res[c] = (x[i0:i1], fs, t0 + i0 / fs); break
        return res


SW.Raw = DistRaw


def task(arg):
    return SW.task(arg)


def pair_stats(Wd, P, Wc, sp, rng):
    tau = S.TAU
    out = dict(n_sample_data=int(np.isfinite(Wd[:, 0, 0]).sum()) if len(Wd) else 0,
               n_sample_c1=int(np.isfinite(Wc[:, 0, 0]).sum()) if len(Wc) else 0)
    def stats(G, prefix):
        e = np.zeros(S.NW); k = 0
        for c in range(3):
            if np.isfinite(G[c]).all():
                e += env(G[c]) ** 2; k += 1
        if k == 0:
            return None
        m = np.flatnonzero((tau >= -SEARCH) & (tau <= SEARCH))
        pk = float(tau[m[np.argmax(e[m])]])
        pre = (tau >= max(pk - sp - 9.0, S.T0W + 0.5)) & (tau <= pk - sp - 1.0)
        cw = (tau >= pk + 2.0) & (tau <= pk + 5.0)
        r = {f"{prefix}pick_shift": pk}
        if pre.sum() < int(2 * S.FS):
            r.update({f"{prefix}pick_snr": np.nan, f"{prefix}causal3c": np.nan, f"{prefix}causalZ": np.nan}); return r
        r[f"{prefix}pick_snr"] = float(np.sqrt(e[np.argmin(np.abs(tau - pk))]) / np.median(np.sqrt(e[pre])))
        rc = [np.mean(G[c][cw] ** 2) / np.mean(G[c][pre] ** 2) for c in range(3) if np.isfinite(G[c]).all()]
        r[f"{prefix}causal3c"] = float(np.sqrt(np.mean(rc)))
        r[f"{prefix}causalZ"] = float(np.sqrt(np.mean(G[0][cw] ** 2) / np.mean(G[0][pre] ** 2))) if np.isfinite(G[0]).all() else np.nan
        return r
    G = None
    if out["n_sample_data"] >= 3:
        G = np.stack([wmean(Wd[:, c], P)[0] for c in range(3)])
        r = stats(G, "")
        if r:
            out.update(r)
    if out["n_sample_c1"] >= 3:
        Gc = np.stack([wmean(Wc[:, c], np.ones(len(Wc)))[0] for c in range(3)])
        r = stats(Gc, "ctrl_")
        if r:
            out.update(r)
    # cc20: random 20-window stack vs the P-weighted stack of the rest, joint 3C cc 2-5 s past the pick
    if G is not None and "pick_shift" in out and out["n_sample_data"] >= 2 * NSUB:
        ok = np.flatnonzero(np.isfinite(Wd[:, 0, 0])); pk = out["pick_shift"]
        cw = (tau >= pk + 2.0) & (tau <= pk + 5.0); ccs = []
        for _ in range(NDRAW):
            sel = rng.choice(ok, NSUB, replace=False); rest = np.setdiff1d(ok, sel)
            v = []
            for c in range(3):
                a_ = Wd[sel, c][:, cw]; b_ = Wd[rest, c][:, cw]
                fa = np.isfinite(a_[:, 0]); fb = np.isfinite(b_[:, 0])
                if fa.sum() < 5 or fb.sum() < 5:
                    continue
                x = (P[sel][fa, None] * a_[fa]).sum(0) / P[sel][fa].sum(); y = (P[rest][fb, None] * b_[fb]).sum(0) / P[rest][fb].sum()
                x = x - x.mean(); y = y - y.mean(); v.append(float(x @ y / np.sqrt((x @ x) * (y @ y))))
            if v:
                ccs.append(np.mean(v))
        out["cc20"] = float(np.median(ccs)) if ccs else np.nan
    return out, G


def main():
    global _RTS, _AV
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--max-stations", type=int, default=0, help="test only")
    a = ap.parse_args()
    D.ensure(); log = C.Log(f"{D.DLOG}/inventory.log"); t00 = time.time()
    log(f"region {D.REG}; grid {S.T0W}..{S.T1W} s; velocity {os.path.basename(D.VEL)}")
    NT = D.near_tables(); F, X, O, A = NT["F"], NT["X"], NT["O"], NT["A"]
    Fi = F.set_index("fam"); XY = D.station_coords(); tt = D.TT()
    log(f"families {len(F)}, rows {len(X)} (data {(X.series == 'data').sum()}, c1 {(X.series == 'c1').sum()}), "
        f"near offsets {len(O)}, active family-days {len(A)} on {A.date.nunique()} dates")
    # ---------------- origin offsets from the near template S offsets
    O = O[O.netsta.isin(XY.index)].copy()
    epi = D.km(Fi.lat[O.fam].to_numpy(), Fi.lon[O.fam].to_numpy(), XY.lat[O.netsta].to_numpy(), XY.lon[O.netsta].to_numpy())
    O["epi"] = epi; O["tts"] = [tt.first(x, Fi.depth_km[f], "S") for x, f in zip(O.epi, O.fam)]
    O["oo"] = O.offset - O.tts
    g = O.groupby("fam").oo
    orig = g.median(); spread = g.agg(lambda x: float(np.median(np.abs(x - x.median()))))
    log(f"origin offset (template S - TT_S) per family: spread (MAD over near stations) median {spread.median():.2f} s, "
        f"90 % {spread.quantile(.9):.2f} s; near TT-vs-template residual vs epi: corr {np.corrcoef(O.epi, O.oo - O.fam.map(orig))[0, 1]:.2f}")
    # ---------------- candidate stations and far pairs
    dates = sorted(A.date.unique())
    _RTS = D.roots()
    have_dir = {ns for r in _RTS if os.path.isdir(r) for ns in os.listdir(r)}
    rows = []
    for ns, r in XY.iterrows():
        if ns.split(".")[0] in D.TEMP_NETS or ns not in have_dir:
            continue
        e = D.km(F.lat.to_numpy(), F.lon.to_numpy(), r.lat, r.lon); h = np.hypot(e, F.depth_km.to_numpy())
        for f, ee, hh in zip(F.fam, e, h):
            if D.HMIN < hh <= D.HMAX and (int(f), ns) not in NT["near_pairs"]:
                rows.append(dict(fam=int(f), netsta=ns, epi_km=float(ee), hypo_km=float(hh)))
    FP = pd.DataFrame(rows)
    log(f"geometric far pairs {len(FP)} at {FP.netsta.nunique()} permanent stations with a raw folder")
    av = []
    for ns in sorted(FP.netsta.unique()):
        for d in dates:
            p = D.raw_path(ns, d, _RTS)
            if p:
                av.append((ns, d, p))
    AV = pd.DataFrame(av, columns=["netsta", "date", "path"])
    nav = AV.groupby("netsta").size()
    keep = nav[nav >= D.MIN_DATES].index
    log(f"stations with files on >= {D.MIN_DATES} of {len(dates)} active dates: {len(keep)} -> {sorted(keep)}; "
        f"dropped {sorted(set(FP.netsta) - set(keep))[:80]}")
    FP = FP[FP.netsta.isin(keep)].reset_index(drop=True)
    if a.max_stations:
        keep = sorted(FP.groupby("netsta").size().sort_values(ascending=False).index[:a.max_stations]); FP = FP[FP.netsta.isin(keep)]
    AV = AV[AV.netsta.isin(keep)]
    AV.to_parquet(f"{D.DTB}/avail.parquet", index=False)
    _AV = {(r.netsta, r.date): r.path for r in AV.itertuples()}
    FP["tt_s"] = [tt.first(x, Fi.depth_km[f], "S") for x, f in zip(FP.epi_km, FP.fam)]
    FP["tt_s_direct"] = [tt.direct(x, Fi.depth_km[f], "S") for x, f in zip(FP.epi_km, FP.fam)]
    FP["tt_p"] = [tt.first(x, Fi.depth_km[f], "P") for x, f in zip(FP.epi_km, FP.fam)]
    FP["sp"] = FP.tt_s - FP.tt_p
    FP["orig_off"] = FP.fam.map(orig); FP["orig_spread"] = FP.fam.map(spread)
    FP = FP[np.isfinite(FP.orig_off)].reset_index(drop=True)
    FP["off_pred"] = FP.orig_off + FP.tt_s
    XYk = XY.loc[sorted(FP.netsta.unique())].reset_index()
    dep = D.sensor_depths(list(XYk.netsta), log)
    XYk["sensor_depth_m"] = XYk.netsta.map(dep); XYk["n_dates_avail"] = XYk.netsta.map(nav)
    XYk["n_far_fams"] = XYk.netsta.map(FP.groupby("netsta").size())
    XYk.to_csv(f"{D.DTB}/stations_far.csv", index=False)
    log(f"far pairs {len(FP)} at {FP.netsta.nunique()} stations; hypo median {FP.hypo_km.median():.0f} km; S-P "
        f"{FP.sp.min():.1f}-{FP.sp.max():.1f} s")
    # ---------------- sample dates
    Xd = X[X.series == "data"]
    cnt = Xd.groupby("date").size().reset_index(name="n"); cnt["yr"] = cnt.date.str[:4]
    SD = cnt.sort_values("n", ascending=False).groupby("yr").head(PER_YEAR).sort_values("date")
    SD.to_csv(f"{D.DTB}/sample_dates.csv", index=False)
    log(f"sample dates {len(SD)} ({PER_YEAR} per year, {SD.yr.nunique()} years), data rows on them {int(SD.n.sum())}")
    offk = FP.set_index(["fam", "netsta"]).off_pred
    Xs = X[X.date.isin(set(SD.date))]
    import multiprocessing as mp
    rng = np.random.default_rng(C.SEED + 41)
    res = []; stacks = {}
    for ns in sorted(FP.netsta.unique()):
        t0 = time.time()
        fams = set(FP.fam[FP.netsta == ns])
        R = Xs[Xs.fam.isin(fams)].copy()
        R = R[[(ns, d) in _AV for d in R.date]]
        if not len(R):
            log(f"{ns}: no sample rows with files"); continue
        R["tS"] = R.t.to_numpy() + np.array([offk[(f, ns)] for f in R.fam])
        R["det_id"] = -1
        tasks = [(ns, d, g.reset_index(drop=True)) for d, g in R.groupby("date")]
        W = {}; nerr = 0
        with mp.get_context("fork").Pool(a.workers, maxtasksperchild=50) as p:
            for (_, d, rid, Wx, src, rr, err) in p.imap_unordered(task, tasks, chunksize=1):
                if err:
                    nerr += 1; continue
                for i_, w_ in zip(rid, Wx):
                    W[int(i_)] = w_
        R = R[R.row_id.isin(W.keys())]
        for f, g in R.groupby("fam"):
            gd = g[g.series == "data"]; gc = g[g.series == "c1"]
            Wd = np.array([W[i] for i in gd.row_id]).reshape(-1, 3, S.NW); Wc = np.array([W[i] for i in gc.row_id]).reshape(-1, 3, S.NW)
            sp = float(FP.sp[(FP.fam == f) & (FP.netsta == ns)].iloc[0])
            st, G = pair_stats(Wd, gd.P.to_numpy(float), Wc, sp, rng)
            st.update(fam=int(f), netsta=ns); res.append(st)
            if G is not None:
                stacks[f"{f}|{ns}"] = G.astype(np.float32)
        log(f"{ns}: {len(tasks)} station-days, {len(R)} windows ({nerr} day errors), {R.fam.nunique()} families, {time.time() - t0:.0f}s")
        del W
    ST = pd.DataFrame(res)
    FP = FP.merge(ST, on=["fam", "netsta"], how="left")
    np.savez_compressed(f"{D.DTB}/sample_stacks.npz", keys=np.array(sorted(stacks)), G=np.array([stacks[k] for k in sorted(stacks)]),
                        tau=S.TAU)
    # ---------------- usable rule
    thr = float(np.nanpercentile(FP.ctrl_pick_snr, 99))
    FP["usable"] = (FP.causal3c > 1.5) & (FP.pick_snr > thr) & (FP.pick_shift.abs() < SEARCH - 0.5)
    log(f"control pick-SNR p99 = {thr:.2f} (n {FP.ctrl_pick_snr.notna().sum()}); data pick-SNR median {FP.pick_snr.median():.2f}; "
        f"causality>1.5 {int((FP.causal3c > 1.5).sum())}, SNR>thr {int((FP.pick_snr > thr).sum())}, usable {int(FP.usable.sum())} of {len(FP)}")
    # ---------------- expected usable cells
    ep_map = dict(zip(zip(A.fam, A.date), A.episode))
    dead = D.dead_days()
    pot = []
    for r in FP[FP.usable].itertuples():
        g = Xd[Xd.fam == r.fam]
        g = g[[(r.netsta, d) in _AV and (r.netsta, d) not in dead for d in g.date]]
        if not len(g):
            continue
        ep = D.epoch_of(D.epochs_table(r.netsta), g.date)
        e_ = np.array([ep_map.get((r.fam, d)) for d in g.date], dtype=object)
        q = pd.DataFrame(dict(epoch=ep, episode=e_)); q = q[(q.epoch != "") & q.episode.notna()]
        n = q.groupby(["epoch", "episode"]).size()
        n = n[n >= 20]
        rho20 = max(r.cc20, 0.01) if np.isfinite(r.cc20) else np.nan
        ccn = np.sqrt(1.0 / (1.0 + (1.0 / rho20 - 1.0) * 20.0 / n.to_numpy())) if np.isfinite(rho20) else np.full(len(n), np.nan)
        for (epc, epi_), nn, cc in zip(n.index, n.to_numpy(), ccn):
            pot.append(dict(fam=r.fam, netsta=r.netsta, epoch=epc, episode=epi_, n=int(nn), cc_pred=float(cc)))
    PC = pd.DataFrame(pot, columns=["fam", "netsta", "epoch", "episode", "n", "cc_pred"])
    PC["expected_usable"] = PC.cc_pred >= 0.35
    PC.to_csv(f"{D.DTB}/potential_cells.csv", index=False)
    FP["n_pot_cells"] = FP.set_index(["fam", "netsta"]).index.map(PC.groupby(["fam", "netsta"]).size()).fillna(0).astype(int)
    FP["n_exp_cells"] = FP.set_index(["fam", "netsta"]).index.map(PC[PC.expected_usable].groupby(["fam", "netsta"]).size()).fillna(0).astype(int)
    FP.to_csv(f"{D.DTB}/far_pairs.csv", index=False)
    nexp = int(PC.expected_usable.sum())
    bysta = FP.groupby("netsta").agg(hypo=("hypo_km", "median"), n_pairs=("fam", "size"), usable=("usable", "sum"),
                                     snr_med=("pick_snr", "median"), causal_med=("causal3c", "median"),
                                     shift_med=("pick_shift", "median"), cc20_med=("cc20", "median"), exp_cells=("n_exp_cells", "sum"))
    bysta.to_csv(f"{D.DTB}/far_stations_summary.csv")
    log("per station:\n" + bysta.round(2).to_string())
    summ = dict(region=D.REG, far_pairs=len(FP), stations=int(FP.netsta.nunique()), ctrl_snr_p99=thr, usable_pairs=int(FP.usable.sum()),
                potential_cells=len(PC), expected_usable_cells=nexp, usable_stations=int(FP[FP.usable].netsta.nunique()),
                min_cells_rule=MIN_CELLS, stop=bool(nexp < MIN_CELLS), written=C.pt_now(), secs=round(time.time() - t00))
    json.dump(summ, open(f"{D.DTB}/inventory_summary.json", "w"), indent=1)
    log(f"SUMMARY {summ}")
    if nexp < MIN_CELLS and not a.max_stations:
        open(f"{D.DOUT}/DIST_STOP", "w").write(f"{C.pt_now()} expected usable cells {nexp} < {MIN_CELLS}\n")
        log("STOP rule: expected usable cells < 100 -> DIST_STOP")
    if not a.max_stations:
        open(f"{D.DOUT}/INVENTORY_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
