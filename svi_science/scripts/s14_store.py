#!/usr/bin/env python
"""s14_store -- S14 step B: from the S11 day-stack stores rebuild each episode stack and its leave-episode-out reference,
(1) reproduce cc3c_25 of episode_raw (rule: median |dcc| <= 0.005), (2) CC of the direct-S window, (3) CC between
successive episode stacks, (4) per-day noise ratio and day CC vs the series grand stack (all windows + direct S, c2 too).
    [S11_REGION=WA_NW] python s14_store.py [--pilot] [--workers 4]
"""
import os, sys, argparse, time
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import pandas as pd
import s14_common as C
import s10_common as S
from s11_measure import Ref3
from s11_r3_inject import load_store_file

REGION = os.environ.get("S11_REGION", "SVI").upper()
STORE = {"SVI": f"{C.S11}/tables/store", "WA_NW": f"{C.S11}/wa_nw/tables/store"}[REGION]
ACT = {"SVI": [f"{C.RES}/s10_dvv/multi/tables/active_days.csv", f"{C.S11}/tables/ext/active_days.csv"],
       "WA_NW": [f"{C.S11}/wa_nw/tables/inv/active_days.csv"]}[REGION]
W = {"S": (-0.3, 1.5), "25": (2.0, 5.0), "2035": (2.0, 3.5), "3560": (3.5, 6.0)}
_G = {}


def rms(x, m):
    return float(np.sqrt(np.mean(x[m] ** 2)))


def series_job(k):
    st = _G["store"][k]; amap = _G["amap"]
    fam = int(k.split("|")[0]); ns = k.split("|")[1]; ep = k.split("|")[2]
    St, Ws, days, tS, sp = st["St"].astype(np.float64), st["Ws"], np.asarray(st["days"]).astype(str), st["tS"], st["sp"]
    okc = np.isfinite(st["G1"]).all(1)
    if not okc.any() or len(days) < 2:
        return None
    tr = S.TAU - tS
    mend = -((sp if np.isfinite(sp) else 6.0) + 1.0)
    pre = (S.TAU >= -13.5) & (tr <= mend)
    cw = (tr >= 2.0) & (tr <= 5.0)
    Gs = np.einsum("dc,dcn->cn", Ws, St); Wt = Ws.sum(0)
    G = np.where(Wt[:, None] > 0, Gs / np.where(Wt[:, None] > 0, Wt[:, None], 1), np.nan); G[~okc] = np.nan
    # ---- day level against the grand stack
    dayrows = []
    R = {w: Ref3(G, tS, win) for w, win in W.items()}
    meas = {w: R[w].measure(St, (0, 1, 2)) for w in W}
    for i, d in enumerate(days):
        r = dict(series=k, fam=fam, netsta=ns, epoch=ep, date=d, kind="data", ws=float(Ws[i].mean()),
                 r_day=rms(St[i, 0], pre) / max(rms(St[i, 0], cw), 1e-30), episode=amap.get((fam, d), None))
        for w in W:
            r[f"eps_{w}"], r[f"cc_{w}"] = float(meas[w][0][i]), float(meas[w][1][i])
        # Amendment 1: leave-day-out reference (removes the self-inclusion bias of heavy days)
        wl = Wt - Ws[i]
        if np.all(wl[okc] > 0):
            Gl = np.stack([(Gs[c] - Ws[i, c] * St[i, c]) / wl[c] if wl[c] > 0 else np.full(S.NW, np.nan) for c in range(3)])
            Gl[~okc] = np.nan
            for w, win in W.items():
                e, c = Ref3(Gl, tS, win).measure(St[i][None], (0, 1, 2))
                r[f"lo_eps_{w}"], r[f"lo_cc_{w}"] = float(e[0]), float(c[0])
        dayrows.append(r)
    if len(st["Sc"]):
        Sc = st["Sc"].astype(np.float64)
        mc = {w: R[w].measure(Sc, (0, 1, 2)) for w in W}
        for j, di in enumerate(st["c2_day"]):
            d = days[di]
            r = dict(series=k, fam=fam, netsta=ns, epoch=ep, date=d, kind="ctrl", ws=float(st["Cn"][j].mean()),
                     r_day=rms(Sc[j, 0], pre) / max(rms(Sc[j, 0], cw), 1e-30), episode=amap.get((fam, d), None))
            for w in W:
                r[f"eps_{w}"], r[f"cc_{w}"] = float(mc[w][0][j]), float(mc[w][1][j])
                r[f"lo_eps_{w}"], r[f"lo_cc_{w}"] = r[f"eps_{w}"], r[f"cc_{w}"]     # control days are not in the reference
            dayrows.append(r)
    # ---- episode level: rebuild stacks, leave-episode-out reference
    epi = np.array([amap.get((fam, d), None) for d in days], object)
    eprows, stacks = [], []
    for e_ in pd.unique(epi):
        if e_ is None:
            continue
        m = epi == e_
        ws = Ws[m].sum(0)
        if np.any((Wt - ws)[okc] <= 0):
            continue
        Sd = np.where(ws[:, None] > 0, np.einsum("dc,dcn->cn", Ws[m], St[m]) / np.where(ws[:, None] > 0, ws[:, None], 1), np.nan)
        ref = np.stack([(Gs[c] - ws[c] * np.nan_to_num(Sd[c])) / (Wt[c] - ws[c]) if Wt[c] - ws[c] > 0 else np.full(S.NW, np.nan)
                        for c in range(3)])
        ref[~okc] = np.nan
        r = dict(fam=fam, netsta=ns, epoch=ep, episode=e_, ndays_store=int(m.sum()), wsum=float(ws.mean()),
                 date_first=str(sorted(days[m])[0]))
        for w, win in W.items():
            e, c = Ref3(ref, tS, win).measure(Sd[None], (0, 1, 2))
            r[f"rec_eps_{w}"], r[f"rec_cc_{w}"] = float(e[0]), float(c[0])
        eprows.append(r); stacks.append(Sd)
    # successive episodes (date order): stack k is the reference for stack k+1
    order = np.argsort([r["date_first"] for r in eprows])
    for a, b in zip(order[:-1], order[1:]):
        ra, rb = eprows[a], eprows[b]
        e, c = Ref3(stacks[a], tS, W["25"]).measure(stacks[b][None], (0, 1, 2))
        rb.update(prev_episode=ra["episode"], prev_wsum=ra["wsum"], next_eps_25=float(e[0]), next_cc_25=float(c[0]),
                  gap_days=(pd.Timestamp(rb["date_first"]) - pd.Timestamp(ra["date_first"])).days)
    return dayrows, eprows


def run_station(path, keys, amap, workers):
    import multiprocessing as mp
    store = {k: v for k, v in load_store_file(path).items() if keys is None or k in keys}
    _G.clear(); _G.update(store=store, amap=amap)
    D, E = [], []
    with mp.get_context("fork").Pool(workers) as p:
        for o in p.imap_unordered(series_job, sorted(store), chunksize=2):
            if o is not None:
                D += o[0]; E += o[1]
    return D, E


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--pilot", action="store_true"); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    tag = f"{REGION}{'_pilot' if a.pilot else ''}"
    t0 = time.time()
    A = pd.concat([pd.read_csv(f) for f in ACT if os.path.exists(f)], ignore_index=True).drop_duplicates(["fam", "date"])
    amap = {(int(f), str(d)): e for f, d, e in zip(A.fam, A.date, A.episode)}
    files = sorted(f for f in os.listdir(STORE) if f.endswith(".npz"))
    if a.pilot:
        files = [f"{s}.npz" for s in C.PILOT_STAS]
        keys = {f"{fa}|{s}|{e}" for fa in C.PILOT_FAMS for s in C.PILOT_STAS for e in ("E1", "E2", "E3", "E4", "E5", "E6")}
    else:
        keys = None
    D, E = [], []
    for f in files:
        d, e = run_station(f"{STORE}/{f}", keys, amap, a.workers); D += d; E += e
        print(f"{C.pt_now()} {f}: {len(e)} episode rows, {len(d)} day rows ({time.time() - t0:.0f} s)", flush=True)
    D, E = pd.DataFrame(D), pd.DataFrame(E)
    D.to_parquet(f"{C.TB}/store_day_{tag}.parquet", index=False)
    E.to_parquet(f"{C.TB}/store_episode_{tag}.parquet", index=False)
    # reproduction check against episode_raw
    R = pd.read_parquet(C.EP_FILE[REGION])[["fam", "netsta", "epoch", "episode", "cc3c_25", "eps3c_25"]]
    M = E.merge(R, on=["fam", "netsta", "epoch", "episode"], how="inner")
    dcc = (M.rec_cc_25 - M.cc3c_25).abs(); dep = (M.rec_eps_25 - M.eps3c_25).abs()
    msg = (f"{tag}: matched {len(M)} of {len(E)} rebuilt / {len(R)} table cells; median |dcc| {dcc.median():.5f}, p95 {dcc.quantile(.95):.4f}; "
           f"median |deps| {dep.median():.2e}; rule (median |dcc| <= 0.005): {'PASS' if dcc.median() <= 0.005 else 'FAIL'}")
    print(msg)
    open(f"{C.TB}/store_check_{tag}.txt", "w").write(msg + "\n")
    open(f"{C.OUT}/STORE_{tag}_DONE", "w").write(C.pt_now() + "\n")


if __name__ == "__main__":
    main()
