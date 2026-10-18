"""Q8 RTR: run the frozen detector on all catalogues, data + scrambled controls (C1 time-shuffled, C2 position-permuted)
+ Bostock thinned to our density. 8 worker processes.

Output: results/q8_rtr/tables/streaks_all.parquet (one row per kept streak, all classes, S >= 2)
        results/q8_rtr/tables/members_data.parquet (streak uid -> event row index in the loader frame; data/thin runs)
        results/q8_rtr/tables/runs.csv
"""
import sys
import time
import numpy as np
import pandas as pd
from multiprocessing import Pool

sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import q8_rtr_data as qd
import q8_rtr_lib as L
import q8_rtr_backzone as BZ

OUT = f"{L.OUT}/tables"
N_HEAD_CTRL = 5
N_VAR_CTRL = 2
N_THIN = 10
_EV = {}
DET = "d1"
TREM = {}


def get(cat, var):
    key = (cat, var)
    if key not in _EV:
        if cat == "bostock":
            _EV[key] = qd.bostock()
        elif cat == "svi":
            _EV[key] = qd.ours_svi(var)
        elif cat == "nw":
            _EV[key] = qd.ours_nw(var)
        elif cat.startswith("ac_"):
            v, s = var.split("|")
            reg = cat[3:]
            if reg == "svi":
                _EV[key] = qd.allcand_svi(v, s)
            elif reg == "back":
                _EV[key] = qd.back_2010(v, s)
            else:
                _EV[key] = qd.allcand_wa(reg.upper(), v, s)
        elif cat.startswith("q8e_"):
            _EV[key] = qd.q8e_nw(cat[4:], var.split("|")[1])
        elif cat.startswith("pnsn_"):
            _EV[key] = qd.pnsn(cat.split("_")[1].upper())
    return _EV[key]


def one(task):
    cat, var, mode, rep, ep, n_thin = task
    ev, E = get(cat, var)
    r = E[E.ep == ep].iloc[0]
    g = ev[ev.ep == ep]
    rows = g.index.values
    rng = np.random.default_rng(abs(hash((cat, var, mode, rep, ep))) % (2 ** 32))
    t = ((g.time - r.w0).dt.total_seconds() / 3600).values
    x = g.x.values.astype(float); unit = g.unit.values
    T = (r.w1 - r.w0).total_seconds() / 3600
    if mode in ("pw", "pw_c2", "pw_c1"):
        keep = rng.random(len(t)) < np.nan_to_num(g.P.values.astype(float), nan=0.0)
        t, x, unit, rows = t[keep], x[keep], unit[keep], rows[keep]
        mode2 = {"pw": "data", "pw_c2": "c2", "pw_c1": "c1"}[mode]
    elif mode.startswith("thin"):
        if len(t) > n_thin:
            keep = np.sort(rng.choice(len(t), n_thin, replace=False))
            t, x, unit, rows = t[keep], x[keep], unit[keep], rows[keep]
        mode2 = "thin" if mode.startswith("thin@") else (mode.split("_")[1] if "_" in mode else "data")
    else:
        mode2 = mode
    day = np.floor(t / 24.0)
    if mode2 == "c1":
        t = L.control_c1(t, day, rng)
    elif mode2 == "c2":
        x = L.control_c2(x, unit, day, rng)
    o = np.argsort(t, kind="stable"); t, x, rows, unit = t[o], x[o], rows[o], unit[o]
    if len(t) < 20:
        return [], []
    t1 = time.time()
    if DET == "d2":
        if cat in ("q8e_F1", "q8e_F3"):
            on_ovr, md = TREM[ep]
            K = BZ.detect(t, x, T, unit=unit, main_dir=md, onset_override=on_ovr)
        else:
            K = BZ.detect(t, x, T, unit=unit, main_dir=BZ.main_direction(t, x))
    else:
        xq = np.percentile(g.x.values, [2, 98])          # range from the real episode, identical for controls
        K = L.detect(t, x, T, xq[0], xq[1], unit=unit)
    out, mem = [], []
    for i, k in enumerate(K):
        uid = f"{cat}|{var}|{mode}|{rep}|{ep}|{i}"
        d = {a: b for a, b in k.items() if a != "idx"}
        d.update(uid=uid, cat=cat, var=var, mode=mode, rep=rep, ep=ep, n_ep=len(t), sec=time.time() - t1)
        out.append(d)
        if mode2 == "data":
            mem.append(pd.DataFrame(dict(uid=uid, row=rows[k["idx"]])))
    return out, mem


def _init(det):
    global DET
    DET = det


def tasks():
    T = []
    svi_head = qd.ours_svi("head")[1]
    n_thin = int(np.median(svi_head.n))
    for cat, vars_ in [("bostock", ["all"]), ("svi", ["head", "p050", "raw075", "scr075"]), ("nw", ["head", "p050", "raw075"])]:
        for var in vars_:
            E = get(cat, var)[1]
            nc = N_HEAD_CTRL if var in ("head", "all") else (0 if var == "scr075" else N_VAR_CTRL)
            for ep in E.ep:
                T.append((cat, var, "data", 0, ep, 0))
                for rep in range(nc):
                    T.append((cat, var, "c1", rep, ep, 0)); T.append((cat, var, "c2", rep, ep, 0))
                if cat == "bostock":
                    for rep in range(N_THIN):
                        T.append((cat, var, "thin", rep, ep, n_thin))
                        T.append((cat, var, "thin_c1", rep, ep, n_thin))
                        T.append((cat, var, "thin_c2", rep, ep, n_thin))
    return T, n_thin


def main():
    T, n_thin = tasks()
    tag = sys.argv[1] if len(sys.argv) > 1 else "all"
    global DET
    DET = sys.argv[2] if len(sys.argv) > 2 else "d1"
    if tag == "pnsn":
        T = []
        for cat in ["pnsn_svi", "pnsn_nw", "pnsn_sw"]:
            for ep in get(cat, "all")[1].ep:
                T.append((cat, "all", "data", 0, ep, 0))
                for rep in range(N_HEAD_CTRL):
                    T.append((cat, "all", "c1", rep, ep, 0)); T.append((cat, "all", "c2", rep, ep, 0))
    if tag == "bostock_only":
        T = [z for z in T if z[0] == "bostock" and not z[2].startswith("thin")]
    suffix = ("" if tag == "all" else "_" + tag) + ("_d2" if DET == "d2" else "")
    print(len(T), "tasks; n_thin", n_thin, flush=True)
    # longest first (bostock dense)
    T.sort(key=lambda z: (z[0] != "bostock" or z[2].startswith("thin"),))
    t0 = time.time()
    allrows, allmem = [], []
    with Pool(8, initializer=_init, initargs=(DET,)) as p:
        for i, (o, m) in enumerate(p.imap_unordered(one, T, chunksize=1)):
            allrows += o; allmem += m
            if i % 50 == 0:
                print(i, len(allrows), f"{time.time() - t0:.0f}s", flush=True)
    S = pd.DataFrame(allrows)
    S.to_parquet(f"{OUT}/streaks_all{suffix}.parquet", index=False)
    pd.concat(allmem).to_parquet(f"{OUT}/members_data{suffix}.parquet", index=False)
    pd.DataFrame(T, columns=["cat", "var", "mode", "rep", "ep", "n_thin"]).to_csv(f"{OUT}/runs{suffix}.csv", index=False)
    print("done", len(S), f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
