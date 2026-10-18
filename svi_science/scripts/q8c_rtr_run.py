"""Q8c runner: D2 on all candidates (variants raw/once/P-weighted; data + scrambled-control series) + controls + thinned Bostock."""
import sys, time
import numpy as np, pandas as pd
import multiprocessing as mp
sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import q8_rtr_run as RR

OUT = RR.OUT


def main():
    RR.DET = "d2"
    T = []
    for cat in ["ac_svi", "ac_nw", "ac_sw", "ac_back"]:
        for v in ["raw", "once"]:
            for s in ["data", "c1"]:
                var = f"{v}|{s}"
                E = RR.get(cat, var)[1]
                for ep in E.ep:
                    T.append((cat, var, "data", 0, ep, 0))
                    if s == "data":
                        for rep in range(5):
                            T.append((cat, var, "c1", rep, ep, 0)); T.append((cat, var, "c2", rep, ep, 0))
                    if v == "raw" and cat in ("ac_svi", "ac_nw", "ac_sw") and (s == "data" or cat == "ac_svi"):
                        for rep in range(5):
                            T.append((cat, var, "pw", rep, ep, 0))
    # thinned Bostock at the new densities
    RR.get("bostock", "all")
    dens = {"svi_raw": int(RR.get("ac_svi", "raw|data")[1].n.median()), "svi_once": int(RR.get("ac_svi", "once|data")[1].n.median()),
            "sw_raw": int(RR.get("ac_sw", "raw|data")[1].n.median())}
    ev = RR.get("ac_svi", "raw|data")[0]
    dens["svi_pw"] = int(ev.groupby("ep").P.sum().median())
    for lab, n in dens.items():
        for ep in RR.get("bostock", "all")[1].ep:
            for rep in range(10):
                T.append(("bostock", "all", f"thin@{lab}", rep, ep, n))
    print(len(T), "tasks", dens, flush=True)
    t0 = time.time(); rows, mem = [], []
    ctx = mp.get_context("fork")
    with ctx.Pool(8, initializer=RR._init, initargs=("d2",)) as p:
        for i, (o, m) in enumerate(p.imap_unordered(RR.one, T, chunksize=2)):
            rows += o; mem += m
            if i % 200 == 0:
                print(i, len(rows), f"{time.time() - t0:.0f}s", flush=True)
    pd.DataFrame(rows).to_parquet(f"{OUT}/streaks_all_q8c_d2.parquet", index=False)
    pd.concat(mem).to_parquet(f"{OUT}/members_data_q8c_d2.parquet", index=False)
    pd.DataFrame(T, columns=["cat", "var", "mode", "rep", "ep", "n_thin"]).to_csv(f"{OUT}/runs_q8c_d2.csv", index=False)
    pd.Series(dens).to_csv(f"{OUT}/q8c_thin_densities.csv")
    print("done", len(rows), f"{time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
