"""Q8e runner: WA-NW variants, P-weighted draws for data + scrambled, position-permuted and time-shuffled on draws."""
import sys
import pandas as pd
import multiprocessing as mp
sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import q8_rtr_run as RR, q8_rtr_data as qd
RR.DET = "d2"
RR.TREM.update(qd.q8e_tremor_overrides())
VARS = sys.argv[1].split(",") if len(sys.argv) > 1 else ["V0", "Pa", "Pb", "Pc", "F1", "F2", "F3", "A1"]
SUF = sys.argv[2] if len(sys.argv) > 2 else ""
T = []
for v in VARS:
    cat = f"q8e_{v}"
    for var in ["raw|data", "raw|c1"]:
        E = RR.get(cat, var)[1]
        for ep in E.ep:
            for rep in range(5):
                T.append((cat, var, "pw", rep, ep, 0))
                if var == "raw|data":
                    T.append((cat, var, "pw_c2", rep, ep, 0)); T.append((cat, var, "pw_c1", rep, ep, 0))
print(len(T), "tasks", flush=True)
rows, mem = [], []
with mp.get_context("fork").Pool(6, initializer=RR._init, initargs=("d2",)) as p:
    for o, m in p.imap_unordered(RR.one, T, chunksize=2):
        rows += o; mem += m
pd.DataFrame(rows).to_parquet(f"{RR.OUT}/streaks_all_q8e{SUF}_d2.parquet", index=False)
pd.concat(mem).to_parquet(f"{RR.OUT}/members_data_q8e{SUF}_d2.parquet", index=False)
pd.DataFrame(T, columns=["cat", "var", "mode", "rep", "ep", "n_thin"]).to_csv(f"{RR.OUT}/runs_q8e{SUF}_d2.csv", index=False)
print("done", len(rows), flush=True)
