import os, numpy as np, pandas as pd
ROOT="/home/jovyan/tremorferometry"
WS=f"{ROOT}/lfe_south"
RES=f"{WS}/results"; TR=f"{WS}/traces"; CACHE=f"{WS}/cache20hz"; LOGS=f"{WS}/logs"
FS=20.0
PROV={'NC':'NCEDC','BK':'NCEDC','PB':'IRIS','UW':'IRIS','CC':'IRIS','US':'IRIS','UO':'IRIS'}
def day_path(root,net,sta,d):
    d=pd.Timestamp(d)
    return os.path.join(root,f"{net}.{sta}",str(d.year),f"{d.dayofyear:03d}.mseed")
def npz_path(net,sta,d):
    d=pd.Timestamp(d)
    return os.path.join(CACHE,f"{net}.{sta}",str(d.year),f"{d.dayofyear:03d}.npz")
def log(m):
    import time,sys
    print(f"[{time.strftime('%H:%M:%S')}] {m}",flush=True)
