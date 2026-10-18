"""Q1 scrambled control (SVI c1 series = stations shifted to different days), P >= 0.5 counted once: same budget tests.
Pre-registered in PLAN.md §5 (run because every big episode has >= 30 control events at P >= 0.5)."""
import json
import numpy as np
import q1_tests as q

CL, PA, EP, TR, DT, S3, KIN, GEO = q.load()
D = dict(CL=CL, PA=PA, EP=EP, TR=TR, DT=DT, S3=S3, KIN=KIN, GEO=GEO)
out = {}
for name, var in [("SVI_scrambled_c1_once050", "c1_once050"), ("SVI_data_once050_same", "once050")]:
    des = q.make(D, "SVI", var, "group", "big", None, False)
    res, df0 = q.run_budget(des, nperm=2000)
    out[name] = {**res, "n_censored": int(df0.cens.sum()), "n_rows": int(len(df0)), "N_targets": des.N[des.tslots].tolist()}
    cov = q.covariates(des, D, "SVI")
    kt, K = q.kin_tests(des, cov, nperm=2000)
    out[name + "|kin"] = kt
    U = (des.pL ** 2).sum(0) / (len(des.patches) - 1)
    out[name + "|U"] = U.tolist()
    print(name, json.dumps(out[name], default=float)[:600], kt, flush=True)
json.dump(out, open(f"{q.OUT}/scrambled_control.json", "w"), indent=1, default=float)
