#!/usr/bin/env python
"""bw_p6_apply -- P6: apply the confidence models to EVERY candidate (data and scrambled control) of every backward run
(rules: results/backward/p6/P6.md).  Output p6/P/<run>_P.parquet: run, date, row, series, fam, t, net_mad, P (= P_same),
P_same, P_any, P_same_foldsd, P_any_foldsd, model, net_set (+ P_same_frozenS15 for CN4/B1 rows, P_era_full / P_era_noamp for
era rows).  row = 0-based row of <date>.parquet (the p4b sidecar key).  Read-only on every feature file.

  python bw_p6_apply.py [run ...]       (default: all runs)
"""
from __future__ import annotations
import os, sys, glob, time, warnings
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
warnings.filterwarnings("ignore")

BW = "/home/jovyan/tremorferometry/svi_science/results/backward"
P6 = f"{BW}/p6"
FWD = "/home/jovyan/tremorferometry/lfe_extension/results/refine/conf/models.joblib"
BOUND = [("2005-07-01", "B1"), ("2005-11-02", "B2"), ("2007-11-13", "B3"), ("2009-08-27", "S15")]
KELB_GAIN_END = "2003-07-31"
_M = {}


def set_of_date(d):
    s = "CN4"
    for b, n in BOUND:
        if d >= b: s = n
    return s


def runs_all():
    R = []
    for p in sorted(glob.glob(f"{BW}/mf/feat/*")):
        R.append((os.path.basename(p), p))
    for p in sorted(glob.glob(f"{BW}/p4b/mf/feat/lad_*")):
        R.append((os.path.basename(p), p))
    for p in sorted(glob.glob(f"{BW}/p4/mf/feat/p4e0_*")):
        R.append((os.path.basename(p), p))
    return R


def models():
    if not _M:
        import joblib
        sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_extension/scripts")
        import svi_ref_fit as F
        bw = joblib.load(f"{P6}/models/bw_models.joblib"); fw = joblib.load(FWD)
        _M.update(bw=bw, fw=fw, WAVE=list(F.WAVE))
    return _M


def predict(x, key):
    """-> (P_same, P_any, sd_same, sd_any) for model key CN4/B1/B2/ERA/ERA_NOAMP/FWD."""
    M = models()
    out = []
    for lab in ("same", "any"):
        if key == "FWD":
            feats = M["WAVE"]; m = M["fw"]["models"][f"WAVE_{lab}"]; fm = M["fw"]["folds"].get(f"WAVE_{lab}", [])
        else:
            feats = M["bw"]["feats"][f"{key}_{lab}"]; m = M["bw"]["models"][f"{key}_{lab}"]; fm = M["bw"]["folds"][f"{key}_{lab}"]
        X = x[feats].to_numpy(float)
        p = m.predict_proba(X)[:, 1]
        sd = np.stack([f.predict_proba(X)[:, 1] for f in fm]).std(0, ddof=1) if len(fm) > 1 else np.full(len(x), np.nan)
        out += [p, sd]
    return out[0], out[2], out[1], out[3]


def one_day(args):
    run, path, ds = args
    p = f"{path}/{ds}.parquet"
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return None
    x = pd.read_parquet(p)
    if not len(x):
        return None
    M = models()
    for c in M["WAVE"]:
        if c not in x: x[c] = np.nan
    for c in ("has_LRIV", "has_SNB", "keep_xfam8"):
        if c in x: x[c] = x[c].astype(float)
    o = pd.DataFrame(dict(run=run, date=x.date.astype(str).to_numpy(), row=np.arange(len(x)), series=x.series.astype(str).to_numpy(),
                          fam=x.fam.to_numpy().astype(int), t=x.t.to_numpy(float), net_mad=x.net_mad.to_numpy(float)))
    era = run.startswith("p5e_")
    if era:
        a = predict(x, "ERA"); b = predict(x, "ERA_NOAMP")
        use = b if ds < KELB_GAIN_END else a
        o["P_era_full"] = a[0]; o["P_era_noamp"] = b[0]
        key = "ERA_NOAMP" if ds < KELB_GAIN_END else "ERA"; ns = "ERA"
    else:
        if run.startswith("lad_"):
            ns = run.split("_")[1]
        elif run.startswith("p4e0_"):
            ns = run.split("_")[1]; ns = "S15" if ns == "S15fix" else ns
        else:
            ns = set_of_date(ds)
        key = ns if ns in ("CN4", "B1", "B2") else "FWD"
        use = predict(x, key)
        if ns in ("CN4", "B1"):
            o["P_same_frozenS15"] = predict(x, "FWD")[0]
    o["P_same"], o["P_any"], o["P_same_foldsd"], o["P_any_foldsd"] = use
    o["P"] = o.P_same; o["model"] = key; o["net_set"] = ns
    return o


def main():
    os.makedirs(f"{P6}/P", exist_ok=True)
    want = sys.argv[1:]
    for run, path in runs_all():
        if want and run not in want: continue
        t0 = time.time()
        ds = sorted(os.path.basename(p)[:10] for p in glob.glob(f"{path}/*.json"))
        with ProcessPoolExecutor(8) as ex:
            res = [r for r in ex.map(one_day, [(run, path, d) for d in ds], chunksize=4) if r is not None]
        A = pd.concat(res, ignore_index=True) if res else pd.DataFrame()
        A.to_parquet(f"{P6}/P/{run}_P.parquet", index=False)
        print(f"[apply] {run}: {len(A)} rows, {len(ds)} days, models {A.model.value_counts().to_dict() if len(A) else {}}, "
              f"{time.time() - t0:.0f}s", flush=True)
    open(f"{P6}/APPLY_DONE", "w").write(pd.Timestamp.now("UTC").isoformat() + "\n")


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    main()
