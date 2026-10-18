#!/usr/bin/env python
"""s11_svdpath -- S11 step 4, data path (ii) = the user's recipe (SVD-Wiener 3C day stacks, 2-5 s stretch anchored at the
picked S, days CC > 0.8, CC-weighted episode mean), from the day-stack store; scrambled-control cells through the same
path; and the per-series waveform-level STEP-INJECTION RECOVERY of both paths (headline rule, PLAN §1):
step at the episode boundary nearest the middle of each series, eps steps +0.1 % and +0.3 %, recovery =
[(post - pre)_injected - (post - pre)_natural] / a over episode values (equal weights). eps = -dv/v.

    python s11_svdpath.py [--workers 12] [--tag ""]
Outputs: tables/episode_svd.parquet, tables/episode_svd_ctrl.parquet, tables/recovery.csv, tables/headline.json, SVD_DONE.
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_paths as P

AMPS = (0.001, 0.003)
_ST = {}


def _series(key):
    st = _ST[key]; unit = _ST["__unit__"]
    fam, ns, ep = key.split("|")
    tag = dict(fam=int(fam), netsta=ns, epoch=ep)
    Dd = P.svd_days(st, unit)
    Ec = P.svd_cells(Dd)
    Dc = P.svd_days(st, unit, control=True)
    Ecc = P.svd_cells(Dc, control=True)
    rec = []
    # ---- step-injection recovery (both paths, headline window, 3C)
    days = st["days"]; ep_day = np.array([unit.get(d, None) for d in days], dtype=object)
    eps_order = [e for e in pd.unique(ep_day) if e is not None]
    if len(eps_order) >= 4:
        mid_day = str(days[len(days) // 2])
        bounds = [i for i in range(1, len(eps_order))]
        first_day = {e: str(days[ep_day == e][0]) for e in eps_order}
        b = min(bounds, key=lambda i: abs(pd.Timestamp(first_day[eps_order[i]]) - pd.Timestamp(mid_day)))
        post_eps = set(eps_order[b:])
        post = np.array([e in post_eps for e in ep_day])
        nat_raw = P.raw_cells(st, unit, wins=[C.HEAD]).set_index("episode")[f"eps3c_{C.HEAD}"]
        nat_svd = Ec.set_index("episode")[f"eps3c_{C.HEAD}"] if len(Ec) else pd.Series(dtype=float)
        for a in AMPS:
            amp = np.where(post, a, 0.0)
            inj_raw = P.raw_cells(st, unit, amp=amp, wins=[C.HEAD]).set_index("episode")[f"eps3c_{C.HEAD}"]
            Ei = P.svd_cells(P.svd_days(st, unit, amp=amp, wins=[C.HEAD]), wins=[C.HEAD])
            inj_svd = Ei.set_index("episode")[f"eps3c_{C.HEAD}"] if len(Ei) else pd.Series(dtype=float)
            for path, nat, inj in (("raw", nat_raw, inj_raw), ("svd", nat_svd, inj_svd)):
                j = pd.concat([nat.rename("n"), inj.rename("i")], axis=1).dropna()
                isp = j.index.isin(post_eps)
                if isp.sum() >= 1 and (~isp).sum() >= 1:
                    r = ((j.i[isp].mean() - j.i[~isp].mean()) - (j.n[isp].mean() - j.n[~isp].mean())) / a
                    lk = (j.i[~isp] - j.n[~isp]).mean() / a
                else:
                    r = lk = np.nan
                rec.append(dict(**tag, path=path, a=a, recovery=r, pre_shift=lk, n_pre=int((~isp).sum()), n_post=int(isp.sum()),
                                boundary=first_day[eps_order[b]]))
    for D in (Ec, Ecc):
        for k, v in tag.items():
            if len(D):
                D[k] = v
    return Ec, Ecc, pd.DataFrame(rec)


def _safe(key):
    try:
        return _series(key)
    except Exception as e:
        import traceback
        print(f"ERROR {key}: {e} {traceback.format_exc()[-400:]}", flush=True)
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--tag", default="")
    a = ap.parse_args()
    log = C.Log(f"{C.LOGD}/svdpath{a.tag}.log")
    unit = P.unit_map(); t0 = time.time()
    E, EC, R = [], [], []
    import multiprocessing as mp
    for ns in C.STATIONS:
        p = f"{C.TB}/store{a.tag}/{ns}.npz"
        if not os.path.exists(p):
            continue
        _ST.clear(); _ST.update(P.load_store(ns, a.tag)); _ST["__unit__"] = unit
        keys = [k for k in _ST if k != "__unit__"]
        with mp.get_context("fork").Pool(a.workers) as pool:
            for k, (Ec, Ecc, rec) in zip(keys, pool.imap(_safe, keys)):
                E.append(Ec); EC.append(Ecc); R.append(rec)
        log(f"{ns}: {len(keys)} series {time.time() - t0:.0f}s")
    E = pd.concat([e for e in E if len(e)], ignore_index=True); EC = pd.concat([e for e in EC if len(e)], ignore_index=True)
    R = pd.concat([r for r in R if len(r)], ignore_index=True)
    E.to_parquet(f"{C.TB}/episode_svd{a.tag}.parquet", index=False); EC.to_parquet(f"{C.TB}/episode_svd_ctrl{a.tag}.parquet", index=False)
    R.to_csv(f"{C.TB}/recovery{a.tag}.csv", index=False)
    # ---- headline rule (PLAN §1)
    s = R.groupby(["path", "a"]).recovery.agg(["median", "count", lambda x: float(np.mean((x >= 0.9) & (x <= 1.1)))])
    s.columns = ["median", "n_series", "frac_in_0.9_1.1"]
    log("step recovery (3C 2-5 s, per series):\n" + s.round(3).to_string())
    med = {aa: float(s.loc[("svd", aa), "median"]) for aa in AMPS}
    svd_ok = all(0.9 <= v <= 1.1 for v in med.values())
    head = dict(headline="svd" if svd_ok else "raw", svd_median_recovery={str(k): v for k, v in med.items()},
                raw_median_recovery={str(aa): float(s.loc[("raw", aa), "median"]) for aa in AMPS},
                rule="SVD headline iff median SVD step recovery in 0.9-1.1 at both 0.1 % and 0.3 %; else raw",
                written=C.pt_now())
    json.dump(head, open(f"{C.TB}/headline{a.tag}.json", "w"), indent=1)
    log(f"HEADLINE PATH: {head['headline']}  {head}")
    open(f"{C.OUT}/SVD{a.tag}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
