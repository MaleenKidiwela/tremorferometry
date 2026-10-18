#!/usr/bin/env python
"""s11_inject -- S11 step 5, injection suite D7(b), differential, waveform level, full chain, BOTH data paths.
Each day stack of a series is stretched about tS by eps_day = -(injected dv/v seen by that series) BEFORE any episode
stacking / SVD; the path is re-run (raw: episode stacks + leave-episode-out reference; svd: SVD-Wiener + grand stack +
days CC > 0.8 + CC-weighted means); the headline cells (same cells, same inflated sigma, nominal kernels) are inverted;
result = injected minus natural.
  A  deep step +0.3 % from 2015-01-01: dv/v = 0.3 % x a_d            -> D recovery 0.7-1.3, leak into each S < 30 %
  B  common shallow annual 0.1 % (peak day 233) at all sites: x a_s  -> D annual amplitude / 0.1 % < 25 % (> 50 % -> "deep +
     regional shallow")
  C  one-site shallow step +0.3 % from 2015-01-01 (each site in turn) -> |D step| / 0.3 % < 20 %
  D  as B with a_s x 1.5 or / 1.5 at random per (family, station)     -> D annual amplitude / 0.1 % < 30 %
eps = -dv/v.   python s11_inject.py [--workers 12]   -> tables/inject_*.csv, INJECT_DONE
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time, json, zlib
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_paths as P
import s11_invert as I

STEP_DATE = "2015-01-01"
PEAK = 233.0
_G = {}


def specs():
    out = {"nat": None, "A": ("deep_step", None), "B": ("shallow_annual", None), "D": ("shallow_annual_kerr", None)}
    for s in C.SITES:
        out[f"C_{s}"] = ("site_step", s)
    return out


def amp_for(kind, site_arg, days, a_s, a_d, site, fam, ns):
    t = pd.to_datetime(pd.Series(days))
    post = (t >= STEP_DATE).to_numpy()
    ann = np.cos(2 * np.pi * (t.dt.dayofyear.to_numpy() - PEAK) / 365.25)
    if kind == "deep_step":
        dvv = 0.003 * a_d * post
    elif kind == "shallow_annual":
        dvv = 0.001 * a_s * ann
    elif kind == "shallow_annual_kerr":
        f = 1.5 if (zlib.crc32(f"{fam}|{ns}".encode()) + C.SEED) % 2 == 0 else 1 / 1.5
        dvv = 0.001 * a_s * f * ann
    elif kind == "site_step":
        dvv = 0.003 * a_s * post * (site == site_arg)
    return -dvv                                                      # eps = -dv/v


def _series(key):
    st = _G["store"][key]; unit = _G["unit"]; kc = _G["kc"]
    fam, ns, ep = key.split("|")
    if key not in kc.index:
        return key, {}
    a_s, a_d = float(kc.loc[key, "a_s"]), float(kc.loc[key, "a_d"])
    out = {}
    for run, sp in specs().items():
        amp = None if sp is None else amp_for(sp[0], sp[1], st["days"], a_s, a_d, C.SITE[ns], fam, ns)
        if amp is not None and not np.any(amp):
            out[run] = None; continue                                 # identical to natural
        r = P.raw_cells(st, unit, amp=amp, wins=[C.HEAD])
        s = P.svd_cells(P.svd_days(st, unit, amp=amp, wins=[C.HEAD]), wins=[C.HEAD])
        out[run] = (r[["episode", f"eps3c_{C.HEAD}"]] if len(r) else r, s[["episode", f"eps3c_{C.HEAD}"]] if len(s) else s)
    return key, out


def _safe(key):
    try:
        return _series(key)
    except Exception as e:
        import traceback
        print(f"ERROR {key}: {e} {traceback.format_exc()[-400:]}", flush=True)
        return key, {}


def annual_amp(par):
    D = par[par.kind == "D"]
    return C.annual_fit(D.mid, D.value.to_numpy(), D.sig.to_numpy())


def step(par, kind="D", site=""):
    q = par[(par.kind == kind) & ((par.site == site) if kind == "S" else True)]
    post = pd.to_datetime(q.mid) >= STEP_DATE
    if post.sum() == 0 or (~post).sum() == 0:
        return np.nan
    w = 1 / q.sig ** 2
    return float(np.average(q.value[post], weights=w[post]) - np.average(q.value[~post], weights=w[~post]))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--tag", default="")
    ap.add_argument("--ccmin", type=float, default=C.CC_HEAD); ap.add_argument("--label", default="")
    ap.add_argument("--paths", default="raw,svd")
    a = ap.parse_args()
    log = C.Log(f"{C.LOGD}/inject{a.tag}.log"); t0 = time.time()
    unit = P.unit_map()
    K = I.kernels("ell40"); K = K[K.win == C.HEAD].set_index("series")
    vals = {}
    import multiprocessing as mp
    for ns in C.STATIONS:
        p = f"{C.TB}/store{a.tag}/{ns}.npz"
        if not os.path.exists(p):
            continue
        _G.clear(); _G.update(store=P.load_store(ns, a.tag), unit=unit, kc=K)
        keys = list(_G["store"])
        with mp.get_context("fork").Pool(a.workers) as pool:
            for key, out in pool.imap_unordered(_safe, keys):
                vals[key] = out
        log(f"{ns}: {len(keys)} series {time.time() - t0:.0f}s")
    rows, deltas = [], []
    for path in a.paths.split(","):
        cells = I.prepare(path, C.HEAD, a.tag, ccmin=a.ccmin)
        ckey = cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch
        pi = 0 if path == "raw" else 1

        def dvec(run):
            d = np.full(len(cells), np.nan)
            for i, (k, e) in enumerate(zip(ckey, cells.episode)):
                o = vals.get(k, {})
                v = o.get(run, None)
                if v is None:
                    v = o.get("nat", None)
                if v is None or not len(v[pi]):
                    continue
                q = v[pi]; m = q.episode == e
                if m.any():
                    d[i] = -100 * float(q[f"eps3c_{C.HEAD}"][m].iloc[0])
            return d
        dn = dvec("nat")
        log(f"{path}: natural-from-store vs measured cells: max |diff| {np.nanmax(np.abs(dn - cells.d)):.2e} %, "
            f"missing {int(np.isnan(dn).sum())}/{len(dn)}")
        for run, sp in specs().items():
            if sp is None:
                continue
            di = dvec(run)
            ok = np.isfinite(di) & np.isfinite(dn)
            cc = cells[ok].reset_index(drop=True)
            pn, _ = I.invert(cc, d_override=dn[ok], extras=False)
            pj, _ = I.invert(cc, d_override=di[ok], extras=False)
            dp = pn[["kind", "site", "episode", "mid", "sig"]].copy(); dp["value"] = pj.value - pn.value
            dp["run"] = run; dp["path"] = path; deltas.append(dp)
            r = dict(path=path, run=run, n_cells=int(ok.sum()), n_lost=int((~ok).sum()))
            if run == "A":
                r["D_recovery"] = step(dp, "D") / 0.3
                for s in C.SITES:
                    r[f"S_{s}_leak"] = step(dp, "S", s) / 0.3
                r["max_S_leak"] = float(np.nanmax(np.abs([r[f"S_{s}_leak"] for s in C.SITES])))
                r["pass"] = bool(0.7 <= r["D_recovery"] <= 1.3 and r["max_S_leak"] < 0.30)
            elif run in ("B", "D"):
                f = annual_amp(dp); r["D_annual_leak"] = f["amp"] / 0.1; r["D_annual_leak_sig"] = f["amp_sig"] / 0.1
                for s in C.SITES:
                    q = dp[(dp.kind == "S") & (dp.site == s)]
                    fs = C.annual_fit(q.mid, q.value.to_numpy(), q.sig.to_numpy()); r[f"S_{s}_annual_rec"] = fs["amp"] / 0.1
                r["pass"] = bool(r["D_annual_leak"] < (0.25 if run == "B" else 0.30))
                if run == "B":
                    r["headline_quantity"] = "deep + regional shallow" if r["D_annual_leak"] > 0.5 else "deep"
            else:
                site = run.split("_", 1)[1]
                r["D_step_leak"] = step(dp, "D") / 0.3; r[f"S_{site}_recovery"] = step(dp, "S", site) / 0.3
                r["pass"] = bool(abs(r["D_step_leak"]) < 0.20)
            rows.append(r)
            log(f"{path} {run}: " + ", ".join(f"{k} {v:.3f}" if isinstance(v, float) else f"{k} {v}" for k, v in r.items()
                                              if k not in ("path", "run")))
    pd.DataFrame(rows).to_csv(f"{C.TB}/inject_summary{a.tag}{a.label}.csv", index=False)
    pd.concat(deltas, ignore_index=True).to_csv(f"{C.TB}/inject_deltas{a.tag}{a.label}.csv", index=False)
    log(f"inject done {time.time() - t0:.0f}s")
    if not a.label:
        open(f"{C.OUT}/INJECT{a.tag}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
