#!/usr/bin/env python
"""Quiet-day waveform test, step 1: cut the station windows of NW core-network candidates exactly where the production
matched filter correlated them (scripts/wa_prod_mf.py run_day/features):

  data, random   station position j = i + round(offset_s * 20)                (template moveout)
  rev            j = i + d, compared with the time-REVERSED template (the control) and with the forward template
  c1 (scrambled) j = ((i - shift_site) % NV) + d, shift_site = run_<date>.json shifts['core|c1|<site>']
                 (= wa_prod_common.site_shift(date, site, NV, 1), fixed per date and site)

Stations = the day's core members (run json, after the drop-dead rule) with a template for (family, station, tkey) and
>= 1 channel live in template and cache; site weights = production site merge B (coef).  Waveforms = the 20 Hz 1-8 Hz
cache (wa_prod_common.find_cache, first hit wins = the files the matched filter read).

Per window: station cc with the forward / reversed template (production formula: mean over channels of
dot(tpl_c, x_c) / |x_c|), the best forward cc within +-15 s, window energy / day level, envelope peak / median,
kurtosis; the waveform [j-20 s, j+32 s] normalised per channel by the day's robust std (float16) for the stacks.
Read-only on production outputs; writes only under results/wa_prod/quietwave/work/."""
import os, sys, json, time
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_extension_wa/scripts")
import wa_prod_common as P
W = P.W
import wa_prod_mf as MF

OD = f"{P.OUT}/quietwave"
WK = f"{OD}/work"
n = W.NT
NV = W.NDAY - n + 1
FS = W.FS
PRE = int(20 * FS); LEN = int(52 * FS)          # stack window [j-20 s, j+32 s]
LAG = int(15 * FS)
TPL = None


def init():
    global TPL
    TPL = MF.load_library(rich=False)


def coefs(ss, f, site, tk):
    wm = {s: float(ch.sum()) * TPL[(f, s, tk[s])]["w"] for s, ch in ss}
    groups = {}
    for s, _ in ss:
        groups.setdefault(site[s], []).append(s)
    Ws = {q: float(np.sqrt(sum(wm[s] ** 2 for s in mem))) for q, mem in groups.items()}
    Wt = sum(Ws.values())
    c = {}
    for q, mem in groups.items():
        wsum = sum(wm[s] for s in mem)
        for s in mem:
            c[s] = Ws[q] * wm[s] / wsum / Wt
    return c


def ccwin(x, t, ch):
    """production station cc: mean over live channels of dot(t_c, x_c)/|x_c| (template L2-normalised per channel)."""
    num = (x * t).sum(-1); den = np.sqrt((x * x).sum(-1))
    with np.errstate(invalid="ignore", divide="ignore"):
        v = num / den
    v = v[ch]
    return float(np.mean(v)) if np.all(np.isfinite(v)) else np.nan


def do_day(args):
    ds, R = args                                   # R: rows of that date (row_id, fid, series, i, group)
    out_p = f"{WK}/{ds}.npz"
    if os.path.exists(out_p):
        return ds, "exists"
    meta = json.load(open(f"{P.OUT}/mf/prod/run_{ds}.json"))
    mem = meta["members"].get("core", [])
    site = meta["sites"]["core"]; tk = meta["tkeys"]
    sh = {q: int(v) for k, v in meta["shifts"].items() for net, ser, q in [k.split("|")] if net == "core" and ser == "c1"}
    X = {}
    for s in mem:
        p = P.find_cache(s, ds)
        if p is None:
            continue
        z = np.load(p)
        x = z["x"].astype(np.float32); live = z["live"].astype(bool)
        lev = np.array([1.4826 * np.median(np.abs(x[c][x[c] != 0])) if np.any(x[c] != 0) else np.nan for c in range(3)])
        X[s] = (x, live, lev)
    recs = []; wav = []; rowinfo = []
    for r in R.itertuples():
        f = int(r.fid)
        ss = []
        for s in mem:
            if s not in X or (f, s, tk.get(s)) not in TPL:
                continue
            ch = TPL[(f, s, tk[s])]["live"] & X[s][1]
            if ch.any():
                ss.append((s, ch))
        if len({site[s] for s, _ in ss}) < 2:
            rowinfo.append((r.row_id, 0, np.nan)); continue
        cf = coefs(ss, f, site, tk)
        net_re = 0.0
        for s, ch in ss:
            t = TPL[(f, s, tk[s])]
            d = int(round(t["off"] * FS))
            i = int(r.i)
            if r.series == "c1":
                j = ((i - sh[site[s]]) % NV) + d
            else:
                j = i + d
            x, live, lev = X[s]
            tpl = t["tpl"]
            if j < 0 or j + n > W.NDAY:
                cf_ = np.nan; cr_ = np.nan
            else:
                seg = x[:, j:j + n].astype(np.float64)
                cf_ = ccwin(seg, tpl, ch); cr_ = ccwin(seg, tpl[:, ::-1], ch)
            if np.isfinite(cf_) and r.series != "rev":
                net_re += cf[s] * cf_
            elif np.isfinite(cr_) and r.series == "rev":
                net_re += cf[s] * cr_
            # lag scan (forward template), energy, envelope
            a0, b0 = j - LAG, j + LAG + n
            best, blag = np.nan, np.nan
            if a0 >= 0 and b0 <= W.NDAY:
                seg2 = x[:, a0:b0].astype(np.float64)
                vals = []
                for c in np.flatnonzero(ch):
                    num = np.correlate(seg2[c], tpl[c], "valid")
                    e = np.cumsum(np.r_[0, seg2[c] ** 2]); den = np.sqrt(e[n:] - e[:-n])
                    with np.errstate(invalid="ignore", divide="ignore"):
                        vals.append(num / den)
                v = np.mean(vals, 0)
                if np.all(np.isfinite(v)):
                    k = int(np.argmax(v)); best = float(v[k]); blag = (k - LAG) / FS
            w0, w1 = j - PRE, j - PRE + LEN
            lv = np.where(ch, lev, 1.0)
            if w0 >= 0 and w1 <= W.NDAY and np.all(np.isfinite(lv)) and np.all(lv > 0):
                wv = (x[:, w0:w1] / lv[:, None]).astype(np.float32)
                wv[~ch] = 0
                sig = wv[ch][:, PRE:PRE + n]; noi = wv[ch][:, :PRE - int(8 * FS)]
                erat = float(np.sqrt((sig ** 2).mean()))
                nrat = float(np.sqrt((noi ** 2).mean()))
                envs = np.sqrt((wv[ch] ** 2).sum(0))
                from scipy.signal import hilbert
                env = np.abs(hilbert(wv[ch], axis=1)).sum(0)
                pk = float(env.max() / np.median(env))
                w12 = wv[ch][:, PRE:PRE + n].ravel()
                kur = float(((w12 - w12.mean()) ** 4).mean() / max(w12.var(), 1e-12) ** 2)
                widx = len(wav); wav.append(wv.astype(np.float16))
            else:
                erat = nrat = pk = kur = np.nan; widx = -1
            recs.append((r.row_id, s, site[s], cf[s], cf_, cr_, best, blag, erat, nrat, pk, kur, widx, j))
        rowinfo.append((r.row_id, len(ss), net_re))
    Rc = pd.DataFrame(recs, columns=["row_id", "sta", "site", "coef", "cc_fwd", "cc_rev", "cc_fwd_best", "lag_best",
                                     "rms_sig", "rms_pre", "env_pk", "kurt", "widx", "j"])
    Ri = pd.DataFrame(rowinfo, columns=["row_id", "n_sta_used", "net_re"])
    np.savez(out_p + ".tmp.npz", wav=np.stack(wav) if wav else np.zeros((0, 3, LEN), np.float16))
    os.replace(out_p + ".tmp.npz", out_p)
    Rc.to_parquet(f"{WK}/{ds}_win.parquet", index=False)
    Ri.to_parquet(f"{WK}/{ds}_row.parquet", index=False)
    return ds, len(Rc)


def build_rows():
    N = pd.read_parquet(f"{OD}/nw_core_candidates.parquet",
                        columns=["date", "series", "fid", "i", "t", "net_cc", "net_mad", "P_own_WAVE", "quiet", "tremor30", "hfr_rel"])
    rng = np.random.default_rng(20261006)
    q = N[N.quiet & (N.P_own_WAVE >= 0.75)].copy(); q["group"] = "quiet"
    t = N[N.tremor30 & (N.P_own_WAVE >= 0.75)].copy(); t["group"] = "tremor"
    # tremor-day reference: all control rows; data rows on 200 sampled tremor days
    td = np.sort(t[t.series == "data"].date.unique())
    keep = set(rng.choice(td, size=min(200, len(td)), replace=False))
    t = t[(t.series != "data") | t.date.isin(keep)]
    S = pd.concat([q, t], ignore_index=True)
    # random origins: for every data row, one uniform random origin of the same family and day (pure-noise reference)
    dd = S[S.series == "data"]
    Rr = dd[["date", "fid", "group"]].copy()
    Rr["i"] = rng.integers(0, NV, len(Rr)); Rr["series"] = "random"
    Rr["t"] = pd.to_datetime(Rr.date).values.astype("datetime64[ns]").astype(np.int64) / 1e9 + Rr.i / FS
    S = pd.concat([S, Rr], ignore_index=True)
    S["row_id"] = np.arange(len(S))
    S.to_parquet(f"{OD}/rows.parquet", index=False)
    return S


if __name__ == "__main__":
    os.makedirs(WK, exist_ok=True)
    S = build_rows()
    print(S.groupby(["group", "series"]).size())
    jobs = [(ds, g[["row_id", "fid", "series", "i", "group"]]) for ds, g in S.groupby("date")]
    jobs.sort(key=lambda x: -len(x[1]))
    t0 = time.time()
    with Pool(int(os.environ.get("QW_NPROC", "24")), initializer=init) as pool:
        for k, (ds, res) in enumerate(pool.imap_unordered(do_day, jobs)):
            if k % 100 == 0:
                print(f"{k}/{len(jobs)} {ds} {res} {time.time() - t0:.0f}s", flush=True)
    print("done", time.time() - t0)
