#!/usr/bin/env python
"""ALL-DAYS DETECTION PASS for the Royer & Bostock (2014) reproduction.  Prefix `ab_alldays`.

USER REQUIREMENT (2026-09-25): detections are NOT gated to tremor / ETS windows.  LFEs may occur
outside slow-slip windows or after ETS, and that could be a finding.  Template building (growth,
ab_grow.py) stays literal on tremor days; this pass then scans EVERY day on which a final
template's stations recorded, in the array lifetime (ab_alldays_fetch.LIFETIME).  No day selection.

What is scanned
  * FINAL templates of a subarray = `tpl_it<final_iter>.npz` of the chain's growth
    (results/ab_royer/<SUB>/grow/summary.json -> final_iter), each with the SAME 7-10 station
    subset (ranked by SNR) growth used in that iteration; every template in the file is scanned.
  * `--kind nullK` = the chain's whole-chain site-shift null realisation K
    (results/ab_royer/<SUB>/null<K>/grow).  Identical treatment: the null's final templates are run
    through the SAME pass over the SAME real (unshifted) continuous data, AB_NULL is forced to 0.
    So real vs null compares real-data families with noise-built families on identical days.

Detector = ab_grow.matched_filter, unchanged in every number: per-CHANNEL normalised cc (masked:
a channel-sample whose window touches an invalid sample, or has zero energy, contributes nothing;
the denominator is +inf there, never clamped on zero-filled gaps -- results/WA_BUG_AUDIT.md), the
network function = mean over valid channels, samples with < 4 valid template stations excluded
from the day's pool, threshold median + 8 x MAD over the 24-hour day, dedupe = keep a sample only
if it is the maximum within +-15 s (the template length) of the same template.  The only change is
WHERE the zero-padded template array is built (on the GPU instead of a 2 GB host array per batch
and station): the same float64 values enter the same torch FFTs.  `--verify` checks the output is
identical to ab_grow.matched_filter on the given days.

Output (never filtered; labels are COLUMNS):
  results/ab_royer/alldays/<SUB>/<kind>/days/<date>.parquet   per-day shards (resumable)
  results/ab_royer/alldays/<SUB>/<kind>/days.csv              per scanned day: stations, detections
  results/ab_royer/alldays/<SUB>/detections_<kind>.parquet    family, time, cc, nmad, n_channels,
      + t (epoch s), cc_mean, n_sta, thr, date, n_det_growth (family's detections in the final
      growth iteration), and the LABELS tremor_day, label_source, in_ets_window,
      days_since_last_ets_end, days_to_next_ets_start, ets_phase, ets_rel_days
  results/ab_royer/alldays/<SUB>/ets_windows.csv, day_labels.csv

TREMOR-DAY LABELS (per subarray, same rule as ab_tremor_days.py)
  * from 2009-08-06: PNSN catalogue day = >= 24 epicentres (2.5 min each -> >= 1 h of tremor)
    within 25 km of any station of the subarray's pool (label_source "pnsn");
  * before 2009-08-06: our Wech & Creager (2008)-style detector's day labels
    (results/ab_royer/tremor/wecc_days_<SUB>.csv, calibrated on the PNSN era; "wecc");
  * NVI (2005-2007): no catalogue and no WECC run; the paper's two sec.2.2 episodes are the only
    known tremor days ("paper_sec2.2"), every other NVI day is UNKNOWN;
  * anything else: UNKNOWN (NaN) -- never counted as a non-tremor day.
ETS WINDOWS (documented rule, `ets_windows()`): tremor days are merged into runs when separated by
<= GAP_DAYS (3) non-tremor days; a run is an ETS window when it holds >= MIN_TREMOR_DAYS (5)
tremor days.  Where no day labels exist, the PNSN published ETS windows used by
ab_tremor_days.ETS_WINDOWS are added (source "pnsn_published"), and for NVI the paper's episodes
(source "paper_sec2.2").  A window runs from 00:00 of its first to 24:00 of its last day.
  in_ets_window           True inside a window; False on a LABELLED day outside every window;
                          NaN on an unlabelled day outside every window
  days_since_last_ets_end (t - end of the latest window that ENDED before t) in days; NaN if none
  days_to_next_ets_start  (start of the first window that STARTS after t - t) in days; NaN if none
  ets_phase               during | post (closer to the previous end) | pre (closer to the next
                          start) | unknown
  ets_rel_days            0 during; +days_since_last_ets_end for post; -days_to_next_ets_start for pre

  python scripts/ab_alldays.py run  --sub NW --kind real
  python scripts/ab_alldays.py run  --sub NW --kind null1
  python scripts/ab_alldays.py run  --sub NW --tpl smoke/grow/tpl_it1.npz --out <dir> --dates a,b
  python scripts/ab_alldays.py run  ... --verify           # identical to ab_grow.matched_filter?
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

os.environ["AB_NULL"] = "0"          # the all-days pass always reads REAL, unshifted data

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402

OUT = f"{A.RES}/alldays"
LIFETIME = {"NW": ("2006-07-21", "2012-05-31"), "SW": ("2006-07-09", "2012-05-31"),
            "NVI": ("2005-06-01", "2007-08-31")}            # = ab_alldays_fetch.LIFETIME
TPL_LEN = 15.0                  # = ab_grow.TPL_LEN (the chain's default --tpl-len)
MIN_STA_VALID = 4               # = ab_grow.MIN_STA_VALID
PNSN_START = "2009-08-06"
RADIUS_KM, MIN_EPI = 25.0, 24   # = ab_tremor_days.RADIUS_KM / MIN_EPI
GAP_DAYS, MIN_TREMOR_DAYS = 3, 5


# =============================================================== templates ======================
def final_dir(sub, kind):
    return f"{A.RES}/{sub}/grow" if kind == "real" else f"{A.RES}/{sub}/{kind}/grow"


def final_tpl_path(sub, kind):
    d = final_dir(sub, kind)
    sj = f"{d}/summary.json"
    if not os.path.exists(sj):
        return None, None
    S = json.load(open(sj))
    it = S.get("final_iter", -1)
    if it is None or it < 0 or not os.path.exists(f"{d}/tpl_it{it}.npz"):
        its = sorted(int(os.path.basename(p)[6:-4]) for p in glob.glob(f"{d}/tpl_it*.npz")
                     if os.path.exists(p.replace("tpl_it", "det_it").replace(".npz", ".parquet")))
        if not its:
            return None, None
        it = its[-1]
    return f"{d}/tpl_it{it}.npz", f"{d}/det_it{it}.parquet"


def load_tpl(npz):
    """{g: {sta: (tpl (3,n) float64, off_s, snr)}} exactly as ab_grow.cut_templates produced it."""
    z = np.load(npz, allow_pickle=True)
    meta = json.loads(str(z["meta"]))
    TPL = {}
    for k in z.files:
        if k == "meta":
            continue
        g, s = k.split("_", 1)
        off, snr = meta[k]
        TPL.setdefault(int(g), {})[s] = (z[k].astype(np.float64), float(off), float(snr))
    return TPL


# =============================================================== detector =======================
def mf_day(D, TPL, date, nmad=A.NMAD, net="mean", min_sta_valid=MIN_STA_VALID, decl_s=TPL_LEN,
           tpl_len=TPL_LEN, dev=None, batch=24, chunk=64):
    """ab_grow.matched_filter for ONE date, given its loaded {sta: (t0, X, m)}.  Every line is the
    chain's; only the zero-padded template array M is built on the device."""
    import torch
    dev = dev or A.torch_dev()
    torch.backends.cuda.matmul.allow_tf32 = False
    nt = int(round(tpl_len * A.FS))
    gs_all = sorted(TPL)
    stas = sorted({s for g in TPL for s in TPL[g]})
    have = [s for s in stas if s in D]
    empty = pd.DataFrame(columns=["grp", "t", "cc", "cc_mean", "nmad", "n_sta", "n_ch", "thr",
                                  "date"])
    if len(have) < min_sta_valid:
        return empty, len(have)
    N = min(D[s][1].shape[1] for s in have)
    NV = N - nt + 1
    Xf, DEN, OK, T0 = {}, {}, {}, {}
    for s in have:
        t0, X, m = D[s]
        X = X[:, :N]
        m = m[:N]
        ew, vw = A.window_stats(X, m, nt, 1)
        Xt = torch.from_numpy(np.ascontiguousarray(X, np.float64)).to(dev)
        Xf[s] = torch.fft.rfft(Xt, n=N)
        e = torch.from_numpy(ew).to(dev)
        okc = torch.from_numpy(vw).to(dev)[None] & (e > 0)
        DEN[s] = torch.where(okc, e.clamp_min(0).sqrt(),
                             torch.full_like(e, float("inf"))).float()
        OK[s] = okc
        T0[s] = float(t0)
        del Xt, e
    # template tensors once per day (float64, identical values to ab_grow's host array)
    TT = {(g, s): torch.from_numpy(np.ascontiguousarray(TPL[g][s][0], np.float64)).to(dev)
          for g in gs_all for s in TPL[g] if s in D}
    out = []
    for c0 in range(0, len(gs_all), chunk):
        gs = gs_all[c0:c0 + chunk]
        acc = torch.zeros((len(gs), NV), device=dev)
        cch = torch.zeros((len(gs), NV), device=dev)
        csta = torch.zeros((len(gs), NV), device=dev)
        for b0 in range(0, len(gs), batch):
            sub = gs[b0:b0 + batch]
            for s in have:
                idx = [i for i, g in enumerate(sub) if s in TPL[g]]
                if not idx:
                    continue
                M = torch.zeros((len(idx), 3, N), dtype=torch.float64, device=dev)
                offs, nrm = [], []
                for r, i in enumerate(idx):
                    tpl, off, _ = TPL[sub[i]][s]
                    M[r, :, :nt] = TT[(sub[i], s)]
                    offs.append(off)
                    nrm.append(np.sqrt((tpl ** 2).sum(1)))
                Tf = torch.fft.rfft(M, n=N)
                del M
                cc = torch.fft.irfft(Xf[s][None] * Tf.conj(), n=N)[..., :NV].float()
                del Tf
                nrmt = torch.tensor(np.array(nrm), device=dev, dtype=torch.float32)
                tok = (nrmt > 0)[:, :, None] & OK[s][None]
                cc = cc / (nrmt[:, :, None] * DEN[s][None]).clamp_min(1e-30)
                cc = torch.where(tok, cc, torch.zeros_like(cc))
                ccs = cc.sum(1)
                nch = tok.float().sum(1)
                del cc, tok
                for r, i in enumerate(idx):
                    k = int(round(offs[r] * A.FS))
                    lo, hi = max(-k, 0), min(NV - k, NV)
                    if hi <= lo:
                        continue
                    acc[b0 + i, lo:hi] += ccs[r, lo + k:hi + k]
                    cch[b0 + i, lo:hi] += nch[r, lo + k:hi + k]
                    csta[b0 + i, lo:hi] += (nch[r, lo + k:hi + k] > 0).float()
                del ccs, nch
            if dev == "cuda":
                torch.cuda.empty_cache()
        enough = csta >= float(min_sta_valid)
        fn = acc if net == "sum" else acc / cch.clamp_min(1.0)
        fn = torch.where(enough, fn, torch.full_like(fn, float("nan")))
        med = torch.nanmedian(fn, dim=1, keepdim=True).values
        madv = torch.nanmedian((fn - med).abs(), dim=1, keepdim=True).values
        thr = med + nmad * madv.clamp_min(1e-9)
        dec = int(round(decl_s * A.FS))
        hit = torch.nan_to_num(fn, nan=-1e9) >= thr
        vals = torch.where(hit, fn, torch.full_like(fn, -1e9))
        pooled = torch.nn.functional.max_pool1d(vals[None], 2 * dec + 1, stride=1, padding=dec)[0]
        keep = hit & (vals >= pooled - 1e-9)
        ij = torch.nonzero(keep, as_tuple=False)
        if ij.numel():
            gi = ij[:, 0].cpu().numpy()
            ti = ij[:, 1].cpu().numpy()
            v = fn[ij[:, 0], ij[:, 1]].cpu().numpy()
            nm = ((fn[ij[:, 0], ij[:, 1]] - med[ij[:, 0], 0])
                  / madv[ij[:, 0], 0].clamp_min(1e-9)).cpu().numpy()
            ns = csta[ij[:, 0], ij[:, 1]].cpu().numpy().astype(np.int16)
            nc = cch[ij[:, 0], ij[:, 1]].cpu().numpy().astype(np.int16)
            ccm = (acc[ij[:, 0], ij[:, 1]] / cch[ij[:, 0], ij[:, 1]].clamp_min(1.0)).cpu().numpy()
            out.append(pd.DataFrame(dict(grp=[gs[i] for i in gi], t=T0[have[0]] + ti / A.FS,
                                         cc=v, cc_mean=ccm, nmad=nm, n_sta=ns, n_ch=nc,
                                         thr=thr[ij[:, 0], 0].cpu().numpy(), date=date)))
        del acc, cch, csta, fn, hit, vals, pooled, keep, ij
    del Xf, DEN, OK, TT
    if dev == "cuda":
        torch.cuda.empty_cache()
    return (pd.concat(out, ignore_index=True) if out else empty), len(have)


# =============================================================== labels =========================
def _pnsn_days(sub):
    import ab_tremor_days as TD
    T = TD.catalogue()
    T = T[T.time >= PNSN_START]
    la = np.array([A.COORD[s][0] for s in A.POOL[sub]])
    lo = np.array([A.COORD[s][1] for s in A.POOL[sub]])
    dk = A.epi_km(T.lat.to_numpy()[:, None], T.lon.to_numpy()[:, None], la, lo).min(1)
    n = T[dk <= RADIUS_KM].groupby(T.time[dk <= RADIUS_KM].dt.strftime("%Y-%m-%d")).size()
    return n, str(T.time.max())[:10]


def day_labels(sub, a, b):
    """DataFrame date, tremor_day (1/0/NaN), label_source, n_epi, wecc_windows."""
    D = pd.DataFrame(dict(date=pd.date_range(a, b).strftime("%Y-%m-%d")))
    D["tremor_day"] = np.nan
    D["label_source"] = "none"
    D["n_epi"] = np.nan
    D["wecc_windows"] = np.nan
    if sub in ("NW", "SW"):
        wp = f"{A.RES}/tremor/wecc_days_{sub}.csv"
        if os.path.exists(wp):
            W = pd.read_csv(wp)
            W["date"] = W.date.astype(str)
            W = W[W.date < PNSN_START].set_index("date")
            m = D.date.isin(W.index)
            D.loc[m, "tremor_day"] = W.loc[D.date[m], "tremor_day"].astype(float).to_numpy()
            D.loc[m, "label_source"] = "wecc"
            if "n_tremor_windows" in W:
                D.loc[m, "wecc_windows"] = W.loc[D.date[m], "n_tremor_windows"].to_numpy()
        n, cat_end = _pnsn_days(sub)
        m = (D.date >= PNSN_START) & (D.date <= cat_end)
        D.loc[m, "n_epi"] = D.date[m].map(n).fillna(0).to_numpy()
        D.loc[m, "tremor_day"] = (D.loc[m, "n_epi"] >= MIN_EPI).astype(float)
        D.loc[m, "label_source"] = "pnsn"
    else:
        ep = {x.strftime("%Y-%m-%d") for s_, e_ in A.SUBARRAY[sub]["episodes"]
              for x in pd.date_range(s_, e_)}
        m = D.date.isin(ep)
        D.loc[m, "tremor_day"] = 1.0
        D.loc[m, "label_source"] = "paper_sec2.2"
    return D


def ets_windows(sub, L):
    """ETS windows from the tremor-day series L (day_labels), rule in the module docstring."""
    td = L[L.tremor_day == 1].date.tolist()
    W = []
    if sub in ("NW", "SW"):
        run = []
        for d in td:
            if run and (pd.Timestamp(d) - pd.Timestamp(run[-1])).days - 1 > GAP_DAYS:
                W.append(run)
                run = []
            run.append(d)
        if run:
            W.append(run)
        W = [dict(start=r[0], end=r[-1], n_tremor_days=len(r), source=L.set_index("date")
                  .loc[r[0], "label_source"]) for r in W if len(r) >= MIN_TREMOR_DAYS]
        import ab_tremor_days as TD
        lab = set(L[L.tremor_day.notna()].date)
        for a, b in TD.ETS_WINDOWS.get(sub, []):
            dd = pd.date_range(a, b).strftime("%Y-%m-%d")
            if not any(x in lab for x in dd) and any(x in set(L.date) for x in dd):
                W.append(dict(start=a, end=b, n_tremor_days=np.nan, source="pnsn_published"))
    else:
        for a, b in A.SUBARRAY[sub]["episodes"]:
            W.append(dict(start=a, end=b, n_tremor_days=len(pd.date_range(a, b)),
                          source="paper_sec2.2"))
    W = pd.DataFrame(W, columns=["start", "end", "n_tremor_days", "source"]).sort_values("start")
    W["duration_days"] = [(pd.Timestamp(e) - pd.Timestamp(s)).days + 1
                          for s, e in zip(W.start, W.end)]
    return W.reset_index(drop=True)


def label_times(t, dates, L, W):
    """Label epoch-second times t (with their UTC dates) against day labels L and windows W."""
    t = np.asarray(t, float)
    ws = np.array([pd.Timestamp(s).timestamp() for s in W.start]) if len(W) else np.zeros(0)
    we = np.array([(pd.Timestamp(e) + pd.Timedelta("1D")).timestamp() for e in W.end]) \
        if len(W) else np.zeros(0)
    Li = L.set_index("date")
    td = pd.Series(dates).map(Li.tremor_day).to_numpy(float)
    src = pd.Series(dates).map(Li.label_source).fillna("none").to_numpy()
    inside = ((t[:, None] >= ws[None]) & (t[:, None] < we[None])).any(1) if len(W) else \
        np.zeros(t.size, bool)
    since = np.full(t.size, np.nan)
    to = np.full(t.size, np.nan)
    if len(W):
        E = np.where(we[None] <= t[:, None], t[:, None] - we[None], np.inf).min(1)
        S = np.where(ws[None] > t[:, None], ws[None] - t[:, None], np.inf).min(1)
        since = np.where(np.isfinite(E), E / 86400.0, np.nan)
        to = np.where(np.isfinite(S), S / 86400.0, np.nan)
    in_ets = np.where(inside, 1.0, np.where(np.isfinite(td), 0.0, np.nan))
    phase = np.where(inside, "during", "unknown").astype(object)
    rel = np.where(inside, 0.0, np.nan)
    post = ~inside & np.isfinite(since) & (~np.isfinite(to) | (since <= to))
    pre = ~inside & np.isfinite(to) & (~np.isfinite(since) | (to < since))
    phase[post], rel[post] = "post", since[post]
    phase[pre], rel[pre] = "pre", -to[pre]
    return dict(tremor_day=td, label_source=src, in_ets_window=in_ets,
                days_since_last_ets_end=since, days_to_next_ets_start=to,
                ets_phase=phase, ets_rel_days=rel)


def labels_for_sub(sub):
    a, b = LIFETIME[sub]
    L = day_labels(sub, a, b)
    W = ets_windows(sub, L)
    os.makedirs(f"{OUT}/{sub}", exist_ok=True)
    L.to_csv(f"{OUT}/{sub}/day_labels.csv", index=False)
    W.to_csv(f"{OUT}/{sub}/ets_windows.csv", index=False)
    return L, W


# =============================================================== driver =========================
def scan_dates(stas, a, b, min_sta=MIN_STA_VALID):
    rows = []
    for d in pd.date_range(a, b).strftime("%Y-%m-%d"):
        n = sum(os.path.exists(A.cache_path(s, d)) for s in stas)
        rows.append((d, n))
    return [d for d, n in rows if n >= min_sta], dict(rows)


def assemble(sub, kind, odir, det_growth=None, label=True):
    shards = sorted(glob.glob(f"{odir}/days/*.parquet"))
    parts = [pd.read_parquet(p) for p in shards]
    parts = [p for p in parts if len(p)]
    cols = ["grp", "t", "cc", "cc_mean", "nmad", "n_sta", "n_ch", "thr", "date"]
    Dt = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=cols)
    Dt = Dt.rename(columns={"grp": "family", "n_ch": "n_channels"})
    Dt["time"] = pd.to_datetime(Dt.t.astype(float), unit="s").astype("datetime64[ns]")
    if det_growth is not None and os.path.exists(det_growth):
        G = pd.read_parquet(det_growth, columns=["grp"]).groupby("grp").size()
        Dt["n_det_growth"] = Dt.family.map(G).fillna(0).astype(int)
    if label and sub in LIFETIME:
        L, W = labels_for_sub(sub)
        lab = label_times(Dt.t.to_numpy(float), Dt.date.astype(str).to_numpy(), L, W)
        for k, v in lab.items():
            Dt[k] = v
    Dt = Dt[["family", "time", "cc", "nmad", "n_channels"]
            + [c for c in Dt.columns if c not in ("family", "time", "cc", "nmad", "n_channels")]]
    return Dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["run", "assemble", "labels"])
    ap.add_argument("--sub", required=True, choices=["NW", "SW", "NVI"])
    ap.add_argument("--kind", default="real", help="real | null1 | null2 | null3")
    ap.add_argument("--tpl", default="", help="override: templates npz (smoke)")
    ap.add_argument("--out", default="", help="override: output dir (smoke)")
    ap.add_argument("--dates", default="", help="override: comma list (smoke)")
    ap.add_argument("--verify", action="store_true",
                    help="also run ab_grow.matched_filter on each date and compare")
    a = ap.parse_args()
    if a.cmd == "labels":
        L, W = labels_for_sub(a.sub)
        A.log(f"[labels {a.sub}] days {len(L)}: " + L.label_source.value_counts().to_dict().__repr__()
              + f"; tremor days {int((L.tremor_day == 1).sum())}; ETS windows {len(W)}")
        A.log("\n" + W.to_string())
        return
    odir = a.out or f"{OUT}/{a.sub}/{a.kind}"
    os.makedirs(f"{odir}/days", exist_ok=True)
    if a.tpl:
        tpl_npz, det_growth = a.tpl, a.tpl.replace("tpl_it", "det_it").replace(".npz", ".parquet")
    else:
        tpl_npz, det_growth = final_tpl_path(a.sub, a.kind)
        if tpl_npz is None:
            raise SystemExit(f"no final templates for {a.sub} {a.kind} ({final_dir(a.sub, a.kind)})")
    final_out = f"{odir}/detections.parquet" if a.out else f"{OUT}/{a.sub}/detections_{a.kind}.parquet"
    if a.cmd == "run":
        TPL = load_tpl(tpl_npz)
        stas = sorted({s for g in TPL for s in TPL[g]})
        lo, hi = LIFETIME[a.sub]
        if a.dates:
            dates = a.dates.split(",")
            ncache = {d: sum(os.path.exists(A.cache_path(s, d)) for s in stas) for d in dates}
        else:
            dates, ncache = scan_dates(stas, lo, hi)
        json.dump(dict(sub=a.sub, kind=a.kind, tpl=tpl_npz, n_templates=len(TPL), stations=stas,
                       lifetime=[lo, hi], n_dates=len(dates), nmad=A.NMAD, tpl_len=TPL_LEN,
                       min_sta_valid=MIN_STA_VALID, net="mean", AB_NULL=os.environ["AB_NULL"],
                       started=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())),
                  open(f"{odir}/_params.json", "w"), indent=1)
        A.log(f"[alldays {a.sub} {a.kind}] {len(TPL)} templates, {len(stas)} stations, "
              f"{len(dates)} dates ({lo}..{hi}) -> {odir}")
        dc = A.DayCache(stas)
        dayrows_p = f"{odir}/days.csv"
        done = {os.path.basename(p)[:10] for p in glob.glob(f"{odir}/days/*.parquet")}
        t_start = time.time()
        for i, d in enumerate(dates, 1):
            if d in done and not a.verify:
                continue
            t0 = time.time()
            D = dc.day(d)
            det, nhave = mf_day(D, TPL, d)
            det.to_parquet(f"{odir}/days/{d}.parquet.tmp", index=False)
            os.replace(f"{odir}/days/{d}.parquet.tmp", f"{odir}/days/{d}.parquet")
            row = dict(date=d, n_sta_cached=ncache.get(d, np.nan), n_sta_have=nhave,
                       scanned=bool(nhave >= MIN_STA_VALID), n_det=len(det),
                       secs=round(time.time() - t0, 1))
            if a.verify:
                import ab_grow
                ref = ab_grow.matched_filter(dc, TPL, [d])
                same = (len(ref) == len(det)) and (len(det) == 0 or (
                    np.array_equal(ref.grp.to_numpy(), det.grp.to_numpy())
                    and np.allclose(ref.t.to_numpy(float), det.t.to_numpy(float))
                    and np.allclose(ref.cc.to_numpy(float), det.cc.to_numpy(float),
                                    rtol=0, atol=1e-6)))
                row.update(verify_ref_n=len(ref), verify_identical=bool(same))
                A.log(f"[verify] {d}: alldays {len(det)} vs ab_grow {len(ref)} identical={same}")
            pd.DataFrame([row]).to_csv(dayrows_p, mode="a", header=not os.path.exists(dayrows_p),
                                       index=False)
            if i % 10 == 0 or i == len(dates) or a.dates:
                el = time.time() - t_start
                A.log(f"  {d} ({i}/{len(dates)}) {nhave} sta -> {len(det)} det, "
                      f"{row['secs']}s; elapsed {el/3600:.2f} h")
    Dt = assemble(a.sub, a.kind, odir, det_growth)
    Dt.to_parquet(final_out + ".tmp", index=False)
    os.replace(final_out + ".tmp", final_out)
    A.log(f"[alldays {a.sub} {a.kind}] {len(Dt)} detections -> {final_out}")


if __name__ == "__main__":
    main()
