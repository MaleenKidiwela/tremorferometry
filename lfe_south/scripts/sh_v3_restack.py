#!/usr/bin/env python
"""SH_V3 template stacking per Shelly (2025) Text S1 (2026-09-29).  New file; v2 files untouched.

Text S1, literally:
  * stack the 100 detections with the highest correlation sum; if a channel lacks data for some of
    them, take further detections in order of strength (max 500) until 100 have data ON THAT CHANNEL
  * each waveform may shift by <= 0.1 s against the stack; only the 50 (of 100) events best
    correlated with the stack enter each channel's final stack (Shelly & Hardebeck 2010)
  * templates are 6 s long, beginning 1 s before the S arrival picked on the stack
  * station set: his 15 detection stations + BK.RVIT, BK.HUNT (restack) -- only those in our cache
Our implementation (deviations flagged in SH_V3_LOG.md):
  * S "pick" = automatic, on a 25 s stack of the top-100 (no shifts), horizontal channels (all channels
    for 1-C): 0.5 s smoothed energy envelope, peak within +-3 s of the predicted S (iasp91 tt_S at
    his Table S1 location + the family's catalogue-time constant = median stage-1 tau), onset = last
    sample before the peak below noise + 0.25 (peak - noise), noise = median envelope 12..7 s before
    the predicted S.  pick_q = peak / noise is written (rank, don't discard: every station kept).
  * 20 sps cache: 0.1 s = +-2 samples; per-channel shifts (his "waveforms"); leave-one-out cc for the
    best-50 ranking; shifts re-estimated twice.
  * every window is read through sh_v2_common.DayPrep (v2 cleaning: dead floors, burst muting, line
    notching) and must pass window_ok (gap guard) on its 6.2 s span to count as "having data".

  --mode his : seeds = his published events on the family's OWN build days, ranked by his ccsum;
               S picked here (iteration 0)
  --mode det : seeds = our build-day detections (t11_mf_v2 parquet) at >= --thr x MAD and summed cc
               >= --floor, ranked by summed cc; S picks taken from --picks (the it0 table)
"""
import os, sys, argparse, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sh_v2_common import DayPrep, coords, epi_km, window_ok, FS, log, GUARD_REL
from sh_v3_common import *
from t4_stack_test import taup_S

S1_STATIONS = ["BK.BLAS", "BK.BRIC", "BK.GHOP", "BK.HALS", "BK.JEWT", "BK.PRDS", "BK.SPRL", "NC.KMPB", "NC.KMR",
               "NC.KPR", "NC.KBU", "PB.B049", "PB.B932", "PB.B934", "PB.B933", "BK.RVIT", "BK.HUNT"]
W_PRE = 12.0; W_LEN = 29.0            # extraction window around the predicted / picked S
LW = int(round(W_LEN * FS)); I_S = int(round(W_PRE * FS))
_G = {}
def _init(G): _G.update(G)


def _day(d):
    warnings.filterwarnings("ignore")
    E, OFF, stas = _G["E"], _G["OFF"], _G["stas"]
    out = []
    P = DayPrep(d)
    if P.z is None: return out
    Ed = E[E.day == d]
    for ns in stas:
        fams = [f for f in Ed.fam.unique() if (f, ns) in OFF]
        if not fams: continue
        v = P.load(ns)
        if v is None: continue
        t0, X, mc, ref = v
        try:
            ch = [str(c) for c in np.load(npz_path(*ns.split("."), d))["chans"]]
        except Exception:
            ch = [str(i) for i in range(X.shape[0])]
        for f in fams:
            Ef = Ed[Ed.fam == f]
            tt = Ef.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
            i0 = np.round((tt + OFF[(f, ns)] - W_PRE - t0) * FS).astype(int)
            for r, a in zip(Ef["rank"].values, i0):
                if a < 0 or a + LW > X.shape[1]: continue
                out.append((f, ns, int(r), X[:, a:a + LW].astype(np.float32), mc[:, a:a + LW].copy(), ref.copy(), tuple(ch)))
    return out


def smooth(x, n):
    k = np.ones(n) / n
    return np.convolve(x, k, mode="same")


def pick_S(stk, chans):
    """stk (C, LW) stack around the predicted S (index I_S).  Returns (pick index, q)."""
    hz = [i for i, c in enumerate(chans) if not str(c).endswith("Z")]
    use = hz if hz else list(range(stk.shape[0]))
    env = smooth((stk[use] ** 2).sum(0), int(0.5 * FS))
    noise = float(np.median(env[I_S - int(12 * FS):I_S - int(7 * FS)])) + 1e-30
    a, b = I_S - int(3 * FS), I_S + int(3 * FS)
    pk = a + int(np.argmax(env[a:b])); peak = float(env[pk])
    lev = noise + 0.25 * (peak - noise)
    j = pk
    lo = I_S - int(4 * FS)
    while j > lo and env[j - 1] >= lev: j -= 1
    return j, peak / noise


def ncc(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float((a * b).sum() / (np.sqrt((a * a).sum() * (b * b).sum()) + 1e-30))


def s1_stack(W, ok, s0, L, m, nstack, nbest):
    """W (n, LW) one channel, candidates in rank order; ok (n,) has data on [s0-m, s0+L+m).
    Returns (template (L,), n_used, mean LOO cc of the kept, median cc of all)."""
    idx = np.flatnonzero(ok)[:nstack]
    if len(idx) < 4: return None
    seg = W[idx, s0 - m:s0 + L + m].astype(np.float64)
    seg = seg / (np.sqrt((seg[:, m:m + L] ** 2).mean(1, keepdims=True)) + 1e-30)
    sh = np.zeros(len(idx), int)
    cur = lambda: np.stack([seg[i, m + sh[i]:m + sh[i] + L] for i in range(len(idx))])
    for _ in range(2):
        C = cur(); stk = C.mean(0)
        for i in range(len(idx)):
            best = -2; bs = 0
            for s in range(-m, m + 1):
                c = ncc(seg[i, m + s:m + s + L], stk)
                if c > best: best, bs = c, s
            sh[i] = bs
    C = cur(); tot = C.sum(0); n = len(idx)
    cc = np.array([ncc(C[i], (tot - C[i]) / max(n - 1, 1)) for i in range(n)])
    keep = np.argsort(-cc)[:max(min(nbest, n // 2), 2)]
    return C[keep].mean(0), n, float(cc[keep].mean()), float(np.median(cc))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["his", "det"]); ap.add_argument("--det", default="")
    ap.add_argument("--picks", default="", help="it0 template csv (det mode): S offsets per fam/station")
    ap.add_argument("--days-file", default=BUILD); ap.add_argument("--out", required=True)
    ap.add_argument("--nstack", type=int, default=100); ap.add_argument("--nmax", type=int, default=500)
    ap.add_argument("--nbest", type=int, default=50); ap.add_argument("--shift", type=float, default=0.1)
    ap.add_argument("--thr", type=float, default=8.0); ap.add_argument("--floor", type=float, default=CCSUM_MIN)
    ap.add_argument("--workers", type=int, default=4); ap.add_argument("--max-days", type=int, default=0)
    ap.add_argument("--v2base", default=f"{RES}/sh_v2/tpl_it0.csv")
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    CO = coords()
    stas = [s for s in S1_STATIONS if s in CO and os.path.isdir(f"{CACHE}/{s}")]
    log(f"[restack-v3] S1 stations in cache: {len(stas)}/{len(S1_STATIONS)} {stas}; "
        f"missing {[s for s in S1_STATIONS if s not in stas]}")
    FAM = pd.read_csv(f"{RES}/families_south.csv"); FAM = FAM[FAM.source == SRC].copy(); FAM["fam"] = FAM.fam.astype(str)
    V2 = pd.read_csv(a.v2base); V2["fam"] = V2.fam.astype(str)
    TAU = V2.groupby("fam").best_tau.median().to_dict()                 # his catalogue-time constant
    ST1 = {(r.fam, r.netsta): float(r.tt_s) + float(r.best_tau) for r in V2.itertuples()}
    rows0 = []
    for r in FAM.itertuples():
        for ns in stas:
            la, lo = CO[ns]; e = float(epi_km(r.lat, r.lon, la, lo))
            tt = taup_S(float(r.depth), e)
            rows0.append(dict(fam=r.fam, netsta=ns, tt_s=round(tt, 3), hypo_km=round(float(np.hypot(e, r.depth)), 1),
                              pred=tt + TAU[r.fam], stage1_off=ST1.get((r.fam, ns), np.nan)))
    G0 = pd.DataFrame(rows0)
    if a.mode == "det":
        PK = pd.read_csv(a.picks); PK["fam"] = PK.fam.astype(str)
        PKd = {(r.fam, r.netsta): float(r.tt_s) + float(r.best_tau) for r in PK.itertuples()}
        G0["pred"] = [PKd.get((f, s), np.nan) for f, s in zip(G0.fam, G0.netsta)]
        G0 = G0[np.isfinite(G0.pred)]
    OFF = {(r.fam, r.netsta): float(r.pred) for r in G0.itertuples()}
    FD = family_days(a.days_file)
    if a.mode == "his":
        S = his_catalogue(); E = S[["fam", "t", "ccsum"]].copy()
    else:
        D = pd.read_parquet(a.det); D["fam"] = D.fam.astype(str)
        E = D[(D.net_mad >= a.thr) & (D.ccsum >= a.floor)][["fam", "t", "ccsum"]].copy()
    E["day"] = E.t.dt.floor("D")
    E = E[[f in FD and d in FD[f] for f, d in zip(E.fam, E.day)]]
    E = E.sort_values("ccsum", ascending=False).groupby("fam").head(a.nmax).copy()
    E["rank"] = E.groupby("fam").cumcount()
    days = sorted(E.day.unique())
    if a.max_days: days = days[:a.max_days]; E = E[E.day.isin(days)]
    log(f"[restack-v3] mode={a.mode}: {len(E)} candidate events (<= {a.nmax}/family, own build days) on {len(days)} days; "
        f"per family min/median {E.groupby('fam').size().min()}/{E.groupby('fam').size().median():.0f}")
    from concurrent.futures import ProcessPoolExecutor
    WIN = {}
    with ProcessPoolExecutor(a.workers, initializer=_init, initargs=(dict(E=E, OFF=OFF, stas=stas),)) as ex:
        for i, res in enumerate(ex.map(_day, days, chunksize=2)):
            for f, ns, r, W, M, ref, ch in res:
                WIN.setdefault((f, ns), []).append((r, W, M, ref, ch))
            if i % 50 == 0: log(f"  {i + 1}/{len(days)} days, pairs {len(WIN)}")
    L = int(round(6.0 * FS)); m = int(round(a.shift * FS)); PRE = 1.0
    TPL = {}; rows = []
    for g in G0.itertuples():
        k = (g.fam, g.netsta); lst = WIN.get(k)
        row = dict(source=SRC, fam=g.fam, net=g.netsta.split(".")[0], sta=g.netsta.split(".")[1], netsta=g.netsta,
                   key=f"{SRC}|{g.fam}|{g.netsta}", tt_s=g.tt_s, hypo_km=g.hypo_km, pre_s=PRE, len_s=6.0,
                   stage1_off=g.stage1_off, pred_off=g.pred)
        if not lst:
            row.update(status="no_data"); rows.append(row); continue
        nc = pd.Series([w[1].shape[0] for w in lst]).mode().iloc[0]
        lst = sorted([w for w in lst if w[1].shape[0] == nc], key=lambda w: w[0])
        W = np.stack([w[1] for w in lst]); M = np.stack([w[2] for w in lst]); R = np.stack([w[3] for w in lst])
        chans = lst[0][4]
        if a.mode == "his":
            # long stack (top nstack with data over -10..+15 s, no shifts) -> S pick
            a0, a1 = I_S - int(10 * FS), I_S + int(15 * FS)
            stk = np.zeros((nc, LW)); q_ = []
            for c in range(nc):
                okc = M[:, c, a0:a1].all(1) & (R[:, c] > 0)
                ii = np.flatnonzero(okc)[:a.nstack]
                if len(ii) < 4: continue
                w = W[ii, c].astype(np.float64); w = w / (np.sqrt((w[:, a0:a1] ** 2).mean(1, keepdims=True)) + 1e-30)
                stk[c] = w.mean(0)
            if not np.any(stk):
                row.update(status="no_long_stack"); rows.append(row); continue
            j, q = pick_S(stk, chans)
            off = g.pred + (j - I_S) / FS; row.update(pick_dt=(j - I_S) / FS, pick_q=q)
        else:
            j = I_S; off = g.pred
        s0 = j - int(PRE * FS)
        if s0 - m < 0 or s0 + L + m > LW:
            row.update(status="pick_at_edge"); rows.append(row); continue
        T = np.zeros((nc, L)); ninfo = []; ccinfo = []
        for c in range(nc):
            seg_ok = M[:, c, s0 - m:s0 + L + m].all(1) & (R[:, c] > 0)
            e = (W[:, c, s0:s0 + L].astype(np.float64) ** 2).sum(1)
            seg_ok &= e > GUARD_REL * R[:, c]
            res = s1_stack(W[:, c], seg_ok, s0, L, m, a.nstack, a.nbest)
            if res is None: ninfo.append(0); ccinfo.append(np.nan); continue
            T[c], n, cck, ccm = res; ninfo.append(n); ccinfo.append(cck)
        if not np.any(T):
            row.update(status="no_template"); rows.append(row); continue
        TPL[row["key"]] = T.astype(np.float32)
        row.update(status="ok", ncomp=int(nc), chans=",".join(map(str, chans)), best_tau=off - g.tt_s,
                   n_stack=int(max(ninfo)), n_stack_ch=",".join(map(str, ninfo)),
                   cc50_ch=",".join(f"{x:.3f}" for x in ccinfo), cc50=float(np.nanmean(ccinfo)),
                   n_cand=len(lst), src="stacked")
        rows.append(row)
    T = pd.DataFrame(rows)
    T.to_csv(a.out.replace(".npz", "_all.csv"), index=False)
    Tok = T[T.status == "ok"].copy(); Tok.to_csv(a.out.replace(".npz", ".csv"), index=False)
    np.savez_compressed(a.out, **TPL)
    log(f"[restack-v3] {len(TPL)} templates over {Tok.fam.nunique()} families, {Tok.netsta.nunique()} stations; "
        f"status {T.status.value_counts().to_dict()}; median n_stack {Tok.n_stack.median():.0f}, median cc50 {Tok.cc50.median():.3f}"
        + (f"; |S pick - stage-1 offset| median {np.nanmedian(np.abs(Tok.best_tau + Tok.tt_s - Tok.stage1_off)):.2f} s" if a.mode == "his" else ""))
if __name__ == "__main__": main()
