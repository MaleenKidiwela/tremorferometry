#!/usr/bin/env python
"""DUC_FIX re-templating of NEW station epochs (status 'new' in station_epochs.csv).

Rule (his own recipe, compute_new_templates.py: best LFEs by cc, max_LFEs = 150, RMS-normalised
linear stack, window = [detection time, + duration] -- his templates' first sample IS the detection
time): for a new epoch E of station S, the detections are pass-1 detections of each family that
names S, made while E had NO template -- i.e. WITHOUT S -- inside the epoch's template window
W = [E.start, min(E.end, E.start + 24 months)).  Candidates: the family's catalogue detections
(cc x nchannel >= his threshold_perm; for his six eliminated families, threshold 0, all 8xMAD
detections), strongest cc first; the station's cleaned window must be fully valid.
Keep / drop, per (family, channel): N >= 20 windows AND the stack coherence
S = rms(stack) * sqrt(N) (1 for incoherent noise) above every one of 20 SCRAMBLED CONTROLS (the same
detection times, the station's data taken from a different day: whole-day shifts of +-3..40 days,
one shift per control, inside the epoch).  Dropped (family, channel) -> the station is absent from
that family for epoch E.  Columns kept whatever the decision.

Finalisation (called after each pass-1 month): an epoch is finalised when W has ended, or when every
family naming S already has >= 150 candidates.  `--final` finalises everything still pending.
Months from E.start to the finalisation month are appended to rerun_months.txt (pass 2 re-scans them
with the new templates).
  python scripts/duc_fix_retemplate.py --upto 2013-12 [--final] [--workers 4]
"""
import os, sys, json, time, argparse, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from duc_fix_common import *
import duc_fix_clean as CL

NMAX = 150; NMIN = 20; NCTRL = 20; WIN_MONTHS = 24
STATE = f"{FIX}/retemplate_state.csv"


def mf1_dets(m0, m1, fams):
    out = []
    for m in pd.period_range(m0, m1, freq="M"):
        p = f"{FIX}/mf1/{m}.parquet"
        if os.path.exists(p):
            D = pd.read_parquet(p, columns=["fam", "t", "cc", "nchannel", "pass_his", "thr_perm"])
            D["fam"] = D.fam.astype(str); out.append(D[D.fam.isin(fams)])
    return pd.concat(out) if out else pd.DataFrame(columns=["fam", "t", "cc", "nchannel", "pass_his", "thr_perm"])


def candidates(D, ep_start, ep_end):
    D = D[(D.t >= ep_start) & (D.t < ep_end)]
    D = D[(D.pass_his == 1) | (D.thr_perm <= 0)]
    return D.sort_values("cc", ascending=False)


def build_epoch(args):
    ep, famlist, dur, wend, his_tpl, C = args
    warnings.filterwarnings("ignore")
    ns = f"{ep['net']}.{ep['sta']}"
    t0e, t1e = pd.Timestamp(ep["start"]), pd.Timestamp(ep["end"])
    chans = ep["chans"].split(",")
    import zlib
    rng = np.random.default_rng(zlib.crc32(ep["epoch"].encode()))
    shifts = []
    while len(shifts) < NCTRL:
        k = int(rng.integers(3, 41)) * int(rng.choice([-1, 1]))
        if k not in shifts: shifts.append(k)
    cache = {}
    def day(d):
        d = pd.Timestamp(d).normalize()
        if d not in cache:
            if len(cache) > 64: cache.pop(next(iter(cache)))
            v = None
            if t0e <= d < t1e:
                v = CL.DayPrep(d).load(ns)
            cache[d] = v
        return cache[d]
    def window(v, t, L, ci):
        t0, X, mc, ref, chs, rec = v
        if chans[ci] not in chs: return None
        j = chs.index(chans[ci]); i0 = int(round((t - t0) * FS))
        if i0 < 0 or i0 + L > X.shape[1] or not mc[j, i0:i0 + L].all() or ref[j] <= 0: return None
        w = X[j, i0:i0 + L].astype(np.float64); r = np.sqrt(np.mean(w * w))
        return w / r if r > 0 else None
    rows = []; tpl = {}
    for f in famlist:
        L = int(round(dur[f] * FS)) + 1
        Cf = C[C.fam == f].head(3 * NMAX)
        for ci, cha in enumerate(chans):
            W = []; T_ = []
            for r in Cf.itertuples():
                te = r.t.value / 1e9
                v = day(r.t)
                if v is None: continue
                w = window(v, te, L, ci)
                if w is None: continue
                W.append(w); T_.append(te)
                if len(W) >= NMAX: break
            N = len(W)
            row = dict(fam=f, epoch=ep["epoch"], sta=ep["sta"], cha=cha, key=f"{f}|{ep['epoch']}|{cha}",
                       win_start=str(t0e.date()), win_end=str(min(t1e, wend).date()), n_cand=int(len(C[C.fam == f])),
                       N=N, S=np.nan, ctrl_p50=np.nan, ctrl_max=np.nan, cc_his=np.nan, keep=0)
            if N >= 2:
                st = np.mean(W, 0); Sr = float(np.sqrt(np.mean(st * st)) * np.sqrt(N))
                Sc = []
                for k in shifts:
                    Wc = []
                    for te in T_:
                        tc = te + 86400.0 * k
                        v = day(pd.Timestamp(tc, unit="s"))
                        if v is None: continue
                        w = window(v, tc, L, ci)
                        if w is not None: Wc.append(w)
                    if len(Wc) >= 2:
                        sc = np.mean(Wc, 0); Sc.append(float(np.sqrt(np.mean(sc * sc)) * np.sqrt(len(Wc))))
                row.update(S=Sr, ctrl_p50=float(np.median(Sc)) if Sc else np.nan,
                           ctrl_max=float(np.max(Sc)) if Sc else np.nan, n_ctrl=len(Sc))
                h = his_tpl.get((f, cha[-1]))
                if h is not None and h.size == st.size:
                    a = (st - st.mean()) / (st.std() + 1e-30); b = (h - h.mean()) / (h.std() + 1e-30)
                    row["cc_his"] = float(max(np.dot(np.roll(a, s), b) / a.size for s in range(-20, 21)))
                row["keep"] = int(N >= NMIN and len(Sc) >= 10 and Sr > max(Sc))
                if row["keep"]: tpl[row["key"]] = st
            rows.append(row)
    return rows, tpl


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--upto", required=True)
    ap.add_argument("--final", action="store_true"); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    E = epochs(); fp = families(); thr = thresholds()
    dur = {r.family: float(r.duration) for r in fp.itertuples()}
    famst = {r.family: r.stations.split(",") for r in fp.itertuples()}
    HT = his_templates()
    E = E[(E.status == "new") & E.sta.isin(set(",".join(fp.stations).split(",")))]
    st = pd.read_csv(STATE) if os.path.exists(STATE) else pd.DataFrame(columns=["epoch", "state", "final_month"])
    done = set(st[st.state == "done"].epoch)
    mend = pd.Timestamp(a.upto + "-01") + pd.offsets.MonthBegin(1)
    todo = []
    pend = [r for r in E.itertuples() if r.epoch not in done and r.start < mend]
    if not pend:
        log(f"[retemplate {a.upto}] nothing pending"); return
    allf = sorted({f for r in pend for f, s in famst.items() if r.sta in s})
    DA = mf1_dets(min(r.start for r in pend).to_period("M"), (mend - pd.Timedelta(days=1)).to_period("M"), allf)
    for r in pend:
        wend = min(r.end, r.start + pd.DateOffset(months=WIN_MONTHS))
        fl = [f for f, s in famst.items() if r.sta in s]
        C = candidates(DA[DA.fam.isin(fl)], r.start, min(wend, mend))
        cnt = C.groupby("fam").size()
        ready = a.final or mend >= wend or all(cnt.get(f, 0) >= NMAX for f in fl)
        if ready:
            ht = {}
            for f in fl:
                for (s2, c), v in HT.get(f, {}).items():
                    if s2 == r.sta: ht[(f, c[-1])] = v
            Cs = pd.concat([C[C.fam == f].head(3 * NMAX) for f in fl]) if len(C) else C
            todo.append((dict(r._asdict()), fl, dur, min(wend, mend), ht, Cs))
    if not todo:
        log(f"[retemplate {a.upto}] nothing to finalise"); return
    t0 = time.time()
    import multiprocessing as mp
    with mp.Pool(min(a.workers, len(todo)), maxtasksperchild=1) as pool:
        res = pool.map(build_epoch, todo)
    Mp = f"{FIX}/templates_new.csv"; Zp = f"{FIX}/templates_new.npz"
    M = pd.read_csv(Mp, dtype={"fam": str}) if os.path.exists(Mp) else pd.DataFrame()
    Z = dict(np.load(Zp)) if os.path.exists(Zp) else {}
    new_rows = []
    for rows, tpl in res:
        new_rows += rows; Z.update(tpl)
    M = pd.concat([M, pd.DataFrame(new_rows)], ignore_index=True)
    M.to_csv(Mp, index=False)
    np.savez(Zp + ".tmp.npz", **Z); os.replace(Zp + ".tmp.npz", Zp)
    stn = pd.DataFrame([dict(epoch=t[0]["epoch"], state="done", final_month=a.upto) for t in todo])
    pd.concat([st, stn], ignore_index=True).to_csv(STATE, index=False)
    with open(f"{FIX}/rerun_months.txt", "a") as fh:
        for t, (rows, tpl) in zip(todo, res):
            if tpl:
                for m in pd.period_range(pd.Timestamp(t[0]["start"]).to_period("M"), pd.Period(a.upto, "M"), freq="M"):
                    fh.write(f"{m}\n")
    R = pd.DataFrame(new_rows)
    log(f"[retemplate {a.upto}] finalised {len(todo)} epochs: {int(R.keep.sum())}/{len(R)} (family,channel) "
        f"templates kept, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
