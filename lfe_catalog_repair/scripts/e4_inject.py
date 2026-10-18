#!/usr/bin/env python
"""E4 step 3 — corrected injection + a RE-ANCHORED real recall curve (Merlin E3-audit step 2).

Three fixes to E3 step 3:
 1. components. e3_inject.py inserted the stack as Z,E,N into Z,N,E data (E4 step 0). Fixed.
 2. the x-axis anchor. E3 put the real events on the injected curve's axis using
    s/n = 0.078 x 10^(1.5 dM), where 0.078 is rms(stack)/rms(null stack) on a ~130-event stack --
    a COHERENT amplitude, which under-states a single event's incoherent s/n. Here the per-event
    INCOHERENT s/n is measured DIRECTLY from the data with no stacking at all:
        s/n^2 = mean over the 3 components of [ rms^2(S window) / rms^2(the same window 20 s
                earlier) ] - 1,     1-8 Hz, S window = S-1 .. S+5 s,
    for every Bostock event with epi_PGC < 40 km, and the per-magnitude-bin median is the anchor.
 3. the injected waveform. Besides the published SeisBench stack, OWN N-event stacks are built at
    PGC from the Bostock times (N = 5 / 20 / 100 / all) and injected at matched snr_rms, so the
    "a stack is a pessimistic proxy for a single event" explanation is tested rather than asserted.

Outputs results/e4_inject_snr_incoherent.csv, e4_inject_recall.csv, e4_inject_vs_real.csv
"""
from __future__ import annotations
import argparse, os, sys, time
import numpy as np, pandas as pd, h5py
from scipy.signal import butter, sosfiltfilt, resample_poly

ROOT = "/home/jovyan/tremorferometry/lfe_catalog_repair"
RES = f"{ROOT}/results"
CACHE = f"{ROOT}/cache20hz"
SD = f"{RES}/external/other/seisbench_lfestacks_cascadia"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e2_pick as EP          # noqa: E402
import e3_inject as EI        # noqa: E402
import e4_chanorder as chanorder  # noqa: E402
import d3_vel1d as v1d        # noqa: E402

FS = 20.0
MB = [(-9, 1.5, "M<1.5"), (1.5, 1.58, "M1.5-1.58"), (1.58, 1.66, "M1.58-1.66"),
      (1.66, 1.74, "M1.66-1.74"), (1.74, 1.82, "M1.74-1.82"), (1.82, 1.95, "M1.82-1.95"),
      (1.95, 2.1, "M1.95-2.1"), (2.1, 9, "M>=2.1")]
QUIET = ["2010-05-20", "2010-05-02", "2010-06-03", "2010-05-30", "2010-06-13", "2010-05-01"]
NSTACK = [5, 20, 100, 0]          # 0 = all available
LEVELS = (0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0)


def bp(x, fs, f1=1.0, f2=8.0):
    sos = butter(4, [f1 / (fs / 2), f2 / (fs / 2)], btype="band", output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def bostock_pgc():
    bo = pd.read_parquet(f"{RES}/e2_gt_bostock.parquet")
    fo = pd.read_csv(f"{RES}/e3_bostock_family_offsets.csv")
    acc = fo[fo.accept.astype(bool)]
    om = dict(zip(acc.bkey.astype(int), acc.offset_s.astype(float)))
    bo = bo[bo.bkey.astype(int).isin(om)].copy()
    bo["time"] = bo.time.astype("datetime64[ns]")
    bo["t_origin"] = bo.time.astype("int64") / 1e9 + bo.bkey.astype(int).map(om).to_numpy()
    bo["t_S_PGC"] = bo.t_origin + bo.ttS_PGC
    return bo[bo.epi_PGC < 40].copy()


def incoherent_snr(bo):
    """Per-event incoherent s/n at PGC from the cached 20 Hz, 1-8 Hz data. No stacking."""
    n = int(round(6 * FS))
    rows = []
    cache = {}
    for r in bo.itertuples():
        d = pd.Timestamp(r.t_S_PGC, unit="s")
        key = d.dayofyear
        if key not in cache:
            p = f"{CACHE}/CN.PGC/2010/{key:03d}.npz"
            cache[key] = (float(np.load(p)["t0"]), np.load(p)["x"]) if os.path.exists(p) else None
        v = cache[key]
        if v is None:
            continue
        t0, X = v
        i = int(round((r.t_S_PGC - 1.0 - t0) * FS))
        j = i - int(round(20 * FS))
        if j < 0 or i + n > X.shape[1]:
            continue
        sg = X[:, i:i + n].astype(np.float64)
        nz = X[:, j:j + n].astype(np.float64)
        ps, pn = (sg ** 2).mean(), (nz ** 2).mean()
        if pn <= 0:
            continue
        rows.append(dict(t=r.t_S_PGC, mag=r.mag, epi=r.epi_PGC, fam=r.fam_id,
                         ratio=ps / pn, snr_inc=np.sqrt(max(ps / pn - 1.0, 0.0))))
    # CONTROL, added before any number is quoted: the same statistic at (a) times drawn uniformly
    # on the same days and (b) the Bostock times shifted by a random 300-3000 s.  Bostock events sit
    # inside tremor bursts, so a raw "S window vs 20 s earlier" ratio can be elevated by the burst
    # envelope and by neighbouring LFEs rather than by the event itself.
    rng = np.random.default_rng(21)
    ctrl = {}
    for name, tt in (("random", rng.uniform(bo.t_S_PGC.min(), bo.t_S_PGC.max(), 4000)),
                     ("shifted", bo.t_S_PGC.to_numpy()
                      + rng.uniform(300, 3000, len(bo)) * rng.choice([-1, 1], len(bo)))):
        rr = []
        for t in tt:
            d = pd.Timestamp(t, unit="s")
            key = d.dayofyear
            if key not in cache:
                p_ = f"{CACHE}/CN.PGC/2010/{key:03d}.npz"
                cache[key] = ((float(np.load(p_)["t0"]), np.load(p_)["x"])
                              if os.path.exists(p_) else None)
            v = cache[key]
            if v is None:
                continue
            t0, X = v
            i = int(round((t - 1.0 - t0) * FS)); j = i - int(round(20 * FS))
            if j < 0 or i + n > X.shape[1]:
                continue
            ps = (X[:, i:i + n].astype(np.float64) ** 2).mean()
            pn = (X[:, j:j + n].astype(np.float64) ** 2).mean()
            if pn > 0:
                rr.append(ps / pn)
        ctrl[name] = np.array(rr)
        print(f"[ctrl] {name}: n={len(rr)} power-ratio mean {np.mean(rr):.4f} "
              f"median {np.median(rr):.4f} trimmed-mean(5-95%) "
              f"{np.mean(np.clip(rr, *np.percentile(rr, [5, 95]))):.4f}", flush=True)

    D = pd.DataFrame(rows)
    out = []
    for lo, hi, lab in MB:
        h = D[(D.mag > lo) & (D.mag <= hi)]
        if len(h) < 20:
            continue
        # the ESTIMATOR: the bin's mean power ratio minus 1, square-rooted (a per-event median of
        # sqrt(ratio-1) is biased low because half the events have ratio < 1 by noise alone)
        # trimmed mean (5-95 %) rather than the raw mean: one glitchy window puts the M>=2.1 bin
        # at ratio 11 with n=60.  The control's trimmed mean is subtracted as the baseline.
        lo5, hi95 = np.percentile(h.ratio, [5, 95])
        r_mean = float(np.mean(np.clip(h.ratio, lo5, hi95)))
        c = ctrl["random"]
        c_mean = float(np.mean(np.clip(c, *np.percentile(c, [5, 95]))))
        out.append(dict(magbin=lab, m_lo=lo, m_hi=hi, n=len(h),
                        mag_median=float(h.mag.median()),
                        ratio_mean=r_mean, ratio_ctrl=c_mean,
                        snr_inc=float(np.sqrt(max(r_mean - c_mean, 0))),
                        snr_inc_vs1=float(np.sqrt(max(r_mean - 1, 0))),
                        ratio_median=float(h.ratio.median()),
                        snr_inc_median=float(np.sqrt(max(h.ratio.median() - 1, 0))),
                        snr_inc_sem=float(h.ratio.std() / np.sqrt(len(h)) / 2 /
                                          max(np.sqrt(max(r_mean - 1, 1e-6)), 1e-6))))
    O = pd.DataFrame(out)
    O.to_csv(f"{RES}/e4_inject_snr_incoherent.csv", index=False)
    D.to_parquet(f"{RES}/e4_inject_snr_events.parquet", index=False)
    print("\n=== MEASURED per-event INCOHERENT s/n at PGC (epi < 40 km, no stacking) ===")
    print(O.round(4).to_string(index=False))
    return O, D


def own_stacks(bo, nlist):
    """Own N-event stacks at PGC, per family, aligned on the predicted S. -> {(N, fam): (3, nt)}"""
    n = int(round(30 * FS))
    pre = int(round(10 * FS))
    cache = {}
    per = {}
    for r in bo.itertuples():
        d = pd.Timestamp(r.t_S_PGC, unit="s")
        key = d.dayofyear
        if key not in cache:
            p = f"{CACHE}/CN.PGC/2010/{key:03d}.npz"
            cache[key] = (float(np.load(p)["t0"]), np.load(p)["x"]) if os.path.exists(p) else None
        v = cache[key]
        if v is None:
            continue
        t0, X = v
        i = int(round((r.t_S_PGC - t0) * FS)) - pre
        if i < 0 or i + n > X.shape[1]:
            continue
        sg = X[:, i:i + n].astype(np.float64)
        rr = np.sqrt((sg ** 2).mean())
        if rr <= 0:
            continue
        per.setdefault(r.fam_id, []).append(sg / rr)
    out = {}
    rng = np.random.default_rng(5)
    for fam, lst in per.items():
        if len(lst) < 20:
            continue
        A = np.stack(lst)
        for N in nlist:
            k = len(A) if N == 0 else N
            if k > len(A):
                continue
            idx = rng.choice(len(A), k, replace=False)
            out[(N, fam)] = (A[idx].mean(0), pre)
    print(f"[own] {len(out)} own-stack templates over {len(per)} families "
          f"(N in {nlist}; N=0 means all, median all = "
          f"{int(np.median([len(v) for v in per.values()]))})", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-day", type=int, default=140)
    ap.add_argument("--days", type=int, default=6)
    ap.add_argument("--model", default="cascadia")
    ap.add_argument("--torch-threads", type=int, default=4)
    ap.add_argument("--skip-pick", action="store_true")
    ap.add_argument("--nstack", default="20,0")
    a = ap.parse_args()

    nst_pre = None
    bo = bostock_pgc()
    print(f"[evt] {len(bo)} Bostock events with epi_PGC < 40 km", flush=True)
    SNR, _ = incoherent_snr(bo)
    nst = [int(x) for x in a.nstack.split(",")]
    own = own_stacks(bo, nst)
    sb = EI.load_stacks("PGC")          # the published SeisBench stacks at PGC
    print(f"[sb] {len(sb)} SeisBench PGC stacks", flush=True)
    if a.skip_pick:
        return

    import torch
    torch.set_num_threads(a.torch_threads)
    from obspy import read
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    mdl = EP.load_model(a.model, dev)
    print(f"[env] {dev}; model {a.model}", flush=True)

    rng = np.random.default_rng(9)
    rows, plan_rows = [], []
    for dstr in QUIET[:a.days]:
        d = pd.Timestamp(dstr)
        p = f"{ROOT}/traces/CN.PGC/{d.year}/{d.dayofyear:03d}.mseed"
        if not os.path.exists(p):
            continue
        st = read(p); st.merge(method=1, fill_value=0)
        chans = sorted({tr.stats.channel for tr in st})
        order = chanorder.resolve_order(chans, station="CN.PGC")
        fs = float(st[0].stats.sampling_rate)
        npts = min(tr.stats.npts for tr in st.select(channel=order[0]))
        X = np.vstack([st.select(channel=c)[0].data[:npts].astype(float) for c in order])
        Xf = bp(X, fs)
        t0 = st.select(channel=order[0])[0].stats.starttime.timestamp
        print(f"[day {dstr}] fs={fs:.0f}, order {order}", flush=True)

        kinds = [("seisbench", None)] + [("own", N) for N in nst]
        for kind, N in kinds:
            n_inj = a.per_day
            slots = np.arange(n_inj) * (npts / fs / n_inj) + 300.0 + rng.uniform(-60, 60, n_inj)
            lev = np.array([LEVELS[i % len(LEVELS)] for i in range(n_inj)])
            rng.shuffle(lev)
            Xi = X.copy()
            plan = []
            up = int(round(fs / FS))
            for ts, L in zip(slots, lev):
                if kind == "seisbench":
                    w, s_samp, fid = sb[rng.integers(len(sb))]
                    wf = bp(w, 20.0)
                else:
                    keys = [k for k in own if k[0] == N]
                    if not keys:
                        continue
                    k = keys[rng.integers(len(keys))]
                    wf, s_samp = own[k]
                    fid = k[1]
                wr = resample_poly(wf, up, 1, axis=-1) if up > 1 else wf
                s_r = s_samp * up
                i_s = int(round(ts * fs))
                aa, bb = max(i_s - int(15 * fs), 0), min(i_s + int(15 * fs), npts)
                rms_host = float(np.sqrt((Xf[:, aa:bb] ** 2).mean()))
                ra, rb = int(s_r - 1 * fs), int(s_r + 5 * fs)
                ra, rb = max(ra, 0), min(rb, wr.shape[1])
                rms_sig = float(np.sqrt((wr[:, ra:rb] ** 2).mean()))
                if rms_host <= 0 or rms_sig <= 0:
                    continue
                g = L * rms_host / rms_sig          # scale so that snr_rms == L exactly
                wr = wr * g
                i0 = i_s - int(round(s_r))
                j0, j1 = max(i0, 0), min(i0 + wr.shape[1], npts)
                if j1 - j0 < wr.shape[1] * 0.9:
                    continue
                Xi[:, j0:j1] += wr[:, j0 - i0:j1 - i0]
                plan.append(dict(day=dstr, kind=kind, nstack=(-1 if N is None else N),
                                 t_S=t0 + ts, snr_rms=float(L), fam=fid))
            plan = pd.DataFrame(plan)
            plan_rows.append(plan)
            picks = {}
            for tagname, DATA in (("base", X), ("inj", Xi)):
                if tagname == "base" and (kind, N) != kinds[0]:
                    picks["base"] = picks.get("base")
                    if picks["base"] is not None:
                        continue
                stw = st.copy()
                for c, row in zip(order, DATA):
                    tr = stw.select(channel=c)[0]
                    tr.data = np.ascontiguousarray(row[:npts], dtype=np.float32)
                    tr.stats.npts = npts
                pk = pd.DataFrame(EP.picks_for_day(mdl, stw, 1024, 600))
                pk = pk[pk.phase == "S"] if len(pk) else pk
                if len(pk):
                    picks[tagname] = (pd.to_datetime(pk.time).astype("datetime64[ns]")
                                      .astype("int64").to_numpy() / 1e9,
                                      pk.thr_first.to_numpy(), pk.prob.to_numpy())
                else:
                    picks[tagname] = (np.zeros(0), np.zeros(0), np.zeros(0))
            for thr in (0.1, 0.2, 0.3, 0.5):
                for lo in sorted(set(plan.snr_rms)):
                    h = plan[plan.snr_rms == lo]
                    fr = {}
                    for tagname in ("inj", "base"):
                        tt, th, pr = picks[tagname]
                        m = (th <= thr) & (pr >= thr)
                        tt2 = np.sort(tt[m])
                        if len(tt2) == 0:
                            fr[tagname] = 0.0
                            continue
                        i = np.searchsorted(tt2, h.t_S.to_numpy())
                        dm = np.full(len(h), np.inf)
                        for k in (-1, 0):
                            j = np.clip(i + k, 0, len(tt2) - 1)
                            dm = np.minimum(dm, np.abs(tt2[j] - h.t_S.to_numpy()))
                        fr[tagname] = float((dm <= 1.0).mean())
                    rows.append(dict(day=dstr, kind=kind, nstack=(-1 if N is None else N),
                                     threshold=thr, snr_rms=lo, n=len(h),
                                     frac_inj=fr["inj"], frac_base=fr["base"],
                                     recall=max(fr["inj"] - fr["base"], 0.0)))
            print(f"  [{kind} N={N}] {len(plan)} injections", flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(f"{RES}/e4_inject_recall_raw.csv", index=False)
    G = R.groupby(["kind", "nstack", "threshold", "snr_rms"]).agg(
        recall=("recall", "mean"), frac_inj=("frac_inj", "mean"),
        frac_base=("frac_base", "mean"), n=("n", "sum")).reset_index()
    G.to_csv(f"{RES}/e4_inject_recall.csv", index=False)
    pd.concat(plan_rows, ignore_index=True).to_parquet(f"{RES}/e4_inject_plan.parquet", index=False)
    print("\n=== injected recall vs snr_rms (thr 0.1) ===")
    print(G[G.threshold == 0.1].round(4).to_string(index=False))
    print(f"[write] {RES}/e4_inject_recall.csv")


if __name__ == "__main__":
    main()
