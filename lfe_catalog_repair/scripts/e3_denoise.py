#!/usr/bin/env python
"""E3 step 3b (user request 2026-09-11) — does DENOISING before picking help?

Three arms, everything else held identical (same days, same injections, same band-pass -- the one
LFEDetect applies itself -- and the SAME sampling rate in every arm, 100 Hz, which is
DeepDenoiser's native rate; the 40 Hz hosts PGC/LZB are resampled once, up front, for all arms, so
no arm sees a different rate):
  raw     the 100 Hz stream, straight into LFEDetect                         (the control)
  dn      DeepDenoiser `original` (Zhu et al. 2019, MIT, weights already in ~/.seisbench -- no
          download was needed) applied per component, then LFEDetect
  wiener  a non-ML baseline: band-limited spectral Wiener / soft-masking with the noise PSD taken
          from the quietest decile of the day's own 30 s windows, then LFEDetect

Measurements
 (i)   recall(SNR) on the SAME injections as step 3 (`e3_inject_plan.parquet` replayed: same
       station, day, S time, family stack and level; the gain is recomputed from the 100 Hz host
       RMS, which is within a per-cent of the 40 Hz one), with the un-injected baseline of the SAME
       arm subtracted, plus pick-time bias and scatter;
 (ii)  recall against Lin-hq and Bostock on REAL dense-ground-truth PGC 2010 days, by the +-1 s
       peak-excess method of `e2_resid_hist.py` (the only chance-free estimator at these pick
       rates);
 (iii) the false-alarm rate F per station-day per threshold, per arm;
 (iv)  WHAT THE DENOISER DOES TO THE LFE ITSELF: with D = arm(host + signal) - arm(host) the
       recovered signal is exact, so rms(D)/rms(signal) in the S window is the amplitude retained,
       corr(D, signal) the waveform fidelity, and |FFT(D)|/|FFT(signal)| in 1 Hz bands the
       spectral shape -- this is what says whether the denoiser strips emergent LFE energy.

Outputs results/e3_denoise_{recall,falsealarm,realdays,fidelity}.csv
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
TRACES = "lfe_catalog_repair/traces"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e2_pick as EP        # noqa: E402
import e2_resid_hist as RH  # noqa: E402
import e3_inject as EI      # noqa: E402
import e4_chanorder as chanorder  # noqa: E402

LEGACY_CHAN = bool(int(os.environ.get("E3_LEGACY_CHAN_ORDER", "0")))  # 1 = reproduce the pre-fix bug
_ORDER_SEEN = {}

FS = 100.0
SWEEP = EP.SWEEP


def to_fs(X, fs_in, fs_out=FS):
    from scipy.signal import resample_poly
    if fs_in == fs_out:
        return X
    from math import gcd
    g = gcd(int(fs_out), int(fs_in))
    return resample_poly(X, int(fs_out) // g, int(fs_in) // g, axis=-1)


def wiener_denoise(X, fs=FS, nper=3000, f1=1.0, f2=8.0, quiet_q=0.10):
    """Band-limited spectral soft-mask: noise PSD = the quietest `quiet_q` quantile of the day's
    own window PSDs; mask = S/(S+N) per time-frequency cell (a Wiener gain), applied in 1-8 Hz and
    passing everything else through unchanged so the picker sees the same band structure."""
    from scipy.signal import stft, istft
    out = np.empty_like(X)
    for c in range(X.shape[0]):
        f, t, Z = stft(X[c], fs=fs, nperseg=nper, noverlap=nper // 2)
        P = np.abs(Z) ** 2
        band = (f >= f1) & (f <= f2)
        Nf = np.quantile(P[:, :], quiet_q, axis=1, keepdims=True)     # per-frequency noise floor
        G = np.ones_like(P)
        num = np.clip(P[band] - Nf[band], 0.0, None)
        G[band] = num / np.clip(P[band], 1e-30, None)
        _, y = istft(Z * G, fs=fs, nperseg=nper, noverlap=nper // 2)
        out[c] = y[:X.shape[1]] if len(y) >= X.shape[1] else np.pad(y, (0, X.shape[1] - len(y)))
    return out


def deepdenoise(X, dnmodel, fs=FS, starttime=None, chans=None):
    """DeepDenoiser per component, via SeisBench annotate on an obspy Stream."""
    from obspy import Stream, Trace, UTCDateTime
    st = Stream()
    for c, row in zip(chans, X):
        tr = Trace(np.ascontiguousarray(row, dtype=np.float32))
        tr.stats.sampling_rate = fs
        tr.stats.starttime = UTCDateTime(starttime)
        tr.stats.network, tr.stats.station, tr.stats.channel = "XX", "DN", c
        st += tr
    ann = dnmodel.annotate(st)
    out = np.zeros_like(X)
    for i, c in enumerate(chans):
        sel = [tr for tr in ann if tr.stats.channel.endswith(c[-1])] or \
              [tr for tr in ann if c[-1] in tr.stats.channel]
        if not sel:
            sel = [ann[i]] if len(ann) > i else []
        if sel:
            d = sel[0].data.astype(float)
            n = min(len(d), X.shape[1])
            out[i, :n] = d[:n]
    return out


def stream_from(X, chans, fs, starttime, net, sta):
    from obspy import Stream, Trace, UTCDateTime
    st = Stream()
    for c, row in zip(chans, X):
        tr = Trace(np.ascontiguousarray(row, dtype=np.float32))
        tr.stats.sampling_rate = fs
        tr.stats.starttime = UTCDateTime(starttime)
        tr.stats.network, tr.stats.station, tr.stats.channel = net, sta, c
        st += tr
    return st


def s_picks(mdl, st):
    """S pick times with their thr_first AND their peak prob.  A pick exists at threshold t iff
    thr_first <= t <= prob (e2_pick docstring); filtering on thr_first alone -- the first version
    of this script, fixed 2026-09-11 -- makes every threshold identical to the 0.1 pick set."""
    df = pd.DataFrame(EP.picks_for_day(mdl, st, 1024, 600))
    if len(df) == 0:
        return np.zeros(0), np.zeros(0), np.zeros(0)
    df = df[df.phase == "S"]
    t = (pd.to_datetime(df.time).astype("datetime64[ns]").astype("int64").to_numpy() / 1e9)
    o = np.argsort(t)
    return t[o], df.thr_first.to_numpy()[o], df.prob.to_numpy()[o]


def at_thr(t, th, pr, thr):
    return np.sort(t[(th <= thr) & (pr >= thr)])


def nearest_dt(q, s):
    if len(s) == 0:
        return np.full(len(q), np.inf)
    i = np.searchsorted(s, q)
    best = np.full(len(q), np.inf)
    for off in (-1, 0):
        j = np.clip(i + off, 0, len(s) - 1)
        d = s[j] - q
        m = np.abs(d) < np.abs(best)
        best = np.where(m, d, best)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", default="raw,dn,wiener")
    ap.add_argument("--models", default="cascadia,jcms")
    ap.add_argument("--stations", default="PGC,LZB,YOUB")
    ap.add_argument("--real-days", type=int, default=4,
                    help="how many dense-ground-truth PGC 2010 days for part (ii)")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--torch-threads", type=int, default=4)
    args = ap.parse_args()

    import torch
    torch.set_num_threads(args.torch_threads)
    from obspy import read
    import seisbench.models as sbm
    dev = args.device if (args.device != "cuda" or torch.cuda.is_available()) else "cpu"
    print(f"[env] {dev} {torch.cuda.get_device_name(0) if dev=='cuda' else ''}", flush=True)
    pickers = {mn: EP.load_model(mn, dev) for mn in args.models.split(",")}
    arms = args.arms.split(",")
    dn = None
    if "dn" in arms:
        dn = sbm.DeepDenoiser.from_pretrained("original")
        dn.to(dev); dn.eval()
        print(f"[dn] DeepDenoiser original, {dn.sampling_rate} Hz, "
              f"{dn.in_samples} samples/window", flush=True)

    plan = pd.read_parquet(f"{RES}/e3_inject_plan.parquet")
    rec_rows, fa_rows, fid_rows, real_rows = [], [], [], []

    # ================= part (i), (iii), (iv): the injected quiet days =================
    for stn in args.stations.split(","):
        netsta, _ = EI.HOSTS[stn]
        stacks = {fid: (w, s) for w, s, fid in EI.load_stacks(stn)}
        for dstr, g in plan[plan.station == stn].groupby("day"):
            d = pd.Timestamp(dstr)
            p = f"{TRACES}/{netsta}/{d.year}/{d.dayofyear:03d}.mseed"
            if not os.path.exists(p):
                continue
            st0 = read(p); st0.merge(method=1, fill_value=0)
            fs0 = float(st0[0].stats.sampling_rate)
            t0 = st0[0].stats.starttime.timestamp
            npts0 = min(tr.stats.npts for tr in st0)
            chans = sorted({tr.stats.channel for tr in st0})
            # FIX 2026-09-11 (E4 step 0): was [Z]+sorted(non-Z) = Z,E,N on CN stations while the
            # injected stack is Z,N,E (1=N, e4_polarization_check.py: 7.9 vs 37.3 deg).
            order = chanorder.resolve_order(chans, station=f"{stn}/{d.date()}", legacy=LEGACY_CHAN)
            if stn not in _ORDER_SEEN:
                print(f"  [chanorder] {stn}: {order}", flush=True); _ORDER_SEEN[stn] = order
            X0 = np.vstack([st0.select(channel=c)[0].data[:npts0].astype(float) for c in order])
            X = to_fs(X0, fs0, FS)
            npts = X.shape[1]
            Xf = EI.bandpass(X, FS)
            # replay the SAME injections at 100 Hz
            S = np.zeros_like(X)
            keep = []
            for r in g.itertuples():
                if r.fam_id not in stacks:
                    continue
                w, s_samp = stacks[r.fam_id]
                wf = EI.bandpass(w, 20.0)
                wr = to_fs(wf, 20.0, FS)
                s_r = s_samp * (FS / 20.0)
                i_s = int(round((r.t_S - t0) * FS))
                a, b = max(i_s - int(15 * FS), 0), min(i_s + int(15 * FS), npts)
                if b <= a:
                    continue
                rms_host = float(np.sqrt((Xf[:, a:b] ** 2).mean()))
                sa, sb = max(int(s_r - 2 * FS), 0), min(int(s_r + 8 * FS), wr.shape[1])
                pk = float(np.abs(wr[:, sa:sb]).max())
                if pk <= 0 or rms_host <= 0 or not np.isfinite(rms_host):
                    continue
                wr = wr * (r.level * rms_host / pk)
                i0 = i_s - int(round(s_r))
                j0, j1 = max(i0, 0), min(i0 + wr.shape[1], npts)
                if j1 - j0 < wr.shape[1] * 0.9:
                    continue
                S[:, j0:j1] += wr[:, j0 - i0:j1 - i0]
                keep.append(dict(t_S=r.t_S, level=r.level, snr_rms=r.snr_rms, fam_id=r.fam_id,
                                 i_s=i_s, s_r=int(round(s_r)), rms_host=rms_host))
            K = pd.DataFrame(keep)
            Xi = X + S
            print(f"  [{stn} {dstr}] {fs0:.0f}->{FS:.0f} Hz, {len(K)} injections replayed",
                  flush=True)

            for arm in arms:
                t1 = time.time()
                if arm == "raw":
                    Ab, Ai = X, Xi
                elif arm == "wiener":
                    Ab, Ai = wiener_denoise(X), wiener_denoise(Xi)
                else:
                    Ab = deepdenoise(X, dn, FS, t0, order)
                    Ai = deepdenoise(Xi, dn, FS, t0, order)
                # (iv) fidelity: D = arm(host+sig) - arm(host) vs the true injected signal
                D = Ai - Ab
                for r in K.itertuples():
                    a, b = max(r.i_s - int(1 * FS), 0), min(r.i_s + int(5 * FS), npts)
                    sig = S[:, a:b]; rec = D[:, a:b]
                    rs = float(np.sqrt((sig ** 2).mean()))
                    rr = float(np.sqrt((rec ** 2).mean()))
                    cc = float(np.corrcoef(sig.ravel(), rec.ravel())[0, 1]) if rs > 0 else np.nan
                    row = dict(station=stn, day=dstr, arm=arm, level=r.level,
                               amp_ratio=rr / rs if rs > 0 else np.nan, corr=cc)
                    if rs > 0:
                        fs_sig = np.fft.rfft(sig, axis=-1); fr = np.fft.rfftfreq(b - a, 1 / FS)
                        fs_rec = np.fft.rfft(rec, axis=-1)
                        for lo, hi in ((1, 2), (2, 4), (4, 6), (6, 8)):
                            mb = (fr >= lo) & (fr < hi)
                            den = np.abs(fs_sig[:, mb]).mean()
                            row[f"specratio_{lo}_{hi}"] = float(
                                np.abs(fs_rec[:, mb]).mean() / den) if den > 0 else np.nan
                    fid_rows.append(row)
                for mn, mdl in pickers.items():
                    tb, hb, pb = s_picks(mdl, stream_from(Ab, order, FS, t0, "XX", stn))
                    ti, hi_, pi = s_picks(mdl, stream_from(Ai, order, FS, t0, "XX", stn))
                    for thr in SWEEP:
                        sb, si = at_thr(tb, hb, pb, thr), at_thr(ti, hi_, pi, thr)
                        fa_rows.append(dict(station=stn, day=dstr, arm=arm, model=mn,
                                            threshold=thr, n_S_base=len(sb), n_S_inj=len(si),
                                            day_seconds=npts / FS))
                        for tol in (1.0, 3.0):
                            db = nearest_dt(K.t_S.to_numpy(), sb)
                            di = nearest_dt(K.t_S.to_numpy(), si)
                            for lev, gg in K.groupby("level"):
                                idx = K.level.to_numpy() == lev
                                rec_rows.append(dict(
                                    station=stn, day=dstr, arm=arm, model=mn, threshold=thr,
                                    tol=tol, level=float(lev), n=int(idx.sum()),
                                    frac_inj=float((np.abs(di[idx]) <= tol).mean()),
                                    frac_base=float((np.abs(db[idx]) <= tol).mean()),
                                    snr_rms=float(gg.snr_rms.median()),
                                    bias_s=float(np.nanmedian(di[idx][np.abs(di[idx]) <= tol]))
                                    if (np.abs(di[idx]) <= tol).sum() > 5 else np.nan,
                                    scatter_s=float((np.nanpercentile(
                                        di[idx][np.abs(di[idx]) <= tol], 84)
                                        - np.nanpercentile(di[idx][np.abs(di[idx]) <= tol], 16)) / 2)
                                    if (np.abs(di[idx]) <= tol).sum() > 10 else np.nan))
                print(f"    arm={arm} done in {time.time()-t1:.0f} s", flush=True)

    R = pd.DataFrame(rec_rows)
    R["recall"] = (R.frac_inj - R.frac_base) / (1 - R.frac_base).clip(lower=1e-9)
    R.to_csv(f"{RES}/e3_denoise_recall.csv", index=False)
    pd.DataFrame(fa_rows).to_csv(f"{RES}/e3_denoise_falsealarm.csv", index=False)
    pd.DataFrame(fid_rows).to_csv(f"{RES}/e3_denoise_fidelity.csv", index=False)

    # ================= part (ii): REAL dense-ground-truth PGC 2010 days =================
    bo = pd.read_parquet(f"{RES}/e2_gt_bostock.parquet")
    lin = pd.read_parquet(f"{RES}/e2_gt_lin.parquet")
    linq = lin[(lin.residual < 0.5) & (lin.N >= 4)].copy()
    fo = pd.read_csv(f"{RES}/e3_bostock_family_offsets.csv")
    accm = dict(zip(fo[fo.accept.astype(bool)].bkey.astype(int),
                    fo[fo.accept.astype(bool)].offset_s))
    medoff = float(fo[fo.accept.astype(bool)].offset_s.median())
    bo["t_pred"] = (bo.time.astype("datetime64[ns]").astype("int64") / 1e9
                    + bo.bkey.astype(int).map(accm).fillna(medoff).to_numpy() + bo.ttS_PGC)
    linq["t_pred"] = (linq.OT.astype("datetime64[ns]").astype("int64") / 1e9 + linq.ttS_PGC)
    cnt = linq.OT.dt.normalize().value_counts()
    days = [d for d in cnt.index
            if os.path.exists(f"{TRACES}/CN.PGC/{d.year}/{d.dayofyear:03d}.mseed")][:args.real_days]
    print(f"[real] PGC ground-truth days: "
          f"{[str(d.date()) for d in days]} (Lin-hq {[int(cnt[d]) for d in days]})", flush=True)
    for d in days:
        p = f"{TRACES}/CN.PGC/{d.year}/{d.dayofyear:03d}.mseed"
        st0 = read(p); st0.merge(method=1, fill_value=0)
        fs0 = float(st0[0].stats.sampling_rate); t0 = st0[0].stats.starttime.timestamp
        npts0 = min(tr.stats.npts for tr in st0)
        chans = sorted({tr.stats.channel for tr in st0})
        # FIX 2026-09-11 (E4 step 0): see the note above -- Z,N,E not Z,E,N.
        order = chanorder.resolve_order(chans, station="CN.PGC", legacy=LEGACY_CHAN)
        if "realday" not in _ORDER_SEEN:
            print(f"  [chanorder] CN.PGC real days: {order}", flush=True); _ORDER_SEEN["realday"] = order
        X0 = np.vstack([st0.select(channel=c)[0].data[:npts0].astype(float) for c in order])
        X = to_fs(X0, fs0, FS); npts = X.shape[1]
        t1e = t0 + npts / FS
        gts = {"lin_hq": linq.t_pred.to_numpy()[(linq.t_pred >= t0) & (linq.t_pred < t1e)],
               "bostock": bo.t_pred.to_numpy()[(bo.t_pred >= t0) & (bo.t_pred < t1e)]}
        for arm in arms:
            A = X if arm == "raw" else (wiener_denoise(X) if arm == "wiener"
                                        else deepdenoise(X, dn, FS, t0, order))
            for mn, mdl in pickers.items():
                tt, hh, pp = s_picks(mdl, stream_from(A, order, FS, t0, "CN", "PGC"))
                for thr in SWEEP:
                    s = at_thr(tt, hh, pp, thr)
                    for gk, gt in gts.items():
                        if len(gt) < 20 or len(s) == 0:
                            continue
                        dd = RH.pairs_within(np.sort(gt), np.sort(s), 30.0)
                        e = np.arange(-30, 30.001, 0.25); h, _ = np.histogram(dd, bins=e)
                        c = 0.5 * (e[1:] + e[:-1])
                        far = np.abs(c) > 10
                        bg = float(np.median(h[far]))
                        m1 = np.abs(c) <= 1.0
                        exc = float(h[m1].sum() - bg * m1.sum())
                        real_rows.append(dict(day=str(d.date()), arm=arm, model=mn, threshold=thr,
                                              gt=gk, n_gt=len(gt), n_picks=int(len(s)),
                                              excess=exc, recall_excess=exc / len(gt),
                                              precision_excess=exc / max(len(s), 1),
                                              sigma=exc / np.sqrt(max(h[m1].sum(), 1))))
        print(f"  [real {d.date()}] done", flush=True)
    pd.DataFrame(real_rows).to_csv(f"{RES}/e3_denoise_realdays.csv", index=False)

    # ---- printed summary
    print("\n[3b (i) recall, tol 1 s, thr 0.3, cascadia, pooled over days]")
    q = R[(R.tol == 1.0) & (R.threshold == 0.3) & (R.model == "cascadia")]
    print(q.groupby(["station", "arm", "level"]).agg(
        n=("n", "sum"), recall=("recall", "mean"), base=("frac_base", "mean"),
        bias=("bias_s", "mean")).round(3).to_string())
    print("\n[3b (iii) false alarms, S picks/station-day]")
    F = pd.DataFrame(fa_rows)
    print(F[F.threshold.isin([0.1, 0.3, 0.5])].groupby(
        ["station", "arm", "model", "threshold"]).n_S_base.mean().round(1).to_string())
    print("\n[3b (iv) what the denoiser does to the injected LFE]")
    Fi = pd.DataFrame(fid_rows)
    print(Fi.groupby(["arm", "level"])[
        [c for c in Fi.columns if c.startswith(("amp_ratio", "corr", "specratio"))]]
        .mean().round(3).to_string())
    print("\n[3b (ii) real-day peak excess]")
    RD = pd.DataFrame(real_rows)
    if len(RD):
        print(RD[RD.threshold.isin([0.1, 0.3])].groupby(
            ["gt", "arm", "model", "threshold"]).agg(
            excess=("excess", "sum"), n_gt=("n_gt", "sum"), n_picks=("n_picks", "sum"))
            .assign(recall=lambda t: t.excess / t.n_gt,
                    precision=lambda t: t.excess / t.n_picks).round(4).to_string())
    print(f"\n[write] {RES}/e3_denoise_{{recall,falsealarm,fidelity,realdays}}.csv")


if __name__ == "__main__":
    main()
