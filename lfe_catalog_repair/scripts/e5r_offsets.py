#!/usr/bin/env python
"""E5 step 5 STEP B(ii) -- the per-family Bostock TIME-BASE OFFSET re-derived on the 2012 BUILD half at PGC.

Two estimators, both at PGC, both on the build half of the provisional split (alternate 2-h blocks of DB time + 2010 offset,
60-s guard; the split is rebuilt on the final offsets afterwards by e5r_groundtruth.py --offsets):
 (a) STACK-ALIGNMENT RESIDUAL HISTOGRAM -- the E3 step-5 method (peak_offset copied verbatim from e3_bostock_timebase.py):
     for every build-half event, the lag of the maximum 3-C normalised cc between the family's SeisBench PGC stack
     (S-1 .. S+7 s) and the PGC 2012 data within +-40 s of the ZERO-offset prediction DB + ttP(Crosson) + (S-P)_stack - 1 s;
     the residual distribution is scanned with a 2-s window stepped 0.25 s and the peak accepted at Poisson z >= 4, mu >= 1.
     Only for the 126 families with a PGC stack.  Convention identical to E3/E4 (origin + ttP + S-P_stack = the stack's S).
 (b) FROM-SCRATCH ENVELOPE PEAK -- linear stack of the rms-normalised PGC 3-C data at DB + ttS(Crosson) over +-40 s, NO
     alignment; the 0.5-s-smoothed 3-C envelope peak within +-20 s of the nominal S; offset_b = t_peak - t_nominal; its
     strength = envelope peak / median envelope (dB), compared with the same statistic on a same-N random-time stack.
     All 130 families; the working offset for the families without a PGC stack.
Drift = offset_a(2012) - offset(2010) for the families with an accepted E3 2010 offset.
Outputs results/e5r_family_offsets_2012.csv (bkey, sb_fam_id, n_build, offset_a, accept_a, sigma_a, peak_over_bg_a, n_pairs_a,
        offset_b, env_db_b, env_db_null_b, accept_b, offset_2010, drift_a, offset_s [working], src, accept) ; figures/e5r_offsets.png
"""
from __future__ import annotations
import glob, os, sys
import numpy as np, pandas as pd, h5py
from scipy.signal import butter, sosfiltfilt, hilbert
ROOT = "/home/jovyan/tremorferometry/lfe_catalog_repair"; RES = f"{ROOT}/results"; CACHE = f"{ROOT}/cache20hz"
SD = f"{RES}/external/other/seisbench_lfestacks_cascadia"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d3_vel1d as v1d  # noqa
FS = 20.0; R = 111.195; MAXLAG = 40.0
PGC = (48.6498, -123.4521)


def bp20(x):
    sos = butter(4, [1.0 / (FS / 2), 8.0 / (FS / 2)], btype="band", output="sos"); return sosfiltfilt(sos, x, axis=-1)


def epi_km(lat1, lon1, lat2, lon2):
    cs = np.cos(np.radians(48.6)); return np.hypot((lon1 - lon2) * R * cs, (lat1 - lat2) * R)


def peak_offset(d, maxlag=MAXLAG, half=1.0, step=0.25, bg_excl=5.0, nsigma=4.0):
    """VERBATIM from scripts/e3_bostock_timebase.py (E3 step 5, bug-fixed version)."""
    out = dict(offset_s=np.nan, peak=0, bg_per_s=np.nan, mu=np.nan, peak_over_bg=np.nan, sigma=np.nan, n_pairs=int(len(d)), excess=np.nan, accept=False)
    if len(d) < 20:
        return out
    d = np.sort(d)
    cen = np.arange(-maxlag + half, maxlag - half + step / 2, step)
    lo = np.searchsorted(d, cen - half); hi = np.searchsorted(d, cen + half); n = hi - lo
    k = int(np.argmax(n))
    n_far = int(((d < cen[k] - bg_excl) | (d > cen[k] + bg_excl)).sum())
    span = 2 * maxlag - 2 * bg_excl
    bg_per_s = n_far / span if span > 0 else np.nan
    mu = bg_per_s * 2 * half
    z = (n[k] - mu) / np.sqrt(mu) if mu > 0 else np.nan
    out.update(offset_s=float(cen[k]), peak=int(n[k]), bg_per_s=float(bg_per_s), mu=float(mu), peak_over_bg=float(n[k] / mu) if mu > 0 else np.nan,
               sigma=float(z) if np.isfinite(z) else np.nan, excess=float(n[k] - mu), accept=bool(np.isfinite(z) and z >= nsigma and mu >= 1.0))
    return out


def load_days(netsta, year=2012):
    D = {}
    for f in sorted(glob.glob(f"{CACHE}/{netsta}/{year}/*.npz")):
        z = np.load(f); D[int(os.path.basename(f)[:3])] = (float(z["t0"]), z["x"])
    return D


def cut(D, starts, n, normalise=True):
    segs, keep = [], []
    for j, ts in enumerate(starts):
        v = D.get(pd.Timestamp(ts, unit="s").dayofyear)
        if v is None:
            continue
        t0, X = v; i0 = int(round((ts - t0) * FS))
        if i0 < 0 or i0 + n > X.shape[1]:
            continue
        seg = X[:, i0:i0 + n].astype(np.float64); r = np.sqrt((seg ** 2).mean())
        if not np.isfinite(r) or r <= 0:
            continue
        segs.append(seg / r if normalise else seg); keep.append(j)
    return (np.stack(segs), np.asarray(keep)) if segs else (None, np.zeros(0, int))


def batch_cc_lags(S, tpl):
    """S (N,3,nw) windows, tpl (3,nt) per-component L2-normalised -> per event the lag (samples) of the max 3-C mean normalised cc."""
    N, nc, nw = S.shape; nt = tpl.shape[1]; nl = nw - nt + 1
    NF = int(2 ** np.ceil(np.log2(nw + nt)))
    num = np.fft.irfft(np.fft.rfft(S, n=NF, axis=-1) * np.conj(np.fft.rfft(tpl, n=NF, axis=-1))[None], n=NF, axis=-1)[:, :, :nl]
    c = np.concatenate([np.zeros((N, nc, 1)), np.cumsum(S * S, axis=-1)], axis=-1)
    e = np.sqrt(np.maximum(c[:, :, nt:nt + nl] - c[:, :, :nl], 1e-12 * c[:, :, -1:]))
    cc = (num / e).mean(1)
    return cc.argmax(1), cc.max(1)


def main():
    D = load_days("CN.PGC")
    print(f"[data] PGC 2012: {len(D)} cached days", flush=True)
    gt = pd.read_parquet(f"{RES}/e5r_gt_bostock2012.parquet")
    F = pd.read_csv(f"{RES}/e5r_families_2012.csv")
    B = gt[(gt.half == "build") & (~gt.guard)].copy()
    B["t_db"] = B.time.astype("datetime64[ns]").astype("int64") / 1e9          # DB time (no offset)
    fo = pd.read_csv(f"{RES}/e3_bostock_family_offsets.csv"); acc = fo[fo.accept.astype(bool)]
    off2010 = dict(zip(acc.fam_id.astype(float), acc.offset_s.astype(float)))
    # SeisBench PGC stacks
    m = pd.read_csv(f"{SD}/metadata.csv"); m = m[(m.station_code == "PGC") & m.trace_S_arrival_sample.notna()]
    SB = {}
    with h5py.File(f"{SD}/waveforms.hdf5", "r") as f:
        for r in m.itertuples():
            b, idx = r.trace_name.split("$"); p = idx.split(",")
            w = bp20(np.asarray(f["data"][b][int(p[0]), :3, :int(p[2].lstrip(":"))], float))
            SB[float(r.source_family_id)] = (w, float(r.trace_S_arrival_sample))
    rng = np.random.default_rng(7)
    nw = int(round(2 * MAXLAG * FS)); nt = int(8 * FS); idx_nom = nw // 2
    tmin, tmax = B.t_db.min(), B.t_db.max()
    rows = []
    for fr in F.itertuples():
        g = B[B.bkey == fr.bkey]
        dx = epi_km(fr.lat, fr.lon, *PGC); ttP = float(v1d.tt(dx, fr.dep, "P")); ttS = float(v1d.tt(dx, fr.dep, "S"))
        row = dict(bkey=fr.bkey, sb_fam_id=fr.sb_fam_id, n_build=len(g), epi_km=dx, offset_2010=off2010.get(fr.sb_fam_id, np.nan))
        # (a) alignment-residual histogram against the SeisBench PGC stack
        ra = dict(offset_s=np.nan, accept=False, sigma=np.nan, peak_over_bg=np.nan, n_pairs=0)
        if np.isfinite(fr.sb_fam_id) and fr.sb_fam_id in SB and len(g) >= 20:
            w, s_samp = SB[fr.sb_fam_id]
            s0 = int(round(s_samp - 1.0 * FS)); tpl = w[:, s0:s0 + nt]; tpl = tpl / np.linalg.norm(tpl, axis=1, keepdims=True)
            sp = (s_samp - 800.0) / FS
            starts = g.t_db.to_numpy() + ttP + sp - 1.0 - MAXLAG
            S, keep = cut(D, starts, nw + nt - 1, normalise=False)
            if S is not None and len(S) >= 20:
                lag, ccmax = batch_cc_lags(S, tpl)
                d = lag / FS - MAXLAG
                ra = peak_offset(d)
                row["ccmax_median_a"] = float(np.median(ccmax))
        row.update(offset_a=ra["offset_s"], accept_a=ra["accept"], sigma_a=ra["sigma"], peak_over_bg_a=ra["peak_over_bg"], n_pairs_a=ra["n_pairs"])
        # (b) from-scratch envelope peak of the un-aligned stack at DB + ttS
        S, keep = cut(D, g.t_db.to_numpy() + ttS - MAXLAG, nw)
        ob, edb, edbn = np.nan, np.nan, np.nan
        if S is not None and len(S) >= 20:
            def envpk(st):
                env = np.sqrt((np.abs(hilbert(st, axis=-1)) ** 2).sum(0)); k = int(0.5 * FS)
                env = np.convolve(env, np.ones(k) / k, mode="same"); wq = int(20 * FS)
                a, b = idx_nom - wq, idx_nom + wq; j = a + int(np.argmax(env[a:b]))
                return j, float(20 * np.log10(env[j] / np.median(env)))
            j, edb = envpk(S.mean(0)); ob = (j - idx_nom) / FS
            Sn, _ = cut(D, rng.uniform(tmin, tmax, len(S)) + ttS - MAXLAG, nw)
            if Sn is not None and len(Sn) >= 10:
                edbn = envpk(Sn.mean(0))[1]
        row.update(offset_b=ob, env_db_b=edb, env_db_null_b=edbn, accept_b=bool(np.isfinite(edb) and np.isfinite(edbn) and edb >= edbn + 3.0))
        rows.append(row)
        print(f"  bkey {fr.bkey:3d} sb {fr.sb_fam_id if np.isfinite(fr.sb_fam_id) else -1:5.0f} n {len(g):3d}  a: {ra['offset_s']!s:>6} z {ra['sigma']!s:>6}  "
              f"b: {ob:+6.2f} ({edb:5.1f} vs null {edbn:5.1f} dB)  2010: {row['offset_2010']}", flush=True)
    T = pd.DataFrame(rows)
    T["drift_a"] = T.offset_a - T.offset_2010
    T["a_minus_b"] = T.offset_a - T.offset_b
    # working offset: (a) if accepted; else (b) if accepted; else 2010; else median of the accepted (a)
    med_a = float(T.loc[T.accept_a, "offset_a"].median())
    # RULE (set before the templates were built, run 2): (a) depends on the SeisBench stack, which is known to be bad for some
    # families (E4/E5: family 15 = bkey 59, S pick 0.9-2.2 s off, cc < 0.3) -- run 1 gave it offset_a 11.0 (z 7.6) against
    # offset_b 5.5 and the 2010 value 5.25.  So: (a) only when (b) is not accepted or agrees within 2 s; otherwise (b).
    use_a = T.accept_a & (~T.accept_b | (T.a_minus_b.abs() <= 2.0))
    use_b = ~use_a & T.accept_b
    src = np.where(use_a, "2012_a", np.where(use_b, "2012_b_envelope", np.where(T.offset_2010.notna(), "2010", "median_2012_a")))
    T["offset_s"] = np.select([use_a, use_b, T.offset_2010.notna()], [T.offset_a, T.offset_b, T.offset_2010], med_a)
    T["src"] = src; T["accept"] = True
    T.to_csv(f"{RES}/e5r_family_offsets_2012.csv", index=False)
    ok = T[T.accept_a]
    both = T[T.accept_a & T.offset_2010.notna()]
    print(f"\n[a] accepted {len(ok)} of {int(T.sb_fam_id.notna().sum())} families with a PGC stack; offset median {ok.offset_a.median():+.2f} s, "
          f"IQR {ok.offset_a.quantile(.25):+.2f}..{ok.offset_a.quantile(.75):+.2f}, range {ok.offset_a.min():+.2f}..{ok.offset_a.max():+.2f}; "
          f"sigma median {ok.sigma_a.median():.1f}, peak/bg median {ok.peak_over_bg_a.median():.1f}")
    print(f"[drift] 2012(a) - 2010 over {len(both)} families: median {both.drift_a.median():+.2f} s, MAD {1.4826*np.median(np.abs(both.drift_a-both.drift_a.median())):.2f}, "
          f"|drift| <= 0.5 s: {int((both.drift_a.abs() <= 0.5).sum())}, <= 1 s: {int((both.drift_a.abs() <= 1.0).sum())}, > 2 s: {int((both.drift_a.abs() > 2.0).sum())}")
    ab = T[T.accept_a & T.accept_b]
    print(f"[a vs b] over {len(ab)} families: a - b median {ab.a_minus_b.median():+.2f} s, IQR {ab.a_minus_b.quantile(.25):+.2f}..{ab.a_minus_b.quantile(.75):+.2f}; "
          f"(b) accepted {int(T.accept_b.sum())} of {len(T)}")
    print(f"[working] sources: {T.src.value_counts().to_dict()}")
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].scatter(both.offset_2010, both.offset_a, c="#0072B2", s=18); lim = [min(both.offset_2010.min(), both.offset_a.min()) - 1, max(both.offset_2010.max(), both.offset_a.max()) + 1]
    ax[0].plot(lim, lim, "k--", lw=0.8); ax[0].set_xlabel("2010 offset at PGC (s, E3 step 5)"); ax[0].set_ylabel("2012 build-half offset (s, method a)"); ax[0].set_title(f"per-family time base, {len(both)} families")
    ax[1].hist(both.drift_a.dropna(), bins=np.arange(-3, 3.01, 0.25), color="#E69F00"); ax[1].set_xlabel("drift 2012 - 2010 (s)"); ax[1].set_ylabel("families"); ax[1].set_title(f"median {both.drift_a.median():+.2f} s")
    ax[2].scatter(ab.offset_b, ab.offset_a, c="#009E73", s=18); ax[2].plot(lim, lim, "k--", lw=0.8); ax[2].set_xlabel("offset (b): envelope peak at DB + ttS (s)"); ax[2].set_ylabel("offset (a) (s)"); ax[2].set_title(f"a - b median {ab.a_minus_b.median():+.2f} s")
    for x in ax: x.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{ROOT}/figures/e5r_offsets.png", dpi=130)
    print(f"[write] results/e5r_family_offsets_2012.csv, figures/e5r_offsets.png")


if __name__ == "__main__":
    main()
