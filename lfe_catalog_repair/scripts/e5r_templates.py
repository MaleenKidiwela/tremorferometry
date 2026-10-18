#!/usr/bin/env python
"""E5 step 5 STEP C -- TEMPLATES FROM SCRATCH for all 15 stations of the 2012 replication set, from the 2012 BUILD half only.

For every Bostock 2012 family with >= 30 build-half events (alternate 2-h blocks, 60-s guard; e5r_gt_bostock2012.parquet) and
every station of results/e5r_stations.csv:
  * stack the station's 3-C, 1-8 Hz, 20-Hz data (cache20hz) at the PREDICTED S arrival = t_origin + ttS(Crosson 1-D, family
    hypocentre -> station), window +-20 s, every event rms-normalised, linear stack;
  * refine the per-(family, station) alignment by cc of every individual event against the running stack, +-1.5 s, 2 iterations;
  * pick S on the 0.5-s-smoothed 3-C envelope of the stack within +-4 s of the nominal S;
  * stack S-SNR (dB) = energy S-0.5..S+3.5 over the same window 15 s earlier, vs a SAME-N RANDOM-TIME null stack put through the
    identical alignment + pick (stricter than E4's un-aligned null); S-SNR of nested random sub-stacks at N = 30/60/100/200;
  * q-control at PGC, VGZ, SNB: cc of our 2012 re-stack against the SeisBench (Bostock) stack of the same family over the 8 s
    from S-1, +-1 s relative shift (expect >= 0.6 at N >= 200, as E4/E5);
  * templates 8 s and 12 s from S-1 s (S-1..S+7 / S-1..S+11), each component L2-normalised; a template is USABLE when
    S-SNR >= null S-SNR + 3 dB; only usable templates enter the network.  NO SeisBench stack is used as a template.
Also the per-station AMPLITUDE REFERENCE for the catalog's relative-amplitude proxy: the median sliding 3-C window energy
(8 s and 12 s, 1-s stride) over the pre-SSE quiet days 2012-08-15..08-29 (doy 228-242).
Outputs results/e5r_templates_2012.npz (tpl {(bkey, sta): {8.0, 12.0}}, sp -- placement S-P for e4_template_test.py, stacks
        "<sta>|<bkey>"), results/e5r_template_qc_2012.csv, results/e5r_weights_snr_2012.csv (usable pairs: w = 10^(S-SNR/20),
        capped at 3x the median, median-normalised -- BUILD half only), results/e5r_amp_ref_2012.csv (sta, L, e3_ref)
PRE bookkeeping: e4_template_test.py places a template at origin + ttP + sp - 0.5 s; our templates start at S - 1.0 s, so
        sp = (ttS - ttP) + (S pick - nominal) - 0.5  ->  start = origin + ttS + pick refinement - 1.0.  (Commented at the site.)
"""
from __future__ import annotations
import argparse, glob, os, sys, time
from multiprocessing import Pool
import numpy as np, pandas as pd, h5py
from scipy.signal import butter, sosfiltfilt, hilbert
ROOT = "/home/jovyan/tremorferometry/lfe_catalog_repair"; RES = f"{ROOT}/results"; CACHE = f"{ROOT}/cache20hz"
SD = f"{RES}/external/other/seisbench_lfestacks_cascadia"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d3_vel1d as v1d  # noqa
FS = 20.0; HALF_W = 20.0; PRE_TPL = 1.0; PRE_TT = 0.5; R_KM_DEG = 111.195
QUIET_DOYS = range(228, 243)
LENS = (8.0, 12.0)
QC_STAS = ("PGC", "VGZ", "SNB")


def bp20(x):
    sos = butter(4, [1.0 / (FS / 2), 8.0 / (FS / 2)], btype="band", output="sos"); return sosfiltfilt(sos, x, axis=-1)


def epi_km(lat1, lon1, lat2, lon2):
    cs = np.cos(np.radians(48.6)); return float(np.hypot((lon1 - lon2) * R_KM_DEG * cs, (lat1 - lat2) * R_KM_DEG))


def load_days(netsta, year=2012):
    D = {}
    for f in sorted(glob.glob(f"{CACHE}/{netsta}/{year}/*.npz")):
        z = np.load(f); D[int(os.path.basename(f)[:3])] = (float(z["t0"]), z["x"])
    return D


def cut_windows(D, starts, n):
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
        segs.append(seg / r); keep.append(j)
    return (np.stack(segs), np.asarray(keep)) if segs else (None, np.zeros(0, int))


def align_stack(S, iters=2, lag=30):
    """rms-normalised windows (N,3,n) -> cc-aligned linear stack (3,n), the aligned windows, per-event shifts (samples).
    roll(sg, k)[t] = sg[t-k]; cc(k) = sum(roll(sg,k) * ref) = sum(sg * roll(ref,-k))  (same convention as e4_transfer.stack_family)."""
    N, nc, n = S.shape; tot = np.zeros(N, int); st = S.mean(0); ar = np.arange(n)
    for _ in range(iters):
        ref = st - st.mean(1, keepdims=True)
        ccs = np.empty((2 * lag + 1, N))
        for kk, k in enumerate(range(-lag, lag + 1)):
            ccs[kk] = np.einsum("ijk,jk->i", S, np.roll(ref, -k, axis=1))
        best = ccs.argmax(0) - lag
        idx = (ar[None, :] - best[:, None]) % n
        S = S[np.arange(N)[:, None, None], np.arange(nc)[None, :, None], idx[:, None, :]]
        tot += best; st = S.mean(0)
    return st, S, tot


def envelope(st, smooth_s=0.5):
    env = np.sqrt((np.abs(hilbert(st, axis=-1)) ** 2).sum(0)); k = int(round(smooth_s * FS))
    return np.convolve(env, np.ones(k) / k, mode="same")


def pick_s(st, nominal_idx, search_s=4.0):
    env = envelope(st); w = int(round(search_s * FS)); a, b = max(nominal_idx - w, 0), min(nominal_idx + w, len(env))
    return int(a + np.argmax(env[a:b]))


def s_snr_db(st, s_idx, win_s=4.0, noise_off_s=15.0):
    n = int(round(win_s * FS)); a = max(s_idx - int(0.5 * FS), 0); sig = st[:, a:a + n]
    b = max(a - int(noise_off_s * FS), 0); noi = st[:, b:b + n]
    if sig.size == 0 or noi.size == 0 or (noi ** 2).mean() <= 0:
        return np.nan
    return float(10 * np.log10((sig ** 2).mean() / (noi ** 2).mean()))


def do_station(args):
    sta, netsta, lat_s, lon_s, klass, fam_tbl, ev, SB, seed = args
    t_st = time.time()
    D = load_days(netsta)
    if not D:
        return sta, [], {}, {}, {}, {}
    # amplitude reference: median sliding 3-C window energy over the quiet days
    ref = {}
    for L in LENS:
        n = int(L * FS); meds = []
        for doy in QUIET_DOYS:
            if doy not in D:
                continue
            X = D[doy][1].astype(np.float64)
            c = np.concatenate([np.zeros((3, 1)), np.cumsum(X * X, axis=1)], axis=1)
            ew = (c[:, n:] - c[:, :-n]).sum(0)[::int(FS)]
            ew = ew[ew > 0]
            if len(ew):
                meds.append(float(np.median(ew)))
        ref[L] = float(np.median(meds)) if meds else np.nan
    rng = np.random.default_rng(seed)
    n = int(round(2 * HALF_W * FS)); idx_nom = int(round(HALF_W * FS))
    tmin, tmax = float(ev.t_origin.min()), float(ev.t_origin.max())
    rows, TPL, SP, STK = [], {}, {}, {}
    for fr in fam_tbl:
        g = ev[ev.bkey == fr["bkey"]]
        if len(g) < 30:
            continue
        dx = epi_km(fr["lat"], fr["lon"], lat_s, lon_s)
        ttP = float(v1d.tt(dx, fr["dep"], "P")); ttS = float(v1d.tt(dx, fr["dep"], "S"))
        S, keep = cut_windows(D, g.t_origin.to_numpy() + ttS - HALF_W, n)
        if S is None or len(S) < 20:
            continue
        st, Sal, sh = align_stack(S, iters=2)
        sp_idx = pick_s(st, idx_nom); snr = s_snr_db(st, sp_idx)
        Sn, _ = cut_windows(D, rng.uniform(tmin, tmax, len(S)) + ttS - HALF_W, n)
        snr_n = np.nan
        if Sn is not None and len(Sn) >= 10:
            stn, _, _ = align_stack(Sn, iters=2); snr_n = s_snr_db(stn, pick_s(stn, idx_nom))
        sub = {}
        for Nk in (30, 60, 100, 200):
            sub[f"snr_N{Nk}"] = s_snr_db(Sal[rng.choice(len(Sal), Nk, replace=False)].mean(0), sp_idx) if Nk < len(Sal) else np.nan
        ccsb = np.nan
        if SB and np.isfinite(fr["sb_fam_id"]) and fr["sb_fam_id"] in SB:
            sb, s_samp = SB[fr["sb_fam_id"]]
            i_sb = int(round(s_samp - 1.0 * FS)); n8 = int(8 * FS); i_re = sp_idx - int(1.0 * FS)
            if i_sb >= 0 and i_sb + n8 <= sb.shape[1] and i_re >= 0 and i_re + n8 <= st.shape[1]:
                A = sb[:, i_sb:i_sb + n8]; A = A - A.mean(1, keepdims=True)
                Bm = st[:, i_re:i_re + n8]; Bm = Bm - Bm.mean(1, keepdims=True)
                best = -9.0
                for kk in range(-int(FS), int(FS) + 1):
                    Bs = np.roll(Bm, kk, axis=1)
                    best = max(best, float((A * Bs).sum() / max(np.sqrt((A ** 2).sum() * (Bs ** 2).sum()), 1e-20)))
                ccsb = best
        usable = bool(np.isfinite(snr) and np.isfinite(snr_n) and snr >= snr_n + 3.0)
        d, ok = {}, True
        for L in LENS:
            nL = int(round(L * FS)); i0 = sp_idx - int(round(PRE_TPL * FS))
            if i0 < 0 or i0 + nL > st.shape[1]:
                ok = False; break
            seg = st[:, i0:i0 + nL].astype(np.float32); nn = np.linalg.norm(seg, axis=1)
            if np.any(nn <= 0):
                ok = False; break
            d[L] = (seg / nn[:, None]).astype(np.float32)
        if ok:
            key = (float(fr["bkey"]), sta)
            TPL[key] = d
            # e4_template_test.py: start = origin + ttP + sp - 0.5; we want start = origin + ttS + (pick - nominal) - 1.0
            SP[key] = (ttS - ttP) + (sp_idx - idx_nom) / FS - (PRE_TPL - PRE_TT)
            STK[f"{sta}|{int(fr['bkey'])}"] = st.astype(np.float32)
        rows.append(dict(sta=sta, klass=klass, bkey=int(fr["bkey"]), sb_fam_id=fr["sb_fam_id"], n_build=len(g), n_events=len(S), epi_km=dx,
                         ttP=ttP, ttS=ttS, align_shift_s=float(np.median(sh)) / FS, s_pick_s=(sp_idx - idx_nom) / FS, snr_db=snr,
                         snr_null_db=snr_n, usable=usable, template=ok, cc_vs_seisbench=ccsb, **sub))
    print(f"  [{sta}] {len(rows)} stacks, {sum(r['usable'] and r['template'] for r in rows)} usable, {time.time()-t_st:.0f} s", flush=True)
    return sta, rows, TPL, SP, STK, ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-events", type=int, default=30)
    ap.add_argument("--tag", default="2012")
    ap.add_argument("--workers", type=int, default=15)
    a = ap.parse_args()
    S = pd.read_csv(f"{RES}/e5r_stations.csv")
    F = pd.read_csv(f"{RES}/e5r_families_2012.csv")
    gt = pd.read_parquet(f"{RES}/e5r_gt_bostock2012.parquet")
    ev = gt[(gt.half == "build") & (~gt.guard)][["bkey", "t_origin"]].reset_index(drop=True)
    fam_tbl = F[F.n_build >= a.min_events][["bkey", "sb_fam_id", "lat", "lon", "dep", "n_build"]].to_dict("records")
    print(f"[fam] {len(fam_tbl)} families with >= {a.min_events} build-half events ({len(ev)} build events); {len(S)} stations", flush=True)
    m = pd.read_csv(f"{SD}/metadata.csv"); m = m[m.station_code.isin(QC_STAS) & m.trace_S_arrival_sample.notna()]
    SBall = {s: {} for s in QC_STAS}
    with h5py.File(f"{SD}/waveforms.hdf5", "r") as f:
        for r in m.itertuples():
            b, idx = r.trace_name.split("$"); p = idx.split(",")
            w = bp20(np.asarray(f["data"][b][int(p[0]), :3, :int(p[2].lstrip(":"))], float))
            SBall[r.station_code][float(r.source_family_id)] = (w, float(r.trace_S_arrival_sample))
    jobs = [(r.sta, r.netsta, float(r.lat), float(r.lon), r.klass, fam_tbl, ev, SBall.get(r.sta, {}), 100 + i) for i, r in enumerate(S.itertuples())]
    with Pool(min(a.workers, len(jobs))) as pool:
        out = pool.map(do_station, jobs)
    rows, TPL, SP, STK, REF = [], {}, {}, {}, []
    for sta, rw, T, P, K, ref in out:
        rows += rw; TPL.update(T); SP.update(P); STK.update(K)
        for L, v in ref.items():
            REF.append(dict(sta=sta, L=L, e3_ref=v))
    Q = pd.DataFrame(rows); Q.to_csv(f"{RES}/e5r_template_qc_{a.tag}.csv", index=False)
    pd.DataFrame(REF).to_csv(f"{RES}/e5r_amp_ref_{a.tag}.csv", index=False)
    use = Q[Q.usable & Q.template]
    TPLu = {k: v for k, v in TPL.items() if ((use.bkey == int(k[0])) & (use.sta == k[1])).any()}
    SPu = {k: SP[k] for k in TPLu}
    np.savez_compressed(f"{RES}/e5r_templates_{a.tag}.npz", tpl=np.array(TPLu, dtype=object), sp=np.array(SPu, dtype=object), **STK)
    w = 10 ** (use.snr_db / 20.0); w = np.minimum(w, 3 * w.median()); w = w / w.median()
    W = pd.DataFrame(dict(fam=use.bkey.astype(float), sta=use.sta, w=w.to_numpy(), snr_db=use.snr_db, n_events=use.n_events))
    W.to_csv(f"{RES}/e5r_weights_snr_{a.tag}.csv", index=False)
    pd.set_option("display.width", 250)
    print("\n[coverage] per station: stacks / usable templates / median S-SNR (dB) / median null / median cc vs SeisBench")
    print(Q.groupby("sta").agg(klass=("klass", "first"), n=("bkey", "size"), usable=("usable", "sum"), nev=("n_events", "median"), snr=("snr_db", "median"),
                               snr_null=("snr_null_db", "median"), ccsb=("cc_vs_seisbench", "median"), spick=("s_pick_s", "median"), epi=("epi_km", "median")).round(2).to_string())
    per_fam = use.groupby("bkey").size()
    print(f"\n[coverage] {len(TPLu)} usable (family, station) templates; families with >= 1 usable template {len(per_fam)}; stations per family "
          f"p10/50/90 = {np.percentile(per_fam, [10, 50, 90])}; families with >= 5 stations {int((per_fam >= 5).sum())}, >= 8: {int((per_fam >= 8).sum())}")
    qc = Q[Q.sta.isin(QC_STAS) & Q.cc_vs_seisbench.notna()]
    for s in QC_STAS:
        h = qc[qc.sta == s]; big = h[h.n_events >= 200]
        print(f"[q-control {s}] {len(h)} pairs: cc median {h.cc_vs_seisbench.median():.2f}, >= 0.6: {int((h.cc_vs_seisbench >= 0.6).sum())}; "
              f"N >= 200 ({len(big)}): median {big.cc_vs_seisbench.median() if len(big) else float('nan'):.2f}, >= 0.6: {int((big.cc_vs_seisbench >= 0.6).sum())}; "
              f"N >= 100: median {h[h.n_events >= 100].cc_vs_seisbench.median():.2f}")
    print(f"[weights] {len(W)} pairs; per-station median w: {W.groupby('sta').w.median().round(2).to_dict()}")
    print(f"[write] results/e5r_templates_{a.tag}.npz, e5r_template_qc_{a.tag}.csv, e5r_weights_snr_{a.tag}.csv, e5r_amp_ref_{a.tag}.csv")


if __name__ == "__main__":
    main()
