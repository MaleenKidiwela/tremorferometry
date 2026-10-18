#!/usr/bin/env python
"""E4 step 2 — (a) the q-CONTROL and (b) TEMPLATE TRANSFER to stations Bostock never had.

Both legs stack the 2010 data at Bostock event times, on a BUILD half, and are scored on a held-out
TEST half.  The hold-out is ALTERNATE 2-HOUR BLOCKS with a 60-s guard at every block edge (Merlin
vet F2: the Aug-build / Sep-test split of the original brief is impossible -- the 2010 ETS is over
by 1 September, the top families have 235-327 August events and 1-35 September events, and splitting
by event parity is invalid because 35 % of events have a same-family neighbour within 30 s).

(a) q-CONTROL, where truth exists.  Re-stack each RING station's own 2010 data at the build-half
    Bostock times of its own families and compare with the published SeisBench stack:
      (i)  cc(2010 re-stack, SeisBench stack) per (family, station) vs the number of events stacked.
           Merlin's bar: >= 0.6 for N >= 200.  cc < 0.3 at N >= 200 would mean the time-base offset
           sign or the S-P application is wrong -- STOP and report.
      (ii) the held-out registration ratio (2010 template / SeisBench template on the SAME test-half
           events) IS the measured template quality q, to be compared with the prediction
           q = 1/sqrt(1 + 1/(r1^2 N)).
(b) TRANSFER to B011, B926, LRIV (3-C) and STW (Z only).  NLLB (6 families within 40 km), OCP (0)
    and GOBB (EHZ only, and no families) are dropped on Merlin's F8 geometry count.  Stack the
    build-half events at the Crosson-predicted S arrival, +-20 s, iterate the alignment by cc
    against the running stack, pick S on the stack envelope, and cut 6- and 8-s templates.
    Reported: stack S-SNR vs the number of events stacked, and vs a same-N random-time null stack.

Outputs
  results/e4_transfer_stacks_<tag>.npz    the stacks and the templates (fed to e4_template_test.py
                                          through --extra-tpl)
  results/e4_transfer_stackqc_<tag>.csv   per (family, station): N, S-SNR, null S-SNR, cc vs the
                                          SeisBench stack, the S pick, the alignment shift
"""
from __future__ import annotations
import argparse, os, sys, time
import numpy as np, pandas as pd, h5py
from scipy.signal import butter, sosfiltfilt, hilbert

ROOT = "/home/jovyan/tremorferometry/lfe_catalog_repair"
RES = f"{ROOT}/results"
CACHE = f"{ROOT}/cache20hz"
SD = f"{RES}/external/other/seisbench_lfestacks_cascadia"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d3_vel1d as v1d  # noqa: E402

FS = 20.0
HALF_W = 20.0                 # +- 20 s window around the predicted S
R_KM_DEG = 111.195
NETSTA = {"PGC": "CN.PGC", "LZB": "CN.LZB", "VGZ": "CN.VGZ", "SNB": "CN.SNB", "YOUB": "CN.YOUB",
          "PFB": "CN.PFB", "B011": "PB.B011", "B926": "PB.B926", "LRIV": "UW.LRIV",
          "STW": "UW.STW", "NLLB": "CN.NLLB", "OCP": "UW.OCP"}
RING = ["PGC", "LZB", "VGZ", "SNB", "YOUB", "PFB"]


def bp20(x):
    sos = butter(4, [1.0 / (FS / 2), 8.0 / (FS / 2)], btype="band", output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def epi_km(lat1, lon1, lat2, lon2):
    cs = np.cos(np.radians(48.6))
    return np.hypot((np.asarray(lon1) - lon2) * R_KM_DEG * cs,
                    (np.asarray(lat1) - lat2) * R_KM_DEG)


def seisbench_stacks():
    m = pd.read_csv(f"{SD}/metadata.csv")
    m = m[m.trace_S_arrival_sample.notna()]
    W, SPs, SN = {}, {}, {}
    with h5py.File(f"{SD}/waveforms.hdf5", "r") as f:
        for r in m.itertuples():
            b, idx = r.trace_name.split("$")
            p = idx.split(",")
            w = np.asarray(f["data"][b][int(p[0]), :3, :int(p[2].lstrip(":"))], float)
            W[(float(r.source_family_id), r.station_code)] = bp20(w)
            SPs[(float(r.source_family_id), r.station_code)] = \
                float(r.trace_S_arrival_sample - 800.0) / FS
            SN[(float(r.source_family_id), r.station_code)] = \
                (float(r.trace_S_snr_db), float(r.trace_stack_count))
    fam = m.groupby("source_family_id")[["source_latitude_deg", "source_longitude_deg",
                                         "source_depth_km"]].first()
    return W, SPs, SN, fam


def load_day(netsta, year, doy, cache):
    k = (netsta, doy)
    if k in cache:
        return cache[k]
    p = f"{CACHE}/{netsta}/{year}/{doy:03d}.npz"
    v = None
    if os.path.exists(p):
        z = np.load(p)
        v = (float(z["t0"]), z["x"])
    cache[k] = v
    return v


def stack_family(times, netsta, year, t_pred_fn, cache, iters=3, align_s=1.5):
    """Linear stack of rms-normalised +-HALF_W windows, iteratively cc-aligned. -> (nc, N), n_used."""
    n = int(round(2 * HALF_W * FS))
    segs = []
    for t in times:
        ts = t + t_pred_fn(t) - HALF_W
        d = pd.Timestamp(ts, unit="s")
        v = load_day(netsta, year, d.dayofyear, cache)
        if v is None:
            continue
        t0, X = v
        i0 = int(round((ts - t0) * FS))
        if i0 < 0 or i0 + n > X.shape[1]:
            continue
        seg = X[:, i0:i0 + n].astype(np.float64)
        r = np.sqrt((seg ** 2).mean())
        if not np.isfinite(r) or r <= 0:
            continue
        segs.append(seg / r)
    if len(segs) < 5:
        return None, 0, 0.0
    S = np.stack(segs)                          # (n_ev, nc, n)
    st = S.mean(0)
    sh_tot = 0.0
    lag = int(round(align_s * FS))
    for _ in range(iters):
        ref = st - st.mean(1, keepdims=True)
        rn = np.sqrt((ref ** 2).sum())
        best = np.zeros(len(S), int)
        for i, sg in enumerate(S):
            cc = np.zeros(2 * lag + 1)
            for k in range(-lag, lag + 1):
                a = np.roll(sg, k, axis=1)
                cc[k + lag] = float((a * ref).sum())
            best[i] = int(np.argmax(cc)) - lag
        S = np.stack([np.roll(sg, k, axis=1) for sg, k in zip(S, best)])
        st = S.mean(0)
        sh_tot += float(np.median(best)) / FS
    return st, len(S), sh_tot


def s_snr_db(st, s_idx, win_s=4.0, noise_off_s=15.0):
    n = int(round(win_s * FS))
    a = max(s_idx - int(0.5 * FS), 0)
    sig = st[:, a:a + n]
    b = max(a - int(noise_off_s * FS), 0)
    noi = st[:, b:b + n]
    if sig.size == 0 or noi.size == 0 or (noi ** 2).mean() <= 0:
        return np.nan
    return float(10 * np.log10((sig ** 2).mean() / (noi ** 2).mean()))


def pick_s(st, nominal_idx, search_s=3.0):
    """S pick = peak of the smoothed 3-C envelope within +-search_s of the nominal position."""
    env = np.sqrt((np.abs(hilbert(st, axis=-1)) ** 2).sum(0))
    k = int(round(0.5 * FS))
    env = np.convolve(env, np.ones(k) / k, mode="same")
    w = int(round(search_s * FS))
    a, b = max(nominal_idx - w, 0), min(nominal_idx + w, len(env))
    return int(a + np.argmax(env[a:b]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", default="PGC,LZB,VGZ,SNB,YOUB,PFB,B011,B926,LRIV,STW")
    ap.add_argument("--year", type=int, default=2010)
    ap.add_argument("--lengths", default="4,6,8")
    ap.add_argument("--min-events", type=int, default=100)
    ap.add_argument("--block-s", type=float, default=7200.0)
    ap.add_argument("--guard-s", type=float, default=60.0)
    ap.add_argument("--half", default="build")
    ap.add_argument("--tag", default="t2")
    a = ap.parse_args()
    STAS = a.stations.split(",")
    LENS = [float(x) for x in a.lengths.split(",")]
    PRE = 0.5

    W, SPs, SN, famxy = seisbench_stacks()
    stas = pd.read_csv(f"{RES}/e2_stations.csv").set_index("sta")
    bo = pd.read_parquet(f"{RES}/e2_gt_bostock.parquet")
    fo = pd.read_csv(f"{RES}/e3_bostock_family_offsets.csv")
    acc = fo[fo.accept.astype(bool)]
    omap = dict(zip(acc.bkey.astype(int), acc.offset_s.astype(float)))
    bo = bo[bo.bkey.astype(int).isin(omap)].copy()
    bo["time"] = bo.time.astype("datetime64[ns]")
    bo["t_origin"] = bo.time.astype("int64") / 1e9 + bo.bkey.astype(int).map(omap).to_numpy()
    blk = np.floor(bo.t_origin / a.block_s).astype("int64")
    pos = bo.t_origin - blk * a.block_s
    bo["half"] = np.where(blk % 2 == 0, "build", "test")
    bo["guard"] = (pos < a.guard_s) | (pos > a.block_s - a.guard_s)
    B = bo[(bo.half == a.half) & (~bo.guard) & bo.fam_id.isin(famxy.index)].copy()
    print(f"[evt] {len(B)} {a.half}-half Bostock events over {B.fam_id.nunique()} families "
          f"(2-h blocks, {a.guard_s:.0f}-s guard)", flush=True)

    cnt = B.groupby("fam_id").size()
    fams = sorted(cnt[cnt >= a.min_events].index)
    print(f"[fam] {len(fams)} families with >= {a.min_events} build-half events", flush=True)

    rng = np.random.default_rng(12)
    cache = {}
    rows, TPL, SPOUT, STACKS = [], {}, {}, {}
    t_start = time.time()
    for s in STAS:
        ns = NETSTA[s]
        if not os.path.isdir(f"{CACHE}/{ns}/{a.year}"):
            print(f"[warn] no cached {a.year} data for {ns}; skipped", flush=True)
            continue
        lat_s, lon_s = float(stas.loc[s, "lat"]), float(stas.loc[s, "lon"])
        for fid in fams:
            g = B[B.fam_id == fid]
            flat, flon, fdep = (float(famxy.loc[fid, "source_latitude_deg"]),
                                float(famxy.loc[fid, "source_longitude_deg"]),
                                float(famxy.loc[fid, "source_depth_km"]))
            dx = epi_km(flat, flon, lat_s, lon_s)
            ttP = float(v1d.tt(dx, fdep, "P"))
            ttS = float(v1d.tt(dx, fdep, "S"))
            has_sb = (fid, s) in W
            # RING (q-control): place the window on Bostock's OWN S-P from the stack pick
            # TRANSFER: place it on the Crosson-predicted S
            s_nom = (ttP + SPs[(fid, s)]) if has_sb else ttS
            st, nused, sh = stack_family(g.t_origin.to_numpy(), ns, a.year,
                                         lambda t, v=s_nom: v, cache)
            if st is None:
                continue
            idx_nom = int(round(HALF_W * FS))
            sp_idx = pick_s(st, idx_nom)
            snr = s_snr_db(st, sp_idx)
            # same-N random-time null stack
            tr = rng.uniform(B.t_origin.min(), B.t_origin.max(), nused)
            stn, nn, _ = stack_family(tr, ns, a.year, lambda t, v=s_nom: v, cache, iters=1)
            snr_n = s_snr_db(stn, idx_nom) if stn is not None else np.nan
            ccsb = np.nan
            if has_sb:
                # cc of the 2010 re-stack against the SeisBench stack over the 8 s from S-1 s
                sb = W[(fid, s)]
                i_sb = int(round(800 + SPs[(fid, s)] * FS - 1.0 * FS))
                n8 = int(8 * FS)
                i_re = sp_idx - int(1.0 * FS)
                if i_sb >= 0 and i_sb + n8 <= sb.shape[1] and i_re >= 0 and i_re + n8 <= st.shape[1]:
                    A = sb[:, i_sb:i_sb + n8]; Bm = st[:st.shape[0], i_re:i_re + n8]
                    nc = min(A.shape[0], Bm.shape[0])
                    A = A[:nc] - A[:nc].mean(1, keepdims=True)
                    Bm = Bm[:nc] - Bm[:nc].mean(1, keepdims=True)
                    # allow +-1 s of relative shift (our S pick vs Bostock's)
                    best = -9
                    for k in range(-int(FS), int(FS) + 1):
                        Bs = np.roll(Bm, k, axis=1)
                        d = float((A * Bs).sum() /
                                  max(np.sqrt((A ** 2).sum() * (Bs ** 2).sum()), 1e-20))
                        best = max(best, d)
                    ccsb = best
            r1 = np.nan
            if (fid, s) in SN:
                db, nst = SN[(fid, s)]
                r1 = 10 ** (db / 20.0) / np.sqrt(max(nst, 1))
            rows.append(dict(sta=s, fam=fid, n_events=nused, align_shift_s=sh,
                             s_pick_s=(sp_idx - idx_nom) / FS, snr_db=snr, snr_null_db=snr_n,
                             snr_ratio=10 ** ((snr - snr_n) / 20.0) if np.isfinite(snr_n) else np.nan,
                             cc_vs_seisbench=ccsb, epi_km=dx, has_seisbench=has_sb, r1=r1))
            # templates
            d = {}
            ok = True
            for L in LENS:
                n = int(round(L * FS))
                i0 = sp_idx - int(round(PRE * FS))
                if i0 < 0 or i0 + n > st.shape[1]:
                    ok = False
                    break
                seg = st[:, i0:i0 + n].astype(np.float32)
                nn_ = np.linalg.norm(seg, axis=1)
                if np.any(nn_ <= 0):
                    ok = False
                    break
                d[L] = (seg / nn_[:, None]).astype(np.float32)
            if ok:
                TPL[(fid, s)] = d
                # the offset formula in e4_template_test.py is ttP + sp - PRE, so store
                # sp = (ttS_pred - ttP) + the S-pick refinement
                SPOUT[(fid, s)] = (s_nom - ttP) + (sp_idx - idx_nom) / FS
                STACKS[f"{s}|{fid}"] = st.astype(np.float32)
        print(f"  [{s}] {sum(1 for k in TPL if k[1]==s)} templates built, "
              f"{time.time()-t_start:.0f} s", flush=True)

    Q = pd.DataFrame(rows)
    Q.to_csv(f"{RES}/e4_transfer_stackqc_{a.tag}.csv", index=False)
    np.savez_compressed(f"{RES}/e4_transfer_stacks_{a.tag}.npz",
                        tpl=np.array(TPL, dtype=object), sp=np.array(SPOUT, dtype=object),
                        **STACKS)
    print(f"[write] results/e4_transfer_stackqc_{a.tag}.csv, e4_transfer_stacks_{a.tag}.npz "
          f"({len(TPL)} templates)")
    pd.set_option("display.width", 200)
    print(Q.groupby("sta").agg(n=("fam", "size"), nev=("n_events", "median"),
                               snr=("snr_db", "median"), snr0=("snr_null_db", "median"),
                               ccsb=("cc_vs_seisbench", "median"),
                               spick=("s_pick_s", "median")).round(3).to_string())


if __name__ == "__main__":
    main()
