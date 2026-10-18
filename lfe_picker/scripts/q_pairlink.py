#!/usr/bin/env python
"""P3 STEP 1a — PAIR-LINKABILITY KILL TEST (template-free, multi-station, positive control).

Question: with Bostock's own 2012 event times handed to us as candidates (the best candidate set any
trigger could ever produce), can a template-free multi-station cross-correlation LINK the events of
one family to each other, on the 13-station / 7-site discovery network (S15 minus SNB and B926)?

Recipe (Merlin, P3 step 1a), verbatim:
  * candidates = every Bostock 2012 ETS-day event, kept AND hidden; family labels are used only
    afterwards, for scoring;
  * at each station cut a 24-s 3-C 20 Hz window starting at t_origin + c_sta, where c_sta is the
    Crosson S travel time from the GRID CENTRE (one constant per station - no family knowledge);
    t_origin is on our time base (Bostock DB time + frozen 2012 family offset, as labels_v0);
  * normalised 3-C cc for every candidate pair at every station over lags +-3 s;
  * MAD-normalise per station-day (here: per station x (day of i, day of j) block, the two-sided
    version of "per station-day" for a pair statistic);
  * merge the two co-located clusters into sites by the MEAN, sum over the 7 sites, max over lag
    -> the network pair cc;
  * NULL: identical computation with each station's continuous record circularly shifted by
    (station index x 1800 s) - destroys inter-station coherence, leaves every single-station
    statistic (including single-station repeaters) intact;
  * threshold = the null's 1e-6 quantile.

Run the NULL first (it sets the threshold), then the TRUE run:

    python scripts/q_pairlink.py --mode null --tag null
    python scripts/q_pairlink.py --mode true --tag true --thr-from results/p3/pairlink_null.json
"""
from __future__ import annotations
import os, sys, json, time, argparse
import numpy as np
import pandas as pd
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P

P3 = f"{P.RES}/p3"
FS = 20.0
WIN_S, LAG_S = 24.0, 3.0
NW = int(round(WIN_S * FS))          # 480
NL = int(round(LAG_S * FS))          # 60
NB = NW + 2 * NL                     # 600
NLAG = 2 * NL + 1                    # 121

WITHHELD = ("SNB", "B926")
STAS = [s for s in P.S15 if s not in WITHHELD]
SITE_OF = {s: P.SITE_OF[s] for s in STAS}
SITES = sorted(set(SITE_OF.values()))
NMEM = {q: sum(1 for s in STAS if SITE_OF[s] == q) for q in SITES}
GRID_CENTRE = (0.5 * (P.GRID_LAT[0] + P.GRID_LAT[1]), 0.5 * (P.GRID_LON[0] + P.GRID_LON[1]))
HIST_LO, HIST_HI, HIST_N = -30.0, 150.0, 3600


def candidates():
    B = P.bostock()
    B = B[B.year == 2012]
    ed = json.load(open(f"{P.RES}/p2/score_p2nnpu_2012.json"))["ets_days"]
    return B[B.doy.isin(ed)].sort_values(["doy", "t_origin"]).reset_index(drop=True)


def station_constants():
    import p_vel1d as V
    ST, CR = P.station_table(), P.station_corr()
    dep = float(np.atleast_1d(P.slab2_interp()(*GRID_CENTRE))[0])
    out = {}
    for s in STAS:
        d = float(P.epi_km(GRID_CENTRE[0], GRID_CENTRE[1], ST.lat[s], ST.lon[s]))
        out[s] = float(V.tt(d, dep, "S")) + CR[s]
    return out, dep


def cut_windows(B, cst, mode, store):
    """(N, nsta, 3, NB) float16 + validity (N, nsta) bool.

    Validity is decided on the TRUE window in both modes, so the null and the true run see exactly
    the same event set.  In null mode the DATA come from the circularly shifted record."""
    N, S = len(B), len(STAS)
    W = np.zeros((N, S, 3, NB), np.float16)
    ok = np.zeros((N, S), bool)
    for si, s in enumerate(STAS):
        shift = int(round(si * 1800.0 * FS)) if mode == "null" else 0
        for doy, g in B.groupby("doy"):
            if not store.has(s, 2012, int(doy)):
                continue
            k = store.idx(s, 2012, int(doy))
            n = int(store.nvalid[k])
            if n < 4 * NB:
                continue
            X = np.asarray(store.get_i(k)[:, :n])
            t0 = float(store.t0[k])
            i0 = np.rint((g.t_origin.to_numpy() + cst[s] - LAG_S - t0) * FS).astype(np.int64)
            for row, aa in zip(g.index.to_numpy(), i0):
                if aa < 0 or aa + NB > n:
                    continue
                true_seg = X[:, aa:aa + NB]
                if not np.square(true_seg.astype(np.float32)).sum() > 0:
                    continue
                seg = X[:, (np.arange(aa, aa + NB) + shift) % n] if shift else true_seg
                if not np.square(seg.astype(np.float32)).sum() > 0:
                    continue
                W[row, si] = seg
                ok[row, si] = True
        print(f"  [cut {mode}] {s}: {ok[:, si].sum()}/{N}", flush=True)
    return W, ok


def unfold_unit(Wd):
    """(n, 3, NB) fp16 GPU -> (n*NLAG, 3*NW) fp16, each row unit Frobenius norm."""
    n = Wd.shape[0]
    U = Wd.unfold(2, NW, 1).permute(0, 2, 1, 3).reshape(n * NLAG, 3 * NW).float()
    U = U / U.pow(2).sum(1, keepdim=True).sqrt().clamp_min(1e-12)
    return U.half()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["true", "null"], required=True)
    ap.add_argument("--tag", default="")
    ap.add_argument("--chunk", type=int, default=768)
    ap.add_argument("--thr-from", default="")
    ap.add_argument("--subsample", type=int, default=4_000_000)
    ap.add_argument("--topk", type=int, default=50000)
    ap.add_argument("--limit-days", type=int, default=0)
    ap.add_argument("--max-edges", type=int, default=60_000_000)
    ap.add_argument("--edge-above", type=float, default=-1e30)
    ap.add_argument("--min-dt", type=float, default=WIN_S + LAG_S,
                    help="pairs closer in time than this SHARE SAMPLES and are excluded "
                         "(their cc is trivially ~1 in the true run AND in the null)")
    a = ap.parse_args()
    tag = a.tag or a.mode
    os.makedirs(P3, exist_ok=True)
    dev = "cuda"
    t_start = time.time()

    B = candidates()
    if a.limit_days:
        B = B[B.doy.isin(sorted(B.doy.unique())[:a.limit_days])].reset_index(drop=True)
    cst, cdep = station_constants()
    print(f"[setup] {len(B)} candidates, {len(STAS)} stations, {len(SITES)} sites "
          f"{SITES}, grid centre {GRID_CENTRE} dep {cdep:.1f} km")
    print("[setup] c_sta = " + "  ".join(f"{s}:{cst[s]:.2f}" for s in STAS))

    store = P.DayStore()
    W, ok = cut_windows(B, cst, a.mode, store)
    use = ok.all(1)                                   # every station of the discovery net present
    print(f"[setup] usable events {use.sum()}/{len(B)} (all {len(STAS)} stations valid)")
    B = B[use].reset_index(drop=True)
    W, ok = W[use], ok[use]
    N = len(B)
    del store

    Wd = torch.from_numpy(W).to(dev)
    del W
    A = {}
    for si, s in enumerate(STAS):
        x = Wd[:, si, :, NL:NL + NW].reshape(N, 3 * NW).float()
        A[s] = (x / x.pow(2).sum(1, keepdim=True).sqrt().clamp_min(1e-12)).half()
        del x
    print(f"[gpu] windows {Wd.numel()*2/1e9:.2f} GB  A {13*N*3*NW*2/1e9:.2f} GB", flush=True)

    days = sorted(B.doy.unique())
    didx = {d: np.flatnonzero(B.doy.to_numpy() == d) for d in days}
    bkey = torch.from_numpy(B.bkey.to_numpy().astype(np.int64)).to(dev)
    torig = torch.from_numpy(B.t_origin.to_numpy()).to(dev)
    magd = torch.from_numpy(B.mag.to_numpy().astype(np.float32)).to(dev)

    FAM = P.families()
    fl = FAM[["lat", "lon"]].to_numpy()
    D5 = np.zeros((130, 130), bool)
    for i in range(len(fl)):
        D5[i] = P.epi_km(fl[:, 0], fl[:, 1], fl[i, 0], fl[i, 1]) < 5.0
    D5d = torch.from_numpy(D5).to(dev)

    thr = thr_lo = None
    THRN = ["q_0.001", "q_0.0001", "q_1e-5", "q_1e-6"]
    THRV = None
    if a.thr_from:
        z = json.load(open(a.thr_from))
        thr, thr_lo = float(z["q_1e-6"]), float(z["q_1e-5"])
        THRV = torch.tensor([float(z[k]) for k in THRN], device=dev)
        print(f"[thr] null quantiles " + "  ".join(f"{k}={float(z[k]):.2f}" for k in THRN))

    hist = torch.zeros(HIST_N, dtype=torch.float64, device=dev)
    top = torch.full((0,), -1e30, device=dev)
    n_pairs = 0
    samp = []
    NT = 4
    MAGB = torch.tensor([-9.0, 1.4, 1.7, 2.0, 9.0], device=dev)
    mag_same = torch.zeros(NT + 1, 4, dtype=torch.int64, device=dev)
    mag_far = torch.zeros(NT + 1, 4, dtype=torch.int64, device=dev)
    same_hitT = torch.zeros(NT, 130, dtype=torch.int64, device=dev)
    ali_hitT = torch.zeros(NT, dtype=torch.int64, device=dev)
    far_hitT = torch.zeros(NT, dtype=torch.int64, device=dev)
    same_hit = torch.zeros(130, dtype=torch.int64, device=dev)
    same_tot = torch.zeros(130, dtype=torch.int64, device=dev)
    same_sum = torch.zeros(130, dtype=torch.float64, device=dev)
    same_max = torch.full((130,), -1e30, device=dev)
    hist_same = torch.zeros(HIST_N, dtype=torch.float64, device=dev)
    hist_near = torch.zeros(HIST_N, dtype=torch.float64, device=dev)
    hist_far = torch.zeros(HIST_N, dtype=torch.float64, device=dev)
    ali_hit = torch.zeros((), dtype=torch.int64, device=dev)
    ali_tot = torch.zeros((), dtype=torch.int64, device=dev)
    far_hit = torch.zeros((), dtype=torch.int64, device=dev)
    far_tot = torch.zeros((), dtype=torch.int64, device=dev)
    edges, n_edge = [], 0
    rng = torch.Generator(device=dev); rng.manual_seed(7)

    for dj in days:
        jjn = didx[dj]
        jj = torch.from_numpy(jjn).to(dev)
        nj = len(jjn)
        Bl = {s: unfold_unit(Wd[jj, si]) for si, s in enumerate(STAS)}
        for di in [d for d in days if d <= dj]:
            iin = didx[di]
            ns = min(256, len(iin))
            isub = torch.from_numpy(iin[np.linspace(0, len(iin) - 1, ns).astype(int)]).to(dev)
            med, mad = {}, {}
            for s in STAS:
                v = torch.mm(A[s][isub], Bl[s].T).float().reshape(-1)
                v = v[torch.randint(0, v.numel(), (min(400_000, v.numel()),),
                                    device=dev, generator=rng)]
                m = v.median()
                med[s] = m
                mad[s] = ((v - m).abs().median().clamp_min(1e-6) * 1.4826)
                del v
            for c0 in range(0, len(iin), a.chunk):
                ic = torch.from_numpy(iin[c0:c0 + a.chunk]).to(dev)
                ni = len(ic)
                acc = torch.zeros(ni, nj * NLAG, device=dev)
                for s in STAS:
                    r = torch.mm(A[s][ic], Bl[s].T).float()
                    acc += (r - med[s]) / (mad[s] * NMEM[SITE_OF[s]])
                    del r
                acc = acc.reshape(ni, nj, NLAG).amax(-1)
                gi = ic[:, None].expand(ni, nj)
                gj = jj[None, :].expand(ni, nj)
                keep = (gi < gj) & ((torig[gi] - torig[gj]).abs() >= a.min_dt)
                vals = acc[keep]
                del acc
                if vals.numel() == 0:
                    continue
                n_pairs += int(vals.numel())
                hist += torch.histc(vals, HIST_N, HIST_LO, HIST_HI).double()
                if len(samp) * 200_000 < a.subsample:
                    k = min(200_000, vals.numel())
                    samp.append(vals[torch.randint(0, vals.numel(), (k,), device=dev,
                                                   generator=rng)].cpu().numpy())
                if a.mode == "null":
                    top = torch.topk(torch.cat([top, vals]),
                                     min(a.topk, top.numel() + vals.numel()))[0]
                    if a.edge_above > -1e29:
                        em = vals >= a.edge_above
                        if em.any():
                            edges.append(np.stack([gi[keep][em].to(torch.int32).cpu().numpy(),
                                                   gj[keep][em].to(torch.int32).cpu().numpy(),
                                                   vals[em].cpu().numpy()], 1))
                            n_edge += int(em.sum())
                else:
                    gik, gjk = gi[keep], gj[keep]
                    bi, bj = bkey[gik], bkey[gjk]
                    mmin = torch.minimum(magd[gik], magd[gjk])
                    mb = torch.bucketize(mmin, MAGB[1:-1])
                    sf = bi == bj
                    near = D5d[bi, bj] & ~sf
                    hit = vals >= thr
                    same_hit += torch.bincount(bi[sf & hit], minlength=130)
                    same_tot += torch.bincount(bi[sf], minlength=130)
                    same_sum += torch.bincount(bi[sf], weights=vals[sf].double(), minlength=130)
                    same_max = torch.maximum(
                        same_max, torch.full((130,), -1e30, device=dev).scatter_reduce(
                            0, bi[sf], vals[sf], reduce="amax"))
                    hist_same += torch.histc(vals[sf], HIST_N, HIST_LO, HIST_HI).double()
                    hist_near += torch.histc(vals[near], HIST_N, HIST_LO, HIST_HI).double()
                    ff = (~near) & (~sf)
                    hist_far += torch.histc(vals[ff], HIST_N, HIST_LO, HIST_HI).double()
                    ali_hit += (near & hit).sum(); ali_tot += near.sum()
                    far_hit += ((~near) & (~sf) & hit).sum(); far_tot += ((~near) & (~sf)).sum()
                    ffm = (~near) & (~sf)
                    for ti in range(NT):
                        ht = vals >= THRV[ti]
                        same_hitT[ti] += torch.bincount(bi[sf & ht], minlength=130)
                        ali_hitT[ti] += (near & ht).sum()
                        far_hitT[ti] += (ffm & ht).sum()
                        mag_same[ti] += torch.bincount(mb[sf & ht], minlength=4)
                        mag_far[ti] += torch.bincount(mb[ffm & ht], minlength=4)
                    mag_same[NT] += torch.bincount(mb[sf], minlength=4)
                    mag_far[NT] += torch.bincount(mb[ffm], minlength=4)
                    em = vals >= thr_lo
                    ne = int(em.sum())
                    if ne and n_edge < a.max_edges:
                        edges.append(np.stack([gi[keep][em].to(torch.int32).cpu().numpy(),
                                               gj[keep][em].to(torch.int32).cpu().numpy(),
                                               vals[em].cpu().numpy()], 1))
                        n_edge += ne
                    del bi, bj, sf, near, hit, em
                del vals, keep, gi, gj
            del med, mad
        del Bl
        torch.cuda.empty_cache()
        print(f"[day {dj}] pairs {n_pairs/1e6:.1f}M  edges {n_edge/1e6:.2f}M  "
              f"({time.time()-t_start:.0f}s)", flush=True)

    out = dict(mode=a.mode, tag=tag, n_events=int(N), n_pairs=int(n_pairs),
               stations=STAS, sites=SITES, nmem=NMEM, withheld=list(WITHHELD),
               c_sta=cst, grid_centre=GRID_CENTRE, grid_dep=cdep,
               win_s=WIN_S, lag_s=LAG_S, min_dt=a.min_dt, hist_lo=HIST_LO, hist_hi=HIST_HI, hist_n=HIST_N,
               runtime_s=time.time() - t_start)
    np.savez(f"{P3}/pairlink_{tag}_hist.npz", hist=hist.cpu().numpy(),
             lo=HIST_LO, hi=HIST_HI, samp=np.concatenate(samp) if samp else np.zeros(0))
    if a.mode == "null":
        t = torch.sort(top, descending=True)[0].cpu().numpy()
        for q, name in [(1e-6, "q_1e-6"), (1e-5, "q_1e-5"), (1e-4, "q_1e-4")]:
            r = int(max(1, round(q * n_pairs)))
            out[name] = float(t[min(r, len(t)) - 1])
        out["null_max"] = float(t[0]); out["topk_stored"] = int(len(t))
        np.save(f"{P3}/pairlink_{tag}_top.npy", t)
        if edges:
            E = np.concatenate(edges)
            np.savez(f"{P3}/pairlink_{tag}_edges.npz", i=E[:, 0].astype(np.int32),
                     j=E[:, 1].astype(np.int32), v=E[:, 2].astype(np.float32))
            out["n_edges_stored"] = int(n_edge)
    else:
        out.update(threshold=thr, edge_floor=thr_lo,
                   alias_hit=int(ali_hit), alias_tot=int(ali_tot),
                   alias_rate=float(ali_hit) / max(1, int(ali_tot)),
                   far_hit=int(far_hit), far_tot=int(far_tot),
                   far_rate=float(far_hit) / max(1, int(far_tot)),
                   n_edges_stored=int(n_edge))
        st = same_tot.cpu().numpy()
        pd.DataFrame(dict(bkey=np.arange(130), n_same_pairs=st,
                          n_same_above=same_hit.cpu().numpy(),
                          mean_same=same_sum.cpu().numpy() / np.maximum(1, st),
                          max_same=same_max.cpu().numpy())).to_csv(
            f"{P3}/pairlink_{tag}_by_family.csv", index=False)
        pd.DataFrame(same_hitT.cpu().numpy().T, columns=THRN).assign(
            bkey=np.arange(130), n_same_pairs=st).to_csv(
            f"{P3}/pairlink_{tag}_by_family_thresholds.csv", index=False)
        out["alias_hit_by_thr"] = {k: int(v) for k, v in zip(THRN, ali_hitT.cpu().numpy())}
        out["far_hit_by_thr"] = {k: int(v) for k, v in zip(THRN, far_hitT.cpu().numpy())}
        out["same_hit_by_thr"] = {k: int(v) for k, v in zip(THRN, same_hitT.sum(1).cpu().numpy())}
        MS, MF = mag_same.cpu().numpy(), mag_far.cpu().numpy()
        pd.DataFrame(dict(magbin=["<1.4", "1.4-1.7", "1.7-2.0", ">2.0"],
                          **{f"same_{k}": MS[i] for i, k in enumerate(THRN)},
                          same_total=MS[NT],
                          **{f"far_{k}": MF[i] for i, k in enumerate(THRN)},
                          far_total=MF[NT])).to_csv(
            f"{P3}/pairlink_{tag}_by_magnitude.csv", index=False)
        np.savez(f"{P3}/pairlink_{tag}_hist_class.npz", same=hist_same.cpu().numpy(),
                 near=hist_near.cpu().numpy(), far=hist_far.cpu().numpy(),
                 lo=HIST_LO, hi=HIST_HI)
        if edges:
            E = np.concatenate(edges)
            np.savez(f"{P3}/pairlink_{tag}_edges.npz", i=E[:, 0].astype(np.int32),
                     j=E[:, 1].astype(np.int32), v=E[:, 2].astype(np.float32))
    B[["bkey", "doy", "t_origin", "hidden", "mag"]].to_parquet(
        f"{P3}/pairlink_{tag}_events.parquet", index=False)
    json.dump(out, open(f"{P3}/pairlink_{tag}.json", "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if k != "c_sta"}, indent=1, default=str))


if __name__ == "__main__":
    main()
