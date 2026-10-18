#!/usr/bin/env python
"""P3D STEP 3 -- pair-rate POSITIVE CONTROL and the discovery PREDICTION (GPU).

The P3 Step-1a statistic (`q_pairlink2.py`), re-run on the 2019-2026 dense network:
  candidates  hp events (Track E p_real >= 0.8) of the KEPT (non-hidden) families on the SSE-window
              days of the CALIBRATION episodes S01-S13
  statistic   per-site MAD-z of the 24-s 3-C normalised cross-correlation, max over +-3 s lag,
              summed over the sites of a configuration; 27-s overlap exclusion
  null        the SAME chain with station si's data circularly shifted by si x 1800 s
  configs     (i) old7   = the seven P3 discovery sites, frozen stations only
              (ii) dense = the a priori geometry rule (median hypocentral distance to the family
                   field <= 60 km), withheld sites dropped
              plus the single-site and drop-one ablations that say where the change comes from
  fit         delta (the same-family shift in units of the null spread) from the tail quantiles,
              exactly the P3 RESULTS 5.5 method, for both configurations

  python scripts/s_pairlink3.py --mode null --tag cal_null
  python scripts/s_pairlink3.py --mode true --tag cal_true --thr-from results/p3d/pairlink3_cal_null.json
"""
from __future__ import annotations
import os, sys, json, time, argparse
import numpy as np
import pandas as pd
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                        # noqa: E402

OUT = S.P3D
FS = 20.0
WIN_S, LAG_S = 24.0, 3.0
NW = int(round(WIN_S * FS))
NL = int(round(LAG_S * FS))
NB = NW + 2 * NL
NLAG = 2 * NL + 1
THRN = ["q_0.001", "q_0.0001", "q_1e-5", "q_1e-6"]
YEARS = list(range(2019, 2027))

# frozen stations of the seven P3 discovery sites (co-located NEW stations deliberately excluded
# from `old7`, so the two configurations differ only by the genuinely new geometry)
OLD7_STAS = ["PGC", "B009", "B010", "B011", "VGZ", "NLLB",
             "LRIV", "B005", "B006", "B007", "B001", "B003", "B004"]


def config_stations(ST, dense_sites):
    """station lists for the two pre-registered configurations, and the union to cut."""
    old7 = [s for s in OLD7_STAS]
    dense = [s for s in ST.index if ST.site[s] in dense_sites and ST.summable[s]
             and not ST.withheld[s]]
    return old7, dense, sorted(set(old7) | set(dense))


def candidates(episodes, stas, idx, hidden):
    H = S.hp_events(0.8, episodes=episodes)
    H = H[~H.bkey.isin(hidden)].copy()
    yr, doy = S.yd_of(H.t.values)
    H["year"], H["doy"] = yr, doy
    have = {(y, d) for (s, y, d) in idx if s == stas[0]}
    for s in stas[1:]:
        have &= {(y, d) for (ss, y, d) in idx if ss == s}
    H = H[[(y, d) in have for y, d in zip(H.year, H.doy)]]
    return H.sort_values("t").reset_index(drop=True), len(have)


def station_constants(stas, ST):
    import p_vel1d as V
    import p_common as P
    ctr = (0.5 * (P.GRID_LAT[0] + P.GRID_LAT[1]), 0.5 * (P.GRID_LON[0] + P.GRID_LON[1]))
    dep = float(np.atleast_1d(P.slab2_interp()(*ctr))[0])
    CR = S.station_corr(stas)
    return {s: float(V.tt(float(S.epi_km(ctr[0], ctr[1], ST.lat[s], ST.lon[s])), dep, "S")) + CR[s]
            for s in stas}, ctr, dep


def cut_windows(B, cst, mode, dc, stas, null_shift=1800.0):
    N, NS = len(B), len(stas)
    W = np.zeros((N, NS, 3, NB), np.float16)
    ok = np.zeros((N, NS), bool)
    for si, s in enumerate(stas):
        shift = int(round(si * null_shift * FS)) if mode == "null" else 0
        for (y, doy), g in B.groupby(["year", "doy"]):
            v = dc.get(s, int(y), int(doy))
            if v is None:
                continue
            t0, X = v
            n = X.shape[1]
            if n < 4 * NB:
                continue
            i0 = np.rint((g.t.to_numpy() + cst[s] - LAG_S - t0) * FS).astype(np.int64)
            for row, aa in zip(g.index.to_numpy(), i0):
                if aa < 0 or aa + NB > n:
                    continue
                if not np.square(X[:, aa:aa + NB], dtype=np.float64).sum() > 0:
                    continue
                seg = X[:, (np.arange(aa, aa + NB) + shift) % n] if shift else X[:, aa:aa + NB]
                e = float(np.square(seg, dtype=np.float64).sum())
                if not (np.isfinite(e) and e > 0):
                    continue
                W[row, si] = (seg / np.sqrt(e / seg.size)).astype(np.float16)
                ok[row, si] = bool(np.isfinite(W[row, si]).all())
        print(f"  [cut {mode}] {s}: {ok[:, si].sum()}/{N}", flush=True)
    return W, ok


def unfold_unit(Wd):
    n = Wd.shape[0]
    U = Wd.unfold(2, NW, 1).permute(0, 2, 1, 3).reshape(n * NLAG, 3 * NW).float()
    U = U / U.pow(2).sum(1, keepdim=True).sqrt().clamp_min(1e-12)
    return U.half()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["true", "null"], required=True)
    ap.add_argument("--tag", default="")
    ap.add_argument("--episodes", default=",".join(S.CAL_EPISODES))
    ap.add_argument("--chunk", type=int, default=512)
    ap.add_argument("--thr-from", default="")
    ap.add_argument("--topk", type=int, default=200000)
    ap.add_argument("--min-dt", type=float, default=WIN_S + LAG_S)
    ap.add_argument("--variant", default="v1", choices=["v1", "v2"])
    ap.add_argument("--null-shift", type=float, default=1800.0,
                    help="the index x SHIFT seconds circular shift of the null; several values give "
                         "independent null REALISATIONS of the same statistic")
    ap.add_argument("--min-avail", type=float, default=0.3,
                    help="drop a station delivering fewer than this fraction of the episode days; "
                         "the candidate set is the INTERSECTION of the stations' cached days, so a "
                         "thin station would gut it")
    a = ap.parse_args()
    tag = a.tag or f"cal_{a.mode}"
    eps = a.episodes.split(",")
    t_start = time.time()
    ST = (S.station_table() if a.variant == "v1"
          else pd.read_csv(f"{S.P3D}/stations_p3d_v2.csv").set_index("sta"))
    SI = S.site_table(ST)
    hidden, iso, _ = S.holdouts()
    dense_sites = S.dense_sites(SI)
    old7, dense, stas = config_stations(ST, dense_sites)
    site_of = {s: ST.site[s] for s in stas}
    sites = sorted(set(site_of.values()))
    nmem = {q: sum(1 for s in stas if site_of[s] == q) for q in sites}
    idx = S.cache_index(years=YEARS, stas=stas)
    # a station that does not exist over these episodes cannot be in either configuration: the
    # candidate set is the intersection of the stations' cached days, so one absent station would
    # empty it.  BUCKS/TURTL/ODUC are post-2023 and drop out of the CALIBRATION run by this rule.
    W = S.sse_windows()
    epdays = {(int(t.year), int(t.dayofyear)) for r in W[W.id.isin(eps)].itertuples()
              for t in pd.date_range(r.start, r.end, freq="D")}
    navail = {s: sum(1 for (ss, y, d) in idx if ss == s and (y, d) in epdays) for s in stas}
    drop = [s for s in stas if navail[s] < a.min_avail * len(epdays)]
    if drop:
        print(f"[setup] dropped (no data over these episodes): "
              f"{ {s: navail[s] for s in drop} } of {len(epdays)} episode days", flush=True)
    stas = [s for s in stas if s not in drop]
    old7 = [s for s in old7 if s not in drop]
    dense = [s for s in dense if s not in drop]
    site_of = {s: ST.site[s] for s in stas}
    sites = sorted(set(site_of.values()))
    nmem = {q: sum(1 for s in stas if site_of[s] == q) for q in sites}
    B, nday = candidates(eps, stas, idx, hidden)
    cst, ctr, cdep = station_constants(stas, ST)
    print(f"[setup] episodes {eps}: {len(B)} hp candidates of kept families on {nday} days with "
          f"all {len(stas)} stations cached", flush=True)
    print(f"[setup] old7 sites {sorted({site_of[s] for s in old7})}; "
          f"dense (<=60 km rule) sites {sorted({site_of[s] for s in dense})}", flush=True)

    CFG = {"old7": sorted({site_of[s] for s in old7}),
           "dense": sorted({site_of[s] for s in dense})}
    for q in sites:                                    # single-site ablations
        CFG[f"site_{q}"] = [q]
    CFG["dense_minus_SYMB"] = [q for q in CFG["dense"] if q != "SYMB"]
    CFG["old7_plus_SYMB"] = sorted(set(CFG["old7"]) | ({"SYMB"} if "SYMB" in sites else set()))
    CFG = {k: v for k, v in CFG.items() if v}
    cnames = list(CFG)

    dc = S.DayCache(maxn=4, years=YEARS, stas=stas)
    W, ok = cut_windows(B, cst, a.mode, dc, stas, a.null_shift)
    use = ok.all(1)
    B = B[use].reset_index(drop=True)
    W = W[use]
    N = len(B)
    print(f"[setup] usable {N}/{len(use)}", flush=True)
    del dc

    dev = "cuda"
    Wd = torch.from_numpy(W).to(dev)
    del W
    A = {}
    for si, s in enumerate(stas):
        x = Wd[:, si, :, NL:NL + NW].reshape(N, 3 * NW).float()
        A[s] = (x / x.pow(2).sum(1, keepdim=True).sqrt().clamp_min(1e-12)).half()
        del x
    print(f"[gpu] windows {Wd.numel()*2/1e9:.2f} GB", flush=True)

    ka = np.array([f"{y}.{d}" for y, d in zip(B.year, B.doy)])
    days = sorted(set(zip(B.year.to_numpy(), B.doy.to_numpy())))
    didx = {k: np.flatnonzero(ka == f"{k[0]}.{k[1]}") for k in days}
    bkey = torch.from_numpy(B.bkey.to_numpy().astype(np.int64)).to(dev)
    tt = torch.from_numpy(B.t.to_numpy()).to(dev)
    F = S.families()
    fl = F[["lat", "lon"]].to_numpy()
    D5 = np.zeros((130, 130), bool)
    for i in range(len(fl)):
        D5[i] = S.epi_km(fl[:, 0], fl[:, 1], fl[i, 0], fl[i, 1]) < 5.0
    D5d = torch.from_numpy(D5).to(dev)

    THRV = None
    if a.thr_from:
        z = json.load(open(a.thr_from))
        THRV = {c: torch.tensor([float(z["thresholds"][c][k]) for k in THRN], device=dev)
                for c in cnames}

    nc = len(cnames)
    tops = {c: torch.full((0,), -1e30, device=dev) for c in cnames}
    same_hit = torch.zeros(nc, 4, 130, dtype=torch.int64, device=dev)
    same_tot = torch.zeros(130, dtype=torch.int64, device=dev)
    ali_hit = torch.zeros(nc, 4, dtype=torch.int64, device=dev)
    far_hit = torch.zeros(nc, 4, dtype=torch.int64, device=dev)
    ali_tot = torch.zeros((), dtype=torch.int64, device=dev)
    far_tot = torch.zeros((), dtype=torch.int64, device=dev)
    n_pairs = 0
    HIST_LO, HIST_HI, HIST_N = -30.0, 150.0, 3600
    hist = {c: {k: torch.zeros(HIST_N, dtype=torch.float64, device=dev)
                for k in ("all", "same", "near", "far")} for c in ("old7", "dense")}
    rng = torch.Generator(device=dev)
    rng.manual_seed(7)

    for dj in days:
        jjn = didx[dj]
        jj = torch.from_numpy(jjn).to(dev)
        nj = len(jjn)
        Bl = {s: unfold_unit(Wd[jj, si]) for si, s in enumerate(stas)}
        for di in [d for d in days if d <= dj]:
            iin = didx[di]
            ns = min(256, len(iin))
            isub = torch.from_numpy(iin[np.linspace(0, len(iin) - 1, ns).astype(int)]).to(dev)
            med, mad = {}, {}
            for s in stas:
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
                gi = ic[:, None].expand(ni, nj)
                gj = jj[None, :].expand(ni, nj)
                keep = (gi < gj) & ((tt[gi] - tt[gj]).abs() >= a.min_dt)
                if not bool(keep.any()):
                    continue
                Z = {q: torch.zeros(ni, nj * NLAG, device=dev) for q in sites}
                for s in stas:
                    r = torch.mm(A[s][ic], Bl[s].T).float()
                    Z[site_of[s]] += (r - med[s]) / (mad[s] * nmem[site_of[s]])
                    del r
                gik, gjk = gi[keep], gj[keep]
                bi, bj = bkey[gik], bkey[gjk]
                sf = bi == bj
                near = D5d[bi, bj] & ~sf
                ff = (~near) & (~sf)
                if a.mode == "true":
                    same_tot += torch.bincount(bi[sf], minlength=130)
                    ali_tot += near.sum()
                    far_tot += ff.sum()
                n_pairs += int(keep.sum())
                for ci, c in enumerate(cnames):
                    sl = CFG[c]
                    if not sl:
                        continue
                    acc = Z[sl[0]].clone()
                    for q in sl[1:]:
                        acc += Z[q]
                    vals = acc.reshape(ni, nj, NLAG).amax(-1)[keep]
                    del acc
                    if c in hist:
                        hist[c]["all"] += torch.histc(vals, HIST_N, HIST_LO, HIST_HI).double()
                        if a.mode == "true":
                            hist[c]["same"] += torch.histc(vals[sf], HIST_N, HIST_LO, HIST_HI).double()
                            hist[c]["near"] += torch.histc(vals[near], HIST_N, HIST_LO, HIST_HI).double()
                            hist[c]["far"] += torch.histc(vals[ff], HIST_N, HIST_LO, HIST_HI).double()
                    if a.mode == "null":
                        tops[c] = torch.topk(torch.cat([tops[c], vals]),
                                             min(a.topk, tops[c].numel() + vals.numel()))[0]
                    else:
                        for ti in range(4):
                            ht = vals >= THRV[c][ti]
                            same_hit[ci, ti] += torch.bincount(bi[sf & ht], minlength=130)
                            ali_hit[ci, ti] += (near & ht).sum()
                            far_hit[ci, ti] += (ff & ht).sum()
                    del vals
                del Z, gi, gj, keep, bi, bj, sf, near, ff
        del Bl
        torch.cuda.empty_cache()
    print(f"[pairs] {n_pairs/1e6:.1f}M  ({time.time()-t_start:.0f}s)", flush=True)

    out = dict(mode=a.mode, tag=tag, episodes=eps, n_events=int(N), n_pairs=int(n_pairs),
               stations=stas, sites=sites, nmem=nmem, min_dt=a.min_dt,
               variant=a.variant, null_shift=a.null_shift,
               configs={c: CFG[c] for c in cnames}, c_sta=cst, grid_centre=ctr, grid_dep=cdep,
               runtime_s=time.time() - t_start)
    np.savez(f"{OUT}/pairlink3_{tag}_hist.npz", lo=HIST_LO, hi=HIST_HI,
             **{f"{c}_{k}": v.cpu().numpy() for c, d in hist.items() for k, v in d.items()})
    B[["bkey", "year", "doy", "t", "p_real"]].to_parquet(
        f"{OUT}/pairlink3_{tag}_events.parquet", index=False)
    if a.mode == "null":
        thr = {}
        for c in cnames:
            t = torch.sort(tops[c], descending=True)[0].cpu().numpy()
            thr[c] = {}
            for q, name in [(1e-3, "q_0.001"), (1e-4, "q_0.0001"), (1e-5, "q_1e-5"),
                            (1e-6, "q_1e-6")]:
                r = int(max(1, round(q * n_pairs)))
                thr[c][name] = float(t[min(r, len(t)) - 1])
            thr[c]["max"] = float(t[0])
            thr[c]["n_kept"] = int(len(t))
        out["thresholds"] = thr
        print(json.dumps({k: thr[k] for k in ("old7", "dense")}, indent=1))
    else:
        rows = []
        st = same_tot.cpu().numpy()
        SH = same_hit.cpu().numpy()
        for ci, c in enumerate(cnames):
            for b in range(130):
                if st[b] == 0:
                    continue
                rows.append(dict(config=c, bkey=b, n_same_pairs=int(st[b]),
                                 **{k: int(SH[ci, ti, b]) for ti, k in enumerate(THRN)}))
        pd.DataFrame(rows).to_csv(f"{OUT}/pairlink3_{tag}_by_family.csv", index=False)
        out["alias_tot"] = int(ali_tot)
        out["far_tot"] = int(far_tot)
        out["same_tot"] = int(st.sum())
        out["alias_hit"] = {c: {k: int(v) for k, v in zip(THRN, ali_hit[ci].cpu().numpy())}
                            for ci, c in enumerate(cnames)}
        out["far_hit"] = {c: {k: int(v) for k, v in zip(THRN, far_hit[ci].cpu().numpy())}
                          for ci, c in enumerate(cnames)}
        out["same_hit"] = {c: {k: int(SH[ci, ti].sum()) for ti, k in enumerate(THRN)}
                           for ci, c in enumerate(cnames)}
    json.dump(out, open(f"{OUT}/pairlink3_{tag}.json", "w"), indent=1, default=str)
    print(f"[done] {time.time()-t_start:.0f}s -> {OUT}/pairlink3_{tag}.json")


if __name__ == "__main__":
    main()
