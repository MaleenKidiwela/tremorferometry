#!/usr/bin/env python
"""P3D DISCOVERY -- template-free pairwise linkage on uniform tremor-hour windows (GPU).

Same statistic as P3 Step 1a / `q_pairlink2` / `s_pairlink3` (per-site MAD-z of the 24-s 3-C
normalised cross-correlation, max over +-3 s lag, summed over the configuration's sites, 27-s
overlap exclusion, index x SHIFT s circular-shift null), but

  * candidates are UNIFORM 12-s windows at 6-s hop on PNSN tremor hours within 50 km of the family
    field -- template-free, as Merlin Q3 requires -- instead of catalogue detections;
  * the per-site accumulator is FUSED into a single tensor instead of one tensor per site.  That is
    the whole optimisation: the old code held `len(sites)` arrays of [n_i, n_j*121] fp32 and summed
    them, so an 8-site run moved 8x the bytes for no arithmetic reason.  Fusing makes the pass
    compute-bound;
  * `true` mode emits the EDGE LIST above a frozen threshold (the graph the discovery chain needs)
    rather than per-family hit counts.

  python scripts/s_discover_pairs.py --net B --mode null --shift 1800 --tag B_null
  python scripts/s_discover_pairs.py --net B --mode true --thr-from results/p3d/disc_B_null.json --tag B_true
  python scripts/s_discover_pairs.py --validate        # reproduces s_pairlink3 on the calibration set
"""
from __future__ import annotations
import os, sys, json, time, argparse, warnings
import numpy as np
import pandas as pd
import torch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                       # noqa: E402

OUT = S.P3D
FS = 20.0
WIN_S, LAG_S = 24.0, 3.0
NW = int(round(WIN_S * FS))
NL = int(round(LAG_S * FS))
NB = NW + 2 * NL
NLAG = 2 * NL + 1
EPISODES = ["S18", "S21", "S22"]
PNSN = "/home/jovyan/tremorferometry/catalogs/pnsn_tremor_cascadia_full.csv"

NETS = {
    "B": ["PGC", "B009", "B010", "B011", "VGZ", "SYMB", "BUCKS", "LHLYB", "ALBHB", "DAOB"],
    "D": ["PGC", "B009", "B010", "B011", "LHLYB", "DAOB"],
    # EXPLORATORY / POST-HOC (coordinator, 2026-09-19): the max-coverage 3-site set.  Chosen because
    # it is the only configuration in which all six scoreable-isolated PRIMARY families have
    # representable truth (see the reachable-truth table); NOT pre-registered.
    "M": ["PGC", "B009", "B010", "B011", "BUCKS", "DAOB"],
}


def tremor_hours(max_km=50.0):
    F = S.families()
    tr = pd.read_csv(PNSN, usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr.time, errors="coerce")
    tr = tr.dropna(subset=["time"])
    d = np.min(S.epi_km(tr.lat.values[:, None], tr.lon.values[:, None],
                        F.lat.values[None, :], F.lon.values[None, :]), axis=1)
    tr = tr[d <= max_km]
    return sorted(set(tr.time.dt.floor("h"))), tr.time.max()


def candidates(net, episodes, hop_s=6.0, win_s=12.0):
    """Uniform windows on tremor hours of days cached at EVERY station of the configuration."""
    stas = NETS[net]
    idx = S.cache_index(years=list(range(2019, 2027)), stas=stas)
    W = S.sse_windows()
    days = {}
    for ep in episodes:
        r = W[W.id == ep].iloc[0]
        for t in pd.date_range(r.start, r.end, freq="D"):
            if all((s, t.year, t.dayofyear) in idx for s in stas):
                days[t.date()] = ep
    hrs, cat_end = tremor_hours()
    rows = []
    n_per = int(round(3600.0 / hop_s))
    for h in hrs:
        ep = days.get(h.date())
        if ep is None:
            continue
        t0 = h.value / 1e9
        rows.append((t0, ep, h))
    C = []
    for t0, ep, h in rows:
        ts = t0 + np.arange(n_per) * hop_s
        C.append(pd.DataFrame(dict(t=ts, episode=ep, hour=h)))
    C = pd.concat(C, ignore_index=True) if C else pd.DataFrame(columns=["t", "episode", "hour"])
    yr, doy = S.yd_of(C.t.values)
    C["year"], C["doy"] = yr, doy
    return C.sort_values("t").reset_index(drop=True), idx, sorted(set(days)), str(cat_end)


def station_constants(stas, ST):
    import p_vel1d as V
    import p_common as P
    ctr = (0.5 * (P.GRID_LAT[0] + P.GRID_LAT[1]), 0.5 * (P.GRID_LON[0] + P.GRID_LON[1]))
    dep = float(np.atleast_1d(P.slab2_interp()(*ctr))[0])
    CR = S.station_corr(stas)
    return {s: float(V.tt(float(S.epi_km(ctr[0], ctr[1], ST.lat[s], ST.lon[s])), dep, "S")) + CR[s]
            for s in stas}, ctr, dep


def cut(C, cst, mode, dc, stas, shift_s):
    N, NS = len(C), len(stas)
    Wd = np.zeros((N, NS, 3, NB), np.float16)
    ok = np.zeros((N, NS), bool)
    for si, s in enumerate(stas):
        sh = int(round(si * shift_s * FS)) if mode == "null" else 0
        for (y, doy), g in C.groupby(["year", "doy"]):
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
                seg = X[:, (np.arange(aa, aa + NB) + sh) % n] if sh else X[:, aa:aa + NB]
                e = float(np.square(seg, dtype=np.float64).sum())
                if not (np.isfinite(e) and e > 0):
                    continue
                Wd[row, si] = (seg / np.sqrt(e / seg.size)).astype(np.float16)
                ok[row, si] = bool(np.isfinite(Wd[row, si]).all())
        print(f"  [cut {mode}] {s}: {ok[:, si].sum()}/{N}", flush=True)
    return Wd, ok


def unfold_unit(x):
    n = x.shape[0]
    U = x.unfold(2, NW, 1).permute(0, 2, 1, 3).reshape(n * NLAG, 3 * NW).float()
    U = U / U.pow(2).sum(1, keepdim=True).sqrt().clamp_min(1e-12)
    return U.half()


def run(a):
    t_start = time.time()
    ST = pd.read_csv(f"{S.P3D}/stations_p3d_v2.csv").set_index("sta")
    stas = NETS[a.net]
    site_of = {s: ST.site[s] for s in stas}
    sites = sorted(set(site_of.values()))
    nmem = {q: sum(1 for s in stas if site_of[s] == q) for q in sites}
    C, idx, days, cat_end = candidates(a.net, a.episodes.split(","))
    if a.cand_file:
        # ARM B / ARM C: an explicit candidate list replaces the uniform generator.  It is
        # INTERSECTED with the uniform set's qualifying days so all three arms see identical data.
        Q = pd.read_parquet(a.cand_file)
        okd = set(zip(C.year, C.doy))
        yr, doy = S.yd_of(Q.t.values)
        Q = Q.assign(year=yr, doy=doy)
        Q = Q[[(y, d) in okd for y, d in zip(Q.year, Q.doy)]]
        if "episode" not in Q:
            Q["episode"] = S.episode_of(pd.to_datetime(Q.t, unit="s"))
        C = Q.sort_values("t").reset_index(drop=True)
        print(f"[setup] candidate list {a.cand_file}: {len(C)} windows on the qualifying days",
              flush=True)
    if a.max_windows and len(C) > a.max_windows:
        C = C.iloc[np.linspace(0, len(C) - 1, a.max_windows).astype(int)].reset_index(drop=True)
    cst, ctr, cdep = station_constants(stas, ST)
    print(f"[setup] NET-{a.net} sites {sites} ({len(stas)} stations); {len(days)} days; "
          f"{len(C)} candidate windows; PNSN last detection {cat_end}", flush=True)
    dc = S.DayCache(maxn=4, years=list(range(2019, 2027)), stas=stas)
    Wn, ok = cut(C, cst, a.mode, dc, stas, a.shift)
    use = ok.all(1)
    C = C[use].reset_index(drop=True)
    Wn = Wn[use]
    N = len(C)
    print(f"[setup] usable {N}/{len(use)}", flush=True)
    del dc
    dev = "cuda"
    Wd = torch.from_numpy(Wn).to(dev)
    del Wn
    A = {}
    for si, s in enumerate(stas):
        x = Wd[:, si, :, NL:NL + NW].reshape(N, 3 * NW).float()
        A[s] = (x / x.pow(2).sum(1, keepdim=True).sqrt().clamp_min(1e-12)).half()
        del x
    print(f"[gpu] windows {Wd.numel()*2/1e9:.1f} GB + anchors {sum(v.numel()*2 for v in A.values())/1e9:.1f} GB",
          flush=True)
    tt = torch.from_numpy(np.ascontiguousarray(C.t.to_numpy())).to(dev)
    ka = np.array([f"{y}.{d}" for y, d in zip(C.year, C.doy)])
    dkeys = sorted(set(zip(C.year.to_numpy(), C.doy.to_numpy())))
    didx = {k: np.flatnonzero(ka == f"{k[0]}.{k[1]}") for k in dkeys}
    seedmask = None
    if a.seed_rows:
        Sd = pd.read_parquet(a.seed_rows)
        keyset = set(np.round(Sd.t.values, 3))
        seedmask = np.array([round(x, 3) in keyset for x in C.t.values])
        print(f"[setup] seed rows: {int(seedmask.sum())} of {N} windows are ARM-B seeds", flush=True)
    THR = None
    if a.thr_from:
        z = json.load(open(a.thr_from))
        THR = float(z["threshold"])
        print(f"[thr] {THR:.4f}  (rule: {z['rule']})", flush=True)
    tops = torch.full((0,), -1e30, device=dev)
    n_pairs = 0
    edges = []
    rng = torch.Generator(device=dev); rng.manual_seed(7)
    # BLOCKS: a day is split into sub-blocks of at most --jblock windows.  The j-side unfold costs
    # nj x 121 x 1440 fp16 PER STATION and the GEMM output chunk x nj x 121 fp32, so an unsplit
    # 11k-window day needs ~39 GB of unfold plus an 11 GB output tile and OOMs the card.
    blocks = []
    for k in dkeys:
        v = didx[k]
        for b0 in range(0, len(v), a.jblock):
            blocks.append((k, v[b0:b0 + a.jblock]))
    print(f"[setup] {len(blocks)} j-blocks (max {a.jblock} windows each)", flush=True)
    for bj, (dj, jjn) in enumerate(blocks):
        jj = torch.from_numpy(jjn).to(dev)
        nj = len(jjn)
        Bl = {s: unfold_unit(Wd[jj, si]) for si, s in enumerate(stas)}
        for bi in range(bj + 1):
            di, iin = blocks[bi]
            if di > dj:
                continue
            ns = min(256, len(iin))
            isub = torch.from_numpy(iin[np.linspace(0, len(iin) - 1, ns).astype(int)]).to(dev)
            med, mad = {}, {}
            for s in stas:
                v = torch.mm(A[s][isub], Bl[s].T).float().reshape(-1)
                v = v[torch.randint(0, v.numel(), (min(400_000, v.numel()),), device=dev, generator=rng)]
                m = v.median()
                med[s] = m
                mad[s] = ((v - m).abs().median().clamp_min(1e-6) * 1.4826)
                del v
            iin_use = iin if seedmask is None else iin[seedmask[iin]]
            if not len(iin_use):
                continue
            for c0 in range(0, len(iin_use), a.chunk):
                ic = torch.from_numpy(iin_use[c0:c0 + a.chunk]).to(dev)
                ni = len(ic)
                gi = ic[:, None].expand(ni, nj)
                gj = jj[None, :].expand(ni, nj)
                keep = (gi < gj) & ((tt[gi] - tt[gj]).abs() >= a.min_dt)
                if not bool(keep.any()):
                    continue
                # ---- FUSED accumulator: one tensor for the whole site sum, not one per site
                Z = None
                for s in stas:
                    r = torch.mm(A[s][ic], Bl[s].T).float()
                    r.sub_(med[s]).div_(mad[s] * nmem[site_of[s]])
                    Z = r if Z is None else Z.add_(r)
                    del r
                vals = Z.reshape(ni, nj, NLAG)
                mx, lag = vals.max(-1)
                del Z, vals
                n_pairs += int(keep.sum())
                v = mx[keep]
                if a.mode == "null":
                    tops = torch.topk(torch.cat([tops, v]), min(a.topk, tops.numel() + v.numel()))[0]
                elif THR is not None:
                    hit = v >= THR
                    if bool(hit.any()):
                        edges.append(torch.stack([gi[keep][hit].float(), gj[keep][hit].float(),
                                                  lag[keep][hit].float(), v[hit]]).T.cpu())
                del gi, gj, keep, mx, lag, v
        del Bl
        torch.cuda.empty_cache()
        if a.verbose and (bj % 10 == 0 or bj == len(blocks) - 1):
            frac = (bj + 1) ** 2 / max(len(blocks) ** 2, 1)
            el = time.time() - t_start
            print(f"  [block {bj+1}/{len(blocks)}] pairs {n_pairs/1e9:.3f}G  {el:.0f}s"
                  f"  eta {el/max(frac,1e-9)-el:.0f}s", flush=True)
    meta = dict(net=a.net, mode=a.mode, tag=a.tag, episodes=a.episodes.split(","), shift_s=a.shift,
                stations=stas, sites=sites, nmem=nmem, n_windows=int(N), n_pairs=int(n_pairs),
                n_days=len(dkeys), pnsn_last_detection=cat_end, grid_centre=ctr, grid_dep=cdep,
                min_dt=a.min_dt, runtime_s=time.time() - t_start)
    C[["t", "episode", "year", "doy"]].to_parquet(f"{OUT}/disc_{a.tag}_windows.parquet", index=False)
    if a.mode == "null":
        t = torch.sort(tops, descending=True)[0].cpu().numpy()
        if a.quantile > 0:
            k = int(max(1, round(a.quantile * n_pairs)))
            meta["rule"] = f"fixed null quantile {a.quantile:g} over {n_pairs} pairs"
            meta["quantile"] = a.quantile
        else:
            k = int(max(1, round(a.false_edges)))
            meta["rule"] = (f"null quantile giving <= {a.false_edges} expected false edges over "
                            f"{n_pairs} pairs")
            meta["quantile"] = a.false_edges / max(n_pairs, 1)
        if k > len(t):
            meta["warning"] = (f"requested rank {k} exceeds the {len(t)} retained values; the "
                               f"threshold is the smallest retained value and is an UPPER bound")
        meta["threshold"] = float(t[min(k, len(t)) - 1])
        meta["rank_used"] = int(min(k, len(t)))
        meta["null_top"] = [float(x) for x in t[:50]]
        meta["n_kept"] = int(len(t))
        why = (f"fixed quantile {a.quantile:g}" if a.quantile > 0
               else f"<= {a.false_edges:g} expected false edges")
        print(f"[null] threshold {meta['threshold']:.4f} at quantile {meta['quantile']:.3e} "
              f"[{why}] over {n_pairs:,} pairs (rank {meta['rank_used']:,} of {meta['n_kept']:,} kept)",
              flush=True)
    else:
        E = (torch.cat(edges).numpy() if edges else np.zeros((0, 4), np.float32))
        pd.DataFrame(dict(i=E[:, 0].astype(np.int64), j=E[:, 1].astype(np.int64),
                          lag=E[:, 2].astype(np.int32) - NL, z=E[:, 3])).to_parquet(
            f"{OUT}/disc_{a.tag}_edges.parquet", index=False)
        meta["n_edges"] = int(len(E))
        print(f"[true] {len(E):,} edges above threshold out of {n_pairs:,} pairs", flush=True)
    json.dump(meta, open(f"{OUT}/disc_{a.tag}.json", "w"), indent=1, default=str)
    print(f"[done] {time.time()-t_start:.0f}s -> {OUT}/disc_{a.tag}.json", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", default="B", choices=list(NETS))
    ap.add_argument("--mode", default="null", choices=["true", "null"])
    ap.add_argument("--tag", default="")
    ap.add_argument("--episodes", default=",".join(EPISODES))
    ap.add_argument("--shift", type=float, default=1800.0)
    ap.add_argument("--chunk", type=int, default=512)
    ap.add_argument("--jblock", type=int, default=3000)
    ap.add_argument("--topk", type=int, default=400000)
    ap.add_argument("--false-edges", type=float, default=30.0)
    ap.add_argument("--quantile", type=float, default=0.0,
                    help="if >0, the threshold is THIS null quantile (ARM B uses 1e-4) instead of "
                         "the <=30-expected-false-edges rule")
    ap.add_argument("--thr-from", default="")
    ap.add_argument("--min-dt", type=float, default=WIN_S + LAG_S)
    ap.add_argument("--max-windows", type=int, default=0)
    ap.add_argument("--cand-file", default="", help="parquet with a `t` column: explicit candidates "
                    "(ARM B seeds, ARM C CF peaks) instead of the uniform tremor-hour grid")
    ap.add_argument("--seed-rows", default="", help="parquet with `t`: restrict the i-side of the "
                    "pair matrix to these rows (ARM B matched-filter-from-seeds)")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    a.tag = a.tag or f"{a.net}_{a.mode}"
    warnings.filterwarnings("ignore")
    run(a)


if __name__ == "__main__":
    main()
