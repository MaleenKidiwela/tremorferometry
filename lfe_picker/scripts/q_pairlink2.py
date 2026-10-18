#!/usr/bin/env python
"""P3 STEP 1a, POOLED REPEAT — 2010 + 2012 Bostock ETS-day candidates, and a site-subset ablation.

Why.  The 2012-only kill test (scripts/q_pairlink.py, results/P3_RESULTS.md § 2) found 0 of the 10
isolated hidden families above the null's 1e-6 quantile.  The one lever that raises the number of
LINKABLE pairs quadratically rather than linearly is pooling ETS episodes: a family with n events in
one ETS has n(n-1)/2 same-family pairs, so 2010 + 2012 together give ~4x the pairs of 2012 alone at
the same per-pair detection rate.  This script runs exactly the Step-1a recipe over the pooled
candidate set (dev years only; 2011 and 2013 are not touched), with the same null construction, the
same overlap exclusion, and the 1e-6 threshold RE-DERIVED on the pooled null.

Second question, in the same pass: do the FAR sites dilute the pair statistic?  Every per-site
lag-resolved z is kept, so an arbitrary number of SITE SUBSETS can be scored from one set of
matmuls.  Three families of configuration are scored:

  all7      the 7 discovery sites (S15 minus SNB and B926)                    [the protocol]
  n4_*      the 4 sites nearest to each family, chosen from those 7           [oracle site choice]
  o4_*      PGC-cluster + B926 + VGZ + the family's nearest far borehole      [oracle, and B926 is a
                                                                               WITHHELD station]

n4_* and o4_* use the family's own position to pick the sites, which a template-free search does not
know, and o4_* additionally puts a withheld station back in.  Both are DIAGNOSTICS of dilution, not
candidate detectors, and are labelled as such everywhere.

Data note: the picker's memory-mapped store has no B001-B007 for 2010, so this script reads the
20 Hz caches directly (lfe_picker/cache20hz, then lfe_catalog_repair, then lfe_extension), in both
years, for every station.  The caches are the frozen 1-8 Hz / 20 Hz chain and are NOT rotated to ZNE;
that is immaterial here because a common orthogonal rotation of the three components leaves the 3-C
inner product of two windows at the SAME station unchanged.  The 2010 borehole caches were built by
scripts/q_cache2010.py from the read-only ../lfe_extension/traces with o_cache_month.cache_one.

    python scripts/q_pairlink2.py --mode null --tag pool_null
    python scripts/q_pairlink2.py --mode true --tag pool_true --thr-from results/p3/pairlink2_pool_null.json
"""
from __future__ import annotations
import os, sys, glob, json, time, argparse
import numpy as np
import pandas as pd
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P

P3 = f"{P.RES}/p3"
FS = 20.0
WIN_S, LAG_S = 24.0, 3.0
NW = int(round(WIN_S * FS))
NL = int(round(LAG_S * FS))
NB = NW + 2 * NL
NLAG = 2 * NL + 1

DISC = [s for s in P.S15 if s not in ("SNB", "B926")]     # the 13-station discovery network
ORACLE_EXTRA = ["B926"]                                    # only ever used by the o4_* diagnostics
STAS = DISC + ORACLE_EXTRA
SITE_OF = {s: P.SITE_OF[s] for s in STAS}
SITES = sorted(set(SITE_OF.values()))
NMEM = {q: sum(1 for s in STAS if SITE_OF[s] == q) for q in SITES}
DISC_SITES = sorted({P.SITE_OF[s] for s in DISC})
GRID_CENTRE = (0.5 * (P.GRID_LAT[0] + P.GRID_LAT[1]), 0.5 * (P.GRID_LON[0] + P.GRID_LON[1]))
CACHES = [f"{P.PICK}/cache20hz", f"{P.REP}/cache20hz", f"{P.EXT}/cache20hz"]
THRN = ["q_0.001", "q_0.0001", "q_1e-5", "q_1e-6"]


# ----------------------------------------------------------------------------- inputs
class CacheStore:
    """(sta, year, doy) -> (t0, x (3, N) float32) from the 20 Hz caches, first root wins."""

    def __init__(self, stas, years):
        self.idx = {}
        for s in stas:
            ns = P.NETSTA[s]
            for y in years:
                for root in CACHES:
                    for f in glob.glob(f"{root}/{ns}/{y}/*.npz"):
                        self.idx.setdefault((s, y, int(os.path.basename(f)[:3])), f)

    def has(self, s, y, d):
        return (s, y, d) in self.idx

    def get(self, s, y, d):
        z = np.load(self.idx[(s, y, d)])
        return float(z["t0"]), z["x"]


def candidates(years, stas, store, min_events=100):
    """Bostock ETS-day events of `years` on days cached at EVERY station of `stas`."""
    B = P.bostock()
    keep = []
    for y in years:
        b = B[B.year == y]
        cnt = b.groupby("doy").size()
        ets = set(cnt[cnt >= min_events].index)
        have = {d for d in ets if all(store.has(s, y, int(d)) for s in stas)}
        keep.append(b[b.doy.isin(have)])
        print(f"[cand] {y}: {len(ets)} ETS days, {len(have)} with all {len(stas)} stations cached, "
              f"{int(b.doy.isin(have).sum())} events")
    C = pd.concat(keep, ignore_index=True)
    return C.sort_values(["year", "doy", "t_origin"]).reset_index(drop=True)


def station_constants(stas):
    import p_vel1d as V
    ST, CR = P.station_table(), P.station_corr()
    dep = float(np.atleast_1d(P.slab2_interp()(*GRID_CENTRE))[0])
    out = {}
    for s in stas:
        d = float(P.epi_km(GRID_CENTRE[0], GRID_CENTRE[1], ST.lat[s], ST.lon[s]))
        out[s] = float(V.tt(d, dep, "S")) + CR[s]
    return out, dep


def cut_windows(B, cst, mode, store, stas):
    N, S = len(B), len(stas)
    W = np.zeros((N, S, 3, NB), np.float16)
    ok = np.zeros((N, S), bool)
    for si, s in enumerate(stas):
        shift = int(round(si * 1800.0 * FS)) if mode == "null" else 0
        for (y, doy), g in B.groupby(["year", "doy"]):
            if not store.has(s, int(y), int(doy)):
                continue
            t0, X = store.get(s, int(y), int(doy))
            n = X.shape[1]
            if n < 4 * NB:
                continue
            i0 = np.rint((g.t_origin.to_numpy() + cst[s] - LAG_S - t0) * FS).astype(np.int64)
            for row, aa in zip(g.index.to_numpy(), i0):
                if aa < 0 or aa + NB > n:
                    continue
                true_seg = X[:, aa:aa + NB]
                if not np.square(true_seg, dtype=np.float64).sum() > 0:
                    continue
                seg = X[:, (np.arange(aa, aa + NB) + shift) % n] if shift else true_seg
                e = float(np.square(seg, dtype=np.float64).sum())
                if not (np.isfinite(e) and e > 0):
                    continue
                # the 20 Hz caches hold RAW counts (up to ~1e6 at the boreholes), which OVERFLOW
                # float16 and poison the whole statistic with inf/NaN.  Every downstream
                # normalisation is per row, so a per-window rms scaling here is free.
                W[row, si] = (seg / np.sqrt(e / seg.size)).astype(np.float16)
                ok[row, si] = bool(np.isfinite(W[row, si]).all())
        print(f"  [cut {mode}] {s}: {ok[:, si].sum()}/{N}", flush=True)
    return W, ok


def site_configs():
    """{name: (list of sites, set of bkeys it is the chosen subset for, kind)}"""
    ST, F = P.station_table(), P.families()
    loc = {q: (np.mean([ST.lat[s] for s in STAS if SITE_OF[s] == q]),
               np.mean([ST.lon[s] for s in STAS if SITE_OF[s] == q])) for q in SITES}
    cfg = {"all7": (DISC_SITES, set(F.index.astype(int)), "protocol")}
    far = [q for q in DISC_SITES if q not in ("PGCq", "VGZ")]
    for b in F.index:
        d = {q: float(P.epi_km(F.lat[b], F.lon[b], *loc[q])) for q in SITES}
        n4 = tuple(sorted(sorted(DISC_SITES, key=lambda q: d[q])[:4]))
        o4 = tuple(sorted(["PGCq", "B926", "VGZ", min(far, key=lambda q: d[q])]))
        for pre, k, kind in (("n4", n4, "oracle-sites"), ("o4", o4, "oracle-sites+withheld")):
            name = f"{pre}_" + "+".join(k)
            cfg.setdefault(name, (list(k), set(), kind))[1].add(int(b))
    return cfg


def unfold_unit(Wd):
    n = Wd.shape[0]
    U = Wd.unfold(2, NW, 1).permute(0, 2, 1, 3).reshape(n * NLAG, 3 * NW).float()
    U = U / U.pow(2).sum(1, keepdim=True).sqrt().clamp_min(1e-12)
    return U.half()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["true", "null"], required=True)
    ap.add_argument("--tag", default="")
    ap.add_argument("--years", default="2010,2012")
    ap.add_argument("--chunk", type=int, default=512)
    ap.add_argument("--thr-from", default="")
    ap.add_argument("--topk", type=int, default=50000)
    ap.add_argument("--min-dt", type=float, default=WIN_S + LAG_S)
    a = ap.parse_args()
    tag = a.tag or f"pool_{a.mode}"
    years = [int(y) for y in a.years.split(",")]
    os.makedirs(P3, exist_ok=True)
    dev = "cuda"
    t_start = time.time()

    store = CacheStore(STAS, years)
    B = candidates(years, STAS, store)
    cst, cdep = station_constants(STAS)
    CFG = site_configs()
    cnames = list(CFG)
    print(f"[setup] {len(B)} pooled candidates, stations {STAS}, sites {SITES}, "
          f"{len(cnames)} site configurations")

    W, ok = cut_windows(B, cst, a.mode, store, STAS)
    use = ok.all(1)
    assert np.isfinite(W[use]).all(), "non-finite windows survived the cut"
    print(f"[setup] usable events {use.sum()}/{len(B)} (all {len(STAS)} stations valid)")
    B = B[use].reset_index(drop=True)
    W = W[use]
    N = len(B)
    print(B.groupby("year").size().to_string())

    Wd = torch.from_numpy(W).to(dev)
    del W
    A = {}
    for si, s in enumerate(STAS):
        x = Wd[:, si, :, NL:NL + NW].reshape(N, 3 * NW).float()
        A[s] = (x / x.pow(2).sum(1, keepdim=True).sqrt().clamp_min(1e-12)).half()
        del x
    print(f"[gpu] windows {Wd.numel()*2/1e9:.2f} GB", flush=True)

    key = list(zip(B.year.to_numpy(), B.doy.to_numpy()))
    days = sorted(set(key))
    ka = np.array([f"{y}.{d}" for y, d in key])
    didx = {k: np.flatnonzero(ka == f"{k[0]}.{k[1]}") for k in days}
    bkey = torch.from_numpy(B.bkey.to_numpy().astype(np.int64)).to(dev)
    torig = torch.from_numpy(B.t_origin.to_numpy()).to(dev)

    F = P.families()
    fl = F[["lat", "lon"]].to_numpy()
    D5 = np.zeros((130, 130), bool)
    for i in range(len(fl)):
        D5[i] = P.epi_km(fl[:, 0], fl[:, 1], fl[i, 0], fl[i, 1]) < 5.0
    D5d = torch.from_numpy(D5).to(dev)

    THRV = None
    if a.thr_from:
        z = json.load(open(a.thr_from))
        THRV = {c: torch.tensor([float(z["thresholds"][c][k]) for k in THRN], device=dev)
                for c in cnames}
        print("[thr] all7 null quantiles " +
              "  ".join(f"{k}={float(z['thresholds']['all7'][k]):.2f}" for k in THRN))

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
    hist = {k: torch.zeros(HIST_N, dtype=torch.float64, device=dev)
            for k in ("all", "same", "near", "far")}
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
                gi = ic[:, None].expand(ni, nj)
                gj = jj[None, :].expand(ni, nj)
                keep = (gi < gj) & ((torig[gi] - torig[gj]).abs() >= a.min_dt)
                if not bool(keep.any()):
                    continue
                Z = {q: torch.zeros(ni, nj * NLAG, device=dev) for q in SITES}
                for s in STAS:
                    r = torch.mm(A[s][ic], Bl[s].T).float()
                    Z[SITE_OF[s]] += (r - med[s]) / (mad[s] * NMEM[SITE_OF[s]])
                    del r
                gik, gjk = gi[keep], gj[keep]
                bi, bj = bkey[gik], bkey[gjk]
                sf = bi == bj
                near = D5d[bi, bj] & ~sf
                ff = (~near) & (~sf)
                if a.mode == "true":
                    same_tot += torch.bincount(bi[sf], minlength=130)
                    ali_tot += near.sum(); far_tot += ff.sum()
                n_pairs += int(keep.sum())
                for ci, c in enumerate(cnames):
                    sl = CFG[c][0]
                    acc = Z[sl[0]].clone()
                    for q in sl[1:]:
                        acc += Z[q]
                    vals = acc.reshape(ni, nj, NLAG).amax(-1)[keep]
                    del acc
                    if c == "all7":
                        hist["all"] += torch.histc(vals, HIST_N, HIST_LO, HIST_HI).double()
                        if a.mode == "true":
                            hist["same"] += torch.histc(vals[sf], HIST_N, HIST_LO, HIST_HI).double()
                            hist["near"] += torch.histc(vals[near], HIST_N, HIST_LO, HIST_HI).double()
                            hist["far"] += torch.histc(vals[ff], HIST_N, HIST_LO, HIST_HI).double()
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
        print(f"[day {dj}] pairs {n_pairs/1e6:.1f}M  ({time.time()-t_start:.0f}s)", flush=True)

    out = dict(mode=a.mode, tag=tag, years=years, n_events=int(N), n_pairs=int(n_pairs),
               stations=STAS, sites=SITES, nmem=NMEM, min_dt=a.min_dt,
               n_by_year={str(k): int(v) for k, v in B.groupby("year").size().items()},
               configs={c: dict(sites=CFG[c][0], kind=CFG[c][2], n_families=len(CFG[c][1]))
                        for c in cnames},
               c_sta=cst, grid_centre=GRID_CENTRE, grid_dep=cdep,
               runtime_s=time.time() - t_start)
    np.savez(f"{P3}/pairlink2_{tag}_hist.npz", lo=HIST_LO, hi=HIST_HI,
             **{k: v.cpu().numpy() for k, v in hist.items()})
    B[["bkey", "year", "doy", "t_origin", "hidden", "mag"]].to_parquet(
        f"{P3}/pairlink2_{tag}_events.parquet", index=False)
    if a.mode == "null":
        thr = {}
        for c in cnames:
            t = torch.sort(tops[c], descending=True)[0].cpu().numpy()
            thr[c] = {}
            for q, name in [(1e-3, "q_0.001"), (1e-4, "q_0.0001"), (1e-5, "q_1e-5"), (1e-6, "q_1e-6")]:
                r = int(max(1, round(q * n_pairs)))
                thr[c][name] = float(t[min(r, len(t)) - 1])
            thr[c]["max"] = float(t[0])
        out["thresholds"] = thr
        print(json.dumps(thr, indent=1))
    else:
        rows = []
        st = same_tot.cpu().numpy()
        SH = same_hit.cpu().numpy()
        for ci, c in enumerate(cnames):
            fams = CFG[c][1]
            for b in range(130):
                if st[b] == 0 or b not in fams:
                    continue
                rows.append(dict(config=c, kind=CFG[c][2], bkey=b, n_same_pairs=int(st[b]),
                                 **{k: int(SH[ci, ti, b]) for ti, k in enumerate(THRN)}))
        pd.DataFrame(rows).to_csv(f"{P3}/pairlink2_{tag}_by_family.csv", index=False)
        out["alias_tot"] = int(ali_tot); out["far_tot"] = int(far_tot)
        out["same_tot"] = int(st.sum())
        out["alias_hit"] = {c: {k: int(v) for k, v in zip(THRN, ali_hit[ci].cpu().numpy())}
                            for ci, c in enumerate(cnames)}
        out["far_hit"] = {c: {k: int(v) for k, v in zip(THRN, far_hit[ci].cpu().numpy())}
                          for ci, c in enumerate(cnames)}
    json.dump(out, open(f"{P3}/pairlink2_{tag}.json", "w"), indent=1, default=str)
    print(f"[done] {time.time()-t_start:.0f}s -> {P3}/pairlink2_{tag}.json")


if __name__ == "__main__":
    main()
