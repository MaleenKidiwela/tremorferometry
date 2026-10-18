#!/usr/bin/env python
"""P2 core — explicit physics back-projection of the per-station learned CF onto the Slab2 grid.

For one day:
  1. read each station's CF (p2_cf.py output), turn it into log-odds, and MAD-normalise it over that
     station-day's valid samples (the oracle's convention, so stations are comparable);
  2. merge co-located stations into 9 SITES by the mean of their MAD units;
  3. for every grid cell c and every origin time t on a 0.25 s lattice,
         B(c, t) = sum_sites  L_site( t + tt(c, site) )                       -- no learned location;
  4. MAD-normalise B per cell per day;
  5. detections = declustered peaks of max_c Z(c, t); the argmax cell is the location.

Nulls: --shuffle moveout permutes the site -> travel-time assignment (the same cells, the wrong
geometry), which is the pre-registered moveout-shuffle null.

  python scripts/p2_backproj.py --cf results/p2/cf_p2core_2012 --year 2012 --tag p2core_2012
"""
from __future__ import annotations
import os, sys, time, json, argparse
import numpy as np
import pandas as pd
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P

OUT = f"{P.RES}/p2"
FS = P.FS
DT = 0.25                      # back-projection origin-time step, s
STEP = int(round(DT * FS))     # 5 samples


def logodds(p, lo=1e-4):
    p = np.clip(p.astype(np.float32), lo, 1 - lo)
    return np.log(p / (1 - p))


def mad_norm(v, valid):
    x = v[valid]
    if x.size < 1000:
        return None
    md = np.median(x)
    sd = np.median(np.abs(x - md)) * 1.4826
    if not np.isfinite(sd) or sd <= 0:
        return None
    return (v - md) / sd


class BackProjector:
    def __init__(self, tt_npz=f"{OUT}/grid_tt.npz", dev="cuda", cell_chunk=96,
                 site_delay=None, shuffle=False, seed=0, stations=None, norm="lowhalf"):
        z = np.load(tt_npz)
        self.lat, self.lon, self.dep = z["lat"], z["lon"], z["dep"]
        self.shape = tuple(z["shape"])
        stas = [str(s) for s in z["stas"]]
        keep = [i for i, s in enumerate(stas) if stations is None or s in stations]
        self.stas = [stas[i] for i in keep]
        TT = z["tt"][:, keep]
        self.sites = sorted({P.SITE_OF[s] for s in self.stas})
        sm = np.array([self.sites.index(P.SITE_OF[s]) for s in self.stas])
        B = np.zeros((TT.shape[0], len(self.sites)))
        for j in range(len(self.sites)):
            B[:, j] = TT[:, sm == j].mean(1)
        if site_delay is not None:              # dev-fitted per-site CF lag, s
            B = B + np.array([site_delay.get(s, 0.0) for s in self.sites])[None, :]
        if shuffle:
            rng = np.random.default_rng(seed)
            B = B[:, rng.permutation(len(self.sites))]
        self.shift = np.rint(B * FS).astype(np.int64)       # (ncell, nsite) samples
        self.dev = dev
        self.cell_chunk = cell_chunk
        self.ncell = B.shape[0]
        self.norm = norm

    def run(self, Lsite, nvalid, want_maps=False):
        """Lsite: (nsite, N) float32 MAD-normalised site series (0 where the site is missing).

        Returns (zmax (nt,), cell_arg (nt,), maps or None)."""
        dev = self.dev
        N = Lsite.shape[1]
        nt = (N - int(self.shift.max()) - 1) // STEP
        L = torch.from_numpy(Lsite).to(dev)
        base = torch.arange(nt, device=dev, dtype=torch.long) * STEP
        zmax = torch.full((nt,), -1e30, device=dev)
        carg = torch.zeros(nt, dtype=torch.long, device=dev)
        maps = torch.zeros((self.ncell, nt), dtype=torch.float16, device=dev) if want_maps else None
        SH = torch.from_numpy(self.shift).to(dev)
        for c0 in range(0, self.ncell, self.cell_chunk):
            c1 = min(self.ncell, c0 + self.cell_chunk)
            acc = torch.zeros((c1 - c0, nt), device=dev)
            for j in range(L.shape[0]):
                idx = (SH[c0:c1, j:j + 1] + base[None, :]).clamp_(0, N - 1)
                acc += L[j][idx]
            sv, _ = torch.sort(acc, dim=1)
            md = sv[:, nt // 2:nt // 2 + 1]
            if self.norm == "lowhalf":
                # NOISE-side scale.  A two-sided MAD is not a noise estimator when a large part of
                # an ETS day carries tremor: it is inflated by the signal itself, which compresses
                # the whole day towards z = 0 and makes a quiet-day threshold untransferable.
                # median(|x-md| | x < md) == md - q25, and for a Gaussian that is 0.6745 sigma.
                sd = (md - sv[:, nt // 4:nt // 4 + 1]) * 1.4826
            else:
                sd = (acc - md).abs().median(dim=1, keepdim=True).values * 1.4826
            z = (acc - md) / sd.clamp_min(1e-6)
            if want_maps:
                maps[c0:c1] = z.half()
            m, i = z.max(0)
            upd = m > zmax
            carg = torch.where(upd, i + c0, carg)
            zmax = torch.where(upd, m, zmax)
        return (zmax.cpu().numpy(), carg.cpu().numpy(),
                None if maps is None else maps.cpu().numpy())


def load_site_series(cfdir, sta_list, year, doy, store, T0):
    """(nsite, NDAY) MAD-normalised log-odds, and the per-site availability."""
    sites = sorted({P.SITE_OF[s] for s in sta_list})
    acc = np.zeros((len(sites), P.NDAY), np.float32)
    cnt = np.zeros(len(sites), np.int32)
    valid_any = np.zeros(P.NDAY, bool)
    used = []
    for s in sta_list:
        f = f"{cfdir}/{s}_{year}_{doy}.npy"
        if not os.path.exists(f) or not store.has(s, year, doy):
            continue
        i = store.idx(s, year, doy)
        cf = np.load(f).astype(np.float32)
        off = int(round((store.t0[i] - T0) * FS))
        v = np.zeros(P.NDAY, np.float32)
        lo = max(0, off)
        hi = min(P.NDAY, off + len(cf))
        if hi <= lo:
            continue
        v[lo:hi] = cf[lo - off:hi - off]
        valid = np.zeros(P.NDAY, bool)
        valid[lo:hi] = True
        valid &= (v > 0)
        z = mad_norm(logodds(v), valid)
        if z is None:
            continue
        z = np.where(valid, z, 0.0).astype(np.float32)
        j = sites.index(P.SITE_OF[s])
        acc[j] += z
        cnt[j] += 1
        valid_any |= valid
        used.append(s)
    ok = cnt > 0
    acc[ok] /= cnt[ok][:, None]
    return acc[ok], [s for s, o in zip(sites, ok) if o], used, valid_any


def decluster_t(v, sep_s=2.0):
    n = max(1, int(round(sep_s / DT)))
    t = torch.from_numpy(v)[None, None]
    mx = torch.nn.functional.max_pool1d(t, kernel_size=2 * n + 1, stride=1, padding=n)[0, 0].numpy()
    return np.flatnonzero(v >= mx - 1e-9)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cf", required=True, help="directory of per-station-day CF .npy")
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--stations", default="", help="comma list; default = all S15 with a CF")
    ap.add_argument("--shuffle", action="store_true", help="moveout-shuffle null")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--site-delay", default="", help="json {site: seconds} dev-fitted CF lag")
    ap.add_argument("--keep-top", type=int, default=20000, help="peaks stored per day")
    ap.add_argument("--maps-days", default="", help="comma list of doys to also store full Z maps")
    ap.add_argument("--norm", default="lowhalf", choices=["lowhalf", "mad"])
    a = ap.parse_args()
    assert a.year != 2013, "2013 is sealed"
    os.makedirs(OUT, exist_ok=True)
    stations = [s for s in a.stations.split(",") if s] or None
    sd = json.load(open(a.site_delay)) if a.site_delay else None
    bp = BackProjector(shuffle=a.shuffle, seed=a.seed, site_delay=sd, stations=stations,
                       norm=a.norm)
    store = P.DayStore()
    di = store.di[store.di.year == a.year]
    doys = sorted(di.doy.unique())
    sta_list = stations or P.S15
    mapdays = {int(x) for x in a.maps_days.split(",") if x}
    print(f"[bp] {a.year}: {len(doys)} days, {bp.ncell} cells, sites={bp.sites}, "
          f"shuffle={a.shuffle}", flush=True)

    rows, t0 = [], time.time()
    for k, doy in enumerate(doys):
        T0 = (pd.Timestamp(f"{a.year}-01-01") + pd.Timedelta(days=int(doy) - 1)).timestamp()
        L, sites_ok, used, valid = load_site_series(a.cf, sta_list, a.year, doy, store, T0)
        if len(sites_ok) < 3:
            print(f"  doy {doy}: only {len(sites_ok)} sites, skipped", flush=True)
            continue
        Lfull = np.zeros((len(bp.sites), P.NDAY), np.float32)
        for j, s in enumerate(sites_ok):
            Lfull[bp.sites.index(s)] = L[j]
        zmax, carg, maps = bp.run(Lfull, valid.sum(), want_maps=doy in mapdays)
        pk = decluster_t(zmax)
        if len(pk) > a.keep_top:
            pk = pk[np.argsort(zmax[pk])[::-1][:a.keep_top]]
            pk.sort()
        rows.append(pd.DataFrame(dict(
            year=a.year, doy=doy, t_origin=T0 + pk * DT, z=zmax[pk], cell=carg[pk],
            lat=bp.lat[carg[pk]], lon=bp.lon[carg[pk]], dep=bp.dep[carg[pk]],
            n_sites=len(sites_ok))))
        if maps is not None:
            np.savez_compressed(f"{OUT}/maps_{a.tag}_{doy}.npz", z=maps.astype(np.float16),
                                shape=np.array(bp.shape), t0=T0, dt=DT,
                                lat=bp.lat, lon=bp.lon)
        if (k + 1) % 5 == 0:
            print(f"  {k+1}/{len(doys)} doy {doy} sites={len(sites_ok)} peaks={len(pk)} "
                  f"zmax={zmax.max():.1f}  {time.time()-t0:.0f}s", flush=True)
    D = pd.concat(rows, ignore_index=True)
    D.to_parquet(f"{OUT}/det_{a.tag}.parquet", index=False)
    json.dump(dict(tag=a.tag, year=a.year, shuffle=bool(a.shuffle), norm=a.norm, stations=sta_list,
                   sites=bp.sites, ncell=int(bp.ncell), dt=DT, n_days=len(doys),
                   n_peaks=int(len(D)), cf=a.cf, site_delay=sd),
              open(f"{OUT}/det_{a.tag}_meta.json", "w"), indent=1)
    print(f"[done] {len(D)} peaks over {len(doys)} days in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
