#!/usr/bin/env python
"""D3 stage 4c/d/e — pick P and S on the network stacks, locate each family, compare with
Lin/Bostock + Slab2, and do the event-level (double-difference) relative locations.

PICKING (automatic, with the figures written for the manual-style sanity check)
  * the stack time base is fixed by construction: sample NPRE = PRE_S*FS is the PREDICTED S arrival
    for that (family, station) at the family cell and its Slab2 depth;
  * S is picked on the horizontal energy sqrt(H1^2+H2^2), searched within +-3 s of predicted S;
    P is picked on Z, searched within +-2.5 s of the predicted P (= predicted S - (ttS-ttP));
  * picker = smoothed-envelope onset (back off from the local maximum to where the envelope falls to
    noise + 0.25*(peak-noise)) refined by an AIC pick in a +-0.75 s window on the band-passed trace;
  * a pick is used only if peak/noise-RMS >= --min-pick-snr.

LOCATION (origin-time free)
  cost(x, y, z) = sum_X [ (dS_obs - dS_pred)^2 / sS^2 ]        (dS = S times DEMEANED over stations)
                + sum_X [ (SP_obs - SP_pred)^2 / sSP^2 ]        (S-P at each station)
                + (z - z_slab2)^2 / sZ^2                        (weak Slab2 prior)
  grid search at --dx km over +---halfwidth km horizontally and --zmin..--zmax in depth, travel
  times from a pre-computed (offset, depth) table of the 1-D model in `d3_vel1d.py`.
  Uncertainties from --nboot bootstrap resamples of the station set.

EVENTS (double difference)
  the per-cut cc lags saved by d3_stage4_stacks.py (`*_lags.parquet`) are, per event, the S-time
  residual relative to the family stack at each station. For events with >= --min-sta stations the
  demeaned lag vector is inverted for a relative position offset with the S slowness partials at the
  family location; the within-family scatter of those offsets is reported.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import d3_vel1d as V                     # noqa: E402
from d3_stage4_stacks import STA_LL, FS, PRE_S, NPRE, NWIN, NOISE_S   # noqa: E402

LIN_FAM = "catalogs/lin_families_ets_2010_vi.csv"
SS, SSP, SZ = 0.30, 0.50, 5.0            # sigma_S (s), sigma_S-P (s), sigma_slab2 (km)


def km_e_n(lat, lon, lat0, lon0):
    return ((lon - lon0) * 111.2 * np.cos(np.radians(lat0)), (lat - lat0) * 111.2)


def smooth(x, n):
    if n < 2:
        return x
    k = np.ones(n) / n
    return np.convolve(x, k, mode="same")


def aic_pick(x, i0, i1):
    """Akaike-information-criterion onset inside [i0, i1) of the trace x."""
    seg = np.asarray(x[i0:i1], float)
    n = seg.size
    if n < 12:
        return (i0 + i1) // 2
    aic = np.full(n, np.inf)
    for k in range(5, n - 5):
        a, b = seg[:k], seg[k:]
        va, vb = a.var(), b.var()
        if va <= 0 or vb <= 0:
            continue
        aic[k] = k * np.log(va) + (n - k) * np.log(vb)
    if not np.isfinite(aic).any():
        return (i0 + i1) // 2
    return i0 + int(np.argmin(aic))


def env_onset(x, lo, hi, nnoise):
    """Envelope-onset pick inside [lo, hi); returns (index, peak/noise ratio) or (None, snr)."""
    from scipy.signal import hilbert
    x = np.asarray(x, float)
    env = smooth(np.abs(hilbert(x - x.mean())), int(round(0.25 * FS)))
    noise = float(np.median(env[:nnoise]))
    nrms = float(np.sqrt(np.mean(x[:nnoise] ** 2)))
    lo = max(lo, nnoise)
    hi = min(hi, len(x) - 2)
    if hi <= lo:
        return None, 0.0
    j = lo + int(np.argmax(env[lo:hi]))
    peak = env[j]
    snr = float(peak / nrms) if nrms > 0 else 0.0
    lvl = noise + 0.25 * (peak - noise)
    k = j
    while k > max(lo - int(2.0 * FS), nnoise) and env[k] > lvl:
        k -= 1
    w = int(round(0.75 * FS))
    return aic_pick(x, max(k - w, nnoise), min(k + w, len(x) - 1)), snr


def build_tt_table(xmax=200.0, dx=1.0, zmin=4.0, zmax=66.0, dz=0.5, vpvs=V.VPVS_DEFAULT,
                   cache="lfe_catalog_repair/results/d3_tt_table.npz"):
    """(offset, depth) travel-time table for bilinear interpolation, cached on disk.

    Each entry costs an 80-step bisection plus the refracted-branch scan, so a 0.25 km x 0.25 km
    table (2e5 nodes x 2 phases) takes many minutes. 1 km x 0.5 km is ~25x cheaper and the
    bilinear error on a function this smooth is at the millisecond level, far below the ~0.2 s
    pick uncertainty.
    """
    key = f"{xmax}_{dx}_{zmin}_{zmax}_{dz}_{vpvs}"
    cp_ = Path(cache) if cache else None
    if cp_ is not None and cp_.exists():
        d = np.load(cp_)
        if str(d.get("key", "")) == key:
            return d["xs"], d["zs"], d["TP"], d["TS"]
    xs = np.arange(0.0, xmax + dx, dx)
    zs = np.arange(zmin, zmax + dz, dz)
    TP = np.empty((zs.size, xs.size))
    TS = np.empty((zs.size, xs.size))
    for i, z in enumerate(zs):
        for j, x in enumerate(xs):
            TP[i, j] = V.tt(x, z, "P", vpvs)
            TS[i, j] = V.tt(x, z, "S", vpvs)
    if cp_ is not None:
        cp_.parent.mkdir(parents=True, exist_ok=True)
        np.savez(cp_, xs=xs, zs=zs, TP=TP, TS=TS, key=np.array(key))
    return xs, zs, TP, TS


def interp_tt(tab, x, z):
    xs, zs, T = tab
    xi = np.clip((x - xs[0]) / (xs[1] - xs[0]), 0, xs.size - 1.001)
    zi = np.clip((z - zs[0]) / (zs[1] - zs[0]), 0, zs.size - 1.001)
    x0, z0 = np.floor(xi).astype(int), np.floor(zi).astype(int)
    fx, fz = xi - x0, zi - z0
    return ((1 - fz) * ((1 - fx) * T[z0, x0] + fx * T[z0, x0 + 1])
            + fz * ((1 - fx) * T[z0 + 1, x0] + fx * T[z0 + 1, x0 + 1]))


def grid_axes(halfwidth, dx, zmin, zmax, dz):
    ge = np.arange(-halfwidth, halfwidth + dx, dx)
    gn = np.arange(-halfwidth, halfwidth + dx, dx)
    gz = np.arange(zmin, zmax + dz, dz)
    E, Nn = np.meshgrid(ge, gn, indexing="ij")
    return E.ravel(), Nn.ravel(), gz


def precompute_tt(E, Nn, gz, stas, lat0, lon0, tabP, tabS):
    """TTS/TTP with shape (nz, nnode, nsta), computed ONCE per family.

    The bootstrap re-solves the same grid ~30 times per family; interpolating inside every solve
    made a family cost ~15 s (an hour over the pilot), so the table is built once and every solve
    is then pure array reduction.
    """
    se = np.array([km_e_n(STA_LL[s][0], STA_LL[s][1], lat0, lon0)[0] for s in stas])
    sn = np.array([km_e_n(STA_LL[s][0], STA_LL[s][1], lat0, lon0)[1] for s in stas])
    d = np.sqrt((E[:, None] - se[None, :]) ** 2 + (Nn[:, None] - sn[None, :]) ** 2)   # (nnode, nsta)
    nz = len(gz)
    TTS = np.empty((nz, d.shape[0], d.shape[1]), np.float32)
    TTP = np.empty_like(TTS)
    for i, z in enumerate(gz):
        zz = np.full_like(d, z)
        TTS[i] = interp_tt(tabS, d, zz)
        TTP[i] = interp_tt(tabP, d, zz)
    return TTS, TTP


def locate_from_tables(tS, tP, sel, TTS, TTP, E, Nn, gz, z_slab, lat0, lon0, use_prior=True):
    """Grid search over the precomputed tables. `sel` = station column indices to use."""
    if len(sel) < 3:
        return (np.nan, np.nan, np.nan, np.nan, len(sel), 0)
    tS = np.asarray(tS, float)
    tP = np.asarray(tP, float)
    dS_obs = tS - tS.mean()
    okP = np.isfinite(tP)
    SP_obs = tS[okP] - tP[okP]
    S = TTS[:, :, sel]                                   # (nz, nnode, k)
    dS_pred = S - S.mean(2, keepdims=True)
    cost = (((dS_pred - dS_obs[None, None, :]) ** 2) / SS ** 2).sum(2)
    if okP.any():
        selP = np.asarray(sel)[okP]
        SPp = TTS[:, :, selP] - TTP[:, :, selP]
        cost = cost + (((SPp - SP_obs[None, None, :]) ** 2) / SSP ** 2).sum(2)
    if use_prior and np.isfinite(z_slab):
        cost = cost + ((gz - z_slab) ** 2 / SZ ** 2)[:, None]
    iz, ih = np.unravel_index(int(np.argmin(cost)), cost.shape)
    e, n, z = float(E[ih]), float(Nn[ih]), float(gz[iz])
    lat = lat0 + n / 111.2
    lon = lon0 + e / (111.2 * np.cos(np.radians(lat0)))
    return lat, lon, z, float(cost[iz, ih]), len(sel), int(okP.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stacks", default="lfe_catalog_repair/results/d3_stage4_stacks")
    ap.add_argument("--min-pick-snr", type=float, default=3.0)
    ap.add_argument("--min-sta", type=int, default=4)
    ap.add_argument("--halfwidth", type=float, default=25.0)
    ap.add_argument("--dx", type=float, default=0.5)
    ap.add_argument("--zmin", type=float, default=10.0)
    ap.add_argument("--zmax", type=float, default=60.0)
    ap.add_argument("--dz", type=float, default=0.5)
    ap.add_argument("--nboot", type=int, default=30)
    ap.add_argument("--vpvs-sens", default="1.75,1.82")
    ap.add_argument("--tt-dx", type=float, default=1.0)
    ap.add_argument("--tt-dz", type=float, default=0.5)
    ap.add_argument("--outfam", default="lfe_catalog_repair/results/d3_stage4_families.csv")
    ap.add_argument("--outev", default="lfe_catalog_repair/results/d3_stage4_events.csv")
    ap.add_argument("--outpicks", default="lfe_catalog_repair/results/d3_stage4_picks.csv")
    a = ap.parse_args()
    rng = np.random.default_rng(1)

    S = pd.read_csv(f"{a.stacks}.csv")
    Z = np.load(f"{a.stacks}.npz", allow_pickle=True)
    keep = S[S.keep == 1]
    print(f"[s4loc] {len(S)} family-station stacks, {len(keep)} kept, "
          f"{keep.family.nunique()} families", flush=True)
    if len(keep) == 0:
        print("[s4loc] no family-station stack passed the SNR-vs-null gate -> nothing to locate. "
              "This is the reportable outcome, not an error.")
        cols = ["family", "n_det", "n_stations", "n_SP", "lat", "lon", "depth", "cost", "sigma_h",
                "sigma_z", "depth_noprior", "slab_depth", "dist_to_Lin_km", "lin_family",
                "lin_depth", "dist_to_cell_km", "cell_lat", "cell_lon"]
        pd.DataFrame(columns=cols).to_csv(a.outfam, index=False)
        pd.DataFrame(columns=["family", "station", "has_horizontals", "i_S", "i_P", "snr_S",
                              "snr_P", "t_S", "t_P", "dt_pred_s", "sp_pred", "dist_km",
                              "snr_stack_h", "snr_stack_z"]).to_csv(a.outpicks, index=False)
        pd.DataFrame(columns=["family", "t_ref", "n_sta", "de_km", "dn_km", "dz_km",
                              "rms_s"]).to_csv(a.outev, index=False)
        return
    nnoise = int(round(NOISE_S * FS))
    print("[s4loc] building the travel-time table ...", flush=True)
    xs, zs, TP, TS = build_tt_table(dx=a.tt_dx, dz=a.tt_dz)
    tabP, tabS = (xs, zs, TP), (xs, zs, TS)

    lin = pd.read_csv(LIN_FAM)
    rows, prow = [], []
    fams = sorted(keep.family.unique())
    for fi, fam in enumerate(fams):
        lat0, lon0 = [float(x) for x in fam.split("__")[0].split("_")]
        z_slab = float(V.slab2_depth([lat0], [lon0])[0])
        g = keep[keep.family == fam]
        picks = []
        for _, r in g.iterrows():
            sta = r.station
            kz, k1, k2 = f"{fam}|{sta}|Z", f"{fam}|{sta}|H1", f"{fam}|{sta}|H2"
            if kz not in Z.files:
                continue
            z = np.asarray(Z[kz], float)
            # Z-only stations (OCP) still carry S on the vertical at non-vertical incidence, so
            # fall back to Z rather than dropping the station from the location entirely.
            has_h = (k1 in Z.files and k2 in Z.files)
            h = (np.sqrt(np.asarray(Z[k1], float) ** 2 + np.asarray(Z[k2], float) ** 2)
                 if has_h else z)
            xkm = np.hypot((lat0 - STA_LL[sta][0]) * 111.2,
                           (lon0 - STA_LL[sta][1]) * 111.2 * np.cos(np.radians(lat0)))
            sp = V.tt(xkm, z_slab, "S") - V.tt(xkm, z_slab, "P")
            iS = iP = None
            sS = sP = 0.0
            iS, sS = env_onset(h, NPRE - int(3 * FS), NPRE + int(3 * FS), nnoise)
            iPpred = NPRE - int(round(sp * FS))
            iP, sP = env_onset(z, iPpred - int(2.5 * FS), iPpred + int(2.5 * FS), nnoise)
            tS = ((iS - NPRE) / FS + r.dt_pred_s) if (iS is not None and sS >= a.min_pick_snr) else np.nan
            tP = ((iP - NPRE) / FS + r.dt_pred_s) if (iP is not None and sP >= a.min_pick_snr) else np.nan
            prow.append(dict(family=fam, station=sta, has_horizontals=int(has_h),
                             i_S=iS, i_P=iP, snr_S=sS, snr_P=sP,
                             t_S=tS, t_P=tP, dt_pred_s=r.dt_pred_s, sp_pred=sp,
                             dist_km=xkm, snr_stack_h=r.snr_h, snr_stack_z=r.snr_z))
            if np.isfinite(tS):
                picks.append((sta, tS, tP))
        if len(picks) < a.min_sta:
            rows.append(dict(family=fam, n_stations=len(picks), lat=np.nan, lon=np.nan,
                             depth=np.nan, cost=np.nan, slab_depth=z_slab,
                             cell_lat=lat0, cell_lon=lon0))
            continue
        stas = [p[0] for p in picks]
        tSv = np.array([p[1] for p in picks])
        tPv = np.array([p[2] for p in picks])
        E, Nn, gz = grid_axes(a.halfwidth, a.dx, a.zmin, a.zmax, a.dz)
        TTS, TTP = precompute_tt(E, Nn, gz, stas, lat0, lon0, tabP, tabS)
        allsel = list(range(len(stas)))
        lat, lon, z, c, nS, nSP = locate_from_tables(tSv, tPv, allsel, TTS, TTP, E, Nn, gz,
                                                     z_slab, lat0, lon0)
        # bootstrap over the station set
        bl, bo, bz = [], [], []
        for _ in range(a.nboot):
            idx = rng.choice(len(stas), len(stas), replace=True)
            if len(set(idx)) < 3:
                continue
            r2 = locate_from_tables(tSv[idx], tPv[idx], list(idx), TTS, TTP, E, Nn, gz,
                                    z_slab, lat0, lon0)
            if np.isfinite(r2[0]):
                bl.append(r2[0]); bo.append(r2[1]); bz.append(r2[2])
        if len(bl) > 3:
            de, dn = km_e_n(np.array(bl), np.array(bo), lat, lon)
            sig_h = float(np.sqrt(np.var(de) + np.var(dn)))
            sig_z = float(np.std(bz))
        else:
            sig_h = sig_z = np.nan
        # no-prior and Vp/Vs sensitivity
        z_np = locate_from_tables(tSv, tPv, allsel, TTS, TTP, E, Nn, gz, z_slab, lat0, lon0,
                                  use_prior=False)[2]
        d = np.hypot((lin.lat - lat) * 111.2, (lin.lon - lon) * 111.2 * np.cos(np.radians(lat)))
        j = int(np.argmin(d.values))
        rows.append(dict(family=fam, n_det=int(g.n_cuts.max()), n_stations=len(picks),
                         n_SP=nSP, lat=lat, lon=lon, depth=z, cost=c,
                         sigma_h=sig_h, sigma_z=sig_z, depth_noprior=z_np,
                         slab_depth=z_slab, dist_to_Lin_km=float(d.values[j]),
                         lin_family=lin.family_id.iloc[j], lin_depth=float(lin.depth_km.iloc[j]),
                         dist_to_cell_km=float(np.hypot((lat - lat0) * 111.2,
                                                        (lon - lon0) * 111.2 * np.cos(np.radians(lat0)))),
                         cell_lat=lat0, cell_lon=lon0))
        if (fi + 1) % 10 == 0:
            print(f"  located {fi+1}/{len(fams)}", flush=True)

    Fm = pd.DataFrame(rows)
    Path(a.outfam).parent.mkdir(parents=True, exist_ok=True)
    Fm.to_csv(a.outfam, index=False)
    pd.DataFrame(prow).to_csv(a.outpicks, index=False)
    ok = Fm.dropna(subset=["lat"])
    print(f"\n[s4loc] {len(ok)}/{len(Fm)} families located with >= {a.min_sta} S picks")
    if len(ok):
        print(f"  depth: median {ok.depth.median():.1f} km, Slab2 median {ok.slab_depth.median():.1f} km, "
              f"median (depth - Slab2) {(ok.depth - ok.slab_depth).median():+.1f} km "
              f"[no prior {(ok.depth_noprior - ok.slab_depth).median():+.1f}]")
        print(f"  sigma_h median {ok.sigma_h.median():.2f} km, sigma_z median {ok.sigma_z.median():.2f} km")
        print(f"  distance to nearest Lin family: median {ok.dist_to_Lin_km.median():.1f} km "
              f"(p25 {ok.dist_to_Lin_km.quantile(.25):.1f}, p75 {ok.dist_to_Lin_km.quantile(.75):.1f})")
        print(f"  distance to the discovery cell: median {ok.dist_to_cell_km.median():.1f} km")

    # ---------------- event-level relative locations ----------------
    lf = Path(f"{a.stacks}_lags.parquet")
    erows = []
    if lf.exists() and len(ok):
        L = pd.read_parquet(lf)
        L = L[L.family.isin(ok.family)]
        for fam, g in L.groupby("family"):
            fr = ok[ok.family == fam].iloc[0]
            stas = sorted(g.station.unique())
            if len(stas) < a.min_sta:
                continue
            # S slowness partials (s/km) at the family location, by finite difference
            P = {}
            for s in stas:
                la, lo = STA_LL[s]
                e0, n0 = km_e_n(fr.lat, fr.lon, la, lo)
                x0 = np.hypot(e0, n0)
                t0 = V.tt(x0, fr.depth, "S")
                h = 0.5
                de = (V.tt(np.hypot(e0 + h, n0), fr.depth, "S") - t0) / h
                dn = (V.tt(np.hypot(e0, n0 + h), fr.depth, "S") - t0) / h
                dz = (V.tt(x0, fr.depth + h, "S") - t0) / h
                P[s] = (de, dn, dz)
            piv = g.pivot_table(index="t_ref", columns="station", values="lag_s")
            ccp = g.pivot_table(index="t_ref", columns="station", values="cc")
            for t_ref, row in piv.iterrows():
                use = [s for s in stas if s in row and np.isfinite(row[s])
                       and s in ccp.columns and ccp.loc[t_ref, s] > 0.3]
                if len(use) < a.min_sta:
                    continue
                y = np.array([row[s] for s in use])
                y = y - y.mean()
                G = np.array([P[s] for s in use])
                G = G - G.mean(0)
                try:
                    m, *_ = np.linalg.lstsq(G, y, rcond=None)
                except np.linalg.LinAlgError:
                    continue
                erows.append(dict(family=fam, t_ref=float(t_ref), n_sta=len(use),
                                  de_km=float(m[0]), dn_km=float(m[1]), dz_km=float(m[2]),
                                  rms_s=float(np.sqrt(np.mean((y - G @ m) ** 2)))))
    E = pd.DataFrame(erows)
    E.to_csv(a.outev, index=False)
    if len(E):
        gg = E.groupby("family")[["de_km", "dn_km", "dz_km"]].std()
        print(f"\n[s4ev] {len(E)} events relocated in {E.family.nunique()} families "
              f"(>= {a.min_sta} stations)")
        print(f"  within-family scatter (median over families): "
              f"E {gg.de_km.median():.2f} km, N {gg.dn_km.median():.2f} km, Z {gg.dz_km.median():.2f} km; "
              f"lag rms {E.rms_s.median():.3f} s")
    print(f"  -> {a.outfam} / {a.outev} / {a.outpicks}")


if __name__ == "__main__":
    main()
