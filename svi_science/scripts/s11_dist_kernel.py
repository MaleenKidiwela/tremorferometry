#!/usr/bin/env python
"""s11_dist_kernel -- S11 distant-station lever (distant/PLAN.md §5): S-pinned time-residence kernels for the far pairs with
the UNCHANGED s11_kernel_mc physics (simulate / coef / window_weights imported; the source region's CVM via S11_REGION),
sensor depths (s11 tables / FDSN), lapse weights from the far store grand stacks on their own tau grid. Unit test
|sum a - 1| <= 0.02; stop if a_d < 0 for > 10 % of far pairs (ell40, 2-5 s, 3C).
Writes distant/<region>/tables/{kernel_mc.npz, kernel_coef.parquet, kernel_unit_test.csv, kernel_summary.csv};
marker KERNEL_DONE or KERNEL_STOP.
    S11_REGION=SVI python s11_dist_kernel.py --workers 10 --npart 100000
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time, glob
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_dist_common as D
import s11_kernel_mc as K

assert os.path.basename(K.C.VEL) == os.path.basename(D.VEL), (K.C.VEL, D.VEL)


def lapse_weights(tag=""):
    from scipy.signal import hilbert
    W = {}
    for p in sorted(glob.glob(f"{D.DTB}/store{tag}/*.npz")):
        z = np.load(p); tau = z["tau"]
        for key, G, tS in zip(z["keys"], z["G1"], z["tS"]):
            e = np.abs(hilbert(np.nan_to_num(G), axis=1)) ** 2
            tr = tau - tS
            Eb = np.array([[e[c][(tr >= K.BIN_EDGES[b]) & (tr < K.BIN_EDGES[b + 1])].mean() for b in range(K.NB)] for c in range(3)])
            Eb[~np.isfinite(G).all(1)] = np.nan
            W[str(key)] = Eb
    return W


def _pair(args):
    return K._pair(args)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--npart", type=int, default=100000); ap.add_argument("--pairs", type=int, default=0)
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    D.ensure(); log = C.Log(f"{D.DLOG}/kernel{a.tag}.log"); t0 = time.time()
    NT = D.near_tables(); F = NT["F"]
    ser = pd.read_csv(f"{D.DTB}/series{a.tag}.csv"); ser = ser[(ser.n_data > 0) & (ser.causal25 > 1.5)]
    pairs = sorted(set(zip(ser.fam.astype(int), ser.netsta)))
    XY = D.station_coords(); dep = D.sensor_depths(sorted({p[1] for p in pairs}), log)
    lat0, lon0 = float(F.lat.mean()), float(F.lon.mean())
    pj = lambda la, lo: np.array([(lo - lon0) * 111.19 * np.cos(np.radians(lat0)), (la - lat0) * 111.19])
    Fx = {int(r.fam): np.r_[pj(r.lat, r.lon), r.depth_km] for r in F.itertuples()}
    Sx = {ns: np.r_[pj(XY.lat[ns], XY.lon[ns]), dep[ns] / 1000.0] for ns in {p[1] for p in pairs}}
    log(f"region {D.REG}, velocity {os.path.basename(K.C.VEL)}; far pairs with a causal series {len(pairs)}; "
        f"{a.npart} particles x {len(K.VARIANTS)} variants, {a.workers} workers")
    jobs = [(2000 + i, f, ns, Fx[f], Sx[ns], a.npart) for i, (f, ns) in enumerate(pairs)]   # seeds distinct from near runs
    if a.pairs:
        jobs = jobs[:a.pairs]
    import multiprocessing as mp
    res = {}
    with mp.get_context("fork").Pool(a.workers) as p:
        for k, (i, fam, ns, out) in enumerate(p.imap_unordered(_pair, jobs)):
            res[(fam, ns)] = out
            if (k + 1) % 20 == 0 or k + 1 == len(jobs):
                log(f"  {k + 1}/{len(jobs)} pairs {time.time() - t0:.0f}s")
    keys = sorted(res); V = list(K.VARIANTS)
    np.savez_compressed(f"{D.DTB}/kernel_mc{a.tag}.npz", fam=np.array([k[0] for k in keys]), netsta=np.array([k[1] for k in keys]),
                        variants=np.array(V), cats=np.array(C.CATS), bin_edges=K.BIN_EDGES,
                        **{f: np.array([[res[k][v][f] for v in V] for k in keys]) for f in ("Rb", "taub", "sumc", "neff", "rS", "tS", "legs", "lenchk")})
    W = lapse_weights(a.tag)
    rows = []
    for key, Eb in W.items():
        fam, ns, ep = key.split("|"); fam = int(fam)
        if (fam, ns) not in res:
            continue
        for v in V:
            k = res[(fam, ns)][v]
            for wn, win in C.WINS.items():
                for comp in ("3c", "Z", "H"):
                    if comp != "3c" and wn != C.HEAD:
                        continue
                    wb = K.window_weights(Eb, win, comp)
                    if wb.sum() <= 0:
                        continue
                    av = K.coef(k, wb)
                    r = dict(series=key, fam=fam, netsta=ns, epoch=ep, variant=v, win=wn, comp=comp, tS_mc=k["tS"],
                             unit_sum=float(np.nansum(av)), a_s=float(C.group_coef(av, "s")), a_m=float(C.group_coef(av, "m")),
                             a_d=float(C.group_coef(av, "d")), a_b=float(C.group_coef(av, "b")),
                             neff_min_win=float(np.nanmin(np.where(wb > 0, k["neff"], np.inf))),
                             hypo_km=float(np.linalg.norm(Sx[ns] - Fx[fam])), epi_km=float(np.linalg.norm(Sx[ns][:2] - Fx[fam][:2])),
                             fam_depth=float(Fx[fam][2]), sensor_depth_km=float(Sx[ns][2]))
                    a3 = av.reshape(len(C.LAYER_NAMES), 3)
                    r["a_s_sta_share"] = float(a3[C.GROUPS["s"], 1].sum() / r["a_s"]) if r["a_s"] != 0 else np.nan
                    r["a_d_src_share"] = float(a3[C.GROUPS["d"], 0].sum() / r["a_d"]) if r["a_d"] != 0 else np.nan
                    for ci, cn in enumerate(C.CATS):
                        r[f"c_{cn}"] = float(av[ci])
                    rows.append(r)
    KC = pd.DataFrame(rows); KC.to_parquet(f"{D.DTB}/kernel_coef{a.tag}.parquet", index=False)
    UT = KC.assign(dev=(KC.unit_sum - 1).abs())
    ut = UT.groupby(["variant", "win", "comp"]).agg(n=("dev", "size"), max_dev=("dev", "max"),
                                                    fail=("dev", lambda x: int((x > 0.02).sum()))).reset_index()
    ut.to_csv(f"{D.DTB}/kernel_unit_test{a.tag}.csv", index=False)
    prim = KC[(KC.variant == "ell40") & (KC.win == C.HEAD) & (KC.comp == "3c")]
    pp = prim.groupby(["fam", "netsta"]).first().reset_index()
    neg = float((pp.a_d < 0).mean())
    sm = []
    for (v, wn, comp), g in KC.groupby(["variant", "win", "comp"]):
        sm.append(dict(variant=v, win=wn, comp=comp, n=len(g), a_s_med=g.a_s.median(), a_m_med=g.a_m.median(), a_d_med=g.a_d.median(),
                       a_d_q10=g.a_d.quantile(.1), a_d_q90=g.a_d.quantile(.9), a_b_med=g.a_b.median(), a_d_neg_frac=float((g.a_d < 0).mean()),
                       neff_min_med=g.neff_min_win.median(), unit_maxdev=float((g.unit_sum - 1).abs().max())))
    pd.DataFrame(sm).to_csv(f"{D.DTB}/kernel_summary{a.tag}.csv", index=False)
    log(f"unit test: max |sum-1| {UT.dev.max():.2e}, failures {int((UT.dev > 0.02).sum())}")
    log(f"far primary (ell40, 2-5 s, 3C): a_s {pp.a_s.median():.3f}, a_m {pp.a_m.median():.3f}, a_d {pp.a_d.median():.3f} "
        f"(10-90 % {pp.a_d.quantile(.1):.3f}..{pp.a_d.quantile(.9):.3f}), a_b {pp.a_b.median():.3f} [forecast 0.23 / -0.34 / 0.91 / 0.20]; "
        f"a_d < 0 for {100 * neg:.1f} % of pairs; neff (min over window bins) median {pp.neff_min_win.median():.0f}")
    stop = bool((UT.dev > 0.02).any() or neg > 0.10)
    if not a.pairs:
        mk = f"{D.DOUT}/KERNEL{a.tag}_{'STOP' if stop else 'DONE'}"
        open(mk, "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n") + f"unit_maxdev {UT.dev.max():.3e} a_d_neg {neg:.3f}\n")
        log(f"wrote {mk} ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
