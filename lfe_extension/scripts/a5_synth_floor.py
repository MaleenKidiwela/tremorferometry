#!/usr/bin/env python
"""Track A-2 step 5(a) — THE SYNTHETIC FLOOR, run FIRST so the cheap test can kill the 5 km bar.

Forward-model Bostock's 130 family hypocentres through the SAME Crosson 1-D model (d3_vel1d) to the
ACTUAL permanent-network site set (PGC cluster = ONE site, LRIV cluster = ONE site -- Merlin F10),
add 0.3 s S / 0.2 s P Gaussian pick noise and a per-site static ("station-correction error")
uniform in +-0.5 s, and invert with the same grid search used for the real relocation.

Because forward and inverse share the velocity model, this measures ONLY the geometry + noise limit.
If the median horizontal error of THIS test exceeds 5 km, the pre-registered 5 km bar is untestable
on this geometry and must be reported as such (the real relocation is still run, for the record).

Variants reported: Run A (free depth) and Run B (Slab2 prior, sigma 5 km); S+P and S-only; the
noise-free control (grid-resolution check); and a matched run restricted, per family, to the sites
that actually yielded picks in step 5(c) when --picks is given.
Outputs: results/a5_synth_floor.csv (per family x realisation summary stats), results/a5_synth_floor.json,
         figures/a5_synthetic_floor.png
"""
from __future__ import annotations
import argparse, json, os, sys
from multiprocessing import Pool
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import a5_locate as A                                          # noqa: E402
import d3_vel1d as V                                           # noqa: E402

REP = A.REP
OUT = os.path.dirname(HERE)
SITES = ["PGC", "VGZ", "SNB", "NLLB", "LRIV", "B926", "B001", "B003", "B004"]

_G = {}


def _init(E, Nn, gz, TTS, TTP):
    _G.update(E=E, Nn=Nn, gz=gz, TTS=TTS, TTP=TTP)


def _job(t):
    (bkey, lat, lon, dep, zslab, tS0, tP0, sel, nreal, seed, sig_s, sig_p, sig_stat) = t
    rng = np.random.default_rng(seed)
    E, Nn, gz, TTS, TTP = _G["E"], _G["Nn"], _G["gz"], _G["TTS"], _G["TTP"]
    rows = []
    for ireal in range(nreal + 1):                 # realisation 0 = noise-free control
        if ireal == 0:
            tS, tP = tS0.copy(), tP0.copy()
        else:
            stat = rng.uniform(-sig_stat, sig_stat, len(tS0))
            tS = tS0 + rng.normal(0, sig_s, len(tS0)) + stat
            tP = tP0 + rng.normal(0, sig_p, len(tP0)) + stat / V.VPVS_DEFAULT
        for prior in (False, True):
            for mode in ("SP", "Sonly"):
                tPu = tP if mode == "SP" else np.full(len(tP), np.nan)
                r = A.locate(tS, tPu, sel, TTS, TTP, E, Nn, gz, zslab, use_prior=prior)
                rows.append(dict(bkey=bkey, ireal=ireal, prior=int(prior), mode=mode,
                                 lat=r[0], lon=r[1], z=r[2], n_sta=r[4], n_sp=r[5], at_bound=r[6],
                                 dh=A.epi_km(r[0], r[1], lat, lon), dz=r[2] - dep))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nreal", type=int, default=50)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--sig-s", type=float, default=0.3)
    ap.add_argument("--sig-p", type=float, default=0.2)
    ap.add_argument("--sig-stat", type=float, default=0.5)
    ap.add_argument("--picks", default="", help="a5_picks.csv -> matched per-family site sets")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    tag = ("_" + a.tag) if a.tag else ""

    F = pd.read_csv(f"{REP}/results/e5r_families_2012.csv")
    sll = A.site_ll()
    print(f"[synth] {len(F)} Bostock families, {len(SITES)} sites {SITES}", flush=True)

    xs, zs, TP, TS = A.build_tt_table()
    E, Nn, gz = A.grid_axes()
    TTS, TTP = A.precompute_tt(E, Nn, gz, SITES, sll, (xs, zs, TP), (xs, zs, TS))
    print(f"[synth] grid {len(E)} nodes x {len(gz)} depths, tables {TTS.nbytes/1e6:.0f} MB", flush=True)

    # --- unit check: the lean grid search reproduces d3_stage4_locate.locate_from_tables ----------
    try:
        sys.path.insert(0, f"{REP}/scripts")
        import importlib.util
        spec = importlib.util.spec_from_file_location("d3loc", f"{REP}/scripts/d3_stage4_locate.py")
        d3 = importlib.util.module_from_spec(spec)
        os.chdir("/home/jovyan/tremorferometry")
        spec.loader.exec_module(d3)
        ok = []
        rng = np.random.default_rng(0)
        for _ in range(5):
            r = F.sample(1, random_state=int(rng.integers(1e6))).iloc[0]
            d = np.array([A.epi_km(r.lat, r.lon, *sll[s]) for s in SITES])
            tS = np.array([V.tt(dd, r.dep, "S") for dd in d]) + rng.normal(0, .3, len(d))
            tP = np.array([V.tt(dd, r.dep, "P") for dd in d]) + rng.normal(0, .2, len(d))
            zsl = float(A.slab2_depth([r.lat], [r.lon])[0])
            m = A.locate(tS, tP, range(len(SITES)), TTS, TTP, E, Nn, gz, zsl)
            o = d3.locate_from_tables(tS, tP, list(range(len(SITES))), TTS, TTP, E, Nn, gz, zsl, A.LAT0, A.LON0)
            ok.append(A.epi_km(m[0], m[1], o[0], o[1]) < 1e-6 and abs(m[2] - o[2]) < 1e-9)
        print(f"[check] lean grid search == d3_stage4_locate.locate_from_tables on 5 random cases: "
              f"{sum(ok)}/5", flush=True)
    except Exception as e:
        print(f"[check] could not run the d3 cross-check: {type(e).__name__}: {e}", flush=True)
    os.chdir(OUT)

    PK = None
    if a.picks:
        PK = pd.read_csv(a.picks)
    zslabs = A.slab2_depth(F.lat.values, F.lon.values)
    jobs = []
    for i, r in enumerate(F.itertuples()):
        d = np.array([A.epi_km(r.lat, r.lon, *sll[s]) for s in SITES])
        tS = np.array([V.tt(dd, r.dep, "S") for dd in d])
        tP = np.array([V.tt(dd, r.dep, "P") for dd in d])
        if PK is not None:
            have = set(PK[(PK.bkey == r.bkey) & PK.t_S.notna()].site)
            sel = [j for j, s in enumerate(SITES) if s in have]
        else:
            sel = list(range(len(SITES)))
        if len(sel) < 3:
            continue
        jobs.append((int(r.bkey), float(r.lat), float(r.lon), float(r.dep), float(zslabs[i]),
                     tS, tP, sel, a.nreal, 1000 + i, a.sig_s, a.sig_p, a.sig_stat))
    with Pool(a.workers, initializer=_init, initargs=(E, Nn, gz, TTS, TTP)) as p:
        out = [x for rs in p.map(_job, jobs) for x in rs]
    R = pd.DataFrame(out)
    R.to_csv(f"{OUT}/results/a5_synth_floor{tag}.csv", index=False)

    # --- geometry table (the ACCEPT subset of step 5c is defined from it) -------------------------
    geo = []
    for r in F.itertuples():
        near = [s for s in SITES if A.epi_km(r.lat, r.lon, *sll[s]) <= 60.0]
        geo.append(dict(bkey=int(r.bkey), n_sites_60km=len(near),
                        gap_deg=A.azimuth_gap(r.lat, r.lon, near, sll) if near else 360.0,
                        gap_all_deg=A.azimuth_gap(r.lat, r.lon, SITES, sll)))
    G = pd.DataFrame(geo)
    G.to_csv(f"{OUT}/results/a5_geometry{tag}.csv", index=False)

    rep = {"n_families": int(F.bkey.nunique()), "sites": SITES, "nreal": a.nreal,
           "sig_s": a.sig_s, "sig_p": a.sig_p, "sig_stat": a.sig_stat,
           "picks_matched": bool(a.picks),
           "accept_subset_n": int(((G.n_sites_60km >= 4) & (G.gap_deg <= 180)).sum())}
    for prior in (0, 1):
        for mode in ("SP", "Sonly"):
            for lab, sub in (("noisefree", R[R.ireal == 0]), ("noisy", R[R.ireal > 0])):
                q = sub[(sub.prior == prior) & (sub["mode"] == mode)]
                if not len(q):
                    continue
                acc = q[q.bkey.isin(G[(G.n_sites_60km >= 4) & (G.gap_deg <= 180)].bkey)]
                rep[f"{'B_prior' if prior else 'A_free'}_{mode}_{lab}"] = dict(
                    n=int(len(q)), med_dh=float(q.dh.median()), p90_dh=float(q.dh.quantile(.9)),
                    med_absdz=float(q.dz.abs().median()), med_dz=float(q.dz.median()),
                    frac_at_bound=float(q.at_bound.mean()),
                    med_dh_accept=float(acc.dh.median()) if len(acc) else None,
                    med_absdz_accept=float(acc.dz.abs().median()) if len(acc) else None)
    json.dump(rep, open(f"{OUT}/results/a5_synth_floor{tag}.json", "w"), indent=1)
    print(json.dumps(rep, indent=1))

    # --- figure ----------------------------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    nz = R[R.ireal > 0]
    for prior, c in ((0, "tab:blue"), (1, "tab:red")):
        for mode, ls in (("SP", "-"), ("Sonly", "--")):
            q = nz[(nz.prior == prior) & (nz["mode"] == mode)]
            if not len(q):
                continue
            lab = f"{'Run B (Slab2 prior)' if prior else 'Run A (free depth)'} {mode}"
            x = np.sort(q.dh.values); ax[0].plot(x, np.linspace(0, 1, len(x)), ls, color=c, label=lab)
            x = np.sort(q.dz.abs().values); ax[1].plot(x, np.linspace(0, 1, len(x)), ls, color=c, label=lab)
    ax[0].axvline(5.0, color="k", lw=1, ls=":"); ax[0].set_xlim(0, 30)
    ax[0].set_xlabel("horizontal error (km)"); ax[0].set_ylabel("cumulative fraction")
    ax[0].set_title(f"synthetic floor, {a.nreal} realisations x {len(jobs)} families\n"
                    f"noise {a.sig_s} s S / {a.sig_p} s P + $\\pm${a.sig_stat} s site statics")
    ax[0].legend(fontsize=7); ax[0].grid(alpha=.3)
    ax[1].set_xlim(0, 30); ax[1].set_xlabel("|depth error| (km)"); ax[1].grid(alpha=.3)
    ax[1].set_title("depth error"); ax[1].legend(fontsize=7)
    q = nz[(nz.prior == 0) & (nz["mode"] == "SP")].groupby("bkey").dh.median()
    F2 = F.set_index("bkey")
    sc = ax[2].scatter(F2.loc[q.index, "lon"], F2.loc[q.index, "lat"], c=q.values, s=26,
                       cmap="viridis", vmin=0, vmax=max(10, float(np.percentile(q.values, 95))))
    sll2 = A.site_ll()
    ax[2].scatter([sll2[s][1] for s in SITES], [sll2[s][0] for s in SITES], marker="^", c="r", s=60)
    for s in SITES:
        ax[2].annotate(s, (sll2[s][1], sll2[s][0]), fontsize=6)
    plt.colorbar(sc, ax=ax[2], label="median horiz. error (km)")
    ax[2].set_xlabel("lon"); ax[2].set_ylabel("lat"); ax[2].set_title("Run A (free depth), S+P")
    plt.tight_layout()
    plt.savefig(f"{OUT}/figures/a5_synthetic_floor{tag}.png", dpi=130)
    print(f"[write] results/a5_synth_floor{tag}.csv/.json, results/a5_geometry{tag}.csv, "
          f"figures/a5_synthetic_floor{tag}.png")


if __name__ == "__main__":
    main()
