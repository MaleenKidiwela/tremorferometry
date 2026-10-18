#!/usr/bin/env python
"""Track A-2 step 5(c) — FAMILY RELOCATION from OUR OWN 2012 stacks, with the four controls.

Picks: on the from-scratch 2012 per-(family, station) stacks in e5r_templates_2012.npz (3-C, 20 Hz,
1-8 Hz, +-20 s around the Crosson-predicted S at the Bostock hypocentre, sample 400 = predicted S).
  * S on sqrt(H1^2 + H2^2) within +-3 s of the predicted S;
  * P on Z within +-2.5 s of the predicted P (= predicted S - (ttS - ttP));
  * envelope-onset + AIC, accepted at peak / noise-rms >= --min-pick-snr (noise = first 10 s, which is
    pre-P for every (family, station) pair here).
The two co-located CLUSTERS COLLAPSE TO ONE SITE (Merlin F10(i)): the site's arrival is the MEDIAN over
its members of (member pick - that member's earthquake station correction), so the PGC cluster
contributes ONE travel time, not four.  Corrections come from a5_station_corrections.csv (intraslab
earthquakes, never from Bostock's hypocentres).

Location: origin-time-free grid search (a5_locate.locate; identical cost function to
d3_stage4_locate.locate_from_tables) — Run A free depth, Run B with the Slab2 prior (sigma 5 km).

Controls: (i) SHUFFLED-PICK null (the pick vector of a random OTHER family, >= --nnull realisations,
distance measured to the family the picks were NOT taken from); (ii) spread along the strike/dip axes
of Bostock's own cloud; (iii) fraction of depths at the grid bounds; (iv) the synthetic floor of
a5_synth_floor.py, re-run matched to the per-family pick sets.

Outputs: results/a5_relocation.csv, results/a5_picks.csv, results/a5_controls.json,
         figures/a5_relocation_map.png
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
from a5_eq_corrections import env_onset                        # noqa: E402

EXT = os.path.dirname(HERE)
REP = A.REP
FS = 20.0
IDX_NOM = 400
NOISE_N = int(10 * FS)
SITES = ["PGC", "VGZ", "SNB", "NLLB", "LRIV", "B926", "B001", "B003", "B004"]

_G = {}


def _init(E, Nn, gz, TTS, TTP):
    _G.update(E=E, Nn=Nn, gz=gz, TTS=TTS, TTP=TTP)


def _solve(t):
    bkey, tS, tP, sel, zslab, lat_b, lon_b, dep_b = t
    E, Nn, gz, TTS, TTP = _G["E"], _G["Nn"], _G["gz"], _G["TTS"], _G["TTP"]
    out = {}
    for lab, prior in (("A", False), ("B", True)):
        r = A.locate(tS, tP, sel, TTS, TTP, E, Nn, gz, zslab, use_prior=prior)
        out[lab] = r
    a_, b_ = out["A"], out["B"]
    return dict(bkey=bkey, lat_A=a_[0], lon_A=a_[1], z_A=a_[2], bound_A=a_[6],
                lat_B=b_[0], lon_B=b_[1], z_B=b_[2], bound_B=b_[6],
                n_sites=a_[4], n_SP=a_[5],
                dh_A=A.epi_km(a_[0], a_[1], lat_b, lon_b), dz_A=a_[2] - dep_b,
                dh_B=A.epi_km(b_[0], b_[1], lat_b, lon_b), dz_B=b_[2] - dep_b)


def _null(t):
    """One shuffled-pick realisation set: picks of family j used for family i."""
    seed, pairs = t
    E, Nn, gz, TTS, TTP = _G["E"], _G["Nn"], _G["gz"], _G["TTS"], _G["TTP"]
    out = []
    for (bkey, tS, tP, sel, zslab, lat_b, lon_b, dep_b) in pairs:
        r = A.locate(tS, tP, sel, TTS, TTP, E, Nn, gz, zslab, use_prior=False)
        rB = A.locate(tS, tP, sel, TTS, TTP, E, Nn, gz, zslab, use_prior=True)
        out.append((bkey, A.epi_km(r[0], r[1], lat_b, lon_b), r[2] - dep_b,
                    A.epi_km(rB[0], rB[1], lat_b, lon_b), rB[2] - dep_b, r[0], r[1]))
    return out


def pick_stacks(min_snr, use_corr=True):
    Z = np.load(f"{REP}/results/e5r_templates_2012.npz", allow_pickle=True)
    QC = pd.read_csv(f"{REP}/results/e5r_template_qc_2012.csv")
    QC = QC[QC.usable & QC.template]
    F = pd.read_csv(f"{REP}/results/e5r_families_2012.csv").set_index("bkey")
    stab = A.station_table()
    C = pd.read_csv(f"{EXT}/results/a5_station_corrections.csv").set_index("sta") \
        if os.path.exists(f"{EXT}/results/a5_station_corrections.csv") else None
    rows = []
    for r in QC.itertuples():
        key = f"{r.sta}|{int(r.bkey)}"
        if key not in Z.files or int(r.bkey) not in F.index:
            continue
        st = np.asarray(Z[key], float)
        fr = F.loc[int(r.bkey)]
        d = A.epi_km(float(fr.lat), float(fr.lon), *stab[r.sta])
        ttP, ttS = float(V.tt(d, float(fr.dep), "P")), float(V.tt(d, float(fr.dep), "S"))
        H = np.sqrt(st[1] ** 2 + st[2] ** 2)
        iP0 = IDX_NOM - int(round((ttS - ttP) * FS))
        # noise window = the start of the stack up to 3 s BEFORE the predicted P (S-P is only 4-8 s
        # here, so a fixed 10-s noise window would straddle P and deflate every SNR)
        nnoise = max(60, iP0 - int(3.0 * FS))
        iS, snrS = env_onset(H, IDX_NOM - int(3 * FS), IDX_NOM + int(3 * FS), nnoise, FS)
        iP, snrP = env_onset(st[0], iP0 - int(2.5 * FS), iP0 + int(2.5 * FS), nnoise, FS)
        cS = float(C.loc[r.sta, "corr_S"]) if (use_corr and C is not None and r.sta in C.index) else 0.0
        cP = float(C.loc[r.sta, "corr_P"]) if (use_corr and C is not None and r.sta in C.index) else 0.0
        cS = 0.0 if not np.isfinite(cS) else cS
        cP = 0.0 if not np.isfinite(cP) else cP
        rows.append(dict(bkey=int(r.bkey), sta=r.sta, site=A.STA2SITE.get(r.sta, r.sta),
                         dist_km=d, ttP=ttP, ttS=ttS, snrS=snrS, snrP=snrP,
                         t_S=(ttS + (iS - IDX_NOM) / FS - cS) if (iS is not None and snrS >= min_snr) else np.nan,
                         t_P=(ttP + (iP - iP0) / FS - cP) if (iP is not None and snrP >= min_snr) else np.nan,
                         corr_S=cS, corr_P=cP, n_events=r.n_events))
    P = pd.DataFrame(rows)
    # site collapse: median over the members that were picked
    g = P.groupby(["bkey", "site"])
    S = g.agg(t_S=("t_S", "median"), t_P=("t_P", "median"), n_mem=("sta", "size"),
              n_mem_S=("t_S", "count"), n_mem_P=("t_P", "count"),
              dist_km=("dist_km", "median"), snrS=("snrS", "max")).reset_index()
    return P, S


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-pick-snr", type=float, default=3.0)
    ap.add_argument("--nnull", type=int, default=100)
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--no-corrections", action="store_true")
    ap.add_argument("--no-p", action="store_true", help="S-only sensitivity: the earthquake P residuals "
                    "grow with epicentral distance (Crosson model bias), so a constant per-station P "
                    "term cannot absorb them; this run drops the S-P constraint entirely")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    tag = ("_" + a.tag) if a.tag else ""

    F = pd.read_csv(f"{REP}/results/e5r_families_2012.csv")
    Fi = F.set_index("bkey")
    sll = A.site_ll()
    P, S = pick_stacks(a.min_pick_snr, use_corr=not a.no_corrections)
    P.to_csv(f"{EXT}/results/a5_picks{tag}.csv", index=False)
    S.to_csv(f"{EXT}/results/a5_picks_site{tag}.csv", index=False)
    print(f"[pick] {len(P)} (family, station) stacks, S accepted {int(P.t_S.notna().sum())}, "
          f"P accepted {int(P.t_P.notna().sum())}; site-collapsed rows {len(S)}, "
          f"sites with S {int(S.t_S.notna().sum())}", flush=True)

    xs, zs, TP, TS = A.build_tt_table()
    E, Nn, gz = A.grid_axes()
    TTS, TTP = A.precompute_tt(E, Nn, gz, SITES, sll, (xs, zs, TP), (xs, zs, TS))
    zslab = dict(zip(F.bkey, A.slab2_depth(F.lat.values, F.lon.values)))

    jobs, meta = [], {}
    for bkey in sorted(S.bkey.unique()):
        g = S[S.bkey == bkey].set_index("site")
        tS = np.array([g.t_S.get(s, np.nan) for s in SITES])
        tP = np.array([g.t_P.get(s, np.nan) for s in SITES])
        if a.no_p:
            tP = np.full(len(SITES), np.nan)
        sel = [i for i, s in enumerate(SITES) if np.isfinite(tS[i])]
        if len(sel) < 3:
            continue
        fr = Fi.loc[bkey]
        picked = [SITES[i] for i in sel]
        near = [s for s in picked if A.epi_km(fr.lat, fr.lon, *sll[s]) <= 60.0]
        meta[bkey] = dict(n_sites_60km=len(near), gap_deg=A.azimuth_gap(fr.lat, fr.lon, picked, sll),
                          gap60_deg=A.azimuth_gap(fr.lat, fr.lon, near, sll) if len(near) >= 2 else 360.0,
                          sites=",".join(picked))
        jobs.append((int(bkey), tS, tP, sel, float(zslab[bkey]),
                     float(fr.lat), float(fr.lon), float(fr.dep)))
    print(f"[loc] {len(jobs)} families with >= 3 picked sites", flush=True)

    with Pool(a.workers, initializer=_init, initargs=(E, Nn, gz, TTS, TTP)) as p:
        res = p.map(_solve, jobs)
        R = pd.DataFrame(res)
        for k in ("n_sites_60km", "gap_deg", "gap60_deg", "sites"):
            R[k] = [meta[b][k] for b in R.bkey]
        R["lat_b"] = [float(Fi.loc[b, "lat"]) for b in R.bkey]
        R["lon_b"] = [float(Fi.loc[b, "lon"]) for b in R.bkey]
        R["dep_b"] = [float(Fi.loc[b, "dep"]) for b in R.bkey]
        R["accept_geom"] = (R.n_sites_60km >= 4) & (R.gap_deg <= 180)

        # ---- shuffled-pick null -------------------------------------------------------------
        rng = np.random.default_rng(11)
        nulljobs = []
        for r_ in range(a.nnull):
            perm = rng.permutation(len(jobs))
            while np.any(perm == np.arange(len(jobs))):          # no family keeps its own picks
                perm = rng.permutation(len(jobs))
            pairs = [(jobs[i][0], jobs[perm[i]][1], jobs[perm[i]][2], jobs[perm[i]][3],
                      jobs[i][4], jobs[i][5], jobs[i][6], jobs[i][7]) for i in range(len(jobs))]
            nulljobs.append((r_, pairs))
        NU = [x for chunk in p.map(_null, nulljobs) for x in chunk]
    N = pd.DataFrame(NU, columns=["bkey", "dh_A", "dz_A", "dh_B", "dz_B", "lat_A", "lon_A"])
    N.to_csv(f"{EXT}/results/a5_null{tag}.csv", index=False)
    R.to_csv(f"{EXT}/results/a5_relocation{tag}.csv", index=False)

    # ---- controls ----------------------------------------------------------------------------
    acc = R[R.accept_geom]
    naccb = set(acc.bkey)
    Na = N[N.bkey.isin(naccb)]

    def spread(lat, lon):
        """std along the principal (strike) and secondary (dip) axes of BOSTOCK's own cloud."""
        eb, nb = A.km_e_n(Fi.loc[list(naccb), "lat"].values, Fi.loc[list(naccb), "lon"].values)
        Xb = np.column_stack([eb, nb]) - np.array([eb.mean(), nb.mean()])
        u, s, vt = np.linalg.svd(Xb, full_matrices=False)
        e, n = A.km_e_n(np.asarray(lat), np.asarray(lon))
        Xo = np.column_stack([e, n]) - np.array([eb.mean(), nb.mean()])
        pr = Xo @ vt.T
        return float(pr[:, 0].std()), float(pr[:, 1].std()), vt

    sb_strike, sb_dip, _ = spread(Fi.loc[list(naccb), "lat"].values, Fi.loc[list(naccb), "lon"].values)
    so_strike, so_dip, _ = spread(acc.lat_A.values, acc.lon_A.values)
    soB_strike, soB_dip, _ = spread(acc.lat_B.values, acc.lon_B.values)

    ctl = dict(
        n_families_located=int(len(R)), n_accept_geom=int(len(acc)),
        min_pick_snr=a.min_pick_snr, corrections_applied=bool(not a.no_corrections),
        nnull=a.nnull,
        runA=dict(med_dh=float(acc.dh_A.median()), p90_dh=float(acc.dh_A.quantile(.9)),
                  med_absdz=float(acc.dz_A.abs().median()), med_dz=float(acc.dz_A.median()),
                  frac_at_bound=float(acc.bound_A.mean())),
        runB=dict(med_dh=float(acc.dh_B.median()), p90_dh=float(acc.dh_B.quantile(.9)),
                  med_absdz=float(acc.dz_B.abs().median()), med_dz=float(acc.dz_B.median()),
                  frac_at_bound=float(acc.bound_B.mean())),
        null=dict(med_dh_A=float(Na.dh_A.median()), med_dh_B=float(Na.dh_B.median()),
                  med_absdz_A=float(Na.dz_A.abs().median())),
        ratio_obs_over_null_A=float(acc.dh_A.median() / Na.dh_A.median()),
        ratio_obs_over_null_B=float(acc.dh_B.median() / Na.dh_B.median()),
        spread=dict(bostock_strike_km=sb_strike, bostock_dip_km=sb_dip,
                    ours_A_strike_km=so_strike, ours_A_dip_km=so_dip,
                    ours_B_strike_km=soB_strike, ours_B_dip_km=soB_dip,
                    ratio_A_strike=so_strike / sb_strike, ratio_A_dip=so_dip / sb_dip,
                    ratio_B_strike=soB_strike / sb_strike, ratio_B_dip=soB_dip / sb_dip),
        verdict=dict(
            horizontal_5km_A=bool(acc.dh_A.median() <= 5.0),
            horizontal_5km_B=bool(acc.dh_B.median() <= 5.0),
            beats_null_A=bool(acc.dh_A.median() < 0.5 * Na.dh_A.median()),
            beats_null_B=bool(acc.dh_B.median() < 0.5 * Na.dh_B.median()),
            spread_A=bool(min(so_strike / sb_strike, so_dip / sb_dip) >= 0.7),
            spread_B=bool(min(soB_strike / sb_strike, soB_dip / sb_dip) >= 0.7),
            bounds_A=bool(acc.bound_A.mean() <= 0.05), bounds_B=bool(acc.bound_B.mean() <= 0.05)))
    json.dump(ctl, open(f"{EXT}/results/a5_controls{tag}.json", "w"), indent=1)
    print(json.dumps(ctl, indent=1))

    # ---- figure ------------------------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(16, 5.2))
    for lab, col, k in (("Run A (free depth)", "tab:blue", "A"), ("Run B (Slab2 prior)", "tab:red", "B")):
        j = 0 if k == "A" else 1
        ax[j].scatter(acc.lon_b, acc.lat_b, s=26, facecolor="none", edgecolor="k", label="Bostock")
        ax[j].scatter(acc[f"lon_{k}"], acc[f"lat_{k}"], s=16, color=col, label="ours")
        for r_ in acc.itertuples():
            ax[j].plot([r_.lon_b, getattr(r_, f"lon_{k}")], [r_.lat_b, getattr(r_, f"lat_{k}")],
                       color=col, lw=.6, alpha=.6)
        ax[j].scatter([sll[s][1] for s in SITES], [sll[s][0] for s in SITES], marker="^", c="g", s=70)
        for s in SITES:
            ax[j].annotate(s, (sll[s][1], sll[s][0]), fontsize=6)
        ax[j].set_title(f"{lab}: median $\\Delta$h {ctl['run' + k]['med_dh']:.1f} km "
                        f"(null {ctl['null']['med_dh_' + k]:.1f} km), n={len(acc)}")
        ax[j].set_xlabel("lon"); ax[j].set_ylabel("lat"); ax[j].legend(fontsize=7); ax[j].grid(alpha=.3)
    x = np.sort(acc.dh_A.values); ax[2].plot(x, np.linspace(0, 1, len(x)), color="tab:blue", label="Run A")
    x = np.sort(acc.dh_B.values); ax[2].plot(x, np.linspace(0, 1, len(x)), color="tab:red", label="Run B")
    x = np.sort(Na.dh_A.values); ax[2].plot(x, np.linspace(0, 1, len(x)), "k--", label="shuffled-pick null")
    ax[2].axvline(5, color="k", lw=1, ls=":")
    ax[2].set_xlim(0, 60); ax[2].set_xlabel("distance to OWN Bostock family (km)")
    ax[2].set_ylabel("cumulative fraction"); ax[2].legend(fontsize=8); ax[2].grid(alpha=.3)
    ax[2].set_title("observed vs shuffled-pick null")
    plt.tight_layout(); plt.savefig(f"{EXT}/figures/a5_relocation_map{tag}.png", dpi=130)
    print(f"[write] results/a5_relocation{tag}.csv, a5_controls{tag}.json, "
          f"figures/a5_relocation_map{tag}.png")


if __name__ == "__main__":
    main()
