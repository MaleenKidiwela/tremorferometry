#!/usr/bin/env python
"""s1_geodetic (2026-10-07) -- SVI science study S1: episode size from GNSS (Tier A, model-free) and the
slip-model gate (Tier B).  Frozen rules: svi_science/results/s1/PLAN_S1.md (written before any number).

    python svi_science/scripts/s1_geodetic.py

Read-only use of /home/jovyan/cascadia-sse-strain (package import only) and of the slip_compare outputs of
lfe_extension/scripts/slip_build.py (corrected dates, thrust-positive).  CPU only, < 5 GB.
Outputs: svi_science/results/s1/{episode_geodetic.csv, tierA_station_offsets.csv, tierA_control.csv,
tierB_gate.json, tierB_model_vs_obs.csv, tierA_summary.json, s1.log}; figures svi_science/figures/s1_*.png
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
sys.path.insert(0, f"{ROOT}/lfe_extension/scripts")
from sse_strain import load_fault_ascii, load_solution          # noqa: E402
from sse_strain.geodesy import llh2localxy                       # noqa: E402
from sse_strain.matio import load_mat                            # noqa: E402
import slip_build as SB                                          # noqa: E402
import slip_resolution as SR                                     # noqa: E402

OUT = f"{ROOT}/svi_science/results/s1"
FIG = f"{ROOT}/svi_science/figures"
SC = f"{ROOT}/lfe_extension/results/slip_compare"
EP_SVI = f"{ROOT}/lfe_extension_wa/results/episode_rates/episodes_svi_final.csv"
EP_NW = f"{ROOT}/lfe_extension_wa/results/episode_rates/episodes_nw_final.csv"
FRAME = json.load(open(f"{ROOT}/svi_science/results/s0/frame_groups.json"))
FAMS = f"{ROOT}/svi_science/results/s0/families.csv"

SVI14 = ["ALBH", "PGC5", "SC02", "SC04", "NEAH", "PTRF", "PTAA", "NANO", "P435", "P436", "P439", "CLRS", "COUP", "SC03"]
WA = ["CHCM", "KTBW", "P426", "ELSR", "P437", "CUSH", "P419", "PFLD", "SSHO", "SEAT", "SMAI", "P399"]
B3_STA = ["ALBH", "PGC5", "SC04"]
FIX_AZ = 236.0
SLAB_AZ = (FRAME["mean_downdip_az"] + 180) % 360
MU = 40e9
NBOOT = 1000
RNG = np.random.default_rng(20261007)
CTRL_START = np.datetime64("2014-01-15")
LOG = open(f"{OUT}/s1.log", "a")


def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


def along_strike(lat, lon):
    x, y = llh2localxy(np.asarray(lat), np.asarray(lon), (FRAME["origin_lon"], FRAME["origin_lat"]))
    s = np.deg2rad(FRAME["strike_az"])
    return x * np.sin(s) + y * np.cos(s)


# ------------------------------------------------------------------------------------------------ data
class GNSS:
    def __init__(self):
        X = load_mat(f"{SB.SOLDIR}/X.mat")["Xd"]
        self.dates_raw = SB.fixed_dates(np.asarray(X["timeline"], float).ravel())
        names = np.asarray(X["name"]).astype(str)
        ts = np.asarray(X["ts"], float); m = np.asarray(X["data_mask"]).astype(bool)
        # full daily calendar (one day, 2009-02-24, missing in the solution -> NaN)
        self.cal = np.arange(self.dates_raw[0], self.dates_raw[-1] + 1)
        ix = (self.dates_raw - self.cal[0]).astype(int)
        self.ts = np.full((ts.shape[0], self.cal.size), np.nan); self.ts[:, ix] = np.where(m, ts, np.nan)
        self.mask = np.zeros((ts.shape[0], self.cal.size), bool); self.mask[:, ix] = m
        self.ix_raw = ix
        self.names = names
        self.sta = np.array([n[:-1] for n in names[::3]])
        llh = np.asarray(X["llh"], float)[::3]
        self.lon, self.lat = llh[:, 0], llh[:, 1]
        self.M5 = self.median5(self.ts)

    @staticmethod
    def median5(Y):
        """M5[:, i] = median over calendar days i..i+4 (>= 3 valid), else NaN."""
        n = Y.shape[1]
        W = np.lib.stride_tricks.sliding_window_view(np.concatenate([Y, np.full((Y.shape[0], 4), np.nan)], 1), 5, axis=1)
        cnt = np.isfinite(W).sum(-1)
        with np.errstate(all="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                med = np.nanmedian(W, axis=-1)
        med[cnt < 3] = np.nan
        return med[:, :n]

    def di(self, d):
        return int((np.datetime64(pd.Timestamp(d).date(), "D") - self.cal[0]).astype(int))

    def rows(self, stas):
        idx = [int(np.flatnonzero(self.sta == s)[0]) for s in stas]
        return idx


def offsets(M5, rows_e, rows_n, i_start, i_end):
    """Offsets (n_sta, 2) [E, N] for scalar indices or (n_sta, n_win, 2) for index arrays."""
    pre, post = np.asarray(i_start) - 6, np.asarray(i_end) + 2
    E = M5[rows_e][:, post] - M5[rows_e][:, pre]
    N = M5[rows_n][:, post] - M5[rows_n][:, pre]
    return np.stack([E, N], -1)


def unit(az):
    a = np.deg2rad(az)
    return np.array([np.sin(a), np.cos(a)])


# ------------------------------------------------------------------------------------------------ Tier A
def tier_a(g, stas, E, exclude, top_names, label):
    """E: episode table (start, end datetime64).  exclude: list of (a, b) datetime64 intervals.
    Returns per-episode DataFrame, station offsets long table, control long table, directions."""
    r = np.array(g.rows(stas)); re, rn = 3 * r, 3 * r + 1
    i0 = np.array([g.di(s) for s in E.start]); i1 = np.array([g.di(e) for e in E.end])
    O = offsets(g.M5, re, rn, i0, i1)                      # (n_sta, n_ep, 2)
    ntop = E.episode.isin(top_names).to_numpy()
    top_idx = np.flatnonzero(ntop)
    # own direction: mean offset over top episodes (LOO for those)
    def dirs(excl=None):
        use = [k for k in top_idx if k != excl]
        v = np.nanmean(O[:, use, :], axis=1)
        return v / np.linalg.norm(v, axis=1, keepdims=True)
    U_all = dirs()
    az_own = (np.degrees(np.arctan2(U_all[:, 0], U_all[:, 1])) + 360) % 360
    # template pattern (projected, own direction, all top)
    P_tmp = np.nanmean((O[:, top_idx, :] * U_all[:, None, :]).sum(-1), axis=1)
    # exclusion -> valid pseudo-start mask per span
    cal = g.cal
    bad = np.zeros(cal.size + 64, bool)
    for a, b in exclude:
        ia, ib = max(g.di(a), 0), min(g.di(b), cal.size - 1)
        if ib >= ia:
            bad[ia:ib + 1] = True
    cs = np.concatenate([[0], np.cumsum(bad)])
    doy = (cal - cal.astype("datetime64[Y]")).astype(int)
    rows, srows, crows = [], [], []
    for k, e in enumerate(E.itertuples()):
        Uk = dirs(k) if ntop[k] else U_all
        variants = {"own": Uk, "fix": np.tile(unit(FIX_AZ), (len(stas), 1)), "slab": np.tile(unit(SLAB_AZ), (len(stas), 1))}
        span = i1[k] - i0[k]
        # control starts
        s_all = np.arange(g.di(CTRL_START) + 6, cal.size - span - 7)
        lo, hi = s_all - 6, s_all + span + 6
        clean = (cs[hi + 1] - cs[lo]) == 0
        dd = np.abs(doy[s_all] - doy[i0[k]]); dd = np.minimum(dd, 365 - dd)
        widened = "no"
        sel = clean & (dd <= 60)
        if sel.sum() < 30:
            sel = clean & (dd <= 120); widened = "doy120"
        if sel.sum() < 30:
            sel = clean; widened = "noseason"
        sc = s_all[sel]
        Oc = offsets(g.M5, re, rn, sc, sc + span)           # (n_sta, n_c, 2)
        rec = dict(episode=e.episode, start=str(e.start)[:10], end=str(e.end)[:10], n_ctrl=int(sel.sum()),
                   ctrl_widened=widened)
        ok = np.isfinite(O[:, k, :]).all(-1)
        rec["n_sta"] = int(ok.sum())
        for vn, U in variants.items():
            x = (O[:, k, :] * U).sum(-1)                     # projected (n_sta,)
            D = np.nanmean(x) if ok.any() else np.nan
            xc = (Oc * U[:, None, :]).sum(-1)               # (n_sta, n_c)
            Dc = np.nanmean(np.where(np.isfinite(xc), xc, np.nan), axis=0)
            # control uses the same station set as the episode (missing stations dropped likewise)
            xc_same = np.where(ok[:, None], xc, np.nan)
            Dc = np.nanmean(xc_same, axis=0) if ok.any() else Dc
            Dc = Dc[np.isfinite(Dc)]
            bias, sig = float(np.median(Dc)), float(np.std(Dc))
            sfx = "" if vn == "own" else f"_{vn}"
            rec[f"D_raw{sfx}"] = D
            rec[f"bias{sfx}"] = bias
            rec[f"sigma{sfx}"] = sig
            rec[f"D{sfx}"] = D - bias
            rec[f"SN{sfx}"] = (D - bias) / sig
            rec[f"SN_uncorr{sfx}"] = D / sig
            if vn == "own":
                xs = x[ok]
                if xs.size:
                    bs = xs[RNG.integers(0, xs.size, (NBOOT, xs.size))].mean(1)
                    rec["sigma_boot"] = float(bs.std())
                rec["ctrl_p99"] = float(np.percentile(Dc - bias, 99))
                pt = P_tmp[ok]
                rec["amp_template"] = float((pt * xs).sum() / (pt ** 2).sum()) if xs.size else np.nan
                for s, xv, okv, u in zip(stas, x, ok, U):
                    srows.append(dict(set=label, episode=e.episode, station=s, E_mm=O[stas.index(s), k, 0],
                                      N_mm=O[stas.index(s), k, 1], proj_own_mm=xv, used=bool(okv),
                                      dir_own_az=float((np.degrees(np.arctan2(u[0], u[1])) + 360) % 360)))
                crows.append(dict(set=label, episode=e.episode, n_ctrl=int(Dc.size), bias=bias, sigma=sig,
                                  q01=float(np.percentile(Dc, 1)), q99=float(np.percentile(Dc, 99))))
                rec["_ctrl_mean_series"] = Dc
        rows.append(rec)
    R = pd.DataFrame(rows)
    return R, pd.DataFrame(srows), pd.DataFrame(crows), dict(zip(stas, np.round(az_own, 1)))


def subset_SN(g, stas, E, exclude, top_names, label):
    R, _, _, _ = tier_a(g, stas, E, exclude, top_names, label)
    return R[["episode", "D", "sigma", "SN", "n_sta"]]


# ------------------------------------------------------------------------------------------------ main
def main():
    os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)
    log("=== S1 start ===")
    g = GNSS()
    log(f"GNSS: {g.sta.size} stations, calendar {g.cal[0]}..{g.cal[-1]} ({g.cal.size} d)")
    for s in SVI14 + WA:
        assert s in g.sta, s

    ES = pd.read_csv(EP_SVI, parse_dates=["start", "end"])
    EN = pd.read_csv(EP_NW, parse_dates=["start", "end"])
    for D in (ES, EN):
        D["start"] = D.start.values.astype("datetime64[D]"); D["end"] = D.end.values.astype("datetime64[D]")
    pad = np.timedelta64(10, "D")
    excl = [(a - pad, b + pad) for a, b in zip(ES.start.values, ES.end.values)] + \
           [(a - pad, b + pad) for a, b in zip(EN[EN.major].start.values, EN[EN.major].end.values)]
    fm = ES[ES.in_frozen_major]
    top_svi = fm.sort_values("family_exposure_days", ascending=False).episode.head(10).tolist()
    log(f"10 largest frozen-major SVI episodes (by family_exposure_days): {top_svi}")

    # station geometry
    st_as = dict(zip(g.sta, along_strike(g.lat, g.lon)))
    north = [s for s in SVI14 if st_as[s] > 0]; south = [s for s in SVI14 if st_as[s] <= 0]
    log(f"north subset {north}; south subset {south}")

    # ---------------- Tier A (SVI-14)
    R, S, C, az = tier_a(g, SVI14, ES, excl, top_svi, "SVI14")
    log(f"own ETS azimuths: {az}")
    Rn, _, _, _ = tier_a(g, north, ES, excl, top_svi, "north")
    Rs, _, _, _ = tier_a(g, south, ES, excl, top_svi, "south")
    for nm, X in (("north", Rn), ("south", Rs)):
        R[f"D_{nm}"] = X.D.values; R[f"sigma_{nm}"] = X.sigma.values; R[f"SN_{nm}"] = X.SN.values
        R[f"n_sta_{nm}"] = X.n_sta.values
    ctrl_series = dict(zip(R.episode, R["_ctrl_mean_series"]))
    R = R.drop(columns="_ctrl_mean_series")
    R = ES[["episode", "major", "in_frozen_major", "span_days", "family_exposure_days"]].merge(R, on="episode")
    # close pairs
    st_, en_ = ES.start.values, ES.end.values
    gp = (st_[1:] - en_[:-1]).astype("timedelta64[D]").astype(int)
    gap_prev = np.r_[np.inf, gp]
    gap_next = np.r_[gp, np.inf]
    R["gap_prev_d"] = gap_prev; R["gap_next_d"] = gap_next
    R["close_pair"] = (gap_prev < 30) | (gap_next < 30)
    R["window_overlaps_other_ep"] = (gap_prev < 7) | (gap_next < 7)
    R["dir_top10"] = R.episode.isin(top_svi)
    S.to_csv(f"{OUT}/tierA_station_offsets.csv", index=False)
    C.to_csv(f"{OUT}/tierA_control.csv", index=False)

    fmR = R[R.in_frozen_major]
    nA = int((fmR.SN > 3).sum()); nA_fix = int((fmR.SN_fix > 3).sum()); nA_slab = int((fmR.SN_slab > 3).sum())
    def verdict(n):
        return "ACCEPT" if n >= 15 else ("FAIL" if n < 8 else "PARTIAL")
    vA = verdict(nA); vA_fix = verdict(nA_fix)
    governing = vA if vA == vA_fix else vA_fix
    sumA = dict(n_frozen_major=int(len(fmR)), detected_own=nA, detected_fixed236=nA_fix, detected_slab231=nA_slab,
                verdict_own=vA, verdict_fixed=vA_fix, verdict_governing=governing,
                detected_all_major_own=int((R[R.major].SN > 3).sum()), n_major=int(R.major.sum()),
                detected_minor_own=int((R[~R.major].SN > 3).sum()), n_minor=int((~R.major).sum()),
                detected_minor_fixed=int((R[~R.major].SN_fix > 3).sum()),
                median_sigma_mm=float(R.sigma.median()), median_bias_mm=float(R.bias.median()),
                median_sigma_boot_frozen_major_mm=float(fmR.sigma_boot.median()),
                frozen_major_D_mm_range=[float(fmR.D.min()), float(fmR.D.median()), float(fmR.D.max())],
                own_azimuths=az, top10=top_svi, north=north, south=south,
                ctrl_widened=R.ctrl_widened.value_counts().to_dict())
    log(f"TIER A: {json.dumps(sumA, default=str)}")

    # ---------------- leakage: WA set
    EN14 = EN[(EN.kind == "pnsn") & (EN.start >= np.datetime64("2014-01-01"))]
    top_nw = EN14[EN14.in_frozen_major].sort_values("family_exposure_days", ascending=False).episode.head(10).tolist()
    log(f"10 largest NW frozen-major 2014+: {top_nw}")
    RW_nw, _, _, az_wa = tier_a(g, WA, EN14[EN14.in_frozen_major].reset_index(drop=True), excl, top_nw, "WA_on_NW")
    norm_wa = float(RW_nw[RW_nw.SN > 3].D.median())
    norm_svi = float(fmR[fmR.SN > 3].D.median())
    log(f"WA own azimuths {az_wa}; norm WA (NW frozen majors S/N>3, n={int((RW_nw.SN > 3).sum())}) = {norm_wa:.2f} mm; "
        f"norm SVI = {norm_svi:.2f} mm")
    # WA Tier A on SVI episode windows: directions from NW top10 -> pass a combined table so dirs come from NW eps
    Ecomb = pd.concat([EN14[EN14.in_frozen_major][["episode", "start", "end"]].assign(episode=lambda d: "NW:" + d.episode),
                       ES[["episode", "start", "end"]]], ignore_index=True)
    RW, _, _, _ = tier_a(g, WA, Ecomb, excl, ["NW:" + t for t in top_nw], "WA_on_SVI")
    RW = RW[~RW.episode.str.startswith("NW:")].set_index("episode")
    ov = []
    for e in ES.itertuples():
        m = (EN.start.values <= e.end) & (EN.end.values >= e.start)
        ov.append(";".join(EN.episode[m].tolist()))
    R["nw_overlap"] = ov
    R["D_WA"] = RW.D.reindex(R.episode).values; R["sigma_WA"] = RW.sigma.reindex(R.episode).values
    R["SN_WA"] = RW.SN.reindex(R.episode).values; R["n_sta_WA"] = RW.n_sta.reindex(R.episode).values
    R["Dnorm_SVI"] = R.D / norm_svi; R["Dnorm_WA"] = R.D_WA / norm_wa
    has = R.nw_overlap != ""
    R["wa_leak_possible"] = has & (R.SN_WA >= 3) & (R.Dnorm_WA >= 2 * R.Dnorm_SVI)
    R["south_only"] = (R.SN >= 3) & (R.SN_north < 2) & (R.SN_south >= 3)

    # ---------------- Tier B
    C8 = np.load(f"{SC}/slip_cache.npz")
    G, Lm, Ls, Ld, A = C8["G"], C8["L"], C8["Ls"], C8["Ld"], C8["A"]
    mdates = C8["dates"].astype("datetime64[D]")
    assert (mdates == g.dates_raw).all(), "model dates differ from X.mat fixed dates"
    fault = load_fault_ascii(SB.MESH)
    sol = load_solution(SB.SOLDIR, fault); sol.L = Lm
    n_p = fault.n_patch
    meta = pd.read_csv(f"{SC}/slip_svi_family_meta.csv")
    Wz = np.load(f"{SC}/slip_svi_weights.npz")
    w_ic = SR.ic_weights(sol, Ls, Ld, A, meta, Wz["r25"])
    RR, D3 = SR.resolution_mats(sol, G)
    Rbar = sum(wk * r["R"] for wk, r in zip(w_ic, RR))
    Rt = Rbar[n_p:, n_p:]
    del RR
    clon, clat, cdep = fault.centroid_llh[:, 0], fault.centroid_llh[:, 1], -fault.centroid_llh[:, 2]
    F0 = pd.read_csv(FAMS)
    fx, fy = llh2localxy(F0.lat.values, F0.lon.values, (FRAME["origin_lon"], FRAME["origin_lat"]))
    px, py = llh2localxy(clat, clon, (FRAME["origin_lon"], FRAME["origin_lat"]))
    dmin = np.sqrt((px[:, None] - fx[None]) ** 2 + (py[:, None] - fy[None]) ** 2).min(1)
    setA = np.array(sorted(meta.patch.unique()))
    setB = np.flatnonzero(dmin <= 15.0)
    p_as = along_strike(clat, clon)
    area = fault.area
    trA, trB = float(np.diag(Rt)[setA].sum()), float(np.diag(Rt)[setB].sum())
    trB_ssds = float(np.diag(Rbar)[setB].sum() + np.diag(Rbar)[n_p + setB].sum())
    halves = {}
    med_as = np.median(p_as[setB]); med_z = np.median(cdep[setB])
    halves["north"] = setB[p_as[setB] > med_as]; halves["south"] = setB[p_as[setB] <= med_as]
    halves["updip"] = setB[cdep[setB] <= med_z]; halves["downdip"] = setB[cdep[setB] > med_z]
    other = {"north": "south", "south": "north", "updip": "downdip", "downdip": "updip"}
    spikes, rec_maps = {}, {}
    for h, P in halves.items():
        m = np.zeros(n_p); m[P] = 1.0
        rec = Rt @ m
        rec_maps[h] = rec
        Q = halves[other[h]]
        r_in = float((rec[P] * area[P]).sum() / area[P].sum())
        r_out = float((rec[Q] * area[Q]).sum() / area[Q].sum())
        spikes[h] = dict(n_patches=int(P.size), recovery=round(r_in, 3), other_half_mean=round(r_out, 3),
                         crosstalk=round(r_out / r_in, 3),
                         potency_recovered_total=round(float((rec * area).sum() / (m * area).sum()), 3))
    ns_ok = all(spikes[h]["recovery"] >= 0.5 and spikes[h]["crosstalk"] < 0.3 for h in ("north", "south"))
    log(f"TIER B: trace setA({setA.size})={trA:.3f} setB({setB.size})={trB:.3f} (ss+ds {trB_ssds:.3f}); spikes {spikes}")

    # B3 modelled vs observed at ALBH/PGC5/SC04 for top10
    GL = G @ Lm                                               # mm per unit amplitude (data unit)
    Dmod = np.full((G.shape[0], g.cal.size), np.nan)
    Dmod[:, g.ix_raw] = GL @ A.T
    rows3 = np.array(g.rows(B3_STA))
    sub = np.sort(np.r_[3 * rows3, 3 * rows3 + 1])
    Dmod_m = np.where(g.mask[sub], Dmod[sub], np.nan)
    M5o = g.M5[sub]; M5m = GNSS.median5(Dmod_m)
    rowmap = {r: i for i, r in enumerate(sub)}
    b3 = []
    Sown = S.set_index(["episode", "station"])
    for ep in top_svi:
        e = ES[ES.episode == ep].iloc[0]
        a0, a1 = g.di(e.start) - 6, g.di(e.end) + 2
        for s, r in zip(B3_STA, rows3):
            ie, inn = rowmap[3 * r], rowmap[3 * r + 1]
            oE, oN = M5o[ie, a1] - M5o[ie, a0], M5o[inn, a1] - M5o[inn, a0]
            mE, mN = M5m[ie, a1] - M5m[ie, a0], M5m[inn, a1] - M5m[inn, a0]
            azs = Sown.loc[(ep, s), "dir_own_az"]; u = unit(azs)
            b3.append(dict(episode=ep, station=s, obs_E=oE, obs_N=oN, mod_E=mE, mod_N=mN,
                           obs_proj=oE * u[0] + oN * u[1], mod_proj=mE * u[0] + mN * u[1], dir_az=azs))
    B3 = pd.DataFrame(b3)
    B3.to_csv(f"{OUT}/tierB_model_vs_obs.csv", index=False)
    ok3 = B3.dropna()
    slope = float((ok3.mod_proj * ok3.obs_proj).sum() / (ok3.obs_proj ** 2).sum())
    oEN = np.r_[ok3.obs_E, ok3.obs_N]; mEN = np.r_[ok3.mod_E, ok3.mod_N]
    slope_EN = float((mEN * oEN).sum() / (oEN ** 2).sum())
    r_proj = float(np.corrcoef(ok3.obs_proj, ok3.mod_proj)[0, 1])
    slope_sta = {s: float((d.mod_proj * d.obs_proj).sum() / (d.obs_proj ** 2).sum()) for s, d in ok3.groupby("station")}
    # bootstrap over episodes for slope uncertainty
    eps_u = ok3.episode.unique(); bsl = []
    for _ in range(1000):
        pick = RNG.choice(eps_u, eps_u.size)
        d = pd.concat([ok3[ok3.episode == p] for p in pick])
        bsl.append((d.mod_proj * d.obs_proj).sum() / (d.obs_proj ** 2).sum())
    slope_ci = [float(np.percentile(bsl, 2.5)), float(np.percentile(bsl, 97.5))]
    slope_ok = 0.7 <= slope <= 1.3
    log(f"B3 slope proj {slope:.3f} CI {slope_ci} (EN {slope_EN:.3f}; r {r_proj:.3f}) per station {slope_sta}")

    # B4 per-episode model slip / M0 (patch thrust, signed)
    T = np.full((n_p, g.cal.size), np.nan)
    T[:, g.ix_raw] = -(Ld @ A.T) * 1e3                         # thrust-positive mm
    T = pd.DataFrame(T.T).ffill().to_numpy().T                # fill the one missing day
    M5T = GNSS.median5(T)
    nc = (clat >= 47.0) & (clat <= 50.5)
    nP, sP = halves["north"], halves["south"]
    # family-mean r25 along-rake (corrected parquet)
    FD = pd.read_parquet(f"{SC}/slip_svi_family_daily.parquet", columns=["fid", "date", "slip_mm_r25"])
    FD["date"] = FD.date.astype("datetime64[ns]")
    Y = FD.pivot(index="date", columns="fid", values="slip_mm_r25")
    Y = Y.reindex(pd.DatetimeIndex(g.cal.astype("datetime64[ns]"))).ffill()
    fam_as = F0.set_index("fid").along_strike_km.reindex(Y.columns)
    yreg = Y.mean(1).to_numpy()[None]; yN = Y.loc[:, (fam_as > 0).values].mean(1).to_numpy()[None]
    yS = Y.loc[:, (fam_as <= 0).values].mean(1).to_numpy()[None]
    M5y = {k: GNSS.median5(v)[0] for k, v in (("reg", yreg), ("N", yN), ("S", yS))}
    b4 = []
    for e in ES.itertuples():
        a0, a1 = g.di(e.start) - 6, g.di(e.end) + 2
        inc = M5T[:, a1] - M5T[:, a0]                            # mm
        def m0(P):
            return float(MU * (inc[P] * 1e-3 * area[P] * 1e6).sum())
        M0b, M0nc = m0(setB), m0(np.flatnonzero(nc))
        M0N, M0S = m0(nP), m0(sP)
        mw = lambda M: 2 / 3 * (np.log10(M) - 9.1) if M > 0 else np.nan
        b4.append(dict(episode=e.episode,
                       tierB_slip_mm=float((inc[setB] * area[setB]).sum() / area[setB].sum()),
                       tierB_M0_Nm=M0b, tierB_Mw=mw(M0b), tierB_M0_ncasc_Nm=M0nc, tierB_Mw_ncasc=mw(M0nc),
                       tierB_north_share=(M0N / (M0N + M0S)) if (M0N > 0 and M0S > 0) else np.nan,
                       tierB_M0_north=M0N, tierB_M0_south=M0S,
                       tierB_fam_r25_mm=float(M5y["reg"][a1] - M5y["reg"][a0]),
                       tierB_fam_r25_north_mm=float(M5y["N"][a1] - M5y["N"][a0]),
                       tierB_fam_r25_south_mm=float(M5y["S"][a1] - M5y["S"][a0])))
    B4 = pd.DataFrame(b4)
    R = R.merge(B4, on="episode")
    # Tier A vs Tier B agreement (not a gate): correlation over frozen majors
    fmm = R[R.in_frozen_major]
    rho_AB = float(fmm[["D", "tierB_slip_mm"]].corr(method="spearman").iloc[0, 1])
    rho_AB_all = float(R[["D", "tierB_slip_mm"]].corr(method="spearman").iloc[0, 1])

    # B5 IC26 / IC17 share 2024-2026
    th = np.deg2rad(meta.rake_deg.to_numpy())[:, None]
    pal = (Wz["r25"] @ Ls) * np.cos(th) + (-(Wz["r25"] @ Ld)) * np.sin(th)       # (n_fam, K) m/amp
    pal_reg = pal.mean(0) * 1e3                                                   # mm/amp
    Cc = A * pal_reg[None]                                                        # (n_t, K) contribution
    chk = float(np.nanmax(np.abs(Cc.sum(1) - Y.mean(1).to_numpy()[g.ix_raw])))
    tt = g.dates_raw
    w24 = tt >= np.datetime64("2024-01-01")
    tot = Cc[w24].sum(1)
    ics = [f"IC{int(i) + 1}" for i in sol.component_index]
    cov_share = {ic: float(np.cov(Cc[w24, k], tot)[0, 1] / tot.var()) for k, ic in enumerate(ics)}
    var_share = {ic: float(Cc[w24, k].var() / Cc[w24].var(0).sum()) for k, ic in enumerate(ics)}
    cov_share_all = {ic: float(np.cov(Cc[:, k], Cc.sum(1))[0, 1] / Cc.sum(1).var()) for k, ic in enumerate(ics)}
    # increments of each IC contribution in 2024+ episodes
    Ccal = np.full((Cc.shape[1], g.cal.size), np.nan); Ccal[:, g.ix_raw] = Cc.T
    M5c = GNSS.median5(Ccal)
    inc_ic = []
    for e in ES[ES.start >= np.datetime64("2024-01-01")].itertuples():
        a0, a1 = g.di(e.start) - 6, g.di(e.end) + 2
        inc_ic.append(M5c[:, a1] - M5c[:, a0])
    inc_ic = np.array(inc_ic)
    abs_share = {ic: float(np.nansum(np.abs(inc_ic[:, k])) / np.nansum(np.abs(inc_ic))) for k, ic in enumerate(ics)}
    ep24 = ES[ES.start >= np.datetime64("2024-01-01")].episode.tolist()
    ic26_inc = dict(zip(ep24, np.round(inc_ic[:, ics.index("IC26")], 2).tolist()))
    tot_inc = dict(zip(ep24, np.round(np.nansum(inc_ic, 1), 2).tolist()))
    log(f"B5 check (sum of IC contributions vs parquet family mean) max |diff| {chk:.2e} mm")
    log(f"B5 cov share 2024+: {json.dumps({k: round(v, 3) for k, v in cov_share.items()})}")

    gateB = dict(
        B1_trace_setB=round(trB, 3), B1_trace_setA=round(trA, 3), B1_trace_setB_ss_plus_ds=round(trB_ssds, 3),
        n_setA=int(setA.size), n_setB=int(setB.size),
        setB_rule="patch centroid within 15 km (horizontal) of any SVI family",
        B2_spikes=spikes, B2_split=dict(along_strike_median_km=round(float(med_as), 2), depth_median_km=round(float(med_z), 2)),
        B3_slope_proj=round(slope, 3), B3_slope_proj_CI95_episode_boot=[round(x, 3) for x in slope_ci],
        B3_slope_EN=round(slope_EN, 3), B3_corr_proj=round(r_proj, 3),
        B3_slope_per_station={k: round(v, 3) for k, v in slope_sta.items()}, B3_n_points=int(len(ok3)),
        B3_pass=bool(slope_ok),
        B4_note="per-episode model Mw in episode_geodetic.csv (tierB_Mw = SVI patch set B, tierB_Mw_ncasc = 47-50.5 N)",
        B5_ic26_cov_share_2024_2026=round(cov_share["IC26"], 3), B5_ic17_cov_share_2024_2026=round(cov_share["IC17"], 3),
        B5_cov_share_2024_2026={k: round(v, 3) for k, v in cov_share.items()},
        B5_var_share_2024_2026={k: round(v, 3) for k, v in var_share.items()},
        B5_cov_share_2007_2026={k: round(v, 3) for k, v in cov_share_all.items()},
        B5_abs_increment_share_2024_eps={k: round(v, 3) for k, v in abs_share.items()},
        B5_ic26_increment_mm_by_2024_episode=ic26_inc, B5_total_increment_mm_by_2024_episode=tot_inc,
        B5_reconstruction_check_max_abs_mm=chk,
        tierA_vs_tierB_spearman_frozen_major=round(rho_AB, 3), tierA_vs_tierB_spearman_all79=round(rho_AB_all, 3),
        verdict_regional="PASS" if slope_ok else "FAIL",
        verdict_along_strike="PASS" if (slope_ok and trB >= 2 and ns_ok) else "FAIL",
        along_strike_components=dict(trace_ge_2=bool(trB >= 2), ns_spikes_ok=bool(ns_ok), slope_ok=bool(slope_ok)))
    R["tierB_usable_regional"] = gateB["verdict_regional"] == "PASS"
    R["tierB_usable_along_strike"] = gateB["verdict_along_strike"] == "PASS"
    if gateB["verdict_along_strike"] != "PASS":
        R["tierB_north_share_gated"] = np.nan
    else:
        R["tierB_north_share_gated"] = R.tierB_north_share
    log(f"TIER B gate: {json.dumps({k: v for k, v in gateB.items() if not k.startswith('B5_')}, default=str)}")

    # flags
    fl = []
    for r in R.itertuples():
        f = []
        if r.close_pair: f.append("close_pair")
        if r.window_overlaps_other_ep: f.append("window_overlaps_other_episode")
        if r.dir_top10: f.append("direction_episode_LOO")
        if r.ctrl_widened != "no": f.append(f"ctrl_{r.ctrl_widened}")
        if r.wa_leak_possible: f.append("wa_leak_possible")
        if r.south_only: f.append("south_only")
        if r.nw_overlap: f.append("overlaps_NW")
        if r.episode.startswith("2017-03"): f.append("detection_problem_2017-03")
        if r.episode.startswith("2025-09") or r.episode.startswith("2025-10"): f.append("detection_problem_2025-10")
        if (r.SN > 3) != (r.SN_fix > 3): f.append("own_vs_fixed_disagree")
        if r.n_sta < 10: f.append("n_sta_lt_10")
        fl.append(";".join(f))
    R["flags"] = fl
    R["detected_own"] = R.SN > 3; R["detected_fixed"] = R.SN_fix > 3
    R.to_csv(f"{OUT}/episode_geodetic.csv", index=False, float_format="%.4g")
    leak = R[R.nw_overlap != ""][["episode", "major", "in_frozen_major", "nw_overlap", "D", "SN", "SN_north", "SN_south",
                                   "D_WA", "SN_WA", "Dnorm_SVI", "Dnorm_WA", "wa_leak_possible", "south_only"]]
    leak.to_csv(f"{OUT}/leakage_svi_vs_wa.csv", index=False, float_format="%.3g")
    sumA.update(dict(norm_wa_mm=norm_wa, norm_svi_mm=norm_svi, wa_azimuths=az_wa, top10_nw=top_nw,
                     n_overlap_NW=int((R.nw_overlap != "").sum()),
                     n_wa_leak_possible=int(R.wa_leak_possible.sum()), wa_leak_episodes=R[R.wa_leak_possible].episode.tolist(),
                     n_south_only=int(R.south_only.sum()), south_only_episodes=R[R.south_only].episode.tolist(),
                     detected_frozen_north=int((fmR.episode.isin(R[R.SN_north > 3].episode)).sum()),
                     detected_frozen_south=int((fmR.episode.isin(R[R.SN_south > 3].episode)).sum()),
                     wa_nw_frozen_major_detected=int((RW_nw.SN > 3).sum()), wa_nw_frozen_major_n=int(len(RW_nw))))
    json.dump(sumA, open(f"{OUT}/tierA_summary.json", "w"), indent=1, default=str)
    json.dump(gateB, open(f"{OUT}/tierB_gate.json", "w"), indent=1, default=str)
    np.savez_compressed(f"{OUT}/tierB_spike_maps.npz", setA=setA, setB=setB, diag=np.diag(Rt),
                        **{f"rec_{h}": v for h, v in rec_maps.items()}, **{f"in_{h}": v for h, v in halves.items()})

    figures(g, R, B3, gateB, fault, rec_maps, halves, setB, F0, ctrl_series)
    log("=== S1 done ===")


# ------------------------------------------------------------------------------------------------ figures
def figures(g, R, B3, gateB, fault, rec_maps, halves, setB, F0, ctrl_series):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    C_MAJ, C_OTH, C_MIN = "#1f5fa8", "#6aa3d8", "#9a9a9a"
    t = pd.to_datetime(R.start) + (pd.to_datetime(R.end) - pd.to_datetime(R.start)) / 2

    # Fig 1: D per episode with sigma over time
    fig, axs = plt.subplots(2, 1, figsize=(11, 6.6), sharex=True, gridspec_kw=dict(height_ratios=[2.2, 1]))
    ax = axs[0]
    for m, c, lab, ms in ((R.in_frozen_major, C_MAJ, "frozen major (acceptance set)", 6),
                          (R.major & ~R.in_frozen_major, C_OTH, "other major", 5), (~R.major, C_MIN, "minor", 4)):
        ax.errorbar(t[m], R.D[m], yerr=R.sigma[m], fmt="o", color=c, ms=ms, lw=1, capsize=0, label=lab)
    ax.scatter(t, 3 * R.sigma, marker="_", s=60, color="#c0392b", lw=1.2, label="3σ of random-window control", zorder=5)
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_ylabel("ETS-direction offset D (mm)\nmean of SVI-14, drift-corrected")
    ax.legend(frameon=False, ncol=2, loc="upper left", fontsize=8)
    nA = int((R[R.in_frozen_major].SN > 3).sum())
    ax.set_title(f"S1 Tier A — GNSS episode size at the 14 SVI stations: {nA} of {int(R.in_frozen_major.sum())} "
                 f"frozen-major episodes above 3σ (rule: ≥ 15 accept, < 8 fail)", fontsize=10, loc="left")
    ax = axs[1]
    for m, c in ((R.in_frozen_major, C_MAJ), (R.major & ~R.in_frozen_major, C_OTH), (~R.major, C_MIN)):
        ax.plot(t[m], R.SN[m], "o", color=c, ms=4)
    ax.axhline(3, color="#c0392b", lw=1, ls="--")
    ax.set_ylabel("S/N (drift-corrected)"); ax.set_yscale("symlog", linthresh=5)
    ax.set_yticks([-3, 0, 3, 5, 10, 20]); ax.set_yticklabels(["−3", "0", "3", "5", "10", "20"])
    fig.savefig(f"{FIG}/s1_episode_size.png", dpi=140, bbox_inches="tight"); plt.close(fig)

    # Fig 2: modelled vs observed
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.6))
    cols = {"ALBH": "#1f5fa8", "PGC5": "#d9822b", "SC04": "#2a9d6f"}
    ax = axs[0]
    for s, d in B3.groupby("station"):
        ax.plot(d.obs_proj, d.mod_proj, "o", color=cols[s], label=s, ms=6)
    lim = np.nanmax(np.abs(np.r_[B3.obs_proj, B3.mod_proj])) * 1.1
    xx = np.array([-lim, lim])
    ax.fill_between(xx, 0.7 * xx, 1.3 * xx, color="0.85", lw=0, label="slope 0.7–1.3")
    ax.plot(xx, xx, color="0.4", lw=0.8)
    ax.plot(xx, gateB["B3_slope_proj"] * xx, color="#c0392b", lw=1.2, label=f"fit slope {gateB['B3_slope_proj']:.2f}")
    ax.set_xlim(-lim * 0.2, lim); ax.set_ylim(-lim * 0.2, lim)
    ax.set_xlabel("observed offset, own ETS direction (mm)"); ax.set_ylabel("modelled offset (mm)")
    ax.legend(frameon=False, fontsize=8); ax.set_title("10 largest SVI episodes, projected", fontsize=10, loc="left")
    ax = axs[1]
    for s, d in B3.groupby("station"):
        ax.plot(d.obs_E, d.mod_E, "o", color=cols[s], ms=5)
        ax.plot(d.obs_N, d.mod_N, "s", mfc="none", color=cols[s], ms=5)
    lim2 = np.nanmax(np.abs(np.r_[B3.obs_E, B3.obs_N, B3.mod_E, B3.mod_N])) * 1.1
    ax.plot([-lim2, lim2], [-lim2, lim2], color="0.4", lw=0.8)
    ax.set_xlabel("observed (mm): ● east, □ north"); ax.set_ylabel("modelled (mm)")
    ax.set_title(f"components, slope {gateB['B3_slope_EN']:.2f}", fontsize=10, loc="left")
    fig.savefig(f"{FIG}/s1_model_vs_obs.png", dpi=140, bbox_inches="tight"); plt.close(fig)

    # Fig 3: spike maps
    n_p = fault.n_patch
    verts = np.column_stack([fault.lon.ravel(), fault.lat.ravel()])
    uniq, inv = np.unique(np.round(verts, 6), axis=0, return_inverse=True)
    trg = mtri.Triangulation(uniq[:, 0], uniq[:, 1], inv.reshape(n_p, 3))
    fig, axs = plt.subplots(1, 4, figsize=(18, 5.6), sharey=True)
    sp = gateB["B2_spikes"]
    for ax, h in zip(axs, ("north", "south", "updip", "downdip")):
        tp = ax.tripcolor(trg, facecolors=rec_maps[h], cmap="RdBu_r", vmin=-1, vmax=1)
        P = halves[h]
        for p in P:
            ax.fill(fault.lon[p], fault.lat[p], fill=False, ec="k", lw=0.8)
        ax.plot(F0.lon, F0.lat, ".", ms=2.5, color="0.15")
        ax.plot(g.lon, g.lat, "^", ms=5, mfc="none", mec="#2a9d6f")
        ax.set_xlim(-125.3, -122.2); ax.set_ylim(47.6, 49.4)
        ax.set_aspect(1 / np.cos(np.deg2rad(48.5)))
        ax.set_title(f"input: {h} half ({sp[h]['n_patches']} patches)\n"
                     f"recovered {sp[h]['recovery']:.2f}, cross-talk {sp[h]['crosstalk']:.2f}", fontsize=9)
    fig.colorbar(tp, ax=axs, shrink=0.75, label="recovered thrust per unit input (IC-weighted R)")
    fig.suptitle(f"S1 Tier B spike tests under SVI (gate set: {len(setB)} patches; resolution trace "
                 f"{gateB['B1_trace_setB']:.2f}; rule trace ≥ 2, N/S recovery ≥ 0.5 with cross-talk < 0.3)\n"
                 "unit thrust on the outlined patches; colour = recovered thrust; dots = LFE families; triangles = GNSS", fontsize=10, y=0.98)
    fig.savefig(f"{FIG}/s1_spike_maps.png", dpi=130, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    main()
