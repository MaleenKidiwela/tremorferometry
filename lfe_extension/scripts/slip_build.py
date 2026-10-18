#!/usr/bin/env python
"""slip_build (2026-10-07) -- daily modelled slow slip at every SVI (and WA-NW) LFE family from the
Gualandi (2025) vbICA solution 2026-08-18, rebuilt with the read-only cascadia-sse-strain package.

    python lfe_extension/scripts/slip_build.py

Slip on patch p at epoch t = sum_k L[p,k] * S_k V[t,k]  (13 retained ICs; L from the per-component
Tarantola-Valette inversion with the exponential 21.09 km prior, sigma0 per IC from misfit_comps CV).
Slip is RELATIVE TO THE SECULAR TREND (no loading term), mean-centred by the ICA -> only differences in
time are meaningful.

Sign convention (verified here, see SIGN in meta.json): the Slab2 mesh in the repo has downward-pointing
triangle normals, so cutde's +dip-slip is NORMAL sense on every patch.  We report ds_thrust_mm = -ds and
along-rake slip with thrust = positive.

Outputs (results/slip_compare/):
  slip_cache.npz                    G, L, Ls, Ld, amplitudes, smoothed rates, fixed dates
  slip_svi_family_meta.csv          per family: patch, depths, rake, n patches in 10/25 km
  slip_svi_family_daily.parquet     fid, date, slip/rate/components for patch, r10, r25
                                    slip_sd_mm = 1-sigma from the vbICA posterior var_V of the time functions
                                    only (slip patterns held fixed) -> a LOWER bound on the daily uncertainty
  slip_wa_nw_family_meta.csv / slip_wa_nw_family_daily.parquet   same for the 132 WA-NW families
  family_id_mapping.json            explicit fid <-> bkey <-> famday.fam <-> families_v1_1 check
  sanity_repo_reproduction.json     reproduction of the repo's run_real.log numbers + episode Mw
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
import pandas as pd

sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
from sse_strain import load_fault_ascii, load_solution                       # noqa: E402
from sse_strain.geodesy import llh2local_xyz                                  # noqa: E402

REPO = "/home/jovyan/cascadia-sse-strain"
SOLDIR = f"{REPO}/data/2026-08-18"
MESH = f"{REPO}/data/cascadia_mesh.txt"
ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/lfe_extension/results/slip_compare"
FAMS12 = f"{ROOT}/lfe_catalog_repair/results/e5r_families_2012.csv"
CAT = f"{ROOT}/data/catalogue/families_v1_1.csv"
COUNTS_SVI = f"{ROOT}/lfe_extension_wa/results/episode_rates/counts_svi_fid_date.parquet"
FAMDAY = f"{ROOT}/lfe_extension/results/refine/rates/famday.parquet"
PBC = f"{ROOT}/lfe_extension/results/refine/conf/P_by_candidate.parquet"
WA_FAM = f"{ROOT}/lfe_extension_wa/results/wa_prod/setup/families.csv"
COUNTS_WA = f"{ROOT}/lfe_extension_wa/results/episode_rates/counts_wa_nw_fid_date_half.parquet"
RADII = (10.0, 25.0)
MU = 40e9                     # repo's nominal fig09 rigidity


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def fixed_dates(t):
    """Gualandi's timeline is year + (doy - 0.5)/ndays (day centres; 2007.00137 = Jan 1).  The repo's
    decimal_year_to_datetime64 ROUNDS frac*ndays (banker's rounding on .5) -> duplicate / skipped days.
    floor() gives the true calendar (UTC) day."""
    t = np.asarray(t, float).ravel()
    yr = np.floor(t).astype(int)
    st = np.array([f"{y:04d}-01-01" for y in yr], dtype="datetime64[D]")
    en = np.array([f"{y + 1:04d}-01-01" for y in yr], dtype="datetime64[D]")
    nd = (en - st).astype(int)
    return st + np.floor((t - yr) * nd + 1e-6).astype(int).astype("timedelta64[D]")


def build_model():
    fault = load_fault_ascii(MESH)
    sol = load_solution(SOLDIR, fault)
    log(f"fault {fault.n_patch} patches; {sol.timeline.size} epochs; ICs {sol.component_index + 1}")
    G = sol.displacement_gf()
    sol.invert(G=G, verbose=False)
    mis = sol.misfit_check(G=G)
    log(f"misfit {np.round(mis, 4)}  sigma0 source: {sol.sigma0_source}")
    Ls, Ld = sol.component_slip_patterns()           # m per unit amplitude
    A = sol.component_amplitudes()                   # (n_t, K)
    rate_amp, t_rate = sol.amplitude_rates(window=7)  # per yr, centred, n_t-6 rows (repo fig07/fig08 rate)
    return fault, sol, G, mis, Ls, Ld, A, rate_amp


def sign_check(fault):
    """Surface response directly above each patch to unit +dip-slip: subsidence + down-dip motion = normal."""
    import cutde.halfspace as HS
    V = fault.vertices_m()
    n = np.cross(V[:, 1] - V[:, 0], V[:, 2] - V[:, 0])
    c = V.mean(1)
    up = np.empty(len(c))
    for i in range(len(c)):
        M = HS.disp_matrix(np.array([[c[i, 0], c[i, 1], 0.0]]), V[i:i + 1], 0.25)
        up[i] = M[0, 2, 0, 1]
    return dict(normals_pointing_down_frac=float((n[:, 2] < 0).mean()),
                plus_ds_gives_subsidence_frac=float((up < 0).mean()),
                convention="cutde +dip-slip on this mesh = NORMAL sense; ds_thrust = -ds")


def locate_on_mesh(fault, lat, lon):
    """-> containing patch (or nearest if outside), barycentric mesh depth, interface point xyz (km)."""
    x, y, _ = llh2local_xyz(np.atleast_1d(lat), np.atleast_1d(lon), np.zeros(1), fault.origin)
    x, y = float(x[0]), float(y[0])
    X, Y, Z = fault.xE, fault.yN, fault.zV
    d = (X[:, 1] - X[:, 0]) * (Y[:, 2] - Y[:, 0]) - (X[:, 2] - X[:, 0]) * (Y[:, 1] - Y[:, 0])
    l1 = ((X[:, 1] - x) * (Y[:, 2] - y) - (X[:, 2] - x) * (Y[:, 1] - y)) / d
    l2 = ((X[:, 2] - x) * (Y[:, 0] - y) - (X[:, 0] - x) * (Y[:, 2] - y)) / d
    l3 = 1 - l1 - l2
    inside = (l1 >= -1e-9) & (l2 >= -1e-9) & (l3 >= -1e-9)
    if inside.any():
        p = int(np.flatnonzero(inside)[0]); ok = True
    else:
        cx, cy = X.mean(1), Y.mean(1)
        p = int(np.argmin((cx - x) ** 2 + (cy - y) ** 2)); ok = False
    z = l1[p] * Z[p, 0] + l2[p] * Z[p, 1] + l3[p] * Z[p, 2]
    return p, ok, float(-z), np.array([x, y, z])


def family_weights(fault, fams):
    """fams: DataFrame fid, lat, lon -> meta DataFrame, weights dict name -> (n_fam, n_patch)."""
    cxyz = fault.centroid_xyz_km
    rows, W = [], {"patch": [], "r10": [], "r25": []}
    for r in fams.itertuples():
        p, ok, zmesh, P = locate_on_mesh(fault, r.lat, r.lon)
        dist = np.sqrt(((cxyz - P) ** 2).sum(1))
        w0 = np.zeros(fault.n_patch); w0[p] = 1.0
        W["patch"].append(w0)
        rec = dict(fid=int(r.fid), lat=r.lat, lon=r.lon, patch=p, inside_mesh=ok, mesh_depth_at_epi_km=round(zmesh, 2),
                   patch_centroid_depth_km=round(float(-fault.centroid_llh[p, 2]), 2),
                   dist_to_patch_centroid_km=round(float(dist[p]), 2))
        for R in RADII:
            m = dist <= R
            w = np.where(m, fault.area * np.exp(-0.5 * (dist / (R / 2.0)) ** 2), 0.0)
            if w.sum() == 0:
                w = w0.copy()
            W[f"r{int(R)}"].append(w / w.sum())
            rec[f"n_patch_r{int(R)}"] = int((w > 0).sum())
        rows.append(rec)
    return pd.DataFrame(rows), {k: np.array(v) for k, v in W.items()}


def family_series(Wd, Ls, Ld, A, rate_amp, n_t):
    """-> dict variant -> (ss, ds_thrust, rate_ss, rate_dsthr) arrays (n_fam, n_t), metres / m per yr."""
    out = {}
    half = 3
    for k, W in Wd.items():
        ps, pd_ = W @ Ls, W @ Ld                      # (n_fam, K)
        ss, ds = ps @ A.T, pd_ @ A.T
        rs = np.full((W.shape[0], n_t), np.nan); rd = rs.copy()
        rs[:, half:half + rate_amp.shape[0]] = ps @ rate_amp.T
        rd[:, half:half + rate_amp.shape[0]] = pd_ @ rate_amp.T
        out[k] = (ss, -ds, rs, -rd, ps, -pd_)
    return out


def rake_from_series(ss, dt):
    """Principal slip direction of the (ss, ds_thrust) cloud; sign so the thrust component is positive."""
    X = np.column_stack([ss - ss.mean(), dt - dt.mean()])
    w, v = np.linalg.eigh(X.T @ X)
    u = v[:, -1]
    if u[1] < 0:
        u = -u
    return float(np.degrees(np.arctan2(u[1], u[0]))), float(w[-1] / w.sum())


def assemble(meta, S, dates, t_days, region, varA):
    """Long table fid x date."""
    n_f, n_t = S["patch"][0].shape
    rakes, frac = [], []
    for i in range(n_f):
        r, f = rake_from_series(S["r25"][0][i], S["r25"][1][i])
        rakes.append(r); frac.append(f)
    meta["rake_deg"] = np.round(rakes, 1)
    meta["rake_var_frac"] = np.round(frac, 3)
    th = np.deg2rad(np.array(rakes))[:, None]
    cols = {"fid": np.repeat(meta.fid.to_numpy(), n_t), "date": np.tile(dates.astype(str), n_f)}
    for k, (ss, dt, rs, rd, ps, pt) in S.items():
        along = ss * np.cos(th) + dt * np.sin(th)
        pal = ps * np.cos(th) + pt * np.sin(th)                    # (n_fam, K) along-rake pattern, m / amplitude
        sd = np.sqrt((pal ** 2) @ varA.T)                          # 1-sigma from var_V only (patterns fixed)
        ralong = rs * np.cos(th) + rd * np.sin(th)
        raw = np.gradient(along, t_days, axis=1)
        sfx = "" if k == "patch" else f"_{k}"
        cols[f"slip_mm{sfx}"] = (1e3 * along).ravel().astype("float32")
        cols[f"slip_sd_mm{sfx}"] = (1e3 * sd).ravel().astype("float32")
        cols[f"slip_rate_mm_per_day{sfx}"] = (1e3 * ralong / 365.25).ravel().astype("float32")
        cols[f"slip_rate_raw_mm_per_day{sfx}"] = (1e3 * raw).ravel().astype("float32")
        cols[f"ss_mm{sfx}"] = (1e3 * ss).ravel().astype("float32")
        cols[f"ds_thrust_mm{sfx}"] = (1e3 * dt).ravel().astype("float32")
    D = pd.DataFrame(cols)
    D["date"] = pd.to_datetime(D.date)
    log(f"[{region}] {len(D)} rows, {n_f} families x {n_t} days")
    return D


def svi_families():
    F = pd.read_csv(FAMS12)
    fid_counts = set(pd.read_parquet(COUNTS_SVI, columns=["fid"]).fid.unique())
    fam_famday = set(pd.read_parquet(FAMDAY, columns=["fam"]).fam.unique())
    fam_pbc = set(pd.read_parquet(PBC, columns=["fam"]).fam.unique())
    C = pd.read_csv(CAT, low_memory=False)
    C = C[C.region == "S_Vancouver_Is"].copy()
    C["family_id"] = C.family_id.astype(float).astype(int)
    J = F.merge(C[["family_id", "family_uid", "lat", "lon", "depth_km"]], left_on="bkey", right_on="family_id",
                how="left", suffixes=("", "_cat"))
    rep = dict(
        e5r_families_2012_rows=len(F), counts_svi_fid_n=len(fid_counts), famday_fam_n=len(fam_famday),
        P_by_candidate_fam_n=len(fam_pbc), catalogue_S_Vancouver_Is_rows=len(C),
        fid_eq_famday_fam=fid_counts == fam_famday, fid_eq_P_by_candidate_fam=fid_counts == fam_pbc,
        fid_subset_of_bkey=fid_counts <= set(F.bkey), catalogue_ids_eq_fid=set(C.family_id.astype(int)) == fid_counts,
        bkey_not_in_production=sorted(set(F.bkey) - fid_counts),
        catalogue_latlon_match_bkey=bool(np.allclose(J.dropna(subset=["family_id"]).lat,
                                                     J.dropna(subset=["family_id"]).lat_cat, atol=1e-6) and
                                         np.allclose(J.dropna(subset=["family_id"]).lon,
                                                     J.dropna(subset=["family_id"]).lon_cat, atol=1e-6)),
        catalogue_depth_match_bkey_dep=bool(np.allclose(J.dropna(subset=["depth_km"]).dep,
                                                        J.dropna(subset=["depth_km"]).depth_km, atol=1e-6)),
        rule="fid (counts_svi_fid_date) == fam (famday, P_by_candidate) == bkey (e5r_families_2012.csv) == "
             "family_id of families_v1_1 rows with region S_Vancouver_Is (family_uid 'SVI:<bkey>'); hypocentre = "
             "Bostock 2015 template location carried in e5r_families_2012 (lat, lon, dep), the same rows svi_ref_common."
             "bostock() joins his catalogue to by exact (lat, lon, dep)")
    F = F[F.bkey.isin(fid_counts)].rename(columns={"bkey": "fid", "dep": "bostock_depth_km"})
    return F[["fid", "lat", "lon", "bostock_depth_km"]].reset_index(drop=True), rep


def wa_families():
    W = pd.read_csv(WA_FAM)
    W = W[W.reg == "NW"]
    cf = set(pd.read_parquet(COUNTS_WA, columns=["fid"]).fid.unique())
    rep = dict(wa_nw_families=len(W), counts_wa_fid_n=len(cf), fid_eq=set(W.fid.astype(int)) == cf)
    W = W.rename(columns={"lat_use": "lat", "lon_use": "lon", "depth_km": "cat_depth_km"})
    return W[["fid", "family_uid", "lat", "lon", "cat_depth_km", "located"]].reset_index(drop=True), rep


def sanity(fault, sol, mis, Ls, Ld, A, rate_amp, dates, G):
    """Reproduce run_real.log (03_slip.py printouts) and the moment of known northern episodes."""
    ss, ds = Ls @ A.T, Ld @ A.T
    mag = np.hypot(ss, ds)
    rate_mag = np.hypot(Ls @ rate_amp.T, Ld @ rate_amp.T)
    mom = MU * (rate_mag * fault.area[:, None] * 1e6).sum(0)
    rep = dict(
        repo_log=dict(misfit_first5=[0.9393, 0.5195, 0.672, 0.5631, 0.9064], peak_slip_cm=3.1,
                      peak_rate_cm_per_yr=108.7, peak_moment_rate=2.338e21),
        mine=dict(misfit_first5=[round(float(x), 4) for x in mis[:5]], peak_slip_cm=round(float(mag.max() * 100), 2),
                  peak_rate_cm_per_yr=round(float(rate_mag.max() * 100), 1), peak_moment_rate=float(f"{mom.max():.4g}")))
    # episode moments: thrust slip increment (mean of 5 days after end - 5 days before start), lat > 47 N
    lat = fault.centroid_llh[:, 1]
    north = lat > 47.0
    eps = [("2010-08-08", "2010-09-05"), ("2011-08-05", "2011-09-10"), ("2012-08-25", "2012-10-10"),
           ("2013-09-01", "2013-10-20"), ("2014-11-14", "2014-12-14"), ("2016-01-25", "2016-02-25"),
           ("2017-01-25", "2017-03-25")]
    dser = pd.Series(np.arange(len(dates)), index=pd.DatetimeIndex(dates))
    out = []
    for a, b in eps:
        i0 = dser[(dser.index >= pd.Timestamp(a) - pd.Timedelta(days=5)) & (dser.index < pd.Timestamp(a))].to_numpy()
        i1 = dser[(dser.index > pd.Timestamp(b)) & (dser.index <= pd.Timestamp(b) + pd.Timedelta(days=5))].to_numpy()
        dss = ss[:, i1].mean(1) - ss[:, i0].mean(1)
        dds = -(ds[:, i1].mean(1) - ds[:, i0].mean(1))
        thr = np.clip(dds, 0, None)
        M0 = MU * (np.hypot(dss, dds) * (dds > 0) * fault.area * 1e6)[north].sum()
        M0t = MU * (thr * fault.area * 1e6)[north].sum()
        out.append(dict(window=f"{a}..{b}", Mw_hypot_thrust_patches=round(2 / 3 * (np.log10(M0) - 9.1), 2),
                        Mw_thrust_component=round(2 / 3 * (np.log10(M0t) - 9.1), 2),
                        peak_thrust_slip_mm=round(float(thr[north].max() * 1e3), 1)))
    rep["episodes_lat_gt_47"] = out
    rep["README_claim"] = "five largest northern episodes Mw 6.6-6.8 with 1-2 cm of slip over 2-5 weeks"
    return rep


def main():
    os.makedirs(OUT, exist_ok=True)
    fault, sol, G, mis, Ls, Ld, A, rate_amp = build_model()
    dates = fixed_dates(sol.timeline)
    rd = sol.dates
    idx = sol.component_index
    varA = np.asarray(sol.ICA["var_V"], float)[:, idx] * sol.S[idx][None, :] ** 2   # (n_t, K) amplitude variance
    assert np.all(np.diff(dates).astype(int) >= 1)
    t_days = (dates - dates[0]).astype(float)
    sg = sign_check(fault)
    log(f"sign: {sg}")
    np.savez_compressed(f"{OUT}/slip_cache.npz", G=G, L=sol.L, Ls=Ls, Ld=Ld, A=A, rate_amp=rate_amp,
                        dates=dates.astype(str), component_index=sol.component_index, mis=mis,
                        station_llh=sol.station_llh,
                        station_names=np.asarray(sol.X["name"]).ravel()[::3].astype(str))
    rep_dates = dict(n_epochs=int(dates.size), first=str(dates[0]), last=str(dates[-1]),
                     gaps_days=[str(dates[i]) for i in np.flatnonzero(np.diff(dates).astype(int) > 1)],
                     repo_dates_duplicates=int((np.diff(rd).astype(int) == 0).sum()),
                     repo_dates_minus_true_days=dict(pd.Series((rd - dates).astype(int)).value_counts().items()))
    log(f"dates {rep_dates}")

    F, idrep = svi_families()
    log(f"id mapping {json.dumps({k: v for k, v in idrep.items() if k != 'rule'}, default=str)}")
    meta, Wd = family_weights(fault, F[["fid", "lat", "lon"]])
    meta = meta.merge(F[["fid", "bostock_depth_km"]], on="fid")
    meta["mesh_minus_bostock_km"] = (meta.mesh_depth_at_epi_km - meta.bostock_depth_km).round(2)
    S = family_series(Wd, Ls, Ld, A, rate_amp, dates.size)
    D = assemble(meta, S, dates, t_days, "SVI", varA)
    D.to_parquet(f"{OUT}/slip_svi_family_daily.parquet", index=False)
    meta.to_csv(f"{OUT}/slip_svi_family_meta.csv", index=False)
    np.savez_compressed(f"{OUT}/slip_svi_weights.npz", fid=meta.fid.to_numpy(), **Wd)

    Wf, warep = wa_families()
    wmeta, WWd = family_weights(fault, Wf[["fid", "lat", "lon"]])
    wmeta = wmeta.merge(Wf[["fid", "family_uid", "cat_depth_km", "located"]], on="fid")
    wmeta["mesh_minus_cat_km"] = (wmeta.mesh_depth_at_epi_km - wmeta.cat_depth_km).round(2)
    WS = family_series(WWd, Ls, Ld, A, rate_amp, dates.size)
    WD = assemble(wmeta, WS, dates, t_days, "WA-NW", varA)
    WD.to_parquet(f"{OUT}/slip_wa_nw_family_daily.parquet", index=False)
    wmeta.to_csv(f"{OUT}/slip_wa_nw_family_meta.csv", index=False)

    json.dump(dict(svi=idrep, wa_nw=warep), open(f"{OUT}/family_id_mapping.json", "w"), indent=1, default=str)
    san = sanity(fault, sol, mis, Ls, Ld, A, rate_amp, dates, G)
    san["sign"] = sg
    san["dates"] = rep_dates
    san["per_family_check"] = "slip_mm_patch rebuilt as W@L@A equals sol.slip() row for that patch (tested below)"
    ss_full, ds_full = sol.slip()
    p = int(meta.patch.iloc[0])
    san["per_family_patch_max_abs_diff_m"] = float(np.abs(ss_full[p] - S["patch"][0][0]).max() +
                                                   np.abs(-ds_full[p] - S["patch"][1][0]).max())
    json.dump(san, open(f"{OUT}/sanity_repo_reproduction.json", "w"), indent=1, default=str)
    log(json.dumps(san, default=str)[:3000])


if __name__ == "__main__":
    main()
