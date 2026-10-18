#!/usr/bin/env python
"""q7_strain (2026-10-07 PT) -- borehole-strain stack on deep between-episode LFE bursts.

Rules frozen in svi_science/results/q7_strainmeter/PLAN.md before any strain number.  Needs q7_process.py output in
data/bsm/processed/.  CPU only, a few GB.

    python -P svi_science/scripts/q7_strain.py
"""
from __future__ import annotations
import os, sys, json, glob, time, warnings
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/q7_strainmeter"
FIG = f"{ROOT}/svi_science/figures"
PROC = f"{ROOT}/data/bsm/processed"
EPR = f"{ROOT}/lfe_extension_wa/results/episode_rates"
REPO = "/home/jovyan/cascadia-sse-strain"
T0 = np.datetime64("2009-01-01T00:00", "ns")
NDRAW = 1000
NINJ = 200
MU = 40e9
NU = 0.25
A_CPL, B_CPL = 1.5, 3.0
W_BURST = 48          # hours, primary half-window
W_ETS = 240           # hours
RMAX_KM = 200.0
MW_GRID = [4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0]
RNG = np.random.default_rng(20261007)
if os.environ.get("Q7_DEBUG"):
    NDRAW, NINJ = 60, 20
    OUT, FIG = "/tmp/q7dbg", "/tmp/q7dbg"
    os.makedirs(OUT, exist_ok=True)
LOG = open(f"{OUT}/q7.log", "a")
RES: dict = {}


def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"
    print(s, flush=True); LOG.write(s + "\n"); LOG.flush()


def mw(m0):
    return (2 / 3) * (np.log10(m0) - 9.1)


def m0_of(mw_):
    return 10 ** (1.5 * mw_ + 9.1)


def hav_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.deg2rad(lat1), np.deg2rad(lat2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.deg2rad(lon2 - lon1) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def hour_index(t):
    t = pd.to_datetime(pd.Series(t), format="ISO8601").values.astype("datetime64[ns]")
    return (t - T0) / np.timedelta64(1, "h")


def mad(x):
    x = x[np.isfinite(x)]
    return float(np.median(np.abs(x - np.median(x)))) * 1.4826 if x.size else np.nan


# ============================================================================ stations
def load_stations():
    from obspy import read_inventory
    S = []
    for f in sorted(glob.glob(f"{PROC}/*_hourly.npz")):
        sta = os.path.basename(f).split("_")[0]
        z = np.load(f)
        az = z["az"]
        if (not np.isfinite(az).all()) or np.allclose(az, 0):
            log(f"drop {sta}: gauge azimuths unknown ({az})"); continue
        inv = read_inventory(f"{ROOT}/data/bsm/stationxml/{sta}.xml")
        ch = inv[0][0]
        Y = z["Y"].astype(float)
        if np.isfinite(Y).sum() < 24 * 365 * 4:
            log(f"drop {sta}: < 1 yr of data"); continue
        S.append(dict(sta=sta, lat=ch.latitude, lon=ch.longitude, az=az, Y=Y, thr=z["thr"], years=z["years"]))
    log(f"stations used: {[s['sta'] for s in S]}")
    return S


def gauge_matrix(az):
    """(4, 6) map from strain [exx eyy ezz exy exz eyz] (ENU) to the 4 gauge readings, nominal isotropic coupling."""
    R = np.zeros((4, 6))
    for i, a in enumerate(np.deg2rad(az)):
        c2, s2 = np.cos(2 * a), np.sin(2 * a)
        # e = 1/2 [A (exx+eyy) + B (-(exx-eyy) cos2phi + 2 exy sin2phi)]
        R[i, 0] = 0.5 * (A_CPL - B_CPL * c2)
        R[i, 1] = 0.5 * (A_CPL + B_CPL * c2)
        R[i, 3] = 0.5 * (2 * B_CPL * s2)
    return R


# ============================================================================ Green's functions
_F: dict = {}


def fault_setup(S):
    sys.path.insert(0, f"{REPO}/src"); sys.path.insert(0, f"{ROOT}/lfe_extension/scripts")
    from sse_strain import load_fault_ascii
    from sse_strain.geodesy import llh2localxy
    import slip_build as SB
    import cutde.halfspace as HS
    fault = load_fault_ascii(SB.MESH)
    tris = fault.vertices_m()
    lat = np.array([s["lat"] for s in S]); lon = np.array([s["lon"] for s in S])
    x, y = llh2localxy(lat, lon, fault.origin)
    obs = np.column_stack([x * 1e3, y * 1e3, np.full(lat.size, -200.0)])
    Sm = HS.strain_matrix(np.ascontiguousarray(obs), np.ascontiguousarray(tris), NU)   # (nsta, 6, np, 3)
    sc = SB.sign_check(fault)
    log(f"sign check: {sc}")
    Rg = np.stack([gauge_matrix(s["az"]) for s in S])                                     # (nsta, 4, 6)
    # gauge response per patch per unit slip: strike (0), thrust (= -dip, 1)
    Gss = np.einsum("sgc,scp->sgp", Rg, Sm[:, :, :, 0])
    Gth = -np.einsum("sgc,scp->sgp", Rg, Sm[:, :, :, 1])
    cxy = fault.centroid_xyz_km
    _F.update(fault=fault, SB=SB, Gss=Gss, Gth=Gth, Rg=Rg, Sm=Sm, cxy=cxy, area=np.asarray(fault.area) * 1e6,
              llh2localxy=llh2localxy, slat=lat, slon=lon)
    # displacement sanity: unit thrust on the patch under B009 -> uplift above it?
    return _F


def template_point(lat, lon, R_km=10.0):
    """(nsta, 4) gauge strain per N.m of uniform thrust on mesh patches within R_km of the interface point."""
    F = _F
    p, ok, zmesh, P = F["SB"].locate_on_mesh(F["fault"], lat, lon)
    d = np.sqrt(((F["cxy"] - P[None, :]) ** 2).sum(1))
    sel = np.flatnonzero(d <= R_km)
    if sel.size == 0:
        sel = np.array([p])
    g = F["Gth"][:, :, sel].sum(-1)                      # 1 m slip on each
    m0 = MU * F["area"][sel].sum()
    return g / m0, dict(mesh_depth_km=zmesh, inside=bool(ok), n_patch=int(sel.size))


def template_tremor(tlat, tlon):
    """Slip moment proportional to tremor epicentre counts per (nearest) patch; per N.m of total moment."""
    F = _F
    x, y = F["llh2localxy"](np.asarray(tlat), np.asarray(tlon), F["fault"].origin)
    from scipy.spatial import cKDTree
    kd = cKDTree(F["cxy"][:, :2])
    _, idx = kd.query(np.column_stack([x, y]))
    cnt = np.bincount(idx, minlength=F["cxy"].shape[0]).astype(float)
    mom = cnt / cnt.sum()                                # N.m per N.m total
    slip = mom / (MU * F["area"])                         # m of slip per patch
    return np.einsum("sgp,p->sg", F["Gth"], slip)


def template_gualandi(dA):
    """Predicted gauge steps (nsta, 4) from Gualandi slip change; dA (K,) amplitude change."""
    F = _F
    z = F["slip"]
    ss = z["Ls"] @ dA; ds = z["Ld"] @ dA                  # cutde convention (ds + = normal)
    eps = np.einsum("scp,p->sc", F["Sm"][:, :, :, 0], ss) + np.einsum("scp,p->sc", F["Sm"][:, :, :, 1], ds)
    return np.einsum("sgc,sc->sg", F["Rg"], eps)


# ============================================================================ step estimator
def steps(Y, a, D, W, return_resid=False):
    """Y (NH,4); a (n,) start hour (int); D (n,) or scalar excluded duration (h, >=1); W half-window (h).
    -> s (n,4) step, ok (n,4).  Model a + b t + s H on [a-W, a) U [a+D, a+D+W)."""
    a = np.asarray(a, int); D = np.broadcast_to(np.asarray(D, int), a.shape)
    NH = Y.shape[0]
    rel = np.r_[np.arange(-W, 0), np.arange(0, W)]
    post = np.r_[np.zeros(W, bool), np.ones(W, bool)]
    idx = a[:, None] + rel[None, :] + np.where(post, 1, 0)[None, :] * D[:, None]
    inb = (idx >= 0) & (idx < NH)
    idxc = np.clip(idx, 0, NH - 1)
    tt = (idx - (a + D / 2)[:, None]) / W                 # scaled time
    s = np.full((a.size, 4), np.nan); ok = np.zeros((a.size, 4), bool)
    res = {} if return_resid else None
    for k in range(4):
        y = Y[idxc, k]
        m = inb & np.isfinite(y)
        cov_pre = m[:, :W].mean(1); cov_post = m[:, W:].mean(1)
        good = (cov_pre >= 0.7) & (cov_post >= 0.7)
        X = np.stack([np.ones_like(tt), tt, np.broadcast_to(post, tt.shape).astype(float)], -1)   # (n, L, 3)
        Xm = X * m[..., None]
        yz = np.where(m, y, 0.0)
        XtX = np.einsum("nli,nlj->nij", Xm, X)
        Xty = np.einsum("nli,nl->ni", Xm, yz)
        g = np.flatnonzero(good)
        if g.size:
            c = np.linalg.solve(XtX[g] + 1e-30 * np.eye(3), Xty[g][..., None])[..., 0]
            s[g, k] = c[:, 2]; ok[g, k] = True
            if return_resid:
                res[k] = (g, c)
    return (s, ok, res) if return_resid else (s, ok)


# ============================================================================ combination
def cov_by_year(S, steps_list, years_list):
    """Robust 4x4 covariance per station-year from control steps.  steps_list[j]: (n, nsta, 4); years_list[j]: (n,)."""
    C = {}
    allS = np.concatenate(steps_list, 0); allY = np.concatenate(years_list, 0)
    for j in range(len(S)):
        X = allS[:, j, :]
        fin = np.isfinite(X).all(1)
        for y in np.unique(allY):
            sel = fin & (allY == y)
            Z = X[sel] if sel.sum() >= 30 else X[fin]
            if Z.shape[0] < 30:
                continue
            Z = Z - np.median(Z, 0)
            m = np.array([mad(Z[:, k]) for k in range(4)])
            Z = np.clip(Z, -5 * m, 5 * m)
            c = np.cov(Z.T)
            c += 1e-3 * np.trace(c) / 4 * np.eye(4)
            c += (2e-9) ** 2 * np.eye(4)            # floor: 2 nanostrain
            C[(j, int(y))] = c
    return C


_CIN: dict = {}


def _cin(C, j, y, code):
    key = (id(C), j, y, code)
    if key not in _CIN:
        k = np.array([(code >> b) & 1 for b in range(4)], bool)
        _CIN[key] = (k, np.linalg.inv(C[(j, y)][np.ix_(k, k)]))
    return _CIN[key]


def gls(s, ok, g, yrs, C, use, geff=None):
    """s (n, nsta, 4) steps; ok same; g (nsta, 4) template per N.m; yrs (n,) year; use (nsta,) bool.
    -> m (n,), sig (n,) in N.m.  With geff (n, nsta, 4): returns the recovery factor sum g'C^-1 geff / sum g'C^-1 g."""
    n = s.shape[0] if s is not None else ok.shape[0]
    num = np.zeros(n); den = np.zeros(n)
    yrs = np.asarray(yrs).astype(int)
    for j in np.flatnonzero(use):
        code = (ok[:, j, :] * np.array([1, 2, 4, 8])).sum(1)
        for y in np.unique(yrs):
            if (j, int(y)) not in C:
                continue
            for cd in np.unique(code[yrs == y]):
                if bin(int(cd)).count("1") < 2:
                    continue
                rows = np.flatnonzero((yrs == y) & (code == cd))
                k, Cin = _cin(C, j, int(y), int(cd))
                w = Cin @ g[j, k]
                if geff is None:
                    num[rows] += s[rows][:, j, k] @ w
                else:
                    num[rows] += geff[rows][:, j, k] @ w
                den[rows] += g[j, k] @ w
    with np.errstate(invalid="ignore", divide="ignore"):
        m = num / den; sig = 1 / np.sqrt(den)
    m[den <= 0] = np.nan; sig[den <= 0] = np.nan
    return (m, sig) if geff is None else m


def gls_eff(ok, g, geff, yrs, C, use):
    return gls(None, ok, g, yrs, C, use, geff=geff)


# ============================================================================ events
def load_episodes():
    E = []
    for f, reg in [("episodes_svi_final.csv", "SVI"), ("episodes_nw_final.csv", "WA")]:
        e = pd.read_csv(f"{EPR}/{f}", parse_dates=["start", "end"]); e["region"] = reg; E.append(e)
    E = pd.concat(E, ignore_index=True)
    E["h0"] = hour_index(E.start); E["h1"] = hour_index(E.end) + 24
    return E


def forbidden_mask(NH, spans_h, pad_h):
    m = np.zeros(NH + 1, int)
    for a, b in spans_h:
        i0 = max(int(np.floor(a - pad_h)), 0); i1 = min(int(np.ceil(b + pad_h)), NH)
        if i1 > i0:
            m[i0] += 1; m[i1] -= 1
    return np.cumsum(m)[:NH] > 0


def load_bursts(E):
    B = pd.read_csv(f"{ROOT}/svi_science/results/q5_bursts/bursts.csv")
    B["kind"] = "BE"
    B["h0"] = hour_index(B.start_utc)
    B["h1"] = B.h0 + B.duration_h
    B["set"] = B["set"].map({"SVI primary": "SVI BE", "WA primary": "WA BE", "WA E1/E2 secondary": "WA BE 2009-14"})
    B["trem"] = B.n_trem_100km_pm6h > 0
    out = [B[["set", "kind", "region", "start_PT", "h0", "h1", "duration_h", "n", "n_fam", "lat", "lon", "dep", "trem"]]]
    # broad bursts (between episodes, local tremor allowed)
    F1 = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv").set_index("fid")
    F2 = pd.read_csv(f"{ROOT}/svi_science/results/a_wa_between/families_wa.csv").set_index("fid")
    for reg, f, F, col in [("SVI", f"{ROOT}/svi_science/results/s7/between_broad_events.parquet", F1, "fam"),
                           ("WA", f"{ROOT}/svi_science/results/a_wa_between/between_broad_events.parquet", F2, "fid")]:
        ev = pd.read_parquet(f); ev = ev[ev.series == "data"].copy()
        ev["fid"] = ev[col]
        ev["lat"] = F.loc[ev.fid, "lat"].values; ev["lon"] = F.loc[ev.fid, "lon"].values; ev["dep"] = F.loc[ev.fid, "dep"].values
        ev = ev.sort_values("time_utc").reset_index(drop=True)
        h = hour_index(ev.time_utc)
        ev["h"] = h
        ev["b"] = np.r_[0, np.cumsum(np.diff(h) >= 24)]
        g = ev.groupby("b")
        bb = pd.DataFrame(dict(set=f"{reg} broad", kind="broad", region=reg, h0=g.h.min(), h1=g.h.max(),
                               n=g.size(), n_fam=g.fid.nunique(), lat=g.lat.mean(), lon=g.lon.mean(), dep=g.dep.mean(),
                               trem=np.nan))
        bb["duration_h"] = bb.h1 - bb.h0
        bb["start_PT"] = (T0 + (bb.h0.values * 3.6e12).astype("timedelta64[ns]")).astype("datetime64[ns]")
        bb["start_PT"] = pd.to_datetime(bb.start_PT).dt.tz_localize("UTC").dt.tz_convert("America/Los_Angeles").dt.strftime("%Y-%m-%d %H:%M")
        out.append(bb.reset_index(drop=True))
    B = pd.concat(out, ignore_index=True)
    B["a"] = np.floor(B.h0).astype(int)
    B["D"] = np.maximum(np.ceil(B.h1).astype(int) - B.a, 1)
    return B


# ============================================================================ draws
def draw_controls(ev_a, ev_D, W, allowed_start, doy_of_hour, span, match_doy=True, tod=True, n=NDRAW * 3):
    """Candidate start hours for one event: same time of day, |doy - doy0| <= 30, within span, window allowed."""
    lo, hi = span
    hrs = np.arange(lo, hi)
    if tod:
        hrs = hrs[(hrs % 24) == (ev_a % 24)]
    if match_doy:
        d = np.abs(doy_of_hour[hrs] - doy_of_hour[ev_a]); d = np.minimum(d, 365 - d)
        hrs = hrs[d <= 30]
    hrs = hrs[allowed_start[np.clip(hrs, 0, allowed_start.size - 1)]]
    hrs = hrs[np.abs(hrs - ev_a) > 24 * 10]
    if hrs.size == 0:
        return hrs
    return RNG.choice(hrs, size=n, replace=True)


def window_allowed(forbid, W, D):
    """allowed[a] True if [a-W, a+D+W) has no forbidden hour."""
    c = np.r_[0, np.cumsum(forbid.astype(int))]
    NH = forbid.size
    a = np.arange(NH)
    lo = np.clip(a - W, 0, NH); hi = np.clip(a + D + W, 0, NH)
    return ((c[hi] - c[lo]) == 0) & (a - W >= 0) & (a + D + W < NH)


# ============================================================================ main analysis of one event set
def run_set(name, ev, S, forbid, span, W, Cstore=None, match_doy=True, tod=True, templates=None, keep_draw_steps=False):
    """ev rows with a, D, template key; returns per-event m, sig, control matrix (n_ev, NDRAW)."""
    NH = S[0]["Y"].shape[0]
    nsta = len(S)
    hours = np.arange(NH)
    doy = ((T0 + (hours * 3.6e12).astype("timedelta64[ns]")).astype("datetime64[D]") -
           (T0 + (hours * 3.6e12).astype("timedelta64[ns]")).astype("datetime64[Y]")).astype(int)
    yr_of_hour = (T0 + (hours * 3.6e12).astype("timedelta64[ns]")).astype("datetime64[Y]").astype(int) + 1970
    allowed_cache = {}
    rows = []
    for i, r in enumerate(ev.itertuples()):
        D = int(r.D)
        if D not in allowed_cache:
            allowed_cache[D] = window_allowed(forbid, W, D)
        cand = draw_controls(int(r.a), D, W, allowed_cache[D], doy, span, match_doy, tod)
        a_all = np.r_[int(r.a), cand]
        st = np.full((a_all.size, nsta, 4), np.nan); ok = np.zeros((a_all.size, nsta, 4), bool)
        for j, s in enumerate(S):
            st[:, j], ok[:, j] = steps(s["Y"], a_all, D, W)
        rows.append(dict(a_all=a_all, st=st, ok=ok, yrs=yr_of_hour[np.clip(a_all, 0, NH - 1)]))
    return rows


def combine(name, ev, rows, S, C, templates):
    """GLS per event and draw; sigma matching; returns dict."""
    nsta = len(S)
    m_obs = np.full(len(ev), np.nan); s_obs = np.full(len(ev), np.nan)
    CTRL = np.full((len(ev), NDRAW), np.nan); CSIG = np.full((len(ev), NDRAW), np.nan)
    CIDX = np.full((len(ev), NDRAW), -1, int)
    for i, (r, R) in enumerate(zip(ev.itertuples(), rows)):
        g = templates[i]
        dist = hav_km(r.lat, r.lon, np.array([s["lat"] for s in S]), np.array([s["lon"] for s in S]))
        use = dist <= RMAX_KM
        m, sg = gls(R["st"], R["ok"], g, R["yrs"], C, use)
        m_obs[i], s_obs[i] = m[0], sg[0]
        if not np.isfinite(sg[0]):
            continue
        good = np.flatnonzero(np.isfinite(m[1:]) & (sg[1:] <= 1.5 * sg[0])) + 1
        if good.size < NDRAW // 2:
            log(f"  {name} event {i}: only {good.size} sigma-matched draws")
        if good.size < min(100, NDRAW // 2):        # too few sigma-matched random days -> event not used
            continue
        pick = good[:NDRAW] if good.size >= NDRAW else RNG.choice(good, NDRAW)
        CTRL[i] = m[pick]; CSIG[i] = sg[pick]; CIDX[i] = pick
    return dict(m=m_obs, sig=s_obs, ctrl=CTRL, csig=CSIG, cidx=CIDX)


def stack_stats(m, sig, ctrl, csig, sel=None):
    if sel is None:
        sel = np.ones(m.size, bool)
    sel = sel & np.isfinite(m) & np.isfinite(sig) & np.isfinite(ctrl).all(1)
    if sel.sum() == 0:
        return None
    w = 1 / sig[sel] ** 2
    obs = float((m[sel] * w).sum() / w.sum())
    wc = 1 / csig[sel] ** 2
    cs = (ctrl[sel] * wc).sum(0) / wc.sum(0)
    mu, sd = float(cs.mean()), float(cs.std())
    z = (obs - mu) / sd
    bound = max(obs - mu, 0) + 3 * sd
    # single-event z vs own control
    zi = (m[sel] - ctrl[sel].mean(1)) / ctrl[sel].std(1)
    # DEVIATION (declared in RESULTS): fixed empirical inverse-variance weights, the same for data and every control
    we = 1 / ctrl[sel].std(1) ** 2
    obs_e = float((m[sel] * we).sum() / we.sum())
    cs_e = (ctrl[sel] * we[:, None]).sum(0) / we.sum()
    mu_e, sd_e = float(cs_e.mean()), float(cs_e.std())
    bound_e = max(obs_e - mu_e, 0) + 3 * sd_e
    zbar = float(zi.mean()); zbar_ctrl = ((ctrl[sel] - ctrl[sel].mean(1, keepdims=True)) / ctrl[sel].std(1, keepdims=True)).mean(0)
    return dict(n=int(sel.sum()), stack_Nm=obs, ctrl_mean=mu, ctrl_sd=sd, z=float(z), p_hi=float((cs >= obs).mean()),
                bound_Nm=float(bound), bound_Mw=float(mw(bound)), stack_Mw_if_pos=float(mw(obs)) if obs > 0 else None,
                n_single_gt2=int((zi > 2).sum()), n_single_gt3=int((zi > 3).sum()),
                frac_single_gt2=float((zi > 2).mean()), median_single_sd_Mw=float(mw(np.median(ctrl[sel].std(1)))),
                emp_stack_Nm=obs_e, emp_ctrl_mean=mu_e, emp_ctrl_sd=sd_e, emp_z=float((obs_e - mu_e) / sd_e),
                emp_p_hi=float((cs_e >= obs_e).mean()), emp_bound_Nm=float(bound_e), emp_bound_Mw=float(mw(bound_e)),
                emp_stack_Mw_if_pos=float(mw(obs_e)) if obs_e > 0 else None,
                mean_single_z=zbar, mean_single_z_ctrl_sd=float(zbar_ctrl.std()),
                mean_single_z_z=float(zbar / zbar_ctrl.std())), cs


# ============================================================================
def main():
    log("=== Q7 start ===")
    S = load_stations()
    NH = S[0]["Y"].shape[0]
    fault_setup(S)
    z = np.load(f"{ROOT}/lfe_extension/results/slip_compare/slip_cache.npz")
    _F["slip"] = dict(Ls=z["Ls"], Ld=z["Ld"], A=z["A"], dates=pd.to_datetime(pd.Series(z["dates"]).str[:10]).values)
    E = load_episodes()
    B = load_bursts(E)
    # exclusion vs episodes +/- 2 d: event window [a-W, a+D+W]
    def overlaps(r, reg, W):
        e = E[E.region == reg]
        lo, hi = r.a - W, r.a + r.D + W
        return bool(((e.h0 - 48 < hi) & (e.h1 + 48 > lo)).any())
    B["excluded"] = [overlaps(r, r.region, W_BURST) for r in B.itertuples()]
    B["has_data"] = True
    log(f"bursts per set: {B.groupby('set').size().to_dict()}; excluded (episode overlap): "
        f"{B[B.excluded].groupby('set').size().to_dict()}")
    Bu = B[~B.excluded].reset_index(drop=True)
    # templates
    tmpl, tinfo = [], []
    for r in Bu.itertuples():
        g, inf = template_point(r.lat, r.lon); tmpl.append(g); tinfo.append(inf)
    Bu["mesh_depth_km"] = [t["mesh_depth_km"] for t in tinfo]
    # expected size: per-station max gauge strain for Mw 5.2 (ne) and nearest-station distance
    Bu["max_gauge_ne_Mw52"] = [float(np.nanmax(np.abs(g)) * m0_of(5.2) * 1e9) for g in tmpl]
    # forbidden hours per region: all episodes +/- 2 d, all bursts +/- 3 d
    forb = {}
    for reg in ["SVI", "WA"]:
        sp = [(a - 48, b + 48) for a, b in E[E.region == reg][["h0", "h1"]].values]
        sp += [(a - 72, b + 72) for a, b in B[B.region == reg][["h0", "h1"]].values]
        forb[reg] = forbidden_mask(NH, sp, 0)
    spans = {"SVI BE": ("2014-01-15", "2026-05-27"), "SVI broad": ("2014-01-15", "2026-05-27"),
             "WA BE": ("2015-01-15", "2026-05-27"), "WA broad": ("2015-01-15", "2026-05-27"),
             "WA BE 2009-14": ("2009-07-01", "2014-12-31")}
    rows_all = {}
    t_start = time.time()
    for name, g in Bu.groupby("set"):
        sp = tuple(int(hour_index([x])[0]) for x in spans[name])
        idx = g.index.values
        rows_all[name] = (idx, run_set(name, g, S, forb[g.region.iloc[0]], sp, W_BURST))
        log(f"steps done {name}: {len(idx)} events ({time.time() - t_start:.0f} s)")
    # covariance per station-year from ALL burst-control draws
    C = cov_by_year(S, [R["st"][1:] for _, rows in rows_all.values() for R in rows],
                    [R["yrs"][1:] for _, rows in rows_all.values() for R in rows])
    log(f"covariances: {len(C)} station-years")
    # GLS
    M = np.full(len(Bu), np.nan); SG = np.full(len(Bu), np.nan)
    CT = np.full((len(Bu), NDRAW), np.nan); CS = np.full((len(Bu), NDRAW), np.nan); CI = np.full((len(Bu), NDRAW), -1)
    ROWS = [None] * len(Bu)
    for name, (idx, rows) in rows_all.items():
        out = combine(name, Bu.loc[idx], rows, S, C, [tmpl[i] for i in idx])
        M[idx], SG[idx], CT[idx], CS[idx], CI[idx] = out["m"], out["sig"], out["ctrl"], out["csig"], out["cidx"]
        for k, i in enumerate(idx):
            ROWS[i] = rows[k]
    Bu["m_Nm"] = M; Bu["sig_gls_Nm"] = SG
    Bu["ctrl_sd_Nm"] = np.nanstd(CT, 1); Bu["z_single"] = (M - np.nanmean(CT, 1)) / Bu.ctrl_sd_Nm
    Bu["single_3sd_Mw"] = mw(3 * Bu.ctrl_sd_Nm)
    Bu["n_sta_used"] = [int((R["ok"][0].sum(1) >= 2).sum()) if R is not None else 0 for R in ROWS]
    # ---------------- stacks
    sets = {
        "PRIMARY: SVI BE + WA BE": Bu.set.isin(["SVI BE", "WA BE"]).values,
        "SVI BE": (Bu.set == "SVI BE").values, "WA BE": (Bu.set == "WA BE").values,
        "WA BE 2009-14": (Bu.set == "WA BE 2009-14").values,
        "all BE (2009-26)": Bu.kind.eq("BE").values,
        "BE big (>=3 LFEs)": Bu.kind.eq("BE").values & (Bu.n >= 3).values,
        "BE with tremor <=100 km +/-6 h": Bu.kind.eq("BE").values & (Bu.trem == True).values,
        "BE without tremor <=100 km": Bu.kind.eq("BE").values & (Bu.trem == False).values,
        "SVI broad": (Bu.set == "SVI broad").values, "WA broad": (Bu.set == "WA broad").values,
        "broad both": Bu.kind.eq("broad").values,
        "broad >=5 LFEs": Bu.kind.eq("broad").values & (Bu.n >= 5).values,
    }
    ST = {}; CSTK = {}
    for k, sel in sets.items():
        r = stack_stats(M, SG, CT, CS, sel)
        if r is None:
            continue
        ST[k], CSTK[k] = r
        log(f"STACK {k}: n={r[0]['n']} stack={r[0]['stack_Nm']:.2e} ctrl {r[0]['ctrl_mean']:.2e}+-{r[0]['ctrl_sd']:.2e} "
            f"z={r[0]['z']:+.2f} bound Mw {r[0]['bound_Mw']:.2f}  single 3sd Mw median {r[0]['median_single_sd_Mw']:.2f} | "
            f"EMP stack {r[0]['emp_stack_Nm']:.2e} z={r[0]['emp_z']:+.2f} bound Mw {r[0]['emp_bound_Mw']:.2f} | mean single z "
            f"{r[0]['mean_single_z']:+.2f} (z {r[0]['mean_single_z_z']:+.2f})")
    RES["stacks"] = ST
    # ---------------- injection (primary set)
    INJ = {}
    for key in ["PRIMARY: SVI BE + WA BE", "all BE (2009-26)", "broad both"]:
        sel = np.flatnonzero(sets[key] & np.isfinite(M) & np.isfinite(CT).all(1))
        st0 = ST[key]
        cur = []
        for mwv in MW_GRID:
            m0 = m0_of(mwv)
            # recovery factor per (event, draw) incl. de-step removal of large hourly increments
            fac = {}
            for i in sel:
                R = ROWS[i]; g = tmpl[i]
                cid = CI[i]
                okd = R["ok"][cid]; yrs = R["yrs"][cid]
                D = int(Bu.D[i])
                inc = np.abs(g) * m0 / max(D, 1)                      # per-hour increment (nsta,4)
                thr = np.stack([[s["thr"][np.clip(y - 2009, 0, len(s["years"]) - 1)] for y in yrs] for s in S], 1)  # (ndraw,nsta,4)
                hit = inc[None] > np.nan_to_num(thr, nan=np.inf)
                if not hit.any():
                    fac[i] = np.ones(cid.size); continue
                geff = np.where(hit, 0.0, g[None])
                dist = hav_km(Bu.lat[i], Bu.lon[i], _F["slat"], _F["slon"]); use = dist <= RMAX_KM
                fac[i] = gls_eff(okd, g, geff, yrs, C, use)
            rec_stack, rec_emp = [], []
            we = 1 / np.nanstd(CT[sel], 1) ** 2
            for rep in range(NINJ):
                pk = RNG.integers(0, NDRAW, size=sel.size)
                mi = np.array([CT[i, p] + m0 * fac[i][p] for i, p in zip(sel, pk)])
                si = np.array([CS[i, p] for i, p in zip(sel, pk)])
                w = 1 / si ** 2
                stk = (mi * w).sum() / w.sum()
                rec_stack.append((stk - st0["ctrl_mean"]) / st0["ctrl_sd"])
                rec_emp.append(((mi * we).sum() / we.sum() - st0["emp_ctrl_mean"]) / st0["emp_ctrl_sd"])
            rec_stack = np.array(rec_stack); rec_emp = np.array(rec_emp)
            zi = np.concatenate([(CT[i] + m0 * fac[i] - np.nanmean(CT[i])) / np.nanstd(CT[i]) for i in sel])
            cur.append(dict(Mw=mwv, stack_recovery=float((rec_stack > 3).mean()), stack_z_median=float(np.median(rec_stack)),
                            emp_stack_recovery=float((rec_emp > 3).mean()), emp_stack_z_median=float(np.median(rec_emp)),
                            single_recovery=float((zi > 3).mean()),
                            mean_recovery_factor=float(np.mean([np.nanmean(fac[i]) for i in sel]))))
            log(f"INJ {key} Mw {mwv}: stack recovered {cur[-1]['stack_recovery']:.2f} (median z {cur[-1]['stack_z_median']:.1f}), "
                f"EMP {cur[-1]['emp_stack_recovery']:.2f} (z {cur[-1]['emp_stack_z_median']:.1f}), "
                f"single {cur[-1]['single_recovery']:.2f}, factor {cur[-1]['mean_recovery_factor']:.2f}")
        INJ[key] = cur
    RES["injection"] = INJ
    # descriptive (not pre-registered): inject each burst with its own slow-earthquake-scaling moment
    # M0 = T x 10^(12.5 +/- 0.5) N.m/s (Ide et al. 2007), T = max(duration, 1 h)
    sel = np.flatnonzero(sets["PRIMARY: SVI BE + WA BE"] & np.isfinite(M) & np.isfinite(CT).all(1))
    st0 = ST["PRIMARY: SVI BE + WA BE"]
    T_s = np.maximum(Bu.duration_h.values[sel], 1.0) * 3600
    scal = []
    for dex in (-0.5, 0.0, 0.5):
        m0i = T_s * 10 ** (12.5 + dex)
        we = 1 / np.nanstd(CT[sel], 1) ** 2
        zs, ze = [], []
        for rep in range(NINJ):
            pk = RNG.integers(0, NDRAW, size=sel.size)
            mi = CT[sel, pk] + m0i               # recovery factor ~1 below Mw 5.8 (see INJ)
            si = CS[sel, pk]; w = 1 / si ** 2
            zs.append(((mi * w).sum() / w.sum() - st0["ctrl_mean"]) / st0["ctrl_sd"])
            ze.append(((mi * we).sum() / we.sum() - st0["emp_ctrl_mean"]) / st0["emp_ctrl_sd"])
        scal.append(dict(log10_rate=12.5 + dex, median_Mw=float(np.median(mw(m0i))), mean_M0=float(m0i.mean()),
                         stack_recovery=float((np.array(zs) > 3).mean()), stack_z_median=float(np.median(zs)),
                         emp_recovery=float((np.array(ze) > 3).mean()), emp_z_median=float(np.median(ze))))
        log(f"SCALING injection 10^{12.5 + dex:.1f} N.m/s: median Mw {scal[-1]['median_Mw']:.2f}, stack recovered "
            f"{scal[-1]['stack_recovery']:.2f} (z {scal[-1]['stack_z_median']:.1f}), EMP {scal[-1]['emp_recovery']:.2f}")
    RES["scaling_injection"] = scal
    # descriptive: single big bursts
    big = Bu[(Bu.kind == "BE") & (Bu.n >= 10)]
    RES["big_bursts"] = big[["set", "start_PT", "n", "duration_h", "m_Nm", "ctrl_sd_Nm", "z_single", "single_3sd_Mw",
                             "n_sta_used"]].to_dict("records")
    log("big bursts:\n" + big[["set", "start_PT", "n", "duration_h", "m_Nm", "z_single", "single_3sd_Mw", "n_sta_used"]].to_string())
    # numerical injection check: add ramp to de-stepped-again raw series for a few (event, draw) pairs
    RES["injection_numeric_check"] = injection_numeric_check(S, Bu, ROWS, CI, tmpl, C, sets["PRIMARY: SVI BE + WA BE"])
    # ---------------- secondary windows W = 1 d, 5 d (primary set only, descriptive)
    RES["alt_windows"] = {}
    for Wa in (24, 120):
        Ba = B[B.set.isin(["SVI BE", "WA BE"])].copy()
        Ba["excluded"] = [overlaps(r, r.region, Wa) for r in Ba.itertuples()]
        Ba = Ba[~Ba.excluded].reset_index(drop=True)
        tm_a = [template_point(r.lat, r.lon)[0] for r in Ba.itertuples()]
        rr = {}
        for name, g in Ba.groupby("set"):
            sp = tuple(int(hour_index([x])[0]) for x in spans[name])
            rr[name] = (g.index.values, run_set(name, g, S, forb[g.region.iloc[0]], sp, Wa))
        Ca = cov_by_year(S, [R["st"][1:] for _, rows in rr.values() for R in rows],
                         [R["yrs"][1:] for _, rows in rr.values() for R in rows])
        Ma = np.full(len(Ba), np.nan); Sa = Ma.copy(); CTa = np.full((len(Ba), NDRAW), np.nan); CSa = CTa.copy()
        for name, (idx, rows) in rr.items():
            o = combine(name, Ba.loc[idx], rows, S, Ca, [tm_a[i] for i in idx])
            Ma[idx], Sa[idx], CTa[idx], CSa[idx] = o["m"], o["sig"], o["ctrl"], o["csig"]
        r = stack_stats(Ma, Sa, CTa, CSa)
        RES["alt_windows"][f"W={Wa // 24} d"] = r[0]
        log(f"ALT W={Wa // 24} d PRIMARY: n={r[0]['n']} z={r[0]['z']:+.2f} bound Mw {r[0]['bound_Mw']:.2f} | EMP z={r[0]['emp_z']:+.2f} "
            f"bound Mw {r[0]['emp_bound_Mw']:.2f}")
    # ---------------- positive control: ETS
    pos = positive_control(S, E, B, C, NH)
    RES["positive_control"] = pos["summary"]
    # ---------------- verdict
    prim = ST["PRIMARY: SVI BE + WA BE"]
    inj52 = [c for c in INJ["PRIMARY: SVI BE + WA BE"] if abs(c["Mw"] - 5.2) < 1e-6][0]
    rec80 = [c["Mw"] for c in INJ["PRIMARY: SVI BE + WA BE"] if c["stack_recovery"] >= 0.8]
    verdict = dict(i_stack_z_gt3=prim["z"] > 3, ii_positive_control=pos["summary"]["P1_pass"],
                   iii_inj52_gt80=inj52["stack_recovery"] > 0.8,
                   smallest_Mw_recovered_80pct=min(rec80) if rec80 else None,
                   bound_Mw_nominal=prim["bound_Mw"])
    gain = pos["summary"].get("checkC_gain")
    if gain and gain > 0:
        verdict["bound_Mw_empirical_gain"] = float(mw(prim["bound_Nm"] / gain))
    verdict["all_pass"] = bool(verdict["i_stack_z_gt3"] and verdict["ii_positive_control"] and verdict["iii_inj52_gt80"])
    rec80e = [c["Mw"] for c in INJ["PRIMARY: SVI BE + WA BE"] if c["emp_stack_recovery"] >= 0.8]
    verdict.update(EMP_i_stack_z_gt3=prim["emp_z"] > 3, EMP_ii_positive_control=pos["summary"]["P1_pass_emp"],
                   EMP_iii_inj52_gt80=inj52["emp_stack_recovery"] > 0.8,
                   EMP_smallest_Mw_recovered_80pct=min(rec80e) if rec80e else None, EMP_bound_Mw_nominal=prim["emp_bound_Mw"])
    verdict["EMP_all_pass"] = bool(verdict["EMP_i_stack_z_gt3"] and verdict["EMP_ii_positive_control"] and verdict["EMP_iii_inj52_gt80"])
    gw = pos["summary"].get("checkC_scale_GLS_median_WA")
    if gw and gw > 0:
        verdict["EMP_bound_Mw_WA_ETS_gain"] = float(mw(prim["emp_bound_Nm"] / gw))
    RES["verdict"] = verdict
    log(f"VERDICT {verdict}")
    # ---------------- save tables
    Bu.drop(columns=[]).to_csv(f"{OUT}/bursts_strain.csv", index=False)
    np.savez_compressed(f"{OUT}/controls.npz", ctrl=CT, csig=CS, m=M, sig=SG,
                        **{f"stack_ctrl_{k.split(':')[0].replace(' ', '_')}": v for k, v in CSTK.items()})
    pd.DataFrame([dict(set=k, **{kk: vv for kk, vv in v.items()}) for k, v in ST.items()]).to_csv(f"{OUT}/stacks.csv", index=False)
    pd.DataFrame([dict(set=k, **c) for k, v in INJ.items() for c in v]).to_csv(f"{OUT}/injection.csv", index=False)
    json.dump(RES, open(f"{OUT}/q7_summary.json", "w"), indent=1, default=float)
    figures(S, Bu, ROWS, CT, CS, CI, tmpl, C, ST, CSTK, INJ, pos, sets)
    log("=== Q7 done ===")


def injection_numeric_check(S, Bu, ROWS, CI, tmpl, C, sel_mask, n=12, mwv=5.6):
    """Add the ramp to the NOT-de-stepped series, re-apply the de-step rule inside the window, refit; compare to
    the analytic (m_ctrl + M0 * factor)."""
    out = []
    sel = np.flatnonzero(sel_mask & np.isfinite(Bu.m_Nm.values))
    m0 = m0_of(mwv)
    for i in RNG.choice(sel, size=min(n, sel.size), replace=False):
        R = ROWS[i]; g = tmpl[i]; D = int(Bu.D[i])
        p = int(CI[i, 0]); a = int(R["a_all"][p]); yr = int(R["yrs"][p])
        st_inj = np.full((1, len(S), 4), np.nan); ok_inj = np.zeros((1, len(S), 4), bool)
        for j, s in enumerate(S):
            raw = np.load(f"{PROC}/{s['sta']}_hourly.npz")["Yraw"].astype(float)
            lo, hi = a - W_BURST - 1, a + D + W_BURST + 1
            seg = raw[lo:hi].copy()
            tt = np.arange(lo, hi)
            ramp = np.clip((tt - a) / max(D, 1), 0, 1)[:, None] * (g[j] * m0)[None, :]
            seg = seg + ramp
            # de-step inside window
            thr = s["thr"][np.clip(yr - 2009, 0, len(s["years"]) - 1)]
            Yw = np.full((raw.shape[0], 4), np.nan)
            for k in range(4):
                x = seg[:, k]; v = np.flatnonzero(np.isfinite(x))
                if v.size < 2:
                    continue
                d = np.diff(x[v]); d[np.abs(d) > thr[k]] = 0.0
                y = np.r_[x[v[0]], x[v[0]] + np.cumsum(d)]
                Yw[lo + v, k] = y
            s1, o1 = steps(Yw, np.array([a]), D, W_BURST)
            st_inj[0, j], ok_inj[0, j] = s1[0], o1[0]
        dist = hav_km(Bu.lat[i], Bu.lon[i], _F["slat"], _F["slon"]); use = dist <= RMAX_KM
        m_num, _ = gls(st_inj, ok_inj, g, np.array([yr]), C, use)
        m_ctrl, _ = gls(R["st"][[p]], R["ok"][[p]], g, np.array([yr]), C, use)
        out.append(dict(event=int(i), injected_Nm=m0, recovered_minus_ctrl_Nm=float(m_num[0] - m_ctrl[0]),
                        ratio=float((m_num[0] - m_ctrl[0]) / m0)))
    log(f"numeric injection check (Mw {mwv}): ratios {[round(o['ratio'], 2) for o in out]}")
    return out


def positive_control(S, E, B, C_burst, NH):
    tr = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    tr["h"] = hour_index(tr.time)
    boxes = {}
    for reg, f in [("SVI", f"{ROOT}/svi_science/results/s0/families.csv"),
                   ("WA", f"{ROOT}/svi_science/results/a_wa_between/families_wa.csv")]:
        F = pd.read_csv(f)
        boxes[reg] = (F.lat.min() - 0.3, F.lat.max() + 0.3, F.lon.min() - 0.4, F.lon.max() + 0.4)
    log(f"tremor boxes: {boxes}")
    Ex = E[(E.start >= "2010-01-15")].copy()
    Ex["P"] = np.where(Ex.in_frozen_major.fillna(False).astype(bool) & (Ex.span_days >= 20), "P1 big ETS",
                       np.where(~Ex.major.astype(bool), "P2 minor", "other major"))
    Ex = Ex[Ex.P != "other major"].reset_index(drop=True)
    Ex["a"] = np.floor(Ex.h0).astype(int); Ex["D"] = (Ex.h1 - Ex.h0).astype(int)
    tm, keep = [], []
    for r in Ex.itertuples():
        b = boxes[r.region]
        t = tr[(tr.h >= r.h0) & (tr.h < r.h1) & tr.lat.between(b[0], b[1]) & tr.lon.between(b[2], b[3])]
        if len(t) < 10:
            keep.append(False); tm.append(None); continue
        keep.append(True); tm.append(template_tremor(t.lat.values, t.lon.values))
    Ex["n_trem"] = [0 if t is None else 1 for t in tm]
    Ex = Ex[np.array(keep)].reset_index(drop=True); tm = [t for t in tm if t is not None]
    Ex["lat"] = 0.0; Ex["lon"] = 0.0
    for i, r in enumerate(Ex.itertuples()):
        b = boxes[r.region]; Ex.loc[i, "lat"] = (b[0] + b[1]) / 2; Ex.loc[i, "lon"] = (b[2] + b[3]) / 2
    # forbidden: P1 controls outside majors (both regions) +/- 2 d; P2 controls outside all episodes +/- 2 d
    majors = E[E.major.astype(bool)]
    f1 = forbidden_mask(NH, [(a - 48, b + 48) for a, b in majors[["h0", "h1"]].values], 0)
    f2 = forbidden_mask(NH, [(a - 48, b + 48) for a, b in E[["h0", "h1"]].values], 0)
    span = (int(hour_index(["2010-01-15"])[0]), int(hour_index(["2026-10-01"])[0]))
    out = {}
    summ = {}
    allrows = {}
    for P, forbid in [("P1 big ETS", f1), ("P2 minor", f2)]:
        sel = np.flatnonzero(Ex.P == P)
        ev = Ex.loc[sel]
        rows = run_set(P, ev, S, forbid, span, W_ETS, match_doy=False, tod=True)
        allrows[P] = (sel, rows)
    C = cov_by_year(S, [R["st"][1:] for _, rows in allrows.values() for R in rows],
                    [R["yrs"][1:] for _, rows in allrows.values() for R in rows])
    for P, (sel, rows) in allrows.items():
        ev = Ex.loc[sel]
        # templates: stations within 200 km of box centre -> use all (box centre distances)
        res = combine(P, ev, rows, S, C, [tm[i] for i in sel])
        r = stack_stats(res["m"], res["sig"], res["ctrl"], res["csig"])
        if r is None:
            log(f"{P}: no usable episodes"); continue
        st, cs = r
        Ex.loc[sel, "m_Nm"] = res["m"]; Ex.loc[sel, "z_single"] = (res["m"] - np.nanmean(res["ctrl"], 1)) / np.nanstd(res["ctrl"], 1)
        Ex.loc[sel, "Mw_strain"] = [float(mw(x)) if x > 0 else np.nan for x in res["m"]]
        summ[P] = st
        out[P] = dict(sel=sel, rows=rows, res=res, cs=cs)
        log(f"POS {P}: n={st['n']} stack M0 {st['stack_Nm']:.2e} (Mw {st['stack_Mw_if_pos']}) z={st['z']:+.2f} "
            f"single>2sd {st['n_single_gt2']}/{st['n']} | EMP z={st['emp_z']:+.2f} Mw {st['emp_stack_Mw_if_pos']} | "
            f"mean single z {st['mean_single_z']:+.2f} (z {st['mean_single_z_z']:+.2f})")
        for reg in ["SVI", "WA"]:
            rr = stack_stats(res["m"], res["sig"], res["ctrl"], res["csig"], (ev.region == reg).values)
            if rr is not None:
                summ[f"{P} {reg}"] = rr[0]
                log(f"   {P} {reg}: n={rr[0]['n']} z={rr[0]['z']:+.2f} EMP z={rr[0]['emp_z']:+.2f} "
                    f"Mw {rr[0]['emp_stack_Mw_if_pos']} single>2sd {rr[0]['n_single_gt2']}")
    p1 = summ.get("P1 big ETS")
    summ["P1_pass"] = bool(p1 and p1["z"] > 3 and p1["frac_single_gt2"] >= 0.5)
    summ["P1_pass_emp"] = bool(p1 and p1["emp_z"] > 3 and p1["frac_single_gt2"] >= 0.5)
    # check C: Gualandi predicted vs observed gauge steps on P1
    dates = _F["slip"]["dates"]; A = _F["slip"]["A"]
    obs_l, pred_l, w_l, sta_l, reg_l = [], [], [], [], []
    m_g = []
    if "P1 big ETS" in out:
        sel, rows = out["P1 big ETS"]["sel"], out["P1 big ETS"]["rows"]
        for i, R in zip(sel, rows):
            r = Ex.loc[i]
            d0 = np.datetime64(pd.Timestamp(r.start).date()) - np.timedelta64(1, "D")
            d1 = np.datetime64(pd.Timestamp(r.end).date()) + np.timedelta64(1, "D")
            i0 = np.searchsorted(dates.astype("datetime64[D]"), d0); i1 = np.searchsorted(dates.astype("datetime64[D]"), d1)
            if i1 >= len(dates) or i0 <= 0:
                continue
            gp = template_gualandi(A[i1] - A[i0])                            # (nsta, 4) predicted strain
            yr = int(R["yrs"][0])
            for j in range(len(S)):
                if (j, yr) not in C:
                    continue
                for k in range(4):
                    if R["ok"][0, j, k]:
                        obs_l.append(R["st"][0, j, k]); pred_l.append(gp[j, k]); w_l.append(1 / C[(j, yr)][k, k])
                        sta_l.append(S[j]["sta"]); reg_l.append(r.region)
            mm, _ = gls(R["st"][[0]], R["ok"][[0]], gp / 1.0, np.array([yr]), C, np.ones(len(S), bool))
            m_g.append((r.region, float(mm[0])))
        o, p_, w = np.array(obs_l), np.array(pred_l), np.array(w_l)
        if o.size:
            gain = float((w * o * p_).sum() / (w * p_ * p_).sum())
            corr = float(np.corrcoef(o, p_)[0, 1])
            sign_agree = float(np.mean(np.sign(o) == np.sign(p_)))
            big = np.abs(p_) > 20e-9
            per_sta = pd.DataFrame(dict(sta=sta_l, o=o, p=p_, w=w)).groupby("sta").apply(
                lambda d: pd.Series(dict(n=len(d), gain=(d.w * d.o * d.p).sum() / (d.w * d.p * d.p).sum(),
                                         corr=np.corrcoef(d.o, d.p)[0, 1] if len(d) > 2 else np.nan,
                                         pred_rms_ne=np.sqrt((d.p ** 2).mean()) * 1e9)))
            per_sta.to_csv(f"{OUT}/checkC_per_station.csv")
            summ.update(checkC_gain=gain, checkC_corr=corr, checkC_sign_agree=sign_agree,
                        checkC_sign_agree_pred_gt20ne=float(np.mean(np.sign(o[big]) == np.sign(p_[big]))) if big.any() else None,
                        checkC_n=int(o.size), checkC_scale_GLS_median=float(np.median([x for _, x in m_g])) if m_g else None,
                        checkC_scale_GLS_median_WA=float(np.median([x for rg, x in m_g if rg == "WA"])),
                        checkC_scale_GLS_median_SVI=float(np.median([x for rg, x in m_g if rg == "SVI"])))
            Ro = np.array(reg_l)
            for rg in ["WA", "SVI"]:
                q = Ro == rg
                summ[f"checkC_gain_{rg}"] = float((w[q] * o[q] * p_[q]).sum() / (w[q] * p_[q] ** 2).sum())
                summ[f"checkC_corr_{rg}"] = float(np.corrcoef(o[q], p_[q])[0, 1])
            pd.DataFrame(dict(sta=sta_l, region=reg_l, obs=o, pred=p_, w=w)).to_csv(f"{OUT}/checkC_obs_pred.csv", index=False)
            log(f"check C: gain {gain:.2f}, corr {corr:.2f}, sign agree {sign_agree:.2f} (n {o.size}); "
                f"GLS scale median {summ['checkC_scale_GLS_median']} WA {summ['checkC_scale_GLS_median_WA']:.2f} "
                f"SVI {summ['checkC_scale_GLS_median_SVI']:.2f}; gain WA {summ['checkC_gain_WA']:.2f} SVI {summ['checkC_gain_SVI']:.2f}")
    Ex.drop(columns=[]).to_csv(f"{OUT}/positive_control_episodes.csv", index=False)
    return dict(summary=summ, out=out, Ex=Ex, C=C, tm=tm)


# ============================================================================ figures
def projected_series(S, R, p, g, C, a, D, W, use):
    """Template-projected strain (N.m equivalent) around one window: residual after the fitted pre/post line
    without the step; returns rel hours (excluded interval collapsed to 0) and series."""
    rel = np.r_[np.arange(-W, 0), np.arange(0, W)]
    num = np.zeros(rel.size); den = 0.0
    yr = int(R["yrs"][p])
    for j, s in enumerate(S):
        if not use[j] or (j, yr) not in C or R["ok"][p, j].sum() < 4:
            continue
        idx = np.r_[np.arange(a - W, a), np.arange(a + D, a + D + W)]
        y = s["Y"][idx]
        if not np.isfinite(y).all():
            continue
        tt = (idx - (a + D / 2)) / W
        Ci = np.linalg.inv(C[(j, yr)]); w = Ci @ g[j]
        pre = rel < 0
        # remove line fitted with step (so the step remains)
        X = np.column_stack([np.ones(rel.size), tt, (~pre).astype(float)])
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        yy = y - X[:, :2] @ c[:2]
        yy = yy - yy[pre].mean(0)
        num += yy @ w; den += g[j] @ w
    return rel, (num / den if den > 0 else np.full(rel.size, np.nan))


def figures(S, Bu, ROWS, CT, CS, CI, tmpl, C, ST, CSTK, INJ, pos, sets):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    # --- fig 1: stacked projected strain, bursts vs random days vs ETS
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    sel = np.flatnonzero(sets["PRIMARY: SVI BE + WA BE"] & np.isfinite(Bu.m_Nm.values))
    def stack_series(pick):
        ser, wts = [], []
        for i in sel:
            R = ROWS[i]; p = pick(i)
            a = int(R["a_all"][p]); D = int(Bu.D[i])
            dist = hav_km(Bu.lat[i], Bu.lon[i], _F["slat"], _F["slon"])
            rel, y = projected_series(S, R, p, tmpl[i], C, a, D, W_BURST, dist <= RMAX_KM)
            if np.isfinite(y).all():
                ser.append(y); wts.append(1 / CS[i, 0] ** 2 if np.isfinite(CS[i, 0]) else 0)
        ser = np.array(ser)
        return rel, np.nanmean(ser, 0) if len(ser) else None
    rel, obs = stack_series(lambda i: 0)
    ctrl = []
    for k in range(40):
        _, c = stack_series(lambda i, k=k: int(CI[i, k]))
        ctrl.append(c)
    ctrl = np.array(ctrl)
    m52 = m0_of(5.2)
    ax[0].fill_between(rel, np.percentile(ctrl, 2.5, 0) / m52, np.percentile(ctrl, 97.5, 0) / m52, color="0.8",
                       label="random days (40 stacks), 95 %")
    ax[0].plot(rel, obs / m52, "C3", lw=1.5, label=f"BE bursts (n={sel.size})")
    ax[0].axvline(0, color="k", lw=0.5)
    ax[0].set_xlabel("hours from burst (burst interval collapsed to 0)")
    ax[0].set_ylabel("stacked strain projected on thrust pattern\n(units of one Mw 5.2 slip)")
    ax[0].set_title("Q7 deep bursts: strain stack vs random days"); ax[0].legend(fontsize=8)
    # panel 2: control distribution and observed for sets
    keys = ["PRIMARY: SVI BE + WA BE", "SVI BE", "WA BE", "WA BE 2009-14", "broad both", "BE big (>=3 LFEs)"]
    for n_, k in enumerate(keys):
        if k not in ST:
            continue
        cs = CSTK[k]
        ax[1].violinplot([cs / m52], positions=[n_], showextrema=False)
        ax[1].plot(n_, ST[k]["stack_Nm"] / m52, "o", color="C3")
        ax[1].plot(n_, (ST[k]["ctrl_mean"] + 3 * ST[k]["ctrl_sd"]) / m52, "_", color="k", ms=14)
    ax[1].set_xticks(range(len(keys))); ax[1].set_xticklabels([k.replace("PRIMARY: ", "").replace(" (", "\n(") for k in keys],
                                                               fontsize=7, rotation=20)
    ax[1].axhline(1, color="C0", lw=0.8, ls="--", label="Mw 5.2 per burst")
    ax[1].set_ylabel("stacked moment per burst / M0(Mw 5.2)"); ax[1].legend(fontsize=8)
    ax[1].set_title("red = bursts, grey = random days, bar = control +3σ")
    # panel 3: ETS positive control series
    P = pos["out"].get("P1 big ETS")
    if P is not None:
        Ex = pos["Ex"]
        for i_, R in zip(P["sel"], P["rows"]):
            r = Ex.loc[i_]
            rel3, y3 = projected_series(S, R, 0, pos["tm"][i_], pos["C"], int(r.a), int(r.D), W_ETS, np.ones(len(S), bool))
            if np.isfinite(y3).all():
                ax[2].plot(rel3 / 24, y3 / 1e18, color="C1", alpha=0.35, lw=0.8)
        ax[2].axvline(0, color="k", lw=0.5)
        ax[2].set_xlabel("days from ETS (episode collapsed to 0)")
        ax[2].set_ylabel("projected moment (10^18 N.m)")
        p1 = pos["summary"]["P1 big ETS"]
        ax[2].set_title(f"positive control: big ETS (n={p1['n']}, stack z={p1['z']:.1f})")
    fig.tight_layout(); fig.savefig(f"{FIG}/q7_strain_stack.png", dpi=140); plt.close(fig)
    # --- fig 2: injection recovery
    fig, ax = plt.subplots(figsize=(6, 4.2))
    for k, c in zip(INJ, ["C3", "C0", "C2"]):
        d = pd.DataFrame(INJ[k])
        ax.plot(d.Mw, d.stack_recovery, "o-", color=c, label=f"stack: {k.replace('PRIMARY: ', '')}")
        ax.plot(d.Mw, d.single_recovery, "s--", color=c, alpha=0.6, ms=4, label=f"single burst: {k.replace('PRIMARY: ', '')}")
    ax.axhline(0.8, color="k", lw=0.6, ls=":"); ax.axvline(5.2, color="k", lw=0.6, ls=":")
    ax.set_xlabel("injected Mw per burst (thrust at burst location)"); ax.set_ylabel("fraction recovered (z > 3)")
    ax.set_title("Q7 injection test"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{FIG}/q7_injection.png", dpi=140); plt.close(fig)
    # --- fig 3: map: stations, bursts, single-burst 3σ Mw
    fig, ax = plt.subplots(figsize=(6.5, 6))
    b = Bu[np.isfinite(Bu.single_3sd_Mw)]
    sc = ax.scatter(b.lon, b.lat, c=b.single_3sd_Mw, cmap="viridis", s=12 + 4 * np.sqrt(b.n), vmin=5.0, vmax=6.2)
    plt.colorbar(sc, ax=ax, label="single-burst 3σ detection threshold (Mw)")
    for s in S:
        ax.plot(s["lon"], s["lat"], "^", color="C3", ms=7); ax.text(s["lon"] + 0.03, s["lat"] + 0.03, s["sta"], fontsize=7)
    ax.set_xlim(-125.8, -121.8); ax.set_ylim(46.2, 49.4); ax.set_aspect(1 / np.cos(np.deg2rad(48)))
    ax.set_title("Q7 strainmeters (triangles) and bursts"); ax.set_xlabel("lon"); ax.set_ylabel("lat")
    fig.tight_layout(); fig.savefig(f"{FIG}/q7_map_sensitivity.png", dpi=140); plt.close(fig)


if __name__ == "__main__":
    warnings.simplefilter("ignore")
    main()
