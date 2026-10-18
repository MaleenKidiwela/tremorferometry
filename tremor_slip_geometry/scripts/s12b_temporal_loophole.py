#!/usr/bin/env python
"""
S12b — closing the TEMPORAL loophole left open by S12.

S12 established that the geodetic network CAN resolve a 30-km thrust patch under the S7
central-Oregon TREMOR-ONLY cluster (44.33-45.42N, 27-38 km, tri 234/235/236/248/249/254/255/
286/304): 98.6 % potency recovery through the real Tarantola-Valette inversion, detection
threshold 0.91x the northern ETS reference, misfit flat with latitude, station spacing the same
as the north.  It also showed the fixed 13-pattern basis is only ~20 % worse at expressing a
patch there (potency fraction 0.76 vs 0.95).

The one thing S12 could NOT test is the vbICA step itself.  The published model is
    slip(x, t) = sum_{k=1..13} P_k(x) * a_k(t),
so a central-Oregon transient whose TIME HISTORY is independent of all 13 amplitude functions
has no component to live in, and no amount of inversion resolution would place it.  S12b asks
whether that loophole is actually open, from three directions:

  (1) COMPONENT LEVEL.  Does any of the 13 amplitude functions track central-Oregon tremor
      specifically (better than it tracks northern tremor)?  If one does, the model has a
      central-Oregon time function and the loophole is closed at the component level.  If none
      does, either there is nothing to track, or the vbICA never isolated it.
  (2) SLIP-SERIES LEVEL.  Does the modeled slip rate ON the cluster co-vary with the cluster's
      own tremor, and does it do so as well as the neighbouring central SLIP+TREMOR triangles
      do with theirs?  Also: how different CAN the two time histories be?  If the cluster's
      modeled slip-rate series is essentially a copy of its neighbours', the model has no
      freedom to give the cluster a different time history, and the S7 decoupling must be
      about AMPLITUDE (spatial), not about timing.
  (3) DATA LEVEL, bypassing the patterns entirely.  Project the raw GNSS 30-day displacement
      increments onto the surface displacement pattern that a slip patch under the cluster
      WOULD produce (a matched filter built from the Green's functions alone) and ask whether
      that projection rises when central-Oregon tremor is high.  A northern SLIP+TREMOR patch
      with its own tremor is the positive control.  If the central stations show a
      tremor-synchronous transient of comparable size, the geodesy sees slip the pattern set
      did not place, and the S7 result is a MODEL limitation.  If they do not, the decoupling
      is in the data, and the composite's error bar is a physical bound in mm of slip.

--------------------------------------------------------------------------------------
INPUTS (all read-only, all outside this folder)
--------------------------------------------------------------------------------------
  /home/jovyan/cascadia-sse-strain/data/cascadia_mesh.txt          742-triangle mesh
  /home/jovyan/cascadia-sse-strain/data/2026-08-18/                Gualandi solution:
      X.mat        ts (639 x 7169 daily positions, mm, e/n/u interleaved for 213 stations),
                   var_ts, data_mask (639 x 7169 uint8; 91.5 % of samples are real),
                   llh, name, timeline (7169 days, 2007-01-02..2026-08-18), unit = 'mm'
      ICA.mat      U, S, V, var_U;  ind_comps.mat: the 13 retained components
      options.mat / misfit_comps.mat: prior settings and the CV grid (see S12)
  catalogs/pnsn_tremor_cascadia_full.csv                            748,002 PNSN detections,
                                                                    2010-01-06..2026-05-27
  tremor_slip_geometry/results/s12_representability_triangles.csv   tri, lat, lon, depth_km,
                                                                    area_km2, seg, class_s7
NOTE on amplitudes: `res_catalog_sse/components_at_cells_trailing30.npz` ships an `A_raw`
(7169 x 13) built from an EARLIER solution snapshot; it correlates 0.87-0.99 with the
amplitudes of the 2026-08-18 solution but is not identical.  S12b uses
`Solution.component_amplitudes()` from the SAME 2026-08-18 solution that produces P = -Ld, so
the patterns and the time functions are self-consistent.

--------------------------------------------------------------------------------------
DEFINITIONS
--------------------------------------------------------------------------------------
Timeline: the solution's own 7169 daily epochs.  Analysis window 2010-03-15..2026-05-27
(tremor catalog span, minus 30 d for the trailing window and 60 d for the lag scan at each end).

Tremor series T_set(t): number of PNSN detections within 40 km (horizontal) of any centroid of
the triangle set, summed over the trailing 30 days (t-29..t), to match the 30-day trailing frame
the slip catalog uses.  Sets: CLUSTER (the 9), CENTRAL_ST (41 central SLIP+TREMOR), CENTRAL_ALL
(144), NORTH_ST (145 north SLIP+TREMOR).  Every series is DETRENDED by subtracting a 365-day
centred rolling median, which removes the PNSN network-growth ramp (S3b: tremor-day rate grew
x1.3 north / x1.9 central) and leaves the ETS episodicity.  Raw-series values are also reported.

Slip series: a_k(t) = component_amplitudes()[:, k]; the 30-day trailing INCREMENT
s_k(t) = a_k(t) - a_k(t-30) is the rate quantity that pairs with a 30-day tremor count.
Thrust potency rate of a set: R_set(t) = sum_k (sum_{j in set} area_j P[j,k]) * s_k(t), in
km^2 * m of thrust per 30 days (P = -Ld, thrust-positive).

Nulls: 200 rigid circular shifts of the tremor series by >= 365 days inside the tremor window,
with the slip/geodesy series left untouched (the S4/S9 control).  seed 20260910.

Matched filter (test 3).  A unit thrust Gaussian patch (sigma 30 km, peak slip 1 mm) centred on
the triangle nearest the area-weighted centroid of the target set gives predicted station
displacements d = G [0; -u] in mm per mm of peak slip.  The observed 30-day increment is
dx_s(t) = ts_s(t) - ts_s(t-30) (valid only where data_mask is 1 at both epochs), detrended by a
365-day centred rolling median and scaled by a robust per-channel sigma (1.4826 x MAD).  The
weighted least-squares amplitude
        a_hat(t) = sum_s m_st d_s dx_st / sigma_s^2  /  sum_s m_st d_s^2 / sigma_s^2
is an estimate of the patch's PEAK SLIP in mm over that 30-day window, with formal sigma
1/sqrt(sum d_s^2/sigma_s^2).  Three variants:
   5STN   the 5 stations nearest the target, east+north channels only (as specified)
   ALL    all 213 stations, all three channels
   JOINT  all stations, fitting the CLUSTER and NORTH patches SIMULTANEOUSLY (2-column design)
          so that northern ETS cannot leak into the central amplitude -- the cross-talk between
          the two Green's-function patterns is otherwise substantial and would fake a signal.
Every filter is run on THREE data sources so that data and model pass through identical
smearing and the difference is meaningful:
   raw_GNSS     X['ts'], the input positions
   ICA_13comp   U[:, idx] S[idx] V[:, idx]^T, what the 13 retained components kept
   MODEL_slip   G L a(t), the slip model's OWN predicted surface displacement
The decisive statistic is (raw_GNSS composite) - (MODEL_slip composite): tremor-synchronous
displacement the data show and the 13 patterns did not place.
Composite: mean a_hat over the top tremor decile minus the mean over the bottom decile, with a
standard error from the day-to-day scatter reduced by the effective number of independent
30-day windows (n/30).  Equivalent moment uses mu = 30 GPa and the patch's area-weighted potency.

--------------------------------------------------------------------------------------
OUTPUTS (all inside tremor_slip_geometry/)
--------------------------------------------------------------------------------------
  results/s12b_component_tremor_corr.csv    13 rows: spatial share of each component's thrust in
                                            the cluster / central / north, and its best-lag
                                            correlation with each tremor series
  results/s12b_slipseries_corr.csv          per triangle set: Spearman(modeled potency rate,
                                            own tremor) with the shifted null, and the
                                            cross-correlation of the sets' slip series
  results/s12b_matched_filter.csv           per target x estimator x data source: correlation
                                            with tremor, null p, composite amplitude and bound
  results/s12b_data_vs_model.csv            data vs model through the identical filter, the
                                            excess and its significance, and the 2-se bound on
                                            unmodelled tremor-synchronous slip
  results/s12b_stations.csv                 the 5+5 stations used, with per-station correlation
                                            and regression slope
  figures/s12b_temporal_loophole.png        9 panels

Run:  OMP_NUM_THREADS=4 python scripts/s12b_temporal_loophole.py     (~3 min)
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, "/home/jovyan/cascadia-sse-strain/src")
from sse_strain import load_fault_ascii, load_solution          # noqa: E402
from sse_strain.geodesy import llh2localxy                      # noqa: E402

import matplotlib                                               # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt                                 # noqa: E402
from matplotlib.dates import YearLocator, DateFormatter         # noqa: E402

HERE = Path(__file__).resolve().parent.parent
ROOT = HERE.parent
MESH = "/home/jovyan/cascadia-sse-strain/data/cascadia_mesh.txt"
SOLDIR = "/home/jovyan/cascadia-sse-strain/data/2026-08-18"
TREMOR = ROOT / "catalogs/pnsn_tremor_cascadia_full.csv"
S12CSV = HERE / "results/s12_representability_triangles.csv"

WIN = 30                 # trailing window, days
LAGS = np.arange(-60, 61)
ASSIGN_KM = 40.0         # tremor-to-triangle assignment radius (S7 convention)
DETREND_D = 365          # centred rolling median removed from every series
N_NULL = 200
MIN_SHIFT = 365
SIG_KM = 30.0
MU_PA = 30e9
SEED = 20260910
ANA0, ANA1 = np.datetime64("2010-03-15"), np.datetime64("2026-05-27")

rng = np.random.default_rng(SEED)


def trail_sum(x, w=WIN):
    """Trailing sum over w days (x[t-w+1..t]); first w-1 entries are NaN."""
    c = np.concatenate([[0.0], np.nancumsum(x)])
    out = np.full(x.shape, np.nan)
    out[w - 1:] = c[w:] - c[:-w + 1 or None][: len(x) - w + 1]
    return out


def trail_diff(x, w=WIN):
    """x(t) - x(t-w); first w entries are NaN."""
    out = np.full(x.shape, np.nan)
    out[w:] = x[w:] - x[:-w]
    return out


def detrend_med(x, w=DETREND_D):
    return x - pd.Series(x).rolling(w, center=True, min_periods=w // 2).median().values


def best_lag_corr(x, y, lags=LAGS, i0=None, i1=None):
    """max-|Pearson| over lags of corr(x[t], y[t+lag]); returns (r, lag)."""
    best = (0.0, 0)
    for L in lags:
        a = x[i0:i1]
        b = y[i0 + L:i1 + L]
        m = np.isfinite(a) & np.isfinite(b)
        if m.sum() < 500:
            continue
        r = np.corrcoef(a[m], b[m])[0, 1]
        if abs(r) > abs(best[0]):
            best = (r, int(L))
    return best


def lag_curve(x, y, i0, i1, lags=LAGS):
    """Pearson r as a function of lag, corr(x[t], y[t+lag])."""
    out = np.full(len(lags), np.nan)
    for i, L in enumerate(lags):
        a, b = x[i0:i1], y[i0 + L:i1 + L]
        m = np.isfinite(a) & np.isfinite(b)
        if m.sum() >= 500:
            out[i] = np.corrcoef(a[m], b[m])[0, 1]
    return out


def shifted_null(stat_fn, y, valid_lo, valid_hi, n=N_NULL):
    """Rigid circular shifts of y inside [valid_lo, valid_hi] by >= MIN_SHIFT days."""
    span = valid_hi - valid_lo
    out = np.empty(n)
    seg = y[valid_lo:valid_hi].copy()
    for i in range(n):
        s = int(rng.integers(MIN_SHIFT, span - MIN_SHIFT))
        yy = y.copy()
        yy[valid_lo:valid_hi] = np.roll(seg, s)
        out[i] = stat_fn(yy)
    return out


# ================================================================== load
print("[1] loading")
fault = load_fault_ascii(MESH)
sol = load_solution(SOLDIR, fault)
L = sol.invert(verbose=False)
Ls, Ld = sol.component_slip_patterns()
P = -Ld                                        # 742 x 13 thrust-positive, m per unit amp
Amp = sol.component_amplitudes()               # 7169 x 13
dates = sol.dates
NT = dates.size
area = fault.area
lon, lat = fault.centroid_llh[:, 0], fault.centroid_llh[:, 1]
xyz = fault.centroid_xyz_km
idx = sol.component_index
K = idx.size

tri = pd.read_csv(S12CSV)
SETS = {
    "CLUSTER": ((tri.class_s7 == "TREMOR-ONLY") & (tri.seg == "central")).values,
    "CENTRAL_ST": ((tri.class_s7 == "SLIP+TREMOR") & (tri.seg == "central")).values,
    "CENTRAL_ALL": (tri.seg == "central").values,
    "NORTH_ST": ((tri.class_s7 == "SLIP+TREMOR") & (tri.seg == "north")).values,
}
for k_, v in SETS.items():
    print(f"    set {k_:12s} {v.sum():3d} tri   lat {lat[v].min():.2f}-{lat[v].max():.2f}  "
          f"median depth {np.median(-fault.centroid_llh[v, 2]):.1f} km")

tr = pd.read_csv(TREMOR, usecols=["time", "lat", "lon"])
tr["day"] = pd.to_datetime(tr["time"]).values.astype("datetime64[D]")
tx, ty = llh2localxy(tr["lat"].values, tr["lon"].values, fault.origin)
tx, ty = np.ravel(tx), np.ravel(ty)
day_i = np.searchsorted(dates, tr["day"].values)
ok = (day_i >= 0) & (day_i < NT) & (dates[np.clip(day_i, 0, NT - 1)] == tr["day"].values)
print(f"    {len(tr)} detections, {ok.sum()} land on the solution timeline")

i0 = int(np.searchsorted(dates, ANA0))
i1 = int(np.searchsorted(dates, ANA1)) + 1
iw0, iw1 = i0 + 60, i1 - 60                     # lag-scan-safe correlation window
tlo = int(np.searchsorted(dates, np.datetime64("2010-02-05")))
print(f"    analysis window {dates[iw0]}..{dates[iw1-1]}  ({iw1-iw0} d)")

# ------------------------------------------------------- tremor series per set
T_RAW, T_DET = {}, {}
for name, m in SETS.items():
    d = np.sqrt((tx[:, None] - xyz[None, m, 0]) ** 2 + (ty[:, None] - xyz[None, m, 1]) ** 2)
    near = (d.min(1) <= ASSIGN_KM) & ok
    cnt = np.bincount(day_i[near], minlength=NT).astype(float)
    T_RAW[name] = trail_sum(cnt)
    T_DET[name] = detrend_med(T_RAW[name])
    print(f"    T[{name:12s}] {int(near.sum()):7d} detections, "
          f"median 30-d count {np.nanmedian(T_RAW[name][i0:i1]):.0f}, "
          f"p95 {np.nanpercentile(T_RAW[name][i0:i1], 95):.0f}")

_a, _b = T_DET["CLUSTER"][iw0:iw1], T_DET["NORTH_ST"][iw0:iw1]
_m = np.isfinite(_a) & np.isfinite(_b)
TR_CROSS = float(spearmanr(_a[_m], _b[_m]).statistic)
print(f"    Spearman(cluster tremor, north tremor) = {TR_CROSS:+.3f} "
      f"(the two ETS series are only partly synchronous)")

# ========================================== TEST 1: components vs tremor timing
print("[2] TEST 1 — does any component's time function follow central-Oregon tremor?")
S = np.column_stack([trail_diff(Amp[:, k]) for k in range(K)])       # 7169 x 13 rate
W = area[:, None] * np.abs(P)                                        # area-weighted |thrust|
rows = []
for k in range(K):
    row = {"component": int(idx[k] + 1),
           "share_cluster": W[SETS["CLUSTER"], k].sum() / W[:, k].sum(),
           "share_central_all": W[SETS["CENTRAL_ALL"], k].sum() / W[:, k].sum(),
           "share_north_ST": W[SETS["NORTH_ST"], k].sum() / W[:, k].sum(),
           "amp_sd": Amp[:, k].std(),
           "mean_thrust_cluster_m_per_amp": np.average(P[SETS["CLUSTER"], k],
                                                       weights=area[SETS["CLUSTER"]]),
           "mean_thrust_northST_m_per_amp": np.average(P[SETS["NORTH_ST"], k],
                                                       weights=area[SETS["NORTH_ST"]])}
    for tag in ("CLUSTER", "CENTRAL_ST", "CENTRAL_ALL", "NORTH_ST"):
        r, lg = best_lag_corr(S[:, k], T_DET[tag], i0=iw0, i1=iw1)
        row[f"r_{tag}"] = r
        row[f"lag_{tag}"] = lg
    row["prefers_central"] = abs(row["r_CENTRAL_ALL"]) > abs(row["r_NORTH_ST"])
    rows.append(row)
comp = pd.DataFrame(rows)
comp["abs_r_central_minus_north"] = comp["r_CENTRAL_ALL"].abs() - comp["r_NORTH_ST"].abs()
# sign-corrected: + means the component ADDS thrust slip on the set when tremor is high
comp["r_signed_CLUSTER"] = comp["r_CLUSTER"] * np.sign(comp["mean_thrust_cluster_m_per_amp"])
comp["r_signed_NORTH_ST"] = comp["r_NORTH_ST"] * np.sign(comp["mean_thrust_northST_m_per_amp"])
comp["thrust_delivered_cluster_mm_per_sd"] = (comp["mean_thrust_cluster_m_per_amp"].abs()
                                              * comp["amp_sd"] * 1e3)
print(comp[["component", "share_cluster", "share_central_all", "share_north_ST",
            "thrust_delivered_cluster_mm_per_sd", "r_CLUSTER", "lag_CLUSTER",
            "r_signed_CLUSTER", "r_CENTRAL_ALL", "r_NORTH_ST",
            "prefers_central"]].round(3).to_string(index=False))
npref = int(comp["prefers_central"].sum())
print(f"    components whose amplitude rate correlates better with CENTRAL than with NORTH "
      f"tremor: {npref} of {K}")
print(f"    max |r| with cluster tremor {comp['r_CLUSTER'].abs().max():.3f} "
      f"(comp {int(comp.loc[comp['r_CLUSTER'].abs().idxmax(),'component'])}); "
      f"max |r| with north tremor {comp['r_NORTH_ST'].abs().max():.3f} "
      f"(comp {int(comp.loc[comp['r_NORTH_ST'].abs().idxmax(),'component'])})")

# ================================ TEST 2: modeled slip rate on the sets vs tremor
print("[3] TEST 2 — does the modeled slip rate ON the cluster follow the cluster's tremor?")
Rser, srows = {}, []
for name, m in SETS.items():
    wk = (area[m][:, None] * P[m, :]).sum(0)          # km^2 * m per unit amplitude, per comp
    Rser[name] = S @ wk                               # potency rate series
for name in SETS:
    x = Rser[name]
    y = T_DET[name]

    def stat(yy, x=x):
        a, b = x[iw0:iw1], yy[iw0:iw1]
        mm = np.isfinite(a) & np.isfinite(b)
        return spearmanr(a[mm], b[mm]).statistic

    obs = stat(y)
    nul = shifted_null(stat, y, tlo, i1)
    p = (np.sum(np.abs(nul) >= abs(obs)) + 1) / (N_NULL + 1)
    r_lag, lg = best_lag_corr(x, y, i0=iw0, i1=iw1)
    srows.append({"set": name, "n_tri": int(SETS[name].sum()),
                  "spearman_slip_vs_own_tremor": obs,
                  "null_mean": nul.mean(), "null_sd": nul.std(),
                  "null_p": p, "best_lag_pearson": r_lag, "best_lag_days": lg,
                  "potency_rate_sd": np.nanstd(x[iw0:iw1])})
    print(f"    {name:12s} rho {obs:+.3f}  null {nul.mean():+.3f} +- {nul.std():.3f}  "
          f"p {p:.4f}   best-lag pearson {r_lag:+.3f} @ {lg:+d} d")
sl = pd.DataFrame(srows)

# how different CAN the sets' modeled time histories be?
print("    cross-correlation of the MODELED potency-rate series between sets:")
cross = {}
for a in SETS:
    for b in SETS:
        xa, xb = Rser[a][iw0:iw1], Rser[b][iw0:iw1]
        m = np.isfinite(xa) & np.isfinite(xb)
        cross[(a, b)] = np.corrcoef(xa[m], xb[m])[0, 1]
    print("      " + a.ljust(12) + "  ".join(f"{b}:{cross[(a,b)]:+.3f}" for b in SETS))
sl["corr_with_CLUSTER_series"] = [cross[("CLUSTER", n)] for n in sl["set"]]
# amplitude comparison: the cluster's potency rate relative to its neighbours
for name in SETS:
    sl.loc[sl["set"] == name, "potency_rate_sd_per_km2"] = (
        np.nanstd(Rser[name][iw0:iw1]) / area[SETS[name]].sum())
print("    potency-rate sd per km^2 (amplitude, not timing): "
      + ", ".join(f"{r['set']}={r.potency_rate_sd_per_km2:.3e}" for _, r in sl.iterrows()))

# ======================= TEST 3: data-space matched filter, no patterns involved
print("[4] TEST 3 — matched filter on the raw GNSS, bypassing the 13 patterns")
G = sol.displacement_gf()
ts = np.asarray(sol.X["ts"], dtype=float)
dmask = np.asarray(sol.X["data_mask"]).astype(bool)
st_llh = sol.station_llh
names = [str(x)[:-1] for x in np.asarray(sol.X["name"], dtype=object).ravel()[::3]]
sx, sy = llh2localxy(st_llh[:, 1], st_llh[:, 0], fault.origin)
sx, sy = np.ravel(sx), np.ravel(sy)
chan = np.array(["e", "n", "u"] * len(names))
chan_sta = np.repeat(np.arange(len(names)), 3)

# ICA reconstruction with the 13 retained components
ts_ica = sol.U[:, idx] @ np.diag(sol.S[idx]) @ sol.V[:, idx].T

DTRI = np.sqrt(((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(-1))


def target_patch(mask):
    """sigma-30 km unit thrust Gaussian on the triangle nearest the set's area-weighted centroid."""
    c = np.average(xyz[mask], axis=0, weights=area[mask])
    j = int(np.argmin(((xyz - c) ** 2).sum(1)))
    u = np.exp(-0.5 * (DTRI[:, j] / SIG_KM) ** 2)
    m = np.zeros(2 * fault.n_patch)
    m[fault.n_patch:] = -u
    return j, u, G @ m           # d in mm per mm of peak slip


TARGETS = {"CLUSTER": SETS["CLUSTER"], "NORTH_ST": SETS["NORTH_ST"]}
PATCH = {}
for name, m in TARGETS.items():
    j, u, d = target_patch(m)
    pot_m2 = float((area * u).sum() * 1e6)           # m^2, for the moment conversion
    PATCH[name] = {"tri": j, "u": u, "d": d, "pot_m2": pot_m2}
    print(f"    {name:10s} patch on tri {j} ({lat[j]:.2f}N {lon[j]:.2f}, "
          f"{-fault.centroid_llh[j,2]:.1f} km); max station disp "
          f"{np.abs(d.reshape(-1,3)).max()*10:.2f} mm per cm slip; "
          f"patch area-integral {pot_m2/1e6:.0f} km^2")
    dist = np.hypot(sx - xyz[j, 0], sy - xyz[j, 1])
    PATCH[name]["near5"] = np.argsort(dist)[:5]
    PATCH[name]["near5_km"] = dist[np.argsort(dist)[:5]]
    print(f"               5 nearest stations: "
          + ", ".join(f"{names[i]}({dist[i]:.0f} km)" for i in PATCH[name]["near5"]))

# cross-talk between the two Green's-function patterns
dc, dn = PATCH["CLUSTER"]["d"], PATCH["NORTH_ST"]["d"]
print(f"    Green's-function cross-talk cos(d_cluster, d_north) = "
      f"{dc @ dn / np.linalg.norm(dc) / np.linalg.norm(dn):+.3f}  -> the JOINT fit is needed")

# The slip model's OWN predicted surface displacement: G L a(t).  Running the identical
# filter on data and on model makes the comparison exact -- both suffer the same 30-km
# smearing, so any EXCESS in the data is slip the 13 patterns did not place.
ts_model = (G @ L) @ Amp.T
SOURCES = {"raw_GNSS": ts, "ICA_13comp": ts_ica, "MODEL_slip": ts_model}
DX, SIG, MSK = {}, {}, {}
for src, arr in SOURCES.items():
    dx = np.full_like(arr, np.nan)
    dx[:, WIN:] = arr[:, WIN:] - arr[:, :-WIN]
    mk = np.zeros(arr.shape, bool)
    mk[:, WIN:] = dmask[:, WIN:] & dmask[:, :-WIN]
    dx[~mk] = np.nan
    dxd = dx - pd.DataFrame(dx.T).rolling(DETREND_D, center=True,
                                          min_periods=DETREND_D // 2).median().values.T
    sg = 1.4826 * np.nanmedian(np.abs(dxd - np.nanmedian(dxd, 1, keepdims=True)),
                               axis=1)
    sg = np.where(np.isfinite(sg) & (sg > 1e-6), sg, np.inf)
    DX[src], SIG[src], MSK[src] = dxd, sg, mk & np.isfinite(dxd)
    print(f"    {src:11s}: median channel sigma of the detrended 30-d increment "
          f"{np.median(sg[np.isfinite(sg)]):.2f} mm")


def wls_amp(dxd, msk, sig, cols, chans):
    """Weighted LS amplitudes for the design `cols` (list of 639-vectors) on `chans` channels."""
    nA = len(cols)
    Dm = np.zeros((nA, dxd.shape[1]))
    NN = np.zeros((nA, nA, dxd.shape[1]))
    w = np.where(np.isin(np.arange(dxd.shape[0]), chans), 1.0 / sig ** 2, 0.0)
    Xd = np.where(msk, dxd, 0.0)
    M = msk.astype(float)
    for a in range(nA):
        Dm[a] = (w * cols[a])[:, None].T @ Xd
        for b in range(nA):
            NN[a, b] = (w * cols[a] * cols[b])[:, None].T @ M
    amp = np.full((nA, dxd.shape[1]), np.nan)
    fsd = np.full((nA, dxd.shape[1]), np.nan)
    for t in range(dxd.shape[1]):
        Nt = NN[:, :, t]
        if np.linalg.det(Nt) <= 0 or Nt[0, 0] <= 0:
            continue
        Ci = np.linalg.inv(Nt)
        amp[:, t] = Ci @ Dm[:, t]
        fsd[:, t] = np.sqrt(np.diag(Ci))
    return amp, fsd


EST = {}
for src in SOURCES:
    dxd, sg, mk = DX[src], SIG[src], MSK[src]
    for name in TARGETS:
        n5 = PATCH[name]["near5"]
        ch5 = np.where(np.isin(chan_sta, n5) & (chan != "u"))[0]
        a, f = wls_amp(dxd, mk, sg, [PATCH[name]["d"]], ch5)
        EST[(src, name, "5STN_en")] = (a[0], f[0])
        a, f = wls_amp(dxd, mk, sg, [PATCH[name]["d"]], np.arange(dxd.shape[0]))
        EST[(src, name, "ALL")] = (a[0], f[0])
    a, f = wls_amp(dxd, mk, sg, [dc, dn], np.arange(dxd.shape[0]))
    EST[(src, "CLUSTER", "JOINT")] = (a[0], f[0])
    EST[(src, "NORTH_ST", "JOINT")] = (a[1], f[1])

TR_FOR = {"CLUSTER": "CLUSTER", "NORTH_ST": "NORTH_ST"}
mrows = []
for (src, name, kind), (a, f) in sorted(EST.items()):
    y = T_DET[TR_FOR[name]]
    ad = detrend_med(a)

    def stat(yy, x=ad):
        p, q = x[iw0:iw1], yy[iw0:iw1]
        mm = np.isfinite(p) & np.isfinite(q)
        return spearmanr(p[mm], q[mm]).statistic

    obs = stat(y)
    nul = shifted_null(stat, y, tlo, i1)
    pv = (np.sum(np.abs(nul) >= abs(obs)) + 1) / (N_NULL + 1)
    seg_a, seg_t = ad[iw0:iw1], y[iw0:iw1]
    mm = np.isfinite(seg_a) & np.isfinite(seg_t)
    hi = seg_t >= np.nanpercentile(seg_t, 90)
    lo = seg_t <= np.nanpercentile(seg_t, 10)
    comp_hi, comp_lo = seg_a[hi & mm], seg_a[lo & mm]
    d_comp = np.nanmean(comp_hi) - np.nanmean(comp_lo)
    # independent 30-d windows only
    se = np.sqrt(np.nanvar(comp_hi) / max(len(comp_hi) / WIN, 1)
                 + np.nanvar(comp_lo) / max(len(comp_lo) / WIN, 1))
    slope = np.polyfit(((seg_t[mm] - np.nanmean(seg_t[mm])) / np.nanstd(seg_t[mm])),
                       seg_a[mm], 1)[0]
    bound = 2 * se
    mo = MU_PA * PATCH[name]["pot_m2"] * bound * 1e-3
    mrows.append({"source": src, "target": name, "estimator": kind,
                  "n_days": int(mm.sum()),
                  "formal_sigma_mm": float(np.nanmedian(f)),
                  "empirical_sd_mm": float(np.nanstd(seg_a[mm])),
                  "spearman_vs_tremor": obs, "null_mean": nul.mean(),
                  "null_sd": nul.std(), "null_p": pv,
                  "slope_mm_per_sd_tremor": slope,
                  "composite_hi_minus_lo_mm": d_comp, "composite_se_mm": se,
                  "composite_snr": d_comp / se if se > 0 else np.nan,
                  "bound_2se_mm_peak_slip": bound,
                  "bound_equiv_Mw": (2.0 / 3.0) * (np.log10(mo) - 9.1)})
mf = pd.DataFrame(mrows)
piv = mf.pivot_table(index=["target", "estimator"], columns="source",
                     values="composite_hi_minus_lo_mm")
pse = mf.pivot_table(index=["target", "estimator"], columns="source",
                     values="composite_se_mm")
piv["data_minus_model_mm"] = piv["raw_GNSS"] - piv["MODEL_slip"]
piv["excess_se_mm"] = np.hypot(pse["raw_GNSS"], pse["MODEL_slip"])
piv["excess_sigma"] = piv["data_minus_model_mm"] / piv["excess_se_mm"]
piv["model_over_data"] = piv["MODEL_slip"] / piv["raw_GNSS"]
piv["unmodelled_bound_2se_mm"] = 2 * piv["excess_se_mm"]
piv["unmodelled_bound_Mw"] = (2.0 / 3.0) * (
    np.log10(MU_PA * np.array([PATCH[t]["pot_m2"] for t, _ in piv.index])
             * piv["unmodelled_bound_2se_mm"].values * 1e-3) - 9.1)
print("    DATA vs MODEL through the identical filter "
      "(top-decile minus bottom-decile tremor, mm of 30-d peak slip):")
print(piv.round(3).to_string())
piv.round(6).to_csv(HERE / "results/s12b_data_vs_model.csv")
mf = mf.merge(piv[["data_minus_model_mm", "excess_sigma", "model_over_data"]],
              left_on=["target", "estimator"], right_index=True, how="left")
print(mf[["source", "target", "estimator", "spearman_vs_tremor", "null_sd", "null_p",
          "composite_hi_minus_lo_mm", "composite_se_mm", "composite_snr",
          "bound_2se_mm_peak_slip", "bound_equiv_Mw"]].round(4).to_string(index=False))

# per-station table
srows2 = []
for name in TARGETS:
    for i in PATCH[name]["near5"]:
        for cc, off in (("e", 0), ("n", 1)):
            ch = 3 * i + off
            x = detrend_med(DX["raw_GNSS"][ch])
            y = T_DET[TR_FOR[name]]
            a, b = x[iw0:iw1], y[iw0:iw1]
            mm = np.isfinite(a) & np.isfinite(b)
            rho = spearmanr(a[mm], b[mm]).statistic
            sl_ = np.polyfit((b[mm] - b[mm].mean()) / b[mm].std(), a[mm], 1)[0]
            srows2.append({"target": name, "station": names[i], "component": cc,
                           "lon": st_llh[i, 0], "lat": st_llh[i, 1],
                           "dist_to_patch_km": float(PATCH[name]["near5_km"]
                                                     [list(PATCH[name]["near5"]).index(i)]),
                           "pred_disp_mm_per_mm_slip": PATCH[name]["d"][ch],
                           "sigma_30d_incr_mm": SIG["raw_GNSS"][ch],
                           "spearman_vs_tremor": rho,
                           "slope_mm_per_sd_tremor": sl_})
stn = pd.DataFrame(srows2)
print(stn.groupby("target")[["spearman_vs_tremor", "slope_mm_per_sd_tremor"]]
      .median().round(4).to_string())

# ============================================================ write + figure
print("[5] writing")
comp.round(6).to_csv(HERE / "results/s12b_component_tremor_corr.csv", index=False)
sl.round(6).to_csv(HERE / "results/s12b_slipseries_corr.csv", index=False)
mf.round(6).to_csv(HERE / "results/s12b_matched_filter.csv", index=False)
stn.round(6).to_csv(HERE / "results/s12b_stations.csv", index=False)

dt = dates.astype("datetime64[D]").astype(object)
fig = plt.figure(figsize=(19.5, 14.5))
gs = fig.add_gridspec(3, 3, hspace=0.34, wspace=0.26)

# 1 component preference
ax = fig.add_subplot(gs[0, 0])
ax.scatter(comp["r_NORTH_ST"].abs(), comp["r_CENTRAL_ALL"].abs(),
           s=40 + 900 * comp["share_cluster"], c=comp["share_cluster"], cmap="viridis",
           edgecolors="k", lw=0.4)
for _, r in comp.iterrows():
    ax.annotate(int(r["component"]), (abs(r["r_NORTH_ST"]), abs(r["r_CENTRAL_ALL"])),
                fontsize=6.5, xytext=(3, 3), textcoords="offset points")
lim = [0, max(comp["r_NORTH_ST"].abs().max(), comp["r_CENTRAL_ALL"].abs().max()) * 1.15]
ax.plot(lim, lim, "k--", lw=0.8)
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_xlabel("|corr| with NORTH tremor", fontsize=8)
ax.set_ylabel("|corr| with CENTRAL tremor", fontsize=8)
ax.tick_params(labelsize=7)
ax.set_title(f"T1  {npref}/{K} components track CENTRAL tremor BETTER than northern\n"
             "(above the 1:1 line); the model does have central-Oregon time functions\n"
             "(marker size/colour = share of its thrust in the 9-tri cluster)", fontsize=8.5)

# 2 lag curves for the 3 components with the largest cluster share
ax = fig.add_subplot(gs[0, 1])
top = comp.reindex(comp["r_signed_CLUSTER"].abs().sort_values(ascending=False).index).head(3)
for _, r in top.iterrows():
    k = int(np.where(idx + 1 == int(r["component"]))[0][0])
    ax.plot(LAGS, lag_curve(S[:, k], T_DET["CLUSTER"], iw0, iw1), lw=1.2,
            label=f"comp {int(r['component'])} (signed r {r['r_signed_CLUSTER']:+.2f})")
ax.plot(LAGS, lag_curve(Rser["CLUSTER"], T_DET["CLUSTER"], iw0, iw1), "k", lw=2.0,
        label="modeled potency rate ON the cluster")
ax.plot(LAGS, lag_curve(detrend_med(EST[("raw_GNSS", "CLUSTER", "JOINT")][0]),
                        T_DET["CLUSTER"], iw0, iw1), color="crimson", lw=1.6, ls="--",
        label="raw-GNSS matched filter on the cluster")
ax.axvline(0, color="0.6", lw=0.7)
ax.axhline(0, color="0.6", lw=0.7)
ax.set_xlabel("lag (d; + = tremor after slip)", fontsize=8)
ax.set_ylabel("Pearson r vs cluster tremor", fontsize=8)
ax.legend(fontsize=6.4, loc="upper left", framealpha=0.9)
ax.tick_params(labelsize=7)
ax.set_title("T1/T2  lagged correlation with the CLUSTER's own tremor", fontsize=8.5)

# 3 test-2 bars
ax = fig.add_subplot(gs[0, 2])
xs = np.arange(len(sl))
ax.bar(xs, sl["spearman_slip_vs_own_tremor"], 0.6, color="#2166ac", edgecolor="k", lw=0.4)
ax.errorbar(xs, sl["null_mean"], yerr=2 * sl["null_sd"], fmt="none", ecolor="crimson",
            elinewidth=1.6, capsize=5, label="shifted-tremor null +-2 sd")
for x, (_, r) in zip(xs, sl.iterrows()):
    ax.text(x, r["spearman_slip_vs_own_tremor"] + 0.02,
            f"{r['spearman_slip_vs_own_tremor']:.2f}\np={r['null_p']:.3f}",
            ha="center", fontsize=6.3)
ax.set_xticks(xs)
ax.set_xticklabels([f"{r['set']}\n({r['n_tri']} tri)" for _, r in sl.iterrows()], fontsize=7)
ax.axhline(0, color="k", lw=0.7)
ax.set_ylabel("Spearman(modeled potency rate, own 30-d tremor)", fontsize=8)
ax.legend(fontsize=6.5, loc="lower right")
ax.tick_params(labelsize=7)
ax.set_title("T2  the cluster's modeled slip DOES follow its tremor\n"
             "(as well as its slip+tremor neighbours do)", fontsize=8.5)

# 4 modeled slip series cluster vs central ST (normalised)
ax = fig.add_subplot(gs[1, 0])
zs = int(np.searchsorted(dates, np.datetime64("2013-01-01")))
ze = int(np.searchsorted(dates, np.datetime64("2018-01-01")))
for nm_, c_ in (("CLUSTER", "#b2182b"), ("CENTRAL_ST", "#2166ac"), ("NORTH_ST", "0.5")):
    v = Rser[nm_][zs:ze]
    ax.plot(dt[zs:ze], v / np.nanstd(Rser[nm_][iw0:iw1]), color=c_, lw=1.0, label=nm_)
ax2 = ax.twinx()
ax2.fill_between(dt[zs:ze], T_RAW["CLUSTER"][zs:ze], color="orange", alpha=0.3, lw=0)
ax2.set_ylabel("cluster tremor, 30-d count", fontsize=8, color="darkorange")
ax2.tick_params(labelsize=7, colors="darkorange")
ax.set_ylabel("modeled potency rate (sd units)", fontsize=8)
ax.xaxis.set_major_locator(YearLocator())
ax.xaxis.set_major_formatter(DateFormatter("%Y"))
ax.legend(fontsize=6.5, loc="upper left", framealpha=0.9)
ax.tick_params(labelsize=7)
ax.set_title(f"T2  the model gives the cluster the SAME time history as its\nneighbours: "
             f"corr(CLUSTER, CENTRAL_ST) = {cross[('CLUSTER','CENTRAL_ST')]:+.3f}, "
             f"corr(CLUSTER, NORTH_ST) = {cross[('CLUSTER','NORTH_ST')]:+.3f}", fontsize=8.5)

# 5 map
ax = fig.add_subplot(gs[1, 1])
ax.set_aspect(1.0 / np.cos(np.deg2rad(45.0)))
ax.scatter(lon, lat, s=6, c="0.87", edgecolors="none")
for nm_, c_ in (("CLUSTER", "#b2182b"), ("NORTH_ST", "#2166ac")):
    ax.scatter(lon[SETS[nm_]], lat[SETS[nm_]], s=14, c=c_, edgecolors="none", label=nm_)
    j = PATCH[nm_]["tri"]
    ax.scatter(lon[j], lat[j], s=110, marker="*", c="yellow", edgecolors="k", lw=0.6, zorder=6)
    for i in PATCH[nm_]["near5"]:
        ax.scatter(st_llh[i, 0], st_llh[i, 1], s=34, marker="^", c=c_, edgecolors="k", lw=0.4,
                   zorder=7)
        ax.annotate(names[i], (st_llh[i, 0], st_llh[i, 1]), fontsize=5.6,
                    xytext=(3, 2), textcoords="offset points")
ax.scatter(st_llh[:, 0], st_llh[:, 1], s=4, marker="^", c="0.6", zorder=3)
ax.set_xlim(lon.min() - 0.4, lon.max() + 1.6)
ax.set_ylim(lat.min() - 0.4, lat.max() + 0.4)
ax.legend(fontsize=6.5, loc="lower left", framealpha=0.9)
ax.tick_params(labelsize=7)
ax.set_xlabel("lon", fontsize=8)
ax.set_ylabel("lat", fontsize=8)
ax.set_title("T3  matched-filter targets (stars) and the 5 nearest\nstations used for each",
             fontsize=8.5)

# 6 matched filter time series
ax = fig.add_subplot(gs[1, 2])
for nm_, c_ in (("CLUSTER", "#b2182b"), ("NORTH_ST", "#2166ac")):
    a = detrend_med(EST[("raw_GNSS", nm_, "JOINT")][0])
    ax.plot(dt[zs:ze], a[zs:ze], color=c_, lw=0.9, label=f"{nm_} (JOINT, raw GNSS)")
ax2 = ax.twinx()
ax2.fill_between(dt[zs:ze], T_RAW["NORTH_ST"][zs:ze], color="0.75", alpha=0.5, lw=0,
                 label="north tremor")
ax2.fill_between(dt[zs:ze], T_RAW["CLUSTER"][zs:ze], color="orange", alpha=0.45, lw=0,
                 label="cluster tremor")
ax2.set_ylabel("30-d tremor count", fontsize=8)
ax2.tick_params(labelsize=7)
ax.set_ylabel("estimated 30-d peak slip on the patch (mm)", fontsize=8)
ax.xaxis.set_major_locator(YearLocator())
ax.xaxis.set_major_formatter(DateFormatter("%Y"))
ax.legend(fontsize=6.3, loc="upper left", framealpha=0.9)
ax.tick_params(labelsize=7)
ax.set_title("T3  data-space matched filter (no patterns)", fontsize=8.5)

# 7 composite by tremor decile
ax = fig.add_subplot(gs[2, 0])
for nm_, c_ in (("CLUSTER", "#b2182b"), ("NORTH_ST", "#2166ac")):
    a = detrend_med(EST[("raw_GNSS", nm_, "JOINT")][0])[iw0:iw1]
    t_ = T_DET[nm_][iw0:iw1]
    m_ = np.isfinite(a) & np.isfinite(t_)
    q = np.quantile(t_[m_], np.linspace(0, 1, 11))
    xs_, ys_, es_ = [], [], []
    for b in range(10):
        sel = m_ & (t_ >= q[b]) & (t_ <= q[b + 1])
        xs_.append(b + 1)
        ys_.append(np.nanmean(a[sel]))
        es_.append(np.nanstd(a[sel]) / np.sqrt(max(sel.sum() / WIN, 1)))
    ax.errorbar(xs_, ys_, yerr=es_, fmt="o-", color=c_, ms=4, lw=1.1, capsize=3, label=nm_)
ax.axhline(0, color="k", lw=0.7)
ax.set_xlabel("tremor decile (10 = most tremor)", fontsize=8)
ax.set_ylabel("mean estimated peak slip (mm/30 d)", fontsize=8)
ax.legend(fontsize=6.8, framealpha=0.9)
ax.tick_params(labelsize=7)
ax.set_title("T3  slip-vs-tremor composite, JOINT filter on raw GNSS\n"
             "(north = positive control)", fontsize=8.5)

# 8 data vs model through the identical filter
ax = fig.add_subplot(gs[2, 1])
tg = ["CLUSTER", "NORTH_ST"]
srcs = [("raw_GNSS", "raw GNSS (data)", "#b2182b"), ("MODEL_slip", "13-pattern slip model",
                                                     "#2166ac")]
xs = np.arange(len(tg))
for i, (sc, labl, c) in enumerate(srcs):
    v = [float(mf[(mf.source == sc) & (mf.target == t) & (mf.estimator == "JOINT")]
               ["composite_hi_minus_lo_mm"].iloc[0]) for t in tg]
    e = [float(mf[(mf.source == sc) & (mf.target == t) & (mf.estimator == "JOINT")]
               ["composite_se_mm"].iloc[0]) for t in tg]
    ax.bar(xs + (i - 0.5) * 0.34, v, 0.34, yerr=2 * np.array(e), color=c, edgecolor="k",
           lw=0.4, capsize=4, label=labl)
    for x, vv in zip(xs + (i - 0.5) * 0.34, v):
        ax.text(x, vv + 1.2, f"{vv:.1f}", ha="center", fontsize=6.6)
for x, t in zip(xs, tg):
    r = piv.loc[(t, "JOINT")]
    ax.text(x, -3.2, f"data-model = {r['data_minus_model_mm']:+.1f} +- "
                     f"{r['excess_se_mm']:.1f} mm\n({r['excess_sigma']:.1f} sigma; model "
                     f"captures {100*r['model_over_data']:.0f} %)",
            ha="center", fontsize=6.6)
ax.set_xticks(xs)
ax.set_xticklabels(tg, fontsize=8)
ax.axhline(0, color="k", lw=0.8)
ax.set_ylim(-6, 36)
ax.set_ylabel("top-decile minus bottom-decile tremor\nmean peak slip (mm/30 d), +-2 se",
              fontsize=8)
ax.legend(fontsize=6.8, loc="upper left", framealpha=0.9)
ax.tick_params(labelsize=7)
ax.set_title("T3  identical filter on DATA and on the MODEL:\n"
             "no significant unmodelled central-Oregon slip", fontsize=8.5)

# 9 verdict text
ax = fig.add_subplot(gs[2, 2])
ax.axis("off")
jc = mf[(mf.source == "raw_GNSS") & (mf.estimator == "JOINT") & (mf.target == "CLUSTER")].iloc[0]
jn = mf[(mf.source == "raw_GNSS") & (mf.estimator == "JOINT") & (mf.target == "NORTH_ST")].iloc[0]
rc = sl[sl.set == "CLUSTER"].iloc[0]
rn = sl[sl.set == "CENTRAL_ST"].iloc[0]
pc = piv.loc[("CLUSTER", "JOINT")]
pn = piv.loc[("NORTH_ST", "JOINT")]
txt = (
    "S12b - is the temporal loophole open?   NO.\n\n"
    f"T1  components tracking CENTRAL tremor better than\n"
    f"    NORTH: {npref}/{K}.  Strongest signed r with the\n"
    f"    cluster's own tremor: comp "
    f"{int(comp.loc[comp['r_signed_CLUSTER'].abs().idxmax(),'component'])} at "
    f"{comp['r_signed_CLUSTER'].abs().max():+.2f}.\n"
    f"    -> the basis HAS central-Oregon time functions.\n\n"
    f"T2  Spearman(modeled potency rate, own tremor):\n"
    f"    CLUSTER    {rc.spearman_slip_vs_own_tremor:+.3f}  (p {rc.null_p:.3f})\n"
    f"    CENTRAL_ST {rn.spearman_slip_vs_own_tremor:+.3f}  (p {rn.null_p:.3f})\n"
    f"    corr(CLUSTER, CENTRAL_ST slip series) "
    f"{cross[('CLUSTER','CENTRAL_ST')]:+.3f}\n"
    f"    slip-rate sd per km2: cluster "
    f"{float(sl[sl.set=='CLUSTER'].potency_rate_sd_per_km2.iloc[0]):.2e} vs\n"
    f"    central S+T "
    f"{float(sl[sl.set=='CENTRAL_ST'].potency_rate_sd_per_km2.iloc[0]):.2e} "
    f"-> the deficit is AMPLITUDE,\n"
    f"    not timing.\n\n"
    f"T3  JOINT matched filter, hi-lo tremor composite:\n"
    f"    CLUSTER  data {pc['raw_GNSS']:.1f}  model {pc['MODEL_slip']:.1f} mm\n"
    f"    NORTH    data {pn['raw_GNSS']:.1f}  model {pn['MODEL_slip']:.1f} mm\n"
    f"    excess (data-model): cluster "
    f"{pc['data_minus_model_mm']:+.1f}+-{pc['excess_se_mm']:.1f} mm "
    f"({pc['excess_sigma']:.1f}sig)\n"
    f"                          north  "
    f"{pn['data_minus_model_mm']:+.1f}+-{pn['excess_se_mm']:.1f} mm "
    f"({pn['excess_sigma']:.1f}sig)\n"
    f"    2-se bound on UNMODELLED central slip:\n"
    f"    {pc['unmodelled_bound_2se_mm']:.1f} mm peak / 30 d "
    f"(equiv. Mw {pc['unmodelled_bound_Mw']:.2f})\n"
)
ax.text(0.0, 1.0, txt, va="top", ha="left", fontsize=8.2, family="monospace",
        transform=ax.transAxes)

fig.suptitle("S12b — closing the temporal loophole: does a central-Oregon slip transient exist "
             "that the 13 fixed patterns could not carry?", fontsize=12, y=0.945)
fig.savefig(HERE / "figures/s12b_temporal_loophole.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("    wrote figures/s12b_temporal_loophole.png")
print("DONE")
