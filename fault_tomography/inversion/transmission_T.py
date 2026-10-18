"""B4 transmission normalization, steps 1-3 of Merlin's prescription (2026-08-22 audit, design note).

T(t, patch) = projection coefficient <recovered, pattern>/<pattern, pattern> of a noise-free +1%
Gaussian-taper blob injected through the REAL per-day estimator (day_solve: joint [m,s], production
lam, capture-weighted Gc, c=0) — i.e. per-day transmission incl. site-term absorption. Computed for
EVERY gate-passing day (no sampling), patches multiplexed odd/even per depth band (4 solves/day).

Also: (step 2) null / 2x-linearity / single-vs-multiplex checks on 3 days; (step 3) full patch-to-
patch transfer matrices T_qp at ~2010/~2016/~2023 (single-patch injections).

SCOPE (Merlin item 2b): T covers ESTIMATOR transmission only, not measurement-stage stretch railing
(fine at the 0.05-0.2% class; >0.5% anomalies would rail upstream). Per item 5 the c/BCD channel is
NOT covered here — per-day T removes roster-growth time-dependence; amplitude claims need the
separate end-to-end BCD injection (step 4, not this script).

Outputs (RESDIR): transmission_T.npz  (T [T x P], days, patch table, Tqp [3 x P x P], checks)
                  transmission_T.png  (per-patch smoothed T + Tqp heatmaps)
"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import invert_dvv_4d_daily as inv                    # module import runs the full data load

RD = inv.OUT
WORKERS = int(os.environ.get("WORKERS", "24"))
SEG_KM = 125.0                                       # along-strike patch width
SIG_KM = SEG_KM / 2.355                              # blob FWHM = patch width
DEEP_KM = 30.0
AMP = 1.0                                            # inject +1% (true dbeta/beta sign)

zprod = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True)
LAMF, LAMS = float(zprod["lam_f"]), float(zprod["lam_s"])
dep = np.asarray(zprod["depth"], float)
cxy = inv.cxy
T = inv.T

# ---- patches: PCA along-strike bins x depth band ----
a = cxy - cxy.mean(0)
u = np.linalg.eigh(a.T @ a)[1][:, -1]                # principal (along-strike) axis
s = a @ u
sbin = np.floor((s - s.min()) / SEG_KM).astype(int)
deepc = dep >= DEEP_KM
patches, patterns = [], []
for band, bmask in (("shallow", ~deepc), ("deep", deepc)):
    for b in np.unique(sbin[bmask]):
        cells = np.where(bmask & (sbin == b))[0]
        if len(cells) < 5:
            continue
        s_c = s[cells].mean()
        pat = np.zeros(inv.Mf)
        pat[bmask] = AMP * np.exp(-0.5 * ((s[bmask] - s_c) / SIG_KM) ** 2)
        pat[pat < 0.01 * AMP] = 0.0                  # compact support
        patches.append(dict(band=band, sbin=int(b), s_c=float(s_c), n=len(cells),
                            parity=int(b % 2), lat=float(zprod["clat"][cells].mean()),
                            depth=float(dep[cells].mean())))
        patterns.append(pat)
patterns = np.array(patterns); P = len(patches)
pt = pd.DataFrame(patches)
selfdot = (patterns ** 2).sum(1)
print(f"patches: {P} ({(pt.band=='deep').sum()} deep, {(pt.band=='shallow').sum()} shallow), "
      f"seg {SEG_KM:.0f} km, lam_f {LAMF:.3f}")

# v2 (Merlin round-2, 2026-08-22): SINGLE-PATCH injections for the daily table (multiplex retired —
# 9% same-parity leakage at weak-diagonal patches; the multiplexed table is archived as
# transmission_T_multiplexed.npz). Daily diagonal is now definitionally consistent with diag(T_qp).


def _recover(t, m_inj):
    idx = inv.BY_DAY[t]
    d = inv.Gc[inv.ROW[idx]] @ m_inj
    m, _, _, _ = inv.day_solve(t, np.zeros(inv.P), LAMF, LAMS, d_override=d)
    return m


def _day_T(t):
    out = np.full(P, np.nan)
    for p in range(P):
        m = _recover(t, patterns[p])
        out[p] = (patterns[p] @ m) / selfdot[p]
    return out


if __name__ == "__main__":
    days = inv.days
    REP_DATES = ("2009-06-15", "2010-06-15", "2016-06-15", "2023-06-15")   # 2009 added (Merlin r2)
    rep = [int(days.get_indexer([pd.Timestamp(x)], method="nearest")[0]) for x in REP_DATES]

    # ---- null / linearity checks on the representative days ----
    checks = {}
    for t in rep:
        m0 = _recover(t, np.zeros(inv.Mf))
        checks[f"null_{t}"] = float(np.abs(m0).max())
        allpat = patterns.sum(0)                      # any-coverage pattern (2009 has no far-N rows)
        m1 = _recover(t, allpat)
        m2 = _recover(t, 2.0 * allpat)
        nz = np.abs(m1) > 1e-9
        checks[f"lin_{t}"] = float(np.nanmax(np.abs(m2[nz] / m1[nz] - 2.0))) if nz.any() else np.nan
    print("CHECK null max|m|:", [f"{checks[f'null_{t}']:.2e}" for t in rep])
    print("CHECK 2x-linearity max dev:", [f"{checks[f'lin_{t}']:.2e}" for t in rep])

    # ---- transfer matrices T_qp (single-patch injections) at the rep epochs ----
    Tqp = np.zeros((len(rep), P, P))
    for k, t in enumerate(rep):
        for p in range(P):
            m = _recover(t, patterns[p])
            Tqp[k, :, p] = (patterns @ m) / selfdot

    # ---- step 1: T for every day ----
    with ProcessPoolExecutor(WORKERS) as ex:
        rows = list(ex.map(_day_T, range(T), chunksize=16))
    Tday = np.array(rows)                             # (T, P)

    np.savez_compressed(f"{RD}/transmission_T.npz",
                        T=Tday.astype(np.float32), days=days.astype(str).values,
                        patch_band=pt.band.values.astype(object), patch_lat=pt.lat.values,
                        patch_depth=pt.depth.values, patch_sc=pt.s_c.values,
                        patch_parity=pt.parity.values, patterns=patterns.astype(np.float32),
                        Tqp=Tqp.astype(np.float32), rep_days=days[rep].astype(str).values,
                        seg_km=SEG_KM, lam_f=LAMF, lam_s=LAMS, amp=AMP,
                        checks=np.array([f"{k}={v:.3e}" for k, v in checks.items()], object))

    # ---- figure ----
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(13, 9))
    Ts = pd.DataFrame(Tday, index=days).rolling("90D", min_periods=10).median()
    for i, bd in enumerate(("shallow", "deep")):
        ax = fig.add_subplot(2, 3, (1, 2) if i == 0 else (4, 5))
        qi = np.where(pt.band == bd)[0]
        latv = pt.lat.values
        cm = plt.cm.viridis((latv[qi] - latv.min()) / (latv.max() - latv.min() + 1e-9))
        for q, c in zip(qi, cm):
            ax.plot(Ts.index, Ts[q], lw=0.9, color=c)
        ax.set_ylabel(f"T {bd}"); ax.set_ylim(0, 1); ax.axhline(0.1, color="r", lw=0.5, ls="--")
        ax.set_title(f"{bd}: per-patch daily T (90-d median), color=latitude", fontsize=10)
    # heatmaps in their own figure
    fig2, ax2 = plt.subplots(1, len(rep), figsize=(4.3 * len(rep), 4))
    order = np.argsort(pt.band.values.astype(str) + pt.s_c.map("{:09.1f}".format).values)
    for k in range(len(rep)):
        im = ax2[k].imshow(Tqp[k][np.ix_(order, order)], vmin=0, vmax=0.8, cmap="magma")
        ax2[k].set_title(f"T_qp {str(days[rep[k]].date())}", fontsize=9)
        nsh = (pt.band == "shallow").sum()
        ax2[k].axhline(nsh - 0.5, color="w", lw=0.5); ax2[k].axvline(nsh - 0.5, color="w", lw=0.5)
    fig2.colorbar(im, ax=ax2, shrink=0.8, label="T_qp (rows=recovered, cols=injected)")
    fig.savefig(f"{RD}/transmission_T.png", dpi=110, bbox_inches="tight")
    fig2.savefig(f"{RD}/transmission_Tqp.png", dpi=110, bbox_inches="tight")

    dd = pt.band == "deep"
    yr = pd.Series(Tday[:, np.where(dd)[0]].mean(1), index=days).resample("YE").mean()
    print("deep-band mean T by year:")
    print(yr.round(3).to_string())
    print(f"wrote {RD}/transmission_T.npz + transmission_T.png + transmission_Tqp.png")
