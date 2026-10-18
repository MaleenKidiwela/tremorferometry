"""Checkerboard resolution of the DAILY estimator (Merlin-audited redo, 2026-08-23 — supersedes the
July study FOR THE DAILY INTERFACE PRODUCT ONLY; the July 3-D volume verdict (0.85-0.86) and the
old-estimator record stand; the sensitivity atlas is marked STALE-FOR-DAILY, ticketed separately).

Design (all rulings in notes/2026-08-22_daily_joint_inversion_design.md):
- scales = July ladder {44,70,100,140,180,250} km, amplitude 1% true dbeta/beta, bar corr>=0.7.
- estimator = production day_solve (lam_f 0.215, ratio 0.125, c=0); linearity lets signal and noise
  be solved separately and summed.
- injection configs per scale: shallow-only, deep-only (HEADLINE per band), joint-opposite-sign
  (adversarial: down-dip leak fights the truth), joint-same-sign (demo of leak inflation ONLY —
  a map-view checker gives both bands the same sign, so same-sign joint FLATTERS deep recovery).
- noise = iid 1/sqrt(w) + station-common (I5b; NOTE SIG_COMMON is a std over time of station-mean
  dvv and so contains real seasonality — red common-mode that averages slower than 1/sqrt(Neff);
  covered by the tau bracket below). R=8 draws/epoch, seeded, reused across scales/configs.
- temporal averaging W in {1 frame, 90 d, 365 d, 1825 d} via NOISE SCALING 1/sqrt(Neff),
  Neff = max(W/tau, 1), bracketed tau in {30, 75, 150} d (75 = project convention, UNVALIDATED —
  if the W=365 verdict flips across the bracket, the split-half empirical measurement gates these
  rows). Mean response for W=365/1825 = average of noise-free recoveries over 3/5 days spanning the
  window (captures roster drift exactly); W<=90 uses the epoch day alone.
- epochs 2011/2013/2016/2019/2023/2026 (2009-10 excluded, user ruling; 2011 = earliest characterized).
- CAVEATS (printed on figure+CSV): corr is amplitude-blind (amplitude fidelity = B4 T/deconv product);
  the 1% amplitude defines the SNR — at the realistic 0.1-0.2% class the noise tier is 5-10x higher;
  long-W rows are geometry+noise only — persistent-field amplitude additionally shrinks via the
  family-era offset channel (B4 end-to-end factor); W=1825 not usable standalone.
- BRIDGE ASSERTION: deep-only, W=365, tau=75, epoch 2023 must land within one scale step of the July
  70-140 km verdict, else HALT and report (do not publish a "degraded" verdict; suspect the noise
  scaling first).

Outputs: res_catalog/checkerboard_daily.{npz,csv,png}.
"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import invert_dvv_4d_daily as inv

RD = inv.OUT
zp = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True)
LAMF, LAMS = float(zp["lam_f"]), float(zp["lam_s"])
dep = np.asarray(zp["depth"], float)
deepc = dep >= 30.0
bands = {"shallow": ~deepc, "deep": deepc}
cx, cy = inv.cxy[:, 0], inv.cxy[:, 1]

SCALES = [44, 70, 100, 140, 180, 250]
TAUS = [30, 75, 150]
WINS = [("frame", 1, 1), ("90d", 90, 1), ("1y", 365, 3), ("5y", 1825, 5)]   # (label, days, n_sample)
EPOCH_DATES = ["2011-06-15", "2013-06-15", "2016-06-15", "2019-06-15", "2023-06-15", "2026-06-01"]
NDRAW = 8
BAR = 0.7

days = inv.days
rep = [int(days.get_indexer([pd.Timestamp(x)], method="nearest")[0]) for x in EPOCH_DATES]


def checker(L):
    return np.sign(np.sin(np.pi * cx / L) * np.sin(np.pi * cy / L))


def build_inj(L, config):
    ck = checker(L)
    if config == "shallow_only": return ck * bands["shallow"]
    if config == "deep_only":    return ck * bands["deep"]
    if config == "joint_same":   return ck.copy()
    if config == "joint_opp":    return ck * bands["shallow"] - ck * bands["deep"]
    raise ValueError(config)


def solve_day(t, dvec):
    m, _, _, _ = inv.day_solve(t, np.zeros(inv.P), LAMF, LAMS, d_override=dvec)
    return m


def _mean_resp(argt):
    """noise-free mean response of injection m_inj over sampled days."""
    tds, m_inj = argt
    acc = np.zeros(inv.Mf)
    for td in tds:
        idx = inv.BY_DAY[td]
        acc += solve_day(td, inv.Gc[inv.ROW[idx]] @ m_inj)
    return acc / len(tds)


def _noise_resp(argt):
    t, seed = argt
    rng = np.random.RandomState(seed)
    idx = inv.BY_DAY[t]
    n = rng.randn(len(idx)) / np.sqrt(inv.Wr[idx]) + \
        (inv.SIG_COMMON * rng.randn(inv.NS))[inv.STJ[idx]]
    return solve_day(t, n)


def window_days(t, wdays, k):
    if k == 1: return [t]
    c = days[t]
    lo, hi = c - pd.Timedelta(days=wdays // 2), c + pd.Timedelta(days=wdays // 2)
    lo = max(lo, days.min()); hi = min(hi, days.max())
    samp = pd.date_range(lo, hi, periods=k)
    return sorted(set(int(days.get_indexer([s], method="nearest")[0]) for s in samp))


if __name__ == "__main__":
    CONFIGS = ["shallow_only", "deep_only", "joint_same", "joint_opp"]
    # ---- noise responses: (epoch, draw), reused across scales & configs ----
    ntasks = [(t, 1000 * k + r) for k, t in enumerate(rep) for r in range(NDRAW)]
    with ProcessPoolExecutor(24) as ex:
        nres = list(ex.map(_noise_resp, ntasks, chunksize=1))
    NOISE = {}  # (epoch_idx, r) -> m
    for (t, seed), m in zip(ntasks, nres):
        k = rep.index(t); NOISE[(k, seed % 1000)] = m

    # ---- mean responses: (epoch, scale, config, window) ----
    mtasks, keys = [], []
    for k, t in enumerate(rep):
        for L in SCALES:
            for cfg in CONFIGS:
                m_inj = build_inj(L, cfg)
                for wl, wd, ns in WINS:
                    tds = window_days(t, wd, ns)
                    keys.append((k, L, cfg, wl))
                    mtasks.append((tds, m_inj))
    # dedupe identical (tds, injection) pairs: frame/90d share tds=[t]
    with ProcessPoolExecutor(24) as ex:
        mres = list(ex.map(_mean_resp, mtasks, chunksize=2))
    MEAN = dict(zip(keys, mres))

    # ---- assemble table ----
    rows = []
    for k, t in enumerate(rep):
        idx = inv.BY_DAY[t]
        ill = np.zeros(inv.Mf, bool); ill[np.unique(inv.pcell_idx[inv.ROW[idx]])] = True
        for L in SCALES:
            ck = checker(L)
            for cfg in CONFIGS:
                m_inj = build_inj(L, cfg)
                for wl, wd, ns in WINS:
                    mr = MEAN[(k, L, cfg, wl)]
                    for tau in TAUS:
                        neff = max(wd / tau, 1.0)
                        for band, bm in bands.items():
                            truth = m_inj
                            sel = bm & ill & (truth != 0)
                            if sel.sum() < 12: continue
                            cs = []
                            for r in range(NDRAW):
                                mh = mr + NOISE[(k, r)] / np.sqrt(neff)
                                cs.append(np.corrcoef(truth[sel], mh[sel])[0, 1])
                            rows.append(dict(epoch=str(days[t].date())[:7], scale=L, config=cfg,
                                             band=band, win=wl, tau=tau,
                                             corr=float(np.mean(cs)), corr_sd=float(np.std(cs))))
    df = pd.DataFrame(rows)
    df.to_csv(f"{RD}/checkerboard_daily.csv", index=False)

    def crossing(sub):
        s = sub.sort_values("scale")
        ok = s[s.corr_ >= BAR] if "corr_" in s else s[s["corr"] >= BAR]
        return int(ok.scale.min()) if len(ok) else -1

    # headline crossings: per-band configs
    print("\nRESOLVABLE SCALE (first scale with mean corr >= 0.7; -1 = never), per-band injection:")
    for band, cfg in (("deep", "deep_only"), ("shallow", "shallow_only")):
        for wl, _, _ in WINS:
            line = []
            for k in range(len(rep)):
                sub = df[(df.band == band) & (df.config == cfg) & (df.win == wl) &
                         (df.tau == 75) & (df.epoch == str(days[rep[k]].date())[:7])]
                line.append(crossing(sub))
            print(f"  {band:8s} W={wl:6s} tau=75: {dict(zip([e[:7] for e in EPOCH_DATES], line))}")
    # tau sensitivity at W=1y, deep
    for tau in TAUS:
        line = []
        for k in range(len(rep)):
            sub = df[(df.band == "deep") & (df.config == "deep_only") & (df.win == "1y") &
                     (df.tau == tau) & (df.epoch == str(days[rep[k]].date())[:7])]
            line.append(crossing(sub))
        print(f"  tau-bracket deep W=1y tau={tau}: {line}")
    # adversarial vs flattering joint, deep, 2023, W=1y
    for cfg in ("deep_only", "joint_opp", "joint_same"):
        sub = df[(df.band == "deep") & (df.config == cfg) & (df.win == "1y") & (df.tau == 75) &
                 (df.epoch == "2023-06")]
        print(f"  deep 2023 W=1y {cfg}: crossing {crossing(sub)} | corr@100km "
              f"{float(sub[sub.scale == 100]['corr'].iloc[0]):.2f}")

    # data-space anchor (July convention): signal std of 1% checker at 100 km, deep-only, 2023
    t = rep[EPOCH_DATES.index('2023-06-15')]; idx = inv.BY_DAY[t]
    sig_ds = float(np.std(inv.Gc[inv.ROW[idx]] @ build_inj(100, "deep_only")))
    print(f"\nanchor: 1% deep 100-km checker -> data-space signal std {sig_ds:.4f}% "
          f"(0.043% class => rescale corr SNR by ~{0.043/max(sig_ds,1e-6):.1f}x amplitude equivalence; "
          f"realistic slow-field 0.1-0.2% => noise tier 5-10x the 1% table)")

    # BRIDGE ASSERTION
    sub = df[(df.band == "deep") & (df.config == "deep_only") & (df.win == "1y") &
             (df.tau == 75) & (df.epoch == "2023-06")]
    br = crossing(sub)
    ladder = {44: 0, 70: 1, 100: 2, 140: 3, 180: 4, 250: 5}
    okset = {44, 70, 100, 140, 180}          # within one step of the July 70-140 verdict
    verdict = "PASS" if br in okset else "FAIL - HALT, audit noise scaling before publishing"
    print(f"\nBRIDGE (deep, W=1y, tau=75, 2023): crossing {br} km vs July 70-140 -> {verdict}")

    np.savez_compressed(f"{RD}/checkerboard_daily.npz",
                        table=df.to_records(index=False), scales=SCALES, taus=TAUS,
                        epochs=[str(days[t].date()) for t in rep], bar=BAR, ndraw=NDRAW,
                        bridge_crossing=br, anchor_ds=sig_ds,
                        caveats="corr is amplitude-blind (B4 governs amplitude); long-W rows are "
                                "geometry+noise only (c-channel shrinkage excluded, see B4 end-to-end); "
                                "interface product only - July 3-D volume verdict untouched")

    # ---- figure: deep-band curves per epoch (per-band injection), tau=75, by window ----
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, len(WINS), figsize=(4.2 * len(WINS), 3.6), sharey=True)
    cm = plt.cm.viridis(np.linspace(0, 1, len(rep)))
    for ax, (wl, _, _) in zip(axs, WINS):
        for k, col in zip(range(len(rep)), cm):
            sub = df[(df.band == "deep") & (df.config == "deep_only") & (df.win == wl) &
                     (df.tau == 75) & (df.epoch == str(days[rep[k]].date())[:7])].sort_values("scale")
            ax.errorbar(sub.scale, sub["corr"], yerr=sub.corr_sd, color=col, lw=1.3,
                        label=str(days[rep[k]].year))
        ax.axhline(BAR, color="r", lw=0.7, ls="--")
        ax.set_title(f"W = {wl}", fontsize=10); ax.set_xlabel("checker scale (km)")
        ax.set_ylim(-0.1, 1.05)
    axs[0].set_ylabel("recovery corr (deep band)"); axs[0].legend(fontsize=7)
    fig.suptitle("DAILY-estimator deep checkerboard (per-band injection, tau=75 d, 1% amplitude; "
                 "corr is amplitude-blind — B4 governs amplitudes; 5-y window excludes c-channel shrinkage)",
                 fontsize=9)
    fig.savefig(f"{RD}/checkerboard_daily.png", dpi=110, bbox_inches="tight")
    print(f"wrote {RD}/checkerboard_daily.npz/.csv/.png")
