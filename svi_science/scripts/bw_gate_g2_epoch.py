#!/usr/bin/env python
"""bw_gate_g2_epoch -- G2: do the 2012 templates still fit VGZ and SNB on their OLDER instruments?

Track A-2 step-4 protocol (lfe_extension/results/TRACK_A2_RESULTS.md), with its own functions unchanged
(a4_epoch_stack.best_cc / registration, e5r_templates.cut_windows / align_stack / pick_s / s_snr_db, d3_vel1d.tt):
for every (family, station) with a usable 2012 template, cut the station's 3-C 20 Hz / 1-8 Hz data at the predicted S
of the reference events, rms-normalise, cc-align, stack; compare with the 2012 stack over S-1..S+7 s (signed 3-C cc
maximised over +-2 s -> cc, S shift dS); registration = median leave-one-out per-event cc.
Reference events = Bostock's catalogue events (truth time + frozen 2012 per-family offset), largest magnitude first,
<= 400 per family, >= 20 needed.  Control leg = the same station on 2012 Bostock events in the TEST blocks only (the
templates were built on the build blocks), each family capped at the new leg's count (matched N).
No-change controls = stations on the same instrument as in 2012 (PGC, NLLB, PB; VGZ in 2010).

  python bw_gate_g2_epoch.py run [--workers 12]     -> OUT/tables/g2_family_rows.csv
  python bw_gate_g2_epoch.py report                 -> OUT/tables/g2_summary.csv, g2_decision.json, figures/g2_*.png
  python bw_gate_g2_epoch.py rebuild --keys ...     -> build/test-block rebuild validation for failing epochs
"""
from __future__ import annotations
import os, sys, json, argparse, warnings
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bw_gate_common as G
sys.path.insert(0, G.LFE_SCRIPTS)
warnings.filterwarnings("ignore")

TAB = f"{G.OUT}/tables"
ROOT = G.ROOT
CACHES = {2004: [G.CACHE_BW], 2005: [G.CACHE_BW, f"{ROOT}/lfe_array_discovery/cache20hz"],
          2008: [G.CACHE_BW, f"{ROOT}/lfe_picker/cache20hz"],
          2010: [f"{G.REP}/cache20hz", f"{G.LFE}/cache20hz"], 2012: [f"{G.REP}/cache20hz"]}
WIN = {"2004-07": ("2004-07-02", "2004-07-27"), "2005-09": ("2005-09-03", "2005-09-25"),
       "2008-05": ("2008-05-10", "2008-06-02"), "2010-08": ("2010-08-16", "2010-09-08"),
       "2012-09": ("2012-08-30", "2012-09-25")}
# (window, station, role, epoch label)
LEGS = []
for w in ("2004-07", "2005-09"):
    LEGS += [(w, "VGZ", "target", "VGZ BH CMG3ESP/GD1 (2003-02..2007-02)"), (w, "SNB", "target", "SNB BH CMG40T/GD1 (..2006-04)"),
             (w, "PGC", "control", "PGC BH CMG3T/GD1 (no change)"), (w, "NLLB", "control", "NLLB BH CMG3ESP/GD2_v12 (no change)")]
LEGS += [("2005-09", s, "control", f"{s} EH (no change)") for s in G.OLY]
LEGS += [("2008-05", "VGZ", "target", "VGZ BH CMG3T/GD2_v21 (2007-02..2009-07)"), ("2008-05", "SNB", "target", "SNB BH CMG3ESP/GD2_v12 (2006-04..2008-12)"),
         ("2008-05", "PGC", "control", "PGC BH CMG3T/GD1 (no change)"), ("2008-05", "NLLB", "control", "NLLB BH CMG3ESP/GD2_v12 (no change)")]
LEGS += [("2008-05", s, "control", f"{s} EH (no change)") for s in G.OLY + G.PGCB + ["B926"]]
LEGS += [("2010-08", "SNB", "target", "SNB BH CMG3ESP/GD1 (2008-12..2011-03)"), ("2010-08", "PGC", "control", "PGC BH CMG3T/GD1 (no change)"),
         ("2010-08", "NLLB", "control", "NLLB BH CMG3ESP/GD2_v12 (no change)"), ("2010-08", "VGZ", "control", "VGZ HH CMG3ESP/GD2_v21 (2012 epoch)")]
LEGS += [("2010-08", s, "control", f"{s} EH (no change)") for s in G.PGCB + ["B926", "B001", "B003", "B004", "B005", "B006", "B007"]]

_STATE = {}


def setup():
    if _STATE: return _STATE
    import e5r_templates as T
    import d3_vel1d as V
    import a4_epoch_stack as A4
    stab = pd.read_csv(f"{G.REP}/results/e5r_stations.csv").set_index("sta")
    F = pd.read_csv(f"{G.REP}/results/e5r_families_2012.csv").set_index("bkey")
    QC = pd.read_csv(f"{G.REP}/results/e5r_template_qc_2012.csv")
    B = G.bostock_all()
    import svi_ref_common as R
    _STATE.update(T=T, V=V, A4=A4, stab=stab, F=F, QC=QC, B=B, R=R)
    return _STATE


def events(win, test_blocks=False, build_blocks=False):
    S = setup(); B = S["B"]; a, b = WIN[win]
    E = B[(B.date >= a) & (B.date <= b)].copy()
    if test_blocks or build_blocks:
        role = S["R"].block_role(E.t_origin.to_numpy())
        E = E[role == ("fit" if test_blocks else "insample")]
    return E


def stack_leg(win, sta, leg, caps=None, ev=None, max_events=400, min_events=20, keep_stacks=False):
    """one (window, station, leg) -> rows per family (a4_epoch_stack.main logic, Bostock events as reference)."""
    S = setup(); T, V, A4 = S["T"], S["V"], S["A4"]
    year = int(WIN[win][0][:4])
    T.CACHES = CACHES[year]; T.YEAR = year
    ns = G.NETSTA[sta]
    days = T.load_days(ns, year)
    lat_s, lon_s = float(S["stab"].loc[sta, "lat"]), float(S["stab"].loc[sta, "lon"])
    QC = S["QC"]; QC = QC[(QC.sta == sta) & QC.usable & QC.template].set_index("bkey")
    Z = np.load(f"{G.REP}/results/e5r_templates_2012.npz", allow_pickle=True)
    E = ev if ev is not None else events(win, test_blocks=(win == "2012-09"))
    rows, STK = [], {}
    nwin = int(round(2 * A4.HALF_W * A4.FS))
    rng = np.random.default_rng(7)
    for bkey in sorted(QC.index):
        g = E[E.fam == int(bkey)]
        cap = max_events if caps is None else min(max_events, caps.get(int(bkey), 0))
        if caps is not None and cap < min_events:
            continue
        base = dict(window=win, sta=sta, leg=leg, family=int(bkey), n_cat=len(g), n_days_cache=len(days))
        if len(g) < min_events or not days:
            rows.append(dict(base, n_events=len(g), status="too_few_events" if days else "no_cache")); continue
        g = g.sort_values("mag", ascending=False).head(cap)
        fr = S["F"].loc[bkey]
        dx = T.epi_km(float(fr.lat), float(fr.lon), lat_s, lon_s)
        ttS = float(V.tt(dx, float(fr.dep), "S"))
        Sw, keep = T.cut_windows(days, g.t_origin.to_numpy() + ttS - A4.HALF_W, nwin)
        if Sw is None or len(Sw) < min_events:
            rows.append(dict(base, n_events=0 if Sw is None else len(Sw), status="no_data")); continue
        st, Sal, sh = T.align_stack(Sw, iters=2)
        sp_idx = T.pick_s(st, A4.IDX_NOM)
        snr = T.s_snr_db(st, sp_idx)
        ref = np.asarray(Z[f"{sta}|{int(bkey)}"], float)
        i_ref = A4.IDX_NOM + int(round(float(QC.loc[bkey, "s_pick_s"]) * A4.FS)) - int(round(A4.PRE_TPL * A4.FS))
        cc, dS = A4.best_cc(ref, st, i_ref)
        reg = A4.registration(Sal, st, i_ref + int(round(dS * A4.FS))) if np.isfinite(dS) else np.nan
        rows.append(dict(base, n_events=int(len(Sw)), epi_km=dx, cc=cc, dS=dS, registration=reg, snr_db=snr,
                         s_pick_s=(sp_idx - A4.IDX_NOM) / A4.FS, s_pick_2012=float(QC.loc[bkey, "s_pick_s"]), status="ok"))
        if keep_stacks:
            STK[int(bkey)] = st.astype(np.float32)
    return (rows, STK) if keep_stacks else rows


def _job(args):
    win, sta, leg, caps = args
    try:
        return stack_leg(win, sta, leg, caps)
    except Exception as e:
        import traceback
        return [dict(window=win, sta=sta, leg=leg, status=f"ERR {type(e).__name__}: {str(e)[:120]}", tb=traceback.format_exc()[-400:])]


def run(a):
    from multiprocessing import Pool
    setup()
    G.log(f"[g2] {len(LEGS)} new legs")
    with Pool(a.workers) as pool:
        res = pool.map(_job, [(w, s, "new", None) for w, s, role, lab in LEGS], chunksize=1)
    N = pd.DataFrame([r for rr in res for r in rr])
    jobs = []
    for (w, s), g in N[N.status == "ok"].groupby(["window", "sta"]):
        jobs.append((f"2012-09", s, f"ctl_for_{w}", {int(r.family): int(r.n_events) for r in g.itertuples()}))
    G.log(f"[g2] {len(jobs)} matched-N 2012 control legs")
    with Pool(a.workers) as pool:
        res2 = pool.map(_job, jobs, chunksize=1)
    Cc = pd.DataFrame([r for rr in res2 for r in rr])
    A = pd.concat([N, Cc], ignore_index=True)
    lab = {(w, s): (role, l) for w, s, role, l in LEGS}
    A["role"] = [lab.get((w, s), ("", ""))[0] if leg == "new" else lab.get((leg.replace("ctl_for_", ""), s), ("", ""))[0]
                 for w, s, leg in zip(A.window, A.sta, A.leg)]
    A["epoch_label"] = [lab.get((w, s), ("", ""))[1] if leg == "new" else "2012 build epoch (control leg, test blocks)"
                        for w, s, leg in zip(A.window, A.sta, A.leg)]
    A.to_csv(f"{TAB}/g2_family_rows.csv", index=False)
    G.log(f"[g2] wrote {len(A)} rows; status {A.status.value_counts().to_dict()}")


def pair_rows(A):
    """new leg joined with its matched-N control leg per (window, station, family)."""
    N = A[(A.leg == "new") & (A.status == "ok")]
    Cc = A[A.leg.str.startswith("ctl_for_") & (A.status == "ok")].copy()
    Cc["window"] = Cc.leg.str.replace("ctl_for_", "")
    M = N.merge(Cc[["window", "sta", "family", "cc", "dS", "registration", "n_events"]], on=["window", "sta", "family"],
                suffixes=("", "_ctl"))
    M["reg_ratio"] = M.registration / M.registration_ctl
    return M


def station_summary(M):
    rows = []
    for (w, s), g in M.groupby(["window", "sta"]):
        rows.append(dict(window=w, sta=s, role=g.role.iloc[0], epoch=g.epoch_label.iloc[0], n_fam=len(g),
                         cc_med=g.cc.median(), cc_ctl_med=g.cc_ctl.median(), cc_ratio=g.cc.median() / g.cc_ctl.median(),
                         reg_ratio=g.reg_ratio.median(), dS_med=g.dS.median(), dS_mad=(g.dS - g.dS.median()).abs().median(),
                         flips=int((g.cc <= -0.6).sum()), flip_frac=float((g.cc <= -0.6).mean()),
                         n_events_med=g.n_events.median()))
    return pd.DataFrame(rows)


def ratios(g):
    return g.cc.median() / g.cc_ctl.median(), g.reg_ratio.median()


def boot_diff(gt, gc, nb=1000, seed=5):
    rs = np.random.default_rng(seed); out = []
    for _ in range(nb):
        a = gt.sample(len(gt), replace=True, random_state=int(rs.integers(1e9)))
        b = gc.sample(len(gc), replace=True, random_state=int(rs.integers(1e9)))
        ra, rb = ratios(a), ratios(b)
        out.append((ra[0] - rb[0], ra[1] - rb[1]))
    o = np.array(out)
    return [float(np.percentile(o[:, 0], 5)), float(np.percentile(o[:, 0], 95)),
            float(np.percentile(o[:, 1], 5)), float(np.percentile(o[:, 1], 95))]


def report(a):
    A = pd.read_csv(f"{TAB}/g2_family_rows.csv")
    M = pair_rows(A)
    M.to_csv(f"{TAB}/g2_pairs.csv", index=False)
    SS = station_summary(M); SS.to_csv(f"{TAB}/g2_station_summary.csv", index=False)
    rows, dec = [], {}
    tests = [("VGZ", ["2008-05"], "VGZ BH CMG3T/GD2_v21 2007-02..2009-07"),
             ("VGZ", ["2004-07"], "VGZ BH CMG3ESP/GD1 2003-02..2007-02"), ("VGZ", ["2005-09"], "VGZ BH CMG3ESP/GD1 2003-02..2007-02"),
             ("VGZ", ["2004-07", "2005-09"], "VGZ BH CMG3ESP/GD1 2003-02..2007-02"),
             ("SNB", ["2010-08"], "SNB BH CMG3ESP/GD1 2008-12..2011-03"), ("SNB", ["2008-05"], "SNB BH CMG3ESP/GD2_v12 2006-04..2008-12"),
             ("SNB", ["2004-07"], "SNB BH CMG40T/GD1 ..2006-04"), ("SNB", ["2005-09"], "SNB BH CMG40T/GD1 ..2006-04"),
             ("SNB", ["2004-07", "2005-09"], "SNB BH CMG40T/GD1 ..2006-04")]
    for sta, wins, ep in tests:
        gt = M[(M.sta == sta) & M.window.isin(wins) & (M.role == "target")]
        gc = M[M.window.isin(wins) & (M.role == "control")]
        if not len(gt) or not len(gc):
            continue
        rt, rc = ratios(gt), ratios(gc)
        ci = boot_diff(gt, gc)
        ok_cc = rt[0] >= rc[0] - 0.05; ok_reg = rt[1] >= rc[1] - 0.05
        ok_ds = abs(gt.dS.median()) <= 0.10; ok_fl = (gt.cc <= -0.6).mean() <= 0.10
        r = dict(station=sta, epoch=ep, windows="+".join(wins), n_fam=len(gt), n_ctl_rows=len(gc),
                 ctl_stations=",".join(sorted(gc.sta.unique())), cc_ratio=rt[0], cc_ratio_nochange=rc[0],
                 reg_ratio=rt[1], reg_ratio_nochange=rc[1], d_cc_ratio_90ci=ci[:2], d_reg_ratio_90ci=ci[2:],
                 cc_med=gt.cc.median(), cc_ctl_med=gt.cc_ctl.median(), dS_med=gt.dS.median(),
                 dS_mad=float((gt.dS - gt.dS.median()).abs().median()), flip_frac=float((gt.cc <= -0.6).mean()),
                 pass_cc=bool(ok_cc), pass_reg=bool(ok_reg), pass_dS=bool(ok_ds), pass_flips=bool(ok_fl))
        r["transfer_pass"] = bool(ok_cc and ok_reg and ok_ds and ok_fl)
        rows.append(r)
        G.log(f"[g2] {sta} {'+'.join(wins)}: cc ratio {rt[0]:.3f} vs no-change {rc[0]:.3f}; reg ratio {rt[1]:.3f} vs {rc[1]:.3f}; "
              f"dS {gt.dS.median():+.2f}; flips {r['flip_frac']:.2f} -> {'PASS' if r['transfer_pass'] else 'FAIL'}")
    T = pd.DataFrame(rows); T.to_csv(f"{TAB}/g2_summary.csv", index=False)
    json.dump(T.to_dict("records"), open(f"{TAB}/g2_decision.json", "w"), indent=1, default=float)
    figure(M, SS)


def figure(M, SS):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    wins = ["2004-07", "2005-09", "2008-05", "2010-08"]
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
    for k, (col, lab) in enumerate((("cc_ratio", "cc ratio (old epoch / 2012 control leg)"), ("reg_ratio", "registration ratio"))):
        ax = axs[k]
        for i, w in enumerate(wins):
            g = SS[SS.window == w]
            c = g[g.role == "control"]; t = g[g.role == "target"]
            ax.scatter(np.full(len(c), i) + np.linspace(-.15, .15, max(len(c), 1))[:len(c)], c[col], s=22, c="0.55",
                       label="no-change stations" if i == 0 else None)
            for _, r in t.iterrows():
                ax.scatter(i + (0.25 if r.sta == "SNB" else -0.25), r[col], s=70, marker="D" if r.sta == "SNB" else "s",
                           c="#c0392b" if r.sta == "SNB" else "#2471a3", label=r.sta if i == wins.index("2008-05") else None)
        ax.set_xticks(range(len(wins))); ax.set_xticklabels(wins); ax.set_ylabel(lab); ax.axhline(1, c="k", lw=.6)
        ax.set_title(lab); ax.legend(fontsize=8)
    plt.tight_layout(); plt.savefig(f"{G.OUT}/figures/g2_transfer.png", dpi=130); plt.close()


# ---------------------------------------------------------------------------------------------------- rebuild
def rebuild(a):
    """for each failing (station, window): templates rebuilt from Bostock BUILD-block events of that window,
    validated on its TEST blocks (cc of the rebuilt build stack vs the test-block stack, and test-block registration),
    against the same numbers for the 2012 templates on the 2012 build/test split."""
    S = setup(); A4 = S["A4"]
    out = []
    for key in a.keys.split(","):
        sta, win = key.split("@")
        res = {}
        for w, tag in ((win, "old"), ("2012-09", "2012")):
            rb, Sb = stack_leg(w, sta, f"build_{tag}", ev=events(w, build_blocks=True), keep_stacks=True)
            rt, St = stack_leg(w, sta, f"test_{tag}", ev=events(w, test_blocks=True), keep_stacks=True)
            regt = {r["family"]: r.get("registration", np.nan) for r in rt if r.get("status") == "ok"}
            for f in sorted(set(Sb) & set(St)):
                c, l = A4.best_cc(Sb[f].astype(float), St[f].astype(float), A4.IDX_NOM - int(round(A4.PRE_TPL * A4.FS)))
                res.setdefault(f, {})[tag] = (c, l, regt.get(f, np.nan))
        for f, d in res.items():
            if "old" in d and "2012" in d:
                out.append(dict(sta=sta, window=win, family=f, cc_build_test_old=d["old"][0], dS_old=d["old"][1],
                                reg_test_old=d["old"][2], cc_build_test_2012=d["2012"][0], reg_test_2012=d["2012"][2]))
    Rb = pd.DataFrame(out); Rb.to_csv(f"{TAB}/g2_rebuild.csv", index=False)
    summ = []
    for (s, w), g in Rb.groupby(["sta", "window"]):
        d = dict(sta=s, window=w, n_fam=len(g), cc_old=g.cc_build_test_old.median(), cc_2012=g.cc_build_test_2012.median(),
                 reg_old=g.reg_test_old.median(), reg_2012=g.reg_test_2012.median())
        d["rebuild_validates"] = bool(d["cc_old"] >= d["cc_2012"] - 0.05 and d["reg_old"] >= d["reg_2012"] - 0.05)
        summ.append(d); G.log(f"[g2 rebuild] {d}")
    pd.DataFrame(summ).to_csv(f"{TAB}/g2_rebuild_summary.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what", choices=["run", "report", "rebuild"])
    ap.add_argument("--workers", type=int, default=12); ap.add_argument("--keys", default="")
    a = ap.parse_args()
    dict(run=run, report=report, rebuild=rebuild)[a.what](a)
