#!/usr/bin/env python
"""Quiet-day waveform test, step 2: stacks (test 1), single windows (test 2), timing + earthquakes (test 3).
Reads rows.parquet + work/<date>{.npz,_win.parquet,_row.parquet}; writes tables + figures under quietwave/."""
import os, sys, glob, json
import numpy as np, pandas as pd
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_extension_wa/scripts")
import wa_prod_common as P
W = P.W
import wa_prod_mf as MF

OD = f"{P.OUT}/quietwave"; WK = f"{OD}/work"; FIG = f"{OD}/figures"; TB = f"{OD}/tables"
os.makedirs(FIG, exist_ok=True); os.makedirs(TB, exist_ok=True)
n = W.NT; NV = W.NDAY - n + 1; FS = W.FS
PRE = int(20 * FS); LEN = int(52 * FS)
CUTS = [0.75, 0.85, 0.9]
SER = ["data", "c1", "rev", "random"]
SNAME = {"data": "data", "c1": "scrambled control", "rev": "reversed control", "random": "random times"}
CAT0, CAT1 = "2009-08-06", "2026-05-27"            # PNSN tremor catalogue coverage used by day_class
rng = np.random.default_rng(7)

# ------------------------------------------------------------------ load
S = pd.read_parquet(f"{OD}/rows.parquet")
# random rows inherit the P of the data row they were drawn for (same order per build_rows)
for g in ("quiet", "tremor"):
    d = S[(S.series == "data") & (S.group == g)].sort_values("row_id")
    r = S[(S.series == "random") & (S.group == g)].sort_values("row_id")
    assert len(d) == len(r)
    S.loc[r.index, "P_own_WAVE"] = d.P_own_WAVE.to_numpy()
S["catalogue_era"] = (S.date >= CAT0) & (S.date <= CAT1)
WIN, ROW = [], []
for p in sorted(glob.glob(f"{WK}/*_win.parquet")):
    ds = os.path.basename(p)[:10]
    a = pd.read_parquet(p); a["date"] = ds; WIN.append(a)
    b = pd.read_parquet(f"{WK}/{ds}_row.parquet"); ROW.append(b)
WIN = pd.concat(WIN, ignore_index=True); ROW = pd.concat(ROW, ignore_index=True)
S = S.merge(ROW, on="row_id", how="left")
chk = S[S.series != "random"]
rep = dict(rows=len(S), rows_cut=int(S.n_sta_used.notna().sum()),
           max_abs_netre_minus_netcc=float(np.nanmax(np.abs(chk.net_re - chk.net_cc))),
           frac_netre_match_1e3=float((np.abs(chk.net_re - chk.net_cc) < 1e-3).mean()))
print(rep)
WIN = WIN.merge(S[["row_id", "group", "series", "fid", "P_own_WAVE", "catalogue_era"]], on="row_id")

# tkey per (date, station) and shifts
META = {}
for ds in S.date.unique():
    m = json.load(open(f"{P.OUT}/mf/prod/run_{ds}.json"))
    META[ds] = dict(tk=m["tkeys"], sh={k.split("|")[2]: int(v) for k, v in m["shifts"].items() if k.startswith("core|c1|")},
                    site=m["sites"]["core"])
WIN["tkey"] = [META[d]["tk"].get(s) for d, s in zip(WIN.date, WIN.sta)]
TPL = MF.load_library(rich=False)


def sel(df, group, cut, series, era=None):
    m = (df.group == group) & (df.series == series) & (df.P_own_WAVE >= cut)
    if era is not None:
        m &= df.catalogue_era == era
    return df[m]


# ------------------------------------------------------------------ waveform access (lazy per date)
_WV = {}


def wave(ds, idx):
    if ds not in _WV:
        _WV[ds] = np.load(f"{WK}/{ds}.npz")["wav"]
    return _WV[ds][idx].astype(np.float32)


def stack_metrics(W_, tpl, ch):
    """W_: (N, 3, LEN) normalised windows (non-finite dropped). returns cc_fwd, cc_rev, snr of the linear stack."""
    st = W_.mean(0)
    seg = st[:, PRE:PRE + n].astype(np.float64)
    num = (seg * tpl).sum(-1); den = np.sqrt((seg * seg).sum(-1))
    numr = (seg * tpl[:, ::-1]).sum(-1)
    cf = float(np.mean((num / den)[ch])); cr = float(np.mean((numr / den)[ch]))
    sig = np.sqrt((st[ch][:, PRE:PRE + n] ** 2).mean()); noi = np.sqrt((st[ch][:, :PRE - int(8 * FS)] ** 2).mean())
    return cf, cr, float(sig / noi), st


# ------------------------------------------------------------------ TEST 1: stacks
REFST = {}          # independent long reference: tremor-day data stack per (fid, sta, tkey), P >= 0.85


def refday(ds):
    """tremor-day split: even day-of-month parity days build the independent long reference; tremor-group stacks use the others"""
    return int(ds[-2:]) % 2 == 0


def build_ref():
    d = sel(WIN, "tremor", 0.85, "data")
    d = d[(d.widx >= 0).to_numpy() & np.array([refday(x) for x in d.date], bool)]
    for (f, s, tk), g in d.groupby(["fid", "sta", "tkey"]):
        if len(g) < 5:
            continue
        A = np.stack([wave(ds, i) for ds, i in zip(g.date, g.widx)])
        A = A[np.isfinite(A).all((1, 2))]
        if len(A) >= 5:
            REFST[(f, s, tk)] = A.mean(0)


def outside_cc(st, ref, ch):
    """correlation of the stack OUTSIDE the template window with the tremor-day reference stack:
    pre-S part [-8 s, 0) (holds the P wave) and coda [12, 24) s, channels concatenated."""
    segs = [slice(PRE - int(8 * FS), PRE), slice(PRE + n, PRE + n + int(12 * FS))]
    a = np.concatenate([st[ch][:, s].ravel() for s in segs]); b = np.concatenate([ref[ch][:, s].ravel() for s in segs])
    return float(np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b)))


def test1():
    build_ref()
    out = []
    keepst = {}
    for group in ("quiet", "tremor"):
        for cut in CUTS:
            for ser in SER:
                for era in ([None, True, False] if group == "quiet" else [None]):
                    h = sel(WIN, group, cut, ser, era)
                    h = h[h.widx >= 0]
                    if group == "tremor" and len(h):
                        h = h[~np.array([refday(x) for x in h.date], bool)]
                    nev = h.groupby("fid").row_id.nunique()
                    fams = nev[nev >= 5].index
                    for f in fams:
                        hf = h[h.fid == f]
                        for (s, tk), g in hf.groupby(["sta", "tkey"]):
                            if len(g) < 5 or (f, s, tk) not in TPL:
                                continue
                            A = np.stack([wave(ds, i) for ds, i in zip(g.date, g.widx)])
                            ok = np.isfinite(A).all((1, 2))
                            A = A[ok]
                            if len(A) < 5:
                                continue
                            t = TPL[(f, s, tk)]; ch = t["live"] & (np.abs(A).sum((0, 2)) > 0)
                            if not ch.any():
                                continue
                            cf, cr, snr, st = stack_metrics(A, t["tpl"], ch)
                            # fixed-size stack (5 random events) for a fair comparison across series
                            sub = A[rng.choice(len(A), 5, replace=False)]
                            cf5, cr5, snr5, _ = stack_metrics(sub, t["tpl"], ch)
                            oc = outside_cc(st, REFST[(f, s, tk)], ch) if (f, s, tk) in REFST else np.nan
                            out.append(dict(group=group, cut=cut, series=ser, era={None: "all", True: "catalogue", False: "no_catalogue"}[era],
                                            fid=f, sta=s, tkey=tk, n_ev=len(A), n_spike=int((~ok).sum()), coef=float(g.coef.mean()),
                                            stack_cc_fwd=cf, stack_cc_rev=cr, stack_snr=snr, cc_fwd_n5=cf5, cc_rev_n5=cr5,
                                            snr_n5=snr5, outside_cc_ref=oc))
                            if era is None:
                                keepst[(group, cut, ser, f, s, tk)] = (st, len(A))
    T = pd.DataFrame(out)
    T.to_csv(f"{TB}/t1_stacks_per_family_station.csv", index=False)
    # summary: median over (family, station); also per family the coef-weighted mean, then median over families
    rows = []
    for k, g in T.groupby(["group", "era", "cut", "series"]):
        fam = g.groupby("fid").apply(lambda x: pd.Series(dict(
            cc=np.average(x.stack_cc_fwd, weights=x.coef), snr=np.average(x.stack_snr, weights=x.coef),
            cc5=np.average(x.cc_fwd_n5, weights=x.coef), snr5=np.average(x.snr_n5, weights=x.coef),
            ccr=np.average(x.stack_cc_rev, weights=x.coef))))
        rows.append(dict(zip(["group", "era", "cut", "series"], k), n_families=g.fid.nunique(), n_fam_sta=len(g),
                         median_events=float(g.groupby("fid").n_ev.max().median()),
                         stack_cc_fwd_med_famsta=g.stack_cc_fwd.median(), stack_snr_med_famsta=g.stack_snr.median(),
                         stack_cc_fwd_med_fam_netw=fam.cc.median(), stack_snr_med_fam_netw=fam.snr.median(),
                         stack_cc_rev_med_fam_netw=fam.ccr.median(),
                         n5_cc_fwd_med_fam_netw=fam.cc5.median(), n5_snr_med_fam_netw=fam.snr5.median(),
                         outside_cc_ref_med=g.outside_cc_ref.median(), n_outside=int(g.outside_cc_ref.notna().sum())))
    Sm = pd.DataFrame(rows)
    Sm.to_csv(f"{TB}/t1_stacks_summary.csv", index=False)
    return T, Sm, keepst


# ------------------------------------------------------------------ TEST 2: single windows
def test2():
    rows = []
    # network mean per candidate (plain mean and production coef-weighted)
    WIN["fin"] = np.isfinite(WIN.cc_fwd)
    w = WIN[WIN.fin]
    net = w.groupby("row_id").apply(lambda x: pd.Series(dict(
        net_fwd_mean=x.cc_fwd.mean(), net_fwd_w=np.average(x.cc_fwd, weights=x.coef),
        net_rev_w=np.average(x.cc_rev.fillna(0), weights=x.coef),
        frac_sta_cc02=(x.cc_fwd >= 0.2).mean(), sig_med=x.rms_sig.median(), snrw_med=(x.rms_sig / x.rms_pre).median(),
        envpk_max=x.env_pk.max(), kurt_max=x["kurt"].max(),
        spike=bool((~np.isfinite(x.env_pk) & (x.widx >= 0)).any() or (x.rms_sig > 50).any()))))
    net = net.reset_index().merge(S[["row_id", "group", "series", "P_own_WAVE", "catalogue_era", "hfr_rel", "date", "fid", "t", "i"]], on="row_id")
    net.to_parquet(f"{TB}/t2_network_per_candidate.parquet", index=False)
    q = lambda v, p: float(np.nanquantile(v, p)) if len(v) else np.nan
    for group in ("quiet", "tremor"):
        for era in (["all", "catalogue", "no_catalogue"] if group == "quiet" else ["all"]):
            for cut in CUTS:
                for ser in SER:
                    m = (net.group == group) & (net.series == ser) & (net.P_own_WAVE >= cut)
                    mw = (WIN.group == group) & (WIN.series == ser) & (WIN.P_own_WAVE >= cut) & WIN.fin
                    if era != "all":
                        m &= net.catalogue_era == (era == "catalogue"); mw &= WIN.catalogue_era == (era == "catalogue")
                    a = net[m]; b = WIN[mw]
                    rows.append(dict(group=group, era=era, cut=cut, series=ser, n_cand=len(a), n_win=len(b),
                                     sta_cc_fwd_p25=q(b.cc_fwd, .25), sta_cc_fwd_med=q(b.cc_fwd, .5), sta_cc_fwd_p75=q(b.cc_fwd, .75),
                                     sta_frac_cc_ge02=float((b.cc_fwd >= 0.2).mean()) if len(b) else np.nan,
                                     net_fwd_mean_med=q(a.net_fwd_mean, .5), net_fwd_w_med=q(a.net_fwd_w, .5),
                                     net_fwd_w_p10=q(a.net_fwd_w, .1), net_fwd_w_p90=q(a.net_fwd_w, .9),
                                     net_rev_w_med=q(a.net_rev_w, .5),
                                     win_energy_med=q(a.sig_med, .5), win_snr_med=q(a.snrw_med, .5),
                                     frac_win_snr_gt12=float((a.snrw_med > 1.2).mean()) if len(a) else np.nan,
                                     envpk_max_med=q(a.envpk_max, .5), frac_spike=float(a.spike.mean()) if len(a) else np.nan,
                                     hfr_rel_med=q(a.hfr_rel, .5)))
    T2 = pd.DataFrame(rows)
    T2.to_csv(f"{TB}/t2_single_windows.csv", index=False)
    return net, T2


# ------------------------------------------------------------------ TEST 3: timing + earthquakes
def test3(net):
    A = pd.read_parquet(f"{OD}/core_data_all_regions.parquet")
    A = A[A.date.isin(set(S.date))]
    SETS = {}
    for lab, m in [("same_fam_P50", None), ("anyfam_P50", A.P_own_WAVE >= 0.5), ("NWfam_P50", (A.P_own_WAVE >= 0.5) & (A.reg == "NW")),
                   ("NWfam_P85", (A.P_own_WAVE >= 0.85) & (A.reg == "NW")), ("NWfam_all", A.reg == "NW")]:
        if m is None:
            continue
        SETS[lab] = {ds: np.sort(g.t.to_numpy()) for ds, g in A[m].groupby("date")}
    SAME = {(ds, f): np.sort(g.t.to_numpy()) for (ds, f), g in A[A.P_own_WAVE >= 0.5].groupby(["date", "fid"])}

    def near(arr, t, tol):
        if arr is None or len(arr) == 0:
            return False
        k = np.searchsorted(arr, t)
        return any(abs(arr[j] - t) <= tol for j in (k - 1, k) if 0 <= j < len(arr))
    res = []
    W2 = WIN[np.isfinite(WIN.cc_fwd)]
    site_score = W2.assign(c=W2.coef * W2.cc_fwd).groupby(["row_id", "site"]).c.sum().reset_index()
    SSD = {k: g for k, g in site_score.groupby("row_id")}
    for r in S.itertuples():
        if r.series not in ("c1", "rev", "random") or not (r.n_sta_used > 0):
            continue
        D0 = pd.Timestamp(r.date).timestamp()
        out = dict(row_id=r.row_id, group=r.group, series=r.series, P=r.P_own_WAVE, era=r.catalogue_era)
        if r.series in ("rev", "random"):
            for tol in (2, 30):
                out[f"same_fam_P50_{tol}"] = near(SAME.get((r.date, r.fid)), r.t, tol)
                for lab, Sd in SETS.items():
                    out[f"{lab}_{tol}"] = near(Sd.get(r.date), r.t, tol)
        if r.series in ("c1", "random"):
            # the real time each site's window came from: origin (i - shift_site) % NV (for random rows the same
            # day's shifts are applied to the random origin = the chance baseline of a scrambled sum)
            ss = SSD.get(r.row_id)
            if ss is not None and len(ss):
                sh = META[r.date]["sh"]
                rt = {q: D0 + ((int(r.i) - sh[q]) % NV) / FS for q in ss.site}
                dom = ss.site.iloc[int(np.argmax(ss.c.to_numpy()))]
                out["c1_dom_same_fam_P50_2"] = near(SAME.get((r.date, r.fid)), rt[dom], 2)
                for lab in ("anyfam_P50", "NWfam_P50", "NWfam_P85", "NWfam_all"):
                    out[f"c1_dom_{lab}_2"] = near(SETS[lab].get(r.date), rt[dom], 2)
                    out[f"c1_nsites_{lab}_2"] = sum(near(SETS[lab].get(r.date), rt[q], 2) for q in rt)
                out["c1_dom_share"] = float(ss.c.max() / ss.c.clip(lower=0).sum()) if ss.c.clip(lower=0).sum() > 0 else np.nan
        res.append(out)
    T3 = pd.DataFrame(res)
    T3.to_parquet(f"{TB}/t3_timing_per_candidate.parquet", index=False)
    # earthquakes: McCrory 2012 ANSS hypocentres (1975-2009, M >= 1); NW box; origin within [t-60 s, t+30 s]
    E = pd.concat([pd.read_csv(f"{P.ROOT}/data/slab_models/mccrory2012/ANSS_hypos_{x}_slab_GCS_WGS84.txt") for x in ("above", "below")])
    E = E[(E.LAT > 46.0) & (E.LAT < 49.5) & (E.LONG > -125.0) & (E.LONG < -121.0)]
    te = np.sort((pd.to_datetime(dict(year=E.YEAR, month=E.MONTH, day=E.DAY, hour=E.HOUR, minute=E.MIN))
                  .values.astype("datetime64[ns]").astype(np.int64) / 1e9 + E.SEC.to_numpy()))
    eqmax = te.max()
    Sx = S[(S.t <= eqmax) & (S.date >= "2007-01-01") & (S.n_sta_used > 0)].copy()
    k0 = np.searchsorted(te, Sx.t.to_numpy() - 60); k1 = np.searchsorted(te, Sx.t.to_numpy() + 30)
    Sx["eq_near"] = k1 > k0
    eq = Sx.groupby(["group", "series"]).apply(lambda g: pd.Series({f"P>={c}": f"{g[g.P_own_WAVE >= c].eq_near.mean():.3f} (n={int((g.P_own_WAVE >= c).sum())})" for c in CUTS}))
    eq.to_csv(f"{TB}/t3_earthquake_coincidence_2007_2009.csv")
    # summaries
    rows = []
    for (group, ser), g in T3.groupby(["group", "series"]):
        for era in (["all", "catalogue", "no_catalogue"] if group == "quiet" else ["all"]):
            for cut in CUTS:
                h = g[g.P >= cut]
                if era != "all":
                    h = h[h.era == (era == "catalogue")]
                d = dict(group=group, era=era, series=ser, cut=cut, n=len(h))
                for c in h.columns:
                    if c.endswith("_2") or c.endswith("_30"):
                        d[c] = float(h[c].astype(float).mean()) if len(h) else np.nan
                if "c1_nsites_NWfam_P50_2" in h:
                    d["c1_ge2sites_NWfam_P50"] = float((h.c1_nsites_NWfam_P50_2 >= 2).mean()) if len(h) else np.nan
                    d["c1_ge2sites_NWfam_all"] = float((h.c1_nsites_NWfam_all_2 >= 2).mean()) if len(h) else np.nan
                    d["c1_dom_share_med"] = float(h.c1_dom_share.median()) if len(h) else np.nan
                rows.append(d)
    T3s = pd.DataFrame(rows)
    T3s.to_csv(f"{TB}/t3_timing_summary.csv", index=False)
    # activity on the day: data detections P >= 0.85 of the same family that day
    act = A[(A.P_own_WAVE >= 0.85)].groupby(["date", "fid"]).size()
    S["fam_day_P85"] = [int(act.get((d, f), 0)) for d, f in zip(S.date, S.fid)]
    actn = A[(A.P_own_WAVE >= 0.85) & (A.reg == "NW")].groupby("date").size()
    S["nw_day_P85"] = [int(actn.get(d, 0)) for d in S.date]
    ad = []
    for (group, ser), g in S[S.n_sta_used > 0].groupby(["group", "series"]):
        for cut in CUTS:
            h = g[g.P_own_WAVE >= cut]
            ad.append(dict(group=group, series=ser, cut=cut, n=len(h), med_nw_day_P85=h.nw_day_P85.median(),
                           frac_day_nw_P85_ge10=(h.nw_day_P85 >= 10).mean(), frac_famday_has_other_P85=(h.fam_day_P85 - (h.series == "data") * (h.P_own_WAVE >= 0.85) >= 1).mean()))
    pd.DataFrame(ad).to_csv(f"{TB}/t3_day_activity.csv", index=False)
    return T3, T3s, eq


# ------------------------------------------------------------------ figures
def figures(T, keepst, net):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    COL = {"data": "#2a78d6", "c1": "#eb6834", "rev": "#1baf7a", "random": "#8a8984"}
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    # example families: quiet P >= 0.75 with stacks in data and at least one control
    q = T[(T.group == "quiet") & (T.cut == 0.75) & (T.era == "all")]
    fam_ser = q.groupby("fid").series.apply(set)
    cand = [f for f, s in fam_ser.items() if "data" in s and ("c1" in s or "rev" in s)]
    cand = sorted(cand, key=lambda f: -q[(q.fid == f) & (q.series != "data") & (q.series != "random")].n_ev.max())[:4]
    for f in cand:
        stas = q[(q.fid == f) & (q.series == "data")].sort_values("coef", ascending=False).drop_duplicates("sta").head(4)
        fig, axs = plt.subplots(len(stas), 1, figsize=(9, 2.8 * len(stas) + 0.8), sharex=True, squeeze=False)
        tt = (np.arange(LEN) - PRE) / FS
        for ax, r in zip(axs[:, 0], stas.itertuples()):
            t = TPL[(f, r.sta, r.tkey)]; c = int(np.flatnonzero(t["live"])[0])
            tp = t["tpl"][c] / np.abs(t["tpl"][c]).max()
            ax.plot(np.arange(n) / FS, 2.5 * tp, color="#0b0b0b", lw=1.0, label="forward template")
            for off, ser in zip((-8, -16, -24, -32), ("data", "c1", "rev", "random")):
                k = ("quiet", 0.75, ser, f, r.sta, r.tkey)
                if k in keepst:
                    st, nn = keepst[k]
                    y = st[c] / (np.sqrt((st[c][:PRE - int(8 * FS)] ** 2).mean()) + 1e-12)   # units of pre-window noise rms
                    ax.plot(tt, np.clip(y, -6, 6) + off, color=COL[ser], lw=0.8, label=f"{SNAME[ser]} stack (n={nn})")
            ax.axvspan(0, n / FS, color="#e9e8e4", zorder=0)
            ax.set_ylabel(f"{r.sta.split('.')[1]} ch{c}"); ax.set_yticks([])
        axs[0, 0].legend(loc="upper left", fontsize=7, ncol=3, frameon=False, bbox_to_anchor=(0, 1.45))
        axs[-1, 0].set_xlabel("time from template start at the station (s); grey = 12 s template window")
        fig.suptitle(f"Family {f}: quiet days, P >= 0.75; stacks in units of their pre-window noise (clipped at 6)", y=0.995, fontsize=10)
        fig.tight_layout(); fig.savefig(f"{FIG}/example_family_{f}.png", dpi=130); plt.close(fig)
    # single-window distributions
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.4))
    for ax, cut in zip(axs, CUTS):
        for ser in SER:
            a = net[(net.group == "quiet") & (net.series == ser) & (net.P_own_WAVE >= cut)].net_fwd_w.dropna()
            if len(a) >= 5:
                ax.hist(a, bins=np.linspace(-0.1, 0.4, 51), histtype="step", density=True, color=COL[ser], lw=1.4,
                        label=f"{SNAME[ser]} (n={len(a)})")
        ax.set_title(f"quiet days, P >= {cut}"); ax.set_xlabel("network cc with FORWARD template (weighted)")
        ax.legend(fontsize=7, frameon=False)
    axs[0].set_ylabel("density")
    fig.tight_layout(); fig.savefig(f"{FIG}/single_window_netcc_quiet.png", dpi=130); plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.4))
    for ax, grp in zip(axs, ("quiet", "tremor")):
        for ser in SER:
            a = net[(net.group == grp) & (net.series == ser) & (net.P_own_WAVE >= 0.85)].snrw_med.dropna()
            if len(a) >= 5:
                ax.hist(np.clip(a, 0.5, 3), bins=np.linspace(0.5, 3, 51), histtype="step", density=True, color=COL[ser], lw=1.4,
                        label=f"{SNAME[ser]} (n={len(a)})")
        ax.set_title(f"{grp} days, P >= 0.85"); ax.set_xlabel("window energy / pre-window energy (station median, rms ratio)")
        ax.legend(fontsize=7, frameon=False)
    fig.tight_layout(); fig.savefig(f"{FIG}/single_window_energy_ratio.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    T, Sm, keepst = test1()
    print(Sm.to_string())
    net, T2 = test2()
    print(T2.to_string())
    T3, T3s, eq = test3(net)
    print(T3s.to_string()); print(eq.to_string())
    figures(T, keepst, net)
    json.dump(rep, open(f"{TB}/replication_check.json", "w"), indent=1)
