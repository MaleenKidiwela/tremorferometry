#!/usr/bin/env python
"""WA PRODUCTION W9 (member probability), W10 (false share by cut-off), W11 (rates with step checks).

Tremor day classes (per family and day): 'tremor30' / 'tremor60' = >= 1 PNSN tremor epicentre within 30 / 60 km of the
family position that UTC day (catalogs/pnsn_tremor_cascadia_full.csv from 2010-01-06 + the 2009 API pull from 2009-08-06);
before 2009-08-06 no catalogue is public: days inside the WECC ETS windows (2007-01-12..02-02, 2008-04-30..06-02,
2009-04-28..06-03; lfe_array_discovery/results/ab_royer/POOL_RECONSTRUCTION.md sec. 2) count as 'tremor' for every family;
'quiet' = neither.

W9   member   for every data candidate of a family with <= 3-km neighbours: the neighbours' candidates of the same network
              and series within +-2 s; P(member) = softmax(beta * net_mad) over them; beta by maximum likelihood on the
              truth fit rows that matched exactly one member's array detection.  A column; it says WHICH member only if the
              repeat is real (the step-2 finding: equally confident on the scrambled control).
W10  cutoff   per region x label x network x day class x cut-off (0.5 .. 0.9): data rows >= cut vs the TIME-REVERSED control
              (main reference) and the scrambled control; false share = control / data; families covered.
W11  rates    headline = the fixed core network of the era (E3 2015-2026); per family and day: count at P >= 0.75 (label own
              for NW / SW / Arm B, nbr3 for Arm B as a second column) minus the time-reversed control count; Sum P minus the
              control Sum P as a second column; k = truth / net count on the truth days with the same station set, shrunk
              to the region median (log space, pseudo-count 20), day-bootstrap 90 % interval; Olympic Arm A: per <= 3-km
              POSITION GROUP only (candidates of the group's families merged at +-2 s; no per-family Olympic A rate).
              Step checks at every change of the core station set (>= 7 days), node change and era boundary: quiet-day rate
              and rate per PNSN tremor hour, +-21 d and +-45 d, against the same calendar windows in other years (|z| > 2
              flagged).  Inter-ETS rates are a control-subtracted bracket (the SVI withheld-station validation was not run).

  python scripts/wa_prod_post.py tremor | member | cutoff | rates
"""
from __future__ import annotations
import os, sys, json, glob, time, argparse, warnings
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_prod_common as P
W = P.W

OPOST = f"{P.OUT}/post"
WECC = [("2007-01-12", "2007-02-02"), ("2008-04-30", "2008-06-02"), ("2009-04-28", "2009-06-03")]
LAB_REG = {"NW": "own", "SW": "own", "ArmB": "own", "OlyA": "nbr3"}


def log(m):
    P.log(m)


# ------------------------------------------------------------------------------------------------ tremor classes
def tremor(a):
    P.ensure(OPOST)
    F = P.families()
    T = pd.read_csv(f"{P.ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    T["date"] = T.time.astype(str).str[:10]
    p9 = f"{P.ROOT}/lfe_array_discovery/results/ab_royer/tremor/pnsn_tremor_2009.csv"
    if os.path.exists(p9):
        T9 = pd.read_csv(p9, skipinitialspace=True)
        T9.columns = [c.strip() for c in T9.columns]
        T9["date"] = T9.starttime.astype(str).str.strip().str[:10]
        T9 = T9[T9.date < "2010-01-06"][["date", "lat", "lon"]]
        T = pd.concat([T[["date", "lat", "lon"]], T9], ignore_index=True)
    T = T[(T.lat > 45.5) & (T.lat < 49.5) & (T.lon > -125.5) & (T.lon < -121.0)]
    la, lo = T.lat.to_numpy(float), T.lon.to_numpy(float)
    rows = []
    for r in F.itertuples():
        d = W.km(la, lo, r.lat_use, r.lon_use)
        for km_, tag in ((30, "n30"), (60, "n60")):
            c = T[d <= km_].groupby("date").size()
            for dd, v in c.items():
                rows.append((int(r.fid), dd, tag, int(v)))
    X = pd.DataFrame(rows, columns=["fid", "date", "tag", "n"]).pivot_table(index=["fid", "date"], columns="tag", values="n",
                                                                         fill_value=0).reset_index()
    X.to_parquet(f"{OPOST}/tremor_family_day.parquet", index=False)
    log(f"[tremor] {len(T)} epicentres; family-days with tremor <= 60 km {int((X.get('n60', 0) > 0).sum())}")


def day_class(D):
    """adds tremor30 / tremor60 / quiet to rows with fid, date."""
    X = pd.read_parquet(f"{OPOST}/tremor_family_day.parquet")
    k = pd.MultiIndex.from_arrays([D.fid.astype(int), D.date])
    Xi = X.set_index(["fid", "date"])
    n30 = Xi.reindex(k)["n30"].fillna(0).to_numpy() if "n30" in Xi else np.zeros(len(D))
    n60 = Xi.reindex(k)["n60"].fillna(0).to_numpy() if "n60" in Xi else np.zeros(len(D))
    wecc = np.zeros(len(D), bool)
    for a0, b0 in WECC:
        wecc |= (D.date >= a0).to_numpy() & (D.date <= b0).to_numpy()
    pre = (D.date < "2009-08-06").to_numpy()
    D["tremor30"] = np.where(pre, wecc, n30 > 0)
    D["tremor60"] = np.where(pre, wecc, n60 > 0)
    D["n_trem60"] = np.where(pre, np.nan, n60)
    D["quiet"] = ~D.tremor60
    return D


# ------------------------------------------------------------------------------------------------ W9
def member_beta():
    A = pd.read_parquet(f"{P.OUT}/conf/truth_rows.parquet", columns=["date", "net", "series", "fid", "t", "net_mad", "m_own", "row_role"])
    F = P.families()
    nb = P.nbr_map()
    D = A[(A.series == "data") & (A.row_role == "fit")]
    groups = []
    for (net, ds), g in D.groupby(["net", "date"]):
        g = g.sort_values("t")
        t = g.t.to_numpy(); f = g.fid.to_numpy().astype(int); m = g.net_mad.to_numpy(); y = g.m_own.to_numpy()
        for i in np.flatnonzero(y == 1):
            nbs = set(nb.get(f[i], [f[i]]))
            if len(nbs) < 2:
                continue
            j = np.flatnonzero((np.abs(t - t[i]) <= W.MATCH_S) & np.isin(f, list(nbs)))
            if len(j) < 2 or (y[j] == 1).sum() != 1:
                continue
            groups.append((m[j], int(np.flatnonzero(j == i)[0])))
    best = (1.0, -np.inf)
    for beta in np.linspace(0.0, 3.0, 61):
        ll = 0.0
        for mm, ii in groups:
            z = beta * mm; z -= z.max()
            ll += z[ii] - np.log(np.exp(z).sum())
        if ll > best[1]:
            best = (float(beta), float(ll))
    return dict(beta=best[0], loglik=best[1], n_groups=len(groups))


def member(a):
    P.ensure(f"{OPOST}/member")
    mb = member_beta()
    json.dump(mb, open(f"{OPOST}/member_beta.json", "w"), indent=1)
    F = P.families()
    nb = P.nbr_map()
    beta = mb["beta"]
    for run in ("prod", "rich"):
        for p in sorted(glob.glob(f"{P.OUT}/conf/P/{run}/P_*.parquet")):
            o = f"{OPOST}/member/{run}_{os.path.basename(p)}"
            if os.path.exists(o):
                continue
            D = pd.read_parquet(p, columns=["date", "net", "series", "fid", "k", "t", "net_mad", "P_nbr3_WAVE"])
            D = D[D.P_nbr3_WAVE >= 0.3].reset_index(drop=True)      # membership matters only for candidates that may be real
            pm = np.full(len(D), np.nan); nm = np.zeros(len(D), np.int16)
            for _, gi in D.groupby(["date", "net", "series"]).indices.items():
                g = D.iloc[gi]
                order = np.argsort(g.t.to_numpy())
                t = g.t.to_numpy()[order]; f = g.fid.to_numpy().astype(int)[order]; m = g.net_mad.to_numpy()[order]
                lo = np.searchsorted(t, t - W.MATCH_S, "left"); hi = np.searchsorted(t, t + W.MATCH_S, "right")
                for q in range(len(t)):
                    nbs = nb.get(f[q], [f[q]])
                    if len(nbs) < 2:
                        continue
                    j = np.arange(lo[q], hi[q]); j = j[np.isin(f[j], nbs)]
                    z = beta * m[j]; z -= z.max(); e = np.exp(z)
                    pm[gi[order[q]]] = e[j == q].sum() / e.sum(); nm[gi[order[q]]] = len(j)
            pd.DataFrame(dict(date=D.date, k=D.k, p_member=pm, n_members_2s=nm)).to_parquet(o, index=False)
        log(f"[member] {run} done (beta {beta:.2f}, {mb['n_groups']} truth groups)")


# ------------------------------------------------------------------------------------------------ W10
def load_P(run, cols):
    parts = []
    for p in sorted(glob.glob(f"{P.OUT}/conf/P/{run}/P_*.parquet")):
        parts.append(pd.read_parquet(p, columns=cols))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=cols)


def cutoff(a):
    od = f"{P.OUT}/cutoff_check"; P.ensure(od)
    acc = {}; fams = {}; ndays = {}
    for run in ("prod", "rich"):
        for p in sorted(glob.glob(f"{P.OUT}/conf/P/{run}/P_*.parquet")):
            D = pd.read_parquet(p, columns=["date", "net", "series", "fid", "reg", "role", "P_own_WAVE", "P_nbr3_WAVE", "P_any_WAVE"])
            D = D[~D.role.astype(str).str.startswith("sealed")]
            if not len(D):
                continue
            D = day_class(D)
            D["all_days"] = True
            for cls in ("all_days", "tremor30", "tremor60", "quiet"):
                h = D[D[cls].astype(bool)]
                for (net, reg), g in h.groupby(["net", "reg"]):
                    ndays.setdefault((run, net, reg, cls), set()).update(g.date.unique())
                    for lab in ("own", "nbr3", "any"):
                        pc = f"P_{lab}_WAVE"
                        for c in P.CUTS:
                            sel = g[g[pc] >= c]
                            key = (run, net, reg, lab, cls, c)
                            v = acc.setdefault(key, [0, 0, 0])
                            vc = sel.series.value_counts()
                            v[0] += int(vc.get("data", 0)); v[1] += int(vc.get("rev", 0)); v[2] += int(vc.get("c1", 0))
                            fams.setdefault(key, set()).update(sel[sel.series == "data"].fid.unique())
        log(f"[cutoff] {run} done")
    rows = []
    for (run, net, reg, lab, cls, c), (nd, nr, n1) in acc.items():
        rows.append(dict(run=run, net=net, region=reg, label=lab, day_class=cls, cut=c, n_days=len(ndays.get((run, net, reg, cls), ())),
                         data=nd, rev_control=nr, scrambled_control=n1, false_share_rev=nr / nd if nd else np.nan,
                         false_share_scrambled=n1 / nd if nd else np.nan, families_covered=len(fams[(run, net, reg, lab, cls, c)])))
    R = pd.DataFrame(rows).sort_values(["run", "net", "region", "label", "day_class", "cut"])
    R.to_csv(f"{od}/false_share_by_cutoff.csv", index=False)
    h = R[(R.net == "core") & (R.cut == P.P_CUT) & (R.run == "prod") & (R.label == R.region.map(LAB_REG))]
    log("[cutoff] headline (core, label per region, P >= 0.75):\n" + h[["region", "label", "day_class", "data", "rev_control",
        "scrambled_control", "false_share_rev", "false_share_scrambled", "families_covered"]].round(3).to_string())


# ------------------------------------------------------------------------------------------------ W11
def _merged(t):
    t = np.sort(np.asarray(t, float))
    return int((np.diff(t, prepend=-1e18) > W.MATCH_S).sum())


def unit_day_counts(D, pc_fam="P_own_WAVE", pc_grp="P_nbr3_WAVE"):
    """per (unit_kind, unit, date, series): count at P >= 0.75 (Olympic A groups: members merged at +-2 s) and Sum P."""
    F = P.families()
    grp = dict(zip(F.fid, F.pos_group.fillna("")))
    D = D.copy()
    D["olyA"] = D.fid.map(lambda f: grp.get(f, "") != "")
    out = []
    Fm = D[~D.olyA]
    for lab, pc in (("own", "P_own_WAVE"), ("nbr3", "P_nbr3_WAVE")):
        c = Fm[Fm[pc] >= P.P_CUT].groupby(["fid", "date", "series"]).size().rename("n")
        sp = Fm.groupby(["fid", "date", "series"])[pc].sum().rename("sumP")
        X = pd.concat([c, sp], axis=1).fillna(0).reset_index()
        X["unit"] = X.fid.astype(int).astype(str); X["unit_kind"] = "family"; X["label"] = lab
        out.append(X.drop(columns=["fid"]))
    Om = D[D.olyA].copy()
    if len(Om):
        Om["unit"] = Om.fid.map(grp)
        hi = Om[Om.P_nbr3_WAVE >= P.P_CUT]
        c = hi.groupby(["unit", "date", "series"]).t.apply(_merged).rename("n")
        sp = Om.groupby(["unit", "date", "series"]).P_nbr3_WAVE.sum().rename("sumP")
        X = pd.concat([c, sp], axis=1).fillna(0).reset_index()
        X["unit_kind"] = "olyA_group"; X["label"] = "nbr3"
        out.append(X)
    return pd.concat(out, ignore_index=True)


def k_factors():
    """k per (unit, era set) from the truth run: truth events / (data - reversed control) at P >= 0.75 (out-of-fold P),
    on fit days the unit's array run scanned; shrunk to the region median in log space (pseudo-count 20); 100 day-bootstraps."""
    A = pd.read_parquet(f"{P.OUT}/conf/truth_rows.parquet", columns=["date", "net", "series", "fid", "t", "reg", "row_role",
                                                                     "P_own_WAVE", "P_nbr3_WAVE"])
    A = A[(A.row_role == "fit") & A.net.isin(["E1", "E2", "E3"])]
    T = pd.read_parquet(f"{P.OUT}/setup/truth.parquet")
    T = T[T.role == "fit"]
    SC = pd.read_csv(f"{P.OUT}/setup/scan_days.csv")
    F = P.families()
    src = dict(zip(F.fid, F.det_source)); grp = dict(zip(F.fid, F.pos_group.fillna(""))); regm = dict(zip(F.fid, F.reg))
    scanned = set(zip(SC.det_source, SC.date))
    T["unit"] = [grp[f] if grp.get(f, "") else str(int(f)) for f in T.fid]
    T["unit_kind"] = ["olyA_group" if grp.get(f, "") else "family" for f in T.fid]
    TR = T.groupby(["unit_kind", "unit", "date"]).t.apply(lambda x: _merged(x)).rename("n_truth").reset_index()
    rs = np.random.default_rng(7)
    out = []
    for net, An in A.groupby("net"):
        U = unit_day_counts(An)
        U = U[((U.unit_kind == "family") & (U.label == "own")) | (U.unit_kind == "olyA_group")]
        piv = U.pivot_table(index=["unit_kind", "unit", "date"], columns="series", values="n", fill_value=0)
        for c in ("data", "rev"):
            if c not in piv:
                piv[c] = 0
        piv["net_c"] = piv["data"] - piv["rev"]
        pvs = U.pivot_table(index=["unit_kind", "unit", "date"], columns="series", values="sumP", fill_value=0.0)
        for c in ("data", "rev"):
            if c not in pvs:
                pvs[c] = 0.0
        piv["net_s"] = (pvs["data"] - pvs["rev"]).reindex(piv.index).fillna(0.0)
        piv = piv.reset_index()[["unit_kind", "unit", "date", "net_c", "net_s"]]
        days_net = sorted(An.date.unique())
        for (uk, u), g in TR.groupby(["unit_kind", "unit"]):
            f0 = int(u) if uk == "family" else int(F[F.pos_group == u].fid.iloc[0])
            sd = [d for d in days_net if (src[f0], d) in scanned]
            if not sd:
                continue
            tr_d = g.set_index("date").n_truth.reindex(sd).fillna(0).to_numpy()
            pu = piv[(piv.unit_kind == uk) & (piv.unit == u)].set_index("date")
            pv = pu.net_c.reindex(sd).fillna(0).to_numpy()
            ps_ = pu.net_s.reindex(sd).fillna(0).to_numpy()
            n_tr, net_c, net_s = float(tr_d.sum()), float(pv.sum()), float(ps_.sum())
            bs = []
            for _ in range(100):
                ii = rs.integers(0, len(sd), len(sd))
                bs.append((float(tr_d[ii].sum()), float(pv[ii].sum())))
            out.append(dict(net=net, unit_kind=uk, unit=u, reg=regm[f0], n_truth=n_tr, net_count=net_c, n_days=len(sd),
                            net_sumP=net_s, k_sumP_raw=n_tr / net_s if net_s > 0 else np.nan,
                            k_raw=n_tr / net_c if net_c > 0 else np.nan, boot_pairs=bs))
    K = pd.DataFrame(out)
    if len(K):
        K["logk"] = np.log(K.k_raw.where(K.k_raw > 0))
        med = K.groupby(["net", "reg"]).logk.median()
        n0 = 20.0
        def shrunk(t_, c_, m_):
            lk = np.log(t_ / c_) if (c_ > 0 and t_ > 0) else m_
            return float(np.exp((max(c_, 0) * lk + n0 * m_) / (max(c_, 0) + n0)))
        K["k"] = [shrunk(r.n_truth, r.net_count, med.get((r.net, r.reg), np.nan)) for r in K.itertuples()]
        lo, hi = [], []
        for r in K.itertuples():             # the day bootstrap of the SHRUNK estimate (same region median)
            v = [shrunk(t_, c_, med.get((r.net, r.reg), np.nan)) for t_, c_ in r.boot_pairs]
            v = [x for x in v if np.isfinite(x)]
            lo.append(float(np.percentile(v, 5)) if len(v) >= 20 else np.nan); hi.append(float(np.percentile(v, 95)) if len(v) >= 20 else np.nan)
        K["k_boot_lo"] = lo; K["k_boot_hi"] = hi
        K["logks"] = np.log(K.k_sumP_raw.where(K.k_sumP_raw > 0))
        meds = K.groupby(["net", "reg"]).logks.median()
        K["k_sumP"] = [np.exp((max(r.net_sumP, 0) * (r.logks if np.isfinite(r.logks) else meds.get((r.net, r.reg), np.nan))
                               + n0 * meds.get((r.net, r.reg), np.nan)) / (max(r.net_sumP, 0) + n0)) for r in K.itertuples()]
        K = K.drop(columns=["boot_pairs"])
        K["shrink_note"] = "log k shrunk to the (era set, region) median with pseudo-count 20 net detections (k_sumP: 20 units of net Sum P)"
    return K


def rates(a):
    od = f"{P.OUT}/rates"; P.ensure(od)
    K = k_factors()
    K.to_csv(f"{od}/k_factors.csv", index=False)
    parts = []
    for p in sorted(glob.glob(f"{P.OUT}/conf/P/prod/P_*.parquet")):
        D = pd.read_parquet(p, columns=["date", "net", "series", "fid", "t", "P_own_WAVE", "P_nbr3_WAVE"])
        D = D[(D.net == "core") & D.series.isin(["data", "rev", "c1"])]
        if len(D):
            parts.append(unit_day_counts(D))
    U = pd.concat(parts, ignore_index=True)
    R = U.pivot_table(index=["unit_kind", "unit", "label", "date"], columns="series", values=["n", "sumP"], fill_value=0)
    R.columns = [f"{a_}_{b_}" for a_, b_ in R.columns]
    R = R.reset_index()
    for c in ("n_data", "n_rev", "n_c1", "sumP_data", "sumP_rev", "sumP_c1"):
        if c not in R:
            R[c] = 0.0
    R["era"] = R.date.map(P.era_of)
    R["net_count"] = R.n_data - R.n_rev
    R["net_sumP"] = R.sumP_data - R.sumP_rev
    Ki = {(r.net, r.unit_kind, r.unit): (r.k, r.k_boot_lo, r.k_boot_hi, r.k_sumP) for r in K.itertuples()} if len(K) else {}
    kk = [Ki.get((e, uk, u), (np.nan,) * 4) for e, uk, u in zip(R.era, R.unit_kind, R.unit)]
    R["k"] = [x[0] for x in kk]; R["k_lo"] = [x[1] for x in kk]; R["k_hi"] = [x[2] for x in kk]; R["k_sumP"] = [x[3] for x in kk]
    R.loc[~((R.unit_kind == "olyA_group") | (R.label == "own")), ["k", "k_lo", "k_hi", "k_sumP"]] = np.nan   # k only for the headline label
    R["rate_scaled"] = R.k * R.net_count
    R["rate_scaled_lo"] = R.k_lo * R.net_count; R["rate_scaled_hi"] = R.k_hi * R.net_count
    R["sumP_scaled"] = R.k_sumP * R.net_sumP                       # second column: Sum P (minus the reversed control) scaled
    F = P.families()
    first = F[F.pos_group.fillna("") != ""].groupby("pos_group").fid.first().to_dict()
    R["fid"] = [int(u) if uk == "family" else int(first[u]) for uk, u in zip(R.unit_kind, R.unit)]
    R = day_class(R)
    R.to_parquet(f"{od}/rates_daily.parquet", index=False)
    M = R.assign(month=R.date.str[:7]).groupby(["unit_kind", "unit", "label", "month"])[
        ["n_data", "n_rev", "n_c1", "net_count", "sumP_data", "sumP_rev", "net_sumP", "rate_scaled", "rate_scaled_lo",
         "rate_scaled_hi", "sumP_scaled"]].sum(min_count=1).reset_index()
    M.to_csv(f"{od}/rates_monthly.csv", index=False)
    step_checks(R[(R.label == "own") | (R.unit_kind == "olyA_group")], od)
    log(f"[rates] {R.unit.nunique()} units, {len(R)} unit-label-days; k factors {len(K)}")


def core_changes():
    """dates where the loaded core set changes for >= 7 days, + node changes of core stations, + era boundaries."""
    ch = []
    runs = sorted(glob.glob(f"{P.OUT}/mf/prod/run_*.json"))
    prev = None; hist = []
    for p in runs:
        J = json.load(open(p))
        mem = frozenset(J.get("members", {}).get("core", []))
        hist.append((J["date"], mem, J.get("nodes", {})))
    for i in range(1, len(hist)):
        d, m, nd = hist[i]; d0, m0, nd0 = hist[i - 1]
        for s in (m ^ m0):
            nxt = [x for x in hist[i:i + 7]]
            if all((s in x[1]) == (s in m) for x in nxt):
                ch.append((d, s, "joins" if s in m else "leaves"))
        for s in m & m0:
            if nd.get(s) != nd0.get(s):
                ch.append((d, s, f"node {nd0.get(s)} -> {nd.get(s)}"))
    for e, (a0, _) in P.ERAS.items():
        ch.append((a0, "*", f"era {e} starts"))
    return pd.DataFrame(ch, columns=["date", "netsta", "change"]).drop_duplicates()


def step_checks(R, od):
    C = core_changes()
    C.to_csv(f"{od}/core_set_changes.csv", index=False)
    rows = []
    F = P.families()
    R = R.assign(reg=R.fid.map(dict(zip(F.fid, F.reg))))
    for reg_kind, Rg in R.groupby("reg"):
        tot = Rg.groupby("date").agg(net_count=("net_count", "sum"), quiet=("quiet", "mean"),
                                     n_trem=("n_trem60", "sum")).reset_index()
        tot["dt"] = pd.to_datetime(tot.date)
        tot = tot.set_index("dt")
        for r in C.itertuples():
            d = pd.Timestamp(r.date)
            for wdays in (21, 45):
                def stat(center):
                    b = tot[(tot.index >= center - pd.Timedelta(days=wdays)) & (tot.index < center)]
                    a_ = tot[(tot.index > center) & (tot.index <= center + pd.Timedelta(days=wdays))]
                    q = lambda x: x[x.quiet >= 0.5].net_count.mean() if (x.quiet >= 0.5).sum() >= 3 else np.nan
                    th = lambda x: x.net_count.sum() / (x.n_trem.sum() * 2.5 / 60) if x.n_trem.sum() > 24 else np.nan
                    return q(a_) - q(b), th(a_) - th(b)
                dq, dt_ = stat(d)
                ref_q, ref_t = [], []
                for y in range(2007, 2027):
                    if y == d.year:
                        continue
                    try:
                        c2 = d.replace(year=y)
                    except ValueError:
                        continue
                    x, z = stat(c2)
                    ref_q.append(x); ref_t.append(z)
                rq = np.array(ref_q, float); rt = np.array(ref_t, float)
                zq = (dq - np.nanmean(rq)) / np.nanstd(rq) if np.isfinite(rq).sum() >= 4 and np.nanstd(rq) > 0 else np.nan
                zt = (dt_ - np.nanmean(rt)) / np.nanstd(rt) if np.isfinite(rt).sum() >= 4 and np.nanstd(rt) > 0 else np.nan
                rows.append(dict(region=reg_kind, date=r.date, netsta=r.netsta, change=r.change, window_days=wdays,
                                 d_quiet_rate=dq, z_quiet=zq, d_rate_per_tremor_hour=dt_, z_tremor_hour=zt,
                                 flag=bool((np.isfinite(zq) and abs(zq) > 2) or (np.isfinite(zt) and abs(zt) > 2))))
    S = pd.DataFrame(rows)
    S.to_csv(f"{od}/step_checks.csv", index=False)
    log(f"[steps] {len(C)} changes; flagged {int(S.flag.sum()) if len(S) else 0} of {len(S)} checks")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["tremor", "member", "cutoff", "rates"])
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    dict(tremor=tremor, member=member, cutoff=cutoff, rates=rates)[a.what](a)
