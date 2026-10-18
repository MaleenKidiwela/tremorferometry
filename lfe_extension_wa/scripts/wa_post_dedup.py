#!/usr/bin/env python
"""WA post_v2 'counted once' columns: cross-family duplicate FLAGS on the production repeats (user-approved 2026-10-07).

Nothing is deleted and no row is dropped anywhere: duplicates are flag columns (dup_A / dup_B / dup_C with the winner's
family, k and time offset) and extra count columns; the raw counts stay the headline.
Does NOT edit or re-run wa_prod_post_v2.py / wa_prod_common.py / wa_common.py (imported read-only). Reads the production
candidates (conf/P/prod/P_*.parquet: the rows catalogue_v2 and rates_daily_v2 were built from; catalogue_v2 holds only the
data series and no net_cc, so it cannot carry the rule) and the post_v2 outputs; writes NEW files only under
results/wa_prod/post_v2/dedup/.

Rule (the SVI keep_xfam8 rule, svi_science/scripts/s0_s0b_dedup.py): per day x series (data, time-reversed, scrambled control
alike), over all families present in the core-net rows, a candidate is FLAGGED iff a candidate of a station-sharing family
(|A & B| / min(|A|, |B|) >= 0.5 over tpl_library station sets ∩ the era core set) lies within +-2 s with a higher raw net_cc
(ties: the earlier one wins). Candidates at any P count as winners. The winner recorded = the highest-cc such candidate.
  A strict (upper bound): any station-sharing core-net family; known to also catch separate events on busy days.
  B neighbours only:      only a station-sharing family within 10 km (great-circle, setup/families.csv lat_use / lon_use).
  D A + same arrivals:    A with the arrival condition of C (any distance).
  C B + same arrivals:    B, and the two detections predict the same arrival at their best shared station:
                          |(t_i - t_w) - (off_w(s*) - off_i(s*))| <= 0.5 s, s* = shared era-core station with the largest
                          min(template weight), off = tpl_library offset_s of the family's highest-weight node at s*.
Arrival-identity check (all variants, reported): the same mismatch for every flagged pair, against the same statistic for
chance pairs (winner pool shifted by +-60..600 s).
Olympic A position groups already merge their members at +-2 s; a group candidate is flagged only by a family OUTSIDE its group
(members never flag each other); kept member candidates are then merged at +-2 s as before.

  python wa_post_dedup.py keep      # pass over conf/P/prod -> per-candidate flag files, unit-day counts, chance / identity stats
  python wa_post_dedup.py daily     # rates_daily_v2_dedup.parquet (all rows of rates_daily_v2 headline; raw == rates_daily_v2)
  python wa_post_dedup.py kfit      # k refit on flagged-out truth + detections (k_factors_v2 logic, imported)
  python wa_post_dedup.py summary   # rates_summary_v2_dedup.csv + headline table
"""
from __future__ import annotations
import os, sys, json, glob, warnings
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wa_prod_post_v2 as M          # read-only import (no side effects at import)
P = M.P
sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
from numba import njit

O2 = f"{P.OUT}/post_v2"
OD = f"{O2}/dedup"
LOGF = f"{OD}/dedup.log"
W_S = 2.0
NEAR_KM = 10.0
ID_TOL = 0.5
TOLS = (0.25, 0.5, 1.0)
SHIFTS = (60, -60, 120, -120, 300, -300, 600, -600)
RCUTS = M.RCUTS
NPROC = int(os.environ.get("DEDUP_NPROC", "8"))
VARS = ("A", "B", "C", "D")


def vsel(Sm, var):
    """variant -> (share matrix, arrival tolerance; -1 = no arrival condition)."""
    return (Sm["A"] if var in ("A", "D") else Sm["B"]), (ID_TOL if var in ("C", "D") else -1.0)


def log(*a):
    os.makedirs(OD, exist_ok=True)
    P.W.log(LOGF, *a)


# ------------------------------------------------------------------------------------------------ rule pieces
def share_mats():
    """per era: dict A, B (bool, fid-indexed), Doff (best-shared-station offset difference off_b - off_a, NaN = none), km."""
    from s0_s0b_dedup import _wa_share            # same function as the SVI / NW S0b check
    F = P.families()
    T = pd.read_csv(f"{P.OUT}/setup/tpl_library.csv")
    core = json.load(open(f"{P.OUT}/setup/core_sets.json"))
    SA = _wa_share(F, T, core)
    nb = next(iter(SA.values())).shape[0]
    la = np.full(nb, np.nan); lo = np.full(nb, np.nan)
    la[F.fid.astype(int)] = F.lat_use; lo[F.fid.astype(int)] = F.lon_use
    p1, p2 = np.radians(la)[:, None], np.radians(la)[None]
    dl = np.radians(lo)[:, None] - np.radians(lo)[None]
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    KMD = 2 * 6371.0 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))
    near = np.nan_to_num(KMD, nan=1e9) <= NEAR_KM
    best = T.sort_values("w", ascending=False).drop_duplicates(["fid", "netsta"])
    out = {}
    for era, S in SA.items():
        st = list(core[era])
        B = best[best.netsta.isin(st)]
        Wt = np.zeros((nb, len(st))); Of = np.full((nb, len(st)), np.nan)
        si = {s: i for i, s in enumerate(st)}
        Wt[B.fid.astype(int).to_numpy(), B.netsta.map(si).to_numpy()] = B.w.to_numpy()
        Of[B.fid.astype(int).to_numpy(), B.netsta.map(si).to_numpy()] = B.offset_s.to_numpy()
        mm = np.minimum(Wt[:, None, :], Wt[None, :, :])
        ix = mm.argmax(2)
        ok = np.take_along_axis(mm, ix[..., None], 2)[..., 0] > 0
        oa = np.take_along_axis(np.broadcast_to(Of[:, None, :], mm.shape), ix[..., None], 2)[..., 0]
        ob = np.take_along_axis(np.broadcast_to(Of[None, :, :], mm.shape), ix[..., None], 2)[..., 0]
        Doff = np.where(ok, ob - oa, np.nan)
        out[era] = dict(A=S, B=S & near, Doff=Doff, km=KMD)
    return out


def unit_codes(F):
    """fid -> unit code: Olympic A members share their group's code (never flag each other); every other family its own."""
    grp = dict(zip(F.fid.astype(int), F.pos_group.fillna("")))
    gl = sorted({g for g in grp.values() if g})
    nb = int(F.fid.max()) + 2
    u = np.arange(nb, dtype=np.int64)
    for f, g in grp.items():
        if g:
            u[f] = 100000 + gl.index(g)
    return u, grp


@njit(cache=True)
def _keep_u(ts, vs, fs, us, S, Doff, tol, w):
    """time-sorted input; flagged iff j within +-w, different unit, S[f_i, f_j], v_j > v_i or (v_j == v_i and j < i), and
    (tol < 0, or |t_i - t_j - Doff[f_i, f_j]| <= tol). Returns keep and the (sorted) index of the highest-v such j."""
    n = len(ts)
    keep = np.ones(n, np.bool_)
    win = np.full(n, -1, np.int64)
    lo = 0
    for i in range(n):
        while ts[lo] < ts[i] - w:
            lo += 1
        best = -1.0e300
        j = lo
        while j < n and ts[j] <= ts[i] + w:
            if j != i and us[j] != us[i] and S[fs[i], fs[j]]:
                if vs[j] > vs[i] or (vs[j] == vs[i] and j < i):
                    okid = True
                    if tol >= 0.0:
                        d = Doff[fs[i], fs[j]]
                        okid = (d == d) and abs(ts[i] - ts[j] - d) <= tol
                    if okid:
                        keep[i] = False
                        if vs[j] > best:
                            best = vs[j]; win[i] = j
            j += 1
    return keep, win


@njit(cache=True)
def _removed_by_pool(tq, vq, fq, uq, tp, vp, fp, up, S, Doff, tol, w):
    """query vs a (shifted) pool, tp sorted: flagged iff a pool row within +-w of a different unit, S-sharing, higher v (and
    arrival-consistent if tol >= 0). Returns flagged and the identity mismatch to the highest-v such row (NaN if none)."""
    out = np.zeros(len(tq), np.bool_)
    mis = np.full(len(tq), np.nan)
    for i in range(len(tq)):
        lo = np.searchsorted(tp, tq[i] - w)
        hi = np.searchsorted(tp, tq[i] + w, side="right")
        best = -1.0e300
        for j in range(lo, hi):
            if up[j] != uq[i] and S[fq[i], fp[j]] and vp[j] > vq[i]:
                d = Doff[fq[i], fp[j]]
                if tol >= 0.0 and not ((d == d) and abs(tq[i] - tp[j] - d) <= tol):
                    continue
                out[i] = True
                if vp[j] > best:
                    best = vp[j]; mis[i] = tq[i] - tp[j] - d
    return out, mis


def keep_flags(t, v, f, u, S, Doff=None, tol=-1.0):
    o = np.argsort(t, kind="stable")
    if Doff is None:
        Doff = np.zeros(S.shape)
    k, ws = _keep_u(t[o].astype(np.float64), v[o].astype(np.float64), f[o].astype(np.int64), u[f[o]].astype(np.int64), S,
                    Doff, float(tol), W_S)
    keep = np.zeros(len(t), bool); keep[o] = k
    win = np.full(len(t), -1, np.int64); win[o] = np.where(ws >= 0, o[np.clip(ws, 0, None)], -1)
    return keep, win


# ------------------------------------------------------------------------------------------------ pass over P tables
_G = {}


def _init():
    F = P.families()
    _G["F"] = F
    _G["S"] = share_mats()
    _G["U"], _G["grp"] = unit_codes(F)
    _G["span"] = M.catalogue_span()


def _month(path):
    if not _G:
        _init()
    S_, U, grp = _G["S"], _G["U"], _G["grp"]
    cols = ["date", "net", "series", "fid", "k", "t", "net_cc", "P_own_WAVE_PNET", "P_nbr3_WAVE", "reg"]
    D = pd.read_parquet(path, columns=cols)
    D = D[D.net == "core"].drop(columns=["net"]).reset_index(drop=True)
    if not len(D):
        return None
    D["era"] = D.date.map(P.era_of)
    t_all = D.t.to_numpy(float); f_all = D.fid.to_numpy().astype(np.int64); k_all = D.k.to_numpy()
    for var in VARS:
        D[f"dup_{var}"] = False; D[f"win_fid_{var}"] = -1; D[f"win_k_{var}"] = -1
        D[f"win_dt_{var}"] = np.nan; D[f"win_km_{var}"] = np.nan; D[f"arrival_mismatch_{var}"] = np.nan
    for (dt, ser), ix in D.groupby(["date", "series"]).indices.items():
        Sm = S_[P.era_of(dt)]
        t = t_all[ix]; v = D.net_cc.to_numpy(float)[ix]; f = f_all[ix]
        for var in VARS:
            S, tol = vsel(Sm, var)
            kp, wn = keep_flags(t, v, f, U, S, Sm["Doff"], tol)
            fl = ~kp
            gi = ix[fl]; wi = ix[wn[fl]]
            D.loc[gi, f"dup_{var}"] = True
            D.loc[gi, f"win_fid_{var}"] = f_all[wi]; D.loc[gi, f"win_k_{var}"] = k_all[wi]
            D.loc[gi, f"win_dt_{var}"] = t_all[gi] - t_all[wi]
            D.loc[gi, f"win_km_{var}"] = Sm["km"][f_all[gi], f_all[wi]]
            D.loc[gi, f"arrival_mismatch_{var}"] = t_all[gi] - t_all[wi] - Sm["Doff"][f_all[gi], f_all[wi]]
    # ---- per-candidate flag file (every core row kept; det_id = catalogue_v2 det_id for the data series)
    ymd = D.date.str.replace("-", "", regex=False).astype(np.int64).to_numpy()
    D["det_id"] = ymd * 10_000_000 + D.k.to_numpy().astype(np.int64)
    os.makedirs(f"{OD}/flags", exist_ok=True)
    fc = ["det_id", "date", "series", "fid", "k", "t", "net_cc", "era"] + \
         [f"{p}_{v}" for v in VARS for p in ("dup", "win_fid", "win_k", "win_dt", "win_km", "arrival_mismatch")]
    D[fc].to_parquet(f"{OD}/flags/flags_{os.path.basename(path)[2:9]}.parquet", index=False)
    # ---- unit-day counts (same layout as unit_day_counts_v2, headline labels only); 'raw' = no flag applied
    oly = D.fid.map(lambda x: grp.get(int(x), "") != "").to_numpy()
    out = []
    for var in ("raw",) + VARS:
        m = np.ones(len(D), bool) if var == "raw" else ~D[f"dup_{var}"].to_numpy()
        X = D[m]; xo = oly[m]
        Fm = X[~xo]
        parts = [Fm[Fm.P_own_WAVE_PNET >= c].groupby(["fid", "date", "series"]).size().rename(f"n{tg}") for tg, c in RCUTS.items()]
        Y = pd.concat(parts, axis=1).fillna(0).reset_index()
        Y["unit"] = Y.fid.astype(int).astype(str); Y["unit_kind"] = "family"; Y["label"] = "own_pnet"
        Y = Y.drop(columns=["fid"])
        Om = X[xo].copy()
        if len(Om):
            Om["unit"] = Om.fid.map(lambda x: grp[int(x)])
            parts = [Om[Om.P_nbr3_WAVE >= c].groupby(["unit", "date", "series"]).t.apply(M.V1._merged).rename(f"n{tg}")
                     for tg, c in RCUTS.items()]
            Z = pd.concat(parts, axis=1).fillna(0).reset_index()
            Z["unit_kind"] = "olyA_group"; Z["label"] = "nbr3"
            Y = pd.concat([Y, Z], ignore_index=True)
        Y["variant"] = var
        out.append(Y)
    Cnt = pd.concat(out, ignore_index=True)
    # ---- chance coincidence + arrival identity, headline-selected candidates (each series) vs the same series' pool
    D = M.day_class_v2(D, _G["span"])
    reg = D.reg.astype(str).to_numpy()
    hcut = M.cut_headline(reg, D.day_class.to_numpy())
    pl = np.where(reg == "OlyA", D.P_nbr3_WAVE.to_numpy(), D.P_own_WAVE_PNET.to_numpy())
    D["sel"] = np.isfinite(hcut) & (pl >= np.nan_to_num(hcut, nan=9.0))
    ch = []
    for (ser, era), g in D.groupby(["series", "era"]):
        q = g[g.sel]
        if not len(q):
            continue
        Sm = S_[era]
        o = np.argsort(g.t.to_numpy(float))
        tp = g.t.to_numpy(float)[o]; vp = g.net_cc.to_numpy(float)[o]; fp = g.fid.to_numpy().astype(np.int64)[o]; up = U[fp]
        tq = q.t.to_numpy(float); vq = q.net_cc.to_numpy(float); fq = q.fid.to_numpy().astype(np.int64); uq = U[fq]
        for var in VARS:
            S, tol = vsel(Sm, var)
            fl = q[f"dup_{var}"].to_numpy(); mis = q[f"arrival_mismatch_{var}"].to_numpy()
            for (rg, cl), qq in q.groupby([q.reg.astype(str), q.day_class]).indices.items():
                r = dict(series=ser, era=era, variant=var, region=rg, day_class=cl, n_selected=len(qq), n_flagged=int(fl[qq].sum()))
                for tl in TOLS:
                    r[f"flagged_pass{tl}"] = int((np.abs(mis[qq]) <= tl).sum())
                r["chance_pairs"] = 0
                for tl in TOLS:
                    r[f"chance_pass{tl}"] = 0
                for sh in SHIFTS:
                    rem, cm = _removed_by_pool(tq[qq], vq[qq], fq[qq], uq[qq], tp + sh, vp, fp, up, S, Sm["Doff"], tol, W_S)
                    r[f"rm_shift{sh:+d}"] = int(rem.sum())
                    r["chance_pairs"] += int(rem.sum())
                    for tl in TOLS:
                        r[f"chance_pass{tl}"] += int((np.abs(cm[rem]) <= tl).sum())
                ch.append(r)
    return Cnt, pd.DataFrame(ch)


def keep(a):
    from multiprocessing import Pool
    os.makedirs(OD, exist_ok=True)
    fs = sorted(glob.glob(f"{P.OUT}/conf/P/prod/P_*.parquet"))
    with Pool(NPROC, initializer=_init) as pool:
        res = [r for r in pool.imap(_month, fs) if r is not None]
    C = pd.concat([r[0] for r in res], ignore_index=True)
    C.to_parquet(f"{OD}/unit_day_counts_dedup.parquet", index=False)
    H = pd.concat([r[1] for r in res], ignore_index=True)
    key = ["series", "era", "variant", "region", "day_class"]
    H = H.groupby(key).sum(numeric_only=True).reset_index()
    sc = [c for c in H.columns if c.startswith("rm_shift")]
    H["dup_share"] = H.n_flagged / H.n_selected
    H["chance_flagged_mean_shift"] = H[sc].mean(axis=1)
    H["chance_share_of_flagged"] = H.chance_flagged_mean_shift / H.n_flagged.where(H.n_flagged > 0)
    for tl in TOLS:
        H[f"flagged_pass_frac{tl}"] = H[f"flagged_pass{tl}"] / H.n_flagged.where(H.n_flagged > 0)
        H[f"chance_pass_frac{tl}"] = H[f"chance_pass{tl}"] / H.chance_pairs.where(H.chance_pairs > 0)
        cp = H[f"chance_pass_frac{tl}"]
        H[f"same_event_frac_est{tl}"] = (H[f"flagged_pass_frac{tl}"] - cp) / (1 - cp)
    H.to_csv(f"{OD}/chance_identity_dedup.csv", index=False)
    log(f"[keep] {len(fs)} months; unit-day count rows {len(C)}; chance/identity table {len(H)} rows")


# ------------------------------------------------------------------------------------------------ daily table
KEY = ["unit_kind", "unit", "label", "date"]


def daily(a):
    import pyarrow.parquet as pq
    rc = KEY + ["fid", "region", "era", "day_class", "cut_v2", "n_v2_data", "n_v2_rev", "n_v2_c1", "net_count_v2", "net_count_v2_scr",
                "k_v2", "rate_v2"] + [f"n{t}_{s}" for t in RCUTS for s in ("data", "rev", "c1")]
    R = pq.read_table(f"{O2}/rates_daily_v2.parquet", columns=rc,
                      filters=[("label", "in", ["own_pnet", "nbr3"])]).to_pandas()
    R = R[((R.unit_kind == "family") & (R.label == "own_pnet")) | (R.unit_kind == "olyA_group")].reset_index(drop=True)
    C = pd.read_parquet(f"{OD}/unit_day_counts_dedup.parquet")
    vals = [f"n{t}" for t in RCUTS]
    Wd = C.pivot_table(index=KEY, columns=["variant", "series"], values=vals, aggfunc="sum", fill_value=0)
    Wd.columns = [f"{v}_{s}_{var}" for v, var, s in Wd.columns]
    Wd = Wd.reset_index()
    J = R.merge(Wd, on=KEY, how="left", indicator=True)
    extra = int(len(Wd) - (J._merge == "both").sum())
    J = J.drop(columns="_merge")
    for v in vals:
        for s in ("data", "rev", "c1"):
            for var in ("raw",) + VARS:
                c = f"{v}_{s}_{var}"
                J[c] = J[c].fillna(0).astype(np.float64) if c in J else 0.0
    # acceptance: raw rebuild == rates_daily_v2 at every cut and series
    acc = {f"{v}_{s}": int((J[f"{v}_{s}_raw"] != J[f"{v}_{s}"]).sum()) for v in vals for s in ("data", "rev", "c1")}
    acc["count_rows_not_in_rates_daily"] = extra
    tag = np.select([J.cut_v2 == 0.5, J.cut_v2 == 0.75, J.cut_v2 == 0.85], ["050", "075", "085"], "")
    for var in ("raw",) + VARS:
        for s in ("data", "rev", "c1"):
            x = np.full(len(J), np.nan)
            for tg in ("050", "075", "085"):
                m = tag == tg
                x[m] = J[f"n{tg}_{s}_{var}"].to_numpy()[m]
            J[f"n_v2_{s}_{var}"] = x
        J[f"net_count_v2_{var}"] = J[f"n_v2_data_{var}"] - J[f"n_v2_rev_{var}"]
        J[f"net_count_v2_scr_{var}"] = J[f"n_v2_data_{var}"] - J[f"n_v2_c1_{var}"]
    acc["n_v2_raw_vs_rates_daily"] = {s: int(((J[f"n_v2_{s}_raw"] != J[f"n_v2_{s}"]) & J.cut_v2.notna()).sum()) for s in ("data", "rev", "c1")}
    keep_cols = KEY + ["fid", "region", "era", "day_class", "cut_v2", "n_v2_data", "n_v2_rev", "n_v2_c1", "net_count_v2",
                       "net_count_v2_scr", "k_v2", "rate_v2"] + \
                [f"{p}_{var}" for var in VARS for p in ("n_v2_data", "n_v2_rev", "n_v2_c1", "net_count_v2", "net_count_v2_scr")] + \
                [f"n{t}_{s}_{var}" for var in VARS for t in RCUTS for s in ("data", "rev", "c1")]
    J = J[keep_cols]
    J.to_parquet(f"{OD}/rates_daily_v2_dedup.parquet", index=False)
    json.dump(dict(rows=len(J), rows_rates_daily_headline=len(R), acceptance_raw_rebuild=acc), open(f"{OD}/acceptance_dedup.json", "w"), indent=1)
    log("[daily]", len(J), "rows; acceptance", json.dumps(acc))


# ------------------------------------------------------------------------------------------------ k refit
class _PDproxy:
    """pandas proxy handed to wa_prod_post_v2.k_factors_v2: truth_rows come back with only the kept (de-duplicated) candidates
    and truth.parquet with cross-family duplicates removed; everything else is plain pandas."""

    def __init__(self, var):
        self.var = var

    def __getattr__(self, n):
        return getattr(pd, n)

    def read_parquet(self, path, columns=None, **kw):
        if path.endswith("conf/truth_rows.parquet"):
            return _truth_rows_kept(self.var, columns)
        if path.endswith("setup/truth.parquet"):
            return _truth_dedup(self.var)
        return pd.read_parquet(path, columns=columns, **kw)


_TR = {}


def _truth_rows_all():
    if "A" not in _TR:
        cols = ["date", "net", "series", "fid", "t", "reg", "row_role", "P_own_WAVE", "P_nbr3_WAVE", "P_own_WITH_PNET", "net_cc"]
        A = pd.read_parquet(f"{P.OUT}/conf/truth_rows.parquet", columns=cols)
        A = A[A.net.isin(["E1", "E2", "E3"])].reset_index(drop=True)
        F = P.families(); S_ = share_mats(); U, _ = unit_codes(F)
        ks = {v: np.ones(len(A), bool) for v in VARS}
        tA = A.t.to_numpy(float); vA = A.net_cc.to_numpy(float); fA = A.fid.to_numpy().astype(np.int64)
        for (net, dt, ser), ix in A.groupby(["net", "date", "series"]).indices.items():
            Sm = S_[net]
            for v in VARS:
                S, tol = vsel(Sm, v)
                ks[v][ix] = keep_flags(tA[ix], vA[ix], fA[ix], U, S, Sm["Doff"], tol)[0]
        for v in VARS:
            A[f"k{v}"] = ks[v]
        _TR["A"] = A
        log(f"[kfit] truth_rows E1-E3 {len(A):,}: flagged " + ", ".join(f"{v} {int((~ks[v]).sum()):,}" for v in VARS))
    return _TR["A"]


def _truth_rows_kept(var, columns):
    A = _truth_rows_all()
    X = A[A[f"k{var}"]]
    return X[[c for c in (columns or X.columns) if c in X.columns]].copy()


def _truth_dedup(var):
    """truth events (setup/truth.parquet; the file is not changed): within +-2 s, events of linked families (same variant rule:
    station sharing of the event's era, B/C also <= 10 km, C/D also the arrival condition) count once. The winner is the
    family whose own data candidate (truth_rows, any P, any truth net) within +-2 s has the highest net_cc (none: -1e30;
    ties: the earlier event). Olympic A members are never linked to each other (their group merge counts them once)."""
    T = pd.read_parquet(f"{P.OUT}/setup/truth.parquet")
    A = _truth_rows_all()
    d = A[A.series == "data"][["fid", "t", "net_cc"]].sort_values("t")
    td = d.t.to_numpy(float); fd = d.fid.to_numpy().astype(np.int64); vd = d.net_cc.to_numpy(float)
    F = P.families(); S_ = share_mats(); U, _ = unit_codes(F)
    tt = T.t.to_numpy(float); ff = T.fid.to_numpy().astype(np.int64)
    score = np.full(len(T), -1e30)
    by_f = {}
    for f in np.unique(fd):
        m = fd == f
        by_f[f] = (td[m], vd[m])
    for i in range(len(T)):
        x = by_f.get(ff[i])
        if x is None:
            continue
        lo = np.searchsorted(x[0], tt[i] - W_S); hi = np.searchsorted(x[0], tt[i] + W_S, "right")
        if hi > lo:
            score[i] = x[1][lo:hi].max()
    k = np.ones(len(T), bool)
    era = T.date.astype(str).map(P.era_of).fillna("E1").to_numpy()
    for e in np.unique(era):
        ix = np.flatnonzero(era == e)
        S, tol = vsel(S_[e], var)
        k[ix] = keep_flags(tt[ix], score[ix], ff[ix], U, S, S_[e]["Doff"], tol)[0]
    log(f"[kfit] truth events {len(T):,}: variant {var} counts once {int(k.sum()):,} (flagged {int((~k).sum()):,})")
    return T[k].reset_index(drop=True)


def kfit(a):
    os.makedirs(OD, exist_ok=True)
    real = M.pd
    try:
        for var in VARS:
            M.pd = _PDproxy(var)
            K = M.k_factors_v2()
            M.pd = real
            K.to_csv(f"{OD}/rates_k_factors_v2_dedup_{var}.csv", index=False)
            log(f"[kfit] {var}: {len(K)} k rows")
    finally:
        M.pd = real


# ------------------------------------------------------------------------------------------------ summary
def _boot(s, n, rs):
    if not len(s) or n.sum() == 0:
        return np.nan, np.nan, np.nan
    bs = []
    for _ in range(1000):
        ii = rs.integers(0, len(s), len(s))
        bs.append(s[ii].sum() / n[ii].sum())
    return s.sum() / n.sum(), float(np.percentile(bs, 5)), float(np.percentile(bs, 95))


def summary(a):
    J = pd.read_parquet(f"{OD}/rates_daily_v2_dedup.parquet")
    J = J[J.cut_v2.notna()].copy()
    for c in ("n_v2_data", "n_v2_rev", "n_v2_c1", "net_count_v2", "net_count_v2_scr"):
        J[f"{c}_raw"] = J[c]
    # k refit (if present): rate_A / rate_B = k_dedup x de-duplicated net count; k at the row's cut and era (E1/E2 own_pnet NaN)
    for var in VARS:
        p = f"{OD}/rates_k_factors_v2_dedup_{var}.csv"
        if os.path.exists(p):
            K = pd.read_csv(p, dtype={"unit": str, "cut_tag": str})
            K = K.set_index(["net", "cut_tag", "label", "unit_kind", "unit"]).k
            tg = np.where(J.cut_v2 == 0.5, "050", "075")
            kk = K.reindex(pd.MultiIndex.from_arrays([J.era, tg, J.label, J.unit_kind, J.unit.astype(str)])).to_numpy()
            J[f"k_v2_{var}"] = kk
            J[f"rate_v2_{var}"] = kk * J[f"net_count_v2_{var}"]
    J["rate_v2_raw"] = J.rate_v2
    rs = np.random.default_rng(11)
    rows = []
    for (reg, era, cls), g in J.groupby(["region", "era", "day_class"]):
        for var in ("raw",) + VARS:
            rc = f"rate_v2_{var}"
            gg = g.assign(one=1.0, _r=g[rc] if rc in g else np.nan)
            dd = gg.groupby("date").agg(net=(f"net_count_v2_{var}", "sum"), n=("one", "sum"),
                                        rate=("_r", lambda x: x.sum(min_count=1)), k_ok=("_r", lambda x: x.notna().all()))
            n = dd.n.to_numpy()
            ne, nlo, nhi = _boot(dd.net.to_numpy(float), n, rs)
            one = np.ones(len(dd))
            ev, elo, ehi = _boot(dd.net.to_numpy(float), one, rs)          # regional events per day (sum over units)
            kok = bool(dd.k_ok.all())
            re_, rlo, rhi = _boot(dd.rate.to_numpy(float), n, rs) if kok else (np.nan, np.nan, np.nan)
            rows.append(dict(view="headline", region=reg, era=era, day_class=cls, column="rule", count_variant=var,
                             units=g.unit.nunique(), unit_days=int(n.sum()), days=len(dd),
                             data=float(g[f"n_v2_data_{var}"].sum()), rev_control=float(g[f"n_v2_rev_{var}"].sum()),
                             scrambled_control=float(g[f"n_v2_c1_{var}"].sum()),
                             net_per_unit_day=ne, net_boot90_lo=nlo, net_boot90_hi=nhi,
                             regional_events_per_day=ev, regional_boot90_lo=elo, regional_boot90_hi=ehi,
                             k_available=kok, rate_per_unit_day=re_, rate_boot90_lo=rlo, rate_boot90_hi=rhi,
                             k_source=("existing k (non-de-duplicated)" if var == "raw" else "refit on de-duplicated truth + detections") if kok else "",
                             subtracted="time-reversed control", status="headline rule" if var == "raw" else "column (counted once)"))
    S = pd.DataFrame(rows)
    S.to_csv(f"{OD}/rates_summary_v2_dedup.csv", index=False)
    # check: raw rows equal rates_summary_v2 headline 'rule' net per unit-day
    O = pd.read_csv(f"{O2}/rates_summary_v2.csv")
    O = O[(O.view == "headline") & (O.column == "rule")].set_index(["region", "era", "day_class"])
    r = S[S.count_variant == "raw"].set_index(["region", "era", "day_class"])
    jj = r.join(O[["net_per_unit_day", "unit_days"]], rsuffix="_v2", how="inner")
    chk = dict(rows=len(jj), max_abs_diff_net=float((jj.net_per_unit_day - jj.net_per_unit_day_v2).abs().max()),
               unit_days_differ=int((jj.unit_days != jj.unit_days_v2).sum()))
    json.dump(chk, open(f"{OD}/summary_check_vs_rates_summary_v2.json", "w"), indent=1)
    # wide headline table
    Wt = S.pivot_table(index=["region", "era", "day_class"], columns="count_variant",
                       values=["net_per_unit_day", "regional_events_per_day", "rate_per_unit_day"])
    Wt.columns = [f"{a_}_{b_}" for a_, b_ in Wt.columns]
    Wt = Wt.reset_index()
    for var in VARS:
        Wt[f"ratio_{var}_raw"] = Wt[f"net_per_unit_day_{var}"] / Wt.net_per_unit_day_raw
    Wt.to_csv(f"{OD}/headline_table_dedup.csv", index=False)
    log("[summary]", len(S), "rows; check vs rates_summary_v2", json.dumps(chk))
    pd.set_option("display.width", 250)
    print(Wt.round(4).to_string())


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("what", choices=["keep", "daily", "kfit", "summary"])
    a = ap.parse_args()
    globals()[a.what](a)
