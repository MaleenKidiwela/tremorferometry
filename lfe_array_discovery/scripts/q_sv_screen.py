#!/usr/bin/env python
"""Apply the FIVE-CLAUSE Washington screen (scripts/n_screen.py, WA_2008_RESULTS.md sec.4),
clause by clause and exactly as coded, to the 2005 POLARIS bed where the ground truth is known.

Populations (results/repro/classify15.csv, the K_REPRO3 15-station run's 1,460 grown templates):
  KEPT    177 templates matching one of the 58 rediscovered Bostock development families
  HIDDEN   64 templates matching one of the 18 located sealed families
  CAND     18 candidate-new templates
  NOISE  1201 unlocatable / unmatched templates

Clause definitions are lifted from n_screen.py (FA_MAX 0.1, FA_RATIO 3, DZ_MAX 8, MIN_DAYS 3,
STAB_KM 10, N_BASE 200) and applied unchanged; only the bed adapter differs.

  python scripts/q_sv_screen.py
"""
from __future__ import annotations
import json, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import k_common as K                                                        # noqa: E402

FA_MAX, FA_RATIO, DZ_MAX, MIN_DAYS, STAB_KM, N_BASE = 0.1, 3.0, 8.0, 3, 10.0, 200
OUT = "results/q_screenval"
LATS = np.arange(47.9, 49.31, 0.03)
LONS = np.arange(-124.9, -122.79, 0.04)
DEPS = np.arange(18.0, 56.1, 2.0)


def hyp_km(la1, lo1, dp1, la2, lo2, dp2):
    import s_common as S
    d = S.epi_km(np.array([la1]), np.array([lo1]), la2, lo2)[0]
    return float(np.sqrt(d ** 2 + (dp1 - dp2) ** 2))


def is_edge(la, lo, dp, t=1e-6):
    return bool(abs(la - LATS[0]) < t or abs(la - LATS[-1]) < t
                or abs(lo - LONS[0]) < t or abs(lo - LONS[-1]) < t
                or abs(dp - DEPS[0]) < t or abs(dp - DEPS[-1]) < t)


def main():
    os.makedirs(OUT, exist_ok=True)
    C = pd.read_csv("results/repro/classify15.csv")
    D = pd.read_parquet("results/repro/grow15/det_it0.parquet")
    tdays = sorted(set(D.date)); n_t = len(tdays)

    # ---- quiet-day detections (clause a) --------------------------------------------------
    qp = f"{OUT}/det_quiet15.parquet"
    Q = pd.read_parquet(qp) if os.path.exists(qp) else None
    if Q is not None:
        QD = pd.read_csv(f"{OUT}/quiet_day_stations.csv")
        qdays = sorted(QD[QD.used].date); n_q = len(qdays)
    else:
        qdays, n_q = [], 0

    # ---- the 2005 tremor-hour reference (clause b) ----------------------------------------
    B = K.bostock(2005).copy()
    B["date"] = B.time.dt.strftime("%Y-%m-%d")
    B = B[B.date.isin(set(tdays))]
    B["ep"] = (B.t_origin.to_numpy(np.float64) // 3600).astype(np.int64)
    # hour-epoch -> per-family event counts, so a family can be left out of its own reference
    cnt_all = B.groupby("ep").size()
    cnt_fam = B.groupby(["ep", "bkey"]).size()
    EPS = np.array(sorted({int(pd.Timestamp(f"{d} 00:00:00").timestamp()) // 3600 + h
                           for d in tdays for h in range(24)}), dtype=np.int64)
    base_all = cnt_all.reindex(EPS).fillna(0).to_numpy()
    # the LOUD-hour variant, matched to the WA 2008 tremor-hour base rate of 0.24
    q_thr = float(np.quantile(base_all, 1 - 0.24))

    def hourset(drop_bkey=None):
        c = base_all.copy()
        if drop_bkey is not None and drop_bkey >= 0:
            sub = cnt_fam.xs(drop_bkey, level="bkey") if drop_bkey in \
                cnt_fam.index.get_level_values("bkey") else None
            if sub is not None:
                c = c - pd.Series(sub).reindex(EPS).fillna(0).to_numpy()
        return set(EPS[c >= 1]), set(EPS[c > q_thr])

    ANY_ALL, LOUD_ALL = hourset(None)
    rng = np.random.default_rng(17)
    BYG = {int(g): gg for g, gg in D.groupby("grp")}
    nq = Q.groupby("grp").size() if Q is not None and len(Q) else pd.Series(dtype=int)

    def cotime(d, HS):
        tv = d.t.to_numpy(np.float64)
        dsec = np.array([pd.Timestamp(x).timestamp() for x in d.date], dtype=np.float64)

        def frac(ts):
            return float(np.mean([int(x) in HS for x in np.floor_divide(
                ts.astype(np.int64), 3600)]))
        obs = frac(tv)
        bs = np.array([frac(dsec + rng.random(len(dsec)) * 86400.0) for _ in range(N_BASE)])
        return obs, float(bs.mean()), float((bs >= obs).mean())

    rows = []
    for r in C.itertuples():
        g = int(r.grp)
        d = BYG.get(g, D.iloc[0:0])
        cls = str(r.cls)
        fam = int(r.bkey) if r.bkey >= 0 else (int(r.hid_bkey) if r.hid_bkey >= 0 else -1)
        rec = dict(grp=g, cls=cls, bkey=int(r.bkey), hid_bkey=int(r.hid_bkey), fam=fam,
                   n_sta=int(r.n_sta), n_det=len(d),
                   n_days=int(d.date.nunique()) if len(d) else 0,
                   n_det_quiet=int(nq.get(g, 0)) if Q is not None else -1)
        # --- (a) -----------------------------------------------------------------------
        if Q is not None:
            rec["fa_quiet_per_tpl_day"] = rec["n_det_quiet"] / max(n_q, 1)
            rec["rate_tremor_per_tpl_day"] = rec["n_det"] / max(n_t, 1)
            rec["c_a"] = bool(rec["fa_quiet_per_tpl_day"] <= FA_MAX
                              and rec["rate_tremor_per_tpl_day"]
                              >= FA_RATIO * rec["fa_quiet_per_tpl_day"] and rec["n_det"] > 0)
        else:
            rec["fa_quiet_per_tpl_day"] = np.nan
            rec["rate_tremor_per_tpl_day"] = rec["n_det"] / max(n_t, 1)
            rec["c_a"] = None
        # --- (b) -----------------------------------------------------------------------
        if len(d):
            HS_any, HS_loud = hourset(fam) if fam >= 0 else (ANY_ALL, LOUD_ALL)
            o, b, p = cotime(d, HS_any)
            rec.update(co_any=o, co_any_base=b, co_any_p=p, c_b=bool(o > b))
            o2, b2, p2 = cotime(d, HS_loud)
            rec.update(co_loud=o2, co_loud_base=b2, co_loud_p=p2,
                       c_b_loud=bool(o2 > b2), c_b_loud_p=bool(p2 < 0.05))
        else:
            rec.update(co_any=np.nan, co_any_base=np.nan, co_any_p=np.nan, c_b=False,
                       co_loud=np.nan, co_loud_base=np.nan, co_loud_p=np.nan,
                       c_b_loud=False, c_b_loud_p=False)
        # --- (c) -----------------------------------------------------------------------
        loc = bool(r.located)
        rec["c_c1"] = loc
        rec["c_c2"] = bool(r.ctrl_pass) if loc else False
        rec["c_c3"] = (not is_edge(r.lat, r.lon, r.dep)) if loc else False
        rec["c_c4"] = bool(np.isfinite(r.dz) and abs(float(r.dz)) <= DZ_MAX) if loc else False
        rec.update(lat=r.lat, lon=r.lon, dep=r.dep, dz=r.dz, rms=r.rms,
                   rms_shuf_p05=r.rms_shuf_p05)
        # --- (d) -----------------------------------------------------------------------
        rec["c_d"] = bool(rec["n_days"] >= MIN_DAYS)
        rows.append(rec)
    V = pd.DataFrame(rows)

    # --- (e) -------------------------------------------------------------------------------
    for pref in ("h1", "h2"):
        p = f"{OUT}/reloc_{pref}.csv"
        if not os.path.exists(p):
            V[f"{pref}_km"] = np.nan
            continue
        R = pd.read_csv(p).set_index("grp")
        km, la, lo, dp, ed = [], [], [], [], []
        for r in V.itertuples():
            if r.grp in R.index and bool(R.located.get(r.grp, False)) and np.isfinite(r.lat):
                q = R.loc[r.grp]
                km.append(hyp_km(r.lat, r.lon, r.dep, float(q.lat), float(q.lon), float(q.dep)))
                la.append(float(q.lat)); lo.append(float(q.lon)); dp.append(float(q.dep))
                ed.append(bool(q.edge))
            else:
                km.append(np.nan); la.append(np.nan); lo.append(np.nan); dp.append(np.nan)
                ed.append(None)
        V[f"{pref}_km"] = km
        V[f"{pref}_lat"] = la; V[f"{pref}_lon"] = lo; V[f"{pref}_dep"] = dp
        V[f"{pref}_edge"] = ed
    V["c_e"] = [bool(np.isfinite(x) and x <= STAB_KM) for x in V.h2_km]
    if "h1_lat" in V.columns:
        V["h1_h2_km"] = [hyp_km(r.h1_lat, r.h1_lon, r.h1_dep, r.h2_lat, r.h2_lon, r.h2_dep)
                         if np.isfinite(getattr(r, "h1_lat", np.nan))
                         and np.isfinite(getattr(r, "h2_lat", np.nan)) else np.nan
                         for r in V.itertuples()]
        V["c_e_hh"] = [bool(np.isfinite(x) and x <= STAB_KM) for x in V.h1_h2_km]

    # --- (a) data-driven reformulation: day-label shuffle -----------------------------------
    if Q is not None:
        alld = list(tdays) + list(qdays)
        lab = np.array([1] * n_t + [0] * n_q)
        cnt_t = {g: gg.groupby("date").size() for g, gg in D.groupby("grp")}
        cnt_q = {g: gg.groupby("date").size() for g, gg in Q.groupby("grp")} if len(Q) else {}
        rs = np.random.default_rng(101)
        nd_all = n_t + n_q
        PM = np.array([rs.permutation(nd_all) for _ in range(2000)])   # shared label shuffles
        ps = []
        for g in V.grp:
            c = np.array([float(cnt_t.get(g, pd.Series(dtype=int)).get(d, 0)) for d in tdays]
                         + [float(cnt_q.get(g, pd.Series(dtype=int)).get(d, 0)) for d in qdays])
            if c.sum() == 0:
                ps.append(1.0); continue
            obs = c[lab == 1].mean() - c[lab == 0].mean()
            CP = c[PM]
            nul = CP[:, :n_t].mean(1) - CP[:, n_t:].mean(1)
            ps.append(float((nul >= obs).mean()))
        V["a_shuf_p"] = ps
        V["c_a_shuf"] = V.a_shuf_p < 0.05
    V.to_csv(f"{OUT}/screen_north.csv", index=False)
    print(f"wrote {OUT}/screen_north.csv  {len(V)} templates; "
          f"tremor days {n_t}, quiet days {n_q}; "
          f"hour base any {len(ANY_ALL)/len(EPS):.3f} loud {len(LOUD_ALL)/len(EPS):.3f}")


if __name__ == "__main__":
    main()
