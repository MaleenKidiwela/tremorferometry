#!/usr/bin/env python
"""E5 step 2 -- the CORRECTED prediction model (Merlin E4-audit F2, F9, F11; prescription 2).

What was wrong in scripts/e4_prediction_table.py (lines 139-140, 171-173): the EXCESS statistic
(median-at-events minus the random-time max-over-scan median, 2.73) was compared with the detection
threshold, but detection is on the RAW statistic: raw = excess + background.  Every Mc was
log10(T/(T-2.73))/alpha too high (0.48 mag at 8 MAD, 1.31 at 4 MAD).

MODEL (every input measured in E4 step 1/2 on the 2010 data, file named next to it)
  background   bg = median raw network statistic at random times (e4_step1_network.csv net_rand_bg, L = 8)
  registration e_s(class, d) = per-station 3-C excess at M1.8-2.1 for PUBLISHED-quality templates:
               surface broadband < 40 km / 40-60 km from the ring (e4_step1_registration.csv); borehole
               from the held-out transfers B011/B926 (e4_step2_registration.csv, leg C_test) divided by the
               measured template-quality ratio of a ~90-event stack (e4_step2_q.csv), i.e. brought to
               published quality; > 60 km: far value x (60/d)^2 (assumption, stated).
  template     q(N) = 1/sqrt(1 + 1/(r1_eff^2 N)) with r1_eff CALIBRATED on the measured held-out ratio
               ours(N~90)/published (e4_step2_q.csv) -- the metadata r1 understates it ~1.9x (E4 step 1).
  length       e(12 s)/e(8 s) = median per-station ratio from e4_step1_registration.csv.
  magnitude    e_s(M) = e_s(1.95) 10^(alpha (M - 1.95)); alpha = 0.38 measured WITHIN family (E4); the
               incoherent-power slope 0.46-0.58 is a sensitivity (F9: 0.38 is compressed by the
               max-over-scan nonlinearity and survivorship).  Recall is reported as CURVES, not one Mc.
  network      equal-weight  e_net = g_eq sqrt(n)^-1 sum_s e_s     g_eq = measured 0.86 (e4_step2_networkgain.csv)
               SNR-weighted  e_net = g_w  sqrt(sum_s e_s^2)         g_w = 0.86 (same degradation) and 1/1.38 (Merlin, conservative)
  detection    raw_med(M) = bg + e_net(M); recall(T, M) = P(raw >= T) with log10(raw) ~ N(log10 raw_med, sigma),
               sigma = the MEASURED log-spread of the raw statistic at fixed magnitude (within-family, L = 8).
IN-SAMPLE ACCEPTANCE: the model with each 2010 family's ACTUAL ring station set must give the ring raw
median 5.9 +- 0.5 and frac >= 8 of 0.21 +- 0.05 at M1.8-2.1 (measured: e4_tt_events_ring_p2.parquet).
Then: Bostock's 21-station array Mc(8) (target 1.6-1.9), the 7-station 2010-13 array vs his modal 1.75,
the 15-station 2012 set (user decision) with and without LZB/YOUB/PFB, and the region-year table at 8 MAD.
Outputs results/e5_prediction_{spread,insample,arrays,2012_families,2012_recall,regionyear}.csv + e5_prediction_model.md
"""
from __future__ import annotations
import argparse, json
import numpy as np, pandas as pd
from scipy.stats import norm

RES = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
SD = f"{RES}/external/other/seisbench_lfestacks_cascadia"
R_KM_DEG = 111.195
MREF = 1.95
Q_PUB = 0.89
RING = ["PGC", "LZB", "VGZ", "SNB", "YOUB", "PFB"]
SET15 = ["PGC", "VGZ", "NLLB", "SNB", "LRIV", "B001", "B003", "B004", "B005", "B006", "B007", "B009", "B010",
         "B011", "B926"]
EXCL3 = ["LZB", "YOUB", "PFB"]
MB12 = [(-9, 1.6, "M<1.6"), (1.6, 1.8, "M1.6-1.8"), (1.8, 2.0, "M1.8-2.0"), (2.0, 9, "M>=2.0")]


def epi_km(lat1, lon1, lat2, lon2):
    cs = np.cos(np.radians(48.6))
    return np.hypot((np.asarray(lon1) - lon2) * R_KM_DEG * cs, (np.asarray(lat1) - lat2) * R_KM_DEG)


def recall(T, raw_med, sig):
    raw_med = np.maximum(np.asarray(raw_med, float), 1e-3)
    return 1.0 - norm.cdf((np.log10(T) - np.log10(raw_med)) / sig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--thr", default="8,9,10", help="detection thresholds (MAD); 8 + the step-1 operating point")
    ap.add_argument("--alphas", default="0.38,0.46,0.58,0.0,0.2", help="0.0 = magnitude-FLAT bracket (the pooled 2010 recall is flat in M)")
    ap.add_argument("--min-build-n", type=int, default=20, help="min 2012 build-half events for a template")
    a = ap.parse_args()
    THR = [float(x) for x in a.thr.split(",")]
    ALPHAS = [float(x) for x in a.alphas.split(",")]
    md = []

    # ------------------------------------------------------------ measured inputs
    E = pd.read_parquet(f"{RES}/e4_tt_events_ring_p2.parquet")
    NET = pd.read_csv(f"{RES}/e4_step1_network.csv")
    bg = float(NET[(NET.L == 8.0) & (NET.magbin == "all")].net_rand_bg.iloc[0])
    REG = pd.read_csv(f"{RES}/e4_step1_registration.csv")
    S2 = pd.read_csv(f"{RES}/e4_step2_registration.csv")
    Q = pd.read_csv(f"{RES}/e4_step2_q.csv")
    G = pd.read_csv(f"{RES}/e4_step2_networkgain.csv")
    QC = pd.read_csv(f"{RES}/e4_transfer_stackqc_t2.csv")
    g_eq = float(G[(G.leg == "A_test") & (G.L == 8.0)].ratio_vs_equalweight.iloc[0])

    # spread of the RAW statistic at fixed magnitude (L = 8 and 12, M1.8-2.1)
    sp_rows = []
    for L in (8.0, 12.0):
        h = E[(E.L == L) & (E.mag > 1.8) & (E.mag <= 2.1) & (E.net_mad > 0.3)].copy()
        lg = np.log10(h.net_mad)
        p10, p50, p90 = np.percentile(h.net_mad, [10, 50, 90])
        sig_pooled = float((np.log10(p90) - np.log10(p10)) / (2 * norm.ppf(0.9)))
        fm = lg.groupby(h.fam).transform("median")
        res = (lg - fm)
        sig_within = float(1.4826 * np.median(np.abs(res - res.median())))
        sig_within_sd = float(res.std())
        sp_rows.append(dict(L=L, n=len(h), p10=p10, p50=p50, p90=p90, sigma_log_pooled=sig_pooled,
                            sigma_log_within_family_robust=sig_within, sigma_log_within_family_sd=sig_within_sd,
                            frac_ge8=float((h.net_mad >= 8).mean()), frac_ge10=float((h.net_mad >= 10).mean())))
    SP = pd.DataFrame(sp_rows); SP.to_csv(f"{RES}/e5_prediction_spread.csv", index=False)
    SIGW = float(SP[SP.L == 8.0].sigma_log_within_family_robust.iloc[0])
    SIGP = float(SP[SP.L == 8.0].sigma_log_pooled.iloc[0])
    # PRIMARY spread = the POOLED log-spread at fixed magnitude (Merlin's p10/50/90 prescription): it carries
    # the between-family registration variance that a station-class model cannot (in-sample the within-family
    # spread under-predicts frac>=8 by 1.8x, the pooled one reproduces it within the bar); within = sensitivity.
    SIG = SIGP
    md.append(f"bg = {bg:.3f}; raw p10/50/90 at M1.8-2.1, L=8: {sp_rows[0]['p10']:.2f}/{sp_rows[0]['p50']:.2f}/"
              f"{sp_rows[0]['p90']:.2f} -> sigma_log pooled {SIGP:.3f} (PRIMARY), within-family robust {SIGW:.3f} "
              f"(sd {sp_rows[0]['sigma_log_within_family_sd']:.3f}); g_eq = {g_eq:.3f}")

    # per-station registration (published quality, L = 8, M1.8-2.1)
    r8 = REG[(REG.L == 8.0) & (REG.magbin == "M1.8-2.1")].set_index("sta")
    r12 = REG[(REG.L == 12.0) & (REG.magbin == "M1.8-2.1")].set_index("sta")
    L12 = float(np.median([r12.excess_meas[s] / r8.excess_meas[s] for s in RING]))
    qratio = Q[Q.L == 8.0].set_index("sta")
    q_meas_ratio = float(qratio.ratio_2010_over_seisbench.median())        # ours(N~90) / published
    n_q = float(qratio.n_events_stacked.median())
    q_meas = q_meas_ratio * Q_PUB
    r1_eff = float(np.sqrt(q_meas ** 2 / ((1 - q_meas ** 2) * n_q)))
    def q_of(N):
        N = np.maximum(np.asarray(N, float), 1.0)
        return 1.0 / np.sqrt(1.0 + 1.0 / (r1_eff ** 2 * N))
    c8 = S2[(S2.leg == "C_test") & (S2.L == 8.0) & (S2.magbin == "M1.8-2.1")].set_index("sta")
    dist = QC.groupby("sta").epi_km.median()
    sta_e = {s: float(r8.excess_meas[s]) for s in RING}
    for s in ("B011", "B926", "LRIV"):                          # transfers -> published-quality equivalent
        sta_e[s] = float(c8.excess3[s]) / q_meas_ratio
    sta_d = {s: float(r8.epi_km[s]) for s in RING}
    for s in ("B011", "B926", "LRIV"):
        sta_d[s] = float(dist[s])
    klass = {s: "surface_bb" for s in RING + ["LRIV"]}; klass.update(B011="borehole", B926="borehole")
    P = pd.DataFrame([dict(sta=s, epi_km=sta_d[s], excess_pubq=sta_e[s], klass=klass[s]) for s in sta_e])
    def med(sel):
        v = P[sel].excess_pubq; return float(v.median()) if len(v) else np.nan
    EMP = {("surface_bb", "near"): med((P.klass == "surface_bb") & (P.epi_km < 40)),
           ("surface_bb", "far"): med((P.klass == "surface_bb") & (P.epi_km >= 40)),
           ("borehole", "near"): med(P.klass == "borehole"), ("borehole", "far"): med(P.klass == "borehole")}
    def e_class(d, kl):
        base = EMP[(kl, "near" if d < 40 else "far")]
        return base if d <= 60 else base * (60.0 / d) ** 2
    md.append(f"per-station published-quality excess (L=8, M1.8-2.1): " +
              ", ".join(f"{r.sta} {r.epi_km:.0f} km {r.excess_pubq:.2f}" for r in P.sort_values('epi_km').itertuples()))
    md.append(f"class model: surface <40 km {EMP[('surface_bb','near')]:.2f}, 40-60 km {EMP[('surface_bb','far')]:.2f}; "
              f"borehole {EMP[('borehole','near')]:.2f} (B011/B926 at 41-47 km, both bands); >60 km: far x (60/d)^2. "
              f"e(12 s)/e(8 s) = {L12:.3f}. q-ratio ours(N={n_q:.0f})/published = {q_meas_ratio:.2f} -> r1_eff = {r1_eff:.3f}; "
              f"q(25/50/100/200/500) = " + "/".join(f"{q_of(n):.2f}" for n in (25, 50, 100, 200, 500)))

    def e_net(es, mode, g_w):
        es = np.asarray(es, float)
        if len(es) == 0:
            return 0.0
        return g_eq * es.sum() / np.sqrt(len(es)) if mode == "equal" else g_w * np.sqrt((es ** 2).sum())

    # ------------------------------------------------------------ IN-SAMPLE check (2010 ring, actual station sets)
    ins = []
    E8 = E[(E.L == 8.0) & (E.mag > 1.8) & (E.mag <= 2.1)]
    for src in ("station_measured", "class_model"):
        for sig_lab, sig in (("within", SIGW), ("pooled", SIGP)):
            pred_med, pred_rec = [], []
            for fam, h in E8.groupby("fam"):
                ss = [s for s in RING if h[f"sta_{s}"].notna().mean() > 0.5]
                if not ss:
                    continue
                es = [sta_e[s] if src == "station_measured" else e_class(sta_d[s], "surface_bb") for s in ss]
                rm = bg + e_net(es, "equal", g_eq)
                pred_med += [rm] * len(h); pred_rec += list(recall(8.0, np.full(len(h), rm), sig))
            ins.append(dict(source=src, sigma=sig_lab, n_ev=len(pred_med), pred_raw_median=float(np.median(pred_med)),
                            meas_raw_median=float(E8.net_mad.median()), pred_frac_ge8=float(np.mean(pred_rec)),
                            meas_frac_ge8=float((E8.net_mad >= 8).mean())))
    INS = pd.DataFrame(ins); INS.to_csv(f"{RES}/e5_prediction_insample.csv", index=False)
    ok = INS[(INS.source == "station_measured") & (INS.sigma == "pooled")].iloc[0]
    md.append(f"IN-SAMPLE (actual station sets, station-measured e_s, equal weight): predicted raw median "
              f"{ok.pred_raw_median:.2f} vs measured {ok.meas_raw_median:.2f} (bar 5.9+-0.5); frac>=8 {ok.pred_frac_ge8:.3f} vs "
              f"{ok.meas_frac_ge8:.3f} (bar 0.21+-0.05) -> {'PASS' if abs(ok.pred_raw_median-5.9)<=0.5 and abs(ok.pred_frac_ge8-0.21)<=0.05 else 'FAIL'}")
    ok2 = INS[(INS.source == "class_model") & (INS.sigma == "pooled")].iloc[0]
    ok3 = INS[(INS.source == "station_measured") & (INS.sigma == "within")].iloc[0]
    CAL = (ok.meas_raw_median - bg) / (ok.pred_raw_median - bg)
    md.append(f"IN-SAMPLE with the CLASS model instead of station values: raw {ok2.pred_raw_median:.2f}, frac>=8 {ok2.pred_frac_ge8:.3f}; "
              f"with the WITHIN-family spread ({SIGW:.3f}) frac>=8 = {ok3.pred_frac_ge8:.3f} (under-predicts the tail: between-family variance). "
              f"Excess calibration factor measured/predicted = {CAL:.2f} (reported as a labelled bracket, NOT applied to the primary numbers)")

    # measured recall-by-magnitude curve (for the figure) at L = 8 / 12
    cur = []
    for L in (8.0, 12.0):
        h = E[E.L == L]
        for lo in np.arange(1.2, 2.5, 0.1):
            k = (h.mag > lo) & (h.mag <= lo + 0.1)
            if k.sum() >= 30:
                cur.append(dict(L=L, m_lo=lo, m_mid=lo + 0.05, n=int(k.sum()),
                                **{f"rec{int(T)}": float((h.net_mad[k] >= T).mean()) for T in (6, 7, 8, 9, 10)}))
    pd.DataFrame(cur).to_csv(f"{RES}/e5_prediction_measured_curve.csv", index=False)
    # WITHIN-family measured curve: magnitude centred on each family's median (removes the distance confound)
    cw = []
    for L in (8.0, 12.0):
        h = E[E.L == L].copy(); h["dM"] = h.mag - h.groupby("fam").mag.transform("median")
        for lo in np.arange(-0.6, 0.6, 0.1):
            k = (h.dM > lo) & (h.dM <= lo + 0.1)
            if k.sum() >= 100:
                cw.append(dict(L=L, dM_lo=lo, dM_mid=lo + 0.05, n=int(k.sum()),
                               **{f"rec{int(T)}": float((h.net_mad[k] >= T).mean()) for T in (6, 8, 10)}))
    CW = pd.DataFrame(cw); CW.to_csv(f"{RES}/e5_prediction_measured_curve_within.csv", index=False)
    g8 = CW[CW.L == 8.0]
    if len(g8) > 3:
        sl = np.polyfit(g8.dM_mid, np.log10(np.maximum(g8.rec8, 1e-3)), 1)[0]
        md.append(f"MEASURED recall vs magnitude, 2010 ring L=8, >=8 MAD: POOLED is flat/non-monotonic "
                  f"({', '.join(f'M{r.m_mid:.2f}:{r.rec8:.2f}' for r in pd.DataFrame(cur).query('L==8.0').itertuples() if 1.3 <= r.m_mid <= 2.3)}); "
                  f"WITHIN-family (centred) rises {g8.rec8.iloc[0]:.3f} at dM {g8.dM_mid.iloc[0]:+.2f} -> {g8.rec8.iloc[-1]:.3f} at dM {g8.dM_mid.iloc[-1]:+.2f} "
                  f"(d log10 recall / dM = {sl:.2f}). The per-bin 2012 numbers inherit the model slope; alpha = 0 is the flat bracket.")

    # ------------------------------------------------------------ Bostock's arrays
    meta = pd.read_csv(f"{SD}/metadata.csv")
    meta = meta[meta.trace_S_arrival_sample.notna()]
    cat = pd.read_csv(f"{RES}/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv")
    cat["time"] = pd.to_datetime(dict(year=cat.year, month=cat.month, day=cat.day, hour=cat.hour, minute=cat["min"],
                                      second=cat.sec.astype(int))).astype("datetime64[ns]")
    w12 = cat[(cat.time >= "2012-08-15") & (cat.time < "2012-10-16")].copy()
    w10 = cat[(cat.time >= "2010-08-10") & (cat.time < "2010-09-15")].copy()
    fam12 = w12.groupby(["lat", "lon", "dep"]).agg(n=("mag", "size"), mag_med=("mag", "median")).reset_index()
    fam12["fid"] = np.arange(len(fam12))
    w12 = w12.merge(fam12[["lat", "lon", "dep", "fid"]], on=["lat", "lon", "dep"])
    def mc_of(enet, T, alpha):
        return MREF + np.log10(max(T - bg, 1e-6) / enet) / alpha if enet > 0 else np.nan
    arr = []
    sta_xy = meta.groupby("station_code")[["station_latitude_deg", "station_longitude_deg"]].first()
    fam_xy = meta.groupby("source_family_id")[["source_latitude_deg", "source_longitude_deg"]].first()
    R7 = ["PGC", "LZB", "VGZ", "SNB", "YOUB", "PFB", "SOKB"]
    for name, stas, allwithin in (("Bostock_21_templated", None, False), ("Bostock_21_all_within60", None, True),
                                  ("ring7_2010-13_templated", R7, False), ("ring7_2010-13_all_within60", R7, True)):
        for fid, g in meta.groupby("source_family_id"):
            if allwithin:      # every station of the array within 60 km of the family, template assumed at each
                cand = sta_xy if stas is None else sta_xy.loc[[s_ for s_ in stas if s_ in sta_xy.index]]
                dd = epi_km(fam_xy.loc[fid, "source_latitude_deg"], fam_xy.loc[fid, "source_longitude_deg"],
                            cand.station_latitude_deg.to_numpy(), cand.station_longitude_deg.to_numpy())
                dd = dd[dd <= 60]
            else:
                if stas is not None:
                    g = g[g.station_code.isin(stas)]
                dd = g.path_ep_distance_km.to_numpy()
            if not len(dd):
                continue
            es = [e_class(d, "surface_bb") for d in dd]   # his templates: q ~ 1 (1,400-2,000 events)
            for mode, gw in (("equal", g_eq), ("snrw", g_eq)):
                en = e_net(es, mode, gw)
                arr.append(dict(array=name, fam=fid, n_sta=len(dd), mode=mode, e_net=en,
                                **{f"Mc{int(T)}": mc_of(en, T, 0.38) for T in THR}))
    ARR = pd.DataFrame(arr)
    summ = ARR.groupby(["array", "mode"]).agg(n_fam=("fam", "size"), n_sta_med=("n_sta", "median"), e_net_med=("e_net", "median"),
                                              **{f"Mc{int(T)}_med": (f"Mc{int(T)}", "median") for T in THR}).reset_index()
    summ.to_csv(f"{RES}/e5_prediction_arrays.csv", index=False)
    for nm in ("Bostock_21_templated", "Bostock_21_all_within60"):
        b21 = summ[(summ["array"] == nm) & (summ["mode"] == "equal")].iloc[0]
        md.append(f"{nm}: median {b21.n_sta_med:.0f} stations/family, e_net {b21.e_net_med:.2f} -> Mc(8) = {b21.Mc8_med:.2f} "
                  f"(target 1.6-1.9; his 2003-05 modal M 1.6-1.7)")
    for nm in ("ring7_2010-13_templated", "ring7_2010-13_all_within60"):
        r7 = summ[(summ["array"] == nm) & (summ["mode"] == "equal")].iloc[0]
        md.append(f"{nm} (ring + SOKB): {r7.n_sta_med:.0f} stations/family, e_net {r7.e_net_med:.2f} -> Mc(8) = {r7.Mc8_med:.2f} "
                  f"vs his 2010-13 modal 1.75: implementation-vs-Bostock deficit {r7.Mc8_med - 1.75:+.2f} mag = "
                  f"{10 ** (0.38 * (r7.Mc8_med - 1.75)):.2f}x in registration (alpha 0.38); i.e. his 7-station matcher needs "
                  f"{(8 - bg) / r7.e_net_med:.2f}x our predicted network excess to put M1.75 at 8 MAD")

    # ------------------------------------------------------------ the 15-station 2012 set
    st = pd.read_csv(f"{RES}/f1_stations.csv").drop_duplicates("sta").set_index("sta")
    def sta_class(s):
        return "borehole" if bool(st.loc[s, "is_borehole_any"]) else "surface_bb"
    famrows, recrows = [], []
    for setname, stas in (("set15", SET15), ("set18_with_LZB_YOUB_PFB", SET15 + EXCL3)):
        for L, Lf in ((8.0, 1.0), (12.0, L12)):
            for alpha in ALPHAS:
                per_fam = {}
                for r in fam12.itertuples():
                    nb = r.n / 2.0
                    if nb < a.min_build_n:
                        per_fam[r.fid] = dict(es=[], nsta=0, q=0.0); continue
                    es, names = [], []
                    for s in stas:
                        d = float(epi_km(r.lat, r.lon, st.loc[s, "lat"], st.loc[s, "lon"]))
                        if d <= 60:
                            es.append(e_class(d, sta_class(s)) * (q_of(nb) / Q_PUB) * Lf); names.append(f"{s}:{d:.0f}")
                    per_fam[r.fid] = dict(es=es, nsta=len(es), q=float(q_of(nb)), names=";".join(names))
                    if L == 8.0 and alpha == ALPHAS[0]:
                        famrows.append(dict(set=setname, fid=r.fid, lat=r.lat, lon=r.lon, dep=r.dep, n_2012=r.n, n_build=nb,
                                            q=float(q_of(nb)), n_sta60=len(es), stations=";".join(names),
                                            e_net_equal=e_net(es, "equal", g_eq), e_net_snrw=e_net(es, "snrw", g_eq),
                                            Mc8_equal=mc_of(e_net(es, "equal", g_eq), 8.0, alpha),
                                            Mc8_snrw=mc_of(e_net(es, "snrw", g_eq), 8.0, alpha)))
                for mode, gw, glab in (("equal", g_eq, "g0.86"), ("snrw", g_eq, "g0.86"), ("snrw", 1 / 1.38, "g1/1.38"),
                                       ("equal", g_eq * CAL, f"g0.86xCAL{CAL:.2f}"), ("snrw", g_eq * CAL, f"g0.86xCAL{CAL:.2f}")):
                    en = w12.fid.map({k: e_net(v["es"], mode, gw) for k, v in per_fam.items()}).to_numpy()
                    rawm = bg + en * 10 ** (alpha * (w12.mag.to_numpy() - MREF))
                    for T in THR:
                        rc = recall(T, rawm, SIG); rc[en <= 0] = 0.0
                        row = dict(set=setname, L=L, alpha=alpha, mode=mode, g=glab, T=T, n_events=len(w12),
                                   recall_all=float(rc.mean()))
                        for lo, hi, lab in MB12:
                            k = (w12.mag >= lo) & (w12.mag < hi)
                            row[f"recall_{lab}"] = float(rc[k].mean()); row[f"n_{lab}"] = int(k.sum())
                        row["fam_with_template"] = int(sum(v["nsta"] > 0 for v in per_fam.values()))
                        row["fam_ge50pct_recall"] = int(sum(1 for k, v in per_fam.items() if v["nsta"] and
                                                            recall(T, bg + e_net(v["es"], mode, gw), SIG) >= 0.5))
                        recrows.append(row)
    FAM = pd.DataFrame(famrows); FAM.to_csv(f"{RES}/e5_prediction_2012_families.csv", index=False)
    REC = pd.DataFrame(recrows); REC.to_csv(f"{RES}/e5_prediction_2012_recall.csv", index=False)
    f15 = FAM[FAM["set"] == "set15"]
    md.append(f"2012 set15: {len(fam12)} families active 2012-08-15..10-15 (per-family n p10/50/90 {np.percentile(fam12.n,[10,50,90])}); "
              f"{int((f15.n_sta60>0).sum())} get >= 1 template (build-half n >= {a.min_build_n}); stations within 60 km per family "
              f"median {f15.n_sta60.median():.0f} (p10 {f15.n_sta60.quantile(.1):.0f}, p90 {f15.n_sta60.quantile(.9):.0f}); "
              f"q(build) median {f15.q.median():.2f}; e_net equal median {f15.e_net_equal.median():.2f}, snrw {f15.e_net_snrw.median():.2f}; "
              f"Mc(8) equal median {f15.Mc8_equal.median():.2f}, snrw {f15.Mc8_snrw.median():.2f}")
    key = REC[(REC.alpha == ALPHAS[0]) & (REC["g"] == "g0.86") & (REC["T"] == 8.0)]
    for rd in REC[(REC["set"] == "set15") & (REC.L == 8.0) & (REC["mode"] == "snrw") & (REC["T"] == 8.0)].to_dict("records"):
        md.append(f"  sensitivity set15 L=8 snrw T=8 alpha={rd['alpha']} g={rd['g']}: recall all {rd['recall_all']:.3f} | "
                  + " ".join(f"{lab} {rd['recall_'+lab]:.3f}" for _, _, lab in MB12))
    for rd in key.to_dict("records"):
        md.append(f"  {rd['set']:24s} L={rd['L']:>4.0f}s {rd['mode']:5s} T={rd['T']:>4.0f}: recall all {rd['recall_all']:.3f} | "
                  + " ".join(f"{lab} {rd['recall_'+lab]:.3f}" for _, _, lab in MB12))

    # ------------------------------------------------------------ region-year table at 8 MAD (+ operating T)
    cd = pd.read_csv(f"{RES}/f1_cell_density.csv")
    mags12 = w12.mag.to_numpy()
    rows = []
    for (region, year, roster), g in cd.groupby(["region", "year", "roster"]):
        if year not in (2010, 2012, 2019, 2021, 2023, 2025):
            continue
        per = {T: [] for T in THR}; recs = {T: [] for T in THR}; nst = []
        for rr in g.itertuples():
            n60 = int(rr.n60)
            if n60 < 2 or not np.isfinite(rr.nearest_km) or rr.nearest_km <= 0:
                nst.append(n60); continue
            nbb = int(getattr(rr, "n60_bb", n60))
            es = [e_class(min(rr.nearest_km * np.sqrt(k), 60.0), "surface_bb" if k <= nbb else "borehole") * (q_of(96) / Q_PUB)
                  for k in range(1, n60 + 1)]
            en_eq, en_w = e_net(es, "equal", g_eq), e_net(es, "snrw", g_eq)
            for T in THR:
                per[T].append((mc_of(en_eq, T, 0.38), mc_of(en_w, T, 0.38)))
                recs[T].append((float(recall(T, bg + en_eq * 10 ** (0.38 * (mags12 - MREF)), SIG).mean()),
                                float(recall(T, bg + en_w * 10 ** (0.38 * (mags12 - MREF)), SIG).mean())))
            nst.append(n60)
        row = dict(region=region, year=year, roster=roster, n_cells=len(g), med_n60=float(np.median(nst)),
                   med_nearest_km=float(g.nearest_km.median()))
        for T in THR:
            if per[T]:
                arrp = np.array(per[T]); arrr = np.array(recs[T])
                row[f"Mc{int(T)}_equal_med"] = float(np.nanmedian(arrp[:, 0])); row[f"Mc{int(T)}_snrw_med"] = float(np.nanmedian(arrp[:, 1]))
                row[f"recall2012dist_T{int(T)}_equal"] = float(np.median(arrr[:, 0])); row[f"recall2012dist_T{int(T)}_snrw"] = float(np.median(arrr[:, 1]))
        rows.append(row)
    RY = pd.DataFrame(rows); RY.to_csv(f"{RES}/e5_prediction_regionyear.csv", index=False)
    md.append("region-year (q at N=96 chain-start, alpha 0.38, sigma pooled): see e5_prediction_regionyear.csv")

    # ------------------------------------------------------------ download volume, 15 stations, 2012-08-15..10-05 at 40 Hz
    days = (pd.Timestamp("2012-10-05") - pd.Timestamp("2012-08-15")).days
    # bytes/sample MEASURED in F1 (results/f1_fleet_inventory.md: STEIM2, 600 day files): 1.416 B/sample at
    # 40 Hz, 0.932 B/sample at 100 Hz; E4 step 6's 0.42 MB per channel-hour at 100 Hz = 1.17 B/sample.
    B40, B100 = 1.416, 0.932
    gb40 = 15 * 3 * 40 * 86400 * days * B40 / 1e9
    sy = pd.read_csv(f"{RES}/f1_station_years.csv"); sy = sy[(sy.year == 2012) & sy.sta.isin(SET15)]
    gbnat = float(sum(3 * fs * 86400 * days * (B40 if fs == 40 else B100) for fs in sy.fs) / 1e9)
    md.append(f"DOWNLOAD (no fetch): 15 stations x 3 ch x {days} d (2012-08-15..10-05) at 40 Hz = {gb40:.1f} GB "
              f"(F1-measured {B40} B/sample at 40 Hz; E4's 1.17 B/sample would give {gb40*1.17/B40:.1f} GB); native rates "
              f"({int((sy.fs==100).sum())} stations at 100 Hz, {int((sy.fs==40).sum())} at 40 Hz) = {gbnat:.1f} GB")
    open(f"{RES}/e5_prediction_model.md", "w").write("# E5 step 2 -- corrected prediction model\n\n" + "\n".join("* " + x for x in md) + "\n")
    print("\n".join(md))
    pd.set_option("display.width", 250)
    print(INS.round(3).to_string(index=False)); print(summ.round(3).to_string(index=False))
    print(REC[(REC.alpha == 0.38) & (REC["g"] == "g0.86")].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
