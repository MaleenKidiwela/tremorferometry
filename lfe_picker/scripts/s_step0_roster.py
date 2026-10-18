#!/usr/bin/env python
"""P3D STEP 0 -- roster, site merge, the PREMISE-CORRECTION site table, the Q2(iii) prediction
table, and the hashed freeze.  Nothing downstream may be opened before this has run.

Writes (all under results/p3d/):
  stations_p3d.csv        the P3D roster (frozen S15 in STAS_ORDER + 10 new 3-C + 2 Z-only helpers)
  sites_p3d.csv           one row per SITE after the co-located merge, with the <=60 km geometry rule
  premise_family_sites.csv per-family site counts within 45/55/60 km BEFORE vs AFTER the ten stations
  pred_q2iii.csv          predicted S-SNR + weight for every (family, new station)
  pred_q2iii_recall.csv   predicted per-family recall at matched FA, old network vs dense network
  freeze.json             the hashed pre-registration (hidden families, withheld sites, split,
                          discovery threshold rule, prediction-table hash)

The prediction model (Merlin Q2(iii)):
  * S-SNR vs HYPOCENTRAL distance, fitted per STATION CLASS on e5r_template_qc_2012.csv at MATCHED
    stack size (snr_N30 and snr_N60 columns), with the per-station residual (site term) reported as
    the honest uncertainty -- the class model cannot know a new site's noise.
  * weight w = 10^(S-SNR/20), then the FROZEN 3x-median clip and median normalisation
    (e5r_templates weight rule).
  * recall through the E5 prediction model (`e5_prediction_model.py`, read as a library of formulae:
    bg + e_net(e_s) with e_s = e_class(epi, klass), recall(T, raw_med, sigma_pooled)).  Matched FA is
    read as matched T in family-day MAD units: both network statistics are MAD-normalised, which is
    exactly the assumption the E5 model already makes.  Stated, not hidden.
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                      # noqa: E402

RES = S.IN                                                                # lfe_catalog_repair/results
OUT = S.P3D
MREF = 1.95
Q_PUB = 0.89
RING = ["PGC", "LZB", "VGZ", "SNB", "YOUB", "PFB"]


# ------------------------------------------------------------------ E5 model, as a library -----
def e5_model():
    """The measured inputs of e5_prediction_model.py, re-read (never re-fitted) from its own files."""
    NET = pd.read_csv(f"{RES}/e4_step1_network.csv")
    bg = float(NET[(NET.L == 8.0) & (NET.magbin == "all")].net_rand_bg.iloc[0])
    REG = pd.read_csv(f"{RES}/e4_step1_registration.csv")
    S2 = pd.read_csv(f"{RES}/e4_step2_registration.csv")
    Qq = pd.read_csv(f"{RES}/e4_step2_q.csv")
    G = pd.read_csv(f"{RES}/e4_step2_networkgain.csv")
    QC = pd.read_csv(f"{RES}/e5r_template_qc_2012.csv")
    E = pd.read_parquet(f"{RES}/e4_tt_events_ring_p2.parquet")
    g_eq = float(G[(G.leg == "A_test") & (G.L == 8.0)].ratio_vs_equalweight.iloc[0])
    h = E[(E.L == 8.0) & (E.mag > 1.8) & (E.mag <= 2.1) & (E.net_mad > 0.3)]
    p10, p90 = np.percentile(h.net_mad, [10, 90])
    sig = float((np.log10(p90) - np.log10(p10)) / (2 * norm.ppf(0.9)))     # POOLED, the E5 primary
    r8 = REG[(REG.L == 8.0) & (REG.magbin == "M1.8-2.1")].set_index("sta")
    qratio = Qq[Qq.L == 8.0].set_index("sta")
    q_meas_ratio = float(qratio.ratio_2010_over_seisbench.median())
    c8 = S2[(S2.leg == "C_test") & (S2.L == 8.0) & (S2.magbin == "M1.8-2.1")].set_index("sta")
    dist = QC.groupby("sta").epi_km.median()
    sta_e = {s: float(r8.excess_meas[s]) for s in RING}
    sta_d = {s: float(r8.epi_km[s]) for s in RING}
    for s in ("B011", "B926", "LRIV"):
        sta_e[s] = float(c8.excess3[s]) / q_meas_ratio
        sta_d[s] = float(dist[s])
    kl = {s: "surface_bb" for s in RING + ["LRIV"]}
    kl.update(B011="borehole", B926="borehole")
    P = pd.DataFrame([dict(sta=s, epi_km=sta_d[s], e=sta_e[s], klass=kl[s]) for s in sta_e])
    def med(sel):
        v = P[sel].e
        return float(v.median()) if len(v) else np.nan
    EMP = {("surface_bb", "near"): med((P.klass == "surface_bb") & (P.epi_km < 40)),
           ("surface_bb", "far"): med((P.klass == "surface_bb") & (P.epi_km >= 40)),
           ("borehole", "near"): med(P.klass == "borehole"),
           ("borehole", "far"): med(P.klass == "borehole")}

    def e_class(d, k):
        base = EMP[(k, "near" if d < 40 else "far")]
        return base if d <= 60 else base * (60.0 / d) ** 2

    def e_net(es, mode="weighted", g=0.86):
        es = np.asarray(es, float)
        if es.size == 0:
            return 0.0
        return g * es.sum() / np.sqrt(len(es)) if mode == "equal" else g * np.sqrt((es ** 2).sum())

    def recall(T, raw_med):
        return float(1.0 - norm.cdf((np.log10(T) - np.log10(max(raw_med, 1e-3))) / sig))

    return dict(bg=bg, sig=sig, g_eq=g_eq, EMP=EMP, e_class=e_class, e_net=e_net, recall=recall)


# ------------------------------------------------------------------ S-SNR vs distance -----------
def fit_snr_distance():
    """S-SNR (dB) vs log10 hypocentral distance, per station class, at matched stack size.
    Returns (fit dict, per-station residual site terms, the fitted rows)."""
    QC = pd.read_csv(f"{RES}/e5r_template_qc_2012.csv")
    F = S.families()
    QC["dep"] = QC.bkey.map(F.dep)
    QC["hyp_km"] = np.hypot(QC.epi_km, QC.dep)
    u = QC[QC.usable & QC.template].copy()
    fits, site_terms = {}, []
    for col in ("snr_N30", "snr_N60", "snr_db"):
        for kl, g in u.groupby("klass"):
            h = g.dropna(subset=[col])
            x = np.log10(h.hyp_km.values)
            b, a = np.polyfit(x, h[col].values, 1)
            pred = a + b * x
            res = h[col].values - pred
            fits[(kl, col)] = dict(a=float(a), b=float(b), n=int(len(h)),
                                   r=float(np.corrcoef(x, h[col].values)[0, 1]),
                                   rms_resid_db=float(np.sqrt((res ** 2).mean())))
            for sta, gg in pd.DataFrame(dict(sta=h.sta.values, res=res)).groupby("sta"):
                site_terms.append(dict(col=col, klass=kl, sta=sta, n=len(gg),
                                       site_term_db=float(np.median(gg.res))))
    return fits, pd.DataFrame(site_terms), u


def main():
    os.makedirs(OUT, exist_ok=True)
    ST = S.build_station_table()
    ST.reset_index().to_csv(f"{OUT}/stations_p3d.csv", index=False)
    SI = S.site_table(ST)
    SI.to_csv(f"{OUT}/sites_p3d.csv", index=False)
    F = S.families()
    hid, iso, ho = S.holdouts()

    # ---------------------------------------------------- premise correction, per family
    # BEFORE = every site that already existed (role frozen OR mixed: B926q and B001q are the SAME
    # SITES as before -- CLRS is 25 m from B926 and OSQM 17 m from B001, so they add NO geometry).
    # AFTER  = those plus the genuinely NEW sites (role == "new").  This is Merlin's flaw 3.
    old_sites = SI[SI.role.isin(["frozen", "mixed"])]
    new_sites = SI                                              # after = every site, merged
    genuinely_new = set(SI[SI.role == "new"].site)
    rows = []
    for b in F.index:
        d_old = S.hyp_km(F.lat[b], F.lon[b], F.dep[b], old_sites.lat.values, old_sites.lon.values)
        d_new = S.hyp_km(F.lat[b], F.lon[b], F.dep[b], new_sites.lat.values, new_sites.lon.values)
        gained = [s for s, d in zip(new_sites.site, d_new)
                  if d <= 45 and s in genuinely_new]
        gained55 = [s for s, d in zip(new_sites.site, d_new)
                    if d <= 55 and s in genuinely_new]
        rows.append(dict(bkey=int(b), lat=F.lat[b], lon=F.lon[b], dep=F.dep[b],
                         hidden=int(b) in hid, isolated=int(b) in iso,
                         n45_before=int((d_old <= 45).sum()), n45_after=int((d_new <= 45).sum()),
                         n55_before=int((d_old <= 55).sum()), n55_after=int((d_new <= 55).sum()),
                         n60_before=int((d_old <= 60).sum()), n60_after=int((d_new <= 60).sum()),
                         new_sites_45="+".join(sorted(gained)),
                         new_sites_55="+".join(sorted(gained55)),
                         n_new45=len(gained), n_new55=len(gained55),
                         d_symb=float(S.hyp_km(F.lat[b], F.lon[b], F.dep[b],
                                               ST.lat["SYMB"], ST.lon["SYMB"]))))
    PR = pd.DataFrame(rows)
    PR["d45"] = PR.n_new45                 # GENUINELY new sites within 45 km (not co-located pairs)
    PR["d55"] = PR.n_new55
    PR.to_csv(f"{OUT}/premise_family_sites.csv", index=False)

    # ---------------------------------------------------- Q2(iii) prediction table
    fits, site_terms, u = fit_snr_distance()
    site_terms.to_csv(f"{OUT}/pred_site_terms.csv", index=False)
    M = e5_model()
    W2012 = pd.read_csv(f"{RES}/e5r_weights_snr_2012.csv")
    w_med_ref = float(W2012.w.median())
    targets = [s for s in ST.index if ST.transfer_target[s] and ST["mode"][s] == "3C"]
    prows = []
    for b in F.index:
        for s in targets:
            d = float(S.hyp_km(F.lat[b], F.lon[b], F.dep[b], ST.lat[s], ST.lon[s]))
            kl = ST.klass[s]
            rec = {}
            for col in ("snr_N30", "snr_N60", "snr_db"):
                f = fits[(kl, col)]
                rec[f"pred_{col}"] = f["a"] + f["b"] * np.log10(d)
            w_raw = 10 ** (rec["pred_snr_N60"] / 20.0)
            prows.append(dict(bkey=int(b), sta=s, site=ST.site[s], klass=kl, hyp_km=d,
                              epi_km=float(S.epi_km(F.lat[b], F.lon[b], ST.lat[s], ST.lon[s])),
                              **rec, w_raw=w_raw,
                              pred_e_excess=float(M["e_class"](
                                  float(S.epi_km(F.lat[b], F.lon[b], ST.lat[s], ST.lon[s])), kl))))
    P = pd.DataFrame(prows)
    wmed = float(P.w_raw.median())
    P["w_pred"] = np.minimum(P.w_raw, 3.0 * wmed) / wmed                  # frozen 3x-median clip
    P["pred_rms_resid_db"] = P.klass.map(
        {k: fits[(k, "snr_N60")]["rms_resid_db"] for k in P.klass.unique()})
    P.to_csv(f"{OUT}/pred_q2iii.csv", index=False)

    # ---------------------------------------------------- predicted recall, old vs dense
    # network for the prediction: SITES, one e_s per site (the merge), at matched T in MAD units.
    fz = json.load(open(f"{S.EXT}/results/frozen_operating_point.json"))
    T_core = float(fz["thresholds"]["T_core"]) if "thresholds" in fz else 8.75
    old_net = [s for s in SI.site if s in S.old7_sites()]
    # the dense set is EPOCH-dependent: BUCKS/TURTL/ODUC do not exist before 2023.
    epoch_dense = {}
    for eid, a, b_ in S.NET_EPOCHS:
        alive = {s for s in ST.index if str(ST.start[s]) <= str(b_) and str(ST.end[s]) >= str(a)}
        d = [s for s in S.dense_sites(SI)
             if any(m in alive for m in SI[SI.site == s].members.iloc[0].split("+"))]
        epoch_dense[eid] = sorted(set(old_net) | set(d))

    def es_of(b, sites):
        out = []
        for si in sites:
            g = SI[SI.site == si].iloc[0]
            kl = ST.klass[g.members.split("+")[0]]
            dd = float(S.epi_km(F.lat[b], F.lon[b], g.lat, g.lon))
            out.append(M["e_class"](dd, kl))
        return out

    rrows = []
    PRi = PR.set_index("bkey")
    for b in F.index:
        r_old = M["recall"](T_core, M["bg"] + M["e_net"](es_of(b, old_net)))
        row = dict(bkey=int(b), hidden=int(b) in hid, isolated=int(b) in iso,
                   n_sites_old=len(old_net), recall_old=r_old,
                   d45=int(PRi.d45[int(b)]), d55=int(PRi.d55[int(b)]))
        for eid, _a, _b in S.NET_EPOCHS:
            r_den = M["recall"](T_core, M["bg"] + M["e_net"](es_of(b, epoch_dense[eid])))
            row[f"recall_dense_{eid}"] = r_den
            row[f"gain_frac_{eid}"] = (r_den - r_old) / max(r_old, 1e-6)
            row[f"n_sites_dense_{eid}"] = len(epoch_dense[eid])
        rrows.append(row)
    RC = pd.DataFrame(rrows)
    RC.to_csv(f"{OUT}/pred_q2iii_recall.csv", index=False)
    sel = RC[RC.d45 >= 1]
    gain_by_epoch = {eid: float((sel[f"recall_dense_{eid}"].sum() - sel.recall_old.sum())
                                / max(sel.recall_old.sum(), 1e-9)) for eid, _a, _b in S.NET_EPOCHS}
    gain_pooled = float(np.mean(list(gain_by_epoch.values())))
    dense = epoch_dense["D3"]
    dense_all = dense

    # ---------------------------------------------------- the freeze
    fr = dict(
        created_utc=pd.Timestamp.utcnow().isoformat(),
        spec="results/merlin_review_p3d_2026-09-19.md",
        hidden_families=sorted(hid),
        hidden_isolated=sorted(iso),
        holdouts_sha256=ho["sha256_of_hidden_family_list"],
        withheld_sites={s: SI[SI.site == s].members.iloc[0] for s in S.WITHHELD_SITES},
        withheld_site_rule=("never enters any network sum in any P3D arm; used only as Gate-1 "
                            "confirmers and to certify marginal detections (a3b stack test)"),
        site_merge={k: v for k, v in S.SITE_CLUSTERS.items()},
        stas_order=S.STAS_FROZEN + S.NEW_3C + S.HELPERS_Z,
        episode_split=dict(calibration=S.CAL_EPISODES, test=S.TEST_EPISODES,
                           late_stations=S.LATE_STATIONS, late_build=S.LATE_BUILD_EPISODE,
                           late_test=S.LATE_TEST_EPISODES,
                           rule=("stations alive before 2023 calibrate on S01-S13 and are tested on "
                                 "S14/S18/S21/S22; BUCKS/TURTL/ODUC build on S18 and are tested on "
                                 "S21/S22 only; 2011 and 2013 are not opened")),
        network_epochs=[dict(id=e, start=a, end=b) for e, a, b in S.NET_EPOCHS],
        quiet_months=S.QUIET_MONTHS,
        geometry_rule=dict(kind="median hypocentral distance to the 130-family field",
                           threshold_km=S.GEOM_RULE_KM,
                           dense_sites=S.dense_sites(SI), dense_plus_old=dense_all, old7=old_net,
                           med_hyp_km={r.site: round(r.med_hyp_km, 2) for r in SI.itertuples()}),
        discovery_threshold_rule=dict(
            statistic="q_pairlink2 network site-sum z, index x 1800 s shift null, 27 s overlap excl.",
            rule="null quantile giving <= 30 EXPECTED FALSE EDGES over the pooled candidate set",
            note=("with ~3.4e10 pairs this is the ~1e-9 quantile; the exact quantile is recomputed "
                  "from the realised pair count of the run and reported, the RULE is what is frozen"),
            expected_false_edges=30),
        step6_licence=dict(rule="predicted delta_full >= 1.5 sigma, i.e. delta_dense/delta_old >= 4",
                           delta_old_p3="0.38 sigma at 7 sites (P3 RESULTS 5.5)"),
        prediction_table=dict(
            file="results/p3d/pred_q2iii.csv", n_rows=int(len(P)),
            recall_file="results/p3d/pred_q2iii_recall.csv",
            model=("S-SNR = a + b log10(hypocentral km) per station class, fitted on "
                   "e5r_template_qc_2012.csv at matched stack size; w = 10^(SNR/20) with the frozen "
                   "3x-median clip; recall through the E5 prediction model "
                   "(bg + e_net(e_class(epi, klass)), pooled log-spread), matched FA read as matched "
                   "T in family-day MAD units"),
            fits={f"{k[0]}|{k[1]}": v for k, v in fits.items()},
            T_core=T_core, bg=M["bg"], sigma_log=M["sig"],
            pred_symb_median_snr_N60_db=float(P[P.sta == "SYMB"].pred_snr_N60.median()),
            pred_symb_median_w=float(P[P.sta == "SYMB"].w_pred.median()),
            PREDICTED_RECALL_GAIN_d45ge1_by_epoch=gain_by_epoch,
            PREDICTED_POOLED_RECALL_GAIN_d45ge1=gain_pooled,
            predicted_gain_by_d45={str(k): float(v) for k, v in
                                   RC.groupby("d45").gain_frac_D3.median().items()},
            epoch_dense_sites=epoch_dense),
        premise_correction=dict(
            file="results/p3d/premise_family_sites.csv",
            median_n45_before=float(PR.n45_before.median()),
            median_n45_after=float(PR.n45_after.median()),
            median_n55_before=float(PR.n55_before.median()),
            median_n55_after=float(PR.n55_after.median()),
            n_families_gaining_a_new_site_le45=int((PR.d45 >= 1).sum())),
        tables_sha256=dict(
            stations=S.sha256_of(pd.read_csv(f"{OUT}/stations_p3d.csv").to_dict("records")),
            sites=S.sha256_of(SI.to_dict("records")),
            premise=S.sha256_of(PR.to_dict("records")),
            prediction=S.sha256_of(P.round(6).to_dict("records")),
            prediction_recall=S.sha256_of(RC.round(6).to_dict("records"))),
    )
    fr["freeze_sha256"] = S.sha256_of(fr)
    with open(f"{OUT}/freeze.json", "w") as f:
        json.dump(fr, f, indent=1, default=str)

    print(f"[step0] roster {len(ST)} stations, {len(SI)} sites")
    print(f"[step0] GEOMETRY RULE (median hypocentral <=60 km, withheld out) = {S.dense_sites(SI)}")
    print(f"[step0] epoch dense networks (old7 + rule-passing new): {epoch_dense}")
    print(f"[step0] premise: median sites <=45 km {PR.n45_before.median():.0f} -> "
          f"{PR.n45_after.median():.0f}; <=55 km {PR.n55_before.median():.0f} -> "
          f"{PR.n55_after.median():.0f}; families gaining a NEW site <=45 km: {(PR.d45>=1).sum()}/130")
    print(f"[step0] predicted SYMB S-SNR(N60) median {P[P.sta=='SYMB'].pred_snr_N60.median():.2f} dB, "
          f"w {P[P.sta=='SYMB'].w_pred.median():.2f}")
    print(f"[step0] PREDICTED recall gain (families with >=1 GENUINELY NEW site <=45 km) by epoch: " +
          ", ".join(f"{k} {v*100:+.1f} %" for k, v in gain_by_epoch.items()) +
          f"  | mean {gain_pooled*100:+.1f} %")
    print(f"[step0] freeze sha256 {fr['freeze_sha256']}")


if __name__ == "__main__":
    main()
