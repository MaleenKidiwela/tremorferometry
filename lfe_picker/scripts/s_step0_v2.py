#!/usr/bin/env python
"""P3D STEP 0 AMENDMENT (coordinator, 2026-09-19) -- two pre-registered NETWORK VARIANTS.

V1 = the frozen P3D roster (already hashed in freeze.json).
V2 = V1 + the four PQ-network broadbands on southern Vancouver Island found by the all-network
     inventory (PQ.LHLYB, PQ.JRDRB, PQ.ALBHB, PQ.DAOB) + UW.OHOH / UW.ODEER IF they deliver.

This script rebuilds, for BOTH variants:
  * the per-family site table that corrects the premise (sites <=45/55/60 km, before vs after,
    with and without PQ), and
  * the Q2(iii) S-SNR / weight / recall prediction table,
and appends the result to results/p3d/freeze.json as `amendment_v2`, hashed.  It is a PRE-ANALYSIS
addition: it is written BEFORE any PQ transfer stack is opened, and it records the delivery state
of the PQ stations at the time of writing so that the claim is checkable.

Geometry comes from the inventory's own station table (results/p3d_puget_geometry/
station_ranking.csv, read-only); nothing in that workspace is written.
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                       # noqa: E402
import s_step0_roster as R0                                                # noqa: E402

INV = f"{S.RES}/p3d_puget_geometry/station_ranking.csv"
OUT = S.P3D


def v2_station_rows():
    """The V2-only stations from the downloader's CORRECTED station_epochs.csv (StationXML
    coordinates, elevations, component dips and the MEASURED data extents), falling back to the
    geometry inventory only for the per-family counts."""
    E = pd.read_csv(f"{S.P3D}/station_epochs.csv").set_index("sta")
    I = pd.read_csv(INV).set_index("site")
    rows = []
    for s in S.V2_3C:
        if s not in E.index:
            print(f"  [v2] {s} not in station_epochs.csv -- skipped")
            continue
        r = E.loc[s]
        rows.append(dict(sta=s, net=str(r.net), netsta=f"{r.net}.{s}",
                         band=str(r.band), fs_native=int(r.fs), lat=float(r.lat), lon=float(r.lon),
                         elev_m=float(r.elev), klass="surface_bb", mode="3C", role="new_v2",
                         site=s, start=str(r.start), end=str(r.end),
                         az1=S.AZ12[s][0], az2=S.AZ12[s][1], zsign=S.ZSIGN.get(s, 1.0),
                         first_data_month=str(getattr(r, "first_data_month", "")),
                         months_with_data=getattr(r, "months_with_data", np.nan),
                         months_probed=getattr(r, "months_probed", np.nan),
                         fam45_inventory=int(I.fam45[s]) if s in I.index else -1,
                         fam55_inventory=int(I.fam55[s]) if s in I.index else -1))
    return pd.DataFrame(rows)


def delivery_state():
    """What the downloader has actually delivered for the V2 stations, right now."""
    p = f"{OUT}/download_manifest.csv"
    out = {}
    if os.path.exists(p):
        M = pd.read_csv(p, low_memory=False)
        for s in S.V2_3C:
            g = M[M.sta == s]
            out[s] = dict(manifest_rows=int(len(g)),
                          data=int((g.status == "data").sum()) if len(g) else 0)
    idx = S.cache_index(years=list(range(2019, 2027)), stas=S.V2_3C)
    for s in S.V2_3C:
        out.setdefault(s, {})["cached_days"] = sum(1 for (ss, _y, _d) in idx if ss == s)
    return out


def main():
    ap = argparse.ArgumentParser()
    a = ap.parse_args()
    ST1 = S.station_table()
    add = v2_station_rows()
    ST2 = ST1.copy()
    if len(add):
        add = add.set_index("sta")
        add["withheld"] = False
        add["summable"] = True
        add["transfer_target"] = True
        add["band_p3d"] = "HH"
        add["order"] = np.arange(len(add)) + int(ST1.order.max()) + 1
        ST2 = pd.concat([ST1, add[[c for c in ST1.columns if c in add.columns]
                                  + [c for c in add.columns if c not in ST1.columns]]])
        for c in ST1.columns:
            if c not in ST2.columns:
                ST2[c] = np.nan
    ST2.reset_index().rename(columns={"index": "sta"}).to_csv(f"{OUT}/stations_p3d_v2.csv",
                                                              index=False)
    SI1 = S.site_table(ST1)
    SI2 = S.site_table(ST2)
    SI2.to_csv(f"{OUT}/sites_p3d_v2.csv", index=False)
    F = S.families()
    hid, iso, _ = S.holdouts()

    # ---------------- premise table, with and without PQ
    old_sites = SI1[SI1.role.isin(["frozen", "mixed"])]
    new_v1 = set(SI1[SI1.role == "new"].site)
    new_v2 = set(SI2[SI2.role == "new"].site) - new_v1
    assert new_v2 or not len(add), "V2 stations did not reach the site table"
    rows = []
    for b in F.index:
        d_old = S.hyp_km(F.lat[b], F.lon[b], F.dep[b], old_sites.lat.values, old_sites.lon.values)
        d_all = S.hyp_km(F.lat[b], F.lon[b], F.dep[b], SI2.lat.values, SI2.lon.values)
        g1 = [s for s, d in zip(SI2.site, d_all) if d <= 45 and s in new_v1]
        g2 = [s for s, d in zip(SI2.site, d_all) if d <= 45 and s in new_v2]
        g1_55 = [s for s, d in zip(SI2.site, d_all) if d <= 55 and s in new_v1]
        g2_55 = [s for s, d in zip(SI2.site, d_all) if d <= 55 and s in new_v2]
        n45_v1 = int((d_old <= 45).sum()) + len(g1)
        rows.append(dict(bkey=int(b), hidden=int(b) in hid, isolated=int(b) in iso,
                         n45_before=int((d_old <= 45).sum()),
                         n45_after_V1=n45_v1, n45_after_V2=n45_v1 + len(g2),
                         n55_before=int((d_old <= 55).sum()),
                         n55_after_V1=int((d_old <= 55).sum()) + len(g1_55),
                         n55_after_V2=int((d_old <= 55).sum()) + len(g1_55) + len(g2_55),
                         d45_V1=len(g1), d45_V2=len(g1) + len(g2),
                         d55_V1=len(g1_55), d55_V2=len(g1_55) + len(g2_55),
                         new_sites_45_V1="+".join(sorted(g1)),
                         new_sites_45_PQ="+".join(sorted(g2))))
    PR = pd.DataFrame(rows)
    PR.to_csv(f"{OUT}/premise_family_sites_v2.csv", index=False)

    # ---------------- Q2(iii) prediction for the V2-only stations + the recall prediction
    fits, site_terms, _u = R0.fit_snr_distance()
    M = R0.e5_model()
    prows = []
    for b in F.index:
        for s in (list(new_v2) if len(add) else []):
            g = SI2[SI2.site == s].iloc[0]
            d = float(S.hyp_km(F.lat[b], F.lon[b], F.dep[b], g.lat, g.lon))
            rec = {f"pred_{c}": fits[("surface_bb", c)]["a"] + fits[("surface_bb", c)]["b"] * np.log10(d)
                   for c in ("snr_N30", "snr_N60", "snr_db")}
            prows.append(dict(bkey=int(b), sta=s, site=s, klass="surface_bb", hyp_km=d,
                              epi_km=float(S.epi_km(F.lat[b], F.lon[b], g.lat, g.lon)),
                              **rec, w_raw=10 ** (rec["pred_snr_N60"] / 20.0),
                              variant="V2_only"))
    P2 = pd.DataFrame(prows)
    if len(P2):
        wm = float(P2.w_raw.median())
        P2["w_pred"] = np.minimum(P2.w_raw, 3.0 * wm) / wm
    P1 = pd.read_csv(f"{OUT}/pred_q2iii.csv")
    P1["variant"] = "V1"
    pd.concat([P1, P2], ignore_index=True).to_csv(f"{OUT}/pred_q2iii_v2.csv", index=False)

    fz = json.load(open(f"{S.EXT}/results/frozen_operating_point.json"))
    T_core = float(fz["T_core"])
    old_net = [s for s in SI1.site if s in S.old7_sites()]

    def es_of(b, sites, SI):
        out = []
        for si in sites:
            g = SI[SI.site == si].iloc[0]
            m0 = g.members.split("+")[0]
            kl = ("borehole" if m0.startswith("B") and m0[1:].isdigit() else "surface_bb")
            out.append(M["e_class"](float(S.epi_km(F.lat[b], F.lon[b], g.lat, g.lon)), kl))
        return out

    dense1 = S.dense_sites(SI1)
    dense2 = [s for s in SI2[SI2.geom_le60 & ~SI2.withheld].site]
    rrows = []
    for b in F.index:
        r_old = M["recall"](T_core, M["bg"] + M["e_net"](es_of(b, old_net, SI1)))
        r_v1 = M["recall"](T_core, M["bg"] + M["e_net"](es_of(b, sorted(set(old_net) | set(dense1)), SI1)))
        r_v2 = M["recall"](T_core, M["bg"] + M["e_net"](es_of(b, sorted(set(old_net) | set(dense2)), SI2)))
        rrows.append(dict(bkey=int(b), recall_old=r_old, recall_V1=r_v1, recall_V2=r_v2,
                          gain_V1=(r_v1 - r_old) / max(r_old, 1e-9),
                          gain_V2=(r_v2 - r_old) / max(r_old, 1e-9),
                          d45_V1=int(PR.set_index("bkey").d45_V1[int(b)]),
                          d45_V2=int(PR.set_index("bkey").d45_V2[int(b)])))
    RC = pd.DataFrame(rrows)
    RC.to_csv(f"{OUT}/pred_q2iii_recall_v2.csv", index=False)
    sel1 = RC[RC.d45_V1 >= 1]
    sel2 = RC[RC.d45_V2 >= 1]
    g1 = float((sel1.recall_V1.sum() - sel1.recall_old.sum()) / max(sel1.recall_old.sum(), 1e-9))
    g2 = float((sel2.recall_V2.sum() - sel2.recall_old.sum()) / max(sel2.recall_old.sum(), 1e-9))

    dl = delivery_state()
    am = dict(
        created_utc=pd.Timestamp.now("UTC").isoformat(),
        authority="coordinator scope amendment, 2026-09-19; PRE-ANALYSIS, written before any PQ "
                  "transfer stack was opened",
        variants=dict(
            V1=dict(stations=S.STAS_FROZEN + S.NEW_3C, dense_sites=dense1,
                    role="reported alongside; the already-hashed network"),
            V2=dict(stations=S.STAS_FROZEN + S.NEW_3C + [s for s in S.V2_3C],
                    v2_only=[s for s in new_v2], dense_sites=dense2,
                    role="the DECISION network for the Step-3 pre-registered rule")),
        pq_time_split=dict(calibrate={s: S.cal_episodes_for(s) for s in S.V2_3C},
                           test={s: S.test_episodes_for(s) for s in S.V2_3C}),
        delivery_state_at_hash_time=dl,
        premise=dict(file="results/p3d/premise_family_sites_v2.csv",
                     median_n45_before=float(PR.n45_before.median()),
                     median_n45_after_V1=float(PR.n45_after_V1.median()),
                     median_n45_after_V2=float(PR.n45_after_V2.median()),
                     median_n55_after_V1=float(PR.n55_after_V1.median()),
                     median_n55_after_V2=float(PR.n55_after_V2.median()),
                     n_fam_gaining_new_le45_V1=int((PR.d45_V1 >= 1).sum()),
                     n_fam_gaining_new_le45_V2=int((PR.d45_V2 >= 1).sum()),
                     isolated_hidden=PR[PR.isolated][["bkey", "n45_before", "n45_after_V1",
                                                      "n45_after_V2", "new_sites_45_V1",
                                                      "new_sites_45_PQ"]].to_dict("records")),
        prediction=dict(file="results/p3d/pred_q2iii_v2.csv",
                        recall_file="results/p3d/pred_q2iii_recall_v2.csv",
                        PREDICTED_RECALL_GAIN_V1=g1, PREDICTED_RECALL_GAIN_V2=g2,
                        pq_median_pred_snr_N60={s: float(P2[P2.sta == s].pred_snr_N60.median())
                                                for s in new_v2} if len(P2) else {}),
        step3_rule="the pre-registered delta_full >= 1.5 sigma rule is applied to V2; V1 reported",
        tables_sha256=dict(stations_v2=S.sha256_of(ST2.reset_index().round(6).to_dict("records")),
                           sites_v2=S.sha256_of(SI2.round(6).to_dict("records")),
                           premise_v2=S.sha256_of(PR.to_dict("records")),
                           prediction_v2=S.sha256_of(
                               pd.concat([P1, P2], ignore_index=True).round(6).to_dict("records")),
                           prediction_recall_v2=S.sha256_of(RC.round(6).to_dict("records"))))
    am["amendment_sha256"] = S.sha256_of(am)
    am["correction"] = ("2026-09-19b: rebuilt from the downloader's CORRECTED station_epochs.csv "
                        "(StationXML coordinates/elevations, component DIPS -- all four PQ stations "
                        "declare HHZ dip +90, positive DOWN -- and the MEASURED data extents). "
                        "amendment_v2 is kept verbatim; this is amendment_v2b.")
    am["z_sign"] = {s: S.ZSIGN.get(s, 1.0) for s in S.STAS_P3D if s in S.ZSIGN}
    am["measured_extents"] = {s: str(v) for s, v in
                              pd.read_csv(f"{S.P3D}/station_epochs.csv").set_index("sta")
                              .first_data_month.dropna().items()}
    fr = json.load(open(f"{OUT}/freeze.json"))
    key = "amendment_v2" if "amendment_v2" not in fr else "amendment_v2b"
    fr[key] = am
    json.dump(fr, open(f"{OUT}/freeze.json", "w"), indent=1, default=str)

    print(f"[v2] V2-only sites: {sorted(new_v2)}")
    print(f"[v2] premise: median sites <=45 km  before {PR.n45_before.median():.0f} -> "
          f"V1 {PR.n45_after_V1.median():.0f} -> V2 {PR.n45_after_V2.median():.0f}; "
          f"<=55 km {PR.n55_before.median():.0f} -> {PR.n55_after_V1.median():.0f} -> "
          f"{PR.n55_after_V2.median():.0f}")
    print(f"[v2] families gaining >=1 genuinely new site <=45 km: V1 {(PR.d45_V1>=1).sum()}/130, "
          f"V2 {(PR.d45_V2>=1).sum()}/130")
    print(f"[v2] dense (<=60 km rule) V1 {dense1}  ->  V2 {dense2}")
    print(f"[v2] PREDICTED recall gain: V1 {g1*100:+.1f} %, V2 {g2*100:+.1f} %")
    print(f"[v2] PQ delivery at hash time: {json.dumps(dl)}")
    print(f"[v2] amendment sha256 {am['amendment_sha256']}")
    print(PR[PR.isolated][["bkey", "n45_before", "n45_after_V1", "n45_after_V2",
                           "new_sites_45_V1", "new_sites_45_PQ"]].to_string(index=False))


if __name__ == "__main__":
    main()
