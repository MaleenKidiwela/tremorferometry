#!/usr/bin/env python
"""bw_scope_coverage -- SVI backward scope, step 4: station availability 2003-2013 and family coverage.

Inputs  samples_raw.csv (+ samples_retry.csv), inventory_stations.csv / inventory_channels.csv (this folder),
        svi_science/results/s0/families.csv (115 families, forward template carriers = usable & template).
Rules
  delivered(sample)   any centre returned >= 1 trace for that 10-min request (main pass or serial retry)
  live(station, date) delivered on that sample date
  carrier             forward rule: the family has a usable 2012 template on that station (families.csv carriers)
  site                forward site merge: PGC+B009+B010+B011 = PGC site; LRIV+B005+B006+B007 = LRIV site
  detectable          >= 2 live carrier SITES (the detector's min_sites = 2)
  epoch_same          the station-day is on the same instrument epoch as the 2012 template build
                      (VGZ before 2009-07-31 and SNB before 2011-03-17 are 'epoch differs'; NLLB 2007-01-27..02-27
                      'epoch differs' (one month CMG3T); PGC, NLLB, PB, LRIV unchanged)
Out     availability_station_year.csv (long), availability_matrix.csv (delivered/requested per year),
        coverage_by_date.csv, coverage_by_year.csv, coverage_forward_reference.csv, candidate_new_carriers.csv
"""
import numpy as np, pandas as pd, os

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/backward_scope"
SITES = {"PGC": ["PGC", "B009", "B010", "B011"], "LRIV": ["LRIV", "B005", "B006", "B007"]}
SITE_OF = {s: q for q, m in SITES.items() for s in m}
EPOCH_DIFF = {"VGZ": [("1990-01-01", "2009-07-31")], "SNB": [("1990-01-01", "2011-03-17")],
              "NLLB": [("2007-01-27", "2007-02-28")]}


def hav(lat1, lon1, lat2, lon2):
    r = np.pi / 180
    a = np.sin((lat2 - lat1) * r / 2) ** 2 + np.cos(lat1 * r) * np.cos(lat2 * r) * np.sin((lon2 - lon1) * r / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def epoch_same(sta, d):
    return not any(a <= d < b for a, b in EPOCH_DIFF.get(sta, []))


def main():
    S = pd.read_csv(os.environ.get("BW_SAMPLES", f"{OUT}/samples_raw.csv"), dtype={"http": str})
    S["ok"] = S.http == "200"
    if os.path.exists(f"{OUT}/samples_retry.csv"):
        R = pd.read_csv(f"{OUT}/samples_retry.csv")
        R = R[R.retry_ok]
        k = set(zip(R.centre, R.net, R.sta, R.date))
        S["ok_retry"] = [(c, n, s, d) in k for c, n, s, d in zip(S.centre, S.net, S.sta, S.date)]
    else:
        S["ok_retry"] = False
    S["deliv"] = S.ok | S.ok_retry
    S["year"] = S.date.str[:4].astype(int)
    req = S[S.http != "no_epoch"]
    # per station-date: delivered at any centre
    D = (S.groupby(["net", "sta", "date", "year"]).agg(
        deliv=("deliv", "any"), in_epoch=("http", lambda h: (h != "no_epoch").any()),
        centres_ok=("centre", lambda c: ""), ).reset_index())
    okc = S[S.deliv].groupby(["net", "sta", "date"]).centre.apply(lambda c: ";".join(sorted(set(c))))
    D["centres_ok"] = [okc.get((n, s, d), "") for n, s, d in zip(D.net, D.sta, D.date)]
    # per station-year
    first_pass = req.groupby(["net", "sta", "year", "centre"]).ok.sum().unstack(fill_value=0)
    first_pass.columns = [f"first_pass_ok_{c}" for c in first_pass.columns]
    A = D.groupby(["net", "sta", "year"]).agg(n_dates=("date", "size"), n_in_epoch=("in_epoch", "sum"),
                                              n_delivered=("deliv", "sum")).reset_index()
    cen = D[D.deliv].groupby(["net", "sta", "year"]).centres_ok.apply(lambda c: ";".join(sorted(set(";".join(c).split(";")) - {""})))
    A["centres_delivering"] = [cen.get((n, s, y), "") for n, s, y in zip(A.net, A.sta, A.year)]
    ch = S[S.ok].groupby(["net", "sta", "year"]).agg(chans=("chans", lambda c: ";".join(sorted(set(";".join(c.dropna()).split(";")) - {""}))),
                                                     fs=("fs", lambda c: ";".join(sorted(set(str(x) for x in c.dropna())))))
    A = A.join(ch, on=["net", "sta", "year"]).join(first_pass, on=["net", "sta", "year"])
    I = pd.read_csv(f"{OUT}/inventory_stations.csv")
    A = A.merge(I[["net", "sta", "dist_min_km", "dist_med_km", "template_station"]], on=["net", "sta"], how="left")
    A.to_csv(f"{OUT}/availability_station_year.csv", index=False)
    M = A.assign(v=np.where(A.n_in_epoch > 0, A.n_delivered.astype(str) + "/" + A.n_in_epoch.astype(str), "-")).pivot_table(index=["net", "sta"], columns="year", values="v", aggfunc="first")
    M = M.join(I.set_index(["net", "sta"])[["dist_med_km", "template_station"]]).sort_values(["template_station", "dist_med_km"], ascending=[False, True])
    M.to_csv(f"{OUT}/availability_matrix.csv")

    # ---------------- family coverage per sample date
    F = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv")
    car = {r.fid: [x for x in str(r.carriers).split(";") if x and x != "nan"] for r in F.itertuples()}
    live = D[D.deliv].groupby("date").sta.apply(set).to_dict()
    rows = []
    for d in sorted(D.date.unique()):
        L = live.get(d, set())
        for variant in ("all", "epoch_same"):
            Lv = {s for s in L if variant == "all" or epoch_same(s, d)}
            nst = np.array([len([s for s in car[f] if s in Lv]) for f in F.fid])
            nsi = np.array([len({SITE_OF.get(s, s) for s in car[f] if s in Lv}) for f in F.fid])
            rows.append(dict(date=d, year=int(d[:4]), variant=variant,
                             live_template_stations=";".join(sorted(s for s in Lv if any(s in car[f] for f in F.fid))),
                             n_live_template_stations=len([s for s in Lv if any(s in car[f] for f in F.fid)]),
                             fam_ge1=int((nst >= 1).sum()), fam_ge2=int((nst >= 2).sum()), fam_ge3=int((nst >= 3).sum()),
                             fam_ge2_sites=int((nsi >= 2).sum()), fam_ge3_sites=int((nsi >= 3).sum()),
                             median_carriers=float(np.median(nst)), median_sites=float(np.median(nsi))))
    C = pd.DataFrame(rows)
    C.to_csv(f"{OUT}/coverage_by_date.csv", index=False)
    Y = C.groupby(["year", "variant"]).agg(fam_ge1=("fam_ge1", "median"), fam_ge2=("fam_ge2", "median"), fam_ge3=("fam_ge3", "median"),
                                           fam_ge2_sites=("fam_ge2_sites", "median"), fam_ge3_sites=("fam_ge3_sites", "median"),
                                           median_carriers=("median_carriers", "median"), median_sites=("median_sites", "median"),
                                           n_live_template_stations=("n_live_template_stations", "median"),
                                           min_fam_ge2_sites=("fam_ge2_sites", "min")).reset_index()
    Y.to_csv(f"{OUT}/coverage_by_year.csv", index=False)

    # ---------------- forward reference (the forward station sets, same rule)
    S15 = ["PGC", "VGZ", "NLLB", "SNB", "LRIV", "B001", "B003", "B004", "B005", "B006", "B007", "B009", "B010", "B011", "B926"]
    sets = {"S15 (2014-01..2019-05)": S15, "S14 (2019-05..2021-03, no LRIV)": [s for s in S15 if s != "LRIV"],
            "S13 (2021-04.., no LRIV, SNB)": [s for s in S15 if s not in ("LRIV", "SNB")],
            "S9 (2010 truth run)": ["PGC", "VGZ", "NLLB", "SNB", "LRIV", "B009", "B010", "B011", "B926"]}
    fr = []
    for nm, Ls in sets.items():
        nst = np.array([len([s for s in car[f] if s in Ls]) for f in F.fid])
        nsi = np.array([len({SITE_OF.get(s, s) for s in car[f] if s in Ls}) for f in F.fid])
        fr.append(dict(set=nm, n_stations=len(Ls), fam_ge1=int((nst >= 1).sum()), fam_ge2=int((nst >= 2).sum()), fam_ge3=int((nst >= 3).sum()),
                       fam_ge2_sites=int((nsi >= 2).sum()), fam_ge3_sites=int((nsi >= 3).sum()),
                       median_carriers=float(np.median(nst)), median_sites=float(np.median(nsi))))
    # backward network eras from the FDSN metadata epochs (all stations assumed live; delivery is checked above)
    CN4 = ["PGC", "VGZ", "NLLB", "SNB"]
    OLY = ["B001", "B003", "B004", "B005", "B006", "B007"]
    eras = {"B0 2003-01..2005-06 (CN4 only)": CN4,
            "B1 2005-07..2005-10 (CN4 + Olympic PB)": CN4 + OLY,
            "B2 2005-11..2007-10 (+ B009-B011 at PGC site)": CN4 + OLY + ["B009", "B010", "B011"],
            "B3 2007-11..2009-08 (+ B926 = forward S14 set)": CN4 + OLY + ["B009", "B010", "B011", "B926"],
            "B4 2009-08..2013-12 (+ LRIV = forward S15 set)": S15,
            "B3 minus PB at PGC site (sensitivity)": CN4 + OLY + ["B926"]}
    for nm, Ls in eras.items():
        nst = np.array([len([s for s in car[f] if s in Ls]) for f in F.fid])
        nsi = np.array([len({SITE_OF.get(s, s) for s in car[f] if s in Ls}) for f in F.fid])
        fr.append(dict(set=nm, n_stations=len(Ls), fam_ge1=int((nst >= 1).sum()), fam_ge2=int((nst >= 2).sum()), fam_ge3=int((nst >= 3).sum()),
                       fam_ge2_sites=int((nsi >= 2).sum()), fam_ge3_sites=int((nsi >= 3).sum()),
                       median_carriers=float(np.median(nst)), median_sites=float(np.median(nsi))))
    pd.DataFrame(fr).to_csv(f"{OUT}/coverage_forward_reference.csv", index=False)

    # ---------------- candidate new carriers (no 2012 template), within the radius existing carriers reach
    stn = pd.read_csv(f"{OUT}/inventory_stations.csv").set_index(["net", "sta"])
    st = stn[stn.index.get_level_values(0).isin(["CN", "PB", "UW"])].reset_index().drop_duplicates("sta").set_index("sta")
    dd = [hav(st.loc[s].lat, st.loc[s].lon, r.lat, r.lon) for r in F.itertuples() for s in car[r.fid] if s in st.index]
    r90 = float(np.percentile(dd, 90)); r50 = float(np.percentile(dd, 50))
    cand = []
    for (n, s), g in A.groupby(["net", "sta"]):
        if s in S15:
            continue
        yrs = g[g.n_delivered >= 2].year.tolist()
        if not yrs:
            continue
        d = hav(stn.loc[(n, s)].lat, stn.loc[(n, s)].lon, F.lat.values, F.lon.values)
        cand.append(dict(net=n, sta=s, years_live=";".join(map(str, yrs)), n_years_live=len(yrs),
                         fam_within_r50=int((d <= r50).sum()), fam_within_r90=int((d <= r90).sum()),
                         dist_med_km=round(float(np.median(d)), 1)))
    CA = pd.DataFrame(cand).sort_values("fam_within_r50", ascending=False)
    CA.attrs["r50"] = r50
    CA.to_csv(f"{OUT}/candidate_new_carriers.csv", index=False)
    print(f"carrier distance r50 {r50:.1f} km r90 {r90:.1f} km")
    print(M.to_string()); print(Y.to_string()); print(pd.DataFrame(fr).to_string()); print(CA.to_string())


if __name__ == "__main__":
    main()
