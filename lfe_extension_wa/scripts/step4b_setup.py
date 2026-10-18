#!/usr/bin/env python
"""Step 4b stage 0: ADD STATIONS so the NW and Olympic families can reach >= 3 permanent-era sites.

Candidates = every station of data/extension_stations/{NW_WA,OLY2010}_stations.csv that
  * has an HH/BH/EH epoch at >= 40 sps overlapping one of the step-4b windows (V2010, A2011, T2012, T2013,
    T2014: 2010-2014, inside the requested 2011-2016 span + the 2010 validation window),
  * has VERIFIED data: a 10-min sample returned data (sample_availability.csv status ok / ok_retest) in a
    year of a window it is live in,
  * lies <= 75 km from at least one NW / Olympic family (pos_lat/lon; approx_lat/lon if unlocated),
  * is not one of the 42 trial stations.
Excluded (a column `excluded`, rows kept): the discovery-era arrays XU XG ZK Z3 YW ZU (their stations ARE
the array era; XU ends 2012-05), PH5-only ZZ YB 3F, GS (no data before 2024), PO (data end 2007),
stations < 3 km from a trial site (the 3-km site rule merges them into that site: no new site, e.g. the
C8 PA array at PGC).  TA and 7D are tier 2: used only for families that the tier-1 set does not
bring to P(>= 3 sites) >= 0.5 (request: "TA/7D allowed if they are the only option").

Projection (a planning number, never a gate): P(clear template | distance, n_comp) = the empirical
clear-arrival rate of the 42-station build templates with >= 30 events (tpl_qc_build.csv) in distance
bins; a new site is templatable for a family only if it is live on days with that family's source
detections: NW = V2010 build half or A2011 array detections; Olympic = V2010 array detections;
later-only sites (live only 2012-06..2014) only for families projected to separate (P >= 0.5 of >= 3
sites from array-era sites), through the bridge rule.  P(>= 3 sites) = Poisson-binomial.
Selection: candidate sites ranked by their expected clear templates over families below 3 sites (the
rank is a column); ALL eligible tier-1 sites are taken (the planned 95 %-saturation stop ended after a
few sites because P(clear) is low beyond 20 km, and taking all costs little); tier 2 (TA/7D) sites are
ranked on the families the tier-1 set leaves below P(>= 3) >= 0.5 and taken when they add >= 1
expected clear template to those families and are live in a scored window (V2010/T2012/T2013/T2014).

  -> results/step4b/setup/{families,stations,station_days,truth,scan_days}.csv|parquet
  -> results/step4b/setup/extra_nodes.csv       node epochs of the added stations (lenient bridge table)
  -> results/step4b/candidates.csv              every candidate with every column (distance, liveness,
                                                verification, P, rank, tier, selected, excluded)
  -> results/step4b/projection_by_family.csv    per family: sites before, expected added, P(>= 3)
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

OUT = f"{W.RES}/step4b"
SET = f"{OUT}/setup"
os.makedirs(SET, exist_ok=True)
WINS = ["V2010", "A2011", "T2012", "T2013", "T2014"]
REG4B = ["Royer_NW_WA", "Olympic_2010", "Olympic_2010_armB"]
EXCL_NETS = {"XU": "discovery-era array (ends 2012-05)", "XG": "discovery-era array", "ZK": "discovery-era array",
             "Z3": "discovery-era array", "YW": "discovery-era array", "ZU": "short 2011 campaign",
             "ZZ": "PH5-only", "YB": "PH5-only", "3F": "PH5-only", "GS": "no data before 2024",
             "PO": "data end 2007-08"}
TIER2 = {"TA", "7D"}
MAXKM = 75.0

OLD = pd.read_csv(f"{W.RES}/setup/stations.csv")
OLDF = pd.read_csv(f"{W.RES}/setup/families.csv")
Fall = pd.read_csv(f"{W.ROOT}/data/catalogue/families_v1_1.csv")
Fall = Fall[Fall.region.isin(W.REGIONS) & (Fall.representative != False)].copy()   # noqa: E712
assert len(Fall) == len(OLDF)
F = OLDF[OLDF.region.isin(REG4B)].copy()
ap = Fall.set_index("family_uid")
key = F.family_uid
F["pos_from"] = np.where(F.lat_use.notna(), "pos", "approx")
F.loc[F.lat_use.isna(), "lat_use"] = key[F.lat_use.isna()].map(ap.approx_lat)
F.loc[F.lon_use.isna(), "lon_use"] = key[F.lon_use.isna()].map(ap.approx_lon)
F["val_windows"] = "V2010"
F["flag"] = np.where(F.family_uid == "OLY:938", "runaway ~200 det/day in data AND scrambled control (B013/B943 pair); kept, flagged", "")
F.to_csv(f"{SET}/families.csv", index=False)
print("families", len(F), F.region.value_counts().to_dict(), "positions from approx:", (F.pos_from == "approx").sum())

# ---------------------------------------------------------------- window day spans
wspan = {w: (min(d for d, _ in W.window_days(w)), max(d for d, _ in W.window_days(w))) for w in WINS}
wyear = {w: int(wspan[w][0][:4]) for w in WINS}

# ---------------------------------------------------------------- candidate epochs
ST = pd.concat([pd.read_csv(f"{W.EXT}/{r}_stations.csv") for r in ["NW_WA", "OLY2010"]]).drop_duplicates(
    ["net", "sta", "loc", "band", "fs", "start", "end"])
ST["ns"] = ST.net + "." + ST.sta
ST = ST[ST.band.isin(["HH", "BH", "EH"]) & (ST.fs >= 40)]
old_ns = set(OLD.netsta)
ST = ST[~ST.ns.isin(old_ns)]
AV = pd.read_csv(f"{W.EXT}/sample_availability.csv")
AV["ns"] = AV.net + "." + AV.sta
okyr = AV[AV.status.isin(["ok", "ok_retest"])].groupby("ns").year.apply(set).to_dict()

rows = []
for ns, g in ST.groupby("ns"):
    net, sta = ns.split(".")
    live = {w: bool(((g.start <= wspan[w][1]) & (g.end > wspan[w][0])).any()) for w in WINS}
    ver = {w: live[w] and (wyear[w] in okyr.get(ns, set())) for w in WINS}
    r0 = g.iloc[0]
    rows.append(dict(netsta=ns, net=net, sta=sta, lat=r0.lat, lon=r0.lon, n_comp=int(g.n_comp.max()),
                     bands=",".join(sorted(set(g.band))), fs=",".join(f"{x:g}" for x in sorted(set(g.fs))),
                     start=g.start.min(), end=g.end.max(),
                     **{f"live_{w}": live[w] for w in WINS}, **{f"ver_{w}": ver[w] for w in WINS}))
C = pd.DataFrame(rows)
C["n_win_verified"] = C[[f"ver_{w}" for w in WINS]].sum(axis=1)
d = W.km(C.lat.to_numpy()[:, None], C.lon.to_numpy()[:, None], F.lat_use.to_numpy()[None], F.lon_use.to_numpy()[None])
d = np.where(np.isfinite(d), d, np.inf)
C["dmin_km"] = d.min(1)
C["n_fam_20km"] = (d <= 20).sum(1); C["n_fam_30km"] = (d <= 30).sum(1); C["n_fam_40km"] = (d <= 40).sum(1)
C["dmed_km"] = np.median(d, 1)
dold = W.km(C.lat.to_numpy()[:, None], C.lon.to_numpy()[:, None], OLD.lat.to_numpy()[None], OLD.lon.to_numpy()[None])
C["nearest_trial_site"] = OLD.site.to_numpy()[dold.argmin(1)]
C["km_to_trial_site"] = dold.min(1)
ex = []
for r in C.itertuples():
    e = EXCL_NETS.get(r.net, "")
    if not e and r.km_to_trial_site < W.SITE_KM:
        e = f"< 3 km from trial site {r.nearest_trial_site} (merges, no new site)"
    if not e and r.dmin_km > MAXKM:
        e = "> 75 km from every NW/Olympic family"
    if not e and r.n_win_verified == 0:
        e = "no verified data in a live step-4b window" if any(getattr(r, f"live_{w}") for w in WINS) else \
            "not live in a step-4b window (2010-2014)"
    ex.append(e)
C["excluded"] = ex
C["tier"] = np.where(C.net.isin(TIER2), 2, 1)

# ---------------------------------------------------------------- empirical P(clear | d, n_comp)
Q = pd.read_csv(f"{W.RES}/templates/tpl_qc_build.csv")
Q = Q[Q.n_events >= W.MIN_EVENTS].merge(OLDF[["fid", "lat_use", "lon_use"]], on="fid").merge(
    OLD[["netsta", "lat", "lon"]], on="netsta")
Q = Q[Q.lat_use.notna()]
Q["d"] = W.km(Q.lat_use, Q.lon_use, Q.lat, Q.lon)
BINS = [0, 5, 10, 15, 20, 30, 40, 60, 100, 1e9]
Q["b"] = pd.cut(Q.d, BINS, labels=False)
Q["c3"] = (Q.nch >= 3).astype(int)
PT = Q.assign(cl=Q.clear.astype(bool)).groupby(["c3", "b"]).cl.mean().unstack().reindex(columns=range(len(BINS) - 1)).fillna(0.0)
PT.to_csv(f"{OUT}/p_clear_table.csv")
print("P(clear | d-bin, 3c):\n", PT.round(3).to_string())


def pclear(dist, c3):
    b = np.clip(np.searchsorted(BINS, dist, side="right") - 1, 0, len(BINS) - 2)
    return PT.loc[int(c3)].to_numpy()[b]


# current clear sites per family (build + t2012new templates of the trial)
QB = pd.concat([pd.read_csv(f"{W.RES}/templates/tpl_qc_{t}.csv") for t in ["build", "t2012new"]])
QB = QB[QB.clear.fillna(False).astype(bool)]
cur = QB.groupby("fid").site.nunique()
F["sites_before"] = F.fid.map(cur).fillna(0).astype(int)

# ---------------------------------------------------------------- candidate SITES (3-km merge among candidates)
K = C[C.excluded == ""].reset_index(drop=True)
lab = list(range(len(K)))
for i in range(len(K)):
    for j in range(i + 1, len(K)):
        if W.km(K.lat[i], K.lon[i], K.lat[j], K.lon[j]) < W.SITE_KM:
            a, b = lab[i], lab[j]; lab = [a if x == b else x for x in lab]
K["cand_site"] = [K.sta[[k for k in range(len(K)) if lab[k] == lab[i]]].str.cat(sep="+") for i in range(len(K))]
# per (family, candidate station) P; templatable by era
isNW = (F.region == "Royer_NW_WA").to_numpy()
arr_live = np.where(isNW[None], (K.ver_V2010 | K.ver_A2011).to_numpy()[:, None], K.ver_V2010.to_numpy()[:, None])
later_live = (K.ver_T2012 | K.ver_T2013 | K.ver_T2014).to_numpy()[:, None] & ~arr_live
dK = W.km(K.lat.to_numpy()[:, None], K.lon.to_numpy()[:, None], F.lat_use.to_numpy()[None], F.lon_use.to_numpy()[None])
dK = np.where(np.isfinite(dK), dK, np.inf)
P = np.stack([pclear(dK[i], K.n_comp[i] >= 3) for i in range(len(K))])
P_arr = P * arr_live
P_late = P * later_live


def site_p(Pm, idx):
    """P(site clear) per candidate site = 1 - prod(1 - p) over its stations; idx = station rows used"""
    out = {}
    for i in idx:
        s = K.cand_site[i]
        out.setdefault(s, np.ones(len(F)))
        out[s] = out[s] * (1 - Pm[i])
    return {s: 1 - v for s, v in out.items()}


def p_ge3(base, plist):
    """Poisson-binomial P(base + sum Bernoulli(p) >= 3) per family"""
    need = np.clip(3 - base, 0, 3)
    dist = np.zeros((len(base), 4)); dist[:, 0] = 1
    for p in plist:
        nd = dist * (1 - p)[:, None]
        nd[:, 1:] += dist[:, :-1] * p[:, None]
        nd[:, 3] += dist[:, 3] * p
        dist = nd
    cum = np.cumsum(dist[:, ::-1], 1)[:, ::-1]   # P(X >= k)
    return np.where(need == 0, 1.0, cum[np.arange(len(base)), need.astype(int)])


base = F.sites_before.to_numpy()


def project(sel_idx):
    SA = site_p(P_arr, sel_idx)
    p3a = p_ge3(base, list(SA.values()))
    SL = site_p(P_late * (p3a >= 0.5)[None], sel_idx)
    exp_add = sum(SA.values()) + sum(SL.values()) if SA or SL else np.zeros(len(F))
    return p3a, exp_add, SA


def rank_and_take(pool, base_sel):
    below = base < 3
    score = {}
    for i in pool:
        score.setdefault(K.cand_site[i], 0.0)
        score[K.cand_site[i]] += float((np.maximum(P_arr[i], P_late[i]) * below).sum())
    order = sorted(score, key=lambda s: -score[s])
    full = project(base_sel + pool)[0]
    target = 0.95 * (full >= 0.5).sum()
    take = []
    for s in order:
        take += [i for i in pool if K.cand_site[i] == s]
        if (project(base_sel + take)[0] >= 0.5).sum() >= target and score[s] >= 0:
            break
    return take, order, score


pool1 = [i for i in range(len(K)) if K.tier[i] == 1]
_, order1, score1 = rank_and_take(pool1, [])
take1 = pool1          # ALL eligible tier-1 stations (the 95 %-saturation rule stops after a few sites because
                       # P(clear) is low beyond 20 km; the added volume of taking all is small)
p3_1 = project(take1)[0]
# tier 2 (TA/7D): only for families the tier-1 set leaves below P >= 0.5
low = p3_1 < 0.5
pool2 = [i for i in range(len(K)) if K.tier[i] == 2 and (np.maximum(P_arr[i], P_late[i]) * low).sum() > 0.05]
_, order2, score2 = rank_and_take(pool2, take1) if pool2 else ([], [], {})
# tier-2 rule (the 95 % stop again ends after one site): take a TA/7D site when it adds >= 1 expected clear
# template over the families the tier-1 set leaves below P(>= 3) >= 0.5 AND it is live in a SCORED window
# (7D FN*A live only in A2011 would template but never be scored)
scored = (K.ver_V2010 | K.ver_T2012 | K.ver_T2013 | K.ver_T2014).to_numpy()
take2 = [i for i in pool2 if score2.get(K.cand_site[i], 0) >= 1.0 and scored[i]]
sel = take1 + take2
K["site_score"] = [score1.get(s, score2.get(s, 0.0)) for s in K.cand_site]
K["selected"] = [i in sel for i in range(len(K))]
C = C.merge(K[["netsta", "cand_site", "site_score", "selected"]], on="netsta", how="left")
C["selected"] = C.selected.fillna(False).astype(bool)
C["expected_clear_arrayera"] = C.netsta.map(dict(zip(K.netsta, P_arr.sum(1))))
C["expected_clear_later"] = C.netsta.map(dict(zip(K.netsta, P_late.sum(1))))
C = C.sort_values(["selected", "site_score"], ascending=[False, False])
C.to_csv(f"{OUT}/candidates.csv", index=False)

p3_all, add_all, _ = project(sel)
p3_none = p_ge3(base, [])
F["p_ge3_before"] = p3_none
F["p_ge3_after_tier1"] = p3_1
F["p_ge3_after"] = p3_all
F["expected_sites_added"] = add_all
F["projected_ge3"] = p3_all >= 0.5
F.drop(columns=[]).to_csv(f"{OUT}/projection_by_family.csv", index=False)
F.to_csv(f"{SET}/families.csv", index=False)
print("selected stations", len(sel), "tier1", len(take1), "tier2", len(take2))
print(C[C.selected][["netsta", "tier", "cand_site", "dmin_km", "n_fam_30km", "n_comp", "bands", "start", "end",
                    "site_score"] + [f"ver_{w}" for w in WINS]].to_string())
for rg, g in F.groupby("region"):
    print(rg, "families", len(g), ">=3 sites before", int((g.sites_before >= 3).sum()),
          "projected P>=0.5 tier1", int((g.p_ge3_after_tier1 >= 0.5).sum()), "all", int((g.p_ge3_after >= 0.5).sum()),
          "expected (sum P)", round(float(g.p_ge3_after.sum()), 1))

# ---------------------------------------------------------------- stations = 42 trial + selected
NEW = C[C.selected][["netsta", "net", "sta", "lat", "lon"]].copy()
NEW["why"] = ["step4b tier" + str(t) for t in C[C.selected].tier]
S = pd.concat([OLD[["netsta", "net", "sta", "lat", "lon", "why"]], NEW], ignore_index=True)
site = list(range(len(S)))
for i in range(len(S)):
    for j in range(i + 1, len(S)):
        if W.km(S.lat[i], S.lon[i], S.lat[j], S.lon[j]) < W.SITE_KM:
            a, b = site[i], site[j]; site = [a if x == b else x for x in site]
S["site_ix"] = site
S["site"] = S.groupby("site_ix").sta.transform(lambda v: "+".join(sorted(v)))
S = S.drop(columns="site_ix")
chk = S.iloc[:len(OLD)].merge(OLD[["netsta", "site"]], on="netsta", suffixes=("", "_old"))
assert (chk.site == chk.site_old).all(), "adding stations changed a trial site"
S.to_csv(f"{SET}/stations.csv", index=False)
print("stations", len(S), "sites", S.site.nunique())

# node epochs of the added stations: the lenient bridge table (merge rule of 05_bridge.py)
LN = pd.concat([pd.read_csv(f"{W.EXT}/{r}_bridge_nodes_lenient.csv") for r in ["NW_WA", "OLY2010"]]).drop_duplicates("node")
LN = LN[(LN.net + "." + LN.sta).isin(set(NEW.netsta))]
LN.to_csv(f"{SET}/extra_nodes.csv", index=False)
missing = set(NEW.netsta) - set(LN.net + "." + LN.sta)
print("extra nodes", len(LN), "stations without a node row:", sorted(missing))

# ---------------------------------------------------------------- station-days
SN = pd.concat([pd.read_csv(f"{W.EXT}/{r}_bridge_nodes_strict.csv") for r in ["NW_WA", "SW_WA", "OLY2010"]]).drop_duplicates("node")
SN = pd.concat([SN[~(SN.net + "." + SN.sta).isin(set(NEW.netsta))], LN]).drop_duplicates("node")
Qb = pd.read_csv(f"{W.RES}/templates/tpl_qc_build.csv")
build_nodes = set(Qb.node)
rows = []
for w in WINS:
    for dd, k in W.window_days(w):
        for r in S.itertuples():
            n = SN[(SN.net == r.net) & (SN.sta == r.sta) & (SN.start <= dd) & (SN.end > dd)]
            need = True
            if w == "A2011" and r.netsta in old_ns:
                # A2011 is templating only: a trial station is needed only on an epoch with no build stack
                need = len(n) > 0 and not set(n.node) <= build_nodes
            rows.append(dict(window=w, date=dd, kind=k, netsta=r.netsta, meta_live=bool(len(n) > 0 and need),
                             meta_nodes=";".join(n.node), raw=W.raw_path(r.netsta, dd) or ""))
SD = pd.DataFrame(rows)
SD.to_csv(f"{SET}/station_days.csv", index=False)
lv = SD[SD.meta_live]
print("station-days", len(SD), "live", len(lv), "on disk", int((lv.raw != "").sum()), "to fetch", int((lv.raw == "").sum()))
print(lv.assign(new=~lv.netsta.isin(old_ns)).groupby(["window", "new"]).size().unstack().to_string())

# ---------------------------------------------------------------- truth (V2010 + A2011 array detections)
days = {dd for w in ["V2010", "A2011"] for dd, _ in W.window_days(w)}
T, SC = [], []
for src, g in F.groupby("det_source"):
    D = pd.read_parquet(f"{W.ROOT}/{src}", columns=["grp", "t", "date"])
    D["date"] = D.date.astype(str)
    SC.append(pd.DataFrame(dict(det_source=src, date=sorted(set(D.date)))))
    D = D[D.date.isin(days)]
    m = dict(zip(g.family_id, g.fid))
    D = D[D.grp.isin(m)].copy()
    D["fid"] = D.grp.map(m)
    D["det_source"] = src
    T.append(D[["fid", "t", "date", "det_source"]])
T = pd.concat(T, ignore_index=True)
T["half"] = W.half_of(T.t.to_numpy())
T["window"] = np.where(T.date < "2011-01-01", "V2010", "A2011")
T.to_parquet(f"{SET}/truth.parquet", index=False)
pd.concat(SC).to_csv(f"{SET}/scan_days.csv", index=False)
print("truth rows", T.groupby(["window", "half"]).size().to_dict())
