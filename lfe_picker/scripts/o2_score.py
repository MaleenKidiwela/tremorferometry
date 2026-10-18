#!/usr/bin/env python
"""P0-O item 2 -- P1a BASELINE: score the FROZEN template configuration on 2011 and 2013.

Scoring is the Track D recipe verbatim (scripts/d1_score_epochs.py): own family, +-2 s, ONE-TO-ONE
greedy by time, HF veto applied, recall read off at the frozen T_core / T_ext.  The time base is the
FROZEN 2012 per-family offsets -- never refit on the referee years (results/referee_protocol.json).

False alarms:
  * quiet family-day = Track C's `quiet_flag` (no PNSN tremor in the family box that day AND >= 3 days
    from any tremor day) and `observable`;
  * FA per template-day = detections on quiet family-days / quiet observable family-days;
  * NETWORK-WIDE false events/day = the same rate x the median number of observable families per day,
    computed on the cross-family-assigned (catalog) set as well as on the raw set.

  python scripts/o2_score.py --year 2011
"""
from __future__ import annotations
import os, sys, json, glob, argparse
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PICK = os.path.dirname(HERE)
RES = f"{PICK}/results"
TRO = f"{RES}/trackO"
IN = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
BOST = f"{IN}/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv"

WINDOWS = {2011: dict(months=["2011-08", "2011-09"], data=("2011-08-05", "2011-09-26"),
                      referee=("2011-08-19", "2011-09-10")),
           2013: dict(months=["2013-08", "2013-09", "2013-10"], data=("2013-08-30", "2013-10-29"),
                      referee=("2013-09-13", "2013-10-14")),
           2012: dict(months=[], data=("2012-08-15", "2012-10-06"), referee=("2012-08-15", "2012-10-06"))}


def one2one(a, b_sorted, tol=2.0):
    a = np.asarray(a, float); out = np.full(len(a), np.inf)
    if len(a) == 0 or len(b_sorted) == 0:
        return out
    j = 0; nb = len(b_sorted)
    for i in np.argsort(a, kind="stable"):
        ti = a[i]
        while j < nb and b_sorted[j] < ti - tol:
            j += 1
        if j < nb and abs(b_sorted[j] - ti) <= tol:
            out[i] = abs(b_sorted[j] - ti); j += 1
    return out


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def gt_year(year, t0, t1):
    """Bostock events in [t0, t1) mapped to bkey by EXACT hypocentre match, with the FROZEN 2012
    per-family offsets applied (t_origin = DB time + offset_s)."""
    FAM = pd.read_csv(f"{IN}/e5r_families_2012.csv")
    OFF = pd.read_csv(f"{IN}/e5r_family_offsets_2012.csv")
    omap = {int(r.bkey): float(r.offset_s) for r in OFF.itertuples() if bool(r.accept)}
    B = pd.read_csv(BOST)
    B["time"] = pd.to_datetime(dict(year=B.year, month=B.month, day=B.day, hour=B.hour,
                                    minute=B["min"], second=B.sec.astype(int)), errors="coerce") \
        + pd.to_timedelta(B.sec % 1.0, unit="s")
    key = {(round(r.lat, 4), round(r.lon, 4), round(r.dep, 4)): int(r.bkey) for r in FAM.itertuples()}
    B["bkey"] = [key.get((round(la, 4), round(lo, 4), round(dp, 4)), -1)
                 for la, lo, dp in zip(B.lat, B.lon, B.dep)]
    Y = B[(B.time.dt.year == year)]
    n_year = len(Y); n_unmapped = int((Y.bkey < 0).sum())
    g = Y[(Y.bkey >= 0) & (Y.time >= t0) & (Y.time < t1)].copy()
    g["t_origin"] = g.time.astype("datetime64[ns]").astype("int64") / 1e9 + g.bkey.map(omap)
    g = g.dropna(subset=["t_origin"]).reset_index(drop=True)
    return g, dict(n_events_in_year=n_year, n_unmapped_families=n_unmapped,
                   n_events_in_data_window=int(len(g)),
                   n_events_outside_data_window=int(n_year - n_unmapped - len(g)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--t-core", type=float, default=8.75)
    ap.add_argument("--t-ext", type=float, default=8.25)
    a = ap.parse_args()
    W = WINDOWS[a.year]
    det = pd.concat([pd.read_parquet(f"{TRO}/detections_{m}.parquet") for m in W["months"]],
                    ignore_index=True)
    fd = pd.concat([pd.read_parquet(f"{TRO}/family_days_{m}.parquet") for m in W["months"]],
                   ignore_index=True)
    det["dt"] = pd.to_datetime(det.date)
    fd["dt"] = pd.to_datetime(fd.date)
    det = det[(det.dt >= W["data"][0]) & (det.dt < W["data"][1])].reset_index(drop=True)
    fd = fd[(fd.dt >= W["data"][0]) & (fd.dt < W["data"][1])].reset_index(drop=True)
    gt, gmeta = gt_year(a.year, W["data"][0], W["data"][1])
    print(f"[gt] {a.year}: {gmeta}")

    out = dict(year=a.year, window=W, gt=gmeta,
               n_detections_recorded=int(len(det)),
               n_hf_vetoed=int((~det.hf_pass).sum()),
               n_family_days=int(len(fd)), n_observable=int(fd.observable.sum()))
    obs = fd[fd.observable]
    nfam_per_day = obs.groupby("date").fam.nunique()
    out["median_observable_families_per_day"] = float(nfam_per_day.median())
    D = det[det.hf_pass].reset_index(drop=True)          # the frozen point adopted the veto

    # ---------------- false alarms on quiet family-days ----------------
    q = obs[obs.quiet_flag]
    out["n_quiet_family_days"] = int(len(q))
    out["n_quiet_calendar_days"] = int(q.date.nunique())
    qkey = set(zip(q.fam.astype(int), q.date))
    dq = D[[ (f, d) in qkey for f, d in zip(D.fam.astype(int), D.date) ]]
    rows_fa = []
    for nm, T, kcol in (("T_ext", a.t_ext, "keep_xfam_ext"), ("T_core", a.t_core, "keep_xfam_core")):
        sub = dq[dq.net_mad >= T]
        fa_tpl = len(sub) / max(len(q), 1)
        cat = sub[sub[kcol]] if kcol in sub.columns else sub
        fa_tpl_cat = len(cat) / max(len(q), 1)
        rows_fa.append(dict(tier=nm, T=T, n_quiet_famdays=int(len(q)),
                            n_det_quiet_raw=int(len(sub)), n_det_quiet_catalog=int(len(cat)),
                            fa_per_template_day_raw=fa_tpl, fa_per_template_day_catalog=fa_tpl_cat,
                            network_fa_per_day_raw=fa_tpl * out["median_observable_families_per_day"],
                            network_fa_per_day_catalog=fa_tpl_cat * out["median_observable_families_per_day"]))
    out["false_alarms"] = rows_fa

    # ---------------- own recall, one-to-one, +-2 s ----------------
    rec_rows, fam_rows = [], []
    for scope, sel in (("data_window", np.ones(len(gt), bool)),
                       ("referee_window", (gt.time >= W["referee"][0]).to_numpy() &
                                          (gt.time < W["referee"][1]).to_numpy())):
        ev = gt[sel]
        for nm, T in (("T_ext", a.t_ext), ("T_core", a.t_core)):
            Dt = D[D.net_mad >= T]
            nh = nt = 0
            for bk, g in ev.groupby("bkey"):
                td = np.sort(Dt[Dt.fam == float(bk)].t.to_numpy())
                hit = one2one(g.t_origin.to_numpy(), td) <= 2.0
                nh += int(hit.sum()); nt += len(g)
                if scope == "data_window":
                    lo, hi = wilson(int(hit.sum()), len(g))
                    fam_rows.append(dict(year=a.year, tier=nm, T=T, family_id=int(bk),
                                         n_events=len(g), n_recalled=int(hit.sum()),
                                         recall=float(hit.mean()), recall_lo=lo, recall_hi=hi,
                                         median_mag=float(np.median(g.mag))))
            lo, hi = wilson(nh, nt)
            rec_rows.append(dict(year=a.year, scope=scope, tier=nm, T=T, n_events=int(nt),
                                 n_recalled=int(nh), recall=nh / max(nt, 1),
                                 recall_lo=lo, recall_hi=hi,
                                 n_families=int(ev.bkey.nunique())))
    out["recall"] = rec_rows
    pd.DataFrame(fam_rows).to_csv(f"{RES}/o2_recall_family_{a.year}.csv", index=False)
    pd.DataFrame(rec_rows).to_csv(f"{RES}/o2_recall_{a.year}.csv", index=False)
    pd.DataFrame(rows_fa).to_csv(f"{RES}/o2_fa_{a.year}.csv", index=False)
    json.dump(out, open(f"{RES}/o2_p1a_{a.year}.json", "w"), indent=1, default=str)
    print(pd.DataFrame(rec_rows).to_string(index=False))
    print(pd.DataFrame(rows_fa).to_string(index=False))


if __name__ == "__main__":
    main()
