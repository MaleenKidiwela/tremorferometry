#!/usr/bin/env python
"""WA pre-production step 4 report: scrambled-control clean-up and the time-reversed-template control.

Day class (PNSN tremor-window tables of the WA regions, windows within 30 km of a region's families):
  quiet   0 windows that day in every WA region table (NW_WA, SW_WA, OLY2010)
  tremor  >= 1 window in the family's region table that day;  ets = >= 100 windows in the family's region table
Tests fixed before looking (main-session acceptance + the bar chosen here for "~ 0"):
  T1  contaminated share of the scrambled-control candidates (>= 8.75 MAD) on quiet days ~ 0   (pass: <= 0.05)
  T2  the contaminated candidates explain >= 70 % of the control excess on tremor days:
      explained = 1 - (clean_trem - clean_quiet) / (ctl_trem - ctl_quiet), rates per observable family-day
  T3  clean control and time-reversed-template control agree within 1.5x on ETS days (60-day set only)
Rows per run x region x variant; the spec variant ("contam": any family's data candidate at the dominant site) is the
headline, the others are columns.

  python scripts/pre_ctl_report.py  -> results/wa_pre/ctl/ctl_summary.csv, ctl_tests.csv, rev_compare.csv
"""
import os, sys, glob, json
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pre_common as P
W = P.W

OD = f"{P.PRE}/ctl"
RUNS = {"core21_n60": (f"{OD}/res/mf_n60", f"{OD}/contam_n60.parquet"),
        "core21_gonogo": (f"{P.STUDY}/results/wa_gonogo/mf_core21", f"{OD}/contam_core21.parquet"),
        "t1_42sta_mt1all": (f"{P.STUDY}/results/mf_mt1_all", f"{OD}/contam_mt1all.parquet")}
VARS = ["contam", "contam_coh", "contam_other", "contam_origin", "contam_origin_coh", "chance_p97", "chance_m97",
        "chance_coh_p97", "chance_other_p97", "flip_sign"]
REGS = {"NW": ["Royer_NW_WA"], "SW": ["Royer_SW_WA"], "OlympicA": ["Olympic_2010"], "ArmB": ["Olympic_2010_armB"],
        "ALL": ["Royer_NW_WA", "Royer_SW_WA", "Olympic_2010", "Olympic_2010_armB"]}


def day_class(td, region, ds):
    if ds < P.PNSN_START:
        return "no_pnsn"
    allq = all(td[r].get(ds, 0) == 0 for r in td)
    n = td[P.REG_TREM[region]].get(ds, 0)
    if allq:
        return "quiet"
    return "ets" if n >= 100 else ("tremor" if n >= 1 else "other_region_tremor")


def famdays(mfdir, F, td):
    rows = []
    for p in sorted(glob.glob(f"{mfdir}/famday_*.parquet")):
        X = pd.read_parquet(p, columns=["fid", "date", "observable"])
        rows.append(X[X.observable])
    FD = pd.concat(rows, ignore_index=True).merge(F[["fid", "region"]], on="fid")
    FD["cls"] = [day_class(td, r, d) for r, d in zip(FD.region, FD.date)]
    return FD


def main():
    F = pd.read_csv(f"{P.STUDY}/results/wa_gonogo/setup/families.csv")
    td = P.tremor_days()
    out, tests = [], []
    for run, (mfdir, cp) in RUNS.items():
        if not os.path.exists(cp):
            print("missing", cp); continue
        X = pd.read_parquet(cp)
        X = X[X.net_mad >= W.T_CORE].merge(F[["fid", "region"]], on="fid")
        X["cls"] = [day_class(td, r, d) for r, d in zip(X.region, X.date)]
        nser = X.series.nunique()
        FD = famdays(mfdir, F, td)
        days = set(X.date)
        FD = FD[FD.date.isin(days)]
        # data candidates (>= 8.75) per family-day
        dd = []
        for p in sorted(glob.glob(f"{mfdir}/det_*.parquet")):
            if os.path.basename(p)[4:14] not in days:
                continue
            d = pd.read_parquet(p, columns=["fid", "t", "net_mad", "series", "date"])
            dd.append(d[(d.series == "data") & (d.net_mad >= W.T_CORE)])
        DD = pd.concat(dd).merge(F[["fid", "region"]], on="fid")
        DD["cls"] = [day_class(td, r, d) for r, d in zip(DD.region, DD.date)]
        for rk, rl in REGS.items():
            for cls in ["quiet", "tremor", "ets", "no_pnsn"]:
                cl = ["tremor", "ets"] if cls == "tremor" else [cls]
                x = X[X.region.isin(rl) & X.cls.isin(cl)]
                fd = FD[FD.region.isin(rl) & FD.cls.isin(cl)]
                n_fd = len(fd)
                if not n_fd:
                    continue
                r = dict(run=run, region=rk, cls=cls if cls != "tremor" else "tremor(incl.ets)", n_days=fd.date.nunique(), famdays=n_fd,
                         n_ctl_series=nser, data_per_famday=((DD.region.isin(rl) & DD.cls.isin(cl)).sum()) / n_fd,
                         ctl_per_famday=len(x) / nser / n_fd, n_ctl=len(x))
                for v in VARS:
                    if v in x:
                        r[f"share_{v}"] = float(x[v].mean()) if len(x) else np.nan
                        r[f"clean_{v}_per_famday"] = float((~x[v]).sum()) / nser / n_fd
                dom = x.dom_site.value_counts(normalize=True)
                r["top_dom_site"] = dom.index[0] if len(dom) else ""
                r["top_dom_share"] = float(dom.iloc[0]) if len(dom) else np.nan
                out.append(r)
    S = pd.DataFrame(out)
    S.to_csv(f"{OD}/ctl_summary.csv", index=False)
    # tests
    for (run, rk), g in S.groupby(["run", "region"]):
        q = g[g.cls == "quiet"]; t = g[g.cls == "tremor(incl.ets)"]
        if not len(q) or not len(t):
            continue
        q = q.iloc[0]; t = t.iloc[0]
        row = dict(run=run, region=rk, quiet_famdays=q.famdays, tremor_famdays=t.famdays,
                   ctl_quiet=q.ctl_per_famday, ctl_trem=t.ctl_per_famday, data_quiet=q.data_per_famday, data_trem=t.data_per_famday)
        for v in ["contam", "contam_coh", "contam_other", "contam_origin_coh"]:
            cq, ct = q[f"clean_{v}_per_famday"], t[f"clean_{v}_per_famday"]
            ex = t.ctl_per_famday - q.ctl_per_famday
            row[f"{v}_share_quiet"] = q[f"share_{v}"]
            row[f"{v}_share_trem"] = t[f"share_{v}"]
            row[f"{v}_clean_quiet"] = cq; row[f"{v}_clean_trem"] = ct
            row[f"{v}_explained"] = 1 - (ct - cq) / ex if ex > 0 else np.nan
            row[f"{v}_T1_pass"] = bool(q[f"share_{v}"] <= 0.05)
            row[f"{v}_T2_pass"] = bool(ex > 0 and (1 - (ct - cq) / ex) >= 0.70)
        row["chance_share_quiet"] = q.share_chance_p97; row["chance_share_trem"] = t.share_chance_p97
        row["chance_coh_share_quiet"] = q.share_chance_coh_p97
        row["flip_sign_share_quiet"] = q.share_flip_sign; row["flip_sign_share_trem"] = t.share_flip_sign
        row["excess_ctl_per_famday"] = t.ctl_per_famday - q.ctl_per_famday
        tests.append(row)
    TT = pd.DataFrame(tests)
    TT.to_csv(f"{OD}/ctl_tests.csv", index=False)
    # time-reversed control on the 60-day set
    rv = []
    days60 = pd.read_csv(f"{OD}/days60.csv")
    rdir = f"{OD}/res/mf_rev60"
    if os.path.isdir(rdir):
        FDr = famdays(rdir, F, td)
        RR = []
        for p in sorted(glob.glob(f"{rdir}/det_*.parquet")):
            d = pd.read_parquet(p, columns=["fid", "t", "net_mad", "series", "date"])
            RR.append(d[(d.series == "data") & (d.net_mad >= W.T_CORE)])
        RR = pd.concat(RR).merge(F[["fid", "region"]], on="fid")
        cls60 = dict(zip(days60.date, days60.cls))
        RR["cls60"] = RR.date.map(cls60); FDr["cls60"] = FDr.date.map(cls60)
        X = pd.read_parquet(f"{OD}/contam_n60.parquet")
        X = X[X.net_mad >= W.T_CORE].merge(F[["fid", "region"]], on="fid"); X["cls60"] = X.date.map(cls60)
        FDn = famdays(f"{OD}/res/mf_n60", F, td); FDn["cls60"] = FDn.date.map(cls60)
        nser = X.series.nunique()
        DDn = []
        for p in sorted(glob.glob(f"{OD}/res/mf_n60/det_*.parquet")):
            d = pd.read_parquet(p, columns=["fid", "net_mad", "series", "date"])
            DDn.append(d[(d.series == "data") & (d.net_mad >= W.T_CORE)])
        DDn = pd.concat(DDn).merge(F[["fid", "region"]], on="fid"); DDn["cls60"] = DDn.date.map(cls60)
        for rk, rl in REGS.items():
            for c in ["quiet", "moderate", "ets"]:
                nf = int((FDn.region.isin(rl) & (FDn.cls60 == c)).sum()); nfr = int((FDr.region.isin(rl) & (FDr.cls60 == c)).sum())
                if not nf:
                    continue
                x = X[X.region.isin(rl) & (X.cls60 == c)]
                r = dict(region=rk, cls=c, famdays=nf, famdays_rev=nfr,
                         data_per_famday=int((DDn.region.isin(rl) & (DDn.cls60 == c)).sum()) / nf,
                         rev_per_famday=int((RR.region.isin(rl) & (RR.cls60 == c)).sum()) / max(nfr, 1),
                         ctl_per_famday=len(x) / nser / nf)
                for v in ["contam", "contam_coh", "contam_other"]:
                    r[f"clean_{v}_per_famday"] = float((~x[v]).sum()) / nser / nf
                    r[f"rev_over_clean_{v}"] = r["rev_per_famday"] / r[f"clean_{v}_per_famday"] if r[f"clean_{v}_per_famday"] > 0 else np.nan
                    r[f"agree_1p5_{v}"] = bool(np.isfinite(r[f"rev_over_clean_{v}"]) and 1 / 1.5 <= r[f"rev_over_clean_{v}"] <= 1.5)
                r["rev_over_raw_ctl"] = r["rev_per_famday"] / r["ctl_per_famday"] if r["ctl_per_famday"] > 0 else np.nan
                rv.append(r)
    RV = pd.DataFrame(rv)
    RV.to_csv(f"{OD}/rev_compare.csv", index=False)
    pd.set_option("display.width", 250)
    print(S.to_string()); print(TT.T.to_string()); print(RV.to_string())


if __name__ == "__main__":
    main()
