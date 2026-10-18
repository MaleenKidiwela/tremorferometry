#!/usr/bin/env python
"""WA pre-production step 3 report: Olympic Arm A dominance-cap test (definitions: results/wa_pre/cap/PREREGISTERED.md).

  python scripts/pre_cap_report.py
  -> results/wa_pre/cap/cap_comparison.csv, cap_verdict.json
"""
import os, sys, glob, json
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pre_common as P
W = P.W

OD = f"{P.PRE}/cap"
RUNS = {"frozen": (f"{P.STUDY}/results/wa_gonogo/mf_core21", f"{OD}/score_frozen"),
        "cap30": (f"{OD}/mf_cap30", f"{OD}/score_cap30"),
        "lodo": (f"{OD}/mf_lodo", f"{OD}/score_lodo")}


def groups(F, region):
    f = F[F.region == region].reset_index(drop=True)
    la, lo = f.lat_use.to_numpy(float), f.lon_use.to_numpy(float)
    d = W.km(la[:, None], lo[:, None], la[None], lo[None])
    A = np.nan_to_num(d, nan=99) <= 3.0
    lab = -np.ones(len(f), int); g = 0
    for i in range(len(f)):
        if lab[i] >= 0:
            continue
        st = [i]; lab[i] = g
        while st:
            k = st.pop()
            for j in np.flatnonzero(A[k] & (lab < 0)):
                lab[j] = g; st.append(j)
        g += 1
    return dict(zip(f.fid, lab))


def dedupe_count(t):
    t = np.sort(t); n = 0; last = -1e18
    for x in t:
        if x - last > W.MATCH_S:
            n += 1; last = x
    return n


def counts(mfdir, F, region):
    G = groups(F, region)
    fams = set(G)
    out = {"data": [0, 0], "control": [0, 0], "control2": [0, 0]}
    famdays = 0
    for p in sorted(glob.glob(f"{mfdir}/det_2010-*.parquet")):
        ds = os.path.basename(p)[4:14]
        FD = pd.read_parquet(f"{mfdir}/famday_{ds}.parquet")
        FD = FD[FD.observable & (FD.kind == "episode") & FD.fid.isin(fams)]
        if not len(FD):
            continue
        ok = set(FD.fid); famdays += len(ok)
        D = pd.read_parquet(p)
        D = D[D.fid.isin(ok) & (D.net_mad >= W.T_CORE)]
        D = D[W.half_of(D.t.to_numpy()) == 1]
        for s, g in D.groupby("series"):
            out[s][0] += len(g)
            out[s][1] += sum(dedupe_count(gg.t.to_numpy()) for _, gg in g.groupby(g.fid.map(G)))
    c_raw = (out["control"][0] + out["control2"][0]) / 2; c_grp = (out["control"][1] + out["control2"][1]) / 2
    return dict(famdays_test_episode=famdays, data_cands=out["data"][0], ctl_cands_mean=c_raw,
                ratio_family_level=out["data"][0] / c_raw if c_raw else np.nan,
                data_group_events=out["data"][1], ctl_group_events_mean=c_grp,
                ratio_nbr3_group=out["data"][1] / c_grp if c_grp else np.nan)


def main():
    F = pd.read_csv(f"{P.STUDY}/results/wa_gonogo/setup/families.csv")
    rows = []
    for run, (mfdir, sd) in RUNS.items():
        if not os.path.exists(f"{sd}/summary.csv"):
            print("missing", sd); continue
        S = pd.read_csv(f"{sd}/summary.csv")
        S = S[S.window == "V2010"]
        for rg in ["Royer_NW_WA", "Olympic_2010", "Olympic_2010_armB"]:
            s = S[S.region == rg]
            if not len(s):
                continue
            s = s.iloc[0]
            r = dict(run=run, region=rg, n_fam_observable=s.n_fam_observable, n_truth_test=s.n_truth_test,
                     n_ours_cat_core=s.n_ours_cat_core, prec_own=s.prec_own_cat_core, prec_nbr3=s.prec_nbr3_cat_core,
                     recall_own=s.recall_own_cat_core, n_raw_core=s.n_ours_raw_core, prec_nbr3_raw=s.prec_nbr3_raw_core,
                     n_ctl_core=s.n_ours_ctl_core, prec_own_ctl=s.prec_own_ctl_core, prec_nbr3_ctl=s.prec_nbr3_ctl_core,
                     data_per_famday=s.data_per_famday_core, ctl_per_famday=s.ctl_per_famday_core,
                     sepboth2_famwin=s.sepboth2_core_data, sep2_ctl_vs_ctl=s.sep2_core_ctl)
            r["truth_matched_ratio_nbr3"] = (s.n_ours_raw_core * s.prec_nbr3_raw_core) / (s.n_ours_ctl_core * s.prec_nbr3_ctl_core) \
                if s.n_ours_ctl_core * s.prec_nbr3_ctl_core > 0 else np.nan
            if rg in ("Olympic_2010", "Royer_NW_WA"):
                r.update(counts(mfdir, F, rg))
            rows.append(r)
        # weight shares
        fds = sorted(glob.glob(f"{mfdir}/famday_2010-*.parquet"))
        FDa = pd.concat([pd.read_parquet(p) for p in fds])
        FDa = FDa[FDa.observable].merge(F[["fid", "region"]], on="fid")
        for rg, g in FDa.groupby("region"):
            for r in rows:
                if r["run"] == run and r["region"] == rg:
                    r["share_max_frozen_median"] = float(g.share_max_frozen.median()) if "share_max_frozen" in g else np.nan
                    r["share_max_used_median"] = float(g.share_max_used.median()) if "share_max_used" in g else np.nan
                    r["cap_infeasible_frac"] = float(g.cap_infeasible.fillna(False).astype(bool).mean()) if "cap_infeasible" in g else np.nan
                    r["obs_famdays"] = len(g)
    C = pd.DataFrame(rows)
    C.to_csv(f"{OD}/cap_comparison.csv", index=False)
    print(C.to_string())
    v = {}
    base_nw = C[(C.run == "frozen") & (C.region == "Royer_NW_WA")].prec_own.iloc[0]
    for run in ("cap30", "lodo"):
        o = C[(C.run == run) & (C.region == "Olympic_2010")]
        n = C[(C.run == run) & (C.region == "Royer_NW_WA")]
        if not len(o) or not len(n):
            continue
        ra = float(o.ratio_nbr3_group.iloc[0]); dn = float(n.prec_own.iloc[0] - base_nw)
        v[run] = dict(olyA_ratio_nbr3_group=ra, nw_prec_own=float(n.prec_own.iloc[0]), nw_prec_own_frozen=float(base_nw),
                      nw_delta=dn, pass_a=bool(ra >= 3), pass_b=bool(abs(dn) <= 0.02), adopt=bool(ra >= 3 and abs(dn) <= 0.02))
    fo = C[(C.run == "frozen") & (C.region == "Olympic_2010")]
    v["frozen_olyA_ratio_nbr3_group"] = float(fo.ratio_nbr3_group.iloc[0]) if len(fo) else np.nan
    v["decision"] = "adopt cap30" if v.get("cap30", {}).get("adopt") else "proceed uncapped; Olympic A -> stacks / group rates only"
    json.dump(v, open(f"{OD}/cap_verdict.json", "w"), indent=1)
    print(json.dumps(v, indent=1))


if __name__ == "__main__":
    main()
