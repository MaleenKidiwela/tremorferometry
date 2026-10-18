#!/usr/bin/env python
"""p3_bostock_episodes (2026-10-07) -- P3 step 2: Bostock-era episodes from Bostock's own family exposure with the
final (ladder) episode rule of lfe_extension_wa/results/episode_rates/PLAN.md (find_episodes copied verbatim).
Rules frozen in svi_science/results/p3_bostock/PLAN_P3.md.
    python svi_science/scripts/p3_bostock_episodes.py
"""
import os, json
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/p3_bostock"
B_EP = f"{ROOT}/svi_science/results/b_ets_record/episodes_b.csv"
RULE_FINAL = dict(frac=0.05, gap=5, minlen=5, major=0.15)
RULE_FROZEN = dict(frac=0.10, gap=3, minlen=5, major=0.25)
TARGETS = ["2003-02", "2004-07", "2005-09", "2007-01", "2008-05", "2009-05", "2010-08", "2011-08", "2012-09", "2013-09"]


def find_episodes(Xe, nfam, t0, t1, rule):
    """verbatim lfe_extension_wa/scripts/episode_rates.py find_episodes (n30 = Bostock events that day)."""
    days = pd.date_range(t0, t1, freq="D")
    ds = days.strftime("%Y-%m-%d")
    nexp = Xe[Xe.n30 > 0].groupby("date").fid.nunique().reindex(ds, fill_value=0).to_numpy()
    q = nexp >= rule["frac"] * nfam
    idx = np.flatnonzero(q)
    eps = []
    if len(idx):
        s = idx[0]; prev = idx[0]
        for i in idx[1:]:
            if i - prev - 1 <= rule["gap"]:
                prev = i; continue
            eps.append((s, prev)); s = i; prev = i
        eps.append((s, prev))
    rows = []
    for s, e in eps:
        L = e - s + 1
        if L < rule["minlen"]:
            continue
        fed = int(nexp[s:e + 1].sum())
        rows.append(dict(start=ds[s], end=ds[e], span_days=L, qual_days=int(q[s:e + 1].sum()), family_exposure_days=fed,
                         max_families_exposed=int(nexp[s:e + 1].max()), major_threshold=rule["major"] * nfam * 14,
                         major=fed >= rule["major"] * nfam * 14, kind="bostock"))
    return pd.DataFrame(rows)


def main():
    M = pd.read_parquet(f"{OUT}/bostock_events_mapped.parquet")
    nfam = 115
    X = M.groupby(["fid", "date"]).size().rename("n30").reset_index()
    R = {}
    out = {}
    for tag, rule in (("final", RULE_FINAL), ("frozen", RULE_FROZEN)):
        E = find_episodes(X, nfam, "2003-01-01", "2013-12-31", rule)
        E["episode"] = E.start.str[:7]
        E["episode"] = E.groupby("episode").cumcount().map(lambda k: "" if k == 0 else f"b{k}").radd(E.episode)
        E["n_events"] = [int(((M.date >= a) & (M.date <= b)).sum()) for a, b in zip(E.start, E.end)]
        out[tag] = E
        acc = []
        for t in TARGETS:
            m0 = pd.Timestamp(t + "-01"); m1 = m0 + pd.offsets.MonthEnd(0)
            hit = E[E.major & (pd.to_datetime(E.start) <= m1 + pd.Timedelta(days=15)) & (pd.to_datetime(E.end) >= m0 - pd.Timedelta(days=15))]
            acc.append(dict(target=t, recovered=len(hit) > 0, n_major=len(hit), episodes=";".join(f"{a}..{b}" for a, b in zip(hit.start, hit.end))))
        A = pd.DataFrame(acc)
        R[tag] = dict(rule=rule, n_episodes=int(len(E)), n_major=int(E.major.sum()), n_recovered=int(A.recovered.sum()),
                      splits=A[A.n_major > 1].target.tolist(), PASS=bool(A.recovered.all()))
        A.to_csv(f"{OUT}/step2_acceptance_{tag}.csv", index=False)
        E.to_csv(f"{OUT}/episodes_bostock_{tag}.csv", index=False)
        print(f"--- rule {tag} {rule}\n{E.to_string(index=False)}\n{A.to_string(index=False)}")
    # cross-check vs PNSN-based list 2010-2013
    Bp = pd.read_csv(B_EP)
    Bp = Bp[(Bp.region == "SVI") & (Bp.start >= "2010-01-01") & (Bp.start <= "2013-12-31") & (Bp.major == True)]
    E = out["final"]; E = E[E.major]
    xc = []
    for r in E[E.start >= "2010-01-01"].itertuples():
        h = Bp[(Bp.start <= r.end) & (Bp.end >= r.start)]
        for q in h.itertuples():
            xc.append(dict(bostock_episode=r.episode, b_start=r.start, b_end=r.end, pnsn_episode=q.episode, p_start=q.start,
                           p_end=q.end, start_diff_d=(pd.Timestamp(r.start) - pd.Timestamp(q.start)).days,
                           end_diff_d=(pd.Timestamp(r.end) - pd.Timestamp(q.end)).days))
    XC = pd.DataFrame(xc); XC.to_csv(f"{OUT}/step2_vs_pnsn_2010_2013.csv", index=False)
    print(XC.to_string(index=False))
    R["vs_pnsn_2010_2013"] = XC.to_dict("records")
    R["PASS"] = R["final"]["PASS"]
    json.dump(R, open(f"{OUT}/step2_checks.json", "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in R.items() if k != "vs_pnsn_2010_2013"}, indent=1, default=str))


if __name__ == "__main__":
    main()
