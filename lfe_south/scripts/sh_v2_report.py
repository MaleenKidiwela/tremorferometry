#!/usr/bin/env python
"""SH_V2 report: per iteration (tpl_it0..it3) on the held-out TEST days -- pooled recall / precision
at 6..12 x MAD (columns), per-family rows, quiet vs active family-day rates, separation curve, and a
mute check (fraction of his arrivals muted vs the station-day muted fraction).  Appends a dated block
to results/SH_V2_LOG.md between markers (re-written each call)."""
import os, sys, json, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sh_v2_common import *
TH = [6, 7, 8, 9, 10, 11, 12]
LOG = f"{RES}/SH_V2_LOG.md"; B0 = "<!-- SH_V2_RESULTS_BEGIN -->"; B1 = "<!-- SH_V2_RESULTS_END -->"


def block(it):
    out = []
    dg = f"{V2}/mf_it{it}_test_diag.csv"; sc = f"{V2}/mf_it{it}_test_score.csv"
    scf = f"{V2}/mf_it{it}_test_score_fam.csv"; ts = f"{V2}/mf_it{it}_test_truthstat.csv"
    if not os.path.exists(dg): return out
    D = pd.read_csv(dg)
    out.append(f"\n### iteration {it} (templates `sh_v2/tpl_it{it}`), test days\n")
    out.append(f"family-days {len(D)}; |cc|>1 samples {int(D.n_ccgt1.sum())}; family-days with mean cc > 1: "
               f"{int((D.mean_max > 1).sum())}; median live channels {D.live_med.median():.0f}; "
               f"max det/family-day at 9x {int(D.n9.max())}\n")
    q = D.n_theirs == 0
    out.append("| xMAD | quiet det/fam-day (n=%d) | active det/fam-day (n=%d) | Shelly active | max det/fam-day |" % (q.sum(), (~q).sum()))
    out.append("|---|---|---|---|---|")
    for t in TH:
        out.append(f"| {t} | {D[q][f'n{t}'].mean():.2f} | {D[~q][f'n{t}'].mean():.2f} | {D[~q].n_theirs.mean():.2f} | {int(D[f'n{t}'].max())} |")
    for lab, p in [("source", sc), ("fam (own test days only)", scf)]:
        if not os.path.exists(p): continue
        S = pd.read_csv(p)
        out.append(f"\nscore, days-mode {lab} (one-to-one +-2 s):\n")
        out.append("| xMAD | his | ours | matched | recall | precision | fam >=1 match | fam precision >=0.3 |")
        out.append("|---|---|---|---|---|---|---|---|")
        for t in TH:
            g = S[S.tmin == t]
            if not len(g): continue
            th, ou, ma = g.n_theirs.sum(), g.n_ours.sum(), g.n_match.sum()
            out.append(f"| {t} | {th} | {ou} | {ma} | {ma / max(th, 1):.3f} | {ma / max(ou, 1):.3f} | "
                       f"{int((g.n_match > 0).sum())}/{len(g)} | {int((g.precision >= 0.3).sum())} |")
        g = S[S.tmin == 9].sort_values("fam", key=lambda s: s.astype(int))
        out.append(f"\nper family at 9x ({lab}): " + "; ".join(
            f"{r.fam}: R {r.recall_all:.2f} P {r.precision:.2f} ({int(r.n_match)}/{int(r.n_theirs)}/{int(r.n_ours)})" for r in g.itertuples()))
    if os.path.exists(ts):
        T = pd.read_csv(ts)
        out.append(f"\nseparation: statistic at his {len(T)} test-day times (best within +-2 s): p10/p50/p90 "
                   f"{T.nm_pm2s.quantile(.1):.2f}/{T.nm_pm2s.median():.2f}/{T.nm_pm2s.quantile(.9):.2f}; "
                   f"fraction >= 6..12: " + ", ".join(f"{t}:{(T.nm_pm2s >= t).mean():.3f}" for t in TH)
                   + f"; dt median {T.dt_best.median():.2f} s; noise p99.9 median {D.nm_p999.median():.2f}")
    return out


def mute_check():
    p = f"{V2}/prep_summary.csv"
    if not os.path.exists(p): return []
    P = pd.read_csv(p)
    out = [f"\n### prep (mute / lines / dead floors)\n",
           f"station-days {len(P)}; muted fraction mean {P.muted_frac.mean():.4f} (median {P.muted_frac.median():.4f}); "
           f"station-days with a dead channel {int((P.dead_ch > 0).sum())}; dead stations {int(P.station_dead.sum())}; "
           f"lines notched {int(P.n_lines.sum())} (station-days with lines {int((P.n_lines > 0).sum())})",
           "most-muted stations: " + ", ".join(f"{k} {v:.3f}" for k, v in P.groupby('netsta').muted_frac.mean().sort_values(ascending=False).head(8).items()),
           "dead-channel station-days by station: " + ", ".join(f"{k} {v}" for k, v in P[P.dead_ch > 0].groupby('netsta').size().sort_values(ascending=False).head(10).items())]
    # his arrivals muted?
    try:
        M = pd.read_csv(f"{V2}/tpl_it0.csv"); M["fam"] = M.fam.astype(str)
        S = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/external/other/LFEcat_MTJ_2018-2024_shelly2025.csv")
        S.columns = [c.strip('%') for c in S.columns]
        S["t"] = (pd.to_datetime(dict(year=S.year, month=S.month, day=S.day, hour=S.hour, minute=S.minute,
                  second=S.second.astype(int))).astype("datetime64[ns]") + pd.to_timedelta(S.second % 1, unit="s"))
        S["fam"] = S.idnum.astype(str); S["day"] = S.t.dt.floor("D")
        days = sorted(S.day.unique()); rng = np.random.default_rng(0)
        days = [days[i] for i in rng.choice(len(days), min(60, len(days)), replace=False)]
        hit = []; base = []
        for d in days:
            pp = f"{PREP}/{pd.Timestamp(d).date()}.npz"
            if not os.path.exists(pp): continue
            z = np.load(pp); st = set(str(s) for s in z["stas"])
            Sd = S[S.day == d]; mid = pd.Timestamp(d).timestamp()
            for r in M.itertuples():
                if r.netsta not in st: continue
                mu = z[f"{r.netsta}|mute"]
                tt = Sd[Sd.fam == r.fam].t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
                k = ((tt + r.tt_s + r.best_tau - mid) / WIN_S).astype(int); k = k[(k >= 0) & (k < mu.size)]
                hit += list(mu[k]); base.append(mu.mean())
        out.append(f"his S-arrival windows muted: {np.mean(hit):.4f} (n={len(hit)}) vs station-day muted fraction {np.mean(base):.4f} "
                   f"({len(days)} random catalogue days)")
    except Exception as e:
        out.append(f"(mute check failed: {e})")
    return out


def main():
    L = ["", B0, f"## SH_V2 results (auto, sh_v2_report.py, {pd.Timestamp.utcnow():%Y-%m-%d %H:%MZ})"]
    L += mute_check()
    for it in range(4): L += block(it)
    L += [B1, ""]
    s = open(LOG).read() if os.path.exists(LOG) else ""
    if B0 in s:
        s = s[:s.index(B0)].rstrip() + "\n" + s[s.index(B1) + len(B1):].lstrip()
    open(LOG, "w").write(s.rstrip() + "\n" + "\n".join(L))
    log(f"[report] {LOG}")
if __name__ == "__main__": main()
