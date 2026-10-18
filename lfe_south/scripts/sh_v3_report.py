#!/usr/bin/env python
"""SH_V3 report: template/pick summary + Text S1 scoring tables (rescore_*.csv) for it0, it1, mix;
written between markers in results/SH_V3_LOG.md (re-written each call)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sh_v3_common import *
LOG = os.environ.get("SH_V3_LOG", f"{RES}/SH_V3_LOG.md")   # v3b: overridable
BASE = os.environ.get("SH_V3_BASE", "")                      # v3b: baseline dir for the v3 vs v3b tables
NAME = os.environ.get("SH_V3_NAME", "SH_V3"); B0 = "<!-- SH_V3_RESULTS_BEGIN -->"; B1 = "<!-- SH_V3_RESULTS_END -->"
TH = [6, 7, 8, 9, 10, 11, 12]


def tpl_block(it):
    p = f"{V3}/tpl_it{it}_all.csv"
    if not os.path.exists(p): return []
    T = pd.read_csv(p); ok = T[T.status == "ok"]
    out = [f"\n### templates it{it}\n",
           f"pairs {len(T)}; ok {len(ok)} over {ok.fam.nunique()} families / {ok.netsta.nunique()} stations; status "
           f"{T.status.value_counts().to_dict()}; median n_stack {ok.n_stack.median():.0f}; median cc50 {ok.cc50.median():.3f}"]
    if "pick_q" in ok and ok.pick_q.notna().any():
        d1 = ok.best_tau + ok.tt_s - ok.stage1_off
        out.append(f"S pick: pick_dt (pick - predicted) median {ok.pick_dt.median():.2f} s, IQR {ok.pick_dt.quantile(.25):.2f}..{ok.pick_dt.quantile(.75):.2f}; "
                   f"pick_q median {ok.pick_q.median():.1f} (q<3: {(ok.pick_q < 3).sum()}); pick - stage-1 offset: median {np.nanmedian(d1):.2f} s, "
                   f"|.|>1.5 s {(np.abs(d1) > 1.5).sum()}/{d1.notna().sum()}, |.|>5 s {(np.abs(d1) > 5).sum()}")
    out.append("per station: " + "; ".join(f"{s} n={len(g)} cc50={g.cc50.median():.2f}" + (f" q={g.pick_q.median():.1f}" if "pick_q" in g else "")
                                             for s, g in ok.groupby("netsta")))
    return out


def score_block(tag):
    p = f"{V3}/rescore_{tag}.csv"
    if not os.path.exists(p): return []
    R = pd.read_csv(p)
    out = [f"\n### {tag}: Text S1 rules, held-out test days (one-to-one +-2 s)\n",
           "| variant | xMAD | mode | ours | his | recall (fam) | precision (fam) | recall (any fam) | precision (any fam) | same-family share of any-fam matches | fam P>=0.3 |",
           "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in R[R.xmad.isin([8, 9, 10, 12]) & R.variant.isin(["v2", "+floor", "S1"])].itertuples():
        out.append(f"| {r.variant} | {r.xmad} | {r.mode} | {r.ours} | {r.his} | {r.recall_fam:.3f} | {r.prec_fam:.3f} | "
                   f"{r.recall_any:.3f} | {r.prec_any:.3f} | {r.same_fam / max(r.match_any, 1):.2f} | {r.fam_p30} |")
    ts = f"{V3}/mf_{tag}_test_truthstat.csv"
    if os.path.exists(ts):
        T = pd.read_csv(ts)
        out.append(f"\nseparation at his {len(T)} test-day times (best within +-2 s): p50 {T.nm_pm2s.median():.2f}; fraction >= 8/9/10: "
                   + "/".join(f"{(T.nm_pm2s >= t).mean():.3f}" for t in (8, 9, 10)))
    return out


def compare_block():
    """v3 (BASE) vs this run: S1 rules, 8/9/10 x MAD, source + fam, family-specific + family-agnostic;
    plus the per-family recall change (S1 9x, both modes)."""
    if not BASE: return []
    out = [f"\n### {NAME} vs baseline ({os.path.relpath(BASE, RES)}): Text S1 rules, held-out test days\n",
           "| iter | xMAD | mode | ours base | ours new | recall fam base | recall fam new | prec fam base | prec fam new | "
           "recall any base | recall any new | prec any base | prec any new |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for tag in ("it0", "it1", "mix"):
        pa, pb = f"{BASE}/rescore_{tag}.csv", f"{V3}/rescore_{tag}.csv"
        if not (os.path.exists(pa) and os.path.exists(pb)): continue
        A = pd.read_csv(pa); B = pd.read_csv(pb)
        for k in (8, 9, 10):
            for m in ("source", "fam"):
                a = A[(A.variant == "S1") & (A.xmad == k) & (A["mode"] == m)]
                b = B[(B.variant == "S1") & (B.xmad == k) & (B["mode"] == m)]
                if not (len(a) and len(b)): continue
                a = a.iloc[0]; b = b.iloc[0]
                out.append(f"| {tag} | {k} | {m} | {a.ours} | {b.ours} | {a.recall_fam:.3f} | {b.recall_fam:.3f} | {a.prec_fam:.3f} | "
                           f"{b.prec_fam:.3f} | {a.recall_any:.3f} | {b.recall_any:.3f} | {a.prec_any:.3f} | {b.prec_any:.3f} |")
    for tag in ("mix", "it0"):
        pa, pb = f"{BASE}/rescore_{tag}_fam.csv", f"{V3}/rescore_{tag}_fam.csv"
        if not (os.path.exists(pa) and os.path.exists(pb)): continue
        A = pd.read_csv(pa, dtype={"fam": str}); B = pd.read_csv(pb, dtype={"fam": str})
        sel = lambda X, m: X[(X.variant == "S1") & (X.xmad == 9) & (X["mode"] == m)].set_index("fam")
        As, Bs, Af, Bf = sel(A, "source"), sel(B, "source"), sel(A, "fam"), sel(B, "fam")
        fams = sorted(set(As.index) | set(Bs.index), key=lambda s: int(s) if s.isdigit() else 999)
        out += [f"\nper-family recall change, {tag}, S1 9x (family-specific one-to-one; his test-day events per family)\n",
                "| fam | his (source) | recall base | recall new | delta | prec base | prec new | his (fam days) | recall base (fam) | recall new (fam) | delta (fam) |",
                "|---|---|---|---|---|---|---|---|---|---|---|"]
        g = lambda X, f, c: float(X.loc[f, c]) if f in X.index else np.nan
        d_all = []
        for f in fams:
            ra, rb = g(As, f, "recall_all"), g(Bs, f, "recall_all"); fa, fb = g(Af, f, "recall_all"), g(Bf, f, "recall_all")
            d_all.append(rb - ra)
            out.append(f"| {f} | {g(Bs, f, 'n_theirs'):.0f} | {ra:.3f} | {rb:.3f} | {rb - ra:+.3f} | {g(As, f, 'precision'):.3f} | "
                       f"{g(Bs, f, 'precision'):.3f} | {g(Bf, f, 'n_theirs'):.0f} | {fa:.3f} | {fb:.3f} | {fb - fa:+.3f} |")
        d = np.array(d_all, float)
        out.append(f"\n{tag}: families with recall up / down / same (source, S1 9x): {(d > 0.005).sum()} / {(d < -0.005).sum()} / "
                   f"{(np.abs(d) <= 0.005).sum()}; median delta {np.nanmedian(d):+.3f}")
    return out


def main():
    L = ["", B0, f"## {NAME} results (auto, sh_v3_report.py, {pd.Timestamp.utcnow():%Y-%m-%d %H:%MZ})"]
    for it in (0, 1): L += tpl_block(it)
    if os.path.exists(f"{V3}/mix_families.csv"):
        m = pd.read_csv(f"{V3}/mix_families.csv")
        L.append(f"\nmix: it1 for the 6 lowest-cc50 families {list(m[m.it == 1].fam)} (build-day stack quality), it0 for the rest")
    for tag in ("it0", "it1", "mix"): L += score_block(tag)
    L += compare_block()
    L += [B1, ""]
    s = open(LOG).read() if os.path.exists(LOG) else ""
    if B0 in s:
        s = s[:s.index(B0)].rstrip() + "\n" + s[s.index(B1) + len(B1):].lstrip()
    open(LOG, "w").write(s.rstrip() + "\n" + "\n".join(L))
    log(f"[report] {LOG}")
if __name__ == "__main__": main()
