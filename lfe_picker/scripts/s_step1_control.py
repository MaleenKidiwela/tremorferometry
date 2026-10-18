#!/usr/bin/env python
"""P3D STEP 1 -- co-located MACHINERY positive controls (Merlin: run FIRST, stop if it fails).

CLRS is 25 m from B926 and OSQM is 17 m from B001: the same SITE, so a transfer stack built at the
new station on 2019-2022 hp events must reproduce the co-located borehole's 2012 template.  The
borehole templates are Z/EH1/EH2 (NOT rotated) and the new stations are ZNE, so the frozen template
is ROTATED to ZNE first (Merlin flaw 6) and re-normalised per component.

Statistic: `g_common.cc_template` -- exactly the statistic whose per-station cross-family p95 is
tabulated in results/gates/g_crossfam_cc.csv (B926 0.348, B001 0.219), measured against the SAME
events at the borehole (the control leg) and at the new station.

control-OK  = the borehole's OWN 2019-2022 leg on those events reaches its own cross-family p95
              (i.e. the test can say something at all: the year and the event population have not
              already destroyed the identity)
match       = the new station's leg reaches the same bar
BAR: >= 70 % of control-OK families match, else the machinery is broken -> STOP.
"""
from __future__ import annotations
import os, sys, json, argparse, warnings
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                        # noqa: E402
import s_transfer as T                                                      # noqa: E402
import g_common as G                                                        # noqa: E402

PAIRS = [("CLRS", "B926", 0.34804172188397436), ("OSQM", "B001", 0.2189784058758827)]


def rotate_tpl(tpl, sta):
    """Borehole 12-s template (Z, EH1, EH2) -> ZNE, re-normalised per component."""
    t = np.asarray(tpl, float).copy()
    t[1:3] = S.rot_matrix(sta) @ t[1:3]
    n = np.linalg.norm(t, axis=1, keepdims=True)
    return t / np.where(n > 0, n, 1.0)


def events_for(H, HF, b, alive_days, min_n=30, min_n_fb=60):
    e = H[(H.bkey == b)]
    e = e[[(pd.Timestamp(t, unit="s").year, pd.Timestamp(t, unit="s").dayofyear) in alive_days
           for t in e.t]]
    if len(e) >= min_n:
        return e.t.to_numpy(), "hp", len(e)
    f = HF[(HF.bkey == b)]
    f = f[[(pd.Timestamp(t, unit="s").year, pd.Timestamp(t, unit="s").dayofyear) in alive_days
           for t in f.t]]
    if len(f) >= min_n_fb:
        return f.t.to_numpy(), "fallback", len(f)
    return None, "insufficient", max(len(e), len(f))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-events", type=int, default=400)
    ap.add_argument("--out", default=f"{S.P3D}/step1_colocated.csv")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    ST = S.station_table()
    F = S.families()
    Z = np.load(f"{S.IN}/e5r_templates_2012.npz", allow_pickle=True)
    TPL = Z["tpl"].item()
    QC = pd.read_csv(f"{S.IN}/e5r_template_qc_2012.csv")
    QC = QC[QC.usable & QC.template]
    H = S.hp_events(0.8, episodes=S.CAL_EPISODES)
    HF = S.hp_events(0.5, episodes=S.CAL_EPISODES, column="cat_core")
    rng = np.random.default_rng(T.RNG_SEED)
    rows = []
    for new, ref, bar in PAIRS:
        tts = S.fam_tt([new, ref], ST)
        idx = S.cache_index(years=list(range(2019, 2027)), stas=[new, ref])
        days_new = {(y, d) for (s, y, d) in idx if s == new}
        days_ref = {(y, d) for (s, y, d) in idx if s == ref}
        both = days_new & days_ref
        dc = S.DayCache(maxn=1400, years=list(range(2019, 2027)), stas=[new, ref])
        fams = sorted(int(b) for b in QC[QC.sta == ref].bkey.unique())
        print(f"[step1] {new} vs {ref}: {len(fams)} families with a {ref} 2012 template, "
              f"{len(both)} shared station-days", flush=True)
        for b in fams:
            key = (float(b), ref)
            if key not in TPL or 12.0 not in TPL[key]:
                continue
            t_ev, kind, n_av = events_for(H, HF, b, both)
            if t_ev is None:
                rows.append(dict(new=new, ref=ref, bkey=b, kind="insufficient", n_events=n_av))
                continue
            if len(t_ev) > a.max_events:
                t_ev = np.sort(rng.choice(t_ev, a.max_events, replace=False))
            tpl_z = rotate_tpl(TPL[key][12.0], ref)
            out = dict(new=new, ref=ref, bkey=b, kind=kind, n_events=len(t_ev))
            for who in (ref, new):
                L = T.leg(dc, who, b, t_ev, tts[(b, who)], ST, rng, want_nulls=False)
                if L.get("status") == "no_data":
                    out[f"cc_{'ctl' if who == ref else 'new'}"] = np.nan
                    continue
                idx_s = T.NOM + int(round(L["s_pick_s"] * S.FS))
                cc = G.cc_template(L["stack_pm20"].astype(float), tpl_z, idx_s=idx_s,
                                   pre_tpl=1.0, lag_s=1.0)
                tag = "ctl" if who == ref else "new"
                out[f"cc_{tag}"] = cc
                out[f"n_{tag}"] = L.get("n_stacked", 0)
                out[f"snr_{tag}"] = L.get("snr_db", np.nan)
                out[f"rho_{tag}"] = L.get("rho_loo", np.nan)
                out[f"stack_{tag}"] = L["stack_pm20"]
                out[f"sidx_{tag}"] = idx_s
            rows.append(out)
        # chance level of the SAME statistic at the NEW station: its stacks vs OTHER families'
        # 2012 borehole templates (the honest null for the bar).
        done = [r for r in rows if r["new"] == new and "stack_new" in r]
        cross = []
        for r in done:
            for r2 in done:
                if r2["bkey"] == r["bkey"]:
                    continue
                k2 = (float(r2["bkey"]), ref)
                cross.append(abs(G.cc_template(r["stack_new"].astype(float), rotate_tpl(TPL[k2][12.0], ref),
                                               idx_s=r["sidx_new"], pre_tpl=1.0, lag_s=1.0)))
        if cross:
            print(f"[step1] {new} cross-family chance: n={len(cross)} med {np.median(cross):.3f} "
                  f"p95 {np.quantile(cross, .95):.3f}", flush=True)
            for r in rows:
                if r["new"] == new:
                    r["crossfam_p95_new"] = float(np.quantile(cross, .95))
    D = pd.DataFrame([{k: v for k, v in r.items() if not k.startswith("stack_")} for r in rows])
    D.to_csv(a.out, index=False)
    summ = []
    for new, ref, bar in PAIRS:
        g = D[(D.new == new) & D.cc_ctl.notna()]
        ok = g[g.cc_ctl >= bar]
        m = ok[ok.cc_new >= bar]
        summ.append(dict(new=new, ref=ref, bar=bar, n_tested=len(g), n_control_ok=len(ok),
                         n_match=len(m), frac=len(m) / max(len(ok), 1),
                         cc_ctl_med=float(g.cc_ctl.median()) if len(g) else np.nan,
                         cc_new_med=float(g.cc_new.median()) if len(g) else np.nan,
                         crossfam_p95_new=float(g.crossfam_p95_new.iloc[0])
                         if "crossfam_p95_new" in g and len(g) else np.nan,
                         verdict="PASS" if len(ok) and len(m) / len(ok) >= 0.70 else "FAIL"))
    SU = pd.DataFrame(summ)
    SU.to_csv(f"{S.P3D}/step1_summary.csv", index=False)
    print(SU.to_string(index=False))
    json.dump(summ, open(f"{S.P3D}/step1_summary.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
