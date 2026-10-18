#!/usr/bin/env python
"""Markdown fragments for results/SCREEN_VALIDATION_NORTH.md.  Prefix `q_sv_`."""
from __future__ import annotations
import json, os, sys
import numpy as np, pandas as pd

OUT = "results/q_screenval"
DZ19, MINDAYS, FA25 = 19.0, 3, 0.25
CL = [("c_a", "(a) quiet-day FA &le; 0.1/tpl-day AND tremor rate &ge; 3&times; it"),
      ("c_b", "(b) tremor co-timing > same-day base rate (any-LFE hours, base 0.77)"),
      ("c_c1", "(c1) located on &ge; 5 stations"),
      ("c_c2", "(c2) shuffled-station control passes"),
      ("c_c3", "(c3) INTERIOR to the search grid"),
      ("c_c4", "(c4) &#124;&Delta;z&#124; &le; 8 km of Slab2"),
      ("c_d", "(d) &ge; 3 distinct tremor days"),
      ("c_e", "(e) 2nd-half re-stack relocates within 10 km")]
ALT = [("c_a25", "(a&prime;) quiet-day FA &le; 0.25/tpl-day AND tremor rate &ge; 3&times; it"),
       ("c_a_shuf", "(a&Prime;) tremor- vs non-tremor-day rate, day-label shuffle p &lt; 0.05"),
       ("c_b_loud", "(b&prime;) co-timing > base, base-rate-matched tremor hours (base 0.23)"),
       ("c_b_loud_p", "(b&Prime;) same, at p &lt; 0.05"),
       ("c_c2m", "(c2&prime;) control loosened to rms &lt; shuffled median"),
       ("c_c4_19", "(c4&prime;) &#124;&Delta;z&#124; &le; 19 km of Slab2"),
       ("c_e_hh", "(e&prime;) FIRST-half vs SECOND-half re-stack within 10 km")]


def load():
    V = pd.read_csv(f"{OUT}/screen_north.csv")
    V = V[[c for c in V.columns if not c.startswith("rms_shuf_med")]]
    C = pd.read_csv("results/repro/classify15.csv")[["grp", "rms_shuf_med"]]
    V = V.merge(C, on="grp", how="left")
    V["c_a25"] = ((V.fa_quiet_per_tpl_day <= FA25)
                  & (V.rate_tremor_per_tpl_day >= 3 * V.fa_quiet_per_tpl_day) & (V.n_det > 0))
    V["c_c4_19"] = V.c_c1.astype(bool) & (V.dz.abs() <= DZ19)
    V["c_c2m"] = V.c_c1.astype(bool) & (V.rms < V.rms_shuf_med)
    for c, _ in CL + ALT:
        if c not in V.columns:
            V[c] = False
        V[c] = V[c].fillna(False).astype(bool)
    return V


def row(V, c, name):
    K = V[V.cls == "kept_family"]; H = V[V.cls == "hidden_family"]
    N = V[V.cls == "unlocatable"]; CN = V[V.cls == "candidate_new"]
    kf = K.groupby("bkey")[c].max(); hf = H.groupby("hid_bkey")[c].max()
    return dict(clause=name, kept_fam=kf.mean(), kept_fam_s=f"**{int(kf.sum())}/58 = {kf.mean():.0%}**",
                kept_tpl=f"{int(K[c].sum())}/177 = {K[c].mean():.0%}",
                hid_fam=f"{int(hf.sum())}/18 = {hf.mean():.0%}",
                hid_tpl=f"{int(H[c].sum())}/64 = {H[c].mean():.0%}",
                cand=f"{int(CN[c].sum())}/18",
                noise=f"{int(N[c].sum())}/1201 = {N[c].mean():.0%}",
                noise_rej=1 - N[c].mean())


def md(rows, cols, hdr):
    out = ["| " + " | ".join(hdr) + " |", "|" + "|".join(["---"] * len(hdr)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(out)


def main():
    V = load()
    os.makedirs(OUT, exist_ok=True)
    frag = {}
    R = [row(V, c, n) for c, n in CL]
    frag["clause_main"] = md(R, ["clause", "kept_fam_s", "kept_tpl", "hid_fam", "hid_tpl",
                                 "cand", "noise"],
                             ["clause (as coded)", "TRUE families (58)", "their 177 templates",
                              "hidden families (18)", "their 64 templates",
                              "18 candidate-new", "noise population (1,201)"])
    RA = [row(V, c, n) for c, n in ALT]
    frag["clause_alt"] = md(RA, ["clause", "kept_fam_s", "kept_tpl", "hid_fam", "cand", "noise"],
                            ["reformulated clause", "TRUE families (58)", "their 177 templates",
                             "hidden families (18)", "18 candidate-new",
                             "noise population (1,201)"])

    def verdict(r):
        if r["kept_fam"] < 0.85:
            return f"**UNJUSTIFIED as coded** (rejects {1-r['kept_fam']:.0%} of the true families)"
        if r["noise_rej"] >= 0.50:
            return f"**JUSTIFIED** ({r['kept_fam']:.0%} pass, rejects {r['noise_rej']:.0%} of noise)"
        return (f"passes {r['kept_fam']:.0%} but rejects only {r['noise_rej']:.0%} of the noise "
                f"population &mdash; **retain, no power**")
    frag["verdicts"] = md([dict(clause=r["clause"], v=verdict(r)) for r in R],
                          ["clause", "v"], ["clause", "verdict by the pre-stated rule"])

    # cumulative, in n_screen's own order
    K = V[V.cls == "kept_family"]; H = V[V.cls == "hidden_family"]
    N = V[V.cls == "unlocatable"]; CN = V[V.cls == "candidate_new"]
    cum = pd.Series(True, index=V.index)
    cr = []
    for c, n in [CL[2], CL[3], CL[4], CL[5], CL[6], CL[0], CL[1], CL[7]]:
        cum = cum & V[c]
        kf = int(K.assign(x=cum.loc[K.index].values).groupby("bkey").x.max().sum())
        hf = int(H.assign(x=cum.loc[H.index].values).groupby("hid_bkey").x.max().sum())
        cr.append(dict(clause=n, kf=f"**{kf}/58**", kt=int(cum.loc[K.index].sum()),
                       hf=f"{hf}/18", ht=int(cum.loc[H.index].sum()),
                       cn=int(cum.loc[CN.index].sum()), no=int(cum.loc[N.index].sum())))
    frag["cumulative"] = md(cr, ["clause", "kf", "kt", "hf", "ht", "cn", "no"],
                            ["cumulative, n_screen's order", "TRUE families left",
                             "their templates", "hidden families", "their templates",
                             "candidate-new", "noise templates"])

    # screen combinations
    def comb(nm, cols):
        m = np.ones(len(V), bool)
        for c in cols:
            m &= V[c].to_numpy(bool)
        kf = int(K.assign(x=m[K.index]).groupby("bkey").x.max().sum())
        hf = int(H.assign(x=m[H.index]).groupby("hid_bkey").x.max().sum())
        return dict(screen=nm, kf=f"**{kf}/58 = {kf/58:.0%}**", kt=int(m[K.index].sum()),
                    hf=f"{hf}/18", cn=f"{int(m[CN.index].sum())}/18",
                    no=f"{int(m[N.index].sum())}/1201 = {m[N.index].mean():.1%}")
    combos = [
        ("the screen as coded, all five clauses",
         ["c_a", "c_b", "c_c1", "c_c2", "c_c3", "c_c4", "c_d", "c_e"]),
        ("as coded without (e)", ["c_a", "c_b", "c_c1", "c_c2", "c_c3", "c_c4", "c_d"]),
        ("as coded without (c3) and (e)", ["c_a", "c_b", "c_c1", "c_c2", "c_c4", "c_d"]),
        ("justified clause only: (c1)", ["c_c1"]),
        ("(c1)+(c4&prime;) only", ["c_c1", "c_c4_19"]),
        ("(c1)+(c4&prime;)+(d)", ["c_c1", "c_c4_19", "c_d"]),
        ("**recommended tier A** = (c1)+(c4&prime;)+(d)+(a&prime;)",
         ["c_c1", "c_c4_19", "c_d", "c_a25"]),
        ("**recommended tier B** = tier A + control (c2)",
         ["c_c1", "c_c4_19", "c_d", "c_a25", "c_c2"]),
        ("**recommended tier C** = tier B + co-timing (b&prime;)",
         ["c_c1", "c_c4_19", "c_d", "c_a25", "c_c2", "c_b_loud"]),
    ]
    frag["combos"] = md([comb(n, c) for n, c in combos],
                        ["screen", "kf", "kt", "hf", "cn", "no"],
                        ["screen", "TRUE families recovered", "their templates",
                         "hidden families", "candidate-new", "noise templates"])

    # sweeps
    Kl = K[K.c_c1]; Nl = N
    sw = []
    for w in (5, 8, 10, 12, 15, 18, 19, 20, 25, 30):
        kp = Kl.assign(ok=Kl.dz.abs() <= w).groupby("bkey").ok.max()
        keep = (N.c_c1.astype(bool) & (N.dz.abs() <= w))
        sw.append(dict(w=f"{w} km", kf=f"{int(kp.sum())}/50 located = {kp.sum()/50:.0%}",
                       kall=f"{int(kp.sum())}/58 = {kp.sum()/58:.0%}",
                       nr=f"{1-keep.mean():.0%}"))
    frag["dz_sweep"] = md(sw, ["w", "kf", "kall", "nr"],
                          ["&#124;&Delta;z&#124; width", "of the 50 LOCATED true families",
                           "of all 58", "noise population rejected"])
    sw = []
    for t in (2, 3, 4, 5, 6):
        kp = K.groupby("bkey").n_days.max() >= t
        sw.append(dict(t=f"&ge; {t} days", kf=f"{int(kp.sum())}/58 = {kp.mean():.0%}",
                       nr=f"{1-(N.n_days>=t).mean():.0%}"))
    frag["days_sweep"] = md(sw, ["t", "kf", "nr"],
                            ["clause (d)", "TRUE families passing", "noise rejected"])
    sw = []
    for x in (0.1, 0.25, 0.5, 1.0, 1e9):
        ok = ((K.fa_quiet_per_tpl_day <= x)
              & (K.rate_tremor_per_tpl_day >= 3 * K.fa_quiet_per_tpl_day) & (K.n_det > 0))
        okn = ((N.fa_quiet_per_tpl_day <= x)
               & (N.rate_tremor_per_tpl_day >= 3 * N.fa_quiet_per_tpl_day) & (N.n_det > 0))
        kf = K.assign(o=ok).groupby("bkey").o.max()
        sw.append(dict(x=("no cap, ratio only" if x > 1e8 else f"&le; {x}/tpl-day"),
                       kf=f"{int(kf.sum())}/58 = {kf.mean():.0%}",
                       nr=f"{1-okn.mean():.0%}"))
    frag["fa_sweep"] = md(sw, ["x", "kf", "nr"],
                          ["clause (a) cap (ratio 3&times; kept throughout)",
                           "TRUE families passing", "noise rejected"])

    # per-family table for the 58
    g = K.groupby("bkey")
    PF = pd.DataFrame(dict(n_tpl=g.size(), n_det=g.n_det.max(), n_days=g.n_days.max(),
                           frac_best=g.n_det.max() * 0))
    for c, _ in CL:
        PF[c] = g[c].max()
    PF["dz_best"] = g.apply(lambda d: (d.dz.abs().min() if d.c_c1.any() else np.nan),
                            include_groups=False)
    PF["located"] = g.c_c1.max()
    PF["n_pass"] = PF[[c for c, _ in CL]].sum(1)
    PF = PF.drop(columns=["frac_best"]).reset_index()
    PF.to_csv(f"{OUT}/per_family_58.csv", index=False)
    lines = ["| bkey | templates | max n_det | max days | " +
             " | ".join(c for c, _ in CL) + " | best &#124;&Delta;z&#124; | clauses passed |",
             "|" + "|".join(["---"] * (5 + len(CL) + 1)) + "|"]
    y = {True: "Y", False: "&middot;"}
    for r in PF.itertuples():
        lines.append("| " + str(int(r.bkey)) + " | " + str(int(r.n_tpl)) + " | "
                     + str(int(r.n_det)) + " | " + str(int(r.n_days)) + " | "
                     + " | ".join(y[bool(getattr(r, c))] for c, _ in CL) + " | "
                     + (f"{r.dz_best:.1f}" if np.isfinite(r.dz_best) else "&mdash;") + " | "
                     + f"{int(r.n_pass)}/8 |")
    frag["per_family"] = "\n".join(lines)
    json.dump(frag, open(f"{OUT}/fragments.json", "w"), indent=1)
    V.to_csv(f"{OUT}/screen_north.csv", index=False)
    for k in frag:
        print(f"\n===== {k} =====\n{frag[k][:4000]}")


if __name__ == "__main__":
    main()
