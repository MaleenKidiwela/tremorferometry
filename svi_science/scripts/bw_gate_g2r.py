#!/usr/bin/env python
"""bw_gate_g2r -- G2, third reference variant: Bostock's events RE-TIMED onto the frozen axis by the other CN stations
(bw_gate_retime.retime with CN4 minus the station under test; truncated 2004/2005 seconds handled there).  Then the
Track A-2 stack comparison (bw_gate_g2_epoch.stack_leg) and the pre-registered G2 rule against the no-change PGC + NLLB.
Control leg = 2012 test-block events re-timed with the same stations, matched N.  Also the rebuild check on re-timed
events (build-block stack vs test-block stack, against the 2012 split).

Out: OUT/tables/g2r_family_rows.csv, g2r_pairs.csv, g2r_station_summary.csv, g2r_summary.csv, g2r_rebuild.csv; figures/g2r_transfer.png
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bw_gate_g2_epoch as E
import bw_gate_retime as RT
G = E.G
TAB = E.TAB
WINS = ["2004-07", "2005-09", "2008-05"]


def ev_rt(win, x, test=False, build=False):
    yr = int(E.WIN[win][0][:4])
    e = E.events(win, test_blocks=test, build_blocks=build)
    return RT.retime(e, yr, E.CACHES[yr], [s for s in G.CN4 if s != x])


def main():
    rows, rb = [], []
    for x in G.CN4:
        ctl_all = ev_rt("2012-09", x)
        role12 = E.setup()["R"].block_role(ctl_all.t_cat.to_numpy())
        for w in WINS:
            en = ev_rt(w, x)
            N = pd.DataFrame(E.stack_leg(w, x, "new", ev=en))
            caps = {int(r.family): int(r.n_events) for r in N[N.status == "ok"].itertuples()}
            Cc = pd.DataFrame(E.stack_leg("2012-09", x, f"ctl_for_{w}", caps=caps, ev=ctl_all[role12 == "fit"]))
            role = "target" if x in ("VGZ", "SNB") and w in WINS else "control"
            N["role"] = Cc["role"] = role; N["epoch_label"] = f"{x} {w}"; Cc["epoch_label"] = "2012 test blocks"
            N["retime_lag_med"] = float(en.lag.median()); N["retime_netcc_med"] = float(en.net_cc.median())
            rows += [N, Cc]
            G.log(f"[g2r] {w} {x}: {len(en)} events, retime lag median {en.lag.median():+.2f} s, net cc {en.net_cc.median():.3f}; ok families {len(caps)}")
            # rebuild on re-timed events: build-block stack vs test-block stack (old window) and the same on 2012
            if w in ("2004-07", "2005-09"):
                yr = int(E.WIN[w][0][:4])
                role_w = E.setup()["R"].block_role(en.t_cat.to_numpy())
                res = {}
                for tag, ewin, eb, et in ((w, w, en[role_w == "insample"], en[role_w == "fit"]),
                                          ("2012", "2012-09", ctl_all[role12 == "insample"], ctl_all[role12 == "fit"])):
                    _, Sb = E.stack_leg(ewin, x, "b", ev=eb, keep_stacks=True)
                    rt_, St = E.stack_leg(ewin, x, "t", ev=et, keep_stacks=True)
                    regt = {r["family"]: r.get("registration", np.nan) for r in rt_ if r.get("status") == "ok"}
                    for f in set(Sb) & set(St):
                        c, l = E.setup()["A4"].best_cc(Sb[f].astype(float), St[f].astype(float), 380)
                        res.setdefault(f, {})[tag] = (c, regt.get(f, np.nan))
                for f, d in res.items():
                    if len(d) == 2:
                        rb.append(dict(sta=x, window=w, family=f, cc_bt_old=d[w][0], reg_t_old=d[w][1], cc_bt_2012=d["2012"][0], reg_t_2012=d["2012"][1]))
    A = pd.concat(rows, ignore_index=True); A.to_csv(f"{TAB}/g2r_family_rows.csv", index=False)
    M = E.pair_rows(A); M.to_csv(f"{TAB}/g2r_pairs.csv", index=False)
    SS = E.station_summary(M); SS.to_csv(f"{TAB}/g2r_station_summary.csv", index=False)
    out = []
    tests = [("VGZ", ["2008-05"], "VGZ BH CMG3T/GD2_v21 2007-02..2009-07"), ("SNB", ["2008-05"], "SNB BH CMG3ESP/GD2_v12 2006-04..2008-12"),
             ("VGZ", ["2004-07"], "VGZ BH CMG3ESP/GD1 2003-02..2007-02"), ("VGZ", ["2005-09"], "VGZ BH CMG3ESP/GD1 2003-02..2007-02"),
             ("VGZ", ["2004-07", "2005-09"], "VGZ BH CMG3ESP/GD1 2003-02..2007-02"),
             ("SNB", ["2004-07"], "SNB BH CMG40T/GD1 ..2006-04"), ("SNB", ["2005-09"], "SNB BH CMG40T/GD1 ..2006-04"),
             ("SNB", ["2004-07", "2005-09"], "SNB BH CMG40T/GD1 ..2006-04")]
    for sta, wins, ep in tests:
        gt = M[(M.sta == sta) & M.window.isin(wins)]; gc = M[M.window.isin(wins) & M.sta.isin(["PGC", "NLLB"])]
        if not len(gt) or not len(gc): continue
        rt, rc = E.ratios(gt), E.ratios(gc); ci = E.boot_diff(gt, gc)
        r = dict(station=sta, epoch=ep, windows="+".join(wins), n_fam=len(gt), n_ctl_rows=len(gc), cc_med=gt.cc.median(),
                 cc_ctl_med=gt.cc_ctl.median(), nochange_cc_med=gc.cc.median(), nochange_cc_ctl_med=gc.cc_ctl.median(),
                 cc_ratio=rt[0], cc_ratio_nochange=rc[0], reg_ratio=rt[1], reg_ratio_nochange=rc[1],
                 d_cc_ratio_90ci=ci[:2], d_reg_ratio_90ci=ci[2:], dS_med=gt.dS.median(), nochange_dS_med=gc.dS.median(),
                 dS_mad=float((gt.dS - gt.dS.median()).abs().median()), flip_frac=float((gt.cc <= -0.6).mean()),
                 pass_cc=bool(rt[0] >= rc[0] - 0.05), pass_reg=bool(rt[1] >= rc[1] - 0.05),
                 pass_dS=bool(abs(gt.dS.median()) <= 0.10), pass_flips=bool((gt.cc <= -0.6).mean() <= 0.10))
        r["transfer_pass"] = bool(r["pass_cc"] and r["pass_reg"] and r["pass_dS"] and r["pass_flips"])
        out.append(r)
        G.log(f"[g2r] {sta} {'+'.join(wins)}: cc {gt.cc.median():.3f} (2012 leg {gt.cc_ctl.median():.3f}) ratio {rt[0]:.3f} vs no-change "
              f"{rc[0]:.3f}; reg {rt[1]:.3f} vs {rc[1]:.3f}; dS {gt.dS.median():+.2f}; flips {r['flip_frac']:.2f} -> {'PASS' if r['transfer_pass'] else 'FAIL'}")
    pd.DataFrame(out).to_csv(f"{TAB}/g2r_summary.csv", index=False)
    Rb = pd.DataFrame(rb); Rb.to_csv(f"{TAB}/g2r_rebuild.csv", index=False)
    if len(Rb):
        s = Rb.groupby(["sta", "window"])[["cc_bt_old", "cc_bt_2012", "reg_t_old", "reg_t_2012"]].median().round(3)
        s["validates"] = (s.cc_bt_old >= s.cc_bt_2012 - 0.05) & (s.reg_t_old >= s.reg_t_2012 - 0.05)
        s.to_csv(f"{TAB}/g2r_rebuild_summary.csv"); G.log(f"[g2r rebuild]\n{s}")
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
    for k, (col, ctl, lab) in enumerate((("cc_med", "cc_ctl_med", "median cc with the 2012 stack"), ("reg_ratio", None, "registration ratio (old / 2012 leg)"))):
        ax = axs[k]
        for i, w in enumerate(WINS):
            for _, r in SS[SS.window == w].iterrows():
                dx = {"PGC": -.21, "NLLB": -.07, "VGZ": .07, "SNB": .21}[r.sta]
                ax.scatter(i + dx, r[col], marker={"PGC": "o", "NLLB": "^", "VGZ": "s", "SNB": "D"}[r.sta], s=60,
                           c="0.55" if r.sta in ("PGC", "NLLB") else ("#2471a3" if r.sta == "VGZ" else "#c0392b"),
                           label=f"{r.sta} ({'no change' if r.sta in ('PGC', 'NLLB') else 'old instrument'})" if i == 0 else None)
                if ctl:
                    ax.scatter(i + dx, r[ctl], marker="_", s=120, c="k", label="same station, 2012 leg" if (i == 0 and r.sta == "PGC") else None)
        ax.set_xticks(range(len(WINS))); ax.set_xticklabels(WINS); ax.set_title(lab); ax.legend(fontsize=7)
    fig.suptitle("G2: 2012 templates on the old VGZ / SNB instruments; Bostock events re-timed by the other CN stations")
    plt.tight_layout(); plt.savefig(f"{G.OUT}/figures/g2r_transfer.png", dpi=130); plt.close()


if __name__ == "__main__":
    main()
