#!/usr/bin/env python
"""bw_gate_g2w_stack -- G2 with the A-2 ORIGINAL reference events: withheld-station detections of the forward detector
(bw_gate_g2w_run.sh), net_mad >= 8.75 (T_core) and the cross-family local maximum (keep_xfam8), strongest first,
<= 400 per family.  Control leg = the same withheld network on the 2012 truth fit days, TEST blocks only, matched N.
Stacking / comparison = bw_gate_g2_epoch.stack_leg (Track A-2 step-4 functions).  Same pass rule as GATES.md G2.

Deviation from the pre-registered G2 (recorded in GATES.md): added after the Bostock-referenced run showed the
never-changed PGC and NLLB failing in 2004/2005 exactly like VGZ/SNB (catalogue time base, bw_gate_g2_lagdiag.py).

Out: OUT/tables/g2w_family_rows.csv, g2w_pairs.csv, g2w_station_summary.csv, g2w_summary.csv; figures/g2w_transfer.png
"""
from __future__ import annotations
import os, sys, glob, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bw_gate_g2_epoch as E
G = E.G
TAB = E.TAB
WINS = {"2004-07": ("2004", "CN4"), "2005-09": ("2005", "B1"), "2008-05": ("2008", "B3")}
TARGET = {"2004-07": ["VGZ", "SNB"], "2005-09": ["VGZ", "SNB"], "2008-05": ["VGZ", "SNB"]}


def dets(run, test_blocks=False):
    fs = [p for p in glob.glob(f"{G.OUT}/mf/feat/{run}/*.parquet") if not p.endswith("_fd.parquet") and os.path.getsize(p) > 0]
    parts = []
    for p in fs:
        try:
            A = pd.read_parquet(p, columns=["series", "fam", "t", "net_mad", "keep_xfam8", "n_sites"])
        except Exception:
            continue
        parts.append(A[(A.series == "data") & (A.net_mad >= 8.75) & A.keep_xfam8.astype(bool)])
    D = pd.concat(parts, ignore_index=True)
    if test_blocks:
        D = D[E.setup()["R"].block_role(D.t.to_numpy()) == "fit"]
    return pd.DataFrame(dict(t_origin=D.t.to_numpy(), fam=D.fam.astype(int).to_numpy(), mag=D.net_mad.to_numpy()))


def main():
    rows = []
    for w, (yr, base) in WINS.items():
        for x in ["PGC", "VGZ", "NLLB", "SNB"]:
            rn, rc = f"g2w{yr}_{base}no{x}", f"g2w2012_{base}no{x}"
            if not os.path.isdir(f"{G.OUT}/mf/feat/{rn}") or not os.path.isdir(f"{G.OUT}/mf/feat/{rc}"):
                G.log(f"[g2w] missing {rn} / {rc}"); continue
            Dn = dets(rn); Dc = dets(rc, test_blocks=True)
            N = pd.DataFrame(E.stack_leg(w, x, "new", ev=Dn))
            caps = {int(r.family): int(r.n_events) for r in N[N.status == "ok"].itertuples()}
            Cc = pd.DataFrame(E.stack_leg("2012-09", x, f"ctl_for_{w}", caps=caps, ev=Dc))
            role = "target" if x in TARGET[w] else "control"
            for T_ in (N, Cc):
                T_["role"] = role; T_["epoch_label"] = f"{x} {w}" if T_ is N else "2012 (control leg, test blocks)"
                T_["n_det_run"] = len(Dn) if T_ is N else len(Dc)
            rows += [N, Cc]
            G.log(f"[g2w] {w} {x}: {len(Dn)} detections new, {len(Dc)} 2012 test-block; ok families {len(caps)}")
    A = pd.concat(rows, ignore_index=True); A.to_csv(f"{TAB}/g2w_family_rows.csv", index=False)
    M = E.pair_rows(A); M.to_csv(f"{TAB}/g2w_pairs.csv", index=False)
    SS = E.station_summary(M); SS.to_csv(f"{TAB}/g2w_station_summary.csv", index=False)
    out = []
    tests = [("VGZ", ["2008-05"], "VGZ BH CMG3T/GD2_v21 2007-02..2009-07"), ("SNB", ["2008-05"], "SNB BH CMG3ESP/GD2_v12 2006-04..2008-12"),
             ("VGZ", ["2004-07"], "VGZ BH CMG3ESP/GD1 2003-02..2007-02"), ("VGZ", ["2005-09"], "VGZ BH CMG3ESP/GD1 2003-02..2007-02"),
             ("VGZ", ["2004-07", "2005-09"], "VGZ BH CMG3ESP/GD1 2003-02..2007-02"),
             ("SNB", ["2004-07"], "SNB BH CMG40T/GD1 ..2006-04"), ("SNB", ["2005-09"], "SNB BH CMG40T/GD1 ..2006-04"),
             ("SNB", ["2004-07", "2005-09"], "SNB BH CMG40T/GD1 ..2006-04")]
    for sta, wins, ep in tests:
        gt = M[(M.sta == sta) & M.window.isin(wins)]
        gc = M[M.window.isin(wins) & M.sta.isin(["PGC", "NLLB"])]
        if not len(gt) or not len(gc): continue
        rt, rc = E.ratios(gt), E.ratios(gc); ci = E.boot_diff(gt, gc)
        r = dict(station=sta, epoch=ep, windows="+".join(wins), n_fam=len(gt), n_ctl_rows=len(gc), ctl_stations="PGC,NLLB",
                 cc_ratio=rt[0], cc_ratio_nochange=rc[0], reg_ratio=rt[1], reg_ratio_nochange=rc[1],
                 d_cc_ratio_90ci=ci[:2], d_reg_ratio_90ci=ci[2:], cc_med=gt.cc.median(), cc_ctl_med=gt.cc_ctl.median(),
                 nochange_cc_med=gc.cc.median(), dS_med=gt.dS.median(), dS_mad=float((gt.dS - gt.dS.median()).abs().median()),
                 nochange_dS_med=gc.dS.median(), flip_frac=float((gt.cc <= -0.6).mean()),
                 pass_cc=bool(rt[0] >= rc[0] - 0.05), pass_reg=bool(rt[1] >= rc[1] - 0.05),
                 pass_dS=bool(abs(gt.dS.median()) <= 0.10), pass_flips=bool((gt.cc <= -0.6).mean() <= 0.10))
        r["transfer_pass"] = bool(r["pass_cc"] and r["pass_reg"] and r["pass_dS"] and r["pass_flips"])
        out.append(r)
        G.log(f"[g2w] {sta} {'+'.join(wins)}: cc ratio {rt[0]:.3f} vs no-change {rc[0]:.3f}; reg {rt[1]:.3f} vs {rc[1]:.3f}; "
              f"dS {gt.dS.median():+.2f} (no-change {gc.dS.median():+.2f}); flips {r['flip_frac']:.2f} -> {'PASS' if r['transfer_pass'] else 'FAIL'}")
    pd.DataFrame(out).to_csv(f"{TAB}/g2w_summary.csv", index=False)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    wl = list(WINS)
    for k, (col, lab) in enumerate((("cc_med", "median cc with the 2012 stack"), ("cc_ratio", "cc ratio (old / 2012 control leg)"),
                                    ("reg_ratio", "registration ratio"))):
        ax = axs[k]
        for i, w in enumerate(wl):
            g = SS[SS.window == w]
            for _, r in g.iterrows():
                mk = {"PGC": "o", "NLLB": "^", "VGZ": "s", "SNB": "D"}[r.sta]
                ax.scatter(i + {"PGC": -.21, "NLLB": -.07, "VGZ": .07, "SNB": .21}[r.sta], r[col], marker=mk, s=60,
                           c="0.55" if r.sta in ("PGC", "NLLB") else ("#2471a3" if r.sta == "VGZ" else "#c0392b"),
                           label=f"{r.sta} ({'no change' if r.sta in ('PGC', 'NLLB') else 'old instrument'})" if i == 0 else None)
            if col == "cc_med":
                ax.scatter(np.full(len(g), i + .3), g.cc_ctl_med, marker="_", s=80, c="k", label="2012 control leg" if i == 0 else None)
        ax.set_xticks(range(len(wl))); ax.set_xticklabels(wl); ax.set_title(lab); ax.axhline(1 if col != "cc_med" else 0, c="k", lw=.5)
        if k == 0: ax.legend(fontsize=7)
    fig.suptitle("G2: 2012 templates on the old VGZ / SNB instruments; reference = withheld-station detections (Track A-2 protocol)")
    plt.tight_layout(); plt.savefig(f"{G.OUT}/figures/g2w_transfer.png", dpi=130); plt.close()


if __name__ == "__main__":
    main()
