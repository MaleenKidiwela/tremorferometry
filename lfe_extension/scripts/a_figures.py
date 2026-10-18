#!/usr/bin/env python
"""Track A figures: (1) matched-FA ROC curves, (2) FA vs threshold with the frozen T_ext/T_core, (3) HF-ratio
distributions with the adopted cut, (4) 2010-vs-2012 recall by magnitude bin (step 2)."""
from __future__ import annotations
import json, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = "/home/jovyan/tremorferometry/lfe_extension/results"
FIG = "/home/jovyan/tremorferometry/lfe_extension/figures"
os.makedirs(FIG, exist_ok=True)

# ---------- 1. ROC
G = pd.read_csv(f"{R}/e6_fa_recall_grid.csv")
fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
for tag, g in G.groupby("tag"):
    g = g.sort_values("quiet_fa")
    lw = 2.4 if tag == "e6_mergeB" else 1.2
    ax[0].plot(g.quiet_fa, g.recall_own, "o-", lw=lw, ms=4, label=tag.replace("e6_", ""))
try:
    G12 = pd.read_csv(f"{R}/a1_fa_recall_grid_L12.csv").sort_values("quiet_fa")
    ax[0].plot(G12.quiet_fa, G12.recall_own, "s-", color="k", lw=2.4, ms=5, label="mergeB 12 s (ADOPTED)")
except FileNotFoundError:
    pass
ax[0].set_xscale("log"); ax[0].set_xlabel("quiet-day false alarms / template / day")
ax[0].set_ylabel("own-family recall, Bostock 2012 test half")
ax[0].axvline(1.0, color="0.6", ls=":"); ax[0].axvline(0.1, color="0.6", ls=":")
ax[0].set_title("Matched-FA ROC (2012, 15 stations)"); ax[0].legend(fontsize=7.5); ax[0].grid(alpha=.3)
for tag, g in G.groupby("tag"):
    g = g.sort_values("quiet_fa")
    ax[1].plot(g.quiet_fa, g.ets_precision_own, "o-", lw=2.4 if tag == "e6_mergeB" else 1.2, ms=4, label=tag.replace("e6_", ""))
ax[1].set_xscale("log"); ax[1].set_xlabel("quiet-day false alarms / template / day")
ax[1].set_ylabel("own-family precision on ETS days"); ax[1].grid(alpha=.3)
ax[1].set_title("Precision at matched FA")
fig.tight_layout(); fig.savefig(f"{FIG}/a0_matched_fa_roc.png", dpi=140); plt.close(fig)

# ---------- 2. FA vs T
fig, ax = plt.subplots(figsize=(7.2, 4.8))
for lbl, fn, st in (("mergeB 12 s, no veto", f"{R}/a1_operating_fa_vs_T_a1_mergeB_L12.csv", "--"),
                    ("mergeB 12 s + HF veto", f"{R}/a1_frozen_fa_vs_T_a1_mergeB_L12.csv", "-")):
    if not os.path.exists(fn):
        continue
    g = pd.read_csv(fn)
    ax.semilogy(g["T"], g.quiet_fa, st, color="C0", label=f"quiet day, {lbl}")
    ax.semilogy(g["T"], g.post_fa, st, color="C3", label=f"post-SSE, {lbl}")
fr = json.load(open(f"{R}/frozen_operating_point.json"))
ax.axhline(1.0, color="C0", lw=.8, alpha=.5); ax.axhline(1.5, color="C3", lw=.8, alpha=.5); ax.axhline(0.1, color="C0", lw=.8, alpha=.5)
ax.axvline(fr["T_ext"], color="k", ls="-", lw=1.2); ax.axvline(fr["T_core"], color="k", ls=":", lw=1.2)
ax.text(fr["T_ext"], 40, f" T_ext {fr['T_ext']}", fontsize=8); ax.text(fr["T_core"], 12, f" T_core {fr['T_core']}", fontsize=8)
ax.set_xlabel("network threshold T (MAD)"); ax.set_ylabel("detections / template / day")
ax.set_title("False alarms vs threshold; thresholds set by FA only"); ax.legend(fontsize=7.5); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(f"{FIG}/a1_fa_vs_T.png", dpi=140); plt.close(fig)

# ---------- 3. HF ratio
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
for k, (tag, L) in enumerate((("e6_mergeB", 8.0), ("a1_mergeB_L12", 12.0))):
    fn = f"{R}/a1_hf_ratio_{tag}.parquet"
    if not os.path.exists(fn):
        continue
    d = pd.read_parquet(fn); d = d[d.L == L]
    dec = pd.read_csv(f"{R}/a1_hf_ratio_{tag}_decision.csv"); cut = float(dec[dec.sta == "PGC"].cut_p99.iloc[0])
    v = d.hfr_PGC
    b = np.logspace(-3, 3, 70)
    ax[k].hist(v[d.match_own & (d.net_mad >= 8)].clip(1e-3, 1e3), bins=b, alpha=.6, density=True, label="matched to own-family Bostock event")
    ax[k].hist(v[d.doy.isin(range(228, 243)) & (d.net_mad >= 8)].clip(1e-3, 1e3), bins=b, alpha=.6, density=True, label="quiet-day detections")
    ax[k].axvline(cut, color="k", lw=1.5); ax[k].text(cut, ax[k].get_ylim()[1] * .85, f" cut {cut:.2f}", fontsize=8)
    ax[k].set_xscale("log"); ax[k].set_xlabel("E(10-19 Hz) / E(2-8 Hz) at PGC over the scored window")
    ax[k].set_title(f"{tag} (L = {L:g} s)"); ax[k].legend(fontsize=7.5); ax[k].grid(alpha=.3)
ax[0].set_ylabel("density")
fig.tight_layout(); fig.savefig(f"{FIG}/a1_hf_ratio.png", dpi=140); plt.close(fig)

# ---------- 4. 2010 vs 2012 recall by magnitude
fn = f"{R}/a2_transfer_table.csv"
if os.path.exists(fn):
    Tb = pd.read_csv(fn)
    bins = ["M<1.6", "M1.6-1.8", "M1.8-2.0", "M>=2.0"]
    fig, ax = plt.subplots(1, len(Tb), figsize=(5.6 * len(Tb), 4.4), squeeze=False)
    for i, r in Tb.iterrows():
        x = np.arange(len(bins)); w = .38
        ax[0][i].bar(x - w / 2, [r[f"rec2012_{b}"] for b in bins], w, label="2012 ablation (7 persistent sites)")
        ax[0][i].bar(x + w / 2, [r[f"rec2010_{b}"] for b in bins], w, label="2010 frozen transfer")
        ax[0][i].set_xticks(x); ax[0][i].set_xticklabels(bins, fontsize=8)
        ax[0][i].set_title(f"{r.threshold} = {r['T']:g} MAD   ratio {r.recall_ratio:.2f}")
        ax[0][i].set_ylabel("own-family recall"); ax[0][i].legend(fontsize=8); ax[0][i].grid(alpha=.3, axis="y")
    fig.tight_layout(); fig.savefig(f"{FIG}/a2_transfer_recall_by_mag.png", dpi=140); plt.close(fig)
# ---------- 5. step 3b: withheld-station stack test
fn = f"{R}/a3b_quietday_stack_test.csv"
if os.path.exists(fn):
    Rt = pd.read_csv(fn)
    q = Rt[(Rt.set == "quiet") & (Rt["T"] == 8.25) & Rt.n_stacked.notna()]
    c = Rt[(Rt.set == "ets_control") & (Rt["T"] == 8.25) & Rt.n_stacked.notna()]
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.4))
    b = np.linspace(-6, 12, 40)
    ax[0].hist(q.excess_db, bins=b, alpha=.65, label=f"quiet-day detections (n fam {len(q)})")
    ax[0].hist(c.excess_db, bins=b, alpha=.65, label=f"ETS-day own-matched (positive control, n fam {len(c)})")
    ax[0].axvline(3, color="k", ls="--", lw=1); ax[0].text(3.1, ax[0].get_ylim()[1]*.9, "PASS bar", fontsize=8)
    ax[0].set_xlabel("withheld-station S-SNR excess over the matched random-time null (dB)")
    ax[0].set_ylabel("families"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
    ax[0].set_title("T_ext = 8.25: quiet-day detections do not register")
    ax[1].scatter(q.cc_null_mean, q.cc_template, s=18, label="quiet-day")
    ax[1].scatter(c.cc_null_mean, c.cc_template, s=18, label="ETS control")
    lim = [0, max(0.9, float(np.nanmax(c.cc_template)) * 1.1)]
    ax[1].plot(lim, lim, "k--", lw=.8); ax[1].set_xlim(lim); ax[1].set_ylim(lim)
    ax[1].set_xlabel("cc of the random-time stack to the template"); ax[1].set_ylabel("cc of the detection stack to the template")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.3); ax[1].set_title("Waveform identity at the withheld station")
    zp = f"{R}/a3b_quietday_stack_test_stacks.npz"
    if os.path.exists(zp):
        Z = np.load(zp)
        best = c.sort_values("excess_db", ascending=False).iloc[0]
        worst = q.sort_values("excess_db").iloc[len(q) // 2]
        t = np.arange(int(2 * 20.0 * 20.0)) / 20.0 - 20.0
        for lbl, row, kind, col in (("ETS control (registers)", best, "ets_control", "C1"),
                                    ("quiet-day (does not)", worst, "quiet", "C0")):
            k = f"{row.sta}|{int(row.fam)}|8.25|{kind}"
            if k in Z:
                w = Z[k]; e = np.sqrt((w ** 2).sum(0)); e = e / e.max()
                ax[2].plot(t, e + (0 if kind == "quiet" else 1.2), col,
                           label=f"{lbl}: fam {int(row.fam)} @ {row.sta}, N {int(row.n_stacked)}, {row.excess_db:+.1f} dB")
        ax[2].axvline(0, color="k", lw=.8); ax[2].set_xlim(-10, 10)
        ax[2].set_xlabel("time from predicted S (s)"); ax[2].set_ylabel("normalised 3-C envelope (offset)")
        ax[2].legend(fontsize=7.5); ax[2].grid(alpha=.3); ax[2].set_title("Withheld-station stacks")
    fig.tight_layout(); fig.savefig(f"{FIG}/a3b_quietday_stack_test.png", dpi=140); plt.close(fig)

# ---------- 6. step 3: what template iteration does
fn = f"{R}/a3_reject_tests.csv"
if os.path.exists(fn):
    D = pd.read_csv(fn)
    Q = pd.read_csv(f"{R}/e5r_template_qc_a3_it1.csv")[["bkey", "sta", "source"]]
    M = D.merge(Q, on=["bkey", "sta"], how="left")
    it = M[(M.source == "iterated") & (~M.new)]
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.4))
    g = it.groupby("sta").dsnr.median().sort_values()
    col = ["C3" if v < 0 else "C0" for v in g]
    ax[0].barh(g.index, g.values, color=col)
    ax[0].axvline(0, color="k", lw=.8); ax[0].axvline(-1, color="0.5", ls="--", lw=.8)
    ax[0].set_xlabel("median change in stack S-SNR, iterated - frozen (dB)")
    ax[0].set_title("Iteration gains at far stations, loses at near ones"); ax[0].grid(alpha=.3, axis="x")
    ax[1].scatter(it.ccx_old, it.ccx_new, s=10, alpha=.5)
    lim = [0, 1]; ax[1].plot(lim, lim, "k--", lw=.8); ax[1].set_xlim(lim); ax[1].set_ylim(lim)
    ax[1].axhline(0.7, color="C3", lw=.8)
    ax[1].set_xlabel("max cc to another family, frozen"); ax[1].set_ylabel("iterated")
    ax[1].set_title("No cross-family collapse"); ax[1].grid(alpha=.3)
    nw = M[M.new].sta.value_counts()
    ax[2].bar(nw.index, nw.values, color="C2")
    ax[2].set_ylabel("new (family, station) templates"); ax[2].set_title(f"Iteration adds {int(M.new.sum())} templates the seed set lacked")
    ax[2].tick_params(axis="x", rotation=60); ax[2].grid(alpha=.3, axis="y")
    fig.tight_layout(); fig.savefig(f"{FIG}/a3_iteration.png", dpi=140); plt.close(fig)

# ---------- 7. step 3 held-out recall at matched FA
fn = f"{R}/a3_heldout_matched_fa_grid.csv"
if os.path.exists(fn):
    G = pd.read_csv(fn)
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    for lbl, g in G.groupby("set"):
        g = g.sort_values("quiet_fa")
        ax.plot(g.quiet_fa, g.recall, "o-", label=f"{lbl} templates")
    ax.set_xscale("log"); ax.set_xlabel("2010 quiet-day false alarms / template / day")
    ax.set_ylabel("own-family recall, Bostock 2010")
    ax.set_title("Held-out year (2010) at matched false-alarm rate"); ax.legend(); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(f"{FIG}/a3_heldout_matched_fa.png", dpi=140); plt.close(fig)

print("[figures] written to", FIG, sorted(os.listdir(FIG)))
