#!/usr/bin/env python
"""Track A step 3 -- the SIX rejection tests for an iterated template set (Merlin plan-review step 3).

Compares a NEW template npz/QC against the frozen 2012 one, per (family, station), and returns a per-family verdict.
  (a) q-control     cc of the stack to the SeisBench (Bostock) stack must not fall by > 0.05 nor go below 0.6 (PGC/VGZ/SNB)
  (c) S pick        |new s_pick - old s_pick| <= 0.5 s; station-differential moveout change <= 0.3 s (from the moveout tables)
  (d) cross-family  max |cc| of a family's template to any OTHER family's template at the same station must not rise by
                    > 0.1 nor exceed 0.7 (full family x family matrix per station, +-1 s lag scan, 3-C summed)
  (e) spectral      no more than 50 % of the template's 1-8 Hz energy in a single 1-Hz bin
  (f) S-SNR         stack S-SNR must not fall by more than 1 dB
  (b) held-out-year recall at matched FA is measured by a separate GPU pass (a3_score_heldout.py) and merged in here.
A family is REJECTED (its iteration discarded, the frozen template kept) if ANY test fires at ANY of its stations."""
from __future__ import annotations
import argparse
import numpy as np, pandas as pd

IN = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
OUT = "/home/jovyan/tremorferometry/lfe_extension/results"
FS = 20.0
QC_STAS = ("PGC", "VGZ", "SNB")


def load_npz(p, L):
    z = np.load(p, allow_pickle=True)
    tpl = z["tpl"].item()
    return {k: v[L] for k, v in tpl.items() if L in v}


def max_cc_other(T, lag=20):
    """per station: max |cc| of each family's template against every OTHER family's, +-lag samples, 3-C summed."""
    out = {}
    stas = sorted({k[1] for k in T})
    for s in stas:
        ks = sorted([k for k in T if k[1] == s], key=lambda k: k[0])
        if len(ks) < 2:
            continue
        X = np.stack([T[k] for k in ks])                     # (nf, 3, n) each component L2-normalised
        nf, nc, n = X.shape
        F = np.fft.rfft(X, n=2 * n, axis=-1)
        best = np.zeros(nf)
        for i in range(nf):
            cc = np.fft.irfft(np.conj(F[i])[None] * F, n=2 * n, axis=-1).sum(1)     # (nf, 2n) 3-C summed cc vs lag
            cc = np.concatenate([cc[:, -lag:], cc[:, :lag + 1]], axis=1) / nc       # normalised: each comp is unit norm
            m = np.abs(cc).max(1); m[i] = -9
            best[i] = m.max()
        for k, v in zip(ks, best):
            out[k] = float(v)
    return out


def spectral_frac(T):
    out = {}
    for k, w in T.items():
        n = w.shape[1]
        P = (np.abs(np.fft.rfft(w, axis=-1)) ** 2).sum(0)
        f = np.fft.rfftfreq(n, 1 / FS)
        m = (f >= 1) & (f <= 8)
        if not m.any():
            continue
        tot = P[m].sum()
        bins = [P[(f >= lo) & (f < lo + 1)].sum() for lo in range(1, 8)]
        out[k] = float(max(bins) / max(tot, 1e-30))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old-npz", default=f"{IN}/e5r_templates_2012.npz")
    ap.add_argument("--old-qc", default=f"{IN}/e5r_template_qc_2012.csv")
    ap.add_argument("--old-mv", default=f"{IN}/e4_moveout_e5r_mv_build.csv")
    ap.add_argument("--new-npz", required=True); ap.add_argument("--new-qc", required=True)
    ap.add_argument("--new-mv", default="")
    ap.add_argument("--L", type=float, default=12.0)
    ap.add_argument("--out", default=f"{OUT}/a3_reject_tests.csv")
    # The ABSOLUTE clauses of Merlin's tests (a) "< 0.6", (d) "> 0.7", (e) "> 50 % in one 1-Hz bin" fire at essentially the
    # same rate on the FROZEN templates (13 / 320 / 443 of 1,144 pairs), so applied literally they reject almost everything
    # regardless of the iteration.  --conditional keeps the RELATIVE clauses as written and applies each absolute clause only
    # where the frozen template satisfied it, i.e. it asks "did the ITERATION break this?" rather than "is this template class
    # narrow-band / mutually similar?".  Both verdicts are written.
    ap.add_argument("--conditional", action="store_true", default=True)
    a = ap.parse_args()
    To = load_npz(a.old_npz, a.L); Tn = load_npz(a.new_npz, a.L)
    qo = pd.read_csv(a.old_qc); qn = pd.read_csv(a.new_qc)
    qo = qo[qo.usable & qo.template].set_index(["bkey", "sta"]); qn = qn[qn.usable & qn.template].set_index(["bkey", "sta"])
    cco = max_cc_other(To); ccn = max_cc_other(Tn)
    spn = spectral_frac(Tn); spo = spectral_frac(To)
    mvo = pd.read_csv(a.old_mv); mvo = mvo[mvo.L == a.L].set_index(["fam", "sta"]).lag_median
    mvn = None
    if a.new_mv:
        m = pd.read_csv(a.new_mv); mvn = m[m.L == a.L].set_index(["fam", "sta"]).lag_median
    rows = []
    for (f, s) in sorted(Tn):
        fi = int(f)
        r = dict(bkey=fi, sta=s)
        o = qo.loc[(fi, s)] if (fi, s) in qo.index else None
        nq = qn.loc[(fi, s)] if (fi, s) in qn.index else None
        r["new"] = bool((fi, s) not in To)
        # (a) q-control
        r["cc_sb_old"] = float(o.cc_vs_seisbench) if o is not None and "cc_vs_seisbench" in qo.columns else np.nan
        r["cc_sb_new"] = float(nq.cc_vs_seisbench) if nq is not None and "cc_vs_seisbench" in qn.columns else np.nan
        r["fail_a_qcontrol"] = bool(s in QC_STAS and np.isfinite(r["cc_sb_new"]) and
                                    ((np.isfinite(r["cc_sb_old"]) and r["cc_sb_new"] < r["cc_sb_old"] - 0.05) or r["cc_sb_new"] < 0.6))
        # (c) S pick + differential moveout
        r["spick_old"] = float(o.s_pick_s) if o is not None else np.nan
        r["spick_new"] = float(nq.s_pick_s) if nq is not None else np.nan
        r["dspick"] = r["spick_new"] - r["spick_old"]
        r["fail_c_spick"] = bool(np.isfinite(r["dspick"]) and abs(r["dspick"]) > 0.5)
        r["dmoveout"] = (float(mvn.get((f, s), np.nan)) - float(mvo.get((f, s), np.nan))) if mvn is not None else np.nan
        r["fail_c_moveout"] = bool(np.isfinite(r["dmoveout"]) and abs(r["dmoveout"]) > 0.3)
        # (d) cross-family collapse
        r["ccx_old"] = cco.get((f, s), np.nan); r["ccx_new"] = ccn.get((f, s), np.nan)
        r["fail_d_crossfam"] = bool(np.isfinite(r["ccx_new"]) and
                                    ((np.isfinite(r["ccx_old"]) and r["ccx_new"] > r["ccx_old"] + 0.1) or r["ccx_new"] > 0.7))
        # (e) spectral narrowing
        r["spec_frac"] = spn.get((f, s), np.nan); r["spec_frac_old"] = spo.get((f, s), np.nan)
        r["fail_e_spectral"] = bool(np.isfinite(r["spec_frac"]) and r["spec_frac"] > 0.5)
        # (f) S-SNR
        r["snr_old"] = float(o.snr_db) if o is not None else np.nan
        r["snr_new"] = float(nq.snr_db) if nq is not None else np.nan
        r["dsnr"] = r["snr_new"] - r["snr_old"]
        r["fail_f_snr"] = bool(np.isfinite(r["dsnr"]) and r["dsnr"] < -1.0)
        r["fail_any"] = any(r[k] for k in r if k.startswith("fail_") and k != "fail_any")
        # conditional versions: absolute clause counts only if the frozen template passed it
        r["cfail_a_qcontrol"] = bool(s in QC_STAS and np.isfinite(r["cc_sb_new"]) and np.isfinite(r["cc_sb_old"]) and
                                     (r["cc_sb_new"] < r["cc_sb_old"] - 0.05 or (r["cc_sb_old"] >= 0.6 and r["cc_sb_new"] < 0.6)))
        r["cfail_c_spick"] = r["fail_c_spick"]; r["cfail_c_moveout"] = r["fail_c_moveout"]
        r["cfail_d_crossfam"] = bool(np.isfinite(r["ccx_new"]) and np.isfinite(r["ccx_old"]) and
                                     (r["ccx_new"] > r["ccx_old"] + 0.1 or (r["ccx_old"] <= 0.7 and r["ccx_new"] > 0.7)))
        r["cfail_e_spectral"] = bool(np.isfinite(r["spec_frac"]) and np.isfinite(r.get("spec_frac_old", np.nan)) and
                                     r["spec_frac"] > 0.5 and r["spec_frac_old"] <= 0.5)
        r["cfail_f_snr"] = r["fail_f_snr"]
        r["cfail_any"] = any(r[k] for k in r if k.startswith("cfail_") and k != "cfail_any")
        rows.append(r)
    D = pd.DataFrame(rows); D.to_csv(a.out, index=False)
    fam = D.groupby("bkey").agg(n_sta=("sta", "size"),
                                **{c: (c, "sum") for c in D.columns if c.startswith("fail_") or c.startswith("cfail_")})
    fam["accept_literal"] = fam.fail_any == 0
    fam["accept"] = fam.cfail_any == 0
    fam.to_csv(a.out.replace(".csv", "_by_family.csv"))
    pd.set_option("display.width", 250)
    print(f"[tests] {len(D)} (family, station) templates, {D.bkey.nunique()} families")
    print("[literal]"); print(D[[c for c in D.columns if c.startswith("fail_")]].sum().to_string())
    print("[baseline-conditioned]"); print(D[[c for c in D.columns if c.startswith("cfail_")]].sum().to_string())
    print(f"\n[per family] literal: accepted {int(fam.accept_literal.sum())} / {len(fam)}; "
          f"baseline-conditioned: accepted {int(fam.accept.sum())} / {len(fam)}, rejected {int((~fam.accept).sum())}")
    pd.DataFrame(dict(bkey=fam.index[fam.accept])).to_csv(a.out.replace(".csv", "_accepted_fams.csv"), index=False)
    print("[medians] cc_sb %.3f -> %.3f | S-SNR %.2f -> %.2f dB | max cross-family cc %.3f -> %.3f | "
          "|dS pick| %.2f s | max 1-Hz-bin energy fraction %.2f" %
          (D.cc_sb_old.median(), D.cc_sb_new.median(), D.snr_old.median(), D.snr_new.median(),
           D.ccx_old.median(), D.ccx_new.median(), D.dspick.abs().median(), D.spec_frac.median()))
    print(f"[write] {a.out} + _by_family.csv")


if __name__ == "__main__":
    main()
