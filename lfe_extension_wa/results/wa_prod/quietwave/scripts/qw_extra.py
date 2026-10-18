#!/usr/bin/env python
"""Quiet-day waveform test, step 3: selection-free evidence per candidate + a rough real-fraction estimate.

The in-template correlation cannot tell real from false (every candidate, real or not, was SELECTED for matching the
template there).  Two measures the selection does not touch:
  out_cc   correlation of the single window OUTSIDE the template window ([-8, 0) s: the P wave; [12, 24) s: coda) with an
           independent long reference stack (tremor-day data P >= 0.85, even-parity days only); coef-weighted over
           the candidate's stations.  Mean over a set is linear in the share of real events.
  e_snr    window energy / pre-window energy ([0, 12) s vs [-20, -8) s), station median, minus 1.
Real fraction (rough) = family-matched ratio of the quiet-set mean to the tremor-day data mean (odd-parity days), with
a day-block bootstrap; both measures assume quiet-day LFEs are as strong as tremor-day LFEs of the same family, so they
are rough: weaker inter-ETS events push both DOWN, quieter pre-windows on quiet days push e_snr UP.
Also: the share of candidates with a transient-like window (network-max envelope peak/median > 12 or a spike)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qw_analyze as Q

Q.build_ref()
REF = Q.REFST
FS = Q.FS; PRE = Q.PRE; n = Q.n
segs = [slice(PRE - int(8 * FS), PRE), slice(PRE + n, PRE + n + int(12 * FS))]
W = Q.WIN[(Q.WIN.widx >= 0)].copy()
W = W[[(f, s, t) in REF for f, s, t in zip(W.fid, W.sta, W.tkey)]]
oc = np.full(len(W), np.nan)
for k, r in enumerate(W.itertuples()):
    x = Q.wave(r.date, r.widx)
    if not np.isfinite(x).all():
        continue
    ref = REF[(r.fid, r.sta, r.tkey)]
    ch = Q.TPL[(r.fid, r.sta, r.tkey)]["live"] & (np.abs(x).sum(1) > 0)
    a = np.concatenate([x[ch][:, s].ravel() for s in segs]); b = np.concatenate([ref[ch][:, s].ravel() for s in segs])
    den = np.sqrt(np.dot(a, a) * np.dot(b, b))
    oc[k] = np.dot(a, b) / den if den > 0 else np.nan
W["out_cc"] = oc
W = W[np.isfinite(W.out_cc)]
C = W.groupby("row_id").apply(lambda x: pd.Series(dict(out_cc=np.average(x.out_cc, weights=x.coef), n_ref_sta=len(x)))).reset_index()
net = pd.read_parquet(f"{Q.TB}/t2_network_per_candidate.parquet")
C = C.merge(net[["row_id", "group", "series", "P_own_WAVE", "catalogue_era", "date", "fid", "snrw_med", "envpk_max", "spike"]], on="row_id")
C["e_snr"] = C.snrw_med.clip(0, 3) - 1            # clipped: gaps make a few pre-windows ~0
C["transient"] = (C.envpk_max > 12) | C.spike.astype(bool)
C.to_parquet(f"{Q.TB}/t4_selection_free_per_candidate.parquet", index=False)

odd = ~np.array([Q.refday(d) for d in C.date], bool)
TM = {}
for cut in Q.CUTS:                     # the real-event reference at the SAME cut (stronger events score higher P)
    trem = C[(C.group == "tremor") & (C.series == "data") & (C.P_own_WAVE >= cut) & odd]
    tcount = trem.groupby("fid").size()
    TM[cut] = trem.groupby("fid")[["out_cc", "e_snr"]].mean()[tcount >= 10]
rng = np.random.default_rng(11)


def frac(h, col, tm):
    h = h[h.fid.isin(tm.index)]
    if len(h) < 5:
        return np.nan, np.nan, np.nan, len(h)
    den = tm.loc[h.fid, col].to_numpy()
    num = h[col].to_numpy()
    est = num.sum() / den.sum()
    days = h.date.to_numpy(); ud = np.unique(days); idx = {d: np.flatnonzero(days == d) for d in ud}
    bs = []
    for _ in range(500):
        pick = np.concatenate([idx[d] for d in rng.choice(ud, len(ud))])
        bs.append(num[pick].sum() / den[pick].sum())
    return est, np.quantile(bs, 0.05), np.quantile(bs, 0.95), len(h)


rows = []
for group in ("quiet", "tremor"):
    for era in (["all", "catalogue", "no_catalogue"] if group == "quiet" else ["all"]):
        for cut in Q.CUTS:
            for ser in Q.SER:
                m = (C.group == group) & (C.series == ser) & (C.P_own_WAVE >= cut)
                if group == "tremor":
                    m &= odd
                if era != "all":
                    m &= C.catalogue_era == (era == "catalogue")
                h = C[m]
                f1 = frac(h, "out_cc", TM[cut]); f2 = frac(h, "e_snr", TM[cut])
                rows.append(dict(group=group, era=era, cut=cut, series=ser, n=len(h), out_cc_mean=h.out_cc.mean(),
                                 out_cc_sem=h.out_cc.std() / np.sqrt(max(len(h), 1)), e_snr_mean=h.e_snr.mean(),
                                 transient_share=h.transient.mean(),
                                 real_frac_outcc=f1[0], real_frac_outcc_p05=f1[1], real_frac_outcc_p95=f1[2],
                                 real_frac_esnr=f2[0], real_frac_esnr_p05=f2[1], real_frac_esnr_p95=f2[2], n_fam_matched=f1[3]))
R = pd.DataFrame(rows)
R.to_csv(f"{Q.TB}/t4_selection_free_summary.csv", index=False)
# transient-like data candidates: do they carry LFE evidence outside the template window?
rows = []
for group in ("quiet", "tremor"):
    for era in (["all", "catalogue", "no_catalogue"] if group == "quiet" else ["all"]):
        for cut in Q.CUTS:
            for tr in (False, True):
                m = (C.group == group) & (C.series == "data") & (C.P_own_WAVE >= cut) & (C.transient == tr)
                if group == "tremor":
                    m &= odd
                if era != "all":
                    m &= C.catalogue_era == (era == "catalogue")
                h = C[m]
                f1 = frac(h, "out_cc", TM[cut])
                rows.append(dict(group=group, era=era, cut=cut, transient=tr, n=len(h), out_cc_mean=h.out_cc.mean(),
                                 real_frac_outcc=f1[0], p05=f1[1], p95=f1[2]))
TT = pd.DataFrame(rows)
TT.to_csv(f"{Q.TB}/t4_transient_split.csv", index=False)
print(TT.round(3).to_string())
pd.set_option("display.width", 250)
print(R.round(3).to_string())
