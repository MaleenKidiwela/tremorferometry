#!/usr/bin/env python
"""DUC_FIX stage 1 -- KNOWN ANSWER FIRST.  His templates, his stations (bed A), his rules (hourly
8 x MAD, complete-hour channels, 1-s decluster, + his cc x nchannel catalogue filter) on the FIXED
data path, 2008-04-01..2008-05-10, scored one-to-one at +-2 s against his catalogue exactly as the
0.836 run was (t8_score.one_to_one / score).

Two variants, both reported (burst muting is a SCREEN, so it must pass the known answer):
  full      all fixes (flat/gap guard, dead floors, burst muting, line notches, exact sums, |cc|>1 zero)
  nomute    all fixes except burst muting
GATE: the variant carried forward must reach mean recall_all >= 0.836 (the old run).  If both pass,
the one with the higher pooled F1 under his cc x nchannel catalogue rule (ours vs his, both filtered) is
carried forward.  If neither passes,
results/duc_fix/STAGE1_FAILED is written and the queue stops.
"""
import os, sys, json, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from duc_fix_common import *
import duc_fix_clean as CL
import duc_fix_mf as MF
from t8_score import score

T0, T1 = "2008-04-01", "2008-05-11"          # days 04-01 .. 05-10 inclusive (as the 0.836 run)
OLD_RECALL = 0.836


def score_df(O, T):
    rows = []
    for f in sorted(set(T.fam) | set(O.fam)):
        oo = O[O.fam == f].t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
        tt = T[T.fam == f]
        r = score(oo, tt.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9, tt.cc.values, 2.0)
        r.update(fam=f); rows.append(r)
    D = pd.DataFrame(rows); D["bostock30"] = (D.precision >= 0.30).astype(int)
    return D


def summary(D, name):
    n_t, n_o, n_m = int(D.n_theirs.sum()), int(D.n_ours.sum()), int(D.n_match.sum())
    pr, pp = n_m / max(n_t, 1), n_m / max(n_o, 1)
    Dt = D[D.n_theirs > 0]
    return dict(variant=name, families=len(Dt), mean_recall_all=float(Dt.recall_all.mean()),
                median_recall_all=float(Dt.recall_all.median()), mean_recall_indep=float(Dt.recall_indep.mean()),
                mean_recall_loud=float(Dt.recall_loud.mean()), mean_precision=float(Dt.precision.mean()),
                pooled_recall=pr, pooled_precision=pp, pooled_f1=2 * pr * pp / max(pr + pp, 1e-9),
                n_theirs=n_t, n_ours=n_o, n_match=n_m, fam_bostock30=int(Dt.bostock30.sum()))


def main():
    os.makedirs(FIX, exist_ok=True)
    thr = thresholds()
    T = his_catalogue(); T = T[(T.t >= T0) & (T.t < T1)].copy()
    T["thr"] = T.fam.map(thr).fillna(0.0)
    Tf = T[(T.thr > 0) & (T.cc * T.nch >= T.thr)]
    days = pd.date_range(T0, T1, freq="D", inclusive="left")
    S = []
    for name, mute in [("full", True), ("nomute", False)]:
        CL.MUTE_ON = mute
        out = f"{FIX}/repro_{name}_mf.parquet"
        O = pd.read_parquet(out) if os.path.exists(out) else MF.run(days, out)
        O["fam"] = O.fam.astype(str); O = O[(O.t >= T0) & (O.t < T1)]
        D = score_df(O, T); D.to_csv(f"{FIX}/repro_{name}_score.csv", index=False)
        S.append(summary(D, name))
        Of = O[O.pass_his == 1]
        Df = score_df(Of, Tf); Df.to_csv(f"{FIX}/repro_{name}_score_hisrule.csv", index=False)
        S.append(summary(Df, name + "+his_cc_x_nch_rule"))
    # the old 12-h-MAD, unfixed run, re-scored by the same code for the table
    Oo = pd.read_parquet(f"{RES}/duc_repro_bedA_detections.parquet"); Oo["fam"] = Oo.fam.astype(str)
    Oo = Oo[(Oo.t >= T0) & (Oo.t < T1)]
    S.append(summary(score_df(Oo, T), "old_12h_unfixed"))
    Oo["thr"] = Oo.fam.map(thr).fillna(0.0)
    S.append(summary(score_df(Oo[(Oo.thr > 0) & (Oo.cc * Oo.nchannel >= Oo.thr)], Tf), "old_12h_unfixed+his_cc_x_nch_rule"))
    SS = pd.DataFrame(S); SS.to_csv(f"{FIX}/repro_summary.csv", index=False)
    print(SS.to_string(index=False))
    ok = SS[SS.variant.isin(["full", "nomute"]) & (SS.mean_recall_all >= OLD_RECALL)]
    if not len(ok):
        open(f"{FIX}/STAGE1_FAILED", "w").write(SS.to_string(index=False))
        print("[stage1] FAILED: no variant reaches mean recall_all >= 0.836"); sys.exit(3)
    # choose by the pooled F1 of the catalogue we deliver (his cc x nchannel rule on both sides): the raw pooled
    # numbers are dominated by one of his six ELIMINATED families (080328.09.029, 10,031 raw detections in 40 days)
    hr = SS.set_index("variant")
    ok = ok.assign(f1_rule=[hr.loc[v + "+his_cc_x_nch_rule", "pooled_f1"] for v in ok.variant])
    best = ok.sort_values("f1_rule", ascending=False).iloc[0]
    ch = dict(variant=best.variant, mute=int(best.variant == "full"), mean_recall_all=float(best.mean_recall_all),
              pooled_f1=float(best.pooled_f1), pooled_precision=float(best.pooled_precision), f1_his_rule=float(best.f1_rule),
              rule="carry the variant with mean recall_all >= 0.836 (raw, as the old run) and the higher pooled F1 under his catalogue rule",
              written=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    json.dump(ch, open(f"{FIX}/detector_choice.json", "w"), indent=1)
    print("[stage1] PASS ->", ch)
    loo_check(ch["mute"] == 1)


LOO = ["LBP", "KCS", "KHBB", "B039", "WDC", "LSF"]
def loo_check(mute):
    """Known-answer check of the RE-TEMPLATING procedure (not a gate): run the carried detector over
    the same 40 days WITHOUT six of his stations, re-template them from those detections exactly as a
    new epoch would be, and compare with HIS template on the same channel (cc_his) and with 20
    scrambled controls (S vs ctrl_max)."""
    import duc_fix_retemplate as RT
    CL.MUTE_ON = mute
    out = f"{FIX}/loo_mf.parquet"
    days = pd.date_range(T0, T1, freq="D", inclusive="left")
    O = pd.read_parquet(out) if os.path.exists(out) else MF.run(days, out, exclude=tuple(LOO))
    O["fam"] = O.fam.astype(str)
    E = epochs(); fp = families(); HT = his_templates()
    dur = {r.family: float(r.duration) for r in fp.itertuples()}
    famst = {r.family: r.stations.split(",") for r in fp.itertuples()}
    rows = []
    for sta in LOO:
        ep = E[(E.sta == sta) & (E.start <= pd.Timestamp(T0)) & (E.end > pd.Timestamp(T0))].iloc[0].to_dict()
        ep["start"], ep["end"] = pd.Timestamp(T0), pd.Timestamp(T1)
        fl = [f for f, s in famst.items() if sta in s]
        C = RT.candidates(O[O.fam.isin(fl)], ep["start"], ep["end"])
        Cs = pd.concat([C[C.fam == f].head(3 * RT.NMAX) for f in fl]) if len(C) else C
        ht = {(f, c[-1]): v for f in fl for (s2, c), v in HT.get(f, {}).items() if s2 == sta}
        r, _ = RT.build_epoch((ep, fl, dur, ep["end"], ht, Cs))
        rows += r
    R = pd.DataFrame(rows); R.to_csv(f"{FIX}/loo_retemplate_check.csv", index=False)
    Rk = R[R.N >= RT.NMIN]
    print(f"[stage1 loo] {len(R)} (family,channel): N>={RT.NMIN}: {len(Rk)}; kept {int(R.keep.sum())}; "
          f"cc_his median {Rk.cc_his.median():.3f} (kept {R[R.keep == 1].cc_his.median():.3f}, "
          f"dropped {R[(R.keep == 0) & (R.N >= RT.NMIN)].cc_his.median():.3f}); S/ctrl_max median {(Rk.S / Rk.ctrl_max).median():.2f}")


if __name__ == "__main__":
    main()
