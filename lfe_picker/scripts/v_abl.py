#!/usr/bin/env python
"""P4-0b SENSITIVITY ABLATION -- NOT THE FROZEN RESULT.

The frozen P4-0b given-seed arm returned 0/18 assembled and THAT VERDICT STANDS.  This module asks
one mechanical question about it: the pre-registered survival floor of >= 20 independent hits was
set for a multi-year record and was applied to a 203-day one, where pass-0 single-event templates
give a median of 6 independent hits -- so the floor removed 49 of 54 seeds BEFORE any template was
iterated.  Two lower floors are run here to bracket how much of the 0/18 is the floor and how much
is the construction.

  A1  pass-0 survival floor >= 8 independent hits (>= 3 days, >= 2 episodes, <= 20 % near seed,
      N-moveout ratio >= 3 all kept) -> iterate twice -> dedupe -> score with the IDENTICAL
      relative clauses of the frozen run.
  A2  floor >= 5 and the day clause relaxed to >= 2 days -- the lowest floor at which the
      N-moveout ratio can still be evaluated -- to bracket from below.
  A3  DIAGNOSTICS ONLY, no GPU: per pass-0 template, the fraction of its independent hits that
      match its OWN family's hp events (one-to-one +-2 s), so a starved-but-pure template (floor
      problem) can be told from an already-impure one (seed/statistic problem); plus the
      per-iteration purity trajectory of the five frozen survivors.
  A4  the reference check: recall/purity of the truth-stack references at the frozen operating
      point's T_core = 8.75 as well as at the 8.0 floor the frozen run used, so the low reference
      recall (0.11-0.73) is understood.

NOTHING under results/p4/p4_0b/ is written or modified.  Everything lands in
results/p4/p4_0b_abl/.  Seeds, templates, T rule, iteration count (exactly two), dedupe rule and
the relative scoring clauses are the frozen ones, imported unchanged.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import v_common as V                                                        # noqa: E402
import s_common as S                                                        # noqa: E402
import v_iterate as IT                                                      # noqa: E402
import v_score as SC                                                        # noqa: E402

FRZ = V.OUT                                   # results/p4/p4_0b -- READ ONLY
ABL = f"{V.RES}/p4/p4_0b_abl"
PREREG = f"{V.RES}/p4/p4_0b_abl_prereg.json"
LABEL = "SENSITIVITY ABLATION -- not the frozen result"
CONF = dict(A1=dict(floor=8, min_days=3), A2=dict(floor=5, min_days=2))


# ------------------------------------------------------------------ survival --------------------
def _indep_of(H, tids_col="tid", thr=None):
    """{tid: DataFrame of independent hits} from a hits table, optionally re-thresholded."""
    out = {}
    if len(H) == 0:
        return out
    if thr is not None:
        H = H[H.net_mad >= thr]
    for tid, g in H.groupby(tids_col):
        g = g.sort_values("t")
        out[tid] = g[V.independent(g.t.to_numpy())].copy()
    return out


def survive(label, level):
    cfg = CONF[label]
    ED = V.episode_days()
    Sd = pd.read_parquet(f"{FRZ}/seeds_given.parquet").set_index("seed")
    if level == 0:
        H = pd.read_parquet(f"{FRZ}/run_given0_hits.parquet")
        NH = pd.read_parquet(f"{FRZ}/run_given_nm_hits.parquet")
        NM = {t: int(len(g)) for t, g in _indep_of(NH).items()}
        out = f"{ABL}/{label}_run0"
    else:
        H = pd.read_parquet(f"{ABL}/{label}_run_it{level}_hits.parquet")
        NM = None
        out = f"{ABL}/{label}_run_it{level}"
    IX = _indep_of(H)
    tids = sorted(set(H.tid)) if len(H) else []
    rows, keepdf = [], []
    for tid in tids:
        gi = IX[tid]
        eps = {ED.get(pd.Timestamp(d), "") for d in gi.date} - {""}
        t_seed = float(Sd.t[tid]) if tid in Sd.index else np.nan
        near = float(np.mean(np.abs(gi.t.to_numpy() - t_seed) <= V.NEAR_SEED_S)) \
            if np.isfinite(t_seed) else np.nan
        nm = NM.get(tid, 0) if NM is not None else np.nan
        ratio = (len(gi) / max(nm, 1)) if NM is not None else np.nan
        surv = bool(len(gi) >= cfg["floor"] and gi.date.nunique() >= cfg["min_days"]
                    and len(eps) >= V.SURV_MIN_EPISODES
                    and (not np.isfinite(near) or near <= V.SURV_MAX_NEAR_SEED)
                    and (NM is None or ratio >= V.NMOVEOUT_RATIO))
        rows.append(dict(tpl=tids.index(tid), tid=tid, n_hits=int((H.tid == tid).sum()),
                         n_indep=int(len(gi)), n_days=int(gi.date.nunique()),
                         n_episodes=int(len(eps)), frac_near_seed=near,
                         n_moveout=nm, moveout_ratio=ratio, survives=surv))
        gi = gi.copy(); gi["tid"] = tid
        keepdf.append(gi)
    R = pd.DataFrame(rows)
    R.to_csv(out + "_survival.csv", index=False)
    if keepdf:
        pd.concat(keepdf, ignore_index=True).to_parquet(out + "_indep.parquet", index=False)
    V.log(f"[{label}/L{level}] floor>={cfg['floor']} days>={cfg['min_days']}: "
          f"{int(R.survives.sum()) if len(R) else 0}/{len(R)} survive "
          f"({R.tid.map(lambda t: t.split('_')[0]).nunique() if len(R) else 0} seed-families seen)")
    return R


# ------------------------------------------------------------------ iteration -------------------
def iterate_abl(bank, indep, survivors, out_bank, out_diag, use_all, min_indep, min_stack):
    """VERBATIM copy of v_iterate.iterate with the stack-count floor exposed as `min_stack`.

    v_iterate.py is hashed in the FROZEN freeze and is not edited; the only change here is that the
    hard-coded `len(g) < 10` stacking floor becomes a parameter, because A2's survival floor is 5."""
    tids, stas, W, off, meta = V.load_bank(bank)
    Hs = pd.read_parquet(indep)
    SV = pd.read_csv(survivors)
    keep = (SV[SV.n_indep >= min_indep].tid.tolist() if use_all else SV[SV.survives].tid.tolist())
    tid2i = {t: i for i, t in enumerate(tids)}
    n_cut = int(round(IT.CUT_LEN_S * V.FS))
    nom = int(round(IT.CUT_PRE_S * V.FS))
    n = int(round(V.TPL1_LEN_S * V.FS))
    idx = V.cache_idx()
    pre_in = float(meta.get("pre_s", V.TPL0_PRE_S))
    Wn = np.full((len(keep), len(stas), 3, n), np.nan, np.float32)
    offn = np.full((len(keep), len(stas)), np.nan)
    diag, tids_out = [], []
    for tid in keep:
        src = tid2i[tid]
        g = Hs[Hs.tid == tid].sort_values("net_mad", ascending=False).head(IT.MAX_STACK)
        if len(g) < min_stack:
            continue
        times = g.t.to_numpy(float)
        stacks, picks, snr = {}, {}, {}
        for si, s in enumerate(stas):
            o = off[src, si]
            if not np.isfinite(o):
                continue
            segs = IT._cut_hits(s, times, o + pre_in - IT.CUT_PRE_S, idx, n_cut, IT.CUT_PRE_S)
            if segs is None or len(segs) < min_stack:
                continue
            st, _, _ = S.align_stack(segs, iters=V.ALIGN_ITERS, lag=V.ALIGN_LAG)
            i_s = S.pick_s(st, nom, search_s=4.0)
            stacks[si], picks[si] = st, i_s
            snr[si] = float(S.s_snr_db(st, i_s))
            diag.append(dict(tid=tid, sta=stas[si], n_stacked=len(segs), s_snr_db=snr[si],
                             s_excess_db=float(S.s_excess_db(st, i_s)),
                             pick_shift_s=float((i_s - nom) / V.FS)))
        if not stacks:
            continue
        best = max(snr, key=lambda k: (snr[k] if np.isfinite(snr[k]) else -1e9))
        dt = (picks[best] - nom) / V.FS
        bi2 = len(tids_out)
        got = 0
        for si, st in stacks.items():
            a = nom + int(round(dt * V.FS)) - int(round(V.TPL1_PRE_S * V.FS))
            if a < 0 or a + n > st.shape[1]:
                continue
            Wn[bi2, si] = V.norm_channels(st[:, a:a + n])
            offn[bi2, si] = off[src, si] + pre_in + dt - V.TPL1_PRE_S
            got += 1
        if got == 0:
            continue
        tids_out.append(tid)
        diag.append(dict(tid=tid, sta="__global__", n_stacked=len(g), s_snr_db=snr[best],
                         s_excess_db=np.nan, pick_shift_s=float(dt)))
    Wn, offn = Wn[:len(tids_out)], offn[:len(tids_out)]
    ok = np.isfinite(Wn).all(axis=(2, 3))
    nsite = np.array([len({V.SITE_OF[stas[i]] for i in np.nonzero(ok[b])[0]})
                      for b in range(len(tids_out))])
    sel = np.nonzero(nsite >= V.MIN_SITES)[0]
    pd.DataFrame(diag).to_csv(out_diag, index=False)
    V.save_bank(out_bank, [tids_out[i] for i in sel], stas, Wn[sel], offn[sel],
                dict(kind="iterated_ABLATION", label=LABEL, len_s=V.TPL1_LEN_S,
                     pre_s=V.TPL1_PRE_S, parent=os.path.basename(bank), n_in=len(keep)))
    V.log(f"[iterate/ABL] {len(keep)} in -> {len(sel)} iterated -> {out_bank}")


# ------------------------------------------------------------------ scoring ---------------------
def _reference(thr=None):
    """{bkey: (recall_ref, purity_ref, n_ref_indep)} from the FROZEN truth-stack reference run."""
    keep, EV = SC._hp_by_family()
    H = pd.read_parquet(f"{FRZ}/run_ref_hits.parquet")
    IX = _indep_of(H, thr=thr)
    out = {}
    for b in keep:
        rt = np.sort(IX.get(f"REF{b}", pd.DataFrame(dict(t=[]))).t.to_numpy(float))
        ev = EV.get(int(b), np.zeros(0))
        mr, _ = V.one2one(ev, rt, V.MATCH_TOL_S)
        mp, _ = V.one2one(rt, ev, V.MATCH_TOL_S)
        out[int(b)] = (float(mr.mean()) if len(ev) else np.nan,
                       float(mp.mean()) if len(rt) else np.nan, int(len(rt)))
    return keep, EV, out


def _score_level(indep_path, keep, EV, REF, level, label, dedupe=None):
    if not os.path.exists(indep_path):
        return pd.DataFrame()
    G = pd.read_parquet(indep_path)
    if dedupe is not None and len(dedupe):
        G = G[G.tid.isin(set(dedupe[dedupe.kept].tid))]
    rows = []
    for b in keep:
        ev = EV.get(int(b), np.zeros(0))
        rec_ref, pur_ref, n_ref = REF[int(b)]
        no_ref = (not np.isfinite(rec_ref)) or rec_ref <= 1e-9
        best = None
        for tid in sorted(set(G.tid)):
            if not str(tid).startswith(f"G{b}_"):
                continue
            ht = np.sort(G[G.tid == tid].t.to_numpy(float))
            if len(ht) == 0:
                continue
            m1, _ = V.one2one(ev, ht, V.MATCH_TOL_S)
            m2, _ = V.one2one(ht, ev, V.MATCH_TOL_S)
            rec, pur = float(m1.mean()) if len(ev) else 0.0, float(m2.mean())
            ok = (not no_ref) and rec >= V.REL_RECALL * rec_ref and pur >= V.REL_PURITY * pur_ref
            r = dict(tid=tid, n_hits=len(ht), recall=rec, purity=pur, assembled=ok)
            if best is None or (ok, rec) > (best["assembled"], best["recall"]):
                best = r
        row = dict(label=label, level=level, bkey=int(b), n_hp=len(ev), n_ref_indep=n_ref,
                   recall_ref=rec_ref, purity_ref=pur_ref, no_reference=bool(no_ref),
                   need_recall=V.REL_RECALL * rec_ref, need_purity=V.REL_PURITY * pur_ref)
        row.update(best or dict(tid=None, n_hits=0, recall=0.0, purity=np.nan, assembled=False))
        rows.append(row)
    return pd.DataFrame(rows)


def score(label):
    keep, EV, REF = _reference()
    dd = pd.read_csv(f"{ABL}/{label}_dedupe.csv") if os.path.exists(f"{ABL}/{label}_dedupe.csv") \
        else None
    parts = [_score_level(f"{ABL}/{label}_run0_indep.parquet", keep, EV, REF, 0, label),
             _score_level(f"{ABL}/{label}_run_it1_indep.parquet", keep, EV, REF, 1, label),
             _score_level(f"{ABL}/{label}_run_it2_indep.parquet", keep, EV, REF, 2, label, dd)]
    T = pd.concat([p for p in parts if len(p)], ignore_index=True)
    T.to_csv(f"{ABL}/{label}_score_by_iteration.csv", index=False)
    fin = T[T.level == 2]
    summ = dict(label=LABEL, ablation=label, clauses=CONF[label],
                n_eligible=len(keep),
                assembled_by_level={int(l): int(g.assembled.sum()) for l, g in T.groupby("level")},
                families_reaching_iteration=int(T[T.level == 1].n_hits.gt(0).sum()),
                n_assembled=int(fin.assembled.sum()) if len(fin) else 0,
                frozen_result="0/18 assembled -- THE FROZEN VERDICT, unchanged by this ablation",
                near_miss=[dict(bkey=int(r.bkey), recall=float(r.recall),
                                need_recall=float(r.need_recall), purity=float(r.purity)
                                if np.isfinite(r.purity) else None,
                                need_purity=float(r.need_purity))
                           for r in fin.itertuples()
                           if np.isfinite(r.recall) and r.recall >= r.need_recall
                           and not r.assembled] if len(fin) else [])
    json.dump(summ, open(f"{ABL}/{label}_score.json", "w"), indent=1, default=str)
    V.log(f"[{label}] assembled by level {summ['assembled_by_level']} "
          f"(frozen: 0/18); near-misses {len(summ['near_miss'])}")
    return summ


# ------------------------------------------------------------------ A3 / A4 ---------------------
def diag():
    """A3 -- no GPU.  Purity of every pass-0 template against its OWN family, and the per-iteration
    purity trajectory of the five frozen survivors."""
    keep, EV = SC._hp_by_family()
    rows = []
    for lvl, p in ((0, f"{FRZ}/run_given0_indep.parquet"), (1, f"{FRZ}/run_given_it1_indep.parquet"),
                   (2, f"{FRZ}/run_given_it2_indep.parquet")):
        if not os.path.exists(p):
            continue
        G = pd.read_parquet(p)
        for tid, g in G.groupby("tid"):
            b = int(str(tid)[1:].split("_")[0])
            ev = EV.get(b, np.zeros(0))
            ht = np.sort(g.t.to_numpy(float))
            m2, _ = V.one2one(ht, ev, V.MATCH_TOL_S)
            m1, _ = V.one2one(ev, ht, V.MATCH_TOL_S)
            rows.append(dict(level=lvl, tid=tid, bkey=b, n_indep=len(ht),
                             n_matched_own=int(m2.sum()), purity_own=float(m2.mean()),
                             recall_own=float(m1.mean()) if len(ev) else np.nan,
                             n_hp=len(ev)))
    D = pd.DataFrame(rows)
    D.to_csv(f"{ABL}/A3_template_purity.csv", index=False)
    P0 = D[D.level == 0]
    frozen_surv = set(pd.read_csv(f"{FRZ}/run_given0_survival.csv").query("survives").tid)
    traj = D[D.tid.isin(frozen_surv)].pivot_table(index="tid", columns="level",
                                                  values=["n_indep", "purity_own", "recall_own"])
    traj.to_csv(f"{ABL}/A3_frozen_survivor_trajectory.csv")
    bins = [(0, 4), (5, 7), (8, 12), (13, 19), (20, 10 ** 6)]
    by = []
    for lo, hi in bins:
        g = P0[(P0.n_indep >= lo) & (P0.n_indep <= hi)]
        by.append(dict(n_indep_bin=f"{lo}-{hi if hi < 10**6 else '+'}", n_templates=len(g),
                       median_purity_own=float(g.purity_own.median()) if len(g) else np.nan,
                       mean_purity_own=float(g.purity_own.mean()) if len(g) else np.nan,
                       n_families=int(g.bkey.nunique())))
    B = pd.DataFrame(by)
    B.to_csv(f"{ABL}/A3_purity_by_hitcount.csv", index=False)
    summ = dict(label=LABEL, question="were the low-hit pass-0 templates PURE but starved (a floor "
                                      "problem) or already IMPURE (a seed/statistic problem)?",
                n_pass0_templates=int(len(P0)),
                median_n_indep=float(P0.n_indep.median()),
                p75_n_indep=float(P0.n_indep.quantile(0.75)),
                max_n_indep=int(P0.n_indep.max()),
                median_purity_own_all=float(P0.purity_own.median()),
                purity_by_hitcount=by,
                spearman_purity_vs_hits=float(P0[["n_indep", "purity_own"]].corr(
                    method="spearman").iloc[0, 1]) if len(P0) > 3 else None)
    json.dump(summ, open(f"{ABL}/A3_diagnostics.json", "w"), indent=1, default=str)
    V.log(f"[A3] {len(P0)} pass-0 templates; median purity_own "
          f"{summ['median_purity_own_all']:.3f}; Spearman(purity, n_indep) "
          f"{summ['spearman_purity_vs_hits']}")
    print(B.to_string(index=False))
    return summ


def refcheck():
    """A4 -- the reference at the frozen operating point's T_core = 8.75 as well as at 8.0."""
    fz = json.load(open("/home/jovyan/tremorferometry/lfe_extension/results/"
                        "frozen_operating_point.json"))
    tcore, text = float(fz["T_core"]), float(fz["T_ext"])
    rows = []
    for thr, name in ((None, f"T=8.0 (the frozen P4-0b FA-only value)"),
                      (text, f"T_ext={text}"), (tcore, f"T_core={tcore}")):
        keep, EV, REF = _reference(thr=thr)
        for b in keep:
            r, p, n = REF[int(b)]
            rows.append(dict(threshold=name, T=thr if thr else 8.0, bkey=int(b),
                             n_ref_indep=n, n_hp=len(EV.get(int(b), [])),
                             recall_ref=r, purity_ref=p,
                             need_recall=V.REL_RECALL * r, need_purity=V.REL_PURITY * p))
    R = pd.DataFrame(rows)
    R.to_csv(f"{ABL}/A4_reference_by_threshold.csv", index=False)
    summ = dict(label=LABEL,
                question="the frozen reference recall spans 0.11-0.73; is that the reference being "
                         "weak, or the 8.0 floor admitting reference false alarms?",
                by_threshold={str(t): dict(median_recall=float(g.recall_ref.median()),
                                           median_purity=float(g.purity_ref.median()),
                                           median_n_indep=float(g.n_ref_indep.median()))
                              for t, g in R.groupby("threshold")},
                frozen_operating_point=dict(T_core=tcore, T_ext=text,
                                            own_recall_at_FA1=fz["T_ext_report"]["recall_own_test"],
                                            own_precision=fz["T_ext_report"]["precision_own_all"]))
    json.dump(summ, open(f"{ABL}/A4_reference.json", "w"), indent=1, default=str)
    V.log(f"[A4] {json.dumps(summ['by_threshold'], default=str)}")
    return summ


# ------------------------------------------------------------------ prereg / report -------------
def prereg():
    os.makedirs(ABL, exist_ok=True)
    frozen = json.load(open(f"{V.RES}/p4/p4_0b_freeze.json"))
    here = {os.path.basename(p): V.sha256_file(p) for p in V.v_scripts()}
    unchanged = {k: (here.get(k) == v) for k, v in frozen["scripts"].items()}
    d = dict(
        created_utc=str(pd.Timestamp.now("UTC")),
        label=LABEL,
        status="SENSITIVITY ABLATION of the frozen P4-0b given-seed arm.  The FROZEN VERDICT "
               "(0/18 assembled under the pre-registered clauses) STANDS and is reported first; "
               "nothing here can change it.  This run exists only to separate the effect of ONE "
               "pre-registered parameter -- the >= 20-independent-hit survival floor, set for a "
               "multi-year record and applied to a 203-day one -- from the construction itself.",
        motivation=dict(
            frozen_result="0/18 assembled",
            mechanics="pass-0 single-event templates give median 6 independent hits (p75 11.5, "
                      "max 34) over the 203-day record, so the >= 20 floor removed 49/54 seeds "
                      "before any iteration; of the 5 that iterated, G44 and G75 reached "
                      "recall >= 0.8x the truth-stack reference and failed only PURITY "
                      "(G44 0.298 vs 0.341 needed)"),
        unchanged_from_frozen=["seeds (results/p4/p4_0b/seeds_given.parquet)",
                               "pass-0 templates and their moveouts (bank_given0.npz)",
                               "pass-0 and N-moveout hit tables (no re-growth at level 0)",
                               "the T rule (FA-only, <= 0.1 per template-day on the 20 quiet days)",
                               "iteration count: exactly two, no third",
                               "the dedupe rule",
                               "the relative scoring clauses (0.5x recall, 0.8x purity)",
                               "the truth-stack reference run (run_ref_*)"],
        changed=dict(
            A1=dict(survival_floor_independent_hits=8, min_days=3, min_episodes=2,
                    max_frac_near_seed=0.20, nmoveout_ratio=3.0,
                    note="only the >= 20 floor moves; every other survival clause is the frozen one"),
            A2=dict(survival_floor_independent_hits=5, min_days=2, min_episodes=2,
                    max_frac_near_seed=0.20, nmoveout_ratio=3.0,
                    note="the lowest floor at which the N-moveout ratio can still be evaluated; "
                         "brackets from below"),
            stacking_floor="v_iterate's hard-coded `len(hits) < 10` stack floor is exposed as a "
                           "parameter in v_abl.iterate_abl (a verbatim copy) and set to the "
                           "ablation's own floor; v_iterate.py itself is NOT edited"),
        A3=dict(kind="diagnostics only, no GPU",
                what="per pass-0 template, the fraction of its independent hits matching its OWN "
                     "family's hp events (one-to-one +-2 s), binned by hit count, plus the "
                     "per-iteration purity trajectory of the five frozen survivors",
                reads="whether the low-hit templates were PURE but starved (floor problem) or "
                      "already IMPURE (seed/statistic problem)"),
        A4=dict(kind="post-processing only, no GPU",
                what="recall/purity of the truth-stack references at the frozen operating point's "
                     "T_ext = 8.25 and T_core = 8.75 as well as at the 8.0 floor the frozen run "
                     "used, to explain the 0.11-0.73 spread in reference recall"),
        outputs=dict(dir=ABL, done=f"{V.RES}/p4/P4_0B_ABL_DONE",
                     failed=f"{V.RES}/p4/P4_0B_ABL_FAILED", log="logs/p4_0b_abl_queue.log"),
        isolation="writes ONLY under results/p4/p4_0b_abl/ and the two markers; results/p4/p4_0b/ "
                  "is opened read-only; the frozen free-arm queue keeps running and the GPU is "
                  "shared through v_engine.wait_gpu",
        scripts=here,
        frozen_scripts_unchanged=unchanged,
        frozen_freeze_sha256=V.sha256_file(f"{V.RES}/p4/p4_0b_freeze.json"))
    json.dump(d, open(PREREG, "w"), indent=1, default=str)
    print(f"[prereg] -> {PREREG}")
    print("  frozen v_ scripts unchanged: " + str(all(unchanged.values())))


def report():
    out = dict(label=LABEL, frozen_verdict="0/18 assembled (stands)")
    for f in sorted(os.listdir(ABL)):
        if f.endswith(".json"):
            out[f[:-5]] = json.load(open(f"{ABL}/{f}"))
    json.dump(out, open(f"{ABL}/P4_0B_ABL_SUMMARY.json", "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in out.items() if not k.startswith("A3_template")},
                     indent=1, default=str)[:5000])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["prereg", "survive", "iterate", "score", "diag", "refcheck",
                                    "report"])
    ap.add_argument("--label", default="A1")
    ap.add_argument("--level", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(ABL, exist_ok=True)
    if a.cmd == "prereg":
        prereg()
    elif a.cmd == "survive":
        survive(a.label, a.level)
    elif a.cmd == "iterate":
        cfg = CONF[a.label]
        if a.level == 1:
            iterate_abl(f"{FRZ}/bank_given0.npz", f"{ABL}/{a.label}_run0_indep.parquet",
                        f"{ABL}/{a.label}_run0_survival.csv",
                        f"{ABL}/{a.label}_bank_it1.npz", f"{ABL}/{a.label}_it1_diag.csv",
                        use_all=False, min_indep=cfg["floor"], min_stack=cfg["floor"])
        else:
            iterate_abl(f"{ABL}/{a.label}_bank_it1.npz",
                        f"{ABL}/{a.label}_run_it1_indep.parquet",
                        f"{ABL}/{a.label}_run_it1_survival.csv",
                        f"{ABL}/{a.label}_bank_it2.npz", f"{ABL}/{a.label}_it2_diag.csv",
                        use_all=True, min_indep=cfg["floor"], min_stack=cfg["floor"])
    elif a.cmd == "score":
        score(a.label)
    elif a.cmd == "diag":
        diag()
    elif a.cmd == "refcheck":
        refcheck()
    else:
        report()


if __name__ == "__main__":
    main()
