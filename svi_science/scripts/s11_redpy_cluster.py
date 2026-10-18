#!/usr/bin/env python
"""s11_redpy_cluster -- PLAN §1: single-event 3C cc matrix, threshold from scrambled-control cc, Path A (average linkage,
REDPy-like) and Path B (stack assignment, CV-chosen k) sub-families; persistence; episode-locked membership.

    python s11_redpy_cluster.py --fam 10 [--workers 4]
"""
from __future__ import annotations
import os, sys, argparse, zlib
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_redpy_common as R

STATIONS = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"]
NMIN = 20


def cramers_v(lab, epi):
    t = pd.crosstab(lab, epi).to_numpy().astype(float)
    if t.shape[0] < 2 or t.shape[1] < 2:
        return np.nan
    e = t.sum(1, keepdims=True) * t.sum(0, keepdims=True) / t.sum()
    chi2 = ((t - e) ** 2 / np.maximum(e, 1e-12)).sum()
    return float(np.sqrt(chi2 / t.sum() / (min(t.shape) - 1)))


def assign(F, labs, iters=10, minsize=NMIN):
    """stack-assignment iterations: each event -> leave-self-out cluster stack with highest cc"""
    for _ in range(iters):
        ks = [k for k in np.unique(labs) if (labs == k).sum() >= 2]
        S_ = np.array([F[labs == k].sum(0) for k in ks])
        C = F @ S_.T                                                   # event . sum-stack
        own = np.array([ks.index(l) if l in ks else -1 for l in labs])
        r = np.arange(len(F))
        # leave-self-out: subtract own feature from own stack
        num = C.copy(); nrm = np.repeat(np.linalg.norm(S_, axis=1)[None], len(F), 0)
        m = own >= 0
        num[r[m], own[m]] -= (F[m] * F[m]).sum(1)
        so = S_[own[m]] - F[m]; nrm[r[m], own[m]] = np.linalg.norm(so, axis=1)
        cc = num / np.maximum(nrm, 1e-12) / np.maximum(np.linalg.norm(F, axis=1)[:, None], 1e-12)
        new = np.array(ks)[np.argmax(cc, 1)]
        # dissolve clusters below minsize into their next-best
        for k in ks:
            if 0 < (new == k).sum() < minsize:
                j = ks.index(k); cc2 = cc.copy(); cc2[:, j] = -9
                new[new == k] = np.array(ks)[np.argmax(cc2[new == k], 1)]
        if (new == labs).all():
            break
        labs = new
    return labs


def tree_init(Z, k, n):
    from scipy.cluster.hierarchy import fcluster
    for mc in range(k, 400):
        lab = fcluster(Z, mc, "maxclust")
        big = [c for c, m in zip(*np.unique(lab, return_counts=True)) if m >= NMIN]
        if len(big) >= k:
            big = sorted(big, key=lambda c: -(lab == c).sum())[:k]
            return np.where(np.isin(lab, big), lab, -1), big
    return None, None


def path_b(F, Z, kmax=5):
    """k = 1..kmax, CV (odd/even) mean assigned cc; returns chosen k, cv table, labels for chosen k"""
    from scipy.cluster.hierarchy import linkage
    n = len(F); halves = [np.arange(n) % 2 == 0, np.arange(n) % 2 == 1]
    cv = []
    fits = {}
    for k in range(1, kmax + 1):
        vals = []
        for h in range(2):
            A, B = halves[h], halves[1 - h]
            FA = F[A]
            if k == 1:
                labA = np.zeros(A.sum(), int)
            else:
                ZA = linkage(np.clip(1 - (FA @ FA.T)[np.triu_indices(len(FA), 1)], 0, 2), "average")
                lab0, big = tree_init(ZA, k, len(FA))
                if lab0 is None:
                    vals.append(np.nan); continue
                # unassigned -> nearest big stack
                Sb = np.array([FA[lab0 == c].sum(0) for c in big]); Sb /= np.linalg.norm(Sb, axis=1, keepdims=True)
                lab0 = np.where(lab0 < 0, np.array(big)[np.argmax(FA @ Sb.T, 1)], lab0)
                labA = assign(FA, lab0)
            ks = np.unique(labA)
            Sa = np.array([FA[labA == c].sum(0) for c in ks]); Sa /= np.linalg.norm(Sa, axis=1, keepdims=True)
            cB = (F[B] @ Sa.T) / np.linalg.norm(F[B], axis=1, keepdims=True).clip(1e-12)
            vals.append(float(cB.max(1).mean()))
        cv.append(dict(k=k, cv_cc=np.nanmean(vals) if np.isfinite(vals).any() else np.nan))
    cv = pd.DataFrame(cv)
    kk = 1
    for k in range(2, kmax + 1):
        g = cv.cv_cc.iloc[k - 1] - cv.cv_cc.iloc[kk - 1]
        if np.isfinite(g) and g >= 0.005:
            kk = k
        else:
            break
    if kk == 1:
        return 1, cv, np.zeros(n, int)
    lab0, big = tree_init(Z, kk, n)
    Sb = np.array([F[lab0 == c].sum(0) for c in big]); Sb /= np.linalg.norm(Sb, axis=1, keepdims=True)
    lab0 = np.where(lab0 < 0, np.array(big)[np.argmax(F @ Sb.T, 1)], lab0)
    lab = assign(F, lab0)
    m = {c: i + 1 for i, c in enumerate(sorted(np.unique(lab), key=lambda c: -(lab == c).sum()))}
    return kk, cv, np.array([m[c] for c in lab])


def run(args):
    fam, ns = args
    from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    D = R.load_series(fam, ns)
    summ, labs_out, cont, cvs, dist = [], [], [], [], []
    for ep, s in D.items():
        g, W, tS = s["g"], s["W"], s["tS"]
        F = R.feat(W, tS); n = len(F)
        C = F @ F.T; iu = np.triu_indices(n, 1)
        Wctl = s["Wca"]
        if len(Wctl) < 30:                                            # pool the station's controls across epochs
            Wctl = np.concatenate([q["Wca"] for q in D.values() if len(q["Wca"])], 0)
        Fc = R.feat(Wctl, tS); X = Fc @ F.T
        thr = float(np.percentile(X, 99))
        dist.append(dict(fam=fam, netsta=ns, epoch=ep, n=n, n_ctrl=len(Fc), thr=thr,
                         **{f"data_p{p}": float(np.percentile(C[iu], p)) for p in (50, 90, 99)},
                         **{f"ctrl_p{p}": float(np.percentile(X, p)) for p in (50, 90, 99)},
                         frac_data_pairs_above_thr=float((C[iu] > thr).mean())))
        dvec = np.clip(1 - C[iu], 0, 2)
        Z = linkage(dvec, "average")
        res = {}
        for meth in ("average", "single", "complete"):
            ZZ = Z if meth == "average" else linkage(dvec, meth)
            lab = fcluster(ZZ, 1 - thr, "distance")
            u, m = np.unique(lab, return_counts=True); big = u[m >= NMIN]
            res[meth] = (lab, big, m)
        labA_raw, bigA, _ = res["average"]
        mapA = {c: i + 1 for i, c in enumerate(sorted(bigA, key=lambda c: -(labA_raw == c).sum()))}
        labA = np.array([mapA.get(c, 0) for c in labA_raw])           # 0 = unclustered
        fracA = float((labA > 0).mean())
        kB, cv, labB = path_b(F, Z)
        cv["fam"], cv["netsta"], cv["epoch"] = fam, ns, ep; cvs.append(cv)
        use = "A" if fracA >= 0.30 else "B"
        lab = labA if use == "A" else labB
        g = g.assign(labA=labA, labB=labB, sub=lab, path=use, netsta=ns, epoch=ep)
        # leave-self-out cc with own sub-family stack (or series stack) -> selective-stacking weights
        keys = np.unique(lab); ccown = np.zeros(n)
        for k in keys:
            mk = (lab == k) if (k > 0 or use == "B") else np.ones(n, bool)
            Sk = F[mk].sum(0)
            idx = np.flatnonzero(lab == k)
            so = Sk[None] - F[idx]
            ccown[idx] = (F[idx] * so).sum(1) / np.maximum(np.linalg.norm(so, axis=1) * np.linalg.norm(F[idx], axis=1), 1e-12)
        g["cc_own"] = ccown
        labs_out.append(g[["row_id", "fam", "netsta", "epoch", "date", "t", "P", "cluster", "labA", "labB", "sub", "path", "cc_own"]])
        yrs = g.date.dt.year
        for k in [k for k in np.unique(lab) if k > 0]:
            yk = yrs[lab == k].value_counts()
            summ.append(dict(fam=fam, netsta=ns, epoch=ep, path=use, sub=int(k), size=int((lab == k).sum()),
                             years_ge3=int((yk >= 3).sum()), years_any=int(yrs.nunique()),
                             first=str(g.date[lab == k].min().date()), last=str(g.date[lab == k].max().date()),
                             mean_cc_within=float(C[np.ix_(lab == k, lab == k)][np.triu_indices((lab == k).sum(), 1)].mean())))
        # episode-locked membership
        ce = g.cluster.value_counts(); epis = ce.index[ce >= 20]
        mm = g.cluster.isin(epis).to_numpy() & (lab > 0 if use == "A" else np.ones(n, bool))
        rng = np.random.default_rng(zlib.crc32(f"{fam}{ns}{ep}".encode()))
        v0 = cramers_v(lab[mm], g.cluster.to_numpy()[mm]) if mm.sum() > 20 else np.nan
        perm = [cramers_v(rng.permutation(lab[mm]), g.cluster.to_numpy()[mm]) for _ in range(2000)] if np.isfinite(v0) else []
        p = float((np.sum(np.array(perm) >= v0) + 1) / (len(perm) + 1)) if len(perm) else np.nan
        tab = pd.crosstab(g.cluster[mm], lab[mm]); tab.index.name = "episode"
        tabl = tab.reset_index().melt(id_vars="episode", var_name="sub", value_name="count"); tabl["netsta"], tabl["epoch"], tabl["fam"] = ns, ep, fam
        cont.append(tabl)
        nA = len(bigA); nS = len(res["single"][1]); nC = len(res["complete"][1])
        dist[-1].update(path=use, nsub_A=nA, frac_A=fracA, sizes_A=";".join(str((labA == k).sum()) for k in range(1, nA + 1)),
                        nsub_single=nS, frac_single=float(np.isin(res["single"][0], res["single"][1]).mean()),
                        nsub_complete=nC, frac_complete=float(np.isin(res["complete"][0], res["complete"][1]).mean()),
                        k_B=kB, sizes_B=";".join(str((labB == k).sum()) for k in range(1, kB + 1)),
                        cramersV=v0, cramersV_p=p, cramersV_perm95=float(np.percentile(perm, 95)) if len(perm) else np.nan,
                        n_episodes_ge20=len(epis))
        # ---- figure: cc matrix ordered by dendrogram + dendrogram
        fig = plt.figure(figsize=(11, 6.5))
        ax0 = fig.add_axes([0.06, 0.73, 0.55, 0.22]); ax1 = fig.add_axes([0.06, 0.08, 0.55, 0.62]); ax2 = fig.add_axes([0.63, 0.08, 0.015, 0.62])
        ax3 = fig.add_axes([0.72, 0.55, 0.26, 0.38]); ax4 = fig.add_axes([0.72, 0.08, 0.26, 0.36])
        dn = dendrogram(Z, ax=ax0, no_labels=True, color_threshold=1 - thr, above_threshold_color="0.6", link_color_func=None)
        ax0.axhline(1 - thr, color="r", lw=0.8); ax0.set_ylabel("1 − cc"); ax0.set_xticks([])
        o = np.array(dn["leaves"])
        im = ax1.imshow(C[np.ix_(o, o)], cmap="magma", vmin=0, vmax=0.5, aspect="auto", interpolation="nearest")
        ax1.set_xlabel("events (dendrogram order)"); ax1.set_ylabel("events")
        ax2.imshow(lab[o][:, None], aspect="auto", cmap="tab10", vmin=0, vmax=9, interpolation="nearest"); ax2.set_xticks([]); ax2.set_yticks([])
        fig.colorbar(im, ax=ax4, orientation="horizontal", fraction=0.05, pad=0.25, label="single-event 3C cc")
        bins = np.linspace(-0.3, 0.8, 111)
        ax3.hist(C[iu], bins, density=True, histtype="step", color="k", label="data–data pairs")
        ax3.hist(X.ravel(), bins, density=True, histtype="step", color="tab:orange", label="scrambled control–data")
        ax3.axvline(thr, color="r", lw=0.8, label=f"threshold {thr:.3f}"); ax3.legend(fontsize=7); ax3.set_xlabel("cc (S−1…S+5 s, 3C)")
        ax3.set_yscale("log")
        for k in np.unique(lab):
            ax4.plot(g.date[lab == k], np.full((lab == k).sum(), k), "|", ms=8, label=f"{k}: {(lab == k).sum()}")
        ax4.set_ylabel("sub-family (0 = unclustered)" if use == "A" else "sub-family (path B)"); ax4.set_title("membership in time", fontsize=8)
        fig.suptitle(f"Family {fam} {ns} {ep}: {n} events P≥0.5; path {use}; A: {nA} sub-families ≥{NMIN} ({100*fracA:.0f} % of events); B: k={kB}", fontsize=10)
        fig.savefig(f"{R.FIG}/cluster_fam{fam}_{ns}_{ep}.png", dpi=110); plt.close(fig)
    return pd.DataFrame(dist), pd.concat(labs_out), pd.concat(cont) if cont else pd.DataFrame(), pd.concat(cvs), pd.DataFrame(summ)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--fam", type=int, default=10); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args(); R.ensure()
    from multiprocessing import Pool
    with Pool(min(a.workers, 8)) as p:
        res = p.map(run, [(a.fam, ns) for ns in STATIONS])
    pd.concat([r[0] for r in res]).to_csv(f"{R.TB}/subfam_summary_fam{a.fam}.csv", index=False)
    pd.concat([r[1] for r in res]).to_parquet(f"{R.TB}/subfam_labels_fam{a.fam}.parquet", index=False)
    pd.concat([r[2] for r in res]).to_csv(f"{R.TB}/subfam_episode_fam{a.fam}.csv", index=False)
    pd.concat([r[3] for r in res]).to_csv(f"{R.TB}/subfam_pathB_cv_fam{a.fam}.csv", index=False)
    pd.concat([r[4] for r in res]).to_csv(f"{R.TB}/subfam_members_fam{a.fam}.csv", index=False)
    print(pd.concat([r[4] for r in res]).to_string())
    print(pd.concat([r[0] for r in res]).to_string())


if __name__ == "__main__":
    main()
